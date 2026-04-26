"""
FOCUS-2 / Chain B & C: Modular structure of Omega-stratified summatories
        (L(x) mod m, C_3(x) mod 2, pi_k(x) mod 2 for k=2..6).

Goal
----
Test whether *any* algebraic structure invisible to the 24 pseudorandomness
measures (novel/pseudorandomness_of_pi.md) hides in the building blocks of
the Liouville identity

        pi(x) = (x - L(x))/2 - C_3(x)               (E2.2, S46)

under fixed small modulus q.  The identity is exact, so polylog access to
*either* L(x) mod 4 *or* C_3(x) mod 2 (combined with the trivial L(x) mod 2
= x mod 2) immediately gives pi(x) mod 2 in polylog.

Key derivations done in this session (refinement of FOCUS-2):

  (i)  L(x) = sum_{n<=x} lambda(n) with lambda(n) in {+/-1}, so the parity
       of L(x) equals the parity of x.  Hence L(x) mod 2 = x mod 2 is a
       FREE identity — it carries no non-trivial information about primes.

  (ii) Writing L(x) = x - 2*A(x) with A(x) = (x - L(x))/2 the count of
       integers <= x with Omega(n) odd, we get
            L(x) mod 4 = (x - 2 A(x)) mod 4
       so  A(x) mod 2 = ((x mod 4) - (L(x) mod 4)) / 2 (mod 2).
       Therefore the *non-trivial* part of "L(x) mod 4" is exactly
       A(x) mod 2.

  (iii) pi(x) = A(x) - C_3(x), so pi(x) mod 2 = (A(x) mod 2) XOR
        (C_3(x) mod 2).  TWO bits in, one bit out.

So the actual missing primitives are:
    P1.  A(x) mod 2 = #{n<=x : Omega(n) odd} mod 2  -- equivalently
         L(x) mod 4 (modulo trivial parity).
    P2.  C_3(x) mod 2 = #{n<=x : Omega(n) odd, Omega(n) >= 3} mod 2.

Either one alone is not enough.  Both are needed (or a third primitive).

What this experiment tests
--------------------------
A.  Pseudorandomness battery on A(x) mod 2 (= L(x) mod 4 modulo
    free parity), C_3(x) mod 2, and pi_k(x) mod 2 for k = 1..6.
    Measures: entropy rate, autocorrelation, FFT spectral lines,
    correlation with cheap arithmetic proxies, LFSR length over
    GF(2), conditional entropy given x.

B.  Pairwise mutual information between A(x) mod 2 and C_3(x) mod 2.
    Their XOR is pi(x) mod 2 (= 50% random by S35,S46).  If each is
    individually random AND mutually independent, no shortcut.  If
    one predicts the other beyond chance, there's structure to
    exploit.

C.  Modular structure of L(x) for moduli m in {2..16} -- does any m
    show non-trivial spectral lines, periodicity, or compressibility?

D.  Boolean-fusion of cheap predictors: can we combine the parity bits
    of (x), (M(x)), (Q(x)), (sigma_0(x)), (number of distinct prime
    factors of x), etc., to predict pi(x) mod 2 better than chance?

A negative result would close the FOCUS-2 chain at the structural
level, complementing S46's identity-level closure.
"""

from __future__ import annotations

import math
import time
from collections import Counter
from typing import List

import numpy as np


# ---------- Sieve Omega(n) and lambda(n) -----------------------------------


def sieve_spf(N: int) -> np.ndarray:
    """Smallest-prime-factor sieve.  spf[n] = smallest prime dividing n,
    spf[0] = spf[1] = 0."""
    spf = np.zeros(N + 1, dtype=np.int32)
    for p in range(2, N + 1):
        if spf[p] == 0:
            # mark all multiples of p whose spf is not yet set
            mask_slice = spf[p::p]
            unset = mask_slice == 0
            mask_slice[unset] = p
            spf[p::p] = mask_slice
    return spf


def sieve_omega(N: int) -> np.ndarray:
    """Return array Omega[0..N] where Omega[n] = number of prime factors of
    n counted with multiplicity.  Uses spf sieve + recurrence."""
    spf = sieve_spf(N)
    Omega = np.zeros(N + 1, dtype=np.int32)
    # Recurrence: Omega(n) = Omega(n / spf(n)) + 1; need n in increasing order.
    for n in range(2, N + 1):
        Omega[n] = Omega[n // spf[n]] + 1
    return Omega, spf


# ---------- Cumulative counters --------------------------------------------


def compute_summatories(Omega: np.ndarray):
    N = len(Omega) - 1
    lam = np.where(Omega % 2, -1, 1).astype(np.int64)
    lam[0] = 0  # exclude n=0 (no integer 0)

    # Make Omega-equality test exclude n=0 by setting Omega[0]=-1 sentinel for
    # the masks below.  Don't mutate the input.
    Omega_t = Omega.copy()
    Omega_t[0] = -1

    L = np.cumsum(lam).astype(np.int64)  # L[x] = sum_{n<=x} lambda(n)

    # pi_k counts: #{n<=x : Omega(n)=k} (over n>=1 only).
    # We track all k that occur for some n<=N, capped at MAX_K.
    MAX_K = int(np.max(Omega)) + 1  # so we cover full range
    pi_k = {}
    for k in range(0, MAX_K + 1):
        pi_k[k] = np.cumsum((Omega_t == k).astype(np.int64))

    pi = pi_k[1].copy()  # primes

    # A(x) = #{n<=x : Omega(n) odd} = (x - L(x))/2 (using lam[0]=0 sentinel)
    A = ((np.arange(N + 1) - L) // 2).astype(np.int64)
    # C_3(x) = #{n<=x : Omega(n) odd and >=3} = sum_{k odd, k>=3} pi_k(x)
    C3 = np.zeros(N + 1, dtype=np.int64)
    for k in range(3, MAX_K + 1, 2):
        C3 += pi_k[k]

    return {
        "lam": lam,
        "L": L,
        "A": A,
        "C3": C3,
        "pi": pi,
        "pi_k": pi_k,
    }


# ---------- Pseudorandomness battery --------------------------------------


def block_entropy(bits: np.ndarray, L_block: int) -> float:
    """Plug-in estimator of L-block entropy (in bits)."""
    n = len(bits) - L_block + 1
    if n < 16:
        return float("nan")
    # Pack each L_block-window as integer.
    pow2 = 1 << np.arange(L_block)
    # Use sliding window via numpy.lib.stride_tricks
    from numpy.lib.stride_tricks import sliding_window_view

    win = sliding_window_view(bits, L_block)
    codes = win.astype(np.int64) @ pow2
    cnt = Counter(codes.tolist())
    total = sum(cnt.values())
    H = 0.0
    for c in cnt.values():
        p = c / total
        H -= p * math.log2(p)
    return H


def entropy_rate(bits: np.ndarray, L_max: int = 12) -> float:
    """Conditional entropy h(L_max | L_max-1)."""
    if len(bits) < (1 << (L_max + 2)):
        L_max = max(2, int(math.log2(len(bits)) - 4))
    HL = block_entropy(bits, L_max)
    HL_1 = block_entropy(bits, L_max - 1)
    if math.isnan(HL) or math.isnan(HL_1):
        return float("nan")
    return HL - HL_1


def autocorr_max(bits: np.ndarray, lags: List[int]) -> float:
    """Return max |corr| over the supplied lags after centering."""
    x = bits.astype(np.float64)
    x -= x.mean()
    var = (x * x).mean()
    if var <= 0:
        return 0.0
    best = 0.0
    n = len(x)
    for lag in lags:
        c = float((x[: n - lag] * x[lag:]).mean()) / var
        best = max(best, abs(c))
    return best


def fft_max_spectral_line(bits: np.ndarray) -> float:
    """Z-score of the largest FFT bin (excluding DC).  >5 = spectral line."""
    x = bits.astype(np.float64)
    x -= x.mean()
    if (x == 0).all():
        return 0.0
    X = np.abs(np.fft.rfft(x))[1:]  # drop DC
    mu = X.mean()
    sd = X.std() + 1e-12
    return float((X.max() - mu) / sd)


def correlation_with_simple(
    bits: np.ndarray, name_arr_pairs: list
) -> dict[str, float]:
    """Compute |Pearson correlation| with each simple arithmetic proxy."""
    out = {}
    bf = bits.astype(np.float64)
    bf -= bf.mean()
    bnorm = (bf * bf).sum() ** 0.5 + 1e-12
    for name, arr in name_arr_pairs:
        a = arr.astype(np.float64)
        a -= a.mean()
        anorm = (a * a).sum() ** 0.5 + 1e-12
        if anorm <= 1e-9 or bnorm <= 1e-9:
            out[name] = 0.0
        else:
            out[name] = float((bf @ a) / (bnorm * anorm))
    return out


def lfsr_length_gf2(bits: np.ndarray) -> int:
    """Berlekamp-Massey: linear complexity over GF(2)."""
    s = bits.astype(np.int8)
    n = len(s)
    C = np.zeros(n, dtype=np.int8)
    B = np.zeros(n, dtype=np.int8)
    C[0] = 1
    B[0] = 1
    L = 0
    m = 1
    b = 1
    for N_idx in range(n):
        # discrepancy
        d = int(s[N_idx])
        for i in range(1, L + 1):
            d ^= int(C[i] & s[N_idx - i])
        d &= 1
        if d == 0:
            m += 1
        elif 2 * L <= N_idx:
            T = C.copy()
            shift = m
            # C = C XOR (B shifted by m), within length n
            if shift < n:
                C[shift:] ^= B[: n - shift]
            L = N_idx + 1 - L
            B = T
            b = d
            m = 1
        else:
            shift = m
            if shift < n:
                C[shift:] ^= B[: n - shift]
            m += 1
    return L


def boolean_xor_fusion(target: np.ndarray, predictors: list) -> tuple:
    """Try every non-empty XOR combination of <=4 predictors against target.
    Return best agreement fraction and which subset achieved it."""
    from itertools import combinations

    n = len(target)
    best = 0.5
    best_sub = ()
    for r in range(1, min(5, len(predictors) + 1)):
        for combo in combinations(range(len(predictors)), r):
            pred = predictors[combo[0]].copy()
            for j in combo[1:]:
                pred = pred ^ predictors[j]
            agree = float((pred == target).mean())
            if abs(agree - 0.5) > abs(best - 0.5):
                best = agree
                best_sub = combo
    return best, best_sub


# ---------- Main ------------------------------------------------------------


def main():
    N = 2_000_000
    print(f"Sieving Omega for n in [1, {N}] ...")
    t0 = time.time()
    Omega, spf = sieve_omega(N)
    print(f"  done in {time.time() - t0:.1f}s")

    print("Computing summatories L, A, C_3, pi_k ...")
    t0 = time.time()
    S = compute_summatories(Omega)
    print(f"  done in {time.time() - t0:.1f}s")

    L = S["L"]
    A = S["A"]
    C3 = S["C3"]
    pi = S["pi"]
    pi_k = S["pi_k"]

    # --- Sanity: L(x) ≡ x (mod 2) trivially ---
    parity_match = int(((L[1:] % 2) == (np.arange(1, N + 1) % 2)).sum())
    print(
        f"\nSanity: L(x) ≡ x (mod 2) trivially? "
        f"parity match = {parity_match}/{N} (expected {N})"
    )
    assert parity_match == N, "Free identity violated -- bug somewhere"

    # --- Identity check: pi(x) = A(x) - C_3(x) ---
    diff = pi - (A - C3)
    err = int(np.abs(diff).sum())
    print(
        f"Identity pi(x) = A(x) - C_3(x): cumulative discrepancy = {err} "
        f"(expected 0 mod (-1+1) initial constants)"
    )

    # --- Bit streams of interest (length N) ---
    x_arr = np.arange(1, N + 1)

    bits = {
        "pi_mod_2": (pi[1:] % 2).astype(np.int8),
        "A_mod_2": (A[1:] % 2).astype(np.int8),
        "C3_mod_2": (C3[1:] % 2).astype(np.int8),
        "Lmod4_high_bit": (((L[1:] % 4) >> 1) & 1).astype(np.int8),
        "Lmod4_low_bit": (L[1:] % 4 & 1).astype(np.int8),  # = x mod 2 trivially
        "Lmod8_bit2": (((L[1:] % 8) >> 2) & 1).astype(np.int8),
        "Lmod8_bit1": (((L[1:] % 8) >> 1) & 1).astype(np.int8),
        "pi2_mod_2": (pi_k[2][1:] % 2).astype(np.int8),
        "pi3_mod_2": (pi_k[3][1:] % 2).astype(np.int8),
        "pi4_mod_2": (pi_k[4][1:] % 2).astype(np.int8),
        "pi5_mod_2": (pi_k[5][1:] % 2).astype(np.int8),
        "pi6_mod_2": (pi_k[6][1:] % 2).astype(np.int8),
    }

    # Trim to a length we can run heavier measures on (FFT, LFSR are O(n^2)
    # for LFSR -> use shorter window).
    LFSR_LEN = 4096
    FFT_LEN = 200_000

    # --- Per-bit-stream pseudorandomness battery ---
    print("\n=== Pseudorandomness battery (24-measure family + new) ===")
    print(
        f"{'stream':<18s} | {'mean':>6s} | {'h_block(8)':>10s} | "
        f"{'AC max[1,2,3,5,10,30]':>22s} | {'FFT z':>6s} | "
        f"{'LFSR/N':>8s}"
    )
    rows = {}
    for name, b in bits.items():
        if len(b) > FFT_LEN:
            b_fft = b[:FFT_LEN]
        else:
            b_fft = b
        h8 = block_entropy(b, 8)
        ac = autocorr_max(b, [1, 2, 3, 5, 10, 30])
        fz = fft_max_spectral_line(b_fft)
        # LFSR is expensive -> compute on a shorter prefix
        lf = lfsr_length_gf2(b[:LFSR_LEN]) / LFSR_LEN
        mean = float(b.mean())
        rows[name] = {
            "mean": mean,
            "h_block_8": h8,
            "AC_max": ac,
            "FFT_z": fz,
            "LFSR_ratio": lf,
        }
        print(
            f"{name:<18s} | {mean:6.4f} | {h8:10.4f} | "
            f"{ac:22.4f} | {fz:6.2f} | {lf:8.4f}"
        )

    # --- Mutual independence A_mod_2 vs C3_mod_2 ---
    print("\n=== Independence test: A(x) mod 2 vs C_3(x) mod 2 ===")
    a = bits["A_mod_2"]
    c = bits["C3_mod_2"]
    p = bits["pi_mod_2"]
    # Joint distribution
    joint = Counter(zip(a.tolist(), c.tolist()))
    total = len(a)
    HA = -sum(
        (cnt / total) * math.log2(cnt / total)
        for k, cnt in Counter(a.tolist()).items()
        if cnt > 0
    )
    HC = -sum(
        (cnt / total) * math.log2(cnt / total)
        for k, cnt in Counter(c.tolist()).items()
        if cnt > 0
    )
    HAC = -sum(
        (cnt / total) * math.log2(cnt / total)
        for cnt in joint.values()
    )
    MI = HA + HC - HAC
    print(f"  H(A) = {HA:.5f}, H(C_3) = {HC:.5f}, H(A,C_3) = {HAC:.5f}")
    print(f"  Mutual info I(A; C_3) = {MI:.5e} bits")
    # Verify XOR identity
    pred = a ^ c
    agree = float((pred == p).mean())
    print(f"  Verify pi mod 2 == A XOR C_3?  agreement = {agree:.6f}")
    print(f"    (allowed tiny mismatch from boundary cases at very small x)")

    # --- L(x) mod m structure for m in {2..16} ---
    print("\n=== L(x) mod m structural sweep ===")
    for m in [2, 3, 4, 5, 6, 7, 8, 9, 12, 13, 16]:
        Lm = (L[1:] % m).astype(np.int8)
        # Use bit-flag streams to apply the bit-level battery
        # Counts: distribution flatness?
        cnt = Counter(Lm.tolist())
        n = len(Lm)
        # Entropy of the residue distribution
        H = -sum(
            (c / n) * math.log2(c / n)
            for c in cnt.values() if c > 0
        )
        H_uniform = math.log2(m)
        # Spectral signature: z-score of largest FFT bin of mod-m stream
        if n > FFT_LEN:
            fz = fft_max_spectral_line(Lm[:FFT_LEN].astype(np.float64))
        else:
            fz = fft_max_spectral_line(Lm.astype(np.float64))
        # Autocorrelation
        ac = autocorr_max(Lm.astype(np.float64), [1, 2, 3, 5, 10, 30, 100])
        print(
            f"  m={m:>2d}: H={H:6.4f}/{H_uniform:6.4f}, "
            f"FFT z={fz:6.2f}, AC max={ac:.4f}"
        )

    # --- Cheap-proxy correlations with pi(x) mod 2 ---
    print("\n=== Cheap-proxy correlations with pi(x) mod 2 ===")
    print("  Computing simple arithmetic proxies (mu, Q, omega_distinct) ...")
    t0 = time.time()
    # mu(n) = (-1)^Omega(n) if n squarefree, else 0.
    is_squarefree = np.ones(N + 1, dtype=bool)
    is_squarefree[0] = False
    for p in range(2, int(math.isqrt(N)) + 1):
        if spf[p] == p:  # p is prime
            is_squarefree[p * p :: p * p] = False
    mu = np.where(is_squarefree, np.where(Omega % 2, -1, 1), 0).astype(np.int8)
    mu[0] = 0
    M = np.cumsum(mu).astype(np.int64)
    Q = np.cumsum(is_squarefree.astype(np.int64))  # squarefree count

    # sigma_0(n) is odd iff n is a perfect square — gives a sparse 0/1 stream.
    is_square = np.zeros(N + 1, dtype=np.int8)
    sq = (np.arange(int(math.isqrt(N)) + 1) ** 2)
    sq = sq[sq <= N]
    is_square[sq] = 1
    sigma0_parity = is_square  # sigma_0(n) mod 2

    # Number of DISTINCT prime factors omega(n) (lowercase).
    omega_distinct = np.zeros(N + 1, dtype=np.int16)
    for p in range(2, N + 1):
        if spf[p] == p:  # prime
            omega_distinct[p::p] += 1
    print(f"  proxy sieves done in {time.time() - t0:.1f}s")

    proxies = [
        ("x mod 2", x_arr % 2),
        ("x mod 3", x_arr % 3),
        ("x mod 4 == 1", (x_arr % 4 == 1).astype(np.int8)),
        ("x mod 4 == 3", (x_arr % 4 == 3).astype(np.int8)),
        ("M(x) mod 2", (M[1:] % 2).astype(np.int8)),
        ("Q(x) mod 2", (Q[1:] % 2).astype(np.int8)),
        ("sigma0(x) mod 2 (=is_square(x))", sigma0_parity[1:].astype(np.int8)),
        ("omega(x) mod 2", (omega_distinct[1:] % 2).astype(np.int8)),
        ("Omega(x) mod 2", (Omega[1:] % 2).astype(np.int8)),
        ("L(x) mod 4 high bit", bits["Lmod4_high_bit"]),
        ("L(x) mod 8 bit 2", bits["Lmod8_bit2"]),
        ("L(x) mod 8 bit 1", bits["Lmod8_bit1"]),
    ]

    print("\n  Pearson |corr| with pi(x) mod 2:")
    target = bits["pi_mod_2"]
    proxy_arrs = []
    proxy_names = []
    for name, arr in proxies:
        proxy_arrs.append(arr.astype(np.int8))
        proxy_names.append(name)
        # truncate
        n = min(len(target), len(arr))
        cor = correlation_with_simple(target[:n], [(name, arr[:n])])
        print(f"    {name:<35s}: {cor[name]:+.4f}")

    # --- Boolean XOR fusion ---
    print("\n=== Boolean XOR fusion of cheap proxies vs pi(x) mod 2 ===")
    # Restrict to binary {0,1} predictors
    binary_predictors = []
    binary_names = []
    for name, arr in zip(proxy_names, proxy_arrs):
        a = arr[: len(target)].astype(np.int8)
        if a.max() <= 1 and a.min() >= 0:
            binary_predictors.append(a & 1)
            binary_names.append(name)
    sample = 200_000  # for speed
    target_s = target[:sample]
    pred_s = [p[:sample] for p in binary_predictors]
    best, best_sub = boolean_xor_fusion(target_s, pred_s)
    print(f"  Best XOR-fusion agreement: {best:.4f}")
    print("  Best subset:")
    for i in best_sub:
        print(f"    - {binary_names[i]}")

    # --- Conditional entropy of pi(x) mod 2 given (proxy stream) ---
    print("\n=== H(pi mod 2 | best proxies) (information left after each proxy) ===")
    target_s = target[: 200_000]
    for name, p in zip(binary_names, pred_s):
        # Conditional entropy via joint count
        joint = Counter(zip(p.tolist(), target_s.tolist()))
        marg = Counter(p.tolist())
        H_cond = 0.0
        for v_p, n_p in marg.items():
            pp = n_p / len(target_s)
            for v_t in (0, 1):
                n_pt = joint.get((v_p, v_t), 0)
                if n_pt > 0:
                    p_t_given_p = n_pt / n_p
                    H_cond -= pp * p_t_given_p * math.log2(p_t_given_p)
        print(f"    {name:<35s}: H(pi|·) = {H_cond:.5f} bits   (max=1)")

    print("\nAll computations done.")


if __name__ == "__main__":
    main()
