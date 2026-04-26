"""
FOCUS-2 extension (S56): Character-twisted Liouville sums L_chi(x) for
non-trivial Dirichlet characters chi mod q in {3, 5, 7, 11, 13}.

Background
----------
S55 closed FOCUS-2 for q = 2 by:
  (i) the FREE identity  L(x) mod 2 = x mod 2  (no prime info);
  (ii) the structural pseudorandomness of the "next bit"
       A(x) mod 2 = (x - L(x))/2 mod 2;
  (iii) mutual independence of A(x) mod 2 and C_3(x) mod 2.

The Chain B target requires polylog access to pi(x) mod q for q in
{3, 5, 7, 11, 13}.  TODO.md amendment: test the natural generalization,
character-twisted Liouville sums

        L_chi(x) := sum_{n<=x} lambda(n) * chi(n)

for non-trivial Dirichlet character chi mod q.  Two questions:

  Q1.  Is there a free identity for the parity of L_chi(x), generalizing
       L(x) mod 2 = x mod 2?

  Q2.  If yes, is the "next bit" structurally pseudorandom (so the
       attack collapses for the same reason as q = 2)?

Analytic precheck
-----------------
For ANY Dirichlet character chi mod q:

  lambda(n) in {-1, +1} for every n >= 1.
  chi(n)   in {0} U {complex roots of unity}.
  lambda(n) * chi(n) is in {0} U {roots of unity}.

For real characters (Legendre symbol type, chi(n) in {-1, 0, +1}):

  L_chi(x) is an integer summing terms in {-1, 0, +1}.
  L_chi(x) mod 2 = (#nonzero summands) mod 2
                 = #{n <= x : gcd(n, q) = 1} mod 2
                 = TRIVIAL FLOOR-COUNT (computable in O(polylog) trivially,
                                        carries no prime info).

For complex characters of order k >= 3, decompose
  chi = chi_re + i * chi_im, both real-valued in {-1, 0, +1, ...}
  Re(L_chi(x)) = sum lambda(n) chi_re(n)  ∈ Z
  Im(L_chi(x)) = sum lambda(n) chi_im(n)  ∈ Z

Their parities collapse similarly:
  Re(L_chi(x)) mod 2 = #{n : chi_re(n) is odd} mod 2
                     = #{n <= x : chi(n) ∈ {real *odd*-magnitude unit}} mod 2

For Dirichlet characters whose values lie in {0, +/-1, +/-i, ...} this
remains a trivial floor count.  In general, chi(n) at gcd(n,q)=1 is a
phi(q)-th root of unity; for k = order of chi, chi(n) takes values in
the cyclic group of k-th roots; the (mod 2) parity of any *integer*
linear functional of L_chi(x) collapses to a count over a fixed coset
union mod q.

So we EXPECT the analog of the S55 free identity for every chi.  We
verify this empirically, then test pseudorandomness of the next bit.

Experiment
----------
For each prime q in {3, 5, 7, 11, 13}:
  enumerate Dirichlet characters chi mod q (phi(q) = q - 1 of them).
  Compute L_chi(x) for x in [1, N] by sieve.
  For each chi:
    a. Verify the predicted free identity for parity (test Re and Im).
    b. Compute A_chi(x) = #{n <= x : lambda(n) chi(n) = +1}, the
       "next-bit" sequence; test A_chi(x) mod 2 with the
       pseudorandomness battery (block entropy, autocorrelation,
       FFT spectral line, LFSR length).
    c. Compute mutual information I(A_chi mod 2 ; pi(x; q, a) mod 2)
       for each a in (Z/q)^*.

If each free identity holds AND each next-bit is pseudorandom AND
MI is near zero, the character-twisted Liouville attack is closed for
q in {3, 5, 7, 11, 13}, completing FOCUS-2 in tandem with S55.

Outputs go to character_twisted_liouville_results.md alongside.
"""

from __future__ import annotations

import math
import time
from collections import Counter
from typing import List, Dict, Tuple

import numpy as np


# ---------- Sieve helpers (same shape as liouville_modular_structure.py) ---


def sieve_spf(N: int) -> np.ndarray:
    spf = np.zeros(N + 1, dtype=np.int32)
    for p in range(2, N + 1):
        if spf[p] == 0:
            mask_slice = spf[p::p]
            unset = mask_slice == 0
            mask_slice[unset] = p
            spf[p::p] = mask_slice
    return spf


def sieve_omega(N: int):
    spf = sieve_spf(N)
    Omega = np.zeros(N + 1, dtype=np.int32)
    for n in range(2, N + 1):
        Omega[n] = Omega[n // spf[n]] + 1
    return Omega, spf


def sieve_primes(N: int) -> np.ndarray:
    """Return boolean array prime[0..N]."""
    is_prime = np.ones(N + 1, dtype=bool)
    is_prime[0] = is_prime[1] = False
    for p in range(2, int(math.isqrt(N)) + 1):
        if is_prime[p]:
            is_prime[p * p :: p] = False
    return is_prime


# ---------- Dirichlet characters mod prime q ------------------------------


def primitive_root(q: int) -> int:
    """Smallest primitive root mod prime q."""
    phi = q - 1
    factors = []
    n = phi
    d = 2
    while d * d <= n:
        if n % d == 0:
            factors.append(d)
            while n % d == 0:
                n //= d
        d += 1
    if n > 1:
        factors.append(n)
    for g in range(2, q):
        if all(pow(g, phi // p, q) != 1 for p in factors):
            return g
    raise ValueError(f"no primitive root mod {q}")


def char_table(q: int) -> List[np.ndarray]:
    """Return list of (q-1) Dirichlet character functions chi[0..q-1] as
    complex arrays.  chi[0] is the principal character (= 1 for gcd=1, 0
    otherwise); chi[k] is e^{2 pi i k log_g(n) / (q-1)} for primitive root g.

    Index k = 0 .. q-2.
    """
    g = primitive_root(q)
    phi = q - 1
    # log_g(n) for n in (Z/qZ)^*
    log_g = {1: 0}
    cur = 1
    for j in range(1, phi):
        cur = (cur * g) % q
        log_g[cur] = j
    chars = []
    for k in range(phi):
        chi = np.zeros(q, dtype=np.complex128)
        for n in range(1, q):
            chi[n] = np.exp(2j * np.pi * k * log_g[n] / phi)
        chars.append(chi)
    return chars, log_g, g


def lift_chi_to_N(chi_q: np.ndarray, N: int) -> np.ndarray:
    """chi(n) for n in 0..N by chi(n) = chi_q[n mod q]."""
    q = len(chi_q)
    idx = np.arange(N + 1) % q
    return chi_q[idx]


# ---------- Pseudorandomness battery (subset of S55) ----------------------


def block_entropy(bits: np.ndarray, L_block: int) -> float:
    n = len(bits) - L_block + 1
    if n < 16:
        return float("nan")
    pow2 = 1 << np.arange(L_block)
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


def autocorr_max(bits: np.ndarray, lags: List[int]) -> float:
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
    x = bits.astype(np.float64)
    x -= x.mean()
    if (x == 0).all():
        return 0.0
    X = np.abs(np.fft.rfft(x))[1:]
    mu = X.mean()
    sd = X.std() + 1e-12
    return float((X.max() - mu) / sd)


def lfsr_length_gf2(bits: np.ndarray) -> int:
    s = bits.astype(np.int8)
    n = len(s)
    C = np.zeros(n, dtype=np.int8)
    B = np.zeros(n, dtype=np.int8)
    C[0] = 1
    B[0] = 1
    L = 0
    m = 1
    for N_idx in range(n):
        d = int(s[N_idx])
        for i in range(1, L + 1):
            d ^= int(C[i] & s[N_idx - i])
        d &= 1
        if d == 0:
            m += 1
        elif 2 * L <= N_idx:
            T = C.copy()
            shift = m
            if shift < n:
                C[shift:] ^= B[: n - shift]
            L = N_idx + 1 - L
            B = T
            m = 1
        else:
            shift = m
            if shift < n:
                C[shift:] ^= B[: n - shift]
            m += 1
    return L


def mutual_info_bits(a: np.ndarray, b: np.ndarray) -> float:
    """MI in bits between two binary streams of equal length."""
    n = len(a)
    joint = Counter(zip(a.tolist(), b.tolist()))
    ma = Counter(a.tolist())
    mb = Counter(b.tolist())
    HA = -sum((c / n) * math.log2(c / n) for c in ma.values() if c > 0)
    HB = -sum((c / n) * math.log2(c / n) for c in mb.values() if c > 0)
    HAB = -sum((c / n) * math.log2(c / n) for c in joint.values() if c > 0)
    return HA + HB - HAB


# ---------- Free-identity prediction --------------------------------------
#
# For ANY Dirichlet character chi mod q of any order d, the value of
# L_chi(x) lies in the ring of integers Z[zeta_d].  Its mod-2 reduction
# (in that ring) is determined by:
#
#   L_chi(x) mod 2 = sum_{r in (Z/q)*} chi(r) * (S_r(x) mod 2)         (*)
#
# where S_r(x) = sum_{n<=x, n=r mod q} lambda(n).  Since each summand
# of S_r(x) is +/-1 (so =1 mod 2), we have
#
#   S_r(x) mod 2 = #{n<=x : n=r mod q} mod 2,
#
# a TRIVIAL floor-count quantity.  So (*) gives a closed-form prediction
# for L_chi(x) mod 2 that uses no prime information.
#
# We verify the free identity by computing both sides of (*) at every x
# in [1,N] and checking actual = predicted (in the appropriate ring).
# In practice we work in C and check that Re(actual - predicted) and
# Im(actual - predicted) are always even integers (resp. even half-
# integers in the order-3/6 case where the basis half-step is allowed).


def predict_chi_mod2(chi_q: np.ndarray, N: int) -> np.ndarray:
    """Predicted L_chi(x) mod 2 (as a complex number) under the free
    identity (*) above.  Returned shape: (N,) complex array indexed by
    x in [1..N].  Each entry is sum_{r in (Z/q)*} chi(r) * (count_r(x) mod 2).
    """
    q = len(chi_q)
    # count_r(x) mod 2 for each r in 1..q-1
    flags = []
    for r in range(1, q):
        # in_progression[n] = 1 iff n mod q == r and n >= 1
        in_prog = (np.arange(1, N + 1) % q == r).astype(np.int64)
        cnt = np.cumsum(in_prog) % 2  # count_r(x) mod 2
        flags.append(cnt.astype(np.int8))
    # predicted = sum_r chi(r) * (count_r(x) mod 2)
    pred = np.zeros(N, dtype=np.complex128)
    for idx, r in enumerate(range(1, q)):
        pred += chi_q[r] * flags[idx]
    return pred


def cyclotomic_reduction(d: int) -> np.ndarray:
    """Return the integer matrix M of shape (d, phi(d)) such that
    zeta_d^k = M[k] . [1, zeta_d, ..., zeta_d^{phi(d)-1}]  (as Z-module elt).

    Uses the cyclotomic polynomial Phi_d to reduce.  We compute Phi_d
    by polynomial division of x^d - 1 by all lower Phi_e (e | d, e < d).
    """
    # Get all divisors of d.
    divs = [e for e in range(1, d + 1) if d % e == 0]
    # Recursively compute Phi_e for each divisor.
    phi_polys = {}
    for e in divs:
        # Start with x^e - 1
        num = [0] * (e + 1)
        num[0] = -1
        num[e] = 1
        # Divide by all Phi_f for f | e, f < e (already computed).
        for f in divs:
            if f >= e:
                break
            if e % f == 0:
                num = poly_div_int(num, phi_polys[f])
        phi_polys[e] = num
    Phi_d = phi_polys[d]
    phi_d = len(Phi_d) - 1
    # Build M[k] for k in 0..d-1.
    # We start with [zeta^0, ..., zeta^{phi_d-1}] = identity rows.
    M = np.zeros((d, phi_d), dtype=np.int64)
    for j in range(phi_d):
        M[j, j] = 1
    # For k >= phi_d, reduce zeta^k = -(Phi_d/lead - x^{phi_d}) lower-degree.
    # I.e., zeta^{phi_d} = - Phi_d_low_part (since Phi_d is monic with leading 1).
    # In general, zeta^k = zeta^{k - phi_d} * zeta^{phi_d}, then expand.
    # Simpler: maintain M[k] by recurrence M[k] = shift(M[k-1]) reduced.
    for k in range(phi_d, d):
        prev = M[k - 1]
        new = np.zeros(phi_d, dtype=np.int64)
        new[1:] = prev[: phi_d - 1]
        # The dropped coefficient was prev[phi_d - 1] (coefficient of zeta^{phi_d - 1});
        # multiplying by zeta sends this to coefficient of zeta^{phi_d}, which we
        # reduce using zeta^{phi_d} = -sum_{i<phi_d} Phi_d[i] * zeta^i.
        carry = prev[phi_d - 1]
        if carry != 0:
            for i in range(phi_d):
                new[i] -= carry * Phi_d[i]
        M[k] = new
    return M


def poly_div_int(num: list, den: list) -> list:
    """Polynomial division (exact) of integer-coeff polynomials.  Returns
    quotient as a list of coefficients (low-degree first).  Asserts zero
    remainder.  Polynomials are stored low-degree first."""
    num = list(num)
    den = list(den)
    while len(den) > 1 and den[-1] == 0:
        den.pop()
    deg_d = len(den) - 1
    if deg_d < 0:
        raise ValueError
    quot = [0] * (len(num) - deg_d)
    for k in range(len(quot) - 1, -1, -1):
        if num[k + deg_d] == 0:
            continue
        coef = num[k + deg_d] // den[deg_d]
        if coef * den[deg_d] != num[k + deg_d]:
            # Should not happen for cyclotomic reductions.
            raise ValueError("non-exact polynomial division")
        quot[k] = coef
        for j in range(deg_d + 1):
            num[k + j] -= coef * den[j]
    # Remainder check
    for v in num[:deg_d]:
        if v != 0:
            raise ValueError("non-zero remainder in poly_div_int")
    return quot


def free_identity_holds_exact(chi_q: np.ndarray, lam: np.ndarray,
                              N: int, n_sample: int = 2000) -> tuple[bool, int, float]:
    """EXACT-arithmetic check that the free identity

         L_chi(x) ≡ sum_r chi(r) count_r(x) (mod 2 Z[zeta_d])

    holds at n_sample evenly-spaced x-values in [1..N].

    Implementation: represent each chi(r) as a Z^{phi(d)}-vector in the
    basis {1, zeta_d, ..., zeta_d^{phi(d)-1}}.  Then both sides are
    vectors in Z^{phi(d)}; check componentwise mod-2 equality.

    Returns (all_match, n_fail, max_residual_in_basis_decomposition).
    """
    q = len(chi_q)
    # Determine character order d
    d = 1
    cur = chi_q.copy()
    chi_pri = np.zeros_like(chi_q)
    chi_pri[1:] = 1.0
    while not np.allclose(cur, chi_pri, atol=1e-10):
        cur = cur * chi_q
        d += 1
        if d > 60:
            raise RuntimeError("character order not found")

    # Phi_d basis reduction
    M = cyclotomic_reduction(d)  # shape (d, phi_d)
    phi_d = M.shape[1]
    # For each r in (Z/q)*, find e_r such that chi(r) = zeta_d^{e_r}.
    e_r_map = {}
    max_resid = 0.0
    for r in range(1, q):
        if abs(chi_q[r]) < 1e-10:
            continue
        # Find e in 0..d-1 with zeta_d^e = chi(r)
        for e in range(d):
            zd_e = np.exp(2j * np.pi * e / d)
            if abs(chi_q[r] - zd_e) < 1e-9:
                e_r_map[r] = e
                break
        else:
            raise RuntimeError(f"chi({r}) is not a d-th root of unity, d={d}")
    # chi(r) as int vec in Z^{phi_d} basis
    chi_r_vecs = {r: M[e_r_map[r]].copy() for r in e_r_map}

    # Sample x values
    sample = np.linspace(1, N, n_sample, dtype=np.int64)
    sample = np.unique(sample)

    n_fail = 0
    for x in sample:
        # Compute LHS = L_chi(x) as Z-vector in basis
        lhs = np.zeros(phi_d, dtype=np.int64)
        # Compute RHS = sum_r chi(r) * count_r(x) as Z-vector in basis
        rhs = np.zeros(phi_d, dtype=np.int64)
        for r, vec in chi_r_vecs.items():
            # S_r(x) = sum_{n<=x, n=r mod q} lambda(n)
            n_vals = np.arange(r, x + 1, q)
            S_r = int(lam[n_vals].sum())
            cnt_r = len(n_vals)
            lhs += S_r * vec
            rhs += cnt_r * vec
        diff = lhs - rhs
        if np.any(diff % 2 != 0):
            n_fail += 1
    return n_fail == 0, n_fail, max_resid


def free_identity_holds(chi_q: np.ndarray, actual_L_chi: np.ndarray,
                        N: int) -> tuple[bool, int, float]:
    """Check that actual_L_chi(x) - predicted L_chi(x) lies in 2 * Z[zeta_d]
    for every x in [1..N].

    Strategy.  Both actual and predicted live in Z[zeta_d] subset of C.
    The ring Z[zeta_d] is a free Z-module of rank phi(d).  We don't need
    its full basis: it suffices to check that 'actual - predicted' has
    rational coordinates in any chosen Z-basis of Z[zeta_d] all of whose
    parities are 0.

    Two cheap sufficient checks:
      (A) For order in {1, 2}: chi takes integer values, so L_chi(x) is
          an integer.  Check Re(actual - predicted) is an even integer.
      (B) For order in {4}: Z[zeta_4] = Z[i].  Check Re/Im differences
          are even integers.
      (C) For order in {3, 6}: Z[zeta_d] = Z[omega] where omega = zeta_3.
          Basis {1, omega}: a + b*omega has Re = a - b/2, Im = b*sqrt(3)/2.
          So a = Re + Im/sqrt(3), b = 2*Im/sqrt(3).  Compute these and
          verify both are even integers.
      (D) For order in {5, 10}: Z[zeta_5] is rank 4.  Use {1, zeta, zeta^2,
          zeta^3} basis.  Solve a 4x4 system per x -- expensive but only
          needs to be done sparsely (subsample).
      (E) For order in {12}: similar.

    Returns (holds_globally, num_x_failed, max_residual_after_proj).
    """
    diff = actual_L_chi - predict_chi_mod2(chi_q, N)
    # Determine character order (smallest d with chi^d = principal)
    d = 1
    cur = chi_q.copy()
    chi0 = np.zeros_like(chi_q)
    chi0[1:] = 1.0  # principal character on (Z/q)*
    while True:
        if np.allclose(cur, chi0, atol=1e-9):
            break
        cur = cur * chi_q
        d += 1
        if d > 60:
            raise RuntimeError("character order not found")

    # Coordinate extraction
    if d in (1, 2):
        # diff should be integer; check Re mod 2.
        a = np.round(diff.real).astype(np.int64)
        ok = bool(np.all(a % 2 == 0))
        n_fail = int(np.sum(a % 2 != 0))
        return ok, n_fail, float(np.max(np.abs(diff.imag)))
    if d == 4:
        a = np.round(diff.real).astype(np.int64)
        b = np.round(diff.imag).astype(np.int64)
        ok = bool(np.all((a % 2 == 0) & (b % 2 == 0)))
        n_fail = int(np.sum((a % 2 != 0) | (b % 2 != 0)))
        return ok, n_fail, max(float(np.max(np.abs(diff.real - a))),
                               float(np.max(np.abs(diff.imag - b))))
    if d in (3, 6):
        # Basis {1, omega} with omega = e^{2 pi i / 3}, Re=-1/2, Im=sqrt(3)/2.
        # diff = a + b*omega = (a - b/2) + b*sqrt(3)/2 * i.
        # So b = 2 * Im / sqrt(3), a = Re + b/2.
        b_est = 2.0 * diff.imag / math.sqrt(3.0)
        a_est = diff.real + 0.5 * b_est
        a = np.round(a_est).astype(np.int64)
        b = np.round(b_est).astype(np.int64)
        # Sanity: residual
        res_a = float(np.max(np.abs(a_est - a)))
        res_b = float(np.max(np.abs(b_est - b)))
        ok = bool(np.all((a % 2 == 0) & (b % 2 == 0)))
        n_fail = int(np.sum((a % 2 != 0) | (b % 2 != 0)))
        return ok, n_fail, max(res_a, res_b)
    # For d in {5, 10, 12} use subsampled basis-coord projection.
    # Pick a sample of x indices.
    sample = np.arange(0, N, max(1, N // 1000))
    rho = np.exp(2j * np.pi / d)
    basis = np.array([rho ** j for j in range(_phi_d(d))],
                      dtype=np.complex128)
    A = np.column_stack([basis.real, basis.imag])  # shape (phi_d, 2)
    # diff[x] = sum_j coeff_j * basis[j], so diff[x] (as 2-vector Re, Im)
    # = A.T @ coeffs.  Solve coeffs = pinv(A.T) @ [Re, Im].
    pinv = np.linalg.pinv(A.T)  # shape (phi_d, 2)
    n_fail = 0
    max_res = 0.0
    for x_idx in sample:
        v = np.array([diff[x_idx].real, diff[x_idx].imag])
        coeffs_real = pinv @ v
        coeffs = np.round(coeffs_real).astype(np.int64)
        res = float(np.max(np.abs(coeffs_real - coeffs)))
        max_res = max(max_res, res)
        if np.any(coeffs % 2 != 0):
            n_fail += 1
    return n_fail == 0, n_fail, max_res


def _phi_d(d: int) -> int:
    """Euler phi of d (used for Z[zeta_d] rank)."""
    cnt = 0
    for k in range(1, d + 1):
        if math.gcd(k, d) == 1:
            cnt += 1
    return cnt


# ---------- Main experiment -----------------------------------------------


def experiment_for_q(q: int, N: int, Omega: np.ndarray,
                     primes_mask: np.ndarray) -> Dict:
    """Run the full battery for prime q.  Returns a dict of summary numbers
    that the results-writer dumps to markdown."""
    print(f"\n========== q = {q} ==========")
    chars, log_g, g = char_table(q)
    print(f"  primitive root g = {g}, phi(q) = {q - 1}")

    # Liouville lambda
    lam = np.where(Omega % 2, -1, 1).astype(np.int64)
    lam[0] = 0

    # pi(x; q, a) for each a in (Z/q)^*
    pi_qa = {}
    for a in range(1, q):
        # Primes <= x with p ≡ a mod q
        flag = (primes_mask & (np.arange(len(Omega)) % q == a)).astype(np.int64)
        pi_qa[a] = np.cumsum(flag)

    summary = {
        "q": q,
        "phi_q": q - 1,
        "primitive_root": g,
        "characters": [],
    }

    for k, chi_q in enumerate(chars):
        # Lift chi(n) for n in 0..N
        chi_n = lift_chi_to_N(chi_q, len(Omega) - 1)
        # f(n) = lambda(n) * chi(n)
        f = lam * chi_n
        # Cumulative L_chi(x) = sum_{n<=x} f(n) (complex array indexed 0..N)
        L_full = np.cumsum(f)  # complex
        # Restricted to x = 1..N
        actual_L_chi = L_full[1:]
        is_real_chi = float(np.max(np.abs(chi_q.imag))) < 1e-12

        # Full free-identity test (handles any character order)
        ok_id, n_fail_id, res_id = free_identity_holds(
            chi_q, actual_L_chi, len(actual_L_chi)
        )
        # Exact integer-arithmetic check on a sample (definitive for orders >= 5)
        ok_id_exact, n_fail_exact, _ = free_identity_holds_exact(
            chi_q, lam, len(actual_L_chi), n_sample=2000
        )

        # For backward compat with prior printing code we keep "Re/Im match" labels.
        re_match = (len(actual_L_chi) - n_fail_id) if ok_id else (
            len(actual_L_chi) - n_fail_id
        )
        im_match = -1 if is_real_chi else re_match
        re_resid = res_id
        im_resid = res_id

        # "Next bit" - the parity of the count of +1 contributions of Re part.
        # For real chi: Re(f(n)) in {-1, 0, +1}; A_chi(x) = #{n<=x : Re(f(n))=+1}.
        # parity = (#nonzero + Re_L)/2 mod 2  iff Re takes only {-1,0,+1}.
        # We compute it directly from the streams.
        nonzero_re = (np.abs(f.real) > 0.5).astype(np.int64)
        # number of +1 contributions in the real part:
        pos_re = ((f.real > 0.5).astype(np.int64))
        A_re = np.cumsum(pos_re)
        A_re_par = (A_re[1:] % 2).astype(np.int8)

        # Pseudorandomness battery on A_re mod 2
        h8 = block_entropy(A_re_par, 8)
        ac = autocorr_max(A_re_par, [1, 2, 3, 5, 10, 30])
        # Use max FFT length for spectral line
        FFT_LEN = min(200_000, len(A_re_par))
        fz = fft_max_spectral_line(A_re_par[:FFT_LEN])
        LFSR_LEN = min(4096, len(A_re_par))
        lf = lfsr_length_gf2(A_re_par[:LFSR_LEN]) / LFSR_LEN

        # MI with pi(x; q, a) mod 2 for each a
        mi_list = {}
        for a in range(1, q):
            target = (pi_qa[a][1:] % 2).astype(np.int8)
            n_use = min(len(target), len(A_re_par))
            mi = mutual_info_bits(A_re_par[:n_use], target[:n_use])
            mi_list[a] = mi

        order_chi = (q - 1) // math.gcd(q - 1, k) if k > 0 else 1
        char_info = {
            "k": k,
            "order": order_chi,
            "is_principal": (k == 0),
            "is_real": is_real_chi,
            "re_resid_FP": re_resid,
            "im_resid_FP": im_resid,
            "free_id_holds_FP": ok_id,
            "free_id_n_fail_FP": n_fail_id,
            "free_id_residual_FP": res_id,
            "free_id_holds_exact": ok_id_exact,
            "free_id_n_fail_exact_of_sample": n_fail_exact,
            "free_id_re_match": re_match,
            "free_id_re_total": len(actual_L_chi),
            "free_id_im_match": im_match,
            "next_bit_h_block_8": h8,
            "next_bit_AC_max": ac,
            "next_bit_FFT_z": fz,
            "next_bit_LFSR_ratio": lf,
            "MI_with_pi_qa_mod2": mi_list,
        }
        summary["characters"].append(char_info)

        marker = "PRINCIPAL" if k == 0 else ("REAL" if is_real_chi else "COMPLEX")
        print(f"  chi[{k:>2d}]  ({marker:9s}, order={order_chi:>2d}):")
        print(f"     FREE IDENTITY in Z[zeta_d]/2:")
        print(f"       FP-projection check: holds={ok_id}, "
              f"n_fail={n_fail_id}/{len(actual_L_chi)}, residual={res_id:.2e}")
        print(f"       EXACT integer check: holds={ok_id_exact}, "
              f"n_fail={n_fail_exact}/2000 sampled x")
        print(f"     next-bit (A_chi mod 2):  h8={h8:.4f}  AC={ac:.4f}  "
              f"FFTz={fz:.2f}  LFSR/N={lf:.4f}")
        max_mi = max(mi_list.values()) if mi_list else 0.0
        print(f"     max MI with pi(x;q,a) mod 2 over a in (Z/q)*: {max_mi:.5e} bits")

    return summary


def main():
    N = 1_000_000
    print(f"Sieving Omega for n in [1, {N}] ...")
    t0 = time.time()
    Omega, spf = sieve_omega(N)
    print(f"  done in {time.time() - t0:.1f}s")
    print("Sieving prime mask ...")
    t0 = time.time()
    primes_mask = sieve_primes(N)
    print(f"  done in {time.time() - t0:.1f}s")

    all_summaries = []
    for q in (3, 5, 7, 11, 13):
        s = experiment_for_q(q, N, Omega, primes_mask)
        all_summaries.append(s)

    # Final compact dump
    print("\n========== COMPACT TABLE (max over chi mod q) ==========")
    print(f"{'q':>3s} | {'#chi':>5s} | {'all free-id?':>12s} | "
          f"{'max h8':>8s} | {'max AC':>8s} | {'max FFTz':>8s} | "
          f"{'min LFSR/N':>10s} | {'max MI bits':>12s}")
    for s in all_summaries:
        all_free_id = all(c["free_id_holds_exact"] for c in s["characters"])
        max_h8 = max(c["next_bit_h_block_8"] for c in s["characters"])
        max_ac = max(c["next_bit_AC_max"] for c in s["characters"])
        max_fz = max(c["next_bit_FFT_z"] for c in s["characters"])
        min_lf = min(c["next_bit_LFSR_ratio"] for c in s["characters"])
        max_mi = max(
            (max(c["MI_with_pi_qa_mod2"].values()) if c["MI_with_pi_qa_mod2"]
             else 0.0)
            for c in s["characters"]
        )
        print(f"{s['q']:>3d} | {s['phi_q']:>5d} | {str(all_free_id):>12s} | "
              f"{max_h8:8.4f} | {max_ac:8.4f} | {max_fz:8.2f} | "
              f"{min_lf:10.4f} | {max_mi:12.5e}")

    return all_summaries


if __name__ == "__main__":
    main()
