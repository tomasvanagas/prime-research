"""S247 — L^p hierarchy of the prime-indicator polynomial f_N(z) = Σ χ_P(n) z^n.

Tests the paradigm-shift composition prediction:

    ‖f_N‖_p^p / [Li(N)^p / N]  →  G_p · S_{p-1}^HL

for even integer p ∈ {4, 6, 8} where:
    G_p = (1/π) ∫ |sin u/u|^p du     (sinc geometric constant)
    S_{p-1}^HL = ∏_{prime r} (1 + 1/(r-1)^{p-1})  (HL singular series)

The same HL singular series governs both L^∞ (E2.21) and L²-spike (S168);
the conjecture is that all higher even L^p inherit it via the major-arc
spike-product-sinc-window structure.

Also tests:
    F3 (cumulative-q-strip closure): ratio of ‖f_N⊥V_q^prim_{q≤Q}‖_p^p
    to ‖f_N‖_p^p across Q ∈ {0, 2, 6, 30, 210}, predicted to match
    the partial Selberg-Delange sum over q sqfree ≤ Q.

    F4 (matched-density Bernoulli null): ‖f^B_N‖_p^p / [Li(N)^p / N]
    should NOT show HL singular-series modulation (Bernoulli has no
    arithmetic spikes).

Run modes:
    --test           N ∈ {2^12} only, fast smoke test.
    --full           N ∈ {2^12, 2^14, 2^16, 2^18}, all p, all strips.
    --F3-only        only F3 closure-fraction.
    --F4-only        only F4 Bernoulli null.

Output: lp_hierarchy_chi_p_results.json plus a printed summary table.
"""

import argparse
import json
import math
import time
from pathlib import Path

import numpy as np

# ----------------------------------------------------------------------
# 1. Closed-form constants (project-internal, no external lookup)
# ----------------------------------------------------------------------

# Borwein-Borwein closed forms for ∫|sin u / u|^p du / π.
# Computed once and verified against high-precision numerical integration.
G_FROM_TABLE = {
    2: 1.0,                # ∫ sinc^2 = π
    4: 2.0 / 3.0,          # ∫ sinc^4 = 2π/3
    6: 11.0 / 20.0,        # ∫ sinc^6 = 11π/20
    8: 151.0 / 315.0,      # ∫ sinc^8 = 151π/315
}


def primes_upto(M: int) -> np.ndarray:
    """Sieve of Eratosthenes, returns primes ≤ M as np.uint32."""
    if M < 2:
        return np.array([], dtype=np.uint32)
    bs = np.ones(M + 1, dtype=bool)
    bs[:2] = False
    for k in range(2, int(M**0.5) + 1):
        if bs[k]:
            bs[k * k :: k] = False
    return np.nonzero(bs)[0].astype(np.uint32)


def hl_singular_series(s: int, P_max: int = 10**6) -> float:
    """Compute S_s^HL = ∏_{prime r ≤ P_max} (1 + 1/(r-1)^s)."""
    primes = primes_upto(P_max)
    log_prod = 0.0
    for r in primes:
        log_prod += math.log1p(1.0 / (int(r) - 1) ** s)
    # Tail estimate: Σ_{r > P_max} 1/(r-1)^s ≈ ∫_{P_max} 1/(t-1)^s dt
    # = 1/((s-1) (P_max-1)^{s-1}) for s ≥ 2; negligible at P_max = 10^6.
    return math.exp(log_prod)


def sinc_integral_numerical(p: int) -> float:
    """G_p numerically via trapezoid rule on a fine grid (sanity check)."""
    # ∫_{-∞}^{∞} sinc^p du = ∫_{-A}^{A} (sin u/u)^p du for u≠0, p>1.
    A = 4000.0
    nu = 800001  # odd → grid lands on u=0
    u = np.linspace(-A, A, nu)
    du = u[1] - u[0]
    with np.errstate(invalid="ignore"):
        val = (np.sin(u) / u) ** p
    val[len(u) // 2] = 1.0  # sinc(0) = 1
    return float(np.sum(val) * du / math.pi)


# ----------------------------------------------------------------------
# 2. L^p norm of f_N via FFT
# ----------------------------------------------------------------------

def f_N_lp_norms(coeffs: np.ndarray, ps: list[int], M_fft: int) -> dict[int, float]:
    """Compute ‖f‖_p^p := (1/2π) ∫_0^{2π} |f(e^{iθ})|^p dθ via FFT.

    Uses an FFT length M_fft >= len(coeffs) to evaluate f on a fine
    angular grid, then approximates the integral by the trapezoid rule
    (= mean over the grid points × 2π, divided by 2π → just the mean).
    """
    M_fft = max(M_fft, len(coeffs))
    if M_fft != int(2 ** round(math.log2(M_fft))):
        M_fft = 1 << int(math.ceil(math.log2(M_fft)))

    padded = np.zeros(M_fft, dtype=np.float64)
    padded[: len(coeffs)] = coeffs
    F = np.fft.fft(padded)
    abs_F = np.abs(F)
    out = {}
    for p in ps:
        out[p] = float(np.mean(abs_F ** p))
    return out


def chi_P_indicator(N: int) -> np.ndarray:
    """1-indexed array of length N+1 with arr[n] = 1 iff n prime, n ≤ N."""
    arr = np.zeros(N + 1, dtype=np.float64)
    primes = primes_upto(N)
    arr[primes] = 1.0
    return arr


# ----------------------------------------------------------------------
# 3. Strip V_q^prim subspaces (F3)
# ----------------------------------------------------------------------

def strip_Vq_prim_subspaces(coeffs: np.ndarray, Q: int) -> np.ndarray:
    """Project coeffs orthogonally onto complement of ⊕_{q sqf ≤ Q} V_q^prim.

    V_q^prim is the span of e^{2πi a n / q} for (a, q) = 1, embedded as
    column vectors of length N+1. The projection uses the (additive)
    discrete Fourier basis.

    Operationally: replace `coeffs` with `coeffs - Σ_{q sqf ≤ Q} coeffs_q`
    where coeffs_q is the orthogonal projection of `coeffs` onto V_q^prim.

    This is implemented in O(N · |Q-additive-orbits|) time using the
    additive-Fourier representation: for each (a, q) with (a, q)=1 and
    q sqfree ≤ Q, project onto the 1-D subspace span(e^{2πi a n/q})_n.
    """
    N = len(coeffs) - 1
    out = coeffs.copy()
    n = np.arange(N + 1, dtype=np.float64)

    # Enumerate squarefree integers q in [2, Q].
    sqf_qs = []
    for q in range(2, Q + 1):
        m = q
        is_sqf = True
        for r in (2, 3, 5, 7, 11, 13, 17, 19, 23):
            if r * r > m:
                break
            if m % (r * r) == 0:
                is_sqf = False
                break
        if is_sqf:
            # also need to check higher prime squares
            r = 29
            while r * r <= q:
                if q % (r * r) == 0:
                    is_sqf = False
                    break
                r += 2
            if is_sqf:
                sqf_qs.append(q)

    # Project off span{e^{2πi a n / q}} for each (a, q) with q sqf in [2, Q],
    # gcd(a, q) = 1. We pair a with q-a (complex conjugate) and project onto
    # the 2D real subspace span{cos, sin} of the corresponding character.
    # For q=2 the subspace collapses to 1D (a=1, the (-1)^n vector).
    for q in sqf_qs:
        if q == 2:
            v = (-1.0) ** n
            c = float(out @ v) / float(v @ v)
            out = out - c * v
            continue
        a_max = q // 2
        for a in range(1, a_max + 1):
            if math.gcd(a, q) != 1:
                continue
            phase = 2.0 * math.pi * a / q
            v_re = np.cos(phase * n)
            v_im = np.sin(phase * n)
            r_norm2 = float(v_re @ v_re)
            i_norm2 = float(v_im @ v_im)
            ri_inner = float(v_re @ v_im)
            det = r_norm2 * i_norm2 - ri_inner ** 2
            if abs(det) < 1e-12 * max(r_norm2, i_norm2, 1.0):
                continue
            ov_re = float(out @ v_re)
            ov_im = float(out @ v_im)
            c_re = (i_norm2 * ov_re - ri_inner * ov_im) / det
            c_im = (-ri_inner * ov_re + r_norm2 * ov_im) / det
            out = out - c_re * v_re - c_im * v_im
    return out


# ----------------------------------------------------------------------
# 4. Closed-form prediction
# ----------------------------------------------------------------------

def predict_norm_p_via_pi(pi_N: int, N: int, p: int) -> float:
    """Closed-form ‖f_N‖_p^p ≈ G_p · S_{p-1} · π(N)^p / N (paradigm-shift conjecture, refined to use π(N) explicitly).

    The G_p constant arises from the dirichlet-kernel L^p norm:
        ‖D_N‖_p^p / N^{p-1} → G_p = (1/π) ∫|sinc|^p du.

    For each sqfree q, the major arc near a/q with (a, q)=1 contributes
    (μ²(q)/φ(q)^p) · π(N)^p · |D_N(δ)|^p · weight_q(N), and summing over
    a's gives a φ(q) factor → φ(q) · μ²(q) / φ(q)^p = μ²(q)/φ(q)^{p-1}.
    Hence:
        ‖f_N‖_p^p ≈ G_p · π(N)^p / N · Σ_{q sqf} μ²(q) / φ(q)^{p-1}
                  = G_p · π(N)^p / N · ∏_{prime r} (1 + 1/(r-1)^{p-1})
    """
    G_p = G_FROM_TABLE[p]
    S = hl_singular_series(p - 1)
    return G_p * S * (pi_N ** p) / N


def pi_of(N: int) -> int:
    return int(np.sum(chi_P_indicator(N) > 0))


# ----------------------------------------------------------------------
# 5. Main driver
# ----------------------------------------------------------------------

def run(Ns: list[int], ps: list[int], strip_Qs: list[int], with_bern: bool, M_fft_pad: int = 8) -> dict:
    """Run the L^p hierarchy experiment for prime indicator and (optionally) Bernoulli null."""
    results: dict = {
        "ps": ps,
        "G_p": {p: G_FROM_TABLE[p] for p in ps},
        "G_p_numerical_check": {p: sinc_integral_numerical(p) for p in ps if p <= 8},
        "S_HL_p_minus_1": {p: hl_singular_series(p - 1) for p in ps},
        "Ns": Ns,
        "pi_of_N": {},
        "strip_Qs": strip_Qs,
        "predictions": {},
        "empirical_chi_P": {},
        "ratio_chi_P": {},
        "empirical_chi_P_stripped": {},
        "ratio_chi_P_stripped": {},
        "empirical_bern": {},
        "ratio_bern": {},
    }
    pi_lookup: dict[int, int] = {}
    for N in Ns:
        pi_lookup[N] = int(chi_P_indicator(N).sum())
        results["pi_of_N"][N] = pi_lookup[N]
    for p in ps:
        for N in Ns:
            key = f"p={p}_N={N}"
            results["predictions"][key] = predict_norm_p_via_pi(pi_lookup[N], N, p)

    rng = np.random.default_rng(2026)
    for N in Ns:
        t0 = time.time()
        coeffs = chi_P_indicator(N)
        pi_N = int(coeffs.sum())
        # Use π(N) directly in the denominator: ratio = ‖f‖_p^p · N / π(N)^p.
        # This removes the Li(N)/π(N) finite-N correction at p ≥ 4.
        norm_factor_per_p = {p: (pi_N ** p) / N for p in ps}

        M_fft = M_fft_pad * (N + 1)
        norms = f_N_lp_norms(coeffs, ps, M_fft)
        for p in ps:
            key = f"p={p}_N={N}"
            results["empirical_chi_P"][key] = norms[p]
            results["ratio_chi_P"][key] = norms[p] / norm_factor_per_p[p]

        # Strips
        for Q in strip_Qs:
            if Q == 0:
                stripped = coeffs
            else:
                stripped = strip_Vq_prim_subspaces(coeffs, Q)
            norms_strip = f_N_lp_norms(stripped, ps, M_fft)
            for p in ps:
                key = f"p={p}_N={N}_Q={Q}"
                results["empirical_chi_P_stripped"][key] = norms_strip[p]
                results["ratio_chi_P_stripped"][key] = norms_strip[p] / norm_factor_per_p[p]

        # Bernoulli null
        if with_bern:
            p_density = pi_N / N
            bern_arr = np.zeros(N + 1, dtype=np.float64)
            # B1: i.i.d. Bernoulli(p_density) on [2, N], 0 at index 0, 1.
            bern_arr[2:] = (rng.random(N - 1) < p_density).astype(np.float64)
            norms_bern = f_N_lp_norms(bern_arr, ps, M_fft)
            for p in ps:
                key = f"p={p}_N={N}"
                results["empirical_bern"][key] = norms_bern[p]
                results["ratio_bern"][key] = norms_bern[p] / norm_factor_per_p[p]

        elapsed = time.time() - t0
        print(f"  N = {N:>10d}  pi(N) = {pi_N}  elapsed = {elapsed:.1f}s")

    return results


def print_summary(results: dict, Ns: list[int], ps: list[int]) -> None:
    G_p = results["G_p"]
    S_p = results["S_HL_p_minus_1"]
    print("\n=== L^p hierarchy of f_N — summary ===\n")
    print(f"{'p':>3} {'G_p':>10} {'G_p (num)':>12} {'S_{p-1}':>10} {'predicted ratio':>18}")
    for p in ps:
        gp = G_p[p]
        gp_num = results["G_p_numerical_check"].get(p, float("nan"))
        sp = S_p[p]
        print(f"{p:>3} {gp:>10.6f} {gp_num:>12.6f} {sp:>10.6f} {gp*sp:>18.6f}")

    print("\n--- Empirical chi_P (raw) ---")
    print(f"{'N':>10} | " + " | ".join(f"p={p:<2} ratio" for p in ps))
    for N in Ns:
        cells = []
        for p in ps:
            key = f"p={p}_N={N}"
            r = results["ratio_chi_P"][key]
            cells.append(f"{r:>10.4f}")
        print(f"{N:>10} | " + " | ".join(cells))

    if results.get("ratio_bern"):
        print("\n--- Empirical Bernoulli (matched density) ---")
        print(f"{'N':>10} | " + " | ".join(f"p={p:<2} ratio" for p in ps))
        for N in Ns:
            cells = []
            for p in ps:
                key = f"p={p}_N={N}"
                r = results["ratio_bern"].get(key, float("nan"))
                cells.append(f"{r:>10.4f}")
            print(f"{N:>10} | " + " | ".join(cells))

    if results.get("ratio_chi_P_stripped"):
        print("\n--- chi_P stripped V_q^prim (sqfree q ≤ Q) ---")
        print("(p=4 only)")
        Qs = sorted({int(k.split("Q=")[1]) for k in results["ratio_chi_P_stripped"].keys() if "p=4_" in k})
        print(f"{'N':>10} | " + " | ".join(f"Q={Q:>3}" for Q in Qs))
        for N in Ns:
            cells = []
            for Q in Qs:
                key = f"p=4_N={N}_Q={Q}"
                r = results["ratio_chi_P_stripped"].get(key, float("nan"))
                cells.append(f"{r:>10.4f}")
            print(f"{N:>10} | " + " | ".join(cells))

    print("\n--- Closed-form prediction (paradigm-shift conjecture) ---")
    print("Ratio should converge to G_p · S_{p-1} as N → ∞.")
    print("If empirical raw ratio ≠ predicted ratio at all N, identify the prefactor.")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--mode", choices=["test", "small", "full", "convergence", "F3-only", "F4-only"], default="small")
    ap.add_argument("--out", type=Path, default=Path(__file__).parent / "lp_hierarchy_chi_p_results.json")
    args = ap.parse_args()

    if args.mode == "test":
        Ns = [1 << 12]
        ps = [2, 4, 6, 8]
        strip_Qs = [0, 2]
        with_bern = True
    elif args.mode == "small":
        Ns = [1 << 12, 1 << 14]
        ps = [2, 4, 6, 8]
        strip_Qs = [0, 2, 6, 30]
        with_bern = True
    elif args.mode == "full":
        # Two-phase: large N for ratio convergence (no strip);
        # smaller N up to 2^16 for the cumulative-Q closure.
        Ns = [1 << 12, 1 << 14, 1 << 16]
        ps = [2, 4, 6, 8]
        strip_Qs = [0, 2, 6, 30, 210]
        with_bern = True
    elif args.mode == "convergence":
        # Larger N, no strip (so each N runs in seconds).
        Ns = [1 << 14, 1 << 16, 1 << 18, 1 << 20]
        ps = [2, 4, 6, 8]
        strip_Qs = []
        with_bern = True
    elif args.mode == "F3-only":
        Ns = [1 << 14]
        ps = [4]
        strip_Qs = [0, 2, 6, 30, 210]
        with_bern = False
    elif args.mode == "F4-only":
        Ns = [1 << 14, 1 << 16]
        ps = [2, 4, 6, 8]
        strip_Qs = []
        with_bern = True

    print("Running L^p hierarchy experiment.")
    print(f"  Ns: {Ns}\n  ps: {ps}\n  strip_Qs: {strip_Qs}\n  with_bern: {with_bern}")
    t0 = time.time()
    results = run(Ns, ps, strip_Qs, with_bern)
    elapsed = time.time() - t0
    print(f"\nTotal elapsed: {elapsed:.1f}s")

    args.out.write_text(json.dumps(results, indent=2))
    print(f"Wrote {args.out}")
    print_summary(results, Ns, ps)


if __name__ == "__main__":
    main()
