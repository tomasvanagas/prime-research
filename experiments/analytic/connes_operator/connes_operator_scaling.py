"""
FOCUS-1: Connes-Consani-Moscovici (arXiv:2511.22755) operator scaling experiment.

Question: For the rank-one perturbation of the scaling operator on L^2([lambda^{-1}, lambda]),
how does the number K of Riemann zeros recovered to within rounding precision (err < 0.5)
scale with the prime budget B?

If K_accurate(B) ~ polylog(B): Chain A is a real polylog architecture (HUGE).
If K_accurate(B) ~ B: Chain A closes via Equivalence (rank-one has B params, fits ~B eigenvalues).

CONSTRUCTION
------------
Use the logarithmic coordinate u = log(t), t in [lambda^{-1}, lambda] => u in [-L, L], L = log(lambda).
The scaling operator -i t d/dt becomes D = -i d/du on L^2([-L, L]).
Discretize on Fourier basis with antiperiodic-style symmetric grid:
    phi_n(u) = (1/sqrt(2L)) exp(i pi n u / L)  for n = -N..N
=> D is diagonal with entries pi n / L.

Perturbation: rank-one V = c |v><v|, with v built from primes <= B:
    v(u) = sum_{p <= B} (log p) / sqrt(p) * delta_supported_kernel(u - log p)
faithful continuous form would use the prime-counting Mellin kernel; we implement the
discretized projection v_n = <phi_n, v> = sum_{p<=B} (log p)/sqrt(p) * phi_n(log p)*  conj.
We choose the Hermitian rank-one V_{mn} = c (v_m conj(v_n) + v_n conj(v_m))/2 -> standard
self-adjoint c |v><v| since |v><v| is already self-adjoint. We tune the coupling c.

BASELINE TO COMPARE
-------------------
True non-trivial Riemann zeta zero ordinates gamma_K from data/zeta_zeros_2000.txt.
We compare the eigenvalue mu_K closest to gamma_K via a stable monotone matching
on the positive eigenvalue branch.

OUTPUT
------
- err(K, B) for K = 1..50, B = 1..9 primes
- K_accurate(B) := largest K with err(K, B) < 0.5
- Fits: K_accurate(B) vs B linear / power law / log
"""
from __future__ import annotations
import os
import sys
import math
import numpy as np
from numpy.linalg import eigvalsh

REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", ".."))
ZEROS_FILE = os.path.join(REPO, "data", "zeta_zeros_2000.txt")

PRIMES = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31]


def load_true_zeros(path: str = ZEROS_FILE, K: int = 60) -> np.ndarray:
    out = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            out.append(float(line))
            if len(out) >= K:
                break
    return np.array(out)


def build_operator(primes: list[int], N: int, L: float, coupling: float,
                   variant: str = "delta_comb", use_prime_powers: bool = True) -> np.ndarray:
    """
    Discretize the scaling operator D = -i d/du on [-L, L] in Fourier basis phi_n(u) = e^{i pi n u/L}/sqrt(2L),
    n = -N..N (size 2N+1). D-eigenvalues d_n = pi n / L.

    Self-adjoint rank-one perturbation V = c |v><v|.

    Variants for v (all encode primes <= max(primes)):
      - "delta_comb":   v(u) = sum_{p^k <= P_max} (log p)/sqrt(p^k) * delta(u - log p^k)
                        i.e. the von-Mangoldt-weighted comb on log space.
      - "exp_kernel":   v(u) = sum_{p} (log p)/sqrt(p) * 1[u <= log p] / sqrt(p)   -- Mertens pairing
      - "psi_step":     v(u) = sum_{p^k <= e^u} (log p) -- Chebyshev psi cumulative;
                        this is the natural functional whose Mellin transform is -zeta'/zeta.
    """
    sz = 2 * N + 1
    ns = np.arange(-N, N + 1, dtype=np.float64)
    d = (np.pi / L) * ns

    # Build prime-power list with weights (von Mangoldt Lambda(n) is the natural weight)
    if use_prime_powers:
        items = []  # (n_value, weight log p, position log n)
        for p in primes:
            pk = p
            while math.log(pk) <= L:
                items.append((pk, math.log(p)))
                pk *= p
    else:
        items = [(p, math.log(p)) for p in primes if math.log(p) <= L]

    v = np.zeros(sz, dtype=np.complex128)

    if variant == "delta_comb":
        for (npp, wt) in items:
            u_pos = math.log(npp)
            # Fourier coeff of delta(u - u_pos):
            v += (wt / math.sqrt(npp)) * (1.0 / math.sqrt(2 * L)) * np.exp(-1j * np.pi * ns * u_pos / L)

    elif variant == "psi_step":
        # v is sum of indicator functions [log p^k, L] weighted by log p / sqrt(p^k)
        # Fourier coeff of 1[a, L]: int_a^L e^{-i pi n u/L}/sqrt(2L) du
        # n != 0: (-L / (i pi n)) * (e^{-i pi n} - e^{-i pi n a/L}) / sqrt(2L)
        # n == 0: (L - a) / sqrt(2L)
        zero_mask = (ns == 0)
        nz = ns.copy()
        nz[zero_mask] = 1.0  # avoid div by zero
        for (npp, wt) in items:
            a = math.log(npp)
            factor = -L / (1j * np.pi * nz)
            coeff = factor * (np.exp(-1j * np.pi * nz) - np.exp(-1j * np.pi * nz * a / L)) / math.sqrt(2 * L)
            coeff[zero_mask] = (L - a) / math.sqrt(2 * L)
            v += (wt / math.sqrt(npp)) * coeff

    elif variant == "exp_kernel":
        # v(u) = sum_p (log p)/sqrt(p) * exp(-(u - log p)^2 / sigma^2) (Gaussian-localized)
        # use sigma = pi/(2 max gamma) ~ small; here just a moderate localizer
        sigma = 0.1
        # Fourier coeff: integral over [-L, L] of Gaussian centered at log p
        u_grid = np.linspace(-L, L, max(2 * N + 1, 4001))
        du = u_grid[1] - u_grid[0]
        v_real = np.zeros_like(u_grid)
        for (npp, wt) in items:
            a = math.log(npp)
            v_real += (wt / math.sqrt(npp)) * np.exp(-((u_grid - a) ** 2) / sigma ** 2)
        # Project onto Fourier basis
        for k_idx, n in enumerate(ns):
            phi = np.exp(-1j * np.pi * n * u_grid / L) / math.sqrt(2 * L)
            v[k_idx] = np.sum(v_real * phi) * du

    else:
        raise ValueError(f"unknown variant {variant}")

    # H = diag(d) + c * v v^H
    H = np.diag(d).astype(np.complex128)
    H += coupling * np.outer(v, np.conjugate(v))
    H = 0.5 * (H + H.conj().T)
    return H


def positive_eigs(H: np.ndarray) -> np.ndarray:
    e = eigvalsh(H)
    e = np.real(e)
    e.sort()
    return e[e > 1e-9]


def measure_error(eigs_pos: np.ndarray, true_zeros: np.ndarray) -> np.ndarray:
    """For each true zero gamma_k, find the closest eigenvalue (greedy, with no reuse)."""
    used = np.zeros(len(eigs_pos), dtype=bool)
    errs = np.zeros(len(true_zeros))
    for k, gamma in enumerate(true_zeros):
        diffs = np.abs(eigs_pos - gamma)
        diffs[used] = np.inf
        j = int(np.argmin(diffs))
        errs[k] = float(diffs[j])
        used[j] = True
    return errs


def measure_error_monotone(eigs_pos: np.ndarray, true_zeros: np.ndarray) -> np.ndarray:
    """STRICT test: eigenvalue ranked K-th must match true zero K-th in order.

    This is the correct test for spectral realization (CCM's claim is that the
    K-th eigenvalue equals gamma_K, not that some eigenvalue is near gamma_K).
    """
    K = len(true_zeros)
    if len(eigs_pos) < K:
        # pad with infinity
        e = np.concatenate([eigs_pos, np.full(K - len(eigs_pos), np.inf)])
    else:
        e = eigs_pos[:K]
    return np.abs(e - true_zeros)


def k_accurate(errs: np.ndarray, threshold: float = 0.5) -> int:
    """Largest K such that err[0..K-1] < threshold (contiguous from start)."""
    K = 0
    for e in errs:
        if e < threshold:
            K += 1
        else:
            break
    return K


def fit_law(B_arr: np.ndarray, K_arr: np.ndarray) -> dict:
    out = {}
    mask = K_arr > 0
    if mask.sum() < 2:
        return {"linear": None, "power": None, "log": None}
    Bm = B_arr[mask].astype(float)
    Km = K_arr[mask].astype(float)
    # K = a * B + b
    a, b = np.polyfit(Bm, Km, 1)
    pred = a * Bm + b
    ss_res = ((Km - pred) ** 2).sum()
    ss_tot = ((Km - Km.mean()) ** 2).sum() if Km.var() > 0 else 1.0
    out["linear"] = {"slope": float(a), "intercept": float(b),
                     "R2": float(1 - ss_res / ss_tot) if ss_tot > 0 else 1.0}
    # K = a * B^alpha
    if (Bm > 0).all() and (Km > 0).all():
        la, lb = np.polyfit(np.log(Bm), np.log(Km), 1)
        out["power"] = {"alpha": float(la), "C": float(math.exp(lb))}
    # K = a + b log B
    if (Bm > 0).all():
        a2, b2 = np.polyfit(np.log(Bm), Km, 1)
        pred2 = a2 * np.log(Bm) + b2
        ss_res2 = ((Km - pred2) ** 2).sum()
        out["log"] = {"slope": float(a2), "intercept": float(b2),
                      "R2": float(1 - ss_res2 / ss_tot) if ss_tot > 0 else 1.0}
    return out


def run_one(B_idx: int, N: int, L: float, coupling: float,
            true_zeros: np.ndarray, variant: str = "delta_comb",
            monotone: bool = True) -> tuple[np.ndarray, np.ndarray]:
    primes = PRIMES[:B_idx]
    H = build_operator(primes, N, L, coupling, variant=variant)
    eigs_pos = positive_eigs(H)
    if monotone:
        errs = measure_error_monotone(eigs_pos, true_zeros)
    else:
        errs = measure_error(eigs_pos, true_zeros)
    return eigs_pos, errs


def main():
    # Cutoff x = lambda^2 = max prime budget * safety; CCM uses x just above max prime.
    # We pick lambda generously: lambda^2 = 50 (so lambda ~ 7.07, log lambda ~ 1.96).
    # Larger L pushes finer eigenvalue spacing and lets us probe more zeros.
    # We pick L so that the unperturbed scaling-operator eigenvalues pi*n/L cover the
    # range of the first K_TRUE Riemann zeros [14, ~270]. With pi/L spacing ~ 0.5,
    # we need N ~ 540 to reach gamma_80 = 270.
    LAMBDA2 = 10000.0
    L = 0.5 * math.log(LAMBDA2)  # ~ 4.6; spacing pi/L ~ 0.68
    N = 1200  # mode count: 2401-dim matrix; max diagonal eig ~ 818
    K_TRUE = 50  # CCM published numerics target
    true_zeros = load_true_zeros(K=K_TRUE)
    print(f"L = {L:.4f}  N = {N}  matrix size = {2*N+1}")
    print(f"Eigenvalue spacing pi/L = {math.pi/L:.4f}")
    print(f"max diagonal eigenvalue = pi N / L = {math.pi*N/L:.2f}")

    # Variant + coupling tuning at B=6 against CCM's anchor of K=50.
    # We use 'matched' (CCM's view: nearest-neighbour matching) for tuning,
    # but report both 'matched' and 'monotone' (rank-K to gamma_K) for honesty.
    variants = ["delta_comb", "psi_step"]
    couplings = [-4.0, -2.0, -1.0, -0.5, 0.5, 1.0, 2.0, 4.0]
    best = None  # (variant, coupling, median_err_matched)
    print("\nTuning at B=6 (CCM anchor: K=50, errors 1e-55 to 1e-3):")
    for variant in variants:
        for c in couplings:
            try:
                eigs_pos, errs = run_one(6, N, L, c, true_zeros, variant=variant, monotone=False)
            except Exception as exc:  # noqa: BLE001
                print(f"  {variant:12s}  c = {c:6.2f}  ERROR: {exc}")
                continue
            median_err = float(np.median(errs[:20]))
            mean_err = float(np.mean(errs[:20]))
            print(f"  {variant:12s}  c = {c:6.2f}  median(matched err first 20) = {median_err:.4g}  mean = {mean_err:.4g}")
            if best is None or median_err < best[2]:
                best = (variant, c, median_err)
    variant, coupling, _ = best
    print(f"=> using variant = {variant}, coupling = {coupling}")

    # CONTROL: the unperturbed scaling operator (B=0, no primes).
    # Eigenvalues are pi*n/L on a uniform grid.
    eigs0 = (np.pi / L) * np.arange(1, N + 1, dtype=float)
    errs0_mono = measure_error_monotone(eigs0, true_zeros)
    errs0_match = measure_error(eigs0, true_zeros)
    print("\n=== CONTROL: B=0 (no prime perturbation) ===")
    print(f"  monotone:   K_accurate(<0.5) = {k_accurate(errs0_mono, 0.5)}, "
          f"K(<0.1) = {k_accurate(errs0_mono, 0.1)}, "
          f"K(<0.01) = {k_accurate(errs0_mono, 0.01)}")
    print(f"  matched:    K_accurate(<0.5) = {k_accurate(errs0_match, 0.5)}, "
          f"K(<0.1) = {k_accurate(errs0_match, 0.1)}, "
          f"K(<0.01) = {k_accurate(errs0_match, 0.01)}")
    print(f"  err first 5 (monotone): {errs0_mono[:5]}")

    # Sweep B = 1..9 primes
    B_max = 9
    K_arr = np.zeros(B_max, dtype=int)
    err_table = np.zeros((B_max, K_TRUE))
    K_arr_05 = np.zeros(B_max, dtype=int)
    K_arr_01 = np.zeros(B_max, dtype=int)
    K_arr_001 = np.zeros(B_max, dtype=int)
    median_err_arr = np.zeros(B_max)
    print("\n=== MATCHED MATCHING (greedy nearest, CCM-style) ===")
    for b_idx in range(1, B_max + 1):
        eigs_pos, errs = run_one(b_idx, N, L, coupling, true_zeros, variant=variant, monotone=False)
        err_table[b_idx - 1, :] = errs
        K_arr[b_idx - 1] = k_accurate(errs, 0.5)
        K_arr_05[b_idx - 1] = k_accurate(errs, 0.5)
        K_arr_01[b_idx - 1] = k_accurate(errs, 0.1)
        K_arr_001[b_idx - 1] = k_accurate(errs, 0.01)
        primes_used = PRIMES[:b_idx]
        median = float(np.median(errs[:50]))
        median_err_arr[b_idx - 1] = median
        print(f"B={b_idx} (p<={primes_used[-1]:2d}) "
              f"K(<0.5)={K_arr_05[b_idx-1]:3d}  K(<0.1)={K_arr_01[b_idx-1]:3d}  K(<0.01)={K_arr_001[b_idx-1]:3d}  "
              f"med_err(50)={median:.4g}  err[1]={errs[0]:.3g}  err[10]={errs[9]:.3g}  err[50]={errs[49]:.3g}")

    # Also report monotone matching (mu_K vs gamma_K, the architecturally correct test)
    print("\n=== MONOTONE MATCHING (mu_K vs gamma_K, structural test) ===")
    K_arr_mono = np.zeros(B_max, dtype=int)
    median_err_mono = np.zeros(B_max)
    for b_idx in range(1, B_max + 1):
        eigs_pos, errs = run_one(b_idx, N, L, coupling, true_zeros, variant=variant, monotone=True)
        K_arr_mono[b_idx - 1] = k_accurate(errs, 0.5)
        median_err_mono[b_idx - 1] = float(np.median(errs[:50]))
        primes_used = PRIMES[:b_idx]
        print(f"B={b_idx} (p<={primes_used[-1]:2d}) K(<0.5)={K_arr_mono[b_idx-1]:3d}  "
              f"med_err(50)={median_err_mono[b_idx-1]:.4g}  err[1]={errs[0]:.3g}  err[10]={errs[9]:.3g}")

    # Fit
    B_arr = np.arange(1, B_max + 1)
    fit = fit_law(B_arr, K_arr)
    print("\n=== FIT RESULTS ===")
    print(f"Linear:  K = {fit['linear']['slope']:.3f} * B + {fit['linear']['intercept']:.3f}, R^2 = {fit['linear']['R2']:.4f}")
    if fit.get("power"):
        print(f"Power:   K = {fit['power']['C']:.3f} * B^{fit['power']['alpha']:.3f}")
    if fit.get("log"):
        print(f"Log:     K = {fit['log']['slope']:.3f} * log(B) + {fit['log']['intercept']:.3f}, R^2 = {fit['log']['R2']:.4f}")

    # Extrapolation
    target_K_at_x_1e100 = math.sqrt(1e100)
    # number of primes up to B equals B by index here (we increment by adding next prime)
    # Convert: if K = a*B + b, to reach K=1e50 we need B ~ 1e50 / a primes
    if fit["linear"]:
        B_needed_linear = (target_K_at_x_1e100 - fit["linear"]["intercept"]) / max(fit["linear"]["slope"], 1e-9)
    else:
        B_needed_linear = float("inf")
    if fit.get("log") and fit["log"]["slope"] > 0:
        try:
            B_needed_log = math.exp((target_K_at_x_1e100 - fit["log"]["intercept"]) / max(fit["log"]["slope"], 1e-9))
        except OverflowError:
            B_needed_log = float("inf")
    else:
        B_needed_log = float("inf")

    # Save numbers
    out_path = os.path.join(os.path.dirname(__file__), "connes_operator_scaling_data.csv")
    with open(out_path, "w") as f:
        f.write("B_idx,prime_max,K_accurate," + ",".join(f"err{k+1}" for k in range(K_TRUE)) + "\n")
        for i in range(B_max):
            f.write(f"{i+1},{PRIMES[i]},{K_arr[i]},")
            f.write(",".join(f"{e:.6g}" for e in err_table[i]))
            f.write("\n")
    print(f"\nWrote: {out_path}")

    # Print verdict-relevant numbers for downstream MD
    print("\n=== VERDICT INPUTS ===")
    print(f"K_accurate(B=2) [<0.5 matched]   = {K_arr_05[1]}")
    print(f"K_accurate(B=6) [<0.5 matched]   = {K_arr_05[5]}  (CCM anchor; their construction reaches K=50)")
    print(f"K_accurate(B=9) [<0.5 matched]   = {K_arr_05[8]}")
    print(f"K_accurate(B=2..9) [<0.1 matched]: {list(K_arr_01)}")
    print(f"K_accurate(B=2..9) [<0.01 matched]: {list(K_arr_001)}")
    print(f"K_accurate(B=2..9) [<0.5 monotone]: {list(K_arr_mono)}")
    print(f"Median(err[:50]) at B=1..9 [matched]:  {[f'{x:.3f}' for x in median_err_arr]}")
    print(f"Median(err[:50]) at B=1..9 [monotone]: {[f'{x:.3f}' for x in median_err_mono]}")
    print(f"\nLinear slope (matched <0.5): {fit['linear']['slope']:.4f}  (zeros recovered per added prime)")
    print(f"To reach K = sqrt(1e100) = 1e50 zeros, linear fit needs B ~ {B_needed_linear:.3g} primes")
    print(f"Log fit needs B ~ {B_needed_log:.3g} primes")

    # CCM fidelity check
    print(f"\nCCM fidelity at B=6: their err[1]=2.5e-55, err[50]~1e-3.")
    print(f"  Our proxy err[1]={err_table[5][0]:.3g}, err[50]={err_table[5][49]:.3g}")
    print(f"  => Proxy does NOT reproduce CCM's quantitative accuracy (off by ~50 orders of magnitude on err[1]).")


if __name__ == "__main__":
    main()
