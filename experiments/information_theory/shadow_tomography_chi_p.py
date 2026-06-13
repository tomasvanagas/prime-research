"""
D11 — Classical-shadow query complexity of pi(x).

Cross-domain ingredient: classical shadow protocol of Huang-Kueng-Preskill 2020
(Nature Physics 16, 1050; arXiv:2002.08953), adapted to a CLASSICAL random-
Rademacher-shadow query model on the prime-indicator signal chi_P in {0, 1}^N.

Question (ATTACK_VECTORS.md D11):
    Under the random-shadow query model -- each query is y_j = <sigma_j, chi_P>
    for sigma_j ~ Uniform({-1, +1}^N) -- how many queries K suffice to estimate
    pi(M) for ALL M in [1, N] simultaneously to within +/- 1?

Theoretical prediction (this file's Theorem 1, derived from random-Rademacher
classical-shadow variance calculation, structurally identical to HKP shadow-
norm bound for global observables):

    pi_hat(M; K) = (1/K) sum_{j=1..K} y_j * <sigma_j, 1_{[1, M]}>

    E[pi_hat(M; K)] = pi(M)                            # unbiased
    Var[pi_hat(M; K)] = (M*pi(N) - pi(M) + pi(M)^2 - pi(M)) / K
                     ~~ M * pi(N) / K                  # leading order

    For epsilon = 1 simultaneous accuracy at all M up to N:
        K = Omega(N * pi(N))                            # NOT polylog

This experiment empirically verifies the M*pi(N)/K scaling and confirms the
closure: shadow protocol on global cumulative-window observables is no better
than direct sample-mean for whole-profile estimation. K_polylog is impossible.

The closure adds EDGE: sample/query complexity of pi(x) under random-mask
oracle is poly(N), distinct from time-complexity (sieve / Galway / Aggarwal)
and information-theoretic (E1.5) bounds.
"""

import numpy as np
from numpy.typing import NDArray
import time
import json
from pathlib import Path

RNG = np.random.default_rng(42)


def sieve(N: int) -> NDArray[np.int8]:
    """Return chi_P as int8 array indexed 0..N-1; chi_P[n-1] = 1 iff n is prime."""
    sieve_arr = np.ones(N + 1, dtype=bool)
    sieve_arr[0] = sieve_arr[1] = False
    for p in range(2, int(N**0.5) + 1):
        if sieve_arr[p]:
            sieve_arr[p * p :: p] = False
    return sieve_arr[1 : N + 1].astype(np.int8)  # chi_P at positions 1..N


def cumulative_pi(chi_p: NDArray[np.int8]) -> NDArray[np.int64]:
    """pi(M) for M = 1..N as cumulative sum."""
    return np.cumsum(chi_p, dtype=np.int64)


def rademacher_shadow_estimate(
    chi_p: NDArray[np.int8],
    M_targets: NDArray[np.int64],
    K: int,
    rng: np.random.Generator,
    block_size: int = 256,
) -> NDArray[np.float64]:
    """
    Compute pi_hat(M; K) = (1/K) sum_j y_j * s_j(M)
    where y_j = <sigma_j, chi_P>  (scalar)
          s_j(M) = <sigma_j, 1_{[1, M]}>  (scalar per M)
          sigma_j in {-1, +1}^N i.i.d. Rademacher.

    Returns array shape (len(M_targets),).
    """
    N = chi_p.size
    # chi_p_signed isn't needed; chi_P is {0,1} so y = sum sigma * chi_P over prime n.
    # For variance/efficiency we use chi_P float dtype.
    chi_pf = chi_p.astype(np.float64)

    # Accumulate
    pi_hat = np.zeros(len(M_targets), dtype=np.float64)

    # Process in blocks to bound memory: K x N = K * 2^20 * 1 byte = up to ~10^11 / block
    # We use block_size shadows per chunk.
    n_blocks = (K + block_size - 1) // block_size
    for b in range(n_blocks):
        kb = min(block_size, K - b * block_size)
        # Sigma block: shape (kb, N), values in {-1, +1}
        # int8 to save memory: kb * N bytes; for kb=256, N=2^20: 256MB (OK)
        sigma = (rng.integers(0, 2, size=(kb, N), dtype=np.int8) * 2 - 1).astype(np.int8)

        # y_j = <sigma_j, chi_P> for each j: shape (kb,)
        y = sigma.astype(np.float64) @ chi_pf  # (kb,)

        # For each M, compute s_j(M) = sum_{n <= M} sigma(n)
        # Using cumulative sum along axis=1 indexed at M_targets:
        sigma_cumsum = np.cumsum(sigma.astype(np.int32), axis=1)  # shape (kb, N)
        # s[j, t] = sigma_cumsum[j, M_t - 1] (since 1-indexed)
        s = sigma_cumsum[:, M_targets - 1].astype(np.float64)  # shape (kb, len(M_targets))

        # pi_hat contributions: sum_j y_j * s_j(M)
        pi_hat += y @ s  # broadcasts: (kb,) . (kb, len(M)) -> (len(M),)

        del sigma, sigma_cumsum, s
    pi_hat /= K
    return pi_hat


def trial_errors(
    chi_p: NDArray[np.int8],
    pi_true: NDArray[np.int64],
    M_targets: NDArray[np.int64],
    K: int,
    n_trials: int,
    rng: np.random.Generator,
    block_size: int = 256,
) -> dict:
    """Run n_trials independent shadow experiments at fixed K. Returns variance stats."""
    errors_per_trial: list[NDArray[np.float64]] = []
    t0 = time.time()
    for t in range(n_trials):
        pi_hat = rademacher_shadow_estimate(chi_p, M_targets, K, rng, block_size)
        err = pi_hat - pi_true[M_targets - 1]
        errors_per_trial.append(err)
    elapsed = time.time() - t0

    errs = np.array(errors_per_trial)  # (n_trials, len(M))
    return {
        "K": K,
        "n_trials": n_trials,
        "mean_err": errs.mean(axis=0).tolist(),  # bias estimate
        "std_err": errs.std(axis=0, ddof=1).tolist(),  # empirical std
        "L_inf_per_trial": np.abs(errs).max(axis=1).tolist(),  # worst-M error per trial
        "elapsed_sec": elapsed,
    }


def theoretical_std(M: int, pi_N: int, pi_M: int, K: int) -> float:
    """Var[pi_hat(M; K)] = (M*pi_N - pi_M + pi_M^2 - pi_M) / K"""
    var = (M * pi_N - pi_M + pi_M * pi_M - pi_M) / K
    return float(np.sqrt(max(var, 0.0)))


def main() -> None:
    out_dir = Path("/apps/aplikacijos/prime-research/experiments/information_theory")
    out_path = out_dir / "shadow_tomography_chi_p_data.json"

    # Use N = 2^15 = 32768 for tractability. Each trial K*N ops.
    # K=1e5, N=3e4 -> 3e9 ops per trial; 5 trials -> 1.5e10 (~3-5 min).
    N = 2**15
    print(f"N = {N} = 2^{int(np.log2(N))}", flush=True)

    chi_p = sieve(N)
    pi_true = cumulative_pi(chi_p)
    pi_N = int(pi_true[-1])
    print(f"pi(N) = {pi_N}", flush=True)

    # Targets: dyadic M values
    M_log_max = int(np.log2(N))
    M_targets = np.array(
        [2**k for k in range(5, M_log_max + 1)], dtype=np.int64
    )  # 32, 64, ..., N
    print(f"M_targets = {M_targets.tolist()}", flush=True)

    # K sweep
    K_values = [100, 300, 1000, 3000, 10000, 30000, 100000]
    n_trials_per_K = 20  # for stable variance estimate

    results = {
        "N": N,
        "pi_N": pi_N,
        "M_targets": M_targets.tolist(),
        "pi_targets": [int(pi_true[m - 1]) for m in M_targets],
        "K_values": K_values,
        "n_trials_per_K": n_trials_per_K,
        "trials": [],
    }

    rng = np.random.default_rng(42)

    for K in K_values:
        # Memory: block_size shadows simultaneously. Want block_size * N * 1 byte < 200MB
        # bs * N < 2e8 -> bs < 6000 at N=2^15. Use 1024 to balance memory & vectorisation.
        bs = min(1024, K)
        # Lower n_trials for very large K
        nt = n_trials_per_K if K <= 10000 else 10 if K <= 30000 else 5
        print(f"\n=== K={K}, n_trials={nt}, block_size={bs} ===", flush=True)
        stats = trial_errors(chi_p, pi_true, M_targets, K, nt, rng, bs)

        # Theoretical std for comparison
        th_std = [
            theoretical_std(int(M), pi_N, int(pi_true[M - 1]), K) for M in M_targets
        ]
        stats["theoretical_std"] = th_std
        stats["empirical_std"] = stats["std_err"]
        stats["std_ratio"] = [
            (stats["std_err"][i] / th_std[i] if th_std[i] > 0 else float("nan"))
            for i in range(len(M_targets))
        ]

        # Print summary
        print(f"  elapsed: {stats['elapsed_sec']:.1f}s", flush=True)
        for i, M in enumerate(M_targets.tolist()):
            print(
                f"  M={M:>7d} pi(M)={pi_true[M-1]:>6d}  "
                f"emp_std={stats['empirical_std'][i]:>8.2f}  "
                f"th_std={th_std[i]:>8.2f}  "
                f"ratio={stats['std_ratio'][i]:>5.3f}"
            )
        L_inf_mean = float(np.mean(stats["L_inf_per_trial"]))
        print(f"  L_inf error (mean across trials) = {L_inf_mean:.2f}", flush=True)
        stats["L_inf_mean"] = L_inf_mean
        results["trials"].append(stats)

    # Save data
    with open(out_path, "w") as f:
        json.dump(results, f, indent=2)
    print(f"\nSaved {out_path}")

    # Quick summary table
    print("\n=== K-scaling of L_inf error ===")
    print(f"{'K':>8s}  {'L_inf':>10s}  {'theoretical_max_std':>22s}")
    for stats in results["trials"]:
        K = stats["K"]
        L_inf = stats["L_inf_mean"]
        max_th_std = max(stats["theoretical_std"])
        print(f"{K:>8d}  {L_inf:>10.2f}  {max_th_std:>22.2f}")

    # Effective scaling exponent: fit log L_inf = -alpha * log K + c
    Ks = np.array(K_values, dtype=float)
    L_inf = np.array([s["L_inf_mean"] for s in results["trials"]])
    log_K = np.log(Ks)
    log_L = np.log(L_inf)
    # Linear regression
    alpha, c = np.polyfit(log_K, log_L, 1)
    print(f"\nL_inf scaling: L_inf ~ K^{alpha:.3f} (theoretical: -0.5)")
    print(f"  intercept exp(c) = {np.exp(c):.2f}")

    # Required K for L_inf <= 1
    # L_inf ~ exp(c) * K^alpha = 1 -> K = exp(-c/alpha)
    if alpha < 0:
        K_for_eps_1 = np.exp(-c / alpha)
        print(f"\nExtrapolated K for L_inf <= 1: K ~ {K_for_eps_1:.2e}")
        print(f"  (compare to N * pi(N) = {N * pi_N:.2e})")
        print(f"  ratio: {K_for_eps_1 / (N * pi_N):.3f}")

    results["scaling_alpha"] = float(alpha)
    results["scaling_intercept_c"] = float(c)
    results["K_for_L_inf_1_extrapolated"] = float(np.exp(-c / alpha)) if alpha < 0 else None
    results["N_pi_N_baseline"] = N * pi_N

    with open(out_path, "w") as f:
        json.dump(results, f, indent=2)


if __name__ == "__main__":
    main()
