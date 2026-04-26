"""
Proposal B5: Self-correcting recursive halving p(n) -> p(n/2).

Define theta(n) := p(n) - 2*p(floor(n/2)).  By PNT,
    p(n)   ~ n log n,
    2*p(n/2) ~ 2 (n/2) log(n/2) = n log n - n log 2,
so theta(n) ~ n log 2 + (lower-order chaotic).

Question: is the chaotic part of theta(n) *smaller* than the chaotic
part of p(n) alone?  If yes, recursive halving can reduce the per-step
unpredictability and yield a polylog algorithm via depth-log-n recursion.

This script:
  1. computes theta(n) for n = 100, 200, ..., n_max.
  2. fits theta(n) = c_0 + c_1 n + c_2 n log n + c_3 n log log n
                   + c_4 sqrt(n) + (residual).
  3. measures residual scaling as a function of n (does residual / n^alpha
     vanish for some alpha < 1?).

Decision rule:
  - residual ~ const  ->  theta is *smoothly* predictable (WIN)
  - residual ~ n^{0.5}  ->  reduces to original problem (CLOSED)
  - residual ~ n^a, 0 < a < 0.5  ->  ambiguous (refine the experiment)

Usage:
    python3 halving_recursion.py
    python3 halving_recursion.py --n_max 50000 --step 100
"""
from __future__ import annotations
import argparse
import json
import time
import numpy as np
import sympy


def prime_table(n_max: int) -> np.ndarray:
    """p(n) for n in [0, n_max], indexed so out[n] = nth prime."""
    print(f"[setup] computing first {n_max} primes")
    t0 = time.time()
    out = np.zeros(n_max + 1, dtype=np.int64)
    p = 1
    for n in range(1, n_max + 1):
        p = sympy.nextprime(p)
        out[n] = p
    print(f"[setup] done in {time.time()-t0:.2f}s")
    return out


def fit_residual(n_grid: np.ndarray, theta: np.ndarray,
                 basis: list[str]) -> tuple[np.ndarray, np.ndarray]:
    """
    Fit theta = sum_k c_k * basis_k(n) by least squares.
    Returns (coefficients, residual array).
    """
    n = n_grid.astype(float)
    cols = []
    for term in basis:
        if term == "1":
            cols.append(np.ones_like(n))
        elif term == "n":
            cols.append(n)
        elif term == "n_logn":
            cols.append(n * np.log(n))
        elif term == "n_log_logn":
            cols.append(n * np.log(np.log(n)))
        elif term == "sqrt_n":
            cols.append(np.sqrt(n))
        elif term == "logn":
            cols.append(np.log(n))
        elif term == "log2n":
            cols.append(np.log(n) ** 2)
        elif term == "n_logn_squared":
            cols.append(n * np.log(n) ** 2)
        else:
            raise ValueError(f"unknown basis term {term}")
    A = np.column_stack(cols)
    coefs, *_ = np.linalg.lstsq(A, theta, rcond=None)
    pred = A @ coefs
    residual = theta - pred
    return coefs, residual


def power_law_fit(x: np.ndarray, y: np.ndarray) -> tuple[float, float, float]:
    """
    Fit |y| ~ a * x^b in log-log.  Returns (a, b, r2).
    """
    mask = (x > 0) & (np.abs(y) > 0)
    lx = np.log(x[mask])
    ly = np.log(np.abs(y[mask]))
    A = np.column_stack([lx, np.ones_like(lx)])
    sol, *_ = np.linalg.lstsq(A, ly, rcond=None)
    b, log_a = sol
    a = float(np.exp(log_a))
    pred = lx * b + log_a
    ss_res = float(np.sum((ly - pred) ** 2))
    ss_tot = float(np.sum((ly - np.mean(ly)) ** 2))
    r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else 0.0
    return a, float(b), r2


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--n_max", type=int, default=10000)
    parser.add_argument("--step", type=int, default=10)
    args = parser.parse_args()

    primes = prime_table(args.n_max)

    n_grid = np.arange(args.step * 2, args.n_max + 1, args.step)
    p = primes[n_grid]
    p_half = primes[n_grid // 2]
    theta = (p - 2 * p_half).astype(float)

    n_log2 = n_grid * np.log(2.0)
    leading = theta - n_log2

    # try basis combinations
    bases = [
        ["1", "n", "n_logn", "n_log_logn", "logn", "log2n"],
        ["1", "n", "n_logn", "n_log_logn", "logn", "log2n", "sqrt_n"],
    ]

    print("\n[fit] theta(n) regressions on n in "
          f"[{n_grid[0]}, {n_grid[-1]}] step {args.step}")
    print(f"      mean theta(n) = {theta.mean():.3e}, "
          f"std theta(n) = {theta.std():.3e}")
    print(f"      mean (theta - n log 2) = {leading.mean():.3e}")

    fit_results = {}
    for basis in bases:
        key = ",".join(basis)
        coefs, residual = fit_residual(n_grid.astype(float), theta, basis)
        a_pl, b_pl, r2 = power_law_fit(n_grid.astype(float),
                                       residual.astype(float))
        rmse = float(np.sqrt(np.mean(residual ** 2)))
        max_abs = float(np.max(np.abs(residual)))
        print(f"\n  basis: {basis}")
        for term, c in zip(basis, coefs):
            print(f"    coef[{term:>14s}] = {c:+.6e}")
        print(f"    RMSE residual = {rmse:.3e}, max |residual| = {max_abs:.3e}")
        print(f"    power-law fit |residual| ~ "
              f"{a_pl:.3e} * n^{b_pl:.3f}   (R^2 = {r2:.3f})")
        fit_results[key] = {
            "basis": basis,
            "coefs": [float(c) for c in coefs],
            "rmse": rmse,
            "max_abs": max_abs,
            "power_law_a": a_pl,
            "power_law_b": b_pl,
            "power_law_r2": r2,
        }

    # decision
    last = fit_results[",".join(bases[-1])]
    print("\n[decision]")
    if last["power_law_b"] < 0.3:
        verdict = "STRUCTURAL: residual exponent < 0.3 (sub-sqrt!)"
    elif 0.3 <= last["power_law_b"] < 0.55:
        verdict = "MAYBE: residual exponent in [0.3, 0.55] -- borderline"
    else:
        verdict = (f"CLOSED: residual exponent {last['power_law_b']:.2f} "
                   f"~ 0.5 -- reduces to original problem")
    print(f"  exponent b = {last['power_law_b']:.3f}, R^2 = "
          f"{last['power_law_r2']:.3f}")
    print(f"  verdict: {verdict}")

    out = {
        "args": vars(args),
        "n_grid_first10": n_grid[:10].tolist(),
        "n_grid_last10": n_grid[-10:].tolist(),
        "theta_first10": theta[:10].tolist(),
        "theta_last10": theta[-10:].tolist(),
        "fit_results": fit_results,
        "verdict": verdict,
    }
    out_path = (
        "/apps/aplikacijos/prime-research/experiments/proposals/"
        "halving_recursion_data.json"
    )
    with open(out_path, "w") as f:
        json.dump(out, f, indent=2, default=float)
    print(f"\n[save] wrote {out_path}")


if __name__ == "__main__":
    main()
