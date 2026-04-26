"""
P3: Random-matrix-theory moment-matching for delta(n) = pi(x) - li(x).

Hypothesis (Montgomery-Dyson + Cramer-style reasoning):
  Under the GUE conjecture for zeta zeros, the sequence Delta(x) = pi(x) - li(x),
  viewed as a function of x with x in a window [X-H, X+H], has predictable
  low-order moments. In particular,
      var(Delta) ~ x / (log x)^2     (Cramer)
      higher moments determined by zero correlations.

If the hypothesis is correct in a SHARP form, then knowing
  Delta(X - H), ..., Delta(X + H)
should let us *predict* the value of Delta(X) (or equivalently pi(X)) within
+/- 0.5, even though X is the central point we want.

We test the simplest Bayesian estimator: "Delta(X) ≈ mean of nearby
Delta(x)". If the hypothesis is right (i.e., delta is smooth on
small scales), this works. If delta is essentially uncorrelated noise,
it fails.

Methodology:
  1. For each test X in {1000, 5000, 10000, 50000}, compute
     Delta(x) for x in [X - H, X + H] using a fast sieve / exact pi.
  2. Estimator A: arithmetic mean of {Delta(x) : x != X}.
  3. Estimator B: weighted mean with weights 1/|x - X|^alpha (alpha = 0, 0.5, 1).
  4. Estimator C: linear regression on Delta vs li(x), use residual structure.
  5. Compare to true Delta(X).

Pass criterion: estimator gets Delta(X) right within +/- 0.5
(integer accuracy needed for pi(X) integer).

Note: even if estimator is good, this gives pi(x) ONLY at the center of a
window — but if H = polylog(X) primes-tested, the cost is polylog. The
question is whether the prediction error is small enough.
"""

import sys
import numpy as np
from mpmath import mp, mpf, li, log
from sympy.ntheory import isprime

mp.prec = 60


def pi_via_test(x_max, x_min=2):
    """Compute (x, pi(x)) for x in [x_min, x_max] using an incremental sieve.
       Returns dict x -> pi(x).
    """
    # Use direct sieve
    sieve = np.ones(int(x_max) + 1, dtype=bool)
    sieve[:2] = False
    for i in range(2, int(x_max ** 0.5) + 1):
        if sieve[i]:
            sieve[i * i::i] = False
    cum = np.cumsum(sieve)
    return {int(x): int(cum[x]) for x in range(x_min, x_max + 1)}


def main():
    print("=" * 70)
    print("P3: RMT moment-matching for Delta(x) = pi(x) - li(x)")
    print("=" * 70)

    test_Xs = [1000, 5000, 10000, 50000]
    H = 100  # half-window

    rmsm_per_X = []  # collect average estimator skill
    print(f"\nWindow half-width H = {H}.")
    print("\nFor each X, compute Delta(x) in [X-H, X+H] and use the OTHER points")
    print("to predict Delta(X). Compare to truth.\n")

    for X in test_Xs:
        x_max = X + H
        pi_dict = pi_via_test(x_max)
        deltas_in_window = {}
        for x in range(X - H, X + H + 1):
            li_x = float(li(mpf(x)))
            deltas_in_window[x] = pi_dict[x] - li_x

        true_delta_X = deltas_in_window[X]

        # Estimator A: arithmetic mean of all OTHER points in window
        others = [d for x, d in deltas_in_window.items() if x != X]
        est_mean = float(np.mean(others))

        # Estimator B: distance-weighted mean
        weights = []
        vals = []
        for x, d in deltas_in_window.items():
            if x == X:
                continue
            weights.append(1.0 / abs(x - X))  # alpha=1
            vals.append(d)
        weights = np.array(weights)
        vals = np.array(vals)
        est_weighted = float(np.average(vals, weights=weights))

        # Estimator C: median (robust to outliers)
        est_median = float(np.median(others))

        # Estimator D: linear regression Delta on (x - X), predict at x = X
        xs = np.array([x for x in deltas_in_window if x != X])
        ds = np.array([deltas_in_window[x] for x in xs])
        slope, intercept = np.polyfit(xs, ds, 1)
        est_linreg = slope * X + intercept

        # baseline: assume Delta(X) = 0 (i.e., pi(X) = li(X))
        est_zero = 0.0

        print(f"X = {X:>6}: true Delta(X) = {true_delta_X:+.4f}")
        print(f"  zero-baseline error:    {est_zero - true_delta_X:+.4f}")
        print(f"  arithmetic mean error:  {est_mean - true_delta_X:+.4f}")
        print(f"  weighted mean error:    {est_weighted - true_delta_X:+.4f}")
        print(f"  median error:           {est_median - true_delta_X:+.4f}")
        print(f"  linreg error:           {est_linreg - true_delta_X:+.4f}")
        print(f"  std(Delta) in window:   {np.std(others):.4f}")
        print()

        rmsm_per_X.append({
            "X": X,
            "true": true_delta_X,
            "zero_err": est_zero - true_delta_X,
            "mean_err": est_mean - true_delta_X,
            "weighted_err": est_weighted - true_delta_X,
            "median_err": est_median - true_delta_X,
            "linreg_err": est_linreg - true_delta_X,
            "std_window": np.std(others),
        })

    # Aggregate skill
    print("=" * 70)
    print("Summary: RMSE of each estimator across test Xs")
    methods = ["zero_err", "mean_err", "weighted_err", "median_err", "linreg_err"]
    for m in methods:
        rmse = np.sqrt(np.mean([r[m] ** 2 for r in rmsm_per_X]))
        print(f"  {m:>15}: RMSE = {rmse:.4f}")

    # Pass criterion: error < 0.5
    best_method = min(methods, key=lambda m: np.sqrt(np.mean([r[m] ** 2 for r in rmsm_per_X])))
    best_rmse = np.sqrt(np.mean([r[best_method] ** 2 for r in rmsm_per_X]))
    print(f"\nBest estimator: {best_method} (RMSE = {best_rmse:.4f})")
    if best_rmse < 0.5:
        print("PASS: best estimator gets Delta(X) within +/- 0.5 on average.")
        return 0
    else:
        print("FAIL: no estimator beats 0.5 RMSE; Delta(X) is NOT predictable")
        print("from neighbors via local moments alone.")
        return 1


if __name__ == "__main__":
    sys.exit(main())
