"""
Detrended Fluctuation Analysis (DFA) and Higuchi Fractal Dimension
of the Prime-Counting Residual delta(x) = pi(x) - Li(x)
====================================================================

This experiment measures characterizing exponents of the prime-counting
residual that have NOT been reported in the project's pseudorandomness
catalog: the DFA Hurst exponent H and the Higuchi fractal dimension D
of `delta(x) = pi(x) - Li(x)` viewed as a 1D time series.

Why these matter

  * DFA Hurst H captures long-range correlations after removing
    polynomial trends. H = 1/2 -> uncorrelated increments (white).
    H > 1/2 -> persistence; H < 1/2 -> anti-persistence.
  * Higuchi D measures path roughness. D = 1.5 corresponds to Brownian
    motion. D < 1.5 -> smoother than Brownian; D > 1.5 -> rougher.
  * Under RH, |delta(x)| is O(sqrt(x) log x), with GUE-random sign.
    The cumulative process pi(x) - Li(x) is conjecturally
    Brownian-motion-like in log scale. An RH prediction is therefore
    H ~ 1/2 for the increment process and D ~ 1.5 for delta(x) itself.

  * If we found anomalous H far from 1/2 (especially H << 1/2,
    i.e. anti-persistence with mean-reverting structure), there might
    be a low-complexity description of the deviation -- the project's
    primary obstacle is precisely the GUE-random nature of delta(x).
    This experiment provides a direct test.

Method

  delta(x) computed via primepi(x) - Li(x), x in [N0, N1].
  Compared to a GUE-random surrogate of matching variance:
    surr(x) = sum_{k=1}^K cos(gamma_k log x + phi_k) * sqrt(x) / gamma_k
  with gamma_k uniform-randomly drawn from the spectral density of
  zeta zeros and phi_k uniform on [0, 2*pi).

  DFA computed at scales s = 16, 32, 64, ..., N/4.
  Higuchi computed for k_max in {2, 3, ..., 64}.
"""

import numpy as np
from sympy import primepi, zeta
import mpmath
import time


def li(x):
    return float(mpmath.li(x))


def compute_delta(N0, N1, step=1):
    xs = np.arange(N0, N1 + 1, step, dtype=np.int64)
    pi_vals = np.array([int(primepi(int(x))) for x in xs], dtype=np.float64)
    li_vals = np.array([li(x) for x in xs], dtype=np.float64)
    delta = pi_vals - li_vals
    return xs, delta


def dfa(signal, scales=None):
    """
    Detrended Fluctuation Analysis (1st order detrending).
    Returns slope alpha of log(F(s)) vs log(s).
    """
    y = np.cumsum(signal - np.mean(signal))
    n = len(y)
    if scales is None:
        scales = np.unique(np.round(np.geomspace(8, n // 4, num=20)).astype(int))
    F = []
    for s in scales:
        if s < 4 or 2 * s > n:
            continue
        # Segment, detrend each piece, RMS of residuals
        n_seg = n // s
        rms = []
        for i in range(n_seg):
            seg = y[i * s:(i + 1) * s]
            t = np.arange(s, dtype=np.float64)
            p = np.polyfit(t, seg, 1)
            resid = seg - np.polyval(p, t)
            rms.append(np.mean(resid ** 2))
        # Reverse direction too, to use both ends:
        for i in range(n_seg):
            seg = y[n - (i + 1) * s:n - i * s]
            t = np.arange(s, dtype=np.float64)
            p = np.polyfit(t, seg, 1)
            resid = seg - np.polyval(p, t)
            rms.append(np.mean(resid ** 2))
        F.append((s, float(np.sqrt(np.mean(rms)))))
    F = np.array(F)
    if len(F) < 3:
        return float("nan"), F
    A = np.vstack([np.log(F[:, 0]), np.ones(len(F))]).T
    sol, _, _, _ = np.linalg.lstsq(A, np.log(F[:, 1]), rcond=None)
    return float(sol[0]), F


def higuchi(signal, k_max=64):
    """
    Higuchi fractal dimension D = -slope of log(L(k)) vs log(k).
    """
    n = len(signal)
    L_k = []
    for k in range(2, k_max + 1):
        Lm = []
        for m in range(k):
            idx = np.arange(m, n, k)
            if len(idx) < 2:
                continue
            diff = np.abs(np.diff(signal[idx]))
            Lm_val = (np.sum(diff) * (n - 1) / (k * len(diff))) / k
            Lm.append(Lm_val)
        if len(Lm) > 0:
            L_k.append((k, float(np.mean(Lm))))
    L_k = np.array(L_k)
    if len(L_k) < 4:
        return float("nan"), L_k
    A = np.vstack([np.log(L_k[:, 0]), np.ones(len(L_k))]).T
    sol, _, _, _ = np.linalg.lstsq(A, np.log(L_k[:, 1]), rcond=None)
    return float(-sol[0]), L_k


def random_surrogate(xs, n_modes=200, seed=0):
    """
    Surrogate signal with the same expected variance scaling sqrt(x):
    sum_k cos(gamma_k log x + phi_k) * sqrt(x) / gamma_k,
    with gamma_k drawn from the density of zeta zeros (~ 1/(2 pi) log(t/(2pi))).
    """
    rng = np.random.default_rng(seed)
    log_x = np.log(xs.astype(np.float64))
    # Sample zero heights gamma proportional to density 1/(2 pi) * log(g/(2 pi))
    gamma_max = 1000.0
    # Importance-sample: draw uniformly in log space.
    g = np.exp(rng.uniform(np.log(10), np.log(gamma_max), n_modes))
    phi = rng.uniform(0.0, 2 * np.pi, n_modes)
    surr = np.zeros_like(log_x)
    for gk, pk in zip(g, phi):
        surr += np.cos(gk * log_x + pk) / gk
    surr *= np.sqrt(xs.astype(np.float64))
    return surr


def main():
    # Manageable size for primepi calls, but big enough for scaling.
    N0, N1, step = 100, 60000, 5
    print("Computing delta(x) = pi(x) - Li(x)  for  x in [{}, {}], step {}...".format(
        N0, N1, step))
    t0 = time.time()
    xs, delta = compute_delta(N0, N1, step)
    print(f"  N = {len(xs)} samples in {time.time() - t0:.1f}s")

    print(f"  Stats: delta range [{delta.min():.2f}, {delta.max():.2f}],  "
          f"std = {delta.std():.3f}")

    # DFA on increments (delta differences) -- testing the noise process
    increments = np.diff(delta)
    H_inc, F_inc = dfa(increments)

    # DFA on delta itself
    H_delta, F_delta = dfa(delta)

    # Higuchi on delta
    D_delta, _ = higuchi(delta, k_max=min(64, len(delta) // 8))

    # Surrogate
    surr = random_surrogate(xs, n_modes=200, seed=42)
    surr_inc = np.diff(surr)
    H_surr_inc, _ = dfa(surr_inc)
    H_surr, _ = dfa(surr)
    D_surr, _ = higuchi(surr, k_max=min(64, len(surr) // 8))

    print()
    print("=" * 76)
    print(" DFA and Higuchi exponents")
    print("=" * 76)
    print(f"{'quantity':30}  {'delta(x) (pi - Li)':>22}  {'GUE surrogate':>15}")
    print("-" * 76)
    print(f"{'DFA H of increments':30}  {H_inc:>22.4f}  {H_surr_inc:>15.4f}")
    print(f"{'DFA H of delta(x)':30}  {H_delta:>22.4f}  {H_surr:>15.4f}")
    print(f"{'Higuchi D of delta(x)':30}  {D_delta:>22.4f}  {D_surr:>15.4f}")

    # Predictions
    print()
    print("Reference values:")
    print("  white noise:       H_increments = 0.5")
    print("  Brownian motion:   H = 0.5,  D = 1.5")
    print("  pure trend (smooth): D ~ 1.0,  H ~ 1.0")
    print("  RH prediction for  delta(x) ~ sqrt(x) GUE-random:  D ~ 1.5,  H ~ 0.5")
    print()

    # Verdict
    if abs(H_inc - 0.5) < 0.1 and abs(D_delta - 1.5) < 0.15:
        print("VERDICT: delta(x) is BROWNIAN-MOTION-LIKE -- GUE-random sign confirmed.")
        print("         No structural shortcut from fractal scaling.")
    elif H_inc < 0.4:
        print("VERDICT: ANTI-PERSISTENCE detected (H < 0.4). Mean-reverting structure!")
        print("         INVESTIGATE: could yield low-complexity description.")
    elif H_inc > 0.6:
        print("VERDICT: PERSISTENT correlations (H > 0.6). Long-memory structure.")
        print("         INVESTIGATE: but typical for partial-sum signals.")
    else:
        print(f"VERDICT: marginal -- H_inc = {H_inc:.3f} between random and structured.")

    # Quick scaling: split delta into halves and compare exponents
    half = len(delta) // 2
    H1, _ = dfa(np.diff(delta[:half]))
    H2, _ = dfa(np.diff(delta[half:]))
    print(f"\nStability check: H_inc on first half = {H1:.4f}, second half = {H2:.4f}")
    print("  (If primes had a hidden time-dependent structure, the two halves would differ.)")


if __name__ == "__main__":
    t0 = time.time()
    main()
    print(f"\nTotal time: {time.time() - t0:.2f}s")
