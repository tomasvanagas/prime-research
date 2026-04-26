"""
Reservoir computing test for delta(n) = pi(n) - R^{-1}(n).

Idea (fresh perspective, session 58):
   Echo-state networks (ESN) and other reservoir computers route input through
   a high-dimensional, fixed, chaotic dynamical system, then learn a *linear*
   readout. Empirically they capture long-range temporal correlations that
   explicit feature engineering misses. Many "irreducible" sequences turn out
   to have structure that an ESN can predict.
   We feed the binary expansion of n into a fixed random ESN, read out the
   reservoir state, and try to predict delta(n) by linear regression.
   The bar is: held-out MSE substantially below var(delta).

Pass/fail
   FAIL  -> reservoir captures no information beyond mean (R^2 near 0 on test).
            Adds a 23rd pseudorandomness measure: delta is incompressible by
            generic chaotic dynamics + linear readout.
   PASS  -> R^2 > 0.1 on held-out.  Worth pursuing: the reservoir state
            is a polylog-sized vector, so any predictive power means
            polylog encoding exists in some form.

Implementation
   - Compute pi(n) for n in [0, N) via segmented sieve.
   - Compute R^{-1}(n) via Lambert W on the prime number theorem inverse.
   - Drive an ESN of size D=200 with the bits of n (LSB first, fixed length).
   - Train ridge regression on first 70%, test on last 30%.
   - Report MSE_test, R^2_test, and compare to var(delta) baseline.
"""

import numpy as np
from scipy.special import wright_bessel  # not used; just availability check
from numpy.linalg import lstsq

# ------------------------------------------------------------------ helpers

def sieve(N: int) -> np.ndarray:
    """Boolean array, sieve[k] = True iff k is prime, for 0 <= k < N."""
    sieve = np.ones(N, dtype=bool)
    sieve[:2] = False
    for p in range(2, int(N**0.5) + 1):
        if sieve[p]:
            sieve[p * p :: p] = False
    return sieve


def riemann_R_inv(n: int) -> float:
    """Inverse of Riemann's R(x) at integer n, evaluated by Newton on R."""
    # R(x) = sum_{k>=1} mu(k)/k * Li(x^{1/k}). Use a simple truncation (k<=8).
    from sympy import mobius
    from mpmath import li, mpf, mp

    mp.dps = 25

    def R(x):
        s = mpf(0)
        xx = mpf(x)
        for k in range(1, 9):
            mu = mobius(k)
            if mu == 0:
                continue
            s += mpf(mu) / k * li(xx ** (mpf(1) / k))
        return s

    # Newton's method.  Seed with x = n * log n (PNT inverse).
    x = max(2.0, n * np.log(max(n, 3)))
    for _ in range(40):
        f = R(x) - n
        # derivative dR/dx ~ 1/log(x)
        dR = mpf(1) / mp.log(x)
        step = float(f / dR)
        x = max(2.0, x - step)
        if abs(step) < 1e-6:
            break
    return float(x)


def riemann_R(x: float, K: int = 12) -> float:
    """Riemann's R(x), used in reverse: pi(n) ~ R(n) gives R^{-1}(n) ~= n.
    Actually we want delta(n) = pi(n) - R^{-1}(n); easier path: define
    delta(n) = pi(n) - approx(n) where approx(n) is the smooth part
    R(n) of pi inverted at n.  We approximate by computing R^{-1}(n)
    from a precomputed table of (R(x), x) values."""
    from sympy import mobius
    from mpmath import li, mpf, mp

    mp.dps = 25
    s = mpf(0)
    for k in range(1, K + 1):
        mu = mobius(k)
        if mu == 0:
            continue
        s += mpf(mu) / k * li(mpf(x) ** (mpf(1) / k))
    return float(s)


def build_R_inverse_table(N: int) -> np.ndarray:
    """For n in [0, N), return R^{-1}(n) by inverting tabulated R via interp.

    R^{-1}(N) ~ N * log N for large N.  Build a logarithmic grid up to
    a safe upper bound `x_hi` such that R(x_hi) >> N.
    """
    import math
    x_hi = max(20.0, 4.0 * N * (math.log(max(N, 4))))
    xs = np.geomspace(1.5, x_hi, 5000)
    ys = np.array([riemann_R(float(x), K=10) for x in xs])
    # Make ys strictly monotone (numerical noise can break this).
    for i in range(1, len(ys)):
        if ys[i] <= ys[i - 1]:
            ys[i] = ys[i - 1] + 1e-9
    ns = np.arange(N, dtype=float)
    Rinv = np.interp(ns, ys, xs, left=xs[0], right=xs[-1])
    # Sanity-check: R(Rinv[k]) should be close to k.
    return Rinv


# ------------------------------------------------------------- reservoir

def make_reservoir(D: int, in_dim: int, spectral_radius: float = 0.9, seed: int = 0):
    rng = np.random.default_rng(seed)
    W = rng.standard_normal((D, D)) / np.sqrt(D)
    eigvals = np.linalg.eigvals(W)
    W *= spectral_radius / max(abs(eigvals).max(), 1e-9)
    Win = rng.standard_normal((D, in_dim)) * 0.5
    return W, Win


def run_reservoir(W, Win, U, leak: float = 0.3):
    """U is (T, in_dim).  Returns states (T, D)."""
    D = W.shape[0]
    T = U.shape[0]
    state = np.zeros(D)
    states = np.empty((T, D))
    for t in range(T):
        pre = W @ state + Win @ U[t]
        state = (1 - leak) * state + leak * np.tanh(pre)
        states[t] = state
    return states


def n_to_bits(n: int, width: int) -> np.ndarray:
    """LSB-first bit array of length `width`."""
    return np.array([(n >> i) & 1 for i in range(width)], dtype=float)


# ----------------------------------------------------------------- main

def main():
    N = 60_000
    D = 200          # reservoir size
    bit_width = 20   # log2(60000) ~= 16, give some headroom
    train_frac = 0.7
    ridge = 1e-3

    print(f"# reservoir test, N={N}, D={D}, bits={bit_width}")

    # Targets.  We define
    #   delta(n) := pi(n) - R(n)
    # i.e. the residual of pi from Riemann's (polylog-computable) smooth
    # approximation R(n) = sum_k mu(k)/k * li(n^{1/k}).
    # |delta(n)| is bounded by O(sqrt(n) log n) under RH, and is conjecturally
    # GUE-distributed.  This is the quantity whose compressibility we want
    # to test.
    pi_arr = np.cumsum(sieve(N).astype(np.int64)).astype(np.float64)
    print("# pi computed")
    print("# computing R(n) for n in [2, N)...")
    # R(n) for integer n; R(0)=R(1)=0 conventionally.
    R_arr = np.zeros(N, dtype=np.float64)
    # Build via interpolation from a logarithmic grid (R is smooth, monotone).
    xs = np.geomspace(1.5, max(N, 100) * 1.2, 4000)
    ys = np.array([riemann_R(float(x), K=10) for x in xs])
    ns = np.arange(2, N, dtype=float)
    R_arr[2:] = np.interp(ns, xs, ys, left=ys[0], right=ys[-1])
    delta = pi_arr - R_arr
    var_delta = float(np.var(delta[2:]))
    print(f"# var(delta) = {var_delta:.4f}, |delta| max = {np.max(np.abs(delta[2:])):.2f}")

    # Build random train / test split so trend / drift doesn't bias R^2.
    rng_split = np.random.default_rng(7)
    indices = np.arange(2, N)
    rng_split.shuffle(indices)
    n_train = int(len(indices) * train_frac)
    train_idx = np.sort(indices[:n_train])
    test_idx = np.sort(indices[n_train:])

    def fit_eval(feats):
        feats_b = np.concatenate([feats, np.ones((N, 1), dtype=np.float32)], axis=1)
        Xtr, ytr = feats_b[train_idx], delta[train_idx]
        Xte, yte = feats_b[test_idx], delta[test_idx]
        A = Xtr.T @ Xtr + ridge * np.eye(Xtr.shape[1], dtype=np.float32)
        w = np.linalg.solve(A, Xtr.T @ ytr.astype(np.float32))
        pred_te = Xte @ w
        return 1.0 - float(np.mean((pred_te - yte) ** 2)) / float(np.var(yte))

    # ---- Reservoir features over multiple seeds ----
    print("# running reservoirs (3 seeds)...")
    r2_list = []
    for seed in (1, 2, 3):
        W, Win = make_reservoir(D, in_dim=1, seed=seed)
        feats = np.empty((N, D), dtype=np.float32)
        for n in range(N):
            bits = n_to_bits(n, bit_width).reshape(-1, 1)
            states = run_reservoir(W, Win, bits)
            feats[n] = states[-1].astype(np.float32)
        r2 = fit_eval(feats)
        r2_list.append(r2)
        print(f"#   seed={seed}: test R^2 = {r2:.4f}")
    r2_res = float(np.mean(r2_list))

    # ---- Baseline 1: linear on raw bits ----
    bits_feats = np.empty((N, bit_width), dtype=np.float32)
    for n in range(N):
        bits_feats[n] = n_to_bits(n, bit_width)
    r2_bits = fit_eval(bits_feats)

    # ---- Baseline 2: linear on log n features ----
    log_feats = np.zeros((N, 4), dtype=np.float32)
    for n in range(N):
        if n >= 2:
            ln = np.log(n)
            log_feats[n] = [ln, np.log(ln), 1.0 / ln, n / max(ln, 1e-9)]
    r2_log = fit_eval(log_feats)

    # ---- Control: shuffled delta on reservoir features ----
    rng = np.random.default_rng(0)
    shuf = rng.permutation(delta.copy())
    real_delta = delta.copy()
    delta[:] = shuf
    W, Win = make_reservoir(D, in_dim=1, seed=1)
    feats_shuf = np.empty((N, D), dtype=np.float32)
    for n in range(N):
        bits = n_to_bits(n, bit_width).reshape(-1, 1)
        states = run_reservoir(W, Win, bits)
        feats_shuf[n] = states[-1].astype(np.float32)
    r2_shuf = fit_eval(feats_shuf)
    delta[:] = real_delta

    print()
    print(f"# reservoir mean test R^2 = {r2_res:.4f} (over 3 seeds)")
    print(f"# raw-bits  test R^2      = {r2_bits:.4f}")
    print(f"# log-feats test R^2      = {r2_log:.4f}")
    print(f"# shuffled  test R^2      = {r2_shuf:.4f}  (control)")

    # Verdict
    print()
    sig_threshold = max(r2_shuf + 0.02, 0.01)
    if r2_res > 0.10:
        print("VERDICT: PASS. Reservoir captures meaningful structure in delta.")
    elif r2_res > sig_threshold:
        print("VERDICT: WEAK. Reservoir slightly above shuffled control.")
    else:
        print("VERDICT: FAIL. Reservoir gives no significant signal on delta.")
    print(f"# significance threshold (shuffled+0.02) = {sig_threshold:.4f}")


if __name__ == "__main__":
    main()
