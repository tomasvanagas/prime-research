"""
PSLQ search for closed-form expressions of delta(n) = p(n) - R^{-1}(n)
on sparse subsequences:
  - dyadic n = 2^k
  - Fibonacci indices n_k = F_k
  - prime indices n_k = p_k (i.e., delta evaluated at prime indices)

Hypothesis: even though delta is incompressible globally, on a
polylog-indexable subsequence it might be expressible in fundamental
constants.
"""
from mpmath import mp, mpf, mpc, log, sqrt, exp, pi, euler, zeta, pslq
from sympy import prime, fibonacci

mp.dps = 80


def riemann_R(x):
    """Riemann's R(x) = sum_{k=1}^inf mu(k)/k * li(x^{1/k})."""
    from sympy import mobius
    s = mpf(0)
    x = mpf(x)
    for k in range(1, 60):
        m = mobius(k)
        if m == 0:
            continue
        # li(x^{1/k}) using mpmath li
        from mpmath import li as mp_li
        try:
            term = mpf(m) / k * mp_li(x ** (mpf(1) / k))
        except Exception:
            break
        s += term
        if abs(term) < mp.mpf(10) ** (-mp.dps + 5):
            break
    return s


def R_inv(n, x_lo=None, x_hi=None):
    """Solve R(x) = n by bisection."""
    n = mpf(n)
    if x_lo is None:
        # rough seed: p_n ~ n log n
        x_lo = max(mpf(2), n * log(max(n, 2)) * mpf("0.5"))
    if x_hi is None:
        x_hi = n * log(max(n, 2)) * mpf("3.0") + mpf(100)
    # ensure bracket
    while riemann_R(x_lo) > n:
        x_lo /= 2
    while riemann_R(x_hi) < n:
        x_hi *= 2
    for _ in range(200):
        mid = (x_lo + x_hi) / 2
        v = riemann_R(mid)
        if abs(v - n) < mpf(10) ** (-mp.dps + 5):
            return mid
        if v < n:
            x_lo = mid
        else:
            x_hi = mid
    return mid


def delta_at(n_idx):
    """delta(n) = p_n - R^{-1}(n)."""
    pn = mpf(int(prime(n_idx)))
    rn = R_inv(n_idx)
    return pn - rn


def try_pslq(target, dictionary, label):
    vec = [target] + dictionary
    rel = pslq(vec, tol=mpf("1e-50"), maxcoeff=10 ** 12)
    if rel is None:
        return None
    # filter trivial
    if rel[0] == 0:
        return None
    if max(abs(r) for r in rel) > 10 ** 8:
        return None
    return rel


def main():
    # Dictionary of fundamental constants
    one = mpf(1)
    z2 = zeta(2)
    z3 = zeta(3)
    g = euler
    log2 = log(2)
    log3 = log(3)
    logpi = log(pi)
    logtwopi = log(2 * pi)

    base_dict = [one, z2, z3, g, log2, logpi, logtwopi]
    print("Dictionary:", ["1", "zeta(2)", "zeta(3)", "gamma", "log 2",
                          "log pi", "log(2 pi)"])

    # Subsequence A: n = 2^k
    print("\n=== Dyadic subsequence: n = 2^k ===")
    deltas_dyadic = []
    for k in range(1, 12):
        n = 2 ** k
        d = delta_at(n)
        deltas_dyadic.append((n, d))
        print(f"  k={k:2d}, n=2^{k}={n:>5d}, p_n={int(prime(n)):>6d}, delta={float(d):+.6f}")

    print("\n  Per-point PSLQ tests:")
    for n, d in deltas_dyadic:
        for scale_name, scale in [("1", one), ("log n", log(n)),
                                   ("sqrt n", sqrt(n)), ("log log n", log(log(n)))]:
            target = d * scale
            rel = try_pslq(target, base_dict,
                           f"n=2^{int(round(log(n)/log2))}, scale={scale_name}")
            if rel is not None:
                print(f"  HIT n={n}, scale={scale_name}: {rel}")

    # Try cross-point relations: do delta(2^k) for k=1..K satisfy a
    # linear recurrence with constant coefficients in the dictionary?
    print("\n=== Cross-point: delta(2^k) for k=1..10 vs dictionary ===")
    deltas_only = [d for _, d in deltas_dyadic]
    # PSLQ on (delta_1, delta_2, ..., delta_10, basis) — too many unknowns
    # Instead check: is there a generating-function style relation
    # Lambda_k * delta(2^k) = sum_j c_kj * basis_j ?
    # Try first three points jointly: maybe delta(2^k) = a + b * k + c * 2^{-k} * something
    # Simpler: check if delta(2^k) / log(2^k) tends to a constant
    print("  delta(n)/log n at n=2^k:")
    for k in range(1, len(deltas_dyadic) + 1):
        n, d = deltas_dyadic[k - 1]
        ratio = d / log(n)
        print(f"    k={k:2d}: {float(ratio):+.6f}")

    # Subsequence B: Fibonacci indices
    print("\n=== Fibonacci index subsequence: n = F_k ===")
    fibs = []
    for k in range(3, 15):
        fk = int(fibonacci(k))
        if fk < 2:
            continue
        try:
            d = delta_at(fk)
            fibs.append((fk, d))
            print(f"  k={k:2d}, F_k={fk:>5d}, delta={float(d):+.6f}")
        except Exception as e:
            print(f"  k={k}: failed: {e}")
            break

    for n, d in fibs:
        for scale_name, scale in [("1", one), ("log n", log(n))]:
            target = d * scale
            rel = try_pslq(target, base_dict, f"n=F_, scale={scale_name}")
            if rel is not None:
                print(f"  HIT n={n}, scale={scale_name}: {rel}")

    # Subsequence C: prime indices
    print("\n=== Prime-index subsequence: n = p_k ===")
    pidx = []
    for k in range(1, 15):
        pk = int(prime(k))
        try:
            d = delta_at(pk)
            pidx.append((pk, d))
            print(f"  k={k:2d}, p_k={pk:>5d}, delta={float(d):+.6f}")
        except Exception as e:
            print(f"  k={k}: failed: {e}")
            break

    for n, d in pidx:
        for scale_name, scale in [("1", one), ("log n", log(n))]:
            target = d * scale
            rel = try_pslq(target, base_dict, f"n=p_k, scale={scale_name}")
            if rel is not None:
                print(f"  HIT n={n}, scale={scale_name}: {rel}")


if __name__ == "__main__":
    main()
