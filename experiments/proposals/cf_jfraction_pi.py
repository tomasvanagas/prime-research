"""F2: Continued-Fraction (J-fraction) Expansion of Prime Generating Function.

Conjecture: pi(n) is P-recursive, i.e. satisfies a finite-order linear recurrence
with polynomial-in-n coefficients. Equivalently, the generating function
G(z) = sum_{n>=1} (pi(n)/n) z^n is D-finite (holonomic).

Test: Compute Hankel determinants H_k = det([m_{i+j}]_{0<=i,j<=k}) of the moment
sequence m_n = pi(n)/n. If pi is D-finite, log|H_k| grows polynomially in k.
If pi is "random-like", log|H_k| grows like k^2 (typical of generic sequences).

We also extract the J-fraction coefficients (a_k, b_k) and check whether they
satisfy any low-order polynomial pattern.

Reference: Krattenthaler (1998) "Advanced determinant calculus".
"""

import sys
import time

import numpy as np
from sympy import primerange, Rational, Matrix, simplify, nsimplify


def pi_n_table(N):
    """pi(n) for n=1..N."""
    primes = list(primerange(2, N + 2))
    pi = [0] * (N + 1)
    j = 0
    for n in range(1, N + 1):
        while j < len(primes) and primes[j] <= n:
            j += 1
        pi[n] = j
    return pi[1:]  # length N


def hankel_det(moments, k):
    """det of (k+1)x(k+1) Hankel matrix M_{ij} = moments[i+j], i,j=0..k.

    Uses sympy for exact rationals.
    """
    n = k + 1
    M = Matrix(n, n, lambda i, j: moments[i + j])
    return M.det()


def jfraction_coeffs(moments, K):
    """Return (a_0, b_0), ..., (a_{K-1}, b_{K-1}) for the J-fraction expansion.

    Uses Hankel determinant ratios. Need 2K+1 moments.
    Formulas (Krattenthaler):
        b_0  = m_0
        D_k  = det Hankel of size k+1 on (m_0, ..., m_{2k})
        D'_k = det shifted Hankel of size k+1 on (m_1, ..., m_{2k+1})
        b_k  = D_k * D_{k-2} / D_{k-1}^2
        a_k  = D'_k / D_k - D'_{k-1} / D_{k-1}
    """
    if len(moments) < 2 * K + 1:
        raise ValueError(f"Need {2*K+1} moments, got {len(moments)}")
    D = [Rational(1)]  # D_{-1} = 1
    Dp = [Rational(1)]  # D'_{-1} = 1
    for k in range(K + 1):
        # Hankel size (k+1)x(k+1) starting at index 0
        M = Matrix(k + 1, k + 1, lambda i, j: moments[i + j])
        D.append(M.det())
        Mp = Matrix(k + 1, k + 1, lambda i, j: moments[i + j + 1])
        Dp.append(Mp.det())
    # D[0] = 1 (D_{-1}), D[1] = D_0, ..., D[K+1] = D_K
    a = []
    b = [moments[0]]  # b_0
    for k in range(1, K + 1):
        bk = D[k + 1] * D[k - 1] / (D[k] ** 2) if D[k] != 0 else None
        b.append(bk)
    for k in range(K):
        ak = Dp[k + 1] / D[k + 1] - Dp[k] / D[k] if D[k + 1] != 0 and D[k] != 0 else None
        a.append(ak)
    return a, b, D


def main():
    N_max = 80  # number of moments. K-th Hankel needs 2K+1, so K up to ~38
    K = 12  # how many J-fraction coefficients to extract

    print(f"Computing pi(n) for n=1..{N_max}")
    pi_vals = pi_n_table(N_max)

    # Use moments m_n = pi(n+1) / (n+1) for n=0,1,2,...
    # so m_0 = pi(1)/1 = 0; that's degenerate. Use shifted: m_n = pi(n+2)/(n+2)
    # m_0 = pi(2)/2 = 1/2.
    moments = [Rational(pi_vals[n + 1], n + 2) for n in range(N_max - 1)]

    print(f"\nFirst 10 moments (pi(n+2)/(n+2)):")
    for i, m in enumerate(moments[:10]):
        print(f"  m_{i} = {m} = {float(m):.6f}")

    # Hankel determinants
    print(f"\nHankel determinants log10|H_k|:")
    log_dets = []
    for k in range(0, min(K + 5, len(moments) // 2)):
        try:
            d = hankel_det(moments, k)
            ld = float(d).__abs__()
            if ld == 0:
                print(f"  H_{k} = 0  (perfect rank deficiency!)")
                log_dets.append(-float("inf"))
            else:
                lg = np.log10(ld) if ld > 0 else np.log10(abs(float(d)))
                log_dets.append(lg)
                print(f"  log10|H_{k}| = {lg:.3f}    (det={float(d):.4e})")
        except Exception as e:
            print(f"  H_{k}: error {e}")

    # Fit growth: log|H_k| vs k (linear) and vs k^2 (quadratic)
    valid = [(k, ld) for k, ld in enumerate(log_dets) if np.isfinite(ld)]
    if len(valid) >= 4:
        ks = np.array([v[0] for v in valid])
        lds = np.array([v[1] for v in valid])
        # Linear fit log|H_k| = A * k + B
        A_lin, B_lin = np.polyfit(ks, lds, 1)
        # Quadratic fit log|H_k| = A * k^2 + B * k + C
        if len(valid) >= 6:
            coefs = np.polyfit(ks, lds, 2)
            A_q, B_q, C_q = coefs
        else:
            A_q = None

        print(f"\nFits to log10|H_k|:")
        print(f"  Linear:   A*k + B,    A = {A_lin:.4f}, B = {B_lin:.4f}")
        if A_q is not None:
            print(f"  Quadratic A*k^2 + B*k + C, A = {A_q:.4f}, B = {B_q:.4f}, C = {C_q:.4f}")

        # If A_q is large compared to A_lin, growth is quadratic (random-like).
        # If A_q ~ 0 and A_lin > 0, growth is linear -> pi might be P-recursive!
        if A_q is not None:
            if abs(A_q) < 0.05 and abs(A_lin) > 0.1:
                print(f"\n*** SUGGESTIVE: linear growth dominates -> hint of D-finiteness! ***")
            elif A_q > 0.1:
                print(f"\nQuadratic growth -> consistent with NON-D-finite (random-like).")

    # J-fraction coefficients
    print(f"\nJ-fraction coefficients (first {K}):")
    try:
        K_eff = min(K, len(moments) // 2 - 1)
        a, b, _ = jfraction_coeffs(moments, K_eff)
        for k in range(K_eff):
            try:
                ak_f = float(a[k]) if a[k] is not None else None
                bk_f = float(b[k + 1]) if b[k + 1] is not None else None
                print(f"  a_{k} = {ak_f}     b_{k+1} = {bk_f}")
            except Exception:
                print(f"  a_{k} = {a[k]}     b_{k+1} = {b[k+1]}")
    except Exception as e:
        print(f"  J-fraction extraction failed: {e}")

    # Check for OEIS-style patterns in (a_k, b_k)
    # Simple check: do a_k or b_k stabilize to rationals with small denominators?
    print(f"\nLooking for arithmetic patterns in J-fraction coefficients:")
    try:
        a_floats = [float(ak) for ak in a if ak is not None]
        b_floats = [float(bk) for bk in b[1:] if bk is not None]
        a_diff = np.diff(a_floats) if len(a_floats) > 1 else []
        b_ratio = [b_floats[i + 1] / b_floats[i] if b_floats[i] != 0 else None
                   for i in range(len(b_floats) - 1)]
        print(f"  a_k differences: {[f'{d:.4f}' for d in a_diff[:8]]}")
        print(f"  b_k ratios:      {[f'{r:.4f}' if r else 'None' for r in b_ratio[:8]]}")
    except Exception as e:
        print(f"  pattern check error: {e}")


if __name__ == "__main__":
    t0 = time.time()
    main()
    print(f"\nTotal time: {time.time() - t0:.2f} s")
