"""
Fast-multipole-style expansion of the Riemann explicit-formula zero sum.

Analogy: in N-body simulation, Σ q_j/(x - y_j) for many sources can be
evaluated in O(log(1/eps) * (M + N)) instead of O(MN) by clustering sources
and expanding in multipoles.

The zero-sum  S(x) = Σ_k x^{rho_k}/rho_k  has a similar structure:
  - Each γ_k is a "source" on the line Re(s) = 1/2
  - x is a "probe" parameter
  - The kernel  x^s/s = exp(s log x)/s

Cluster zeros into B groups by γ-value. Within each cluster, expand the
kernel around the cluster center γ_c by Taylor series to order P:
  x^{1/2 + iγ}/(1/2 + iγ)  ≈  (kernel at γ_c) * Σ_{p=0}^P (i(γ - γ_c))^p * D_p(γ_c, x)

If the per-cluster expansion converges to ε with order P = O(log 1/ε), then
the full sum costs O(B * P) instead of O(K). For K ≈ x/log x zeros, this is
O((K/B) * log(1/eps) + B * polylog) — minimized at B = sqrt(K), giving
O(sqrt(K) * polylog) — a SQRT speedup over direct summation.

Test:
- For x = 100, 1000, 10000 use first K=200 zeros.
- Cluster into B groups by γ. For each P ∈ {2, 5, 10, 20}, compute the
  multipole approximation. Compare to the exact partial sum.
- Plot error vs (B, P). If error ~ exp(-c*P) at fixed B, the analogy works.
"""

import math
import mpmath
from mpmath import mp, mpf, mpc

mp.dps = 30


def load_zeros(path, K):
    with open(path) as f:
        return [mpf(line.strip()) for line in f if line.strip()][:K]


def exact_zero_sum(x, gammas):
    """Sum over zeros (and conjugates) of x^rho/rho."""
    x_mp = mpf(x)
    log_x = mp.log(x_mp)
    sqrt_x = mp.sqrt(x_mp)
    s = mpf(0)
    for g in gammas:
        # 2 Re(x^rho / rho) where rho = 1/2 + i g
        # = 2 sqrt(x) * Re(exp(i g log x) / (1/2 + i g))
        c = mp.cos(g * log_x)
        si = mp.sin(g * log_x)
        denom = mpf("0.25") + g * g
        re_inv = mpf("0.5") / denom
        im_inv = -g / denom
        s += 2 * sqrt_x * (c * re_inv - si * im_inv)
    return s


def multipole_zero_sum(x, gammas, B, P):
    """Multipole approximation: cluster gammas into B contiguous groups,
    Taylor-expand kernel around cluster mean to order P."""
    x_mp = mpf(x)
    log_x = mp.log(x_mp)
    sqrt_x = mp.sqrt(x_mp)
    K = len(gammas)
    if B > K:
        B = K
    chunk = (K + B - 1) // B
    total = mpf(0)
    for b in range(B):
        cluster = gammas[b * chunk : (b + 1) * chunk]
        if not cluster:
            continue
        gc = sum(cluster) / len(cluster)
        # Kernel f(g) = exp(i g log x) / (1/2 + i g)
        # Taylor in delta = g - gc:  f(gc + delta) = sum_p f^{(p)}(gc)/p! * delta^p
        # f'(g) = (i log x * exp(i g log x))/(1/2+ig) - exp(i g log x)*i/(1/2+ig)^2
        # We compute coefficients f^{(p)}(gc)/p! by recursion on the product
        # exp(i g log x) and 1/(1/2 + i g) separately, then convolve.
        rho_c = mpc(mpf("0.5"), gc)
        # Series for exp(i g log x) = exp(i (gc + delta) log x) = exp(i gc log x) * exp(i delta log x)
        # coefficients: a_p = (i log x)^p / p!
        a_coefs = []
        a = mpc(1, 0)
        for p in range(P + 1):
            if p > 0:
                a = a * mpc(0, log_x) / p
            a_coefs.append(a)
        # rotation factor exp(i gc log x)
        rot = mpc(mp.cos(gc * log_x), mp.sin(gc * log_x))
        # Series for 1/(rho_c + i delta) = (1/rho_c) sum_p (-i delta / rho_c)^p
        # b_p = (-i)^p / rho_c^{p+1}
        b_coefs = []
        for p in range(P + 1):
            b_coefs.append(((mpc(0, -1)) ** p) / (rho_c ** (p + 1)))
        # Convolution: c_p = sum_{j=0..p} a_j * b_{p-j}, then multiply by rot
        c_coefs = []
        for p in range(P + 1):
            cp = mpc(0, 0)
            for j in range(p + 1):
                cp += a_coefs[j] * b_coefs[p - j]
            c_coefs.append(rot * cp)
        # Compute moments M_p = sum_{g in cluster} (g - gc)^p
        moments = []
        for p in range(P + 1):
            mom = mpf(0)
            for g in cluster:
                d = g - gc
                mom += d ** p
            moments.append(mom)
        # Cluster contribution: sum_p c_p * M_p (in complex)
        # We want 2 Re of (sqrt(x) * cluster contribution), since each rho contributes
        # x^rho/rho + x^rhobar/rhobar = 2 Re(x^rho/rho), and x^rho = sqrt(x) * exp(i g log x).
        contrib = mpc(0, 0)
        for p in range(P + 1):
            contrib += c_coefs[p] * moments[p]
        total += 2 * sqrt_x * contrib.real
    return total


def main():
    data = "/apps/aplikacijos/prime-research/data/zeta_zeros_1000.txt"
    K = 200
    gammas = load_zeros(data, K)
    print(f"Loaded {K} zeros, range γ ∈ [{float(gammas[0]):.2f}, {float(gammas[-1]):.2f}]")

    for x in [100, 1000, 10000]:
        print(f"\n{'='*60}\nx = {x}")
        exact = exact_zero_sum(x, gammas)
        print(f"  exact sum  = {mp.nstr(exact, 10)}")

        for B in [1, 2, 5, 10, 20, 50]:
            for P in [2, 5, 10, 15]:
                approx = multipole_zero_sum(x, gammas, B, P)
                err = float(abs(approx - exact))
                relerr = err / max(abs(float(exact)), 1e-30)
                print(f"  B={B:3d}, P={P:2d}: approx={mp.nstr(approx, 8):>20} err={err:.4e}  rel={relerr:.2e}")
            print()

        # Now check: at fixed B, how does error scale with P?
        print("  Convergence in P at B=10:")
        prev = None
        for P in range(1, 25):
            approx = multipole_zero_sum(x, gammas, 10, P)
            err = float(abs(approx - exact))
            ratio = "" if prev is None else f"  ratio={err/max(prev,1e-30):.3f}"
            print(f"    P={P:2d}: err={err:.4e}{ratio}")
            prev = err


if __name__ == "__main__":
    main()
