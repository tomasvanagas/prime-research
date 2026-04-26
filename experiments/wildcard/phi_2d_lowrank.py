"""
phi_2d_lowrank.py
==================

FRESH ANGLE (Session 65, fresh-perspective mode).

Question: viewed as a 2D matrix M[i,j] = phi(x_i, a_j), does the Meissel phi
function have low numerical rank or low effective sparsity in some basis?

If yes (rank r polylog in input), phi(x,a) — bottleneck of Meissel-Lehmer —
could in principle be reconstructed from O(r) entries via low-rank
interpolation, giving sub-x^{2/3} pi(x) computation.

We test framings:
  A:  phi(x_i, a_j) raw
  B:  phi - x*prod(1-1/p)  (residual from Mertens smooth)
  C:  phi(x,a)/x  (normalised)
  D:  phi(x,a) - phi(x,a-1)  (col-difference)

phi computed via direct inclusion-exclusion:
   phi(x, a) = sum_{d squarefree, all prime factors <= p_a} mu(d) * floor(x/d)
"""

import time
import numpy as np
from sympy import primerange


def primes_upto(n):
    return list(primerange(2, n + 1))


def phi_direct(x, a, primes):
    """phi(x, a) via inclusion-exclusion over smooth squarefree divisors.
    Only enumerate d <= x with all factors <= p_a."""
    if a == 0:
        return x
    if x == 0:
        return 0
    ps = primes[:a]
    total = x
    # DFS: pick subset of primes, multiply, with sign
    def dfs(idx, prod, sign):
        nonlocal total
        for j in range(idx, len(ps)):
            np_ = prod * ps[j]
            if np_ > x:
                break  # ps sorted, further only larger
            total += sign * (x // np_)
            dfs(j + 1, np_, -sign)
    dfs(0, 1, -1)
    return total


def phi_smooth(x, a, primes):
    prod = 1.0
    for p in primes[:a]:
        prod *= (1.0 - 1.0 / p)
    return x * prod


def build_matrix(framing, K, x_grid, a_grid, primes):
    M = np.zeros((K, K))
    if framing == "A":
        for i, x in enumerate(x_grid):
            for j, a in enumerate(a_grid):
                M[i, j] = phi_direct(x, a, primes)
    elif framing == "B":
        for i, x in enumerate(x_grid):
            for j, a in enumerate(a_grid):
                M[i, j] = phi_direct(x, a, primes) - phi_smooth(x, a, primes)
    elif framing == "C":
        for i, x in enumerate(x_grid):
            for j, a in enumerate(a_grid):
                if x > 0:
                    M[i, j] = phi_direct(x, a, primes) / x
    elif framing == "D":
        for i, x in enumerate(x_grid):
            prev = x  # phi(x, 0) = x
            for j, a in enumerate(a_grid):
                cur = phi_direct(x, a, primes)
                M[i, j] = cur - prev
                prev = cur
    return M


def analyze(M, name, log):
    U, S, Vt = np.linalg.svd(M, full_matrices=False)
    fro = np.linalg.norm(M)
    if S[0] == 0:
        log(f"\n[{name}] zero matrix, skip.")
        return None
    Sn = S / S[0]
    log(f"\n[{name}] shape {M.shape}  ||M||_F={fro:.4g}  sigma_1={S[0]:.4g}")
    log(f"  top 8 normalised: {[f'{v:.3e}' for v in Sn[:8]]}")
    eff = int((Sn > 1e-10).sum())
    log(f"  effective rank (>1e-10): {eff}")

    Kfit = min(20, len(S) - 1)
    if Kfit >= 3:
        ks = np.arange(1, Kfit + 1)
        ys = np.log(Sn[1:Kfit + 1] + 1e-30)
        b_exp, a_exp = np.polyfit(ks, ys, 1)
        alpha, _ = np.polyfit(np.log(ks), ys, 1)
        log(f"  exp fit:   sigma_k ~ exp({a_exp:.3g} {b_exp:+.3g}*k)")
        log(f"  power fit: sigma_k ~ k^{alpha:.3g}")
        if -b_exp > 0.3:
            verdict = "EXP-DECAY (low effective rank, COMPRESSIBLE)"
        elif alpha < -1.5:
            verdict = "FAST-POWER (compressible-ish)"
        elif alpha < -0.7:
            verdict = "MODERATE-POWER (mild structure)"
        else:
            verdict = "SLOW/FLAT (incompressible / pseudorandom-like)"
        log(f"  VERDICT: {verdict}")
        return {"sigma": Sn, "alpha": alpha, "b_exp": b_exp, "verdict": verdict}
    return {"sigma": Sn}


def main():
    out = []
    log = lambda s: (print(s), out.append(s))[1]

    log("=" * 60)
    log("PHI 2D LOW-RANK / SVD TEST  (fresh-perspective S65)")
    log("=" * 60)

    # Test at multiple scales to detect crossover from "small structure" to
    # "incompressible regime"
    scales = [
        ("small (K=18)", 18, 1.6, 8),
        ("medium (K=40)", 40, 1.35, 16),
        ("large (K=60)", 60, 1.25, 32),
    ]
    primes = primes_upto(20000)

    all_results = {}
    t0 = time.time()
    for scale_name, K, ratio, x0 in scales:
        x_grid = sorted(set(int(round(x0 * (ratio ** i))) for i in range(K)))
        while len(x_grid) < K:
            x_grid.append(x_grid[-1] + max(1, x_grid[-1] // 7))
        x_grid = x_grid[:K]
        a_grid = list(range(1, K + 1))
        log("\n" + "#" * 60)
        log(f"### SCALE {scale_name}: x_grid {x_grid[0]}..{x_grid[-1]}, a 1..{K}, p_K={primes[K-1]}")
        log("#" * 60)

        scale_res = {}
        for framing in ["A", "B", "C", "D"]:
            t1 = time.time()
            M = build_matrix(framing, K, x_grid, a_grid, primes)
            log(f"\n>> Framing {framing} built in {time.time()-t1:.2f}s")
            scale_res[framing] = analyze(M, f"{scale_name} | Framing {framing}", log)
        all_results[scale_name] = scale_res

    log("\n" + "=" * 60)
    log("CROSS-SCALE SUMMARY")
    log("=" * 60)
    log(f"{'scale':<18s} {'frame':<6s} {'b_exp':>8s} {'alpha':>8s}  verdict")
    for scale_name, res in all_results.items():
        for f, r in res.items():
            if r and 'b_exp' in r:
                log(f"{scale_name:<18s} {f:<6s} {r['b_exp']:>8.3f} {r['alpha']:>8.3f}  {r['verdict']}")
    log(f"\nTotal time: {time.time()-t0:.2f}s")

    log("""
INTERPRETATION:
  Framing A includes the rank-1 dominant smooth piece x*prod(1-1/p), so it
  trivially decays. The interesting question is the RESIDUAL spectrum in
  Framing B (phi minus the smooth Mertens approximation).

  If Framing B's spectrum decays EXPONENTIALLY, polylog reconstruction of
  phi(x, a) via low-rank interpolation is plausible -> direct attack on
  Meissel-Lehmer's O(x^{2/3}) bottleneck.

  If Framing B is flat / power-law with small exponent, the residual is
  numerically incompressible -- a structural obstruction.

  Framings C and D are diagnostic: C kills the trivial rank-1, D probes
  the column-incremental structure (related to phi(x/p_a, a-1) recursion).
""")
    return "\n".join(out)


if __name__ == "__main__":
    main()
