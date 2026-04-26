"""
Wirtinger-flow recovery of zeta zero ordinates from π(x) values.

Reverse the usual question: instead of "given zeros, compute π(x)", ask:
"given π(x) at moderate x, can we *recover* the first K zeros γ_1, ..., γ_K
via nonconvex optimization?"

If recovery is sub-exponential in K, then by duality we get a fast π-counter
for large x: precompute K = polylog(x) zeros once, then evaluate explicit
formula in O(K).

The objective is roughly:
    L(γ) = Σ_{x in S}  (π_explicit(x; γ) − π_observed(x))^2

This is highly nonconvex (each γ_k contributes a sinusoid in log x). We use
gradient descent (Wirtinger flow's real analog) from a "warm" initialization
and check if the recovered γ̂ approaches the true γ.

Test:
- True γ_1, ..., γ_K from data file (K ∈ {3, 5, 10}).
- Generate π̂(x) = R(x) − Σ x^ρ/ρ + ... (oracle) at S = {x_1, ..., x_M}.
- Initialize γ̂ at perturbed truth (test convergence basin width) AND from random.
- Run gradient descent on L(γ̂).
- Report: final loss, distance to truth, basin radius.

If the loss landscape has many bad local minima even from near-true init,
the inverse problem is hard → recovery from π values is at least as hard as
direct computation.
"""

import numpy as np
import math
import os

# We work in float for speed. Higher precision not needed for landscape diagnostic.

def explicit_psi_from_gammas(x, gammas):
    """psi(x) ≈ x − Σ 2 sqrt(x) Re(exp(i γ log x)/(1/2 + i γ)) − log(2π) − ½ log(1−x⁻²)"""
    log_x = np.log(x)
    sqrt_x = np.sqrt(x)
    s = 0.0
    for g in gammas:
        c = np.cos(g * log_x)
        si = np.sin(g * log_x)
        denom = 0.25 + g * g
        re_inv = 0.5 / denom
        im_inv = -g / denom
        s += 2 * sqrt_x * (c * re_inv - si * im_inv)
    return x - s - np.log(2 * np.pi) - 0.5 * np.log(1 - 1 / (x * x))


def loss(gammas_hat, xs, observed_psi):
    preds = np.array([explicit_psi_from_gammas(x, gammas_hat) for x in xs])
    return np.mean((preds - observed_psi) ** 2)


def grad_loss(gammas_hat, xs, observed_psi, h=1e-6):
    """Numerical gradient (avoid hand-derivation)."""
    g = np.zeros(len(gammas_hat))
    base = loss(gammas_hat, xs, observed_psi)
    for i in range(len(gammas_hat)):
        gp = gammas_hat.copy()
        gp[i] += h
        g[i] = (loss(gp, xs, observed_psi) - base) / h
    return g


def gd_recover(init, xs, observed, lr=0.01, max_iter=200, verbose=False):
    g = np.array(init, dtype=float)
    L_history = []
    for it in range(max_iter):
        L = loss(g, xs, observed)
        L_history.append(L)
        gr = grad_loss(g, xs, observed)
        # gradient clip
        gnorm = np.linalg.norm(gr)
        if gnorm > 100:
            gr = gr / gnorm * 100
        g_new = g - lr * gr
        # Backtrack if loss increases
        if loss(g_new, xs, observed) > L:
            lr *= 0.5
            if lr < 1e-8:
                break
            continue
        g = g_new
        if verbose and it % 20 == 0:
            print(f"    iter {it}: loss = {L:.6e}")
        if L < 1e-20:
            break
    return g, L_history


def main():
    data = "/apps/aplikacijos/prime-research/data/zeta_zeros_1000.txt"
    with open(data) as f:
        true_gammas_all = [float(l) for l in f if l.strip()]

    # Test x values
    xs = np.array([20, 30, 50, 70, 100, 150, 200, 300, 500, 700, 1000])

    for K in [3, 5, 10]:
        true_gammas = true_gammas_all[:K]
        observed = np.array([explicit_psi_from_gammas(x, true_gammas) for x in xs])
        print("\n" + "=" * 60)
        print(f"K = {K}, true γ = {[round(g, 3) for g in true_gammas]}")

        # Test 1: init perturbed by epsilon
        for eps in [0.01, 0.1, 0.5, 1.0, 2.0]:
            np.random.seed(42)
            init = np.array(true_gammas) + eps * np.random.randn(K)
            g_hat, hist = gd_recover(init, xs, observed, lr=0.01, max_iter=200)
            dist = np.linalg.norm(g_hat - np.array(true_gammas))
            final_loss = hist[-1] if hist else float("inf")
            print(f"  init perturb ε={eps}: dist={dist:.4f}, final_loss={final_loss:.3e}")

        # Test 2: init random in correct range
        print("  Random initialization tests:")
        successes = 0
        trials = 5
        for seed in range(trials):
            np.random.seed(seed)
            # init in [γ_1 - 5, γ_K + 5]
            init = np.random.uniform(true_gammas[0] - 5, true_gammas[-1] + 5, K)
            init.sort()
            g_hat, hist = gd_recover(init, xs, observed, lr=0.005, max_iter=300)
            g_hat_sorted = np.sort(g_hat)
            dist = np.linalg.norm(g_hat_sorted - np.array(sorted(true_gammas)))
            converged = dist < 0.1
            if converged:
                successes += 1
            print(f"    seed={seed}: init={[round(g,2) for g in init]}, recovered={[round(g,2) for g in g_hat_sorted]}, dist={dist:.3f}, conv={converged}")
        print(f"  K={K} random init success rate: {successes}/{trials}")


if __name__ == "__main__":
    main()
