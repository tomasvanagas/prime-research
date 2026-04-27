"""
ATTACK_VECTORS §D.D30 — Pollicott-Ruelle resonances of χ_P-weighted
Gauss-map transfer operator.

Question: build the χ_P-weighted Gauss-map transfer operator
    (L_h f)(x) = Σ_{n=1}^{n_max} h(n) f(1/(x+n)) / (x+n)^2
for several arithmetic weights h ∈ {1, χ_P, λ, Λ} and matched-density
Bernoulli controls. Compute its top-K eigenvalues. Test:
  (i)  unweighted h=1 reproduces Mayer 1991 / Gauss-Kuzmin-Wirsing
       spectrum (λ_0=1, λ_1 ≈ -0.30366..., λ_2 ≈ 0.10088..., ...).
  (ii) χ_P-weighted spectral radius and gap vs matched-density
       Bernoulli controls — does a non-trivial isolated arithmetic
       resonance exist?
  (iii) Stability under discretisation refinement (M_grid, n_max) —
        true Pollicott-Ruelle resonances are stable; spurious are not.

Cross-domain technique: Pollicott-Ruelle resonance theory + Mayer 1991
arithmetic transfer-operator approach to ζ via Gauss map. Channelled
mathematician: Ruelle / Baladi (transfer operator spectral theory),
Mayer (arithmetic transfer operators).

Pre-stated falsification criteria (results.md):
  F-A  unweighted h=1: top-5 eigenvalues match published GKW spectrum
       (λ_0=1.0 ± 1e-6, λ_1 = -0.30366300289... ± 1e-4) — sanity check.
  F-B  χ_P spectral radius/gap match matched-density Bernoulli mean
       within 2σ across 30 controls — closure mode E (discretisation-
       stable spectral bound is the new B-grade content).
  F-C  χ_P top eigenvalue or first sub-leading eigenvalue deviates from
       Bernoulli baseline by ≥ 3σ on at least one (M_grid, n_max)
       AND deviation persists under 2× refinement of M_grid AND n_max
       — claim mode I, B-grade negative-shape edge candidate.
  F-D  χ_P resonance has closed-form interpretation tying it to a known
       arithmetic invariant (Mertens product, prime sum, etc.) — A-grade.

Author: session 140.
"""

import argparse
import json
import os
import sys
import time
from pathlib import Path

import numpy as np
from sympy import isprime, primerange, primepi, factorint


def chebyshev_lobatto(M):
    """Chebyshev-Lobatto nodes on [0,1] and barycentric weights."""
    k = np.arange(M + 1)
    x = 0.5 * (1 - np.cos(np.pi * k / M))
    w = np.ones(M + 1)
    w[1::2] = -1.0
    w[0] *= 0.5
    w[-1] *= 0.5
    return x, w


def barycentric_weights(y, x_nodes, w_bary):
    """
    Lagrange barycentric basis coefficients at point y for nodes x_nodes
    with given weights w_bary. Returns array c of length len(x_nodes)
    with f(y) = sum_j c[j] * f(x_nodes[j]).
    """
    diffs = y - x_nodes
    # Exact-match check
    if np.any(np.abs(diffs) < 1e-15):
        j = int(np.argmin(np.abs(diffs)))
        out = np.zeros_like(x_nodes)
        out[j] = 1.0
        return out
    terms = w_bary / diffs
    return terms / terms.sum()


def transfer_operator_matrix(weights, M_grid=80, n_max=200):
    """
    Build the matrix of the Gauss-map transfer operator
       (L_h f)(x) = Σ_{n=1}^{n_max} h(n) f(1/(x+n)) / (x+n)^2
    on a Chebyshev-Lobatto grid of M_grid+1 nodes via barycentric
    interpolation.

    weights : 1D array, weights[n] for n=0..n_max (weights[0] ignored).
    """
    x_nodes, w_bary = chebyshev_lobatto(M_grid)
    N = M_grid + 1
    L = np.zeros((N, N), dtype=np.float64)

    for i, x in enumerate(x_nodes):
        for n in range(1, n_max + 1):
            wn = weights[n]
            if wn == 0:
                continue
            y = 1.0 / (x + n)
            jac = 1.0 / (x + n) ** 2
            c = barycentric_weights(y, x_nodes, w_bary)
            L[i, :] += wn * jac * c
    return L, x_nodes


def top_eigenvalues(L, K=30):
    """Top-K eigenvalues by magnitude. Returns sorted array (complex)."""
    e = np.linalg.eigvals(L)
    order = np.argsort(-np.abs(e))
    return e[order[:K]]


def eigenpairs_top(L, K=10):
    """Top-K (eigenvalue, eigenvector) pairs by |λ|."""
    e, V = np.linalg.eig(L)
    order = np.argsort(-np.abs(e))
    return e[order[:K]], V[:, order[:K]]


def chi_P_weights(n_max):
    w = np.zeros(n_max + 1, dtype=np.float64)
    for p in primerange(2, n_max + 1):
        w[p] = 1.0
    return w


def lambda_weights(n_max):
    """Liouville function λ(n) = (-1)^Ω(n)."""
    w = np.zeros(n_max + 1, dtype=np.float64)
    for n in range(1, n_max + 1):
        Omega_n = sum(factorint(n).values())
        w[n] = (-1) ** Omega_n
    return w


def Lambda_weights(n_max):
    """von Mangoldt Λ(n) = log p if n=p^k else 0."""
    w = np.zeros(n_max + 1, dtype=np.float64)
    for n in range(2, n_max + 1):
        f = factorint(n)
        if len(f) == 1:
            p = next(iter(f.keys()))
            w[n] = float(np.log(p))
    return w


def bernoulli_weights(n_max, density, rng):
    """{0,1} Bernoulli weights, density-matched to χ_P (NAIVE baseline)."""
    return (rng.random(n_max + 1) < density).astype(np.float64)


def support_matched_weights(n_max, k, rng, support_lo=2):
    """
    Cardinality-matched random subset of size k drawn uniformly from
    [support_lo, n_max]. Forces w[0] = w[1] = ... = w[support_lo-1] = 0,
    matching χ_P's exact lack of weight at n=1 (and at n=0).
    """
    pool = np.arange(support_lo, n_max + 1)
    sel = rng.choice(pool, size=k, replace=False)
    w = np.zeros(n_max + 1, dtype=np.float64)
    w[sel] = 1.0
    return w


def odd_parity_matched_weights(n_max, k_odd, has_two, rng):
    """
    Random subset matching χ_P's parity structure: exactly 1 even (n=2)
    if has_two else 0 evens, plus k_odd odd indices ≥ 3.
    """
    odds = np.arange(3, n_max + 1, 2)
    sel_odd = rng.choice(odds, size=k_odd, replace=False)
    w = np.zeros(n_max + 1, dtype=np.float64)
    w[sel_odd] = 1.0
    if has_two:
        w[2] = 1.0
    return w


def random_signed_weights(n_max, rng):
    """{-1, +1} Rademacher weights — Liouville baseline control."""
    w = (rng.integers(0, 2, size=n_max + 1) * 2 - 1).astype(np.float64)
    w[0] = 0
    return w


def random_signed_weights_match_lam_at_1(n_max, rng):
    """
    Rademacher ±1 with w[1] = +1 forced (matches λ(1) = +1, the dominant
    cylinder kernel). Eliminates the noise of the random sign at n=1.
    """
    w = random_signed_weights(n_max, rng)
    w[1] = 1.0
    return w


def cramer_weights(n_max, rng, support_lo=2):
    """
    Cramér model: w[n] ~ Bernoulli(1/log n) for n ≥ support_lo.
    The natural prime-density-matched random control with the same
    log-density profile as χ_P.
    """
    probs = np.zeros(n_max + 1, dtype=np.float64)
    for n in range(support_lo, n_max + 1):
        probs[n] = 1.0 / max(np.log(n), 1.0)
    w = (rng.random(n_max + 1) < probs).astype(np.float64)
    w[:support_lo] = 0.0
    return w


def cramer_odd_weights(n_max, rng, support_lo=3):
    """
    Cramér model restricted to odd integers ≥ 3, matching χ_P's parity
    structure: only one even prime (n=2), all others odd. Includes
    n=2 with probability 1.
    """
    probs = np.zeros(n_max + 1, dtype=np.float64)
    for n in range(support_lo, n_max + 1, 2):
        probs[n] = 2.0 / max(np.log(n), 1.0)  # 2× since restricted to odds
    w = (rng.random(n_max + 1) < probs).astype(np.float64)
    w[:support_lo] = 0.0
    w[2] = 1.0  # always include the unique even prime
    # zero out evens beyond 2
    w[4::2] = 0.0
    return w


def density_of(weights):
    return float(np.mean(weights[1:] != 0))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--M_grid", type=int, default=80,
                    help="Chebyshev grid size (matrix dim = M_grid+1)")
    ap.add_argument("--n_max", type=int, default=200,
                    help="Cylinder truncation for transfer operator sum")
    ap.add_argument("--n_seeds", type=int, default=30,
                    help="Number of Bernoulli matched-density controls")
    ap.add_argument("--K_eig", type=int, default=30,
                    help="Top-K eigenvalues to retain")
    ap.add_argument("--out", type=str, default="results.json")
    ap.add_argument("--include_signed", action="store_true",
                    help="Also include Rademacher controls (for λ baseline)")
    args = ap.parse_args()

    out_path = Path(__file__).parent / args.out

    print(f"[D30] Pollicott-Ruelle resonances of χ_P-weighted Gauss-map "
          f"transfer operator")
    print(f"  M_grid = {args.M_grid}, n_max = {args.n_max}, "
          f"n_seeds = {args.n_seeds}")

    rng = np.random.default_rng(20260427)

    # === Build weights ===
    chi_P = chi_P_weights(args.n_max)
    lam = lambda_weights(args.n_max)
    Lam = Lambda_weights(args.n_max)
    unw = np.zeros(args.n_max + 1)
    unw[1:] = 1.0  # h=1 unweighted Gauss map
    rho_chi = density_of(chi_P)
    rho_lam = density_of(lam)  # nominally 1.0 for Liouville
    print(f"  density(χ_P up to {args.n_max}) = {rho_chi:.5f}  "
          f"(π({args.n_max}) = {int(primepi(args.n_max))})")

    # === Compute eigenvalues for the named operators ===
    cases = {}
    t0 = time.time()
    L_unw, x_nodes = transfer_operator_matrix(unw, args.M_grid, args.n_max)
    e_unw = top_eigenvalues(L_unw, args.K_eig)
    cases["unweighted"] = e_unw
    print(f"[unweighted] top-5 eigenvalues by |λ|:")
    for k, ev in enumerate(e_unw[:5]):
        print(f"   λ_{k} = {ev.real:+.10f}{ev.imag:+.4e}j  "
              f"|λ| = {abs(ev):.10f}")
    print(f"  elapsed: {time.time()-t0:.2f}s")

    t0 = time.time()
    L_chi, _ = transfer_operator_matrix(chi_P, args.M_grid, args.n_max)
    e_chi = top_eigenvalues(L_chi, args.K_eig)
    cases["chi_P"] = e_chi
    print(f"[χ_P] top-5 eigenvalues by |λ|:")
    for k, ev in enumerate(e_chi[:5]):
        print(f"   λ_{k} = {ev.real:+.10f}{ev.imag:+.4e}j  "
              f"|λ| = {abs(ev):.10f}")
    print(f"  elapsed: {time.time()-t0:.2f}s")

    t0 = time.time()
    L_lam, _ = transfer_operator_matrix(lam, args.M_grid, args.n_max)
    e_lam = top_eigenvalues(L_lam, args.K_eig)
    cases["lambda"] = e_lam
    print(f"[λ] top-5 eigenvalues by |λ|:")
    for k, ev in enumerate(e_lam[:5]):
        print(f"   λ_{k} = {ev.real:+.10f}{ev.imag:+.4e}j  "
              f"|λ| = {abs(ev):.10f}")

    t0 = time.time()
    L_Lam, _ = transfer_operator_matrix(Lam, args.M_grid, args.n_max)
    e_Lam = top_eigenvalues(L_Lam, args.K_eig)
    cases["vonMangoldt"] = e_Lam
    print(f"[Λ] top-5 eigenvalues by |λ|:")
    for k, ev in enumerate(e_Lam[:5]):
        print(f"   λ_{k} = {ev.real:+.10f}{ev.imag:+.4e}j  "
              f"|λ| = {abs(ev):.10f}")

    # === Three baselines: Bernoulli (naive) / support-matched / parity-matched ===
    n_primes = int(np.sum(chi_P[1:] != 0))
    has_two = bool(chi_P[2] != 0)
    n_odd_primes = int(np.sum(chi_P[3::2] != 0))
    print(f"\n[control baselines] n_primes = {n_primes}, "
          f"has_two = {has_two}, n_odd_primes = {n_odd_primes}")

    def run_baseline_ensemble(label, sampler):
        print(f"\n[{label}] {args.n_seeds} controls...")
        t0_ = time.time()
        eigs = np.zeros((args.n_seeds, args.K_eig), dtype=np.complex128)
        for s in range(args.n_seeds):
            wb = sampler()
            L_b, _ = transfer_operator_matrix(wb, args.M_grid, args.n_max)
            eigs[s] = top_eigenvalues(L_b, args.K_eig)
        print(f"  elapsed: {time.time()-t0_:.2f}s")
        return eigs

    bern_eigs = run_baseline_ensemble(
        "B_naive (Bernoulli rho)",
        lambda: bernoulli_weights(args.n_max, rho_chi, rng))
    supp_eigs = run_baseline_ensemble(
        "B_supp (support [2..n_max], cardinality-matched)",
        lambda: support_matched_weights(args.n_max, n_primes, rng))
    par_eigs = run_baseline_ensemble(
        "B_par (parity-matched: 1 even=2 + odd primes count)",
        lambda: odd_parity_matched_weights(args.n_max, n_odd_primes, has_two, rng))
    cra_eigs = run_baseline_ensemble(
        "B_cra (Cramér model 1/log n)",
        lambda: cramer_weights(args.n_max, rng))
    crao_eigs = run_baseline_ensemble(
        "B_crao (Cramér + odd parity, n=2 fixed)",
        lambda: cramer_odd_weights(args.n_max, rng))

    def stats(label, control_eigs):
        ctop1 = np.abs(control_eigs[:, 0])
        ctop2 = np.abs(control_eigs[:, 1])
        cgap = ctop2 / np.maximum(ctop1, 1e-12)
        chi_top1 = abs(e_chi[0])
        chi_top2 = abs(e_chi[1])
        chi_gap = chi_top2 / max(chi_top1, 1e-12)
        z_top1 = (chi_top1 - ctop1.mean()) / max(ctop1.std(), 1e-12)
        z_top2 = (chi_top2 - ctop2.mean()) / max(ctop2.std(), 1e-12)
        z_gap = (chi_gap - cgap.mean()) / max(cgap.std(), 1e-12)
        print(f"\n[stats: χ_P vs {label}]")
        print(f"  χ_P |λ_0| = {chi_top1:.6f};  ctrl mean ± std = "
              f"{ctop1.mean():.6f} ± {ctop1.std():.6f};  z = {z_top1:+.3f}")
        print(f"  χ_P |λ_1| = {chi_top2:.6f};  ctrl mean ± std = "
              f"{ctop2.mean():.6f} ± {ctop2.std():.6f};  z = {z_top2:+.3f}")
        print(f"  χ_P gap   = {chi_gap:.6f};  ctrl mean ± std = "
              f"{cgap.mean():.6f} ± {cgap.std():.6f};  z = {z_gap:+.3f}")
        return {
            "chi_top1_abs": float(chi_top1),
            "chi_top2_abs": float(chi_top2),
            "chi_gap": float(chi_gap),
            "ctrl_top1_mean": float(ctop1.mean()),
            "ctrl_top1_std": float(ctop1.std()),
            "ctrl_top2_mean": float(ctop2.mean()),
            "ctrl_top2_std": float(ctop2.std()),
            "ctrl_gap_mean": float(cgap.mean()),
            "ctrl_gap_std": float(cgap.std()),
            "z_top1": float(z_top1),
            "z_top2": float(z_top2),
            "z_gap": float(z_gap),
        }

    stat_naive = stats("B_naive", bern_eigs)
    stat_supp = stats("B_supp", supp_eigs)
    stat_par = stats("B_par", par_eigs)
    stat_cra = stats("B_cra (Cramér)", cra_eigs)
    stat_crao = stats("B_crao (Cramér+odd)", crao_eigs)

    # === Liouville-baseline: λ vs Rademacher ±1 controls ===
    rade_eigs = run_baseline_ensemble(
        "Rademacher (±1)",
        lambda: random_signed_weights(args.n_max, rng))
    rade1_eigs = run_baseline_ensemble(
        "Rademacher_w1=1 (matches λ(1)=+1)",
        lambda: random_signed_weights_match_lam_at_1(args.n_max, rng))

    def lam_stats(label, ctrl_eigs):
        ctop1 = np.abs(ctrl_eigs[:, 0])
        lam_top1 = abs(e_lam[0])
        z = (lam_top1 - ctop1.mean()) / max(ctop1.std(), 1e-12)
        print(f"[stats: λ vs {label}]")
        print(f"  λ |λ_0| = {lam_top1:.6f};  ctrl mean ± std = "
              f"{ctop1.mean():.6f} ± {ctop1.std():.6f};  z = {z:+.3f}")
        return {
            "lam_top1_abs": float(lam_top1),
            "ctrl_top1_mean": float(ctop1.mean()),
            "ctrl_top1_std": float(ctop1.std()),
            "z_lam_top1": float(z),
        }

    print(f"\n[Liouville comparisons]")
    stat_lam_rade = lam_stats("Rademacher", rade_eigs)
    stat_lam_rade1 = lam_stats("Rademacher_w1=1", rade1_eigs)

    # === Save artefacts ===
    record = {
        "M_grid": args.M_grid,
        "n_max": args.n_max,
        "n_seeds": args.n_seeds,
        "K_eig": args.K_eig,
        "rho_chi_P": rho_chi,
        "primes_le_nmax": int(primepi(args.n_max)),
        "n_primes": n_primes,
        "has_two": has_two,
        "n_odd_primes": n_odd_primes,
        "cases": {
            k: {"real": [float(x.real) for x in v],
                "imag": [float(x.imag) for x in v]}
            for k, v in cases.items()
        },
        "ensemble_top_eigs": {
            "B_naive": {"real": bern_eigs.real.tolist(),
                         "imag": bern_eigs.imag.tolist()},
            "B_supp": {"real": supp_eigs.real.tolist(),
                       "imag": supp_eigs.imag.tolist()},
            "B_par": {"real": par_eigs.real.tolist(),
                      "imag": par_eigs.imag.tolist()},
            "B_cra": {"real": cra_eigs.real.tolist(),
                      "imag": cra_eigs.imag.tolist()},
            "B_crao": {"real": crao_eigs.real.tolist(),
                       "imag": crao_eigs.imag.tolist()},
            "Rademacher": {"real": rade_eigs.real.tolist(),
                           "imag": rade_eigs.imag.tolist()},
            "Rademacher_w1": {"real": rade1_eigs.real.tolist(),
                               "imag": rade1_eigs.imag.tolist()},
        },
        "stats": {
            "B_naive": stat_naive,
            "B_supp": stat_supp,
            "B_par": stat_par,
            "B_cra": stat_cra,
            "B_crao": stat_crao,
            "lambda_vs_Rademacher": stat_lam_rade,
            "lambda_vs_Rademacher_w1": stat_lam_rade1,
        },
    }
    with open(out_path, "w") as fp:
        json.dump(record, fp, indent=2)
    print(f"\n[saved] {out_path}")


if __name__ == "__main__":
    main()
