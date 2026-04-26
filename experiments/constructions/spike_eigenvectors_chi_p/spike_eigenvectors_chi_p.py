#!/usr/bin/env python3
"""S82 / C2 sub-arc — Characterization of spike singular vectors of chi_P MPS unfolding.

S74 (free_cumulants_chi_p) found that the chi_P MPS unfolding M^(j) of
shape (W^j, W^{d-j}) has a three-component spectrum: a finite structural
peak, a band of O(N^{0.42}) "spike" eigenvalues absent from the matched
random baseline, and an MP bulk at rate c = phi(W)/W.

The spike COUNT is empirically polynomial; the next question is what the
spike EIGENVECTORS encode arithmetically. Three candidate origins:
  (a) Residue-class biases pi(N; q, a) for small q;
  (b) Selberg / explicit-formula low-lying zeta-zero contributions;
  (c) Cross-block correlation modes (i.e., pure linear-algebraic noise).

This script tests (a) directly. For each top singular pair (u_k, v_k)
of M, reconstruct f_k = sigma_k * outer(u_k, v_k) flattened to length N,
then for each modulus q in q_list compute the L^2 fraction of f_k that
sits in the subspace of "functions of n mod q". A small Q* such that the
spikes are entirely absorbed by V_{<= Q*} would identify the spike
subspace as the residue-class character subspace; a "diffuse" pattern
(no q captures > 50% of f_k for k > 1) would falsify (a) and push us
toward (b) or (c).

Pre-stated falsifiers:
  PR1: For W=2 d=20 j=10, EACH of the top k_*=26 spike eigenvectors
       (excluding the leading rank-1 mean component) achieves
       max_q (residue-energy fraction) > 0.5 for some q in {2..200}.
  PR2: For the matched-active-Bernoulli baseline at the same shape, NO
       singular vector beyond k=1 achieves max_q (residue-energy
       fraction) > 0.2 for any q in {2..200}.
  PR3: The total dimension of residue-class subspaces visited by spike
       eigenvectors satisfies sum_{q in visited} (q-1) ~ k_*(N) within
       a factor of 4. (Order-of-magnitude check on the dimension count.)

If PR1 fails (e.g., spikes are diffuse over residue classes), the
spike origin is (b) or (c) and a different characterization is required.

Cross-domain reference:
  Iwaniec & Kowalski, "Analytic Number Theory" (AMS Colloq. 53), Ch. 3
    on Dirichlet characters and orthogonality of residue-class projections.
  Mingo & Speicher 2017 (Free Probability and Random Matrices), Ch. 1
    on free-cumulant standardization (carry-over from S74).
"""

import argparse
import json
import time
from pathlib import Path

import numpy as np


def sieve(N):
    is_p = np.ones(N + 1, dtype=bool)
    is_p[:2] = False
    for i in range(2, int(N**0.5) + 1):
        if is_p[i]:
            is_p[i * i :: i] = False
    return is_p


def chi_P_vec(N):
    return sieve(N)[1 : N + 1].astype(np.float64)


def euler_phi(n):
    r = n
    p = 2
    nn = n
    while p * p <= nn:
        if nn % p == 0:
            while nn % p == 0:
                nn //= p
            r -= r // p
        p += 1
    if nn > 1:
        r -= r // nn
    return r


def precompute_residues(N, q_list):
    """For each q in q_list, return array residues_q of length N where
    residues_q[i] = (i+1) mod q (so n = 1..N, residue in 0..q-1)."""
    n_arr = np.arange(1, N + 1, dtype=np.int64)
    return {q: (n_arr % q).astype(np.int32) for q in q_list}


def residue_energy_fraction(f, residues, q):
    """Fraction of L^2 energy of f that sits in the 'function of n mod q' subspace.

    Computes ||P_q f||^2 / ||f||^2 where (P_q f)(n) = mean of f over residue
    class containing n.
    """
    f_sq = float(np.dot(f, f))
    if f_sq <= 0:
        return 0.0
    class_sums = np.bincount(residues, weights=f, minlength=q)
    class_counts = np.bincount(residues, minlength=q)
    nz = class_counts > 0
    means = np.zeros_like(class_sums, dtype=np.float64)
    means[nz] = class_sums[nz] / class_counts[nz]
    energy = float(np.sum(class_counts * means**2))
    return energy / f_sq


def residue_energy_fraction_nonprincipal(f, residues, q):
    """Fraction in residue-mod-q subspace MINUS the constant component.

    Subtracts the global mean before computing the residue projection so
    that the principal character (constant function on residues) doesn't
    contribute. This isolates the genuinely q-arithmetic content.
    """
    f_centered = f - f.mean()
    return residue_energy_fraction(f_centered, residues, q)


def analyze_spike_vectors(M, N, q_list, K_top, label="chi_P"):
    """Compute the top K_top singular vectors of M and their residue energies.

    Returns a list of K_top dicts, each with sigma, residue energies per q,
    dominant q, minimum-conductor q (smallest q whose centered residue
    energy is within 5% of the maximum across q_list), and small-prime
    energy decomposition.
    """
    print(f"  [{label}] SVD of {M.shape}...")
    t0 = time.time()
    U, S, Vt = np.linalg.svd(M, full_matrices=False)
    elapsed_svd = time.time() - t0
    print(f"  [{label}] SVD done in {elapsed_svd:.2f}s; top-30 sv = {S[:30].round(3).tolist()}")

    K = min(K_top, len(S))
    residues_map = precompute_residues(N, q_list)

    # Specifically track these "small-prime / small-prime-power" moduli.
    small_q = [q for q in [2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16,
                            18, 20, 21, 24, 25, 27, 28, 30, 32, 33, 35, 36, 42,
                            45, 48, 49, 50, 56, 60, 64, 72, 84, 90, 96, 105,
                            120, 144, 168, 180, 210] if q in q_list]

    out = []
    print(f"  [{label}] computing residue energies for top-{K} sv pairs...")
    t0 = time.time()
    for k in range(K):
        sigma_k = float(S[k])
        # Reconstruct rank-1 contribution as length-N vector.
        f_k = sigma_k * np.outer(U[:, k], Vt[k, :]).reshape(-1)
        energies_centered = {}
        for q in q_list:
            energies_centered[q] = residue_energy_fraction_nonprincipal(
                f_k, residues_map[q], q
            )
        # Identify dominant q (by centered energy).
        cand = [(q, e) for q, e in energies_centered.items() if q > 1]
        cand.sort(key=lambda x: x[1], reverse=True)
        dom_q, dom_e = cand[0] if cand else (None, 0.0)
        top5 = cand[:5]

        # Find the SMALLEST q whose centered energy is within 5% (relative)
        # of the maximum dom_e. This is the "conductor" of the spike: the
        # smallest residue-class subspace that captures the spike's content.
        min_q_conductor = None
        if dom_e > 0:
            threshold = 0.95 * dom_e
            # Sort q in increasing order, find first that crosses threshold.
            for q in sorted(energies_centered.keys()):
                if energies_centered[q] >= threshold:
                    min_q_conductor = q
                    break

        # Small-prime / small-q energy panel (for diagnostic / pretty-print).
        small_q_energy = {q: float(energies_centered[q]) for q in small_q}

        out.append({
            "k": k,
            "sigma": sigma_k,
            "residue_energy_fraction_centered": {str(q): float(e) for q, e in energies_centered.items()},
            "dom_q": int(dom_q) if dom_q is not None else None,
            "dom_q_centered_energy": float(dom_e),
            "min_q_conductor": int(min_q_conductor) if min_q_conductor is not None else None,
            "top5_q": [(int(q), float(e)) for q, e in top5],
            "small_q_energy": small_q_energy,
        })
    elapsed_proj = time.time() - t0
    print(f"  [{label}] projections done in {elapsed_proj:.2f}s")
    return out, S


def make_chi_p_matrix(W, d, j=None):
    if j is None:
        j = d // 2
    N = W**d
    chi = chi_P_vec(N)
    M = chi.reshape(W**j, W ** (d - j))
    return M, N, j


def make_active_bernoulli_matrix(W, d, j=None, density=None, seed=0):
    """Random Bernoulli baseline matched to the active-submatrix shape.

    Mimics the S74 matched baseline. Note: this is NOT the full
    (W^j, W^{d-j}) matrix; it's the active sub-block of shape
    (W^j - 1, phi(W)*W^{d-j-1}). For comparison-of-spectrum-shape with
    chi_P, we embed it back into (W^j, W^{d-j}) by zero-padding so
    indices line up.
    """
    if j is None:
        j = d // 2
    N = W**d
    phiW = euler_phi(W)
    m_act = W**j
    n_act = W ** (d - j)
    n_live = phiW * W ** (d - j - 1)

    if density is None:
        # Use the actual chi_P density on the active submatrix.
        chi = chi_P_vec(N)
        chi_M = chi.reshape(W**j, W ** (d - j))
        # Active columns: those with (j' + 1) coprime to W.
        col_idx = np.arange(W ** (d - j))
        n_at_col = col_idx + 1
        gcd = np.gcd(n_at_col, W * np.ones_like(n_at_col))
        live_col_mask = (gcd == 1)
        n_live_actual = int(live_col_mask.sum())
        # Density on live columns.
        density = float(chi_M[:, live_col_mask].mean())
    else:
        live_col_mask = None

    rng = np.random.default_rng(seed)
    # Build full (m_act, n_act) baseline:
    # - on live columns: iid Bernoulli(density)
    # - on dead columns: zeros (matching chi_P's structural zeros).
    M = np.zeros((m_act, n_act), dtype=np.float64)
    if live_col_mask is None:
        col_idx = np.arange(n_act)
        n_at_col = col_idx + 1
        gcd = np.gcd(n_at_col, W * np.ones_like(n_at_col))
        live_col_mask = (gcd == 1)
    M[:, live_col_mask] = rng.binomial(1, density, size=(m_act, int(live_col_mask.sum()))).astype(np.float64)
    return M, N, j


def fmt_e(x, prec=4):
    return f"{x:.{prec}f}"


def summarize_run(spike_data, k_star, label, show_K=None):
    """Print a compact summary of dominant q + conductor q for top spike eigenvectors."""
    if show_K is None:
        show_K = min(max(k_star + 5, 20), len(spike_data))
    else:
        show_K = min(show_K, len(spike_data))
    print(f"\n  [{label}] Top {show_K} singular vectors (* = within k_star={k_star}):")
    print(f"  {'k':>3}  {'sigma':>10}  {'min_q_cond':>10}  {'dom_q':>6}  {'dom_e':>8}  top3 (q,e)")
    for s in spike_data[:show_K]:
        top5 = s["top5_q"][:3]
        top3_str = "  ".join(f"({q},{e:.2f})" for q, e in top5)
        marker = " *" if s["k"] < k_star else ""
        cond = s.get("min_q_conductor", "?")
        print(f"  {s['k']:>3}  {s['sigma']:>10.3f}  {str(cond):>10}  "
              f"{s['dom_q']:>6}  {s['dom_q_centered_energy']:>8.3f}  {top3_str}{marker}")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--W", type=int, default=2)
    parser.add_argument("--d", type=int, default=20)
    parser.add_argument("--j", type=int, default=None,
                        help="cut position; defaults to d/2")
    parser.add_argument("--K_top", type=int, default=60,
                        help="number of top singular vectors to analyze")
    parser.add_argument("--k_star", type=int, default=26,
                        help="empirical spike count from S74 (for marking)")
    parser.add_argument("--q_max", type=int, default=200)
    parser.add_argument("--with_baseline", action="store_true",
                        help="also analyze the matched active-Bernoulli baseline")
    parser.add_argument("--also", type=str, default="",
                        help="comma-separated additional configs as W:d (e.g. '2:18,3:12')")
    parser.add_argument("--out", type=str, default="spike_eigenvectors_chi_p_results.json")
    args = parser.parse_args()

    # Build q_list: 2..q_max (all integers).
    q_list = list(range(2, args.q_max + 1))

    configs = [(args.W, args.d, args.j)]
    if args.also:
        for tok in args.also.split(","):
            tok = tok.strip()
            if not tok:
                continue
            w_, d_ = tok.split(":")
            configs.append((int(w_), int(d_), None))

    all_results = []

    for W, d, j_arg in configs:
        N = W**d
        if j_arg is None:
            j_use = d // 2
        else:
            j_use = j_arg
        print(f"\n=== Config: W={W} d={d} j={j_use} N={N} ===")

        M_chi, N_chi, j_chi = make_chi_p_matrix(W, d, j_use)
        spike_chi, S_chi = analyze_spike_vectors(
            M_chi, N_chi, q_list, args.K_top, label="chi_P")
        summarize_run(spike_chi, args.k_star, label=f"chi_P W={W} d={d}")

        run_record = {
            "W": W, "d": d, "j": j_use, "N": N,
            "shape": list(M_chi.shape),
            "k_star_assumed": args.k_star,
            "chi_p_spikes": spike_chi,
            "S_chi_top60": [float(s) for s in S_chi[:60]],
        }

        if args.with_baseline:
            M_b, N_b, j_b = make_active_bernoulli_matrix(W, d, j_use, seed=42)
            spike_b, S_b = analyze_spike_vectors(
                M_b, N_b, q_list, args.K_top, label="active_bernoulli")
            summarize_run(spike_b, args.k_star, label=f"active_bernoulli W={W} d={d}")
            run_record["baseline_spikes"] = spike_b
            run_record["S_baseline_top60"] = [float(s) for s in S_b[:60]]

        all_results.append(run_record)

    out_path = Path(__file__).parent / args.out
    out_path.write_text(json.dumps(all_results, indent=2))
    print(f"\nWrote {out_path}")

    # Print PR1 / PR2 / PR3 verdict:
    print("\n=== Pre-stated falsifier verdicts ===")
    for rec in all_results:
        W, d, j_ = rec["W"], rec["d"], rec["j"]
        chi_spikes = rec["chi_p_spikes"]
        k_star = rec["k_star_assumed"]
        # PR1: each of top k_star spike vectors has max_q centered_energy > 0.5
        pr1_passes = 0
        pr1_total = min(k_star, len(chi_spikes))
        for s in chi_spikes[1:pr1_total]:  # skip k=0 (leading mean)
            if s["dom_q_centered_energy"] > 0.5:
                pr1_passes += 1
        # PR3: which q's show up in dom_q for k=1..k_star
        dom_qs = set()
        for s in chi_spikes[1:pr1_total]:
            if s["dom_q"] is not None:
                dom_qs.add(s["dom_q"])
        sum_phi_visited = sum(euler_phi(q) for q in dom_qs)
        print(f"  W={W} d={d} j={j_}: PR1 passes {pr1_passes}/{pr1_total - 1} "
              f"(target: >50% per spike); visited q={sorted(dom_qs)} "
              f"sum_phi={sum_phi_visited} k*={k_star}")
        if "baseline_spikes" in rec:
            b_spikes = rec["baseline_spikes"]
            pr2_violations = 0
            for s in b_spikes[1:pr1_total]:
                if s["dom_q_centered_energy"] > 0.2:
                    pr2_violations += 1
            print(f"    PR2: baseline violations {pr2_violations}/{pr1_total - 1} "
                  f"(target: 0)")


if __name__ == "__main__":
    main()
