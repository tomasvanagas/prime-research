"""
D22 — Higher-order Hodge Laplacian spectrum of the coprimality flag complex.

Builds K_N = the flag (clique) complex of the coprimality graph
G_N = ([2..N], {{i,j} : gcd(i,j)=1}). Computes the spectra of the Hodge
Laplacians L_0, L_1 (and at small N, L_2):

    L_k = B_k^T B_k + B_{k+1} B_{k+1}^T

with B_k the boundary matrix from k-simplices to (k-1)-simplices, oriented
ascending by sorted vertex tuple. Hodge decomposition gives
dim ker(L_k) = beta_k(K_N, R) (Eckmann 1944; Friedman 1996; Lim 2020).

Comparison: a matched-density Erdős-Rényi flag complex Y(n, p) with
n = N - 1, p = empirical edge density of G_N. For each N, S random seeds
build the ER 1-skeleton, take its flag complex, and report spectra.
Z-score is (beta_k(coprime) - mean_ER) / std_ER.

The coprimality flag complex's L_k spectrum at k >= 1 has not (to the
author's knowledge) been computed. CLOSED_PATHS lines 356, 387 only
address L_0 (graph Laplacian = Ramanujan sums = Meissel-Lehmer cost).

Cross-domain ingredient: combinatorial Hodge theory.
- Eckmann 1944, "Harmonische Funktionen und Randwertaufgaben in einem
  Komplex," Comment. Math. Helv. 17, 240.
- Friedman 1996, "Computing Betti numbers via combinatorial Laplacians,"
  Algorithmica 21, 331.
- Horak-Jost 2013, "Spectra of combinatorial Laplace operators on
  simplicial complexes," Adv. Math. 244, 303.
- Lim 2020, "Hodge Laplacians on graphs," SIAM Review 62(3), 685
  (= arXiv:1507.05379). UNUSED in CROSS_DOMAIN_TECHNIQUES §1
  (PROPOSED for D22) — first quantitative use here.
- Random flag complexes: Kahle 2009, "Topology of random clique
  complexes," Discrete Math 309, 1658; Kahle 2014, "Sharp vanishing
  thresholds for cohomology of random flag complexes," Annals 179, 1085.

Cites edges: E2.1 (MPS bond-dim φ(W)/W mechanism), E2.13 (Gowers U^k
W-trick erasure), E2.14 (Anderson chi_P), E2.16 (DPP), E2.17 (PH on
gaps), E2.19 (subword complexity), E7.12 (Cayley spectrum), E7.16
(prime-Cayley Friedman). Distinct from CLOSED 356/387 (L_0 only).
"""
from __future__ import annotations

import argparse
import json
import math
import sys
import time
from itertools import combinations
from math import gcd

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla


# ---------------------------------------------------------------------------
# Coprimality graph and flag complex
# ---------------------------------------------------------------------------

def coprime_edges(N: int) -> list[tuple[int, int]]:
    """Edges of the coprimality graph on vertex set [2..N]."""
    verts = list(range(2, N + 1))
    edges = []
    for i in range(len(verts)):
        a = verts[i]
        for j in range(i + 1, len(verts)):
            b = verts[j]
            if gcd(a, b) == 1:
                edges.append((a, b))
    return edges


def er_edges(n_vertices: int, p: float, rng: np.random.Generator) -> list[tuple[int, int]]:
    """Erdős-Rényi G(n_vertices, p) edge list on labels [0..n_vertices-1]."""
    coin = rng.random((n_vertices, n_vertices))
    edges = []
    for i in range(n_vertices):
        for j in range(i + 1, n_vertices):
            if coin[i, j] < p:
                edges.append((i, j))
    return edges


def neighbors_dict(edges: list[tuple[int, int]]) -> dict[int, set[int]]:
    nb: dict[int, set[int]] = {}
    for a, b in edges:
        nb.setdefault(a, set()).add(b)
        nb.setdefault(b, set()).add(a)
    return nb


def enumerate_triangles(edges: list[tuple[int, int]]) -> list[tuple[int, int, int]]:
    """All triangles (i<j<k) in the graph."""
    nb = neighbors_dict(edges)
    tris = []
    for a, b in edges:
        # iterate over vertices > b adjacent to both a and b
        if a not in nb or b not in nb:
            continue
        common = nb[a] & nb[b]
        for c in common:
            if c > b:
                tris.append((a, b, c))
    return tris


def enumerate_tetrahedra(edges: list[tuple[int, int]]) -> list[tuple[int, int, int, int]]:
    """All 4-cliques (i<j<k<l)."""
    nb = neighbors_dict(edges)
    tets = []
    tris = enumerate_triangles(edges)
    for a, b, c in tris:
        common = nb[a] & nb[b] & nb[c]
        for d in common:
            if d > c:
                tets.append((a, b, c, d))
    return tets


# ---------------------------------------------------------------------------
# Boundary matrices
# ---------------------------------------------------------------------------

def build_B1(verts: list[int], edges: list[tuple[int, int]]) -> sp.csr_matrix:
    """B_1 : C_1 -> C_0. shape = (|V|, |E|), sign +1 for head, -1 for tail."""
    v_idx = {v: i for i, v in enumerate(verts)}
    rows, cols, data = [], [], []
    for ek, (a, b) in enumerate(edges):
        rows.append(v_idx[a]); cols.append(ek); data.append(-1.0)  # ∂[a,b] = b - a
        rows.append(v_idx[b]); cols.append(ek); data.append(+1.0)
    return sp.csr_matrix((data, (rows, cols)), shape=(len(verts), len(edges)), dtype=np.float64)


def build_B2(edges: list[tuple[int, int]], tris: list[tuple[int, int, int]]) -> sp.csr_matrix:
    """B_2 : C_2 -> C_1. shape = (|E|, |T|).
    ∂[a,b,c] = [b,c] - [a,c] + [a,b] for a<b<c.
    """
    e_idx = {tuple(sorted(e)): i for i, e in enumerate(edges)}
    rows, cols, data = [], [], []
    for tk, (a, b, c) in enumerate(tris):
        rows.append(e_idx[(b, c)]); cols.append(tk); data.append(+1.0)
        rows.append(e_idx[(a, c)]); cols.append(tk); data.append(-1.0)
        rows.append(e_idx[(a, b)]); cols.append(tk); data.append(+1.0)
    return sp.csr_matrix((data, (rows, cols)), shape=(len(edges), len(tris)), dtype=np.float64)


def build_B3(tris: list[tuple[int, int, int]],
             tets: list[tuple[int, int, int, int]]) -> sp.csr_matrix:
    """B_3 : C_3 -> C_2. shape = (|T|, |Tet|).
    ∂[a,b,c,d] = [b,c,d] - [a,c,d] + [a,b,d] - [a,b,c].
    """
    t_idx = {tuple(sorted(t)): i for i, t in enumerate(tris)}
    rows, cols, data = [], [], []
    for tk, (a, b, c, d) in enumerate(tets):
        rows.append(t_idx[(b, c, d)]); cols.append(tk); data.append(+1.0)
        rows.append(t_idx[(a, c, d)]); cols.append(tk); data.append(-1.0)
        rows.append(t_idx[(a, b, d)]); cols.append(tk); data.append(+1.0)
        rows.append(t_idx[(a, b, c)]); cols.append(tk); data.append(-1.0)
    return sp.csr_matrix((data, (rows, cols)), shape=(len(tris), len(tets)), dtype=np.float64)


# ---------------------------------------------------------------------------
# Spectra and Betti
# ---------------------------------------------------------------------------

def hodge_spectra(B1, B2, B3, max_eig: int = 50, dense_threshold: int = 600) -> dict:
    """Compute spectra of L_0, L_1, L_2 = B_k^T B_k + B_{k+1} B_{k+1}^T.

    L_0 = B_1 B_1^T (no L_0^down)
    L_1 = B_1^T B_1 + B_2 B_2^T
    L_2 = B_2^T B_2 + B_3 B_3^T
    """
    out = {}

    # L_0
    L0 = (B1 @ B1.T).toarray() if B1.shape[0] <= dense_threshold else (B1 @ B1.T)
    if isinstance(L0, np.ndarray):
        eig0 = np.linalg.eigvalsh(L0)
    else:
        eig0 = np.sort(np.linalg.eigvalsh(L0.toarray()))
    eig0 = np.real_if_close(np.sort(eig0))
    eig0[eig0 < 1e-9] = 0.0
    out["L0_eigs"] = eig0.tolist()
    out["beta_0"] = int(np.sum(eig0 < 1e-7))

    # L_1
    nE = B1.shape[1]
    nT = B2.shape[1]
    L1_up = B2 @ B2.T  # (|E|, |E|)
    L1_dn = B1.T @ B1  # (|E|, |E|)
    L1 = L1_up + L1_dn
    if nE <= dense_threshold:
        eig1 = np.linalg.eigvalsh(L1.toarray())
        eig1 = np.real_if_close(np.sort(eig1))
        eig1[eig1 < 1e-9] = 0.0
        out["L1_eigs"] = eig1.tolist()
        out["beta_1"] = int(np.sum(eig1 < 1e-7))
        out["L1_dense"] = True
    else:
        # Sparse: report top-k and bottom-k via shift-invert
        L1s = sp.csr_matrix(L1)
        # Top eigenvalues
        if nE > 1:
            try:
                top = spla.eigsh(L1s, k=min(max_eig, nE - 1), which="LA",
                                  return_eigenvectors=False, tol=1e-5, maxiter=10000)
                top = np.sort(top.real)[::-1]
            except spla.ArpackNoConvergence as exc:
                top = np.sort(exc.eigenvalues.real)[::-1] if exc.eigenvalues.size else np.array([])
            try:
                bot = spla.eigsh(L1s, k=min(max_eig, nE - 1), which="SM",
                                  return_eigenvectors=False, tol=1e-5, maxiter=10000)
                bot = np.sort(bot.real)
            except spla.ArpackNoConvergence as exc:
                bot = np.sort(exc.eigenvalues.real) if exc.eigenvalues.size else np.array([])
        else:
            top, bot = np.array([]), np.array([])
        out["L1_top_eigs"] = top.tolist()
        out["L1_bot_eigs"] = bot.tolist()
        # beta_1 = nE - rank(B1) - rank(B2) (formal)
        # but to report a Betti we use Z-zero count among bot eigs (tolerant)
        out["beta_1"] = int(np.sum(bot < 1e-7))
        out["L1_dense"] = False

    # L_2 (only if reasonable size)
    if B3 is not None and nT > 0 and nT <= dense_threshold:
        L2_up = B3 @ B3.T
        L2_dn = B2.T @ B2
        L2 = (L2_up + L2_dn).toarray()
        eig2 = np.linalg.eigvalsh(L2)
        eig2 = np.real_if_close(np.sort(eig2))
        eig2[eig2 < 1e-9] = 0.0
        out["L2_eigs"] = eig2.tolist()
        out["beta_2"] = int(np.sum(eig2 < 1e-7))
        out["L2_dense"] = True
    elif B3 is not None and nT > 0:
        out["L2_skipped_too_large"] = True
        out["L2_size"] = nT
    else:
        out["L2_skipped_no_triangles"] = True

    return out


# ---------------------------------------------------------------------------
# Driver
# ---------------------------------------------------------------------------

def run_one(N: int, edges: list[tuple[int, int]], compute_L2: bool, max_eig: int,
            dense_threshold: int) -> dict:
    """Build B_1, B_2, optional B_3, compute spectra and Betti."""
    verts = sorted(set(v for e in edges for v in e))
    if not verts:
        verts = [0]
    nV = len(verts)
    nE = len(edges)
    tris = enumerate_triangles(edges) if edges else []
    nT = len(tris)
    if compute_L2 and tris:
        tets = enumerate_tetrahedra(edges)
    else:
        tets = []
    nTet = len(tets)

    B1 = build_B1(verts, edges) if edges else sp.csr_matrix((nV, 0), dtype=np.float64)
    B2 = (build_B2(edges, tris) if tris else sp.csr_matrix((nE, 0), dtype=np.float64))
    if compute_L2 and tets:
        B3 = build_B3(tris, tets)
    else:
        B3 = sp.csr_matrix((nT, 0), dtype=np.float64)

    spec = hodge_spectra(B1, B2, B3, max_eig=max_eig, dense_threshold=dense_threshold)
    spec.update({
        "N_label": N,
        "n_vertices": nV,
        "n_edges": nE,
        "n_triangles": nT,
        "n_tetrahedra": nTet,
        "edge_density": (2 * nE) / (nV * (nV - 1)) if nV > 1 else 0.0,
    })
    return spec


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--N-list", type=int, nargs="+", default=[32, 64, 128, 256])
    parser.add_argument("--seeds", type=int, default=20)
    parser.add_argument("--max-eig", type=int, default=50,
                        help="Number of leading/trailing eigenvalues if sparse.")
    parser.add_argument("--dense-threshold", type=int, default=900,
                        help="If matrix dim ≤ this, do dense diag.")
    parser.add_argument("--out", type=str, required=True)
    args = parser.parse_args()

    rng = np.random.default_rng(20260427)
    results = {"args": vars(args), "by_N": {}}

    for N in args.N_list:
        t0 = time.time()
        # Coprimality flag complex
        cp_edges = coprime_edges(N)
        nV_cp = N - 1
        # Decide whether to compute L_2 (need |triangles| ≤ dense_threshold for full)
        # Always enumerate triangles; only build B_3 if reasonable
        cp_result = run_one(N, cp_edges, compute_L2=True,
                            max_eig=args.max_eig, dense_threshold=args.dense_threshold)
        cp_result["wall_time_sec"] = time.time() - t0
        cp_result["edge_density"] = (2 * len(cp_edges)) / (nV_cp * (nV_cp - 1)) if nV_cp > 1 else 0.0

        p_match = cp_result["edge_density"]

        # ER controls
        ctrl_results = []
        for s in range(args.seeds):
            t1 = time.time()
            er_e = er_edges(nV_cp, p_match, rng)
            ctrl = run_one(-1, er_e, compute_L2=True,
                           max_eig=args.max_eig, dense_threshold=args.dense_threshold)
            ctrl["wall_time_sec"] = time.time() - t1
            ctrl["seed"] = s
            ctrl_results.append(ctrl)
            # Write progress for long-running runs
            if (s + 1) % 5 == 0:
                print(f"  N={N} seed {s+1}/{args.seeds} done", flush=True)

        # Aggregate
        record = {
            "coprime": cp_result,
            "controls_p": p_match,
            "controls": ctrl_results,
        }
        # Z-scores for beta_0, beta_1, beta_2
        for k in (0, 1, 2):
            key = f"beta_{k}"
            cp_val = cp_result.get(key, None)
            ctrl_vals = [c.get(key, None) for c in ctrl_results if c.get(key) is not None]
            if cp_val is None or len(ctrl_vals) < 2:
                continue
            mu = float(np.mean(ctrl_vals))
            sd = float(np.std(ctrl_vals, ddof=1))
            record[f"z_{key}"] = {
                "cp": cp_val,
                "mean_er": mu,
                "std_er": sd,
                "z": (cp_val - mu) / sd if sd > 0 else float("nan"),
                "ctrl_n": len(ctrl_vals),
            }
        # Spectral comparison (eigenvalue moments)
        # KS distance between coprime L_1 spectrum and pooled-ER L_1 spectrum
        cp_eigs1 = cp_result.get("L1_eigs", None)
        if cp_eigs1 is not None:
            er_eigs1 = []
            for c in ctrl_results:
                e1 = c.get("L1_eigs", None)
                if e1 is not None:
                    er_eigs1.extend(e1)
            if er_eigs1:
                from scipy.stats import ks_2samp
                ks = ks_2samp(cp_eigs1, er_eigs1)
                record["L1_ks_stat"] = float(ks.statistic)
                record["L1_ks_pval"] = float(ks.pvalue)
                # Moment comparison
                cp_mean = float(np.mean(cp_eigs1))
                cp_std = float(np.std(cp_eigs1))
                er_mean = float(np.mean([np.mean(c["L1_eigs"]) for c in ctrl_results
                                          if "L1_eigs" in c]))
                er_std_of_means = float(np.std([np.mean(c["L1_eigs"]) for c in ctrl_results
                                                  if "L1_eigs" in c], ddof=1))
                record["L1_mean_cp"] = cp_mean
                record["L1_mean_er"] = er_mean
                record["L1_std_of_er_means"] = er_std_of_means
                record["L1_z_mean"] = (
                    (cp_mean - er_mean) / er_std_of_means if er_std_of_means > 0 else float("nan")
                )

        results["by_N"][str(N)] = record
        print(f"N={N} done in {time.time()-t0:.1f}s. "
              f"|V|={cp_result['n_vertices']} |E|={cp_result['n_edges']} "
              f"|T|={cp_result['n_triangles']} "
              f"beta_0={cp_result['beta_0']} beta_1={cp_result['beta_1']} "
              f"beta_2={cp_result.get('beta_2')} dens={cp_result['edge_density']:.4f}",
              flush=True)
        # incremental write
        with open(args.out, "w") as f:
            json.dump(results, f, indent=1)

    print(f"Wrote {args.out}", flush=True)


if __name__ == "__main__":
    main()
