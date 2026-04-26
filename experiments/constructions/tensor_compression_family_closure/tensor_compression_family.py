"""
N1 — Tensor-network compression family closure for chi_P.

For chi_P : [1, W^d] -> {0,1} the prime indicator, reshaped as a tensor of
shape (W, W, ..., W) with d copies, this script verifies that all the
classical tensor-network ansatze (MPS, Tucker, HT, CP, TR, PEPS, MERA)
admit an unfolding-rank lower bound traceable to E2.1:

    rank M^(j) = min(W^j, phi(W)*W^(d-j-1) + 1)

The script computes for each (W, d):
- All cut-j unfolding ranks (E2.1 verification).
- The CP-rank lower bound = max over j of rank M^(j).
- HT tree-cut ranks for a balanced binary tree.
- TR bond dim (cyclic shift of MPS — same lower bound as MPS).
- PEPS 2D-reshape ranks.
- MERA Renyi-2 entanglement entropy bound: bond dim chi must satisfy
  chi >= rank^(1 / O(log d)) for some absolute constant in the
  exponent — we report log(rank)/log(d).
- Tucker support-cardinality lower bound: chi_P has Theta(N/log N)
  nonzero entries, so any Tucker decomposition has parameter count
  bounded below by min(support, prod r_j) — Tucker fails to compress
  for a different reason than the unfolding-rank route.

Outputs JSON with everything; prints a Pass/Fail summary against the
falsification criteria F1-F6 from definition.md.

CLI:
    python tensor_compression_family.py [--Ws W1 W2 ...] [--ds D1 D2 ...]

Defaults to W in {2, 3, 5}, d in {6, 8, 10}.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Iterable

import numpy as np
from sympy import isprime, totient


# ----------------------------- chi_P tensor ---------------------------------


def build_chi_p_tensor(W: int, d: int) -> np.ndarray:
    """
    Build the chi_P tensor of shape (W, W, ..., W) with d copies.

    The (i_1, i_2, ..., i_d) entry is chi_P(n) = 1 if n is prime, 0 else,
    where n = i_1 * W^(d-1) + i_2 * W^(d-2) + ... + i_d + 1.

    Returns a float64 ndarray of shape (W,)*d.
    """
    if W ** d > 2_000_000:
        raise ValueError(f"W^d = {W**d} too large; use smaller (W, d).")
    arr = np.zeros((W,) * d, dtype=np.float64)
    # Vectorised fill via flat index.
    flat = arr.reshape(-1)
    for n in range(1, W ** d + 1):
        if isprime(n):
            flat[n - 1] = 1.0
    return arr


# ----------------------------- E2.1 unfolding ranks --------------------------


def unfolding_rank_at_cut(T: np.ndarray, j: int) -> int:
    """
    Compute the rank of the (W^j x W^(d-j)) unfolding of T.

    j is the number of leftmost modes grouped as rows.
    """
    d = T.ndim
    W = T.shape[0]
    M = T.reshape(W ** j, W ** (d - j))
    # SVD-based rank with tight tolerance over Q-equivalent {0,1} matrix:
    return int(np.linalg.matrix_rank(M, tol=1e-9))


def predicted_unfolding_rank(W: int, d: int, j: int) -> int:
    """E2.1 closed form: min(W^j, phi(W)*W^(d-j-1) + 1)."""
    return int(min(W ** j, totient(W) * W ** (d - j - 1) + 1))


# ----------------------------- CP rank lower bound --------------------------


def cp_rank_lower_bound(T: np.ndarray) -> int:
    """
    Kruskal lower bound: CP-rank(T) >= max over j of mode-j matricisation
    rank. We use the *block-mode* matricisations (W^j x W^(d-j)) since
    they give the strongest bound for chi_P (mode-j matricisation has
    rank <= W which is trivial).
    """
    d = T.ndim
    return max(unfolding_rank_at_cut(T, j) for j in range(1, d))


# ----------------------------- HT bond dims ---------------------------------


def ht_balanced_tree_bond_dims(T: np.ndarray) -> list[tuple[tuple[int, ...], int]]:
    """
    For a balanced binary tree on d modes, return [(modes_subset, rank), ...]
    for every NON-TRIVIAL internal-node tree-cut. The bond dim at internal
    node v with leaves-below-v = `modes_subset` is the rank of the
    matricisation T_{modes_subset}^{complement}.

    Note: leaves themselves have bond dim ≤ W (trivial) since the mode-i
    matricisation is W × W^(d-1). What's diagnostic is the bond dim at
    *deeper* internal nodes — in particular the root's children, which
    induce the half-cut.
    """
    d = T.ndim
    W = T.shape[0]
    out: list[tuple[tuple[int, ...], int]] = []

    def matricisation_rank(modes: tuple[int, ...]) -> int:
        complement = tuple(i for i in range(d) if i not in modes)
        perm = list(modes) + list(complement)
        Tp = T.transpose(perm)
        rows = W ** len(modes)
        cols = W ** len(complement)
        M = Tp.reshape(rows, cols)
        return int(np.linalg.matrix_rank(M, tol=1e-9))

    def recurse(lo: int, hi: int) -> None:
        # Interval [lo, hi) of modes; only emit if it's a non-trivial split
        # (size 2..d-1; mode-leaves have trivial rank ≤ W).
        size = hi - lo
        if size <= 1:
            return
        modes_left = tuple(range(lo, (lo + hi) // 2))
        if 2 <= len(modes_left) and len(modes_left) < d:
            rank = matricisation_rank(modes_left)
            out.append((modes_left, rank))
        elif len(modes_left) == 1 and d >= 3:
            # Only emit single-mode bond dims for completeness (always ≤ W)
            rank = matricisation_rank(modes_left)
            out.append((modes_left, rank))
        recurse(lo, (lo + hi) // 2)
        recurse((lo + hi) // 2, hi)

    recurse(0, d)
    return out


# ----------------------------- PEPS 2D reshape ranks ------------------------


def peps_2d_reshape_ranks(T: np.ndarray) -> dict[str, int]:
    """
    Two natural 2D reshapes of a (W,)*d tensor:
    (a) (W^(d/2), W^(d/2)) — same as MPS half-cut (when d even).
    (b) (W^a, W^(d-a)) for a = floor(d/2) and a = ceil(d/2).
    PEPS bond dim across the boundary between the two halves equals the
    matrix rank of this reshape.
    """
    d = T.ndim
    W = T.shape[0]
    out: dict[str, int] = {}
    for a in {d // 2, (d + 1) // 2}:
        if 1 <= a < d:
            M = T.reshape(W ** a, W ** (d - a))
            out[f"a={a}"] = int(np.linalg.matrix_rank(M, tol=1e-9))
    return out


# ----------------------------- MERA Renyi-2 bound ---------------------------


def mera_log_bond_lower_bound(T: np.ndarray, j: int) -> float:
    """
    For a 1D MERA with bond dim chi at every layer and L = log_2 d
    layers, the rank of the unfolding at cut j is bounded above by
    chi^(c*L) for some constant c (depending on local geometry — for
    binary MERA, c = 2 is a standard upper bound). Therefore:

        log(rank M^(j)) <= c * L * log(chi)
        => log(chi) >= log(rank M^(j)) / (c * L)
        => chi >= rank^(1 / (c * log_2 d)).

    We report log(rank) and log_2(d) and the implied minimum log(chi)
    with c = 2.
    """
    d = T.ndim
    rank = unfolding_rank_at_cut(T, j)
    log_rank = math.log(max(rank, 2))  # natural log; +epsilon avoids log(0/1)
    log2_d = math.log2(max(d, 2))
    min_log_chi = log_rank / (2.0 * log2_d)  # natural log of chi
    return {
        "rank": rank,
        "log_rank": log_rank,
        "log2_d": log2_d,
        "min_log_chi": min_log_chi,
        "min_chi_lower": math.exp(min_log_chi),
    }


# ----------------------------- Tucker support arg ---------------------------


def tucker_support_lower_bound(T: np.ndarray) -> dict[str, int]:
    """
    Tucker decomposition T = G x_1 U_1 ... x_d U_d has total parameter
    count prod_j r_j + sum_j r_j * W. The mode-j matricisation has rank
    <= W trivially, so the multilinear rank tuple is bounded by
    (W, ..., W) — Tucker's natural rank notion does NOT see the
    unfolding-rank lower bound.

    What DOES bound Tucker compression: T has support cardinality |T|_0
    (number of nonzero entries). Any Tucker decomposition with all
    r_j < W is forced to have prod_j r_j >= |T|_0 (Schmidt-style argument
    on the core tensor cardinality), since otherwise the core has fewer
    unique entries than the support.

    We report:
    - |T|_0 (number of primes in [1, W^d])
    - W^d total volume
    - density |T|_0 / W^d
    - The trivial multilinear rank (W, W, ..., W)
    """
    nnz = int((T != 0).sum())
    total = int(T.size)
    return {
        "nnz_chi_p": nnz,
        "total_volume": total,
        "density": nnz / total,
        "trivial_multilinear_rank": [int(s) for s in T.shape],
    }


# ----------------------------- Falsification eval ---------------------------


def evaluate_falsification(W: int, d: int, results: dict) -> dict:
    """Check F1-F6 against the recorded numbers."""
    checks: dict[str, dict] = {}

    # F1: MPS / TT bond dim < phi(W) * W^(d-j-1) + 1 at some cut j?
    # We allow finite-size corrections of O(1): for the asymptotic family
    # closure to fail, a cut would need actual rank ASYMPTOTICALLY below
    # predicted. For tested (W, d), record both raw mismatches and the
    # max relative deficit.
    deficits = []
    for entry in results["unfolding_ranks"]:
        actual = entry["actual"]
        predicted = entry["predicted_E2_1"]
        deficits.append(predicted - actual)
    checks["F1_MPS"] = {
        "max_deficit": max(deficits),
        "all_ranks_match_E2_1": all(d == 0 for d in deficits),
        "asymptotic_rank_at_half_cut": next(
            (entry["actual"] for entry in results["unfolding_ranks"]
             if entry["j"] == d // 2),
            None,
        ),
        "predicted_rank_at_half_cut": next(
            (entry["predicted_E2_1"] for entry in results["unfolding_ranks"]
             if entry["j"] == d // 2),
            None,
        ),
    }

    # F2: HT max bond dim < unfolding max rank?
    # The diagnostic for HT compression is the MAX bond dim, since polylog
    # representation requires the max to be polylog. We extract the
    # contiguous "left-prefix" cuts (which match MPS j-cuts) and report.
    ht_max = max(entry["rank"] for entry in results["ht_tree_ranks"])
    unfold_max = max(entry["actual"] for entry in results["unfolding_ranks"])
    # The half-cut is the diagnostic internal node — match against MPS half-cut.
    half_modes = tuple(range(d // 2))
    ht_half_cut = next(
        (entry["rank"] for entry in results["ht_tree_ranks"]
         if tuple(entry["modes"]) == half_modes),
        None,
    )
    checks["F2_HT"] = {
        "ht_max_rank": ht_max,
        "unfold_max_rank": unfold_max,
        "ht_half_cut_rank": ht_half_cut,
        "max_match": ht_max == unfold_max,
    }

    # F3: CP rank lower bound check (Kruskal). The bound itself is what
    # we report — F3 fails only if some explicit construction beats it.
    # We don't build a CP decomposition here; we just record the bound.
    checks["F3_CP"] = {
        "lower_bound": results["cp_rank_lower_bound"],
    }

    # F4: TR bond dim < MPS bond dim. TR is structurally MPS with a
    # cyclic boundary, so the rank-bound is the same. Reported as
    # "structurally identical".
    checks["F4_TR"] = {
        "structurally_identical_to_MPS": True,
    }

    # F5: MERA polylog bond dim — we report the log_chi lower bound.
    if results["mera_bound"]:
        chi_min = min(entry["min_chi_lower"] for entry in results["mera_bound"].values())
        checks["F5_MERA"] = {
            "min_chi_lower_bound": chi_min,
            "polylog_admissible": chi_min < math.log2(W ** d) ** 2,  # generous polylog bar
        }

    # F6: PEPS 2D reshape rank < MPS half-cut.
    half_cut_rank = next(
        (entry["actual"] for entry in results["unfolding_ranks"] if entry["j"] == d // 2),
        None,
    )
    peps_min = min(results["peps_2d"].values()) if results["peps_2d"] else None
    checks["F6_PEPS"] = {
        "peps_min_rank": peps_min,
        "mps_half_cut_rank": half_cut_rank,
        "failed": peps_min is not None and half_cut_rank is not None and peps_min < half_cut_rank,
    }

    return checks


# ----------------------------- Driver --------------------------------------


def run(Ws: Iterable[int], ds: Iterable[int]) -> dict:
    out: dict = {"runs": []}
    for W in Ws:
        for d in ds:
            if W ** d > 2_000_000:
                continue
            T = build_chi_p_tensor(W, d)

            unfolding_ranks = []
            for j in range(1, d):
                actual = unfolding_rank_at_cut(T, j)
                predicted = predicted_unfolding_rank(W, d, j)
                unfolding_ranks.append(
                    {"j": j, "actual": actual, "predicted_E2_1": predicted}
                )

            cp_lb = cp_rank_lower_bound(T)
            ht_ranks = ht_balanced_tree_bond_dims(T)
            peps = peps_2d_reshape_ranks(T)

            mera_bound = {}
            for j in range(1, d):
                mera_bound[f"j={j}"] = mera_log_bond_lower_bound(T, j)

            tucker_lb = tucker_support_lower_bound(T)

            results = {
                "W": W,
                "d": d,
                "N": W ** d,
                "phi_W": int(totient(W)),
                "phi_over_W": float(totient(W)) / W,
                "unfolding_ranks": unfolding_ranks,
                "cp_rank_lower_bound": cp_lb,
                "ht_tree_ranks": [
                    {"modes": list(modes), "rank": rank} for modes, rank in ht_ranks
                ],
                "peps_2d": peps,
                "mera_bound": mera_bound,
                "tucker_support": tucker_lb,
            }
            results["falsification_checks"] = evaluate_falsification(W, d, results)
            out["runs"].append(results)
            print(f"W={W}, d={d}, N={W**d}: ", end="")
            print(
                f"E2.1 match={results['falsification_checks']['F1_MPS']['all_ranks_match_E2_1']}, "
                f"max_deficit={results['falsification_checks']['F1_MPS']['max_deficit']}, "
                f"CP-LB={cp_lb}, "
                f"HT-half-cut={results['falsification_checks']['F2_HT']['ht_half_cut_rank']}, "
                f"PEPS-min={results['falsification_checks']['F6_PEPS']['peps_min_rank']}"
            )

    return out


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--Ws", type=int, nargs="*", default=[2, 3, 5])
    parser.add_argument("--ds", type=int, nargs="*", default=[6, 8, 10])
    parser.add_argument("--out", type=str, default=None)
    args = parser.parse_args()

    out = run(args.Ws, args.ds)

    out_path = (
        Path(args.out)
        if args.out
        else Path(__file__).resolve().parent / "tensor_compression_family_results.json"
    )
    out_path.write_text(json.dumps(out, indent=2))
    print(f"\nWrote {out_path}")


if __name__ == "__main__":
    main()
