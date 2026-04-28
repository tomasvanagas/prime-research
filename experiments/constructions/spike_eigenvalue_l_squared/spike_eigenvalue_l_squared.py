#!/usr/bin/env python3
"""S82 commit-thread step 1 — eigenvalue ~ |L(1, chi)|^2 prediction test.

S82 conjectures: spike eigenvectors of the chi_P MPS-Gram matrix M^T M
are residue-class character vectors with eigenvalue ~ |L(1, chi)|^2.

This script tests the eigenvalue claim numerically. It does NOT prove the
theorem; it asks whether the empirical spike sigma_k^2 values from S82's
JSON dumps are consistent with a formula

    sigma_chi^2 ~ K(N) * |L(1, chi)|^2

across (a) different conductors at fixed d, and (b) the same conductor at
different d. K(N) is allowed to depend on N but not chi.

For each conductor q = 2*p (p odd prime, p ∤ W=2), the spike block has
phi(p) singular values per S82's empirical finding. We total their
sigma^2's and compare with sum_{chi primitive mod p} |L(1, chi)|^2 plus
the "indicator-of-(p ∤ n)" mode contribution.

Key reference for L(1, chi) formula:
  Davenport, "Multiplicative Number Theory" (3rd ed, GTM 74), Ch. 6.
    L(1, chi) for primitive chi mod q given by
      L(1, chi) = -(1/q) * sum_{a=1}^{q-1} chi(a) * psi(a/q)
    where psi is the digamma function (psi(x) = Gamma'(x)/Gamma(x)).
"""

import argparse
import cmath
import json
import math
from pathlib import Path

import numpy as np
from scipy.special import digamma


# ----- Dirichlet character enumeration -----

def primitive_root(p):
    """Smallest primitive root mod p (p prime)."""
    for g in range(2, p):
        seen = set()
        x = 1
        for _ in range(p - 1):
            x = (x * g) % p
            seen.add(x)
        if len(seen) == p - 1:
            return g
    raise ValueError(f"no primitive root for {p}")


def characters_mod_p(p):
    """All Dirichlet characters mod p (p prime), as a dict a -> chi(a).

    Returns list of (label, chi) where chi is dict {0: 0, 1: chi(1), ...}.
    Includes the principal character (chi_0).
    """
    g = primitive_root(p)
    log_g = {1: 0}
    x = 1
    for k in range(1, p - 1):
        x = (x * g) % p
        log_g[x] = k
    log_g[0] = None  # 0 has no log

    chars = []
    for j in range(p - 1):
        chi = {0: 0}
        for a in range(1, p):
            k = log_g[a]
            chi[a] = cmath.exp(2j * math.pi * j * k / (p - 1))
        label = f"chi_{j}_mod_{p}"
        chars.append((label, chi, j))
    return chars  # j=0 is principal


def L1_of_char_mod_p(chi, p):
    """Compute L(1, chi) for chi non-principal mod p (p prime).

    Uses the closed-form identity for primitive chi mod q:
      L(1, chi) = -(1/q) * sum_{a=1}^{q-1} chi(a) * psi(a/q)
    where psi(x) = digamma(x). For chi mod prime p, chi is primitive iff
    non-principal.

    Returns the complex value L(1, chi).
    """
    val = 0.0 + 0.0j
    for a in range(1, p):
        val += chi[a] * digamma(a / p)
    return -val / p


def is_principal(chi, p):
    """True if chi is the principal character mod p."""
    for a in range(1, p):
        if abs(chi[a] - 1.0) > 1e-10:
            return False
    return True


# ----- Predicted eigenvalue contributions -----

def L_squared_contributions(p):
    """For prime p, return (sum of |L(1, chi)|^2 over non-principal chi mod p,
    list of individual |L(1, chi)|^2 values).

    These are the per-character contributions to spike sigma^2 in the
    "conductor 2p" block (under the S82 hypothesis with eigenvalue ~|L(1,chi)|^2).
    """
    chars = characters_mod_p(p)
    contribs = []
    for label, chi, j in chars:
        if j == 0:
            continue  # skip principal
        L1 = L1_of_char_mod_p(chi, p)
        contribs.append((label, abs(L1) ** 2, L1))
    total = sum(c[1] for c in contribs)
    return total, contribs


# ----- Empirical extraction from S82 JSON -----

def empirical_block_totals(json_path, conductor_groups):
    """Load S82 JSON and total sigma^2 by conductor block.

    conductor_groups: dict mapping group name to list of conductors that
    fall in that group, e.g. {"mod_3": [6], "mod_5": [10], ...}.

    Returns dict: group_name -> {"total_sigma_sq": ..., "n_spikes": ...,
                                 "sigmas": [...], "ks": [...]}
    """
    with open(json_path) as f:
        data = json.load(f)
    rec = data[0]
    spikes = rec["chi_p_spikes"]

    out = {}
    for group, qs in conductor_groups.items():
        total = 0.0
        sigmas = []
        ks = []
        for s in spikes:
            cond = s.get("min_q_conductor")
            if cond in qs:
                total += s["sigma"] ** 2
                sigmas.append(s["sigma"])
                ks.append(s["k"])
        out[group] = {
            "total_sigma_sq": total,
            "n_spikes": len(sigmas),
            "sigmas": sigmas,
            "ks": ks,
        }
    out["_meta"] = {"N": rec["N"], "d": rec["d"], "W": rec["W"], "shape": rec["shape"]}
    return out


# ----- Main -----

def run_q_sweep(json_paths):
    """Run the test across multiple JSON dumps."""
    print("=== Predicted |L(1, chi)|^2 totals per conductor-2p block ===\n")

    L_table = {}
    for p in [3, 5, 7, 11, 13]:
        total_L2, contribs = L_squared_contributions(p)
        L_table[p] = {"total": total_L2, "contribs": contribs}
        print(f"p={p}: phi(p)={p-1} non-principal chars, sum |L(1,chi)|^2 = {total_L2:.6f}")
        for label, L2, L1 in contribs:
            print(f"    {label}: L(1,chi) = {L1:.5f},  |L|^2 = {L2:.5f}")
        print()

    print("=== Empirical block totals from S82 JSON dumps ===\n")

    # Conductor groupings: 2p-conductor + cross conductors.
    cgroups = {
        "mod_3 (q=6)": [6],
        "mod_5 (q=10)": [10],
        "mod_7 (q=14)": [14],
        "cross_3_5 (q=30)": [30],
        "mod_11 (q=22 or 88)": [22, 88, 44],
        "cross_3_7 (q=42)": [42],
        "cross_5_7 (q=70)": [70],
    }

    rows = []  # (d, group, total_sigma_sq, n_spikes, predicted_L2_sum)
    for jp in json_paths:
        emp = empirical_block_totals(jp, cgroups)
        d = emp["_meta"]["d"]
        N = emp["_meta"]["N"]
        print(f"\n--- d={d}  N={N}  shape={emp['_meta']['shape']} ---")
        for group, info in emp.items():
            if group.startswith("_"):
                continue
            n = info["n_spikes"]
            tot = info["total_sigma_sq"]
            sigmas = ", ".join(f"{s:.2f}" for s in info["sigmas"][:6])
            if n > 6:
                sigmas += ", ..."
            print(f"  {group}: n={n}  sum sigma^2 = {tot:.2f}  sigmas=[{sigmas}]")
            rows.append((d, N, group, tot, n, info["sigmas"]))

    # Compare ratios across blocks at fixed d.
    print("\n\n=== Block-ratio test (does sigma^2_q / sigma^2_{q'} match L^2 ratios?) ===\n")

    # For each d, compute ratio (block sigma^2 sum) / (predicted L^2 sum) per p.
    by_d = {}
    for (d, N, group, tot, n, _) in rows:
        by_d.setdefault((d, N), {})[group] = (tot, n)

    p_groups = {3: "mod_3 (q=6)", 5: "mod_5 (q=10)",
                7: "mod_7 (q=14)", 11: "mod_11 (q=22 or 88)"}

    print(f"{'d':>4} {'N':>10} {'p':>4} {'n_spikes':>10} {'sum_sigma^2':>14} "
          f"{'sum |L|^2':>12} {'ratio = sigma^2/|L|^2':>22}")
    for (d, N), groups in sorted(by_d.items()):
        for p, gname in p_groups.items():
            if gname not in groups:
                continue
            tot, n = groups[gname]
            L2_sum = L_table[p]["total"]
            ratio = tot / L2_sum if L2_sum > 0 else float("nan")
            print(f"{d:>4} {N:>10} {p:>4} {n:>10} {tot:>14.2f} "
                  f"{L2_sum:>12.5f} {ratio:>22.2f}")

    # If the prediction sigma^2 ~ K(N) * |L|^2 holds, the "ratio" column
    # should be approximately constant ACROSS p at FIXED d.

    print("\n\n=== Constancy-of-ratio diagnostic ===\n")
    print("If sigma^2 = K(N) * |L|^2 holds, ratio (sigma^2 / |L|^2) is constant across p at fixed d.")
    print("Compute the std/mean of the ratio across p={3,5,7,11} for each d:\n")

    print(f"{'d':>4} {'N':>10}  ratios (p=3,5,7,11)                       mean    std    cv")
    for (d, N), groups in sorted(by_d.items()):
        ratios = []
        for p, gname in p_groups.items():
            if gname not in groups:
                continue
            tot, n = groups[gname]
            L2_sum = L_table[p]["total"]
            if L2_sum > 0 and n >= max(p - 2, 1):  # only blocks with full saturation
                ratios.append(tot / L2_sum)
        if not ratios:
            continue
        arr = np.array(ratios)
        mean = arr.mean()
        std = arr.std()
        cv = std / mean if mean > 0 else float("nan")
        rstr = " ".join(f"{r:>7.1f}" for r in ratios)
        print(f"{d:>4} {N:>10}  [{rstr}]  {mean:>7.1f} {std:>7.1f} {cv:>5.3f}")

    print("\nLow CV (<0.2) => |L(1,chi)|^2 scaling holds as predicted by S82.")
    print("High CV (>0.5) => the eigenvalue is NOT simply proportional to |L|^2; the")
    print("                  S82 prediction is wrong or needs a per-q correction.")

    # ----- Alternative scaling test: sigma^2 ~ N / (p * log^2 N) -----
    print("\n\n=== Alternative scaling test: sum_sigma^2 ~ K * N / (p * log^2 N) ===\n")
    print(f"{'d':>4} {'N':>10} {'p':>4} {'sum_sigma^2':>14} "
          f"{'K = sigma^2 * p * log^2N / N':>32}")
    K_table = []
    for (d, N), groups in sorted(by_d.items()):
        log2N = math.log(N) ** 2
        for p, gname in p_groups.items():
            if p == 11:
                continue  # mod-11 has cross contamination + incomplete saturation
            if gname not in groups:
                continue
            tot, n = groups[gname]
            if n < p - 2:
                continue  # need saturated block
            K = tot * p * log2N / N
            K_table.append((d, N, p, K))
            print(f"{d:>4} {N:>10} {p:>4} {tot:>14.2f} {K:>32.4f}")

    Karr = np.array([k for (_, _, _, k) in K_table])
    print(f"\n  Mean K = {Karr.mean():.4f},  std K = {Karr.std():.4f},  CV = {Karr.std()/Karr.mean():.3f}")
    print(f"  Equivalent: sum_sigma^2 ~ {Karr.mean():.3f} * N / (p * log^2 N)")
    print(f"            ~ {Karr.mean():.3f} * pi(N)^2 / (N * p)  by PNT")
    print(f"            ~ {Karr.mean():.3f} * pi(N) / (p * log N)")
    print("\nA constant K (low CV across (d, p)) means the eigenvalue scales as a")
    print("PNT-in-arithmetic-progressions variance, not as |L(1, chi)|^2.")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--json_dir", type=str,
                        default="../spike_eigenvectors_chi_p")
    args = parser.parse_args()

    here = Path(__file__).parent
    json_dir = (here / args.json_dir).resolve()
    json_paths = sorted([p for p in json_dir.glob("spike_d*.json")
                         if "chi_only" not in p.name and "w6" not in p.name])
    print(f"JSON files to test: {[p.name for p in json_paths]}")

    run_q_sweep(json_paths)


if __name__ == "__main__":
    main()
