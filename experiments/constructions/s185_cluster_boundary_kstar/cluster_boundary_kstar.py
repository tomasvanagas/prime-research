"""
S185 (verify-15 of S169): cluster-boundary k_* probe.

Adversarial test: the 21% spike-block fraction in S169/S183 is computed
at k_* extrapolated linearly from S82's three-point fit. Is this k_*
choice structurally privileged, or only one of many reasonable k_*
selections each giving a different fraction?

We compute the spike block / pi(N) ratio at MULTIPLE k_* selection rules
applied to the same SVD data:

    R0  S82-extrapolated (the canonical S169/S183 choice)
    R1  Marchenko-Pastur upper edge: k_* = #{sigma_k > 2*sqrt(M*p_N)}
    R2  Spectral elbow: k_* = argmax_k of (sigma_k - sigma_{k+1})/sigma_k
    R3  Largest plateau-end gap: k_* at the deepest gap below k <= 2*M^{1/2}
    R4  Sector-completion (cumulative phi(q) up to small primes)

If the 21% asymptote is a real structural fact, its value should be
similar across these rules. If it depends on the canonical k_* choice,
the spread reveals model-dependence.

Saved sigma data:
    d in {14, 18, 20}: spike_eigenvectors_chi_p/spike_d{d}_results.json
                       (top 25, 40, 50 sigmas respectively)
    d = 24:            s183_d24_svd_verify/d24_svd_verify_results.json
                       (top 100 sigmas)

For d in {14, 18, 20} we also recompute the FULL sigma spectrum from
chi_P to validate our spike-block sums against S169's reported numbers.

What would falsify our finding: the cluster-boundary fractions land in
[0.20, 0.22] across all four rules, confirming that 0.21 is robust to
k_* choice. What would confirm our finding: a spread > 0.04 across
natural rules, indicating that the "0.21" choice is one in a family.
"""

import json
import math
import os
import sys

import numpy as np
from sympy import primerange

DATA_ROOT = "/apps/aplikacijos/prime-research/experiments/constructions"
S82_DIR = os.path.join(DATA_ROOT, "spike_eigenvectors_chi_p")
S183_FILE = os.path.join(DATA_ROOT, "s183_d24_svd_verify", "d24_svd_verify_results.json")


def load_s82_sigmas():
    out = {}
    for d in (14, 18, 20):
        with open(os.path.join(S82_DIR, f"spike_d{d}_results.json")) as f:
            payload = json.load(f)[0]
        spikes = payload["chi_p_spikes"]
        sigmas = [s["sigma"] for s in spikes]
        out[d] = {
            "N": payload["N"],
            "k_star_assumed": payload["k_star_assumed"],
            "sigmas": sigmas,
            "shape": payload["shape"],
            # only stored sigmas; full SVD will be recomputed below
        }
    return out


def load_s183_sigmas():
    with open(S183_FILE) as f:
        payload = json.load(f)
    return {
        24: {
            "N": payload["N"],
            "pi_N": payload["pi_N"],
            "shape": payload["shape"],
            "sigmas": payload["sigmas_top100"],
            "frobenius_sq": payload["frobenius_sq"],
            "sigma_0_sq": payload["sigma_0_sq"],
        }
    }


def chi_P(N):
    arr = np.zeros(N, dtype=np.float64)
    for p in primerange(2, N + 1):
        arr[p - 1] = 1.0
    return arr


def full_svd_sigmas(d):
    """Recompute full sigma spectrum from chi_P at scale d."""
    N = 1 << d
    cp = chi_P(N)
    M = 1 << (d // 2)
    # MPS unfolding: 2^d = 2^{d/2} * 2^{d/2} reshape (S82 / S169 convention)
    A = cp.reshape(M, M)
    sigmas = np.linalg.svd(A, compute_uv=False)
    return N, M, sigmas


def cum_block(sigmas, k_star, exclude_k0=True):
    """Sum sigma^2 from k=1..k_star inclusive (excluding rank-1 mean k=0)."""
    if exclude_k0:
        return float(np.sum(np.asarray(sigmas[1 : k_star + 1]) ** 2))
    return float(np.sum(np.asarray(sigmas[: k_star + 1]) ** 2))


def k_star_mp_edge(sigmas, M, pi_N):
    """Marchenko-Pastur upper edge for chi_P after rank-1 mean removal."""
    p_N = pi_N / (M * M)
    # variance per entry approx p_N (since chi is 0/1, after centering ~ p(1-p))
    edge = 2.0 * math.sqrt(M * p_N * (1.0 - p_N))
    # count sigmas strictly above edge, excluding sigma_0 (rank-1 mean)
    s = np.asarray(sigmas[1:])
    return int(np.sum(s > edge)), edge


def k_star_spectral_elbow(sigmas, k_min=2, k_max=None):
    """k where (sigma_k - sigma_{k+1})/sigma_k is maximised in [k_min, k_max].

    We exclude k=0,1,2 (always huge gaps for V_3 cluster) and look in the
    bulk-transition region where the spike block ends. To target the
    *bulk-transition* elbow rather than the V_3/V_5 plateau gaps that
    dominate the early spectrum, we use a k_min that skips past the
    deterministic prime-character clusters.
    """
    s = np.asarray(sigmas, dtype=np.float64)
    if k_max is None:
        k_max = len(s) - 2
    rel = (s[k_min:k_max] - s[k_min + 1 : k_max + 1]) / s[k_min:k_max]
    k_idx = int(np.argmax(rel)) + k_min
    return k_idx, float(rel[k_idx - k_min])


def k_star_spectral_elbow_bulk(sigmas, k_min):
    """Spectral elbow targeted at the bulk transition (not the V_q plateaus).

    Returns the k beyond k_min where the relative gap is largest. This
    captures the "spike block ends" boundary if there is one.
    """
    s = np.asarray(sigmas, dtype=np.float64)
    if k_min >= len(s) - 2:
        return None, None
    rel = (s[k_min:-1] - s[k_min + 1 :]) / s[k_min:-1]
    k_idx = int(np.argmax(rel)) + k_min
    return k_idx, float(rel[k_idx - k_min])


def k_star_first_bulk_plateau(sigmas, plateau_tol=0.01, min_plateau_len=10):
    """k where the relative variation flattens to a plateau for >=min_plateau_len.

    Captures the bulk noise-floor onset.
    """
    s = np.asarray(sigmas, dtype=np.float64)
    for k in range(2, len(s) - min_plateau_len):
        window = s[k : k + min_plateau_len]
        rel_spread = (window.max() - window.min()) / window.mean()
        if rel_spread < plateau_tol:
            return k - 1
    return len(s) - 1


def k_star_sector_completion(sigmas, primes=(3, 5, 7, 11, 13)):
    """Cumulative phi(q) up to a given prime q (excluding q=2 which is in sigma_0)."""
    cum = 0
    sectors = {}
    for q in primes:
        cum += q - 1  # phi(q) = q-1 for prime q
        sectors[q] = cum
    return sectors


def analyse_d(d, payload, mp_pi_N=None):
    """Run all k_* rules on a given d's sigma spectrum."""
    sigmas = payload["sigmas"]
    M = payload["shape"][0]
    pi_N = mp_pi_N if mp_pi_N is not None else int(np.round(np.sum(np.asarray(sigmas) ** 2)))
    print(f"\n=== d={d}, N={payload['N']}, M={M}, pi(N)~{pi_N} ===")

    if "frobenius_sq" in payload:
        print(f"  saved Frobenius_sq = {payload['frobenius_sq']:.1f}")

    # R0 — canonical
    k0 = payload.get("k_star_assumed") or payload.get("k_star_canonical")
    if k0 is None:
        # synthesise via S82 rule for d=24: k_*(d) = round(exp(0.275*d - 2.24))
        k0 = int(round(math.exp(0.275 * d - 2.24)))
    block0 = cum_block(sigmas, k0)
    print(f"  R0 canonical k_* = {k0}: spike_block = {block0:.1f}, frac = {block0 / pi_N:.4f}")

    # R1 — MP edge
    k1, edge = k_star_mp_edge(sigmas, M, pi_N)
    block1 = cum_block(sigmas, k1)
    print(f"  R1 MP edge {edge:.2f} -> k_* = {k1}: spike_block = {block1:.1f}, frac = {block1 / pi_N:.4f}")

    # R2 — spectral elbow in [k_min..k_max], k_min set to skip V_3/V_5 plateaus
    if len(sigmas) >= 5:
        # Skip past first three clusters: V_3 (2 chars), V_5 (4 chars), V_7 (6 chars)
        k_min_bulk = 2 + 4 + 6  # = 12
        k2, gap = k_star_spectral_elbow_bulk(sigmas, k_min=k_min_bulk)
        if k2 is not None:
            block2 = cum_block(sigmas, k2)
            print(f"  R2 bulk-transition elbow gap {gap:.4f} at k = {k2}: spike_block = {block2:.1f}, frac = {block2 / pi_N:.4f}")
        else:
            block2 = None
    else:
        k2, block2, gap = None, None, None

    # R3 — first plateau onset
    if len(sigmas) >= 15:
        k3 = k_star_first_bulk_plateau(sigmas, plateau_tol=0.02, min_plateau_len=8)
        block3 = cum_block(sigmas, k3)
        print(f"  R3 first plateau onset k = {k3}: spike_block = {block3:.1f}, frac = {block3 / pi_N:.4f}")
    else:
        k3, block3 = None, None

    # R4 — sector completion (cumulative phi(q) up to fixed primes)
    sectors = k_star_sector_completion(sigmas)
    print(f"  R4 sector completion (cumulative phi):")
    sector_results = {}
    for q, k4 in sectors.items():
        if k4 < len(sigmas):
            block4 = cum_block(sigmas, k4)
            sector_results[q] = {"k_star": k4, "spike_block": block4, "frac": block4 / pi_N}
            print(f"    up to q={q}: k_* = {k4}, spike_block = {block4:.1f}, frac = {block4 / pi_N:.4f}")

    return {
        "N": payload["N"],
        "M": M,
        "pi_N": pi_N,
        "R0": {"k_star": k0, "spike_block": block0, "frac": block0 / pi_N},
        "R1": {"k_star": k1, "edge": edge, "spike_block": block1, "frac": block1 / pi_N},
        "R2": ({"k_star": k2, "gap": gap, "spike_block": block2, "frac": block2 / pi_N} if k2 is not None else None),
        "R3": ({"k_star": k3, "spike_block": block3, "frac": block3 / pi_N} if k3 is not None else None),
        "R4": sector_results,
    }


def main():
    out = {"by_d": {}, "headline": {}}

    print("Loading S82 saved sigmas (d=14, 18, 20)...")
    s82 = load_s82_sigmas()
    print("Loading S183 saved sigmas (d=24)...")
    s183 = load_s183_sigmas()

    # For d in {14, 18, 20}, recompute FULL sigmas to get a proper spectrum
    # (the saved data only has the top 25/40/50). This lets MP-edge, elbow,
    # and plateau rules see beyond the spike region.
    for d in (14, 18, 20):
        print(f"\nRecomputing full SVD at d={d}...")
        N, M, full_sigmas = full_svd_sigmas(d)
        # compose payload using full sigma list but keep S82's k_star_assumed
        s82[d]["sigmas_full"] = full_sigmas.tolist()
        s82[d]["pi_N"] = int(np.sum(np.asarray([1 for p in primerange(2, N + 1)])))
        # Use full sigmas for analysis
        full_payload = {
            "N": N,
            "shape": [M, M],
            "sigmas": full_sigmas.tolist(),
            "k_star_assumed": s82[d]["k_star_assumed"],
            "frobenius_sq": float(np.sum(full_sigmas ** 2)),
        }
        out["by_d"][d] = analyse_d(d, full_payload, mp_pi_N=s82[d]["pi_N"])

    # d=24 from S183 — explicit pi_N from saved Frobenius_sq
    out["by_d"][24] = analyse_d(24, s183[24], mp_pi_N=s183[24]["pi_N"])

    # Summary across rules
    print("\n\n=== SUMMARY: spike_block / pi(N) by rule and d ===")
    print(f"{'d':<4} {'R0_canon':<10} {'R1_MP':<10} {'R2_elbow':<10} {'R3_plateau':<11} {'R4_q11':<10}")
    for d in (14, 18, 20, 24):
        row = out["by_d"][d]
        cells = [
            f"{row['R0']['frac']:.4f}",
            f"{row['R1']['frac']:.4f}",
            f"{row['R2']['frac']:.4f}" if row.get("R2") else "—",
            f"{row['R3']['frac']:.4f}" if row.get("R3") else "—",
            f"{row['R4'].get(11, {}).get('frac', '—'):.4f}" if isinstance(row['R4'].get(11), dict) else "—",
        ]
        print(f"{d:<4} " + " ".join(f"{c:<10}" for c in cells))

    # Compute spread per d
    print("\n=== SPREAD across rules (max - min frac) ===")
    for d in (14, 18, 20, 24):
        row = out["by_d"][d]
        fracs = [row["R0"]["frac"], row["R1"]["frac"]]
        if row.get("R2"): fracs.append(row["R2"]["frac"])
        if row.get("R3"): fracs.append(row["R3"]["frac"])
        for q, val in row["R4"].items():
            fracs.append(val["frac"])
        spread = max(fracs) - min(fracs)
        print(f"d={d}: rule fracs = [{min(fracs):.4f}, {max(fracs):.4f}], spread = {spread:.4f}")

    out_file = os.path.join(
        os.path.dirname(__file__), "cluster_boundary_kstar_results.json"
    )
    # Drop very large arrays from output JSON
    for d in (14, 18, 20):
        if "sigmas_full" in s82[d]:
            del s82[d]["sigmas_full"]
    with open(out_file, "w") as f:
        json.dump(out, f, indent=2, default=float)
    print(f"\nWrote results to {out_file}")


if __name__ == "__main__":
    main()
