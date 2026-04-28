"""
S186 (verify-16 of S169): character-alignment k_* probe.

Adversarial direction not tried in any of S170-S185:

The S82 thesis is that spike eigenvectors of chi_P MPS unfolding ARE
Dirichlet-character vectors of squarefree q. The S169 spike-block fraction
counts the top k_* sigmas, with k_* chosen by canonical S82 extrapolation.
S185 added a Marchenko-Pastur-edge alternative rule, finding rule
divergence with ~0.10 spread in extrapolated asymptotes.

But neither S169's canonical rule nor S185's MP-edge rule USES the
S82-thesis criterion directly: "the eigenvector is character-aligned."

Direct rule: k_char(d) = max k such that the k-th right-singular-vector
has dominant centered residue energy (over Z/qZ for some q) > τ.

If S82's character-vector identification is intrinsic to the spike block,
then k_char ≈ k_canonical and the fractions should agree. If they
disagree substantially, then S169 is counting additional sigmas that are
NOT cleanly character-aligned, and the "spike block = character-vector
energy" thesis weakens at the canonical boundary.

The saved S82 data (`spike_d{14,18,20}_results.json`) gives us exactly
this — `dom_q_centered_energy` for each spike's right singular vector.

This script:
  1. Reads dom_q_centered_energy for each spike at d=14, 18, 20.
  2. For thresholds τ ∈ {0.01, 0.02, 0.05}, finds k_char(d, τ).
  3. Computes spike block sum and frac = block / π(N) at each (d, τ).
  4. Compares to S169 canonical R0 and S185 MP-edge R1.

WHAT WOULD FALSIFY THIS DIRECTION (i.e. confirm S169):
  k_char(d, τ=0.02) ≈ k_canonical(d) for all τ ∈ [0.01, 0.05] and the
  spike-block fractions agree to within 0.01.

WHAT WOULD REFUTE S169 (or narrow its scope):
  k_char(d) << k_canonical(d). The "character cliff" sits below the
  canonical k_*, meaning S169 counts sigmas that are NOT cleanly
  character-aligned. The 21% prediction holds only for the *combined*
  sum (character + post-cliff), not for the cleanly-character-aligned
  block alone.
"""

import json
import os

import numpy as np
from sympy import primerange

DATA_ROOT = "/apps/aplikacijos/prime-research/experiments/constructions"
S82_DIR = os.path.join(DATA_ROOT, "spike_eigenvectors_chi_p")


def pi_N(N):
    return sum(1 for _ in primerange(2, N + 1))


def load_spike_payload(d):
    with open(os.path.join(S82_DIR, f"spike_d{d}_results.json")) as f:
        payload = json.load(f)
    if isinstance(payload, list):
        payload = payload[0]
    return payload


def k_char_threshold(spikes, tau, k_min=1):
    """
    Largest k such that ALL spikes in [k_min, k] have dom_q_E > tau.
    Equivalently: the last "run" of character-aligned spikes from k_min.
    Stops at the first spike with dom_q_E <= tau (the "character cliff").
    """
    last_aligned = k_min - 1
    for k in range(k_min, len(spikes)):
        s = spikes[k]
        if s["dom_q_centered_energy"] > tau:
            last_aligned = k
        else:
            break
    return last_aligned


def k_char_count(spikes, tau, k_min=1):
    """
    Permissive: k_char = (k_min - 1) + count of spikes in [k_min, K] with
    dom_q_E > tau. Allows gaps; reports total count of aligned spikes.
    """
    cnt = sum(1 for s in spikes[k_min:] if s["dom_q_centered_energy"] > tau)
    return k_min - 1 + cnt


def spike_block_sum(sigmas, k_star, exclude_k0=True):
    if exclude_k0:
        return float(np.sum(np.asarray(sigmas[1 : k_star + 1]) ** 2))
    return float(np.sum(np.asarray(sigmas[: k_star + 1]) ** 2))


def main():
    out = {"by_d": {}}

    for d in (14, 18, 20):
        payload = load_spike_payload(d)
        N = payload["N"]
        spikes = payload["chi_p_spikes"]
        sigmas = [s["sigma"] for s in spikes]
        k_canon = payload["k_star_assumed"]
        pN = pi_N(N)

        print(f"\n=== d={d}, N={N}, π(N)={pN}, k_canonical={k_canon} ===")
        print(f"  Total spikes saved: {len(spikes)}")

        # Print the character-cliff structure
        print("  Spike-by-spike alignment:")
        for k in range(min(len(spikes), 30)):
            s = spikes[k]
            mark = " <-- CANONICAL" if k == k_canon else ""
            print(
                f"    k={k:3d}, σ={s['sigma']:7.3f}, σ²={s['sigma']**2:9.2f}, "
                f"dom_q={s['dom_q']:3d}, dom_q_E={s['dom_q_centered_energy']:.4f}, "
                f"min_q_cond={s['min_q_conductor']:3d}{mark}"
            )

        # Compute fractions at canonical k_*
        block_canon = spike_block_sum(sigmas, k_canon)
        frac_canon = block_canon / pN
        print(f"\n  R0 canonical k_*={k_canon}: spike_block={block_canon:.2f}, frac={frac_canon:.4f}")

        # Compute fractions at character-cliff cutoffs (consecutive-run rule)
        d_results = {
            "N": N,
            "pi_N": pN,
            "k_canonical": k_canon,
            "frac_R0_canonical": frac_canon,
            "block_R0_canonical": block_canon,
            "thresholds": {},
        }
        for tau in (0.005, 0.01, 0.02, 0.05, 0.10):
            k_run = k_char_threshold(spikes, tau)
            k_cnt = k_char_count(spikes, tau)
            block_run = spike_block_sum(sigmas, k_run) if k_run >= 1 else 0.0
            block_cnt = spike_block_sum(sigmas, k_cnt) if k_cnt >= 1 else 0.0
            frac_run = block_run / pN
            frac_cnt = block_cnt / pN
            print(
                f"  τ={tau:.3f}: k_run={k_run} (consecutive cliff), frac={frac_run:.4f} | "
                f"k_cnt={k_cnt} (count rule), frac={frac_cnt:.4f}"
            )
            d_results["thresholds"][f"{tau:.3f}"] = {
                "tau": tau,
                "k_consecutive_run": k_run,
                "frac_consecutive_run": frac_run,
                "block_consecutive_run": block_run,
                "k_count": k_cnt,
                "frac_count": frac_cnt,
                "block_count": block_cnt,
            }

        # Distance from canonical
        print(f"\n  Distance from canonical:")
        for tau in (0.01, 0.02, 0.05):
            k_run = d_results["thresholds"][f"{tau:.3f}"]["k_consecutive_run"]
            f_run = d_results["thresholds"][f"{tau:.3f}"]["frac_consecutive_run"]
            print(
                f"    τ={tau}: k_char(run)={k_run} vs k_canon={k_canon}, "
                f"Δk={k_run - k_canon}, Δfrac={f_run - frac_canon:+.4f}"
            )

        out["by_d"][d] = d_results

    # Cross-d summary
    print("\n\n=== SUMMARY: spike-block fraction by rule ===")
    print(f"{'d':<4} {'k_canon':<8} {'frac_canon':<12} {'k_τ=0.02':<10} {'frac_τ=0.02':<12} {'k_τ=0.05':<10} {'frac_τ=0.05':<12}")
    for d in (14, 18, 20):
        r = out["by_d"][d]
        kc = r["k_canonical"]
        fc = r["frac_R0_canonical"]
        k02 = r["thresholds"]["0.020"]["k_consecutive_run"]
        f02 = r["thresholds"]["0.020"]["frac_consecutive_run"]
        k05 = r["thresholds"]["0.050"]["k_consecutive_run"]
        f05 = r["thresholds"]["0.050"]["frac_consecutive_run"]
        print(f"{d:<4} {kc:<8} {fc:<12.4f} {k02:<10} {f02:<12.4f} {k05:<10} {f05:<12.4f}")

    out_file = os.path.join(
        os.path.dirname(__file__), "character_alignment_kstar_results.json"
    )
    with open(out_file, "w") as f:
        json.dump(out, f, indent=2, default=float)
    print(f"\nWrote {out_file}")


if __name__ == "__main__":
    main()
