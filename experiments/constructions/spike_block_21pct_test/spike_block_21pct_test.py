#!/usr/bin/env python3
"""S169 — Empirical test of the S168 21% spike-block-fraction prediction.

S168 prediction (paired with S74's k_*(N) ~ N^{0.42}):

    Sum over squarefree q in [2, Q*] of E(q, N) ~ 0.21 * pi(N),  Q* ~ N^{0.21}.

Equivalently, the chi_P MPS spike-block energy fraction (centered, after the
rank-1 mean) approaches 0.21 of pi(N) as N grows.

This script tests the prediction in two ways:

(A) ANALYTIC vs ANALYTIC.  Compute  cum(Q) := sum_{sqf q in [2, Q]} E(q, N)
    by direct Fourier sieve, and ask: at Q = round(N^{0.21}), is cum(Q)
    close to 0.21 * pi(N)?  Also: at what Q does cum hit 0.21*pi(N)?

(B) ANALYTIC vs SVD SPIKE BLOCK.  Load saved sigma values from
    spike_d{14,18,20}_results.json (S82), compute the spike-block sum
    sum_{k=1..k_*} sigma_k^2, and compare to cum(Q) at the matching Q*.
    The gap, if any, is the "leakage" identified in S166.

Reads:
    ../spike_eigenvectors_chi_p/spike_d14_results.json
    ../spike_eigenvectors_chi_p/spike_d18_results.json
    ../spike_eigenvectors_chi_p/spike_d20_results.json

Outputs:
    spike_block_21pct_test_results.md  (table + narrative)
    spike_block_21pct_test_results.json (machine-readable)

Edges composed: S168 (squarefree extension), S82 (spike eigenvectors),
                S74 (spike count), S166 (q in {p,2p} formula).

Falsification: the prediction FAILS if the SVD spike block sum is not
within 20% of 0.21*pi(N) for d in {14, 18, 20} (after correcting for
finite-N leakage), AND if the analytic cum(Q) at Q = round(N^{0.21}) is
not within 20% of 0.21*pi(N).
"""

import json
import math
import os
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent
SPIKE_DIR = HERE.parent / "spike_eigenvectors_chi_p"


def sieve(N):
    is_p = np.ones(N + 1, dtype=bool)
    is_p[:2] = False
    for i in range(2, int(N**0.5) + 1):
        if is_p[i]:
            is_p[i * i :: i] = False
    return is_p


def primes_up_to(N):
    return np.nonzero(sieve(N))[0]


def factor(q):
    f = []
    n = q
    p = 2
    while p * p <= n:
        if n % p == 0:
            e = 0
            while n % p == 0:
                n //= p
                e += 1
            f.append((p, e))
        p += 1
    if n > 1:
        f.append((n, 1))
    return f


def is_squarefree(q):
    return q == 1 or all(e == 1 for _, e in factor(q))


def euler_phi(q):
    if q == 1:
        return 1
    res = q
    for p, _ in factor(q):
        res = res // p * (p - 1)
    return res


def num_distinct_prime_factors(q):
    return len(factor(q))


def fourier_energy_at_q(prime_arr, q, N):
    """E(q, N) = (1/N) sum_{k coprime to q} |sum_{p' prime <= N} omega_q^{k p'}|^2."""
    energy = 0.0
    for k in range(1, q):
        if math.gcd(k, q) != 1:
            continue
        angle = 2 * np.pi * k * prime_arr / q
        c = float(np.sum(np.cos(angle)))
        s = float(np.sum(np.sin(angle)))
        energy += (c * c + s * s) / N
    return energy


def predicted_main_energy(q, pi_N, N):
    """S168 main term: mu(q)^2 * (pi(N) - r(q))^2 / (phi(q) * N)."""
    if not is_squarefree(q):
        return 0.0
    r = num_distinct_prime_factors(q)
    return (pi_N - r) ** 2 / (euler_phi(q) * N)


def load_svd_block(d):
    fn = SPIKE_DIR / f"spike_d{d}_results.json"
    if not fn.exists():
        return None
    data = json.loads(fn.read_text())[0]
    sigmas = [s["sigma"] for s in data["chi_p_spikes"]]
    return {
        "d": d,
        "N": data["N"],
        "k_star": data["k_star_assumed"],
        "n_recorded": len(sigmas),
        "sigma_0_sq": sigmas[0] ** 2,
        "spike_block_sum": float(sum(s * s for s in sigmas[1 : 1 + data["k_star_assumed"]])),
        "all_recorded_sum_excluding_mean": float(sum(s * s for s in sigmas[1:])),
        "all_recorded_sum_including_mean": float(sum(s * s for s in sigmas)),
        "sigma_per_k": sigmas,
    }


def main():
    out_md = []
    out_json = {"per_d": {}, "summary": {}}

    out_md.append("# S169 — empirical test of S168's 21% spike-block-fraction prediction\n")
    out_md.append(
        "S168 prediction (paired with S74 k_* ~ N^{0.42}):\n"
        "    sum_{sqf q in [2, Q*]} E(q, N) ~ 0.21 * pi(N),   Q* ~ N^{0.21}.\n\n"
        "This script measures the LHS for d in {14, 16, 18, 20} and compares\n"
        "to (i) 0.21*pi(N) and (ii) the SVD spike-block sum from S82.\n"
    )

    ds = [14, 16, 18, 20, 22, 24]

    out_md.append("## (A) Analytic cumulative sum vs 0.21*pi(N)\n")
    out_md.append(
        "Computes  cum(Q) = sum_{sqf q in [2, Q]} E(q, N)  by direct Fourier\n"
        "sieve (this is the LHS of S168's prediction).  Q* := round(N^{0.21}).\n"
    )
    out_md.append(
        "| d | N | pi(N) | Q*=round(N^.21) | 0.21*pi(N) | cum(Q*) emp | "
        "cum(Q*) main-term | ratio cum_emp / 0.21pi | ratio main / 0.21pi |"
    )
    out_md.append("|---|---|---|---|---|---|---|---|---|")

    cum_data_per_d = {}

    for d in ds:
        N = 2 ** d
        prime_arr = primes_up_to(N)
        pi_N = int(len(prime_arr))
        Q_max = max(60, int(round(N ** 0.30)))
        sqf = [q for q in range(2, Q_max + 1) if is_squarefree(q)]

        # Compute E(q, N) empirically for every squarefree q in [2, Q_max]
        per_q = []
        cum_emp = 0.0
        cum_main = 0.0
        for q in sqf:
            E_emp = fourier_energy_at_q(prime_arr, q, N)
            E_main = predicted_main_energy(q, pi_N, N)
            cum_emp += E_emp
            cum_main += E_main
            per_q.append((q, E_emp, E_main, cum_emp, cum_main))
        cum_data_per_d[d] = {"per_q": per_q, "pi_N": pi_N, "N": N, "Q_max": Q_max}

        # Q* = N^0.21 rounded
        Q_star = max(2, int(round(N ** 0.21)))
        cum_emp_Qstar = 0.0
        cum_main_Qstar = 0.0
        for q, E_emp, E_main, ce, cm in per_q:
            if q <= Q_star:
                cum_emp_Qstar = ce
                cum_main_Qstar = cm
        target_21 = 0.21 * pi_N
        ratio_emp = cum_emp_Qstar / target_21 if target_21 > 0 else float("nan")
        ratio_main = cum_main_Qstar / target_21 if target_21 > 0 else float("nan")
        out_md.append(
            f"| {d} | {N} | {pi_N} | {Q_star} | {target_21:.2f} | "
            f"{cum_emp_Qstar:.2f} | {cum_main_Qstar:.2f} | "
            f"{ratio_emp:.4f} | {ratio_main:.4f} |"
        )

        out_json["per_d"][d] = {
            "N": N,
            "pi_N": pi_N,
            "Q_star": Q_star,
            "target_21": target_21,
            "cum_emp_Qstar": cum_emp_Qstar,
            "cum_main_Qstar": cum_main_Qstar,
            "ratio_cum_emp_to_21": ratio_emp,
            "ratio_cum_main_to_21": ratio_main,
            "Q_max_swept": Q_max,
        }

    out_md.append("")

    # Cumulative table at canonical Q exponents
    out_md.append("## (A.1) cum(Q) / pi(N) at canonical Q = round(N^a) for a in {0.10, 0.15, 0.20, 0.21, 0.25, 0.30}\n")
    out_md.append("| d | N | pi(N) | a=0.10 | a=0.15 | a=0.20 | a=0.21 | a=0.25 | a=0.30 |")
    out_md.append("|---|---|---|---|---|---|---|---|---|")
    canonical_table = {}
    for d in ds:
        N = 2 ** d
        pi_N = cum_data_per_d[d]["pi_N"]
        per_q = cum_data_per_d[d]["per_q"]
        row_parts = [f"| {d} | {N} | {pi_N} |"]
        canonical_table[d] = {}
        for a in [0.10, 0.15, 0.20, 0.21, 0.25, 0.30]:
            Q = max(2, int(round(N ** a)))
            cum = 0.0
            for q, E_emp, E_main, ce, cm in per_q:
                if q <= Q:
                    cum = ce
            frac = cum / pi_N
            row_parts.append(f" {frac:.4f} |")
            canonical_table[d][a] = {"Q": Q, "cum": cum, "frac": frac}
        out_md.append("".join(row_parts))
    out_md.append("")
    out_json["canonical_table"] = canonical_table

    # Find the Q at which cum/pi(N) crosses 0.21
    out_md.append("## (A.2) crossing Q at which cum(Q) / pi(N) reaches 0.21\n")
    out_md.append("| d | N | pi(N) | Q_cross | log Q_cross / log N |")
    out_md.append("|---|---|---|---|---|")
    crossing = {}
    for d in ds:
        N = 2 ** d
        pi_N = cum_data_per_d[d]["pi_N"]
        per_q = cum_data_per_d[d]["per_q"]
        target = 0.21 * pi_N
        Q_cross = None
        for q, E_emp, E_main, ce, cm in per_q:
            if ce >= target:
                Q_cross = q
                break
        exp = math.log(Q_cross) / math.log(N) if Q_cross else float("nan")
        crossing[d] = {"Q_cross": Q_cross, "exp": exp}
        Q_cross_str = str(Q_cross) if Q_cross is not None else "(>Q_max)"
        exp_str = f"{exp:.4f}" if Q_cross else "n/a"
        out_md.append(f"| {d} | {N} | {pi_N} | {Q_cross_str} | {exp_str} |")
    out_md.append("")
    out_json["crossing"] = crossing

    # (B) SVD spike block comparison
    out_md.append("## (B) SVD spike block sum vs 0.21*pi(N) and vs analytic cum(Q*)\n")
    out_md.append(
        "From S82 saved JSONs: spike block = sum_{k=1..k_*} sigma_k^2 (excluding rank-1 mean).\n"
    )
    out_md.append(
        "| d | k_* | spike block | 0.21*pi(N) | cum(Q*) emp | "
        "spike/0.21pi | spike/cum(Q*) emp | leakage % |"
    )
    out_md.append("|---|---|---|---|---|---|---|---|")
    svd_results = {}
    for d in ds:
        svd = load_svd_block(d)
        if svd is None:
            out_md.append(f"| {d} | (no JSON) | - | - | - | - | - | - |")
            continue
        N = 2 ** d
        pi_N = cum_data_per_d[d]["pi_N"]
        target_21 = 0.21 * pi_N
        cum_emp = out_json["per_d"][d]["cum_emp_Qstar"]
        spike = svd["spike_block_sum"]
        ratio_to_21 = spike / target_21 if target_21 > 0 else float("nan")
        ratio_to_analytic = spike / cum_emp if cum_emp > 0 else float("nan")
        leakage_pct = (
            100.0 * (spike - cum_emp) / cum_emp if cum_emp > 0 else float("nan")
        )
        svd_results[d] = {
            "k_star": svd["k_star"],
            "spike_block_sum": spike,
            "ratio_to_21pi": ratio_to_21,
            "ratio_to_analytic_Qstar": ratio_to_analytic,
            "leakage_pct": leakage_pct,
        }
        out_md.append(
            f"| {d} | {svd['k_star']} | {spike:.2f} | {target_21:.2f} | "
            f"{cum_emp:.2f} | {ratio_to_21:.4f} | {ratio_to_analytic:.4f} | "
            f"{leakage_pct:+.2f}% |"
        )
    out_md.append("")
    out_json["svd_block"] = svd_results

    # (B.1) Test the "right" comparison: at what Q does cum(Q) match the SVD spike block?
    out_md.append("## (B.1) Q at which cum(Q) matches the SVD spike block\n")
    out_md.append(
        "If the spike block is the V_q^prim direct sum for some `effective Q*`,\n"
        "then cum(Q_eff) should match the SVD spike block.  Solve for Q_eff.\n"
    )
    out_md.append("| d | spike block | Q_eff | log Q_eff / log N | predicted = 0.21? |")
    out_md.append("|---|---|---|---|---|")
    Q_eff_table = {}
    for d in ds:
        if d not in svd_results:
            continue
        spike = svd_results[d]["spike_block_sum"]
        per_q = cum_data_per_d[d]["per_q"]
        N = 2 ** d
        Q_eff = None
        for q, E_emp, E_main, ce, cm in per_q:
            if ce >= spike:
                Q_eff = q
                break
        exp_eff = math.log(Q_eff) / math.log(N) if Q_eff else float("nan")
        match_021 = abs(exp_eff - 0.21) < 0.05 if Q_eff else False
        Q_eff_table[d] = {"Q_eff": Q_eff, "exp": exp_eff}
        Q_eff_str = str(Q_eff) if Q_eff is not None else "(>Q_max)"
        exp_eff_str = f"{exp_eff:.4f}" if Q_eff else "n/a"
        out_md.append(
            f"| {d} | {spike:.2f} | {Q_eff_str} | {exp_eff_str} | "
            f"{'YES' if match_021 else 'NO'} |"
        )
    out_md.append("")
    out_json["Q_eff_table"] = Q_eff_table

    # (B.2) Cumulative trajectory, save to log
    out_md.append("## (B.2) Full cum(Q) trajectory at d=20 (last 30 rows)\n")
    out_md.append("| q | E_emp | cum_emp | cum_emp / pi(N) | leadup |")
    out_md.append("|---|---|---|---|---|")
    per_q = cum_data_per_d[20]["per_q"]
    pi_N = cum_data_per_d[20]["pi_N"]
    for q, E_emp, E_main, ce, cm in per_q[-30:]:
        out_md.append(f"| {q} | {E_emp:.4f} | {ce:.4f} | {ce/pi_N:.4f} | |")
    out_md.append("")

    # Save
    (HERE / "spike_block_21pct_test_results.md").write_text("\n".join(out_md))
    (HERE / "spike_block_21pct_test_results.json").write_text(
        json.dumps(out_json, indent=2, default=lambda o: float(o) if isinstance(o, np.floating) else None)
    )

    # Stdout summary
    print("\n".join(out_md))


if __name__ == "__main__":
    main()
