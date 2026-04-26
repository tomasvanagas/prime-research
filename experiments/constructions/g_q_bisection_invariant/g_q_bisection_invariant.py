"""
Construction C1 (composing E1.6 + E1.5):
    g_q(x) := ( A(x) mod q, C_3(x) mod q ) where
      A(x)   = (x - L(x)) / 2 = #{n <= x : Omega(n) odd}
      C_3(x) = A(x) - pi(x)   = #{n <= x : Omega(n) odd, >= 3}.

Three coupled predictions are tested:

    PR1 (per-component per-step rate)
        H( F mod q | F(x-1) mod q ) = h_2( F(X) / X ) + O( 1 / F(X) )
        for F in { A, C_3, pi }, q in {2, 3, 5, 7, 11, 13}, X up to 2 * 10^6.

    PR2 (joint per-step rate)
        H( g_q(x) | g_q(x-1) ) = H_3( 1 - A(X)/X, pi(X)/X, C_3(X)/X )
                                  + O( 1 / pi(X) )
        in regime q^2 << pi(X) (state-space coverage).

    PR3 (q-stable marginal independence, generalising E1.6)
        I( A mod q ; C_3 mod q )  small for all q in {2, 3, 5, 7, 11, 13}.

See `definition.md` and `g_q_bisection_invariant_results.md` for detailed
falsification criteria.

Usage
-----
python g_q_bisection_invariant.py [--N 2_000_000] [--quick]

Quick mode runs N = 200_000 in <1 s for fast iteration. Default N = 2 * 10^6
takes ~10 s.
"""

from __future__ import annotations

import argparse
import json
import math
import time

import numpy as np


# ---------- Sieves -------------------------------------------------------


def sieve_omega(N: int):
    """Return Omega[0..N]: Omega[n] = number of prime factors of n with
    multiplicity. Omega[0] = Omega[1] = 0. Uses smallest-prime-factor sieve."""
    spf = np.zeros(N + 1, dtype=np.int32)
    for p in range(2, N + 1):
        if spf[p] == 0:
            slc = spf[p::p]
            slc[slc == 0] = p
            spf[p::p] = slc
    Omega = np.zeros(N + 1, dtype=np.int32)
    for n in range(2, N + 1):
        Omega[n] = Omega[n // spf[n]] + 1
    return Omega


def build_summatories(Omega: np.ndarray):
    """Return dict { 'A', 'C3', 'pi', 'L' }: cumulative arrays indexed 0..N."""
    N = len(Omega) - 1
    lam = np.where(Omega % 2 == 1, -1, 1).astype(np.int64)
    lam[0] = 0  # n=0 excluded
    L = np.cumsum(lam).astype(np.int64)

    A = ((np.arange(N + 1, dtype=np.int64) - L) // 2).astype(np.int64)

    pi = np.zeros(N + 1, dtype=np.int64)
    pi_step = (Omega == 1).astype(np.int64)
    pi_step[0] = 0
    pi = np.cumsum(pi_step)

    C3 = A - pi

    return {"A": A, "C3": C3, "pi": pi, "L": L, "lam": lam}


# ---------- Information measures -----------------------------------------


def conditional_entropy_pair(prev: np.ndarray, curr: np.ndarray, q: int) -> float:
    """H( curr mod q | prev mod q ), bits.
    `prev` and `curr` are integer arrays of equal length n; we compute
    the empirical conditional entropy treating each (prev[i], curr[i])
    pair as a sample.
    """
    p = (prev % q).astype(np.int64)
    c = (curr % q).astype(np.int64)
    joint = p * q + c  # encode pair as single integer in [0, q^2)
    counts_joint = np.bincount(joint, minlength=q * q)
    counts_prev = counts_joint.reshape(q, q).sum(axis=1)
    n = counts_joint.sum()
    H = 0.0
    for s in range(q):
        if counts_prev[s] == 0:
            continue
        p_s = counts_prev[s] / n
        # P(c | s) for each c with count > 0
        for ck in range(q):
            j = s * q + ck
            cnt = counts_joint[j]
            if cnt == 0:
                continue
            p_cs = cnt / counts_prev[s]
            H -= p_s * p_cs * math.log2(p_cs)
    return H


def joint_conditional_entropy(
    A: np.ndarray, C3: np.ndarray, q: int
) -> float:
    """H( (A, C_3)(x) mod q | (A, C_3)(x-1) mod q ), bits."""
    n = len(A) - 1
    a_curr = (A[1 : n + 1] % q).astype(np.int64)
    c_curr = (C3[1 : n + 1] % q).astype(np.int64)
    a_prev = (A[0:n] % q).astype(np.int64)
    c_prev = (C3[0:n] % q).astype(np.int64)
    # Encode pair as integer in [0, q^2)
    state_curr = a_curr * q + c_curr
    state_prev = a_prev * q + c_prev
    Q = q * q
    joint = state_prev * Q + state_curr  # in [0, q^4)
    counts_joint = np.bincount(joint, minlength=Q * Q)
    counts_prev = counts_joint.reshape(Q, Q).sum(axis=1)
    total = counts_joint.sum()
    H = 0.0
    for s in range(Q):
        if counts_prev[s] == 0:
            continue
        p_s = counts_prev[s] / total
        for ck in range(Q):
            cnt = counts_joint[s * Q + ck]
            if cnt == 0:
                continue
            p_cs = cnt / counts_prev[s]
            H -= p_s * p_cs * math.log2(p_cs)
    return H


def marginal_mutual_information(
    A: np.ndarray, C3: np.ndarray, q: int
) -> float:
    """I( A mod q ; C_3 mod q ) over x in [1, N], bits."""
    a = (A[1:] % q).astype(np.int64)
    c = (C3[1:] % q).astype(np.int64)
    joint = a * q + c
    counts = np.bincount(joint, minlength=q * q).astype(np.float64)
    total = counts.sum()
    P_joint = counts.reshape(q, q) / total
    P_a = P_joint.sum(axis=1)
    P_c = P_joint.sum(axis=0)
    eps = 1e-300
    H_a = -np.sum(P_a * np.log2(P_a + eps))
    H_c = -np.sum(P_c * np.log2(P_c + eps))
    H_joint = -np.sum(P_joint * np.log2(P_joint + eps))
    return float(H_a + H_c - H_joint)


def h_binary(p: float) -> float:
    if p <= 0 or p >= 1:
        return 0.0
    return -p * math.log2(p) - (1 - p) * math.log2(1 - p)


def h_three(p1: float, p2: float, p3: float) -> float:
    s = 0.0
    for p in (p1, p2, p3):
        if p > 0:
            s -= p * math.log2(p)
    return s


# ---------- Experiment driver --------------------------------------------


def run_experiment(N: int, q_list, X_list):
    print(f"Sieving Omega for N = {N:,} ...")
    t0 = time.time()
    Omega = sieve_omega(N)
    print(f"  done in {time.time() - t0:.2f} s")

    summ = build_summatories(Omega)
    A, C3, pi = summ["A"], summ["C3"], summ["pi"]

    # Sanity: pi(x) = A(x) - C_3(x) for all x
    assert np.array_equal(pi, A - C3), "Bisection identity broken!"
    print("  bisection identity pi = A - C_3 verified bit-exact.")

    results = {
        "N": N,
        "q_list": q_list,
        "X_list": X_list,
        "summary": {},
        "PR1_per_step": {},
        "PR2_joint_step": {},
        "PR3_marginal_MI": {},
    }

    for X in X_list:
        if X > N:
            continue
        pi_X = int(pi[X])
        A_X = int(A[X])
        C3_X = int(C3[X])
        rho_pi = pi_X / X
        rho_A = A_X / X
        rho_C3 = C3_X / X
        rho_even = 1.0 - rho_A
        h_pi = h_binary(rho_pi)
        h_A = h_binary(rho_A)
        h_C3 = h_binary(rho_C3)
        h_joint_pred = h_three(rho_even, rho_pi, rho_C3)

        results["summary"][str(X)] = {
            "pi_X": pi_X,
            "A_X": A_X,
            "C3_X": C3_X,
            "pi_density": rho_pi,
            "A_density": rho_A,
            "C3_density": rho_C3,
            "Omega_even_density": rho_even,
            "h_pi_pred_bits": h_pi,
            "h_A_pred_bits": h_A,
            "h_C3_pred_bits": h_C3,
            "h_joint_pred_bits": h_joint_pred,
        }

    # PR1: per-component per-step conditional entropy.
    # We compute H over the prefix x in [1, X] for each X in X_list.
    print("\nPR1: per-component per-step conditional entropies")
    print("  (H_emp - H_pred should be ~0 in regime q << density * X)")
    for X in X_list:
        if X > N:
            continue
        results["PR1_per_step"][str(X)] = {}
        s = results["summary"][str(X)]
        for q in q_list:
            entry = {}
            # Compute H( F(x) mod q | F(x-1) mod q ) for x in [1, X],
            # using prev = F[0..X-1], curr = F[1..X].
            for F_name, F in (("pi", pi), ("A", A), ("C3", C3)):
                H_emp = conditional_entropy_pair(F[0:X], F[1 : X + 1], q)
                H_pred = s[f"h_{F_name}_pred_bits"]
                entry[F_name] = {
                    "H_emp": H_emp,
                    "H_pred": H_pred,
                    "diff": H_emp - H_pred,
                }
            results["PR1_per_step"][str(X)][str(q)] = entry

    # PR2: joint per-step conditional entropy of (A, C3) mod q.
    print("PR2: joint per-step conditional entropy H(g_q(x)|g_q(x-1))")
    for X in X_list:
        if X > N:
            continue
        results["PR2_joint_step"][str(X)] = {}
        s = results["summary"][str(X)]
        for q in q_list:
            H_emp = joint_conditional_entropy(A[0 : X + 1], C3[0 : X + 1], q)
            H_pred = s["h_joint_pred_bits"]
            results["PR2_joint_step"][str(X)][str(q)] = {
                "H_emp": H_emp,
                "H_pred": H_pred,
                "diff": H_emp - H_pred,
            }

    # PR3: marginal mutual information I(A mod q ; C_3 mod q) on x in [1, X].
    print("PR3: marginal mutual information I(A mod q ; C_3 mod q)")
    for X in X_list:
        if X > N:
            continue
        results["PR3_marginal_MI"][str(X)] = {}
        for q in q_list:
            I = marginal_mutual_information(A[: X + 1], C3[: X + 1], q)
            results["PR3_marginal_MI"][str(X)][str(q)] = I

    return results


# ---------- Reporting -----------------------------------------------------


def print_report(results):
    qs = results["q_list"]
    Xs = [int(x) for x in results["X_list"] if int(x) <= results["N"]]

    print("\n=====  SUMMARY (densities and closed-form predictions)  =====")
    print(
        f"{'X':>10} {'pi_X':>8} {'A_X':>10} {'C3_X':>10} "
        f"{'rho_pi':>9} {'rho_A':>8} {'rho_C3':>8} "
        f"{'h_pi':>7} {'h_A':>7} {'h_C3':>7} {'h_joint':>8}"
    )
    for X in Xs:
        s = results["summary"][str(X)]
        print(
            f"{X:>10} {s['pi_X']:>8} {s['A_X']:>10} {s['C3_X']:>10} "
            f"{s['pi_density']:>9.5f} {s['A_density']:>8.5f} "
            f"{s['C3_density']:>8.5f} "
            f"{s['h_pi_pred_bits']:>7.4f} {s['h_A_pred_bits']:>7.4f} "
            f"{s['h_C3_pred_bits']:>7.4f} {s['h_joint_pred_bits']:>8.4f}"
        )

    print("\n=====  PR1: |H_emp - H_pred| for each (X, q, F)  =====")
    print("(Falsification: any |diff| > 5e-3 in regime q << pi(X)/100 fails)")
    for X in Xs:
        print(f"\nX = {X}")
        print(
            f"  {'q':>3}    "
            f"{'pi diff':>11}  {'A diff':>11}  {'C3 diff':>11}"
        )
        for q in qs:
            e = results["PR1_per_step"][str(X)][str(q)]
            print(
                f"  {q:>3}    "
                f"{e['pi']['diff']:>+11.6f}  "
                f"{e['A']['diff']:>+11.6f}  "
                f"{e['C3']['diff']:>+11.6f}"
            )

    print("\n=====  PR2: |H_emp - H_pred| for joint (A,C3) mod q  =====")
    print("(Falsification: |diff| > 5e-3 in regime q^2 << pi(X) fails)")
    for X in Xs:
        print(f"\nX = {X}   pi(X)/q^2 cutoff at q^2 << {results['summary'][str(X)]['pi_X']}")
        print(f"  {'q':>3}  {'q^2':>5}  {'H_emp':>9}  {'H_pred':>9}  {'diff':>11}")
        for q in qs:
            e = results["PR2_joint_step"][str(X)][str(q)]
            print(
                f"  {q:>3}  {q*q:>5}  "
                f"{e['H_emp']:>9.5f}  {e['H_pred']:>9.5f}  "
                f"{e['diff']:>+11.6f}"
            )

    print("\n=====  PR3: I(A mod q ; C_3 mod q) marginal in bits  =====")
    print("(Falsification: I > 0.01 bits at X = 2*10^6 for any q in {2..13})")
    print(f"  {'q':>3}  ", end="")
    for X in Xs:
        print(f"X={X:<10}", end=" ")
    print()
    for q in qs:
        print(f"  {q:>3}  ", end="")
        for X in Xs:
            I = results["PR3_marginal_MI"][str(X)][str(q)]
            print(f"{I:>10.6f}  ", end="")
        print()

    # Verdicts
    print("\n=====  VERDICTS  =====")
    pr1_max = 0.0
    pr1_max_loc = None
    for X in Xs:
        # restrict to q with q^2 << pi(X), but for PR1 we just need q << pi(X)
        if results["summary"][str(X)]["pi_X"] < 100:
            continue
        for q in qs:
            if q * 100 > results["summary"][str(X)]["pi_X"]:
                continue
            for F in ("pi", "A", "C3"):
                d = abs(
                    results["PR1_per_step"][str(X)][str(q)][F]["diff"]
                )
                if d > pr1_max:
                    pr1_max = d
                    pr1_max_loc = (X, q, F)
    print(f"PR1 worst |diff| (in regime q*100 < pi(X)): {pr1_max:.6f}  at  {pr1_max_loc}")
    print(f"  PR1 verdict: {'PASS' if pr1_max < 5e-3 else 'FAIL'} (threshold 5e-3)")

    pr2_max = 0.0
    pr2_max_loc = None
    for X in Xs:
        if results["summary"][str(X)]["pi_X"] < 100:
            continue
        for q in qs:
            if q * q * 100 > results["summary"][str(X)]["pi_X"]:
                continue
            d = abs(results["PR2_joint_step"][str(X)][str(q)]["diff"])
            if d > pr2_max:
                pr2_max = d
                pr2_max_loc = (X, q)
    print(f"PR2 worst |diff| (in regime q^2*100 < pi(X)): {pr2_max:.6f}  at  {pr2_max_loc}")
    print(f"  PR2 verdict: {'PASS' if pr2_max < 5e-3 else 'FAIL'} (threshold 5e-3)")

    X_top = max(Xs)
    pr3_max = 0.0
    pr3_max_q = None
    for q in qs:
        I = results["PR3_marginal_MI"][str(X_top)][str(q)]
        if I > pr3_max:
            pr3_max = I
            pr3_max_q = q
    print(f"PR3 worst I (at X={X_top}, all q in {qs}): {pr3_max:.6f}  at  q={pr3_max_q}")
    print(f"  PR3 verdict: {'PASS' if pr3_max < 0.01 else 'FAIL'} (threshold 0.01 bits)")


# ---------- Main ----------------------------------------------------------


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--N", type=int, default=2_000_000)
    p.add_argument("--quick", action="store_true",
                   help="Use N = 200_000 for fast iteration.")
    p.add_argument("--out", type=str, default="g_q_bisection_invariant_data.json")
    args = p.parse_args()

    if args.quick:
        N = 200_000
    else:
        N = args.N

    q_list = [2, 3, 5, 7, 11, 13]
    X_list = [10_000, 100_000, 1_000_000, 2_000_000]
    if N < 2_000_000:
        X_list = [x for x in X_list if x <= N]
        if N not in X_list:
            X_list.append(N)
        X_list.sort()

    t0 = time.time()
    results = run_experiment(N, q_list, X_list)
    elapsed = time.time() - t0
    print(f"\nElapsed: {elapsed:.2f} s")
    print_report(results)

    with open(args.out, "w") as f:
        json.dump(results, f, indent=2)
    print(f"\nWrote {args.out}")


if __name__ == "__main__":
    main()
