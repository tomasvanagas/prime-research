"""
S198 — JOINT conditional entropy of (pi(x) mod m_1, ..., pi(x) mod m_k).

Adversarial probe of E1.5's closure of "CRT reconstruction cannot win".

E1.5 measured *per-modulus* H(pi(x) mod m | pi(x-1) mod m) = h_2(pi(X)/X)
across m. The *joint* conditional entropy across k moduli simultaneously
was never measured directly. Two competing hypotheses:

  (H1) JOINT conditional entropy = sum of marginals = k * h_2(pi(X)/X).
       (Would be true if the residues were independent across moduli.)
  (H2) JOINT conditional entropy = h_2(pi(X)/X), independent of k.
       (Would be true if the residues are perfectly correlated.)

The mechanism predicts (H2): the joint transition
  (s_1, ..., s_k) -> (s_1 + b, ..., s_k + b)  with  b = 1[x prime]
is determined by a SINGLE bit b. So the JOINT new info per step is at most
h_2(P[x prime]) = h_2(pi(X)/X), regardless of k.

If (H2) holds empirically, then combining many moduli gives the SAME bits/step
as a single modulus -- this is a STRONGER statement than E1.5 currently makes.
If (H1) held, CRT-incremental could in principle gain k-fold per step, which
would be algorithmically interesting.

The outcome falsifies (H1) sharply and confirms (H2). This is the precise
information-theoretic reason that "incremental CRT cannot win": k coprime
residue chains share the SAME single-bit randomness source (the prime
indicator), not k independent bit sources.

Edges touched: E1.5 (sharpened: per-step joint info is CONSTANT in k, not
linear). Closes-under: CLOSED_PATHS row 243.
"""

from __future__ import annotations

import math
import sys
from itertools import combinations

import numpy as np


def sieve_pi_sequence(X: int) -> np.ndarray:
    """Return numpy array pi[0..X] where pi[x] = #primes <= x."""
    if X < 2:
        return np.zeros(X + 1, dtype=np.int64)
    is_prime = np.ones(X + 1, dtype=bool)
    is_prime[:2] = False
    for p in range(2, int(math.isqrt(X)) + 1):
        if is_prime[p]:
            is_prime[p * p :: p] = False
    return np.cumsum(is_prime, dtype=np.int64)


def joint_state_entropy(pi_seq: np.ndarray, moduli: tuple[int, ...]) -> float:
    """
    Compute H( joint(x) | joint(x-1) ) where
       joint(x) = (pi(x) mod m_1, ..., pi(x) mod m_k).

    Using the identity that pi(x) - pi(x-1) in {0, 1}, the joint transition
    is determined by a single bit b = 1[x prime]: joint(x) = joint(x-1) + b
    (componentwise mod m_i). So we can compute the joint conditional entropy
    by, for each visited joint state s, looking at the empirical
    P[b = 1 | joint(x-1) = s] and averaging the binary entropy h_2(.) of that.

    Returns H_joint_cond in bits/step.
    """
    # Encode joint state as a single integer: state = sum_i s_i * prod_{j<i} m_j.
    # Using the CRT-mixed-radix representation. (We don't actually need CRT;
    # any injective encoding works. Use a Python tuple key for safety with
    # arbitrary moduli.)
    n = len(pi_seq)

    # Compute joint states: (pi[x] mod m_1, ..., pi[x] mod m_k) for each x.
    # We'll bucket by the previous joint state s and count whether b=1 (prime
    # at x) or b=0 (not prime at x).
    cols = [pi_seq % m for m in moduli]
    # Hash each joint tuple to a single integer (since moduli are bounded,
    # we can use a positional encoding).
    M = 1
    state_int = np.zeros(n, dtype=np.int64)
    for col, m in zip(cols, moduli):
        state_int += col * M
        M *= m

    # b[x] = 1 if x is prime (i.e., pi[x] != pi[x-1]).
    # Step at index x in [1, n-1].
    prev_state = state_int[:-1]
    b = (pi_seq[1:] - pi_seq[:-1]).astype(np.int64)  # in {0, 1}
    assert b.min() >= 0 and b.max() <= 1, "pi grows by 0 or 1 per step"

    # For each (prev_state, b), tally counts.
    # Use np.unique to find buckets cheaply.
    n_states = M  # may be large; use a dict if too big
    # Heuristic: if M <= 10^7, use direct array; else use dict.
    if M <= 10_000_000:
        zeros = np.zeros(n_states, dtype=np.int64)
        ones = np.zeros(n_states, dtype=np.int64)
        np.add.at(zeros, prev_state[b == 0], 1)
        np.add.at(ones, prev_state[b == 1], 1)
        total = zeros + ones
        nz = total > 0
        p_state = total[nz] / total[nz].sum()
        p_one = ones[nz] / total[nz]
        # Entropy contribution
        h_per_state = np.where(
            (p_one > 0) & (p_one < 1),
            -p_one * np.log2(np.where(p_one > 0, p_one, 1)) - (1 - p_one) * np.log2(np.where(p_one < 1, 1 - p_one, 1)),
            0.0,
        )
        H_cond = float((p_state * h_per_state).sum())
    else:
        from collections import Counter

        zeros: Counter[int] = Counter()
        ones: Counter[int] = Counter()
        for s, bb in zip(prev_state.tolist(), b.tolist()):
            if bb == 1:
                ones[s] += 1
            else:
                zeros[s] += 1
        total_pairs = len(b)
        H_cond = 0.0
        for s in set(zeros) | set(ones):
            n0 = zeros[s]
            n1 = ones[s]
            tot = n0 + n1
            p_state = tot / total_pairs
            p1 = n1 / tot
            if 0 < p1 < 1:
                h = -p1 * math.log2(p1) - (1 - p1) * math.log2(1 - p1)
            else:
                h = 0.0
            H_cond += p_state * h
    return H_cond


def per_modulus_entropy(pi_seq: np.ndarray, m: int) -> float:
    """H(pi(x) mod m | pi(x-1) mod m) by the same method but k=1."""
    return joint_state_entropy(pi_seq, (m,))


def h2(p: float) -> float:
    if p <= 0 or p >= 1:
        return 0.0
    return -p * math.log2(p) - (1 - p) * math.log2(1 - p)


def report_test(pi_seq: np.ndarray, label: str, moduli: tuple[int, ...]) -> dict:
    X = len(pi_seq) - 1
    n_primes = int(pi_seq[-1])
    p_dens = n_primes / X
    h2_pred = h2(p_dens)
    H_marginal_sum = sum(per_modulus_entropy(pi_seq, m) for m in moduli)
    H_joint = joint_state_entropy(pi_seq, moduli)
    k = len(moduli)
    print(
        f"  [{label}]  X={X}  k={k}  moduli={moduli}\n"
        f"    h_2(pi/X)      = {h2_pred:.6f}\n"
        f"    sum marginal H = {H_marginal_sum:.6f}  (= k * {H_marginal_sum/k:.6f})\n"
        f"    JOINT H(.|.)   = {H_joint:.6f}\n"
        f"    H_joint / h_2  = {H_joint / h2_pred:.4f}"
        f"  (predicted 1.0 under H2; predicted {k} under H1)\n"
        f"    H_joint / sum  = {H_joint / H_marginal_sum:.4f}"
        f"  (predicted 1/k = {1/k:.4f} under H2; 1.0 under H1)"
    )
    return {
        "label": label,
        "X": X,
        "k": k,
        "moduli": list(moduli),
        "h2_pred": h2_pred,
        "H_marginal_sum": H_marginal_sum,
        "H_joint": H_joint,
        "ratio_joint_to_h2": H_joint / h2_pred,
        "ratio_joint_to_marginal_sum": H_joint / H_marginal_sum,
    }


def main():
    quick = "--quick" in sys.argv
    X_values = [10**4, 10**5, 10**6] if not quick else [10**4]

    # Various k-tuples of moduli. To probe (H1) vs (H2) we want k > 1 coprime
    # moduli, and we want the *product* small enough to fit a state-count array.
    moduli_sets = [
        (2,),  # k=1 baseline
        (2, 3),  # k=2 small primes
        (2, 3, 5),  # k=3
        (2, 3, 5, 7),  # k=4
        (2, 3, 5, 7, 11),  # k=5
        (2, 3, 5, 7, 11, 13),  # k=6 (product 30030)
        (4, 9, 25),  # k=3 prime powers
        (8, 9, 5, 7, 11),  # k=5 mixed
    ]

    results = []
    for X in X_values:
        print(f"\n[sieve] X = {X}")
        pi_seq = sieve_pi_sequence(X)
        for moduli in moduli_sets:
            label = f"X={X} k={len(moduli)}"
            res = report_test(pi_seq, label, moduli)
            results.append(res)

    # Summary: tabulate H_joint vs k * h_2 vs h_2.
    print("\n# Summary table\n")
    print("| X | k | moduli | h_2 | sum marginal | JOINT | ratio J/h_2 | ratio J/sum |")
    print("|---|---|--------|-----|-------------|-------|-------------|-------------|")
    for r in results:
        print(
            f"| {r['X']} | {r['k']} | {r['moduli']} | {r['h2_pred']:.4f} | "
            f"{r['H_marginal_sum']:.4f} | {r['H_joint']:.4f} | "
            f"{r['ratio_joint_to_h2']:.4f} | {r['ratio_joint_to_marginal_sum']:.4f} |"
        )

    # Falsification.
    print("\n# Falsification check\n")
    h1_evidence_count = 0  # (H1) would show ratio J/h_2 ~ k
    h2_evidence_count = 0  # (H2) would show ratio J/h_2 ~ 1
    for r in results:
        if r["k"] == 1:
            continue
        if abs(r["ratio_joint_to_h2"] - 1.0) < 0.005:
            h2_evidence_count += 1
        elif abs(r["ratio_joint_to_h2"] - r["k"]) < 0.05:
            h1_evidence_count += 1

    if h2_evidence_count > 0 and h1_evidence_count == 0:
        print(
            f"PASS: (H2) confirmed in {h2_evidence_count}/{len(results)-len(X_values)} k>1 cells "
            "— joint conditional entropy is CONSTANT in k, equal to h_2(pi(X)/X)."
        )
        print(
            "Implication: combining k moduli gives ZERO new bits over a single modulus.\n"
            "The k coordinates of joint(x) are perfectly correlated via the SAME "
            "single-bit prime indicator b = 1[x prime].\n"
            "This sharpens E1.5: per-step joint info is BOUNDED by h_2 from above,\n"
            "INDEPENDENT of how many moduli are combined."
        )
    elif h1_evidence_count > 0:
        print(f"FAIL: (H1) evidence in {h1_evidence_count} cells — joint info scales LINEARLY in k.")
    else:
        print(f"INCONCLUSIVE: H2 cells {h2_evidence_count}, H1 cells {h1_evidence_count}.")


if __name__ == "__main__":
    main()
