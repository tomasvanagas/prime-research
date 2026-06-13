"""
Adversarial probe of E2.13 closure (re-verify-closure mode, S237).

CLAIM ATTACKED: E2.13's closure says U^k(chi_P) = HL singular series gives
"no algorithmic content beyond HL". The S205/S191 spike approximator
    T_Q(n) := (pi(N)/N) * sum_{q sqf <= Q} mu(q)/phi(q) * c_q(n)
is polylog-evaluable and reproduces HL 2-pt correlation within 0.6%. The
re-verify probe asks: as a *pointwise* primality predictor, is T_Q(n)
better than the trivial wheel sieve W_Q at the same Q?

Key structural facts (this session's clarification):
  - T_Q(n) is a sum over SQUAREFREE q in [1, Q]. The sum does NOT factorize
    as the full Euler product over primes <= Q unless Q >= primorial of all
    primes <= Q, which fails as soon as a prime p satisfies p_2 such that
    p * p' > Q for some smaller prime p'.
  - At Q = primorial(k) (i.e., Q = 2, 6, 30, 210, 2310, ...), the truncated
    sum INCLUDES all squarefree products of primes <= p_k AND additional
    primes p_{k+1}, p_{k+2}, ..., q_max alone (without products).
  - The wheel-W identity at Q = primorial(k):
        T_Q(n) | gcd(n, primorial(k)) = 1  ~~  (pi/N) * primorial(k)/phi(primorial(k))
                                             + corrections from primes > p_k.
    S191 noted this; this session asks whether the corrections cause
    T_Q to discriminate primes BETTER than the wheel.

Falsifier:
  (F-A) AUC(T_Q) - AUC(W_Q) > 0.05 absolute at any Q       -> A-grade reopen.
  (F-B) AUC(T_Q) ~~ AUC(W_Q) within +/-0.005             -> closure confirmed.
  (F-C) T_Q better at large Q only                       -> partial reopen.
"""

import math
import time
import json
import os
from sympy import primerange, isprime, mobius, totient

OUT_DIR = os.path.dirname(os.path.abspath(__file__))


def squarefree_le(Q):
    """Return list of (q, mu(q), phi(q)) for squarefree q in [1, Q]."""
    if Q < 1:
        return []
    mu_arr = [1] * (Q + 1)
    sqf = [True] * (Q + 1)
    for p in primerange(2, Q + 1):
        for k in range(p, Q + 1, p):
            mu_arr[k] = -mu_arr[k]
        for k in range(p * p, Q + 1, p * p):
            sqf[k] = False
    phi_arr = list(range(Q + 1))
    for p in primerange(2, Q + 1):
        for k in range(p, Q + 1, p):
            phi_arr[k] -= phi_arr[k] // p
    out = []
    for q in range(1, Q + 1):
        if sqf[q]:
            out.append((q, mu_arr[q], phi_arr[q]))
    return out


def precompute_mu_phi(Q):
    """For ALL q in [1, Q], precompute mu and phi (used for q/gcd values)."""
    mu_arr = [1] * (Q + 1)
    sqf = [True] * (Q + 1)
    for p in primerange(2, Q + 1):
        for k in range(p, Q + 1, p):
            mu_arr[k] = -mu_arr[k]
        for k in range(p * p, Q + 1, p * p):
            sqf[k] = False
    for q in range(Q + 1):
        if not sqf[q]:
            mu_arr[q] = 0
    phi_arr = list(range(Q + 1))
    for p in primerange(2, Q + 1):
        for k in range(p, Q + 1, p):
            phi_arr[k] -= phi_arr[k] // p
    phi_arr[0] = 0
    return mu_arr, phi_arr


def compute_T_Q(N, Q, density):
    """T_Q(n) for n in [1, N] via direct Ramanujan-Fourier sum.
    Uses the closed form c_q(n) = mu(q/g) * phi(q) / phi(q/g) where g = gcd(n, q)
    (Hardy-Wright). Cost: O(N * #{sqf q <= Q})."""
    mu_arr, phi_arr = precompute_mu_phi(Q)
    sqf_data = [(q, mu_arr[q], phi_arr[q]) for q in range(1, Q + 1) if mu_arr[q] != 0]
    T = [0.0] * (N + 1)
    for q, mu_q, phi_q in sqf_data:
        if q == 1:
            for n in range(1, N + 1):
                T[n] += density * 1.0  # mu(1)/phi(1) * c_1(n) = 1 * 1 = 1
            continue
        # mu(q)/phi(q) * c_q(n) = mu(q/g) / phi(q/g) where g = gcd(n,q)
        # via the identity c_q(n) = mu(q/g) phi(q) / phi(q/g) and the
        # mu(q) cancellation.
        # Actually mu(q)/phi(q) * c_q(n) = mu(q) * mu(q/g) / phi(q/g)
        # (cancelling phi(q)). Since mu(q) = mu(g) * mu(q/g) for sqf q with
        # g = gcd(n,q) and g | q sqf so g sqf. Thus mu(q)*mu(q/g) = mu(g).
        # Therefore: mu(q)/phi(q) * c_q(n) = mu(g) / phi(q/g).
        for n in range(1, N + 1):
            g = math.gcd(n, q)
            qg = q // g
            mu_g = mu_arr[g]
            if mu_g == 0:
                continue
            T[n] += density * mu_g / phi_arr[qg]
    return T


def wheel_indicator(N, Q):
    """W_Q(n) = 1 if gcd(n, prod_{p<=Q} p) = 1, else 0."""
    primes_le_Q = list(primerange(2, Q + 1))
    W = [1] * (N + 1)
    W[0] = 0
    for p in primes_le_Q:
        for k in range(p, N + 1, p):
            W[k] = 0
    return W


def auc_score(score, label):
    """AUC of score vs binary label, ties handled with average rank."""
    paired = sorted(zip(score, label), key=lambda x: x[0])
    n_pos = sum(label)
    n_neg = len(label) - n_pos
    if n_pos == 0 or n_neg == 0:
        return float("nan")
    rank_sum_pos = 0.0
    i = 0
    rank_pos = 1
    n = len(paired)
    while i < n:
        j = i
        while j + 1 < n and paired[j + 1][0] == paired[i][0]:
            j += 1
        avg_rank = (rank_pos + (rank_pos + (j - i))) / 2.0
        for k in range(i, j + 1):
            if paired[k][1] == 1:
                rank_sum_pos += avg_rank
        rank_pos += (j - i + 1)
        i = j + 1
    U = rank_sum_pos - n_pos * (n_pos + 1) / 2.0
    return U / (n_pos * n_neg)


def precision_at_top_k(score, label, k):
    paired = sorted(zip(score, label), key=lambda x: -x[0])
    top = paired[:k]
    return sum(1 for _, l in top if l == 1) / k if k > 0 else float("nan")


def main():
    N = 30_000
    print(f"# E2.13 closure adversarial probe — S237")
    print(f"# T_Q (true Ramanujan-Fourier sum) vs W_Q (wheel sieve)")
    print(f"# N = {N}")
    print(f"# AUC(T_Q) vs AUC(W_Q): is T_Q a better primality discriminator?")

    primes_set = set(primerange(2, N + 1))
    label = [1 if n in primes_set else 0 for n in range(N + 1)]
    pi_N = len(primes_set)
    density = pi_N / N
    print(f"pi(N) = {pi_N}, density = {density:.6f}")

    results = []
    for Q in [6, 30, 210, 2310]:
        print(f"\n--- Q = {Q} ---")
        t0 = time.time()
        T = compute_T_Q(N, Q, density)
        t1 = time.time()
        W = wheel_indicator(N, Q)
        t2 = time.time()

        # Restrict to n in [Q+1, N] so the small-prime n's are NOT
        # boundary-special (n = p for small p satisfies p | n and is the
        # only prime divisible by p, sqewing wheel).
        # Actually the wheel sieve filters out all multiples of small
        # primes including the primes themselves (e.g., W_30(7) = 1 since
        # 7 is coprime to 2,3,5; but wait, 7 ≤ 30. Is 7 |  prod(p<=30)?
        # The wheel uses prod of primes <=Q, so 7 | prod, and W_30(7) = 0
        # since 7 % 7 = 0. That's the wheel filtering its own primes.).
        # T_30(7) however is non-zero because the truncated sum includes
        # contributions from q's not divisible by 7. So T_Q DOES retain
        # information that W_Q discards: the small primes themselves.
        # For a clean comparison, restrict to n > Q.
        n_start = Q + 1
        T_score = T[n_start:N + 1]
        W_score = [float(w) for w in W[n_start:N + 1]]
        labels = label[n_start:N + 1]
        n_test = len(labels)
        n_pos_test = sum(labels)

        # Also report unrestricted (n in [2, N])
        T_score_all = T[2:N + 1]
        W_score_all = [float(w) for w in W[2:N + 1]]
        labels_all = label[2:N + 1]

        auc_T = auc_score(T_score, labels)
        auc_W = auc_score(W_score, labels)
        auc_T_all = auc_score(T_score_all, labels_all)
        auc_W_all = auc_score(W_score_all, labels_all)

        # Check structural identity: is T_score perfectly correlated with W?
        # If T_Q = const * W on n > Q, then T's continuous score reduces
        # to W's binary indicator. Check unique levels of T.
        uniq_T_levels = sorted(set(round(t / density, 8) for t in T_score))
        n_uniq = len(uniq_T_levels)

        # Top-K precision at K = pi(N) (the "expected primes" count)
        prec_T_topK = precision_at_top_k(T_score, labels, n_pos_test)
        prec_W_topK = precision_at_top_k(W_score, labels, n_pos_test)

        T_primes = [t for t, l in zip(T_score, labels) if l == 1]
        T_comp = [t for t, l in zip(T_score, labels) if l == 0]

        out = {
            "Q": Q,
            "compute_time_T_sec": t1 - t0,
            "compute_time_W_sec": t2 - t1,
            "n_test_window_(Q+1,N)": n_test,
            "n_primes_in_window": n_pos_test,
            "auc_T_Q_window": auc_T,
            "auc_W_Q_window": auc_W,
            "auc_diff_window": auc_T - auc_W,
            "auc_T_Q_all_n_ge_2": auc_T_all,
            "auc_W_Q_all_n_ge_2": auc_W_all,
            "auc_diff_all": auc_T_all - auc_W_all,
            "T_unique_levels": n_uniq,
            "T_min": min(T_score),
            "T_max": max(T_score),
            "T_mean_on_primes": sum(T_primes) / len(T_primes) if T_primes else 0.0,
            "T_mean_on_composites": sum(T_comp) / len(T_comp) if T_comp else 0.0,
            "prec_topK_T": prec_T_topK,
            "prec_topK_W": prec_W_topK,
            "n_uniq_levels_sample": uniq_T_levels[:5] + (["..."] if n_uniq > 10 else []) + uniq_T_levels[-5:],
        }
        results.append(out)
        print(json.dumps(out, indent=2, default=str))

    with open(os.path.join(OUT_DIR, "results.json"), "w") as f:
        json.dump({"N": N, "results": results}, f, indent=2, default=str)
    print("\nWrote results.json")

    print("\n=== Verdict ===")
    diffs = [r["auc_diff_window"] for r in results]
    diffs_all = [r["auc_diff_all"] for r in results]
    max_abs_diff = max(abs(d) for d in diffs + diffs_all)
    print(f"max |AUC(T_Q) - AUC(W_Q)| over Q in [6, 2310] = {max_abs_diff:.6f}")
    if max_abs_diff < 0.005:
        print("  -> F-B holds: closure confirmed (T_Q ~ W_Q discrimination).")
    elif max_abs_diff < 0.05:
        print("  -> Marginal: re-verify with larger N or k-th moment exploitation.")
    else:
        print("  -> F-A holds: A-grade reopen needed.")


if __name__ == "__main__":
    main()
