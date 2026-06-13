"""
mobius_bisection.py — Möbius-side bisection of pi(x) and the 4-cell
{squarefree × Omega-parity} decomposition.

Composes:
    E1.6  (parity bisection pi(x) = A(x) XOR C_3(x) mod 2)
    E2.2  (Liouville bisection pi(x) = (x - L(x))/2 - C_3(x))
    E2.10 (parity trap L(x) mod 2 = x mod 2)

with new content:

    Möbius bisection (B):
        pi(x) = (Q(x) - M(x))/2 - C_3*(x)
    where
        Q(x) = #{n <= x : squarefree}
        M(x) = sum_{n <= x} mu(n)         (Mertens function)
        C_3*(x) = #{n <= x : sqfree, composite, omega(n) odd, omega(n) >= 3}

    Bridge identity (closed-form for non-squarefree Omega-odd count):
        NS_o(x) := #{n <= x : not sqfree, Omega(n) odd}
                = (x - Q(x) - L(x) + M(x)) / 2.

    4-cell decomposition (Omega-parity x sqfree):
        x = S_e + S_o + NS_e + NS_o,
        S_e = (Q + M)/2,  S_o = (Q - M)/2,
        NS_e = ((x - Q) + (L - M))/2,  NS_o = ((x - Q) - (L - M))/2.

Pre-stated falsifiers (must all PASS for any x <= N for the construction
to count as integer-exact):
    F1: pi(x) = (Q(x) - M(x))/2 - C_3*(x)   for all x in [1, N].
    F2: NS_o(x) = (x - Q - L + M)/2          for all x in [1, N].
    F3: C_3(x) - C_3*(x) = NS_o(x)           for all x in [1, N].
    F4: A(x) XOR B(x) = NS_o(x) mod 2        (where A=(x-L)/2, B=(Q-M)/2).
    F5: pi(x) mod 2 = B(x) mod 2 XOR C_3*(x) mod 2 for all x.

If any falsifier fires, the identities are wrong (or the implementation is).

CLI:
    python3 mobius_bisection.py --N 1000000

Outputs JSON to mobius_bisection_results.json + summary printed to stdout.
"""

from __future__ import annotations

import argparse
import json
import math
import time
from collections import Counter
from pathlib import Path


def sieve_omega_Omega_squarefree(N: int):
    """Return arrays omega[1..N], Omega[1..N], is_sqfree[1..N], lpf[1..N]."""
    omega = [0] * (N + 1)
    Omega = [0] * (N + 1)
    is_sqfree = [True] * (N + 1)
    is_sqfree[0] = False
    # smallest prime factor sieve to also get distinct prime exponents
    spf = [0] * (N + 1)
    for i in range(2, N + 1):
        if spf[i] == 0:
            for j in range(i, N + 1, i):
                if spf[j] == 0:
                    spf[j] = i
    # factor each n via spf
    for n in range(2, N + 1):
        m = n
        last_p = 0
        while m > 1:
            p = spf[m]
            e = 0
            while m % p == 0:
                m //= p
                e += 1
            if p != last_p:
                omega[n] += 1
                last_p = p
            Omega[n] += e
            if e >= 2:
                is_sqfree[n] = False
    is_sqfree[1] = True  # convention: 1 is squarefree
    return omega, Omega, is_sqfree


def compute_indicators(omega, Omega, is_sqfree, N: int):
    """Per-n indicators and cumulative summatories."""
    # mu(n) and lambda(n)
    mu = [0] * (N + 1)
    lam = [0] * (N + 1)
    mu[1] = 1
    lam[1] = 1
    for n in range(2, N + 1):
        lam[n] = 1 if Omega[n] % 2 == 0 else -1
        if is_sqfree[n]:
            mu[n] = 1 if omega[n] % 2 == 0 else -1
        else:
            mu[n] = 0
    return mu, lam


def cumsum(seq):
    out = [0] * len(seq)
    s = 0
    for i, v in enumerate(seq):
        s += v
        out[i] = s
    return out


def main(N: int) -> dict:
    t0 = time.time()
    omega, Omega, is_sqfree = sieve_omega_Omega_squarefree(N)
    mu, lam = compute_indicators(omega, Omega, is_sqfree, N)

    # Per-n indicators of the 4-cell decomposition.
    # Note: 1 is squarefree, omega(1)=0 (even). So 1 is in S_e.
    is_prime = [False] * (N + 1)
    is_S_e = [False] * (N + 1)
    is_S_o = [False] * (N + 1)
    is_NS_e = [False] * (N + 1)
    is_NS_o = [False] * (N + 1)
    is_C3 = [False] * (N + 1)         # composite with Omega odd, Omega >= 3
    is_C3_star = [False] * (N + 1)    # sqfree composite, omega odd, omega >= 3

    for n in range(1, N + 1):
        if is_sqfree[n]:
            if omega[n] % 2 == 0:
                is_S_e[n] = True
            else:
                is_S_o[n] = True
                # composite if omega >= 3 (omega = 1 means n = prime; for sqfree)
                if omega[n] >= 3:
                    is_C3_star[n] = True
        else:
            if Omega[n] % 2 == 0:
                is_NS_e[n] = True
            else:
                is_NS_o[n] = True
        # primality: omega == 1 AND Omega == 1 (then n is prime; only sqfree)
        if omega[n] == 1 and Omega[n] == 1:
            is_prime[n] = True
        # C_3 (Liouville-side): composite with Omega odd, Omega >= 3
        if (not is_prime[n]) and n != 1 and Omega[n] % 2 == 1 and Omega[n] >= 3:
            is_C3[n] = True

    # Cumulative summatories.
    pi_cum = cumsum([1 if is_prime[n] else 0 for n in range(N + 1)])
    L_cum = cumsum(lam)
    M_cum = cumsum(mu)
    Q_cum = cumsum([1 if is_sqfree[n] else 0 for n in range(N + 1)])
    C3_cum = cumsum([1 if is_C3[n] else 0 for n in range(N + 1)])
    C3star_cum = cumsum([1 if is_C3_star[n] else 0 for n in range(N + 1)])
    NSo_cum = cumsum([1 if is_NS_o[n] else 0 for n in range(N + 1)])
    NSe_cum = cumsum([1 if is_NS_e[n] else 0 for n in range(N + 1)])
    Se_cum = cumsum([1 if is_S_e[n] else 0 for n in range(N + 1)])
    So_cum = cumsum([1 if is_S_o[n] else 0 for n in range(N + 1)])

    # ============= Falsifier checks =============
    failures: dict[str, int] = {}

    def check_all(name: str, predicate):
        first_fail = -1
        for x in range(1, N + 1):
            if not predicate(x):
                first_fail = x
                break
        if first_fail >= 0:
            failures[name] = first_fail

    # F1: pi(x) == (Q(x) - M(x))/2 - C3*(x).
    check_all(
        "F1_mobius_bisection",
        lambda x: pi_cum[x] == (Q_cum[x] - M_cum[x]) // 2 - C3star_cum[x]
        and (Q_cum[x] - M_cum[x]) % 2 == 0,
    )

    # F2: NS_o(x) == (x - Q - L + M)/2.
    check_all(
        "F2_NSo_closed_form",
        lambda x: NSo_cum[x] == (x - Q_cum[x] - L_cum[x] + M_cum[x]) // 2
        and (x - Q_cum[x] - L_cum[x] + M_cum[x]) % 2 == 0,
    )

    # F3: C_3(x) - C_3*(x) == NS_o(x).
    check_all(
        "F3_C3_diff_equals_NSo",
        lambda x: C3_cum[x] - C3star_cum[x] == NSo_cum[x],
    )

    # F4: A XOR B = NS_o mod 2 (parity).
    # A(x) = (x - L(x))/2,  B(x) = (Q(x) - M(x))/2
    def F4_pred(x):
        A = (x - L_cum[x]) // 2
        B = (Q_cum[x] - M_cum[x]) // 2
        return (A ^ B) & 1 == NSo_cum[x] & 1

    check_all("F4_AB_xor_NSo", F4_pred)

    # F5: pi(x) mod 2 = B mod 2 XOR C_3* mod 2.
    def F5_pred(x):
        B = (Q_cum[x] - M_cum[x]) // 2
        return (pi_cum[x] & 1) == ((B & 1) ^ (C3star_cum[x] & 1))

    check_all("F5_pi_parity_via_mobius", F5_pred)

    # F6: original Liouville bisection (sanity, E2.2).
    def F6_pred(x):
        A = (x - L_cum[x]) // 2
        return pi_cum[x] == A - C3_cum[x]

    check_all("F6_liouville_bisection", F6_pred)

    # ============= Pseudorandomness measures (parity of B, C3*, A, C3) =============
    # Sequences over x in [1, N]
    A_par = [(((x - L_cum[x]) // 2) & 1) for x in range(1, N + 1)]
    B_par = [(((Q_cum[x] - M_cum[x]) // 2) & 1) for x in range(1, N + 1)]
    C3_par = [(C3_cum[x] & 1) for x in range(1, N + 1)]
    C3s_par = [(C3star_cum[x] & 1) for x in range(1, N + 1)]
    pi_par = [(pi_cum[x] & 1) for x in range(1, N + 1)]
    NSo_par = [(NSo_cum[x] & 1) for x in range(1, N + 1)]

    def H(p1: float) -> float:
        if p1 <= 0 or p1 >= 1:
            return 0.0
        return -(p1 * math.log2(p1) + (1 - p1) * math.log2(1 - p1))

    def freq1(seq) -> float:
        return sum(seq) / len(seq)

    def lag1_ac(seq) -> float:
        n = len(seq)
        s = sum(seq)
        ssq = sum(v * v for v in seq)
        if ssq * n - s * s == 0:
            return 0.0
        cross = sum(seq[i] * seq[i + 1] for i in range(n - 1))
        # Pearson correlation between seq[:-1] and seq[1:]
        m1 = sum(seq[:-1]) / (n - 1)
        m2 = sum(seq[1:]) / (n - 1)
        num = cross / (n - 1) - m1 * m2
        var1 = sum((v - m1) ** 2 for v in seq[:-1]) / (n - 1)
        var2 = sum((v - m2) ** 2 for v in seq[1:]) / (n - 1)
        if var1 <= 0 or var2 <= 0:
            return 0.0
        return num / math.sqrt(var1 * var2)

    def joint_MI_bits(s1, s2) -> float:
        # Both binary. MI in bits.
        n = len(s1)
        c = Counter(zip(s1, s2))
        m_a = Counter(s1)
        m_b = Counter(s2)
        mi = 0.0
        for (a, b), cnt in c.items():
            p_ab = cnt / n
            p_a = m_a[a] / n
            p_b = m_b[b] / n
            if p_ab > 0:
                mi += p_ab * math.log2(p_ab / (p_a * p_b))
        return mi

    stats = {
        "freq1_A_par": freq1(A_par),
        "freq1_B_par": freq1(B_par),
        "freq1_C3_par": freq1(C3_par),
        "freq1_C3s_par": freq1(C3s_par),
        "freq1_pi_par": freq1(pi_par),
        "freq1_NSo_par": freq1(NSo_par),
        "lag1_ac_A_par": lag1_ac(A_par),
        "lag1_ac_B_par": lag1_ac(B_par),
        "lag1_ac_C3_par": lag1_ac(C3_par),
        "lag1_ac_C3s_par": lag1_ac(C3s_par),
        "lag1_ac_NSo_par": lag1_ac(NSo_par),
        "lag1_ac_pi_par": lag1_ac(pi_par),
        "MI_A_B": joint_MI_bits(A_par, B_par),
        "MI_C3_C3s": joint_MI_bits(C3_par, C3s_par),
        "MI_A_C3": joint_MI_bits(A_par, C3_par),
        "MI_B_C3s": joint_MI_bits(B_par, C3s_par),
        "MI_pi_NSo": joint_MI_bits(pi_par, NSo_par),
    }

    # Final aggregate sizes
    aggregate = {
        "N": N,
        "pi_N": pi_cum[N],
        "L_N": L_cum[N],
        "M_N": M_cum[N],
        "Q_N": Q_cum[N],
        "C3_N": C3_cum[N],
        "C3s_N": C3star_cum[N],
        "NSo_N": NSo_cum[N],
        "NSe_N": NSe_cum[N],
        "Se_N": Se_cum[N],
        "So_N": So_cum[N],
        "wall_seconds": time.time() - t0,
    }

    out = {
        "aggregate": aggregate,
        "failures": failures,  # empty dict means all 6 falsifiers PASS
        "stats": stats,
    }
    return out


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--N", type=int, default=200000,
                        help="Upper bound for x (and for sieve).")
    parser.add_argument("--out", type=str, default=None,
                        help="JSON output path. Defaults to "
                             "<script_dir>/mobius_bisection_results.json")
    args = parser.parse_args()

    res = main(args.N)
    print(json.dumps(res, indent=2))

    out_path = (
        Path(args.out)
        if args.out
        else Path(__file__).parent / "mobius_bisection_results.json"
    )
    out_path.write_text(json.dumps(res, indent=2))
    print(f"\n[wrote {out_path}]")
