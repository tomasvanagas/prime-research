"""Verify the conjectured S425.1 identity:

  q-reduce(D_omega1, q=1, H_N) = pi(N) * delta_1 + sum_{p^k, k>=2, p^k<=N} delta_{p^k}.

Print the full non-zero chip distribution at multiple N.
"""

import argparse
from sympy import isprime, factorint

from baker_norine_chi_p import q_reduce
from s425_inverse_chipfiring_probe import (
    hasse_cover_graph,
    divisor_omega1,
)


def is_prime_power(n):
    if n < 2:
        return False
    f = factorint(n)
    return len(f) == 1


def prime_powers_le(N, k_min=2):
    return [n for n in range(2, N + 1) if is_prime_power(n) and list(factorint(n).values())[0] >= k_min]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--N", nargs="+", type=int, default=[16, 32, 64, 128])
    args = ap.parse_args()

    for N in args.N:
        V, adj = hasse_cover_graph(N)
        D = divisor_omega1(N)
        D_red = q_reduce(D, q=1, V=V, adj=adj)
        nonzero = sorted([(v, c) for v, c in D_red.items() if c != 0])
        sink = D_red[1]
        pi_N = sum(1 for k in range(2, N + 1) if isprime(k))
        pp_ge2 = prime_powers_le(N, k_min=2)
        pp_chips = [(v, D_red[v]) for v in pp_ge2]
        other = sorted([(v, c) for v, c in D_red.items() if v != 1 and c != 0 and v not in pp_ge2])

        match_sink = (sink == pi_N)
        match_pp = all(D_red[v] == 1 for v in pp_ge2)
        match_other = all(c == 0 for v, c in D_red.items() if v != 1 and v not in pp_ge2)

        print(f"\n=== N = {N} ===")
        print(f"  pi(N) = {pi_N}")
        print(f"  prime powers k>=2: {pp_ge2}")
        print(f"  D'(sink=1) = {sink}; pi(N) match: {match_sink}")
        print(f"  Chips at prime powers k>=2: {pp_chips}; all == 1: {match_pp}")
        print(f"  Other non-zero chips off-sink: {other}; all zero: {match_other}")

        all_match = match_sink and match_pp and match_other
        print(f"  IDENTITY S425.1 holds at N={N}: {all_match}")


if __name__ == "__main__":
    main()
