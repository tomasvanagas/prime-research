"""
Full experiment for §D45 wild-swing: Baker-Norine graph Riemann-Roch / chip-
firing rank of the prime divisor.

Computes:
- Closed-form q-reduced identity for D_P on (a) divisibility graph and
  (b) Hasse cover graph, for N in {32, 64, 128, 256, 512}.
- Empirical comparison vs matched-degree random divisors.
- Generalised identity to other arithmetic divisors (sqfree, lambda, mu).
- Rank computation up to k=2 where tractable.
"""

import json
import time
import random
import statistics
from sympy import isprime, primepi, factorint, mobius

from baker_norine_chi_p import (
    divisibility_graph, q_reduce, is_winnable, rank,
    prime_divisor, random_matched_divisor, graph_stats, is_connected,
)
from baker_norine_hasse import hasse_graph


def liouville(n):
    if n == 1: return 1
    return (-1) ** sum(factorint(n).values())


def squarefree(n):
    if n == 1: return 1
    return 1 if max(factorint(n).values()) == 1 else 0


def signature(D, q, V, adj):
    Dr = q_reduce(D, q, V, adj)
    return {
        "D_red_at_q": Dr[q],
        "support": [v for v in V if v != q and Dr[v] > 0],
        "max_chip": max((Dr[v] for v in V if v != q), default=0),
        "total_nonq_chips": sum(Dr[v] for v in V if v != q),
    }


def run_identity_test():
    """Verify the closed-form q-reduced identities for various N."""
    print("=" * 70)
    print("PART 1: Closed-form q-reduced identity verification")
    print("=" * 70)
    out = {"divisibility_identity": [], "hasse_identity": []}
    for N in [16, 32, 64, 128, 256, 512]:
        # Divisibility graph identity: D' = (pi(N) - pi(N/2)) * delta_1 + sum_{p<=N/2} delta_p
        V, adj = divisibility_graph(N)
        DP = prime_divisor(V)
        Dr = q_reduce(DP, 1, V, adj)
        primes_low = [p for p in V if isprime(p) and p <= N // 2]
        primes_high = [p for p in V if isprime(p) and p > N // 2]
        pred_q = len(primes_high)
        pred_supp = set(primes_low)
        actual_supp = {v for v in V if v != 1 and Dr[v] > 0}
        match = (Dr[1] == pred_q) and (actual_supp == pred_supp) and \
                all(Dr[p] == 1 for p in pred_supp)
        out["divisibility_identity"].append({
            "N": N, "pi_N": int(primepi(N)), "pi_N_over_2": int(primepi(N // 2)),
            "D_red_at_q": Dr[1], "predicted_q": pred_q,
            "support_size": len(actual_supp), "predicted_support_size": len(pred_supp),
            "match": match,
        })
        print(f"  Γ_{N:4d}: pi(N)={int(primepi(N)):3d}, pi(N/2)={int(primepi(N//2)):3d}, "
              f"D'(q)={Dr[1]:3d} (pred={pred_q}), |supp|={len(actual_supp):3d} "
              f"(pred={len(pred_supp)}), match={match}")

        # Hasse identity: D' = pi(N) * delta_1
        V, adj = hasse_graph(N)
        DP = prime_divisor(V)
        Dr = q_reduce(DP, 1, V, adj)
        nonq_supp = [v for v in V if v != 1 and Dr[v] > 0]
        match_h = (Dr[1] == int(primepi(N))) and len(nonq_supp) == 0
        out["hasse_identity"].append({
            "N": N, "pi_N": int(primepi(N)),
            "D_red_at_q": Dr[1], "nonq_supp_size": len(nonq_supp),
            "match": match_h,
        })
        print(f"  H_{N:4d}: pi(N)={int(primepi(N)):3d}, "
              f"D'(q)={Dr[1]:3d}, |non-q supp|={len(nonq_supp)}, match={match_h}")
    return out


def run_random_comparison():
    """Compare D_P signature to random matched-degree divisors."""
    print()
    print("=" * 70)
    print("PART 2: Z-scores vs matched-density random divisors")
    print("=" * 70)
    out = []
    for graph_name, graph_fn in [("divisibility", divisibility_graph),
                                 ("hasse", hasse_graph)]:
        print(f"\n--- Graph: {graph_name} ---")
        for N in [32, 64, 128]:
            V, adj = graph_fn(N)
            n, m, g = graph_stats(V, adj)
            DP = prime_divisor(V)
            deg = sum(DP.values())
            sig_P = signature(DP, 1, V, adj)
            rng = random.Random(42)
            r_data = []
            for i in range(20):
                DR = random_matched_divisor(V, deg, rng)
                sig_R = signature(DR, 1, V, adj)
                r_data.append(sig_R)
            for key in ["D_red_at_q", "max_chip", "total_nonq_chips"]:
                vals = [d[key] for d in r_data]
                m_ = statistics.mean(vals)
                s_ = statistics.stdev(vals) if len(vals) > 1 else 0
                z = (sig_P[key] - m_) / s_ if s_ > 0 else (
                    0.0 if sig_P[key] == m_ else float('inf'))
                print(f"  N={N:4d} {graph_name:14s}: {key:20s}: "
                      f"D_P={sig_P[key]:5d}  R={m_:5.2f}±{s_:5.2f}  z={z:+.2f}")
            out.append({
                "graph": graph_name, "N": N, "n": n, "m": m, "g": g, "deg_DP": deg,
                "sig_DP": sig_P, "n_random": len(r_data),
                "rand_means": {k: statistics.mean([d[k] for d in r_data])
                               for k in ["D_red_at_q", "max_chip", "total_nonq_chips"]},
                "rand_sds": {k: (statistics.stdev([d[k] for d in r_data])
                                 if len(r_data) > 1 else 0)
                             for k in ["D_red_at_q", "max_chip", "total_nonq_chips"]},
            })
    return out


def run_other_divisors():
    """Test the Hasse identity generalisation for other divisors."""
    print()
    print("=" * 70)
    print("PART 3: Generalised Hasse identity for other arithmetic divisors")
    print("=" * 70)
    print()
    print("  Theorem: D'(1) = sum_{p prime, p<=N} D(p) on Hasse(N), q=1")
    print()
    out = []
    for N in [32, 64, 128]:
        V, adj = hasse_graph(N)
        print(f"\n--- N={N} ---")
        divisors = {
            "D_P": lambda v: 1 if isprime(v) else 0,
            "D_sqfree": lambda v: squarefree(v) if v != 1 else 0,
            "D_lam_pos": lambda v: 1 if v != 1 and liouville(v) == 1 else 0,
            "D_mu_pos": lambda v: 1 if v != 1 and mobius(v) == 1 else 0,
            "D_Omega2": lambda v: 1 if v != 1 and sum(factorint(v).values()) == 2 else 0,
        }
        for name, f in divisors.items():
            D = {v: f(v) for v in V}
            Dr = q_reduce(D, 1, V, adj)
            chips_on_primes = sum(D[v] for v in V if isprime(v))
            chips_on_composites = sum(D[v] for v in V if v != 1 and not isprime(v))
            ident_holds = (Dr[1] == chips_on_primes)
            comp_match = all(
                Dr[v] == D[v] for v in V if v != 1 and not isprime(v)
            )
            prime_zeroed = all(
                Dr[v] == 0 for v in V if v != 1 and isprime(v)
            )
            print(f"  {name:12s}: chips_on_primes={chips_on_primes:3d}, "
                  f"D'(q)={Dr[1]:3d}, comp_match={comp_match}, "
                  f"prime_zeroed={prime_zeroed}, identity={ident_holds and comp_match and prime_zeroed}")
            out.append({
                "N": N, "name": name, "chips_on_primes": chips_on_primes,
                "chips_on_composites": chips_on_composites,
                "D_red_at_q": Dr[1],
                "comp_match": comp_match, "prime_zeroed": prime_zeroed,
                "identity_holds": ident_holds and comp_match and prime_zeroed,
            })
    return out


def main():
    out = {}
    out["identity"] = run_identity_test()
    out["random_comparison"] = run_random_comparison()
    out["other_divisors"] = run_other_divisors()
    with open("full_results.json", "w") as f:
        json.dump(out, f, indent=2, default=str)
    print("\n=== Saved to full_results.json ===")


if __name__ == "__main__":
    main()
