#!/usr/bin/env python3
"""
Critical analysis: Is the STRONG Lucas test actually in TC^0?

The strong Lucas test has a loop: check V_{d*2^r} ≡ 0 for r = 0, ..., s-1
where n+1 = d * 2^s.

This loop has s iterations, and s can be up to O(log n) = O(N).
Each iteration squares V and subtracts 2*Q^k.

Is this O(N)-iteration loop compatible with TC^0 (constant depth)?

Author: Session 13 sub-agent
Date: 2026-04-04
"""

import math
from sympy import isprime, jacobi_symbol
import time


def analyze_strong_lucas_depth():
    """
    The strong Lucas test:
    1. Compute d, s where n+1 = d * 2^s  → TC^0 (bit operations)
    2. Compute U_d, V_d mod n            → 2x2 MPOW → TC^0
    3. Check U_d ≡ 0 (mod n)             → comparison → TC^0
    4. For r = 0 to s-1:
       Check V_{d*2^r} ≡ 0 (mod n)
       V_{d*2^{r+1}} = V_{d*2^r}^2 - 2*Q^{d*2^r}

    Step 4 has s iterations. Is this a problem?

    KEY INSIGHT: We don't actually need to iterate!
    V_{d*2^r} can be computed DIRECTLY as a 2x2 matrix power:
    V_{d*2^r} comes from M^{d*2^r} where M = [[P,-Q],[1,0]]

    So we need to compute M^{d*2^r} mod n for r = 0, 1, ..., s-1.
    These are s different matrix powers of the SAME 2x2 matrix.
    Each is independently computable in TC^0 (MPOW_2).

    But s can be O(N) = O(log n). So we need O(N) independent MPOW_2
    computations, each of which is TC^0.

    O(N) independent TC^0 computations in parallel → still TC^0!
    (TC^0 circuits have polynomial size, and O(N) copies of a poly(N) circuit
    is O(N * poly(N)) = poly(N).)

    The OR over s results is a single OR gate → TC^0.

    CONCLUSION: The strong Lucas test IS in TC^0.
    """
    print("=" * 70)
    print("STRONG LUCAS TEST: TC^0 DEPTH ANALYSIS")
    print("=" * 70)
    print()
    print("The concern: strong Lucas has a loop of s iterations (s up to O(N))")
    print()
    print("Resolution: Each V_{d*2^r} is an INDEPENDENT 2x2 MPOW computation.")
    print("We compute all s values in PARALLEL, each in TC^0.")
    print("Then take OR (is any of them 0?) → single threshold gate.")
    print()
    print("Detailed circuit structure for strong Lucas PRP(n):")
    print("  Layer 1: Extract d, s from n+1 = d * 2^s  [bit ops, TC^0]")
    print("  Layer 2: Compute Jacobi (D/n) = delta       [TC^0, via GCD]")
    print("  Layer 3: In parallel, for r = 0, ..., s-1:")
    print("           Compute M^{d*2^r} mod n             [MPOW_2, TC^0 each]")
    print("  Layer 4: Extract U_d (from r=0) and V_{d*2^r} from each")
    print("  Layer 5: Check U_d == 0 OR any V_{d*2^r} == 0 [OR gate, TC^0]")
    print()
    print("  Total: O(N) parallel MPOW_2 computations + O(1) extra layers")
    print("  Circuit size: O(N * poly(N)) = poly(N)")
    print("  Circuit depth: O(1)  [each MPOW_2 is constant depth]")
    print("  → STRONG LUCAS IS IN TC^0. ✓")
    print()

    # Verify: what is typical s for various n?
    print("Distribution of s (valuation of 2 in n+1):")
    s_values = []
    for n in range(3, 100000, 2):
        if not isprime(n):
            continue
        d = n + 1
        s = 0
        while d % 2 == 0:
            d //= 2
            s += 1
        s_values.append(s)

    max_s = max(s_values)
    avg_s = sum(s_values) / len(s_values)
    print(f"  For primes in [3, 100000]:")
    print(f"    Max s: {max_s}")
    print(f"    Avg s: {avg_s:.2f}")
    print(f"    s distribution: {dict(sorted([(s, s_values.count(s)) for s in set(s_values)]))}")
    print()

    # Even if s = N = log2(n), we just need N parallel MPOW_2 instances
    print(f"  Even for s = N = log2(n), we need N = {math.ceil(math.log2(100000))} parallel")
    print(f"  MPOW_2 instances. Each is poly(N) size. Total: O(N * poly(N)) = poly(N).")
    print()

    # Now: full BPSW in TC^0
    print("=" * 70)
    print("FULL BPSW CIRCUIT IN TC^0")
    print("=" * 70)
    print()
    print("BPSW = Miller-Rabin(base 2) + Strong Lucas(Selfridge)")
    print()
    print("Miller-Rabin(2):")
    print("  1. Write n-1 = d * 2^s                    [bit ops, TC^0]")
    print("  2. Compute 2^d mod n                       [scalar powering, TC^0]")
    print("  3. For r = 0, ..., s-1:")
    print("     Compute 2^{d*2^r} mod n                 [s independent scalar pows, TC^0]")
    print("  4. Check conditions                        [comparisons + OR, TC^0]")
    print("  Total: O(N) parallel scalar powering + O(1) extra → TC^0")
    print()
    print("Strong Lucas(Selfridge):")
    print("  1. Find D: first in {5,-7,9,-11,...} with (D/n) = -1")
    print("     This requires testing O(log n) Jacobi symbols in worst case")
    print("     Each Jacobi is TC^0; test all in parallel, take first = -1")
    print("     → TC^0 (poly(N) parallel Jacobi computations + MUX)")
    print("  2. P = 1, Q = (1-D)/4")
    print("  3. n+1 = d * 2^s")
    print("  4. Compute U_d, V_d, V_{d*2^r} as above → TC^0")
    print("  5. Check conditions → TC^0")
    print("  Total: TC^0")
    print()
    print("BPSW output = MR(2) passes AND Strong Lucas passes → AND gate")
    print("  → BPSW IS IN TC^0.")
    print()

    # The critical assumption
    print("=" * 70)
    print("THE GAP: CORRECTNESS")
    print("=" * 70)
    print()
    print("BPSW is in TC^0 as a COMPUTATION. The open question:")
    print("Does BPSW correctly identify ALL primes and ALL composites?")
    print()
    print("Known status:")
    print("  - Verified correct for all n < 2^64 (exhaustive computation)")
    print("  - No BPSW pseudoprime is known")
    print("  - Pomerance heuristic: BPSW PSPs probably exist (very rare)")
    print("  - No PROOF that BPSW is correct for all n")
    print()
    print("Conditional results:")
    print("  - GRH → Miller's test (O(log^2 n) bases) is correct → PRIMES in TC^0")
    print("  - BPSW correct → PRIMES in TC^0 (immediate, with O(1) TC^0 tests)")
    print()
    print("The weakest sufficient condition for PRIMES in TC^0:")
    print("  'There exist O(poly(N)) TC^0-computable witness values")
    print("   such that checking the witness is TC^0 and the combined")
    print("   test is always correct.'")
    print()
    print("  Under GRH: witnesses = MR bases {2, 3, ..., O(log^2 n)}")
    print("    → O(log^2 n) = O(N^2) scalar powerings, each TC^0")
    print("    → total circuit size O(N^2 * poly(N)) = poly(N) → TC^0")
    print()
    print("  Under BPSW: witnesses = {MR base 2, Lucas (Selfridge D)}")
    print("    → O(1) tests, each TC^0 → TC^0")
    print()
    print("  Unconditionally: UNKNOWN whether any O(poly(N)) witness set exists")
    print("    that makes a combined TC^0 test always correct.")

    # Final status table
    print()
    print("=" * 70)
    print("SUMMARY TABLE")
    print("=" * 70)
    print()
    print("  Test                      | In TC^0? | Deterministic? | Condition")
    print("  " + "-" * 72)
    print("  MR(2)                     | YES      | NO             | -")
    print("  MR(2,3,...,37)            | YES      | YES for n<3e24 | unconditional")
    print("  MR(2,...,2ln^2(n))        | YES      | YES for all n  | GRH")
    print("  Strong Lucas (Selfridge)  | YES      | NO             | -")
    print("  BPSW                      | YES      | YES for n<2^64 | unconditional")
    print("  BPSW                      | YES      | YES for all n  | BPSW conjecture")
    print("  QFT (Grantham)            | YES      | NO             | -")
    print("  AKS                       | UNKNOWN  | YES            | unconditional")
    print("  Trial division            | NO       | YES            | unconditional")
    print("  Wilson's                  | NO       | YES            | unconditional")


if __name__ == "__main__":
    t0 = time.time()
    analyze_strong_lucas_depth()
    print(f"\nTotal time: {time.time()-t0:.1f}s")
