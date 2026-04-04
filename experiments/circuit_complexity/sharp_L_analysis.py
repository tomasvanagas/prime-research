"""
Session 14: Critical analysis of the #L chain for pi(x).

CLAIMED CHAIN:
  BPSW correct → PRIMES ∈ TC^0 → PRIMES ∈ L → pi(x) ∈ #L → pi(x) ∈ NC^2

Each step is analyzed for correctness.

FINDING: The chain BREAKS at step "pi(x) ∈ #L".

The issue: #L counts accepting paths of NL machines with O(log N) work space,
where N is the input length. For pi(x), the input is x (N = log(x) bits).
The natural NL machine needs to handle candidate n (log(x) = N bits) as
work data, requiring O(N) work space — but #L only allows O(log N) = O(log log x).

This is the WORK SPACE MISMATCH: primality testing uses logspace relative to n,
but pi(x) counting needs n as intermediate data, which is N bits relative to x.
"""

print("""
Analysis of the #L chain for pi(x)
===================================

Step 1: BPSW correct → PRIMES ∈ TC^0
  VALID. Session 13 proved: BPSW test (MR base 2 + Strong Lucas)
  consists entirely of TC^0 operations (scalar pow + 2x2 MPOW + GCD).
  If BPSW has no pseudoprimes, this TC^0 circuit decides PRIMES.

Step 2: TC^0 ⊆ NC^1 ⊆ L
  VALID. Standard inclusions in circuit/space complexity.
  TC^0 → NC^1: replace MAJORITY gates with O(log n)-depth adder+compare.
  NC^1 → L: Borodin 1977, logspace can evaluate O(log n)-depth circuits.

Step 3: PRIMES ∈ L → pi(x) ∈ #L ???
  ANALYSIS: This step FAILS.

  #L = {f : f(x) = #accepting_paths(M, x) for some NL machine M}
  where M has O(log(|x|)) work space on input x.

  For pi(x): input is x, with |x| = N = log2(x) bits.
  The NL machine for pi(x) would need O(log N) = O(log log x) work space.

  Natural approach: nondeterministically guess n, check n ≤ x and PRIMES(n).
  Problem: n has N = log(x) bits, but work space is only O(log N) bits.
  Cannot store n on work tape.

  Alternative: generate n bit by bit via nondeterministic choices.
  Problem: after generating N bits of n, need to check primality.
  PRIMES(n) requires access to ALL bits of n simultaneously
  (e.g., compute 2^{n-1} mod n needs random access to bits of n-1).
  With O(log N) work space and no stored copy of n, this is impossible.

  The logspace algorithm for PRIMES(n) treats n as INPUT (on read-only tape).
  But in the NL machine for pi(x), n is NOT on the input tape — only x is.
  n must be stored as intermediate data, requiring O(N) work space.

  CORRECT CLASSIFICATION:
  - pi(x) ∈ #SPACE(N) with binary input: YES
    (NL machine with O(N) work space has 2^{O(N)} = poly(x) configurations)
    This gives circuits of size poly(x) = 2^{O(N)} — EXPONENTIAL in N.
  - pi(x) ∈ #SPACE(log N) = #L with binary input: UNKNOWN
    Would require O(log N) = O(log log x) work space — extremely tight.
    No known way to check primality of N-bit n in O(log N) work space
    when n is not on the input tape.

Step 4: #L ⊆ GapL ⊆ NC^2
  VALID (but moot since Step 3 fails).
  GapL → NC^2 via Berkowitz algorithm for determinant.

CONCLUSION:
  The chain breaks at Step 3.
  PRIMES ∈ L does NOT imply pi(x) ∈ #L (binary input model).
  The work space mismatch is fundamental: O(log N) ≠ O(N).

  "Is pi(x) ∈ #L?" is EQUIVALENT to "Is pi(x) ∈ GapL?" (from novel/gapl_question.md).
  Both ask for a poly(N)-size structure (matrix/graph) encoding pi(x).
  This remains OPEN and is the sharpest formulation of our target.

KEY INSIGHT:
  The natural connection between PRIMES ∈ L and pi(x) ∈ NC^2 DOES NOT EXIST.
  The problem of COMPUTING pi(x) is fundamentally harder than TESTING primality,
  even in the logspace/NC world. The gap comes from the SUMMATION of 2^N
  logspace computations, which cannot be done in logspace.

  Compare:
  - PRIMES(n): one evaluation, n on input tape, O(log log n) work space
  - pi(x) = SUM_{n=2}^{x} PRIMES(n): 2^N evaluations, each n NOT on tape
  - The SUM over 2^N terms requires generating each n as intermediate data
  - Storing n requires O(N) bits, not O(log N) bits

  This is analogous to the unary/binary distinction:
  - Unary input (size x): pi(x) ∈ FP (sublinear, O(x^{2/3}))
  - Binary input (size N = log x): pi(x) ∈ ? (no poly(N) algorithm known)
""")

# Verify the work space claims numerically
print("\nNumerical illustration of the work space mismatch:")
print("=" * 60)
for N in [10, 20, 50, 100, 333]:
    x = 2**N
    logN = N.bit_length()  # bits needed for log(N)
    print(f"N = {N} bits (x ≈ 10^{N*0.301:.0f}):")
    print(f"  Bits to store candidate n: {N}")
    print(f"  #L work space (O(log N)):  {logN}")
    print(f"  Ratio n/workspace: {N/logN:.1f}x")
    print(f"  #L configurations: O({N**2}) ≈ {N**2}")
    print(f"  #SPACE(N) configs: O(2^{N}) ≈ 10^{N*0.301:.0f}")
