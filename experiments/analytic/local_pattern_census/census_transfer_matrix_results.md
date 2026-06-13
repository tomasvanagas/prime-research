# Counting A_k without enumeration: the transfer-matrix state space is Θ(A_k) (S519)

## Question (open item 4, PROGRAM.md NEXT ACTION)

Count the admissible-pattern census number

> A_k = #{ S ⊆ allowed(k) : for every prime q ≤ k with q ∤ k, the occupied
> residues {o mod q : o ∈ S} do not cover all of Z_q }

(allowed(k) = offsets in [0,k) divisible by no prime dividing k) **without
enumerating all A_k patterns**, so as to push the entropy law
A_k = e^{(1+o(1))·π(k)} from its DFS reach (k ≤ 80, A₈₀ ≈ 5.8·10⁶) toward
k = 128 (A₁₂₈ ≈ e²⁵ ≈ 10¹¹). The DFS in `local_pattern_census.py` visits
~A_k nodes, so it is enumeration-bound. The NEXT ACTION asked for a
transfer-matrix / DP over positions carrying, per small prime, the occupied
residues — *or* a proof that the transfer-matrix state space is itself
exponential (closing the method).

**Result: the transfer-matrix state space is Θ(A_k) — exponential with the
SAME entropy exponent as A_k itself. The DP does not escape the
enumeration barrier (it is in fact strictly worse than DFS: Θ(A_k) memory
vs Θ(k)). The method is closed.** En route we obtained two clean
structural facts that DO scale: an exact **active-prime reduction**
(A_k depends only on the primes ≤ B(k), with B(k) ≈ ρ*(k) ≪ k) and the
identity **W(B(k)) = ρ*(k) = the classical maximal admissible-tuple size**,
which explains why ln A_k scales as π(k).

Code: `census_transfer_matrix.py` (one script, `--selftest`, `--k`,
`--scan KMIN KMAX [--powers]`, `--reduce K`). 9 selftest groups; the TM
reproduces every known A_k bit-for-bit.

## Method

**Transfer-matrix DP (method "tm").** Process the allowed offsets in order;
the state is the *coverage tuple* — for each enforced prime q, the bitmask
of residues mod q occupied so far (a Python int, q-bit fields concatenated).
Each offset is included/excluded; on inclusion the masks OR in the offset's
residues, and a state where *any* enforced prime becomes fully covered is
dropped immediately (coverage only grows — a covered prime can never
recover, so the pattern is permanently inadmissible). The DP collapses all
subsets sharing a coverage tuple into one `{state: count}` entry. A_k is the
sum of the counts of surviving states (each misses a class in every enforced
prime). Counts are exact Python big-ints — no floating-point cancellation.

**Active-prime reduction (rigorous, two independent proofs).** Enforcing a
prime is a *monotone filter* on subsets — adding a prime can only keep or
decrease the count. Two ways to find the minimal sufficient prime set B(k):

- *count stabilisation* (`minimal_B`): enforce primes {q ≤ B}; if the count
  equals the count for {q ≤ next prime}, no pattern was removed by the primes
  in between, so they never bind and the count is exact;
- *weight bound* (`reduce_primes`, no count needed): W(B) = |allowed| −
  minkill(B), where minkill(B) is the minimum (over choosing one residue
  class per prime in B to leave empty) of the size of the union of those
  classes — computed exactly by branch-and-bound (`minkill_bb`, never
  iterating the Π q drop-tuples). Any S that misses a class mod every prime
  in B has weight ≤ W(B); if W(B) < the next enforceable prime p, then S has
  < p elements and so cannot cover Z_p (covering needs ≥ p elements) — hence
  S misses a class mod *every* enforceable prime, i.e. S is admissible.
  Therefore A_k = #{S misses a class mod all q ∈ B}, and the larger primes
  are provably irrelevant.

The two agree (selftest 8): the W-based B(k) equals the count-based B(k) at
k = 8,16,32,64. The W-based one needs no count, so it gives B(128), B(256)
even though A₁₂₈, A₂₅₆ are unreachable.

## Measured — the state space (clean family k = 2^m)

```
   k          A_k        states  A_k/states   lnA    k/lnk   ratio
   8           13             7      1.857    2.565   3.847  0.6667
  16          106            64      1.656    4.663   5.771  0.8081
  32         3573          1505      2.374    8.181   9.233  0.8861
  64      1581920       1020333      1.550   14.274  15.389  0.9276
```

- **States track A_k to a bounded constant factor (1.55–2.37, not growing).**
  Slope of ln(states) vs k/ln k = **1.0202**; slope of ln A_k vs k/ln k =
  **1.0105**. The two exponents coincide ⇒ **#states = Θ(A_k)**. The
  many-to-one collapse of the DP is only a ~2× constant, because the large
  enforced primes (11, 13 at k = 64) give each admissible pattern a nearly
  distinct coverage tuple. The transfer matrix is therefore exponential and
  cannot reach k = 128 (it would need ~5·10¹⁰ live states ≈ A₁₂₈/2). It is
  strictly worse than the DFS, which uses Θ(k) memory for the same Θ(A_k)
  time. **The transfer-matrix / position-DP route is closed.**
- The entropy-law trend is reconfirmed independently of the DFS: slope
  **1.01 → conjecture 1.0**, ratio climbing 0.667 → 0.928 (geometric
  convergence). No new clean-family exact point is reachable (k = 128 is the
  next power of two and is past both DFS and TM).

## Measured — the active-prime reduction B(k) and the weight identity

```
  k     |allowed|  enforceable  B(k) = primes ≤   W(B)=ρ*(k)   #primes ignored
  8        4         {3}            3                 3            0
 16        8         {3,5,7}        5                 5            1
 32       16         {3,5,7,11,13}  7                 9            2
 64       32         10 primes      13               16           5
128       64         17 primes      23               28           9
256      128         30 primes      47               52          16
```

Two clean facts:

1. **B(k) ≪ (enforceable primes), growing like ρ*(k) ~ k/ln k.** At k = 128,
   A₁₂₈ depends on only the **8 primes {3,5,7,11,13,17,19,23}**; the 9 larger
   enforceable primes {29,…,61} provably never bind. (This shrinks the
   problem but does not break the exponential state wall — see above.)

2. **W(B(k)) = ρ*(k) = the maximal admissible-tuple size in a window of k.**
   The reduction's weight bound is *tight*: a weight-W(B) subset missing a
   class mod every prime in B is itself admissible (it cannot cover any
   larger prime), so W(B(k)) equals the maximum weight of any admissible
   pattern, ρ*(k). The measured ρ*(k) = 3,5,9,16,28,52 at k = 8,16,32,64,
   128,256:
   - **k ≤ 32 reproduce the values already cross-validated** in
     `local_pattern_census_results.md` (max realised/admissible weights
     3,5,9, matching narrowest-tuple H(9)=30 ≤ 31 < H(10)=32, H() = minimal
     admissible-tuple diameter, OEIS A008407).
   - **k = 64,128 extend the cross-check exactly** against H(): a window of
     k integers admits an m-tuple iff H(m) ≤ k−1. ρ*(64)=16 because
     H(16)=60 ≤ 63 < H(17)=66; ρ*(128)=28 because H(28)=126 ≤ 127 < H(29)=
     132. (k=256: ρ*=52, consistent with the H()-growth, not separately
     table-checked here.) The growth ρ*(k) ≈ (1+o(1))·k/ln k
     (ρ*/(k/ln k) = 0.78,0.87,0.97,1.04,1.06,1.13 over k = 8…256) is the
     Hensley–Richards densest-admissible-tuple law.

   This is the structural reason for the entropy law: admissible patterns are
   exactly the subsets that fit inside the ρ*(k)-bounded admissible geometry,
   and ρ*(k) = Θ(k/ln k) = Θ(π(k)), so ln A_k = Θ(π(k)). (The lower bound
   A_k ≥ 2^{ρ*(k)} gives ratio ≥ 0.69·ρ*/(k/ln k); the measured ratio
   exceeds it, since A_k counts more than the subsets of a single maximal
   tuple — the route to the *exact* constant 1 is a generating-function /
   cluster-expansion asymptotic analysis, a separate open avenue from exact
   counting.)

## What this closes / what stays open

- **Closed:** the transfer-matrix / position-DP method for *exact* A_k at
  k = 128. Its state space is Θ(A_k), measured (slope 1.02 ≈ A_k's 1.01;
  collapse factor a bounded ~2×). No offset ordering helps — the final live
  state count is order-independent and equals ~A_k/2, so memory is Θ(A_k).
- **Delivered (scales):** the exact active-prime reduction (A_k uses only
  primes ≤ B(k) ~ ρ*(k)) and the W(B(k)) = ρ*(k) identity, both rigorous and
  validated; B(128), B(256) computed though A₁₂₈, A₂₅₆ are uncountable here.
- **Still open:** the entropy *constant* → 1 (needs an analytic /
  generating-function argument, not enumeration); an exact A₁₂₈ (would need a
  genuinely sub-A_k algorithm — none found; the natural transfer matrix is
  exponential, and inclusion–exclusion over the coverage constraints has the
  catastrophic cancellation / term-count blow-up noted in PROGRAM open
  item 4).

## What would falsify

- Any A_k from the TM disagreeing with DFS / the known exact value (k ≤ 80) —
  a bug (selftests 2,3 assert bit-for-bit equality, brute == DFS == TM).
- The stabilisation check passing (count(B) == count(next)) while a larger
  enforced set later changes the count — impossible for a monotone filter,
  so a sign of a state-collision bug (selftest 4 asserts monotonicity +
  stabilisation; selftest 8 asserts W-based B == count-based B).
- The collapse factor A_k/states growing with k (would reopen the method) —
  not observed (1.55–2.37, flat); ln(states) and ln A_k slopes within 0.01.
- minkill_bb returning a non-minimal union (would make W(B) wrong and could
  drop a binding prime) — selftest 9 checks it against a full brute over
  drop-tuples; selftest 8 checks W(B) ≥ the brute-enumerated max admissible
  weight.

## Files

`census_transfer_matrix.py`, this file. Cross-refs:
`local_pattern_census.py` / `_results.md` (the DFS counter and the
A_k=13,106,3573,… data this validates against and explains),
`experiments/constructions/automaton_width_dichotomy/` (the width profile
whose cliff height is A_k+1), Hardy–Littlewood Conjecture B,
Hensley–Richards densest-admissible-tuple bound, narrowest-tuple H(k) tables.
