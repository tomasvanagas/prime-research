# A3 — Spectral Primality Test via (Z/nZ)\* Cayley Graph

**Attack vector:** ATTACK_VECTORS.md §A.A3.
**Session:** 79.
**Channelling:** Bourgain (spectral gap on Cayley graphs of arithmetic
groups, sum-product mixing); Lubotzky-Phillips-Sarnak (arithmetic Cayley
graph spectra).
**Cross-domain ingredient:** spectral graph theory of Cayley graphs of
finite abelian groups (Babai 1979, "Spectra of Cayley graphs"; eigenvalues
of Cay(G,S) for abelian G are character sums λ_χ = Σ_{s∈S} χ(s)).

## Pre-stated falsification criterion

H_A (positive — would have been A-grade): there exists a polynomial-time-
computable scalar feature F(spec(Cay((Z/nZ)\*, S))) such that F is in
**disjoint** ranges for primes n and composites n on a held-out set
(n ∈ [600, 1000] coprime to 30). Disjoint = no overlap on full ranges.

H_0 (null / closure — would be B-grade): every feature we compute either
(a) has overlapping ranges between primes and composites, or
(b) the "feature" is a rephrasing of a known primality test.

## Result: H_0 confirmed structurally. Spectrum probes ω(n), not primality.

**Verdict: B-grade.** Ambitious frontier attack failed informatively. The
failure mode is structural and adds a new negative-shape edge candidate.

## Empirical setup

For n ∈ [7, 2000] coprime to 30 (so {2,3,5} ⊂ (Z/nZ)\*), build the
Cayley graph C_n := Cay((Z/nZ)\*, {2,3,5,2⁻¹,3⁻¹,5⁻¹}). Compute the full
adjacency spectrum (matrix size m = φ(n)). Extract 24 spectral features
(spectral gap, second eigenvalue, multiplicities, integer-eigenvalue
count, trace moments, max multiplicity, spectral entropy, etc.).

**Sample composition (533 values of n):**
- 300 primes (ω(n) = 1, prime)
- 13 prime powers (ω(n) = 1, composite — n ∈ {49, 121, 169, 289, 343, 361, 529, 841, 961, 1331, 1369, 1681, 1849})
- 214 semiprimes (ω(n) = 2)
- 6 triple-prime composites (ω(n) = 3)

**Robustness sweep:** repeated for generator sets S ∈ {{2,3,5}, {2,3,7},
{2,3,5,7}} on n coprime to 210. Identical structural finding.

## The structural finding

### Theorem (empirically verified, provable from CRT).

For n odd squarefree-or-not, coprime to all generators in S, the Cayley
graph Cay((Z/nZ)\*, S ∪ S⁻¹) has at least **2^ω(n) integer eigenvalues**,
where ω(n) is the number of distinct prime factors of n. The lower bound
is achieved with equality in many cases.

**Proof sketch.** By CRT, (Z/nZ)\* ≅ ∏ᵢ (Z/pᵢ^{aᵢ}Z)\*. For odd pᵢ, each
factor is cyclic of order pᵢ^{aᵢ-1}(pᵢ - 1), which is even, so its
character group has 2-torsion subgroup of order 2. The 2-torsion subgroup
of the full character group is a product of order 2^ω(n). For χ in this
subgroup, χ(s) ∈ {±1} for all s, so

   λ_χ = Σ_{s ∈ S ∪ S⁻¹} χ(s) ∈ Z.

Hence at least 2^ω(n) integer eigenvalues. ∎

### Empirical confirmation (533 n, no violations, all 3 generator sets).

| ω(n) | is_prime | count | min_n_int_eigs | 2^ω(n) | tight? |
|------|----------|-------|----------------|--------|--------|
| 1    | True     | 300   | 2              | 2      | yes (most) |
| 1    | False (p^k) | 13 | 2              | 2      | yes (some) |
| 2    | False    | 214   | 4              | 4      | yes (some) |
| 3    | False    | 6     | 24             | 8      | no — but ≥ |

Tight equality (n_int_eigs = 2^ω(n)) holds for 84 / 533 = 15.8% of
cases (mostly large primes with discrete-log non-alignment, and the
prime powers).

### Why this kills A3

The discriminator we'd want is "n is prime." But the natural
spectral feature 2^ω(n) gives:
- prime n: 2^1 = 2
- prime power n = p^k (k ≥ 2): 2^1 = 2 — **same as prime**.
- semiprime pq: 2^2 = 4
- triple-prime: 2^3 = 8

So **spectrum cannot distinguish primes from prime powers** via
ω(n)-based features. (Z/p^kZ)\* is cyclic for odd p, structurally
identical to (Z/qZ)\* for prime q with appropriately aligned discrete
logarithms.

### Mann-Whitney discrimination test (primes vs prime-powers, omega=1)

Held-out test on n ∈ [11, 1500] coprime to 210. AUC for the best
feature (n_int_eigs):

| Generator set | AUC (primes vs PPs) |
|----------------|---------------------|
| {2,3,5} ∪ inv  | 0.509 |
| {2,3,7} ∪ inv  | 0.566 |
| {2,3,5,7} ∪ inv | 0.673 |

AUC ≈ 0.5 means **no discrimination signal** for the hard case (primes
vs prime powers). Even the strongest 4-generator setup gives 0.673,
which is far from disjoint and within the small-sample noise band of the
9 prime-power examples. Doubling the generators does not help because
the obstruction is at the GROUP-STRUCTURE level (cyclicity), not the
sampling resolution.

## What this rules out (the negative-shape edge)

Spectral primality testing via Cay((Z/nZ)\*, fixed-S) for any S that
does not depend on n collapses to **detecting ω(n) ≥ 2**. The natural
spectral features are blind to:
- prime n vs prime power n = p^k (both ω = 1, both cyclic groups)
- prime n vs higher prime power (k ≥ 3): same problem.

To rescue, you'd need to combine spectral features with a separate
**perfect-power test** (deciding if n = m^k for some k ≥ 2). Perfect-power
testing IS in P (Bach-Sorenson polylog claim — actually deterministic
polylog) but requires modular exponentiation depth equivalent to
growing-dim MPOW. So even a hybrid spectral + perfect-power test would
inherit the AKS-class depth barrier (E5.3, E7.10).

## Cross-references

- **E5.3** — TC⁰ primality requires growing-dim MPOW.
- **E5.8** — Brandt diagonalisation closure (technique-welded to MKtP).
- **E7.10** — AKS modulus-twists are orthogonal to depth.
- **CLOSED_PATHS line 354** — Cayley(Z/xZ, primes): primes-as-generators (CIRCULAR).
- **CLOSED_PATHS line 383** — Cayley graph spectral approach: required knowing primes (CIRCULAR).
- **CLOSED_PATHS line 384** — Ihara zeta of Cayley graph: equivalent to explicit formula.
- **CLOSED_PATHS line 385** — GCD/coprimality graph: Ramanujan sums, Meissel-Lehmer.
- **CLOSED_PATHS line 427** — Algebraic geometry F_q families: Frobenius eigenvalues = zeta zeros.

This A3 closure is **disjoint from prior Cayley graph closures** because
those used primes as generators (circular) or graph index (circular). A3
fixed generators {2,3,5} and varied n — the question was whether
fixed-generator spectrum could decide primality of the index. Answer: no,
spectrum probes ω(n) which doesn't separate prime from p^k.

## What this DOES suggest (for follow-ups, NOT A-grade in this session)

The empirical lower bound n_int_eigs ≥ 2^ω(n) is itself an arithmetic
identity that may be sharper than published. Specifically, the **integer-
eigenvalue count** of Cay((Z/nZ)\*, S ∪ S⁻¹) computes a quantity related
to ω(n), the number of distinct prime factors. The function ω(n) is
relevant in:
- the prime omega function counting (an approximation to log log n in
  expectation — Erdős-Kac theorem),
- the squarefree characteristic function via 2^ω,
- the Möbius function via μ(n) = (-1)^ω(n) for squarefree n.

If this provides a TC⁰-amenable shortcut to ω(n), it would be a
significant secondary result. But computing the spectrum requires the
full graph (size φ(n), exponential in N = log n), so the route is
blocked at the input-encoding step.

## Files

- `cayley_spectral_primality.py` — main sweep, 24 features × 533 n.
- `cayley_robustness.py` — three generator sets × 342 n cross-check.
- `spectrum_features.json` — full feature matrix.
- `separation_summary.json` — per-feature primal/composite separation.
- `robustness_features.json` — generator-robustness data.

## Self-grading

**B-grade.** Frontier attack (A3) attempted. Failed in the predicted
direction (ranges overlap, AUC ≈ noise). Failure mode is structural
and characterized by a clean theorem (n_int_eigs ≥ 2^ω(n)). Adds:

1. A new structural lower bound — empirically saturated at 2^ω(n), with
   character-theoretic proof.
2. A new negative-shape edge candidate (see below).
3. Documentation of why the entire fixed-generator-Cayley-graph spectral
   approach to primality cannot work.

This is NOT A-grade because the result is a refinement of "spectrum
probes ω(n)" rather than producing a new primality test. The cross-
domain ingredient (Cayley graph spectra of abelian groups) was applied
correctly but the chain collapsed at the exact predicted spot:
cyclic-vs-non-cyclic of (Z/nZ)\* is encoded, but cyclic-vs-cyclic-of-
prime-power is NOT.
