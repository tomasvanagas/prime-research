# Derandomization and Prime Counting: A Rigorous Analysis

**Date:** 2026-04-04 (Session 18)
**Goal:** Analyze whether derandomization theory provides new algorithmic
approaches to computing pi(x), or proves impossibility.

---

## 1. Impagliazzo-Wigderson (1997): BPP = P Under Hardness

### The Theorem
**IW97:** If there exists a language L in E = DTIME(2^{O(n)}) that requires
circuits of size 2^{Omega(n)}, then BPP = P.

More precisely: if E requires exponential-size circuits, then every BPP language
has a deterministic polynomial-time algorithm. The proof constructs a
pseudorandom generator G: {0,1}^{O(log n)} -> {0,1}^n that fools polynomial-size
circuits.

### Application to pi(x)

**Question:** Can computing pi(x) be cast as a BPP problem?

**Analysis:** Consider the decision version: "Is pi(x) >= k?" where x and k are
given in binary (N = log x bits each). This is a well-defined language.

**Key observation:** We do NOT know that pi(x) is in BPP (binary input model).
In fact, all known algorithms for pi(x) run in time 2^{Omega(N^{1/2})}, which is
super-polynomial in N. For pi(x) to be in BPP, we would need a randomized
polynomial-time (in N = log x) algorithm -- this is EXACTLY the open question.

**What IW97 gives us:** If we could show pi(x) is in BPP (under some hardness
assumption or unconditionally), then IW97 + assumed circuit hardness would give
pi(x) in P. But we CANNOT show pi(x) is in BPP without already having essentially
solved the problem.

**Converse direction:** Could IW97 help REFUTE pi(x) in BPP? No. IW97 is a
positive result (derandomization), not a separation result. It says BPP collapses
to P under hardness assumptions, but this doesn't tell us anything about whether
pi(x) is in BPP or P in the first place.

**Conditional argument:** If BPSW is unconditionally correct, then PRIMES is in
TC^0 (Session 13). Then pi(x) = #TC^0 query. The question "#TC^0 in BPP?" is
subsumed by "#TC^0 in NC?" -- our existing open question. Since TC^0 ⊂ P, we
have #TC^0 ⊂ #P, and computing #P functions in BPP would imply NP ⊂ BPP ⊂ P/poly
(Adleman), which doesn't help directly.

**Verdict: IW97 does NOT provide a new route.** The bottleneck is not randomness
removal -- it is that no polynomial-time (randomized or deterministic) algorithm
is known in the binary input model. IW97 presupposes the existence of an efficient
randomized algorithm, which we lack.

---

## 2. Nisan-Wigderson Generators

### The Framework
**NW94 (Nisan-Wigderson 1994):** Constructs a PRG G: {0,1}^{O(log^2 n)} -> {0,1}^n
that fools logarithmic-space computations (i.e., is epsilon-biased against L/poly).
Under stronger hardness assumptions, the seed length drops to O(log n) and fools
polynomial-size circuits.

**Nisan (1992):** For bounded-space computation, Nisan's PRG fools space-S
computations with seed length O(S * log n). For S = O(log n), seed = O(log^2 n).

### Application: PRG Based on Primes?

**Idea:** Could a pseudorandom generator whose output distribution "encodes" prime
positions allow reconstruction of pi(x)?

**Analysis:** This conflates two things:
1. PRGs as derandomization tools (replacing random bits in BPP algorithms)
2. PRGs as structured distributions that encode combinatorial information

For (1): As established above, we don't have a BPP algorithm to derandomize.

For (2): Consider constructing a distribution D on {0,1}^n such that
sum_i D(i) * 1_prime(i) = pi(n)/2^n. If D could be generated from a short
seed, then summing over all 2^{O(log n)} seeds would give pi(x).

**The problem:** For this to work, D must distinguish primes from composites.
But a PRG that fools bounded computations explicitly CANNOT distinguish things
that bounded computations cannot. If pi(x) is not computable in the bounded
class, the PRG is useless for reconstructing pi(x).

**More precisely:** If G fools circuits of size s, then for any circuit C of
size s: |Pr[C(G(U_k))=1] - Pr[C(U_n)=1]| < epsilon. This means G cannot help
compute functions that require circuits larger than s.

**A different angle -- primes as pseudorandom source:** The Fermat residue
function f(n) = 2^{n-1} mod n is pseudorandom-looking for consecutive n
(Session 16: zero autocorrelation at all lags). If this function were a good PRG,
it might fool the "count primes" computation. But:
- f(n) = 1 iff n passes MR(2), which is NOT equivalent to primality
- The correlation between f(n) and 1_prime(n) is high but NOT exact (pseudoprimes exist)
- Even if it were a perfect PRG, "counting" through a PRG requires evaluating
  the discriminator on exponentially many inputs

**Verdict: NW generators do NOT provide a new route.** PRGs are tools for
removing randomness from efficient algorithms; they cannot create efficiency
where none exists. The prime indicator function has cryptographically hard
structure that prevents any short-description generation.

---

## 3. Hardness of Approximation: Smooth Part to Exact

### The Question
R(x) approximates pi(x) to within O(sqrt(x)/ln(x)) (assuming RH). This is
~50% of the digits. Is there a derandomization argument that "amplifies"
approximation to exactness?

### Relevant Results

**Sipser (1983), Lautemann (1983):** BPP ⊂ Sigma_2^P ∩ Pi_2^P. A BPP algorithm
with error 1/3 can be amplified to error 2^{-poly(n)} by majority vote. This is
about error PROBABILITY, not approximation quality.

**Approximate counting (Stockmeyer 1983):** Given an NP oracle, one can
APPROXIMATE #SAT to within a factor (1 + epsilon) in polynomial time. Combined
with Toda's theorem (#P ⊂ P^{PP}), this is part of the counting hierarchy.

**Our situation:** R(x) gives pi(x) ± O(sqrt(x)). To get pi(x) exactly, we need
to determine the correction delta(x) with |delta(x)| = O(sqrt(x)).

**Is there a "rounding" trick?** In some cases, if you know a quantity is an
integer and you have an approximation within 1/2, you can round. Here, the
approximation error is O(sqrt(x)), not O(1). We'd need to improve the approximation
from sqrt(x) error to <1 error.

**Key negative result (Session 17):** The communication matrix rank of the
oscillatory residual pi(x) - R(x) is EXACTLY 2^{N/2-1}. This means:
- The residual has FULL information-theoretic complexity
- No polynomial number of "corrections" to R(x) can determine it
- Specifically, D(pi(x) - R(x)) = Omega(N/2) communication bits

**Goldreich-Wigderson (1999) on search-to-decision:** For NP-complete problems,
search reduces to decision. But pi(x) is not NP-complete (it's in P for unary
input), so this reduction framework doesn't apply.

**Self-reducibility:** If f is downward self-reducible (computing f(n) reduces
to f(m) for m < n), then approximate counting can yield exact counting via
binary search on the error. Is pi(x) self-reducible?

pi(x) = pi(x-1) + 1_prime(x). But this trivial recursion requires x steps.
The Meissel-Lehmer recursion pi(x) = Phi(x, pi(sqrt(x))) - pi(sqrt(x)) + 1 is
a self-reduction but doesn't reduce the complexity class (the "small" instances
still have exponential binary input).

**Verdict: No derandomization argument can amplify R(x) to exact pi(x).**
The approximation gap is O(sqrt(x)), which carries 2^{N/2-1} bits of
information (proven by communication complexity, Session 17). Closing this gap
requires computing those bits, which IS the hard part of the problem.

---

## 4. Kabanets-Impagliazzo (2004): PIT and Circuit Lower Bounds

### The Theorem
**KI04:** If polynomial identity testing (PIT) requires super-polynomial
circuits (i.e., PIT is not in P/poly), then either:
(a) NEXP is not in P/poly, or
(b) the permanent does not have polynomial-size arithmetic circuits.

Equivalently: if NEXP ⊂ P/poly and VNP = VP, then PIT ∈ P.

### Connection to pi(x)

**Session 17 result:** The multilinear polynomial pi_N(x_1,...,x_N) representing
pi(x) as a function of the bits of x has determinantal complexity dc(pi_N) >=
2^{N/2-1} + 2. This means pi_N cannot be expressed as det(M) for any polynomial-
size matrix M with affine entries in the x_i variables.

**Does KI04 help?** The connection is:
1. KI04 connects PIT derandomization to arithmetic circuit lower bounds
2. We HAVE arithmetic complexity lower bounds for pi_N (exponential dc)
3. Does this yield PIT derandomization? Or vice versa?

**Analysis:** KI04 works in the algebraic (arithmetic circuit) model. The relevant
complexity classes are VP (polynomial-size arithmetic circuits) and VNP (permanent-
like). The key connection is:

- If the permanent has polynomial arithmetic circuits (VP = VNP), then PIT ∈ P
  (conditional on NEXP ⊂ P/poly).
- If pi_N has exponential dc, this means pi_N is NOT in VBP (algebraic branching
  programs). Since VBP ⊂ VP, pi_N might or might not be in VP.

**Critical issue:** pi_N is defined over the integers/rationals, but the
exponential dc lower bound means it's not in VBP over ANY field. However:
- pi_N is NOT a family of polynomials in the usual algebraic complexity sense
  (it's one polynomial per N, not a uniform family)
- pi_N has degree N and N variables -- it's "tiny" by algebraic complexity standards
  (the hard part is that N itself grows as log x)
- KI04's implications are about uniform families of polynomial-size circuits

**The real gap:** KI04 says that if we can prove lower bounds for EXPLICIT
polynomial families (like the permanent), we get PIT derandomization. We can prove
lower bounds for pi_N. But pi_N is not "explicit" in the right sense -- it's not
a VNP-family (we can't even define it algebraically without reference to primes).

Moreover, even if KI04 gave us PIT ∈ P, this wouldn't help compute pi(x).
PIT solves "is this polynomial identically zero?" -- a fundamentally different
question from evaluating pi_N at a specific point.

**A more nuanced attempt:** Could we express "pi(x) = k" as a polynomial identity
and use PIT? We'd need a polynomial P(x,k) that vanishes iff pi(x) = k. Such a
polynomial would need to encode primality testing for all numbers up to x, which
is circular.

**Verdict: KI04 is tangentially related but does NOT help.** The exponential dc
of pi_N lives in a different regime than the polynomial families KI04 addresses.
KI04 connects PIT to VNP vs VP, while our lower bound is about a specific
non-algebraic function. The frameworks don't interact productively.

---

## 5. Circuit-to-Formula Conversion and Lower Bounds

### Relevant Results

**Valiant (1976):** Any formula (circuit where every gate has fan-out 1) for a
function f: {0,1}^n -> {0,1} computed by circuits of size s can be converted to
a formula of size O(s^2). Conversely, a formula lower bound of Omega(n^2/log n)
is known for an explicit function (Andreev, Hastad-Boppana).

**Nechiporuk (1966):** For certain explicit functions, formula size >=
Omega(n^2/log^2 n). This remains the best known formula lower bound for
EXPLICIT functions.

**Karchmer-Wigderson (1990):** Formula depth (= formula size^{1/2+o(1)} via
balanced formulas) equals the communication complexity of the KW game associated
with f. Specifically, D(KW_f) = formula depth of f.

### Application to pi(x)

**Our result (Session 17):** rank(pi_N) = 2^{N/2-1} + 2 for the balanced
partition. This gives:

1. **Formula size lower bound:**
   The log-rank of the communication matrix gives a lower bound on communication
   complexity: D(pi_N) >= log2(rank) = N/2 - O(1).

   For the KW relation associated with pi(x) (Alice gets an input where pi(x) >= k,
   Bob gets an input where pi(x) < k, they must find a bit position where the
   inputs differ), D(KW_pi) >= N/2 - O(1).

   By KW theorem: formula depth >= N/2 - O(1), hence formula size >= 2^{N/2 - O(1)}.

   **This matches our communication complexity lower bound exactly.**

2. **Comparison with Valiant's conversion:**
   If pi(x) had circuits of size s, it would have formulas of size O(s^2).
   Our formula lower bound 2^{N/2} implies circuit size >= 2^{N/4} (via Valiant).
   This is WEAKER than the direct communication complexity bound of 2^{N/2} on
   TC^0 circuit size (Session 16).

3. **Comparison with Nechiporuk:**
   Nechiporuk's method partitions variables and uses the number of distinct
   subfunctions. For pi_N with balanced partition, we have 2^{N/2-1} + 2 distinct
   subfunctions (= rank). This gives formula size >= sum_blocks (log distinct
   subfunctions)^2 / N. With a single balanced block:
   >= (N/2)^2 / N = N/4.
   With N/2 blocks of size 2 each: sum >= (N/2) * c for small c.
   Nechiporuk gives at most O(N^2/log N) which is polynomial -- far weaker than
   our exponential bound.

   **Key distinction:** Nechiporuk's bound is polynomial because it works for ALL
   explicit functions. Our exponential bound is specific to pi(x) and comes from
   the exponential RANK of the communication matrix, which Nechiporuk's method
   cannot capture directly.

4. **Relationship between the bounds:**
   - Communication complexity D(pi_N) = Theta(N/2) bits (balanced partition)
   - Formula depth >= N/2 - O(1) (via KW)
   - Formula size >= 2^{N/2 - O(1)} (exponential)
   - Circuit size >= 2^{N/4 - O(1)} (via Valiant)
   - TC^0 circuit size >= 2^{N/2}/poly(N) (via communication complexity, Session 16)

   The TC^0 lower bound is the STRONGEST, but applies only to constant-depth circuits.
   The formula lower bound is exponential for unbounded-depth formulas.
   The general circuit lower bound (2^{N/4}) is the weakest.

**Verdict:** The formula complexity picture is CONSISTENT with but does not improve
upon our existing barriers. The KW connection gives a clean reformulation: pi(x)
formula complexity >= 2^{N/2}. But this does not yield new algorithmic insight --
it merely confirms that the sqrt(x) barrier manifests in ALL computational models.

---

## 6. Natural Proofs Barrier (Razborov-Rudich 1997)

### The Theorem
**RR97:** A "natural proof" against a circuit class C is a property P of Boolean
functions such that:
1. **Constructivity:** P can be decided in time 2^{O(n)} (polynomial in the truth table)
2. **Largeness:** P is satisfied by >= 2^{-O(n)} fraction of all functions
3. **Usefulness:** No function in C satisfies P

If secure pseudorandom function generators exist (a standard cryptographic
assumption), then no natural proof can prove super-polynomial lower bounds against
TC^0 (or any class containing such generators).

### Application to pi(x) and TC^0

**The question:** Can we prove pi(x) is not computable by polynomial-size TC^0
circuits? If not, the natural proofs barrier explains why.

**Analysis:**

1. **Our Session 17 Fourier analysis results are "natural":**
   We showed the prime indicator has near-random Fourier spectrum, noise sensitivity,
   and total influence. These are ALL "natural" properties in the Razborov-Rudich
   sense: they are constructive (computable from truth table in 2^{O(N)} time),
   large (satisfied by ~50% of functions), and we hoped they would be useful
   (distinguishing pi_N from TC^0-computable functions).

   **But RR97 says:** If one-way functions exist, such properties CANNOT
   distinguish pi_N from TC^0 functions. A TC^0-computable pseudorandom function
   would have the SAME Fourier profile as a random function, so any "natural"
   measure that shows pi_N looks random also shows the PRF looks random.

2. **The communication complexity bound IS potentially non-natural:**
   Our rank(pi_N) = 2^{N/2-1} + 2 result is a SPECIFIC structural property of pi_N,
   not a random-like property. The rank is NOT large -- a random function has rank
   2^{N/2} with high probability, and pi_N has rank 2^{N/2-1} + 2, which is HALF
   of random. This suggests the bound might evade the natural proofs barrier.

   However: the rank is still EXPONENTIAL, which is the same order as random.
   The factor-of-2 savings is a structural consequence (even numbers are not prime)
   but doesn't change the asymptotics. A TC^0-computable PRF would also likely have
   rank close to 2^{N/2}, so our rank bound cannot distinguish pi_N from a PRF.

   **Therefore:** Our communication complexity lower bound, while correct, is
   consistent with pi(x) being in TC^0. It shows that the MULTILINEAR POLYNOMIAL
   representation requires exponential resources, but TC^0 circuits are NOT
   constrained to work through multilinear polynomials.

3. **What RR97 means for our project:**

   **Scenario A: One-way functions exist (standard assumption).**
   Then no natural proof can show pi(x) requires super-polynomial TC^0 circuits.
   Any proof that pi(x) is not in TC^0 must use non-natural methods, such as:
   - Diagonalization (but this can't prove circuit lower bounds -- Baker-Gill-Solovay)
   - Algebraization (but Aaronson-Wigderson 2009 shows this also has barriers)
   - Indirect/structural arguments specific to pi(x)

   This means: **we almost certainly CANNOT prove pi(x) is not in TC^0 using
   currently known techniques**, assuming standard cryptographic assumptions.

   **Scenario B: One-way functions do NOT exist (unlikely but possible).**
   Then natural proofs could work, and our Fourier/communication results might
   suffice. But if OWFs don't exist, then P = BPP unconditionally (Impagliazzo-
   Wigderson), and the entire landscape changes.

   **Scenario C: pi(x) IS in TC^0 (our goal scenario).**
   Then natural proofs are irrelevant -- we're trying to prove membership, not
   separation. The barrier only blocks impossibility proofs.

4. **Critical implication for the project:**
   The natural proofs barrier means we should focus on UPPER bounds (finding
   algorithms), not lower bounds (proving impossibility). We cannot hope to
   prove pi(x) is not in TC^0 or NC using standard methods. But we CAN hope
   to find an algorithm -- nothing blocks that direction.

   **This is GOOD NEWS for our project:** The absence of a proof that pi(x)
   is not in NC is not because we haven't tried hard enough -- it's because
   such a proof would require techniques beyond current complexity theory.
   The problem remains genuinely open in both directions.

**Verdict: The natural proofs barrier prevents proving pi(x) is not in TC^0
using known techniques.** This is important meta-information: we should not
spend effort trying to prove impossibility (which is blocked by RR97), but
rather search for algorithms (which is not blocked). However, RR97 does NOT
say an algorithm exists -- it merely says we can't prove one doesn't.

---

## 7. Synthesis: What Derandomization Theory Tells Us

### Summary Table

| Direction | Relevant? | New Route? | Why/Why Not |
|-----------|-----------|------------|-------------|
| IW97 (BPP=P) | No | No | No BPP algorithm known for pi(x) |
| NW generators | No | No | PRGs can't create efficiency from nothing |
| Approx -> Exact | No | No | Gap is sqrt(x) = 2^{N/2-1} bits of info |
| KI04 (PIT) | Tangential | No | pi_N not in right algebraic framework |
| Formula bounds | Confirmatory | No | Matches existing sqrt(x) barrier |
| Natural proofs | YES | METARESULT | Blocks impossibility proofs |

### The Key Insight from This Analysis

**Derandomization theory assumes the existence of efficient algorithms and
asks whether randomness helps.** Our problem is more fundamental: we don't
know whether ANY efficient algorithm exists (randomized or deterministic).

The derandomization paradigm is:
  HARD FUNCTION EXISTS → PRG EXISTS → BPP = P

For our problem, we need:
  ??? → EFFICIENT pi(x) ALGORITHM EXISTS

Derandomization can't fill in the "???". The chain goes in the wrong direction.

### What WOULD Help (Speculative)

The one genuinely interesting connection is through **hardness amplification**
in reverse. Normally:
- Worst-case hard function → average-case hard function → PRG

But what if we could run this in reverse for pi(x)?
- pi(x) - R(x) LOOKS random (Session 17 Fourier analysis)
- Could we prove that if the residual is pseudorandom, then pi(x) is efficiently
  computable? (i.e., if the hard part looks random, maybe it can be replaced by
  a PRG?)

**This fails because:** The residual pi(x) - R(x) is not pseudorandom -- it IS
random-looking but is a DETERMINISTIC function. A PRG that generates something
with the same distribution would give a random variable, not the exact value.
For exact computation, we need the exact residual, not a distributional match.

### Connection to #TC^0 ⊆ NC?

The most relevant derandomization question is actually:

**Can the "randomness" in BPSW's primality decisions across consecutive integers
be derandomized into a structured batch computation?**

BPSW performs 2^{n-1} mod n for each n, and these residues look cryptographically
random (Session 16). If we could show that the function n -> (2^{n-1} mod n)
has small circuit complexity as a batch computation, we could potentially count
primes efficiently.

But Session 16 showed: the Fermat residue coupling (exponent and modulus both
depend on n) makes this function have zero autocorrelation at all lags. This is
precisely the signature of a function that CANNOT be derandomized -- it behaves
like a PRF keyed by n.

### Broader Complexity-Theoretic Picture

The situation for pi(x) mirrors the broader landscape:
- We cannot prove super-polynomial circuit lower bounds for explicit functions
  (natural proofs barrier + relativization + algebraization)
- We cannot prove BPP != P unconditionally
- We cannot prove pi(x) is or isn't in NC/TC^0
- The best we can do is conditional results and structural understanding

**Our contribution:** The exact communication complexity formula rank(pi_N) =
2^{N/2-1} + 2 is a CONCRETE structural result about pi(x) that sits alongside
(but does not resolve) these meta-questions. It shows that any efficient algorithm
must work in a fundamentally different representation than multilinear polynomials
over the bits of x.

---

## 8. New Closed Paths from This Analysis

| Path | Failure Mode | Session |
|------|-------------|---------|
| IW97 BPP=P derandomization for pi(x) | No BPP algorithm exists to derandomize | 18 |
| NW/Nisan PRG-based prime reconstruction | PRGs can't create computational efficiency | 18 |
| Approximation amplification R(x) -> pi(x) via derandomization | Info barrier: 2^{N/2-1} bits in residual | 18 |
| KI04 PIT connection via dc(pi_N) | Wrong algebraic framework; pi_N not explicit VNP family | 18 |
| Reverse hardness amplification (random residual -> PRG) | Residual is deterministic, not distributional | 18 |
| Batch Fermat residue derandomization | Zero autocorrelation = PRF-like; cannot derandomize | 18 |

## 9. Impact on Open Directions

**No new algorithmic direction found.** All six derandomization connections
analyzed fail to provide a route to efficient pi(x) computation.

**One important metaresult confirmed:** The natural proofs barrier (RR97)
explains WHY we cannot prove pi(x) is not in TC^0. This is not a failure of our
analysis -- it is a fundamental limitation of known proof techniques. This
reinforces that the problem is genuinely open: neither algorithms nor impossibility
proofs are within reach of current methods.

**Updated assessment:** The derandomization lens confirms the existing barrier
structure but does not change it. The sqrt(x) wall remains universal. The only
paths that survive are:
1. #TC^0 ⊆ NC? (pure complexity theory)
2. Novel number-theoretic identity (serendipity)
3. BPSW correctness proof (number theory)
4. Entirely new computational paradigm (unknown unknowns)

---

## References

- Impagliazzo, R. & Wigderson, A. (1997). "P = BPP if E requires exponential circuits." STOC 1997.
- Nisan, N. & Wigderson, A. (1994). "Hardness vs randomness." JCSS 49(2), 149-167.
- Nisan, N. (1992). "Pseudorandom generators for space-bounded computations." Combinatorica 12(4), 449-461.
- Kabanets, V. & Impagliazzo, R. (2004). "Derandomizing polynomial identity tests means proving circuit lower bounds." Computational Complexity 13, 1-46.
- Razborov, A.A. & Rudich, S. (1997). "Natural proofs." JCSS 55(1), 24-35.
- Valiant, L.G. (1976). "Graph-theoretic arguments in low-level complexity." MFCS 1976.
- Karchmer, M. & Wigderson, A. (1990). "Monotone circuits for connectivity require super-logarithmic depth." SIAM J. Disc. Math. 3(2), 255-265.
- Nechiporuk, E.I. (1966). "On a Boolean function." Soviet Mathematics Doklady 7, 999-1000.
- Goldreich, O. & Wigderson, A. (1999). "Tiny families of functions with random properties." Lecture Notes in Computer Science.
- Stockmeyer, L. (1983). "The complexity of approximate counting." STOC 1983.
- Sipser, M. (1983). "A complexity-theoretic approach to randomness." STOC 1983.
- Aaronson, S. & Wigderson, A. (2009). "Algebrization: A new barrier in complexity theory." TOCT 1(1).
- Allender, E. (1999). "The permanent requires large uniform threshold circuits." Chicago J. Theor. Comp. Sci.
- Mereghetti, C. & Palano, B. (2000). "Threshold circuits for iterated matrix product and powering." RAIRO 34(1), 39-46.
- Hesse, W., Allender, E. & Barrington, D.A.M. (2002). "Uniform constant-depth threshold circuits for division and iterated multiplication." JCSS 65(4), 695-716.
- Session 16-17 results: communication complexity rank(pi_N) = 2^{N/2-1} + 2, Fourier analysis of prime indicator.

---

Generated: 2026-04-04 (Session 18)
