# Session-by-Session Key Insights

Detailed findings from each research session. Read the relevant session
if you need deep context on a specific topic. For a quick overview,
see the Status section of CLAUDE.md.

---

## Session 12
- "Is pi(x) in NC?" is EQUIVALENT to our target.
  All known approaches produce circuits of size 2^{Theta(N)} (exponential in input).

## Session 13
(a) BPSW IS computable in TC^0 (MR=scalar pow, Strong Lucas=2x2 MPOW, Jacobi=GCD).
    PRIMES in TC^0 iff BPSW (or similar) is unconditionally correct.
    Verified correct to 2^64. GRH also suffices (Miller's test = O(N^2) scalar pows).
(b) Prime indicator ANF degree = Theta(N) over GF(2), 50% sparsity (random-like).
    No GF(2) algebraic shortcut for counting primes.
(c) Zeta zero minimum K_min ~ 0.35 * x^{0.27}: power law, no reordering helps.
(d) Spectral graph approaches all circular or equivalent to Meissel-Lehmer.
(e) No new algorithmic breakthroughs in 2025-2026 literature.

## Session 14
(a) **PRIMES in L and pi(x) in NC are INDEPENDENT questions.** The chain
    PRIMES in L -> pi(x) in #L -> pi(x) in NC^2 BREAKS due to workspace mismatch:
    NL machine needs O(N) bits for candidate n, but #L allows only O(log N).
(b) **I-E fractional parts carry O(2^k) independent bits** -> no determinant
    smaller than 2^{Theta(sqrt(x)/log(x))} can encode the Legendre sieve.
    A GapL algorithm MUST avoid floor functions entirely.
(c) **Lucy DP matrices have NO algebraic structure**: unipotent, displacement
    rank 50-60% of dimension, full-rank product. No compression possible.
(d) **Nonlinear sieve breaks parity in theory but NOT in efficiency**: nonlinear
    ops on floor values CAN distinguish primes from semiprimes but cost >= O(x^{2/3}).
(e) **All algebraic variety approaches fail**: low-dim can't encode pi(x),
    high-dim has slow point counting, Frobenius eigenvalues = zeta zeros.
(f) **pi(x) mod m is NOT a linear recurrence** for any m (Mauduit-Rivat consequence).

## Session 15
(a) **Determinantal complexity connection**: pi(x) as degree-N multilinear polynomial
    in bits. Found NxN det representations for N=2,3,4. "Is pi(x) in GapL?" <=>
    "polynomial determinantal complexity?" For N>=10, GENERIC polynomials don't fit.
(b) **#TC^0 subset NC? is THE complexity question**: If BPSW in TC^0 (conditional),
    pi(x) in NC iff #TC^0 subset NC. Fermat residue coupling prevents batch counting.
(c) **Uniformity is the true barrier**: Nonuniform circuits trivially poly(N).
    Hard part = UNIFORM construction. Natural proofs barrier blocks lower bounds.
(d) **ALL randomized approaches fail**: zeta zero sampling (100% needed), probabilistic
    sieve (10^6x worse), hash counting, quantum counting -- all rigorously excluded.
(e) **Divide-and-conquer fails**: error accumulates O(sqrt(x)) regardless of depth.
(f) **Monotone complexity inapplicable**: [pi(x)>=k] = [x>=p(k)] trivially O(N).
(g) **Arithmetic circuits don't help**: VP=VNC^2 (depth free), but pi(x) not a
    polynomial over fields. Tau conjecture orthogonal.
(h) **Systematic analysis of 8 intermediate quantity families**: residues, polynomial
    evals, matrix eigenvalues, topology, representation theory, entropy, recursive,
    physical -- ALL route back to floor values or zeta zeros.
(i) **No 2026 breakthroughs**: Chen-Tal-Wang STOC 2026 closest to TC^0 frontier.

## Session 16
(a) **15 intermediate quantity families now CLOSED** (8 from S15 + 4 GapL + 7 novel):
    class numbers h(-d), L-function L(1,chi), elliptic curve a_p, regulators,
    additive combinatorics/sumsets, ergodic theory, model theory, tropical geometry,
    sufficient statistics, algebraic geometry/F_q, representation theory S_n/GL_n.
    ALL route to primes (C), zeta zeros (E), or lose information (I).
(b) **Three "pillars" are the ONLY exact encodings of pi(x)**: prime positions,
    zeta zeros, floor values {floor(x/k)}. These are informationally equivalent.
    No fourth encoding found across 15+ candidate families.
(c) **TC^0 batch counting has 5 closed routes**: MAJORITY fan-in, divide-and-conquer,
    batch CRT, Carmichael structure, period exploitation.
    Communication complexity gives 2^{N/2}/poly(N) lower bound on TC^0 circuit size.
(d) **H-T signed cancellation transfer to pi(x) is CLOSED via 6 routes**:
    POSITIVITY of prime indicator prevents cancellation. Converting M(x)->pi(x)
    always costs O(x^{2/3}).
(e) **Lambert W error is structurally random**: delta(n) uncorrelated with gaps,
    uniform mod m, ~sqrt(x) magnitude. No exploitable structure.
(f) **HKM 2023 achieves O~(sqrt(x)) for pi(x) elementarily** (Math. Comp. 2024).
(g) **Aggarwal 2025 gives O(sqrt(n)*log^4(n)) for p(n)** via binary search + HKM.

## Session 17
(a) **EXACT communication complexity of pi(x)**: rank(pi_N) = 2^{N/2-1} + 2 for
    balanced bit partition (verified N=2..20). This gives dc(pi_N) >= Omega(sqrt(x)).
    The multilinear polynomial route to GapL is DEFINITIVELY CLOSED.
(b) **Boolean Fourier analysis**: Prime indicator has ~30% excess low-degree Fourier
    weight vs random (from parity/mod-4 structure), but noise sensitivity is near-random
    (~0.9x). No evidence of low-depth circuit structure. Spectral profile near-random.
(c) **Ono partition characterization (PNAS 2024) CLOSED**: requires divisors = circular.
(d) **Chen-Tal-Wang ECCC 2026**: THR-of-THR lower bounds, not number-theoretic.
(e) **sqrt(x) barrier is UNIVERSAL**: communication complexity, Fourier analysis,
    determinantal complexity, substitution rank -- ALL converge to sqrt(x).

## Session 18
(a) **Derandomization theory CLOSED (6 routes):** IW97, NW generators, approximation
    amplification, KI04/PIT, reverse hardness amplification, batch Fermat derandomization.
(b) **Natural proofs barrier (RR97) CONFIRMED relevant:** Cannot prove pi(x) not in
    TC^0 via "natural" methods. All our Fourier/rank measures are natural.
(c) **META-INSIGHT:** Derandomization addresses removing randomness from EXISTING
    efficient algorithms. Our problem is more fundamental: no efficient algorithm exists.
(d) **Formula complexity:** KW theorem gives formula size >= 2^{N/2-O(1)} for pi(x).
    Via Valiant, general circuit size >= 2^{N/4}.
    TC^0 bound 2^{N/2}/poly(N) remains strongest for constant-depth.
