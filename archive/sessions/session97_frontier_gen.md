# Session 97: Frontier Generation — 5 New Cross-Domain Attack Vectors

**Mode:** frontier_gen (auto-fired by run.sh).
**Date:** 2026-04-27.
**Self-grade:** **B**.

## Summary

Produced 5 new ATTACK_VECTORS.md entries (§D9–§D13), each grounded in
a distinct UNUSED cross-domain technique from CROSS_DOMAIN_TECHNIQUES.md.
All five have a single-session concrete first step, a pre-stated
falsification criterion, an A-grade success target, and a B-grade
fallback. CROSS_DOMAIN_TECHNIQUES.md updated to mark all five
techniques PROPOSED with survey references.

## The 5 new vectors

### D9 — Sum-product gain on the prime set in F_p (Bourgain-Glibichuk-Konyagin)

Tests whether the prime set `A_prime ⊂ F_p` saturates the
BGK 2006 sum-product lower bound `max(|A+A|, |A·A|) ≥ |A|^{15/14}`
or has *strictly* greater joint additive-multiplicative gain than
matched-density random subsets. UNUSED in §7. **A-grade if** primes
have a gain that grows with `p` and admits an HL-mod-q
interpretation; **B-grade if** matched (38th pseudorandomness
measure of joint-structure type — strongest such category).

### D10 — Mahler measure m(f_N) of f_N(z) = Σ_{n≤N} χ_P(n) z^n

Tests whether the prime generating polynomial has near-cyclotomic
factorisation (`m(f_N) = O(log^c N)`) or Lehmer-typical scaling
(`m(f_N) ~ √N`). UNUSED in §2 algebraic-height. **A-grade if**
near-cyclotomic — gives polylog π(N) via roots-of-unity sampling;
**B-grade if** scaling distinguishable from random matched-density
0/1 polynomial baseline (first project measure detecting algebraic-
height structure of χ_P).

### D11 — Shadow tomography sample complexity for π(x) extraction

Aaronson 2018 / Huang-Kueng-Preskill 2020 give `poly(log M)` query
complexity for many properties of a quantum state. Question: under
the random-Rademacher-shadow query model on `χ_P`, how many queries
suffice to estimate `π(M)` for ALL `M ≤ N` simultaneously? **A-grade
if** `K = poly(log N)` (polylog QUERY complexity for π — new
computational-model result); **B-grade if** explicit `K ≥ Ω(N^β)`
lower bound (strengthens information-theoretic bounds with a query-
complexity bound). UNUSED in §8.

### D12 — Compressed sensing of χ_P in arithmetic-progression dictionary

Tests whether `χ_P ∈ {0, 1}^N` is `K`-sparse in a structured
AP dictionary (`{1_{aZ + b}}` for sieve moduli `a ∈ {2, 3, 5, 7, ...,
30030}`). **A-grade if** OMP / L1-min finds `K = poly(log N)` —
gives a finite-`N` polylog representation of χ_P, hence polylog
π(M) at the project's target scale; **B-grade if** `K_ε = Θ(π(N))`
rigorously, closing this dictionary family. UNUSED in §6.

### D13 — Subword complexity p_w(n) of χ_P as binary infinite word

Computes `p_w(n)` for `w = (χ_P(n))_{n≥2}` at `n ∈ [1, 30]` via
sliding-window factor counting at `N = 10^7`. **A-grade if**
`p_w(n)` provably bounded by `2^n · exp(-c√n)` matching an HL
prediction; **B-grade if** statistically distinguishable from
matched-density Bernoulli (first topological-dynamics deviation of
χ_P from random). UNUSED in §5 symbolic dynamics.

## Cross-domain literature consulted

- **Sum-product (D9):** WebSearch retrieved Bourgain-Glibichuk-
  Konyagin 2006 J. London Math. Soc. (the explicit `|A|^{15/14}/(log
  |A|)^{2/7}` exponent for `|A| ≤ p^{7/13}`), Garaev 2008 arXiv:math/
  0702780 ("An explicit sum-product estimate"), Iosevich 2009
  arXiv:0904.2075 ("Sum-product phenomena in F_p"). Tao-Vu 2006 textbook
  for general framework.

- **Shadow tomography (D11):** WebFetched Aaronson 2018 arXiv:1711.01053
  abstract — `Õ(ε^{-4} log^4 M · log D)` copies for M observables.
  Huang-Kueng-Preskill 2020 Nature Physics 16, 1050 (classical shadows)
  cited in CROSS_DOMAIN_TECHNIQUES.md as the classical analogue.

- **Mahler measure (D10):** WebFetched Wikipedia article (definition
  via Jensen integral, Kronecker's theorem, Lehmer's conjecture,
  Weil height connection). Smyth 2008 survey URL added to refs.
  Boyd 1981 / Lehmer 1933 cited.

- **Compressed sensing (D12):** WebFetched Wikipedia RIP article
  (verified Gaussian/Bernoulli/partial-Fourier RIP, no general
  dictionary RIP for AP-structured matrices). Candes-Tao 2006 / Foucart-
  Rauhut 2013 cited.

- **Subword complexity (D13):** WebFetched Wikipedia "Complexity
  function" article (Morse-Hedlund classification, Sturmian p(n) =
  n+1, topological entropy). Lind-Marcus 1995 / Cassaigne-Nicolas
  2010 cited.

## Cross-checks against CLOSED_PATHS.md

Every candidate was checked against ~750 closures. **Key avoidances:**

- **Selberg trace formula** (lines 200, 347, 520, 593, 653, 656, 716,
  739) — explicitly closed as "isomorphic to explicit formula".
  Excluded.
- **Sheaf cohomology** (line 202) — closed E. Excluded.
- **Tropical geometry** (lines 199, 318, 423, 424) — closed C+E+I,
  full battery. Excluded.
- **Goldreich-Levin / list decoding** (line 473) — closed for π(x)
  on R(x)-low-bit decomposition. Distinct from D11 (shadow tomography
  is not GL — different query model and different bounding technique).
- **FRACTRAN / transfer operators / Furstenberg** (lines 179, 421,
  593, 656) — closed E+C as "spectral theory = zeta zeros".
  D13 (subword complexity) is purely combinatorial-on-words
  (topological), distinct from measure-theoretic transfer-operator
  spectral theory. Explicitly noted in vector D13 text.
- **Trace formula moment method** (lines 520, 653, 656) — closed.
  Mahler measure (D10) is height-theoretic, not spectral / trace —
  orthogonal mathematical category.
- **p-adic interpolation** (lines 8, 10) — closed as "Mahler
  coefficients grow super-exponentially". This is a SPECIFIC
  closure of the von-Mangoldt-side p-adic Mahler series. D10 is the
  COMPLEX Mahler MEASURE of the prime generating polynomial (Jensen
  integral) — completely distinct invariant in a different complete
  field.
- **Fourier dictionary** (line 47) — closed. D12 uses STRUCTURED
  AP dictionary (sieve moduli), not Fourier — distinct as discussed.

## Self-grade reasoning

**B not A:** none of the five vectors is a confident "≥ 2 will
produce A-grade" bet. Each has a *plausible* pathway to A-grade
(D9 if HL gives sum-product bias; D10 if cyclotomic factors abound;
D11 if random-Rademacher shadow inverts globally; D12 if AP-sparse
representation exists; D13 if HL constrains forbidden factors), but
each has equally plausible E-mode closure (38th pseudorandomness
measure type). Honest expectation: ~3-4 will close as B-grade
"new pseudorandomness measure in NEW mathematical category", 1
will close as duplicate-plus, 0-1 might surprise to A-grade.

The 5 vectors collectively SPAN 5 distinct mathematical categories
not represented in the existing 35-measure battery:
1. **Joint additive-multiplicative combinatorics** (D9).
2. **Algebraic height / Mahler measure** (D10).
3. **Quantum-info-style sample complexity** (D11).
4. **Structured-dictionary compressive sensing** (D12).
5. **Topological / symbolic dynamics on infinite words** (D13).

Even all-B closures would advance the project: the pseudorandomness
battery would expand to 40 measures across 9-10 mathematical
categories rather than 35 across 5, with each new category being
harder to dismiss as "obviously HL-isomorphic". The S87 / S88 / S92 /
S95 / S96 pattern (each new mathematical category producing a fresh
HL-detection edge) suggests every new orthogonal category is a real
project step forward — no longer pure duplicate.

**Not C:** the proposed vectors cite genuinely-unused techniques
(none has been imported in any prior session per CROSS_DOMAIN_TECHNIQUES.md
status). Each has a single concrete first step (no "study X" non-vectors).
Each has a pre-stated falsification criterion. Five WebFetched survey
references (one failed initially, replaced by alternative). No vectors
proposed without grounding.

**Not A:** none has the "I expect this to produce paper-grade fresh
result" confidence. The rotation should fire D11 (shadow tomography)
first if the harness picks at random — it has the cleanest A-grade
pathway (genuinely new computational model) and the lowest pre-test
duplicate risk. D9 (sum-product) is the lowest-risk B-grade candidate;
D13 (subword complexity) is the cheapest to execute first.

## Suggested rotation ordering

If the harness offers a choice in the next 5 production-mode sessions:
1. **D11 (shadow tomography)** — newest computational model; 1-2 sessions.
2. **D9 (sum-product on F_p)** — cheapest substantive computation; 1 session.
3. **D13 (subword complexity)** — single-session, single-script, B-grade likely.
4. **D10 (Mahler measure)** — single-session Jensen-FFT; testable scaling law.
5. **D12 (compressed sensing AP)** — most likely to require 2 sessions.

## Next-action

The next production-mode session should either pick one of D9-D13 or
attempt §G1 / §G2 (Liouville-side multiplicative regime, also
single-session). Both classes are now well-stocked. The project's
ATTACK_VECTORS frontier is no longer thin — `frontier_gen` should NOT
fire again until 2+ A-grade-misses have been recorded among D9-D13 +
G1-G3 + the existing open vectors.
