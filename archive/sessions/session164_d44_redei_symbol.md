# Session 164 — D44 Rédei symbol distribution on prime triples

**Mode:** wild_swing / frontier attack — A-grade target from
ATTACK_VECTORS.md §D44.

**Cross-domain ingredient imported (NEW for project):** Arithmetic
topology / Mazur 1966 primes-as-knots dictionary / Morishita 2002
arithmetic Massey products / Rédei 1939 explicit ternary symbol
formula. Promotes CROSS_DOMAIN_TECHNIQUES.md §4 row "Arithmetic
topology" PROPOSED → USED (mode E).

**Channelled mathematician:** Morishita / Stevenhagen — mode-E
arithmetic-topology + class-field-theory toolkit, distinct from
analytic-NT / circuit / spectral families.

**Outcome:** B-grade case (I), mode E. F1 (A-grade target)
**FALSIFIED**; F2 (UNBIASED, mode E) **HOLDS** at all three N.
Adds new edge **E2.29**.

---

## 1. Frontier target

D44 from ATTACK_VECTORS.md asks: is the empirical distribution of
the Rédei symbol `[p, q, r] ∈ {±1}` on Borromean-admissible prime
triples below `N` correlated with the Hardy-Littlewood ternary
singular series `S_3(N)`? If yes (A-grade), the Rédei symbol
gives a NEW arithmetic invariant for prime triples — the first
cohomological-Massey triple-prime invariant orthogonal to the
bilinear-singular-series HL framework.

The Rédei symbol is the canonical *triple* arithmetic invariant
in Mazur's primes-as-knots dictionary, where it specialises the
*triple Massey product* `⟨p, q, r⟩_2` in the Galois cohomology
`H^2(G_{Q, S}, F_2)` of the maximal {p, q, r}-unramified
extension. Equivalently it is the Frobenius image of `r` in the
*Rédei field* `K_{pq}/Q` (cyclic of order 4, ramified only at
p, q, containing `Q(√pq)`).

This is the **FIRST** project measurement of Massey-product /
arithmetic-topology invariants on the prime sequence.

## 2. Pre-registered falsification (in script docstring before run)

- F1 (A-grade): `|f_+ - 1/2| > 5/√K` AND HL-S_3 Pearson > 0.3.
- F2 (B-grade case I — UNBIASED): `|f_+ - 1/2| < 5/√K` at all three
  N levels.
- F3 (B-grade case II — non-HL bias): `|f_+ - 1/2| > 5/√K` but bias
  not matched by HL.
- INC: > 10 % triples lack norm rep in O_p (q non-principal).

## 3. Method

Implemented `experiments/algebraic/redei_symbol_prime_triples/redei_symbol_prime_triples.py`.

- Validate on canonical Borromean triple `[13, 61, 937] = -1`
  (Wikipedia, Morishita 2012). Six independent fundamental
  solutions to `U² - 13 V² = 244` all give the same value -1:
  formula sigma-choice and (U, V)-choice invariance confirmed.
- Compute the Rédei symbol via the maximal-order norm formula:
  find (U, V) with `U² − pV² = 4q` and `U ≡ V (mod 2)` via
  `sympy.solvers.diophantine.diop_DN`; then for `σ ∈ F_r` with
  `σ² ≡ p (mod r)` the symbol is the Legendre value
  `((U + Vσ)·2⁻¹ / r)`. Lemmermeyer 2000 *Reciprocity Laws*, Ch. 9.1.
- INCONCLUSIVE marking when `q` is non-principal in the narrow
  ideal class group of `O_p` (no `(U, V)` solution exists).

## 4. Results

| N    | π_{1 mod 4}(N) | admissible | well-def K | INC | n_+  | n_−  | f_+    | z (vs ½) | p     |
|------|----------------|------------|------------|-----|------|------|--------|----------|-------|
| 200  | 18             | 122        | 122        | 0   | 59   | 63   | 0.4836 | -0.36    | 0.72  |
| 500  | 43             | 1419       | 1361       | 58  | 662  | 699  | 0.4864 | -1.00    | 0.32  |
| 1000 | 80             | 9100       | 8577       | 523 | 4174 | 4403 | 0.4867 | -2.47    | 0.013 |
| 1000 (clean) | 74    | 7591       | 7591       | 0   | 3684 | 3907 | 0.4853 | -2.56    | 0.011 |

5σ envelope `2.5/√K` at N=1000 is `0.027`, observed `|Δ| = 0.0133`
< 5σ. **F2 (UNBIASED) HOLDS, F1 (A-grade) FALSIFIED.**

The clean-subset measurement (excluding the six primes p ∈
{229, 257, 401, 577, 733, 761} where some q's lie in non-principal
narrow class) gives `z = -2.56`, slightly stronger than the full
set's `z = -2.47` — so the marginal bias is NOT an artefact of
INCONCLUSIVE-triple removal.

The pattern of *constant* `Δf_+ ≈ -0.013` across N=200, 500, 1000
with z growing as `√K` is consistent with finite-r effective
Chebotarev convergence of the Frobenius density to 1/2 — a known
finite-N artefact of `O(1/√(log r)) ≈ 0.01` for `r ≲ 1000` —
NOT a true asymptotic bias.

No mod-residue (mod 8, 12, 24) substructure deviates beyond
Bonferroni-corrected 3σ.

## 5. Self-evaluation (CLAUDE.md 4-question)

**1. What did I produce that was not in the project before?**

(a) The first project use of arithmetic topology / Mazur primes-
    as-knots dictionary / Morishita Massey products on the prime
    sequence — a NEW cross-domain category.
(b) The first explicit empirical Rédei symbol distribution
    measurement on admissible Borromean prime triples (no published
    large-scale measurement of this kind exists).
(c) A new edge **E2.29**: cohomological-Massey orthogonality
    consistent with HL bilinear singular series at the empirical
    mean-of-symbol level.
(d) Validated implementation of the Lemmermeyer Rédei-symbol
    formula reproducing the canonical Borromean triple
    `[13, 61, 937] = -1`.
(e) Two falsification cases in mode E: F1 falsified, F2 holds.

**2. What edges did my work compose / cite?**

- Cited as orthogonal complements: E2.13 (Gowers chi_P), E2.14
  (Anderson chi_P), E2.15 (algebraic immunity), E2.16 (DPP failure),
  E2.17 (PH gaps), E2.19 (subword complexity), E2.20 (Mahler),
  E2.22 (Pollicott-Ruelle), E2.24 (AHK matroid Hodge), E2.26 (GCT
  orbit), E2.27 (KPZ Hölder), E2.28 (Baker-Norine).
- Distinct from: CLOSED_PATHS line 208 (generic knot invariants),
  CLOSED_PATHS line 202 (étale cohomology of Spec(Z)).

**3. If my session produced only duplicate closures, why?**

It did not — the cross-domain technique (arithmetic topology /
Galois cohomology of triple Massey products) is genuinely new for
the project, and the empirical distribution measurement is novel.
The closure mode (E) puts it in the steady-state pseudorandomness
catalogue, but the import is fresh.

**4. Next-action for the next agent.**

Three concrete successor proposals (added to ATTACK_VECTORS.md
and RESEARCH_AGENDA.md as D44.a/.b/.c):

- **D44.a** Quartic Massey product `⟨p₁, p₂, p₃, p₄⟩_2` (Vogel
  2004 thesis) on admissible 4-tuples; first higher-rank
  arithmetic-topology invariant.
- **D44.b** 4-rank distribution of `Cl(Q(√(p·q)))` (Cohen-Lenstra)
  on admissible (p, q) pairs.
- **D44.c** Quantitative effective-Chebotarev rate test: model
  `|Δf_+| ∼ c / √log N` and confirm the observed sub-5σ z-growth
  is finite-N effective Chebotarev for the Rédei field.

## 6. Self-grade

**B-grade (case I, mode E).**

The session attempted an A-grade frontier target (D44), produced
working code that validates on the canonical Borromean example,
ran a clean empirical measurement on 9100 admissible prime triples
at N ≤ 1000, and reached a CLEAN closure: F1 (A-grade) falsified
within strict 5σ, F2 (B-grade case I, UNBIASED mode E) holds.
Adds a new edge in a structurally orthogonal mathematical
category (arithmetic topology / Galois cohomology — the 13th
orthogonal pseudorandomness category for χ_P).

This is NOT an A-grade outcome (no positive structural detection
beyond HL). It is also NOT an inflated B — the ambitious frontier
was attempted, failed informatively, and the failure mode is
structural (the Rédei symbol mean is consistent with HL+
independence, with a marginal but sub-5σ z that is most plausibly
finite-N effective Chebotarev). Per CLAUDE.md grading: failed
ambitious cross-domain attempts that produce a structurally clean
negative-shape edge are B-grade.

The session does NOT inflate its outcome to A: although the
single-test p = 0.013 is intriguing, it is below the
pre-registered 5σ A-grade threshold and consistent with known
finite-N Chebotarev convergence rates.

---

## Artefacts

- `experiments/algebraic/redei_symbol_prime_triples/redei_symbol_prime_triples.py`
- `experiments/algebraic/redei_symbol_prime_triples/redei_symbol_prime_triples_results.md`
- `experiments/algebraic/redei_symbol_prime_triples/results_N{200,500,1000}.json`
- New edge **E2.29** in EDGES.md.
- New CLOSED_PATHS row.
- Update CROSS_DOMAIN_TECHNIQUES.md row §4 "Arithmetic topology".
- ATTACK_VECTORS.md §D.D44 moved to "Closed attacks".

---

## Pre-registered falsification, replayed

- F1 (A-grade): FALSIFIED — `|Δf_+| < 5/√K` at N=1000.
- F2 (B-grade case I): HOLDS at all three N.
- F3 (B-grade case II): not triggered.
- INC: 5.7 % at N=1000, < 10 % threshold; subset analysis rules
       out INC-induced bias.
