# sharpP_probe.py — results

**Question (open item 3, the COMPUTATIONAL lower-bound face).** S511 closed the
*information* half: the Θ(√x) certificate is information-forced for the
sieve-reconstruction CLASS (the K=π(√x) checkpoints carry Θ(√x) joint hard bits).
The explicit remaining lever S511 named is the COMPUTATIONAL route — a *different*
(non-sieve) witness for `π(x)=c` is not ruled out by information (cf. factoring's
tiny witness for a hard answer). The decisive formal question is whether
`L_π = {(x,c) : π(x)=c}` is in NP at all, and where it sits in the counting
hierarchy (the #P side). This is a STANDALONE measurement + complexity anchor, not
a protocol.

Code: `sharpP_probe.py` (`--selftest`, 40 checks; `--kmax K` sets the ladder
top). Self-contained number theory (spf sieve, primality predicate, Legendre
inclusion-exclusion, n-th prime); no chain imports.

## Headline

1. **Upper bound (folklore, demonstrated exactly):** `π(x)` and the Legendre
   partial sieve `φ(x,a)` are **#P functions of the binary input x** — count the
   N=log x bit integers `n≤x` satisfying a poly(N)-time predicate (primality /
   coprimality). Hence `L_π ∈ C_=P` (exact counting) and
   `{(x,c):π(x)≥c} ∈ PP ⊆ P^#P`. (Toda: `PH ⊆ P^#P`, the *other* direction;
   `PP ⊆ PH` is not known — would collapse PH — so not claimed.)
2. **NP-completeness obstruction (provable):** `π` is a function, so
   `L_π ∈ NP ⇒ L_π ∈ NP∩coNP` (to certify `π(x)≠c`, exhibit the true `c'` with
   its NP-cert; `c'≠c`). Therefore **`L_π` NP-complete ⇒ NP=coNP.** The live
   question is mere NP-MEMBERSHIP (a polylog(x) witness), not NP-completeness.
3. **Witness-size ladder (measured):** every NATURAL certificate family for
   `π(x)=c` sits at one of three rungs (leading-power exponent α, model
   `log bits = α·log x + δ·log log x + c`):

   | witness family | α (leading power) | δ (polylog order) | rung |
   |---|---|---|---|
   | enumeration NP-cert (list π(x) primes w/ Pratt + every other n composite w/ a factor) | **0.985** | +1.16 | Θ(x) |
   | sieve transcript (the chain; S509 real naive=0.473) | **0.490** | +1.03 | Θ(√x) |
   | S511 joint-info floor | **0.486** | +0.12 | Θ(√x) |
   | zeta-zero / explicit-formula (Galway K~c·√x·log²x, **cited** EDGES Thread 3) | **0.500** | +3.00 | Θ(√x) |
   | hypothetical polylog witness (RULED OUT, S511) | **−0.000** | +2.00 | poly(N) |

   **THREE independent natural families (sieve, analytic/zeta, info floor)
   converge at the √x rung** (max |α−0.5| = 0.014 over k=10..20); enumeration sits
   a full power higher (Θ(x)); polylog is the rung S511 rules out for the sieve
   class. **No natural witness reaches poly(N); every one is ≥ 2^{N/2}.**
4. **Parsimonious-reduction (combinatorial #P-hardness) obstruction (precise +
   exact):** a parsimonious reduction `#A → π` maps `w ↦ x(w)` with
   `π(x(w))=#A(w)`. Realizing target count `c` forces `x ∈ [p_c, p_{c+1}-1]`, so
   the reduction map `c ↦ x(c)` **IS the inverse-prime `p(·)` = the project goal.**
   The reduction's "easy direction" is therefore as hard as π itself —
   combinatorial #P-hardness via target-count realization is **CIRCULAR** (closure
   mode C). And the sieve/φ "instance" is lattice-structured (sets are forced to
   be `{multiples of d}`); it has no instance richness to embed an arbitrary
   #P-complete count. **#P-hardness of π, if true, cannot come from
   instance-embedding — it would have to come from value-incompressibility (the
   S511 route).**
5. **Correction of CLOSED_PATHS row 175** ("…exact pi(x) is #P-hard"): that
   parenthetical is an **unsubstantiated S7-era assertion.** The true upper-bound
   statement is **π ∈ #P**; #P-HARDNESS is open and the natural combinatorial
   reductions fail (4). The BQP-exclusion conclusion row 175 reaches should rest
   on the √x information floor (S511) + no polylog quantum algorithm, NOT on an
   unproven #P-hardness.

## Detail

### (A) #P membership, constructive and exact
`π(200)` as a predicate-count = 46 = sieve π(200); `φ(200,4)` as a coprimality
count = 47 = the 2⁴-term inclusion-exclusion `Σ_{S⊆{2,3,5,7}}(−1)^|S|⌊200/∏S⌋`.
Both `π` and `φ` are #P (count the ≤2^N candidates `n≤x` under a poly(N)
predicate). The catch making this an *upper bound only*: the #P "subset witness"
for `φ(x,a)` has 2^a terms — for `a=π(√x)` that is `2^{Θ(√x/log x)}`, **doubly
exponential in N**. #P membership gives no short witness.

### (B) NP ⇒ NP∩coNP
Standard function-graph argument, illustrated at x=1000 (π=168): a poly(N) NP-cert
for `π(x)=c` immediately yields a coNP-cert for `π(x)≠c` (exhibit `c'=π(x)` with
its cert), so `L_π ∈ NP ⇒ NP∩coNP`, and NP-completeness would collapse NP=coNP.

### (C) the witness ladder
Bit-sizes vs x=2^k, k=10..20 (`--kmax 20`):

```
  k          x    pi(x)  K=pi(vx)          enum      sieve       info        zeta    poly
 10       1024      172        11         15493       3520        110       32000     100
 16      65536     6542        54       1608512      43200        864     1048576     256
 20    1048576    82025       172      32079562     213280       3440     8192000     400
```
- **enum** is computed from a real spf sieve: `Σ_{p≤x}` Pratt-cert bits (model
  `⌈log₂p⌉²`) `+ Σ_{composite n≤x}` smallest-factor bits — exact, dominated by the
  composite list (`≈ x·log x`), α=0.985.
- **sieve / info / zeta** are Θ(√x·polylog) proxies (K=π(√x) layers/checkpoints ×
  the appropriate polylog), continuous-log to avoid integer-bit_length step
  artifacts; their leading power α is the claim (the constants are illustrative).
- The leading-power fit `fit_power_log` separates α from the polylog order δ; the
  naive single-slope (shown in the report) is polylog-inflated over a short window
  (e.g. log²x gives naive slope 0.21 but α=0 — selftest case 6 verifies exact
  recovery on closed-form controls).

### (D) the reduction circularity
For target counts c∈{1,10,100,500}: p_c∈{2,29,541,3571}, with π(p_c)=c and
π(p_c−1)=c−1 exactly. The realizing window `[p_c, p_{c+1}-1]` is pinned by the
prime gaps; computing any x in it from a #SAT instance (without already knowing
the count) is the inverse-prime problem `p(n)`.

## What would falsify each claim

- **(A)** the predicate-count ≠ sieve π(x), or ≠ the inclusion-exclusion φ(x,a),
  at any tested x,a. (Not observed; selftest cases 3–4.)
- **(C)** the enumeration exponent not ~1.0; OR any of the three nontrivial
  natural families not fitting α~0.5 over the window; OR — the strong falsifier —
  ANY natural family yielding a polylog (α~0) exact witness, which would put
  `L_π ∈ NP` and contradict the S511 floor for the sieve class. (Selftest case 7:
  α_enum=0.98, α_sieve=0.49, α_info=0.49, α_zeta=0.50, α_poly=0.00; convergence
  max|α−0.5|<0.08.)
- **(D)** exhibiting a poly(log c)-time map `c ↦ x` with π(x)=c — that map is
  exactly the project's `p(n)` goal, so exhibiting it is a breakthrough, not a
  refutation of the obstruction's logic.

## Honest scope (what this is NOT)

- The ladder enumerates the **natural** witness families. A universal "`L_π ∉ NP`"
  needs ruling out ALL witnesses — that is the open #P-hardness question, NOT
  settled here. This probe shows the natural families are all ≥ √x and the
  combinatorial-reduction route is circular: a **precisely-stated obstruction**,
  not a non-membership proof.
- The zeta-zero rung is **cited** (Galway `K~c·√x·log²x`; EDGES Thread 3,
  S195/196 conditional, S434–436 worst-case-of-N strengthening), not recomputed —
  it is an established edge, included to show the analytic witness lands on the
  SAME √x rung as the sieve and the info floor. The probe does not re-measure it.
- The sieve / info / zeta bit-counts are **proxies** for the leading power, not
  the literal certificate constants. The exact chain certificate size was measured
  on the real protocol in S509 (naive exponent 0.473, dominated by `comm_outer`);
  this probe's contribution is the *ladder placement and the NP/#P framing*, not a
  remeasurement of the chain.
- "π ∈ #P" is folklore; the contribution is the **NP∩coNP reframing**, the
  **witness ladder**, the **reduction-circularity obstruction**, and the
  **correction of row 175**.

## Bottom line for open item 3

The two halves now agree and point the same way without closing the universal
question:
- **information side (S511):** the √x is forced for any sieve-reconstructing
  verifier;
- **computational side (this cycle):** `L_π ∈ #P` (upper), cannot be NP-complete
  unless NP=coNP, every natural witness family is ≥ √x (no poly(N) one known), and
  the combinatorial #P-hardness route is circular (so hardness, if real, is a
  value-incompressibility statement, i.e. the S511 face).

The remaining open sub-question (filed, not closed): a genuine #P-hardness proof
(or a non-sieve sub-√x witness) for π — neither is delivered by any natural route
examined here. `L_π` is, on present evidence, an **NP-intermediate-flavoured
counting problem with a √x certificate and no proven poly(N) one** — exactly the
status the verification program's Õ(√x) construction (S491–S509) realizes from
above and S511 pins from below.
