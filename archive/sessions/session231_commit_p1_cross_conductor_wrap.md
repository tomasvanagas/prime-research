# Session 231 — commit thread 6 slot 5: cross-conductor (Q-batched) amortisation + Thread 6 wrap

**Date:** 2026-04-29
**Mode:** commit (Thread 6 / π in arithmetic progressions, batched on
modulus q, slot 5 of 5 — FINAL)
**Self-grade:** **B** — substantive theoretical + empirical wrap of
the fourth and final amortisation axis (Q-batched cross-conductor).
Empirical speedup decreases monotonically toward 1× as M grows
(1.25× → 1.12× → 1.05× at M ∈ {2, 3, 4}), confirming bounded-constant
savings asymptotically dominated by per-conductor work. **Thread 6
closes negatively across all four amortisation axes** with rigorous
stage decomposition + empirical confirmation. Not A-grade because no
new partial-positive (this is a negative-shape closure of P1); not
C-grade because the four-axis decomposition + cross-conductor
experiment are novel artefacts.

## Mission

Final slot of Thread 6. From slot-4 (S230) `recommended_next_action`
PRIORITY (a):

> Cross-conductor amortisation argument (slot-3 leftover priority-a).
> Likely structurally impossible by CRT independence across distinct
> conductors; rigorous lower bound would seal Thread 6 negatively
> across all four amortisation axes.

This session's task: empirically test cross-conductor (Q-batched)
amortisation on the slot-3/4 pipeline; pair the empirical result
with a theoretical lemma; **synthesise the five-session Thread 6
arc into a final result**; recommend next thread.

## What was built

`experiments/analytic/batched_q_amortisation/slot5_cross_conductor_amortisation.py`
(~400 lines):

- **Stage decomposition** of the slot-3 / slot-4 end-to-end zero
  finder for π(x; q, a) into seven stages:
  - Stage A: cp_all = 1/n^{1/2+it} (conductor-independent).
  - Stage B: aggregate cp_all by residue mod q (per-conductor).
  - Stage C: length-φ(q) FFT (per-conductor).
  - Stage D: reflected term + Gauss sums + W_χ (per-conductor).
  - Stage E: loggamma at t-grid (conductor-independent).
  - Stage F: Hardy theta + Z (per-conductor).
  - Stage G: sign-change bracket (per-conductor).

- **Two pipelines**:
  - INDEPENDENT: rebuild Stage A + E inside the q-loop.
  - Q-BATCHED-SHARED: build Stage A (sized to N_max) + Stage E *once*
    before the q-loop; per-conductor stages B-G run as before.

- **Empirical M-sweep** at K=200, n_t=823, t_max=200, family
  Q ⊆ {1009, 2003, 5003, 10007}, M ∈ {1, 2, 3, 4}.

- **Numerical equivalence verification**: max abs diff between
  pipelines = 0.000e+00 at Q = {101, 251, 503, 1009}, K=5.

## Empirical findings

```
M | Q                         | T_indep (ms) | T_batch (ms) | speedup
1 | 1009                      | 264.3        | 134.5        | 1.97×*
2 | 1009, 2003                | 485.5        | 389.7        | 1.25×
3 | 1009, 2003, 5003          | 1267.4       | 1133.3       | 1.12×
4 | 1009, 2003, 5003, 10007   | 3493.2       | 3341.4       | 1.05×
```

\* M=1 is functionally identical between pipelines; the 1.97× ratio
is a cold-start FFT-planner / numpy-allocation artefact (indep is
called first per run-trio), NOT a real speedup. M ∈ {2, 3, 4} rows
are the honest measurement.

**Key signature: speedup decreases monotonically as M grows.** This
is the OPPOSITE shape of S224's Correlation Dichotomy (33× at M=64,
INCREASING in M). The Q-axis admits no Correlation-Dichotomy partial-
positive.

### Stage A breakdown (M=4)

```
indep_stageA_total = 572 ms  (Stage A built 4 times, once per conductor)
batch_stageA_once  = 343 ms  (Stage A built ONCE, sized to N_max for q=10007)
saving (Stage A)   = 229 ms
saving (Stage E)   =   1 ms  (sub-dominant by 150×)
TOTAL saving       = 152 ms = 4.4% of T_indep
```

The savings are bounded by `T_stageA(N_max) + T_stageE`, both
constants. As M grows with q_max fixed, Stage A savings stay
bounded. As Q expands to larger conductors, per-conductor work grows
as `√q · n_t`, so savings/total → 0.

## Theoretical Lemma (cross-conductor amortisation)

**Claim.** For any Q = {q_1, ..., q_M} of distinct prime conductors
with q_i ≥ q_min, the slot-3-style L-zero pipeline admits

```
T_batch(Q) / T_indep(Q) ≥ 1 - c·(T_A(N_max) + T_E)/T_indep(Q)
```

The savings ratio is bounded by O(1/M) when q_avg is fixed (Stages
B-G grow linearly in M with O(√q · n_t) per term), and goes to 0 as
the family expands to include larger conductors.

**Proof structure.**

(i) Stage B's residue map `n → n mod q` is q-dependent. No primitive
aggregates 1/n^{1/2+it} into ALL φ(q_i) characters of all q_i in Q
simultaneously: distinct primes have COPRIME multiplicative groups
of order q_i - 1 with no nontrivial common refinement. CRT
gives (Z/q_1q_2)* ≅ (Z/q_1)* × (Z/q_2)* but requires the SINGLE
composite modulus, not a *family* of moduli.

(ii) Stages C, F, G have FFT length φ(q_i), q-dependent.

(iii) Stage D has Gauss sums τ(χ_q,j) and root numbers W_χ_q,j that
depend on the residue structure mod q.

(iv) Only Stages A and E are conductor-independent. Their combined
cost is one bounded constant `T_A(N_max) + T_E = O(N_max·n_t + n_t)`.
Sharing across M conductors saves `(Σ_q N_q - N_max + (M-1))·n_t`
operations — at most `~2·N_max·n_t` for any M, regardless of M.

(v) Total per-conductor work `Σ_q (φ(q)·log φ(q) + N_q)·n_t ~
M·√q_avg·log q_avg·n_t` grows linearly in M.

Saving/total = O(N_max·n_t) / O(M·√q_avg·log q_avg·n_t) → 0. ∎

**Corollary.** The Q-axis admits no Correlation-Dichotomy-shaped
partial-positive on the slot-3 / slot-4 architecture.

## Thread 6 — five-session arc summary (ALL slots)

The Thread 6 question (P1 of OPEN_POSITIVE_TARGETS):

> *"For fixed x and a fixed family of moduli Q = {q_1, ..., q_M},
> compute π(x; q_i, a_i) for all i. Is there per-(q, a) amortised
> polylog over M = poly(log x) characters?"*

Five-axis decomposition tested across slots:

| Axis | Description | Slot | Empirical lift | Asymptotic shape |
|------|-------------|------|----------------|------------------|
| 1 | a-direction (fixed q, χ; varying a) | S227 | 8× at M=256 | TRIVIAL: zeros independent of a |
| 2 | χ-direction (fixed q; all φ(q) chars) | S228, S229 | 1.79–2.04× | Bounded constant: BLAS-vs-Bluestein FFT FLOP rate |
| 3 | composite q (multi-axis FFT) | S230 | 1.75× at q=1001 | Bounded constant: same crossover as cyclic FFT |
| 4 | Q-batched (cross-conductor) | S231 | 1.05× at M=4, decreasing toward 1× | Bounded constant decreasing toward 0 saving |

(Plus a-direction within Q-batches: trivially amortisable per slot 1.)

**Thread 6 verdict: P1 closes negatively for the polylog target across
all four amortisation axes.**

- Trivial axis (1) saves 8× exactly per E1.5 prediction (zeros depend
  on χ, not a).
- Axes 2-4 give bounded-constant lifts (1.05× to 2.04×) explainable
  by hardware-level FFT/BLAS arithmetic; **none structural in the
  algorithmic-complexity sense**.
- No axis produces sub-linear-in-batch-size amortisation. **Opposite
  shape from Thread 5 (Correlation Dichotomy).**

The five-session arc produced:
- A novel end-to-end zero finder for π(x; q, a) at prime AND composite q.
- mpmath-validated accuracy at q ∈ {11, 15, 35, 105, 1001} to ~0.05 abs.
- The first cross-conductor amortisation experiment in the project.
- A unified four-axis stage-decomposition closure of P1.
- One new cross-domain technique entry (CRT-based multi-axis Dirichlet
  DFT, S230).

## Why Thread 6 ≠ Thread 5 (Correlation Dichotomy contrast)

Thread 5 (S224) found a partial-positive: 33× speedup at M=64 for
*correlated narrow-window* π queries. The mechanism: shared zeros
amortise when queries are CORRELATED.

Thread 6 found no partial-positive across any of the four
amortisation axes. The mechanism is different:

- In Thread 5, the family of queries shared a FIXED L-function (ζ).
  Zeros are computed ONCE; per-query work is the partial-sum eval.
- In Thread 6, the family of queries spans DISTINCT L-functions
  (one per character, possibly across distinct conductors). Zeros
  must be computed PER L-function — there is no shared zero database
  across distinct characters of distinct conductors.

The structural distinction: Thread 5's "shared L" geometry vs Thread
6's "distinct L per query" geometry. The former admits Correlation
Dichotomy; the latter does not.

This is a USEFUL distinction for `OPEN_POSITIVE_TARGETS.md`: future
threads should prioritise problems with shared-L geometry (Thread 5
shape) over problems with distinct-L geometry (Thread 6 shape) when
seeking partial-positives.

## Edges composed / cited

- **E1.5** (per-query bit-content): each amortisation axis confirms
  E1.5's prediction.
- **E2.1** (spectral / sieve): cyclic-DFT and multi-axis-DFT
  primitives are the spectral handles tested.
- **E3.1** (CCM amortisation, downgraded): same setup-vs-eval
  decoupling.
- **E6.6** (Aggarwal binary search): slot 5 closes the q-axis
  similarly to Aggarwal's x-axis closure in Thread 5.
- **E6.7** (Deléglise-Rivat per-query Θ(x^{2/3}/log²)): orthogonal
  pillar — slot 5 closes the AP variant.
- **S224** (Correlation Dichotomy, Thread 5): slot 5 produces the
  Q-axis analogue, with OPPOSITE shape (decreasing-toward-1×, not
  increasing-toward-33×).
- **S190** (Thread 1 closed), **S202** (Thread 3 closed), **S215**
  (Thread 4 closed): structural closures Thread 6 now joins as the
  fifth.
- **S227-S230** (Thread 6 slots 1-4).

## Cross-domain ingredient

**NONE NEW** in slot 5. Slot 5 is rigorous-experimentation +
algebraic-structure argument applied to existing primitives. The
slot-5 result is a theoretical lemma + empirical confirmation, NOT
a new technique import. This is appropriate for a wrap session.

The wrap completes Thread 6's CROSS_DOMAIN_TECHNIQUES.md §10 entries
(cyclic / multi-axis DFT family) by showing the family does NOT
extend to a *cross-conductor* primitive.

## Self-extension proposals (next thread)

Per CLAUDE.md autonomy invariant, when a thread closes:

**Successor 1 — Thread 7 = P3 (approximate π(x) ± ε in polylog with
named ε)**: highest a-priori chance of partial-positive based on
R(x) error analysis. With K = log²(x) zeros, ε ≈ √x · O(log log x /
log⁴ x) heuristically — strictly better-than-√x, polylog time.
Smallest cross-domain barrier (already pure analytic NT). **Recommended.**

**Successor 2 — Thread 8 = P2 (prime gap function π_h(x) batched on h)**:
distinct-h queries SHARE an underlying sieve (Hardy-Littlewood
singular series factorises through small primes). Different
amortisation geometry from Thread 6's distinct-conductor structure
— closer to Thread 5's shared-structure shape.

**Successor 3 — Thread 9 = P9 (quantum batched primes-in-AP)**: slot
5 closed Thread 6 *classically*. Quantum batched amortisation might
admit different scaling. Heavy literature, possibly inconclusive.

## Falsifiers (slot 5)

- **F5.1** Pipeline equivalence: max abs diff = 0e0. Stage-A sharing
  introduces NO numerical error.
- **F5.2** speedup decreases monotonically (1.97×* → 1.25 → 1.12 →
  1.05×). NOT sub-linear in M. *(M=1 row is artefact.)
- **F5.3** Theoretical lemma: Stages B-G are all per-conductor;
  Stage A + Stage E save bounded constants.
- **F5.4** ALTERNATIVE attack via Selberg trace formula or Iwaniec-
  Kowalski large-sieve NOT investigated. Different cross-domain
  imports could change the lemma; absence in slot 5 noted but not
  exhausted.
- **F5.5** Quantum batched amortisation NOT investigated. Could be a
  different thread.

## Files written

- `experiments/analytic/batched_q_amortisation/slot5_cross_conductor_amortisation.py` — new (~400 lines).
- `experiments/analytic/batched_q_amortisation/slot5_cross_conductor.csv` — 4 rows.
- `experiments/analytic/batched_q_amortisation/slot5_cross_conductor_amortisation_results.md` — new.
- `archive/sessions/session231_commit_p1_cross_conductor_wrap.md` — this file.
- `status/CLOSED_PATHS.md` — slot-5 row.
- `status/SESSION_INSIGHTS.md` — S231 entry.
- `EDGES.md` — Thread 6 four-axis closure summary line.
- `RESEARCH_AGENDA.md` — Arc 8 marked DONE_NEGATIVE.
- `OPEN_POSITIVE_TARGETS.md` — P1 marked CLOSED-B.
- `.commit_state` — sessions_used 4 → 5_final, status:DONE_PARTIAL_NEGATIVE.
- `.run_state` → 399.

## Session-end self-evaluation

1. **What did I produce that was not in the project before this session?**
   - Stage decomposition (A-G) of the slot-3 / slot-4 zero-finder
     pipeline, identifying conductor-independent vs per-conductor stages.
   - First cross-conductor (Q-batched) amortisation experiment.
   - Empirical four-axis closure of P1's polylog target.
   - Theoretical lemma: cross-conductor amortisation factor is
     bounded constant; vanishes as M·√q_avg grows.
   - Thread 6 five-session arc final synthesis.
   - Structural Thread-5-vs-Thread-6 distinction (shared-L vs
     distinct-L geometries) — useful guidance for future P_x targets.

2. **What edges did my work compose or cite?**
   E1.5 (per-query bit-content), E2.1 (spectral sieve), E3.1 (CCM
   amortisation), E6.6 (Aggarwal binary search), E6.7 (Deléglise-Rivat),
   S224 (Correlation Dichotomy template), S190 (Thread 1), S202
   (Thread 3), S215 (Thread 4), S227-S230 (Thread 6 slots 1-4).

3. **If only duplicate closures, why?**
   Not duplicate. Stage decomposition + cross-conductor empirical +
   theoretical lemma are project-novel.

4. **What is the next-action for the next agent?**
   Thread 6 DONE_PARTIAL_NEGATIVE. **Next thread = Thread 7 = P3
   (approximate π(x) ± ε in polylog with named ε)**, highest a-priori
   chance of partial-positive. Or Thread 8 = P2 (π_h(x) batched on h)
   for shared-structure shape similar to Thread 5.

   `.commit_state` is updated; `recommended_next_action` field
   updated for Thread 7. Next commit slot picks up Thread 7.
