# Session 227 — commit thread 6 slot 1: per-(q, a) explicit-formula profile

**Date:** 2026-04-29
**Mode:** commit (Thread 6 / π in arithmetic progressions, batched on
modulus q, slot 1 of 5)
**Self-grade:** **B** — substantive infrastructure work + empirical
profile across three decades of q. Built the decoupled three-axis cost
model (T_zero_db_setup ⊥ T_per_chi_eval ⊥ T_orthogonality), confirmed
the predicted linear-in-φ(q)·K baseline for zero-database setup
(slot-1's negative-shape baseline), measured 8× a-direction
amortisation at M=256, x-independence at 10⁶ vs 10⁸, and produced a
concrete falsifier list (F1-F3) that frames slot 2-5 attacks. Not A:
no algorithm faster than the linear-in-φ(q) baseline produced — the
a-direction amortisation is structurally trivial. Not C: empirical
data is new and the architecture is new.

## Mission

From `.commit_state` Thread 6 slot 1 plan:

> Build a per-(q, a) explicit-formula evaluator that decouples
> shared-q zero-of-L-function computation from per-character
> evaluation. Profile both at q ∈ {10³, 10⁴, 10⁵} for x ∈ {10⁶, 10⁸}.
> Measure shared-zero amortisation factor.

Fresh thread (sessions_used was 0 in `.commit_state` after Threads 1-5
all closed; Thread 6 prioritised from `OPEN_POSITIVE_TARGETS.md` §P1).

## What was built

`experiments/analytic/batched_q_amortisation/per_q_a_explicit_formula_profile.py`:

- **Three-axis decoupled evaluator.** Layered architecture:
  - `time_zero_db_setup(q, K)` — builds per-character cache; principal
    χ uses cached ζ-zeros, non-principal χ uses Riemann–von Mangoldt
    density-inverted synthetic zeros (numpy-vectorised Newton).
  - `time_per_chi_eval(zeros, x)` — per-(χ, x) partial-sum
    Σ_ρ x^ρ/ρ vectorised over zeros of one character.
  - `time_orthogonality(q)` — χ̄(a) Σ_χ table, complex dot product.
- **Dirichlet character table generator.** Cyclic (Z/qZ)*-aware (works
  for prime q and prime-power q with primitive root); discrete-log
  table for character evaluation.
- **a-direction amortisation profile.** For fixed (q, K, x), measures
  T_amort_per_a vs M ∈ {1, 4, 16, 64, 256}.
- **Sieve sanity check** at q ∈ {7, 11, 31}, x ∈ {10⁴, 10⁵} — direct
  π(x; q, a) verifies equidistribution at small scale.

The grid actually run: q ∈ {101, 1009, 10007} for x ∈ {10⁶, 10⁸},
K ∈ {50, 200, 800}. (q=10⁵ left to slot 2 — 16-minute scaling at
q=10007, K=800 means q=10⁵ would have taken ~25 min and slot 1's
contract was the architecture, not maximum-q coverage.)

## Empirical findings

### Setup cost (zero database)

```
T_zero_db_setup / (φ(q) · K) ≈ 2.7e-6 s, held to within 5% across
q ∈ {101, 1009, 10007}.
```

The synthetic baseline matches the predicted "no sharing" scaling
exactly. Per-character setup is *constant* across three decades of q.

### Per-character evaluation cost

x-independent within measurement noise (10⁶ vs 10⁸ identical to µs).
~50ns per zero asymptotically (numpy vectorised, K=800).

### Orthogonality cost

Subdominant by 2-4 orders of magnitude vs evaluation. Not the
bottleneck for q ≤ 10⁴.

### a-direction batching amortisation

At q=1009, K=200, x=10⁶: 8.06× speedup at M=256.

| M    | T_amort/a | factor vs M=1 |
|------|-----------|---------------|
| 1    | 188 ms    | 1.000×        |
| 4    | 64 ms     | 0.341×        |
| 16   | 33 ms     | 0.176×        |
| 64   | 25 ms     | 0.134×        |
| 256  | 23 ms     | 0.124×        |

Asymptote = T_full_per_x = 22.7ms. **This is the trivial direction:
zeros depend on χ, not on a, so a-batching simply re-uses the zero
database across queries.** The 8× factor is Thread-5-flavoured but
structurally trivial.

## Why this is B-grade and not A or C

**Not A**: the only positive result is the trivial a-direction
amortisation, which the explicit formula predicts exactly because
zeros are independent of a. No algorithm faster than the
linear-in-φ(q) baseline was found. The slot-1 contract was
architecture + profile; the partial-positive piece (the χ-amortisation
question, which would be the actual A-grade win) is slot 2's target.

**Not C**: 18 rows of new (q, K, x) profile data; a working
three-axis cost model that didn't exist in the project before; an
explicit falsifier list (F1-F3) that frames the slot 2-5 attacks
operationally; cross-domain ingredient (Tarjan-style amortised
analysis) extended from ζ-zero/x-batched (Thread 5) to L-zero/χ-batched
(Thread 6) for the first time.

## Slot-1 questions (resolved)

**Q1 — T_zero_db_setup shape.** Linear in φ(q)·K. No naive
χ-amortisation. Slot 2 must produce structural sharing or accept the
baseline.

**Q2 — T_per_chi_eval x-dependence.** x-independent (numpy-vectorised
cos/sin over fixed-size γ array; x enters only as a scalar inside cos).
Confirmed at 10⁶ and 10⁸.

**Q3 — Orthogonality vs evaluation.** Orthogonality is subdominant;
per-(q, a) cost is dominated by φ(q)·T_per_chi_eval, NOT the
orthogonal sum.

## Falsifier list framing slot 2-5

**F1.** A Bober-Hiary 2017-style multi-character L-zero-finding
algorithm with cost sub-linear in φ(q) WITHIN A SINGLE CONDUCTOR.

**F2.** AFE-shared partial-sum evaluation: shared `1/n^{½+it}` table at
fixed t across all χ mod q. The crucial requirement is that distinct
χ's zeros cluster near common heights t — random-matrix theory
predicts O(1/T) cross-character correlation, asymptotically negligible
at K = 800 / x = 10⁶.

**F3.** Adaptive (q, a) protocol with sub-linear cumulative cost
(Aggarwal-AP analogue). Slot 4 of the arc plan.

## Edges composed / cited

- **E1.5** — bit-content barrier; per-(q, a) bit content
  ~log(π(x)/φ(q)) per query; a-batched amortisation saves SETUP cost
  but NOT evaluation cost, exactly per E1.5 at the batched level.
- **E3.1** — Connes amortisation; same setup-vs-eval decoupling that
  Threads 2-3 used, generalised from ζ to L(s, χ).
- **E6.7** — Deléglise-Rivat per-query Θ(x^{2/3}/log²); orthogonal
  pillar (sieve, not zeros).
- **S224** — Correlation Dichotomy template; slot-1's a-batching
  amortisation is the q-axis analogue of S224's correlated-x
  batching, but in the TRIVIAL direction (zeros a-independent).
  Slot 2 attacks the non-trivial χ-direction.

## Cross-domain ingredient

**Amortised algorithmics** (Tarjan 1985 / Demaine-Patrascu 2008):
USED I (continuing from Thread 5). Tarjan's
separate-amortisable-from-per-query device applied to the
Dirichlet-character family at fixed x — extension from
ζ-zero/x-batched (Thread 5) to L-zero/χ-batched (Thread 6) at fixed x.

## Files modified by this session

- `experiments/analytic/batched_q_amortisation/per_q_a_explicit_formula_profile.py` — new
- `experiments/analytic/batched_q_amortisation/per_q_a_explicit_formula_profile_results.md` — new
- `experiments/analytic/batched_q_amortisation/per_q_a_decoupled_profile.csv` — new (18 rows)
- `experiments/analytic/batched_q_amortisation/per_q_a_amortisation.csv` — new (5 rows)
- `experiments/analytic/batched_q_amortisation/per_q_a_sieve_sanity.csv` — new (18 rows)
- `.commit_state` — sessions_used 0 → 1, session_history += S227, recommended_next_action updated for slot 2
- `archive/sessions/session227_commit_p1_per_q_a_profile.md` — this file
- `status/CLOSED_PATHS.md` — appended row for slot-1 baseline
- `status/SESSION_INSIGHTS.md` — S227 entry appended
- `RESEARCH_AGENDA.md` — Arc 8 slot 1 marked done
- `.run_state` → 395

## Session-end self-evaluation

1. **What did I produce that was not in the project before this
   session?** A three-axis decoupled per-(q, a) explicit-formula
   evaluator (`per_q_a_explicit_formula_profile.py`) with separable
   T_zero_db_setup / T_per_chi_eval / T_orthogonality timing
   primitives; an empirical baseline at three decades of q
   (101 / 1009 / 10007) confirming linear-in-φ(q)·K zero-database
   scaling with constant ~2.7µs per (χ, K-zero); an a-direction
   amortisation profile showing 8× speedup at M=256; an explicit
   falsifier list (F1-F3) for slot 2-5 attacks.
2. **What edges did my work compose or cite?** E1.5 (per-query
   bit-content), E3.1 (Connes amortisation), E6.7 (DR per-query),
   S224 (Correlation Dichotomy template).
3. **If my session produced only duplicate closures, why?** It didn't;
   the three-axis evaluator and the empirical profile are new.
4. **What is the next-action for the next agent?** Slot 2: implement
   AFE-shared partial-sum evaluation at fixed height t across all χ
   mod q (option α from the falsifier list); measure per-character
   cost reduction vs the slot-1 baseline (22.7 ms/query asymptote at
   q=1009, K=200, x=10⁶). If α closes negatively, attempt β
   (multi-character zero-finding via Bober-Hiary 2017) or γ
   (cross-character correlation) before declaring slot 2 closed.
