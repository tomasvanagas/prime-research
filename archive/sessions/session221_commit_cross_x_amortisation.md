# Session 221 — commit thread 5 step 2: cross-x amortisation, multipoint falsifier

**Date:** 2026-04-29
**Mode:** commit (Thread 5 / cross-x amortisation, slot 2 of 5)
**Prior:** S220 (slot 1) decoupled `K_zeros_setup` from
`K_per_x_evaluation`, established per-term Ei cost flooring at ~600 ns
by K ≈ 1600 (so `T_eval = Θ(K)` asymptotically), and structurally closed
the explicit-formula pillar of cross-x amortisation conditional on
Thread 3's Montgomery closure — pending one live falsifier:
**Schönhage 1990 / Odlyzko-Schönhage multipoint zeta evaluation.**
**Self-grade:** **B** — substantive empirical refinement: built the
Taylor-P multipoint evaluator, established empirically that Taylor-2
gives a constant-factor wall-clock speedup (~10× at K=200, M=16) but
NO asymptotic α reduction; closed the slot-2 falsifier under the
quantitative cluster-width bound `Δx ≤ x · (1/(2K))^{1/(P+1)} / γ_K`.
The explicit-formula pillar of Thread 5 is now structurally closed
(modulo Montgomery, inherited from Thread 3).

## Mission (slot-2, from `.commit_state` recommended_next_action)

> "Empirically measure per-x amortised cost when M = log x, log² x,
> x^{1/4}, x^{1/2} queries share K precomputed zeros. Build a multi-x
> evaluator that times the partial sum at M correlated x_i's. Plot
> per-x amortised cost vs M and confirm slot-1's forecast that
> amortised cost saturates at the per-term floor (~600 ns × K)
> independent of M. Live falsifier: Schönhage 1990 / Odlyzko-Schönhage
> multipoint zeta evaluation."

## What was built

`experiments/analytic/cross_x_amortisation/cross_x_batched_evaluator.py`:

- `batched_direct(xs, K, gammas)`: baseline M-batched evaluator. For
  each `x_i` in the cluster, evaluates `R(x_i) - 2 Σ_{j=1}^K Re
  R(x_i^{ρ_j})` using the slot-1 `R_at_rho` primitive.
- `batched_taylor(xs, K, gammas, x0, P)`: the slot-2 Schönhage falsifier
  — multipoint evaluation via Taylor-P expansion of `f_j(x) =
  R(x^{ρ_j})` around `x_0`. Setup computes `f_j(x_0), f_j'(x_0),
  f_j''(x_0)` for j = 1..K *exactly in closed form* via the identity

  ```
  d/dx Ei(ρ_j ln x / n) = exp(ρ_j ln x / n) / (x ln x)
  ```

  (slot-2 derivation; the `1/(x ln x)` factor is `n`-independent — see
  `f_j_value_and_derivs` docstring). Per-x evaluation reduces to K · P
  arithmetic ops (vectorised numpy), no `Ei` call.
- `cluster_around(x, M, mode)`: integer / log-spread / uniform-eta
  cluster generators.
- `policy_K_values` / `policy_M_values`: the slot-prescribed
  `{log²x, log³x, x^{1/4}, x^{1/2}}` × `{1, log x, log²x, x^{1/4}}`
  sweep.
- Q1 saturation curve, Q2 cross-x sharing check (direct), Q3 Taylor-P
  speedup measurement.

`experiments/analytic/cross_x_amortisation/cross_x_taylor_scaling.py`:
supplementary K-scaling check confirming the Taylor speedup formula is
K-independent at fixed M.

The evaluator reuses `R_at_rho`, `riemann_R`, `get_zeros`, `mobius_sieve`
from `experiments/analytic/connes_amortisation/connes_amortisation.py`
— the same primitive that backed S195 / S202 / S220.

## Empirical findings

### Q1 — saturation of T_per_x_amortised vs M

Direct M-batched evaluator at (x = 10⁶, K-policy = log³x = 191,
integer cluster), measured `T_per_x_amortised(K, M) = T_setup(K)/M +
T_eval(K, x)`:

| M           | T_eval/x (s) | T_amort_cached (s) | gain over M=1 |
|-------------|--------------|--------------------|---------------|
| 1           | 0.237        | 0.237              | 0.1%          |
| 14 (log x)  | 0.217        | 0.217              | 0.0%          |
| 32 (x^{1/4})| 0.209        | 0.209              | 0.0%          |
| 191 (log²x) | 0.206        | 0.206              | 0.0%          |

Setup amortisation gain is **0.0–0.1%** at all M because cached zero
loading is microseconds, vs ~200 ms per-x. At Hiary 2011 production
scale (`K^{17/13}` ops at 1 ns/op): for K=200 setup = 0.74 ms = 0.3% of
per-x cost. **T_per_x_amortised saturates at T_eval(K) immediately**,
exactly per slot-1's forecast.

Same pattern at x ∈ {10⁵, 10⁶} and K ∈ {log²x, log³x, x^{1/4},
x^{1/2}}. (Q1 confirmed.)

### Q3 — Schönhage / Odlyzko-Schönhage Taylor-P falsifier check

At (x = 10⁶, K = 200, P = 2, integer cluster):

| M  | T_direct (s) | T_taylor (s) | speedup | max err  |
|----|--------------|--------------|---------|----------|
| 4  | 0.752        | 0.317        | 2.37×   | 1.86e-10 |
| 8  | 1.779        | 0.316        | 5.63×   | 2.37e-09 |
| 16 | 3.462        | 0.311        | 11.15×  | 2.34e-08 |

Taylor-2 setup cost = 0.31 s (independent of M; one-time per cluster).
Taylor-2 eval per x = ~25 µs (K · P vectorised arithmetic ops).
Direct per x = ~200 ms.

**Asymptotic speedup formula**:

```
speedup(K, M) = M · a / (a' + M · b · P)  → a/(b · P)  as M → ∞
```

with `a` = per-zero per-x Ei cost (~1 ms), `a'` = per-zero Taylor setup
cost (~Ei + exp ≈ 1.5 ms), `b` = per-zero per-x arithmetic-op cost
(~25 ns / 200 zeros ≈ 0.125 µs). With P = 2: `a/(b·P) ≈ 4000×`. The
asymptotic speedup is **K-independent** (K cancels from numerator and
denominator).

### Cluster-width bound — the structural ceiling on M

Taylor-P truncation error per zero:

```
err_per_zero ≤ (γ_K · Δx / x)^{P+1} / (P+1)! ≤ 1/(2K)
⇒ Δx ≤ x · (1/(2K))^{1/(P+1)} / γ_K
```

With γ_K ≈ K and K = x^{1/2} (Galway 2004 explicit-formula
requirement under Thread 3):

| P        | Δx_max     | M_max (integer) |
|----------|------------|-----------------|
| 2        | x^{1/3}    | x^{1/3}         |
| 4        | x^{2/5}    | x^{2/5}         |
| 8        | x^{4/9}    | x^{4/9}         |
| O(log x) | ≈ x^{1/2}  | ≈ x^{1/2}       |

For P = O(1) and γ_K = K = x^{1/2}, per-x Taylor cost (cluster of
clusters tiling [1, x_max]):

```
Total_per_x_taylor = K · (P+1) · a' / Δx + K · P · b
                   = K^{2 + 1/(P+1)} / x · a' + K · b
                   = x^{1/(2(P+1))} · a' + x^{1/2} · b
                   = Θ̃(x^{1/2})  for P = O(1)
```

P = O(log x) gives Θ̃(x^{1/2} log x). **Same asymptotic order as direct
evaluator**. The Taylor / Schönhage approach gives only constant-factor
speedup.

## Structural conclusion (slot-2 final closure of explicit-formula pillar)

Combining slot 1 (`T_eval = Θ(K)` asymptotically, with per-zero floor
~600 ns) with slot 2:

1. **Direct M-batched per-x amortised cost saturates at T_eval(K)** for
   any M, because setup is microseconds while per-x is milliseconds.
2. **Taylor-P / Schönhage-Odlyzko multipoint evaluation is constant-
   factor only**, with cluster-width-bounded M ≤ x^{1/2}. Cluster-
   stitched per-x cost remains Θ̃(x^{1/2}).
3. **Combined with Thread 3 (S195+S196+S202): K* = Θ̃(x).** Therefore
   `T_per_x_amortised = Θ̃(x)` for any (M, cluster, multipoint scheme).

**The explicit-formula pillar of Thread 5 admits no batch-polylog
algorithm.** The closure inherits Thread 3's Montgomery pair-correlation
random-phase heuristic.

## Falsifiers (open for slot 5 to address)

The slot-2 closure is falsified only by:

1. **A multipoint primitive beating the cluster-width bound** —
   conformal-map continuation, FFT-based bandlimited reconstruction,
   etc. — with per-zero accuracy 1/(2K). Not known to exist for
   `R(x^{ρ_j})`.
2. **Per-zero accuracy is the wrong metric** — if per-zero truncation
   errors cancel systematically across `j`, the per-zero `1/(2K)`
   pointwise bound could relax. Phases of `R(x^{ρ_j})` are GUE-random
   under Montgomery (E2.10), so the variance of the partial-sum
   truncation error is controlled by RMS not pointwise — but the
   slot-1 / S195 / S202 measurements already used this implicitly
   (the in-distribution closure under Montgomery already has
   pair-correlation built in).
3. **An evaluator not expressible as `R(x) - Σ_j R(x^{ρ_j})`** with
   shared zeros. Galway 2004 / Platt 2017 are this form. Slot 3 (HKM /
   Meissel-Lehmer) is structurally distinct and outside scope.

## Slot 2 does NOT close

- **Slot 3 (HKM / Meissel-Lehmer cross-x state-sharing).** The
  combinatorial pillar uses `φ(x, a)` recursion, not partial sums of
  zero contributions; structurally distinct.
- **Slot 4 (Aggarwal binary search).** O(log x) sub-queries at
  binary-search-correlated x values. Slot-2 cluster-width bound applies
  but Aggarwal's analysis is more nuanced.
- **Slot 5 (theoretical wrap).** Rigorous lower bound matching the
  empirical per-x curve; rigorous α=1 proof.

## Edges composed / cited

- **E1.5** (information-theoretic per-query barrier): the Taylor
  cluster-width bound `Δx ≤ x · K^{-1-1/(P+1)}` quantitatively realises
  the bit-content barrier at the multipoint level — even with the most
  aggressive multipoint scheme the per-x asymptotic stays Θ̃(x^{1/2}).
- **E3.1** (Connes setup-cost dominance): Hiary `K^{17/13}` setup
  amortises to `K^{17/13}/M`, but per-x cost is the bottleneck for any M.
- **E6.6** (Aggarwal binary search): the slot-2 cluster-width and
  speedup formulae will feed slot 4's analysis directly.
- **S195 row 816 / S202 row** (Thread 3 closure under Montgomery):
  combined with slot-2 to give per-x amortised = Θ̃(x).
- **Schönhage 1990** (multipoint zeta evaluation): the live falsifier
  identified at slot 1; closed empirically as constant-factor only.
- **E2.10 / Montgomery pair correlation**: the GUE phases enable the
  in-distribution lift S195/S202 use; mentioned for completeness.

## Cross-domain ingredient (CROSS_DOMAIN_TECHNIQUES.md §8)

**Amortised algorithmics** (Tarjan 1985 / Iacono 2008 / Demaine-
Patrascu 2008): USED PARTIAL → USED I.

The slot-2 contribution is the *Schönhage-Odlyzko-Taylor multipoint
analysis* of the explicit-formula evaluator. Tarjan's framework
(separate amortisable from per-query work) is the structural device.
Schönhage 1990 / Odlyzko-Schönhage gives the multipoint primitive,
imported as Taylor-P expansion. Slot 2 imports both and reduces to a
structural negative result: cross-x amortisation gives constant-factor
speedup only — mode I closure.

## Why B-grade and not A or C

**Not A**: no batch-polylog algorithm achieved; the explicit-formula
pillar closes structurally as Θ̃(x^{1/2}) per-x asymptotic, conditional
on Thread 3's Montgomery closure. The closure is a refinement of
Thread 3, not a new positive result.

**Not C**: a new evaluator was built (Taylor-P multipoint, with
closed-form derivative computation via the `1/(x ln x)` n-independence
identity); the empirical Taylor speedup formula `a/(b·P)` was
quantified to ~13× at K=200 and predicted ~4000× asymptotically; the
structural cluster-width-and-stitch argument transformed the empirical
constant-factor speedup into an asymptotic equivalence. The argument
closes the live falsifier identified at slot 1 — substantive
refinement, not housekeeping.

## .commit_state changes

Slot 2 advances `sessions_used` 1 → 2, `session_history` += S221,
`recommended_next_action` updated to slot 3 (HKM / Meissel-Lehmer
cross-x state-sharing).

```
thread:cross_x_amortisation
sessions_used:2
status:ACTIVE
session_history:S220,S221
last_synthesis:archive/sessions/session221_commit_cross_x_amortisation.md
recommended_next_action: Slot 3 of Thread 5 — apply the cross-x
    amortisation decoupling to Meissel-Lehmer / HKM combinatorial-
    pillar π(x). Quantify the shareable-state fraction of HKM's
    T·S = x^{0.866} when querying π(x_1), ..., π(x_M) at evenly-
    spaced x_i. Slot-2 closes the explicit-formula pillar as
    constant-factor Schönhage / O(1)-asymptotic; slots 3-4 open the
    structurally distinct combinatorial and Aggarwal-binary-search
    pillars.
```

## Files modified by this session

- `experiments/analytic/cross_x_amortisation/cross_x_batched_evaluator.py` — new
- `experiments/analytic/cross_x_amortisation/cross_x_batched_evaluator_results.md` — new
- `experiments/analytic/cross_x_amortisation/cross_x_taylor_scaling.py` — new
- `experiments/analytic/cross_x_amortisation/cross_x_batched_direct.csv` — new
- `experiments/analytic/cross_x_amortisation/cross_x_batched_amortised.csv` — new
- `experiments/analytic/cross_x_amortisation/cross_x_batched_taylor.csv` — new
- `experiments/analytic/cross_x_amortisation/cross_x_taylor_scaling.csv` — new (supplementary)
- `.commit_state` — sessions_used 1 → 2, session_history += S221
- `RESEARCH_AGENDA.md` — Arc 7 slot 2 marked [x]
- `status/CLOSED_PATHS.md` — appended S221 row (§H.H1 slot 2)
- `status/SESSION_INSIGHTS.md` — S221 entry appended
- `CROSS_DOMAIN_TECHNIQUES.md` — §8 amortised-algorithmics row USED
  PARTIAL → USED I
- `archive/sessions/session221_commit_cross_x_amortisation.md` — this file
- `.run_state` → 223

## Session-end self-evaluation (CLAUDE.md 4-question)

1. **What did I produce that was not in the project before this
   session?** A multipoint Taylor-P evaluator with closed-form
   derivatives (via the n-independent `1/(x ln x)` factor identity);
   empirical speedup measurements at three (K, M) regimes; the
   cluster-width bound `Δx ≤ x · (1/(2K))^{1/(P+1)} / γ_K`; the
   structural argument that the cluster-stitched per-x Taylor cost is
   Θ̃(x^{1/2}) — same as direct, only constant-factor improvement.
   None of this existed in the project prior to S221. The slot-2
   closure of the explicit-formula pillar of Thread 5 is the new
   structural content.

2. **What edges did my work compose or cite?** E1.5 (per-query barrier
   — multipoint cluster-width bound is its quantitative realisation);
   E3.1 (Connes setup-cost dominance via Hiary `K^{17/13}` — amortises
   away but irrelevant); E6.6 (Aggarwal binary search — slot 4 will
   directly use the cluster-width and speedup formulae); S195 row 816
   / S202 row (Thread 3 K* = Θ̃(x) closure under Montgomery — combined
   with slot-2 to yield per-x = Θ̃(x)); E2.10 / Montgomery pair
   correlation (GUE phases enabling in-distribution lift); Schönhage
   1990 / Odlyzko-Schönhage (multipoint zeta evaluation, the live
   slot-1 falsifier — closed as constant-factor only); Tarjan 1985 /
   Iacono 2008 / Demaine-Patrascu 2008 (amortised-algorithmics
   framework).

3. **If my session produced only duplicate closures, why?** It didn't.
   The Taylor-P multipoint evaluator with closed-form derivatives is
   new code; the cluster-width bound combined with cluster-stitching is
   a new structural argument; the slot-2 closure of the explicit-
   formula pillar is the explicit follow-up to S202 WRAP's "open
   falsifier" listing — exactly what the slot-2 commitment promised.

4. **What is the next-action for the next agent?** Slot 3 of Thread 5
   (3 commit slots remaining). Apply the cross-x amortisation
   decoupling to Meissel-Lehmer / HKM combinatorial-pillar π(x). Build
   an instrumented HKM / Meissel-Lehmer evaluator that exposes which
   intermediate state can be reused across M batched x_i. Measure
   shareable fraction of T·S = x^{0.866}. The combinatorial pillar is
   structurally distinct from the explicit-formula pillar (slot-2's
   closure does NOT carry over): HKM's `φ(x, a)` recursion has a
   different cross-x behaviour. A-grade probability bumps slightly
   because the HKM tables are larger and reuse may be more substantive
   than zero precomputation.
