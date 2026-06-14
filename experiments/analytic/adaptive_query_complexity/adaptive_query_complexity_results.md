# P10 — adaptive query complexity of computing p(n): CLOSED-NEGATIVE with mechanism (S532)

**Target:** OPEN_POSITIVE_TARGETS **P10** — "π(x₁) computed first, then x₂
chosen adaptively as a function of π(x₁). … other adaptive strategies
(adversarial, information-theoretic-optimal) might amortise better [than
Aggarwal binary search]." A-grade target: an adaptive strategy *strictly
better* than Aggarwal on a concrete benchmark.

**Verdict: CLOSED-NEGATIVE.** Adaptive π-querying gives **no asymptotic
benefit** — neither in query *cost* nor (in the zone that matters) query
*count* — over the **non-adaptive** "analytic-predict + local-sieve" recipe.
The reason is the project's own **SMOOTH + RANDOM** decomposition: the
smooth far zone is where adaptivity helps and is already *free*; the random
√x hard zone is where the work lives and admits *no* adaptive shortcut.

Code: `adaptive_query_complexity.py` (`--selftest`, `--measure`,
`--window-demo`). Log: `run_measure.log`. 9 selftest groups (R≈π within
3√x; Newton round-trip + iter bound; exact p(n)/π boundary values; landing
window < √x; both searches recover p(n); `_bracket` repair; within-window
residual is random; R monotone).

This **builds on and sharpens** `guess_comparison_oracle/` (S491), which
established the *geography* (far zone polylog-decidable, hard-zone width
~√p(n) with median 0.13/max 0.53, R-guess comparison = GUE fair coin) but
never made the *query-complexity* statement nor measured the Newton-iter /
window-content laws. It also makes precise the CLOSED_PATHS row "Binary
search + local sieve — problem relocated, not solved" (the *cost* face for
a single oracle) as a full count-vs-cost dichotomy.

---

## The algorithm being defended (non-adaptive, Õ(√x))

To compute p(n):
1. **Predict** x̂ = R⁻¹(n) by Newton on the polylog Riemann R (smooth far
   zone). Lands within Θ(√x) of p(n).
2. **Sieve** the window [x̂ − W, x̂ + W], W = Θ(√x), once; read off p(n).

Total Õ(√x), **zero adaptive π-queries**. P10 asks whether *adaptive*
π-querying can beat this. It cannot — measured below.

## Measurements (`--measure --xmax 2·10⁷`, π up to 1.27M, p(n)≤1.9·10⁷)

**M1 — smooth zone is free & non-adaptive (Newton iters O(log log x)).**
R⁻¹(n) by Newton from x₀=n·ln n converges in **2–3 iterations, flat across
22 scales** (n = 199 … 1.2M, ~4 orders of magnitude in x). min/mean/max =
2/2.95/3. The leading ~half of p(n)'s digits are non-adaptively, polylog
computable; no querying involved.

**M2 — hard zone has width Θ(√x) (the part you must read).** Over 240
log-spaced n in 8 scale-bins:
- RMS|p(n) − R⁻¹(n)| ~ **x^0.488** (predict 0.50).
- err/√(p(n)): median **0.145**, mean 0.178, max **0.721** (guess-oracle
  S491 over 30 pts: 0.13 / 0.53; denser sampling here catches a fatter tail
  — sized at ~0.72 at this scale, RH worst case √x·log x).
- trend of err/√x vs x: slope **−0.007 ≈ 0** → the ratio is **bounded, no
  upward drift** ⇒ window = Θ(√x). (F2 not triggered.)

**M3 — the hard zone carries no smooth signal to interpolate.** Within the
±√x window, fit π(x) by polynomials and take the RMS residual:

| n | p(n) | window count | √count | rms deg1 | rms deg3 | rms/√count |
|---|---|---|---|---|---|---|
| 1 586 | 13 367 | 23 | 4.8 | 0.70 | 0.65 | 0.135 |
| 12 613 | 135 409 | 69 | 8.3 | 1.40 | 0.96 | 0.115 |
| 100 285 | 1 303 597 | 161 | 12.7 | 2.31 | 1.48 | 0.117 |
| 797 355 | 12 151 679 | 422 | 20.5 | 4.55 | 2.45 | 0.119 |

The residual sits at the **√count ~ x^{1/4} Poisson scale** (rms/√count flat
~0.12) and **does not drop** going deg1→deg3 (~0.6×, not orders) — a
higher-degree smooth model captures essentially nothing. No interpolable
signal. (F3 not triggered.)

**M3 corollary — inside the random window interpolation ≈ binary (count).**
Bracketing p(n) inside the ±√x window with an *exact* π oracle:

| n | p(n) | window=2√x | binary_q | interp_q |
|---|---|---|---|---|
| 1 586 | 13 367 | 230 | 8 | 9 |
| 12 613 | 135 409 | 734 | 10 | 9 |
| 100 285 | 1 303 597 | 2 282 | 12 | 13 |
| 797 355 | 12 151 679 | 6 970 | 13 | 13 |

`interp_q ≈ binary_q ≈ ½·log₂(x)` — interpolation gets **no traction**
inside the random window (it would need a smooth value-vs-position signal;
M3 shows there is none). Interpolation's textbook loglog advantage exists
only in the *smooth* zone — which M1 already does for free.

## The closure (count vs cost)

- **Smooth far zone:** adaptivity (interpolation/Newton) gives loglog
  convergence — but this zone is already *free* (polylog R, M1). No
  π-oracle queries needed; R⁻¹(n) lands directly in the √x window.
- **Random √x hard zone:** width Θ(√x) (M2), content is √x/ln x near-Poisson
  primes with no smooth structure (M3). To pin p(n) you must *read* this
  window. In the **count** metric interpolation does not beat binary here
  (M3 corollary) — ½log₂x probes either way. In the **cost** metric each
  hard-zone π-eval costs ~√x (best-known), so ½log₂x of them = √x·polylog =
  the cost of **one sieve of the window**, which a single non-adaptive sweep
  performs while answering *all* positions at once.

⇒ The **cost-optimal** p(n) algorithm uses **zero adaptive π-queries**
(predict + sieve, Õ(√x)). Every adaptive π-query strategy (Aggarwal binary
search, interpolation search, …) is at best equal and only "relocates" the
√x into ½log₂x hard-zone π-evals. **Adaptivity cannot beat non-adaptive;
the √x random width is irreducible.** The SMOOTH+RANDOM split is exactly the
adaptivity boundary.

This is consistent with — and the query-complexity face of — the
information barrier (S511 √x info floor; S518 the analytic witness is √x)
and the guess-comparison geography (S491).

## Scope / honesty

- **Negative, adjacent.** This closes an *adjacent* target (a query model on
  top of a π oracle); it does **not** advance the goal (polylog exact π
  stays blocked). It is a precisely-closed question with a mechanism, not a
  novelty claim.
- The "no adaptive strategy beats it" argument is a **mechanism backed by
  measurement on natural strategies** (binary, interpolation), not a formal
  lower bound over *all* conceivable adaptive query strategies. A formal
  universal adaptive lower bound would need the information-theoretic /
  #P-hardness machinery (PROGRAM item 3). What is rigorous here: R⁻¹(n)+sieve
  is Õ(√x) and non-adaptive; the measured window is Θ(√x); within-window
  structure shows no interpolable signal.
- Ground truth is an honest Eratosthenes sieve; p(n)/π values exact.

## Falsifiers (pre-stated)

- **F1** Newton iters growing like log x (not ~const / loglog) → smooth zone
  not polylog-cheap. **Not seen** (2–3, flat).
- **F2** |p(n) − R⁻¹(n)| growing slower than √x (e.g. polylog) → hard zone
  shrinks, p(n) polylog-sieveable. **Not seen** (exp 0.488, ratio trend
  −0.007). (Faster than √x·polylog would violate RH.)
- **F3** within-window π residual dropping well below √count for some
  low-degree model → an interpolable smooth signal. **Not seen** (flat
  rms/√count ~0.12, no deg1→deg3 drop).
- **F4** an adaptive strategy resolving p(n) in o(√x) total *cost* on the
  benchmark. **Not seen** (every probe is a hard-zone √x eval; window read
  dominates).

Reproduce: `python3 adaptive_query_complexity.py --selftest` then
`--measure --xmax 20000000` (≈ a few seconds; RAM-light, ~30 MB sieve).
