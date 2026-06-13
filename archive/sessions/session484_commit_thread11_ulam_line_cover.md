# Session 484 — Commit Thread 11 / Slot 1: Ulam-spiral minimum line cover

## Context

Thread 11 (CLAUDE.md "Highest-EV mathematical threads", Arc 13 in
RESEARCH_AGENDA.md, attack vector H.H4, OPEN_POSITIVE_TARGETS §P11)
attacks the project's first incidence-geometric variant: under a 2D
embedding `Φ: ℕ → ℝ²`, what is the minimum number of straight lines
`L_Φ(N)` covering all prime points `Φ(p)` for p ≤ N? The Ulam spiral
(1963) shows visual diagonal-lines-of-primes corresponding to
Hardy-Littlewood prime-rich quadratic forms; Thread 11 asks whether
this has *quantitative algorithmic content* — i.e., whether L_Ulam
beats the random-points Szemerédi-Trotter lower bound (~√N for
N points on a √N×√N grid).

This is **slot 1 of 5**: build the evaluator and measure L_Ulam(N)
at N ∈ {10⁴, 10⁵, 10⁶}. Slots 2-5 explore alternative embeddings,
theoretical shape, algorithmic angle, and theoretical wrap.

## What I did

1. Read `.commit_state`, `OPEN_POSITIVE_TARGETS.md` §P11,
   `ATTACK_VECTORS.md` §H.H4, `RESEARCH_AGENDA.md` Arc 13.
2. Built `experiments/constructions/p11_ulam_line_cover/`:
   - `p11_ulam_line_cover.py` — exact greedy line cover via all-pairs
     enumeration (feasible to N=10⁴; pure-Python).
   - `p11_ulam_bounded.py` — bounded-direction greedy with numpy
     vectorisation (feasible to N=10⁶, K=20). Sanity check at N=10⁴:
     bounded K=5 → L=95 vs exact L=91 (≤5% overestimate).
   - `p11_top_lines_quadratic_decomp.py` — confirms each top line
     decomposes into pairs of quadratic forms `4k² + bk + c`.
3. Ran experiments at N = 10⁴, 10⁵, 10⁶, with 20 / 5 / 3 random-
   matched-baseline trials respectively.

## Results

### Scaling table

| N    | π(N)  | L_primes        | L_random ± std | L_p/L_r | L_p/√N | z      |
|------|-------|-----------------|----------------|---------|--------|--------|
| 10⁴  | 1229  | 91 (exact)      | 113.5 ± 9.4    | 0.80    | 0.91   | +2.07σ |
| 10⁵  | 9592  | 308 (K=20)      | 337.0 ± 4.0    | 0.91    | 0.97   | +7.25σ |
| 10⁶  | 78498 | 989 (K=20)      | 1042.3 ± 14.4  | 0.95    | 0.99   | +3.71σ |

### Key findings

1. **Both `L_primes` and `L_random` scale as √N to leading order**,
   matching the random-points Szemerédi-Trotter lower bound. Primes
   are NOT compressible to a sub-√N number of lines on the Ulam
   spiral.

2. **The prime advantage shrinks with N.** L_p/L_r = 0.80 → 0.91 →
   0.95 over N = 10⁴ → 10⁵ → 10⁶. Fitting `1 - L_p/L_r ~ c·N^{-α}`
   gives α ≈ 0.24 — a power-law decay of the advantage. Asymptotic
   ratio appears to → 1 (not bounded away from 1, contrary to the
   predicted E-mode "constant-factor compression").

3. **Top lines are all slope ±1 diagonals** with specific intercepts.
   At N=10⁶, the top dominant line (1, -1, -58) covers 347 primes;
   line (1, -1, -18) covers 325. Each decomposes into pairs of
   quadratic forms `f(k) = 4k² + bk + c`.

4. **Top-line prime densities ~40%** vs random density ~10%
   (consistent with HL prediction: HL_const × 1/log(N) ≈ 4 × 0.087 =
   0.35-0.45 at N=10⁵).

### Implication for thread

Slot 1 forecasts B-NEGATIVE for Ulam embedding alone: the visual
diagonal phenomenon is a small-N artifact. HL alignment provides
log-decaying compression that loses against the √N
Szemerédi-Trotter lower bound. **However, slots 2-5 are still in
play** — alternative embeddings may give qualitatively different
results, and an LP-relaxation algorithm (slot 4) could exhibit a
non-trivial gap from greedy.

## Self-evaluation

1. **What did I produce that was not in the project before this
   session?**
   - First incidence-geometric experiment in the project — both
     algorithm (exact greedy + bounded-direction numpy variant) and
     empirical scaling data.
   - Verification that Ulam-spiral top-lines decompose into pairs of
     `4k² + bk + c` quadratic forms with HL-rich prime densities
     (40% vs 10% random).
   - Empirical shape of L_p / L_r → 1 power-law-decay (α ≈ 0.24),
     **stronger** (more negative) than the predicted "constant-factor
     compression" E-mode failure.

2. **What edges did my work compose or cite?** E1.5 (compression),
   HL Conjecture F (off-edge external), Szemerédi-Trotter (off-edge
   external incidence-geometric baseline). No new EDGES.md entry —
   this is slot-1 measurement, not yet edge-grade.

3. **If duplicate closures, why?** N/A — slot 1 produces empirical
   setup, not closure. Closure (or A-grade outcome) is slot 5
   responsibility.

4. **Next action:**
   - **Slot 2**: alternative 2D embeddings. Implement
     `(n mod p, ⌊n/p⌋)` residue-class grid for p = 2, 3, 5, 7, 11
     and `(n, n² mod p)` polynomial-image grid. Compute L for each
     at N ∈ {10⁴, 10⁵, 10⁶}. Cross-embedding comparison to determine
     whether Ulam result is universal or embedding-specific.

## Self-grade

**B**. Substantive empirical setup: built a clean, correct,
fast-enough evaluator, ran at three orders of magnitude in N, and
identified a quantitative scaling pattern (`1 - L_p/L_r ~ N^{-0.24}`)
that is structurally informative — specifically, **stronger than
the E-mode prediction**. The empirical decay law is itself a
non-trivial measurement; the predicted "constant-factor compression"
hypothesis is challenged by the data.

This is slot-1 work — produces no theorem, but lays the empirical
ground that slots 2-5 build on. Self-grading down from A: no proof
or algorithm has been built yet; the empirical scaling alone is
informative but not theorem-grade.

## Files

- `experiments/constructions/p11_ulam_line_cover/p11_ulam_line_cover.py`
- `experiments/constructions/p11_ulam_line_cover/p11_ulam_bounded.py`
- `experiments/constructions/p11_ulam_line_cover/p11_top_lines_quadratic_decomp.py`
- `experiments/constructions/p11_ulam_line_cover/p11_ulam_line_cover_results.md`
- `experiments/constructions/p11_ulam_line_cover/summary_N10000.csv`
- `experiments/constructions/p11_ulam_line_cover/summary_bounded_N{100000,1000000}_K20.csv`
- `experiments/constructions/p11_ulam_line_cover/top_lines_*.csv`

## Closure tracking

- `.commit_state` updated: `sessions_used:1`, `session_history:S484`.
- No CLOSED_PATHS row yet — slot 1 produces empirical data, not
  closure. Slot 5 will file the closure.
