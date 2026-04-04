# Prime Research: Computing p(n) Exactly Without Bruteforcing

## Goal
Find an O(polylog) algorithm to compute the nth prime p(n) exactly.
Target: p(10^100) in <1 second, 100% accurate.

## Status (April 2026)
- **520+ approaches tested** across 19 sessions, 135+ sub-agents
- **All known paths closed** but no proof that polylog is impossible
- **Problem is genuinely open** -- no unconditional lower bound beyond Omega(log x)
- Best exact: `algorithms/v10_c_accelerated.py` -- O(p(n)^{2/3}), p(10^9) in 0.175s
- Best approximate: R^{-1}(n) -- O(polylog), ~50% digits correct

## The Barrier (One Paragraph)
p(n) = SMOOTH(n) + RANDOM(n). The smooth part R^{-1}(n) is O(polylog) and gives
~50% of digits. The remaining ~50% encode oscillatory contributions of ~10^48
Riemann zeta zeros with GUE-random phases, information-theoretically incompressible.
Best known: O(x^{2/3}) combinatorial, O(x^{1/2+epsilon}) analytic.

## Where to Find Things

```
status/
  CLOSED_PATHS.md      <-- 500+ tested approaches. SEARCH before proposing anything.
  OPEN_PROBLEMS.md     <-- The ONLY viable research directions. Start here.
  BEST_ALGORITHMS.md   <-- Working implementations with benchmarks.
  SESSION_INSIGHTS.md  <-- Detailed per-session findings (Sessions 12-18).
proven/
  barriers.md          <-- Mathematically proven impossibility results.
  complexity.md        <-- Upper/lower bounds, circuit complexity.
  information.md       <-- Entropy, Kolmogorov complexity, entanglement.
  quantum.md           <-- Why quantum computing doesn't help.
novel/
  info_computation_gap.md    <-- Best original insight (delta(n) framing).
  failure_taxonomy.md        <-- Three failure modes: Circularity / Equivalence / Info Loss.
  session16_synthesis.md     <-- 15 intermediate families closed; three "pillars".
  session17_synthesis.md     <-- Exact communication rank formula.
  critique_latest.md         <-- Latest adversarial critique results (if exists).
  proposals_session.md       <-- Latest wildcard/propose session ideas (if exists).
algorithms/                  <-- ONLY working, tested code.
literature/
  references.md              <-- All citations.
  state_of_art_2026.md       <-- Latest published results.
experiments/                 <-- ALL experiments organized by topic.
  analytic/ algebraic/ quantum/ ml/ information_theory/
  dynamical/ topological/ sieve/ circuit_complexity/ other/
  wildcard/ proposals/
data/                        <-- Zeta zeros (200/300/500/1000).
archive/                     <-- Session logs and visualizations.
FOCUS_QUEUE.md               <-- Deep-dive tasks for focused sessions.
```

## Rules for AI Agents
1. **On startup:** check if `TODO.md` exists. If so, work through its items first.
2. **Read `status/OPEN_PROBLEMS.md`** for viable directions.
3. **Search `status/CLOSED_PATHS.md`** before proposing ANY approach.
4. Save experiments to `experiments/<topic>/` with descriptive filenames.
5. Update `status/CLOSED_PATHS.md` when closing an approach.
6. Novel findings go to `novel/` with evidence.
7. New literature goes to `literature/`.
8. If you beat the current best algorithm, save to `algorithms/` with benchmarks.
9. DO NOT modify `run.sh` or `FOCUS_QUEUE.md`.
10. Use sub-agents to save context window.
11. **Context management:** when context is filling up, write remaining work to
    `TODO.md` so the next session can continue. Then halt.
12. **Update this file** only for significant status changes (new best algorithm,
    major barrier proven, goal change). Do NOT add session-by-session details here --
    those go to `status/SESSION_INSIGHTS.md`.
13. If you find the breakthrough, respond with exactly: I FOUND IT!!!
