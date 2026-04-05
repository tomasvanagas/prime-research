# Prime Research: Computing p(n) Exactly Without Bruteforcing

## Goal
Find an O(polylog) algorithm to compute the nth prime p(n) exactly.
Target: p(10^100) in <1 second, 100% accurate.

## Status (April 2026)
- **649+ approaches tested** across 35 sessions, 192+ sub-agents
- **All known paths closed** but no proof that polylog is impossible
- **Problem is genuinely open** -- no unconditional lower bound beyond Omega(log x)
- **Proving impossibility faces Natural Proofs barrier** -- as hard as P != NP (S23)
- **pi(x) mod 2 is random-like in 20+ structural measures** (S35)
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
  CLOSED_PATHS.md      <-- 633+ tested approaches. SEARCH before proposing anything.
  OPEN_PROBLEMS.md     <-- The ONLY viable research directions. Start here.
  BEST_ALGORITHMS.md   <-- Working implementations with benchmarks.
  SESSION_INSIGHTS.md  <-- Detailed per-session findings (Sessions 12-31).
proven/
  barriers.md          <-- Mathematically proven impossibility results.
  complexity.md        <-- Upper/lower bounds, circuit complexity.
  information.md       <-- Entropy, Kolmogorov complexity, entanglement.
  quantum.md           <-- Why quantum computing doesn't help.
  *_barrier.md         <-- Individual barrier proofs (circuit size, uniformity, etc.)
novel/                       <-- ONLY genuinely original findings not in published literature.
  info_computation_gap.md    <-- Best original insight (delta(n) framing).
  failure_taxonomy.md        <-- Three failure modes: Circularity / Equivalence / Info Loss.
  approx_degree_prime.md     <-- adeg(chi_P) = N/2 (novel measurement).
  determinantal_complexity.md <-- dc(pi_N) connection to GapL.
  (+ other novel findings with evidence)
algorithms/                  <-- ONLY working, tested code.
literature/
  references.md              <-- All citations.
  state_of_art_2026.md       <-- Latest published results.
experiments/                 <-- ALL experiments organized by topic.
  analytic/ algebraic/ quantum/ ml/ information_theory/
  dynamical/ topological/ sieve/ circuit_complexity/ other/
  wildcard/ proposals/
  Each .py MUST have a companion <name>_results.md alongside it.
data/                        <-- Zeta zeros (200/300/500/1000).
archive/
  sessions/            <-- ALL session summaries/syntheses (session11-31+).
  ephemeral/           <-- Mode outputs (overwritten each cycle). DO NOT put these elsewhere.
    proposals_latest.md, proposals_session.md, critique_latest.md, wildcard_findings.md
  CLAUDE_OUTPUTS/      <-- Raw run logs (JSON + human-readable).
  visualizations/      <-- Generated plots and charts.
FOCUS_QUEUE.md               <-- Deep-dive tasks for focused sessions.
```

## Rules for AI Agents

### Workflow
1. **On startup:** check if `TODO.md` exists. If so, work through its items first.
2. **Read `status/OPEN_PROBLEMS.md`** for viable directions.
3. **Search `status/CLOSED_PATHS.md`** before proposing ANY approach.
4. Use sub-agents to save context window.
5. **Context management:** when context is filling up, write remaining work to
   `TODO.md` so the next session can continue. Then halt.
6. **Update this file** only for significant status changes (new best algorithm,
   major barrier proven, goal change). Do NOT add session-by-session details here --
   those go to `status/SESSION_INSIGHTS.md`.
7. DO NOT modify `run.sh` or `FOCUS_QUEUE.md`.
8. If you find the breakthrough, respond with exactly: I FOUND IT!!!

### File Placement (STRICT — read carefully)
9.  **Experiments** go to `experiments/<topic>/` with descriptive filenames.
10. **Every .py script MUST have a companion `<name>_results.md`** saved alongside
    it with the experiment's findings, verdict, and key numbers. Do NOT just capture
    results in your context — persist them to disk.
11. **Do NOT create multiple versions** of the same script (e.g., `foo.py`,
    `foo_v2.py`, `foo_quick.py`, `foo_small.py`). Refactor the original or use
    command-line arguments/flags. One script per experiment.
12. **Results format is `.md` only.** No `.txt` or `.json` for human-readable results.
    Raw data (if needed) goes in `.json`/`.csv` but must have a `.md` summary too.
13. **`novel/`** is ONLY for genuinely original findings not in published literature.
    - YES: new formulas, novel measurements, original theoretical connections
    - NO: session syntheses, barrier proofs, ephemeral mode outputs, literature surveys
14. **Session syntheses** (per-session summaries) go to `archive/sessions/`.
15. **Proven barriers** go to `proven/`, not `novel/`.
16. **Ephemeral mode outputs** (proposals, critiques, wildcard brainstorms) go to `archive/ephemeral/`.
17. **New literature** goes to `literature/`.
18. **Working algorithms** go to `algorithms/` with benchmarks.
19. Update `status/CLOSED_PATHS.md` when closing an approach.

### Cleanup (before finishing each session)
20. Delete any `__pycache__` directories you created.
21. Verify every `.py` you wrote has a companion `_results.md`.
22. Do NOT leave orphaned or duplicate scripts.
