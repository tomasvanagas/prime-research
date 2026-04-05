# Prime Research: Computing p(n) Exactly Without Bruteforcing

## Goal
Find an O(polylog) algorithm to compute the nth prime p(n) exactly.
Target: p(10^100) in <1 second, 100% accurate.

## Status (April 2026)
- **525+ approaches tested** across 37 sessions, 197+ sub-agents
- **All known paths closed** but no proof that polylog is impossible
- **Problem is genuinely open** -- no unconditional lower bound beyond Omega(log x)
- **Proving impossibility faces Natural Proofs barrier** -- as hard as P != NP (S23)
- **pi(x) mod 2 is random-like in 20+ structural measures** (S35)
- Best exact: `algorithms/v10_c_accelerated.py` -- O(p(n)^{2/3}), p(10^9) in 0.175s
- Best approximate: R^{-1}(n) -- O(polylog), ~50% digits correct
- **TECH DEBT RESOLVED (S37):** All 351 experiment scripts now have _results.md companions.

## The Barrier (One Paragraph)
p(n) = SMOOTH(n) + RANDOM(n). The smooth part R^{-1}(n) is O(polylog) and gives
~50% of digits. The remaining ~50% encode oscillatory contributions of ~10^48
Riemann zeta zeros with GUE-random phases, information-theoretically incompressible.
Best known: O(x^{2/3}) combinatorial, O(x^{1/2+epsilon}) analytic.

## Where to Find Things

```
status/
  CLOSED_PATHS.md      <-- 525+ tested approaches. SEARCH before proposing anything.
  OPEN_PROBLEMS.md     <-- The ONLY viable research directions. Start here.
  BEST_ALGORITHMS.md   <-- Working implementations with benchmarks.
  SESSION_INSIGHTS.md  <-- Detailed per-session findings (Sessions 12-37).
proven/
  barriers.md          <-- Mathematically proven impossibility results.
  complexity.md        <-- Upper/lower bounds, circuit complexity.
  information.md       <-- Entropy, Kolmogorov complexity, entanglement.
  quantum.md           <-- Why quantum computing doesn't help.
  *_barrier.md         <-- Individual barrier proofs (circuit size, uniformity, etc.)
novel/                       <-- ONLY genuinely original findings not in published literature.
  pseudorandomness_of_pi.md  <-- STRONGEST FINDING: 21 measures show pi(x) mod 2 is random-like.
  info_computation_gap.md    <-- Best original insight (delta(n) framing).
  failure_taxonomy.md        <-- Three failure modes: Circularity / Equivalence / Info Loss.
  approx_degree_prime.md     <-- adeg(chi_P) = N/2 (novel measurement).
  bpsw_tc0_reduction.md      <-- BPSW correct => PRIMES in TC^0 (conditional).
  delta_spectrum.md          <-- delta(n) has 1/f^1.69 spectrum (no shortcut).
  determinantal_complexity.md <-- dc(pi_N) >= 2^{N/2-1}+2 (CLOSED).
  (+ entropy, ANF, GapL, algebraic geometry findings)
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
  sessions/            <-- ALL session summaries/syntheses (session11-37+).
  ephemeral/           <-- Mode outputs (overwritten each cycle). DO NOT put these elsewhere.
    proposals_latest.md, proposals_session.md, critique_latest.md, wildcard_findings.md
  CLAUDE_OUTPUTS/      <-- Raw run logs (JSON + human-readable).
  visualizations/      <-- Generated plots and charts.
FOCUS_QUEUE.md               <-- Deep-dive tasks for focused sessions.
```

## Rules for AI Agents

### Workflow (follow this EXACT order on startup)
1. **Read `TODO.md`.** If it has [HOUSEKEEPING] items, complete them first.
   If the only remaining housekeeping is file deletion/merging that you cannot
   safely do (e.g., duplicate cleanup requiring human review), skip it and
   note in TODO.md that it's blocked. Do NOT let blocked housekeeping prevent
   all research.
2. **Read `status/OPEN_PROBLEMS.md`** for viable research directions.
   As of Session 37, only circuit complexity of pi(x) remains genuinely open.
   If no viable experiment exists, check `literature/state_of_art_2026.md` for
   new publications that might open a direction. If nothing new, say so and stop.
   Do NOT re-run experiments on closed paths.
3. **Search `status/CLOSED_PATHS.md`** before proposing ANY approach.
4. **Read `novel/pseudorandomness_of_pi.md`** — this is the project's strongest
   finding (21 measures showing pi(x) mod 2 is random-like). Any new approach
   must explain how it circumvents this.
5. Use sub-agents to save context window.
6. **Context management:** when context is filling up, write remaining work to
   `TODO.md` so the next session can continue. Then halt.
7. **Update this file** only for significant status changes (new best algorithm,
   major barrier proven, goal change). Do NOT add session-by-session details here --
   those go to `status/SESSION_INSIGHTS.md`.
8. DO NOT modify `run.sh`. Update `FOCUS_QUEUE.md` only to mark tasks
   COMPLETED or add new tasks — do not delete completed task descriptions.
9. If you find the breakthrough, respond with exactly: I FOUND IT!!!

### What to do when all paths seem closed
The project is in a mature state. Most sessions will NOT produce breakthroughs.
Productive work includes:
- **Literature monitoring:** check for new 2026 publications on pi(x) complexity,
  TC^0 lower bounds, or zeta zero computation. Update `literature/state_of_art_2026.md`.
- **Engineering:** improve `algorithms/v10_c_accelerated.py` (Gourdon variant,
  segmented sieve, SIMD — see BEST_ALGORITHMS.md comparison table).
- **Theoretical sharpening:** tighten existing novel results (e.g., extend
  pseudorandomness measures to larger N, prove the N/2 threshold rigorously).
- **Do NOT** re-run completed experiments, propose approaches already in CLOSED_PATHS,
  or generate speculative proposals without concrete experiments to run.

### File Placement (STRICT — read carefully)
- **Experiments** go to `experiments/<topic>/` with descriptive filenames.
- **Every .py script MUST have a companion `<name>_results.md`** saved alongside
    it with the experiment's findings, verdict, and key numbers. Do NOT just capture
    results in your context — persist them to disk.
    - **Write the _results.md IMMEDIATELY after running the script.** Do not batch
      them up. Do not move on to the next experiment until the results file exists.
    - **If you cannot run a script** (missing deps, too slow), write a _results.md
      anyway noting: what it attempts, why it couldn't run, and your best assessment
      from reading the code.
    - **ENFORCEMENT:** Before ending ANY session, run this check:
      `find experiments/ -name "*.py" | while read f; do r="${f%.py}_results.md"; [ ! -f "$r" ] && echo "MISSING: $r"; done`
      If any are missing, write them before stopping.
- **Do NOT create multiple versions** of the same script (e.g., `foo.py`,
    `foo_v2.py`, `foo_quick.py`, `foo_small.py`). Refactor the original or use
    command-line arguments/flags. One script per experiment.
- **Results format is `.md` only.** No `.txt` or `.json` for human-readable results.
    Raw data (if needed) goes in `.json`/`.csv` but must have a `.md` summary too.
- **`novel/`** is ONLY for genuinely original findings not in published literature.
    - YES: new formulas, novel measurements, original theoretical connections
    - NO: session syntheses, barrier proofs, ephemeral mode outputs, literature surveys
- **Session syntheses** (per-session summaries) go to `archive/sessions/`.
- **Proven barriers** go to `proven/`, not `novel/`.
- **Ephemeral mode outputs** (proposals, critiques, wildcard brainstorms) go to `archive/ephemeral/`.
- **New literature** goes to `literature/`.
- **Working algorithms** go to `algorithms/` with benchmarks.
- Update `status/CLOSED_PATHS.md` when closing an approach.

### Status File Hygiene (STRICT — previous sessions violated this)
- **When an experiment closes a path:** add it to CLOSED_PATHS.md in the SAME session.
  Do not leave it for "later." The entry needs: approach name, verdict, failure mode
  (C/E/I), one-line key finding, session number.
- **When an open problem is resolved:** update OPEN_PROBLEMS.md in the SAME session.
  Mark it CLOSED with the evidence citation. Do not leave stale "open" problems.
- **When a completed experiment is labeled "pending" in ephemeral docs:** correct
  the label immediately. Stale "pending" labels create confusion for future sessions.
- **Every approximation formula or algorithm** must appear in BEST_ALGORITHMS.md,
  even if it's not exact. Separate sections for "Exact" and "Approximate" methods.
- **Novel findings that unify multiple session results** (e.g., "20+ measures show
  pseudorandomness") deserve their own document in `novel/`. Do not leave cross-session
  synthesis undone — it's arguably more valuable than any single experiment.

### Cleanup (before finishing each session)
- Delete any `__pycache__` directories you created.
- Verify every `.py` you wrote has a companion `_results.md` (run the find
  command from the File Placement section — zero tolerance for missing results files).
- Do NOT leave orphaned or duplicate scripts.
- Verify every experiment you ran has its verdict in `status/CLOSED_PATHS.md`.
- Verify no "pending" labels remain in ephemeral docs for completed work.
- Update `status/SESSION_INSIGHTS.md` with this session's key findings.
- Create `archive/sessions/sessionNN_<topic>.md` for this session.
