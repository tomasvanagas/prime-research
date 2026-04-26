# Prime Research: Computing p(n) Exactly Without Bruteforcing

## Goal
Find an O(polylog) algorithm to compute the nth prime p(n) exactly.
Target: p(10^100) in <1 second, 100% accurate.

## Status (April 2026)
- **533+ approaches tested** across 41 sessions, 197+ sub-agents
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
EDGES.md                <-- 66 real mathematical edges across 10 sections,
                            tagged with IDs (E1.x, ...). Cite by ID when
                            classifying closures or building chains. READ
                            EARLY in every session. Extend with new edges
                            discovered during the session.
status/
  CLOSED_PATHS.md      <-- 720+ tested approaches. SEARCH before proposing anything.
  OPEN_PROBLEMS.md     <-- The ONLY viable research directions. Start here.
  BEST_ALGORITHMS.md   <-- Working implementations with benchmarks.
  SESSION_INSIGHTS.md  <-- Detailed per-session findings (Sessions 12-66+).
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
  constructions/             <-- NEW (S60): purpose-built mathematical objects
                                 (custom circuits, rings, transforms, representations,
                                 algorithms). Each subdir has <name>.py +
                                 <name>_results.md, optionally definition.md.
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
2. **Read `EDGES.md`.** This is the catalogue of every real mathematical edge
   surfaced across 60+ sessions, grouped into surviving chains and tagged with
   IDs (E1.x, E2.x, E3.x, E4.x, E5.x, E6.x, E7.x). Cite edge IDs by name when
   classifying closures, building chains, or arguing why an approach is
   constrained. The §0 frame (T1-T6) and §10 single-target objectives tell you
   what the chain has to fit through and where uncovered attack surface
   actually lives. Do NOT propose a chain whose endpoint reduces to a known
   negative-shape edge (E7.x).
3. **Read `status/OPEN_PROBLEMS.md`** for viable research directions.
   As of Session 66, only circuit complexity of pi(x) remains genuinely open,
   and Chain E (the AKS-family route to TC^0) is "computationally cornered"
   per E7.10 — modulus / coefficient-ring / gcd-strengthening twists are
   orthogonal to depth. If no viable experiment exists, check
   `literature/state_of_art_2026.md` for new publications that might open a
   direction. If nothing new, say so and stop. Do NOT re-run experiments on
   closed paths.
4. **Search `status/CLOSED_PATHS.md`** before proposing ANY approach.
5. **Read `novel/pseudorandomness_of_pi.md`** — this is the project's strongest
   finding (22+ measures showing pi(x) mod 2 is random-like). Any new approach
   must explain how it circumvents this.
6. Use sub-agents to save context window.
7. **Context management:** when context is filling up, write remaining work to
   `TODO.md` so the next session can continue. Then halt.
8. **Update this file** only for significant status changes (new best algorithm,
   major barrier proven, goal change). Do NOT add session-by-session details here --
   those go to `status/SESSION_INSIGHTS.md`.
9. DO NOT modify `run.sh`. Update `FOCUS_QUEUE.md` only to mark tasks
   COMPLETED or add new tasks — do not delete completed task descriptions.
10. **When closing or building, cite EDGES.md edge IDs in your CLOSED_PATHS
    entry and session synthesis.** A closure that reduces to E5.3 should say
    so by name; a chain that exploits E1.6 (A⊕C₃ bisection) or E2.1 (MPS
    bond-dim identity) should reference the source edge. This keeps EDGES.md
    cross-linked with closures and prevents drift.
11. **When you discover a new edge** — a verified mathematical fact, identity,
    or measurable structural deviation not already in EDGES.md — add it
    inline in the appropriate section (§1 information, §2 algebraic, §3
    analytic, §4 sieve, §5 conditional/TC⁰, §6 computational, §7 negative-
    shape) with a fresh ID, EVS rating (H/M/L/shape), and a one-line
    "why this is an edge" justification. Update the closing footer to record
    the addition.
12. If you find the breakthrough, respond with exactly: I FOUND IT!!!

### Construction is encouraged, not just measurement (S60)

**Strategic note:** the project is heavy on measurement (test, batter,
PSLQ, close) and light on construction. After 58+ sessions, most work has
been "measure existing object X and report whether it deviates from random."
Genuine construction work — *building* a new circuit, algebraic object,
transform, representation, or algorithm — is welcome and may now be the
highest-leverage path forward.

**You have explicit permission to invent new mathematical objects** if
you see how they could help. Examples of what counts as construction:

- **New circuits:** prototype a TC^0 / NC^1 / AC^0 circuit for a specific
  primitive (e.g., the three FOCUS-1 sub-attacks: Bernstein 2003 smaller-r
  AKS, non-cyclotomic ring AKS, Healy-Viola Frobenius transplant — none
  of which has yet been built).
- **New algebraic objects:** define a custom ring, polynomial system,
  module, or representation tailored to pi(x). E.g., a non-cyclotomic
  quotient `Z_n[x]/f(x)` with a designed splitting structure; a custom
  L-function whose zero set encodes pi(x) less wastefully than zeta;
  a custom group action whose orbits separate primes from composites.
- **New transforms:** invent an integral / discrete transform purpose-built
  for the prime indicator (not just DCT/Fourier/wavelet). The bar is that
  it must have a concrete forward+inverse computation, not be hand-wavy.
- **New representations:** alternative encodings of pi(x) — tensor
  networks, spectral triples, sheaves, simplicial complexes — provided
  you build an evaluator that produces a number, not just a definition.
- **New algorithms:** combinations of existing primitives in non-obvious
  ways (Newton + bisection + Miller-Rabin walks etc. is the model).
- **New lower-bound techniques:** apply Brandt-MKtP-style diagonalisation,
  pebbling arguments on non-DAG models, or invent your own.

**Construction rules (so it stays disciplined):**

1. **Each construction MUST produce code that runs**, even if the object
   is "theoretical." A circuit prototype: write a Python simulator that
   evaluates it on small inputs. An algebraic object: write its operations
   and test associativity / cardinality on small parameters. A transform:
   compute it on at least one concrete input and verify the inverse.
2. **Each construction MUST be filed under
   `experiments/constructions/<descriptive_name>/`** with a
   `<name>.py`, `<name>_results.md`, and (for algebraic objects only) a
   short `definition.md` giving the object's signature and its expected
   relationship to pi(x).
3. **Failed constructions are valuable.** A built-and-tried-and-broken
   object generates a CLOSED_PATHS entry of much higher information
   content than yet another statistical battery. File the failure mode
   (C/E/I or "construction-incoherent" if the object turned out
   ill-defined) with the same discipline as any other experiment.
4. **Genuinely novel objects that turn out to work in some non-trivial
   way** (even partially) get a `novel/<name>.md` entry describing the
   object, its construction, what it does, and what it doesn't do. The
   bar for `novel/` placement remains "not in published literature."
5. **You still cannot rebuild objects that are explicit re-runs of
   closed paths.** E.g., do not "construct" yet another DCT-of-delta
   sparsifier (E6.3 caveat closed it). But you CAN construct a novel
   non-DCT basis if you have a concrete reason to think it might
   compress where DCT didn't.

**The bar for speculation:** if your construction is more than ~200
lines of prose without any code or concrete object, you are speculating.
Either turn it into code or save it as an `archive/ephemeral/` brainstorm
and move on.

### What to do when all paths seem closed
The project is in a mature state. Most sessions will NOT produce breakthroughs.
Productive work includes (roughly in order of expected leverage):
- **Construction work** as described above — currently the most underweighted
  category and the most likely path to a result. The three FOCUS-1
  sub-attacks in `TODO.md` are explicit construction tasks none of which
  has been attempted; FOCUS-3 (Brandt) is another.
- **Theoretical sharpening:** tighten existing novel results (e.g., extend
  pseudorandomness measures to larger N, prove the N/2 threshold rigorously,
  formalise the MPS bond-dim theorem in Lean).
- **Engineering:** improve `algorithms/v10_c_accelerated.py` (Gourdon variant,
  segmented sieve, SIMD — see BEST_ALGORITHMS.md comparison table).
- **Literature monitoring:** check for new 2026 publications on pi(x) complexity,
  TC^0 lower bounds, or zeta zero computation. Update `literature/state_of_art_2026.md`.
- **Do NOT** re-run completed experiments, OR propose approaches already in
  CLOSED_PATHS without a *new* angle on them. (The "no speculative proposals"
  rule from earlier sessions is **relaxed**: speculation that lands as
  buildable code in `experiments/constructions/` is now welcome.)

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
    - YES: new formulas, novel measurements, original theoretical connections,
      **new mathematical objects you constructed that turned out to do something
      non-trivial** (S60 expansion).
    - NO: session syntheses, barrier proofs, ephemeral mode outputs, literature surveys.
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
