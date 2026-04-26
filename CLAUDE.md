# Prime Research: Computing p(n) Exactly Without Bruteforcing

## Goal
Find an O(polylog) algorithm to compute the nth prime p(n) exactly.
Target: p(10^100) in <1 second, 100% accurate.

## Status (April 2026, after 67+ sessions)

The project is in a **post-closure mature state**. Every known attack family
on the only open problem (circuit complexity of pi(x)) has been
systematically closed:

- **AKS family closed** (E7.10, S61/S64/S66) — modulus / ring / gcd twists
  are orthogonal to depth.
- **Brandt MKtP / diagonalisation-via-meta-complexity family closed**
  (E5.8, S51) — technique is structurally welded to MKtP itself.
- **Convergence-acceleration / variance-reduction family closed** (E7.11) —
  13+ transforms tested, all strictly worse than partial sum.
- **Sieve route asymptotically tight** (E6.7 + E7.6) within x^{0.034}.
- **35+ pseudorandomness measures** confirm pi(x) mod 2 is random-like.
- **730+ approaches tested** in CLOSED_PATHS.md.

## What this means for the workflow

The project's needs have shifted. **Disciplined-duplicate-closure has hit
diminishing returns.** Most "new" angles map to existing CLOSED_PATHS lines
within 30 minutes of work; the marginal information per session is dropping.

Going forward, sessions are evaluated on **novelty production**, not on
discipline alone. Every session must answer: *what mathematical object,
identity, structure, or proof did this session produce that did not exist
in the project before?* If the answer is "another duplicate closure" —
that is a failed session, even if every CLAUDE.md rule was followed.

Five things genuinely move the project forward at this stage. Sessions
should pick exactly one of these as their primary target:

1. **Composition.** Build a mathematical object that exploits ≥ 2 existing
   edges from EDGES.md simultaneously (see `NOVELTY_CHALLENGES.md` §1).
2. **Frame-shift.** Attack a research-grade open question that is NOT
   "polylog π(x)" but is in the same neighbourhood (see
   `NOVELTY_CHALLENGES.md` §2).
3. **Lean formalisation.** Pick an EDGES.md edge or a `novel/` result
   and produce a Lean 4 proof (see `NOVELTY_CHALLENGES.md` §3).
4. **Synthesis.** Write a publishable-quality unification across multiple
   `novel/` results (see `RESEARCH_AGENDA.md` arc 1).
5. **Multi-session arc work.** Continue an arc from `RESEARCH_AGENDA.md`.

Measurement-only sessions (run a battery, close a path, file a row) are
no longer the default. They are now a fallback option, used only when
none of the above five fits the session's available time.

## The Barrier (One Paragraph)

p(n) = SMOOTH(n) + RANDOM(n). The smooth part R^{-1}(n) is O(polylog) and
gives ~50% of digits. The remaining ~50% encode oscillatory contributions
of ~10^48 Riemann zeta zeros with GUE-random phases, information-
theoretically incompressible. Best known: O(x^{2/3}) combinatorial,
O(x^{1/2+epsilon}) analytic. The polylog gap is real but every standard
attack family is now closed.

## Where to Find Things

```
NOVELTY_CHALLENGES.md       <-- The active research targets. READ FIRST.
                                Composition challenges, frame-shift questions,
                                Lean formalisation queue.
RESEARCH_AGENDA.md          <-- Long-horizon multi-session arcs that survive
                                across sessions. Each arc has a state and a
                                next-action. Update on every visit.
TODO.md                     <-- Session-level tasks. Now thin: housekeeping
                                + recurring work + pointer to the above.
EDGES.md                    <-- 68 real mathematical edges across 10 sections
                                tagged with IDs. Cite by ID.
status/
  CLOSED_PATHS.md           <-- 730+ tested approaches.
  OPEN_PROBLEMS.md          <-- Maintained but mostly closed.
  BEST_ALGORITHMS.md
  SESSION_INSIGHTS.md
proven/                     <-- Mathematically proven results.
novel/                      <-- ONLY genuinely original findings.
algorithms/                 <-- Working tested code.
literature/
experiments/
  analytic/ algebraic/ quantum/ ml/ information_theory/
  dynamical/ topological/ sieve/ circuit_complexity/ other/
  wildcard/ proposals/
  constructions/            <-- Purpose-built mathematical objects.
                                <name>.py + <name>_results.md (+ definition.md
                                for algebraic objects).
  formalisations/           <-- NEW: Lean 4 proofs of edges and novel results.
                                <name>.lean + <name>_notes.md.
data/
archive/
  sessions/
  ephemeral/
  CLAUDE_OUTPUTS/
  visualizations/
```

## Rules for AI Agents

### Workflow (the ordering matters)

1. **Read `NOVELTY_CHALLENGES.md`.** Pick one challenge that fits the
   session's available time. This is your primary target.
2. **Read `RESEARCH_AGENDA.md`** to see if any in-flight arc is waiting
   for the kind of work you're about to do. If so, prefer continuing that
   arc over starting fresh.
3. **Read `EDGES.md`** to ground your work in existing mathematical
   facts. You will cite edge IDs.
4. **Read `status/OPEN_PROBLEMS.md`** if your target involves the only
   open problem (circuit complexity of pi(x)).
5. **Read `TODO.md`** for housekeeping items.
6. **Search `status/CLOSED_PATHS.md`** before claiming any result is
   novel. The "novel" bar is HIGH at this stage of the project.
7. **Read `novel/pseudorandomness_of_pi.md`** if your target involves
   pi(x) mod 2 structural claims. Any new approach to this function must
   either circumvent the 35+ measures or explain why they don't apply.

### The Novelty Bar (THIS IS THE CORE CRITERION)

A session is **successful** if and only if it produces at least one of:

- **A mathematical object** (circuit, ring, transform, representation,
  function, or algorithm) that is not equivalent to anything in CLOSED_PATHS,
  with code that runs and a results file.
- **A composition** of ≥ 2 existing edges into a single new statement
  or object, where the composition itself is the novel content.
- **A Lean 4 proof** of an existing edge or novel result, with the proof
  type-checking under Lean.
- **A formal-statement-grade unification** across multiple `novel/`
  results, written as `novel/<name>.md` with theorem/lemma structure.
- **A negative-shape edge** (E7.x family) that closes an entire technique
  family at the structural level (not a single approach).

A session that produces only "DUPLICATE-PLUS" closures of fresh-perspective
brainstorms — without any of the above — is a **failed session**, even if
every CLAUDE.md rule was followed. File the closures and note in the
synthesis: "session produced no novel artifact; primary failure mode was
duplicate ideation."

### Construction Discipline (was the S60 expansion, now the default)

Every construction MUST:

1. **Produce code that runs.** A Python simulator on small inputs counts.
   A Lean file that type-checks counts. A 200-line prose essay does NOT.
2. **Live under a topical experiments/ subdirectory** with `<name>.py` +
   `<name>_results.md`. For new mathematical objects use
   `experiments/constructions/`. For Lean proofs use
   `experiments/formalisations/`.
3. **Carry a "what would falsify this" statement** in the results file.
   If your construction is a non-falsifiable claim, it is not yet a
   construction — sharpen it.
4. **Cite the edges it composes / uses / contradicts**, by ID. A
   composition of E1.5 + E2.1 + E5.8 should say so explicitly.

### Honest Failure Reporting

The single biggest mistake in the project's history has been agents
inflating measurement results into "novel" findings. The policy is now
strict:

- **Confirming a known wall via a new measurement is a measurement, not
  novelty.** Adding a 36th pseudorandomness measure that lands at the
  noise floor is a CLOSED_PATHS entry, not a `novel/` entry.
- **Reproducing an existing edge with a new method is verification, not
  novelty.** Re-deriving E2.1 in a different basis goes in the same edge,
  not a new one.
- **An object you built that turned out to equal an existing one is a
  duplicate.** File honestly as such.

The bar for `novel/` placement: *"a published paper-grade number theorist
or complexity theorist could not, after one careful read of the prior
literature and CLOSED_PATHS, produce this."*

### When you discover a new edge

A new EDGES.md entry requires:

- A verifiable mathematical fact (theorem, identity, exact formula, or
  empirically-saturated bound at meaningful scale).
- An EVS rating (H / M / L / shape).
- A "why this is an edge" justification — what does it constrain or
  enable?
- Cross-references to adjacent edges by ID.

If your finding doesn't meet this bar, it goes in CLOSED_PATHS.md as a
refining row, not in EDGES.md.

### File Placement (STRICT)

- **Experiments** → `experiments/<topic>/` with `<name>.py` +
  `<name>_results.md`. ALWAYS.
- **Constructions** → `experiments/constructions/<descriptive_name>/`
  with `<name>.py` + `<name>_results.md` + optional `definition.md`.
- **Lean formalisations** → `experiments/formalisations/<edge_id>/`
  with `<name>.lean` + `<name>_notes.md`.
- **Novel results** → `novel/<descriptive_name>.md`. Bar: original to
  this project, paper-grade.
- **Session syntheses** → `archive/sessions/sessionNN_<topic>.md`.
- **Proven barriers** → `proven/`.
- **Ephemeral brainstorms** → `archive/ephemeral/`.
- **Literature** → `literature/`.
- **Working algorithms** → `algorithms/` with benchmarks.
- **NEVER** create `<name>_v2.py` / `_quick.py` / `_small.py` variants.
  One script per experiment, parameterised by CLI args.
- **NEVER** put session-by-session details in CLAUDE.md. They go in
  `status/SESSION_INSIGHTS.md`.

### Status File Hygiene (STRICT)

- When an experiment closes a path: add it to CLOSED_PATHS.md in the
  SAME session, with mode (C/E/I) and edge IDs cited.
- When an open problem is resolved: update OPEN_PROBLEMS.md in the SAME
  session.
- When a path / arc / open question is closed: cross-reference it in
  EDGES.md and / or RESEARCH_AGENDA.md.
- Update `status/SESSION_INSIGHTS.md` once per session.

### Cleanup before halting

```
find experiments/ -name "*.py" | while read f; do
  r="${f%.py}_results.md"; [ ! -f "$r" ] && echo "MISSING: $r"
done
```

Run this. Zero tolerance for missing results files. Then:

- Delete any `__pycache__` directories you created.
- Verify every closure has a CLOSED_PATHS row with edge IDs.
- Verify no "pending" labels remain in ephemeral docs for completed work.
- Update `RESEARCH_AGENDA.md` if you advanced or closed an arc.
- Create `archive/sessions/sessionNN_<topic>.md` for this session.

### Session-end self-evaluation (NEW — read before you halt)

Answer these in your session synthesis:

1. **What did I produce that was not in the project before this session?**
   (Object, composition, proof, edge, synthesis, arc-step.)
2. **What edges did my work compose or cite?**
3. **If my session produced only duplicate closures, why?** (Was the
   challenge ill-posed? Was the input frame stale? Was there a tool /
   piece of context I lacked?)
4. **What is the next-action for the next agent?** (One sentence,
   actionable, written into NOVELTY_CHALLENGES.md or RESEARCH_AGENDA.md.)

If question 1's answer is "nothing," do not write a triumphant synthesis.
Write a 5-line note saying so and what blocked you. Honest failure
reporting beats inflated success any day at this stage.

### Rules unchanged from earlier policy

- DO NOT modify `run.sh`. Update `FOCUS_QUEUE.md` only to mark tasks
  COMPLETED or to add new ones.
- If you find the breakthrough, respond with exactly: **I FOUND IT!!!**
- Use sub-agents to save context window when the task is decomposable.
- Context management: when context is filling up, write the in-progress
  state to `RESEARCH_AGENDA.md` (or `TODO.md` for short tasks) so the
  next session can continue, then halt.
- DO NOT add session-by-session details to this file. Put them in
  `status/SESSION_INSIGHTS.md`.
