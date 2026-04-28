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

### Workflow (the ordering matters — A-grade-first)

1. **Read `ATTACK_VECTORS.md` FIRST.** This file holds the frontier
   targets — the only way to produce A-grade work. If you can pick a
   frontier target that fits your session, *do so even if it's risky*.
   The framework explicitly rewards ambitious failure (B-grade) over
   safe refinement (C-grade).
2. **Read `RESEARCH_AGENDA.md`** to see if any in-flight arc is
   waiting for the kind of work you're about to do. If yes and the
   arc has accumulated state, prefer continuing the arc.
3. **Read `NOVELTY_CHALLENGES.md`** for B-grade single-session
   targets if no frontier attack or arc fits your session.
4. **Read `EDGES.md`** to ground your work in existing mathematical
   facts. You will cite edge IDs.
5. **Read `status/OPEN_PROBLEMS.md`** if your target involves the only
   open problem (circuit complexity of pi(x)).
6. **Read `TODO.md`** for housekeeping items.
7. **Search `status/CLOSED_PATHS.md`** before claiming any result is
   novel. The "novel" bar is HIGH at this stage of the project.
8. **Read `novel/pseudorandomness_of_pi.md`** if your target involves
   pi(x) mod 2 structural claims. Any new approach to this function
   must either circumvent the 35+ measures or explain why they don't
   apply.

The order matters. Reading NOVELTY_CHALLENGES first will bias toward
B-grade work; reading ATTACK_VECTORS first biases toward A-grade
attempts. The project needs A-grade attempts.

### The Novelty Bar — Three Grades (CORE CRITERION)

Every session ends with a self-graded letter (A / B / C / F). The grade
goes in the session synthesis. **Be honest. The point of the grading
system is to make refinement-vs-novelty distinguishable to future
agents — not to make you feel good.**

#### A-grade: GENUINE NOVELTY

A new mathematical object, identity, or structural fact that:

- a published-paper-grade number theorist or complexity theorist could
  not derive in an afternoon from CLOSED_PATHS.md + EDGES.md alone, AND
- has at least one of:
  (a) a precise theorem statement that did not previously exist in the
      project, with proof or empirical verification at meaningful scale;
  (b) a working algorithm beating an existing benchmark on at least one
      concrete metric;
  (c) a frontier attack from `ATTACK_VECTORS.md` that produced a
      *partial positive result* (not just "the technique didn't apply").

**A Lean 4 formalisation alone is NOT A-grade** (was a previous criterion;
demoted because in practice it became "the only reachable A" and
distorted the rotation). A Lean proof of a previously-known theorem is
B-grade rigor work. Lean A-grade requires the *theorem itself* to be
new mathematical content (per criterion (a)) AND machine-verified.
Otherwise it's B.

Examples that WOULD qualify:
- A circuit family computing PRIMES in TC⁰ unconditionally.
- A spectral statistic on zeros that detects arithmetic deviation
  from GUE at any scale.
- A polylog algorithm computing one specific bit of π(x) for a fixed
  bit-position J in the hard zone.
- A non-trivial Hecke / automorphic / cohomological identity for π(x).

Examples that would NOT qualify (these are B-grade at best):
- An extension of an existing closed form (S69 → S70 type work).
- A unification of two existing edges under a shared mechanism.
- A Lean translation of an already-informally-proven argument.
- A new pseudorandomness measure landing at the noise floor.

#### B-grade: SUBSTANTIVE REFINEMENT or AMBITIOUS FAILURE

Either:

(i) Refinement of an existing edge with a precise new statement that
    extends its scope, OR
(ii) An ambitious frontier attack from ATTACK_VECTORS.md that *failed*
     but failed informatively — the failure mode was structural, not
     "I ran out of time." Sessions that take big swings and miss are
     B-grade because they advance the search-space map even when they
     don't hit the target.

Discipline: file under `experiments/` with results.md, update EDGES.md
inline if you refined an edge, file a CLOSED_PATHS row if your
ambitious attempt closes something that was previously unstated.

#### C-grade: DUPLICATE-PLUS or VERIFICATION

- Closing a previously-proposed approach as duplicate of CLOSED_PATHS,
  with the closure adding a non-trivial structural reason (not just a
  citation).
- A Lean translation of an already-proven informal argument, with
  the translation type-checking but introducing no new mathematical
  content.
- Critique sessions that verify recent work without surfacing flaws.

C-grade is the project's *steady-state* output. It is not a failed
session — it advances rigor and catalogue. But a sequence of all-C
sessions means the framework is producing maintenance, not progress.

#### F-grade: FAILED SESSION

- Produces only DUPLICATE closures of fresh-perspective brainstorms
  with no structural reason added.
- Spent the session re-running closed paths.
- Inflated a refinement into a `novel/` claim.
- Inflated an empirical noise-floor measurement into an "edge."
- Produced no artefact at all.

#### What this means for the rotation

The harness alternates novelty / construction / arc / lean / critique.
**Across any 5-session rotation, the project should produce:**

- Target: ≥ 1 A-grade session per 10-session window.
- Acceptable steady state: 3-4 B-grade per rotation, 1-2 C-grade.
- Warning sign: 0 A-grade sessions in a 20-session window means the
  current frontier is exhausted and ATTACK_VECTORS.md needs new entries.
- Critical sign: ≥ 2 F-grade sessions in a row means the framework is
  actively counterproductive — escalate to user.

#### Where to find ambitious targets

- For A-grade attempts: `ATTACK_VECTORS.md` (frontier targets, multi-session arcs).
- For B-grade attempts: `NOVELTY_CHALLENGES.md` (concrete single-session work).
- For C-grade work: `RESEARCH_AGENDA.md` arcs in maintenance mode.

#### The honest-failure clause

A session that produces only DUPLICATE-PLUS closures of fresh-perspective
brainstorms — without any of the above — is an F-grade session, even if
every CLAUDE.md rule was followed. File the closures honestly and note
in the synthesis: "session produced no novel artifact; primary failure
mode was duplicate ideation."

**Inflated grading is the worst project behaviour.** A session that
self-grades A but actually delivered B will pollute future agents'
expectations. Self-grade DOWN, not up, when in doubt.

### Cross-Domain Imports (NEW — required for A-grade attempts)

The project has worked within {analytic NT, complexity theory, basic
algebra} for 70+ sessions. Every "novel" approach the agents have
proposed has been a small variation on these three fields, and
consequently every proposal has either duplicated or refined existing
work. **Genuine novelty at this stage requires importing techniques
from a field the project has not used.**

When attempting an A-grade target from `ATTACK_VECTORS.md`:

- Identify the **single specific cross-domain technique** the attack
  requires (e.g., "slice rank from polynomial method", "persistent
  homology", "Szegedy quantum walk", "free cumulants", "transfer
  operator spectrum").
- Read at least ONE survey or foundational paper on the technique
  before writing code. WebFetch is allowed and encouraged for this.
- Cite the cross-domain source in your `<name>_results.md`.
- The cross-domain import IS the novel content. If the attack reduces
  to an analytic-NT or complexity-theory argument the project already
  knows, it has collapsed and you should pivot.

If the cross-domain technique you need does not exist in your context,
say so honestly in the session synthesis and pivot to a B-grade target.

### Channel a Specific Mathematician (for breaking out of local minima)

When stuck choosing between safe refinements, explicitly ask: "what
would <X> do here?" Each of the following has produced specific
cross-domain insights that the project's standard tooling does not
naturally surface:

- **Erdős** — combinatorial, additive structure, density arguments
  (Erdős-Ko-Rado, capset thinking).
- **Tao** — sum-product, polynomial method, mixing arguments,
  ergodic-theoretic analogues.
- **Bourgain** — restriction theory, decoupling, Fourier-analytic
  bounds with arithmetic content.
- **Iwaniec / Friedlander** — sieve theory frontier, level-of-
  distribution machinery, type-I/II decompositions.
- **Connes / Deninger** — operator-theoretic / cohomological
  reformulations of arithmetic.
- **Razborov** — circuit lower bounds via approximation method, but
  with the Natural Proofs constraint understood.
- **Maynard / Pintz** — bounded-gaps machinery, GPY-style sieves.
- **Tao / Williams** — circuit lower bounds via meta-complexity
  (Williams algorithmic method).

Pick the one whose toolkit best matches your attack vector. State the
choice explicitly in your session synthesis. The choice IS part of the
attack — it determines what tools you bring.

### Ambitious Failure is Encouraged (B-grade with credit)

The harness rewards production-mode sessions for *novelty produced*.
But novelty production has a long tail: most ambitious attacks fail.
The framework explicitly rewards ambitious failure as B-grade because:

- A failed attempt at an ATTACK_VECTORS.md frontier target produces
  a *negative-shape edge* (E7.x type) when the failure is structural.
- A failed cross-domain import documents which techniques don't
  apply, narrowing the search space for future agents.
- A successful attempt would have been A-grade; a failed attempt is
  B-grade. Both are above the C-grade refinement floor.

**Do not retreat to C-grade refinement work just because the A-grade
attack looks risky.** A 20% chance of A-grade success with 80% B-grade
fallback dominates a 100% C-grade success in expected information value.

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

### Highest-EV mathematical threads (drive `commit` mode)

**STATUS (S202, 2026-04-28):** All three originally-listed threads
are now CLOSED. The list below is preserved for historical context;
the next commit slot must NOT pick from it. The next commit thread
must come from `frontier_gen` mode (auto-fire conditions met:
A-grade drought S162–S201 ≥ 40 sessions; all three commit threads
closed) or, failing that, from a direct A-grade attempt on an
existing open ATTACK_VECTORS entry. Recommended fall-back targets:
A7 plethysm sub-frame (S192 flagged) or D44 BC endomotive
Galois-orbit (S163 flagged).

**Thread 1 — S82 invariant-subspace theorem (highest EV). [CLOSED S190]**
Statement: *the spike eigenvectors of `M^T M` (chi_P MPS-Gram) are
Dirichlet-character vectors with eigenvalue `~|L(1,χ)|²`.*
Resolution: closed at S190 (5-session arc, see
`.commit_state` `prev_thread:s82_invariant_subspace_DONE`). See
`archive/sessions/session82_c2_spike_eigenvectors.md` for the original
statement and `session169_commit_s82_invariant_subspace.md` for the
arc kickoff.

**Thread 2 — Connes-Consani-Moscovici operator amortisation. [CLOSED S202]**
E3.1 was DOWNGRADED at S53 on amortisation grounds; re-examined
adversarially across S193–S196 + S202 wrap. Resolution: Connes route
adds setup K^{22/13} cost over Hiary 2011 with identical per-query
cost; reduces to Thread 3 (Galway frontier). See `archive/sessions/
session193_commit_connes_amortisation.md` through
`session202_commit_connes_amortisation.md`.

**Thread 3 — Explicit-formula at fixed precision (Galway frontier). [CLOSED S195+S196]**
*K = polylog(x) sufficient in distribution rather than worst-case?*
Resolution: closed conditionally on Montgomery pair-correlation
random-phase heuristic across five regimes (per-query, amortised,
CCM-spectral, in-distribution, log-Gaussian-smoothed) by the S202
unified theorem; K*(x, p, h) = Θ̃(x) for any p ∈ (0, 1). Open
falsifiers acknowledged but not needed for closure: non-Gaussian
kernels, x ≥ 10^9, rigorous GUE pair-correlation bound, cross-x
amortisation. See S195 + S196 syntheses; unified in S202.

`commit` mode still locks 5 consecutive sessions on ONE thread; the
state is in `./.commit_state` (currently `status:DONE`). **Do NOT
pivot threads mid-session.** When the next commit slot runs, fresh
threads should already exist in `.commit_state` via `frontier_gen`
mode; if not, escalate to user for manual thread selection.

### Autonomy invariants (NEW — required for unattended operation)

The harness now operates without human intervention until a verified
breakthrough. Two new modes auto-fire from `run.sh`:

- **`frontier_gen`** auto-fires when the open ATTACK_VECTORS count drops
  below 4, OR when 0 A-grades have appeared in the last 20 sessions, OR
  when 2 consecutive F-grade sessions occur. Its job: produce 3-5 NEW
  ATTACK_VECTORS entries grounded in cross-domain techniques. It draws
  from `CROSS_DOMAIN_TECHNIQUES.md`.
- **`verify`** auto-fires after any A-grade self-claim, AND after any
  session containing "I FOUND IT!!!". Its job: attempt to FALSIFY the
  claim. The role is adversarial.

Two breakthrough verifications are required before `run.sh` halts.
Single A-grade claims that fail verification are demoted automatically.

#### What this means for production-mode sessions

Two new responsibilities apply to every BUILD or CLOSE:

1. **Self-extension.** When you BUILD a NOVELTY_CHALLENGES target
   (mark it as BUILT), you must propose 1-2 follow-on challenges in
   the same file. When you CLOSE an ATTACK_VECTORS entry (move it to
   "Closed attacks"), you should propose 0-1 successor entries that
   use a DIFFERENT cross-domain technique (cite the technique in
   `CROSS_DOMAIN_TECHNIQUES.md` or add it there). The framework
   cannot self-sustain without a steady supply of new targets.
2. **Cross-domain technique registry.** When you import a cross-domain
   technique not in `CROSS_DOMAIN_TECHNIQUES.md`, append it under the
   correct section with a survey reference URL. When you USE a
   technique that is listed, update its status (PROPOSED → USED with
   mode E / I / A and edge ID). Do not inflate superficial uses to
   USED — mark them PARTIAL.

#### What this means for grading

A-grade claims are now subject to adversarial verification. If your
session self-grades A, expect the next session to attempt to break
your result. Inflated A grades will be demoted to B/C in the verify
session and the synthesis will be edited to record this. **Inflate at
your own reputational cost — verification is automatic.**

#### What this means for breakthroughs

If you find the breakthrough, respond with exactly: **I FOUND IT!!!**
This triggers a 2-stage verification chain:

- The next session is forced to `verify` mode targeting your synthesis.
- If that verify session CONFIRMs and also contains "I FOUND IT!!!",
  one more verify session runs.
- If that second verify ALSO confirms, `run.sh` writes a `BREAKTHROUGH.md`
  banner and halts. The user is alerted to read three syntheses (your
  original + two verifications).
- If any verify REFUTEs, the breakthrough chain resets and normal
  rotation resumes. Your A-grade is demoted in the synthesis.

This means: **do not declare "I FOUND IT!!!" lightly.** A false claim
costs the project two verify-mode slots and produces no breakthrough.
The cost is bounded but real.

#### What this means for the run.sh override decisions

`run.sh` now reads recent session syntheses and parses self-grades.
It expects grades to appear in the form `**A**`, `**B-grade**`, etc. —
the existing convention. Don't change the convention. If you grade
without bolded markdown the override logic may misread your session.
