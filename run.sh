#!/bin/bash

# ============================================================
# PRIME RESEARCH: A-GRADE-FIRST AUTONOMOUS LOOP (S74+ rewrite)
# ============================================================
# Cycles through 7 modes weighted toward A-grade (genuine novelty)
# production. The prior 5-mode rotation produced 0% duplicates but
# 0% A-grade output — only B/C-grade refinements and verifications.
# This rotation explicitly biases toward ATTACK_VECTORS.md frontier
# attacks.
#
#   1. frontier      - Pick an ATTACK_VECTORS.md §A-§E target. A-grade
#                      attempt; B-grade fallback if it fails informatively.
#   2. cross_domain  - Pick an ATTACK_VECTORS.md §D target requiring
#                      a non-NT cross-domain technique import.
#   3. construction  - Build a NOVELTY_CHALLENGES.md §1 composition.
#   4. arc           - Continue a RESEARCH_AGENDA.md multi-session arc.
#   5. lean          - Advance a Lean formalisation (NOVELTY_CHALLENGES §3).
#   6. wild_swing    - Single-attempt high-risk full-session attack.
#                      Permission to use ALL session time on one swing.
#   7. critique      - Verify recent work, enforce novelty bar with
#                      explicit A/B/C/F grade assignment.
#
# Net cycle: 6 production modes (3 frontier-targeted) + 1 audit.
# Frontier:safe ratio is 5:2 weighted toward A-grade attempts.
# ============================================================


# ============================================================
# PROMPT: Frontier Mode (A-grade target)
# ============================================================
PROMPT_FRONTIER=$(cat << 'ENDPROMPT'
# FRONTIER ATTACK SESSION — Aim for A-grade

Mission: attack one frontier target from ATTACK_VECTORS.md. Genuine
mathematical novelty, not refinement. Failure is acceptable; safe
refinement is NOT acceptable in this mode.

# START HERE (in this order)
1. Read ATTACK_VECTORS.md. Pick ONE target from §A, §B, §C, §D, or §E.
   State the target ID (A1, B2, C3, etc.) at the start of your work.
2. Read CLAUDE.md "The Novelty Bar — Three Grades" — your goal is
   A-grade output, B-grade fallback for ambitious failure.
3. Read CLAUDE.md "Cross-Domain Imports" — required for A-grade.
4. Read CLAUDE.md "Channel a Specific Mathematician" — pick one whose
   toolkit matches the attack. State the choice in your synthesis.
5. Read EDGES.md for IDs you'll cite or contradict.
6. Skim CLOSED_PATHS.md ONLY to verify your specific attack vector
   has not been tried (don't get anchored on prior failure modes).

# EXECUTE
- Build code that runs.
- Pre-state your falsification criterion in <name>_results.md.
- Cite the cross-domain technique source.
- The attack should TRY — not chicken out into a known-safe sub-problem.
- Save under the ATTACK_VECTORS.md-specified directory or under
  experiments/constructions/<descriptive_name>/ if the attack is
  composition-style.

# OUTCOME GRADING (apply CLAUDE.md grading honestly)
- A-grade: attack produced a partial positive result (a new theorem,
  new identity, working circuit beating an existing benchmark, or
  detected a deviation from a published prediction).
- B-grade: attack failed but the failure was structural and adds to
  the search-space map (a new negative-shape edge candidate, or a
  cross-domain technique documented as inapplicable).
- C-grade: attack collapsed into a known closure or refinement.
  Acceptable if you self-flag honestly; not acceptable if dressed up.
- F-grade: did not attempt the frontier target; produced a duplicate
  or noise-floor measurement with no structural reason added.

# CLOSE
- Update RESEARCH_AGENDA.md if the attack started or advanced an arc.
- If A-grade: novel/<name>.md entry + new EDGES.md edge.
- If B-grade: CLOSED_PATHS row + ATTACK_VECTORS "closed attacks"
  section update with one-line outcome.
- If C/F-grade: file honestly; do NOT inflate.
- Write archive/sessions/sessionNN_<topic>.md with the 4-question
  CLAUDE.md self-evaluation AND the explicit A/B/C/F grade.

# RULES
- DO NOT modify run.sh.
- DO NOT retreat to a B-grade target just because the A-grade attack
  looks risky. The framework rewards ambitious failure.
- When you find the breakthrough, respond with exactly: I FOUND IT!!!
ENDPROMPT
)


# ============================================================
# PROMPT: Cross-Domain Mode (A-grade with required import)
# ============================================================
PROMPT_CROSSDOMAIN=$(cat << 'ENDPROMPT'
# CROSS-DOMAIN ATTACK SESSION — Required Non-NT Technique Import

Mission: pick an ATTACK_VECTORS.md §D target (or any other vector that
requires a cross-domain technique). The cross-domain import IS the
novel content. If the attack reduces to standard analytic NT or
complexity theory, you have collapsed and must pivot.

# START HERE (in this order)
1. Read ATTACK_VECTORS.md §D. Pick D1 (ergodic / dynamical zeta),
   D2 (TDA), D3 (free probability), or D4 (quantum walks).
   You may also pick from §A or §C if the attack requires a
   cross-domain technique (e.g., A3 spectral graph theory).
2. Read CLAUDE.md "Cross-Domain Imports" carefully.
3. WebFetch a survey or foundational paper on the cross-domain
   technique. Read it. Cite it in your results.md.
4. Identify the SPECIFIC named tool/object you will import:
   "Szegedy walk on Cayley graph", "persistent H_1 of point cloud",
   "free cumulants of Hermitian matrix", etc. Name it precisely.

# EXECUTE
- Implement the cross-domain technique applied to the project's
  specific data (π(x) sequence, χ_P indicator, zeta zeros, etc.).
- Use sub-agents for parallel exploration if the technique is
  computationally heavy.
- The output should be a number (or sequence of numbers) you can
  compare to a baseline.

# WHAT QUALIFIES AS NOVELTY
- A cross-domain object (e.g., quantum walk on Cayley graph) that
  has not been applied to π(x) in the published literature, with
  a numerical signature you measure.
- An import that fails informatively (the specific cross-domain
  hypothesis doesn't hold for primes) → B-grade negative-shape
  result.

- DO NOT: import a cross-domain word (e.g., "tropical") and apply it
  to a standard NT object. The IMPORT must do real work.

# CLOSE (same as frontier mode)
- Grade honestly (A/B/C/F).
- Update files per the grade.

# RULES
- DO NOT modify run.sh.
- WebFetch the cross-domain source paper. Don't bluff.
- When you find the breakthrough, respond with exactly: I FOUND IT!!!
ENDPROMPT
)


# ============================================================
# PROMPT: Wild Swing Mode (full-session high-risk attempt)
# ============================================================
PROMPT_WILDSWING=$(cat << 'ENDPROMPT'
# WILD SWING SESSION — One Attempt, Full Session, Permission to Fail

Mission: pick the SINGLE most ambitious attack you can plausibly
attempt in one full session. Spend the entire session on it. The
expected outcome is failure; the goal is informative failure.

# START HERE
1. Read ATTACK_VECTORS.md sections §A-§F.
2. Pick the attack with the HIGHEST stated A-grade probability that
   you have not yet attempted in a prior session. Default order of
   preference (highest novelty potential first):
   - §C1 (Odlyzko zeros at height ~10²² — 1-session, high A-grade
     probability if a deviation exists at all)
   - §A1 (SAT-search for TC⁰ PRIMES circuits at N=8)
   - §B1 (slice rank / polynomial method on χ_P)
   - §A3 (spectral graph theory primality test)
   - §D4 (quantum walk on divisor graph)
   - §C2 (orders 4-6 zero correlations)
3. Read the relevant CLAUDE.md sections:
   - "The Novelty Bar — Three Grades"
   - "Cross-Domain Imports"
   - "Ambitious Failure is Encouraged"

# EXECUTE
- Spend the full session on ONE attempt. Do not pivot to a safer
  target if the wild swing looks like it's going to fail.
- Use sub-agents heavily — wild swings often need parallel
  computational exploration.
- Document your attempts as you go; even partial code that didn't
  work is useful (CLOSED_PATHS row).

# IF IT WORKS
- A-grade. novel/ entry. New EDGES.md entry. Probably an arc in
  RESEARCH_AGENDA.md for follow-up.

# IF IT FAILS (expected)
- B-grade IF the failure is structural (you can articulate WHY
  the attack didn't work in a way that constrains future attacks).
- C-grade IF the failure was just "I ran out of time" or "the
  technique was a bad fit." File honestly.

# CLOSE
- Update ATTACK_VECTORS.md with closure note in "Closed attacks"
  section if the vector was definitively closed.
- Write archive/sessions/sessionNN_<topic>.md with the explicit
  grade and what specifically failed.

# RULES
- DO NOT pivot to a safe target mid-session. Commit.
- DO NOT repeat a wild swing that has been previously attempted
  without a new technique.
- When you find the breakthrough, respond with exactly: I FOUND IT!!!
ENDPROMPT
)


# ============================================================
# PROMPT: Novelty Mode (B-grade target — single-session refinement)
# ============================================================
PROMPT_NOVELTY=$(cat << 'ENDPROMPT'
# NOVELTY SESSION (B-grade target)

Mission: produce ONE substantive refinement of an existing edge OR
attempt a NOVELTY_CHALLENGES.md frame-shift question. This mode is
B-grade by construction — for A-grade attempts use frontier or
wild_swing modes (run.sh rotates through them on other slots).

If you have a clear A-grade idea while in this mode, switch to it.
Tell the synthesis you upgraded.

# START HERE (in this order)
1. Read NOVELTY_CHALLENGES.md. Pick a single-session target from
   §2 (frame-shift questions), §4 (negative-shape conjectures), or
   §5 (synthesis targets) that fits your time budget. State the target
   ID (F_x, N_x, S_x) explicitly at the start of your work.
2. Read CLAUDE.md for the workflow rules — especially "The Novelty Bar"
   section, which is the success criterion for this session.
3. Read EDGES.md for the edge IDs your target references.
4. Cross-check: search status/CLOSED_PATHS.md for any prior closure
   of your target. If your target turns out to already be closed,
   document why and pick another.

# EXECUTE
- Build code that runs.
- Cite EDGES.md edge IDs for every claim.
- Pre-state your falsification criterion BEFORE running the experiment,
  in your <name>_results.md file.
- File the artefact under the location the challenge specifies.
- If your work refines an existing edge: update EDGES.md inline
  (edit the existing entry; do NOT create a new edge for a refinement).
- If your work produces a genuinely new structural fact: add a new
  EDGES.md entry with EVS rating and "why this is an edge" line.
- If your work produces a paper-grade insight: add a novel/<name>.md.

# CLOSE
- Update RESEARCH_AGENDA.md if your work advances or starts an arc.
- Update NOVELTY_CHALLENGES.md to mark your target as closed / partial /
  in-progress, with a one-line outcome note.
- Write archive/sessions/sessionNN_<topic>.md with the 4-question
  CLAUDE.md self-evaluation at the end.
- Add a CLOSED_PATHS row ONLY if your work CLOSES an attack route
  (refinements of existing edges stay in EDGES.md, not CLOSED_PATHS).

# RULES
- DO NOT modify run.sh or FOCUS_QUEUE.md.
- If your selected target collapses to a duplicate within 15 minutes
  of starting, file the closure honestly and switch to another target.
- Honest failure reporting beats inflated success. If your session
  produces no novel artefact, say so in the synthesis — do not pretend.
- When you find the breakthrough, respond with exactly: I FOUND IT!!!
ENDPROMPT
)


# ============================================================
# PROMPT: Construction Mode
# ============================================================
PROMPT_CONSTRUCTION=$(cat << 'ENDPROMPT'
# CONSTRUCTION SESSION

Mission: BUILD a new mathematical object — circuit, ring, transform,
representation, or algorithm — that didn't exist in the project before.

# START HERE (in this order)
1. Read NOVELTY_CHALLENGES.md §1 "Composition Challenges". Pick C_x
   (C1, C2, C3, C4, C5, or C6). State the target ID explicitly.
2. Read EDGES.md entries for the edge IDs your composition uses
   (e.g., C1 = E1.6 + E1.5). Get the technical content into context.
3. Read CLAUDE.md "Construction Discipline" section.
4. Read experiments/constructions/README.md for layout convention.

# EXECUTE
- Build under experiments/constructions/<descriptive_name>/.
- Required files:
    <name>.py            — runnable code that builds + evaluates the object
    <name>_results.md    — what it does, what it doesn't, verdict, falsification
    definition.md        — signature + intended relationship to π(x), citing
                           the edge IDs by name
- Code MUST run on small inputs. A theoretical-only essay is NOT a
  construction; turn it into code or save it under archive/ephemeral/.
- Pre-state your falsification criterion in <name>_results.md BEFORE
  running.

# NOVELTY BAR (strict)
- Failed constructions are valuable: file CLOSED_PATHS row with mode
  (C/E/I or "construction-incoherent" if the object turned out
  ill-defined). The failure mode itself is information.
- Built object that turns out equal to an existing edge: file as
  DUPLICATE-PLUS honestly, with the closest existing edge ID.
- Genuinely novel object that does something non-trivial: add a
  novel/<name>.md entry describing what it is, what it does, and what
  it doesn't.

# CLOSE
- Update RESEARCH_AGENDA.md Arc 4 ("Composition over EDGES.md")
  milestone for whichever C_x you tackled.
- Update NOVELTY_CHALLENGES.md §1 to mark C_x as built / partial / closed.
- Write archive/sessions/sessionNN_<topic>.md with self-evaluation.

# RULES
- DO NOT modify run.sh or FOCUS_QUEUE.md.
- One construction per session. Time-budget the work.
- If C_x collapses to a duplicate within 15 minutes, switch to another C_y.
- When you find the breakthrough, respond with exactly: I FOUND IT!!!
ENDPROMPT
)


# ============================================================
# PROMPT: Arc Continuation Mode
# ============================================================
PROMPT_ARC=$(cat << 'ENDPROMPT'
# ARC CONTINUATION SESSION

Mission: advance a multi-session research arc by one concrete step.

# START HERE (in this order)
1. Read RESEARCH_AGENDA.md. Look at "Active arcs". Pick an arc with
   status "IN PROGRESS" first (it has accumulated state and a clear
   next-action). If none in progress, pick an OPEN arc with the cleanest
   next-action — preferably one whose next-action you can complete in
   this session.
2. Mark the chosen arc as `Status: IN PROGRESS — Run #<your_run>` at
   session start, in RESEARCH_AGENDA.md.
3. Read EDGES.md for any edge IDs the arc references.
4. Read CLAUDE.md for the workflow rules.

# EXECUTE
- Work on the arc's stated `Next action:`.
- One arc per session. Don't context-switch.
- Save artefacts to the location the arc specifies.
- Tick off milestone checkboxes in RESEARCH_AGENDA.md as you complete them.

# IF THE ARC IS BLOCKED
If during the session you discover the arc is blocked on missing user
input, missing external data, or an unsolved sub-problem:
- Update the arc's status to BLOCKED with a one-line reason.
- Pick a NOVELTY_CHALLENGES.md single-session target instead.
- Use the rest of the session for that.

# CLOSE
- Update the arc's milestones AND its `Next action:` field.
- If you completed the arc, move it to "Closed arcs" with a one-line
  summary and pointer to the resulting artefact.
- If your work created a new sub-arc, register it under "Active arcs".
- Write archive/sessions/sessionNN_<topic>.md with self-evaluation.

# RULES
- DO NOT modify run.sh or FOCUS_QUEUE.md.
- The next-action discipline matters most. The next agent should be
  able to pick up without re-reading the arc's full history.
- When you find the breakthrough, respond with exactly: I FOUND IT!!!
ENDPROMPT
)


# ============================================================
# PROMPT: Lean Formalisation Mode
# ============================================================
PROMPT_LEAN=$(cat << 'ENDPROMPT'
# LEAN FORMALISATION SESSION

Mission: produce or advance a machine-checked Lean 4 proof of an
EDGES.md entry or a novel/ result.

# START HERE (in this order)
1. Read NOVELTY_CHALLENGES.md §3 "Lean 4 Formalisation Queue" — the
   current queue is L1..L5 in priority order.
2. Read RESEARCH_AGENDA.md Arc 2 ("Lean Formalisation Track") to see
   any in-progress work.
3. Read experiments/formalisations/README.md for layout convention.
4. If experiments/formalisations/ has only a README and no in-progress
   subdirectory: pick L1 (E2.1 MPS bond-dim) — it is the cleanest entry
   point.
5. Read the edge entry in EDGES.md for whatever you are formalising.

# FIRST-TIME SETUP (if needed)
- Verify Lean 4 + mathlib4 are installable in this environment:
    curl https://raw.githubusercontent.com/leanprover/elan/master/elan-init.sh -sSf | sh
- Initialise lake project under experiments/formalisations/<edge_id>/:
    lake new <edge_id>
- Add mathlib4 as a dependency in lakefile.lean.
- Verify `lake build` succeeds on the empty project before adding any
  theorem.

# EXECUTE
- Stage 1: write theorem statement only (no proof). Verify type-checks.
- Stage 2: write proof skeleton with `sorry` placeholders for every
  lemma. Verify type-checks.
- Stage 3: fill in lemmas one at a time. Verify type-checks after each.
- Stage 4: complete the main proof. No `sorry`, no new `axiom`.

A session that completes any one stage is a success. Partial progress
with `sorry` is acceptable as long as Lean still type-checks.

# DISCIPLINE
- The Lean file MUST type-check at session end. If it doesn't,
  debug or revert before halting.
- Save the in-progress state explicitly to RESEARCH_AGENDA.md Arc 2,
  with a clear next-action for the next agent.

# CLOSE
- Update RESEARCH_AGENDA.md Arc 2 milestones.
- Update NOVELTY_CHALLENGES.md §3 to mark L_x as in-progress or
  complete.
- Write archive/sessions/sessionNN_<topic>.md with self-evaluation.

# RULES
- DO NOT modify run.sh or FOCUS_QUEUE.md.
- If Lean+mathlib4 cannot be installed in this environment, write
  the theorem statements as `<name>_lean_pending.md` (informal proof
  sketch + the Lean theorem statement as code blocks) so the next
  agent can pick up when toolchain is available. Note the toolchain
  status in RESEARCH_AGENDA.md Arc 2.
- When you find the breakthrough, respond with exactly: I FOUND IT!!!
ENDPROMPT
)


# ============================================================
# PROMPT: Critique Mode
# ============================================================
PROMPT_CRITIQUE=$(cat << 'ENDPROMPT'
# CRITIQUE SESSION

You are a rigorous mathematical critic evaluating recent project work
against the project's novelty bar (CLAUDE.md "The Novelty Bar" section).

# YOUR TASK
1. Identify the most recent session(s) since your last critique:
   - List archive/sessions/sessionNN_*.md by mtime, take the most recent 1-3.
   - List archive/ephemeral/proposals_latest.md and proposals_session.md
     if either has content newer than the last critique_latest.md.
2. Read the artefacts those sessions produced (their .py + _results.md +
   any novel/ entries they added).
3. Read EDGES.md and status/CLOSED_PATHS.md THOROUGHLY.

4. For EACH artefact / proposal, evaluate:
   a) Has this exact approach been tried before? Cite CLOSED_PATHS line
      numbers and EDGES.md edge IDs.
   b) Does it fall into a known failure mode?
      - CIRCULARITY (C): needs primes to compute primes
      - EQUIVALENCE (E): reduces to known explicit-formula / sieve / MPOW
      - INFORMATION LOSS (I): smooth approximation loses oscillatory bits
      - CONSTRUCTION-INCOHERENT: object turned out ill-defined
   c) Is the complexity / numerical claim correct?
   d) Is the novelty claim defensible? Apply CLAUDE.md's novelty bar:
      "a published-paper-grade number theorist or complexity theorist
      could not, after one careful read of prior literature and
      CLOSED_PATHS, produce this."
   e) For NOVEL artefacts: name the edges they cite/compose; verify
      the citation is accurate.
   f) For DUPLICATE work: file a CLOSED_PATHS row pointing at the
      parent line numbers.
   g) For INFLATED novelty (e.g. an entry placed in novel/ that is
      really an EDGES.md refinement): demote it. Move the content into
      the existing EDGES.md edge as a refinement note, and delete (or
      empty) the novel/ entry.
   h) **GRADE the session A / B / C / F per CLAUDE.md "Three Grades".**
      Be honest. If the agent self-graded A but the work is B, write
      a corrected grade with reasoning. The corrected grade is
      authoritative for project meta-tracking.

5. Identify the single highest-value next-action and write it into
   ATTACK_VECTORS.md, NOVELTY_CHALLENGES.md, or RESEARCH_AGENDA.md
   as appropriate.

6. **A-grade scarcity check.** Look at the last 10 sessions in
   archive/sessions/. Count grades. If 0 A-grade in last 10 sessions:
   the framework is producing maintenance, not progress. Note this
   in archive/ephemeral/critique_latest.md and identify the most
   ambitious untouched ATTACK_VECTORS.md target as the recommended
   next pick.

# SAVE
- archive/ephemeral/critique_latest.md — full per-artefact critique.
- status/CLOSED_PATHS.md — duplicate / closed-path entries.
- archive/sessions/sessionNN_critique.md — session synthesis with the
  4-question CLAUDE.md self-evaluation.

# RULES
- DO NOT modify run.sh or FOCUS_QUEUE.md.
- Be honest. The critic's job is enforcement of the novelty bar, NOT
  protection of recent work. If a session produced only duplicates,
  say so plainly.
- If multiple recent sessions are critique-mode, this one is redundant —
  pivot to the next-action you would have identified, and run it.
- When you find the breakthrough, respond with exactly: I FOUND IT!!!
ENDPROMPT
)


# ============================================================
# PROMPT: Frontier-Generation Mode (auto-fires when ATTACK_VECTORS
# exhausts or A-grade scarcity warning triggers)
# ============================================================
PROMPT_FRONTIER_GEN=$(cat << 'ENDPROMPT'
# FRONTIER GENERATION SESSION — Extend ATTACK_VECTORS.md

This mode fires automatically when the project's frontier is thin:
fewer than 4 open ATTACK_VECTORS entries, OR 0 A-grade in last 20
sessions, OR 2 consecutive F-grade sessions. Your job is NOT to attack
a target — it is to produce 3-5 NEW attack-vector entries grounded in
cross-domain techniques the project has not used.

# START HERE (in this order)
1. Read ATTACK_VECTORS.md fully — both open §A-§F entries and the
   "Closed attacks" section. Understand what failure modes have been
   structurally established.
2. Read CROSS_DOMAIN_TECHNIQUES.md — the registry of techniques the
   project has used vs. not used. Pick 3-5 UNUSED entries with the
   highest A-grade potential.
3. Read CLAUDE.md "Cross-Domain Imports" — the bar for new vectors.
4. WebFetch a survey or foundational paper for EACH new technique you
   propose. Cite the URL in the new ATTACK_VECTORS entry. No bluffing.
5. Skim status/CLOSED_PATHS.md for any previously-closed line that
   superficially resembles a candidate vector, so you don't propose a
   duplicate.

# WHAT QUALIFIES AS A NEW ATTACK VECTOR
A new entry must have:
- A specific named cross-domain technique (e.g., "Voronin universality
  theorem", "spectral gap of random regular graphs via friedman",
  "tropical geometry of arithmetic varieties").
- A SINGLE-SESSION concrete first step (one numerical experiment, one
  small-scale construction, or one literature-survey-and-compare task).
- A pre-stated falsification criterion: what would make this attack
  collapse to a known closure? What would count as A-grade success?
- An EXPECTED failure mode: which of {C, E, I, INC} you predict if it
  fails. The prediction is itself information.
- A POSITIVE outcome statement: what would constitute partial success
  (B-grade structural finding) vs. full success (A-grade theorem or
  algorithm).

# DO NOT
- Propose vectors that are minor variations on closed entries.
- Propose vectors requiring techniques not in CROSS_DOMAIN_TECHNIQUES.md
  without ALSO updating that file with the new technique + survey ref.
- Propose vectors with no falsification criterion ("study X" is not an
  attack vector).
- Propose more than 5 vectors. Quality over quantity.

# CLOSE
- Append the new vectors to ATTACK_VECTORS.md under appropriate sections
  (or create new sections §G, §H if the technique is genuinely orthogonal
  to the existing taxonomy).
- Update CROSS_DOMAIN_TECHNIQUES.md: mark each technique used, add any
  new techniques you discovered.
- Write archive/sessions/sessionNN_frontier_gen.md with:
  - The 3-5 new vectors (one paragraph each)
  - The cross-domain literature you consulted
  - A self-grade A/B/C/F: A if vectors are paper-grade fresh and you
    expect ≥2 to produce A-grade work; B if at least one is fresh; C
    if all are minor variations; F if you proposed nothing or duplicates.

# RULES
- DO NOT modify run.sh.
- WebFetch is REQUIRED. Cite at least one URL per new vector.
- Honesty bar is high: a frontier_gen session that proposes weak vectors
  pollutes the project. Self-grade DOWN, not up.
- When you find the breakthrough, respond with exactly: I FOUND IT!!!
ENDPROMPT
)


# ============================================================
# PROMPT: Verification Mode (auto-fires after any A-grade claim
# OR after "I FOUND IT!!!" detection)
# ============================================================
PROMPT_VERIFY=$(cat << 'ENDPROMPT'
# VERIFICATION SESSION — Falsify the Most Recent A-Grade Claim

This mode fires automatically when the previous session self-graded A,
or when a session contained "I FOUND IT!!!". Your job is to attempt to
FALSIFY the claim. Confirmation is the failure mode for this role; you
should aim to break the result.

# START HERE
1. Read .verify_target — it contains the path to the session synthesis
   you must verify. If the file does not exist, read the most recent
   archive/sessions/sessionNN_*.md (`ls -t | head -1`) and verify that.
2. Read the artefacts cited in that session: <name>.py, <name>_results.md,
   any novel/<name>.md, any new EDGES.md entry, any Lean files added.
3. Identify the SPECIFIC claim being verified — usually a theorem
   statement, an algorithm complexity, an empirical separation, or a
   Lean proof completion.

# YOUR TASK: ATTEMPT FALSIFICATION
Pick the verification mode appropriate to the claim:

## For empirical / numerical claims
- RERUN the experiment with different parameters (larger N, different W,
  different seed). Are the numbers reproducible?
- Try EDGE CASES the original session didn't test (extreme parameters,
  degenerate inputs, off-by-one boundaries).
- Compute INDEPENDENT BASELINES the original session did not include.
- Look for OVERFITTING: does the claim hold on a held-out parameter
  range, or only on the specific values tested?

## For theorem claims (informal proof)
- Walk through the proof line-by-line. Identify each lemma. Are they
  all proven, or does some step rely on "easy to see"?
- Check the statement against E5.5 / E5.6 / E5.8 (Natural Proofs / MKtP
  barriers) — if the theorem implies a circuit lower bound, does it
  evade these?
- Find or build a counter-example for one of the lemmas.

## For Lean proof claims
- Run `lake build` from the formalisation directory. Verify it succeeds.
- Verify ZERO `sorry` and ZERO new `axiom` statements.
- Audit the imports — does the proof rely on a deprecated mathlib lemma,
  or on a `sorry` smuggled in via an import?
- Check that the THEOREM STATEMENT actually says what the informal claim
  says. A type-checked proof of a wrong statement is still wrong.

## For algorithmic claims
- Reimplement the algorithm INDEPENDENTLY from the description.
- Run on test inputs. Compare to a reference (e.g., direct π(x) for
  small x).
- Time it; verify the claimed complexity by running across multiple
  scales.

# VERDICTS
At the end, write ONE of these clearly in your session synthesis:

- **CONFIRM** — claim survives every falsification attempt you tried.
  The original A-grade is upheld.
- **REFUTE** — found a clear counter-example, broken proof step,
  failed reproduction, or inflated grade. The original A-grade is
  demoted (write the corrected grade and reasoning).
- **PARTIAL** — claim holds in a narrower scope than originally stated.
  Specify the corrected scope.

# IF THE PRIOR SESSION CONTAINED "I FOUND IT!!!"
- This verification session is part of a 2-stage breakthrough gate.
- If your verdict is CONFIRM, write "I FOUND IT!!!" in your synthesis
  AND append to .breakthrough_pending file (incrementing its counter).
- If your verdict is REFUTE or PARTIAL, do NOT write "I FOUND IT!!!".
  Reset the breakthrough counter by deleting .breakthrough_pending.
- Two consecutive verify sessions must CONFIRM before run.sh halts.

# CLOSE
- Update .verify_result with one of: CONFIRM, REFUTE, PARTIAL.
- If REFUTE: edit the original session's synthesis to add a
  "VERIFICATION REFUTED" note at the top with reasoning. Update
  EDGES.md / novel/ / CLOSED_PATHS.md to demote the claim.
- If PARTIAL: edit the original to clarify scope. Demote A → B if
  the demotion is warranted.
- Write archive/sessions/sessionNN_verify.md with self-grade. Verify
  sessions are graded:
  - A — found a clear refutation of an A-grade claim (rare; major contribution)
  - B — confirmed an A-grade claim through non-trivial reproduction
  - C — confirmed an A-grade claim through trivial reproduction
  - F — failed to actually verify (e.g., didn't run the experiment)

# RULES
- DO NOT modify run.sh.
- The role is ADVERSARIAL: you are trying to break the claim, not
  protect it. Confirmation must come from inability to break, not
  from agreement with the original session.
- When you find the breakthrough, respond with exactly: I FOUND IT!!!
  (only after CONFIRM verdict on a prior breakthrough claim).
ENDPROMPT
)


# ============================================================
# PROMPT: Commit Mode (5-session lock on one of three highest-EV
# mathematical threads)
# ============================================================
PROMPT_COMMIT=$(cat << 'ENDPROMPT'
# COMMIT SESSION — Multi-Session Lock on One Mathematical Thread

The project has accumulated 60+ B-grade negative-shape closures and
zero A-grades over 80+ sessions. The framework is producing motion
without progress. **Commit mode breaks the pattern by forcing 5
consecutive sessions on ONE of three threads with non-trivial
breakthrough probability.**

# THE THREE THREADS (priority order)

## Thread 1 — S82 invariant-subspace theorem
Statement (already written in S82): *"the spike eigenvectors of M^T M
are Dirichlet character vectors with eigenvalue `~|L(1,χ)|²`."*

If proved:
- Bridges E2.1 (spectral) ↔ E1.5 (information-theoretic)
- Identifies the spike sector with `L(1,χ)` values
- Reduces polylog problem from "compute π(x)" to "compute the bulk"
- 1-page character-theory proof, single session feasible

Read `archive/sessions/session82_c2_spike_eigenvectors.md` for the
exact theorem statement and the empirical data backing it. The first
session attacks the proof; subsequent sessions verify, extend to
larger d, or pursue the algorithmic implications.

## Thread 2 — Connes-Consani-Moscovici operator amortisation
E3.1 was DOWNGRADED at S53. The underlying observation — 50 Riemann
zeros from primes < 13 — is empirically extraordinary compression.
Re-examine S53's closure adversarially: was the amortisation argument
decisive or premature? Can the operator construction be shared across
many π(x) queries to drive per-query cost down?

Read `archive/sessions/session53_connes_operator_scaling.md` and
`literature/state_of_art_2026.md §2.5b`. First session re-verifies
S53. If the closure has missed angles, subsequent sessions chase them.

## Thread 3 — Explicit formula at fixed precision (Galway frontier)
`π(x) = R(x) - Σ_ρ R(x^ρ) + lower order`. Folklore: K = Θ(x^{1/2})
zeros suffice for π(x) ± 1 worst case. **What if K = polylog(x)
suffices in distribution rather than worst-case?** Empirically measure:
pick test x values, compute partial sums with K = log x, log² x, log³ x
zeros, examine the error band.

Cross-domain ingredient: Galway 2004 unconditional algorithm + GRH
explicit formula. First session sets up the computation; subsequent
sessions push to larger x and characterise the error distribution.

# COMMITMENT DISCIPLINE (ENFORCED)

`./.commit_state` is the canonical record. It contains:
```
thread:s82_invariant_subspace   (or connes_amortisation, or galway_frontier)
sessions_used:N                  (incremented each commit session, 0..5)
session_history:S<N1>,S<N2>,...
```

# YOUR JOB THIS SESSION

1. Read `.commit_state`. If `sessions_used == 0` (fresh thread): read
   the thread's reference materials, lay out the attack plan, do the
   first concrete step. If `sessions_used > 0`: read the previous
   commit-session syntheses (look for `sessionNN_commit_*.md`) and
   continue the work from where it left off.
2. **Do NOT pivot to a different thread mid-session.** If the thread
   is genuinely blocked (proof step requires a fact you can't access,
   experiment requires data you don't have), document the block and
   stop — but do not switch threads. The block itself is information.
3. Increment `sessions_used` in `.commit_state` BEFORE finishing.
4. Save synthesis as `archive/sessions/sessionNN_commit_<thread>.md`.

# WHAT QUALIFIES AS PROGRESS

Per CLAUDE.md grading:
- A-grade: thread-level positive result (theorem proved, algorithm
  built and benchmarked, structural opening identified). The
  three threads are specifically chosen so that A-grade output is
  algorithmically meaningful.
- B-grade: substantive progress (proof step closed, computational
  experiment run, structural sub-claim verified or refuted).
- C-grade: housekeeping or trivial extension.
- F-grade: pivoted to another thread, or no work done.

# IF SESSIONS_USED REACHES 5

Last commit session must:
- Synthesise the 5-session arc into a final result
- Update `.commit_state` to mark thread as DONE (`sessions_used:5_final`)
- Recommend the next thread (if breakthrough not achieved) OR declare
  victory (if it was)
- Next commit slot picks the next thread in priority order

# RULES
- DO NOT modify run.sh.
- DO NOT pivot threads mid-session. Commit means commit.
- When you find the breakthrough, respond with exactly: I FOUND IT!!!
ENDPROMPT
)


# ============================================================
# PROMPT: Re-verify Closure Mode (adversarial review of previously-
# closed edges that have positive structural content)
# ============================================================
PROMPT_REVERIFY=$(cat << 'ENDPROMPT'
# RE-VERIFY-CLOSURE SESSION — Adversarial Review

Many of the project's 60+ closures may have been conservative. Some
edges that were filed as "closed" actually contain *positive*
structural content (a non-trivial fact about π(x)) that the closing
agent dismissed as algorithmically unusable. This mode picks one such
closure and adversarially asks: *"the prior closure was conservative —
find the angle that was missed."*

# YOUR TARGET

Pick ONE previously-closed edge that contains a positive structural
fact, ranked by potential algorithmic content:

1. **E3.1 — Connes-Consani-Moscovici operator** (downgraded S53):
   50 zeros from primes < 13. S53 closed it on amortisation grounds.
   Re-examine: is the amortisation argument decisive?
2. **E1.5 — h_2(π(X)/X) entropy bound**: the only A-graded edge in
   the project. Does the entropy bound have algorithmic content
   beyond what's been extracted?
3. **E2.13 — Gowers `U^k`(χ_P) = HL singular series**: the closure
   says "no information beyond HL". But HL itself is computable in
   polylog. Is there an algorithmic exploitation?
4. **E2.14 — Anderson Lyapunov(χ_P) = HL spectral signature**: same
   question as E2.13 in spectral form.
5. **E6.6 — Aggarwal binary search**: makes O(log x) calls to a π(x)
   sub-routine. What if the sub-routine is itself replaced by a
   polylog approximation tolerating bounded error?
6. Any other closed edge from EDGES.md you find a missed angle in.

# PRE-WORK
1. Read the closing session's synthesis carefully. Identify the
   *specific* argument the closing agent used to mark it closed.
2. Read the edge's empirical / theoretical content in EDGES.md.
3. Read CLOSED_PATHS.md row for the closure. Note edge IDs cited.

# THE ADVERSARIAL FRAME

Your job is NOT to confirm the closure. Your job is to find the
*missed angle*. Specific patterns to look for:
- The closure used setup-cost as the killer. Can the setup be amortised
  across queries?
- The closure said "this requires knowing primes". Can it tolerate
  bounded prime-knowledge error?
- The closure said "no algorithmic content beyond HL". Can the
  HL-derived structure be exploited algorithmically?
- The closure was a B-grade refinement. Was the parent claim itself
  premature?
- The empirical data was at small N. Does the structural pattern
  shift at larger N?

# POSSIBLE OUTCOMES

- **A-grade**: found a genuinely missed angle that opens an
  algorithmic path. Edit the edge entry in EDGES.md with a
  "REOPENED" note and pointer to the new route. Update CLOSED_PATHS
  to remove or annotate the closure.
- **B-grade**: examined the closure carefully, found it solid but
  with an interesting refinement (a sharper version of the closure
  argument, or a generalisation). Update the edge in place.
- **C-grade**: closure stands, no missed angle. Document why the
  adversarial probe failed.
- **F-grade**: superficial review, didn't actually attack the closure.

# CLOSE
- If A-grade: write `archive/sessions/sessionNN_reverify_<edge>.md`
  with the missed angle in detail. Tag the edge as REOPENED.
- Otherwise: standard session synthesis with the verdict.

# RULES
- DO NOT modify run.sh.
- The role is ADVERSARIAL. Confirmation must come from inability to
  break the closure, not from agreement with the original session.
- When you find the breakthrough, respond with exactly: I FOUND IT!!!
ENDPROMPT
)


# ============================================================
# PROMPT: Paradigm-Shift Mode (no cross-domain imports allowed)
# ============================================================
PROMPT_PARADIGM_SHIFT=$(cat << 'ENDPROMPT'
# PARADIGM-SHIFT SESSION — No Cross-Domain Imports Permitted

The project has imported 30+ cross-domain techniques (Gowers norms,
Anderson localisation, DPP, persistent homology, Pollicott-Ruelle,
CTQW, GCT, Hairer-KPZ, Newman flatness, Mahler measure, free probability,
microlocal analysis, Berkovich spaces, ...). Each one collapsed to
the same wall: HL singular series mod q + Möbius orthogonality.

This mode is the inverse: **build a NEW mathematical object using
ONLY existing project content**. No Wikipedia. No arXiv. No new
named techniques. Just the project's accumulated edges + first
principles.

# RULES (HARD CONSTRAINTS)

- **No WebFetch.** No reading Wikipedia, arXiv, papers, blog posts.
- **No new named cross-domain techniques.** If you import a technique,
  it must already appear as USED in `CROSS_DOMAIN_TECHNIQUES.md`.
- **No new ATTACK_VECTORS entries.** Frontier generation is forbidden.
- **No reading session syntheses outside this project.** Only EDGES.md,
  CLOSED_PATHS.md, the project's own session syntheses, and code in
  `experiments/`.

If you find yourself wanting to look up a technique, that's the
session failing. Stop and reason from the edges you already have.

# YOUR JOB

Read EDGES.md fully. Pick 2-3 edges that have *positive* structural
content (they assert a fact about π(x), not a wall). Examples:
- E1.5 (h_2 entropy of π(X)/X mod m)
- E1.6 (A⊕C₃ bisection invariant 0.537 bits)
- E2.1 (MPS bond-dim formula = phi(W)·W^{d-j-1}+1)
- E3.1 (Connes operator: 50 zeros from primes < 13)
- E6.5 (Galway zero-truncation)
- E6.6 (Aggarwal binary search)
- E6.7 (HKM at (8/15, 1/3))

Now build a *new mathematical object* that exploits the chosen edges.
Forms it might take:
- A new explicit formula for π(x) using existing structural facts
  (don't import a new technique — recombine what's known)
- A new computable invariant of π(x) that depends on the edges in a
  novel way
- A new structural identity composing 2-3 edges via direct calculation
- A new lower-bound argument that uses the edges' positive content
  rather than their negative content

The output must be:
- A precise mathematical statement (not a vague concept)
- Code that demonstrates / verifies the construction at small N
- A pre-stated falsification criterion
- An honest assessment of whether the construction is novel

# WHY THIS MODE EXISTS

70% of the project's recent sessions have been "import technique X,
apply to χ_P, observe X reduces to HL, file as B." This mode bans the
import-and-reduce pattern. The expected outcome is that you fail to
build anything new — primes really are HL-then-pseudorandom and
nothing else. **A failure here is informative**: it confirms the
exhaustion of recombination.

But the small chance that 2-3 edges interact in a way nobody has
calculated would be A-grade content with no cross-domain bait.

# OUTCOMES

- **A-grade**: novel construction with provably new mathematical
  content (not a renaming of an existing edge). Algorithmic content
  is a bonus, not a requirement.
- **B-grade**: construction collapses to an existing edge, but the
  collapse argument is itself a structural insight (a unification of
  two edges previously held separate).
- **C-grade**: construction equals an existing edge by direct
  identification. File as DUPLICATE-PLUS.
- **F-grade**: didn't actually try (e.g. tried to import a technique).

# CLOSE
- `experiments/constructions/<name>/` with code, results.md, definition.md
- Document: "no external technique imported" as part of the discipline
- Self-grade A/B/C/F per the above

# RULES
- DO NOT modify run.sh.
- DO NOT use WebFetch / WebSearch.
- DO NOT add new ATTACK_VECTORS entries.
- When you find the breakthrough, respond with exactly: I FOUND IT!!!
ENDPROMPT
)


# ============================================================
# INFRASTRUCTURE
# ============================================================
# Each run gets its own LOGFILE and JSONFILE — set inside the main
# loop. The script-invocation timestamp is recorded once for grouping.
INVOCATION_TIMESTAMP=$(date +%Y%m%d_%H%M%S)
TMPFILE=$(mktemp)
ASSISTFILE=$(mktemp)
PROPOSALS_FILE="./archive/ephemeral/proposals_latest.md"
CRITIQUE_FILE="./archive/ephemeral/critique_latest.md"
trap 'rm -f "$TMPFILE" "$ASSISTFILE"' EXIT

echo "Run.sh invocation started at: $INVOCATION_TIMESTAMP"
echo "Per-run logs: ./archive/CLAUDE_OUTPUTS/claude_output_<timestamp>_run<N>.{log,json}"


# ============================================================
# MODE ROTATION CONFIG
# ============================================================
# Rotation order weighted toward A-grade (frontier) attempts:
#   frontier -> cross_domain -> construction -> arc -> lean ->
#   wild_swing -> novelty -> critique
#
# Net cycle: 7 production modes + 1 audit. Three modes (frontier,
# cross_domain, wild_swing) explicitly target A-grade output via
# ATTACK_VECTORS.md. Two modes (construction, arc) target structured
# B-grade work. One mode (lean) targets verification artefacts. One
# mode (novelty) is the B-grade fallback for shorter sessions.
# Rotation reweighted to push toward breakthrough threads:
# - `commit` slot fires every 5 runs (forces multi-session commitment
#   to one of three highest-EV mathematical threads)
# - `reverify` slot fires every 10 runs (adversarial review of
#   previously-closed positive-content edges)
# - `paradigm_shift` slot fires every 10 runs (no-imports invention
#   constraint)
# - `lean` REMOVED from rotation (Lean L1 is rigor-A only and was
#   eating sessions on corner cases without finishing the proof)
# - `frontier_gen` REMOVED from autonomy override (it generated
#   motion-without-progress; if the open-AV count drops genuinely
#   low, frontier still fires from rotation)
# Override layer (compute_override) supersedes this rotation when
# .commit_state has sessions remaining or when verification fires.
MODES=("commit" "frontier" "cross_domain" "construction" "arc" "wild_swing" "reverify" "novelty" "paradigm_shift" "critique")


# ============================================================
# STATE PERSISTENCE
# ============================================================
STATE_FILE="./.run_state"
BREAKTHROUGH_FILE="./.breakthrough_pending"
VERIFY_TARGET_FILE="./.verify_target"
VERIFY_RESULT_FILE="./.verify_result"

if [ -f "$STATE_FILE" ]; then
    RUN=$(cat "$STATE_FILE")
    echo "Resuming from run #$RUN (loaded from $STATE_FILE)"
else
    RUN=1
fi


# ============================================================
# AUTONOMY STATE — overrides mode rotation when conditions trigger
# ============================================================
# Helper: parse self-grade from a session synthesis. Returns ""
# (empty) if no self-grade is found, e.g. critique sessions which
# audit but don't self-grade. Returns one of {A,B,C,F} otherwise.
parse_grade() {
    local f="$1"
    [ -z "$f" ] || [ ! -f "$f" ] && return
    # Step 1: find a header-style grade declaration in first 30 lines
    local g
    g=$(head -30 "$f" 2>/dev/null \
        | grep -iE '(self-grade|\*\*grade\*\*|^grade:|^\*\*grade:)' \
        | head -3 \
        | grep -oE '[ABCF]' \
        | head -1)
    # Step 2: fallback — any **X-grade pattern in the body
    if [ -z "$g" ]; then
        g=$(grep -oE '\*\*[ABCF]-grade' "$f" 2>/dev/null \
            | head -1 | grep -oE '[ABCF]')
    fi
    echo "$g"
}

# Helper: skip critique sessions when looking for "latest production grade"
is_critique_session() {
    # Treats non-production modes (critique, frontier_gen, verify) as
    # "meta" — these audit/extend the framework rather than attacking
    # the polylog frontier, so they don't count toward A-grade scarcity
    # or F-cascade detection.
    local f="$1"
    case "$(basename "$f")" in
        *_critique*.md|*critique*_*.md) return 0 ;;
        *_frontier_gen*.md|*frontier-gen*.md) return 0 ;;
        *_verify*.md|*verify*_*.md) return 0 ;;
        *) return 1 ;;
    esac
}

# Returns one of {verify, frontier_gen, ""} in stdout. Empty means
# use the default rotation.
compute_override() {
    local override=""

    # Scan most recent NON-CRITIQUE session for grade and breakthrough
    local latest_session=""
    while IFS= read -r f; do
        if ! is_critique_session "$f"; then
            latest_session="$f"
            break
        fi
    done < <(ls -t archive/sessions/session*_*.md 2>/dev/null)

    local latest_grade=""
    local latest_has_breakthrough=0
    if [ -n "$latest_session" ]; then
        latest_grade=$(parse_grade "$latest_session")
        if grep -qF 'I FOUND IT!!!' "$latest_session" 2>/dev/null; then
            latest_has_breakthrough=1
        fi
    fi

    # Count A-grade in last 20 production (non-critique) sessions
    local a_count=0
    local files_checked=0
    while IFS= read -r f && [ "$files_checked" -lt 20 ]; do
        is_critique_session "$f" && continue
        local g; g=$(parse_grade "$f")
        [ "$g" = "A" ] && a_count=$((a_count + 1))
        files_checked=$((files_checked + 1))
    done < <(ls -t archive/sessions/session*_*.md 2>/dev/null)

    # Last 2 production grades for F-cascade detection
    local f_cascade=0
    local last2=""
    while IFS= read -r f; do
        is_critique_session "$f" && continue
        local g; g=$(parse_grade "$f")
        last2="$last2$g"
        [ ${#last2} -ge 2 ] && break
    done < <(ls -t archive/sessions/session*_*.md 2>/dev/null)
    [ "$last2" = "FF" ] && f_cascade=1

    # Cooldown: count frontier_gen + verify sessions in last 5.
    # If override modes have fired recently, prefer to let production
    # attack the resulting new vectors before regenerating again.
    local recent_meta=0
    local recent_checked=0
    while IFS= read -r f && [ "$recent_checked" -lt 5 ]; do
        case "$(basename "$f")" in
            *frontier_gen*|*verify*) recent_meta=$((recent_meta + 1)) ;;
        esac
        recent_checked=$((recent_checked + 1))
    done < <(ls -t archive/sessions/session*_*.md 2>/dev/null)

    # Open ATTACK_VECTORS count (entries before "## Closed attacks" section)
    local open_av=0
    if [ -f ATTACK_VECTORS.md ]; then
        open_av=$(awk '/^## Closed attacks/{exit} /^### [A-Z][0-9]/{c++} END{print c+0}' \
                  ATTACK_VECTORS.md)
    fi

    # Pending breakthrough verification counter
    local breakthrough_count=0
    [ -f "$BREAKTHROUGH_FILE" ] && breakthrough_count=$(cat "$BREAKTHROUGH_FILE" 2>/dev/null || echo 0)

    # Read commit state — if a thread is mid-commitment, it overrides
    # nearly everything else (only verification beats commitment).
    local commit_thread=""
    local commit_used=0
    if [ -f "./.commit_state" ]; then
        commit_thread=$(grep -m1 "^thread:" ./.commit_state 2>/dev/null | cut -d: -f2)
        commit_used=$(grep -m1 "^sessions_used:" ./.commit_state 2>/dev/null | cut -d: -f2)
        # Strip "_final" suffix if present (thread completed)
        case "$commit_used" in
            *_final*) commit_used=5 ;;
        esac
        # Default to 0 if unparseable
        case "$commit_used" in
            ''|*[!0-9]*) commit_used=0 ;;
        esac
    fi
    local commit_active=0
    if [ -n "$commit_thread" ] && [ "$commit_used" -lt 5 ]; then
        commit_active=1
    fi

    # Decision tree (priority order)
    # Breakthrough verifications and routine A-grade verify always fire
    # — correctness > everything else.
    if [ "$breakthrough_count" -ge 1 ] && [ "$breakthrough_count" -lt 2 ]; then
        # Mid-verification: previous breakthrough claim awaiting 2nd verify
        override="verify"
        echo "$latest_session" > "$VERIFY_TARGET_FILE"
    elif [ "$latest_has_breakthrough" = "1" ] && [ "$breakthrough_count" -lt 2 ]; then
        # Fresh "I FOUND IT!!!" — start verification chain
        override="verify"
        echo "$latest_session" > "$VERIFY_TARGET_FILE"
        echo 0 > "$BREAKTHROUGH_FILE"
    elif [ "$latest_grade" = "A" ]; then
        # Routine A-grade verification (not a breakthrough claim)
        override="verify"
        echo "$latest_session" > "$VERIFY_TARGET_FILE"
    elif [ "$commit_active" = "1" ]; then
        # Active commit thread: force commit mode, no other override
        override="commit"
    fi
    # NOTE: frontier_gen scarcity / F-cascade triggers REMOVED.
    # frontier_gen still appears in the rotation for default-cycle
    # firing, but autonomy no longer forces it on A-grade scarcity —
    # that just produces motion-without-progress. The commit thread
    # mechanism replaces it.

    # Telemetry to stderr (won't pollute stdout return value)
    {
        echo "[autonomy] open_AV=$open_av  A_in_20=$a_count  last2=$last2  breakthrough_pending=$breakthrough_count  recent_meta=$recent_meta"
        echo "[autonomy] commit_thread='$commit_thread'  sessions_used=$commit_used  commit_active=$commit_active"
        if [ -n "$override" ]; then
            echo "[autonomy] OVERRIDE → $override"
        fi
    } 1>&2

    echo "$override"
}


# ============================================================
# MAIN LOOP
# ============================================================
while true; do
    # Per-run log files — each run gets its own pair so logs don't
    # accumulate into one giant file across the whole script invocation.
    RUN_TIMESTAMP=$(date +%Y%m%d_%H%M%S)
    LOGFILE="./archive/CLAUDE_OUTPUTS/claude_output_${RUN_TIMESTAMP}_run${RUN}.log"
    JSONFILE="./archive/CLAUDE_OUTPUTS/claude_output_${RUN_TIMESTAMP}_run${RUN}.json"

    MODE_IDX=$(( (RUN - 1) % ${#MODES[@]} ))
    MODE=${MODES[$MODE_IDX]}

    # Auto-pivot: state-based mode override
    OVERRIDE_MODE=$(compute_override)
    if [ -n "$OVERRIDE_MODE" ]; then
        echo "[autonomy] Default mode '$MODE' overridden by '$OVERRIDE_MODE'" | tee -a "$LOGFILE"
        MODE="$OVERRIDE_MODE"
    fi

    # Select prompt based on mode
    case $MODE in
        frontier)
            CURRENT_PROMPT="$PROMPT_FRONTIER"
            ;;
        cross_domain)
            CURRENT_PROMPT="$PROMPT_CROSSDOMAIN"
            ;;
        construction)
            CURRENT_PROMPT="$PROMPT_CONSTRUCTION"
            ;;
        arc)
            CURRENT_PROMPT="$PROMPT_ARC"
            ;;
        lean)
            CURRENT_PROMPT="$PROMPT_LEAN"
            ;;
        wild_swing)
            CURRENT_PROMPT="$PROMPT_WILDSWING"
            ;;
        novelty)
            CURRENT_PROMPT="$PROMPT_NOVELTY"
            ;;
        critique)
            CURRENT_PROMPT="$PROMPT_CRITIQUE"
            ;;
        frontier_gen)
            CURRENT_PROMPT="$PROMPT_FRONTIER_GEN"
            ;;
        verify)
            CURRENT_PROMPT="$PROMPT_VERIFY"
            ;;
        commit)
            CURRENT_PROMPT="$PROMPT_COMMIT"
            ;;
        reverify)
            CURRENT_PROMPT="$PROMPT_REVERIFY"
            ;;
        paradigm_shift)
            CURRENT_PROMPT="$PROMPT_PARADIGM_SHIFT"
            ;;
    esac

    # Tell the agent to persist the next run number before it finishes
    CURRENT_PROMPT="$CURRENT_PROMPT

# STATE: When you finish your work, run this command BEFORE your final message:
#   echo $((RUN + 1)) > .run_state
# This saves progress so the next ./run.sh invocation resumes at the correct run."

    echo ""
    echo "============================================================"
    echo "=== Run #$RUN — Mode: $MODE — $(date) ==="
    echo "    Log:  $LOGFILE"
    echo "    JSON: $JSONFILE"
    echo "============================================================"
    echo "=== Run #$RUN — Mode: $MODE — $(date) ===" | tee -a "$LOGFILE"
    echo '{"run":'"$RUN"',"mode":"'"$MODE"'","timestamp":"'"$(date -Iseconds)"'"}' >> "$JSONFILE"
    : > "$TMPFILE"
    : > "$ASSISTFILE"

    claude -p "$CURRENT_PROMPT" --output-format stream-json --verbose 2>&1 \
        | python3 -u -c "
import sys, json

json_f = open(sys.argv[1], 'a')
log_f  = open(sys.argv[2], 'a')
tmp_f  = open(sys.argv[3], 'w')
assist_f = open(sys.argv[4], 'w')

def out(text):
    print(text, flush=True)
    log_f.write(text + '\n')
    log_f.flush()

while True:
    line = sys.stdin.readline()
    if not line:
        break
    json_f.write(line)
    json_f.flush()
    tmp_f.write(line)
    tmp_f.flush()
    stripped = line.strip()
    if not stripped:
        continue
    try:
        obj = json.loads(stripped)
    except:
        out(stripped)
        continue
    if not isinstance(obj, dict):
        continue
    try:
        t = obj.get('type', '')
        if t == 'assistant':
            for block in obj.get('message', {}).get('content', []):
                if not isinstance(block, dict):
                    continue
                if block.get('type') == 'text':
                    out(block['text'])
                    assist_f.write(block['text'] + '\n')
                    assist_f.flush()
                elif block.get('type') == 'tool_use':
                    name = block.get('name', '')
                    inp = block.get('input', {})
                    cmd = inp.get('command', inp.get('file_path', inp.get('pattern', '')))
                    out(f'  [{name}] {cmd}')
        elif t == 'user':
            for block in obj.get('message', {}).get('content', []):
                if not isinstance(block, dict):
                    continue
                if block.get('type') == 'tool_result':
                    content = block.get('content', '')
                    if isinstance(content, str) and content:
                        lines = content.split('\n')
                        for l in lines[:30]:
                            out(f'    {l}')
                        if len(lines) > 30:
                            out(f'    ... ({len(lines) - 30} more lines)')
                    elif isinstance(content, list):
                        for item in content:
                            if isinstance(item, dict):
                                text = item.get('text', '')
                                if text:
                                    for l in text.split('\n')[:30]:
                                        out(f'    {l}')
        elif t == 'result':
            out('--- Session ended ---')
    except Exception:
        pass

json_f.close()
log_f.close()
tmp_f.close()
assist_f.close()
" "$JSONFILE" "$LOGFILE" "$TMPFILE" "$ASSISTFILE"

    # If a session generated freeform proposals and saved them, surface them
    # for the next critique run. (Most production modes won't, but if any
    # session leaves novel proposals in proposals_session.md, copy to latest.)
    if [ -f "./archive/ephemeral/proposals_session.md" ]; then
        if [ ! -f "$PROPOSALS_FILE" ] || [ "./archive/ephemeral/proposals_session.md" -nt "$PROPOSALS_FILE" ]; then
            cp "./archive/ephemeral/proposals_session.md" "$PROPOSALS_FILE"
        fi
    fi

    # Breakthrough handling: requires TWO consecutive verifications.
    # - "I FOUND IT!!!" in a non-verify session starts the chain (counter=0)
    # - Each subsequent verify session that ALSO contains "I FOUND IT!!!"
    #   increments counter
    # - When counter reaches 2, halt the loop with a banner file
    # - A verify session WITHOUT "I FOUND IT!!!" resets counter, continues
    if grep -qF 'I FOUND IT!!!' "$ASSISTFILE"; then
        if [ "$MODE" = "verify" ]; then
            CURRENT_COUNT=0
            [ -f "$BREAKTHROUGH_FILE" ] && CURRENT_COUNT=$(cat "$BREAKTHROUGH_FILE" 2>/dev/null || echo 0)
            CURRENT_COUNT=$((CURRENT_COUNT + 1))
            echo "$CURRENT_COUNT" > "$BREAKTHROUGH_FILE"
            echo "[breakthrough] Verification confirmed (count=$CURRENT_COUNT/2)" | tee -a "$LOGFILE"
            if [ "$CURRENT_COUNT" -ge 2 ]; then
                echo "" | tee -a "$LOGFILE"
                echo "============================================================" | tee -a "$LOGFILE"
                echo "=== BREAKTHROUGH VERIFIED (2/2) — Run #$RUN ===" | tee -a "$LOGFILE"
                echo "============================================================" | tee -a "$LOGFILE"
                {
                    echo "# BREAKTHROUGH DETECTED AND DOUBLE-VERIFIED"
                    echo ""
                    echo "**Detection run:** see counter file history"
                    echo "**Final verification run:** $RUN ($(date -Iseconds))"
                    echo "**Mode at final verify:** $MODE"
                    echo ""
                    echo "## What was verified"
                    echo "See \`.verify_target\` for the path of the synthesis file"
                    echo "containing the original breakthrough claim."
                    echo ""
                    echo "Latest verify_target: $(cat "$VERIFY_TARGET_FILE" 2>/dev/null)"
                    echo ""
                    echo "## Next steps"
                    echo "1. Read the original breakthrough session synthesis."
                    echo "2. Read both verify session syntheses (most recent two)."
                    echo "3. If genuine: write up novel/breakthrough.md and announce."
                    echo "4. If still suspect: invoke a third independent verify"
                    echo "   manually before publishing."
                } > "./BREAKTHROUGH.md"
                echo "Wrote BREAKTHROUGH.md banner. Halting loop." | tee -a "$LOGFILE"
                break
            fi
        else
            # Fresh "I FOUND IT!!!" outside verify mode — start chain
            echo 0 > "$BREAKTHROUGH_FILE"
            echo "[breakthrough] Claim detected in mode '$MODE'. Forcing 2-step verification." | tee -a "$LOGFILE"
        fi
    else
        # No "I FOUND IT!!!" in this session
        if [ "$MODE" = "verify" ] && [ -f "$BREAKTHROUGH_FILE" ]; then
            echo "[breakthrough] Verify session did NOT confirm. Resetting counter." | tee -a "$LOGFILE"
            rm -f "$BREAKTHROUGH_FILE"
        fi
    fi

    RUN=$((RUN + 1))
    echo "$RUN" > "$STATE_FILE"
    echo "Sleeping 60 seconds before next run..." | tee -a "$LOGFILE"
    sleep 60
done
