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
# INFRASTRUCTURE
# ============================================================
TIMESTAMP=$(date +%Y%m%d_%H%M%S)
LOGFILE="./archive/CLAUDE_OUTPUTS/claude_output_${TIMESTAMP}.log"
JSONFILE="./archive/CLAUDE_OUTPUTS/claude_output_${TIMESTAMP}.json"
TMPFILE=$(mktemp)
ASSISTFILE=$(mktemp)
PROPOSALS_FILE="./archive/ephemeral/proposals_latest.md"
CRITIQUE_FILE="./archive/ephemeral/critique_latest.md"
trap 'rm -f "$TMPFILE" "$ASSISTFILE"' EXIT

echo "Human-readable log: $LOGFILE"
echo "Raw JSON log:       $JSONFILE"


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
MODES=("frontier" "cross_domain" "construction" "arc" "lean" "wild_swing" "novelty" "critique")


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

    # Decision tree (priority order)
    # Breakthrough verifications and routine A-grade verify always fire
    # — these are ungated by cooldown because correctness > cooldown.
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
    # frontier_gen triggers gated by cooldown to avoid flapping.
    # If meta sessions (frontier_gen / verify) have fired in last 5,
    # prefer to let production sessions attack the resulting state.
    elif [ "$open_av" -lt 4 ]; then
        # Frontier exhaustion (always fire — frontier truly bare)
        override="frontier_gen"
    elif [ "$f_cascade" = "1" ] && [ "$recent_meta" -eq 0 ]; then
        # Two F's in a row, no recent meta intervention
        override="frontier_gen"
    elif [ "$a_count" -eq 0 ] && [ "$RUN" -gt 20 ] && [ "$recent_meta" -eq 0 ]; then
        # 20-session A-grade scarcity, no recent intervention
        override="frontier_gen"
    fi

    # Telemetry to stderr (won't pollute stdout return value)
    {
        echo "[autonomy] open_AV=$open_av  A_in_20=$a_count  last2=$last2  breakthrough_pending=$breakthrough_count  recent_meta=$recent_meta"
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
    esac

    # Tell the agent to persist the next run number before it finishes
    CURRENT_PROMPT="$CURRENT_PROMPT

# STATE: When you finish your work, run this command BEFORE your final message:
#   echo $((RUN + 1)) > .run_state
# This saves progress so the next ./run.sh invocation resumes at the correct run."

    echo ""
    echo "============================================================"
    echo "=== Run #$RUN — Mode: $MODE — $(date) ==="
    echo "============================================================"
    echo "" | tee -a "$LOGFILE"
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
    echo "Sleeping 5 seconds before next run..." | tee -a "$LOGFILE"
    sleep 5
done
