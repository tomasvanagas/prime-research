#!/bin/bash

# ============================================================
# PRIME RESEARCH: NOVELTY-FIRST AUTONOMOUS LOOP (S70+ rewrite)
# ============================================================
# Cycles through 5 modes optimised for *genuine novelty production*
# rather than disciplined-duplicate-closure:
#
#   1. novelty       - Pick a NOVELTY_CHALLENGES.md §2 / §4 / §5 target
#   2. construction  - Pick a NOVELTY_CHALLENGES.md §1 composition target,
#                      build it under experiments/constructions/
#   3. arc           - Continue a RESEARCH_AGENDA.md multi-session arc
#   4. lean          - Work on a NOVELTY_CHALLENGES.md §3 Lean formalisation
#   5. critique      - Verify recent work, enforce novelty bar
#
# Replaces the prior {normal, wildcard, focused, propose, critique}
# rotation. The "no-anchor" propose and wildcard prompts produced ~70%
# duplicate-shape outputs at the project's current saturation (per S60
# critic's diagnosis), and `focused` mode rotated through a stale
# FOCUS_QUEUE.md (all tasks marked COMPLETED). All four are now retargeted
# to read NOVELTY_CHALLENGES.md and RESEARCH_AGENDA.md instead.
#
# Critique mode is preserved largely as-is — it was working correctly.
# ============================================================


# ============================================================
# PROMPT: Novelty Mode
# ============================================================
PROMPT_NOVELTY=$(cat << 'ENDPROMPT'
# NOVELTY SESSION

Mission: produce ONE new mathematical artefact for the project this session.
"Artefact" = an object, identity, proof, refinement, or composition that did
not exist in the project before this session.

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

5. Identify the single highest-value next-action and write it into
   NOVELTY_CHALLENGES.md or RESEARCH_AGENDA.md as appropriate.

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
# Rotation order: novelty -> construction -> arc -> lean -> critique
# The first four are content-producing modes; critique enforces the
# novelty bar on what they produced. Net cycle: 4 production + 1 audit.
MODES=("novelty" "construction" "arc" "lean" "critique")


# ============================================================
# STATE PERSISTENCE
# ============================================================
STATE_FILE="./.run_state"
if [ -f "$STATE_FILE" ]; then
    RUN=$(cat "$STATE_FILE")
    echo "Resuming from run #$RUN (loaded from $STATE_FILE)"
else
    RUN=1
fi


# ============================================================
# MAIN LOOP
# ============================================================
while true; do
    MODE_IDX=$(( (RUN - 1) % ${#MODES[@]} ))
    MODE=${MODES[$MODE_IDX]}

    # Select prompt based on mode
    case $MODE in
        novelty)
            CURRENT_PROMPT="$PROMPT_NOVELTY"
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
        critique)
            CURRENT_PROMPT="$PROMPT_CRITIQUE"
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

    # Check for breakthrough
    if grep -qF 'I FOUND IT!!!' "$ASSISTFILE"; then
        echo ""
        echo "============================================================"
        echo "=== BREAKTHROUGH DETECTED — Run #$RUN — Mode: $MODE ==="
        echo "============================================================"
        echo "Detected 'I FOUND IT!!!' — stopping." | tee -a "$LOGFILE"
        break
    fi

    RUN=$((RUN + 1))
    echo "$RUN" > "$STATE_FILE"
    echo "Sleeping 5 seconds before next run..." | tee -a "$LOGFILE"
    sleep 5
done
