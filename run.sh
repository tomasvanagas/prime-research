#!/bin/bash
# ============================================================
# PRIME RESEARCH — CONTINUOUS LOOP (v2, 2026-06-13)
# ============================================================
# One mode. Each cycle is a fresh Claude session (Fable 5, xhigh)
# working directly toward the goal the way the S491 cycles did:
# read PROGRAM.md, take the next action (or justify a better one),
# build code that runs, measure, report honestly, update PROGRAM.md.
#
# The v1 multi-stage rotation framework (modes, grading, commit
# threads, verify chains) is retired: archive/framework_v1/.
#
# The loop halts only when BREAKTHROUGH.md exists — which a cycle may
# create ONLY with machine-verified, reproducible evidence of the
# actual goal.
# ============================================================

set -u
cd "$(dirname "$0")"

MODEL="claude-opus-4-8[1m]"
EFFORT="xhigh"
OUTDIR="archive/CLAUDE_OUTPUTS"
CYCLE_FILE=".cycle"
SLEEP_BETWEEN=600

mkdir -p "$OUTDIR"
[ -f "$CYCLE_FILE" ] || echo 0 > "$CYCLE_FILE"

PROMPT=$(cat <<'ENDPROMPT'
# RESEARCH CYCLE — work directly toward the goal

You are one cycle of a continuous research loop on computing the n-th
prime in O(polylog) time (project root: this directory). There are no
modes, no grades, no rotation. Work plainly and honestly.

## Start here (briefly)
1. Read PROGRAM.md — the live program. It records what is done, what
   is open, and ONE concrete next action.
2. If continuing the active line, read the doc it points to
   (novel/succinct_verification_of_pi.md) and the relevant
   experiments/ results files.
3. Before calling anything new, search status/CLOSED_PATHS.md (730+
   closed routes) and EDGES.md. CLAUDE.md has the conventions.

## The contract (how the S491 cycles worked — match them)
- Pick ONE item: the recorded NEXT ACTION, or a better idea you can
  explicitly justify against PROGRAM.md and the catalogue. Commit to
  it for the whole cycle. Multi-cycle builds are normal: record the
  design and partial state in PROGRAM.md so the next cycle continues
  mid-build instead of restarting.
- Build code that runs: experiments/<topic>/<name>.py (one script,
  CLI-parameterised, with --selftest) + <name>_results.md. Every
  boundary case you debug goes into the selftest. Measure; use
  matched baselines/controls where they mean something. Every results
  file states what would falsify the result.
- Honest reporting is absolute. Failures are filed as failures with
  the structural mechanism. Your own (or earlier cycles') claims get
  corrected when wrong — that is progress, not embarrassment. A clean
  negative with a mechanism beats a vague positive. Never inflate.
- Do not redo what PROGRAM.md records as done. Do not re-derive
  closed paths. Do not rebuild working artifacts.
- The goal may be unreachable with current mathematics. The
  deliverable is real progress: working algorithms/protocols, new
  structural facts with measurements, precisely-closed questions,
  adjacent-problem wins (verification, batching, certification).

## Before you stop
- Update PROGRAM.md: one entry for what this cycle produced, the
  current open list, and a single concrete NEXT ACTION.
- Every new .py has its _results.md. Remove __pycache__ directories.
- BREAKTHROUGH.md: write it ONLY if you have machine-verified,
  reproducible evidence of the actual goal (exact pi(x)/p(n) at large
  x in polylog time — verified output equality, measured scaling,
  reproduction instructions). The loop halts on its existence. Never
  write it on heuristics, approximations, or conditional results.
ENDPROMPT
)

echo "Continuous research loop v2 — model=$MODEL effort=$EFFORT"
echo "Halt condition: BREAKTHROUGH.md exists. Ctrl-C to stop manually."

while true; do
    if [ -f BREAKTHROUGH.md ]; then
        echo "BREAKTHROUGH.md present — halting loop. Read it, then the"
        echo "cycle log in $OUTDIR that produced it."
        exit 0
    fi

    N=$(($(cat "$CYCLE_FILE") + 1))
    echo "$N" > "$CYCLE_FILE"
    TS=$(date +%Y%m%d_%H%M%S)
    LOG="$OUTDIR/cycle_${TS}_n${N}.log"
    JSON="$OUTDIR/cycle_${TS}_n${N}.json"

    echo ""
    echo "============================================================"
    echo "=== Cycle $N — $(date) ==="
    echo "    Log:  $LOG"
    echo "============================================================"
    echo "=== Cycle $N — $(date) ===" >> "$LOG"

    claude -p "$PROMPT" \
        --model "$MODEL" \
        --effort "$EFFORT" \
        --output-format stream-json --verbose 2>&1 \
        | python3 -u -c "
import sys, json
json_f = open(sys.argv[1], 'a')
log_f = open(sys.argv[2], 'a')
def out(text):
    print(text, flush=True)
    log_f.write(text + '\n'); log_f.flush()
for line in sys.stdin:
    json_f.write(line); json_f.flush()
    s = line.strip()
    if not s:
        continue
    try:
        obj = json.loads(s)
    except Exception:
        out(s); continue
    if not isinstance(obj, dict):
        continue
    t = obj.get('type', '')
    if t == 'assistant':
        for blk in obj.get('message', {}).get('content', []):
            if blk.get('type') == 'text':
                out(blk['text'])
            elif blk.get('type') == 'tool_use':
                name = blk.get('name', '?')
                inp = blk.get('input', {})
                hint = inp.get('command') or inp.get('file_path') or ''
                out(f'  [tool: {name}] {str(hint)[:120]}')
    elif t == 'result':
        out(f\"--- cycle result: {obj.get('subtype','')} \"
            f\"(turns={obj.get('num_turns','?')}, \"
            f\"duration={obj.get('duration_ms','?')}ms) ---\")
" "$JSON" "$LOG"

    echo "=== Cycle $N finished — $(date) ===" >> "$LOG"
    sleep "$SLEEP_BETWEEN"
done
