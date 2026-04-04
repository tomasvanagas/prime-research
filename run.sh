#!/bin/bash

# ============================================================
# PRIME RESEARCH: MODE-ROTATING AUTONOMOUS LOOP
# ============================================================
# Cycles through 5 research modes to break the plateau:
#   1. normal    - Read CLAUDE.md, explore open problems
#   2. wildcard  - Fresh thinking, NO closed paths, first principles
#   3. focused   - Deep dive on one task from FOCUS_QUEUE.md
#   4. propose   - Generate novel approaches (without seeing failures)
#   5. critique  - Evaluate proposals against 500+ closed paths
# ============================================================


# ============================================================
# PROMPT: Normal Mode
# ============================================================
PROMPT_NORMAL=$(cat << 'ENDPROMPT'
# MISSION
Make a breakthrough in computing the nth prime p(n) EXACTLY without bruteforcing/sieving/enumeration.
Target: p(10^100) in <1 second, 100% accurate. Current best is O(x^{2/3}) -- you need O(polylog).

# FIRST: READ CLAUDE.md
Start by reading CLAUDE.md. It contains the project status, directory layout,
and rules. It points to status/ files for closed paths and open problems.

# PROJECT STRUCTURE
- status/CLOSED_PATHS.md  -- 500+ approaches already tried. SEARCH before proposing anything.
- status/OPEN_PROBLEMS.md -- The ONLY viable research directions. Focus here.
- status/BEST_ALGORITHMS.md -- Current best working code.
- proven/                 -- Mathematically proven barriers with citations.
- novel/                  -- Our original insights.
- algorithms/             -- Working, tested code only.
- literature/             -- All references and 2026 state of the art.
- experiments/<topic>/    -- Past experiments organized by topic.
- data/                   -- Zeta zeros for explicit formula work.

# HOW TO WORK
1. Read CLAUDE.md and status/OPEN_PROBLEMS.md
2. Pick a direction from OPEN_PROBLEMS.md (or propose a genuinely new one)
3. Check status/CLOSED_PATHS.md to confirm it has not been tried
4. Spin sub-agents to explore in parallel (save context!)
5. Save experiments to experiments/<topic>/
6. Update status/CLOSED_PATHS.md when you close an approach
7. Novel findings go to novel/ with evidence

# ALSO CHECK
- Read novel/critique_latest.md if it exists -- it may contain NOVEL proposals
  that survived critique and are worth implementing.

# RULES
- DO NOT touch run.sh or FOCUS_QUEUE.md
- Update CLAUDE.md only for significant status changes (new best algorithm, major barrier)
- Use sub-agents to save context window
- When you find the breakthrough, respond with exactly: I FOUND IT!!!
ENDPROMPT
)


# ============================================================
# PROMPT: Wildcard Mode (fresh thinking, no anchoring)
# ============================================================
PROMPT_WILDCARD=$(cat << 'ENDPROMPT'
# FRESH PERSPECTIVE SESSION

You are attacking one of the great open problems in computational number theory:
computing p(n) -- the nth prime number -- EXACTLY in O(polylog(n)) time.

Current best: O(x^{2/3}) via Meissel-Lehmer combinatorial sieve.
Target: p(10^100) in <1 second, 100% accuracy.

The smooth approximation R^{-1}(n) gives ~50% of digits in O(polylog) time.
The remaining ~50% encode oscillatory contributions from Riemann zeta zeros.

## YOUR MISSION
Do NOT start by reading CLAUDE.md or status/CLOSED_PATHS.md.
This is a FRESH THINKING session. Start from first principles.

## APPROACH
1. Think about analogies from other fields where "impossible" barriers were broken:
   - Shor's algorithm broke factoring via quantum Fourier transform
   - Compressed sensing broke Nyquist via sparsity assumptions
   - AlphaFold broke protein folding via learned energy landscapes
   - Fast multipole method broke N-body from O(N^2) to O(N) via hierarchy
   - Candes-Tao broke matrix completion via nuclear norm relaxation
   What is the analogous structural insight for prime counting?

2. Brainstorm at least 5 genuinely unconventional ideas. Spirit examples:
   - Can primes be characterized by a DYNAMICAL SYSTEM with fast-forwardable orbits?
   - Is there a PROBABILISTIC identity that gives exact answers?
   - Can pi(x) be encoded as a LINEAR ALGEBRA problem over a clever ring/field?
   - Can AUTOMATED THEOREM PROVING search for novel identities?
   - Is there a HIERARCHICAL DECOMPOSITION of the sieve that enables recursion?
   - Can the zeta zero sum be replaced by a SPECTRAL SHORTCUT?
   - Does the ADELIC perspective give a fast algorithm via local-global?

3. For each idea: write pseudocode, test on small cases (n < 10000), analyze complexity.

4. Search the internet for latest 2025-2026 papers. Look in unexpected places:
   quantum information, algebraic topology, machine learning theory, cryptography.

## SAVING
- Save experiments to experiments/wildcard/
- Promising findings go to novel/wildcard_findings.md
- Better algorithms go to algorithms/
- DO NOT read CLAUDE.md or status/CLOSED_PATHS.md
- DO NOT spend time proving things impossible
- When you find the breakthrough, respond with exactly: I FOUND IT!!!
ENDPROMPT
)


# ============================================================
# PROMPT: Focused Mode (deep dive on one task)
# ============================================================
PROMPT_FOCUSED=$(cat << 'ENDPROMPT'
# DEEP FOCUS SESSION — Task #__TASK_NUM__

This is a DEEP FOCUS session. Work on exactly ONE research direction
for the entire session. No broad exploration, no topic switching.

## YOUR TASK
Read FOCUS_QUEUE.md in the project root. Work on Task #__TASK_NUM__.
Follow the task description: run every experiment listed, write code,
test hypotheses, and record results.

## CONTEXT
- Read CLAUDE.md for project status and existing knowledge
- Read status/OPEN_PROBLEMS.md for background on your task
- Check experiments/ for any prior work on this topic

## RULES
- Stay focused on Task #__TASK_NUM__. Do NOT explore other directions.
- Write code and run experiments. This is NOT a theory-only session.
- Save all experiments to the directory specified in FOCUS_QUEUE.md.
- Significant discoveries: update CLAUDE.md and novel/
- Closed directions: update status/CLOSED_PATHS.md with evidence
- Use sub-agents for parallel experiments
- DO NOT touch run.sh or FOCUS_QUEUE.md
- When you find the breakthrough, respond with exactly: I FOUND IT!!!
ENDPROMPT
)


# ============================================================
# PROMPT: Propose Mode (generate ideas without seeing failures)
# ============================================================
PROMPT_PROPOSE=$(cat << 'ENDPROMPT'
# PROPOSAL SESSION

You are a mathematician attacking the problem of computing p(n) -- the nth prime --
exactly in O(polylog(n)) time.

## BACKGROUND
- Prime number theorem: p(n) ~ n*ln(n)
- R^{-1}(n) gives ~50% of digits in O(polylog) time
- Riemann explicit formula: pi(x) = R(x) - sum_rho R(x^rho) - ...
  Requires O(sqrt(x)) zeta zeros for exact results
- Meissel-Lehmer: O(x^{2/3}) via inclusion-exclusion on floor(x/k) values
- The correction delta(n) = p(n) - R^{-1}(n) has O(log n) bits but costs O(x^{2/3}) to compute

## YOUR MISSION
Propose AT LEAST 3 concrete approaches. For each:
1. State the mathematical idea clearly
2. Write pseudocode (runnable Python preferred)
3. Analyze the time complexity
4. Identify the key assumption or conjecture it relies on
5. Design a computational test for n < 10000

IMPORTANT: Do NOT read CLAUDE.md or status/CLOSED_PATHS.md for this session.
We want fresh ideas unconstrained by prior analysis.

Search the internet for latest papers, tools, and techniques that might help.
Look in unexpected fields: quantum computing, cryptography, compressed sensing,
algebraic topology, representation theory, information theory.

## CREATIVE DIRECTIONS TO CONSIDER
- Algebraic geometry / etale cohomology shortcuts for counting
- Quantum-inspired classical algorithms (dequantization)
- Machine learning on the RESIDUAL delta(n), not the primes
- Automated conjecture generation via PSLQ / LLL / OEIS
- Connections to other counting problems with known fast algorithms
- Modular forms / automorphic L-functions
- Operator-theoretic approaches (trace formula shortcuts)
- Information-theoretic: prove O(polylog) is achievable via entropy arguments

## SAVING
- Save detailed proposals to novel/proposals_session.md
- If any idea works on small cases, save code to experiments/proposals/
- When you find the breakthrough, respond with exactly: I FOUND IT!!!
ENDPROMPT
)


# ============================================================
# PROMPT: Critique Mode (evaluate proposals rigorously)
# ============================================================
PROMPT_CRITIQUE=$(cat << 'ENDPROMPT'
# CRITIQUE SESSION

You are a rigorous mathematical critic evaluating proposed approaches
to computing p(n) in O(polylog(n)) time.

## YOUR TASK
1. Read archive/proposals_latest.md for the latest proposals.
   If missing/empty, try novel/proposals_session.md instead.
   If neither exists, read CLAUDE.md and status/OPEN_PROBLEMS.md for context,
   then search the internet for any 2025-2026 papers proposing new prime algorithms.

2. Read status/CLOSED_PATHS.md THOROUGHLY. It contains 500+ tested approaches.

3. For EACH proposal, evaluate:
   a) Has this exact approach been tried? (check CLOSED_PATHS.md)
   b) Does it fall into a known failure mode?
      - CIRCULARITY: needs primes to compute primes
      - EQUIVALENCE: reduces to explicit formula / zeta zero sum
      - INFORMATION LOSS: smooth approximation loses ~170 bits
   c) Is the complexity analysis correct?
   d) What is the specific mathematical obstacle?
   e) VERDICT: NOVEL / DUPLICATE / FLAWED

4. For NOVEL proposals: design a concrete experiment to test the core hypothesis.
   Use sub-agents to actually RUN the experiment and report results.

5. For PARTIALLY novel proposals: extract the novel component and test it separately.

## SAVING
- Save your critique to novel/critique_latest.md
- NOVEL proposals: add to status/OPEN_PROBLEMS.md
- DUPLICATE proposals: note which CLOSED_PATHS entry they match
- Update CLAUDE.md if your analysis reveals new insights
- DO NOT touch run.sh or FOCUS_QUEUE.md
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
PROPOSALS_FILE="./archive/proposals_latest.md"
trap 'rm -f "$TMPFILE" "$ASSISTFILE"' EXIT

echo "Human-readable log: $LOGFILE"
echo "Raw JSON log:       $JSONFILE"


# ============================================================
# MODE ROTATION CONFIG
# ============================================================
MODES=("normal" "wildcard" "focused" "propose" "critique")
NUM_FOCUS_TASKS=4


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

    # For focused mode: compute which task to work on (cycles 1-4)
    CYCLE=$(( (RUN - 1) / ${#MODES[@]} ))
    TASK_NUM=$(( (CYCLE % NUM_FOCUS_TASKS) + 1 ))

    # Select prompt based on mode
    case $MODE in
        normal)
            CURRENT_PROMPT="$PROMPT_NORMAL"
            ;;
        wildcard)
            CURRENT_PROMPT="$PROMPT_WILDCARD"
            ;;
        focused)
            CURRENT_PROMPT="${PROMPT_FOCUSED//__TASK_NUM__/$TASK_NUM}"
            ;;
        propose)
            CURRENT_PROMPT="$PROMPT_PROPOSE"
            ;;
        critique)
            CURRENT_PROMPT="$PROMPT_CRITIQUE"
            ;;
    esac

    echo ""
    echo "============================================================"
    echo "=== Run #$RUN — Mode: $MODE — $(date) ==="
    if [ "$MODE" = "focused" ]; then
        echo "=== Focus Task: #$TASK_NUM ==="
    fi
    echo "============================================================"
    echo "" | tee -a "$LOGFILE"
    echo "=== Run #$RUN — Mode: $MODE — $(date) ===" | tee -a "$LOGFILE"
    echo '{"run":'"$RUN"',"mode":"'"$MODE"'","task_num":'"$TASK_NUM"',"timestamp":"'"$(date -Iseconds)"'"}' >> "$JSONFILE"
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

    # For propose mode: save proposals for the next critique run
    if [ "$MODE" = "propose" ]; then
        cp "$ASSISTFILE" "$PROPOSALS_FILE"
        echo "Proposals saved to: $PROPOSALS_FILE" | tee -a "$LOGFILE"
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
    echo "Sleeping 600 seconds before next run..." | tee -a "$LOGFILE"
    sleep 600
done
