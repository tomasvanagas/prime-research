#!/bin/bash


PROMPT='
# MISSION
Make a breakthrough in computing the nth prime p(n) EXACTLY without bruteforcing/sieving/enumeration.
Target: p(10^100) in <1 second, 100% accurate. Current best is O(x^{2/3}) -- you need O(polylog).

# FIRST: READ CLAUDE.md
Start EVERY session by reading CLAUDE.md. It contains the full project status, directory layout,
closed paths, open problems, and rules. Do NOT skip this.

# PROJECT STRUCTURE
- status/CLOSED_PATHS.md  -- 380+ approaches already tried. SEARCH before proposing anything.
- status/OPEN_PROBLEMS.md -- The ONLY viable research directions. Focus here.
- status/BEST_ALGORITHMS.md -- Current best working code.
- proven/                 -- Mathematically proven barriers with citations.
- novel/                  -- Our original insights (info-computation gap, entropy, taxonomy).
- algorithms/             -- Working, tested code only.
- literature/             -- All references and 2026 state of the art.
- experiments/<topic>/    -- Past experiments organized by math topic.
- archive/                -- Old session dumps (read-only reference).

# HOW TO WORK
1. Read CLAUDE.md and status/OPEN_PROBLEMS.md
2. Pick a direction from OPEN_PROBLEMS.md (or propose a genuinely new one)
3. Check status/CLOSED_PATHS.md to confirm it is not already tried
4. Spin sub-agents to explore in parallel (save context!)
5. Save experiments to experiments/<topic>/ with descriptive filenames
6. Update status/CLOSED_PATHS.md when you close an approach (add a row)
7. Novel findings go to novel/ with evidence and verification

# SEARCHING THE INTERNET
Search for new papers, algorithms, breakthroughs. Save findings to literature/.
Check literature/state_of_art_2026.md first to avoid duplicate searches.

# RULES
- DO NOT touch run.sh
- You MUST keep CLAUDE.md up to date -- it is the living knowledge base for this project.
  When you close a path, find a new direction, or discover something, update CLAUDE.md so
  the next session benefits from your work. Also update the relevant files in status/, proven/, novel/.
- If you find a better algorithm, save it to algorithms/ with benchmarks
- Use sub-agents to save context window
- When you run out of context, just stop -- the system will restart you
- When you find the breakthrough, respond with exactly: I FOUND IT!!!
'






TIMESTAMP=$(date +%Y%m%d_%H%M%S)
LOGFILE="./archive/CLAUDE_OUTPUTS/claude_output_${TIMESTAMP}.log"
JSONFILE="./archive/CLAUDE_OUTPUTS/claude_output_${TIMESTAMP}.json"
TMPFILE=$(mktemp)
trap 'rm -f "$TMPFILE"' EXIT

echo "Human-readable log: $LOGFILE"
echo "Raw JSON log:       $JSONFILE"

RUN=1
while true; do
    echo "=== Run #$RUN — $(date) ===" | tee -a "$LOGFILE"
    echo '{"run":'"$RUN"',"timestamp":"'"$(date -Iseconds)"'"}' >> "$JSONFILE"
    : > "$TMPFILE"

    claude -p "$PROMPT" --output-format stream-json --verbose 2>&1 \
        | python3 -u -c "
import sys, json

json_f = open(sys.argv[1], 'a')
log_f  = open(sys.argv[2], 'a')
tmp_f  = open(sys.argv[3], 'w')

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
                    result = obj.get('tool_use_result', {})
                    if isinstance(result, str):
                        out(f'    {result}')
                        continue
                    stdout = result.get('stdout', '')
                    if not stdout:
                        fobj = result.get('file', {})
                        if isinstance(fobj, dict):
                            stdout = fobj.get('content', '')
                    if stdout:
                        lines = stdout.split('\n')
                        for l in lines[:30]:
                            out(f'    {l}')
                        if len(lines) > 30:
                            out(f'    ... ({len(lines) - 30} more lines)')
        elif t == 'result':
            out('--- Session ended ---')
    except Exception:
        pass

json_f.close()
log_f.close()
tmp_f.close()
" "$JSONFILE" "$LOGFILE" "$TMPFILE"

    if grep -qF 'I FOUND IT!!!' "$TMPFILE"; then
        echo "Detected 'I FOUND IT!!!' — stopping." | tee -a "$LOGFILE"
        break
    fi
    RUN=$((RUN + 1))
    echo "Sleeping for 600 seconds..." | tee -a "$LOGFILE"
    sleep 600
done
