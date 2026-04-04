# Prime Research

Autonomous AI research project to find an O(polylog) algorithm for computing
the nth prime number p(n) exactly. Target: p(10^100) in <1 second.

<br/>

## Status

**500+ approaches tested** across 18 sessions with 130+ AI sub-agents.
All known paths are closed, but no proof of impossibility exists.
The problem remains genuinely open.

**Best result:** `algorithms/v10_c_accelerated.py` -- p(10^9) in 0.175s, 100% exact, O(x^{2/3}).

<br/>

## Running

```bash
./run.sh
```

Launches a mode-rotating autonomous loop that cycles through 5 research strategies:
normal exploration, wildcard (fresh thinking), focused deep-dives, proposal generation,
and adversarial critique. Logs go to `archive/CLAUDE_OUTPUTS/`.

<br/>

## Project Structure

```
CLAUDE.md              -- AI agent entry point (concise, read this first)
FOCUS_QUEUE.md         -- Deep-dive tasks for focused research sessions
run.sh                 -- Mode-rotating autonomous research loop
status/
  CLOSED_PATHS.md      -- 500+ tested approaches (searchable lookup table)
  OPEN_PROBLEMS.md     -- Viable research directions
  BEST_ALGORITHMS.md   -- Working implementations with benchmarks
  SESSION_INSIGHTS.md  -- Per-session detailed findings (Sessions 12-18)
proven/
  barriers.md          -- Proven impossibility results
  complexity.md        -- Upper/lower bounds, circuit complexity
  information.md       -- Entropy, Kolmogorov, entanglement
  quantum.md           -- Why quantum doesn't help
novel/
  info_computation_gap.md  -- Best original insight
  failure_taxonomy.md      -- Three failure modes classification
  session16_synthesis.md   -- 15 families closed, three "pillars"
  session17_synthesis.md   -- Exact communication rank formula
algorithms/            -- Working, tested code only
literature/
  references.md        -- All citations
  state_of_art_2026.md -- Latest published results
experiments/           -- Past experiments by topic
  analytic/ algebraic/ quantum/ ml/ information_theory/
  dynamical/ topological/ sieve/ circuit_complexity/ other/
  wildcard/ proposals/
data/                  -- Zeta zeros (200/300/500/1000)
archive/               -- Session logs and visualizations
```

<br/>

## Key Insight

The correction delta(n) = p(n) - R^{-1}(n) has only O(log n) bits of information,
but extracting those bits requires O(x^{2/3}) computation. No unconditional lower
bound beyond Omega(log x) exists. The circuit complexity of pi(x) is unstudied.

<br/>

## References

See `literature/references.md` for the full bibliography.
Key papers: Aggarwal (2025), Guth-Maynard (2024), Latorre-Sierra (2014),
Kolpakov-Rocke (2023), Deleglise-Rivat (1996), Lagarias-Odlyzko (1987).
