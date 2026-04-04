# Prime Research

Autonomous AI research project to find an O(polylog) algorithm for computing
the nth prime number p(n) exactly. Target: p(10^100) in <1 second.

<br/>

## Status

**380+ approaches tested** across 10 sessions with 89+ AI sub-agents.
All known paths are closed, but no proof of impossibility exists.
The problem remains genuinely open.

**Best result:** `algorithms/v10_c_accelerated.py` -- p(10^9) in 0.175s, 100% exact, O(x^{2/3}).

<br/>

## Running

```bash
./run.sh
```

Launches an autonomous loop: Claude reads CLAUDE.md, explores open problems,
spins sub-agents, saves results, and restarts when context runs out.
Logs go to `archive/CLAUDE_OUTPUTS/`.

<br/>

## Project Structure

```
CLAUDE.md              -- AI agent entry point (read this first)
run.sh                 -- Autonomous research loop
status/
  CLOSED_PATHS.md      -- 380+ tested approaches (searchable lookup table)
  OPEN_PROBLEMS.md     -- Viable research directions
  BEST_ALGORITHMS.md   -- Working implementations with benchmarks
proven/
  barriers.md          -- Proven impossibility results
  complexity.md        -- Upper/lower bounds, circuit complexity
  information.md       -- Entropy, Kolmogorov, entanglement
  quantum.md           -- Why quantum doesn't help
novel/
  info_computation_gap.md  -- Best original insight
  entropy_measurements.md  -- Empirical measurements
  failure_taxonomy.md      -- Three failure modes classification
algorithms/            -- Working, tested code only
literature/
  references.md        -- All citations
  state_of_art_2026.md -- Latest published results
experiments/           -- Past experiments by topic
  analytic/ algebraic/ quantum/ ml/ information_theory/
  dynamical/ topological/ sieve/ other/
archive/               -- Session history (read-only)
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
