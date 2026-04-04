# Prime Research: Computing p(n) Exactly Without Bruteforcing

## Project Goal
Find an efficient algorithm to compute the nth prime number p(n) without
enumeration, sieving, or brute force. Target: p(10^100) in <1 second, 100% exact.

## Status (April 2026)
- **380+ approaches tested** across 10 sessions, 89+ sub-agents
- **ALL KNOWN PATHS CLOSED** but no proof that polylog is impossible
- **Problem is GENUINELY OPEN** -- no unconditional lower bound beyond Omega(log x)
- Best exact: `algorithms/v10_c_accelerated.py` -- O(p(n)^{2/3}), p(10^9) in 0.175s
- Best approximate: R^{-1}(n) -- O(polylog), ~47% digits correct

## Directory Layout

```
status/
  CLOSED_PATHS.md    <-- SEARCH HERE before proposing ANY approach (380+ entries)
  OPEN_PROBLEMS.md   <-- The ONLY viable research directions
  BEST_ALGORITHMS.md <-- Working implementations with benchmarks
proven/
  barriers.md        <-- Mathematically proven impossibility results
  complexity.md      <-- Upper/lower bounds, circuit complexity
  information.md     <-- Entropy, Kolmogorov complexity, entanglement
  quantum.md         <-- Why quantum computing doesn't help
novel/
  info_computation_gap.md  <-- Best original insight (delta(n) framing)
  entropy_measurements.md  <-- Empirical: 5.04 bits/prime, gap analysis
  failure_taxonomy.md      <-- Three failure modes classification
algorithms/              <-- ONLY working, tested code
literature/
  references.md          <-- All citations in one place
  state_of_art_2026.md   <-- Latest published results
experiments/             <-- ALL experiments organized by topic
  analytic/ algebraic/ quantum/ ml/ information_theory/
  dynamical/ topological/ sieve/ other/
archive/                 <-- Old session dumps (read-only reference)
```

## The Barrier in One Paragraph
p(n) = SMOOTH(n) + RANDOM(n). The smooth part R^{-1}(n) is O(polylog) and gives
~50% of digits. The remaining ~50% encode oscillatory contributions of ~10^48
Riemann zeta zeros with GUE-random phases, information-theoretically incompressible.
Best known: O(x^{2/3}) combinatorial, O(x^{1/2+epsilon}) analytic. For p(10^100):
minimum ~10^49 operations. At 10^15 ops/sec = 10^34 seconds = 10^26 years.

## Three Failure Modes (Every Failed Approach Falls Into One)
1. **Circularity** -- needs primes to compute primes (Mills, Chebotarev, EGF)
2. **Equivalence** -- reduces to explicit formula / zeta zero sum (all spectral methods)
3. **Information Loss** -- smooth approximations lose ~170 bits distinguishing p(n)

## Do NOT Re-explore These (Thoroughly Closed)
- Zeta zero summation (always O(sqrt(x)) terms, GUE prevents compression)
- ML/neural (1.1% best, fundamentally limited by Prime Coding Theorem)
- Modular/CRT reconstruction (circular or same complexity as pi(x))
- Quantum algorithms (volume-law entanglement, no speedup beyond O(x^{1/3}))
- Exotic algebra (topos/motivic/prismatic/condensed/perfectoid -- all closed)
- Dynamical systems (FRACTRAN/CA/transfer operators -- all closed)
- p-adic/adelic (not q-adically continuous, Mahler coefficients diverge)
- Interpolation/regression on delta(n) (random walk, zero generalization)

## Viable Research Directions
1. **Circuit complexity of pi(x)** -- TC^0? NC^1? (genuinely unstudied)
2. **Time-bounded Kolmogorov complexity of delta(n)** (open, connects to circuit bounds)
3. **Zeta zero compressibility** -- if zeros have exploitable structure, barrier falls
4. **Concrete Berry-Keating Hamiltonian** -- Hilbert-Polya conjecture (quantum, open)
5. **Novel number-theoretic identity** relating sum_rho R(x^rho) to computable function

## Rules for AI Agents
1. **Read `status/CLOSED_PATHS.md` before proposing ANY approach**
2. Save experiments to `experiments/<topic>/` with descriptive filenames
3. Update `status/CLOSED_PATHS.md` when closing an approach (add a row)
4. Genuinely novel findings go to `novel/` with evidence
5. New literature goes to `literature/`
6. **UPDATE this CLAUDE.md** when you discover something that changes the picture:
   - New barrier proven? Update "The Barrier" and "Do NOT Re-explore" sections
   - New viable direction found? Add to "Viable Research Directions"
   - Closed a direction from OPEN_PROBLEMS.md? Move it to "Do NOT Re-explore"
   - Found a better algorithm? Update "Status" section and `status/BEST_ALGORITHMS.md`
   - Keep this file accurate and current -- future sessions depend on it
7. If you beat the current best algorithm, save it to `algorithms/` with benchmarks
8. DO NOT modify `run.sh`
9. Use sub-agents to save context window
10. When you run out of context, just stop -- the system will restart you
11. If you find the breakthrough, respond with exactly: I FOUND IT!!!
