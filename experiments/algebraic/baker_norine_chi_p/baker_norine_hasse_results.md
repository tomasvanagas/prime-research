# Hasse cover graph Baker-Norine experiments

This script extends `baker_norine_chi_p.py` to the Hasse cover graph
`H_N := ([1, N], {(a, pa) : p prime, pa ≤ N})` and to ratio-graph
variants `(a, ka)` for k ∈ {2, 3, ...}.

**Status:** built and used in S161. See main writeup
`baker_norine_chi_p_results.md` for full results, theorems, falsifier
verdicts, and self-graded synthesis.

**Key finding:** on `H_N` with sink q=1, the q-reduced form of the
prime divisor is `D'_P^N = π(N) · δ_1` (Theorem E2.28-2). Verified for
N ∈ {16, 32, 64, 128, 256, 512}.

**Output:** `hasse_results.json` (raw z-scores), consolidated into
`full_results.json` by `run_full.py`.
