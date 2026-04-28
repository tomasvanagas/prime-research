# Q-reduced form signature extraction battery

This script extracts richer signals from the q-reduced form of a divisor
(beyond the rank): D'(q) margin, support size, max chip, total non-q
chips, and the "winnable subtraction set" `W(D) := {v : D − δ_v is
winnable}`.

**Status:** built and used in S161. See main writeup
`baker_norine_chi_p_results.md` for full results.

**Key finding:** D'(q) for the prime divisor on Hasse(N) equals π(N)
exactly while matched-density random divisors give E ≈ π(N)²/N — z-score
grows as `√π(N) ~ √(N/log N)`. The W-set restricted to primes has all
primes for D_P (trivially, since D_P(p) = 1) but reduces to a tautological
"primes are in supp(D_P)" observation under analysis.

**Output:** `sig_N16_32.json` and consolidated into `full_results.json`.
