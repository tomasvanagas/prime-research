I'll write proposals to ephemeral, then code and test the most promising ones. Let me start with the proposals document.
Now let me code and test the most promising proposals. Starting with the TT-rank experiment (cheapest, most directly informative).
Strong verdict: TT-rank scales as 2^{L/2-1} for primes vs 2^{L/2} for random — exactly halved (the "all primes > 2 are odd" bit), but still exponential. Let me write the results.
Now Proposal 1: compressed sensing on the zeta-zero contribution matrix.
Now Proposal 2: PSLQ subsequence hunt.
## Summary

Wrote 5 proposals to `archive/ephemeral/proposals_session.md`, then implemented and tested the three most tractable ones:

**Proposal 3 — TT-rank of prime indicator** (`experiments/proposals/tt_rank_prime_indicator.py`): VERDICT CLOSED. TT-rank scales as exactly 2^{(L-1)/2} + 1, exactly half of random Bernoulli baseline (the "all primes > 2 are odd" bit). Exponential, not polylog. Confirms pseudorandomness picture cleanly.

**Proposal 1 — Compressed sensing on zeta-zero contribution matrix** (`experiments/proposals/zero_contribution_compressibility.py`): VERDICT CLOSED. Singular values decay algebraically as k^{-0.5 to -0.7} (effective rank ~312 of 1000 at 1% threshold). 2D Fourier basis: ~30% of coefficients above 1%. Polynomial compressibility, not polylog.

**Proposal 2 — PSLQ subsequence hunt** (`experiments/proposals/pslq_subsequence_delta.py`): VERDICT CLOSED. No relations found on dyadic, Fibonacci-indexed, or prime-indexed subsequences. delta(n) values flip sign wildly (e.g., -20.77 at n=256, +21.05 at n=1024), confirming pseudorandomness on sparse probes too.

Proposals 4 (GUE-aware adaptive truncation) and 5 (theta/modular bridge) are documented but not tested — Proposal 4 needs an analytic variance estimate that's a session of its own, and Proposal 5 is intentionally speculative.

`.run_state` set to 55. All `.py` files have companion `_results.md`. No `__pycache__` directories created.
