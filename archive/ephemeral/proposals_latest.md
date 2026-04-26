I'll dive into this proposal session. Let me start by understanding the workspace and then develop concrete proposals with runnable code.
Good, the environment is set up. Now I'll write up the proposals document and prepare runnable experiments. Let me develop 4 concrete proposals.
Now let me write the four experiments. I'll start with the most important ones.
Proposal A is partially encouraging — works for x=100, 1000 but fails at x=10000. Let me write proposal B.
The R(x^rho) computation in Proposal A blew up because the Gram series doesn't converge for complex arguments with large imaginary part. Let me fix it using a proper li/ei evaluation for complex arguments.
Proposal A path closes cleanly: zero-sum partial sums oscillate without smooth-envelope structure, so Wynn extrapolation fails catastrophically. Let me finish with cleanup.
No `__pycache__` to clean, no missing results files. Now write the run state.
Session complete. Four concrete proposals filed in `archive/ephemeral/proposals_session.md`, each with a runnable experiment and `_results.md` companion in `experiments/proposals/`:

| # | Proposal | Verdict |
|---|----------|---------|
| A | Pade/Wynn extrapolation of psi(x) zero-sum | **Closed** — partial sums oscillate, no smooth envelope; Wynn fails |
| B | Borel-resummation of Cipolla asymptotic | **Closed** — Borel-Pade does not improve over best truncation |
| C | Hecke-eigenvalue fingerprint primality oracle | **Partially open** but blocked: tau(n) for composite n requires factoring |
| D | Zero-aware control variate for MC pi(x) | **Closed** — empirical decay slope ~-0.25, far worse than RH-bound -1; PNT control variate gives only constant-factor variance reduction |

`.run_state` set to 45.
