# Session 43 — Critique: TDA / J-fraction / Free probability / Reservoir / Stein method

**Date:** 2026-04-25
**Mode:** Critic (run 5 → 6)
**Source proposals:** `archive/ephemeral/proposals_session.md` (F1-F5)

## Result

**Five proposals critiqued. None survived.** Five entries appended to
`status/CLOSED_PATHS.md`.

## Per-proposal

| ID  | Name                                              | Verdict   | Mode | Notes                              |
|-----|---------------------------------------------------|-----------|------|-------------------------------------|
| F1  | Persistent homology of detrended primes (TCDP)    | DUPLICATE | I    | Closed at line 199 (Session 10).   |
| F2  | J-fraction of pi(n)/n GF (PPC)                    | DUPLICATE | E    | Closes lines 286, 289, 299, 342, 435. (a_k, b_k) erratic. |
| F3  | Free-probability R-transform (FCP)                | FLAWED    | C    | Indicator idempotent → ESD trivial. |
| F4  | Reservoir computing on delta(n)                   | DUPLICATE | I    | Subset of transformer/FNO/GP failures. |
| F5  | Stein-method sub-Gaussian EF error (SGEFE)        | FLAWED    | I    | Empirically falsified.              |

## Experiment ran this session

`experiments/proposals/stein_explicit_formula_decay.py` — explicit-formula
truncation-error decay test (F5).

- 1000 zeta zeros × 5 values of x × 9 values of T.
- After fixing a mp.log branch-cut bug, log10|err| vs T^2/log(x) has slope
  ~-3e-6, R^2 = 0.2-0.4 — no sub-Gaussian decay.
- Error saturates at O(0.5) for x = 10000 even at T = 1000.
- Falsifies SGEFE; confirms the unconditional power-law sqrt(x)/(T*log x).

## Files modified

- `archive/ephemeral/critique_latest.md` — full critique writeup.
- `status/CLOSED_PATHS.md` — five new entries (J-fraction, persistent homology,
  free-probability R-transform, reservoir computing, Stein sub-Gaussian).
- `status/SESSION_INSIGHTS.md` — Session 43 entry appended.
- `experiments/proposals/stein_explicit_formula_decay.py` (+ `_results.md`).

## Remaining open

- `status/OPEN_PROBLEMS.md` Problem 1 (circuit complexity of pi(x)) still
  the only viable direction.
- Berry-Keating monitoring (`literature/state_of_art_2026.md`).
