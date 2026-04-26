# Session 45 — Wildcard: 5 Fresh-perspective experiments

Date: 2026-04-25
Mode: wildcard (no CLAUDE.md or CLOSED_PATHS read; first-principles brainstorm)

## What was tested
Five new wildcard experiments, each picking an analogy from another field and
checking if it transplants to the prime-counting / explicit-formula problem:

1. **Borel resummation of the explicit formula** (`borel_explicit_formula.py`)
2. **Padé approximants of δ(n) = p(n) − R⁻¹(n)** (`pade_residual.py`)
3. **Fast-multipole expansion of the zero sum** (`multipole_zero_sum.py`)
4. **Continued-fraction expansion of π(x)/li(x)** (`cf_pi_over_li.py`)
5. **Wirtinger-flow recovery of γ_k from π(x)** (`wirtinger_zero_recovery.py`)

All have companion `_results.md` files in `experiments/wildcard/`.
Brainstorm doc: `session45_brainstorm.md`.

## Results — all five close

| # | Approach              | Best evidence                                | Failure mode |
|---|-----------------------|----------------------------------------------|--------------|
| 1 | Borel-Padé            | err 0.22 at x=500 (vs raw 2.04) but doesn't break √x scaling — same γ_max needed | E |
| 2 | Padé of δ(n)          | RMSE 3.5 vs baseline-mean RMSE 4.7; marginal | I |
| 3 | FMM multipole         | Diverges (P=24 error 10³⁷); kernel oscillation defeats local Taylor | I |
| 4 | CF of π(x)/li(x)      | Geo mean 2.52 ≈ Khinchin 2.685; max a_i unbounded | I |
| 5 | Wirtinger recovery    | Basin radius ≲ 0.1 around truth; 0/5 random init at K=10 | E |

## Mildly novel diagnostics produced

- **Multipole impossibility quantified**: cluster Taylor needs P ~ e·R·log(x);
  any clustering with B clusters needs P · B ≥ K · log(x), so no √K speedup
  achievable — *first explicit constant in the prior pessimistic claim*.
- **Wirtinger basin-radius vs K table**: documents how fast the recovery problem
  becomes intractable as we ask for more zeros.
- **CF Khinchin test as 22nd pseudorandomness measure**: formerly untested in
  wildcard list. π(x)/li(x) and the normalized residual both look like generic
  irrationals.

None of these is independently a "novel result" worth promoting to `novel/` —
they reinforce the existing pseudorandomness/oscillatory-incompressibility
picture by closing three additional named techniques.

## Open follow-up (not pursued in this session)

- The Borel-Padé regularization at moderate x is real (10× error reduction at
  x=500). A targeted study of the optimal Borel order vs x might give a useful
  *constant-factor* improvement for analytic prime counting libraries.
- Riemann–Siegel transplant (mentioned in brainstorm but not implemented) is
  the next natural angle: their split into main + correction terms is the only
  technique that has ever broken a √-barrier in this neighborhood.

## State
`.run_state` set to 8 (advance to next mode in run.sh rotation).
