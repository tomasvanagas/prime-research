# Critique Session — S51 (critic of S50-fresh proposals)

Date: 2026-04-26
Mode: critique
Source proposals: `archive/ephemeral/proposals_session.md` (S50-fresh, four proposals A/B/C/D)
Source experiments: `experiments/proposals/proposal_b_cf_prime_constant_results.md`,
                    `experiments/proposals/proposal_c_zero_algebraicity_results.md`

## Outcome

All four proposals classified as duplicates of pre-existing CLOSED_PATHS
entries:

| proposal | object | verdict | matches |
|---|---|---|---|
| A | tropical Hankel of π | DUPLICATE | 199, 318, 423, 424, 652 |
| B | regular CF of α=Σ 2^{-p} | DUPLICATE-PLUS (object new) | 288, 298, 588, 615 |
| C | PSLQ on γ_k − γ_1 | DUPLICATE | 23 / E7.1 |
| D | reservoir fit on π(n) | DUPLICATE | 683 |

## Actions taken

1. **Critique written** to `archive/ephemeral/critique_latest.md` with
   per-proposal failure-mode classification and CLOSED_PATHS line-number
   references.

2. **CLOSED_PATHS.md** got one new entry for B (regular CF of the prime
   constant α): the *object* was new even though the *mechanism* was a
   strict consequence of E7.7 (representation-changing bijections preserve
   entropy). Khintchine-typical signature in 5/5 measures (K_0, Lévy,
   Gauss-Kuzmin, autocorrelation, k-automaticity, spectral slope).

3. **`novel/pseudorandomness_of_pi.md`** updated to 33 measures: added
   measure 32 (Dirichlet-character spectrum of χ_P, which the index had
   referenced but the table had not listed) and measure 33 (regular-CF
   quotients of α). Header bumped from "31 Measures" to "33 Measures".

4. **No new entries** to OPEN_PROBLEMS or `novel/` standalone — every
   proposal closed under existing edges.

## Cross-cutting observation

The proposer in S50-fresh deliberately avoided reading CLAUDE.md /
CLOSED_PATHS before generating proposals (per the file's own
"discipline: do not look at CLAUDE.md / CLOSED_PATHS.md before
generating" line). This produced four proposals all of which fall
directly under existing closures. The two proposals that ran (B, C)
each consumed compute (PSLQ at 50 dp on 199 zero differences;
mpmath-precision 16800-bit CF expansion to depth 1500) to confirm
predictions already on file.

The "fresh, no prior context" methodology is no longer net-positive at
the project's current maturity (720+ closed paths). The duplicate-
detection step should run **first** in future proposal sessions.

## Underweighted directions (per CLAUDE.md S60)

The proposer touched none of:
- Construction work (FOCUS-1 sub-attacks, FOCUS-3 Brandt MKtP)
- Lower-bound techniques (pebbling, MKtP-style diagonalisation)
- Theoretical sharpening (Lean formalisation of MPS theorem,
  rigorous N/2 threshold proof)

These remain the highest-leverage open work per project frontier.

## Files touched

- `archive/ephemeral/critique_latest.md` (rewritten)
- `status/CLOSED_PATHS.md` (+1 entry, line 734)
- `novel/pseudorandomness_of_pi.md` (header bumped 31→33; +2 table rows)
- `archive/sessions/session_critique50.md` (this file)

## Files NOT touched (intentionally)

- `EDGES.md` — no new edges discovered
- `status/OPEN_PROBLEMS.md` — no resolutions, no new opens
- `algorithms/` — no algorithmic change
- `run.sh`, `FOCUS_QUEUE.md`, `TODO.md` — per protocol
- `CLAUDE.md` — no significant status change
