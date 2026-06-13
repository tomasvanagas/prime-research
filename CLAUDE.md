# Prime Research: Computing p(n) Exactly Without Bruteforcing

## Goal

Find an O(polylog) algorithm to compute the nth prime p(n)
(equivalently π(x)) exactly. Target: p(10^100) in <1 second, 100%
accurate.

## How this project runs (v2, since 2026-06-13)

A continuous loop (`./run.sh`) starts fresh research cycles (Fable 5,
xhigh effort). There are **no modes, no session grades, no rotation,
no commit threads** — the v1 multi-stage framework is retired and
archived at `archive/framework_v1/`. Each cycle:

1. reads `PROGRAM.md` — the single live state document (what is done,
   what is open, ONE concrete next action);
2. does real work toward the goal;
3. updates `PROGRAM.md` before stopping.

## Honest status

Every known *computation* route is blocked: 730+ approaches catalogued
in `status/CLOSED_PATHS.md`, four obstruction theorems, and the
information-theoretic barrier — p(n) = SMOOTH + RANDOM, where the
smooth part R⁻¹(n) is polylog and gives ~50% of digits, and the rest
encodes ~√x zeta-zero phases with GUE-random behaviour. Best known:
O(x^{2/3}) combinatorial, O(x^{1/2+ε}) analytic. The goal may require
mathematics that does not exist yet.

The active line (see `PROGRAM.md`) exploits what is NOT blocked:
verification, certification, structural measurement, adjacent-problem
wins. The deliverable of a cycle is honest progress — a working
algorithm or protocol, a new structural fact with measurements, a
precisely-closed question with its mechanism — never a faked
breakthrough.

## The working contract

- **Pick ONE item per cycle**: the NEXT ACTION in `PROGRAM.md`, or a
  better idea explicitly justified against `PROGRAM.md` and the
  catalogues. Multi-cycle builds are normal — record design and
  partial state in `PROGRAM.md` so the next cycle continues mid-build.
- **Build code that runs**: `experiments/<topic>/<name>.py` — one
  script per experiment, CLI-parameterised (never `_v2.py` /
  `_quick.py` variants) — plus `<name>_results.md`. Include
  `--selftest`; every boundary case you debug becomes a selftest case.
- **Measure**, with matched baselines/controls where meaningful.
- **Every results file states what would falsify the result.**
- **Search before claiming novelty**: `status/CLOSED_PATHS.md` and
  `EDGES.md` first. The bar for `novel/` is: a paper-grade number
  theorist or complexity theorist, after one careful read of the
  prior literature and the catalogue, could not produce it.
- **Honest reporting is absolute.** Failures are filed as failures
  with the structural mechanism. Earlier claims (including your own
  cycles') get corrected when wrong — that is progress. A clean
  negative with a mechanism beats a vague positive. Confirming a
  known wall with a new measurement is a measurement, not novelty.
  Inflated claims poison future cycles.

## File placement

- Live state → `PROGRAM.md` (keep short; one entry per cycle).
- Experiments → `experiments/<topic>/<name>.py` + `<name>_results.md`.
- Genuinely original findings → `novel/<descriptive_name>.md`.
- Proven results / barriers → `proven/`.
- Working tested algorithms → `algorithms/`.
- Cycle transcripts → `archive/CLAUDE_OUTPUTS/` (written by run.sh).
- v1-era session syntheses → `archive/sessions/` (read-only history).

## Reference catalogues (read-only unless you add real content)

- `status/CLOSED_PATHS.md` — 730+ tested-and-closed approaches. Add a
  row when you close something; cite the mechanism.
- `EDGES.md` — verified mathematical edges, cited by ID (e.g. E2.1).
- `OPEN_POSITIVE_TARGETS.md` — adjacent-problem targets (P6 carries a
  duplication warning — read it before spending anything there).
- `ATTACK_VECTORS.md`, `NOVELTY_CHALLENGES.md`, `RESEARCH_AGENDA.md`,
  `CROSS_DOMAIN_TECHNIQUES.md`, `status/SESSION_INSIGHTS.md` —
  v1-era catalogues; useful as idea sources and history.
- `novel/succinct_verification_of_pi.md` — the active line's
  canonical document (protocol stack, theory, open program).

## Cleanup before stopping

- `PROGRAM.md` updated (produced / open / single NEXT ACTION).
- Every new `.py` has its `_results.md`:
  `find experiments/ -path "*/.lake" -prune -o -name "*.py" -print |
   while read f; do r="${f%.py}_results.md"; [ ! -f "$r" ] && echo
   "MISSING: $r"; done`
- No `__pycache__` directories left behind.

## Breakthrough protocol

Write `BREAKTHROUGH.md` ONLY with machine-verified, reproducible
evidence of the actual goal: exact π(x)/p(n) at large x in polylog
time — verified output equality against ground truth, measured
scaling, and reproduction instructions. The loop halts when the file
exists. Never write it on heuristics, approximations, or conditional
results. A false breakthrough is the worst possible output of this
project.
