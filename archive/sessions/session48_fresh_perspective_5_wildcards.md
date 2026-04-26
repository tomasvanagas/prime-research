# Session 48 -- Fresh-Perspective Wildcards (5 angles)

**Date:** 2026-04-26
**Mode:** fresh-perspective (no CLOSED_PATHS pre-read)
**Net status:** all 5 angles closed; literature scan returned nothing new.

## What was attempted

Five fresh-angle wildcard experiments selected after a brief audit
showed they had not been explicitly tried among the existing 142
wildcard experiments:

1. **`newton_diffs_decay.py`** -- Newton forward-difference series
   for pi(x). Test whether truncated Newton series gives polylog eval.
2. **`multilinear_ext_pi.py`** -- Multilinear extension of chi_P
   over the boolean hypercube, in view of GKR/sumcheck protocols.
3. **`parallel_sieve_matrix_exp.py`** -- Sieve as a linear dynamical
   system; CP/Tucker rank decomposition of alive-vector via wheel.
4. **`dirichlet_char_decomp.py`** -- chi_P spectrum in the Dirichlet
   character basis modulo prime q.
5. **`orbit_resum_pi.py`** -- Semiclassical block-resummation of zeta
   zero contributions in the explicit formula.

## Verdicts

| # | Angle | Verdict | Mode | One-line key finding |
|---|-------|---------|------|----------------------|
| 1 | Newton series | CLOSED | I | log2(Delta^k pi(0)) ~ 0.995k - 3.17 (pure 2^k growth); divergence + Omega(x) bit complexity |
| 2 | Multilinear extension | CLOSED | I | nonzero-coef ratio converges to ~0.71; max coef grows as 2^{n/2}; matches existing rank(pi_N)=2^{N/2-1}+2 (E2.1) |
| 3 | Linear sieve dynamics | CLOSED | E | Rank-1 wheel factorization stops at primorial(K) > N; reshape rank Theta(sqrt N) (matches MPS bond-dim) |
| 4 | Dirichlet char spectrum | CLOSED | E | L1/L2 ~ 0.85 sqrt(phi(q)) -- white-noise limit; equivalent to Polya-Vinogradov barrier |
| 5 | Block-resummation of zeros | CLOSED | E | Per-block contribution grows like 2^{i/2}; log-block grouping gives no advantage over flat truncation |

(Modes: C = circular, E = equivalent to known barrier, I = info loss.)

All five collapse to existing project barriers. None opens new attack
surface. None produces a partial / approximate algorithm worth
keeping.

## Quantitative highlights

The ONE genuinely useful number obtained from this session is the
linear fit on Newton differences:

    log2 |Delta^k pi(0)| ≈ 0.995 * k - 3.17    (for k = 1..200)

This is a **clean empirical confirmation that pi(x) is, in the
forward-difference sense, an order-1 entire function** -- the
sharpest possible "non-smooth" behavior. Combined with the multilinear
test (max coef ~ 2^{n/2}) we now have two independent measurements
of the same Omega(2^{N/2}) coefficient-magnitude wall.

These confirm but do not add fundamentally new structure to E2.x
(algebraic edges) and the existing determinantal complexity result.

## Literature scan (parallel sub-agent)

Targeted topics: pi(x) algorithms, GKR/sumcheck for #P counting,
Barvinok short generating functions, MLE/sumcheck for prime indicators,
explicit-formula evaluation, Gutzwiller / periodic-orbit truncation,
Newton-series for arithmetic functions, tensor networks for arithmetic
indicators, holographic algorithms, quantum (non-Shor) for pi(x).

**Nothing new.** Three already-cataloged 2025-2026 contributions
(Aggarwal 2510.16285, Cully-Hugill-Lee 2402.04272, Inoue 2604.05733)
are the only relevant items; all three are in
`literature/state_of_art_2026.md` already. Topics 2, 3, 4, 7, 9, 10
returned **zero relevant publications**.

The cached `literature/state_of_art_2026.md` is current; no update
required.

## Status-file impact

CLOSED_PATHS.md additions (5 entries to be added):
- Newton forward-difference series for pi(x) -- I-mode, divergent
- Multilinear extension of chi_P with sumcheck -- I-mode, dense MLE
- Linear sieve as fast-forwardable matrix -- E-mode, wheel barrier
- Dirichlet character spectrum of chi_P -- E-mode, Polya-Vinogradov
- Block-resummation of zeta zeros -- E-mode, sqrt(x) bound

OPEN_PROBLEMS.md: no change. Chain E (AKS-family TC^0) still open.

EDGES.md: no new edges. Findings are quantitative refinements of
E2.1 (MPS bond-dim), E3.x (analytic), E5.x (sieve barrier).

## What's still worth trying

Per CLAUDE.md guidance ("construction is encouraged"):
- The three FOCUS-1 sub-attacks (Bernstein 2003 smaller-r AKS,
  non-cyclotomic ring AKS, Healy-Viola Frobenius transplant) remain
  the highest-leverage CONSTRUCTION tasks and were not touched here.
- FOCUS-3 (Brandt-MKtP) also still open.

## Files

Code + results (5 each, all in `experiments/wildcard/`):
- `newton_diffs_decay.py` + `_results.md`
- `multilinear_ext_pi.py` + `_results.md`
- `parallel_sieve_matrix_exp.py` + `_results.md`
- `dirichlet_char_decomp.py` + `_results.md`
- `orbit_resum_pi.py` + `_results.md`
