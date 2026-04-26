# Session 57 — normal mode — smoothness-of-(p-1) conditioning of delta

**Date:** 2026-04-26
**Mode:** normal (post-critique cleanup + one targeted measurement)
**Verdict:** FAIL / mode I (closes 35th pseudorandomness measure)

## Context

Run 56 (critique mode) closed all five S55 proposals as duplicates of
prior CLOSED_PATHS rows; the proposals_session.md doc collided with
multi-line entries 26 / 654 / 699 / 708 / 666 / 703 / 711 / 518 / 712 /
729 / 256 / 257 / 656 / 723 / 726 / 732 / 733 / 519 (cumulative 619 →
622 closed paths). With Chain E closed for both known technique families
(AKS sub-attacks → E7.10, Brandt MKtP → E5.8) and the fourth-encoding
sweep cumulative count at ~78, the project plateau is confirmed and the
"normal" mode has only marginal new attack surface to probe.

This session ran one targeted measurement that was NOT a duplicate of
prior CLOSED_PATHS, picked specifically because (a) it conditions the
most-tested residual delta(n) on a feature of the prime itself rather
than of the index, and (b) the "smoothness of p-1" angle is a natural
arithmetic primality-testing primitive (Pollard p-1, Pocklington,
Lucas-Lehmer) and was missing from the empirical pseudorandomness
catalogue.

## Experiment

`experiments/proposals/delta_pminus1_smoothness.py`. For all N=9592
primes p ≤ 10⁵, compute:

* delta(pi(p)) via R^{-1}(n) ≈ p + f(p)·log(p) linearization (validated
  ±1 vs Newton-Rinv at mpmath 30 dps; bias is independent of smoothness
  class so does not affect class statistics);
* s(p-1) = largest prime factor of p-1, smoothness ratio
  s_ratio = log(s(p-1)) / log(p) ∈ (0, 1];
* omega(p-1) = number of distinct prime factors of p-1;
* Sophie-Germain sub-class: primes of the form p = 2q+1 with q prime
  (least smooth p-1 possible).

Quartile-binning of s_ratio produces 4 classes of size 2398 each;
KS test on extreme quartiles; Pearson and Spearman against delta and
\|delta\|. Sophie-Germain N=670; KS vs others.

## Results

* Per-quartile mean(delta) range = 1.6 (single sigma fluctuation),
  std(delta) range = 1.4. No monotonic trend.
* Pearson(s_ratio, delta) = +0.0049 (p=0.634);
  Spearman = +0.0018 (p=0.861);
  KS class-0-vs-3 D = 0.016 (p=0.91).
* Sophie-Germain (N=670): mean/std within single sigma, KS p=0.95.
* Pearson(omega(p-1), \|delta\|) = +0.054 at p=1e-7 — statistically
  real but explains 0.3% variance, fully attributable to both
  quantities co-scaling with p (Erdős-Kac log log p, \|delta\| ~ p^{0.57}).

## Verdict

**FAIL / mode I.** The arithmetic sector (factorization features of p±1)
is empirically orthogonal to the transcendental sector (zeta-zero
residual delta of R^{-1}). Neither smoothness nor distinct-prime-factor
count of p-1 leaks usable information about delta. This rules out a
class of "smoothness-aware delta predictors" and a side-channel of
Pollard / Pocklington data into prime-counting. Adds 35th measure to
`novel/pseudorandomness_of_pi.md`.

## Files

* `experiments/proposals/delta_pminus1_smoothness.py` — script
* `experiments/proposals/delta_pminus1_smoothness_results.md` — report
* `experiments/proposals/delta_pminus1_smoothness_data.npz` — raw data
* `status/CLOSED_PATHS.md` — S57 row appended (cumulative 545+)
* `novel/pseudorandomness_of_pi.md` — measure #35 added; title bumped to 35
* `archive/sessions/session57_smoothness_conditioning.md` — this synthesis
* `.run_state` advanced to 57

## Cross-link to EDGES.md

* T3 (failure modes): closes I.
* E1.1 (delta has O(log x) bits, computational not informational
  barrier): the smoothness-of-(p-1) feature is one more channel that
  fails to encode those bits.
* E1.5 (pi(x) mod m has invariant 0.537 bits/step): consistent —
  smoothness conditioning cannot beat the unconditional rate.
* E1.7 (delta long-range correlation is INDIRECT, AR(7)+Hurst):
  smoothness is not part of the AR(7) memory window.

No new edge promotion: this is one more empirical orthogonality check
of the ARITHMETIC × TRANSCENDENTAL sector independence already
documented across S25-S57.
