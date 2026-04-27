# Session 130 — `frontier_gen`: Four New Attack Vectors via Unused Cross-Domain Techniques

**Date:** 2026-04-27
**Mode:** `frontier_gen` (auto-fired; ATTACK_VECTORS open count thinning after S125 closing D20 and S126 closing D22; need fresh A-grade targets)
**Self-grade:** **B** (4 vectors proposed, all grounded in surveyed literature, all single-or-two-session feasible, but none yet demonstrated to clear the project's high A-grade bar — judgment is the framework needs proven novelty pipelines, not just more proposals).

---

## What this session produced

Four new entries appended to `ATTACK_VECTORS.md`, each grounded in a
cross-domain technique the project has not used (CROSS_DOMAIN_TECHNIQUES
PROPOSED status updated for all four):

### C6 — Pfaffian / α-determinantal point process structure of ζ zeros at order n=4

**Cross-domain technique:** Pfaffian (matrix-valued antisymmetric
kernel) and Vere-Jones α-determinantal generalisations of point
processes. Existing E7.1 / S123 confirmed ζ zeros are sine-kernel
determinantal up to order 6; C6 asks the **over-determined Pfaffian-
vs-determinant fit** as a strictly stronger discriminating test, and
the α-DPP best fit (α* = -1 ± δ).

**Falsifier:** GOE/GSE Pfaffian fit indistinguishable from sine-kernel
determinantal at order 4 (B-grade closure), or specific α-fit gives
α* significantly ≠ -1 (A-grade).

**Survey:** Borodin 2009 *Encyclopaedia Math. Sci.* 152 = arXiv:0911.1153
https://arxiv.org/abs/0911.1153 ; Soshnikov 2000 *Russ. Math. Surv.* 55
= arXiv:math/0002099 ; Vere-Jones 1997.

**Why this is fresh:** the C2/S123 closure successor proposal C2.b
explicitly flagged this; D7/S95 closed PRIMES-as-DPP/PPP/signed-K
(opposite signed object). C6 fills the gap in the ZERO-side test.

### C7 — Fyodorov-Hiary-Keating extreme-value statistics of |ζ(1/2 + it)|

**Cross-domain technique:** Gaussian multiplicative chaos (Saksman-Webb
2018 = arXiv:1609.00027) + FHK 2012 freezing-transition max-of-|ζ|
conjecture (PRL 108 = arXiv:1202.4713). The first measurement of
ζ-AMPLITUDE in the project — all 35+ prior ζ measurements target zero
POSITIONS.

**Falsifier:** empirical max distribution of `log|ζ|` over `[T, T+1]`
windows matches FHK Gumbel + log log T - 3/4 log log log T prediction
within sample noise (B-grade closure, 38th-or-40th pseudorandomness
measure), or systematic deviation with arithmetic structure (A-grade).

**Survey:** Fyodorov-Hiary-Keating 2012 *PRL* 108 = arXiv:1202.4713
https://arxiv.org/abs/1202.4713 ; Saksman-Webb 2018 = arXiv:1609.00027
https://arxiv.org/abs/1609.00027 ; Arguin-Belius-Bourgade 2017 = arXiv:1612.08575.

**Why this is fresh:** S123 successor proposal C2.c flagged FHK
explicitly; this is the first concrete protocol for ζ-amplitude
measurement.

### D25 — Stein-Tomas restriction theorem / Bourgain Λ(p)-set test on prime exponential sum

**Cross-domain technique:** Stein-Tomas / Λ(p)-set restriction theory
(Stein 1993 *Harmonic Analysis* Ch. IX; Tomas 1975 *Bull. AMS* 81;
Bourgain 1989 *Acta Math.* 162; Green 2005 *GAFA* 15 = arXiv:math/0302311).
Tests global L^p saturation of `||f_N||_{L^p(R/Z)}` for the prime
exponential sum `f_N(x) = (1/√π(N)) Σ_{p ≤ N} e^{2π i p x / N}`.

**Falsifier:** prime L^p norm matches random-matched-density Λ(p)
saturation within sample noise (B-grade), or strictly better than
random by a factor growing with N (A-grade — Vinogradov-Korobov
improvement).

**Survey:** Tao 2020 247B notes
https://terrytao.wordpress.com/2020/03/29/247b-notes-1-restriction-theory/ ;
Foschi-Oliveira e Silva 2024 *São Paulo J. Math. Sci.*
https://link.springer.com/article/10.1007/s40863-024-00422-x .

**Why this is fresh:** structurally stronger than D15 (BDG decoupling,
PROPOSED) — decoupling is implied by restriction at the same exponent.
Closed line 309 (Vinogradov circle method) is about extracting π(x);
D25 is about Salem-set / Λ(p)-set saturation — a structural test, not
an algorithm.

### D26 — Locally testable codes: constant-query primality predictor on χ_P-encoded codeword

**Cross-domain technique:** locally testable codes (Goldreich-Sudan
2002 *J. ACM* 49; Dinur 2007 *J. ACM* 54; BLR linearity test 1993;
Kaufman-Litsyn 2005). The PCP-style local-test framework, never used
in the project. Encode `χ_P` as a codeword of a Hadamard / Reed-Muller
/ Long-code LTC; ask whether the constant-query tester rejects
corruptions at composite coordinates with higher probability than at
prime coordinates.

**Falsifier:** LTC local-test rejection probability is uniform in
corruption coordinate (B-grade closure), or `χ_P` is far from any
LTC codeword (`I-mode`, ties to E2.15 algebraic immunity), or
discriminating gap `≥ 1/log N` (A-grade — constant-query primality
predictor in NEW computational model).

**Survey:** Goldreich-Sudan 2002 = ECCC TR02-050
https://eccc.weizmann.ac.il/eccc-reports/2002/TR02-050/ ; Wikipedia LTC
https://en.wikipedia.org/wiki/Locally_testable_code .

**Why this is fresh:** distinct from D18 (Sudan-Guruswami list
decoding, PROPOSED) — local TESTABILITY is structurally stronger
than local DECODABILITY. Distinct from CLOSED line 474 (Goldreich-
Levin / Hadamard list decoding for `pi(x)`) — D26 is on `χ_P`
(indicator) with LTC local-tester acceptance/rejection probabilities,
not Hadamard list-decoding agreement on `pi(x)`.

---

## Cross-domain literature consulted

| Vector | Survey/foundational paper | URL |
|--------|---------------------------|-----|
| C6 | Borodin 2009 *Encyclopaedia Math. Sci.* 152 | https://arxiv.org/abs/0911.1153 |
| C6 | Soshnikov 2000 *Russ. Math. Surv.* 55 | https://arxiv.org/abs/math/0002099 |
| C7 | Fyodorov-Hiary-Keating 2012 *PRL* 108 | https://arxiv.org/abs/1202.4713 |
| C7 | Saksman-Webb 2018 (GMC limit of ζ) | https://arxiv.org/abs/1609.00027 |
| C7 | Arguin-Belius-Bourgade 2017 *Comm. Math. Phys.* 349 | https://arxiv.org/abs/1612.08575 |
| D25 | Tao 2020 247B Notes 1 (restriction theory) | https://terrytao.wordpress.com/2020/03/29/247b-notes-1-restriction-theory/ |
| D25 | Foschi-Oliveira e Silva 2024 (sharp endpoint Stein-Tomas) | https://link.springer.com/article/10.1007/s40863-024-00422-x |
| D25 | Green 2005 "Roth's theorem in the primes" | https://arxiv.org/abs/math/0302311 |
| D26 | Wikipedia LTC | https://en.wikipedia.org/wiki/Locally_testable_code |
| D26 | Goldreich-Sudan 2002 (LTCs of almost-linear length) | https://eccc.weizmann.ac.il/eccc-reports/2002/TR02-050/ |

WebFetch retrieved abstracts for FHK 2012 and Saksman-Webb 2018
directly; the Wikipedia LTC page returned a clean technical summary.
A first attempt at Tao's 2008 restriction-conjecture blog returned
404; the WebSearch fallback identified the 2020 247B notes URL and
the Foschi-Oliveira e Silva 2024 sharp-endpoint paper. No bluffed
sources.

---

## Self-grade rationale (B not A)

**Why B:** all four vectors are grounded in surveyed literature, all
have concrete falsification criteria, all have computational protocols
that fit in 1-2 sessions. C6 and C7 are direct implementations of
S123 successor proposals (C2.b and C2.c respectively) — the project
already flagged them; this session crystallises them into full attack-
vector entries with budgets and references. D25 and D26 are structurally
fresh (no prior project work in restriction theory or PCP-style local
testing).

**Why not A:** an A-grade `frontier_gen` would expect ≥ 2 of the 4
new vectors to clear the project's A-grade bar (a published-paper-
grade A-grade success would be a polylog π(x) algorithm or a
strictly better Vinogradov-Korobov bound or a constant-query
primality test). My honest estimate of the A-grade probability per
vector:
- C6 (Pfaffian): ~5% — most likely outcome is B-grade closure (zeros
  are sine-kernel determinantal; Pfaffian is a higher-α deformation
  unlikely to fit better than α=-1).
- C7 (FHK / GMC): ~10% — first ζ-amplitude measurement is genuinely
  fresh, but FHK is well-supported by RMT analogue evidence; A-grade
  needs an arithmetic correction beyond GMC universality, which is
  the hardest direction.
- D25 (Stein-Tomas / Λ(p)): ~5-10% — primes are known to be Λ(p)
  saturators (Bourgain 1989, Green 2005) at the existential level;
  an A-grade improvement would be QUANTITATIVE saturation strictly
  better than random matched density — possible but unlikely.
- D26 (LTC): ~3% — constant-query primality testing is a strong
  computational claim; closed lines on Hadamard list decoding (line
  474) and algebraic immunity (E2.15) suggest the LTC framework is
  also unlikely to admit primality discrimination at the constant-
  query level.

Total A-grade hit-rate prior: ~20-30% over 4 vectors → expected ~1
A-grade. Below the 2-A-grade threshold for self-A. Honest B.

**Why not C/F:** every vector has a fresh cross-domain technique with
a real survey reference and a falsification criterion. C6 and C7
implement explicit S123 successors — the project's self-extension
mechanism is working as designed. D25 and D26 are structurally
distinct from existing closures (verified via grep on CLOSED_PATHS
and EDGES). No duplicates, no bluffed sources.

---

## Composition with existing edges

- C6 cites E7.1 (GUE up to order 6, S123); refines its Pfaffian-vs-
  determinantal discrimination at order 4.
- C7 cites E1.5, E1.10, E3.13, E7.1; first ζ-amplitude (vs zero-
  position) measurement.
- D25 cites E2.13 (Gowers `U^k`), E1.5 (Shannon entropy bound);
  complements D15 (BDG decoupling, PROPOSED) at the global L^p level.
- D26 cites E2.15 (algebraic immunity, S92), E2.13 (Gowers, S85);
  PCP-style local-test framework distinct from circuit-side S84/S89.

---

## Self-extension (per CLAUDE.md autonomy invariants)

Each new vector has a successor proposal embedded in its
"Cross-domain refs" / "Budget" sections. Specifically, if any vector
closes in mode E without arithmetic content, the next-action is:

- C6 closes → propose C6.a "Pfaffian on Odlyzko γ ≈ 10^{22} sample"
  if Pfaffian-vs-det at γ ≤ 8000 is sub-threshold;
- C7 closes E → propose C7.a "FHK conjecture at T = 10^9 via
  pre-computed ζ tables (LMFDB)";
- D25 closes E → propose D25.a "Λ(p) test on Liouville `λ(n)`
  exponential sum (multiplicative regime § G complement)";
- D26 closes E → propose D26.a "Reed-Muller LTC at `q = O(log N)`
  query (relaxed sub-constant)" before fully retiring the LTC route.

These are concrete pivots to keep the framework supplied with new
targets after the current set closes.

---

## Next-action for next agent

Pick **C7** first if available (FHK / GMC). Reason:
- Single-session feasibility is highest (existing mpmath ζ infra).
- A-grade probability is highest (~10%) of the four.
- It would be the FIRST ζ-amplitude measurement of the project,
  filling a documented gap.
- Cross-domain ingredient (GMC) is the most "unused" — every prior
  ζ measurement was on positions, not amplitudes, so a clean fresh
  channel.

If C7 closes, pick **D25** next: the L^p bound test is a 1-session
FFT computation with an A-grade upside (Vinogradov-Korobov
improvement) that would be publishable on its own.

Vectors C6, D26 are slightly heavier (need careful Pfaffian formula
implementation; LTC encoding & subsampled tester implementation
respectively); they are second-tier picks.
