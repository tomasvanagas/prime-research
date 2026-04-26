# TODO

> **For agents starting a session:** read this file first, then `EDGES.md`,
> then `status/OPEN_PROBLEMS.md`. EDGES catalogues every real mathematical
> edge found across 58 sessions and groups them into the two surviving
> chains (E and F). Every other chain (A, B, C, D, G) has been closed.
>
> **Critical priority is FOCUS-1 below.** It is the only direction whose
> resolution would substantively change project status. FOCUS-2 / FOCUS-3
> are sharpening targets; FOCUS-4 is a recurring lightweight watch; the
> housekeeping and benchmarking items are non-blocking.
>
> Recently closed and removed from this file: Connes operator scaling
> (S53), pi(x) mod q for q in {2..13} via Liouville (S55+S56), Borel+Fejer
> hybrid (S52), 3-point GUE correlation (S57), N=2000 zeta structure (S45),
> Helfgott-Thompson M(x) benchmark (S60, sieve route shown tight via E6.7+E7.6).
>
> S60 EDGES extension added 10 new edges (E1.7, E1.8, E2.11, E2.12, E3.12,
> E6.7, E6.8, E6.9, E7.8, E7.9). Effect on this file: FOCUS-2 gains an
> E2.11 pre-disqualification test (5s vs 30s per candidate); FOCUS-5 gains
> level-of-distribution as a tracked literature category; the H-T benchmark
> is removed. FOCUS-1 / FOCUS-3 / FOCUS-4 unchanged.

---

# === CRITICAL PATH ===

## [FOCUS-1] AKS growing-dim MPOW in TC^0  [Chain E, the only open frontier]

**Why this is THE critical item.** After S50-S58 closed Chains A/B/C/D/G,
Chain E (TC^0 via growing-dimension matrix powering) is the only
remaining EVS-H surviving chain (`EDGES.md` §8). The missing primitive is
either:

1. resolve `growing-dim MPOW in TC^0?` (positive → PRIMES in TC^0),
2. or find an alternative TC^0 counting primitive (4th encoding,
   subsumed by FOCUS-2).

S47 closed cyclotomic-CRT splitting (AKS-prescribed r is prime in 21/22
sampled n). Three uncovered sub-attacks remain — each tractable for a
focused-mode session:

### Sub-attack 1: Smaller-r AKS variant (Bernstein 2003)

Lenstra-Pomerance proved AKS works at r = O(log^4 n). Bernstein 2003
strengthened the gcd condition. **Question:** can r = O(log^2 n) be
forced under any sufficient condition that itself sits in TC^0?

**Concrete next step:** read Bernstein 2003 (cr.yp.to/papers/aks.pdf),
identify the exact gcd condition, prototype the matrix-powering circuit
for the strengthened test, measure dim/r ratio over the same 22-sample
n grid as S47.

Save under `experiments/circuit_complexity/aks_alternative/bernstein_smaller_r/`.

### Sub-attack 2: Non-cyclotomic ring decomposition

AKS works in `Z_n[x]/(x^r - 1)`. Replace the modulus by an irreducible
`f(x)` of degree exactly polylog(n) over `Z_n`. The cyclotomic-CRT
split (S47) closes because `x^r - 1` factors trivially when r is prime;
choosing `f` irreducible bypasses that.

**Concrete next step:** prototype with `f(x) = x^d + a` (Eisenstein-style)
for d = ceil(log^2 n), check (a) irreducibility certification cost in
TC^0, (b) whether the AKS congruence still characterises primality
when the ring is `Z_n[x]/f(x)`.

Save under `experiments/circuit_complexity/aks_alternative/non_cyclotomic_ring/`.

### Sub-attack 3: Healy-Viola Frobenius transplant

Healy-Viola 2006 showed `F_{2^n}` exponentiation is in TC^0 by
exploiting that Frobenius `x -> x^2` is the identity on `F_2`-coefficients.
**Question:** does any partial Frobenius-like map exist on `Z_n[x]/(x^r - 1)`
that gives the same depth collapse?

Specifically: if we work modulo a prime `q` dividing `n - 1`, the
`q`-th-power Frobenius on `(Z/qZ)[x]/(x^r - 1)` IS a ring homomorphism.
Whether AKS-correctness lifts back to `Z_n` from this quotient is the
question.

**Concrete next step:** compute the full power-Frobenius orbit structure
on `(Z/qZ)[x]/(x^r - 1)` for small (n, q, r) triples, see if AKS
verification can be staged through a TC^0-friendly chain.

Save under `experiments/circuit_complexity/aks_alternative/frobenius_transplant/`.

### Backup objective if all three close

Then the sub-attack space is exhausted and the only escape is **Brandt
MKtP** or a fundamentally new lower-bound technique. Mark Chain E as
"computationally cornered" in `status/OPEN_PROBLEMS.md` and shift
remaining sessions to FOCUS-2 / FOCUS-3.

---

# === SECONDARY (sharpening) ===

## [FOCUS-2] 4th encoding of pi(x) — concrete candidate sweep

E7.7 (three-pillars meta-theorem) says only 3 informationally-complete
encodings of pi(x) are known: prime positions, zeta zeros, floor values.
15+ candidate intermediate-quantity families have been individually
closed (S15, S16). The *space of additive number-theoretic functions*
has not been enumerated systematically.

### Pre-test (fail-fastest, **NEW from S60 EDGES extension**)

E2.11 pins R(x) as the **exact** smooth/random separator: iterated finite
differences of `f(x) = pi(x) - R(x)` grow as white noise with RMS ratio
→ 2.0 per Δ-application; Hankel rank is full (250/250). So:

> **Pre-disqualification:** for each candidate `T_i`, compute its known
> leading asymptotic `A_i(x)`. Iterate `Delta^k (T_i - A_i)` for k = 1..7
> on x ∈ [10^4, 10^6]. **If RMS ratio Delta^{k+1}/Delta^k → 2.0**, the
> residual is GUE-type and `T_i` is just another encoding of the same
> information — close immediately, no PSLQ run needed (~5 seconds, not 30).
> Only candidates whose residual difference operator is **non-white**
> (ratio bounded away from 2, or annihilated at finite order) advance to
> the PSLQ stage.

Expected outcome from the pre-test: candidates 1, 2, 4, 5, 6 below all have
closed-form leading terms whose residuals empirically look GUE-type — they
should pre-disqualify in seconds. Only candidate 3 (Psi(x, log^c x))
plausibly survives because the smooth/random split for *smooth-number
counts* is genuinely different from R(x)'s split.

### Concrete candidates (test in this order, fail-fast)

1. **Sum of log Gamma fractional part:** `T_1(x) = sum_{n<=x} {log Gamma(n)}`.
   Distinct from psi(x); ties to Stirling's series. Known asymptotic but
   no polylog inversion attempt logged.
2. **Harmonic-number summatories:** `T_2(x) = sum_{n<=x} H_n`,
   `T_3(x) = sum_{n<=x} H_n^2`. Both have closed-form leading terms;
   residuals never PSLQ'd against pi(x).
3. **Smooth-number count Psi(x, B) for varying B:** S26 partially tested
   B = sqrt(x); B = log^c(x) regime untested. **Most likely to survive
   the E2.11 pre-test** — the Dickman/Buchstab piece is structurally
   different from R(x).
4. **Divisor function sigma_k summatories (k = 1, 2, 3):** the
   Dirichlet series structures are different from zeta(s); cross-PSLQ
   never run against pi(x) - R(x).
5. **Squarefree-count Q(x):** ratio Q(x)/x converges to 6/pi^2; the
   residual `Q(x) - 6 x / pi^2` was tested (closed S21) but not against
   pi(x) directly via Mobius inversion.
6. **Jacobi totient phi(n) summatory:** `T_6(x) = sum_{n<=x} phi(n)`.
   Closed-form leading 3 x^2 / pi^2; residual periodicity untested.

For each candidate that passes the pre-test:
- Compute `T_i(x)` exactly for x up to 10^6.
- PSLQ at 60 dps against `{pi(x), pi(x)/log(x), li(x), R(x), x, sqrt(x)}`
  with cross-validation at 4 distinct x values.
- If any relation survives cross-validation: deep dive.
- If none does: log closure with three-pillars verdict.

Save under `experiments/algebraic/fourth_encoding_search/<candidate>/`.

**Each pre-test is ~5 seconds; full PSLQ is ~30 seconds.** With the E2.11
pre-test, six candidates take well under one focused session. Most likely
outcome: 5/6 pre-disqualify, 1/6 (Psi) reaches PSLQ and closes there,
but the search-space narrowing is itself valuable.

## [FOCUS-3] Brandt MKtP framework deep dive

S30 flagged Brandt (TCC 2024, MKtP not in DTIME[O(n^2)]) as "the only
known technique that could lead to unconditional superpolynomial lower
bounds without hitting the barriers." S39 confirmed no follow-up papers
through April 2026. **Project has never engaged with the technique
beyond logging it.**

### Concrete actions

1. Read Brandt 2024 (TCC, alternatively arXiv listing) carefully.
2. Identify whether the diagonalization argument can be applied to
   pi(x) mod 2 specifically — the function is total, computable, has
   sub-exponential time bound (O(x^{2/3})), and is conjectured outside
   TC^0.
3. Even a *conditional* Brandt-style lower bound on pi(x) mod 2 would
   be the first non-trivial circuit lower bound for our problem.
4. Save analysis at `experiments/circuit_complexity/brandt_mktp/`
   (theory dump + any concrete computations).

This is exploratory — most likely closure mode is "Brandt's hypothesis
doesn't extend to natural functions", but the upside is unique among
remaining options.

## [FOCUS-4] Higher-order zeta-zero correlation at N >= 10^4

S57 closed 3-point correlation at N=2000 to GUE noise floor, but
explicitly noted: *"arithmetic corrections [Conrey-Snaith] become
visible only at large lag in the Fourier-conjugate variable, not here."*
The S57 windows-disjoint test at L=32 gave c_3 ~ 10^-3, but the
disjoint-window count was only ~62 windows.

### Concrete action

1. Use Odlyzko's tabulated zeros (online) to extend N from 2000 to
   10^5 (factor 50 in sample).
2. Re-run S57's pair / triple / cumulant battery with proper disjoint
   windows.
3. Specifically test the **Conrey-Snaith arithmetic correction** for
   pair correlation:
   `R_2_arith(s) = R_2_GUE(s) + (1/2) (log(gamma_n / 2 pi))^{-2} A_HL(s)`
   where `A_HL` is a Hardy-Littlewood-style correction.
4. If the arithmetic correction is empirically confirmed at this scale,
   it is genuinely novel structural info (not in current literature
   for pi(x)) and could feed back into Chain D-style attacks if anyone
   ever revisits them.

Save at `experiments/analytic/zeta_structure/large_n_correlations/`.

## [FOCUS-5] Literature watch (recurring, lightweight)

Last run: S58 (window 2026-04-05 -> 2026-04-26, NO-DELTA).
See `archive/sessions/session58_literature_watch.md`.

Re-run every 2-3 weeks. Update `literature/state_of_art_2026.md` only if
the asymptotic landscape changes.

Sources to scan:
- arXiv math.NT recent submissions for "pi(x)", "p_n", "zeta zero"
- ECCC TR2026 for TC^0/NC^1 separations and **Brandt MKtP follow-ups**
- Connes / Yakaboylu / van Ittersum / Ono author streams
- primecount / kim-walisch GitHub for new releases
- **Unconditional level-of-distribution past x^{5/8} (NEW from S60 EDGES
  extension).** Pascadi 2025 (E3.12) reached x^{5/8} = 0.625 unconditionally,
  within 0.04 of the algorithmic Meissel-Lehmer threshold x^{2/3} ≈ 0.667.
  Any unconditional improvement past x^{2/3+epsilon} would convert several
  conditional pi(x) algorithms (Lagarias-Odlyzko under GRH, etc.) into
  unconditional ones — this is a publishable status-change for our project.
  Watch for: Pascadi follow-ups, Maynard streams, Bombieri-Vinogradov
  improvements, Friedlander-Iwaniec extensions.

Note: most months produce zero algorithmic deltas. Absence of news is
itself information.

---

# === HOUSEKEEPING & BENCHMARKING (non-blocking) ===

## [HOUSEKEEPING] Flag duplicate scripts (human review required)

Rule 11 says no duplicate scripts. Suspects flagged S39, awaiting human review:

- [ ] `experiments/analytic/weil_optimized.py` vs `weil_optimized_v2.py` vs `weil_optimized_v3.py` (v3 same size as v1, possible revert)
- [ ] `experiments/sieve/ht_transfer_attempt.py` vs `ht_signed_transfer_v2.py`
- [ ] `experiments/circuit_complexity/k_party_nof.py` vs `k_party_nof_v2.py`
- [ ] `experiments/circuit_complexity/approx_degree.py` + `approx_degree_quick.py` + `approx_degree_small.py` + `approx_degree_prime.py` + `approx_degree_counting.py`
- [ ] `experiments/information_theory/information_shortcut.py` vs `information_shortcut_v2.py`
- [ ] `experiments/other/breakthrough_attempt_v2.py` (no v1 exists)
- [ ] `experiments/other/self_correcting_v2.py` (no v1 exists)

Each pair has a companion `_results.md`. Do NOT delete without human approval.

> **Removed S60:** `[BENCHMARK] Helfgott-Thompson O(x^{3/5}) for M(x)`
> dropped from this file. E6.7 (HKM time-space curve O~(N^{8/15}) /
> O~(N^{1/3})) combined with E7.6 (Lucy-DAG pebbling lower bound
> Omega(x^{5/6}/ln x)) shows the sieve route is asymptotically tight to
> within x^{0.034}. H-T sits inside this Pareto frontier and benchmarking
> it would not change project status. The script `experiments/sieve/
> mertens_speedup.py` remains in place for archival reference.
