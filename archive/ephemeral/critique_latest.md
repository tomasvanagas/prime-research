# Critique — Session 51 — 2026-04-26

Proposals reviewed: A (Tropical Hankel rank of π), B (CF of prime
constant), C (Algebraicity of zero differences), D (Reservoir computer
fit on π(n)). Both experiments executed by the proposer (B, C) returned
clean nulls. After cross-referencing against `status/CLOSED_PATHS.md`
(733 lines, 720+ tested approaches) and `EDGES.md`, the verdict on
each is below.

---

## Proposal A — Tropical Hankel Rank of π

### Match in CLOSED_PATHS

- **Line 199** (S4): "Tropical geometry — FAIL/C — wrong compression"
- **Line 318** (S16): "Tropical geometry for π(x) — FAIL/I+E —
  Tropicalization loses all info (min=smallest prime); floor values ARE
  tropical objects; no new content"
- **Line 423** (S16): "Tropical geometry of P(s) — FAIL/I+E — tropical
  convolution = standard convolution"
- **Line 424** (S25+): "Tropical/min-plus full battery (7 tests) —
  FAIL/I — tropical Dirichlet conv, Mertens, sieve, Euler product,
  p-adic val of x!, **tropical det**, tropical conv — all fail. Core
  reason: min-plus is OPTIMIZATION not COUNTING; val(a+b) ≥
  min(val(a), val(b)) destroys sum info."
- **Line 652** (S32): "Tropical/min-plus fast-forward via gap structure
  — FAIL/I — Hankel rank 418/500 (ratio to random: 0.98). Min-plus
  accuracy 7-11%. **Algebraic obstruction: p(n)=Σ gaps needs (+,×) not
  (min,+)**."

### Verdict: **DUPLICATE**

The proposal frames it as Hankel-of-π in (max,+), which is a permitted
narrow specialisation, but:

1. The general "tropical battery" (line 424) already includes tropical
   determinant — and tropical rank is a relaxation of the Hankel-tropical-
   determinant test that closure already ran.
2. Line 652 already measured the Hankel rank of the prime gaps (418/500
   ≈ random) and closed the min-plus fast-forward route. Hankel of π is
   the cumulative integral of Hankel of g (gaps); rank can only differ
   by 1 (adding a rank-1 prefix-sum operator).
3. The structural obstruction in line 424 — "tropical = optimization,
   not counting" — applies to the proposed rank in the same way: a
   max-plus decomposition $\pi(i+j) = \max_k(u_{ik}+v_{kj})$ gives a
   piecewise-linear MAX, not a count, and does not yield a primality
   evaluator.

Re-running CGB tropical rank on H[i,j]=π(i+j) at m=2000 would almost
certainly reproduce the line-652 finding (rank ratio ~ 1.0 vs random),
adding nothing to the catalogue.

### Failure mode

**E (equivalence)** for the rank measurement, **C (circularity)** for
the proposed evaluator (building H requires π up to 2N-2 already).

---

## Proposal B — Continued Fraction of Prime Constant

### Match in CLOSED_PATHS

- **Line 615** (S29): "CF / Stern-Brocot structure of p(n)/n — FAIL/I —
  CF partial quotients follow Gauss-Kuzmin (random)"
- **Line 288** (S8): "Continued fraction of p(n)/n — FAIL/I — Higher CF
  terms unpredictable"
- **Line 298** (S4): "BBP digit extraction — FAIL — No BBP formula for
  prime constant"
- **Line 588** (S24): "k-automatic prime indicator — FAIL/I — 2-kernel
  has 38 growing subsequences (vs 6 for Thue-Morse). NOT k-automatic
  for any k"

The specific object — **regular CF of α = Σ_{p prime} 2^{-p}** — is
slightly different from "CF of p(n)/n" (line 615) and does not appear
verbatim in CLOSED_PATHS. So the *measurement* was novel-in-object even
though it was an obvious extension.

### Verdict: **DUPLICATE-PLUS** (object new, mechanism predicted)

The author's own experiment confirms what line 615 predicted at the
sequence level: regular CF preserves entropy and exposes Gauss-Kuzmin
statistics with no compressible deviation:

| measurement | empirical | theoretical | comment |
|---|---|---|---|
| K₀ (Khintchine geomean) | 2.7336 | 2.6854 | within 1.8% |
| Lévy log q_n / n | 1.2056 | 1.1866 | within 1.6% |
| max autocorr lag 1..10 | 0.013 | 0 | noise floor |
| 1/f^β slope | −0.091 | 0 | white |
| k-automaticity ratio (k=10) | 0.011 | 0 | noise floor |

This is a textbook **Khintchine-typical** signature. The negative
result is high-information at the *measurement* level (it adds α to
the list of computable reals empirically indistinguishable from generic
in CF statistics, a 33rd pseudorandomness measure for the project) but
it does NOT reopen anything.

### Failure mode

**E (equivalence)** — switching from binary to regular CF is a measure-
preserving bijection on irrationals; entropy is preserved; CF-shortcut
to a_n in polylog requires automaticity which is empirically absent.
Adjacent to the **k-automatic prime indicator** closure (line 588, S24,
"2-kernel has 38 growing subsequences vs 6 for Thue-Morse — NOT
k-automatic for any k").

### Recommendation

Add as a new CLOSED_PATHS entry (object is new even though mechanism
isn't), and append the CF result to `novel/pseudorandomness_of_pi.md`
as a CF-side variant of the existing k-automaticity closure. Do NOT
escalate to standalone `novel/` — the result is the absence of
structure, not the discovery of one.

---

## Proposal C — Algebraicity of Riemann Zero Differences

### Match in CLOSED_PATHS

- **Line 23** (S25): "PSLQ linear relations among zeros — FAIL/I —
  **13,000+ tests at 60-digit precision**: zero relations among
  subsets of 3-5 zeros with {1, π, log(2π)}. **1,225 pairwise tests
  negative**. Large-K hits are lattice artifacts (random baseline
  confirms). **Zeros linearly independent over Q.**"
- **EDGES E7.1** (S25, S45, S20, S57): "PSLQ on first 1000-2000 zeros:
  **0/1225 pairs, 0/4060 triples, 0/400 cross-block** integer relations
  at 30-60 digit precision."

### Verdict: **DUPLICATE**

The proposal's experiment runs PSLQ on $\gamma_k - \gamma_1$ for
k = 2..200 against a 12-element basis of arithmetic constants
{1, π, log 2, log 3, log 5, log 7, log 11, log 13, log Γ(1/4), log Γ(1/3),
log Γ(1/6), log π} at 50 dp with height cap 10^15. Result: 0/199 hits.

This is a **strict subset** of the line-23 / E7.1 coverage:
- E7.1 ran 1225 pairwise + 4060 triple + 400 cross-block tests
  (≈ 5685 tests) at higher precision (60 dp) on a smaller arithmetic-
  constants basis.
- The only "novelty" is enlarging the basis from {1, π, log(2π)} to a
  12-element set — but Γ-values at small denominators and log p_k for
  p_k ≤ 13 are already covered under the E7.1 footprint.

E7.1's verdict is rated **shape-revealing** in EDGES, meaning *every
zero is independent information*. A 200-test, 12-basis re-run at lower
precision cannot move that needle.

### Failure mode

**I** — information is irreducibly distributed across zeros, not
concentrated in a finite-rank algebraic shadow.

### Recommendation

Do NOT add a new CLOSED_PATHS entry. The result strictly subsumes
under line 23 and E7.1. The author's own results file already notes
the verdict is consistent with folklore.

---

## Proposal D — Reservoir Echo of π(x) via Liquid State Machines

### Match in CLOSED_PATHS

- **Line 683** (S43): "Reservoir computing on δ(n) (Dynamical
  Compressibility) — FAIL/I — ESN is a strict subset of recurrent
  dynamical systems; transformer (line 149, 1.1% accuracy), FNO, GP,
  PSLQ symbolic regression on δ(n) all already failed. δ(n) Kt entropy
  5.58 bits/value (line 555). **Polylog reservoir state cannot encode
  irreducible info.** Symbolic-regression readout duplicates line 146."

### Verdict: **DUPLICATE**

Line 683 closed reservoir computing on δ(n) = p(n) − R⁻¹(n). Proposal
D fits a reservoir to π(n) directly (not the residual). This is a
weaker target — the smooth part R(x) is a $C^∞$ function the reservoir
*can* fit by polynomial/PNT-baseline interpolation, so a small-L
reservoir might pass the proposed MAE = 0 test on n ∈ [5001, 10000]
purely by tracking the smooth envelope. The interesting question is
the residual, which is exactly what line 683 tested and closed.

If the experiment were re-run on π(n) itself it would either:
1. **Pass at large L** by memorising — exactly the failure mode line
   683 calls out for transformers (1.1% accuracy on δ at line 149
   scaled up).
2. **Pass at L = poly(log n)** by tracking only R(x) — already known
   information-free (R⁻¹ gets ~50% of digits, see CLAUDE.md).
3. **Fail extrapolation** — reproducing line 683.

None of these is a new outcome.

### Failure mode

**I (information loss).** Polylog reservoir state has $O(\log^c n)$
bits of memory but π(n) requires the same ≥ 50% of $n$ irreducible
bits as δ(n) does (the smooth half is just R(x) and is removable in
polylog).

### Recommendation

Do NOT escalate. Do NOT add to CLOSED_PATHS (already covered by line
683). Remove from "next experiments" list.

---

## Summary table

| proposal | object | verdict | matching CLOSED_PATHS | action |
|---|---|---|---|---|
| A | tropical Hankel of π | **DUPLICATE** | 199, 318, 423, 424, 652 | none |
| B | regular CF of α=Σ 2^{-p} | **DUPLICATE-PLUS** (object new) | 288, 298, 588, 615 | add CLOSED_PATHS entry, append measure to `novel/pseudorandomness_of_pi.md` |
| C | PSLQ on γ_k − γ_1 | **DUPLICATE** | 23 / E7.1 | none (strict subset) |
| D | reservoir fit on π(n) | **DUPLICATE** | 683 | none |

## Meta-observation

The proposer flagged C as the "highest signal yes/no" but it was the
*most-replicated* closure in the pre-existing literature (E7.1 +
13,000+ tests). The cleanest novel-in-object result is B, but its
mechanism — that representation-changing bijections preserve entropy —
is the **three-pillars meta-theorem** (E7.7). Proposals A and D
re-traverse closed analytical paths (tropical, reservoir) without
addressing the structural obstructions on file.

**Underweighted directions** the proposer did not touch and which the
project explicitly requests in CLAUDE.md S60 onward:
- Construction of new circuits / algebraic objects / transforms (the
  three FOCUS-1 sub-attacks; only one — Bernstein smaller-r — is
  already built and closed at line 722; FOCUS-3 Brandt unbuilt).
- Lower-bound techniques — Brandt MKtP-style diagonalisation against
  TC^0; pebbling on non-DAG models.
- Theoretical sharpening of existing novel results (rigorous N/2
  threshold proof, Lean formalisation of MPS bond-dim theorem).

Future "fresh, no prior context" proposal sessions should be cross-
referenced against CLOSED_PATHS sections **before** the experiment is
designed, not after — the project is mature enough that the duplicate-
detection step is the highest-value activity.
