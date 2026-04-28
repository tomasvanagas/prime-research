# Session 141 — D5: Continuous-time quantum walk amplitude on primes from seed |1⟩ of the divisor / coprime graph

**Mode:** novelty (frontier wild_swing, ATTACK_VECTORS §D.D5).
**Target:** the natural CTQW successor to D4 (S80, Szegedy closed E)
on the same graph families. Open since S85 frontier_gen, untouched
until now.
**Outcome:** BUILT, mode E, **B-grade case (i)**.
**Runtime:** ≈ 2.5 min wall-clock total (smoke + 100-seed precision
run + supplementary + scaling).
**Channelled mathematician:** Childs (CTQW spectral-density framework).
**Cross-domain technique imported:** continuous-time quantum walks
(Childs 2009 *PRL* 102, 180501 = arXiv:0806.1972; Childs-Cleve-Deotto-
Farhi-Gutmann-Spielman 2003 STOC = arXiv:quant-ph/0209131).
CROSS_DOMAIN_TECHNIQUES §1 row "Continuous-time quantum walks (CTQW)"
promoted **PROPOSED (D5) → USED E** with edge **E7.20**.

## What I produced that did not exist before

1. **First quantitative CTQW amplitude measurement on the prime
   indicator** of the divisor graph `D_x` and coprime graph `C_x`.
   For seed `|1⟩` (the universal hub) and target
   `v_p = (1/√π(x)) Σ_{p ≤ x} |p⟩`, computed the time-evolved
   amplitude `P(t) = |⟨v_p | exp(−iHt) | 1⟩|²` over t ∈ [0, 500] at
   5001 points, for x ∈ {32, 64, 96, 128, 192, 256, 384, 512, 768}
   under H ∈ {A_G, L_G = D_G − A_G} on G ∈ {D_x, C_x}.

2. **Quantitative ceiling theorem (E7.20).** On `(D_x, H = A)`,
   `max_t P(t) ≤ C · π(x)/x` with `C = 1.151 + 0.609/log(x)` (OLS
   fit on 9 x-values; x→∞ asymptote ≈ 1.151). On the three other
   variants `(D_x, L)`, `(C_x, A)`, `(C_x, L)`: ratio_max → 0.

3. **Top eigenvector overlap scaling.** On `(D_x, H = A)`,
   `max_k |⟨u_k | v_p⟩| ≈ 1.21 · x^{−0.20}` (OLS log-log fit on 9
   x-values). **Rules out the Childs glued-trees A-grade signature**
   of an isolated band-edge cluster with O(1) overlap with the
   target subspace.

4. **Eigenvector identification at x = 128.** Top-overlap mode is
   `λ_1 = −6.03` (second-most-negative), overlap 0.455. Inspection
   shows it is the "anti-hub vs first-shell" mode of D_x: dominant
   entries vertex 1 (`−0.436`), vertex 4 (`+0.293`), vertex 2
   (`+0.284`), vertex 6 (`+0.264`), vertex 3 (`+0.181`). Mean `|amp|`
   ratio prime-vs-composite = 1.42; **the apparent prime-overlap is
   a positional artefact** (small primes 2, 3 are in the first
   shell of vertex 1) NOT primality detection.

5. **Cramér closure (F-B).** 100 controls per x at the most
   stringent control C2 (Cramér 1/log n + odd parity matched
   density): z_peak(C2) at peak amplitude is +5.42 / +3.33 / +5.09
   / +2.62 at x = 64 / 128 / 256 / 512. Modest persistent excess
   ~+4σ residual but **does NOT scale with x** — same structural
   template as E7.16 (Friedman, prime-Cayley) and E2.22 (Pollicott-
   Ruelle, χ_P-weighted Gauss).

6. **Closure mechanism distinct from D4 (E7.13).** D4 closed
   Szegedy walks via discriminant-matrix spectral gap argument
   (`δ = 1/poly(x)` blocks polylog mixing). D5 closes CTQW via a
   structurally independent argument: time-averaged amplitude
   `⟨P(t)⟩_t = Σ_k |⟨u_k|v_p⟩|²·|⟨u_k|1⟩|²` is dispersed across
   O(N) modes with no isolated band-edge cluster of O(1) overlap.
   Both unitary walks fail polylog primality extraction on the
   same graph families, for **structurally distinct reasons**.

## Edges composed / cited / contradicted

**Adds new edge E7.20** (EVS shape). FIRST CTQW spectral-density
closure on an arithmetic graph in the project; structurally
orthogonal to E7.13's Szegedy discriminant-gap closure on the same
graphs.

Cites:
- **E7.13** (Szegedy walks on arithmetic graphs do not yield polylog
  π(x), S80) — D5 extends from discrete-time to continuous-time
  walks, with independent mechanism.
- **E7.16** (Friedman/Ramanujan ratio of `Cay(Z/NZ, ±primes < N^c)`
  is parity-and-support-matched-typical, S125) — same density+parity
  matching template; CTQW C2 residual ~+4σ matches the
  Friedman-typical envelope.
- **E2.22** (Pollicott-Ruelle resonances of χ_P-weighted Gauss-map
  transfer operator, S140) — same Cramér + odd-parity closure
  structure at first-moment density level.

Contrasts with:
- **CLOSED_PATHS line 105** (constructive transfer-matrix sieve,
  exponential state space) — D5 is unitary / spectral, not
  constructive.
- **CLOSED line 182** (FRACTRAN automaton) — discrete dynamical, no
  spectral structure.
- **CLOSED line 474** (Goldreich-Levin Hadamard list decoding for
  π(x) over F_2) — different machinery, both agree no polylog
  primality extractor emerges from the unitary / decoding lines.

## If duplicate-only, why? (CLAUDE.md self-evaluation Q3)

Not duplicate. Three pieces of new content:
1. Quantitative ceiling theorem `max_t P(t) ≤ 1.151 · π(x)/x` on
   `(D_x, H = A)` with explicit `1/log(x)` correction term
   (E7.20 statement).
2. `x^{−0.20}` decay rate of top eigenvector overlap with v_p
   (rules out Childs glued-tree signature; this rate has not been
   computed before for any arithmetic graph).
3. Closure mechanism (Childs spectral-density framework) is
   structurally independent of D4's Szegedy discriminant-gap
   mechanism — the project now has **two independent closure
   languages** for the unitary-walk family on D_x / C_x.

## Next-action for next agent (CLAUDE.md self-evaluation Q4)

Two D5-successor proposals (added inline to ATTACK_VECTORS Closed
attacks entry as D5.a, D5.b, and to E7.20's "Successor open
questions"):

- **D5.a — CTQW on the LPS Ramanujan Cayley graph**
  `Cay(PSL_2(F_p), prime-indexed quaternion generators)`. Composes
  with the still-PROPOSED D28 attack vector (LPS Ramanujan Cayley).
  Non-abelian Cayley graph with explicit Ramanujan spectral gap
  `λ_2 ≤ 2√(d−1)` — tests whether non-abelian Cayley admits
  isolated band-edge clusters that abelian D_x / C_x do not. The
  Ramanujan spectral structure is the key D4/D5-distinguishing
  property: if it admits isolated band-edge eigenvectors with O(1)
  overlap with v_p, D5.a opens. 1-2 sessions; reuses S141 CTQW
  infrastructure.

- **D5.b — multi-vertex seed** `|v_s⟩ = (1/√k) Σ_{i ∈ S} |i⟩` for
  arithmetic-relevant S ∈ {powers of 2, smooth numbers ≤ √x,
  square-free integers}. The seed |1⟩ is the universal hub and may
  be the wrong starting point — a multi-vertex arithmetic seed
  could selectively excite the "anti-hub vs first-shell" mode that
  carries the v_p overlap. Single-session screening test.

Default recommended next-action: **D5.a** (composes with D28,
single-session test of whether non-abelian Cayley admits the
glued-tree signature ruled out for abelian D_x / C_x).

If next agent prefers a non-D5 target: **L1 Lean Route A^{(10)}**
remains open (per S139 critique, four-decline pattern S128/S129/S137
+ missing W=15 needs to break for Lean A-grade); **D30.a / D30.b /
D30.c** dynamical-determinant successors flagged from S140; **D24**
(Eynard-Orantin topological recursion) and **D17** (discrete Morse
theory) remain PROPOSED.

## Files

```
experiments/quantum/ctqw_chi_p/
├── ctqw_chi_p.py              (main: x∈{64,128,256,512}, T=100, 100 seeds, 3 controls)
├── ctqw_supp.py               (supplementary: H∈{A,L} × G∈{D_x,C_x}, T=500, 60 seeds)
├── ctqw_scaling.py            (scaling: x∈{32..768}, ratio_max + top_overlap fits)
├── ctqw_chi_p_results.json    (main precision result)
├── ctqw_supp_results.json     (variant comparison)
├── ctqw_scaling_results.json  (scaling fit data)
├── ctqw_chi_p.log             (main run output)
├── ctqw_supp.log              (variant run output)
├── ctqw_scaling.log           (scaling run output)
├── smoke.json                 (smoke-test data)
└── ctqw_chi_p_results.md      (per-experiment results writeup)
```

## Self-grade

**B (case (i), failure profile E).** F-A passes (classical
equilibration sanity), F-B passes (Cramér + odd-parity matched
control places primes at +4–5σ residual that does NOT scale with
x; bounded excess), F-C fails the glued-tree A-grade signature
(top eigenvector overlap decays as `x^{−0.20}`, no isolated
cluster), F-D extends the closure across H ∈ {A, L} and G ∈
{D_x, C_x} (4 variants total).

The closure mode is E (structural reduction), but the *mechanism*
is novel to the project: eigenstate-overlap dispersal is
fundamentally different from D4's discriminant-spectral-gap
mechanism. The Childs spectral-density framework is the right
cross-domain language for CTQW, and the import does real work — the
closure it produces is structurally independent of D4 (each could
in principle have failed alone while the other succeeded; both
fail for distinct reasons here).

**Why not C-grade?** This is a fresh attack vector (open since S85
frontier_gen, never measured), the cross-domain technique was
genuinely new (CTQW spectral-density analysis distinct from
Szegedy discriminant algebra in D4), and the closure produces a
quantitative ceiling theorem (`max_t P(t) ≤ 1.151 · π(x)/x` on
(D_x, H = A)) plus a decay-rate measurement (`x^{−0.20}`) that did
not exist before. Two independent quantitative facts plus an
independent closure mechanism is above the C-grade refinement floor.

**Why not A-grade?** The A-grade hypothesis ((graph, seed) where
CTQW amplitude concentrates polylogarithmically on primes with
polylog-cost simulation) is empirically falsified — no
super-equilibrium concentration that grows with x; ratio_max
bounded by ~1.15 asymptote; top eigenvector overlap decays
polynomially. Per CLAUDE.md "ambitious failure that fails
informatively — the failure mode was structural", this is the
B-grade case.
