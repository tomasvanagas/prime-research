# Novelty Challenges — Active Research Targets

This file holds the active research questions for sessions targeting
**B-grade work** (substantive refinement or ambitious failure of a
single-session attack — see CLAUDE.md "Novelty Bar").

**For A-grade work, read `ATTACK_VECTORS.md` first.** That file holds
the frontier targets — the only way to produce genuinely novel
mathematics. The challenges below are ranked B-grade because they are
single-session-tractable; the trade-off is that even when they succeed,
they produce refinements rather than discoveries.

**Pick one challenge per session.** Time-budget the work. If you can't
finish, write the partial state into RESEARCH_AGENDA.md and halt.

The framework's grading discipline (CLAUDE.md): a session that
completes a NOVELTY_CHALLENGES.md target with a clean refinement is
B-grade. A session that *attempts* an ATTACK_VECTORS.md target and
fails informatively is also B-grade. A session that completes 5
NOVELTY_CHALLENGES targets in one sitting is *not* A-grade — it's 5×B,
and a single A-grade attempt would have been more valuable.

---

## §0. Pointer to Frontier Attacks (A-grade work)

**`ATTACK_VECTORS.md` §A through §F** holds the project's frontier
targets. If your session has 1+ hour available and you are not in
the middle of an arc, pick from there instead. The frontier sections:

- **§A** TC⁰ Primality, Beyond AKS (3 attacks)
- **§B** Spectral Identity Searches Beyond Standard Bases (3 attacks)
- **§C** Break GUE Universality of Zeta Zeros at Some Scale (3 attacks)
- **§D** Cross-Domain Imports the Project Has Never Tried (4 attacks)
- **§E** Meta-Analysis of CLOSED_PATHS as a Data Source (2 attacks)
- **§F** Synthesis Targets That Would Be Publishable (3 multi-session arcs)

**S147 critique recommended next pick (highest A-prior open vector):**
**ATTACK_VECTORS.md §D.D31 — Adiprasito-Huh-Katz combinatorial Hodge
theory of an arithmetic prime-matroid.** S142 self-stated 25%
A-grade prior, AHK 2018 *Annals* 188 (= arXiv:1511.02888) machinery
fully mature, Heron-Rota-Welsh log-concavity of `|μ_{M_P^N}(F)|` is
an explicit testable inequality, flat enumeration feasible at N ≤ 32
with Bernoulli-matched-density control. If D31 closes E or I, pivot
to D32 (ACSV Pemantle-Wilson). Backup synthesis: a unifying
`novel/cramer_parity_envelope.md` covering E2.18 / E2.20 / E2.21 /
E2.22 / E2.23 / E7.16 / E7.20 — ≥ 8 measurements of the same
structural fact (primes are density+parity Cramér-typical at first
moment) is a defensible synthesis target.

**§D6 — CLOSED (S85, mode E, B-grade).** Gowers U^k norms of χ_P
empirically match the Hardy-Littlewood {0,1}^k-cube singular series
S_2 = 2.301 (Q^2(χ_P) = 2.153 at N=2^18, monotonically converging
from below) and S_3 = 54.12 (Q^3 ≈ 35.5 at N=2^13, slow finite-N
convergence). W=210 W-trick restores Gowers uniformity at U^2 within
0.1%. Outcome (b) of the original challenge — the 36th-37th
pseudorandomness measure with closed-form HL prediction. New EDGE
**E2.13**. See `experiments/information_theory/gowers_uk_chi_p_results.md`,
`archive/sessions/session85_d6_gowers_uk_chi_p.md`.

**Successor challenges (proposed in S85):**

**§D6.a — Gowers U^4 of χ_P at N ≤ 2^12.** {0,1}^4-cube has 16
forms; α_p(k=4) requires `(Z/p)^5` enumeration, borderline tractable
for p ≤ 23 with vectorisation. Predicted S_4 ~ 10^3 to 10^4 from
the per-prime factor pattern. Empirical Q^4(χ_P) verifying ~ S_4
would extend the Gowers-norm fingerprint chain by one level
(B-grade refinement). Cost: 1 session.

**§D6.b — Λ vs χ_P U^k comparison — CLOSED (S93, mode E, B-grade).**
Side-by-side measurement of Q^k(χ_P) and Q^k(Λ) on Z/NZ for N up to
2^17 (U^2) and 2^12 (U^3). Result: **the Hardy-Littlewood
{0,1}^k-cube singular series S_k is invariant under log-weighting** —
Q^k(Λ) and Q^k(χ_P) agree to within 0.4% at U^2 and 2.5% at U^3, both
converging to S_k from below at similar rates. After W-trick at
W = 210: Q^2(χ_W) = Q^2(Λ_W) = 1.0029 to four decimals (the W-trick
also erases the log-weight discrepancy). The residual +0.3% offset
Q^k(Λ) > Q^k(χ_P) is identified as the prime-power weight contribution
(Λ counts n = p^k for k ≥ 2), vanishing as π(√N)/π(N) → 0. **Refines
E2.13** to a statement about all natural prime-supported weightings,
not just the bare indicator. See
`experiments/information_theory/lambda_vs_chi_p_uk/lambda_vs_chi_p_uk_results.md`,
`archive/sessions/session93_d6b_lambda_vs_chi_p_uk.md`.

**Successor challenges (proposed in S93):**

**§D6.c — CLOSED (S101, mode E, B-grade).** μ-weighted χ_P at U^k.
The literal target collapses: μ(p) = −1 on every prime, so
`μ·χ_P ≡ −χ_P` pointwise on ℤ and `Q²(μ·χ_P) = Q²(χ_P) → S_2`. We
pivoted to the broader question "does Möbius cancellation propagate to
indicator level sets?" and tabulated 11 multiplicatively-defined
indicators on N ≤ 2^17. Sharp dichotomy:
(i) signed μ, λ are Gowers-uniform (Q²_norm → 1 ≡ Bernoulli baseline);
(ii) density-1/2 indicators 1[λ=±1] inherit Gowers-uniformity (Q² → 1.0000);
(iii) asymmetric-density indicators 1[μ=±1] (density 3/π²), sqfree
(density 6/π²), 1[Ω=2] retain HL structure (Q² ≈ 1.038 / 1.038 / 1.116);
(iv) W=210 W-trick uniformly collapses every indicator to Q² ≤ 1.0041.
**Refines E2.13** to family-level statement: HL-singular-series
structure is killable by W=210 W-trick across the whole multiplicative
indicator family, not just χ_P. New empirical constant
`S_2^{sqfree} ≈ 1.0384` for the squarefree singular series at U^2.
S87's Liouville-uniformity is now structurally explained as a
density-1/2 phenomenon, not a Liouville-specific one. See
`experiments/information_theory/mu_weighted_chi_p_uk/mu_weighted_chi_p_uk_results.md`,
`archive/sessions/session101_d6c_mu_weighted_chi_p.md`.

**§B1 CLOSED (S92, mode E, B-grade).** Polynomial method / slice rank
on χ_P closed via algebraic immunity AI(χ_P, N) = 2 with explicit
annihilator (1 + x_0)(1 + x_1) encoding the trivial mod-4 sieve fact;
W-trick at W ≥ 6 fully removes the deviation. New EDGE **E2.15**.
Third leg of the W-trick triple alongside E2.13 (Gowers) and E2.14
(Anderson). See `experiments/algebraic/algebraic_immunity_chi_p/`,
`archive/sessions/session92_b1_algebraic_immunity_chi_p.md`.

**Highest-leverage attempt now:** §C5 — Stein's method for quantitative
finite-x Gaussianity of `(π(x) - Li(x))/(√x/log x)` (S91 frontier_gen).
Single-session viable: precompute π(x) on `[10^6, 10^7]` via segmented
sieve, compute `D(x_k) = (π(x_k) - Li(x_k)) log(x_k) / √x_k` at K=1000
log-uniform anchors, compute empirical Wasserstein-1 distance
`W_1(D, N(0, σ̂²))` via sorted-CDF integration. Outcome shapes:
(a) `W_1 ≥ c > 0` even as K → ∞ → **A-grade**: FIRST quantitative
finite-x non-Gaussianity result for π(x) - Li(x);
(b) `W_1 ≈ 1/√K` matching Gaussian baseline → **B-grade**: 38th
pseudorandomness measure of strongest type. Cross-domain ingredient:
Stein's method (Stein 1986; Ross 2011 *Probability Surveys* 8, 210;
Chen-Goldstein-Shao 2011 Springer) — never applied to π(x) deviations
in the project (the existing "Stein" CLOSED_PATHS hits are about a
different Stein — modular form decay analysis).

**§D2 — CLOSED (S96, mode I, B-grade).** Persistent homology of
Takens-embedded normalised prime gaps deviates ≥ 5σ from BOTH IID
Exp(1) AND from gap-marginal-permuted control across d ∈ {2, 3, 4}
and M ∈ {500..4000}. PRIMES T0 (total H_0 persistence) z(B2) =
−8.70, T1 z(B2) = −11.99 at M=2000, d=3, x ≈ 10^6 (rank 0/20 in K=20
bootstrap each). Sign: PRIMES is MORE clustered and has FEWER
persistent loops than random — geometric self-similarity from HL
k-tuple admissibility. New EDGE **E2.17**; fifth orthogonal HL-
detection category after E2.13/E2.14/E2.15/E2.16, in algebraic-
topological / metric-space geometry. NOT an algorithmic opening —
VR-PH is O(M^3). See
`experiments/topological/persistent_homology_chi_p/persistent_homology_chi_p_results.md`,
`archive/sessions/session96_d2_persistent_homology_chi_p.md`.

**Successor challenges (proposed in S96):**

**§D2.a — CLOSED (S117, mode E, B-grade refinement of E2.17).**
W=210 W-trick on prime-gap sequence: filter to `q ≡ b (mod 210)` for
b ∈ {1, 11, 13}, Cramér-normalise as `x_n = g_n / (φ(W) log q_n)`,
re-run Takens+ripser PH. **Outcome: F_a holds for the serial-
correlation component.** z(W-tricked PRIMES; B2)_T0 = −1.99 / T1 =
−0.67 (M=1000, d=3, x≈10⁶, three residues pooled), versus S96
unconditioned z(B2)_T0 = −7.45 / T1 = −4.05 — collapse by 3.7×–6×
across all (M, d, x) cells tested. At x ≈ 5·10⁶ even cleaner:
|z(B2)| ≤ 1 on both summaries. The B1 (IID Exp(1)) deficit is
preserved/amplified because the W-tricked gap MARGINAL is non-Exp(1)
(gaps are multiples of W, giving a discrete spectrum). **E2.17
decomposes** as serial-correlation (HL k-tuple) + marginal-
distribution (gap-quantisation) components: W-trick kills the first
and structurally amplifies the second. **E2.17 refined inline** —
sixth leg of the W-trick HL-fingerprint family alongside E2.13 /
E2.14 / E2.15 / E2.16 / E2.20. See
`experiments/topological/persistent_homology_w_trick/persistent_homology_w_trick_results.md`,
`archive/sessions/session117_d2a_ph_w_trick.md`.

**Successor challenges (proposed in S117):**

**§D2.a.1 — CLOSED (S124, mode E, B-grade refinement of E2.17).**
PH on the residual marginal-distribution component. Adding B3 = IID
inverse-transform samples from the W-tricked empirical CDF (Devroye
1986 §II.2.1; *continuous* envelope without discreteness) to the
S117 baselines, at M=1000, d=3, x≈10⁶, three residues pooled:
z(P_W; B3)_T0 = −0.05, z(P_W; B3)_T1 = +0.46 (vs S117 B1: −9.07 /
+5.56; B2: −1.99 / −0.67). Robust across d ∈ {2,3,4}: |z(B3)| ≤ 0.65
on every cell. **Refines E2.17** to a three-component decomposition:
(i) **marginal-envelope component** (~7–9σ on T0; the W-tricked
marginal variance ≈ 0.55 vs Exp(1) variance 1, support shifted off
zero by gap-quantum 0.318) — dominant; absorbed by B3.
(ii) **discreteness component** (~1–3σ on T0; B2 mean > B3 mean,
discrete-grid permutation has higher T0 than continuous-envelope IID)
— sub-leading.
(iii) **residual serial-correlation component** (~1–2σ on T0;
S117's W-tricked gap-correlation residual) — sub-leading.
Components (ii) and (iii) partially cancel on the (PRIMES_W vs B3)
comparison, which is why z(B3) ≈ 0. Pre-stated F_a HOLDS on absolute
thresholds, partially fails on relative |Δz| ≤ 1 condition for d=2,3
(the failure IS the new content). See
`experiments/topological/persistent_homology_w_trick_marginal_b3/persistent_homology_w_trick_marginal_b3_results.md`,
`archive/sessions/session124_d2a1_ph_marginal_b3.md`.

**§D2.a.2 — CLOSED (S138, mode E, B-grade refinement of E2.17).**
Tabulated `z(B2; T0)` across W ∈ {2, 6, 30, 210, 2310} at M ∈ {500,
1000} pooling first min(3, φ(W)) coprime residues, otherwise matched
to S117 protocol. Clean monotone decay at matched M=500: 4.89 →
1.52 → 0.93 → 0.99 → 0.30. The HL closed-form fit
`r(W) = ∏_{p|W, p>2} (1 - α/p)` gives `α ≈ 2.07`, matching the
Hardy-Littlewood twin-prime per-prime local factor `1 - 2/p` exactly
on the W=6 → W=30 cell-pair (residuals 0.001 / 0.008). **The p=3
factor alone accounts for 70% of the deficit; by W=6 the serial
component is at the K=20 noise floor.** S117's W=210 anchor was in
the saturation regime, not the HL-active regime. M=1000 W=2310
rebound to z = -3.04 identified as finite-size window non-stationarity
(at φ(2310)=480 the M=1000 window spans log range 2.13; M=500 collapses
the window and the rebound disappears, +0.30). **E2.17 refined inline**
with per-prime decay rate `(1 - 2/p)`. Sixth leg of the W-trick HL
fingerprint family alongside E2.13/E2.14/E2.15/E2.16/E2.20 retains
its status; the new content is the **explicit per-prime decay rate
matching the twin-prime Hardy-Littlewood constant**. See
`experiments/topological/persistent_homology_w_scan/persistent_homology_w_scan_results.md`,
`archive/sessions/session138_d2a2_ph_w_scan.md`.

**Successor challenges (proposed in S138):**

**§D2.a.2.i — Higher-K noise floor at W ∈ {30, 210, 2310}.** At
K=20 baselines the noise floor is ~ 1/√K ≈ 0.22σ on z, ~0.045 on r,
which prevents tight α-fitting at high W. Re-run W ∈ {30, 210, 2310}
with K=200 baselines (~10× cost, ~3× tighter floor). If r(W=2310)
remains in [0.04, 0.08] this would falsify the closed-form prediction
`r(W=2310) ≈ 0.10` at α=2.07 (would suggest the per-prime decay rate
slows for large primes). Cost: 1 session.

**§D2.a.2.ii — Matched-physical-window protocol.** Hold q_end - q_start
fixed (say 1.5·10⁶ matching W=210 at M=1000) and let M shrink with W.
At W=2310 this gives M ≈ 200 (too few for stable PH); the protocol
breaks down. Alternative: hold log q_n drift ≤ 1 across all W. Re-run
with M_W = round(q_target / (φ(W) · ⟨log q⟩)) for q_target = 10⁶ above
start_x = 10⁶. Confirms the matched-physical-span scan reproduces
matched-M=500 cleanly. Cost: 1 session.

**Successor challenges (proposed in S124):**

**§D2.a.1.i — CLOSED (S131, mode E, B-grade refinement of E2.17).**
B4 = IID with replacement from the W-tricked empirical PMF, run
alongside B1/B2/B3 at d ∈ {2, 3, 4}, M=1000, K=20, three residues.
Three of four pre-stated falsifiers (F_i.2 discreteness direction,
F_i.3 T1 consistency, F_i.4 duplicate-count) PASS at all d. F_i.1
T0 IID-vs-permutation **MARGINAL FAIL** at d=3 (|Δz| = 2.11 vs
threshold 2.0); the failure is fully explained by the duplicate-
point cloud-compression artifact of IID-with-replacement sampling
(0.368M duplicate values per draw create zero-distance pairs in
the Takens cloud, which contribute zero-length H_0 bars but do not
affect H_1 loops — F_i.3 PASS at |Δz|_T1 = 0.08 confirms the H_0-
specific localisation). **NEW T0 component Δ_duplicate ≈ 2-3
mean-gap** (B2 vs B4) isolated; technical book-keeping for any
future PH-on-IID study. Δ_serial_residual after disentanglement is
bounded by ≤ 1σ on T0 (S124's "1-2σ" tightened). E2.17 refined
inline with four-way decomposition.
See `experiments/topological/persistent_homology_w_trick_discrete_b4/persistent_homology_w_trick_discrete_b4_results.md`,
`archive/sessions/session131_d2a1i_ph_discrete_b4.md`.

**Successor challenges (proposed in S131):**

**§D2.a.1.iii — H_0-persistence-without-zero-bars renormalisation.**
Re-compute T0 over only the **non-degenerate** H_0 bars (bars with
finite, nonzero death time). With ≈ 0.368M zero-length bars
filtered out of B4, B4 should match B2 within 0.5σ on the
renormalised T0 — an alternative discreteness probe that bypasses
the duplicate-compression term. Cost: 1 session.

**§D2.a.1.iv — Stratified IID baseline B6.** Sample WITHOUT
replacement M' = 0.632M values uniformly from x, then permute. This
has fewer effective points (matching B4's effective unique count
≈ 632) but no duplicates. The (B6 − B2) gap isolates the *fewer-
points-effect* from the *duplicate-compression-effect*. Cost: 1
session.

**§D2.a.1.ii — Sliding-bandwidth KDE B5(σ).** Vary smoothing
bandwidth σ from 0 (≈ B2 discrete) to a large value (≈ near-Exp(1)
smooth) and trace z(P_W; B5(σ)) on T0 / T1. Predicted: a sigmoidal
crossover at σ ≈ grid-spacing / 2 ≈ 0.16 where the discreteness
component switches off. Cost: 1 session.

**§D2.b — Persistence-image vector classifier.** Replace scalar
T0/T1/L0/L1 with the persistence-image vector (Adams et al. 2017
JMLR 18) and fit a linear classifier on {primes-window, B1-window,
B2-window}. ROC-AUC reports total information content of the PH
descriptor; the classifier's discriminating axis localises which
persistence band carries the signal (interpretable HL fingerprint).
Predicted AUC > 0.95 vs both baselines. Cost: 1 session.

**§D2.c — Sliding-window χ_P indicator embedding.** Embed the
indicator itself rather than gaps: y_n = (χ_P(n), χ_P(n+1), ...,
χ_P(n+d-1)) ∈ {0,1}^d, run PH on Hamming distance. Point cloud
restricted to cube vertices with multiplicities — Gabriel /
witness-complex PH may give a cleaner signature of HL admissible
patterns. Cost: 1 session.

**§D7 — CLOSED (S95, mode I, B-grade).** Determinantal / permanantal
point process fit to χ_P. Tested four progressively flexible kernel
hypotheses at N = 10^7: real DPP (F1: K^2_DPP < 0 at all 15 even t),
real PPP (F2: K^2_PPP < 0 at all 14 odd t > 1), real-signed-K PPP on
all-even sub-lattice (F3+F4: max gap 79.16% on 3-AP triples; required
sign cross-term σ_req ∈ (-0.541, +0.769) never ±1), complex-Hermitian-K
(F5: best LS phase residual 0.0746 vs 0.01 noise floor). All five
falsifiers HOLD. New EDGE **E2.16** — first 3-point structural
confirmation that χ_P deviation = HL singular series (complementing
2-point closures at E2.13 / E2.14 / E2.15). HL prime-by-prime
factorisation cannot reduce to pair-level kernel because triple
admissibility (ν_p saturation) is a multi-body fact. See
`experiments/constructions/primes_dpp_ppp_fit/primes_dpp_ppp_fit_results.md`,
`archive/sessions/session95_d7_dpp_ppp_fit.md`.

**Successor challenges (proposed in S95):**

**§D7.b — Pfaffian point process fit on W-tricked χ_P.** A Pfaffian
PP has matrix-valued kernel `K(t) ∈ R^{2x2}` (block-anti-symmetric)
with k-point correlations as Pfaffians of a 2k×2k matrix. Strictly
extends DPP/PPP. Test on W=210 W-tricked χ_P (where pair correlations
S(0, t)_W ≈ 1, kernel near-zero) and on bare χ_P. Predicted: same
structural failure (admissibility is multi-body, not pairwise) but
the falsification mode is different — produces a structurally
distinct edge candidate. Cost: 1-2 sessions.

**§D7.c — α-determinantal generalisation (Vere-Jones 1997)**: allow
`R_2(t) = ρ^2 + α(t) K(t)^2` with `α(t) ∈ ℝ` offset-dependent. Setting
`α(odd) = +1, α(even) = -1` immediately accommodates the parity sign
flip. The question: do the 3-point α-identities (which generalise
both DPP and PPP) match HL across admissible all-even triples? Build
an α-permanent solver, fit α(t) numerically to minimise 3-point
residual. **A-grade if** consistent fit across all triples (would be
a NEW arithmetic invariant of χ_P). **B-grade if** structural failure
shows the obstruction is prime-vs-pair factorisation independent of
α(t) flexibility. Cost: 1-2 sessions.

**Tertiary backup:** L1 Lean — close `exists_invertible_submatrix` for
the corner case R ≤ 2 (Route A', single session) using `Nat.bertrand`
once at `n = W^(d-1)` plus the existing `chiP_two_eq_one` row-0 entry
to exhibit a 2×2 submatrix with determinant ±1. Note S90 audit:
Route A (full closure) requires Hoheisel-grade prime-density not in
mathlib — multi-session arc on its own.

**§D43 PARTIAL CLOSURE (S157, mode E, B-grade case (i)).** KPZ
universality / Hairer regularity structures on `D(x) := (π(x) − Li(x))
· log(x) / √x` on KPZ-spaced grid `x_k = X/2 + k · ⌊X^{1/3}⌋`. Wavelet
Hölder α(D) ≈ 0.85 stable across logX ∈ {18..24} with linear-fit
r²>0.998 — far above KPZ ceiling 1/2. Pre-stated falsifiers F1–F4
reject KPZ; F2 TW2 skewness rejected at |z|≥5 in 5/7 windows.
Wide-range FFT (x ∈ [10⁴, 10⁷]) confirms D's top peaks at first 12
non-trivial Riemann γ_k (γ_1=14.135 peak/median 770000, γ_2=21.022
ratio 760000) — D is deterministic almost-periodic, NOT stochastic
SPDE. Adds **edge E2.27**. CROSS_DOMAIN_TECHNIQUES §3 Hairer/KPZ/TW2
PROPOSED → USED PARTIAL. See
`experiments/analytic/kpz_pi_li_d43/d43_kpz_pi_li_results.md`,
`archive/sessions/session157_d43_kpz_pi_li.md`.

**Successor challenges (proposed in S157):**

**§D43.b — logX-extension to 2^28.** Scale wavelet Hölder measurement
to logX = 28 (X = 2²⁸ = 268M, ~5x current cost). The CLT-scale theory
of Bombieri/Heath-Brown predicts `Var(D) ~ log log x` (slowly growing)
with Gauss-asymptotic shape; but does Hölder α(D) drift toward 1/2 with
X (asymptotic KPZ-creep — would be a partial positive A-grade signal),
or stay flat at 0.85 (conclusively reject KPZ)? Either outcome is
informative. Cost: 1 session at logX=28; 2 sessions if also testing
logX=26, 27 to extract α(logX) trend slope. Cross-domain ingredient
(unchanged): Hairer regularity / wavelet Hölder.

**§D43.c — CLOSED (S160, mode E, B-grade refinement of E2.27).**
K-truncated residual `R_K(x) := (π(x) − Li(x)) + Σ_{k≤K} 2 Re Li(x^{ρ_k})`
(corrected sign — challenge spec had it reversed; the explicit formula
gives `π − Li ≈ −Σ 2 Re Li(x^{ρ_k}) − log 2` so to cancel we ADD).
Tested at logX = 22, db4 wavelet, K ∈ {0, 1, 5, 10, 25, 50, 100, 200,
500, 1000, 2000, 4000} using the 8000-zero Odlyzko table. Asymptotic
`Ei(z) ~ e^z/z · Σ n!/z^n` is exact on `Re Ei` to rel_err 3·10⁻¹¹
(the iπ Stokes term is purely imaginary; cancels in 2 Re).

**Outcome: A-grade conjecture REFUTED.** Full-band α(D_K) drops below
1/2 at K* ∈ [200, 500] but with crashing r² (0.998 → 0.005); fine-band
fit (only wavelet levels above γ_K cutoff) gives U-shape α: smooth at
K = 0..25, dipping to −0.40 at K=1000 (r²=0.48), rebounding to +0.74
at K=4000 (r²=0.70). No (K, fit) cell achieves both α<0.5 AND r²>0.5
— no clean KPZ regime.

**B-grade content (new):** (i) variance reduction curve
var(R_K)/var(R_0) = 1.00, 0.85, 0.68, 0.54, 0.48, 0.35, 0.27, 0.23, 0.21,
0.19, 0.18, 0.18 at K = 0,1,5,10,25,50,100,200,500,1000,2000,4000 —
first 50 zeros account for ~65% of var(π−Li); first 4000 for ~82%;
(ii) Cramér control with same subtraction is structurally inert:
var(R_K^C)/var(R_0^C) = 1.012–1.022 flat, α_full(C) ≈ 0.98 at all K,
α_fine(C) ≥ 0.89 at all K — the explicit formula has nothing to remove
from Cramér; (iii) **π-vs-Cramér α_fine gap ≤ −0.2 from K=50 onwards**
is a new structural measurement: the explicit formula structurally
describes π in a way invisible on a Bernoulli prime model.
**Refines E2.27 inline.** See
`experiments/analytic/d43c_k_truncated_residual/d43c_k_truncated_residual_results.md`,
`archive/sessions/session160_d43c_k_truncated_residual.md`.

**Successor challenges (proposed in S160):**

**§D43.c.i — Higher logX (logX=26) confirmation of variance curve.**
At logX=22 the variance ratio plateaus at ~0.18 by K=2000. Does the
plateau drop at higher logX (more variance to extract from more
zeros)? Predicted: var(R_K)/var(R_0) at fixed K decreases as logX
grows. If the plateau drops below 0.10 at logX=26, the structural-
description claim sharpens. Cost: 1 session (sieve at 2^26 ≈ 67M is
~30s; full sweep ~5 min).

**§D43.c.ii — Asymptotic 1/2 ceiling test.** The infinite-tail residual
has Hölder regularity exactly α* = 1/2 − ε. To detect it empirically,
need K ~ √x to reach amplitudes ~ x^{1/2}/γ_K large enough to escape
the Cramér noise floor. This requires K ≫ 10⁵, beyond the available
8000-zero Odlyzko table. Path: rigorously compute high zeros via
mpmath.zetazero up to γ ~ 10⁵ (cost: hours per zero) AND scale logX
to ~30. Multi-session arc (probably out of practical reach).

§C1 (BK Odlyzko zeros) is now **closed** by S71-redux with a quantitative
L⁴ obstruction (N_required ≥ 0.81 L⁴ for κ=3 detection); no further
attempts there.

---

## §1. Composition Challenges

A composition challenge says: build a single mathematical object whose
*definition* uses two or more existing edges from EDGES.md, and whose
*evaluator* produces a number you can compare to a baseline. Why this
matters: edges are catalogued individually; their *interactions* are
mostly unexplored. Composition is the cheapest source of genuinely new
mathematical objects.

**Discipline:** every composition target must produce code that runs.
File under `experiments/constructions/<name>/`. Cite the edges by ID
in `definition.md`. Failure mode is "duplicate of [composed-edge-X]"
or "construction-incoherent" if the object isn't well-defined.

### C1 — A⊕C₃ bisection × 0.537-bits invariant — **BUILT (S70)**
**Edges:** E1.6 + E1.5
**Object:** A function `g_q : ℕ → 𝔽_q × 𝔽_q` with `g_q(x) = (A(x) mod q, C_3(x) mod q)` (the paired bisection; q=2 recovers E1.6's two parities). Test whether the per-step entropy closed form (S69-sharpened E1.5) holds for `g_q` for q ∈ {3, 5, 7, 11, 13}. If yes, a *family* of E1.5-style invariants where the original was a single statement.
**Why novel:** E1.5/S69 was for π(x) only. Testing on A and C_3 (and their joint state) is new.
**Falsification:** PR1 = |H_emp(F mod q | prev) − h_2(F(X)/X)| < 5e-3 for F ∈ {π, A, C_3}; PR2 = joint H(g_q(x) | prev) − H_3(1−ρ_A, ρ_π, ρ_C3) < 5e-3; PR3 = I(A mod q ; C_3 mod q) < 0.01 bits at X = 2·10⁶.
**Outcome (S70):** All three PASS at full N = 2·10⁶. Worst PR1 |diff| 1.6e-3 (at small X), 5e-6 at X = 2·10⁶; worst PR2 |diff| 8e-4 in strict regime q²·100 ≤ π(X); worst PR3 I = 1.17e-4 bits at q = 13. Three new artifacts: (i) per-component closed form extended from π to A, C_3 (mechanism applies to all monotone Ω-stratified summatories); (ii) **new joint closed form** `H(g_q | prev) = H_3(1−ρ_A, ρ_π, ρ_C3)`; (iii) q-stable strengthening of E1.6 marginal independence to q ≤ 13. Composition reveals destructive-interference picture: A, C_3 each carry 1 bit/step asymptotically, π carries 0 bits/step, marginally q-stably independent.
**Status:** BUILT, no polylog opening. Closure row in CLOSED_PATHS.md (S70). E1.5 and E1.6 annotated in EDGES.md.
**Save under:** `experiments/constructions/g_q_bisection_invariant/`

### C2 — MPS bond-dim × free-probability moments — **BUILT (S74)**
**Edges:** E2.1 + (free probability framework, partially explored S58)
**Object:** Treat χ_P as a non-commutative random variable in the algebra generated by base-W shift operators. Compute its first 4 free cumulants via the MPS reshape. The MPS bond-dim formula `min(W^j, φ(W)·W^{d-j-1}+1)` should constrain which free moment patterns are achievable. Build the cumulant evaluator and check whether the resulting sequence `(κ_1, κ_2, κ_3, κ_4)` matches a known free distribution (semicircle, free Poisson, Marchenko-Pastur).
**Why novel:** No one has computed free cumulants of χ_P in this framework. If they hit a known distribution, you've found a *non-commutative analogue* of the wheel-W theorem. If they don't, you've added a 36th pseudorandomness measure.
**Falsification:** Free cumulants land within or outside known distributions.
**Outcome (S74):** χ_P's MPS unfolding spectrum has a **three-part structure**: (a) a finite "structural" peak (rank-1 mean + φ(W)−1 residue indicators) reproduced by the matched-active-Bernoulli baseline; (b) a **spike band of O(N^{0.42}) outlier eigenvalues** present in χ_P but absent from the baseline; (c) an MP bulk that **matches Marchenko-Pastur with rate `c = φ(W)/W = ∏_{p≤W}(1 − 1/p)` (the Mertens product)** within 5–10% relative on κ_2, κ_3, κ_4. Cross-domain import: free probability (Mingo-Speicher 2017). New empirical regularity: outlier count `k*(W, d) ≈ R^{0.85}` where `R` is the rank — fitting `α ≈ 0.85` over W=2 d=14..22. **Algorithmic implication:** any spectral compression of χ_P faithful at the second-moment level needs rank Ω(N^{0.42}) — recovers a polynomial-in-N barrier from a free-probabilistic angle, parallel to (not stronger than) Lagarias-Odlyzko. No polylog opening.
**Status:** BUILT. Three new artifacts: (i) **MP-bulk identification** — the wheel-W sieve density `φ(W)/W` IS the free-Poisson rate of χ_P's MPS bulk; (ii) **spike-count regularity** — `k*(W, d) ∝ R^{0.85}`, a new structural empirical fact beyond E2.1's rank statement; (iii) verified MP(c) standardized free cumulant identity `κ_r = c^{r-1}` and cross-domain free-cumulant evaluator.
**Save under:** `experiments/constructions/free_cumulants_chi_p/`

### C3 — Brandt obstructions × per-bit difficulty — **BUILT (S105)**
**Edges:** E5.8 + E1.3
**Object:** E5.8's obstruction (O1) says Brandt's hard string is oracle-dependent, not fixed. E1.3 says π(x) has a sharp 4-bit difficulty boundary. Compose: define `π_J(x)` = "the J-th bit of π(x) for J = 0.7N" (within E1.3's hard zone). Ask whether **`π_J` is a fixed natural function** for which a Brandt-style oracle-aware argument might still apply. The point is to find the *minimal weakening* of E5.8's O1 that still admits Brandt's TRAVERSE technique.
**Why novel:** E5.8 closed Brandt for π(x) mod 2 wholesale. A per-bit version might side-step O1 if J is parameter-controlled.
**Falsification:** Either O2/O3/O4 also obstruct π_J (closes per-bit too), or one of them survives (genuinely new partial result).
**Outcome (S105):** **F2 holds in the saturation regime (N ≥ 4); empirical sharpening of F2.** Built per-bit family `{s_J^(N) : J = 0..N-1}` for N ∈ {3..7} and computed `Kt_b(s_J^(N))` under brandt_mktp's 3-bit stack VM (L_MAX=12, INF=61). Bounded-Kt cut location: `J*(N) := ⌈log₂(π(2^N - 1) + 1)⌉` — for J < J*(N), all bits saturate (Kt_b = 61), including E1.3's "easy zone" J ∈ [⌈0.5N⌉, J*); for J ≥ J*(N), `π_J ≡ 0` and Kt_b ∈ {5..8}. The bounded VM is **blind to E1.3's smooth/oscillatory transition** — it only sees the trivially-zero boundary, materially higher than `0.5N` (gap ≈ 0.5N − log₂ N). Structurally, all four Brandt obstructions O1-O4 still apply to every fixed `π_J`: O1 (each π_J is a fixed function — per-bit decomposition is parameter-controlled but each π_J itself is a fixed total Boolean function); O2 (π_J answers "is bit J set", not "is z complicated"); O3 (per-bit family supplies J*(N) ≈ N − log₂ N strings, but no traversal-path-dependent fresh prefixes); O4 (uniform DTIME, not circuits). The "minimal weakening of O1 admitting TRAVERSE" the spec asked about does **not** exist on `{π_J}`. **Status:** BUILT, no polylog opening. Closure mode E (DUPLICATE-PLUS of E5.8) at structural level + empirical refinement of E1.3 at bounded-Kt resolution.
**Successor C3.a (S105 → BUILT S150):** arithmetic-primitive bounded-Kt VM. Built 4-bit-per-op extended VM with all four R⁻¹-kernel primitives {LOG2, LI_APPROX, DIV_LOG, GEO_SUM} added to base 8-op stack VM; T_MAX = 4096; C inner-loop simulator (`sim.c → sim.so`) with batch evaluator. Scanned (N=3..6, J=0..N-1) at L_max ∈ {24, 28}. **Verdict: F3 (intermediate hierarchy)**. The cut shifts to E1.3's `⌈N/2⌉` boundary fully for N ∈ {4, 5} at L_max=24 (matched to target_lens 16, 32) — `bit_2(π)` at N=4 and `bit_3(π)` at N=5 BOTH compress, with explicit programs `ADD, LI, LI, EMIT_S, PUSH_N, PUSH_N` (24 bits, computes `bit_0(LI(LI(2x)))` matching `bit_2(π(x))` exactly on [0,16)) and `EMIT_S, LI, DUP, INC, PUSH_N, PUSH0` (24 bits). At N=6 the easy zone splits at L_max=28: `J=4` (closer to `J*(N)=5`) compresses via triple-LI program `EMIT_S, PUSH_N, LI, LI, LI, DUP, EMIT_S`; `J=3` (closer to `⌈N/2⌉`) remains saturated. **Within-easy-zone J-monotone hierarchy identified**: bits closer to `J*(N)` compress at smaller L_max — bounded-Kt cost of `bit_J(π)` is empirically monotone in N − J. Hard-zone bits saturate at L_max ≤ target_len for all N ≥ 4. Iterated LI applications (LI∘LI, LI∘LI∘LI) are the dominant compression mechanism. **Refines E1.3** with explicit VM-richness × N-dependent cut hierarchy (annotated EDGES.md). **E5.8 unchanged** — empirical compressibility of some `s_J^(N)` does not enable a Brandt-style diagonalisation; obstructions O1-O4 are independent of VM choice. CLOSED_PATHS row added (S150). See `experiments/constructions/brandt_per_bit_arith_vm/`.

**Successor challenges (proposed in S150):**

**C3.a.i — L_max = 32 sweep.** At L_max=32 (4.3B programs, ~5h with current C simulator; doable with multi-core), test whether `(N=6, J=3)` finally compresses. If yes, the F1 cut would extend to N=6 and the hypothesis "L_max needs only ~ target_len bits to reach E1.3's `⌈N/2⌉` boundary" gets one more confirmation. Cost: 1 session with parallelism.

**C3.a.ii — Random-program sampling at larger N.** Run N=7 (target_len=128) with random sampling of K=10⁹ random 128-bit programs (full enumeration intractable). Test whether the within-easy-zone J-monotone hierarchy persists at larger N. Cost: 1 session.

**C3.a.iii — VM-budget vs RH-scale.** S146 found the bit-level RH-scale anti-correlation valley sits at `J* = ⌊log₂(p(N))/2⌋` (Skewes-direction error scale). Test whether the VM-budget threshold `L_max(J)` for compressing bit J of π scales with the RH error magnitude. If `L_max(J) ∝ √π(2^N) · log` near the RH-scale J, this would be the first project link between **bounded-Kt complexity** and **RH-shadow bit phenomena** (S146). Cost: 1-2 sessions.

**C3.a.iv — Drop / swap arithmetic primitives — BUILT (S158, mode E, B-grade refinement of E1.3).** Six ablation conditions on `{LOG2, LI, DIV_LOG, GEO_SUM}` at `L_max=28`, `N ∈ {3, 4, 5, 6}`: baseline + four single-drops + only_LI. **Verdict: F4 — no single arithmetic primitive is strictly necessary for the easy-zone cut shift.** Every easy-zone cell `{(3,2), (4,2), (5,3), (6,4)}` that compresses in baseline also compresses under every single-drop AND under only_LI. drop_LI matches baseline at L≤24 cells with alternative programs (`(N=4,J=2)`: GEO_SUM+DIV_LOG; `(N=5,J=3)`: LOG2 alone), and is +1 bit at `(N=6, J=4)`. **Refutes S150's narrowest reading** that LI is the dominant compression mechanism — it is one of four substitutable mechanisms. Hard-zone primitive sensitivity differs: meaningful cell `(N=5, J=2)` requires LI ∧ DIV_LOG (either-alone saturates), an orthogonal observation. **Refines E1.3** with primitive-class-robustness: the cut shift is driven by the FAMILY of slow-growing integer-function primitives, not by LI specifically. CLOSED_PATHS row (S158); EDGES.md annotated. See `experiments/constructions/brandt_per_bit_arith_vm_ablation/`, `archive/sessions/session158_c3aiv_arith_vm_ablation.md`.

**Successor challenges (proposed in S158):**

**C3.a.iv.α — LI-vs-non-LI gap scaling at large N.** The +1-bit
`(N=6, J=4)` `drop_LI` gap is the only cell where LI is strictly
shorter than non-LI alternatives at `L_max=28`. At larger N the
iterated-LI mechanism may pull strictly ahead of GEO_SUM/LOG2-based
alternatives (the slow-growth match to `Li(x)` is asymptotically
sharper than to alternative slow-growers). Run N=7 with random
sampling at `L_max=32` under (baseline) and (drop_LI) to estimate
gap(N) := Kt_b'(drop_LI) − Kt_b'(baseline) at the easy-zone cell
closest to `J*(N)`. Cost: 1 session with multi-core sampling.

**C3.a.iv.β — Successor primitive set: LN, INV_LI, ZETA_K.**
Replace one of the four S150 primitives with an alternative
slow-growing primitive (LN truncating-natural-log; INV_LI = R⁻¹;
ZETA_K computing `floor(zeta(2 + k))` for small k). Re-run S158's F4
test with the new primitive. Question: is the F4 family-robustness
specific to {LOG2, LI, DIV_LOG, GEO_SUM} or does it survive any
4-primitive set whose elements grow at sub-linear rates? Cost: 1-2
sessions.

**Save under:** `experiments/constructions/brandt_per_bit/` (C3); `experiments/constructions/brandt_per_bit_arith_vm/` (C3.a); `experiments/constructions/brandt_per_bit_arith_vm_ablation/` (C3.a.iv).

### C4 — Aggarwal binary search × Dusart bracket × BPSW oracle — **BUILT (S120)**
**Edges:** E6.6 + E6.8 + E5.1
**Object:** Build a *unified* p(n) library that combines all three: Aggarwal's `O(sqrt(n) log^4 n)` binary-search reduction, Dusart's narrow bracket of width n, and BPSW as a TC⁰ primality oracle. The composition is non-trivial because BPSW is conditional. Build it, benchmark vs `algorithms/v10_c_accelerated.py`. The novel content is whether the conditional BPSW step propagates correctly through Aggarwal's wrapper.
**Why novel:** No published primecount-style library combines all three. Even the integration is publishable.
**Falsification:** Wall-clock benchmarks at p(10^k) for k = 6, 9, 12, 15.
**Outcome (S120):** Three modes built (`agg`, `bpsw`, `hybrid`); F1 (sympy.prime agreement) HOLDS at n ∈ {10⁴..10⁷}; F3 (hybrid 10×-faster than bpsw) HOLDS at 21×/34×/53× over n ∈ {10⁴..10⁶}; F2 (hybrid 1.5×-faster than agg) HOLDS only with C-accelerated pi (1.64× at n = 10⁷); fails in pure Python (1.30×); F4 (U-shape K-curve) HOLDS only with fast pi oracle: optimal `K*` shifts as `K* ~ pi_cost / bpsw_cost`. **Net new content**: (a) **K* depends on the pi/bpsw cost ratio** — Python `pi_lucy`: K* = width (composition collapses to E6.8 + E5.1, no Aggarwal narrowing); C-Lucy: K* ≈ 16K ≈ √width at n = 10⁷; HKM projection: K* small constant. This cost-ratio knob is invisible in Aggarwal 2025's asymptotic analysis. (b) **BPSW conditionality propagates 1-to-1 through Aggarwal's wrapper**, not amplified — a single BPSW pseudoprime in the final bracket shifts the answer by exactly one prime; the wrapper does not compound conditionals. (c) **Dusart bracket alone is worth ~50× over naive BPSW-from-2** (predicted `2 log p_n`, observed 21-53×). **Closure mode E** (constant-factor re-arrangement; no asymptotic improvement). The asymptotic falsifier on p(10^k) for k ≥ 9 is bottlenecked by Lucy DP / HKM in pi(x), not by the composition itself. Cites E5.1, E6.6, E6.8 (annotated S120). See `algorithms/aggarwal_dusart_bpsw/`.
**Status:** BUILT, no polylog opening. Closure row in CLOSED_PATHS.md (S120). E5.1, E6.6, E6.8 annotated in EDGES.md.
**Save under:** `algorithms/aggarwal_dusart_bpsw/` (working algorithm goes in `algorithms/`, not `experiments/constructions/`).

### C4.a — Replace `pi_lucy` with HKM/primecount and re-locate K* (NEW, follow-up to S120)
**Why:** S120's K-sweep showed U-shape only with C-accelerated Lucy DP (K* ≈ √width at n = 10⁷). The natural extension: replace the Lucy DP oracle with primecount or a tractable HKM port (the asymptotically O~(√x) algorithm from Hirsch-Kessler-Mendlovic 2024) and trace the optimal K. **Predicted outcome**: K* drops further toward small constant ≪ √width, and the hybrid scheme's gain over `agg` shrinks toward zero in the asymptotic limit (Aggarwal-pure dominates). If TRUE, the C4 composition is a *finite-x*-only improvement; if K* stays at √width even with HKM, BPSW retains constant-factor value asymptotically — a cleaner phrasing of "BPSW saves the trailing pi-call regime even at the asymptotic frontier."
**Save under:** `algorithms/aggarwal_dusart_bpsw/hkm_extension/`. Cost: 1-2 sessions.

### C5 — N/2 universality × non-Boolean function — **BUILT (S71)**
**Edges:** E1.4 + E2.5
**Object:** The N/2 universality (E1.4) is verified for π(x) mod 2 across 6 measures. Pick another natural number-theoretic Boolean function — e.g., `is_squarefree(n)`, `is_squarefree_AND_n_mod_3==1`, `mu(n) == +1`, `Liouville(n) == +1`. Run the same 6 N/2 measurements. If N/2 universality holds across multiple functions, you've found a *meta-theorem* about a class of NT functions; if it doesn't, you've found a function where π(x) mod 2 is special (which is also genuinely interesting).
**Why novel:** The N/2 universality has only been measured on one function.
**Falsification:** Per-function table of the 6 measurements.
**Outcome (S71):** Ran cheap 4-measure subset (M1 comm-rank, M2 BM-LFSR, M3 approx-degree, M4 PTF) at N up to 14 on six functions {f_pi, f_sqfree, f_mu_pos, f_lam_pos, f_sqfree3, density-matched PRF}. **N/2 universality is NOT universal**: it holds tightly only for the parity-of-Omega family `{chi_P, mu_pos, lam_pos}` at adeg/PTF, but BREAKS for sqfree (below), sqfree3 (above), and PRF (below). M1 produces three distinct **closed-form rank regimes**: rank(M_chi_P) = 2^{N/2-1}+1 exactly; rank(M_sqfree) = rank(M_mu_pos) = 3·2^{N/2-1} exactly; rank(M_lam_pos) = 2^{N/2} (full). All three explained by the same structural identity: rank(M_f) ≤ (1−ρ_f)·2^{N/2} where ρ_f is the density of lower-half columns forced to zero by f's smallest non-witness modulus (chi_P: ρ=1/2 from "x even ⇒ x non-prime"; sqfree, mu_pos: ρ=1/4 from "4 | x ⇒ x non-squarefree"; lam_pos: ρ=0). M2 BM/2^N ≈ 0.50 universally — the LFSR test alone is non-distinguishing. **Status:** BUILT, no polylog opening. This **refines E1.4** (scope = parity-of-Omega family) and **unifies E2.7 + E2.8** by a single column-zero density principle.
**Status:** BUILT, no polylog opening. Closure row in CLOSED_PATHS.md (S71). E1.4 and E2.7 annotated in EDGES.md.
**Save under:** `experiments/constructions/n_over_2_universality_class/`

### C6 — Three-pillars × HKM time-space curve — **BUILT (S81)**
**Edges:** E7.7 + E6.7
**Object:** E7.7 says only 3 informationally-complete encodings of π(x) exist (prime positions, zeta zeros, floor values). E6.7 puts HKM at the Pareto-frontier `(time = N^{8/15}, space = N^{1/3})`. Build a "three-pillar tradeoff diagram" — for each pillar, compute the achievable time-space curve, and ask whether HKM is on the *zeta-zero* pillar's frontier or the *floor-values* pillar's frontier. The novel content is the structural placement.
**Why novel:** No one has classified the algorithmic tradeoffs by which pillar they live on.
**Falsification:** Tradeoff diagrams for at least 2 pillars.
**Outcome (S81):** Built (alpha, beta) catalog of 14 algorithms (4 prime / 4 zeta / 6 floor). All four pre-stated falsifiers (F1: HKM dominated; F2: HKM cross-pillar dominated; F3: identical pillar frontiers; F4: floor T*S < 5/6) PASS. **Three structural findings:** (i) HKM is on the **floor pillar** Pareto frontier and dominates every other floor-pillar entry elementwise; (ii) **HKM uniqueness**: no zero-pillar or prime-pillar algorithm achieves both `T ≤ N^{8/15}` AND `S ≤ N^{1/3}` simultaneously — HKM's `(8/15, 1/3)` point is uniquely accessible on the floor pillar; (iii) **non-overlapping pillar regions**: time-only minimum is shared by prime+zeta pillars at α=1/2, space-only minimum is unique to floor at β=1/3, T*S minimum is unique to floor at 13/15 ≈ 0.867 (saturating E7.6 to within N^{0.034}). Aggarwal (E6.6) is observed as a meta-algorithm whose effective placement migrates with the chosen pi(x). No polylog opening — refinement of E7.7 + E6.7. CLOSED_PATHS row added (S81). E6.7 and E7.7 annotated.
**Status:** BUILT, no polylog opening.
**Save under:** `experiments/constructions/pillar_tradeoff_diagram/`

### C7 — Calibrated 1-bit-bias random control for the S84 PRIMES-vs-random depth-2 sign-threshold gap — **BUILT (S89)**
**Edges:** E1.10 / E3.13 + S84 result.
**Object:** Class-conditional matched-density random Boolean function `f_cal` on {0..63}: P(f_cal(x)=1 | x odd) = 17/32, P(f_cal(x)=1 | x even) = 1/32 (matching PRIMES at N=6). Two modes: STRATIFIED (exact 17 odd + 1 even, weight always 18, bit_0_acc always 0.75) and BERNOULLI (independent draws). For each, run S84's `enum_d2_smart` ILP harness over K=1458 W=1 candidates and find min depth-2 sign-threshold size M.
**Why novel:** No previous project session has measured the depth-2 sign-threshold size of any class-conditional matched random function (calibrated or otherwise). The S84 result was a single-mechanism conjecture; C7 tests it directly.
**Falsification (pre-stated):** F1 ≥1 stratified > 6 → residual; F2 all stratified ≤ 6 → gap is oddness; F3 ≥1 stratified < 6 → calibration overshoots PRIMES; F4 stratified mean < unbiased mean → bit_0 primary driver.
**Outcome (S89):** **F2 + F3 + F4 hold; F1 fails.** Stratified (n=20, weight=18 always, bit_0_acc=0.75 always): histogram = {5: 4, 6: 16}, mean = 5.80, max = 6. **0/20 stratified samples need M > 6** — full S84 deviation absorbed. Bernoulli (n=20, variable weight, bit_0_acc 0.625-0.81): histogram = {5: 7, 6: 11, 7: 2}; both M=7 cases have bit_0_acc < 0.75, empirically confirming the proposed monotonic relationship. Calibrated mean (5.80) is BELOW PRIMES (6.00); 4/20 stratified samples are *strictly easier* than PRIMES (M=5). PRIMES sits at the +0.5σ position of the calibrated distribution; under the calibrated null `P(M ≤ 6) = 1.0` (vs the unbiased null where `P(M ≤ 6) = 0/10`). **The S84 gap REDUCES to elementary parity ("π(x) ≈ 1 iff x odd, for x > 2"), no PRIMES-specific structure beyond oddness.** Recommends a footnote on `novel/pseudorandomness_of_pi.md`: the 36th measure deviates from unbiased random but the deviation is fully captured by the trivial 1-bit oddness fact; pseudorandomness thesis stands "modulo the obvious mod-2 bias."
**Status:** BUILT, no polylog opening. Closure row in CLOSED_PATHS.md (S89). Refines E1.10 / E3.13 / E1.6 reading.
**Save under:** `experiments/circuit_complexity/sat_tc0_primes_n8_calibrated/`

### C7.a — Calibrated control at N=8 (NEW from S89, follow-up)
**Why:** S89 closed C7 at N=6. The natural extension: redo the calibrated baseline at N=8 (where 53 of 54 primes are odd, P(f=1|odd) = 53/128, P(f=1|even) = 1/128). The unbiased N=8 search in S84 didn't terminate at M ≤ 16 (k_max=4); a calibrated N=8 study would be marginal at the same time-budget.
**Predicted outcome:** Same — calibration absorbs the gap; PRIMES sits in the calibrated distribution. If FALSIFIED at N=8 (calibrated needs M much larger than PRIMES), the bit_0 explanation breaks at higher N and a residual mechanism enters — A-grade material.
**Save under:** `experiments/circuit_complexity/sat_tc0_primes_n8_calibrated/n8_extension/`

### C8 — Depth-2 sign-threshold weight-vs-size tradeoff for PRIMES — **BUILT (S127)**
**Edges:** E5.3 + S84 framework + E1.6 / C7-S89.
**Pivot from spec:** the N=8 ILP did not terminate within session budget at any cell beyond M ≥ 4 even at W=1 (S84 already reported this). N=6 was used as the headline scale; the N=4 cross-N analog was added as validation.
**Outcome (S127):** First measured weight-vs-size tradeoff curve for any natural NT Boolean function in depth-2 sign-threshold class.
- **N=6 curve: M\*(W) = (6, 4, 3, 3, 3, 3, 3, 3) at W ∈ {1, 2, 3, 4, 8, 16, 32, 64}**, with `M*(W=1) = 6` from S84 column enumeration, all other cells from this session. M=2 UNSAT proven up to W=64 (4-5 s each); M=3 SAT proven at W ∈ {3, 4, 8, 16, 32, 64} (6-65 s, witnesses verified 64/64); M=3 UNSAT proven at W=2 (277 s).
- **Structural floor at M\*=3 reached at W=3 and held through W=64** — 16× weight increase yields zero gate reduction.
- **N=4 cross-N analog:** PRIMES `M* ∈ {3, 2, 2, 2}`, random `M* ∈ {4, 3, 3, 3}` at `W ∈ {1, 2, 4, 8}`; F4 (PRIMES easier than random) holds with Δ=1 gate at every W.
- **N=6 random control (seed 1):** even at 600 s budget per cell the `(W=4, M=3)` random cell is UNKNOWN, vs PRIMES (W=4, M=3) SAT in 113 s. Time-asymmetry is itself a structural gap (in search-problem complexity, not necessarily circuit complexity).

**Falsifiers verdict:** F1 (flat plateau) FAILS as predicted — `M*(N=6)` decreases. F2 (geometric decay) HOLDS at W=1→2 and W=2→4 then FAILS (saturates) from W=4. F3 (M=1 collapse at finite W) FAILS up to W=64. F4 (PRIMES easier than random) HOLDS at N=4 with Δ=1 gate at every W; partial evidence at N=6 via time-asymmetry. F5 (N=4 closed-form `M*(W ≥ 3)=1`) FAILS — `M*(N=4)=2` even at W=8.

**Status:** BUILT, no polylog opening. **Refines E5.3** with quantitative `(W, M*)` curve at N=6, the first such measurement for any natural NT function. **Refines C7-S89** PRIMES-vs-random gap from a single-W observation to a curve. CLOSED_PATHS row added (S127). E5.3 annotated.
**Save under:** `experiments/constructions/d2_sign_threshold_w_m_tradeoff/`

### C8.a — N=8 partial scan (NEW from S127, follow-up)
**Why:** The C8 spec originally targeted N=8 but the N=8 ILP scaling pushed the work to N=6. The N=8 column at non-trivial W is genuinely unmeasured. S84 reports the W=1 column is `M ≥ 17` (column-enum partial bound, k_max=5).
**Predicted shape:** by the N=4 → N=6 progression (`M*(W=∞) ≈ 2 → 3`), expect `M*(N=8, W=∞) ∈ {4, 5, 6}` — a structural floor in single-digits even at unbounded weight. The transition `W` at which the floor is reached should also grow with N.
**Save under:** `experiments/constructions/d2_sign_threshold_w_m_tradeoff/n8_extension/`. Cost: 1-2 sessions (each cell allow ≥ 30-min CBC, parallelise across W).

### C8.b — Random-control F4 resolution at N=6 — **BUILT (S135)**
**Edges composed:** E5.3 + S84 framework + E1.6 / C7-S89 + C8/S127.
**Object:** Extended column enumeration `Θ(N, W) = {distinct sign-threshold truth tables on N inputs with weights in {-W..W}}` to W ≥ 2 (catalog sizes K(W=1) = 1458, K(W=2) = 30,898, K(W=3) = 218,066). Run S84's `depth2_search` on PRIMES vs density-matched random N=6 across multiple seeds.
**Why novel:** S127's joint (alpha-bilinear) ILP returned UNKNOWN on every random `(W ≥ 4, M=3)` cell within 600 s. The cross-encoding shift (column-enum) eliminates the bilinear constraints and resolves the M=3 cells in 147–230 s.
**Falsification:** F0/F0' sanity (PRIMES reproduces S84/S127); F1 (random easier at W=2 → REJECTED); F2 (Δ=0 tie at W=2); F3 (Δ ≥ 1 strict, requires M=4 UNSAT proof); F4 (cross-seed robust).
**Outcome (S135):** **F0, F0' HOLD; F1 REJECTED on every seed; F2 partially established; F3 UNRESOLVED (random M=4 UNKNOWN at 600 s); F4 HOLDS at three seeds {1, 5, 42}.** Random N=6 W=2 M=2,3 are both UNSAT for every tested seed (M=2 in 130–196 s, M=3 in 147–230 s); PRIMES W=2 M=3 UNSAT in 157 s, M=4 SAT in 181 s. **`M*(rand_s; W=2) ≥ 4 = M*(PRIMES; W=2)` for s ∈ {1, 5, 42}**, robust direction-confirmation of F4. Magnitude (Δ = 0 vs Δ ≥ 1) remains open at W=2 because random M=4 cells exceed CBC's 600 s budget at K=30898. Cross-encoding finding for the S84 framework: column-enum proves W=2 M=3 UNSAT 1.8× faster than joint-ILP on PRIMES; intractable joint-ILP M=3 random cells become tractable in column-enum.
**Status:** BUILT, no polylog opening. **Refines E5.3** (PRIMES `M*(W=2) = 4` not breakable by density-matched random) and **C7-S89/E1.6** (W=2 PRIMES-vs-random gap holds in direction even outside calibrated-oddness regime).
**Save under:** `experiments/constructions/d2_sign_threshold_w_m_tradeoff/random_n6_resolution/`

**Successor challenges (proposed in S135):**

**C8.b.i — Random M=4 W=2 SAT search via greedy/LNS.** CBC ILP UNSAT proof on (random, W=2, M=4) is intractable; a greedy triple-enum or local search may find a SAT witness if one exists. If no SAT in ~1 h CPU, empirical evidence for `M*(rand; W=2) ≥ 5` (Δ ≥ 1 strict). Cost: 1 session. Save under `random_n6_resolution/m4_search/`.

**C8.b.ii — Tighter F4 via W=3 column-enum on random.** K=218,066 at W=3 borderline. PRIMES W=3 M=3 SAT (S127, 65 s); does random M=3 stay UNSAT at W=3? If yes across multiple seeds, F4 strengthens by one weight-doubling step. Cost: 1-2 sessions (each cell ~600-1200 s).

**C8.b.iii — Seed-distribution of random M\*.** Run N=6 W=2 M=3 across 100+ seeds in parallel. Empirical histogram of M\*(rand) tests "random is always ≥ 4" vs "random is sometimes 3 on a non-trivial seed fraction". With 32-core parallelism + ~150 s/cell, ≈ 10 min wall-clock. Cost: 1 session.

### C9 — Pointwise Ramanujan-spike approximator T_Q(n) — **BUILT (S191, paradigm-shift mode)**
**Edges:** E2.1 + E1.5 + E1.6 + E2.2 + S168.
**Object:** `T_Q(n) := (π(N)/N) · Σ_{q sqf ≤ Q} mu(gcd(q,n)) / phi(q/gcd(q,n))` — pointwise dual of S168's energy-level spike formula, computable in `O(Q·ω(n))` per evaluation. Falsifiers: (PR1) `||T_Q − const||² / π(N)` at Q=N^{0.185}, d=20 in [0.18, 0.26] (matches S169 SVD spike block). (PR2) Pearson(chi_P, T_Q) monotone-increasing in Q. (PR3) Precision-at-π(N) > Q · π(N)/N. (PR4) Hölder identity `mu(q) c_q(n)/phi(q) = mu(gcd)/phi(q/gcd)` exact for sqf q ≤ 30, n ≤ 60.
**Outcome (S191):** All four falsifiers PASS. At d=20, Q=13: ratio = 0.2229 (matches S169's 0.220 to 1.4%). Precision lift 5.6× at Q=N^{0.185}, scaling to 12.8× at Q=√N. **Wheel-Mertens identity**: at primorial Q=W and n coprime to W, `T_W(n) = (π(N)/N) · W/φ(W)`. **Status:** BUILT. No polylog opening (cost is `Q · ω(n) = Ω(N^{0.185})`); structural picture only. EDGES.md E2.1 annotated, CLOSED_PATHS row added. **Save under:** `experiments/constructions/ramanujan_spike_pointwise/`.

**Successor challenges (proposed in S191):**

**C9.a — Divisor-only restriction at primorial Q — BUILT (S208, mode E, B-grade refinement of E2.1 / S191).**
**Edges:** E2.1 + E2.2 + S191 / C9.
**Object:** `M_W^{div}(n) := Σ_{q | W, q sqf} mu(gcd(q,n)) / phi(q/gcd(q,n))`
and `T_W^{div}(n) := (pi(N)/N) · M_W^{div}(n)`; the divisor-only
truncation of S191's full pointwise spike `M_Q`.
**Closed-form identity (proved analytically and verified pointwise to
machine epsilon):** for every squarefree `W ≥ 1`,
```
   M_W^{div}(n)  =  [gcd(n, W) = 1] · W / phi(W).
```
For `W` not squarefree, the same identity holds with `W` replaced by
`rad(W)`. The proof is a one-line Möbius collapse:
`Σ_{A ⊆ P_n}(-1)^|A| = [P_n = ∅]` factorising the summation across the
squarefree-divisor lattice of `W`.
**Outcome (S208):** All six pre-stated falsifiers PASS. F1 (5
primorials), F2 (squarefree non-primorials `15, 105`), F3 (non-squarefree
`12, 60, 420` reducing to radical) all PASS pointwise to abs_err
≤ 8.88·10⁻¹⁶ on `N = 10⁶`. F4 (mean = 1) PASS to ≤ 3·10⁻⁶. F5 (L²
closed form `||T_W^{div} − chi_P||²/N` matches explicit prime-distribution
count) PASS to rel_err < 10⁻⁹ on every cell. F6 (Pearson gap
`Pearson(chi_P, T_W) > Pearson(chi_P, T_W^{div})`) PASS at every
`W ∈ {6, 30, 210, 2310}`; equality at `W = 2` is structural (the only
squarefree integers `≤ 2` are `{1, 2}`, both divisors). The Pearson gap
grows monotonically: `0.028 → 0.089 → 0.177 → 0.281` at `W ∈ {6, 30, 210,
2310}` (N=10⁵), so the divisor-only restriction captures only **69% of
T_W's full Pearson with chi_P at W=2310**; the missing signal is in the
non-divisor squarefree conductors `q ≤ W`.
**Net new content:** (i) **closed-form pointwise identity on the full
domain** of `n`, generalising S191's coprime-only corollary; (ii)
**radical reduction** for non-squarefree `W`; (iii) explicit L² closed
form via prime-distribution boundary count; (iv) localisation of the
prime-discriminating signal to the *non-divisor squarefree conductors*
in the full `M_Q`. **Refines E2.1** with a pointwise identity (the
`phi(W)/W` factor is now the value of an explicit scalar field, not
merely an asymptotic rank quotient).
**Status:** BUILT, no polylog opening (the divisor-only `T_W^{div}` is
information-theoretically a wheel sieve). CLOSED_PATHS row added (S208).
EDGES.md E2.1 annotated.
**Save under:** `experiments/constructions/spike_divisor_only_wheel/`.

**Successor challenges (proposed in S208):**

**C9.a.i — Off-divisor squarefree expansion at primorial W.** Build the
residual `M_W − M_W^{div}` at primorial W = 2310 and decompose the
0.281 Pearson gap by conductor `q`. Predicted: each prime
`13 ≤ p ≤ √W` contributes `~ 1/φ(p)²` to the Pearson² gap (Selberg-
Delange–style). Cost: 1 session.

**C9.a.ii — Lean formalisation of the divisor-only identity.** The
one-line Möbius collapse plus the Hölder-Möbius reduction (L6 queue) is
a clean Lean target. Estimated 1 session. Save under
`experiments/formalisations/L6_holder_normalised_ramanujan/`.

**C9.a.iii — Λ-modulated divisor-only sum.** Replace `mu(q) c_q(n)/phi(q)`
with `μ²(q) λ(q) c_q(n)/phi(q)` (or similar parity-twisted weight) and
ask whether the divisor-only sum is again a closed-form scalar on the
residue lattice mod `W`. Cost: 1 session.

**C9.b — Higher-moment composition: T_Q correlations and Hardy-Littlewood — BUILT (S205).**
**Edges:** E2.1 + E2.2 + E2.13 + E1.6 + C9 / S191.
**Object:** Two-point function `R_h(Q, N) := ⟨T_Q(n) T_Q(n+h)⟩` with
connected variant `R_h^conn := R_h − ⟨T_Q⟩²`. Predicted closed-form
identity: `R_h^conn(Q, N) = (π(N)/N)² · (S_Q(h) − 1)` where
`S_Q(h) = Σ_{q sqf ≤ Q} μ²(q)/φ²(q) c_q(h)` is the truncated
Hardy-Littlewood twin-prime singular series.
**Falsifiers:** F1 identity ratio in [0.85, 1.15]; F2 convergence to
full HL at Q=√N; F3 odd-h asymptote → −(π/N)²; F4 h=0 self-consistency
recovers S191 L²; F5 prime correlation matches HL.
**Outcome (S205):** All five falsifiers PASS, with F1 holding to
**< 0.6 %** uniformly across 14 shifts × 8 conductors at d = 20 (two
orders of magnitude tighter than the pre-stated band). F2 (HL recovery
at Q = √N): even-h ratios within 0.2 %; F3 (odd-h asymptote): within
0.03 %; F4 (h=0 variance): within 0.01 %; F5 (prime correlation): within
~6 % (standard finite-N for primes at d=20). The identity composes
S191's pointwise spike with E2.13's HL signature: `T_Q` is a closed-form
finite-conductor approximator of HL whose two-point function reproduces
the truncated singular series exactly. The construction also recovers
C9 / S191's single-point L² statement as the h = 0 diagonal.
**Status:** BUILT, no polylog opening (cost is `O(N · |H|)`). Refines
E2.13 (cube → two-point shift case), E2.1 (single-point spike → two-point
spike correlation), C9 / S191 (single-point → two-point). CLOSED_PATHS
row added (S205). EDGES.md E2.13 / E2.1 annotated.
**Save under:** `experiments/constructions/spike_pointwise_HL_correlation/`.

**Successor challenges (proposed in S205):**

**C9.b.i — Cross-conductor off-diagonal explicit form.** F1 ratios
drift by ~0.6 % at the highest tested conductor Q = 2310, d = 20. Build
the explicit cross-conductor (gcd(q1, q2) > 1, both squarefree) sum and
verify it matches the predicted `O(Q · log Q / N)` leakage. Cost: 1
session. **Save under:** `experiments/constructions/spike_pointwise_HL_correlation/cross_conductor/`.

**C9.b.ii — Triple correlation `<T_Q · shift_{h1} · shift_{h2}>` —
BUILT (S209, mode I, B-grade refinement of E2.1 / S205 / S208).**
Test whether the three-point function reproduces the truncated HL
triple-prime singular series `S(0, h1, h2)`. This bridges S205's
two-point identity to E2.13's full `U^2` cube content via an explicit
sequence of `k`-point identities. **Save under:**
`experiments/constructions/spike_pointwise_HL_triple/`.

**Outcome (S209):** all six falsifiers PASS; F1 (primorial-W identity)
within **0.01% at d=20** (100× tighter than the pre-stated 0.5% band).
**Theorem (S209)**: for every squarefree primorial `W ≥ 1` and integers
`(h_1, h_2)`,
`⟨T_W^{div}(n) T_W^{div}(n+h_1) T_W^{div}(n+h_2)⟩_n =
(π(N)/N)^3 · ∏_{p|W} (p − ν_p(0, h_1, h_2)) · p² / (p−1)^3 =
(π(N)/N)^3 · S_HL^{(W)}(0, h_1, h_2)`. **Two independent proofs**:
(a) direct cube of S208's wheel-indicator collapse + prime-factor
inclusion-exclusion; (b) Ramanujan-Fourier prime-by-prime expansion
giving `G_p = 1 + (1/(p-1)²)[c_p(h_1) + c_p(h_2) + c_p(h_2 - h_1)] −
f_p/(p-1)³` with `f_p ∈ {(p-1)(p-2), -(p-2), +2}` for `ν_p ∈ {1, 2, 3}`,
algebraically equal to `(p − ν_p) p²/(p−1)³` across 70/70 prime ×
shift cells. F3 (HL recovery for general T_Q at Q ≈ √N) within 1%
across all 6 admissible (h_1, h_2) tested at d=20. F4 / F5 (h=0 and
h_2=0 self-consistency) within 0.001% at d=20. F6 (prime triple density
finite-N) within 5% standard. **Composes** S205 (2-point), S208 (1-point
divisor collapse), E2.1 (spike), E2.13 (cube), E2.16 (3-point HL factors
over primes — turned negative-shape into positive identity), E1.6 / E2.2
(parity sign at q=2). **No polylog opening** (cost `O(N · |H|²)`).
CLOSED_PATHS row added (S209). EDGES.md E2.1 inline annotated.

**Successor challenges (proposed in S209):**

**C9.b.iv — k = 4 four-point identity at primorial W — BUILT (S234,
mode E, B-grade refinement of E2.1 / S205 / S208 / S209).**
**Theorem (S234, proven analytically; verified to 0.06% across
admissible cells at d=20).** For every squarefree W ≥ 1 and integers
(h_1, h_2, h_3),
`<T_W^{div}(n) T_W^{div}(n+h_1) T_W^{div}(n+h_2) T_W^{div}(n+h_3)>_n
= (π(N)/N)^4 · ∏_{p|W} (p − ν_p(0, h_1, h_2, h_3)) · p^3 / (p−1)^4
= (π(N)/N)^4 · S_HL^{(W)}(0, h_1, h_2, h_3)`. Two independent proofs:
(a) S208 pointwise collapse + CRT density factorisation; (b)
Ramanujan-Fourier 4-cumulant expansion algebraically equal to the
closed form across 78/78 small-prime cells (p ∈ {2,3,5,7,11,13} × 13
shift triples) at abs_err ≤ 4.44e-16. Reduction lemmas:
h_1=h_2=h_3=0 → `<T_W^4>=(π/N)^4(W/φ(W))^3` (F4); h_3=0 → (W/φ(W))×
S209-3pt (F5). F1 (primorial-W identity) PASS at 100× tighter than
pre-stated band; F2 / F3 partial fails (1-2 of 6 admissible cells outside
pre-stated 2.5% / 3% bands at d=20 Q=2310; cross-conductor leakage at
k=4 amplified by 4-6× vs S205 k=2 — 3 disconnected pair contractions).
**Composes** S205 / S208 / S209 (k=1,2,3 → k=4 hierarchy now closed) +
E2.13 (`U^2` cube subsumed at (h_1, h_2, h_1+h_2)) + E2.16 (DPP
3-point negative shape extended to 4-tuples as positive identity) +
E1.6 / E2.2 (q=2 parity admissibility). General-k form
`G_p^{(k)} = (p − ν_p) p^{k-1} / (p−1)^k` follows by induction. No
polylog opening. CLOSED_PATHS row added (S234). EDGES.md E2.1 inline
annotated. See `experiments/constructions/spike_pointwise_HL_quad/`,
`archive/sessions/session234_c9biv_4point_HL_identity.md`.

**Successor challenges (proposed in S234):**

**C9.b.iv.α — Cross-conductor leakage at k=4: closed-form bound.**
S234's F2/F3 partial fails (worst 5.7% at (6,10,12)) reflect 3
disconnected pair contractions at the 4-point level. Predicted scaling
`O(Q · log^j Q / N)` for some j ≤ 3 (refining S205's `O(Q log Q / N)`
at the 2-point level). Quantify analytically and verify at d ∈ {18, 20,
22} for general Q. Cost: 1 session. **Save under:**
`experiments/constructions/spike_pointwise_HL_quad/cross_conductor_k4/`.

**C9.b.iv.β — Lean formalisation of the 4-point identity at primorial
W.** Proof is `S208 + CRT + factor identity`. Add to L1 / L6 queue
after C9.b.iii (S205 2-point) and C9.b.vi (S209 3-point). Cost: 1-2
sessions.

**C9.b.v — General-k closed form `G_p^{(k)} = (p−ν_p) p^{k-1}/(p−1)^k`.**
With k=1,2,3,4 verified analytically (S208/S205/S209/S234), the inductive
proof via S208 collapse on the k-fold indicator product is essentially
mechanical. Verify with k=5 at one (h_1,h_2,h_3,h_4) tuple as a
check + write up the inductive theorem. Cost: 0.5 session.

**C9.b.vii — Hoffman / triple-MZV interpretation of `f_p^{(4)}` per
multiplicity profile.** The closed-form values per ν_p and multiplicity
split — (p-1)[1+(p-1)^3]/p for (4,); (p-3) for (2,1,1); -3 for
(1,1,1,1); (1/p)[(p-2)+2(p-1)^2] for (2,2); etc. — are explicit small
polynomials in p. Possibly bridges to D38 (S233) prime-MZV antisymmetric
A(s,t) at weight 4 via Hoffman algebra. Cost: 1 session.

**C9.b.v — General-Q off-diagonal calibration.** The F2 over-count
at non-primorial Q (predicted `prod_{p sqf primes ≤ Q} G_p` minus
empirical) is due to subsets of primes whose product exceeds Q being
absent from the squarefree-cap sum. Quantify analytically. Cost: 1
session. **Save under:** `experiments/constructions/spike_pointwise_HL_correlation/general_Q_calibration/`.

**C9.b.vi — Lean formalisation of the 3-point identity.** ~30-line
target via S208's Möbius collapse + 3-coprimality factorisation. Add to
L1 / L6 queue after L1.a (S205 2-point) target. Cost: 1-2 sessions.

**C9.b.iii — Lean formalisation of the two-point identity.** The
identity `<T_Q(n) T_Q(n+h)>_n = (π(N)/N)² · S_Q(h)` is a one-step
character-theoretic computation: Ramanujan-Fourier expansion +
diagonal-conductor selection. Tractable Lean 4 target; add to L1 / L6
queue. Cost: 2-3 sessions.

### C10 — Wheel-coprime cumulative Liouville identity — **BUILT (S219, paradigm-shift mode, B-grade refinement of E1.6 / E2.10)**
**Edges:** E1.6 + E2.10 + E2.1 / E4.1 / S208.
**Object:** `L_W(x) := Σ_{n ≤ x, gcd(n,W)=1} λ(n)` — the cumulative
Liouville sum restricted to wheel-W-coprime integers. The cumulative-
side analogue of S208's pointwise indicator `T_W^{div}(n) = (π/N) ·
(W/φ(W)) · [gcd(n,W)=1]`.
**Closed-form identities (proved analytically, integer-exact verified):**
* (T1) `L_W(x) = Σ_{d | rad(W)} L(⌊x/d⌋)` for every W ≥ 1.
* (T2) `L_W(x) mod 2 = (Σ_{d | rad(W)} ⌊x/d⌋) mod 2` (combines T1 + E2.10).
* (T3) Wheel-graded prime bisection lift: `π_W(x) = (1/2) Σ_{d | rad(W)}
  (μ(d)⌊x/d⌋ − L(⌊x/d⌋)) − C_{3,W}(x)`, where π_W(x) := π(x) − π_W^{div}(x).
**Outcome (S219):** All six pre-stated falsifiers PASS integer-exactly
at 16 W (sqfree primorial / sqfree non-primorial / non-squarefree) ×
~200 x grid points at N = 10⁴. F1 (T1 pointwise), F2 (radical reduction
`L_W = L_{rad(W)}`), F3 (T2 mod-2), F4 (T3 wheel-bisection), F5 (RHS
divisor count = `2^{ω(rad(W))}`), F6 (mod-q lift at q ∈ {2, 3, 4, 8}).
**Algorithmic content:** `2^{ω(W)}`-call oracle reduction `L_W → L`,
polylog in W for primorial W (since `2^{ω(W)} = W^{O(1/log log W)}`).
Mod-2 unconditionally polylog (T2). Does NOT bypass the unconditional
`O(x^{2/3})` ceiling for `L(x)` — the bottleneck stays on `L`.
**Symmetric counterpart to S208** (indicator-side wheel-collapse):
S208's `T_W^{div}` is pointwise-side, polylog-cheap, finite-state
on coprime cosets; S219's `L_W` is cumulative-side, retains the full
oscillatory content of `λ` decomposed into `2^{ω(W)}` scaled copies.
**Status:** BUILT, no polylog opening. Refines E1.6 (bisection) and
E2.10 (parity) inline. CLOSED_PATHS row added (S219).
**Save under:** `experiments/constructions/wheel_coprime_liouville_identity/`.

**Successor challenges (proposed in S219):**

**C10.a — Selberg-Delange decomposition of the L_W − (φ(W)/W) L(x) heuristic.**
The heuristic `L_W(x) ~ (φ(W)/W) L(x)` is folklore in PNT-in-AP. T1
gives an EXACT rewriting:
`L_W(x) − (φ(W)/W) L(x) = Σ_{d | rad(W), d > 1} (L(⌊x/d⌋) − ⌊x/d⌋ · L(x)/x)`.
Each summand is `o(x/d)` by PNT for λ. Quantify the magnitude of
`Σ_{d | rad(W), d > 1}` empirically at scales `N ∈ {10⁵, 10⁶, 10⁷}` and
across primorial W ∈ {2, 6, 30, 210, 2310}. Predicted: heuristic error is
dominated by the `d = 2` term `L(⌊x/2⌋) − ⌊x/2⌋ · L(x)/x` for primorial
W. Cost: 1 session. **Save under:**
`experiments/constructions/wheel_coprime_liouville_identity/selberg_delange/`.

**C10.b — Mod-4 lift of L_W.** T2 gives the closed form for L_W mod 2
(via E2.10). Mod 4 is the E1.5 / T6 hard-zone bit boundary.
Apply T1 mod 4 directly: `L_W(x) mod 4 = Σ_{d | rad(W)} L(⌊x/d⌋) mod 4`.
Empirically test whether at primorial W, the `2^{ω(W)}`-fold "blurring"
of the L mod 4 hard-zone bit gives any algorithmic leverage at small
scale. Predicted: no leverage (E1.5 saturation per E1.5 / T6), but
worth empirical confirmation. Cost: 1 session. **Save under:**
`experiments/constructions/wheel_coprime_liouville_identity/mod4_lift/`.

**C10.c — Lean formalisation of T1 and T3.** T1 is a one-line Möbius
inversion + completely-multiplicative collapse. T3 substitutes T1 into
the bisection. Tractable Lean 4 target; add to L1 / L6 queue. Cost:
1-2 sessions.

### C11 — D17.b: Discrete Morse on the squarefree-only divisibility Hasse diagram — **BUILT (S422, mode E, B-grade refinement of D17 / S232)**
**Edges/sources:** D17 closed at S232 (mode E, B-grade) with sharp
identity `collapses(H_N) = (π(N) − π(N/2)) + Π_pow(N) + 1` for the
full divisibility Hasse 1-skeleton. The squarefree-only restriction
is a distinct combinatorial object: vertices `{n ≤ N : μ(n) ≠ 0}`;
covering relations `(m, mp)` with `m, mp` squarefree, `p` prime,
`p ∤ m`, `mp ≤ N`.
**Outcome (S422):** Greedy random Morse collapse halts in **exactly
one wave** at every `N ∈ {64, 128, 256, 512, 1024, 2048, 4096, 8192}`.
Sharp identity (proven analytically + verified pointwise):
> **Theorem (S422).** `collapses(H_N^sqfree) = π(N) − π(N/2)` exactly,
> with `ε(N) ≡ 0`. Equivalently `m_0(H_N^sqfree) = |V_sqfree(N)| − (π(N) − π(N/2))`,
> `m_1(H_N^sqfree) = |E(H_N^sqfree)| − (π(N) − π(N/2))`.

Two-line proof: vertex-degree case-analysis shows wave-0 leaves are
*exactly* primes `p ∈ (N/2, N]` (each has unique edge to vertex 1);
peeling drops `deg(1)` from `π(N)` to `π(N/2) ≥ 2` for `N ≥ 6`; no
other vertex's degree is modified, so wave 1 is empty. Leaves
mutually non-adjacent, so Morse rigidity is automatic.

Sharper than D17 by removing BOTH the `Π_pow(N) ~ Θ(√N/log N)`
prime-power term AND the constant `+1` chained-collapse residual.
F1 (polylog A-grade) FAILS: `m_0/|V| → 1` monotone. F2 (ER baseline
match) FAILS: divisibility-Hasse is ≈ 4 % more collapsible than
matched-density ER (vs D17's 0.5–2 %); squarefree restriction
*amplifies* the divisibility-vs-random distinguishability gap. F3
HOLDS exactly. F4 (Morse rigidity) HOLDS analytically across 200
seeds at every measured N. **Status:** BUILT, no polylog opening
(circular: `m_0` reduces to `π(N) − π(N/2)`, no easier than `π(N)`).
**Refines D17 inline.** CLOSED_PATHS row added (S422). See
`experiments/topological/d17b_squarefree_morse/`.
**Save under:** `experiments/topological/d17b_squarefree_morse/`.

**Successor challenges (proposed in S422):**

**C11.a — D17.b.i full order complex.** Lift from 1-skeleton to the
full chain complex of the squarefree-poset order complex. Cell count
`~ Σ_{n sqf ≤ N} ω(n)!`, tractable for `N ≤ 256`. The unrestricted
Boolean lattice's order complex is shellable with one critical top
cell; the truncation `∏ p ≤ N` may reduce chain-level `m_0` below
`Θ(N)`. **A-grade target if** `m_0(chain) = O(polylog N)` despite the
1-skeleton failure. Cost: 1-2 sessions. Save under
`experiments/topological/d17b_squarefree_morse/order_complex/`.

**C11.b — D17.b.ii multiplicative-indicator generalisation.** Replace
`μ²(n)` by other multiplicative indicators (k-free numbers,
k-almost-primes, smooth-y-rough numbers). Each defines a sub-Hasse of
`(Z, |) ∩ [1, N]`. Predicted: same `(π(N) − π(N/2))`-leading-term
identity with a different second-order correction matching the
indicator's level-2 structure. Cost: 1 session per indicator family.
Save under
`experiments/topological/d17b_squarefree_morse/multiplicative_indicators/`.

**C11.c — D17.b.iii Lean formalisation.** Two-line proof
(vertex-degree case-analysis + wave-1 emptiness lemma). Tractable
Lean 4 target; add to L1 / L6 queue. Cost: 1-2 sessions. Save under
`experiments/formalisations/L_d17b_morse_identity/`.

### C12 — D38.a: Prime-MZV antisymmetric A(s,t) at weight ≥ 8 with Hoffman-irreducible MZV basis (NEW from S233 follow-up)
**Edges/sources:** E7.23 (D38 closed S233, mode I, B-grade) — PSLQ
NO_RELATION on `A(s,t) := ζ_P(s,t) − ζ_P(t,s)` against the COMPLETE
weight-`(s+t)` basis at `(s,t) ∈ {(2,3), (2,4), (2,5), (3,4)}` (weights
5/6/7). The weight-8 result was flagged INC because the basis missed
the irreducible double zeta `ζ(3,5)` and the Hoffman-depth-3 generators
`ζ(3,3,2), ζ(2,3,3), ζ(3,2,3)`. At weight 9 the gap widens further
(`ζ(3,3,3)`, `ζ(2,2,5)`, etc.).
**Question:** Compute the depth-3 Hoffman generators `ζ(3,3,2)`,
`ζ(2,3,3)`, `ζ(3,2,3)` numerically to ≥ 50 digits via mpmath nested
sums (or use Bigotte/Goncharov shuffle/stuffle reductions). Add to
the weight-8 basis. Re-run PSLQ on `A(2,6), A(3,5)`. At weight 9 add
`ζ(3,3,3)` and re-run PSLQ on `A(2,7), A(3,6), A(4,5)`.
**Why this is a real successor (not a duplicate of D38 closure):** the
weight-8 basis was the INC-flagged tier in S233. Adding the missing
irreducible MZV either:
- **Closes the gap** → primes ARE in Brown's algebra at higher depth
  via Hoffman generators (would require an explicit reduction to
  re-open D38 as E-mode); or
- **Confirms** the structural negative at higher weights →
  strengthens E7.23 by ruling out the depth-≥3 saturation hope.
**A-grade target:** an explicit closed-form
`A(2,6) = α_1 ζ(3,3,2) + α_2 ζ(2,3,3) + α_3 ζ(3,2,3)
         + Σ β_j (Mertens monomial of weight 8)`
with rational `α_j, β_j` would partially rehabilitate D38 and yield
polylog evaluators for prime-MZVs at infinite depth via Brown's
recursive Hoffman reduction → polylog access to certain prime
correlation sums.
**B-grade outcome:** PSLQ NO_RELATION at weight 8 (and 9) with the
expanded depth-3 basis — confirms E7.23 universally, lifts the INC
flag.
**First step:** Implement nested-sum evaluator for `ζ(s_1, s_2, s_3)`
in mpmath (`Σ_{n_1>n_2>n_3≥1} 1/(n_1^{s_1} n_2^{s_2} n_3^{s_3})`)
with tail-correction. Compute `ζ(3,3,2), ζ(2,3,3), ζ(3,2,3),
ζ(3,3,3)` to 50 dps. Re-run d38_prime_mzv.py PSLQ with expanded
basis at weights 8, 9.
**Save under:** `experiments/algebraic/d38_prime_mzv/` (extend the
existing experiment; add `d38_higher_weight_results.md`).
**Cost:** 1-2 sessions.

---

## §2. Frame-Shift Questions

These are research-grade open questions that *re-frame* the project away
from "polylog π(x)" while remaining in the same neighbourhood. The point
is to escape the local minimum: 67+ sessions on one frame have been
exhaustive, but adjacent frames are essentially untouched.

### F1 — Per-bit polylog extraction — **PARTIALLY CLOSED (S146, mode E, B-grade refinement of E1.3)**
**Question:** Find J such that the J-th bit of `p(n)` is computable in polylog *for fixed J independent of N*. E1.3 says bits 0..0.6N match `round(R^{-1}(n))`. So J = 0 is trivially polylog. What is the LARGEST J for which `bit_J(p(n))` is provably polylog?
**Why this is unblocked:** The whole-number question is closed; the per-bit question has an obvious YES answer for small J and an unknown answer for medium J. The boundary is the research target.
**First step:** Read `novel/carry_propagation_boundary.md`. Build a per-bit accuracy curve as a function of J for N up to 30.

**S146 outcome (refinement, NOT closure):** Direct bit-position
measurement on N = π(2·10⁸) ≈ 1.1·10⁷ primes refines E1.3 with
three new structural facts:
- **bit_0(p(n)) is trivially polylog** (= 1 for n ≥ 2). E1.3's
  R⁻¹-agreement metric mis-classifies it as "hard" (ag = 0.500)
  because Li⁻¹'s LSB is uniform mod 2, but the bit IS deterministic.
- **RH-scale anti-correlation valley**: ag_Li(J) drops to ~0.36
  at J* = ⌊log₂(p(N))/2⌋ — the bit-level signature of Skewes-
  direction error sign in `Li(x) > π(x)`. Position shifts with
  N as predicted; magnitude stable across N ∈ {10⁷, 5·10⁷, 2·10⁸}.
- **PNT first-order ≠ Li⁻¹ at bit level**: top-bit ag_PNT < 0.95
  while ag_Li → 1.0; full PNT closed form is required.
The original F1 question (largest J for provable polylog) remains
**OPEN** — the new content is a structural map (LSB triviality,
RH-scale valley, PNT-vs-Li bit-level inequivalence) that constrains
where to look for non-trivial polylog opportunities. See
`experiments/wildcard/bit_J_pn_polylog_map/bit_J_pn_polylog_map_results.md`,
`archive/sessions/session146_f1_bit_J_pn_polylog_map.md`.

**Successor challenges (proposed in S146):**

**§F1.a — CLOSED (S199, mode E, B-grade refinement of E1.3).**
Cross-modulus generalisation of the S146 RH-shadow valley to bases
m ∈ {3, 5, 6, 30, 210} at L = 2·10⁸. **Hypothesis CONFIRMED**: the
dip position `J*(m) = argmin_J ag_Li(m, J)` matches `⌊log_m(L)/2⌋`
**exactly (Δ = 0)** for all 5 cross-modulus bases. Five F-verdicts:
F1.a-1 HOLDS unanimously; F1.a-2 (sub-baseline by ≥ 5%) HOLDS
unanimously; F1.a-3 (random regime far from J*) HOLDS; F1.a-4
(modal shift at +1) REJECTED, **revised to F1.a-4'**: the modal
shift at J*(m) is `s* ≈ ⌊⟨e⟩/m^J*⌋ mod m` (S146's "+1 mod 2" is
the special case `⟨e⟩/2^J* = 0.66 < 1`); F1.a-5 HOLDS exactly
(`H_p(J=0) = log_2(φ(m))` for primorial m). **Net new content:**
(a) the RH-shadow valley is *m-adic universal*, not a base-2
artefact; (b) relative dip `rel(m) = ag·m` deepens monotonically
with the conductor — 0.722 (m=2) → 0.543 (m=3) → 0.521 (m=6) →
0.035 (m=30) → 0.010 (m=210), so the Li⁻¹ predictor is essentially
deterministic-wrong at the m=210 half-conductor digit; (c) m=5
is structurally shallow (mid-wrap modal shift); (d) primorial-m
digit-0 entropy = log_2(φ(m)) exactly. **Refines E1.3 inline**.
See `experiments/wildcard/bit_J_pn_cross_modulus/`,
`archive/sessions/session199_f1a_cross_modulus.md`.

**Successor challenges (proposed in S199):**

**§F1.a.i — CLOSED (S218, mode E, B-grade refinement of E1.3).**
Tabulated `rel_emp(m)` across `m ∈ {2..30}` at L = 2·10⁸ and tested
the proposed Gaussian-RH-shadow closed-form against three concrete
predictions: (Y) Gaussian-Y, (R) Empirical-e + uniform-r, (GR)
Gaussian-e + uniform-r. **Outcome: A-grade refuted, B-grade
established.**

**EXACT identity** (verified to 0.04 % mean abs error across 29
moduli): under uniformly-distributed `L_n mod m^J*`,
`ag_emp(m, J*) = E_n[(1−frac_n)·𝟙[q_n ≡ 0 (mod m)] +
frac_n·𝟙[q_n ≡ m−1 (mod m)]]` with `q_n = ⌊e_n/m^J*⌋`,
`frac_n = (e_n mod m^J*)/m^J*`.

**Gaussian-Y is APPROXIMATE.** Mean abs error 0.0962, holds to
6–15 % for `σ_Y ≤ 1.5` cells but fails by 30 %–530 % on
`σ_Y ≥ 2` cells (`m ∈ {12, 13, 27, 28, 29, 30}`) because empirical
e is bounded `[0, 21648]` while Gaussian e leaks mass outside.

**NEW PEAK structure.** `rel_emp(m)` non-monotone, `0.025 ≤
rel_emp(m) ≤ 6.32`, peaking at `m = 24` (μ_Y = 1.28, σ_Y = 0.42).
The "RH-shadow valley" reframes as **"RH-shadow phase alignment"**:
agreement ratio is m·P[⌊Y⌋ ≡ 0 mod m], controlled by μ_Y(m) mod m
relative to the integer lattice with σ_Y as a diffusion parameter.
F5/F6 (m=5 mid-wrap) PASS; F7 (composite=parent) PARTIAL (3/6).
**Refines E1.3 inline.** See
`experiments/wildcard/bit_J_pn_dip_scaling/`,
`archive/sessions/session218_f1ai_dip_scaling.md`.

**Successor challenges (proposed in S218):**

**§F1.a.i.α — Sub-Gaussian tail correction.** Replace the Gaussian
e model with a truncated Gaussian on `[0, e_max]` (or beta-Cramér
mixture matched to the empirical skew = -0.108 and the |e| < 22000
bound). Re-fit the closed form. Predicted: closed-form error falls
from O(10²) % to O(1) % across `m ∈ {2..30}`. **A-grade if** the
tail-corrected closed form matches empirical to within 1 % on every
m. Cost: 1 session.

**§F1.a.i.β — Cramér-asymptotic prediction.** Under the Cramér model
`e ~ √p_N · N(c, σ²)` with `c, σ` constants, derive the asymptotic
`rel(m, p_N → ∞)` as a function of m only (J* absorbed via
`m^J* = √p_N`). Test cross-scale at L ∈ {10⁷, 5·10⁷, 2·10⁸, 10⁹}
(sympy/primecount for L = 10⁹). The Cramér prediction is L-independent
in the limit; cross-scale stability would confirm. Cost: 2-3 sessions.

**§F1.a.i.γ — CLOSED (S238, mode E, B-grade refinement of E1.3).**
Phase-diagram sweep across `m ∈ {2..100, 110, 120, 140, 170, 200, 250,
300, 400, 500, 700, 1000, 1500}` (111 cells, both J*=2 and J*=1) at
L = 2·10⁸ in `(α, σ_norm)` coordinates produces four new structural
facts: (a) wrapped-Gaussian density regime is σ_Y ≥ 2 valid (mean
|Δrel| = 0.21 over 35 cells), σ_Y < 2 fails (3.47 over 76 cells);
(b) U-shape against α decile-binned (max=4.67 at decile 0, min=0.05
at decile 4) HOLDS; (c) peak ridge along α ≈ 0 reaches `rel_emp(110)
= 22.37` (3.5× S218's m=24 max); (d) J*=1 bounded-support trough
mechanism — deepest trough rel_emp(170)=0.005 at α=0.376 (NOT α=0.5),
caused by truncation of Gaussian-Y's right-tail at the bounded-e
support [0, 21648]; sharpens S218's Gaussian-Y downgrade to 17×
worst-case at J*=1 (m=200). Refines E1.3 inline. F-verdicts: F_β
HOLDS, F_ε HOLDS, F_α/F_δ/F_ζ/F_η FAIL informatively. See
`experiments/wildcard/bit_J_pn_phase_diagram/`,
`archive/sessions/session238_f1aig_phase_diagram.md`.

**Successor challenges (proposed in S238):**

**§F1.a.i.γ.i — Truncated-Gaussian closed form for J*=1 cells.**
Replace Gaussian e by truncated normal on [0, e_max] (or beta-Cramér
mixture matched to the empirical e support [0, 21648] and skew
−0.108). Re-derive the closed form. Predicted: J*=1 trough errors
collapse from 17× (m=200) to ≤ 30 %. **A-grade if** the tail-
corrected closed form matches empirical to within 1 % across all
m ∈ [2, 1500]. Cost: 1 session. (Subsumes §F1.a.i.α.)

**§F1.a.i.γ.ii — Asymptotic peak height.**
Conjecture: `max_m rel_emp(m, J*=2) → 1/(σ_norm·√(2π))` as m
approaches the J*=2 ↔ J*=1 boundary at `m ≈ p_N^{1/4}`. Verify
peak ridge scaling at `L ∈ {10⁷, 5·10⁷, 2·10⁸, 10⁹}`. Predicted:
peak height grows as L^{1/8} or similar. Cost: 1-2 sessions.

**§F1.a.i.γ.iii — Symmetry restoration via L-scaling.**
F_δ failed because the natural α(m) distribution at fixed L is
α-skewed. Test whether scanning `L ∈ {10⁶..10⁹}` produces mirror
α-coordinates for the same m, and whether `rel_emp(m, L)` and
`rel_emp(m, L')` with `α(m, L) ≈ 1 − α(m, L')` match within 30%.
Cost: 1 session.

**§F1.a.ii — Higher m primorial extension.** Extend to m ∈ {2310,
30030} where `J*(m) ∈ {0, 1}` at L = 2·10⁸. Test whether the dip
is resolvable at J = 0 (where the predictor and truth are *both*
on the full conductor scale) — a regime where E1.3's per-bit
decomposition breaks down because the "digit" is the entire
conductor's reduced-residue class. Cost: 1 session.

**§F1.a.iii — Cross-zero-truncation.** Compose with §D43.c (S160
K-truncated explicit-formula residual). Hypothesis: subtracting
the first K explicit-formula zero contributions from Li⁻¹(n) should
*cancel* the RH-shadow valley at every base m once K is large
enough. Empirical test: re-measure m-ary dip on R_K predictor
agreement for K ∈ {0, 1, 5, 25, 100, 1000}. Direct cross-validation
of S195/S160's explicit-formula picture at the cross-modulus bit-
shadow level. Cost: 1-2 sessions.

**§F1.b — Magnitude scaling.** The dip magnitude (ag_Li(J*) ≈ 0.36)
is stable at L ∈ {10⁷, 5·10⁷, 2·10⁸}. Test up to L = 10⁹ (sympy or
primecount) to see whether the magnitude tracks (or saturates at)
the Skewes-region transition. Predicted: magnitude stable until
the first sign-change of `Li(x) - π(x)` at x ≈ 10^316 (Bays-Hudson)
which is unreachable. So the valley should be **a fixed empirical
constant** in the practical regime. Verifying this constancy at
larger L is direct and tractable. Cost: 1 session.

### F2 — π_+(x) := π(x) mod p^k for fixed prime p, growing k — **PARTIALLY CLOSED (S69)**
**Question:** π(x) mod 2 is pseudorandom (35 measures). Is π(x) mod 4 also pseudorandom in 35 measures? Mod 8? At what point does the "free-bit" argument break, if ever? Is the *information rate* of π(x) mod 2^k the constant 0.537·k bits per step (E1.5 says yes) or does it saturate?
**Why novel:** E1.5 was tested only at small q. The saturation behaviour is not in the literature.
**First step:** Extend E1.5's measurement to q = 2^k for k = 1..10.

**S69 outcome (information-rate side, CLOSED with closed form):**
`experiments/information_theory/pi_mod_2k_saturation/pi_mod_2k_saturation_results.md`.
For m = 2^k with k = 1..10 at X up to 10^7:
`H(π(x) mod m | π(x-1) mod m) = h_2(π(X)/X) + O(1/π(X))` exactly
(verified to 7 decimals), m-independent for m ≪ π(X). This **falsifies the
linear-in-k framing** and **sharpens E1.5** with a closed form. The
"0.537" was the X = 10^4 value of `h_2(π(X)/X)`; the constant decays as
X grows, asymptotically to 0.
**Still open (other side of F2):** the 35-measure pseudorandomness battery
applied to π(x) mod 2^k for k > 1 is not yet run. If a future agent picks
this up, the information-rate question has been settled — focus on
non-information measures.

### F3 — Algorithms for π(x) given an oracle
**Question:** Suppose you have a TC⁰ primality oracle (E5.1, conditional on BPSW). What is the minimum complexity of π(x) given this oracle? Aggarwal-style binary search gives O(sqrt(x) log^4 x). Is sub-sqrt achievable with the oracle?
**Why novel:** Most circuit-complexity work treats π(x) and primality together. Separating them — assuming primality is free, asking about counting — is a clean sub-problem nobody has formalised.
**First step:** Define the oracle Turing-machine model precisely. Re-derive Aggarwal in this model. Look for sublinear improvements.

### F4 — Smooth-π(x) and its complexity
**Question:** Define π_smooth(x, B) = "number of B-smooth integers ≤ x with no large prime factor." Its complexity is well-studied (Buchstab function). Now define π_BD(x) = π(x) - π_smooth(x, x^{1/3}) = "primes minus B-smooth count." Is π_BD(x) more or less compressible than π(x)? Does subtracting a "computable noise" change the residual's spectral / information-theoretic structure?
**Why novel:** The decomposition π = R(x) + (π - R(x)) has been beaten to death. The decomposition π = π_smooth + π_BD has not.
**First step:** Compute π_BD(x) for x ≤ 10^7. Apply the 35-measure pseudorandomness battery.

### F5 — Find a function f with π-comparable statistics that IS in TC⁰
**Question:** The pseudorandomness of π(x) mod 2 (35 measures) is taken as evidence that it's outside TC⁰. But pseudorandomness alone doesn't prove circuit-hardness (cryptographic PRGs are pseudorandom AND in TC⁰). Find a concrete function f: ℕ → {0,1} that:
  (a) passes all 35 pseudorandomness measures applied to π(x) mod 2;
  (b) is provably in TC⁰.
If such f exists, it weakens the pseudorandomness-→-hardness argument. If you can argue no such f exists at this level of pseudorandomness, you've sharpened the argument.
**Why novel:** This is the dual of "is π in TC⁰?" — and is more tractable.
**First step:** AES output stream restricted to {0,1}, sampled at the same density as π(n) mod 2. Apply the 35 measures.

### F6 — Compute π(x) for x in a parametric family — **PARTIALLY CLOSED (S246, mode E, B-NEGATIVE on dyadic per-query)**
**Question:** Existing algorithms compute π(x) for arbitrary x. What is the complexity of computing `π(2^k)` for k = 1..N, all at once? The shared structure (powers of 2) might admit polylog amortised cost.
**Why novel:** Batch / parametric / structured-input variants of π(x) are essentially unstudied.
**First step:** Compute π(2^k) for k = 1..40 using primecount. Look for differences `π(2^{k+1}) - π(2^k)` and see if they admit a closed form or a polylog recurrence.

**S246 outcome (PER-QUERY side, B-NEGATIVE).** Three independent
structural tests on `(π(2^k))_{k=1..56}` (OEIS A007053) all return
random-baseline values within Monte Carlo noise:
- **F1 — sign-sequence linear complexity.** BM(sign(r_R)) = 28 against
  MC random-shuffle null `mean=28.25, std=1.01, [p05=27, p95=30]`.
  Sign sequence of `r_R(k) = π(2^k) - R(2^k)` is statistically
  identical to random.
- **F2 — autocorrelation of normalised residual.** `c_R(k) =
  r_R(k)/2^{k/2}`, max `|ac|` over lags 1..10 is **0.283** at lag 1
  vs MC iid-Gaussian `p999 = 0.437`. Bonferroni-corrected p ≈ 0.4.
  Not significant. (NB: the v1 Li-baseline gave spurious lag-1 = 0.87
  driven by the smooth `-Li(√x)/2` Möbius-series correction; corrected
  by switching to R-baseline.)
- **F3 — HKM dyadic speedup.** `sympy.primepi(2^k) / primepi(2^k ± 1)`
  cold-start subprocess timing ratio = **1.013** averaged over
  k ∈ {28, 30, 32}, range [0.998, 1.047]. Meissel-Lehmer cost is
  x-magnitude-driven, not x-form-driven.

**Net new content.** Empirical B-NEGATIVE closure of the per-query
dyadic question with two independent mechanisms documented:
(i) Weyl equidistribution of `γ_n · k log 2 mod 2π` makes phases
random across k regardless of dyadic structure; (ii) Lucy / Meissel-
Lehmer outer loop processes `{n ≤ √x}` and `{n smooth up to x^{2/3}}`
which carry no binary-representation structure. **Refines E1.1
inline.** See `experiments/wildcard/dyadic_pi_structure/`,
`archive/sessions/session246_f6_dyadic_pi_structure.md`. Successor
challenges below.

**Successor challenges (proposed in S246):**

**§F6.a — CLOSED (S426, mode E, B-grade refinement of E1.1).**
Cross-modulus structural test on `x_k = m^k` at m ∈ {3, 5, 6, 10, 30}.
**B-NEGATIVE shape carries over universally**: F1 (BM ≤ MC.p05) and F2
(max|ac| ≥ MC.p999) FAIL for every m. K-budgets {22, 15, 14, 28, 8};
total wall 30.6s.

Bonus γ_1-cosine ceiling: `|emp_lag1_ac(m)| ≤ |cos(γ_1·log(m) mod 2π)|
+ 0.05` empirically saturates at 5/5 cross-modulus cases. m=6 identified
as near-resonance (`φ_1(6) = 0.193 rad`, ceiling cos = +0.981, empirical
lag-1 = 0.529 = 54% of ceiling) — even this near-saturation case stays
below MC p999 = 0.687. Sign agreement of `cos(φ_1)` vs `emp_lag1`:
3/5 m. **Refines E1.1 inline**. See
`experiments/wildcard/cross_modulus_pi_structure/`,
`archive/sessions/session426_f6a_cross_modulus.md`.

**Successor challenges (proposed in S426):**

**§F6.a.i — γ_1 ceiling tightness scan.** Extend the ceiling
diagnostic to m ∈ {2..50}; conjecture: m^* ∈ {6, 14, 36}
(`e^{2π·j/γ_1}` integers) are strongest near-resonances. Test
each up to K_m = 14 anchors. Cost: 1 session.

**§F6.a.ii — γ_2 decoherence rate.** Predicted `|emp_lag1_ac| → 0
as K_m → ∞` at rate `~K_m^{-1/2}`. Test via sub-windowing m=10
into rolling K=14 windows; or extend to K_m=40 via primecount.
Cost: 1-2 sessions.

**§F6.a.iii — Higher-lag ceiling breakdown.** m=6 lag-1 = 0.529
matches γ_1 ceiling, but lag-2 = -0.049 vs γ_1-prediction +0.926
fails — γ_1-only model collapses at lag ≥ 2. Map breakdown
boundary in (m, ℓ) cells. Cost: 1 session.

**§F6.b — Batched-on-k amortisation (Thread 5 shape transposed).**
For `x_1 = 2^{k_1}, ..., x_M = 2^{k_M}`, the explicit-formula zero
table `(γ_n)_{n=1..K}` is independent of k. Computing all `Re(R(x_i^ρ_n))`
for i = 1..M and n = 1..K costs `O(M·K)`; per-query amortised
`O(K)`. For `K = polylog(max x_i)`, per-query = polylog. **Exactly
Thread 5 / S224 Correlation Dichotomy shape transposed to the
dyadic setting.** Cost: 2-session arc.

**§F6.c — Lower bound from output bit-budget.** Outputting all
π(2^k) for k = 1..N requires `Σ_k log₂ π(2^k) = Σ_k k − O(log k) =
Ω(N²)` bits. Make rigorous as a CLOSED_PATHS row under "lower
bounds via output size". Cost: 0.5 session.

---

## §3. Lean 4 Formalisation Queue

Goal: produce verifiable proofs of the project's main results. Lean
forces precision — vague claims fail to type-check, and proving an edge
formally often surfaces hidden structure. Each formalisation is a
permanent project asset and a publishable artifact.

**File layout:** `experiments/formalisations/<edge_id_or_result>/<name>.lean`
plus `<name>_notes.md`.

**Discipline:** the Lean file must type-check (run `lake build` or
`lean --run`). If it doesn't, the session is an in-progress formalisation,
not a closure. Save the in-progress state to `RESEARCH_AGENDA.md`.

### L1 — E2.1: MPS bond-dim identity (HIGHEST priority) — **IN PROGRESS (S76, S83, S98, S99, S106, S107, S117, S122, S128, S129, S137, S143, S144, S151, S152)**
**Statement:** For χ_P : [1, W^d] → {0,1} the prime indicator reshaped in
base W ≥ 2, for every cut 1 ≤ j < d:
```
   rank M^{(j)} = min( W^j , φ(W)·W^{d-j-1} + 1 )
```
where M^{(j)} is the (W^j × W^{d-j}) unfolding.
**Why this first:** clean clean statement, proof in `novel/mps_bond_dimension.md`
is short, no deep number theory needed (just CRT + φ definition).
**Estimated effort:** 1-2 sessions.
**Save under:** `experiments/formalisations/E2_1_mps_bond_dim/`
**S76 progress:** Lake project + mathlib `v4.30.0-rc2` installed; theorem
statement + 5-lemma skeleton typechecks (`lake build` succeeds);
`rank_le_min_dim`, `row_support_coprime`, `live_columns_count`,
`upper_bound` fully proved (no `sorry`, no `axiom`); main theorem
`mps_bond_dim` reduced to `Nat.le_antisymm` of the lemmas (3-line
term-mode proof, no `sorry`).
**S83 progress:** Closed `lower_bound` itself (sorry-free 6-line proof)
by introducing a new declaration `exists_invertible_submatrix` that
isolates the prime-density content as `∃ ρ σ, IsUnit (submatrix ρ σ)`.
The reduction uses mathlib's `Matrix.rank_of_isUnit` +
`Matrix.rank_submatrix_le`. Net: the only remaining `sorry` is now on
the pure prime-existence existential, not on the rank-bound logic.
**S98 progress:** Closed corner case `(W=2, j=1)` unconditionally via
Bertrand. **S99 progress:** Closed the orthogonal corner case
`(W=2, d=j+1)` *without even Bertrand* — only `Nat.prime_two`,
`Nat.prime_three`, decidable non-primality of 1 and 4. The matrix has
only 2 columns; rows {0, 1} with column-swap give the identity matrix.
Together S98+S99 give unconditional Lean proofs of `mps_bond_dim` on
the entire `(j, d - j)` boundary for `W = 2`. **S106 progress:**
Extended Route A''' to `W = 3`: `mps_bond_dim_W_eq_3_d_eq_j_plus_1 :
(unfolding 3 (j+1) j).rank = 3` for every `j ≥ 1`, sorry-free, via a
`3 × 3` determinant computation with `Matrix.det_fin_three`. **S107
progress:** Extended Route A'''' to `W = 4`:
`mps_bond_dim_W_eq_4_d_eq_j_plus_1 : (unfolding 4 (j+1) j).rank = 3`
for every `j ≥ 1`, sorry-free. The matrix has 4 columns but `R = 3`
(column 3 is `chiP` at multiples of 4, all zeros), so this is the
**first orthogonal corner where `rank_le_width` is not tight** —
required citing the general `upper_bound` lemma for the upper direction.
**S117 progress:** Extended Route A''''' to `W = 5` —
`mps_bond_dim_W_eq_5_d_eq_j_plus_1 : (unfolding 5 (j+1) j).rank = 5`
for every `j ≥ 1`, sorry-free. First instance with `R = W` (all `W`
columns retained); first instance using `Matrix.det_of_upperTriangular`
(BlockTriangular id route) rather than `Matrix.det_fin_three`.
**S122 progress:** Extended Route A'''''' to `W = 6` —
`mps_bond_dim_W_eq_6_d_eq_j_plus_1 : (unfolding 6 (j+1) j).rank = 3`
for every `j ≥ 1`, sorry-free. **First orthogonal-corner instance where
the working row set is not `{0, 1, ..., R-1}`**: rows `{1, 2, 3}` of
the `6 × 6` slab are linearly dependent (each has identical support
pattern `(1, 0, 0, 0, 1, 0)` from primes `{7, 11}, {13, 17}, {19, 23}`),
so the construction uses rows `{0, 1, 5}` with `chiP 31` (the smallest
"row-5 prime") providing the third linearly independent row. Returns
to `det_fin_three` (since `R = 3`) plus `upper_bound` (since
`rank_le_width` gives only `rank ≤ 6`). Sets the template for higher-
`W` corners where the first `R` rows of the `W × W` slab are LD.
**S128 progress:** Extended Route A''''''' to `W = 8` —
`mps_bond_dim_W_eq_8_d_eq_j_plus_1 : (unfolding 8 (j+1) j).rank = 5`
for every `j ≥ 1`, sorry-free, via the BlockTriangular route with
permutation `ρ ↦ (2, 0, 1, 3, 4)`, `σ ↦ (0, 1, 2, 6, 4)`, diagonal
primes `{17, 2, 11, 31, 37}`, and `det_of_upperTriangular` +
`Fin.prod_univ_five`. **Seventh unconditional instance.** **W=7
deferred** (no leading-row triangulation; needs `det_fin_seven`
or Laplace expansion).
**S129 progress:** Extended Route A^{(8)} to `W = 12`, **skipping
the structurally-obstructed `W ∈ {7, 9, 10, 11}` corners** —
`mps_bond_dim_W_eq_12_d_eq_j_plus_1 : (unfolding 12 (j+1) j).rank = 5`
for every `j ≥ 1`, sorry-free. Permutation `ρ ↦ (0, 9, 10, 4, 7)`,
`σ ↦ (1, 0, 6, 10, 4)`, diagonal primes `{2, 109, 127, 59, 89}`,
below-diagonal composites `{49, 50, 55, 85, 86, 91, 95, 110, 121, 122}`.
**First instance using FOUR non-leading rows** — extends S122's
single-non-leading-row trick to the maximally non-leading regime.
**Eighth unconditional instance; sixth instance over a wheel `W ≥ 3`;
third instance using `det_of_upperTriangular`.** Skip-rationale: W=9
has clean block-diagonal structure (1+3+3) but neither 3×3 block
admits standalone upper-triangulation; closing W=9 cleanly requires
`Matrix.det_of_blockTriangular` (a new technique not yet used in
the file). W=7, W=10, W=11 hit similar obstructions.
**S137 progress:** Extended Route A^{(9)} to `W = 18` —
`mps_bond_dim_W_eq_18_d_eq_j_plus_1 : (unfolding 18 (j+1) j).rank = 7`
for every `j ≥ 1`, sorry-free. Permutation `ρ ↦ (0, 2, 9, 1, 11, 6, 16)`,
`σ ↦ (1, 6, 16, 10, 12, 0, 4)`, diagonal primes
`{2, 43, 179, 29, 211, 109, 293}`, 21 below-diagonal composites.
**First instance with `R = 7`** — uses `Fin.prod_univ_seven` (mathlib).
**First instance using `norm_num` for primality** — `decide` hits
Lean's `maxRecDepth` for primes ≥ 150 (e.g. 179, 211, 293). Pre-search
showed **W = 14 (also R=7) is structurally obstructed**: rows 2 and 5
of the 14×14 j=1 slab have identical support pattern at the chosen 7
cols, joining `{7, 9, 10, 11}` as needs-`det_of_blockTriangular`.
**Ninth unconditional `mps_bond_dim` instance; seventh instance over a
wheel `W ≥ 3`; fourth instance using `det_of_upperTriangular`.**
**S143 progress:** Extended Route A^{(10)} to `W = 20` —
`mps_bond_dim_W_eq_20_d_eq_j_plus_1 : (unfolding 20 (j+1) j).rank = 9`
for every `j ≥ 1`, sorry-free. **First instance with `R = 9`** and
**first heartbeat-bumped declaration** (`maxHeartbeats 2000000`). Adds
local `prod_univ_nine'` lemma. **Tenth unconditional instance.**
**S144 progress:** Extended Route A^{(11)} to `W = 10`, **refuting
S128/S129's "structurally obstructed" claim** — `mps_bond_dim_W_eq_10_d_eq_j_plus_1
: (unfolding 10 (j+1) j).rank = 5`, sorry-free. Permutation
`ρ ↦ (1, 0, 4, 3, 9)` uses **row 9** (not in {0..4}); the earlier
S128/S129 search restricted to row prefixes {0..R-1}. S144 also
performed a DEFINITIVE DP-based enumeration over `W ∈ [2, 72]`,
`R ≤ 22`: leading-row + dead-col upper-triangulation closes
**EXACTLY** `W ∈ {2, 3, 4, 5, 6, 8, 10, 12, 18, 20}` and is
structurally obstructed for every other W in that range.
**Eleventh unconditional instance.**
**S151 progress:** Pre-search for W=9 BlockTriangular route — Python
script `w9_blocktriangular_search.py` enumerated 32 valid (1+3+3)
block-DIAGONAL decompositions; minimum-new-helpers candidate uses
`ρ ↦ (0, 1, 3, 5, 2, 4, 6), σ ↦ (2, 1, 3, 7, 0, 4, 6)` and adds chiP
helpers `{13, 41, 53, 61}`. Documented S151 LESSON: `Mexp + Matrix.ext
+ fin_cases <;> rw [h_sub]` shortcut hits "motive not type-correct"
errors due to dependent-type proof terms.
**S152 progress:** Closed Route A^{(12)} `W = 9` corner —
`mps_bond_dim_W_eq_9_d_eq_j_plus_1 : (unfolding 9 (j+1) j).rank = 7`,
sorry-free. **FIRST closure of an S128/S129/S144 "block-triangular-
required" wheel; FIRST use of `Matrix.det_fromBlocks_zero₂₁`** —
orthogonal to the previous nine corner closures (all
`det_of_upperTriangular`). Used a NESTED `det_fromBlocks_zero₂₁`
decomposition: outer 1+6 split (1×1 block via `det_fin_one`, gives `det
= 1`) plus inner 3+3 split of the 6×6 (two 3×3 blocks via
`det_fin_three`, each gives `det = -1`). Total `det = 1 · (-1) · (-1)
= 1`. **Crucial design choice:** the 1+6 outer split (vs the S151-
proposed 4+3) avoids any 4×4 det computation, sidestepping the
mathlib-no-`det_fin_four` issue. ~610 Lean lines (49 entry-lemmas +
85 fromBlocks reindex case checks via `rcases ... <;> fin_cases <;>
rfl` + structural assembly). **Eleventh unconditional instance; tenth
over a wheel `W ≥ 3`.** Closed-W set: `{2, 3, 4, 5, 6, 8, 9, 10, 12,
18, 20}`.

**Next action (post-S152):** The S152 nested-`det_fromBlocks_zero₂₁`
template is now the canonical pattern for the remaining "block-
triangular-required" wheels `{7, 11, 14, 15, 16, 24, 30, ...}`. Each
needs a Python pre-search analogous to S151's `w9_blocktriangular_search.py`
to find the row/col permutation that makes one (or more) off-diagonal
block zero. The Lean assembly is reusable: 49+ entry-lemmas in W=20
style + nested fromBlocks reindex + per-block `det_fin_one` /
`det_fin_three` (or `det_fin_two` for 2×2 blocks). Cleanest next
single-session targets: W=7 (R=5) and W=11 (R=11). For W=11 the larger
`R` may need 1+(3+3+3+1) or similar partition. Or pivot to (b) Route C
(mathlib PNT for the low-density regime — leaves saturating half-cut
open). See
`experiments/formalisations/E2_1_mps_bond_dim/mps_bond_dim_notes.md`.

### L2 — E1.5: 0.537-bits invariant
**Statement:** For every modulus m ≥ 2, H(π(x) mod m | π(x-1) mod m) → c
where c = 0.537... bits as x → ∞ (or precisely: the constant is `1 - H_2(1/log 2)` or whatever the exact form is).
**Why this is interesting:** the constant has a closed form related to the prime density 1/log x; making the constant exact in Lean would force precision about its derivation.
**Estimated effort:** 2-3 sessions.
**Save under:** `experiments/formalisations/E1_5_invariant_information_rate/`

### L3 — E2.7: Communication rank exactly 2^{N/2-1}+2
**Statement:** For the prime indicator on N bits with balanced 2-party split,
the communication matrix has rank exactly 2^{N/2-1} + 2.
**Why this:** verified empirically to N=20; the proof exists informally in
`novel/determinantal_complexity.md`. Formalising it is the next step.
**Estimated effort:** 2-3 sessions.
**Save under:** `experiments/formalisations/E2_7_comm_rank_plus_2/`

### L4 — E5.1: BPSW correctness ⇒ PRIMES in TC⁰
**Statement:** A formal Lean construction of the TC⁰ circuit family.
**Why this:** ties together Hesse-Allender-Barrington (scalar powering in TC⁰),
Mereghetti-Palano (2×2 matrix powering in TC⁰), Jacobi symbol in TC⁰. The
formalisation has educational value plus locks in the conditional claim.
**Estimated effort:** 4-5 sessions (substantial circuit-class library work).
**Save under:** `experiments/formalisations/E5_1_bpsw_tc0/`

### L5 — E5.8: Brandt 4-obstruction argument
**Statement:** Formally state and prove that Brandt's TRAVERSE technique
relies on 1-Kt-randomness of Chaitin Ω prefixes (Lemma 2 in Brandt 2024)
and that this property has no analog for any fixed total Boolean function.
**Why this:** the structural argument is the project's freshest closure;
formalising it would be the first formal verification of a Brandt-MKtP
extension claim.
**Estimated effort:** 5-8 sessions (Kt complexity in Lean is non-trivial).
**Save under:** `experiments/formalisations/E5_8_brandt_obstruction/`

### L6 — Hölder simplification of normalised Ramanujan sum (S191 byproduct)
**Statement:** For squarefree q and any integer n with `d := gcd(q, n)`,
```
   mu(q) · c_q(n) / phi(q)  =  mu(d) / phi(q/d).
```
**Why this:** Single-step character-theoretic identity, fully visible
from definitions. Used pointwise in S191's `T_Q` construction and in
S168's main-term simplification. Estimated 1 session of Lean. Mathlib
has `Nat.ArithmeticFunction.moebius`, `Nat.totient`, and the
Ramanujan-sum definition follows from a simple `Finset` over coprime
units. **Bonus value:** isolating the identity enables a future Lean
proof of the *S191 wheel-Mertens identity*: at primorial Q = W,
`M_W(coprime n) = W/φ(W)`.
**Estimated effort:** 1 session.
**Save under:** `experiments/formalisations/L6_holder_normalised_ramanujan/`

---

## §4. Negative-Shape Edges — Family-Level Closures

The strongest project edges are **family-level closures**: E7.10 closes the
AKS family, E7.11 closes the convergence-acceleration family, E5.8 closes
the Brandt family. Each is a single argument that closes dozens of
individual approaches. Producing more of these is high-leverage.

### N1 — Tensor-network compression family — **BUILT (S77)**
**Hypothesis:** No tensor network of polylog bond dimension can represent
χ_P (for any reshape, any tensor-train ordering, any matrix product
structure). E2.1 is the W=primorial special case; E1.9 is the φ(x,a) 2D case;
E6.3 is the wavelet/DCT case. Conjecture a unified statement: "any
factorisation of χ_P with bond dimension polylog(N) requires
exponentially-many factors."
**Why novel:** the individual cases are each known; the family-level
statement is not.
**First step:** State the conjecture precisely. Test it on 3-4 new tensor
factorisations (CP, Tucker, hierarchical Tucker, MERA).
**Outcome (S77):** Family-level closure verified empirically across **22
(W, d) pairs** spanning W ∈ {2,3,4,5,6,7}, d ∈ {4,6,8,10,12}. Five
classical ansätze (MPS, Hierarchical Tucker root-children, PEPS 2D
boundary, Tensor Ring, CP-rank Kruskal LB) are reduced to E2.1's
unfolding-rank lower bound — half-cut bond dim is **identical** across
all five and equals `min(W^j, φ(W)·W^(d-j-1)+1)` (exact saturation in
21/22 cases; one finite-size deficit of 1 at W=5, d=4). Tucker and MERA
close by orthogonal mechanisms (mode-slice independence; parameter
count). Asymptotic ratio `bond_dim / sqrt(N) → φ(W)/W` (the **Mertens
product**) — verified at d=12 for W∈{2,3,4} and at d=8 for W=6: 0.515
vs 0.500; 0.668 vs 0.667; 0.501 vs 0.500; 0.334 vs 0.333. **Net new
content**: (i) E2.1 lifted from single-ansatz to family-level scope —
unfolding-rank bound is the *universal* mechanism; (ii) Mertens product
identified as the universal asymptotic compression ratio across the
family; (iii) the closure subsumes prior single-decomposition
CLOSED_PATHS rows (lines 171, 185, 517, 518, 600). No polylog opening.
**Status:** BUILT. CLOSED_PATHS row added (S77). E2.1 annotated in
EDGES.md.
**Save under:** `experiments/constructions/tensor_compression_family_closure/`

### N2 — Linear-functional family on zero-sums
**Hypothesis:** Any linear functional `Σ w(γ) · h(γ, x)` with `Σ |w(γ)| ≥ c·N(T)`
that approximates ψ(x) has tail error `Ω(sqrt(x))`. E7.11 is the
acceleration-family special case; E3.7 / E3.5 are special instances.
**Why novel:** the tail bound is folklore but a clean theorem statement
covering ALL linear functionals would be publishable.
**First step:** Write the precise statement. Cite Iwaniec-Kowalski Ch. 5
for closest published bound.

### N3 — Reduction-to-MKtP family
**Hypothesis:** Any Brandt-style technique that uses 1-Kt-randomness of
oracle prefixes cannot be applied to fixed total Boolean functions.
E5.8 is the special case for π(x) mod 2.
**Why novel:** this is a meta-theorem about a *class* of complexity
arguments, not just one.
**First step:** Catalogue Brandt-class arguments in the literature
(Brandt 2024, plus any predecessors / extensions). Identify the common
structural ingredient. State the hypothesis precisely.

---

## §5. Synthesis Targets (publishable-paper grade)

### S1 — "Four structural barriers to polylog π(x)" (renamed post-S116)
**Content:** Unify E7.6 + E7.10 + E5.8 + E7.11 + **E7.14** into a single
negative-result paper. Each is a family-level closure of a major
technique. Together they cover sieve-pebbling, AKS modulus-twist,
Brandt-meta, convergence-acceleration, and **Maynard multidim sieve**.
Audience: analytic number theorists + complexity theorists. Likely
venue: arxiv preprint targeting cstheory.
**Save under:** `novel/four_structural_barriers.md` (draft) →
`literature/preprint_four_barriers.md` (final).
**Estimated length:** 18-30 pages.

**S116 update (2026-04-27):** Maynard sieve closed as the fifth
explicit barrier (or fourth, depending on whether E7.6 sieve-pebbling
is counted as a separate family from E7.14 Maynard sieve-weight).
This synthesis now has the natural structure: (i) AKS-family
[E7.10, E5.3]; (ii) Brandt MKtP [E5.8]; (iii) convergence-
acceleration [E7.11]; (iv) Sieve-route [E7.6 ∧ E7.14].

### S2 — "Pseudorandomness of π(x) mod m: a 35-measure battery"
**Content:** Catalogue the 35 measures with definitions, methodology, and
results. Cite project sessions. Frame as "the largest reported empirical
pseudorandomness battery on a single natural number-theoretic function."
**Save under:** `novel/pseudorandomness_battery_paper.md`.

### S3 — "Information-computation gap for π(x): O(log x) information,
O(x^{2/3}) computation"
**Content:** Formalise the project's `novel/info_computation_gap.md`
into a paper-grade statement. Cross-reference E1.1, E1.3, E2.6, E2.7.

---

## How to use this file

When you start a session:

1. Skim the section that matches your time budget:
   - 30-60 minutes: pick a §1 composition or §3 small Lean L_i.
   - 1-2 hours: pick a §2 frame-shift or §3 medium Lean.
   - 2+ hours and clear head: pick a §4 negative-shape or §5 synthesis.

2. **Mark your pick in `RESEARCH_AGENDA.md`** as in-progress. If the
   challenge has multi-session structure, register it as an arc.

3. Build, test, file. Cite edge IDs.

4. **Update this file at session end:** if you completed a challenge,
   move it to `archive/sessions/sessionNN_*.md` and add a closure note
   here. If you produced a *new* challenge through your work, add it
   to the appropriate section.

5. Honest failure: if your selected challenge collapsed to an existing
   closure, file it in CLOSED_PATHS but ALSO note in the challenge entry
   why. The next agent can adjust the challenge or remove it.

---

## What this file is NOT

- It is not TODO.md (housekeeping + recurring tasks live there).
- It is not RESEARCH_AGENDA.md (long-horizon multi-session arcs live there).
- It is not EDGES.md (verified mathematical facts live there).
- It is not CLOSED_PATHS.md (closed approaches live there).

It is a **menu of concrete novelty targets**. Updated when a challenge is
completed or when a new one is identified mid-session.
