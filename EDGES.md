# EDGES.md — Mathematical Cracks in the sqrt(x) Wall

Vulnerability-research-style catalogue of every **real** mathematical edge,
identity, structural deviation, or partial signal that AI agents have
surfaced over 49 sessions / 690+ closed paths. Listed individually with
exact provenance, then chained into hypothetical attack scenarios.

Conventions:
- **Edge** = a verified mathematical fact, identity, or measurable statistical
  deviation from what a "fully random" object would give. Even if it has
  already failed in isolation, it is recorded here because its value in a
  *chain* is what matters.
- **EVS** (Edge-Value Score, qualitative): H = high (could be a primary
  primitive), M = medium (useful in conjunction), L = low (probably noise but
  recorded for completeness).
- Sources cited as `file:line` or `Sxx` for session number; consult
  `status/CLOSED_PATHS.md` for individual closure entries.

---

## 0. The frame the chain has to fit through

Every chain must address all four targets simultaneously, because a chain that
misses any one of them lands in a known closed path:

T1. **Smooth target (~50% of bits):** R^{-1}(n) gives the top ~60% of the
    bit-length of p(n) in O(polylog) (carry-propagation boundary, S40,
    `novel/carry_propagation_boundary.md`). Any chain that only delivers
    these bits is a duplicate of R^{-1}(n).
T2. **Oscillatory bits (~40%):** the remaining bits encode the GUE-random
    contribution of ~x^{1/2} zeta zeros (`novel/info_computation_gap.md`,
    S36 line 212 of CLOSED_PATHS).
T3. **Avoid the three failure modes:**
    C = "needs primes to compute primes",
    E = "reduces to the explicit formula / O(sqrt x) zero sum",
    I = "smooth approximation loses the oscillatory bits"
    (`novel/failure_taxonomy.md`).
T4. **Avoid the Natural Proofs barrier when reasoning about
    impossibility** — every "natural" structural test we tried gives the
    same answer for chi_P and for a uniformly random Boolean function
    (`novel/pseudorandomness_of_pi.md`, 24 measures).

If a chain ends in "and now compute pi(x;q,a) for q < 13" → C-circular.
If a chain ends in "and now sum sqrt(x) zeta zeros" → E-equivalent.
If a chain ends in "round R(x) + smooth correction" → I-info-loss.

T5. **Three-pillars constraint (E7.7):** any new informationally-complete
    encoding of pi(x) must either (a) be a 4th informationally-complete
    primitive distinct from prime positions / zeta zeros / floor values
    (FOCUS-6), or (b) admit polylog conversion to one of the three. 15+
    candidate intermediate-quantity families have been individually
    closed against this constraint (S15, S16).

T6. **CRT-mod-m information rate (E1.5):** every modulus m << pi(X)
    gives `H(pi(x) mod m | pi(x-1) mod m) = h_2(pi(X)/X) + O(1/pi(X))`,
    independent of m (S69 closed-form sharpening of S12's "0.537 bits";
    that constant is the X = 10^4 special case). Combining moduli adds
    linearly, never multiplicatively — CRT-style chains cannot win.
    The asymptotic per-step rate goes to 0 as X to infinity, since
    pi(X)/X to 0.

---

## 1. Information-content edges  (these are why the problem might still be open)

### E1.1 — E(x) = pi(x) - Li(x) carries only O(log x) bits        [EVS H]

`novel/info_computation_gap.md`; S36;
`experiments/wildcard/local_global_reconstruction_results.md`.

For x = 10^100, |E(x)| ~ 10^48, so E uses ~160 bits — under 0.15% of the
polylog budget (log2 x)^2 ~ 110,000. **There is NO information-theoretic
barrier to polylog.** The barrier is purely computational: extracting
160 bits from a sum of ~10^50 GUE-random terms with massive cancellation.

> Why this is an "edge": the Razborov-Rudich pseudorandomness picture is
> output-side. On the *input* side, the answer is small. This is the
> primary asymmetry we want to exploit.

### E1.2 — Adjacent-x autocorrelation of the error is 0.996      [EVS M]

`experiments/wildcard/local_global_reconstruction_results.md` Exp 3;
S36; `archive/sessions/session36_synthesis` style block.

`Corr(E(x), E(x+1)) = 0.996` and `std(E(x+1)-E(x)) = 0.326`.  Helpless for
*single-point* p(n) computation, but **very useful for batch / sliding /
oracle-amortised algorithms**: a single anchor pi(x_0) plus 0.326-std
local increments determines pi at neighbouring points.

### E1.3 — Per-bit difficulty of p(n) is a sharp 4-bit sigmoid    [EVS H]

`novel/carry_propagation_boundary.md`; S40
`experiments/wildcard/carry_propagation_boundary.py`.

Bit positions 0-7 (MSB) match round(R^{-1}(n)) at ≥95%. Positions 8-9
are 80-90%. From bit 10 onward, agreement is 50% (coin flip). Width of
the "hard zone" is *only ~4 bits* and located at ~60% from MSB. For
N-bit p(n), the hard region is bits 0.6·N .. N (asymptotic).

> Why this is an edge: the hard region has a sharp boundary, not a long
> tail. A method that produces *any* extra bit beyond R^{-1}(n)'s ceiling
> is enough — you don't have to break the whole barrier, only the first
> bit past the cliff.

**Bounded-Kt refinement (S105, C3):** Under the 3-bit-per-op stack VM
of `experiments/constructions/brandt_mktp/` (L_MAX = 12, T_MAX = 4096),
the per-bit family `{s_J^(N) : J = 0..N-1}` for N ∈ {4..7} exhibits
its bounded-Kt cut at `J*(N) := ⌈log₂( π(2^N - 1) + 1 )⌉ ≈ N - log₂ N`,
**not** at the 0.5N boundary above. For J ≥ J*(N), `π_J ≡ 0` and
Kt_b ∈ {5..8}; for J < J*(N) — including E1.3's "easy zone"
J ∈ [⌈0.5N⌉, J*) — Kt_b saturates at INF = 61. Stack-VM-class programs
resolve only the *trivially-zero* boundary; the smooth/oscillatory
transition E1.3 detects via Fourier weight is invisible to bounded-Kt
at this resolution. Detecting it would require arithmetic primitives
in the VM (LOG2, LI_APPROX, the R^{-1} kernel) — proposed C3.a. See
`experiments/constructions/brandt_per_bit/`.

### E1.5 — pi(x) mod m saturates at h_2(pi(X)/X) per step                        [EVS H]

`proven/complexity.md` §"Circuit Complexity"; S12;
`proven/circuit_size_barrier.md`;
**closed-form sharpening:** `experiments/information_theory/pi_mod_2k_saturation/pi_mod_2k_saturation_results.md` (S69, F2);
**generalisation to bisection components and joint state:**
`experiments/constructions/g_q_bisection_invariant/` (S70, C1) — same
mechanism (binary increment + conditional independence of indicator from
state) applies to *any* monotone Ω-stratified summatory; verified for
A(x) and C_3(x), and for the joint state `(A mod q, C_3 mod q)` which
satisfies `H(g_q(x) | g_q(x-1)) = H_3(1 - rho_A, rho_pi, rho_C3) +
O(1/pi(X))` in regime q^2 << pi(X).

**Sharpened statement.** For every modulus m and every X with
`m <= pi(X)/100`,
```
   H( pi(x) mod m | pi(x-1) mod m, x in [1, X] )
       =  h_2( pi(X)/X )  +  O( 1 / pi(X) )
```
where `h_2(p) = -p log_2 p - (1-p) log_2 (1-p)` is the binary entropy.
The "0.537 bits" of S12 is specifically the value at X = 10^4 (where
`h_2(pi(10^4)/10^4) = h_2(0.1229) = 0.5376`); the constant tracks prime
density and **decays as X grows**:

| X | pi(X)/X | h_2(pi(X)/X) |
|---|---------|---------------|
| 10^4 | 0.1229  | 0.5376  |
| 10^5 | 0.0959  | 0.4559  |
| 10^6 | 0.0785  | 0.3969  |
| 10^7 | 0.0665  | 0.3526  |

Verified to 7-decimal precision at X = 10^7 across moduli m in {2..1024}
(powers of 2 plus cross-checks). For m approaching pi(X), the conditional
entropy collapses below the closed form (the state pi(x-1) mod m starts
encoding x itself).

**Asymptotic implication.** As X to infinity, `pi(X)/X to 0` and hence
the per-step rate `h_2(pi(X)/X) to 0`. The information rate of pi(x) mod
m **vanishes** in the limit, regardless of m.

> Why this is an edge: the constancy across m (in regime m << pi(X))
> proves that **CRT reconstruction cannot win** — knowing pi(x) mod m
> for any m gives exactly h_2(pi(X)/X) bits of new info per step, so
> combining k moduli scales linearly (not multiplicatively) and offers
> no compression. It also says pi(x) mod 2 is *exactly as hard* as
> pi(x) mod m for any other m. The mechanism is closed-form: since
> `pi(x) - pi(x-1) = 1[x prime] in {0,1}`, the conditional entropy is
> `h_2(P[x prime]) = h_2(pi(X)/X)` whenever the prime indicator is
> conditionally independent of the state — which holds for m << pi(X).

### E1.6 — pi(x) mod 2 bisects into two statistically-independent bits  [EVS M]

S55; `experiments/sieve/pi_mod_q_fixed/liouville_modular_structure_results.md`.
**q-stable extension:** `experiments/constructions/g_q_bisection_invariant/` (S70, C1).

Via the verified identity `pi(x) = (x - L(x))/2 - C_3(x)`,

```
   pi(x) mod 2  =  A(x) mod 2  XOR  C_3(x) mod 2
```

where `A(x) = (x - L(x))/2`. **Both bits are uniform on {0,1}** (H(A) =
H(C_3) = 1.00000 bits) and **statistically independent** (mutual info
1.94 × 10⁻⁵ bits, verified on x ∈ [1, 2 × 10⁶]).

**S70 q-stable extension (composition C1):** the marginal independence
extends from q = 2 to all prime q in {3, 5, 7, 11, 13}, with mutual
information `I( A(x) mod q ; C_3(x) mod q ) <= 1.2 * 10^-4` bits on
x ∈ [1, 2 * 10⁶] (worst case q = 13). Combined with the per-component
S70 finding that A and C_3 each carry asymptotically 1 bit/step
(`rho_A, rho_C3 -> 1/2`) while their integer difference π carries 0
bits/step (`rho_pi -> 0`), the bisection has a *destructive-interference*
character: two high-rate channels cancel into one low-rate channel.

> Why this is an edge: the parity question is split into two equal halves
> with no shared structure to exploit. A polylog algorithm for *one* side
> still needs the other. Crucially: if (and only if) some future technique
> handles A mod 2 cheaply, then C_3 mod 2 is the entire residual problem
> (not "a fraction of the residual"). The split is informationally clean.
> The q-stable extension means the "no shared structure" property is not
> a q = 2 artefact — it persists across all prime moduli probed.

### E1.7 — delta(n) long-range correlation is INDIRECT (AR(7) + Hurst crossover) [EVS M]

S20, S36; `experiments/information_theory/kt_complexity/SYNTHESIS.md`,
`experiments/information_theory/kt_complexity/kt_deep_analysis_results.md`.

The 1/f^{1.69} spectrum of delta(n) (E3-side characterisation, see also
`novel/delta_spectrum.md`) initially looks like genuine long-memory. PACF
analysis at N = 100 000 corrects this:

* **PACF drops to 0.056 at lag 2** despite ACF(200) = 0.84.
* BIC model-selects **AR(7)** as direct memory order; PACF ~ k^{-1.33}
  (alpha > 1, summable).
* DFA Hurst exponent **H = 1.31 with crossover at scale ~572**: H_small =
  1.41 (scales < 572), H_large = 1.19 (scales > 572) — two regimes.
* Total Kt(1..N) ~ 5.58 N + 0.023 N log N — **EXTENSIVE** (linear),
  negligible log correction.
* Transfer entropy: **TE(n -> delta) = 0.013 bits, TE(delta -> delta) =
  0.051 bits.** Past delta values are 4x more informative than n itself;
  n mod 30 contributes only 0.003 bits.

> Why this is an edge: the apparent long memory is a *short-range AR
> process viewed at distance*. Indirect correlation is **not** the
> signature of an exploitable hidden state — only ~7 lags of direct
> dependence, then independence. Combined with TE(n -> delta) ≈ 0,
> this rules out *any* recurrence-based shortcut from n to delta(n)
> (a recurrence of order > 7 would have non-zero PACF at higher lags).

### E1.10 — Gap-shuffled-zeros: the only valid null for prime-frequency probes of zeta-zero D(s)  [EVS shape]

S49; `experiments/analytic/zeta_structure/large_n_correlations/large_n_battery_results.md`;
`archive/sessions/session49_focus4_large_n_zeta.md`.

For any test that Fourier-decomposes the pair-correlation residual
D(s) = R_2_emp(s) − R_2_GUE(s) at frequencies tied to primes
(e.g. f_p = log(p)/L for the Bogomolny-Keating arithmetic
correction), the **uniform-random-frequency null is wrong**: it
produces false-positive p < 0.05 detections even on data with no
arithmetic content beyond local GUE statistics.

The *correct* null is the **gap-shuffled** sequence: take the
unfolded zeros, compute the gaps, randomly permute them, regenerate
the sequence by cumulative summing the shuffled gaps. Gap-shuffling
preserves the local GUE pair correlation (the gap distribution) but
destroys long-range arithmetic structure beyond nearest-neighbour.

S49 measurement at N=8000:

```
Statistic                  Zeta      Gap-shuffled null    z-score
Pearson(D, BK_template)   +0.111     +0.487 ± 0.035       −10.85σ
Phase coherence ⟨cos(φ-π)⟩ +0.544     +0.624 ± 0.036       −2.20σ
```

Both statistics on real zeta are **below** the gap-shuffled null —
zeta is *more* GUE-like than the matched null, the opposite of what
a real BK arithmetic signature would produce.

> Why this is an edge: shape-edge, methodological. Any future zeta-
> zero structural test that uses Fourier-prime amplitude or phase
> coherence MUST calibrate against gap-shuffled null, not uniform-
> random frequencies. Without this calibration, one would falsely
> 'detect' BK at p < 0.05 in any GUE-like sequence.  This narrows
> the failure modes of all future Fourier-prime probes (Chain D
> extensions, FOCUS-4-style retests at higher N, conditional-test
> n-correlation extensions in S57's Conrey-Snaith framework).

### E1.9 — phi(x,a) 2D rank: 22nd pseudorandomness measure       [EVS L]

S65; `experiments/wildcard/phi_2d_lowrank.py` and `_results.md`;
`novel/pseudorandomness_of_pi.md` measure #22.

Viewing the Meissel `phi(x, a)` function as a 2D matrix `M[i, j] =
phi(x_i, a_j)` under four framings (raw, Mertens-residual,
normalised, column-difference), the relative SVD spectrum decays as
`exp(-0.33 * k)` — apparent exponential low-rank — but `||M||_F` grows
linearly in x, so

```
   rank required for integer-precision (+/- 0.5) recovery
        scales LINEARLY in K
        (12 at K=18, ~35 at K=60)
```

not polylog. The relative compressibility cancels exactly against the
absolute scale, the same pattern that breaks DCT/wavelet sparsity at
asymptotic x (E6.3 caveat, wavelet `N^{0.75}` scaling).

> Why this is an edge: adds a clean *combinatorial* pseudorandomness
> measure to a battery previously dominated by algebraic / spectral /
> communication tests. Confirms that even when chi_P is replaced by
> the smoother sieve auxiliary phi, the linear-rank wall persists.

### E1.8 — Spectral reconstruction of delta needs 82% of all modes [EVS M]

S36; `experiments/information_theory/kt_complexity/spectral_algebraic_structure_results.md`.

For exact recovery (RMSE < 1) of delta(n) at N = 100 000, the Fourier
basis needs **41 182 of 50 001 modes (~82%)**; max error < 1 requires
**ALL modes**. The spectrum is a smooth continuum — no discrete spectral
lines, no algebraic relations among coefficients (PSLQ on top peaks: zero).
Spectral peaks weakly correlate with zeta zeros (4/10 top peaks match at
5% error, well within chance given 1000 candidate zeros).

> Why this is an edge: it is the **frequency-domain analogue of E2.7**
> (communication rank deficiency = +2). The "+2" of E2.7 corresponds to
> the ~18% of modes that are NOT needed; the remaining 82% is irreducible.
> Closes the "Fourier sparsity in delta" attack route at scale.

### E1.4 — N/2 universality                                        [EVS M]

`novel/pseudorandomness_of_pi.md` table; `novel/approx_degree_prime.md`;
S28 synthesis. Five independent measures land on the same threshold:

| Measure                             | Boundary       | Source |
|-------------------------------------|----------------|--------|
| Approximate degree of chi_P at eps=0.49 | ceil(N/2)  | S28 |
| Communication matrix rank           | 2^{N/2-1}+2    | S17 |
| PTF (polynomial-threshold) degree   | N/2 (LP-verified, N=4..12) | S35 |
| LFSR length over GF(p) for delta(n) | N/2            | S24,S26 |
| Per-bit R-correlation crossover     | ~bit N/2       | S28 |
| Real tensor rank of chi_P (3-way)   | ~d^{1.5}=2^{N/2} | S31 |

Empirically: bits 0..N/2 are oscillatory, bits N/2..N are smooth.
**This is not a barrier per se; it's a bisection of the problem into a
solved half (smooth) and an unsolved half (oscillatory).**

**S71 scope refinement (C5):** the N/2 boundary holds tightly for the
**parity-of-Omega family** `{chi_P, f_mu_pos, f_lam_pos}` (where adeg
and PTF degree match `ceil(N/2)` at every tested N up to 10), but
NOT for all natural NT Boolean functions: `f_sqfree` falls below
(adeg = 4 < 5 at N = 10 because density 6/pi^2 is far from the
worst-case 1/2), `f_sqfree3` exceeds (adeg = 4 > 3 at N = 6, 6 > 5
at N = 10), and density-matched PRF falls below. The universality is
real but smaller in scope than the original framing suggested. See
`experiments/constructions/n_over_2_universality_class/` and
EDGES.md E2.7 cross-reference.

---

## 2. Algebraic edges  (small but non-zero structural deviations from random)

### E2.1 — TT (MPS) bond-dim identity at every primorial cut       [EVS H]

`novel/mps_bond_dimension.md`; S41/S42; `experiments/wildcard/mps_bond_dimension*.md`.

For the prime indicator chi_P : [1, W^d] → {0,1} reshaped in base W:

```
   rank M^{(j)} = min( W^j , phi(W) * W^{d-j-1} + 1 )
```

at *every* cut 1 ≤ j < d, *every* primorial W ≥ 2. Empirically saturated to
strict equality for W ∈ {2, 6, 30, 210}, d up to 20.

The +1 row is exactly row-0 (the trivial "all integers > W are coprime to
W constraint"). The rest of the rank reduction is **exactly the Mertens
product prod_{p ≤ W}(1 - 1/p)** in disguise.

> Why this is an edge: it is the **only known closed-form identity of any
> kind** for the entanglement-theoretic compressibility of chi_P. The
> savings is only constant-factor phi(W)/W, so on its own it doesn't break
> sqrt(x), but it gives an exact algebraic handle on the wheel-W sieve in
> tensor-network language.

**S74 free-probability refinement (C2):** the squared-singular-value
distribution of `M^(j)` decomposes as **MP(c = φ(W)/W) bulk + O(N^{0.42})
spike outliers + finite structural peak**. The MP-bulk rate equals
`φ(W)/W = ∏_{p≤W}(1 − 1/p)` (the Mertens product). After projecting out
the spike band (`k* ∝ R^{0.85}` on W=2 d=14..22), the bulk free cumulants
match the MP standardized identity `κ_r = c^{r-1}` within 5–10% across
W ∈ {2, 3, 5, 6, 30}. Matched-active iid Bernoulli baseline matches MP
from drop_1; chi_P needs `O(R^{0.85})` drops — the difference is the new
content. **Algorithmic implication:** any spectral compression of `M^(j)`
faithful at the second-moment level (i.e., reproducing κ_2) needs rank
Ω(N^{0.42}). Polynomial-in-N barrier recovered from a free-probabilistic
angle independent of the explicit-formula machinery.
See `experiments/constructions/free_cumulants_chi_p/`.

**S82 spike-eigenvector identification (C2 sub-arc):** the S74 spike
band IS the **residue-class character subspace at small odd primes
coprime to W**. Empirically, the top `k_*` singular vectors of `M^(j)`
(after the leading mean) decompose into per-prime sectors of dimension
exactly `phi(p)` indexed by primes `p ≤ P*(N)` not dividing `W`, plus
cross-sectors of dimension `phi(p_1)·phi(p_2)`. Conductor of each
spike is `2 · p` (or `2 · p_1 · p_2`). Verified at (W=2, d ∈ {14, 16,
18, 20}) and (W=6, d=6); matched-Bernoulli baseline shows centered
residue energy ≤ `5×10⁻⁴` while chi_P spikes ≥ `10⁻²` (≥ 20× ratio).
**Sharpened algorithmic implication:** the polynomial-in-N spike
content IS exactly the small-modulus residue-class biases
`pi(N; q, a)` for `q ≤ Q*(N) ≈ N^{0.21}`, so the C2 spectral compression
barrier is the same object as the E1.5 / T6 saturation barrier viewed
from the spectral side — a textbook **C-circular** collapse.
See `experiments/constructions/spike_eigenvectors_chi_p/`.

**S83 / S98 / S99 / S106 / S107 Lean formalisation status (L1):** in
`experiments/formalisations/E2_1_mps_bond_dim/MPSBondDim/MPSBondDim/Basic.lean`,
`mps_bond_dim` (the main theorem), `upper_bound`, `rank_le_min_dim`,
`row_support_coprime`, `live_columns_count`, and `lower_bound` are all
closed `sorry`-free over Lean 4 + mathlib `v4.30.0-rc2`. The single
remaining obligation is the existential `exists_invertible_submatrix`
isolating the prime-density content. Net: the *rank-bound* logic of
E2.1 is fully formalised; only the prime-existence exhibit remains.
**S98 closed the corner case `(W = 2, j = 1)` unconditionally** via
two new sorry-free theorems
(`exists_invertible_submatrix_W_eq_2_j_eq_1` and
`mps_bond_dim_W_eq_2_j_eq_1 : (unfolding 2 d 1).rank = 2`) using
Bertrand's postulate alone — the first formally verified instance of
`mps_bond_dim` for any concrete `(W, j)`. **S99 closed the orthogonal
corner case `(W = 2, d = j + 1)` without even Bertrand** via
`exists_invertible_submatrix_W_eq_2_d_eq_j_plus_1` and
`mps_bond_dim_W_eq_2_d_eq_j_plus_1 : (unfolding 2 (j+1) j).rank = 2`,
using only `Nat.prime_two` and `Nat.prime_three`. The two corners
together cover the **entire `(j, d - j)` boundary** (i.e. the full
`R = 2` regime) of the parameter grid for `W = 2`. **S106 extended
the orthogonal-corner argument to `W = 3`** —
`exists_invertible_submatrix_W_eq_3_d_eq_j_plus_1` and
`mps_bond_dim_W_eq_3_d_eq_j_plus_1 : (unfolding 3 (j+1) j).rank = 3`
for every `j ≥ 1`, sorry-free, using `Matrix.det_fin_three` and
the explicit primes `2, 3, 5, 7`. **First unconditional `mps_bond_dim`
instance over a wheel `W ≥ 3`** — three corners are now formally
verified. **S107 extended the orthogonal-corner argument to `W = 4`** —
`exists_invertible_submatrix_W_eq_4_d_eq_j_plus_1` and
`mps_bond_dim_W_eq_4_d_eq_j_plus_1 : (unfolding 4 (j+1) j).rank = 3`
for every `j ≥ 1`, sorry-free, using `Matrix.det_fin_three` and the
explicit primes `2, 3, 5, 7, 11`. **The first orthogonal-corner instance
where `rank_le_width` is not tight** — column `3` is `chiP` at multiples
of `4` (all zeros), so `R = 3 = φ(4)·4^0 + 1 < 4` and the upper bound
must cite the general `upper_bound` lemma rather than the trivial column
count. Four corners are now formally verified. See `mps_bond_dim_notes.md`
for the routes still open on the general case.

**S77 family-level scope refinement (N1):** the unfolding-rank lower
bound `≥ φ(W)·W^(d-j-1)+1` is the *universal* half-cut bond dim across
classical polynomial-spatial-locality tensor-network ansätze:
**MPS / Tensor-Train, Hierarchical Tucker (root-children bond dim),
Tensor Ring, PEPS (2D-reshape boundary), and CP-rank Kruskal lower
bound**. Verified at 22 (W, d) pairs spanning W ∈ {2, 3, 4, 5, 6, 7},
d ∈ {4, 6, 8, 10, 12}: half-cut bond dim is **identical** across all
five ansätze in 21/22 cases, with a single finite-size deficit of 1 at
(W=5, d=4, N=625). Tucker and MERA close by orthogonal mechanisms
(mode-slice independence and parameter-count, respectively) — the
conclusion (no polylog ansatz exists) is uniform. The unifying
mechanism is the unfolding-rank lower bound ITSELF; E2.1 is no longer
an MPS-only statement. See
`experiments/constructions/tensor_compression_family_closure/`.

### E2.2 — Liouville/parity identity                                [EVS M (was H, downgraded S55)]

S46, S55; `experiments/proposals/proposal25_liouville_parity_check.py`
and `experiments/sieve/pi_mod_q_fixed/liouville_modular_structure.py`.

```
   pi(x) = (x - L(x)) / 2 - C_3(x)
```

where `L(x) = sum_{n<=x} (-1)^{Omega(n)}` is the Liouville summatory and
`C_3(x) = #{n <= x : Omega(n) odd, Omega(n) >= 3}`. **Verified exactly for
all x in [1, 2 000 000]** in S55 (extending the S46 range).

The identity is structural, not a closed-form shortcut — `C_3(x)` grows
~3·pi(x) and parity-tests pseudorandom against simple proxies.

**Free identity discovered S55:** `L(x) mod 2 = x mod 2` *trivially*
(since L is a sum of x terms each in {+1,-1}, parity = count parity).
So L(x) mod 2 carries ZERO non-trivial prime information.  The S46
quote "polylog L(x) mod 2 → polylog pi(x) mod 2" is vacuously satisfied
without yielding a new bit.

The actual missing primitives via E2.2 are:

* `L(x) mod 4` (equivalently `A(x) mod 2 = (x - L(x))/2 mod 2`), AND
* `C_3(x) mod 2`,

with `pi(x) mod 2 = A(x) mod 2 XOR C_3(x) mod 2` (verified bit-exact on
x ∈ [1, 2 000 000]).

S55 measured both components on x up to 2 × 10⁶ and found:

* A(x) mod 2 has block-entropy 7.9999/8 bits, AC-max[1..30] = 0.0010,
  FFT z = 5.5 (no spectral line), GF(2) LFSR length = N/2 — *more*
  pseudorandom than pi(x) mod 2 itself.
* C_3(x) mod 2 has block-entropy 7.88/8, AC = 0.148 (density bias only),
  GF(2) LFSR length = N/2.
* Mutual information I(A mod 2; C_3 mod 2) ≈ 2 × 10⁻⁵ bits — the two
  are statistically independent.
* 11 cheap arithmetic proxies (M(x), Q(x), is_square(x), x mod 4, ...)
  all show |corr| < 0.002 with pi(x) mod 2; best XOR-fusion of any
  4-subset achieves 0.4951 agreement (worse than chance).

> Edge value DOWNGRADED to M (from H) by S55. The identity is exact and
> beautiful, but it splits pi(x) mod 2 into two components which are
> *each* structurally pseudorandom (LFSR, AC, FFT, MI) and *mutually
> independent*. So a polylog algorithm for either bit alone is not in
> any current attack surface, and the "MISSING" primitives are now
> known to be at-least-as-hard as pi(x) mod 2 itself under the 26-measure
> pseudorandomness battery (`novel/pseudorandomness_of_pi.md`).

### E2.3 — Ono / van Ittersum / Craig partition identity            [EVS L]

`literature/state_of_art_2026.md` §4.1, S17, S23, S36 closures.

n is prime iff certain MacMahon partition equations hold:
  `(n^2 - 3n + 2) sigma_1(n) - 8 M_2(n) = 0` (and infinitely many similar).

Verified, but uses sigma_1 (divisor sum, factoring) → mode C. Adds nothing
to the chain unless `M_k(n)` admits a polylog evaluation that doesn't
factor through divisors, which it doesn't (Ono himself: O(n^{3/2}) at best).

### E2.4 — Cloitre 2025 effective analytic recurrence               [EVS L]

`literature/state_of_art_2026.md` §4.2; arXiv:2508.02690.

```
   p(n+1) = ceil( ( -1 + zeta(2 p_n) * prod_{j=1..n} (1 - 1/p_j^{2 p_n}) )^{-1/(2 p_n)} )
```

EXACT for all n ≥ 1, no sieving. But (a) the product is over primes
(circular), (b) the zeta evaluation needed is at extreme precision —
trades sieving for equally expensive zeta computation. Recorded because
it is the ONLY known fully-analytic exact recurrence for p(n).

### E2.5 — pi(x) is a multilinear polynomial of degree N             [EVS M]

S15; `novel/determinantal_complexity.md`;
`experiments/circuit_complexity/det_perm_encoding.py` and `_results.md`.

For N-bit input, pi(x) is a multilinear polynomial of degree exactly N
over R / Z (verified empirically to N=10). Determinantal representations
(N x N matrices with affine-linear entries whose determinant equals pi)
**have been numerically constructed for N=2,3,4** but NOT for N=5,6 —
optimization landscape is highly non-convex.

> Edge value: GapL ⊂ NC^2. If an N x N determinantal representation
> exists for all N, pi(x) is in NC^2. Generic degree-N polynomials in
> N variables don't fit (parameter count drops below 2^N at N >= 10), so
> any fit must come from special prime structure.

### E2.6 — Higher-order finite differences of pi(x) are tiny         [EVS M]

`novel/prime_indicator_anf_structure.md`.

For the integer multilinear representation of pi(x), order-7 differences
have max = 15 at N=14, vs random expectation 2^7 = 128. **pi(x) has very
low higher-order bit interactions** — it is *almost* an additive function
of bit positions, with the residue concentrated in the top-degree term
(coefficient grows as ~e^{0.38 N} ~ x^{0.55}).

> Edge value: the additive part is exactly R(x). The "residue" is exactly
> the explicit-formula correction. This is a third witness (after E1.4 and
> E2.1) of the smooth/oscillatory bisection.

**S48-fresh refinement:** A *fourth* independent witness of the 2^{N/2}
coefficient-magnitude wall comes from forward-differences of pi at 0:

```
   log_2 |Delta^k pi(0)|  ≈  0.995 k - 3.17    (k = 1..200)
```

Cleanly linear with slope 1, confirming pi(x) is **order-1 entire** in
the forward-difference sense (sharpest non-smooth behaviour possible).
S48-fresh's parallel multilinear-extension test on (Z/2)^n gave
max-coefficient `~ 2^{n/2}` and nonzero-coefficient ratio `~ 0.71` — both
matching E2.1's TT bond-dim 2^{N/2} and E2.7's communication rank
2^{N/2-1}+2. Three independent operations (TT reshape, 2-party split,
forward differences) on three different objects (chi_P, comm matrix,
pi at 0) all reproduce the same Theta(2^{N/2}) wall. Not a separate
edge — just the cleanest empirical confirmation that the wall is
*operation-invariant*.
See `archive/sessions/session48_fresh_perspective_5_wildcards.md`.

### E2.7 — Rank deficiency of communication matrix is *exactly* +2  [EVS M]

S17, S19; `experiments/circuit_complexity/min_circuit_size.py`.

For the prime indicator on N bits with balanced 2-party split, the
communication matrix has rank exactly `2^{N/2-1} + 2` (verified N=2..20).
Universal formula extends to unbalanced splits: rank `= 2^{min(k,N-k)-1}+2`
when min(k, N-k) ≥ 3.

The "+2" is **the entire signal** distinguishing chi_P from random
(which would saturate at 2^{N/2}). It comes from the
*monotonicity-against-x* structure of chi_P (only 2 is even).

> Edge value: minimal but non-trivial. The "+2" is the only piece of
> exploitable algebraic structure in the matrix, and it is exactly what
> makes pi(x) special.

**S71 indicator analogue + structural unification (C5):**
For the **indicator** chi_P (vs the counting function pi(x) used here),
the analogous balanced split gives `rank(M_chi_P) = 2^{N/2-1} + 1`
exactly (verified N=6,8,10,12,14). The "+1 vs +2" is not a discrepancy:
M_pi adds a constant column from the cumulative integral, M_chi_P does
not. Both are special cases of the **column-zero density bound**:
`rank(M_f^{balanced}) <= (1 - rho_f) * 2^{N/2}`, where rho_f is the
fraction of lower-half indices b for which `f(a*2^{N/2} + b) = 0`
for every row a. For chi_P, even-b columns are zero (x even => x
non-prime for x > 2), so rho_chi_P = 1/2; for sqfree and mu_pos,
b mod 4 == 0 columns are zero (4 | x => x non-squarefree), so
rho = 1/4 and rank = 3 * 2^{N/2-1} exactly; for lam_pos, no column is
forced to zero (lambda(4) = +1), so rho = 0 and rank = 2^{N/2}. This
is a **structural unification** with E2.8's 25-35% tensor rank
deficiency: same column-zero principle, different complexity measure.
See `experiments/constructions/n_over_2_universality_class/`.

### E2.8 — chi_P is 25-35% simpler than random in tensor rank       [EVS L]

S31; `experiments/circuit_complexity/tensor_rank_robust.py`.

3-way tensor rank of chi_P scales as ~d^{1.5} = 2^{N/2} but is
consistently 25-35% below a random function at the same density:

```
  N=6  (4x4x4)   rank = 5,    random ~ 7
  N=9  (8x8x8)   rank ≤ 19,   random ~ 29
  N=12 (16x16x16) rank ≈ 67,  random > 100
```

> Edge value: confirms E2.7's "small but non-zero" structural deviation
> at a different complexity measure. Constant-factor savings, not polylog.

### E2.10 — Free identity: L(x) mod 2 = x mod 2                     [EVS L]

S55; `experiments/sieve/pi_mod_q_fixed/liouville_modular_structure.py`.

Since `L(x) = sum_{n=1..x} lambda(n)` with `lambda(n) ∈ {-1, +1}`,
parity of L equals parity of the count: `L(x) mod 2 = x mod 2`. Bit-exact
on all 2 × 10⁶ tested x.

> Edge value: low, but it's a **trap warning**. Any chain whose endpoint
> is "compute L(x) mod 2 in polylog" is vacuously satisfied (it's just
> x mod 2) without yielding a single non-trivial bit. The actual missing
> Liouville-route primitives are L(x) mod 4 (= A(x) mod 2) and C_3(x)
> mod 2 — see E1.6.

### E2.11 — Higher-order finite differences of f(x)=pi(x)-R(x) GROW like white noise [EVS M]

S29; `experiments/algebraic/identity_search/RESULTS_SUMMARY.md`,
`experiments/algebraic/identity_search/wz_definite_sum.py`.

Iterated finite-difference operator Delta^k applied to f(x) = pi(x) - R(x)
on x ∈ [2, 100 000]:

* RMS of Delta^k f **grows monotonically** with k, ratio between
  successive RMS values converging to **exactly 2.0** as k increases.
  This is the signature of an i.i.d. Gaussian-like sequence: each Delta^k
  doubles variance.
* No finite-order linear difference operator annihilates f.
* Hankel matrix of K(x) = f(x)·x·log(x) has **full rank 250/250**
  (incompressible — no Padé / continued-fraction representation).
* Wilf-Zeilberger definite-sum certificate fit: total variance explained = 0.

> Why this is an edge: contrasts with **E2.6** (low higher-order ANF
> diffs of pi(x)'s integer multilinear representation, max=15 at order
> 7). E2.6 says pi(x) the *function-on-bits* is ANF-sparse; E2.11 says
> the residual after subtracting the smooth part R(x) is *fully white*
> in the difference operator. So the smooth part R(x) absorbs all the
> ANF-sparse piece exactly. **This pins R(x) as the precise smooth/random
> separator** — there is no third "intermediate" piece to peel off.

**S67 extension (FOCUS-2 closure):** the E2.11 finite-difference
signature was applied as a *pre-test* against six concrete candidate
fourth-encoding intermediate quantities (T_1=Σ{logΓ}, T_2=ΣH_n,
T_3=ΣH_n², Ψ(x, log^c x) for c∈{2,3,4}, Σσ_2, Σσ_3, Q(x), Σφ(n)).
Result: all 9 close as mode I in two flavours — WHITE-A
(residual std ~ √x, ratio → 2.0 like π−R) for 7 candidates, WHITE-B
(residual at f64 precision, function entirely smooth) for T_2, T_3.
**The pre-test rules out fourth-encoding candidates ~150x faster
than S64's ρ-correlation method** (0.2s vs 30s/candidate) and is a
strictly stronger filter, since it tests structural finite-difference
equivalence rather than pairwise correlation. See
`experiments/algebraic/fourth_encoding_search/e211_pretest_focus2_results.md`.

### E2.12 — Chebyshev psi link captures 91% of f(x) variance       [EVS L]

S29; `experiments/algebraic/identity_search/algebraic_relations.py`.

Define `g(x) = psi(x) - x` (Chebyshev's prime-power summatory minus its
PNT mean). Then `g(x)/log(x)` and `f(x) = pi(x) - R(x)` satisfy

```
   Pearson correlation r(g(x)/log(x), f(x)) = 0.996,
   variance of f explained by g/log = 91%.
```

This is the **partial-summation identity** between pi and psi made
quantitative — but psi(x) costs O(x) to compute exactly, so the link
yields no algorithmic shortcut. Recorded because it is the *only*
non-trivial correlation found across 7 algebraic-relation experiments
(Bernoulli numbers, zeta(2..7), Dirichlet L(1,chi), Ramanujan tau,
LLL minimal polynomials, ODE/integral-equation fits all returned
zero correlation or x-specific spurious fits).

> Why this is an edge: it singles out g(x)/log(x) as the *unique*
> universal predictor of f(x) in the project's tested basis. Any
> fundamentally new identity must either improve on 91% (very hard,
> the residual is GUE-noise) or beat the O(x) computation cost of psi —
> the latter is the actual barrier.

### E2.9 — F_2 Fourier weight is BELOW random for d ≥ 2             [EVS L]

S31; `experiments/circuit_complexity/f2_correlation_profile.py`.

After accounting for the parity term (W(0) + W(1) = 47-68% of weight
from "primes are odd"), the remaining Fourier weight at degree d ≥ 2 is
*more pseudorandom than random* — z-scores -2 to -30 (i.e., chi_P sits
2-30 standard deviations *below* the average random function in higher-
degree weight). Combined with E2.7, the algebraic story is:
"trivially-explained on the surface, anti-structured underneath."

### E2.13 — Gowers U^k norm of chi_P matches Hardy-Littlewood {0,1}^k-cube singular series; universal across multiplicative-indicator family  [EVS M]

S85, refined S93, S101; ATTACK_VECTORS §D6, NOVELTY_CHALLENGES §D6.b/§D6.c;
`experiments/information_theory/gowers_uk_chi_p.py`,
`experiments/information_theory/hardy_littlewood_box.py`,
`experiments/information_theory/gowers_uk_chi_p_results.md`,
`experiments/information_theory/lambda_vs_chi_p_uk/lambda_vs_chi_p_uk.py`,
`experiments/information_theory/lambda_vs_chi_p_uk/lambda_vs_chi_p_uk_results.md`,
`experiments/information_theory/mu_weighted_chi_p_uk/mu_weighted_chi_p_uk.py`,
`experiments/information_theory/mu_weighted_chi_p_uk/mu_weighted_chi_p_uk_results.md`;
CLOSED_PATHS row at session 85.

Define the normalized Gowers uniformity norm `Q^k(f) := ||f||_{U^k}^{2^k} / E[f]^{2^k}`
on `Z/NZ`. Hardy-Littlewood for the {0,1}^k-cube prime k-tuple
configuration predicts

   `Q^k(chi_P) -> S_k    as N -> infinity,`
   `S_k = product_p alpha_p(k) / (1 - 1/p)^{2^k},`
   `alpha_p(k) = #{(x, h_1, ..., h_k) in (Z/p)^{k+1} : x + epsilon . h != 0 for all epsilon in {0,1}^k} / p^{k+1}`.

Numerically (direct enumeration in `(Z/p)^{k+1}`, product truncated):

  `S_2 = 2.300938`   (truncation error <= 1e-5 at P_max = 100)
  `S_3 = 54.116088`  (truncation error ~ 0.2% at P_max = 47)

Empirically (cyclic Z/N, FFT for U^2, recursion for U^3):

  `Q^2(chi_P)`: 2.103 (N=2^10) → 2.132 (N=2^12) → 2.146 (N=2^14)
                → 2.149 (N=2^16) → 2.153 (N=2^18). Trends monotonically to S_2.
  `Q^3(chi_P)`: 35.5 +- 0.1 across N ∈ [2^10, 2^13]; far above the
                random-Bernoulli baseline (which → 1) but slow to
                approach S_3 = 54 due to compounded finite-N corrections.

The W-trick `chi_{W,b}(n) = chi_P(W*n + b)` for `(b, W) = 1`:

  `S^(W)_k = product_{p does not divide W} alpha_p(k)/(1-1/p)^{2^k}`.

  `S^(W=210)_2 = 1.0023`, empirical `Q^2(chi_{W=210, b=1}) = 1.003`
  at N=2^12 — match within 0.1%, confirming W-tricked chi_P is
  Gowers-uniform at U^2.

**Universality across log-weighting (S93, D6.b):** the same statement
`Q^k(f) → S_k` holds for **f = Λ (von Mangoldt log-weighted)** as well
as for f = chi_P. Side-by-side measurement at N up to 2^17:

  `Q^2(chi_P)`: 2.103 (N=2^10) → 2.149 (N=2^16) → 2.151 (N=2^17). Diff
  `Q^2(Λ)`:    2.107 (N=2^10) → 2.155 (N=2^16) → 2.155 (N=2^17). Diff
  `|Q^2(Λ) − Q^2(chi_P)| / Q^2(chi_P)` ≤ 0.4 % across all N tested.

After W-trick at W = 210, b = 1, N = 2^12: Q^2(chi_W) = Q^2(Λ_W) =
1.0029 to four decimals. The W-trick erases not just the bulk HL
structure but ALSO the residual log-weight-vs-bare discrepancy. The
small persistent positive offset Q^k(Λ) > Q^k(chi_P) (≤ 0.4 % at U^2,
≤ 2.5 % at U^3) is identifiable with the prime-power weight in Λ
(Λ counts n = p^k for k ≥ 2), a contribution that decays as
π(√N) / π(N) → 0.

Refined statement: **S_k is the universal Gowers fingerprint of the
{0,1}^k-cube prime correlation, independent of the weighting scheme
applied to detect it (bare indicator chi_P or log-weighted Λ).**

**Family-level extension (S101, D6.c):** the same U^2 / W-trick
instrument applied to a panel of 11 multiplicatively-defined
arithmetic indicators on Z/NZ produces a sharp dichotomy between
"HL-structured" and "Gowers-uniform" indicators, governed by support
density:

  | indicator class                                      | density   | Q^2_inf at large N |
  |------------------------------------------------------|-----------|--------------------|
  | chi_P (primes, Omega = 1, sqfree)                    | -> 0      | 2.301 = S_2        |
  | mu^2 (squarefree)                                    | 6/pi^2    | **1.0384** stable  |
  | 1[mu = +1] (sqfree, even omega)                      | 3/pi^2    | 1.0384             |
  | 1[mu = -1] (sqfree, odd omega)                       | 3/pi^2    | 1.0384             |
  | 1[Omega = 2] (semi-primes)                           | -> 0      | 1.05 -> 1.12 slow  |
  | 1[lambda = +1] (Liouville-positive)                  | -> 1/2    | **1.0000** stable  |
  | 1[lambda = -1] (Liouville-negative)                  | -> 1/2    | 1.0000             |
  | mu (signed Mobius)                                   | mean -> 0 | Q^2_norm -> 1.0    |
  | lambda (signed Liouville)                            | mean -> 0 | Q^2_norm -> 1.0    |

Two structural facts:

  1. **Mobius cancellation** propagates from signed mu, lambda to
     **density-1/2 indicators** 1[lambda = +/-1] (which inherit
     Gowers-uniformity Q^2 -> 1.0000), but NOT to **asymmetric-density
     indicators** 1[mu = +/-1] / mu^2 (which retain Q^2 ~ 1.0384).
     New empirical constant `S_2^{sqfree} ≈ 1.0384` (the
     squarefree-mod-p^2 Hardy-Littlewood product).

  2. **W-trick at W = 210 collapses the entire family** to
     Q^2 ≈ 1 ± 0.005: chi_P → 1.004, sqfree → 1.000, 1[mu = ±1] →
     1.001, 1[Omega = 2] → 1.001. The W-trick is therefore not a
     prime-indicator-specific tool but a family-level HL-killing
     operation on multiplicatively-defined arithmetic indicators.

The S87 "Liouville-positive is Gowers-uniform" result is now
structurally explained as a *density-1/2 phenomenon*, not a
Liouville-specific one.

> Why this is an edge: **first project measure where chi_P
> deviates from random by a CLOSED-FORM, predicted, constant
> multiplicative factor (the Hardy-Littlewood singular series)**, in
> contrast to the 35+ prior pseudorandomness measures that all
> returned "indistinguishable from random". The deviation is exactly
> the prime k-tuple correlation captured by HL — there is no
> additional structure beyond HL. Combined with E2.9 (degree-d
> Fourier weight BELOW random for d >= 2): chi_P's higher-order
> structure is "exactly Hardy-Littlewood, no more". This narrows the
> "what extra structure could chi_P have" search space.

> Why this is **not** an algorithmic opening: extracting Q^k(chi_P)
> empirically requires summing over (x, h) ∈ (Z/N)^{k+1}, cost
> Theta(N^2 log N) at U^2 — already non-polylog. The structural
> identity Q^k(chi_P) = S_k provides only the same information that
> HL already gives, in compressed form.

> Cross-domain reference: Gowers norms (Gowers 2001), Green-Tao
> arXiv:math/0606088 (linear equations in primes), Green-Tao
> arXiv:0807.1736 (Mobius / nilsequences), Green-Tao-Ziegler
> arXiv:1009.3998 (U^{s+1} inverse theorem), Hardy-Littlewood 1923.

### E2.14 — Anderson Lyapunov gamma(E) of chi_P-driven Schrödinger operator: deviation cascade matches W-trick  [EVS M]

S88; ATTACK_VECTORS §C4;
`experiments/dynamical/anderson_localisation_chi_p/anderson_localisation_chi_p.py`,
`experiments/dynamical/anderson_localisation_chi_p/parity_control.py`,
`experiments/dynamical/anderson_localisation_chi_p/wtrick_control.py`,
`experiments/dynamical/anderson_localisation_chi_p/anderson_localisation_chi_p_results.md`;
CLOSED_PATHS row at session 88.

Define the discrete 1D Schrödinger operator
`H psi(n) = -psi(n+1) - psi(n-1) + V(n) psi(n)`, transfer matrix
`T_n(E) = [[V(n) - E, -1], [1, 0]] in SL(2, R)`, and Lyapunov exponent
`gamma(E) = lim_N (1/N) log ||T_N ... T_1||`. For sparse `V in {0, 1}`
at density `rho` Pastur-Figotin gives `gamma(E) ~ rho(1-rho)/(8 sin^2(k))`
inside the band `E = -2 cos(k) in (-2, 2)`.

Empirically (`N = 10^5`, 50 random seeds, 51 energies):

   Baseline (Bernoulli at matched density rho = pi(N)/N):
       max |z|(gamma_prime - gamma_bern) = 88.5 sigma at E = 0.108

The deviation reduces monotonically when the random control sieves
out small-modulus residue classes (W-trick):

   W = 2  (parity matched)         max |z| = 32.7   at E = 1.088
   W = 6  (mod 2, 3 matched)       max |z| = 11.9   at E = 0.010 (N=2*10^5)
   W = 30                           max |z| = 6.3
   W = 210                          max |z| = 6.1
   W = 2310                         max |z| = 4.0

The peak at `E = 1.088 ~ -2 cos(2 pi / 3) = +1` is the **mod-3
resonance**: the free transfer matrix rotates by `2 pi / 3` per step,
maximally coupling to mod-3-periodic potentials. The peak at
`E ~ 0` (`k = pi / 2`) is the parity resonance.

> Why this is an edge: chi_P deviates from random on a *spectral
> global* statistic (Lyapunov exponent of the associated Schrödinger
> operator), in a category orthogonal to all 36 prior pseudorandomness
> measures (which are local correlation / entropy / Fourier / tensor
> rank, plus E2.13 Gowers norms). Its deviation is captured by the
> W-trick to within ~ 4 sigma at W = 2310, giving a SECOND
> independent confirmation of the "chi_P structure = Hardy-Littlewood
> equidistribution mod q" picture (E2.13 was the first). This rules
> out hypothetical "spectral-only" structure of chi_P that would not
> show up in additive-combinatorics measures.

> Why this is **not** an algorithmic opening: gamma(E) extraction
> requires N transfer-matrix multiplications => Theta(N) cost per
> energy, no polylog gain. The cascade pattern matches HL exactly,
> producing no new bit beyond what equidistribution already gives.

> Composes with: E2.13 (Gowers norms -> HL singular series); E1.10
> (chi_P locally indistinguishable from random); E7.x family
> (negative-shape bound on chi_P-as-deviation-from-random).

> Cross-domain references: Aizenman-Warzel 2015 *Random Operators*
> (AMS GSM 168), chs 6-7; Furstenberg-Kifer 1983 Israel J. Math. 46
> (Lyapunov exponents of SL(2,R) products); Pastur-Figotin 1992
> *Spectra of Random and Almost-Periodic Operators*; Green-Tao
> arXiv:math/0606088 (W-trick).

### E2.15 — Algebraic immunity AI_F_2(chi_P) = 2, encoded by (1+x_0)(1+x_1)        [EVS M]

S92; ATTACK_VECTORS §B1;
`experiments/algebraic/algebraic_immunity_chi_p/algebraic_immunity_chi_p.py`,
`experiments/algebraic/algebraic_immunity_chi_p/extract_annihilator.py`,
`experiments/algebraic/algebraic_immunity_chi_p/wtrick_AI.py`,
`experiments/algebraic/algebraic_immunity_chi_p/algebraic_immunity_chi_p_results.md`;
CLOSED_PATHS row at session 92.

The algebraic immunity (Courtois-Meier 2003) of `chi_P : F_2^N -> F_2`
is the smallest degree `d` such that some non-zero polynomial
`g in F_2[x_0, ..., x_{N-1}]` of total degree `<= d` annihilates
chi_P or `1 + chi_P`:

    AI(f) = min{ deg(g) : g != 0 and (g*f == 0 OR g*(1+f) == 0) }.

**Empirical:** `AI(chi_P, N) = 2` for ALL `N in [4, 13]`. The unique
smallest-degree annihilator is the SAME polynomial for every N >= 5:

    g(x) = 1 + x_0 + x_1 + x_{0,1} = (1 + x_0)(1 + x_1).

`g(n) = 1` iff `n ≡ 0 mod 4`. The annihilation is structural:
- composite n ≡ 0 mod 4: g(n) = 1, chi_P(n) = 0, product = 0;
- odd n or n ≡ 1, 2, 3 mod 4: g(n) = 0, product = 0;
- the only prime with bit_0 = 0 is n = 2, and bit_1(2) = 1, so g(2) = 0.

In contrast, AI(matched-density Bernoulli) grows as
`Theta(log_2(1/rho))` (Faugere-Ars 2003 heuristic), reaching 4 at
N = 11 with zero std (every random seed strictly higher than chi_P).

**Liouville+** (1[lambda(n) = +1]) does NOT have this annihilator
(lambda(4) = +1, so 4 is in support); AI(Liouville+) grows linearly,
matching random.

**Mobius!=0** (squarefree indicator) DOES inherit the same annihilator
because no n divisible by 4 is squarefree.

**W-trick correction**: `chi_P_{W,b}(n) := chi_P(W*n + b)` with
gcd(b, W) = 1 removes the small-modulus residue structure.
Empirically:

    W = 1 (no sieve):                   AI(chi_P, N=11) = 2
    W = 2 (parity):                     AI(chi_P, N=11) = 4
    W = 6 (b in {1, 5}):                AI(chi_P, N=11) = 5
    W = 30 (b in {1, 7, 11}):           AI(chi_P, N=11) = 5

vs AI(random matched-density, N=11) = 4 to 5 with zero std at the
W-tricked densities. **The deviation is fully removed by W >= 6.**

> Why this is an edge: chi_P deviates from random on the
> polynomial-method invariant central to the Croot-Lev-Pach / cap-set
> machinery (Tao 2016). The deviation is sharp (zero std at most N >= 8)
> and persistent (constant 2 across all N tested). It is the
> POLYNOMIAL-METHOD encoding of the trivial mod-4 sieve fact, in a
> mathematical category (Boolean polynomial method / algebraic
> cryptanalysis) ORTHOGONAL to E2.13 (additive combinatorics) and
> E2.14 (spectral / Anderson localisation). Three independent
> confirmations of the "chi_P structure = Hardy-Littlewood
> equidistribution mod q" picture in three orthogonal categories.

> Why this is **not** an algorithmic opening: AI extraction is
> exponential time in N (LP over a `2^N x sum_i C(N,i)` F_2 matrix).
> The annihilator is the trivial mod-4 fact, fully captured by the
> W-trick at W >= 6. F_p multilinear ANF degree of chi_P is
> near-saturated for all tested (p, k); slice rank brackets are
> non-informative at p = 2 and match random at p = 3, k >= 3.

> Composes with: E2.13 (Gowers U^k -> HL singular series, S85);
> E2.14 (Anderson Lyapunov, S88); together with E2.15 these form a
> triple of independent confirmations of mod-q sieve structure.

> Cross-domain references: Croot-Lev-Pach arXiv:1605.01506
> "Progression-free sets in Z_4^n are exponentially small";
> Ellenberg-Gijswijt 2017 *Annals* 185, 339; Tao 2016 blog
> (slice rank); Courtois-Meier 2003 Eurocrypt LNCS 2656 (algebraic
> immunity); Faugere-Ars 2003 (heuristic AI scaling); Green-Tao
> arXiv:math/0606088 (W-trick).

### E2.16 — Primes are NOT a translation-invariant DPP / PPP / signed-K / complex-Hermitian-K point process  [EVS M]

S95; ATTACK_VECTORS §D7;
`experiments/constructions/primes_dpp_ppp_fit/primes_dpp_ppp_fit.py`,
`experiments/constructions/primes_dpp_ppp_fit/primes_dpp_ppp_fit_results.md`;
CLOSED_PATHS row at session 95.

A translation-invariant simple point process on `Z` with intensity
`rho` and real-symmetric kernel `K(t) = K(-t)` has 2-point inclusion
probability `R_2(t) = rho^2 - K(t)^2` (DPP) or `rho^2 + K(t)^2` (PPP).
For prime `chi_P`:

  - **Pair level (F1):** `K^2_DPP(t) = rho^2 - R_2(t) < 0` for all 15
    admissible even t in [2, 30] at N = 10^7 (Hardy-Littlewood
    `S(0, t) > 1` ⇒ `R_2 > rho^2`). DPP infeasible.
  - **Pair level (F2):** `K^2_PPP(t) = R_2(t) - rho^2 < 0` for all 14
    odd t > 1 (`R_2 ≈ 0` ⇒ `R_2 < rho^2`). PPP infeasible.
  - **3-point (F3):** restricting to all-even sub-lattice where
    `K^2_PPP > 0`, the permanent prediction
    `R_3^PPP = perm[K]` overshoots the HL singular series prediction
    `R_3^HL = rho^3 S(0, t_1, t_2)` by 7.5–79.2% across 19 admissible
    triples. Maximum gap on 3-AP triples (0, 6, 12), (0, 12, 18),
    (0, 18, 24).
  - **3-point (F4):** real-signed K with magnitudes
    `|K(t)| = rho sqrt(S(0, t) - 1)` and arbitrary signs `s : t → ±1`
    fails at every triple — required cross-term ratio
    `sigma_req ∈ (-0.541, +0.769)` is never `±1`.
  - **3-point (F5):** complex Hermitian K of the same magnitudes with
    free phases `phi : t → [0, 2π)` admits no globally consistent
    assignment; least-squares fit over 200 random starts produces
    best max residual 0.0746 ≫ 0.01 sample-noise floor.

**Mechanism:** the HL singular series `S(0, t_1, …, t_k)` factorises
over PRIMES (`prod_p alpha_p`), with `alpha_p` depending on
`nu_p({0, t_1, …, t_k})` = number of distinct residues mod p.
DPP/PPP correlations factorise over PAIRS. Pairwise admissibility
does NOT imply triple admissibility: e.g., `(0, 4, 14)` is admissible
in every pair mod 3 but inadmissible as a triple (covers all of
`Z/3Z`), giving `R_3^HL = 0` while PPP predicts a non-zero value.

> Why this is an edge: this is the FIRST 3-point structural statement
> ruling out a kernel-factorisation representation of chi_P,
> complementing the 2-point closures E2.13 (Gowers norm), E2.14
> (Anderson Lyapunov), E2.15 (algebraic immunity). The DPP/PPP
> framework is canonical in random-matrix theory and the GUE eigenvalue
> ensemble IS a sine-kernel DPP — the contrast highlights that the
> "GUE pair correlation" of zeta zeros (E3.1, E3.13) does NOT transfer
> to a kernel structure on chi_P itself. Adds a fourth-leg confirmation
> in a new mathematical category (point-process theory).

> Why this is **not** an algorithmic opening: even if a kernel existed,
> `K(0) = rho = 1/log N` gives only `pi(N) ≈ N/log N` (PNT). DPP/PPP
> 2-point evaluation requires `K(t)` for t up to log² N anyway —
> already the singular-series cost.

> Composes with: E2.13 (Gowers U^k); E2.14 (Anderson Lyapunov);
> E2.15 (algebraic immunity); E1.10 / E3.13 (GUE pair correlation
> on zeros, NOT on chi_P).

> Cross-domain references: Hough-Krishnapur-Peres-Virag 2009 *Zeros
> of Gaussian Analytic Functions and Determinantal Point Processes*
> (AMS ULect 51); Soshnikov 2000 "Determinantal random point fields"
> Russ. Math. Surv. 55, 923 (arXiv:math/0002099); Vere-Jones 1997
> (alpha-determinantal generalisation); Hardy-Littlewood 1923
> (singular series).

### E2.17 — Persistent homology of Takens-embedded normalised prime gaps deviates ≥ 5σ from IID Exp(1) AND from gap-marginal-permuted control; W=210 W-trick erases the serial-correlation component  [EVS M]

S96, refined S117 (D2.a, two-component), S124 (D2.a.1, three-
component), and S131 (D2.a.1.i, four-baseline disentanglement);
ATTACK_VECTORS §D2; NOVELTY_CHALLENGES §D2.a, §D2.a.1, §D2.a.1.i;
`experiments/topological/persistent_homology_chi_p/persistent_homology_chi_p.py`,
`experiments/topological/persistent_homology_chi_p/persistent_homology_chi_p_results.md`,
`experiments/topological/persistent_homology_w_trick/persistent_homology_w_trick.py`,
`experiments/topological/persistent_homology_w_trick/persistent_homology_w_trick_results.md`,
`experiments/topological/persistent_homology_w_trick_marginal_b3/persistent_homology_w_trick_marginal_b3.py`,
`experiments/topological/persistent_homology_w_trick_marginal_b3/persistent_homology_w_trick_marginal_b3_results.md`,
`experiments/topological/persistent_homology_w_trick_discrete_b4/persistent_homology_w_trick_discrete_b4.py`,
`experiments/topological/persistent_homology_w_trick_discrete_b4/persistent_homology_w_trick_discrete_b4_results.md`;
CLOSED_PATHS rows at sessions 96, 117, 124, 131.

For a window of M = 2000 consecutive normalised gaps
`x_n = (p_{n+1} - p_n) / log(p_n)` starting near `p ≈ 10^6`, Takens-
embedded at delay 1 in dimension d ∈ {2, 3, 4}, the Vietoris-Rips
persistence diagram (computed via ripser, `thresh = 4`) has total
H_0 persistence `T0` and total H_1 persistence `T1` BOTH lower than
random baselines:

  - **vs B1 (IID Exp(1)):** `T0` z-score = -10.31, `T1` z-score = -4.20
    at d=3, M=2000 (rank 0/20 in K=20 baseline replicates, both
    summaries).
  - **vs B2 (gap-permuted, marginal preserved):** `T0` z = -8.70,
    `T1` z = -11.99 (rank 0/20 both).
  - **Robust:** signature persists at d ∈ {2, 3, 4} (T0 z(B2) ∈
    [-8.7, -5.1]) and at a different window x ≈ 5·10^6 (T0 z(B2) =
    -7.58).
  - **Scaling:** z-scores grow super-linearly with M
    (M=500 ⇒ T0 z(B1) = -4.2; M=4000 ⇒ T0 z(B1) = -17.8). Signal is
    at least linear in window size, not a finite-N artifact.

**Mechanism:** Hardy-Littlewood k-tuple admissibility constrains
consecutive gaps to repeat residue patterns more often than random,
creating geometric self-similarity in the delay-embedding cloud
(small `T0` = clusters merge faster) and suppressing random
"out-and-back" loops in delay space (small `T1`). The B2 control
preserves the gap MARGINAL but destroys serial correlation; the
T0/T1 deficit relative to B2 isolates the *serial-structure*
component.

> Why this is an edge: it is the first quantitative persistent-
> homology measurement on a prime sequence in this project, and the
> signed deviation (T0, T1 *below* random) is genuinely informative —
> primes are *more clustered and less loop-rich* than Poisson on a
> sliding gap window. Joins E2.13 (Gowers U^k), E2.14 (Anderson
> Lyapunov), E2.15 (algebraic immunity), E2.16 (DPP failure) as a
> fifth-leg confirmation of HL serial structure, in a new
> mathematical category (algebraic-topological / metric-space
> geometry).

> Why this is **not** an algorithmic opening: persistence on a
> Vietoris-Rips filtration of an M-point cloud costs O(M^2) distance
> computations and O(M^3) worst-case PH. No closed-form polylog
> evaluator is suggested; PH is a measurement instrument here, not
> an algorithm.

> Composes with: E2.13 (Gowers); E2.14 (Anderson); E2.15 (alg.
> immunity); E2.16 (DPP failure). Detects the same singular-series
> structure but via a metric-topological observable rather than
> Fourier or kernel.

**S117 refinement (D2.a, W-trick decomposition):** filtering primes
to a single residue class `q ≡ b (mod W)` with `W = 210, b ∈ {1, 11,
13}`, then re-Cramér-normalising as `x_n = g_n / (φ(W) log q_n)`
with `φ(210) = 48`, gives a sub-sequence whose Takens-embedded PH
diagram **matches its own gap-marginal-permuted shuffle to within
Gaussian noise**. Concretely (M=1000, d=3, x≈10⁶, K=20, pooled
across the three residues):

  - z(W-tricked PRIMES; B2)_T0 = -1.99 (S96 unconditioned: -7.45)
  - z(W-tricked PRIMES; B2)_T1 = -0.67 (S96 unconditioned: -4.05)

At a different window (x ≈ 5·10⁶) and at d ∈ {2, 4} the same
collapse holds (|z(B2)| ≤ 2.5 across every (M, d, x) cell vs ≥ 4
unconditioned). At the original S96 anchor (M=2000, d=3, x≈10⁶,
single residue b=1) the reduction is 3.0×–5.8× on T0/T1
(z(B2)_T0 = -2.87, z(B2)_T1 = -2.08).

The IID Exp(1) control B1 retains a large deficit (z up to ≈ 12,
with T1 sign-flipping positive at d ∈ {3, 4}, M = 1000) because
the W-tricked gap MARGINAL itself differs from Exp(1) — gaps are
multiples of W = 210, giving a discrete Cramér-normalised spectrum
on a quasi-grid of spacing ≈ 0.318. This marginal-distribution
deviation from Exp(1) is structurally distinct from HL serial
structure and is preserved by construction in B2.

**Decomposition (S117 — two-component):** E2.17 PH deficit on bare
prime gaps = (serial-correlation component) + (marginal-distribution
component).

  - Serial component → killed by W=210 W-trick (z vs B2 collapses).
  - Marginal component → preserved/amplified (gap-quantisation).

This places E2.17 as the **sixth orthogonal leg** of the W-trick-
erases-everything HL fingerprint family alongside E2.13 (Gowers
U^k), E2.14 (Anderson Lyapunov), E2.15 (algebraic immunity),
E2.16 (DPP failure), E2.20 (subword complexity). PH is a metric-
space topological observable on a delay embedding; its serial-
correlation deviation from random is exactly the HL singular series
at W=210.

**S124 refinement (D2.a.1 — three-component decomposition).** The
S117 marginal-distribution component splits further. Adding a third
baseline B3 = IID samples from the **continuous** distribution
matched to the W-tricked empirical marginal (linearly-interpolated
empirical CDF; Devroye 1986 *Non-Uniform Random Variate Generation*
§II.2.1) gives, at M=1000, d=3, x≈10⁶, three residues pooled:

  - z(W-tricked PRIMES; B3)_T0 = −0.05  (B1: −9.64; B2: −1.99)
  - z(W-tricked PRIMES; B3)_T1 = +0.46  (B1: +5.71; B2: −0.74)

Across d ∈ {2, 3, 4}: |z(B3)| ≤ 0.65 on every (d, T0/T1) cell —
the continuous-envelope baseline absorbs the entire S117 residual
to within Gaussian noise. Hence:

  E2.17 deficit = Δ_envelope (≈ 7-9σ on T0; the W-tricked marginal
                  variance ≈ 0.55, support on (0.318, ∞), differs
                  from Exp(1) variance 1)
                + Δ_discreteness (1-3σ on T0, sign: B2 mean > B3
                  mean — discrete-grid permutation has higher T0
                  than continuous-envelope IID)
                + Δ_serial_residual (~1-2σ on T0 post-W-trick;
                  S117's residual gap correlation among primes
                  p > 7 within a residue class).

The discreteness and serial-residual components partially cancel on
the (PRIMES_W vs B3) comparison: B2 has T0 ≈ +5 above B3 (discrete-
grid PH lattice effect raises T0); PRIMES_W has T0 ≈ −5 below B2
(serial correlation lowers T0). They null on B3 → z(B3) ≈ 0 even
though both individual components are 1-3σ. Marginal-envelope
remains the sole dominant component.

This refines E2.17's status as a singular-series fingerprint to a
specifically-decomposed statement: PH detects the HL singular
series **primarily through the marginal CDF shape**, not through
the discrete-grid quantisation or through residual gap-correlation.

**S131 refinement (D2.a.1.i — four-baseline disentanglement).** Adding
a fourth baseline B4 = IID with replacement from the empirical PMF
of the W-tricked normalised gaps separates S124's "discreteness" sub-
component from a previously-unisolated H_0-specific cloud-geometry
artifact. At M=1000, K=20, three residues pooled, x ≈ 10⁶:

  - z(P_W; B4)_T0 = −0.67 / −0.56 / +0.07 across d ∈ {2, 3, 4}.
  - z(P_W; B4)_T1 = −1.23 / −0.51 / +0.45 across d ∈ {2, 3, 4}.
  - (B4 − B3)_T0 = +1.74 / +2.87 / +3.95 — same sign as (B2 − B3)_T0
    = +3.61 / +5.42 / +6.86 across d ∈ {2, 3, 4}. Both discrete-
    support baselines lift T0 above continuous-support B3.
  - (B2 − B4)_T0 = +1.87 / +2.56 / +2.91 across d ∈ {2, 3, 4} —
    strictly positive, monotone in d. **NEW T0 component**: the
    duplicate-compression effect.
  - (B4 − B2)_T1 = −0.11 / +0.20 / −0.98 — small at all d, consistent
    with H_0-specific (not H_1) localisation of the new component.
  - B4 duplicate counts: 366 / 368 / 371 of M=1000 across the three
    runs, matching theory `M(1 − (1−1/M)^M) ≈ 0.368M` to 0.5%.

**Mechanism of the new component:** B4's IID-with-replacement draws
generate ≈ 0.368M duplicate values per replicate. In the Takens delay
embedding, duplicate underlying values produce zero-distance pairs in
the cloud whenever the duplicates align in any d consecutive coordinates
(or, more generally, induce small-distance pairs by partial overlap).
In the H_0 filtration, zero-distance pairs merge at edge weight 0,
contributing **zero** to T0 (bar `(0, 0)` has length 0). The B4 T0
mean thus reflects H_0 persistence over the ≈ 0.632M *unique* delay
points only, against B2's full M points. B2 has each empirical value
exactly once and therefore no zero-distance pairs.

The H_1 (loop) statistic is **insensitive** to a single duplicate
pair because a loop requires ≥ 3 vertices and an interior. Hence
T1(B4) ≈ T1(B2), as observed.

**Four-way decomposition (S131-refined):**

  E2.17 deficit on bare W-tricked primes
   = Δ_envelope        ≈ 7-9σ on T0    [B1 vs B3, S117]
   + Δ_discreteness    ≈ 3-7 mean-gap on T0
                        ≈ 0.5-1.5 mean-gap on T1
                        [B3 vs {B2, B4}; S124 + S131 confirms baseline-
                         independence: BOTH discrete-support baselines
                         exhibit the lift]
   + Δ_duplicate       ≈ 2-3 mean-gap on T0, NULL on T1
                        [B4 vs B2; **NEW S131**; H_0-specific
                         cloud-geometry artifact of IID-with-replacement
                         sampling, NOT a primes-structural feature]
   + Δ_serial_residual ≤ 1σ on T0, ≤ 1.2σ on T1
                        [P_W vs B4; **tightened from S124's 1-2σ**]

Δ_duplicate is technical book-keeping for any future PH-on-empirical-
IID study (any natural number-theoretic indicator with a non-degenerate
empirical PMF will exhibit the same artifact under IID-with-replacement
baselines). Δ_serial_residual after this disentanglement is bounded
above by ≤ 1σ on T0 — the W-tricked primes' residual gap-correlation is
**at the noise floor** of the PH instrument once duplicate-compression
is factored out.

> Cross-domain references: Carlsson 2009 "Topology and data" Bull.
> AMS 46(2); Edelsbrunner-Harer 2010 *Computational Topology: An
> Introduction*; Bauer 2021 "Ripser" J. Appl. Comput. Topol. 5
> (arXiv:1908.02518); Cohen-Steiner-Edelsbrunner-Harer 2007 (PH
> stability); Green-Tao 2008 *Annals* 167 (arXiv:math/0404188) —
> origin of W-trick; Devroye 1986 *Non-Uniform Random Variate
> Generation* Springer §II.2.1 — inverse-transform sampling;
> coupon-collector / birthday-problem theory for the 0.368M
> duplicate-count formula.

### E2.18 — Anderson Lyapunov gamma(E) of Liouville-driven Schrödinger operator (V = lambda in {-1, +1}) is INDISTINGUISHABLE from i.i.d. Rademacher across [-3, 3] **without W-trick**  [EVS M (multiplicative regime, paired with E2.14)]

S100; ATTACK_VECTORS §G1;
`experiments/dynamical/liouville_anderson_lyapunov/liouville_anderson_lyapunov.py`,
`experiments/dynamical/liouville_anderson_lyapunov/analyze_g1.py`,
`experiments/dynamical/liouville_anderson_lyapunov/liouville_anderson_lyapunov_results.md`;
CLOSED_PATHS row at session 100.

The discrete 1D Schrödinger operator
`H psi(n) = -psi(n+1) - psi(n-1) + V(n) psi(n)` with the **centered
multiplicative potential** `V(n) = lambda(n) in {-1, +1}` (Liouville
function) gives a Lyapunov exponent `gamma_lambda(E)` matching an
i.i.d. Rademacher baseline **within seed noise** at all 51 energies
in `E ∈ [-2.95, 2.95]`, across three orders of N magnitude:

```
N         seeds  max |z|   argmax E    chi^2 / K   L^2 rank of lambda
10^5      50     1.78      -0.236      0.63        21 / 50
3*10^5    50     2.16      +0.118      0.49         7 / 50
10^6      100    2.04      -2.006      0.69        41 / 100
```

`max |z|` is **flat in N** (1.8-2.2, well below the 51-energy
Bonferroni threshold z = 3.16) and the argmax energy **wanders**
between runs — statistical, not arithmetic. **No W-trick is needed.**
Pastur-Figotin agreement: `gamma_lambda / gamma_PF = 0.9317 ± 0.32`
inside the band, identical to `gamma_Rademacher / gamma_PF = 0.9309`
to 4 decimals.

Independent two-point Chowla aggregate at N = 10^6:
`Σ_{h=1..16} z_h^2 = 4.77` (vs Rademacher chi^2_16 mean 16, std
sqrt(32) ≈ 5.7). Lambda is *more Rademacher-like than Rademacher*
on this independent off-spectral test.

**Stark contrast with E2.14 (chi_P, S88):** at the same N=10^5 grid,
chi_P had max |z| = 88.5 sigma at E ≈ 0.108 (parity resonance) and
required W = 2310 sieve (Green-Tao W-trick) to reduce to ~ 4 sigma.
Liouville needs no W-trick at all.

> Why this is an edge: it is the **first non-W-tricked spectral
> measurement** of a prime-related sequence in this project's 38-
> measure pseudorandomness battery to land at noise floor without
> any sieving. It pairs with E2.14 to **isolate** the source of
> chi_P's Anderson-Lyapunov deviation as exclusively HL singular-
> series mod-q resonance: the multiplicative-regime companion of
> chi_P (Liouville, density 1/2, centered) carries no such
> resonance and therefore no spectral signal at the scales
> measured. Spectral analogue of Green-Tao's Möbius/nilsequence
> orthogonality theorem (arXiv:0807.1736) and Tao's logarithmic-
> Chowla theorem (arXiv:1509.05422). First project use of Mobius
> orthogonality as a CROSS_DOMAIN_TECHNIQUES entry.

> Why this is **not** an algorithmic opening: the absence of
> spectral structure for lambda **closes** the intended polylog
> opening "if Liouville-side Lyapunov has detectable structure, the
> explicit formula μ-side `Σ Li(x^ρ) μ_ρ` could be partial-summed in
> polylog time." Lambda is featureless, so the explicit-formula
> Möbius side is irreducible at this measurement scale.

> Composes with: E2.14 (chi_P Anderson Lyapunov, direct paired
> contrast); E2.13 (Gowers U^k of chi_P, additive analogue closes
> via W-trick); E1.10, E3.13 (chi_P pseudorandomness battery).

> Cross-domain references: Aizenman-Warzel 2015 *Random Operators*
> AMS GSM 168 chs 6-7; Furstenberg-Kifer 1983 Israel J. Math. 46;
> Pastur-Figotin 1992; **Green-Tao 2012 Annals 175 = arXiv:0807.1736
> (Möbius nilsequence orthogonality);** Sarnak 2010 IAS lectures
> (Möbius randomness conjecture); **Tao 2016 Forum Math Pi 4 e8 =
> arXiv:1509.05422 (logarithmically averaged Chowla);** Tao-
> Teräväinen 2017 arXiv:1708.02610.

### E2.19 — Subword complexity p_w(n) of chi_P deviates from matched-density Bernoulli on a clean monotone W-trick cascade {RAW > ODD > W6 > W30 > W210}; W=210 erases the deviation to ≤ 8.4 sigma at L = 2.4·10^4  [EVS M]

S104; ATTACK_VECTORS §D13;
`experiments/dynamical/subword_complexity_chi_p/subword_complexity_chi_p.py`,
`experiments/dynamical/subword_complexity_chi_p/subword_complexity_chi_p_results.md`,
`experiments/dynamical/subword_complexity_chi_p/results.json`;
CLOSED_PATHS row at session 104.

For a binary infinite word `w` (Wikipedia "Complexity function"; Lind-
Marcus 1995 *An Introduction to Symbolic Dynamics and Coding* CUP;
Cassaigne-Nicolas 2010 "Factor complexity" in CANT vol. 135), the
**subword complexity**

  `p_w(n) := #{distinct length-n factors of w}`

with topological entropy `h_w := lim_n log_2 p_w(n) / n` is the
canonical combinatorial-on-words invariant. Morse-Hedlund 1938
*Amer. J. Math.* 60: `p_w(n) ≤ n` for some `n` iff `w` is ultimately
periodic.

For five chi_P-derived streams at N = 5 · 10^6, n_max = 22, K = 20
permutation and Bernoulli matched-density baselines:

| Stream | density rho | L         | max\|z_perm\| | at n |
|--------|------------:|----------:|--------------:|-----:|
| RAW    | 0.0697 | 4 999 999 | 132.7 | 22 |
| ODD    | 0.1394 | 2 499 999 | 277.1 |  7 |
| W6     | 0.2090 |   833 334 | 120.5 |  8 |
| W30    | 0.2611 |   166 667 |  24.8 | 17 |
| W210   | 0.3053 |    23 810 |   8.4 | 12 |

`RAW` is the bare chi_P(2..N) stream; `ODD` keeps only chi_P(2k+1);
`W{q}` = Green-Tao W-trick at primorial q with residue r = 1.
`p_chi(22) / p_perm_mean(22)` for the same five streams: 0.018, 0.216,
0.806, 0.994, 1.011.

**Sign:** chi_P has FEWER distinct n-grams than matched-density random
on RAW/ODD/W6/W30 (restricted alphabet from divisibility constraints
mod small primes). On W=210 the sign flips at n ∈ [17, 22]: chi_P
slightly above the permutation mean, consistent with finite-L
saturation noise (`h_eff = log_2 L / n = 0.661` at n=22).

**Effective topological entropy match (W=210):** at n=22, h_chi =
0.6581 vs h_perm = 0.6574 — agreement to ≤ 0.001, indistinguishable
from Bernoulli matched-density across n ∈ [1, 22].

The cascade `|z_perm|: 132 → 277 → 120 → 25 → 8` is monotone in the
HL framework: each W-trick step removes the corresponding mod-q
admissibility constraints (primes ≤ p_k coprime to W = primorial(p_k));
the remaining residual at W=210 is `O(1/log N)` consistent with the
HL singular-series tail.

> Why this is an edge: **first quantitative subword-complexity
> measurement of chi_P** in the project (CLOSED_PATHS line 181 was an
> informal placeholder with no scale, no baseline, no W-trick). The
> 38th pseudorandomness measure in the project's battery, and the
> SIXTH orthogonal HL-detection family — first in symbolic-dynamics /
> factor-complexity / topological-entropy category. Joins E2.13
> (Gowers U^k), E2.14 (Anderson Lyapunov), E2.15 (algebraic
> immunity), E2.16 (DPP failure), E2.17 (persistent homology) under
> the same W-trick-erases-everything fingerprint — a sixth-leg
> confirmation that chi_P deviation from random is exactly Hardy-
> Littlewood mod q.

> Why this is **not** an algorithmic opening: computing `p_w(n)` for
> a given window is `O(L log L)` (rolling encode + sort), no polylog
> gain. The W-trick fingerprint just instruments the same HL physics
> in a new mathematical category — new probe, same signal.

> Composes with: E2.13 (Gowers norms — additive-combinatorics
> analogue at the same W-trick limit); E2.14 (Anderson Lyapunov —
> spectral analogue with same W-trick cascade); E2.15 (algebraic
> immunity); E2.16 (DPP failure); E2.17 (PH); E7.7 (three-pillars).

> Cross-domain references: Cassaigne-Nicolas 2010 "Factor complexity"
> in *Combinatorics, Automata and Number Theory* CANT vol. 135
> (Cambridge); Lind-Marcus 1995 *An Introduction to Symbolic Dynamics
> and Coding* (Cambridge); Morse-Hedlund 1938-40 "Symbolic Dynamics"
> *Amer. J. Math.* 60.

### E2.20 — Mahler measure of the prime indicator polynomial f_N(z) = Σ_{n≤N} χ_P(n) z^n is below density-matched Bernoulli baseline by a CONSTANT asymptotic deficit Δ_∞ ≈ −0.307 nat; f_N / z² irreducible over Q[z] for N ∈ {64, 128, 256} (no cyclotomic share)  [EVS M]

S134, mode I, B-grade negative-shape edge.

**Statement.** Let `f_N(z) := Σ_{n=1}^{N} χ_P(n) z^n ∈ Z[z]` and define
`m(f_N) = exp(∫₀¹ log|f_N(e^{2π i θ})| dθ)` (Jensen). Empirically at
`N ∈ {2^{10}, 2^{12}, 2^{14}, 2^{16}, 2^{17}, 2^{18}}`, with M = 2^{18}
or 2^{19} Jensen-FFT sample points and 50–100 random-Bernoulli /
matched-cardinality controls per N:

  `Δ(N) := log m(f_PRIMES,N) − log m(f_BERN(d_N),N)`

is monotone in N and plateaus at **`Δ_∞ = −0.307 ± 0.001 nat`** by
`N = 2^{16}`. z-score `z(N=2^{18})` is `−110σ` against Bernoulli, `−337σ`
against matched-cardinality. The shared exponent fit
`log m ≈ α log N + β` gives `α_PRIMES = 0.4566, α_BERN = 0.4577`
(indistinguishable) but distinct intercepts `β_PRIMES − β_BERN ≈ −0.30`.

**Cyclotomic structure (small N).** Factoring `f_N(z)` over Q[z] via
sympy at `N ∈ {64, 128, 256}`:
- `f_N(z) = z² · g_N(z)` where `g_N` is **irreducible** of degree
  `N − 2 − 1{N is composite} = π(N) − 1` (just enumerated as 59,
  125, 249 respectively).
- **Zero cyclotomic factors.** No `Φ_n(z) | f_N`. The A-grade
  cyclotomic-compressibility hypothesis `m(f_N) = O((log N)^c)` is
  decisively **falsified**.
- Max root modulus drops with N: `|α|_max(64) = 1.1249`,
  `|α|_max(128) = 1.0666`. Roots cluster near the unit circle but
  not on it.

**Why this is an edge.** This is the **6th orthogonal pseudorandomness
measure category** for χ_P, after E2.13 (Gowers — additive-combinatorial),
E2.14 (Anderson Lyapunov — spectral), E2.15 (algebraic immunity —
F_2-algebraic), E2.16 (DPP failure — point-process / determinantal),
E2.17 (persistent homology — metric / topological). The Mahler measure
is a **multiplicative-height / log-Weil-height** invariant from
algebraic number theory and adds a fundamentally distinct probe.

The deficit direction (`< 0`) is **same as E2.17** (PRIMES persistence
T₀, T₁ smaller than Poisson) — primes are *more constrained* than
random, simultaneously in topology AND algebraic height. This is a
coherent under-randomness pattern that does NOT translate to polylog
extraction (`m(f_N) = Θ(√N)` is still polynomial, and `f_N` is
irreducible — no cyclotomic-evaluator gateway).

> Numerical signature: `Δ_∞ ≈ −0.307 nat`, `m_PRIMES / m_BERN ≈ 0.736`.

> Falsifies / constrains: F1 (Lehmer-typical) — primes are NOT
> indistinguishable from random Bernoulli at this metric. F2
> (cyclotomic / A-grade) — `f_N` has zero cyclotomic share.

> Composes with: E2.13 (Gowers); E2.14 (Anderson); E2.15 (algebraic
> immunity); E2.16 (DPP failure); E2.17 (PH); E7.7 (three-pillars). The
> Hardy-Littlewood singular series (in major-arc form) and the
> determinantal-failure edge E2.16 are the most plausible mechanisms;
> reduction to either would close E2.20 as duplicate.

> Open: derive `Δ_∞ = −0.307` from H-L singular series via Cesàro-summed
> major-arc Jensen contributions. If matches, E2.20 reduces to E2.13/
> E2.16. If not, E2.20 is an algebraic-height invariant orthogonal to HL.

> Cross-domain references: Smyth 2008 "The Mahler measure of algebraic
> numbers: a survey" (CUP); Lehmer 1933 *Ann. Math.* 34; Boyd 1981
> *Canad. Math. Bull.* 24; Dobrowolski 1979 lower bound; Wikipedia:
> Mahler measure / Lehmer's conjecture.

> Files: `experiments/algebraic/mahler_measure_chi_p/` — `.py`,
> `_results.md`, `results/main.json`, `results/scale.json`,
> `results/cyclo_N{64,128,256}.json`, `results/roots_N{64,128}.json`.
> See `archive/sessions/session134_d10_mahler_chi_p.md`.

### E2.21 — L^∞ norm of the prime-indicator polynomial f_N(z) = Σ χ_P(n) z^n at every rational point a/q matches the Hardy-Littlewood density factor μ(q)²/φ(q)·√π(N); ‖f_N‖_∞ = (1+o(1))·√π(N), attained at z = −1 (parity major arc)  [EVS M]

S138, mode I (HL singular series imprint), B-grade negative-shape edge.

**Statement.** Let `f_N(z) := Σ_{n=2}^{N} χ_P(n) z^n` (length N+1, weight
M = π(N)). For every rational `a/q` with `(a, q) = 1` and `q ≤ Q`,

  `R(a/q; primes) := |f_N(e^{2π i a/q})| / √π(N)
                   = √π(N) · μ(q)² / φ(q) + O(√(2 log N))`

verified empirically at `N ∈ {2^{10}, 2^{12}, 2^{14}, 2^{16}, 2^{18}}`
for `1 ≤ q ≤ 58` (limited only by FFT resolution at small N). The
fluctuation `O(√(2 log N))` is Salem-Zygmund / Vinogradov-typical noise.

**Consequence.** `‖f_N‖_∞ = (1 + o(1)) · √π(N) = (1 + o(1)) · M^{1/2}`,
attained at `z = e^{i π} = −1` (the q=2 parity major arc). Equivalently:
`R_N := ‖f_N‖_∞ / √π(N) → 1` (because the dominant rational point
contributes √π(N) · μ²(2)/φ(2) = √π(N)).

**Why this is an edge.** The L^∞-norm endpoint (p = ∞) of the L^p
Fourier-side characterisation of χ_P. This is **the OPPOSITE extreme of
Erdős/Newman flatness**: a flat polynomial would have `R = O(1)` (e.g.,
Rudin-Shapiro `R = √2`); χ_P is far from flat — the parity peak alone
forces `R = √π(N) → ∞`. Complements:
- E2.20 (Mahler = log integral / geometric mean): D27 makes the parity
  major-arc contribution to the Jensen integral explicit, partially
  accounting for Δ_∞ ≈ −0.307 nat via `(1/2π)∫ log|1 − e^{iθ}| dθ`-type
  contributions concentrated near z = −1.
- E2.13 (Gowers U^k matches HL singular series): same μ²/φ structure
  governs both U^k and L^∞ at major arcs.
- E1.5 (explicit formula): the major-arc value
  `Σ_{p ≤ N} e^{2π i p a/q} ≈ μ(q)/φ(q) · π(N)` is the HL/Vinogradov
  exponential-sum identity.

**Q_max scan (minor-arc residual).** Excluding all rationals up to
Q_max, the residual R(prime; minor) ≈ √π(N)/φ(q*) where q* is the
smallest squarefree integer > Q_max. As Q_max → ∞, R(prime; minor)
matches matched-density Bernoulli noise floor (z ≤ 2.5σ at Q_max=32,
N=2^{16}).

**Negative-shape consequence for algorithms.** Any sup-norm-based
compression of f_N (e.g., L^∞-Chebyshev approximation) cannot beat
√π(N) ≈ √(N/log N) storage. The parity peak alone forces this lower
bound, ruling out a "flat-polynomial" representation of χ_P as a
polylog primality witness.

> Numerical signature at N = 2^{16} (π(N) = 6542):
> `√π(N) = 80.88`; `R(prime, q=2) = 80.85`; HL prediction = 80.88.
> Empirical = HL ± 0.04 (≈ 0.05% relative). Across q ≤ 58: max
> relative deviation 6%, mean 1.5%.

> Falsifies / constrains: (a) Erdős-extremal flatness for prime
> indicator. (b) Salem-Zygmund-typicality of bare ‖f_N‖_∞.

> Composes with: E1.5, E1.10, E2.13, E2.20.

> Open: D27.a — Vinogradov minor-arc supremum at N = 2^{18}, 10^4
> random rationals a/q with q ∈ (N^{2/5}, N^{3/5}); does the empirical
> max R(a/q) match `C · N^{1/4} (log N)^A` for explicit C, A?
> D27.b — Liouville L^∞ at major arcs (Möbius/nilsequence prediction
> R^λ(a/q) = o(1) for ALL q including q=2). D27.c — twin-prime Newman
> polynomial L^∞ uses HL twin-prime singular series — different q-pattern.

> Cross-domain references: Erdős 1957 *Mich. Math. J.* 4; Newman 1976
> *Proc. AMS* 51; Kahane 1985 *Some Random Series of Functions* CUP §6
> (Salem-Zygmund); Bourgain 1988 *Acta Math.* 161 (Λ(p) sets);
> Balister-Bollobás-Morris-Sahasrabudhe 2020 *Annals* 192, 977
> (resolution of Erdős for ±1).

> Files: `experiments/analytic/newman_linfty_chi_p/` — `.py`,
> `_results.md`, `results.json`. See
> `archive/sessions/session138_d27_newman_linfty_chi_p.md`.

---

## 3. Analytic edges  (zeta-zero and explicit-formula side)

### E3.1 — Connes-Consani-Moscovici operator: 50 zeros from primes < 13  [EVS L (was H, downgraded S53)]

arXiv:2511.22755 (Nov 2025); `literature/state_of_art_2026.md` §2.5b.
**Polylog claim closed S53.** See `archive/sessions/session53_connes_operator_scaling.md`.

A spectral-triple operator built using ONLY primes p ≤ 13 produces the
first 50 nontrivial zeta zeros with errors ranging from
**2.5 x 10^{-55} (first zero) to ~10^{-3} (50th zero)**. Regularized
determinant converges to the Riemann Xi function.

> Edge value DOWNGRADED to L (from H) by S53. Three independent arguments
> close the polylog claim: (i) rank-one perturbation has B parameters,
> Cauchy interlacing limits eigenvalue placement to ~B substantial
> shifts; (ii) diagonalization of N×N truncation costs O(N^3) = O(K^3)
> for K eigenvalues — strictly worse than O(x^{1/2+eps}) summation when
> K = sqrt(x); (iii) CCM's published B=6 → K=50 with err range
> 10^{-55..-3} implies geometric per-zero error growth ~11.3x, yielding
> K_accurate(B=6) ≈ 53 even at face value, with no multi-B data point
> in the literature establishing super-linear scaling. CCM's
> construction is a *one-shot fit*, not an algorithm. See
> `experiments/analytic/connes_operator/` and CLOSED_PATHS line 696.

### E3.2 — Explicit formula needs **only ~1200 zeros at x with 10^8
digits** (DEBUNKED but the count is real)                          [EVS L]

Kilictas-Alpay arXiv:2506.22634; debunked S12, S30.

The paper *claims* exact pi(x) with O(log x) zeros via a Truncated
Gaussian kernel. Actually wrong (uncertainty principle, missing
sqrt(x) factor, AI-generated). **But the underlying observation that
the truncated explicit formula at x ~ 10^{10^8} converges with ~1200
zeros under their (incorrect) kernel is empirically reproducible at
small scales.** The arithmetic complexity is not in zero count but in
operating on enormous-digit x.

> Edge value: low because debunked, but flagged because the rumor that
> "you only need O(log x) zeros" recurs in the literature and chains
> sometimes silently invoke it. The *correct* count is and remains
> O(sqrt(x)/log(x)).

### E3.3 — Riemann's R(x) gets ~50% of the digits, exactly           [EVS H]

S29 measurements; `archive/sessions/session29_fresh_perspective.md`.

```
  x         pi(x)    |pi-li|   |pi-R|
  10^4      1229     17.14     2.15
  10^6      78498   129.55    29.35
```

R(x) converges in **only 10-20 Möbius terms** (it's polylog).  At
x = 10^{200}, R^{-1}(n) gives ~99 of 203 correct digits in 1.15s
(see `status/BEST_ALGORITHMS.md`). **The smooth approximation is NOT
the bottleneck.** Anything cheap that gives even one bit beyond the
R-ceiling is interesting.

### E3.4 — Spectral truncation gives *unique* prime in tiny window  [EVS H]

S21 finding;
`experiments/proposals/proposal2_crt_modular_results.md`,
`experiments/proposals/proposal1_spectral_sieve_results.md`.

Combining (a) R^{-1}(n) (smooth target), (b) ~K zeros gives an interval
of width O(x/K · log^2 x). With K = sqrt(x/log^3 x) zeros, the interval
is O(sqrt(x) · polylog), and after a CRT filter with M=30030 (six
small primes) **every n ≤ 10^4 has a UNIQUE prime in the filtered
interval** (7,462x search compression).

> Edge value: this is the cleanest "spectral + sieve" reduction.
> Bottleneck: computing pi(x) mod q for small q is mode-C circular at
> O(x^{1/2+epsilon}) cost via L-functions. **If pi(x) mod q for q ≤ 13
> were polylog, p(n) is polylog.** This is the single clearest "missing
> primitive" the project has identified.

### E3.5 — Borel-Padé regularises the divergent zero sum at moderate x [EVS L]

S45; `experiments/wildcard/borel_explicit_formula.py` and `_results.md`.

At x=500 with K=200 zeros and Borel-Padé (M,N)=(5,5):
- raw partial-sum error: -2.04
- Borel-Padé error: +0.22  (10x improvement)

Does NOT change asymptotic scaling (still O(sqrt(x)) zeros), but is the
*only* convergence-acceleration method tested that consistently
*improves* the error rather than worsening it (Richardson, Levin,
Weniger all WORSE than plain summation, see
`proven/convergence_acceleration_barrier.md`).

> Edge value: low but non-zero. Adaptive Borel order vs. x might give a
> useful constant-factor improvement to analytic prime-counting libraries.

### E3.6 — Numerical bug fix worth flagging                          [EVS M]

S43, S46; `experiments/proposals/stein_explicit_formula_decay.py`,
`experiments/proposals/proposal26_fejer_damped_explicit.py`.

`mpmath.li(x**rho)` is **wrong** for rho = 1/2 + i*gamma with large
gamma — it silently uses the principal branch of complex log, discarding
~ gamma·log x / 2*pi full rotations. Always use
`mpmath.ei(rho * log(x))` instead.

This invalidates several published-style numerical experiments that
"prove" sub-sqrt(x) convergence using `li(x^rho)`. **Any future chain
that includes a zero-sum numerical check must handle the branch cut
correctly.**

### E3.7 — Cesaro-Fejer window: 67% recovery at 1/3 the zeros        [EVS L]

S46; `experiments/proposals/proposal26_fejer_damped_explicit_results.md`.

Replacing the sharp truncation `sum_{|gamma|<T}` by the Fejer kernel
`(sin(gamma/2T)/(gamma/2T))^2` gives:
- T=100, sharp: 53% exact recovery on x ∈ [100, 3000]
- T=100, Fejer: 67% exact recovery
- T=300: 70% → 80%
- T=1000: 97% → 97% (saturates)

Constant-factor savings: Fejer recovers 2/3 of primes with ~1/3 of the
zeros. **Failures cluster on (a) x near a prime — Gibbs jump — and (b)
smooth x (small lpf)**, both circular features.

**S51 follow-up:** Stacking E3.5 (Borel-Padé) on top of E3.7 (Fejér) is
strictly worse than Fejér alone. Fejér + Borel-Padé recovery rate is flat at
3/8 across T ∈ {50, 100, 300, 1000} for x ∈ [10³, 10⁵], vs Fejér-alone 6/8
at T=1000. Borel-Padé on Fejér-windowed increments locks into a T-independent
asymptote (at x=10000: S = 1229.70 to 4 decimals across all T because the
Borel transform's 1/k! suppresses information from added zeros). The two
acceleration tricks **cannot be stacked**; each re-parametrises the same
√x-bounded zero-sum information. CLOSED in
`experiments/proposals/borel_fejer_hybrid_results.md`.

### E3.8 — Schoenfeld interval under RH                              [EVS L]

S33; `experiments/analytic/conditional/schoenfeld_cramer_results.md`.

Under RH, `|pi(x) - li(x)| < sqrt(x) ln(x) / (8 pi)` (Schoenfeld
explicit). Sieving the resulting Schoenfeld interval costs
O(sqrt(x) · ln x · ln ln x) = O(x^{1/2 + epsilon}), strictly better
asymptotically than Meissel-Lehmer O(x^{2/3}), but the constants are
prohibitive in practice and it does not approach polylog.

### E3.9 — pi(x) mod 2 sign predictability via Chebyshev bias        [EVS L]

S22; Rubinstein-Sarnak prediction.

`E(x; 4) > 0` for **99.52%** of x (vs theoretical 99.59%). Sign of
the prime race in mod-4 progressions is predictable; magnitude is not.
Provides O(1) bits, not O(log n). Used as a sanity check in chains
involving CRT modular reconstruction.

### E3.11 — K_min ~ 0.35 * x^{0.27} under per-x optimal zero ordering   [EVS M]

S13; `experiments/analytic/zero_compression*`;
`status/CLOSED_PATHS.md` line 380.

For each x, the **minimum number of zeta zeros** needed to achieve
|truncation error| < 0.5 — measured under per-x greedy-optimal selection
(not the natural ordering by gamma) — scales as

```
   K_min(x)  ~  0.35 · x^{0.27}
```

vs the **universal-ordering** floor of x^{0.5}. So per-x oracle reordering
strictly beats sqrt(x). The catch: the optimal ordering is a different
permutation of zeros for every x, and **no universal pre-permutation
achieves this**. Trying to identify the optimal ordering for x without
knowing pi(x) is circular.

> Why this is an edge: it draws a line between "information needed
> (x^{0.27})" and "information needed via universal algorithms (x^{0.5})".
> Bridging them would require a structural rule for **which** zeros matter
> at which x — a problem effectively as hard as pi(x) itself, but
> conceptually cleaner than zero-summation as a target.

**S62 caveat (NEW):** the universal-ordering K_min(x) measured at single
x values is **wildly non-monotonic** and dominated by GUE phase-luck.
Concrete (`experiments/analytic/conditional/k_min_extended/`):

| x       | K_min | K_min* (50-window stable) | sqrt(x) |
|---------|------:|--------------------------:|--------:|
| 10^3    |     0 |                        81 |    31.6 |
| 10^4    |     1 |                      1250 |   100.0 |
| 10^5    |     3 |                       572 |   316.2 |
| 10^6    |   N/A |    N/A (still off ±0.5 at K=2000) |  1000.0 |
| 10^7    |   563 |                      1912 |  3162.3 |

10^4 needs more zeros than 10^5; 10^6 hasn't settled at K=2000 even though
10^7 has. **A single-x K_min measurement is not a reliable estimator of
the asymptotic curve** — the residual oscillates around 0 with x-dependent
phase, and whether the rounding-window crossing happens "early" or "late"
is phase-luck. The asymptotic exponent stays at 1/2 (universal ordering).
Linear regression on log K* vs log x giving x^{0.275} from four points
is meaningless, two of them phase-luck-dominated. Future K_min
experiments should use the median over an x-interval, not a single x.

### E3.12 — Pascadi 2025 unconditional x^{5/8} equidistribution     [EVS M]

`literature/aggarwal_2025_analysis.md`; arXiv 2025 (Pascadi).

The unconditional level-of-distribution exponent for primes in arithmetic
progressions has improved from x^{1/2} (Bombieri-Vinogradov) toward the
conjectural x^{1/2+epsilon} ceiling. **Pascadi 2025: x^{5/8}**, achieved
without assuming the Selberg eigenvalue conjecture (which previously was
needed for any improvement past Bombieri-Vinogradov via Friedlander-Iwaniec).

> Why this is an edge: x^{5/8} = 0.625, vs the algorithmic threshold of
> x^{2/3} ≈ 0.667 — analytic level-of-distribution is now within 0.04
> of the algorithmic Meissel-Lehmer barrier. If unconditional level-of-
> distribution reaches x^{2/3+epsilon}, several conditional pi(x)
> algorithms (Lagarias-Odlyzko under GRH, etc.) would become unconditional.
> Not a polylog edge, but the **direction of analytic progress** for the
> first time creeps toward the algorithmic boundary.

### E3.10 — Cully-Hugill–Lee improved error term                     [EVS L]

`literature/state_of_art_2026.md` §6 (April S29 search update);
arXiv:2402.04272.

Riemann–von Mangoldt error improved to O(x/T) without log factor via
*averaged truncated Perron formula*. Still O(sqrt x) when T = sqrt x,
but cleaner constants.

### E3.13 — Bogomolny-Keating arithmetic correction is empirically absent at N=8000 [EVS M]

S49; `experiments/analytic/zeta_structure/large_n_correlations/large_n_battery.py`
and `_results.md`; `archive/sessions/session49_focus4_large_n_zeta.md`.
Companion methodology edge: E1.10.

The Bogomolny-Keating semiclassical conjecture predicts a non-universal
prime-arithmetic departure from GUE pair correlation:

```
   D(s) = R_2_emp(s) - R_2_GUE(s)
   D_BK(s; T) ≈ -(2/L²) sum_{p, k≥1} ((log p)² / pᵏ) · cos(2π s · k log p / L)
```

with `L = log(T/(2π))`. If detectable empirically, this would be the
**first non-universal structural deviation of zeta zeros from random-matrix
predictions** in the project — and would constitute exploitable arithmetic
information accessible from local pair-correlation data alone.

**S49 measurement at N=8000 (zero-height range up to ~T ≈ 9878):** zeta
sits **below** the gap-shuffled null on every BK statistic. Concretely:

| Statistic | Zeta | Gap-shuffled null mean ± std | z-score |
|-----------|------|------------------------------|---------|
| Pearson(D, BK_template), s>0.5 | +0.111 | +0.487 ± 0.035 | **−10.85σ** |
| Phase coherence ⟨cos(φ_p − π)⟩, p ≤ 50 | +0.544 | +0.624 ± 0.036 | −2.20σ |

Per-prime phase pattern matches the gap-shuffled null pattern (small
primes 2,3,5 have phases far from π; larger primes 29..47 have phases
near π). This is exactly the signature of *no real BK content* — purely
local-GUE behaviour with the same finite-N artefacts as a matched null.

> Why this is an edge: it is a **quantified empirical falsification** of
> a specific arxiv-grade theoretical prediction (BK arithmetic correction
> visible at this scale). Any chain that depends on BK-style residuals
> being recoverable from pair-correlation data is provably empty up to
> N=8000. **Combined with E7.1 (zeros linearly independent over Q in
> every measurable sense), the structural-arithmetic content of zero
> spacings is now empirically constrained to be at most O(1/√N) at every
> measurable scale.** Future Chain-D-style attacks that try to extract
> arithmetic from zero correlations need to find a different statistic
> than BK or operate at N >> 10⁴ — and even there, the gap-shuffled
> null calibration of E1.10 is mandatory.

This is a *negative* analytic edge: it doesn't open a chain, it closes
the BK-correction-extraction route up to current scale. Pairs with E7.1.

**S71 extension to Odlyzko's high-height tables (T ~ 10²⁰-10²¹):**
re-ran the BK probe at zero indices 10²¹+1..10²¹+10⁴ (`zeros4`,
L=44.6) and 10²²+1..10²²+10⁴ (`zeros5`, L=46.8). With a **proper
random-prime null** (template-shape-matched but with frequencies drawn
uniformly from [2,50]), the empirical BK Pearson is statistically
indistinguishable from the null:

| Block | L | empirical Pearson | random-prime null μ ± σ | z |
|-------|---|--------------------|----------------------------|---|
| zeros4_T1e20 | 44.6 | +0.0628 | +0.0630 ± 0.00021 | −0.94σ |
| zeros5_T1e21 | 46.8 | +0.0349 | +0.0345 ± 0.00037 | +0.93σ |

Direct sanity check: prime-frequency |Fourier amplitude| medians
(0.0141, 0.0100) are **not enhanced** over random-frequency amplitudes
in the same band (0.0148, 0.0101).

The previously-large gap-shuffled-null z-scores (~+30σ) at high height
are an artefact of E1.10's known bias inflated even further at large
L (Poisson-leakage in the gap-shuffled long-range tail anti-correlates
strongly with any oscillatory template — the inflation grows with L).

**Quantitative obstruction (E7.1 / E1.10 / E3.13 sharpening, S71):**
Empirical fit confirms `|BK_pred|_max · L² ≈ 13.6` (invariant across
L=44.6, 46.8) and `pair_rms ≈ 4/√N` (verified across N ∈ {2000, 8000,
10000}). Detection threshold for κ-σ Pearson:

```
  N_required(L; κ) ≥ (4κ / 13.6)² · L⁴ ≈ 0.09 κ² L⁴
```

At κ=3 detection: L=10 → N ≥ 10⁴; L=20 → N ≥ 1.5·10⁵;
L=44.6 → N ≥ 3.5·10⁵ (Odlyzko gives 10⁴, **35× short**); L=80 →
N ≥ 3·10⁷ (hopeless for any tabulated source). **L⁴ scaling is
structural** — doubling height requires 16× more zeros; the asymptotic
regime suppresses the BK signal faster than data accumulation can
compensate.

> Strengthens this edge from "empirically absent at N=8000" to
> "empirically absent at all heights through T ≈ 10²¹, AND requires
> N ≥ 0.81 L⁴ zeros for any future detection — a hard scaling barrier
> independent of computational budget."

This is the first quantitative version of the "BK-asymptotically-too-
weak" obstruction. Closes ATTACK_VECTORS §C1.

---

## 4. Sieve / combinatorial edges

### E4.1 — phi(W)/W bond-dim ratio = wheel-W gain in Mertens form    [EVS M]

(Same as E2.1.) Worth restating: every prime p added to the wheel
multiplies the bond-dim ratio by `(p-1)/p`. So the wheel-W sieve and
the base-W MPS reshape are *identical* compression mechanisms — one
is the Mertens product on the prime side, the other its tensor-network
shadow. **Useful as a translation table.**

### E4.2 — Lucy DP transition matrices are unipotent                 [EVS L]

S14; `experiments/circuit_complexity/lucy_dp_matrix_structure*.md`.

All eigenvalues of Lucy DP step matrices are 1; matrices are unipotent
with displacement rank 50-60% of dimension. Product is full-rank but
the unipotent structure means **no rotation/eigen tricks can speed up
the iteration**. The structure is "lower triangular with all-1 diagonal" —
free to invert, but not free to exponentiate.

### E4.3 — Hyperbola method floor-value count                        [EVS L]

Folklore + S24 measurement (`experiments/wildcard/ntt_sieve_results.md`).

Exactly **2*sqrt(x) - 1 distinct values** of `floor(x/n)` for n ≤ x.
This is the irreducible state-space size that any DP / sieve / Möbius
algorithm has to traverse. It IS the sqrt(x) barrier in disguise.

> Edge value: clarifies that "compress the floor-value set" = "break
> sqrt x". So the floor set itself is not the lever.

### E4.4 — Telescoped (Meissel-Lehmer) parallel depth = pi(x^{1/3})  [EVS L]

S12; `experiments/circuit_complexity/parallel_depth_scaling.py`,
`proven/circuit_size_barrier.md`.

The telescoped recursion has depth `O(log log x)` (≈ 1.6 log N with
constant = log(3/2)) but width O(x^{2/3}). **Free parallel depth, but
exponential in N=log x size.** Useful as a target: any circuit with
depth O(log log x) can match Meissel-Lehmer. To beat it we need to also
shrink the width.

### E4.5 — Lucy DP wheel-30 + prefetch: 2-35% wall-clock              [EVS L]

S41; `experiments/sieve/segmented_lucy.py`.

Cheap engineering edge for the project's V10 implementation. Confirms
that within-Lucy-DP optimisations are marginal; the 100-1000x gap
to primecount is in *infrastructure* (segmented sieve, OpenMP,
AVX-512), not in the algorithm.

### E4.6 — Prime gaps under residue + history features              [EVS L]

`proven/information.md`.

Conditional entropy of g(k):
- unconditional: 3.94 bits
- mod 30030 + 2 prev gaps: **0.54 bits irreducible**

86.3% of the gap entropy is reducible; the residue-decoupled 0.54 bits
is the irreducible piece (carried by larger primes 17, 19, 23, ...
not in the wheel). Capturing them all → wheel of primorial(sqrt(x)) =
sieving.

> Edge value: shows the smooth/oscillatory split also at the gap level.
> The 0.54 bits/gap is the gap analogue of the ~40% hard bits of p(n).

### E4.7 — Prime parity stream has h_mu deficit ~0.020 nats/step    [EVS L]

S48; `experiments/wildcard/causal_state_complexity.py`,
`novel/pseudorandomness_of_pi.md` measure #24.

For q(n) = 1{p(n) ≡ 5 mod 6}, the entropy rate is
**h_mu = 0.97686 nats/step** vs i.i.d. = 0.99708, a 0.0202-nat deficit.
This is the *first* genuine, measurable stochastic-memory deviation
from random — coming from Hardy-Littlewood pair correlations of
consecutive primes' residues mod 6 (e.g. twin pair (5,1)).

Excess entropy E plateaus at +0.026 nats above the Bernoulli null,
implying an epsilon-machine of size exp(0.026) ≈ 1.026 — basically
Markov-1.

> Edge value: tiny but *not* zero. Twenty-three earlier
> pseudorandomness measures rounded to "indistinguishable from random";
> this one is the first that doesn't.

---

## 5. Conditional / TC^0 edges

### E5.1 — BPSW correctness ⇒ PRIMES in TC^0                          [EVS H]

`novel/bpsw_tc0_reduction.md`; S13.

If BPSW (verified to 2^{64}, conjecture: no pseudoprimes) is
unconditionally correct then PRIMES is in DLOGTIME-uniform TC^0:
- MR(2): scalar pow x^k mod m is in TC^0 (Hesse-Allender-Barrington)
- Strong Lucas: 2x2 matrix powering, in TC^0 (Mereghetti-Palano)
- Jacobi symbol: GCD-style, in TC^0
- BPSW = AND of these → TC^0.

This is **strictly weaker than GRH** (BPSW is a specific computational
statement, not an analytic conjecture). Conditional on it, p(n) -> p(n)
becomes a TC^0 primality test in a polylog window from R^{-1}(n).
*Doesn't* immediately give pi(x) in polylog (counting is the bottleneck,
T1+T2 unaffected), but eliminates the testing-vs-counting confusion in
chain construction.

**S120 annotation (C4):** BPSW conditionality propagates **1-to-1**
through Aggarwal's binary-search wrapper, not amplified. A single
BPSW pseudoprime in the final bracket `[L, R]` (after Aggarwal
narrowing) shifts the answer by exactly one prime. Aggarwal narrowing
itself runs on `pi_lucy(x)`, which uses no primality testing — only
divisibility and small-prime sieving. Hence the BPSW conditional
enters only at the trailing walk, and propagates linearly. This is
structurally cleaner than naïve "BPSW everywhere" approaches that
would compound conditionals through a sieve. See
`algorithms/aggarwal_dusart_bpsw/`.

### E5.2 — QFT (Grantham) in TC^0                                    [EVS M]

S13. Quadratic Frobenius test is in TC^0 (operates in 2-D algebra =
2x2 MPOW). Error < 1/7710 per parameter set. Deterministic correctness
under GRH.

> Edge value: a second TC^0 primality oracle, redundant with E5.1 but
> useful if BPSW counterexample ever appears. Independent assumption.

### E5.3 — Growing-dim MPOW in TC^0 = the only OPEN frontier         [EVS H]

S11, S39, S47; `literature/matrix_powering_tc0.md`;
`status/CLOSED_PATHS.md` line 232.
**S127 quantitative refinement (C8):** depth-2 sign-threshold size for
PRIMES at N=6 is *exactly* `M* = 3` for every weight bound `W ≥ 3`
(plateau verified at W ∈ {3, 4, 8, 16, 32, 64}; all `M=2` cells UNSAT,
all `M=3` cells SAT, witness verified 64/64). At W=2, `M*(N=6) = 4`
(M=3 UNSAT proven 277 s, M=4 SAT 17 s). Combined with S84 column
enumeration `M*(W=1) = 6`, the curve is `(6, 4, 3, 3, 3, 3, 3, 3)` at
`W ∈ {1, 2, 3, 4, 8, 16, 32, 64}` — **structural floor reached at
W=3, not asymptotic**. Doubling W from 4 to 64 (16×) yields zero
gate reduction. Refines the prior "depth-2 sign-threshold restricted-
weight complexity" sub-question: weight is *not* a fixable axis.
See `experiments/constructions/d2_sign_threshold_w_m_tradeoff/`.

**S135 random-control sharpening (C8.b):** extended column enumeration
`Θ(N=6, W=2) = 30,898` distinct sign-threshold truth tables. Across
three independent random seeds {1, 5, 42}, density-matched random N=6
satisfies `M*(rand; W=2) ≥ 4 = M*(PRIMES; W=2)` (M=2,3 both UNSAT in
130–230 s per cell), confirming the direction of S127's F4 (PRIMES not
strictly easier than random at W=2) outside the C7-S89 oddness-
calibration regime. Strict-magnitude question (Δ = 0 vs Δ ≥ 1)
remains open: random M=4 W=2 cells UNKNOWN at CBC's 600 s budget.
Methodological side-finding: column-enum proves W=2 M=3 UNSAT 1.8×
faster than joint-ILP on PRIMES (157 s vs 277 s) and dominates joint-
ILP entirely for the random-target M=3 case (which joint-ILP could
not resolve). See `experiments/constructions/d2_sign_threshold_w_m_tradeoff/random_n6_resolution/`.

AKS in TC^0 ↔ matrix powering of polylog-dimensional matrices in TC^0.
Fixed-k MPOW IS in TC^0 (Mereghetti-Palano 2000) but k=polylog(n)
GROWING is **genuinely open**. Healy-Viola 2006 settles F_{2^n} via
Frobenius — has no analog over Z_n.

S47 closed the cyclotomic CRT sub-attack: AKS-prescribed r is prime
in 21/22 sampled n ∈ {10^2..10^6}, so x^r-1 = (x-1)·Phi_r(x) and the
non-trivial factor still has dimension r-1 (max_dim/r >= 0.94).

> Edge value: this is the **ONE genuine open question** at the
> TC^0/NC^1 boundary. A positive resolution puts PRIMES in TC^0; a
> negative resolution would separate TC^0 from NC^1 — both major
> complexity-theory results.

**S66 side-finding (folklore-grade):** the gcd-of-residual-coefficient
test, applied to composite n that fail Bernstein 2003's strengthened
AKS, **extracts a non-trivial factor of n in 13/13 sampled composite
failures including all 7 Carmichaels in [101, 410041]**. Mechanism:
residual coefficient at degree k is `≡ a^{n-k} · binom(n, k) mod n`,
and `gcd(binom(n,k), n)` is non-trivial whenever k is a base-p
sub-string of n for some prime p | n (Lucas + CRT). This is implicit
in the AKS proof and not novel-grade, but the empirical 100% extraction
rate is worth knowing — anyone proposing "use AKS residuals for
factoring" should start from this S66 number, not re-measure.

### E5.5 — Karchmer-Wigderson formula lower bound: 2^{N/2-O(1)}      [EVS M]

S18; `proven/circuit_size_barrier.md`.

The Karchmer-Wigderson theorem combined with the communication rank
2^{N/2-1} + 2 (E2.7) gives an **unconditional formula-size lower bound**
for pi(x):

```
   formula_size(pi(x))  >=  2^{N/2 - O(1)}
```

By Valiant's depth-reduction trick, this transfers to a general
circuit-size lower bound

```
   circuit_size(pi(x))  >=  2^{N/4}
```

These are the **only known unconditional super-linear lower bounds** for
pi(x) circuits (the Omega(log x) input-reading bound is trivial). They
fall short of separating pi(x) from TC^0 (which needs poly(N) lower
bounds), but they are real, proven, and not subject to Natural Proofs.

### E5.6 — PRIMES is in NONUNIFORM TC^0 unconditionally               [EVS M]

S13; `proven/complexity.md` table; Sorenson-Webster 2016.

For any **fixed** input length N, the deterministic Miller-Rabin base
set {2, 3, 5, 7, ..., 37} suffices for n < 3.317 × 10²⁴ (12 bases).
Each base is a scalar modular powering, in TC^0 (Hesse-Allender-
Barrington 2002). So:

```
   PRIMES (length N)  in  nonuniform TC^0  unconditionally
```

> Why this is an edge: it isolates the *uniformity* part of the question.
> All the "hardness" of PRIMES in TC^0 is in *uniform* construction
> (knowing which 12 bases work without testing), not in the existence
> of a depth-O(1) circuit. Combined with E5.1: BPSW correctness moves
> nonuniform to uniform.

### E5.7 — #AC^0 ⊆ TC^0 known; #TC^0 ⊆ NC open                       [EVS M]

`proven/complexity.md` §"Determinantal Complexity"; Allender-Gore 1993.

The counting-class inclusion `#AC^0 ⊆ TC^0` is a **theorem** (Allender-
Gore 1993). The next-level inclusion `#TC^0 ⊆ NC` is the central
remaining open question for our problem: under BPSW, it is the **only**
remaining gate between PRIMES in TC^0 and pi(x) in NC.

```
   pi(x) in NC  iff  #TC^0 in NC      (assuming BPSW correctness)
```

> Why this is an edge: it pinpoints the precise complexity-theoretic
> open problem. If somebody proves #TC^0 ⊆ NC, our problem becomes
> conditional only on BPSW (already verified to 2^{64}). If somebody
> proves the opposite, pi(x) in NC is impossible.

### E5.4 — Sorenson-Webster deterministic MR bases                    [EVS L]

`status/CLOSED_PATHS.md` line 253.

MR({2, 3, 5, 7, ..., 37}) is deterministic for n < 3.317 x 10^{24}
(Sorenson-Webster 2016). 12 bases × scalar pow each → uniform TC^0
**for fixed input length**, not asymptotic. Useful in any chain with
explicit small-x precomputation steps.

### E5.8 — Brandt 2024 is structurally welded to MKtP                 [EVS shape]

S51; `experiments/constructions/brandt_mktp/`;
`status/CLOSED_PATHS.md` last row (S51 Brandt).

Brandt 2024 (TCC, IACR ePrint 2024/687) proves `MKtP ∉ DTIME[O(n)]`
unconditionally via a length-monotonic depth-first traversal that
diagonalises against any candidate linear-time Kt-decider, using
1-Kt-randomness of Chaitin's Ω prefixes as the unique ingredient
that bypasses the black-box barrier (page 5 of the paper). The
proof relativizes and **does** thread the Natural Proofs barrier (T4).

> However it does NOT extend to fixed natural functions like
> pi(x) mod 2.  The hard string z is an oracle-dependent Kt-random
> prefix, not a fixed function; the contradiction uses self-referential
> Kt on both sides (`Kt(z) ≥ |z|` vs `Kt(z) ≤ |M| + log₂ t`); the
> Chaitin-Ω density argument has no analog for fixed totally-defined
> functions; and the bound is on uniform time (linear; conditional
> super-linear), not on circuits.  Brandt explicitly contrasts himself
> (page 4) with the Williams/Hirahara algorithmic-method approach that
> *would* yield circuit bounds and *is* subject to Natural Proofs on
> stronger classes — i.e., the price of relativisation is no circuit
> lower bound out of this technique.

Why this is shape-revealing: it is a **family-level closure on
diagonalisation-via-meta-complexity** approaches to Chain E (E5.3),
parallel to E7.10 (AKS modulus-twist orthogonality, AKS family) and
E7.11 (linear post-processing of zero sums, analytic family).  With
E7.10 + E5.8, Chain E is now closed for **both** known technique
families: AKS-style and diagonalisation-via-meta-complexity.  The
remaining unconstrained levers on E5.3 are non-AKS TC⁰ primality
tests or entirely-new lower-bound techniques.

**Per-bit extension confirmed (S105, C3).** The per-bit family
`{π_J(x) : J = 0..N-1}` (where `π_J(x) := bit J of π(x)`) inherits
all four obstructions O1-O4 wholesale. (O1) Each `π_J` is itself a
fixed total Boolean function — parameter-control of J is not the
oracle-dependence Brandt's TRAVERSE requires; the hard string of the
proof is generated by the run, not by the function-family index. (O2)
The π_J-decider answers "is bit J of π(x) set", not "is z complicated"
— no self-referential Kt inequality. (O3) The per-bit family has at
most `J*(N) ≈ N − log₂ N` non-trivial members per N (others are
identically zero), but they are pre-determined fixed strings, not
fresh Kt-random prefixes generated by traversal-path-dependent oracle
queries — Lemma 2 still fails. (O4) Uniform DTIME bound is unaffected
by bit-decomposition. The "minimal weakening of O1 admitting Brandt's
TRAVERSE" the C3 spec asked about does not exist on `{π_J}`. See
`experiments/constructions/brandt_per_bit/`.

---

## 6. Computational / oracle-style edges

### E6.1 — R^{-1}(n) gives ~99 of 203 digits in 1.15s at x = 10^200   [EVS M]

`status/BEST_ALGORITHMS.md`.

Concrete numbers: R^{-1}(n) is fast in *practice* well beyond the
theoretical "polylog gives ~50% of digits" prediction. At x = 10^{200},
50% of digits correct in ~1s. **The smooth-side primitive is engineering-
ready.**

### E6.2 — V10 hybrid: O(x^{2/3}) but with Newton-bisection-walk-
verify chain                                                         [EVS M]

`algorithms/v10_c_accelerated.py`; `status/BEST_ALGORITHMS.md`.

```
1. x0 = R^{-1}(n)                   [O(polylog)]
2. Newton: x0 += (n - C_Lucy(x0))*ln(x0)   [O(x^{2/3}) per call, ~4 steps]
3. Bisection in ~70-int bracket     [~6 steps]
4. Walk down to nearest prime       [O(ln^2(x)) trial divisions / Miller-Rabin]
```

Each layer is a primitive of independent value: layer 1 is polylog,
layer 4 is polylog (under BPSW), layer 3 is polylog. Only layer 2 is
exponential. **Replacing or eliminating layer 2 would break the barrier.**

### E6.3 — DCT of delta(n) is 99% sparse                              [EVS L]

S24, S35;
`experiments/wildcard/spectral_compression.py`,
`experiments/wildcard/wavelet_zero_compression.py`,
`experiments/wildcard/sublinear_correction.py`.

In the DCT basis, 99% of delta(n) energy lives in 10.4% of coefficients.
**But**: identifying *which* coefficients matter requires having computed
delta(n) (mode C), and the dominant indices change with N. Wavelet 99%
sparsity scales as N^{0.75} — sublinear, but not polylog.

**S24 follow-up (CRITICAL caveat):** the DCT compressibility is a
**finite-size effect that vanishes asymptotically**:

* At N = 5000, exact rounding of delta(n) needs 17.3% of coefficients
  (867/4999), with index spread spanning the full spectrum (median
  index 597).
* Hybrid K zeros + DCT: K=10 reduces needed coefficients from 867 to
  442; K=50 makes it WORSE (510). Mean residual stays ~2.5 regardless
  of K.
* FRI sampling density: **14.8% at small scale → 100% at x ~ 10^{100}**.
* Cross-scale correlation of DCT coefficients: only 0.742, NOT predictable.

So the 90/10 sparsity "edge" is illusory at our target scale x ~ 10^{100}.
This invalidates compressed-sensing / DCT-recovery chains for delta(n).

### E6.4 — 90/10 split in sieve update spectrum                      [EVS L]

S20; `experiments/wildcard/sieve_function_compression.py`.

Sieve update operations have 90% energy in the rank-1 component, while
log-Fourier of the sieve function has 90% energy in 10 modes. This is
the same smooth-vs-oscillatory split surfacing again — cleanly visible
in the **multiplicative-to-additive change of variables** (log-x → mod
2*pi/log_x).

### E6.5 — Hirsch-Kessler-Mendlovic O~(sqrt(N)) elementary algorithm [EVS M]

`literature/state_of_art_2026.md` §1.2; arXiv:2212.09857.

First elementary (non-analytic) algorithm matching analytic
O(x^{1/2+eps}) — uses Dirichlet convolutions via NTT to speed up the
combinatorial side. ~400x slower than primecount in practice but **matches
analytic-side asymptotics with no zeta zeros**. A chain that combines
HKM with E5.1 (BPSW-TC^0) and E1.3 (carry-boundary) is the closest the
project gets to a believable polylog architecture.

### E6.7 — HKM time-space tradeoff curve O~(N^{8/15}) / O~(N^{1/3}) [EVS M]

`literature/aggarwal_2025_analysis.md`; arXiv:2212.09857 (Helfgott-Kessler-
Mendlovic);
**pillar placement:** `experiments/constructions/pillar_tradeoff_diagram/` (S81, C6).

HKM is not just a single point at O~(sqrt(x)) (E6.5); it is a tradeoff curve.
With NTT-based fast multiplication on Dirichlet series:

```
   time = O~(N^{8/15})  using  space = O~(N^{1/3}).
```

This is the **only known sub-x^{2/3} time achievable simultaneously with
sub-x^{1/2} space** for pi(x). Beats Meissel-Lehmer (x^{2/3} time / x^{1/3}
space) on time AND beats classical analytic (x^{1/2+eps} time / x^{1/2}
space) on space.

> Why this is an edge: it is the **first elementary algorithm to do strictly
> better than both classical sieve and explicit-formula on a Pareto axis**.
> The 8/15 = 0.533 exponent is between 1/2 and 2/3. Combined with the
> Meissel-Lehmer pebbling lower bound T*S >= Omega(x^{5/6}/ln x) (E7.6):
> HKM achieves T*S ~ x^{8/15+1/3} = x^{13/15} = x^{0.867}, almost matching
> the lower bound x^{5/6} = x^{0.833}. **The sieve route is asymptotically
> tight up to x^{0.034}.** Polylog cannot come from any further
> sieve-pebbling improvement.

**S81 pillar-placement extension (composition C6):** HKM's `(8/15, 1/3)`
point is achievable on the **floor pillar only** (E7.7). Per-pillar Pareto
classification (catalog of 14 algorithms across the three pillars) shows:
(i) HKM dominates every other floor-pillar entry elementwise (Lucy DP,
Lucy+Fenwick, Meissel-Lehmer classical, Gourdon, primecount); (ii) HKM has
no cross-pillar dominators — no zero-pillar or prime-pillar algorithm
achieves both `T <= N^{8/15}` AND `S <= N^{1/3}` simultaneously; (iii)
pillar dominance regions are non-overlapping: prime/zeta pillars share
the time-only minimum at `alpha = 1/2`, the floor pillar uniquely accesses
the space-only minimum at `beta = 1/3`, and the floor pillar uniquely
achieves `T*S = 13/15`. See
`experiments/constructions/pillar_tradeoff_diagram/`.

### E6.8 — Dusart bracketing: p_n in interval of width n            [EVS M]

`literature/aggarwal_2025_analysis.md`; Dusart 2010, 2018.

For n >= 6:

```
   n (log n + log log n - 1)  <=  p_n  <=  n (log n + log log n).
```

The bracket width is exactly **n** (logarithmically narrow on the scale
of p_n). This is what allows Aggarwal-style binary search (E6.6) to take
**only O(log n) calls to pi(x)** rather than O(log p_n) ~ O(log x) calls,
saving an extra log factor.

> Why this is an edge: it is the only published explicit bracket of width
> not depending on log x. Combined with E1.3 (4-bit sigmoid carry boundary)
> and the carry-propagation theorem: the *bit-level* structure of p_n
> within the Dusart bracket is what determines the polylog gap. A polylog
> p(n) algorithm needs to *navigate* this n-wide window in polylog time;
> Dusart guarantees the window itself is polylog-describable, the navigation
> is the bottleneck.

**S120 annotation (C4):** quantitative leverage measured. A pure
BPSW-walk from 2 needs ~p_n BPSW tests, while a Dusart-bracketed walk
from L_dusart needs only ~width / 2 = n / 2 tests on average. Ratio
2 p_n / n ~ 2 log p_n predicts a **30-50× speedup** at n = 10⁴-10⁷;
empirically observed 21×/34×/53× (within factor 2 of prediction). The
Dusart bracket is the cheap value-add of the C4 composition.

### E6.9 — primecount asymptotic hierarchy                          [EVS L]

`literature/primecount_analysis.md`; Lucy 2008, Walisch primecount.

The state-of-art combinatorial sieve sequence:

| Algorithm | Time | Space | Engineering |
|-----------|------|-------|-------------|
| Lucy DP basic | O(x^{3/4}) | O(sqrt x) | direct DP |
| Lucy + Fenwick | O(x^{2/3} log^{1/3} x) | O(sqrt x) | log update |
| Gourdon variant | O(x^{2/3} / log^2 x) | O(sqrt x) | two-param family |
| primecount production | same Big-O | O(sqrt x) | POPCNT counter, mod-240 wheel, multi-level (O(log log x) per event) |

Combined ~10^5 asymptotic constant-factor improvement over naive Lucy DP.
Mod-240 wheel encoding (8 bits/byte for residues mod 30) is the **concrete
realisation of E2.1 / E4.1's wheel-W = phi(W)/W bond-dim theorem** — a 4x
storage compression equal to phi(30)/30 * 30/8 = 8/30 * 30/8 = 1, with the
factor coming from 8 residues per primorial.

> Why this is an edge: it is the **engineering ceiling of the entire
> sieve route**. primecount/kim-walisch is within the polylog factor
> log^{1/3}(x) of HKM (E6.7); breaking past requires NOT a sieve speedup
> (those have been wrung dry) but a fundamentally different algorithmic
> structure.

### E6.6 — Aggarwal 2025 binary search optimality                     [EVS H]

`literature/state_of_art_2026.md` §1.4; arXiv:2510.16285.

p(n) computable in O(sqrt(n) (log n)^4) unconditionally via binary search
on pi(x). **Proves sieve methods CANNOT asymptotically beat
binary-search-on-pi(x).** This is *not* polylog, but it is the tightest
known upper bound and **explicitly identifies pi(x) as the bottleneck**.

> Edge value: **frames the entire research target.** Any polylog p(n)
> algorithm reduces (via Aggarwal) to a polylog pi(x) algorithm.

**S120 annotation (C4 composition E6.6 + E6.8 + E5.1):** the
Aggarwal-pure binary search has a constant-factor inefficiency in the
trailing `O(log K)` pi-calls when narrowing the bracket to width K.
Replacing the trailing narrowing with a BPSW walk over `O(K)`
candidates is strictly faster when `pi(x)` cost exceeds `K * BPSW`
cost. Empirical optimum K satisfies `K* ~ pi_cost / bpsw_cost`:
pure-Python Lucy regime gives `K* = width` (no Aggarwal narrowing
useful), C-accelerated Lucy gives `K* ~ sqrt(width)` (mid-regime),
HKM projection gives `K* ~ const` (Aggarwal-pure dominates). The
asymptotic O(sqrt(n) log^4 n) bound is preserved; the composition is
a **practical-constant tightening only**. See
`algorithms/aggarwal_dusart_bpsw/`.

---

## 7. Negative-but-shape-revealing edges

These are not "edges" in the affirmative sense, but they constrain the chain.

### E7.1 — Each zeta zero contributes "independent" information     [EVS shape]

S25, S45, S20, S57, S123; `experiments/wildcard/zero_sum_acceleration.py`,
`experiments/analytic/zeta_structure/triple_correlation.py`,
`experiments/analytic/zeta_structure/n_correlations_4_5_6/`.

PSLQ on first 1000-2000 zeros: **0/1225 pairs, 0/4060 triples, 0/400
cross-block** integer relations at 30-60 digit precision. Pair-correlation
matches GUE to RMS deviation 0.0864. **Zeros are linearly independent
over Q in every measurable sense.**

S57 extension to **3-point correlation** at N=2000: bulk RMS deviation
from GUE sine-kernel determinant `rho_3(s_1, s_2) = 1 - K(s_1)^2 -
K(s_2-s_1)^2 - K(s_2)^2 + 2 K(s_1) K(s_2-s_1) K(s_2)` is 0.0875 — same
noise floor as the pair-correlation test. Cumulant rigidity:
`c_3` of zero count in disjoint windows of length L ∈ {1..32} stays at
~10⁻³ while a Poisson process would give c_3 ranging from 1 to 32 — a
factor-10⁴ separation, fully consistent with GUE rigidity.

**S123 extension to orders 4, 5, 6 at N=8000** along three independent
probes:

- **R_n along equally-spaced slices** `R_n(0, s, 2s, ..., (n-1)s) =
  det[K((j-i)s)]_{i,j}`: max |z_vs_theory| = 2.36σ at n=4, s=2.0;
  n=5 raw 6σ at s=2.0 reproduced by matched-finite-N GUE Monte Carlo
  (K=20 batches × 1200 central evs each), so attributable to Poisson
  shot noise on `n_ref · tol^{n-1} ≈ 12` expected tuples — NOT
  arithmetic content.
- **k-th nearest-neighbor spacings P_k(s)** for k ∈ {0..5} (probing
  up to 7-point correlation): rms vs GUE pool ∈ [0.09, 0.22] growing
  monotonically with k; per-bin |z_vs_GUE_pool| ≤ 1.5× sample-noise
  scale at every k.
- **Higher cumulants κ_n(L)** for n ∈ {3..6} at L ≥ 8: all
  |z_vs_GUE_batches| < 2.1σ. The variance `k_2` retains its ~9σ
  rigidity z-score at L=32 (the GUE-rigidity signal). At L=1, 2
  cumulants k_4, k_6 deviate strongly between Riemann-vM unfolded
  zeta and semicircle-unfolded GUE — finite-N unfolding mismatch,
  not arithmetic content (eliminated at L ≥ 8).

**No higher-order cluster structure exists at any tested order up to 6
within the resolution allowed by N=8000 / γ ≤ 8148.** Conrey-Snaith's
conjectural higher-order arithmetic correction scales as `1/L²` where
L = log(γ/2π) ≈ 6.5; the empirical noise floor for n-correlations at
this N is `1/√(n_tuples)` which exceeds the predicted correction by
1.5+ orders of magnitude. Detection demands Odlyzko-tabulated zeros
at γ ≥ 10⁶ (cf. S71 closure of §C1; same `N ≥ 0.81 L⁴`-style scaling
barrier).

> Composes with: E1.10 (gap-shuffled is the right null for prime-
> frequency probes; S123 confirms it is the WRONG null for higher-
> order GUE-vs-arithmetic discrimination — destroys GUE rigidity);
> E3.13 (BK arithmetic correction below noise floor); S71 (Odlyzko-
> height closure of §C1 shares the same `1/L²` scaling barrier).

> Cross-domain references: Mehta 2004 *Random Matrices* (3rd ed.) §6.2;
> Hough-Krishnapur-Peres-Virág 2009 *Zeros of GAFs and DPPs* §1.2;
> Conrey-Snaith 2007 arXiv:0710.5304 (target conjecture).

### E7.2 — Floor-division is not a ring homomorphism                [EVS shape]

S20, S35;
`experiments/wildcard/adelic_prime_count_results.md`,
`experiments/proposals/proposal17_padic_lifting_results.md`.

`floor(x/d) mod q != floor((x mod q)/(d mod q))` — kills ALL purely
modular decompositions of pi(x). Every CRT-based attack ultimately
trips on this.

### E7.3 — Holographic projection: bits are entangled                [EVS shape]

S20, S38; `archive/ephemeral/wildcard_findings.md`.

For x = 10^{100}, the ~173 "hard" bits of pi(x) are a holographic
projection of ~10^{50} zeta-zero contributions. **No single hard bit is
independently computable** — every bit needs O(sqrt x) zeros. So per-bit
attacks (E1.3 sigmoid) cannot reduce to "compute one bit each, get O(1)
per bit."

### E7.4 — Inversion iteration is non-contracting                    [EVS shape]

S44; `experiments/sieve/inversion_search_results.md`.

Fixed-point iteration `p_{k+1} = R^{-1}(n + sum_rho R(p_k^rho))` has
|f'| oscillating ~ |sum_rho R'(p^rho)|, which has GUE-random phases,
so |f'| is not bounded < 1. Iterating does NOT converge geometrically.
**This rules out *any* attempt to "self-correct" R^{-1}(n) using its
own truncated zero sum.**

### E7.5 — Communication rank caps lower bounds at log^2 x          [EVS shape]

S20, S23; `experiments/circuit_complexity/space_time_tradeoff_results.md`.

D(pi) ≤ N = log x, so all communication-complexity-based time-space
tradeoffs give T*S ≥ Omega(log^2 x), polylog. **The communication route
to proving polylog impossible is closed.** Only circuit lower bounds
could prove it, and those face Natural Proofs (T4).

### E7.6 — Meissel-Lehmer DAG pebbling: T*S ≥ Ω(x^{5/6}/ln x)       [EVS shape]

S20, S23; `experiments/circuit_complexity/space_time_tradeoff_results.md`;
`proven/circuit_size_barrier.md`.

Lucy DP DAG: depth = pi(sqrt(x)), width = O(sqrt(x)). Black-pebbling
gives the time-space tradeoff `T*S >= Omega(x^{5/6}/ln x)`. This rules
out simultaneous polylog time AND polylog space for **any** Lucy/Meissel-
Lehmer-style sieve algorithm.

> Why this is shape-revealing: it is the **strongest known lower bound
> for any specific pi(x) algorithm**, but **algorithm-specific** — a
> non-sieve approach could in principle bypass this DAG entirely.
> Establishes that the sieve route is provably stuck at sqrt(x)-class.

### E7.8 — Curve-encoding cost is at least O(g^2) Frobenius        [EVS shape]

S24; `experiments/wildcard/etale_counting_results.md`.

To represent p(n) as a Frobenius eigenvalue of a curve over F_2, the curve
must have genus

```
   g >= p(n) / (2 sqrt 2)             (Hasse-Weil bound saturation)
```

For x = 10^{100}, g ~ 10^{102}; computing Frobenius eigenvalues of a
genus-g curve costs O(g^2) by Schoof-Pila or O(g log g) by Couveignes-Lercier
in the best case. **Either way, exceeds the trivial sqrt(x) by many orders
of magnitude.**

> Why this is shape-revealing: closes the *entire algebraic-geometry route*
> (via Weil conjectures / etale cohomology / curve point counting) by a
> simple counting argument. Any curve big enough to encode pi(x) has too
> many points to count quickly. **Algebraic geometry is at least as hard
> as direct sieving for our problem.** Pairs with E7.2 (floor-division
> not a ring homomorphism) to close most "structured arithmetic" routes.

### E7.9 — Hierarchical Phi recursion calls/x ratio is constant     [EVS shape]

S24; `experiments/wildcard/hierarchical_sieve_results.md`.

For the Meissel-Lehmer recursion `Phi(x, a) = Phi(x, a-1) - Phi(x/p_a, a-1)`,
the number of recursive calls per unit x is **0.03-0.04, constant in x** for
x up to 10^9. This is the *fast-multipole hope* (analytic far-field +
exact near-field) made empirical: there is **no fast-multipole-style
constant-factor compressibility beyond what the recursion already exploits**.
The sub-tree sparsity (7.5% non-trivial Möbius terms at N=1000) **does NOT
scale**: Psi(10^100, 10^50) ≈ 10^100 — the smooth-part shortcut works only
when the smoothness parameter y << x.

> Why this is shape-revealing: rules out the Greengard-Rokhlin "FMM-for-
> primes" analogy at the **constant scaling level**, not just on technical
> grounds. The recursion is Theta(x), not o(x), in calls.

### E7.10 — AKS modulus-twist orthogonality theorem               [EVS shape]

S61, S64, S66 (combined evidence);
`status/CLOSED_PATHS.md` lines 714, 719, ~722;
`archive/sessions/session{61,64,66}_*.md`.

The three FOCUS-1 sub-attacks closed in 2026-04-26 produced a clean
structural conjecture, well-evidenced if not formally proved:

> **Modulus / coefficient-ring / gcd-strengthening twists of AKS are
> orthogonal to circuit depth.** Every AKS variant tested preserves
> the matrix-power dimension `r = polylog(n)` and therefore reduces
> to E5.3 (growing-dim MPOW in TC⁰), independently of the surrounding
> algebraic choices.

Concrete twist-by-twist evidence:

| Variant | What changed | What stayed | Source |
|---------|--------------|-------------|--------|
| Non-cyclotomic ring `Z_n[x]/(x^d+a)` | Modulus structure (Eisenstein vs cyclotomic-flavour) | r×r MPOW over Z_n | S61, line 714 |
| Frobenius transplant `(Z/qZ)[x]/(x^r-1)`, `q ∣ n−1` | Coefficient ring (Z_n → F_q) | r×r MPOW over F_q | S64, line 719 |
| Bernstein 2003 strengthened gcd | Adds gcd-of-residual certificate | r×r MPOW over Z_n + integer gcd | S66, line ~722 |

In all three cases the "interesting" structural change (Eisenstein
twist breaking Carmichael alignment, q-power Frobenius giving binomial
decomposition, gcd extraction recovering 13/13 Carmichael factors) is
*orthogonal* to the depth question. The bottleneck primitive is the
matrix-power chain itself, not any algebraic surrounding.

> Why this is shape-revealing: it is a **search-space constraint on the
> entire AKS family of TC⁰ approaches.** Any future AKS-style
> construction whose only novelty is a modulus / ring / gcd choice is
> orthogonal to E5.3 and provably cannot resolve Chain E. **Chain E is
> "computationally cornered"** within the AKS family per TODO.md's
> backup objective. The remaining levers on the only open problem are
> non-AKS in flavour: Brandt MKtP (FOCUS-3), a fundamentally new
> lower-bound technique, or a non-AKS TC⁰ primality test.

This complements E7.6 (Lucy-DAG pebbling closes the *sieve* family
asymptotically): E7.6 + E7.10 together close two of the three known
informationally-complete encodings (E7.7) at the algorithmic level,
leaving only zero-summation routes — themselves blocked by E7.1
(zeros are linearly independent) and E7.4 (inversion non-contracting).

### E7.12 — Cayley graph spectrum probes ω(n), not primality          [EVS shape]

S79; ATTACK_VECTORS §A.A3;
`experiments/circuit_complexity/cayley_spectral_primality/`;
CLOSED_PATHS row at session 79.

For n coprime to the generator set S, the adjacency spectrum of the
Cayley graph Cay((Z/nZ)*, S ∪ S⁻¹) has at least **2^ω(n) integer
eigenvalues**, where ω(n) is the number of distinct prime factors of n.

Proof: by CRT, (Z/nZ)\* ≅ ∏ (Z/p_i^{a_i}Z)\*. For odd p_i each factor
has even order so its character group has 2-torsion of order 2; the
2-torsion of the full character group has size 2^ω(n). For χ in the
2-torsion subgroup, χ(s) ∈ {±1} for all s, so the eigenvalue
λ_χ = Σ_{s∈S∪S⁻¹} χ(s) is an integer. ∎

Empirically verified (533 values of n × 3 generator sets, 0 violations).
Bound is sharp in 84/533 cases.

> Why this is shape-revealing: any **fixed-generator abelian Cayley
> graph spectral feature** is a function of (Z/nZ)\*'s decomposition
> into prime-power factors, hence detects ω(n), not primality. Since
> (Z/p^kZ)\* is cyclic for odd p (Gauss), prime n and prime power
> n=p^k both give cyclic unit groups with ω=1; they are
> spectrally indistinguishable up to discrete-log alignment. Distinguishing
> primes from prime powers requires a separate perfect-power test
> which inherits AKS-class growing-dim MPOW depth (E5.3, E7.10). Adds
> a search-space constraint to the spectral-primality family — closing
> A3 from ATTACK_VECTORS.md and ruling out the entire fixed-generator
> abelian-Cayley spectral approach.

Distinct from CLOSED_PATHS lines 354, 383, 384, 385 (which used
primes-as-generators or graph index — circular). E7.12 is the
**index-fixed, generator-fixed** version of the failure: the graph
itself is buildable without primes, but the spectrum doesn't decide
primality.

**S80 quantum extension:** Szegedy quantisation
of these Cayley walks inherits the same poly-mixing barrier. The
discriminant matrix `D(x,y) = sqrt(P(x,y) P(y,x))` has the same
spectrum as the (lazy) walk for reversible chains, so eigenvalues
`cos(θ_k) → e^{±2iθ_k}` give phase gap `2 arcsin(sqrt(δ))`. Empirical
sweep over 14 primes N ∈ {31, …, 1009} fits classical mixing
`t_C(N) ~ N^{0.789}` and Szegedy mixing `t_Q(N) ~ N^{0.414}` —
quadratic speedup confirmed but not polylog. The structural barrier
(spectrum probes ω(n)) is unchanged by quantisation.

### E7.15 — L(s, Δ) Hecke partial-sum subspace is ~3× more obstructed from spanning π(x) − Li(x) than matched-multiplicative random Sato-Tate ensembles [EVS M]

S118; ATTACK_VECTORS §B2;
`experiments/algebraic/automorphic_l_function_basis/`;
CLOSED_PATHS row at session 118.

Empirical OLS-projection ratio at finite N:

> **Claim.** For y(x) := π(x) − Li(x) sampled at M = 200 log-uniform
> anchors x_i ∈ [N/50, N−1], and feature matrix F_τ[i, k] := Σ_{n≤x_i}
> a(n) cos/sin(t_k log n) where a(n) = τ(n)/n^{11/2} (normalised
> Ramanujan tau, level-1 weight-12 cusp form Δ) and t_k log-uniform in
> [1, 50] at K twist frequencies, the 5-fold cross-validated test rms
> residual rms(y − F_τ ĉ_τ) is **≈ 2.83× larger** than for matched-
> multiplicative random Sato-Tate ensembles F_mult (a_mult(p) drawn
> from Sato-Tate density (2/π) sin² θ at primes, Hecke-recursion at
> prime powers, multiplicativity at composites). Empirically:
>
> - rms(F_τ) / √N = 0.0415, 0.0424, 0.0406 at N = 5k, 10k, 20k.
> - rms(F_mult) / √N = 0.0183, 0.0157, 0.0146.
> - **Ratio ρ := rms(F_τ) / rms(F_mult) = 2.83 ± 0.02 across (N, K) ∈
>   {5k, 10k, 20k} × {4, 6, 8, 10, 12} (15 cells).**
> - Z-score per cell: 17–58σ.
> - β_τ ≈ 0.05 (Hecke fit captures √x scaling — residuals are flat in
>   x), but the constant residual is ~3× larger than F_random.
> - Effective rank of F_τ at K=8 is 15/16 = 0.94, matching F_iid
>   (14.9/16) and F_mult (15.5/16) — NOT a rank-deficiency artifact.

Mechanism (closure mode E):

> By Mellin-Perron, F_τ_k(x) is dominated by L(s, Δ) zero contributions
> at heights γ_l^Δ. By the explicit formula, y oscillates at heights
> γ_l^ζ. Sato-Tate equidistribution / Katz-Sarnak GUE statistics imply
> {γ_l^Δ} and {γ_l^ζ} are GUE-distributed independently. So F_τ is a
> **narrow-band basis at the wrong heights**; F_random_mult is a
> **broadband basis covering all heights**. Both have the same K and
> similar effective rank, but F_τ spans a "wrong subspace" of feature
> space. The empirical 3× obstruction ratio is consistent with
> (effective spectral coverage)^{−1/2} for K ≈ 8 narrow- vs broad-band
> bases.

> Why this is shape-revealing: §B2 was the recommended next frontier
> attack per the post-S117 critique — the only major function-field /
> spectral candidate untouched. The classical "Sato-Tate kills
> automorphic attacks on π(x)" intuition is now numerically anchored:
> the Hecke partial sum subspace at K=8 frequencies recovers only
> ~89% of y energy (R² in-sample), versus ~98.6% for matched random
> ensembles, with Z = 30σ obstruction. This is the **first
> automorphic-spectral measurement of the project**, complementary to
> the Riemann-side measurements E7.1 / E1.10 / E3.13.

Cross-domain reference: Sarnak 2008 "Letter to Bombieri on the
Riemann hypothesis" (heuristics for {γ^Δ} vs {γ^ζ} independence);
Katz-Sarnak 1999 *Random Matrices, Frobenius Eigenvalues, and
Monodromy* (AMS Coll. 45) for the GUE-independence framework;
Iwaniec, *Spectral Methods of Automorphic Forms* (AMS, 2002) for
the partial-sum / Mellin-Perron machinery.

### E7.16 — Friedman / Ramanujan ratio of `Cay(Z/NZ, ±primes < N^c)` is indistinguishable from a parity-and-support-matched random subset of [3, N^c) within ±2σ [EVS shape]

S125; ATTACK_VECTORS §D.D20;
`experiments/algebraic/friedman_ramanujan_prime_cayley/`;
CLOSED_PATHS row at session 125.

Empirical Cayley-spectral measurement at finite N:

> **Claim.** For `N` prime in `{509, 1009, 4001, 16001, 65537}`,
> `c ∈ {1/2, 2/3}`, the abelian Cayley graph
> `G_N = Cay(Z/NZ, S_N)` with `S_N = {±p : p prime, p < N^c}` of degree
> `d = 2 |{p < N^c}|` has FFT-computable second eigenvalue
> `λ_2 = max_{k ≠ 0} |Σ_{p < N^c} 2 cos(2π pk/N)|` and Ramanujan ratio
> `r_N = λ_2 / (2 √(d-1))`. Empirically:
>
> - `r_N(prime)` ranges 2.05 (N=509, c=0.5) → 11.30 (N=65537, c=2/3),
>   sub-Ramanujan by orders of magnitude (Friedman 2008 reference is
>   r ≈ 1).
> - vs UNIFORM random subset of `Z/NZ` (Friedman reference): Z = +5 to
>   +66 — TRIVIAL bounded-support FFT spike at `k = 1`.
> - vs SUPPORT-MATCHED random subset of `[2, N^c)`: **Z = +0.68 to
>   +1.87 across all 10 cells (within ±2σ noise)**.
> - vs PARITY-MATCHED random odd subset of `[3, N^c)`, MINOR-arc band
>   `k ∈ [N/4, 3N/4]`: Z = -31 to -15622 — entirely the single-prime
>   `p = 2` parity-spike artefact.
> - vs PARITY-MATCHED, with `p = 2` removed from prime set:
>   **Z = +0.51 to +2.07 across all 10 cells (within ±2σ noise)**;
>   sign positive in 10/10 (sign-test p ≈ 0.001) but no individual cell
>   exceeds Bonferroni-3.4σ.

Mechanism (closure mode E — explicit reduction to two trivial controlled
effects):

> The two non-trivial spectral spikes are at `k = 1` (low-frequency
> major-arc) and `k ≈ (N-1)/2` (parity frequency). The first is a
> deterministic bounded-support FFT spike: any subset of `[2, N^c)`
> gives `λ_1 ≈ |S| (1 - O(M²/N²))` because `cos(2π p / N) ≈ 1` for
> `p ≪ N`. Vinogradov's prime-exp-sum bound `|Σ_p e(pα)| ≪ M (log M)^A
> / √q` does not apply at `α = 1/N` since `q = N` is not bounded by a
> fractional power of `M`. The second spike at the parity frequency is
> the all-odd alignment: for primes > 2, `(-1)^p = -1`, giving
> `|λ_{(N-1)/2}| ≈ d - 4` (the `-4` is the single even prime `p = 2`
> contributing `(-1)^2 = +1` instead of `-1`). Random odd-only subsets
> achieve `|λ_{(N-1)/2}| ≈ d` — a 4-unit fixed gap explains the
> thousand-σ Z scores at large N (random-control std → 0 because
> parity alignment is deterministic up to negligible cosine
> corrections). After both controls, the residual Hardy-Littlewood
> mod-q signal is invisible at the scales tested (`p = 2` artefact
> dominates over any mod-3 / mod-5 sieve-cancellation correction).

> Why this is shape-revealing: it adds a **NEGATIVE-shape** edge in the
> abelian Cayley spectral category, distinct from line 752 / E7.12
> (S79: fixed-generator (Z/nZ)\* spectrum probes ω(n), not primality)
> and line 754 / E7.13 (S80: Szegedy quantum walks on arithmetic
> graphs do not yield polylog mixing). E7.16 closes the *prime-as-
> generator* abelian Cayley spectral gap as carrying no detectable
> arithmetic content beyond support concentration and parity. Joins
> the HL-detection family — E2.13 (Gowers), E2.14 (Anderson),
> E2.15 (alg. immunity), E2.16 (DPP), E2.17 (PH), E2.18 (Liouville
> Anderson), E2.19 (subword complexity) — but the HL singular-series
> mod-3 / mod-5 component is here invisible because the bounded-
> support and parity dominators are an order of magnitude larger.
> Composes with: E2.13, E2.14, E2.18, E7.12, E7.15.

Cross-domain reference: Friedman 2008 "A proof of Alon's second
eigenvalue conjecture and related problems" Memoirs AMS 195 =
arXiv:cs/0405020. Hoory-Linial-Wigderson 2006 "Expander graphs and
their applications" Bull. AMS 43, 439. CROSS_DOMAIN_TECHNIQUES §1
row "Random regular graph spectral gap (Friedman)" promoted
PROPOSED → USED-E with this edge.

### E7.17 — Hodge L_1 max-eigenvalue multiplicity of the coprimality flag complex equals the Bertrand-prime triangular count [EVS M]

S126; ATTACK_VECTORS §D.D22;
`experiments/topological/hodge_coprimality_flag/`;
CLOSED_PATHS row at session 126.

Empirical Hodge-spectral identity for an arithmetic flag complex:

> **Claim (verified at 9 values of N spanning a 16× range).**
> For the coprimality flag complex `K_N := \{σ ⊆ [2, N] : σ is pairwise
> coprime\}` with N ≥ 3, the combinatorial Hodge Laplacian
> `L_1 = B_1^T B_1 + B_2 B_2^T` has top eigenvalue
>
>     λ_max(L_1, K_N) = |V| = N − 1
>
> exactly (numerical residual ≤ 1e-13), with multiplicity
>
>     mult(λ_max = N − 1) = C(k + 1, 2) = k(k + 1) / 2,
>
> where `k = π(N) − π(N/2)` is the Bertrand-prime count in `(N/2, N]`.
>
> Verification:
>
> | N    | k  | C(k+1, 2) | empirical mult |
> |------|----|-----------|----------------|
> | 8    | 2  | 3         | 3              |
> | 12   | 2  | 3         | 3              |
> | 16   | 2  | 3         | 3              |
> | 24   | 4  | 10        | 10             |
> | 32   | 5  | 15        | 15             |
> | 48   | 6  | 21        | 21             |
> | 64   | 7  | 28        | 28             |
> | 96   | 9  | 45        | 45             |
> | 128  | 13 | 91        | 91             |

Mechanism (closure mode E — single-page Friedman/Horak-Jost exercise):

> By Bertrand's postulate, for every N ≥ 3 there exists a prime
> `p ∈ (N/2, N]`. Such a prime is **universal** in `G_N`: for any
> `v ∈ [2, N] \setminus \{p\}`, `gcd(p, v) = 1` (since `p > N/2 ≥ v`
> in most cases and `p ∤ v` always except when `v = p`). The set `U`
> of all such Bertrand primes is itself pairwise coprime (distinct
> primes), so `U` induces `K_k` in `G_N`. Furthermore `U` is in the
> *full join* with `V \setminus U`: every non-Bertrand-prime vertex is
> adjacent to every `p ∈ U`. Hence the flag complex factors as a join
>
>     K_N = Δ^{k-1} \ast F(H),
>
> where `Δ^{k-1}` is the (k-1)-simplex on `U` and `H = G_N - U`. The
> Hodge spectrum of a join (Goff 2009; Horak-Jost 2013) gives
> contributions from each factor shifted by the other's vertex count;
> the eigenvalue `|V| = k + |V(H)|` arises from the **full simplex
> factor `Δ^{k-1}`** with multiplicity equal to the number of 0- and
> 1-faces of `Δ^{k-1}`, namely `k + C(k, 2) = C(k+1, 2)`. This matches
> empirical multiplicity exactly across all 9 tested N.

Why this is shape-revealing: it adds a **NEGATIVE-shape** edge in the
higher-order Hodge / arithmetic-flag-complex category, structurally
distinct from CLOSED 356/387/E7.12/E7.16 (all `L_0`-only). The L_1
spectrum's top-eigenvalue multiplicity is the FIRST Hodge fingerprint
in the project tying *higher-order* spectral data to a prime-counting
quantity (`π(N) − π(N/2)`). Composes with: triple-coprime singular
series (E2.13 family — the L_1 mean shift `mean(L_1, K_N) − mean(L_1,
ER) = 3(T_cp − T_ER) / |E|` reduces to `T_cp / T_ER → ∏_p (1 − 3/p² +
2/p³) / (6/π²)³ ≈ 1.27628`).

Auxiliary findings:

> - **β_k(K_N) = δ_{k,0} for all N ≥ 3.** Bertrand-prime universal
>   vertex makes K_N a cone, hence contractible. Hodge kernel is
>   deterministically trivial; the entire signal lives in the nonzero
>   spectrum.
> - **L_1 mean shift Z-score grows.** Z[mean(L_1, K_N) vs ER] grows
>   3.04σ → 4.40σ → 5.82σ → 9.48σ → 18.33σ across N = 32..128
>   (≈ √N scaling, from the |T| ~ N³ vs σ_T ~ N^{3/2} ratio).
> - **KS distance grows.** KS p-value drops from 1.9e-11 (N=32) to
>   < 1e-300 (N=128).

Composes with: E2.13 (Gowers / triple-coprime), E2.14 (Anderson chi_P),
E2.16 (DPP), E2.17 (PH on gaps), E2.19 (subword complexity), E7.12
(Cayley (Z/nZ)^* spectrum), E7.16 (prime-Cayley Friedman). Distinct
from CLOSED 387 (`L_0` only).

Cross-domain reference: Eckmann 1944 *Comment. Math. Helv.* 17, 240
(discrete Hodge theorem); Friedman 1996 *Algorithmica* 21, 331
(combinatorial Laplacian Betti algorithm); Horak-Jost 2013 *Adv. Math.*
244, 303 (join spectrum); Lim 2020 *SIAM Review* 62(3), 685 =
arXiv:1507.05379 (survey, primary cite); Goff 2009 *J. Algebraic
Combin.* 30, 215 (join-spectrum machinery); Kahle 2009/2014 (random
flag complexes). CROSS_DOMAIN_TECHNIQUES §1 row "Hodge / Laplacian on
simplicial complexes (higher-order L_k, k ≥ 1)" promoted PROPOSED (D22)
→ USED E with this edge.

### E7.18 — FHK-renormalised window max M_T = max log|ζ(1/2 + it)| − log log T + (3/4) log log log T has empirical mean −0.66 ± 0.05 that is T-independent across T ∈ {10⁴, 10⁵, 10⁶}, BUT distribution shape at finite T ≤ 10⁶ is approximately Gaussian (var ≈ 0.60 ≈ 1.47× FHK Gumbel(1/2) prediction π²/24, skew ≈ 0.10 vs Gumbel +1.14, excess kurt ≈ −0.4 vs Gumbel +2.4); FHK Gumbel limit shape is NOT finite-T-detectable [EVS M]

S133; ATTACK_VECTORS §C.C7;
`experiments/analytic/zeta_structure/fhk_amplitude_max/`;
CLOSED_PATHS row at session 133.

**FIRST ζ-amplitude (vs zero-position) measurement of the project**
— in a category structurally orthogonal to the 35+ existing
pseudorandomness measures and to the position-side family E7.1 / E1.10
/ E3.13 / E7.15.

Empirical FHK-mean-T-independence + FHK-shape-finite-T-undetectability:

> **Claim (verified at three T anchors with K = 100 unit windows
> per anchor, M = 200 evenly-spaced samples per window).** Define
>
>     M_T := max_{t ∈ [T, T+1]} log|ζ(1/2 + it)|
>             − log log T + (3/4) log log log T.
>
> (FHK 2012 PRL 108 conjectures: M_T → randomly shifted Gumbel(loc,
> scale = 1/2) as T → ∞; mean = μ_∞ + γ/2; variance = π²/24 ≈ 0.4112.)
>
> Empirically:
>
> | T_base | M_T mean ± sem | M_T var | KS to FHK Gumbel(1/2) | KS to free Gauss |
> |--------|----------------|---------|------------------------|-------------------|
> | 10⁴    | −0.699 ± 0.067 |  0.452  |        0.088           |      0.061        |
> | 10⁵    | −0.632 ± 0.083 |  0.692  |        0.169           |      0.062        |
> | 10⁶    | −0.641 ± 0.082 |  0.677  |        0.128           |      0.050        |
> | pooled | **−0.657 ± 0.045** |  0.604 |          —             |       —           |
>
> 1. **Mean T-INDEPENDENCE (FHK normalisation works):** pairwise
>    Z(M_T mean) ≤ 0.7 across all three pairs. Pooled mean
>    M_∞-mean = −0.657 ± 0.045 over 300 windows.
> 2. **Shape Gaussian-NOT-Gumbel (FHK shape NOT finite-T-detectable):**
>    KS to free Gaussian < KS to FHK Gumbel(1/2) at all 3 T (ratio
>    0.4-0.7); skewness 0.02-0.15 vs Gumbel +1.14; excess kurtosis
>    −0.85 to −0.14 vs Gumbel +2.4; Vuong z (Gauss vs Gumbel)
>    ∈ {−1.79, −1.43, −1.58} (joint Z ≈ −2.8, Gauss preferred).
> 3. **Variance 1.47× too large** vs FHK prediction 0.4112; bootstrap
>    95% CI [0.50, 0.85] at T = 10⁶ excludes the FHK value.

**Mechanism (closure mode I).** Saksman-Webb 2018 proved ζ(1/2 + it)
on a *mesoscopic* scale converges to a Gaussian multiplicative chaos
measure; FHK Gumbel-limit is a refined consequence applying log-
correlation theory to the asymptotic GMC. **The convergence rate at
finite T is NOT addressed in the published literature.** Empirically
the FHK *mean* convergence is FAST (T-stable at T = 10⁴ already) but
the *shape* convergence is SLOW (still ~Gaussian at T = 10⁶).
Plausible explanation: pre-freezing log-correlated noise is
approximately Gaussian (standard CLT-style argument on the scale-
summed log-correlation kernel); the freezing transition that produces
the heavy-tailed Gumbel structure has not yet activated at finite T.

> Why this is shape-revealing: complements the **position-side**
> measurements (E7.1 GUE-statistical zeros, E1.10 gap-shuffled
> control, E3.13 BK arithmetic correction below noise) with the FIRST
> **amplitude-side** project measurement. The amplitude SHAPE
> convergence to FHK is provably-empirically slow; this is the FIRST
> quantitative bound on the FHK convergence rate at finite T in any
> setting (project-internal or published).
>
> Quantitative content for future agents: the empirical FHK universal
> intercept M_∞-mean = −0.657 ± 0.045 ⇒ Gumbel-loc μ ≈ −0.946 (under
> the FHK Gumbel(loc, 1/2) form), giving GMC moment-generating
> constant `c ≈ 0.151` (under the FHK relation
> M_∞-mean = (γ + log c)/2). This empirical c is below the random-
> matrix-side prediction `c = π² G^4 / 4! ≈ 0.79` (Bourgade-Kuan 2014;
> G = Barnes G-function) by a factor 5×; the discrepancy is the
> finite-T finite-shape correction not yet absorbed.

> Composes with: E7.1 (zero-position GUE); E1.10 (gap-shuffled null
> for prime-frequency probes); E3.13 (BK below noise); E7.15
> (automorphic Hecke L(s, Δ) basis obstructed); pseudorandomness
> battery as the position-side complement to this amplitude-side
> measurement.
>
> Cross-domain references: **Fyodorov-Hiary-Keating 2012** *PRL* 108 =
> arXiv:1202.4713; **Saksman-Webb 2018** arXiv:1609.00027 (GMC limit
> of ζ on mesoscopic scale); **Arguin-Belius-Bourgade 2017** *CMP* 349
> = arXiv:1612.08575 (RMT-side Gumbel proof); **Bourgade-Kuan 2014**
> *CPAM* 67; **Selberg 1946** (Selberg CLT for log|ζ|).
> CROSS_DOMAIN_TECHNIQUES §3 row "Gaussian multiplicative chaos / FHK
> extreme-value statistics" promoted PROPOSED (C7) → USED I with
> this edge.

Successor open questions (S133 self-extension): (a) **C7.a mesoscopic
window** — repeat at window length δ = (log T)^{1/2}, where Saksman-
Webb proved sharp GMC convergence; does the Gumbel shape emerge at
this scale? (b) **C7.b joint argmax × prime alignment** — argmax
is uniform within sample noise (KS 0.16-0.20); does conditioning on
zero proximity reveal arithmetic structure? (c) **C7.c higher-order
Keating-Snaith joint moments** of (max, second_max) at λ < 1.

### E7.14 — Maynard 2015 multidimensional sieve weight is not a TC⁰ primality witness [EVS shape]

S116; ATTACK_VECTORS §A5;
`experiments/sieve/maynard_weight_pointwise/`;
CLOSED_PATHS row at session 116.

Two compositional obstructions show that Maynard's optimal k-dim
weighted Selberg sieve weight `w(n) = (Σ_{d_i|n+h_i, gcd=1, ∏d_i≤R}
μ(d_1)…μ(d_k) F(log d_i/log R))^2` is **not** a polylog primality test:

- **(1) Aggregate-not-pointwise.** Maynard's main theorem says
  `Σ_{N≤n<2N} w(n) χ_P(n+h_i) ≥ c(k) Σ w(n)` for some `i`. Empirically
  at single n: AUC restricted to odd n stays in `[0.66, 0.69]` across
  θ ∈ {0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40} and across four F
  choices. Best F1 = 0.62 with F = `(1-Σx_i)^2`. The "high AUC at
  low θ" signal (≤ 0.90) is parity detection (for H = {0, 2, 6} the
  three coordinates share parity), NOT sieve content; restricted to
  odd n it disappears.
- **(2) Divisor-enumeration cost.** Mean coprime tuple count scales as
  N^{0.10–0.12} for θ ∈ [0.20, 0.40]; worst-case tuple count grows
  proportionally. Listing squarefree divisors `d ≤ R = N^θ` of `n+h_i`
  requires Ω(R/log R) work without precomputation. The sub-routine
  reduces to growing-dim modular powering / divisor enumeration —
  exactly E5.3.

> Why this is shape-revealing: this closes the **most refined
> explicit prime-detection sieve in modern analytic number theory**
> (Maynard 2015, Polymath8b 2014, Pascadi 2025) as a candidate for
> polylog π(x) at single-n resolution. Combined with E6.7 (HKM time-
> space), E5.3 (MPOW), and E7.10 (AKS modulus-twist), Maynard
> sieve is now the FOURTH family closed in the explicit
> "construction-side" attack space: its information content is
> aggregate, not pointwise; its evaluation cost is sub-poly, not
> polylog. The pre-stated A-grade target (PRIMES ∈ TC⁰ via Maynard)
> is empirically falsified at single-n resolution and theoretically
> closed by reduction to E5.3.

Cross-domain reference: Maynard 2015 "Small gaps between primes"
arXiv:1311.4600; Polymath8b 2014 arXiv:1407.4897.

### E7.13 — Szegedy walks on arithmetic graphs do not yield polylog π(x) [EVS shape]

S80; ATTACK_VECTORS §D.D4;
`experiments/quantum/szegedy_walk_prime_graphs/`;
CLOSED_PATHS row at session 80.

For Szegedy quantum walks (Szegedy 2004, arxiv quant-ph/0401053) on
reversible Markov chains over `[1..x]` whose transition graph is
defined by a simple arithmetic relation, polylog primality extraction
requires *jointly* (i) classical spectral gap `δ = Ω(1/polylog(x))`
AND (ii) a single eigenvector of the discriminant matrix with
primality-localised mass. Three test families demonstrate (i) and (ii)
are **incompatible** in this setting:

- **Cayley `(Z/NZ)*`:** `δ → 0` as `N^{-c}`, c ≈ 0.87 — Szegedy
  mixing `~N^{0.41}` is poly, not polylog.
- **Coprime graph `C_x = ([1..x], gcd=1)`:** `δ = Ω(1)` (asymptotically
  `δ → 0.4166...`); only the trivial Perron eigenvector carries
  primality information (stationary mass on primes ≈ `(π²/6)/x`,
  exactly Mertens density `ζ(2)`).
- **Divisor graph `D_x`:** `δ = Ω(1)`; high-prime-mass eigenvectors
  exist but localise on prime-centered divisor clusters (one specific
  prime `p` + its multiples per eigenvector), not on global primality.
  Same degree-class probing mechanism as E7.12.

Verified at full sweep: 14 Cayley primes, 8 coprime sizes (x ≤ 1000),
8 divisor sizes (x ≤ 1000); cross-domain quadratic speedup
`t_Q ≈ sqrt(t_C)` empirically realised on the Cayley sweep.

> Why this is shape-revealing: this is the **quantum-walk family
> closure**. The Szegedy discriminant theorem cleanly converts the
> classical spectral structure of an arithmetic chain into a
> quantum-walk operator, but the conversion is gap-preserving: it can
> at best square-root the mixing time, never log-compress it. So any
> arithmetic graph whose classical chain mixes in poly(x) cannot mix
> in polylog under Szegedy walks; and any graph that mixes in O(1)
> classically already fails on the spectral-localisation side.
> Combined with E7.12 this rules out a substantial swath of
> "graph-based polylog primality test" attempts.

Cross-domain reference: Szegedy, M. *Quantum Speed-Up of Markov Chain
Based Algorithms*, FOCS 2004, arxiv quant-ph/0401053.

### E1.7 — Quantitative finite-x Wasserstein plateau for D(x)=(π(x)-Li(x))log(x)/√x  [EVS L]

S108; refined by S109/S110/S111 (CONFIRM); refined by S112 (PARTIAL);
refined by S113 (PARTIAL — plateau is universal across non-Gaussian
distributions, magnitude predicted by kurtosis alone);
refined by S114 (PARTIAL — S108 numeric reproduces to 0.3% via
three independent W_1 routines; S113 kurtosis-only fit refuted on
Beta(α,α) at matched kurt = -0.41, which gives W_1/σ = 0.0328 vs
S113-predicted 0.0426 vs D_emp 0.0376; universality is qualitative,
not tightly quantitative under kurtosis alone);
refined by S115 (CONFIRM — sub-window r=0.906 IS Riemann-φ-specific:
random-phase null on same Riemann γ_k gives r = -0.04 ± 0.39, 0/200
reach r ≥ 0.906; structural matching survives but content = E1.5).
`experiments/analytic/stein_wasserstein_pi/`;
`novel/finite_x_wasserstein_plateau.md`.

For K log-uniform anchors `x_k ∈ [X, eX]` with `X = 10^6, K ≤ 10⁴`, the
empirical Wasserstein-1 distance `W_1(D̂_K, N(μ̂(K), σ̂(K)²))` to the
sample-fitted Gaussian *plateaus* at a positive value
`c(X) ≈ 0.0083 (= 0.038 σ̂)`, while the i.i.d.-Gaussian-fluctuation
control follows the Stein-CLT rate `~ c_G(σ̂)/√K → 0`. Inflation
factor at K=10000: 4.5×, z-score = 15.34.

**Pointwise structural matching (CONFIRMED, S110):** at n=1000 zeros,
`corr(D_emp, D_th(n)) = 0.98` and `W_1(D_th)/W_1(D_emp) = 0.98`.
Sub-window correlation `r=0.906` (S108), `r=0.9154` (S110 disjoint
windows). These confirm D_emp ≈ low-zero explicit-formula sum
pointwise — i.e., re-confirm E1.5.

**W_1 magnitude is NOT Riemann-specific (S112 PARTIAL).** The
"indistinguishable from random-phase variant" test originally cited
as evidence for Riemann structural origin is a generic oscillatory-
sum property: at K=5000, n_modes=50, n_trials=60 per family,
D_emp's W_1 = 0.00863 is indistinguishable from
- random-phase Riemann (z=-0.93, mean=0.0120, std=0.0036);
- random-phase non-Riemann uniform [10,145] (z=-1.26, mean=0.0204,
  std=0.0094);
- random-phase non-Riemann equispaced (z=-1.55, mean=0.0175, std=0.0057).

Riemann ensemble vs non-Riemann ensemble IS distinguishable
(KS p<10⁻⁴ at n=60 trials), but a single empirical observation lacks
resolving power. So the W_1 *magnitude* is consistent with any low-
frequency oscillatory sum; the §C5 verbatim criterion clause
"the gap…ties to a specific zeta-zero contribution" fails.

**Excess kurtosis** = -0.41 (95% CI [-0.46, -0.36] at K=10000) —
sub-Gaussian, sourced from arcsine-distributed individual cosine
modes `cos(γ_k log x − arctan(2 γ_k))`. The sub-Gaussian signature
is window-width dependent (S111: kurt flips sign at narrow logw=0.5).

Why an edge: this is the FIRST quantitative finite-x Wasserstein-shape
bound for `π(x) - Li(x)`. The bound itself is a generic-oscillatory-
sum quantity that happens to be realised by D(x) via E1.5; it is not
a Riemann-specific quantitative statistic. EVS L (was M; demoted by
S112 because the quantitative magnitude is not arithmetic-specific).

**S113 universality strengthening:** the plateau is NOT EVEN
oscillatory-sum-specific. Across 9 distribution families (uniform,
single arcsine, sums of arcsines with various weights, two-Gaussian
mixture, t df=10, Laplace, analytic low-zero sum, Gaussian control),
*every* distribution with non-zero kurtosis plateaus at K=10000,
while pure Gaussian decays as 1/√K (ratio 6.42 ≈ √50, matching the
Stein-CLT rate). The plateau magnitude W_1/σ tracks |kurtosis(D)|
monotonically; linear interpolation at kurt=-0.41 (D_emp's kurtosis)
predicts W_1/σ ≈ 0.042, vs the observed 0.038 (within 10%). So the
W_1 plateau is a generic kurtosis-driven W_1(P, N(μ_P, σ_P²))
quantity — well-defined positive whenever P is not Gaussian — and
its magnitude is fully derivable from D's kurtosis under log-uniform
x. The cross-domain Stein technique import survives but is not
load-bearing for the conclusion.

Why **not** an attack route: the plateau reduces *quantitatively* to
the explicit-formula low-zero contribution (E1.5) — pointwise. It
joins the GUE-sieve-circuit closure family (E7.1, E7.6, E7.11). The
W_1 magnitude carries no information beyond "oscillatory sum on a
log-uniform grid"; the cross-domain Stein machinery did not surface
arithmetic structure orthogonal to E1.5.

Sub-window range dependence: at `x ∈ [10^7, 10^8]` (K=1000),
`c(10^7) ≈ 0.0067 < c(10^6) ≈ 0.0087` — consistent with asymptotic
Hejhal CLT prediction `c(X) → 0` as `X → ∞`.

Cross-domain reference: Chen-Goldstein-Shao 2011 *Normal Approximation
by Stein's Method* (Springer); Ross 2011 *Probability Surveys* 8, 210.
First application of Stein's method to π(x) - Li(x) in this project.

**S115 sub-window correlation IS Riemann-φ-specific (CONFIRM).** S110
reported `r=0.9154` on disjoint windows using ACTUAL zero phases.
S115 ran a 200-trial random-phase null preserving the same 50 Riemann
γ_k: mean r = -0.04 ± 0.39 (matches the n=10 noise-floor 1/√10), with
0/200 trials reaching r ≥ 0.906. Random-phase + non-Riemann γ ∈
[10, 145] and pure-noise controls give statistically identical
distributions (means within 0.02; stds within 0.01). The actual zero
phases ARE necessary and sufficient to reproduce V_emp pointwise
across sub-windows — exactly as E1.5 predicts. No grade change: the
B demotion was about W_1 *magnitude* generic-ness, not the pointwise
match. The structural-matching claim survives but its mathematical
content is exactly E1.5.

### E7.11 — Convergence-acceleration / variance-reduction family is exhausted [EVS shape]

Critique-46 synthesis (`archive/sessions/session_critique46.md`);
sessions S5, S6, S10, S15, S25, S32, S43-S46, S48, S51, S63, critique-46.

Across 13+ sessions the project has tested every standard
convergence-acceleration and variance-reduction transform applied to the
truncated explicit formula `psi_T(x) = sum_{|gamma|<=T}` (or to the
analogous Cipolla / Lagarias-Odlyzko zero sums). **None beats the raw
partial sum at x large enough for the rounding test.** Many are strictly
worse:

| Transform | Outcome | Worst measured ratio vs partial sum |
|-----------|---------|-------------------------------------|
| Padé / Wynn-eps (= Shanks) | strictly worse | **1900–4100x WORSE** (P-46 A at x in {100,1000,10000}) |
| Aitken Δ², Richardson, Cesàro | strictly worse | order-of-magnitude regressions |
| Fejér window | constant gain only | 67% recovery at K=100 vs 53% sharp (E3.7) |
| Borel (single-dir), Borel-Padé | mixed → equivalent | locks into T-independent asymptote (E3.5 + E3.7 follow-up) |
| Cesàro × Borel-Padé hybrid | strictly worse | 3/8 vs 6/8 Fejér-alone (line 711) |
| Borel-Padé applied to Cipolla for p(n) | strictly worse | **1.4–14x worse** than raw Cipolla (P-46 B) |
| Mellin-Barnes contour | equivalent | (line 49) |
| Hermite / Gaussian / Riesz mollification | strictly worse | 54x worse at x=100, K=800 (line 693) |
| Selberg Dirichlet-poly mollifier | strictly worse | best 1.72 vs sharp 0.44 at x=10000 (S63 P2, line 716) |
| PNT zero-aware control variate (MC) | trivial | **1.006x variance reduction** (P-46 D) |
| Harper random-multiplicative | equivalent | (line 257) |
| Random K-subset of zeros | strictly worse | 4x worse than first-K (S54, line 700) |
| Randomized I-E sieve | equivalent | (S54) |

The unifying mechanism is the **linear-functional bound** captured in
line 693:

> Any kernel `w(gamma)` with `sum |w(gamma)| >= c · N(T)` inherits the
> tail bound `Omega(sqrt(x) sqrt(log T / log x))` from the underlying
> explicit-formula tail. Re-weighting cannot beat this — only NONLINEAR
> operations on zeros could.

Borel-Padé and the hybrid stacks fail by a related but distinct
mechanism: the `1/k!` Borel transform exponentially suppresses
information from zeros beyond the optimal-truncation point, locking the
estimate into a leading-term asymptote that doesn't move with K.

> Why this is shape-revealing: it is a **search-space constraint on the
> entire family of linear post-processing of zero sums.** Any future
> proposal in this family — re-weighting, smoothing, mollifying,
> transforming, control-variate — is provably orthogonal to the sqrt(x)
> wall. The only escape would be a *nonlinear* operation on zeros (e.g.,
> a non-trivial functional that uses products / quotients / max-plus over
> the zero spectrum, not weighted sums). No such operation has been
> identified; E7.1 (zeros are linearly independent over Q) and E7.4
> (inversion non-contracting) further constrain what nonlinear operations
> could even leverage.

Pairs with E7.10 (AKS modulus-twist orthogonality) as the second major
*family-level closure* — between them, two of the largest unexplored
attack families have been systematically exhausted, in the sieve+analytic
and TC⁰-primality directions respectively.

### E7.7 — Three-pillars meta-theorem (informational closure)        [EVS shape]

S15, S16; `archive/sessions/session16_synthesis.md`;
`proven/circuit_size_barrier.md`;
**algorithmic projection (Pareto map):** `experiments/constructions/pillar_tradeoff_diagram/` (S81, C6).

Across **15+ candidate intermediate-quantity families** systematically
tested (class numbers h(-d), L-values L(1,chi), elliptic curve a_p,
regulators, additive sumsets, ergodic theory, model theory, tropical
geometry, sufficient statistics, F_q point counting, S_n/GL_n
representation theory, and 4 more), every one routes back to one of:

```
   prime positions  |  zeta zeros  |  floor values {floor(x/k)}
```

These three are the **only known informationally complete encodings**
of pi(x), and they are pairwise informationally equivalent:
prime positions ↔ explicit formula ↔ Möbius/Lucy DP each costs
O(x^{1/2 .. 2/3}) to convert.

> Why this is shape-revealing: it is a **search-space constraint** on
> all future chain construction. Any new chain that doesn't reduce to
> one of these three encodings is suspect of either circularity (it's
> secretly going through primes), equivalence (it's secretly the
> explicit formula or sieve), or information loss. FOCUS-6 in TODO.md
> is the explicit project to find a fourth pillar.

---

## 8. Putative chain templates

A chain is a sequence of edges that together produce a polylog
architecture for p(n). Each surviving chain has at least one **missing
primitive** (clearly labelled), without which it would close. The
purpose of catalogueing them is to direct future work at the missing
primitives, not to claim a result.

> **Chains A, B, C, D, G have been removed.** Each was empirically
> closed by S52-S57 (Connes operator scaling, character-twisted
> Liouville pseudorandomness, free-identity vacuity, Borel-Fejer
> regression, per-bit reduction to A+D). See
> `status/CLOSED_PATHS.md` and `archive/sessions/session5{2,3,5,6,7}_*.md`
> for full closure evidence. The two remaining chains below are the
> only ones with *any* uncovered attack surface.

### Chain E — TC^0 via growing-dimension MPOW                        [EVS H]

Primitives:
1. E5.3: PRIMES in TC^0 ↔ growing-dim MPOW in TC^0.
2. E5.1 / E5.2: BPSW or QFT in TC^0 *bypasses* the growing-dim issue
   for primality, but not for counting.
3. **MISSING**: either resolve growing-dim MPOW in TC^0 (positive),
   or find an alternative primitive that turns sieve-style counting
   into TC^0.

If primes are in TC^0, sum_{k <= x} chi_P(k) is an `x`-input MAJORITY-
style threshold counting circuit of size poly(x), depth O(1). To
reach polylog *time*, we further need #TC^0 ⊆ NC, which is not known
either way. **This is the deepest open frontier of the project.**

S47 closed the cyclotomic-CRT sub-attack (AKS-prescribed r is prime in
21/22 sampled cases). FOCUS-5 in `TODO.md` lists three untested
sub-attacks (smaller-r workaround, non-cyclotomic ring decompositions,
partial-Frobenius transplant) — these are where any new work belongs.

### Chain F — Aggarwal-style binary search amplification             [EVS L]

Primitives:
1. E6.6: p(n) reduces to O(log p(n)) calls to pi(·).
2. E6.5: pi at one point in O~(sqrt x) elementarily (HKM).
3. **MISSING**: any one polylog speedup at *any* step of the chain
   — none currently known.

This is the "shortest path from current state-of-art to polylog," but
each step needs an independent breakthrough. Listed for completeness;
contains no fresh attack surface beyond what Chain E already covers.

---

## 9. Tactical primitives ready for re-use in chains

| Primitive | Cost | Source |
|-----------|------|--------|
| R^{-1}(n) — Newton inversion of R(x), 50% of digits | O(polylog) | E6.1, V10 algorithm |
| Lambert-W Pade approximation, single closed form, 0.019% error | O(1) | `algorithms/v1_pade_approximation.py` |
| BPSW primality test (verified to 2^64) | O(polylog) | E5.1 |
| Strong Lucas (2x2 MPOW), Jacobi (TC^0) | O(polylog) | E5.1 |
| Liouville parity identity verified to x = 10^4 | exact | E2.2 |
| Connes operator giving 50 zeros from primes ≤ 13 to 10^{-55..-3} (one-shot fit, not algorithm; see S53) | poly(prime budget) | E3.1 |
| Borel-Padé regularisation of explicit formula | constant-factor | E3.5 |
| Cesaro-Fejer kernel: 67% recovery at K=100 (vs 53% sharp) | constant-factor | E3.7 |
| MPS-base-W identity: rank = min(W^j, phi(W) W^{d-j-1} + 1) | exact | E2.1 |
| Wheel-30 + prefetch on Lucy DP | 2-35% | E4.5 |
| HKM elementary O~(sqrt x) for pi(x) | O(x^{1/2+eps}) | E6.5 |
| Aggarwal binary search: p(n) from pi(x) calls | O(log x) calls | E6.6 |
| Bug fix: use mpmath.ei(rho * log x), never li(x**rho) | numerical | E3.6 |
| pi(x) mod 2 sign predictable by Chebyshev bias 99.5% | O(1) bits | E3.9 |

---

## 10. Single-target objectives ranked by "if you cracked this, the
   rest of the chain is trivial"

Only items with *uncovered* attack surface are listed. Objectives whose
last attack sub-route was structurally closed in S52-S57 have been
removed.

1. **Growing-dim MPOW in TC^0** (Chain E). — Pure complexity-theory
   open problem. Cyclotomic-CRT sub-attack closed S47; three untested
   sub-attacks (smaller-r, non-cyclotomic ring, partial-Frobenius) are
   listed in `TODO.md` FOCUS-5.
2. **A non-AKS TC^0 primality test** (`status/OPEN_PROBLEMS.md`). —
   Eliminates the AKS bottleneck of E5.3 if found. No identified
   attack surface beyond Chain E sub-attacks.
3. **A 4th encoding of pi(x)** (S16, FOCUS-6). — Three known: prime
   positions, zeta zeros, floor values. A fourth, info-equivalent
   encoding could sit in TC^0. 15+ candidate families closed S15/S16
   but the search was not exhaustive.

---

## 11. Reading-list of supporting documents

- **The big picture:** `novel/info_computation_gap.md`,
  `novel/failure_taxonomy.md`, `novel/pseudorandomness_of_pi.md`.
- **delta(n) deep structure (E1.7, E1.8, E2.11):**
  `experiments/information_theory/kt_complexity/SYNTHESIS.md`,
  `experiments/algebraic/identity_search/RESULTS_SUMMARY.md`.
- **HKM time-space curve & Aggarwal-Dusart bracket (E6.7, E6.8):**
  `literature/aggarwal_2025_analysis.md`.
- **Engineering ceiling of sieve route (E6.9):**
  `literature/primecount_analysis.md`.
- **What's proven hard:** `proven/circuit_size_barrier.md`,
  `proven/uniformity_barrier.md`,
  `proven/convergence_acceleration_barrier.md`,
  `proven/workspace_mismatch_barrier.md`.
- **The only open frontier:** `status/OPEN_PROBLEMS.md`,
  `novel/tc0_matrix_powering_reduction.md`,
  `novel/bpsw_tc0_reduction.md`,
  `literature/matrix_powering_tc0.md`.
- **The bisection at N/2:** `novel/approx_degree_prime.md`,
  `novel/carry_propagation_boundary.md`, `novel/delta_spectrum.md`.
- **The MPS theorem:** `novel/mps_bond_dimension.md`.
- **The Liouville identity:** `experiments/proposals/proposal25_liouville_parity_check.py` and `_results.md`.
- **The Connes spectral triple:** `literature/state_of_art_2026.md` §2.5b.
- **Bug to avoid:** `experiments/proposals/proposal26_fejer_damped_explicit.py` (mpmath.ei vs li branch cut).

---

## 12. Open composition spaces

Edges are catalogued individually (§1-§7); their *interactions* are
mostly unexplored. Compositions are the cheapest source of genuinely
new mathematical objects at the project's current maturity. See
`NOVELTY_CHALLENGES.md` §1 for the active composition challenges.

The pair-space `EDGES.md × EDGES.md` is large (~68² ≈ 4600 entries);
most pairs are uninteresting (e.g. composing two analytic edges from §3
typically just stacks errors). The interesting compositions cross
sections — algebraic × analytic, information × computational,
negative-shape × positive — and use complementary structure.

### Cross-section pair grid (which composition has been built?)

| Pair | Status | Where |
|------|--------|-------|
| E1.6 (A⊕C₃ bisection) × E1.5 (0.537-bit invariant) | OPEN | C1 |
| E2.1 (MPS bond-dim) × free-probability | OPEN | C2 |
| E5.8 (Brandt obstructions) × E1.3 (per-bit difficulty) | OPEN | C3 |
| E6.6 (Aggarwal binary search) × E6.8 (Dusart bracket) × E5.1 (BPSW) | OPEN | C4 |
| E1.4 (N/2 universality) × E2.5 (multilinear poly) on non-π | BUILT (S71) | C5 — refined E1.4, unified E2.7+E2.8 by column-zero density |
| E7.7 (three pillars) × E6.7 (HKM tradeoff) | OPEN | C6 |

### Triple-edge compositions worth scoping

- **E1.5 + E1.6 + E2.10** — three "modular structure" edges. Conjecture: the
  free identity L mod 2 = x mod 2 (E2.10) extends to a *family* of free
  identities indexed by the bisection (E1.6) at rate 0.537 bits/step (E1.5).
  If true, this would be a single statement subsuming three edges.
- **E2.1 + E2.7 + E2.8** — three "structural rank" edges (MPS bond-dim,
  communication rank +2, tensor rank 25-35% below random). The ratios
  φ(W)/W (E2.1), 1+2/2^{N/2} (E2.7), 0.65-0.75 (E2.8) might be facets of
  a single underlying rank-deficiency invariant. Worth a unified
  statement attempt.
- **E5.5 + E5.6 + E5.7** — three "TC⁰ frontier" edges (formula bound
  2^{N/2-O(1)}, nonuniform TC⁰ unconditional, #TC⁰ ⊆ NC open). Composition:
  is there a function in nonuniform TC⁰ that requires formula-size
  2^{N/2-O(1)} AND whose counting class falls outside NC? If so, π(x)
  is consistent with this spec — but does the spec FORCE π(x) outside
  uniform TC⁰?

### Unexplored cross-section frontier

These edges have NEVER been composed with any other edge:

- E2.4 Cloitre 2025 effective analytic recurrence (× anything)
- E3.10 Cully-Hugill-Lee improved error term (× any algorithmic edge)
- E4.2 Lucy DP transition matrices unipotent (× spectral edges)
- E4.7 Prime parity stream h_μ deficit 0.020 nats (× pseudorandomness battery)
- E6.4 90/10 split in sieve update spectrum (× compression edges)

Each of these is a one-session composition target.

---

*Compiled: 2026-04-25, scanning Sessions 1-49 (690+ closed paths).
Extended: 2026-04-26 deep file-by-file review — added E1.7, E1.8, E2.11,
E2.12, E3.12, E6.7, E6.8, E6.9, E7.8, E7.9 plus the finite-size DCT
caveat to E6.3.
Extended: 2026-04-26 post-S67 — E2.11 pre-test closure of FOCUS-2 (9
candidates, mode I, 8 NEW); methodological note added under E2.11 about
its use as a fourth-encoding pre-disqualifier.
Extended: 2026-04-26 post-S61..S66 (FOCUS-1 sub-attack closures) — added
E1.9 (phi 2D rank, 22nd pseudorandomness measure), E7.10 (AKS modulus-
twist orthogonality theorem, the structural meta-finding from S61/S64/
S66), K_min non-monotonicity caveat to E3.11, AKS-gcd-extraction
folklore note under E5.3.
Extended: 2026-04-26 post-critique-46 — added E7.11 (convergence-
acceleration / variance-reduction family systematically exhausted across
13+ sessions; pairs with E7.10 as the second family-level closure).
Extended: 2026-04-26 post-S49 (FOCUS-4 closure) — added E1.10 (gap-
shuffled-zeros null methodology, mandatory for prime-frequency Fourier
probes of pair-correlation residuals).
Extended: 2026-04-26 audit (post-S49 content gap) — added E3.13 (BK
arithmetic correction empirically absent at N=8000, z=-10.85σ) as the
content companion to E1.10's methodology entry; refined E2.6 with the
S48-fresh forward-difference signature as the fourth independent
witness of the 2^{N/2} wall.
Extended: 2026-04-26 post-S51 (FOCUS-3 closure) — added E5.8 (Brandt 2024
TCC MKtP-diagonalisation is structurally welded to MKtP and does not
extend to natural functions like π(x) mod 2). With E7.10 + E5.8,
**Chain E is closed for both known technique families** (AKS-style and
diagonalisation-via-meta-complexity).
Extended: 2026-04-26 post-S67 mature-state rebalance — added §12 "Open
composition spaces" listing the unexplored cross-section edge pairs as
the new high-leverage direction. Active research targets are now in
`NOVELTY_CHALLENGES.md` and long-horizon arcs in `RESEARCH_AGENDA.md`.
This document is alive — extend it with new edges as they appear.*
