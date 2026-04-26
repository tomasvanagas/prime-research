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

### E5.2 — QFT (Grantham) in TC^0                                    [EVS M]

S13. Quadratic Frobenius test is in TC^0 (operates in 2-D algebra =
2x2 MPOW). Error < 1/7710 per parameter set. Deterministic correctness
under GRH.

> Edge value: a second TC^0 primality oracle, redundant with E5.1 but
> useful if BPSW counterexample ever appears. Independent assumption.

### E5.3 — Growing-dim MPOW in TC^0 = the only OPEN frontier         [EVS H]

S11, S39, S47; `literature/matrix_powering_tc0.md`;
`status/CLOSED_PATHS.md` line 232.

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
Mendlovic).

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

---

## 7. Negative-but-shape-revealing edges

These are not "edges" in the affirmative sense, but they constrain the chain.

### E7.1 — Each zeta zero contributes "independent" information     [EVS shape]

S25, S45, S20, S57; `experiments/wildcard/zero_sum_acceleration.py`,
`experiments/analytic/zeta_structure/triple_correlation.py`.

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
factor-10⁴ separation, fully consistent with GUE rigidity. **No higher-
order cluster structure exists at any tested order.**

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
`proven/circuit_size_barrier.md`.

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
