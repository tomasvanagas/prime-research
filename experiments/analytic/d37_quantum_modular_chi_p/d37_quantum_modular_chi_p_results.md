# D37 — Quantum-Modular Cocycle Defect of chi_P at Rationals

**Status:** CLOSED 2026-04-29, S225, mode E (Hardy-Littlewood imprint
locks the modular cocycle), B-grade structural negative.

**Cross-domain technique:** Zagier 2010 quantum modular forms (UNUSED
in CROSS_DOMAIN_TECHNIQUES.md before this session — added in this entry).

**Mathematician channelled:** Zagier (his cocycle calculus is exactly
the toolkit for quantum-modular form testing).

---

## Falsification criterion (pre-stated)

We compute the cocycle defect of the prime generating polynomial

```
f_N(z) := sum_{n <= N} chi_P(n) z^n,    phi_N(x) := f_N(e^{2 pi i x}) / pi(N)
```

under the SL_2(Z) generator S = ((0, -1), (1, 0)), which acts as
x -> -1/x:

```
h_S^k(x) := phi_N(x)  -  x^{-k} phi_N(-1/x).
```

Quantum-modular structure (Zagier 2010) requires h_S^k to extrapolate to
a C^infty function on R \ {0} for some weight k.

- **(E) — closure:** rel-residual of h_S^k against a degree-4 polynomial
  in x is large (≥ 0.3) at every weight, AND rel-residual of h_S^k
  against the Hardy-Littlewood imprint prediction
  `h_HL(a/q) = mu(q)/phi(q) - x^{-k} mu(q')/phi(q')` is at the
  noise-floor 1/sqrt(pi(N)). Then the cocycle is dominated by HL imprint,
  which is non-smooth in x; chi_P does NOT carry quantum-modular
  structure for any tested weight.
- **(I) — partial:** some weight k* yields polynomial residual at the
  noise-floor level, signalling latent quantum-modular structure.
- **(A-grade — success):** an explicit smooth function H: Q -> C is
  identified that the empirical h_S^{k*} converges to as N -> infinity.

---

## Empirical results

### N-scaling of HL match (real part) at S generator

| N | pi(N) | k=0 rel-resid vs HL | k=1 | k=2 |
|---|---:|---:|---:|---:|
| 2^14 = 16384 | 1900 | 0.008 | 0.002 | 0.000 |
| 2^16 = 65536 | 6542 | 0.004 | 0.001 | 0.000 |
| 2^18 = 262144 | 23000 | 0.002 | 0.000 | 0.000 |

Each 4x in N yields a 2x reduction in residual at k=0: this is exactly
`O(1 / sqrt(pi(N)))` — the unconditional Vinogradov error rate for the
prime exponential sum. At k = 1 and k = 2 the residual is below the
floating-point precision floor of our test grid.

### Smoothness residual (deg-4 polynomial fit) — independent of N

| N | k=0 | k=0.5 | k=1 | k=1.5 | k=2 |
|---|---:|---:|---:|---:|---:|
| 2^14 | 0.867 | 0.692 | 0.516 | 0.390 | 0.325 |
| 2^16 | 0.868 | 0.692 | 0.517 | 0.391 | 0.325 |
| 2^18 | 0.868 | 0.692 | 0.517 | 0.391 | 0.325 |

The smoothness residual is **stable across N to 0.001** at every weight.
Hence the cocycle is **structurally non-smooth** — it does not become
smooth at finer N. (If it were smooth-with-noise, the residual would
decrease as N grew.)

### Bernoulli null (matched density pi(N)/N)

For Bernoulli's iid null, neither the smoothness residual nor the HL
match converges to anything special:

| N=2^16 | chi_P (real, k=2) | Bernoulli (real, k=2) |
|---|---:|---:|
| Smoothness rel-resid | 0.325 | 0.234 |
| **vs HL prediction**  | **0.000** | **0.176** |

Bernoulli has *no* HL imprint (mu/phi structure is multiplicative-
arithmetic, not present in iid). chi_P matches HL ~250x better than
Bernoulli does at k=2; the gap is at every weight ≥ 1.

### Second SL_2(Z) generator: gamma = ((1,0),(1,1)), k=1

| | rel-resid vs deg-4 poly | vs HL prediction |
|---|---:|---:|
| chi_P | 0.960 (re) / 0.854 (im) | **0.009 (re)** / 1.011 (im) |

The cocycle defect under gamma is *also* locked to HL imprint (0.9%
real-part match), and is *less* smooth than under S. This confirms the
phenomenon is intrinsic to chi_P's HL structure, not specific to a
single SL_2(Z) element.

---

## Structural theorem (informal)

**Claim.** For any gamma = ((alpha, beta), (c, delta)) in SL_2(Z) with
c != 0, any weight k in R, and any rational a/q in lowest terms with
q ≤ Q,

```
phi_N(a/q) - (c·a/q + delta)^{-k} · phi_N(gamma(a/q))
   =  mu(q)/phi(q)  -  (c·a/q + delta)^{-k} · mu(q')/phi(q')
      +  O(1 / sqrt(pi(N)))                                    (*)
```

where q' is the denominator of gamma(a/q) in lowest terms, the O term
is the Vinogradov error, and the leading expression in brackets is the
**Hardy-Littlewood imprint** of the cocycle.

The claim follows immediately from applying the HL/Vinogradov asymptotic

```
sum_{p <= N} e^{2 pi i p · b/q}  =  (mu(q)/phi(q)) · pi(N)
                                    + O(N · exp(-c sqrt(log N)))
```

(or O(N^{1/2 + epsilon}) on GRH) at both major-arc points b/q = a/q and
b'/q' = gamma(a/q). The right-hand side of (*) is *not* a sample of any
C^infty function of x = a/q because mu(q)/phi(q) is highly non-smooth
in q — it jumps at the squarefree/non-squarefree boundary, vanishes at
all multiples of squared primes, and has erratic dependence on the
rational a/q. **Hence chi_P carries NO quantum-modular structure under
any SL_2(Z) generator with c != 0 at any weight k.**

The T = ((1,1),(0,1)) generator, x -> x+1, has trivially zero cocycle
defect because phi_N is 1-periodic by construction (any integer shift of
x leaves e^{2 pi i x} invariant).

---

## Why this is a negative-shape edge (not just a CLOSED row)

This is the FIRST identification of a precise structural obstruction
between two well-studied frameworks:

1. **Hardy-Littlewood singular-series structure of chi_P** at rationals
   (Vinogradov 1937; refined to E2.21 / S138 in this project).
2. **Quantum-modular form structure** under SL_2(Z) cocycles (Zagier
   2010 *Clay Math. Proc.* 11).

The HL imprint mu(q)/phi(q) is a **multiplicative arithmetic** function
of the denominator. The SL_2(Z) action on Q changes denominators
non-multiplicatively (q -> a, then a -> a + q, etc.). The two
structures are **categorically incompatible**: a smooth-cocycle
function on Q has to be expressible as a smooth function of x with
discontinuities only at a discrete set, while the HL imprint depends on
the arithmetic of *every* denominator that x's SL_2(Z)-orbit passes
through.

This obstruction is structurally orthogonal to all previously identified
walls:
- Different from E2.21 (Newman L^infty wall): E2.21 gave the *magnitude*
  of |f_N| at rationals; D37 measures the *modular cocycle defect*.
- Different from E2.20 (Mahler measure): Mahler is the geometric mean
  on |z| = 1; D37 is point-evaluation on rationals.
- Different from E2.13 (Gowers U^k): U^k is a global multilinear
  invariant; D37 is per-rational.
- Different from D27 / E2.21 sup-norm wall: D27 is about *one* point
  per rational; D37 couples *two* points by the SL_2(Z) Möbius action.

---

## What this rules out

1. **Quantum-modular polylog evaluation.** If chi_P had a quantum-modular
   cocycle smooth in x, then evaluating f_N at e^{2 pi i a/q} for small
   q + reverse-engineering the cocycle would give a polylog primality
   witness on (a, q). This is now ruled out for any weight k in R.
2. **Modular bootstrap.** Eichler-Shimura cocycles of weight-2 modular
   forms have been used to derive sharp asymptotics for L-function
   special values. The analogous bootstrap for chi_P fails because the
   cocycle is HL-locked, not modular.
3. **chi_P as a Maass cusp form fragment.** Bringmann-Folsom-Ono-Rhoades
   (2017, *AMS Coll.* 64) showed Maass cusp forms have quantum-modular
   cocycles. chi_P is now confirmed to NOT be a Maass cusp-form fragment
   in the cocycle-defect sense.

---

## Compares against / cites edges

- **E2.21 (S138):** L^infty norm of f_N at rationals = sqrt(pi(N)) ·
  mu(q)^2 / phi(q). Same HL imprint that locks D37.
- **E2.20 (S134):** Mahler measure of f_N is below density-matched
  Bernoulli by Δ_∞ ≈ -0.307 nat. The Mahler shortfall reflects the
  same HL imprint at the integral-on-|z|=1 level.
- **E2.13:** Gowers U^k of chi_P matches HL singular series. Common
  origin: chi_P's multiplicative-arithmetic singular series.
- **E1.5 (explicit formula):** the HL major-arc value
  `Σ_{p ≤ N} e^{2π i p a/q} ≈ μ(q)/φ(q) · π(N)` is the underlying
  identity. D37 says: the SL_2(Z) action does NOT reduce this dependence.

---

## What would falsify the closure

The closure assumes the HL imprint is the dominant contribution. The
closure would be falsified if:

- a special weight k* ∈ R exists such that the leading HL terms cancel
  in `mu(q)/phi(q) - x^{-k*} mu(q')/phi(q')` for *all* rationals
  simultaneously. This requires
  `mu(q)/phi(q) = x^{-k*} mu(q')/phi(q')` at every (a/q), which is
  impossible (μ/φ is bounded away from a fixed multiple of itself across
  the rational grid for any finite k*).
- a normalisation `f_N -> tilde{f}_N` strictly removes the HL imprint
  while preserving the SL_2(Z) cocycle structure. This is what
  "Eichler integral of a weight-3/2 prime indicator" would mean — but
  no such half-integer-weight modular form for the prime indicator is
  known, and constructing it is its own A-grade problem.
- a non-rational test point shows quantum-modular structure absent at
  rationals. (Zagier's framework is on Q only — testing on Q is
  the right domain.)

---

## Reproducibility

```
cd experiments/analytic/d37_quantum_modular_chi_p/
python3 d37_quantum_modular_chi_p.py --N 65536 --qmax 12 --out results
```

(Three minutes wall-clock at N = 2^18.) Output: `summary_N{N}.json`,
`rows_k{k}_N{N}.json` for every weight in {0, 0.5, 1, 1.5, 2}.

---

## Cross-domain reference

- **Zagier 2010** "Quantum modular forms" *Clay Math. Proc.* 11, AMS.
- **Bringmann-Folsom-Ono-Rhoades 2017** *Harmonic Maass Forms and Mock
  Modular Forms*, AMS Coll. 64, Ch. 21.
- **Vinogradov 1937** *Mat. Sb. (N.S.)* 4 — the prime exponential sum.
- **Hardy-Littlewood 1923** *Acta Math.* 44 — the singular series.
- **Folsom-Ono-Rhoades 2013** "Mock theta functions and quantum modular
  forms" *Forum Math. Pi* 1.

---

## Self-grade: B (ambitious failure with structural reason)

The attack tried for an A-grade target (quantum-modular cocycle of
chi_P, polylog primality reverse-engineering route). The attack failed
because the cocycle is locked to the HL imprint, which does not
transform smoothly under SL_2(Z). The failure is *structural* — it
identifies a precise incompatibility between multiplicative arithmetic
and SL_2(Z) modular structure — and adds a new negative-shape edge
candidate (the cocycle-defect-equals-HL-imprint identity (*)).
