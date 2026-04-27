# D22 — Higher-order Hodge Laplacian on the coprimality flag complex

**Session 126 (2026-04-27).**
**Cross-domain attack target:** §D.D22 (ATTACK_VECTORS.md).
**Cross-domain technique:** Combinatorial Hodge Laplacian on simplicial
complexes (Eckmann 1944; Friedman 1996; Horak-Jost 2013; Lim 2020).
**Cross-domain status:** marked PROPOSED in `CROSS_DOMAIN_TECHNIQUES.md`
§1 ("Hodge / Laplacian on simplicial complexes (higher-order L_k, k ≥ 1)").
This is the first quantitative use; promote PROPOSED → USED-E.
**Channelled mathematician:** Friedman (combinatorial Laplacians /
simplicial-complex Betti via spectral methods).

## What is computed

Build the **coprimality flag complex** `K_N := \{σ ⊆ [2..N] : σ is
pairwise-coprime\}` (a `(k+1)`-simplex is a pairwise-coprime
`(k+2)`-tuple). Compute boundary matrices `B_1` (edges → vertices),
`B_2` (triangles → edges), `B_3` (tetrahedra → triangles), and the
combinatorial Hodge Laplacians

    L_0 = B_1 B_1^T                       (graph Laplacian)
    L_1 = B_1^T B_1 + B_2 B_2^T            (edge Laplacian)
    L_2 = B_2^T B_2 + B_3 B_3^T            (triangle Laplacian, where feasible).

Hodge decomposition: `dim ker(L_k) = β_k(K_N, R)` (Eckmann 1944).
**Crucially, `L_1` and `L_2` have never previously been computed for
arithmetic flag complexes** — CLOSED_PATHS lines 356, 387 only address
`L_0` (Cayley / GCD-graph Laplacian = Ramanujan sums).

Compare against a matched-density Erdős-Rényi flag complex `Y(n, p)` with
`n = N - 1`, `p = empirical edge density of G_N` (≈ 6/π² + finite-N).

## Setup

| field | value |
|---|---|
| `N`s tested | `{8, 12, 16, 24, 32, 48, 64, 96, 128}` |
| seeds per N | 30 (N=32, 48, 64), 15 (N=96, 128), 5 (N≤24) |
| spectral routine | `numpy.linalg.eigvalsh` on dense `|E|×|E|` matrix |
| dense threshold | up to 4894 (N=128 `\|E\|`) |
| KS test | scipy.stats.ks_2samp on coprime vs pooled-ER `L_1` spectra |

Code: `hodge_coprimality_flag.py`. Raw data: `pilot.json`, `pilot32.json`,
`small.json`, `main_N32_64.json`, `main_N96_128.json`.

## Pre-stated falsifiers

The §D22 ATTACK_VECTORS entry pre-stated three failure profiles:
- **(E)** spectra match within KS noise floor; β_k within 2σ of ER;
  closes as 41st pseudorandomness category. **B-grade.**
- **(I)** spectra match BUT β_1 deviates by > 3σ — Möbius-singular-
  series fingerprint at edge level. **B-grade.**
- **(A)** L_1 has a uniform spectral gap λ_2^{(1)} > c > 0 in N (or
  reverse for ER); persistent Hodge spectral gap encodes new
  arithmetic constraint. **A-grade IF combined with TC⁰-computable
  boundary matrix.**

## Result summary

```
   N  |V|    |E|    |T|   k  C(k+1,2)  mult  Lmean_cp  Lmean_ER  Z_mean   Z_T     KS         KSp
  32   31    292   1376   5        15    15   16.137    13.209    3.04   2.11  0.210   1.92e-11
  48   47    664   4721   6        21    21   23.330    19.106    4.40   2.89  0.205   3.27e-24
  64   63   1196  11418   7        28    28   30.641    24.808    5.82   3.92  0.229   9.35e-54
  96   95   2710  38955   9        45    45   45.124    36.179    9.48   6.45  0.288  1.30e-187
 128  127   4894  96412  13        91    91   61.100    48.544   18.33  12.34  0.326   < 1e-300
```

`k = π(N) - π(N/2)` is the count of **Bertrand primes** in `(N/2, N]`.
`mult` is the multiplicity of `λ = |V|` in the empirical `L_1` spectrum.
`Lmean_cp - Lmean_ER` reduces (analytically) to `3 (T_cp - T_ER) / |E|`.

## Findings

### F1. Hodge KERNEL is deterministically trivial.

For all `N ≥ 3`, `β_0(K_N) = 1` and `β_k(K_N) = 0` for `k ≥ 1`.

**Why:** by Bertrand's postulate, there is a prime `p ∈ (N/2, N]`. Since
`p > N/2 ≥ |v|/2` for every `v ∈ [2, N]`, `p` does not divide any other
`v`, hence `gcd(p, v) = 1` for all `v ≠ p`. So `p` is a **universal
vertex** in the coprimality graph, and the flag complex `K_N` is the cone
on the rest, hence contractible — collapsing all `β_k` for `k ≥ 1`.

This forces **complete absence of harmonic forms in the Hodge spectrum**.
Every nonzero eigenvalue carries the structural information.

### F2. λ_max(L_1) saturates the universal bound `|V|`.

For all 9 tested `N`, `max(spec L_1(K_N)) = |V| = N - 1` to numerical
precision (≤ 1e-13 absolute residual).

The matched-density ER controls have `λ_max(L_1) ≈ p \cdot |V| + O(\sqrt{|V|})`,
strictly less than `|V|`. Empirical: across 30 ER seeds at N=64,
`max λ_max(L_1) = 51.7 < 63 = |V|_cp`.

### F3. **Multiplicity identity (main empirical theorem).**

> **Conjecture (verified at 9 values of N spanning 16-fold range)**:
>
> `mult(λ_max(L_1, K_N) = |V|) = C(k+1, 2) = k(k+1)/2`,
>
> where `k = π(N) - π(N/2)` is the Bertrand-prime count in `(N/2, N]`.

| N    | k  | C(k+1, 2) | empirical mult | match |
|------|----|-----------|----------------|-------|
| 8    | 2  | 3         | 3              | ✓     |
| 12   | 2  | 3         | 3              | ✓     |
| 16   | 2  | 3         | 3              | ✓     |
| 24   | 4  | 10        | 10             | ✓     |
| 32   | 5  | 15        | 15             | ✓     |
| 48   | 6  | 21        | 21             | ✓     |
| 64   | 7  | 28        | 28             | ✓     |
| 96   | 9  | 45        | 45             | ✓     |
| 128  | 13 | 91        | 91             | ✓     |

**Mechanism (proof sketch):** since every Bertrand prime `p ∈ (N/2, N]`
is universal in `G_N` and Bertrand primes are pairwise coprime (distinct
primes), the set `U` of Bertrand primes induces `K_k` in `G_N`. Moreover
`U` is in the *full join* with the rest: every non-Bertrand-prime vertex
is adjacent to every `p ∈ U`. Hence the flag complex factors as a join
`K_N = Δ^{k-1} \ast F(H)`, where `Δ^{k-1}` is the `(k-1)`-simplex on
`U` and `H = G_N - U` is the induced subgraph on `[2..N] \setminus U`.

The Hodge spectrum of a join (Goff 2009; Horak-Jost 2013) gives
contributions from each factor shifted by the other's vertex count. The
`L_1` eigenvalue `|V| = k + |V(H)|` arises with multiplicity equal to
`#(0-faces of Δ^{k-1}) + #(1-faces of Δ^{k-1}) = k + C(k, 2) = C(k+1, 2)`,
matching empirical exactly. Full join-spectrum derivation is a
single-page exercise in Friedman / Horak-Jost machinery; we report the
empirical identity here and conjecture the closed-form theorem.

**Why this is structurally novel:** identity F3 connects a **purely
spectral quantity** (the multiplicity of a specific Hodge-Laplacian
eigenvalue) to a **prime-counting fingerprint** (`π(N) - π(N/2)`) via
the explicit formula `k(k+1)/2`. CLOSED_PATHS lines 356/387/E7.12/E7.16
all addressed `L_0` only and produced `λ` formulas tied to Ramanujan
sums or character sums — never to `π(N) - π(N/2)`. This is the first
spectral identity in the project tying *higher-order* Hodge data to
Bertrand-prime counts.

### F4. **Mean shift = triple-coprime singular series (closure mode E).**

For any flag complex,
`mean(spec L_1) = trace(L_1) / |E| = 2 + 3 |T| / |E|`.

Hence `mean(L_1, K_N) - mean(L_1, ER) = 3 (T_cp - T_ER) / |E|`.

Triangle counts asymptotically:
- ER at density `p`: `T_ER ~ p^3 \cdot \binom{N}{3}`.
- Coprime: `T_cp ~ \binom{N}{3} \cdot \prod_p (1 - 3/p^2 + 2/p^3)`
  (probability of pairwise coprimality of a random triple).

So `T_cp / T_ER \to \prod_p (1 - 3/p^2 + 2/p^3) / (1 - 1/p^2)^3` as
`N \to \infty`. Numerically:
- `\prod_p (1 - 3/p^2 + 2/p^3) ≈ 0.286747`
- `(6/π²)^3 ≈ 0.224675`
- **Asymptotic ratio: 1.27628**

Empirical at our scales:
| N   | T_cp/T_ER | asy ratio |
|-----|-----------|-----------|
| 32  | 1.267     | 1.276     |
| 48  | 1.241     | 1.276     |
| 64  | 1.257     | 1.276     |
| 96  | 1.262     | 1.276     |
| 128 | 1.273     | 1.276     |

`Z[mean(L_1)]` grows from 3.0σ at N=32 to **18.3σ at N=128**. Scaling is
roughly `Z ~ √N` (since `|T| ~ N^3`, `σ(T_ER) ~ N^{3/2}`, `\Delta T ~ 0.276
\cdot N^3`).

This is **closure mode E (refinement of E2.13 / Gowers `U^k` family)**:
the L_1-mean fingerprint reduces analytically to triple-coprime singular
series, a known arithmetic identity. Adds a new *empirical category*
(Hodge `L_1`-mean) to the W-trick / singular-series saturation family,
but no new arithmetic content.

### F5. KS distance and growing rejection.

The full L_1 spectrum (not just its mean) deviates from matched ER:

| N   | KS stat | p-value      |
|-----|---------|--------------|
| 32  | 0.210   | 1.9e-11      |
| 48  | 0.205   | 3.3e-24      |
| 64  | 0.229   | 9.4e-54      |
| 96  | 0.288   | 1.3e-187     |
| 128 | 0.326   | < 1e-300     |

KS distance grows with N. The rejection is overwhelming.

### F6. β_2 ≈ 0 in both coprime and ER (Z[β_2] ≈ -0.6σ at N=32).

`β_2(K_N) = 0` deterministically by F1 (contractibility). ER controls
sometimes have small `β_2 ∈ {0, 1, 2}` from random clique gaps. Z-scores
are -0.6 to -1.7σ at N≤32 (suppressed) but the ER `β_2` rapidly drops to
0 as ER density grows, so the Z degenerates beyond N=32.

## Pre-stated falsifier outcomes

| Falsifier | Outcome |
|-----------|---------|
| (E) spectra match within KS noise; β_k within 2σ | **FALSIFIED** — KS stat 0.21–0.33, p < 1e-11 across all N |
| (I) β_1 ≠ 0 with σ-significance | **TRIVIALLY FALSIFIED** — β_1 = 0 deterministically (F1) |
| (A) uniform `λ_2^{(1)} > c > 0` gap distinct from ER | **NOT INSTANTIATED** — `λ_min(L_1)` Z-score is +5.5σ at N=48 but only +2.5σ at N=64; not a uniform-in-N gap |

The pre-stated A-grade success path is NOT instantiated. The result is
**B-grade**: a sharp empirical identity (F3) plus a known-mechanism
shift (F4) plus failure of the L_1-spectral-gap A-target.

## Edge implications (proposed)

### EDGE E7.17 (proposed) — Hodge L_1 max-multiplicity = Bertrand-prime triangular count

**Statement:** For the coprimality flag complex `K_N` on `[2, N]` with
`N ≥ 3`, the Hodge Laplacian `L_1 = B_1^T B_1 + B_2 B_2^T` has
`λ_max(L_1) = |V| = N - 1` exactly, with multiplicity

    mult(λ = |V|) = (k(k+1))/2,   where k = π(N) - π(N/2).

**Verification:** empirical at 9 values of N ∈ {8, 12, 16, 24, 32, 48,
64, 96, 128}; perfect match in every cell.

**Mechanism (proof sketch):** Bertrand's postulate plus universal-prime
join decomposition `K_N = Δ^{k-1} * F(H)`. Multiplicity follows from
Horak-Jost 2013 join-spectrum formula. Full proof: single-page exercise.

**Why an edge:** first identity in the project tying a higher-order
Hodge-spectral quantity to a prime-counting fingerprint. Distinct from
E7.12 (Cayley `(Z/nZ)^*` spectrum), E7.16 (prime-Cayley Friedman), and
all other graph-Laplacian (`L_0`-only) closures. Negative-shape: the
identity *closes* the question "does L_1 detect arithmetic structure"
with a clean closed-form answer that derives from elementary Bertrand
plus join machinery, not from new arithmetic content.

**EVS rating: M** (medium edge-value-score). The identity is sharp and
new, but its derivation is short and the underlying mechanism (Bertrand
universal vertices + join formula) is classical.

### Refinement of E2.13 (Gowers / triple-coprime singular series)

F4 adds the **L_1-mean fingerprint** to the family of triple-pairwise-
coprime singular-series detectors:

    mean(L_1, K_N) - mean(L_1, ER_p) = 3 (T_cp - T_ER) / |E|
                                    → 3 \cdot (\prod_p(1-3p^{-2}+2p^{-3}) / (6/π²)^3 - 1) \cdot p \cdot (N-2) / 6
                                    ≈ 0.276 \cdot p \cdot N as N → ∞.

This is a 39th saturation measure in the project's HL-singular-series
family.

## Grade and closure

**Self-grade: B (substantive refinement + new empirical identity).**

Reasoning:
- Cross-domain Hodge Laplacian was previously PROPOSED in
  `CROSS_DOMAIN_TECHNIQUES.md` §1; this is its first quantitative use;
  the import is real and produces work, not a name.
- F3 is a clean new identity tying spectral multiplicity to Bertrand-
  prime counts; verified at 9 N values; explained by a known
  decomposition (join formula).
- F4 is duplicate-plus of triple-coprime singular series (E2.13 / line
  387 family).
- The pre-stated A-target ("uniform `λ_2^{(1)}` gap") is NOT met;
  spectral gap Z-score is non-monotone in N.
- No polylog opening, no circuit-complexity content, no new
  pseudorandomness measure beyond singular-series family.

**Closure mode:** E (closes the line "L_1 spectrum of coprimality flag
complex" as: structurally determined by Bertrand-prime universal-
vertex count + triple-coprime singular series; no further arithmetic
content). Refines E2.13 / 387 family.

## What would falsify F3

The conjecture `mult(λ = |V|) = k(k+1)/2` would be falsified by:
- A counterexample at any specific N (none seen across N ∈ {8..128}).
- A break of Bertrand's postulate at large N (impossible: Bertrand
  proven; in fact stronger Hoheisel-type results give primes in shorter
  intervals).
- Floating-point degeneracy at large N — but the integer equality
  `λ = |V|` survived to 1e-13 absolute precision, ruling this out.

The conjecture would *strengthen* to a theorem if Horak-Jost 2013
join-spectrum machinery is applied formally; this is feasible in 1-2
sessions of careful linear algebra.

## What did NOT close

- Whether `mult(λ = |V|)` matches `C(k+1, 2)` for `N` in the millions
  (only verified to 128).
- Whether the SECOND-largest eigenvalue (or smallest nonzero) carries
  additional arithmetic content beyond join-decomposition.
- `L_2` and `L_3` spectra: tractable only at `N ≤ 24` due to triangle
  count `O(N^3)`; brief peek at N=24 shows `β_2 = 0` (consistent with
  F1) and `L_2` mean shift behavior similar to L_1.

## Successors / follow-on attack vectors

### D22.a — Higher Hodge Laplacian on the 2-skeleton at sieved N

Restrict the coprimality flag complex to a sub-complex obtained by
removing vertices below `N^{1/2}` (which removes the small-prime hubs).
Test whether `L_1` still has integer-saturating `λ_max = |V|`. **Why
interesting:** Bertrand-prime count `k = π(N) - π(N/2)` would still
control multiplicity if join structure persists. **Single-session.**

### D22.b — `β_2` and `β_3` on a SPARSER flag complex

The full coprimality flag complex is contractible. A *sparser* arithmetic
flag complex (e.g., `K^{(2)}_N := \{σ : pairwise coprime AND no element
of σ exceeds N/2\}`, removing universal Bertrand primes by hand) might
have non-trivial `β_2, β_3`. Hodge identities at higher dimensions are
the cleanest entrypoint. **Single-session at N ≤ 64.**

### D22.c — Frame-shift to sub-Laplacian fingerprints of `chi_P` directly

Skip the flag complex; compute Hodge Laplacian of the `chi_P`-INDUCED
simplex structure (only include `n` with `chi_P(n) = 1` as vertices).
**Single-session.** Tests whether the prime sub-graph of the coprimality
graph has non-trivial Hodge content.

## Cross-domain refs (cited)

- **Eckmann 1944** "Harmonische Funktionen und Randwertaufgaben in
  einem Komplex" Comment. Math. Helv. 17, 240. (Discrete Hodge theorem.)
- **Friedman 1996** "Computing Betti numbers via combinatorial
  Laplacians" Algorithmica 21, 331. (Algorithm.)
- **Horak-Jost 2013** "Spectra of combinatorial Laplace operators on
  simplicial complexes" Adv. Math. 244, 303. (Join-spectrum machinery.)
- **Lim 2020** "Hodge Laplacians on graphs" SIAM Review 62(3), 685
  = arXiv:1507.05379. (Survey; first quantitative project use.)
- **Goff 2009** "Spectra of cycle and path families of orthogonal
  subgraphs" J. Algebraic Combin. (Used for join-spectrum identities.)
- **Kahle 2009, 2014** "Topology of random clique complexes" Discrete
  Math 309 / "Sharp vanishing thresholds for cohomology of random flag
  complexes" Annals 179. (ER flag-complex baseline.)

## Edges cited

E2.13 (Gowers `U^k` triple-coprime), E2.14 (Anderson `chi_P`), E2.16
(DPP), E2.17 (PH on gaps), E2.19 (subword complexity), E7.12 (Cayley
`(Z/nZ)^*` spectrum), E7.16 (prime-Cayley Friedman). Distinct from
CLOSED 387 (`L_0` only).
