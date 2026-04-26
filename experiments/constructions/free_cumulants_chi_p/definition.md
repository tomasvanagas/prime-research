# Construction Definition: Free cumulants of the chi_P MPS unfolding operator

**Composition target:** C2 (NOVELTY_CHALLENGES.md §1).
**Composes edges:** E2.1 (TT/MPS bond-dim identity) × free probability framework.
**Cross-domain import:** Voiculescu free probability (Mingo & Speicher,
"Free Probability and Random Matrices", Fields Inst. Monographs 35, 2017).

## Object

Let `chi_P : [1, W^d] → {0, 1}` be the prime indicator. View it as a vector in
`(C^W)^{⊗d}` via the natural base-W reshape. For a cut `1 ≤ j < d` the
unfolding `M^(j)` is a rectangular matrix of shape `(W^j, W^{d-j})` with

```
M^(j)[i, k] = chi_P(i · W^{d-j} + k + 1).
```

Define the positive self-adjoint operator
```
A_chi^pos = M^(j) M^(j)^T   ∈   M_{W^j}(R)
```
with spectral measure under the tracial state `τ(X) = (1/W^j) tr(X)`. The
empirical spectral distribution of `A_chi^pos` is the law of its eigenvalues,
which equal the squared singular values of `M^(j)`.

The construction's evaluator is the free-cumulant transform applied to the
empirical squared-singular-value distribution, standardized to mean 1.

## Free cumulants

Given empirical raw moments `(m_1, m_2, m_3, m_4)`, the first four free
cumulants `(κ_1, κ_2, κ_3, κ_4)` are the unique solution of the moment-
cumulant relation summed over non-crossing partitions:

```
κ_1 = m_1
κ_2 = m_2 − m_1^2
κ_3 = m_3 − 3 m_1 m_2 + 2 m_1^3
κ_4 = m_4 − 4 m_1 m_3 − 2 m_2^2 + 10 m_1^2 m_2 − 5 m_1^4
```

(The κ_4 formula differs from the classical cumulant κ_4^{cl} = m_4 −
4 m_1 m_3 − 3 m_2^2 + … by the −2 versus −3 coefficient on `m_2^2`. This is
where free / classical first separate.)

## Reference free distributions

- **Wigner semicircle** (squared = MP(1)): support [0, 4·σ²], free cumulants
  `κ_2 = σ², κ_n = 0 for n ≠ 2`. Standardized: `κ_r ∈ {1, 1, 0, 0}` is
  diagnostic for n ≥ 3.
- **Marchenko–Pastur MP(c)** (= free Poisson rate-1/c jump-c): the limiting
  spectral distribution of `(1/m) X X^T` for `X` an `m × n` iid matrix of
  variance 1, `c = n/m`. Mean = 1, variance = c. After standardization to
  mean 1, the empirical free cumulants are
  ```
  κ_r(MP(c)) = c^{r−1}.
  ```
  Verified against Gaussian iid simulation (m=4000, n=1000, c=0.25) to 3
  decimals for r ≤ 4.

## Relation to π(x)

The squared singular values of `M^(j)` carry strictly more information than
their multiset cardinality (= rank), which is the content of E2.1:
```
rank M^(j) = min(W^j, φ(W) · W^{d-j-1} + 1).
```
The free-cumulant transform reduces the multiset to four real numbers
that characterize the *bulk* spectral structure. If chi_P's bulk matches
MP(c) for some `c`, then chi_P's MPS spectrum is "free-Poisson-like" up to
finitely many outliers — a strictly stronger statement than rank deficit.

Conversely, if chi_P's bulk free cumulants do not match MP(c) at the
expected `c = φ(W)/W` after a finite outlier projection, the deviation is a
new pseudorandomness signature beyond the 35-measure battery in
`novel/pseudorandomness_of_pi.md`.

## Outputs

The evaluator produces, for each `(W, d)` configuration with `j = d/2`:

1. The empirical raw moments `(m_1, m_2, m_3, m_4)` of squared singular
   values, under several standardizations (full spectrum, nonzero-only,
   bulk after dropping `k` leading singular values).
2. The free cumulants `(κ_1, κ_2, κ_3, κ_4)` for each standardization.
3. The MP(c=φ(W)/W) reference cumulants for comparison.
4. The "outlier count" `k_*(W, d) = min{ k : the drop-k bulk free
   cumulants match MP(c) within relative tolerance ε }` for ε ∈ {0.1,
   0.2, 0.5}.
5. Two iid Bernoulli baselines for sanity:
   - `bernoulli` — full unfolding shape, density matched to π(N)/N.
   - `active_bernoulli` — restricted to the active block of shape
     `(W^j − 1, φ(W) · W^{d-j-1})`, density matched to coprime-residue
     prime density `π(N) · W / (N · φ(W))`. This is the matched random
     surrogate for the non-trivial part of `M^(j)`.

## Why this is original to this project

The S53 file `experiments/wildcard/free_probability_delta.py` measured
free moments of the *scalar* fluctuation `δ(x)/√x` and found classical
Gaussian (CLT, not free CLT). That probed `δ(x)` as a 1D process.

This construction probes the *operator* `A_chi^pos` derived from the
*MPS unfolding* of `chi_P` itself. The two are different objects:
the S53 measure tests "is the partial-sum fluctuation freely
distributed?" and answered NO. This construction tests "is the
matrix-spectrum of `chi_P` freely Poisson-distributed?" — a
strictly different question whose answer is partially YES (bulk
matches MP) with a quantifiable deviation (outlier count).

E2.1 already establishes that `M^(j)` has a deterministic rank deficit.
What was previously unknown: the *singular value distribution within
the rank* and how it compares to a matched random ensemble.
