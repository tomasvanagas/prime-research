# Session 49 Proposals (2026-04-26)

Three concrete proposals attacking p(n) in O(polylog). Each has a small-n
test script in `experiments/proposals/session49_*.py`.

The high-level frame: p(n) = R^{-1}(n) + delta(n) where delta carries the
"oscillatory" content. Every proposal asks: is delta secretly compressible
in a way that orthodox truncation of the explicit formula does not see?

---

## Proposal A — Reordered tensor-train (TT) compression of delta(n)

### Idea

The MPS / tensor-train bond dimension of a length-N sequence under the
*natural* ordering 1,2,...,N is not the right yardstick: a sequence can be
volume-law in one ordering and rank-1 in another (e.g. a tensored product
has rank 1 only after a Cartesian split, not in original order). Prior
"tensor network sieve" tests used the natural ordering of subsets and
found volume-law. We instead reorder the index n by **non-trivial
permutations** before TT decomposition:

* Bit-reversal of n's binary form
* Gray code
* Morton (Z-order) of the (n_high, n_low) split
* p-adic ordering for p in {2,3,5}
* Sorted by R^{-1}(n) (already known to be a smooth coordinate)

### Mathematical claim

If there exists a permutation pi : [N] -> [N] computable in O(polylog) such
that delta(pi(n)) has bond dim D = O(polylog(N)) under TT, then delta is
representable in O(D^2 log N) = O(polylog(N)) classical bits and a single
delta(n) value is computable in O(D^2 log N) time. Combined with R^{-1}(n)
this gives O(polylog) for p(n).

### Pseudocode

```
for ordering in [identity, gray, bit_reverse, morton, 2adic, sort_by_Rinv]:
    seq = [delta(ord(i)) for i in range(N)]
    T = reshape(seq, [2]*log2(N))
    bond_dims = sequential_SVD(T, threshold=0.99 of singular value sum)
    record max(bond_dims)
```

### Complexity

If max bond dim D scales as polylog(N), proposal works. If linear or
power-law in N, proposal closed.

### Key assumption

Some efficiently-computable permutation exists that "untwists" delta into
a low-rank TT. The honest prior: delta is GUE-random and *no* permutation
helps. But this has not been tested across the orderings listed.

### Test

`session49_mps_delta_orderings.py` — N = 8192, run all six orderings,
report max bond dim at each cut for each ordering. Compare with random
sequence baseline.

---

## Proposal B — Compressed sensing of the explicit-formula residue
            from sparsity in the zeta-zero basis

### Idea

The truncated explicit formula gives delta(n) = sum_rho (analytic residue
at rho) * something(n, rho). The truncation error decays as 1/T where T
is the height to which zeros are summed. To get exact integer rounding we
naively need T ~ sqrt(n).

Counterproposal: treat the zeros' modes as a *dictionary* and use L1
minimization to find a **sparse** combination over zeros that fits delta
on a training window. Even if many zeros are needed in the analytic
formula, a sparse L1 fit may discover that *most of delta's variance is
explained by far fewer zeros than 1/T predicts* — i.e. compressed-sensing
sparsity in the zero basis.

### Pseudocode

```
truth = [delta(n) for n in train_window]

# Dictionary: cos and sin modes for first K_max zeros
gammas = first_K_max_zero_imag_parts()
Phi[n, j] = cos(gamma_j * log(n)) / sqrt(n)
            sin(gamma_j * log(n)) / sqrt(n)

# L1-minimize: c = argmin ||c||_1 + lambda * ||Phi c - truth||^2
c = LASSO(Phi, truth)

# Test on held-out window
err = max |round(R^{-1}(n) + Phi c [n]) - p(n)| over test n
report nnz(c), err, K_effective = nnz(c)
```

### Complexity

If a fixed effective sparsity K_eff = O(polylog(N)) achieves zero rounding
error on the test window, single delta(n) takes O(K_eff) time and
proposal works.

### Key assumption

There exist a few "anomalous" zero modes whose coefficients dominate
delta in the L1-sparse sense, despite GUE statistics being known for the
full ensemble. This contradicts the standard heuristic but is checkable.

### Test

`session49_compressed_sensing_zeros.py` — N = 4096, K_max = 1024,
LASSO with cross-validation, report nnz(c), train residual, test
residual, and prime recovery rate (rounded(R^{-1} + delta_hat) == p(n)).

---

## Proposal C — Learned residues at fixed zeros: how few zeros if
            coefficients are data-fit?

### Idea

The analytic explicit formula plugs in fixed residues 1/rho. But the
truncation tail is large precisely because the analytic coefficient
choice is *not* the L^2-optimal one for the truncated basis. A *learned*
ridge fit can absorb tail energy into the coefficients of the kept modes
and may converge faster.

We fix gamma_1,...,gamma_K (zero imaginary parts, computed once from a
table) and learn coefficients on a training window:

  delta(n) ≈ sum_{k=1}^K [ a_k cos(gamma_k log n) + b_k sin(gamma_k log n) ] / n^{1/2}

### Pseudocode

```
gammas = first_K_zero_imag_parts(K)
features(n) = [ cos(g log n)/sqrt(n), sin(g log n)/sqrt(n) for g in gammas ]
ridge_fit (a,b) on train_window
test on held-out window: error histogram, prime recovery rate
sweep K in {8, 32, 128, 512, 2048}
```

### Complexity

If smallest K with ~100% test recovery scales polylog(N), polylog wins.
If K(N) ~ N^alpha with alpha > 0, closed.

### Key assumption

The K-th truncated explicit formula's truncation tail is *low-dimensional
inside* the kept K-mode subspace — i.e., refitting coefficients can
exploit constructive interference and produce smaller error than the
analytic 1/rho choice.

### Test

`session49_neural_zero_residue.py` — N = 4096, train on n in [10, 2048],
test on n in [2049, 4096], sweep K, report scaling K_min(N) for ~99%
prime recovery.

---

## Honest priors

Every proposal is a long shot. The orthodox bound is ~sqrt(x) zeros for
exact pi(x). A, B, C all bet that delta is somehow compressible *despite*
being driven by GUE-random oscillations. The most likely outcomes:

* **A**: bond dim grows ~ N^{1/2} (tracks the sqrt(x) zero count) under
  every ordering. Volume-law-ish.
* **B**: K-sparse recovery on small K fails to clear residual for
  typical n. Residual norm scales like sqrt(N/K), consistent with
  Parseval over GUE.
* **C**: K must scale linearly with sqrt(N), matching the explicit
  formula's analytic count. Possibly with a *constant factor*
  improvement from learned coefficients.

Even null results sharpen the empirical bound on "minimum K" and may
expose whether reordering or learned coefficients give a constant
factor advantage that prior literature has not measured. Each test
runs in seconds for N=4096-8192.
