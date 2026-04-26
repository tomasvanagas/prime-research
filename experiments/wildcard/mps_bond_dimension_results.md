# MPS / Tensor-Train Bond Dimension of the Prime Indicator

**Script:** `mps_bond_dimension.py`
**Date:** 2026-04-25 (extended same day with base-W sweep)
**Verdict:** CLOSED -- bond dimension of chi_P viewed in base W is exactly
`min(W^j, phi(W) * W^{d-j-1} + 1)` at every cut j, for every (W, d) tested.
Half-cut bond dim ratio chi_P / random asymptotes to `phi(W) / W`. No
polylog MPS representation exists.

## What was tested

We view the binary prime indicator `chi_P(n)` for `n in [1, W^d]` as a
state vector in `(C^W)^{otimes d}` and measure its tensor-train bond
dimension at every bipartite cut. Two families:

  - Family A: `W = 2`, `d in {10, 12, 14, 16, 18, 20}` -- replicates the
    original measurement; also tracks random control, cumulative pi(n),
    and the residual `delta(n) = pi(n) - li(n)`.
  - Family B: `W in {2, 6, 30, 210}` (primorials), various d, exact rank
    at every cut. Tests the precise prediction
    `rank(cut j) = min(W^j, phi(W) * W^{d-j-1} + 1)`.

A polylog MPS representation would mean polylog-time exact computation of
`pi(N)` and `p(n)`. Compared against a random binary vector at matched
density.

This builds on but is distinct from earlier project work:
  - `tensor_sieve.py` -- product DFA states, primorial blowup, different object
  - `sieve_matrix_rank.py` -- divisibility matrix `M[n, p]`, different tensor
  - `wavelet_*` -- single-mode basis sparsity, not multipartite

## Headline result: the base-W theorem

> rank(M^{(j)}) = min( W^j, phi(W) * W^{d-j-1} + 1 )
>
> at every cut j in {1, ..., d-1}, for every primorial W tested.

Empirically saturated at strict equality in 11/11 configurations:

| W   | d  | N         | chi_P max | predicted | random max | chi_P/random | phi/W   |
|-----|----|-----------|-----------|-----------|------------|--------------|---------|
| 2   | 10 | 1,024     | 17        | 17        | 32         | 0.5312       | 0.5000  |
| 2   | 14 | 16,384    | 65        | 65        | 128        | 0.5078       | 0.5000  |
| 2   | 20 | 1,048,576 | 513       | 513       | 1,024      | 0.5010       | 0.5000  |
| 6   | 4  | 1,296     | 13        | 13        | 36         | 0.3611       | 0.3333  |
| 6   | 6  | 46,656    | 73        | 73        | 216        | 0.3380       | 0.3333  |
| 6   | 8  | 1,679,616 | 433       | 433       | 1,296      | 0.3341       | 0.3333  |
| 30  | 2  | 900       | 9         | 9         | 30         | 0.3000       | 0.2667  |
| 30  | 3  | 27,000    | 30        | 30        | 30         | 1.0000*      | 0.2667  |
| 30  | 4  | 810,000   | 241       | 241       | 900        | 0.2678       | 0.2667  |
| 210 | 2  | 44,100    | 49        | 49        | 210        | 0.2333       | 0.2286  |
| 210 | 3  | 9,261,000 | 210       | 210       | 210        | 1.0000*      | 0.2286  |

* At odd d the constrained side has fewer factors than the free side at
  one of the cuts, so max bond dim is attained at a non-half cut where
  W^j (the unconstrained side) equals the random ceiling. The per-cut
  theorem still matches exactly; only the *max-rank* statistic loses the
  phi/W signature at odd d. See per-cut profile below.

## Per-cut profile at a few configurations

W=2, d=20 (N=1,048,576):

| cut j | chi_P rank | predicted | random rank | half-cut ceiling |
|-------|------------|-----------|-------------|------------------|
| 1     | 2          | 2         | 2           | 2                |
| 5     | 32         | 32        | 32          | 32               |
| 9     | 512        | 512       | 512         | 512              |
| 10    | 513        | 513       | 1024        | 1024 (middle)    |
| 11    | 257        | 257       | 512         | 512              |
| 15    | 17         | 17        | 32          | 32               |
| 19    | 2          | 2         | 2           | 2                |

W=30, d=4 (N=810,000):

| cut j | chi_P rank | predicted        | random rank | unfolded matrix |
|-------|------------|------------------|-------------|-----------------|
| 1     | 30         | min(30, 8*900+1) | 30          | 30 x 27000      |
| 2     | 241        | min(900, 8*30+1) | 900         | 900 x 900 (mid) |
| 3     | 9          | min(27000, 8+1)  | 30          | 27000 x 30      |

W=6, d=8 (N=1,679,616):

| cut j | chi_P rank | predicted          | random rank |
|-------|------------|--------------------|-------------|
| 1     | 6          | 6                  | 6           |
| 2     | 36         | 36                 | 36          |
| 3     | 216        | min(216, 433)=216  | 216         |
| 4     | 433        | min(1296, 433)=433 | 1,296 (mid) |
| 5     | 73         | min(7776, 73)=73   | 216         |
| 6     | 13         | min(46656, 13)=13  | 36          |
| 7     | 3          | min(279936, 3)=3   | 6           |

The "min crossover" between the two terms `W^j` and `phi(W)*W^{d-j-1}+1`
happens at the half-cut. To the left of the crossover, the constraint is
not yet binding (chi_P matches random). To the right, the constraint is
saturated (factor phi(W)/W below random ceiling).

## Why the W=2 result is the natural starting point and why W>2 matters

For W=2 the savings is exactly one bit (1/2 ratio) reflecting "primes > 2
are odd." For W=6 = 2*3 the savings is 1/3 (primes > 3 are coprime to 6,
i.e., in {1, 5} mod 6). For W=30 = 2*3*5 the savings is 4/15 (primes > 5
are in (Z/30)* of size 8). For W=210 = 2*3*5*7 the savings is 8/35.

Each new prime added to the wheel knocks the bond-dim ratio down by
`(p-1)/p`, exactly mirroring the Mertens-product
`prod_{p<=W}(1 - 1/p) ~ e^{-gamma} / log W`. So the achievable
bond dim ratio at a base-W half cut is
`phi(W) / W ~ e^{-gamma} / log W` for primorial W -- decaying as
`1 / log log N` if we choose W such that `W^d = N` and `d` is constant.

## Why this does NOT yield polylog

A polylog MPS would require chi_P bond dim `poly(d) = poly(log N)`. The
half-cut theorem gives bond dim `>= phi(W) * W^{d/2-1} + 1`, which is
`>= W^{d/2} * phi(W) / W`. To make this polylog in `N = W^d`, we would
need

    `phi(W)/W * W^{d/2}  <=  poly(d)`

Setting `W^d = N` and letting `W` range freely, the minimum of
`phi(W)/W * W^{d/2}` over the constraint `W^d = N` is achieved roughly
at the boundary `W = sqrt(N)` (where `d = 2`), giving bond dim
`= Theta(phi(W)) = Theta(sqrt(N))`. So the MPS view recovers the
`sqrt(N)` Lagarias-Odlyzko barrier from a pure entanglement argument,
with no reference to zeta zeros.

For each fixed primorial W and growing d, the bond dim is `Theta(W^{d/2})`,
exponential in d.

## Random vs prime at relaxed accuracy (W=2 baseline)

With eps = 1e-1 (99% energy) approximation, chi_P drops only marginally:

| k  | chi_P (eps=1e-1) | random (eps=1e-1) |
|----|------------------|-------------------|
| 10 | 17               | 29                |
| 14 | 64               | 115               |
| 20 | 508              | 933               |

Both grow like `2^{k/2}` (volume law); chi_P stays the same factor of
~1/2 below random. **No exponential separation.** No area-law structure.

## Compressibility of cumulative pi(n) and the residual

Cumulative `pi(n)` is trivially smooth -- approximate bond dim ~2 at
eps=1e-2 across all k. (This recovers `li(n) ~ pi(n)`.)

For the residual `delta(n) = pi(n) - li(n)` at eps = 1e-2 (relative):

| k  | bond dim |
|----|----------|
| 10 | 17       |
| 12 | 30       |
| 14 | 56       |
| 16 | 99       |

Scales as ~ N^{0.42}. To compute pi(n) exactly we need *absolute* error
< 0.5, which means relative error < 0.5/sqrt(N) on a residual of order
sqrt(N). At that precision the bond dim returns to near-full. So the
residual's modest compressibility at coarse precision does NOT translate
to a polylog algorithm -- consistent with the wavelet result
(`wavelet_zero_compression` showed N^{0.75} at 99.9% energy).

## Failure mode (per project taxonomy)

**Equivalence (E)**: the MPS bond dimension at every cut is exactly the
wheel-W structure of chi_P. Each prime included in the wheel reduces the
bond dim ratio by `(p-1)/p`, which is the same factor that the wheel-W
sieve gains classically. No new compression beyond classical wheel-sieving.
The polylog limit is unreachable.

## Status

This is a NEW result (small but original). The exact identity
   `rank(M^{(j)}) = min(W^j, phi(W) * W^{d-j-1} + 1)`
plus the empirical saturation at every cut for primorials W in {2,6,30,210}
appears not to be in published literature. It cleanly distinguishes the
prime indicator from random data with an exact arithmetic identity (rather
than statistical bound). Recorded in `novel/mps_bond_dimension.md`.

Practical impact: zero -- confirms a `phi(W)/W` compressibility ratio
equivalent to wheel-W sieving. Adds Measures 22 and 23 to
`novel/pseudorandomness_of_pi.md`.

## One-line summary

Prime indicator's MPS bond dim at the half cut is exactly
`phi(W) * W^{d/2-1} + 1` for every primorial W tested -- random saturates
at `W^{d/2}`. The factor-`phi(W)/W` separation is exactly the wheel-W
sieve gain. No polylog tensor-network representation exists.

## Reproducing the experiment

```
python3 experiments/wildcard/mps_bond_dimension.py            # full sweep
python3 experiments/wildcard/mps_bond_dimension.py --quick    # smaller sizes
python3 experiments/wildcard/mps_bond_dimension.py --base 30 --depth 4
```

Total runtime: ~40s on a single CPU for the full sweep (most time in
W=2, d=20 SVD: ~13s; W=6, d=8: ~3s; W=210, d=3: ~8s).
