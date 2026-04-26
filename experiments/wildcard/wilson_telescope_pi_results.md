# Wilson Telescoping for pi(x) — Results

## Question
Wilson's theorem gives a primality indicator: `(n-1)! ≡ -1 (mod n)` iff n is prime. Can we compute the *sequence* `{(n-1)! mod n}_{n=2..x}` collectively, using product-tree or telescoping structure, faster than `O(x · M(log x))`?

## Setup
- Verified `pi(200) = 46` via Wilson formula (correct).
- Compared two methods for computing `(n-1)! mod n` per n in `[2..N]`:
  - **Naive**: Sequential modular multiplications, O(n) per n → O(N²) total.
  - **Balanced product tree**: Recursive halving, O((log n)² · M(log n)) per n.

## Numbers
| N    | Naive (s) | Tree (s) | Both correct? |
|------|-----------|----------|---------------|
| 200  | 0.001     | 0.002    | yes (46)      |
| 1000 | 0.019     | 0.074    | yes (168)     |
| 3000 | 0.168     | 0.723    | yes (430)     |

In Python's bignum, the tree approach is *slower* than naive (overhead of recursion + Python integer multiplication is faster naive than tree for these sizes). Asymptotically both are O(N polylog N) — worse than the basic sieve `O(N log log N)`.

## Key observation: barrier structure
The fundamental obstruction is **incompatible moduli**:
- `(n-1)! mod n` is computed in ring `Z/nZ`. For different n, these rings are unrelated.
- Sharing work would require a single ring containing all `Z/nZ` for `n ≤ N`. The only such ring is `Z/L` where `L = lcm(1..N)`. By PNT, `log L = ψ(N) ~ N` — exponential size.
- Even computing `N! mod L` once would need O(N log L · M(log L)) = O(N² polylog) bit operations.

So Wilson's theorem provides a clean *characterization* of primes, but its computational complexity is dominated by:
1. Either O(N) per modulus (no sharing), totalling O(N² polylog), or
2. A single huge modulus L = lcm(1..N) of size 2^N.

Either way, **information-theoretically blocked from polylog**.

## Verdict
**CLOSED — failure mode (E) Equivalence.** Wilson's theorem reduces prime counting to factorial-modular computation, but each `(n-1)! mod n` requires Ω(log n) bit work and has no algebraic structure shared across different n. The sequence has the same information content as π(x) itself. No telescoping or product-tree shortcut bridges the gap.

## What would change the verdict
A *single* polylog-computable integer Q (not a sequence) such that `pi(N) = f(Q)` with f also polylog. Wilson's theorem doesn't provide such a Q — it provides a *family* of factorial residues, one per modulus.
