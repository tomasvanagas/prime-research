# Recursive Prime Counting: FMM-Inspired Aggressive Decomposition

## Hypothesis

Push the Meissel-Lehmer recursion aggressively (FMM-inspired):
- pi(x) decomposes into pi(x^{1/k}) subproblems plus correction terms
- Recursion depth is O(log log x) (since x -> x^{1/2} -> x^{1/4} -> ...)
- If each level costs O(polylog), total is O(polylog * log log x) = O(polylog)

## Prior Closed Paths

- **#368**: Recursive identity pi(x) via pi(x/d) — correction encodes primes in interval
- **#560**: Hierarchical sieve (FMM Phi) — Phi calls/x = 0.03-0.04 constant (linear)
- **#605**: Primorial decomposition — optimal c=3 gives O(x^{2/3}), IS Meissel-Lehmer

## Experiments Run

### Part 1: Meissel-Lehmer with explicit work tracking

Tracked arithmetic operations at each level of the phi recursion tree.

| x     | Total ops | ops/x^{2/3} | ops/log(x)^3 | Max depth |
|-------|-----------|-------------|---------------|-----------|
| 10^3  | 82        | 0.82        | 0.2           | 4         |
| 10^4  | 552       | 1.19        | 0.7           | 8         |
| 10^5  | 3,103     | 1.44        | 2.0           | 15        |
| 10^6  | 13,819    | 1.38        | 5.2           | 25        |
| 10^7  | 74,067    | 1.60        | 17.7          | 47        |

**Key finding**: ops/x^{2/3} converges to constant (~1.4). ops/log(x)^3 is growing (0.2 -> 17.7). Work is **Theta(x^{2/3})**, not polylog.

The phi recursion depth grows as pi(x^{1/3}) = O(x^{1/3}/ln x), which is exponential in log x. Work per level grows toward the leaves, with the deepest levels accumulating the most operations.

### Part 2: Aggressive recursive decomposition

Tried maximizing recursive decomposition: pi(x) -> pi(x^{1/3}) -> pi(x^{1/9}) -> ...

| x     | Total work | work/x^{2/3} | Recursion depth |
|-------|-----------|-------------|-----------------|
| 10^3  | 93        | 0.93        | 1               |
| 10^4  | 582       | 1.25        | 1               |
| 10^5  | 3,064     | 1.42        | 1               |
| 10^6  | 13,198    | 1.32        | 1               |
| 10^7  | 67,460    | 1.45        | 2               |

**Key finding**: 99.9%+ of all work is at depth 0. The top-level phi computation dominates entirely. Deeper recursion reduces sub-problem sizes to trivial (x^{1/3} is small quickly) but does NOT reduce the top-level work.

Depth 0 alone: 67,415 ops out of 67,460 total for x=10^7.

### Part 3: Hierarchical sieve compression

Separated sieve into small primes (p <= x^eps) and large primes.

At x = 10^7, eps = 1/3 (threshold = 216):
- Small primes: 47, Large: 400
- Total large prime removals: 3,957,923
- Average removals per large prime: 9,895

The large primes collectively make O(x^{2/3}) removals at irregular positions. These positions depend on the primes themselves — no periodic or compressible structure.

### Part 4: Theoretical recursion depth analysis

At x = 10^100 (the target), the recursion chain has 7 levels:

| Depth | log(y) | phi_work at this level |
|-------|--------|----------------------|
| 0     | 230.3  | 2.0 * 10^64          |
| 1     | 115.1  | 1.9 * 10^31          |
| 2     | 57.6   | 8.1 * 10^14          |
| 3     | 28.8   | 7.5 * 10^6           |
| 4     | 14.4   | 1,020                |
| 5     | 7.2    | 17                   |
| 6     | 3.6    | 3                    |

**Depth 0 dominates by a factor of 10^{33} over depth 1.** The recursion depth IS O(log log x) = 7, but the work at the top level is O(x^{2/3}) = O(10^{66.7}), utterly non-polylog.

## Why the FMM Analogy Fails

The Fast Multipole Method works because:
1. The 1/r kernel is **smooth** away from the singularity
2. Smooth functions compress hierarchically (multipole expansions)
3. Each level of the hierarchy does O(N) work with bounded interaction lists

For prime counting:
1. The prime indicator function has **no smooth kernel**
2. The phi recursion tree branches on floor(x/product-of-primes) values
3. These floor values are **not compressible** — they encode the very primes we seek
4. The tree has O(x^{2/3}/ln x) leaves, each requiring O(1) work = total O(x^{2/3})

The "correction terms" at each recursion level encode the irregularity of primes. This irregularity is exactly what makes the problem hard.

## Why Hierarchical Sieve Compression Fails

- Small primes (p <= x^eps): periodic sieving pattern, efficiently computable
- Large primes (p > x^eps): each removes few elements, but collectively O(x^{2/3}) removals
- The positions of these removals are **at multiples of primes we don't yet know**
- No way to batch or compress these without knowing the primes themselves (circularity)

## Verdict: FAIL

**This reduces to Meissel-Lehmer (confirms paths #368, #560, #605).**

The aggressive recursion achieves O(log log x) depth as hoped, but the work is catastrophically concentrated at the top level: depth 0 accounts for >99.9% of all operations, and that work is Theta(x^{2/3}/ln x).

The fundamental issue: the phi(x, a) computation at the top level IS the Meissel-Lehmer bottleneck. No amount of recursive decomposition of sub-problems helps because the sub-problems are already tiny — the hard part is the combinatorial explosion at the top.

For x = 10^100: total work ~ 10^{66.7}, versus the polylog target of ~10^3. The gap is 63+ orders of magnitude.

**Close as: equivalent to Meissel-Lehmer, O(x^{2/3}/ln x) confirmed.**
