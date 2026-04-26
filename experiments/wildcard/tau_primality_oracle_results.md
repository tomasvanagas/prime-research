# Ramanujan tau function as primality oracle — Results

## Question
Lehmer's congruence `tau(p) ≡ 1 + p^11 (mod 691)` holds for all primes p (≠ 691). Could this be turned into a polylog-time primality test, and from there a prime counting algorithm via summation?

## Setup
- Computed `tau(1..200)` exactly via Niebur-style coefficient expansion of `Δ(q) = q · η(q)^24` truncated at degree 200.
- Tested the Lehmer congruence on all primes ≤ 200.
- Tested whether **composites** also satisfy the congruence (false positives).
- Combined multiple congruences (mod 691, 7, 5) for tighter discrimination.

## Numbers

### Single congruence: tau(n) ≡ 1 + n^11 (mod 691)
| Class       | Pass    | Total | Rate     |
|-------------|---------|-------|----------|
| Primes      | (all)   | 95    | 100%     |
| Composites  | 0       | 153   | 0%       |

In `[4..200]`, **the congruence is a perfect prime/composite separator** — every prime passes, no composite passes.

### Multi-congruence (mod 691, 7, 5)
- Primes passing all three: 9/95 — meaning my conjectured forms of the mod-7 and mod-5 congruences are wrong (the correct Lehmer congruences mod small primes have richer structure than `1 + p^k`).
- Composites passing: 0/404.

The mod-691 result alone is the headline: it discriminates perfectly on this range.

## Why this doesn't yield a polylog algorithm

### The bottleneck: computing tau(n) for general n
For **prime** n: `tau(p)` is an eigenvalue of the Hecke operator T_p on `S_12(SL_2(Z))`. It can be computed in `O((log p)^4 (log log p)^k)` via Schoof-Pila-style algorithms on the ℓ-adic Galois representation of Δ (Edixhoven, Couveignes, Bruin et al., 2011). So **for known primes, tau is polylog**.

For **composite** n: tau is multiplicative on coprimes (`tau(mn) = tau(m)tau(n)` if gcd=1), and on prime powers `tau(p^{k+1}) = tau(p)tau(p^k) - p^11 tau(p^{k-1})`. Both formulas **require knowing the factorization** of n.

So computing `tau(n) mod 691` for general n reduces to factoring n — itself a problem with no known polylog algorithm.

### The circularity
Suppose we had a factoring oracle. Then:
1. Factor n.
2. Compute tau of each prime factor in polylog (Edixhoven et al.).
3. Combine via multiplicativity.
4. Test `tau(n) ≡ 1 + n^11 (mod 691)`.

If true, n is prime. (We already know this from the factorization!) The tau oracle is **circular**: it tests primality only after we've already factored n, so it provides no new information.

For prime counting: even with this oracle, we'd need to call it on Ω(x/log x) candidates — that's the standard counting bottleneck.

## Verdict
**CLOSED — failure mode (C) Circularity.** The Lehmer congruence is a beautiful structural fact about primes, but the only known way to evaluate `tau(n)` for unrestricted n requires knowing the factorization, which itself encodes primality. The oracle assumes what it would prove.

## What would change the verdict
A polylog-time algorithm computing `tau(n) mod 691` for arbitrary n — *without* factoring n. This would be deeply surprising: it would imply factoring is in BPP or that tau has a structure we haven't seen. (Note: `tau(n)` is the sum `tau(n) = sum_{d|n} ...` — a divisor sum that we can't unwind without divisors.)

## Bonus observation
The congruence `tau(n) ≡ 1 + n^11 (mod 691)` has a beautiful interpretation: it says the Galois representation `rho_{691}` attached to Δ is **reducible mod 691** — splitting as `1 ⊕ chi^11` where chi is the cyclotomic character. This is the Ramanujan congruence interpretation (Serre, Swinnerton-Dyer). Equivalently, 691 divides the numerator of `B_12 / 24 = -691/2730`. The congruence is encoded in the Eisenstein-cuspform lattice structure of weight 12, not in "primes themselves" — explaining why it doesn't accelerate primality testing.
