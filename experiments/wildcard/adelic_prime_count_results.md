# Adelic Local-Global Prime Counting — Results

**Experiment:** `experiments/wildcard/adelic_prime_count.py`
**Date:** 2026-04-05
**Verdict:** CLOSED — No speedup from adelic/CRT decomposition of pi(x)

## Idea

Use the adelic (local-global) principle: compute pi(x) mod q for several small primes q, then reconstruct pi(x) via CRT. The key question is whether computing pi(x) mod q is cheaper than computing pi(x) itself.

## Prior Closed Paths

This is the 6th CRT-based attempt:
- S24: CRT Prime Locator — pi(x;q,a) as hard as pi(x)
- S24: Adelic local-global reconstruction — L-function zeros needed per modulus
- S29: CRT from modular periods — p(n) mod m NOT periodic for m >= 5
- S33: p-adic lifting via CRT — floor division not a ring homomorphism

## Results

### Test 1: Liouville/Mertens Parity vs pi(x) mod 2

L(x) mod 2 and M(x) mod 2 match pi(x) mod 2 approximately 57% of the time (7 test points), consistent with random (50% baseline). No predictive relationship.

### Test 2: Smooth Approximation R(x) mod q

| x | pi(x) | round(R(x)) | error | mods correct (of 6) |
|---|-------|-------------|-------|---------------------|
| 100 | 25 | 26 | +1 | 0/6 |
| 500 | 95 | 94 | -1 | 0/6 |
| 1,000 | 168 | 168 | 0 | 6/6 |
| 5,000 | 669 | 669 | 0 | 6/6 |
| 10,000 | 1,229 | 1,227 | -2 | 1/6 |
| 50,000 | 5,133 | 5,133 | 0 | 6/6 |
| 100,000 | 9,592 | 9,587 | -5 | 1/6 |
| 500,000 | 41,538 | 41,529 | -9 | 1/6 |
| 1,000,000 | 78,498 | 78,527 | +29 | 0/6 |

R(x) gives correct mod-q results only when the absolute error happens to be < 1 (i.e., when R(x) rounds to the exact answer). The error grows as O(sqrt(x)), so correctness **degrades** with x. Taking mod q provides zero benefit over rounding.

### Test 3: Timing — pi(x) vs AP Decomposition

Computing pi(x) mod q via arithmetic progression decomposition is **100-1000x SLOWER** than computing pi(x) directly via sympy's sieve. The AP approach enumerates primes per residue class, which is strictly more work.

### Test 4: Pattern Analysis of pi(x) mod q

- **Autocorrelation:** ~0.877 for all q (= 1 - prime_density near x=10000). pi(x) mod q changes only at primes.
- **Residue frequencies:** Approximately uniform (chi-squared reasonable for q=2,3; slightly elevated for q=5,11,13 but not extreme).
- **No exploitable pattern:** The sequence pi(x) mod q changes at unpredictable positions (primes), and the jump direction depends on which residue class the prime falls in.

### Test 5: CRT Reconstruction

| x | pi(x) | Min moduli needed | Product |
|---|-------|-------------------|---------|
| 100 | 25 | 3: {2,3,5} | 30 |
| 1,000 | 168 | 4: {2,3,5,7} | 210 |
| 10,000 | 1,229 | 5: {2,3,5,7,11} | 2,310 |
| 100,000 | 9,592 | 6: {2,3,5,7,11,13} | 30,030 |

CRT reconstruction works correctly (verified). Only O(log log x) moduli are needed since primorial grows super-exponentially. But the bottleneck is computing **each** pi(x) mod q, which costs as much as pi(x) itself.

### Test 6: Error Growth

| x | |R(x) - pi(x)| |
|---|----------------|
| 1,000 | 0.34 |
| 10,000 | 2.09 |
| 100,000 | 4.58 |
| 1,000,000 | 29.39 |

Error grows roughly as sqrt(x)/log(x), consistent with known estimates. For pi(x) mod q to work from R(x), would need error < 0.5, which requires the oscillatory correction from zeta zeros — the same bottleneck as exact pi(x).

## Why This Fails

1. **Computing pi(x) mod q is NOT easier than pi(x):** No known method computes any non-trivial function of pi(x) without essentially computing pi(x).

2. **Floor division breaks modularity:** The key sieve identity pi(x) = ... involves floor(x/d) terms. floor(x/d) mod q != floor((x mod q)/(d mod q)) because floor division is not a ring homomorphism. This kills all attempts to "work mod q."

3. **Error analysis confirms theory:** Truncating the explicit formula at N zeros gives error ~ Cx/(N log x). For exact pi(x), need error < 0.5. For pi(x) mod q, still need error < 0.5 (because we round first, then take mod). No improvement.

4. **Information-theoretic barrier:** pi(x) mod q for all q up to some bound encodes pi(x) itself via CRT. So computing all the mod-q values is at least as hard as pi(x).

## Verdict

**CLOSED.** The adelic/CRT approach provides no speedup. This is the 6th independent confirmation that CRT-based decomposition of prime counting fails. The fundamental reason is that the arithmetic of prime counting (sieving, floor division) does not decompose modularly.
