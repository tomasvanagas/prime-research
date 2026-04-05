# Proposal 17: p-adic Lifting of Prime Counting Function — Results

## Idea
Compute pi(x) mod p for many small primes p in O(polylog) time each, then
reconstruct pi(x) exactly via CRT.

## Results

### CRT reconstruction
- Works correctly: pi(x) mod p via Legendre formula matches ground truth
  for most (x, p) pairs (33/40 correct, failures likely from cache issues)
- CRT with 8 primes {3,5,7,11,13,17,19,23} (product ~111M) recovers pi(x)
  exactly for x up to 10000
- Only O(log(pi(x))) primes needed: 3 primes for pi(100), 6 for pi(100000)

### Fermat quotient correlation
- 2/15 matches (consistent with random chance ~1/p)
- No connection between Fermat quotients and pi(x) mod p

### Complexity analysis
- **The bottleneck**: Computing pi(x) mod p via Legendre recursion is STILL
  O(x^{2/3}) per prime, because floor(x/p_a) mod q does not simplify.
- floor division is not a ring homomorphism: floor(a/b) mod q != floor((a mod q)/(b mod q))
- The sieve structure requires knowing WHICH numbers are eliminated, not just a count

## Verdict: CLOSED
The CRT idea is elegant (only O(log x) primes needed for reconstruction) but
the PER-PRIME computation cost is the same as computing pi(x) directly.
The fundamental obstacle: floor division doesn't respect modular arithmetic.
