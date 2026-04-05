# Proposal Session: Fresh Approaches to O(polylog) p(n)
## Session 33 — April 5, 2026

**Methodology:** 4 proposals developed, each with runnable Python code + computational tests on n < 10000.
**Internet search:** Searched for 2025-2026 papers on dequantization, Ramanujan Library (ICLR 2025), p-adic L-functions, spectral methods, compressed sensing, and trace formula shortcuts.

---

## Proposal 16: Spectral Truncation with Adaptive Zero Selection

**Idea**: Select O(polylog(n)) "resonant" zeta zeros (those whose
gamma * log(x) is near k*pi) instead of all O(sqrt(x)) zeros.

**Key result**: For n < 100, just 1-2 zeros suffice. For n >= 500,
even 1000 zeros are insufficient.

**Why it fails**: Which zeros are "resonant" depends on x in an unpredictable
way. The number of significant zeros grows as O(x^{1/2}) by the explicit
formula's error term. Adaptive selection provides ~40% improvement over
sequential but cannot reduce the O(x^{1/2}) requirement to O(polylog).

**Conjecture status**: REFUTED by computation.

---

## Proposal 17: p-adic Lifting via CRT

**Idea**: Compute pi(x) mod p for O(log x) small primes, reconstruct via CRT.

**Key result**: CRT works perfectly — only 3-6 primes needed for reconstruction.
But computing pi(x) mod p costs O(x^{2/3}) per prime.

**Why it fails**: floor(x/d) mod q != floor((x mod q)/(d mod q)).
Floor division is not a ring homomorphism, so modular reduction
doesn't simplify the Legendre recursion.

**Novel insight**: The information-theoretic cost is low (O(log x) bits via
O(log x) primes), but the COMPUTATIONAL cost per bit is O(x^{2/3}).
This cleanly separates information complexity from computational complexity.

---

## Proposal 18: Dequantized Grover Counting

**Idea**: Apply classical dequantization (Tang 2024, Chia et al. 2025) to
the Mobius-inversion formulation of pi(x).

**Key result**: Truncated Mobius sums diverge badly. Sampling-based estimates
have O(x^2/S) variance. Hyperbola method gives O(sqrt(x)) but is already known.

**Why it fails**: Dequantization requires low-rank structure in the input.
The Mobius function mu(d) is pseudorandom with no exploitable structure.
Prime indicator has approximate degree N/2 — maximum possible — confirming
full-rank resistance to polynomial approximation.

---

## Proposal 19: Ramanujan Library / PSLQ on delta(n)

**Idea**: Apply automated conjecture generation (PSLQ/LLL, inspired by the
Ramanujan Library ICLR 2025) to discover structure in delta(n) = p(n) - R^{-1}(n).

**Key result**: delta(n) has range [-130, 102], std ~33, strong short-range
autocorrelation (0.95 at lag 1). BUT:
- No linear recurrence (orders 1-6, all errors ~4)
- No modular patterns (chi-squared all non-significant)
- PSLQ finds different integer relations for each n (no universal formula)
- Best basis decomposition RMSE = 31 (too large for rounding to exact)

**Why it fails**: delta(n) is determined by the phases of O(sqrt(p(n))) zeta
zeros, which are GUE-distributed. No finite basis of known functions captures
this pseudorandom structure.

---

## Proposal 20: Étale Cohomology Point-Counting

**Idea**: Encode pi(x) as point-counting on an algebraic variety, exploiting
the Grothendieck-Lefschetz trace formula which runs in O(polylog) time for
fixed varieties.

**Key result**: Character sums correctly compute pi(x, q, a), and Hasse-Weil
constraints provide enough information bits. BUT computing each constraint
costs O(x^{2/3}).

**Why it fails**: Point-counting on a FIXED variety is O(polylog), but
the primality predicate over [1,x] has UNBOUNDED algebraic complexity.
No fixed-dimensional variety encodes pi(x) for all x.

---

## Cross-cutting Insights

All four proposals encounter the SAME fundamental barrier from different angles:

1. **Spectral (P16)**: Need O(sqrt(x)) zeros → too many
2. **Modular (P17)**: Floor division breaks homomorphism → can't lift mod p
3. **Sampling (P18)**: mu(d) is full-rank → no low-rank shortcut
4. **Algebraic (P20)**: Primality has unbounded algebraic complexity → no fixed variety

These are all manifestations of the **information-computation gap**: the BITS of
information in delta(n) are O(log n), but they are ENCODED in the interaction of
O(sqrt(x)) independent pieces (zeta zeros, sieve values, Mobius function values).

No proposal found a way to extract O(log n) bits without processing
O(x^{alpha}) pieces for some alpha > 0.

## Files
- `experiments/proposals/proposal16_spectral_truncation_bound.py` + `_results.md`
- `experiments/proposals/proposal17_padic_lifting.py` + `_results.md`
- `experiments/proposals/proposal18_dequantized_grover_count.py` + `_results.md`
- `experiments/proposals/proposal19_ramanujan_library_delta.py` + `_results.md`
- `experiments/proposals/proposal20_etale_cohomology_count.py` + `_results.md`
