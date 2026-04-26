# Proposal B — Continued Fraction of Prime Constant (results)

**Run date:** 2026-04-26
**Script:** `proposal_b_cf_prime_constant.py`

## What it tests

The prime constant `α = Σ_{p prime} 2^{-p}` encodes π(·) in its binary
digits. We compute its regular continued fraction
`α = [a_0; a_1, a_2, …]` to depth 1500 (using 12000 binary digits, i.e.
all primes up to 12000 at mpmath precision 16800 bits) and probe `{a_n}`
for any deviation from Khintchine-typical statistics.

Tests applied:
- Khintchine geometric-mean K₀ (theoretical 2.6854)
- Lévy constant log q_n / n (theoretical 1.1866)
- Gauss-Kuzmin distribution fit for a_n ∈ {1,…,6}
- Autocorrelation at lags 1..10
- Power-spectrum slope (1/f^β fit on log a_n)
- k-automaticity inter/intra-variance proxy for k ∈ {2,3,5,7,10}
- Shannon entropy of bucketed a_n distribution

## Result

| Test | Empirical | Theoretical | Δ |
|---|---|---|---|
| Khintchine geomean K₀ | 2.7336 | 2.6854 | +0.05 |
| Lévy constant | 1.2056 | 1.1866 | +0.02 |
| GK match (a_n=1) | 0.4327 | 0.4150 | +0.018 |
| GK match (a_n=2) | 0.1607 | 0.1699 | −0.009 |
| GK match (a_n=3) | 0.0787 | 0.0931 | −0.014 |
| max |autocorrelation| (lag 1..10) | 0.013 | 0 | ~0 |
| power-spectrum β | −0.091 | 0 (white) | ~white |
| k-automaticity ratio (max k=10) | 0.011 | 0 | ~0 |
| Shannon entropy (bucketed) | 3.258 b | n/a | — |

First 30 partial quotients (for the record):
`[0, 2, 2, 2, 3, 12, 131, 1, 7, 1, 2, 1, 3, 3, 1, 2, 5, 39, 2, 1, 169,
 2, 2, 2, 1, 1, 2, 5, 11, 1, 2]`

## Verdict

**KHINTCHINE-TYPICAL.** The CF of the prime constant is statistically
indistinguishable from a "generic" real number in every test we ran.

Specifically, `a_n` looks like:
- a typical Gauss map orbit (GK distribution within sampling error),
- with no autocorrelation structure (max |ρ| at lags 1..10 is 0.013),
- a roughly white-noise power spectrum (β ≈ −0.09),
- no k-automaticity signal at k ≤ 10 (inter/intra variance ratios all
  < 0.012),
- entropy 3.26 bits per partial quotient (bucketed at 100), close to
  what GK predicts.

## Failure mode classification

**Equivalence (E).** Switching representations from binary digits
(primality indicator) to regular CF preserved entropy without exposing
any compressible structure — the CF is just another faithful
high-entropy encoding of the same data. No polylog shortcut available
through this representation.

## Caveats / what might change the verdict

- **Depth.** 1500 partial quotients is small. A k-automatic structure
  with period > 1500 could be invisible.
- **Specific CF.** Regular CF is one choice; Engel, Stern-Brocot, and
  Lüroth expansions are different and untested.
- **Conditional statistics.** Conditioning on a_n given e.g. n's
  divisibility class might expose structure invisible in the marginal.

None of these is a strong lead. Recording the negative result and
moving on.
