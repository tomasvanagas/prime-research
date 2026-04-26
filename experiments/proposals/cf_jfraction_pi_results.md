# Results: J-fraction Expansion of Prime Generating Function (F2)

**Date:** 2026-04-25
**Script:** `cf_jfraction_pi.py`
**Conjecture tested:** PPC — pi(n) is P-recursive (satisfies a finite-order linear
recurrence with polynomial-in-n coefficients), equivalently the moment-generating
function of pi(n)/n is D-finite.

## Setup
- Computed pi(n) for n=1..80
- Used moments m_n = pi(n+2)/(n+2), n=0..78
- Computed Hankel determinants H_k of size (k+1)x(k+1) for k=0..16
- Extracted J-fraction coefficients (a_k, b_k) for k=0..11

## Key numbers
| k  | log10\|H_k\| |
|----|-------------|
| 0  | -0.301      |
| 4  | -4.665      |
| 8  | -8.516      |
| 12 | -15.853     |
| 16 | -21.192     |

Linear fit: log10|H_k| ≈ -1.32 k + 0.68
Quadratic fit: A_2 = -0.016 (essentially 0), so growth is **linear in k**, not quadratic.

| k  | a_k        | b_{k+1}     |
|----|------------|-------------|
| 0  |  0.333     | -0.778      |
| 1  | -2.105     | -0.029      |
| 2  | -0.600     | -3.000      |
| 3  | -2.790     | -8.283      |
| 4  |  9.729     | -2.142      |
| 5  | -7.062     |  0.140      |
| 8  | -0.217     | -0.034      |
| 11 |  0.885     | **-407.4**  |

## Verdict: NEGATIVE

Although log|H_k| decays linearly in k, this is a **trivially expected** consequence
of the moments being bounded in [0,1] — a uniform random sequence in [0,1] gives the
same. The linear log-decay does NOT imply D-finiteness; it would imply it only if
log|H_k| stayed BOUNDED or grew polynomially in *positive* direction.

The decisive evidence is in the **J-fraction coefficients themselves**:
- b_3 = -3, b_4 = -8.28, b_10 = -17.9, **b_12 = -407.4**
- a_k oscillates between -7 and +10 with no apparent pattern

For a D-finite sequence, the (a_k, b_k) themselves form a P-recursive sequence
(possibly periodic). Here they are clearly erratic, with one coefficient blowing up
two orders of magnitude. **pi(n) is not D-finite by any low-order recurrence
detectable from the first 80 values.**

## Failure mode

C — Circularity / Information Loss is not the issue. This is direct **Equivalence
failure**: the J-fraction representation captures pi(n) only by encoding all of its
random-like structure into the (a_k, b_k) sequence one-for-one — no compression.

## What would change the verdict

A clear pattern in (a_k, b_k) — e.g., a_k = (-1)^k * polynomial(k) / polynomial(k) —
would imply D-finiteness and a polylog algorithm. None observed. To rule out
D-finiteness rigorously one would need to extend to ~200 moments and apply
gfun-style detection (Maple); negligible expected return given current evidence.

## Next steps (if pursued)

- Try alternative moment families: m_n = pi(2^n), m_n = pi(n)^2/n, m_n = log p_n.
- Check Hankel determinants of (delta(n)) instead of pi(n)/n; if delta(n) is
  D-finite even though pi(n) is not, this would give a polylog algorithm.
- Cross-reference (a_k, b_k) with OEIS to detect known sequences.
- Compare to a baseline: same computation on a Cramér model of "pseudo-primes"
  to verify the chaotic (a_k, b_k) is specific to true primes vs. their model.
