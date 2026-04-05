# Short-Interval Explicit Formula Analysis: Results

**Script:** short_interval_explicit.py

## What Was Tested
Cost analysis of counting primes in short intervals [x-W, x+W]: how many zeta zeros K are needed vs interval width W, and whether the scaling matches theoretical predictions. Tests whether short-interval counting has different zero requirements than full pi(x).

## Key Findings
- Full pi(x): needs K ~ sqrt(x) zeros for error < 1
- Short interval width W: needs K ~ x/W zeros (from sinc-function decay in Fourier domain)
- Empirically confirmed for x up to 10^6: K_min scales linearly with x/W
- Hybrid optimum: W = sqrt(x) gives K = sqrt(x), total cost O(sqrt(x) * polylog) -- same as direct
- Iteration does NOT help: each round needs approximately the same number of zeros
- Short-interval counting is cheaper per interval but the total work for identifying p(n) is unchanged
- The hybrid method (short contour + sieve) gives O(x^{1/2+eps}) with Huxley improvements, but not polylog
- No combination of short intervals can reduce total zero requirement below O(sqrt(x))

## Verdict
**CLOSED**
**Failure Mode:** E (Equivalence -- short-interval method trades interval width for zeros; total cost O(sqrt(x)) unchanged)

## One-Line Summary
Short-interval explicit formula: K ~ x/W zeros needed; hybrid optimum W=sqrt(x) gives O(sqrt(x)) total cost, same as direct.
