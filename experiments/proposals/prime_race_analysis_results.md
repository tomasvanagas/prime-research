# Prime Race Analysis: Information Content of E(x;q) -- Results

**Script:** prime_race_analysis.py

## What Was Tested
Can the prime race difference E(x;q) = pi(x;q,a) - pi(x;q,b) be computed more cheaply than pi(x)? Experiments: (a) compute E(x;4) for x up to 10^6; (b) spectral analysis comparing smoothness of E(x;4) vs pi(x)-Li(x); (c) Chebyshev bias |E(x;4)|/sqrt(x); (d) information content of p(n) mod q for q=3,5,7,11,13.

## Key Findings
- E(x;4) = pi(x;4,3) - pi(x;4,1) computed exactly via sieve up to 10^6. Exhibits the classic Chebyshev bias (E > 0 most of the time).
- Spectral analysis: E(x;4) has the SAME spectral complexity as pi(x) - Li(x). The Fourier spectrum shows no additional sparsity -- both are dominated by zeta/L-function zeros.
- Chebyshev bias: |E(x;4)|/sqrt(x) ~ O(1/log(x)), consistent with the Rubinstein-Sarnak analysis. The bias is a STATISTICAL property, not exploitable for exact computation.
- Information content of p(n) mod q: each mod-q value carries ~log2(phi(q)) bits. For q=3: ~1 bit, q=5: ~2 bits, etc. Total across many q matches the ~log(n) bits needed.
- Computing p(n) mod q to certainty requires L-function zeros for chi mod q at O(sqrt(x)) cost -- no savings over direct pi(x).
- The prime race is a CONSEQUENCE of L-function zeros, not a shortcut around them.

## Verdict
**CLOSED**
**Failure Mode:** Equivalence (E) -- E(x;q) has the same spectral complexity as pi(x); prime races are consequences of L-function zeros, not shortcuts.

## One-Line Summary
Prime race E(x;q): same spectral complexity as pi(x)-Li(x); Chebyshev bias is statistical, not computationally exploitable.
