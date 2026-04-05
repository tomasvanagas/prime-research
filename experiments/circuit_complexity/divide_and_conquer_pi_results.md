# Divide-and-Conquer Structure of pi(x): Results

**Script:** `divide_and_conquer_pi.py`
**Session:** 15

## What Was Tested
Whether pi(x) has a recursive decomposition with polylog overhead. Analyzed binary recursion pi(x) = pi(x/2) + delta(x), Buchstab-like identities, and mutual information between pi(x) and pi(x/2).

## Key Findings
- delta(x) = pi(x) - pi(x/2) has ~50% unpredictable bits relative to PNT or Li predictions
- The error in delta(x) relative to Li(x) - Li(x/2) encodes the same oscillatory zeta-zero information
- Buchstab identity unfolds into a tree of depth O(log log x) but total work is still O(x^{2/3})
- Mutual information I(pi(x); pi(x/2)) is high but the conditional entropy H(pi(x) | pi(x/2)) grows as sqrt(x)/log(x)
- No recursive scheme achieves polylog overhead per level

## Verdict
**CLOSED**
**Failure Mode:** Information loss (each recursion level requires sqrt-scale new information from zeta zeros)

## One-Line Summary
Divide-and-conquer for pi(x) fails: each recursion level requires O(sqrt(x)/log x) new bits of zeta-zero information.
