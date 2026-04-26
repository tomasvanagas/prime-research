# Newton forward-difference series for pi(x)

## Hypothesis
If the forward-difference table Delta^k pi(0) decays or has compact
structure, then pi(x) = sum_{k>=0} C(x, k) * Delta^k pi(0) gives a
fast (polylog) evaluation by truncating to the first few terms.

## Setup
- Computed pi(0..256), then exact integer Delta^k pi(0) for k = 0..200.
- Tested truncated reconstruction at x = 10, 50, 100, 200 with K_use in
  {5, 10, 20, 40, 80, 160}.

## Findings

**Delta^k pi(0) GROWS EXPONENTIALLY** in k.
Linear fit of log2 |Delta^k pi(0)| vs k (over k = 1..200, omitting zeros):

    log2 |Delta^k pi(0)| ~ 0.9948 * k - 3.17

so |Delta^k pi(0)| ~ 1.99^k -- almost exactly 2^k. Concretely:
- k=8:  56
- k=16: 9918
- k=64: 2.16e18
- k=128: 4.05e37
- k=200: 2.04e59

**Truncation error is enormous.** Even with K_use = 160 terms, the
reconstruction error at x = 200 is ~ 10^91. Newton series for pi
**diverges** under truncation -- the alternating-sum cancellation that
yields the small integer pi(x) requires retaining all O(x) terms with
exact >= x bits of precision per coefficient.

## Verdict (failure mode I -- Information loss / wrong basis)

**Closed.** Newton-forward-difference series is a useless basis for
pi(x). Bit complexity is Omega(x): both the number of needed terms
and the bits per term grow linearly. This is **strictly worse** than
the trivial Eratosthenes O(x log log x) bound.

## Theoretical reason
A function with exponentially growing forward differences is, by a
standard result, the restriction to integers of an entire function of
exponential type >= log 2. The prime indicator chi_P, viewed as a
{0,1}-valued sequence, has no smooth analytic continuation -- it is
a "noise" function in the sense of difference theory.

## Connection to existing project results
This is consistent with E1.x (information-theoretic edges): pi has
Omega(x/log x) bits of incompressible content. The Newton basis is
agnostic to structure and inherits this bound directly.

## Files
- `newton_diffs_decay.py` -- experiment driver
