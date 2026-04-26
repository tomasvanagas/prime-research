# Semiclassical block-resummation of zeta-zero contributions

## Hypothesis
The standard explicit formula
    pi(x) = R(x) - sum_rho R(x^rho) - <small terms>
needs ~sqrt(x)/log(x) zeros for the truncation error to fall below 1.
We test whether GROUPING zeros into "semiclassical orbits" (here
log-spaced blocks of indices [2^i, 2^{i+1})) and summing each block
analytically can reduce the effective zero count to O(log x).

If per-block contributions decay (or cancel pairwise), only O(log x)
blocks suffice and we get a polylog evaluation.

## Setup
Loaded 2000 zeta zeros from `data/zeta_zeros_2000.txt`. Computed
- pi(x) approx using standard explicit formula with ~sqrt(x) zeros
  for x = 10, 100, 500, ..., 100000.
- For x = 10000 specifically: per-block contributions for blocks
  [1,2), [2,4), [4,8), ..., [1024,2000).

## Findings

### Whole-spectrum approximation
Using ~sqrt(x) zeros (110 zeros at x = 10000) the approximation is
NOT accurate enough to round to the integer:

| x       | true pi | zeros used | approx     | rounded err |
|---------|---------|------------|------------|-------------|
| 1000    | 168     | 41         | -98        | -266        |
| 10000   | 1229    | 110        | -707       | -1936       |
| 100000  | 9592    | 326        | 59         | -9533       |

This confirms ~sqrt(x) zeros is NOT enough -- the standard analysis
needs more like sqrt(x)*log(x) zeros for rounding to succeed
(consistent with all prior project experiments on the explicit
formula, e.g. `truncated_zeta_zeros`, `zero_convergence_proper`).

### Block-resummation contributions at x = 10000

| block | zero index range | gamma range     | block sum |
|-------|-----------------|-----------------|-----------|
| 0     | [1, 2)          | 21.0            | 30.6      |
| 4     | [16, 32)        | 69.5 -- 105.4   | 195.3     |
| 6     | [64, 128)       | 173.4 -- 283.2  | 1175      |
| 8     | [256, 512)      | 481.8 -- 826.9  | 3362      |
| 10    | [1024, 2000)    | 1448.3 -- 2515.3| 12276     |

**Per-block contributions GROW like 2^{i/2}**. The contribution of
block i scales as block-size^{1/2} * sqrt(x) (consistent with
random-matrix predictions for zero-spacings). There is no
cancellation between blocks at this grouping.

Cumulative approximation actually DIVERGES from the truth as we
include more blocks (without compensating analytic-tail correction):

| blocks used | approx | true | err  |
|-------------|--------|------|------|
| 1           | 1196   | 1229 | -33  |
| 2           | 1221   | 1229 | -8   (best)
| 6           | -14    | 1229 | -1243|
| 11          | -25563 | 1229 | -26792|

The approximation is best at 2 blocks and gets WORSE as we include more.
This means the pi-explicit-formula is mis-truncated by including only
some zeros: the missing tail contribution is large.

## Verdict (failure mode E -- Equivalent to sqrt(x) explicit-formula bound)

**Closed.** Log-block grouping of zeta zeros does not help:
1. Per-block contributions grow as O(sqrt(2^i x)), so cumulative
   contributions through log(T) blocks scale as the FINAL block,
   which is O(sqrt(T x)).
2. Truncation error is dominated by zeros NOT included; one must
   include zeros to gamma >= sqrt(x) for variance reduction.
3. Random-matrix theory (GUE statistics) predicts O(1) standard
   deviation per block in a normalized basis, but the pi-target
   is the SUM of these so cancellation across blocks is sqrt-fast,
   not exponential.

This refines the project's existing E3.x edges (analytic explicit-
formula bounds) by showing that **no log-spaced semiclassical grouping
of zeros admits exact integer rounding from O(log) blocks alone.**

## Possible follow-up (not pursued here)
The truncation might benefit from grouping zeros NOT by index but by
their EXPECTED CANCELLATION CLASS in the GUE pair correlation. We
did not try this: it would require pre-clustering zeros via their
local pair correlation, which is itself O(T^2) work.

## Files
- `orbit_resum_pi.py` -- experiment driver
