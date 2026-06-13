# Results: p11_top_lines_quadratic_decomp.py

Quadratic decomposition tool for the top-K dominant Ulam-spiral
lines. For each input canonical line key (a, b, c), enumerates the
integers v ∈ [1, N] whose Ulam coords (X(v), Y(v)) lie on that line,
and verifies that the resulting integer sequence decomposes into
pairs of quadratic forms `f(k) = 4k² + bk + c` (one per Ulam side).

## Confirmed quadratic decomposition (N=10⁵)

For top line (1, -1, -18) → x+y=18 line on the Ulam plane:
- 299 total Ulam-points on this line; 121 of them prime → density 0.405
- Sorted integer sequence: 307, 379, 383, 459, 467, 547, 559, 643, ...
- Splits into 2 interleaved quadratics:
  - Even-indexed: 1st-differences 76, 84, 92, 100, 108 (2nd-diff = 8)
    → leading coef 4 → quadratic `4k² + bk + c`
  - Odd-indexed: 1st-differences 80, 88, 96, 104, 112 (2nd-diff = 8)
    → leading coef 4 → another quadratic `4k² + b'k + c'`

All three top lines (intercepts -18, -58, +40) verified to decompose
into 2 quadratics with leading coefficient 4 and second-differences
exactly 8. This is structurally consistent with the standard fact
that each Ulam-spiral line is a union of polynomial sequences of
degree 2.

The two quadratics per line correspond to the two "sides" of the
Ulam ring intersecting the line. For x+y=18 specifically, these are
the east side (`x = k, y = 18-k` for k ≥ 9) and the north side
(`y = k, x = 18-k` for k ≥ 10), giving formulas like
`f₁(k) = 4k² - 4k + 19` and `f₂(k) = 4k² - 17`.

## Cross-reference

See `p11_ulam_line_cover_results.md` for the connection to
Hardy-Littlewood prime-rich quadratic forms. Top-line densities
~40% match HL Conjecture F predictions for small-discriminant
class-number-1 quadratics.

## What would falsify

If a top-K dominant line **fails** the `4k² + bk + c` decomposition,
that would indicate either (a) a bug in the Ulam coordinate
function, or (b) a previously-unknown structure on the Ulam spiral
(unlikely — this is well-established). Both top-3 lines verified at
N=10⁵.
