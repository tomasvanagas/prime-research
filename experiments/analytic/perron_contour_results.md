# Perron's Formula / Contour Integral Approach: Results

**Script:** perron_contour.py

## What Was Tested
Numerical contour integration of Perron's formula pi(x) = (1/2pi*i) integral of log(zeta(s)) * x^s/s ds. Methods: (1) Gauss-Legendre quadrature, (2) tanh-sinh (double exponential) quadrature, (3) saddle-point / steepest descent contour, (4) pole-subtracted Hankel contour, (5) residue sum (explicit formula) for comparison.

## Key Findings
- Gauss-Legendre: converges with O(T) quadrature points for height-T contour; total cost O(T^{4/3+eps}) using fast zeta evaluation
- Tanh-sinh: exponential convergence in number of quadrature points but contour height T still required
- Saddle-point contour: optimal contour shape reduces quadrature points ~2x but same T requirement
- Pole-subtracted: removes the s=1 pole contribution analytically; cleaner numerics but same asymptotics
- Residue sum (explicit formula): mathematically EQUIVALENT to the contour integral -- residues at zeros
- For error < 1: need contour height T ~ x^{1/2} regardless of quadrature method
- Total cost: O(x^{2/3+eps}) via Lagarias-Odlyzko, matching combinatorial methods
- Proved: contour integral and residue sum are mathematically identical; no "numerical shortcut" exists

## Verdict
**CLOSED**
**Failure Mode:** E (Equivalence -- Perron contour integral IS the explicit formula; residue expansion and contour integration are mathematically identical)

## One-Line Summary
Perron contour integral with various quadrature methods: mathematically equivalent to explicit formula; O(x^{2/3}) cost unavoidable.
