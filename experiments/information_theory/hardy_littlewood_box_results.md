# Hardy-Littlewood {0,1}^k-cube Singular Series — Module Results

This module is a utility called by `gowers_uk_chi_p.py` and
`gowers_uk_chi_p_analyze.py`. It computes the singular series

  `S_k = product_p alpha_p(k) / (1 - 1/p)^{2^k}`

for the Hardy-Littlewood prime k-tuple conjecture applied to the
{0,1}^k-cube configuration `{x + epsilon . h : epsilon in {0,1}^k}`
in `Z^{k+1}`. `alpha_p(k)` is computed by direct numpy-vectorised
enumeration in `(Z/p)^{k+1}` checking the 2^k cube-form constraints.

## Numerical results

  `S_2 = 2.300938`   (P_max = 100, truncation error <= 1e-5)
  `S_3 = 54.116088`  (P_max = 47,  truncation error ~ 0.2%)

Convergence per prime is rapid; factor at p decays like 1 + O(1/p^{k+1}).
The dominant contribution at all k comes from p=2 (special case at the
{0,1}^k cube due to the cube embedding into F_2^k):

  `factor_2(k=2) = 2.0`        (singular series jump from p=2)
  `factor_2(k=3) = 16.0`

For full empirical comparison and falsification statement, see
`gowers_uk_chi_p_results.md`.
