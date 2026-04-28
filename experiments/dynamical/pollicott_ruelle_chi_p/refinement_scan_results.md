# refinement_scan.py — Pollicott-Ruelle resonance stability scan

Sweeps `(M_grid, n_max) ∈ {30, 60, 90, 120, 160} × {100, 200, 400, 800}`
for each weight `h ∈ {unweighted, χ_P, λ, Λ}`. Reports top-5 eigenvalues
per cell and the `L^∞` stability gap across cells. Output:
`refinement_scan.json`.

**Headline:** Top-5 eigenvalues stable to <1% coefficient-of-variation
across the full 20-cell sweep for `n_max ≥ 100`, confirming real
Pollicott-Ruelle resonances (not spurious-discretisation eigenvalues).

For full discussion + numerics see `pollicott_ruelle_chi_p_results.md`
(parent results file).
