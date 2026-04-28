# eigenfunction_analysis.py — leading-eigenvector / trace identity analysis

For each `h ∈ {unweighted, χ_P, λ, Λ}`, extracts the leading right
eigenvector of `L_h` and computes:
1. cosine overlap with the unweighted Gauss-Kuzmin invariant density
   `g(x) = 1/((1+x) log 2)`;
2. traces `tr(L_h^k)` for `k = 1..6`;
3. fixed-point trace prediction `Σ_n h(n)/(n²+n-1)` (Mayer-style);
4. arithmetic invariants for χ_P (prime zeta `Σ_p 1/p²`, prime-fixed-
   point sum `Σ_p 1/(p²+p-1)`, kernel-sup sum `Σ_p 1/(p(p+1))`).

**Headline:** leading right eigenvector overlaps Gauss-Kuzmin density
at cosine 1.000 (unweighted, sanity ✓), 0.99853 (χ_P), 0.99191 (λ),
0.99525 (Λ). **The arithmetic content of all three weighted operators
lives entirely in the eigenvalue, not in the eigenfunction.** This
justifies the Rayleigh-quotient closed-form prediction implemented in
`closed_form_prediction.py`.

For full discussion + numerics see `pollicott_ruelle_chi_p_results.md`.
