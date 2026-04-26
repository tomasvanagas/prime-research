# Eigenvector inspection — pointer

Auxiliary inspection script for the §D.D4 Szegedy-walk experiment.
The main results writeup is in `szegedy_walk_prime_graphs_results.md`.
Output of this specific script is in `eigvec_inspection.json` and
`eigvec_inspection.log`.

What this script confirmed:
- Top discriminant eigenvectors of the divisor graph at x=250 with
  eigenvalues in [0.70, 0.77] do NOT carry global primality content;
  each is concentrated on a 1-2-prime "cluster" `{p, 2p, 3p, …, 1}`.
- Specifically: k=14 (eigenvalue 0.768) localises 100% on
  {43, 47, 86, 94, 172, 188, …, 1}.
- This is the structural-localization result quoted in the main
  results.md and CLOSED_PATHS row.
