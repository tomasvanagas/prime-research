# Degree-class check — pointer

Auxiliary script for the §D.D4 Szegedy-walk experiment.
Main writeup is in `szegedy_walk_prime_graphs_results.md`.
Output of this specific script is in `degree_class.log`.

What this script confirmed:
- High-ratio eigenvectors of the divisor graph (eigenvalues 0.70-0.77)
  concentrate mass on specific degree classes (degree-3 to degree-5
  vertices), NOT on primality globally.
- For x=250, k=14 eigenvector: 60% mass on degree-5 vertices, 20% on
  degree-4 vertices. The "high prime ratio" emerges because each
  cluster contains ONE prime + its multiples, so primes carry an
  outsized share of cluster mass — but every cluster picks its own
  small set of primes.
- This matches E7.12's mechanism (Cayley spectrum probes ω(n) /
  degree class, not primality).
