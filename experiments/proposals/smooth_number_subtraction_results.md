# Smooth Number Counting + Subtraction for pi(x) -- Results

**Script:** smooth_number_subtraction.py

## What Was Tested
Can counting B-smooth numbers Psi(x, B) be done faster than counting primes? If so, reconstruct pi(x) from Psi values via linear combination or Buchstab identity. Tests: exact Psi(x, B) computation, linear combination fitting, timing comparison, Buchstab recursion tree size.

## Key Findings
- Psi(x, B) computed exactly for various B and x up to 10000. Exact computation via factoring each integer is O(x * sqrt(x)) -- slower than sieve.
- Sieve-based Psi(x, B) computation: mark multiples of primes > B. Cost O(x * pi(x)/log(x)) -- comparable to Eratosthenes.
- Linear combination: pi(x) = sum c_i * Psi(x, B_i) has no exact finite solution. The Buchstab identity gives pi(x) in terms of Psi values, but the recursion tree has 2^{pi(sqrt(x))} nodes.
- Memoized Buchstab recursion IS the Lucy DP algorithm, costing O(x^{2/3}). No shortcut from the smooth-number angle.
- Dickman's rho function gives Psi(x, B) ~ x * rho(u) where u = log(x)/log(B), but this is an APPROXIMATION with error O(x/log(x)).
- The correction from approximate Psi to exact pi(x) encodes the same zeta-zero information.
- Previously closed: "Smooth number counts" (S6), "Buchstab tree pruning" (S11), "Buchstab signed identity" (S16), "Recursive Dickman DDE shortcut" (S24).

## Verdict
**CLOSED**
**Failure Mode:** Equivalence (E) -- memoized Buchstab recursion IS Lucy DP at O(x^{2/3}); smooth number counting is not fundamentally easier.

## One-Line Summary
Smooth number subtraction: Buchstab recursion tree is exponential; memoized version IS Lucy DP at O(x^{2/3}); no advantage.
