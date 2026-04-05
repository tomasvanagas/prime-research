# Quantum-Inspired Classical Algorithms: Results

**Script:** quantum_inspired.py

## What Was Tested
Six quantum-inspired classical approaches to computing p(n): (1) amplitude amplification analog via importance sampling near R_inv(n), (2) QFT analog via FFT on prime distribution / zeta zero quasi-periodicity, (3) phase estimation analog via power iteration on prime-counting operator, (4) tensor network / MPS as sieve via low-bond-dimension tensor contraction, (5) variational approach with parameterized function optimized via gradient descent, (6) adiabatic / homotopy continuation from li_inv(n) to p(n).

## Key Findings
- Importance sampling: concentrating samples near R_inv(n) does not help because delta(n) has no structure to exploit
- FFT on prime distribution: recovers known spectral structure (zeta zero frequencies) but cannot be inverted without O(sqrt(x)) terms
- Power iteration: dominant eigenvalue of prime-counting operator gives PNT; subdominant eigenvalues correspond to zeta zeros (same barrier)
- MPS sieve: bond dimension grows polynomially; cannot compress sieve into polylog representation
- Variational optimization: gradient descent on parameterized correction converges to R_inv + constant offset; cannot capture oscillations
- Homotopy continuation: path from li_inv(n) to p(n) passes through singularities at each prime gap; path-tracking requires O(n) steps

## Verdict
**CLOSED**
**Failure Mode:** E (Equivalence -- all quantum-inspired methods reduce to known classical methods: PNT approximation, explicit formula, or sieve)

## One-Line Summary
Six quantum-inspired classical methods (importance sampling, FFT, power iteration, MPS, variational, homotopy) all reduce to known barriers.
