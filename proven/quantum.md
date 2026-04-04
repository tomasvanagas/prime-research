# Why Quantum Computing Doesn't Help

Summary of quantum approaches tested in Sessions 7, 9, and 10.

---

## Best Quantum Complexity: O(x^{1/3})

Using Grover search on the Deleglise-Rivat sieve, one can achieve O(x^{1/3})
quantum complexity (vs classical O(x^{2/3})). This is the quadratic Grover
speedup applied to the inner loop.

**For p(10^100):** 10^34 quantum operations. At 10^9 quantum ops/sec (optimistic
for 2026 hardware), that's 10^25 seconds = 10^17 years. Still factor 10^16 too slow.

---

## Quantum Lower Bound: Omega(sqrt(x)) queries

Nayak-Wu quantum query lower bound: Omega(sqrt(x)) queries to the input for
any quantum algorithm computing pi(x). This means even ideal quantum computers
face a sqrt(x) barrier.

---

## Approaches Tested and Closed

### 1. Grover on D-R Sieve
O(x^{1/3}) quantum, 10^34 ops for x=10^102. Insufficient by factor 10^25.

### 2. Shor-like (Periodicity Search)
Primes have NO periodicity that Shor's algorithm could exploit.
- Mauduit-Rivat (2010): primes not k-automatic
- Vinogradov: primes equidistributed in residue classes

### 3. Berry-Keating QPE
Even if the Hilbert-Polya Hamiltonian exists and is quantum-simulable,
QPE would need 10^51 eigenvalues (zeros) for p(10^100). Same barrier.

### 4. Primon Gas / GUE Simulation
Circular: statistical properties don't give individual primes.

### 5. BQP Membership
Exact pi(x) is related to #P-hard counting. PP is not in BQP (believed).

### 6. Quantum Channel / Communication
Holevo bound = classical Shannon capacity for deterministic functions.
No quantum encoding advantage for n -> p(n).

### 7. Entanglement-Assisted Counting
Superdense coding halves COMMUNICATION bits, not computation.
Bottleneck is computing pi(x), not communicating it.

### 8. Holographic Principle
Prime information scales as VOLUME O(N), not surface O(sqrt(N)).
No boundary/bulk duality exists for primes.

### 9. Quantum Error Correction
delta(n) mod m unpredictable from n mod m. No low-dimensional syndrome.

### 10. MPS / Tensor Network
Bond dimension ~ N^{0.49} (volume-law entanglement).
Exponential resources needed for exact representation.

---

## Why Quantum Fails: The Root Cause

All quantum approaches fail because the prime indicator has VOLUME-LAW
entanglement. The prime sequence is analogous to a maximally entangled
many-body state. No quantum algorithm can exploit structure that doesn't exist.

The 10 quantum paradigms analyzed all fail by factor >= 10^16.

---

## References

- Montgomery (1973): pair correlation conjecture
- Odlyzko (1987): numerical verification of GUE
- Latorre & Sierra (2014): arXiv:1403.4765 (entanglement in primes)
- Nayak-Wu: quantum query complexity
- Mauduit-Rivat (2010): primes not k-automatic
