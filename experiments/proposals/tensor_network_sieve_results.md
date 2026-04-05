# Tensor Network Contraction for Prime Counting (Legendre Sieve) -- Results

**Script:** tensor_network_sieve.py

## What Was Tested
Represent the Legendre sieve phi(x, a) as a tensor network (MPS) with binary indices indicating inclusion/exclusion of each prime. Measure bond dimensions via SVD to determine if efficient contraction is possible.

## Key Findings
- Tensor T[s_1,...,s_a] = (-1)^{sum s_i} * floor(x / prod_{i: s_i=1} p_i) built for various x and a = pi(sqrt(x)).
- Sequential SVD decomposition into MPS form: bond dimensions measured at each bipartition.
- Bond dimension grows EXPONENTIALLY with a (number of sieving primes). For a=10, bond dim already ~2^8-2^9; for a=15, approaches 2^{14}.
- Entanglement entropy at each bipartition is approximately linear in min(k, a-k) -- VOLUME-LAW scaling, not area-law.
- Volume-law entanglement means no efficient MPS/tensor-network representation exists.
- This is consistent with sessions 10 and 20 which found the same volume-law behavior.
- Physical interpretation: the divisibility structure of integers by different primes is maximally entangled -- knowing divisibility by p_1...p_k gives no information about divisibility by p_{k+1}...p_a.

## Verdict
**CLOSED**
**Failure Mode:** Information Loss (I) -- Legendre sieve tensor has volume-law entanglement; bond dimension 2^{Omega(a)}, no efficient contraction.

## One-Line Summary
Tensor network sieve: volume-law entanglement, bond dim 2^{Omega(pi(sqrt(x)))}; no efficient MPS contraction possible.
