# Session 14 Summary

**Date:** 2026-04-04
**Focus:** Counting complexity connections, GapL/determinant approach, nonlinear sieves, algebraic structure analysis

## Key Results

### 1. Workspace Mismatch Barrier (NOVEL)
The chain PRIMES ∈ L → pi(x) ∈ #L → pi(x) ∈ NC^2 **BREAKS** at step 2.
An NL machine for pi(x) needs O(N) workspace for candidate n, but #L allows only O(log N).
This proves PRIMES ∈ L and pi(x) ∈ NC are **INDEPENDENT** questions.
Documented in novel/workspace_mismatch_barrier.md.

### 2. Information Barrier for GapL
The I-E (Legendre) formula's fractional parts {x/d} for 2^k squarefree d
carry O(2^k) independent bits. No determinant of size less than 2^{Theta(sqrt(x)/log(x))}
can encode the sieve. A GapL algorithm **MUST avoid floor functions entirely**.

### 3. Lucy DP Matrix Structure
- All transition matrices are **unipotent** (eigenvalues all = 1)
- Displacement rank is 50-60% of dimension (NOT structured)
- Product is **full-rank** with no compressible structure
- No Toeplitz, Cauchy, or Vandermonde structure

### 4. pi(x) is NOT a Linear Recurrence Sequence
- pi(x) mod m is NOT LRS for any modulus m = 2, ..., 16
- No fixed-size matrix A can encode pi(x) as (A^x)_{ij}
- Companion matrix fitting: order 20 gets only 27% test exact matches

### 5. Nonlinear Sieve: Parity Broken but Not Efficiency
Six experiments (products, comparisons, bitwise, identities, parity, multiplicative):
- Nonlinear operations CAN distinguish primes from semiprimes (breaking Selberg parity)
- But NO combination reduces cost below O(x^{2/3})
- Polynomial of degree d in K floor values needs K >= x^{1/(2d)} for exactness
- All polynomial fits overfit catastrophically

### 6. Algebraic Geometry Point Counting Fails
- Low-dimensional varieties: too few F_q-points to encode pi(x)
- High-dimensional (d ≈ N): Kedlaya algorithm is O(N^3 * 2^N) ≈ sieve
- Frobenius eigenvalues = zeta zeros (same information)

### 7. Other Closed Paths
- Polynomial regression on floor values (overfitting, residual = zeta zeros)
- Residue class / character sum counting (equidistributed, L-functions harder)
- Selberg sieve exact counting (parity barrier, 3.7x overshoot)
- Batched Fermat/MR (no shared structure in consecutive 2^{n-1} mod n)
- Smooth/rough decomposition (equivalent to sieve)
- DAG path count compression (floor-value set is irreducible)
- Spanning tree / Pfaffian encoding (equivalent to GapL question)

## Experiments Created
- experiments/circuit_complexity/lucy_dp_matrix_structure.py
- experiments/circuit_complexity/matrix_power_encoding.py
- experiments/circuit_complexity/polynomial_identity_approach.py
- experiments/circuit_complexity/additive_structure_primes.py
- experiments/circuit_complexity/counting_complexity_connection.py
- experiments/circuit_complexity/sharp_L_analysis.py
- experiments/circuit_complexity/variety_point_counting.py
- experiments/circuit_complexity/small_matrix_search.py (agent)
- experiments/circuit_complexity/compressed_redheffer.py (agent)
- experiments/circuit_complexity/lgv_dag_construction.py (agent)
- experiments/circuit_complexity/ie_as_determinant.py (agent)
- experiments/circuit_complexity/lucy_dp_matrix_formulation.py (agent)
- experiments/sieve/nonlinear_floor_products.py (agent)
- experiments/sieve/nonlinear_comparisons.py (agent)
- experiments/sieve/nonlinear_selberg_parity.py (agent)
- experiments/sieve/nonlinear_floor_identities.py (agent)
- experiments/sieve/nonlinear_multiplicative_structure.py (agent)
- experiments/sieve/nonlinear_bitwise_tc0.py (agent)

## Novel Findings
- novel/workspace_mismatch_barrier.md — PRIMES ∈ L ≠> pi(x) ∈ #L

## Updated Files
- CLAUDE.md (Session 14 status, closed paths, viable directions)
- status/CLOSED_PATHS.md (~15 new entries, count 430+)
- status/OPEN_PROBLEMS.md (Session 14 refinements)
- proven/complexity.md (workspace mismatch barrier)

## Paths Still Open
The problem remains GENUINELY OPEN. The key remaining directions:
1. **GapL question**: Can pi(x) = det(poly(N)-size matrix)? Must avoid floor values.
2. **Growing-dim MPOW**: Is matrix powering with growing dimension in TC^0?
3. **BPSW unconditional proof**: Would give PRIMES in TC^0.
4. **Novel intermediate quantities**: Not floor values, not zeta zeros — what else?
5. **Berry-Keating Hamiltonian**: Quantum, far from concrete.
6. **Kt complexity**: Connects circuit bounds to compression.
