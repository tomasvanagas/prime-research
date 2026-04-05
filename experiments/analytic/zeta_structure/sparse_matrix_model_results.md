# Sparse Matrix Model for Zeta Zeros: Results

**Script:** sparse_matrix_model.py

## What Was Tested
Whether zeta zeros can be modeled as eigenvalues of a SPARSE Hermitian matrix with O(N) entries instead of O(N^2): (A) tridiagonal (Jacobi) matrix, (B) banded matrix with bandwidth k, (C) structured sparse (circulant+diagonal, Toeplitz, prime-based sparsity).

## Key Findings
- Tridiagonal (2N-1 parameters): per-eigenvalue relative error ~5% for N=100 -- poor
- Bandwidth 5: relative error ~1% -- better but grows with N
- Bandwidth 10: relative error ~0.3% for N=100, but ~1% for N=500 -- bandwidth must grow
- Minimum bandwidth for <1% error: grows as ~sqrt(N), giving O(N^{3/2}) parameters -- not sparse enough
- Toeplitz structure: ~2% error with N parameters; similar scaling issues
- Circulant+diagonal: captures mean spacing well but misses GUE fluctuations
- The GUE fluctuations require O(N^2) information to specify precisely (eigenvalue repulsion is non-local)
- No sparse matrix model achieves <0.1% per-eigenvalue error with O(N) parameters

## Verdict
**CLOSED**
**Failure Mode:** I (Information Loss -- sparse models lose GUE fluctuation information; bandwidth must grow as sqrt(N) for accuracy)

## One-Line Summary
Sparse matrix model for zeros: minimum bandwidth for <1% error grows as sqrt(N); GUE fluctuations require O(N^2) information.
