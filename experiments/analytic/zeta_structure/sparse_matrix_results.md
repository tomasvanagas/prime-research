# Sparse Matrix Model for Zeta Zeros

Date: 2026-04-05

Zeros loaded: 1000


## Approach A: Tridiagonal (Jacobi) Matrix

Parameters: 2N-1 (diagonal + off-diagonal)

- N=50: params=99, mean_rel_err=9.233910e-05, max_rel_err=7.007685e-04, rms=0.0180, converged=True, time=0.2s

- N=100: params=199, mean_rel_err=7.776646e-05, max_rel_err=8.794714e-04, rms=0.0246, converged=True, time=2.5s

- N=200: params=399, mean_rel_err=1.122295e-04, max_rel_err=4.380733e-04, rms=0.0424, converged=True, time=57.3s


## Approach B: Banded Matrix (varying bandwidth k)


### N = 50

- k=1: params=99, mean_rel_err=9.233910e-05, max_rel_err=7.007685e-04, rms=0.0180, converged=True, time=0.3s

- k=2: params=147, mean_rel_err=6.767566e-05, max_rel_err=3.840647e-04, rms=0.0122, converged=True, time=0.4s

- k=3: params=194, mean_rel_err=6.431410e-05, max_rel_err=2.942770e-04, rms=0.0115, converged=True, time=0.6s

- k=5: params=285, mean_rel_err=6.133813e-05, max_rel_err=3.530804e-04, rms=0.0119, converged=True, time=0.9s

- k=10: params=495, mean_rel_err=6.393095e-05, max_rel_err=3.676090e-04, rms=0.0121, converged=True, time=2.2s


### N = 100

- k=1: params=199, mean_rel_err=7.776646e-05, max_rel_err=8.794714e-04, rms=0.0246, converged=True, time=1.8s

- k=2: params=297, mean_rel_err=9.118409e-05, max_rel_err=1.172965e-03, rms=0.0315, converged=True, time=2.8s

- k=3: params=394, mean_rel_err=8.058719e-05, max_rel_err=1.191966e-03, rms=0.0311, converged=True, time=4.1s

- k=5: params=585, mean_rel_err=7.499008e-05, max_rel_err=1.194453e-03, rms=0.0305, converged=True, time=7.2s

- k=10: params=1045, mean_rel_err=9.857535e-05, max_rel_err=1.296899e-03, rms=0.0340, converged=False, time=10.2s


### N = 200

- k=1: params=399, mean_rel_err=1.122295e-04, max_rel_err=4.380733e-04, rms=0.0424, converged=True, time=23.5s

- k=2: params=597, mean_rel_err=7.820848e-05, max_rel_err=6.800513e-04, rms=0.0350, converged=True, time=43.6s

- k=3: params=794, mean_rel_err=1.051702e-04, max_rel_err=9.270158e-04, rms=0.0462, converged=True, time=37.0s

- k=5: params=1185, mean_rel_err=1.553818e-04, max_rel_err=1.213466e-03, rms=0.0626, converged=False, time=45.4s

- k=10: params=2145, mean_rel_err=2.853948e-04, max_rel_err=1.443806e-03, rms=0.0895, converged=False, time=46.4s


## Approach C: Structured Sparse


### C1: Circulant + Diagonal (2N params)

- N=50: params=100, mean_rel_err=0.000000e+00, max_rel_err=0.000000e+00, rms=0.0000, converged=True, time=0.1s

- N=100: params=200, mean_rel_err=0.000000e+00, max_rel_err=0.000000e+00, rms=0.0000, converged=True, time=0.4s


### C2: Toeplitz (N params)

- N=50: params=50, mean_rel_err=8.429653e-05, max_rel_err=3.055258e-04, rms=0.0118, converged=True, time=1.4s

- N=100: params=100, mean_rel_err=1.892630e-03, max_rel_err=3.575029e-02, rms=0.2181, converged=False, time=7.0s


### C3: Prime-based M_{ij} = a*log(pi*pj) + b*log(gcd) + c*delta*log(pi) + shift (4 params)

- N=50: params=4, mean_rel_err=1.482148e-01, rms=16.5496, converged=True, time=2.7s

  Best params: a=-0.1212, b=-1.1512, c=21.3571, shift=-3.1357


## Prediction Test

Fit on first N zeros, check fit quality on zeros N+1..N+50

- Jacobi N=50, predict next 50: train_rel_err=8.171396e-05, test_rel_err=7.381897e-05, time=1.5s

  (Note: this fits all N+50 zeros jointly, not a true extrapolation test.)

- Jacobi N=100, predict next 50: train_rel_err=7.710772e-05, test_rel_err=1.436003e-04, time=5.0s

  (Note: this fits all N+50 zeros jointly, not a true extrapolation test.)


### True Extrapolation (Jacobi)

Fit N zeros -> get 2N-1 params. Extend to (N+50) matrix by

extrapolating diagonal/off-diagonal trends, then check eigenvalues.

- N=50: extrapolated mean_rel_err=2.601347e-02, max_rel_err=5.477088e-02

  First 5 predicted vs actual:

    gamma_51: predicted=145.690007, actual=146.000982, rel_err=2.129952e-03

    gamma_52: predicted=147.892922, actual=147.422765, rel_err=3.189170e-03

    gamma_53: predicted=149.996637, actual=150.053520, rel_err=3.790843e-04

    gamma_54: predicted=152.094493, actual=150.925258, rel_err=7.747116e-03

    gamma_55: predicted=154.192208, actual=153.024694, rel_err=7.629581e-03

- N=100: extrapolated mean_rel_err=7.831206e-03, max_rel_err=1.923897e-02

  First 5 predicted vs actual:

    gamma_101: predicted=237.896493, actual=237.769820, rel_err=5.327515e-04

    gamma_102: predicted=239.603449, actual=239.555478, rel_err=2.002532e-04

    gamma_103: predicted=241.366833, actual=241.049158, rel_err=1.317884e-03

    gamma_104: predicted=243.134478, actual=242.823272, rel_err=1.281616e-03

    gamma_105: predicted=244.902252, actual=244.070898, rel_err=3.406195e-03


## Key Analysis: Minimum Bandwidth for <1% Error

- N=50: minimum bandwidth k=1 achieves mean_rel_err=9.233910e-05 (<1%)
  Parameters: 99 (vs N^2=2500 for dense)
  Compression ratio: 25.3x

- N=100: minimum bandwidth k=1 achieves mean_rel_err=7.776646e-05 (<1%)
  Parameters: 199 (vs N^2=10000 for dense)
  Compression ratio: 50.3x

- N=200: minimum bandwidth k=1 achieves mean_rel_err=1.122295e-04 (<1%)
  Parameters: 399 (vs N^2=40000 for dense)
  Compression ratio: 100.3x


## Analysis: Structure of Fitted Parameters


### N=50 Jacobi parameter structure

- Diagonal a_i vs gamma_i (trivial match): mean_rel_err=2.197530e-03

- Diagonal a_i ~ linear(i): slope=2.4792, intercept=26.8934, mean_rel_err=5.745744e-02

- Off-diagonal b_i ~ constant(1.3271): mean_rel_err=3.279195e-01

- Off-diagonal b_i ~ linear(i): slope=-0.021862, intercept=1.8518, mean_rel_err=2.438890e-01

- Off-diagonal stats: mean=1.3271, std=0.4744, min=0.4186, max=2.5850


### N=100 Jacobi parameter structure

- Diagonal a_i vs gamma_i (trivial match): mean_rel_err=1.186873e-03

- Diagonal a_i ~ linear(i): slope=2.1025, intercept=35.6459, mean_rel_err=6.523695e-02

- Off-diagonal b_i ~ constant(1.1406): mean_rel_err=3.369901e-01

- Off-diagonal b_i ~ linear(i): slope=-0.008768, intercept=1.5702, mean_rel_err=2.846172e-01

- Off-diagonal stats: mean=1.1406, std=0.4286, min=0.4035, max=2.5898


## Conclusions

### 1. Tridiagonal (Jacobi) fits are excellent -- but illusory

A tridiagonal Jacobi matrix with 2N-1 parameters achieves <0.05% mean relative
error for all tested N (50, 100, 200). Bandwidth k=1 is already sufficient for
<1% error. Increasing bandwidth beyond k=1 gives only marginal improvement.

However, this is **not compression**. The Jacobi matrix has 2N-1 free parameters
to fit N eigenvalues -- more parameters than targets. The optimizer has enough
degrees of freedom to fit any N numbers, not just zeta zeros. This is
overfitting, not structure discovery.

### 2. Circulant + Diagonal achieves zero error -- trivially

The circulant+diagonal model (2N params for N targets) achieves machine-zero
error because the diagonal alone can encode all N zeros. The circulant part
converges to zero. This confirms the overfitting concern.

### 3. Toeplitz (N params for N targets) is the real test

Toeplitz uses exactly N parameters for N eigenvalues -- the minimum for a
general set. It works well for N=50 (rel_err ~ 8e-5) but degrades for N=100
(rel_err ~ 0.2%, did not converge). The zeros become harder to fit with a
single structured family as N grows.

### 4. Prime-based matrix (4 params) fails completely

Mean relative error ~15%. A simple prime-based structure cannot capture zero
locations. Not surprising -- the connection between primes and zeros is through
the Euler product / explicit formula, not through a simple matrix entry formula.

### 5. True extrapolation fails

When we fit N zeros and extrapolate to predict the next 50:
- N=50: 2.6% mean error (max 5.5%)
- N=100: 0.78% mean error (max 1.9%)

The extrapolation improves with N but remains significantly worse than the
fitting error (~0.01%). The off-diagonal parameters b_i have std/mean ~ 35%
and do NOT follow a simple pattern (linear fit has ~25% residual error).
This means the b_i encode genuinely irregular information about zero spacings.

### 6. Key answer: bandwidth does not need to grow, but parameters do

Bandwidth k=1 (tridiagonal) suffices for <1% error at all N tested. But this
uses 2N-1 free parameters -- more than N. The "sparsity" is structural (O(N)
nonzeros vs O(N^2) in a dense matrix) but not informational (still O(N) free
numbers to store).

**The off-diagonal entries b_i are irregular with ~35% relative variation and
cannot be compressed to fewer than O(N) numbers.** They encode the same
GUE-random fluctuations that make the zeros incompressible. Repackaging the
zeros as Jacobi matrix parameters (a_i, b_i) does not reduce the information
content -- it merely changes the representation.

### 7. Implications for the prime research goal

This experiment confirms the information-theoretic barrier from a new angle:
- The zeros CAN be represented as eigenvalues of a sparse (tridiagonal) matrix.
- But the matrix entries themselves require O(N) bits of information.
- No low-parameter structured matrix family captures the zeros.
- This is consistent with delta(n) incompressibility and the known barriers.

**Status: CLOSED.** Sparse matrix representation does not provide a shortcut.
