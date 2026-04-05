# GF(2) ANF Structure and SLP Complexity of pi(x) mod 2 — Results

## Experiment (Session 35)

Deep analysis of the Algebraic Normal Form (ANF) over GF(2) for pi(x) mod 2,
going beyond Session 13's initial finding of degree Θ(N) and 50% sparsity.

## Key Findings

### 1. ANF Sparsity is Exactly 0.50 (Random)

| N | Monomials | Total possible | Sparsity |
|---|-----------|---------------|----------|
| 4 | 8 | 16 | 0.5000 |
| 6 | 32 | 64 | 0.5000 |
| 8 | 128 | 256 | 0.5000 |
| 10 | 512 | 1024 | 0.5000 |
| 12 | 2048 | 4096 | 0.5000 |
| 14 | 8192 | 16384 | 0.5000 |

**Exactly 50% of all possible monomials appear.** This is the expected value for a
random Boolean function. No deviation at any N tested.

### 2. Degree-N Monomial Has a Parity Pattern

| N | Max degree | Top monomial present? |
|---|-----------|----------------------|
| 4 | 4 (=N) | YES |
| 5 | 4 (=N-1) | NO |
| 6 | 6 (=N) | YES |
| 7 | 6 (=N-1) | NO |
| ... | ... | ... |
| 13 | 13 (=N) | YES |
| 14 | 14 (=N) | YES |

For even N: degree = N (the all-variables monomial is present).
For odd N ≤ 11: degree = N-1. For N = 13: degree = N.

The all-variables monomial coefficient equals Σ_x f(x) mod 2 = parity of |{x : pi(x) is odd}|.
This has a weak parity pattern for small N but breaks at N = 13.

### 3. CSE (Common Subexpression Elimination) — Same as Random

| N | Monomials | Naive ops | CSE ops | Savings | Random savings |
|---|-----------|-----------|---------|---------|---------------|
| 4 | 8 | 18 | 14 | 22.2% | ~22% |
| 6 | 32 | 98 | 62 | 36.7% | ~37% |
| 8 | 128 | 524 | 261 | 50.2% | 49.5% |
| 10 | 512 | 2597 | 1051 | 59.5% | 58.5% |
| 12 | 2048 | 12276 | 4206 | 65.7% | ~65% |

**CSE savings for pi(x) mod 2 match random functions within 1%.** There is NO
additional compressibility from number-theoretic structure. The SLP length after
CSE is Θ(2^N) — exponential in N.

**Direct comparison (N=10):**
- pi(x) mod 2: CSE = 1051 ops (savings 59.5%)
- Random function: CSE = 989 ops (savings 58.5%)
- MAJORITY: CSE = 566 ops (savings 65.0%)
- Inner Product: CSE = 54 ops (savings 0%)

MAJORITY has slightly BETTER CSE than random (because it has fewer monomials = 25% sparsity
and more structure). Inner Product has degree 2 so no CSE benefit. pi(x) mod 2 matches random exactly.

### 4. Degree-Truncated Accuracy Confirms adeg = N/2

| Degree cutoff | N=6 | N=8 | N=10 | N=12 |
|--------------|-----|-----|------|------|
| d ≤ 0 | 0.55 | 0.56 | 0.53 | 0.50 |
| d ≤ 1 | 0.55 | 0.51 | 0.51 | 0.50 |
| d ≤ 2 | 0.67 | 0.58 | 0.54 | 0.52 |
| d ≤ 3 | 0.86 | 0.68 | 0.58 | 0.53 |
| d ≤ N/2 | 0.95 | 0.82 | 0.81 | 0.81 |

At degree N/2: accuracy ≈ 0.80 (20% wrong). This is consistent with adeg(chi_P) = N/2
from Session 28. Half the Fourier weight is in degrees > N/2.

### 5. Variable Frequencies are Uniform

Mean variable frequency ≈ (monomials * degree) / N. Standard deviation / mean < 1%.
No variable (= no bit position of x) is more "important" than others in the ANF.

### 6. Monomial Patterns are Random

At each degree d, the fraction of present monomials approaches 0.50. The average
variable index of present monomials matches the expected value (N-1)/2. No systematic
bias toward high bits or low bits.

## Implications

1. **SLP complexity is Θ(2^N)**: Same as random, no compression from number-theoretic structure.
2. **No algebraic shortcut in GF(2)**: The polynomial structure is completely generic.
3. **Consistent with all prior measurements**: ANF sparsity (S13), degree (S13), 
   approximate degree (S28), and now SLP length all match random functions.

## Verdict

CLOSED: GF(2) SLP structure analysis. pi(x) mod 2 has NO exploitable GF(2) algebraic
structure. Its ANF is statistically indistinguishable from a random Boolean function
in all tested metrics: sparsity, degree, variable frequency, CSE savings, and monomial patterns.
