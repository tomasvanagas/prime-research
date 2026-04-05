# Sieve Matrix Rank: Results

**Script:** sieve_matrix_rank.py

## What Was Tested
Rank structure of the binary sieve matrix M[n][j] = 1 if p_j divides n: rank over R and GF(2), SVD spectrum, whether zero-row count (= primes + 1) can be extracted from rank/SVD structure without examining all rows, and transfer matrix approach.

## Key Findings
- Rank over R: equals pi(y) (number of sieving primes) -- full column rank, no deficiency to exploit
- Rank over GF(2): also full column rank = pi(y) for tested sizes
- SVD spectrum: singular values do NOT decay fast; no low-rank approximation captures the zero-row structure
- Zero-row count cannot be deduced from rank alone -- need explicit null-space enumeration which costs O(x * pi(y))
- Transfer matrix approach: per-prime projection matrices have size equal to primorial(y) -- exponential in pi(y)
- The sieve matrix is "maximally unstructured" for the purpose of counting zero rows

## Verdict
**CLOSED**
**Failure Mode:** Equivalence (E) -- sieve matrix has full column rank, non-decaying SVs, and exponential-size transfer matrices; no shortcut to zero-row counting.

## One-Line Summary
Sieve matrix rank analysis: full rank = pi(y), SVs don't decay, transfer matrix size = primorial -- no algebraic shortcut to prime counting.
