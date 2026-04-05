# Commutative Matrix Powering (AKS Special Case): Results

**Script:** `commutative_matrix_powering.py`
**Session:** 11

## What Was Tested
Whether the commutativity of the AKS polynomial ring Z_n[x]/(x^r - 1) allows matrix powering to be done in TC^0. The AKS primality test requires (x+a)^n mod (x^r - 1, n), which is matrix powering in a commutative ring. Scalar powering is TC^0 (Allender 1999); does commutativity extend this?

## Key Findings
- For scalar powering, the key trick is reducing the exponent mod phi(m) -- works because the group order is known
- For commutative ring powering, the group order |R*| depends on the factorization of x^r - 1 over Z_n, which depends on whether n is prime -- circular!
- The "reduce exponent" step requires knowing if n is prime, which is what AKS is trying to determine
- Even with commutativity, the ring dimension r = O(log^c n) makes the matrix exponential in the relevant parameter

## Verdict
**CLOSED**
**Failure Mode:** Circularity (computing the group order of Z_n[x]/(x^r-1) requires knowing if n is prime)

## One-Line Summary
Commutative ring powering for AKS is circular: reducing the exponent mod |R*| requires knowing primality of n.
