# Class Field Tower Prime Reconstruction -- Results

**Script:** class_field_splitting.py

## What Was Tested
Use splitting patterns of primes across quadratic fields Q(sqrt(d)) (via Legendre symbols) to determine p(n) mod M for large M. If pi_split(x, K) (counting primes with specific splitting pattern across K fields) is easier than pi(x), this could enable CRT reconstruction.

## Key Findings
- A prime p splits in Q(sqrt(d)) iff Legendre(d, p) = 1. Each quadratic field gives ~1 bit of information about p.
- With m fundamental discriminants, the splitting pattern determines p mod M for M ~ 2^m (each field halves the candidates).
- However, computing pi_split(x, K) -- counting primes up to x with a specific splitting pattern across K fields -- is NOT easier than pi(x). It requires the Chebotarev density theorem, which involves L-function zeros.
- Each quadratic field contributes a Hecke L-function L(s, chi_d) with its own zeros. Computing the split-count to unit accuracy requires summing over sqrt(x) zeros of each L-function.
- Previously closed: "Chebotarev density theorem" (session 4, failure C) and "Number field algebraic decomposition" (session 7, failure E).
- The splitting pattern provides information about p mod discriminants, but extracting that information computationally has the same L-function zero cost.

## Verdict
**CLOSED**
**Failure Mode:** Equivalence (E) -- pi_split(x, K) requires L-function zeros at O(sqrt(x)) cost per field; no easier than pi(x).

## One-Line Summary
Class field splitting patterns: each field gives 1 bit but costs O(sqrt(x)) via L-function zeros; no faster than direct pi(x).
