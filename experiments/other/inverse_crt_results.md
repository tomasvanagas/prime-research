# Inverse CRT -- Constructing p(n) from Residues: Results

**Script:** inverse_crt.py

## What Was Tested
Whether p(n) mod q for small primes q can be computed without knowing p(n), enabling CRT reconstruction. Analyzed information content of p(n) mod q, entropy of residue distributions, and total bits needed.

## Key Findings
- p(n) mod 2 = 1 for n > 1 (free, 0 bits)
- p(n) mod q for q = 3, 5, 7, ... has near-maximal entropy (uniform over phi(q) residues)
- Total entropy from q=2..31 is ~32 bits; need q up to ~600 for 170 bits total
- No shortcut exists for computing p(n) mod q without knowing p(n) or pi(x)
- CRT needs ~52 primes as moduli (primorial > 10^102), each residue carrying ~log2(q-1) independent bits
- The total information content matches the barrier: ~170 bits needed, ~170 bits unavoidable

## Verdict
**CLOSED** -- Failure Mode: C (Circularity)

## One-Line Summary
Computing p(n) mod q requires knowing p(n); CRT reconstruction just redistributes the same information.
