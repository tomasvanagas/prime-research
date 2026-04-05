# Algebraic Number Theory Approaches: Results

**Script:** algebraic_nt.py

## What Was Tested

Five algebraic number theory approaches to computing p(n): (1) CRT via splitting types in cyclotomic fields, (2) Artin L-function Frobenius extraction, (3) elliptic curve point counts (a_p recovery), (4) Hecke eigenvalue approach, (5) Iwasawa theory connection. Each analyzed for whether it avoids the O(x^{2/3}) or O(x^{1/2+eps}) barrier.

## Key Findings

- Cyclotomic CRT: Frobenius Frob_p at p in Gal(Q(zeta_q)/Q) is "multiplication by p mod q" -- requires knowing p first (circular)
- Artin L-functions: computing Frobenius elements requires factoring ideals, which requires knowing p
- Elliptic curve point counts: computing a_p = p+1-#E(F_p) requires knowing p
- Hecke eigenvalues: the eigenvalue of T_p on a newform requires knowing p to even define the operator
- Iwasawa theory: connects to p-adic L-functions via the explicit formula, which requires O(sqrt(x)) terms
- The Langlands program unifies all these L-functions, and the explicit formula barrier applies universally
- All 5 approaches are circular, equivalent to known methods, or encode only local/p-adic info

## Verdict

**CLOSED** -- Failure Mode: Circularity (C) / Equivalence (E). Every algebraic NT approach either requires knowing p first or reduces to the explicit formula with O(x^{1/2+eps}) cost.

## One-Line Summary

All 5 algebraic NT approaches (cyclotomic CRT, Artin L, elliptic curves, Hecke, Iwasawa) are circular or equivalent to the explicit formula barrier.
