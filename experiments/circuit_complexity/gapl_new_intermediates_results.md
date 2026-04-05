# GapL via New Intermediate Quantities: Results

**Script:** `gapl_new_intermediates.py`
**Session:** 16

## What Was Tested
Whether intermediate quantities OTHER than floor values and zeta zeros can serve as matrix entries in a poly-size det = pi(x) construction. Tested four families: class numbers h(-d), L-function values L(1, chi), elliptic curve point counts a_p, and number field regulators.

## Key Findings
- Class numbers h(-d): computing them requires knowing primes (circularity)
- L-function values L(1, chi): reduces to sums over primes/zeta zeros (equivalence)
- Elliptic curve a_p: point counts a_p = p+1 - #E(F_p) encode individual prime information; summing them gives L-function values (equivalence)
- Regulators: require factoring and computing units in number fields (circularity)
- ALL four families route back to circularity, equivalence to L-functions, or insufficient encoding

## Verdict
**CLOSED**
**Failure Mode:** Circularity / Equivalence (all four families reduce to known barriers)

## One-Line Summary
Class numbers, L-values, elliptic curve counts, and regulators all reduce to circularity or L-function equivalence for the GapL question.
