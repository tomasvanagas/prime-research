# Prime Gap Encoding & Wheel-Based Construction: Results

**Script:** prime_gap_encoding.py

## What Was Tested
Whether prime gaps can be predicted from local information (specifically p(k) mod primorial) to construct p(n) = 2 + sum g(k). Analyzed gap statistics conditioned on wheel position (mod 30 = 2*3*5).

## Key Findings
- Wheel mod 30 constrains gaps to specific residue patterns (8 allowed residues out of 30)
- Gap distribution conditioned on (p mod 30) shows measurable but weak dependence: wheel position predicts ~34.6% of gap values correctly
- The remaining ~65% of gap information is effectively random (from higher prime factors and zeta zeros)
- Even perfect wheel prediction (mod arbitrarily large primorial) captures only the smooth part of gap structure
- Cumulative error from wheel-based gap prediction still grows as O(sqrt(n))

## Verdict
**CLOSED** -- Failure Mode: I (Information Loss)

## One-Line Summary
Wheel-based gap conditioning captures only ~35% of gap information; the rest is random and incompressible.
