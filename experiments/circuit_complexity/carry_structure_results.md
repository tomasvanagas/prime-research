# Binary Carry Structure of pi(x): Results

**Script:** `carry_structure.py`
**Session:** 17

## What Was Tested
Five analyses of the binary arithmetic structure of pi(x) = sum 1_P(k): carry propagation patterns, per-bit Boolean complexity, communication matrix rank for each bit, threshold [pi(x) >= t] complexity, and comparison with random subsets.

## Key Findings
- Carry propagation is irregular and unpredictable; no exploitable pattern in carry chains
- Individual bits of pi(x) have communication rank close to full (exponential)
- Threshold functions [pi(x) >= t] are monotone but have exponential monotone circuit complexity
- All measures match random functions of the same density within statistical error
- No shortcut visible from the binary accumulation perspective

## Verdict
**CLOSED**
**Failure Mode:** Information loss (binary carry structure matches random functions)

## One-Line Summary
Binary carry/accumulation structure of pi(x) is random-like; no exploitable pattern in carry propagation or per-bit complexity.
