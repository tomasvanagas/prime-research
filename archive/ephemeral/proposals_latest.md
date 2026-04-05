I'll search for recent papers and techniques in parallel, then develop concrete proposals.
Good research done. Now let me develop the proposals and write test code.
Let me check what's in the data directory for zeta zeros, and look at the existing best algorithm for reference.
Now let me write all four proposals as testable Python code. I'll create them in parallel.
All four proposal scripts written. Now let me run them all in parallel.
The zeros aren't loading. Let me debug the file path issue and fix the scripts.
The file has no index column — just one value per line. Let me re-read to confirm.
The file has one value per line (no index column — the Read tool adds line numbers). The parser expects 2 columns but finds 1. Let me fix all scripts.
Now run all four in parallel:
Good results from proposals 21 and 22. Let me fix the bugs in 23 and run 23+24.
All four proposals ran successfully. Very interesting results. Let me now write the results files and the session proposal document.
Now let me write the session proposals document.
Now clean up and save run state.
## Session 36 Complete — 4 Proposals Tested

All four proposals were developed, coded, tested on n < 10,000, and documented. Here's the executive summary:

### Results

| # | Proposal | Verdict |
|---|---------|---------|
| **21** | Zero Clustering Truncation | CLOSED — tail contributions unpredictable; each zero carries independent info |
| **22** | Compressed Sensing on delta(n) | CLOSED — delta is NOT sparse in any tested basis; 11-16% extrapolation accuracy |
| **23** | PSLQ Integer Relations | CLOSED — normalized delta passes Gaussian normality test; no universal formula |
| **24** | Dequantized Grover Sieve | CLOSED — primality operator rank ~ N^{0.365}, not polylog |

### Key New Insights

1. **10 zeros give O(1) accuracy** but ~1.3 average error — not enough to pin down the exact prime. The "last mile" is the hard part.

2. **Normalized delta(n) is Gaussian** (p=0.19 normality test). This is strong evidence the residual is genuinely random, not just hard to compute.

3. **Sieve matrix rank scales as N^{0.365}** — any rank-based approach is bounded away from polylog.

4. **delta(n) is dense in all tested bases** — Fourier needs 50+ coefficients for 87% energy, zeta zero basis is 69% dense at 100 terms.

### Files Created
- 4 experiment scripts in `experiments/proposals/` (proposals 21-24)
- 4 companion `_results.md` files
- Session summary in `archive/ephemeral/proposals_session.md`
