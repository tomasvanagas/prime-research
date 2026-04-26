I'll propose fresh approaches to computing p(n) in O(polylog(n)) time, save the detailed proposals, and build at least one runnable experiment to test the most promising idea.
Many proposals exist, so I need genuinely novel angles. Let me think fresh and propose ideas not in the existing list.
Now I'll write up 4 fresh proposals to `archive/ephemeral/proposals_session.md`, then build runnable experiments for the most testable two.
Now let me build the most testable experiments — Proposals C and B — and run them.
Proposal C: clean negative — 0/199 zero differences admit integer relations of bounded height in the chosen basis. That's a strong null. Now Proposal B.
Both proposals returned clean nulls:
- **C:** 0/199 zero differences admit relations.
- **B:** CF of prime constant is Khintchine-typical across every test (K₀=2.73 vs 2.69, Lévy=1.21 vs 1.19, Gauss-Kuzmin matches, autocorrelation < 0.013, flat spectrum, k-automaticity ratios all ≪ 1).

Now writing the companion result files.
## Summary

**Proposed 4 fresh approaches** in `archive/ephemeral/proposals_session.md`:

- **A — Tropical Hankel Rank of π** (max-plus rank of `H[i,j]=π(i+j)`)
- **B — CF of prime constant** (regular continued fraction of `α=Σ 2^{-p}`)
- **C — Algebraicity of Riemann zero differences** (PSLQ on `γ_k − γ_1`)
- **D — Reservoir computer fit on π(n)** (fading-memory test)

**Built and ran two experiments** (n < 10000 regime):

| Proposal | Verdict |
|---|---|
| C — zero algebraicity PSLQ (200 zeros, 12-element basis, height 10¹⁵, 50 dp) | **0/199 hits → NEGATIVE** |
| B — CF of α (1500 partial quotients, 12000 binary digits, 16800-bit precision) | **KHINTCHINE-TYPICAL** (K₀=2.73 vs 2.69; Lévy 1.21 vs 1.19; flat spectrum; no autocorrelation) |

Both close cleanly with high-information null results: zeros have no detectable Q-linear structure, and the CF representation of the prime constant carries the same entropy as the binary one with no exposed compression.

**Files created:**
- `archive/ephemeral/proposals_session.md`
- `experiments/proposals/proposal_c_zero_algebraicity.{py,_results.md}`
- `experiments/proposals/proposal_b_cf_prime_constant.{py,_results.md}`
- `.run_state` set to `50`.
