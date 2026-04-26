I'll generate proposals for fast prime computation. Let me start by loading needed tools and checking the directory structure.
Existing proposals already cover compressed sensing, PSLQ, and continued fractions. To stay fresh, I'll pivot to less-explored angles: hierarchical-matrix compression of the explicit-formula kernel, reservoir-computing prediction of delta-bits, determinant-representation hunt, and a modular-form / Hecke-trace sketch.

Let me write the proposals file first.
I'll append my fresh proposals to the existing session file rather than overwrite.
Now let me write the four runnable experiments. I'll start with B1 (H-matrix kernel), then B3 (determinant hunt), B5 (halving), B2 (ESN), in order of expected runtime.
The result is striking: kernel rank ≈ block width (linear fit beats log fit by 5x). H-matrix compression does NOT work on this kernel. Let me also test with more rows than columns (where H-matrix structure usually appears).
Crystal-clear negative: kernel rank = block width exactly. The H-matrix approach is fundamentally blocked because the columns are independent characters x^{1/2+i*gamma_k} (Vandermonde structure → no compressibility). Let me record results and move to B3.
Top candidates are essentially rediscovering `li(x) - log(x)·loglog(x)` — refined li, not a polylog representation. Let me probe extrapolation to larger x to confirm.
The exponent 0.40 is borderline interesting but R²=0.11 is weak. Let me extend the experiment with larger n and windowed RMS to get a cleaner exponent.
Asymptotic behavior of theta-residual is misleading: upper-window scaling is n^0.61, fully sqrt-chaos. The polylog look at low n is fitting artifact. Let me write up B5 honestly and move to B2.
ESN trained 64% but tested 50% — pure overfitting, no generalization. Let me confirm with a larger run.
All four experiments ran, results files written, no `__pycache__` to delete, `.run_state` set to 60. Final summary.

## Session summary

Five fresh proposals drafted in `archive/ephemeral/proposals_session.md` (appended as second batch — prior batch from earlier session preserved). Four runnable experiments implemented under `experiments/proposals/`:

| # | Idea | Verdict | Result |
|---|------|---------|--------|
| **B1** | H-matrix compression of explicit-formula kernel | **CLOSED** (structural) | Off-diagonal ranks = block width *exactly* (zero residual). The kernel `x^{1/2+i*gamma}` is Vandermonde-like; columns are independent characters. No off-diagonal cancellation possible. |
| **B3** | Determinant representation `det M(x) = pi(x)` | **CLOSED** (equivalence) | Best k=2 candidate over 12-feature library is `li(x) - log(x)·loglog(x)`, with sqrt(x)-scale error tracking li(x) - π(x) exactly. |
| **B5** | Self-correcting halving recursion `p(n) → p(n/2)` | **CLOSED** (info loss) | Asymptotic exponent of theta(n) residual is ~0.6 ≥ 0.5 — chaotic phases at scales n and n/2 add in RMS rather than cancel. Halving propagates, doesn't shrink, the entropy. |
| **B2** | ESN on bit-0 of delta(n) | **CLOSED** (no signal) | Train acc 67%, test acc 45-50% — pure overfitting. Logistic on same features also at chance. Confirms entropy barrier on delta parity. |
| **B4** | Hecke-trace coupling | sketch only (no code) | Left as future-arc proposal. |

Honest self-assessment: this session produced **4 clean negative closures** with concrete structural reasons, but **0 novel artifacts** (objects, compositions, proofs). By the project's own success criterion this is a "duplicate-shape" session — the fresh-frame proposals each map onto a known structural barrier within ~30 minutes of measurement. The B1 H-matrix result is the most pedagogically useful (clear Vandermonde argument explains why kernel-side compression fails, complementing prior signal-side compressed-sensing closures).

`.run_state` set to 60.
