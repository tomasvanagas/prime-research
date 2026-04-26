# Session 48 Wildcard Brainstorm — 5 angles not in the existing 57 wildcard scripts

Date: 2026-04-25
Constraint: must (a) not appear by name among existing wildcard .py scripts, (b) testable on n < 10⁵, (c) have a pass/fail criterion before running.

## 1. Causal-state / excess-entropy complexity of prime parity stream
The Berlekamp-Massey linear complexity of 1{n prime} over GF(2) is already
known to be N/2 (maximal). But linear complexity is a *deterministic*
measure: the "best" finite-state model is allowed to be linear feedback only.
The **excess entropy** E = lim_{L→∞} (H_L − L·h_μ) measures the *stochastic*
memory, i.e. mutual information between past and future. A finite-state
stochastic generator (epsilon-machine of finite size) ⇔ finite E.
- Pass: E plateaus as L grows → small stochastic state model exists.
- Fail: E grows linearly with L → no finite hidden Markov / generative model.
- Cost: O(N·2^L) for histograms, ≤ 30 s for N=10⁵, L=18.
- Status: this is the test I run.

## 2. Cumulant expansion of δ(n) = p(n) − R⁻¹(n)
δ(n) is "random-like" (21 measures). But Gaussianity is a 2-cumulant
statement. Compute κ₁,…,κ₆ for windows of n, check if higher cumulants
follow Gaussian scaling (κ_k = 0 for k≥3) or grow. If e.g. κ₄ scales
non-trivially in n, δ(n) is non-Gaussian and the 4th-cumulant has
arithmetic content extractable in polylog.
- Pass: κ_k for k≥3 has clean asymptotic in n.
- Fail: κ_k consistent with Gaussian + sampling noise.
- Why fresh: 21 measures focus on autocorrelation / spectrum / digit. Higher
  free cumulants of the marginal distribution have not been done.

## 3. Phase retrieval of pi(x) from |ζ(½+it)|² samples only
The magnitude |ζ(½+it)|² determines zero locations (zeros are zeros of |ζ|²)
but not phases. Phase retrieval theory (e.g. Wirtinger flow, sparsity-priors)
recovers a function from |F̂|² under genericity. Question: from O(polylog)
samples of |ζ|² on the critical line, can we approximate the zero-counting
function N(T) well enough for the explicit formula? Test on small T and
benchmark.
- Pass: NRMSE on N(T) decays sub-exponentially with #samples.
- Fail: requires Ω(T) samples (no improvement over direct enumeration).

## 4. Quasi-modular Eichler integral pairing for log ζ
log ζ(s) has known Eichler integrals relating it to weight-2 mock modular
forms (Bringmann-Folsom-Ono-Rolen 2017+). These have **algorithmically
computable** Fourier coefficients via Hecke operators. Question: does any
linear combination of mock-modular completions of log ζ admit a closed-form
that gives ψ(x) at integer x with no zero sum?
- Pass: identify a finite-dim space of mock modular forms whose
  evaluations encode π(n) at integer n via a fast algorithm.
- Fail: completions reduce back to the same explicit-formula sum.
- Note: highly speculative; would need lit search before coding.

## 5. Fast-multipole on **prime side** instead of zero side
The Riemann–Weil dual sum is over PRIMES, but with a fixed test function f.
For f localized near x, the sum Σ_p f(log p) determines π(x) up to small
correction. Group primes p₁<p₂<… into geometric clusters [2^k, 2^(k+1)).
Each cluster contributes Σ_{p in cluster} f(log p). Multipole expansion of f
around log(cluster center) compresses the cluster's contribution to O(d)
moments where d = degree of multipole. Question: for f a smooth cutoff, can
we evaluate Σ_p f(log p) in O(polylog · d) given we **already** had O(x^{2/3})
machinery to enumerate clusters?
- Pass: single-cluster evaluation reduces multi-prime sum to O(d) terms with
  d=O(log(1/ε)).
- Fail: still need to enumerate primes inside each cluster — nothing saved.
- Note: this is the *prime-side* mirror of `multipole_zero_sum.py`. Even if
  it fails, the **failure mode** comparison between zero-side and prime-side
  multipoles is informative.

## Selected for this session
**Idea 1** (causal-state / excess entropy). Cleanest pass/fail, fastest to
run, sharpens the pseudorandomness narrative if E grows linearly, and would
be a *new finding* if E plateaus.
