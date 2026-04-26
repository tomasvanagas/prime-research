# Proposal Session — Fresh Ideas (2026-04-26, session 63)

Four concrete proposals for computing p(n) in O(polylog). Each is distinct from
prior closed paths to my best knowledge (I deliberately did not consult
CLOSED_PATHS.md, per instructions). Each has a runnable test in
`experiments/proposals/session63fresh_*.py`.

---

## P1 — D-finite (Apéry-style) recurrence hunt for δ(n)

### Idea

Apéry famously showed ζ(3) admits a 2-term recurrence with polynomial
coefficients. PSLQ-style searches (already tried) look for *algebraic*
relations on δ(n) = p(n) − R⁻¹(n). A *D-finite* (holonomic) recurrence is
strictly more general:

    Σ_{j=0..L} P_j(n) · δ(n+j) = 0

where each P_j ∈ ℤ[n] has bounded degree. Holonomic sequences include all
hypergeometric, algebraic, and many transcendental sequences. If δ(n)
is D-finite of order L with degree-d coefficients, only O(L·d) integers
specify it — giving O(polylog(n)) evaluation via fast holonomic
extrapolation (Bostan-Salvy-Schost, ISSAC 2014).

### Pseudocode

```python
# Build a Krattenthaler-style ansatz matrix over ℤ:
# Rows: equations from forcing the recurrence at n=N0..N1
# Cols: unknowns are coefficients of P_0(n)..P_L(n) up to degree d
# Solve over ℚ; non-trivial null vector = recurrence
A = []
for n in range(N0, N1):
    row = []
    for j in range(L+1):
        for k in range(d+1):
            row.append(n**k * delta(n+j))
    A.append(row)
# rank-deficiency in A => D-finite recurrence
```

### Complexity & evaluation

If found at order L, degree d: evaluation of δ(n) is **O(L·d·M(log n))**
via Bostan's binary-splitting shift — genuinely polylog. Combined with
R⁻¹(n) (also polylog), p(n) lands in polylog.

### Key assumption

δ(n) is D-finite with small (L, d). Probably false (because the GUE-like
randomness inherited from ζ-zeros), but the test costs only K² linear
algebra over ℤ for K data points and is decisive.

### Test (n ≤ 10000)

Compute δ(n) exactly for n = 1..10000 via sieve. For each (L ∈ {1..4},
d ∈ {1..4}), build the ansatz matrix and check rank. Verdict if a
recurrence is found: extrapolate to n = 10001..10100 and verify against
sieve. **Strong falsification:** the matrix is full-rank at every (L,d).

---

## P2 — Mollifier-corrected explicit formula

### Idea

Selberg's *mollifier* trick: replace ζ(s) by ζ(s)·M(s) where M(s) is a
Dirichlet polynomial of length Y designed to cancel ζ's zeros in a
target window. The zeros of ζ·M include all zeros of ζ but with reduced
"amplitude" in the explicit formula. If we pick

    M(s) = Σ_{n ≤ Y} μ(n) P(log(Y/n)/log Y) n^{-s}

with a polynomial P (Conrey-Iwaniec-style), then the contribution of
zeros γ < T to π(x) is *attenuated* by a factor depending on M(½+iγ).

If we choose M to make M(½+iγ) ≈ 0 for the *first K zeros*, then we only
need to sum over zeros γ > T_K to compute π(x) accurately. The cost: the
mollifier is a one-time precomputation of size O(K log K).

### Pseudocode

```python
# Step 1: precompute first K zeros γ_1..γ_K
# Step 2: solve linear system for length-Y mollifier coefficients
#         to make |M(0.5 + i γ_j)|² minimal for j=1..K
# Step 3: explicit formula
#   π(x) ≈ R(x) - Σ_{γ > γ_K} R(x^ρ) · M(ρ)/M(1) + boundary terms
# Step 4: with M zeroing first K zeros, the remaining sum can be much shorter
```

### Complexity

If M zeros K zeros via a length-Y mollifier (Y ≤ x^θ), the effective
truncation is T_eff ≈ T₀ / K^α (Conrey-style mean-square estimate).
This *might* drop the zero-count from O(√x) to O(x^{1/2−ε}).
Polylog only if K can be made polylog while T_eff drops to polylog —
this is the conjecture.

### Key assumption

Cancellation extends *uniformly* in K with Y growing polynomially.
Conrey-style results establish this for K up to log T but not exponential.

### Test (n ≤ 10000)

For x = 100, 1000, 10000:
1. Compute baseline sum-over-zeros up to T = 1000.
2. Build mollifier with K = 5, 10, 20 zeroed and Y = 50, 100, 200.
3. Measure how much T can be reduced while keeping |π(x) − S_T(x)| < 0.5.
   If T = 100 suffices with K = 20, that's a 10× savings — encouraging.

---

## P3 — Random-Matrix-Theory moment-matching for δ(n)

### Idea

Conjecturally, ζ-zero spacings follow GUE statistics (Montgomery-Dyson).
The random variable Δ(x) = π(x) − li(x) over a random window [x, x+H]
has predictable *moments* (Cramér / Heath-Brown / Soundararajan-Young
predictions).

Hypothesis: by computing low-order moments of Δ(x) over a SHORT window
(cheap), we can localize π(x) within that window using Bayesian inference
*without* summing O(√x) zeros. Compute Δ(x) for x ∈ [X−H, X+H] using only
O(H · polylog X) primality tests (cheap if H = polylog(X)), and match to
GUE-predicted moments to recover the *expected value* of Δ(X).

### Pseudocode

```python
H = log(X)**3
# probe in window using fast primality testing (Miller-Rabin)
window_x = list(range(X - H, X + H))
li_vals = [li(x) for x in window_x]
pi_vals = [pi_local(x) for x in window_x]  # via per-integer prime test
delta_samples = [pi_vals[i] - li_vals[i] for i in range(len(window_x))]
# match to GUE-predicted moments
predicted_delta = bayesian_estimator(delta_samples, GUE_moments(X, H))
pi_X_estimate = li(X) + predicted_delta
```

### Complexity

Each prime test: O(log² X) via Miller-Rabin. H = polylog → polylog total.

### Key assumption

GUE predictions hold *pointwise* (not just in average). Soundararajan's
predictions are for the *tail* not pointwise extraction. Empirically,
GUE matches moments extremely well, providing a strong prior.

### Test (n ≤ 10000)

For X = 1000, 5000, 10000:
1. Compute Δ(x) for x in window of size H = 100.
2. Compute first 4 sample moments.
3. Use Wishart-style Bayesian update to predict Δ(X).
4. Compare to true π(X) − li(X).
5. Pass criterion: prediction within ±0.5 of true delta.

---

## P4 — Iterated Newton with progressive zero-budget

### Idea

p(n) is the inverse of π. Newton iteration converges quadratically:

    x_{k+1} = x_k − (π(x_k) − n) / π'(x_k)

If we use a *cheap approximate* π_K(x) using only K zeros, error is
~x/log(x) · e^{−cK/log x}. Newton converges quadratically, so digits
double per iteration. Hypothesis: at iteration k we only need
K_k = O(2^k) zeros to maintain error below the current Newton-residual
scale. Total zero work: K_total = Σ 2^k = 2 · K_final ≈ O(log n).

### Pseudocode

```python
x = R_inv(n)  # initial guess (polylog)
zeros = []  # zero cache, fetched on demand
for k in range(ceil(log2(log(n)))):
    needed_K = 2 ** (k+1)
    while len(zeros) < needed_K:
        zeros.append(next_zero())
    pi_approx = li(x) - sum(li_complex(x, rho) for rho in zeros[:needed_K])
    pi_prime = 1/log(x) - sum(...) # derivative
    x = x - (pi_approx - n) / pi_prime
return round(x)  # = p(n)
```

### Complexity

If `next_zero()` costs polylog AND zero-count needed is polylog, total
is polylog. The first is conditional on Riemann-Siegel formula running
in polylog per zero (currently O((zero_index)^{1/2}) — sub-polynomial
but not polylog).

### Key assumption

(a) Newton residual at iteration k is dominated by the *truncated zero
    tail*, which decays like K^{-α}.
(b) Polylog zeros suffice at the final precision.

(a) is reasonable: the explicit formula's tail behaves like 1/(γ log γ).

### Test (n ≤ 10000)

For n = 100, 1000, 5000, 10000:
1. Start at x = n · log(n).
2. At each Newton iteration, use K_k = 2, 4, 8, 16, ... zeros.
3. Track |x − p(n)| per iteration.
4. Pass: residual halves per iteration AND total zeros used is O(log² n).

---

## Summary of testability

| P | Test | Cost | Falsification |
|---|------|------|----------------|
| P1 | Linear algebra over ℤ on δ(1..10000) | sec | rank always full |
| P2 | Mollified zero sum at x=100,1000,10000 | min | T cannot be reduced |
| P3 | Bayesian update on H=100 samples | sec | prediction error ≈ full Δ |
| P4 | Newton with growing K | sec | residual stagnates |

Each test resolves in seconds-to-minutes. None requires GPU or cluster.
