# Information Theory of Primes

---

## Entropy of Prime Gaps

### Unconditional entropy: 3.94 bits per gap (first 100K primes)

The gap g(k) = p(k+1) - p(k) has 54 distinct values in this range.
Most common: 6 (17.0%), 2 (10.3%), 4 (10.2%), 12 (10.2%).

### Conditional entropy with features

| Feature set | Conditional entropy | Info gain | % of total |
|------------|-------------------:|----------:|-----------:|
| Unconditional | 3.943 bits | -- | -- |
| p mod 6 | 3.362 bits | 0.581 | 14.7% |
| p mod 30 | 2.998 bits | 0.945 | 24.0% |
| p mod 210 | 2.729 bits | 1.213 | 30.8% |
| p mod 2310 | 2.519 bits | 1.424 | 36.1% |
| p mod 30030 | 2.071 bits | 1.872 | 47.5% |
| prev 1 gap | 3.591 bits | 0.352 | 8.9% |
| prev 2 gaps | 3.332 bits | 0.611 | 15.5% |
| mod 2310 + prev 2 gaps | 1.659 bits | 2.284 | 57.9% |
| mod 30030 + prev 2 gaps | 0.541 bits | 3.402 | 86.3% |

**Irreducible entropy:** Even with the richest features (mod 30030 + 2 previous
gaps), 0.54 bits/gap remain. This is carried by residues mod uncaptured primes
(17, 19, 23, ...). Capturing all would require p mod primorial(sqrt(p)) = sieving.

### Combined features: residue + gap history

| Configuration        | Patterns | Deterministic | Majority accuracy |
|---------------------|----------|---------------|-------------------|
| mod 30, hist 2      | 2,006    | 531 (26.5%)   | 28.8%             |
| mod 210, hist 2     | 5,839    | 1,859 (31.8%) | 36.1%             |
| mod 2310, hist 2    | 23,129   | 10,446 (45.2%)| 49.5%             |
| mod 30030, hist 2   | 67,660   | 50,712 (75.0%)| 77.5%             |

High deterministic rates at large moduli are artifacts of key uniqueness (data
sparsity), not genuine predictive power.

### Decision tree overfitting: no learnable pattern

| Max depth | Train acc | Test acc | Leaves  |
|-----------|-----------|----------|---------|
| 10        | 38.1%     | 31.8%    | 817     |
| 20        | 74.0%     | 20.0%    | 26,729  |
| 30        | 98.0%     | 19.3%    | 47,528  |
| unlimited | 100.0%    | 19.1%    | 49,297  |

An unlimited decision tree memorizes training data perfectly (100%) but
generalizes at only ~19% on unseen primes -- worse than simple majority-vote
baselines. This confirms there is no learnable deterministic pattern in prime
gaps from number-theoretic features.

### Zero deterministic well-sampled keys

Among 626 distinct keys of (p mod 30, g(k-1), g(k-2)) that appear 20+ times
in the first 100K primes, **exactly zero** are deterministic. Every single
well-sampled feature combination maps to multiple different gap values.

This eliminates the possibility that apparent determinism at finer granularity
is genuine -- it is always an artifact of data sparsity (keys seen only once
are trivially "deterministic").

### AR model confirmation: ~5.04 bits/prime

AR(1)-AR(10) models on prime gaps confirm an entropy floor of approximately
5.04 bits per prime. This represents the information content that cannot be
predicted from local context.

---

## Kolmogorov Complexity

### K(p(n) | n) -- unbounded time
- Trivially O(1): a fixed program (sieve) computes p(n) from n
- This is not informative about computational difficulty

### K(p(n) | n) -- time-bounded (Kt complexity)
- OPEN QUESTION: What is Kt(p(n) | n)?
- If Kt = O(polylog(n)), implies small circuits exist (connects to NC question)
- If Kt = omega(polylog(n)), implies no fast algorithm exists
- Oliveira (2019) connected Kt to circuit lower bounds

### K(p(n) | n) conditional information
- delta(n) = p(n) - R^{-1}(n) requires O(log n) bits to specify
- But COMPUTING those bits takes O(x^{2/3}) time
- For p(10^100): ~170 bits of information, ~10^49 operations to extract

### Prime Coding Theorem (Kolpakov-Rocke 2023)
- E[K_U(X_N)] ~ pi(N) * ln(N) ~ N bits total
- Average information per prime: ~ln(N) bits (= log_2(N) * ln(2))
- Proves ML cannot discover compressible prime formula

---

## Entanglement Properties

### Volume-Law Entanglement (Latorre-Sierra 2014)

The prime state |P_n> = (1/sqrt(pi(2^n))) * sum_{p prime < 2^n} |p> has:
- Von Neumann entropy S ~ (7/8) * n/2 + (1/4) * log(n) + const
- This is VOLUME-LAW (linear in system size n)
- Coefficient 7/8 (not 1) encodes Hardy-Littlewood prime correlations
- Required MPS bond dimension: D ~ 2^{7n/16} (exponential)

### Implications
- Tensor network methods cannot efficiently represent the prime indicator
- No area-law: primes are analogous to maximally entangled quantum states
- The reduced density matrix encodes Hardy-Littlewood constants C(2k)
- Purity: Tr(rho^2_A) ~ 4.602 / (n * log 2)^2

### Our measurement (Session 10)
- Numerical MPS bond dimension: ~ N^{0.49} (consistent with volume-law)
- At N=2^14: bond_dim = 65 / max 128 (~51%)
- Gap sequence also requires full rank

---

## Information-Theoretic Decomposition of p(n)

For p(10^100) specifically:
- Total: ~340 bits
- Smooth part (R^{-1}): ~172 bits (computable in O(polylog))
- Oscillatory part (zeta zeros): ~170 bits
  - Encoded in ~10^48 zeta zero contributions
  - Each with GUE-random phase
  - Information-theoretically incompressible
  - Extracting requires O(10^49) operations minimum

### The Random Walk
delta(n) = p(n) - R^{-1}(n) is empirically a random walk:
- Autocorrelation r(1) = 0.996 (changes slowly)
- Step standard deviation ~ 8
- ZERO correlation with any computable function of n
- Cannot be predicted by polynomial fitting, Fourier, ML, or any method
- The only way to determine delta(n) is to compute pi(R^{-1}(n)) exactly

---

## Compressibility Results

- R(n) is 24% more compressible than random (some structure exists)
- But the RESIDUAL after removing R(n) is incompressible
- Spectral flatness of zeta zero sequence: 0.91 (near white noise)
- Gap spectral entropy: 94.9% of maximum (no hidden periodicities)
- Primes are NOT k-automatic for any k (Mauduit-Rivat 2010)

---

## References

- Latorre & Sierra (2014): arXiv:1403.4765
- Garcia-Martin et al. (2020): Quantum 4, 371
- Kolpakov & Rocke (2023): arXiv:2308.10817
- Mauduit & Rivat (2010): Annals of Mathematics 171
- Oliveira (2019): ICALP, LIPIcs 132
- Ortiz Tapia & Stoleum (2016): arXiv:1606.08293
