# Local-to-Global Reconstruction for pi(x)

## Experiments Run
5 experiments testing whether local/cheap computations can reconstruct pi(x).

## Key Results

### Exp 1: Partial Sieve Residual
Sieving by primes ≤ B leaves false positives. For x = 1000:
- B = 7: 63 false positives (37% of pi(x))
- B = 23: 3 false positives
- B = 31 (= √1000): 0 false positives → exact
Need B ≈ √x always. **No shortcut.**

### Exp 2: Randomized Inclusion-Exclusion
Random subset sampling of IE requires ~370,000 samples for x = 500 (k = 8 sieve primes).
Full IE has only 256 terms. **Randomization makes it WORSE.**

### Exp 3: Error Autocorrelation ★
E(x) = pi(x) - Li(x) has lag-1 autocorrelation = **0.996** — extremely correlated at nearby points!
- std(E(x+1) - E(x)) = 0.326 (very predictable locally)
- But only helps for adjacent queries, not single-point computation

### Exp 4: Information Content ★★★ KEY FINDING
| x | |E(x)| | bits(E) | (log₂x)² |
|---|--------|---------|-----------|
| 100 | 4 | 3.2 | 44 |
| 1000 | 9 | 4.2 | 99 |
| 10000 | 16 | 5.1 | 177 |

**E(x) has MUCH LESS information than polylog!** Only O(log x) bits.

Scaling for larger x (on RH):
- x = 10^100: |E| ≈ 10^48, needs ~160 bits. (log₂x)² ≈ 110,000.
- Ratio: 160 / 110,000 = 0.15% of polylog budget

**There is NO information-theoretic barrier to polylog computation of pi(x).**

### Exp 5: Multi-Resolution with Zeta Zeros
30 zeros insufficient for any tested x (crude Li(x^ρ) approximation used).
Individual zeros contribute < 1 bit each; massive cancellation among ~√x terms
produces only O(log x) bits of net information.

## The Critical Insight

The error term E(x) = pi(x) - Li(x) has only **O(log x) bits** of information content,
which is well within polylog(x) bounds. The barrier to polylog computation is:

**COMPUTATIONAL, NOT INFORMATION-THEORETIC.**

We need to compute O(log x) bits of a function that's defined as a sum of ~x^{1/2} 
oscillating terms (zeta zero contributions). Each term has O(log x) bits, but their
sum (after massive cancellation) has only O(log x) bits.

Analogy: computing Σ_{k=1}^N cos(kθ) requires O(N) naively but has a closed form 
computable in O(1). The zero sum Σ_ρ R(x^ρ) is similar — a sum of many oscillating 
terms with a small net result — but NO known closed form exists.

The problem is structurally similar to a one-way function: the output has O(log x) bits
but seems to require O(x^{1/2}) computation to produce. Whether this hardness is 
fundamental or merely a gap in our knowledge is the key open question.

## Verdict
- Information-theoretically: **polylog is possible** (only O(log x) bits needed)
- Computationally: **no known method** extracts those bits in less than O(x^{1/2}) time
- The breakthrough (if it exists) must find a way to compute the NET effect of ~x^{1/2}
  zeta zero contributions without computing them individually

## Implication for Research Direction
The most promising direction is finding algebraic/analytic structure in the zero sum
that enables bulk evaluation — similar to how geometric series collapse sums of 
exponentials into closed forms.
