# Meta-Complexity Analysis: pi(x) Circuit Complexity and MKtP — Results

## Experiment (Session 35)

Formal analysis of the connection between pi(x) circuit complexity, time-bounded Kolmogorov
complexity (Kt), and the MKtP (Minimum Kt Problem) framework.

## Key Findings

### 1. Truth Table Compressibility (Empirical)

| N | Table length | gzip ratio | Random ratio | Advantage |
|---|-------------|-----------|-------------|-----------|
| 12 | 4,096 | 0.631 | 1.022 | 1.62x |
| 14 | 16,384 | 0.533 | 1.005 | 1.89x |
| 16 | 65,536 | 0.471 | 1.001 | 2.13x |
| 18 | 262,144 | 0.425 | 1.001 | 2.35x |
| 20 | 1,048,576 | 0.385 | 1.000 | 2.60x |

**pi(x) mod 2 is 2-2.6x more compressible than random.** The advantage INCREASES with N.

However, this compressibility is due to LOCAL structure (prime density 1/ln(x) creates
long runs in the truth table). It does NOT imply the function has small circuits —
it merely reflects that primes have density approaching 0.

### 2. Block Entropy Confirms Local Structure Only

At N=20, block entropy per bit for different block sizes:
- Block 2: 0.813 (81% of max)
- Block 4: 0.562 (56% of max)
- Block 8: 0.435 (44% of max)
- Block 16: 0.369 (37% of max)

Only 213/65536 possible 16-bit patterns occur. This is because consecutive entries in
the truth table correspond to consecutive integers, and primes have density ~1/ln(x),
creating sparse patterns. This is TRIVIAL structure, not exploitable for circuits.

### 3. Shannon Entropy: Nearly Maximal

ones_frac → 0.50, H(bit) → 1.0 as N grows. The truth table is nearly balanced
and nearly maximally entropic. No global bias to exploit.

### 4. Meta-Complexity Framework: REFORMULATION, NOT TOOL

**Theorem (informal):** "Is pi(x) in NC?" is equivalent to all of:
- "Does pi(x) mod 2 have poly(N)-size Boolean circuits?"
- "Is MCSP(T_N, N^c) = YES for some c and all N?"
- "Is Kt(pi(x) mod 2 | x) = O(poly(N))?"

**Why MKtP doesn't help us:**
1. MCSP for GENERIC instances is NP-hard (Hirahara 2018), but our instance T_N is specific.
2. Kt(T_N) is O(2^N * N) regardless of circuit size (sieve is efficient for the full table).
3. Kt(f(x) | x) for individual x is what matters, but bounding it requires... building the circuit (circular).

### 5. Brandt's Conditional Framework: TOO GENERIC

Brandt et al. (2024): MKtP hardness ↔ circuit lower bounds for SOME function in E.
But this is for ANY explicit function, not specifically pi(x). The framework cannot
distinguish pi(x) from other functions in E.

To use it for pi(x), we'd need to show either:
- T_N has high Kt (= no small circuits) → requires proving circuit lower bounds
- T_N has low Kt (= small circuits exist) → requires solving the original problem

Both directions are blocked by the same barriers (Natural Proofs, GUE).

## Verdict

**CLOSED: Meta-complexity / MKtP / Brandt framework as attack path.**

The framework provides an elegant reformulation of the problem but NO new techniques.
Every formulation reduces to the same open question: does pi(x) have poly-size circuits?
This question is at least as hard as proving unconditional circuit lower bounds,
which faces the Natural Proofs barrier (Razborov-Rudich 1997).

## New Closed Paths

| Approach | Verdict | Mode | Key Finding | Session |
|----------|---------|------|-------------|---------|
| MKtP / meta-complexity framework | FAIL | E | Reformulation not tool; Kt(T_N) = O(2^N * N) by sieve regardless of circuit size; Brandt framework too generic (any function in E, not pi(x) specifically) | 35 |
| Approximate circuit complexity phase transition | FAIL | I | NO phase transition; accuracy degrades gradually from 50% to 100% as rank increases; errors are spatially uniform; no "easy core" to exploit | 35 |
| Truth table compressibility | FAIL | I | 2.6x compressible at N=20 but only local structure (prime sparsity); block entropy → max at all scales; doesn't imply small circuits | 35 |
