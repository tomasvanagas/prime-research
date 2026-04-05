# Automata-Theoretic Prime Counting -- Results

## Experiment Date: 2026-04-05

## Idea
Build product DFA for "n has no factor <= B" by composing individual divisibility DFAs for each prime p <= B. Count accepted binary strings up to x via digit DP over DFA states. This gives O(M * log x) time where M = primorial(B).

## Key Results

### 1. Correctness: VERIFIED
Digit DP over product DFA correctly counts B-rough numbers for all tested (B, x) pairs. Matched brute-force counts exactly.

| B | primorial(B) | x=100 | x=1000 | x=10000 |
|---|---|---|---|---|
| 7 | 210 | 21 | 227 | 2284 |
| 11 | 2,310 | 20 | 206 | 2076 |
| 13 | 30,030 | 19 | 189 | 1916 |

### 2. The Gap: B-rough count vs prime count
At x=100,000, pi(x)=9,592:

| B | primorial | rough_count | gap | gap/pi(x) | extra_composites |
|---|---|---|---|---|---|
| 2 | 2 | 49,999 | 40,407 | 4.21 | 40,408 |
| 5 | 30 | 26,665 | 17,073 | 1.78 | 17,076 |
| 7 | 210 | 22,856 | 13,264 | 1.38 | 13,268 |
| 11 | 2,310 | 20,778 | 11,186 | 1.17 | 11,191 |
| 13 | 30,030 | 19,180 | 9,588 | 1.00 | 9,594 |

Even at B=13, the number of extra composites (9,594) nearly equals the number of primes (9,592). The gap decreases slowly with B.

### 3. Primorial Growth (The Barrier)

| B | primorial(B) | log2 |
|---|---|---|
| 7 | 210 | 7.7 |
| 13 | 30,030 | 14.9 |
| 19 | 9,699,690 | 23.2 |
| 29 | 6,469,693,230 | 32.6 |
| 31 | 200,560,490,130 | 37.6 |

To eliminate the gap entirely: need B >= sqrt(x). For x=10^100, need B >= 10^50, giving primorial(10^50) ~ e^{10^50} states. This is super-exponential, far worse than O(x).

### 4. Timing (Digit DP)

| B | M | x=10^6 | x=10^7 | x=10^8 |
|---|---|---|---|---|
| 7 | 210 | 0.006s | 0.009s | 0.013s |
| 11 | 2,310 | 0.039s | 0.076s | 0.114s |
| 13 | 30,030 | 0.220s | 0.533s | 0.922s |

Timing scales linearly with M and logarithmically with x, confirming O(M * log x). But M itself is the problem.

### 5. Matrix Structure
- Transition matrices have exactly M nonzero entries (one per column)
- They are NOT permutation matrices (the map s -> 2s mod M is not bijective when gcd(2,M)>1, but M is always odd since primorial(B>=2) includes factor 2... actually M includes 2 so 2s mod M maps both s and s+M/2 to same target)
- Sparsity is 1/M, extreme sparsity but doesn't help since we need all M states

### 6. DFA Minimization
The product DFA is already minimal (or near-minimal). States are residues mod M = primorial(B). Distinct residues give distinct acceptance behavior since accept = "nonzero mod every p <= B". Therefore all M states are distinguishable and necessary.

## Fundamental Barrier
1. **To count primes exactly**: need B >= sqrt(x), giving M = primorial(sqrt(x)) ~ e^{sqrt(x)} states
2. **The DFA is already minimal**: no compression below primorial(B) states
3. **Sub-sqrt sieving leaves a massive gap**: extra composites with all factors > B cannot be identified by the DFA
4. **Even with digit DP** (O(M log x) instead of matrix exponentiation O(M^omega log x)), M = e^{sqrt(x)} is exponential

## Verdict: CLOSED

The automata-theoretic approach is a reformulation of the classical sieve with worse complexity. The digit DP method is an interesting algorithmic technique (O(M * log x) for counting B-rough numbers), but the required state space M = primorial(sqrt(x)) ~ e^{sqrt(x)} makes it exponential. The DFA is provably minimal, so no compression is possible.

This confirms the existing CLOSED_PATHS.md entries:
- "DFA product automaton sieve: Minimized states = primorial(p_k) ~ e^{sqrt(x)}, worse than O(x)"
- "Automata-theoretic sieve (DFA product): Minimized DFA has exactly primorial(y) states; exponential"

No new insight toward O(polylog) prime computation.
