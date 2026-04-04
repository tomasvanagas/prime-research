# Focus Queue: Deep Dive Research Tasks

Each focused session works on exactly ONE task from this list.
The task number rotates automatically across sessions.

---

## Task 1: Kt Complexity of delta(n)

**Goal:** Empirically estimate the time-bounded Kolmogorov complexity Kt(delta(n) | n)
for small n values. Is there hidden structure?

**Experiments to run:**
1. Compute delta(n) = p(n) - round(R^{-1}(n)) for n = 1..100000
2. For each delta(n), find the shortest program that outputs it given n,
   using brute-force search over small programs (BFS over instruction sequences)
3. Plot Kt(delta(n)) vs log(n) -- is it O(polylog)? O(sqrt(n))? O(n^{1/3})?
4. Check if delta(n) values cluster or have patterns when viewed as binary strings
5. Test compressibility: gzip/bz2/lzma compression ratio of delta sequence vs random
6. Autocorrelation analysis of the delta sequence at multiple lags

**Key question:** If Kt(delta(n)|n) grows slower than expected, it suggests
small circuits exist, which would be strong evidence for a fast algorithm.

**Save to:** experiments/information_theory/kt_complexity/

---

## Task 2: Zeta Zero Structural Patterns

**Goal:** Search for algebraic or arithmetic relations among Riemann zeta zeros
that could enable fast summation of the explicit formula.

**Experiments to run:**
1. Load the 1000 zeros from data/ directory
2. Test all pairwise ratios gamma_i/gamma_j -- are any close to simple rationals?
3. Use PSLQ (integer relation algorithm) to search for linear relations
   among small subsets of zeros with algebraic coefficients
4. Compute the discrete Fourier transform of {gamma_1, ..., gamma_N} --
   is there spectral structure beyond GUE?
5. Test if partial sums S_K = sum_{k=1}^{K} R(x^{rho_k}) satisfy any recurrence
6. Search for patterns in gamma_n mod 1, gamma_n mod pi, gamma_n mod log(2*pi)
7. Test if the zero sequence can be modeled as eigenvalues of a SPARSE matrix

**Key question:** Any exploitable global structure in the zeros would collapse
the sqrt(x) barrier for the explicit formula.

**Save to:** experiments/analytic/zeta_structure/

---

## Task 3: Novel Identity Search

**Goal:** Use computational algebra to search for identities that relate
the oscillatory part of pi(x) to computable functions.

**Experiments to run:**
1. Define f(x) = pi(x) - R(x) (the oscillatory residual)
2. Compute f(x) for x = 2..100000 using sieve + mpmath R(x)
3. Use PSLQ to search for relations: does f(x) satisfy any identity
   involving elementary functions of x? (log, sqrt, trig, etc.)
4. Test Wilf-Zeilberger method: can f(x) be expressed as a definite sum
   with a closed-form certificate?
5. Search for algebraic relations between f(x) and:
   - Bernoulli numbers B_k
   - Values of zeta at integers zeta(2), zeta(3), ...
   - Dirichlet L-function values
   - Modular form coefficients
6. Lattice reduction (LLL) on vectors of f(x) values to find minimal polynomials
7. Test if f(x) satisfies a DIFFERENTIAL equation in some regularized sense

**Key question:** A computable identity for the oscillatory part would directly
give an O(polylog) algorithm.

**Save to:** experiments/algebraic/identity_search/

---

## Task 4: Conditional Algorithms

**Goal:** Assuming GRH or other standard conjectures, what is the best possible
exact algorithm for p(n)?

**Experiments to run:**
1. Under GRH: Miller's test is O(log^2(n)) per number. What is the best
   algorithm for pi(x) that exploits this? Can batch Miller testing help?
2. Under GRH: the explicit formula error with T zeros is O(x*log(x)/T).
   What is the optimal T for error < 1? Is it sublinear in x?
3. Under Elliott-Halberstam: Goldston-Pintz-Yildirim gives bounded gaps.
   Does this help with counting?
4. Under strong forms of the twin prime conjecture: can gap structure
   accelerate counting?
5. Under the Riemann Hypothesis: implement Schoenfeld's explicit bounds
   and test if they enable a faster binary search for p(n)
6. Under Cramer's conjecture (gaps < C*log^2(p)): what is the best
   search algorithm for p(n) given R^{-1}(n)?
7. Implement and benchmark the best CONDITIONAL algorithm

**Key question:** Even a conditional O(x^{1/3}) or O(x^{1/4}) algorithm
would be a significant theoretical contribution.

**Save to:** experiments/analytic/conditional/
