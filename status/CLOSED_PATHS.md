# Closed Paths: Master Lookup (655+ Approaches)

**SEARCH THIS FILE before proposing any approach.** Use grep/ctrl-F.

Last updated: 2026-04-05 (Sessions 1-36, 197+ sub-agents)

## Failure Modes
- **C** = Circularity (needs primes to compute primes)
- **E** = Equivalence (reduces to zeta zero sum / explicit formula)
- **I** = Information Loss (smooth approximation loses critical bits)

---

## Analytic / Explicit Formula / Zeta Zeros

| Approach | Verdict | Mode | Key Finding | Session |
|----------|---------|------|-------------|---------|
| Sparse matrix model for zeta zeros | FAIL | I | Tridiagonal Jacobi fits zeros with <0.05% error but uses 2N-1 params for N zeros (no compression). Off-diagonal entries b_i have 35% variation, no simple pattern. Toeplitz (N params) degrades at N=100. True extrapolation: 0.8-2.6% error. Information content not reduced. | 25 |
| Explicit formula partial sums recurrence | FAIL | I | No linear (r<=20), nonlinear (deg<=3), or difference recurrence. Residuals 2-6%. GUE phases make each term unpredictable from predecessors | 25 |
| Zeta zero pairwise rational relations | FAIL | I | 499,500 ratios tested; closest to p/q (q<=100) matches are no better than random (KS p=1.3e-6 but wrong direction — zeros FARTHER from rationals). GUE repulsion effect. | 25 |
| PSLQ linear relations among zeros | FAIL | I | 13,000+ tests at 60-digit precision: zero relations among subsets of 3-5 zeros with {1,pi,log(2pi)}. 1,225 pairwise tests negative. Large-K hits are lattice artifacts (random baseline confirms). Zeros linearly independent over Q. | 25 |
| DFT spectral structure of zeros | FAIL | I | Power spectrum matches GUE (corr=0.9999). Spectral flatness 0.93-0.999 at high freq (white noise). Pair correlation matches GUE 1-(sin(pi*r)/(pi*r))^2. Only faint p=2 signal (12x median). Number variance logarithmic (GUE). O(N) bits incompressible. | 25 |
| Zeta zeros mod constants (equidistribution) | FAIL | I | gamma_n mod m uniform for all 10 moduli tested (1, pi, log(2pi), 2pi, e, log(2,3,5,7)). KS p-values all >0.4. Weyl sums at 1/sqrt(N). Discrepancy BELOW random (GUE repulsion). Joint (mod 1, mod pi) independent. No arithmetic structure. | 25 |
| Convergence acceleration of zero sum | FAIL | I | Tested Richardson, Aitken Δ², Shanks, Padé, Cesàro, Euler-Maclaurin on partial zero sums for x=10^4..10^7. Errors GROW as N^{+1.0} (random walk). Best method (Shanks) gives O(1) improvement only. GUE-random phases make each term independent — no structured error to accelerate. | 32 |
| Self-correcting explicit formula + integer constraints | FAIL | I | Use fewer zeros + rounding/monotonicity/primality constraints. Truncation error is systematic bias ~O(N), not random noise. Constraints only help when error < 1, but error >> 1 with few zeros. Autocorr > 0.98. No triangulation possible. | 35,36 |
| Explicit formula proper convergence (mpmath R(x^ρ)) | FAIL | I | Naive R(x^ρ) summation DIVERGES: error grows from 3.5 (K=0) to 2076 (K=100) at x=10^4. Complex li branch cuts cause numerical instability. Only Lagarias-Odlyzko contour method is numerically stable. | 36 |
| Explicit formula + few zeros | FAIL | E | K^{-0.01} convergence | 5 |
| Explicit formula + 50 zeros | FAIL | I | More zeros can make WORSE | 6 |
| Explicit formula + 30 zeros | PARTIAL | E | 10% exact, best analytic | 10 |
| Weil explicit + Gaussian kernel | FAIL | I | Uncertainty principle: resolution * bandwidth >= const | 4 |
| Weil explicit + Beurling-Selberg | FAIL | I | Same uncertainty barrier | 4 |
| Mobius inversion of Pi(x) | PARTIAL | E | 10-30x better but still O(sqrt(x)/T) | 4 |
| Smooth zero sum integral | PARTIAL | E | Better than discrete but O(1) errors | 4 |
| Heuristic candidate generation + AKS | FAIL | I | Generate candidates near R^{-1}(n), filter with mod30/zeta/sieve. Candidates scale as n^0.577 (=sqrt(n)), not polylog. Interval width O(sqrt(p(n))*log(p(n))) is irreducible without computing zeta zero sum. With C=3.0 multiplier, 100% containment for n<=1000. Sieve reduces to ~4% of interval but still O(sqrt(n)) candidates. | 25 |
| Smoothed explicit + deconvolution | FAIL | I | Uncertainty principle blocks | 4 |
| Perron shifted contour (sigma>1) | FAIL | E | Numerically unstable, 20-40% error | 6 |
| Perron contour integral | FAIL | E | PROVEN equivalent to explicit formula | 5 |
| Cipolla Pade resummation | FAIL | E | Divergent, Pade can't generalize | 5 |
| Richardson/Aitken extrapolation | FAIL | E | No acceleration | 5 |
| Carlson's theorem | FAIL | E | Analytic continuation = explicit formula | 6 |
| Ramanujan-type series / Pade | FAIL | E | Cipolla DIVERGES, Pade gives 1372 mean error | 6 |
| Resurgent / trans-series | FAIL | E | Non-perturbative sectors ARE zeta zeros | 6 |
| Borel summation | FAIL | E | NOT APPLICABLE (power-law not factorial) | 6 |
| Riemann-Siegel analog | FAIL | E | = Lagarias-Odlyzko method | 6 |
| Zeta zero compression/FMM | FAIL | I | Spectral flatness 0.91, phase catastrophe | 9 |
| Connes-Weil quadratic form | FAIL | E | Finds zeros but can't sum them fast | 9 |
| Resummation (Borel/Pade/Shanks/Richardson/Mellin-Barnes) | FAIL | E | GUE statistics prevent partial sum prediction | 10 |
| FFT acceleration of explicit formula | FAIL | E | Constant factor only, not asymptotic | 7 |
| Gram points + Riemann-Siegel | PARTIAL | E | O(x^{1/2}*polylog) -- already known best | 8 |
| Analytic continuation (7 methods) | FAIL | E | Natural boundary at Re=0 | 8 |
| LLL/PSLQ integer relations | FAIL | I | Coefficients grow as O(n^{1/3}) | 6 |
| Selberg deconvolution | FAIL | E | 7-64% error, = Legendre sieve | 9 |
| Higher-order Richardson (orders 2-10) | FAIL | E | Error oscillatory not smooth; condition number explodes; diverges at high orders | 11 |
| Levin u/t-transform on zero sum | FAIL | E | WORSE than plain for large x; error not alternating | 11 |
| Weniger delta transform on zero sum | FAIL | E | Same failure; oscillatory error defeats all transforms | 11 |
| Smoothed explicit formula (psi_k/ln) | FAIL | I | psi(x)/ln(x) loses prime power corrections; error ~10% | 11 |
| Zero-count scaling analysis | CLOSED | E | N_min scales as x^{0.25-0.37} (power law), not polylog | 11 |
| Zeta oracle query complexity | CLOSED | E+I | M(x)=Theta(sqrt(x)*log(x)) evals needed; cond(sigma=1/2)=sqrt(x); N(T) is O(1) but positions cost Omega(N(T)); oracle cannot bypass sqrt(x) barrier | 21 |

## Algebraic / Number Theory / p-adic

| Approach | Verdict | Mode | Key Finding | Session |
|----------|---------|------|-------------|---------|
| CRT/modular construction | FAIL | C | Fundamentally circular | 3,7,9 |
| Modular CRT for pi(x) | FAIL | E | pi(x) mod m costs same as pi(x) | 4,5,7 |
| Chebotarev density theorem | FAIL | C | Need p to compute Frob_p | 4 |
| Artin L-functions | FAIL | E | = Riemann explicit formula | 4 |
| Elliptic curve point counts | FAIL | C | Can't invert a_p to p | 4 |
| p-adic counting | FAIL | E | Encode class numbers, not pi(x) | 4 |
| Modular forms approach | FAIL | C | Euler products need primes | 5 |
| p-adic interpolation | FAIL | I | Mahler coefficients grow super-exponentially | 8 |
| p-adic analysis (5 experiments) | FAIL | I | Not q-adically continuous, equidistributed | 10 |
| Modular forms (7 experiments) | FAIL | C | Encode primes individually not collectively | 10 |
| Perfectoid tilting | FAIL | E | Tilting preserves complexity | 9 |
| Prismatic cohomology | FAIL | E | Formula exists but >= O(x^{1/2+eps}) | 9 |
| Condensed mathematics | FAIL | E | Categorical, not computational | 9 |
| Geometric Langlands | FAIL | E | Same complexity as explicit formula | 9 |
| Motivic integration | FAIL | E | Structural obstruction | 9 |
| Number field algebraic decomposition | FAIL | E | Splitting doesn't reduce total cost | 7 |
| Class field tower prime reconstruction | FAIL | C+E | Splitting in m quadratic fields gives ~1 bit each (Chebotarev ~50/50), O(log n) fields suffice for CRT. But pi_split(x,d) = (pi(x) + character sums)/2, so computing it requires pi(x) (C) or zeros of L(s,chi_d) which are ADDITIONAL to zeta zeros (E). Net cost O(x^{1/2+eps}*log n), strictly worse. See experiments/proposals/class_field_splitting.py | 25 |
| PSLQ/LLL exhaustive identity search for f(x)=pi(x)-R(x) | FAIL | I | All 6 relation types (linear, polynomial deg 2-4, recurrence, modular, functional, discrete derivative) tested x=2..10000; all single-point PSLQ relations spurious (cross-validation residuals ~10^4); no recurrence (rms~0.31); Fourier confirms zeta-zero dominance | 18 |
| Iwasawa theory | FAIL | I | lambda,mu invariants give no shortcut | 7 |
| Cloitre analytic recurrence (2025) | EXISTS | E | Exact but uses zeta evaluation, not competitive | lit |
| Farey fractions | FAIL | C | Require phi(k) -> factoring | 8 |
| Arithmetic derivative | FAIL | C | Requires factoring, O(n^{1/4}) | 8 |
| NFS-type L[1/3] for pi(x) (4 sub-approaches) | FAIL | E+C | Norm sieve = residue class counting (E); Chebotarev = circular + O(x^{1/2}) error (C+I); class groups = Euler products (E); Artin L-fns = explicit formula (E). NFS exploits multiplicative structure of single N; pi(x) is additive-global with no analog. See experiments/algebraic/lthird_analysis.md | 19 |
| LLL lattice reduction for algebraic relations in f(x)=pi(x)-R(x) | FAIL | I | Minimal polynomial search deg 2-8 at x=100,1000,10000,100000: coefficient norms follow Dirichlet bounds for generic transcendentals (no anomalous structure). Multi-point polynomial search (k=2,3, deg<=3): all residuals at float64 machine epsilon, no genuine relations. Algebraic independence test on 6-point vector: min ||c||=201, consistent with random reals. Scaled g(x)=f(x)/sqrt(log(x)) search: no integer relations. Polynomial-in-log(x) models: RMSE~1.35 (explains <10% variance). See experiments/algebraic/identity_search/lll_results.md | 29 |

## Sieve / Combinatorial / Counting

| Approach | Verdict | Mode | Key Finding | Session |
|----------|---------|------|-------------|---------|
| Sieve weights formula | FAIL | I | Parity barrier (Selberg) | 4 |
| Batch Möbius sieve (recursive halving phi) | FAIL | E | Split prime set in half at each recursion level. Tree depth O(log(pi(√x))) but branching 2^{half_size}. Total leaves = 2^{pi(√x)} same as unpruned standard tree. At x=10^4, 30x slower than Lucy_Hedgehog despite pruning 64% of terms. Branching at x=10^6 would be 2^84 ~ 10^25. Mathematically identical to standard IE, just reordered. | 36 |
| Automata digit DP prime counting | FAIL | I | Product DFA for B-rough numbers: O(M·log x) via digit DP where M=primorial(B). Timing: M=210 at x=10^6 in 6ms, M=30030 at x=10^8 in 0.9s. DFA provably minimal (all M states distinguishable). Gap: at B=13, false positives ≈ true primes. Need B≥√x → M=e^{√x} states. | 36 |
| Randomized inclusion-exclusion for pi(x) | FAIL | I | Random subset sampling: std ~2^k/√(samples). For x=500 (k=8), need ~370K samples vs 256 full IE terms. Variance reduction impossible — randomization makes IE strictly worse. | 36 |
| Convolution sieve via polynomial NTT | FAIL | E | Π_{p≤B}(1-z^p) mod z^{x+1}: each factor sparse (2 terms) but π(B) factors give O(π(B)·x log x) total ~ O(x^{3/2}/log x), worse than Meissel-Lehmer. Evaluation at roots of unity gives Ramanujan sums, no new info. | 36 |
| Redheffer/Mertens matrix | FAIL | E | O(x^{2/3}) same as pi(x) | 8 |
| Wheel factorization | PARTIAL | - | Constant factor only, NOT complexity class | 8 |
| Hierarchical decomposition | FAIL | E | Legendre tree depth = pi(sqrt(x)) | 7 |
| Buchstab tree pruning/memo | FAIL | E | Memoized Buchstab IS Lucy DP; O(x^{3/4}) distinct args, O(x^{2/3}) with DR | 11 |
| Vaughan identity exact computation | FAIL | E | Type I/II/bilinear all O(x^{2/3}) individually; designed for bounds not computation | 11 |
| Legendre sieve FFT/NTT | FAIL | E | floor(x/d) breaks multiplicativity; cannot use convolution acceleration | 11 |
| Alternative DP formulations (4 variants) | FAIL | E | Largest/smallest prime factor DP, P_k DP, Mobius DP all reduce to O(x^{2/3}) | 11 |
| Dirichlet series extrapolation to s=0 | FAIL | E | Pole at s=1 makes extrapolation unstable; IS Lagarias-Odlyzko method | 11 |
| Double hyperbola decomposition | FAIL | E | No way to split pi(x) into two independently O(x^{1/2})-computable pieces | 11 |
| Mertens function shortcuts | FAIL | E | O(x^{2/3}) -- same as pi(x); H-T gives O(x^{3/5}) for M(x) but doesn't transfer to pi(x) | 6,11 |
| M(x) hyperbola factorization | FAIL | E | No mu=f*g with both partial sums easy; self-referential structure fundamental | 11 |
| M(x)->pi(x) transfer | FAIL | E | H-T O(x^{3/5}) for M(x) exploits signed cancellation; pi(x) is positive count, no analog | 11 |
| Dirichlet series / prime zeta | FAIL | E | Perron integral = explicit formula | 6 |
| Prime zeta P(s) Mobius extraction | FAIL | E | P(s) poles at rho, rho/2, rho/3... MORE terms than explicit formula; Vandermonde recovery cond~10^26 | 13 |
| CRT reconstruction of pi(x) | FAIL | E | pi(x) mod q costs same O(x^{2/3}) as pi(x); CRT multiplies cost by k moduli; entropy ratio ~1.0 | 13 |
| Prime race E(x;q) for CRT | FAIL | E | E(x;4) spectral flatness/power-law same as pi(x)-Li(x); 40 L-function zeros give 0% exactness; K_exact~sqrt(x); p(n)mod q has near-max entropy; L-zeros INDEPENDENT of zeta zeros so CRT costs MORE | 20 |
| Recursive pi(x/d) identity search | FAIL | E | Best: pi(x)=pi(x/3)+pi(x/5)+2*pi(x/7)+g(x), but g(x) grows as x^{0.22}; correction encodes primes in interval | 13 |
| Smooth number elimination | FAIL | I | 12% reduction, insufficient | 5 |
| Bit-by-bit construction | WORKS | E | log2(p(n)) pi(x) calls, each O(x^{2/3}) | 5 |
| Recursive pi(x) via pi(x/2) | FAIL | I | Error O(sqrt(x)/ln^2(x)), too large | 7 |
| Recursive identity (8 methods) | FAIL | E | All reduce to Meissel-Lehmer or explicit formula | 9 |
| GCD/coprimality matrix | FAIL | I | Eigenvalues != primes | 8 |
| Newton's identities (power sums) | WORKS | E | O(x^{2/3}) per sum, unstable deg>50 | 7 |
| Nonlinear sieve: floor products | FAIL | E | delta(a,b,x) doesn't separate primes; polynomial fits overfit catastrophically; K^d >= sqrt(x) required | 14 |
| Nonlinear sieve: comparisons/thresholds | FAIL | E | Parity of floor(x/k) uncorrelated with primality; gap function misses 96%+ of primes; mod/XOR no signal | 14 |
| Nonlinear sieve: multiplicative structure | FAIL | E | GCD, collisions, quadratic-sieve analog all fail; M(x/k) products: max_err=44 on test set | 14 |
| Nonlinear sieve: bitwise/TC^0 | FAIL | E | Bit features + AND: 85% exact but no generalization; degree d in K values needs K >= x^{1/(2d)} | 14 |
| Nonlinear sieve: parity barrier bypass | FAIL | E | Nonlinear ops CAN distinguish primes from semiprimes but cost O(x^{3/2}) -- worse than Meissel-Lehmer | 14 |
| Nonlinear sieve: floor identities | FAIL | E | No exact polynomial identity; recursive pi(x/k) formulas have O(1) residual; overfitting scales with degree | 14 |

## Machine Learning / Statistical / Fitting

| Approach | Verdict | Mode | Key Finding | Session |
|----------|---------|------|-------------|---------|
| ML correction fitting | FAIL | I | 0.6% exact, no generalization | 3 |
| Oscillatory formula search | FAIL | I | 1.1% exact | 3 |
| Simultaneous least-squares (8 features) | FAIL | I | 1.7% train / 0.3% test | 6 |
| Correction lattice + zeta oscillatory | FAIL | I | 10.5% train / 5.0% test | 6 |
| Optimal 10-param formula | FAIL | I | 1.5% train / 0.8% test (overfits) | 6 |
| AR(1)-AR(10) correction | FAIL | I | 4-5% exact, 5.04 bits/prime entropy floor | 6 |
| Deep ML (ridge/AR) | FAIL | I | 5.4% test exact -- fundamentally limited | 8 |
| Polynomial in ln(n) | FAIL | I | 1.2% exact, 0.3% generalization | 3 |
| Fourier correction | FAIL | I | 2.9% overfit, 0% generalization | 3 |
| Symbolic regression/PSLQ | FAIL | I | No pattern in delta(n) | 5 |
| Polynomial hash (deg 2-50) | FAIL | I | 0% test accuracy at ALL degrees | 6 |
| Decision tree (unlimited depth) | FAIL | I | 100% train, 19.1% test -- memorizes | 6 |
| Transformer neural | FAIL | I | 1.1% -- best ML result, still useless | 10 |
| FNO (Fourier Neural Operator) | FAIL | I | 0.5% | 10 |
| Gaussian Process | FAIL | I | 0.45% | 10 |
| MLP | FAIL | I | 0.4% | 10 |
| KAN (Kolmogorov-Arnold Network) | FAIL | I | 0.35% | 10 |
| Botkin partition model | FAIL | I | 10/10000 exact, model not formula | 8 |

## Quantum Computing / Quantum Information

| Approach | Verdict | Mode | Key Finding | Session |
|----------|---------|------|-------------|---------|
| Grover on Deleglise-Rivat sieve | FAIL | E | O(x^{1/3}) quantum, still 10^34 ops for 10^100 | 7 |
| Shor-like (periodicity search) | FAIL | I | Primes have no periodicity (Mauduit-Rivat) | 7 |
| Berry-Keating QPE | FAIL | E | Even if Hamiltonian exists, need 10^51 zeros | 7 |
| Primon gas / GUE simulation | FAIL | C | Circular or statistical-only | 7 |
| BQP membership of exact pi(x) | FAIL | - | Almost certainly NO (exact pi(x) is #P-hard) | 7 |
| Quantum channel capacity | FAIL | I | Holevo bound = classical, no advantage | 10 |
| Entanglement-assisted counting | FAIL | I | Superdense coding halves comms, not computation | 10 |
| Holographic principle for primes | FAIL | I | Volume-law O(N), not surface O(sqrt(N)) | 10 |
| Quantum error correction / stabilizer | FAIL | I | No low-dimensional error syndrome | 10 |
| MPS / tensor network | FAIL | I | Bond dimension ~ N^{0.49}, volume-law entanglement | 10 |
| Black hole information analogy | FAIL | I | No Page curve, no monotonic recovery | 10 |

## Dynamical Systems / Automata

| Approach | Verdict | Mode | Key Finding | Session |
|----------|---------|------|-------------|---------|
| Conway PRIMEGAME (FRACTRAN) | FAIL | E | O(p(n)^2) per prime, 37K steps to p=29 | 4,6 |
| FRACTRAN transfer operator | FAIL | E | O(2^p) per prime | 10 |
| Bost-Connes quantum mechanics | FAIL | E | Reduces to zeta(s) | 4 |
| Symbolic dynamics | FAIL | I | Near-random block complexity | 4 |
| All 256 elementary CA rules | FAIL | I | 0 generate primes; best 75.8% primality | 7 |
| Cellular automata (general) | FAIL | E | O(sqrt(N)) minimum | 10 |
| DFA product automaton sieve | FAIL | E | Minimized states = primorial(p_k) ~ e^{sqrt(x)}, worse than O(x) | 20 |
| Tensor sieve (MPS of divisibility DFAs) | FAIL | I | Transfer matrices full rank, volume-law SV entropy, no compression | 20 |
| Iterative maps (x+ln(x)+...) | FAIL | I | 4% exact at best | 6 |
| Rowland acceleration | FAIL | E | O(p^2) | 4 |
| Deterministic gap recurrence (7 variants) | FAIL | I | <25% accuracy, 5.04 bits/prime entropy | 7 |

## Topological / Exotic Mathematics

| Approach | Verdict | Mode | Key Finding | Session |
|----------|---------|------|-------------|---------|
| K-theory of Z | FAIL | I | K_0(Z)=Z, K_1(Z)=Z/2Z -- almost no prime info | 7 |
| Motivic cohomology | FAIL | C | H^1_mot(Spec(Z)) = Z^inf -- isomorphic to primes | 7 |
| Etale cohomology of Spec(Z) | FAIL | C | Recovers Euler product = knowing all primes | 7 |
| Topos theory / Spec(Z) sheaves | FAIL | E | Recovers explicit formula, no new content | 7 |
| Tropical geometry | FAIL | C | Wrong compression | 4 |
| Selberg trace / geodesics | FAIL | E | Isomorphic to explicit formula | 4 |
| TDA (topological data analysis) | FAIL | I | Detects structure but not computable | 10 |
| Sheaf cohomology | FAIL | E | Same barrier | 10 |
| Knot invariants | FAIL | C | No connection to prime counting | 10 |
| Non-standard analysis | FAIL | E | Transfer principle preserves complexity | 10 |
| Surreal numbers | FAIL | E | Same | 10 |
| IUT (Inter-Universal Teichmuller) | FAIL | E | No computational content | 10 |

## Information Theory / Complexity Theory

| Approach | Verdict | Mode | Key Finding | Session |
|----------|---------|------|-------------|---------|
| Levin universal search (AIT) | FAIL | I | K(p(n)|n) >= 0.5*log2(n) bits | 7 |
| Local-to-global reconstruction of pi(x) | FAIL | I | E(x)=pi(x)-Li(x) has only O(log x) bits — NO info-theoretic barrier to polylog. But O(log x) bits encoded across ~x^{1/2} GUE-random zero contributions via massive cancellation. Error autocorr=0.996 at lag 1 (helps adjacent queries only). Barrier is COMPUTATIONAL not informational. | 36 |
| Busy Beaver connection | FAIL | - | Wrong hierarchy (computability vs complexity) | 7 |
| Curry-Howard type theory | FAIL | - | Proof = algorithm, no gap to exploit | 7 |
| SNARG/STARK verification | PARTIAL | - | 50K ops verify but 10^68 to COMPUTE | 7 |
| Arthur-Merlin protocols | FAIL | - | No natural polylog hint exists | 7 |
| Reverse mathematics | FAIL | - | p(n) provable in RCA_0, barrier is computational | 7 |
| Hypercomputation / real oracles | FAIL | C | All constants tautological or need sqrt(p(n)) ops | 7 |
| CRT novel encoding (178 bits) | FAIL | - | Identifies info but no way to obtain it | 7 |
| PRG/Derandomization | FAIL | E | delta(n) appears easy for circuits | 10 |
| Randomness Extractors | FAIL | E | Seed-finding as hard as original problem | 10 |
| Crypto hardness reduction | FAIL | - | No reduction to/from known hard problems | 10 |
| One-way function connection | FAIL | - | Both directions equally hard | 10 |
| Expander graphs | FAIL | C | Need primes to build the graph | 10 |
| Pfaffian/FKT planar sieve encoding | FAIL | I | Three graph encodings tested (sieve bipartite, Lucy DP recurrence, prime dependency). ALL non-planar for x>=80. Treewidth grows as x^0.36 (sieve), x^0.51 (Lucy DP), x^0.43 (prime dep). Prime dep graph is complete (density=1.0). Unbounded treewidth rules out FPT/Pfaffian approach. | 26 |
| Derandomization (NW/IW) | FAIL | E | R(n) is deterministic, not random | 9 |
| TC^0 via Legendre sieve | FAIL | E | 2^{pi(sqrt(x))} terms, super-polynomial | 11 |
| TC^0 via Meissel-Lehmer | FAIL | E | Super-polynomial special leaves | 11 |
| TC^0 via Lucy DP | FAIL | E | Sequential depth O(sqrt(x)/ln(x)), not constant | 11 |
| TC^0 via Möbius function | FAIL | C | Requires factoring, not known in TC^0 | 11 |
| TC^0 via partial sieve (constant primes) | FAIL | I | Error Theta(x/ln(x)), same order as pi(x) | 11 |
| AKS in TC^0 via matrix powering | OPEN | - | Fixed-k MPOW IS in TC^0 (RAIRO 2000), but AKS needs k=polylog(n) GROWING; combination step needs depth O(log log n). IMP_k (k≥3) NOT in TC^0 unless TC^0=NC^1. Growing-dim MPOW is genuinely open at TC^0/NC^1 boundary | 11 |
| pi(x) mod 2 shortcut | FAIL | E | Parity as hard as full problem; pseudo-random, no pattern | 11 |
| Compressed pi(x) | PARTIAL | I | 5 bits/prime irreducible entropy | 8 |
| Lucy DP parallel depth | CLOSED | E | DAG depth = pi(sqrt(x)) exactly; fundamentally sequential for S(x). Telescoped (Meissel-Lehmer) gives depth pi(x^{1/3}) but both exponential in N=log(x) | 12 |
| Lucy DP floor-value mapping commutativity | FAIL | E | Mapping matrices NON-commutative (0 commuting pairs); cannot reorder sieve steps | 12 |
| pi(x) mod m shortcut (m>2) | FAIL | E | Conditional entropy H(Y|X) = 0.537 bits INVARIANT across all moduli m; no modular shortcut | 12 |
| Floor-value linear algebra | FAIL | E | Transformation is full-rank (~80% nonzero coeffs), massive cancellation (|coeffs| >> pi(x)); no compression | 12 |
| Meissel-Lehmer as NC circuit | FAIL | E | Depth O(log log x) = O(log N) but width O(x^{2/3}) = O(2^{2N/3}); exponential size in input | 12 |
| TG kernel for exact pi(x) | FAIL | I | arXiv:2506.22634 DEBUNKED: violates uncertainty principle, omits x^{1/2} factor, kernel too wide to resolve primes. AI-generated paper with errors. Authors (Kilictas/Alpay) have no number theory background; Alpay publicly admits uploading arxiv papers to influence LLM training data (Medium, July 2025). Zero citations, no peer review. | 12,30 |
| Andrews-Wigderson for pi(x) | FAIL | E | FOCS 2024 constant-depth arithmetic circuits: wrong model (fields not rings), wrong bottleneck (GCD not matrix powering), wrong circuit model (arithmetic not Boolean) | 12 |
| Succinct/lattice point counting | FAIL | E | Barvinok requires fixed dim; encoding coprimality needs dim ~ sqrt(x)/log(x), exponential in N | 12 |
| Permanent/determinant encoding | FAIL | E | Would resolve NC question; no natural matching interpretation for primes; circular (matrix entries encode primes) | 12 |
| H-T cancellation transfer quantified | CONFIRMED | E | M(x) cancels 99.9+% of terms (|M(x)|~sqrt(x) vs 0.6x nonzero); pi(x) has ZERO cancellation (all +1). Definitive. | 12 |
| Wilson's theorem in TC^0 | FAIL | E | (n-1)! mod n requires 2^N multiplications; Shamir's GapL gives NC^2 not TC^0 | 13 |
| Sum-of-two-squares primality in TC^0 | FAIL | C | Counting reps needs factoring; Cornacchia has O(log n) sequential depth; only handles p≡1(4) | 13 |
| BPSW in TC^0 (computation) | WORKS | - | BPSW IS computable in TC^0: MR(2) = scalar pow, strong Lucas = 2x2 MPOW (Mereghetti-Palano), Jacobi = TC^0 via GCD. Correctness for all n conditional on BPSW conjecture (verified to 2^64) | 13 |
| QFT (Grantham) in TC^0 | WORKS | - | Quadratic Frobenius test IS in TC^0: operates in 2D algebra = 2x2 MPOW. Error < 1/7710 per param. Deterministic only with GRH. 4 QFT pseudoprimes below 50000 | 13 |
| Strong Lucas standalone | PARTIAL | - | 12 pseudoprimes below 100000 (all caught by second param set). No unconditional correctness proof | 13 |
| Miller-Rabin fixed bases in TC^0 | WORKS | - | MR({2,3,...,37}) deterministic for n<3.317×10^24 (Sorenson-Webster 2016). Each base = scalar pow → TC^0. Gives nonuniform TC^0 for any fixed input length | 13 |
| Divide-and-conquer pi(x) | FAIL | E | pi(x)=pi(x/2)+delta: error in each interval O(sqrt(x_level)), accumulates to O(sqrt(x)). Buchstab tree has O(x^{2/3+eps}) nodes. No recursion gives sublinear error accumulation | 15 |
| Randomized zeta zero sampling | FAIL | E | Must use 100% of zeros (variance analysis); no sampling gain possible | 15 |
| Probabilistic sieve (I-E sampling) | FAIL | E | Variance 10^6-10^9x WORSE than exhaustive due to massive cancellation in I-E | 15 |
| Hash-based prime counting sketch | FAIL | I | Exact sketch needs O(x^2/ln^2(x)) bits — worse than storing all primes | 15 |
| Quantum counting (structured) | OPEN | E | Black-box: O(sqrt(x)). Lower bound: Omega(sqrt(ln x)). Gap requires exploiting PRIMES structure = original problem | 15 |
| Arithmetic circuit complexity path | FAIL | E | VP=VNC^2 (VSBR): depth is FREE in algebraic model. But pi(x) is NOT a polynomial over fields. Tau conjecture orthogonal. Does not offer new viable path | 15 |
| Monotone complexity lower bounds | FAIL | - | [pi(x)>=k]=[x>=p(k)] has TRIVIAL O(N) monotone complexity. Individual bits of pi(x) NOT monotone. Approach doesn't capture hardness of pi(x) | 15 |
| #TC^0 counting for BPSW | OPEN | E | If BPSW∈TC^0, pi(x)∈NC iff #TC^0⊆NC. Fermat residue coupling (2^{n-1} mod n) prevents batch counting. Equivalent difficulty to original problem | 15 |
| TC^0 MAJORITY fan-in for pi(x) | FAIL | E | MAJORITY can aggregate x bits in O(1) depth, but generating x=2^N input bits costs 2^N*poly(N). Bottleneck is INPUT GENERATION not aggregation | 16 |
| Batch modular exp via CRT | FAIL | E | f(k)=2^{k-1} mod k: CRT decomposition requires factoring all k (O(x log log x)); zero autocorrelation; decoupling exponent from modulus loses primality testing | 16 |
| Algebraic batch via Carmichael lambda | FAIL | E | f(k)=2^{(k-1) mod lambda(k)} mod k; computing lambda(k) for all k requires complete factorization sieve O(x); no aggregation shortcut | 16 |
| Period/Fourier exploitation of f(k) | FAIL | E | g_p(k)=2^{k-1} mod p has period lcm(ord_p(2),p), but k has varying factorization; mutual info between BPSW(n), BPSW(n+k) near zero for all lags | 16 |
| Novel intermediate quantity families | FAIL | E | Systematic analysis of 8 families (residues, polynomial evals, matrix eigenvalues, topology, representation theory, entropy, recursive, physical). ALL route back to floor values or zeta zeros. Only #TC^0 counting not immediately closed | 15 |
| Determinantal complexity of pi(x) | OPEN | - | pi(x) = degree-N multilinear polynomial in bits. Found N×N det reps for N=2,3,4. dc(pi_N)=poly(N) iff pi(x)∈GapL. For N>=10, generic polynomials DON'T have N×N det reps — pi(x) would need special structure | 15 |
| Class numbers h(-d) as det entries | FAIL | E | h(-d) = (w*sqrt(d))/(2*pi)*L(1,chi_d) — equivalent to L-values; constants with no x-dependence; cannot be affine-linear matrix entries | 16 |
| L-function values L(1,chi) as det entries | FAIL | C+E | Partial L-functions need primes (circularity); full values are global aggregates; inverting = explicit formula (equivalence to zeta zeros) | 16 |
| Elliptic curve a_p as det entries | FAIL | C | a_p defined only at primes; a_n at composites = products of a_p (needs factoring); constants with no x-dependence | 16 |
| Regulators as det entries | FAIL | E | h*R = f(L(1,chi_d)) by class number formula — equivalent to L-values; transcendental, incompatible with GapL integer entries | 16 |
| Hybrid (class num + L-val + a_p + char) | FAIL | C+E | Characters periodic (residue class info, not primality); linear prediction POOR; products of characters = higher modulus characters; all roads back to zeta zeros/floor values | 16 |
| Binary carry structure of pi(x) | FAIL | I | Carry chains identical to generic counter (avg ~2.0). Communication matrix rank for bit_j(pi(x)) = 2^{N/2-1}+2 EXACTLY (verified N=8..16) = exponential. Spectral weight spread across all Fourier levels. Primes indistinguishable from random subsets. No shortcut from binary arithmetic. | 17 |
| NC^1 branching programs for primality | FAIL | E | OBP lengths for is_prime match random functions at all widths (2,3,5) and N (4,6,8,10). Communication matrix rank = 2^{N/2-1}+1 (half dimension + 1, explained by even/odd structure; odd columns have FULL rank). Fourier L1 norm ratio to random: 0.94->0.67 (decreasing, approaching random). Sensitivity = N (maximal). No NC^1-exploitable structure beyond TC^0 threshold gates. | 17 |
| IW97 BPP=P derandomization for pi(x) | FAIL | E | No BPP algorithm known for pi(x) in binary input model; IW97 presupposes efficient randomized algorithm. Cannot create efficiency from nothing. | 18 |
| NW/Nisan PRG-based prime reconstruction | FAIL | E | PRGs fool bounded computations but cannot compute functions outside the bounded class. If pi(x) not in the class, PRG useless. Prime indicator has PRF-like structure. | 18 |
| Approximation amplification R(x)->pi(x) | FAIL | I | Residual pi(x)-R(x) has communication rank 2^{N/2-1} (Session 17). Gap is sqrt(x) = 2^{N/2-1} independent bits. No derandomization/rounding/self-reduction can bridge this. | 18 |
| KI04 PIT connection via dc(pi_N) | FAIL | E | dc(pi_N) >= 2^{N/2-1}+2 (exponential) but pi_N is not an explicit VNP family. KI04 connects PIT to VP vs VNP, not to specific non-algebraic functions. PIT solves identity testing, not evaluation. | 18 |
| Reverse hardness amplification (random-looking residual -> PRG -> exact pi) | FAIL | I | Residual pi(x)-R(x) is DETERMINISTIC, not distributional. A PRG matching its distribution gives random samples, not exact values. | 18 |
| Batch Fermat residue derandomization for counting | FAIL | E | f(n)=2^{n-1} mod n has zero autocorrelation at all lags (Session 16). Cryptographically pseudorandom across consecutive n. Cannot be batch-computed or derandomized. | 18 |

## Encoding / Novel Representations

| Approach | Verdict | Mode | Key Finding | Session |
|----------|---------|------|-------------|---------|
| Super-convergent series | FAIL | I | 1.4% exact (R^{-1} ceiling) | 3 |
| Continued fraction of p(n)/n | FAIL | I | Higher CF terms unpredictable | 8 |
| Prime constant C_p | FAIL | C | Circular to compute | 8 |
| Integer rounding R(x) | PARTIAL | I | 26.7% exact for x in [1000,5000] | 8 |
| Prony method on P(s) | PARTIAL | E | 6/9 primes, Hankel ill-conditioned | 8 |
| Novel encoding (6 methods) | FAIL | E | All analytically equivalent | 8 |
| Self-referential formulas | FAIL | I | All smooth -> error O(sqrt(p)) | 8 |
| Prunescu-Shunia arithmetic term | EXISTS | - | Exact but 10^78913-digit intermediates | 5 |
| Willans/Wilson formula | EXISTS | - | Exact but (2^n)! needed | 1 |
| Mills' constant | FAIL | C | Circular -- constant encodes all primes | 1 |
| Euler product inversion | FAIL | C | Needs primes to compute | 5 |
| BBP digit extraction | FAIL | - | No BBP formula for prime constant | 4 |
| DPRM polynomial witness | FAIL | E | 10^102 bit witness for p(10^100) | 9 |
| JSWW polynomial inversion | FAIL | - | Exponential witness sizes | 7 |
| Generating function (Cauchy integral) | FAIL | C | Need primes to evaluate EGF | 7 |

## Additive Combinatorics / Goldbach / Partitions

| Approach | Verdict | Mode | Key Finding | Session |
|----------|---------|------|-------------|---------|
| Green-Tao AP interpolation | FAIL | C | Circular, tower-of-exponentials bounds | 4 |
| Circle method / Vinogradov | FAIL | I | Asymptotic only | 4 |
| Goldbach FFT deconvolution | FAIL | I | Phase loss destroys positional info | 6 |
| Goldbach phase retrieval | FAIL | I | 0% recall, ill-posed | 9 |
| Ono-Craig partition detection | FAIL | E | O(n^2) per test -- worse than Miller-Rabin | 4 |
| Partition-based (Ono-Craig) | FAIL | - | Detects primes but CANNOT be inverted | 6 |
| Additive structure (7 methods) | FAIL | I | No exploitable structure | 9 |
| Sumset P+P / representation r_2(n) | FAIL | E+C | r_2(n) defined in terms of primes; HL circle method error = zeta zeros; smooth sum gives 20% error | 16 |
| Ergodic theory / orbit complexity | FAIL | E+C | Block complexity maximal; entropy 5.04 bits/prime irreducible; transfer ops = zeta zeros; Gauss map unrelated | 16 |
| Model theory / o-minimality | FAIL | - | pi(x) not definable in o-minimal structures (step function); adding Z = undecidable; cell count O(N/lnN); orthogonal to computation | 16 |
| Tropical geometry for pi(x) | FAIL | I+E | Tropicalization loses all info (min=smallest prime); floor values ARE tropical objects; no new content | 16 |
| Sufficient statistics of floor values | FAIL | E | poly(logx)-bit statistic EXISTS but computing it from floor values IS the Meissel-Lehmer problem; no compression shortcut | 16 |
| Curve families over F_p for pi(x) | FAIL | E+C | a_p encodes L(E,s) not pi(x); AKS variety costs O(x*polylog); Frobenius eigenvalues = zeta zeros | 16 |
| S_n/GL_n representation theory | FAIL | C+E | Cycle structure gives sum 1/p ~ loglogn (Mertens, not pi); P(k) recovery cond~10^26; Fourier coeffs O(sqrt(pi(N))) random | 16 |
| Green-Tao nilsystem correlation | FAIL | I+E | 1-step nilseqs (Fourier) R^2=0.52, 2-step (bracket quadratics) R^2=0.37, combined CV R^2=0.25. Top freqs = 1/2,1/3,1/6 (sieve bias). p(n) prediction: 0.84 correct digits (CV) vs li^{-1} 1.51 digits. Nilseqs capture AP correlations not individual primes. Residual non-normal, spectral peaks remain. Reduces to Fourier fitting (Sessions 3,7,14). | 25 |

## Gap Prediction / Interpolation

| Approach | Verdict | Mode | Key Finding | Session |
|----------|---------|------|-------------|---------|
| Recursive gap prediction | FAIL | I | Drifts catastrophically | 3 |
| Gap history (1-5 previous) | FAIL | I | Best 71.4% (sparse artifact), 0.54 bits irreducible | 6 |
| Wheel transition matrix | FAIL | I | 34.6% prediction (vs 12.5% random) | 6 |
| Nearest-prime decoding | FAIL | I | 0.7% raw, 15.2% after bias correction | 6 |
| Ensemble voting (3 approx) | FAIL | I | 0% accuracy | 6 |
| Fractal/self-similar structure | FAIL | I | Zeta frequencies but NOT periodic | 7 |
| Scaling hypothesis | FAIL | I | FALSE: incommensurate frequencies never repeat | 7 |
| Interpolation (5 methods) | FAIL | I | Random walk delta(n) is irreducible | 9 |
| Multi-approximation cascade | FAIL | I | Errors too correlated (rho=0.31-0.90) | 5 |
| Self-correcting R^{-1} | PARTIAL | I | 75% correct small n, error grows | 5 |
| Recursive halving p(2n) from p(n) | FAIL | I | Error O(sqrt(n)*ln(n)) >> gap | 5 |

## Spectral Methods

| Approach | Verdict | Mode | Key Finding | Session |
|----------|---------|------|-------------|---------|
| Spectral/eigenvalue approach | FAIL | E | No viable formula found | 3 |
| Jacobi matrix approach | FAIL | E | Parameters not closed-form | 3 |
| Heat kernel / inverse Laplace | FAIL | E | Numerically unstable | 3 |
| Selberg trace formula | FAIL | E | Same barrier as explicit formula | 3 |
| Graph spectrum approach | FAIL | C | Requires knowing primes | 3 |
| MUSIC/ESPRIT spectral estimation | FAIL | E | Requires O(e^T) signal | 3 |
| Inverse Stieltjes | FAIL | E | Spectral resolution barrier | 8 |
| Compressed sensing on K zeros | PARTIAL | I | Lucky cancellations small x, doesn't scale | 7 |
| Matching pursuit (K zeros) | PARTIAL | I | K=2 error<0.5 at x=1000, fails large x | 7 |
| RMT/GUE approximation | FAIL | I | Correct statistics, wrong values | 8,10 |
| Cayley graph / Ihara zeta / spectral graph | FAIL | C+E | Cayley(Z/xZ, primes): circular (lambda_0=pi(x)); Ihara zeta = Selberg trace analog; GCD graph: eigenvalues = Ramanujan sums c_q(n), same ingredients as Meissel-Lehmer O(x^{2/3}); expander mixing error >> 1 | 13 |
| Trace formula / contour integral shortcut | FAIL | E | Contour ∮ f(s)ζ'/ζ(s)ds needs O(T^{1/2}) per eval (Riemann-Siegel); total O(T^{3/2} log T) = same as direct zero sum. Heat kernel smoothing trades accuracy for speed (no free lunch). First 50 zeros explain only 22% of oscillatory variance. See experiments/wildcard/spectral_shortcut.py | 20 |

## Other / Miscellaneous

| Approach | Verdict | Mode | Key Finding | Session |
|----------|---------|------|-------------|---------|
| Lambert W Pade | PARTIAL | I | 0.019% error, not exact | 1 |
| LCM/primorial | FAIL | C | Requires factoring | 4 |
| Totient summatory | FAIL | C | IS primality testing | 4 |
| Lucas/Fibonacci entry points | FAIL | C | Requires factoring exp-large numbers | 5 |
| Quadratic residue patterns | FAIL | E | Equivalent to Dirichlet's theorem | 5 |
| Wilson inversion / fast factorial | FAIL | E | O(p^{1/2}) per test | 5 |
| Wheel + CRT survivor enumeration | PARTIAL | - | O(1) per survivor but need too many | 5 |
| Smooth number counts | FAIL | I | 5-15% error, inherently approximate | 6 |
| Deep structure (8 experiments) | FAIL | I | No hidden periodicity | 6 |
| Function field analog (F_q) | FAIL | - | Genus finite -> O(g); for Z, genus = infinity. Quantified: virtual curve needs genus~sqrt(x), giving O(x^{3/2}). N_e(ln x)=PNT in disguise. Z[i]/Z[omega] no advantage. See experiments/wildcard/finite_field_lift.py | 7,19 |
| F_1 / Lattice / CFs | FAIL | I | q->1 degenerate: N_{1+eps}(n)->(n-1)*eps/n->0. No structure | 10,19 |
| Gap encoding | FAIL | I | ~4 bits irreducible per gap | 10 |
| Wavelet analysis (100K primes) | PARTIAL | I | alpha=1.68 spectrum, 1.2 bits/prime irreducible | 5 |
| Cramer model + corrections | FAIL | I | Parity of pi(x) from li: 49.3% (random) | 9 |
| Information-theoretic shortcut | FAIL | E | Zero sum is O(N), no FMM | 9 |
| Novel exact formulas (7 methods) | FAIL | E | All confirm barrier | 9 |
| Radical bypass (5 methods) | FAIL | - | Three failure modes identified | 9 |
| Internet search 2024-2026 | - | - | No new algorithms found | 5,8,9,10,13 |
| Wilson's theorem in TC^0 | FAIL | - | (n-1)! mod n needs 2^N mults, not poly(N). Shamir GapL only gives NC^2 | 13 |
| Sum-of-two-squares primality | FAIL | C | 173 composite counterexamples below 10000; needs factoring | 13 |
| Optimized li(x^{1/k}) coefficients | FAIL | I | Optimal coeffs overfit (test err 2000+); Riemann R(x) is essentially optimal; error grows as x^{0.3} | 13 |
| Zeta zero reordering/weighting | FAIL | E | K_min ~ 0.35 * x^{0.27} (power law). Greedy-optimal ordering varies by x; no universal shortcut below O(x^{1/4}) terms | 13 |
| Cayley graph spectral approach | FAIL | C | Graph requires knowing primes; lambda_0 = pi(x) trivially; eigenvalues = exponential sums over primes | 13 |
| Ihara zeta function of graphs | FAIL | E | For Cayley graph, Ihara zeta factors via character sums = discrete FT of prime indicator ≡ explicit formula | 13 |
| GCD/coprimality graph spectrum | FAIL | E | Constructible without primes BUT eigenvalues = Ramanujan sums c_q(n) involving μ and floor(x/d); same as Meissel-Lehmer O(x^{2/3}) | 13 |
| Expander mixing lemma for primes | FAIL | C+E | Either need primes to build graph (C) or GCD graph gives Meissel-Lehmer equivalence (E) | 13 |
| CRT reconstruction of pi(x) | FAIL | E | Each pi(x) mod q costs same O(x^{2/3}) as pi(x); k moduli needed → k × worse. Entropy invariant 0.537 bits confirmed | 13 |
| Recursive identity pi(x) via pi(x/d) | FAIL | E | For ANY fixed divisors d_1,...,d_k: correction function encodes primes in (x/max(d_i), x], as hard as pi(x). Tree of O(x^{1/3}) subproblems but correction is hard | 13 |
| Prime zeta P(s) Perron extraction | FAIL | E | P(s) = sum mu(k)/k * ln(zeta(ks)); Perron integral has poles at rho, rho/2, rho/3,... MORE terms than standard explicit formula | 13 |
| GF(2) algebraic shortcuts for primes | FAIL | I | Prime indicator ANF degree = Theta(N), 50% sparsity, indistinguishable from random over GF(2). Consistent with PRIMES ∉ AC^0[2] | 13 |
| Algebraic geometry point counting for pi(x) | FAIL | C+E | Genus must be Omega(x/ln(x)) to encode pi(x); H^1_et(Spec(Z)) infinite-dim; Weil duality = explicit formula; F_q[T] works only because genus=0 and zeta rational | 13 |
| Lucy DP matrix displacement rank | FAIL | E | Matrices are unipotent (all eigenvals=1), displacement rank 50-60% of dimension (NOT structured), product is full-rank with no low-rank approximation. No Toeplitz/Cauchy structure | 14 |
| Matrix power encoding of pi(x) | FAIL | E | pi(x) mod m NOT an LRS for any m=2..16; 2x2 matrix search gives only constant outputs; confirms Mauduit-Rivat | 14 |
| Polynomial in floor values | FAIL | I | Degree-2 poly in 20 floor values: 98% train exact but 16% test; catastrophic overfitting. Residual is smooth (autocorr 0.97) = zeta zero contribution | 14 |
| #L chain: PRIMES ∈ L → pi(x) ∈ NC^2 | FAIL | E | Chain breaks: NL machine needs O(N) workspace for candidate n but #L allows only O(log N). Work space mismatch is fundamental. PRIMES ∈ L and pi(x) ∈ NC are INDEPENDENT | 14 |
| Small matrix det = pi(x) search | FAIL | E | Redheffer full-rank, Lucy DP product dense. I-E fractional parts carry O(2^k) independent bits → no det smaller than 2^{sqrt(x)/2logx}. Must avoid floor functions entirely | 14 |
| LGV/DAG path count compression | FAIL | E | Sieve DAG has 2^k paths (exponential); recursive DAG has O(sqrt(x)) nodes; matrix product full-rank throughout. Floor-value set is IRREDUCIBLE state space | 14 |
| Residue class decomposition | FAIL | I | Primes equidistribute mod q (Dirichlet). Each L(s,chi) has own zeros, making problem HARDER. Character sums tautological for exact counting | 14 |
| Selberg sieve exact counting | FAIL | E | Upper bound only (3.7x actual). Parity barrier for linear weights. Nonlinear weights break parity but don't reduce cost | 14 |
| Batched primality for pi(x) | FAIL | E | No shared structure in 2^{n-1} mod n for consecutive n. Fermat count ≈ pi(x) (rare pseudoprimes) but computing it costs O(x) | 14 |
| Algebraic variety F_q point count (Session 14) | FAIL | E | Low-dim: too few points. High-dim (d≈N): Kedlaya O(N^3 * 2^N) = same as sieve. Frobenius eigenvalues = zeta zeros (same info). Constructing variety circular | 14 |
| Companion matrix / LRS fitting | FAIL | I | Order-20 recurrence: 27% test exact, errors grow with extrapolation. pi(x) fundamentally non-recurrent | 14 |
| Smooth/rough decomposition | FAIL | E | When B≥sqrt(x), rough numbers = primes (exact). But computing rough count IS the sieve. No shortcut | 14 |
| QFT deterministic (Grantham) | OPEN | - | QFT IS in TC^0; 4 PSPs below 50000; error < 1/7710. Deterministic correctness unknown without GRH | 13 |
| BPSW unconditional correctness | OPEN | - | BPSW IS in TC^0 as computation. No pseudoprime below 2^64. Proving correctness ⟺ PRIMES in TC^0 | 13 |
| Class number h(-d) as det entry | FAIL | E | h(-d) = (w*sqrt(d))/(2pi) * L(1,chi_d); constants, no x-dependence; equivalent to L-values | 16 |
| L-function L(1,chi) as det entry | FAIL | C+E | Partial L-functions need primes (C); full L-values global, no x-dependence; inversion = zeta zeros (E) | 16 |
| Elliptic curve a_p as det entry | FAIL | C | a_p defined only at primes; a_n for composite needs factoring; #E(Z/nZ) multiplicative, doesn't encode pi(x) | 16 |
| Number field regulators as det entry | FAIL | E | h*R = f(L(1,chi)) by class number formula; transcendental, incompatible with integer GapL entries | 16 |
| TC^0 MAJORITY fan-in for pi(x) | FAIL | E | Inputs (BPSW(k)) must be computed first; MAJORITY helps aggregation, not generation; size 2^N * poly(N) | 16 |
| Batch modular exp via CRT | FAIL | E | 2^{k-1} mod k: exponent-modulus coupling prevents batch. Fixed-base CRT computes different function | 16 |
| Carmichael lambda batch structure | FAIL | C | Identifying common lambda(k) requires factoring all k; no batch shortcut | 16 |
| Period exploitation in BPSW | FAIL | E | Autocorrelation of BPSW(n) ≈ 0 at all lags; no periodic structure | 16 |
| Lambda W error structure exploitation | FAIL | I | delta(n) uncorrelated with gaps (r~0.01); uniform mod m; power law ~sqrt(x); incompressible | 16 |
| Cheaper oracle hierarchy (R -> Lucy DP) | FAIL | I | R(x) error O(sqrt(x)) = same as search range; saves constants only, not asymptotics | 16 |
| H-T transfer via explicit formula zeros | FAIL | E | H-T is combinatorial; zero sum is analytic; incompatible techniques | 16 |
| Weighted prime counting for cancellation | FAIL | E | Multiplicative weights collapse to constants on primes; non-multiplicative lose sieve structure | 16 |
| Buchstab signed identity + H-T | FAIL | E | Signed Buchstab IS M(x); recovering pi(x) needs A_3, A_5 terms at O(x^{2/3}). The conversion M(x) -> pi(x) costs O(x^{2/3}) minimum | 16 |
| M(x) -> pi(x) conversion identity | FAIL | E | pi(x) = sum omega(d)*M(floor(x/d)); omega partial sums cost O(x^{2/3}). All conversion paths cost >= O(x^{2/3}) | 16 |
| Additive combinatorics / sumsets | FAIL | E | r(n) representations via circle method = zeta zeros; Goldbach-type identities don't give pi(x) | 16 |
| Ergodic theory / orbit complexity | FAIL | E+C | Transfer operators have spectral theory = zeta zeros; Furstenberg correspondence is circular | 16 |
| Model theory / o-minimality | FAIL | - | o-minimal structures can't represent step functions like pi(x); definability orthogonal to computation | 16 |
| Tropical geometry of P(s) | FAIL | I+E | Tropicalization loses all but smallest prime factor; tropical convolution = standard convolution | 16 |
| Tropical/min-plus full battery (7 tests) | FAIL | I | Tropical Dirichlet conv, Mertens, sieve, Euler product, p-adic val of x!, tropical det, tropical conv -- all fail. Core reason: min-plus is OPTIMIZATION not COUNTING; val(a+b)>=min(val(a),val(b)) destroys sum info. See experiments/proposals/tropical_sieve.py | 25+ |
| Smooth number Ψ(x,B) subtraction for π(x) | FAIL | E | Ψ(x,B) computation IS the sieve. 50-640x slower than Eratosthenes. Linear combo of 11 Ψ features: train RMSE=3.3, test RMSE=37.3 (0/500 exact). Buchstab tree has 1.8·x^{2/3} distinct args = Lucy DP. Legendre identity confirmed exact. See experiments/proposals/smooth_number_subtraction.py | 26 |
| Sufficient statistics of floor values | FAIL | I | Any T with |T|=poly(log x) determining pi(x) must encode O(N/2) bits; hashing loses exactness | 16 |
| Algebraic geometry F_q families | FAIL | E+C | Frobenius eigenvalues = zeta zeros; constructing the right variety requires knowing primes | 16 |
| Representation theory S_n/GL_n | FAIL | C+E | Character sums at primes = Dirichlet characters; Euler product through rep theory = L-functions | 16 |

## Session 17: Communication Complexity, Fourier Analysis, New Literature

| Approach | Verdict | Mode | Key Finding | Session |
|----------|---------|------|-------------|---------|
| Ono partition characterization | FAIL | C | n prime iff (n²-3n+2)σ₁(n)-8M₂(n)=0; requires divisors (circular); O(n²) per test, O(x³) total | 17 |
| GapL via multilinear polynomial | FAIL | - | dc(pi_N) >= 2^{N/2-1}+2 = Omega(sqrt(x)); substitution rank exponential in N; definitively closed | 17 |
| Boolean Fourier low-degree structure | FAIL | - | 30% excess low-deg weight = parity/mod-4 only; noise sensitivity near-random; no junta structure | 17 |
| Communication matrix rank shortcut | FAIL | - | rank(pi_N) = 2^{N/2-1}+2 exactly; converges to 50% of max; no efficient 2-party protocol | 17 |
| Partition generating function sum | FAIL | C | Sum [f(n)=0] has no GF shortcut; indicator of zeros doesn't factor through GFs | 17 |
| MacMahonesque generalizations | FAIL | C | Higher M_a also circular + more expensive; polynomials in n with periodic corrections | 17 |
| Additive number theory encoding | FAIL | - | Ono shows additive characterization exists but counting is equally hard from both viewpoints | 17 |

---

## Session 18: Preprocessing / Succinct Data Structures

| Approach | Verdict | Mode | Key Finding | Session |
|----------|---------|------|-------------|---------|
| Preprocessing + polylog query | FAIL | E | Reduces exactly to pi(x) complexity; polylog T requires M ~ 10^100 checkpoints (exponential storage); S*T tradeoff: S=T=(M*logx)^{2/5} unconditional | 18 |
| Range tree / skip structure for p(n) | PARTIAL | E | Optimal Delta=(332x)^{3/5}; balanced S=T~10^{41.8} for M=10^100; wins vs plain DR but not polylog | 18 |
| Cell-probe lower bounds for primes | FAIL | - | No prime-specific lower bounds exist; Patrascu rank/select gives only Omega(N) with poly(N) space (trivial) | 18 |
| Comm complexity -> data structure LB | FAIL | - | Session 17 rank=2^{N/2-1}+2 is for computing pi(x), not data structure SELECT; no non-trivial transfer | 18 |
| Non-uniform / geometric checkpoints | FAIL | E | Widest interval dominates; geometric checkpoints give O(N^2) storage but no asymptotic query improvement | 18 |

## Session 18: Streaming / Space Complexity

| Approach | Verdict | Mode | Key Finding | Session |
|----------|---------|------|-------------|---------|
| Streaming exact pi(x) | FAIL | E | Sieve-based streaming needs Omega(sqrt(x)/log(x)) space (must encode primes up to sqrt(x)); one-way CC from rank gives same bound; no non-sieve streaming algorithm known | 18 |
| Online p(n) with small state | PARTIAL | - | O(log n) state suffices (just store last prime, test next candidates with AKS); but time = O(n log n), catastrophic. Space is NOT the bottleneck | 18 |
| Space-time tradeoff for pi(x) (general) | OPEN | - | Communication complexity gives T*S >= Omega(log^2 x) only -- too weak to rule out polylog(x). Proving T >= x^{Omega(1)} requires circuit lower bounds (Natural Proofs barrier). OBDD model: Omega(sqrt(x)) proven. M-L DAG pebbling: T*S >= Omega(x^{5/6}/ln x) algorithm-specific. Nechiporuk: Omega(log x) trivial. See experiments/circuit_complexity/space_time_tradeoff_results.md | 20 |
| Space-time tradeoff via comm complexity | CLOSED | E | D(pi) = Theta(N/2) bits, giving T*S >= Omega(N^2) = Omega(log^2 x). This is polylog(x), so communication complexity CANNOT rule out polylog-time algorithms. The input has only N=log(x) bits, bounding D(f) <= N, so comm-complexity-based tradeoffs are inherently poly(N). | 20 |
| Space-time tradeoff via Nechiporuk | CLOSED | E | Nechiporuk formula bound uses rank(M_pi, k-bit) = 2^{k-1}+2 but optimal block size s*=3 gives only L(pi) >= Omega(N) = Omega(log x). Method inherently limited to O(N^2) for ANY function. | 20 |
| Space-time via OBDD/BDD size | PARTIAL | E | OBDD size ~ 2^{0.79*N} empirically (N=4..12); matches random functions. OBDD width >= 2^{N/2-1} proven from comm rank. But general BPs can be exponentially smaller than OBDDs. | 20 |
| Space-time via M-L DAG pebbling | CLOSED | E | Lucy DP DAG: depth=pi(sqrt(x)), width=O(sqrt(x)). Pebbling gives T*S >= Omega(x^{5/6}/ln x) but ONLY for Meissel-Lehmer; a non-sieve algorithm could bypass this DAG. | 20 |
| Branching program width for pi(x) | FAIL | E | Nechiporuk bound gives BP size >= 2^{N/2} = sqrt(x) from rank computation; confirms no NC^1 shortcut via BPs | 18 |

## Session 18: Self-Correction / Amplification / Boosting

| Approach | Verdict | Mode | Key Finding | Session |
|----------|---------|------|-------------|---------|
| Self-correction for pi(x) (BLR/Lipton) | FAIL | E+I | pi(x) not polynomial over any field; no algebraic identity at random points; Meissel-Lehmer IS the self-reduction but costs O(x^{2/3}) | 18 |
| Random self-reducibility of pi(x) | FAIL | E | Differences pi(x+r)-pi(x) are hard counting problems; no group/field structure over Z; even if RSR, proves hardness not easiness | 18 |
| Goldreich-Levin / list decoding for pi(x) | FAIL | I | R(x) has zero correlation with low ~170 bits of pi(x); comm rank 2^{N/2} blocks sparse Fourier recovery; no linear encoding helps | 18 |
| Sumcheck / arithmetization of pi(x) | FAIL | E | MLE of 1_P at non-Boolean points requires knowing all 1_P values (cost O(x)); varying modulus in primality tests blocks F_p arithmetization; efficient prover requires the circuit complexity breakthrough we seek | 18 |
| Local decodability of pi(x) | FAIL | C | All natural encodings (Dirichlet series, character sums, pi at other points) are hard to compute; encoding IS the bottleneck, not decoding; LDCs cannot help when encoding requires O(x^{2/3}) | 18 |

---

## Session 19: Communication Complexity, SVD Spectral, Identity Search

| Approach | Verdict | Mode | Key Finding | Session |
|----------|---------|------|-------------|---------|
| Unbalanced communication complexity | FAIL | - | rank = 2^{min(k,N-k)-1}+2 for ALL bit partitions k=1..N-1 (verified N=4..20). No polynomial-rank partition exists. Barrier is intrinsic. | 19 |
| 3-party NOF communication complexity | PARTIAL | - | Balanced (N/3,N/3,N/3): max cut rank = 2^{N/3}, NOF = N/3. Best (1,1,N-2): rank = 4 (trivial). Consistent with TC^0 but insufficient to prove/disprove. | 19 |
| k-party NOF (k=2..8) mode-unfolding rank | CLOSED | - | For k>=3: ALL mode-unfoldings have FULL RANK = 2^{sizes[i]}. Max rank = 2^{ceil(N/k)} exactly. Verified N=6..18, k=2..8 (35+ cases). For k=2: rank = 2^{N/2-1}+2 (sub-full). GF(p) rank identical for p=2,3,5,7. SVD: rank-1 dominates (99.9% variance). Residual (pi-smooth) also full rank. Mode-unfolding rank is too coarse to resolve ACC^0 question -- need tensor rank, discrepancy, or polynomial method. | 20 |
| SVD spectral decay of oscillatory part | FAIL | I | Power-law S~i^{-1} (not geometric). 90% osc variance in ~20 SVs but 99% needs ~30% of total. Max osc SV scales as x^{0.66}. Insufficient for exact computation. | 19 |
| SVD ↔ zeta zero correspondence | CONFIRMED | E | Top osc SVs match zeta zeros (corr 0.95 at N=20 for gamma_1). Zeta basis explains 0.12% of variance at N=20 (need more zeros for larger N). SVD IS explicit formula decomposition. | 19 |
| PSLQ/LLL identity search for f(x)=pi(x)-R(x) | FAIL | I | All 6 relation types tested (linear, polynomial deg 2-4, recurrence, modular, functional, discrete derivative) x=2..10000. All single-point PSLQ relations spurious (cross-validation residuals ~10^4). Fourier confirms zeta-zero dominance. | 19 |
| NFS-type L[1/3] for pi(x) (4 sub-approaches) | FAIL | E+C | Norm sieve = residue class counting (E); Chebotarev = circular + O(x^{1/2}) error (C+I); class groups = Euler products (E); Artin L-fns = explicit formula (E). NFS exploits multiplicative structure; pi(x) is additive-global. | 19 |
| Gap-based predictability for pi(x) | FAIL | I | AR(1..50) gives NO improvement over baseline. MI(g_n;g_{n+1})=0.38 bits (10.3%). Gaps 9% more compressible than i.i.d. Near Cramér random model. | 19 |
| Kt complexity of delta(n) empirical | FAIL | I | |delta|~n^{0.57}; bit ratio 0.52; AR(1) R²=0.996 but RMSE=10.5 (innovations random); uniform mod m; 18% compressible; sign run length 38.5 (zeta oscillation); no exploitable structure beyond smoothness | 19 |
| Fast-forwardable dynamical system on gaps | FAIL | I | AR(1..50) R^2<0 (worse than mean). Residue-conditional (g_n, p_n mod m) R^2<0.05. HMM(2-8 states) = baseline. MI(g_n;g_{n+1})=0.35 bits (8.9%). LZ complexity 71% of random. Correlation dim grows with embedding (no attractor). No substitution/morphism structure. p(n) mod m equidistributed, no autocorrelation. | 20 |

---

## Summary Statistics

- **Total approaches tested:** 641+
- **FAIL:** ~452 (confirmed impossible or impractical)
- **PARTIAL:** ~19 (works partially but doesn't meet target)
- **WORKS:** ~8 (correct but O(x^{2/3}) or worse, or TC^0 but conditional)
- **EXISTS:** ~3 (formulas exist but computationally useless)
- **OPEN:** 5 (circuit complexity TC^0/NC, BPSW correctness, QFT deterministic, #TC^0⊆NC?, S*T tradeoff)

### By Failure Mode
- **Circularity (C):** ~82 approaches
- **Equivalence (E):** ~222 approaches
- **Information Loss (I):** ~180 approaches

---

## Session 20: Fresh Perspective (10 wildcard experiments)

| Approach | Verdict | Mode | Key Finding | Session |
|----------|---------|------|-------------|---------|
| CRT modular π(x) mod m reconstruction | FAIL | C | π(x) mod m = random walk with step ~1/ln(x); no structure | 20 |
| Tensor network / MPS sieve | FAIL | I | Bond dim = primorial (volume-law entanglement); exponential | 20 |
| Automata-theoretic sieve (DFA product) | FAIL | I | Minimized DFA has exactly primorial(y) states; exponential | 20 |
| Trace formula / spectral shortcut | FAIL | E | Trace formula IS explicit formula (circular); O(T^{3/2}) | 20 |
| Heat kernel smoothing of zero sum | FAIL | I | Smoothing introduces O(√x) bias; no free lunch | 20 |
| Fast-forwardable dynamical system on gaps | FAIL | I | Gaps ~90% random; MI(g_n;g_{n+1})~0.3 bits; no attractor | 20 |
| HMM / substitution morphism for gaps | FAIL | I | LZ complexity near-random; no finite-state generator | 20 |
| Finite field lifting F_q→Z via q→1 | FAIL | E | q→1 degenerates to 0; virtual curve needs genus ~√x | 20 |
| F₁ (field with one element) for primes | FAIL | E | No meaningful polynomial ring; produces PNT not exact | 20 |
| Z[i]/Z[ω] prime counting shortcut | FAIL | E | Reduces to rational primes + congruence; MORE zeros | 20 |
| Sieve matrix SVD / low-rank | FAIL | I | Matrix has FULL rank = #primes; SVs don't decay fast | 20 |
| Lucy_Hedgehog S(v,p) compression | FAIL | I | Binary step function; Fourier NOT sparse; 90/10 split | 20 |
| Log-Fourier sieve (multiplicative → translation) | FAIL | I | 90% energy in 10 modes BUT 99% needs ~70% of modes | 20 |
| Iterative zero-sum with self-correction | FAIL | E | Sensitivity <1 (converges), but to value ~√x from exact | 20 |
| Wilson's theorem batch computation | FAIL | C | Computing (k-1)! mod k IS primality testing | 20 |
| Determinant/permanent sieve formulation | FAIL | E | Reduces to Möbius inclusion-exclusion with 2^k terms | 20 |
| Cyclotomic polynomial prime encoding | FAIL | C | Φ_n(1)=p iff n=p^k, but evaluation requires factoring | 20 |
| Arithmetic derivative (n'=1 iff prime) | FAIL | C | Perfect encoding but requires factorization | 20 |
| Number field sieve analog | FAIL | E | ζ_K has MORE zeros than ζ_Q; harder not easier | 20 |
| Matrix exponentiation for π(n) | FAIL | I | Non-autonomous system (M_n depends on primality); no fixed M | 20 |
| Selberg sieve as optimization / LP | FAIL | I | Parity barrier: sieves cannot distinguish set from complement | 20 |
| Dirichlet character decomposition | FAIL | E | Deviations O(√x/φ(q)); L-function zeros equally hard | 20 |
| L-function zero convergence rate | FAIL | E | Empirical: pi(x;q,a) error is ROUGHER than pi(x) (TV ratio 1.5-3x for q=3..7); spectral 90% power at K=39-81 vs K=148 for pi(x) but total zeros needed = phi(q) sets; principal char alone gives ~same error as zeta for pi(x); no convergence advantage | 20+ |
| Bit-by-bit computation of Δ=π(x)-Li(x) | FAIL | I | Even MSB requires O(√x) zeros; bits are entangled | 20 |
| Contour integral with Euler product | FAIL | E | Requires O(T^{3/2}) quadrature; worse than sieving | 20 |
| Euler product convergence for ζ | FAIL | E | Need primes up to ~x for accuracy; circular | 20 |

## Kt Complexity / Information Theory (Session 20 Deep Focus)

| Approach | Verdict | Mode | Key Finding | Session |
|----------|---------|------|-------------|---------|
| Predict delta(n) from n via polynomial | FAIL | I | R²<0, worse than naive predict-zero. Poly degree 6: RMSE 11x naive | 20 |
| Predict delta(n) from n via ML (RF, GBM) | FAIL | I | R²=-0.17 (RF), -0.14 (GBM). n is informationally irrelevant to delta | 20 |
| Predict delta(n) from n via Fourier series | FAIL | I | 200 Fourier terms: RMSE=55.5 vs std=176. Need ~10^4 params to match AR | 20 |
| AR + n-features | FAIL | I | Adding n as feature to AR(k) gives exactly 0% improvement | 20 |
| Predict delta(n) mod m from n | FAIL | I | At random baseline for all m=2,3,4,5,6 | 20 |
| Predict sign(delta(n)) from n | MARGINAL | I | Best 55.9% accuracy vs 50% random. Only from long sign runs | 20 |
| Explicit formula partial sums (K zeros) | FAIL | E | Diverges: error grows linearly with K at x<10^6 (conditional convergence) | 20 |
| Kt(delta) via incremental entropy | INCONCLUSIVE | - | Growth ~0.22*log(n) but this is statistical, not algorithmic Kt | 20 |
| Low-dim manifold for delta sequence | FAIL | I | Correlation dimension grows with embedding dim (no saturation) | 20 |
| SVD low-rank approx of delta | FAIL | I | S_k ~ k^{-1.1} power law, need 71% of SVs for 99.9% energy | 20 |
| Conditional entropy H(delta|n mod m) | FAIL | I | Even n mod 1000 only captures 31% of delta entropy | 20 |
| ML (Ridge/kNN) prediction of delta(n) | FAIL | I | 0% exact on 1000 test cases. Ridge RMSE improves 35% over baseline but still off by 190. kNN RMSE=78 but 0.2% exact | 21 |
| Fourier interpolation of S(x) via zero frequencies | FAIL | I | Test RMSE saturates at ~2.2 regardless of zeros used (5→30). Need O(sqrt(x)) zeros to resolve | 21 |
| F_1 extrapolation (q→1 limit of F_q counting) | FAIL | E | Richardson extrapolation DIVERGES. The limit is singular—not computationally transferable | 21 |
| Gaussian test function in Weil explicit formula | FAIL | E | Smoothing that kills zeros also destroys resolution. No optimum exists | 21 |
| Diophantine quasi-periodicity of zero ratios | FAIL | I | Super-period LCM=819 fundamental periods (~10^{158} multiplicatively). Zero ratios are generic irrationals | 21 |
| CRT reconstruction of delta(n) (direct) | PARTIAL | C | Framework works (unique candidate with M=30030) but computing p(n) mod q requires pi(x;q,a) at cost O(x^{2/3})—CIRCULAR | 21 |
| Bisection with R(x)-only oracle | FAIL | E | Error O(sqrt(x)) at convergence. Without zero info, oscillatory correction dominates | 21 |
| Prime autocorrelation mod W for prediction | FAIL | I | Strong lag-1 correlation (chi2/df=159) in p(n) mod 30 but only short-range (decays by lag 5). O(1) bits, not O(log n) | 21 |
| Prime race pi(x;q,a) shortcut for CRT | FAIL | E+C | pi(x;q,a) errors are 2x ROUGHER than pi(x) (TV ratio~2.0); each L-function needs O(sqrt(x)) zeros; phi(q) L-functions multiply cost; Chebyshev bias gives O(1) bits not O(log n). Katz-Sarnak universality: L-function zeros have GUE stats like zeta zeros | 22 |
| L-function zero convergence vs zeta convergence | FAIL | E | Per-character 90% power at K~39-81 vs K=148 for zeta, but phi(q) characters needed → total zeros >= zeta zeros. TV rougher. No faster convergence | 22 |
| k-party NOF communication (k≥3) | FAIL | E | Mode-unfolding rank is FULL (=2^{ceil(N/k)}) for ALL k≥3, indistinguishable from random functions. Mode-unfolding rank too coarse for ACC^0 question. Need true tensor rank or discrepancy | 23 |
| Communication complexity → time-space lower bound | FAIL | E | Communication complexity bounded by input length N=log(x), so can NEVER give super-polylog lower bounds. Max achievable: T*S >= Omega(log^2 x). Fundamental limitation | 23 |
| OBDD for pi(x) | FAIL | E | Size grows as 2^{0.79*N} — comparable to random functions. OBDD lower bound Omega(sqrt(x)) from communication rank, but only for OBDD model | 23 |
| Meissel-Lehmer pebbling space-time | PARTIAL | E | T*S >= Omega(x^{5/6}/ln x) for M-L DAG specifically. But algorithm-specific: different algorithm could bypass entirely | 23 |
| Holonomic (D-finite) recurrence for pi(n) | FAIL | I | NOT holonomic for any order d≤20 with polynomial degree r≤8. Test/random ratio ~1.0-1.7 across all (d,r). Stronger than "not LRS" (Session 14) — polynomial coefficients don't help either | 23 |
| Prime indicator holonomic | FAIL | I | Prime indicator (delta pi(n)) NOT holonomic. Same ratio to random as pi(n) itself | 23 |
| Ono partition characterization (p-adic lift) | FAIL | C+E | M_k(n) DP has O(n^2) ops — WORSE than O(x^{2/3}). Modular version same op count, only bounded values. Ono criterion mod l gives 46-72% accuracy (near random). Ramanujan congruences special cases only | 23 |
| Ono partition criterion (GF approach) | FAIL | E | M₁(n) via divisor-partition convolution: O(n^{3/2}). Direct enumeration: O(exp(π√(2n/3))). Both worse than Meissel-Lehmer O(x^{2/3}). Partition functions inherently require O(n) terms minimum. | 36 |
| Short-interval explicit formula iteration | FAIL | E | Each round needs K_i ~ K_1 zeros (sinc envelope gives x/W cutoff). Iteration does NOT reduce zero count. Hybrid optimum at W=sqrt(x) gives O(sqrt(x)*polylog) — matching known best | 23 |
| Approx degree of isPrime over GF(p) | PARTIAL | E | GF(2/3): trivial deg 0 for N≥8 (prime density < 1/(2p)). GF(5): deg grows ~2N/3 up to N=12, then full N at N=14. GF(7): similar. Linear growth CONSISTENT with ∉ACC^0 but only upper bounds (ANF truncation). Razborov-Smolensky breaks down for sparse functions | 23 |
| Real-valued approx degree of chi_P (polynomial method) | CLOSED | I | adeg_0.49(chi_P) = ceil(N/2) for N=4..11 (LP-verified). Power-law fit: 0.601*N^0.909 = Theta(N). adeg(pi(x)) = adeg(chi_P) exactly. Error = 0.5 for deg < N/2, then exponential drop (phase transition). Promise version (coprime to 2,3,5) reduces by only 1-2 degrees. Quantum lower bound: Q(chi_P) >= N/4 queries. SOS degree = adeg. PARITY-like, not MAJORITY-like. See experiments/circuit_complexity/approx_degree_prime_results.md | 28 |
| Zeta oracle query complexity | FAIL | E | M(x) ~ x^{alpha} evaluations needed, alpha ∈ [0.25, 0.50]. Three arguments: (A) K_min zeros scales as x^alpha; (B) N(sqrt(x)) ~ sqrt(x)*log(x) zeros needed, each O(1) eval; (C) condition number O(sqrt(x)). Oracle model CONFIRMS sqrt(x) barrier | 23 |
| CRT Prime Locator (progression counting) | FAIL | C | π(x;q,a) as hard as π(x). L-function zeros same barrier. Only 4-5 CRT moduli needed but each requires full prime counting | 24 |
| Hierarchical sieve (fast-multipole Φ) | FAIL | E | Φ calls/x = 0.03-0.04 constant (linear, not polylog). Prime indicator Fourier: 99% energy needs 78% of coefficients. Möbius sparsity doesn't scale | 24 |
| Spectral compression of zero sum | FAIL | I | 50 zeros: residual 1.84 at x=10000. DCT best: 99% in 10.4% coefficients. Critical 5% dense in ALL bases. Fourier interpolation: only 2x savings | 24 |
| Prime gap linear recurrence | FAIL | I | R² negative for all orders 1-20 (worse than mean predictor). Next gap ≈ 10±8 regardless of history | 24 |
| k-automatic prime indicator | FAIL | I | 2-kernel has 38 growing subsequences (vs 6 for Thue-Morse). NOT k-automatic for any k | 24 |
| LFSR encoding over finite fields | FAIL | I | L/N = 0.5000 over ALL GF(p) tested (p=2..23). Maximally random-like in every field | 24 |
| Recursive Dickman DDE shortcut | FAIL | E | Correction grows with exponent α≈7.9, dominates signal. 93-100% recursion nodes need exact computation | 24 |
| Neural/ML for δ(n) | FAIL | I | Test RMSE=3.44 (random features), 1.3-1.6 (zeta features). Never <0.5. No generalization | 24 |
| Étale/Weil formula analog | FAIL | E | Needs genus ~10^102 for target. Frobenius eigenvalues cost O(g²) = O(10^204) | 24 |
| Connes trace formula optimization | FAIL | E | σ=0.1 optimal (5 zeros for 99%), but smoothed π(x) error O(√x/log x). No σ gives <0.5 | 24 |
| NTT/Dirichlet convolution for sieve | FAIL | E | Perron needs T~x, O(T) quadrature. Floor grouping: 2√x−1 values, reproduces O(x^{2/3}) | 24 |
| Binary search + local sieve | FAIL | C | Works mechanically but requires π(x) oracle at window boundary. Problem relocated, not solved | 24 |
| Adelic local-global reconstruction | FAIL | C | Only 2 moduli suffice for CRT, but computing p(n) mod p^k needs Dirichlet L-function zeros | 24 |
| Kolmogorov complexity of δ(n) | PARTIAL | I | 5x more compressible than random (0.049 vs 0.256). Improves with length (1.56 bits/symbol at N=10000). Shannon entropy ratio 0.84. Structure exists but insufficient for exactness | 24 |
| Multiplicative Fourier / characters | FAIL | E | L(1,χ) polylog via digamma but gives density not count. Summing L-function zeros same barrier | 24 |
| Sublinear spectral correction (DCT+zeros) | FAIL | C+I | 17.3% of DCT coefficients needed for rounding (867/4999). Coefficients unpredictable across scales (correlation 0.742). Hybrid zeros+DCT circular. FRI sampling → 100% at scale | 24 |
| Tensor network (MPS) Legendre sieve contraction | FAIL | I | Bond dim ~ 2^(0.33-0.43 * a), exponential in #primes. Sieve 5-20x more structured than random but still exponential. Dominant SV is rank-1 (floor(x)), corrections full-rank. O(sqrt(x)) unique floors don't help: non-local assignment to bipartitions. Confirms S10/S20 | 26 |
| Explicit BDD circuit synthesis for pi(x) | PARTIAL | E | ROBDD with multi-ordering: LSB BDD ~ 2^(0.73*N). Better than OBDD 2^(0.79*N) (S20) but worse than sqrt = 2^(0.5*N). Per-bit: MSB~N+1 (trivial), LSB~2^(0.73*N). Influence(LSB)~N/2. BDDs are restricted (branching programs); general circuits could be exponentially smaller. Does NOT prove circuit lower bound. | 28 |
| Ramanujan Library PSLQ on delta(n)=p(n)-R^{-1}(n) | FAIL | I | PSLQ on n=1..2000, validated n=2001..2500. Linear recurrence orders 2-20: all relations unique (spurious). Polynomial deg 2-3 in consecutive deltas: all unique. Mixed delta(n) vs log(n)/sqrt(n): all unique. Modular: brute-force mod 2 order-4 gives 70% but trivial predictor gives 92% (autocorrelation artifact). Berlekamp-Massey: LFSR length L=N/2 for all mod 2,3,5,7,11 (maximally complex). Differenced delta same. Entropy 6.26 bits. Confirms algebraic independence of delta sequence | 26 |
| Group-theoretic cross-modulus prime race shortcut | FAIL | E+C | E(x;q) correlations between moduli are trend artifacts (QR/NQR imbalance); oscillatory parts independent (roughness ratio 0.997). Galois lifts go wrong direction (fine->coarse free, coarse->fine needs new zeros). Linear prediction of p(n) mod q from E(x;q_i): R^2 < 0.03, MI < 3% of max. Character orthogonality makes cross-modulus shortcuts mathematically impossible. CRT via prime races costs 96*O(sqrt(x)) for q<=23, worse than direct pi(x) | 26 |
| Wheel decomposition circuit complexity | FAIL | I | Per-class pi_r(x) circuits 3-4x MORE complex (normalized transitions) than full pi(x). Classes near-independent (I/H < 0.01), preventing shortcuts. Total entropy grows as phi(M). Entropy reduction from mixed-radix conditioning is finite-size effect (shrinks as N grows). Decomposition destroys sequential regularity that makes pi(x) partially structured | 28 |
| Multiplicative-additive circuit structure (Legendre I-E cancellation) | FAIL | I | Signed matrix rank = #distinct floors (full rank, no reduction). 90% of floor values have nonzero net I-E contribution. Effective terms O(sqrt(x)) but no way to evaluate sum sublinearly | 28 |
| Carry propagation structure in pi(x) sum | FAIL | I | Carry chain distribution IDENTICAL to random Bernoulli(1/ln x). Max carry = O(log pi(x)) = O(N). No exploitable additive structure | 28 |
| Monochromatic rectangle partition of pi(x) | FAIL | I | Partition number ~ 2^(0.76*N), exceeding rank ~ 2^(0.41*N). Nondeterministic communication complexity exponential. partition/rank ratio grows as 2^(0.35*N) | 28 |
| Per-bit complexity gradient analysis | DIAGNOSTIC | I | LSB-half influence 2x MSB-half, ratio growing with N. Bit 0 influence ~ N/2, max sensitivity = N. MSB influence O(1). Crossover at bit N/2. All CONSISTENT with polylog circuits | 28 |
| Approximate degree of chi_P at rounding threshold | DIAGNOSTIC | -- | adeg(chi_P, 0.49) = ceil(N/2) for N=4..10. Same for pi(x) mod 2. Counting adds no difficulty. Quantum lower bound Omega(N/4), still polylog | 28 |
| Rounding boundary / frac(R(x)) analysis | FAIL | I | frac(R(x)) perfectly uniform; precision follows geometric distribution P(k bits)=2^{-(k-1)}; no easy subset. R(x) accuracy drops to 13% at N=16 | 28 |
| Legendre I-E cancellation structure | FAIL | E | Signed matrix rank = distinct floors (full rank). 90% nonzero net. Max collision grows but doesn't cancel | 28 |
| Linear dynamical systems over Z/MZ | FAIL | I | p(n) mod M has no linear recurrence dim<=4 for M=3..30. No affine recurrence. Random matrix orbits match at chance level. Only trivial mod 2 (all odd) | 29 |
| GF(2) sieve matrix product fast-forward | FAIL | E | Product of sieve matrices is DIAGONAL (commuting). Encodes Mobius parity = Legendre sieve in disguise. Rank = #primes+1s. No matrix structure shortcut | 29 |
| Polynomial recurrence mod m (deg<=3) | FAIL | I | No degree-1/2/3 polynomial recurrence p(n)=f(p(n-1),p(n-2)) mod M for M=3,5,7. Only trivial mod 2. Extends S24 LFSR to nonlinear | 29 |
| CF / Stern-Brocot structure of p(n)/n | FAIL | I | CF partial quotients follow Gauss-Kuzmin (random). SB path lengths O(log p(n)). p(n)/n behaves as generic real | 29 |
| CRT from modular periods of p(n) | FAIL | I | p(n) mod m NOT periodic for m>=5 (within 500 terms). Only mod 2/3/4/6/12 trivially periodic (all primes odd, Dirichlet). Cannot fast-forward residues | 29 |
| Cipolla residual autoregression | FAIL | I | AR(1) 91.4% reduction but irreducible error grows as O(log n) = O(gap std). AR(20) gives no improvement over AR(1). Cannot achieve O(1) for exact | 29 |
| GUE CLT for zero sum | FAIL | I | CLT gives O(sqrt(log x)) typical size but actual correction is O(sqrt(x)/log(x)). CLT applies to random ensembles, not specific zero config | 29 |
| Zero grouping/clustering | FAIL | I | Group size 2→error 9.65, group 4→0.51 (lucky), group 8→2.65 at x=10K. Rapid oscillation of exp(i*gamma*log(x)) makes each zero position essential | 29 |
| GUE surrogate replacement | FAIL | I | Random matrix surrogates give mean=-0.81±1.18 vs actual -0.25. Right order of magnitude, wrong value. Specific zero config matters | 29 |
| Interpolation of pi correction | FAIL | I | 100 Chebyshev nodes on [2,10000]: max error=7. Fine structure from zeta zeros prevents polynomial interpolation with few points | 29 |
| Polynomial empirical correction to R(x) | FAIL | I | Polynomial in 1/log(x) reduces error 60-80% but cross-validation confirms overfitting. Oscillatory part genuinely non-polynomial | 29 |
| Compressed sensing on prime indicator | FAIL | C | pi(x)-x/log(x) has 99% energy in 1.25% of Fourier components, but those ARE zeta zero frequencies. Circular | 29 |
| Transfer matrix / stat mech for sieve | FAIL | E | State space 2^{pi(sqrt(x))}—exponential. Wheel sieve already exploits small-prime periodicity. No stat-mech shortcut | 29 |
| p-adic interpolation of pi(x) | FAIL | I | pi(x) NOT p-adically continuous: pi(x) mod p^k varies across all residue classes of x mod p^k. Local info doesn't determine global count | 29 |
| Wilson's theorem bulk primality | FAIL | E | Running factorial requires O(x!) precision. Computing (k-1)! mod k for single k is O(k). No bulk shortcut | 29 |
| Goldbach representation extraction | FAIL | C | Computing r_2(n) requires knowing primes. Hardy-Littlewood formula error O(sqrt(n/log(n))) too noisy | 29 |
| Finite field F_q[x] analogy lift | FAIL | I | F_q[x] works because its zeta has ONE zero. Z has infinitely many. No deformation from F_q to Z preserves the counting formula | 29 |
| PSLQ identity search | FAIL | I | No integer relation found among pi(x), li, R, sqrt, log, zeta values with coefficients ≤1000. Known asymptotics are the only relations | 29 |
| Primorial decomposition optimization | FAIL | E | Optimal c=3 in x^{1/c} split → O(x^{2/3}). This IS Meissel-Lehmer. No c gives better tradeoff | 29 |
| Short-interval error structure | FAIL | C | Short-interval Fourier structure from zeta zeros. CS recovery requires zero frequencies a priori | 29 |
| Number-theoretic hash p(n) mod m | FAIL | I | p(n) mod m for m≥3 shows no periodicity, no autocorrelation, pseudorandom. No polylog-computable formula | 29 |
| DE search for f(x)=pi(x)-R(x) | FAIL | I | Tested linear ODEs order<=3 poly deg<=3, Euler-type, nonlinear poly ODEs, Volterra integrals. Best linear ODE residual matches random noise (spurious from smoothing). Euler/nonlinear residuals ~1.0. f encodes infinitely many zeta-zero oscillations incompatible with finite-order DE | 29 |
| Extended PSLQ identity search x=2..100000 | FAIL | I | 14-element basis (log,sqrt,roots,li-variants,zeta zero oscillations). 18 PSLQ tests, 15 relations with nonzero f-coeff -- ALL fail cross-validation (residuals 13-53000). Functional relations f(ax) vs f(x) for a=2,3,4: fail. Shift recurrences f(x)..f(x+10): fail. Extends S19 to 10x range | 29 |
| Wilf-Zeilberger definite sum for f(x) | FAIL | I | Delta_f bimodal (prime indicator). Higher-order differences RMS GROW (ratio->2.0=white noise). Hypergeometric recurrence R^2=0.997 is spurious (trivial autocorrelation). Summation kernel full Hankel rank 250/250 (incompressible). No WZ certificate | 29 |
| f(x) vs Bernoulli numbers | FAIL | I | Bernoulli corrections < 10^{-7} at x=1000, zero correlation with f(x) (r=-0.006). Trivial zeros contribute negligibly | 29 |
| f(x) vs zeta values zeta(2..7) PSLQ | FAIL | I | PSLQ finds x-dependent relations with large coefficients (~10^6). Different coefficients at each x. Perfect-square x gives trivial sqrt relations. No universal algebraic relation | 29 |
| f(x) vs Dirichlet L-values L(1,chi) | FAIL | I | Same as zeta values -- x-dependent large-coefficient relations only. L(1,chi_3)=0.605, L(1,chi_4)=pi/4. No universal relation | 29 |
| f(x) vs Ramanujan tau function | FAIL | I | Pearson r=+0.010, p=0.93. Spearman r=+0.054. Zero correlation. Different mathematical objects (weight-12 modular form vs zeta zero oscillations) | 29 |
| Chebyshev psi(x)/log(x) as f(x) shortcut | FAIL | E | g(x)/log(x) captures 91% of f(x) variance (r=0.996) via partial summation identity. BUT psi(x) requires O(x) computation -- worse than O(x^{2/3}) sieve. Known identity, not shortcut | 29 |
| LLL minimal polynomials for f(x) | FAIL | I | "Candidate" polynomials at every x but DIFFERENT polynomials at each x (x=100: [156396,287371,96731], x=1000: [119569,141263,-381961]). Multi-point relations at float64 precision limit. f(x) values effectively algebraically independent transcendentals | 29 |
| Polynomial-in-log(x) fit for f(x) | FAIL | I | f(x) ~ sum a_j*log(x)^j: best validation RMSE=1.976 (~47% of std). With sqrt(x) factor: no improvement. Degree>4 overfits. f(x) not a polynomial in log(x) | 29 |
| Incremental delta(n) via autocorrelation | FAIL | I | r(1)=0.975 but smooth-trend artifact (after MA-10 detrending: r(1)=0.41). AR(1) RMSE=7.31, AR(5) RMSE=7.27, both >> 0.5 threshold. Only 5.3% of predictions within rounding range. H(delta(n+1)\|delta(n))=4.93 bits/step. RMSE grows with n (4.8->8.1). See experiments/proposals/critique_incremental_delta.py | 30 |
| Verification-prediction separation | FAIL | C | Primality testing is O(polylog) but ordinality (is g the n-th prime?) requires pi(g) = O(x^{2/3}). Claim conflates primality with ordinality. Local sieving gives interval count but not global rank. At x=10^100: primality ~10^9 ops, ordinality ~10^68 ops. See experiments/proposals/critique_verification_separation.py | 30 |
| NOF discrepancy for chi_P (3-party) | FAIL | - | Discrepancy HIGH (~bias), dominated by density imbalance. Trivial cylinder (A=B=C=all) achieves max. No communication lower bound. Degree-1 Walsh correlation 3-1513x random (parity), degree ≥ 2 matches random. Consistent with TC^0. See experiments/circuit_complexity/nof_discrepancy_chi_P.py | 31 |
| Gamma-2 (factorization) norm / sign-rank | FAIL | - | γ₂ ~ 2^{0.186N}, 85-93% of random matrices. Sign-rank = rank for all N=4..10. No hidden low-rank sign structure. SM complexity bound O(log γ₂) = O(N) trivially weak. F_2 rank = rank_R - 1. See experiments/circuit_complexity/gamma2_norm_chi_P.py | 31 |
| True tensor rank of chi_P (3-way) | DIAGNOSTIC | - | Tensor rank ~ d^{1.5} = 2^{N/2} = sqrt(x). N=6: rank=5 (random=7), N=9: rank≤19 (random=29), N=12: rank≈67 (random>100). chi_P 25-35% BELOW random but still EXPONENTIAL. Confirms N/2 universality. Close to generic d²/3 bound. See experiments/circuit_complexity/tensor_rank_robust.py | 31 |
| BDD size / sensitivity / decision tree | FAIL | - | BDD grows as 1.661^N (R²=0.997). chi_P ~30% simpler than random (ratio 0.69-0.71). Sensitivity = block sensitivity = certificate complexity = N for N≥7. Decision tree depth = N. All N variables essential. See experiments/circuit_complexity/min_circuit_size.py | 31 |
| F_2 correlation profile (degree-d polynomials) | DIAGNOSTIC | - | Degree-0,1 dominate: W(0)+W(1) = 47%(N=6) to 68%(N=16). W(1) z-score grows from 3 to 1513 (parity). W(d≥2) BELOW random (z-scores -2 to -30). After bias+parity removal, chi_P more pseudorandom than random. No exploitable low-degree F_2 structure. See experiments/circuit_complexity/f2_correlation_profile.py | 31 |
| Recursive prime counting (FMM-inspired) | FAIL | E | Aggressive recursion achieves O(log log x) depth but work catastrophically concentrated at depth 0: Theta(x^{2/3}/ln x). phi recursion tree has O(x^{2/3}/ln x) leaves. FMM analogy fails: no smooth kernel. Hierarchical sieve: large primes make O(x^{2/3}) irregular removals. Reduces to Meissel-Lehmer. See experiments/wildcard/recursive_prime_counting.py | 32 |
| Tropical/min-plus fast-forward via gap structure | FAIL | I | Hankel rank 418/500 (ratio to random: 0.98). Corr dim grows with embedding (2.1→5.8, no attractor). Linear R²<0.001, quadratic R²<0.002. MI(g_n;g_{n+1})=0.013 bits (0.4% of entropy). Min-plus accuracy 7-11%. Algebraic obstruction: p(n)=sum(gaps) needs (+,×) not (min,+). Gaps indistinguishable from IID. See experiments/wildcard/tropical_prime_gaps.py | 32 |
| Trace formula / moment method for zero sums | FAIL | I | Moment expansion S(x)=Σ(lnx)^k/k!·M_{k-1} DIVERGES: |M_k|~γ_max^k>>k!, radius of convergence x<1.0007. GUE trace analogy breaks: bounded eigenvalues required. Weil geometric side circular (needs all primes to x). Effective rank ~N (no compression). Gaussian damping unfavorable trade-off. 4 independent failures. See experiments/wildcard/trace_formula_approach.py | 33 |
| Contour integral evaluation of zero sum | FAIL | E | Contour ∮(ζ'/ζ)(s)·x^s/s ds equivalent to zero summation by residue theorem. Nyquist: quadrature points ≥ enclosed zeros. Smoothing trades accuracy for fewer zeros (uncertainty principle). SVD of zero-contribution matrix: 99% energy needs 133/500 components. No sparse representation in any basis. See experiments/wildcard/contour_integral_zerosum.py | 32 |
| Hybrid analytic + local sieve | FAIL | E | R⁻¹(n) error O(√p(n)) with ratio 0.05-0.25. Truncated zero sum often WORSE than none (K=0 error 12 vs K=500 error 208 at n=50000). Optimal hybrid complexity O(x^{1/2+ε}) matches Lagarias-Odlyzko. Iterative refinement oscillates. Cannot beat O(√x) zeros required. See experiments/wildcard/hybrid_analytic_sieve.py | 32 |
| Hilbert-Pólya trace / GUE model / functional equation | FAIL | E,I | GUE: correct statistics, wrong values (std≈20-60). Moments: M_k~γ_max^k diverges. Selberg trace: circular (spectral↔geometric duality). Functional equation: doesn't change kernel. δ(x) incompressible: poly deg-50 RMSE=0.21≈deg-2. 90% Fourier power needs 138/500 modes. See experiments/wildcard/hilbert_polya_trace.py | 32 |
| Convergence acceleration of zeta zero sum | FAIL | I | Richardson, Euler-Maclaurin, Padé, Cesàro, Aitken Δ², Shanks tested. Errors GROW as ~N^{0.8-1.0} (random walk). Best: Shanks(3) gives ~10x constant improvement. Root cause: each zero carries independent info (GUE-random phases). No structured error terms to exploit. See experiments/wildcard/zero_sum_acceleration.py | 32 |
| GRH batch Miller testing for pi(x) | FAIL | E | GRH Miller (witnesses ≤2ln²n) correct to 10^5. Batch testing: 1.02x speedup (negligible). Cost O(n·polylog n) — worse than sieve O(n·loglog n) and much worse than Meissel-Lehmer O(x^{2/3}). Each number requires independent modular exponentiation. See experiments/analytic/conditional/grh_miller_batch.py | 33 |
| GRH explicit formula optimal T | FAIL | E | Under GRH, need T=O(√x·log²x) zeros for error<0.5. With 1000 zeros, error≥18 for x≥10^4. Individual zeros at γ~200-400 still contribute >0.5 to pi(10^6). No way to use fewer than O(√x) zeros. Confirms Lagarias-Odlyzko at O(√x·polylog x). See experiments/analytic/conditional/grh_miller_batch.py | 33 |
| Elliott-Halberstam for pi(x) counting | FAIL | E | EH controls distribution in residue classes, not total count. Summing li(x)/φ(q) over coprime residues gives li(x) (same error as direct). Residue-class method 10^3-10^6x SLOWER than direct primepi. EH is about equidistribution, not enumeration. See experiments/analytic/conditional/elliott_halberstam_gaps.py | 33 |
| Gap structure (Cramér) for p(n) | FAIL | I | Cramér model fits well (max gaps 40-60% of ln²p). R⁻¹(n) search interval has O(ln x) primes. BUT identifying which is p(n) requires pi(x) at boundary — the counting bottleneck. At x=10^100: counting 10^57x more expensive than searching. Cramér solves wrong problem. See experiments/analytic/conditional/elliott_halberstam_gaps.py | 33 |
| Schoenfeld interval sieve under RH | FAIL | E | |pi(x)-li(x)|<√x·ln(x)/(8π). Sieve cost O(√x·ln x·ln ln x)=O(x^{1/2+ε}), better than ML O(x^{2/3}) asymptotically but still exponential in input bits. Iterative refinement with 1000 zeros: only 21% error reduction at x=10^6. Need O(√x·log²x) zeros for exactness. See experiments/analytic/conditional/schoenfeld_cramer.py | 33 |
| Cramér search + counting bottleneck | FAIL | E | Walk phase O(ln⁴x) trivial (0-114 steps verified). Counting phase pi(x₀) grows to dominate: 25% at n=5M, heading to 100%. At x=10^100: count=10^66.7 ops vs walk=10^9.5 ops. No conjecture reduces counting below O(x^{1/2+ε}). See experiments/analytic/conditional/schoenfeld_cramer.py | 33 |
| Best conditional algorithm (all conjectures) | FAIL | E | Under RH+Odlyzko-Schönhage: O(x^{1/2+ε}). RH alone with Turing zeros: O(x^{2/3+ε}), WORSE than unconditional. GRH, EH, Cramér all fail to improve beyond O(x^{1/2+ε}). For p(10^100): ~10^51 ops minimum. The √x barrier is fundamental — no standard conjecture breaks it. See experiments/analytic/conditional/best_conditional_algorithm.py | 33 |
| Spectral truncation / adaptive zero selection | FAIL | E | Select "resonant" zeros (γ·log(x) near kπ). For n<100, 1-2 zeros suffice (small-number artifact). For n≥500, 1000 zeros insufficient. Adaptive selection 41% better than sequential but cannot reduce O(√x) requirement. = reordering explicit formula sum. See experiments/proposals/proposal16_spectral_truncation_bound.py | 33 |
| p-adic lifting via CRT for pi(x) | FAIL | C+E | CRT works (3-6 primes suffice) but pi(x) mod p costs O(x^{2/3}) per prime. floor(x/d) mod q ≠ floor((x mod q)/(d mod q)) — floor division not a ring homomorphism. 12th CRT variant tested. See experiments/proposals/proposal17_padic_lifting.py | 33 |
| Dequantized Grover counting for pi(x) | FAIL | I | Tang 2024/Chia 2025 dequantization requires low-rank input. Prime indicator has adeg=N/2, comm rank=2^{N/2-1}+2 — full-rank. Truncated Möbius sums diverge. Sampling variance O(x²/S). Hyperbola=known O(√x). See experiments/proposals/proposal18_dequantized_grover_count.py | 33 |
| PSLQ/LLL on delta(n) (Ramanujan Library) | FAIL | I | delta(n) range [-130,102], std~33, lag-1 autocorr 0.95. PSLQ finds DIFFERENT relations per n (no universal formula). No recurrence order 1-6. No modular patterns. Best basis RMSE=31. 7th PSLQ/identity variant tested. See experiments/proposals/proposal19_ramanujan_library_delta.py | 33 |
| Étale cohomology point-counting for pi(x) | FAIL | C+E | Grothendieck-Lefschetz is polylog for FIXED variety, but encoding pi(x) needs unbounded dimension. Frobenius eigenvalues = zeta zeros. Character sums cost O(x). Hasse-Weil gives enough bits but each costs O(x^{2/3}). EC trace correlation 0.39 — too weak. 3rd algebraic geometry variant tested. See experiments/proposals/proposal20_etale_cohomology_count.py | 33 |
| MKtP / meta-complexity framework | FAIL | E | Kt framework reformulates "Is pi(x) in NC?" as "Is Kt(pi(x)mod2\|x) = O(polylog)?" — equivalent, not easier. Brandt's MKtP↔circuit lower bounds is GENERIC (any function in E, not pi(x) specifically). Kt(T_N) = O(2^N·N) by sieve regardless of circuit size. No new technique beyond circuit complexity theory. See experiments/circuit_complexity/meta_complexity_analysis.py | 35 |
| Approximate circuit complexity phase transition | FAIL | I | NO phase transition in rank-k approximation accuracy. Accuracy degrades GRADUALLY (rank-2: 63%, rank-10: 83%, rank-20: 95% at N=14). Errors spatially uniform (CV=0.22-0.25). No "easy core" to exploit for approx-then-correct strategy. See experiments/circuit_complexity/approx_circuit_complexity.py | 35 |
| Truth table compressibility for pi(x) | FAIL | I | 2.6x more compressible than random at N=20 (gzip 0.385 vs 1.0) but only LOCAL structure (prime sparsity creates runs). Block entropy → max at all scales. Shannon entropy → 1.0 bit/entry. Doesn't imply small circuits. See experiments/circuit_complexity/meta_complexity_analysis.py | 35 |
| GF(2) SLP / common subexpression analysis | FAIL | I | ANF sparsity EXACTLY 0.50 (random). CSE savings 50-66% = SAME as random functions. SLP length Θ(2^N). No GF(2) algebraic compression. Variable frequencies uniform. Monomial patterns at each degree → 0.50. pi(x) mod 2 indistinguishable from random Boolean function in ANF structure. See experiments/circuit_complexity/gf2_slp_structure.py | 35 |
| Smooth approximation for pi(x) parity | FAIL | I | R(x) accuracy for pi(x) mod 2 → 0.50 (random) as N grows. Smooth part provides ZERO parity information. Top-2 SVs capture only 12% of variance at N=14 (decreasing as 1/N). Parity entirely determined by oscillatory (zeta zero) contributions. See experiments/circuit_complexity/approx_circuit_complexity.py | 35 |
| Unbalanced communication matrix rank | DIAGNOSTIC | - | For k_LSB > N/2, actual rank EXCEEDS formula 2^{min(k,N-k)-1}+2 — often full row rank. Formula only tight for k ≤ N/2. SVD spectrum flat (power-law decay, no gap). Each additional SV contributes small equal increment. Polynomial rank for k=2logN but component functions g_i(MSBs) are as hard as pi(x). See experiments/circuit_complexity/approx_circuit_complexity.py | 35 |
| Depth-2 threshold circuits for pi(x) LSB | FAIL | I | PTF degree = N/2 exactly (LP-verified N=4-12, ratio 0.50 stable from N=6). Requires C(N,N/2) ~ 2^N/sqrt(N) monomials = EXPONENTIAL. Single LTF accuracy → 0.50 (random). Depth-2 heuristic with k=64 gates fails for N≥8. Matches random function PTF degree (Gotsman 1994). Does NOT rule out poly-depth TC^0. See experiments/circuit_complexity/threshold_circuit_construction.py | 35 |
| Adelic local-global prime counting (CRT mod q) | FAIL | E | 6th CRT variant. pi(x) mod q via AP decomposition 100-1000x SLOWER than pi(x) direct. R(x) error grows as O(sqrt(x)), mod-q correctness degrades to 0/6 at x=10^6. Liouville/Mertens parity ~50% (random). floor division not a ring homomorphism kills all modular decomposition. 3-6 moduli suffice for CRT but each costs O(x^{2/3}). See experiments/wildcard/adelic_prime_count.py | 35 |
