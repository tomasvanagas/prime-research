# Closed Paths: Master Lookup (445+ Approaches)

**SEARCH THIS FILE before proposing any approach.** Use grep/ctrl-F.

Last updated: 2026-04-04 (Sessions 1-15, 115+ sub-agents)

## Failure Modes
- **C** = Circularity (needs primes to compute primes)
- **E** = Equivalence (reduces to zeta zero sum / explicit formula)
- **I** = Information Loss (smooth approximation loses critical bits)

---

## Analytic / Explicit Formula / Zeta Zeros

| Approach | Verdict | Mode | Key Finding | Session |
|----------|---------|------|-------------|---------|
| Explicit formula + few zeros | FAIL | E | K^{-0.01} convergence | 5 |
| Explicit formula + 50 zeros | FAIL | I | More zeros can make WORSE | 6 |
| Explicit formula + 30 zeros | PARTIAL | E | 10% exact, best analytic | 10 |
| Weil explicit + Gaussian kernel | FAIL | I | Uncertainty principle: resolution * bandwidth >= const | 4 |
| Weil explicit + Beurling-Selberg | FAIL | I | Same uncertainty barrier | 4 |
| Mobius inversion of Pi(x) | PARTIAL | E | 10-30x better but still O(sqrt(x)/T) | 4 |
| Smooth zero sum integral | PARTIAL | E | Better than discrete but O(1) errors | 4 |
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
| Iwasawa theory | FAIL | I | lambda,mu invariants give no shortcut | 7 |
| Cloitre analytic recurrence (2025) | EXISTS | E | Exact but uses zeta evaluation, not competitive | lit |
| Farey fractions | FAIL | C | Require phi(k) -> factoring | 8 |
| Arithmetic derivative | FAIL | C | Requires factoring, O(n^{1/4}) | 8 |

## Sieve / Combinatorial / Counting

| Approach | Verdict | Mode | Key Finding | Session |
|----------|---------|------|-------------|---------|
| Sieve weights formula | FAIL | I | Parity barrier (Selberg) | 4 |
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
| TG kernel for exact pi(x) | FAIL | I | arXiv:2506.22634 DEBUNKED: violates uncertainty principle, omits x^{1/2} factor, kernel too wide to resolve primes. AI-generated paper with errors. | 12 |
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
| Novel intermediate quantity families | FAIL | E | Systematic analysis of 8 families (residues, polynomial evals, matrix eigenvalues, topology, representation theory, entropy, recursive, physical). ALL route back to floor values or zeta zeros. Only #TC^0 counting not immediately closed | 15 |
| Determinantal complexity of pi(x) | OPEN | - | pi(x) = degree-N multilinear polynomial in bits. Found N×N det reps for N=2,3,4. dc(pi_N)=poly(N) iff pi(x)∈GapL. For N>=10, generic polynomials DON'T have N×N det reps — pi(x) would need special structure | 15 |

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
| Function field analog (F_q) | FAIL | - | Genus finite -> O(g); for Z, genus = infinity | 7 |
| F_1 / Lattice / CFs | FAIL | I | q->1 degenerate, no structure | 10 |
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

---

## Summary Statistics

- **Total approaches tested:** 445+
- **FAIL:** ~385 (confirmed impossible or impractical)
- **PARTIAL:** ~16 (works partially but doesn't meet target)
- **WORKS:** ~8 (correct but O(x^{2/3}) or worse, or TC^0 but conditional)
- **EXISTS:** ~3 (formulas exist but computationally useless)
- **OPEN:** 5 (circuit complexity TC^0/NC, BPSW correctness, QFT deterministic, #TC^0⊆NC?, det complexity of pi(x))

### By Failure Mode
- **Circularity (C):** ~65 approaches
- **Equivalence (E):** ~190 approaches
- **Information Loss (I):** ~160 approaches
