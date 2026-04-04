# Closed Paths: Master Lookup (380+ Approaches)

**SEARCH THIS FILE before proposing any approach.** Use grep/ctrl-F.

Last updated: 2026-04-04 (Sessions 1-10, 89+ sub-agents)

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
| Mertens function shortcuts | FAIL | E | O(x^{2/3}) -- same as pi(x) | 6 |
| Dirichlet series / prime zeta | FAIL | E | Perron integral = explicit formula | 6 |
| Smooth number elimination | FAIL | I | 12% reduction, insufficient | 5 |
| Bit-by-bit construction | WORKS | E | log2(p(n)) pi(x) calls, each O(x^{2/3}) | 5 |
| Recursive pi(x) via pi(x/2) | FAIL | I | Error O(sqrt(x)/ln^2(x)), too large | 7 |
| Recursive identity (8 methods) | FAIL | E | All reduce to Meissel-Lehmer or explicit formula | 9 |
| GCD/coprimality matrix | FAIL | I | Eigenvalues != primes | 8 |
| Newton's identities (power sums) | WORKS | E | O(x^{2/3}) per sum, unstable deg>50 | 7 |

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
| Compressed pi(x) | PARTIAL | I | 5 bits/prime irreducible entropy | 8 |

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
| Internet search 2024-2026 | - | - | No new algorithms found | 5,8,9,10 |

---

## Summary Statistics

- **Total approaches tested:** 380+
- **FAIL:** ~350 (confirmed impossible or impractical)
- **PARTIAL:** ~15 (works partially but doesn't meet target)
- **WORKS:** ~5 (correct but O(x^{2/3}) or worse)
- **EXISTS:** ~3 (formulas exist but computationally useless)
- **OPEN:** 1 (circuit complexity -- TC^0/NC^1 status unknown)

### By Failure Mode
- **Circularity (C):** ~60 approaches
- **Equivalence (E):** ~160 approaches
- **Information Loss (I):** ~150 approaches
