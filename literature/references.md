# Consolidated Reference List

All references from the prime-research project, deduplicated and organized by category.
Each entry includes author(s), title, year, source, and relevance to the project.

Last updated: 2026-04-04

---

## 1. Algorithms for pi(x) and p(n)

**Legendre, A.-M.** "Essai sur la theorie des nombres." 1808.
Recursive inclusion-exclusion formula for pi(x), O(x^{3/4}). Foundation of all combinatorial pi(x) methods.

**Meissel, E.** Combinatorial decomposition for pi(x). 1870-1885.
Split contributions by number of prime factors, achieving O(x/(log x)^3). First sublinear pi(x) method.

**Lehmer, D. H.** Extension of Meissel's method with P3 term. 1959.
Improved to O(x/(log x)^4). Basis of all modern combinatorial pi(x) implementations.

**Lagarias, J.; Miller, V.; Odlyzko, A.** "Computing pi(x): The Meissel-Lehmer Method." 1985.
Achieved O(x^{2/3}/log x) via improved combinatorial identities.

**Lagarias, J.; Odlyzko, A.** "Computing pi(x): An Analytic Method." J. reine angew. Math. 387, 1987. [PDF](https://www-users.cse.umn.edu/~odlyzko/doc/arch/analytic.pi.of.x.pdf)
First O(x^{1/2+eps}) algorithm for pi(x) using Perron's formula and zeta zeros. Theoretical best upper bound; never implemented by its authors due to large constants.

**Deleglise, M.; Rivat, J.** "Computing pi(x): The Meissel, Lehmer, Lagarias, Miller, Odlyzko method." Math. Comp. 65, 1996. [PDF](https://cr.yp.to/bib/1996/deleglise.pdf)
Refined combinatorial method to O(x^{2/3}/log^2 x). Current best practical asymptotic for combinatorial pi(x).

**Gourdon, X.** Implementation of an improved Deleglise-Rivat variant. 2001.
O(x^{2/3}/log^2 x) with different space tradeoffs. Used in primecount library.

**Lucy_Hedgehog.** Iterative DP algorithm for pi(x). 2012 (Project Euler community). [Analysis](https://gbroxey.github.io/blog/2023/04/09/lucy-fenwick.html)
O(x^{2/3}) time, O(x^{1/2}) space, ~15 lines of code. Exploits that floor(x/k) takes O(sqrt(x)) distinct values. Widely used in competitive programming.

**Platt, D.; Trudgian, T.** "Computing pi(x) Analytically." 2012. [arXiv:1203.5712](https://arxiv.org/abs/1203.5712)
Rigorous interval-arithmetic implementation of analytic pi(x). Computed pi(10^25) unconditionally.

**Walisch, K.** primecount library. [GitHub](https://github.com/kimwalisch/primecount)
State-of-art implementation of Deleglise-Rivat and Gourdon algorithms. Parallelized, computes pi(10^25) in 330 CPU core hours. Versions 8.0-8.4 (2025-2026) added 128-bit API, SIMD Gourdon D formula, AVX512/ARM SVE support.

**Hirsch, D.; Kessler, I.; Mendlovic, U.** "Computing pi(N): An elementary approach in O~(sqrt(N)) time." arXiv:2212.09857, December 2022; Mathematics of Computation (published November 2024). [arXiv](https://arxiv.org/abs/2212.09857) [AMS](https://www.ams.org/journals/mcom/0000-000-00/S0025-5718-2024-04039-5/) [GitHub](https://github.com/PrimeCounting/PrimeCounting)
First elementary (non-analytic) algorithm achieving O(sqrt(N) * (log N)^{5/2}) for pi(N), matching Lagarias-Odlyzko asymptotics without complex analysis. Uses NTT-based Dirichlet convolutions extending the hyperbola method. Space-time tradeoff: O(N^{8/15}) time with O(N^{1/3}) space. Also improves Mertens function, totient sum, square-free counting. In practice ~400x slower than primecount at N=10^{14} due to constants. Peer-reviewed (Math. Comp.).

**Aggarwal, A.** "A Note on Algorithms for Computing p_n." arXiv:2510.16285, October 2025. [arXiv](https://arxiv.org/abs/2510.16285)
Proves p(n) computable in O(sqrt(n) * (log n)^4) via binary search on pi(x) using HKM algorithm. Conditional: O(sqrt(n) * (log n)^{7/2} * log log n) under RH + Cramer. Proves sieve-based methods cannot beat binary-search-on-pi(x) asymptotically: sublinear sieve is worse by factor ~sqrt(n).

---

## 2. Complexity Theory and Lower Bounds

**Agrawal, N.; Kayal, N.; Saxena, N.** "PRIMES is in P." Annals of Mathematics 160, 2002. [PDF](https://www.cse.iitk.ac.in/users/manindra/algebra/primality_v6.pdf)
AKS deterministic primality test in O~(log^6 n). Proved PRIMES is in P, but does not help with counting primes.

**Allender, E.; Saks, M.; Shparlinski, I.** "A Lower Bound for Primality." J. Comput. Sys. Sci. 62, 2001. [ScienceDirect](https://www.sciencedirect.com/science/article/pii/S0022000000917252)
PRIMES is NOT in AC^0[p] for any prime p. Best known unconditional circuit lower bound for primality.

**Hesse, W.; Allender, E.; Barrington, D.** "Division is in Uniform TC^0." J. Comput. Sys. Sci. 65, 2002. [ScienceDirect](https://www.sciencedirect.com/science/article/pii/S0022000002000259)
Division and iterated multiplication in uniform TC^0. Context for circuit complexity of arithmetic primitives.

**Agrawal, M.; Allender, E. et al.** "On TC0, AC0, and Arithmetic Circuits." [PDF](https://people.cs.rutgers.edu/~allender/papers/arithmetic.pdf)
Framework for understanding arithmetic circuit complexity. Relevant to whether pi(x) can be parallelized (NC question).

**Barrington, D.A.M.** "Bounded-width polynomial-size branching programs recognize exactly those languages in NC^1." JCSS 38(1), 1989. [ScienceDirect](https://www.sciencedirect.com/science/article/pii/0022000089900378)
NC^1 = width-5 polynomial-size branching programs. Foundation for matrix powering in NC^1.

**Barrington, D.A.M.; Therien, D.** "Finite monoids and the fine structure of NC^1." JACM 35(4), 1988. [ACM](https://dl.acm.org/doi/10.1145/48014.63138)
Algebraic characterization of TC^0: word problem over monoid M is in TC^0 iff every simple group dividing M is abelian (solvable).

**Immerman, N.; Landau, S.** "The complexity of iterated multiplication." Information and Computation 116(1), 1995. [PDF](https://people.cs.umass.edu/~immerman/pub/mult.pdf)
Iterated multiplication of n integers is in TC^0 via Chinese Remainder Representation.

**Mereghetti, C.; Palano, B.** "Threshold circuits for iterated matrix product and powering." RAIRO-ITA 34(1), 39-46, 2000. [NUMDAM](https://www.numdam.org/article/ITA_2000__34_1_39_0.pdf)
MPOW_k (M^n for fixed k x k matrix) IS in TC^0 for any fixed k. IMP_k (product of N given k x k matrices) NOT in TC^0 for k >= 3 unless TC^0 = NC^1. Critical for our matrix powering reduction.

**Healy, A.; Viola, E.** "Constant-Depth Circuits for Arithmetic in Finite Fields of Characteristic Two." STACS 2006. [PDF](https://www.ccs.neu.edu/home/viola/papers/FieldOps.pdf)
Exponentiation alpha^k in F_{2^n} is in TC^0 for specific irreducible polynomials. Shows polynomial-ring powering CAN be in TC^0 over F_2 even with growing dimension.

**Andrews, R.; Wigderson, A.** "Constant-Depth Arithmetic Circuits for Linear Algebra Problems." FOCS 2024. [arXiv](https://arxiv.org/abs/2404.10839)
Polynomial GCD, discriminant, resultant, Bezout coefficients in constant-depth arithmetic circuits. Does not include matrix powering.

**Krebs, A.; Lange, K.-J.; Reifferscheid, S.** "Characterizing TC^0 in Terms of Infinite Groups." STACS 2005. [Springer](https://link.springer.com/chapter/10.1007/978-3-540-31856-9_41)
Extended algebraic characterization of TC^0 to infinite groups via typed monoids.

**Gesmundo, F.; Ikenmeyer, C.** "Geometric complexity theory and matrix powering." Differential Geometry and its Applications 55, 2017. [arXiv](https://arxiv.org/abs/1611.00827) [ScienceDirect](https://www.sciencedirect.com/science/article/pii/S0926224517300876)
Replaced determinant with trace of variable matrix power in GCT. No-go result for orbit occurrence obstructions.

**Galby, E.; Ouaknine, J.; Worrell, J.** "On Matrix Powering in Low Dimensions." STACS 2015. [Dagstuhl](https://drops.dagstuhl.de/entities/document/10.4230/LIPIcs.STACS.2015.329)
Matrix powering positivity problem decidable in polynomial time for fixed dimensions 2 and 3 via Baker's theorem.

**Mauduit, C.; Rivat, J.** "Sur un probleme de Gelfond: la somme des chiffres des nombres premiers." Annals of Mathematics 171, 2010.
Proved primes are NOT k-automatic for any k. Rules out finite-automaton-based approaches.

**Raz, R.** "Natural Proofs for Super-Linear Lower Bounds on Linear Functions." ECCC TR26-008, January 2026. [ECCC](https://eccc.weizmann.ac.il/report/2026/008/)
First formalization of natural proof barriers for linear circuit complexity. Shows what approaches cannot work for proving super-linear lower bounds on linear functions over finite fields.

**Chen, L.; Tal, A.; Wang, Y.** "Super-Quadratic Lower Bounds for Depth-2 Linear Threshold Circuits." ECCC TR26-039, March 2026. [ECCC](https://eccc.weizmann.ac.il/report/2026/039/)
Function in E^NP requiring n^{2.5-eps}-size THR-of-THR circuits. Strongest known lower bound for depth-2 linear threshold model.

**Rossman, B.** "Riffle Rank." Preprint, May 2025. [PDF](https://users.cs.duke.edu/~br148/riffle-rank.pdf)
New tensor complexity measure on even-order tensors. Lower bounds on riffle rank of identity tensor would separate VNC^1 from VBP. Algebraic circuits only, not Boolean.

**Chukhin, Kulikov, Mihajlin, Smirnova.** "Conditional Complexity Hardness Results." ECCC TR25-038, April 2025. [ECCC](https://eccc.weizmann.ac.il/report/2025/038/)
Under NSETH assumptions, constructs explicit tensors with superlinear rank and monotone circuit lower bounds. Conditional, not unconditional.

**Kelley, Z.; Lovett, S.; Meka, R.** "Explicit Separations between Randomized and Deterministic Number-on-Forehead Communication." STOC 2024. [ACM](https://dl.acm.org/doi/pdf/10.1145/3618260.3649721)
First explicit 3-player function with O(1)-bit randomized but Omega(log^{1/3}(N)) deterministic NOF complexity. Most methodologically relevant recent advance for NOF lower bounds.

**Improved NOF Sifting.** arXiv:2505.01587, May 2025.
Improves KLM'24 lower bound from Omega(log^{1/3}) to Omega(log^{1/2}) via better grid-norm sifting. Shows weak one-sided pseudorandomness suffices for discrepancy bounds.

**Communication complexity from information causality.** arXiv:2602.10206, February 2026. [arXiv](https://arxiv.org/abs/2602.10206)
New information-theoretic framework for one-way communication complexity lower bounds. Potentially relevant to communication rank of prime counting (Session 17 result: rank(π_N) = 2^{N/2-1}+2).

---

## 3. Exact Formulas for Primes

**Mills, W. H.** "A Prime-Representing Function." Bull. Amer. Math. Soc. 53, 604, 1947.
floor(A^{3^n}) is prime for all n. Circular: constant A encodes the primes.

**Wright, E. M.** "A Prime-Representing Function." Amer. Math. Monthly 58, 616-618, 1951.
Iterated-exponential prime formula. Same circularity as Mills.

**Robinson, J.** Arithmetic terms for number-theoretic functions. 1952.
Showed binomial coefficients have arithmetic term representations: C(a,b) = floor((2^a+1)^a / 2^(ab)) mod 2^a. Foundation for Prunescu-Shunia.

**Willans, C. P.** "A Formula for the Nth Prime Number." Math. Gaz. 48, 413-415, 1964.
Wilson-based formula. O(4^n * n!) -- computationally useless beyond p(6), but proves a formula exists.

**Gandhi, J. M.** "Formulae for the Nth Prime." Proc. Washington State Univ. Conf. Number Theory, 96-106, 1971.
Primorial recurrence using Mobius function. Exact but exponential in primorial; infeasible beyond p(8).

**Jones, J. P.** "Formula for the Nth Prime Number." Canad. Math. Bull. 18, 433-434, 1975.
Compact Wilson variant using monus operator and Bertrand's postulate bound n^2.

**Jones, J. P.; Sato, D.; Wada, H.; Wiens, D.** "Diophantine Representation of the Set of Prime Numbers." Amer. Math. Monthly 83, 449-464, 1976.
Polynomial of degree 25 in 26 variables whose positive values are exactly the primes. Solutions involve numbers ~10^{10^52}; computationally impossible.

**Golomb, S. W.** Zeta-based recurrence for primes. 1976.
Requires knowing prior primes; non-constructive limit.

**Conway, J. H.** "FRACTRAN: A Simple Universal Programming Language for Arithmetic." Open Problems in Communication and Computation, Springer, 1987. [OEIS](https://oeis.org/wiki/Conway's_PRIMEGAME)
14-fraction program producing primes. Turing-complete but O(p^2) per prime.

**Mazzanti, S.** "Plain Bases for Classes of Primitive Recursive Functions." Math. Logic Quarterly 48, 2002.
Proved all Kalmar elementary functions have arithmetic term representations (non-constructive). Theoretical basis for Prunescu-Shunia.

**Rowland, E. S.** "A Natural Prime-Generating Recurrence." J. Integer Seq. 11, 2008.
GCD recurrence a(n) = a(n-1) + gcd(n, a(n-1)). Produces only primes but misses many; 18 distinct primes in 100,000 terms.

**Plouffe, S.** "A formula for primes." arXiv:1901.01849, 2018 (revised 2022). [PDF](http://plouffe.fr/articles/Formula%20for%20primes.pdf)
Exponential formulas {a_0 * r^n} with r=5/4 producing 50 probable primes. Conjectured, not proven. Constant-dependent.

**Fridman, D.; Garbulsky, J.; Glecer, B.; Grime, J.** "A Prime-Representing Constant." Amer. Math. Monthly 126, 70-73, 2019.
Elegant fractional expansion revealing successive primes. Circular: computing the constant requires knowing primes.

**Elsholtz, C.** "Unconditional Prime-Representing Functions." Bull. London Math. Soc. 52, 2020.
Mills-type formula proved without RH. Same constant-dependency problem.

**Prunescu, M.; Sauras-Altuzarra, L.** "An Arithmetic Term for the Factorial Function." J. Integer Seq., 2024. [ScienceDirect](https://www.sciencedirect.com/science/article/pii/S2666657X24000028)
First explicit arithmetic term for n!. Building block for the prime formula.

**Prunescu, M.; Shunia, J. M.** "On arithmetic terms expressing the prime-counting function and the n-th prime." arXiv:2412.14594, December 2024; J. Integer Seq. 28, Article 25.5.3, 2025. [arXiv](https://arxiv.org/abs/2412.14594) [JIS](https://cs.uwaterloo.ca/journals/JIS/VOL28/Prunescu/prunescu3.pdf)
First explicit fixed-length arithmetic terms for pi(n) and p(n) using only +, -, *, div, exp. Tower-of-exponentials intermediate values make evaluation impossible for any non-trivial input.

**Prunescu, M.; Shunia, J. M.** "Computational considerations on arithmetic terms for number-theoretic functions." J. Logic & Computation 35(3), April 2025. [Oxford](https://academic.oup.com/logcom/article-abstract/35/3/exaf012/8071363)
Arithmetic terms for tau(n), sigma(n), phi(n), modular inverse, integer root, discrete logarithm.

**Prunescu, M.; Shunia, J. M.** "Elementary closed-forms for non-trivial divisors." arXiv:2510.26939, October 2025.
Arithmetic terms for smallest/largest prime divisor. O(1) operations but O(2^n) time due to intermediate values.

**Prunescu, M.; Shunia, J. M.** "On the representation of number-theoretic functions by arithmetic terms." arXiv:2407.12928, July 2024. [arXiv](https://arxiv.org/abs/2407.12928)
Arithmetic terms for GCD function.

**Prunescu, M.** "Arithmetic closed forms count the Mersenne primes, the Fermat primes and the twin-prime pairs." arXiv:2512.01680, December 2025. [arXiv](https://arxiv.org/abs/2512.01680)
Extended arithmetic term framework to special prime classes using Lucas-Lehmer, Pepin, Clement tests.

**Semenov, S.** "Smooth Analytical Approximation of Prime Characteristic Function." arXiv:2504.14414, April 2025. [arXiv](https://arxiv.org/abs/2504.14414)
Smooth function P(n) via triple integral with periodic kernel: P(n) -> 1 for primes, P(n) < 1 for composites. Not a finite closed-form; evaluation cost at least as expensive as trial division.

---

## 4. Sieve Methods and Prime Distribution

**Maynard, J.; Tao, T.** Bounded gaps between primes. 2013-2014.
Infinitely many prime pairs with gap <= 246. Sieve-based (weighted Selberg sieve).

**Craig, W.; Ono, K.; van Ittersum, J.-W.** "Integer partitions detect the primes." PNAS, September 2024. arXiv:2405.06451. [arXiv](https://arxiv.org/abs/2405.06451) [PNAS](https://www.pnas.org/doi/10.1073/pnas.2409417121)
Primes are solutions to Diophantine equations in MacMahon partition functions. Infinitely many new prime characterizations. Structural, not computational.

**Lichtman, J. D.** "A modified linear sieve." Algebra & Number Theory 19(1), 2025. arXiv:2109.02851. [ANT](https://msp.org/ant/2025/19-1/ant-v19-n1-p01-p.pdf)
Equidistributes primes in arithmetic progressions to moduli up to x^{10/17}, surpassing the x^{7/12} barrier. Largest improvement to twin prime upper bound since 1986.

**Pascadi, A.** "Primes and smooth numbers to x^{5/8}." arXiv:2505.00653, May 2025. [arXiv](https://arxiv.org/abs/2505.00653)
Equidistribution to moduli x^{5/8} without Selberg eigenvalue conjecture. Uses triply-well-factorable weights.

**Bouras, M.** "A new primes-generating sequence." arXiv:2509.09745, 2025. [arXiv](https://arxiv.org/abs/2509.09745)
GCD-based recurrence producing only 1s and primes. Improvement over Rowland but still incomplete.

**Tao, T.; Teravainen, J.** "Quantitative Correlations on Prime Factors of Consecutive Integers." arXiv:2512.01739, December 2025. [arXiv](https://arxiv.org/abs/2512.01739)
Resolves multiple Erdos conjectures: omega(n+k) <= k infinitely often, irrationality of sum omega(n)/2^n, asymptotic for consecutive integers with equal prime divisor counts.

**Volfson, V.** "Generalization of Hardy-Littlewood Conjecture to Almost-Prime Tuples." arXiv:2603.13416, March 2026. [arXiv](https://arxiv.org/abs/2603.13416)
Extends Hardy-Littlewood k-tuple conjecture to almost-primes with Selberg constant correction factors.

**Hakobyan, T.** "Almost Prime Numbers." arXiv:2603.00679, February 2026. [arXiv](https://arxiv.org/abs/2603.00679)
New generalization of primality; shows composite almost-primes must be Carmichael numbers.

**Li, R.** "Prime-Producing Sieves and Distribution of alpha*p - beta mod 1." arXiv:2504.13195, April 2025. [arXiv](https://arxiv.org/abs/2504.13195)
Infinitely many primes p with ||alpha*p - beta|| < p^{-28/87}. Sharpens Jia (2000) using minimal Type-II information.

**Gensel, B.** "The Prime Number Formula of Gandhi." arXiv:1910.08362, 2019.
Exposition and simplification of Gandhi's primorial formula.

**Jakimczuk, R.** "Generalizations of the Gandhi Formula." 2024.
Extensions of the Gandhi primorial recurrence.

---

## 5. Analytic Number Theory and Zeta Function

**Cloitre, B.** "An effective analytic recurrence for prime numbers." arXiv:2508.02690, 2025. [arXiv](https://arxiv.org/abs/2508.02690)
First proven effective analytic recurrence: p(n+1) via zeta evaluation and Euler product over prior primes. Exact for all n but requires extreme-precision arithmetic.

**Saito, K.** "Mills' constant is irrational." arXiv:2404.19461, April 2024; Mathematika 71(3), 2025. [arXiv](https://arxiv.org/abs/2404.19461)
Resolved a long-standing open problem. Partial transcendence results obtained.

**Visser, M.** "The nth Prime Exponentially." Mathematics (MDPI) 13(11), 1844, 2025. [MDPI](https://www.mdpi.com/2227-7390/13/11/1844)
Exponentially tight explicit bounds on p(n) location. Useful for narrowing binary search.

**Visser, M.** "Effective exponential bounds on the prime gaps." International Mathematical Forum, 2025. [PDF](https://www.m-hikari.com/imf/imf-2025/1-4-2025/p/visserIMF1-4-2025.pdf)
Explicit bounds on g(n) = p(n+1) - p(n) with exponential factors and proven validity ranges.

**Visser, M.** "Behaviour of the sequence theta_n = theta(p_n)." arXiv:2507.14410, July 2025.
Analysis of Chebyshev's theta function at primes.

**Explicit formula for zeta zeros via Hermite polynomials.** arXiv:2312.00108; Ramanujan Journal, 2025. [Springer](https://link.springer.com/article/10.1007/s11139-025-01297-y)
Claims primes determine zero distribution and zeros can be computed without using zeta directly. Could simplify analytic pi(x) pipeline.

**Connes, A.** "The Riemann Hypothesis: Past, Present and a Letter Through Time." arXiv:2602.04022, February 2026. [arXiv](https://arxiv.org/abs/2602.04022)
42-page survey. Novel: approximates first 50 zeta zeros from primes < 13 to 10^{-55} accuracy using Weil quadratic form. Does not provide computational speedup for pi(x).

**Bailleul, A.; Hayani, M.; Untrau, T.** "Wasserstein Approach to Generalized Skewes Numbers." arXiv:2603.20093, March 2026. [arXiv](https://arxiv.org/abs/2603.20093)
Unconditionally disproves Fiorilli conjecture. Conditional upper bounds for generalized Skewes numbers using quantitative Kronecker-Weyl theorem in 1-Wasserstein metric.

**Vartziotis, D.** "Fourier Series from Additive Prime Factor Functions." arXiv:2602.13342, February 2026. [arXiv](https://arxiv.org/abs/2602.13342)
Sparse prime-indexed Fourier series from sopfr(n). Exact decomposition of summatory function B(n) ~ pi^2 x^2/(12 log x). See also companion: arXiv:2602.21270.

**Li, Z.** "Deterministic Fractal Set from Prime Number Sequence." arXiv:2603.00658, February 2026. [arXiv](https://arxiv.org/abs/2603.00658)
Cantor-like fractal from primes mod 16 with universal Hausdorff dimension 1/4, independent of specific prime residue sequence.

**Valley Scanner Algorithm.** arXiv:2512.09960, December 2025. [arXiv](https://arxiv.org/abs/2512.09960)
New geometry-driven method for locating zeta zeros without Gram points or bracketing. Uses real-to-complex conversion and tracks minima of |Z(t)|. Tested up to t ~ 10^20.

---

## 6. 2024-2026 Breakthroughs

**Guth, L.; Maynard, J.** "New large value estimates for Dirichlet polynomials." arXiv:2405.20552, May 2024; Annals of Mathematics, 2025. [arXiv](https://arxiv.org/abs/2405.20552) [Annals](https://annals.math.princeton.edu/articles/22049)
Improved Ingham zero-density bound from 3/5 to 13/25 -- first improvement in 80+ years. PNT in short intervals [x, x+x^{17/30}]. Most significant advance in analytic number theory in 2024.

**Gafni, A.; Tao, T.** "Exceptional intervals for the prime number theorem in short intervals." arXiv:2505.24017, May 2025. [arXiv](https://arxiv.org/abs/2505.24017)
Explicit quantitative relationship between zero density estimates and exceptional PNT failures. Uses Guth-Maynard estimates.

**Le Duc Hieu.** "Arithmetic progressions of primes in short intervals beyond 17/30." arXiv:2509.04883, September 2025. [arXiv](https://arxiv.org/abs/2509.04883)
Combines uniform short-interval PNT at theta > 17/30 with Green-Tao transference principle.

**Chen, L.; Tal, A.; Wang, Y.** "Superquadratic Lower Bounds for Depth-2 Linear Threshold Circuits." STOC 2026; ECCC TR26-039. [ECCC](https://eccc.weizmann.ac.il/report/2026/039/)
Breaks decade-old quadratic barrier for depth-2 THR-of-THR circuits: proves n^{2.5-eps} size lower bound. Closest result to TC^0 frontier but hard function is in E^NP, not number-theoretic.

**Raz, R.** "Natural Proofs for Super-Linear Lower Bounds for Linear Functions." ECCC TR26-008, Jan 2026. [ECCC](https://eccc.weizmann.ac.il/report/2026/008/)
Natural proofs barrier extends to linear functions under crypto assumptions. Relevant barrier for proving pi(x) hard.

**Andrews, R.; Garg, A.; Schost, E.** "Hilbert's Nullstellensatz is in the Counting Hierarchy." ECCC TR26-024 / arXiv:2602.17904, Feb 2026.
Constant-depth arithmetic circuits for multivariate resultant. Places Nullstellensatz in CH (previously only PSPACE).

**Aravind et al.** "Learning Read-Once Determinants and the Principal Minor Assignment Problem." ECCC TR26-035, 2026.
First poly-time algorithm for learning ROD polynomials; first NC algorithm for principal minor equivalence. Relevant to GapL/determinantal complexity question.

**Pain, J.-C.** "Three expressions of the n-th prime number." arXiv:2601.18816, Jan 2026.
Formal/exact expressions via Mobius, von Mangoldt, survival dynamics. NOT computational algorithms.

**Lamzouri, Y.** "An effective Linear Independence conjecture for the zeros of the Riemann zeta function and applications." arXiv:2311.04860, November 2023; revised May 2024. [arXiv](https://arxiv.org/abs/2311.04860)
Formulates effective version of the Linear Independence conjecture for zeta zero ordinates using heuristics from linear forms in logarithms (Lang-Waldschmidt). Under this conjecture: obtains best-possible Omega results for the error in PNT (matching Montgomery's conjecture), and conditionally resolves part of Gonek's conjecture on Mobius summatory function: M(x) = Omega_pm(sqrt(x) * (log log log x)^{5/4}). Relevant: strengthens understanding of oscillatory contributions that encode delta(n); if zeta zeros are linearly independent (over Q), no finite linear combination of zero contributions vanishes, blocking shortcuts via partial zero sums.

**Brandt, N.** "Lower Bounds for Levin-Kolmogorov Complexity." TCC 2024; ePrint 2024/687. [arXiv](https://eprint.iacr.org/2024/687)
Proves unconditionally that MKtP (decisional Levin-Kolmogorov complexity) cannot be decided in deterministic linear time. For RAM models with linear-time universal simulation, obtains quadratic lower bound: MKtP not in DTIME[O(n^2)]. Novel diagonalization technique that bypasses algebrization and natural proofs barriers. Relevant: MKtP is closely connected to prime detection complexity; stronger lower bounds via this approach could eventually constrain pi(x) computation, and the technique itself avoids the barriers that block proving primes hard.

**Elimelech, R.; David, O.; De la Cruz Mengual, C.; Kalisch, R.; Berndt, W.; Shalyt, M.; Silberstein, M.; Hadad, Y.; Kaminer, I.** "Algorithm-assisted discovery of an intrinsic order among mathematical constants." PNAS 121(25), 2024; arXiv:2308.11829. [PNAS](https://www.pnas.org/doi/abs/10.1073/pnas.2321440121) [arXiv](https://arxiv.org/abs/2308.11829)
Massively parallel algorithm discovering unprecedented number of continued fraction formulas for fundamental constants. Unveiled the Conservative Matrix Field (CMF) -- a novel structure relating constants through continued fractions analogous to conservative vector fields. CMFs provide a powerful framework for generating irrationality proofs. Relevant: the CMF framework extends PSLQ-type integer relation detection to structured families of continued fractions; could potentially be applied to search for identities involving zeta-zero-related constants, though no prime-counting connections found.

---

## 7. Quantum Computing

**Nayak, A.; Wu, F.** Quantum lower bound: Omega(sqrt(x)) queries for functions with prime structure.
Proves quantum speedup is at most quadratic over classical for counting-type problems.

**Berry-Keating conjecture.** Hypothetical Hamiltonian whose eigenvalues are zeta zeros.
Even if it exists, would need ~10^51 eigenvalues for p(10^100). No efficient quantum simulation known.

---

## 8. Information Theory and ML

**Prime Coding Theorem.** arXiv:2308.10817, 2023.
Information content of prime sequence exceeds any learnable pattern. Explains why ML cannot discover a prime formula.

**Prime classification with sparse encoding.** arXiv:2402.03363, 2024 (WITHDRAWN).
ResNet + Transformer for prime/non-prime classification. ~99% recall, never 100%. Withdrawn by authors.

**Kolpakov, A.; Rocke, A.** "Machine Learning of the Prime Distribution." arXiv:2403.12588, March 2024; PLOS ONE, September 2024. [arXiv](https://arxiv.org/abs/2403.12588) [PMC](https://pmc.ncbi.nlm.nih.gov/articles/PMC11432912/)
Uses maximum entropy methods to derive theorems in probabilistic number theory (Hardy-Ramanujan Theorem). Provides theoretical argument explaining Yang-Hui He's observations about learnability of primes. Argues Erdos-Kac law would very unlikely be discovered by current ML. Confirms ML learns the smooth part of prime distribution but NOT the oscillatory correction.

**AI in Number Theory: LLMs for Algorithm Generation.** arXiv:2504.19451, April 2025. [arXiv](https://arxiv.org/abs/2504.19451)
Qwen2.5-Math-7B achieves >= 0.95 accuracy on algorithmic number theory tasks with optimal hints. LightGBM classifier achieves >= 93.9% accuracy identifying Dirichlet character modulus from L-function zeros.

**Google DeepMind.** "AlphaEvolve: A Gemini-Powered Coding Agent for Designing Advanced Algorithms." arXiv:2506.13131, 2025. [arXiv](https://arxiv.org/abs/2506.13131)
Evolutionary algorithm discovery: 4x4 matrix multiplication in 48 scalar ops (breaking Strassen), improved solutions on 50+ open problems. Not yet applied to prime computation.

**Google DeepMind.** "Towards Autonomous Mathematics Research" (Aletheia). arXiv:2602.10177, February 2026. [arXiv](https://arxiv.org/abs/2602.10177)
AI research agent (Gemini Deep Think) with Generator/Verifier/Reviser subagents. Autonomous solutions to 4 open Erdos problems. Generated research paper on eigenweights in arithmetic geometry without human intervention. 95.1% on IMO-ProofBench. No application to prime counting algorithms.

**Gauthier, T.; Urban, J.** "Learning Conjecturing from Scratch." arXiv:2503.01389, March 2025. [arXiv](https://arxiv.org/abs/2503.01389)
Self-learning system for generating induction predicates on 16,197 OEIS-derived problems. Neural translator in feedback loop with Z3 prover. Solved 5,565 problems vs 2,265 by CVC5/Vampire/Z3 baselines. Relevant approach: automated conjecture generation from integer sequences could potentially discover prime-related identities, though no prime results reported.

**Beit-Halachmi, I.; Kaminer, I.** "The Ramanujan Library -- Automated Discovery on the Hypergraph of Integer Relations." arXiv:2412.12361, December 2024; ICLR 2025. [arXiv](https://arxiv.org/abs/2412.12361)
Hypergraph representation of mathematical constants (nodes) and formulas (edges). Uses PSLQ algorithm to discover 75 previously unknown connections between constants, including new formulas for pi, ln(2), and relations generalizing century-old Ramanujan pi-e formula. Open-source API. No prime-counting or zeta-related discoveries, but the PSLQ framework could be applied to search for identities involving pi(x) corrections.

**Shalyt, M.; Seligmann, U.; Beit-Halachmi, I.; David, O.; Elimelech, R.; Kaminer, I.** "Unsupervised Discovery of Formulas for Mathematical Constants." arXiv:2412.16818, December 2024; NeurIPS 2024. [arXiv](https://arxiv.org/abs/2412.16818)
Key innovation: convergence-dynamics-based distance metric for formula discovery (traditional numerical precision fails since near-true formulas give no signal). Tested on 1,768,900 Polynomial Continued Fraction formulas; discovered novel formulas for pi, ln(2), Gauss' constant, Lemniscate constant. Revealed patterns enabling generalization to infinite families. Paradigm of using convergence dynamics rather than endpoint values could potentially apply to discovering prime correction formulas.

**Nakasho, K.** "IntSeqBERT: Learning Arithmetic Structure in OEIS via Modulo-Spectrum Embeddings." arXiv:2603.05556, March 2026. [arXiv](https://arxiv.org/abs/2603.05556)
Dual-stream Transformer (91.5M params) encoding OEIS sequences via log-scale magnitude + sin/cos modulo embeddings (moduli 2-101), fused via FiLM. CRT-based solver converts probabilistic predictions to integers with 7.4x improvement over tokenized baselines (Top-1: 19.09% vs 2.59%). Strong negative correlation between information gain and Euler's totient ratio suggests composite moduli capture arithmetic structure more efficiently. Potentially applicable to predicting prime-related sequences modularly.

**Moreira, W.; Stubbs, J.** "Sequencelib: A Computational Platform for Formalizing the OEIS in Lean." arXiv:2601.11757, January 2026. [arXiv](https://arxiv.org/abs/2601.11757)
Formalized 25,000+ OEIS sequences in Lean 4, proved 1.6M theorems about their values. Includes SML-to-Lean transpiler. Infrastructure for automated conjecture verification on integer sequences.

---

## 9. Summatory Functions and Related

**Computation of Totient Summatory Function.** arXiv:2506.07386.
Related O(x^{2/3}) method for Euler totient sum. Same barrier as pi(x).

**Griff's blog on summing multiplicative functions.** [Blog](https://gbroxey.github.io/blog/2023/04/30/mult-sum-1.html)
Practical techniques for summatory functions using Dirichlet hyperbola method.

---

## 10. General References and Surveys

- [Prime-counting function (Wikipedia)](https://en.wikipedia.org/wiki/Prime-counting_function)
- [Formula for primes (Wikipedia)](https://en.wikipedia.org/wiki/Formula_for_primes)
- [Primality test (Wikipedia)](https://en.wikipedia.org/wiki/Primality_test)
- [NC complexity (Wikipedia)](https://en.wikipedia.org/wiki/NC_(complexity))
- [FRACTRAN (Wikipedia)](https://en.wikipedia.org/wiki/FRACTRAN)
- [Mills' constant (Wikipedia)](https://en.wikipedia.org/wiki/Mills'_constant)
- [Meissel-Lehmer algorithm (Wikipedia)](https://en.wikipedia.org/wiki/Meissel%E2%80%93Lehmer_algorithm)
- [Analytic Computation of pi(x) (Uni Bonn)](https://www.math.uni-bonn.de/people/jbuethe/topics/AnalyticPiX.html)
- [Wolfram MathWorld: Prime Formulas](https://mathworld.wolfram.com/PrimeFormulas.html)
- [Wolfram MathWorld: Prime Diophantine Equations](https://mathworld.wolfram.com/PrimeDiophantineEquations.html)
- [Willans' Formula (MathWorld)](https://mathworld.wolfram.com/WillansFormula.html)
- [FAQ: Is there a formula for the nth prime?](https://t5k.org/notes/faq/p_n.html)
- [Quanta: Sensational proof delivers new insights into prime numbers (2024)](https://www.quantamagazine.org/sensational-proof-delivers-new-insights-into-prime-numbers-20240715/)
- [Quanta: Mathematicians uncover a new way to count prime numbers (2024)](https://www.quantamagazine.org/mathematicians-uncover-a-new-way-to-count-prime-numbers-20241211/)
- [Scientific American: Top 10 Math Discoveries of 2025](https://www.scientificamerican.com/article/the-top-10-math-discoveries-of-2025/)
- [AMS Notices on partition primes (June 2025)](https://www.ams.org/journals/notices/202506/noti3198/noti3198.html)
- [Modern Breakthroughs in Prime Gaps (2025, Mathematics Magazine)](https://www.tandfonline.com/doi/full/10.1080/0025570X.2025.2481010)
- [Shunia's blog on prime formula](https://www.josephshunia.com/2024/08/10/an-easy-formula-for-prime-numbers/)
- [Shunia's GitHub](https://github.com/jshunia)

---

**Dawar, A.; Evans, A. T.** "Characterizing NC1 with Typed Monoids." arXiv:2508.11019, August 2025. FSTTCS 2025.
Extends Krebs et al.'s algebraic characterization of TC^0 via typed monoids to NC^1. Proves higher-dimensional monoid quantifiers collapse to unary. New algebraic tools at TC^0/NC^1 boundary but does NOT separate TC^0 from NC^1.

**Hu, J.; Manor, Y.; Oliveira, I.** "Failure of Symmetry of Information for Randomized Computations." ECCC TR26-021, 2026.
Proves unconditionally that Symmetry of Information fails for rKt complexity. First unconditional result for randomized time-bounded Kolmogorov complexity. Relevant to Kt complexity direction.

**Green, B.; Sawhney, M.** "New Way to Count Primes." 2024.
Proved Friedlander-Iwaniec conjecture: infinitely many primes of form p^2 + 4q^2 with p, q prime. Used Gowers norms. Structural result, does not help with computing pi(x).

**Bhattacharjee et al.** "Constant-depth circuits for polynomial GCD over any characteristic." arXiv:2506.23220, June 2025.
Extends Andrews-Wigderson to fields of any sufficiently large characteristic. Still over fields, not applicable to Z_n.

**Tao, T.; Trudgian, T.; Yang, A.** "New exponent pairs, zero density estimates, and zero additive energy estimates: a systematic approach." arXiv:2501.16779, January 2025.
Builds on Guth-Maynard. Creates ANTEDB database (teorth.github.io/expdb/). Obtains 4 new exponent pairs and new zero-density estimates via computational optimization. Significant infrastructure for analytic number theory but does not change pi(x) computational barrier.

**Chiari, M.** "Feasibility of Primality in Bounded Arithmetic." arXiv:2504.17041, April 2025/2026.
Proves AKS correctness in the bounded arithmetic theory VTC^0_2 (= first-order theory of TC^0 with smash function). Shows the PROOF of primality testing is formalizable in TC^0-level reasoning. Does not put PRIMES in TC^0 as a computational class.

**Semenov, S.** "A Smooth Analytical Approximation of the Prime Characteristic Function." arXiv:2504.14414, April 2025.
Constructs smooth P(n) via triple integral with periodic kernel. P(n)→1 for prime n, P(n)<1 for composite. Computationally O(n^2) or worse — no algorithmic value.

**Garg, A.; Oliveira, R.; Saxena, N.** "Primes via Zeros: Interactive Proofs for Testing Primality of Natural Classes of Ideals." arXiv:2503.20071, STOC 2025.
Reduces ideal primality testing to Sigma_3^p ∩ Pi_3^p (under GRH). About IDEALS in polynomial rings, not integer primes. Not relevant to pi(x).

---

**arXiv:2602.20917** (February 2026). Primes in arithmetic progressions to large moduli via refinements of Harman's sieve. Bilinear to x^{9/17}, trilinear to x^{17/32}. New sieve equidistribution bounds but does NOT change computational barrier for pi(x).

**Chiari (Jalali et al.)** "Feasibility of Primality in Bounded Arithmetic." arXiv:2504.17041, April 2025 (revised January 2026). Proves AKS correctness is formalizable in VTC^0_2. Does NOT place PRIMES in TC^0; shows the proof is TC^0-level reasoning.

**arXiv:2505.00730** (April 2025). Circulant matrix primality test. n is prime iff C_n has exactly two irreducible factors over Q. AI-assisted, no complexity analysis, unlikely competitive.

**ECCC TR26-038** (March 2026). Hardness Amplification for Counting over Finite Fields. Extends XOR lemma to sums over F_p. Not about number-theoretic counting.

**ECCC TR26-039** (March 2026). Super-Quadratic THR-of-THR Lower Bounds. Function in E^NP requiring n^{2.5-eps} depth-2 threshold circuits. Strongest known for this model.

---

**Craig, W.; van Ittersum, J.-W.; Ono, K.** "Integer partitions detect the primes." PNAS 121(39), 2024. arXiv:2405.06451.
Proves primes are solutions to infinitely many Diophantine equations in MacMahon partition functions. n>=2 prime iff (n²-3n+2)σ₁(n)-8M₂(n)=0. Runner-up 2025 Cozzarelli Prize. Beautiful but computationally worse than BPSW (requires divisor computation, O(n²) per test).

**Harper, A.; Xu, M. W.; Wang, V.; Soundararajan, K.** "Prime distribution via Gaussian multiplicative chaos and random fractal measures." Announced September 2025 conference.
Harper proved that statistics of zeta function zeros are captured by random fractal measures (Gaussian multiplicative chaos). Xu-Wang (2025) verified Harper's conjecture that GMC gives a better prime counting formula than Riemann's in the "beyond sqrt(x) barrier" transition interval [x, x+y]. In the transition regime, the exact mix of randomness and chaos is computable. The statistics revert to pure randomness for sufficiently small intervals. CAVEAT: This is about distributional statistics of primes, not exact computation. Does not provide a new algorithm for pi(x), but characterizes the information-theoretic structure of the oscillatory correction.

**Connes, A.; Consani, C.; Moscovici, H.** "Zeta zeros and prolate wave operators." arXiv:2310.18423, October 2023 (revised May 2024). [arXiv](https://arxiv.org/abs/2310.18423)
Introduces semilocal analogue of the prolate wave operator for spectral realization of zeta zeros. Positive spectrum captures low-lying zeros; negative spectrum (Sonin space) captures ultraviolet behavior. Operator is sum of squared scaling operator with grading of orthogonal polynomials. Extends to semilocal case via metaplectic representation of SL(2,R). Purely theoretical framework; no computational implications for prime counting. Relevant to Berry-Keating / Hilbert-Polya direction.

**Gafni, A.; Tao, T.** "Rough numbers between consecutive primes." arXiv:2508.06463, August 2025.
Almost all prime gaps contain a rough number with lpf ≥ gap length. Confirms Erdős conjecture (problem #682). Number of exceptions ≤ O(X/log² X). Gap structure, not counting.

**Dey, A.; Guo, Z.** "Debordering Closure Results in Determinantal and Pfaffian Ideals." arXiv:2511.16492, ECCC TR25-189, November 2025.
Deborders Andrews-Forbes STOC 2022 result: for f in determinantal ideal with poly degree, determinant is exactly (not approximately) computed by constant-depth f-oracle circuit. Advances algebraic complexity but doesn't apply to pi(x).

**Vijayaraghavan, T. C.** "The Complexity of Logarithmic Space Bounded Counting Classes." arXiv:2507.23563, December 2025.
Comprehensive textbook on GapL and related classes. Standard GapL theory. No new results on counting primes.

**Dwivedi, P.; Pago, B.; Seppelt, T.** "Lower Bounds in Algebraic Complexity via Symmetry and Homomorphism Polynomials." STOC 2026. [STOC](https://acm-stoc.org/stoc2026/accepted-papers.html)
New algebraic circuit lower bounds via symmetry. Relevant to determinantal/algebraic complexity of pi(x).

**Saraf, S.; Shringi, D.; Varadarajan, N.** "Reconstruction of Depth-3 Arithmetic Circuits with Constant Top Fan-in." STOC 2026.
Arithmetic circuit reconstruction. Relevant to understanding structure of bounded-depth arithmetic circuits.

**Grochow, J. A.** "Polynomial Identity Testing and the Ideal Proof System." STOC 2026.
Connects PIT to algebraic proof complexity. Foundational for algebraic circuit complexity.

**Alman, J.; Duan, R.; et al.** "More Asymmetry Yields Faster Matrix Multiplication." SODA 2025. arXiv:2404.16349. [arXiv](https://arxiv.org/abs/2404.16349)
New bound omega < 2.371339 (from 2.371552). Builds on Duan-Wu-Zhou FOCS 2023 and Vassilevska Williams-Xu-Xu-Zhou SODA 2024. Circumvents Schonhage tau theorem via asymmetric analysis. Relevant to any matrix-based approach to pi(x).

**Tao, T.; Trudgian, T.; Yang, A.** "New exponent pairs and zero density estimates via ANTEDB." arXiv:2501.16779, January 2025. [ANTEDB](https://teorth.github.io/expdb/)
Living database + Python optimizer for automated derivation of exponent pairs, zero density estimates, and zero additive energy estimates. Found 4 new exponent pairs and improved multiple zero density bounds through automated optimization. Most actionable new tool for systematically checking whether cascading known results yields better pi(x) complexity.

**Goldston, D.; Lee, Y.; Schettler, J.; Suriajaya, A.** "PCC implies 100% critical zeros." arXiv:2503.15449 (Part I, March 2025), arXiv:2507.06823 (Part II, July 2025), arXiv:2511.20059 (November 2025).
Using "horizontal multiplicity" technique, PCC alone (without RH) implies asymptotically 100% of zeta zeros are simple and on the critical line. First time PCC yields horizontal distribution. Simplifies conditional analytic pi(x) framework but doesn't help compute zeros.

**Hardy, O.; Xu, M. W.** "Helson's conjecture for smooth numbers." arXiv:2511.03430, November 2025.
Random multiplicative function sums over y-smooth numbers have beyond-sqrt cancellation, uniformly for (log x)^30 ≤ y ≤ x. Structural insight about smooth number sums; potentially relevant to smooth/rough decomposition in analytic pi(x) methods.

**Ferrari, M.; Hainry, E.; Pechoux, R.; Silva, R.** "Quantum Programming in Polylogarithmic Time." arXiv:2507.15415, July 2025. MFCS 2025.
First characterization of FBQPOLYLOG (functions computable in quantum polylog time). Proves FBQPOLYLOG strictly contained in QNC. Formal target class for impossibility results about quantum prime counting.

**Kuznetsov, A.** "Simple and accurate approximations to the Riemann zeta function." arXiv:2503.09519, March 2025.
Elementary function replacements for Riemann-Siegel remainder + precomputed Gaussian quadrature. Reduces constants in zeta evaluation, not asymptotics.

**Orellana Real, C.** "Valley scanner for zeta zeros." arXiv:2512.09960, December 2025.
Mountain-valley geometry approach to Hardy Z-function, replacing Gram-point bracketing. Cloud-based (AWS EC2), reaches height t ~ 10^20. New computational paradigm for individual zero finding.

Generated: 2026-04-05 (updated Session 29 internet search)
Sources: FINDINGS.md, all_known_exact_formulas.md, arithmetic_pi_formula.md, combinatorial_formulas.md, recursive_formulas.md, prunescu_shunia_formula.md, complexity_of_primes.md, latest_2026_breakthroughs.md, latest_research_2026.md, exact_formulas_2026.md, exact_formulas_research.md, april_2026_new_findings.md, session10_internet_search.md, session29_internet_search
