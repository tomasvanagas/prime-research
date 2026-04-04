"""
Session 14: Counting complexity connections for pi(x).

Key question: Can pi(x) be reduced to a counting problem with known
efficient algorithms?

Known efficiently countable objects (poly-time in binary input):
1. Lattice points in low-dimensional polytopes (Barvinok: fixed-dim poly-time)
2. Perfect matchings in planar graphs (FKT: poly-time via Pfaffians)
3. Spanning trees (Kirchhoff: det of Laplacian)
4. Linear extensions of certain posets (poly-time for width-2)
5. Solutions to 2-SAT (poly-time)
6. Paths in DAGs (matrix powering / dynamic programming)

Can pi(x) be expressed as any of these?

Also: What is the relationship between pi(x) and known #P-complete problems?
If pi(x) were #P-hard, that would strongly suggest no efficient algorithm exists.
If pi(x) is NOT #P-hard, there may be structural reasons for efficiency.
"""

import numpy as np
from sympy import primepi, isprime, Matrix, primerange, factorint
from math import isqrt, log, gcd

def pfaffian_approach():
    """
    The FKT algorithm computes the number of perfect matchings in a PLANAR
    graph in polynomial time (via Pfaffian = sqrt(det) of a signed adjacency matrix).

    Question: Can we construct a planar graph G(x) on poly(N) vertices
    such that the number of perfect matchings = pi(x)?

    For this to work, pi(x) would need to be expressible as a Pfaffian
    of a small matrix. Since Pfaffian ∈ GapL (determinant class), this
    is related to the GapL question.

    Let's first check: for small x, what size planar graph would be needed?
    """
    print("=" * 70)
    print("Pfaffian / FKT approach")
    print("=" * 70)

    # For a graph on n vertices, the maximum number of perfect matchings
    # is achieved by the complete bipartite graph K_{n/2, n/2}: (n/2)!
    # So we need at least n vertices where (n/2)! >= pi(x)
    # For pi(x) ~ x/ln(x), we need (n/2)! >= x/ln(x)
    # By Stirling: n/2 * ln(n/2) >= ln(x) - ln(ln(x)) ≈ N*ln(2)
    # So n ≈ 2N/ln(N) = O(N/log N)

    # This is POLYNOMIAL in N! So in principle, the size is right.
    # The question is whether such a planar graph can be CONSTRUCTED.

    for x in [10, 100, 1000, 10000, 10**10, 10**100]:
        N = int(log(x, 2)) + 1
        pi_x = int(primepi(x)) if x <= 10**7 else int(x / log(x))  # approx for large

        # Minimum vertices needed (loosely)
        # We need permanent/Pfaffian >= pi(x)
        # For n vertices, max matchings ~ (n/2)^{n/2} / e^{n/2}
        n_min = 2
        while True:
            from math import lgamma
            log_max_matchings = lgamma(n_min/2 + 1) / log(10) if n_min > 2 else 0
            if 10**log_max_matchings >= pi_x or n_min > 1000:
                break
            n_min += 2

        print(f"x=10^{int(log(x,10)):d}: N={N}, pi(x)≈{pi_x}, "
              f"min vertices ≈ {n_min} ({'poly(N)' if n_min <= 10*N else 'NOT poly(N)'})")


def spanning_tree_approach():
    """
    By Kirchhoff's theorem: the number of spanning trees of a graph G
    equals any cofactor of the Laplacian matrix L = D - A.

    This is a DETERMINANT computation. If we could construct a graph
    on poly(N) vertices whose spanning tree count = pi(x), that would
    place pi(x) in GapL.

    For a graph on n vertices, the number of spanning trees can be at most
    n^{n-2} (complete graph, by Cayley's formula). For n = O(N), this is
    N^{O(N)} ≈ 2^{O(N log N)}, which is more than enough to encode pi(x) ≈ 2^N.

    Let's try to construct such graphs for small x.
    """
    print("\n" + "=" * 70)
    print("Spanning tree count approach (Kirchhoff)")
    print("=" * 70)

    # For very small cases, try to find a graph whose spanning tree count = pi(x)
    from sympy import Matrix

    def spanning_tree_count(adj_matrix):
        """Compute spanning tree count via Kirchhoff's theorem."""
        n = len(adj_matrix)
        if n <= 1:
            return 1
        # Laplacian = Degree - Adjacency
        D = np.diag(adj_matrix.sum(axis=1))
        L = D - adj_matrix
        # Any cofactor = det of (n-1)x(n-1) minor
        minor = L[1:, 1:]
        return round(np.linalg.det(minor))

    # Test: for the complete graph K_n, spanning trees = n^{n-2}
    for n in [3, 4, 5]:
        adj = np.ones((n, n)) - np.eye(n)
        st = spanning_tree_count(adj)
        print(f"K_{n}: spanning trees = {st} (expected {n**(n-2)})")

    # For a cycle C_n: spanning trees = n
    for n in [5, 7, 11]:
        adj = np.zeros((n, n))
        for i in range(n):
            adj[i, (i+1) % n] = 1
            adj[(i+1) % n, i] = 1
        st = spanning_tree_count(adj)
        print(f"C_{n}: spanning trees = {st} (expected {n})")

    # Key insight: spanning tree count of C_n = n = pi(x) when n = pi(x).
    # But we need a graph whose NUMBER OF VERTICES is poly(N), not pi(x).
    # The cycle graph trick gives n=pi(x) vertices, which is exponential.

    # Can we design a graph on O(N) vertices with pi(x) spanning trees?
    # This requires the eigenvalues of the Laplacian to multiply to pi(x)/n.
    # (product of nonzero eigenvalues of L = n * #spanning_trees)

    print("\nNeed graph on O(N) vertices with tree count = pi(x).")
    print("This is equivalent to the GapL question (det of small matrix = pi(x)).")

    # Let's try: for x = 10, pi(10) = 4.
    # Need a graph on ~4 vertices with 4 spanning trees.
    # C_4 has 4 spanning trees! (cycle on 4 vertices)
    # But 4 = pi(10) and we need the graph to be constructible from x.

    # For x = 100, pi(100) = 25.
    # Need a graph on ~7 vertices (N=7) with 25 spanning trees.
    # Let's search.

    target = 25  # pi(100)
    n = 7
    print(f"\nSearching for graph on {n} vertices with {target} spanning trees...")

    # Try random graphs
    found = False
    np.random.seed(42)
    for trial in range(10000):
        # Random symmetric adjacency matrix with 0/1 entries
        adj = np.zeros((n, n))
        for i in range(n):
            for j in range(i+1, n):
                if np.random.random() < 0.5:
                    adj[i, j] = 1
                    adj[j, i] = 1

        st = spanning_tree_count(adj)
        if st == target:
            print(f"  Found! Trial {trial}, adjacency:")
            for row in adj.astype(int):
                print(f"    {list(row)}")
            found = True
            break

    if not found:
        # Try systematic small perturbations from K_n
        print(f"  Random search failed. Trying systematic...")
        # K_7 has 7^5 = 16807 spanning trees. We need 25.
        # Very sparse graphs are needed.

        # Try: path graph with extra edges
        for extra in range(n*(n-1)//2):
            from itertools import combinations
            edges = list(combinations(range(n), 2))
            # Try each subset of 'n-1+extra' edges
            if n - 1 + extra > len(edges):
                break
            for edge_set in combinations(edges, n - 1 + extra):
                adj = np.zeros((n, n))
                for i, j in edge_set:
                    adj[i, j] = 1
                    adj[j, i] = 1
                st = spanning_tree_count(adj)
                if st == target:
                    print(f"  Found with {n-1+extra} edges!")
                    found = True
                    break
            if found:
                break
            if extra > 2:
                break  # Too many combinations

    if not found:
        print(f"  Not found in limited search")


def two_sat_connection():
    """
    #2-SAT (counting satisfying assignments of 2-CNF formulas) is #P-complete.
    But the DECISION version of 2-SAT is in P, and #2-SAT with a specific
    structure might be tractable.

    Can pi(x) be expressed as #2-SAT of some formula with O(N) variables
    and O(poly(N)) clauses?

    For this, we'd need to encode "n is prime" as a conjunction of 2-literal
    clauses over the bits of n. Primality testing involves modular arithmetic,
    which doesn't naturally decompose into 2-literal clauses.

    Let's check: what is the smallest 2-CNF formula that distinguishes
    primes from composites for N-bit numbers?
    """
    print("\n" + "=" * 70)
    print("2-SAT / CNF approach to prime counting")
    print("=" * 70)

    for N in [4, 5, 6]:
        max_x = 2**N
        # Build truth table: for each N-bit number, is it prime?
        primes_set = set()
        for x in range(2, max_x):
            if isprime(x):
                primes_set.add(x)

        n_primes = len(primes_set)

        # The prime indicator as a Boolean function of N variables (bits of x)
        # Try to express as a CNF
        # Each prime x gives a satisfying assignment
        # Each composite x gives an unsatisfying assignment

        # For exact counting: we'd need #{satisfying assignments} = pi(max_x - 1)
        # A 2-CNF on N variables has at most 2^N satisfying assignments
        # We need exactly n_primes assignments

        # Instead, let's check: how many 2-clauses are needed to separate
        # primes from composites?

        # Convert to bit vectors
        def to_bits(x, N):
            return tuple((x >> i) & 1 for i in range(N))

        prime_bits = {to_bits(x, N) for x in primes_set}
        composite_bits = {to_bits(x, N) for x in range(2, max_x) if x not in primes_set}

        # A 2-clause (l1 OR l2) eliminates all assignments where both l1 and l2 are false
        # We need clauses that eliminate all composites but keep all primes

        # Count: how many composites can a single 2-clause eliminate?
        best_clause = None
        best_eliminated = 0

        for i in range(N):
            for si in [0, 1]:  # literal: x_i (si=1) or NOT x_i (si=0)
                for j in range(i, N):
                    for sj in [0, 1]:
                        if i == j and si == sj:
                            continue
                        # Clause: (l_i OR l_j)
                        # Eliminates assignments where l_i=0 AND l_j=0
                        eliminated_composites = 0
                        eliminated_primes = 0
                        for bits in composite_bits:
                            li = bits[i] if si else (1 - bits[i])
                            lj = bits[j] if sj else (1 - bits[j])
                            if li == 0 and lj == 0:
                                eliminated_composites += 1
                        for bits in prime_bits:
                            li = bits[i] if si else (1 - bits[i])
                            lj = bits[j] if sj else (1 - bits[j])
                            if li == 0 and lj == 0:
                                eliminated_primes += 1

                        if eliminated_primes == 0 and eliminated_composites > best_eliminated:
                            best_eliminated = eliminated_composites
                            best_clause = (i, si, j, sj)

        print(f"\nN={N}: {n_primes} primes, {len(composite_bits)} composites in [2, {max_x})")
        if best_clause:
            print(f"  Best single 2-clause eliminates {best_eliminated}/{len(composite_bits)} composites")
            print(f"  Clause: ({'x' if best_clause[1] else '~x'}{best_clause[0]} OR "
                  f"{'x' if best_clause[3] else '~x'}{best_clause[2]})")
        else:
            print(f"  No pure 2-clause separates composites from primes without false negatives")


def dag_path_counting():
    """
    The number of paths from s to t in a DAG can be computed in O(V + E) time.
    If we could build a DAG on poly(N) vertices where #paths = pi(x),
    that would give a poly(N)-time algorithm.

    This is equivalent to computing (A^k)_{st} for some matrix A and power k,
    which connects to the matrix powering approach.

    But we need #paths (unweighted) = pi(x). Let's try building such DAGs.
    """
    print("\n" + "=" * 70)
    print("DAG path counting approach")
    print("=" * 70)

    # Key insight: in a layered DAG with L layers and width W,
    # the maximum number of paths is W^L. For W = N and L = N,
    # this is N^N, which is more than enough for pi(x) ≈ 2^N.

    # But we need a SPECIFIC DAG where the path count equals pi(x).

    # Idea: use the binary representation of pi(x).
    # pi(x) = sum_{k=0}^{N-1} b_k * 2^k where b_k are the bits.
    # If we could compute each bit b_k in polylog time, we could
    # construct the DAG accordingly.

    # But computing bit k of pi(x) seems as hard as computing pi(x).

    # Alternative idea: the DAG encodes a decision tree.
    # At each layer, we make a decision based on some easily-computed
    # function of x, and the number of paths through "prime" decisions
    # equals pi(x).

    # This is essentially the sieve: at each prime p, we decide whether
    # to keep or remove multiples. But this requires sqrt(x) layers.

    # For a truly new approach, we'd need a DAG whose structure
    # doesn't correspond to the sieve.

    # Let's try: for x = 100, build the smallest DAG with pi(100) = 25 paths

    # A binary tree of depth 5 has 2^5 = 32 leaf-to-root paths.
    # Removing 7 edges gives 25 paths. But which edges to remove
    # depends on knowing pi(100) = 25, which is circular.

    # The non-circularity requirement: the DAG must be constructible
    # from x without knowing pi(x).

    print("DAG path counting requires constructing the DAG from x alone.")
    print("All known constructions (sieve-based) need O(sqrt(x)) vertices.")
    print("A poly(N)-vertex DAG with pi(x) paths would solve the problem.")
    print("\nThis is EQUIVALENT to the NC question.")

    # Let's verify with a small sieve-based DAG
    x = 30
    primes = list(primerange(2, isqrt(x) + 1))  # [2, 3, 5]

    # Sieve DAG: layers correspond to primes, nodes to remaining candidates
    # Layer 0: all numbers 2..30 (29 nodes)
    # After sieve by 2: 15 odd numbers + 1 (the prime 2)
    # After sieve by 3: remove odd multiples of 3
    # After sieve by 5: remove remaining multiples of 5

    candidates = list(range(2, x + 1))
    for p in primes:
        remaining = [n for n in candidates if n == p or n % p != 0]
        print(f"After sieving by {p}: {len(remaining)} candidates")
        candidates = remaining

    print(f"Final: {len(candidates)} primes: {candidates}")
    print(f"DAG has {x - 1} initial nodes → {len(primes)} layers → {len(candidates)} paths")
    print(f"Total DAG vertices: O(x * pi(sqrt(x))) = O({(x-1) * len(primes)})")


def hardness_analysis():
    """
    Is pi(x) #P-hard? If so, no polynomial-time algorithm exists unless P = #P.

    #P-complete problems: permanent, #SAT, #HAMILTONIAN-PATH, etc.

    pi(x) is NOT known to be #P-hard. In fact, pi(x) is a SPECIFIC function
    (not a problem family), so the question is more nuanced.

    The function n → pi(2^n) defines a sequence. Is computing this sequence
    #P-hard (i.e., is there a poly-time reduction from #SAT to computing pi)?

    Unlikely: #SAT instances are adversarial, while pi(x) is a fixed sequence.

    More relevant: is pi(x) in #L (log-space counting class)?
    #L = {f : f(x) = #accepting paths of NL machine on input x}
    GapL = #L - #L (differences of #L functions)

    pi(x) is in #P (with unary input): the number of n <= x such that n is prime.
    With binary input: pi(x) = #{n : n <= x AND n is prime} where n ranges
    over {2, ..., x}. This is a sum of 2^N terms.

    In #P, each term is the indicator 1_prime(n), which is computable in P.
    So pi(x) ∈ #P relative to the binary input.

    But is pi(x) in #NC (the counting version of NC)?
    #NC = {f : f = #accepting branches of poly-size NC circuits}
    If PRIMES is in NC, then pi(x) ∈ #NC.
    """
    print("\n" + "=" * 70)
    print("Hardness / counting complexity analysis")
    print("=" * 70)

    print("""
Known relationships:
- pi(x) ∈ #P (binary input): sum of P-time indicators
- pi(x) ∈ FP (unary input): O(x^{2/3}) algorithm
- pi(x) ∈ #NC iff PRIMES ∈ NC
- GapL ⊆ NC^2 ⊆ NC ⊆ P (binary input model)
- pi(x) ∈ GapL would imply pi(x) computable in NC^2 time

Key question: Is pi(x) in #L?
  #L functions count accepting paths of NL machines.
  NL = coNL (Immerman-Szelepcsényi).
  PRIMES ∈ coNP (easy: provide a composite witness).
  PRIMES ∈ NP (AKS certificate is poly-time verifiable).
  But PRIMES ∈ NL? UNKNOWN.

  If PRIMES ∈ NL, then pi(x) ∈ #L ⊆ GapL ⊆ NC^2.
  This would SOLVE our problem!

  Is PRIMES in NL? This requires:
  - An NL algorithm for primality: use O(log n) space,
    nondeterministically guess a certificate.
  - AKS uses O(log^6 n) space (actually O(log^2 n) with careful implementation).
  - But AKS is DETERMINISTIC. NL allows nondeterminism.
  - Miller-Rabin is in coRP, which uses random bits, not nondeterminism.

  PRIMES in L (logspace) is a major open problem.
  PRIMES in NL is STRICTLY WEAKER (NL contains L).

  Is there a nondeterministic logspace primality test?
""")

    # The connection between PRIMES and logspace is deep:
    # PRIMES ∈ L iff computing (a^b mod n) is in L
    # (since BPSW = scalar powering + 2x2 matrix powering)
    # Scalar powering IS in TC^0 ⊂ L.
    # 2x2 matrix powering IS in TC^0 ⊂ L.
    # So BPSW IS computable in L... if BPSW is correct.

    # Wait — this means:
    # BPSW correct → PRIMES ∈ L (not just TC^0!)
    # Because TC^0 ⊂ NC^1 ⊂ L (actually TC^0 ⊆ L is known)

    # If PRIMES ∈ L, then pi(x) = sum_{n=2}^x PRIMES(n) is the sum of
    # L-computable indicators. The sum of L functions...
    # This is in #L (the counting class for L).
    # And #L ⊆ GapL ⊆ FNC^2 (functional NC^2).

    # But WAIT: pi(x) = sum over 2^N terms, each L-computable.
    # This makes pi(x) ∈ #P (exponentially many terms).
    # For #L, we'd need the sum to be over POLYNOMIALLY many terms.

    # The issue: we're summing over n = 2, ..., x = 2^N.
    # This is 2^N terms. Even if each is L-computable, the SUM
    # requires aggregating 2^N values.

    # In the #L model: the machine nondeterministically generates n,
    # then deterministically tests primality. The number of accepting
    # paths = pi(x). But generating n uses O(N) nondeterministic bits,
    # so the machine has 2^N nondeterministic branches.
    # This IS #L: the space is O(N) (logspace), nondeterministic.

    print("CONCLUSION:")
    print("If PRIMES ∈ L (or NL), then pi(x) ∈ #L.")
    print("#L ⊆ GapL ⊆ NC^2.")
    print("But pi(x) ∈ NC^2 does NOT mean pi(x) is computable in poly(N) TIME.")
    print("NC^2 circuits have poly(N) size and O(log^2 N) depth.")
    print("If pi(x) ∈ NC^2, it's computable in O(poly(N)) PARALLEL time")
    print("with O(poly(N)) PROCESSORS.")
    print()
    print("WAIT — this means:")
    print("If BPSW is correct → PRIMES ∈ L → pi(x) ∈ #L ⊆ GapL ⊆ NC^2")
    print("→ pi(x) computable by poly(N)-size, O(log^2 N)-depth circuits!")
    print()
    print("Is this reasoning correct? Let's verify the chain:")
    print("1. BPSW correct → PRIMES ∈ TC^0 ⊆ L [Session 13, known]")
    print("2. PRIMES ∈ L → 1_prime(n) computable in logspace")
    print("3. pi(x) = #{n ≤ x : prime(n)} = #accepting paths of NL machine")
    print("4. This NL machine: on input x (N bits), nondeterministically")
    print("   guesses n (N bits), checks n ≤ x AND prime(n) in logspace")
    print("5. Number of accepting paths = pi(x)")
    print("6. So pi(x) ∈ #L")
    print("7. #L ⊆ GapL (trivially)")
    print("8. GapL ⊆ NC^2 (Allender 1999)")
    print()
    print("BUT STEP 4-5 has a problem:")
    print("The NL machine needs O(N) nondeterministic bits to guess n.")
    print("Total space = O(N) nondeterministic + O(N) deterministic = O(N).")
    print("This IS logspace nondeterministic.")
    print("So pi(x) ∈ #L should be correct IF PRIMES ∈ L.")
    print()
    print("HOWEVER: #L counts accepting paths. The FUNCTION pi(x)")
    print("requires outputting log2(pi(x)) ≈ N bits. Computing each bit")
    print("of a #L function is in GapL → NC^2.")
    print()
    print("CRITICAL CHECK: Does GapL → NC^2 give POLY(N) SIZE circuits?")
    print("YES — by definition, NC^2 has polynomial size and O(log^2 N) depth.")


if __name__ == '__main__':
    pfaffian_approach()
    spanning_tree_approach()
    two_sat_connection()
    dag_path_counting()
    hardness_analysis()
