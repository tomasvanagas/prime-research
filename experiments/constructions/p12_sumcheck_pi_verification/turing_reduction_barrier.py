#!/usr/bin/env python3
"""turing_reduction_barrier.py -- a FINE-GRAINED (SETH) barrier to #P-hardness
of pi, closing the Turing-reduction lever S512 left open (open item 3).

THE QUESTION
------------
S512 (sharpP_probe.py) settled the *parsimonious* (many-one) #P-hardness route:
a reduction #A -> pi with pi(x(w)) = #A(w) must realize the target count c by
landing x in [p_c, p_{c+1}-1], so the reduction map c |-> x(c) IS the inverse
prime p(.) = the project goal -- CIRCULAR. S512 explicitly left open
"direction (i): a Turing (not parsimonious) reduction from a #P-complete count
to a pi-oracle that sidesteps the c|->x(c) circularity."

This cycle closes direction (i) CONDITIONALLY, with a fine-grained-complexity
argument that the parsimonious analysis cannot see: pi has its OWN sub-linear
algorithm, and that algorithm is itself the barrier.

THE BARRIER (conditional theorem)
---------------------------------
Let  alpha = inf { a : pi(x) is computable in time x^{a+o(1)} }.
  - Best PROVEN combinatorial:  a = 2/3   (Lagarias-Miller-Odlyzko 1985;
    Deleglise-Rivat 1996)            -> the implementable/measured rung.
  - Best known (analytic):      a = 1/2   (Lagarias-Odlyzko 1987, x^{1/2+eps})
    -> the strongest rung.
  - The PROJECT GOAL would be    a = 0    (polylog).

THEOREM.  Assume SETH (Impagliazzo-Paturi). Then there is NO poly(N)-time
Turing reduction from #SAT to a pi-FUNCTION-oracle in which every oracle query
has bit-length <= c*n with  c < 1/alpha  (n = number of #SAT variables).
Equivalently: any such reduction must query pi at some x >= 2^{(1/alpha)*n}.

PROOF.  Suppose M^pi computes #SAT in poly(n) time with all queries <= c*n bits,
c < 1/alpha. Pick a' with alpha < a' < 1/c. Replace each pi-oracle call by the
time-x^{a'} algorithm; each costs <= 2^{a'*c*n}. M makes poly(n) queries + poly(n)
extra work, so it computes #SAT(w) in time 2^{a'c*n + o(n)} with a'c < 1. Then
SAT (decide #SAT>0) is solved in 2^{(1-delta)n}, delta>0, for every k -- which
refutes SETH (a 2^{(1-delta)n} #SAT algorithm decides every k-SAT in that time).
QED.

So the critical query-blowup is  c* = 1/alpha:
    measured Lucy-Hedgehog a~=3/4 -> c* ~= 1.33   (this script MEASURES it)
    combinatorial          a =2/3 -> c* =  1.5     (proven algorithm)
    analytic               a =1/2 -> c* =  2       (best known)
A natural (parsimonious-ish) reduction queries pi NEAR the realizing prime
p_c ~ c*ln c, whose bit-length is ~ n + log n, i.e. blowup c -> 1 (MEASURED
below). Since 1 < c* for every sub-linear pi, ALL such reductions are
SETH-forbidden -- the circularity of S512 is now matched by a complexity wall
that also covers adaptive, non-parsimonious Turing reductions.

THE COROLLARY (the project tie-in)
----------------------------------
The same argument with alpha -> 0 gives: if pi has a polylog-time algorithm
(the PROJECT GOAL, a=0, c*=infinity), then NO poly-time Turing reduction of any
polynomial query-blowup makes pi #P-hard (under SETH) -- in fact a polylog pi
with pi #P-hard under poly-time Turing reductions would give P = #P outright.
=> "pi has a polylog algorithm" and "pi is #P-hard (poly-time Turing)" are
MUTUALLY EXCLUSIVE under P != NP. A #P-hardness PROOF for pi would, under P!=NP,
prove THIS PROJECT'S GOAL IMPOSSIBLE; the project's pursuit of the algorithm and
the search for hardness are the two horns of one dichotomy.

WHAT THIS IS / IS NOT
---------------------
This is a STANDALONE measurement + complexity anchor, not a protocol. The
meta-principle ("a sub-linear target + a tight reduction = a sub-linear source")
is folklore fine-grained complexity; the CONTRIBUTION is (1) applying it to pi
with pi's own sqrt(x) algorithm, (2) the quantitative threshold c* = 1/alpha and
the MEASURED parsimonious blowup c -> 1 < c*, (3) the polylog/#P-hardness
mutual-exclusivity corollary tying the barrier to the project goal. It is
CONDITIONAL (SETH) and covers sub-c* -blowup Turing reductions; super-c*-blowup
reductions (querying pi exponentially far from where the answer lives) remain
formally open but are not realized by any natural construction.

Run:
  python3 turing_reduction_barrier.py --selftest      # all checks
  python3 turing_reduction_barrier.py                 # full report
  python3 turing_reduction_barrier.py --pi-bench --xmax 8
  python3 turing_reduction_barrier.py --blowup-nmax 22
"""
import argparse
import math
import time

import numpy as np


# ===========================================================================
#  Number theory (self-contained; no chain imports)
# ===========================================================================
def sieve_bytearray(limit):
    """Boolean primality sieve [0..limit]."""
    if limit < 2:
        return bytearray([0]) * (limit + 1)
    s = bytearray([1]) * (limit + 1)
    s[0] = s[1] = 0
    for i in range(2, int(limit ** 0.5) + 1):
        if s[i]:
            s[i * i::i] = bytearray(len(range(i * i, limit + 1, i)))
    return s


def primes_upto(limit):
    s = sieve_bytearray(limit)
    return [i for i in range(2, limit + 1) if s[i]]


def pi_sieve(x):
    """pi(x) ground truth via sieve."""
    if x < 2:
        return 0
    s = sieve_bytearray(x)
    return int(sum(s))


def nth_prime(n):
    """The n-th prime p_n (1-indexed: p_1=2), exact via sieve with a PNT bound."""
    if n < 1:
        raise ValueError("n>=1")
    if n < 6:
        return [2, 3, 5, 7, 11][n - 1]
    # Rosser: p_n < n(ln n + ln ln n) for n>=6; pad for safety.
    ln = math.log(n)
    bound = int(n * (ln + math.log(ln))) + 10
    s = sieve_bytearray(bound)
    cnt = 0
    for i in range(2, bound + 1):
        if s[i]:
            cnt += 1
            if cnt == n:
                return i
    raise RuntimeError("nth_prime bound too small")  # pragma: no cover


def pi_lucy(x):
    """pi(x) via the Lucy_Hedgehog method: O(x^{3/4}) time, O(sqrt x) space.
    Computes S[v] = pi(v) for the O(sqrt x) distinct values v=floor(x/i). This is
    a genuinely SUB-LINEAR pi algorithm -- the object whose time exponent we
    measure to instantiate the barrier threshold c* = 1/alpha."""
    if x < 2:
        return 0
    r = int(math.isqrt(x))
    # the distinct values of floor(x/i): large block x//i for i=1..r, then the
    # small block 1..x//(r+1) ... assembled as the classic Lucy key list.
    V = [x // i for i in range(1, r + 1)]
    small_top = V[-1] - 1
    V.extend(range(small_top, 0, -1))
    S = {v: v - 1 for v in V}          # start: S[v] = v-1 (all of 2..v "prime")
    for p in range(2, r + 1):
        if S[p] > S[p - 1]:            # p prime  <=>  count jumped at p
            sp = S[p - 1]              # = pi(p-1)
            p2 = p * p
            for v in V:
                if v < p2:
                    break              # V is descending; rest are < p^2
                S[v] -= S[v // p] - sp
    return S[x]


# ===========================================================================
#  #SAT brute counter + the realizing-prime encoding (the worked example)
# ===========================================================================
def count_sat(clauses, n):
    """Exact #satisfying assignments of a CNF over variables 1..n.
    A clause is a list of nonzero ints: +i means x_i, -i means NOT x_i.
    Brute force over 2^n assignments (n small)."""
    cnt = 0
    for a in range(1 << n):
        ok = True
        for cl in clauses:
            sat = False
            for lit in cl:
                v = (a >> (abs(lit) - 1)) & 1
                if (lit > 0 and v == 1) or (lit < 0 and v == 0):
                    sat = True
                    break
            if not sat:
                ok = False
                break
        if ok:
            cnt += 1
    return cnt


def realizing_prime(c):
    """Smallest x with pi(x)=c is x=p_c (the c-th prime); window [p_c, p_{c+1}-1].
    Returns (p_c, p_{c+1}-1, bitlen(p_c))."""
    if c == 0:
        return (0, 1, 1)               # pi(x)=0 for x in [0,1]
    pc = nth_prime(c)
    pc1 = nth_prime(c + 1)
    return (pc, pc1 - 1, pc.bit_length())


# ===========================================================================
#  exponent fitting
# ===========================================================================
def fit_exponent(xs, ys):
    """Least-squares slope of log y vs log x (single power law)."""
    lx = np.log(np.asarray(xs, dtype=float))
    ly = np.log(np.asarray(ys, dtype=float))
    A = np.vstack([lx, np.ones_like(lx)]).T
    slope, intercept = np.linalg.lstsq(A, ly, rcond=None)[0]
    return float(slope), float(intercept)


# ===========================================================================
#  (1) MEASURE pi's sub-linear time exponent alpha
# ===========================================================================
def measure_pi_exponent(x_decades, repeats=1, verbose=True):
    """Time pi_lucy over x = 10^d for d in x_decades; fit the time exponent.
    Returns (rows, alpha_fit). alpha is the slope of log(time) vs log(x): pi_lucy
    runs in x^{3/4}, so alpha should approach ~0.75 as x grows."""
    rows = []
    for d in x_decades:
        x = 10 ** d
        best = math.inf
        res = None
        for _ in range(repeats):
            t0 = time.perf_counter()
            res = pi_lucy(x)
            dt = time.perf_counter() - t0
            best = min(best, dt)
        rows.append(dict(d=d, x=x, pi=res, t=best))
    # fit on the larger x (drop the smallest 1-2 points: Python/dict overhead
    # dominates there and flattens the slope).
    fit_rows = rows[1:] if len(rows) >= 4 else rows
    alpha, _ = fit_exponent([r['x'] for r in fit_rows], [r['t'] for r in fit_rows])
    if verbose:
        print("=" * 70)
        print("(1) MEASURED sub-linear time exponent of a real pi algorithm")
        print("    Lucy_Hedgehog pi(x): O(x^{3/4}) time, O(sqrt x) space")
        print("=" * 70)
        print(f"{'d':>3} {'x':>13} {'pi(x)':>12} {'time(s)':>10} "
              f"{'t/x^.75':>10}")
        for r in rows:
            print(f"{r['d']:>3} {r['x']:>13} {r['pi']:>12} {r['t']:>10.4f} "
                  f"{r['t'] / r['x'] ** 0.75:>10.2e}")
        print("-" * 70)
        print(f"  fitted time exponent alpha = {alpha:.3f}  "
              f"(theory 0.75 for Lucy; -> 2/3, 1/2 for LMO, analytic)")
        print(f"  => barrier threshold c* = 1/alpha = {1.0 / alpha:.3f}  "
              f"(measured); 1.5 (LMO), 2.0 (analytic)")
    return rows, alpha


# ===========================================================================
#  (2) the SETH / blowup accounting
# ===========================================================================
def barrier_threshold(alpha):
    """c* = 1/alpha: query-blowup below which a Turing reduction is SETH-forbidden."""
    return math.inf if alpha <= 0 else 1.0 / alpha


def seth_forbidden(alpha, c):
    """True iff a poly-time Turing reduction #SAT->pi with query-blowup c (max
    query bit-length = c*n) is forbidden under SETH: c < 1/alpha. (alpha=0 ->
    forbidden for every finite c.)"""
    return c < barrier_threshold(alpha)

def implied_sat_exponent(alpha, c):
    """A reduction of blowup c would solve #SAT in time 2^{(implied)*n}."""
    return alpha * c


def report_barrier(alphas, blowups, verbose=True):
    """For each (algorithm exponent alpha, claimed reduction blowup c): the
    implied #SAT time-exponent alpha*c and the SETH verdict (forbidden iff < 1)."""
    rows = []
    for name, a in alphas:
        cstar = barrier_threshold(a)
        for c in blowups:
            e = implied_sat_exponent(a, c)
            rows.append(dict(algo=name, alpha=a, cstar=cstar, c=c,
                             sat_exp=e, forbidden=seth_forbidden(a, c)))
    if verbose:
        print("=" * 70)
        print("(2) SETH BARRIER: a Turing reduction #SAT->pi with query-blowup c")
        print("    plugs pi's x^alpha algorithm in -> #SAT in 2^{alpha*c*n}.")
        print("    SETH forbids alpha*c < 1  <=>  c < c* = 1/alpha.")
        print("=" * 70)
        print(f"{'pi-algorithm':>20} {'alpha':>6} {'c*=1/a':>8} | "
              f"{'blowup c':>9} {'#SAT 2^(.)n':>12} {'verdict':>11}")
        last = None
        for r in rows:
            algo = r['algo'] if r['algo'] != last else ''
            last = r['algo']
            astr = f"{r['alpha']:.3f}" if algo else ''
            cs = ('  inf' if math.isinf(r['cstar']) else f"{r['cstar']:.3f}") \
                if algo else ''
            v = 'FORBIDDEN' if r['forbidden'] else 'allowed'
            print(f"{algo:>20} {astr:>6} {cs:>8} | {r['c']:>9.2f} "
                  f"{r['sat_exp']:>12.3f} {v:>11}")
    return rows


# ===========================================================================
#  (3) MEASURE the parsimonious blowup c -> 1
# ===========================================================================
def parsimonious_blowup_exact(n):
    """For count target C=2^n: the realizing prime p_C has bit-length B; the
    blowup of the parsimonious reduction (query x ~ p_C to encode 'answer=C') is
    c = B / n. EXACT via nth_prime (n small enough to sieve)."""
    C = 1 << n
    pc, _, B = realizing_prime(C)
    return dict(n=n, C=C, p_C=pc, bitlen=B, c=B / n)


def parsimonious_blowup_pnt(n):
    """The same blowup via PNT (p_C ~ C ln C), valid for any n: shows c -> 1.
    With C=2^n: log2(C ln C) = n + log2(ln C) = n + log2(n ln 2), so
    c = 1 + log2(n ln 2)/n  (computed symbolically -> no overflow for large n)."""
    lnC = n * math.log(2.0)          # ln(2^n)
    bitlen = n + math.log2(lnC)      # = log2(2^n * ln C)
    return dict(n=n, bitlen=bitlen, c=bitlen / n)


def report_parsimonious(exact_nmax=22, pnt_ns=(10, 20, 40, 80, 160, 320),
                        verbose=True):
    exact = [parsimonious_blowup_exact(n) for n in range(4, exact_nmax + 1, 2)]
    pnt = [parsimonious_blowup_pnt(n) for n in pnt_ns]
    if verbose:
        print("=" * 70)
        print("(3) the NATURAL (parsimonious) reduction has blowup c -> 1 < c*")
        print("    encode 'count = C=2^n' as pi(x)=C: smallest x is p_C, the")
        print("    C-th prime. blowup c = bitlen(p_C)/n.")
        print("=" * 70)
        print("  EXACT (x = p_C via sieve):")
        print(f"  {'n':>4} {'C=2^n':>12} {'p_C':>14} {'bitlen':>7} {'c=B/n':>7}")
        for r in exact:
            print(f"  {r['n']:>4} {r['C']:>12} {r['p_C']:>14} "
                  f"{r['bitlen']:>7} {r['c']:>7.3f}")
        print("  PNT extrapolation (p_C ~ C ln C), c = 1 + log2(ln C)/n -> 1:")
        print(f"  {'n':>4} {'bitlen(p_C)':>12} {'c':>7}")
        for r in pnt:
            print(f"  {r['n']:>4} {r['bitlen']:>12.1f} {r['c']:>7.3f}")
        print("-" * 70)
        print("  => parsimonious blowup c in (1, ~1.4] for these n, -> 1 as n->inf.")
        print("     1 < c* = 1/alpha for EVERY sub-linear pi (alpha<1), so the")
        print("     parsimonious reduction is SETH-forbidden -- independently of")
        print("     (and matching) the S512 circularity obstruction.")
    return exact, pnt


# ===========================================================================
#  (4) the worked #SAT -> pi encoding example
# ===========================================================================
def worked_example(verbose=True):
    """A concrete small 3-CNF: brute-count it, realize the count as pi(x)=C with
    x=p_C, and display the parsimonious-encoding blowup + the SETH verdict."""
    n = 12
    # a small satisfiable-but-nontrivial 3-CNF (count is some C in [1, 2^12]).
    clauses = [
        [1, -2, 3], [-1, 4, 5], [2, -3, 6], [-4, -5, 7], [3, 6, -8],
        [-6, 9, 10], [7, -8, 11], [-9, -10, 12], [8, 11, -12], [-7, 9, -11],
        [1, -6, 12], [-2, 5, -10], [4, -7, 8], [-3, 6, -9],
    ]
    C = count_sat(clauses, n)
    pc, hi, B = realizing_prime(C)
    c = B / n
    # cross-check: pi(p_C)=C, pi(p_C - 1)=C-1
    assert pi_sieve(pc) == C
    assert pi_sieve(pc - 1) == C - 1
    # the SETH-relevant quantity is the WORST-CASE query length over all
    # n-variable instances: some instance has count ~2^n (e.g. near-tautology),
    # forcing the largest query x ~ p_{2^n}.
    worst = parsimonious_blowup_exact(n)
    if verbose:
        print("=" * 70)
        print("(4) WORKED #SAT -> pi parsimonious encoding")
        print("=" * 70)
        print(f"  3-CNF: n={n} variables, {len(clauses)} clauses")
        print(f"  #SAT (brute over 2^{n}) = C = {C}")
        print(f"  realize as pi(x)=C: x in [p_C, p_(C+1)-1] = [{pc}, {hi}]")
        print(f"  smallest x = p_C = {pc}  (bitlen {B}); verified pi(p_C)=C, "
              f"pi(p_C-1)=C-1")
        print(f"  this instance blowup c = bitlen(x)/n = {B}/{n} = {c:.3f}")
        print(f"  WORST-CASE over n={n}-var instances (count ~ 2^{n}): "
              f"x ~ p_(2^{n}) = {worst['p_C']}, bitlen {worst['bitlen']}, "
              f"c = {worst['c']:.3f}")
        cw = worst['c']
        for nm, a in (('analytic a=1/2', 0.5), ('LMO a=2/3', 2 / 3)):
            cstar = barrier_threshold(a)
            verdict = 'FORBIDDEN' if seth_forbidden(a, cw) else 'allowed'
            print(f"    vs {nm}: c*={cstar:.3f} -> worst-case reduction "
                  f"{verdict} under SETH (#SAT in 2^{{{a * cw:.3f} n}})")
        print("  ALSO circular (S512): computing x=p_C from the formula is the")
        print("  inverse-prime p(.) = the project goal.")
    return dict(n=n, C=C, p_C=pc, bitlen=B, c=c, worst_c=worst['c'])


# ===========================================================================
#  (5) the mutual-exclusivity corollary
# ===========================================================================
def report_corollary(verbose=True):
    if verbose:
        print("=" * 70)
        print("(5) COROLLARY: polylog-pi  XOR  #P-hard-pi   (under P != NP)")
        print("=" * 70)
        print("  alpha -> 0 (the PROJECT GOAL, polylog pi) makes c* = 1/alpha = inf:")
        print("  a polylog pi-oracle answers every query in poly(N) time, so a")
        print("  poly-time Turing reduction from #SAT to pi would put #SAT in")
        print("  poly(n) time => P = #P (=> P = NP). Hence, assuming P != NP:")
        print("    * if pi has a polylog algorithm, pi is NOT #P-hard (poly-T);")
        print("    * a #P-hardness PROOF for pi (poly-time Turing) proves NO")
        print("      polylog algorithm exists -- i.e. proves THIS PROJECT")
        print("      IMPOSSIBLE.")
        print("  The two open questions are the two horns of one dichotomy; the")
        print("  sqrt(x) algorithm sits exactly between them (c* = 2).")
    return True


# ===========================================================================
#  full report
# ===========================================================================
def run_report(pi_xmax=7, blowup_nmax=20):
    decades = list(range(4, pi_xmax + 1))
    measure_pi_exponent(decades)
    print()
    alphas = [('measured Lucy 3/4', 0.75), ('combinatorial LMO 2/3', 2 / 3),
              ('analytic 1/2', 0.5), ('PROJECT GOAL polylog', 0.0)]
    report_barrier(alphas, blowups=[1.0, 1.33, 1.5, 2.0])
    print()
    report_parsimonious(exact_nmax=blowup_nmax)
    print()
    worked_example()
    print()
    report_corollary()
    print()
    print("BOTTOM LINE: pi's own sub-linear algorithm is a fine-grained barrier")
    print("to #P-hardness. Under SETH no poly-time Turing reduction #SAT->pi has")
    print("query-blowup < 1/alpha; the natural (parsimonious) blowup is c->1, so")
    print("every natural Turing reduction is forbidden -- extending S512's")
    print("parsimonious-circularity result to adaptive Turing reductions, and")
    print("(corollary) tying #P-hardness of pi to the impossibility of the goal.")


# ===========================================================================
#  selftest
# ===========================================================================
def selftest():
    rng = np.random.default_rng(0)
    n_ok = 0

    def check(cond, msg):
        nonlocal n_ok
        assert cond, "FAIL: " + msg
        n_ok += 1
        print(f"  ok  {msg}")

    print("-- number theory ground truth --")
    # 1-4: pi_lucy == pi_sieve over a range incl. boundaries
    for x in [0, 1, 2, 3, 10, 100, 1000, 10000, 100003, 222222]:
        check(pi_lucy(x) == pi_sieve(x), f"pi_lucy({x}) == pi_sieve = {pi_sieve(x)}")
    # known values
    check(pi_sieve(10) == 4, "pi(10)=4")
    check(pi_sieve(100) == 25, "pi(100)=25")
    check(pi_sieve(10 ** 6) == 78498, "pi(10^6)=78498")
    # 5: nth_prime exact small + larger
    check([nth_prime(i) for i in range(1, 11)] ==
          [2, 3, 5, 7, 11, 13, 17, 19, 23, 29], "nth_prime 1..10")
    check(nth_prime(168) == 997, "nth_prime(168)=997 (last prime < 1000)")
    check(nth_prime(1000) == 7919, "nth_prime(1000)=7919")
    # 6: realizing prime: pi(p_C)=C and pi(p_C-1)=C-1
    for C in [1, 10, 100, 500, 1229]:
        pc, hi, B = realizing_prime(C)
        check(pi_sieve(pc) == C and pi_sieve(pc - 1) == C - 1,
              f"realize C={C}: pi(p_C={pc})=C, pi(p_C-1)=C-1, window top={hi}")

    print("-- #SAT brute counter --")
    # 7: trivial formulas
    check(count_sat([], 3) == 8, "#SAT(empty, n=3) = 2^3 = 8")
    check(count_sat([[1]], 2) == 2, "#SAT(x1, n=2) = 2")
    check(count_sat([[1], [-1]], 1) == 0, "#SAT(x1 AND NOT x1) = 0")
    check(count_sat([[1, 2]], 2) == 3, "#SAT(x1 OR x2) = 3")
    # 8: brute count matches an independent enumeration on a random CNF
    n = 8
    cls = [[int(rng.integers(1, n + 1)) * (1 if rng.integers(2) else -1)
            for _ in range(3)] for _ in range(10)]
    brute = 0
    for a in range(1 << n):
        if all(any(((a >> (abs(l) - 1)) & 1) == (1 if l > 0 else 0) for l in cl)
               for cl in cls):
            brute += 1
    check(count_sat(cls, n) == brute, f"count_sat == independent enum ({brute})")

    print("-- exponent fitter (closed-form controls) --")
    # 9: fit recovers known power laws
    xs = [10 ** d for d in range(3, 9)]
    for true_a, ys in [(0.5, [x ** 0.5 for x in xs]),
                       (0.75, [x ** 0.75 for x in xs]),
                       (1.0, [x for x in xs]),
                       (2 / 3, [3 * x ** (2 / 3) for x in xs])]:
        a, _ = fit_exponent(xs, ys)
        check(abs(a - true_a) < 1e-6, f"fit_exponent recovers alpha={true_a:.3f}")

    print("-- barrier accounting --")
    # 10: c* = 1/alpha
    check(abs(barrier_threshold(0.5) - 2.0) < 1e-9, "c*(1/2) = 2")
    check(abs(barrier_threshold(2 / 3) - 1.5) < 1e-9, "c*(2/3) = 1.5")
    check(math.isinf(barrier_threshold(0.0)), "c*(0) = inf (polylog pi)")
    # 11: SETH verdict logic (the theorem's boundary)
    check(seth_forbidden(0.5, 1.0), "blowup 1.0 forbidden at alpha=1/2 (c*<2)")
    check(seth_forbidden(0.5, 1.99), "blowup 1.99 forbidden at alpha=1/2")
    check(not seth_forbidden(0.5, 2.01), "blowup 2.01 allowed at alpha=1/2")
    check(not seth_forbidden(0.5, 2.0), "blowup exactly 2.0 not (<2) -> allowed")
    check(seth_forbidden(0.0, 1e9), "polylog pi: every finite blowup forbidden")
    # 12: implied #SAT exponent and the SETH crossing alpha*c<1
    check(abs(implied_sat_exponent(0.5, 1.0) - 0.5) < 1e-9,
          "blowup 1, a=1/2 -> #SAT in 2^{0.5 n} (forbidden, <1)")
    check(implied_sat_exponent(0.5, 2.0) == 1.0,
          "blowup 2, a=1/2 -> #SAT in 2^{1.0 n} (boundary)")

    print("-- parsimonious blowup measurement --")
    # 13: exact blowup matches PNT trend and is in (1, c*) -> forbidden
    for n in (8, 12, 16, 20):
        ex = parsimonious_blowup_exact(n)
        pn = parsimonious_blowup_pnt(n)
        # exact uses integer bit_length (rounds up to <1 bit over log2); the
        # closeness tolerance must absorb that ~1.5/n quantization at small n.
        check(abs(ex['c'] - pn['c']) < 1.5 / n,
              f"n={n}: exact blowup {ex['c']:.3f} ~= PNT {pn['c']:.3f}")
        check(1.0 < ex['c'] < 2.0 and seth_forbidden(0.5, ex['c']),
              f"n={n}: blowup {ex['c']:.3f} in (1,2) -> SETH-forbidden")
    # 14: blowup is DECREASING toward 1 (c->1 as n->inf)
    cs = [parsimonious_blowup_pnt(n)['c'] for n in (10, 40, 160, 640, 2560)]
    check(all(cs[i] > cs[i + 1] for i in range(len(cs) - 1)) and cs[-1] < 1.05,
          f"PNT blowup decreasing to ~1: {[round(c,3) for c in cs]}")

    print("-- worked example & pi-exponent measurement --")
    # 15: worked example self-consistent (asserts pi(p_C)=C internally)
    we = worked_example(verbose=False)
    check(we['C'] >= 1 and we['c'] <= 1.5 and we['worst_c'] <= 1.5
          and seth_forbidden(0.5, we['c']) and seth_forbidden(0.5, we['worst_c']),
          f"worked example: C={we['C']}, blowup {we['c']:.3f}/worst "
          f"{we['worst_c']:.3f} (-> 1), both < c*=2 -> SETH-forbidden")
    # 16: measured pi exponent is sub-linear (alpha < 1 => c* > 1, the barrier
    #     exists for the project's OWN measured algorithm) and broadly ~0.75
    rows, alpha = measure_pi_exponent([4, 5, 6], verbose=False)
    check(0.4 < alpha < 1.0,
          f"measured Lucy time exponent alpha={alpha:.3f} in (0.4,1.0): "
          f"sub-linear => c*={1/alpha:.2f} > 1")
    check(all(r['pi'] == pi_sieve(r['x']) for r in rows),
          "pi-bench values all match sieve ground truth")

    print(f"\nALL {n_ok} CHECKS PASSED")
    return n_ok


# ===========================================================================
def main():
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--selftest', action='store_true', help='run all checks')
    ap.add_argument('--pi-bench', action='store_true',
                    help='only the pi time-exponent measurement')
    ap.add_argument('--xmax', type=int, default=7,
                    help='pi-bench: max decade d (x=10^d), default 7')
    ap.add_argument('--blowup-nmax', type=int, default=20,
                    help='max n (count=2^n) for the EXACT parsimonious blowup')
    args = ap.parse_args()

    if args.selftest:
        selftest()
        return
    if args.pi_bench:
        measure_pi_exponent(list(range(4, args.xmax + 1)), repeats=2)
        return
    run_report(pi_xmax=args.xmax, blowup_nmax=args.blowup_nmax)


if __name__ == '__main__':
    main()
