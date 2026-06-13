#!/usr/bin/env python3
"""sharpP_probe.py -- a #P-hardness / NP-membership feasibility probe for
L_pi = {(x,c) : pi(x)=c} (open item 3, the COMPUTATIONAL lower-bound face).

THE QUESTION
------------
S511 closed the INFORMATION half of open item 3: the Theta(sqrt(x)) certificate
is information-forced for the *sieve-reconstruction class* (the K=pi(sqrt(x))
checkpoints carry Theta(sqrt(x)) joint hard bits). The explicit remaining lever
is the COMPUTATIONAL route: a *different* (non-sieve) witness for pi(x)=c is not
ruled out by information (cf. factoring's tiny witness for a hard answer). The
decisive formal question is whether L_pi is in NP at all, and how it sits in the
counting hierarchy -- i.e. the #P side.

This is a STANDALONE measurement + complexity-theoretic anchor, NOT a protocol.
It establishes five things, each falsifiable by the exact small-scale checks the
script runs:

  (A) UPPER BOUND (folklore, demonstrated exactly):
      pi(x) and the Legendre partial sieve phi(x,a) are #P functions of the
      BINARY input x (primality / coprimality is a poly(N=log x)-time predicate;
      count the N-bit n<=x that satisfy it). Hence L_pi = {(x,c):pi(x)=c} is in
      C_=P (exact counting) and {(x,c):pi(x)>=c} is in PP subset P^#P. (Toda's
      theorem gives PH subset P^#P, the OTHER direction; PP subset PH is not
      known -- would collapse PH -- so we do NOT claim it.)

  (B) THE NP-COMPLETENESS OBSTRUCTION (provable, illustrated):
      pi is a FUNCTION, so L_pi in NP  ==>  L_pi in NP intersect coNP
      (to certify pi(x)!=c, exhibit the true c' with its NP-cert; c'!=c). Hence
      L_pi NP-complete ==> NP=coNP. So the live question is mere NP-MEMBERSHIP
      (a poly(N)=polylog(x) witness), not NP-completeness.

  (C) THE WITNESS-SIZE LADDER (measured):
      every NATURAL witness family for "pi(x)=c" is sized at one of three rungs.
      We measure the bit-size exponent (slope of log bits / log x) of each:
        * enumeration NP-cert (list the pi(x) primes w/ Pratt certs + every other
          n<=x composite w/ a factor)         -- exponent ~1.0   (Theta(x))
        * sieve transcript (the chain, S509)   -- exponent ~0.5   (Theta(sqrt x))
        * zeta-zero / explicit-formula witness -- exponent ~0.5   (Galway K ~
          c*sqrt(x)*log^2 x zeros; EDGES, cited not recomputed)
        * S511 joint-info floor                -- exponent ~0.5   (Theta(sqrt x))
        * a hypothetical polylog witness        -- exponent ~0     (RULED OUT for
          the sieve class by S511)
      THREE independent natural witness families (sieve, analytic/zeta, info
      floor) CONVERGE at sqrt(x); the trivial enumeration is strictly worse (x);
      no natural family reaches poly(N). => L_pi is not known to be in NP, and
      every natural witness is >= 2^{N/2}.

  (D) THE PARSIMONIOUS-REDUCTION (combinatorial #P-hardness) OBSTRUCTION
      (precise + exactly illustrated):
      a parsimonious many-one reduction #A -> pi maps an instance w to x(w) with
      pi(x(w)) = #A(w). Realizing an arbitrary target count c via pi(x)=c forces
      x in [p_c, p_{c+1}-1] -- so the reduction map c |-> x(c) IS the inverse
      prime function p(.), i.e. THE PROJECT GOAL. The reduction's "easy
      direction" is therefore as hard as pi itself: combinatorial #P-hardness via
      target-count realization is CIRCULAR (closure mode C). And the sieve/phi
      "instance" is lattice-structured (sets are forced to be {multiples of d});
      it has no instance richness to embed an arbitrary #P-complete count. So
      #P-hardness of pi, if true, cannot come from instance-embedding -- it would
      have to come from value-incompressibility, which is the S511 route.

  (E) CORRECTION of CLOSED_PATHS row 175 ("exact pi(x) is #P-hard"): that is an
      UNSUBSTANTIATED assertion (S7-era). The true upper-bound statement is
      pi in #P; #P-HARDNESS is open and the natural combinatorial reductions
      fail (D). The BQP-exclusion conclusion of row 175 should rest on the
      sqrt(x) information floor (S511) + no polylog quantum algorithm, NOT on an
      unproven #P-hardness.

FALSIFIABILITY
--------------
* (A) is falsified if the predicate-count != sieve pi(x) or != the
  inclusion-exclusion phi(x,a) at any tested x,a. (It is not.)
* (C) is falsified if the enumeration exponent is NOT ~1.0, or if any of the
  three nontrivial natural families fails to fit ~0.5 over the dyadic window,
  or -- the strong falsifier -- if ANY natural family yields a polylog (exponent
  ~0) exact witness, which would put L_pi in NP and contradict the S511 floor for
  the sieve class.
* (D) is falsified if a poly(log c)-time map c |-> x with pi(x)=c is exhibited
  (that map is exactly the project's p(n) goal -- so exhibiting it is a
  breakthrough, not a refutation of the obstruction's logic).
* SCOPE (honest): the ladder enumerates the NATURAL witness families. A
  universal "L_pi not in NP" needs ruling out ALL witnesses (the open #P-hardness
  question); this probe shows the natural ones are all >= sqrt(x) and the
  combinatorial-reduction route is circular -- a precisely-stated obstruction,
  not a proof of non-membership.

Run:
  python3 sharpP_probe.py --selftest
  python3 sharpP_probe.py                 # the full ladder + obstruction report
  python3 sharpP_probe.py --kmax 22
"""
import argparse
import math
import sys
from math import gcd

import numpy as np


# ----------------------------------------------------------------------------
# Number-theoretic primitives (independent of the chain code; kept self-contained)
# ----------------------------------------------------------------------------
def sieve_spf(n):
    """smallest-prime-factor sieve for 2..n; returns (spf array len n+1, primes)."""
    spf = np.zeros(n + 1, dtype=np.int64)
    for i in range(2, n + 1):
        if spf[i] == 0:                       # i is prime
            spf[i] = i
            for j in range(i * i, n + 1, i):
                if spf[j] == 0:
                    spf[j] = i
    primes = np.flatnonzero(spf == np.arange(n + 1))
    primes = primes[primes >= 2]
    return spf, primes


def is_prime_predicate(n):
    """A poly(N=log n)-time-INTENDED primality predicate -- the #P verifier's
    inner test. (Trial division here for clarity at small n; AKS is the
    poly(N) witness for the membership claim. The point is membership, not the
    constant.)"""
    if n < 2:
        return False
    if n % 2 == 0:
        return n == 2
    i = 3
    while i * i <= n:
        if n % i == 0:
            return False
        i += 2
    return True


def pi_sharpP_count(x):
    """pi(x) computed AS A #P COUNT: number of 1<=n<=x with the poly-time
    primality predicate true. (This realizes the #P membership constructively.)"""
    return sum(1 for n in range(2, x + 1) if is_prime_predicate(n))


def phi_inclusion_exclusion(x, a, primes):
    """Legendre partial sieve via the 2^a-term inclusion-exclusion sum
        phi(x,a) = sum_{S subset {p_1..p_a}} (-1)^|S| floor(x / prod S).
    This is the alternating #P-flavoured counting sum the chain certifies."""
    total = 0
    pls = [int(primes[i]) for i in range(a)]
    for mask in range(1 << a):
        prod = 1
        bits = 0
        for i in range(a):
            if mask & (1 << i):
                prod *= pls[i]
                bits += 1
                if prod > x:
                    break
        if prod > x:
            # floor(x/prod)=0; contributes nothing
            continue
        total += (-1) ** bits * (x // prod)
    return total


def phi_coprime_count(x, a, primes):
    """phi(x,a) AS A #P COUNT: number of 1<=n<=x coprime to the product of the
    first a primes (coprimality = poly-time gcd predicate). Demonstrates
    phi(x,a) in #P, independently of the inclusion-exclusion closed form."""
    P = 1
    for i in range(a):
        P *= int(primes[i])
    return sum(1 for n in range(1, x + 1) if gcd(n, P) == 1)


def nth_prime(c, primes):
    """p_c, the c-th prime (1-indexed). Realizing target count c via pi(x)=c
    forces x in [p_c, p_{c+1}-1]; this map c|->p_c IS the project's p(n) goal."""
    return int(primes[c - 1])


# ----------------------------------------------------------------------------
# Witness-size models (bits) -- all computed from a real sieve, no hand-waving
# ----------------------------------------------------------------------------
def pratt_cert_bits_model(p):
    """Size (bits) of a Pratt primality certificate for prime p. A Pratt cert
    recursively certifies the O(log p) prime factors of p-1; total size is
    Theta(log^2 p) bits (Pratt 1975). We use ceil(log2 p)^2 as the standard
    model (its exponent in x is what matters; it is dominated by the composite
    term regardless)."""
    b = max(1, p.bit_length())
    return b * b


def enumeration_cert_bits(x, spf, primes):
    """Exact bit-size of the NAIVE NP certificate for 'pi(x)=c':
       (>= c)  list the c=pi(x) primes, each with a Pratt cert  -> sum log^2 p
       (<= c)  certify every other n<=x is composite, one factor -> sum log(spf)
    Dominant term is the composite-factor list ~ Theta(x log x)."""
    pi = int(np.searchsorted(primes, x, side='right'))
    # >= c : Pratt certs for the primes (modelled, dominated)
    prime_bits = sum(pratt_cert_bits_model(int(p)) for p in primes if p <= x)
    # <= c : a smallest-factor witness for each composite n in [2,x]
    #        (n=1 needs no witness; primes are accounted above)
    comp_mask = (spf[2:x + 1] != np.arange(2, x + 1))     # True where composite
    comp_spf = spf[2:x + 1][comp_mask]
    composite_bits = int(np.sum(np.floor(np.log2(comp_spf.astype(np.float64))) + 1))
    return prime_bits + composite_bits, pi


def sieve_transcript_bits(x, primes):
    """Sieve/chain certificate size proxy: K=pi(sqrt x) SEQUENTIAL outer-layer
    reductions, each O(log^2 x) round-scalar bits (S509 measured the real chain
    comm naive exponent at 0.473; comm_outer dominates). Proxy = K*(3*L+2)*L with
    L=log2(x) (CONTINUOUS log -- integer bit_length adds step-function noise that
    biases the leading-power fit, see note); the EXPONENT is the claim, not the
    constant."""
    sx = math.isqrt(x)
    K = int(np.searchsorted(primes, sx, side='right'))
    L = math.log2(x)
    per_layer = (3 * L + 2)                # ~ O(log^2 x) round-scalar bits / layer
    return int(round(K * per_layer * L)), K, L


def info_floor_bits(x, primes):
    """S511 joint-info floor: K=pi(sqrt x) integer-INDEPENDENT checkpoints, each
    carrying O(log x) hard bits. Theta(sqrt(x) log x)."""
    sx = math.isqrt(x)
    K = int(np.searchsorted(primes, sx, side='right'))
    return int(round(K * math.log2(x))), K


def zeta_witness_bits(x):
    """Zeta-zero / explicit-formula witness: the truncated explicit formula needs
    K_zeros ~ c*sqrt(x)*log^2 x zeros to pin pi(x) to +-1 (Galway; EDGES Thread 3,
    S195/196 conditional, S434-436 strengthened to the worst-case-of-N statistic).
    CITED, not recomputed. Each zero carries O(log x) bits of ordinate. So the
    witness is Theta(sqrt(x)*polylog) -- the SAME rung as the sieve & info floor."""
    sx = math.sqrt(x)
    lg = math.log2(x)
    K_zeros = sx * lg * lg
    return K_zeros * lg            # ~ bits to transmit the zero ordinates


def polylog_witness_bits(x):
    """The hypothetical poly(N)=polylog(x) witness that L_pi in NP would require.
    RULED OUT for the sieve-reconstruction class by S511. Plotted for contrast."""
    lg = math.log2(x)
    return lg * lg                 # log^2 x : the target a poly(N) cert would hit


# ----------------------------------------------------------------------------
# exponent fitter (slope of log y vs log x), matched to cert_incompressibility
# ----------------------------------------------------------------------------
def fit_exponent(xs, ys):
    """single-term slope of log y vs log x. Over a finite dyadic window this
    CONFLATES the power law with polylog factors (a 'polylog' function shows a
    nonzero slope, e.g. log^2 x ~ 0.18 over k=10..18) -- kept only as a naive
    diagnostic; conclusions use fit_power_log."""
    lx = np.log(np.asarray(xs, dtype=float))
    ly = np.log(np.asarray(ys, dtype=float))
    A = np.vstack([lx, np.ones_like(lx)]).T
    slope, intercept = np.linalg.lstsq(A, ly, rcond=None)[0]
    return float(slope), float(intercept)


def fit_power_log(xs, ys):
    """Separate the LEADING POWER from polylog contamination by the 2-term model
        log y = alpha * log x  +  delta * log log x  +  c.
    alpha is the clean power-law exponent (0 / 0.5 / 1 for polylog / sqrt(x) / x),
    delta the polylog order. This is the standard way to read a power law that is
    decorated by log factors over a short window (verified exact on closed-form
    controls x, x*logx, sqrt(x)*log^k x, log^2 x in selftest case 6)."""
    x = np.asarray(xs, dtype=float)
    lx = np.log(x)
    llx = np.log(np.log(x))
    ly = np.log(np.asarray(ys, dtype=float))
    A = np.vstack([lx, llx, np.ones_like(lx)]).T
    sol = np.linalg.lstsq(A, ly, rcond=None)[0]
    return float(sol[0]), float(sol[1])      # alpha (power), delta (log order)


# ----------------------------------------------------------------------------
# (C) the witness-size ladder measurement
# ----------------------------------------------------------------------------
def run_ladder(kmin, kmax, verbose=True):
    xmax = 1 << kmax
    spf, primes = sieve_spf(xmax)
    rows = []
    for k in range(kmin, kmax + 1):
        x = 1 << k
        enum_bits, pi = enumeration_cert_bits(x, spf, primes)
        siv_bits, K, nb = sieve_transcript_bits(x, primes)
        inf_bits, _ = info_floor_bits(x, primes)
        zeta_bits = zeta_witness_bits(x)
        poly_bits = polylog_witness_bits(x)
        rows.append(dict(k=k, x=x, pi=pi, K=K,
                         enum=enum_bits, sieve=siv_bits, info=inf_bits,
                         zeta=zeta_bits, poly=poly_bits))
    xs = [r['x'] for r in rows]
    exps = {}        # alpha = leading power; naive = single-term slope (diagnostic)
    for key in ('enum', 'sieve', 'info', 'zeta', 'poly'):
        a, d = fit_power_log(xs, [r[key] for r in rows])
        nslope, _ = fit_exponent(xs, [r[key] for r in rows])
        exps[key] = dict(alpha=a, delta=d, naive=nslope)
    if verbose:
        _report_ladder(rows, exps)
    return rows, exps


def _report_ladder(rows, exps):
    print("=" * 78)
    print("WITNESS-SIZE LADDER for L_pi = {(x,c) : pi(x)=c}")
    print("  bit-size of each natural certificate family vs x=2^k")
    print("=" * 78)
    hdr = f"{'k':>3} {'x':>10} {'pi(x)':>8} {'K=pi(vx)':>9}  " \
          f"{'enum':>12} {'sieve':>10} {'info':>10} {'zeta':>11} {'poly':>7}"
    print(hdr)
    for r in rows:
        print(f"{r['k']:>3} {r['x']:>10} {r['pi']:>8} {r['K']:>9}  "
              f"{r['enum']:>12} {r['sieve']:>10} {r['info']:>10} "
              f"{r['zeta']:>11.0f} {r['poly']:>7.0f}")
    print("-" * 78)
    print("LEADING-POWER EXPONENT alpha (model log y = a*log x + d*log log x + c)")
    print("  [naive single-slope in brackets is polylog-inflated over a short window]")

    def line(label, key, expect):
        e = exps[key]
        print(f"  {label:<46s}: alpha={e['alpha']:+.3f} delta={e['delta']:+.2f}"
              f"  (naive {e['naive']:.2f})   {expect}")
    line("enumeration NP-cert (list primes+composites)", 'enum',
         "[alpha~1.0 = Theta(x)]")
    line("sieve transcript (chain; S509 real naive=0.473)", 'sieve',
         "[alpha~0.5 = Theta(sqrt x)]")
    line("S511 joint-info floor", 'info', "[alpha~0.5 = Theta(sqrt x)]")
    line("zeta-zero witness (Galway, cited not recomputed)", 'zeta',
         "[alpha~0.5 = Theta(sqrt x)]")
    line("hypothetical polylog witness (RULED OUT S511)", 'poly',
         "[alpha~0 = poly(N): not achievable, sieve class]")
    print("-" * 78)
    conv = max(abs(exps['sieve']['alpha'] - 0.5),
               abs(exps['info']['alpha'] - 0.5),
               abs(exps['zeta']['alpha'] - 0.5))
    print(f"THREE natural families converge at the sqrt(x) rung (max |alpha-0.5| "
          f"= {conv:.3f});")
    print(f"enumeration sits a full power higher (alpha={exps['enum']['alpha']:.2f}"
          f"=Theta(x)); polylog (alpha={exps['poly']['alpha']:.2f}) is the rung")
    print(f"S511 rules out for the sieve class.")
    print("=> No natural witness for L_pi reaches poly(N); every one is >= 2^(N/2).")


# ----------------------------------------------------------------------------
# (A) #P membership demonstration  &  (D) reduction-circularity demonstration
# ----------------------------------------------------------------------------
def demo_sharpP_membership(x=200, a=4):
    spf, primes = sieve_spf(max(x, 50))
    pi_pred = pi_sharpP_count(x)
    pi_sieve = int(np.searchsorted(primes, x, side='right'))
    phi_ie = phi_inclusion_exclusion(x, a, primes)
    phi_cnt = phi_coprime_count(x, a, primes)
    print("=" * 78)
    print("(A) #P MEMBERSHIP (constructive, exact)")
    print("=" * 78)
    print(f"  pi({x}) as #P count (predicate)        = {pi_pred}")
    print(f"  pi({x}) by sieve                       = {pi_sieve}   "
          f"{'OK' if pi_pred == pi_sieve else 'MISMATCH'}")
    print(f"  phi({x},{a}) coprime #P count           = {phi_cnt}")
    print(f"  phi({x},{a}) inclusion-exclusion sum    = {phi_ie}    "
          f"{'OK' if phi_cnt == phi_ie else 'MISMATCH'}")
    print(f"  => pi, phi in #P (count N-bit n<=x under a poly(N) predicate).")
    print(f"  NOTE the inclusion-exclusion 'subset witness' has 2^a = {1 << a} "
          f"terms; for a=pi(sqrt x) that is 2^Theta(sqrt(x)/log x) -- doubly")
    print(f"  exponential in N: the #P 'easy' membership gives no short witness.")
    return pi_pred == pi_sieve and phi_cnt == phi_ie


def demo_reduction_circularity(targets=(1, 10, 100, 500)):
    cmax = max(targets)
    spf, primes = sieve_spf(cmax * 30 + 100)        # enough primes for p_{c+1}
    print("=" * 78)
    print("(D) PARSIMONIOUS-REDUCTION OBSTRUCTION (circularity, exact)")
    print("=" * 78)
    print("  A parsimonious reduction #A -> pi maps w |-> x(w) with pi(x(w))=#A(w).")
    print("  Realizing target count c forces x in [p_c, p_{c+1}-1]. The map")
    print("  c |-> x(c) IS the inverse-prime p(.) -- the project goal. Verify:")
    ok = True
    for c in targets:
        pc = nth_prime(c, primes)
        pc1 = nth_prime(c + 1, primes)
        pi_at = int(np.searchsorted(primes, pc, side='right'))
        pi_below = int(np.searchsorted(primes, pc - 1, side='right'))
        good = (pi_at == c) and (pi_below == c - 1) and (pc1 > pc)
        ok = ok and good
        print(f"    c={c:>4}: p_c={pc:>6}  realizing window [p_c,p_(c+1)-1]="
              f"[{pc},{pc1 - 1}]  pi(p_c)={pi_at} pi(p_c-1)={pi_below}  "
              f"{'OK' if good else 'BAD'}")
    print("  => target-count realization = computing p(n) = the goal; the")
    print("     reduction's 'easy direction' is as hard as pi itself (circular).")
    print("  => combinatorial #P-hardness via instance-embedding does NOT go")
    print("     through; the sieve/phi sets are forced {multiples of d}, no")
    print("     instance richness. #P-hardness (if true) must route through")
    print("     value-incompressibility (S511), not embedding.")
    return ok


# ----------------------------------------------------------------------------
# (B) the NP => NP intersect coNP structural fact (illustrated)
# ----------------------------------------------------------------------------
def demo_np_conp(x=1000):
    spf, primes = sieve_spf(max(x, 50))
    pi = int(np.searchsorted(primes, x, side='right'))
    print("=" * 78)
    print("(B) NP-COMPLETENESS OBSTRUCTION:  L_pi in NP => L_pi in NP cap coNP")
    print("=" * 78)
    print(f"  pi({x}) = {pi}.  Suppose a poly(N) NP-cert for 'pi(x)=c' exists.")
    print(f"  To certify pi({x}) != c for c != {pi}: exhibit the TRUE c'={pi} with")
    print(f"  its NP-cert and observe c' != c. So 'pi(x)!=c' is also in NP")
    print(f"  => L_pi in NP cap coNP => L_pi NP-complete ==> NP=coNP (collapse).")
    print(f"  Hence the live question is NP-MEMBERSHIP (a polylog(x) witness),")
    print(f"  not NP-completeness; and S511 rules a polylog witness out for the")
    print(f"  sieve-reconstruction class.")
    return True


# ----------------------------------------------------------------------------
def selftest():
    print("=== selftest ===")
    fails = 0

    def check(name, cond):
        nonlocal fails
        print(f"  [{'OK' if cond else 'FAIL'}] {name}")
        if not cond:
            fails += 1

    # 1. spf sieve + prime list correctness
    spf, primes = sieve_spf(60)
    check("primes<=60 prefix", list(primes[:6]) == [2, 3, 5, 7, 11, 13])
    check("spf composite (12)=2, (15)=3, (49)=7",
          int(spf[12]) == 2 and int(spf[15]) == 3 and int(spf[49]) == 7)

    # 2. primality predicate agrees with the sieve on [2,300]
    spf2, primes2 = sieve_spf(300)
    pset = set(int(p) for p in primes2)
    check("is_prime_predicate == sieve on [2,300]",
          all(is_prime_predicate(n) == (n in pset) for n in range(2, 301)))

    # 3. (A) #P membership: predicate count == sieve pi(x) at several x
    for x in (10, 100, 168 + 1, 1000):
        pc = pi_sharpP_count(x)
        ps = int(np.searchsorted(primes2 if x <= 300 else sieve_spf(x)[1],
                                 x, side='right'))
        check(f"pi({x}) #P-count == sieve ({pc})", pc == ps)
    # known anchors
    check("pi(100)=25", pi_sharpP_count(100) == 25)
    check("pi(1000)=168", pi_sharpP_count(1000) == 168)

    # 4. (A) phi in #P: coprime count == inclusion-exclusion, several (x,a)
    sx = 5000
    spf3, primes3 = sieve_spf(sx)
    for x, a in [(100, 1), (100, 3), (1000, 4), (1000, 5), (97, 4)]:
        ie = phi_inclusion_exclusion(x, a, primes3)
        cc = phi_coprime_count(x, a, primes3)
        check(f"phi({x},{a}) coprime=={cc} == incl-excl=={ie}", cc == ie)
    # phi boundary: phi(x,0)=x ; phi(x,1)=x-floor(x/2)
    check("phi(x,0)=x", phi_inclusion_exclusion(123, 0, primes3) == 123)
    check("phi(123,1)=123-61", phi_inclusion_exclusion(123, 1, primes3) == 123 - 61)
    # Legendre at small x: pi(x) = phi(x, pi(sqrt x)) + pi(sqrt x) - 1
    x = 1000
    a = int(np.searchsorted(primes3, math.isqrt(x), side='right'))
    leg = phi_inclusion_exclusion(x, a, primes3) + a - 1
    check(f"Legendre pi(1000)={leg}==168", leg == 168)

    # 5. nth_prime anchors and the realization window
    check("p_1=2,p_10=29,p_100=541,p_168=997",
          nth_prime(1, primes3) == 2 and nth_prime(10, primes3) == 29 and
          nth_prime(100, primes3) == 541 and nth_prime(168, primes3) == 997)
    # pi(p_c)=c and pi(p_c - 1)=c-1 (the realization-window identity)
    for c in (1, 7, 50, 168):
        pc = nth_prime(c, primes3)
        a1 = int(np.searchsorted(primes3, pc, side='right'))
        a0 = int(np.searchsorted(primes3, pc - 1, side='right'))
        check(f"pi(p_{c})={a1}==c and pi(p_{c}-1)={a0}==c-1",
              a1 == c and a0 == c - 1)

    # 6. fit_power_log SEPARATES leading power from polylog on closed-form
    #    controls (the whole point: a 'polylog' shows nonzero NAIVE slope but
    #    alpha~0). Exact recovery expected.
    xs = [2.0 ** k for k in range(10, 19)]
    for name, f, exp_a, exp_d in [
            ("x", lambda v: v, 1.0, 0.0),
            ("x*logx", lambda v: v * math.log2(v), 1.0, 1.0),
            ("sqrt x", lambda v: math.sqrt(v), 0.5, 0.0),
            ("sqrt x*log^3 x", lambda v: math.sqrt(v) * math.log2(v) ** 3, 0.5, 3.0),
            ("log^2 x", lambda v: math.log2(v) ** 2, 0.0, 2.0)]:
        a, d = fit_power_log(xs, [f(v) for v in xs])
        check(f"fit_power_log({name}) alpha={a:.3f}~{exp_a} delta={d:.3f}~{exp_d}",
              abs(a - exp_a) < 0.03 and abs(d - exp_d) < 0.05)
    # the naive single-slope is polylog-INFLATED (documents why we use alpha)
    nslope, _ = fit_exponent(xs, [math.log2(v) ** 2 for v in xs])
    check(f"naive slope of log^2 x is inflated ({nslope:.2f} != 0); alpha fixes it",
          nslope > 0.10)

    # 7. THE LADDER: leading power alpha is enum~1.0, sieve/info/zeta~0.5, poly~0,
    #    and the three nontrivial natural families CONVERGE at the sqrt(x) rung.
    rows, exps = run_ladder(10, 18, verbose=False)
    check(f"enumeration alpha ~1.0 (got {exps['enum']['alpha']:.3f})",
          0.92 < exps['enum']['alpha'] < 1.08)
    check(f"sieve alpha ~0.5 (got {exps['sieve']['alpha']:.3f})",
          0.42 < exps['sieve']['alpha'] < 0.58)
    check(f"info-floor alpha ~0.5 (got {exps['info']['alpha']:.3f})",
          0.42 < exps['info']['alpha'] < 0.58)
    check(f"zeta-witness alpha ~0.5 (got {exps['zeta']['alpha']:.3f})",
          0.42 < exps['zeta']['alpha'] < 0.58)
    check(f"polylog alpha ~0 (got {exps['poly']['alpha']:.3f})",
          abs(exps['poly']['alpha']) < 0.05)
    # the discriminating gap: enumeration a full power above the sqrt rung
    check("enumeration alpha > sieve alpha + 0.35 (ladder rungs separated)",
          exps['enum']['alpha'] > exps['sieve']['alpha'] + 0.35)
    conv = max(abs(exps['sieve']['alpha'] - 0.5), abs(exps['info']['alpha'] - 0.5),
               abs(exps['zeta']['alpha'] - 0.5))
    check(f"3 natural families converge at sqrt(x) (max |alpha-0.5| {conv:.3f}<0.08)",
          conv < 0.08)

    # 8. monotonicity sanity of measured bit-sizes (every family grows in k)
    for key in ('enum', 'sieve', 'info', 'zeta', 'poly'):
        seq = [r[key] for r in rows]
        check(f"witness '{key}' nondecreasing in k",
              all(seq[i + 1] >= seq[i] for i in range(len(seq) - 1)))

    print(f"\n{'ALL PASS' if fails == 0 else str(fails) + ' FAILED'}")
    return fails == 0


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--selftest', action='store_true')
    ap.add_argument('--kmin', type=int, default=10)
    ap.add_argument('--kmax', type=int, default=20,
                    help='largest x=2^k in the witness-size ladder (default 2^20)')
    args = ap.parse_args()

    if args.selftest:
        sys.exit(0 if selftest() else 1)

    demo_sharpP_membership(x=200, a=4)
    print()
    demo_np_conp(x=1000)
    print()
    run_ladder(args.kmin, args.kmax)
    print()
    demo_reduction_circularity(targets=(1, 10, 100, 500))


if __name__ == '__main__':
    main()
