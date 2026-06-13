#!/usr/bin/env python3
"""
compressed_prover_mult_trace.py — S492 (step 1) + S498 (field lift).

The exact-pi(x) interactive proof (lucy_dp_verification.py +
lucy_dp_delegated_wiring.py) has an Õ(sqrt x) VERIFIER but an Õ(x^{3/2})
table-based demo PROVER, because the demo materializes every DP layer over
the full 2^n value space. A production prover must work the O(sqrt x)-size
compressed Lucy state: the small side (values v <= sqrt x, value-addressed)
and the large side (values floor(x/d), d-addressed by the index d).

The large-side layer update needs S_{i-1}(floor(floor(x/d)/p)) =
S_{i-1}(floor(x/(d p))). When d p <= sqrt x this is large-side index d p
(multiply the INDEX by the constant p -> affine wiring, existing
machinery). When d p > sqrt x the target floor(x/(d p)) is a SMALL value
and must be looked up by value in the small side -- which requires the
verifier to be convinced of floor(x/(d p)) WITHOUT computing it. That is a
division of the public constant x by the prover-variable (d p): the
certifying identity

        u_d * (d p) + r_d = x ,   0 <= r_d < d p

has the product u_d * d of two quantities that BOTH range over the batch
(p fixed within a layer), i.e. a genuine variable x variable
multiplication. This script builds and tests, in isolation, the batched
sum-check primitive that certifies a claim of the form

        C = sum_{d in {0,1}^m}  eq(rho, d) * S~( u_d ) ,
        u_d = floor( X / a_d ),   a_d = (dstart + d) * p ,

where S is the (committed) small-side table, rho a verifier point, and the
correctness of every u_d is enforced ENTIRELY by prover-supplied trace
tables checked through low-degree multilinear identities. The verifier
never divides and never touches a quotient.

Construction (all over the d-cube of size D = 2^m unless noted):
  Witness tables: U[d]=u_d, R[d]=r_d, Qv[d]=q_d=a_d-1-r_d, the MSB-first
  bit decompositions Ub_j (Lv bits of u_d), Rb_j / Qb_j (Lr bits of r_d /
  q_d), the wiring table A[d]=a_d, and ONE[d]=1.
  Constraint zero-test (one degree-3 sum-check at a random point tau,
  constraints batched by powers of a random alpha):
    * bitness        T*(1-T)=0           for every bit table
    * recomposition  U = sum 2^j Ub_j , R = sum 2^j Rb_j , Qv = sum 2^j Qb_j
    * multiply+rem   U*a_d + R - X = 0      (a_d via the verifier's own
                                             degree-1 wiring a~(pt)=p*(e+dstart))
    * remainder bnd  a_d - 1 - R - Qv = 0   (q_d>=0, a sum of bits, => r_d<a_d)
  Together with the bit-ranges these pin u_d = floor(X/a_d) UNIQUELY.
  Lookup (the claim C): route C = sum_v S(v)*omega(v) with
  omega(v)=sum_d eq(rho,d)[u_d=v]; verify by two sum-checks --
    B1 (deg 2, v-cube): C = sum_v S~(v) omega~(v) -> claims S~(r_B), omega~(r_B);
    B2 (deg Lv+1, d-cube): omega~(r_B) = sum_d eq(rho,d) prod_j EB_j(d),
       EB_j(d) = eq(r_B[j], Ub_j(d)) the per-bit selector built from the
       SAME committed Ub_j proven above.

THE FIELD LIFT (S498). The soundness of the multiply identity U*a+R-X=0
depends on |F| being a TRUE upper bound on the integers it relates. The
demo prime q = 2^31-1 is sound only while every certified integer (values,
and the product U*a ~ X) stays below q, i.e. n <~ 30. Above that a cheating
prover can supply a *different* Lv-bit quotient U'' and remainder R'' < a
with U''*a + R'' = X + k*q for some k >= 1: this aliases to X mod q (so the
field check passes) yet U'' != floor(X/a). The bit-range tables alone do
NOT stop it -- U'' = floor((X+kq)/a) is a genuine Lv-bit integer whenever
2^Lv has headroom (precisely the regime q < X). This is a real soundness
hole of the demo above its field, not a stylistic caveat.

  *** The fix is to RAISE THE CHARACTERISTIC: a prime field F_Q with
      prime Q > (the largest integer the identity relates). ***

  q = 2^31-1   characteristic 2^31-1, sound for n <~ 30
  Q = 2^61-1   characteristic 2^61-1, sound for n <~ 60 (Mersenne, mulmod
               needs >64-bit intermediates -> Python-int object arrays here)

  A NON-FIX, recorded so it is not retried: the EXTENSION field F_{q^2}.
  |F_{q^2}| = q^2 ~ 2^62 counts the number of ELEMENTS, but its
  CHARACTERISTIC is still q. Integers embed mod the characteristic, so an
  integer k*q embeds as the degree-0 element (k*q mod q, 0) = (0,0) = the
  ZERO of F_{q^2}. The wrap-around forgery differs from the truth by
  exactly k*q, hence is the SAME element of F_{q^2} as the truth -- F_{q^2}
  cannot tell them apart no matter how the challenges are drawn (the larger
  element count only shrinks the Schwartz-Zippel error, an orthogonal
  concern). `--refute-q2` demonstrates this directly. The genuine
  alternative to a larger prime is a small-field schoolbook CARRY trace
  (limb products + carry constraints), which never asks the field to hold
  the product; that is a separate construction, not a field swap.

This file therefore runs the primitive over an arbitrary prime modulus
(--field q | big | small, or --qmod). The sum-check skeleton, the trace
identities and the cost accounting (verifier polylog, prover Õ(√x)) are all
field-agnostic; only the modulus and the numpy dtype (uint64 vs object for
>32-bit primes) change.

Leaf claims about committed witness tables (U, R, bits, S, omega at the
final sum-check points) are closed in this ISOLATED harness by direct MLE
evaluation -- the stand-in for the polynomial-commitment opening / outer-
protocol line-batching that integration supplies.

What would falsify this: an honest run rejected; a cheating prover (wrong
quotient, wrong remainder, non-Boolean bit, corrupted lookup table, wrong
claimed C, OR a wrap-around alias) accepted by a field whose characteristic
exceeds the integers in play; or the verifier work scaling with D = 2^m or
with a quotient size instead of with m and Lv.

Usage:
  python3 compressed_prover_mult_trace.py --selftest
  python3 compressed_prover_mult_trace.py --n 12 --p 7 --field q
  python3 compressed_prover_mult_trace.py --n 34 --p 7 --field big --alias-demo
  python3 compressed_prover_mult_trace.py --refute-q2
  python3 compressed_prover_mult_trace.py --bench
"""

import argparse
import time

import numpy as np

DEFAULT_Q = (1 << 31) - 1            # 2^31-1, the original demo prime
BIG_Q = (1 << 61) - 1                # 2^61-1 Mersenne prime, characteristic > X for n<~60
SMALL_Q = 1000003                    # a small prime (< X for the fast uint64 alias demo)
FIELDS = {"q": DEFAULT_Q, "big": BIG_Q, "small": SMALL_Q}


# ----------------------------------------------------------------------
# fast Mersenne-61 field arithmetic for BIG_Q = 2^61-1 (numpy uint64)
# ----------------------------------------------------------------------
# A field element of F_q for q <= 2^31-1 fits a uint64 and so does the product
# of two of them (< 2^62), so the q / small fields use plain `(a*b) % q`.  For
# BIG_Q the product reaches ~2^122 and overflows uint64; the S498/S499 lift
# therefore fell back to Python-int `object` arrays (exact but ~10-50x slower).
# These helpers compute (a*b) mod (2^61-1) on uint64 numpy arrays using a
# split-limb schoolbook product folded by the Mersenne identity 2^61 == 1
# (mod 2^61-1), keeping every intermediate < 2^64.  Set FAST_BIG to route the
# BIG_Q path through them; default False so importers (compressed_layer,
# lucy_dp_*) keep their object-array behaviour untouched until separately
# converted.  This module's own --field big / --bench-big / selftest enable it
# locally (it calls no importer), and the selftest cross-checks fast==object.

_P61 = (1 << 61) - 1
_P61u = np.uint64(_P61)
_M31 = np.uint64((1 << 31) - 1)
_M30 = np.uint64((1 << 30) - 1)
_S31 = np.uint64(31)
_S30 = np.uint64(30)
_S61 = np.uint64(61)
_U1 = np.uint64(1)
_U0 = np.uint64(0)

FAST_BIG = False   # True => uint64 Mersenne path for BIG_Q (else exact object)


def _reduce61(s):
    """s: uint64 array/scalar with s < 2^64.  Return s mod (2^61-1), canonical
    in [0, 2^61-2].  Two folds (each value drops below 2^61) then map the
    single residue 2^61-1 (== 0 mod p) to 0."""
    s = (s >> _S61) + (s & _P61u)          # < 2^61 + 7
    s = (s >> _S61) + (s & _P61u)          # in [0, 2^61-1]
    return np.where(s == _P61u, _U0, s)


def _mul61(a, b):
    """Element-wise (a*b) mod (2^61-1) for uint64 a,b in [0, 2^61).  Vectorized;
    accepts array x array, array x scalar, scalar x scalar.  All intermediates
    stay < 2^64: split a = a1*2^31 + a0 (a0<2^31, a1<2^30), same for b, so the
    four partials lo,a0b1,a1b0,hi are each < 2^62; fold mid*2^31 and hi*2^62 by
    2^61==1, 2^62==2."""
    a = np.asarray(a, dtype=np.uint64)
    b = np.asarray(b, dtype=np.uint64)
    a0 = a & _M31; a1 = a >> _S31          # a0 < 2^31, a1 < 2^30
    b0 = b & _M31; b1 = b >> _S31
    lo = a0 * b0                           # < 2^62
    mid = a0 * b1 + a1 * b0                # < 2^62
    hi = a1 * b1                           # < 2^60
    mc = mid >> _S30                       # < 2^32   (mid = mc*2^30 + ml)
    ml = mid & _M30                        # < 2^30
    # a*b == lo + mid*2^31 + hi*2^62 == lo + (mc + ml*2^31) + 2*hi   (mod p)
    s = (hi << _U1) + mc + (ml << _S31) + lo      # < 2^63
    return _reduce61(s)


def _sum61(arr):
    """Exact sum (mod 2^61-1) of a uint64 array of reduced residues (< 2^61),
    folding blocks of 4 (each block sum < 2^63) with a reduce between passes."""
    a = np.asarray(arr, dtype=np.uint64).ravel()
    if a.size == 0:
        return 0
    while a.size > 1:
        pad = (-a.size) % 4
        if pad:
            a = np.concatenate([a, np.zeros(pad, dtype=np.uint64)])
        a = _reduce61(a.reshape(-1, 4).sum(axis=1, dtype=np.uint64))
    return int(a[0])


# ----------------------------------------------------------------------
# field-parameterised helpers (uint64 fast path for q <= 2^31-1; BIG_Q via the
# Mersenne path above when FAST_BIG, else exact Python-int object arrays)
# ----------------------------------------------------------------------

def _dt(q):
    """numpy dtype: uint64 keeps a+b of two products below 2^64 only while
    q <= 2^31-1; BIG_Q fits a uint64 too and uses the Mersenne mulmod when
    FAST_BIG; other larger primes fall back to exact Python-int objects."""
    if q <= DEFAULT_Q or (FAST_BIG and q == BIG_Q):
        return np.uint64
    return object


def fast_big(q):
    """True when the uint64 Mersenne path is active for this modulus (reads the
    live FAST_BIG global, so importers see toggles)."""
    return FAST_BIG and q == BIG_Q


def fmul(a, b, q):
    """Element-wise (a*b) mod q.  Mersenne fast path for BIG_Q under FAST_BIG;
    otherwise plain `(a*b) % q` (exact for object arrays and for q <= 2^31-1)."""
    if FAST_BIG and q == BIG_Q:
        return _mul61(a, b)
    return (a * b) % q


def _asum(arr, q):
    """Exact (mod q) sum of a field array, for either dtype."""
    if FAST_BIG and q == BIG_Q and arr.dtype != object:
        return _sum61(arr)
    if arr.dtype == object:
        return int(arr.sum() % q)
    return int(arr.sum(dtype=np.uint64) % q)


def lagrange_eval(ys, r, q=DEFAULT_Q):
    """Evaluate the polynomial through (j, ys[j]), j=0..deg, at r over F_q."""
    k = len(ys) - 1
    r = int(r) % q
    total = 0
    for i in range(k + 1):
        num, den = 1, 1
        for j in range(k + 1):
            if j == i:
                continue
            num = (num * ((r - j) % q)) % q
            den = (den * ((i - j) % q)) % q
        total = (total + int(ys[i]) * num % q * pow(den, q - 2, q)) % q
    return total


def eq_table(z, q=DEFAULT_Q):
    """Table of eq~(z, w) over Boolean w (MSB-first). PROVER-side O(2^n)."""
    t = np.array([1], dtype=_dt(q))
    for zj in reversed([int(v) % q for v in z]):
        a = fmul(t, (q + 1 - zj) % q, q)
        b = fmul(t, zj, q)
        t = np.concatenate([a, b])
    return t


def eq_point(a, b, q=DEFAULT_Q):
    """eq~(a, b) for two field points. VERIFIER: O(n)."""
    acc = 1
    for x, y in zip(a, b):
        x, y = int(x) % q, int(y) % q
        acc = acc * ((x * y + (q + 1 - x) * (q + 1 - y)) % q) % q
    return acc


def mle_eval(table, pt, q=DEFAULT_Q):
    """tilde{table}(pt). PROVER-side: O(2^n)."""
    return _asum(fmul(table, eq_table(pt, q), q), q)


def int_of_point(z, q=DEFAULT_Q):
    """Int(z) = sum_j 2^{n-1-j} z_j mod q."""
    n = len(z)
    return sum((1 << (n - 1 - j)) % q * (int(z[j]) % q) for j in range(n)) % q


def sumcheck(claim, tables, terms, degree, rng, q=DEFAULT_Q):
    """Sum-check for claim = sum_w sum_terms coef * prod tables[name][w] over
    F_q. tables: dict name -> field array (consumed). Returns
    (ok, r_list, final_claim, bound_scalars). Field-agnostic: the dtype of
    the input tables (uint64 or object) carries the modulus size."""
    nvars = int(round(np.log2(len(next(iter(tables.values()))))))
    rs = []
    c = int(claim) % q
    for _ in range(nvars):
        evals = []
        for X in range(degree + 1):
            xf = X % q
            xc = (q + 1 - X) % q
            tot = 0
            for coef, names in terms:
                prod = None
                for nm in names:
                    tb = tables[nm]
                    h = len(tb) >> 1
                    row = (fmul(tb[:h], xc, q) + fmul(tb[h:], xf, q)) % q
                    prod = row if prod is None else fmul(prod, row, q)
                tot = (tot + (coef % q) * _asum(prod, q)) % q
            evals.append(tot)
        if (evals[0] + evals[1]) % q != c:
            return False, rs, None, None
        r = int(rng.integers(0, q))
        rs.append(r)
        c = lagrange_eval(evals, r, q)
        rc = (q + 1 - r) % q
        for nm in list(tables):
            tb = tables[nm]
            h = len(tb) >> 1
            tables[nm] = (fmul(tb[:h], rc, q) + fmul(tb[h:], r, q)) % q
    scalars = {nm: int(tables[nm][0]) for nm in tables}
    return True, rs, c, scalars


# ----------------------------------------------------------------------
# verifier-side O(m) wiring for a_d = (dstart + e) * p  (degree-1 in e-bits)
# ----------------------------------------------------------------------

def a_tilde(pt, p, dstart, q=DEFAULT_Q):
    """MLE of a_e = (dstart + e)*p at a field point pt (e = Int(pt)). O(m)."""
    return (p % q) * ((int_of_point(pt, q) + dstart) % q) % q


# ----------------------------------------------------------------------
# honest prover: the compressed trace over the d-cube
# ----------------------------------------------------------------------

def build_witness(X, p, m, dstart, q=DEFAULT_Q):
    """Materialize the trace tables for a_e = (dstart+e)*p, e in [0, 2^m).
    Integer arrays (u, r, q-complement, a) are kept so the witness can be
    re-materialised over a different prime (`rematerialize`) or forged
    (`forge_alias`)."""
    D = 1 << m
    e = np.arange(D, dtype=np.int64)
    a = (e + dstart) * p
    u = X // a
    r = X - u * a
    qd = a - 1 - r
    assert (r >= 0).all() and (r < a).all() and (qd >= 0).all()
    Lv = max(int(u.max()).bit_length(), 1)
    Lr = max(int(a.max()).bit_length(), 1)
    W = dict(D=D, m=m, Lv=Lv, Lr=Lr, u=u, r=r, q=qd, a=a, p=p,
             dstart=dstart, X=X, mod=q)
    _materialize(W, q)
    return W


def _materialize(W, q):
    """(Re)build the field tables of W over modulus q from its integer
    arrays. Bit tables are 0/1 (field-independent); U/R/Qv/A reduce mod q."""
    dt = _dt(q)
    u, r, qd, a = W["u"], W["r"], W["q"], W["a"]
    D, Lv, Lr = W["D"], W["Lv"], W["Lr"]
    tabs = {
        "U": (u % q).astype(dt),
        "R": (r % q).astype(dt),
        "Qv": (qd % q).astype(dt),
        "A": (a % q).astype(dt),
        "ONE": np.ones(D, dtype=dt),
    }
    for j in range(Lv):
        tabs[f"Ub{j}"] = ((u >> (Lv - 1 - j)) & 1).astype(dt)
    for j in range(Lr):
        tabs[f"Rb{j}"] = ((r >> (Lr - 1 - j)) & 1).astype(dt)
        tabs[f"Qb{j}"] = ((qd >> (Lr - 1 - j)) & 1).astype(dt)
    W["tabs"] = tabs
    W["mod"] = q
    return W


def rematerialize(W, q):
    """Switch a (possibly forged) witness to a different prime field."""
    return _materialize(W, q)


# ----------------------------------------------------------------------
# wrap-around alias forgery (the soundness hole the field lift closes)
# ----------------------------------------------------------------------

def forge_alias(W, alias_mod):
    """Return (W2, info): a copy of W with ONE entry's quotient replaced by a
    wrap-around alias u'' = floor((X + k*alias_mod)/a) that satisfies
    u''*a + r'' = X + k*alias_mod with u'' an Lv-bit integer and 0<=r''<a.
    Hence u''*a + r'' == X (mod alias_mod) but != X over the integers (and
    != X mod any prime > X + k*alias_mod). Searches entries with the most
    Lv-headroom (largest a -> smallest u) first, k>=1 up to 64. Returns
    (None, None) when no alias fits (e.g. alias_mod > X, the safe regime).
    W2 carries forged INTEGER arrays; call rematerialize(W2, test_q) before
    verifying over a chosen field."""
    X, Lv = W["X"], W["Lv"]
    cap = 1 << Lv
    a_arr, u_arr = W["a"], W["u"]
    for e0 in np.argsort(-a_arr):          # large a first => most headroom
        a = int(a_arr[e0])
        u_true = int(u_arr[e0])
        for k in range(1, 65):
            num = X + k * alias_mod
            upp = num // a
            rpp = num - upp * a
            if upp != u_true and upp < cap and 0 <= rpp < a:
                W2 = dict(W)
                W2["u"] = W["u"].copy()
                W2["r"] = W["r"].copy()
                W2["q"] = W["q"].copy()
                W2["u"][e0] = upp
                W2["r"][e0] = rpp
                W2["q"][e0] = a - 1 - rpp
                return W2, dict(e0=int(e0), k=int(k), a=a, u_true=u_true,
                                u_forged=int(upp), r_forged=int(rpp))
    return None, None


# ----------------------------------------------------------------------
# the batched constraint zero-test (degree-3 sum-check)
# ----------------------------------------------------------------------

def build_terms(Lv, Lr, X, alpha, q):
    """Flat (coef,[names]) term list for F = sum_c alpha^c * constraint_c,
    every term carrying the eq-selector EQ. Order MUST match
    constraint_eval below (selftest checks honest equivalence)."""
    ap = [pow(int(alpha), k, q) for k in range(6 + Lv + 2 * Lr)]
    negX = (q - (X % q)) % q
    terms = []
    a1 = ap[1]                                              # c1: U*a + R - X
    terms += [(a1, ["EQ", "U", "A"]), (a1, ["EQ", "R"]),
              (a1 * negX % q, ["EQ", "ONE"])]
    a2 = ap[2]                                              # c2: a - 1 - R - Qv
    terms += [(a2, ["EQ", "A"]), ((q - a2) % q, ["EQ", "ONE"]),
              ((q - a2) % q, ["EQ", "R"]), ((q - a2) % q, ["EQ", "Qv"])]
    a3 = ap[3]                                              # c3: recompU
    terms.append((a3, ["EQ", "U"]))
    for j in range(Lv):
        terms.append(((q - a3 * pow(2, Lv - 1 - j, q)) % q, ["EQ", f"Ub{j}"]))
    a4 = ap[4]                                              # c4: recompR
    terms.append((a4, ["EQ", "R"]))
    for j in range(Lr):
        terms.append(((q - a4 * pow(2, Lr - 1 - j, q)) % q, ["EQ", f"Rb{j}"]))
    a5 = ap[5]                                              # c5: recompQ
    terms.append((a5, ["EQ", "Qv"]))
    for j in range(Lr):
        terms.append(((q - a5 * pow(2, Lr - 1 - j, q)) % q, ["EQ", f"Qb{j}"]))
    c = 6                                                   # c6..: bitness
    for nm in ([f"Ub{j}" for j in range(Lv)]
               + [f"Rb{j}" for j in range(Lr)]
               + [f"Qb{j}" for j in range(Lr)]):
        ac = ap[c]
        c += 1
        terms += [(ac, ["EQ", nm]), ((q - ac) % q, ["EQ", nm, nm])]
    return terms


def constraint_eval(scal, av, X, Lv, Lr, alpha, q):
    """Verifier's recomputation of F/eq at the final point from the bound
    witness scalars (U,R,Qv,bits) and its OWN wiring value av = a~(r). The
    prover's EQ/A/ONE scalars are ignored -- av ties A to the true wiring."""
    ap = [pow(int(alpha), k, q) for k in range(6 + Lv + 2 * Lr)]
    U, R, Qv = scal["U"], scal["R"], scal["Qv"]
    val = ap[1] * ((U * av + R - X) % q) % q
    val = (val + ap[2] * ((av - 1 - R - Qv) % q)) % q
    su = U
    for j in range(Lv):
        su = (su - pow(2, Lv - 1 - j, q) * scal[f"Ub{j}"]) % q
    val = (val + ap[3] * su) % q
    sr = R
    for j in range(Lr):
        sr = (sr - pow(2, Lr - 1 - j, q) * scal[f"Rb{j}"]) % q
    val = (val + ap[4] * sr) % q
    sq = Qv
    for j in range(Lr):
        sq = (sq - pow(2, Lr - 1 - j, q) * scal[f"Qb{j}"]) % q
    val = (val + ap[5] * sq) % q
    c = 6
    for nm in ([f"Ub{j}" for j in range(Lv)]
               + [f"Rb{j}" for j in range(Lr)]
               + [f"Qb{j}" for j in range(Lr)]):
        T = scal[nm]
        val = (val + ap[c] * ((T - T * T) % q)) % q
        c += 1
    return val % q


def verify_constraints(W, rng, stats, q=DEFAULT_Q):
    """Degree-3 zero-test certifying u_d = floor(X/a_d) for the whole batch
    over F_q. Returns ok. (W["tabs"] must already be materialised over q;
    build_witness/rematerialize do this.)"""
    Lv, Lr, m = W["Lv"], W["Lr"], W["m"]
    X, p, dstart = W["X"], W["p"], W["dstart"]
    tau = [int(rng.integers(0, q)) for _ in range(m)]
    alpha = int(rng.integers(2, q))
    t0 = time.perf_counter()
    tabs = dict(W["tabs"])
    tabs["EQ"] = eq_table(tau, q)
    terms = build_terms(Lv, Lr, X, alpha, q)
    stats["t_prover"] += time.perf_counter() - t0
    ok, r_A, finalA, scal = sumcheck(0, tabs, terms, 3, rng, q)
    stats["comm"] += 4 * m                                  # degree-3 round polys
    if not ok:
        return False
    t0 = time.perf_counter()
    av = a_tilde(r_A, p, dstart, q)
    eqv = eq_point(tau, r_A, q)
    expect = eqv * constraint_eval(scal, av, X, Lv, Lr, alpha, q) % q
    res = (finalA % q) == expect
    stats["t_verifier"] += time.perf_counter() - t0
    return res


# ----------------------------------------------------------------------
# the batched lookup  C = sum_d eq(rho,d) S(u_d)
# ----------------------------------------------------------------------

def honest_C(W, S_tab, rho, q):
    eqr = eq_table(rho, q)
    return _asum(fmul(eqr, S_tab[W["u"]], q), q)


def build_omega(W, rho, q):
    eqr = eq_table(rho, q)
    if FAST_BIG and q == BIG_Q:
        # < D summands each < q overflow uint64; accumulate exactly in object
        acc = np.zeros(1 << W["Lv"], dtype=object)
        np.add.at(acc, W["u"], eqr.astype(object))
        return (acc % q).astype(np.uint64)
    omega = np.zeros(1 << W["Lv"], dtype=_dt(q))
    np.add.at(omega, W["u"], eqr)            # < D summands each < q
    return omega % q


def verify_lookup(W, S_tab, rho, claimedC, omega, rng, stats, q):
    """B1 then B2; leaf claims on S / Ub_j closed by direct MLE opening."""
    Lv, m = W["Lv"], W["m"]
    # B1: claimedC = sum_v S(v) omega(v)
    okB1, r_B, finalB1, sB1 = sumcheck(claimedC,
                                       {"S": S_tab.copy(), "W": omega.copy()},
                                       [(1, ["S", "W"])], 2, rng, q)
    stats["comm"] += 3 * Lv
    if not okB1:
        return False
    t0 = time.perf_counter()
    okf = (finalB1 % q) == sB1["S"] * sB1["W"] % q
    okf = okf and sB1["S"] == mle_eval(S_tab, r_B, q)       # opening of S
    stats["t_verifier"] += time.perf_counter() - t0
    if not okf:
        return False
    sW = sB1["W"]                                           # claim for B2
    # B2: sW = sum_d eq(rho,d) prod_j EB_j(d),  EB_j = eq(r_B[j], Ub_j)
    t0 = time.perf_counter()
    tabs = {"EQrho": eq_table(rho, q)}
    ebnames = []
    for j in range(Lv):
        rb = int(r_B[j]) % q
        Ub = W["tabs"][f"Ub{j}"]
        EB = (((1 - rb) % q) + ((2 * rb - 1) % q) * Ub) % q
        nm = f"EB{j}"
        tabs[nm] = EB.astype(_dt(q))
        ebnames.append(nm)
    stats["t_prover"] += time.perf_counter() - t0
    okB2, r_C, finalB2, _ = sumcheck(sW, tabs, [(1, ["EQrho"] + ebnames)],
                                     Lv + 1, rng, q)
    stats["comm"] += (Lv + 2) * m
    if not okB2:
        return False
    t0 = time.perf_counter()
    expect = eq_point(rho, r_C, q)
    for j in range(Lv):
        rb = int(r_B[j]) % q
        ub = mle_eval(W["tabs"][f"Ub{j}"], r_C, q)         # opening of Ub_j
        ebj = (((1 - rb) % q) + ((2 * rb - 1) % q) * ub) % q
        expect = expect * ebj % q
    res = (finalB2 % q) == expect
    stats["t_verifier"] += time.perf_counter() - t0
    return res


# ----------------------------------------------------------------------
# full primitive with cheat injection
# ----------------------------------------------------------------------

def run_primitive(X, p, m, dstart, rng, q=DEFAULT_Q, S_seed=0, cheat=None):
    """Build the trace, the lookup claim, optionally corrupt, run both the
    constraint zero-test and the lookup over F_q. Returns dict(accepted, ...)."""
    W = build_witness(X, p, m, dstart, q)
    Lv = W["Lv"]
    S_tab = (np.random.default_rng(S_seed).integers(0, min(q, 1 << 62),
             size=1 << Lv) % q).astype(_dt(q))
    rho = [int(rng.integers(0, q)) for _ in range(m)]
    C = honest_C(W, S_tab, rho, q)
    omega = build_omega(W, rho, q)
    claimedC = C

    if cheat == "u_consistent":          # self-consistent wrong quotient
        d = 1 % W["D"]
        nu = int(W["u"][d]) + 1
        W["tabs"]["U"] = W["tabs"]["U"].copy(); W["tabs"]["U"][d] = nu % q
        for j in range(Lv):
            W["tabs"][f"Ub{j}"] = W["tabs"][f"Ub{j}"].copy()
            W["tabs"][f"Ub{j}"][d] = (nu >> (Lv - 1 - j)) & 1
    elif cheat == "u_value":             # quotient off, bits untouched
        d = 1 % W["D"]
        W["tabs"]["U"] = W["tabs"]["U"].copy()
        W["tabs"]["U"][d] = (int(W["tabs"]["U"][d]) + 1) % q
    elif cheat == "r_value":             # remainder corrupted
        d = 1 % W["D"]
        W["tabs"]["R"] = W["tabs"]["R"].copy()
        W["tabs"]["R"][d] = (int(W["tabs"]["R"][d]) + 1) % q
    elif cheat == "nonbit":              # a 'bit' set to 2
        W["tabs"]["Ub0"] = W["tabs"]["Ub0"].copy()
        W["tabs"]["Ub0"][0] = 2
    elif cheat == "wrong_C":             # lie about the claimed value
        claimedC = (C + 1) % q
    elif cheat == "omega_route":         # corrupt routing, keep <S,omega>=C
        omega = omega.copy()
        omega[0] = (int(omega[0]) + 1) % q
        omega[-1] = (int(omega[-1]) - 1) % q
    elif cheat == "u_alias":             # wrap-around: u'' aliases X mod q
        W2, info = forge_alias(W, q)
        if W2 is not None:               # only constructible while q < X
            rematerialize(W2, q)
            W = W2
            # honest lookup keeps the true u (only the trace is forged) ->
            # this isolates the constraint test, the seat of the hole
            C = honest_C(W, S_tab, rho, q)
            omega = build_omega(W, rho, q)
            claimedC = C

    stats = {"t_prover": 0.0, "t_verifier": 0.0, "comm": 0}
    okC = verify_constraints(W, rng, stats, q)
    okL = verify_lookup(W, S_tab, rho, claimedC, omega, rng, stats, q)
    return dict(accepted=(okC and okL), okC=okC, okL=okL, D=W["D"],
                Lv=Lv, Lr=W["Lr"], m=m, trueC=C, **stats)


# ----------------------------------------------------------------------
# F_{q^2} refutation: an extension field does NOT lift the integer range
# ----------------------------------------------------------------------

class GF2:
    """F_{q^2} = F_q[t]/(t^2 + 1), valid when q % 4 == 3 (so -1 is a
    quadratic non-residue). Element a + b*t. Characteristic is q."""
    __slots__ = ("a", "b", "q")

    def __init__(self, a, b, q):
        self.a = a % q
        self.b = b % q
        self.q = q

    def __add__(s, o):
        return GF2(s.a + o.a, s.b + o.b, s.q)

    def __sub__(s, o):
        return GF2(s.a - o.a, s.b - o.b, s.q)

    def __mul__(s, o):                            # t^2 = -1
        q = s.q
        return GF2((s.a * o.a - s.b * o.b) % q, (s.a * o.b + s.b * o.a) % q, q)

    def __eq__(s, o):
        return isinstance(o, GF2) and s.a == o.a and s.b == o.b

    def is_zero(s):
        return s.a == 0 and s.b == 0

    def __repr__(s):
        return f"({s.a}+{s.b}t)"

    @staticmethod
    def embed(x, q):
        return GF2(x % q, 0, q)

    @staticmethod
    def rand(rng, q):
        return GF2(int(rng.integers(0, q)), int(rng.integers(0, q)), q)


def _eqw_sum_c1(W, embed, eq_of, zero, tau, X):
    """Compute  sum_e eq(tau, e) * c1(e),  c1(e) = U[e]*A[e] + R[e] - X,
    directly over a field given by (embed, multiply/add via field objects,
    eq_of(tau,e), additive zero). This is EXACTLY the quantity the constraint
    zero-test certifies equals 0 for the multiply identity; the other
    constraints (bitness, recompositions, remainder bound) are satisfied as
    EXACT integers by the alias forgery, hence 0 in every field. So this sum
    is the whole soundness question for the wrap-around."""
    D = W["D"]
    u, a, r = W["u"], W["a"], W["r"]
    eX = embed(X)
    total = zero
    for e in range(D):
        c1 = embed(int(u[e])) * embed(int(a[e])) + embed(int(r[e])) - eX
        total = total + eq_of(tau, e) * c1
    return total


def refute_q2(n=33, p=7, m=4, seed=0, q0=DEFAULT_Q):
    """Demonstrate that the extension field F_{q0^2} CANNOT detect the
    wrap-around alias, while a true larger prime can. Builds an alias forgery
    against the base prime q0 (needs X > q0, so n > 31 for q0 = 2^31-1) and
    evaluates the multiply-identity zero-test over three fields."""
    assert q0 % 4 == 3, "GF2 with t^2=-1 needs q0 % 4 == 3"
    X = (1 << n) - 1
    dstart = max(1, int(X ** 0.5) // p + 1)
    W = build_witness(X, p, m, dstart, q=q0)
    W2, info = forge_alias(W, q0)
    assert W2 is not None, f"no alias fits at n={n}; raise n above log2(q0)"
    e0 = info["e0"]

    rng = np.random.default_rng(seed)
    # the integer the forged entry violates the multiply identity by:
    delta = info["k"] * q0                              # u''*a+r'' - X = k*q0

    out = {"info": info, "delta": delta}

    # --- F_{q0}: the base prime (too small) -- aliases to zero ---
    # (the prime path accumulates unreduced ints; reduce at the end, since
    #  mod is a ring homomorphism -> this is exactly the field sum)
    tau_p = [int(rng.integers(0, q0)) for _ in range(m)]
    sp = _eqw_sum_c1(W2, lambda x: x % q0,
                     lambda t, e: _eq_bool_prime(t, e, m, q0),
                     0, tau_p, X)
    out["base_prime_sum"] = sp % q0                      # 0 -> undetected

    # --- F_{q0^2}: the EXTENSION field -- challenges drawn from all of F_{q0^2} ---
    tau_g = [GF2.rand(rng, q0) for _ in range(m)]
    sg = _eqw_sum_c1(W2, lambda x: GF2.embed(x, q0),
                     lambda t, e: _eq_bool_gf2(t, e, m, q0),
                     GF2(0, 0, q0), tau_g, X)
    out["ext_field_sum"] = sg                           # (0,0) -> undetected

    # --- F_Q: a genuine larger prime (Q > X + k*q0) -- detects it ---
    Q = BIG_Q
    assert Q > X + delta
    tau_b = [int(rng.integers(0, Q)) for _ in range(m)]
    sb = _eqw_sum_c1(W2, lambda x: x % Q,
                     lambda t, e: _eq_bool_prime(t, e, m, Q),
                     0, tau_b, X)
    out["big_prime_sum"] = sb % Q                        # != 0 -> detected

    # honest control over the extension field -> still (0,0)
    sh = _eqw_sum_c1(W, lambda x: GF2.embed(x, q0),
                     lambda t, e: _eq_bool_gf2(t, e, m, q0),
                     GF2(0, 0, q0), tau_g, X)
    out["ext_field_honest"] = sh
    return out


def _eq_bool_prime(tau, e, m, q):
    """eq~(tau, e) for Boolean e (MSB-first), tau a prime-field point."""
    acc = 1
    for j in range(m):
        bit = (e >> (m - 1 - j)) & 1
        tj = int(tau[j]) % q
        acc = acc * (tj if bit else (q + 1 - tj)) % q
    return acc


def _eq_bool_gf2(tau, e, m, q):
    """eq~(tau, e) for Boolean e (MSB-first), tau a GF2 point."""
    acc = GF2(1, 0, q)
    one = GF2(1, 0, q)
    for j in range(m):
        bit = (e >> (m - 1 - j)) & 1
        acc = acc * (tau[j] if bit else (one - tau[j]))
    return acc


# ----------------------------------------------------------------------
# tests / experiments
# ----------------------------------------------------------------------

def auto_params(n, p):
    X = (1 << n) - 1
    m = max(2, n // 2)
    dstart = max(1, int(X ** 0.5) // p + 1)      # a_d > sqrt X => u_d < sqrt X
    return X, m, dstart


def selftest():
    rng = np.random.default_rng(0)
    # 0. field sanity: BIG_Q / SMALL_Q primality + characteristic argument
    assert all(BIG_Q % d for d in range(2, 1000)) and pow(2, BIG_Q - 1, BIG_Q) == 1
    assert all(SMALL_Q % d for d in range(2, int(SMALL_Q ** 0.5) + 1))
    for k in (1, 2, 7):                          # k*q is the field zero of F_{q^2}
        assert (k * DEFAULT_Q) % DEFAULT_Q == 0

    # 1. honest trace relations hold exactly (incl. p=2 and a 1-bit-u case)
    for (n, p) in [(12, 7), (12, 2), (10, 3), (8, 5), (14, 11)]:
        X, m, dstart = auto_params(n, p)
        W = build_witness(X, p, m, dstart)
        assert (W["u"] * W["a"] + W["r"] == X).all()
        assert (W["q"] == W["a"] - 1 - W["r"]).all()
        Lv = W["Lv"]
        ru = np.zeros(W["D"], dtype=np.int64)
        for j in range(Lv):
            ru += (1 << (Lv - 1 - j)) * W["tabs"][f"Ub{j}"].astype(np.int64)
        assert (ru == W["u"]).all()

    # 2. honest_C == <S,omega> (routing identity), default + a big prime
    for (n, p, sd, q) in [(12, 7, 1, DEFAULT_Q), (10, 3, 2, DEFAULT_Q),
                          (12, 7, 3, BIG_Q)]:
        X, m, dstart = auto_params(n, p)
        W = build_witness(X, p, m, dstart, q)
        S = (np.random.default_rng(sd).integers(0, min(q, 1 << 62),
             size=1 << W["Lv"]) % q).astype(_dt(q))
        rho = [int(rng.integers(0, q)) for _ in range(m)]
        C = honest_C(W, S, rho, q)
        om = build_omega(W, rho, q)
        assert _asum((S * om) % q, q) == C

    # 3. honest end-to-end accepts; soundness against each cheat class.
    #    Run over the default prime AND the big prime (object dtype) to prove
    #    the lift's arithmetic is correct in the full machinery.
    cheats = ["u_consistent", "u_value", "r_value", "nonbit",
              "wrong_C", "omega_route"]
    for q in (DEFAULT_Q, BIG_Q):
        for (n, p) in [(12, 7), (10, 3)]:
            X, m, dstart = auto_params(n, p)
            res = run_primitive(X, p, m, dstart, np.random.default_rng(5), q)
            assert res["accepted"], (q, n, p, res)
            for ch in cheats:
                rej = sum(not run_primitive(X, p, m, dstart,
                          np.random.default_rng(100 + t), q, cheat=ch)["accepted"]
                          for t in range(6))
                assert rej == 6, (q, n, p, ch, rej)

    # 4. constraint test alone catches a corrupted trace
    X, m, dstart = auto_params(12, 7)
    W = build_witness(X, p=7, m=m, dstart=dstart)
    st = {"t_prover": 0.0, "t_verifier": 0.0, "comm": 0}
    assert verify_constraints(W, np.random.default_rng(9), st)
    W["tabs"]["R"] = W["tabs"]["R"].copy(); W["tabs"]["R"][2] += 1
    bad = sum(not verify_constraints(W, np.random.default_rng(200 + t), st)
              for t in range(8))
    assert bad == 8

    # 5. THE FIELD LIFT (fast, uint64): a wrap-around alias built against a
    #    too-small prime is ACCEPTED by that prime but REJECTED by a prime
    #    that exceeds the integers in play. n=24 with SMALL_Q < X < DEFAULT_Q.
    n, p = 24, 7
    X, m, dstart = auto_params(n, p)
    Wh = build_witness(X, p, m, dstart, q=SMALL_Q)
    W2, info = forge_alias(Wh, SMALL_Q)
    assert W2 is not None and info["u_forged"] != info["u_true"]
    assert info["u_forged"] < (1 << Wh["Lv"])             # fits Lv bits
    # the forged integer relation: u''*a + r'' = X + k*SMALL_Q
    assert (info["u_forged"] * info["a"] + info["r_forged"]
            == X + info["k"] * SMALL_Q)
    # ACCEPTED over the too-small prime (aliases), deterministically:
    rematerialize(W2, SMALL_Q)
    stA = {"t_prover": 0.0, "t_verifier": 0.0, "comm": 0}
    assert verify_constraints(W2, np.random.default_rng(1), stA, SMALL_Q)
    # REJECTED over the bigger prime (no wrap), every trial:
    rematerialize(W2, DEFAULT_Q)
    rej = sum(not verify_constraints(W2, np.random.default_rng(300 + t),
              {"t_prover": 0.0, "t_verifier": 0.0, "comm": 0}, DEFAULT_Q)
              for t in range(8))
    assert rej == 8, ("alias not rejected by lifted prime", rej)
    # honest witness still accepted over both primes
    assert verify_constraints(build_witness(X, p, m, dstart, SMALL_Q),
                              np.random.default_rng(2),
                              {"t_prover": 0.0, "t_verifier": 0.0, "comm": 0}, SMALL_Q)
    assert verify_constraints(build_witness(X, p, m, dstart, DEFAULT_Q),
                              np.random.default_rng(2),
                              {"t_prover": 0.0, "t_verifier": 0.0, "comm": 0}, DEFAULT_Q)

    # 6. THE PRODUCTION BOUNDARY: at n=33 (X > 2^31) the DEMO prime DEFAULT_Q
    #    itself aliases; BIG_Q (2^61-1) rejects. Small m keeps object arrays fast.
    n2, p2, m2 = 33, 7, 4
    X2 = (1 << n2) - 1
    ds2 = max(1, int(X2 ** 0.5) // p2 + 1)
    Wb = build_witness(X2, p2, m2, ds2, q=DEFAULT_Q)
    Wf, info2 = forge_alias(Wb, DEFAULT_Q)               # alias mod 2^31-1
    assert Wf is not None
    rematerialize(Wf, DEFAULT_Q)
    assert verify_constraints(Wf, np.random.default_rng(7),
                              {"t_prover": 0.0, "t_verifier": 0.0, "comm": 0}, DEFAULT_Q), \
        "demo prime should accept its own alias"
    rematerialize(Wf, BIG_Q)
    rej2 = sum(not verify_constraints(Wf, np.random.default_rng(400 + t),
               {"t_prover": 0.0, "t_verifier": 0.0, "comm": 0}, BIG_Q)
               for t in range(5))
    assert rej2 == 5, ("BIG_Q failed to reject the n=34 alias", rej2)

    # 7. F_{q^2} REFUTATION: the extension field does NOT detect the alias.
    r = refute_q2(n=33, p=7, m=4)
    assert r["base_prime_sum"] == 0                       # too-small prime: undetected
    assert r["ext_field_sum"].is_zero()                   # EXTENSION FIELD: undetected
    assert r["ext_field_honest"].is_zero()                # honest control
    assert r["big_prime_sum"] != 0                        # larger prime: DETECTED

    # 8. FAST MERSENNE-61 PATH (the speed fix): the uint64 mulmod/sum primitives
    #    are bit-for-bit exact, and the WHOLE BIG_Q protocol over the fast path
    #    is identical to the object-array reference (same C, accept/reject, comm).
    global FAST_BIG
    saved = FAST_BIG
    try:
        # 8a. primitives vs Python-int reference (random + edge values)
        prng = np.random.default_rng(11)
        edges = [0, 1, 2, _P61 - 1, _P61 - 2, 1 << 60, (1 << 60) - 1,
                 1 << 31, (1 << 31) - 1, 1 << 30, (1 << 61) - 2]
        ea = np.array(edges, dtype=np.uint64)
        eb = np.array(list(reversed(edges)), dtype=np.uint64)
        assert (_mul61(ea, eb) ==
                np.array([(int(x) * int(y)) % _P61 for x, y in zip(ea, eb)],
                         dtype=np.uint64)).all()
        for _ in range(30):
            k = int(prng.integers(1, 400))
            aa = prng.integers(0, _P61, size=k, dtype=np.uint64)
            bb = prng.integers(0, _P61, size=k, dtype=np.uint64)
            mm = _mul61(aa, bb)
            assert (mm == np.array([(int(x) * int(y)) % _P61
                                    for x, y in zip(aa, bb)],
                                   dtype=np.uint64)).all()
            assert (mm < _P61).all()                      # canonical residues
            assert _sum61(mm) == int(sum(int(x) for x in mm) % _P61)
            sc = int(prng.integers(0, _P61))              # array x scalar
            assert (_mul61(aa, sc) ==
                    np.array([(int(x) * sc) % _P61 for x in aa],
                             dtype=np.uint64)).all()
        # 8b. whole BIG_Q protocol: fast (uint64) == object reference, bit-identical
        cheats8 = ["u_consistent", "u_value", "r_value", "nonbit",
                   "wrong_C", "omega_route"]
        for (n, p) in [(12, 7), (10, 3), (16, 5), (20, 7)]:
            X, m, dstart = auto_params(n, p)
            FAST_BIG = False
            ref = run_primitive(X, p, m, dstart, np.random.default_rng(5), BIG_Q)
            FAST_BIG = True
            fst = run_primitive(X, p, m, dstart, np.random.default_rng(5), BIG_Q)
            assert ref["trueC"] == fst["trueC"], (n, p, ref["trueC"], fst["trueC"])
            assert ref["accepted"] and fst["accepted"]
            assert ref["comm"] == fst["comm"]
            for ch in cheats8:
                FAST_BIG = False
                ro = [run_primitive(X, p, m, dstart, np.random.default_rng(100 + t),
                                    BIG_Q, cheat=ch)["accepted"] for t in range(4)]
                FAST_BIG = True
                rf = [run_primitive(X, p, m, dstart, np.random.default_rng(100 + t),
                                    BIG_Q, cheat=ch)["accepted"] for t in range(4)]
                assert ro == rf == [False, False, False, False], (n, p, ch, ro, rf)
        # 8c. the alias forgery is still rejected by BIG_Q on the FAST path
        n8, p8, m8 = 33, 7, 4
        X8 = (1 << n8) - 1
        ds8 = max(1, int(X8 ** 0.5) // p8 + 1)
        Wb8 = build_witness(X8, p8, m8, ds8, q=DEFAULT_Q)
        Wf8, _ = forge_alias(Wb8, DEFAULT_Q)
        assert Wf8 is not None
        rematerialize(Wf8, BIG_Q)
        rej8 = sum(not verify_constraints(Wf8, np.random.default_rng(400 + t),
                   {"t_prover": 0.0, "t_verifier": 0.0, "comm": 0}, BIG_Q)
                   for t in range(5))
        assert rej8 == 5, ("FAST BIG_Q failed to reject the alias", rej8)
    finally:
        FAST_BIG = saved

    print("selftest OK")


def bench(p=7):
    """Verifier work vs batch size D = 2^m over the default prime."""
    print(f"{'n':>3} {'m':>3} {'D':>7} {'2^Lv':>7} {'Lr':>4} "
          f"{'t_prover_ms':>12} {'t_verifier_ms':>14} {'comm_elems':>11}")
    for n in (8, 10, 12, 14, 16):
        X, m, dstart = auto_params(n, p)
        res = run_primitive(X, p, m, dstart, np.random.default_rng(n))
        assert res["accepted"]
        print(f"{n:>3} {m:>3} {res['D']:>7} {1<<res['Lv']:>7} {res['Lr']:>4} "
              f"{res['t_prover']*1000:>12.2f} {res['t_verifier']*1000:>14.3f} "
              f"{res['comm']:>11}")


def bench_big(p=7, ns=(12, 16, 20, 24)):
    """The speed fix, measured: the full BIG_Q (2^61-1) trace-certifier run over
    the uint64 Mersenne path vs the exact Python-int object-array reference.
    Reports the prover/verifier wall and the fast/object ratio; asserts the two
    paths agree on the certified value C (so the speedup is free of soundness)."""
    global FAST_BIG
    saved = FAST_BIG
    print(f"BIG_Q = 2^61-1 trace certifier: uint64 Mersenne vs object arrays "
          f"(p={p})")
    print(f"{'n':>3} {'D':>7} {'2^Lv':>7} {'obj_ms':>10} {'fast_ms':>10} "
          f"{'speedup':>8}  {'C match':>7}")
    try:
        for n in ns:
            X, m, dstart = auto_params(n, p)
            FAST_BIG = False
            t0 = time.perf_counter()
            ro = run_primitive(X, p, m, dstart, np.random.default_rng(n), BIG_Q)
            t_obj = time.perf_counter() - t0
            FAST_BIG = True
            t0 = time.perf_counter()
            rf = run_primitive(X, p, m, dstart, np.random.default_rng(n), BIG_Q)
            t_fast = time.perf_counter() - t0
            ok = ro["trueC"] == rf["trueC"] and ro["accepted"] and rf["accepted"]
            print(f"{n:>3} {ro['D']:>7} {1<<ro['Lv']:>7} {t_obj*1000:>10.1f} "
                  f"{t_fast*1000:>10.1f} {t_obj/max(t_fast,1e-9):>7.1f}x  "
                  f"{'yes' if ok else 'NO!':>7}")
            assert ok, (n, ro["trueC"], rf["trueC"])
    finally:
        FAST_BIG = saved


def alias_demo(n, p, q):
    """Show the wrap-around alias accepted by a too-small field and rejected
    by a field whose characteristic exceeds the integers in play."""
    X = (1 << n) - 1
    dstart = max(1, int(X ** 0.5) // p + 1)
    W = build_witness(X, p, max(2, n // 2), dstart, q=q)
    W2, info = forge_alias(W, q)
    if W2 is None:
        print(f"no alias fits at n={n} over q={q} (q > X: this field is already "
              f"sound here)")
        return
    print(f"alias forged against q={q}: entry e0={info['e0']}, k={info['k']}, "
          f"a={info['a']}, true u={info['u_true']} -> forged u={info['u_forged']} "
          f"(< 2^Lv={1<<W['Lv']}), r={info['r_forged']}")
    print(f"  u''*a + r'' = {info['u_forged']*info['a']+info['r_forged']} "
          f"= X + {info['k']}*q  (X={X})")
    rematerialize(W2, q)
    st = {"t_prover": 0.0, "t_verifier": 0.0, "comm": 0}
    acc = verify_constraints(W2, np.random.default_rng(1), st, q)
    print(f"  over q={q:>19}: {'ACCEPTED (aliases!)' if acc else 'rejected'}")
    for lift in (DEFAULT_Q, BIG_Q):
        if lift <= q:
            continue
        rematerialize(W2, lift)
        rej = sum(not verify_constraints(W2, np.random.default_rng(500 + t),
                  {"t_prover": 0.0, "t_verifier": 0.0, "comm": 0}, lift)
                  for t in range(8))
        print(f"  over lift={lift:>16}: rejected {rej}/8")


def main(n, p, cheat_trials, seed, field, alias):
    global FAST_BIG
    q = FIELDS[field]
    if field == "big":
        FAST_BIG = True                  # use the uint64 Mersenne speed path
    if alias:
        alias_demo(n, p, q)
        return
    X, m, dstart = auto_params(n, p)
    res = run_primitive(X, p, m, dstart, np.random.default_rng(seed), q)
    print(f"field = {field} (q = {q}), X = 2^{n}-1 = {X}, p = {p}, "
          f"a_d = (d+{dstart})*{p}")
    print(f"batch D = {res['D']} = 2^{m}, small-side 2^Lv = {1<<res['Lv']}, "
          f"Lr = {res['Lr']}")
    print(f"honest run: {'ACCEPTED' if res['accepted'] else 'REJECTED'} "
          f"(constraints={res['okC']}, lookup={res['okL']}), "
          f"verified C = sum_d eq(rho,d) S(floor(X/a_d)) = {res['trueC']}")
    assert res["accepted"]
    print(f"t_prover = {res['t_prover']*1000:.2f} ms   "
          f"t_verifier = {res['t_verifier']*1000:.3f} ms   "
          f"comm = {res['comm']} field elems")
    for ch in ["u_consistent", "u_value", "r_value", "nonbit", "wrong_C",
               "omega_route"]:
        rej = sum(not run_primitive(X, p, m, dstart,
                                    np.random.default_rng(seed + 50 + t),
                                    q, cheat=ch)["accepted"]
                  for t in range(cheat_trials))
        print(f"  cheat {ch:>13}: rejected {rej}/{cheat_trials}")
    Lv = res["Lv"]
    print(f"soundness bound ~ (deg*vars)/q = "
          f"{(3*m + 2*Lv + (Lv+1)*m)/q:.2e}")
    if q < X:
        print(f"  NOTE: q < X here -> wrap-around aliasing is possible at this "
              f"field; use --field big or --alias-demo.")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, default=12)
    ap.add_argument("--p", type=int, default=7)
    ap.add_argument("--field", choices=list(FIELDS), default="q")
    ap.add_argument("--cheat-trials", type=int, default=20)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--alias-demo", action="store_true",
                    help="forge a wrap-around alias and show the field lift")
    ap.add_argument("--refute-q2", action="store_true",
                    help="show the extension field F_{q^2} does NOT fix aliasing")
    ap.add_argument("--selftest", action="store_true")
    ap.add_argument("--bench", action="store_true")
    ap.add_argument("--bench-big", action="store_true",
                    help="BIG_Q uint64 Mersenne vs object-array speedup")
    args = ap.parse_args()
    if args.selftest:
        selftest()
    elif args.bench:
        bench()
    elif args.bench_big:
        bench_big()
    elif args.refute_q2:
        r = refute_q2()
        print(f"alias forged against q=2^31-1 (k={r['info']['k']}, entry "
              f"e0={r['info']['e0']}): the multiply identity is off by "
              f"delta = k*q = {r['delta']}")
        print(f"  zero-test sum over F_q   (q=2^31-1, too small): "
              f"{r['base_prime_sum']}  -> {'UNDETECTED' if r['base_prime_sum']==0 else 'detected'}")
        print(f"  zero-test sum over F_q^2 (extension, |F|~2^62):  "
              f"{r['ext_field_sum']}  -> "
              f"{'UNDETECTED (extension does NOT help)' if r['ext_field_sum'].is_zero() else 'detected'}")
        print(f"  zero-test sum over F_Q   (Q=2^61-1, > X):       "
              f"{r['big_prime_sum']}  -> "
              f"{'detected (correct lift)' if r['big_prime_sum']!=0 else 'UNDETECTED'}")
        print(f"  honest control over F_q^2:                      "
              f"{r['ext_field_honest']}  (zero, as it must be)")
    else:
        main(args.n, args.p, args.cheat_trials, args.seed, args.field,
             args.alias_demo)
