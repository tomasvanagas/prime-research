#!/usr/bin/env python3
"""
merged_layer_comm.py -- S510.  Is the Õ(sqrt x) certificate SIZE inherent to the
LAYERING, or just an artifact of the unbatched K = pi(sqrt x)-layer chain?

CONTEXT.  S509 (certificate_profile.py) measured the full succinct certificate for
"pi(x) = c": size alpha = 0.473 (Theta(sqrt x)), DOMINATED by comm_outer (the K
SEQUENTIAL two-sided layer reductions; every batchable obligation already polylog,
the base committed to x^{1/4}).  It flagged ONE open quantitative question: could
MERGING j consecutive sieve layers into one reduction drop the comm exponent below
0.5?  Open problem 2 (CLOSED NEGATIVE) showed merging blows up the PROVER (a j-layer
merge has 2^j fill per row -- one entry per subset-product floor(v / p_{i1}...p_{is});
tree cost m*K^2/2 ~ x^{3/2}); it never measured the COMM of a merged-layer verifier.
This file does, as a STANDALONE primitive (one merged layer's comm vs j single
layers), not a full-chain rebuild.

THE MERGED OPERATOR.  Work over the SMALL (value-addressed) side, where the sieve
wiring is the pure value->value map floor(v/p) -- the cleanest face of the
floor(v/p) semigroup.  The homogeneous Lucy layer operator for prime p is
    A_p[S](v) = S(v) - g_p(v) * S(floor(v/p)),   g_p(v) = [p^2 <= v <= V].
(The public Lucy correction scalar sp_i = S_{i-1}(p_i - 1) = i-1 adds only a
VERIFIER-COMPUTABLE additive term per layer -- no S-claim -- so it does not touch the
2^j fill-in being measured; this homogeneous model isolates exactly that fill.)
Composing j layers (p_1 applied first .. p_j last) gives the inclusion-exclusion
expansion with 2^j subset terms (verified == direct composition, selftest 1):
    M[S](v) = sum_{T subset {1..j}} (-1)^|T| * GateProd_T(v) * S(floor(v / m_T)),
    m_T = prod_{k in T} p_k,
    GateProd_T(v) = prod_{k in T} g_{p_k}(floor(v / prod_{l in T, l>k} p_l)).

TWO SOUND REDUCTIONS of a claim  M~[S](r) = c  to a claim about the committed input
table S (both verified against ground truth, both cheat-tested):

 DIRECT.  The prover commits the 2^j routed tables ST_T[v] = S[floor(v/m_T)] and the
   2^j PUBLIC gate-product tables GP_T.  ONE degree-3 sum-check over the v-cube
   (nb rounds) reduces c to the 2^j openings {ST_T(r')}; each is routed back to a
   claim S(.) by its own degree-2 wiring sum-check, and the 2^j S-claims are folded
   (leaf_open.open_batch) to one residual closed at the base.  comm ~ 2^j * nb.

 BATCHED.  Stack the 2^j subset terms on a j-bit axis (the batched_trace pattern):
   ONE degree-3 sum-check over the (j+nb)-cube reduces c to a SINGLE opening
   STst(R), routed to one S-claim.  The 2^j lives ENTIRELY in verifier WORK
   (recomputing the public stacked operator MLE, O(2^j 2^nb)) and prover work -- NOT
   in comm.  comm ~ j + nb (polylog).

THE MEASUREMENT (per merged layer, attributed comm_main vs comm_fill; plus the
chain extrapolation total = ceil(K/j) * per-merged, and the prover fill-table size
2^j * 2^nb).  Note comm is PRIME-INDEPENDENT (the transcript size depends only on
nb, j -- primes only set gate/wiring VALUES; selftest 4 asserts this), so
total = ceil(K/j) * c_merged(j) is exact, not a representative-group estimate.

HEADLINE (FALSIFIABLE).  The chain TOTAL comm exponent in x stays ~0.5 for ALL j in
BOTH modes -- the sqrt x is LAYERING-INHERENT, on the verification face exactly as
open problem 2 found it on the prover face:
  * DIRECT: per-merged comm grows ~2^j, so total ceil(K/j)*2^j*nb has exponent 0.5
    (K = Theta(sqrt x)) and a prefactor 2^j/j that GROWS with j -- merging makes it
    worse, never better.
  * BATCHED: per-merged comm is polylog (j+nb), but total ceil(K/j)*(j+nb) -> K = sqrt x
    as j grows -- even the most aggressive comm-batching is FLOORED at sqrt x, because
    the merge-axis sum-check itself costs log(2^j)=j per merged layer and (K/j)*j = K.
  * The PROVER fill-table total ceil(K/j)*2^j*2^nb crosses from Theta(x) (small j)
    toward x^{3/2} as j -> log2 K (2^j ~ K) -- reproducing open problem 2's
    m*K^2/2 ~ x^{3/2} on the SAME construction.
So merging trades a bounded comm-prefactor win for a 2^j prover blow-up, and the comm
exponent never drops below 0.5.  That STRENGTHENS the S509 membership result: the
sieve certificate is Theta(sqrt x) at every merge depth -- the floor(v/p)-semigroup
obstruction is real on the verification face too.

WHAT WOULD FALSIFY THIS (any one):
  * the subset-sum operator disagreeing with direct j-layer composition;
  * an honest merged reduction (either mode) rejected, or its residual != true S~(r);
  * any cheating prover (wrong top claim / corrupted gate table / self-consistent
    false routing / wrong base residual) accepted above the field bound;
  * comm NOT prime-independent at fixed (nb, j);
  * DIRECT per-merged comm NOT growing ~2^j with j (e.g. staying linear);
  * BATCHED per-merged comm growing like 2^j (it should be polylog);
  * the chain TOTAL comm exponent dropping clearly below 0.5 for some fixed j in
    either mode -- THAT would be the genuine breakthrough S509 flagged (a real
    batching of the SEQUENTIAL floor(v/p) reductions); verify it against open
    problem 2's fill mechanism and re-derive the prover cost before believing it;
  * the prover fill-table total NOT trending toward x^{3/2} as j grows.

CLI:
  python3 merged_layer_comm.py --selftest
  python3 merged_layer_comm.py                      # profile (n=8..16, j=1,2,4) + jsweep
  python3 merged_layer_comm.py --field big          # over BIG_Q = 2^61-1
  python3 merged_layer_comm.py --jsweep --n 16      # fixed-n comm/prover vs j curve
"""
import argparse
import math
import time

import numpy as np

from compressed_prover_mult_trace import (DEFAULT_Q as Q, BIG_Q, SMALL_Q, FIELDS,
                                          _dt, _asum, eq_table, eq_point, mle_eval,
                                          sumcheck, lagrange_eval, fmul)
import compressed_layer as cl
import leaf_open as lo


# ----------------------------------------------------------------------
# the homogeneous merged operator (ground truth + the 2^j subset form)
# ----------------------------------------------------------------------

def to_field(int_tab, q):
    """Reduce an integer table (possibly negative) into a field array of dtype _dt."""
    a = np.asarray(int_tab, dtype=object) % q
    return a.astype(_dt(q))


def gate_table_int(p, V, nb):
    """g_p(v) = [p^2 <= v <= V] over v in [0, 2^nb), as a 0/1 int64 array."""
    N = 1 << nb
    g = np.zeros(N, dtype=np.int64)
    p2 = p * p
    hi = min(V, N - 1)
    if p2 <= hi:
        g[p2:hi + 1] = 1
    return g


def apply_layer_int(S, p, V):
    """homogeneous A_p[S](v) = S(v) - [p^2<=v<=V] S(floor(v/p)), integer table."""
    N = len(S)
    out = S.copy()
    p2 = p * p
    hi = min(V, N - 1)
    for v in range(p2, hi + 1):
        out[v] = S[v] - S[v // p]
    return out


def merged_direct_int(S0, primes, V):
    """Apply primes in order (primes[0] FIRST = innermost) -- the ground-truth
    merged table M[S0]."""
    S = list(S0)
    for p in primes:
        S = apply_layer_int(S, p, V)
    return S


def iter_subsets(primes):
    """Yield (mask, T_indices, sign, m_T) over all 2^j subsets, primes increasing."""
    j = len(primes)
    for mask in range(1 << j):
        T = [k for k in range(j) if (mask >> k) & 1]
        sign = -1 if (len(T) & 1) else 1
        mT = 1
        for k in T:
            mT *= primes[k]
        yield mask, T, sign, mT


def route_index(m, nb):
    """The wiring map floor(v/m) over v in [0, 2^nb)  (<= v < 2^nb, always valid)."""
    return np.arange(1 << nb, dtype=np.int64) // m


def gateprod_table_int(primes, T, V, nb):
    """GateProd_T(v) = prod_{k in T} g_{p_k}(floor(v / prod_{l in T,l>k} p_l)),
    0/1 int64 array over v."""
    N = 1 << nb
    GP = np.ones(N, dtype=np.int64)
    arange = np.arange(N, dtype=np.int64)
    for k in T:
        outer = 1
        for l in T:
            if l > k:
                outer *= primes[l]
        g = gate_table_int(primes[k], V, nb)
        GP = GP * g[arange // outer]
    return GP


def build_subsets(S0_int, primes, V, nb, q):
    """For every subset T: (mask, sign mod q, m_T, GP field table, ST field table).
    ST_T[v] = S0[floor(v/m_T)] (the routed input);  GP_T public 0/1."""
    out = []
    for mask, T, sign, mT in iter_subsets(primes):
        idx = route_index(mT, nb)
        ST_int = np.asarray(S0_int, dtype=object)[idx]            # S0[floor(v/m_T)]
        ST = (ST_int % q).astype(_dt(q))
        GP = to_field(gateprod_table_int(primes, T, V, nb), q)
        out.append((mask, sign % q, mT, GP, ST))
    return out


# ----------------------------------------------------------------------
# routing: reduce a claim about a routed table ST_T(rp) to a claim about S
# ----------------------------------------------------------------------

def _scatter_eq(idx, EQ, nb, q):
    """W[u] = sum_{v : idx[v]=u} EQ[v]  over u in [0,2^nb).  Exact mod q for both
    dtypes (object accumulates exactly; uint64 sums of < 2^nb terms < q stay below
    2^64 for q=Q, and we accumulate BIG_Q in object to avoid the 2^61 overflow)."""
    N = 1 << nb
    if EQ.dtype == object:
        acc = np.zeros(N, dtype=object)
        np.add.at(acc, idx, EQ)
        return (acc % q).astype(object)
    if q == BIG_Q:                                                # uint64 but > 2^31
        acc = np.zeros(N, dtype=object)
        np.add.at(acc, idx, EQ.astype(object))
        return (acc % q).astype(_dt(q))
    acc = np.zeros(N, dtype=_dt(q))
    np.add.at(acc, idx, EQ)
    return acc % q


def route_claim(S0, m, nb, rp, claimed, rng, stats, q):
    """Reduce  ST~(rp) == claimed  (ST[v]=S0[floor(v/m)]) to a claim S0~(rr)==sval.
    Identity: claimed = sum_u W[u] S0[u],  W[u] = sum_{v:floor(v/m)=u} eq~(rp,v).
    Degree-2 sum-check over u.  Returns (rr, sval, ok)."""
    t0 = time.perf_counter()
    EQrp = eq_table(rp, q)
    W = _scatter_eq(route_index(m, nb), EQrp, nb, q)
    _vw(stats, 1 << nb)                                           # verifier builds wiring
    stats["t_prover"] += time.perf_counter() - t0
    ok, rr, final, scal = sumcheck(claimed, {"W": W.copy(), "S": S0.copy()},
                                   [(1, ["W", "S"])], 2, rng, q)
    stats["comm"] += 3 * nb                                       # deg-2 round polys
    if not ok:
        return rr, None, False
    t0 = time.perf_counter()
    wv = mle_eval(W, rr, q)                                       # verifier grounds wiring
    _vw(stats, 1 << nb)
    ok = (int(scal["W"]) % q == wv
          and final % q == int(scal["W"]) * int(scal["S"]) % q)
    stats["t_verifier"] += time.perf_counter() - t0
    return rr, int(scal["S"]) % q, ok


# ----------------------------------------------------------------------
# DIRECT mode: one degree-3 sum-check over the v-cube, 2^j subset openings
# ----------------------------------------------------------------------

def reduce_direct(S0, subs, nb, r, claimed, rng, stats, q):
    """Reduce M~[S0](r)==claimed to a base close, via the 2^j subset terms summed
    in ONE sum-check, each opening routed back to S0 and the 2^j S-claims batched."""
    EQ = eq_table(r, q)
    tables = {"EQ": EQ.copy()}
    terms = []
    for (mask, sgn, mT, GP, ST) in subs:
        gn, sn = f"GP{mask}", f"ST{mask}"
        tables[gn] = GP.copy()
        tables[sn] = ST.copy()
        terms.append((sgn, ["EQ", gn, sn]))
    cmain0 = stats["comm"]
    ok, rp, final, scal = sumcheck(claimed, tables, terms, 3, rng, q)
    stats["comm"] += 4 * nb                                       # deg-3 round polys
    stats["comm_main"] = stats.get("comm_main", 0) + (stats["comm"] - cmain0)
    if not ok:
        return False
    t0 = time.perf_counter()
    eqv = eq_point(r, rp, q)
    if int(scal["EQ"]) % q != eqv:
        return False
    acc = 0
    for (mask, sgn, mT, GP, ST) in subs:
        gpv = int(scal[f"GP{mask}"]) % q
        _vw(stats, 1 << nb)                                       # verifier grounds public gate
        if gpv != mle_eval(GP, rp, q):
            stats["t_verifier"] += time.perf_counter() - t0
            return False
        stv = int(scal[f"ST{mask}"]) % q
        acc = (acc + sgn * eqv % q * gpv % q * stv) % q
    ok = (final % q == acc)
    stats["t_verifier"] += time.perf_counter() - t0
    stats["comm"] += len(subs)                                    # the 2^j ST openings
    if not ok:
        return False
    # route every ST_T(rp) -> a claim about S0, then batch
    sclaims = []
    for (mask, sgn, mT, GP, ST) in subs:
        stv = int(scal[f"ST{mask}"]) % q
        rr, sval, okr = route_claim(S0, mT, nb, rp, stv, rng, stats, q)
        if not okr:
            return False
        sclaims.append((rr, sval))
    stats["comm"] += len(subs)                                    # routing residual scalars
    res = lo.open_batch(S0, sclaims, rng, stats, q)               # comm 1+3nb (or 3nb if k=1)
    stats["comm"] += 1                                            # base residual scalar
    _vw(stats, len(S0))                                           # the one-time base close
    return bool(res["ok"] and lo.close_eval(S0, res["r"], res["residual"], q))


# ----------------------------------------------------------------------
# BATCHED mode: one degree-3 sum-check over the (j+nb)-cube, ONE opening
# ----------------------------------------------------------------------

def reduce_batched(S0, subs, primes, nb, r, claimed, rng, stats, q):
    """Stack the 2^j subset terms on a j-bit axis (MSB) and discharge in ONE
    degree-3 sum-check over the (j+nb)-cube -> ONE opening STst(R), routed to one
    S0 claim.  The 2^j lives only in verifier/prover WORK, not comm."""
    j = len(primes)
    Lk = j                                                        # 2^j subsets exactly
    N = 1 << nb
    M = 1 << Lk
    assert M == len(subs)
    t0 = time.perf_counter()
    EQr = eq_table(r, q)
    A = np.zeros(M * N, dtype=_dt(q))                             # A[mask,v]=sgn(mask)eq(r,v)
    GPst = np.zeros(M * N, dtype=_dt(q))
    STst = np.zeros(M * N, dtype=_dt(q))
    mlist = [0] * M
    for (mask, sgn, mT, GP, ST) in subs:
        b = mask * N
        A[b:b + N] = fmul(EQr, sgn, q)
        GPst[b:b + N] = GP
        STst[b:b + N] = ST
        mlist[mask] = mT
    stats["t_prover"] += time.perf_counter() - t0
    cmain0 = stats["comm"]
    ok, R, final, scal = sumcheck(claimed,
                                  {"A": A.copy(), "GP": GPst.copy(), "ST": STst.copy()},
                                  [(1, ["A", "GP", "ST"])], 3, rng, q)
    stats["comm"] += 4 * (Lk + nb)                               # deg-3 round polys
    stats["comm_main"] = stats.get("comm_main", 0) + (stats["comm"] - cmain0)
    if not ok:
        return False
    t0 = time.perf_counter()
    RT, Rv = R[:Lk], R[Lk:]
    sgn_R = 1                                                     # sgn~(R_T)=prod(1-2R_k)
    for rk in RT:
        sgn_R = sgn_R * ((1 - 2 * int(rk)) % q) % q
    Av = sgn_R * eq_point(r, Rv, q) % q                          # A public -> verifier
    gpv = mle_eval(GPst, R, q)                                   # GP public -> verifier
    _vw(stats, M * N)                                            # stacked-operator eval (2^j 2^nb)
    if int(scal["A"]) % q != Av or int(scal["GP"]) % q != gpv:
        stats["t_verifier"] += time.perf_counter() - t0
        return False
    stv = int(scal["ST"]) % q
    ok = (final % q == Av * gpv % q * stv % q)
    stats["t_verifier"] += time.perf_counter() - t0
    stats["comm"] += 1                                            # the one ST opening
    if not ok:
        return False
    # route STst(R) -> one claim about S0:  STst~(R) = sum_u Wbig[u] S0[u]
    t0 = time.perf_counter()
    EQR = eq_table(R, q)
    _vw(stats, M * N)                                            # stacked wiring build (2^j 2^nb)
    Wbig = np.zeros(N, dtype=object if EQR.dtype == object or q == BIG_Q else _dt(q))
    ar = np.arange(N, dtype=np.int64)
    for mask in range(M):
        blk = EQR[mask * N:(mask + 1) * N]
        idx = ar // mlist[mask]
        if Wbig.dtype == object:
            np.add.at(Wbig, idx, (blk.astype(object) if blk.dtype != object else blk))
        else:
            np.add.at(Wbig, idx, blk)
    Wbig = (Wbig % q)
    if Wbig.dtype == object and q != BIG_Q:
        Wbig = Wbig.astype(_dt(q))
    elif Wbig.dtype == object and q == BIG_Q:
        Wbig = Wbig.astype(_dt(q))
    stats["t_prover"] += time.perf_counter() - t0
    okr, rr, finalr, scalr = sumcheck(stv, {"W": Wbig.copy(), "S": S0.copy()},
                                      [(1, ["W", "S"])], 2, rng, q)
    stats["comm"] += 3 * nb
    if not okr:
        return False
    t0 = time.perf_counter()
    wv = mle_eval(Wbig, rr, q)
    _vw(stats, 1 << nb)
    okr = (int(scalr["W"]) % q == wv
           and finalr % q == int(scalr["W"]) * int(scalr["S"]) % q)
    stats["t_verifier"] += time.perf_counter() - t0
    stats["comm"] += 2                                            # routing + base residual scalars
    if not okr:
        return False
    _vw(stats, len(S0))                                           # the one-time base close
    return bool(lo.close_eval(S0, rr, int(scalr["S"]) % q, q))


# ----------------------------------------------------------------------
# the measured run of one merged layer
# ----------------------------------------------------------------------

def _stats():
    return {"t_prover": 0.0, "t_verifier": 0.0, "comm": 0, "comm_main": 0,
            "vwork": 0}


def _vw(stats, size):
    """Account ONE verifier large-table evaluation (mle_eval / wiring scatter) of
    O(size) field ops -- the same metric S507 used for the chain.  In BOTH modes
    the verifier must touch the routed/stacked operator, so vwork carries the 2^j
    fill exactly as the prover does: batching the COMM does not remove the 2^j,
    it only moves it out of the transcript and into the WORK."""
    stats["vwork"] = stats.get("vwork", 0) + int(size)


def run_merged(x, j, mode, rng, q=Q, cheat=None, start=0):
    """Build the depth-j merged layer over primes[start:start+j], reduce a claim
    M~[S0](r)=c, and return (ok, stats, info).  cheat in {None,'c_value','c_gate',
    'c_route'}.  Returns prover fill-table size etc. in info."""
    V = cl.isqrt(x)
    nb = max(1, V.bit_length())
    primes = cl.compressed_lucy(x)[0]
    assert start + j <= len(primes), (x, j, start, len(primes))
    P = primes[start:start + j]
    S0_int = list(cl.compressed_lucy(x)[1][0])                    # real S_0 small side
    S0_int += [0] * ((1 << nb) - len(S0_int))
    S0 = to_field(S0_int, q)

    M_int = merged_direct_int(S0_int, P, V)                      # ground truth
    M_field = to_field(M_int, q)
    subs = build_subsets(S0_int, P, V, nb, q)

    # cheat: corrupt one subset's routed table (self-consistent liar recomputes c)
    if cheat == "c_route":
        # corrupt mask=1 (single smallest prime, gate active) ST at one entry
        subs = [list(t) for t in subs]
        ST = subs[1][4].copy()
        idx = (1 << nb) // 2
        ST[idx] = (int(ST[idx]) + 1) % q
        subs[1][4] = ST
        subs = [tuple(t) for t in subs]
        # rebuild the (corrupted) merged table so the top claim is self-consistent
        Mc = np.zeros(1 << nb, dtype=object)
        for v in range(1 << nb):
            tot = 0
            for (mask, sgn, mT, GP, STt) in subs:
                s = 1 if sgn == 1 else -1
                tot += s * int(GP[v]) * int(STt[v])
            Mc[v] = tot % q
        M_field = Mc.astype(_dt(q))
    elif cheat == "c_gate":
        subs = [list(t) for t in subs]
        GP = subs[1][3].copy()
        idx = (1 << nb) // 2
        GP[idx] = (int(GP[idx]) + 1) % q                          # corrupt a public gate
        subs[1][3] = GP
        subs = [tuple(t) for t in subs]

    r = [int(rng.integers(0, q)) for _ in range(nb)]
    claimed = mle_eval(M_field, r, q)
    if cheat == "c_value":
        claimed = (claimed + 1) % q

    st = _stats()
    if mode == "direct":
        ok = reduce_direct(S0, subs, nb, r, claimed, rng, st, q)
    else:
        ok = reduce_batched(S0, subs, P, nb, r, claimed, rng, st, q)
    st["comm_fill"] = st["comm"] - st["comm_main"]
    info = dict(V=V, nb=nb, primes=P, prover_tab=(1 << j) * (1 << nb))
    return ok, st, info


# ----------------------------------------------------------------------
# measurement helpers
# ----------------------------------------------------------------------

def chain_K(x):
    return len(cl.compressed_lucy(x)[0])


def measure_point(x, j, mode, q, seed=1):
    """Honest merged-layer run -> per-merged comm + the chain extrapolation."""
    ok, st, info = run_merged(x, j, mode, np.random.default_rng(seed), q)
    assert ok, ("honest merged reduction rejected", x, j, mode)
    K = chain_K(x)
    groups = math.ceil(K / j)
    return dict(
        x=x, j=j, mode=mode, nb=info["nb"], K=K, groups=groups,
        per_comm=st["comm"], comm_main=st["comm_main"], comm_fill=st["comm_fill"],
        total_comm=groups * st["comm"],
        vwork=st["vwork"], total_vwork=groups * st["vwork"],
        prover_tab=info["prover_tab"],
        total_prover=groups * info["prover_tab"],
    )


def _fit(ns, vals):
    return cl._fit_slope(list(ns), [math.log2(max(1, v)) for v in vals])


# ----------------------------------------------------------------------
# the profile: total comm exponent ~0.5 for all j, both modes
# ----------------------------------------------------------------------

def profile(ns=(8, 10, 12, 14, 16), js=(1, 2, 4), field="q", seed=1):
    q = FIELDS[field]
    print(f"=== merged-layer comm probe, field={field} ===")
    print("homogeneous small-side operator A_p[S](v)=S(v)-[p^2<=v<=V]S(v//p); "
          "merge j consecutive primes")
    print(f"x = 2^n - 1; K = pi(sqrt x); chain total = ceil(K/j) * per-merged\n")

    data = {}
    for mode in ("direct", "batched"):
        print(f"-- mode = {mode} --")
        print(f"{'n':>3} {'nb':>3} {'K':>4} | " + " ".join(
            f"j={j}:{'per':>5}/{'main':>4}/{'fill':>5}|{'TOTcomm':>8}" for j in js))
        for n in ns:
            x = (1 << n) - 1
            row = []
            for j in js:
                m = measure_point(x, j, mode, q, seed)
                data[(mode, n, j)] = m
                row.append(f"j={j}:{m['per_comm']:>5}/{m['comm_main']:>4}/"
                           f"{m['comm_fill']:>5}|{m['total_comm']:>8}")
            nb = data[(mode, n, js[0])]["nb"]
            K = data[(mode, n, js[0])]["K"]
            print(f"{n:>3} {nb:>3} {K:>4} | " + " ".join(row))
        print()

    # ---- exponent fits: TOTAL comm vs x, per j, per mode (the headline) ----
    last = ns[-4:] if len(ns) >= 4 else ns
    print(f"-- fitted exponents alpha (metric ~ x^alpha; log2 metric vs n, "
          f"n in {tuple(last)}) --")
    print(f"{'mode':>8} {'j':>2} {'TOTcomm_a':>10} {'TOTvwork_a':>11} "
          f"{'TOTprover_a':>12}   expectation")
    fits = {}
    for mode in ("direct", "batched"):
        for j in js:
            tc = _fit(last, [data[(mode, n, j)]["total_comm"] for n in last])
            tv = _fit(last, [data[(mode, n, j)]["total_vwork"] for n in last])
            tp = _fit(last, [data[(mode, n, j)]["total_prover"] for n in last])
            fits[(mode, j)] = dict(tc=tc, tv=tv, tp=tp)
            print(f"{mode:>8} {j:>2} {tc:>10.3f} {tv:>11.3f} {tp:>12.3f}   "
                  f"{'TOTcomm ~0.5 (sqrt x) for ALL j'}")
    print()

    # ---- per-merged growth with j: fill (DIRECT ~2^j) + vwork (BOTH ~2^j) ----
    nstar = ns[-1]
    print(f"-- per-merged growth vs merge depth j at n={nstar}: comm_FILL (DIRECT "
          f"tracks 2^j; BATCHED flat),\n   and vwork (BOTH track 2^j -- batching "
          f"moves the 2^j out of comm but NOT out of the WORK) --")
    print(f"{'j':>2} {'2^j':>5} | {'dir_fill':>8} {'fill/2^j':>8} {'dir_vwork':>9} | "
          f"{'bat_fill':>8} {'bat_vwork':>9}")
    f1 = data[("direct", nstar, js[0])]["comm_fill"]
    for j in js:
        d = data[("direct", nstar, j)]
        b = data[("batched", nstar, j)]
        print(f"{j:>2} {1 << j:>5} | {d['comm_fill']:>8} "
              f"{d['comm_fill'] / (1 << j):>8.1f} {d['vwork']:>9} | "
              f"{b['comm_fill']:>8} {b['vwork']:>9}")
    print()

    # ---- falsifiable assertions ----
    for mode in ("direct", "batched"):
        for j in js:
            a = fits[(mode, j)]["tc"]
            assert 0.40 <= a <= 0.65, \
                (mode, j, "TOTAL comm exponent not ~sqrt x", a)
    # DIRECT comm_FILL grows ~2^j (the inclusion-exclusion fill-in):
    if len(js) >= 2 and js[-1] >= 2:
        f_lo = data[("direct", nstar, js[0])]["comm_fill"]
        f_hi = data[("direct", nstar, js[-1])]["comm_fill"]
        ratio = f_hi / f_lo
        assert ratio > js[-1] / js[0], \
            ("DIRECT comm_fill not growing faster than linear in j", ratio)
        # BATCHED comm fill is FLAT in j (the 2^j left the transcript)
        b_lo = data[("batched", nstar, js[0])]["comm_fill"]
        b_hi = data[("batched", nstar, js[-1])]["comm_fill"]
        assert b_hi <= b_lo + 4, ("BATCHED comm_fill not flat in j", b_lo, b_hi)
        # ...but BATCHED vwork still carries the 2^j (work, not comm)
        vw_lo = data[("batched", nstar, js[0])]["vwork"]
        vw_hi = data[("batched", nstar, js[-1])]["vwork"]
        assert vw_hi / vw_lo > js[-1] / js[0], \
            ("BATCHED vwork escaped the 2^j -- would be the real win", vw_hi / vw_lo)
    # PROVER fill total trends UP with j (toward x^{3/2}); strictly above comm exp.
    tp_dir = [fits[("direct", j)]["tp"] for j in js]
    assert tp_dir[-1] >= tp_dir[0] - 0.05, \
        ("prover-table exponent not non-decreasing in j", tp_dir)
    assert tp_dir[-1] > 0.65, ("prover-table exponent not >> comm (toward x^1.5)", tp_dir[-1])

    print("HEADLINE: chain TOTAL comm exponent ~0.5 (Theta(sqrt x)) for EVERY merge "
          "depth j in BOTH\n  modes -> the sqrt x certificate size is LAYERING-INHERENT. "
          "DIRECT comm_fill grows ~2^j\n  (the inclusion-exclusion fill); BATCHED pushes "
          "the 2^j OUT of comm (fill flat) into\n  verifier WORK (vwork exponent rises "
          "with j, matching the prover's), comm still floored\n  at ceil(K/j)*(j+nb) -> "
          "K = Theta(sqrt x).  The PROVER fill-table total trends toward\n  x^{3/2} as j "
          "grows -- open problem 2's m*K^2/2, now shown on the VERIFICATION face.")
    return data, fits


def jsweep(n=16, field="q", seed=1):
    """Fixed n, j = 1..jmax: the comm/prover-vs-j curve.  Shows comm has a (shallow)
    minimum but never drops below ~K=sqrt x while the prover fill blows up 2^j."""
    q = FIELDS[field]
    x = (1 << n) - 1
    K = chain_K(x)
    nb = max(1, cl.isqrt(x).bit_length())
    jmax = min(6, len(cl.compressed_lucy(x)[0]))
    print(f"=== j-sweep at n={n} (x={x}, sqrt x ~ 2^{nb}={1 << nb}, K=pi(sqrt x)={K}), "
          f"field={field} ===")
    print("merge depth j vs: chain TOTAL comm = ceil(K/j)*per-merged, verifier WORK "
          "total, PROVER fill-table total = ceil(K/j)*2^j*2^nb")
    print(f"{'j':>2} {'groups':>6} | {'dir_TOTAL':>10} {'bat_TOTAL':>10} | "
          f"{'bat_vwkTOT':>11} | {'prover_TOT':>12} {'/x':>7}")
    for j in range(1, jmax + 1):
        groups = math.ceil(K / j)
        md = measure_point(x, j, "direct", q, seed)
        mb = measure_point(x, j, "batched", q, seed)
        ptot = md["total_prover"]
        print(f"{j:>2} {groups:>6} | {md['total_comm']:>10} {mb['total_comm']:>10} | "
              f"{mb['total_vwork']:>11} | {ptot:>12} {ptot / x:>7.2f}")
    print(f"\nfloor: chain total comm cannot drop below ~K = {K} (the merge-axis "
          f"sum-check costs\n  log2(2^j)=j per merged layer, and ceil(K/j)*j = K). "
          f"DIRECT total bottoms out then\n  GROWS (2^j fill); BATCHED total keeps "
          f"dropping by a constant factor toward the K floor,\n  BUT its verifier WORK "
          f"total and the prover_TOT/x both grow ~2^j/j -- reaching x^{{1/2}} extra\n  "
          f"(=> x^{{3/2}} total) as 2^j -> K.  Batching moves the 2^j out of comm, never "
          f"out of the WORK\n  (open problem 2, both faces).  No j drops the comm "
          f"exponent below 0.5.")


# ----------------------------------------------------------------------
# selftest
# ----------------------------------------------------------------------

def selftest():
    print("[selftest] merged_layer_comm (S510)")
    rng = np.random.default_rng

    # 1. the 2^j subset-sum operator == direct j-layer composition (exact ints),
    #    AND its field MLE == the subset-sum's field MLE at random points.
    for n in (8, 10, 12, 14):
        x = (1 << n) - 1
        V = cl.isqrt(x); nb = max(1, V.bit_length())
        S0 = list(cl.compressed_lucy(x)[1][0]); S0 += [0] * ((1 << nb) - len(S0))
        primes = cl.compressed_lucy(x)[0]
        for j in (1, 2, 3, 4):
            if j > len(primes):
                continue
            P = primes[:j]
            d = merged_direct_int(S0[:], P, V)
            # subset sum (int)
            s = [0] * (1 << nb)
            for v in range(1 << nb):
                tot = 0
                for mask, T, sign, mT in iter_subsets(P):
                    gp = 1
                    for k in T:
                        outer = 1
                        for l in T:
                            if l > k:
                                outer *= P[l]
                        gp *= int(gate_table_int(P[k], V, nb)[v // outer])
                    tot += sign * gp * S0[v // mT]
                s[v] = tot
            assert d == s, (n, j, "subset-sum != direct composition")
    print("  1. inclusion-exclusion subset operator == direct composition OK")

    # 2. honest merged reduction ACCEPTS and its base residual is the true S~(r),
    #    over q AND BIG_Q, both modes, j in {1,2,3,4}.
    for q in (Q, BIG_Q):
        for n in (8, 10, 12):
            x = (1 << n) - 1
            primes = cl.compressed_lucy(x)[0]
            for j in (1, 2, 3, 4):
                if j > len(primes):
                    continue
                for mode in ("direct", "batched"):
                    ok, st, info = run_merged(x, j, mode, rng(1), q)
                    assert ok, (q, n, j, mode, "honest rejected")
                    assert st["comm_main"] > 0 and st["comm_fill"] >= 0
    print("  2. honest merged reduction accepts (q & BIG_Q, both modes, j<=4) OK")

    # 3. cheating provers REJECT in both modes, over q and BIG_Q.
    for q in (Q, BIG_Q):
        for n in (10, 12):
            x = (1 << n) - 1
            for j in (2, 3):
                for mode in ("direct", "batched"):
                    for cheat in ("c_value", "c_gate", "c_route"):
                        ok, st, info = run_merged(x, j, mode, rng(2), q, cheat=cheat)
                        assert not ok, (q, n, j, mode, cheat, "cheat accepted")
    print("  3. cheats reject: c_value / c_gate / c_route (q & BIG_Q, both modes) OK")

    # 4. comm is PRIME-INDEPENDENT: different prime windows (same nb,j) -> same comm;
    #    and field-independent (q vs BIG_Q identical comm counts).
    for n in (12, 14):
        x = (1 << n) - 1
        primes = cl.compressed_lucy(x)[0]
        for j in (2, 3):
            for mode in ("direct", "batched"):
                _, s0, _ = run_merged(x, j, mode, rng(3), Q, start=0)
                _, s1, _ = run_merged(x, j, mode, rng(7), Q, start=1)
                assert s0["comm"] == s1["comm"], (n, j, mode, "comm prime-dependent")
                _, sb, _ = run_merged(x, j, mode, rng(3), BIG_Q, start=0)
                assert s0["comm"] == sb["comm"], (n, j, mode, "comm field-dependent")
    print("  4. comm prime-independent AND field-independent OK")

    # 5. scaling sanity (short sweep): TOTAL comm exponent in [0.4,0.65] for all j,
    #    both modes; DIRECT per-merged comm grows faster than linear in j; BATCHED
    #    per-merged < DIRECT at j=4; prover-table exponent clearly above comm.
    ns = (8, 10, 12, 14)
    js = (1, 2, 4)
    cache = {}
    for mode in ("direct", "batched"):
        for n in ns:
            for j in js:
                cache[(mode, n, j)] = measure_point((1 << n) - 1, j, mode, Q, 1)
    for mode in ("direct", "batched"):
        for j in js:
            a = _fit(ns, [cache[(mode, n, j)]["total_comm"] for n in ns])
            assert 0.40 <= a <= 0.70, (mode, j, "total comm exp", a)
    # DIRECT comm_fill ~2^j; BATCHED comm_fill flat but its vwork still ~2^j.
    df1 = cache[("direct", ns[-1], 1)]["comm_fill"]
    df4 = cache[("direct", ns[-1], 4)]["comm_fill"]
    assert df4 / df1 > 4, ("direct comm_fill not ~2^j in j", df4 / df1)
    bf1 = cache[("batched", ns[-1], 1)]["comm_fill"]
    bf4 = cache[("batched", ns[-1], 4)]["comm_fill"]
    assert bf4 <= bf1 + 4, ("batched comm_fill not flat", bf1, bf4)
    bv1 = cache[("batched", ns[-1], 1)]["vwork"]
    bv4 = cache[("batched", ns[-1], 4)]["vwork"]
    assert bv4 / bv1 > 4, ("batched vwork escaped the 2^j", bv4 / bv1)
    ap = _fit(ns, [cache[("direct", n, 4)]["total_prover"] for n in ns])
    ac = _fit(ns, [cache[("direct", n, 4)]["total_comm"] for n in ns])
    assert ap > ac + 0.2, ("prover-table exponent not above comm", ap, ac)
    print(f"  5. scaling OK (total comm exp ~0.5 all j; direct comm_fill "
          f"df4/df1={df4/df1:.1f}~2^j; batched fill flat ({bf1}->{bf4}) but vwork "
          f"bv4/bv1={bv4/bv1:.1f}~2^j; prover exp {ap:.2f} > comm exp {ac:.2f})")

    print("[selftest] PASS")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--selftest", action="store_true")
    ap.add_argument("--jsweep", action="store_true")
    ap.add_argument("--n", type=int, default=16, help="n for --jsweep")
    ap.add_argument("--nmax", type=int, default=16, help="largest n for --profile")
    ap.add_argument("--field", choices=list(FIELDS), default="q")
    ap.add_argument("--seed", type=int, default=1)
    args = ap.parse_args()
    if args.selftest:
        selftest()
    elif args.jsweep:
        jsweep(n=args.n, field=args.field, seed=args.seed)
    else:
        ns = tuple(range(8, args.nmax + 1, 2))
        profile(ns, field=args.field, seed=args.seed)
        print()
        jsweep(n=args.nmax, field=args.field, seed=args.seed)


if __name__ == "__main__":
    main()
