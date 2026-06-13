#!/usr/bin/env python3
"""
lucy_dp_delegated_wiring.py — S491 continuation: the delegated-wiring
upgrade of the exact-pi(x) interactive proof.

The base protocol (lucy_dp_verification.py) has the verifier evaluate
the division-wiring MLE W~(r_v, r_u) = MLE of [u = floor(v/p)] by a
p-state automaton DP: O(n*p) per layer, total O~(x/log x) over all
layers (and the width-dichotomy measurement proved 2p+1 is TIGHT for
explicit automatons). This script removes that bottleneck by DELEGATING
the automaton DP itself to the prover:

  The DP is a chain v_{j+1} = M_j v_j of n matrix-vector products over
  the 2^l-embedded state space (l = ceil(log2 p)). Each M_j is known to
  the verifier symbolically: its nonzero pattern is the AFFINE relation
    s' = 2 s + a - b p,  (a, b) in {0,1}^2,
  weighted by the outer challenge weights for bit a of v and quotient
  bit b of u. The agreed low-degree extension is
    M~_j(sig', sig) = LT~(sig') * LT~(sig) * R~_j(sig', sig),
  R~_j = sum_{a,b} wv_j(a) wu_j(b) * AFF~_{a-bp}(sig, sig'),
  and the verifier evaluates AFF~_c (the MLE of [y = 2x + c] on l-bit
  integers) in O(l) by a WIDTH-4 carry automaton — no p-dependence.
  Each chain layer is reduced by a 2l-variable sum-check with the
  eq-selector inserted (standard GKR form, needed because the layer
  identity only extends multilinearly through Boolean points).

Verifier wiring cost per Lucy layer: O(n * log p) instead of O(n * p).
Total exact-pi(x) verifier: sum_{p <= sqrt x} O(n log p) = O~(sqrt x).
Soundness remains unconditional. Trade-off (measured): communication
rises (the chain transcripts), prover gains O(n p^2)-ish extra work —
both still O~(sqrt x)-compatible overall.

What would falsify: honest delegated run rejected or accepted with a
count differing from the sieve; cheating prover (false wiring claim, or
corrupted DP propagated consistently) accepted above the soundness
bound; wiring-check verifier time failing to flatten in p (must scale
~log p while the automaton baseline scales ~p).

Usage:
  python3 lucy_dp_delegated_wiring.py --selftest
  python3 lucy_dp_delegated_wiring.py --n 16 --cheat-trials 10
  python3 lucy_dp_delegated_wiring.py --wiring-bench
"""

import argparse
import time

import numpy as np

# Field-parameterised helpers (q kwarg, _dt/_asum dtype machinery) come from
# compressed_prover_mult_trace; they default to q=Q (= DEFAULT_Q = 2^31-1), so
# the 2^n-table demo below (run_protocol_delegated) and every existing caller
# stay BIT-IDENTICAL while the compressed chain can now run over a larger prime
# (BIG_Q = 2^61-1, object dtype) -- the field lift threaded end-to-end.
from compressed_prover_mult_trace import (DEFAULT_Q as Q, BIG_Q, _dt, _asum,
                                          lagrange_eval, eq_table, eq_point,
                                          mle_eval, sumcheck, fmul, fast_big)
from lucy_dp_verification import (s0_closed_form, ge_const_eval, w_div_eval,
                                  prover_tables, direct_pi)


# ----------------------------------------------------------------------
# verifier-side O(l) primitives
# ----------------------------------------------------------------------

def lt_const_eval(rho, M, l, q=Q):
    """MLE of [s < M] on l-bit s at point rho. O(l). When M >= 2^l the
    predicate is constantly true (the p = 2, l = 1 boundary case)."""
    if M >= (1 << l):
        return 1
    return (q + 1 - ge_const_eval(rho, M, l, q)) % q


def aff_rel_eval(rho_x, rho_y, c, l, q=Q):
    """MLE of the affine relation [y = 2x + c] over l-bit integers x, y,
    at field points rho_x, rho_y (bits MSB-first). Width-4 carry
    automaton (carry, prev-x-bit), processed LSB-first. O(l).

    c >= 0: y = 2x + c. At LSB-position i: t = x_{i-1} + c_i + carry;
            requires y_i = t & 1; accept iff final x_prev = carry = 0.
    c < 0:  y + d = 2x, d = -c > 0. t = y_i + d_i + carry; requires
            (2x)_i = x_{i-1} = t & 1; accept iff final carry == x_prev.
    """
    xb = list(reversed([int(v) % q for v in rho_x]))
    yb = list(reversed([int(v) % q for v in rho_y]))
    # process enough positions to exhaust the constant and the shift;
    # bits of x, y beyond position l-1 are virtual zeros
    L = max(l, abs(c).bit_length()) + 1
    st = {(0, 0): 1}  # state: (carry, xprev) -> weight
    for i in range(L):
        wx = [(q + 1 - xb[i]) % q, xb[i]] if i < l else [1, 0]
        wy = [(q + 1 - yb[i]) % q, yb[i]] if i < l else [1, 0]
        new = {}
        for (carry, xprev), w in st.items():
            if w == 0:
                continue
            if c >= 0:
                ci = (c >> i) & 1
                t = xprev + ci + carry
                ybit, carry2 = t & 1, t >> 1
                if wy[ybit] == 0:
                    continue
                for xi in (0, 1):
                    if wx[xi] == 0:
                        continue
                    key = (carry2, xi)
                    new[key] = (new.get(key, 0)
                                + w * wy[ybit] % q * wx[xi]) % q
            else:
                d = -c
                di = (d >> i) & 1
                for yi in (0, 1):
                    if wy[yi] == 0:
                        continue
                    t = yi + di + carry
                    if (t & 1) != xprev:
                        continue
                    carry2 = t >> 1
                    for xi in (0, 1):
                        if wx[xi] == 0:
                            continue
                        key = (carry2, xi)
                        new[key] = (new.get(key, 0)
                                    + w * wy[yi] % q * wx[xi]) % q
        st = new
    # all bit positions through L-1 are consumed (x_{L-1} = 0 virtual),
    # so the equation holds iff no carry remains
    total = 0
    for (carry, xprev), w in st.items():
        if carry == 0:
            total = (total + w) % q
    return total


def r_tilde_eval(rho_in, rho_out, p, l, wv, wu, q=Q):
    """Agreed extension R~_j(rho_out, rho_in) =
    sum_{a,b} wv[a] wu[b] AFF~_{a - b p}(x=rho_in, y=rho_out). O(l)."""
    tot = 0
    for a in (0, 1):
        for b in (0, 1):
            tot = (tot + wv[a] * wu[b] % q
                   * aff_rel_eval(rho_in, rho_out, a - b * p, l, q)) % q
    return tot


# ----------------------------------------------------------------------
# inner protocol: delegated evaluation of W~(r_v, r_u)
# ----------------------------------------------------------------------

def bit_weights(r, q=Q):
    r = int(r) % q
    return [(q + 1 - r) % q, r]


def inner_chain_vectors(r_v, r_u, p, l, n, q=Q):
    """Honest prover: chain v_0..v_n over the 2^l-embedded state space.
    Support is kept inside the valid states [0, p) — transitions into
    junk states [p, 2^l) must be excluded so the layer identity
    v_{j+1}(s') = LT(s') * sum_s R(s', s) LT(s) v_j(s) holds on the
    WHOLE Boolean cube (the agreed extension kills junk via LT)."""
    dim = 1 << l
    vs = [np.zeros(dim, dtype=_dt(q))]
    vs[0][0] = 1
    for j in range(n):
        wv, wu = bit_weights(r_v[j], q), bit_weights(r_u[j], q)
        v = vs[-1]
        nv = np.zeros(dim, dtype=_dt(q))
        for s in range(p):
            if v[s] == 0:
                continue
            for a in (0, 1):
                for b in (0, 1):
                    y = 2 * s + a - b * p
                    if 0 <= y < p:
                        nv[y] = (int(nv[y]) + int(v[s]) * wv[a] % q
                                 * wu[b]) % q
        vs.append(nv)
    return vs


def inner_r_table(p, l, wv, wu, q=Q):
    """Prover table of R~_j on the Boolean (sig', sig) cube,
    index = (sig' << l) | sig."""
    dim = 1 << l
    Rt = np.zeros(dim * dim, dtype=_dt(q))
    for s in range(dim):
        for a in (0, 1):
            for b in (0, 1):
                y = 2 * s + a - b * p
                if 0 <= y < dim:
                    Rt[(y << l) | s] = (int(Rt[(y << l) | s])
                                        + wv[a] * wu[b]) % q
    return Rt


# ----------------------------------------------------------------------
# STRUCTURED chain-layer reduction (S496): O(2^l)=O(p) prover, NOT O(p^2)
# ----------------------------------------------------------------------
#
# The dense backward sweep materializes five 2^{2l}-size tables per chain
# layer (np.repeat/np.tile/inner_r_table) and runs a 2l-variable degree-5
# sum-check over them -> O(p^2) per chain layer -> Sum_p n p^2 ~ O~(x^{3/2}),
# reintroducing the x^{3/2} the compressed DP removed.  But the summand
#   F(sig',sig) = E(sig') P(sig') R_j(sig',sig) L(sig) v_j(sig)
# factors: E,P depend only on sig'; L,v depend only on sig; R_j is the
# SPARSE width-4 affine operator (4 nonzeros per column sig, via
# sig'=2sig+a-bp) and v_j has support <= p.  Splitting the 2l-variable
# sum-check at the sig'/sig boundary lets BOTH phases run over O(p)-size
# tables while producing the IDENTICAL transcript (same round polynomials,
# challenges, final claim, scal['V']):
#
#   phase 1  (bind sig' = y):  product  E(y) P(y) G(y),  where
#            G(y) = Sum_sig R_j(y,sig) L(sig) v_j(sig)  is the exact
#            multilinear y-table of the folded sig-sum -- built in O(p)
#            because each sig<p sends weight to <=4 targets y=2sig+a-bp.
#   phase 2  (bind sig = w):   product  (kappa R*(w)) L(w) v_j(w),  where
#            R*(w) = R_j(rho_out, w),  kappa = E(rho_out) P(rho_out),
#            both from O(p) work over eq_table(rho_out).
#
# Why the transcript is bit-identical: (i) sum-check challenges are drawn
# from rng, NEVER from the eval values, and both the dense (one 2l-round
# call) and the split (l+l rounds) draw exactly 2l challenges in the same
# order from the same rng; (ii) the round polynomial is the unique
# multilinear sum-check polynomial of the SAME function F -- G just regroups
# the sig-sum, which commutes with the multilinear extension in y -- so the
# threaded claim, final value and scal['V'] coincide.  The selftest checks
# round-poly equality exhaustively for small p.

def _modmul_arr(arr, scalar, q=Q):
    """(arr * scalar) % q for a field array (entries < q) and a scalar < q.
    BIG_Q routes through fmul (uint64 Mersenne under FAST_BIG, else object %q).
    For q <= 2^31-1 the scalar must be a np.uint64 so numpy keeps a uint64*uint64
    product (Python-int scalars would promote to float64); object arrays multiply
    by a plain Python int, exact at any q."""
    s = int(scalar) % q
    if q == BIG_Q:
        return fmul(arr, s, q)
    if arr.dtype == object:
        return (arr * s) % q
    return (arr * np.uint64(s)) % q


def _sumcheck_rec(claim, tables, terms, degree, rng, rec):
    """Exact copy of lucy_dp_verification.sumcheck that also appends each
    round's eval list to rec (test instrumentation only)."""
    nvars = int(np.log2(len(next(iter(tables.values())))))
    rs = []
    c = int(claim) % Q
    for _ in range(nvars):
        evals = []
        for X in range(degree + 1):
            xf = X % Q
            xc = (Q + 1 - X) % Q
            tot = 0
            for coef, names in terms:
                prod = None
                for nm in names:
                    tb = tables[nm]
                    h = len(tb) >> 1
                    row = (tb[:h] * xc + tb[h:] * xf) % Q
                    prod = row if prod is None else (prod * row) % Q
                tot = (tot + (coef % Q) * int(prod.sum(dtype=np.uint64) % Q)) % Q
            evals.append(tot)
        rec.append(list(evals))
        if (evals[0] + evals[1]) % Q != c:
            return False, rs, None, None
        r = int(rng.integers(0, Q))
        rs.append(r)
        c = lagrange_eval(evals, r)
        rc = (Q + 1 - r) % Q
        for nm in tables:
            tb = tables[nm]
            h = len(tb) >> 1
            tables[nm] = (tb[:h] * rc + tb[h:] * r) % Q
    scalars = {nm: int(tables[nm][0]) for nm in tables}
    return True, rs, c, scalars


def chain_layer_reduce_structured(claim, vj, rho, lt_tab, p, l, wv, wu,
                                  rng, stats, rec=None, q=Q):
    """One backward-sweep chain layer in O(2^l)=O(p) prover work, producing
    the SAME sum-check transcript as the dense product-of-5-tables reduction.
    Returns (ok, rr, final, {'V': v~_j(rho_in)}).  rec, if given, receives
    each round's eval list (so a test can assert round-poly equality with the
    dense reduction; the rec path is default-q test instrumentation).  See the
    block comment above for the factorization."""
    dim = 1 << l
    if rec is not None:
        def sc(cl, tb, tm, dg, rg, qq):
            return _sumcheck_rec(cl, tb, tm, dg, rg, rec)
    else:
        sc = sumcheck

    # ---- phase 1: bind sig'; tables E(y), P(y), G(y) over the y-cube ----
    t0 = time.perf_counter()
    E_y = eq_table(rho, q)                     # eq~(rho, y)
    P_y = lt_tab.copy()                        # [y < p]
    G_y = np.zeros(dim, dtype=_dt(q))          # G(y)=Sum_w R_j(y,w) L(w) v_j(w)
    wbeg = vj[:p].astype(_dt(q))               # L(w)=1 for w<p; v_j zero for w>=p
    wcol = np.arange(p)
    for a in (0, 1):
        for b in (0, 1):
            term = wv[a] * wu[b] % q
            y0 = 2 * wcol + a - b * p
            m = (y0 >= 0) & (y0 < dim)
            np.add.at(G_y, y0[m], _modmul_arr(wbeg[m], term, q))
    G_y %= q
    stats["t_prover"] += time.perf_counter() - t0
    ok1, rs1, c1, sc1 = sc(claim, {"E": E_y, "P": P_y, "G": G_y},
                           [(1, ["E", "P", "G"])], 5, rng, q)
    if not ok1:
        return False, None, None, None
    rho_out = rs1

    # ---- phase 2: bind sig; tables (kappa R*(w)), L(w), v_j(w) over w ----
    t0 = time.perf_counter()
    eqo = eq_table(rho_out, q)                 # eq~(rho_out, .) on the y-cube
    kappa = sc1["E"] * sc1["P"] % q            # E(rho_out) P(rho_out)
    Rstar = np.zeros(dim, dtype=_dt(q))        # R*(w) = R_j(rho_out, w)
    wall = np.arange(dim)
    for a in (0, 1):
        for b in (0, 1):
            term = wv[a] * wu[b] % q
            y0 = 2 * wall + a - b * p
            m = (y0 >= 0) & (y0 < dim)
            Rstar[m] = (Rstar[m] + _modmul_arr(eqo[y0[m]], term, q)) % q
    Rstar = _modmul_arr(Rstar, kappa, q)
    L_w = lt_tab.copy()
    V_w = vj.copy()
    stats["t_prover"] += time.perf_counter() - t0
    ok2, rs2, c2, sc2 = sc(c1, {"R": Rstar, "L": L_w, "V": V_w},
                           [(1, ["R", "L", "V"])], 5, rng, q)
    if not ok2:
        return False, None, None, None
    return True, rs1 + rs2, c2, {"V": sc2["V"]}


def _dense_chain_evals(claim, vj, rho, lt_tab, p, l, wv, wu, rng, rec):
    """Reference: the dense product-of-5-tables backward-sweep reduction,
    recording round polynomials (test only). Mirrors inner_verify_div's loop
    body."""
    dim = 1 << l
    EQb = np.repeat(eq_table(rho), dim)
    LTo = np.repeat(lt_tab, dim)
    Rt = inner_r_table(p, l, wv, wu)
    LTi = np.tile(lt_tab, dim)
    Vb = np.tile(vj, dim)
    return _sumcheck_rec(claim, {"EQ": EQb, "LTo": LTo, "R": Rt,
                                 "LTi": LTi, "V": Vb},
                         [(1, ["EQ", "LTo", "R", "LTi", "V"])], 5, rng, rec)


def inner_verify_div(claim, r_v, r_u, p, n, rng, stats, accept_rem=None,
                     lie=False, structured=False, q=Q):
    """Verify, by GKR over the long-division automaton chain, the value of a
    division-relation wiring MLE at (r_v=dividend bits, r_u=quotient bits).

    accept_rem=None  -> claim = W~(r_v, r_u), W(v,u)=[u=floor(v/p)] (the
                        division wiring: accept summed over ALL remainders).
    accept_rem=m     -> claim = D~(r_v, r_u), D(v,u)=[u=floor(v/p) AND
                        v mod p = m]  (accept only when the final remainder
                        is exactly m).  The affine-image wiring e'=p*e+(p-1)
                        is the m=p-1 case with v<-e', u<-e.

    Only the ACCEPT weighting of sum-check #0 depends on accept_rem; the
    backward chain sweep (which uses LT~ to mask junk states >= p in the
    AGREED EXTENSION, unrelated to acceptance) is identical either way.
    Verifier cost O(n*l), l=ceil(log2 p). Returns ok.

    structured=False (default): the backward sweep builds dense 2^{2l}-size
    tables -> O(n*p^2) PROVER per call (the demo's Õ(x^{3/2}) chain-prover
    bottleneck).  structured=True: chain_layer_reduce_structured runs the
    same reduction in O(n*p) prover with a bit-identical transcript -- the
    verifier checks and rng draws are unchanged, so honest/cheat outcomes
    and all downstream protocols are byte-for-byte the same.  sum-check #0
    is already O(p) and is left as-is."""
    l = max(1, (p - 1).bit_length())
    dim = 1 << l
    t0 = time.perf_counter()
    vs = inner_chain_vectors(r_v, r_u, p, l, n, q)
    if lie:  # prover corrupts the chain consistently from the middle
        vs[n // 2] = vs[n // 2].copy()
        vs[n // 2][0] = (int(vs[n // 2][0]) + 1) % q
        for j in range(n // 2, n):
            wv, wu = bit_weights(r_v[j], q), bit_weights(r_u[j], q)
            nv = np.zeros(dim, dtype=_dt(q))
            for s in range(p):
                for a in (0, 1):
                    for b in (0, 1):
                        y = 2 * s + a - b * p
                        if 0 <= y < p:
                            nv[y] = (int(nv[y]) + int(vs[j][s]) * wv[a] % q
                                     * wu[b]) % q
            vs[j + 1] = nv
    lt_tab = (np.arange(dim) < p).astype(_dt(q))           # junk-state mask
    if accept_rem is None:
        acc_tab = lt_tab.copy()                            # [sigma < p]
    else:
        acc_tab = (np.arange(dim) == accept_rem).astype(_dt(q))
        acc_bits = [(accept_rem >> (l - 1 - i)) & 1 for i in range(l)]
    stats["t_prover"] += time.perf_counter() - t0

    # sum-check #0: claim = sum_sig ACC(sig) * v_n(sig)
    ok, r0, final0, scal0 = sumcheck(claim,
                                     {"ACC": acc_tab,
                                      "V": vs[n].copy()},
                                     [(1, ["ACC", "V"])], 2, rng, q)
    stats["comm"] += l * 3
    if not ok:
        return False
    t0 = time.perf_counter()
    acc_mle = (lt_const_eval(r0, p, l, q) if accept_rem is None
               else eq_point(acc_bits, r0, q))
    if final0 != acc_mle * scal0["V"] % q:
        stats["t_verifier"] += time.perf_counter() - t0
        return False
    stats["t_verifier"] += time.perf_counter() - t0

    rho, claim = r0, scal0["V"]  # claim about v~_n at rho
    for j in range(n - 1, -1, -1):
        wv, wu = bit_weights(r_v[j], q), bit_weights(r_u[j], q)
        # claim = sum_{sig',sig} eq(rho,sig') LT(sig') R_j(sig',sig)
        #                        LT(sig) v~_j(sig)
        if structured:
            ok, rr, final, scal = chain_layer_reduce_structured(
                claim, vs[j], rho, lt_tab, p, l, wv, wu, rng, stats, q=q)
        else:
            t0 = time.perf_counter()
            eqt = eq_table(rho, q)
            EQb = np.repeat(eqt, dim)
            LTo = np.repeat(lt_tab, dim)
            Rt = inner_r_table(p, l, wv, wu, q)
            LTi = np.tile(lt_tab, dim)
            Vb = np.tile(vs[j], dim)
            stats["t_prover"] += time.perf_counter() - t0
            ok, rr, final, scal = sumcheck(claim,
                                           {"EQ": EQb, "LTo": LTo, "R": Rt,
                                            "LTi": LTi, "V": Vb},
                                           [(1, ["EQ", "LTo", "R", "LTi",
                                                 "V"])], 5, rng, q)
        stats["comm"] += 2 * l * 6
        if not ok:
            return False
        rho_out, rho_in = rr[:l], rr[l:]
        t0 = time.perf_counter()
        expect = (eq_point(rho, rho_out, q)
                  * lt_const_eval(rho_out, p, l, q) % q
                  * r_tilde_eval(rho_in, rho_out, p, l, wv, wu, q) % q
                  * lt_const_eval(rho_in, p, l, q) % q
                  * scal["V"]) % q
        okf = final == expect
        stats["t_verifier"] += time.perf_counter() - t0
        if not okf:
            return False
        rho, claim = rho_in, scal["V"]
    # base: v~_0 = MLE of e_0
    t0 = time.perf_counter()
    base = 1
    for ri in rho:
        base = base * ((q + 1 - int(ri)) % q) % q
    okb = claim == base
    stats["t_verifier"] += time.perf_counter() - t0
    return okb


def inner_verify_W(claimW, r_v, r_u, p, n, rng, stats, lie=False,
                   structured=False, q=Q):
    """Division wiring W~(r_v, r_u) = MLE of [u=floor(v/p)]: accept over all
    remainders. Thin wrapper over inner_verify_div(accept_rem=None)."""
    return inner_verify_div(claimW, r_v, r_u, p, n, rng, stats,
                            accept_rem=None, lie=lie, structured=structured,
                            q=q)


# ----------------------------------------------------------------------
# outer protocol with delegated wiring
# ----------------------------------------------------------------------

def verify_layer_delegated(i, p, n, z, claim, S_prev_field, rng, stats,
                           lie_inner=False, structured=False):
    """Same as lucy_dp_verification.verify_layer, but the phase-B final
    wiring value is verified by the inner GKR chain instead of the
    verifier's O(n*p) automaton DP."""
    N = 1 << n
    v = np.arange(N, dtype=np.int64)
    c_i = (i - 1) % Q

    t0 = time.perf_counter()
    EQ = eq_table(z)
    Sf = S_prev_field.copy()
    Bf = (v >= p * p).astype(np.uint64)
    Gf = S_prev_field[v // p]
    stats["t_prover"] += time.perf_counter() - t0

    termsA = [(1, ["EQ", "S"]), (Q - 1, ["EQ", "B", "G"]),
              (c_i, ["EQ", "B"])]
    okA, r_v, finalA, scal = sumcheck(claim, {"EQ": EQ, "S": Sf, "B": Bf,
                                              "G": Gf}, termsA, 3, rng)
    if not okA:
        return False, None, None
    s1, g1 = scal["S"], scal["G"]

    t0 = time.perf_counter()
    eqv = eq_point(z, r_v)
    bv = ge_const_eval(r_v, p * p, n)
    ok = (finalA % Q) == eqv * ((s1 - bv * ((g1 - c_i) % Q)) % Q) % Q
    stats["t_verifier"] += time.perf_counter() - t0
    if not ok:
        return False, None, None

    t0 = time.perf_counter()
    Wt = np.zeros(N, dtype=np.uint64)
    np.add.at(Wt, v // p, eq_table(r_v))
    Wt %= Q
    stats["t_prover"] += time.perf_counter() - t0

    okB, r_u, finalB, scalB = sumcheck(g1, {"W": Wt,
                                            "S": S_prev_field.copy()},
                                       [(1, ["W", "S"])], 2, rng)
    if not okB:
        return False, None, None
    s2 = scalB["S"]

    # prover claims the wiring value (computed by its own automaton DP);
    # verifier checks the algebra, then verifies the claim via the inner
    # chain at O(n log p) instead of evaluating it at O(n p)
    t0 = time.perf_counter()
    wv_claim = w_div_eval(r_v, r_u, p, n)
    if lie_inner:
        wv_claim = (wv_claim + 1) % Q
    stats["t_prover"] += time.perf_counter() - t0

    t0 = time.perf_counter()
    ok = (finalB % Q) == wv_claim * s2 % Q
    stats["t_verifier"] += time.perf_counter() - t0
    if not ok:
        return False, None, None
    if not inner_verify_W(wv_claim, r_v, r_u, p, n, rng, stats,
                          structured=structured):
        return False, None, None

    # line batching s1@r_v, s2@r_u
    t0 = time.perf_counter()
    hs = []
    for t in range(n + 1):
        tf = t % Q
        gpt = [((Q + 1 - tf) * int(a) + tf * int(b)) % Q
               for a, b in zip(r_v, r_u)]
        hs.append(mle_eval(S_prev_field, gpt))
    stats["t_prover"] += time.perf_counter() - t0
    t0 = time.perf_counter()
    if hs[0] != s1 or hs[1] != s2:
        stats["t_verifier"] += time.perf_counter() - t0
        return False, None, None
    tstar = int(rng.integers(0, Q))
    new_claim = lagrange_eval(hs, tstar)
    tf = tstar % Q
    new_z = [((Q + 1 - tf) * int(a) + tf * int(b)) % Q
             for a, b in zip(r_v, r_u)]
    stats["t_verifier"] += time.perf_counter() - t0
    stats["comm"] += (4 * n) + (3 * n) + (n + 1) + 3
    return True, new_z, new_claim


def run_protocol_delegated(n, rng, corrupt_layer=None, lie_inner_at=None,
                           structured=False):
    x = (1 << n) - 1
    primes = [p for p in range(2, int(x ** 0.5) + 1)
              if all(p % d for d in range(2, int(p ** 0.5) + 1))]
    K = len(primes)
    S = prover_tables(n, primes)
    if corrupt_layer is not None:
        i = corrupt_layer
        S[i] = S[i].copy()
        S[i][x] += 1
        v = np.arange(1 << n, dtype=np.int64)
        for j in range(i + 1, K + 1):
            p = primes[j - 1]
            b = v >= p * p
            S[j] = S[j - 1] - b * (S[j - 1][v // p] - (j - 1))
    claimed_pi = int(S[K][x])
    Sfield = [(s % Q).astype(np.uint64) for s in S]
    stats = {"t_prover": 0.0, "t_verifier": 0.0, "comm": 0}
    z = [(x >> (n - 1 - j)) & 1 for j in range(n)]
    claim = claimed_pi % Q
    for i in range(K, 0, -1):
        ok, z, claim = verify_layer_delegated(
            i, primes[i - 1], n, z, claim, Sfield[i - 1], rng, stats,
            lie_inner=(lie_inner_at == i), structured=structured)
        if not ok:
            return dict(accepted=False, claimed=claimed_pi, layers=K,
                        **stats)
    t0 = time.perf_counter()
    ok = s0_closed_form(z) == claim % Q
    stats["t_verifier"] += time.perf_counter() - t0
    return dict(accepted=ok, claimed=claimed_pi, layers=K, **stats)


# ----------------------------------------------------------------------
# tests
# ----------------------------------------------------------------------

def selftest():
    rng = np.random.default_rng(0)
    # 1. affine relation MLE vs brute force, exhaustive over small l and
    #    constants exceeding the l-bit range (the p=2 regression case)
    for l in (1, 2, 3, 4):
        for c in range(-(1 << (l + 1)), 3):
            for x in range(1 << l):
                for y in range(1 << l):
                    bx = [(x >> (l - 1 - i)) & 1 for i in range(l)]
                    by = [(y >> (l - 1 - i)) & 1 for i in range(l)]
                    want = 1 if y == 2 * x + c else 0
                    got = aff_rel_eval(bx, by, c, l)
                    assert got == want, (l, c, x, y, got)
    # 2. R~ table vs r_tilde_eval at random points
    for p in (3, 7, 13):
        l = max(1, (p - 1).bit_length())
        wv, wu = bit_weights(rng.integers(0, Q)), bit_weights(
            rng.integers(0, Q))
        Rt = inner_r_table(p, l, wv, wu)
        for _ in range(10):
            ro = [int(rng.integers(0, Q)) for _ in range(l)]
            ri = [int(rng.integers(0, Q)) for _ in range(l)]
            assert mle_eval(Rt, ro + ri) == r_tilde_eval(ri, ro, p, l,
                                                         wv, wu)
    # 3. inner chain total equals direct automaton MLE; inner protocol
    #    accepts truth and rejects truth+1
    n = 10
    for p in (3, 13, 31):
        l = max(1, (p - 1).bit_length())
        r_v = [int(rng.integers(0, Q)) for _ in range(n)]
        r_u = [int(rng.integers(0, Q)) for _ in range(n)]
        vs = inner_chain_vectors(r_v, r_u, p, l, n)
        chainW = int(sum(int(vs[n][s]) for s in range(p)) % Q)
        trueW = w_div_eval(r_v, r_u, p, n)
        assert chainW == trueW, (p, chainW, trueW)
        st = {"t_prover": 0, "t_verifier": 0, "comm": 0}
        assert inner_verify_W(trueW, r_v, r_u, p, n,
                              np.random.default_rng(1), st)
        assert not inner_verify_W((trueW + 1) % Q, r_v, r_u, p, n,
                                  np.random.default_rng(2), st)
        assert not inner_verify_W(trueW, r_v, r_u, p, n,
                                  np.random.default_rng(3), st, lie=True)
    # 3b. affine-image relation D~(v,u) = [v = p*u + (p-1)] via the
    #     accept_rem = p-1 path, vs brute full-table MLE (small n) and vs the
    #     chain accept weight vs[n][p-1]; rejects truth+1 and a lying chain.
    for (n2, p) in [(6, 3), (6, 5), (7, 7)]:
        l = max(1, (p - 1).bit_length())
        full = np.zeros(1 << (2 * n2), dtype=np.uint64)
        for u in range(1 << n2):
            v = p * u + (p - 1)
            if v < (1 << n2):
                full[(v << n2) | u] = 1          # high n2 bits = v, low = u
        for _ in range(4):
            r_v = [int(rng.integers(0, Q)) for _ in range(n2)]   # dividend v
            r_u = [int(rng.integers(0, Q)) for _ in range(n2)]   # quotient u
            vs = inner_chain_vectors(r_v, r_u, p, l, n2)
            chainD = int(vs[n2][p - 1] % Q)
            bruteD = mle_eval(full, r_v + r_u)
            assert chainD == bruteD, (n2, p, chainD, bruteD)
            st = {"t_prover": 0, "t_verifier": 0, "comm": 0}
            assert inner_verify_div(bruteD, r_v, r_u, p, n2,
                                    np.random.default_rng(11), st,
                                    accept_rem=p - 1)
            assert not inner_verify_div((bruteD + 1) % Q, r_v, r_u, p, n2,
                                        np.random.default_rng(12), st,
                                        accept_rem=p - 1)
            assert not inner_verify_div(bruteD, r_v, r_u, p, n2,
                                        np.random.default_rng(13), st,
                                        accept_rem=p - 1, lie=True)
    # 5. STRUCTURED chain-layer reduction (S496): its per-round polynomials
    #    equal the dense product-of-5-tables ones at EVERY round (exhaustive
    #    small p, random points and chain vectors), and the returned
    #    challenges / final claim / scal['V'] coincide -> bit-identical
    #    transcript.  This is what guarantees structured=True is a drop-in.
    for p in (2, 3, 5, 7, 13, 16, 31):
        l = max(1, (p - 1).bit_length())
        dim = 1 << l
        lt_tab = (np.arange(dim) < p).astype(np.uint64)
        for _ in range(6):
            vj = np.zeros(dim, dtype=np.uint64)
            for s in range(p):                       # support in [0, p)
                vj[s] = int(rng.integers(0, Q))
            rho = [int(rng.integers(0, Q)) for _ in range(l)]
            wv = bit_weights(rng.integers(0, Q))
            wu = bit_weights(rng.integers(0, Q))
            EQb = np.repeat(eq_table(rho), dim)
            Rt = inner_r_table(p, l, wv, wu)
            F = (EQb * np.repeat(lt_tab, dim) % Q * Rt % Q
                 * np.tile(lt_tab, dim) % Q * np.tile(vj, dim) % Q)
            claim_true = int(F.sum(dtype=np.uint64) % Q)
            rd, rs_ = [], []
            okd, rsd, fd, scd = _dense_chain_evals(
                claim_true, vj, rho, lt_tab, p, l, wv, wu,
                np.random.default_rng(77), rd)
            st = {"t_prover": 0.0, "t_verifier": 0.0, "comm": 0}
            oks, rrs, fs, scs = chain_layer_reduce_structured(
                claim_true, vj, rho, lt_tab, p, l, wv, wu,
                np.random.default_rng(77), st, rec=rs_)
            assert okd and oks
            assert rd == rs_, (p, "round polynomials differ")
            assert rsd == rrs and fd == fs and scd["V"] == scs["V"], p
            # a wrong claim makes BOTH reject at round 0, identically
            rd2, rs2 = [], []
            bad = (claim_true + 1) % Q
            okd2, *_ = _dense_chain_evals(bad, vj, rho, lt_tab, p, l, wv, wu,
                                          np.random.default_rng(78), rd2)
            oks2, *_ = chain_layer_reduce_structured(
                bad, vj, rho, lt_tab, p, l, wv, wu,
                np.random.default_rng(78),
                {"t_prover": 0.0, "t_verifier": 0.0, "comm": 0}, rec=rs2)
            assert (not okd2) and (not oks2) and rd2 == rs2
    # 6. structured=True is a drop-in for inner_verify_W / inner_verify_div:
    #    identical accept/reject AND identical comm to the dense path.
    n = 10
    for p in (2, 3, 13, 31):
        r_v = [int(rng.integers(0, Q)) for _ in range(n)]
        r_u = [int(rng.integers(0, Q)) for _ in range(n)]
        trueW = w_div_eval(r_v, r_u, p, n)
        for claimW, lie, want in [(trueW, False, True),
                                  ((trueW + 1) % Q, False, False),
                                  (trueW, True, False)]:
            sd = {"t_prover": 0.0, "t_verifier": 0.0, "comm": 0}
            ss = {"t_prover": 0.0, "t_verifier": 0.0, "comm": 0}
            od = inner_verify_W(claimW, r_v, r_u, p, n,
                                np.random.default_rng(9), sd, lie=lie)
            os_ = inner_verify_W(claimW, r_v, r_u, p, n,
                                 np.random.default_rng(9), ss, lie=lie,
                                 structured=True)
            assert od == os_ == want and sd["comm"] == ss["comm"], (p, claimW)
    # 6b. accept_rem (affine-image) path is also a clean drop-in.
    for (n2, p) in [(6, 3), (7, 7)]:
        l = max(1, (p - 1).bit_length())
        full = np.zeros(1 << (2 * n2), dtype=np.uint64)
        for u in range(1 << n2):
            v = p * u + (p - 1)
            if v < (1 << n2):
                full[(v << n2) | u] = 1
        r_v = [int(rng.integers(0, Q)) for _ in range(n2)]
        r_u = [int(rng.integers(0, Q)) for _ in range(n2)]
        bruteD = mle_eval(full, r_v + r_u)
        for claimD, want in [(bruteD, True), ((bruteD + 1) % Q, False)]:
            od = inner_verify_div(claimD, r_v, r_u, p, n2,
                                  np.random.default_rng(15),
                                  {"t_prover": 0.0, "t_verifier": 0.0,
                                   "comm": 0}, accept_rem=p - 1)
            os_ = inner_verify_div(claimD, r_v, r_u, p, n2,
                                   np.random.default_rng(15),
                                   {"t_prover": 0.0, "t_verifier": 0.0,
                                    "comm": 0}, accept_rem=p - 1,
                                   structured=True)
            assert od == os_ == want
    # 4. tiny end-to-end with delegated wiring (dense and structured prover
    #    give the SAME accepted count == sieve; structured is the O~(x) prover)
    res = run_protocol_delegated(10, np.random.default_rng(4))
    assert res["accepted"] and res["claimed"] == direct_pi((1 << 10) - 1)
    res_s = run_protocol_delegated(10, np.random.default_rng(4),
                                   structured=True)
    assert res_s["accepted"] and res_s["claimed"] == res["claimed"]
    assert res_s["comm"] == res["comm"]
    # structured corrupt-DP and false-wiring liars still rejected
    assert not run_protocol_delegated(
        10, np.random.default_rng(5), corrupt_layer=2, structured=True
    )["accepted"]
    assert not run_protocol_delegated(
        10, np.random.default_rng(6), lie_inner_at=2, structured=True
    )["accepted"]
    print("selftest OK")


def wiring_bench(n=16, trials=5):
    """Verifier wiring-check cost: automaton O(n*p) vs delegated
    O(n*log p), across p."""
    rng = np.random.default_rng(7)
    print(f"{'p':>5} {'automaton_ms':>13} {'delegated_verifier_ms':>22} "
          f"{'delegated_comm(elems)':>22}")
    for p in (3, 13, 31, 61, 127, 251):
        r_v = [int(rng.integers(0, Q)) for _ in range(n)]
        r_u = [int(rng.integers(0, Q)) for _ in range(n)]
        t0 = time.perf_counter()
        for _ in range(trials):
            trueW = w_div_eval(r_v, r_u, p, n)
        t_auto = (time.perf_counter() - t0) / trials
        st = {"t_prover": 0.0, "t_verifier": 0.0, "comm": 0}
        for t in range(trials):
            st_one = {"t_prover": 0.0, "t_verifier": 0.0, "comm": 0}
            assert inner_verify_W(trueW, r_v, r_u, p, n,
                                  np.random.default_rng(100 + t), st_one)
            st["t_verifier"] += st_one["t_verifier"]
            st["comm"] = st_one["comm"]
        print(f"{p:>5} {t_auto*1000:>13.3f} "
              f"{st['t_verifier']/trials*1000:>22.3f} {st['comm']:>22}")


def prover_bench(n=16, trials=2):
    """Chain PROVER wall-clock per delegated wiring evaluation: the dense
    backward sweep (2^{2l}-size tables -> O(n*p^2)) vs the structured
    reduction (O(p)-size tables -> O(n*p)).  Both produce the identical
    transcript (asserted by the selftest); this measures the SCALING in p.
    The win is the exponent: dense ~ p^2, structured ~ p (so the ratio
    grows ~ linearly in p), turning the whole-chain prover Sum_p n p^2 ~
    O~(x^{3/2}) into Sum_p n p ~ O~(x)."""
    rng = np.random.default_rng(7)
    print(f"{'p':>5} {'l':>3} {'dense_ms':>10} {'struct_ms':>10} "
          f"{'ratio':>7} {'dense/p^2':>11} {'struct/p':>10}")
    for p in (7, 13, 31, 61, 127, 251, 509, 1021):
        l = max(1, (p - 1).bit_length())
        r_v = [int(rng.integers(0, Q)) for _ in range(n)]
        r_u = [int(rng.integers(0, Q)) for _ in range(n)]
        trueW = w_div_eval(r_v, r_u, p, n)

        def timed(structured):
            t0 = time.perf_counter()
            for t in range(trials):
                st = {"t_prover": 0.0, "t_verifier": 0.0, "comm": 0}
                ok = inner_verify_W(trueW, r_v, r_u, p, n,
                                    np.random.default_rng(200 + t), st,
                                    structured=structured)
                assert ok
            return (time.perf_counter() - t0) / trials

        td = timed(False)
        ts = timed(True)
        print(f"{p:>5} {l:>3} {td*1000:>10.2f} {ts*1000:>10.2f} "
              f"{td/ts:>7.2f} {td/p**2*1e6:>11.4f} {ts/p*1e6:>10.4f}")


def main(n, cheat_trials, seed):
    x = (1 << n) - 1
    res = run_protocol_delegated(n, np.random.default_rng(seed))
    truth = direct_pi(x)
    print(f"x = 2^{n}-1 = {x}, layers: {res['layers']} (delegated wiring)")
    print(f"honest run: {'ACCEPTED' if res['accepted'] else 'REJECTED'}, "
          f"claimed pi(x) = {res['claimed']}, sieve = {truth}, "
          f"match = {res['claimed'] == truth}")
    assert res["accepted"] and res["claimed"] == truth
    print(f"t_prover = {res['t_prover']:.2f}s  "
          f"t_verifier = {res['t_verifier']*1000:.1f}ms  "
          f"comm = {res['comm']} field elems (~{res['comm']*4/1024:.0f} KB)")
    rej = 0
    for t in range(cheat_trials):
        r = run_protocol_delegated(n, np.random.default_rng(seed + 50 + t),
                                   corrupt_layer=res["layers"] // 2)
        rej += (not r["accepted"])
    print(f"self-consistent DP liar: rejected {rej}/{cheat_trials}")
    rej = 0
    for t in range(cheat_trials):
        r = run_protocol_delegated(n, np.random.default_rng(seed + 90 + t),
                                   lie_inner_at=res["layers"] // 3)
        rej += (not r["accepted"])
    print(f"false wiring-claim liar:  rejected {rej}/{cheat_trials}")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, default=16)
    ap.add_argument("--cheat-trials", type=int, default=10)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--selftest", action="store_true")
    ap.add_argument("--wiring-bench", action="store_true")
    ap.add_argument("--prover-bench", action="store_true")
    args = ap.parse_args()
    if args.selftest:
        selftest()
    elif args.wiring_bench:
        wiring_bench()
    elif args.prover_bench:
        prover_bench()
    else:
        main(args.n, args.cheat_trials, args.seed)
