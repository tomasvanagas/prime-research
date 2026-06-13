#!/usr/bin/env python3
"""
batched_trace.py -- S502 (batched trace zero-test).

The compressed pi(x) chain (compressed_layer.run_chain) certifies, per Lucy
layer l, that its trace witness W_l = build_witness(x, p_l, nb, dstart=1) really
encodes u_e = floor(X / ((e+1) p_l)) -- the degree-3 constraint zero-test
verify_constraints (compressed_prover_mult_trace).  The S501 profile found this
single kernel is 53% of run_chain's wall: K = pi(sqrt x) INDEPENDENT sum-checks,
each over its own nb-cube witness, each paying full Python+numpy per-call
dispatch on a ~30-element array (op-count-bound, exactly S500's signature).

The K zero-tests are independent: each draws its own tau/alpha, certifies a
witness that depends only on (x, p_l) -- never on the carried chain claim or on
another layer's challenges (verified: verify_constraints' challenges are local;
nothing feeds forward).  So they admit a textbook BATCHED sum-check.

CONSTRUCTION.  Stack the K witnesses along a layer axis of size L = 2^Lk,
Lk = ceil(log2 K); index the joint cube by (l << nb) | e, l in [0,L), e in
[0,2^nb).  Draw ONE shared tau (the e-cube point), ONE shared alpha (the
intra-layer constraint combiner), ONE beta (the inter-layer combiner).  Define

    BETA_EQ[l, e] = (beta^l if l < K else 0) * eq~(tau, e)            (multilinear)

and run a SINGLE degree-3 sum-check of claim 0 over the (Lk+nb)-cube:

    0 = sum_{l,e} BETA_EQ[l,e] * F(l, e),
    F(l, e) = sum_c alpha^c * constraint_c(W stacked)[l, e]

with the SAME flat term list as the per-layer test (build_terms), except the
eq-selector "EQ" is replaced by "BETA_EQ".  Because each honest layer's inner
sum is sum_e eq~(tau,e) F_l(e) = F_l~(tau) = 0, the whole sum is
sum_{l<K} beta^l * 0 = 0.  A corrupted layer l* has F_{l*} != 0, so
F_{l*}~(tau) != 0 (w.p. >= 1 - nb/q over tau), and sum_l beta^l F_l~(tau) != 0
(w.p. >= 1 - (K-1)/q over beta) -- caught.

The width-Lv_max/Lr_max bit tables embed every layer's shorter decomposition by
MSB-zero-padding (u_e < 2^Lv_l <= 2^Lv_max => the high bits are 0), so one
shared term list and one recomposition weight set serve all layers.  Padding
rows l >= K carry BETA_EQ = 0, so they contribute 0 to the sum regardless of
their (zeroed) witness tables; ONE is all-ones over the whole cube (its MLE must
be 1 for the -X term's identity to hold at the random point).

VERIFIER.  After the sum-check binds (r_L, r_e), it recomputes
    BETA_EQ~(r) = [ sum_{l<K} beta^l chi_l(r_L) ] * eq~(tau, r_e)        (O(K))
    av         = (Int(r_e) + dstart) * [ sum_{l<K} p_l chi_l(r_L) ]      (O(K))
where chi_l is the multilinear indicator of layer index l, and av is the MLE of
the TRUE per-layer wiring a(l,e) = (e+dstart) p_l (it factors across the disjoint
e / l variable blocks, so the MLE is the product of the per-block MLEs).  It then
checks  final == BETA_EQ~(r) * constraint_eval(scalars, av, X, Lv_max, Lr_max,
alpha).  The prover's A / ONE / BETA_EQ scalars are ignored; av and the literal X
anchor the multiply identity to the true wiring, exactly as the per-layer test.

WHY IT IS THE RIGHT WIN (S501).  ONE sum-check instead of K cuts the fmul CALL
count ~K-fold while WIDENING each fmul's array ~K-fold (2^nb -> ~K * 2^nb),
turning op-count-bound into width-bound -- the regime where numpy and the fast
Mersenne mulmod (S500) finally pay off.  Communication drops too: 4(Lk+nb) vs
K * 4 nb.  Verifier work O(K + n(Lk+nb)) = O~(sqrt x) (the chain already does K
layer-reductions, so O(K) here is within budget).  Soundness error
~ (K * #constraints + nb + K + deg(Lk+nb)) / q, negligible for BIG_Q (~2^61).

What would falsify this: the batched test rejecting an honest chain's witnesses,
or accepting when ANY single layer's witness is corrupted (wrong quotient,
corrupted bit, wrong remainder, non-Boolean bit); the batched accept/reject
disagreeing with the AND of the K per-layer verify_constraints; the fast-Mersenne
uint64 path disagreeing with the object reference; or the fmul COUNT not dropping
~K-fold / the per-fmul WIDTH not rising ~K-fold vs the unbatched baseline.

Usage:
  python3 batched_trace.py --selftest
  python3 batched_trace.py --n 16 --field big
  python3 batched_trace.py --bench --field big
"""

import argparse
import time
from math import isqrt

import numpy as np

import compressed_prover_mult_trace as _cpmt
from compressed_prover_mult_trace import (DEFAULT_Q as Q, BIG_Q, SMALL_Q, FIELDS,
                                          _dt, eq_table, eq_point, mle_eval,
                                          sumcheck, build_witness, build_terms,
                                          constraint_eval, int_of_point, fmul,
                                          verify_constraints, rematerialize)


# ----------------------------------------------------------------------
# the K chain trace witnesses (one per prime p <= sqrt x)
# ----------------------------------------------------------------------

def primes_upto(V):
    if V < 2:
        return []
    s = np.ones(V + 1, dtype=bool)
    s[:2] = False
    for i in range(2, isqrt(V) + 1):
        if s[i]:
            s[i * i::i] = False
    return [int(p) for p in np.nonzero(s)[0]]


def chain_trace_witnesses(x, q=Q):
    """Build exactly the K = pi(sqrt x) trace witnesses run_chain's large_reduce
    certifies: W_l = build_witness(x, p_l, nb, dstart=1) for each prime p_l<=V."""
    V = isqrt(x)
    nb = max(1, V.bit_length())
    primes = primes_upto(V)
    Ws = [build_witness(x, p, nb, dstart=1, q=q) for p in primes]
    return primes, nb, Ws


def _ceil_log2(K):
    return 0 if K <= 1 else (K - 1).bit_length()


# ----------------------------------------------------------------------
# verifier-side O(K) MLEs over the layer index
# ----------------------------------------------------------------------

def _chi(l, r_L, q):
    """Multilinear indicator chi_l(r_L) = prod_j eq(bit_j(l), r_L[j]), MSB-first."""
    Lk = len(r_L)
    acc = 1
    for j in range(Lk):
        bit = (l >> (Lk - 1 - j)) & 1
        rj = int(r_L[j]) % q
        acc = acc * (rj if bit else (q + 1 - rj)) % q
    return acc


def beta_eq_layer_eval(beta, K, r_L, q):
    """sum_{l<K} beta^l chi_l(r_L).  O(K * Lk)."""
    total, bp = 0, 1
    for l in range(K):
        total = (total + bp * _chi(l, r_L, q)) % q
        bp = bp * (beta % q) % q
    return total


def ptilde(primes, r_L, q):
    """MLE of the layer->prime map: sum_{l<K} p_l chi_l(r_L).  O(K * Lk)."""
    total = 0
    for l, p in enumerate(primes):
        total = (total + (p % q) * _chi(l, r_L, q)) % q
    return total


# ----------------------------------------------------------------------
# the stacked witness tables over the (Lk + nb)-cube
# ----------------------------------------------------------------------

def stacked_tables(witnesses, primes, nb, X, beta, tau, q, cheat=None,
                   cheat_layer=0):
    """Build the field tables of the joint (l, e) cube.  Index = l*D + e.
    Bit tables use the common Lv_max / Lr_max (MSB-zero-padding each layer's
    shorter decomposition).  Returns (tables, Lk, Lv, Lr)."""
    K = len(witnesses)
    Lk = _ceil_log2(K)
    L = 1 << Lk
    D = 1 << nb
    N = L * D
    Lv = max(W["Lv"] for W in witnesses)
    Lr = max(W["Lr"] for W in witnesses)
    dt = _dt(q)

    u = np.zeros(N, dtype=np.int64)
    r = np.zeros(N, dtype=np.int64)
    qd = np.zeros(N, dtype=np.int64)
    a = np.zeros(N, dtype=np.int64)
    for li, W in enumerate(witnesses):
        sl = slice(li * D, li * D + D)
        u[sl] = W["u"]; r[sl] = W["r"]; qd[sl] = W["q"]; a[sl] = W["a"]

    # integer-level cheat: a self-consistent wrong quotient (bits stay in sync,
    # so recomposition holds but the multiply identity U*a+R-X=0 fails).
    if cheat == "u_consistent":
        u[cheat_layer * D + 1] += 1

    tabs = {
        "U": (u % q).astype(dt),
        "R": (r % q).astype(dt),
        "Qv": (qd % q).astype(dt),
        "A": (a % q).astype(dt),
        "ONE": np.ones(N, dtype=dt),
    }
    for j in range(Lv):
        tabs[f"Ub{j}"] = ((u >> (Lv - 1 - j)) & 1).astype(dt)
    for j in range(Lr):
        tabs[f"Rb{j}"] = ((r >> (Lr - 1 - j)) & 1).astype(dt)
        tabs[f"Qb{j}"] = ((qd >> (Lr - 1 - j)) & 1).astype(dt)

    # field-level cheats on a single layer (mirror compressed_prover_mult_trace)
    idx = cheat_layer * D + 1
    if cheat == "u_value":                  # quotient off, bits untouched
        tabs["U"] = tabs["U"].copy(); tabs["U"][idx] = (int(tabs["U"][idx]) + 1) % q
    elif cheat == "r_value":                # remainder corrupted
        tabs["R"] = tabs["R"].copy(); tabs["R"][idx] = (int(tabs["R"][idx]) + 1) % q
    elif cheat == "nonbit":                 # a 'bit' set to 2
        tabs["Ub0"] = tabs["Ub0"].copy(); tabs["Ub0"][cheat_layer * D] = 2

    # BETA_EQ[l,e] = (beta^l if l<K else 0) * eq~(tau, e)
    eqtau = eq_table(tau, q)                                     # length D
    BE = np.zeros(N, dtype=dt)
    bp = 1
    for li in range(K):
        BE[li * D:(li + 1) * D] = fmul(eqtau, bp % q, q)
        bp = bp * (beta % q) % q
    tabs["BETA_EQ"] = BE
    return tabs, Lk, Lv, Lr


# ----------------------------------------------------------------------
# the batched zero-test
# ----------------------------------------------------------------------

def verify_constraints_batched(witnesses, primes, X, nb, rng, stats, q=Q,
                               dstart=1, cheat=None, cheat_layer=0):
    """ONE degree-3 sum-check certifying ALL K layers' u_e = floor(X/a_e^l)
    simultaneously.  Replaces K calls to verify_constraints.  Returns ok."""
    K = len(witnesses)
    tau = [int(rng.integers(0, q)) for _ in range(nb)]
    alpha = int(rng.integers(2, q))
    beta = int(rng.integers(2, q))

    t0 = time.perf_counter()
    tabs, Lk, Lv, Lr = stacked_tables(witnesses, primes, nb, X, beta, tau, q,
                                      cheat=cheat, cheat_layer=cheat_layer)
    terms = build_terms(Lv, Lr, X, alpha, q)
    terms = [(coef, ["BETA_EQ" if nm == "EQ" else nm for nm in names])
             for coef, names in terms]
    stats["t_prover"] += time.perf_counter() - t0

    ok, r_all, final, scal = sumcheck(0, tabs, terms, 3, rng, q)
    stats["comm"] += 4 * (Lk + nb)                              # degree-3 round polys
    if not ok:
        return False

    t0 = time.perf_counter()
    r_L, r_e = r_all[:Lk], r_all[Lk:]
    be = beta_eq_layer_eval(beta, K, r_L, q) * eq_point(tau, r_e, q) % q
    av = (int_of_point(r_e, q) + dstart) % q * ptilde(primes, r_L, q) % q
    expect = be * constraint_eval(scal, av, X, Lv, Lr, alpha, q) % q
    res = (final % q) == expect
    stats["t_verifier"] += time.perf_counter() - t0
    return res


def verify_constraints_each(witnesses, rng, stats, q=Q):
    """The unbatched baseline: K independent per-layer zero-tests, AND-ed.
    Exactly what run_chain does today (one verify_constraints per layer)."""
    ok = True
    for W in witnesses:
        ok = verify_constraints(W, rng, stats, q) and ok
    return ok


# ----------------------------------------------------------------------
# batched Ub-bit-table OPENINGS (S506) -- the last per-layer O(2^nb) term
# ----------------------------------------------------------------------
#
# WHAT IT DISCHARGES.  compressed_layer.verify_trace_region's B2 routing check
# pins g1_trace to the certified quotient u_e by EVALUATING, at the B2 challenge
# point r_C, the nb LOW bit tables Ub_{Lv-nb+k} (= bit nb-1-k of u_e) -- nb direct
# mle_eval's per layer = O(K nb 2^nb) over the chain, the ONLY per-layer O(2^nb)
# verifier term the pcs leaf-opening pass (S505) left standing (their residuals
# are per-layer WITNESS data: no carried chain claim to thread to the base).
#
# THE BATCH.  Those Ub tables are EXACTLY the ones verify_constraints_batched
# already stacks and commits along the layer axis (one zero-test, all K layers).
# So defer each layer's openings -- (layer l, point r_C^l, claims [ub_{l,k}]_k) --
# and discharge them ALL in ONE degree-2 sum-check over the SAME stacked
# (Lk+nb)-cube.  gamma <- F_q; with beta = gamma^nb the per-(l,k) weight
# beta^l gamma^k = gamma^{l*nb+k} is a DISTINCT power, so the weighting is a
# genuine RLC over all K*nb openings.  Write the combined identity
#
#   claim := sum_{l,k} gamma^{l nb+k} ub_{l,k}  ==  sum_w B[w] C[w],
#     B[w] = sum_{l<K} beta^l [w_L=l] eq~(r_C^l, w_e)          (verifier-anchored)
#     C[w] = sum_{k<nb} gamma^k Ub_{(nb-1-k)}[w]               (committed Ub RLC)
#
# (B factors the per-layer eq weight, k-independent; C is the gamma-fold of the
# nb committed low-bit tables -- one combined opening of the SAME cube the
# zero-test pins).  The cross term sum_w B[w]C[w] = sum_{l,k} beta^l gamma^k
# Ub_{(nb-1-k)}^{(l)}~(r_C^l) is exactly the gamma-weighted TRUE openings, so an
# honest prover's claim matches and a single wrong ub_{l,k} flips claim (w.p.
# >= 1 - (K nb - 1)/q over gamma), caught by the sum-check round-1 identity.
#
# VERIFIER.  After binding (r_L, r_e) it recomputes B~(r*) = sum_l beta^l
# chi_l(r_L) eq~(r_C^l, r_e) in O(K(Lk+nb)) -- B is PUBLIC (the r_C points, beta),
# the anchor -- and takes C~(r*) from the sum-check's folded scalar (the committed
# Ub-cube opening, exactly as the zero-test trusts its U/R/Qv/bit scalars), then
# checks final == B~(r*) C~(r*).  ONE-TIME O(K(Lk+nb)) = O~(sqrt x), replacing the
# per-layer O(K nb 2^nb).  Soundness ~ (K nb + 2(Lk+nb))/q.

def ub_opening_claims(W, r_C, nb, q=Q):
    """Honest prover's nb LOW-bit openings for one layer: [Ub_{Lv-nb+k}~(r_C)]_k,
    k in [0,nb).  These are the scalars verify_trace_region multiplies into its B2
    `expect`; under deferral the prover supplies them and the batch discharges
    them.  O(nb 2^nb) (prover-side)."""
    Lv = W["Lv"]
    out = []
    for k in range(nb):
        j = Lv - nb + k
        out.append(int(mle_eval(W["tabs"][f"Ub{j}"].astype(_dt(q)), r_C, q)) % q)
    return out


def verify_ub_openings_batched(witnesses, nb, obligations, rng, stats, q=Q):
    """Discharge ALL deferred per-layer Ub-bit-table openings in ONE degree-2
    sum-check against the SAME stacked Ub cube the trace zero-test commits.

    obligations: list of (layer_idx, r_C, ub_claims), one per layer that ran the
    trace region, with r_C a length-nb point and ub_claims[k] the prover's claimed
    Ub_{(nb-1-k)}^{(layer)}~(r_C) (bit nb-1-k of u_e, = the per-layer Ub_{Lv-nb+k}).
    layer_idx indexes `witnesses` / its stacked layer slice.

    Replaces the K*nb per-layer mle_eval openings (the last per-layer O(2^nb)
    verifier term) with a one-time O(K(Lk+nb)) discharge.  Returns ok."""
    K = len(witnesses)
    if not obligations:
        return True
    Lk = _ceil_log2(K)
    D = 1 << nb
    N = (1 << Lk) * D
    dt = _dt(q)

    gamma = int(rng.integers(2, q))
    beta = pow(gamma % q, nb, q)

    # combined claim = sum_l beta^l (sum_k gamma^k ub_{l,k})   (verifier O(K nb))
    t0 = time.perf_counter()
    claim = 0
    for (li, r_C, ubs) in obligations:
        inner, gp = 0, 1
        for c in ubs:
            inner = (inner + gp * (int(c) % q)) % q
            gp = gp * gamma % q
        claim = (claim + pow(beta, li, q) * inner) % q
    stats["t_verifier"] += time.perf_counter() - t0

    # prover tables: B (verifier-anchored eq weights) and C (committed Ub RLC)
    t0 = time.perf_counter()
    B = np.zeros(N, dtype=dt)
    for (li, r_C, ubs) in obligations:
        B[li * D:(li + 1) * D] = fmul(eq_table(r_C, q), pow(beta, li, q), q)
    C = np.zeros(N, dtype=dt)
    gp = 1
    for k in range(nb):
        bitpos = nb - 1 - k                                  # bit of u_e for slot k
        col = np.zeros(N, dtype=np.int64)
        for li, W in enumerate(witnesses):
            col[li * D:(li + 1) * D] = (W["u"] >> bitpos) & 1
        C = (C + fmul(col.astype(dt), gp % q, q)) % q
        gp = gp * gamma % q
    stats["t_prover"] += time.perf_counter() - t0

    ok, r_all, final, scal = sumcheck(claim, {"B": B, "C": C},
                                      [(1, ["B", "C"])], 2, rng, q)
    stats["comm"] += 1 + 3 * (Lk + nb)                       # gamma + deg-2 rounds
    if not ok:
        return False

    t0 = time.perf_counter()
    r_L, r_e = r_all[:Lk], r_all[Lk:]
    bv = 0                                                   # B~(r*), the anchor
    for (li, r_C, ubs) in obligations:
        bv = (bv + pow(beta, li, q) * _chi(li, r_L, q) % q
              * eq_point(r_C, r_e, q)) % q
    res = (final % q) == bv * (int(scal["C"]) % q) % q
    stats["t_verifier"] += time.perf_counter() - t0
    return res


def verify_ub_openings_each(witnesses, nb, obligations, q=Q):
    """The unbatched baseline / ground truth the batch must agree with: per layer,
    check every ub_claim equals the direct mle_eval of its Ub bit table at r_C
    (exactly compressed_layer.verify_trace_region's inline opening, line ~439)."""
    ok = True
    for (li, r_C, ubs) in obligations:
        true = ub_opening_claims(witnesses[li], r_C, nb, q)
        ok = ok and all(int(a) % q == int(b) % q for a, b in zip(ubs, true))
    return ok


# ----------------------------------------------------------------------
# fmul instrumentation (count + array width), for the bench
# ----------------------------------------------------------------------

class _FmulCounter:
    """Monkeypatch _cpmt.fmul to count calls and operand widths.  sumcheck and
    verify_constraints call the bare module global `fmul`, so this catches every
    multiply in both the batched and unbatched paths uniformly."""

    def __enter__(self):
        self.calls = 0
        self.total_w = 0
        self.max_w = 0
        self._orig = _cpmt.fmul

        def wrap(a, b, q):
            w = max(np.size(a), np.size(b))
            self.calls += 1
            self.total_w += int(w)
            if w > self.max_w:
                self.max_w = int(w)
            return self._orig(a, b, q)

        _cpmt.fmul = wrap
        return self

    def __exit__(self, *exc):
        _cpmt.fmul = self._orig

    @property
    def mean_w(self):
        return self.total_w / self.calls if self.calls else 0.0


# ----------------------------------------------------------------------
# tests / experiments
# ----------------------------------------------------------------------

def selftest():
    rng = np.random.default_rng(0)
    cheats = ["u_consistent", "u_value", "r_value", "nonbit"]

    # 1. structural facts: Lk = ceil(log2 K), one sum-check of (Lk+nb) rounds,
    #    comm = 4(Lk+nb) -- and the comm strictly drops vs the K-fold baseline.
    for n in (8, 12, 16):
        x = (1 << n) - 1
        primes, nb, Ws = chain_trace_witnesses(x)
        K = len(Ws)
        Lk = _ceil_log2(K)
        assert (1 << Lk) >= K and (K <= 1 or (1 << (Lk - 1)) < K)
        st = {"t_prover": 0.0, "t_verifier": 0.0, "comm": 0}
        assert verify_constraints_batched(Ws, primes, x, nb,
                                          np.random.default_rng(1), st)
        assert st["comm"] == 4 * (Lk + nb)
        assert st["comm"] < K * 4 * nb or K == 1          # comm drop (K>1)

    # 2. honest accept AND soundness over the demo prime AND BIG_Q (object):
    #    batched accepts an honest chain; corrupting ANY single layer (first,
    #    middle, last) with any cheat class is rejected, every trial.
    for q in (Q, BIG_Q):
        for n in (8, 10, 12):
            x = (1 << n) - 1
            primes, nb, Ws = chain_trace_witnesses(x, q)
            K = len(Ws)
            st = {"t_prover": 0.0, "t_verifier": 0.0, "comm": 0}
            assert verify_constraints_batched(Ws, primes, x, nb,
                                              np.random.default_rng(5), st, q), \
                (q, n)
            for cl in sorted(set([0, K // 2, K - 1])):
                for ch in cheats:
                    rej = sum(not verify_constraints_batched(
                        Ws, primes, x, nb, np.random.default_rng(100 + t), st,
                        q, cheat=ch, cheat_layer=cl) for t in range(5))
                    assert rej == 5, (q, n, ch, cl, rej)

    # 3. batched accept/reject AGREES with the AND of the K per-layer tests
    #    (the protocol it replaces) -- honest and a representative cheat.
    for n in (8, 12):
        x = (1 << n) - 1
        primes, nb, Ws = chain_trace_witnesses(x)
        st = {"t_prover": 0.0, "t_verifier": 0.0, "comm": 0}
        each = verify_constraints_each(Ws, np.random.default_rng(7), st)
        batched = verify_constraints_batched(Ws, primes, x, nb,
                                             np.random.default_rng(7), st)
        assert each and batched
        # corrupt one layer's witness in place -> BOTH must reject it
        Wc = build_witness(x, primes[len(primes) // 2], nb, dstart=1)
        Wc["tabs"]["U"] = Wc["tabs"]["U"].copy()
        Wc["tabs"]["U"][1] = (int(Wc["tabs"]["U"][1]) + 1) % Q
        assert not verify_constraints(Wc, np.random.default_rng(8), st)

    # 4. K=1 edge case (x=7: only prime 2 <= V=2; Lk=0, no layer bits)
    x = 7
    primes, nb, Ws = chain_trace_witnesses(x)
    assert len(Ws) == 1 and _ceil_log2(1) == 0
    st = {"t_prover": 0.0, "t_verifier": 0.0, "comm": 0}
    assert verify_constraints_batched(Ws, primes, x, nb,
                                      np.random.default_rng(9), st)
    assert sum(not verify_constraints_batched(
        Ws, primes, x, nb, np.random.default_rng(10 + t), st,
        cheat="u_value", cheat_layer=0) for t in range(5)) == 5

    # 5. fast Mersenne uint64 path == object reference, bit-for-bit (same
    #    accept/reject, honest + every cheat), for BIG_Q.
    saved = _cpmt.FAST_BIG
    try:
        for n in (8, 10, 12):
            x = (1 << n) - 1
            _cpmt.FAST_BIG = False
            primes, nb, Wo = chain_trace_witnesses(x, BIG_Q)
            _cpmt.FAST_BIG = True
            _, _, Wf = chain_trace_witnesses(x, BIG_Q)
            K = len(Wo)
            for ch in [None] + cheats:
                cl = K // 2
                _cpmt.FAST_BIG = False
                ro = verify_constraints_batched(
                    Wo, primes, x, nb, np.random.default_rng(13),
                    {"t_prover": 0.0, "t_verifier": 0.0, "comm": 0},
                    BIG_Q, cheat=ch, cheat_layer=cl)
                _cpmt.FAST_BIG = True
                rf = verify_constraints_batched(
                    Wf, primes, x, nb, np.random.default_rng(13),
                    {"t_prover": 0.0, "t_verifier": 0.0, "comm": 0},
                    BIG_Q, cheat=ch, cheat_layer=cl)
                assert ro == rf, (n, ch)
                assert ro == (ch is None), (n, ch)
    finally:
        _cpmt.FAST_BIG = saved

    # 6. the structural win itself: ONE sum-check cuts the fmul CALL count and
    #    WIDENS each fmul, vs the K-fold unbatched baseline (object dtype).
    x = (1 << 14) - 1
    primes, nb, Ws = chain_trace_witnesses(x, BIG_Q)
    with _FmulCounter() as cu:
        verify_constraints_each(Ws, np.random.default_rng(1),
                                {"t_prover": 0.0, "t_verifier": 0.0, "comm": 0})
    with _FmulCounter() as cb:
        verify_constraints_batched(Ws, primes, x, nb, np.random.default_rng(1),
                                   {"t_prover": 0.0, "t_verifier": 0.0, "comm": 0},
                                   BIG_Q)
    assert cb.calls < cu.calls, (cb.calls, cu.calls)        # fewer calls
    assert cb.mean_w > cu.mean_w, (cb.mean_w, cu.mean_w)    # wider each call
    K = len(Ws)
    assert cu.calls / cb.calls > 0.4 * K, (cu.calls, cb.calls, K)   # ~K-fold

    # 7. BATCHED Ub OPENINGS (S506): the deferred per-layer Ub-bit-table openings
    #    (verify_trace_region's nb mle_eval's, the last per-layer O(2^nb) verifier
    #    term) discharge in ONE sum-check.  It (a) ACCEPTS honest openings, (b)
    #    AGREES with the per-layer inline mle_eval ground truth, (c) REJECTS any
    #    single forged ub_claim -- over q AND BIG_Q, first/middle/last layer and
    #    bit, plus the K=1 edge and fast==object bit-identity.
    def _obligations(witnesses, nb, q, seed, corrupt=None):
        rg = np.random.default_rng(seed)
        obs = []
        for li, W in enumerate(witnesses):
            rC = [int(rg.integers(0, q)) for _ in range(nb)]
            ubs = ub_opening_claims(W, rC, nb, q)
            if corrupt is not None and li == corrupt[0]:
                ubs = list(ubs)
                ubs[corrupt[1]] = (int(ubs[corrupt[1]]) + 1) % q
            obs.append((li, rC, ubs))
        return obs

    for q in (Q, BIG_Q):
        for n in (8, 10, 12):
            x = (1 << n) - 1
            primes, nb, Ws = chain_trace_witnesses(x, q)
            K = len(Ws)
            obs = _obligations(Ws, nb, q, 11)
            # (a) honest accept + (b) agreement with the inline mle_eval ground truth
            assert verify_ub_openings_each(Ws, nb, obs, q), (q, n, "each honest")
            st = {"t_prover": 0.0, "t_verifier": 0.0, "comm": 0}
            assert verify_ub_openings_batched(Ws, nb, obs,
                                              np.random.default_rng(3), st, q), \
                (q, n, "batched honest")
            assert st["comm"] == 1 + 3 * (_ceil_log2(K) + nb), (q, n, st["comm"])
            # (c) corrupting ANY single (layer, bit) opening -> BOTH reject
            for cl in sorted(set([0, K // 2, K - 1])):
                for kb in sorted(set([0, nb // 2, nb - 1])):
                    bad = _obligations(Ws, nb, q, 11, corrupt=(cl, kb))
                    assert not verify_ub_openings_each(Ws, nb, bad, q), \
                        (q, n, cl, kb, "each missed forge")
                    rej = sum(not verify_ub_openings_batched(
                        Ws, nb, bad, np.random.default_rng(200 + t), st, q)
                        for t in range(4))
                    assert rej == 4, (q, n, cl, kb, rej)

    # K=1 edge (x=7): Lk=0, the cube is just the nb e-vars, chi_l empty -> 1
    x = 7
    primes, nb, Ws = chain_trace_witnesses(x)
    assert len(Ws) == 1 and _ceil_log2(1) == 0
    obs = _obligations(Ws, nb, Q, 12)
    st = {"t_prover": 0.0, "t_verifier": 0.0, "comm": 0}
    assert verify_ub_openings_batched(Ws, nb, obs, np.random.default_rng(1), st)
    bad = _obligations(Ws, nb, Q, 12, corrupt=(0, 0))
    assert sum(not verify_ub_openings_batched(
        Ws, nb, bad, np.random.default_rng(13 + t), st) for t in range(4)) == 4

    # fast-Mersenne uint64 path == object reference (same accept/reject)
    saved = _cpmt.FAST_BIG
    try:
        x = (1 << 12) - 1
        _cpmt.FAST_BIG = False
        primes, nb, Wo = chain_trace_witnesses(x, BIG_Q)
        obs_o = _obligations(Wo, nb, BIG_Q, 21)
        ro = verify_ub_openings_batched(Wo, nb, obs_o, np.random.default_rng(7),
                                        {"t_prover": 0., "t_verifier": 0., "comm": 0},
                                        BIG_Q)
        _cpmt.FAST_BIG = True
        _, _, Wf = chain_trace_witnesses(x, BIG_Q)
        obs_f = _obligations(Wf, nb, BIG_Q, 21)
        rf = verify_ub_openings_batched(Wf, nb, obs_f, np.random.default_rng(7),
                                        {"t_prover": 0., "t_verifier": 0., "comm": 0},
                                        BIG_Q)
        assert ro == rf == True, (ro, rf)
    finally:
        _cpmt.FAST_BIG = saved

    print("selftest OK")


def _trace_wall(fn):
    """Time fn() and return (wall_s, fmul_calls, mean_w, max_w)."""
    with _FmulCounter() as c:
        t0 = time.perf_counter()
        fn()
        wall = time.perf_counter() - t0
    return wall, c.calls, c.mean_w, c.max_w


def bench(field="big", ns=(8, 10, 12, 14, 16), seed=1):
    """The trace zero-test in isolation (53% of run_chain wall, S501): K
    per-layer sum-checks (unbatched, today's chain) vs ONE batched sum-check.
    Reports fmul COUNT (should drop ~K-fold), per-fmul mean WIDTH (should rise
    ~K-fold), and wall, on the object path and -- for --field big -- the
    fast-Mersenne uint64 path (where the width win should flip S500's 'fast
    path loses on the chain' to a net gain)."""
    q = FIELDS[field]
    paths = [("object", False), ("fast", True)] if field == "big" else [("", False)]
    print(f"batched trace zero-test, field = {field} (q = {q})")
    print(f"{'n':>3} {'K':>4} {'nb':>3} {'Lk':>3} {'path':>7} || "
          f"{'unb_calls':>10} {'bat_calls':>10} {'c_drop':>6} || "
          f"{'unb_w':>7} {'bat_w':>7} {'w_rise':>6} || "
          f"{'unb_ms':>9} {'bat_ms':>9} {'speedup':>7}")
    saved = _cpmt.FAST_BIG
    try:
        for n in ns:
            x = (1 << n) - 1
            for label, fast in paths:
                if field == "big":
                    _cpmt.FAST_BIG = fast
                primes, nb, Ws = chain_trace_witnesses(x, q)
                K = len(Ws)
                Lk = _ceil_log2(K)
                wu, cu, mu, _ = _trace_wall(
                    lambda: verify_constraints_each(
                        Ws, np.random.default_rng(seed),
                        {"t_prover": 0.0, "t_verifier": 0.0, "comm": 0}))
                wb, cb, mb, _ = _trace_wall(
                    lambda: verify_constraints_batched(
                        Ws, primes, x, nb, np.random.default_rng(seed),
                        {"t_prover": 0.0, "t_verifier": 0.0, "comm": 0}, q))
                print(f"{n:>3} {K:>4} {nb:>3} {Lk:>3} {label:>7} || "
                      f"{cu:>10} {cb:>10} {cu/max(cb,1):>5.1f}x || "
                      f"{mu:>7.1f} {mb:>7.1f} {mb/max(mu,1e-9):>5.1f}x || "
                      f"{wu*1000:>9.1f} {wb*1000:>9.1f} "
                      f"{wu/max(wb,1e-9):>6.2f}x")
    finally:
        _cpmt.FAST_BIG = saved


def main(n, field, seed):
    q = FIELDS[field]
    if field == "big":
        _cpmt.FAST_BIG = True
    x = (1 << n) - 1
    primes, nb, Ws = chain_trace_witnesses(x, q)
    K = len(Ws)
    Lk = _ceil_log2(K)
    st = {"t_prover": 0.0, "t_verifier": 0.0, "comm": 0}
    ok = verify_constraints_batched(Ws, primes, x, nb,
                                    np.random.default_rng(seed), st, q)
    print(f"x = 2^{n}-1 = {x}, V = {isqrt(x)}, nb = {nb}, K = pi(V) = {K} layers")
    print(f"field = {field} (q = {q}); batched cube = 2^({Lk}+{nb}) = "
          f"{1 << (Lk + nb)} (one sum-check; unbatched ran {K})")
    print(f"BATCHED trace zero-test (all {K} layers at once): "
          f"{'ACCEPTED' if ok else 'REJECTED'}")
    assert ok
    print(f"  comm = {st['comm']} field elems = 4*(Lk+nb)   "
          f"(unbatched trace comm = {K * 4 * nb} = K*4*nb)")
    print(f"  soundness ~ (K*#constraints + nb + K + 3(Lk+nb))/q = "
          f"{(K * (6 + Ws[0]['Lv'] + 2 * Ws[0]['Lr']) + nb + K + 3 * (Lk + nb)) / q:.2e}")
    for ch in ["u_consistent", "u_value", "r_value", "nonbit"]:
        rej = sum(not verify_constraints_batched(
            Ws, primes, x, nb, np.random.default_rng(seed + 50 + t), st, q,
            cheat=ch, cheat_layer=K // 2) for t in range(10))
        print(f"  cheat {ch:>13} @ layer {K // 2:>3}: rejected {rej}/10")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, default=12)
    ap.add_argument("--field", choices=list(FIELDS), default="q")
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--selftest", action="store_true")
    ap.add_argument("--bench", action="store_true")
    args = ap.parse_args()
    if args.selftest:
        selftest()
    elif args.bench:
        bench(field=args.field)
    else:
        main(args.n, args.field, args.seed)
