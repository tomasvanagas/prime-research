#!/usr/bin/env python3
"""
batched_wiring.py -- S503 (batched wiring delegation).

The compressed pi(x) chain (compressed_layer.run_chain), with --delegate, verifies
each Lucy layer's division/affine wiring MLE by the inner GKR chain
inner_verify_div (lucy_dp_delegated_wiring): a degree-2 sum-check #0 over the
2^l state cube, then n backward-sweep layers (each a 2l-var structured
sum-check), then a base check.  The S501 profile found this wiring delegation is
30% of run_chain's wall and 76% of ALL sum-check CALLS (1712): K = pi(sqrt x)
INDEPENDENT inner_verify_div calls, each on a TINY 2^l cube (l = ceil(log2 p_l)),
each paying full Python+numpy per-call dispatch -- op-count-bound, the exact
regime where the fast-Mersenne path (S500) is a NET LOSS.  This is the second
big kernel; the trace batch (S502) widened the first.

Unlike the trace witnesses, the wiring obligations are queried at points
(r_v_l, r_u_l) produced SEQUENTIALLY by each Lucy layer's outer sum-checks, so
they cannot be batched in place.  The restructure (this file) is defer-and-batch:
COLLECT the K obligations (p_l, r_v_l, r_u_l, claim_l, accept_rem_l) and discharge
ALL of them in ONE batched inner protocol.

CONSTRUCTION (data-parallel GKR with a non-uniform, per-layer wiring).  Add
Lk = ceil(log2 K) layer-index variables; pad every layer's 2^{l_l} state cube up
to a common 2^{l_max} (l_max = ceil(log2 sqrt x); MSB-zero-pad -- valid states
stay < p_l, junk stays 0).  The K chains share a single (r_L, rho) trajectory; a
random gamma combines the K claims.  Because the layer index stays a genuine
sum-check variable, the round-poly final factorizes as (product of stacked-table
MLEs), and the per-layer NON-UNIFORM wiring (R_j^l, LT_l depend on p_l and on the
per-layer per-depth bit weights wv_l[j], wu_l[j]) is recomputed by the VERIFIER in
O(K * l) via Sum_l chi_l(r_L) * <per-layer wiring MLE> -- exactly the trace's av
trick generalized.  No dense 2^{2l} table is ever built: each backward layer uses
the S496 STRUCTURED phase split (phase 1 binds sigma', phase 2 binds sigma), and
the layer axis rides along inside both phases, so every prover table is
K * 2^{l_max} WIDE, not 2^{Lk + 2 l_max} dense.

  sum-check #0   : C = Sum_l gamma^l Sum_sig ACC_l(sig) v_n^l(sig)        deg 2
                   over (l, sig);  ACC_l = [sig<p_l] (accept_rem None) or
                   [sig=accept_rem_l].  Reduces C to a claim
                   theta = Sum_l chi_l(r_L0) v_n^l~(rho0).
  backward j     : reduce theta_j = Sum_l w_l v_{j+1}^l~(rho), w_l = chi_l(r_L),
                   to theta_{j-1} = Sum_l chi_l(r_L') v_j^l~(rho').  Two phases:
                   phase 1 over (l, sigma') of W.E.P.G  (W=eq(r_L), E=eq(rho),
                            P=LT_l(sig'), G=Sum_sig R_j^l(sig',sig)LT_l v_j^l);
                   phase 2 over (l, sigma)  of (kappa.chi(r_L1)).Rstar.L.V
                            (Rstar_l(w)=R_j^l~(rho_out,w), L=LT_l, V=v_j^l),
                            kappa = W~ E~ P~ from phase 1.
                   external check ties the phase-2 final to the verifier's
                   recomputed Sum_l chi(r_L2) {R_j^l~, LT_l} and the prover's
                   carried V~ (the absorb-into-next-weight step).
  base           : Sum_l chi_l(r_L) v_0^l~(rho) = (Sum_{l<K} chi_l(r_L)) Prod(1-rho_i)
                   since v_0^l = e_0 for every real layer.

WHY THE WIN (S501).  ONE batched chain instead of K cuts the inner-sum-check CALL
count ~K-fold while WIDENING each fmul's array ~K-fold (2^{l_l} -> ~K 2^{l_max}),
turning op-count-bound into width-bound -- the precondition that makes the fast
Mersenne mulmod (S500) a net win on the chain.  Verifier work stays
O(n * K * l_max) = O~(sqrt x) (the chain already does K layer-reductions, so O(K)
per backward step is within budget).  Soundness ~ (K + n*K*deg*l_max + ...)/q.

What would falsify this: the batched chain rejecting an honest set of wiring
obligations, or ACCEPTING when ANY single obligation is corrupted (wrong claim,
self-consistent lying chain, broken chain vector) -- i.e. the batched accept/reject
disagreeing with the AND of the K independent inner_verify_div; the
fast-Mersenne uint64 path disagreeing with the object reference; or the fmul COUNT
not dropping ~K-fold / per-fmul WIDTH not rising ~K-fold vs the K-independent
baseline.

Usage:
  python3 batched_wiring.py --selftest
  python3 batched_wiring.py --n 14 --field big
  python3 batched_wiring.py --bench --field big
"""

import argparse
import time
from math import isqrt

import numpy as np

import compressed_prover_mult_trace as _cpmt
import lucy_dp_delegated_wiring as _ldw
from compressed_prover_mult_trace import (DEFAULT_Q as Q, BIG_Q, SMALL_Q, FIELDS,
                                          _dt, eq_table, eq_point, sumcheck, fmul)
from lucy_dp_delegated_wiring import (bit_weights, inner_chain_vectors,
                                      r_tilde_eval, lt_const_eval, _modmul_arr,
                                      inner_verify_div)


# ----------------------------------------------------------------------
# small helpers
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


def _ceil_log2(K):
    return 0 if K <= 1 else (K - 1).bit_length()


def _chi(l, r_L, q):
    """Multilinear indicator chi_l(r_L) = prod_j eq(bit_j(l), r_L[j]), MSB-first.
    Matches the binding order of `sumcheck` (most significant variable first)."""
    Lk = len(r_L)
    acc = 1
    for j in range(Lk):
        bit = (l >> (Lk - 1 - j)) & 1
        rj = int(r_L[j]) % q
        acc = acc * (rj if bit else (q + 1 - rj)) % q
    return acc


# ----------------------------------------------------------------------
# the K wiring obligations (what run_chain's large_reduce / small_reduce emit)
# ----------------------------------------------------------------------

def make_obligation(p, n, l_max, r_v, r_u, accept_rem, q=Q, cheat=None):
    """Build one wiring obligation exactly as the delegated chain would query it.
    chain vectors are built on the COMMON l_max cube (MSB-zero-pad: valid states
    stay < p).  The honest claim is read off the chain (= the true wiring MLE):
      accept_rem is None -> W(v,u)=[u=floor(v/p)], value = Sum_{s<p} v_n(s)
      accept_rem = m      -> D(v,u)=[u=floor(v/p) and v mod p = m], value = v_n[m]
    cheat in {None,'claim','lie','vbreak'} corrupts this one obligation."""
    vs = inner_chain_vectors(r_v, r_u, p, l_max, n, q)
    if accept_rem is None:
        claim = sum(int(vs[n][s]) for s in range(p)) % q
    else:
        claim = int(vs[n][accept_rem]) % q

    if cheat == "lie":
        # self-consistent lying chain: corrupt the middle vector and re-propagate
        # forward, so the *claim* still matches v_n (the backward sweep must catch
        # the broken layer identity, not a claim mismatch).
        vs = [v.copy() for v in vs]
        dim = 1 << l_max
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
        if accept_rem is None:
            claim = sum(int(vs[n][s]) for s in range(p)) % q
        else:
            claim = int(vs[n][accept_rem]) % q
    elif cheat == "vbreak":
        # corrupt one final-layer chain entry WITHOUT re-propagating: breaks the
        # last layer identity (caught by the first backward step).
        vs = [v.copy() for v in vs]
        vs[n][0] = (int(vs[n][0]) + 1) % q

    if cheat == "claim":
        claim = (claim + 1) % q          # wrong claim, honest chain

    return dict(p=p, n=n, r_v=list(r_v), r_u=list(r_u), accept_rem=accept_rem,
                claim=claim, vs=vs)


def synth_obligations(x, n, accept_mode="div", q=Q, seed=0, cheat=None,
                      cheat_layer=0):
    """Synthesize the K = pi(sqrt x) wiring obligations a delegated run_chain emits
    (one per prime p <= V).  accept_mode: 'div' (all accept_rem=None, the division
    wiring), 'affine' (all accept_rem=p-1), 'mixed' (alternate).  A single
    obligation `cheat_layer` may be corrupted by `cheat`."""
    V = isqrt(x)
    primes = primes_upto(V)
    l_max = max(max(1, (p - 1).bit_length()) for p in primes)
    rng = np.random.default_rng(seed)
    obs = []
    for li, p in enumerate(primes):
        r_v = [int(rng.integers(0, q)) for _ in range(n)]
        r_u = [int(rng.integers(0, q)) for _ in range(n)]
        if accept_mode == "div":
            ar = None
        elif accept_mode == "affine":
            ar = p - 1
        else:                              # mixed
            ar = None if (li % 2 == 0) else p - 1
        ch = cheat if li == cheat_layer else None
        obs.append(make_obligation(p, n, l_max, r_v, r_u, ar, q, cheat=ch))
    return primes, l_max, obs


# ----------------------------------------------------------------------
# the UNBATCHED baseline: K independent inner_verify_div, AND-ed
# ----------------------------------------------------------------------

def verify_wiring_each(obs, rng, stats, q=Q, structured=True):
    """Exactly what run_chain --delegate does today: one inner_verify_div per
    obligation, AND-ed.  The ground-truth oracle for the batched test."""
    ok = True
    for ob in obs:
        ok = inner_verify_div(ob["claim"], ob["r_v"], ob["r_u"], ob["p"],
                              ob["n"], rng, stats, accept_rem=ob["accept_rem"],
                              structured=structured, q=q) and ok
    return ok


# ----------------------------------------------------------------------
# one batched backward-sweep step (structured phase split, layer axis inside)
# ----------------------------------------------------------------------

def _acc_bits(accept_rem, l_max):
    return [(accept_rem >> (l_max - 1 - i)) & 1 for i in range(l_max)]


def batch_backward_step(obs, r_L_in, rho, claim, j, l_max, Lk, rng, stats, q):
    """Reduce  theta = Sum_l chi_l(r_L_in) v_{j+1}^l~(rho)  to
       theta' = Sum_l chi_l(r_L2)   v_j^l~(rho_in2)  via two O(K 2^l_max) phases.
    Returns (ok, r_L2, rho_in2, claim')."""
    K = len(obs)
    dim = 1 << l_max
    L = 1 << Lk
    N = L * dim
    primes = [ob["p"] for ob in obs]
    wv = [bit_weights(ob["r_v"][j], q) for ob in obs]
    wu = [bit_weights(ob["r_u"][j], q) for ob in obs]
    valid = np.arange(dim)

    # ---- phase 1: bind (l, sigma'); tables W, E, P, G ----
    t0 = time.perf_counter()
    W_lay = eq_table(r_L_in, q)                       # chi_l(r_L_in), length L
    Wt = np.repeat(W_lay, dim)                        # broadcast over sigma'
    Et = np.tile(eq_table(rho, q), L)                 # eq(rho, sigma'), over l
    Pt = np.zeros(N, dtype=_dt(q))
    Gt = np.zeros(N, dtype=_dt(q))
    for li in range(K):
        p = primes[li]
        sl = slice(li * dim, li * dim + dim)
        Pt[sl] = (valid < p).astype(_dt(q))
        G = np.zeros(dim, dtype=_dt(q))
        vbeg = obs[li]["vs"][j][:p]
        wcol = np.arange(p)
        for a in (0, 1):
            for b in (0, 1):
                term = wv[li][a] * wu[li][b] % q
                y0 = 2 * wcol + a - b * p
                m = (y0 >= 0) & (y0 < dim)
                np.add.at(G, y0[m], _modmul_arr(vbeg[m], term, q))
        Gt[sl] = G % q
    stats["t_prover"] += time.perf_counter() - t0
    ok1, rs1, c1, sc1 = sumcheck(claim, {"W": Wt, "E": Et, "P": Pt, "G": Gt},
                                 [(1, ["W", "E", "P", "G"])], 4, rng, q)
    stats["comm"] += 5 * (Lk + l_max)
    if not ok1:
        return False, None, None, None
    r_L1, rho_out = rs1[:Lk], rs1[Lk:]
    kappa = int(sc1["W"]) * int(sc1["E"]) % q * int(sc1["P"]) % q

    # ---- phase 2: bind (l, sigma); tables CHIk=kappa*chi(r_L1), Rstar, L, V ----
    t0 = time.perf_counter()
    chi_rL1 = eq_table(r_L1, q)                        # chi_l(r_L1), length L
    CHIk = np.repeat(_modmul_arr(chi_rL1, kappa, q), dim)
    eqo = eq_table(rho_out, q)                         # length dim
    Rs = np.zeros(N, dtype=_dt(q))
    Lt = np.zeros(N, dtype=_dt(q))
    Vt = np.zeros(N, dtype=_dt(q))
    for li in range(K):
        p = primes[li]
        sl = slice(li * dim, li * dim + dim)
        R = np.zeros(dim, dtype=_dt(q))
        for a in (0, 1):
            for b in (0, 1):
                term = wv[li][a] * wu[li][b] % q
                y0 = 2 * valid + a - b * p
                m = (y0 >= 0) & (y0 < dim)
                R[m] = (R[m] + _modmul_arr(eqo[y0[m]], term, q)) % q
        Rs[sl] = R
        Lt[sl] = (valid < p).astype(_dt(q))
        Vt[sl] = obs[li]["vs"][j]
    stats["t_prover"] += time.perf_counter() - t0
    ok2, rs2, c2, sc2 = sumcheck(c1, {"C": CHIk, "R": Rs, "L": Lt, "V": Vt},
                                 [(1, ["C", "R", "L", "V"])], 4, rng, q)
    stats["comm"] += 5 * (Lk + l_max)
    if not ok2:
        return False, None, None, None
    r_L2, rho_in2 = rs2[:Lk], rs2[Lk:]

    # ---- external check: tie c2 to the verifier-recomputed wiring + carried V ----
    t0 = time.perf_counter()
    W_v = eq_point(r_L_in, r_L1, q)
    E_v = eq_point(rho, rho_out, q)
    P_v = 0
    Rstar_v = 0
    L_v = 0
    for li in range(K):
        p = primes[li]
        c1l = _chi(li, r_L1, q)
        c2l = _chi(li, r_L2, q)
        P_v = (P_v + c1l * lt_const_eval(rho_out, p, l_max, q)) % q
        Rstar_v = (Rstar_v + c2l * r_tilde_eval(rho_in2, rho_out, p, l_max,
                                                wv[li], wu[li], q)) % q
        L_v = (L_v + c2l * lt_const_eval(rho_in2, p, l_max, q)) % q
    chi_eq = eq_point(r_L1, r_L2, q)
    expect = (W_v * E_v % q * P_v % q * chi_eq % q
              * Rstar_v % q * L_v % q * int(sc2["V"])) % q
    ok = (c2 % q) == expect
    stats["t_verifier"] += time.perf_counter() - t0
    if not ok:
        return False, None, None, None
    return True, r_L2, rho_in2, int(sc2["V"])


# ----------------------------------------------------------------------
# the batched wiring protocol (sum-check #0 + n backward steps + base)
# ----------------------------------------------------------------------

def verify_wiring_batched(obs, n, l_max, rng, stats, q=Q):
    """ONE batched inner GKR chain certifying ALL K wiring obligations at once.
    Replaces K calls to inner_verify_div.  Returns ok."""
    K = len(obs)
    Lk = _ceil_log2(K)
    dim = 1 << l_max
    L = 1 << Lk
    Ncube = L * dim
    primes = [ob["p"] for ob in obs]
    valid = np.arange(dim)
    gamma = int(rng.integers(2, q))

    # combined initial claim  C = Sum_l gamma^l claim_l
    C, gp = 0, 1
    for ob in obs:
        C = (C + gp * (ob["claim"] % q)) % q
        gp = gp * gamma % q

    # ---- sum-check #0:  C = Sum_{l,sig} GACC[l,sig] V[l,sig] ----
    t0 = time.perf_counter()
    GACC = np.zeros(Ncube, dtype=_dt(q))
    Vt = np.zeros(Ncube, dtype=_dt(q))
    gp = 1
    for li in range(K):
        p = primes[li]
        ar = obs[li]["accept_rem"]
        sl = slice(li * dim, li * dim + dim)
        acc = (valid < p) if ar is None else (valid == ar)
        GACC[sl] = _modmul_arr(acc.astype(_dt(q)), gp % q, q)
        Vt[sl] = obs[li]["vs"][n]
        gp = gp * gamma % q
    stats["t_prover"] += time.perf_counter() - t0
    ok0, rs0, final0, scal0 = sumcheck(C, {"GACC": GACC, "V": Vt},
                                       [(1, ["GACC", "V"])], 2, rng, q)
    stats["comm"] += 3 * (Lk + l_max)
    if not ok0:
        return False
    r_L0, rho0 = rs0[:Lk], rs0[Lk:]

    t0 = time.perf_counter()
    gacc_v, gp = 0, 1
    for li in range(K):
        p = primes[li]
        ar = obs[li]["accept_rem"]
        acc_mle = (lt_const_eval(rho0, p, l_max, q) if ar is None
                   else eq_point(_acc_bits(ar, l_max), rho0, q))
        gacc_v = (gacc_v + gp * _chi(li, r_L0, q) % q * acc_mle) % q
        gp = gp * gamma % q
    ok = (final0 % q) == gacc_v * int(scal0["V"]) % q
    stats["t_verifier"] += time.perf_counter() - t0
    if not ok:
        return False
    r_L, rho, claim = r_L0, rho0, int(scal0["V"])

    # ---- n backward steps ----
    for j in range(n - 1, -1, -1):
        ok, r_L, rho, claim = batch_backward_step(
            obs, r_L, rho, claim, j, l_max, Lk, rng, stats, q)
        if not ok:
            return False

    # ---- base:  Sum_l chi_l(r_L) v_0^l~(rho), v_0^l = e_0 -> Prod_i (1 - rho_i) ----
    t0 = time.perf_counter()
    base_per = 1
    for ri in rho:
        base_per = base_per * ((q + 1 - int(ri)) % q) % q
    wsum = sum(_chi(li, r_L, q) for li in range(K)) % q
    base = wsum * base_per % q
    res = (claim % q) == base
    stats["t_verifier"] += time.perf_counter() - t0
    return res


# ----------------------------------------------------------------------
# fmul instrumentation (count + array width), for the bench
# ----------------------------------------------------------------------

class _FmulCounter:
    """Count calls and operand widths through BOTH module fmul references (the
    sum-check rounds go through compressed_prover_mult_trace.fmul; the phase-table
    builders' _modmul_arr goes through lucy_dp_delegated_wiring.fmul on BIG_Q)."""

    def __enter__(self):
        self.calls = 0
        self.total_w = 0
        self.max_w = 0
        self._o_cpmt = _cpmt.fmul
        self._o_ldw = _ldw.fmul

        def wrap(a, b, q):
            w = max(np.size(a), np.size(b))
            self.calls += 1
            self.total_w += int(w)
            if w > self.max_w:
                self.max_w = int(w)
            return self._o_cpmt(a, b, q)

        _cpmt.fmul = wrap
        _ldw.fmul = wrap
        return self

    def __exit__(self, *exc):
        _cpmt.fmul = self._o_cpmt
        _ldw.fmul = self._o_ldw

    @property
    def mean_w(self):
        return self.total_w / self.calls if self.calls else 0.0


# ----------------------------------------------------------------------
# tests / experiments
# ----------------------------------------------------------------------

def _new_stats():
    return {"t_prover": 0.0, "t_verifier": 0.0, "comm": 0}


def selftest():
    # 'claim' (wrong wiring value, honest chain) and 'lie' (self-consistent lying
    # chain) are expressible against BOTH the batched protocol and the
    # inner_verify_div baseline (which regenerates its chain from (r_v,r_u), so the
    # corrupted claim mismatches at sum-check #0).  'vbreak' (corrupt a leaf of the
    # PROVER's submitted chain but keep the honest claim) is a BATCHED-ONLY
    # soundness cheat: the baseline never consumes the prover's vectors, so it
    # cannot see it -- it must still be rejected by the batched chain.
    cheats_both = ["claim", "lie"]
    cheats_batched_only = ["vbreak"]
    cheats = cheats_both + cheats_batched_only

    # 1. structural: one batched chain has Lk+l_max-var sum-checks; comm strictly
    #    drops vs the K-fold independent baseline (which is K * O(n l) per layer).
    for n in (8, 10, 12):
        x = (1 << n) - 1
        primes, l_max, obs = synth_obligations(x, n, "div", Q, seed=1)
        K = len(obs)
        Lk = _ceil_log2(K)
        assert (1 << Lk) >= K and (K <= 1 or (1 << (Lk - 1)) < K)
        sb = _new_stats()
        assert verify_wiring_batched(obs, n, l_max, np.random.default_rng(2), sb)
        se = _new_stats()
        assert verify_wiring_each(obs, np.random.default_rng(2), se)
        assert sb["comm"] < se["comm"] or K == 1, (n, sb["comm"], se["comm"])

    # 2. AGREEMENT with the AND of K independent inner_verify_div, over the demo
    #    prime AND BIG_Q (object), for div / affine / mixed wirings:
    #    honest -> both accept; corrupting ANY single obligation (first, middle,
    #    last) with any cheat class -> batched rejects every trial AND the
    #    independent baseline also rejects.
    for q in (Q, BIG_Q):
        for n in (8, 10):
            for mode in ("div", "affine", "mixed"):
                x = (1 << n) - 1
                primes, l_max, obs = synth_obligations(x, n, mode, q, seed=3)
                K = len(obs)
                sb, se = _new_stats(), _new_stats()
                assert verify_wiring_batched(obs, n, l_max,
                                             np.random.default_rng(4), sb, q), \
                    (q, n, mode)
                assert verify_wiring_each(obs, np.random.default_rng(4), se, q), \
                    (q, n, mode)
                for cl in sorted(set([0, K // 2, K - 1])):
                    for ch in cheats:
                        _, lm2, oc = synth_obligations(x, n, mode, q, seed=3,
                                                       cheat=ch, cheat_layer=cl)
                        rej = sum(not verify_wiring_batched(
                            oc, n, lm2, np.random.default_rng(50 + t),
                            _new_stats(), q) for t in range(4))
                        assert rej == 4, (q, n, mode, ch, cl, rej)
                        # the protocol it replaces ALSO rejects, for the cheats it
                        # can see (a corrupted claim mismatches its honest chain).
                        if ch in cheats_both:
                            assert not verify_wiring_each(
                                oc, np.random.default_rng(7), _new_stats(), q), \
                                (q, n, mode, ch, cl)

    # 3. K=1 edge (x=7: only prime 2 <= V=2; Lk=0, no layer bits)
    x = 7
    primes, l_max, obs = synth_obligations(x, 3, "div", Q, seed=5)
    assert len(obs) == 1 and _ceil_log2(1) == 0
    assert verify_wiring_batched(obs, 3, l_max, np.random.default_rng(6),
                                 _new_stats())
    _, lm2, oc = synth_obligations(x, 3, "div", Q, seed=5, cheat="claim",
                                   cheat_layer=0)
    assert sum(not verify_wiring_batched(
        oc, 3, lm2, np.random.default_rng(10 + t), _new_stats())
        for t in range(4)) == 4

    # 4. fast Mersenne uint64 path == object reference, bit-for-bit (same
    #    accept/reject, honest + every cheat), for BIG_Q.
    saved = _cpmt.FAST_BIG
    try:
        for n in (8, 10):
            x = (1 << n) - 1
            for ch in [None] + cheats:
                cl = 1
                _cpmt.FAST_BIG = False
                _, lm, oo = synth_obligations(x, n, "div", BIG_Q, seed=8,
                                              cheat=ch, cheat_layer=cl)
                ro = verify_wiring_batched(oo, n, lm, np.random.default_rng(9),
                                           _new_stats(), BIG_Q)
                _cpmt.FAST_BIG = True
                _, lm, of = synth_obligations(x, n, "div", BIG_Q, seed=8,
                                              cheat=ch, cheat_layer=cl)
                rf = verify_wiring_batched(of, n, lm, np.random.default_rng(9),
                                           _new_stats(), BIG_Q)
                assert ro == rf, (n, ch)
                assert ro == (ch is None), (n, ch)
    finally:
        _cpmt.FAST_BIG = saved

    # 5. the structural win: ONE batched chain cuts the fmul CALL count and
    #    WIDENS each fmul, vs the K-fold independent baseline (object dtype).
    x = (1 << 12) - 1
    primes, l_max, obs = synth_obligations(x, 12, "div", BIG_Q, seed=1)
    K = len(obs)
    with _FmulCounter() as cu:
        verify_wiring_each(obs, np.random.default_rng(1), _new_stats(), BIG_Q)
    with _FmulCounter() as cb:
        verify_wiring_batched(obs, 12, l_max, np.random.default_rng(1),
                              _new_stats(), BIG_Q)
    assert cb.calls < cu.calls, (cb.calls, cu.calls)
    assert cb.mean_w > cu.mean_w, (cb.mean_w, cu.mean_w)

    print("selftest OK")


def _wall(fn):
    with _FmulCounter() as c:
        t0 = time.perf_counter()
        fn()
        wall = time.perf_counter() - t0
    return wall, c.calls, c.mean_w, c.max_w


def bench(field="big", ns=(8, 10, 12, 14), seed=1):
    """The wiring delegation in isolation (30% of run_chain wall, 76% of its
    sum-check calls, S501): K independent inner_verify_div (today's chain) vs ONE
    batched chain.  Reports fmul COUNT (should drop ~K-fold), per-fmul mean WIDTH
    (should rise ~K-fold), and wall, on the object path and -- for --field big --
    the fast-Mersenne uint64 path (where the width win should flip S500's
    'fast path loses on the chain' to a net gain)."""
    q = FIELDS[field]
    paths = [("object", False), ("fast", True)] if field == "big" else [("", False)]
    print(f"batched wiring delegation, field = {field} (q = {q})")
    print(f"{'n':>3} {'K':>4} {'lmax':>4} {'Lk':>3} {'path':>7} || "
          f"{'unb_calls':>10} {'bat_calls':>10} {'c_drop':>6} || "
          f"{'unb_w':>7} {'bat_w':>7} {'w_rise':>6} || "
          f"{'unb_ms':>9} {'bat_ms':>9} {'speedup':>7}")
    save = _cpmt.FAST_BIG
    try:
        for n in ns:
            x = (1 << n) - 1
            for label, fast in paths:
                if field == "big":
                    _cpmt.FAST_BIG = fast
                primes, l_max, obs = synth_obligations(x, n, "div", q, seed=seed)
                K = len(obs)
                Lk = _ceil_log2(K)
                wu_, cu, mu, _ = _wall(lambda: verify_wiring_each(
                    obs, np.random.default_rng(seed), _new_stats(), q))
                wb, cb, mb, _ = _wall(lambda: verify_wiring_batched(
                    obs, n, l_max, np.random.default_rng(seed), _new_stats(), q))
                print(f"{n:>3} {K:>4} {l_max:>4} {Lk:>3} {label:>7} || "
                      f"{cu:>10} {cb:>10} {cu/max(cb,1):>5.1f}x || "
                      f"{mu:>7.1f} {mb:>7.1f} {mb/max(mu,1e-9):>5.1f}x || "
                      f"{wu_*1000:>9.1f} {wb*1000:>9.1f} "
                      f"{wu_/max(wb,1e-9):>6.2f}x")
    finally:
        _cpmt.FAST_BIG = save


def main(n, field, seed):
    q = FIELDS[field]
    if field == "big":
        _cpmt.FAST_BIG = True
    x = (1 << n) - 1
    primes, l_max, obs = synth_obligations(x, n, "div", q, seed=seed)
    K = len(obs)
    Lk = _ceil_log2(n if False else K)
    sb = _new_stats()
    ok = verify_wiring_batched(obs, n, l_max, np.random.default_rng(seed), sb, q)
    print(f"x = 2^{n}-1 = {x}, V = {isqrt(x)}, l_max = {l_max}, "
          f"K = pi(V) = {K} wiring obligations")
    print(f"field = {field} (q = {q}); batched cube = 2^({Lk}+{l_max}) "
          f"per phase (one chain; unbatched ran {K} chains)")
    print(f"BATCHED wiring chain (all {K} obligations at once): "
          f"{'ACCEPTED' if ok else 'REJECTED'}")
    assert ok
    se = _new_stats()
    assert verify_wiring_each(obs, np.random.default_rng(seed), se, q)
    print(f"  comm batched = {sb['comm']} field elems   "
          f"(unbatched comm = {se['comm']} = K independent chains)")
    for ch in ["claim", "lie", "vbreak"]:
        _, lm2, oc = synth_obligations(x, n, "div", q, seed=seed, cheat=ch,
                                       cheat_layer=K // 2)
        rej = sum(not verify_wiring_batched(
            oc, n, lm2, np.random.default_rng(seed + 70 + t), _new_stats(), q)
            for t in range(10))
        print(f"  cheat {ch:>7} @ obligation {K // 2:>3}: rejected {rej}/10")


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
