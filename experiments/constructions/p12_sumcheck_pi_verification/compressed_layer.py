#!/usr/bin/env python3
"""
compressed_layer.py — S493/S494, compressed-prover steps 2 & 3.

Step 1 (compressed_prover_mult_trace.py) built, in isolation, the batched
multiplication-trace sum-check that certifies a small-side value lookup
C = sum_d eq(rho,d) S~(floor(X/(d p))) with a polylog verifier that never
divides. Step 2 (S493) integrated it into ONE large-side Lucy layer over the
compressed O(sqrt x)-size state. Step 3 (S494, THIS file) adds the SMALL-side
layer and CHAINS all K layers into an end-to-end verification of pi(x) over
the compressed state -- the gate to a production Õ(sqrt x) prover (the demo
prover of lucy_dp_verification materializes the full 2^n table, Õ(x^{3/2})).

The compressed Lucy_Hedgehog key-value DP, V = floor(sqrt x), nb = bits(V):
  small[v]  = S(v)              value-addressed,  v in [1, V]   (index = value)
  large[e]  = S(floor(x/(e+1))) d-addressed,      e in [0, V-1] (index e = d-1)
and pi(x) = large[0] after sieving by every prime <= V. Both sides are padded
to the nb-cube of size 2^nb (V < 2^nb <= 2V); padding entries are 0 every
layer.

The layer recurrence S_i(w) = S_{i-1}(w) - [w>=p^2]*(S_{i-1}(floor(w/p))-(i-1)).

LARGE side (key w = floor(x/(e+1)), step 2): floor(w/p)=floor(x/((e+1)p)).
With ep=(e+1)p, ep<=V is an AFFINE large-index map e->ep-1 (long-division
automaton waff_eval, no trace); ep>V is a SMALL value looked up by value
WITHOUT the verifier dividing -- step 1's trace primitive. Phase-A folds
validity AND the p^2 gate into ONE threshold B(e)=[e<=Bcut]; phase-B splits
g1=g1_aff+g1_trace and dispatches per key. Outputs: a large point-claim
(line-batched s1@r_v + s2_aff@r_u) and a SMALL point-claim s_B@r_B (carried).

SMALL side (key w = v, step 3, NEW): floor(w/p)=floor(v/p) stays ENTIRELY
within the small side (value->value) -- the existing division automaton
w_div_eval / lucy_dp_verification.verify_layer phase-B over the nb-cube. The
gate is a BAND B(v)=[p^2<=v<=V]: the upper bound [v<=V] is load-bearing in
general (it keeps the recurrence off the padding columns whose update would
otherwise corrupt the zero padding); when V=2^nb-1 (x=2^n-1, n even) the band
collapses to the bare comparator [v>=p^2] of the NEXT-ACTION note -- so the
"comparator suffices, no validity fold" claim holds exactly in the no-padding
case and the band is its general-V refinement. The small layer emits TWO
point-claims on S_{i-1}^small: s1@r_v (the v itself) and s2@r_u (floor(v/p)).

CHAIN (step 3). Per layer we carry a width-2 claim vector
(large-claim about S_i^large, small-claim about S_i^small). Start: layer K
has only large-claim S_K^large(e=0)=pi(x) (no S_K^small is referenced).
Running layer i:
  * large_reduce consumes the large-claim -> a large-claim about S_{i-1}^large
    AND a small-claim s_B about S_{i-1}^small (its trace region);
  * small_reduce (skipped at i=K) consumes the small-claim -> two small-claims
    about S_{i-1}^small.
The same-table claims about S_{i-1}^small (s_B + the two from small_reduce)
are folded by sequential LINE restriction into one; the large-claim is its own
line-batch. Base: open S_0 on both sides. Soundness: a self-consistent liar
(corrupt one DP entry, re-propagate upward) is caught at the layer where the
corruption meets the honest layer below -- phase-A's round-1 sum on the honest
S_{i-1} tables disagrees with the corrupted carried claim.

Leaf openings of committed tables (incl. the two S_0 bases) are closed here by
direct MLE evaluation -- the stand-in for the polynomial-commitment / outer-
protocol line batching, identical to steps 1-2's scope. The S_0 bases have no
cheap closed form (floor(x/(e+1)) on the large side; padded max(v-1,0) on the
small side), so they are opened directly: a ONE-TIME O(sqrt x), within the
Õ(sqrt x) budget. The reduction logic is real: every per-layer check is a
genuine MLE-consistency test against true wiring (eq, band/threshold
comparators, division & affine automata); honest accepts, every cheating
prover is rejected.

Sizes (the point of the line): small & large cubes 2^nb ~ sqrt x; the trace
constraint test runs on the e-cube with Lv=O(log x) bit tables -> prover
Õ(sqrt x); verifier touches only nb-length points and O(nb)-state automata
(except waff's O(nb*p), delegatable to O(nb*log p)).

What would falsify this: the compressed DP disagreeing with a sieve; an honest
layer or honest chain rejected; output/base claims not equal to the true MLEs;
any cheating prover (corrupted DP claim/table, wrong region split, wrong trace
quotient, wrong affine image, corrupted small-side wiring, tampered cross-side
carry, wrong claimed pi) accepted above the field soundness bound; or verifier
work scaling with 2^n instead of nb.

Usage:
  python3 compressed_layer.py --selftest
  python3 compressed_layer.py --n 12 --cheat-trials 20
  python3 compressed_layer.py --bench
"""

import argparse
import time
from math import isqrt

import numpy as np

# Field-parameterised helpers (q kwarg, _dt/_asum dtype machinery) come from
# compressed_prover_mult_trace and default to q=Q (= DEFAULT_Q = 2^31-1), so the
# whole chain stays BIT-IDENTICAL at the demo prime while every layer/wiring can
# now run over a larger-characteristic prime (BIG_Q = 2^61-1, object dtype) --
# the field lift (S498) threaded end-to-end through run_chain.  ge_const_eval /
# w_div_eval are the q-defaulted automata from lucy_dp_verification.
import compressed_prover_mult_trace as _cpmt
from compressed_prover_mult_trace import (DEFAULT_Q as Q, BIG_Q, SMALL_Q,
                                          FIELDS, _dt, _asum, lagrange_eval,
                                          eq_table, eq_point, mle_eval, sumcheck,
                                          build_witness, verify_constraints,
                                          forge_alias, rematerialize, fmul,
                                          fast_big)
from lucy_dp_verification import ge_const_eval, w_div_eval
from lucy_dp_delegated_wiring import inner_verify_div, inner_chain_vectors
import leaf_open as _lo


# ----------------------------------------------------------------------
# verifier large-table-eval accounting (S507): the op-count curve
# ----------------------------------------------------------------------
# The ONLY verifier operation whose per-call cost is Theta(2^nb) = Theta(sqrt x)
# is the direct MLE evaluation of a committed sqrt x-size table (mle_eval -- a
# "leaf opening").  Every OTHER verifier step in the chain -- the sum-check
# round-poly checks, the eq_point / le_eval / band_eval point evals, the
# delegated wiring carry chain (O(nb log p)), and the batched-discharge MLE
# recomputes (O(K polylog)) -- costs O(polylog) or O(K polylog) and so is bounded
# by Õ(sqrt x) over the whole K = pi(sqrt x)-layer chain.  Therefore the LEADING
# term of the whole-chain verifier cost IS the count of these large-table evals,
# weighted by table size.  _acct_vleaf tallies them, split into the PER-LAYER
# critical path (which scales with K ~ sqrt x / ln) vs the ONE-TIME terms (the
# two S_0 base closes), so a single run exposes whether the leading term is
# Theta(x) (per-layer leaf opens: K * sqrt x) or Theta(sqrt x) (one-time only).
# Keys: vleaf_ops_pl / vleaf_ops_ot (summed table sizes = field-ops touched, up
# to a constant), vleaf_cnt_pl / vleaf_cnt_ot (eval counts).  The wall-clock
# t_verifier is a misleading proxy (the dominant sum-check FOLDING is untimed
# prover work, S501, and vectorised numpy hides the asymptotics, S506) -- this
# op count is the honest headline.

def _acct_vleaf(stats, size, onetime=False):
    """Account ONE verifier large-table (O(size)) evaluation toward the op-count
    curve.  size = len(table), the field-ops a direct mle_eval touches up to a
    constant.  onetime=True buckets it as a chain-level one-time term (the S_0
    base closes); otherwise it is per-layer critical-path work (scales with K)."""
    suf = "ot" if onetime else "pl"
    stats["vleaf_ops_" + suf] = stats.get("vleaf_ops_" + suf, 0) + int(size)
    stats["vleaf_cnt_" + suf] = stats.get("vleaf_cnt_" + suf, 0) + 1


# ----------------------------------------------------------------------
# the compressed Lucy_Hedgehog DP (prover substrate), all layers, exact
# ----------------------------------------------------------------------

def compressed_lucy(x):
    """Return (primes, small_layers, large_layers).
    small_layers[i][v] = S_i(v) for v in [0,V];  small_layers[0] = S_0.
    large_layers[i][e] = S_i(floor(x/(e+1))) for e in [0,V).
    pi(x) = large_layers[-1][0]."""
    V = isqrt(x)
    small = [max(v - 1, 0) for v in range(V + 1)]
    large = [(x // (e + 1)) - 1 for e in range(V)]
    small_layers, large_layers, primes = [small[:]], [large[:]], []
    for p in range(2, V + 1):
        if small[p] - small[p - 1] > 0:                     # p is prime
            primes.append(p)
            sp, p2 = small[p - 1], p * p
            ns, nl = small[:], large[:]
            for v in range(p2, V + 1):
                ns[v] = small[v] - (small[v // p] - sp)
            for e in range(V):
                val = x // (e + 1)
                if val >= p2:
                    ep = (e + 1) * p
                    sval = large[ep - 1] if ep <= V else small[val // p]
                    nl[e] = large[e] - (sval - sp)
            small, large = ns, nl
            small_layers.append(small[:])
            large_layers.append(large[:])
    return primes, small_layers, large_layers


def corrupt_compressed(primes, sm_layers, lg_layers, x, i0):
    """Self-consistent liar: corrupt large_layers[i0][0] (the pi-carrying
    slot, = S_{i0}(x)) and re-propagate layers i0+1..K so the corrupted tables
    are mutually consistent. The small side is unchanged (its recurrence never
    reads the large side), so the corruption is large-only and must be caught
    at layer i0 where it meets the honest layer i0-1 below."""
    V = isqrt(x)
    K = len(primes)
    sm = [s[:] for s in sm_layers]
    lg = [l[:] for l in lg_layers]
    lg[i0][0] = lg[i0][0] + 1
    small, large = sm[i0][:], lg[i0][:]
    for j in range(i0 + 1, K + 1):
        p = primes[j - 1]
        sp, p2 = small[p - 1], p * p
        ns, nl = small[:], large[:]
        for v in range(p2, V + 1):
            ns[v] = small[v] - (small[v // p] - sp)
        for e in range(V):
            val = x // (e + 1)
            if val >= p2:
                ep = (e + 1) * p
                sval = large[ep - 1] if ep <= V else small[val // p]
                nl[e] = large[e] - (sval - sp)
        small, large = ns, nl
        sm[j], lg[j] = small[:], large[:]
    return sm, lg


# ----------------------------------------------------------------------
# verifier-side O(nb) threshold MLEs and the O(nb*p) affine-wiring MLE
# ----------------------------------------------------------------------

def le_eval(r, M, n, q=Q):
    """MLE of [e <= M] on n-bit e at point r. O(n)."""
    if M < 0:
        return 0
    if M >= (1 << n) - 1:
        return 1
    return (q + 1 - ge_const_eval(r, M + 1, n, q)) % q


def band_eval(r, lo, hi, n, q=Q):
    """MLE of [lo <= e <= hi].  [lo<=e<=hi] = [e<=hi] - [e<=lo-1] pointwise,
    and the MLE is linear over that pointwise difference. O(n).  An EMPTY band
    (lo > hi, e.g. the small-side gate when p^2 > V) is the all-zero indicator,
    whose MLE is 0 -- and le(hi)-le(lo-1) is NOT 0 there, so guard it."""
    if lo > hi:
        return 0
    return (le_eval(r, hi, n, q) - le_eval(r, lo - 1, n, q)) % q


def waff_eval(r_v, r_u, p, ecut, n, q=Q):
    """MLE of the affine-region wiring W_aff(e, e') = [e <= ecut] * [e' = p*e
    + (p-1)] at (r_v over e, r_u over e').  e' = p*e+(p-1) <=> floor(e'/p)=e
    AND e' mod p = p-1, so a long-division automaton on e' (state = running
    remainder in [0,p)) whose emitted quotient bits must equal e, accepting
    on final remainder p-1, crossed with a 3-state comparator e<=ecut.
    Bits MSB-first. VERIFIER: O(n*p)."""
    st = {(0, 0): 1}                       # (remainder, cmp) -> weight; cmp 0=eq 1=gt 2=lt
    for j in range(n):
        cbit = (ecut >> (n - 1 - j)) & 1
        wv = [(q + 1 - int(r_v[j])) % q, int(r_v[j]) % q]   # e bit
        wu = [(q + 1 - int(r_u[j])) % q, int(r_u[j]) % q]   # e' bit
        new = {}
        for (rem, cmp), w in st.items():
            if w == 0:
                continue
            for epb in (0, 1):
                t = 2 * rem + epb
                qb, rem2 = divmod(t, p)
                if qb > 1:                 # quotient bit out of range -> dead
                    continue
                cmp2 = (cmp if cmp != 0 else
                        (1 if qb > cbit else 2 if qb < cbit else 0))
                wgt = w * wu[epb] % q * wv[qb] % q          # eb := qb (forced)
                key = (rem2, cmp2)
                new[key] = (new.get(key, 0) + wgt) % q
        st = new
    return sum(w for (rem, cmp), w in st.items()
               if rem == p - 1 and cmp in (0, 2)) % q       # e <= ecut


def verify_waff_value(waffv, r_v, r_u, p, ecut, nb, rng, stats,
                      lie_chain=False, structured=False, q=Q, defer=None):
    """Prove the affine-wiring scalar
        waffv = waff_eval(r_v, r_u) = sum_{e<=ecut} eq(r_v,e) eq(r_u,(e+1)p-1)
    with an O(nb log p) verifier that NEVER divides -- the delegated drop-in
    for the O(nb*p) waff_eval automaton.  Decomposition:
      waffv = sum_e EQ(r_v,e) * C(e) * Phi(e),  C(e)=[e<=ecut],
              Phi(e) = eq(r_u, (e+1)p-1) = the pure affine-image relation.
    An inner sum-check over the e-cube folds the comparator C (its MLE is the
    O(nb) le_eval -- the 3-state comparator KEPT AS-IS) and the eq(r_v,.),
    leaving a single claim Phi~(r_e) = D~(r_e, r_u) for the affine-image
    relation D(e,e')=[e'=p*e+(p-1)].  That claim is delegated to the SAME
    carry-chain template (lucy_dp_delegated_wiring.inner_verify_div) with the
    dividend = e' (point r_u), the quotient = e (point r_e), and acceptance
    restricted to final remainder == p-1.  Returns ok.

    structured (S496 drop-in): when True the delegated inner_verify_div runs
    its backward chain sweep in O(nb*p) prover work (chain_layer_reduce_
    structured) instead of the dense O(nb*p^2) product-of-5-tables sweep, with
    a BIT-IDENTICAL transcript -- accept/reject and comm are unchanged.  Only
    affects the delegated path; the comparator fold above is already O(2^nb).

    defer (S503): when a list is passed, the inner affine-image chain is NOT run
    inline; instead the obligation (p, dividend=r_u, quotient=r_e, accept_rem=p-1,
    claim=phi_claim, lie_chain) is appended and the function returns True after the
    comparator fold.  run_chain(batch_wiring=True) collects these across all layers
    and discharges them in ONE batched_wiring.verify_wiring_obligations call.  The
    O(2^nb) comparator fold is p-independent and stays per-layer."""
    Dml = 1 << nb
    t0 = time.perf_counter()
    eqrv = eq_table(r_v, q)
    eqru = eq_table(r_u, q)
    arange = np.arange(Dml, dtype=np.int64)
    img = (arange + 1) * p - 1                  # image of each e under the map
    vmask = img < Dml
    Phi = np.zeros(Dml, dtype=_dt(q))           # Phi[e] = eq(r_u, (e+1)p-1)
    Phi[arange[vmask]] = eqru[img[vmask]]
    Cind = (arange <= ecut).astype(_dt(q))      # [e <= ecut]
    stats["t_prover"] += time.perf_counter() - t0

    # inner sum-check (deg 3): waffv = sum_e EQrv(e) Cind(e) Phi(e)
    okI, r_e, finalI, scalI = sumcheck(waffv,
                                       {"EQ": eqrv.copy(), "C": Cind.copy(),
                                        "PHI": Phi.copy()},
                                       [(1, ["EQ", "C", "PHI"])], 3, rng, q)
    stats["comm"] += 4 * nb
    if not okI:
        return False
    t0 = time.perf_counter()
    phi_claim = scalI["PHI"]
    ok = (finalI % q) == (eq_point(r_v, r_e, q) * le_eval(r_e, ecut, nb, q) % q
                          * phi_claim) % q
    stats["t_verifier"] += time.perf_counter() - t0
    if not ok:
        return False
    # delegate phi_claim = D~(r_e over quotient e, r_u over dividend e'),
    # D(e',e)=[e=floor(e'/p) and e' mod p = p-1]  <=>  e'=p*e+(p-1).
    if defer is not None:                       # batch_wiring: collect, discharge later
        defer.append((p, list(r_u), list(r_e), p - 1, int(phi_claim) % q,
                      lie_chain))
        return True
    return inner_verify_div(phi_claim, r_u, r_e, p, nb, rng, stats,
                            accept_rem=p - 1, lie=lie_chain,
                            structured=structured, q=q)


# ----------------------------------------------------------------------
# region verifiers (large side)
# ----------------------------------------------------------------------

def verify_affine_region(S_large, r_v, claimed, p, ecut, nb, rng, stats,
                         cheat=None, delegate=False, structured=False, q=Q,
                         defer=None, pcs=False):
    """g1_aff = sum_{e<=ecut} eq(r_v,e) S_large((e+1)p-1).  Routed by image
    index e': g1_aff = sum_{e'} Waff[e'] S_large(e').  The wiring scalar
    Waff~(r_u) = waff_eval(r_v, r_u) is either evaluated by the O(nb*p)
    automaton (delegate=False) or PROVEN by the O(nb*log p) inner protocol
    verify_waff_value (delegate=True).  cheat in
    {None,'waff','waff_forge','waff_chain'}.  structured passes the S496
    O(nb*p) chain prover down to verify_waff_value (delegate path only).
    Returns (ok, r_u, s2)."""
    Dml = 1 << nb
    t0 = time.perf_counter()
    eqrv = eq_table(r_v, q)
    Waff = np.zeros(Dml, dtype=_dt(q))
    for e in range(ecut + 1):
        ep1 = (e + 1) * p - 1
        if ep1 < Dml:
            Waff[ep1] = (int(Waff[ep1]) + int(eqrv[e])) % q
    if cheat in ("waff", "waff_forge"):                     # wrong affine image
        bad = (p) % Dml
        Waff[bad] = (int(Waff[bad]) + 1) % q
    stats["t_prover"] += time.perf_counter() - t0
    okB, r_u, finalB, scal = sumcheck(claimed,
                                      {"W": Waff.copy(), "S": S_large.copy()},
                                      [(1, ["W", "S"])], 2, rng, q)
    stats["comm"] += 3 * nb
    if not okB:
        return False, None, None
    s2 = scal["S"]
    # pcs (S505): s2 = S~_large(r_u) is line-batched into the large carried claim
    # (new_z, new_claim) by the caller and grounded transitively at the base, so
    # the per-layer O(2^nb) close is redundant and skipped.  Default: close it.
    if not pcs:
        t0 = time.perf_counter()
        ok_open = s2 == mle_eval(S_large, r_u, q)
        _acct_vleaf(stats, len(S_large))                    # per-layer leaf open
        stats["t_verifier"] += time.perf_counter() - t0
        if not ok_open:
            return False, None, None

    if not delegate:
        t0 = time.perf_counter()
        waffv = waff_eval(r_v, r_u, p, ecut, nb, q)         # O(nb*p) automaton
        ok = (finalB % q) == waffv * s2 % q
        stats["t_verifier"] += time.perf_counter() - t0
        return (ok, r_u, s2) if ok else (False, None, None)

    # delegate: prover supplies waffv (computed Õ(sqrt x), NO automaton), the
    # verifier binds it to finalB and PROVES it via verify_waff_value.
    t0 = time.perf_counter()
    if cheat == "waff_forge" and int(s2) % q != 0:
        waffv = finalB % q * pow(int(s2) % q, q - 2, q) % q  # forge to pass algebra
    else:
        eqru = eq_table(r_u, q)
        e_idx = np.arange(ecut + 1, dtype=np.int64)
        img = (e_idx + 1) * p - 1
        m = img < Dml
        waffv = _asum(fmul(eqrv[e_idx[m]], eqru[img[m]], q), q)
    stats["t_prover"] += time.perf_counter() - t0
    t0 = time.perf_counter()
    ok = (finalB % q) == waffv * s2 % q
    stats["t_verifier"] += time.perf_counter() - t0
    if not ok:
        return False, None, None
    okw = verify_waff_value(waffv, r_v, r_u, p, ecut, nb, rng, stats,
                            lie_chain=(cheat == "waff_chain"),
                            structured=structured, q=q, defer=defer)
    return (True, r_u, s2) if okw else (False, None, None)


def verify_trace_region(W, S_small, r_v, claimed, lo, hi, nb, rng, stats, q=Q,
                        skip_constraints=False, pcs=False, ub_defer=None,
                        layer_idx=0):
    """g1_trace = sum_{lo<=e<=hi} eq(r_v,e) S_small(floor(x/((e+1)p))).
    Certify every u_e by the step-1 trace zero-test, then route the masked
    small-side lookup. The band mask T(e)=[lo<=e<=hi] zeroes the affine /
    padding columns (whose u_e exceed V) out of the routing.  The lookup is
    over nb-bit values (<= V), taking the LOW nb bits of the certified u_e,
    so omega / S_small stay sqrt x.  Returns (ok, r_B, s_B).

    skip_constraints (S502): the per-layer trace zero-test is omitted because
    the caller has already certified THIS witness in ONE batched zero-test over
    all K layers (batched_trace.verify_constraints_batched) -- the S501 profile's
    53%-of-wall kernel, run K-fold cheaper.  The masked-lookup routing below
    (which still pins g1_trace to the certified u_e via B1/B2) is unchanged."""
    if not skip_constraints and not verify_constraints(W, rng, stats, q):
        return False, None, None
    D, Lv = W["D"], W["Lv"]
    e = np.arange(D, dtype=np.int64)
    t0 = time.perf_counter()
    Tband = ((e >= lo) & (e <= hi)).astype(_dt(q))
    eqr = eq_table(r_v, q)
    weight = (eqr * Tband) % q
    mask = (1 << nb) - 1
    vlow = (W["u"] & mask).astype(np.int64)                 # low nb bits of u_e
    if fast_big(q):                                         # < D summands overflow u64
        acc = np.zeros(1 << nb, dtype=object)
        np.add.at(acc, vlow, weight.astype(object))
        omega = (acc % q).astype(np.uint64)
    else:
        omega = np.zeros(1 << nb, dtype=_dt(q))
        np.add.at(omega, vlow, weight)                      # < D summands < q
        omega %= q
    stats["t_prover"] += time.perf_counter() - t0

    # B1 (deg 2, v-cube): claimed = sum_v S_small(v) omega(v)
    okB1, r_B, finalB1, sB1 = sumcheck(claimed,
                                       {"S": S_small.copy(), "W": omega.copy()},
                                       [(1, ["S", "W"])], 2, rng, q)
    stats["comm"] += 3 * nb
    if not okB1:
        return False, None, None
    # pcs (S505): s_B = sB1["S"] = S~_small(r_B) is folded into the small carried
    # claim by the caller (batch_on_table) and grounded at the base, so its
    # per-layer O(2^nb) close is redundant and skipped; the wiring-consistency
    # check finalB1 == s_B * omega~(r_B) stays (grounds g1_trace to the lookup).
    t0 = time.perf_counter()
    okf = (finalB1 % q) == sB1["S"] * sB1["W"] % q
    if not pcs:
        okf = okf and sB1["S"] == mle_eval(S_small, r_B, q)  # opening of S_small
        _acct_vleaf(stats, len(S_small))                     # per-layer leaf open
    stats["t_verifier"] += time.perf_counter() - t0
    if not okf:
        return False, None, None

    # B2 (deg nb+2, e-cube): sB1[W] = sum_e eq(r_v,e) T(e) prod_k EB_k(e),
    #   EB_k = eq(r_B[k], Ub_{Lv-nb+k})  (the LOW nb bits of the SAME u_e).
    t0 = time.perf_counter()
    tabs = {"EQ": eq_table(r_v, q), "T": Tband.copy()}
    ebn, jidx = [], []
    for k in range(nb):
        j = Lv - nb + k
        jidx.append(j)
        rb = int(r_B[k]) % q
        Ub = W["tabs"][f"Ub{j}"].astype(_dt(q))
        tabs[f"EB{k}"] = (((1 - rb) % q) + ((2 * rb - 1) % q) * Ub) % q
        ebn.append(f"EB{k}")
    stats["t_prover"] += time.perf_counter() - t0
    okB2, r_C, finalB2, _ = sumcheck(sB1["W"], tabs,
                                     [(1, ["EQ", "T"] + ebn)], nb + 2, rng, q)
    stats["comm"] += (nb + 3) * nb
    if not okB2:
        return False, None, None
    if ub_defer is None:
        t0 = time.perf_counter()
        expect = eq_point(r_v, r_C, q) * band_eval(r_C, lo, hi, nb, q) % q
        for k in range(nb):
            rb = int(r_B[k]) % q
            ub = mle_eval(W["tabs"][f"Ub{jidx[k]}"], r_C, q)    # opening of Ub_j
            _acct_vleaf(stats, len(W["tabs"][f"Ub{jidx[k]}"]))  # per-layer leaf open
            expect = expect * (((1 - rb) % q) + ((2 * rb - 1) % q) * ub) % q
        ok = (finalB2 % q) == expect
        stats["t_verifier"] += time.perf_counter() - t0
        stats["ub_leaf_v"] = stats.get("ub_leaf_v", 0) + nb      # O(2^nb) VERIFIER evals
    else:
        # S506: defer the nb Ub openings -- the LAST per-layer O(2^nb) verifier
        # term.  The prover supplies the Ub_{jidx[k]}~(r_C) scalars (t_prover); the
        # verifier folds them into `expect` in O(nb) here and records
        # (layer_idx, r_C, ub_claims) for ONE batched discharge after the layer
        # loop (batched_trace.verify_ub_openings_batched), against the SAME stacked
        # Ub cube the trace zero-test commits.  A forged ub_claim is caught there.
        t0 = time.perf_counter()
        ubs = [int(mle_eval(W["tabs"][f"Ub{jidx[k]}"], r_C, q)) % q
               for k in range(nb)]                          # prover-supplied openings
        stats["t_prover"] += time.perf_counter() - t0
        stats["ub_leaf_p"] = stats.get("ub_leaf_p", 0) + nb     # now PROVER-side
        t0 = time.perf_counter()
        expect = eq_point(r_v, r_C, q) * band_eval(r_C, lo, hi, nb, q) % q
        for k in range(nb):
            rb = int(r_B[k]) % q
            expect = expect * (((1 - rb) % q) + ((2 * rb - 1) % q) * ubs[k]) % q
        ok = (finalB2 % q) == expect
        stats["t_verifier"] += time.perf_counter() - t0
        ub_defer.append((layer_idx, list(r_C), ubs))
    return ok, r_B, sB1["S"]


# ----------------------------------------------------------------------
# claim aggregation: fold same-table point-claims by line restriction
# ----------------------------------------------------------------------

def line_batch_pair(S, pt0, val0, pt1, val1, nb, rng, stats, q=Q, pcs=False):
    """Reduce val0=S~(pt0), val1=S~(pt1) to ONE claim (new_pt, new_claim) via
    the degree-nb line restriction h(t)=S~((1-t)pt0+t pt1). Endpoint checks
    anchor a false input; new_claim is closed by mle_eval (the PCS stand-in).

    pcs (S505): replace the line restriction + its closing mle_eval with a real
    sum-check OPENING (leaf_open.open_batch over the two claims). The returned
    (new_pt, new_claim) is then the sum-check RESIDUAL -- NOT closed here, but
    threaded by the caller as a carried chain claim and discharged transitively
    at the base. Verifier O(nb) instead of the line close's O(2^nb); a wrong
    val0/val1 propagates to a wrong residual and is caught at the base."""
    if pcs:
        res = _lo.open_batch(S, [(list(pt0), val0), (list(pt1), val1)],
                             rng, stats, q)
        return res["r"], res["residual"], res["ok"]
    t0 = time.perf_counter()
    hs = []
    for t in range(nb + 1):
        tf = t % q
        gp = [((q + 1 - tf) * int(a) + tf * int(b)) % q
              for a, b in zip(pt0, pt1)]
        hs.append(mle_eval(S, gp, q))
    stats["t_prover"] += time.perf_counter() - t0
    t0 = time.perf_counter()
    ok = hs[0] == val0 and hs[1] == val1
    tstar = int(rng.integers(0, q))
    new_claim = lagrange_eval(hs, tstar, q)
    new_pt = [((q + 1 - tstar) * int(a) + tstar * int(b)) % q
              for a, b in zip(pt0, pt1)]
    ok = ok and new_claim == mle_eval(S, new_pt, q)         # leaf close
    _acct_vleaf(stats, len(S))                              # per-layer leaf open
    stats["t_verifier"] += time.perf_counter() - t0
    stats["comm"] += nb + 1
    return new_pt, new_claim, ok


def batch_on_table(S, claims, nb, rng, stats, q=Q, pcs=False):
    """Fold k same-table point-claims into one by sequential line restriction.
    claims: list of (point, value). Returns (point, value, ok).

    pcs (S505): fold all k claims in ONE sum-check opening (leaf_open.open_batch)
    and RETURN the residual without closing it -- the caller carries it as the
    chain's small-side claim, discharged at the base. Replaces the (k-1) line
    closes' O(k*2^nb) verifier work with O(k*nb); the single-claim case becomes a
    plain open_eval (no anchoring close)."""
    if pcs:
        res = _lo.open_batch(S, [(list(pt), val) for (pt, val) in claims],
                             rng, stats, q)
        return res["r"], res["residual"], res["ok"]
    pt, val = claims[0]
    if len(claims) == 1:                                    # still anchor it
        t0 = time.perf_counter()
        ok = val == mle_eval(S, pt, q)
        _acct_vleaf(stats, len(S))                          # per-layer leaf open
        stats["t_verifier"] += time.perf_counter() - t0
        return pt, val, ok
    ok = True
    for (pt1, val1) in claims[1:]:
        pt, val, okp = line_batch_pair(S, pt, val, pt1, val1, nb, rng, stats, q)
        ok = ok and okp
    return pt, val, ok


# ----------------------------------------------------------------------
# one compressed layer, each side
# ----------------------------------------------------------------------

def pad(tab, nb, q=Q):
    out = np.zeros(1 << nb, dtype=_dt(q))
    a = (np.asarray(tab, dtype=np.int64) % q).astype(_dt(q))
    out[:len(a)] = a
    return out


def large_reduce(x, p, sp, ecut, nb, S_prev_small, S_prev_large, z, C,
                 rng, stats, cheat=None, delegate=False, structured=False,
                 q=Q, skip_trace=False, defer=None, pcs=False, ub_defer=None,
                 layer_idx=0):
    """Reduce S~_i^large(z)=C to: a large point-claim (new_z,new_claim) on
    S_{i-1}^large (line-batched s1@r_v + s2_aff@r_u) and a small point-claim
    (r_B,s_B) on S_{i-1}^small (trace region). cheat in
    {None,'G','split','u','waff','waff_forge','waff_chain'};
    'wrong_C'/'delta_pi' handled by callers via C.  delegate routes the affine
    wiring scalar through the O(nb log p) inner protocol; structured makes that
    inner protocol's chain prover O(nb*p) (S496 drop-in, delegate path only)."""
    V = isqrt(x)
    ci = sp % q
    lo, hi = ecut + 1, V - 1
    Bcut = min(V - 1, x // (p * p) - 1)
    Dml = 1 << nb

    # prover-committed G table over the e-cube
    t0 = time.perf_counter()
    G = np.zeros(Dml, dtype=_dt(q))
    for e in range(V):
        ep = (e + 1) * p
        G[e] = S_prev_large[ep - 1] if ep <= V else S_prev_small[x // ep]
    if cheat == "G":
        G[0] = (int(G[0]) + 1) % q
    Bf = (np.arange(Dml) <= Bcut).astype(_dt(q))
    stats["t_prover"] += time.perf_counter() - t0

    # ---- phase A ----
    termsA = [(1, ["EQ", "S"]), (q - 1, ["EQ", "B", "G"]), (ci, ["EQ", "B"])]
    okA, r_v, finalA, scal = sumcheck(C, {"EQ": eq_table(z, q),
                                          "S": S_prev_large.copy(),
                                          "B": Bf.copy(), "G": G.copy()},
                                      termsA, 3, rng, q)
    stats["comm"] += 4 * nb
    if not okA:
        return dict(accepted=False)
    s1, g1 = scal["S"], scal["G"]
    t0 = time.perf_counter()
    bv = le_eval(r_v, Bcut, nb, q)
    expect = eq_point(z, r_v, q) * ((s1 - bv * ((g1 - ci) % q)) % q) % q
    okA2 = (finalA % q) == expect
    stats["t_verifier"] += time.perf_counter() - t0
    if not okA2:
        return dict(accepted=False)

    # ---- honest split of g1 ----
    eqrv_full = eq_table(r_v, q)
    g1_aff = 0
    for e in range(ecut + 1):
        g1_aff = (g1_aff + int(eqrv_full[e]) * int(S_prev_large[(e + 1) * p - 1])) % q
    g1_trace = 0
    for e in range(lo, hi + 1):
        g1_trace = (g1_trace + int(eqrv_full[e]) * int(S_prev_small[x // ((e + 1) * p)])) % q
    if cheat == "split":
        g1_aff = (g1_aff + 1) % q
        g1_trace = (g1_trace - 1) % q
    t0 = time.perf_counter()
    ok_sum = (g1_aff + g1_trace) % q == g1 % q
    stats["t_verifier"] += time.perf_counter() - t0
    if not ok_sum:
        return dict(accepted=False)

    # ---- affine region ----
    okAff, r_u, s2_aff = verify_affine_region(
        S_prev_large, r_v, g1_aff, p, ecut, nb, rng, stats,
        cheat=(cheat if cheat in ("waff", "waff_forge", "waff_chain")
               else None),
        delegate=delegate, structured=structured, q=q, defer=defer, pcs=pcs)
    if not okAff:
        return dict(accepted=False)

    # ---- trace region ----
    W = build_witness(x, p, nb, dstart=1, q=q)
    if cheat == "u":                                        # wrong quotient in trace
        d = max(lo, 1) % W["D"]
        W["tabs"]["U"] = W["tabs"]["U"].copy()
        W["tabs"]["U"][d] = (int(W["tabs"]["U"][d]) + 1) % q
    okTr, r_B, s_B = verify_trace_region(
        W, S_prev_small, r_v, g1_trace, lo, hi, nb, rng, stats, q,
        skip_constraints=skip_trace, pcs=pcs, ub_defer=ub_defer,
        layer_idx=layer_idx)
    if not okTr:
        return dict(accepted=False)

    # ---- line-batch s1@r_v and s2_aff@r_u (both on S_prev_large) ----
    new_z, new_claim, okL = line_batch_pair(
        S_prev_large, r_v, s1, r_u, s2_aff, nb, rng, stats, q, pcs=pcs)
    if not okL:
        return dict(accepted=False)

    return dict(accepted=True, s_B=s_B, r_B=r_B,
                new_claim=new_claim, new_z=new_z, r_v=r_v)


def small_reduce(x, p, sp, nb, S_prev_small, z, C, rng, stats, cheat=None,
                 delegate=False, structured=False, q=Q, defer=None):
    """Reduce S~_i^small(z)=C to TWO point-claims on S_{i-1}^small:
    (r_v,s1) [the v itself] and (r_u,s2) [floor(v/p)]. The wiring v->floor(v/p)
    stays entirely small-side (value->value): the existing division automaton.
    Gate B(v)=[p^2<=v<=V] is a BAND (the upper bound keeps the recurrence off
    the padding). cheat in {None,'small_G','small_u','small_forge',
    'small_chain'}; 'small_C' via C.  delegate routes the division wiring
    scalar W~(r_v,r_u) through the SAME carry-chain (inner_verify_div);
    structured makes that chain prover O(nb*p) (S496 drop-in, delegate only).
    defer (S503, delegate only): when a list is passed, the division chain is NOT
    run inline; the obligation (p, r_v, r_u, accept_rem=None, claim=wv_claim,
    lie=small_chain) is appended for ONE batched discharge (run_chain
    batch_wiring=True)."""
    V = isqrt(x)
    ci = sp % q
    p2 = p * p
    Dml = 1 << nb
    arange = np.arange(Dml, dtype=np.int64)

    # ---- phase A: C = sum_v eq(z,v)[S(v) - B(v)(G(v)-ci)] ----
    t0 = time.perf_counter()
    Bf = ((arange >= p2) & (arange <= V)).astype(_dt(q))        # band [p^2<=v<=V]
    Gf = S_prev_small[arange // p].astype(_dt(q)).copy()        # G(v)=S(floor(v/p))
    if cheat == "small_G":
        gi = (p2) % Dml
        Gf[gi] = (int(Gf[gi]) + 1) % q
    stats["t_prover"] += time.perf_counter() - t0
    termsA = [(1, ["EQ", "S"]), (q - 1, ["EQ", "B", "G"]), (ci, ["EQ", "B"])]
    okA, r_v, finalA, scal = sumcheck(C, {"EQ": eq_table(z, q),
                                          "S": S_prev_small.copy(),
                                          "B": Bf.copy(), "G": Gf.copy()},
                                      termsA, 3, rng, q)
    stats["comm"] += 4 * nb
    if not okA:
        return dict(accepted=False)
    s1, g1 = scal["S"], scal["G"]
    t0 = time.perf_counter()
    bv = band_eval(r_v, p2, V, nb, q)
    expect = eq_point(z, r_v, q) * ((s1 - bv * ((g1 - ci) % q)) % q) % q
    okA2 = (finalA % q) == expect
    stats["t_verifier"] += time.perf_counter() - t0
    if not okA2:
        return dict(accepted=False)

    # ---- phase B: g1 = sum_u W(r_v,u) S(u),  W(v,u)=[u=floor(v/p)] ----
    t0 = time.perf_counter()
    eqrv = eq_table(r_v, q)
    if fast_big(q):                                         # < p summands overflow u64
        acc = np.zeros(Dml, dtype=object)
        np.add.at(acc, arange // p, eqrv.astype(object))
        Wt = (acc % q).astype(np.uint64)
    else:
        Wt = np.zeros(Dml, dtype=_dt(q))
        np.add.at(Wt, arange // p, eqrv)                    # < p summands < q
        Wt %= q
    if cheat in ("small_u", "small_forge"):
        Wt[0] = (int(Wt[0]) + 1) % q
    stats["t_prover"] += time.perf_counter() - t0
    okB, r_u, finalB, scalB = sumcheck(g1, {"W": Wt, "S": S_prev_small.copy()},
                                       [(1, ["W", "S"])], 2, rng, q)
    stats["comm"] += 3 * nb
    if not okB:
        return dict(accepted=False)
    s2 = scalB["S"]

    if not delegate:
        t0 = time.perf_counter()
        wv = w_div_eval(r_v, r_u, p, nb, q)                 # O(nb*p) automaton
        okB2 = (finalB % q) == wv * s2 % q
        stats["t_verifier"] += time.perf_counter() - t0
        if not okB2:
            return dict(accepted=False)
        return dict(accepted=True, claims=[(r_v, s1), (r_u, s2)])

    # delegate: prover supplies the wiring scalar (Õ(sqrt x), NO automaton),
    # verifier binds it to finalB then PROVES it via the carry chain.
    t0 = time.perf_counter()
    if cheat == "small_forge" and int(s2) % q != 0:
        wv_claim = finalB % q * pow(int(s2) % q, q - 2, q) % q  # forge to pass algebra
    else:
        eqru = eq_table(r_u, q)
        wv_claim = _asum(fmul(eqrv, eqru[arange // p], q), q)
    stats["t_prover"] += time.perf_counter() - t0
    t0 = time.perf_counter()
    okB2 = (finalB % q) == wv_claim * s2 % q
    stats["t_verifier"] += time.perf_counter() - t0
    if not okB2:
        return dict(accepted=False)
    if defer is not None:                       # batch_wiring: collect, discharge later
        defer.append((p, list(r_v), list(r_u), None, int(wv_claim) % q,
                      cheat == "small_chain"))
        return dict(accepted=True, claims=[(r_v, s1), (r_u, s2)])
    if not inner_verify_div(wv_claim, r_v, r_u, p, nb, rng, stats,
                            accept_rem=None, lie=(cheat == "small_chain"),
                            structured=structured, q=q):
        return dict(accepted=False)
    return dict(accepted=True, claims=[(r_v, s1), (r_u, s2)])


# ----------------------------------------------------------------------
# single-layer test drivers
# ----------------------------------------------------------------------

def _result(accepted, x, p, C, stats, nb, ecut, **extra):
    d = dict(accepted=accepted, x=x, p=p, C=C, nb=nb, ecut=ecut, **stats)
    d.update(extra)
    return d


def run_layer(x, layer, rng, z=None, cheat=None, delegate=False,
              structured=False, q=Q):
    """Verify ONE large-side compressed layer (step-2 single-layer test).
    cheat in {None, wrong_C, G, split, u, waff, waff_forge, waff_chain}."""
    V = isqrt(x)
    nb = max(1, V.bit_length())
    primes, sm_layers, lg_layers = compressed_lucy(x)
    assert 1 <= layer <= len(primes)
    p = primes[layer - 1]
    sp = sm_layers[layer - 1][p - 1]
    ecut = V // p - 1
    S_prev_small = pad(sm_layers[layer - 1], nb, q)
    S_prev_large = pad(lg_layers[layer - 1], nb, q)
    S_i_large = pad(lg_layers[layer], nb, q)
    if z is None:
        z = [int(rng.integers(0, q)) for _ in range(nb)]
    C = mle_eval(S_i_large, z, q)
    if cheat == "wrong_C":
        C = (C + 1) % q
    stats = {"t_prover": 0.0, "t_verifier": 0.0, "comm": 0}
    res = large_reduce(x, p, sp, ecut, nb, S_prev_small, S_prev_large, z, C,
                       rng, stats,
                       cheat=(cheat if cheat in ("G", "split", "u", "waff",
                                                 "waff_forge", "waff_chain")
                              else None),
                       delegate=delegate, structured=structured, q=q)
    out = _result(res["accepted"], x, p, C, stats, nb, ecut)
    if res["accepted"]:
        out.update(s_B=res["s_B"], r_B=res["r_B"], new_claim=res["new_claim"],
                   new_z=res["new_z"], S_prev_small=S_prev_small,
                   S_prev_large=S_prev_large)
    return out


def run_two_sided_layer(x, layer, rng, cheat=None, delegate=False,
                        structured=False, q=Q, batch_wiring=False):
    """Verify ONE full two-sided layer: (large-claim, small-claim) ->
    (large-claim', small-claim'). Small-side claims about S_{i-1}^small (s_B
    from the large trace + the two from the small layer) are folded to one.

    batch_wiring (S503, delegate only): both delegated wiring chains (large affine
    image + small division) are DEFERRED and discharged in ONE
    batched_wiring.verify_wiring_obligations -- the single-layer analogue of
    run_chain(batch_wiring=True), used to confirm the deferred/batched integration
    still rejects the wiring-specific cheats (small_forge/small_chain/waff_chain)."""
    if batch_wiring and not delegate:
        raise ValueError("batch_wiring requires delegate=True")
    V = isqrt(x)
    nb = max(1, V.bit_length())
    primes, sm_layers, lg_layers = compressed_lucy(x)
    assert 1 <= layer <= len(primes)
    p = primes[layer - 1]
    sp = sm_layers[layer - 1][p - 1]
    ecut = V // p - 1
    S_prev_small = pad(sm_layers[layer - 1], nb, q)
    S_prev_large = pad(lg_layers[layer - 1], nb, q)
    S_i_large = pad(lg_layers[layer], nb, q)
    S_i_small = pad(sm_layers[layer], nb, q)
    stats = {"t_prover": 0.0, "t_verifier": 0.0, "comm": 0}

    zL = [int(rng.integers(0, q)) for _ in range(nb)]
    zS = [int(rng.integers(0, q)) for _ in range(nb)]
    CL = mle_eval(S_i_large, zL, q)
    CS = mle_eval(S_i_small, zS, q)
    if cheat == "wrong_C":
        CL = (CL + 1) % q
    if cheat == "small_C":
        CS = (CS + 1) % q

    defer = [] if batch_wiring else None
    lres = large_reduce(x, p, sp, ecut, nb, S_prev_small, S_prev_large, zL, CL,
                        rng, stats,
                        cheat=(cheat if cheat in ("G", "split", "u", "waff",
                                                  "waff_forge", "waff_chain")
                               else None),
                        delegate=delegate, structured=structured, q=q,
                        defer=defer)
    if not lres["accepted"]:
        return dict(accepted=False, x=x, p=p, nb=nb, ecut=ecut, **stats)
    sres = small_reduce(x, p, sp, nb, S_prev_small, zS, CS, rng, stats,
                        cheat=(cheat if cheat in ("small_G", "small_u",
                                                  "small_forge", "small_chain")
                               else None),
                        delegate=delegate, structured=structured, q=q,
                        defer=defer)
    if not sres["accepted"]:
        return dict(accepted=False, x=x, p=p, nb=nb, ecut=ecut, **stats)

    sB = lres["s_B"]
    if cheat == "carry":                                    # tamper cross-side carry
        sB = (sB + 1) % q
    small_claims = [(lres["r_B"], sB)] + sres["claims"]
    zS2, CS2, okb = batch_on_table(S_prev_small, small_claims, nb, rng, stats, q)

    out_ok = (okb
              and CS2 == mle_eval(S_prev_small, zS2, q)
              and lres["new_claim"] == mle_eval(S_prev_large, lres["new_z"], q))
    if batch_wiring and defer:                  # discharge both deferred wirings at once
        import batched_wiring as _bw
        l_max = max(1, (p - 1).bit_length())
        out_ok = out_ok and _bw.verify_wiring_obligations(
            defer, nb, l_max, rng, stats, q)
    return dict(accepted=out_ok, x=x, p=p, nb=nb, ecut=ecut,
                CS2=CS2, zS2=zS2, new_claim=lres["new_claim"],
                new_z=lres["new_z"], S_prev_small=S_prev_small,
                S_prev_large=S_prev_large, **stats)


# ----------------------------------------------------------------------
# the full compressed chain: verify pi(x) end-to-end
# ----------------------------------------------------------------------

def run_chain(x, rng, cheat=None, corrupt_layer=None, delegate=False,
              structured=False, q=Q, batch_trace=False, batch_wiring=False,
              pcs=False, batch_ub=False, commit_base=False, commit_code="rs"):
    """Chain all K layers from S_K^large(e=0)=pi(x) down to S_0 (opened on
    both sides). cheat='delta_pi' lies about pi(x); corrupt_layer=i0 is a
    self-consistent liar. delegate routes BOTH wirings (small division,
    large affine image) through the O(nb log p) inner protocol -- the whole
    per-layer verifier is then O(nb log p), no p-linear automaton term.
    structured (delegate only) makes the inner chain PROVER O(nb*p) per wiring
    (S496 drop-in) -- whole-chain wiring prover Õ(x^{3/2})->Õ(x), transcript
    bit-identical.  q selects the prime field: Q=2^31-1 (demo, uint64) is sound
    only for x < q (n <~ 30); BIG_Q=2^61-1 (object dtype) lifts the
    CHARACTERISTIC so the trace identities stay sound to n <~ 60 -- the
    wrap-around alias the demo prime admits above its field is rejected (S498).

    batch_trace (S502): the K per-layer trace zero-tests (verify_constraints,
    53% of wall by the S501 profile -- K independent nb-cube sum-checks) are
    replaced by ONE batched sum-check over the stacked (K*2^nb) cube
    (batched_trace.verify_constraints_batched), run BEFORE the layer loop, then
    skipped inside each large_reduce.  fmul count drops ~K-fold / width rises
    ~K-fold (the precondition for the fast-Mersenne path to win).  The chain's
    trace witnesses are deterministic in (x, p_l) and independent of the carried
    chain claim, so this is a faithful detach -- the chain's soundness vs
    delta_pi / self-consistent liar lives in phase-A and the base check, not the
    trace test, so accept/claimed are unchanged.

    batch_wiring (S503, delegate only): the K delegated per-layer wiring chains
    (76% of all sum-check CALLS by the S501 profile -- the small division and large
    affine-image inner_verify_div, each on a TINY 2^l cube) are DEFERRED via the
    `defer` accumulator and discharged in ONE batched chain over the stacked
    (2^Lk * 2^l_max) cube (batched_wiring.verify_wiring_obligations) AFTER the layer
    loop.  fmul CALL count drops ~K-fold / per-fmul WIDTH rises ~K-fold -- with
    batch_trace this is the precondition for the fast-Mersenne path (FAST_BIG) to be
    a net end-to-end win.  The per-layer comparator fold of verify_waff_value
    (O(2^nb), p-independent) stays per-layer.  The obligations are deterministic in
    (x, p_l) and the per-layer outer sum-checks pin each claim BEFORE deferral, so
    the detach is faithful: the chain's soundness vs delta_pi / liar lives in
    phase-A and the base check; the wiring batch additionally rejects a forged
    wiring scalar or a self-consistent lying inner chain (small_forge/small_chain/
    waff_chain), exactly as the K inline inner_verify_div did.

    pcs (S505): replace the per-layer carried-claim LEAF OPENINGS -- the
    line_batch_pair / batch_on_table folds and the redundant s2/s_B closes in
    verify_affine_region/verify_trace_region -- with real sum-check OPENINGS
    (leaf_open.open_batch).  Each fold becomes a degree-2 sum-check whose RESIDUAL
    is carried as the next layer's claim and discharged transitively at the base
    (the two S_0 mle_eval closes, a one-time O(sqrt x)).  This removes the
    line-restriction folding's O(k*nb*2^nb) per-layer verifier work; the only
    remaining per-layer O(2^nb) term is verify_trace_region's nb Ub-bit-table
    openings (line 429 -- still mle_eval; their residuals are per-layer witness
    data with no carried claim to thread to, so they want the batched-trace
    integration, a documented follow-on).  UNCONDITIONAL: no commitment -- the
    sum-check soundness rides the existing GKR layer reductions to the closed base.
    delta_pi / self-consistent liar are caught in phase-A round-1 (untouched by the
    leaf openings), so the verdict is preserved.  Default pcs=False keeps every
    prior artifact verbatim (the transcript is NOT bit-identical under pcs).

    batch_ub (S506, requires pcs AND batch_trace): the ONE per-layer O(2^nb) verifier
    term pcs leaves standing -- verify_trace_region's nb Ub-bit-table openings (the B2
    routing's Ub_{Lv-nb+k}~(r_C), per-layer WITNESS data with no carried claim to
    thread to the base) -- is DEFERRED via the `ub_defer` accumulator (the prover
    supplies the scalars; the verifier folds them into its B2 `expect` in O(nb)) and
    discharged in ONE degree-2 sum-check (batched_trace.verify_ub_openings_batched)
    against the SAME stacked Ub cube batch_trace's zero-test commits.  The verifier-
    side O(2^nb) Ub-leaf eval count drops from K*nb to 0 (moved to the prover); the
    per-layer verifier is then leaf-eval-free, so the only O(sqrt x) verifier work
    left is one-time (the two S_0 base closes + the batched discharges) -- the chain
    verifier is honestly Õ(sqrt x) end-to-end.  Verdict unchanged (the same `expect`
    check runs; a forged Ub opening is caught by the batched discharge, the chain
    cheats by phase-A/base as before).

    commit_base (S508): discharge the TWO one-time S_0 base closes -- the last
    verifier work whose cost grows with x (each a direct O(2^nb)=O(sqrt x) mle_eval)
    -- through the tensor / Ligero-Brakedown polynomial commitment (pcs_commit:
    commit/prove/verify).  Its evaluation-proof VERIFIER is SUB-sqrt x
    (O(t*(r+k)) ~ O(x^{1/4}) field-ops, tallied in stats['vcommit_ops']); the full
    2^nb tables stay on the prover side.  With pcs+batch_trace+batch_ub this makes
    the WHOLE-CHAIN verifier sub-sqrt x (polylog above sqrt x) -- the milestone's
    full succinctness -- now under a collision-resistant-hash / random-oracle
    assumption (vs the otherwise-unconditional Õ(sqrt x) verifier).  Verdict
    unchanged: the commitment opens to the SAME S~(z) value the mle_eval close
    checked, so honest accepts and delta_pi / self-consistent liar (caught in
    phase-A / by the carried-claim value reaching the base) are rejected as before.
    commit_code (S514) selects the pcs row code: "rs" (Reed-Solomon, needs the
    codeword length N<=q -- the demo constraint that forced large runs onto BIG_Q)
    or "expander" (Brakedown-style linear-time code, field-size-free, no N<=q
    requirement; same sub-sqrt-x verifier).  Verdict identical either way.
    Returns dict(accepted, claimed, layers, ...)."""
    if batch_wiring and not delegate:
        raise ValueError("batch_wiring requires delegate=True")
    if batch_ub and not (pcs and batch_trace):
        # the Ub openings ride the pcs deferral path and discharge against the
        # SAME stacked Ub cube the trace zero-test (batch_trace) commits.
        raise ValueError("batch_ub requires pcs=True and batch_trace=True")
    V = isqrt(x)
    nb = max(1, V.bit_length())
    primes, sm_layers, lg_layers = compressed_lucy(x)
    K = len(primes)
    if corrupt_layer is not None:
        sm_layers, lg_layers = corrupt_compressed(primes, sm_layers,
                                                  lg_layers, x, corrupt_layer)
    claimed_pi = lg_layers[K][0]
    if cheat == "delta_pi":
        claimed_pi += 1

    stats = {"t_prover": 0.0, "t_verifier": 0.0, "comm": 0,
             "vleaf_ops_pl": 0, "vleaf_cnt_pl": 0,        # per-layer leaf-eval ops
             "vleaf_ops_ot": 0, "vleaf_cnt_ot": 0,        # one-time leaf-eval ops
             "vcommit_ops": 0,                            # commit-verifier field-ops
             # S509: certificate-SIZE attribution -- partition stats["comm"] by
             # source so the Õ(sqrt x) total can be localized.  Purely additive
             # (boundary snapshots of stats["comm"]); the transcript is unchanged.
             "comm_outer": 0,   # K SEQUENTIAL per-layer two-sided reductions
             "comm_bt": 0,      # one-time batched trace zero-test
             "comm_bw": 0,      # one-time batched wiring discharge
             "comm_ub": 0,      # one-time batched Ub-opening discharge
             "comm_base": 0}    # one-time S_0 base close / commit proof
    z_large = [0] * nb                                       # Boolean point e=0
    C_large = claimed_pi % q
    z_small, C_small = None, None

    def fail():
        return dict(accepted=False, claimed=claimed_pi, layers=K, **stats)

    if batch_trace:                       # one batched zero-test for all K layers
        import batched_trace as _bt
        Ws = [build_witness(x, p, nb, dstart=1, q=q) for p in primes]
        _c0 = stats["comm"]
        if not _bt.verify_constraints_batched(Ws, primes, x, nb, rng, stats, q):
            return fail()
        stats["comm_bt"] += stats["comm"] - _c0

    defer = [] if batch_wiring else None   # collect K wiring obligations; discharge once
    ub_defer = [] if batch_ub else None    # collect K*nb Ub openings; discharge once

    for i in range(K, 0, -1):
        p = primes[i - 1]
        sp = sm_layers[i - 1][p - 1]
        ecut = V // p - 1
        S_prev_small = pad(sm_layers[i - 1], nb, q)
        S_prev_large = pad(lg_layers[i - 1], nb, q)
        _c0 = stats["comm"]                 # S509: attribute this layer's outer comm

        lres = large_reduce(x, p, sp, ecut, nb, S_prev_small, S_prev_large,
                            z_large, C_large, rng, stats, delegate=delegate,
                            structured=structured, q=q, skip_trace=batch_trace,
                            defer=defer, pcs=pcs, ub_defer=ub_defer,
                            layer_idx=i - 1)
        if not lres["accepted"]:
            return fail()
        small_claims = [(lres["r_B"], lres["s_B"])]
        if C_small is not None:                             # skipped at i=K
            sres = small_reduce(x, p, sp, nb, S_prev_small, z_small, C_small,
                                rng, stats, delegate=delegate,
                                structured=structured, q=q, defer=defer)
            if not sres["accepted"]:
                return fail()
            small_claims += sres["claims"]
        z_small, C_small, okb = batch_on_table(S_prev_small, small_claims,
                                               nb, rng, stats, q, pcs=pcs)
        if not okb:
            return fail()
        z_large, C_large = lres["new_z"], lres["new_claim"]
        stats["comm_outer"] += stats["comm"] - _c0   # this layer's outer comm

    if batch_wiring and defer:             # ONE batched chain for all 2K-1 wirings
        import batched_wiring as _bw
        l_max = max(max(1, (pp - 1).bit_length()) for pp in primes)
        _c0 = stats["comm"]
        if not _bw.verify_wiring_obligations(defer, nb, l_max, rng, stats, q):
            return fail()
        stats["comm_bw"] += stats["comm"] - _c0

    if batch_ub and ub_defer:              # ONE batched opening for all K*nb Ub leaves
        import batched_trace as _bt        # (Ws built above under batch_trace)
        _c0 = stats["comm"]
        if not _bt.verify_ub_openings_batched(Ws, nb, ub_defer, rng, stats, q):
            return fail()
        stats["comm_ub"] += stats["comm"] - _c0

    # base: open S_0 on both sides.  Default: a direct mle_eval close (the
    # leaf-opening stand-in, a one-time O(sqrt x) = O(2^nb)).  commit_base (S508):
    # discharge each close through the tensor polynomial commitment (pcs_commit),
    # whose evaluation-proof VERIFIER is sub-sqrt x (O(t*(r+k)) ~ O(x^{1/4})), so
    # the last one-time x-scaling verifier term drops below sqrt x -- making the
    # whole-chain verifier polylog under the hash (RO/CRH) assumption.  The full
    # 2^nb tables never touch the verifier; the commit is the prover's one-time
    # O(x^{3/4}) work (within its Õ(x) budget).
    t0 = time.perf_counter()
    _c0 = stats["comm"]                                     # S509: base-close comm
    S0_large = pad(lg_layers[0], nb, q)
    S0_small = pad(sm_layers[0], nb, q)
    if commit_base:
        import pcs_commit as _pcs
        ok = True
        for S0, z, C in ((S0_large, z_large, C_large),
                         (S0_small, z_small, C_small)):
            com = _pcs.commit(S0, q, code=commit_code)
            pf = _pcs.prove(com, z, C)
            okp, vop = _pcs.verify(com["root"], z, C, pf, q=q)
            ok = ok and okp
            _acct_vleaf(stats, vop, onetime=True)           # sub-sqrt-x base open
            stats["vcommit_ops"] += vop
            stats["comm"] += (len(pf["v"]) + len(pf["w"])  # proof: messages + cols
                              + sum(len(cv) for (_, cv, _) in pf["opened"]))
    else:
        ok = (C_large == mle_eval(S0_large, z_large, q)
              and C_small == mle_eval(S0_small, z_small, q))
        _acct_vleaf(stats, len(S0_large), onetime=True)     # one-time base close
        _acct_vleaf(stats, len(S0_small), onetime=True)     # one-time base close
    stats["t_verifier"] += time.perf_counter() - t0
    stats["comm_base"] += stats["comm"] - _c0
    return dict(accepted=ok, claimed=claimed_pi, layers=K, **stats)


# ----------------------------------------------------------------------
# tests / experiments
# ----------------------------------------------------------------------

def sieve_pi(x):
    s = np.ones(x + 1, dtype=bool)
    s[:2] = False
    for q in range(2, isqrt(x) + 1):
        if s[q]:
            s[q * q::q] = False
    return int(s.sum())


def selftest():
    rng = np.random.default_rng(0)
    # 1. compressed DP == sieve
    for n in range(6, 17):
        x = (1 << n) - 1
        _, _, lg = compressed_lucy(x)
        assert lg[-1][0] == sieve_pi(x), n

    # 2. le / band MLEs vs Boolean truth (exhaustive small n)
    for n in (3, 4, 5):
        for M in range(-1, (1 << n) + 1):
            for e in range(1 << n):
                be = [(e >> (n - 1 - j)) & 1 for j in range(n)]
                assert le_eval(be, M, n) == (1 if e <= M else 0), (n, M, e)
        for lo in range(0, 1 << n):
            for hi in range(lo, 1 << n):
                for e in range(1 << n):
                    be = [(e >> (n - 1 - j)) & 1 for j in range(n)]
                    assert band_eval(be, lo, hi, n) == (1 if lo <= e <= hi else 0)
        # empty band (lo>hi, e.g. small-side gate when p^2>V) -> all-zero MLE
        for lo in range(1, 1 << n):
            r = [int(rng.integers(0, Q)) for _ in range(n)]
            assert band_eval(r, lo, lo - 1, n) == 0

    # 3. waff_eval vs Boolean truth and vs full-table MLE at random points
    for (n, p) in [(4, 3), (5, 5), (5, 7), (6, 2)]:
        for ecut in (0, 1, (1 << n) // p):
            full = np.zeros(1 << (2 * n), dtype=np.uint64)
            for e in range(1 << n):
                ep1 = p * e + (p - 1)
                if e <= ecut and ep1 < (1 << n):
                    full[(e << n) | ep1] = 1
            for e in range(1 << n):
                for ep in range(1 << n):
                    be = [(e >> (n - 1 - j)) & 1 for j in range(n)]
                    bp = [(ep >> (n - 1 - j)) & 1 for j in range(n)]
                    want = 1 if (e <= ecut and ep == p * e + (p - 1)) else 0
                    assert waff_eval(be, bp, p, ecut, n) == want, (n, p, ecut, e, ep)
            for _ in range(5):
                rv = [int(rng.integers(0, Q)) for _ in range(n)]
                ru = [int(rng.integers(0, Q)) for _ in range(n)]
                assert waff_eval(rv, ru, p, ecut, n) == mle_eval(full, rv + ru)

    # 4. honest LARGE layer accepts; output claims equal true MLEs of S_{i-1}
    for (n, layer) in [(10, 2), (10, 4), (12, 3), (12, 5)]:
        x = (1 << n) - 1
        res = run_layer(x, layer, np.random.default_rng(7))
        assert res["accepted"], (n, layer)
        assert res["s_B"] == mle_eval(res["S_prev_small"], res["r_B"])
        assert res["new_claim"] == mle_eval(res["S_prev_large"], res["new_z"])

    # 5. honest g1 == g1_aff + g1_trace decomposition (phase-A / phase-B tie)
    x = (1 << 12) - 1
    V = isqrt(x); nb = V.bit_length()
    primes, sm, lg = compressed_lucy(x)
    p = primes[3]; ecut = V // p - 1; lo, hi = ecut + 1, V - 1
    S_prev_small = pad(sm[3], nb); S_prev_large = pad(lg[3], nb)
    r_v = [int(rng.integers(0, Q)) for _ in range(nb)]
    eqrv = eq_table(r_v)
    G = np.zeros(1 << nb, dtype=np.uint64)
    for e in range(V):
        ep = (e + 1) * p
        G[e] = S_prev_large[ep - 1] if ep <= V else S_prev_small[x // ep]
    g1 = mle_eval(G, r_v)
    g1a = sum(int(eqrv[e]) * int(S_prev_large[(e + 1) * p - 1])
              for e in range(ecut + 1)) % Q
    g1t = sum(int(eqrv[e]) * int(S_prev_small[x // ((e + 1) * p)])
              for e in range(lo, hi + 1)) % Q
    assert (g1a + g1t) % Q == g1 % Q

    # 6. LARGE-side soundness: every cheat class rejected across primes/layers
    for (n, layer) in [(10, 3), (12, 4), (12, 6)]:
        x = (1 << n) - 1
        for ch in ["wrong_C", "G", "split", "u", "waff"]:
            rej = sum(not run_layer(x, layer, np.random.default_rng(50 + t),
                                    cheat=ch)["accepted"] for t in range(8))
            assert rej == 8, (n, layer, ch, rej)

    # 7. honest SMALL layer accepts; both output claims equal true MLEs.
    #    (also exercises padding via odd n: V != 2^nb-1)
    for (n, layer) in [(10, 2), (11, 3), (12, 4), (12, 6), (13, 5)]:
        x = (1 << n) - 1
        Vv = isqrt(x); nbv = max(1, Vv.bit_length())
        primes, sm, lg = compressed_lucy(x)
        p = primes[layer - 1]
        sp = sm[layer - 1][p - 1]
        Spsmall = pad(sm[layer - 1], nbv)
        Sismall = pad(sm[layer], nbv)
        zS = [int(rng.integers(0, Q)) for _ in range(nbv)]
        CS = mle_eval(Sismall, zS)
        st = {"t_prover": 0.0, "t_verifier": 0.0, "comm": 0}
        r = small_reduce(x, p, sp, nbv, Spsmall, zS, CS,
                         np.random.default_rng(11 + layer), st)
        assert r["accepted"], (n, layer)
        for (pt, val) in r["claims"]:
            assert val == mle_eval(Spsmall, pt), (n, layer)

    # 8. SMALL-side soundness: small_C / small_G / small_u rejected
    for (n, layer) in [(10, 3), (12, 5), (13, 4)]:
        x = (1 << n) - 1
        Vv = isqrt(x); nbv = max(1, Vv.bit_length())
        primes, sm, lg = compressed_lucy(x)
        p = primes[layer - 1]; sp = sm[layer - 1][p - 1]
        Spsmall = pad(sm[layer - 1], nbv); Sismall = pad(sm[layer], nbv)
        for ch in ["small_C", "small_G", "small_u"]:
            rej = 0
            for t in range(8):
                rr = np.random.default_rng(300 + t)
                zS = [int(rr.integers(0, Q)) for _ in range(nbv)]
                CS = mle_eval(Sismall, zS)
                if ch == "small_C":
                    CS = (CS + 1) % Q
                st = {"t_prover": 0.0, "t_verifier": 0.0, "comm": 0}
                res = small_reduce(x, p, sp, nbv, Spsmall, zS, CS, rr, st,
                                   cheat=(ch if ch != "small_C" else None))
                rej += (not res["accepted"])
            assert rej == 8, (n, layer, ch, rej)

    # 9. honest TWO-SIDED single layer accepts; outputs == true MLEs
    for (n, layer) in [(10, 2), (11, 4), (12, 3), (12, 7)]:
        x = (1 << n) - 1
        res = run_two_sided_layer(x, layer, np.random.default_rng(9 + layer))
        assert res["accepted"], (n, layer)
        assert res["CS2"] == mle_eval(res["S_prev_small"], res["zS2"])
        assert res["new_claim"] == mle_eval(res["S_prev_large"], res["new_z"])

    # 10. TWO-SIDED soundness: every cheat (both sides + carry) rejected
    for (n, layer) in [(10, 3), (12, 5)]:
        x = (1 << n) - 1
        for ch in ["wrong_C", "G", "split", "u", "waff",
                   "small_C", "small_G", "small_u", "carry"]:
            rej = sum(not run_two_sided_layer(x, layer,
                                              np.random.default_rng(400 + t),
                                              cheat=ch)["accepted"]
                      for t in range(6))
            assert rej == 6, (n, layer, ch, rej)

    # 11. honest CHAIN accepts and matches the sieve (incl. odd n / padding)
    for n in (8, 10, 11, 12):
        x = (1 << n) - 1
        res = run_chain(x, np.random.default_rng(5))
        assert res["accepted"] and res["claimed"] == sieve_pi(x), (n, res)

    # 12. CHAIN soundness: delta_pi and self-consistent liar rejected
    for n in (10, 12):
        x = (1 << n) - 1
        K = len(compressed_lucy(x)[0])
        rej = sum(not run_chain(x, np.random.default_rng(600 + t),
                                cheat="delta_pi")["accepted"] for t in range(5))
        assert rej == 5, (n, "delta_pi", rej)
        for i0 in sorted(set([K, max(1, K // 2), 1])):
            rej = sum(not run_chain(x, np.random.default_rng(700 + i0 + t),
                                    corrupt_layer=i0)["accepted"]
                      for t in range(5))
            assert rej == 5, (n, "corrupt", i0, rej)

    # 13. DELEGATED wiring scalars equal the O(nb*p) automaton values, and the
    #     inner protocol accepts truth / rejects truth+1 / rejects a lying
    #     chain -- for BOTH the large affine map (verify_waff_value) and the
    #     small division (inner_verify_div, accept_rem=None == w_div_eval).
    for (n, layer) in [(10, 3), (12, 4), (12, 6)]:
        x = (1 << n) - 1
        Vv = isqrt(x); nbv = max(1, Vv.bit_length())
        primes, _, _ = compressed_lucy(x)
        p = primes[layer - 1]; ecut = Vv // p - 1
        for _ in range(4):
            r_v = [int(rng.integers(0, Q)) for _ in range(nbv)]
            r_u = [int(rng.integers(0, Q)) for _ in range(nbv)]
            st = {"t_prover": 0.0, "t_verifier": 0.0, "comm": 0}
            waffv = waff_eval(r_v, r_u, p, ecut, nbv)        # automaton truth
            assert verify_waff_value(waffv, r_v, r_u, p, ecut, nbv,
                                     np.random.default_rng(1), st)
            assert not verify_waff_value((waffv + 1) % Q, r_v, r_u, p, ecut,
                                         nbv, np.random.default_rng(2), st)
            assert not verify_waff_value(waffv, r_v, r_u, p, ecut, nbv,
                                         np.random.default_rng(3), st,
                                         lie_chain=True)
            divv = w_div_eval(r_v, r_u, p, nbv)              # automaton truth
            assert inner_verify_div(divv, r_v, r_u, p, nbv,
                                    np.random.default_rng(4), st)
            assert not inner_verify_div((divv + 1) % Q, r_v, r_u, p, nbv,
                                        np.random.default_rng(5), st)
            assert not inner_verify_div(divv, r_v, r_u, p, nbv,
                                        np.random.default_rng(6), st, lie=True)

    # 14. DELEGATED end-to-end: honest layer/chain accept & match the sieve;
    #     every cheat (incl. forge / bad-chain on both wirings) rejected. The
    #     verifier never runs an O(nb*p) automaton -- wiring is the carry chain.
    for (n, layer) in [(10, 2), (12, 5)]:
        x = (1 << n) - 1
        res = run_two_sided_layer(x, layer, np.random.default_rng(9 + layer),
                                  delegate=True)
        assert res["accepted"], (n, layer, "deleg")
        assert res["CS2"] == mle_eval(res["S_prev_small"], res["zS2"])
        assert res["new_claim"] == mle_eval(res["S_prev_large"], res["new_z"])
    for (n, layer) in [(10, 3), (12, 5)]:
        x = (1 << n) - 1
        for ch in ["wrong_C", "G", "split", "u", "waff", "waff_forge",
                   "waff_chain", "small_C", "small_G", "small_u",
                   "small_forge", "small_chain", "carry"]:
            rej = sum(not run_two_sided_layer(
                x, layer, np.random.default_rng(800 + t), cheat=ch,
                delegate=True)["accepted"] for t in range(6))
            assert rej == 6, (n, layer, ch, rej, "delegate")
    for n in (8, 10, 11, 12):
        x = (1 << n) - 1
        res = run_chain(x, np.random.default_rng(5), delegate=True)
        assert res["accepted"] and res["claimed"] == sieve_pi(x), (n, "deleg")
    for n in (10, 12):
        x = (1 << n) - 1
        K = len(compressed_lucy(x)[0])
        rej = sum(not run_chain(x, np.random.default_rng(600 + t),
                                cheat="delta_pi", delegate=True)["accepted"]
                  for t in range(5))
        assert rej == 5, (n, "delta_pi", "deleg", rej)
        for i0 in sorted(set([K, max(1, K // 2), 1])):
            rej = sum(not run_chain(x, np.random.default_rng(700 + i0 + t),
                                    corrupt_layer=i0, delegate=True)["accepted"]
                      for t in range(5))
            assert rej == 5, (n, "corrupt", "deleg", i0, rej)

    # 15. STRUCTURED (S496) is a clean drop-in for the DELEGATED compressed
    #     chain: the inner wiring sweep runs in O(nb*p) prover with a BIT-
    #     IDENTICAL transcript, so accept/reject, comm and claimed pi are
    #     UNCHANGED vs the dense O(nb*p^2) prover -- honest AND every cheat.
    #     This is what makes the whole-chain wiring prover Õ(x^{3/2})->Õ(x)
    #     a safe, verifier-untouched swap (mirrors lucy_dp selftest 6/6b).
    #     (a) the two delegated primitives, structured on vs off, agree on
    #         accept/reject AND comm (affine map + small division).
    for (n, layer) in [(10, 3), (12, 4)]:
        x = (1 << n) - 1
        Vv = isqrt(x); nbv = max(1, Vv.bit_length())
        primes, _, _ = compressed_lucy(x)
        p = primes[layer - 1]; ecut = Vv // p - 1
        for _ in range(3):
            r_v = [int(rng.integers(0, Q)) for _ in range(nbv)]
            r_u = [int(rng.integers(0, Q)) for _ in range(nbv)]
            waffv = waff_eval(r_v, r_u, p, ecut, nbv)
            divv = w_div_eval(r_v, r_u, p, nbv)
            for (claim, lie, want) in [(waffv, False, True),
                                       ((waffv + 1) % Q, False, False),
                                       (waffv, True, False)]:
                sd = {"t_prover": 0.0, "t_verifier": 0.0, "comm": 0}
                ss = {"t_prover": 0.0, "t_verifier": 0.0, "comm": 0}
                od = verify_waff_value(claim, r_v, r_u, p, ecut, nbv,
                                       np.random.default_rng(20), sd,
                                       lie_chain=lie)
                os_ = verify_waff_value(claim, r_v, r_u, p, ecut, nbv,
                                        np.random.default_rng(20), ss,
                                        lie_chain=lie, structured=True)
                assert od == os_ == want and sd["comm"] == ss["comm"]
            for (claim, lie, want) in [(divv, False, True),
                                       ((divv + 1) % Q, False, False),
                                       (divv, True, False)]:
                sd = {"t_prover": 0.0, "t_verifier": 0.0, "comm": 0}
                ss = {"t_prover": 0.0, "t_verifier": 0.0, "comm": 0}
                od = inner_verify_div(claim, r_v, r_u, p, nbv,
                                      np.random.default_rng(21), sd, lie=lie)
                os_ = inner_verify_div(claim, r_v, r_u, p, nbv,
                                       np.random.default_rng(21), ss, lie=lie,
                                       structured=True)
                assert od == os_ == want and sd["comm"] == ss["comm"]
    #     (b) full DELEGATED two-sided layer: structured on vs off match
    #         accept/reject AND comm on honest + every cheat class.
    for (n, layer) in [(10, 2), (12, 5)]:
        x = (1 << n) - 1
        for ch in [None, "wrong_C", "G", "waff", "waff_forge", "waff_chain",
                   "small_C", "small_u", "small_forge", "small_chain", "carry"]:
            d = run_two_sided_layer(x, layer, np.random.default_rng(33),
                                    cheat=ch, delegate=True)
            s = run_two_sided_layer(x, layer, np.random.default_rng(33),
                                    cheat=ch, delegate=True, structured=True)
            assert d["accepted"] == s["accepted"] and d["comm"] == s["comm"], \
                (n, layer, ch)
    #     (c) full DELEGATED chain: structured on vs off match accept/reject,
    #         comm AND claimed pi -- honest, delta_pi, self-consistent liar.
    for n in (8, 10, 12):
        x = (1 << n) - 1
        K = len(compressed_lucy(x)[0])
        for kw in [dict(), dict(cheat="delta_pi"),
                   dict(corrupt_layer=max(1, K // 2))]:
            d = run_chain(x, np.random.default_rng(44), delegate=True, **kw)
            s = run_chain(x, np.random.default_rng(44), delegate=True,
                          structured=True, **kw)
            assert (d["accepted"] == s["accepted"] and d["comm"] == s["comm"]
                    and d["claimed"] == s["claimed"]), (n, kw)

    # 16. THE FIELD LIFT THROUGH THE CHAIN (S499): run the WHOLE compressed
    #     chain over BIG_Q = 2^61-1 (characteristic > x, object-dtype exact
    #     arithmetic), proving the q-threading is correct end-to-end -- claimed
    #     pi == sieve in all three modes (automaton / delegated-dense /
    #     delegated-structured) and every chain cheat rejected.  The DEMO prime
    #     Q = 2^31-1 is sound only for x < Q (n <~ 30); BIG_Q lifts it.
    for n in (8, 10):
        x = (1 << n) - 1
        truth = sieve_pi(x)
        for (deleg, struct) in [(False, False), (True, False), (True, True)]:
            r = run_chain(x, np.random.default_rng(5), delegate=deleg,
                          structured=struct, q=BIG_Q)
            assert r["accepted"] and r["claimed"] == truth, (n, deleg, struct)
        K = len(compressed_lucy(x)[0])
        assert all(not run_chain(x, np.random.default_rng(600 + t),
                                 cheat="delta_pi", q=BIG_Q)["accepted"]
                   for t in range(3))
        assert all(not run_chain(x, np.random.default_rng(700 + t),
                                 corrupt_layer=max(1, K // 2),
                                 q=BIG_Q)["accepted"] for t in range(3))
    # one larger cube (n=12) end-to-end over the object-dtype field
    x = (1 << 12) - 1
    assert run_chain(x, np.random.default_rng(5), delegate=True,
                     structured=True, q=BIG_Q)["claimed"] == sieve_pi(x)
    #     The lift also CLOSES the wrap-around aliasing hole in the chain's
    #     EXACT trace-region config (build_witness(x, p, nb, dstart=1)): at an n
    #     with x > SMALL_Q a forged quotient u'' aliases X mod SMALL_Q (so the
    #     too-small prime ACCEPTS it) yet violates U*a+R-X=0 over BIG_Q (which
    #     REJECTS it).  The full chain at n>31 is prover-bound (compressed_lucy
    #     is an O(x/log x) Python DP), so the n>31 soundness is exercised here at
    #     the trace certifier -- the seat of the hole -- not the whole chain.
    xa = (1 << 22) - 1
    Va = isqrt(xa); nba = max(1, Va.bit_length())
    primes_a = compressed_lucy(xa)[0]
    pa = primes_a[len(primes_a) // 2]                       # a representative layer prime
    assert xa > SMALL_Q
    W = build_witness(xa, pa, nba, dstart=1, q=SMALL_Q)     # what large_reduce builds
    W2, info = forge_alias(W, SMALL_Q)
    assert W2 is not None and info["u_forged"] != info["u_true"]
    rematerialize(W2, SMALL_Q)
    assert verify_constraints(W2, np.random.default_rng(1),
                              {"t_prover": 0., "t_verifier": 0., "comm": 0},
                              SMALL_Q), "demo prime should alias-accept"
    rematerialize(W2, BIG_Q)
    assert all(not verify_constraints(
        W2, np.random.default_rng(300 + t),
        {"t_prover": 0., "t_verifier": 0., "comm": 0}, BIG_Q)
        for t in range(5)), "lift should reject the alias"

    # 17. FAST BIG_Q PATH (the speed fix): the uint64 Mersenne mulmod gives a
    #     BIT-IDENTICAL chain to the object-array reference of case 16 -- same
    #     claimed pi, accept/reject, comm -- in all three modes, with the chain
    #     cheats and the wrap-around alias still rejected.  (The arithmetic is
    #     exact either way, so equality is the correctness proof; the win is
    #     speed, measured by --bench-big.)
    saved = _cpmt.FAST_BIG
    try:
        for n in (8, 10, 12):
            x = (1 << n) - 1
            truth = sieve_pi(x)
            for (deleg, struct) in [(False, False), (True, False), (True, True)]:
                _cpmt.FAST_BIG = False
                ro = run_chain(x, np.random.default_rng(5), delegate=deleg,
                               structured=struct, q=BIG_Q)
                _cpmt.FAST_BIG = True
                rf = run_chain(x, np.random.default_rng(5), delegate=deleg,
                               structured=struct, q=BIG_Q)
                assert (rf["accepted"] and rf["claimed"] == truth
                        and rf["claimed"] == ro["claimed"]
                        and rf["accepted"] == ro["accepted"]
                        and rf["comm"] == ro["comm"]), (n, deleg, struct, ro, rf)
            K = len(compressed_lucy(x)[0])
            _cpmt.FAST_BIG = True
            assert all(not run_chain(x, np.random.default_rng(600 + t),
                                     cheat="delta_pi", delegate=True,
                                     structured=True, q=BIG_Q)["accepted"]
                       for t in range(3)), (n, "delta_pi fast")
            assert all(not run_chain(x, np.random.default_rng(700 + t),
                                     corrupt_layer=max(1, K // 2), delegate=True,
                                     structured=True, q=BIG_Q)["accepted"]
                       for t in range(3)), (n, "liar fast")
        # the wrap-around alias is still rejected by BIG_Q on the fast path
        _cpmt.FAST_BIG = True
        rematerialize(W2, BIG_Q)
        assert all(not verify_constraints(
            W2, np.random.default_rng(300 + t),
            {"t_prover": 0., "t_verifier": 0., "comm": 0}, BIG_Q)
            for t in range(5)), "fast lift should reject the alias"
    finally:
        _cpmt.FAST_BIG = saved

    # 18. BATCHED TRACE (S502): replacing the K per-layer trace zero-tests with
    #     ONE batched sum-check (batch_trace=True) leaves the chain's verdict
    #     UNCHANGED -- honest accepts & matches the sieve, delta_pi and the
    #     self-consistent liar are still rejected -- over q AND BIG_Q, automaton
    #     and delegated+structured.  (Transcript is NOT bit-identical: the
    #     batched check consumes the rng before the layer loop, so only the
    #     verdict/claimed pi are asserted, not comm.)
    for n in (8, 10, 12):
        x = (1 << n) - 1
        truth = sieve_pi(x)
        K = len(compressed_lucy(x)[0])
        for qf in (Q, BIG_Q):
            for (deleg, struct) in [(False, False), (True, True)]:
                r = run_chain(x, np.random.default_rng(5), delegate=deleg,
                              structured=struct, q=qf, batch_trace=True)
                assert r["accepted"] and r["claimed"] == truth, \
                    (n, qf, deleg, struct, r)
            assert all(not run_chain(x, np.random.default_rng(600 + t),
                                     cheat="delta_pi", q=qf,
                                     batch_trace=True)["accepted"]
                       for t in range(3)), (n, qf, "delta_pi batch")
            assert all(not run_chain(x, np.random.default_rng(700 + t),
                                     corrupt_layer=max(1, K // 2), q=qf,
                                     batch_trace=True)["accepted"]
                       for t in range(3)), (n, qf, "liar batch")

    # 19. BATCHED WIRING (S503): replacing the K delegated per-layer inner wiring
    #     chains (inner_verify_div, 76% of all sum-check CALLS by the S501 profile)
    #     with ONE batched chain -- DEFERRED via batch_wiring=True and discharged by
    #     batched_wiring.verify_wiring_obligations after the layer loop -- leaves the
    #     chain's verdict UNCHANGED: honest accepts & matches the sieve; delta_pi
    #     and the self-consistent liar still rejected.  Tested over q AND BIG_Q,
    #     structured AND dense, ALONE and COMPOSED with batch_trace (the full target
    #     config).  (Transcript not bit-identical: the batched wiring draws its rng
    #     after the layer loop, so only verdict/claimed pi are asserted, not comm.)
    for n in (8, 10, 12):
        x = (1 << n) - 1
        truth = sieve_pi(x)
        K = len(compressed_lucy(x)[0])
        for qf in (Q, BIG_Q):
            for struct in (False, True):
                for bt in (False, True):
                    r = run_chain(x, np.random.default_rng(5), delegate=True,
                                  structured=struct, q=qf, batch_trace=bt,
                                  batch_wiring=True)
                    assert r["accepted"] and r["claimed"] == truth, \
                        (n, qf, struct, bt, r)
            assert all(not run_chain(x, np.random.default_rng(800 + t),
                                     delegate=True, cheat="delta_pi", q=qf,
                                     batch_trace=True,
                                     batch_wiring=True)["accepted"]
                       for t in range(3)), (n, qf, "delta_pi batch_wiring")
            assert all(not run_chain(x, np.random.default_rng(900 + t),
                                     delegate=True, corrupt_layer=max(1, K // 2),
                                     q=qf, batch_trace=True,
                                     batch_wiring=True)["accepted"]
                       for t in range(3)), (n, qf, "liar batch_wiring")
    # batch_wiring requires delegate (the wirings it batches only exist delegated)
    try:
        run_chain((1 << 8) - 1, np.random.default_rng(0), batch_wiring=True)
        assert False, "batch_wiring without delegate must raise"
    except ValueError:
        pass

    # 19b. The deferred/batched wiring still rejects the WIRING-SPECIFIC liars
    #      (forged wiring scalar / self-consistent lying inner chain), exercised at
    #      the two-sided single layer (run_chain injects only chain-level cheats).
    #      waff_forge is caught in the per-layer comparator fold (before deferral);
    #      small_forge by sum-check #0 of the batch; small_chain / waff_chain by the
    #      batched backward sweep -- all THROUGH verify_wiring_obligations.
    for qf in (Q, BIG_Q):
        x = (1 << 12) - 1
        K = len(compressed_lucy(x)[0])
        layer = max(1, K // 2)
        assert run_two_sided_layer(x, layer, np.random.default_rng(5),
                                   delegate=True, q=qf,
                                   batch_wiring=True)["accepted"], qf
        for ch in ("small_forge", "small_chain", "waff_chain", "waff_forge"):
            assert all(not run_two_sided_layer(
                x, layer, np.random.default_rng(1000 + t), cheat=ch,
                delegate=True, q=qf, batch_wiring=True)["accepted"]
                for t in range(4)), (qf, ch)

    # 20. PCS LEAF OPENINGS (S505): pcs=True replaces the per-layer carried-claim
    #     leaf openings (line_batch_pair / batch_on_table folds + the redundant
    #     s2/s_B closes) with real sum-check openings (leaf_open.open_batch) whose
    #     residuals thread to the base.  The chain VERDICT is unchanged -- honest
    #     accepts & matches the sieve; delta_pi and the self-consistent liar still
    #     rejected -- over q AND BIG_Q, automaton AND delegated+structured, and
    #     COMPOSED with batch_trace + batch_wiring (the full winning config).  The
    #     residual-threading grounds every carried scalar at the S_0 base; the
    #     liar is caught in phase-A round-1 (untouched), so soundness is preserved.
    #     (Transcript not bit-identical -- a sum-check fold replaces a line
    #     restriction -- so only verdict/claimed pi are asserted.)
    for n in (8, 10, 12):
        x = (1 << n) - 1
        truth = sieve_pi(x)
        K = len(compressed_lucy(x)[0])
        for qf in (Q, BIG_Q):
            # honest: automaton, delegated+structured, and the full batched config
            cfgs = [dict(delegate=False),
                    dict(delegate=True, structured=True),
                    dict(delegate=True, structured=True, batch_trace=True,
                         batch_wiring=True)]
            for kw in cfgs:
                r = run_chain(x, np.random.default_rng(5), q=qf, pcs=True, **kw)
                assert r["accepted"] and r["claimed"] == truth, \
                    (n, qf, kw, "pcs honest")
            # delta_pi and self-consistent liar rejected under pcs (full config)
            assert all(not run_chain(x, np.random.default_rng(600 + t),
                                     cheat="delta_pi", delegate=True,
                                     structured=True, q=qf, batch_trace=True,
                                     batch_wiring=True, pcs=True)["accepted"]
                       for t in range(3)), (n, qf, "pcs delta_pi")
            for i0 in sorted(set([K, max(1, K // 2), 1])):
                assert all(not run_chain(x, np.random.default_rng(700 + i0 + t),
                                         corrupt_layer=i0, delegate=True,
                                         structured=True, q=qf, batch_trace=True,
                                         batch_wiring=True, pcs=True)["accepted"]
                           for t in range(3)), (n, qf, "pcs liar", i0)

    # 21. BATCHED Ub OPENINGS (S506): batch_ub=True defers verify_trace_region's nb
    #     Ub-bit-table openings -- the LAST per-layer O(2^nb) verifier term left after
    #     pcs (S505) -- and discharges them in ONE sum-check against the SAME stacked
    #     Ub cube batch_trace commits (batched_trace.verify_ub_openings_batched).  The
    #     chain VERDICT is unchanged: honest accepts & matches the sieve; delta_pi and
    #     the self-consistent liar still rejected -- over q AND BIG_Q.  The asymptotic
    #     fact, MEASURED: the verifier-side O(2^nb) Ub-leaf eval count drops from K*nb
    #     to EXACTLY 0 (moved to the prover, ub_leaf_p == K*nb); the per-layer verifier
    #     is then leaf-eval-free (the only O(sqrt x) verifier ops left are one-time --
    #     the two S_0 base closes + the batched discharges).  Requires pcs+batch_trace.
    for n in (8, 10, 12):
        x = (1 << n) - 1
        truth = sieve_pi(x)
        K = len(compressed_lucy(x)[0])
        nb = max(1, isqrt(x).bit_length())
        for qf in (Q, BIG_Q):
            cfg = dict(delegate=True, structured=True, batch_trace=True,
                       batch_wiring=True, pcs=True)
            off = run_chain(x, np.random.default_rng(5), q=qf, batch_ub=False, **cfg)
            on = run_chain(x, np.random.default_rng(5), q=qf, batch_ub=True, **cfg)
            assert on["accepted"] and on["claimed"] == truth, (n, qf, "ub honest")
            # verdict identical to pcs-without-ub-batch; the op-count shifts
            assert off["accepted"] and off["claimed"] == on["claimed"]
            assert off.get("ub_leaf_v", 0) == K * nb, (n, qf, off.get("ub_leaf_v"))
            assert on.get("ub_leaf_v", 0) == 0, (n, qf, on.get("ub_leaf_v"))
            assert on.get("ub_leaf_p", 0) == K * nb, (n, qf, on.get("ub_leaf_p"))
            # delta_pi and the self-consistent liar still rejected under batch_ub
            assert all(not run_chain(x, np.random.default_rng(600 + t),
                                     cheat="delta_pi", q=qf, batch_ub=True,
                                     **cfg)["accepted"]
                       for t in range(3)), (n, qf, "ub delta_pi")
            for i0 in sorted(set([K, max(1, K // 2), 1])):
                assert all(not run_chain(x, np.random.default_rng(700 + i0 + t),
                                         corrupt_layer=i0, q=qf, batch_ub=True,
                                         **cfg)["accepted"]
                           for t in range(3)), (n, qf, "ub liar", i0)
    # batch_ub requires pcs AND batch_trace (the committed cube it discharges against)
    for bad in (dict(batch_ub=True),
                dict(batch_ub=True, pcs=True),
                dict(batch_ub=True, batch_trace=True, delegate=True)):
        try:
            run_chain((1 << 8) - 1, np.random.default_rng(0), **bad)
            assert False, ("batch_ub guard", bad)
        except ValueError:
            pass

    # 22. VERIFIER OP-COUNT ATTRIBUTION (S507): the per-region large-table-eval
    #     counter (_acct_vleaf) splits the verifier's O(2^nb) leaf opens into the
    #     PER-LAYER critical path vs the ONE-TIME base closes.  Exact count
    #     identities (the falsifiable instrumentation check), and the leading-term
    #     direction: config (c) (pcs+batch_trace+batch_ub) is per-layer leaf-eval-
    #     FREE so its total is the two S_0 closes ONLY -- collapsing the Theta(x)
    #     per-layer term (a)/(b) carry to a Theta(sqrt x) one-time term.  Over q AND
    #     BIG_Q (the op COUNT is field-independent: same tables, same eval count).
    for n in (8, 10, 12):
        x = (1 << n) - 1
        truth = sieve_pi(x)
        K = len(compressed_lucy(x)[0])
        nb = max(1, isqrt(x).bit_length())
        for qf in (Q, BIG_Q):
            ra = run_chain(x, np.random.default_rng(5), q=qf,
                           delegate=True, structured=True)               # (a)
            rb = run_chain(x, np.random.default_rng(5), q=qf,
                           delegate=True, structured=True, pcs=True)       # (b)
            rc = run_chain(x, np.random.default_rng(5), q=qf, delegate=True,
                           structured=True, pcs=True, batch_trace=True,
                           batch_ub=True)                                  # (c)
            for r in (ra, rb, rc):
                assert r["accepted"] and r["claimed"] == truth, (n, qf)
                assert r["vleaf_cnt_ot"] == 2                  # two S_0 base closes
                assert r["vleaf_ops_ot"] == 2 * (1 << nb)
            # exact per-layer leaf-eval COUNTS
            assert ra["vleaf_cnt_pl"] == K * (nb + 5) - 1, (n, qf, ra["vleaf_cnt_pl"])
            assert rb["vleaf_cnt_pl"] == K * nb, (n, qf, rb["vleaf_cnt_pl"])
            assert rc["vleaf_cnt_pl"] == 0, (n, qf, rc["vleaf_cnt_pl"])
            # field-independence of the op count (same table sizes, same evals)
            if qf == BIG_Q:
                assert ra["vleaf_ops_pl"] == ra0_pl and rc["vleaf_ops_pl"] == 0
            else:
                ra0_pl = ra["vleaf_ops_pl"]
            # all per-layer leaves are 2^nb-size, so ops = cnt * 2^nb exactly
            assert ra["vleaf_ops_pl"] == (K * (nb + 5) - 1) * (1 << nb)
            assert rb["vleaf_ops_pl"] == K * nb * (1 << nb)
            # leading-term direction: (c) total is the one-time base ONLY
            tot_a = ra["vleaf_ops_pl"] + ra["vleaf_ops_ot"]
            tot_c = rc["vleaf_ops_pl"] + rc["vleaf_ops_ot"]
            assert tot_c == 2 * (1 << nb)
            assert tot_a > K * tot_c // 4              # per-layer term dominates (a)

    # 23. COMMIT-BASE (S508): the two one-time S_0 base closes discharged through
    #     the tensor polynomial commitment (pcs_commit) -- verdict UNCHANGED vs the
    #     direct mle_eval close (config (c)), per-layer still leaf-eval-free, and the
    #     ONLY one-time x-scaling verifier work is now the sub-sqrt-x commit opens
    #     (vcommit_ops == vleaf_ops_ot, < 2*2^nb asymptotically -- here just the
    #     instrumentation identity + the correctness/soundness equivalence).  Over q
    #     AND BIG_Q (the commitment is field-parameterised through pcs_commit).
    full = dict(delegate=True, structured=True, pcs=True,
                batch_trace=True, batch_ub=True)
    for n in (8, 10, 12):
        x = (1 << n) - 1
        truth = sieve_pi(x)
        K = len(compressed_lucy(x)[0])
        # honest accept + instrumentation identity over q AND BIG_Q (the commitment
        # is field-parameterised through pcs_commit)
        for qf in (Q, BIG_Q):
            cb = run_chain(x, np.random.default_rng(5), q=qf, commit_base=True, **full)
            assert cb["accepted"] and cb["claimed"] == truth, (n, qf, "commit-base")
            assert cb["vleaf_cnt_pl"] == 0                  # per-layer still leaf-free
            assert cb["vleaf_cnt_ot"] == 2                  # two base opens
            assert cb["vcommit_ops"] == cb["vleaf_ops_ot"] > 0
        # cheat rejection THROUGH the committed base (cheat mechanism is
        # field-independent -- already covered over BIG_Q in case 22 -- so run the
        # panel over q to keep the selftest fast)
        assert all(not run_chain(x, np.random.default_rng(600 + t),
                                 cheat="delta_pi", commit_base=True,
                                 **full)["accepted"] for t in range(3)), (n,)
        for i0 in sorted(set([1, max(1, K // 2), K])):
            assert all(not run_chain(x, np.random.default_rng(700 + i0 + t),
                                     corrupt_layer=i0, commit_base=True,
                                     **full)["accepted"]
                       for t in range(3)), (n, "liar", i0)

    # 24. EXPANDER BASE (S514): the field-size-free Brakedown-style row code is a
    #     drop-in for the RS commitment behind the SAME interface -- VERDICT
    #     UNCHANGED (honest pi==sieve, delta_pi + self-consistent liar rejected),
    #     per-layer still leaf-free, base opens still sub-sqrt-x (vcommit_ops ==
    #     vleaf_ops_ot).  Over q AND BIG_Q.
    for n in (8, 10, 12):
        x = (1 << n) - 1
        truth = sieve_pi(x)
        K = len(compressed_lucy(x)[0])
        for qf in (Q, BIG_Q):
            cb = run_chain(x, np.random.default_rng(5), q=qf, commit_base=True,
                           commit_code="expander", **full)
            assert cb["accepted"] and cb["claimed"] == truth, (n, qf, "expander-base")
            assert cb["vleaf_cnt_pl"] == 0
            assert cb["vleaf_cnt_ot"] == 2
            assert cb["vcommit_ops"] == cb["vleaf_ops_ot"] > 0
        assert all(not run_chain(x, np.random.default_rng(800 + t),
                                 cheat="delta_pi", commit_base=True,
                                 commit_code="expander", **full)["accepted"]
                   for t in range(2)), (n, "expander delta_pi")
        i0 = max(1, K // 2)
        assert not run_chain(x, np.random.default_rng(900 + i0),
                             corrupt_layer=i0, commit_base=True,
                             commit_code="expander", **full)["accepted"], \
            (n, "expander liar", i0)

    print("selftest OK")


def bench():
    """Full compressed chain: automaton wiring vs DELEGATED wiring. Headline is
    t_verifier -- in delegate mode it sheds the O(K nb p) automaton term and
    scales ~O(K nb log p), at the cost of more communication (chain transcripts)
    and more prover work (both still Õ(sqrt x)-compatible)."""
    print(f"{'n':>3} {'V':>5} {'nb':>3} {'K':>4} {'pi(x)':>7} "
          f"{'tV_auto_ms':>11} {'tV_deleg_ms':>12} {'tV_ratio':>9} "
          f"{'comm_auto':>10} {'comm_deleg':>11}")
    for n in (8, 10, 12, 14, 16):
        x = (1 << n) - 1
        a = run_chain(x, np.random.default_rng(n), delegate=False)
        d = run_chain(x, np.random.default_rng(n), delegate=True)
        assert a["accepted"] and a["claimed"] == sieve_pi(x)
        assert d["accepted"] and d["claimed"] == sieve_pi(x)
        ratio = a["t_verifier"] / d["t_verifier"] if d["t_verifier"] else 0
        print(f"{n:>3} {isqrt(x):>5} {max(1, isqrt(x).bit_length()):>3} "
              f"{a['layers']:>4} {a['claimed']:>7} "
              f"{a['t_verifier']*1000:>11.3f} {d['t_verifier']*1000:>12.3f} "
              f"{ratio:>9.2f} {a['comm']:>10} {d['comm']:>11}")


def bench_pcs():
    """The S505 measurement: pcs=True replaces the per-layer carried-claim LEAF
    openings (line_batch_pair / batch_on_table folds + the redundant s2/s_B closes)
    with sum-check openings whose residuals thread to the base.  Headline is the
    VERIFIER term -- pcs removes ~5 direct O(2^nb) leaf closes per layer, leaving
    only verify_trace_region's nb Ub-bit-table openings (line 429) as the remaining
    per-layer O(2^nb) cost (the documented batched-trace follow-on).  So this is a
    CONSTANT-FACTOR verifier win, NOT yet the asymptotic Õ(sqrt x): the per-layer
    leaf cost drops from ~(nb+5) to ~nb direct openings.  Comm rises a little (a
    degree-2 sum-check fold has more round scalars than a line poly).  Reports
    delegated chain, pcs off vs on.  Falsifier: pcs t_verifier not below pcs-off,
    or claimed pi changing."""
    print("delegated compressed chain: leaf-opening stand-in (pcs off) vs "
          "sum-check openings (pcs on)")
    print(f"{'n':>3} {'V':>6} {'nb':>3} {'K':>4} {'pi(x)':>7} "
          f"{'tV_off_ms':>10} {'tV_pcs_ms':>10} {'tV_ratio':>9} "
          f"{'comm_off':>9} {'comm_pcs':>9}")
    for n in (8, 10, 12, 14, 16):
        x = (1 << n) - 1
        a = run_chain(x, np.random.default_rng(n), delegate=True, pcs=False)
        b = run_chain(x, np.random.default_rng(n), delegate=True, pcs=True)
        assert a["accepted"] and a["claimed"] == sieve_pi(x)
        assert b["accepted"] and b["claimed"] == sieve_pi(x)
        ratio = a["t_verifier"] / b["t_verifier"] if b["t_verifier"] else 0
        print(f"{n:>3} {isqrt(x):>6} {max(1, isqrt(x).bit_length()):>3} "
              f"{a['layers']:>4} {a['claimed']:>7} "
              f"{a['t_verifier']*1000:>10.3f} {b['t_verifier']*1000:>10.3f} "
              f"{ratio:>9.2f} {a['comm']:>9} {b['comm']:>9}")


def bench_ub():
    """The S506 measurement: batch_ub=True defers verify_trace_region's nb Ub-bit-
    table openings (the LAST per-layer O(2^nb) VERIFIER term after pcs) and discharges
    them in ONE sum-check against the committed stacked Ub cube.  The ASYMPTOTIC
    headline is the COUNT 'ub_leaf_v' of O(2^nb) Ub-leaf evals on the verifier's
    per-layer critical path: K*nb with batch_ub off, EXACTLY 0 with it on (they move
    to the prover, 'ub_leaf_p').  That is the falsifiable signal that the per-layer
    verifier is now leaf-eval-free (Õ(sqrt x) end-to-end -- only one-time S_0 closes +
    batched discharges remain O(sqrt x)).

    Honest scope: the WALL-CLOCK t_verifier drop is a modest ~1.1-1.2x and roughly
    flat in n -- because a vectorised numpy mle_eval on an O(sqrt x)-size array is
    cheap relative to the Python-loop wiring/eq recomputes that dominate the measured
    t_verifier, the asymptotic O(K nb 2^nb)->one-time improvement does NOT surface as
    a growing wall-clock ratio at reachable n (it would past the array-size crossover
    S500 measured).  t_prover rises by the shifted openings.  This CORRECTS the S505
    NEXT-ACTION prediction that the ratio would grow with n.  Falsifier: ub_leaf_v
    not 0 under batch_ub, or claimed pi changing, or t_verifier(on) > t_verifier(off).
    Full config (delegate+structured+batch_trace+batch_wiring+pcs), q field."""
    print("compressed chain, full config: per-layer Ub openings inline (batch_ub "
          "off) vs ONE batched discharge (on)")
    print(f"{'n':>3} {'V':>6} {'nb':>3} {'K':>4} {'pi(x)':>7} "
          f"{'Knb':>6} {'vleaf_off':>9} {'vleaf_on':>8} "
          f"{'tV_off_ms':>10} {'tV_on_ms':>9} {'tVr':>5} "
          f"{'tP_off_ms':>10} {'tP_on_ms':>9} {'comm_on':>8}")
    cfg = dict(delegate=True, structured=True, batch_trace=True,
               batch_wiring=True, pcs=True)
    for n in (8, 10, 12, 14, 16):
        x = (1 << n) - 1
        a = run_chain(x, np.random.default_rng(n), batch_ub=False, **cfg)
        b = run_chain(x, np.random.default_rng(n), batch_ub=True, **cfg)
        assert a["accepted"] and b["accepted"]
        assert a["claimed"] == b["claimed"] == sieve_pi(x)
        nb = max(1, isqrt(x).bit_length())
        K = a["layers"]
        r = a["t_verifier"] / b["t_verifier"] if b["t_verifier"] else 0
        print(f"{n:>3} {isqrt(x):>6} {nb:>3} {K:>4} {a['claimed']:>7} "
              f"{K*nb:>6} {a.get('ub_leaf_v',0):>9} {b.get('ub_leaf_v',0):>8} "
              f"{a['t_verifier']*1000:>10.3f} {b['t_verifier']*1000:>9.3f} "
              f"{r:>4.2f}x {a['t_prover']*1000:>10.3f} {b['t_prover']*1000:>9.3f} "
              f"{b['comm']:>8}")


def _fit_slope(xs, ys):
    """Least-squares slope of ys vs xs (xs = n, ys = log2(metric))."""
    k = len(xs)
    mx = sum(xs) / k
    my = sum(ys) / k
    num = sum((a - mx) * (b - my) for a, b in zip(xs, ys))
    den = sum((a - mx) ** 2 for a in xs)
    return num / den if den else 0.0


def bench_verifier_ops(ns=(8, 10, 12, 14, 16, 18), seed=1):
    """S507 -- the whole-chain VERIFIER op-count CURVE, the honest headline the
    wall-clock t_verifier hides (the dominant sum-check FOLDING is untimed prover
    work, S501; vectorised numpy flattens the asymptotics, S506).

    Metric `vleaf_ops` = the count of VERIFIER large-table (mle_eval, O(2^nb) =
    O(sqrt x)) evaluations, weighted by table size -- the ONLY verifier operation
    whose per-call cost grows with x.  Every OTHER verifier step (sum-check round
    checks, eq/le/band point evals O(nb), the delegated wiring carry chain
    O(nb log p), the batched-discharge recomputes O(K polylog)) is polylog or
    O(K polylog), bounded by Õ(sqrt x); the printed `comm` curve corroborates that
    residual stays Õ(sqrt x) in ALL configs.  So vleaf_ops IS the leading term.

    Four configs (ALL delegate+structured -- so the wiring verifier is already
    polylog and the leaf opens are the SOLE x-scaling verifier term):
      (a) no-pcs              per-layer leaf closes  K*(nb+5)-1 evals  ~ Theta(x)
      (b) pcs                 only Ub openings survive   K*nb evals    ~ Theta(x)
      (c) pcs+batch_trace+ub  per-layer leaf-eval-FREE; one-time 2 evals ~ Th(sqrt x)
      (d) +commit_base        one-time base via tensor PCS (S508) ~ Th(x^{1/4})

    FALSIFIABLE CLAIM: (c)'s whole-chain verifier leaf-op count has leading term
    Theta(sqrt x polylog) -- slope of log2(ops) vs n approaches 0.5 -- while (a)/(b)
    are Theta(x polylog) -- slope -> 1.0.  (d) discharges the one-time base closes
    through the sub-sqrt-x tensor commitment, so its leading term drops BELOW sqrt x
    -- slope -> 0.25 (Theta(x^{1/4})), i.e. the whole-chain verifier is sub-sqrt x
    (succinct under the hash assumption).  Falsifier: (c)/(d) per-layer leaf count
    != 0; (c)'s slope not -> 0.5; (d)'s slope not -> 0.25 (not sub-(c)); (a)/(b)'s
    slope not -> 1.0; or claimed pi(x) differing across configs."""
    cfgs = [("a no-pcs", dict(delegate=True, structured=True)),
            ("b pcs", dict(delegate=True, structured=True, pcs=True)),
            ("c pcs+bt+ub", dict(delegate=True, structured=True, pcs=True,
                                 batch_trace=True, batch_ub=True)),
            ("d +commit", dict(delegate=True, structured=True, pcs=True,
                               batch_trace=True, batch_ub=True, commit_base=True))]
    print("whole-chain VERIFIER large-table-eval op count (field-ops mle_eval "
          "touches), default q")
    print(f"{'n':>3} {'x':>9} {'V':>6} {'nb':>3} {'K':>4} {'pi':>6} | "
          + " | ".join(f"{nm:>11}: {'pl_cnt':>6} {'total':>9} {'x4?':>4}"
                       for nm, _ in cfgs))
    data = {nm: [] for nm, _ in cfgs}                # (n, x, total, comm, r)
    for n in ns:
        x = (1 << n) - 1
        V = isqrt(x); nb = max(1, V.bit_length())
        K = len(compressed_lucy(x)[0]); truth = sieve_pi(x)
        cells = []
        for nm, kw in cfgs:
            r = run_chain(x, np.random.default_rng(seed), **kw)
            assert r["accepted"] and r["claimed"] == truth, (n, nm)
            tot = r["vleaf_ops_pl"] + r["vleaf_ops_ot"]
            # exact per-config leaf-eval COUNT identities (instrumentation check)
            if nm.startswith("a"):
                assert r["vleaf_cnt_pl"] == K * (nb + 5) - 1, (n, r["vleaf_cnt_pl"])
            elif nm.startswith("b"):
                assert r["vleaf_cnt_pl"] == K * nb, (n, r["vleaf_cnt_pl"])
            elif nm.startswith("c"):
                assert r["vleaf_cnt_pl"] == 0, (n, r["vleaf_cnt_pl"])
                assert r["vleaf_ops_ot"] == 2 * (1 << nb), (n, r["vleaf_ops_ot"])
            else:                       # (d) commit_base: per-layer free + sub-sqrt-x base
                assert r["vleaf_cnt_pl"] == 0, (n, r["vleaf_cnt_pl"])
                assert r["vcommit_ops"] == r["vleaf_ops_ot"], (n,)   # only one-time work
            assert r["vleaf_cnt_ot"] == 2, (n, r["vleaf_cnt_ot"])    # 2 S_0 closes
            prev = data[nm][-1][2] if data[nm] else 0
            fac = tot / prev if prev else 0.0
            data[nm].append((n, x, tot, r["comm"], r))
            cells.append(f"{nm:>11}: {r['vleaf_cnt_pl']:>6} {tot:>9} "
                         + (f"{fac:>3.1f}x" if fac else f"{'-':>4}"))
        print(f"{n:>3} {x:>9} {V:>6} {nb:>3} {K:>4} {truth:>6} | "
              + " | ".join(cells))

    import math
    print("\nfitted leading-term exponent  alpha  (total_ops ~ x^alpha; "
          "log2(total) vs n, last 4 pts):")
    print(f"{'config':>14} {'alpha_ops':>10} {'alpha_comm':>11}  expectation")
    expect = {"a no-pcs": "ops~x  (1.0), comm~sqrt x",
              "b pcs": "ops~x  (1.0), comm~sqrt x",
              "c pcs+bt+ub": "ops~sqrt x (0.5), comm~sqrt x",
              "d +commit": "ops~x^1/4 (0.25), comm~sqrt x"}
    for nm, _ in cfgs:
        pts = data[nm][-4:]
        nsl = [p[0] for p in pts]
        a_ops = _fit_slope(nsl, [math.log2(p[2]) for p in pts])
        a_comm = _fit_slope(nsl, [math.log2(p[3]) for p in pts])
        print(f"{nm:>14} {a_ops:>10.3f} {a_comm:>11.3f}  {expect[nm]}")
    print("\n(x = 2^n-1, so alpha = d log2(total)/d n: alpha~1.0 => Theta(x), "
          "alpha~0.5 => Theta(sqrt x).  Per-step factor over dn=2: ~4x for x, "
          "~2x for sqrt x.)")
    # headline assertion: (c) leading exponent is sub-linear and (a)/(b) linear
    a_a = _fit_slope([p[0] for p in data["a no-pcs"][-4:]],
                     [math.log2(p[2]) for p in data["a no-pcs"][-4:]])
    a_c = _fit_slope([p[0] for p in data["c pcs+bt+ub"][-4:]],
                     [math.log2(p[2]) for p in data["c pcs+bt+ub"][-4:]])
    a_d = _fit_slope([p[0] for p in data["d +commit"][-4:]],
                     [math.log2(p[2]) for p in data["d +commit"][-4:]])
    assert a_a > 0.9, ("config a not ~linear", a_a)
    assert a_c < 0.6, ("config c not ~sqrt", a_c)
    assert a_d < 0.4 and a_d < a_c, ("config d not sub-sqrt-x", a_d, a_c)
    print(f"\nHEADLINE: (a) alpha_ops = {a_a:.3f} ~ 1.0 (Theta(x));  "
          f"(c) alpha_ops = {a_c:.3f} ~ 0.5 (Theta(sqrt x));  "
          f"(d) alpha_ops = {a_d:.3f} ~ 0.25 (Theta(x^1/4), SUB-sqrt x).  "
          f"The whole-chain verifier leading term drops Theta(x) -> Theta(sqrt x) "
          f"(per-layer leaf opens removed) -> Theta(x^1/4) (base via tensor PCS).")


def wiring_bench(nb=12, trials=5):
    """The p-isolation measurement: per-layer WIRING-check verifier cost vs p
    at FIXED nb. Automaton (waff_eval / w_div_eval) is O(nb*p); the delegated
    chain (verify_waff_value / inner_verify_div) is O(nb*log p). The delegated
    columns must FLATTEN in p (track l = ceil(log2 p)) while the automaton
    columns grow ~linearly. Falsifier: delegated verifier time scaling with p."""
    rng = np.random.default_rng(7)
    Dml = 1 << nb
    V = Dml - 1                           # mirror a full nb-bit layer (x ~ V^2)
    print(f"nb = {nb}  (compressed cube size 2^nb = {Dml}, V = {V})")
    print(f"{'p':>5} {'l':>3} || {'waff_auto_ms':>13} {'waff_deleg_ms':>14} "
          f"{'ratio':>6} || {'div_auto_ms':>12} {'div_deleg_ms':>13} "
          f"{'ratio':>6} || {'deleg_comm':>10}")
    for p in (3, 7, 13, 31, 61, 127, 251, 509, 1021):
        ecut = V // p - 1                 # the affine region is e <= ecut
        if ecut < 0:                      # p > V: no affine region
            continue
        l = max(1, (p - 1).bit_length())
        ta = td = ta2 = td2 = 0.0
        comm = 0
        for t in range(trials):
            r_v = [int(rng.integers(0, Q)) for _ in range(nb)]
            r_u = [int(rng.integers(0, Q)) for _ in range(nb)]
            t0 = time.perf_counter()
            waffv = waff_eval(r_v, r_u, p, ecut, nb)            # O(nb*p)
            ta += time.perf_counter() - t0
            st = {"t_prover": 0.0, "t_verifier": 0.0, "comm": 0}
            assert verify_waff_value(waffv, r_v, r_u, p, ecut, nb,
                                     np.random.default_rng(100 + t), st)
            td += st["t_verifier"]
            t0 = time.perf_counter()
            divv = w_div_eval(r_v, r_u, p, nb)                  # O(nb*p)
            ta2 += time.perf_counter() - t0
            st2 = {"t_prover": 0.0, "t_verifier": 0.0, "comm": 0}
            assert inner_verify_div(divv, r_v, r_u, p, nb,
                                    np.random.default_rng(200 + t), st2)
            td2 += st2["t_verifier"]
            comm = st["comm"] + st2["comm"]
        r1 = ta / td if td else 0
        r2 = ta2 / td2 if td2 else 0
        print(f"{p:>5} {l:>3} || {ta/trials*1000:>13.3f} "
              f"{td/trials*1000:>14.3f} {r1:>6.2f} || "
              f"{ta2/trials*1000:>12.3f} {td2/trials*1000:>13.3f} {r2:>6.2f} "
              f"|| {comm:>10}")


def prover_bench(nb=14, trials=2):
    """The compressed chain's WIRING-PROVER scaling in p at fixed nb (mirrors
    lucy_dp_delegated_wiring.prover_bench, S496).  The delegated inner sweep,
    dense (five 2^{2l} tables, O(nb*p^2)) vs structured (chain_layer_reduce_
    structured, O(nb*p)).  BOTH wirings the compressed chain delegates are
    measured: the small DIVISION (inner_verify_div accept_rem=None ==
    w_div_eval) and the large AFFINE image (accept_rem=p-1, the verify_waff_
    value inner claim).  Transcript is bit-identical (selftest 15), so this is
    pure scaling: dense/p^2 -> const, struct/p -> const, ratio ~ p.  Whole-chain
    wiring prover Sum_{p<=sqrt x} nb p^2 ~ Õ(x^{3/2})  ->  Sum nb p ~ Õ(x).
    Falsifier: structured column not flattening to ~p (struct/p not ~const)."""
    rng = np.random.default_rng(7)
    print(f"nb = {nb}  (delegated chain wiring prover; both wirings)")
    print(f"{'p':>5} {'l':>3} || {'div_dense':>10} {'div_struct':>11} "
          f"{'r':>6} {'d/p^2':>8} {'s/p':>8} || {'aff_dense':>10} "
          f"{'aff_struct':>11} {'r':>6} {'d/p^2':>8} {'s/p':>8}")
    for p in (7, 13, 31, 61, 127, 251, 509, 1021):
        l = max(1, (p - 1).bit_length())
        r_v = [int(rng.integers(0, Q)) for _ in range(nb)]   # dividend bits
        r_u = [int(rng.integers(0, Q)) for _ in range(nb)]   # quotient bits
        vs = inner_chain_vectors(r_v, r_u, p, l, nb)
        divv = int(sum(int(vs[nb][s]) for s in range(p)) % Q)   # == w_div_eval
        affv = int(vs[nb][p - 1] % Q)                           # accept rem p-1

        def timed(accept_rem, claim, structured):
            t0 = time.perf_counter()
            for t in range(trials):
                st = {"t_prover": 0.0, "t_verifier": 0.0, "comm": 0}
                ok = inner_verify_div(claim, r_v, r_u, p, nb,
                                      np.random.default_rng(200 + t), st,
                                      accept_rem=accept_rem,
                                      structured=structured)
                assert ok, (p, accept_rem, structured)
            return (time.perf_counter() - t0) / trials

        dd = timed(None, divv, False)
        ds = timed(None, divv, True)
        ad = timed(p - 1, affv, False)
        as_ = timed(p - 1, affv, True)
        print(f"{p:>5} {l:>3} || {dd*1000:>10.2f} {ds*1000:>11.2f} "
              f"{dd/ds:>6.2f} {dd/p**2*1e6:>8.4f} {ds/p*1e6:>8.4f} || "
              f"{ad*1000:>10.2f} {as_*1000:>11.2f} {ad/as_:>6.2f} "
              f"{ad/p**2*1e6:>8.4f} {as_/p*1e6:>8.4f}")


def bench_big(seed=1, ns=(8, 10, 12, 14, 16)):
    """The chain speed fix, measured: the full compressed delegated+structured
    pi(x) chain over BIG_Q = 2^61-1 via the uint64 Mersenne path vs the object-
    array reference.  run_chain's wall is field arithmetic (the Lucy DP is
    <0.1%), so this is the headline speedup for a sound (n<~60) large-x run.
    Asserts the two paths agree on claimed pi(x)."""
    saved = _cpmt.FAST_BIG
    print("full compressed chain over BIG_Q = 2^61-1 (delegate+structured): "
          "uint64 Mersenne vs object arrays")
    print(f"{'n':>3} {'x':>10} {'V':>7} {'K':>5} {'obj_ms':>10} {'fast_ms':>10} "
          f"{'speedup':>8}  {'pi==sieve':>9}")
    try:
        for n in ns:
            x = (1 << n) - 1
            V = isqrt(x); K = len(compressed_lucy(x)[0]); truth = sieve_pi(x)
            _cpmt.FAST_BIG = False
            t0 = time.perf_counter()
            ro = run_chain(x, np.random.default_rng(seed), delegate=True,
                           structured=True, q=BIG_Q)
            t_obj = time.perf_counter() - t0
            _cpmt.FAST_BIG = True
            t0 = time.perf_counter()
            rf = run_chain(x, np.random.default_rng(seed), delegate=True,
                           structured=True, q=BIG_Q)
            t_fast = time.perf_counter() - t0
            ok = (ro["claimed"] == rf["claimed"] == truth
                  and ro["accepted"] and rf["accepted"])
            print(f"{n:>3} {x:>10} {V:>7} {K:>5} {t_obj*1000:>10.1f} "
                  f"{t_fast*1000:>10.1f} {t_obj/max(t_fast,1e-9):>7.1f}x  "
                  f"{'yes' if ok else 'NO!':>9}")
            assert ok, (n, ro["claimed"], rf["claimed"], truth)
    finally:
        _cpmt.FAST_BIG = saved


def bench_combined(n=16, seed=1):
    """The S503 headline: end-to-end run_chain over BIG_Q with BOTH big kernels
    batched (batch_trace S502 + batch_wiring S503).  The S501 profile found the
    trace zero-test (53% of wall) and the wiring delegation (30% of wall, 76% of
    sum-check calls) were the two op-count-bound kernels that made FAST_BIG a NET
    LOSS unbatched (S502: 22 s vs 16 s, because the tiny per-layer wiring cubes
    penalised the 24-op Mersenne mulmod).  With BOTH widened, globally enabling the
    fast path should finally win end-to-end.  Reports wall across configs; asserts
    every config returns the correct pi(x).  Falsifier: batch_trace+wiring+FAST
    not beating the baseline, or any config mis-counting pi(x)."""
    x = (1 << n) - 1
    V = isqrt(x); K = len(compressed_lucy(x)[0]); truth = sieve_pi(x)
    saved = _cpmt.FAST_BIG
    configs = [
        ("baseline (no batch, obj)",  dict(batch_trace=False, batch_wiring=False), False),
        ("batch_trace (obj)",         dict(batch_trace=True,  batch_wiring=False), False),
        ("batch_trace+wiring (obj)",  dict(batch_trace=True,  batch_wiring=True),  False),
        ("batch_trace (FAST)",        dict(batch_trace=True,  batch_wiring=False), True),
        ("batch_trace+wiring (FAST)", dict(batch_trace=True,  batch_wiring=True),  True),
    ]
    print(f"end-to-end run_chain over BIG_Q=2^61-1 (delegate+structured), "
          f"n={n}: x=2^{n}-1, V={V}, K={K}, pi={truth}")
    print(f"{'config':>28} {'wall_ms':>10} {'comm':>10} {'pi==sieve':>10}")
    base_ms = None
    try:
        for label, kw, fast in configs:
            _cpmt.FAST_BIG = fast
            t0 = time.perf_counter()
            r = run_chain(x, np.random.default_rng(seed), delegate=True,
                          structured=True, q=BIG_Q, **kw)
            wall = (time.perf_counter() - t0) * 1000
            if base_ms is None:
                base_ms = wall
            ok = r["accepted"] and r["claimed"] == truth
            print(f"{label:>28} {wall:>10.1f} {r['comm']:>10} "
                  f"{('yes' if ok else 'NO!'):>10}  ({base_ms/max(wall,1e-9):.2f}x vs base)")
            assert ok, (label, r["claimed"], truth)
    finally:
        _cpmt.FAST_BIG = saved


def main(n, layer, cheat_trials, seed, delegate, structured=False, field="q",
         fast_big_flag=False, pcs=False, commit_base=False, commit_code="rs"):
    q = FIELDS[field]
    # The uint64 Mersenne path wins only on LARGE arrays; the chain is dominated
    # by many SMALL per-layer cubes, so it is SLOWER at reachable n (see
    # --bench-big and the results file).  Keep BIG_Q on the object reference by
    # default; --fast-big opts into the uint64 path (pays off only at large n,
    # where the sqrt(x)-cube layer sum-checks dominate).
    if field == "big" and fast_big_flag:
        _cpmt.FAST_BIG = True
    x = (1 << n) - 1
    V = isqrt(x); nb = max(1, V.bit_length())
    primes, _, _ = compressed_lucy(x)
    K = len(primes)
    truth = sieve_pi(x)
    mode = "DELEGATED wiring (O(nb log p))" if delegate else "automaton wiring (O(nb p))"
    if delegate:
        mode += (", STRUCTURED prover O(nb p)" if structured
                 else ", dense prover O(nb p^2)")

    # ---- full compressed chain (the headline) ----
    res = run_chain(x, np.random.default_rng(seed), delegate=delegate,
                    structured=structured, q=q, pcs=pcs, commit_base=commit_base,
                    commit_code=commit_code)
    if pcs:
        mode += ", PCS leaf openings (residuals threaded to base)"
    if commit_base:
        mode += (f", S_0 base via tensor PCS ({commit_code} code, "
                 f"sub-sqrt-x verifier)")
    print(f"x = 2^{n}-1 = {x}, V = floor(sqrt x) = {V}, compressed cubes "
          f"2^{nb} = {1 << nb}, layers K = {K}")
    fld = ("  [uint64 Mersenne]" if fast_big(q)
           else "  [object dtype]" if q > Q else "")
    print(f"field = {field} (q = {q}{fld}); mode: {mode}")
    if q < x:
        print(f"  WARNING: q < x -- the trace certifier admits a wrap-around "
              f"alias at this field; use --field big for sound n>~30.")
    print(f"CHAIN honest run: {'ACCEPTED' if res['accepted'] else 'REJECTED'}; "
          f"claimed pi(x) = {res['claimed']}, sieve pi(x) = {truth}, "
          f"match = {res['claimed'] == truth}")
    assert res["accepted"] and res["claimed"] == truth
    print(f"  t_prover = {res['t_prover']*1000:.2f} ms   "
          f"t_verifier = {res['t_verifier']*1000:.3f} ms   "
          f"comm = {res['comm']} field elems (~{res['comm']*4/1024:.1f} KB)")
    rej = sum(not run_chain(x, np.random.default_rng(seed + 10 + t),
                            cheat="delta_pi", delegate=delegate,
                            structured=structured, q=q, pcs=pcs,
                            commit_base=commit_base, commit_code=commit_code)["accepted"]
              for t in range(cheat_trials))
    print(f"  chain cheat delta_pi (claim pi+1):           "
          f"rejected {rej}/{cheat_trials}")
    for i0 in sorted(set([K, max(1, K // 2), 1])):
        rej = sum(not run_chain(x, np.random.default_rng(seed + 100 + i0 + t),
                                corrupt_layer=i0, delegate=delegate,
                                structured=structured, q=q, pcs=pcs,
                                commit_base=commit_base, commit_code=commit_code)["accepted"]
                  for t in range(cheat_trials))
        print(f"  chain self-consistent liar @ layer {i0:>2}/{K}:        "
              f"rejected {rej}/{cheat_trials}")

    # ---- two-sided single layer cheat panel ----
    if layer is None:
        layer = max(1, K // 2)
    tl = run_two_sided_layer(x, layer, np.random.default_rng(seed),
                             delegate=delegate, structured=structured, q=q)
    print(f"two-sided single layer {layer}/{K} (prime p = {primes[layer-1]}): "
          f"{'ACCEPTED' if tl['accepted'] else 'REJECTED'}")
    assert tl["accepted"]
    cheats = ["wrong_C", "G", "split", "u", "waff",
              "small_C", "small_G", "small_u", "carry"]
    if delegate:                       # delegation-specific liars
        cheats += ["waff_forge", "waff_chain", "small_forge", "small_chain"]
    for ch in cheats:
        rej = sum(not run_two_sided_layer(x, layer,
                                          np.random.default_rng(seed + 30 + t),
                                          cheat=ch, delegate=delegate,
                                          structured=structured, q=q)["accepted"]
                  for t in range(cheat_trials))
        print(f"  layer cheat {ch:>11}: rejected {rej}/{cheat_trials}")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, default=12)
    ap.add_argument("--layer", type=int, default=None)
    ap.add_argument("--cheat-trials", type=int, default=20)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--delegate", action="store_true",
                    help="route both wirings through the O(nb log p) chain")
    ap.add_argument("--structured", action="store_true",
                    help="S496 O(nb p) chain prover (delegate only; "
                         "bit-identical transcript)")
    ap.add_argument("--field", choices=list(FIELDS), default="q",
                    help="prime field: q=2^31-1 (demo, uint64), big=2^61-1 "
                         "(uint64 Mersenne fast path, sound for n<~60), "
                         "small=1000003")
    ap.add_argument("--selftest", action="store_true")
    ap.add_argument("--bench", action="store_true")
    ap.add_argument("--bench-big", action="store_true",
                    help="BIG_Q chain: uint64 Mersenne vs object speedup")
    ap.add_argument("--bench-combined", action="store_true",
                    help="end-to-end run_chain over BIG_Q across batch_trace/"
                         "batch_wiring/FAST_BIG configs (S503 headline)")
    ap.add_argument("--fast-big", action="store_true",
                    help="opt into the uint64 Mersenne path for --field big "
                         "(pays off only at large n; slower for small chains)")
    ap.add_argument("--wiring-bench", action="store_true")
    ap.add_argument("--prover-bench", action="store_true")
    ap.add_argument("--bench-pcs", action="store_true",
                    help="leaf-opening stand-in vs sum-check openings (S505)")
    ap.add_argument("--bench-ub", action="store_true",
                    help="per-layer Ub openings inline vs ONE batched discharge "
                         "(S506): verifier O(2^nb) leaf-eval count K*nb -> 0")
    ap.add_argument("--bench-verifier-ops", action="store_true",
                    help="whole-chain verifier large-table-eval op-count CURVE "
                         "across no-pcs / pcs / pcs+batch_trace+batch_ub (S507): "
                         "leading term Theta(x) -> Theta(sqrt x)")
    ap.add_argument("--pcs", action="store_true",
                    help="use real sum-check leaf openings (S505), residuals "
                         "threaded to the base")
    ap.add_argument("--commit-base", action="store_true",
                    help="discharge the two one-time S_0 base closes via the "
                         "tensor polynomial commitment (S508): sub-sqrt-x verifier")
    ap.add_argument("--commit-code", choices=("rs", "expander"), default="rs",
                    help="pcs row code (S514): rs (Reed-Solomon, needs N<=q) or "
                         "expander (Brakedown-style, field-size-free)")
    args = ap.parse_args()
    if args.selftest:
        selftest()
    elif args.bench:
        bench()
    elif args.bench_pcs:
        bench_pcs()
    elif args.bench_ub:
        bench_ub()
    elif args.bench_verifier_ops:
        bench_verifier_ops()
    elif args.bench_big:
        bench_big()
    elif args.bench_combined:
        bench_combined(args.n)
    elif args.wiring_bench:
        wiring_bench()
    elif args.prover_bench:
        prover_bench()
    else:
        main(args.n, args.layer, args.cheat_trials, args.seed, args.delegate,
             args.structured, args.field, args.fast_big, args.pcs, args.commit_base,
             args.commit_code)
