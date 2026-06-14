#!/usr/bin/env python3
"""census_entropy_asymptotic.py — the census entropy CONSTANT without
exact full enumeration, via a mean-field + cluster (Mayer/cumulant)
expansion.

Open item 4, PROGRAM.md NEXT ACTION (the part that survives S519).

S519 proved the transfer-matrix route to the EXACT A_k is closed (state
space Theta(A_k)) and reduced A_k to a count over only the active primes
B(k) ~ rho*(k) ~ k/ln k:

    A_k = #{ S subset allowed(k) : for every prime q in B(k), the occupied
            residues {o mod q : o in S} do NOT cover all of Z_q }.

The remaining question is purely asymptotic: the conjecture
A_k = e^{(1+o(1)) pi(k)}, i.e. ln A_k / (k/ln k) -> 1.  This script
estimates that constant WITHOUT exact counting.

THE IDEA (exact reformulation).  Put the uniform measure on subsets:
each allowed offset is included independently with probability 1/2.  Then

    A_k = 2^N * P[ admissible ],   N = |allowed(k)| = phi(k),

EXACTLY (every subset has probability 2^{-N}).  Let C_q = "S covers all
classes mod q" (an INCREASING event) and \bar C_q = "S misses a class
mod q" (DECREASING).  admissible = intersection over q in B(k) of \bar C_q.

  * MARGINAL (independent within a prime).  Classes mod q partition the
    allowed offsets into disjoint groups of sizes n_{q,0..q-1}, so
        p_q := P[\bar C_q] = 1 - prod_j (1 - 2^{-n_{q,j}}),
    computed exactly from the class sizes (no transfer matrix).

  * MEAN FIELD = product over q (independent-primes leading term).  The
    \bar C_q are DECREASING events, hence positively associated (FKG), so
        P[ intersection ] >= prod_q p_q   ==>   A_k >= 2^N prod_q p_q
    is a RIGOROUS LOWER BOUND.

  * CLUSTER EXPANSION (the corrections = CRT overlaps).  For T subset
    B(k) let P_T = P[ intersection over q in T of \bar C_q ].  Mobius
    inversion of f(T) = ln P_T gives the connected cluster terms
        c_T = sum_{U subset T} (-1)^{|T|-|U|} ln P_U,
    and EXACTLY  ln P_B = sum_{T subset B, T != empty} c_T.
    Truncating at |T| <= r is the order-r cluster approximation (r=1 is
    the mean field).

THE ENGINE (joint cover via nested inclusion-exclusion).  P_T needs the
joint cover probabilities P[ intersection_{q in S} C_q ] for S subset T.
These have a NESTED I-E form: pick the largest prime L in S as the
"product" prime and inclusion-exclude over the missed classes of the
smaller "head" primes,

  P[cover all q in S] = sum_{(U_h subset Z_h)_h} (-1)^{sum|U_h|}
        2^{-(#offsets in any head-class U_h)} * prod_{w in Z_L}
        (1 - 2^{-n_w}),   n_w = #offsets with (mod L = w) and
        (mod h not in U_h for every head h).

Cost = 2^{sum of the all-but-largest primes of S} -- TINY for the small
primes that dominate the entropy, and INDEPENDENT of k (it depends only
on the primes, not on N).  This is exact (validated bit-for-bit against
the S519 transfer matrix count_tm) and FEASIBLE at k = 128, 256 where
count_tm on the full B(k) (Theta(A_k) states, S519) is dead.

WHY THIS BEATS THE FULL TRANSFER MATRIX.  count_tm on the full B(k) has
Theta(A_k) states.  A 3- or 4-prime cluster term needs only a 2^{O(small
primes)} I-E, so the r=3,4 cluster approximation is computable at k=128
(and k=256), validated against exact A_k for k <= 64.

MEASURED (validation, clean family):
  - the mean field is always below exact (FKG holds);
  - the cluster errors ALTERNATE in sign and shrink ~5x per order
    (odd r under, even r over), so [r=3, r=2] empirically BRACKETS A_k;
  - r=3 is within <1% of ln A_k by k=64.

What would falsify a result here:
  - the mean field exceeding exact A_k for any k (FKG/lower-bound bug);
  - the joint-cover I-E disagreeing with count_tm (selftest asserts
    bit-for-bit on pairs and triples over k=32,64);
  - the full cluster sum (r=|B|) NOT equalling ln A_k - N ln2 (Mobius
    identity bug) -- selftest asserts it;
  - the cluster truncation error NOT shrinking with r at validated k
    (then the expansion does not converge and the estimate is void --
    reported honestly, no extrapolation).

Usage:
  python3 census_entropy_asymptotic.py --selftest
  python3 census_entropy_asymptotic.py --k 64               # detail one k
  python3 census_entropy_asymptotic.py --k 128 --rmax 4     # estimate A_128
  python3 census_entropy_asymptotic.py --scan 8 256 --rmax 4
"""

import argparse
import math
import os
import sys
import time
from itertools import combinations, product

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from census_transfer_matrix import (  # noqa: E402
    allowed_offsets, enforceable_primes, reduce_primes, count_tm,
    count_admissible_dfs, KNOWN,
)

try:
    import numpy as np
    _HAVE_NP = hasattr(np, "bitwise_count")     # vectorised popcount (numpy>=2)
except Exception:                                # pragma: no cover
    np = None
    _HAVE_NP = False

# joint-cover engine selector: the numpy-vectorised path (default when
# available) is bit-for-bit equal to the pure-Python recursion (selftest [12])
# but ~50-100x faster, which is what makes the bracket sweep and the k=256
# single-block upper bound feasible.  --no-vec forces the reference path.
USE_VEC = _HAVE_NP

LN2 = math.log(2.0)
_MASK64 = (1 << 64) - 1


# ----------------------------------------------------------------------
# class sizes and the marginal p_q (exact, no transfer matrix)
# ----------------------------------------------------------------------
def class_sizes(allowed, q):
    per = [0] * q
    for o in allowed:
        per[o % q] += 1
    return per


def marginal_lnp(allowed, q):
    """ln p_q = ln P[S misses a class mod q], p=1/2.
    p_q = 1 - prod_j (1 - 2^{-n_j}); computed via log1p/expm1 so we never
    subtract two near-equal floats for tiny p_q."""
    ln_cover = 0.0
    for n in class_sizes(allowed, q):
        ln_cover += math.log1p(-(2.0 ** (-n)))
    pq = -math.expm1(ln_cover)
    if pq <= 0.0:                              # cover==1 to float precision
        pq = sum(2.0 ** (-n) for n in class_sizes(allowed, q))
    return math.log(pq), pq


# ----------------------------------------------------------------------
# joint cover probability via nested inclusion-exclusion (the engine)
# ----------------------------------------------------------------------
def joint_cover_prob(allowed, S, cost_cap=None, memo=None):
    """P[ intersection_{q in S} C_q ] = P[cover every class mod q, all q
    in S] under the uniform (p=1/2) measure.

    Nested I-E with the LARGEST prime as the product prime.  Cost =
    prod over head primes of 2^h = 2^{sum(S)-max(S)}.  Returns None if
    that exceeds cost_cap (caller treats the cluster as unavailable).
    S = [] -> 1.0 (empty cover).  `memo` (keyed by sorted-tuple) avoids
    recomputing the same joint cover (the costly large blocks recur in
    both the upper bound and the lower-bound partition).

    Dispatches to the numpy-vectorised engine when available (bit-for-bit
    equal, selftest [12]); `USE_VEC=False`/`--no-vec` forces the reference
    Python recursion."""
    if USE_VEC:
        return _jc_vec(allowed, S, cost_cap, memo)
    return _jc_py(allowed, S, cost_cap, memo)


def _jc_py(allowed, S, cost_cap=None, memo=None):
    """Reference pure-Python nested-I-E joint cover (the readable spec the
    vectorised engine is validated against)."""
    S = sorted(S)
    key = tuple(S)
    if memo is not None and key in memo:
        return memo[key]
    if not S:
        return 1.0
    heads, L = S[:-1], S[-1]
    leaves = 1
    for h in heads:
        leaves <<= h                            # 2^h
    if cost_cap is not None and leaves > cost_cap:
        if memo is not None:
            memo[key] = None
        return None
    no = len(allowed)
    # per-head per-residue offset bitmasks; per-L-class offset bitmasks
    head_cls = [[0] * h for h in heads]
    for i, o in enumerate(allowed):
        for hi, h in enumerate(heads):
            head_cls[hi][o % h] |= (1 << i)
    Lcls = [0] * L
    for i, o in enumerate(allowed):
        Lcls[o % L] |= (1 << i)

    total = 0.0
    # recurse over head subsets, accumulating the excluded-offset bitmask
    def rec(hi, excl_mask, sign):
        nonlocal total
        if hi == len(heads):
            excl = excl_mask.bit_count()
            prod = 2.0 ** (-excl)
            for w in range(L):
                n_w = (Lcls[w] & ~excl_mask).bit_count()
                prod *= (1.0 - 2.0 ** (-n_w))
            total += sign * prod
            return
        cls = head_cls[hi]
        h = heads[hi]
        for U in range(1 << h):
            m = 0
            uu = U
            r = 0
            while uu:
                if uu & 1:
                    m |= cls[r]
                uu >>= 1
                r += 1
            s = -sign if (U.bit_count() & 1) else sign
            rec(hi + 1, excl_mask | m, s)
    rec(0, 0, 1)
    if memo is not None:
        memo[key] = total
    return total


# ----- numpy-vectorised joint cover (the same nested I-E, batched) ----------
def _nlimbs(N):
    return (N + 63) // 64


def _int_to_limbs(m, nlimbs):
    """Python-int bitmask -> (nlimbs,) uint64 little-endian limb vector."""
    return np.array([(m >> (64 * j)) & _MASK64 for j in range(nlimbs)],
                    dtype=np.uint64)


def _class_int_masks(allowed, q):
    """Per-residue offset bitmask (Python int): bit i set iff allowed[i]%q==r."""
    cls = [0] * q
    for i, o in enumerate(allowed):
        cls[o % q] |= (1 << i)
    return cls


def _subset_or_masks(cls_ints, nlimbs):
    """Build the FULL inclusion-exclusion subset table over a list of
    per-class offset masks (Python ints), by numpy doubling: for each class
    mask c, [masks, masks|c] and [signs, -signs].  After processing all C
    classes the result is (limb masks (nlimbs, 2^C), signs (2^C,)) with
    mask[U] = OR of the selected classes, sign[U] = (-1)^|U|.  This is O(C)
    vectorised numpy ops -- crucially, NO Python loop over the 2^C subsets,
    so a single large head prime (e.g. 23) builds in milliseconds, not 8.4M
    Python iterations."""
    masks = np.zeros((nlimbs, 1), dtype=np.uint64)
    signs = np.array([1.0])
    for c in cls_ints:
        cl = _int_to_limbs(c, nlimbs)[:, None]          # (nlimbs,1)
        masks = np.concatenate([masks, masks | cl], axis=1)
        signs = np.concatenate([signs, -signs])
    return masks, signs


def _jc_vec(allowed, S, cost_cap=None, memo=None, batch_log2=20):
    """Vectorised nested-I-E joint cover.  Identical math to _jc_py: the
    largest prime L is the product prime, the smaller "head" primes are
    inclusion-excluded.  The largest heads (suffix with sum<=batch_log2) are
    batched into ONE numpy array of 2^{sum} excluded-offset masks; the rest
    recurse in Python.  Multi-limb uint64 masks support N>64 (k=256)."""
    S = sorted(S)
    key = tuple(S)
    if memo is not None and key in memo:
        return memo[key]
    if not S:
        return 1.0
    heads, L = S[:-1], S[-1]
    sum_heads = sum(heads)
    if cost_cap is not None and (1 << sum_heads) > cost_cap:
        if memo is not None:
            memo[key] = None
        return None
    N = len(allowed)
    nlimbs = _nlimbs(N)
    Lcls_col = np.stack([_int_to_limbs(m, nlimbs)
                         for m in _class_int_masks(allowed, L)], axis=1)  # (nlimbs,L)

    if not heads:                                # single-prime cover: no I-E
        prod = 1.0
        for w in range(L):
            nw = int(np.bitwise_count(Lcls_col[:, w]).sum())
            prod *= (1.0 - 2.0 ** (-nw))
        if memo is not None:
            memo[key] = prod
        return prod

    # split heads (sorted ascending) into a batched suffix (the LARGEST heads,
    # vectorised) and a recursed prefix (the smallest heads, Python loop).
    # The largest head is ALWAYS batched even if it alone exceeds batch_log2:
    # recursing it would mean 2^h Python iterations, while the doubling builder
    # handles it in O(h) numpy ops.  The Python recursion then runs only over
    # the small remaining heads (2^{sum of those}, tiny).
    batched_heads = []
    s = 0
    for h in reversed(heads):                    # largest first
        if not batched_heads or s + h <= batch_log2:
            batched_heads.append(h)
            s += h
        else:
            break
    batched_set = set(batched_heads)
    outer = [h for h in heads if h not in batched_set]

    batch_cls = []
    for h in batched_heads:
        batch_cls.extend(_class_int_masks(allowed, h))
    bm, bs = _subset_or_masks(batch_cls, nlimbs)

    outer_cls = [_class_int_masks(allowed, h) for h in outer]
    total = 0.0

    def leaf(base_mask, base_sign):
        nonlocal total
        full = bm | base_mask[:, None]                  # (nlimbs, M)
        excl = np.bitwise_count(full).sum(axis=0)       # (M,)
        prodvec = np.exp2(-excl.astype(np.float64))
        notfull = ~full
        for w in range(L):
            aw = Lcls_col[:, w][:, None] & notfull
            nw = np.bitwise_count(aw).sum(axis=0)
            prodvec *= (1.0 - np.exp2(-nw.astype(np.float64)))
        total += base_sign * float(np.dot(bs, prodvec))

    def rec(hi, base_mask, base_sign):
        if hi == len(outer):
            leaf(base_mask, base_sign)
            return
        cls, h = outer_cls[hi], outer[hi]
        for U in range(1 << h):
            m, uu, r = 0, U, 0
            while uu:
                if uu & 1:
                    m |= cls[r]
                uu >>= 1
                r += 1
            sgn = -base_sign if (U.bit_count() & 1) else base_sign
            rec(hi + 1, base_mask | _int_to_limbs(m, nlimbs), sgn)

    rec(0, np.zeros(nlimbs, dtype=np.uint64), 1.0)
    if memo is not None:
        memo[key] = total
    return total


# ----- EXACT integer joint cover (no floating-point cancellation) -----------
# jc(S) = P[cover every class mod q, all q in S] is a DYADIC rational m/2^N:
# each nested-I-E term equals sign * prod_w (2^{n_w}-1) / 2^N exactly (because
# 2^{-excl} * prod_w 2^{n_w}(1-2^{-n_w}) and excl + sum_w n_w = N).  So jc*2^N
# is an EXACT integer -- joint_cover_int returns it.  This matters at large N
# (k>=256): there P[miss] ~ 2^{-Theta(N)} ~ 1e-17, while the I-E sums terms ~1,
# so the FLOAT engine cancels to ZERO significant digits (selftest [14]); the
# integer engine has none of that error (bounded only by the I-E cost cap).
def joint_cover_int(allowed, S, cost_cap=None, memo=None):
    """jc(S) * 2^N as an exact Python integer (jc = this / 2^len(allowed))."""
    S = sorted(S)
    key = tuple(S)
    if memo is not None and key in memo:
        return memo[key]
    N = len(allowed)
    if not S:
        v = 1 << N                                   # jc=1 -> num=2^N
        if memo is not None:
            memo[key] = v
        return v
    heads, L = S[:-1], S[-1]
    if cost_cap is not None and (1 << sum(heads)) > cost_cap:
        if memo is not None:
            memo[key] = None
        return None
    head_cls = [[0] * h for h in heads]
    for i, o in enumerate(allowed):
        for hi, h in enumerate(heads):
            head_cls[hi][o % h] |= (1 << i)
    Lcls = [0] * L
    for i, o in enumerate(allowed):
        Lcls[o % L] |= (1 << i)

    total = 0

    def rec(hi, excl_mask, sign):
        nonlocal total
        if hi == len(heads):
            prod = 1
            for w in range(L):
                n_w = (Lcls[w] & ~excl_mask).bit_count()
                prod *= ((1 << n_w) - 1)
            total += sign * prod
            return
        cls, h = head_cls[hi], heads[hi]
        for U in range(1 << h):
            m, uu, r = 0, U, 0
            while uu:
                if uu & 1:
                    m |= cls[r]
                uu >>= 1
                r += 1
            s = -sign if (U.bit_count() & 1) else sign
            rec(hi + 1, excl_mask | m, s)
    rec(0, 0, 1)
    if memo is not None:
        memo[key] = total
    return total


def lnP_miss_exact(allowed, T, cost_cap=None, memo=None):
    """ln P[miss a class mod every q in T], computed in EXACT integers (no
    cancellation at any N).  None if any joint-cover term exceeds cost_cap.
    `memo` here keys integer numerators (keep it separate from the float
    joint-cover memo)."""
    N = len(allowed)
    tot = 0
    for r in range(len(T) + 1):
        for S in combinations(T, r):
            num = joint_cover_int(allowed, list(S), cost_cap, memo)
            if num is None:
                return None
            tot += ((-1) ** r) * num
    if tot <= 0:                                     # P[miss]<=0 is impossible
        return None
    return math.log(tot) - N * LN2                   # math.log handles big int


def P_miss(allowed, T, cost_cap=None, memo=None, noise_rtol=1e-11):
    """P[ intersection_{q in T} \bar C_q ] (miss a class mod every q in T)
    = sum_{S subset T} (-1)^{|S|} P[cover all q in S], FLOAT engine.  None if
    any needed joint-cover term exceeds cost_cap.

    NOISE GUARD: the I-E sums terms ~1 to reach P[miss] ~ 2^{-Theta(N)}; once
    |sum| falls below noise_rtol * max|term| the float result is dominated by
    rounding error (this happens at N>~100, k>=256).  We then return None
    rather than a garbage value -- callers should use lnP_miss_exact there.
    Set noise_rtol=0 to disable (e.g. when validating against the exact engine
    at small N)."""
    tot = 0.0
    maxabs = 0.0
    for r in range(len(T) + 1):
        for S in combinations(T, r):
            jc = joint_cover_prob(allowed, list(S), cost_cap, memo)
            if jc is None:
                return None
            tot += ((-1) ** r) * jc
            maxabs = max(maxabs, abs(jc))
    if noise_rtol and (tot <= 0.0 or tot < noise_rtol * maxabs):
        return None                                  # below the float noise floor
    return tot


# ----------------------------------------------------------------------
# cluster expansion
# ----------------------------------------------------------------------
def ie_cost_log2(block):
    """log2 of the nested-I-E cost of a block = sum of all but the largest."""
    if len(block) <= 1:
        return 0
    return sum(block) - max(block)


def greedy_blocks(B, cost_cap):
    """Partition the active primes into consecutive blocks, each with I-E
    cost <= cost_cap, greedily packing the SMALLEST primes together (they
    carry the strongest correlations).  Coarser blocks => tighter FKG
    lower bound."""
    blocks, cur = [], []
    for q in B:
        if not cur:
            cur = [q]
            continue
        if (1 << ie_cost_log2(cur + [q])) <= cost_cap:
            cur.append(q)
        else:
            blocks.append(cur)
            cur = [q]
    if cur:
        blocks.append(cur)
    return blocks


def block_lower_bound(allowed, blocks, cost_cap, memo=None):
    """ln of the FKG block lower bound: A_k >= 2^N * prod_blocks P[miss
    within block].  Each block's joint-miss prob is exact (nested I-E).
    Returns (ln_lower, ok)."""
    N = len(allowed)
    tot = N * LN2
    for g in blocks:
        pm = P_miss(allowed, list(g), cost_cap, memo)
        if pm is None or pm <= 0.0:
            return None, False
        tot += math.log(pm)
    return tot, True


def best_single_block(B, cost_cap):
    """The largest prefix of B (smallest primes) whose single-block I-E is
    within the cap -- gives the tightest feasible UPPER bound
    2^N * P[miss within that block] >= A_k."""
    g = []
    for q in B:
        if (1 << ie_cost_log2(g + [q])) <= cost_cap:
            g.append(q)
        else:
            break
    return g


def cluster_term(T, lnP):
    """c_T = sum_{U subset T} (-1)^{|T|-|U|} ln P_U.  None if any ln P_U
    unavailable."""
    Tl = list(T)
    s = 0.0
    for r in range(len(Tl) + 1):
        for U in combinations(Tl, r):
            v = lnP.get(frozenset(U))
            if v is None:
                return None
            s += ((-1) ** (len(Tl) - r)) * v
    return s


def estimate(k, rmax=4, core_cut=None, cost_cap=1 << 20, max_states=6_000_000,
             want_exact=True):
    """Mean-field lower bound + cluster truncations r=1..rmax + rigorous
    count_tm upper bound + (if feasible) exact ln A_k.

    core_cut: only cluster over active primes <= core_cut (larger active
    primes stay at mean-field order; justified because p_q ~ 1 there)."""
    allowed = allowed_offsets(k)
    enf = enforceable_primes(k, allowed)
    B, W, dropped = reduce_primes(k, allowed, enf)
    N = len(allowed)
    base = N * LN2

    lnp, pq = {}, {}
    for q in B:
        lnp[q], pq[q] = marginal_lnp(allowed, q)
    mf_all = sum(lnp[q] for q in B)
    meanfield = base + mf_all                       # all-singletons FKG lower

    memo = {}                                       # shared joint-cover cache

    # rigorous FKG block lower bound (coarser partition => tighter)
    blocks = greedy_blocks(B, cost_cap)
    block_lower, _ = block_lower_bound(allowed, blocks, cost_cap, memo)

    # rigorous single-block upper bound: 2^N P[miss within G] >= A_k
    Gup = best_single_block(B, cost_cap)
    pm_up = P_miss(allowed, Gup, cost_cap, memo)
    upper_ie = base + math.log(pm_up) if (pm_up and pm_up > 0) else None

    core = list(B) if core_cut is None else [q for q in B if q <= core_cut]

    # ln P_T for every subset of core with 1 <= |T| <= rmax (singletons are
    # the exact marginals; |T|>=2 via nested I-E)
    lnP = {frozenset(): 0.0}
    for q in core:
        lnP[frozenset((q,))] = lnp[q]
    skipped = []
    for r in range(2, min(rmax, len(core)) + 1):
        for T in combinations(core, r):
            pm = P_miss(allowed, list(T), cost_cap, memo)
            if pm is None or pm <= 0.0:
                lnP[frozenset(T)] = None
                skipped.append(tuple(T))
            else:
                lnP[frozenset(T)] = math.log(pm)

    # cluster truncations
    trunc = {}
    cluster_contrib = {1: mf_all}
    for r in range(2, rmax + 1):
        if r > len(core):
            break
        s = 0.0
        for T in combinations(core, r):
            cT = cluster_term(T, lnP)
            if cT is not None:
                s += cT
        cluster_contrib[r] = s
    for rr in range(1, rmax + 1):
        tot = base
        for r in range(1, rr + 1):
            tot += cluster_contrib.get(r, 0.0)
        trunc[rr] = tot

    # exact ln A_k via the S519 transfer matrix on the full active set, but
    # ONLY when it is feasible (small k) -- prefixes are nested so we climb
    # and stop at the first abort (never pays the slow full-B abort).
    exact = None
    if want_exact and N <= 40:                 # k <= 80: full B count_tm fits
        last = None
        for i in range(1, len(B) + 1):
            c, st = count_tm(k, list(B[:i]), allowed, max_states=max_states,
                             want_stats=True)
            if c is None:
                last = None
                break
            last = (math.log(c), B[:i])
        if last is not None and last[1] == B:
            exact = last[0]

    return {
        "k": k, "N": N, "B": B, "W": W, "core": core, "base": base,
        "meanfield": meanfield, "block_lower": block_lower, "blocks": blocks,
        "upper": upper_ie, "upper_B": Gup, "trunc": trunc,
        "cluster_contrib": cluster_contrib, "exact": exact, "lnp": lnp,
        "pq": pq, "skipped": skipped,
    }


def ratio(lnval, k):
    return lnval / (k / math.log(k))


# ----------------------------------------------------------------------
# reporting
# ----------------------------------------------------------------------
def report_one(k, rmax, core_cut, cost_cap, max_states, verbose=False):
    t0 = time.time()
    r = estimate(k, rmax, core_cut, cost_cap, max_states, want_exact=True)
    dt = time.time() - t0
    kk = k / math.log(k)
    print(f"\nk={k}  N={r['N']}  active B={r['B']} (W={r['W']})  "
          f"core={r['core']}  k/lnk={kk:.3f}  [{dt:.1f}s]")
    if verbose:
        print("   prime   class sizes        p_q        ln p_q")
        allowed = allowed_offsets(k)
        for q in r["B"]:
            sz = class_sizes(allowed, q)
            print(f"   {q:5d}   n in [{min(sz)},{max(sz)}]      "
                  f"{r['pq'][q]:.4e}  {r['lnp'][q]:9.4f}")
    ex = r["exact"]
    if ex is not None:
        print(f"  EXACT  ln A_k        = {ex:8.4f}   ratio = {ratio(ex,k):.4f}")
    print("  --- rigorous bracket (FKG block lower, single-block I-E upper) ---")
    print(f"  UPPER  block {r['upper_B']} = {r['upper']:8.4f}   "
          f"ratio = {ratio(r['upper'], k):.4f}")
    bl = r["block_lower"]
    if bl is not None:
        print(f"  LOWER  blocks {r['blocks']} = {bl:8.4f}   "
              f"ratio = {ratio(bl, k):.4f}")
    print(f"  (mean field, all-singleton FKG lower = {r['meanfield']:8.4f}   "
          f"ratio = {ratio(r['meanfield'], k):.4f})")
    print("  --- cluster-expansion point estimates (not one-sided) ---")
    for rr in sorted(r["trunc"]):
        if rr == 1:
            continue
        v = r["trunc"][rr]
        err = f"   err {v-ex:+.4f}" if ex is not None else ""
        print(f"  cluster r<={rr}         = {v:8.4f}   "
              f"ratio = {ratio(v, k):.4f}{err}")
    if r["skipped"]:
        print(f"  ({len(r['skipped'])} clusters skipped over cost cap "
              f"2^{int(math.log2(cost_cap))}: e.g. {r['skipped'][:4]})")
    return r


def _ls_slope(pts):
    n = len(pts)
    sx = sum(x for x, _ in pts)
    sy = sum(y for _, y in pts)
    sxx = sum(x * x for x, _ in pts)
    sxy = sum(x * y for x, y in pts)
    return (n * sxy - sx * sy) / (n * sxx - sx * sx)


def report_scan(kmin, kmax, rmax, core_cut, cost_cap, max_states):
    ks = [1 << m for m in range(2, 20) if kmin <= (1 << m) <= kmax]
    print(f"{'k':>6} {'N':>4} {'|B|':>3} {'exact':>9} {'LOWER':>9} "
          f"{'UPPER':>9} {'clus_r4':>9} | {'rt_ex':>7} {'rt_lo':>7} "
          f"{'rt_up':>7} {'rt_cl':>7}")
    recs = []
    for k in ks:
        r = estimate(k, rmax, core_cut, cost_cap, max_states, want_exact=True)
        recs.append(r)

        def fmt(v, w=9):
            return f"{v:{w}.3f}" if v is not None else f"{'--':>{w}}"

        def rt(v, w=7):
            return f"{ratio(v,k):{w}.4f}" if v is not None else f"{'--':>{w}}"
        ex, lo, up = r["exact"], r["block_lower"], r["upper"]
        cl = r["trunc"].get(min(rmax, 4)) or r["trunc"].get(rmax)
        print(f"{k:>6} {r['N']:>4} {len(r['B']):>3} {fmt(ex)} {fmt(lo)} "
              f"{fmt(up)} {fmt(cl)} | {rt(ex)} {rt(lo)} {rt(up)} {rt(cl)}")

    print("\nclean family k=2^m  (slope of ln(.) vs k/ln k -> conjecture 1.0):")
    for name, getter in (("exact", lambda r: r["exact"]),
                         ("FKG block UPPER", lambda r: r["upper"]),
                         ("FKG block LOWER", lambda r: r["block_lower"]),
                         ("mean field (singletons)",
                          lambda r: r["meanfield"])):
        pts = [(r["k"] / math.log(r["k"]), getter(r))
               for r in recs if getter(r) is not None]
        if len(pts) >= 2:
            sl = _ls_slope(pts)
            print(f"  {name:24s}: slope = {sl:.4f}   "
                  f"last ratio = {pts[-1][1] / pts[-1][0]:.4f}")

    print("\nRIGOROUS bracket on the entropy ratio ln A_k / (k/ln k):")
    for r in recs:
        lo, up = r["block_lower"], r["upper"]
        ex = r["exact"]
        tag = f"  (exact {ratio(ex,r['k']):.4f})" if ex is not None else \
              "  [exact unreachable -- bracket is the result]"
        if lo is not None and up is not None:
            print(f"  k={r['k']:>4}:  {ratio(lo,r['k']):.4f}  <=  "
                  f"ratio  <=  {ratio(up,r['k']):.4f}{tag}")


# ----------------------------------------------------------------------
# bracket WIDTH vs cost budget -- does the FKG-block method close to a point?
# ----------------------------------------------------------------------
def prefix_thresholds(B):
    """The discrete cost thresholds (log2) at which the single-block UPPER
    bound admits one more prime: cost of each prefix of B.  Between two
    consecutive thresholds the best feasible upper block -- and hence the
    upper bound -- is CONSTANT (the plateau)."""
    out = []
    for j in range(1, len(B) + 1):
        out.append((B[:j], ie_cost_log2(B[:j])))
    return out


def bracket_sweep(k, caps_log2, max_states=6_000_000, exact_ie=None):
    """For a fixed k, sweep the I-E cost budget C and record the best
    FEASIBLE rigorous bracket at each C:
      UPPER = base + ln P[miss within the largest affordable prefix block],
      LOWER = base + Sum_blocks ln P[miss within block] (greedy coarsest
              partition affordable at C).
    Returns the rows plus the threshold structure and the closing cost
    (ie_cost_log2(B), the single full-B block that would shut the bracket).

    exact_ie: use the EXACT integer joint-cover (no float cancellation).
    Default = auto (exact when N>96, i.e. k>=256, where the float engine
    underflows; float otherwise, where it is validated and far faster)."""
    allowed = allowed_offsets(k)
    enf = enforceable_primes(k, allowed)
    B, W, _ = reduce_primes(k, allowed, enf)
    N = len(allowed)
    base = N * LN2
    if exact_ie is None:
        exact_ie = N > 96
    memo = {}

    def lnpm(T, cap):
        if exact_ie:
            return lnP_miss_exact(allowed, list(T), cap, memo)
        pm = P_miss(allowed, list(T), cap, memo)
        return math.log(pm) if (pm is not None and pm > 0) else None

    exact = None
    if N <= 40:                                  # k<=80: full B count_tm fits
        c = count_tm(k, list(B), allowed, max_states=max_states)
        if c is not None:
            exact = math.log(c)

    rows = []
    for C in caps_log2:
        cap = 1 << C
        Gup = best_single_block(B, cap)
        lu = lnpm(Gup, cap)
        upper = base + lu if lu is not None else None
        blocks = greedy_blocks(B, cap)
        lower = base
        for g in blocks:
            lg = lnpm(g, cap)
            if lg is None:
                lower = None
                break
            lower += lg
        ru = ratio(upper, k) if upper is not None else None
        rl = ratio(lower, k) if lower is not None else None
        width = (ru - rl) if (ru is not None and rl is not None) else None
        rows.append({"C": C, "Gup": Gup, "blocks": blocks, "ru": ru,
                     "rl": rl, "width": width})
    return {"k": k, "N": N, "B": B, "W": W, "exact": exact, "exact_ie": exact_ie,
            "rows": rows, "full_cost": ie_cost_log2(B),
            "thresholds": prefix_thresholds(B)}


def report_bracket_sweep(k, caps_log2, max_states, exact_ie=None):
    t0 = time.time()
    r = bracket_sweep(k, caps_log2, max_states, exact_ie=exact_ie)
    dt = time.time() - t0
    kk = k / math.log(k)
    ex = r["exact"]
    eng = "EXACT integer" if r["exact_ie"] else "float I-E"
    print(f"\n=== bracket-width sweep, k={k}  N={r['N']}  "
          f"active B={r['B']} (|B|={len(r['B'])}, W={r['W']})  "
          f"k/lnk={kk:.3f}  engine={eng}  [{dt:.1f}s] ===")
    if ex is not None:
        print(f"  exact ln A_k = {ex:.4f}   ratio = {ratio(ex,k):.4f}")
    print("  upper-block prefix thresholds (cost to admit one more prime):")
    feas = []
    for blk, c in r["thresholds"]:
        feas.append(c)
        print(f"    block {blk}  ->  cost 2^{c}")
    print(f"  CLOSING cost (single full-B block => bracket shuts) = "
          f"2^{r['full_cost']}")
    print(f"\n  {'C(=log2 cap)':>12} {'upper-block':>26} {'ratio_up':>9} "
          f"{'#blk_lo':>8} {'ratio_lo':>9} {'WIDTH':>8}")
    prev_w = None
    for row in r["rows"]:
        gup = "{" + ",".join(map(str, row["Gup"])) + "}"
        nlo = len(row["blocks"])
        ru = f"{row['ru']:.4f}" if row["ru"] is not None else "--"
        rl = f"{row['rl']:.4f}" if row["rl"] is not None else "--"
        w = row["width"]
        flat = ""
        if w is not None and prev_w is not None and abs(w - prev_w) < 1e-9:
            flat = "  (plateau)"
        if w is not None:
            prev_w = w
        ws = f"{w:.4f}" if w is not None else "--"
        print(f"  {row['C']:>12} {gup:>26} {ru:>9} {nlo:>8} {rl:>9} "
              f"{ws:>8}{flat}")
    # plateau summary
    widths = [row["width"] for row in r["rows"] if row["width"] is not None]
    if widths:
        print(f"\n  width over the swept range: max {max(widths):.4f} -> "
              f"min {min(widths):.4f}")
        # the smallest cap at which the upper block is already the max
        # feasible prefix => any larger feasible C leaves it unchanged
        last_gup = r["rows"][-1]["Gup"]
        nxt = len(last_gup) + 1
        if nxt <= len(r["B"]):
            nxt_cost = r["thresholds"][nxt - 1][1]
            print(f"  NEXT upper improvement needs cost 2^{nxt_cost} "
                  f"(adds prime {r['B'][nxt-1]}); "
                  f"bracket fully shuts only at 2^{r['full_cost']}.")
    return r


# ----------------------------------------------------------------------
# selftest
# ----------------------------------------------------------------------
def selftest():
    print("=== census_entropy_asymptotic selftest ===")
    n_ok = 0

    # 1. marginal p_q (product) == count_tm({q})/2^N
    for k in (16, 32, 64):
        allowed = allowed_offsets(k)
        N = len(allowed)
        for q in enforceable_primes(k, allowed)[:4]:
            lp, p = marginal_lnp(allowed, q)
            lp_tm = math.log(count_tm(k, [q], allowed)) - N * LN2
            assert abs(lp - lp_tm) < 1e-9, (k, q, lp, lp_tm)
    n_ok += 1
    print("[1] marginal p_q (class-size product) == count_tm({q})/2^N")

    # 2. joint cover I-E == count_tm joint miss, bit-for-bit, pairs+triples
    for k in (32, 64):
        allowed = allowed_offsets(k)
        N = len(allowed)
        B, _, _ = reduce_primes(k, allowed)
        tests = [B[:2], B[:3]]
        if len(B) >= 5:
            tests += [list(B[2:5]), list(B[1:4])]
        for T in tests:
            ie = P_miss(allowed, list(T))
            tm = count_tm(k, list(T), allowed) / (2.0 ** N)
            assert abs(ie - tm) <= 1e-12 + 1e-9 * abs(tm), (k, T, ie, tm)
    n_ok += 1
    print("[2] joint-cover nested I-E == count_tm (pairs+triples, k=32,64)")

    # 3. Mobius identity: full cluster sum (r=|B|) == ln A_k - N ln2.
    #    P_T from count_tm (feasible for any subset at k<=64; the full-B
    #    nested-I-E would be 2^{sum of all but max} -- too slow -- which is
    #    exactly why we TRUNCATE; the IE engine itself is checked in [2]).
    for k in (8, 16, 32, 64):
        allowed = allowed_offsets(k)
        B, _, _ = reduce_primes(k, allowed)
        N = len(allowed)
        lnP = {frozenset(): 0.0}
        for r in range(1, len(B) + 1):
            for T in combinations(B, r):
                lnP[frozenset(T)] = math.log(count_tm(k, list(T), allowed)) \
                    - N * LN2
        total = 0.0
        for r in range(1, len(B) + 1):
            for T in combinations(B, r):
                total += cluster_term(T, lnP)
        lnA = math.log(count_tm(k, list(B), allowed))
        assert abs((N * LN2 + total) - lnA) < 1e-7, (k, N * LN2 + total, lnA)
    n_ok += 1
    print("[3] full cluster sum (r=|B|) == ln A_k - N ln2 (Mobius exact)")

    # 4. mean field <= exact A_k (FKG lower bound), k=8..64
    for k in (8, 16, 32, 48, 64):
        r = estimate(k, rmax=1, want_exact=True)
        assert r["exact"] is not None
        assert r["meanfield"] <= r["exact"] + 1e-9, \
            (k, r["meanfield"], r["exact"])
    n_ok += 1
    print("[4] mean field <= exact ln A_k (FKG lower bound), k=8..64")

    # 5. rigorous upper bound: count_tm(prefix) >= A_k, monotone
    for k in (16, 32, 64):
        allowed = allowed_offsets(k)
        B, _, _ = reduce_primes(k, allowed)
        A = count_tm(k, list(B), allowed)
        prev = None
        for i in range(1, len(B) + 1):
            up = count_tm(k, list(B[:i]), allowed)
            assert up >= A, (k, B[:i], up, A)
            if prev is not None:
                assert up <= prev, (k, i, up, prev)   # monotone non-increasing
            prev = up
    n_ok += 1
    print("[5] count_tm(prefix of B) >= A_k, monotone (rigorous upper bound)")

    # 6. cluster error shrinks with r at validated k; r<=3 within 1% at k=64.
    #    (rmax=3: every cluster of order <=3 is below the default cost cap at
    #    k<=64, so the truncations are complete -- no skips.)
    for k, rm in ((32, 3), (64, 3)):
        r = estimate(k, rmax=rm, want_exact=True)
        ex = r["exact"]
        errs = [abs(r["trunc"][rr] - ex) for rr in sorted(r["trunc"])]
        for a, b in zip(errs, errs[1:]):
            assert b <= a + 1e-9, (k, errs)
    r64 = estimate(64, rmax=3, want_exact=True)
    assert abs(r64["trunc"][3] - r64["exact"]) / abs(r64["exact"]) < 0.01, \
        (r64["trunc"][3], r64["exact"])
    n_ok += 1
    print("[6] cluster error shrinks with r; r<=3 within 1% at k=64")

    # 7. alternating sign of the truncation error (odd r under, even over)
    #    at validated k=32,64 for r=1,2,3 -- the empirical-bracket basis
    for k in (32, 64):
        r = estimate(k, rmax=3, want_exact=True)
        ex = r["exact"]
        for rr in (1, 2, 3):
            if rr >= len(r["B"]):
                continue                       # exact at/after full order
            e = r["trunc"][rr] - ex
            if rr % 2 == 1:
                assert e <= 1e-9, (k, rr, e)   # odd: under
            else:
                assert e >= -1e-9, (k, rr, e)  # even: over
    n_ok += 1
    print("[7] truncation error alternates: odd r under, even r over (k=32,64)")

    # 8. exact field == known A_k / DFS
    for k in (8, 16, 32, 64):
        allowed = allowed_offsets(k)
        B, _, _ = reduce_primes(k, allowed)
        assert count_tm(k, list(B), allowed) == KNOWN[k], k
    assert abs(math.exp(estimate(16, rmax=2)["exact"])
               - count_admissible_dfs(16)) < 0.5
    n_ok += 1
    print("[8] count_tm(B) == known A_k {8,16,32,64}; exact == DFS at k=16")

    # 9. cost cap behaves: a cluster over the cap is skipped (None), and an
    #    estimate at a lower cap is still produced (graceful degradation)
    allowed = allowed_offsets(128)
    assert joint_cover_prob(allowed, [11, 13, 17], cost_cap=1 << 18) is None
    assert joint_cover_prob(allowed, [3, 5, 7], cost_cap=1 << 18) is not None
    rr = estimate(128, rmax=3, cost_cap=1 << 16, want_exact=False)
    assert rr["meanfield"] < rr["trunc"][2]            # mean field below r=2
    assert rr["trunc"].get(2) is not None
    n_ok += 1
    print("[9] cost cap skips over-budget clusters; estimate still produced")

    # 10. the rigorous bracket BRACKETS exact, and the block lower bound is a
    #     genuine improvement over the all-singletons mean field, k=8..64
    for k in (16, 32, 48, 64):
        r = estimate(k, rmax=3, want_exact=True)
        ex = r["exact"]
        assert ex is not None
        assert r["block_lower"] <= ex + 1e-9, (k, "block_lower>exact",
                                               r["block_lower"], ex)
        assert r["upper"] >= ex - 1e-9, (k, "upper<exact", r["upper"], ex)
        assert r["meanfield"] <= r["block_lower"] + 1e-9, \
            (k, "mf>block_lower")
    n_ok += 1
    print("[10] FKG bracket: block_lower <= exact <= upper; "
          "mean field <= block_lower (k=16..64)")

    # 11. coarsening the partition monotonically tightens the FKG lower bound:
    #     a single block {3,5,7} >= product of {3,5}*{7} >= product of
    #     {3}*{5}*{7} (all <= exact), at k=64
    allowed = allowed_offsets(64)
    N = len(allowed)
    one = N * LN2 + math.log(P_miss(allowed, [3, 5, 7]))
    two = N * LN2 + math.log(P_miss(allowed, [3, 5])) + \
        math.log(P_miss(allowed, [7]))
    three = N * LN2 + sum(marginal_lnp(allowed, q)[0] for q in (3, 5, 7))
    ex357 = math.log(count_tm(64, [3, 5, 7], allowed))
    assert three <= two + 1e-9 <= one + 1e-9 <= ex357 + 1e-9, \
        (three, two, one, ex357)
    n_ok += 1
    print("[11] coarser FKG partition tightens the lower bound monotonically")

    # 12. the numpy-vectorised engine == the pure-Python reference (bit-for-bit
    #     to float tolerance), over multi-limb N (k=128 N=64; k=256 N=128),
    #     singletons / pairs / triples and the unaffordable-cap None path.
    if _HAVE_NP:
        for k in (32, 64, 128, 256):
            allowed = allowed_offsets(k)
            B, _, _ = reduce_primes(k, allowed)
            tests = [B[:1], B[:2], B[:3], B[1:3]]
            if len(B) >= 5:
                tests += [B[2:5], [B[0], B[4]]]
            for S in tests:
                a = _jc_py(allowed, S, cost_cap=1 << 20)
                b = _jc_vec(allowed, S, cost_cap=1 << 20)
                if a is None or b is None:
                    assert (a is None) == (b is None), (k, S, a, b)
                    continue
                rel = abs(a - b) / (abs(a) + 1e-300)
                assert rel < 1e-10, (k, S, a, b, rel)
            # the cost cap returns None identically on both engines
            assert (_jc_py(allowed, B[:4], cost_cap=1 << 4) is None
                    and _jc_vec(allowed, B[:4], cost_cap=1 << 4) is None), k
        # vectorised P_miss still equals count_tm (independent oracle), k=128
        allowed = allowed_offsets(128)
        Np = len(allowed)
        ie = P_miss(allowed, [3, 5, 7])
        tm = count_tm(128, [3, 5, 7], allowed) / (2.0 ** Np)
        assert abs(ie - tm) <= 1e-12 + 1e-9 * abs(tm), (ie, tm)
        n_ok += 1
        print("[12] numpy-vectorised joint cover == pure-Python reference "
              "(k=32..256, multi-limb) and == count_tm")
    else:                                            # pragma: no cover
        print("[12] SKIPPED (numpy>=2 with bitwise_count not available)")

    # 13. bracket sweep: (a) WIDTH non-increasing in the cost budget and the
    #     upper block CONSTANT (plateau) between prefix-cost thresholds (k=64,
    #     caps <= 2^20 so the expensive 2^26 block is never built); (b) at a
    #     cap >= the closing cost the bracket shuts to the exact value (k=32,
    #     closing cost 2^8, cheap).
    sw = bracket_sweep(64, [6, 9, 16, 20])           # thresholds 0,3,8,15,(26)
    ws = [row["width"] for row in sw["rows"] if row["width"] is not None]
    for a, b in zip(ws, ws[1:]):
        assert b <= a + 1e-9, ("width not monotone", ws)
    g16 = next(r["Gup"] for r in sw["rows"] if r["C"] == 16)
    g20 = next(r["Gup"] for r in sw["rows"] if r["C"] == 20)
    assert g16 == g20 == [3, 5, 7, 11], ("plateau broken", g16, g20)
    assert sw["full_cost"] == 26, sw["full_cost"]    # k=64 closes only at 2^26
    sw32 = bracket_sweep(32, [2, 4, 9])              # closing cost 8
    assert sw32["full_cost"] == 8, sw32["full_cost"]
    r9 = next(r for r in sw32["rows"] if r["C"] == 9)
    assert abs(r9["width"]) < 1e-6, ("bracket should shut at closing cost",
                                     r9["width"])
    assert sw32["exact"] is not None
    assert abs((sw["exact"] / (64 / math.log(64))) - 0.9276) < 0.01
    n_ok += 1
    print("[13] bracket width monotone in cost; upper block constant between "
          "thresholds (k=64); bracket shuts to exact at closing cost (k=32)")

    # 14. EXACT integer joint cover: (a) jc_int/2^N == float at small N;
    #     (b) the I-E numerator == count_tm integer exactly (no cancellation);
    #     (c) the FLOAT PRECISION WALL -- at k=256 (N=128) the raw float P_miss
    #     of a small-prime block UNDERFLOWS to <=0 while the exact engine gives
    #     a positive ~1e-19, and the guarded P_miss returns None (refuses).
    for k in (64, 128):
        allowed = allowed_offsets(k)
        N = len(allowed)
        for S in ([3, 5], [3, 5, 7], [5, 7, 11], [7, 11, 13]):
            ni = joint_cover_int(allowed, S)
            fi = _jc_py(allowed, S)
            assert abs(ni / (2.0 ** N) - fi) / (abs(fi) + 1e-300) < 1e-9, (k, S)
    for k in (64, 128, 256):                          # exact == count_tm integer
        allowed = allowed_offsets(k)
        for T in ([3, 5], [3, 5, 7]):
            num = 0
            for r in range(len(T) + 1):
                for S in combinations(T, r):
                    num += ((-1) ** r) * joint_cover_int(allowed, list(S))
            assert num == count_tm(k, T, allowed, max_states=8_000_000), (k, T)
    allowed = allowed_offsets(256)                    # the precision wall
    lp = lnP_miss_exact(allowed, [3, 5, 7])
    assert lp is not None and lp < -40, lp           # exact ~ ln(1.48e-19)
    assert math.exp(lp) > 0                           # genuinely positive
    raw = P_miss(allowed, [3, 5, 7], noise_rtol=0.0)  # float, guard OFF
    assert raw <= 1e-30, ("float should underflow at N=128", raw)
    assert P_miss(allowed, [3, 5, 7]) is None         # guard ON -> refuses
    n_ok += 1
    print("[14] exact integer joint cover == float (small N) == count_tm; "
          "float underflows at k=256 (N=128), exact + guard handle it")

    print(f"\nALL SELFTESTS PASSED ({n_ok} groups)")


# ----------------------------------------------------------------------
def main():
    global USE_VEC
    ap = argparse.ArgumentParser()
    ap.add_argument("--selftest", action="store_true")
    ap.add_argument("--k", type=int, default=0)
    ap.add_argument("--scan", type=int, nargs=2, default=None,
                    metavar=("KMIN", "KMAX"))
    ap.add_argument("--bracket-sweep", type=int, default=0, metavar="K",
                    help="sweep the I-E cost budget at this k and report the "
                         "rigorous bracket width vs cost (the method-closure "
                         "measurement)")
    ap.add_argument("--caps", type=int, nargs="+", default=None,
                    help="explicit list of log2 cost caps for --bracket-sweep")
    ap.add_argument("--rmax", type=int, default=4)
    ap.add_argument("--core-cut", type=int, default=None)
    ap.add_argument("--cost-cap", type=int, default=20,
                    help="log2 of the joint-cover I-E cost cap (default 20)")
    ap.add_argument("--max-states", type=int, default=6_000_000)
    ap.add_argument("--no-vec", action="store_true",
                    help="force the pure-Python joint-cover engine")
    ap.add_argument("--exact-ie", dest="exact_ie", action="store_true",
                    default=None,
                    help="force the EXACT integer joint cover in the bracket "
                         "sweep (no float cancellation; slower)")
    ap.add_argument("--float-ie", dest="exact_ie", action="store_false",
                    help="force the float joint cover in the bracket sweep")
    ap.add_argument("--verbose", action="store_true")
    args = ap.parse_args()

    if args.no_vec:
        USE_VEC = False

    cap = 1 << args.cost_cap
    if args.selftest:
        selftest()
    elif args.bracket_sweep:
        caps = args.caps if args.caps else list(range(8, args.cost_cap + 1, 2))
        report_bracket_sweep(args.bracket_sweep, caps, args.max_states,
                             exact_ie=args.exact_ie)
    elif args.k:
        report_one(args.k, args.rmax, args.core_cut, cap, args.max_states,
                   args.verbose)
    elif args.scan:
        report_scan(args.scan[0], args.scan[1], args.rmax, args.core_cut,
                    cap, args.max_states)
    else:
        selftest()


if __name__ == "__main__":
    main()
