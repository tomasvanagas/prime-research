#!/usr/bin/env python3
"""
dk_sliding_census.py — S523: the SLIDING-WINDOW local pattern census
D_k^slide(x), the unaligned companion of the aligned census
(local_pattern_census.py / dk_first_occurrence.py).

Aligned census (prior work):   windows [k*m, k*m+k), base n = k*m ≡ 0
(mod k).  Admissibility for the occupancy pattern S ⊆ {0..k-1} has TWO
conditions:
  (i)  q | k  =>  no o in S with q | o   (k*m+o ≡ o, a fixed residue),
  (ii) q ∤ k  =>  {o mod q : o in S} ≠ Z_q.
A_k^aligned = 13, 106, 3573 for k = 8, 16, 32.

Sliding census (THIS script):  windows [n, n+k) with n ranging over ALL
residues.  Because n is now unconstrained mod every prime, condition (i)
DISAPPEARS and admissibility collapses to the CLASSICAL Hardy-Littlewood
prime-constellation criterion on the offset set S:

    S ⊆ {0..k-1} is admissible  <=>  for EVERY prime q <= k,
        {o mod q : o in S} ≠ Z_q    (the offsets miss a class mod q).

The q = 2 clause forces every admissible S to be all-even or all-odd
(consecutive primes share parity), so for k a power of two A_k^slide
splits cleanly into an even-offset family and an odd-offset family.

Structural facts this script establishes (selftest-asserted):
  * A_k^slide = 2*A_k^aligned - 1   (k a power of two): 25, 211, 7145.
    The odd-offset admissible family is EXACTLY the aligned admissible
    set; the even-offset family is its mirror under the reflection
    o ↦ k-1-o (an affine bijection mod every q, so it preserves residue
    coverage ⇒ preserves admissibility), and the two families share only
    the empty pattern.  Hence #even = #odd = A_k^aligned, A_slide =
    2*A_aligned - 1.
  * Max admissible weight (densest constellation) = ρ*(k) = 9 at k=32,
    the SAME H-geometry as aligned — parity-symmetric, so sliding does
    not reach a denser tuple than the aligned census.
  * Aligned == sliding restricted to n ≡ 0 (mod k)  (selftest, the base
    case the NEXT ACTION asked for): an aligned window is literally a
    sliding window anchored at a multiple of k.

Census law (falsifiable): with sliding windows every classical
admissible pattern occurs (prime k-tuple conjecture), and a finite set
of small-n windows passing over the small primes {2,3,5,7,...} realize
inadmissible patterns, so

    D_k^slide(x) -> A_k^slide + E_k^slide   (E_k^slide = small-n
                                             exceptional patterns, finite)

We measure E_k^slide and the approach D_k^slide(x) -> A_k^slide+E with a
constant-memory segmented sieve, and test first occurrences against the
CLASSICAL HL singular series N_S(x) = S(S) * x / (ln x)^w (no /k, no q|k
factors — n runs over all integers, ~x windows up to x not x/k).

What would falsify:
  - an inadmissible pattern realized at LARGE n (impossible by plain
    divisibility — a bug check; none expected beyond the small-n band);
  - an admissible pattern still missing at x far beyond its HL mean
    count N_S(x) >= 5 while denser cousins appear (breaks HL ordering);
  - D_k^slide(x) plateauing strictly below A_k^slide + E_k^slide.

Usage:
  python3 dk_sliding_census.py --k 32 --xmax $((1<<28))
  python3 dk_sliding_census.py --k 8  --xmax $((1<<20)) --validate-hl
  python3 dk_sliding_census.py --selftest
"""

import argparse
import math
from collections import Counter

import numpy as np


# ----------------------------------------------------------------------
# primes / classical admissibility
# ----------------------------------------------------------------------
def primes_upto(limit):
    """All primes <= limit as a python list."""
    if limit < 2:
        return []
    s = np.ones(limit + 1, dtype=bool)
    s[:2] = False
    for p in range(2, int(limit ** 0.5) + 1):
        if s[p]:
            s[p * p::p] = False
    return [int(p) for p in np.nonzero(s)[0]]


def admissible_codes_sliding(k):
    """All CLASSICAL-admissible offset sets S ⊆ {0..k-1}, returned as
    bitmask codes.  DFS over offsets with the standard downward-closed
    pruning: including an offset that completes a full residue class
    mod some prime q kills the whole supertree.  Node count
    ~ O(A_k^slide * k) — the q=2 prune (mixed parity dead) keeps it
    tiny."""
    primes = primes_upto(k)
    full = {q: (1 << q) - 1 for q in primes}
    obits = [{q: 1 << (o % q) for q in primes} for o in range(k)]
    out = []
    # iterative DFS: (offset index, per-prime coverage masks, set bitmask)
    stack = [(0, tuple(0 for _ in primes), 0)]
    while stack:
        idx, cov, S = stack.pop()
        if idx == k:
            out.append(S)
            continue
        # exclude offset idx
        stack.append((idx + 1, cov, S))
        # include offset idx (prune if it completes a class mod some q)
        newcov = []
        dead = False
        for qi, q in enumerate(primes):
            c = cov[qi] | obits[idx][q]
            if c == full[q]:
                dead = True
                break
            newcov.append(c)
        if not dead:
            stack.append((idx + 1, tuple(newcov), S | (1 << idx)))
    return out


def code_of(S, k):
    c = 0
    for o in S:
        c |= 1 << o
    return c


def set_of(code, k):
    return frozenset(o for o in range(k) if (code >> o) & 1)


def parity_split(adm_codes, k):
    """(#even-only, #odd-only) admissible patterns (empty set counts in
    neither bucket here; it is the single shared pattern)."""
    even_mask = sum(1 << o for o in range(0, k, 2))
    odd_mask = sum(1 << o for o in range(1, k, 2))
    neven = sum(1 for c in adm_codes if c != 0 and (c & odd_mask) == 0)
    nodd = sum(1 for c in adm_codes if c != 0 and (c & even_mask) == 0)
    return neven, nodd


# ----------------------------------------------------------------------
# segmented sliding census
# ----------------------------------------------------------------------
def _seg_bounds(xmax, seg_size):
    """Breakpoints on the window-START axis [0, xmax] that include every
    power of two <= xmax (so each power-of-two checkpoint lands on a
    segment boundary) and every multiple of seg_size."""
    brk = {0, xmax}
    p = 1
    while p <= xmax:
        brk.add(p)
        p <<= 1
    m = seg_size
    while m < xmax:
        brk.add(m)
        m += seg_size
    return sorted(brk)


def iter_sliding_codes(k, xmax, seg_size, nmod=None):
    """Yield (n_start, codes) over consecutive window-start segments,
    codes[i] = k-bit pattern of window [n_start+i, n_start+i+k):
    bit o set iff (n_start+i+o) is prime.  Window starts n in [0, xmax).
    Constant memory ~ O(seg_size).  Sieves k-1 positions past each
    segment for the tail windows.

    nmod = (div, rem): keep only window starts n ≡ rem (mod div)
    (used to recover the aligned census as n ≡ 0 mod k)."""
    sieve_hi = xmax + k                      # positions [0, sieve_hi)
    base = primes_upto(int(sieve_hi ** 0.5) + 1)
    shifts = np.arange(k, dtype=np.int64)
    bounds = _seg_bounds(xmax, seg_size)
    for lo, hi in zip(bounds[:-1], bounds[1:]):
        n_count = hi - lo                    # window starts lo..hi-1
        if n_count <= 0:
            continue
        L = n_count + (k - 1)                # positions [lo, lo+L)
        seg = np.ones(L, dtype=bool)
        seg_hi = lo + L
        for p in base:
            start = max(p * p, ((lo + p - 1) // p) * p)
            if start >= seg_hi:
                continue
            seg[start - lo::p] = False
        # 0 and 1 are not prime but have no prime factor, so the sieve
        # leaves them True — mask them explicitly wherever they land
        # (a power-of-two boundary can start a segment AT position 1).
        for pos in (0, 1):
            if lo <= pos < seg_hi:
                seg[pos - lo] = False
        codes = np.zeros(n_count, dtype=np.int64)
        for o in range(k):
            col = seg[o:o + n_count]
            if col.any():
                codes |= col.astype(np.int64) << shifts[o]
        if nmod is not None:
            div, rem = nmod
            # absolute n = lo + i ; keep n % div == rem
            i0 = (rem - lo) % div
            sel = np.arange(i0, n_count, div, dtype=np.int64)
            yield lo, codes[sel] if sel.size else codes[:0]
        else:
            yield lo, codes


def census(k, xmax, seg_size=1 << 24, want_counts=False, progress=False,
           nmod=None):
    """Segmented sliding census. Returns:
      first[code]  -> first-occurrence window start n
      checkpoints  -> [(x, D, n_adm_found, n_missing, n_exceptional)]
      counts[code] -> #occurrences up to xmax (only if want_counts)
      adm_codes, A_k
    Checkpoints are taken at powers of two of the window-start bound x
    (D counts windows with start n < x)."""
    adm_codes = set(admissible_codes_sliding(k))
    A_k = len(adm_codes)

    first = {}
    counts = Counter() if want_counts else None
    checkpoints = []
    pow2 = set()
    p = 1
    while p <= xmax:
        pow2.add(p)
        p <<= 1

    for n_start, codes in iter_sliding_codes(k, xmax, seg_size, nmod=nmod):
        if codes.size:
            uniq, idx = np.unique(codes, return_index=True)
            for c, i in zip(uniq.tolist(), idx.tolist()):
                if c not in first:
                    first[c] = n_start + i
            if want_counts:
                uq, cc = np.unique(codes, return_counts=True)
                for c, n in zip(uq.tolist(), cc.tolist()):
                    counts[c] += n

    # Checkpoints at powers of two of the window-start bound x: a pattern
    # counts toward D(x) iff its first occurrence n < x.  Done in one
    # cheap pass over the (small) first-occurrence list.
    occ = sorted((n, c) for c, n in first.items())
    for x in sorted(pow2 | {xmax}):
        seen = {c for n, c in occ if n < x}
        found = len(adm_codes & seen)
        exc = len(seen - adm_codes)
        checkpoints.append((x, len(seen), found, A_k - found, exc))
        if progress:
            lx = int(round(math.log2(x))) if x > 0 else 0
            print(f"  x=2^{lx:<2} D={len(seen):>6} adm={found:>5}/{A_k} "
                  f"missing={A_k - found:>5} exc={exc:>4}", flush=True)

    return {"first": first, "checkpoints": checkpoints, "counts": counts,
            "adm_codes": adm_codes, "A_k": A_k}


# ----------------------------------------------------------------------
# classical Hardy-Littlewood singular series + first-occurrence estimate
# ----------------------------------------------------------------------
def singular_series(S, pbound=20000):
    """Classical HL singular series for the constellation S (offsets):
        S(S) = prod_q (1 - omega(q)/q) / (1 - 1/q)^w ,
    omega(q) = #{o mod q : o in S}.  Inadmissible (omega=q) -> 0."""
    w = len(S)
    if w == 0:
        return 1.0
    prod = 1.0
    for q in primes_upto(pbound):
        omega = len({o % q for o in S})
        if omega == q:
            return 0.0
        prod *= (1.0 - omega / q) / ((1.0 - 1.0 / q) ** w)
    return prod


def hl_count(S, x, pbound=20000, Sa=None):
    """Expected #sliding windows [n,n+k), n <= x, realizing the (>= S)
    constellation:  N_S(x) = S(S) * x / (ln x)^w  (no /k: n runs over
    all integers)."""
    w = len(S)
    if w == 0:
        return float(x)
    if Sa is None:
        Sa = singular_series(S, pbound)
    if Sa == 0.0:
        return 0.0
    lx = math.log(x)
    return Sa * x / (lx ** w)


def predict_first_occ(S, pbound=20000, Sa=None):
    """Smallest x with N_S(x) >= 1 (HL first-occurrence estimate)."""
    w = len(S)
    if w == 0:
        return 2.0
    if Sa is None:
        Sa = singular_series(S, pbound)
    if Sa == 0.0:
        return math.inf
    lo, hi = 8.0, 1.0e6
    while hl_count(S, hi, Sa=Sa) < 1.0 and hi < 1e305:
        hi *= 10.0
    if hl_count(S, hi, Sa=Sa) < 1.0:
        return math.inf
    for _ in range(300):
        mid = math.sqrt(lo * hi)
        if hl_count(S, mid, Sa=Sa) < 1.0:
            lo = mid
        else:
            hi = mid
    return hi


# ----------------------------------------------------------------------
# reporting
# ----------------------------------------------------------------------
def rho_star(k):
    """Max admissible weight (densest classical constellation in an
    interval of length k) = ρ*(k)."""
    adm = admissible_codes_sliding(k)
    return max(bin(c).count("1") for c in adm)


def report(k, xmax, res, pbound=20000, baseline_pow=0):
    first = res["first"]
    adm_codes = res["adm_codes"]
    A_k = res["A_k"]
    seen = set(first)
    found = adm_codes & seen
    missing = adm_codes - seen
    exceptional = seen - adm_codes

    neven, nodd = parity_split(adm_codes, k)  # both EXCLUDE empty
    A_aligned = nodd + 1  # odd-offset family + empty == aligned set (k=2^j)
    rk = max(bin(c).count("1") for c in adm_codes)

    print(f"\n=== D_{k}^slide(x) census, window starts n < "
          f"2^{int(round(math.log2(xmax)))} ({xmax}) ===")
    print(f"A_{k}^slide = {A_k}   ({neven} even-only + {nodd} odd-only "
          f"+ 1 shared empty)")
    print(f"  odd-offset family + empty == aligned admissible set: "
          f"A_{k}^aligned = {A_aligned};  A_slide = 2*A_aligned - 1 = "
          f"{2 * A_aligned - 1}")
    print(f"  densest admissible weight (ρ*({k})) = {rk}  "
          f"[parity-symmetric: even & odd both reach it]")
    print(f"  target census height: A_slide + E (E = small-n "
          f"exceptional inadmissible patterns, finite)")

    print(f"\n{'x':>6} {'D_slide':>8} {'adm':>6} {'missing':>8} "
          f"{'exc':>5} {'D/(A+E)':>9}")
    A_plus_E = A_k + len(exceptional)
    for x, D, fnd, miss, exc in res["checkpoints"]:
        lx = int(round(math.log2(x))) if x > 0 else 0
        denom = A_plus_E if A_plus_E else 1
        print(f"2^{lx:<4} {D:>8} {fnd:>6} {miss:>8} {exc:>5} "
              f"{D / denom:>9.4f}")

    # exceptional (inadmissible) patterns from the small-n band
    print(f"\nexceptional (inadmissible) patterns realized at small n: "
          f"{len(exceptional)}")
    exc_sorted = sorted(exceptional, key=lambda c: first[c])
    for c in exc_sorted[:20]:
        n = first[c]
        print(f"  n={n:>6} window [{n},{n + k}): offsets {sorted(set_of(c, k))}")
    if len(exc_sorted) > 20:
        print(f"  ... ({len(exc_sorted) - 20} more)")
    if exceptional:
        max_exc_n = max(first[c] for c in exceptional)
        print(f"  last exceptional first-occurrence at n={max_exc_n} "
              f"(beyond this, no inadmissible pattern should appear)")

    # newest dense first occurrences (admissible, weight-sorted high)
    baseline_x = 1 << baseline_pow if baseline_pow else 0
    new = sorted(((first[c], c) for c in found if first[c] >= baseline_x),
                 key=lambda t: (-len(set_of(t[1], k)), t[0]))
    print(f"\ndensest admissible first occurrences (top by weight):")
    print(f"  {'w':>2} {'first-occ n':>12} {'HL x*':>12} {'obs/x*':>8} "
          f"{'N(occ)':>8}  offsets")
    for n, c in new[:25]:
        S = set_of(c, k)
        nocc = max(n, 2)
        Sa = singular_series(S, pbound)
        xpred = predict_first_occ(S, pbound, Sa=Sa)
        ratio = nocc / xpred if math.isfinite(xpred) and xpred else \
            float('nan')
        N = hl_count(S, nocc, Sa=Sa)
        print(f"  {len(S):>2} {nocc:>12} {('%.3e' % xpred):>12} "
              f"{ratio:>8.2f} {N:>8.2f}  {sorted(S)}")

    # still-missing: HL falsification monitor (mean count, not x*<xmax)
    print(f"\nstill MISSING at x=2^{int(round(math.log2(xmax)))}: "
          f"{len(missing)} (by weight: "
          f"{dict(sorted(Counter(len(set_of(c, k)) for c in missing).items()))})")
    miss_info = []
    for c in missing:
        S = set_of(c, k)
        Sa = singular_series(S, pbound)
        N = hl_count(S, xmax, Sa=Sa)
        xstar = predict_first_occ(S, pbound, Sa=Sa)
        miss_info.append((N, len(S), xstar, c))
    miss_info.sort(reverse=True)
    if miss_info:
        print("  most-expected still-missing (largest HL mean N_S(xmax) "
              "first; P(0 occ)=e^{-N}):")
        print(f"  {'w':>2} {'N_S(xmax)':>10} {'P(0 occ)':>9} "
              f"{'HL x*':>12}  offsets")
        for N, w, xstar, c in miss_info[:12]:
            print(f"  {w:>2} {N:>10.3f} {math.exp(-N):>9.4f} "
                  f"{('%.3e' % xstar):>12}  {sorted(set_of(c, k))}")
    THRESH = 5.0
    overdue = [(N, w, c) for N, w, xstar, c in miss_info if N >= THRESH]
    if overdue:
        print(f"\n  [!] {len(overdue)} missing pattern(s) have HL mean "
              f"N_S(xmax) >= {THRESH:.0f} (P(0 occ) < "
              f"{math.exp(-THRESH):.2%}) — STATISTICAL ANOMALY vs HL:")
        for N, w, c in sorted(overdue, reverse=True)[:10]:
            print(f"      w={w} N={N:.2f} offsets {sorted(set_of(c, k))}")
    else:
        mx = miss_info[0][0] if miss_info else 0.0
        print(f"\n  census OK: every still-missing pattern has HL mean "
              f"N_S(xmax) < {THRESH:.0f} (largest {mx:.2f}, "
              f"P(0 occ)={math.exp(-mx):.2%}) — all in the expected HL "
              f"first-occurrence tail.")
    if not missing:
        print(f"\n  FULL SATURATION: D_{k}^slide = A_slide + E = "
              f"{A_k} + {len(exceptional)} = {A_plus_E}, every classical "
              f"admissible pattern realized.")


# ----------------------------------------------------------------------
# HL count validation (occurrence counts vs prediction, per weight)
# ----------------------------------------------------------------------
def validate_hl(k, xmax, res, pbound=20000):
    counts = res["counts"]
    adm_codes = res["adm_codes"]
    print(f"\n=== HL count validation, k={k}, x up to "
          f"2^{int(round(math.log2(xmax)))} ===")
    print(f"{'w':>2} {'#patterns':>10} {'obs_total':>14} "
          f"{'pred_total':>14} {'obs/pred':>9}")
    byw = {}
    for c in adm_codes:
        if c not in counts:
            continue
        S = set_of(c, k)
        w = len(S)
        obs = counts[c]
        pred = hl_count(S, xmax, pbound)
        d = byw.setdefault(w, [0, 0.0, 0])
        d[0] += obs
        d[1] += pred
        d[2] += 1
    for w in sorted(byw):
        obs, pred, n = byw[w]
        r = obs / pred if pred else float('nan')
        print(f"{w:>2} {n:>10} {obs:>14} {pred:>14.1f} {r:>9.3f}")
    print("(obs/pred ~ 1 validates the classical singular series + "
          "normalization; high-w drifts low — HL counts the >=S "
          "superset and dense patterns are pre-asymptotic.)")


# ----------------------------------------------------------------------
# references for selftest
# ----------------------------------------------------------------------
def _sieve_bool(n):
    s = np.ones(n, dtype=bool)
    s[:2] = False
    for p in range(2, int(n ** 0.5) + 1):
        if s[p]:
            s[p * p::p] = False
    return s


def _direct_sliding(k, x):
    """Reference: distinct sliding patterns for window starts n in [0,x),
    via a single full sieve (independent of the segmented path)."""
    s = _sieve_bool(x + k)
    codes = np.zeros(x, dtype=np.int64)
    for o in range(k):
        codes |= s[o:o + x].astype(np.int64) << o
    return set(int(c) for c in np.unique(codes))


def _naive_sliding(k, x):
    """Trial-division reference (validates the bit-packing itself)."""
    hi = x + k
    isp = [n >= 2 and all(n % d for d in range(2, int(n ** 0.5) + 1))
           for n in range(hi)]
    pats = set()
    for n in range(x):
        c = 0
        for o in range(k):
            if isp[n + o]:
                c |= 1 << o
        pats.add(c)
    return pats


def _direct_aligned(k, x):
    """Reference aligned realized set (windows [k*m,k*m+k), k*m<x)."""
    s = _sieve_bool(x)
    nwin = x // k
    win = s[:nwin * k].reshape(nwin, k)
    codes = win.astype(np.int64).dot(1 << np.arange(k, dtype=np.int64))
    return set(int(c) for c in np.unique(codes))


def _aligned_admissible_codes(k):
    """Aligned admissible patterns (the prior-work definition) as codes,
    for the cross-check that odd-offset sliding == aligned."""
    div_primes = [q for q in primes_upto(k) if k % q == 0]
    other_primes = [q for q in primes_upto(k) if k % q != 0]
    allowed = [o for o in range(k) if all(o % q != 0 for q in div_primes)]
    out = set()
    for mask in range(1 << len(allowed)):
        S = [allowed[i] for i in range(len(allowed)) if (mask >> i) & 1]
        ok = all(len({o % q for o in S}) != q for q in other_primes)
        if ok:
            out.add(code_of(S, k))
    return out


# ----------------------------------------------------------------------
# selftest
# ----------------------------------------------------------------------
def selftest():
    print("=== selftest ===")

    # 1. classical admissibility counts + the 2*A_aligned-1 identity
    counts = {}
    for k in (8, 16, 32):
        adm = set(admissible_codes_sliding(k))
        counts[k] = len(adm)
        neven, nodd = parity_split(adm, k)
        # odd-offset sliding admissible == aligned admissible set
        assert {c for c in adm if c != 0 and
                (c & sum(1 << o for o in range(0, k, 2))) == 0} == \
            (_aligned_admissible_codes(k) - {0}), k
        a_aligned = len(_aligned_admissible_codes(k))   # includes empty
        assert nodd == a_aligned - 1, (k, nodd, a_aligned)  # excl. empty
        assert neven == nodd, (k, neven, nodd)        # mirror symmetry
        assert len(adm) == 2 * a_aligned - 1, (k, len(adm), a_aligned)
    assert counts == {8: 25, 16: 211, 32: 7145}, counts
    print(f"[1] A_slide = {counts} ; even==odd==A_aligned-1; "
          f"A_slide = 2*A_aligned-1 OK")

    # 1b. mirror symmetry o->k-1-o bijects even<->odd admissible
    for k in (8, 16):
        adm = set(admissible_codes_sliding(k))
        def mirror(c):
            return code_of({k - 1 - o for o in set_of(c, k)}, k)
        assert {mirror(c) for c in adm} == adm, k
    print("[2] reflection o->k-1-o is an admissibility automorphism "
          "(even<->odd) OK")

    # 2. code/set round-trip
    for S in [frozenset(), frozenset({0, 2, 6}), frozenset({1, 3, 31})]:
        assert set_of(code_of(S, 32), 32) == S
    print("[3] code<->set round-trip OK")

    # 3. segmented == direct full sieve == naive (the bit logic)
    assert _naive_sliding(6, 200) == _direct_sliding(6, 200)
    for k, x in [(8, 1 << 14), (16, 1 << 15), (32, 1 << 16)]:
        direct = _direct_sliding(k, x)
        for seg in (k * 4, 1 << 12, 1 << 15):
            r = census(k, x, seg_size=max(seg, k))
            seen = set(int(c) for c in r["first"])
            assert seen == direct, (k, x, seg, len(seen), len(direct))
    print("[4] segmented census == direct full sieve == naive "
          "(k=6,8,16,32; 3 seg sizes) OK")

    # 4. THE BASE CASE: sliding restricted to n ≡ 0 (mod k) == aligned
    for k, x in [(8, 1 << 16), (16, 1 << 18), (32, 1 << 20)]:
        r = census(k, x, seg_size=max(1 << 15, k), nmod=(k, 0))
        seen = set(int(c) for c in r["first"])
        assert seen == _direct_aligned(k, x), (k, x, len(seen))
    print("[5] sliding census | n≡0 mod k  ==  aligned census "
          "(k=8,16,32) OK")

    # 5. saturation at k=8,16 (full reach), and exceptional count finite
    r8 = census(8, 1 << 18)
    D8, fnd8, miss8, exc8 = (r8["checkpoints"][-1][1],
                             r8["checkpoints"][-1][2],
                             r8["checkpoints"][-1][3],
                             r8["checkpoints"][-1][4])
    assert fnd8 == 25 and miss8 == 0, (D8, fnd8, miss8)   # all admissible
    assert D8 == 25 + exc8, (D8, exc8)
    r16 = census(16, 1 << 20)
    fnd16, miss16, exc16 = (r16["checkpoints"][-1][2],
                            r16["checkpoints"][-1][3],
                            r16["checkpoints"][-1][4])
    assert fnd16 == 211 and miss16 == 0, (fnd16, miss16)
    print(f"[6] saturation: D_8^slide={D8} (25 adm + {exc8} exc) by 2^18, "
          f"D_16^slide all 211 adm found by 2^20 ({exc16} exc) OK")

    # 5b. exceptional patterns: genuinely inadmissible, only at small n,
    #     and the count E_k^slide = 4,6,8 for k=8,16,32 — the windows
    #     n=0..q*(k) over the small-prime cluster (q*=3,5,7), realized by
    #     2^16 already (no new inadmissible pattern appears past n=q*).
    r32e = census(32, 1 << 16, seg_size=1 << 15)
    exc_n = {8: (r8, 3), 16: (r16, 5), 32: (r32e, 7)}
    E_expect = {8: 4, 16: 6, 32: 8}
    for k, (r, nmax_exc) in exc_n.items():
        exc = set(r["first"]) - r["adm_codes"]
        assert len(exc) == E_expect[k], (k, len(exc))
        for c in exc:
            assert singular_series(set_of(c, k)) == 0.0, (k, c)  # inadm
        assert max(r["first"][c] for c in exc) == nmax_exc, k    # n=0..q*
        assert min(r["first"][c] for c in exc) == 0, k
    print(f"[7] exceptional patterns inadmissible (S=0), E_k^slide="
          f"{{8:4,16:6,32:8}}=q*(k)+1, all at n=0..q*(3,5,7) OK")

    # 6. monotonic non-decreasing D
    for r in (r8, r16):
        ds = [c[1] for c in r["checkpoints"]]
        assert all(b >= a for a, b in zip(ds, ds[1:])), ds
    print("[8] D_k^slide(x) monotone non-decreasing OK")

    # 7. classical singular series: admissible>0, inadmissible=0,
    #    weight-1 == 1 (NOT 2 — no q|k aligned factor), twin {0,2}>1
    assert singular_series(frozenset({1})) == 1.0
    assert singular_series(frozenset({0})) == 1.0          # parity-free
    assert singular_series(frozenset({0, 1})) == 0.0       # covers mod 2
    assert singular_series(frozenset({0, 2, 4})) == 0.0    # covers mod 3
    tw = singular_series(frozenset({0, 2}))
    assert abs(tw - 1.3203) < 5e-3, tw                     # 2*C2 twin const
    # even/odd mirror patterns share the SAME singular series
    assert abs(singular_series(frozenset({0, 2, 6})) -
               singular_series(frozenset({1, 5, 7})) ) < 1e-12  # mirror in k=8
    print(f"[9] classical S(S): weight-1==1, mod-2/mod-3 cover==0, "
          f"twin {{0,2}}=={tw:.4f} (2C2), mirror-invariant OK")

    # 8. HL first-occ finite for all admissible; per-weight median grows
    adm32 = [set_of(c, 32) for c in admissible_codes_sliding(32)]
    by_w = {}
    for S in adm32:
        if S:
            by_w.setdefault(len(S), []).append(S)
    med = {}
    for w, lst in by_w.items():
        ps = sorted(predict_first_occ(S) for S in lst)
        assert all(math.isfinite(p) for p in ps), w
        med[w] = ps[len(ps) // 2]
    ws = sorted(med)
    # non-decreasing throughout; strictly increasing for the non-floored
    # weights (low weights all occur immediately -> floor at the search
    # lower bound, so they tie there — that is the honest behaviour).
    assert all(med[a] <= med[b] for a, b in zip(ws, ws[1:])), med
    hi_ws = [w for w in ws if w >= 5]
    assert all(med[a] < med[b] for a, b in zip(hi_ws, hi_ws[1:])), med
    print(f"[10] HL first-occ finite for all admissible; per-weight "
          f"median non-decreasing, strictly increasing for w>=5 "
          f"(w5 {med[5]:.1e} -> w9 {med[9]:.1e}) OK")

    # 9. HL count monotone in x; max weight == rho*(k)
    S = frozenset({0, 2, 6, 8, 12})
    assert hl_count(S, 1 << 26) > hl_count(S, 1 << 20)
    assert rho_star(8) == 3 and rho_star(16) == 5 and rho_star(32) == 9
    print("[11] HL count monotone in x; ρ*(8,16,32)=(3,5,9) OK")

    print("\nall selftests passed.")


# ----------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--k", type=int, default=32)
    ap.add_argument("--xmax", type=int, default=1 << 28,
                    help="upper bound on window start n (exclusive)")
    ap.add_argument("--seg", type=int, default=1 << 24,
                    help="segment size in window starts (memory ~ seg)")
    ap.add_argument("--baseline-pow", type=int, default=0,
                    help="2^p: only first-occurrences at n>=2^p are 'new'")
    ap.add_argument("--pbound", type=int, default=20000,
                    help="prime bound for the singular-series product")
    ap.add_argument("--aligned", action="store_true",
                    help="restrict to n≡0 mod k (recovers aligned census)")
    ap.add_argument("--validate-hl", action="store_true")
    ap.add_argument("--selftest", action="store_true")
    args = ap.parse_args()

    if args.selftest:
        selftest()
        return

    nmod = (args.k, 0) if args.aligned else None
    seg = max(args.seg, args.k)
    res = census(args.k, args.xmax, seg_size=seg,
                 want_counts=args.validate_hl, progress=True, nmod=nmod)
    report(args.k, args.xmax, res, pbound=args.pbound,
           baseline_pow=args.baseline_pow)
    if args.validate_hl:
        validate_hl(args.k, args.xmax, res, pbound=args.pbound)


if __name__ == "__main__":
    main()
