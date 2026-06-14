#!/usr/bin/env python3
"""
dk_first_occurrence.py — S522: the D_k(x) census first-occurrence probe.

Background (see local_pattern_census.py / _results.md): for aligned
k-windows of the prime indicator, D_k(x) = #distinct exact occupancy
patterns realized in windows [k*m, k*m+k) with k*m < x. The saturation
law is

    D_k(x) -> A_k + 1   (A_k = #Hardy-Littlewood admissible patterns,
                          +1 = the m=0 small-prime window),

proven (= divisibility) and verified saturated at k=8 (D=14, by 2^16)
and k=16 (D=107, by 2^18). At k=32 the census is still UNSATURATED at
x=2^26: D_32(2^26)=3385, with **189 admissible patterns unrealized**,
every one of weight >= 6 (23 of wt 6, 115 of 7, 47 of 8, 4 of 9). The
deficit D_32 = A_32 + 1 - D_32(x) is exactly the first-occurrence tail
of the densest admissible constellations.

This script extends the census past 2^26 with a SEGMENTED sieve
(constant memory, CLI --xmax) to:
  (i)  find the first occurrences of the 189 unrealized dense 32-tuples
       and measure how D_32(x) approaches A_32 + 1;
  (ii) compare each first occurrence to the Hardy-Littlewood
       prediction (aligned singular series S_aligned(pattern), expected
       count N_S(x) = S_aligned * (x/k) / (ln x)^w, first occurrence
       x* where N_S(x*) = 1) — a falsifiable HL test;
  (iii) state what would falsify saturation at k=32.

HL singular series for ALIGNED windows. The constellation is the family
of linear forms f_i(m) = k*m + o_i, o_i in S, |S| = w. Standard HL for
linear forms in one variable m: with omega(q) = #{ m mod q : some
f_i(m) ≡ 0 mod q },
    S_aligned = prod_q  (1 - omega(q)/q) / (1 - 1/q)^w .
  - q | k:  f_i(m) ≡ o_i (const). omega = q if some o_i ≡ 0 mod q
            (then count 0, S=0 — but admissibility (i) forbids this),
            else omega = 0  =>  factor 1/(1-1/q)^w  (for q=2: 2^w —
            aligned-even base makes every candidate odd).
  - q ∤ k:  k*m hits every residue, so f_i vanishes for exactly one m
            per offset; omega = #{ o_i mod q } = the admissibility-(ii)
            residue-coverage count. Admissible <=> omega < q for all q.

What would falsify saturation at k=32:
  - an INADMISSIBLE pattern realized at m > 0 (impossible by plain
    divisibility — a bug check; none ever found);
  - a LOW-weight admissible pattern still missing at x far beyond its
    HL prediction x* while higher-weight cousins appear (would break HL
    density ordering / k-tuple uniformity);
  - D_32(x) plateauing strictly below A_32 + 1 = 3574 as x -> infinity
    (an admissible pattern that provably never occurs — would refute
    the prime k-tuple conjecture's qualitative content).

Usage:
  python3 dk_first_occurrence.py --k 32 --xmax $((1<<30))
  python3 dk_first_occurrence.py --selftest
  python3 dk_first_occurrence.py --k 32 --xmax $((1<<28)) --validate-hl
"""

import argparse
import math
from collections import Counter

import numpy as np


# ----------------------------------------------------------------------
# primes / admissibility
# ----------------------------------------------------------------------
def primes_upto(limit):
    """All primes <= limit as a python list (limit small: sqrt(xmax))."""
    if limit < 2:
        return []
    s = np.ones(limit + 1, dtype=bool)
    s[:2] = False
    for p in range(2, int(limit ** 0.5) + 1):
        if s[p]:
            s[p * p::p] = False
    return [int(p) for p in np.nonzero(s)[0]]


def allowed_offsets(k):
    """Offsets that can ever be occupied (condition (i): not divisible
    by any prime q | k)."""
    div_primes = [q for q in primes_upto(k) if k % q == 0]
    return [o for o in range(k) if all(o % q != 0 for q in div_primes)]


def admissible_patterns(k):
    """Exact enumeration of admissible aligned-k patterns (frozensets).
    Brute over allowed offsets — fine through k≈40 (2^|allowed|)."""
    other_primes = [q for q in primes_upto(k) if k % q != 0]
    allowed = allowed_offsets(k)
    adm = []
    for mask in range(1 << len(allowed)):
        S = [allowed[i] for i in range(len(allowed)) if (mask >> i) & 1]
        ok = True
        for q in other_primes:
            if len({o % q for o in S}) == q:
                ok = False
                break
        if ok:
            adm.append(frozenset(S))
    return adm


def code_of(S, k):
    """Bitmask code of a pattern (matches realized-window codes)."""
    c = 0
    for o in S:
        c |= 1 << o
    return c


def set_of(code, k):
    return frozenset(o for o in range(k) if (code >> o) & 1)


# ----------------------------------------------------------------------
# segmented census
# ----------------------------------------------------------------------
def _segment_bounds(k, xmax, seg_size):
    """Break [0, xmax] into k-aligned chunks of size <= seg_size whose
    union of breakpoints includes every power of two in [k, xmax] (so
    every power-of-two checkpoint lands exactly on a chunk boundary)."""
    assert seg_size % k == 0
    brk = {0, xmax}
    p = 1
    while p < k:
        p <<= 1
    while p <= xmax:                       # powers of two >= k
        brk.add(p)
        p <<= 1
    m = seg_size
    while m < xmax:                        # multiples of seg_size
        brk.add(m)
        m += seg_size
    bounds = sorted(b - (b % k) for b in brk)   # k-align (no-op here)
    # dedup after alignment
    out = []
    for b in bounds:
        if not out or b != out[-1]:
            out.append(b)
    return out


def iter_window_codes(k, xmax, seg_size):
    """Yield (m_start, codes) for consecutive k-aligned segments, where
    codes[i] is the bitmask of prime offsets in window
    [k*(m_start+i), k*(m_start+i)+k). Segmented sieve, O(seg_size)
    memory."""
    base = primes_upto(int(xmax ** 0.5) + 1)
    shifts = np.arange(k, dtype=np.int64)
    bounds = _segment_bounds(k, xmax, seg_size)
    for lo, hi in zip(bounds[:-1], bounds[1:]):
        length = hi - lo
        if length < k:
            continue
        seg = np.ones(length, dtype=bool)
        for p in base:
            start = max(p * p, ((lo + p - 1) // p) * p)
            if start >= hi:
                continue
            seg[start - lo::p] = False
        if lo == 0:
            seg[0] = False
            if length > 1:
                seg[1] = False
        nwin = length // k
        win = seg[:nwin * k].reshape(nwin, k)
        # accumulate code low-memory: one int64 pass per offset
        codes = np.zeros(nwin, dtype=np.int64)
        for o in range(k):
            col = win[:, o]
            if col.any():
                codes |= col.astype(np.int64) << shifts[o]
        yield lo // k, codes


def census(k, xmax, seg_size=1 << 24, want_counts=False, progress=False):
    """Run the segmented census. Returns dict with:
      first[code] -> first-occurrence window index m
      checkpoints  -> list of (x, D, n_adm_found, n_missing)
      counts[code] -> #occurrences up to xmax  (only if want_counts)
    """
    adm = admissible_patterns(k)
    adm_codes = {code_of(S, k) for S in adm}
    A_k = len(adm)

    first = {}
    counts = Counter() if want_counts else None
    checkpoints = []
    pow2 = set()
    p = 1
    while p < k:
        p <<= 1
    while p <= xmax:
        pow2.add(p)
        p <<= 1

    for m_start, codes in iter_window_codes(k, xmax, seg_size):
        uniq, idx = np.unique(codes, return_index=True)
        for c, i in zip(uniq.tolist(), idx.tolist()):
            if c not in first:
                first[c] = m_start + i
        if want_counts:
            uq, cc = np.unique(codes, return_counts=True)
            for c, n in zip(uq.tolist(), cc.tolist()):
                counts[c] += n
        # checkpoint if this segment ends on a power of two (or xmax)
        seg_end = (m_start + len(codes)) * k
        if seg_end in pow2 or seg_end == xmax:
            seen = set(first)
            found = len(adm_codes & seen)
            checkpoints.append((seg_end, len(seen), found, A_k - found))
            if progress:
                print(f"  x=2^{int(round(math.log2(seg_end))):<2} "
                      f"D={len(seen):>5} adm={found:>5}/{A_k} "
                      f"missing={A_k - found:>4}", flush=True)

    return {"first": first, "checkpoints": checkpoints,
            "counts": counts, "adm_codes": adm_codes, "A_k": A_k}


# ----------------------------------------------------------------------
# Hardy-Littlewood aligned singular series + first-occurrence prediction
# ----------------------------------------------------------------------
def singular_series_aligned(S, k, pbound=20000):
    """S_aligned for the aligned constellation {k*m + o : o in S}."""
    w = len(S)
    if w == 0:
        return 1.0
    prod = 1.0
    for q in primes_upto(pbound):
        if k % q == 0:
            if any(o % q == 0 for o in S):
                return 0.0
            omega = 0
        else:
            omega = len({o % q for o in S})
        prod *= (1.0 - omega / q) / ((1.0 - 1.0 / q) ** w)
    return prod


def hl_count(S, k, x, pbound=20000, Sa=None):
    """Expected #aligned windows realizing the (>= S) constellation up
    to x:  N_S(x) = S_aligned * (x/k) / (ln x)^w."""
    w = len(S)
    if w == 0:
        return x / k
    if Sa is None:
        Sa = singular_series_aligned(S, k, pbound)
    if Sa == 0.0:
        return 0.0
    lx = math.log(x)
    return Sa * (x / k) / (lx ** w)


def predict_first_occ(S, k, pbound=20000, Sa=None):
    """Smallest x with N_S(x) >= 1 (HL first-occurrence estimate)."""
    w = len(S)
    if w == 0:
        return float(k)
    if Sa is None:
        Sa = singular_series_aligned(S, k, pbound)
    if Sa == 0.0:
        return math.inf
    lo, hi = float(max(k, 8)), 1.0e6
    while hl_count(S, k, hi, Sa=Sa) < 1.0 and hi < 1e305:
        hi *= 10.0
    if hl_count(S, k, hi, Sa=Sa) < 1.0:
        return math.inf
    for _ in range(300):
        mid = math.sqrt(lo * hi)
        if hl_count(S, k, mid, Sa=Sa) < 1.0:
            lo = mid
        else:
            hi = mid
    return hi


# ----------------------------------------------------------------------
# reporting
# ----------------------------------------------------------------------
def report(k, xmax, res, baseline_x=1 << 26, pbound=20000):
    first = res["first"]
    adm_codes = res["adm_codes"]
    A_k = res["A_k"]
    seen = set(first)
    found = adm_codes & seen
    missing = adm_codes - seen

    print(f"\n=== D_{k}(x) census, x up to 2^{int(round(math.log2(xmax)))}"
          f" ({xmax}) ===")
    print(f"A_{k} = {A_k}, target saturation D_{k} = A_{k}+1 = {A_k + 1}")
    print(f"\n{'x':>14} {'D_k(x)':>8} {'adm':>6} {'missing':>8} "
          f"{'D/(A+1)':>9}")
    for x, D, fnd, miss in res["checkpoints"]:
        print(f"2^{int(round(math.log2(x))):<2}{'':>9} {D:>8} {fnd:>6} "
              f"{miss:>8} {D / (A_k + 1):>9.4f}")

    # patterns first realized ABOVE the baseline x (the new finds)
    new = [(first[c], c) for c in found
           if first[c] * k >= baseline_x]
    new.sort()
    print(f"\nadmissible patterns first realized at x >= "
          f"2^{int(round(math.log2(baseline_x)))} (new this run): "
          f"{len(new)} "
          f"(by weight: "
          f"{dict(sorted(Counter(len(set_of(c, k)) for _, c in new).items()))})")
    if new:
        print(f"  {'w':>2} {'first-occ x':>13} {'HL x*':>13} "
              f"{'obs/x*':>8} {'N(occ)':>8}  offsets")
        for m, c in new[:40]:
            S = set_of(c, k)
            xocc = max(m * k, k)
            Sa = singular_series_aligned(S, k, pbound)
            xpred = predict_first_occ(S, k, pbound, Sa=Sa)
            ratio = xocc / xpred if math.isfinite(xpred) and xpred else \
                float('nan')
            nocc = hl_count(S, k, xocc, Sa=Sa)  # expected count by 1st occ
            print(f"  {len(S):>2} {xocc:>13} "
                  f"{('%.3e' % xpred):>13} {ratio:>8.2f} {nocc:>8.2f}  "
                  f"{sorted(S)}")
        if len(new) > 40:
            print(f"  ... ({len(new) - 40} more)")

    # still-missing patterns: rank by EXPECTED count N_S(xmax) (Poisson
    # mean). A pattern is 'overdue' only if its mean is large enough that
    # P(0 occurrences)=e^{-N} is tiny — that, not x* < xmax, is the sound
    # falsification trigger.
    print(f"\nstill MISSING at x=2^{int(round(math.log2(xmax)))}: "
          f"{len(missing)} "
          f"(by weight: "
          f"{dict(sorted(Counter(len(set_of(c, k)) for c in missing).items()))})")
    miss_info = []
    for c in missing:
        S = set_of(c, k)
        Sa = singular_series_aligned(S, k, pbound)
        N = hl_count(S, k, xmax, Sa=Sa)
        xstar = predict_first_occ(S, k, pbound, Sa=Sa)
        miss_info.append((N, len(S), xstar, c))
    miss_info.sort(reverse=True)  # most overdue (largest expected) first
    print("  most-expected still-missing (largest HL mean count "
          "N_S(xmax) first; P(0 occ)=e^{-N}):")
    print(f"  {'w':>2} {'N_S(xmax)':>10} {'P(0 occ)':>9} {'HL x*':>13}"
          f"  offsets")
    for N, w, xstar, c in miss_info[:12]:
        S = set_of(c, k)
        print(f"  {w:>2} {N:>10.3f} {math.exp(-N):>9.4f} "
              f"{('%.3e' % xstar):>13}  {sorted(S)}")

    # falsification check: a missing pattern whose mean count is large
    THRESH = 5.0  # P(0)=e^-5 = 0.7%
    overdue = [(N, w, c) for N, w, xstar, c in miss_info if N >= THRESH]
    if overdue:
        print(f"\n  [!] {len(overdue)} missing pattern(s) have HL mean "
              f"count N_S(xmax) >= {THRESH:.0f} (P(0 occ) < "
              f"{math.exp(-THRESH):.2%}) — STATISTICAL ANOMALY vs HL:")
        for N, w, c in sorted(overdue, reverse=True)[:10]:
            print(f"      w={w} N={N:.2f} offsets {sorted(set_of(c, k))}")
    else:
        mx = miss_info[0][0] if miss_info else 0.0
        print(f"\n  saturation OK: every still-missing pattern has HL "
              f"mean count N_S(xmax) < {THRESH:.0f} (largest {mx:.2f}, "
              f"P(0 occ)={math.exp(-mx):.2%}) — all missing are in the "
              f"expected HL first-occurrence tail, none statistically "
              f"overdue.")


# ----------------------------------------------------------------------
# HL constant validation (occurrence counts vs prediction)
# ----------------------------------------------------------------------
def validate_hl(k, xmax, res, pbound=20000):
    """Compare observed occurrence counts to HL prediction, per weight,
    for realized admissible patterns. Validates S_aligned and the
    overall normalization."""
    counts = res["counts"]
    adm_codes = res["adm_codes"]
    print(f"\n=== HL count validation, k={k}, x up to "
          f"2^{int(round(math.log2(xmax)))} ===")
    print(f"{'w':>2} {'#patterns':>10} {'obs_total':>12} "
          f"{'pred_total':>12} {'obs/pred':>9}")
    byw = {}
    for c in adm_codes:
        if c not in counts:
            continue
        S = set_of(c, k)
        w = len(S)
        obs = counts[c]
        pred = hl_count(S, k, xmax, pbound)
        d = byw.setdefault(w, [0, 0.0, 0])
        d[0] += obs
        d[1] += pred
        d[2] += 1
    for w in sorted(byw):
        obs, pred, n = byw[w]
        r = obs / pred if pred else float('nan')
        print(f"{w:>2} {n:>10} {obs:>12} {pred:>12.1f} {r:>9.3f}")
    print("(obs/pred ~ 1 validates the aligned singular series + "
          "normalization; ratio drifts low for high w because HL "
          "counts the >=S superset and dense patterns are below their "
          "asymptotic regime.)")


# ----------------------------------------------------------------------
# selftest
# ----------------------------------------------------------------------
def _direct_realized(k, x):
    """Reference: distinct exact patterns via a single full sieve."""
    n = x
    s = np.ones(n, dtype=bool)
    s[:2] = False
    for p in range(2, int(n ** 0.5) + 1):
        if s[p]:
            s[p * p::p] = False
    nwin = n // k
    win = s[:nwin * k].reshape(nwin, k)
    codes = win.astype(np.int64).dot(1 << np.arange(k, dtype=np.int64))
    return set(int(c) for c in np.unique(codes))


def selftest():
    print("=== selftest ===")

    # 1. admissibility counts (regression vs known A_k)
    assert len(admissible_patterns(8)) == 13
    assert len(admissible_patterns(16)) == 106
    assert len(admissible_patterns(32)) == 3573
    print("[1] A_8=13, A_16=106, A_32=3573 OK")

    # 2. code/set round-trip
    for S in [frozenset(), frozenset({1, 3, 7}), frozenset({0, 2, 30})]:
        assert set_of(code_of(S, 32), 32) == S
    print("[2] code<->set round-trip OK")

    # 3. segmented census == direct full sieve (codes), several seg sizes
    for k, n in [(8, 1 << 16), (16, 1 << 18), (32, 1 << 20)]:
        direct = _direct_realized(k, n)
        for seg in (k * 1, 1 << 14, 1 << 17):
            r = census(k, n, seg_size=max(seg, k))
            seen = set(int(c) for c in r["first"])
            assert seen == direct, (k, n, seg, len(seen), len(direct))
    print("[3] segmented census == direct full sieve (k=8,16,32; "
          "3 seg sizes) OK")

    # 4. saturation reproduced: D_8=14 by 2^16, D_16=107 by 2^18
    r8 = census(8, 1 << 16)
    assert r8["checkpoints"][-1][1] == 14, r8["checkpoints"][-1]
    assert r8["checkpoints"][-1][2] == 13  # all 13 admissible found
    r16 = census(16, 1 << 18)
    assert r16["checkpoints"][-1][1] == 107, r16["checkpoints"][-1]
    assert r16["checkpoints"][-1][2] == 106
    print("[4] saturation D_8=14 (2^16), D_16=107 (2^18) reproduced OK")

    # 5. the single inadmissible '+1' is exactly the m=0 window
    #    (k=16: primes {2,3,5,7,11,13} in [0,16))
    adm16 = {code_of(S, 16) for S in admissible_patterns(16)}
    seen16 = set(int(c) for c in r16["first"])
    inadm = seen16 - adm16
    assert len(inadm) == 1, inadm
    c0 = next(iter(inadm))
    assert r16["first"][c0] == 0  # window m=0
    assert set_of(c0, 16) == frozenset({2, 3, 5, 7, 11, 13})
    print("[5] the +1 inadmissible pattern is the m=0 small-prime "
          "window, only at m=0 OK")

    # 6. k=32 baseline at 2^26 matches the recorded census (3385, 189)
    r32 = census(32, 1 << 26)
    D, found, miss = (r32["checkpoints"][-1][1], r32["checkpoints"][-1][2],
                      r32["checkpoints"][-1][3])
    assert D == 3385, D
    assert miss == 189, miss
    mw = Counter(len(set_of(c, 32)) for c in (r32["adm_codes"] -
                                              set(r32["first"])))
    assert dict(sorted(mw.items())) == {6: 23, 7: 115, 8: 47, 9: 4}, mw
    print("[6] k=32 baseline at 2^26: D=3385, 189 missing "
          "(23/115/47/4 by weight 6/7/8/9) OK")

    # 7. checkpoints are monotone non-decreasing in D
    for r in (r8, r16, r32):
        ds = [c[1] for c in r["checkpoints"]]
        assert all(b >= a for a, b in zip(ds, ds[1:])), ds
    print("[7] D_k(x) monotone non-decreasing OK")

    # 8. singular series: admissible -> S_aligned > 0; an inadmissible
    #    (covers a prime) -> 0; factor 2^w shape on the q=2|k term
    Sgood = frozenset({1, 3, 7, 9, 13})        # admissible-ish odd set
    assert singular_series_aligned(Sgood, 32) > 0.0
    # inadmissible mod 3: residues {0,1,2} covered, e.g. odd offsets
    Sbad = frozenset({1, 3, 5})                # 1,3,5 mod3 = 1,0,2 cover
    assert singular_series_aligned(Sbad, 32) == 0.0
    # weight-1 single offset: S_aligned = prod_{q|k}1/(1-1/q)^1 *
    #   prod_{q∤k}(1-1/q)/(1-1/q) = 2  (only q=2|32 contributes 1/(1/2))
    s1 = singular_series_aligned(frozenset({1}), 32)
    assert abs(s1 - 2.0) < 1e-9, s1
    print("[8] singular series: admissible>0, inadmissible=0, "
          f"weight-1 == 2 (={s1:.6f}) OK")

    # 9. predicted first-occurrence is finite for every admissible
    #    pattern, and (per-weight medians) grows with weight
    adm32 = admissible_patterns(32)
    by_w = {}
    for S in adm32:
        if S:
            by_w.setdefault(len(S), []).append(S)
    med = {}
    for w, lst in by_w.items():
        ps = [predict_first_occ(S, 32) for S in lst]
        assert all(math.isfinite(p) for p in ps), w  # admissible => finite
        ps.sort()
        med[w] = ps[len(ps) // 2]
    ws = sorted(med)
    assert all(med[a] < med[b] for a, b in zip(ws, ws[1:])), med
    print(f"[9] HL first-occ finite for all admissible, per-weight "
          f"median strictly increasing (w1 {med[1]:.1e} -> "
          f"w9 {med[9]:.1e}) OK")

    # 10. HL count is monotincreasing in x and ordering-consistent:
    #     every pattern realized at 2^26 should have HL x* not absurdly
    #     larger than 2^26 (sanity: median obs<=pred slack)
    S = frozenset({1, 7, 11, 13, 17})
    assert hl_count(S, 32, 1 << 26) > hl_count(S, 32, 1 << 20)
    print("[10] HL count monotone in x OK")

    print("\nall selftests passed.")


# ----------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--k", type=int, default=32)
    ap.add_argument("--xmax", type=int, default=1 << 28)
    ap.add_argument("--seg", type=int, default=1 << 24,
                    help="segment size in integers (memory ~ seg bytes)")
    ap.add_argument("--baseline", type=int, default=1 << 26,
                    help="x above which a first occurrence counts as "
                         "'new this run'")
    ap.add_argument("--pbound", type=int, default=20000,
                    help="prime bound for the singular-series product")
    ap.add_argument("--validate-hl", action="store_true",
                    help="also count occurrences and compare to HL "
                         "prediction per weight")
    ap.add_argument("--selftest", action="store_true")
    args = ap.parse_args()

    if args.selftest:
        selftest()
        return

    seg = args.seg - (args.seg % args.k)
    res = census(args.k, args.xmax, seg_size=max(seg, args.k),
                 want_counts=args.validate_hl, progress=True)
    report(args.k, args.xmax, res, baseline_x=args.baseline,
           pbound=args.pbound)
    if args.validate_hl:
        validate_hl(args.k, args.xmax, res, pbound=args.pbound)


if __name__ == "__main__":
    main()
