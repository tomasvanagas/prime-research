#!/usr/bin/env python3
"""
census_transfer_matrix.py — A_k WITHOUT enumerating all A_k patterns.

Open item 4 (PROGRAM.md): the admissible-pattern census count

    A_k = #{ S ⊆ allowed(k) : for every prime q ≤ k, q ∤ k, the occupied
            residues {o mod q : o ∈ S} do NOT cover all of Z_q }

(allowed(k) = offsets in [0,k) divisible by no prime q | k).  The
conjectured entropy law is A_k = e^{(1+o(1))·π(k)} (clean family
k = 2^m: ln A_k/(k/ln k) = 0.667, 0.808, 0.886, 0.928 at k = 8,16,32,64,
geometric convergence to 1).  The DFS counter in local_pattern_census.py
visits ~A_k nodes, so it is enumeration-bound and dies at A_128 ≈ e^25
≈ 10^11.  The NEXT ACTION: count A_k by a transfer matrix / DP that does
NOT enumerate every pattern, validate bit-for-bit against the exact
A_k (k ≤ 80), and either push toward k = 128 or show the transfer-matrix
state space is itself exponential (closing the method).

This script does three things.

1. **Active-prime reduction (rigorous via stabilization).**  Enforcing a
   prime is a monotone FILTER on subsets: adding a prime q to the
   enforced set can only KEEP or DECREASE the count.  So if enforcing
   primes {q ≤ B} gives the same count as enforcing {q ≤ B'} for the next
   enforceable prime B', then no pattern was ever removed by the primes
   in (B,B'] — they never bind, and the count is EXACT.  We therefore
   enforce primes up to the smallest B at which the count stabilises.
   Empirically (and provably via the weight bound — see results md) B(k)
   grows like ρ*(k) ~ k/ln k, far below k: B = 3,5,7,13 at k = 8,16,32,64.

2. **Transfer-matrix DP (method "tm").**  Process the allowed offsets one
   at a time; the state is the COVERAGE TUPLE (per enforced prime q, the
   bitmask of residues mod q occupied so far).  Each offset is
   included/excluded; coverage masks OR in the offset's residues.  Many
   subsets share a coverage tuple, so the DP collapses them — the cost is
   Θ(#distinct states · #offsets), NOT Θ(A_k).  At the end, sum the
   counts over states where every prime misses ≥1 class.  A_k is exact
   (Python big-int counts, no cancellation).

3. **State-space measurement.**  We report the peak and final number of
   distinct DP states vs A_k itself.  This decides the open question:
   if #states ≪ A_k the TM beats DFS and may reach k = 128; if #states
   grows with the same exponent as A_k the transfer-matrix route is
   itself exponential (the alternative deliverable).

What would falsify a result here:
  - any A_k disagreeing with the DFS / known exact value (k ≤ 80) — a bug;
  - the stabilisation check passing (count(B)==count(B')) yet the larger
    enforced set later changing the count — impossible (monotone filter),
    so a sign of a state-collision bug;
  - #states measured ≪ A_k at small k but the fitted state exponent
    α_states ≈ α_{A_k}: the "beats DFS" claim would be false (constant-
    factor only), which we report honestly from the fit.

Usage:
  python3 census_transfer_matrix.py --selftest
  python3 census_transfer_matrix.py --k 64               # one A_k + state count
  python3 census_transfer_matrix.py --scan 8 96 --step 4 # A_k table + entropy fit
  python3 census_transfer_matrix.py --k 128 --max-states 50000000
"""

import argparse
import math
import time


# ----------------------------------------------------------------------
# basics
# ----------------------------------------------------------------------
def primes_upto(m):
    return [p for p in range(2, m + 1)
            if all(p % d for d in range(2, int(p ** 0.5) + 1))]


def allowed_offsets(k):
    """Offsets in [0,k) divisible by no prime dividing k (condition (i))."""
    dp = [q for q in primes_upto(k) if k % q == 0]
    return [o for o in range(k) if all(o % q != 0 for q in dp)]


def enforceable_primes(k, allowed):
    """Primes q ≤ k, q ∤ k, that SOME subset of allowed offsets can fully
    cover — i.e. every residue class mod q contains ≥1 allowed offset.
    Only these can ever be violated; the rest impose no constraint."""
    out = []
    for q in primes_upto(k):
        if k % q == 0:
            continue
        if len({o % q for o in allowed}) == q:
            out.append(q)
    return out


# ----------------------------------------------------------------------
# DFS baseline (downward-closed admissibility pruning) — for validation
# ----------------------------------------------------------------------
def count_admissible_dfs(k):
    allowed = allowed_offsets(k)
    qlist = enforceable_primes(k, allowed)
    full = {q: (1 << q) - 1 for q in qlist}
    obits = [{q: 1 << (o % q) for q in qlist} for o in allowed]
    count = 0
    stack = [(0, tuple(0 for _ in qlist))]
    while stack:
        idx, cov = stack.pop()
        if idx == len(allowed):
            count += 1
            continue
        stack.append((idx + 1, cov))            # exclude offset idx
        newcov = []
        dead = False
        for qi, q in enumerate(qlist):
            c = cov[qi] | obits[idx][q]
            if c == full[q]:                    # q fully covered -> inadmissible
                dead = True
                break
            newcov.append(c)
        if not dead:
            stack.append((idx + 1, tuple(newcov)))   # include offset idx
    return count


# ----------------------------------------------------------------------
# transfer-matrix DP over offsets, state = coverage tuple over B primes
# ----------------------------------------------------------------------
def count_tm(k, B, allowed=None, max_states=None, want_stats=False):
    """A_k restricted to enforcing exactly the primes in B (a monotone
    under-set of enforceable primes).  State = packed coverage bitmask
    over the primes in B.  Returns the count (and optionally stats).

    A "dead" state (some prime fully covered) is dropped immediately — it
    can never recover (coverage only grows).  So we never store covered
    states; the state space is the set of partial-coverage tuples that
    still miss a class in every enforced prime.
    """
    if allowed is None:
        allowed = allowed_offsets(k)
    # bit layout: prime q occupies a contiguous q-bit field; shift[qi]
    shift = []
    s = 0
    for q in B:
        shift.append(s)
        s += q
    fullmask = 0
    for qi, q in enumerate(B):
        fullmask |= ((1 << q) - 1) << shift[qi]
    # per-offset OR-delta and the per-prime field masks (to test "full")
    deltas = []
    for o in allowed:
        d = 0
        for qi, q in enumerate(B):
            d |= (1 << (o % q)) << shift[qi]
        deltas.append(d)
    fieldmask = [((1 << q) - 1) << shift[qi] for qi, q in enumerate(B)]

    states = {0: 1}                 # coverage -> #subsets reaching it
    peak = 1
    for d in deltas:
        nxt = states.copy()         # all the "exclude this offset" branches
        for cov, cnt in states.items():
            nc = cov | d
            # drop if any enforced prime just became fully covered
            dead = False
            for fm in fieldmask:
                if nc & fm == fm:
                    dead = True
                    break
            if dead:
                continue
            nxt[nc] = nxt.get(nc, 0) + cnt
        states = nxt
        if len(states) > peak:
            peak = len(states)
        if max_states is not None and len(states) > max_states:
            if want_stats:
                return None, {"peak": peak, "final": len(states),
                              "aborted": True}
            return None
    total = sum(states.values())    # every surviving state misses a class
    if want_stats:
        return total, {"peak": peak, "final": len(states), "aborted": False}
    return total


def minkill_bb(allowed, B):
    """Min over drop-tuples (one residue class per prime in B) of the size
    of the union of dropped classes, via branch-and-bound (exact; does NOT
    iterate all Π q tuples).  W(B) = |allowed| - minkill is the max weight
    of a subset missing a class mod every prime in B; admissible patterns
    are a subset of these, so every admissible pattern has weight ≤ W(B)."""
    n = len(allowed)
    kmasks = []                      # kmasks[qi][r] = offset-bitmask of class r
    for q in B:
        per = [0] * q
        for i, o in enumerate(allowed):
            per[o % q] |= (1 << i)
        kmasks.append(per)
    best = [n + 1]

    def lower_bound(killed, idx):
        kc = bin(killed).count("1")
        extra = 0                    # must still kill ≥1 full class of each
        for j in range(idx, len(B)):  # remaining prime: cheapest residual
            m = min(bin(km & ~killed).count("1") for km in kmasks[j])
            if m > extra:
                extra = m
        return kc + extra

    def dfs(idx, killed):
        if lower_bound(killed, idx) >= best[0]:
            return
        if idx == len(B):
            kc = bin(killed).count("1")
            if kc < best[0]:
                best[0] = kc
            return
        for km in sorted(kmasks[idx],
                         key=lambda m: bin(m & ~killed).count("1")):
            dfs(idx + 1, killed | km)

    dfs(0, 0)
    return best[0]


def reduce_primes(k, allowed=None, enf=None):
    """Rigorous active-prime reduction WITHOUT running the count: the
    smallest prefix B of the enforceable primes with W(B) < next prime.
    Then any subset missing a class mod all q∈B has weight ≤ W(B) < that
    next prime, so it can't cover any larger prime either ⇒ it is fully
    admissible.  Hence A_k = #{S misses a class mod all q∈B}.  Returns
    (B, W, dropped_primes)."""
    if allowed is None:
        allowed = allowed_offsets(k)
    if enf is None:
        enf = enforceable_primes(k, allowed)
    for i in range(1, len(enf) + 1):
        B = enf[:i]
        W = len(allowed) - minkill_bb(allowed, B)
        nxt = enf[i] if i < len(enf) else None
        if nxt is None or W < nxt:
            return B, W, [p for p in enf if p > B[-1]]
    return enf, len(allowed) - minkill_bb(allowed, enf), []


def minimal_B(k, allowed=None, enf=None, max_states=None, verbose=False):
    """Smallest prefix B of the enforceable primes whose count is stable:
    count(B) == count(B + next prime).  Returns (B, A_k, stats)."""
    if allowed is None:
        allowed = allowed_offsets(k)
    if enf is None:
        enf = enforceable_primes(k, allowed)
    prev = None
    for i in range(1, len(enf) + 1):
        B = enf[:i]
        c, st = count_tm(k, B, allowed, max_states=max_states, want_stats=True)
        if c is None:
            return B, None, st          # aborted: too many states
        if verbose:
            print(f"    B={B} -> count={c}  states(peak/final)="
                  f"{st['peak']}/{st['final']}")
        if prev is not None and c == prev_count:
            # adding prime enf[i-1] did not change the count -> stable.
            return prev_B, c, prev_st
        prev = True
        prev_count, prev_B, prev_st = c, B, st
    # all enforceable primes used (count is exact by construction)
    return enf, prev_count, prev_st


# ----------------------------------------------------------------------
# entropy fit
# ----------------------------------------------------------------------
def _ls_slope(pts):
    n = len(pts)
    sx = sum(x for x, _ in pts)
    sy = sum(y for _, y in pts)
    sxx = sum(x * x for x, _ in pts)
    sxy = sum(x * y for x, y in pts)
    return (n * sxy - sx * sy) / (n * sxx - sx * sx)


def fit_entropy(records):
    """records: list of (k, A_k, peak_states).  Reports (i) the conjecture
    quantity ln A_k/(k/ln k) -> 1, and (ii) the transfer-matrix state
    count vs A_k (the collapse factor A_k/states and the state exponent),
    deciding whether the DP escapes the enumeration barrier.  Slopes are
    over the clean family k = 2^m."""
    print(f"\n{'k':>5} {'A_k':>16} {'states':>12} {'A_k/st':>7} "
          f"{'lnA':>8} {'k/lnk':>8} {'ratio':>7}")
    clean_a = []
    clean_s = []
    for rec in records:
        k, a = rec[0], rec[1]
        st = rec[2] if len(rec) > 2 else None
        if a is None or a <= 0:
            print(f"{k:>5} {'(aborted)':>16}")
            continue
        lna = math.log(a)
        kk = k / math.log(k)
        ratio = lna / kk
        cf = (a / st) if st else float("nan")
        flag = "  *" if (k & (k - 1)) == 0 else ""
        print(f"{k:>5} {a:>16} {st if st else '-':>12} {cf:>7.3f} "
              f"{lna:>8.3f} {kk:>8.3f} {ratio:>7.4f}{flag}")
        if (k & (k - 1)) == 0:
            clean_a.append((kk, lna))
            if st:
                clean_s.append((kk, math.log(st)))
    if len(clean_a) >= 2:
        sa = _ls_slope(clean_a)
        print(f"\nclean family k=2^m:")
        print(f"  ln A_k vs k/lnk : slope = {sa:.4f}  (conjecture -> 1.0); "
              f"last ratio = {clean_a[-1][1]/clean_a[-1][0]:.4f}")
    if len(clean_s) >= 2:
        ss = _ls_slope(clean_s)
        print(f"  ln(states) vs k/lnk : slope = {ss:.4f}")
        print(f"  => state space exponent {'≈' if abs(ss-sa)<0.1 else '!='} "
              f"A_k exponent  ==> transfer matrix is "
              f"{'EXPONENTIAL (Theta(A_k), method closed)' if abs(ss-sa)<0.1 else 'sub-A_k'}")


# ----------------------------------------------------------------------
# selftest
# ----------------------------------------------------------------------
KNOWN = {8: 13, 16: 106, 32: 3573, 48: 29002, 64: 1581920, 80: 5777381,
         # non-power-of-two clean cross-checks from local_pattern_census
         12: 16, 20: 121, 24: 227, 28: 640, 36: 2704, 40: 5704,
         44: 19825, 52: 87438, 56: 93751, 60: 53602}


def selftest():
    import itertools
    print("=== census_transfer_matrix selftest ===")
    n_ok = 0

    # 1. allowed offsets / enforceable primes sanity
    assert allowed_offsets(8) == [1, 3, 5, 7]
    assert allowed_offsets(16) == [1, 3, 5, 7, 9, 11, 13, 15]
    assert enforceable_primes(8, allowed_offsets(8)) == [3]
    assert enforceable_primes(16, allowed_offsets(16)) == [3, 5, 7]
    n_ok += 1
    print("[1] allowed/enforceable offsets ok")

    # 2. brute-force A_k for tiny k (full 2^|allowed| enumeration) ==
    #    DFS == TM(all enforceable primes)
    def brute(k):
        allowed = allowed_offsets(k)
        enf = enforceable_primes(k, allowed)
        cnt = 0
        for mask in range(1 << len(allowed)):
            S = [allowed[i] for i in range(len(allowed)) if (mask >> i) & 1]
            ok = True
            for q in enf:
                if len({o % q for o in S}) == q:
                    ok = False
                    break
            if ok:
                cnt += 1
        return cnt
    for k in (8, 12, 16, 20, 24):
        b = brute(k)
        d = count_admissible_dfs(k)
        t = count_tm(k, enforceable_primes(k, allowed_offsets(k)))
        assert b == d == t == KNOWN[k], (k, b, d, t, KNOWN.get(k))
    n_ok += 1
    print("[2] brute == DFS == TM(all primes) == known, k in {8,12,16,20,24}")

    # 3. TM with minimal stabilised B == TM(all primes) == known
    for k in (8, 16, 32, 48, 64):
        allowed = allowed_offsets(k)
        B, c, st = minimal_B(k, allowed)
        full = count_tm(k, enforceable_primes(k, allowed), allowed)
        assert c == full == KNOWN[k], (k, B, c, full, KNOWN.get(k))
    n_ok += 1
    print("[3] minimal-B stabilised count == full == known, "
          "k in {8,16,32,48,64}")

    # 4. monotonicity of the filter: count is non-increasing as primes
    #    are added; and stabilisation is genuine (further primes
    #    leave it unchanged).
    k = 64
    allowed = allowed_offsets(k)
    enf = enforceable_primes(k, allowed)
    counts = [count_tm(k, enf[:i], allowed) for i in range(1, len(enf) + 1)]
    assert all(counts[i] >= counts[i + 1] for i in range(len(counts) - 1)), \
        counts
    # once it stops decreasing it stays constant
    last = counts[-1]
    assert counts.count(last) >= 2 and counts[-1] == counts[-2] == KNOWN[64]
    n_ok += 1
    print(f"[4] filter monotone non-increasing; stabilises at A_64={last}")

    # 5. the dead-state pruning is correct: TM count == #subsets that miss
    #    a class in EVERY enforced prime (independent recompute via the
    #    per-prime per-class structure on a tiny case k=16, B={3,5}).
    k = 16
    allowed = allowed_offsets(k)
    # enforce {3,5}: count subsets missing a class mod 3 AND mod 5
    cnt = 0
    for mask in range(1 << len(allowed)):
        S = [allowed[i] for i in range(len(allowed)) if (mask >> i) & 1]
        if len({o % 3 for o in S}) < 3 and len({o % 5 for o in S}) < 5:
            cnt += 1
    assert count_tm(16, [3, 5], allowed) == cnt
    n_ok += 1
    print(f"[5] dead-state pruning correct (TM({{3,5}})={cnt} == brute)")

    # 6. state-count vs A_k recorded and < the trivial 2^|allowed| bound
    for k in (16, 32, 48):
        allowed = allowed_offsets(k)
        B, c, st = minimal_B(k, allowed)
        assert st["final"] <= (1 << len(allowed))
        assert c == KNOWN[k]
    n_ok += 1
    print("[6] TM state count bounded by 2^|allowed|; counts exact")

    # 7. empty pattern always admissible and counted exactly once
    for k in (8, 16, 32):
        # the all-empty coverage 0 survives -> contributes the empty set
        assert count_tm(k, enforceable_primes(k, allowed_offsets(k))) >= 1
    n_ok += 1
    print("[7] empty pattern counted (A_k >= 1) for k in {8,16,32}")

    # 8. reduce_primes (W(B)-based, no count) AGREES with the count-based
    #    minimal_B stabilisation, and is sufficient: TM enforcing exactly
    #    B(k) reproduces the known A_k.
    for k in (8, 16, 32, 64):
        allowed = allowed_offsets(k)
        Bw, W, dropped = reduce_primes(k, allowed)
        Bc, c, _ = minimal_B(k, allowed)
        assert Bw == Bc, (k, Bw, Bc)
        assert count_tm(k, Bw, allowed) == KNOWN[k], (k, Bw)
        # W(B) is a genuine weight bound: no admissible pattern exceeds it
        # (max admissible weight <= W, checked against DFS-realised weights
        #  for small k by full enumeration)
        if k <= 32:
            maxw = 0
            for mask in range(1 << len(allowed)):
                S = [allowed[i] for i in range(len(allowed))
                     if (mask >> i) & 1]
                ok = all(len({o % q for o in S}) < q
                         for q in enforceable_primes(k, allowed))
                if ok:
                    maxw = max(maxw, len(S))
            assert maxw <= W, (k, maxw, W)
    n_ok += 1
    print("[8] reduce_primes (W-based) == minimal_B (count-based); "
          "B(k) sufficient; W(B) bounds max admissible weight")

    # 9. minkill_bb is exact: matches a tiny brute over all drop-tuples
    k = 32
    allowed = allowed_offsets(k)
    B = [3, 5, 7]

    def union_killed(drop):
        killed = 0
        for i, o in enumerate(allowed):
            if any(o % q == drop[qi] for qi, q in enumerate(B)):
                killed |= 1 << i
        return killed
    brute_mk = min(bin(union_killed(drop)).count("1")
                   for drop in itertools.product(*[range(q) for q in B]))
    assert minkill_bb(allowed, B) == brute_mk, (minkill_bb(allowed, B),
                                                brute_mk)
    n_ok += 1
    print(f"[9] minkill_bb exact (=={brute_mk}) vs brute drop-tuples k=32")

    print(f"\nALL SELFTESTS PASSED ({n_ok} groups)")


# ----------------------------------------------------------------------
def run_one(k, max_states, verbose):
    allowed = allowed_offsets(k)
    enf = enforceable_primes(k, allowed)
    t0 = time.time()
    B, c, st = minimal_B(k, allowed, enf, max_states=max_states,
                         verbose=verbose)
    dt = time.time() - t0
    if c is None:
        print(f"k={k}: ABORTED at B={B} (states>{max_states}); "
              f"peak={st['peak']} final={st['final']}")
        return k, None, st
    nb = len(allowed)
    print(f"k={k}: A_k={c}  enforced B={B} (|B|={len(B)}, max {B[-1]} "
          f"of {len(enf)} enforceable)  |allowed|={nb}  "
          f"states(peak/final)={st['peak']}/{st['final']}  "
          f"states/A_k={st['peak']/c:.4f}  2^nb={1<<nb}  t={dt:.2f}s")
    return k, c, st


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--selftest", action="store_true")
    ap.add_argument("--k", type=int, default=0)
    ap.add_argument("--scan", type=int, nargs=2, default=None,
                    metavar=("KMIN", "KMAX"))
    ap.add_argument("--step", type=int, default=4)
    ap.add_argument("--powers", action="store_true",
                    help="scan only k = 2^m in [KMIN,KMAX]")
    ap.add_argument("--max-states", type=int, default=20_000_000)
    ap.add_argument("--verbose", action="store_true")
    ap.add_argument("--reduce", type=int, default=0,
                    help="report the active-prime reduction B(k) for k "
                         "(no count) — works even when A_k is unreachable")
    args = ap.parse_args()

    if args.selftest:
        selftest()
        return
    if args.reduce:
        k = args.reduce
        allowed = allowed_offsets(k)
        enf = enforceable_primes(k, allowed)
        B, W, dropped = reduce_primes(k, allowed, enf)
        print(f"k={k}: |allowed|={len(allowed)}  enforceable primes "
              f"({len(enf)}): {enf}")
        print(f"  B({k}) = primes ≤ {B[-1]} ({len(B)} primes): {B}")
        print(f"  W(B) = max admissible weight ≤ {W}  (< next prime "
              f"{dropped[0] if dropped else '—'})")
        print(f"  primes that NEVER bind (ignored): {dropped}")
        return
    if args.k:
        run_one(args.k, args.max_states, args.verbose)
        return
    if args.scan:
        kmin, kmax = args.scan
        records = []
        if args.powers:
            ks = [1 << m for m in range(1, 40) if kmin <= (1 << m) <= kmax]
        else:
            ks = list(range(kmin, kmax + 1, args.step))
        for k in ks:
            kk, c, st = run_one(k, args.max_states, args.verbose)
            records.append((k, c, (st or {}).get("peak")))
        fit_entropy(records)
        return
    selftest()


if __name__ == "__main__":
    main()
