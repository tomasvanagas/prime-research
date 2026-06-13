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

from lucy_dp_verification import (Q, lagrange_eval, eq_table, eq_point,
                                  mle_eval, s0_closed_form, ge_const_eval,
                                  w_div_eval, sumcheck, prover_tables,
                                  direct_pi)


# ----------------------------------------------------------------------
# verifier-side O(l) primitives
# ----------------------------------------------------------------------

def lt_const_eval(rho, M, l):
    """MLE of [s < M] on l-bit s at point rho. O(l). When M >= 2^l the
    predicate is constantly true (the p = 2, l = 1 boundary case)."""
    if M >= (1 << l):
        return 1
    return (Q + 1 - ge_const_eval(rho, M, l)) % Q


def aff_rel_eval(rho_x, rho_y, c, l):
    """MLE of the affine relation [y = 2x + c] over l-bit integers x, y,
    at field points rho_x, rho_y (bits MSB-first). Width-4 carry
    automaton (carry, prev-x-bit), processed LSB-first. O(l).

    c >= 0: y = 2x + c. At LSB-position i: t = x_{i-1} + c_i + carry;
            requires y_i = t & 1; accept iff final x_prev = carry = 0.
    c < 0:  y + d = 2x, d = -c > 0. t = y_i + d_i + carry; requires
            (2x)_i = x_{i-1} = t & 1; accept iff final carry == x_prev.
    """
    xb = list(reversed([int(v) % Q for v in rho_x]))
    yb = list(reversed([int(v) % Q for v in rho_y]))
    # process enough positions to exhaust the constant and the shift;
    # bits of x, y beyond position l-1 are virtual zeros
    L = max(l, abs(c).bit_length()) + 1
    st = {(0, 0): 1}  # state: (carry, xprev) -> weight
    for i in range(L):
        wx = [(Q + 1 - xb[i]) % Q, xb[i]] if i < l else [1, 0]
        wy = [(Q + 1 - yb[i]) % Q, yb[i]] if i < l else [1, 0]
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
                                + w * wy[ybit] % Q * wx[xi]) % Q
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
                                    + w * wy[yi] % Q * wx[xi]) % Q
        st = new
    # all bit positions through L-1 are consumed (x_{L-1} = 0 virtual),
    # so the equation holds iff no carry remains
    total = 0
    for (carry, xprev), w in st.items():
        if carry == 0:
            total = (total + w) % Q
    return total


def r_tilde_eval(rho_in, rho_out, p, l, wv, wu):
    """Agreed extension R~_j(rho_out, rho_in) =
    sum_{a,b} wv[a] wu[b] AFF~_{a - b p}(x=rho_in, y=rho_out). O(l)."""
    tot = 0
    for a in (0, 1):
        for b in (0, 1):
            tot = (tot + wv[a] * wu[b] % Q
                   * aff_rel_eval(rho_in, rho_out, a - b * p, l)) % Q
    return tot


# ----------------------------------------------------------------------
# inner protocol: delegated evaluation of W~(r_v, r_u)
# ----------------------------------------------------------------------

def bit_weights(r):
    r = int(r) % Q
    return [(Q + 1 - r) % Q, r]


def inner_chain_vectors(r_v, r_u, p, l, n):
    """Honest prover: chain v_0..v_n over the 2^l-embedded state space.
    Support is kept inside the valid states [0, p) — transitions into
    junk states [p, 2^l) must be excluded so the layer identity
    v_{j+1}(s') = LT(s') * sum_s R(s', s) LT(s) v_j(s) holds on the
    WHOLE Boolean cube (the agreed extension kills junk via LT)."""
    dim = 1 << l
    vs = [np.zeros(dim, dtype=np.uint64)]
    vs[0][0] = 1
    for j in range(n):
        wv, wu = bit_weights(r_v[j]), bit_weights(r_u[j])
        v = vs[-1]
        nv = np.zeros(dim, dtype=np.uint64)
        for s in range(p):
            if v[s] == 0:
                continue
            for a in (0, 1):
                for b in (0, 1):
                    y = 2 * s + a - b * p
                    if 0 <= y < p:
                        nv[y] = (nv[y] + int(v[s]) * wv[a] % Q
                                 * wu[b]) % Q
        vs.append(nv)
    return vs


def inner_r_table(p, l, wv, wu):
    """Prover table of R~_j on the Boolean (sig', sig) cube,
    index = (sig' << l) | sig."""
    dim = 1 << l
    Rt = np.zeros(dim * dim, dtype=np.uint64)
    for s in range(dim):
        for a in (0, 1):
            for b in (0, 1):
                y = 2 * s + a - b * p
                if 0 <= y < dim:
                    Rt[(y << l) | s] = (Rt[(y << l) | s]
                                        + wv[a] * wu[b]) % Q
    return Rt


def inner_verify_W(claimW, r_v, r_u, p, n, rng, stats, lie=False):
    """Verify claimW = W~(r_v, r_u) by GKR over the automaton chain.
    Verifier cost O(n*l). Returns ok."""
    l = max(1, (p - 1).bit_length())
    dim = 1 << l
    t0 = time.perf_counter()
    vs = inner_chain_vectors(r_v, r_u, p, l, n)
    if lie:  # prover corrupts the chain consistently from the middle
        vs[n // 2] = vs[n // 2].copy()
        vs[n // 2][0] = (int(vs[n // 2][0]) + 1) % Q
        for j in range(n // 2, n):
            wv, wu = bit_weights(r_v[j]), bit_weights(r_u[j])
            nv = np.zeros(dim, dtype=np.uint64)
            for s in range(p):
                for a in (0, 1):
                    for b in (0, 1):
                        y = 2 * s + a - b * p
                        if 0 <= y < p:
                            nv[y] = (nv[y] + int(vs[j][s]) * wv[a] % Q
                                     * wu[b]) % Q
            vs[j + 1] = nv
    lt_tab = (np.arange(dim) < p).astype(np.uint64)
    stats["t_prover"] += time.perf_counter() - t0

    # sum-check #0: claimW = sum_sig LT(sig) * v_n(sig)
    ok, r0, final0, scal0 = sumcheck(claimW,
                                     {"LT": lt_tab.copy(),
                                      "V": vs[n].copy()},
                                     [(1, ["LT", "V"])], 2, rng)
    stats["comm"] += l * 3
    if not ok:
        return False
    t0 = time.perf_counter()
    if final0 != lt_const_eval(r0, p, l) * scal0["V"] % Q:
        stats["t_verifier"] += time.perf_counter() - t0
        return False
    stats["t_verifier"] += time.perf_counter() - t0

    rho, claim = r0, scal0["V"]  # claim about v~_n at rho
    for j in range(n - 1, -1, -1):
        wv, wu = bit_weights(r_v[j]), bit_weights(r_u[j])
        # claim = sum_{sig',sig} eq(rho,sig') LT(sig') R_j(sig',sig)
        #                        LT(sig) v~_j(sig)
        t0 = time.perf_counter()
        eqt = eq_table(rho)
        EQb = np.repeat(eqt, dim)
        LTo = np.repeat(lt_tab, dim)
        Rt = inner_r_table(p, l, wv, wu)
        LTi = np.tile(lt_tab, dim)
        Vb = np.tile(vs[j], dim)
        stats["t_prover"] += time.perf_counter() - t0
        ok, rr, final, scal = sumcheck(claim,
                                       {"EQ": EQb, "LTo": LTo, "R": Rt,
                                        "LTi": LTi, "V": Vb},
                                       [(1, ["EQ", "LTo", "R", "LTi",
                                             "V"])], 5, rng)
        stats["comm"] += 2 * l * 6
        if not ok:
            return False
        rho_out, rho_in = rr[:l], rr[l:]
        t0 = time.perf_counter()
        expect = (eq_point(rho, rho_out)
                  * lt_const_eval(rho_out, p, l) % Q
                  * r_tilde_eval(rho_in, rho_out, p, l, wv, wu) % Q
                  * lt_const_eval(rho_in, p, l) % Q
                  * scal["V"]) % Q
        okf = final == expect
        stats["t_verifier"] += time.perf_counter() - t0
        if not okf:
            return False
        rho, claim = rho_in, scal["V"]
    # base: v~_0 = MLE of e_0
    t0 = time.perf_counter()
    base = 1
    for ri in rho:
        base = base * ((Q + 1 - int(ri)) % Q) % Q
    okb = claim == base
    stats["t_verifier"] += time.perf_counter() - t0
    return okb


# ----------------------------------------------------------------------
# outer protocol with delegated wiring
# ----------------------------------------------------------------------

def verify_layer_delegated(i, p, n, z, claim, S_prev_field, rng, stats,
                           lie_inner=False):
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
    if not inner_verify_W(wv_claim, r_v, r_u, p, n, rng, stats):
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


def run_protocol_delegated(n, rng, corrupt_layer=None, lie_inner_at=None):
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
            lie_inner=(lie_inner_at == i))
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
    # 4. tiny end-to-end with delegated wiring
    res = run_protocol_delegated(10, np.random.default_rng(4))
    assert res["accepted"] and res["claimed"] == direct_pi((1 << 10) - 1)
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
    args = ap.parse_args()
    if args.selftest:
        selftest()
    elif args.wiring_bench:
        wiring_bench()
    else:
        main(args.n, args.cheat_trials, args.seed)
