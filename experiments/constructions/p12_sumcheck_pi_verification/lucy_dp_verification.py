#!/usr/bin/env python3
"""
lucy_dp_verification.py — interactive proof for the EXACT value of pi(x).

An untrusted prover runs the Lucy_Hedgehog / Legendre sieve DP and then
convinces a verifier of the exact value of pi(x) with UNCONDITIONAL
(information-theoretic) soundness. The verifier checks every layer of
the DP recursion via sum-check without ever sieving, touching only
O(n)-size objects per layer plus an O(n*p) automaton evaluation for the
floor-division wiring.

Protocol (new object, S491). State function after sieving by the first
i primes, defined for ALL v in [0, 2^n):

    S_0(v) = max(v - 1, 0)                       (count of 2..v)
    S_i(v) = S_{i-1}(v) - b_i(v) * (S_{i-1}(floor(v/p_i)) - (i-1))

where b_i(v) = [v >= p_i^2]. The recurrence holds pointwise for every v
(Lucy_Hedgehog), and S_K(x) = pi(x) once K = pi(sqrt(x)). Crucially
S_{i-1}(p_i - 1) = i - 1 is a PUBLIC constant (the primes below p_i),
so each layer's correction scalar needs no extra claim.

Per layer i (claim: tilde{S_i}(z) = C):
  Phase A  (n rounds, degree <= 3): sum-check of
      C = sum_v eq(z,v) * [ S(v) - B(v) * (G(v) - (i-1)) ]
    with G(v) := S_{i-1}(floor(v/p_i)) as an independent table; ends
    with prover claims s1 = tilde{S}(r_v), g1 = tilde{G}(r_v); verifier
    evaluates eq(z, r_v) in O(n) and tilde{B}(r_v) by a 3-state
    comparator automaton in O(n), and checks the combination.
  Phase B  (n rounds, degree <= 2): sum-check of
      g1 = sum_u W(r_v, u) * S(u),   W(v,u) = [u = floor(v/p_i)]
    ends with claim s2 = tilde{S}(r_u); verifier evaluates
    tilde{W}(r_v, r_u) by the LONG-DIVISION AUTOMATON: p states
    (the running remainder), reading bits of v MSB-first and emitting
    the quotient bits of u — O(n*p) per layer. This is the structured-
    wiring primitive that makes the verifier sublinear.
  Line batching: the two claims s1@r_v, s2@r_u about tilde{S_{i-1}}
    reduce to one claim at a random point on the line through them
    (prover sends the degree-<=n restriction).
Base layer: tilde{S_0}(z) = Int(z) - 1 + prod_j (1 - z_j), a closed
form the verifier evaluates in O(n). Int(z) = sum_j 2^{n-1-j} z_j.

Costs: prover O~(K * 2^n) here (table-based demo; a production prover
works on the O(sqrt x)-size compressed Lucy table). Verifier
O(n * sum_{p <= sqrt x} p) field ops with the automaton wiring;
communication O(K * n) field elements. Soundness error <= ~7nK/q,
unconditional. Field q = 2^31 - 1 for numpy-exact arithmetic (a
production protocol lifts to q ~ 2^61 or an extension field).

What would falsify this construction: an accepted honest run whose
output differs from a direct sieve's pi(x); a cheating prover (wrong
claim, adaptive round patching, or a corrupted DP table propagated
consistently downstream) accepted at a rate above the soundness bound;
or verifier work scaling with 2^n instead of n * sum p.

Usage:
  python3 lucy_dp_verification.py --selftest
  python3 lucy_dp_verification.py --n 16 --cheat-trials 25
"""

import argparse
import time

import numpy as np

Q = (1 << 31) - 1


# ----------------------------------------------------------------------
# small field helpers
# ----------------------------------------------------------------------

def lagrange_eval(ys, r):
    """Evaluate the polynomial through (j, ys[j]), j = 0..deg, at r."""
    k = len(ys) - 1
    r = int(r) % Q
    total = 0
    for i in range(k + 1):
        num, den = 1, 1
        for j in range(k + 1):
            if j == i:
                continue
            num = (num * ((r - j) % Q)) % Q
            den = (den * ((i - j) % Q)) % Q
        total = (total + ys[i] * num % Q * pow(den, Q - 2, Q)) % Q
    return total


def eq_table(z):
    """Table of eq~(z, w) over Boolean w; index bit j (MSB-first) pairs
    with z[j]. PROVER-side object: O(2^n)."""
    t = np.array([1], dtype=np.uint64)
    for zj in reversed([int(v) % Q for v in z]):
        a = (t * ((Q + 1 - zj) % Q)) % Q
        b = (t * zj) % Q
        t = np.concatenate([a, b])
    return t


def eq_point(a, b):
    """eq~(a, b) for two field points. VERIFIER: O(n)."""
    acc = 1
    for x, y in zip(a, b):
        x, y = int(x) % Q, int(y) % Q
        acc = acc * ((x * y + (Q + 1 - x) * (Q + 1 - y)) % Q) % Q
    return acc


def mle_eval(table, pt):
    """tilde{table}(pt). PROVER-side: O(2^n)."""
    return int((table * eq_table(pt) % Q).sum(dtype=np.uint64) % Q)


def int_of_point(z):
    """Int(z) = sum_j 2^{n-1-j} z_j mod Q (degree-1 extension of the
    bits-to-integer map)."""
    n = len(z)
    return sum((1 << (n - 1 - j)) % Q * (int(z[j]) % Q) for j in range(n)) % Q


def s0_closed_form(z):
    """tilde{S_0}(z) = Int(z) - 1 + prod_j(1 - z_j). VERIFIER: O(n)."""
    prod = 1
    for zj in z:
        prod = prod * ((Q + 1 - int(zj)) % Q) % Q
    return (int_of_point(z) - 1 + prod) % Q


# ----------------------------------------------------------------------
# verifier-side automaton evaluations (the structured-wiring primitives)
# ----------------------------------------------------------------------

def ge_const_eval(rv, M, n, q=Q):
    """MLE of [v >= M] at point rv. 3-state comparator (equal-so-far /
    greater / less), bits MSB-first. VERIFIER: O(n). Field-parameterised
    (q=Q default keeps every existing caller verbatim)."""
    st = [1, 0, 0]
    for j in range(n):
        mj = (M >> (n - 1 - j)) & 1
        w = [(q + 1 - int(rv[j])) % q, int(rv[j]) % q]
        new = [0, 0, 0]
        for a in (0, 1):
            if a == mj:
                new[0] += st[0] * w[a]
            elif a > mj:
                new[1] += st[0] * w[a]
            else:
                new[2] += st[0] * w[a]
            new[1] += st[1] * w[a]
            new[2] += st[2] * w[a]
        st = [x % q for x in new]
    return (st[0] + st[1]) % q


def w_div_eval(rv, ru, p, n, q=Q):
    """MLE of W(v, u) = [u = floor(v/p)] at (rv, ru). Long-division
    automaton: state = running remainder (< p); reading bit a of v
    MSB-first from state s gives t = 2s + a, emits quotient bit
    floor(t/p) and moves to t mod p. VERIFIER: O(n*p). Field-parameterised
    (q=Q default keeps every existing caller verbatim)."""
    st = [0] * p
    st[0] = 1
    for j in range(n):
        wv = [(q + 1 - int(rv[j])) % q, int(rv[j]) % q]
        wu = [(q + 1 - int(ru[j])) % q, int(ru[j]) % q]
        new = [0] * p
        for s in range(p):
            if st[s] == 0:
                continue
            for a in (0, 1):
                t = 2 * s + a
                b, s2 = divmod(t, p)
                new[s2] = (new[s2] + st[s] * wv[a] % q * wu[b]) % q
        st = new
    return sum(st) % q


# ----------------------------------------------------------------------
# generic sum-check over products of multilinear tables
# ----------------------------------------------------------------------

def sumcheck(claim, tables, terms, degree, rng, patch=False):
    """Sum-check for  claim = sum_w  sum_terms coef * prod tables[name][w].
    tables: dict name -> uint64 array (consumed). Interleaved prover
    (round polynomials by table halving) and verifier (consistency
    checks, random challenges). patch=True simulates an adaptive
    cheating prover that fixes up g(0) each round to mask a false claim.
    Returns (ok, r_list, final_claim, bound_scalars)."""
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
        if patch:
            deficit = (c - (evals[0] + evals[1])) % Q
            evals[0] = (evals[0] + deficit) % Q
        # verifier checks and challenges
        if (evals[0] + evals[1]) % Q != c:
            return False, rs, None, None
        r = int(rng.integers(0, Q))
        rs.append(r)
        c = lagrange_eval(evals, r)
        # bind the variable in all tables
        rc = (Q + 1 - r) % Q
        for nm in tables:
            tb = tables[nm]
            h = len(tb) >> 1
            tables[nm] = (tb[:h] * rc + tb[h:] * r) % Q
    scalars = {nm: int(tables[nm][0]) for nm in tables}
    return True, rs, c, scalars


# ----------------------------------------------------------------------
# the layered protocol
# ----------------------------------------------------------------------

def prover_tables(n, primes):
    """Prover precomputes all DP layers S_0..S_K (exact integers)."""
    N = 1 << n
    v = np.arange(N, dtype=np.int64)
    S = [np.maximum(v - 1, 0)]
    for i, p in enumerate(primes, start=1):
        prev = S[-1]
        b = v >= p * p
        upd = prev - b * (prev[v // p] - (i - 1))
        S.append(upd)
    return S


def verify_layer(i, p, n, z, claim, S_prev_field, rng, stats,
                 patch=False):
    """One layer reduction: claim about tilde{S_i}(z) -> claim about
    tilde{S_{i-1}} at a single new point. Returns (ok, new_z, new_claim)."""
    N = 1 << n
    v = np.arange(N, dtype=np.int64)
    c_i = (i - 1) % Q

    # ---- prover builds phase-A tables (O(2^n)) ----
    t0 = time.perf_counter()
    EQ = eq_table(z)
    Sf = S_prev_field.copy()
    Bf = (v >= p * p).astype(np.uint64)
    Gf = S_prev_field[v // p]
    stats["t_prover"] += time.perf_counter() - t0

    termsA = [(1, ["EQ", "S"]), (Q - 1, ["EQ", "B", "G"]),
              (c_i, ["EQ", "B"])]
    okA, r_v, finalA, scal = sumcheck(claim, {"EQ": EQ, "S": Sf, "B": Bf,
                                              "G": Gf}, termsA, 3, rng,
                                      patch=patch)
    if not okA:
        return False, None, None
    s1, g1 = scal["S"], scal["G"]  # prover claims

    # ---- verifier-side phase-A final check (O(n)) ----
    t0 = time.perf_counter()
    eqv = eq_point(z, r_v)
    bv = ge_const_eval(r_v, p * p, n)
    expect = eqv * ((s1 - bv * ((g1 - c_i) % Q)) % Q) % Q
    ok = (finalA % Q) == expect
    stats["t_verifier"] += time.perf_counter() - t0
    if not ok:
        return False, None, None

    # ---- phase B: verify g1 against the division wiring ----
    t0 = time.perf_counter()
    Wt = np.zeros(N, dtype=np.uint64)
    np.add.at(Wt, v // p, eq_table(r_v))   # sums of < p terms < Q each
    Wt %= Q
    stats["t_prover"] += time.perf_counter() - t0

    okB, r_u, finalB, scalB = sumcheck(g1, {"W": Wt,
                                            "S": S_prev_field.copy()},
                                       [(1, ["W", "S"])], 2, rng,
                                       patch=patch)
    if not okB:
        return False, None, None
    s2 = scalB["S"]

    t0 = time.perf_counter()
    wv = w_div_eval(r_v, r_u, p, n)
    ok = (finalB % Q) == wv * s2 % Q
    stats["t_verifier"] += time.perf_counter() - t0
    if not ok:
        return False, None, None

    # ---- line batching: s1@r_v, s2@r_u -> one claim ----
    t0 = time.perf_counter()
    line_pts = []
    for t in range(n + 1):
        tf = t % Q
        gpt = [((Q + 1 - tf) * int(a) + tf * int(b)) % Q
               for a, b in zip(r_v, r_u)]
        line_pts.append(gpt)
    hs = [mle_eval(S_prev_field, gp) for gp in line_pts]  # prover
    stats["t_prover"] += time.perf_counter() - t0

    t0 = time.perf_counter()
    if hs[0] != s1 or hs[1] != s2:   # verifier checks endpoints
        stats["t_verifier"] += time.perf_counter() - t0
        return False, None, None
    tstar = int(rng.integers(0, Q))
    new_claim = lagrange_eval(hs, tstar)
    tf = tstar % Q
    new_z = [((Q + 1 - tf) * int(a) + tf * int(b)) % Q
             for a, b in zip(r_v, r_u)]
    stats["t_verifier"] += time.perf_counter() - t0
    stats["comm"] += (4 * n) + (3 * n) + (n + 1) + 2  # field elems sent
    return True, new_z, new_claim


def run_protocol(n, rng, patch=False, corrupt_layer=None):
    """Full chain: verify pi(2^n - 1). Returns dict with outcome.
    corrupt_layer=i: prover corrupts one entry of S_i and recomputes
    all later layers consistently (a 'self-consistent liar')."""
    x = (1 << n) - 1
    primes = [p for p in range(2, int(x ** 0.5) + 1)
              if all(p % d for d in range(2, int(p ** 0.5) + 1))]
    K = len(primes)

    S = prover_tables(n, primes)
    if corrupt_layer is not None:
        i = corrupt_layer
        S[i] = S[i].copy()
        S[i][x] += 1  # off-by-one in the DP table at the output value
        v = np.arange(1 << n, dtype=np.int64)
        for j in range(i + 1, K + 1):
            p = primes[j - 1]
            b = v >= p * p
            S[j] = S[j - 1] - b * (S[j - 1][v // p] - (j - 1))

    claimed_pi = int(S[K][x])
    Sfield = [(s % Q).astype(np.uint64) for s in S]

    stats = {"t_prover": 0.0, "t_verifier": 0.0, "comm": 0}
    z = [(x >> (n - 1 - j)) & 1 for j in range(n)]  # Boolean output point
    claim = claimed_pi % Q
    for i in range(K, 0, -1):
        ok, z, claim = verify_layer(i, primes[i - 1], n, z, claim,
                                    Sfield[i - 1], rng, stats, patch=patch)
        if not ok:
            return dict(accepted=False, claimed=claimed_pi, layers=K,
                        **stats)
    # base layer: closed form, O(n)
    t0 = time.perf_counter()
    ok = s0_closed_form(z) == claim % Q
    stats["t_verifier"] += time.perf_counter() - t0
    return dict(accepted=ok, claimed=claimed_pi, layers=K, **stats)


# ----------------------------------------------------------------------
# tests / experiments
# ----------------------------------------------------------------------

def direct_pi(x):
    sieve = np.ones(x + 1, dtype=bool)
    sieve[:2] = False
    for p in range(2, int(x ** 0.5) + 1):
        if sieve[p]:
            sieve[p * p::p] = False
    return int(sieve.sum())


def selftest():
    rng = np.random.default_rng(0)
    n = 8
    # W automaton vs floor division on Boolean points
    for p in (2, 3, 7, 13):
        for _ in range(200):
            v = int(rng.integers(0, 1 << n))
            u = int(rng.integers(0, 1 << n))
            bv = [(v >> (n - 1 - j)) & 1 for j in range(n)]
            bu = [(u >> (n - 1 - j)) & 1 for j in range(n)]
            assert w_div_eval(bv, bu, p, n) == (1 if u == v // p else 0)
    # comparator vs >=
    for M in (0, 1, 37, 200, 255):
        for v in range(1 << n):
            bv = [(v >> (n - 1 - j)) & 1 for j in range(n)]
            assert ge_const_eval(bv, M, n) == (1 if v >= M else 0)
    # S0 closed form vs table MLE at random field points
    v = np.arange(1 << n, dtype=np.int64)
    S0 = (np.maximum(v - 1, 0) % Q).astype(np.uint64)
    for _ in range(20):
        pt = [int(rng.integers(0, Q)) for _ in range(n)]
        assert s0_closed_form(pt) == mle_eval(S0, pt)
    # single-layer reduction consistency at random z
    primes = [2, 3, 5, 7, 11, 13]
    S = prover_tables(n, primes)
    Sf = [(s % Q).astype(np.uint64) for s in S]
    i, p = 4, primes[3]
    pt = [int(rng.integers(0, Q)) for _ in range(n)]
    claim = mle_eval(Sf[i], pt)
    stats = {"t_prover": 0, "t_verifier": 0, "comm": 0}
    ok, z2, c2 = verify_layer(i, p, n, pt, claim, Sf[i - 1],
                              np.random.default_rng(1), stats)
    assert ok and c2 == mle_eval(Sf[i - 1], z2)
    # tiny end-to-end
    res = run_protocol(10, np.random.default_rng(2))
    assert res["accepted"] and res["claimed"] == direct_pi((1 << 10) - 1)
    print("selftest OK")


def main(n, cheat_trials, seed):
    x = (1 << n) - 1
    rng = np.random.default_rng(seed)

    res = run_protocol(n, rng)
    truth = direct_pi(x)
    print(f"x = 2^{n}-1 = {x}, layers (primes <= sqrt x): {res['layers']}")
    print(f"honest run: {'ACCEPTED' if res['accepted'] else 'REJECTED'}, "
          f"claimed pi(x) = {res['claimed']}, direct sieve pi(x) = {truth}, "
          f"match = {res['claimed'] == truth}")
    assert res["accepted"] and res["claimed"] == truth
    print(f"t_prover = {res['t_prover']:.2f}s   "
          f"t_verifier = {res['t_verifier']*1000:.1f}ms   "
          f"communication = {res['comm']} field elems "
          f"(~{res['comm']*4/1024:.1f} KB)")

    rej_patch = 0
    for t in range(cheat_trials):
        r = run_protocol(n, np.random.default_rng(seed + 100 + t),
                         patch=True, corrupt_layer=res["layers"] // 2)
        rej_patch += (not r["accepted"])
    print(f"self-consistent liar + adaptive round patching: "
          f"rejected {rej_patch}/{cheat_trials}")

    rej_cor = 0
    for t in range(cheat_trials):
        r = run_protocol(n, np.random.default_rng(seed + 500 + t),
                         corrupt_layer=res["layers"] // 2)
        rej_cor += (not r["accepted"])
    print(f"self-consistent liar (honest protocol on corrupted DP): "
          f"rejected {rej_cor}/{cheat_trials}")
    nK = n * res["layers"]
    print(f"soundness bound ~ 7nK/q = {7*nK/Q:.2e}")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, default=16)
    ap.add_argument("--cheat-trials", type=int, default=25)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--selftest", action="store_true")
    args = ap.parse_args()
    if args.selftest:
        selftest()
    else:
        main(args.n, args.cheat_trials, args.seed)
