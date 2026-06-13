#!/usr/bin/env python3
"""
sumcheck_wheel_verification.py — S491 proof-of-concept for Thread 12 (P12):
succinct verification of pi(x)-adjacent counts.

Implements an unconditionally-sound interactive proof (LFKN 1992 sum-check)
for the partial-sieve count

    Phi(2^n, P) = #{ 0 <= w < 2^n : p does not divide w, for all p in P }

where P is a set of odd "wheel" primes (w = 0 is divisible by every p, so
it never counts; the count equals the number of integers in [1, 2^n)
coprime to prod(P)).

The VERIFIER runs in O(n * sum_{p in P} p) field operations — polylog in
x = 2^n for any polylog-sized wheel — while the PROVER does the Theta(2^n)
work. Soundness is information-theoretic (no hardness assumptions): a
cheating prover passes with probability <= n*|P|/q over the verifier's
random field challenges.

Key primitive (the project-new part): the multilinear extension (MLE) of
the divisibility indicator delta_p(w) = [w == 0 mod p] is evaluable at an
ARBITRARY field point in O(n*p) time by a transfer-matrix DP over residue
states, because delta_p is recognised by a p-state automaton reading bits
MSB-first. This circumvents the S48 closure (the MLE of the PRIMALITY
indicator chi_P is dense, with 2^{n/2}-size coefficients — true for
primality, NOT for fixed-modulus divisibility, which is an automaton).

Field: F_q with q = 2^31 - 1 (Mersenne prime); numpy uint64 products stay
below 2^62 so all arithmetic is exact. Soundness error <= n*|P|/q ~ 1e-7
at demo sizes (a production protocol would use q ~ 2^128 via extension
fields; orthogonal to the structure demonstrated here).

What would falsify this construction: an accepted run whose output count
differs from the direct sieve count (completeness bug), or a cheating
prover acceptance rate materially above n*k/q (soundness bug), or verifier
cost growing like 2^n (the whole point is that it must not).

Usage:
  python3 sumcheck_wheel_verification.py --selftest
  python3 sumcheck_wheel_verification.py --n 20 --primes 3,5,7,11,13 --trials 100
  python3 sumcheck_wheel_verification.py --scaling
"""

import argparse
import time

import numpy as np

Q = (1 << 31) - 1  # Mersenne prime 2147483647


# ----------------------------------------------------------------------
# MLE of delta_p via residue-automaton transfer DP
# ----------------------------------------------------------------------

def automaton_step(v, z, p):
    """One DP step: feed coordinate value z (an element of F_q) into the
    p-state residue automaton state vector v (length p, entries mod Q).
    Transition s -> 2s+b mod p on bit b; the MLE weights are (1-z, z).
    Both transition maps are bijections on Z_p (p odd), so each target
    state receives exactly two contributions: values < 2Q < 2^32, exact
    in uint64."""
    z = int(z) % Q
    zc = (Q + 1 - z) % Q
    s = np.arange(p, dtype=np.int64)
    t0 = (2 * s) % p
    t1 = (2 * s + 1) % p
    nv = np.zeros(p, dtype=np.uint64)
    nv[t0] = (v * zc) % Q
    nv[t1] = (nv[t1] + (v * z) % Q) % Q
    return nv


def prefix_state(prefix_vals, p):
    """Residue-state vector after feeding the (possibly non-Boolean)
    prefix into the p-automaton MLE DP. O(len(prefix) * p)."""
    v = np.zeros(p, dtype=np.uint64)
    v[0] = 1
    for z in prefix_vals:
        v = automaton_step(v, z, p)
    return v


def delta_mle_eval(point, p):
    """Evaluate the MLE of delta_p (bits MSB-first) at an arbitrary field
    point. O(n*p). This is the verifier's final-check primitive."""
    return int(prefix_state(point, p)[0])


def g_eval(point, primes):
    """g(z) = prod_{p in P} (1 - delta~_p(z)) mod Q at an arbitrary point.
    Per-variable degree |P|. Restricted to Boolean inputs, equals the
    coprime-to-wheel indicator."""
    acc = 1
    for p in primes:
        acc = (acc * ((Q + 1 - delta_mle_eval(point, p)) % Q)) % Q
    return acc


# ----------------------------------------------------------------------
# Prover
# ----------------------------------------------------------------------

def prover_round_poly(r_prefix, n, primes):
    """Round j = len(r_prefix)+1 message: evaluations [g_j(0..k)] of the
    univariate g_j(X) = sum over Boolean suffixes of
    g(r_1..r_{j-1}, X, suffix). Degree <= k = |P|.

    Suffix trick: with a Boolean suffix of integer value t (m = n-j bits),
    the automaton from state s ends at (s*2^m + t) mod p, so
    delta~ = u[s*] with the single state s* = -t * inv(2^m) mod p.
    Cost O(2^m * k) per evaluation point via numpy gather."""
    j = len(r_prefix) + 1
    m = n - j
    k = len(primes)
    nsuf = 1 << m
    evals = []
    # suffix value mod p, per prime (recomputed per round; geometric total)
    t_mod = {p: np.arange(nsuf, dtype=np.int64) % p for p in primes}
    sstar = {p: ((-t_mod[p]) % p * pow(2, -m, p)) % p for p in primes}
    for X in range(k + 1):
        pref = list(r_prefix) + [X]
        prod = np.ones(nsuf, dtype=np.uint64)
        for p in primes:
            u = prefix_state(pref, p)
            delta = u[sstar[p]]
            prod = (prod * ((Q + 1 - delta) % Q)) % Q
        # entries < Q < 2^31; sum of 2^m of them fits uint64 for m <= 33
        evals.append(int(prod.sum(dtype=np.uint64) % Q))
    return evals


# ----------------------------------------------------------------------
# Verifier helpers
# ----------------------------------------------------------------------

def lagrange_eval(ys, r):
    """Evaluate the degree-<=k polynomial through (i, ys[i]), i = 0..k,
    at field point r. O(k^2) — k is the wheel size, tiny."""
    k = len(ys) - 1
    r = int(r) % Q
    total = 0
    for i in range(k + 1):
        num, den = 1, 1
        for jj in range(k + 1):
            if jj == i:
                continue
            num = (num * ((r - jj) % Q)) % Q
            den = (den * ((i - jj) % Q)) % Q
        total = (total + ys[i] * num % Q * pow(den, Q - 2, Q)) % Q
    return total


# ----------------------------------------------------------------------
# Protocol driver
# ----------------------------------------------------------------------

def run_protocol(n, primes, rng, cheat_offset=0):
    """Run the full sum-check. The prover's claim is the honest sum plus
    cheat_offset (mod Q). With cheat_offset != 0 the prover additionally
    patches each round polynomial minimally (additive fix to g_j(0)) so
    every running consistency check passes; the lie must then survive the
    verifier's final independent evaluation of g at the random point —
    probability <= n*k/q.

    Returns dict with: accepted, claimed, t_prover, t_verifier, comm_field_elems.
    """
    t_p = 0.0
    t_v = 0.0

    # Round 1 prover message also determines the honest claim.
    tic = time.perf_counter()
    evals = prover_round_poly([], n, primes)
    t_p += time.perf_counter() - tic
    honest_S = (evals[0] + evals[1]) % Q
    claimed = (honest_S + cheat_offset) % Q

    r_prefix = []
    c = claimed
    comm = 0
    for j in range(1, n + 1):
        if j > 1:
            tic = time.perf_counter()
            evals = prover_round_poly(r_prefix, n, primes)
            t_p += time.perf_counter() - tic
        if cheat_offset:
            deficit = (c - (evals[0] + evals[1])) % Q
            evals = list(evals)
            evals[0] = (evals[0] + deficit) % Q
        comm += len(evals)

        tic = time.perf_counter()
        ok = (evals[0] + evals[1]) % Q == c
        if not ok:
            t_v += time.perf_counter() - tic
            return dict(accepted=False, claimed=claimed, t_prover=t_p,
                        t_verifier=t_v, comm_field_elems=comm)
        r_j = int(rng.integers(0, Q))
        c = lagrange_eval(evals, r_j)
        t_v += time.perf_counter() - tic
        r_prefix.append(r_j)

    tic = time.perf_counter()
    final_ok = g_eval(r_prefix, primes) == c
    t_v += time.perf_counter() - tic
    return dict(accepted=final_ok, claimed=claimed, t_prover=t_p,
                t_verifier=t_v, comm_field_elems=comm)


def direct_count(n, primes):
    """Ground truth by sieve: integers in [0, 2^n) with no factor in P
    (0 is excluded automatically as a multiple of every p)."""
    N = 1 << n
    mask = np.ones(N, dtype=bool)
    for p in primes:
        mask[::p] = False
    return int(mask.sum())


# ----------------------------------------------------------------------
# Self-test and experiments
# ----------------------------------------------------------------------

def selftest():
    rng = np.random.default_rng(0)
    # 1. MLE agrees with the Boolean indicator on the cube.
    n = 8
    for p in (3, 5, 7):
        for w in range(1 << n):
            bits = [(w >> (n - 1 - i)) & 1 for i in range(n)]
            assert delta_mle_eval(bits, p) == (1 if w % p == 0 else 0), (p, w)
    # 2. Multilinearity: affine in each coordinate.
    for p in (3, 7):
        pt = [int(rng.integers(0, Q)) for _ in range(n)]
        for i in range(n):
            a, b, m = pt.copy(), pt.copy(), pt.copy()
            a[i], b[i] = 0, 1
            m[i] = int(rng.integers(0, Q))
            va, vb = delta_mle_eval(a, p), delta_mle_eval(b, p)
            expect = (va * ((Q + 1 - m[i]) % Q) + vb * m[i]) % Q
            assert delta_mle_eval(m, p) == expect, (p, i)
    # 3. Round-poly sums match brute force at small n.
    n, primes = 6, [3, 5]
    evals = prover_round_poly([], n, primes)
    s = (evals[0] + evals[1]) % Q
    assert s == direct_count(n, primes), (s, direct_count(n, primes))
    print("selftest OK")


def experiment(n, primes, trials, seed, skip_direct):
    rng = np.random.default_rng(seed)
    k = len(primes)
    print(f"n={n} (x = 2^n = {1 << n}), wheel P={primes} (k={k}), field q=2^31-1")

    # Honest run + completeness.
    res = run_protocol(n, primes, rng)
    assert res["accepted"], "honest prover was rejected — completeness bug"
    print(f"honest run: ACCEPTED, claimed Phi = {res['claimed']}")
    if not skip_direct:
        tic = time.perf_counter()
        truth = direct_count(n, primes)
        t_direct = time.perf_counter() - tic
        assert res["claimed"] == truth % Q, (res["claimed"], truth)
        print(f"direct sieve count  = {truth}  (matches; t_direct = {t_direct:.3f}s)")
    print(f"t_prover   = {res['t_prover']:.3f}s")
    print(f"t_verifier = {res['t_verifier']*1000:.2f}ms   "
          f"(communication: {res['comm_field_elems']} field elems = "
          f"{res['comm_field_elems']*4} bytes)")

    # Soundness: adaptive additive cheater claiming Phi+1.
    rejected = 0
    for t in range(trials):
        r = run_protocol(n, primes, np.random.default_rng(seed + 1000 + t),
                         cheat_offset=1)
        rejected += (not r["accepted"])
    bound = n * k / Q
    print(f"cheating prover (claim+1, adaptive patch): rejected {rejected}/{trials} "
          f"(soundness error bound n*k/q = {bound:.2e})")
    return res


def scaling(seed):
    print("=== scaling sweep: verifier vs prover vs direct ===")
    rows = []
    for n, primes in [(12, [3, 5, 7]),
                      (16, [3, 5, 7, 11, 13]),
                      (20, [3, 5, 7, 11, 13]),
                      (22, [3, 5, 7, 11, 13, 17, 19, 23, 29, 31])]:
        rng = np.random.default_rng(seed)
        res = run_protocol(n, primes, rng)
        tic = time.perf_counter()
        truth = direct_count(n, primes)
        t_direct = time.perf_counter() - tic
        assert res["accepted"] and res["claimed"] == truth % Q
        rows.append((n, len(primes), res["t_prover"], res["t_verifier"],
                     t_direct, res["comm_field_elems"]))
        print(f"n={n:2d} k={len(primes):2d}  prover={res['t_prover']:8.3f}s  "
              f"verifier={res['t_verifier']*1000:7.2f}ms  "
              f"direct={t_direct:7.3f}s  comm={res['comm_field_elems']*4}B")
    return rows


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, default=20)
    ap.add_argument("--primes", type=str, default="3,5,7,11,13")
    ap.add_argument("--trials", type=int, default=100)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--skip-direct", action="store_true")
    ap.add_argument("--selftest", action="store_true")
    ap.add_argument("--scaling", action="store_true")
    args = ap.parse_args()

    if args.selftest:
        selftest()
    elif args.scaling:
        scaling(args.seed)
    else:
        primes = [int(t) for t in args.primes.split(",")]
        assert all(p % 2 == 1 for p in primes), "wheel primes must be odd (2^m must be invertible mod p)"
        experiment(args.n, primes, args.trials, args.seed, args.skip_direct)
