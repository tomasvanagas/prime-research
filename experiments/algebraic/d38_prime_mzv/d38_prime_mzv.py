"""
D38 — Prime-restricted multiple zeta values vs Brown 2012 motivic Galois algebra.

Defines:
    M_s   := P(s)   = Σ_{p prime} 1/p^s              (prime zeta function)
    ζ_P(s,t) := Σ_{p<q prime} 1/(p^s q^t)            (prime-restricted MZV depth 2)

Symmetric identity (proven, sanity check):
    ζ_P(s,t) + ζ_P(t,s) = M_s · M_t − M_{s+t}        (Euler reflection over primes)
    ζ_P(s,s) = (M_s² − M_{2s}) / 2

Antisymmetric content:
    A(s,t) := ζ_P(s,t) − ζ_P(t,s),    A(s,t) = −A(t,s),  A(s,s) = 0

Question (D38): is A(s,t) a Q-linear combination of weight-(s+t) elements of
    BASIS = { ζ(weight-w MZVs), products M_a · M_b · ..., M_a · ζ(...), 1 }
where standard MZV closed forms are taken as known (Hoffman basis ⊕ products)?

PSLQ-test A(s,t) numerically against this basis. If integer relation found
at high precision and survives prime-cutoff cross-validation, the prime-MZV
algebra reduces to Brown's MZV algebra ⊕ {Mertens constants} (E). If PSLQ
returns norm > 10^something at full precision, A(s,t) generates a strictly
new period — primes carry motivic content beyond Brown's framework (I).

Numerical strategy: high-precision (60 dps) computation with N up to 10^6 or
10^7, with the analytic correction
    ζ_P(s,t) = inner_sum_{p<q≤N} 1/(p^s q^t)
             + prefix_s(N) · tail_t(N)
             + tail_above
where tail_above is bounded by tail_s(N) · tail_t(N) ~ 1/(N^{s+t-2} · log² N).

Run:
  python d38_prime_mzv.py [--N 1000000] [--dps 50] [--weights 5,6,7,8]
"""
import argparse
import json
import time
from pathlib import Path

import numpy as np
from mpmath import mp, mpf, zeta, pslq, log, primezeta, fabs


def primes_up_to(N):
    """Sieve of Eratosthenes; return list of primes p ≤ N."""
    if N < 2:
        return []
    sieve = np.ones(N + 1, dtype=bool)
    sieve[:2] = False
    for i in range(2, int(N**0.5) + 1):
        if sieve[i]:
            sieve[i * i :: i] = False
    return np.where(sieve)[0].tolist()


def prefix_M(primes, s, dps):
    """Compute Σ_{p in primes} 1/p^s exactly to current dps."""
    mp.dps = dps
    total = mpf(0)
    for p in primes:
        total += mpf(1) / mpf(p) ** s
    return total


def prime_zeta_high(s, dps):
    """M_s = P(s) = Σ_p 1/p^s to current dps via mpmath.primezeta."""
    mp.dps = dps
    return primezeta(s)


def zeta_P_double(primes, s, t, M_s, M_t, dps, return_components=False):
    """Compute ζ_P(s, t) = Σ_{p<q ≤ N} 1/(p^s q^t) + prefix_s · tail_t.

    Truncation error is bounded by tail_s · tail_t (when both s, t ≥ 2).
    """
    mp.dps = dps
    inner = mpf(0)
    # Streaming: for each prime p, accumulate partial 1/q^t over q > p.
    # Work in two passes: (a) prefix_t_partial[i] = Σ_{j<i} 1/p_j^t, (b)
    # then for each p_i compute tail = M_t_finite - prefix_t_partial[i] - 1/p_i^t,
    # i.e., Σ_{q>p_i, q ≤ N} 1/q^t.
    # Build prefix_t_partial in mpf:
    inv_p_t = [mpf(1) / mpf(p) ** t for p in primes]
    M_t_finite = mpf(0)
    prefix_t_part = [mpf(0)] * len(primes)
    for i, val in enumerate(inv_p_t):
        prefix_t_part[i] = M_t_finite
        M_t_finite += val
    # M_t_finite now = prefix sum to all primes.
    for i, p in enumerate(primes):
        # Σ_{q>p, q ≤ N} 1/q^t = M_t_finite - prefix_t_part[i] - 1/p^t
        inner_sum_q = M_t_finite - prefix_t_part[i] - inv_p_t[i]
        inner += (mpf(1) / mpf(p) ** s) * inner_sum_q

    # Compute prefix_s and tail_t analytic correction
    prefix_s_finite = mpf(0)
    for p in primes:
        prefix_s_finite += mpf(1) / mpf(p) ** s
    tail_t = M_t - M_t_finite
    correction = prefix_s_finite * tail_t

    val = inner + correction

    if return_components:
        # Bound on tail_above (the q > N, p > N part)
        tail_s = M_s - prefix_s_finite
        tail_above_bound = tail_s * tail_t
        return val, {
            "inner": inner,
            "correction": correction,
            "tail_above_bound": tail_above_bound,
            "prefix_s": prefix_s_finite,
            "M_t_finite": M_t_finite,
            "tail_s": tail_s,
            "tail_t": tail_t,
        }
    return val


def run_zeta_P_panel(primes_list, weights, dps, M_cache):
    """Compute ζ_P(s,t) for all (s,t) with s,t ≥ 2 and s+t ∈ weights."""
    mp.dps = dps
    panel = {}
    for w in weights:
        for s in range(2, w):
            t = w - s
            if t < 2:
                continue
            key = f"({s},{t})"
            t0 = time.time()
            val, info = zeta_P_double(
                primes_list, s, t, M_cache[s], M_cache[t], dps,
                return_components=True,
            )
            elapsed = time.time() - t0
            panel[key] = {
                "s": s, "t": t, "weight": w,
                "value": str(val),
                "tail_above_bound": str(info["tail_above_bound"]),
                "elapsed_sec": elapsed,
            }
            print(
                f"  ζ_P({s},{t}) = {mp.nstr(val, 25)}   "
                f"err≤{mp.nstr(info['tail_above_bound'], 3)}   ({elapsed:.1f}s)"
            )
    return panel


def verify_symmetric_identity(panel, M_cache, weights, dps):
    """Sanity check: ζ_P(s,t) + ζ_P(t,s) =? M_s · M_t − M_{s+t}."""
    mp.dps = dps
    print("\nSanity: Euler reflection ζ_P(s,t) + ζ_P(t,s) = M_s M_t − M_{s+t}")
    out = {}
    for w in weights:
        for s in range(2, w):
            t = w - s
            if t < 2 or s >= t:
                continue
            k1, k2 = f"({s},{t})", f"({t},{s})"
            if k1 not in panel or k2 not in panel:
                continue
            v1 = mpf(panel[k1]["value"])
            v2 = mpf(panel[k2]["value"])
            lhs = v1 + v2
            rhs = M_cache[s] * M_cache[t] - M_cache[s + t]
            diff = fabs(lhs - rhs)
            out[f"{s},{t}"] = {"diff": str(diff)}
            print(
                f"  s={s},t={t}: LHS={mp.nstr(lhs, 20)}  RHS={mp.nstr(rhs, 20)}  diff={mp.nstr(diff, 3)}"
            )
    return out


def _enumerate_compositions(weight, min_part=2):
    """Yield each composition (multiset) of `weight` into parts ≥ min_part as
    a sorted tuple. E.g. weight=5 → [(2,3), (5,)]."""
    def rec(remaining, lower_bound):
        if remaining == 0:
            yield ()
            return
        if remaining < lower_bound:
            return
        for k in range(lower_bound, remaining + 1):
            for tail in rec(remaining - k, k):
                yield (k,) + tail
    yield from rec(weight, min_part)


def build_basis(weight, M_cache, dps):
    """Construct weight-w basis: ALL monomials of total weight w in the
    generators
        ζ(2), ζ(3), ζ(5), ζ(7), ...   (weight = argument)
        M_2, M_3, M_4, M_5, ...        (weight = subscript)

    For each multiset partition of `weight` into parts ≥ 2, label each part
    as either ζ or M and form the product. Even ζ(2k) (k≥2) are
    intentionally OMITTED because ζ(2k) ∈ Q · ζ(2)^k → linearly dependent
    (PSLQ would find a trivial relation). ζ(2) and ζ(odd) for odd ≥ 3 are
    kept as algebraic generators.
    """
    mp.dps = dps
    basis = []
    # Generators of each weight k ≥ 2:
    #   weight 2: ζ(2), M_2  (2 generators)
    #   weight 3: ζ(3), M_3
    #   weight 4: M_4  (ζ(4) reduces)
    #   weight k odd ≥ 3: ζ(k), M_k
    #   weight k even ≥ 4: M_k only
    def gens_of_weight(k):
        if k == 2:
            return [("zeta(2)", zeta(2)), (f"M_{k}", M_cache[k])]
        elif k % 2 == 1 and k >= 3:
            return [(f"zeta({k})", zeta(k)), (f"M_{k}", M_cache[k])]
        else:
            return [(f"M_{k}", M_cache[k])]

    seen_labels = set()
    for parts in _enumerate_compositions(weight, min_part=2):
        # Build a list of (gen_lists) per part, then take all combinations.
        per_part = [gens_of_weight(p) for p in parts]
        # Collapse identical part-positions to avoid double-counting equal
        # multisets: enforce non-decreasing label index within equal-part
        # blocks. For simplicity here we just dedupe via label sort.
        from itertools import product as iprod
        # Group equal parts together in `parts` (already sorted ascending)
        # and within each block of equal parts, enumerate multisets of choices
        # rather than ordered tuples.
        blocks = []
        i = 0
        while i < len(parts):
            j = i
            while j < len(parts) and parts[j] == parts[i]:
                j += 1
            blocks.append((parts[i], j - i))
            i = j
        # For each block (part_value, count), enumerate multisets of size `count`
        # from gens_of_weight(part_value)
        from itertools import combinations_with_replacement as cwr
        block_choices = []
        for pv, cnt in blocks:
            gens = gens_of_weight(pv)
            choices = list(cwr(range(len(gens)), cnt))
            block_choices.append((pv, gens, choices))
        for combo in iprod(*[bc[2] for bc in block_choices]):
            # combo is a tuple of multiset-choice tuples per block
            label_parts = []
            value = mpf(1)
            for (pv, gens, _), choice in zip(block_choices, combo):
                # `choice` is a tuple of indices into `gens`
                for idx in choice:
                    label_parts.append(gens[idx][0])
                    value = value * gens[idx][1]
            label = "*".join(label_parts)
            # Canonicalise label by sorting factors lexicographically
            canon = "*".join(sorted(label.split("*")))
            if canon not in seen_labels:
                seen_labels.add(canon)
                basis.append((canon, value))
    return basis


def pslq_test_target(target, basis, dps, tol_exp=None):
    """PSLQ test whether target ∈ Q-span(basis). Return relation or None.

    If tol_exp is given, declare a relation only if the dot-product of
    integer coefficients · values equals target within tol_exp digits.
    """
    mp.dps = dps
    vals = [target] + [b[1] for b in basis]
    labels = ["TARGET"] + [b[0] for b in basis]
    # Default tol: half the working precision
    if tol_exp is None:
        tol_exp = -int(0.5 * dps)
    try:
        rel = pslq(vals, tol=mpf(10) ** tol_exp, maxcoeff=10**18)
    except Exception as e:
        return {"status": "FAIL", "msg": str(e)}
    if rel is None:
        return {"status": "NO_RELATION", "tol_exp": tol_exp}
    # Check the residual: |Σ rel_i · vals_i|
    residual = mpf(0)
    for r, v in zip(rel, vals):
        residual += mpf(int(r)) * v
    res_str = mp.nstr(fabs(residual), 5)
    out = {
        "status": "RELATION_FOUND",
        "coeffs": [int(r) for r in rel],
        "labels": labels,
        "residual_abs": res_str,
        "tol_exp": tol_exp,
    }
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--N", type=int, default=1000000)
    ap.add_argument("--dps", type=int, default=50)
    ap.add_argument("--weights", default="5,6,7,8")
    ap.add_argument("--out", default="d38_results.json")
    args = ap.parse_args()
    weights = [int(w) for w in args.weights.split(",")]
    dps = args.dps
    mp.dps = dps

    print(f"=== D38 prime-restricted MZV — N={args.N}, dps={dps} ===")
    t0 = time.time()
    primes_list = primes_up_to(args.N)
    print(f"  {len(primes_list)} primes ≤ {args.N}  ({time.time() - t0:.1f}s)")

    # Pre-cache M_s for all needed s
    needed_s = set()
    for w in weights:
        for s in range(2, w):
            t = w - s
            if t >= 2:
                needed_s.add(s)
                needed_s.add(t)
                needed_s.add(s + t)
    needed_s = sorted(needed_s)
    M_cache = {}
    print("  Computing prime-zeta M_s to high precision via mpmath.primezeta...")
    for s in needed_s:
        t1 = time.time()
        M_cache[s] = prime_zeta_high(s, dps)
        print(f"    M_{s} = {mp.nstr(M_cache[s], 30)}   ({time.time()-t1:.2f}s)")

    # Compute ζ_P(s, t)
    print("\n  Computing ζ_P(s, t) for all (s,t), s,t≥2, s+t ∈ weights...")
    panel = run_zeta_P_panel(primes_list, weights, dps, M_cache)

    # Sanity: Euler reflection
    sym_check = verify_symmetric_identity(panel, M_cache, weights, dps)

    # Antisymmetric A(s,t) and PSLQ tests
    print("\n  PSLQ tests on antisymmetric A(s,t) = ζ_P(s,t) − ζ_P(t,s):")
    pslq_results = {}
    for w in weights:
        for s in range(2, w):
            t = w - s
            if t < 2 or s >= t:
                continue
            k1, k2 = f"({s},{t})", f"({t},{s})"
            if k1 not in panel or k2 not in panel:
                continue
            v1 = mpf(panel[k1]["value"])
            v2 = mpf(panel[k2]["value"])
            A = v1 - v2
            print(
                f"\n  --- A({s},{t}) = {mp.nstr(A, 25)} (weight {w}) ---"
            )
            basis = build_basis(w, M_cache, dps)
            print(f"    basis size: {len(basis)} items")
            for lbl, val in basis:
                pass
            res = pslq_test_target(A, basis, dps, tol_exp=-int(0.4 * dps))
            print(f"    PSLQ status: {res.get('status')}")
            if res.get("status") == "RELATION_FOUND":
                # pretty-print the relation: TARGET coefficient is the first
                # entry of res['coeffs']
                tc = res["coeffs"][0]
                if tc != 0:
                    sgn = -1 if tc > 0 else 1
                    print(f"    A({s},{t}) ≈ ", end="")
                    parts = []
                    for r, lbl in zip(res["coeffs"][1:], res["labels"][1:]):
                        if r != 0:
                            num = sgn * r
                            parts.append(f"({num}/{abs(tc)}) · {lbl}")
                    print(" + ".join(parts))
                print(f"    residual |Σ| = {res['residual_abs']}")
            pslq_results[f"A({s},{t})"] = {
                "value": str(A),
                "weight": w,
                "pslq": res,
            }

    out = {
        "args": vars(args),
        "needed_s": needed_s,
        "M_cache": {s: str(v) for s, v in M_cache.items()},
        "panel": panel,
        "symmetric_check": sym_check,
        "pslq_results": pslq_results,
    }
    Path(args.out).write_text(json.dumps(out, indent=2, default=str))
    print(f"\nWrote {args.out}")


if __name__ == "__main__":
    main()
