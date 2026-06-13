"""Wheel-coprime Liouville-Möbius identity verification.

Tests Theorems 1-3 of `definition.md` pointwise on `x ∈ [0, N]` for a
range of W (squarefree primorial, squarefree non-primorial, and
non-squarefree).

Composition: E1.6 (bisection π = (x−L)/2 − C_3) lifted to wheel-graded
form via the Möbius / completely-multiplicative inversion of Liouville
restricted to coprime-to-W subsets.

Output: F1-F6 PASS/FAIL verdicts, JSON results.
"""

import argparse
import json
import math
import time
from pathlib import Path


def sieve(N):
    """Return:
    - lpf[1..N] : smallest prime factor (lpf[1] = 1)
    - is_prime[1..N] : 0/1 prime indicator (is_prime[1] = 0)
    - bigomega[1..N] : Ω(n) = number of prime factors with multiplicity
    - liouville[1..N] : λ(n) = (-1)^Ω(n)
    - mu[1..N] : μ(n)
    """
    lpf = [0] * (N + 1)
    lpf[1] = 1
    is_prime = [0] * (N + 1)
    primes = []
    for i in range(2, N + 1):
        if lpf[i] == 0:
            lpf[i] = i
            is_prime[i] = 1
            primes.append(i)
        for p in primes:
            if p > lpf[i] or i * p > N:
                break
            lpf[i * p] = p
    # Ω, λ, μ
    bigomega = [0] * (N + 1)
    bigomega[1] = 0
    mu = [0] * (N + 1)
    mu[1] = 1
    for n in range(2, N + 1):
        p = lpf[n]
        m = n // p
        bigomega[n] = bigomega[m] + 1
        if m % p == 0:
            mu[n] = 0
        else:
            mu[n] = -mu[m]
    liouville = [1 if (bigomega[n] % 2 == 0) else -1 for n in range(N + 1)]
    return lpf, is_prime, bigomega, liouville, mu


def cum(arr):
    """In-place running prefix sum returning new list."""
    out = [0] * len(arr)
    s = 0
    for i, v in enumerate(arr):
        s += v
        out[i] = s
    return out


def L_partial(liouville_cum, x):
    """L(x) = Σ_{n ≤ x} λ(n)."""
    if x < 0:
        return 0
    if x >= len(liouville_cum):
        x = len(liouville_cum) - 1
    return liouville_cum[x]


def divisors(n):
    """Sorted divisors of n via trial division."""
    out = []
    i = 1
    while i * i <= n:
        if n % i == 0:
            out.append(i)
            if i != n // i:
                out.append(n // i)
        i += 1
    return sorted(out)


def radical(n):
    """rad(n) = product of distinct primes dividing n."""
    if n == 1:
        return 1
    r = 1
    p = 2
    m = n
    while p * p <= m:
        if m % p == 0:
            r *= p
            while m % p == 0:
                m //= p
        p += 1
    if m > 1:
        r *= m
    return r


def omega(n):
    """ω(n) = number of distinct prime divisors."""
    return len([d for d in divisors(radical(n)) if d > 1 and is_prime_int(d)])


def is_prime_int(n):
    if n < 2:
        return False
    if n < 4:
        return True
    if n % 2 == 0:
        return False
    p = 3
    while p * p <= n:
        if n % p == 0:
            return False
        p += 2
    return True


def L_W_direct(liouville, W, x):
    """Direct computation of L_W(x) by enumeration."""
    s = 0
    for n in range(1, x + 1):
        if math.gcd(n, W) == 1:
            s += liouville[n]
    return s


def L_W_via_identity(liouville_cum, W, x):
    """L_W(x) via Theorem 1: Σ_{d | rad(W)} L(⌊x/d⌋)."""
    r = radical(W)
    s = 0
    calls = 0
    for d in divisors(r):
        s += L_partial(liouville_cum, x // d)
        calls += 1
    return s, calls


def n_W_via_mobius(W, x, mu):
    """n_W(x) = Σ_{d | rad(W)} μ(d) ⌊x/d⌋."""
    r = radical(W)
    s = 0
    for d in divisors(r):
        s += mu[d] * (x // d)
    return s


def n_W_direct(W, x):
    """Direct count of n ≤ x with gcd(n, W) = 1."""
    return sum(1 for n in range(1, x + 1) if math.gcd(n, W) == 1)


def A_W_direct(bigomega, W, x):
    """A_W(x) = #{n ≤ x : gcd(n, W) = 1, Ω(n) odd}."""
    return sum(
        1
        for n in range(1, x + 1)
        if math.gcd(n, W) == 1 and bigomega[n] % 2 == 1
    )


def pi_W_direct(is_prime, W, x):
    """π_W(x) = #{p ≤ x : p prime, gcd(p, W) = 1}."""
    return sum(
        1 for n in range(2, x + 1) if is_prime[n] and math.gcd(n, W) == 1
    )


def C3_W_direct(bigomega, is_prime, W, x):
    """C_{3,W}(x) = #{n ≤ x : gcd(n,W)=1, composite, Ω(n) odd}."""
    return sum(
        1
        for n in range(2, x + 1)
        if math.gcd(n, W) == 1 and (not is_prime[n]) and bigomega[n] % 2 == 1
    )


def main(N=10_000, out_path=None):
    print(f"Sieving up to N = {N} ...")
    t0 = time.time()
    lpf, is_prime, bigomega, liouville, mu = sieve(N)
    liouville_cum = cum(liouville[1:])
    # liouville_cum[i] = L(i+1); make 1-indexed accessor:
    #   L(0) = 0; L(x) = liouville_cum[x-1] for x ≥ 1.
    L_arr = [0] * (N + 1)
    for x in range(1, N + 1):
        L_arr[x] = liouville_cum[x - 1]
    print(f"  Sieve done in {time.time() - t0:.2f}s.")

    # Test grid for W
    W_list = [2, 3, 5, 6, 10, 12, 15, 30, 60, 105, 210, 420, 2310]
    W_nonsqfree = [4, 9, 12, 60, 420, 2700]
    W_all = sorted(set(W_list + W_nonsqfree))

    x_grid_F1 = list(range(0, N + 1, max(1, N // 200)))  # ~200 points
    if x_grid_F1[-1] != N:
        x_grid_F1.append(N)

    x_grid_F4 = [N // 4, N // 2, 3 * N // 4, N]

    results = {
        "N": N,
        "F1": {},
        "F2": {},
        "F3": {},
        "F4": {},
        "F5": {},
        "F6": {},
    }

    # ---------------- F1 + F3 + F5 (pointwise) ----------------
    print("\nF1 + F3 + F5 — pointwise Liouville identity, parity, divisor count:")
    for W in W_all:
        rW = radical(W)
        n_div = len(divisors(rW))
        ow = sum(1 for p in [2, 3, 5, 7, 11, 13, 17, 19, 23] if rW % p == 0)
        # Predicted divisor-count: 2^ω(W) iff ω(rW) = number of distinct prime
        # factors. Compute exactly.
        ow_exact = 0
        m = rW
        p = 2
        while p * p <= m:
            if m % p == 0:
                ow_exact += 1
                while m % p == 0:
                    m //= p
            p += 1
        if m > 1:
            ow_exact += 1

        max_abs_F1 = 0
        max_abs_F3 = 0
        calls_used = []
        for x in x_grid_F1:
            ldirect = L_W_direct(liouville, W, x)
            lid, calls = L_W_via_identity(L_arr, W, x)
            calls_used.append(calls)
            diff = ldirect - lid
            if abs(diff) > max_abs_F1:
                max_abs_F1 = abs(diff)
            # F3
            parity_pred = sum(x // d for d in divisors(rW)) % 2
            parity_emp = ldirect % 2
            # python mod is well-defined for negative, but ensure {0,1}:
            parity_emp = parity_emp & 1
            if parity_pred != parity_emp:
                max_abs_F3 += 1

        results["F1"][str(W)] = {
            "W": W,
            "rad_W": rW,
            "omega_W": ow_exact,
            "max_abs_diff": max_abs_F1,
            "n_divisors_used": calls_used[0],
            "predicted_2_omega": 2 ** ow_exact,
            "verdict": "PASS" if max_abs_F1 == 0 else "FAIL",
        }
        results["F3"][str(W)] = {
            "W": W,
            "n_parity_mismatches": max_abs_F3,
            "verdict": "PASS" if max_abs_F3 == 0 else "FAIL",
        }
        # F5: divisor count
        results["F5"][str(W)] = {
            "W": W,
            "n_divisors": calls_used[0],
            "predicted": 2 ** ow_exact,
            "verdict": "PASS" if calls_used[0] == 2 ** ow_exact else "FAIL",
        }
        print(
            f"  W={W:5d} rad={rW:5d} ω={ow_exact} 2^ω={2**ow_exact:3d} "
            f"max|diff|={max_abs_F1} parity_mismatches={max_abs_F3} "
            f"verdict={results['F1'][str(W)]['verdict']}"
        )

    # ---------------- F2 (radical reduction) ----------------
    print("\nF2 — radical reduction L_W = L_rad(W) for non-squarefree W:")
    for W in W_nonsqfree:
        rW = radical(W)
        max_abs = 0
        for x in x_grid_F1:
            l_W = L_W_direct(liouville, W, x)
            l_rW = L_W_direct(liouville, rW, x)
            d = l_W - l_rW
            if abs(d) > max_abs:
                max_abs = abs(d)
        results["F2"][str(W)] = {
            "W": W,
            "rad_W": rW,
            "max_abs_diff": max_abs,
            "verdict": "PASS" if max_abs == 0 else "FAIL",
        }
        print(
            f"  W={W:5d} rad={rW:5d} max|L_W − L_rad(W)|={max_abs} "
            f"verdict={results['F2'][str(W)]['verdict']}"
        )

    # ---------------- F4 (wheel-graded prime bisection) ----------------
    print("\nF4 — wheel-graded bisection π_W = (n_W − L_W)/2 − C_{3,W}:")
    f4_overall = "PASS"
    for W in W_all:
        rW = radical(W)
        max_diff = 0
        rows = []
        for x in x_grid_F4:
            n_W = n_W_direct(W, x)
            n_W_id = n_W_via_mobius(W, x, mu)
            assert n_W == n_W_id, (W, x, n_W, n_W_id)

            l_W = L_W_direct(liouville, W, x)
            a_W = (n_W - l_W) // 2
            a_W_dir = A_W_direct(bigomega, W, x)
            assert a_W == a_W_dir, (W, x, a_W, a_W_dir)

            pi_W = pi_W_direct(is_prime, W, x)
            c3_W = C3_W_direct(bigomega, is_prime, W, x)
            assert a_W - c3_W == pi_W, (W, x, a_W, c3_W, pi_W)

            # Direct identity (3d): π_W = (1/2) Σ_d (μ(d) ⌊x/d⌋ − L(⌊x/d⌋)) − C_{3,W}
            rhs = 0
            for d in divisors(rW):
                rhs += mu[d] * (x // d) - L_partial(L_arr, x // d)
            assert rhs % 2 == 0, (W, x, rhs)
            rhs = rhs // 2 - c3_W
            d = pi_W - rhs
            if abs(d) > max_diff:
                max_diff = abs(d)
            rows.append({"x": x, "pi_W": pi_W, "rhs": rhs, "diff": d})
        verdict = "PASS" if max_diff == 0 else "FAIL"
        if verdict == "FAIL":
            f4_overall = "FAIL"
        results["F4"][str(W)] = {
            "W": W,
            "max_abs_diff": max_diff,
            "rows": rows,
            "verdict": verdict,
        }
        print(
            f"  W={W:5d} rad={rW:5d} max|π_W − rhs|={max_diff} "
            f"verdict={verdict}"
        )

    # ---------------- F6 (mod-q lift; informational, not strict pass/fail) ----------------
    print("\nF6 — mod-q lift L_W mod q = (Σ L(⌊x/d⌋)) mod q:")
    for q in [2, 3, 4, 8]:
        for W in [2, 6, 30, 210, 2310]:
            rW = radical(W)
            mismatches = 0
            for x in x_grid_F4:
                l_W = L_W_direct(liouville, W, x) % q
                l_id = sum(L_partial(L_arr, x // d) for d in divisors(rW)) % q
                if l_W != l_id:
                    mismatches += 1
            key = f"q={q},W={W}"
            results["F6"][key] = {
                "q": q,
                "W": W,
                "mismatches_over_4_x": mismatches,
                "verdict": "PASS" if mismatches == 0 else "FAIL",
            }
            print(
                f"  q={q} W={W:5d} mismatches/4={mismatches} "
                f"verdict={results['F6'][key]['verdict']}"
            )

    # ---------------- Aggregate verdict ----------------
    pass_F1 = all(v["verdict"] == "PASS" for v in results["F1"].values())
    pass_F2 = all(v["verdict"] == "PASS" for v in results["F2"].values())
    pass_F3 = all(v["verdict"] == "PASS" for v in results["F3"].values())
    pass_F4 = f4_overall == "PASS"
    pass_F5 = all(v["verdict"] == "PASS" for v in results["F5"].values())
    pass_F6 = all(v["verdict"] == "PASS" for v in results["F6"].values())

    results["overall"] = {
        "F1": "PASS" if pass_F1 else "FAIL",
        "F2": "PASS" if pass_F2 else "FAIL",
        "F3": "PASS" if pass_F3 else "FAIL",
        "F4": "PASS" if pass_F4 else "FAIL",
        "F5": "PASS" if pass_F5 else "FAIL",
        "F6": "PASS" if pass_F6 else "FAIL",
    }
    print("\n=== AGGREGATE VERDICT ===")
    for k, v in results["overall"].items():
        print(f"  {k}: {v}")

    if out_path is not None:
        Path(out_path).write_text(json.dumps(results, indent=2))
        print(f"\nResults written to {out_path}")
    return results


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--N", type=int, default=10_000)
    p.add_argument(
        "--out",
        type=str,
        default=str(
            Path(__file__).parent / "wheel_coprime_liouville_results.json"
        ),
    )
    args = p.parse_args()
    main(N=args.N, out_path=args.out)
