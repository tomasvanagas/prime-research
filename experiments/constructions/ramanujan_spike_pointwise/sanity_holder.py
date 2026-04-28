"""
Sanity check: Hoelder identity for Ramanujan sums on squarefree q.

For squarefree q with gcd(q,n) = d:
    c_q(n)            = mu(q/d) * phi(d)                  (standard Hoelder)
    mu(q) c_q(n)/phi(q) = mu(d) / phi(q/d)                (project simplification)

This script computes both sides directly for q in [1, 30], n in [0, 60]
and verifies they agree.
"""

import math
import cmath


def gcd(a, b):
    while b:
        a, b = b, a % b
    return a


def mobius(n):
    if n == 1:
        return 1
    x, cnt = n, 0
    p = 2
    while p * p <= x:
        if x % p == 0:
            x //= p
            cnt += 1
            if x % p == 0:
                return 0
        p += 1
    if x > 1:
        cnt += 1
    return -1 if cnt % 2 else 1


def totient(n):
    res = n
    x = n
    p = 2
    while p * p <= x:
        if x % p == 0:
            res -= res // p
            while x % p == 0:
                x //= p
        p += 1
    if x > 1:
        res -= res // x
    return res


def is_squarefree(n):
    return mobius(n) != 0


def ramanujan_sum_direct(q, n):
    s = 0.0 + 0.0j
    for j in range(1, q + 1):
        if gcd(j, q) == 1:
            s += cmath.exp(2j * math.pi * j * n / q)
    return s.real  # purely real for c_q(n)


def main():
    fails = []
    for q in range(1, 31):
        if not is_squarefree(q):
            continue
        muq, phiq = mobius(q), totient(q)
        for n in range(0, 61):
            d = gcd(q, n)
            phi_qd = totient(q // d)
            mud = mobius(d)
            # LHS: project simplification
            lhs = mud / phi_qd
            # RHS: directly compute mu(q) c_q(n)/phi(q)
            cqn = ramanujan_sum_direct(q, n)
            rhs = muq * cqn / phiq
            if abs(lhs - rhs) > 1e-9:
                fails.append((q, n, lhs, rhs, cqn))
    if not fails:
        print("OK: Hoelder simplification mu(q) c_q(n)/phi(q) = mu(gcd) / phi(q/gcd)")
        print("    verified for all squarefree q in [1, 30] and n in [0, 60].")
    else:
        print(f"FAIL: {len(fails)} mismatches.  Sample:")
        for q, n, lhs, rhs, c in fails[:10]:
            print(f"  q={q} n={n}  lhs={lhs:.6f}  rhs={rhs:.6f}  c_q(n)={c:.6f}")


if __name__ == "__main__":
    main()
