"""
sanity_singular_series.py
=========================

Standalone numerical check that the truncated singular series

    S_Q(h) = sum_{q squarefree, 1<=q<=Q}  mu(q)^2 / phi(q)^2 * c_q(h)

converges (as Q -> infinity) to the textbook Hardy-Littlewood twin-prime
singular series

    S_HL(h) =
        2 * C_2 * prod_{p | h, p >= 3} (p-1)/(p-2)   if h even, h != 0,
        0                                            if h odd,

with C_2 = prod_{p >= 3} (1 - 1/(p-1)^2) ~ 0.66016 18158 ...

This validates the closed-form identity used in tq_correlation.py
WITHOUT using any T_Q evaluation, so the comparison is independent.
"""

from math import gcd, prod


def mobius_table(n_max: int) -> list[int]:
    """mu(n) for n in [0, n_max]; mu(0) = 0 by convention."""
    mu = [1] * (n_max + 1)
    mu[0] = 0
    is_prime = [True] * (n_max + 1)
    is_prime[0] = is_prime[1] = False
    for p in range(2, n_max + 1):
        if is_prime[p]:
            for m in range(p, n_max + 1, p):
                if m != p:
                    is_prime[m] = False
                mu[m] = -mu[m]
            p2 = p * p
            for m in range(p2, n_max + 1, p2):
                mu[m] = 0
    return mu


def phi_table(n_max: int) -> list[int]:
    """Euler totient phi(n) for n in [0, n_max]."""
    phi = list(range(n_max + 1))
    for p in range(2, n_max + 1):
        if phi[p] == p:  # p is prime
            for m in range(p, n_max + 1, p):
                phi[m] -= phi[m] // p
    return phi


def ramanujan_sum(q: int, h: int, phi_q: int, mu_table) -> int:
    """c_q(h) = mu(q/d) * phi(d) where d = gcd(q, h). Standard identity."""
    if q == 1:
        return 1
    d = gcd(q, h % q if q > 0 else 0)
    if d == 0:
        d = q  # gcd(q, 0) = q
    return mu_table[q // d] * (d if d == 1 else _phi_of(d))


# We keep a small phi cache to avoid recomputing phi(d) when d != q in
# c_q(h). Use the table directly via closure.

def make_c_q_evaluator(mu, phi):
    def c_q(q: int, h: int) -> int:
        if q == 1:
            return 1
        h_mod = h % q
        d = gcd(q, h_mod) if h_mod != 0 else q
        return mu[q // d] * phi[d]
    return c_q


def trunc_singular_series(Q: int, h: int, mu, phi, c_q) -> float:
    """Sum mu(q)^2 / phi(q)^2 * c_q(h) over squarefree q in [1, Q]."""
    s = 0.0
    for q in range(1, Q + 1):
        if mu[q] == 0:
            continue  # not squarefree
        s += (mu[q] * mu[q]) / (phi[q] * phi[q]) * c_q(q, h)
    return s


def textbook_S(h: int) -> float:
    """Full Hardy-Littlewood twin-prime singular series for shift h."""
    if h == 0:
        return float('inf')
    if h % 2 == 1:
        return 0.0
    # h even: 2 * C_2 * prod_{p | h, p odd} (p-1)/(p-2)
    # Compute C_2 by truncated Euler product over odd primes up to a cutoff.
    PMAX = 100_000
    is_prime = [True] * (PMAX + 1)
    is_prime[0] = is_prime[1] = False
    for p in range(2, int(PMAX ** 0.5) + 1):
        if is_prime[p]:
            for m in range(p * p, PMAX + 1, p):
                is_prime[m] = False
    primes = [p for p in range(2, PMAX + 1) if is_prime[p]]
    C2 = 1.0
    for p in primes:
        if p == 2:
            continue
        C2 *= 1 - 1.0 / (p - 1) ** 2
    # Odd prime factors of h
    h_abs = abs(h)
    n = h_abs
    while n % 2 == 0:
        n //= 2
    factors = []
    p = 3
    while p * p <= n:
        if n % p == 0:
            factors.append(p)
            while n % p == 0:
                n //= p
        p += 2
    if n > 1:
        factors.append(n)
    factor_prod = prod((p - 1) / (p - 2) for p in factors) if factors else 1.0
    return 2 * C2 * factor_prod


def main():
    Q_LIST = [10, 30, 100, 300, 1000, 3000, 10_000, 30_000]
    H_LIST = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 30, 210, 2310]

    Q_MAX = max(Q_LIST)
    print(f"# Building mu, phi tables up to Q_max = {Q_MAX}")
    mu = mobius_table(Q_MAX)
    phi = phi_table(Q_MAX)
    c_q = make_c_q_evaluator(mu, phi)

    print(
        "\n# S_Q(h) = sum_{q sqf <= Q} mu(q)^2 / phi(q)^2 * c_q(h)\n"
        "# Last column: textbook 2*C_2*prod (p-1)/(p-2) over odd p | h.\n"
    )
    header = "  h      " + "".join(f"   Q={Q:<6d}" for Q in Q_LIST) + "      S(h)_HL"
    print(header)
    for h in H_LIST:
        row = f"  h={h:<5d}"
        for Q in Q_LIST:
            S = trunc_singular_series(Q, h, mu, phi, c_q)
            row += f"  {S:>+9.5f}"
        S_HL = textbook_S(h)
        row += f"      {S_HL:>+9.5f}"
        print(row)

    # Spot check c_q identity for a few small (q, h)
    print("\n# Sanity: c_q(h) values for small q, h")
    print("    h\\q   1    2    3    4    5    6    7    8    9   10   12   15   30")
    Q_SAMPLE = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 15, 30]
    for h in [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12]:
        row = f"  h={h:<3d}"
        for q in Q_SAMPLE:
            row += f" {c_q(q, h):>+4d}"
        print(row)


if __name__ == "__main__":
    main()
