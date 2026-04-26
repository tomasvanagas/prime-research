"""
Proposal 25: Liouville-Parity Triangulation

Identity (corrected during this experiment): with
    L(x) = sum_{n<=x} lambda(n)              (Liouville summatory)
    C_3(x) = #{n <= x : Omega(n) is odd and >= 3}
the count #{n <= x : Omega(n) odd} = (x - L(x))/2 since Omega(1)=0 is even,
so
    pi(x) = (x - L(x))/2 - C_3(x)
and in particular
    pi(x) mod 2 == ((x - L(x))/2 - C_3(x)) mod 2.

The original proposal mis-identified C_3 with (Q(x)-1).  This experiment:
  (a) confirms the corrected identity,
  (b) probes whether C_3(x) admits a polylog parity oracle of its own,
  (c) tests truncated-L_K reconstructions.

Run:  python3 proposal25_liouville_parity_check.py
"""

import time
from sympy import factorint, primepi


def Omega(n: int) -> int:
    """Total number of prime factors of n with multiplicity."""
    if n == 1:
        return 0
    return sum(factorint(n).values())


def liouville(n: int) -> int:
    return -1 if Omega(n) % 2 else 1


def squarefree(n: int) -> bool:
    if n == 1:
        return True
    for p, e in factorint(n).items():
        if e >= 2:
            return False
    return True


def main(N: int = 10000) -> None:
    print(f"Computing arithmetic functions for x in [1, {N}]...")
    t0 = time.time()
    omega_arr = [0] * (N + 1)         # Omega(n) with multiplicity
    for n in range(1, N + 1):
        omega_arr[n] = Omega(n)
    lam = [(-1) ** omega_arr[n] for n in range(N + 1)]
    print(f"  done in {time.time() - t0:.1f}s")

    # Cumulative sums
    L = [0] * (N + 1)
    C3 = [0] * (N + 1)         # #{n <= x : Omega(n) odd and >= 3}
    pi_arr = [0] * (N + 1)
    for n in range(1, N + 1):
        L[n] = L[n - 1] + lam[n]
        is_C3 = 1 if (omega_arr[n] >= 3 and omega_arr[n] % 2 == 1) else 0
        C3[n] = C3[n - 1] + is_C3
        pi_arr[n] = pi_arr[n - 1] + (1 if omega_arr[n] == 1 else 0)

    # Verify CORRECTED identity:  pi(x) = (x - L(x))/2 - C_3(x)
    print("\nVerifying CORRECTED identity  pi(x) = (x - L(x))/2 - C_3(x)  ...")
    integer_match = 0
    parity_match = 0
    for x in range(2, N + 1):
        diff = x - L[x]
        if diff % 2 == 0:
            rhs_int = diff // 2 - C3[x]
            if rhs_int == pi_arr[x]:
                integer_match += 1
            if (rhs_int % 2) == (pi_arr[x] % 2):
                parity_match += 1
    total = N - 1
    print(f"  integer match: {integer_match}/{total}")
    print(f"  parity  match: {parity_match}/{total}")

    # Now the algorithmic question: can we compute C_3(x) more cheaply than pi(x)?
    # By definition C_3(x) >= pi(x^{1/3}) (every cube of a prime contributes if
    # Omega = 3 odd, e.g. p^3 has Omega=3); and C_3(x) <= sum of multinomials.
    # But more importantly C_3(x) = #{Omega=3,5,7,...} <= x.
    # Compare growth:
    print("\nGrowth comparison pi(x) vs C_3(x):")
    for x in [10, 100, 1000, 10000]:
        print(f"  x={x:>5}:  pi={pi_arr[x]:>5}, "
              f"C_3={C3[x]:>5}, "
              f"L={L[x]:>5}, "
              f"(x-L)/2={(x - L[x]) // 2:>5}, "
              f"sum check: {(x - L[x]) // 2 - C3[x]} vs pi={pi_arr[x]}")

    # The KEY question for Proposal 25:
    # If we want pi(x) parity, we need both L(x) parity and C_3(x) parity.
    # C_3(x) counts integers with >= 3 prime factors of *odd* total multiplicity.
    # For x ~ 10^4, C_3(x) is ~30% of x; computing it from scratch is no easier
    # than pi(x).  But maybe its parity has structure?
    #
    # Probe: does C_3(x) parity admit a divisor-sum identity?
    # By inclusion-exclusion on Omega = 3,5,7,...:
    #   C_3(x) = sum_{k odd, k>=3} pi_k(x)   where pi_k(x) = #{Omega(n)=k, n<=x}
    # Each pi_k can be expressed via Dirichlet convolutions.  Skip for now;
    # just measure: how often does C_3 mod 2 equal a "fast" approximation?
    print("\nProbing C_3(x) parity vs simple proxies...")
    proxies = {
        "pi(x^(1/3))": [int(primepi(int(x ** (1.0 / 3.0)))) for x in range(N + 1)],
        "pi(x^(1/2))": [int(primepi(int(x ** 0.5))) for x in range(N + 1)],
        "x // 8":    [x // 8 for x in range(N + 1)],
        "L(x)":      L,
    }
    for name, arr in proxies.items():
        agree = sum(1 for x in range(2, N + 1) if (arr[x] % 2) == (C3[x] % 2))
        print(f"  C_3 parity agrees with parity({name}) for {agree}/{total} = "
              f"{agree/total:.3f}")

    # Truncated-L probe.  Define
    #     L_K(x) = sum_{n<=x, Omega(n) <= K} lambda(n)
    print("\nTruncated L_K probe (parity match with full L):")
    for K in [1, 2, 3, 4, 5, 6, 8, 10]:
        L_K = [0] * (N + 1)
        for n in range(1, N + 1):
            if omega_arr[n] <= K:
                L_K[n] = L_K[n - 1] + lam[n]
            else:
                L_K[n] = L_K[n - 1]
        agree_full = sum(1 for x in range(1, N + 1) if (L_K[x] % 2) == (L[x] % 2))
        print(f"  K={K:2d}:  L_K parity agrees with L for {agree_full}/{N} x "
              f"({agree_full/N:.3f})")


if __name__ == "__main__":
    main()
