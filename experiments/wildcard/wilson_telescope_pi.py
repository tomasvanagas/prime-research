"""
Wilson's theorem with telescoping factorial residues.

Wilson: (p-1)! ≡ -1 mod p iff p prime.
=> floor( ((n-1)! mod n) / (n-1) ) = 1 if n prime, 0 otherwise (for n >= 2).
=> pi(x) = sum_{n=2..x} floor((n-1)!  mod n / (n-1)).

Computing (n-1)! mod n for each n separately is O(n) per n, O(x^2) total.

KEY IDEA: Can we compute the SEQUENCE { (n-1)! mod n }_{n=2..x} faster
collectively using telescoping or product-tree structure?

Observation: (n-1)! mod n is NOT related to (n-2)! mod (n-1) in any simple way
(the modulus changes). So naive telescoping fails.

BUT: We can fix a single large modulus M (much greater than x) and work with
factorials mod M. By CRT, knowing M! mod M determines (n-1)! mod n for all
n with n | M, but n ranges over arbitrary integers, so M would need to be
LCM(1..x), which is exp(x).

Alternative: COMPUTE (n-1)! mod n in batches via Wilson PSEUDOPRIME pattern.

This experiment:
1. Verify Wilson formula on small range.
2. Time the naive computation vs batched approaches.
3. Try a "polynomial product tree" approach: build (x-1)(x-2)...(x-(n-1))
   as polynomial mod n, then evaluate at x=0 -- gives (-1)^{n-1} (n-1)!.
   This is still O(n log^2 n) per n.
4. Estimate the COMPLEXITY BARRIER explicitly: any algorithm computing
   (n-1)! mod n must do O(log n) multiplications at minimum... but
   does it need log(n!) ~ n log n bits of work overall?
"""

import time
from sympy import isprime
from math import log, log2


def wilson_count(x):
    """pi(x) via Wilson's formula, naive."""
    count = 0
    fact_mod_n_record = []
    for n in range(2, x + 1):
        f = 1
        for k in range(1, n):
            f = (f * k) % n
        # f is (n-1)! mod n
        # if f == n-1 (== -1 mod n), n is prime
        if f == n - 1:
            count += 1
        fact_mod_n_record.append((n, f))
    return count, fact_mod_n_record


def fast_factorial_mod(n):
    """Compute (n-1)! mod n via product tree in O((log n)^2 M(log n)) bits."""
    if n < 2:
        return 0
    # Build product 1*2*...*(n-1) mod n via balanced product tree
    vals = list(range(1, n))
    while len(vals) > 1:
        new = []
        for i in range(0, len(vals), 2):
            if i + 1 < len(vals):
                new.append((vals[i] * vals[i + 1]) % n)
            else:
                new.append(vals[i])
        vals = new
    return vals[0]


def time_naive_vs_tree(N):
    t0 = time.time()
    count1 = 0
    for n in range(2, N + 1):
        f = 1
        for k in range(1, n):
            f = (f * k) % n
        if f == n - 1:
            count1 += 1
    t1 = time.time()

    t2 = time.time()
    count2 = 0
    for n in range(2, N + 1):
        f = fast_factorial_mod(n)
        if f == n - 1:
            count2 += 1
    t3 = time.time()

    return (t1 - t0, count1), (t3 - t2, count2)


def estimate_bit_complexity_lower_bound(N):
    """
    Each (n-1)! mod n requires reading the modulus n (log n bits) and producing
    log n bits of output. The minimum number of bit operations across all n in [2..N]
    is at least N log N (output size).

    For a polylog algorithm computing pi(N), we'd need to AVOID computing each
    (n-1)! mod n separately. We need a single quantity that captures pi(N) in
    polylog space.

    Question: is there a single integer Q computable in polylog(N) bits such that
    pi(N) = f(Q) with f also polylog?

    The Mertens-like sum sum_{n} mu(n)/n * f(N/n) compresses some info, but
    it requires f(N/n) for many n.
    """
    # Collect Wilson-failure differences (n-1)! mod n - (n-1)
    diffs = []
    for n in range(2, min(N, 200) + 1):
        f = fast_factorial_mod(n)
        diffs.append(((n, f, n - 1, (f - (n - 1)) % n, isprime(n))))
    return diffs


def main():
    print("=" * 60)
    print("Wilson telescoping investigation")
    print("=" * 60)

    print("\nVerify formula on N=200:")
    count, _ = wilson_count(200)
    true_pi_200 = sum(1 for n in range(2, 201) if isprime(n))
    print(f"  Wilson count: {count}, true pi(200): {true_pi_200}, match: {count == true_pi_200}")

    print("\nTiming naive vs balanced product tree:")
    for N in [200, 1000, 3000]:
        (t_naive, c_n), (t_tree, c_t) = time_naive_vs_tree(N)
        print(f"  N={N}: naive={t_naive:.3f}s, tree={t_tree:.3f}s (counts={c_n},{c_t})")

    print("\nNote: product tree per n is O((log n)^2 M(log n)).")
    print(f"Total over [2..N] -> O(N (log N)^2 M(log N)) = O(N polylog N).")
    print(f"This is *worse* than the naive sieve. Wilson is computationally useless")
    print(f"as direct prime counting unless we can SHARE work across different n.")

    print("\nFundamental obstruction:")
    print("- Each modulus n is distinct; no global ring contains all (n-1)! mod n.")
    print("- LCM(1..N) ~ exp(N) by PNT, so a single big modulus is exp-size.")
    print("- A polylog Wilson-based algorithm is information-theoretically blocked")
    print("  unless we can share factorials across moduli.")


if __name__ == "__main__":
    main()
