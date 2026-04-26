"""
P2: Mollifier-corrected explicit formula.

Selberg's mollifier replaces zeta(s) by zeta(s) * M(s), where
    M(s) = sum_{n <= Y} a_n n^{-s}
is a Dirichlet polynomial whose values at the first K critical zeros
{0.5 + i gamma_j : j = 1..K} are forced to be zero (or as small as
possible). The product zeta * M then has those K zeros suppressed in
the residue sum.

Concretely we use the explicit-formula approximation:
    pi(x) ~ R(x) - sum_{|gamma| <= T} Re[ Ei(rho * log x) ]
where rho = 1/2 + i gamma. Multiplying inside the sum by the mollifier
weights M(rho)/M(1) produces a CORRECTED estimate; if M(rho_j) ~ 0 for
j <= K, those terms vanish even though we never actually summed them.

The QUESTION is: at what cost in T can we compute pi(x) to error < 0.5
when the mollifier zeros K zeros?

Plan
----
1. Read first 1000 zeros.
2. Build a length-Y mollifier minimising sum_{j<=K} |M(rho_j)|^2 by
   solving a least-squares problem in the Y mollifier coefficients.
3. For each (K, Y) pair and each x in {100, 1000, 10000}, compare:
      sharp_T(x) = R(x) - sum_{|gamma_j| <= T} Ei(rho_j log x)
   to
      mollified_T(x) = R(x) - sum_{|gamma_j| <= T}
                       Ei(rho_j log x) * (M(rho_j) / M(1))_real
   plus correction terms for 1 <= j <= K which were zeroed.

We measure the smallest T such that the mollified estimate stays within
+/- 0.5 of pi(x).

Falsification: increasing K does NOT reduce the required T.

Important caveat
----------------
Multiplying inside the sum is NOT mathematically equivalent to running
the explicit formula for zeta*M. The explicit formula for zeta*M would
involve the von-Mangoldt-like coefficients of zeta'/zeta + M'/M, which
introduce a primary contribution from M itself. But for the *zero
contribution* alone, the multiplication-by-M-evaluated-at-rho weighting
captures the dominant correction. We are testing whether this
heuristic helps.

The clean version of the experiment: solve for a_n that minimise
|M(rho_j)|^2 + lambda * sum |a_n|^2 (Tikhonov), then look at the
smallest T at which the weighted partial sum is within ±0.5 of pi(x).
"""

import sys
import numpy as np
from mpmath import mp, mpc, mpf, log as mlog, ei, riemannr, log

mp.prec = 100

ZEROS_PATH = "data/zeta_zeros_1000.txt"

def load_zeros():
    with open(ZEROS_PATH) as f:
        return [mpf(line.strip()) for line in f if line.strip()]


def li_complex(rho, log_x):
    """Compute Ei(rho * log x). Use mpmath for accuracy."""
    z = rho * log_x
    return ei(z)


def sharp_partial_sum(zeros, x, T):
    """R(x) - sum_{|gamma_j|<=T} Re[ Ei(rho_j log x) ] (using rho and conjugate)."""
    log_x = log(mpf(x))
    R = riemannr(mpf(x))
    total = mpc(0)
    n_used = 0
    for g in zeros:
        if g > T:
            break
        rho = mpc(mpf("0.5"), g)
        rho_bar = mpc(mpf("0.5"), -g)
        total += li_complex(rho, log_x) + li_complex(rho_bar, log_x)
        n_used += 1
    return float(R - total.real), n_used


def mollifier_evaluate(coeffs_an, n_values, s):
    """Evaluate M(s) = sum_n a_n n^{-s} where coeffs_an[i] is a_{n_values[i]}."""
    return sum(a * mpc(n) ** (-s) for a, n in zip(coeffs_an, n_values))


def build_mollifier(zeros, K, Y, lam=1e-6):
    """Solve for length-Y mollifier coefficients a_1..a_Y minimising
        sum_{j=1..K} |M(0.5+i gamma_j)|^2 + lam * ||a||^2
       subject to a_1 = 1 (normalisation).

       This is a complex least-squares: build the design matrix
       D_{j, n} = n^{-(0.5 + i gamma_j)} for j=1..K, n=1..Y, and
       minimise ||D a||^2 + lam ||a||^2 with a[0] = 1.
    """
    n_values = list(range(1, Y + 1))
    # use double precision floats here for speed; zeros are ~1e-30 accurate
    Df = np.zeros((K, Y), dtype=complex)
    for j in range(K):
        gamma = float(zeros[j])
        for i, n in enumerate(n_values):
            log_n = np.log(n)
            # n^{-(0.5 + i gamma)} = n^{-0.5} * exp(-i gamma log n)
            Df[j, i] = (n ** -0.5) * np.exp(-1j * gamma * log_n)
    # Real-stack the system: split D into [Re(D); Im(D)] of shape (2K, Y)
    A_real = np.vstack([Df.real, Df.imag])  # shape (2K, Y)
    # Constrain a_1 = 1: substitute a_1=1, solve for a_2..a_Y
    # equations: A_real[:, 0] * 1 + A_real[:, 1:] * a_rest = -A_real[:, 0]
    #          → A_real[:, 1:] * a_rest = -A_real[:, 0]
    rhs = -A_real[:, 0]
    A_sub = A_real[:, 1:]
    # add Tikhonov regularisation: minimise ||A_sub @ a_rest - rhs||^2 + lam ||a_rest||^2
    n_unknowns = Y - 1
    A_aug = np.vstack([A_sub, np.sqrt(lam) * np.eye(n_unknowns)])
    rhs_aug = np.concatenate([rhs, np.zeros(n_unknowns)])
    a_rest, *_ = np.linalg.lstsq(A_aug, rhs_aug, rcond=None)
    coeffs = np.concatenate([[1.0], a_rest])
    return coeffs, n_values


def mollified_partial_sum(zeros, x, T, coeffs, n_values):
    """R(x) - sum_{|gamma_j| <= T} (M(rho_j)/M(1)) * Ei(rho_j log x) + conj.
    """
    log_x = log(mpf(x))
    R = riemannr(mpf(x))
    M1 = sum(coeffs)  # M(1) = sum a_n / n^1 ... actually M(1) = sum a_n /n
    M1 = sum(a / n for a, n in zip(coeffs, n_values))
    total = mpc(0)
    n_used = 0
    for g in zeros:
        if g > T:
            break
        rho = mpc(mpf("0.5"), g)
        rho_bar = mpc(mpf("0.5"), -g)
        # M(rho) for double-precision rho (we already computed coeffs in float)
        gamma_f = float(g)
        Mrho = sum(
            a * (n ** -0.5) * np.exp(-1j * gamma_f * np.log(n))
            for a, n in zip(coeffs, n_values)
        )
        weight = Mrho / float(M1)  # complex weight
        # weight * Ei(rho log x) + conj weight * Ei(rho_bar log x)
        # = 2 * Re(weight * Ei(rho log x))
        ei_val = li_complex(rho, log_x)
        # convert weight to mpc
        w = mpc(weight.real, weight.imag)
        total += w * ei_val + mpc(weight.real, -weight.imag) * li_complex(rho_bar, log_x)
        n_used += 1
    return float(R - total.real), n_used


def pi_exact(x):
    """Exact pi(x) via sieve."""
    if x < 2:
        return 0
    sieve = np.ones(int(x) + 1, dtype=bool)
    sieve[:2] = False
    for i in range(2, int(x ** 0.5) + 1):
        if sieve[i]:
            sieve[i * i::i] = False
    return int(sieve.sum())


def main():
    print("=" * 70)
    print("P2: Mollifier-corrected explicit formula for pi(x)")
    print("=" * 70)

    zeros = load_zeros()
    print(f"\nLoaded {len(zeros)} zeta zeros (first at gamma_1 = {float(zeros[0]):.6f}).")

    test_xs = [100, 1000, 10000]
    pi_truth = {x: pi_exact(x) for x in test_xs}
    print(f"\nGround truth pi(x):")
    for x in test_xs:
        print(f"  pi({x:>5}) = {pi_truth[x]}")

    # --- baseline: sharp partial sum at various T ---
    print("\n--- Baseline: sharp explicit-formula partial sums (no mollifier) ---")
    print(f"{'T':>6} | " + " | ".join(f"x={x} (true={pi_truth[x]})".rjust(20) for x in test_xs))
    for T in [50, 100, 200, 500, 1000]:
        row = [f"{T:>6}"]
        for x in test_xs:
            est, n = sharp_partial_sum(zeros, x, T)
            err = est - pi_truth[x]
            row.append(f"{est:>9.3f} ({err:+.2f})".rjust(20))
        print(" | ".join(row))

    # --- mollified ---
    print("\n--- With mollifier (zeroing first K zeros, length Y) ---")
    configs = [(0, 1), (5, 30), (10, 50), (20, 100), (40, 200)]
    for K, Y in configs:
        if K == 0:
            print(f"\n[K={K}, Y={Y}]  (no mollifier baseline)")
            continue
        print(f"\n[K={K}, Y={Y}]  mollifier suppressing first {K} zeros")
        coeffs, n_vals = build_mollifier(zeros, K, Y, lam=1e-8)
        # measure how well M zeroes the first K zeros
        max_M_first_K = 0.0
        for j in range(K):
            gamma = float(zeros[j])
            Mr = sum(a * (n ** -0.5) * np.exp(-1j * gamma * np.log(n)) for a, n in zip(coeffs, n_vals))
            max_M_first_K = max(max_M_first_K, abs(Mr))
        print(f"   max |M(rho_j)| for j=1..{K}: {max_M_first_K:.3e}")
        # M(1)
        M1_val = sum(a / n for a, n in zip(coeffs, n_vals))
        print(f"   M(1) = {M1_val:.4f}")
        print(f"   {'T':>6} | " + " | ".join(f"x={x}".rjust(22) for x in test_xs))
        for T in [50, 100, 200, 500, 1000]:
            row = [f"{T:>6}"]
            for x in test_xs:
                est, n = mollified_partial_sum(zeros, x, T, coeffs, n_vals)
                err = est - pi_truth[x]
                row.append(f"{est:>9.3f} ({err:+.2f})".rjust(22))
            print(" | ".join(row))

    print("\n" + "=" * 70)
    print("Interpretation: if mollifier truly suppresses early zeros' contribution,")
    print("the sharp partial sum (which truncates them) should AGREE with the")
    print("mollified partial sum at much SMALLER T.")
    print("Equivalently: mollified-T-error should track sharp-(T+contribution)-error.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
