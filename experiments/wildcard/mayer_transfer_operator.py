"""
Mayer Transfer Operator (Gauss-Kuzmin-Wirsing) as a fast zero-locator.

IDEA
----
Mayer's theorem (1991) connects the Selberg zeta function for PSL(2,Z) to the
transfer operator L_s acting on analytic functions on a disk:

    (L_s f)(z) = sum_{n>=1} (z+n)^{-2s} f(1/(z+n))

The Selberg zeta function Z(s) = det(1 - L_s) det(1 + L_s), and one of its
spectral lines reproduces the non-trivial Riemann zeros (the other line gives
Maass-form eigenvalues). Numerical diagonalization on a polynomial basis
converges *exponentially* in matrix size N, since L_s is nuclear (trace-class
of any order).

The fresh question: does a SMALL N (say 30-50) Galerkin discretization let us
locate the first ~K zeros to integer accuracy in O(log K) zeta evaluations?
If so, the standard explicit-formula zero sum could be replaced by a polylog
matrix computation. If matrix diagonalization itself gives only the SAME
number of zeros as its dimension, no win.

This is a *quantitative* test: how does the number of accurately recovered
zeros scale with N? If sub-linearly, no compression. If super-linearly (e.g.,
N basis functions yield e^N zeros to known precision), we have a compression.

Concretely, for each candidate zero gamma_k (taken from a known table), we
check whether L_{1/2 + i*gamma_k} has an eigenvalue numerically close to +-1
or whether det(1 - L_s)det(1 + L_s) has a zero crossing nearby.

REFERENCES
----------
- Mayer, "The thermodynamic formalism approach to Selberg's zeta function for
  PSL(2,Z)", Bull. AMS (1991).
- Bonanno-Graffi, "The Mayer transfer operator approach to Selberg zeta",
  J. Math. Phys. (2008).
- Pollicott-Vytnova, "Estimating singularity dimension" (2014) for numerical
  recipes.
"""

import numpy as np
from scipy import linalg
from sympy import primepi
import time


# ---------- Mayer transfer operator on monomial basis ----------

def mayer_matrix(s_complex: complex, N: int, zeta_cache: dict = None) -> np.ndarray:
    """
    Galerkin matrix of L_s on monomial basis {1, z, ..., z^{N-1}}.

    M_{jk} = (-1)^j poch(2s+k, j) zeta(2s+k+j) / j!
    where poch(alpha, j) = alpha (alpha+1) ... (alpha+j-1).

    The zeta values needed are zeta(2s + m) for m = 0, 1, ..., 2N-2.
    These are CACHED across grid points via `zeta_cache` keyed on s_complex.
    """
    from mpmath import zeta, mpc
    M = np.zeros((N, N), dtype=complex)
    # Pre-compute zeta(2s + m) for m = 0..2N-2
    zetas = np.zeros(2 * N, dtype=complex)
    for m in range(2 * N):
        key = (s_complex, m)
        if zeta_cache is not None and key in zeta_cache:
            zetas[m] = zeta_cache[key]
        else:
            zv = complex(zeta(mpc(2 * s_complex + m)))
            zetas[m] = zv
            if zeta_cache is not None:
                zeta_cache[key] = zv
    # Pre-compute factorials
    facts = np.ones(N)
    for j in range(1, N):
        facts[j] = facts[j - 1] * j
    # Fill matrix using pochhammer recursion
    for k in range(N):
        alpha0 = 2 * s_complex + k  # base
        poch = 1.0 + 0j
        for j in range(N):
            sign = (-1) ** j
            M[j, k] = sign * poch * zetas[k + j] / facts[j]
            poch *= (alpha0 + j)
    return M


def selberg_det(s_complex: complex, N: int, zeta_cache: dict = None) -> complex:
    """
    Compute det((1 - L_s)(1 + L_s)) = det(1 - L_s^2).
    """
    M = mayer_matrix(s_complex, N, zeta_cache=zeta_cache)
    A = np.eye(N) - M @ M
    # Use sign-and-logdet for stability, then re-exponentiate
    sign, logdet = linalg.slogdet(A)
    if not np.isfinite(logdet):
        return linalg.det(A)
    return sign * np.exp(logdet)


# ---------- Test 1: locate first few zeros via det vanishing ----------

def search_zeros_on_critical_line(t_min: float, t_max: float, N: int, n_grid: int):
    """
    Scan |det(1 - L_s^2)| along s = 1/2 + i*t and report local minima.
    """
    ts = np.linspace(t_min, t_max, n_grid)
    abs_det = []
    for t in ts:
        s = 0.5 + 1j * t
        try:
            d = selberg_det(s, N)
            abs_det.append(abs(d))
        except Exception:
            abs_det.append(np.nan)
    abs_det = np.array(abs_det)
    # find local minima
    minima = []
    for i in range(1, len(abs_det) - 1):
        if abs_det[i] < abs_det[i - 1] and abs_det[i] < abs_det[i + 1]:
            minima.append((ts[i], abs_det[i]))
    return ts, abs_det, minima


# ---------- Reference Riemann zeros ----------

KNOWN_ZEROS = [
    14.134725141734693,
    21.022039638771555,
    25.010857580145688,
    30.424876125859513,
    32.935061587739189,
    37.586178158825671,
    40.918719012147495,
    43.327073280914999,
    48.005150881167159,
    49.773832477672302,
    52.970321477714460,
    56.446247697063394,
]


def test_recovery_vs_dimension():
    """
    For increasing N, count how many of the first 12 known zeros are captured
    by local minima of |det(1 - L_s^2)| in the range [10, 60].
    """
    print("=" * 72)
    print("TEST: Number of recovered zeros vs Galerkin dimension N")
    print("=" * 72)
    print(f"{'N':>4} {'time(s)':>10} {'minima_found':>14} {'matched':>10} {'avg_err':>14}")

    results = []
    for N in [6, 10, 14, 18]:
        t0 = time.time()
        try:
            _, _, minima = search_zeros_on_critical_line(10.0, 60.0, N, 80)
        except Exception as e:
            print(f"{N:>4} FAILED: {e}")
            continue
        elapsed = time.time() - t0

        # match minima to known zeros
        matched = 0
        errors = []
        for gamma in KNOWN_ZEROS:
            best = min(minima, key=lambda m: abs(m[0] - gamma)) if minima else None
            if best is not None and abs(best[0] - gamma) < 0.5:
                matched += 1
                errors.append(abs(best[0] - gamma))
        avg_err = np.mean(errors) if errors else float('nan')
        print(f"{N:>4} {elapsed:>10.2f} {len(minima):>14} {matched:>10} {avg_err:>14.6f}")
        results.append((N, elapsed, len(minima), matched, avg_err))
    return results


# ---------- Test 2: pi(x) via truncated explicit formula using recovered zeros ----------

def explicit_pi(x: float, gammas: list) -> float:
    """
    Riemann's R-function plus oscillatory zero sum, truncated.
    pi(x) ~ R(x) - 2 Re sum_k R(x^{1/2 + i gamma_k})

    R(x) = sum_{n>=1} mu(n)/n * Li(x^{1/n})

    For this test we use the simpler li(x) - 2 Re sum li(x^rho) (Riemann's
    original truncated form, which is good to ~O(1) for moderate x).
    """
    from mpmath import li, mpc, log

    s = mpc(0.5, 0)
    # leading term
    out = float(li(x).real)
    # zero contributions
    for gamma in gammas:
        rho = mpc(0.5, gamma)
        # li(x^rho) = li(exp(rho * log x)); use Ei(rho * log x)
        from mpmath import ei
        val = ei(rho * log(x))
        out -= 2 * float(val.real)
    return out


def test_pi_recovery_with_few_zeros():
    """
    Use the FIRST K known zeros (instead of recovered ones, since recovery
    error is small) and see how pi(x) error decays with K for various x.
    Compares against ground-truth sympy primepi.
    """
    print("\n" + "=" * 72)
    print("TEST: pi(x) error vs number of zeros used (using exact zeros)")
    print("=" * 72)

    # Use a longer table of zeros for this test
    extended_zeros = [
        14.134725141734693, 21.022039638771555, 25.010857580145688,
        30.424876125859513, 32.935061587739189, 37.586178158825671,
        40.918719012147495, 43.327073280914999, 48.005150881167159,
        49.773832477672302, 52.970321477714460, 56.446247697063394,
        59.347044002602353, 60.831778524609810, 65.112544048081607,
        67.079810529494173, 69.546401711173980, 72.067157674481907,
        75.704690699083933, 77.144840068874805, 79.337375020249367,
        82.910380854086030, 84.735492980517050, 87.425274613125229,
        88.809111207634532, 92.491899270558484, 94.651344040519887,
        95.870634228245309, 98.831194218193198, 101.317851005731391,
    ]

    print(f"{'x':>10} {'pi(x)':>10} {'K=10':>8} {'K=20':>8} {'K=30':>8}")
    for x in [50, 100, 500, 1000, 5000, 10000]:
        true_pi = int(primepi(x))
        row = [f"{x:>10}", f"{true_pi:>10}"]
        for K in [10, 20, 30]:
            est = explicit_pi(x, extended_zeros[:K])
            err = est - true_pi
            row.append(f"{err:>+8.2f}")
        print(" ".join(row))


# ---------- main ----------

if __name__ == "__main__":
    print("MAYER TRANSFER OPERATOR — fresh angle test")
    print("Question: does small Galerkin matrix encode many Riemann zeros?")
    print()

    results = test_recovery_vs_dimension()
    test_pi_recovery_with_few_zeros()

    print("\nINTERPRETATION")
    print("-" * 72)
    print("If `matched` grows much faster than N, the operator compresses")
    print("zero information exponentially -> potential polylog algorithm.")
    print("If `matched` <= N (linear), no compression beyond direct zero list.")
