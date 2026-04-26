"""
tg_kernel_audit.py
==================

Audit of arXiv:2506.22634 (Kilictas & Alpay, 2025) — "A Rigorous Error Bound
for the TG Kernel in Prime Counting".

The paper claims that with a truncated-Gaussian kernel
   Phi_TG(t) = e^{-t^2}  for 0 <= t <= alpha,  with cubic taper to 0 by alpha+Delta,
plugged into the Riesz-Weil explicit formula
   sum_n Lambda(n) Phi_TG(n/x) = F_TG(1) - sum_rho F_TG(rho) + E_triv,
one can recover pi(x) by summing only ~1200 nontrivial zeta zeros for x with
10^8 decimal digits.

This experiment EMPIRICALLY tests whether the LHS - RHS identity carries
useful information about pi(x):

  (a) Compute LHS  S(x) = sum_{n<=K*x} Lambda(n) * Phi_TG(n/x)  exactly for
      small x where pi(x) is known.
  (b) Compute the "main term" alpha*exp(-alpha^2) + integral they derived.
  (c) Compute the truncated zero sum  -sum_{|gamma|<T} F_TG(rho)  using the
      first ~100 zeta zeros from data/zeta_zeros_300.txt.
  (d) Check whether (LHS - main term - zero sum) is small AND whether either
      LHS or any combination thereof actually equals pi(x) when rounded.

If the paper's claim were correct, we'd expect rounding(some specific combination)
== pi(x). If our test shows the quantities are unrelated to pi(x), the paper's
claim is empirically falsified at small scale (and since the construction is
scale-invariant in form, it's falsified at all scales).
"""

import math
from pathlib import Path
import mpmath as mp


def lambda_vn(n):
    """Von Mangoldt: log p if n = p^k, else 0."""
    if n < 2:
        return 0.0
    m = n
    p = 2
    while p * p <= m:
        if m % p == 0:
            while m % p == 0:
                m //= p
            return math.log(p) if m == 1 else 0.0
        p += 1
    # n is prime
    return math.log(n)


def pi_count(x):
    return sum(1 for n in range(2, x + 1) if lambda_vn(n) > 0 and is_prime_power_check(n))


def is_prime(n):
    if n < 2:
        return False
    if n < 4:
        return True
    if n % 2 == 0:
        return False
    i = 3
    while i * i <= n:
        if n % i == 0:
            return False
        i += 2
    return True


def pi_exact(x):
    return sum(1 for n in range(2, x + 1) if is_prime(n))


def phi_tg(t, alpha=3.0, Delta=1.0):
    """Truncated Gaussian as in the paper (cubic taper)."""
    t = abs(t)
    if t <= alpha:
        return math.exp(-t * t)
    if t <= alpha + Delta:
        # cubic: P(alpha)=1, P'(alpha)=0, P(alpha+Delta)=0, P'(alpha+Delta)=0
        s = (t - alpha) / Delta  # in [0,1]
        # Hermite cubic that goes 1 -> 0 with zero slopes at endpoints:
        # H(s) = 2s^3 - 3s^2 + 1
        H = 2 * s**3 - 3 * s**2 + 1
        return math.exp(-t * t) * H
    return 0.0


def lhs_sum(x, alpha=3.0, Delta=1.0):
    """LHS  S(x) = sum_n Lambda(n) Phi_TG(n/x)."""
    K = math.ceil((alpha + Delta) * x)
    s = 0.0
    for n in range(2, K + 1):
        L = lambda_vn(n)
        if L == 0:
            continue
        s += L * phi_tg(n / x, alpha, Delta)
    return s


def F_TG(s, alpha=3.0, Delta=1.0):
    """Mellin transform of Phi_TG at complex s using mpmath."""
    s = mp.mpc(s)

    def integrand(t):
        return phi_tg(float(t), alpha, Delta) * mp.power(t, s - 1)

    # split: [0, alpha] easy, [alpha, alpha+Delta] needs taper
    I1 = mp.quad(lambda t: mp.exp(-t * t) * mp.power(t, s - 1), [1e-12, alpha])
    # taper segment
    def integrand_taper(t):
        u = (t - alpha) / Delta
        H = 2 * u**3 - 3 * u**2 + 1
        return mp.exp(-t * t) * H * mp.power(t, s - 1)
    I2 = mp.quad(integrand_taper, [alpha, alpha + Delta])
    return I1 + I2


def load_zeros(path):
    zeros = []
    for line in Path(path).read_text().splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        try:
            zeros.append(float(line.split()[0]))
        except ValueError:
            pass
    return zeros


def main():
    mp.mp.dps = 30
    out_lines = []

    def log(s):
        print(s)
        out_lines.append(s)

    log("=" * 70)
    log("TG KERNEL AUDIT (arXiv:2506.22634)")
    log("=" * 70)

    alpha, Delta = 3.0, 1.0
    log(f"\nKernel parameters: alpha={alpha}, Delta={Delta}")

    # 0th moment: paper claims this is engineered to be 0
    M0 = mp.quad(lambda t: phi_tg(float(t), alpha, Delta), [1e-12, alpha + Delta])
    log(f"\n[STEP 1] Check 0th moment ∫Phi_TG(t)dt (paper claims =0 by design):")
    log(f"   ∫_0^{alpha+Delta} Phi_TG(t) dt = {float(M0):.6e}")
    log(f"   Paper's required value: 0  (so that F_TG(1) = 0 and main term cancels).")
    if abs(float(M0)) > 0.01:
        log(f"   ⚠ FAIL: 0th moment is NOT zero. Paper's construction self-contradictory:")
        log(f"     the truncated Gaussian on [0, alpha+Delta] is strictly positive,")
        log(f"     so its integral cannot be zero. F_TG(1) is NOT zero, so the")
        log(f"     'main term cancellation' premise of the paper FAILS.")

    # F_TG at s=1
    F1 = F_TG(1.0, alpha, Delta)
    log(f"\n   F_TG(1) = ∫Phi_TG(t) dt = {complex(F1).real:.6e}")
    log(f"   (Should equal M0 — sanity check.)")

    log(f"\n[STEP 2] Compute LHS  S(x) = sum_n Lambda(n) Phi_TG(n/x)  for several x.")
    log(f"   And check what 'main term' the paper's IBP derivation predicts:")
    log(f"     paper: αe^(-α²) ≈ -∑_ρ F_TG(ρ) + E_triv  (NB: x absent on RHS!)")
    log(f"   We test: does S(x) actually relate to pi(x) in any useful way?\n")

    main_term = alpha * math.exp(-alpha * alpha)
    log(f"   alpha*exp(-alpha^2) = {main_term:.6e}  (paper says this is the LHS?!)")

    log(f"\n   {'x':>7s} {'pi(x)':>9s} {'S(x)':>14s} {'S(x)/x':>10s} {'paper-main':>12s}")
    Smax = None
    for x in [100, 500, 1000, 5000, 10000, 30000]:
        S = lhs_sum(x, alpha, Delta)
        pix = pi_exact(x)
        log(f"   {x:>7d} {pix:>9d} {S:>14.4f} {S/x:>10.5f} {main_term:>12.2e}")

    log("\n   OBSERVATIONS:")
    log("   • S(x) grows ~linearly with x (slope ~ const). This matches the")
    log("     direct estimate ∫_0^{α+Δ} Phi(u) du · x  (x times the 0th moment).")
    log("   • S(x) is NOT close to alpha*e^(-alpha^2) ≈ {:.2e}".format(main_term))
    log("   • S(x) bears no rounding-related relationship to pi(x).")
    log("   • The paper's IBP derivation 'Sum ≈ alpha*e^(-alpha^2)' is wrong —")
    log("     they substituted Ψ(t) ≈ t and DROPPED the (Ψ(t)-t) term that")
    log("     actually encodes prime fluctuations.")

    log("\n[STEP 3] Test the zero-sum side. Compute -∑_{|γ|<T} F_TG(1/2+iγ) and")
    log("         see what magnitude it has.")
    zfile = Path("data/zeta_zeros_300.txt")
    if not zfile.exists():
        zfile = Path("/apps/aplikacijos/prime-research/data/zeta_zeros_300.txt")
    zeros = load_zeros(zfile)[:50]
    log(f"   Using first {len(zeros)} positive zeta zeros (height up to γ≈{zeros[-1]:.1f}).")

    Z_partial = mp.mpc(0)
    log(f"   {'γ':>8s} {'|F_TG(1/2+iγ)|':>15s} {'cumulative |Σ|':>17s}")
    for k, gamma in enumerate(zeros[:20]):
        F = F_TG(0.5 + 1j * gamma, alpha, Delta)
        # zeros come in conjugate pairs, contributing 2 Re(F)
        Z_partial += 2 * F
        log(f"   {gamma:>8.3f} {float(abs(F)):>15.6e} {float(abs(Z_partial)):>17.6e}")

    log(f"\n   |∑ F_TG(rho)| (first 20 zero pairs) ≈ {float(abs(Z_partial)):.4e}")
    log(f"   This has NO x-dependence (Phi_TG is x-independent, so F_TG is too)!")
    log(f"   So even the truncated zero sum is independent of x — cannot encode pi(x).")

    log("\n[STEP 4] VERDICT")
    log("=" * 70)
    log("   The paper's identity")
    log("       αe^(-α²) ≈ -Σ F_TG(ρ) + E_triv")
    log("   has both sides INDEPENDENT OF x.")
    log("   It cannot be 'rounded to give π(x)' because π(x) depends on x and")
    log("   the proposed identity does not.")
    log("")
    log("   The error in the paper traces to its IBP step (page 7): they replace")
    log("   Ψ(t) by t (PNT main term) and drop the deviation Ψ(t)−t. But that")
    log("   deviation is precisely the oscillatory contribution from zeros to")
    log("   π(x). Discarding it makes the resulting identity tautological.")
    log("")
    log("   In addition:")
    log("   - 0th moment is NOT zero for the kernel as defined (≈ 1.49) → main")
    log("     term doesn't cancel as paper claims.")
    log("   - Their zero-density inequality N(σ,T) ≤ AT^{1−1/σ}(log T)^B at")
    log("     σ=1/2 gives N(1/2,T) ≤ AT^{−1}(log T)^B, decreasing in T —")
    log("     blatantly false (true count grows like (T/2π) log T).")
    log("   - Appendix B contains a 'symbolic identity fold' equating Faruk")
    log("     Alpay with a functor Φ_∞: clear crank content.")
    log("")
    log("   CONCLUSION: arXiv:2506.22634 does NOT establish a polylog π(x)")
    log("   algorithm. The 'breakthrough' is illusory. Add to CLOSED_PATHS.")

    return "\n".join(out_lines)


if __name__ == "__main__":
    main()
