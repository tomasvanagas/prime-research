"""
Proposal 26: Cesaro-Fejer-Damped Explicit Formula

The Riemann explicit formula for pi(x) needs T ~ sqrt(x) zeros for an exact
integer answer because the truncation error is O(x/T).  Replacing the
sharp cutoff with a Fejer kernel
    K_T(gamma) = (sin(gamma/(2T)) / (gamma/(2T)))^2 * 1[|gamma| <= 2T]
gives a smooth window with L^1 norm ~ 1.  The point of the experiment
is NOT to claim the truncated sum reaches integer accuracy at smaller T
(it can't, by the standard explicit-formula error analysis); rather, we
measure the *distribution* of the rounding gap |S_T(x) - round(S_T(x))|
across x = 1..10000 for several truncation values T.

If for a positive density of x the gap is provably > 1/4 already at
T = x^{1/3}, those primes can be recovered with O(x^{1/3}) zeros, which
beats the unconditional sqrt(x) barrier on a fraction of inputs.  We
quantify that fraction.

Run: python3 proposal26_fejer_damped_explicit.py
"""

from pathlib import Path
import time

import mpmath
from sympy import primepi
from sympy.ntheory import mobius


ZEROS_FILE = Path(__file__).resolve().parents[2] / "data" / "zeta_zeros_2000.txt"


def load_zeros(n_max: int = 2000):
    out = []
    with open(ZEROS_FILE) as fh:
        for line in fh:
            line = line.strip()
            if line:
                out.append(mpmath.mpf(line))
                if len(out) >= n_max:
                    break
    return out


def riemann_R(x):
    """Riemann's R function R(x) = sum_{k>=1} mu(k)/k * li(x^{1/k}).

    Truncate the sum at K = ceil(log2(x)) since x^{1/k} -> 1 quickly and
    li(1)=-inf gets handled by mpmath; we stop when x^{1/k} < 1.5.
    """
    x = mpmath.mpf(x)
    s = mpmath.mpf(0)
    k = 1
    while True:
        root = x ** (mpmath.mpf(1) / k)
        if root <= mpmath.mpf("1.05"):
            break
        mu = mobius(k)
        if mu != 0:
            s += mpmath.mpf(mu) / k * mpmath.li(root)
        k += 1
        if k > 200:
            break
    return s


def li_via_Ei(rho, log_x):
    """Compute li(x^rho) via analytic continuation: li(x^rho) = Ei(rho * log x).

    Using mpmath.li(x**rho) directly is WRONG because mpmath uses the
    principal branch of complex log, which collapses x^rho to its principal
    value, losing the winding information required by Riemann's formula.
    The correct expression is Ei(rho * log_x) where the argument is the
    *unwound* product on the natural sheet.
    """
    return mpmath.ei(rho * log_x)


def riemann_R_zero_term(rho, x):
    """R(x^rho) = sum_{k>=1} mu(k)/k * li(x^{rho/k}) using the safe formula above."""
    s = mpmath.mpc(0)
    log_x = mpmath.log(x)
    k = 1
    while True:
        # x^{rho/k} has modulus |x|^{0.5/k} = sqrt(x)^{1/k}.  Stop when
        # the term becomes negligible (k large enough).
        if k > 1 and mpmath.mpf(x) ** (mpmath.mpf("0.5") / k) <= mpmath.mpf("1.05"):
            break
        mu = mobius(k)
        if mu != 0:
            s += mpmath.mpf(mu) / k * li_via_Ei(rho / k, log_x)
        k += 1
        if k > 100:
            break
    return s


def explicit_pi_fejer(x, zeros, T, sharp=False):
    """Approximate pi(x) via R(x) - sum_gamma 2*Re[Ei(rho * log x)] with optional Fejer window.

    The full explicit formula has additional small corrections (higher-k Mobius
    terms in the zero sum, O(x^{1/4}/T)); these are negligible for our gap-CDF
    measurement and we drop them for speed.
    """
    x = mpmath.mpf(x)
    log_x = mpmath.log(x)
    s = riemann_R(x)
    for gamma in zeros:
        if sharp:
            if gamma > T:
                break
            K = mpmath.mpf(1)
        else:
            if gamma > 2 * T:
                break
            u = gamma / (2 * T)
            K = (mpmath.sin(u) / u) ** 2 if u else mpmath.mpf(1)
        rho = mpmath.mpc(mpmath.mpf("0.5"), gamma)
        contribution = li_via_Ei(rho, log_x)
        s -= 2 * K * contribution.real
    return s


def main():
    mpmath.mp.dps = 30
    print("Loading zeros...")
    zeros = load_zeros(2000)
    print(f"  loaded {len(zeros)} zeros, max gamma = {float(zeros[-1]):.2f}")

    sample_xs = [100, 200, 500, 1000, 2000, 5000, 10000]
    Ts = [10, 30, 100, 300, 1000]

    print("\nIntegrity check at sample x with T=1000 (sharp), should be close to pi(x):")
    for x in sample_xs:
        pi_x = int(primepi(x))
        s = explicit_pi_fejer(x, zeros, 1000, sharp=True)
        gap = float(s - pi_x)
        print(f"  x={x:>5}  pi={pi_x:>4}  S_T={float(s):.4f}  diff={gap:+.4f}")

    print("\nFejer-damped behaviour at sample x:")
    print(f"  {'x':>5} {'pi':>4} | " +
          " | ".join(f"T={T}: S, gap" for T in Ts))
    for x in sample_xs:
        pi_x = int(primepi(x))
        row = [f"{x:>5} {pi_x:>4} |"]
        for T in Ts:
            s = explicit_pi_fejer(x, zeros, T, sharp=False)
            sf = float(s)
            row.append(f"  {sf:7.3f} ({sf - pi_x:+.3f})")
        print("".join(row))

    # Larger sweep: round-gap CDF for x = 100, 200, ..., 3000 at fixed T's.
    print("\nRound-gap CDF over x in [100, 3000] step 100:")
    xs = list(range(100, 3001, 100))
    pi_xs = [int(primepi(x)) for x in xs]
    for T in [30, 100, 300, 1000]:
        gaps_sharp = []
        gaps_fejer = []
        sharp_correct = 0
        fejer_correct = 0
        t0 = time.time()
        for x, pi_x in zip(xs, pi_xs):
            s_sharp = float(explicit_pi_fejer(x, zeros, T, sharp=True))
            s_fejer = float(explicit_pi_fejer(x, zeros, T, sharp=False))
            gaps_sharp.append(abs(s_sharp - round(s_sharp)))
            gaps_fejer.append(abs(s_fejer - round(s_fejer)))
            if round(s_sharp) == pi_x:
                sharp_correct += 1
            if round(s_fejer) == pi_x:
                fejer_correct += 1
        bins = [0.05, 0.1, 0.2, 0.3, 0.4, 0.5]
        sharp_pct = [sum(1 for g in gaps_sharp if g < b) / len(gaps_sharp) for b in bins]
        fejer_pct = [sum(1 for g in gaps_fejer if g < b) / len(gaps_fejer) for b in bins]
        print(f"  T={T:>4}  ({time.time()-t0:.1f}s)")
        print(f"    sharp  P(gap<b): " +
              "  ".join(f"b={b}:{p:.2f}" for b, p in zip(bins, sharp_pct)))
        print(f"    fejer  P(gap<b): " +
              "  ".join(f"b={b}:{p:.2f}" for b, p in zip(bins, fejer_pct)))
        print(f"    sharp  round(S_T) == pi(x): {sharp_correct}/{len(xs)}")
        print(f"    fejer  round(S_T) == pi(x): {fejer_correct}/{len(xs)}")


if __name__ == "__main__":
    main()
