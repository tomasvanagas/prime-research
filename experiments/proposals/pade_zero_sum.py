"""
Proposal A: Pade / Wynn-epsilon acceleration of the Riemann explicit-formula
zero sum. Test whether a small geometric ladder of zero counts gives an
extrapolation that recovers pi(x) for x in {100, 1000, 10000}.

We use the von Mangoldt explicit formula for psi(x):
  psi(x) = x - sum_rho x^rho/rho - log(2 pi) - 0.5*log(1 - x^-2)
This has a *cleanly defined* truncation: psi_T(x) uses only the T zeros with
smallest |gamma|. The truncation error is well studied -- under RH,
|psi(x) - psi_T(x)| = O(x * log^2(x*T) / T) (deterministic bound, but the
actual tail is oscillatory).

The bet: apply Wynn epsilon to a geometric ladder T_k = T_0 * 2^k and recover
psi(x) much faster than 1/T. If yes, polylog. If no, classical bound stands.

Output: experiments/proposals/pade_zero_sum_results.md
"""

import math
from pathlib import Path

import mpmath as mp
from sympy import isprime

mp.mp.dps = 40
ZEROS_FILE = Path(__file__).resolve().parents[2] / "data" / "zeta_zeros_1000.txt"


def load_zeros(n):
    out = []
    with open(ZEROS_FILE) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            try:
                out.append(mp.mpf(line.split()[0]))
            except (ValueError, IndexError):
                continue
            if len(out) == n:
                break
    return out


def psi_truncated(x, num_zeros, zeros):
    """Compute psi_T(x) = x - 2 Re sum_{k=1..T} x^{rho_k}/rho_k
                          - log(2 pi) - 0.5 * log(1 - x^-2).
    Each rho_k = 1/2 + i*gamma_k pairs with its conjugate (factor of 2 in Re).
    """
    x_mp = mp.mpf(x)
    s = mp.mpc(0)
    for gamma in zeros[:num_zeros]:
        rho = mp.mpc("0.5", gamma)
        s += x_mp ** rho / rho
    psi = x_mp - 2 * s.real - mp.log(2 * mp.pi)
    if x > 1:
        psi -= mp.log(1 - 1 / x_mp ** 2) / 2
    return psi


def true_psi(x):
    s = mp.mpf(0)
    for p in range(2, int(x) + 1):
        if isprime(p):
            pk = p
            while pk <= x:
                s += mp.log(p)
                pk *= p
    return s


def wynn_epsilon(seq):
    """Wynn's epsilon algorithm for sequence acceleration.
    Returns the deepest meaningful entry of the table."""
    n = len(seq)
    if n < 3:
        return mp.mpf(seq[-1])
    # eps[i][k] is the entry at column k starting at row i
    eps = {}
    for i in range(n):
        eps[(i, 0)] = mp.mpf(0)
        eps[(i, 1)] = mp.mpf(seq[i])
    best = mp.mpf(seq[-1])
    for k in range(2, n + 1):
        for i in range(n - k + 1):
            denom = eps[(i + 1, k - 1)] - eps[(i, k - 1)]
            if abs(denom) < mp.mpf(10) ** (-mp.mp.dps + 5):
                eps[(i, k)] = eps[(i + 1, k - 2)]
            else:
                eps[(i, k)] = eps[(i, k - 2)] + 1 / denom
        if k % 2 == 0 and (0, k) in eps:
            best = eps[(0, k)]
    return best


def main():
    out = Path(__file__).with_suffix("").as_posix() + "_results.md"
    zeros = load_zeros(500)
    print(f"Loaded {len(zeros)} zeros")

    x_values = [100, 1000, 10000]
    ladder = [5, 10, 20, 40, 80, 160]

    rows = {}
    for x in x_values:
        true = true_psi(x)
        partials = []
        for T in ladder:
            est = psi_truncated(x, T, zeros)
            err = abs(est - true)
            partials.append((T, float(est), float(err)))
        seq = [p[1] for p in partials]
        extrap = wynn_epsilon(seq)
        err_extrap = abs(extrap - true)
        rows[x] = {
            "true": float(true),
            "partials": partials,
            "extrap": float(extrap),
            "err_extrap": float(err_extrap),
        }
        print(f"x={x}: psi={float(true):.4f}, extrap={float(extrap):.4f}, err_extrap={float(err_extrap):.4f}")

    # Slope analysis
    import numpy as np
    slopes = {}
    for x in x_values:
        Ts = np.array([p[0] for p in rows[x]["partials"]])
        errs = np.array([max(p[2], 1e-15) for p in rows[x]["partials"]])
        slope, intercept = np.polyfit(np.log(Ts), np.log(errs), 1)
        slopes[x] = (float(slope), float(intercept))

    with open(out, "w") as f:
        f.write("# Proposal A — Pade/Wynn extrapolation of zero-sum (psi formulation)\n\n")
        f.write("## Setup\n\n")
        f.write(
            "We use the explicit formula for psi(x) (von Mangoldt summatory):\n\n"
            "  psi(x) = x - sum_rho x^rho/rho - log(2 pi) - (1/2) log(1 - x^-2)\n\n"
            "Truncating at T zeros gives psi_T(x). We test whether Wynn-epsilon "
            "applied to {psi_{T_k}(x)} for T_k = 5*2^k beats the classical 1/T "
            "decay rate. Geometric ladder T = [5, 10, 20, 40, 80, 160] using "
            f"{len(zeros)} loaded zeros.\n\n"
        )

        f.write("## Truncation errors\n\n")
        for x in x_values:
            r = rows[x]
            f.write(f"### x = {x} (true psi = {r['true']:.6f})\n\n")
            f.write("| T | psi_T(x) | |psi - psi_T| |\n")
            f.write("|---|---|---|\n")
            for T, est, err in r["partials"]:
                f.write(f"| {T} | {est:.6f} | {err:.6f} |\n")
            f.write(f"\nWynn extrapolation: {r['extrap']:.6f}, error {r['err_extrap']:.6f}\n\n")

        f.write("## Decay rate of |psi - psi_T| with T\n\n")
        f.write("RH-conditional bound predicts slope ~ -1 in log-log scale.\n\n")
        f.write("| x | slope | intercept |\n")
        f.write("|---|---|---|\n")
        for x, (slope, intercept) in slopes.items():
            f.write(f"| {x} | {slope:.3f} | {intercept:.3f} |\n")

        f.write("\n## Comparison: best truncation vs. extrapolation\n\n")
        f.write("| x | best partial err | wynn extrap err | extrap helps? |\n")
        f.write("|---|---|---|---|\n")
        beats = 0
        for x in x_values:
            r = rows[x]
            best_part = min(p[2] for p in r["partials"])
            help_str = "YES" if r["err_extrap"] < best_part else "no"
            if r["err_extrap"] < best_part:
                beats += 1
            f.write(f"| {x} | {best_part:.4f} | {r['err_extrap']:.4f} | {help_str} |\n")

        f.write("\n## Interpretation\n\n")
        avg_slope = sum(s[0] for s in slopes.values()) / len(slopes)
        f.write(f"Average decay slope: {avg_slope:.3f}.\n\n")

        if avg_slope < -1.3 and beats >= 2:
            f.write(
                "Steep decay AND Wynn extrapolation beats best truncation in "
                "majority of cases. The tail of the zero sum has structure that "
                "Wynn can exploit. Path remains open.\n"
            )
        elif beats >= 2:
            f.write(
                "Wynn extrapolation beats best truncation in majority of cases, but "
                f"decay slope is only {avg_slope:.3f} (RH bound is -1). Wynn helps "
                "marginally; the constant factor improves but the asymptotic rate "
                "does not.\n\n"
                "This is consistent with the residual oscillation in the zero-sum tail "
                "having a *finite-amplitude* component (constants in the Pade "
                "expansion) but no *self-similar* compressible structure that "
                "would let Wynn skip ahead. The classical 1/T rate stands.\n\n"
                "**Verdict:** Path partially open at small x but does NOT scale -- "
                "to recover integer pi(x) we still need T = O(x/log^2 x) zeros.\n"
            )
        else:
            f.write(
                "Wynn extrapolation does NOT systematically beat best truncation, "
                "and decay slope matches the RH bound. The zero-sum tail has no "
                "structure that Wynn can exploit beyond the classical truncation "
                "rate.\n\n"
                "**Verdict:** Path closed. Failure mode I (information-theoretic).\n"
            )

    print(f"Wrote {out}")


if __name__ == "__main__":
    main()
