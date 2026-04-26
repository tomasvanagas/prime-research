"""
FOCUS-4: Borel-Padé + Cesàro-Fejér hybrid for π(x).

The Riemann explicit formula π(x) ≈ R(x) − 2·Σ_γ Re(R(x^ρ)) needs O(√x) zeros
for an integer-accurate answer.  Two earlier sessions tested two
*independent* convergence-acceleration tricks:

    S45  Borel-Padé regularisation  (experiments/wildcard/borel_explicit_formula.py)
    S46  Cesàro-Fejér damped window (experiments/proposals/proposal26_*.py)

This script applies them **sequentially**: Fejér-window the per-zero
contributions FIRST, then Borel-Padé the resulting partial-sum sequence.
The two interventions have independent failure modes (Fejér halves the L^1
truncation constant; Borel-Padé reshapes the series into a Borel-summable
asymptotic), so combining them is the natural next test.

Bug to avoid (S46 / EDGES.md §E3.6):
    use mpmath.ei(rho * log(x)),  NEVER mpmath.li(x**rho).

Run:  python3 borel_fejer_hybrid.py
"""

from __future__ import annotations

from pathlib import Path
import argparse
import time

import mpmath
from mpmath import mp, mpf, mpc, polyval, quad
from sympy import primepi
from sympy.ntheory import mobius

ZEROS_FILE = Path(__file__).resolve().parents[2] / "data" / "zeta_zeros_2000.txt"


def load_zeros(n_max: int = 2000) -> list[mpf]:
    out: list[mpf] = []
    with open(ZEROS_FILE) as fh:
        for line in fh:
            line = line.strip()
            if line:
                out.append(mpmath.mpf(line))
                if len(out) >= n_max:
                    break
    return out


def riemann_R(x) -> mpf:
    """Riemann's R(x) = Σ_k μ(k)/k · li(x^{1/k})."""
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


def li_via_Ei(rho: mpc, log_x: mpf) -> mpc:
    """li(x^ρ) via analytic continuation Ei(ρ · log x) — branch-correct version."""
    return mpmath.ei(rho * log_x)


def zero_term_pi(gamma: mpf, log_x: mpf) -> mpf:
    """2·Re(li(x^ρ)) for a single zero ρ = 1/2 + iγ — the per-zero contribution."""
    rho = mpc(mpf("0.5"), gamma)
    return 2 * li_via_Ei(rho, log_x).real


def fejer_weight(gamma: mpf, T: mpf) -> mpf:
    """Cesàro-Fejér kernel:  (sin(γ/(2T)) / (γ/(2T)))²  on |γ| ≤ 2T, else 0."""
    if gamma > 2 * T:
        return mpf(0)
    u = gamma / (2 * T)
    if u == 0:
        return mpf(1)
    s = mpmath.sin(u) / u
    return s * s


def explicit_pi(x, zeros, T, mode: str) -> mpf:
    """Compute π(x) approx via R(x) − Σ_γ w(γ)·2·Re(li(x^ρ)).

    mode: 'sharp'  — w(γ) = 1[γ ≤ T]
          'fejer'  — w(γ) = K_T(γ)
          'borel_fejer' — Fejér-windowed increments + Borel-Padé tail.
    """
    x_mp = mpf(x)
    log_x = mpmath.log(x_mp)
    R = riemann_R(x_mp)

    if mode == "sharp":
        s = mpf(0)
        for gamma in zeros:
            if gamma > T:
                break
            s += zero_term_pi(gamma, log_x)
        return R - s

    if mode == "fejer":
        s = mpf(0)
        for gamma in zeros:
            if gamma > 2 * T:
                break
            s += fejer_weight(gamma, T) * zero_term_pi(gamma, log_x)
        return R - s

    if mode == "borel_fejer":
        increments = []
        for gamma in zeros:
            if gamma > 2 * T:
                break
            increments.append(fejer_weight(gamma, T) * zero_term_pi(gamma, log_x))
        if not increments:
            return R
        s_borel = borel_pade_sum_best(increments)
        if s_borel is None:
            s_borel = sum(increments, mpf(0))
        return R - s_borel

    raise ValueError(f"unknown mode {mode!r}")


def borel_pade_sum_best(increments: list[mpf]) -> mpf | None:
    """Borel-Padé on a finite increment sequence; pick the (M,N) order with
    smallest residual variance among the last few orders.

    Returns the resummed Σ a_k or None if all (M,N) orders fail.
    """
    if len(increments) < 6:
        return sum(increments, mpf(0))

    # Borel transform: B(z) = Σ a_k z^k / k!
    bcoefs = []
    fact = mpf(1)
    for k, a in enumerate(increments):
        if k > 0:
            fact *= k
        bcoefs.append(a / fact)

    # Try a few low/medium Padé orders and aggregate.
    candidates: list[tuple[mpf, tuple[int, int]]] = []
    L = len(bcoefs)
    orders = []
    for total in (10, 14, 20, 28):
        if total + 1 > L:
            continue
        for M in (total // 2,):
            N = total - M
            orders.append((M, N))
    # Always try (5,5) and (10,10) when available.
    for mn in [(5, 5), (10, 10), (15, 15)]:
        if sum(mn) + 1 <= L and mn not in orders:
            orders.append(mn)

    for (M, N) in orders:
        try:
            p, q = mp.pade(bcoefs, M, N)
        except Exception:
            continue

        def integrand(z, p=p, q=q):
            try:
                num = polyval(p[::-1], z)
                den = polyval(q[::-1], z)
                if abs(den) < mpf("1e-30"):
                    return mpf(0)
                return mpmath.exp(-z) * num / den
            except Exception:
                return mpf(0)

        try:
            val = quad(integrand, [0, mp.inf])
        except Exception:
            continue
        if val is None or not mpmath.isfinite(val):
            continue
        candidates.append((val, (M, N)))

    if not candidates:
        return None

    # Median-of-candidates: more robust to one rogue (M,N) divergence than min/mean.
    vals_sorted = sorted(c[0] for c in candidates)
    return vals_sorted[len(vals_sorted) // 2]


def sweep(zeros, xs, Ts, modes, dps_local=30, verbose=True):
    """Return dict[(mode, T)] -> list of (x, pi_x, S, round_correct, gap)."""
    mp.dps = dps_local
    pi_xs = {x: int(primepi(x)) for x in xs}
    results: dict[tuple[str, int], list[tuple[int, int, float, bool, float]]] = {}
    for mode in modes:
        for T in Ts:
            t0 = time.time()
            row = []
            for x in xs:
                S = explicit_pi(x, zeros, mpf(T), mode)
                Sf = float(S)
                rounded = round(Sf)
                gap = abs(Sf - rounded)
                correct = (rounded == pi_xs[x])
                row.append((x, pi_xs[x], Sf, correct, gap))
            results[(mode, T)] = row
            if verbose:
                hits = sum(1 for r in row if r[3])
                avg_gap = sum(r[4] for r in row) / len(row)
                print(
                    f"  {mode:>13}  T={T:>4}  "
                    f"hits {hits}/{len(xs)}   "
                    f"avg-gap {avg_gap:.3f}   "
                    f"({time.time()-t0:.1f}s)"
                )
    return results


def fmt_table(results, xs, Ts, modes):
    header = f"{'x':>7}  {'π(x)':>6}  | " + " | ".join(
        f"{mode}@T={T}" for T in Ts for mode in modes
    )
    rows = [header, "-" * len(header)]
    for i, x in enumerate(xs):
        pi_x = results[(modes[0], Ts[0])][i][1]
        cells = []
        for T in Ts:
            for mode in modes:
                _, _, S, correct, gap = results[(mode, T)][i]
                marker = "✓" if correct else "✗"
                cells.append(f"{S:8.2f}{marker}({gap:.2f})")
        rows.append(f"{x:>7}  {pi_x:>6}  | " + " | ".join(cells))
    return "\n".join(rows)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--quick", action="store_true",
                    help="Reduced sweep for smoke testing.")
    args = ap.parse_args()

    mp.dps = 30
    print("Loading zeta zeros...")
    zeros = load_zeros(2000)
    print(f"  loaded {len(zeros)} zeros (max γ = {float(zeros[-1]):.1f})")

    if args.quick:
        xs = [1000, 5000, 20000]
        Ts = [100, 300]
    else:
        # 8 x-values spanning 10^3..10^5; 2000 zeros enough for T ≤ 1000.
        xs = [1000, 2500, 5000, 10000, 25000, 50000, 75000, 100000]
        Ts = [50, 100, 300, 1000]

    modes = ["sharp", "fejer", "borel_fejer"]

    print(f"\nx-values:   {xs}")
    print(f"T-values:   {Ts}")
    print(f"modes:      {modes}\n")

    print("Sweep:")
    results = sweep(zeros, xs, Ts, modes, dps_local=30)

    print("\nDetailed table (S = approx, ✓/✗ = round matches π(x), gap = |S−round(S)|):\n")
    print(fmt_table(results, xs, Ts, modes))

    # Aggregate hit-rate
    print("\n--- Recovery rate (round(S) == π(x)) ---")
    print(f"{'mode':>14} | " + " | ".join(f"T={T}" for T in Ts))
    for mode in modes:
        cells = []
        for T in Ts:
            row = results[(mode, T)]
            hits = sum(1 for r in row if r[3])
            cells.append(f"{hits:>2}/{len(row):<2}")
        print(f"{mode:>14} | " + " | ".join(cells))

    print("\n--- Mean rounding gap |S − round(S)| ---")
    print(f"{'mode':>14} | " + " | ".join(f"T={T}" for T in Ts))
    for mode in modes:
        cells = []
        for T in Ts:
            row = results[(mode, T)]
            avg = sum(r[4] for r in row) / len(row)
            cells.append(f"{avg:.3f}")
        print(f"{mode:>14} | " + " | ".join(cells))

    # Failure clustering: at the strongest T where any mode wins, list which xs
    # the hybrid uniquely solves.
    print("\n--- Mode-by-mode wins ---")
    for T in Ts:
        sharp = {r[0] for r in results[("sharp", T)] if r[3]}
        fejer = {r[0] for r in results[("fejer", T)] if r[3]}
        bofe = {r[0] for r in results[("borel_fejer", T)] if r[3]}
        only_hybrid = sorted(bofe - fejer - sharp)
        only_fejer = sorted(fejer - sharp - bofe)
        regr_hybrid = sorted(fejer - bofe)  # Fejer solves but hybrid loses
        print(
            f"  T={T:>4}: hybrid-only x={only_hybrid}  "
            f"fejer-only x={only_fejer}  "
            f"hybrid-regressions x={regr_hybrid}"
        )


if __name__ == "__main__":
    main()
