"""
Polylog approximation to pi(x) with named-exponent error.

Thread 7 (commit, slot 1, S240). Target: OPEN_POSITIVE_TARGETS.md §P3.

The partial-sum truncation under RH is

  pi_K(x) = R(x) - 2 sum_{j=1..K} Re R(x^{rho_j}),
  rho_j = 1/2 + i gamma_j.

S195 derived (under the Montgomery random-phase heuristic for {gamma_j log x mod 2 pi})

  sigma(K, x) := stddev(pi(x) - pi_K(x))
              ~  sqrt(x) * log K / (pi * sqrt(2 K) * log x).               (*)

This script (a) builds a clean evaluator using the same R_at_rho machinery
as connes_amortisation.py but parameterised for the broader K-policy set
needed for P3, and (b) measures the empirical |error| at x in
{10^5, 10^6, 10^7, 10^8, 10^10} for K policies
{1, log x, log^2 x, log^3 x, log^4 x, x^{1/4}, x^{1/2 cap 8000}}.

Comparing the empirical |error| to (*) gives the named-exponent corollary

  eps(x, K = log^alpha x) ~ Theta( sqrt(x) * log log x / log^{1 + alpha/2} x )

under the random-phase heuristic. For alpha = 6 this is sqrt(x) * log log x / log^4 x;
for alpha = 2 it is sqrt(x) * log log x / log^2 x. Both are strictly better
than the trivial Riemann smooth bound R(x)-pi(x) = O(sqrt(x)).

The script is parameterised by CLI args:

  --xs           comma-separated list of "10^k" exponents (default: 5,6,7,8,10)
  --policies     comma-separated list of policy names (default: all)
  --kmax         hard cap on K, defaults to min(8000, ceil(sqrt(x)))
  --dps          mpmath dps (default 25)
  --M            Mobius truncation in R_at_rho (default 8)
  --csv          path to write per-(x, policy) rows
  --quiet        suppress per-K progress prints

Hardcoded pi(x) at x in {10^5, 10^6, 10^7, 10^8, 10^9, 10^10, 10^12, 10^15}
uses canonical OEIS A006880 values.
"""

from __future__ import annotations

import argparse
import csv
import math
import os
import sys
import time
from dataclasses import dataclass, field
from typing import List, Tuple

import mpmath as mp


# -----------------------------------------------------------------------------
# Canonical pi(10^k) values (OEIS A006880).
# -----------------------------------------------------------------------------

PI_KNOWN = {
    10 ** 1: 4,
    10 ** 2: 25,
    10 ** 3: 168,
    10 ** 4: 1229,
    10 ** 5: 9592,
    10 ** 6: 78498,
    10 ** 7: 664579,
    10 ** 8: 5761455,
    10 ** 9: 50847534,
    10 ** 10: 455052511,
    10 ** 11: 4118054813,
    10 ** 12: 37607912018,
    10 ** 13: 346065536839,
    10 ** 14: 3204941750802,
    10 ** 15: 29844570422669,
}


# -----------------------------------------------------------------------------
# Zeta zeros.
# -----------------------------------------------------------------------------

def get_zeros(K: int, dps: int = 25) -> List[mp.mpf]:
    """Load the imaginary parts gamma_1..gamma_K from the project zeros files.

    Falls back to mpmath.zetazero for uncached counts, but all our K are <= 8000.
    """
    mp.mp.dps = dps
    here = os.path.dirname(os.path.abspath(__file__))
    data_dir = os.path.normpath(os.path.join(here, "..", "..", "..", "data"))
    for size in (300, 500, 1000, 2000, 8000):
        if size >= K:
            fname = os.path.join(data_dir, f"zeta_zeros_{size}.txt")
            if os.path.exists(fname):
                with open(fname) as fh:
                    zs = [line.strip() for line in fh if line.strip()]
                if len(zs) >= K:
                    return [mp.mpf(z) for z in zs[:K]]
    return [mp.zetazero(k).imag for k in range(1, K + 1)]


# -----------------------------------------------------------------------------
# Riemann's R(x) and R(x^rho).
# -----------------------------------------------------------------------------

def mobius_sieve(N: int) -> List[int]:
    mu = [1] * (N + 1)
    is_prime = [True] * (N + 1)
    if N >= 0:
        is_prime[0] = False
    if N >= 1:
        is_prime[1] = False
    for i in range(2, N + 1):
        if is_prime[i]:
            for j in range(i, N + 1, i):
                if j != i:
                    is_prime[j] = False
                mu[j] *= -1
            sq = i * i
            for j in range(sq, N + 1, sq):
                mu[j] = 0
    if N >= 0:
        mu[0] = 0
    return mu


def riemann_R(x: float, M: int = 20, dps: int = 25) -> mp.mpf:
    """Riemann's R(x) = sum_n mu(n)/n * li(x^{1/n}), truncated at n <= M."""
    mp.mp.dps = dps
    if x <= 1:
        return mp.mpf(0)
    mobius = mobius_sieve(M)
    total = mp.mpf(0)
    x_mp = mp.mpf(x)
    for n in range(1, M + 1):
        mu = mobius[n]
        if mu == 0:
            continue
        try:
            val = mp.li(mp.power(x_mp, mp.mpf(1) / n))
            total += mp.mpf(mu) / n * val
        except Exception:
            continue
    return total


def R_at_rho(x: float, gamma: mp.mpf, M: int = 8, dps: int = 25) -> mp.mpc:
    """R(x^rho) for rho = 1/2 + i*gamma using

       R(x^rho) = sum_n mu(n)/n * Ei(rho * ln(x) / n).

    Uses Ei to avoid principal-log phase wrap of x^rho.
    """
    mp.mp.dps = dps
    if x <= 1:
        return mp.mpc(0)
    ln_x = mp.log(mp.mpf(x))
    rho = mp.mpc(mp.mpf("0.5"), gamma)
    mobius = mobius_sieve(M)
    total = mp.mpc(0)
    for n in range(1, M + 1):
        mu = mobius[n]
        if mu == 0:
            continue
        arg = rho * ln_x / n
        if abs(arg) < mp.mpf("1e-20"):
            break
        total += mp.mpc(mu) / n * mp.ei(arg)
    return total


# -----------------------------------------------------------------------------
# Theoretical sigma prediction (S195).
# -----------------------------------------------------------------------------

def sigma_predicted(x: float, K: int) -> float:
    """sigma(K, x) ~ sqrt(x) * log K / (pi * sqrt(2K) * log x).

    Returns +inf for K <= 1.
    """
    if K <= 1 or x <= 1:
        return float("inf")
    return float(
        mp.sqrt(mp.mpf(x)) * mp.log(K)
        / (mp.pi * mp.sqrt(2 * K) * mp.log(x))
    )


# -----------------------------------------------------------------------------
# K-policies.
# -----------------------------------------------------------------------------

@dataclass
class Policy:
    name: str
    K_of: callable  # x: float -> K: int

    def K(self, x: float) -> int:
        return max(1, int(round(self.K_of(x))))


def all_policies() -> List[Policy]:
    return [
        Policy("K=1",      lambda x: 1),
        Policy("logx",     lambda x: math.log(x)),
        Policy("log2x",    lambda x: math.log(x) ** 2),
        Policy("log3x",    lambda x: math.log(x) ** 3),
        Policy("log4x",    lambda x: math.log(x) ** 4),
        Policy("log6x",    lambda x: math.log(x) ** 6),
        Policy("x14",      lambda x: x ** 0.25),
        Policy("x12",      lambda x: x ** 0.5),
    ]


# -----------------------------------------------------------------------------
# Main evaluator.
# -----------------------------------------------------------------------------

@dataclass
class Measurement:
    x: float
    pi_x: int
    R_x: float
    K: int
    pi_K: float
    err: float
    abs_err: float
    sigma_pred: float
    ratio_emp_over_pred: float
    policy_name: str = ""


def evaluate_x(
    x: float,
    pi_x: int,
    K_max: int,
    gammas: List[mp.mpf],
    policies: List[Policy],
    dps: int = 25,
    M: int = 8,
    quiet: bool = False,
) -> List[Measurement]:
    """Compute pi_K(x) for all K = 1..K_max cumulatively and report at policy K's.

    Time = O(K_max) calls to R_at_rho (~5ms each at dps=25, M=8).
    """
    t0 = time.time()
    R_x = riemann_R(x, M=20, dps=dps)
    if not quiet:
        print(f"  R(x={x:.3g}) = {float(R_x):.6f}, pi(x)={pi_x}, "
              f"R-pi={float(R_x) - pi_x:+.4f}", flush=True)

    K_targets = sorted({min(K_max, p.K(x)) for p in policies})
    K_target_to_policies = {k: [p for p in policies if min(K_max, p.K(x)) == k]
                            for k in K_targets}

    results: List[Measurement] = []
    correction = mp.mpf(0)
    K_done = 0
    for K_target in K_targets:
        for k in range(K_done, K_target):
            Rrho = R_at_rho(x, gammas[k], M=M, dps=dps)
            correction += 2 * Rrho.real
        K_done = K_target

        pi_K = R_x - correction
        err = float(pi_K) - pi_x
        sig = sigma_predicted(x, K_target)
        ratio = abs(err) / sig if (sig > 0 and not math.isinf(sig)) else float("nan")

        for pol in K_target_to_policies[K_target]:
            results.append(Measurement(
                x=x, pi_x=pi_x, R_x=float(R_x),
                K=K_target, pi_K=float(pi_K), err=err, abs_err=abs(err),
                sigma_pred=sig, ratio_emp_over_pred=ratio,
                policy_name=pol.name,
            ))
    if not quiet:
        elapsed = time.time() - t0
        print(f"  done x={x:.3g} in {elapsed:.1f}s "
              f"(K_target_max={K_targets[-1]}).", flush=True)
    return results


# -----------------------------------------------------------------------------
# Driver.
# -----------------------------------------------------------------------------

def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--xs", default="5,6,7,8,10",
                   help="Comma list of log10 x exponents.")
    p.add_argument("--policies", default="all",
                   help="Comma list of policy names or 'all'.")
    p.add_argument("--kmax", type=int, default=0,
                   help="Hard cap on K (0 = auto = min(8000, ceil(sqrt(x_max)))).")
    p.add_argument("--dps", type=int, default=25)
    p.add_argument("--M", type=int, default=8)
    p.add_argument("--csv", default=None)
    p.add_argument("--quiet", action="store_true")
    return p.parse_args()


def select_policies(arg: str) -> List[Policy]:
    pols = all_policies()
    if arg == "all":
        return pols
    names = [s.strip() for s in arg.split(",") if s.strip()]
    name_to_pol = {p.name: p for p in pols}
    out = []
    for n in names:
        if n not in name_to_pol:
            raise SystemExit(f"unknown policy: {n}; available: "
                             f"{', '.join(name_to_pol)}")
        out.append(name_to_pol[n])
    return out


def main():
    args = parse_args()
    log10_xs = [int(s) for s in args.xs.split(",")]
    xs: List[Tuple[float, int]] = []
    for k in log10_xs:
        x = 10.0 ** k
        x_int = 10 ** k
        if x_int not in PI_KNOWN:
            raise SystemExit(f"pi(10^{k}) not in PI_KNOWN; add it.")
        xs.append((x, PI_KNOWN[x_int]))

    policies = select_policies(args.policies)
    K_needed = max(p.K(x) for x, _ in xs for p in policies)
    K_max = args.kmax if args.kmax > 0 else min(8000, K_needed)
    K_needed_capped = min(K_needed, K_max)

    print(f"Loading {K_needed_capped} zeta zeros (mp dps={args.dps})...",
          flush=True)
    t0 = time.time()
    gammas = get_zeros(K_needed_capped, dps=args.dps)
    print(f"  loaded {len(gammas)} zeros in {time.time()-t0:.1f}s "
          f"(gamma_1={float(gammas[0]):.4f}, "
          f"gamma_K={float(gammas[-1]):.4f}).", flush=True)

    rows: List[Measurement] = []
    for x, pi_x in xs:
        print(f"\n[x = {x:.3g}]", flush=True)
        ms = evaluate_x(x, pi_x, K_max, gammas, policies,
                        dps=args.dps, M=args.M, quiet=args.quiet)
        rows.extend(ms)

    # ----------- Print summary table -----------
    print("\n" + "=" * 92)
    print(f"{'x':>10s} {'policy':>8s} {'K':>7s} {'sigma':>10s} "
          f"{'|err|':>10s} {'eps/sx':>11s} {'sgima/sx':>10s} {'r_em':>7s}")
    print("-" * 92)
    for m in rows:
        sx = math.sqrt(m.x)
        print(f"{m.x:>10.2g} {m.policy_name:>8s} {m.K:>7d} "
              f"{m.sigma_pred:>10.3f} {m.abs_err:>10.3f} "
              f"{m.abs_err / sx:>11.3e} {m.sigma_pred / sx:>10.3e} "
              f"{m.ratio_emp_over_pred:>7.3f}")

    # ----------- Print named-exponent corollary -----------
    print("\n" + "=" * 92)
    print("NAMED-EXPONENT COROLLARY (under Montgomery random-phase heuristic):")
    print("  K = log^alpha(x) zeros => sigma(x, K) ~ sqrt(x) * log log x"
          " / log^{1 + alpha/2} x")
    print()
    print(f"{'alpha':>6s} {'K formula':>15s} {'x=1e6':>10s} {'x=1e8':>10s} "
          f"{'x=1e10':>10s} {'x=1e12':>10s} {'x=1e15':>10s}")
    for alpha in [2, 3, 4, 6, 8]:
        K_at = lambda lx: max(1, int(round(lx ** alpha)))
        sigs = []
        for log10x in [6, 8, 10, 12, 15]:
            x_ = 10.0 ** log10x
            lx = math.log(x_)
            K = K_at(lx)
            sig = sigma_predicted(x_, K)
            sigs.append(sig)
        print(f"{alpha:>6d} {'log^'+str(alpha)+' x':>15s} "
              + " ".join(f"{s:>10.3g}" for s in sigs))
    print()
    print("  Concrete sigma / sqrt(x) at K = log^alpha x:")
    print(f"{'alpha':>6s}  {'x=1e6':>10s} {'x=1e8':>10s} {'x=1e10':>10s} "
          f"{'x=1e12':>10s} {'x=1e15':>10s}")
    for alpha in [2, 3, 4, 6, 8]:
        ratios = []
        for log10x in [6, 8, 10, 12, 15]:
            x_ = 10.0 ** log10x
            lx = math.log(x_)
            K = max(1, int(round(lx ** alpha)))
            sig = sigma_predicted(x_, K)
            ratios.append(sig / math.sqrt(x_))
        print(f"{alpha:>6d}  " + " ".join(f"{r:>10.3e}" for r in ratios))

    # ----------- CSV output -----------
    if args.csv:
        csv_path = args.csv if os.path.isabs(args.csv) else os.path.join(
            os.path.dirname(os.path.abspath(__file__)), args.csv)
        with open(csv_path, "w", newline="") as fh:
            w = csv.writer(fh)
            w.writerow([
                "x", "pi_x", "R_x", "policy", "K",
                "pi_K", "err", "abs_err", "sigma_pred",
                "ratio_emp_over_pred", "abs_err_over_sqrt_x",
                "sigma_pred_over_sqrt_x",
            ])
            for m in rows:
                sx = math.sqrt(m.x)
                w.writerow([
                    f"{m.x:.6g}", m.pi_x, f"{m.R_x:.6f}", m.policy_name, m.K,
                    f"{m.pi_K:.6f}", f"{m.err:.6f}", f"{m.abs_err:.6f}",
                    f"{m.sigma_pred:.6f}",
                    f"{m.ratio_emp_over_pred:.6f}",
                    f"{m.abs_err / sx:.6e}",
                    f"{m.sigma_pred / sx:.6e}",
                ])
        print(f"\nWrote CSV: {csv_path}")


if __name__ == "__main__":
    main()
