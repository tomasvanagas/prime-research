"""
EXTENSION OF TASK #4 (Session 33): K_min with 2000 zeros, resolve x=10^4 anomaly
================================================================================

Session 33 reported a persistent +1 error at x = 10^4 with all 1000 zeta zeros
available. The classical truncation-error bound under RH is

    |pi(x) - pi_K(x)|  =  O(x^{1/2} * log(x) / T_K)

so at x = 10^4 with T = 1419 (the 1000-th gamma) the bound is ~0.65, which is
right on the +/-0.5 rounding cliff. Doubling to 2000 zeros (max gamma = 2515)
should drop the bound to ~0.37 -- inside the rounding window.

This experiment:
  (i)   resolves the x = 10^4 persistent error: with 2000 zeros, does the
        rounded analytic value finally settle on 1229?
  (ii)  reports K_min(x) (smallest K with round(pi_K) == pi(x)) and a stronger
        K_min* (first K stable for >= 50 consecutive successes) at
        x in {10^3, 10^4, 10^5, 10^6}.
  (iii) checks the residual at x = 10^7, where K_min(x) is still beyond 2000.

Implementation note. We compute the entire pi_K(x) trajectory for K = 1..2000
in a single O(K_max) pass, accumulating the zero sum incrementally.

Save: alongside this script as k_min_extended_results.md.
"""

import math
import os
import sys
import time

import mpmath


DATA_DIR = "/apps/aplikacijos/prime-research/data"
DPS = 30


def load_zeros():
    path = os.path.join(DATA_DIR, "zeta_zeros_2000.txt")
    with open(path) as f:
        return [mpmath.mpf(line.strip()) for line in f if line.strip()]


_MU = [0, 1, -1, -1, 0, -1, 1, -1, 0, 0, 1, -1, 0, -1, 1, 1, 0, -1, 0, -1, 0,
       1, 1, -1, 0, 0, 1, 0, 0, -1, -1, -1, 0, 1, 1, 1, 0, -1, 1, 1, 0, -1,
       -1, -1, 0, 0, -1, -1, 0, 0, 0, 1, -1, 0, 0, -1, 1, 0, 1, 1, -1, 0, -1,
       1, 0, 0, 1, -1, 0, 1, -1, -1, 0, -1, 1, 0, 0, 1, -1, 0, 0, 0, -1, 1,
       -1, 0, 1, 1, -1, 0, 0, -1, 0, 0, -1, -1, 1, 0, 1, 0, 0]


def R_mp(x, dps=DPS):
    with mpmath.workdps(dps):
        x = mpmath.mpf(x)
        if x <= 1:
            return mpmath.mpf(0)
        s = mpmath.mpf(0)
        for n in range(1, len(_MU)):
            if _MU[n] == 0:
                continue
            xn = mpmath.power(x, mpmath.mpf(1) / n)
            if xn <= mpmath.mpf("1.0001"):
                break
            s += mpmath.mpf(_MU[n]) / n * mpmath.li(xn)
        return s


def pi_trajectory(x, zeros, dps=DPS):
    """Return list pi_K(x) for K = 0, 1, 2, ..., len(zeros), where K=0 is the
    main term R(x) alone (no zero corrections). Uses incremental zero-sum
    accumulation for O(K_max) cost."""
    with mpmath.workdps(dps):
        x_mp = mpmath.mpf(x)
        lnx = mpmath.log(x_mp)
        main = R_mp(x, dps)

        zero_sum = mpmath.mpf(0)
        out = [float(mpmath.re(main))]  # K = 0
        for gamma in zeros:
            rho = mpmath.mpc(0.5, gamma)
            r1 = mpmath.ei(rho * lnx)
            r2 = mpmath.ei(rho / 2 * lnx)
            r3 = mpmath.ei(rho / 3 * lnx)
            r5 = mpmath.ei(rho / 5 * lnx)
            r_rho = r1 - r2 / 2 - r3 / 3 - r5 / 5
            zero_sum += 2 * mpmath.re(r_rho)
            out.append(float(mpmath.re(main - zero_sum)))
        return out


def lucy_pi(x):
    x = int(x)
    if x < 2:
        return 0
    sx = int(x ** 0.5)
    while (sx + 1) ** 2 <= x:
        sx += 1
    while sx * sx > x:
        sx -= 1

    small = list(range(-1, sx + 1))
    large = [0] * (sx + 2)
    for v in range(1, sx + 1):
        large[v] = x // v - 1
    for p in range(2, sx + 1):
        if small[p] == small[p - 1]:
            continue
        pcnt = small[p - 1]
        p2 = p * p
        for v in range(1, min(sx, x // p2) + 1):
            d = v * p
            if d <= sx:
                large[v] -= large[d] - pcnt
            else:
                large[v] -= small[x // d] - pcnt
        for v in range(sx, p2 - 1, -1):
            small[v] -= small[v // p] - pcnt
    return large[1]


def k_min_metrics(traj, exact, run_after=50):
    """Given pi_K trajectory and exact value, return:
        K_min  -- first K with round(pi_K) == exact (or None)
        K_min* -- first K such that round(pi_K')==exact for K' = K..K+run_after
                  (or None)
    K = 0 means main term only (no zeros)."""
    K_min = None
    K_stable = None
    streak_start = None
    streak_len = 0
    for K in range(len(traj)):
        if round(traj[K]) == exact:
            if K_min is None:
                K_min = K
            if streak_start is None:
                streak_start = K
            streak_len += 1
            if streak_len >= run_after and K_stable is None:
                K_stable = streak_start
        else:
            streak_start = None
            streak_len = 0
    return K_min, K_stable


def main():
    print("=" * 80)
    print("EXTENSION OF TASK #4: K_min with 2000 zeros (incremental)")
    print("=" * 80, flush=True)

    t0 = time.time()
    zeros = load_zeros()
    print(f"Loaded {len(zeros)} zeros (max gamma = {float(zeros[-1]):.3f})  "
          f"[{time.time() - t0:.2f}s]", flush=True)

    x_values = [10 ** 3, 10 ** 4, 10 ** 5, 10 ** 6, 10 ** 7]
    print("\n--- Computing exact pi(x) via Lucy DP ---", flush=True)
    exact = {}
    for x in x_values:
        t0 = time.time()
        exact[x] = lucy_pi(x)
        print(f"  pi({x:.0e}) = {exact[x]:>10}   [{time.time() - t0:.3f}s]",
              flush=True)

    # Compute trajectories.
    trajectories = {}
    print("\n--- Computing pi_K trajectories (K = 0..2000) ---", flush=True)
    for x in x_values:
        t0 = time.time()
        trajectories[x] = pi_trajectory(x, zeros, dps=DPS)
        print(f"  trajectory for x = {x:.0e}: {len(trajectories[x])} pts  "
              f"[{time.time() - t0:.1f}s]", flush=True)

    # ------------------------------------------------------------------
    # (i) x = 10^4 persistent-error investigation
    # ------------------------------------------------------------------
    print("\n--- (i) x=10^4 persistent-error investigation ---", flush=True)
    x = 10 ** 4
    traj = trajectories[x]
    for K in [100, 200, 500, 1000, 1500, 2000]:
        if K < len(traj):
            raw = traj[K]
            print(f"  K={K:>4}:  raw = {raw:.6f}  rounded = {round(raw)}  "
                  f"residual = {raw - exact[x]:+.4f}", flush=True)
    K_min, K_stable = k_min_metrics(traj, exact[x], run_after=50)
    print(f"  K_min  (first round-correct): {K_min}", flush=True)
    print(f"  K_min* (>=50 consecutive correct): {K_stable}", flush=True)
    # Also report final state
    last = traj[-1]
    print(f"  Final (K=2000) raw = {last:.6f}  residual = "
          f"{last - exact[x]:+.4f}", flush=True)

    # ------------------------------------------------------------------
    # (ii) K_min, K_min* across x
    # ------------------------------------------------------------------
    print("\n--- (ii) K_min and K_min* across x ---", flush=True)
    print(f"{'x':>10} | {'K_min':>6} | {'K_min*':>7} | "
          f"{'sqrt(x)':>10} | {'sqrt(x)/log(x)':>16} | "
          f"{'sqrt(x)*log(x)':>16} | {'K*/sqrt(x)':>11} | "
          f"{'K*/(sqrt(x)/log(x))':>22}", flush=True)
    print("-" * 130, flush=True)
    rows = []
    for x in x_values:
        traj = trajectories[x]
        K_min, K_stable = k_min_metrics(traj, exact[x], run_after=50)
        sx = math.sqrt(x)
        sxlogx = sx / math.log(x)
        sxLogx = sx * math.log(x)
        rows.append((x, K_min, K_stable, sx, sxlogx, sxLogx))
        kf = str(K_min) if K_min is not None else "N/A"
        ks = str(K_stable) if K_stable is not None else "N/A"
        if K_stable is not None:
            r1 = f"{K_stable / sx:.4f}"
            r2 = f"{K_stable / sxlogx:.4f}"
        else:
            r1 = r2 = "--"
        print(f"{x:>10.0e} | {kf:>6} | {ks:>7} | {sx:>10.1f} | "
              f"{sxlogx:>16.2f} | {sxLogx:>16.2f} | "
              f"{r1:>11} | {r2:>22}", flush=True)

    # Linear regression K_min* ~ x^slope.
    pts = [(x, ks) for (x, _, ks, _, _, _) in rows if ks is not None]
    if len(pts) >= 2:
        lx = [math.log(p[0]) for p in pts]
        ly = [math.log(p[1]) for p in pts]
        n = len(pts)
        mx = sum(lx) / n
        my = sum(ly) / n
        num = sum((a - mx) * (b - my) for a, b in zip(lx, ly))
        den = sum((a - mx) ** 2 for a in lx)
        slope = num / den if den else float("nan")
        print(f"\n  Fitted: K_min* ~ x^{slope:.4f}  "
              f"(RH-conjectural ~ x^0.5 with log factors)",
              flush=True)

    # ------------------------------------------------------------------
    # (iii) x = 10^7 residual at K = 2000
    # ------------------------------------------------------------------
    print("\n--- (iii) Residual at x = 10^7 with K up to 2000 ---", flush=True)
    traj7 = trajectories[10 ** 7]
    for K in [50, 200, 500, 1000, 1500, 2000]:
        if K < len(traj7):
            raw = traj7[K]
            print(f"  K={K:>4}:  raw = {raw:.4f}  residual = "
                  f"{raw - exact[10 ** 7]:+.2f}", flush=True)
    T_min_classical = 2 * math.sqrt(10 ** 7) * math.log(10 ** 7) ** 2
    print(f"  Classical T_min for error<0.5 at x=10^7: ~{T_min_classical:.0f} "
          f"(need this many zeros, vs available 2000).", flush=True)

    # ------------------------------------------------------------------
    # Summary
    # ------------------------------------------------------------------
    print("\n=== SUMMARY ===", flush=True)
    for x in x_values:
        traj = trajectories[x]
        K_min, K_stable = k_min_metrics(traj, exact[x], run_after=50)
        final_err = traj[-1] - exact[x]
        if K_stable is not None:
            print(f"  x={x:.0e}: K_min* = {K_stable}, "
                  f"final residual at K=2000 = {final_err:+.3f}", flush=True)
        else:
            print(f"  x={x:.0e}: K_min* not reached within 2000 zeros, "
                  f"final residual at K=2000 = {final_err:+.3f}", flush=True)


if __name__ == "__main__":
    main()
