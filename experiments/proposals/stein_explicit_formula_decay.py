"""
F5 test — Sub-Gaussian Explicit-Formula Error (SGEFE) conjecture.

Compute pi(x) approximation
    f_T(x) = R(x) - sum_{|gamma_n| < T} 2*Re( R(x^{1/2+i*gamma_n}) )
where R is the Riemann prime-counting function, R(x) = sum_{k>=1} mu(k)/k * li(x^{1/k}).
Then plot log|error|(x,T) versus T, T^2/log(x), sqrt(T), and 1/T to determine
the actual decay rate of the truncation error.

SGEFE predicts: log|err| ~ -T^2 / (c log x) (sub-Gaussian)
Classical / RH: log|err| ~ -log T (power-law sqrt(x)/T)

If the slope vs T^2/log(x) is not significantly negative compared to the slope
vs log T, SGEFE is empirically falsified.
"""

import numpy as np
import mpmath as mp
from sympy import primerange

mp.mp.dps = 30

DATA = "/apps/aplikacijos/prime-research/data/zeta_zeros_1000.txt"


def riemann_R_from_logx(log_x, terms=20):
    """R(y) = sum_{k>=1} mu(k)/k * li(y^{1/k}) = sum_k mu(k)/k * Ei(log(y)/k).

    Pass log(y) directly to avoid branch-cut issues when y = x^rho is complex
    with large imaginary part — the correct analytic continuation uses the
    NON-principal log along ds Schwartz path, which equals rho*log(x) directly.
    """
    from sympy import mobius
    s = mp.mpc(0)
    for k in range(1, terms + 1):
        mu = mobius(k)
        if mu == 0:
            continue
        s += mu * mp.ei(log_x / k) / k
    return s


def riemann_R(x, terms=20):
    return riemann_R_from_logx(mp.log(mp.mpc(x)), terms=terms)


def pi_explicit_truncated(x, gammas, T):
    """f_T(x) = R(x) - sum_{|gamma| < T} 2 Re R_log(rho * log x).
    Uses log-space evaluation to handle complex rho correctly.
    """
    log_x = mp.log(mp.mpf(x))
    main = riemann_R_from_logx(log_x).real
    osc = mp.mpf(0)
    for g in gammas:
        if g >= T:
            break
        rho = mp.mpc(0.5, g)
        osc += 2 * riemann_R_from_logx(rho * log_x).real
    return float(main - osc)


def main():
    # Load zeros
    gammas = []
    with open(DATA) as f:
        for line in f:
            line = line.strip()
            if line:
                gammas.append(mp.mpf(line))

    # Targets
    xs = [100, 500, 1000, 5000, 10000]
    Ts = [10, 30, 50, 100, 200, 300, 500, 700, 1000]

    # True pi(x)
    pi_true = {}
    for x in xs:
        pi_true[x] = sum(1 for _ in primerange(2, x + 1))

    # Compute approximations
    print(f"{'x':>8} {'T':>6} {'pi(x)':>8} {'f_T(x)':>14} {'|err|':>14} {'log10|err|':>10}")
    rows = []
    for x in xs:
        for T in Ts:
            f = pi_explicit_truncated(x, gammas, T)
            err = pi_true[x] - f
            ae = abs(err)
            le = float(mp.log10(ae)) if ae > 0 else -99
            print(f"{x:>8} {T:>6} {pi_true[x]:>8} {f:>14.4f} {ae:>14.6f} {le:>10.4f}")
            rows.append((x, T, pi_true[x], f, err))

    # Fit log|err| vs T, log T, T^2/log x, sqrt(T) for each x
    print()
    print("Decay-rate fits per x (log10|err| vs predictor):")
    print(f"{'x':>8} {'slope T':>12} {'slope logT':>12} {'slope T^2/lx':>14} {'slope sqrtT':>12}")
    fits = {}
    for x in xs:
        lerrs = []
        Ts_use = []
        for (xi, Ti, _, _, e) in rows:
            if xi == x and abs(e) > 0:
                lerrs.append(float(mp.log10(abs(e))))
                Ts_use.append(Ti)
        Ts_arr = np.array(Ts_use, dtype=float)
        lerr_arr = np.array(lerrs)
        if len(Ts_arr) < 3:
            continue
        # Linear regression slopes
        s_T = np.polyfit(Ts_arr, lerr_arr, 1)[0]
        s_logT = np.polyfit(np.log(Ts_arr), lerr_arr, 1)[0]
        s_T2lx = np.polyfit(Ts_arr ** 2 / np.log(x), lerr_arr, 1)[0]
        s_sqrtT = np.polyfit(np.sqrt(Ts_arr), lerr_arr, 1)[0]
        fits[x] = (s_T, s_logT, s_T2lx, s_sqrtT)
        print(f"{x:>8} {s_T:>12.5g} {s_logT:>12.5g} {s_T2lx:>14.5g} {s_sqrtT:>12.5g}")

    # R^2 of each model
    print()
    print("R^2 of each predictor model:")
    print(f"{'x':>8} {'R^2 T':>10} {'R^2 logT':>10} {'R^2 T^2/lx':>12} {'R^2 sqrtT':>10}")
    for x in xs:
        lerrs = []
        Ts_use = []
        for (xi, Ti, _, _, e) in rows:
            if xi == x and abs(e) > 0:
                lerrs.append(float(mp.log10(abs(e))))
                Ts_use.append(Ti)
        Ts_arr = np.array(Ts_use, dtype=float)
        lerr_arr = np.array(lerrs)
        if len(Ts_arr) < 3:
            continue
        def r2(pred):
            p = np.polyfit(pred, lerr_arr, 1)
            yhat = np.polyval(p, pred)
            ss_res = np.sum((lerr_arr - yhat) ** 2)
            ss_tot = np.sum((lerr_arr - lerr_arr.mean()) ** 2)
            return 1 - ss_res / ss_tot if ss_tot > 0 else 0
        r_T = r2(Ts_arr)
        r_logT = r2(np.log(Ts_arr))
        r_T2lx = r2(Ts_arr ** 2 / np.log(x))
        r_sqrtT = r2(np.sqrt(Ts_arr))
        print(f"{x:>8} {r_T:>10.4f} {r_logT:>10.4f} {r_T2lx:>12.4f} {r_sqrtT:>10.4f}")

    # Bottom line: for SGEFE, log|err| should decrease with T^2/log x; check sign of slope.
    print()
    print("Verdict diagnostic:")
    print("- SGEFE predicts slope_(T^2/lx) << 0 with high R^2.")
    print("- Power-law (RH) predicts slope_logT << 0 with high R^2; slope_T near 0.")
    print("- Random-walk / no decay predicts all slopes near 0 with low R^2.")


if __name__ == "__main__":
    main()
