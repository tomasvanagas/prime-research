"""S184 verify — adversarial test of S183's asymptote-pinning claim.

S183 claimed: 5-point fit `ratio(d) = a + b/d` over d in {14,18,20,22,24}
gives `a = 0.2117`, "within 1% of theoretical 0.21".

Adversarial questions:
1. d=16 (datum 0.2132 from S173, k_*=8 by S74's table) was excluded
   from the fit on the basis that "k_* selection at d=16 is uncertain".
   Does including it shift the asymptote materially?
2. The model `a + b/d` showed systematically positive then negative
   residuals (per S183 §"5-point asymptote fit"). With one extra
   parameter, what asymptote do we get?
3. Bootstrap CI: how stable is the asymptote under leave-one-out?
4. The k_* "linear extrapolation" gave 78 at d=24 — but the underlying
   sequence (5, 15, 26, 45, 78) has growing deltas (10, 11, 19, 33),
   inconsistent with a true linear rule. What asymptote do we get if
   k_* is chosen by a different sensible rule?

If the asymptote is fragile to any of (1)-(4), the "0.2117 within 1% of
0.21" framing is overstated; the honest scope is the bootstrap interval.
"""
import numpy as np


def fit(model, ds, ys, p0=None):
    """Linear least-squares. model(d) -> design row. Returns (params, asymptote, resid)."""
    X = np.array([model(d) for d in ds])
    y = np.array(ys)
    params, *_ = np.linalg.lstsq(X, y, rcond=None)
    fit_y = X @ params
    resid = y - fit_y
    asymptote = params[0]  # coefficient of the constant column
    return params, asymptote, resid


def model_a_b_over_d(d):
    return [1.0, 1.0 / d]


def model_a_b_c_over_d_d2(d):
    return [1.0, 1.0 / d, 1.0 / d**2]


def model_a_b_over_logd(d):
    return [1.0, 1.0 / np.log(d)]


def model_a_b_over_sqrtd(d):
    return [1.0, 1.0 / np.sqrt(d)]


def model_pure_power(d):
    # y = a * d^b -> log(y) = log(a) + b * log(d). But we want asymptote of a + b * d^k.
    # Skip; not in same form.
    raise NotImplementedError


# ----------------------------------------------------------------------
# DATA: spike-block / pi(N) at each d.
# Sources: S82 saved JSONs (d=14,18,20), S173 (d=16, k_*=8), S174 (d=22,
# k_*=45 via linear-log extrap), S183 (d=24, k_*=78 via "linear extrap").
# ----------------------------------------------------------------------
data_5pt_S183 = [
    (14, 0.2236),
    (18, 0.2212),
    (20, 0.2198),
    (22, 0.2178),
    (24, 0.2160),
]

# 6-point: include d=16 from S173 (k_*=8, S74's data-derived value).
data_6pt = [
    (14, 0.2236),
    (16, 0.2132),
    (18, 0.2212),
    (20, 0.2198),
    (22, 0.2178),
    (24, 0.2160),
]


def report(name, data, model, model_name):
    ds = [d for d, _ in data]
    ys = [y for _, y in data]
    params, asy, resid = fit(model, ds, ys)
    print(f"\n=== {name}, model = {model_name} ===")
    print(f"  params: {[f'{p:.6f}' for p in params]}")
    print(f"  asymptote a = {asy:.4f}")
    print(f"  residuals: {[f'{r:+.5f}' for r in resid]}")
    print(f"  max |resid|: {np.max(np.abs(resid)):.5f}")
    return asy


def loocv(data, model, model_name):
    """Leave-one-out asymptotes. Returns list of (omitted_d, asymptote)."""
    out = []
    for i in range(len(data)):
        sub = data[:i] + data[i+1:]
        ds = [d for d, _ in sub]
        ys = [y for _, y in sub]
        params, asy, _ = fit(model, ds, ys)
        out.append((data[i][0], asy))
    return out


def main():
    print("S184 — adversarial robustness probe of S183's `a = 0.2117` asymptote claim")
    print("=" * 80)

    # Probe A: reproduce S183's headline (5pt, a + b/d).
    a_5pt_lin = report("(A) 5pt S183 trajectory", data_5pt_S183,
                        model_a_b_over_d, "a + b/d")

    # Probe B: include d=16.
    a_6pt_lin = report("(B) 6pt INCLUDE d=16", data_6pt,
                        model_a_b_over_d, "a + b/d")

    # Probe C: alternative model (a + b/d + c/d^2) on 5pt.
    a_5pt_quad = report("(C) 5pt, quadratic-in-1/d", data_5pt_S183,
                         model_a_b_c_over_d_d2, "a + b/d + c/d^2")

    # Probe D: alternative model on 6pt.
    a_6pt_quad = report("(D) 6pt, quadratic-in-1/d", data_6pt,
                         model_a_b_c_over_d_d2, "a + b/d + c/d^2")

    # Probe E: 1/log(d) model (slower convergence, more honest for log-scale data).
    a_5pt_log = report("(E) 5pt, 1/log(d)", data_5pt_S183,
                        model_a_b_over_logd, "a + b/log(d)")
    a_6pt_log = report("(F) 6pt, 1/log(d)", data_6pt,
                        model_a_b_over_logd, "a + b/log(d)")

    # Probe G: 1/sqrt(d) model.
    a_5pt_sqrt = report("(G) 5pt, 1/sqrt(d)", data_5pt_S183,
                         model_a_b_over_sqrtd, "a + b/sqrt(d)")

    # Probe H: leave-one-out on 5pt with a + b/d.
    print("\n=== (H) 5pt LOO with model a + b/d ===")
    for d_out, asy in loocv(data_5pt_S183, model_a_b_over_d, "a + b/d"):
        print(f"  omit d={d_out}: asymptote = {asy:.4f}")

    # Probe I: LOO on 6pt.
    print("\n=== (I) 6pt LOO with model a + b/d ===")
    loo_6 = loocv(data_6pt, model_a_b_over_d, "a + b/d")
    for d_out, asy in loo_6:
        print(f"  omit d={d_out}: asymptote = {asy:.4f}")
    asys_loo_6 = [a for _, a in loo_6]
    print(f"  range: [{min(asys_loo_6):.4f}, {max(asys_loo_6):.4f}]")

    # Summary table.
    print("\n" + "=" * 80)
    print("SUMMARY: asymptote `a` under different model / data choices")
    print("=" * 80)
    print(f"{'(A) 5pt, a+b/d (S183 hdline)':<45} a = {a_5pt_lin:.4f}")
    print(f"{'(B) 6pt, a+b/d':<45} a = {a_6pt_lin:.4f}")
    print(f"{'(C) 5pt, a+b/d+c/d²':<45} a = {a_5pt_quad:.4f}")
    print(f"{'(D) 6pt, a+b/d+c/d²':<45} a = {a_6pt_quad:.4f}")
    print(f"{'(E) 5pt, a+b/log(d)':<45} a = {a_5pt_log:.4f}")
    print(f"{'(F) 6pt, a+b/log(d)':<45} a = {a_6pt_log:.4f}")
    print(f"{'(G) 5pt, a+b/sqrt(d)':<45} a = {a_5pt_sqrt:.4f}")
    print()
    asys = [a_5pt_lin, a_6pt_lin, a_5pt_quad, a_6pt_quad,
            a_5pt_log, a_6pt_log, a_5pt_sqrt]
    print(f"Range across model+data choices: [{min(asys):.4f}, {max(asys):.4f}]")
    print(f"Spread: {max(asys) - min(asys):.4f}")
    print(f"Theory target: 0.2100")
    print(f"S183 headline: a = 0.2117 (within 1% of 0.21)")

    # Adversarial verdict criterion.
    if max(asys) - min(asys) > 0.05:
        print("\nADVERSARIAL VERDICT: model-choice spread > 0.05 — '1% of 0.21' overstated.")
    elif max(asys) - min(asys) > 0.02:
        print("\nADVERSARIAL VERDICT: model-choice spread 0.02-0.05 — 'within 1%' overstated, asymptote ≈ 0.21 ± 0.02.")
    else:
        print("\nADVERSARIAL VERDICT: model-choice spread < 0.02 — '1% of 0.21' framing robust.")


if __name__ == "__main__":
    main()
