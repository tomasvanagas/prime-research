"""Post-hoc diagnostic: Gamma_1-cosine ceiling for cross-modulus lag-1 autocorr.

For each m, the leading explicit-formula contribution to r_R(k) =
pi(m^k) - R(m^k) is roughly proportional to (1/gamma_1) *
cos(gamma_1 * k * log(m) + phase_1). Higher zeros add Weyl-equidistributed
phases that average out as k grows.

The lag-1 autocorrelation of a pure cosine sequence (cos(omega * k))
with omega = gamma_1 * log(m) is exactly cos(omega) (mod 2*pi).

This script computes:
  phi_1(m) = gamma_1 * log(m) mod 2*pi
  cos(phi_1(m))                  -- ceiling for lag-1 autocorr from gamma_1 mode

and compares to the empirical max-lag-1 autocorr from
cross_modulus_pi_structure.summary.json.
"""
from __future__ import annotations

import json
from pathlib import Path

import mpmath

mpmath.mp.dps = 60

# First few non-trivial zeta zeros (imaginary parts).
GAMMAS = [
    14.134725141734693790457251983562470270784257115699,
    21.022039638771554992628479593896902777334340524903,
    25.010857580145688763213790992562821818659549672558,
    30.424876125859513210311897530584091320181560023715,
    32.935061587739189690662368964074903488812715603517,
]


def main():
    here = Path(__file__).parent
    summary = json.loads((here / "summary.json").read_text())
    per_m = summary["per_m"]

    print(f"{'m':>4} {'phi_1':>10} {'cos(phi_1)':>12} "
          f"{'emp_lag1':>10} {'max_lag':>8} {'max|ac|':>10}")
    print("-" * 64)
    g1 = mpmath.mpf(GAMMAS[0])
    rows = []
    for m_str, info in per_m.items():
        m = int(m_str)
        phi_1 = (g1 * mpmath.log(m)) % (2 * mpmath.pi)
        cphi = mpmath.cos(phi_1)
        # Recover empirical lag-1 autocorr from raw_data
        raw = json.loads((here / "raw_data.json").read_text())
        ac_R = raw[m_str]["ac_R"]
        emp_lag1 = float(ac_R["1"])
        max_lag = info["max_abs_ac_lag"]
        max_ac = info["max_abs_ac"]
        rows.append({
            "m": m, "phi_1": float(phi_1), "cos_phi_1": float(cphi),
            "emp_lag1": emp_lag1, "max_lag": max_lag, "max_abs_ac": max_ac,
        })
        print(f"{m:>4} {float(phi_1):>10.4f} {float(cphi):>+12.4f} "
              f"{emp_lag1:>+10.4f} {max_lag:>8d} {max_ac:>+10.4f}")

    # Sign-agreement test: does sign(cos(phi_1)) match sign(emp_lag1)?
    sign_match = sum(1 for r in rows
                    if (r["cos_phi_1"] > 0) == (r["emp_lag1"] > 0))
    print(f"\nsign(cos(phi_1)) matches sign(emp_lag1) for {sign_match}/{len(rows)} m")

    # Magnitude bound: |emp_lag1| <= |cos(phi_1)|?
    bound_holds = sum(1 for r in rows
                     if abs(r["emp_lag1"]) <= abs(r["cos_phi_1"]) + 0.05)
    print(f"|emp_lag1| <= |cos(phi_1)| + 0.05 for {bound_holds}/{len(rows)} m")

    out = {"rows": rows, "sign_match": sign_match,
           "bound_holds": bound_holds, "n": len(rows)}
    (here / "gamma1_phase_diagnostic.json").write_text(json.dumps(out, indent=2))
    print(f"\nwrote {here / 'gamma1_phase_diagnostic.json'}")


if __name__ == "__main__":
    main()
