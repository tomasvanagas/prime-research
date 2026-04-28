#!/usr/bin/env python3
"""S183 — direct SVD of chi_P at d=24 to test S168/S169 21% prediction.

Background: S169 verify chain (S170-S182) has CONFIRMed/PARTIAL'd the 21%
spike-block claim 12 times. S178, S179 identified d=24 SVD as the most
informative remaining probe — never run because S82's saved JSONs cover
only d ∈ {14, 18, 20}, S173 added d=16, S174 added d=22. d=24 is the
largest scale the analytic side reaches (S169's table runs through d=24).

This script:
  1. Builds chi_P on [0, 2^24).
  2. Reshapes to 4096 x 4096.
  3. Computes the full SVD (~minutes for 4096x4096 dense).
  4. Computes the spike-block sum at the linear-extrapolation k_*(24)=78
     from S174's fit log k_* = 0.275 d - 2.24, plus a sweep over k_*.
  5. Compares to S169's predicted analytic cum(Q*=33) at d=24 (= 264244).
  6. Computes Q_eff (smallest sqf q with cum(q) >= spike_block).

Pre-stated falsifiers:
  PR1. spike_block / pi(N) within 5% of S169 trajectory's 0.218..0.224
       band (d=14..22).  PASS if 0.20 <= spike/pi <= 0.24.
  PR2. log Q_eff / log N within range [0.16, 0.20] (the d=14..22 band).
       PASS if 0.15 <= exp <= 0.20.
  PR3. Asymptote pinning: 5-point fit (d=14..24) gives intercept-at-
       d=infinity within [0.18, 0.24].  PASS if 0.18 <= a <= 0.24.

REFUTE conditions:
  - spike/pi(N) at d=24 lies OUTSIDE [0.18, 0.25] (would refute the 21% trend).
  - 5-point fit asymptote outside [0.17, 0.25] (would refute 0.21 pinning).
  - Q_eff exponent <0.10 or >0.25 (would refute the matching cutoff claim).
"""

import json
import math
import time
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent


def sieve(N):
    is_p = np.ones(N + 1, dtype=bool)
    is_p[:2] = False
    for i in range(2, int(N**0.5) + 1):
        if is_p[i]:
            is_p[i * i :: i] = False
    return is_p[:N]  # length-N indicator on [0, N)


def factor(q):
    f = []
    n = q
    p = 2
    while p * p <= n:
        if n % p == 0:
            e = 0
            while n % p == 0:
                n //= p
                e += 1
            f.append((p, e))
        p += 1
    if n > 1:
        f.append((n, 1))
    return f


def is_squarefree(q):
    return q == 1 or all(e == 1 for _, e in factor(q))


def euler_phi(q):
    if q == 1:
        return 1
    res = q
    for p, _ in factor(q):
        res = res // p * (p - 1)
    return res


def fourier_energy_at_q(prime_arr, q, N):
    """Same definition as S169."""
    energy = 0.0
    for k in range(1, q):
        if math.gcd(k, q) != 1:
            continue
        angle = 2 * np.pi * k * prime_arr / q
        c = float(np.sum(np.cos(angle)))
        s = float(np.sum(np.sin(angle)))
        energy += (c * c + s * s) / N
    return energy


def main():
    out = {}
    d = 24
    N = 1 << d  # 16,777,216
    side = 1 << (d // 2)  # 4096

    print(f"=== S183 d=24 SVD verification ===")
    print(f"N = {N}, matrix shape = {side} x {side}")

    print("Step 1: sieve chi_P ...")
    t0 = time.time()
    chi_P = sieve(N).astype(np.float64)
    pi_N = int(chi_P.sum())
    print(f"  pi(N) = {pi_N}, time = {time.time() - t0:.1f}s")
    out["N"] = N
    out["pi_N"] = pi_N
    out["d"] = d
    out["shape"] = [side, side]

    print("Step 2: reshape and SVD ...")
    t0 = time.time()
    M = chi_P.reshape(side, side)
    # Use full SVD: sigma values only (no U, V).  Faster: compute_uv=False.
    sigmas = np.linalg.svd(M, compute_uv=False)
    print(f"  SVD time = {time.time() - t0:.1f}s")
    print(f"  top 30 sigmas:", " ".join(f"{s:.2f}" for s in sigmas[:30]))
    print(f"  sigma_0^2 = {sigmas[0]**2:.2f}")
    print(f"  Frobenius^2 = sum sigma^2 = {(sigmas**2).sum():.2f}, pi(N) = {pi_N}")

    # Save sigmas
    out["sigmas_top100"] = [float(s) for s in sigmas[:100]]
    out["sigma_0_sq"] = float(sigmas[0] ** 2)
    out["frobenius_sq"] = float((sigmas**2).sum())

    # Step 3: spike-block sweep over k_*
    print("Step 3: spike-block sweep ...")
    target_21 = 0.21 * pi_N
    sigma_sq = sigmas ** 2
    k_star_sweep = list(range(20, 121, 5)) + [78]  # include linear-extrap prediction
    k_star_sweep.sort()

    sweep = {}
    for k_star in k_star_sweep:
        # Spike block = sigmas[1..k_star]^2 (excluding rank-1 mean)
        spike_block = float(sigma_sq[1 : 1 + k_star].sum())
        ratio_to_21pi = spike_block / target_21
        ratio_to_pi = spike_block / pi_N
        sweep[k_star] = {
            "k_star": k_star,
            "spike_block_sum": spike_block,
            "spike_over_pi_N": ratio_to_pi,
            "ratio_to_21pi": ratio_to_21pi,
        }
    out["sweep"] = sweep

    # The headline k_* (linear extrapolation S174)
    k_star_lin = 78
    headline = sweep[k_star_lin]
    out["headline_k_star"] = k_star_lin
    out["headline"] = headline
    print(f"  k_*={k_star_lin} (linear extrap): "
          f"spike/pi = {headline['spike_over_pi_N']:.4f}, "
          f"spike/(0.21 pi) = {headline['ratio_to_21pi']:.4f}")

    # Step 4: matching Q_eff (smallest sqf q with cum(q) >= spike_block).
    print("Step 4: compute analytic cum(q) at d=24 for sqf q in [2, 60] ...")
    t0 = time.time()
    primes = np.where(chi_P > 0)[0]
    print(f"  primes array length = {len(primes)}, time to extract = {time.time() - t0:.1f}s")

    sqf_qs = [q for q in range(2, 61) if is_squarefree(q)]
    cum_table = []
    cum_emp = 0.0
    print(f"  Fourier sieve over {len(sqf_qs)} sqf q ...")
    t0 = time.time()
    for q in sqf_qs:
        E = fourier_energy_at_q(primes, q, N)
        cum_emp += E
        cum_table.append((q, E, cum_emp))
    print(f"  Fourier sieve time = {time.time() - t0:.1f}s")
    out["cum_table"] = [{"q": q, "E": E, "cum": c} for q, E, c in cum_table]

    # Q_eff for headline k_star
    spike_headline = headline["spike_block_sum"]
    Q_eff = None
    for q, E, c in cum_table:
        if c >= spike_headline:
            Q_eff = q
            break
    if Q_eff is None:
        Q_eff_str = ">60"
        exp_eff = float("nan")
    else:
        exp_eff = math.log(Q_eff) / math.log(N)
        Q_eff_str = str(Q_eff)
    print(f"  Q_eff (k_*=78) = {Q_eff_str}, log Q_eff / log N = {exp_eff:.4f}")
    out["Q_eff_headline"] = {"Q_eff": Q_eff, "exp": exp_eff}

    # Q_eff for a sweep of k_*
    Q_eff_table = {}
    for k_star, row in sweep.items():
        sb = row["spike_block_sum"]
        Q = None
        for q, E, c in cum_table:
            if c >= sb:
                Q = q
                break
        Q_eff_table[k_star] = {
            "Q_eff": Q,
            "exp": math.log(Q) / math.log(N) if Q else None,
        }
    out["Q_eff_sweep"] = Q_eff_table

    # Step 5: 5-point fit on spike_over_pi_N vs 1/d using d=14, 18, 20, 22, 24
    # Prior reported values:
    prior = {
        14: 0.2236,  # S169
        16: 0.2132,  # S173 (k_*=8)
        18: 0.2212,  # S169
        20: 0.2198,  # S169
        22: 0.2178,  # S174 (k_*=45 lin extrap)
    }
    full = dict(prior)
    full[24] = headline["spike_over_pi_N"]
    out["all_d_ratios"] = full

    # Linear fit ratio = a + b/d
    ds = sorted(full.keys())
    xs = np.array([1.0 / dd for dd in ds])
    ys = np.array([full[dd] for dd in ds])
    A = np.vstack([np.ones_like(xs), xs]).T
    coef, *_ = np.linalg.lstsq(A, ys, rcond=None)
    a, b = float(coef[0]), float(coef[1])
    out["fit"] = {"a": a, "b": b, "d_used": ds, "ratios": [float(full[dd]) for dd in ds]}
    print(f"  5-point fit: ratio = {a:.4f} + {b:.4f}/d "
          f"(asymptote a = {a:.4f}, slope b = {b:.4f})")

    # Pre-stated falsifiers
    PR1_pass = 0.20 <= headline["spike_over_pi_N"] <= 0.24
    PR2_pass = 0.15 <= exp_eff <= 0.20 if exp_eff == exp_eff else False
    PR3_pass = 0.18 <= a <= 0.24
    out["falsifiers"] = {
        "PR1_spike_fraction_in_band": PR1_pass,
        "PR2_Q_eff_exponent_in_band": PR2_pass,
        "PR3_asymptote_in_band": PR3_pass,
    }
    print(f"  PR1 (spike/pi in [0.20, 0.24]): {'PASS' if PR1_pass else 'FAIL'} (got {headline['spike_over_pi_N']:.4f})")
    print(f"  PR2 (Q_eff exp in [0.15, 0.20]): {'PASS' if PR2_pass else 'FAIL'} (got {exp_eff:.4f})")
    print(f"  PR3 (asymptote in [0.18, 0.24]): {'PASS' if PR3_pass else 'FAIL'} (got {a:.4f})")

    (HERE / "d24_svd_verify_results.json").write_text(json.dumps(out, indent=2))
    print(f"\nSaved: {HERE/'d24_svd_verify_results.json'}")


if __name__ == "__main__":
    main()
