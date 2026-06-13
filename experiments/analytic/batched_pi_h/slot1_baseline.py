"""
slot1_baseline.py — Thread 8 (P2) commit slot 1
================================================

Establish the empirical baseline for batched-on-h prime-gap counting:

    pi_h(x) := #{ p prime, p <= x : p + h prime }

for h in {2, 4, 6, ..., H} simultaneously.

Three measurements:

  (M1) BATCHED SIEVE: sieve [1, x+H] once, then for each prime p <= x and
       each h in {2, 4, ..., H} test is_prime(p+h). Total cost
       O((x+H) log log x) sieve + O(pi(x) * (H/2)) pair-tests.

  (M2) PER-H SIEVE: for each h, sieve [1, x+h] independently and count
       pairs. Total cost M * O((x+h) log log x).
       (We measure the single-h cost and multiply, so we can predict M*T_1.)

  (M3) HARDY-LITTLEWOOD APPROXIMATION:
        HL_h(x) = S_h * li_2(x),   li_2(x) = integral_2^x dt/log^2 t,
        S_h = 2 C_2 prod_{p|h, p odd} (p-1)/(p-2),  C_2 = 0.66016181...
       Per-h cost: factor h (polylog x), multiply by precomputed li_2(x).

Outputs:
  - exact pi_h(x) for h in {2, 4, ..., H} at x in {x_values}
  - HL_h(x) and signed error  pi_h(x) - HL_h(x)
  - timing tables for M1/M2/M3
  - per-h amortised cost shape

Falsification statement:
  If the H/log measurements show batched-sieve per-h cost stays at
  Theta(x/log x) (= one prime-density factor) regardless of M, then
  EXACT pi_h(x) admits NO polylog amortisation and matches the P1
  shape (Thread 6). The partial-positive then must come from the
  approximate (HL) regime: HL_h(x) costs O(polylog x) per-h with
  empirical error bounded by the Hardy-Littlewood deviation.

Edges cited:
  - E1.5 (information density of pi(x))
  - E6.x (sieve hierarchy)
  - S224 / Thread 5 Correlation Dichotomy (the dichotomy pattern this
    slot probes for the h-axis instead of the x-axis)
  - S231 (Thread 6 P1 final wrap; demonstrates the negative shape we
    test against)
"""

from __future__ import annotations

import math
import sys
import time
from pathlib import Path
from typing import Callable

# ----------------------------------------------------------------------
# Constants
# ----------------------------------------------------------------------

# Hardy-Littlewood twin-prime constant
# C_2 = prod_{p >= 3} (1 - 1/(p-1)^2) = 0.660161815846869573927812110...
C_2 = 0.6601618158468695739278121100145557784326233628

# ----------------------------------------------------------------------
# Sieve primitives
# ----------------------------------------------------------------------


def sieve_eratosthenes(N: int) -> bytearray:
    """Return bytearray is_p of length N+1, is_p[i] = 1 iff i is prime."""
    is_p = bytearray(b"\x01") * (N + 1)
    is_p[0] = 0
    is_p[1] = 0
    bound = int(N**0.5) + 1
    for i in range(2, bound):
        if is_p[i]:
            start = i * i
            step = i
            is_p[start : N + 1 : step] = b"\x00" * (((N - start) // step) + 1)
    return is_p


def count_pi(is_p: bytearray, x: int) -> int:
    """Count primes <= x using is_p."""
    return sum(is_p[2 : x + 1])


# ----------------------------------------------------------------------
# Exact pi_h(x) — batched and single-h
# ----------------------------------------------------------------------


def pi_h_batched(x: int, H: int, is_p: bytearray | None = None) -> dict[int, int]:
    """
    Compute pi_h(x) for h in {2, 4, ..., H} simultaneously.
    Single sieve pass for [1, x+H], then for each prime p <= x do H/2 tests.
    """
    assert H % 2 == 0
    if is_p is None:
        is_p = sieve_eratosthenes(x + H)
    counts = {h: 0 for h in range(2, H + 1, 2)}
    # Iterate over primes p <= x via the sieve mask
    # For speed: extract primes once, then loop
    primes = [i for i in range(3, x + 1) if is_p[i]]
    # NOTE: we exclude p=2 because p+h with h even => p+h even; the only
    # prime pair (2, 2+h) requires 2+h prime. Add separately:
    if x >= 2:
        for h in range(2, H + 1, 2):
            if is_p[2 + h]:
                counts[h] += 1
    for p in primes:
        for h in range(2, H + 1, 2):
            if is_p[p + h]:
                counts[h] += 1
    return counts


def pi_h_single(x: int, h: int) -> int:
    """Single-h count via independent sieve of [1, x+h]."""
    is_p = sieve_eratosthenes(x + h)
    c = 0
    for p in range(2, x + 1):
        if is_p[p] and is_p[p + h]:
            c += 1
    return c


# ----------------------------------------------------------------------
# Hardy-Littlewood approximation
# ----------------------------------------------------------------------


def li2(x: float) -> float:
    """
    li_2(x) = integral_2^x dt / log^2(t).
    Use scipy if available, else trapezoidal rule with adaptive N.
    """
    try:
        from scipy.integrate import quad

        v, _ = quad(lambda t: 1.0 / math.log(t) ** 2, 2.0, float(x), limit=200)
        return v
    except Exception:
        # Composite Simpson's rule
        N = max(10000, int(math.log(x) ** 3))
        if N % 2 == 1:
            N += 1
        a, b = 2.0, float(x)
        h = (b - a) / N
        s = 1.0 / math.log(a) ** 2 + 1.0 / math.log(b) ** 2
        for k in range(1, N, 2):
            s += 4.0 / math.log(a + k * h) ** 2
        for k in range(2, N, 2):
            s += 2.0 / math.log(a + k * h) ** 2
        return s * h / 3.0


def singular_series(h: int) -> float:
    """
    S_h = 2 * C_2 * prod_{p|h, p odd} (p - 1)/(p - 2).
    Polylog in h: trial-divide h to find odd prime factors.
    Returns 0 for odd h (twin-prime-style admissible 2-tuple {0, h}
    requires h even, else h+0 covers a residue class mod 2).
    """
    if h % 2 != 0:
        return 0.0
    n = h
    while n % 2 == 0:
        n //= 2
    factor = 1.0
    p = 3
    while p * p <= n:
        if n % p == 0:
            factor *= (p - 1) / (p - 2)
            while n % p == 0:
                n //= p
        p += 2
    if n > 1:
        factor *= (n - 1) / (n - 2)
    return 2.0 * C_2 * factor


# ----------------------------------------------------------------------
# Experiment driver
# ----------------------------------------------------------------------


def run_anchor(x: int, H: int, do_per_h_timing: bool = False) -> dict:
    """
    Run all three measurements at one anchor x.
    Returns a dict suitable for serialising into a markdown table.
    """
    print(f"\n=== ANCHOR x = {x:,}  H = {H} ===", flush=True)

    # M1: BATCHED SIEVE (sieve once, count all h)
    t0 = time.perf_counter()
    is_p = sieve_eratosthenes(x + H)
    t_sieve = time.perf_counter() - t0
    pi_x = count_pi(is_p, x)
    print(f"  pi(x) = {pi_x:,};  sieve cost t_sieve = {t_sieve:.3f}s")

    t0 = time.perf_counter()
    exact = pi_h_batched(x, H, is_p=is_p)
    t_batched_count = time.perf_counter() - t0
    t_batched_total = t_sieve + t_batched_count

    M = H // 2  # number of h-values
    per_h_amort_batched = t_batched_total / M
    print(
        f"  M1 batched: t_sieve+count = {t_batched_total:.3f}s "
        f"({t_sieve:.3f} sieve + {t_batched_count:.3f} count over {M} h-values), "
        f"per-h amortised = {per_h_amort_batched:.4f}s"
    )

    # M2: PER-H SIEVE — measure single-h, predict M*T_1
    if do_per_h_timing:
        t0 = time.perf_counter()
        single_count = pi_h_single(x, 2)
        t_single = time.perf_counter() - t0
        predicted_per_h_total = M * t_single
        print(
            f"  M2 per-h single (h=2): t = {t_single:.3f}s, "
            f"predicted M*T_1 = {predicted_per_h_total:.3f}s"
        )
        speedup = predicted_per_h_total / t_batched_total if t_batched_total > 0 else float("inf")
        # sanity check
        assert single_count == exact[2], (single_count, exact[2])
    else:
        t_single = None
        predicted_per_h_total = None
        speedup = None

    # M3: HARDY-LITTLEWOOD APPROXIMATION
    t0 = time.perf_counter()
    li2_x = li2(x)
    t_li2 = time.perf_counter() - t0
    t0 = time.perf_counter()
    hl = {h: singular_series(h) * li2_x for h in exact}
    t_HL_all = time.perf_counter() - t0
    print(
        f"  M3 HL approx: li_2(x) = {li2_x:,.2f} (t={t_li2:.4f}s), "
        f"all S_h*li_2 over {M} h-values t = {t_HL_all:.6f}s, "
        f"per-h = {t_HL_all/M*1e6:.2f} µs"
    )

    # Residuals
    residuals = {h: exact[h] - hl[h] for h in exact}
    abs_residuals = [abs(r) for r in residuals.values()]
    sqrt_x = math.sqrt(x)
    mean_abs_norm = sum(abs_residuals) / M / sqrt_x
    max_abs_norm = max(abs_residuals) / sqrt_x
    print(
        f"  HL residual: mean|err|/sqrt(x) = {mean_abs_norm:.4f}, "
        f"max|err|/sqrt(x) = {max_abs_norm:.4f}"
    )

    # Sample h values for the table
    sample_h = [2, 4, 6, 8, 10, 12, 30, 60, 90, 120, 150, 180, 198]
    sample_h = [h for h in sample_h if h <= H]
    sample_table = []
    for h in sample_h:
        e = exact[h]
        a = hl[h]
        r = e - a
        sample_table.append((h, e, a, r, r / sqrt_x))

    return {
        "x": x,
        "H": H,
        "M": M,
        "pi_x": pi_x,
        "t_sieve": t_sieve,
        "t_batched_count": t_batched_count,
        "t_batched_total": t_batched_total,
        "per_h_amort_batched": per_h_amort_batched,
        "t_single_h2": t_single,
        "predicted_per_h_total": predicted_per_h_total,
        "speedup_batched_vs_perh": speedup,
        "li2_x": li2_x,
        "t_li2": t_li2,
        "t_HL_all": t_HL_all,
        "per_h_HL": t_HL_all / M,
        "exact": exact,
        "hl": hl,
        "residuals": residuals,
        "mean_abs_residual_norm": mean_abs_norm,
        "max_abs_residual_norm": max_abs_norm,
        "sample_table": sample_table,
    }


def main():
    # Anchors and H values
    # We pick H = 2 * floor(log^2 x / 2) so M = floor(log^2 x).
    # That keeps the polylog regime honest.
    plan = []
    for x in [10**5, 10**6, 10**7]:
        log2x = max(2, int(math.log(x) ** 2))
        # round up to nearest even, cap at 200
        H = min(200, 2 * ((log2x + 1) // 2))
        plan.append((x, H))
    # also a fixed-H run at x=10^6 with several H values
    fixed_x = 10**6
    H_sweep = [20, 50, 100, 200, 400]

    all_results = []

    print("\n========== ANCHORS WITH H = round(log^2 x) ==========")
    for x, H in plan:
        # Per-h single-h timing only at smaller x to avoid quadratic blowup
        do_per_h = x <= 10**6
        r = run_anchor(x, H, do_per_h_timing=do_per_h)
        all_results.append(("plan", r))

    print("\n========== H-SWEEP AT x = 10^6 ==========")
    for H in H_sweep:
        r = run_anchor(fixed_x, H, do_per_h_timing=False)
        all_results.append(("Hsweep", r))

    # Save raw results to a TSV-ish text dump for later inspection
    out_dir = Path(__file__).parent
    raw = out_dir / "slot1_raw.tsv"
    with raw.open("w") as f:
        f.write("# tag\tx\tH\tM\tpi_x\tt_sieve\tt_batched_count\tt_batched_total\t"
                "per_h_amort_batched\tt_single_h2\tspeedup\tt_li2\tt_HL_all\tper_h_HL\t"
                "mean_abs_resid_over_sqrt_x\tmax_abs_resid_over_sqrt_x\n")
        for tag, r in all_results:
            f.write(
                "\t".join(
                    str(v)
                    for v in [
                        tag,
                        r["x"],
                        r["H"],
                        r["M"],
                        r["pi_x"],
                        f'{r["t_sieve"]:.4f}',
                        f'{r["t_batched_count"]:.4f}',
                        f'{r["t_batched_total"]:.4f}',
                        f'{r["per_h_amort_batched"]:.6f}',
                        f'{r["t_single_h2"]:.4f}' if r["t_single_h2"] is not None else "NA",
                        f'{r["speedup_batched_vs_perh"]:.3f}'
                        if r["speedup_batched_vs_perh"] is not None
                        else "NA",
                        f'{r["t_li2"]:.6f}',
                        f'{r["t_HL_all"]:.6f}',
                        f'{r["per_h_HL"]:.8f}',
                        f'{r["mean_abs_residual_norm"]:.6f}',
                        f'{r["max_abs_residual_norm"]:.6f}',
                    ]
                )
                + "\n"
            )
    print(f"\nRaw results written to {raw}")

    # Also dump the sample tables
    samp = out_dir / "slot1_samples.tsv"
    with samp.open("w") as f:
        f.write("# tag\tx\tH\th\texact_pi_h\tHL\tresidual\tresidual_over_sqrt_x\n")
        for tag, r in all_results:
            for h, e, a, res, resn in r["sample_table"]:
                f.write(
                    f"{tag}\t{r['x']}\t{r['H']}\t{h}\t{e}\t{a:.4f}\t{res:.4f}\t{resn:.4f}\n"
                )
    print(f"Sample tables written to {samp}")
    return all_results


if __name__ == "__main__":
    main()
