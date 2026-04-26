"""
Semiclassical / periodic-orbit GROUPED resummation of the explicit
formula for pi(x), with rounding to integer.

Idea: the explicit formula is
    pi(x) = R(x) - sum_rho R(x^rho) - <small terms>

where rho runs over ALL non-trivial zeta zeros. Direct evaluation requires
~x^{1/2} / log(x) zeros to bring error below 1, which kills polylog hope.

Standard approaches truncate by individual zeros. We try a different
truncation: GROUP zeros into "semiclassical orbits" indexed by the
period of an underlying putative chaotic dynamical system and
resum each group analytically.

In the Hilbert-Polya framework, zeros are eigenvalues of a yet-unknown
operator. In the Berry-Keating model, the zeros are quantized energies
of the 1D Hamiltonian xp. The PERIODIC ORBITS of the corresponding
chaotic flow (geodesics on a hyperbolic surface, etc.) give CLOSED-FORM
contributions to the spectral sum after Gutzwiller's trace formula.

We test a SIMPLER question: if zeros are grouped by their ordinate
into BLOCKS of size B, can we approximate the sum over each block
analytically (using the average density 1/(2*pi)*log(t/(2*pi))), and
sum only O(log x) such blocks?

Concretely:
  sum_{rho} li(x^rho) = sum_{blocks i} avg-contribution(block i)
where the average over a block is well approximated by a smooth
integral. The CORRECTION (deviation from average) is what we're
interested in.

If the per-block deviation is bounded by O(1) (deviating from the
integer truth by less than 0.5 cumulatively), we can ROUND to the
nearest integer and get pi(x) exactly with only O(log x) blocks.

We test this on small N where we have all zeros.
"""

import numpy as np
import os
import sys
from sympy import primepi, mobius
import mpmath
mpmath.mp.dps = 30


def load_zeros():
    """Load Riemann zeros from project data directory."""
    candidates = [
        "data/zeta_zeros_2000.txt",
        "data/zeta_zeros_1000.txt",
        "data/zeta_zeros_500.txt",
        "data/zeta_zeros_300.txt",
        "data/zeta_zeros_200.txt",
    ]
    for c in candidates:
        if os.path.exists(c):
            xs = []
            with open(c) as f:
                for line in f:
                    line = line.strip()
                    if not line:
                        continue
                    try:
                        xs.append(float(line.split()[-1]))
                    except Exception:
                        pass
            if xs:
                print(f"Loaded {len(xs)} zeros from {c}")
                return xs
    # Fallback: compute first 100 zeros via mpmath
    print("Computing first 100 zeros via mpmath.zetazero ...")
    return [float(mpmath.zetazero(k).imag) for k in range(1, 101)]


def Rfunc(x, terms=40):
    """Riemann's R function = sum_{n>=1} mu(n)/n * li(x^{1/n})."""
    s = mpmath.mpf(0)
    for n in range(1, terms + 1):
        m = int(mobius(n))
        if m == 0:
            continue
        s += mpmath.mpf(m) / n * mpmath.li(mpmath.mpf(x) ** (mpmath.mpf(1) / n))
    return s


def Rfunc_complex(z, terms=40):
    """Riemann's R extended to complex z: sum_{n>=1} mu(n)/n * li(z^{1/n}).
    Used for the standard explicit formula pi(x) = R(x) - sum_rho R(x^rho)."""
    s = mpmath.mpc(0)
    for n in range(1, terms + 1):
        m = int(mobius(n))
        if m == 0:
            continue
        s += mpmath.mpf(m) / n * mpmath.li(z ** (mpmath.mpf(1) / n))
    return s


def explicit_pi(x, gamma_list, num_zeros):
    """Approximate pi(x) via the (rounded) explicit formula using
    `num_zeros` zeros from the list. Form: pi(x) = R(x) - sum_rho R(x^rho)."""
    x = mpmath.mpf(x)
    R = Rfunc(x)
    s = mpmath.mpf(0)
    for g in gamma_list[:num_zeros]:
        rho = mpmath.mpc("0.5", g)
        s += 2 * mpmath.re(Rfunc_complex(x ** rho))
    return R - s


def main():
    print("=== Semiclassical orbit resummation of explicit formula ===\n")

    gammas = load_zeros()
    print(f"\nUsing {len(gammas)} zeros, max gamma = {gammas[-1]:.2f}\n")

    test_xs = [10, 100, 500, 1000, 5000, 10000, 50000, 100000]
    print(f"{'x':>8} {'true_pi':>10} {'R(x)':>10} {'used_zeros':>10}  {'approx':>15}  {'err':>10}  {'rounded_err':>10}")

    for x in test_xs:
        truth = int(primepi(x))
        R = float(Rfunc(x))
        # Test: use up to N_max zeros and evaluate
        # Heuristic threshold for needed zeros: T ~ sqrt(x)
        need = int(np.sqrt(x)) + 10
        used = min(need, len(gammas))
        approx = float(explicit_pi(x, gammas, used))
        rounded = int(round(approx))
        err = approx - truth
        rerr = rounded - truth
        print(f"{x:>8} {truth:>10} {R:>10.2f} {used:>10}  {approx:>15.4f}  {err:>10.4f}  {rerr:>10}")

    print("\nNow test BLOCK-RESUMMATION: group zeros into blocks and sum each block,")
    print("comparing per-block contributions to the BLOCK-AVERAGE smooth approximation.")

    # Block analysis: split first num_zeros zeros into log2(num) blocks
    # Block i: zeros with index in [2^i, 2^{i+1}).
    # Compute per-block deviation from expected smooth integral.
    x = 10000
    truth = int(primepi(x))
    print(f"\nFor x = {x}, true pi(x) = {truth}.")
    blocks = []
    i = 0
    while (1 << i) < len(gammas):
        a = 1 << i
        b = min(1 << (i + 1), len(gammas))
        blocks.append((a, b))
        i += 1
    contrib_block = []
    for (a, b) in blocks:
        s = mpmath.mpf(0)
        for g in gammas[a:b]:
            rho = mpmath.mpc("0.5", g)
            s += 2 * mpmath.re(Rfunc_complex(mpmath.mpf(x) ** rho))
        contrib_block.append(float(s))

    print(f"\n{'block_idx':>10} {'zeros':>15} {'gamma_range':>20} {'block_sum':>15} {'cumulative':>15}")
    cum = 0.0
    for k, (a, b) in enumerate(blocks):
        cum += contrib_block[k]
        gr = f"[{gammas[a]:.1f},{gammas[b-1]:.1f}]"
        print(f"{k:>10} {f'[{a},{b})':>15} {gr:>20} {contrib_block[k]:>15.4f} {cum:>15.4f}")

    print(f"\nR(x) - cumulative_zero_contribution_after_log2_blocks:")
    for k_use in range(1, len(blocks) + 1):
        s = sum(contrib_block[:k_use])
        approx = float(Rfunc(x)) - s
        last_g = gammas[blocks[k_use - 1][1] - 1]
        print(f"  use {k_use} blocks (zeros up to gamma~{last_g:.1f}):  approx={approx:.4f} round={round(approx)} truth={truth} err={round(approx)-truth}")

    print("\nVerdict heuristic:")
    print("  - If per-block contributions decay -> finite blocks suffice -> polylog OK.")
    print("  - If per-block contributions oscillate as O(1/sqrt(B)) (random walk) ->")
    print("    Omega(sqrt(x)/log) zeros needed and grouping doesn't save anything.")


if __name__ == "__main__":
    main()
