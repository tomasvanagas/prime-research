"""Proposal C — Algebraicity of Riemann zero differences.

Tests whether {gamma_k - gamma_1} admit integer relations of bounded
height against a small Q-basis. A success would dramatically compress
the explicit formula; a clean miss (0/N) is a strong negative.

Prints results and writes to proposal_c_zero_algebraicity_results.md is
the responsibility of the human; this script just dumps to stdout.
"""

from pathlib import Path

from mpmath import mp, mpf, pi, log, gamma, pslq

mp.dps = 50

ZEROS_FILE = Path(__file__).resolve().parents[2] / "data" / "zeta_zeros_8000.txt"
N_PROBES   = 200            # try gamma_2 .. gamma_200
HEIGHT_CAP = 10**15
TOL        = mpf(10) ** -40


def load_zeros(n: int) -> list:
    out = []
    with open(ZEROS_FILE) as f:
        for line in f:
            out.append(mpf(line.strip()))
            if len(out) == n:
                break
    return out


def main() -> None:
    zeros = load_zeros(N_PROBES + 1)
    basis = [
        mpf(1),
        pi,
        log(2),
        log(3),
        log(5),
        log(7),
        log(11),
        log(13),
        log(gamma(mpf(1) / 4)),
        log(gamma(mpf(1) / 3)),
        log(gamma(mpf(1) / 6)),
        log(pi),
    ]
    basis_names = [
        "1",
        "pi",
        "log2",
        "log3",
        "log5",
        "log7",
        "log11",
        "log13",
        "log G(1/4)",
        "log G(1/3)",
        "log G(1/6)",
        "log pi",
    ]
    print(f"# Proposal C — algebraicity of zero differences")
    print(f"basis size = {len(basis)}, height cap = {HEIGHT_CAP}, "
          f"tol = 1e-40, dps = {mp.dps}")
    print()
    print(f"{'k':>4} {'gamma_k-gamma_1':>26} | result")
    print("-" * 70)
    hits = 0
    for k in range(2, N_PROBES + 1):
        delta = zeros[k - 1] - zeros[0]
        try:
            rel = pslq([delta] + basis, tol=TOL, maxcoeff=HEIGHT_CAP)
        except Exception as exc:
            print(f"{k:>4} {str(delta)[:26]:>26} | err: {exc}")
            continue
        if rel is None:
            tag = "."
        else:
            terms = []
            for c, name in zip(rel[1:], basis_names):
                if c != 0:
                    terms.append(f"{c:+d}*{name}")
            tag = f"REL c0={rel[0]:+d}: " + " ".join(terms)
            hits += 1
        print(f"{k:>4} {str(delta)[:26]:>26} | {tag}")
    print()
    print(f"hits: {hits} / {N_PROBES - 1}")
    print(f"verdict: {'CANDIDATE' if hits >= 5 else 'NEGATIVE'}")


if __name__ == "__main__":
    main()
