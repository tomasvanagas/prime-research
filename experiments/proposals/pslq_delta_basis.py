"""Proposal C: PSLQ + structured-basis search on delta(n).

We compute delta(n) = p(n) - R^{-1}(n) at high precision for n = 1..N
and then apply PSLQ to find integer relations between delta(n) and a
basis of "natural atoms":

    1                                            (constant offset)
    log(n), log log(n), 1/log(n)                 (smooth scales)
    sqrt(n)/log(n)                               (PNT correction)
    sin(2 pi gamma_k log n) / sqrt(n)            (zeta-zero oscillations)
    cos(2 pi gamma_k log n) / sqrt(n)            (k = 1..K)

If PSLQ returns a relation of moderate norm (< 10^6) for several values
of n simultaneously, that is evidence delta(n) admits an integer linear
combination of these atoms — i.e. a closed form.

Equivalently we test the *small-norm null vector* hypothesis on the
matrix M[i, j] = atom_j(n_i) augmented by a column of delta(n_i): does
the matrix have a small-norm integer null vector?

Run: python3 pslq_delta_basis.py
"""
from __future__ import annotations

import math
import sys
import time
from pathlib import Path

import numpy as np
from mpmath import mp, mpf, pslq

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "algorithms"))
from v5_pure_analytic import _R, _R_inverse  # type: ignore  # noqa: E402


def load_zeros(K: int) -> list[float]:
    p = REPO / "data" / "zeta_zeros_2000.txt"
    if not p.exists():
        p = REPO / "data" / "zeta_zeros_1000.txt"
    out = []
    with p.open() as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            out.append(float(line))
            if len(out) >= K:
                break
    return out


def sieve_primes(N: int) -> list[int]:
    sv = [True] * (N + 1)
    sv[:2] = [False, False]
    for k in range(2, int(math.isqrt(N)) + 1):
        if sv[k]:
            for j in range(k * k, N + 1, k):
                sv[j] = False
    return [i for i, b in enumerate(sv) if b]


def compute_delta_hp(n_max: int, dps: int = 50) -> tuple[list[mpf], list[int]]:
    mp.dps = dps
    primes = sieve_primes(20 * n_max + 200)
    deltas: list[mpf] = []
    pn: list[int] = []
    for n in range(1, n_max + 1):
        pn.append(primes[n - 1])
        # Use the mp-aware R-inverse: bisect using R(x) - n
        # For small n we just use the plain double-precision R_inverse
        # since 30 digits is plenty here.
        rinv = mpf(_R_inverse(n))
        deltas.append(mpf(primes[n - 1]) - rinv)
    return deltas, pn


def build_basis(n: int, gammas: list[float]) -> list[mpf]:
    n_mp = mpf(n + 1)  # shift to avoid log(1) = 0
    log_n = mp.log(n_mp)
    row: list[mpf] = []
    row.append(mpf(1))
    row.append(log_n)
    row.append(mp.log(log_n) if n >= 2 else mpf(0))
    row.append(1 / log_n)
    row.append(mp.sqrt(n_mp) / log_n)
    sqn = mp.sqrt(n_mp)
    for g in gammas:
        ang = 2 * mp.pi * mpf(g) * log_n / (2 * mp.pi)  # = gamma * log_n
        # Standard zero contribution oscillates as cos(gamma log n) / sqrt(n)
        row.append(mp.cos(mpf(g) * log_n) / sqn)
        row.append(mp.sin(mpf(g) * log_n) / sqn)
    return row


def main() -> None:
    print("# PSLQ on delta(n) with structured basis")
    mp.dps = 40
    N_TARGET = 50  # primes considered
    K_GAMMAS = 5
    gammas = load_zeros(K_GAMMAS)
    print(f"# Using first {K_GAMMAS} zeros: {gammas[:K_GAMMAS]}")

    delta, pn = compute_delta_hp(N_TARGET, dps=40)
    print(f"# delta(1..10) = {[float(d) for d in delta[:10]]}")
    print()

    # ------------------------------------------------------------------
    # PSLQ test: for each n, augment basis with delta(n) and seek
    # integer relations across multiple n. We do this by "row-stacking".
    # PSLQ requires a single vector; we adapt by averaging atoms over n
    # which is a weak test. Below is a stronger test: we form a matrix
    # M of size (N_TARGET, M_basis+1) and ask whether the column vector
    # delta is in the integer span of the basis columns at low height.
    # ------------------------------------------------------------------

    rows = []
    for n in range(1, N_TARGET + 1):
        rows.append(build_basis(n, gammas))
    M = len(rows[0])
    A = np.array([[float(x) for x in r] for r in rows])  # double precision matrix
    y = np.array([float(d) for d in delta])

    # Least-squares fit
    coef, residuals, rank, sv = np.linalg.lstsq(A, y, rcond=None)
    pred = A @ coef
    err = y - pred
    rms = math.sqrt(float(np.mean(err * err)))
    max_abs = float(np.max(np.abs(err)))
    print("## Least-squares fit:")
    print(f"  coefficients (size {M}): {coef.round(4).tolist()}")
    print(f"  RMS error  = {rms:.4f}")
    print(f"  max |err|  = {max_abs:.4f}")
    print(f"  cond(A)    = {np.linalg.cond(A):.2e}")
    print()

    # ------------------------------------------------------------------
    # PSLQ-style: try to find integer relation among (delta_n_avg, basis_avgs).
    # We compute averages of each basis fn over n=1..N to get a single vector,
    # then call mpmath.pslq.
    # ------------------------------------------------------------------
    avg_delta = sum(delta) / N_TARGET
    avg_basis: list[mpf] = []
    for j in range(M):
        col = [rows[i][j] for i in range(N_TARGET)]
        avg_basis.append(sum(col) / N_TARGET)

    vec = [avg_delta] + avg_basis
    print("## PSLQ on column-averaged vector:")
    print(f"  Vector length = {len(vec)}")
    try:
        rel = pslq(vec, tol=1e-15, maxcoeff=10**6)
    except Exception as exc:  # noqa: BLE001
        rel = None
        print(f"  PSLQ exception: {exc}")
    if rel is None:
        print("  No relation of norm <= 10^6 found.")
    else:
        print(f"  Relation: {rel}")
        # Verify magnitude
        s = sum(int(c) * v for c, v in zip(rel, vec))
        print(f"  Residual: {float(s):.3e}")

    # ------------------------------------------------------------------
    # Stronger test: PSLQ across multiple n simultaneously. We seek a
    # single coefficient vector v such that A v = delta. We linearize:
    # for each row i, build (delta[i], A[i,:]) and run PSLQ. Then average
    # the relation found and check stability.
    # ------------------------------------------------------------------
    print("\n## PSLQ on individual rows (looking for stable relation):")
    relations: list[list[int]] = []
    for i in range(0, N_TARGET, 5):
        try:
            v = [delta[i]] + rows[i]
            r = pslq(v, tol=1e-18, maxcoeff=10**5)
        except Exception:  # noqa: BLE001
            r = None
        if r is not None:
            relations.append(list(r))
            print(f"  n={i+1}: {r}")
        else:
            print(f"  n={i+1}: no relation found")

    if len(relations) >= 2:
        # Check whether the same first coefficient appears stably
        first_coefs = [r[0] for r in relations]
        from collections import Counter
        ct = Counter(first_coefs)
        print(f"  delta-coefficients frequency: {ct}")


if __name__ == "__main__":
    main()
