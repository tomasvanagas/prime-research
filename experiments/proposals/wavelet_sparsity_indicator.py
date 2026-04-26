"""Proposal B: Wavelet sparsity of the prime indicator function.

Test whether 1_P[1..N] (prime indicator) admits a sparse representation
in some wavelet basis. Sparsity-K representations would allow
compressed-sensing recovery in O(K^2 polylog N) per query.

We measure k_99: the number of wavelet coefficients (sorted by magnitude)
needed to capture 99% of the signal energy, for several wavelet families.
A sparse signal has k_99 = polylog(N); a dense signal has k_99 ~ N.

Run: python wavelet_sparsity_indicator.py
"""
from __future__ import annotations

import math
import sys

import numpy as np

try:
    import pywt
except ImportError:
    pywt = None


def sieve(N: int) -> np.ndarray:
    """Return float array of length N where entry k is 1 iff k+1 is prime."""
    sv = np.ones(N + 1, dtype=bool)
    sv[:2] = False
    for k in range(2, int(math.isqrt(N)) + 1):
        if sv[k]:
            sv[k * k :: k] = False
    return sv[1 : N + 1].astype(np.float64)


def k99_haar_manual(f: np.ndarray) -> tuple[int, int]:
    """Haar transform via repeated averaging — works without pywt."""
    coefs: list[float] = []
    cur = f.copy()
    while len(cur) > 1:
        even = cur[0::2]
        odd = cur[1::2]
        avg = (even + odd) / math.sqrt(2.0)
        diff = (even - odd) / math.sqrt(2.0)
        coefs.extend(diff.tolist())
        cur = avg
    coefs.extend(cur.tolist())
    arr = np.asarray(coefs)
    energy = arr * arr
    if energy.sum() == 0.0:
        return 0, len(arr)
    idx = np.argsort(energy)[::-1]
    cumulative = np.cumsum(energy[idx]) / energy.sum()
    k99 = int(np.searchsorted(cumulative, 0.99) + 1)
    return k99, len(arr)


def k99_pywt(f: np.ndarray, wavelet: str, level: int) -> tuple[int, int]:
    coeffs = pywt.wavedec(f, wavelet, level=level)
    flat = np.concatenate([c.ravel() for c in coeffs])
    energy = flat * flat
    if energy.sum() == 0.0:
        return 0, len(flat)
    idx = np.argsort(energy)[::-1]
    cumulative = np.cumsum(energy[idx]) / energy.sum()
    k99 = int(np.searchsorted(cumulative, 0.99) + 1)
    return k99, len(flat)


def random_density_matched(f: np.ndarray, seed: int = 0) -> np.ndarray:
    rng = np.random.default_rng(seed)
    g = np.zeros_like(f)
    n_ones = int(f.sum())
    idx = rng.choice(len(f), size=n_ones, replace=False)
    g[idx] = 1.0
    return g


def sin_signal(N: int, freq: float = 0.01) -> np.ndarray:
    x = np.arange(N)
    return 0.5 + 0.5 * np.sin(2 * math.pi * freq * x)


def main() -> None:
    powers = [10, 12, 14]
    print(f"# Wavelet sparsity of prime indicator")
    print(f"# k_99: number of largest-magnitude coefficients holding 99% energy")
    print()
    for p in powers:
        N = 2**p
        f = sieve(N)
        n_primes = int(f.sum())
        print(f"## N = 2^{p} = {N} (pi(N) = {n_primes})")

        k99_h, total_h = k99_haar_manual(f)
        print(f"  Haar (manual)    : k_99 = {k99_h:>6} / {total_h}  ratio = {k99_h/total_h:.4f}")

        # Compare against random density-matched signal
        g = random_density_matched(f)
        k99_rg, total_rg = k99_haar_manual(g)
        print(f"  Haar random-match: k_99 = {k99_rg:>6} / {total_rg}  ratio = {k99_rg/total_rg:.4f}")

        # And against a smooth sinusoid (ground-truth sparse signal)
        s = sin_signal(N)
        k99_s, total_s = k99_haar_manual(s)
        print(f"  Haar sinusoid    : k_99 = {k99_s:>6} / {total_s}  ratio = {k99_s/total_s:.4f}")

        if pywt is not None:
            for wav in ("db4", "db8", "sym8"):
                try:
                    level = min(int(math.log2(N)) - 3, pywt.dwt_max_level(N, wav))
                    k99, total = k99_pywt(f, wav, level)
                    print(f"  {wav:>3} (lvl={level})    : k_99 = {k99:>6} / {total}  ratio = {k99/total:.4f}")
                except Exception as exc:  # noqa: BLE001
                    print(f"  {wav}: skipped ({exc})")
        else:
            print("  pywt not installed; only Haar reported")

        log2N = math.log2(N)
        polylog_target = int(log2N**3)
        print(f"  polylog target ((log2 N)^3 = {polylog_target}) — sparsity benchmark")
        print()


if __name__ == "__main__":
    sys.exit(main())
