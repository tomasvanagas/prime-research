"""D13 — Subword complexity / topological entropy of chi_P as binary infinite word.

For w_N = (chi_P(k))_{k=2..N} in {0,1}^*, compute the subword complexity
    p_w(n) := #{distinct length-n factors of w_N}
for n = 1..n_max, on multiple chi_P-derived streams:

  RAW  : raw chi_P(2..N) — dominated by trivial parity constraint
         (even positions never prime except 2).
  ODD  : odd-position-only stream chi'_P(k) := chi_P(2k+1) for k=1..(N-1)/2
         — removes parity constraint; tests for non-parity structure.
  W{q} : Green-Tao W-trick with W = q (primorial), r coprime to W.
         Stream is (chi_P(W*n + r))_{n=0,1,...} with r=1, fixed.
         Density 1/phi(W) and removes ALL small-prime mod-q constraints
         for primes <= p_k. q in {6, 30, 210}.

Baselines per stream:
  PERM : K random permutations of the stream (preserves density and
         multiset of bits exactly; destroys all >1-gram structure).
  BERN : K independent iid Bernoulli at matched density rho.

These two are the F3-style "matched-density same-resolution" controls
borrowed from S96 (D2 persistent homology). PERM is the strict null
because it preserves density precisely; BERN gives a cross-check that
any density-fluctuation effect doesn't bias the result.

Topological entropy (Lind-Marcus 1995, Cassaigne-Nicolas 2010,
Wikipedia "Complexity function"):
    h_eff(n) := log_2(p_w(n)) / n
    h_w := lim_{n->infty} h_eff(n)   for binary alphabet, h_w in [0, 1].

Pre-registered F3 falsifier:
    "PRIMES > 3 sigma below mean of BOTH PERM AND BERN at every n in
    [n_lo, n_hi]" -> structural deviation; otherwise -> at noise floor.

A non-trivial deviation must hold on AT LEAST ONE of {ODD, W{6}, W{30}}
streams to be considered an arithmetic edge (raw is parity-trivial).

Cross-domain refs:
  - Lind-Marcus 1995 *An Introduction to Symbolic Dynamics and Coding* CUP
  - Cassaigne-Nicolas 2010 "Factor complexity" in CANT vol. 135
  - Morse-Hedlund 1938-40 *Amer. J. Math.* 60
  - Wikipedia: Complexity function
"""
import argparse
import json
import os
import time

import numpy as np


# ---------- Streams ----------

def sieve(N: int) -> np.ndarray:
    """chi_P[k] = 1 iff k prime, for k = 0..N."""
    s = np.ones(N + 1, dtype=np.uint8)
    s[0] = 0
    s[1] = 0
    for p in range(2, int(N**0.5) + 1):
        if s[p]:
            s[p * p::p] = 0
    return s


def stream_raw(chi: np.ndarray) -> np.ndarray:
    """RAW chi_P stream: chi_P(2..N) as uint8."""
    return chi[2:].copy()


def stream_odd(chi: np.ndarray) -> np.ndarray:
    """ODD stream: chi_P(2k+1) for k = 1, 2, ...
    First element corresponds to n=3 (smallest odd >= 3)."""
    N = len(chi) - 1  # chi indexed 0..N
    # take values at indices 3, 5, 7, ..., (largest odd <= N)
    odd_idx = np.arange(3, N + 1, 2)
    return chi[odd_idx].copy()


def stream_wtrick(chi: np.ndarray, W: int, r: int = 1) -> np.ndarray:
    """W-tricked stream: chi_P(W*n + r) for n = 0, 1, 2, ...
    First element corresponds to W*0+r = r."""
    N = len(chi) - 1
    n_idx = np.arange(0, (N - r) // W + 1)
    indices = W * n_idx + r
    return chi[indices].copy()


# ---------- Subword complexity ----------

def factor_count(bits: np.ndarray, n: int) -> int:
    """Return p_w(n) = #{distinct length-n factors of bits} via vectorized
    rolling encoding. Memory: ~8*L bytes for the encoding array."""
    L = len(bits)
    if L < n:
        return 0
    if n > 63:
        raise ValueError("n must be <= 63 to fit in uint64.")
    # Build encoding[i] = sum_{j=0..n-1} bits[i+j] * 2^j
    enc = np.zeros(L - n + 1, dtype=np.uint64)
    bits64 = bits.astype(np.uint64)
    for j in range(n):
        enc += bits64[j:L - n + j + 1] << np.uint64(j)
    return int(np.unique(enc).size)


def topological_entropy_curve(bits: np.ndarray, n_max: int) -> np.ndarray:
    return np.array([factor_count(bits, n) for n in range(1, n_max + 1)],
                    dtype=np.int64)


# ---------- Experiment ----------

def run_stream(name: str, bits: np.ndarray, n_max: int, K: int,
               rng: np.random.Generator) -> dict:
    L = len(bits)
    rho = float(bits.sum()) / L if L > 0 else 0.0
    print(f"\n=== Stream: {name} (L={L}, density rho={rho:.6f}) ===",
          flush=True)
    t0 = time.time()
    p_chi = topological_entropy_curve(bits, n_max)
    print(f"  p_chi computed in {time.time() - t0:.2f}s", flush=True)

    # Permutation controls
    p_perm = np.zeros((K, n_max), dtype=np.int64)
    t0 = time.time()
    for k in range(K):
        perm_bits = rng.permutation(bits)
        p_perm[k] = topological_entropy_curve(perm_bits, n_max)
    print(f"  K={K} permutation controls in {time.time() - t0:.2f}s",
          flush=True)
    perm_mean = p_perm.mean(axis=0)
    perm_std = p_perm.std(axis=0, ddof=1)

    # Bernoulli matched-density controls
    p_bern = np.zeros((K, n_max), dtype=np.int64)
    t0 = time.time()
    for k in range(K):
        bern_bits = (rng.random(L) < rho).astype(np.uint8)
        p_bern[k] = topological_entropy_curve(bern_bits, n_max)
    print(f"  K={K} bern controls in {time.time() - t0:.2f}s", flush=True)
    bern_mean = p_bern.mean(axis=0)
    bern_std = p_bern.std(axis=0, ddof=1)

    # Z-scores (avoid divide-by-zero with small-n constants)
    with np.errstate(divide="ignore", invalid="ignore"):
        z_perm = np.where(perm_std > 0,
                          (p_chi - perm_mean) / perm_std, 0.0)
        z_bern = np.where(bern_std > 0,
                          (p_chi - bern_mean) / bern_std, 0.0)
    h_eff = np.log2(np.maximum(p_chi, 1)) / np.arange(1, n_max + 1)
    h_eff_perm = np.log2(np.maximum(perm_mean, 1)) / np.arange(1, n_max + 1)

    print(f"  {'n':>3} | {'p_chi':>10} | {'p_perm':>12} | {'sig_perm':>9} |"
          f" {'z_perm':>8} | {'p_bern':>12} | {'z_bern':>8} | {'h_chi':>6} | {'h_perm':>6}",
          flush=True)
    for n in range(1, n_max + 1):
        print(f"  {n:>3} | {p_chi[n-1]:>10d} | {perm_mean[n-1]:>12.1f} | "
              f"{perm_std[n-1]:>9.2f} | {z_perm[n-1]:>8.2f} | "
              f"{bern_mean[n-1]:>12.1f} | {z_bern[n-1]:>8.2f} | "
              f"{h_eff[n-1]:.4f} | {h_eff_perm[n-1]:.4f}",
              flush=True)

    return dict(
        name=name, L=L, rho=rho, n_max=n_max, K=K,
        p_chi=p_chi.tolist(),
        p_perm_mean=perm_mean.tolist(),
        p_perm_std=perm_std.tolist(),
        p_perm_all=p_perm.tolist(),
        z_perm=z_perm.tolist(),
        p_bern_mean=bern_mean.tolist(),
        p_bern_std=bern_std.tolist(),
        p_bern_all=p_bern.tolist(),
        z_bern=z_bern.tolist(),
        h_eff=h_eff.tolist(),
        h_eff_perm=h_eff_perm.tolist(),
    )


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--N", type=int, default=10_000_000)
    ap.add_argument("--n_max", type=int, default=24)
    ap.add_argument("--K", type=int, default=20)
    ap.add_argument("--Ws", type=str, default="6,30,210",
                    help="comma-separated W-trick moduli")
    ap.add_argument("--out", type=str,
                    default="/apps/aplikacijos/prime-research/experiments/dynamical/subword_complexity_chi_p/results.json")
    ap.add_argument("--streams", type=str, default="raw,odd,wtrick",
                    help="comma-separated subset of {raw,odd,wtrick}")
    args = ap.parse_args()

    print(f"[t={time.time():.1f}] Building chi_P up to N={args.N}...",
          flush=True)
    t0 = time.time()
    chi = sieve(args.N)
    print(f"  sieve cost: {time.time() - t0:.2f}s", flush=True)

    rng = np.random.default_rng(seed=42)
    streams_to_run = set(args.streams.split(","))

    out = {"N": args.N, "n_max": args.n_max, "K": args.K, "streams": {}}

    if "raw" in streams_to_run:
        out["streams"]["raw"] = run_stream(
            "raw", stream_raw(chi), args.n_max, args.K, rng)

    if "odd" in streams_to_run:
        out["streams"]["odd"] = run_stream(
            "odd", stream_odd(chi), args.n_max, args.K, rng)

    if "wtrick" in streams_to_run:
        for W_str in args.Ws.split(","):
            W = int(W_str.strip())
            out["streams"][f"W{W}"] = run_stream(
                f"W{W}", stream_wtrick(chi, W, r=1),
                args.n_max, args.K, rng)

    os.makedirs(os.path.dirname(args.out), exist_ok=True)
    with open(args.out, "w") as f:
        json.dump(out, f, indent=2)
    print(f"\nWrote {args.out}")


if __name__ == "__main__":
    main()
