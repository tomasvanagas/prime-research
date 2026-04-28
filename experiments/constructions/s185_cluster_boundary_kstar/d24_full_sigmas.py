"""Compute full sigma spectrum at d=24 to enable MP-edge counting beyond top-100."""
import json
import os
import numpy as np
from sympy import primerange

d = 24
N = 1 << d
M = 1 << (d // 2)
print(f"d={d}, N={N}, M={M}")

print("Sieving primes...")
arr = np.zeros(N, dtype=np.float64)
for p in primerange(2, N + 1):
    arr[p - 1] = 1.0

print("Reshaping and SVD (sigmas only)...")
A = arr.reshape(M, M)
sigmas = np.linalg.svd(A, compute_uv=False)
print(f"Total sigmas: {len(sigmas)}")
print(f"Frobenius_sq: {np.sum(sigmas**2):.1f}")
print(f"sigma_0: {sigmas[0]:.2f}")
print(f"sigma_99: {sigmas[99]:.4f}")
print(f"sigma_500: {sigmas[500]:.4f}")
print(f"sigma_1000: {sigmas[1000]:.4f}")
print(f"sigma_-1: {sigmas[-1]:.4f}")

# MP edge at d=24
pi_N = int(np.sum(arr))
p_N = pi_N / N
edge = 2.0 * np.sqrt(M * p_N * (1.0 - p_N))
print(f"\npi(N) = {pi_N}, MP edge = {edge:.4f}")

# Count sigmas above MP edge
above = np.sum(sigmas > edge)
spike_block_mp = float(np.sum(sigmas[1:above] ** 2))
print(f"#sigmas above MP edge: {above}")
print(f"spike_block (MP edge, exc k=0): {spike_block_mp:.1f}")
print(f"spike_block / pi(N): {spike_block_mp / pi_N:.4f}")

# Look at sigma values in transition zone
print("\nSigma values in MP transition zone:")
for k in [50, 75, 100, 150, 200, 300, 500, 700, 1000, 1500, 2000, 3000, 4000]:
    if k < len(sigmas):
        print(f"  sigma_{k} = {sigmas[k]:.4f}")

# Compute spike_block at multiple k_* values
print("\nspike_block / pi(N) at various k_*:")
for k_star in [50, 75, 78, 100, 150, 200, 300, 500, 1000]:
    if k_star < len(sigmas):
        block = float(np.sum(sigmas[1:k_star + 1] ** 2))
        print(f"  k_* = {k_star}: spike_block = {block:.1f}, frac = {block / pi_N:.4f}")

# Save
out = {
    "N": N,
    "M": M,
    "pi_N": pi_N,
    "MP_edge": edge,
    "k_star_MP_edge": int(above),
    "spike_block_MP_edge": spike_block_mp,
    "frac_MP_edge": spike_block_mp / pi_N,
    "sigmas_top500": sigmas[:500].tolist(),
    "frobenius_sq": float(np.sum(sigmas ** 2)),
    "frac_at_k": {
        str(k): float(np.sum(sigmas[1:k+1] ** 2) / pi_N)
        for k in [10, 20, 30, 40, 50, 60, 70, 78, 80, 90, 100, 120, 150, 200, 300, 500, 1000, 2000]
    },
}
with open(os.path.join(os.path.dirname(__file__), "d24_full_sigmas.json"), "w") as f:
    json.dump(out, f, indent=2, default=float)
print(f"\nWrote d24_full_sigmas.json")
