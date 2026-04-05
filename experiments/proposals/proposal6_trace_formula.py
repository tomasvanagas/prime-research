"""
PROPOSAL 6: Trace Formula with Optimal Test Function
"""
import numpy as np
from mpmath import mp, im, re, zetazero
mp.dps = 30

K = 20
zeros = [(float(re(zetazero(k))), float(im(zetazero(k)))) for k in range(1, K+1)]

print("=" * 70)
print("PROPOSAL 6: Trace Formula with Optimal Test Function")
print("=" * 70)

print("\n--- Decay of zero contributions with Gaussian smoothing ---")
for x in [100, 1000, 10000, 100000]:
    sigma = np.log(x)
    print(f"\nx={x}, sigma={sigma:.2f}:")
    for k in range(min(10, K)):
        gamma = zeros[k][1]
        decay = np.exp(-sigma**2 * gamma**2 / 2)
        print(f"  zero {k+1} (gamma={gamma:.2f}): decay = {decay:.2e}")

print("\n\n--- Sharp cutoff: how many zeros for error < 1? ---")
for x in [1000, 10000, 100000, 10**6, 10**9, 10**12, 10**20, 10**50, 10**100]:
    half_root = np.sqrt(float(x))
    T_needed = half_root / np.pi
    if T_needed < 1:
        K_needed = 0
    else:
        K_needed = int(T_needed * np.log(max(2,T_needed)) / (2 * np.pi))
    print(f"  x={x:.0e}: need ~{K_needed:>20} zeros (T_needed={T_needed:.1e})")

print("\n\n--- IDEA: Can we avoid individual zeros entirely? ---")
print("The Riemann-von Mangoldt formula: N(T) = T/(2*pi)*log(T/(2*pi*e)) + O(log T)")
print("Total contribution of all zeros with |gamma| in [T, T+1]:")
print("  ~ sqrt(x) * log(T) / T  (each zero contributes ~ x^{1/2}/|gamma|)")
print("Sum over T from 1 to T_max: ~ sqrt(x) * log^2(T_max)")
print("For error < 1: T_max ~ sqrt(x) * log^2(x)")
print("Number of zeros below T_max: ~ T_max * log(T_max) / (2*pi) ~ sqrt(x) * log^3(x)")
print()
print("CONCLUSION: sqrt(x) * polylog(x) zeros always needed for the explicit formula.")
print("The Gaussian smoothing doesn't help: smoother test function = less resolution.")
print()
print("HOWEVER: What if we DON'T need individual zeros but only STATISTICS?")
print("Random Matrix Theory: the zeros locally behave like GUE eigenvalues.")
print("Their PAIR CORRELATION is universal. Could we use this to shortcut the sum?")
