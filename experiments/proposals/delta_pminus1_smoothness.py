"""
Conditional structure of delta(n) = p(n) - round(R^{-1}(n))
on the SMOOTHNESS of p(n) - 1.

Motivation
----------
Most prior delta(n) tests look at delta as a sequence in n (PSLQ, AR, MPS,
spectral, autocorrelation, conditional entropy on n mod m). All close as
mode-I (incompressible). One angle not yet probed in this exact form:
condition delta on a FACTORIZATION feature of the prime itself.

The natural arithmetic feature of a prime p is the smoothness of p-1
(Pollard p-1, Pocklington primality, Lucas-Lehmer machinery all exploit
this). If primes with B-smooth p-1 had detectably different delta
statistics, that would be a fresh structural signal worth understanding.

Hypothesis: H_0 = "delta(pi(p)) is independent of s(p-1)" where
s(m) = largest prime factor of m. We test by partitioning primes p<=X
into smoothness classes and comparing delta statistics class by class.

Output: per-class mean(delta), std(delta), sign(delta) frequency,
Pearson correlation of |delta| with log(s(p-1))/log(p), KS test
between smoothest and least-smooth quartiles.

Closes either:
  (a) FAIL/I (no smoothness-conditional structure) -- adds 35th measure
      to novel/pseudorandomness_of_pi.md, OR
  (b) PARTIAL (real correlation) -- promote to a follow-up exact-rounding
      test to see if the signal is exploitable for prediction.

Implementation note
-------------------
We use the existing fx_data.npz containing f(x) = pi(x) - R(x) for
x = 2..1e5. Since R(p) = n - f(p) for n = pi(p), and locally
R(x+dx) approx R(p) + dx/log(p), we get
    R^{-1}(n) approx p + f(p) * log(p),
so delta(n) = p - round(R^{-1}(n)) = -round(f(p) * log(p)).
This avoids per-prime Newton iteration in mpmath while remaining
exact-rounding-faithful in the regime where |f(p) log(p)| << p
(verified true for x <= 1e5: max|f(p) log(p)| ~ 11).
"""

import numpy as np
from sympy import sieve, factorint
from scipy.stats import pearsonr, spearmanr, ks_2samp

DATA = np.load("/apps/aplikacijos/prime-research/experiments/algebraic/identity_search/fx_data.npz")
xs = DATA["x"]      # x = 2,3,...,1e5
fs = DATA["f"]      # f(x) = pi(x) - R(x)
pis = DATA["pi"]    # pi(x)
print(f"Loaded f(x) for x in [{xs[0]},{xs[-1]}], N={len(xs)}")

# Find primes <= max x
X = int(xs[-1])
primes = list(sieve.primerange(2, X + 1))
N_p = len(primes)
print(f"Primes in range: {N_p}")

# For each prime, get its f(p) value (xs is contiguous starting at 2 -> index = p-2)
prime_to_idx = {p: p - 2 for p in primes}

# delta(n) = -round(f(p) * log(p)) where n = pi(p)
log_p = np.log(np.array(primes, dtype=np.float64))
f_at_p = np.array([fs[prime_to_idx[p]] for p in primes])
delta_real = -f_at_p * log_p
delta = np.round(delta_real).astype(np.int64)

print(f"delta statistics (N={N_p}):")
print(f"  mean={delta.mean():.4f}, std={delta.std():.4f}, "
      f"min={delta.min()}, max={delta.max()}")
print(f"  P(delta=0) = {(delta==0).mean():.4f}")
print(f"  P(delta>0) = {(delta>0).mean():.4f}")

# Cross-check: validate against direct computation for first 50 primes
import mpmath
mpmath.mp.dps = 30
print("\nValidation: comparing 'delta = -round(f log p)' approximation "
      "to direct Newton-Rinv for first 30 primes...")

def mob(n):
    f = factorint(n)
    if any(v > 1 for v in f.values()):
        return 0
    return (-1) ** len(f)

def R_mp(x):
    s = mpmath.mpf(0)
    for k in range(1, 30):
        m = mob(k)
        if m == 0:
            continue
        s += m / k * mpmath.li(mpmath.power(x, mpmath.mpf(1) / k))
    return s

def Rinv_mp(n, x0):
    x = mpmath.mpf(x0)
    for _ in range(30):
        rx = R_mp(x) - n
        d = 1 / mpmath.log(x)
        dx = rx / d
        x_new = x - dx
        if abs(x_new - x) < mpmath.mpf('1e-12'):
            return x_new
        x = x_new
    return x

mismatches = 0
for i in range(min(30, N_p)):
    p = primes[i]
    n = i + 1
    rinv = Rinv_mp(n, p)
    direct_delta = p - int(mpmath.nint(rinv))
    if direct_delta != delta[i]:
        mismatches += 1
        print(f"  n={n}, p={p}: approx={delta[i]}, direct={direct_delta}, "
              f"f(p)log(p)={delta_real[i]:.3f}")
print(f"Validation mismatches: {mismatches}/30")

# Compute s(p-1) for each prime
print("\nFactoring p-1 ...")
def lpf(m):
    if m <= 1:
        return 1
    return max(factorint(m).keys())

s_pm1 = np.zeros(N_p, dtype=np.int64)
for i in range(N_p):
    if primes[i] == 2:
        s_pm1[i] = 1  # p-1=1 sentinel
    else:
        s_pm1[i] = lpf(primes[i] - 1)
print("Done factoring.")

# Smoothness ratio: log(s(p-1)) / log(p) in (0, 1]
# Larger ratio = "rougher" p-1 (largest prime factor close to p itself)
# Smaller ratio = smoother p-1
mask = np.array([p > 2 for p in primes])
log_s_ratio = np.zeros(N_p)
log_s_ratio[mask] = np.log(s_pm1[mask].astype(np.float64)) / log_p[mask]

# Quartile classes
qs = np.quantile(log_s_ratio[mask], [0.25, 0.5, 0.75])
print(f"\nQuartiles of log(s(p-1))/log(p): {qs}")
bins = np.digitize(log_s_ratio, qs)

print(f"\n{'class':>6} {'count':>8} {'mean(delta)':>12} {'std(delta)':>12} "
      f"{'P(d>0)':>10} {'mean|d|':>10}")
for c in range(4):
    sel = (bins == c) & mask
    if sel.sum() > 0:
        d = delta[sel]
        print(f"{c:>6} {sel.sum():>8} {d.mean():>12.4f} {d.std():>12.4f} "
              f"{(d>0).mean():>10.4f} {np.abs(d).mean():>10.4f}")

# Correlations
r_p, p_p = pearsonr(log_s_ratio[mask], delta[mask])
r_s, p_s = spearmanr(log_s_ratio[mask], delta[mask])
r_pa, p_pa = pearsonr(log_s_ratio[mask], np.abs(delta[mask]).astype(np.float64))
print(f"\nPearson(log_s_ratio, delta):  r={r_p:+.5f}, p={p_p:.3g}")
print(f"Spearman(log_s_ratio, delta): r={r_s:+.5f}, p={p_s:.3g}")
print(f"Pearson(log_s_ratio, |delta|):r={r_pa:+.5f}, p={p_pa:.3g}")

# KS test on smoothest vs least smooth quartiles
sel0 = (bins == 0) & mask
sel3 = (bins == 3) & mask
ks_stat, ks_p = ks_2samp(delta[sel0], delta[sel3])
print(f"\nKS class 0 (smoothest p-1) vs class 3 (largest s(p-1)/p): "
      f"D={ks_stat:.5f}, p={ks_p:.3g}")

# Sophie Germain check
primes_set = set(primes)
sg_idx = [i for i in range(1, N_p)
          if (primes[i] - 1) // 2 in primes_set
          and 2 * ((primes[i] - 1) // 2) + 1 == primes[i]]
print(f"\nSophie-Germain primes p=2q+1 in range: N={len(sg_idx)}")
if len(sg_idx) >= 30:
    sg_d = delta[sg_idx]
    other = [i for i in range(1, N_p) if i not in set(sg_idx)]
    other_d = delta[other]
    ks_stat, ks_p = ks_2samp(sg_d, other_d)
    print(f"  mean(delta)_SG = {sg_d.mean():.4f} vs others {other_d.mean():.4f}")
    print(f"  std(delta)_SG  = {sg_d.std():.4f} vs others {other_d.std():.4f}")
    print(f"  KS p-value     = {ks_p:.3g}")

# Number-of-distinct-prime-factors of p-1 (Erdos-Kac style)
omega_pm1 = np.zeros(N_p, dtype=np.int64)
for i in range(N_p):
    if primes[i] == 2:
        omega_pm1[i] = 0
    else:
        omega_pm1[i] = len(factorint(primes[i] - 1))

r_o, p_o = pearsonr(omega_pm1[mask].astype(np.float64), delta[mask])
r_oa, p_oa = pearsonr(omega_pm1[mask].astype(np.float64),
                      np.abs(delta[mask]).astype(np.float64))
print(f"\nPearson(omega(p-1), delta):  r={r_o:+.5f}, p={p_o:.3g}")
print(f"Pearson(omega(p-1), |delta|):r={r_oa:+.5f}, p={p_oa:.3g}")

# Save data
np.savez_compressed(
    "/apps/aplikacijos/prime-research/experiments/proposals/delta_pminus1_smoothness_data.npz",
    primes=np.array(primes, dtype=np.int64),
    delta=delta,
    s_pm1=s_pm1,
    log_s_ratio=log_s_ratio,
    omega_pm1=omega_pm1,
)
print("\nDone.")
