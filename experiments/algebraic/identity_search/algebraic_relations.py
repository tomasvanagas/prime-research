#!/usr/bin/env python3
"""
Algebraic relations search for f(x) = pi(x) - R(x).

Tests:
1. Bernoulli number relations
2. Zeta value relations via PSLQ
3. Dirichlet L-function value relations via PSLQ
4. Modular form (Ramanujan tau) correlations
5. Von Mangoldt / Chebyshev connection: f(x) vs g(x)/log(x)

Uses mpmath with mp.dps=50 for PSLQ.
"""

import numpy as np
import mpmath
from mpmath import mp, mpf, log as mlog, sqrt as msqrt, zeta as mzeta, bernoulli, pslq
from scipy.stats import pearsonr, spearmanr
import os, json, time

mp.dps = 50

# ── Load data ────────────────────────────────────────────────────────────────
DATA_DIR = os.path.dirname(os.path.abspath(__file__))
d = np.load(os.path.join(DATA_DIR, "fx_data.npz"))
x_arr = d["x"]       # 2..100000
f_arr = d["f"]        # float64, f(x) = pi(x) - R(x)
pi_arr = d["pi"]      # integer pi(x)

results = {}

# ══════════════════════════════════════════════════════════════════════════════
# 1. BERNOULLI NUMBER RELATIONS
# ══════════════════════════════════════════════════════════════════════════════
print("=" * 70)
print("1. BERNOULLI NUMBER RELATIONS")
print("=" * 70)

# Test: f(x) ~ sum_{k=1}^{K} B_{2k} * x^{-2k}
# This comes from the connection zeta(-n) = -B_{n+1}/(n+1).
# The oscillatory part of pi(x) relates to zeta zeros, but the Bernoulli
# correction terms come from the trivial zeros of zeta.

sample_xs = [100, 500, 1000, 5000, 10000, 50000, 100000]
K_values = [1, 2, 3, 5, 10, 20]

bernoulli_results = {}
for K in K_values:
    # Precompute Bernoulli numbers B_{2k} for k=1..K
    B_coeffs = [float(bernoulli(2 * k)) for k in range(1, K + 1)]

    corrections = []
    f_vals = []
    for xv in sample_xs:
        idx = xv - 2  # x_arr starts at 2
        fx = f_arr[idx]
        corr = sum(B_coeffs[k] * float(xv) ** (-2 * (k + 1)) for k in range(K))
        corrections.append(corr)
        f_vals.append(fx)

    # Also compute full-range correlation
    # Use subset of x values for speed (every 100th)
    subset = np.arange(0, len(x_arr), 100)
    f_sub = f_arr[subset]
    x_sub = x_arr[subset].astype(float)
    corr_sub = np.zeros(len(subset))
    for k in range(K):
        corr_sub += B_coeffs[k] * x_sub ** (-2 * (k + 1))

    r_val, p_val = pearsonr(f_sub, corr_sub)

    bernoulli_results[K] = {
        "pearson_r": float(r_val),
        "p_value": float(p_val),
        "sample_corrections": {str(xv): float(c) for xv, c in zip(sample_xs, corrections)},
        "sample_f_values": {str(xv): float(fv) for xv, fv in zip(sample_xs, f_vals)},
    }
    print(f"  K={K:3d}: Pearson r = {r_val:+.6f}, p = {p_val:.2e}")
    print(f"         Sample: f(1000)={f_vals[2]:.4f}, corr(1000)={corrections[2]:.2e}")

results["bernoulli"] = bernoulli_results

# Also test: f(x) * x^{1/2} vs Bernoulli sum (scaling might help)
print("\n  Scaled test: f(x)*sqrt(x) vs Bernoulli sums:")
f_scaled = f_sub * np.sqrt(x_sub)
for K in [1, 5, 10]:
    B_coeffs = [float(bernoulli(2 * k)) for k in range(1, K + 1)]
    corr_sub = np.zeros(len(subset))
    for k in range(K):
        corr_sub += B_coeffs[k] * x_sub ** (-2 * (k + 1))
    corr_scaled = corr_sub * np.sqrt(x_sub)
    r_val, p_val = pearsonr(f_scaled, corr_scaled)
    print(f"  K={K:3d}: Pearson r = {r_val:+.6f}, p = {p_val:.2e}")

# ══════════════════════════════════════════════════════════════════════════════
# 2. ZETA VALUE RELATIONS (PSLQ)
# ══════════════════════════════════════════════════════════════════════════════
print("\n" + "=" * 70)
print("2. ZETA VALUE RELATIONS (PSLQ)")
print("=" * 70)

# Precompute zeta values at high precision
zeta_vals = {}
for s in range(2, 8):
    zeta_vals[s] = mzeta(s)
    print(f"  zeta({s}) = {zeta_vals[s]}")

pslq_sample_xs = [100, 500, 1000, 5000, 10000, 50000]
zeta_pslq_results = {}

for xv in pslq_sample_xs:
    idx = xv - 2
    fx = mpf(str(f_arr[idx]))
    lx = mlog(mpf(xv))
    sx = msqrt(mpf(xv))

    # Vector: [f(x), zeta(2), zeta(3), zeta(4), zeta(5), zeta(6), zeta(7),
    #           log(x)*zeta(2), log(x)*zeta(3), sqrt(x), 1]
    vec = [
        fx,
        zeta_vals[2],
        zeta_vals[3],
        zeta_vals[4],
        zeta_vals[5],
        zeta_vals[6],
        zeta_vals[7],
        lx * zeta_vals[2],
        lx * zeta_vals[3],
        sx,
        mpf(1),
    ]

    rel = pslq(vec, maxcoeff=10**6, maxsteps=5000)
    label_names = [
        "f(x)", "zeta(2)", "zeta(3)", "zeta(4)", "zeta(5)",
        "zeta(6)", "zeta(7)", "log(x)*zeta(2)", "log(x)*zeta(3)",
        "sqrt(x)", "1"
    ]

    if rel is not None:
        # Format the relation
        terms = []
        for coeff, name in zip(rel, label_names):
            if coeff != 0:
                terms.append(f"{coeff}*{name}")
        relation_str = " + ".join(terms) + " = 0"
        print(f"  x={xv:6d}: RELATION FOUND: {relation_str}")
        zeta_pslq_results[str(xv)] = {"relation": rel, "formatted": relation_str}
    else:
        print(f"  x={xv:6d}: No relation found (coeffs < 10^6)")
        zeta_pslq_results[str(xv)] = {"relation": None}

results["zeta_pslq"] = zeta_pslq_results

# ══════════════════════════════════════════════════════════════════════════════
# 3. DIRICHLET L-FUNCTION VALUES (PSLQ)
# ══════════════════════════════════════════════════════════════════════════════
print("\n" + "=" * 70)
print("3. DIRICHLET L-FUNCTION VALUES (PSLQ)")
print("=" * 70)

def dirichlet_L(s_val, chi_func, mod_q, num_terms=50000):
    """Compute L(s, chi) = sum_{n=1}^{N} chi(n)/n^s via mpmath."""
    s_val = mpf(s_val)
    total = mpf(0)
    for n in range(1, num_terms + 1):
        c = chi_func(n, mod_q)
        if c != 0:
            total += mpf(c) / mpf(n) ** s_val
    return total

# Define characters
def chi_4_nontrivial(n, q):
    """Non-principal character mod 4: chi(1)=1, chi(3)=-1, chi(0)=chi(2)=0"""
    r = n % 4
    if r == 1: return 1
    if r == 3: return -1
    return 0

def chi_3_nontrivial(n, q):
    """Non-principal character mod 3: chi(1)=1, chi(2)=-1, chi(0)=0"""
    r = n % 3
    if r == 1: return 1
    if r == 2: return -1
    return 0

def chi_5_2(n, q):
    """A non-principal real character mod 5"""
    # Legendre symbol (n/5)
    r = n % 5
    if r == 0: return 0
    if r in [1, 4]: return 1
    return -1

def chi_7_quadratic(n, q):
    """Quadratic character mod 7 = Legendre (n/7)"""
    r = n % 7
    if r == 0: return 0
    if r in [1, 2, 4]: return 1
    return -1

characters = [
    ("chi_3", chi_3_nontrivial, 3),
    ("chi_4", chi_4_nontrivial, 4),
    ("chi_5", chi_5_2, 5),
    ("chi_7", chi_7_quadratic, 7),
]

print("  Computing L(1, chi) values...")
L_values = {}
for name, chi, q in characters:
    Lv = dirichlet_L(1, chi, q)
    L_values[name] = Lv
    print(f"    L(1, {name}) = {Lv}")

# PSLQ with L-values
dirichlet_pslq_results = {}
for xv in [100, 1000, 10000, 50000]:
    idx = xv - 2
    fx = mpf(str(f_arr[idx]))
    lx = mlog(mpf(xv))
    sx = msqrt(mpf(xv))

    vec = [
        fx,
        L_values["chi_3"],
        L_values["chi_4"],
        L_values["chi_5"],
        L_values["chi_7"],
        lx,
        sx,
        mpf(1),
    ]

    rel = pslq(vec, maxcoeff=10**6, maxsteps=5000)
    label_names = [
        "f(x)", "L(1,chi_3)", "L(1,chi_4)", "L(1,chi_5)",
        "L(1,chi_7)", "log(x)", "sqrt(x)", "1"
    ]

    if rel is not None:
        terms = []
        for coeff, name in zip(rel, label_names):
            if coeff != 0:
                terms.append(f"{coeff}*{name}")
        relation_str = " + ".join(terms) + " = 0"
        print(f"  x={xv:6d}: RELATION FOUND: {relation_str}")
        dirichlet_pslq_results[str(xv)] = {"relation": rel, "formatted": relation_str}
    else:
        print(f"  x={xv:6d}: No relation found")
        dirichlet_pslq_results[str(xv)] = {"relation": None}

results["dirichlet_pslq"] = dirichlet_pslq_results

# ══════════════════════════════════════════════════════════════════════════════
# 4. MODULAR FORM COEFFICIENTS (RAMANUJAN TAU)
# ══════════════════════════════════════════════════════════════════════════════
print("\n" + "=" * 70)
print("4. MODULAR FORM COEFFICIENTS (RAMANUJAN TAU)")
print("=" * 70)

def ramanujan_tau(N_max):
    """Compute tau(n) for n=1..N_max using the product formula for Delta(q).
    Delta(q) = q * prod_{n>=1} (1-q^n)^24 = sum tau(n) q^n.
    We work with power series truncated at N_max+1 terms.
    """
    # Start with coefficients of prod (1-q^n)^24
    # Build up iteratively
    coeffs = [0] * (N_max + 1)
    coeffs[0] = 1  # constant term of the product

    for n in range(1, N_max + 1):
        # Multiply current series by (1 - q^n)^24
        # Use the fact that (1-q^n)^24 has expansion via binomial
        # More efficient: compute prod (1-q^n)^24 step by step
        # Multiply by (1 - q^n) twenty-four times
        for _ in range(24):
            for j in range(N_max, n - 1, -1):
                coeffs[j] -= coeffs[j - n]

    # Delta(q) = q * prod, so tau(n) = coeffs[n-1]
    tau = {}
    for n in range(1, N_max + 1):
        tau[n] = coeffs[n - 1]
    return tau

print("  Computing Ramanujan tau(n) for n=1..100...")
tau = ramanujan_tau(100)
print(f"  tau(1)={tau[1]}, tau(2)={tau[2]}, tau(3)={tau[3]}, tau(4)={tau[4]}")
print(f"  tau(5)={tau[5]}, tau(10)={tau[10]}, tau(11)={tau[11]}")

# Correlation: f(n) vs tau(n) for n=2..100
n_range = range(2, 101)
f_vals = [f_arr[n - 2] for n in n_range]
tau_vals = [tau[n] for n in n_range]

r_direct, p_direct = pearsonr(f_vals, tau_vals)
rho_direct, p_rho = spearmanr(f_vals, tau_vals)
print(f"\n  f(n) vs tau(n), n=2..100:")
print(f"    Pearson  r = {r_direct:+.6f}, p = {p_direct:.2e}")
print(f"    Spearman r = {rho_direct:+.6f}, p = {p_rho:.2e}")

# Correlation: f(n) vs sum_{d|n} tau(d)
def divisor_sum_tau(n, tau_dict):
    s = 0
    for d in range(1, n + 1):
        if n % d == 0 and d in tau_dict:
            s += tau_dict[d]
    return s

tau_div_sums = [divisor_sum_tau(n, tau) for n in n_range]
r_div, p_div = pearsonr(f_vals, tau_div_sums)
rho_div, p_rho_div = spearmanr(f_vals, tau_div_sums)
print(f"\n  f(n) vs sum_{{d|n}} tau(d), n=2..100:")
print(f"    Pearson  r = {r_div:+.6f}, p = {p_div:.2e}")
print(f"    Spearman r = {rho_div:+.6f}, p = {p_rho_div:.2e}")

# Normalized version: tau(n)/n^{11/2} (Deligne bound says |tau(n)| < d(n)*n^{11/2})
tau_norm = [tau[n] / n ** 5.5 for n in n_range]
r_norm, p_norm = pearsonr(f_vals, tau_norm)
print(f"\n  f(n) vs tau(n)/n^{{11/2}}, n=2..100:")
print(f"    Pearson  r = {r_norm:+.6f}, p = {p_norm:.2e}")

results["ramanujan_tau"] = {
    "tau_samples": {str(n): tau[n] for n in [1, 2, 3, 4, 5, 10, 11, 12]},
    "f_vs_tau": {"pearson_r": float(r_direct), "p": float(p_direct),
                 "spearman_r": float(rho_direct)},
    "f_vs_tau_divisor_sum": {"pearson_r": float(r_div), "p": float(p_div),
                             "spearman_r": float(rho_div)},
    "f_vs_tau_normalized": {"pearson_r": float(r_norm), "p": float(p_norm)},
}

# ══════════════════════════════════════════════════════════════════════════════
# 5. VON MANGOLDT / CHEBYSHEV CONNECTION
# ══════════════════════════════════════════════════════════════════════════════
print("\n" + "=" * 70)
print("5. VON MANGOLDT / CHEBYSHEV CONNECTION")
print("=" * 70)

# Compute psi(x) = sum_{n<=x} Lambda(n) where Lambda(n) = log(p) if n=p^k
# First build a sieve for primes up to 100000
print("  Building prime sieve up to 100000...")
sieve_limit = 100001
is_prime = np.ones(sieve_limit, dtype=bool)
is_prime[0] = is_prime[1] = False
for i in range(2, int(sieve_limit ** 0.5) + 1):
    if is_prime[i]:
        is_prime[i * i::i] = False

# Compute Lambda(n) for n=1..100000
print("  Computing von Mangoldt function Lambda(n)...")
Lambda = np.zeros(sieve_limit, dtype=np.float64)
primes = np.where(is_prime)[0]
for p in primes:
    pk = p
    logp = np.log(p)
    while pk < sieve_limit:
        Lambda[pk] = logp
        pk *= p

# Compute psi(x) = cumulative sum of Lambda
psi = np.cumsum(Lambda)  # psi[n] = sum_{k=1}^{n} Lambda(k)

# g(x) = psi(x) - x
g_arr = psi[2:sieve_limit] - np.arange(2, sieve_limit, dtype=np.float64)

# f(x) and g(x)/log(x)
x_float = x_arr.astype(np.float64)
g_over_logx = g_arr[:len(x_arr)] / np.log(x_float)

# Residual: f(x) - g(x)/log(x)
residual = f_arr - g_over_logx

# Statistics
print(f"\n  f(x) statistics:   mean={np.mean(f_arr):.4f}, std={np.std(f_arr):.4f}")
print(f"  g(x)/log(x) stats: mean={np.mean(g_over_logx):.4f}, std={np.std(g_over_logx):.4f}")
print(f"  Residual stats:    mean={np.mean(residual):.4f}, std={np.std(residual):.4f}")
print(f"  |f(x)| mean:      {np.mean(np.abs(f_arr)):.4f}")
print(f"  |residual| mean:  {np.mean(np.abs(residual)):.4f}")

# Is the residual simpler (smaller)?
improvement = 1.0 - np.std(residual) / np.std(f_arr)
print(f"\n  Std reduction: {improvement*100:.2f}%")
print(f"  {'YES' if improvement > 0.01 else 'NO'}: residual {'is' if improvement > 0.01 else 'is NOT'} simpler than f(x)")

# Correlation between f and g/log
r_fg, p_fg = pearsonr(f_arr, g_over_logx)
print(f"\n  Pearson corr f(x) vs g(x)/log(x): r = {r_fg:+.6f}, p = {p_fg:.2e}")

# Check at sample points
print("\n  Sample values:")
print(f"  {'x':>8s} {'f(x)':>12s} {'g(x)/log(x)':>14s} {'residual':>12s}")
for xv in sample_xs:
    idx = xv - 2
    if idx < len(f_arr):
        print(f"  {xv:8d} {f_arr[idx]:12.4f} {g_over_logx[idx]:14.4f} {residual[idx]:12.4f}")

# Test if residual correlates with 1/log(x), Li(x) correction terms
corr_invlog, p_invlog = pearsonr(residual, 1.0 / np.log(x_float))
print(f"\n  Residual vs 1/log(x): r = {corr_invlog:+.6f}, p = {p_invlog:.2e}")

# Check variance ratio at different scales
print("\n  Scale analysis (std ratio |residual|/|f| in windows):")
for lo, hi in [(2, 100), (100, 1000), (1000, 10000), (10000, 100000)]:
    sl = slice(lo - 2, hi - 2)
    sf = np.std(f_arr[sl])
    sr = np.std(residual[sl])
    ratio = sr / sf if sf > 0 else float('inf')
    print(f"    [{lo:6d}, {hi:6d}): std(f)={sf:.4f}, std(res)={sr:.4f}, ratio={ratio:.4f}")

results["chebyshev"] = {
    "f_std": float(np.std(f_arr)),
    "residual_std": float(np.std(residual)),
    "std_reduction_pct": float(improvement * 100),
    "pearson_f_vs_g_over_logx": float(r_fg),
    "residual_simpler": bool(improvement > 0.01),
    "sample_values": {
        str(xv): {
            "f": float(f_arr[xv - 2]),
            "g_over_logx": float(g_over_logx[xv - 2]),
            "residual": float(residual[xv - 2]),
        }
        for xv in sample_xs if xv - 2 < len(f_arr)
    },
}

# ══════════════════════════════════════════════════════════════════════════════
# SUMMARY
# ══════════════════════════════════════════════════════════════════════════════
print("\n" + "=" * 70)
print("SUMMARY")
print("=" * 70)

print("""
1. BERNOULLI: The Bernoulli correction terms sum_k B_{2k} x^{-2k} are tiny
   (order 10^{-6} to 10^{-70}) while f(x) is O(1)-O(10). No correlation.

2. ZETA PSLQ: Searched for integer relations among f(x) and zeta values.
   Any relations found are x-dependent (hence likely spurious/coincidental).

3. DIRICHLET PSLQ: Searched for integer relations with L(1, chi) values.
   Again x-dependent if found.

4. RAMANUJAN TAU: Correlation between f(n) and tau(n) for small n.
   Expected to be negligible (tau grows as n^{11/2}, f is O(sqrt(n)/log(n))).

5. CHEBYSHEV: f(x) - g(x)/log(x) where g(x) = psi(x) - x.
   This is the key partial-summation test. If the residual variance is
   significantly smaller, it means the Chebyshev residual captures a
   substantial part of f(x).
""")

# Save results
out_path = os.path.join(DATA_DIR, "algebraic_relations_results.json")
# Convert non-serializable types
def make_serializable(obj):
    if isinstance(obj, (np.integer,)):
        return int(obj)
    if isinstance(obj, (np.floating,)):
        return float(obj)
    if isinstance(obj, np.ndarray):
        return obj.tolist()
    if isinstance(obj, dict):
        return {k: make_serializable(v) for k, v in obj.items()}
    if isinstance(obj, list):
        return [make_serializable(i) for i in obj]
    return obj

with open(out_path, "w") as fp:
    json.dump(make_serializable(results), fp, indent=2, default=str)
print(f"\nResults saved to {out_path}")
print("Done.")
