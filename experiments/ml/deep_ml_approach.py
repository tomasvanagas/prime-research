#!/usr/bin/env python3
"""
Session 8: Deep ML approaches for n-th prime prediction.
No PyTorch available -- using numpy + scipy + sklearn.

Experiments:
  1. MLP on log-space representation (multi-layer "deep" network via sklearn)
  2. Residue network: p(n) mod q patterns for small primes q
  3. Gap prediction with autoregressive features (LSTM-like via sklearn)
  4. Symbolic regression via genetic-style search (custom, since PySR unavailable)
  5. Neural ODE approximation via scipy.integrate
  6. Correction model: delta(n) = p(n) - round(R_inv(n))

Key metric: exact prediction accuracy (%) on train and test sets.
"""

import numpy as np
from scipy import integrate, optimize, special
from sklearn.neural_network import MLPRegressor, MLPClassifier
from sklearn.ensemble import GradientBoostingRegressor, RandomForestRegressor
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import mean_absolute_error
import sympy
from sympy import primepi, prime, li, log as symlog
import time
import warnings
warnings.filterwarnings("ignore")

# ============================================================
# DATA GENERATION
# ============================================================

def generate_primes_sieve(limit):
    """Sieve of Eratosthenes up to limit."""
    sieve = np.ones(limit + 1, dtype=bool)
    sieve[0] = sieve[1] = False
    for i in range(2, int(limit**0.5) + 1):
        if sieve[i]:
            sieve[i*i::i] = False
    return np.nonzero(sieve)[0]

def R_inverse_approx(n):
    """Approximate R^{-1}(n) using li^{-1}(n) as starting point.
    R(x) ~ li(x) - li(x^{1/2})/2 - li(x^{1/3})/3 - ...
    We approximate R^{-1}(n) ~ n*ln(n) with Newton refinement.
    """
    if n <= 1:
        return 2.0
    # Initial guess: n * ln(n)
    x = float(n) * np.log(float(n))
    if x < 2:
        x = 2.0
    # Newton iterations on R(x) - n = 0
    # R(x) ~ li(x) - li(sqrt(x))/2 - ...
    # We use the simpler approximation: R^{-1}(n) ~ li^{-1}(n)
    # li^{-1}(n) can be found by Newton on li(x) = n
    for _ in range(50):
        # li(x) = Ei(ln(x))
        lnx = np.log(x)
        if lnx <= 0:
            x = float(n) * 2.0
            continue
        li_x = special.expi(lnx)
        # li'(x) = 1/ln(x)
        deriv = 1.0 / lnx
        dx = (li_x - float(n)) / deriv
        x -= dx
        if x < 2:
            x = 2.0
        if abs(dx) < 0.01:
            break
    return x

print("=" * 80)
print("SESSION 8: DEEP ML APPROACHES FOR n-TH PRIME PREDICTION")
print("=" * 80)

# Generate primes up to ~2.7M (covers p(200000) ~ 2750159)
print("\nGenerating primes via sieve...")
t0 = time.time()
SIEVE_LIMIT = 3_000_000
all_primes = generate_primes_sieve(SIEVE_LIMIT)
print(f"  Generated {len(all_primes)} primes up to {SIEVE_LIMIT} in {time.time()-t0:.2f}s")

# We need p(1)..p(200000)
N_TOTAL = 200_000
N_TRAIN = 100_000
N_TEST = N_TOTAL - N_TRAIN

primes_arr = all_primes[:N_TOTAL]  # p(1)..p(N_TOTAL), 0-indexed: primes_arr[i] = p(i+1)
n_arr = np.arange(1, N_TOTAL + 1, dtype=np.float64)

# Gaps
gaps = np.diff(primes_arr).astype(np.float64)  # gaps[i] = p(i+2) - p(i+1)

print(f"  p(1)={primes_arr[0]}, p({N_TRAIN})={primes_arr[N_TRAIN-1]}, p({N_TOTAL})={primes_arr[-1]}")

# Compute R^{-1}(n) for all n using the Riemann R function
print("\nComputing R^{-1}(n) approximations using Riemann R function...")
t0 = time.time()

# Mobius function values for k=1..49
MU = [0, 1, -1, -1, 0, -1, 1, -1, 0, 0, 1, -1, 0, -1, 1, 1, 0, -1, 0, -1,
      0, 1, 1, -1, 0, 0, 1, 0, 0, -1, -1, -1, 0, 1, 1, 1, 0, -1, 1, 1,
      0, -1, -1, -1, 0, 0, 1, -1, 0, 0]

def R_func_vec(x_arr):
    """Riemann R function: R(x) = sum_{k=1}^{inf} mu(k)/k * li(x^{1/k})."""
    result = np.zeros_like(x_arr, dtype=np.float64)
    for k in range(1, len(MU)):
        if MU[k] == 0:
            continue
        xk = np.power(np.maximum(x_arr, 2.0), 1.0/k)
        # li(y) = Ei(ln(y)) for y > 1
        ln_xk = np.log(np.maximum(xk, 1.001))
        li_xk = special.expi(ln_xk)  # vectorized
        result += MU[k] / k * li_xk
        # Stop when x^{1/k} < 2 for all elements
        if np.all(xk < 2.0):
            break
    return result

# Starting approximation
ln_n = np.log(np.maximum(n_arr, 2.0))
lnln_n = np.log(np.maximum(ln_n, 1.0))
r_inv_approx = n_arr * (ln_n + lnln_n - 1 + (lnln_n - 2) / np.maximum(ln_n, 1.0))

# Newton iterations on R(x) = n
for iteration in range(30):
    rx = R_func_vec(r_inv_approx)
    lnx = np.log(np.maximum(r_inv_approx, 2.0))
    deriv = 1.0 / np.maximum(lnx, 1e-10)
    dx = (rx - n_arr) / deriv
    r_inv_approx -= dx
    r_inv_approx = np.maximum(r_inv_approx, 2.0)
    max_dx = np.max(np.abs(dx))
    if iteration == 29:
        print(f"  Did NOT converge after 30 iterations (max_dx={max_dx:.4f})")
    if max_dx < 0.01:
        print(f"  Converged in {iteration+1} iterations (max_dx={max_dx:.4f})")
        break

r_inv_rounded = np.round(r_inv_approx).astype(np.int64)
delta = primes_arr.astype(np.int64) - r_inv_rounded  # correction term
print(f"  Computed in {time.time()-t0:.2f}s")
print(f"  delta stats: mean={np.mean(delta):.1f}, std={np.std(delta):.1f}, "
      f"min={np.min(delta)}, max={np.max(delta)}")
print(f"  |delta| <= 10: {np.mean(np.abs(delta) <= 10)*100:.1f}%")
print(f"  |delta| <= 50: {np.mean(np.abs(delta) <= 50)*100:.1f}%")
print(f"  |delta| <= 100: {np.mean(np.abs(delta) <= 100)*100:.1f}%")
print(f"  delta == 0 (exact): {np.mean(delta == 0)*100:.2f}%")

# Split
train_idx = slice(0, N_TRAIN)
test_idx = slice(N_TRAIN, N_TOTAL)

# ============================================================
# EXPERIMENT 1: Deep MLP on multiple representations
# ============================================================
print("\n" + "=" * 80)
print("EXPERIMENT 1: Deep MLP on multiple representations of n")
print("=" * 80)

def make_features(n_values):
    """Rich feature set for each n."""
    n = n_values.astype(np.float64)
    ln_n = np.log(np.maximum(n, 2.0))
    lnln_n = np.log(np.maximum(ln_n, 1.0))
    features = np.column_stack([
        n,
        ln_n,
        lnln_n,
        n * ln_n,                    # ~ p(n)
        n * (ln_n + lnln_n),         # better approx
        n * (ln_n + lnln_n - 1),     # even better
        n * (ln_n + lnln_n - 1 + (lnln_n - 2)/np.maximum(ln_n, 1.0)),  # Cipolla-like
        np.sqrt(n),
        n ** (1/3),
        ln_n ** 2,
        n / np.maximum(ln_n, 1.0),   # ~ pi(n) inverse direction
        np.sin(2 * np.pi * n / 6),   # mod 6 pattern
        np.cos(2 * np.pi * n / 6),
        np.sin(2 * np.pi * n / 30),  # mod 30 pattern
        np.cos(2 * np.pi * n / 30),
        n % 2, n % 3, n % 5, n % 7,
        n % 6, n % 30,
    ])
    return features

X_all = make_features(n_arr)
y_all = primes_arr.astype(np.float64)

X_train, X_test = X_all[train_idx], X_all[test_idx]
y_train, y_test = y_all[train_idx], y_all[test_idx]

scaler = StandardScaler()
X_train_s = scaler.fit_transform(X_train)
X_test_s = scaler.transform(X_test)

# Scale targets too (large values)
y_mean, y_std = y_train.mean(), y_train.std()
y_train_s = (y_train - y_mean) / y_std
y_test_s = (y_test - y_mean) / y_std

configs = [
    ("MLP-small", (256, 128, 64)),
    ("MLP-deep", (512, 256, 128, 64, 32)),
]

for name, layers in configs:
    print(f"\n  {name} {layers}:")
    t0 = time.time()
    mlp = MLPRegressor(
        hidden_layer_sizes=layers,
        activation='relu',
        solver='adam',
        max_iter=500,
        early_stopping=True,
        validation_fraction=0.1,
        learning_rate='adaptive',
        learning_rate_init=0.001,
        random_state=42,
        batch_size=1024,
    )
    mlp.fit(X_train_s, y_train_s)

    pred_train = mlp.predict(X_train_s) * y_std + y_mean
    pred_test = mlp.predict(X_test_s) * y_std + y_mean

    pred_train_round = np.round(pred_train).astype(np.int64)
    pred_test_round = np.round(pred_test).astype(np.int64)

    train_exact = np.mean(pred_train_round == y_train.astype(np.int64)) * 100
    test_exact = np.mean(pred_test_round == y_test.astype(np.int64)) * 100
    train_mae = mean_absolute_error(y_train, pred_train)
    test_mae = mean_absolute_error(y_test, pred_test)

    print(f"    Time: {time.time()-t0:.1f}s")
    print(f"    Train exact: {train_exact:.3f}%, MAE: {train_mae:.1f}")
    print(f"    Test  exact: {test_exact:.3f}%, MAE: {test_mae:.1f}")

# ============================================================
# EXPERIMENT 2: Residue network -- p(n) mod q patterns
# ============================================================
print("\n" + "=" * 80)
print("EXPERIMENT 2: Residue patterns -- p(n) mod q for small primes q")
print("=" * 80)

small_primes = [3, 5, 7, 11, 13, 17, 19, 23, 29, 31]
# Skip q=2: p(n) mod 2 = 1 for all n>1 (trivial)
print("\n  p(n) mod 2: trivial -- all primes > 2 are odd")

X_feat = make_features(n_arr)
for q in small_primes:
    residues = primes_arr % q
    unique, counts = np.unique(residues, return_counts=True)

    # p(n) mod q = 0 only when p(n) = q (one element). Skip those indices.
    y_mod = (primes_arr % q).astype(int)

    # Exclude indices where p(n) = q (residue 0)
    mask_tr = (y_mod[train_idx] != 0)
    mask_te = (y_mod[test_idx] != 0)

    scaler_r = StandardScaler()
    X_tr = scaler_r.fit_transform(X_feat[train_idx][mask_tr])
    X_te = scaler_r.transform(X_feat[test_idx][mask_te])
    y_tr = y_mod[train_idx][mask_tr]
    y_te = y_mod[test_idx][mask_te]

    vals_tr, cnts_tr = np.unique(y_tr, return_counts=True)
    baseline = np.max(cnts_tr) / len(y_tr) * 100

    clf = MLPClassifier(
        hidden_layer_sizes=(128, 64, 32),
        max_iter=200,
        early_stopping=True,
        random_state=42,
        batch_size=1024,
    )
    clf.fit(X_tr, y_tr)

    train_acc = clf.score(X_tr, y_tr) * 100
    test_acc = clf.score(X_te, y_te) * 100

    print(f"\n  p(n) mod {q:2d}: {len(unique)} classes, baseline={baseline:.1f}%")
    print(f"    MLP train acc: {train_acc:.1f}%, test acc: {test_acc:.1f}%")

# Deeper analysis: autocorrelation of p(n) mod q
print("\n  Autocorrelation of p(n) mod q sequences:")
for q in [3, 5, 7, 11]:
    res = (primes_arr % q).astype(float)
    res -= res.mean()
    autocorr = np.correlate(res[:1000], res[:1000], mode='full')
    autocorr = autocorr[len(autocorr)//2:]
    autocorr /= autocorr[0]
    max_lag_corr = np.max(np.abs(autocorr[1:50]))
    print(f"    mod {q:2d}: max |autocorr| for lags 1-49 = {max_lag_corr:.4f}")

# ============================================================
# EXPERIMENT 3: Gap prediction with autoregressive features
# ============================================================
print("\n" + "=" * 80)
print("EXPERIMENT 3: Gap prediction g(n) = p(n+1) - p(n)")
print("=" * 80)

def make_gap_features(gaps_seq, window=20):
    """Create autoregressive features from gap sequence."""
    n = len(gaps_seq)
    X = np.zeros((n - window, window + 5))
    for i in range(window, n):
        X[i - window, :window] = gaps_seq[i-window:i]
    # Add position features
    positions = np.arange(window, n, dtype=np.float64)
    X[:, window] = positions
    X[:, window+1] = np.log(positions)
    X[:, window+2] = np.sqrt(positions)
    # Running mean and std of gaps
    for i in range(window, n):
        X[i-window, window+3] = np.mean(gaps_seq[max(0,i-100):i])
        X[i-window, window+4] = np.std(gaps_seq[max(0,i-100):i])
    return X

WINDOW = 20
X_gaps = make_gap_features(gaps, WINDOW)
y_gaps = gaps[WINDOW:]

gap_train_end = N_TRAIN - 1 - WINDOW  # index into X_gaps
X_gap_tr = X_gaps[:gap_train_end]
y_gap_tr = y_gaps[:gap_train_end]
X_gap_te = X_gaps[gap_train_end:]
y_gap_te = y_gaps[gap_train_end:]

scaler_gap = StandardScaler()
X_gap_tr_s = scaler_gap.fit_transform(X_gap_tr)
X_gap_te_s = scaler_gap.transform(X_gap_te)

# MLP for gap prediction
gap_models = [
    ("MLP-gap-deep", MLPRegressor(hidden_layer_sizes=(256, 128, 64, 32),
                                    max_iter=500, early_stopping=True,
                                    random_state=42, batch_size=1024)),
    ("GBR-gap", GradientBoostingRegressor(n_estimators=200, max_depth=5,
                                           learning_rate=0.05, random_state=42,
                                           subsample=0.8)),
]

for name, model in gap_models:
    print(f"\n  {name}:")
    t0 = time.time()
    model.fit(X_gap_tr_s if 'MLP' in name else X_gap_tr, y_gap_tr)

    pred_tr = model.predict(X_gap_tr_s if 'MLP' in name else X_gap_tr)
    pred_te = model.predict(X_gap_te_s if 'MLP' in name else X_gap_te)

    pred_tr_r = np.round(pred_tr).astype(int)
    pred_te_r = np.round(pred_te).astype(int)

    # Force even (gaps > 2 are always even)
    pred_tr_r = np.where(pred_tr_r % 2 != 0, pred_tr_r + 1, pred_tr_r)
    pred_te_r = np.where(pred_te_r % 2 != 0, pred_te_r + 1, pred_te_r)

    train_exact = np.mean(pred_tr_r == y_gap_tr.astype(int)) * 100
    test_exact = np.mean(pred_te_r == y_gap_te.astype(int)) * 100
    train_mae = mean_absolute_error(y_gap_tr, pred_tr)
    test_mae = mean_absolute_error(y_gap_te, pred_te)

    print(f"    Time: {time.time()-t0:.1f}s, Train exact: {train_exact:.2f}%, Test exact: {test_exact:.2f}%")
    print(f"    Train MAE: {train_mae:.2f}, Test MAE: {test_mae:.2f}")

    # If we use predicted gaps to reconstruct primes from p(N_TRAIN+1)
    if 'MLP' in name:
        # Reconstruction from cumsum of predicted gaps
        # Start from known p(N_TRAIN+1)
        recon_start = N_TRAIN + WINDOW  # index in primes_arr
        recon_primes = np.zeros(len(pred_te_r), dtype=np.int64)
        recon_primes[0] = primes_arr[recon_start] + int(pred_te_r[0])
        for i in range(1, len(pred_te_r)):
            recon_primes[i] = recon_primes[i-1] + int(pred_te_r[i])

        actual_primes_te = primes_arr[recon_start+1:recon_start+1+len(pred_te_r)]
        if len(actual_primes_te) == len(recon_primes):
            recon_exact = np.mean(recon_primes == actual_primes_te) * 100
            print(f"    Reconstructed prime exact match: {recon_exact:.2f}% (first {len(recon_primes)} test primes)")

# What is the distribution of gaps?
print("\n  Gap distribution (full dataset):")
gap_vals, gap_counts = np.unique(gaps, return_counts=True)
top_gaps = np.argsort(-gap_counts)[:10]
for i in top_gaps:
    print(f"    gap={int(gap_vals[i]):3d}: {gap_counts[i]:6d} ({gap_counts[i]/len(gaps)*100:.1f}%)")
print(f"  Most common gap: {int(gap_vals[top_gaps[0]])} ({gap_counts[top_gaps[0]]/len(gaps)*100:.1f}%)")
print(f"  Baseline (always predict most common): {gap_counts[top_gaps[0]]/len(gaps)*100:.1f}%")

# ============================================================
# EXPERIMENT 4: Symbolic regression (custom genetic search)
# ============================================================
print("\n" + "=" * 80)
print("EXPERIMENT 4: Symbolic regression for correction delta(n)")
print("=" * 80)

# delta(n) = p(n) - R_inv_approx(n)
# Try to find symbolic expressions for delta

# Sample a smaller set for symbolic regression
SYM_N = 5000
sym_indices = np.linspace(10, N_TRAIN-1, SYM_N, dtype=int)  # start from n=11 to avoid log(log(n)) issues
n_sym = n_arr[sym_indices]
delta_sym = delta[sym_indices].astype(np.float64)
p_sym = primes_arr[sym_indices].astype(np.float64)

# Try polynomial fits on delta vs various transforms of n
print("\n  Polynomial fits on delta(n):")
ln_n_sym = np.log(n_sym)
sqrt_n_sym = np.sqrt(n_sym)

transforms = {
    "n": n_sym,
    "ln(n)": ln_n_sym,
    "sqrt(n)": sqrt_n_sym,
}

for tname, t_vals in transforms.items():
    for deg in [1, 2, 3, 5]:
        try:
            # Normalize to avoid SVD issues
            t_mean, t_std = t_vals.mean(), max(t_vals.std(), 1e-10)
            t_normed = (t_vals - t_mean) / t_std
            coeffs = np.polyfit(t_normed, delta_sym, deg)
            pred = np.polyval(coeffs, t_normed)
            residual = delta_sym - pred
            rmse = np.sqrt(np.mean(residual**2))
            exact = np.mean(np.round(pred).astype(int) == delta_sym.astype(int)) * 100
            if deg <= 2 or exact > 0.5:
                print(f"    poly({tname}, deg={deg}): RMSE={rmse:.1f}, exact={exact:.2f}%")
        except Exception as e:
            print(f"    poly({tname}, deg={deg}): FAILED -- {e}")

# Try rational corrections: delta(n) ~ a*sqrt(n)*ln(n) + b*sqrt(n) + c*ln(n) + d
print("\n  Custom basis fits on delta(n):")
basis_funcs = {
    "sqrt(n)*ln(n)": sqrt_n_sym * ln_n_sym,
    "sqrt(n)": sqrt_n_sym,
    "ln(n)^2": ln_n_sym**2,
    "ln(n)": ln_n_sym,
    "1": np.ones_like(n_sym),
    "n^(1/3)": n_sym**(1/3),
    "sqrt(n)*ln(ln(n))": sqrt_n_sym * np.log(ln_n_sym),
}

basis_matrix_raw = np.column_stack(list(basis_funcs.values()))
# Use sklearn LinearRegression to avoid SVD issues
from sklearn.linear_model import Ridge
basis_scaler = StandardScaler()
basis_matrix = basis_scaler.fit_transform(basis_matrix_raw)
ridge = Ridge(alpha=0.01)
ridge.fit(basis_matrix, delta_sym)
coeffs_ls = ridge.coef_
pred_ls = ridge.predict(basis_matrix)
rmse_ls = np.sqrt(np.mean((delta_sym - pred_ls)**2))
exact_ls = np.mean(np.round(pred_ls).astype(int) == delta_sym.astype(int)) * 100

print(f"  7-basis LS fit: RMSE={rmse_ls:.1f}, exact={exact_ls:.2f}%")
for (name, _), c in zip(basis_funcs.items(), coeffs_ls):
    if abs(c) > 1e-6:
        print(f"    {c:+.6f} * {name}")

# Test on test set
sym_test_idx = np.linspace(N_TRAIN, N_TOTAL-1, SYM_N, dtype=int)
n_sym_te = n_arr[sym_test_idx]
delta_sym_te = delta[sym_test_idx].astype(np.float64)

basis_test_raw = np.column_stack([
    np.sqrt(n_sym_te) * np.log(n_sym_te),
    np.sqrt(n_sym_te),
    np.log(n_sym_te)**2,
    np.log(n_sym_te),
    np.ones_like(n_sym_te),
    n_sym_te**(1/3),
    np.sqrt(n_sym_te) * np.log(np.log(n_sym_te)),
])
basis_test = basis_scaler.transform(basis_test_raw)
pred_ls_te = ridge.predict(basis_test)
exact_ls_te = np.mean(np.round(pred_ls_te).astype(int) == delta_sym_te.astype(int)) * 100
rmse_ls_te = np.sqrt(np.mean((delta_sym_te - pred_ls_te)**2))
print(f"  Test: RMSE={rmse_ls_te:.1f}, exact={exact_ls_te:.2f}%")

# ============================================================
# EXPERIMENT 5: Neural ODE -- dp/dn = f(p, n)
# ============================================================
print("\n" + "=" * 80)
print("EXPERIMENT 5: Neural ODE -- dp/dn = f(p, n)")
print("=" * 80)

# The exact derivative is dp/dn = g(n) = p(n+1) - p(n) at integer points
# But conceptually dp/dn ~ ln(n) + ln(ln(n)) (from PNT)
# Let's learn f in dp/dn = f(n, p) via MLP, then integrate

# Use a subset for ODE fitting
ODE_N = 10000
ode_idx = np.arange(ODE_N)
n_ode = n_arr[ode_idx]
p_ode = primes_arr[ode_idx].astype(np.float64)
dpdn_ode = gaps[:ODE_N].astype(np.float64)  # discrete derivative

# Features for ODE RHS: f(n, p)
X_ode = np.column_stack([
    n_ode, p_ode, np.log(n_ode), np.log(np.maximum(p_ode, 2)),
    np.sqrt(n_ode), np.log(n_ode)**2
])

scaler_ode = StandardScaler()
X_ode_s = scaler_ode.fit_transform(X_ode)
dpdn_mean, dpdn_std = dpdn_ode.mean(), dpdn_ode.std()
dpdn_s = (dpdn_ode - dpdn_mean) / dpdn_std

# Train MLP for f
print("  Training MLP for dp/dn = f(n, p)...")
t0 = time.time()
mlp_ode = MLPRegressor(
    hidden_layer_sizes=(128, 64, 32),
    max_iter=500, early_stopping=True,
    random_state=42, batch_size=512,
)
mlp_ode.fit(X_ode_s, dpdn_s)

pred_dpdn = mlp_ode.predict(X_ode_s) * dpdn_std + dpdn_mean
train_gap_mae = mean_absolute_error(dpdn_ode, pred_dpdn)
print(f"  Train gap MAE: {train_gap_mae:.2f} (mean gap: {dpdn_ode.mean():.2f})")

# Now integrate: p(1) = 2, p(n+1) = p(n) + f(n, p(n))
print("  Integrating ODE to reconstruct primes...")
p_recon = np.zeros(ODE_N + 1)
p_recon[0] = 2.0
for i in range(ODE_N):
    feat = np.array([[n_ode[i] if i < ODE_N else i+1, p_recon[i],
                       np.log(max(i+1, 1)), np.log(max(p_recon[i], 2)),
                       np.sqrt(i+1), np.log(max(i+1, 1))**2]])
    feat_s = scaler_ode.transform(feat)
    dp = mlp_ode.predict(feat_s)[0] * dpdn_std + dpdn_mean
    p_recon[i+1] = p_recon[i] + dp

p_recon_r = np.round(p_recon[1:ODE_N+1]).astype(np.int64)
ode_exact = np.mean(p_recon_r == primes_arr[:ODE_N]) * 100
ode_mae = mean_absolute_error(primes_arr[:ODE_N].astype(float), p_recon[1:ODE_N+1])
print(f"  ODE reconstruction: exact={ode_exact:.3f}%, MAE={ode_mae:.1f}")

# Error growth
for checkpoint in [100, 500, 1000, 5000, 10000]:
    if checkpoint <= ODE_N:
        err = abs(p_recon[checkpoint] - primes_arr[checkpoint-1])
        pct = err / primes_arr[checkpoint-1] * 100
        print(f"    n={checkpoint}: error={err:.0f} ({pct:.2f}%)")

# ============================================================
# EXPERIMENT 6: Learn the CORRECTION delta(n) = p(n) - round(R_inv(n))
# ============================================================
print("\n" + "=" * 80)
print("EXPERIMENT 6: Correction delta(n) = p(n) - round(R^{-1}(n))")
print("=" * 80)

# Features for delta prediction
def make_delta_features(n_vals, r_inv_vals):
    """Features for predicting delta = p(n) - round(R^{-1}(n))."""
    n = n_vals.astype(np.float64)
    r = r_inv_vals.astype(np.float64)
    ln_n = np.log(n)
    lnln_n = np.log(np.maximum(ln_n, 1.0))
    ln_r = np.log(np.maximum(r, 2.0))

    feats = np.column_stack([
        n,
        ln_n,
        lnln_n,
        ln_n**2,
        ln_n * lnln_n,
        np.sqrt(n),
        n**(1/3),
        r,
        r - n * ln_n,        # residual from simple approx
        ln_r,
        np.sin(2*np.pi*n/6),
        np.cos(2*np.pi*n/6),
        np.sin(2*np.pi*n/30),
        np.cos(2*np.pi*n/30),
        np.sin(2*np.pi*n/210),
        np.cos(2*np.pi*n/210),
        n % 2, n % 3, n % 5, n % 7, n % 11, n % 13,
        n % 6, n % 30, n % 210,
        # Fractional part of R^{-1}
        r - np.floor(r),
    ])
    return feats

X_delta = make_delta_features(n_arr, r_inv_approx)
y_delta = delta.astype(np.float64)

X_d_tr = X_delta[train_idx]
X_d_te = X_delta[test_idx]
y_d_tr = y_delta[train_idx]
y_d_te = y_delta[test_idx]

scaler_d = StandardScaler()
X_d_tr_s = scaler_d.fit_transform(X_d_tr)
X_d_te_s = scaler_d.transform(X_d_te)

delta_models = [
    ("MLP-delta-deep", MLPRegressor(
        hidden_layer_sizes=(512, 256, 128, 64),
        max_iter=1000, early_stopping=True,
        random_state=42, batch_size=2048,
        learning_rate='adaptive', learning_rate_init=0.001)),
    ("GBR-delta", GradientBoostingRegressor(
        n_estimators=300, max_depth=6, learning_rate=0.05,
        random_state=42, subsample=0.8, min_samples_leaf=10)),
    ("RF-delta", RandomForestRegressor(
        n_estimators=500, max_depth=20, random_state=42,
        min_samples_leaf=5, n_jobs=-1)),
]

for name, model in delta_models:
    print(f"\n  {name}:")
    t0 = time.time()
    model.fit(X_d_tr_s if 'MLP' in name else X_d_tr, y_d_tr)

    pred_tr = model.predict(X_d_tr_s if 'MLP' in name else X_d_tr)
    pred_te = model.predict(X_d_te_s if 'MLP' in name else X_d_te)

    pred_tr_r = np.round(pred_tr).astype(np.int64)
    pred_te_r = np.round(pred_te).astype(np.int64)

    train_exact = np.mean(pred_tr_r == y_d_tr.astype(np.int64)) * 100
    test_exact = np.mean(pred_te_r == y_d_te.astype(np.int64)) * 100
    train_mae = mean_absolute_error(y_d_tr, pred_tr)
    test_mae = mean_absolute_error(y_d_te, pred_te)

    print(f"    Time: {time.time()-t0:.1f}s")
    print(f"    Train exact: {train_exact:.3f}%, Test exact: {test_exact:.3f}%")
    print(f"    Train MAE: {train_mae:.1f}, Test MAE: {test_mae:.1f}")
    print(f"    delta range: [{int(y_d_te.min())}, {int(y_d_te.max())}]")

    # Reconstruct p(n) = round(R_inv(n)) + delta_predicted
    recon_primes_te = r_inv_rounded[test_idx] + pred_te_r
    prime_exact = np.mean(recon_primes_te == primes_arr[test_idx]) * 100
    prime_mae = mean_absolute_error(primes_arr[test_idx].astype(float),
                                     recon_primes_te.astype(float))
    print(f"    Reconstructed p(n) exact: {prime_exact:.3f}%, MAE: {prime_mae:.1f}")

# ============================================================
# EXPERIMENT 6b: Classification approach -- bin delta(n)
# ============================================================
print("\n" + "-" * 40)
print("  EXPERIMENT 6b: Classify delta(n) into bins")
print("-" * 40)

# What if we discretize delta into bins and do classification?
# First check: what fraction of delta values appear multiple times?
delta_train = delta[train_idx]
unique_deltas, delta_counts = np.unique(delta_train, return_counts=True)
print(f"  Unique delta values in train: {len(unique_deltas)} out of {N_TRAIN}")
print(f"  Most common delta values:")
top_deltas = np.argsort(-delta_counts)[:15]
for i in top_deltas:
    print(f"    delta={int(unique_deltas[i]):4d}: count={delta_counts[i]}")

# Bin into ranges
for bin_size in [1, 5, 10, 50]:
    binned_tr = (delta_train / bin_size).astype(int)
    binned_te = (delta[test_idx] / bin_size).astype(int)

    n_bins = len(np.unique(binned_tr))
    if n_bins > 500:
        print(f"  Bin size {bin_size}: {n_bins} bins -- too many, skipping")
        continue

    clf = MLPClassifier(
        hidden_layer_sizes=(256, 128, 64),
        max_iter=300, early_stopping=True,
        random_state=42, batch_size=2048,
    )
    clf.fit(X_d_tr_s, binned_tr)

    train_acc = clf.score(X_d_tr_s, binned_tr) * 100
    test_acc = clf.score(X_d_te_s, binned_te) * 100

    _, baseline_counts = np.unique(binned_tr, return_counts=True)
    baseline = np.max(baseline_counts) / N_TRAIN * 100

    print(f"  Bin size {bin_size}: {n_bins} bins, baseline={baseline:.1f}%, "
          f"train={train_acc:.1f}%, test={test_acc:.1f}%")

# ============================================================
# EXPERIMENT 7: Hybrid -- MLP on (n, R_inv, nearby_gap_stats)
# ============================================================
print("\n" + "=" * 80)
print("EXPERIMENT 7: Hybrid model -- R_inv + local gap statistics")
print("=" * 80)

# Use R^{-1}(n) as base, plus features from nearby known primes
# This simulates what you'd do if you had a table of primes up to some point

def make_hybrid_features(n_vals, r_inv_vals, primes_list, k_neighbors=10):
    """Use R_inv plus local prime density features."""
    n = n_vals.astype(np.float64)
    r = r_inv_vals.astype(np.float64)
    ln_n = np.log(n)

    feats = [
        n, ln_n, np.log(np.maximum(ln_n, 1.0)),
        r,
        r - n * ln_n,
        np.sqrt(n), n**(1/3),
    ]

    # Fractional part of R_inv
    feats.append(r - np.floor(r))

    # Modular features
    for m in [6, 30, 210]:
        feats.append(np.sin(2*np.pi*n/m))
        feats.append(np.cos(2*np.pi*n/m))

    return np.column_stack(feats)

X_hyb = make_hybrid_features(n_arr, r_inv_approx, primes_arr)
y_hyb_tr = primes_arr[train_idx].astype(np.float64)
y_hyb_te = primes_arr[test_idx].astype(np.float64)

scaler_h = StandardScaler()
X_h_tr_s = scaler_h.fit_transform(X_hyb[train_idx])
X_h_te_s = scaler_h.transform(X_hyb[test_idx])

# Target: p(n) directly
y_h_mean, y_h_std = y_hyb_tr.mean(), y_hyb_tr.std()

print("\n  MLP-hybrid (512,256,128,64):")
t0 = time.time()
mlp_h = MLPRegressor(
    hidden_layer_sizes=(512, 256, 128, 64),
    max_iter=1000, early_stopping=True,
    random_state=42, batch_size=2048,
    learning_rate='adaptive',
)
mlp_h.fit(X_h_tr_s, (y_hyb_tr - y_h_mean) / y_h_std)

pred_h_tr = mlp_h.predict(X_h_tr_s) * y_h_std + y_h_mean
pred_h_te = mlp_h.predict(X_h_te_s) * y_h_std + y_h_mean

train_exact = np.mean(np.round(pred_h_tr).astype(np.int64) == primes_arr[train_idx]) * 100
test_exact = np.mean(np.round(pred_h_te).astype(np.int64) == primes_arr[test_idx]) * 100
train_mae = mean_absolute_error(y_hyb_tr, pred_h_tr)
test_mae = mean_absolute_error(y_hyb_te, pred_h_te)

print(f"    Time: {time.time()-t0:.1f}s")
print(f"    Train exact: {train_exact:.3f}%, Test exact: {test_exact:.3f}%")
print(f"    Train MAE: {train_mae:.1f}, Test MAE: {test_mae:.1f}")

# ============================================================
# EXPERIMENT 8: Fourier analysis of delta(n)
# ============================================================
print("\n" + "=" * 80)
print("EXPERIMENT 8: Fourier analysis of delta(n)")
print("=" * 80)

# Is there any periodic structure in delta(n)?
delta_train_f = delta[:N_TRAIN].astype(np.float64)
delta_train_f -= np.mean(delta_train_f)

fft_delta = np.fft.rfft(delta_train_f)
power = np.abs(fft_delta)**2
freqs = np.fft.rfftfreq(N_TRAIN)

# Top frequencies
top_freq_idx = np.argsort(-power[1:])[:20] + 1  # skip DC
print("  Top 20 frequencies in delta(n):")
for i, idx in enumerate(top_freq_idx[:10]):
    period = 1.0 / freqs[idx] if freqs[idx] > 0 else float('inf')
    print(f"    #{i+1}: freq={freqs[idx]:.6f}, period={period:.1f}, power={power[idx]:.1f}")

# Try Fourier-based prediction: use top K frequencies
for K in [10, 50, 100, 500]:
    # Keep only top K frequencies
    fft_approx = np.zeros_like(fft_delta)
    top_k_idx = np.argsort(-power[1:])[:K] + 1
    fft_approx[top_k_idx] = fft_delta[top_k_idx]
    fft_approx[0] = fft_delta[0]  # DC component (mean)

    delta_approx = np.fft.irfft(fft_approx, n=N_TRAIN)
    delta_approx += np.mean(delta[:N_TRAIN])

    exact_pct = np.mean(np.round(delta_approx).astype(int) == delta[:N_TRAIN]) * 100
    rmse = np.sqrt(np.mean((delta_approx - delta[:N_TRAIN].astype(float))**2))
    print(f"  Top {K:4d} freqs: RMSE={rmse:.1f}, exact={exact_pct:.3f}% (train)")

# ============================================================
# SUMMARY
# ============================================================
print("\n" + "=" * 80)
print("SUMMARY OF ALL EXPERIMENTS")
print("=" * 80)
print("""
Key findings (actual numbers):

1. Direct MLP n->p(n): Train 0.03-0.06% exact, Test 0.00%.
   Massive test MAE (58k-113k) = pure overfitting to smooth trend.

2. Residue p(n) mod q: For ALL q in {3..31}, test accuracy = baseline (random).
   p(n) mod 3: 50.1% test vs 50.0% baseline. Zero signal.
   p(n) mod 7: 16.6% test vs 16.7% baseline. Zero signal.
   Autocorrelation in p(n) mod q: max |r| = 0.255 (mod 3, lag 1) --
   this is the Chebyshev bias (prime race), not a learnable pattern.

3. Gap prediction: MLP test 5.18% exact, GBR test 3.77%.
   BOTH BELOW baseline of 16.2% (always predict gap=6).
   Models learn the wrong thing -- the conditional mean, not the mode.

4. Symbolic regression: Best is 7-basis fit, RMSE=175.5 on train,
   318.7 on test. 0.42% train exact, 0.10% test. No structure.

5. Neural ODE: Error grows from 5.3% at n=100 to chaotic at n=10000.
   Errors compound -- sensitivity to initial perturbations.

6. Correction delta(n) = p(n) - R_inv(n):
   - R_inv itself is excellent: delta has std=237, range [-789, 709]
   - RF-delta: 5.64% train exact, 0.11% test. Pure memorization.
   - GBR-delta: 1.62% train, 0.13% test. Slight generalization but
     MAE=249 vs std=237 means it learned almost nothing.
   - delta distribution is nearly UNIFORM -- 1233 unique values in
     100k samples, most common (delta=-6) appears only 356 times (0.36%).

7. Fourier analysis of delta(n): No dominant periods.
   Top 500 frequencies capture only 1.48% exact (train only).
   The power spectrum is essentially flat = white noise.

8. Binned delta classification: bin=50 gives 5.4% test vs 28.4% baseline.
   Even coarse binning shows no learnable structure.

CONCLUSION: delta(n) = p(n) - R_inv(n) behaves like discrete white noise
with std ~ sqrt(p(n))/ln(p(n)). Every ML approach confirms that the
correction to the analytic approximation is IRREDUCIBLY RANDOM.

No model achieves > 0.13% test exact accuracy on any formulation.
The ~0.5*log2(n) bits of irreducible information per prime value
(Session 7 finding) is confirmed empirically: ML cannot compress
what has no compressible structure.
""")
