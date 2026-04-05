"""
PROPOSAL 3: ML Surrogate Model for delta(n)

IDEA: delta(n) = p(n) - round(R^{-1}(n)) encodes the "random" part.
But what if delta(n) has learnable structure that's cheaper to compute
than the Riemann explicit formula?

We train a model on (n, features(n)) -> delta(n) and see if it generalizes.

FEATURES to try:
- n mod small primes (residue class features)
- Fractional parts of n * log(gamma_k) for first few zeta zeros gamma_k
- log(n), sqrt(n), n^{1/3} (scale features)
- Li^{-1}(n) fractional part
- Möbius-related features

KEY CONJECTURE: If delta(n) has a sparse representation in some feature
basis, a linear model (or shallow network) could learn it.

Even partial learning is useful: if ML predicts delta(n) to within ±C,
we reduce the sieve search to 2C candidates.
"""

import numpy as np
from sympy import prime, isprime
from mpmath import mp, mpf, log, li, im, zetazero
import time

mp.dps = 30

# Precompute zeta zeros
print("Computing zeta zeros...")
n_zeros = 20
gamma = [float(im(zetazero(k))) for k in range(1, n_zeros + 1)]
print(f"First 5 zeros: {gamma[:5]}")

def compute_R_inv(n):
    """Fast R^{-1}(n) approximation"""
    x = mpf(n) * log(mpf(n))
    for _ in range(30):
        rx = li(x)
        dx = 1 / log(x)
        x += (mpf(n) - rx) / dx
    return float(x)

def features(n, gamma_list):
    """Compute feature vector for integer n"""
    x = float(n)
    lx = np.log(x)
    feats = [
        lx,
        np.sqrt(x),
        x ** (1/3),
        1.0 / lx,
    ]
    # Residue class features
    for q in [2, 3, 5, 7, 11, 13, 17, 19, 23, 29]:
        feats.append(n % q)
    # Zeta zero features: fractional parts of x * gamma_k / (2*pi)
    for g in gamma_list:
        phase = x * g / (2 * np.pi)
        feats.append(np.sin(2 * np.pi * phase))
        feats.append(np.cos(2 * np.pi * phase))
    # Cross features
    for g in gamma_list[:5]:
        feats.append(np.sin(lx * g))
        feats.append(np.cos(lx * g))
    return np.array(feats)

# Generate training data
print("\nGenerating training data...")
t0 = time.time()

N_train = 3000
N_test = 1000

# Training: n from 10 to 3009
train_ns = list(range(10, 10 + N_train))
train_deltas = []
train_features = []

for n in train_ns:
    p = prime(n)
    r_inv = compute_R_inv(n)
    delta = p - round(r_inv)
    train_deltas.append(delta)
    train_features.append(features(n, gamma))

X_train = np.array(train_features)
y_train = np.array(train_deltas)

# Test: n from 3010 to 4009
test_ns = list(range(10 + N_train, 10 + N_train + N_test))
test_deltas = []
test_features = []

for n in test_ns:
    p = prime(n)
    r_inv = compute_R_inv(n)
    delta = p - round(r_inv)
    test_deltas.append(delta)
    test_features.append(features(n, gamma))

X_test = np.array(test_features)
y_test = np.array(test_deltas)

t1 = time.time()
print(f"Data generation: {t1-t0:.1f}s")

# Model 1: Ridge Regression
from numpy.linalg import lstsq, norm

# Normalize features
mu = X_train.mean(axis=0)
sigma = X_train.std(axis=0) + 1e-10
X_train_n = (X_train - mu) / sigma
X_test_n = (X_test - mu) / sigma

# Ridge regression
lam = 1.0
XtX = X_train_n.T @ X_train_n + lam * np.eye(X_train_n.shape[1])
Xty = X_train_n.T @ y_train
w = np.linalg.solve(XtX, Xty)

y_pred_train = X_train_n @ w
y_pred_test = X_test_n @ w

train_rmse = np.sqrt(np.mean((y_train - y_pred_train)**2))
test_rmse = np.sqrt(np.mean((y_test - y_pred_test)**2))
test_mae = np.mean(np.abs(y_test - y_pred_test))
test_max_err = np.max(np.abs(y_test - y_pred_test))

# Baseline: predict delta = 0
baseline_rmse = np.sqrt(np.mean(y_test**2))
baseline_mae = np.mean(np.abs(y_test))

print("\n" + "=" * 70)
print("PROPOSAL 3: ML Surrogate for delta(n)")
print("=" * 70)

print(f"\nFeature dim: {X_train.shape[1]}")
print(f"Train size: {N_train}, Test size: {N_test}")
print(f"\nBaseline (predict 0): RMSE={baseline_rmse:.2f}, MAE={baseline_mae:.2f}")
print(f"Ridge regression:     RMSE={test_rmse:.2f}, MAE={test_mae:.2f}, MaxErr={test_max_err:.2f}")
print(f"Improvement over baseline: {(1 - test_rmse/baseline_rmse)*100:.1f}%")

# How often would ML + rounding give exact answer?
exact_count = sum(1 for i in range(N_test) if round(y_pred_test[i]) == y_test[i])
print(f"\nExact match (round prediction): {exact_count}/{N_test} = {exact_count/N_test*100:.1f}%")

# How many candidates would we need to check?
for threshold in [1, 2, 5, 10, 20, 50]:
    within = sum(1 for i in range(N_test) if abs(y_pred_test[i] - y_test[i]) <= threshold)
    print(f"  Within ±{threshold:>2}: {within}/{N_test} = {within/N_test*100:.1f}%")

# Feature importance
print("\nTop 10 features by |weight|:")
feat_names = ['log(n)', 'sqrt(n)', 'n^(1/3)', '1/log(n)']
for q in [2, 3, 5, 7, 11, 13, 17, 19, 23, 29]:
    feat_names.append(f'n%{q}')
for k in range(n_zeros):
    feat_names.append(f'sin(n*g{k+1})')
    feat_names.append(f'cos(n*g{k+1})')
for k in range(5):
    feat_names.append(f'sin(ln*g{k+1})')
    feat_names.append(f'cos(ln*g{k+1})')

importance = np.abs(w)
top_idx = np.argsort(importance)[::-1][:10]
for i in top_idx:
    name = feat_names[i] if i < len(feat_names) else f'feat_{i}'
    print(f"  {name:>15}: weight={w[i]:>8.3f}, |w|={importance[i]:.3f}")

# Model 2: k-NN approach (local structure)
print("\n--- k-NN approach ---")
from scipy.spatial.distance import cdist

# Try different k values
for k in [1, 3, 5, 10, 20]:
    dists = cdist(X_test_n, X_train_n)
    knn_preds = []
    for i in range(N_test):
        nearest = np.argsort(dists[i])[:k]
        knn_preds.append(np.mean(y_train[nearest]))
    knn_preds = np.array(knn_preds)
    knn_rmse = np.sqrt(np.mean((y_test - knn_preds)**2))
    knn_exact = sum(1 for i in range(N_test) if round(knn_preds[i]) == y_test[i])
    print(f"  k={k:>2}: RMSE={knn_rmse:.2f}, exact={knn_exact}/{N_test}")

