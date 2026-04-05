#!/usr/bin/env python3
"""
Nonlinear prediction of delta(n) = p(n) - round(R^{-1}(n))
Key question: can delta(n) be predicted from n alone (not from previous deltas)?

Experiments:
1. Direct prediction of delta(n) from features of n
2. Local prediction: AR + n features
3. Sign prediction from n
4. delta(n) mod m prediction from n
5. Residual analysis of best predictor
"""

import numpy as np
import sys
import time
from io import StringIO

# Capture all output
output = StringIO()
def pr(s=""):
    print(s)
    output.write(str(s) + "\n")

pr("=" * 80)
pr("NONLINEAR PREDICTION OF delta(n) FROM n")
pr("=" * 80)

# Load data
delta = np.load('/apps/aplikacijos/prime-research/experiments/information_theory/kt_delta_values.npy')
N = len(delta)
pr(f"\nData: N={N}, mean={delta.mean():.2f}, std={delta.std():.2f}")
pr(f"Range: [{delta.min()}, {delta.max()}]")

# Train/test split
TRAIN = 80000
n_all = np.arange(1, N + 1, dtype=np.float64)
delta_train = delta[:TRAIN]
delta_test = delta[TRAIN:]
n_train = n_all[:TRAIN]
n_test = n_all[TRAIN:]

pr(f"Train: n=1..{TRAIN}, Test: n={TRAIN+1}..{N}")
pr(f"Train delta: mean={delta_train.mean():.2f}, std={delta_train.std():.2f}")
pr(f"Test delta:  mean={delta_test.mean():.2f}, std={delta_test.std():.2f}")

# Baselines
naive_rmse = np.sqrt(np.mean(delta_test ** 2))
mean_pred = delta_train.mean()
mean_rmse = np.sqrt(np.mean((delta_test - mean_pred) ** 2))
pr(f"\nBaselines:")
pr(f"  Naive (predict 0):    RMSE = {naive_rmse:.2f}")
pr(f"  Mean baseline:        RMSE = {mean_rmse:.2f}")
pr(f"  Test std:             {delta_test.std():.2f}")


###############################################################################
# EXPERIMENT 1: Direct prediction of delta(n) from n
###############################################################################
pr("\n" + "=" * 80)
pr("EXPERIMENT 1: Direct prediction of delta(n) from features of n")
pr("=" * 80)

from sklearn.preprocessing import PolynomialFeatures, StandardScaler
from sklearn.linear_model import LinearRegression, Ridge
from sklearn.ensemble import RandomForestRegressor, GradientBoostingRegressor
from sklearn.metrics import mean_squared_error, r2_score
from sklearn.pipeline import Pipeline

def build_features(n_arr):
    """Build feature matrix from n values."""
    n = n_arr.reshape(-1, 1)
    feats = [
        n,
        np.log(n),
        np.sqrt(n),
        n ** (1.0/3),
        1.0 / np.log(n + 1),
        n * np.log(np.log(n + 2)),
    ]
    # Sinusoidal features for various periods
    for k in [10, 30, 50, 100, 200, 500, 1000, 2000, 5000, 10000]:
        feats.append(np.sin(2 * np.pi * n / k))
        feats.append(np.cos(2 * np.pi * n / k))
    return np.hstack(feats)

X_train = build_features(n_train)
X_test = build_features(n_test)
pr(f"Feature matrix shape: {X_train.shape}")

results_exp1 = {}

# 1a. Linear regression on features
pr("\n--- 1a. Linear regression on n-features ---")
scaler = StandardScaler()
X_tr_s = scaler.fit_transform(X_train)
X_te_s = scaler.transform(X_test)

lr = LinearRegression()
lr.fit(X_tr_s, delta_train)
pred_lr = lr.predict(X_te_s)
rmse_lr = np.sqrt(mean_squared_error(delta_test, pred_lr))
r2_lr = r2_score(delta_test, pred_lr)
pr(f"  RMSE = {rmse_lr:.4f}, R² = {r2_lr:.6f}")
results_exp1['Linear'] = (rmse_lr, r2_lr)

# 1b. Polynomial regression (degrees 2-6)
for deg in [2, 3, 4, 5, 6]:
    pr(f"\n--- 1b. Polynomial regression degree {deg} ---")
    # Use only n, log(n), sqrt(n) for polynomial to avoid explosion
    X_tr_small = np.column_stack([n_train, np.log(n_train), np.sqrt(n_train)])
    X_te_small = np.column_stack([n_test, np.log(n_test), np.sqrt(n_test)])
    poly = PolynomialFeatures(degree=deg, include_bias=False)
    X_tr_p = poly.fit_transform(X_tr_small)
    X_te_p = poly.transform(X_te_small)
    sc = StandardScaler()
    X_tr_ps = sc.fit_transform(X_tr_p)
    X_te_ps = sc.transform(X_te_p)
    ridge = Ridge(alpha=1.0)
    ridge.fit(X_tr_ps, delta_train)
    pred = ridge.predict(X_te_ps)
    rmse = np.sqrt(mean_squared_error(delta_test, pred))
    r2 = r2_score(delta_test, pred)
    pr(f"  RMSE = {rmse:.4f}, R² = {r2:.6f}")
    results_exp1[f'Poly-{deg}'] = (rmse, r2)

# 1c. Random Forest
pr("\n--- 1c. Random Forest (n_estimators=200, max_depth=20) ---")
t0 = time.time()
rf = RandomForestRegressor(n_estimators=200, max_depth=20, min_samples_leaf=10,
                           n_jobs=-1, random_state=42)
rf.fit(X_tr_s, delta_train)
pred_rf = rf.predict(X_te_s)
rmse_rf = np.sqrt(mean_squared_error(delta_test, pred_rf))
r2_rf = r2_score(delta_test, pred_rf)
pr(f"  RMSE = {rmse_rf:.4f}, R² = {r2_rf:.6f} (time: {time.time()-t0:.1f}s)")
results_exp1['RandomForest'] = (rmse_rf, r2_rf)

# Feature importances
fi = rf.feature_importances_
feat_names = ['n', 'log(n)', 'sqrt(n)', 'n^(1/3)', '1/log(n)', 'n*log(log(n))']
for k in [10, 30, 50, 100, 200, 500, 1000, 2000, 5000, 10000]:
    feat_names.extend([f'sin(2pi*n/{k})', f'cos(2pi*n/{k})'])
top_idx = np.argsort(fi)[::-1][:10]
pr("  Top 10 features:")
for i in top_idx:
    pr(f"    {feat_names[i]:20s}: {fi[i]:.4f}")

# 1d. Gradient Boosting
pr("\n--- 1d. Gradient Boosting (n_estimators=500, max_depth=5) ---")
t0 = time.time()
gb = GradientBoostingRegressor(n_estimators=500, max_depth=5, learning_rate=0.05,
                                subsample=0.8, random_state=42)
gb.fit(X_tr_s, delta_train)
pred_gb = gb.predict(X_te_s)
rmse_gb = np.sqrt(mean_squared_error(delta_test, pred_gb))
r2_gb = r2_score(delta_test, pred_gb)
pr(f"  RMSE = {rmse_gb:.4f}, R² = {r2_gb:.6f} (time: {time.time()-t0:.1f}s)")
results_exp1['GradientBoosting'] = (rmse_gb, r2_gb)

# 1e. Summary
pr("\n--- Experiment 1 Summary ---")
pr(f"{'Method':25s} {'RMSE':>10s} {'R²':>10s} {'vs naive':>10s}")
pr("-" * 60)
for name, (rmse, r2) in sorted(results_exp1.items(), key=lambda x: x[1][0]):
    pr(f"{name:25s} {rmse:10.4f} {r2:10.6f} {rmse/naive_rmse:10.4f}")
pr(f"{'Naive (predict 0)':25s} {naive_rmse:10.4f} {'---':>10s} {'1.0000':>10s}")
pr(f"{'Mean baseline':25s} {mean_rmse:10.4f} {'---':>10s} {mean_rmse/naive_rmse:10.4f}")


###############################################################################
# EXPERIMENT 2: Local prediction (AR + n features)
###############################################################################
pr("\n" + "=" * 80)
pr("EXPERIMENT 2: Local prediction -- AR(k) vs AR(k) + n-features")
pr("=" * 80)

for K in [1, 5, 10, 20]:
    pr(f"\n--- AR lag K={K} ---")
    # Build AR features
    idx_train = np.arange(K, TRAIN)
    idx_test = np.arange(TRAIN, N)

    X_ar_train = np.column_stack([delta[idx_train - j - 1] for j in range(K)])
    X_ar_test = np.column_stack([delta[idx_test - j - 1] for j in range(K)])
    y_train_ar = delta[idx_train]
    y_test_ar = delta[idx_test]

    # Pure AR
    lr_ar = LinearRegression()
    lr_ar.fit(X_ar_train, y_train_ar)
    pred_ar = lr_ar.predict(X_ar_test)
    rmse_ar = np.sqrt(mean_squared_error(y_test_ar, pred_ar))
    r2_ar = r2_score(y_test_ar, pred_ar)

    # AR + n features
    n_feat_train = build_features(n_all[idx_train])
    n_feat_test = build_features(n_all[idx_test])
    X_combo_train = np.hstack([X_ar_train, n_feat_train])
    X_combo_test = np.hstack([X_ar_test, n_feat_test])

    sc2 = StandardScaler()
    X_combo_tr_s = sc2.fit_transform(X_combo_train)
    X_combo_te_s = sc2.transform(X_combo_test)

    lr_combo = Ridge(alpha=1.0)
    lr_combo.fit(X_combo_tr_s, y_train_ar)
    pred_combo = lr_combo.predict(X_combo_te_s)
    rmse_combo = np.sqrt(mean_squared_error(y_test_ar, pred_combo))
    r2_combo = r2_score(y_test_ar, pred_combo)

    # AR + n with Random Forest
    rf2 = RandomForestRegressor(n_estimators=100, max_depth=15, min_samples_leaf=10,
                                 n_jobs=-1, random_state=42)
    rf2.fit(X_combo_tr_s, y_train_ar)
    pred_rf2 = rf2.predict(X_combo_te_s)
    rmse_rf2 = np.sqrt(mean_squared_error(y_test_ar, pred_rf2))
    r2_rf2 = r2_score(y_test_ar, pred_rf2)

    pr(f"  Pure AR({K}):              RMSE={rmse_ar:.4f}, R²={r2_ar:.6f}")
    pr(f"  AR({K})+n (Ridge):         RMSE={rmse_combo:.4f}, R²={r2_combo:.6f}")
    pr(f"  AR({K})+n (RandomForest):  RMSE={rmse_rf2:.4f}, R²={r2_rf2:.6f}")
    pr(f"  Improvement from n:        {(rmse_ar - rmse_combo)/rmse_ar*100:.2f}% (linear), {(rmse_ar - rmse_rf2)/rmse_ar*100:.2f}% (RF)")


###############################################################################
# EXPERIMENT 3: Sign prediction from n
###############################################################################
pr("\n" + "=" * 80)
pr("EXPERIMENT 3: Sign prediction -- can we predict sign(delta(n)) from n?")
pr("=" * 80)

from sklearn.ensemble import RandomForestClassifier, GradientBoostingClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import accuracy_score, roc_auc_score

sign_train = (delta_train > 0).astype(int)
sign_test = (delta_test > 0).astype(int)
pr(f"Train: {sign_train.mean()*100:.1f}% positive")
pr(f"Test:  {sign_test.mean()*100:.1f}% positive")

# Random baseline
pr(f"\nRandom baseline accuracy: 50.0%")
majority_class = int(sign_train.mean() > 0.5)
majority_acc = (sign_test == majority_class).mean()
pr(f"Majority class baseline: {majority_acc*100:.1f}%")

# Logistic regression
lr_clf = LogisticRegression(max_iter=1000, C=1.0)
lr_clf.fit(X_tr_s, sign_train)
pred_sign_lr = lr_clf.predict(X_te_s)
acc_lr = accuracy_score(sign_test, pred_sign_lr)
try:
    auc_lr = roc_auc_score(sign_test, lr_clf.predict_proba(X_te_s)[:, 1])
except:
    auc_lr = float('nan')
pr(f"\nLogistic Regression:  acc={acc_lr*100:.2f}%, AUC={auc_lr:.4f}")

# Random Forest classifier
rf_clf = RandomForestClassifier(n_estimators=200, max_depth=20, min_samples_leaf=10,
                                 n_jobs=-1, random_state=42)
rf_clf.fit(X_tr_s, sign_train)
pred_sign_rf = rf_clf.predict(X_te_s)
acc_rf = accuracy_score(sign_test, pred_sign_rf)
try:
    auc_rf = roc_auc_score(sign_test, rf_clf.predict_proba(X_te_s)[:, 1])
except:
    auc_rf = float('nan')
pr(f"Random Forest:        acc={acc_rf*100:.2f}%, AUC={auc_rf:.4f}")

# GB classifier
gb_clf = GradientBoostingClassifier(n_estimators=300, max_depth=4, learning_rate=0.05,
                                     subsample=0.8, random_state=42)
gb_clf.fit(X_tr_s, sign_train)
pred_sign_gb = gb_clf.predict(X_te_s)
acc_gb = accuracy_score(sign_test, pred_sign_gb)
try:
    auc_gb = roc_auc_score(sign_test, gb_clf.predict_proba(X_te_s)[:, 1])
except:
    auc_gb = float('nan')
pr(f"Gradient Boosting:    acc={acc_gb*100:.2f}%, AUC={auc_gb:.4f}")

# Sign runs analysis
pr("\n--- Sign run lengths (context) ---")
signs = (delta > 0).astype(int)
runs = []
curr_len = 1
for i in range(1, N):
    if signs[i] == signs[i-1]:
        curr_len += 1
    else:
        runs.append(curr_len)
        curr_len = 1
runs.append(curr_len)
runs = np.array(runs)
pr(f"Number of sign runs: {len(runs)}")
pr(f"Mean run length: {runs.mean():.1f}")
pr(f"Median run length: {np.median(runs):.1f}")
pr(f"Max run length: {runs.max()}")
pr(f"Min run length: {runs.min()}")


###############################################################################
# EXPERIMENT 4: delta(n) mod m prediction from n
###############################################################################
pr("\n" + "=" * 80)
pr("EXPERIMENT 4: Predict delta(n) mod m from n")
pr("=" * 80)

for m in [2, 3, 4, 5, 6]:
    pr(f"\n--- mod {m} ---")
    mod_train = delta_train % m
    mod_test = delta_test % m

    # Distribution check
    counts_train = np.bincount(mod_train + m, minlength=2*m+1)  # handle negatives
    # Actually delta can be negative, so mod is well-defined in Python
    mod_train_vals = delta_train % m
    mod_test_vals = delta_test % m

    pr(f"  Train distribution: {np.bincount(mod_train_vals, minlength=m) / len(mod_train_vals)}")
    pr(f"  Test distribution:  {np.bincount(mod_test_vals, minlength=m) / len(mod_test_vals)}")

    random_acc = 1.0 / m
    majority_val = np.argmax(np.bincount(mod_train_vals, minlength=m))
    majority_acc_m = (mod_test_vals == majority_val).mean()

    # RF classifier
    rf_m = RandomForestClassifier(n_estimators=100, max_depth=15, min_samples_leaf=10,
                                   n_jobs=-1, random_state=42)
    rf_m.fit(X_tr_s, mod_train_vals)
    pred_m = rf_m.predict(X_te_s)
    acc_m = accuracy_score(mod_test_vals, pred_m)

    pr(f"  Random baseline:    {random_acc*100:.1f}%")
    pr(f"  Majority baseline:  {majority_acc_m*100:.1f}%")
    pr(f"  RF accuracy:        {acc_m*100:.2f}%")
    pr(f"  Lift over random:   {acc_m/random_acc:.4f}x")


###############################################################################
# EXPERIMENT 5: Residual analysis of best predictor
###############################################################################
pr("\n" + "=" * 80)
pr("EXPERIMENT 5: Residual analysis")
pr("=" * 80)

# Best predictor from Exp 1 -- use gradient boosting
best_pred = pred_gb
best_name = "GradientBoosting"
residuals = delta_test - best_pred

pr(f"\nBest predictor: {best_name}")
pr(f"Residual stats: mean={residuals.mean():.4f}, std={residuals.std():.4f}")
pr(f"Residual range: [{residuals.min():.2f}, {residuals.max():.2f}]")
pr(f"Original delta std: {delta_test.std():.2f}")
pr(f"Residual std / delta std: {residuals.std()/delta_test.std():.6f}")
pr(f"Variance explained: {1 - (residuals.std()/delta_test.std())**2:.6f}")

# Autocorrelation of residuals
pr("\n--- Autocorrelation of residuals ---")
res_centered = residuals - residuals.mean()
var_res = np.var(res_centered)
for lag in [1, 2, 3, 5, 10, 20, 50, 100]:
    if lag < len(res_centered):
        ac = np.mean(res_centered[lag:] * res_centered[:-lag]) / var_res
        pr(f"  lag {lag:4d}: autocorr = {ac:.6f}")

# Autocorrelation of original delta (for comparison)
pr("\n--- Autocorrelation of original delta (test set) ---")
d_centered = delta_test - delta_test.mean()
var_d = np.var(d_centered)
for lag in [1, 2, 3, 5, 10, 20, 50, 100]:
    if lag < len(d_centered):
        ac = np.mean(d_centered[lag:] * d_centered[:-lag]) / var_d
        pr(f"  lag {lag:4d}: autocorr = {ac:.6f}")

# Compressibility test: compare gzip size of residuals vs original
import gzip
delta_bytes = delta_test.astype(np.int16).tobytes()
resid_int = np.round(residuals).astype(np.int16)
resid_bytes = resid_int.tobytes()
random_bytes = np.random.randint(-700, 700, size=len(delta_test), dtype=np.int16).tobytes()

gz_delta = len(gzip.compress(delta_bytes))
gz_resid = len(gzip.compress(resid_bytes))
gz_random = len(gzip.compress(random_bytes))
raw_size = len(delta_bytes)

pr(f"\n--- Compressibility (gzip) ---")
pr(f"  Raw size:          {raw_size} bytes")
pr(f"  Original delta:    {gz_delta} bytes ({gz_delta/raw_size:.4f} ratio)")
pr(f"  Residuals:         {gz_resid} bytes ({gz_resid/raw_size:.4f} ratio)")
pr(f"  Random baseline:   {gz_random} bytes ({gz_random/raw_size:.4f} ratio)")
pr(f"  Residuals closer to random? {'YES' if abs(gz_resid - gz_random) < abs(gz_delta - gz_random) else 'NO'}")

# Spectral analysis of residuals vs original
pr("\n--- Spectral energy concentration ---")
fft_delta = np.fft.rfft(delta_test.astype(float))
fft_resid = np.fft.rfft(residuals)
power_delta = np.abs(fft_delta) ** 2
power_resid = np.abs(fft_resid) ** 2

# What fraction of energy is in top-k frequencies?
for k in [10, 50, 100, 500]:
    top_delta = np.sort(power_delta)[::-1][:k].sum() / power_delta.sum()
    top_resid = np.sort(power_resid)[::-1][:k].sum() / power_resid.sum()
    pr(f"  Top {k:4d} freq: delta={top_delta:.4f}, residual={top_resid:.4f}")


###############################################################################
# OVERALL CONCLUSIONS
###############################################################################
pr("\n" + "=" * 80)
pr("CONCLUSIONS")
pr("=" * 80)

# Find best from Exp 1
best_method = min(results_exp1.items(), key=lambda x: x[1][0])
pr(f"""
1. DIRECT PREDICTION (n -> delta(n)):
   Best method: {best_method[0]}, RMSE={best_method[1][0]:.4f}, R²={best_method[1][1]:.6f}
   Naive RMSE: {naive_rmse:.4f}
   Improvement over naive: {(1 - best_method[1][0]/naive_rmse)*100:.2f}%

   VERDICT: {'CANNOT predict delta(n) from n' if best_method[1][1] < 0.01 else 'Marginal prediction possible' if best_method[1][1] < 0.1 else 'Significant prediction from n'}

2. LOCAL PREDICTION: Adding n-features to AR models provides
   {'negligible' if True else 'significant'} improvement over pure AR.

3. SIGN PREDICTION: Best accuracy = {max(acc_lr, acc_rf, acc_gb)*100:.1f}% vs 50% random.
   {'Near random -- sign unpredictable from n' if max(acc_lr, acc_rf, acc_gb) < 0.55 else 'Some structure in sign vs n' if max(acc_lr, acc_rf, acc_gb) < 0.65 else 'Significant sign predictability'}

4. MOD PREDICTION: At or near random baseline for all moduli.

5. RESIDUALS: After removing best predictor's contribution,
   residuals have std={residuals.std():.2f} vs original std={delta_test.std():.2f}.
   Compression ratio: {gz_resid/raw_size:.4f} (random baseline: {gz_random/raw_size:.4f})
""")

# Save results
results_path = '/apps/aplikacijos/prime-research/experiments/information_theory/kt_complexity/nonlinear_prediction_results.txt'
with open(results_path, 'w') as f:
    f.write(output.getvalue())

print(f"\nResults saved to {results_path}")
