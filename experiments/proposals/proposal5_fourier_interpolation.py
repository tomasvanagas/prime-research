"""
PROPOSAL 5: Fourier Interpolation of the Oscillatory Part

IDEA: The oscillatory correction to pi(x) is essentially a quasiperiodic 
function when viewed in log-space:
  S(x) = pi(x) - R(x) ~ sum_k c_k * x^{rho_k}

In variable t = log(x), this becomes:
  S(e^t) ~ e^{t/2} * sum_k c_k * e^{i*gamma_k*t}

This is a GENERALIZED trigonometric polynomial in t, with "frequencies" 
gamma_1, gamma_2, ... (the imaginary parts of zeta zeros).

If we can learn the coefficients c_k from samples of S at known points,
we can PREDICT S at any new point in O(K) time.

KEY: Sampling S requires knowing pi(x) exactly at sample points. 
This costs O(x^{2/3}) per sample via Meissel-Lehmer.
But once the K coefficients are learned, prediction is O(K).

COST ANALYSIS:
- Preprocessing: O(K * x_max^{2/3}) to sample at K points  
- Per query: O(K) to evaluate the trig polynomial
- If K = O(polylog(x)), we get O(polylog) per query!
- But K must be large enough to capture all significant zeros...
"""

import numpy as np
from sympy import prime, primepi, mobius
from mpmath import mp, mpf, log, li, im, zetazero
import time

mp.dps = 30

# Precompute zeta zeros
print("Computing zeta zeros...")
K = 30
gammas = [float(im(zetazero(k))) for k in range(1, K+1)]
print(f"First 5 zeros: {[f'{g:.3f}' for g in gammas[:5]]}")

def R_function(x, terms=15):
    """Riemann R(x) = sum mu(k)/k * li(x^{1/k})"""
    x = mpf(x)
    result = mpf(0)
    for k in range(1, terms+1):
        mu_k = int(mobius(k))
        if mu_k != 0:
            result += mpf(mu_k)/k * li(x ** (mpf(1)/k))
    return float(result)

def oscillatory_part(x):
    """S(x) = pi(x) - R(x)"""
    return primepi(int(x)) - R_function(x)

print("\n" + "=" * 70)
print("PROPOSAL 5: Fourier Interpolation of Oscillatory Part")
print("=" * 70)

# Step 1: Sample S(x) at training points in the range [100, 10000]
print("\n--- Step 1: Sampling oscillatory correction ---")
t0 = time.time()

# Use logarithmically spaced points
train_x = np.unique(np.round(np.exp(np.linspace(np.log(100), np.log(8000), 200))).astype(int))
train_S = np.array([oscillatory_part(x) for x in train_x])
t1 = time.time()
print(f"Sampled {len(train_x)} points in {t1-t0:.1f}s")

# Step 2: Fit Fourier model S(x) = sum c_k e^{i gamma_k log(x)} / sqrt(x)
# In real form: S(x) = (1/sqrt(x)) * [a_0 + sum_k (a_k cos(g_k t) + b_k sin(g_k t))]
# where t = log(x)

print("\n--- Step 2: Fourier decomposition with zeta zeros ---")

for n_zeros in [5, 10, 15, 20, 30]:
    t_vals = np.log(train_x.astype(float))
    sqrt_x = np.sqrt(train_x.astype(float))
    y = train_S * sqrt_x  # normalize out the x^{1/2} decay
    
    # Design matrix
    A = np.ones((len(train_x), 2*n_zeros + 1))
    for k in range(n_zeros):
        A[:, 2*k+1] = np.cos(gammas[k] * t_vals)
        A[:, 2*k+2] = np.sin(gammas[k] * t_vals)
    
    # Solve
    coeffs, _, _, _ = np.linalg.lstsq(A, y, rcond=None)
    
    # Train error
    y_pred = A @ coeffs
    train_rmse = np.sqrt(np.mean((y - y_pred)**2))
    S_pred_train = y_pred / sqrt_x
    S_rmse_train = np.sqrt(np.mean((train_S - S_pred_train)**2))
    
    # Test on [8500, 10000]
    test_x = np.arange(8500, 10001, 5)
    test_S = np.array([oscillatory_part(x) for x in test_x])
    
    t_test = np.log(test_x.astype(float))
    sqrt_test = np.sqrt(test_x.astype(float))
    A_test = np.ones((len(test_x), 2*n_zeros + 1))
    for k in range(n_zeros):
        A_test[:, 2*k+1] = np.cos(gammas[k] * t_test)
        A_test[:, 2*k+2] = np.sin(gammas[k] * t_test)
    
    S_pred_test = (A_test @ coeffs) / sqrt_test
    test_rmse = np.sqrt(np.mean((test_S - S_pred_test)**2))
    
    print(f"  {n_zeros:>2} zeros ({2*n_zeros+1} params): train_S_RMSE={S_rmse_train:.3f}, test_S_RMSE={test_rmse:.3f}")

# Step 3: Use best model to predict p(n)
print("\n--- Step 3: Predict p(n) via interpolated S ---")

# Use 30 zeros
n_z = 30
t_vals = np.log(train_x.astype(float))
sqrt_x = np.sqrt(train_x.astype(float))
y = train_S * sqrt_x
A = np.ones((len(train_x), 2*n_z + 1))
for k in range(n_z):
    A[:, 2*k+1] = np.cos(gammas[k] * t_vals)
    A[:, 2*k+2] = np.sin(gammas[k] * t_vals)
coeffs, _, _, _ = np.linalg.lstsq(A, y, rcond=None)

correct = 0
total = 0
abs_errors = []

for n in range(100, 1201):
    p_actual = prime(n)
    
    # Estimate x = R^{-1}(n)
    x_est_mp = mpf(n) * log(mpf(n))
    for _ in range(30):
        rx = li(x_est_mp)
        dx = 1/log(x_est_mp)
        x_est_mp += (mpf(n) - rx) / dx
    x_est = float(x_est_mp)
    
    # Predict S(x_est)
    t_e = np.log(x_est)
    sq_e = np.sqrt(x_est)
    a_vec = np.ones(2*n_z + 1)
    for k in range(n_z):
        a_vec[2*k+1] = np.cos(gammas[k] * t_e)
        a_vec[2*k+2] = np.sin(gammas[k] * t_e)
    S_pred = (a_vec @ coeffs) / sq_e
    
    # pi(x_est) ≈ R(x_est) + S_pred
    # We want x such that pi(x) = n
    # R(x_est) ≈ n (by construction), so pi(x_est) ≈ n + S_pred
    # Need to shift x by about -S_pred * log(x_est) to compensate
    x_corrected = x_est - S_pred * np.log(x_est)
    
    x_int = int(round(x_corrected))
    # Find nearest prime
    from sympy import nextprime, prevprime
    if x_int < 2: x_int = 2
    
    candidates = []
    try: candidates.append(prevprime(x_int))
    except: pass
    if isprime(x_int): candidates.append(x_int)
    candidates.append(nextprime(x_int))
    
    # Pick the one closest to x_corrected
    from sympy import isprime
    best = min(candidates, key=lambda c: abs(c - x_corrected))
    
    error = best - p_actual
    abs_errors.append(abs(error))
    if best == p_actual:
        correct += 1
    total += 1

abs_errors = np.array(abs_errors)
print(f"\nResults for n in [100, 1200]:")
print(f"Exact: {correct}/{total} = {correct/total*100:.1f}%")
print(f"Mean |error|: {abs_errors.mean():.1f}")
print(f"Median |error|: {np.median(abs_errors):.1f}")
for thr in [0, 2, 5, 10, 20, 50, 100]:
    within = (abs_errors <= thr).sum()
    print(f"  |error| <= {thr:>3}: {within}/{total} ({within/total*100:.1f}%)")

