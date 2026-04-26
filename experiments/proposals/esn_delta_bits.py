"""
Proposal B2: Reservoir-computing (Echo State Network) prediction of
delta(n) bit residuals.

delta(n) = p(n) - round(li_inv(n))  has the chaotic ~sqrt(n) part that
no smooth approximant captures.  The conjecture: delta is the orbit of
a finite-dimensional dynamical system whose first coordinate equals
delta(n).  An ESN is the cheapest test: a fixed random recurrent network
with only the readout layer trained.  If an ESN beats the majority-class
baseline on bit-0 of delta(n), there *is* dynamical structure.

This script:
  1. Computes p(n) and li_inv(n) for n in [1, n_max] (n_max = 10_000 default).
  2. Forms delta(n), takes bit-0 (parity).
  3. Builds polylog feature vectors per n: (log n, R^{-1} fractional part,
     n mod p for small primes, prior-bit lookback).
  4. Trains an ESN with reservoir size {64, 128, 256} on the first 70%,
     evaluates on the last 30%.
  5. Compares against (a) majority-class baseline, (b) logistic regression
     on the same features, (c) random-readout baseline.

Win condition: ESN test accuracy > majority + 2sigma on the held-out tail
AND ESN > logistic on same features (to rule out "the features alone
explain the bit").

Usage:
    python3 esn_delta_bits.py
    python3 esn_delta_bits.py --n_max 5000 --reservoir 256 --seed 0
"""
from __future__ import annotations
import argparse
import json
import time
import numpy as np
import sympy
from mpmath import li, findroot, mp


def li_inv_int(n: int) -> float:
    mp.dps = 20
    if n < 2:
        return 2.0
    x0 = max(2.0, float(n) * np.log(max(n, 2)))
    try:
        return float(findroot(lambda y: li(y) - n, x0, tol=1e-4))
    except Exception:
        return x0


def primes_and_delta(n_max: int) -> tuple[np.ndarray, np.ndarray]:
    primes = np.zeros(n_max + 1, dtype=np.int64)
    p = 1
    for n in range(1, n_max + 1):
        p = sympy.nextprime(p)
        primes[n] = p
    print(f"[setup] computing li_inv on {n_max} points (slow)")
    t0 = time.time()
    li_inv_vec = np.zeros(n_max + 1)
    for n in range(2, n_max + 1):
        li_inv_vec[n] = li_inv_int(n)
        if n % 1000 == 0:
            print(f"  ...{n}/{n_max}  ({time.time()-t0:.1f}s)")
    delta = primes.astype(float) - li_inv_vec
    return primes, delta


def build_features(primes: np.ndarray, delta: np.ndarray,
                   small_primes: list[int],
                   lookback: int = 4) -> np.ndarray:
    n_max = len(primes) - 1
    rows = []
    for n in range(2 + lookback, n_max + 1):
        log_n = np.log(n)
        log_log_n = np.log(log_n)
        feats = [log_n, log_log_n,
                 (delta[n - 1] - delta[n - 2]) / max(np.sqrt(n), 1),
                 delta[n - 1] / max(np.sqrt(n), 1)]
        for q in small_primes:
            feats.append((n % q) / q)
        for k in range(1, lookback + 1):
            d_prev = delta[n - k]
            feats.append(int(int(round(d_prev)) & 1) - 0.5)
        rows.append(feats)
    return np.array(rows, dtype=np.float64)


def make_targets(delta: np.ndarray, lookback: int) -> np.ndarray:
    out = []
    for n in range(2 + lookback, len(delta)):
        bit = int(round(delta[n])) & 1
        out.append(bit)
    return np.array(out, dtype=np.int64)


class ESN:
    def __init__(self, in_dim: int, res_dim: int = 128,
                 spec_radius: float = 0.95, sparsity: float = 0.1,
                 leak: float = 0.3, seed: int = 0) -> None:
        rng = np.random.default_rng(seed)
        self.W_in = rng.standard_normal((res_dim, in_dim)) * 0.5
        W = rng.standard_normal((res_dim, res_dim))
        mask = rng.random(W.shape) < sparsity
        W = W * mask
        eigs = np.linalg.eigvals(W)
        max_abs = float(np.max(np.abs(eigs)))
        if max_abs > 0:
            W *= spec_radius / max_abs
        self.W = W
        self.res_dim = res_dim
        self.leak = leak
        self.h = np.zeros(res_dim)
        self.W_out: np.ndarray | None = None
        self.bias_out: float = 0.0

    def reset(self) -> None:
        self.h = np.zeros(self.res_dim)

    def step(self, u: np.ndarray) -> np.ndarray:
        self.h = (1 - self.leak) * self.h + self.leak * np.tanh(
            self.W @ self.h + self.W_in @ u)
        return self.h.copy()

    def collect(self, U: np.ndarray) -> np.ndarray:
        H = np.zeros((U.shape[0], self.res_dim))
        for t in range(U.shape[0]):
            H[t] = self.step(U[t])
        return H

    def fit(self, U: np.ndarray, y: np.ndarray, ridge: float = 1e-4) -> None:
        self.reset()
        H = self.collect(U)
        H_aug = np.column_stack([H, np.ones(len(H))])
        # ridge regression
        A = H_aug.T @ H_aug + ridge * np.eye(H_aug.shape[1])
        b = H_aug.T @ (y.astype(float) - 0.5)
        w = np.linalg.solve(A, b)
        self.W_out = w[:-1]
        self.bias_out = float(w[-1])

    def predict(self, U: np.ndarray) -> np.ndarray:
        # do NOT reset --- carry hidden state from training
        H = self.collect(U)
        if self.W_out is None:
            raise RuntimeError("not fit")
        scores = H @ self.W_out + self.bias_out
        return (scores > 0).astype(np.int64)


def logistic_fit(X: np.ndarray, y: np.ndarray, ridge: float = 1.0,
                 iters: int = 200, lr: float = 0.05) -> tuple[np.ndarray, float]:
    n, d = X.shape
    w = np.zeros(d)
    b = 0.0
    y2 = y.astype(float)
    for _ in range(iters):
        z = X @ w + b
        p = 1.0 / (1.0 + np.exp(-z))
        grad_w = X.T @ (p - y2) / n + ridge * w
        grad_b = float(np.mean(p - y2))
        w -= lr * grad_w
        b -= lr * grad_b
    return w, b


def logistic_predict(X: np.ndarray, w: np.ndarray, b: float) -> np.ndarray:
    z = X @ w + b
    return (z > 0).astype(np.int64)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--n_max", type=int, default=5000)
    parser.add_argument("--reservoir", type=int, default=128)
    parser.add_argument("--leak", type=float, default=0.3)
    parser.add_argument("--spec_radius", type=float, default=0.95)
    parser.add_argument("--lookback", type=int, default=4)
    parser.add_argument("--train_frac", type=float, default=0.7)
    parser.add_argument("--seed", type=int, default=0)
    args = parser.parse_args()

    print(f"[setup] n_max={args.n_max}, reservoir={args.reservoir}")
    primes, delta = primes_and_delta(args.n_max)

    small_primes = [3, 5, 7, 11, 13]
    X = build_features(primes, delta, small_primes, args.lookback)
    y = make_targets(delta, args.lookback)
    assert X.shape[0] == y.shape[0]

    n_train = int(args.train_frac * len(y))
    X_tr, y_tr = X[:n_train], y[:n_train]
    X_te, y_te = X[n_train:], y[n_train:]

    print(f"[setup] train={len(y_tr)} test={len(y_te)} feature_dim={X.shape[1]}")
    print(f"[base]  bit balance train: {y_tr.mean():.4f}, "
          f"test: {y_te.mean():.4f}")

    # majority class baseline
    majority = int(round(y_tr.mean()))
    base_acc = float(np.mean((np.full_like(y_te, majority) == y_te)))
    base_majority = max(y_te.mean(), 1 - y_te.mean())
    n_te = len(y_te)
    sigma_majority = float(np.sqrt(0.25 / n_te))
    print(f"[base]  majority-class test acc = {base_acc:.4f}  "
          f"(unconditional bound = {base_majority:.4f}, "
          f"+/- {2*sigma_majority:.4f} = 2sigma)")

    # logistic baseline on same features
    w_lr, b_lr = logistic_fit(X_tr, y_tr)
    lr_pred = logistic_predict(X_te, w_lr, b_lr)
    lr_acc = float(np.mean(lr_pred == y_te))
    print(f"[base]  logistic regression test acc = {lr_acc:.4f}")

    # ESN
    esn = ESN(in_dim=X.shape[1], res_dim=args.reservoir,
              spec_radius=args.spec_radius, leak=args.leak, seed=args.seed)
    print(f"[esn]   training ridge readout...")
    esn.fit(X_tr, y_tr, ridge=1e-4)
    esn_train_pred = esn.predict(X_tr)
    train_acc = float(np.mean(esn_train_pred == y_tr))
    # for test: continue hidden state from end of train
    esn_te_pred = esn.predict(X_te)
    esn_acc = float(np.mean(esn_te_pred == y_te))
    print(f"[esn]   train acc = {train_acc:.4f}, test acc = {esn_acc:.4f}")

    # random-readout baseline
    rng = np.random.default_rng(123)
    esn_rand = ESN(in_dim=X.shape[1], res_dim=args.reservoir,
                   spec_radius=args.spec_radius, leak=args.leak, seed=999)
    esn_rand.W_out = rng.standard_normal(args.reservoir) * 0.01
    esn_rand.bias_out = 0.0
    rand_pred = esn_rand.predict(X_te)
    rand_acc = float(np.mean(rand_pred == y_te))
    print(f"[base]  random-readout ESN test acc = {rand_acc:.4f}")

    # confidence assessment
    advantage = esn_acc - base_majority
    z_score = advantage / sigma_majority
    print(f"\n[verdict]")
    print(f"  ESN advantage over unconditional majority = {advantage:+.4f}")
    print(f"  z-score = {z_score:+.2f}")
    if z_score > 2 and esn_acc > lr_acc:
        verdict = "TRUE WIN: ESN beats both majority and logistic"
    elif z_score > 2:
        verdict = ("WEAK WIN: ESN beats majority but logistic equally good "
                   "-- features alone explain it (no recurrent dynamics needed)")
    else:
        verdict = "CLOSED: ESN no advantage over baseline"
    print(f"  {verdict}")

    out = {
        "args": vars(args),
        "n_train": int(n_train),
        "n_test": int(len(y_te)),
        "feature_dim": int(X.shape[1]),
        "bit_balance_train": float(y_tr.mean()),
        "bit_balance_test": float(y_te.mean()),
        "majority_acc": base_acc,
        "majority_bound": base_majority,
        "logistic_acc": lr_acc,
        "esn_train_acc": train_acc,
        "esn_test_acc": esn_acc,
        "esn_random_readout_acc": rand_acc,
        "z_score": z_score,
        "verdict": verdict,
    }
    out_path = (
        "/apps/aplikacijos/prime-research/experiments/proposals/"
        "esn_delta_bits_data.json"
    )
    with open(out_path, "w") as f:
        json.dump(out, f, indent=2, default=float)
    print(f"\n[save] wrote {out_path}")


if __name__ == "__main__":
    main()
