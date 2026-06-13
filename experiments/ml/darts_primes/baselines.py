"""
Reference loss baselines for PRIMES at N=12.

Computes BCE loss for several elementary predictors so DARTS final loss
can be interpreted:

  - Constant predictor (entropy floor)
  - 1-bit oddness predictor
  - mod-6 residue predictor
  - mod-30 residue predictor (1, 7, 11, 13, 17, 19, 23, 29 mod 30 = candidate primes)
  - mod-210 residue predictor (excluding multiples of 2, 3, 5, 7)

The DARTS solution is non-trivially better than baseline_k iff its loss is
< baseline_k AND the gap is statistically significant across seeds.

Also reports, for each predictor, the equivalent "discrete accuracy" if
predictions are thresholded at 0.5.
"""
import numpy as np
from sympy import isprime


def bce(y_hat, y, eps=1e-12):
    y_hat = np.clip(y_hat, eps, 1 - eps)
    return float(np.mean(-(y * np.log(y_hat) + (1 - y) * np.log(1 - y_hat))))


def predictor_from_residue_table(table_size, residue_set, density_global):
    """y_hat = (n mod table_size in residue_set ? p_residue : p_off)
    where p_residue = density of primes among integers with that residue,
    p_off = density of primes among other residues.
    """
    def fn(n_arr, y_arr):
        residues = n_arr % table_size
        in_set = np.isin(residues, list(residue_set))
        # Compute conditional densities.
        if in_set.sum() > 0:
            p_in = float(y_arr[in_set].mean())
        else:
            p_in = density_global
        if (~in_set).sum() > 0:
            p_off = float(y_arr[~in_set].mean())
        else:
            p_off = density_global
        y_hat = np.where(in_set, p_in, p_off)
        return y_hat, (p_in, p_off)
    return fn


def main(N=12):
    n_arr = np.arange(2 ** N, dtype=np.int64)
    y = np.array([1.0 if isprime(int(k)) else 0.0 for k in n_arr], dtype=np.float32)
    p = float(y.mean())
    print(f"N={N}, |[0, 2^N)| = {len(y)}, prime count = {int(y.sum())}, density = {p:.4f}")
    print()

    # 1. Constant.
    y_hat_const = np.full_like(y, p)
    print(f"Constant predictor (= density):   loss = {bce(y_hat_const, y):.4f}, acc = {1 - p:.4f}")

    # 2. Oddness.
    fn = predictor_from_residue_table(2, {1}, p)
    y_hat_odd, (p_in, p_off) = fn(n_arr, y)
    print(f"1-bit oddness:                    loss = {bce(y_hat_odd, y):.4f}, density(odd|prime)={p_in:.4f}")

    # 3. mod-6.
    fn = predictor_from_residue_table(6, {1, 5}, p)
    y_hat, (p_in, p_off) = fn(n_arr, y)
    print(f"mod-6 ∈ {{1, 5}}:                  loss = {bce(y_hat, y):.4f}, density(in|prime)={p_in:.4f}")

    # 4. mod-30.
    cand_30 = {1, 7, 11, 13, 17, 19, 23, 29}
    fn = predictor_from_residue_table(30, cand_30, p)
    y_hat, (p_in, p_off) = fn(n_arr, y)
    print(f"mod-30 wheel (8/30):              loss = {bce(y_hat, y):.4f}, density(in|prime)={p_in:.4f}")

    # 5. mod-210.
    primorials = [2, 3, 5, 7]
    coprime_to = lambda r: all(r % p != 0 for p in primorials)
    cand_210 = {r for r in range(210) if coprime_to(r) and r > 0}
    fn = predictor_from_residue_table(210, cand_210, p)
    y_hat, (p_in, p_off) = fn(n_arr, y)
    print(f"mod-210 wheel ({len(cand_210)}/210): loss = {bce(y_hat, y):.4f}, density(in|prime)={p_in:.4f}")

    # 6. Perfect (oracle) predictor on training set.
    print(f"Oracle (perfect on train):        loss ≈ 0.0000, acc = 1.0000")


if __name__ == "__main__":
    main()
