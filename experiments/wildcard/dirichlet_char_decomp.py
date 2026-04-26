"""
Dirichlet-character decomposition of the prime indicator chi_P modulo q.

Idea: any function f: (Z/qZ)^* -> C decomposes as
    f(n) = (1/phi(q)) * sum_{chi mod q} f_hat(chi) * chi(n),
where f_hat(chi) = sum_{a in (Z/qZ)^*} f(a) chi^{-1}(a).

If chi_P (restricted to (Z/qZ)^*) is "sparse" in the character basis
(say, only O(polylog q) nonzero f_hat), then evaluating chi_P at any n
mod q would take only O(polylog q) character evaluations. Combined with
a finer-grained recovery (e.g., over multiple moduli or via CRT), this
could yield a fast approximation/exact recovery of pi(x).

The well-known fact: chi_P restricted to (Z/qZ)^* is "almost balanced"
across residue classes (Dirichlet's theorem) but the small-scale
fluctuations are NOT character-sparse - they are random-like in
character coefficient. This experiment quantifies the sparsity.

We compute, for q in a range of moduli:
- The full Dirichlet-character spectrum |f_hat(chi)|^2 of chi_P on a
  window [N, N + q*M] mod q.
- The number of "significant" coefficients (above noise threshold).
- The L1/L2 ratio (small => sparse, large => dense).

Verdict: if number of significant coeffs grows like q (dense), no
shortcut. If it grows polylog, this is a path forward.
"""

import numpy as np
from math import gcd
from sympy import isprime, totient


def dirichlet_chars(q):
    """Build all Dirichlet characters mod q as a phi(q) x phi(q) matrix.
    Rows are characters, columns are residues (only the units).
    Uses a CRT-based decomposition for cyclic q only; for the test we
    use q with Z/qZ^* cyclic (q = prime, or 2,4,p^k,2p^k).
    """
    # For simplicity use q prime, so (Z/qZ)^* is cyclic of order q-1.
    if not isprime(q):
        raise ValueError("Use prime q only for this test.")
    # Find a primitive root g
    from sympy.ntheory.residue_ntheory import primitive_root
    g = primitive_root(q)
    # Build discrete log table
    units = []
    dlog = {}
    cur = 1
    for k in range(q - 1):
        units.append(cur)
        dlog[cur] = k
        cur = (cur * g) % q
    n = q - 1
    # chi_k(g^j) = exp(2 pi i k j / n)
    K = np.arange(n).reshape(-1, 1)
    J = np.array([dlog[u] for u in units]).reshape(1, -1)
    M = np.exp(2j * np.pi * K * J / n)
    return units, M


def chi_P_mod_q(q, window_start=0, window_size=None):
    """Empirical density of primes in each residue class mod q over
    [window_start, window_start + window_size). Returns vector of length
    phi(q) (only units), values in [0, count]."""
    if window_size is None:
        window_size = q * 50  # average ~50 hits per class on average
    counts = np.zeros(q, dtype=np.int64)
    for n in range(window_start, window_start + window_size):
        if n < 2:
            continue
        if isprime(n):
            counts[n % q] += 1
    return counts


def main():
    print("=== Dirichlet-character spectrum of chi_P ===\n")
    print(f"  q   phi(q)  significant_coefs  L1/L2  max_|coef|/mean  shape")

    for q in [11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97]:
        units, M = dirichlet_chars(q)
        # Get prime counts per unit residue
        ws = q * 200  # window large enough for stable signal
        counts = chi_P_mod_q(q, 0, ws)
        f = np.array([counts[u] for u in units], dtype=np.float64)
        # Subtract mean (this removes the principal-character contribution
        # which is "the smooth part" we already understand).
        f_centered = f - f.mean()
        # Compute Dirichlet transform: f_hat(chi) = sum_a f(a) chi^{-1}(a)
        # M[k, j] = chi_k(units[j]); chi^{-1}(a) = chi(a^{-1}) = conj(chi(a))
        # so f_hat = M.conj() @ f_centered
        f_hat = M.conj() @ f_centered
        spec = np.abs(f_hat) ** 2
        # Skip principal char (k=0)
        nonprincipal = spec[1:]
        total = nonprincipal.sum()
        if total == 0:
            continue
        # "significant" = above noise floor (mean + 2*std of nonprincipal)
        thr = nonprincipal.mean() + 2.0 * nonprincipal.std()
        sig = int((nonprincipal > thr).sum())
        # L1/L2
        l1 = np.abs(f_hat[1:]).sum()
        l2 = np.sqrt(spec[1:].sum())
        ratio = l1 / max(l2, 1e-12)
        m_over_mean = nonprincipal.max() / max(nonprincipal.mean(), 1e-12)
        # Shape: how concentrated? "uniform" if ratio close to sqrt(phi(q)-1)
        ideal = np.sqrt(q - 2)  # uniform random would give ~sqrt(N)
        shape = "uniform" if ratio > 0.85 * ideal else "concentrated"
        print(f"  {q:3d}  {q-1:5d}     {sig:5d}        {ratio:6.2f}  {m_over_mean:6.2f}        {shape} (ideal={ideal:.2f})")

    print("\nVerdict: sparsity of chi_P in character basis directly tracks")
    print("how 'random' chi_P looks across residue classes. If 'significant_coefs'")
    print("grows linearly in phi(q), the spectrum is dense and there is no")
    print("polylog shortcut from this angle. If 'L1/L2 ratio' approaches sqrt(phi(q))")
    print("the characters are equally weighted i.e. white-noise like.")


if __name__ == "__main__":
    main()
