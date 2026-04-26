"""Szegedy quantum walks on prime/divisor/coprime graphs.

Cross-domain import: Szegedy (2004), "Quantum Speed-Up of Markov Chain Based
Algorithms" (arxiv quant-ph/0401053). Discrete-time quantum walk on a
reversible Markov chain P with transition matrix P; the unitary W(P) is
defined as a product of two reflections, with eigenvalues e^{+/- 2 i theta_k}
where cos(theta_k) are the eigenvalues of the discriminant matrix
D(x,y) = sqrt(P(x,y) P(y,x)). For reversible P, D has the same spectrum
as P (up to similarity), so the spectral gap of W(P) is sqrt(spectral
gap of P) -- the celebrated quadratic mixing speedup.

Application target: can a Szegedy walk on a number-theoretic graph
mix in polylog time, providing a polylog quantum algorithm for some
prime-related primitive?

We test three families:
  (A) Divisor graph D_x: vertices [1..x], edges (m,n) iff m|n and m!=n.
  (B) Coprime graph C_x: vertices [1..x], edges (m,n) iff gcd(m,n)=1
      and m!=n.
  (C) Multiplicative-group Cayley graph Cay((Z/N)*, S) for small N
      (compare against the S79 closure -- this checks the QUANTUM walk,
      not classical Laplacian).

For each graph and each x, we compute:
  1. Classical spectral gap delta(P) = 1 - lambda_2 of normalized walk.
  2. Predicted Szegedy gap: 2 * arcsin(sqrt(delta(P))) ~ 2 sqrt(delta).
  3. Empirical Szegedy mixing time: smallest t such that
        || (1/t) sum_{s=1..t} |W^s |psi_0> <psi_0| W^{-s}| - Pi || < 1/4
     where Pi is the stationary projector. (Cesaro-mean mixing.)
  4. Scaling of mixing time with x: fit to (log x)^a, x^b.
  5. Whether any eigenvector of D is supported predominantly on primes
     (= "primality witness" candidate).

Falsification: if any of D_x, C_x, or Cayley graph yields a Szegedy walk
with mixing time O(polylog(x)) AND a primality-correlated eigenvector,
we have a polylog quantum primality-extraction primitive.

Edges cited: E5.3 (TC^0 / circuit limits via AKS), E7.10 (modulus-twist
orthogonality), E7.12 (Cayley spectrum probes omega(n) not primality;
S79).

Negative-shape result if mixing time scales as poly(x) or the primality
correlation is at noise level: contributes E7.x line "Szegedy walks on
divisor / coprime / Cayley graphs do not provide polylog primality
extraction".
"""

from __future__ import annotations

import math
import time
import json
from dataclasses import dataclass, asdict
from pathlib import Path

import numpy as np
from scipy import sparse
from scipy.sparse.linalg import eigsh
from scipy.linalg import eigh
from sympy import isprime, primerange, totient, gcd


# ---------------------------------------------------------------------------
# Graph constructors
# ---------------------------------------------------------------------------

def divisor_graph(x: int) -> np.ndarray:
    """Adjacency matrix of D_x: vertices [1..x], edges (m,n) with m|n.

    The vertex 1 is connected to all -- a "hub" that dominates the
    spectrum. Returns dense float64 adjacency.
    """
    n = x
    A = np.zeros((n, n), dtype=np.float64)
    for m in range(1, n + 1):
        for k in range(2, n // m + 1):
            j = m * k
            A[m - 1, j - 1] = 1.0
            A[j - 1, m - 1] = 1.0
    return A


def coprime_graph(x: int) -> np.ndarray:
    """Adjacency of C_x: vertices [1..x], edges (m,n) iff gcd(m,n)=1, m!=n.

    Densest among the three families (~6/pi^2 fraction of edges).
    """
    n = x
    A = np.zeros((n, n), dtype=np.float64)
    for m in range(1, n + 1):
        for k in range(m + 1, n + 1):
            if math.gcd(m, k) == 1:
                A[m - 1, k - 1] = 1.0
                A[k - 1, m - 1] = 1.0
    return A


def cayley_unit_graph(N: int, generators: list[int]) -> tuple[np.ndarray, list[int]]:
    """Adjacency of Cay((Z/N)*, S) for S = generators (with inverses).

    Vertices: residues coprime to N (in (Z/N)^*).
    Returns (A, vertex_list) where vertex_list[i] is the residue.
    """
    units = [r for r in range(1, N) if math.gcd(r, N) == 1]
    idx = {r: i for i, r in enumerate(units)}
    n = len(units)
    full_gens = set()
    for g in generators:
        if math.gcd(g, N) != 1:
            continue
        full_gens.add(g % N)
        # multiplicative inverse
        for r in units:
            if (g * r) % N == 1:
                full_gens.add(r)
                break
    A = np.zeros((n, n), dtype=np.float64)
    for u in units:
        for g in full_gens:
            v = (u * g) % N
            if v in idx and v != u:
                A[idx[u], idx[v]] = 1.0
                A[idx[v], idx[u]] = 1.0
    return A, units


# ---------------------------------------------------------------------------
# Markov chain + spectral analysis
# ---------------------------------------------------------------------------

def lazy_walk_transition(A: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Lazy random walk: P = (I + D^{-1} A) / 2.

    Lazy walks avoid bipartite parity issues so the spectral gap matches
    mixing time directly. Returns (P, deg).
    """
    deg = A.sum(axis=1)
    n = A.shape[0]
    if (deg == 0).any():
        # add self-loops of weight 1 to isolated vertices
        for i in range(n):
            if deg[i] == 0:
                A = A.copy()
                A[i, i] = 1.0
                deg[i] = 1.0
    Dinv = 1.0 / deg
    P = 0.5 * np.eye(n) + 0.5 * (Dinv[:, None] * A)
    return P, deg


def discriminant_matrix(P: np.ndarray) -> np.ndarray:
    """D(x,y) = sqrt(P(x,y) * P(y,x)). Symmetric for reversible P."""
    return np.sqrt(P * P.T)


def spectral_gap_lazy(D: np.ndarray) -> tuple[float, np.ndarray, np.ndarray]:
    """Eigenvalue gap of D (= eigenvalues of lazy walk).

    For lazy reversible chains, eigenvalues of P (= eigenvalues of D
    for symmetric D) lie in [0,1] with 1 simple. Gap = 1 - lambda_2.

    Returns (gap, eigenvalues sorted descending, eigenvectors columns).
    """
    eigvals, eigvecs = eigh(D)
    order = np.argsort(eigvals)[::-1]
    eigvals = eigvals[order]
    eigvecs = eigvecs[:, order]
    gap = 1.0 - eigvals[1] if len(eigvals) > 1 else 1.0
    return gap, eigvals, eigvecs


def szegedy_phases(eigvals_D: np.ndarray) -> np.ndarray:
    """Szegedy walk eigenphases.

    For each eigenvalue lambda = cos(theta) of D (with |lambda|<=1),
    W(P) has eigenvalues e^{+/- 2 i theta}. The gap of W(P) is
    2 * arcsin(sqrt(1 - |lambda_2|)) ~ 2 sqrt(1 - lambda_2)
    in the small-gap regime -- the quadratic mixing speedup.
    """
    eigvals_clip = np.clip(eigvals_D, -1.0, 1.0)
    theta = np.arccos(eigvals_clip)
    return theta


# ---------------------------------------------------------------------------
# Primality correlation of eigenvectors
# ---------------------------------------------------------------------------

def primality_correlation(eigvecs: np.ndarray, vertex_list: list[int],
                          top_k: int = 30) -> dict:
    """For each of the top-k eigenvectors of D, compute the inner product
    of |v_i|^2 with the prime indicator chi_P (mass on prime vertices).

    Returns dict with arrays of correlations and a "max" score.
    """
    n = len(vertex_list)
    chi_P = np.array([1.0 if isprime(v) else 0.0 for v in vertex_list])
    p_count = chi_P.sum()
    expected = p_count / n  # expected mass for uniform vector

    k = min(top_k, eigvecs.shape[1])
    masses = np.zeros(k)
    for i in range(k):
        v = eigvecs[:, i]
        prob = (v * v)
        prob = prob / prob.sum()
        masses[i] = (prob * chi_P).sum()

    # z-score against expected mass under uniform
    # variance of mass for random unit vector is small; use ratio
    ratios = masses / max(expected, 1e-12)
    return {
        "expected_uniform_mass": float(expected),
        "eigenvector_masses": masses.tolist(),
        "ratios": ratios.tolist(),
        "max_ratio": float(np.max(np.abs(ratios - 1.0))),
        "argmax_eig": int(np.argmax(np.abs(ratios - 1.0))),
    }


# ---------------------------------------------------------------------------
# Mixing time (classical AND Szegedy proxy)
# ---------------------------------------------------------------------------

def classical_mixing_time(P: np.ndarray, eps: float = 0.25) -> int:
    """Smallest t such that ||P^t(x_0, *) - pi||_TV < eps for x_0 = vertex 0.

    Stationary pi is left eigenvector at eigenvalue 1.
    """
    n = P.shape[0]
    # stationary: solve pi P = pi
    eigvals, eigvecs = np.linalg.eig(P.T)
    idx = np.argmin(np.abs(eigvals - 1.0))
    pi = np.real(eigvecs[:, idx])
    pi = pi / pi.sum()

    p_t = np.zeros(n)
    p_t[0] = 1.0
    Pt = P.copy()
    p_t = p_t @ Pt  # one step
    t = 1
    max_t = 4 * n
    while t < max_t:
        tv = 0.5 * np.abs(p_t - pi).sum()
        if tv < eps:
            return t
        p_t = p_t @ P
        t += 1
    return t


def szegedy_mixing_time_estimate(gap_classical: float) -> float:
    """Theoretical Szegedy mixing time = O(1 / sqrt(gap)) (Szegedy 2004).

    Returns 1 / phase_gap where phase_gap = 2 arcsin(sqrt(gap)).
    """
    if gap_classical <= 0:
        return float("inf")
    g = min(gap_classical, 1.0)
    phase_gap = 2.0 * math.asin(math.sqrt(g))
    return 1.0 / phase_gap


# ---------------------------------------------------------------------------
# Single experiment runner
# ---------------------------------------------------------------------------

@dataclass
class Result:
    family: str
    x: int
    n_vertices: int
    n_edges: int
    spectral_gap: float
    classical_mixing_time: int
    szegedy_mixing_estimate: float
    quadratic_speedup_ratio: float
    primality_max_eigenvector_ratio: float
    primality_argmax_eig: int
    primality_count_in_graph: int
    runtime_s: float


def run_one(family: str, x: int, generators: list[int] | None = None) -> tuple[Result, dict]:
    t0 = time.time()
    if family == "divisor":
        A = divisor_graph(x)
        vertex_list = list(range(1, x + 1))
    elif family == "coprime":
        A = coprime_graph(x)
        vertex_list = list(range(1, x + 1))
    elif family == "cayley_unit":
        assert generators is not None
        A, vertex_list = cayley_unit_graph(x, generators)
    else:
        raise ValueError(family)

    P, deg = lazy_walk_transition(A)
    D = discriminant_matrix(P)
    gap, eigvals, eigvecs = spectral_gap_lazy(D)
    cmix = classical_mixing_time(P)
    smix = szegedy_mixing_time_estimate(gap)
    pcorr = primality_correlation(eigvecs, vertex_list, top_k=20)

    res = Result(
        family=family,
        x=x,
        n_vertices=A.shape[0],
        n_edges=int(A.sum() / 2),
        spectral_gap=float(gap),
        classical_mixing_time=int(cmix),
        szegedy_mixing_estimate=float(smix),
        quadratic_speedup_ratio=float(cmix / max(smix, 1e-9)),
        primality_max_eigenvector_ratio=float(pcorr["max_ratio"]),
        primality_argmax_eig=int(pcorr["argmax_eig"]),
        primality_count_in_graph=int(sum(1 for v in vertex_list if isprime(v))),
        runtime_s=time.time() - t0,
    )
    detail = {
        "result": asdict(res),
        "top_eigvals": eigvals[:20].tolist(),
        "primality_correlation": pcorr,
    }
    return res, detail


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

def main(out_dir: Path):
    out_dir.mkdir(parents=True, exist_ok=True)
    log_lines: list[str] = []

    def log(msg: str):
        print(msg, flush=True)
        log_lines.append(msg)

    all_results: list[Result] = []
    all_details: list[dict] = []

    log("=" * 78)
    log("Szegedy quantum walks on prime/divisor/coprime graphs")
    log("Cross-domain import: Szegedy (2004), discriminant matrix theorem.")
    log("=" * 78)

    # --- Family A: divisor graph D_x ---
    log("\n## Family A: Divisor graph D_x = ([1..x], m|n)")
    for x in [30, 50, 80, 120, 180, 250]:
        res, det = run_one("divisor", x)
        log(f"  x={x:4d}  V={res.n_vertices:4d}  E={res.n_edges:6d}  "
            f"gap={res.spectral_gap:.4f}  cmix={res.classical_mixing_time:4d}  "
            f"smix~{res.szegedy_mixing_estimate:6.2f}  "
            f"speedup={res.quadratic_speedup_ratio:5.2f}  "
            f"prime_max_ratio={res.primality_max_eigenvector_ratio:.3f}")
        all_results.append(res)
        all_details.append(det)

    # --- Family B: coprime graph C_x ---
    log("\n## Family B: Coprime graph C_x = ([1..x], gcd=1)")
    for x in [30, 50, 80, 120, 180, 250]:
        res, det = run_one("coprime", x)
        log(f"  x={x:4d}  V={res.n_vertices:4d}  E={res.n_edges:6d}  "
            f"gap={res.spectral_gap:.4f}  cmix={res.classical_mixing_time:4d}  "
            f"smix~{res.szegedy_mixing_estimate:6.2f}  "
            f"speedup={res.quadratic_speedup_ratio:5.2f}  "
            f"prime_max_ratio={res.primality_max_eigenvector_ratio:.3f}")
        all_results.append(res)
        all_details.append(det)

    # --- Family C: Cayley graph on (Z/NZ)^* ---
    log("\n## Family C: Cayley graph Cay((Z/NZ)*, {2,3,5,inv})")
    for N in [31, 61, 127, 211, 307, 499]:
        # Note: when N is prime, (Z/NZ)^* has N-1 elements (cyclic).
        # We use N-1 vertices.
        res, det = run_one("cayley_unit", N, generators=[2, 3, 5])
        log(f"  N={N:4d}  V={res.n_vertices:4d}  E={res.n_edges:6d}  "
            f"gap={res.spectral_gap:.4f}  cmix={res.classical_mixing_time:4d}  "
            f"smix~{res.szegedy_mixing_estimate:6.2f}  "
            f"speedup={res.quadratic_speedup_ratio:5.2f}  "
            f"prime_max_ratio={res.primality_max_eigenvector_ratio:.3f}")
        all_results.append(res)
        all_details.append(det)

    # --- Scaling fits ---
    log("\n" + "=" * 78)
    log("## Scaling laws: fit cmix and smix vs (log x)^a, x^b")
    log("=" * 78)

    for fam in ["divisor", "coprime", "cayley_unit"]:
        sub = [r for r in all_results if r.family == fam]
        if len(sub) < 3:
            continue
        xs = np.array([r.x for r in sub], dtype=float)
        cmix = np.array([r.classical_mixing_time for r in sub], dtype=float)
        smix = np.array([r.szegedy_mixing_estimate for r in sub], dtype=float)
        # log-log fit: cmix ~ x^b
        b_cmix = np.polyfit(np.log(xs), np.log(cmix + 1e-9), 1)[0]
        b_smix = np.polyfit(np.log(xs), np.log(smix + 1e-9), 1)[0]
        # log-loglog fit: cmix ~ (log x)^a
        a_cmix = np.polyfit(np.log(np.log(xs)), np.log(cmix + 1e-9), 1)[0]
        a_smix = np.polyfit(np.log(np.log(xs)), np.log(smix + 1e-9), 1)[0]
        log(f"  [{fam}] cmix ~ x^{b_cmix:.3f} or (log x)^{a_cmix:.3f}")
        log(f"  [{fam}] smix ~ x^{b_smix:.3f} or (log x)^{a_smix:.3f}")

    # --- Verdict on polylog mixing ---
    log("\n## Polylog test")
    log("If smix ~ (log x)^a for small a, Szegedy is polylog-mixing on this graph.")
    log("If smix ~ x^b for b > 0, Szegedy is poly(x)-mixing -- no polylog opening.")

    # --- Persist ---
    with open(out_dir / "results.json", "w") as f:
        json.dump({
            "results": [asdict(r) for r in all_results],
            "details": all_details,
        }, f, indent=2)
    with open(out_dir / "run.log", "w") as f:
        f.write("\n".join(log_lines))

    log(f"\nWrote results.json and run.log to {out_dir}")


if __name__ == "__main__":
    out = Path(__file__).parent
    main(out)
