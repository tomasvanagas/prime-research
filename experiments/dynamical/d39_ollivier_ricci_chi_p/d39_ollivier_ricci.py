"""
D39 — Ollivier-Ricci curvature of prime-indexed Cayley graphs.

Computes κ_O(x, y) := 1 − W_1(m_x, m_y) for every edge (x, y) of
G_P = Cay(Z/NZ, ±primes < N^c) and a matched random Cayley control,
where m_x is the uniform measure on the neighbourhood N(x).

By Cayley translation invariance, κ depends only on the generator
g = y − x ∈ S, so we compute one W_1 per (g, −g) pair: |S|/2 values.

Cross-domain reference: Ollivier 2009, J. Funct. Anal. 256
"Ricci curvature of Markov chains on metric spaces."
W_1 LP solved with scipy.optimize.linprog (HiGHS backend).

Run:
  python d39_ollivier_ricci.py [--Ns 64,128,256] [--cs 0.5,0.7,1.0]
                               [--ncontrols 10] [--seed 0]
"""
import argparse
import json
import time
from collections import deque

import numpy as np
from scipy.optimize import linprog


def primes_lt(M):
    if M < 3:
        return []
    sieve = np.ones(M, dtype=bool)
    sieve[:2] = False
    for i in range(2, int(M**0.5) + 1):
        if sieve[i]:
            sieve[i * i :: i] = False
    return np.where(sieve)[0].tolist()


def cayley_S_primes(N, c):
    """Return sorted generator set S = ±{primes p < N^c} ⊂ Z/NZ \\ {0}."""
    M = max(int(N**c), 2)
    P = primes_lt(M + 1)
    S = set()
    for p in P:
        if p >= 2 and p < M + 1:
            S.add(p % N)
            S.add((-p) % N)
    S.discard(0)
    return sorted(S)


def matched_window_random_S(N, S_template, M_window, rng):
    """Random S_R ⊂ ±{1,...,M_window} ⊂ Z/NZ\\{0}, parity-matched to S_template.

    Each generator a ∈ [1, M_window] contributes the pair {a mod N, (-a) mod N}.
    This is the FAIR control: same window size as S_P = ±{primes < M_window},
    same parity profile, just randomly chosen within the window."""
    pos_reps_template = set()
    for s in S_template:
        a = s if s <= N - s else N - s
        pos_reps_template.add(a)
    n_pairs = len(pos_reps_template)
    even_pairs = sum(1 for a in pos_reps_template if a % 2 == 0)
    odd_pairs = n_pairs - even_pairs

    # Canonical reps: a ≤ M_window AND a < N - a (so each a generates a distinct
    # ± pair of size 2 in S_R)
    upper = min(M_window, (N - 1) // 2)
    even_pool = [a for a in range(2, upper + 1) if a % 2 == 0]
    odd_pool = [a for a in range(1, upper + 1) if a % 2 == 1]

    if even_pairs > len(even_pool) or odd_pairs > len(odd_pool):
        return None  # window too small
    chosen_even = rng.choice(len(even_pool), size=even_pairs, replace=False)
    chosen_odd = rng.choice(len(odd_pool), size=odd_pairs, replace=False)
    chosen = [even_pool[i] for i in chosen_even] + [odd_pool[i] for i in chosen_odd]
    S_R = set()
    for a in chosen:
        S_R.add(a % N)
        S_R.add((-a) % N)
    S_R.discard(0)
    return sorted(S_R)


def matched_random_S(N, S_template, rng):
    """Random S_R ⊂ Z/NZ \\ {0}, |S_R|=|S_template|, closed under negation,
    parity-matched to S_template."""
    even_count = sum(1 for s in S_template if s % 2 == 0)
    odd_count = len(S_template) - even_count

    # Build by sampling unordered ± pairs.
    # Even pairs: a even, -a mod N even (true if N even); pair = {a, N-a}.
    # Odd pairs: similar.
    needed_even_pairs = even_count // 2
    needed_odd_pairs = odd_count // 2

    # Available residues
    even_residues = [a for a in range(1, N) if a % 2 == 0 and a != N - a]
    odd_residues = [a for a in range(1, N) if a % 2 == 1 and a != N - a]
    # Pair canonical reps: for a < N-a, the pair is (a, N-a)
    even_reps = [a for a in even_residues if a < N - a]
    odd_reps = [a for a in odd_residues if a < N - a]

    def sample_pairs(reps, k):
        if k > len(reps):
            raise ValueError(f"asked {k} pairs from {len(reps)} reps")
        idx = rng.choice(len(reps), size=k, replace=False)
        return [reps[i] for i in idx]

    chosen_even = sample_pairs(even_reps, needed_even_pairs)
    chosen_odd = sample_pairs(odd_reps, needed_odd_pairs)

    S_R = set()
    for a in chosen_even + chosen_odd:
        S_R.add(a)
        S_R.add((N - a) % N)

    # Handle self-paired (a = N - a, i.e., a = N/2): if odd_count or even_count odd
    if even_count % 2 == 1 and N % 2 == 0:
        S_R.add(N // 2)
    return sorted(S_R)


def bfs_dist_from_zero(N, S):
    """Single-source BFS from 0 on Cay(Z/NZ, S). Returns dist[v]."""
    dist = np.full(N, -1, dtype=np.int32)
    dist[0] = 0
    q = deque([0])
    while q:
        u = q.popleft()
        du = dist[u]
        for s in S:
            v = (u + s) % N
            if dist[v] == -1:
                dist[v] = du + 1
                q.append(v)
    return dist


def w1_emd(p, q, C):
    """Earth-mover W_1 between discrete measures p, q with cost matrix C.

    Solves LP: min Σ C_ij π_ij  s.t.  Σ_j π_ij = p_i, Σ_i π_ij = q_j, π ≥ 0.
    """
    n, m = len(p), len(q)
    c_obj = C.flatten()
    A_eq = np.zeros((n + m, n * m))
    for i in range(n):
        A_eq[i, i * m : (i + 1) * m] = 1.0
    for j in range(m):
        A_eq[n + j, j::m] = 1.0
    b_eq = np.concatenate([p, q])
    res = linprog(c_obj, A_eq=A_eq, b_eq=b_eq, bounds=(0, None), method="highs")
    if not res.success:
        raise RuntimeError(f"LP failed: {res.message}")
    return float(res.fun)


def cayley_ricci_distribution(N, S):
    """Return list of (g_canonical, kappa) for each generator-pair {g, -g}.
    Also returns mean and std of κ across edges (each generator counted once)."""
    dist = bfs_dist_from_zero(N, S)
    if (dist == -1).any():
        # disconnected — should not happen for prime Cayley with 2 ∈ S
        raise RuntimeError("graph disconnected")

    S_arr = np.asarray(S, dtype=np.int64)
    p_unif = np.full(len(S), 1.0 / len(S))
    seen = set()
    kappa_list = []
    for g in S:
        if g in seen:
            continue
        seen.add(g)
        seen.add((-g) % N)
        # Cost: C[i, j] = dist[(s_j + g - s_i) mod N]
        diff = (S_arr[None, :] + g - S_arr[:, None]) % N
        C = dist[diff].astype(float)
        w = w1_emd(p_unif, p_unif, C)
        kappa = 1.0 - w  # d(0, g) = 1
        kappa_list.append((int(g), float(kappa)))
    kappas = np.array([k for _, k in kappa_list])
    return kappa_list, float(kappas.mean()), float(kappas.std()), float(kappas.min())


def run_pair(N, c, n_controls, seed):
    out = {"N": N, "c": c, "n_controls": n_controls, "seed": seed}
    S_P = cayley_S_primes(N, c)
    out["S_P_size"] = len(S_P)
    if len(S_P) == 0 or len(S_P) > 200:
        out["status"] = "skip (S size)"
        return out
    t0 = time.time()
    klist_P, mean_P, std_P, min_P = cayley_ricci_distribution(N, S_P)
    t_P = time.time() - t0
    out["prime"] = {
        "kappa_mean": mean_P,
        "kappa_std": std_P,
        "kappa_min": min_P,
        "kappa_per_g": klist_P,
        "elapsed_sec": t_P,
    }

    M_window = max(int(N**c), 2)
    rng = np.random.default_rng(seed)

    def collect(make_S, label):
        means, stds, mins = [], [], []
        for _ in range(n_controls):
            S_R = make_S()
            if S_R is None or len(S_R) != len(S_P):
                continue
            try:
                _, mean_R, std_R, min_R = cayley_ricci_distribution(N, S_R)
                means.append(mean_R)
                stds.append(std_R)
                mins.append(min_R)
            except RuntimeError:
                continue
        if not means:
            return {"label": label, "n_realisations": 0}
        block = {
            "label": label,
            "n_realisations": len(means),
            "kappa_mean_mean": float(np.mean(means)),
            "kappa_mean_std": float(np.std(means)),
            "kappa_min_mean": float(np.mean(mins)),
        }
        block["z_score"] = float(
            (mean_P - np.mean(means)) / (np.std(means) + 1e-12)
        )
        return block

    out["random_global"] = collect(
        lambda: matched_random_S(N, S_P, rng), "random_in_ZNZ"
    )
    out["random_window"] = collect(
        lambda: matched_window_random_S(N, S_P, M_window, rng),
        f"random_in_window_{M_window}",
    )
    out["M_window"] = M_window
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--Ns", default="64,128,256")
    ap.add_argument("--cs", default="0.5,0.7,1.0")
    ap.add_argument("--ncontrols", type=int, default=10)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--out", default="d39_results.json")
    args = ap.parse_args()

    Ns = [int(x) for x in args.Ns.split(",")]
    cs = [float(x) for x in args.cs.split(",")]
    runs = []
    for N in Ns:
        for c in cs:
            print(f"=== N={N}, c={c} ===")
            res = run_pair(N, c, args.ncontrols, args.seed + N + int(c * 100))
            runs.append(res)
            if "prime" in res:
                p = res["prime"]
                rg = res.get("random_global", {})
                rw = res.get("random_window", {})
                print(
                    f"  |S|={res['S_P_size']}  M_window={res['M_window']}  "
                    f"prime κ̄={p['kappa_mean']:+.4f} "
                    f"σ={p['kappa_std']:.4f} κ_min={p['kappa_min']:+.4f}  "
                    f"({p['elapsed_sec']:.1f}s)"
                )
                if rg.get("kappa_mean_mean") is not None:
                    print(
                        f"  random_global  κ̄={rg['kappa_mean_mean']:+.4f} "
                        f"± {rg['kappa_mean_std']:.4f}  z={rg['z_score']:+.2f}"
                    )
                if rw.get("kappa_mean_mean") is not None:
                    print(
                        f"  random_window κ̄={rw['kappa_mean_mean']:+.4f} "
                        f"± {rw['kappa_mean_std']:.4f}  z={rw['z_score']:+.2f}"
                    )

    with open(args.out, "w") as f:
        json.dump({"runs": runs, "args": vars(args)}, f, indent=2)
    print(f"\nWrote {args.out}")


if __name__ == "__main__":
    main()
