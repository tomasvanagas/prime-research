#!/usr/bin/env python3
"""
Minimum Boolean Circuit Size for pi(x) mod 2
=============================================

For N-bit inputs x in {0,...,2^N-1}, find the minimum number of gates
in a Boolean circuit computing f(x) = pi(x) mod 2.

Methods:
1. Exact BFS for N<=4: enumerate ALL reachable truth tables layer by layer
2. SAT-based upper bound verification (N<=7, confirm size works)
3. Multi-strategy synthesis for upper bounds (all N)
   - Shannon decomposition (top-down, recursive)
   - Beam-search BFS (bottom-up, Hamming-distance guided)
   - Best-of-many random synthesis attempts

Gate basis: {AND, OR, XOR} with free NOT on inputs.
"""

import sys
import time
import math
import random
from collections import defaultdict
import heapq

# ============================================================
# Truth table computation
# ============================================================

def sieve_primes(limit):
    if limit < 2:
        return []
    is_prime = bytearray(b'\x01') * (limit + 1)
    is_prime[0] = is_prime[1] = 0
    for i in range(2, int(limit**0.5) + 1):
        if is_prime[i]:
            is_prime[i*i::i] = bytearray(len(is_prime[i*i::i]))
    return [i for i in range(limit + 1) if is_prime[i]]

def compute_pi_mod2_table(N):
    limit = (1 << N) - 1
    primes = set(sieve_primes(limit))
    tt, count = [], 0
    for x in range(limit + 1):
        if x in primes: count += 1
        tt.append(count & 1)
    return tt

def compute_isprime_table(N):
    limit = (1 << N) - 1
    primes = set(sieve_primes(limit))
    return [1 if x in primes else 0 for x in range(limit + 1)]

def random_balanced_table(N, weight):
    size = 1 << N
    tt = [0] * size
    for i in random.sample(range(size), min(weight, size)):
        tt[i] = 1
    return tt

def tt_to_int(tt_list):
    val = 0
    for i, b in enumerate(tt_list):
        if b: val |= (1 << i)
    return val

def make_base_signals(N):
    size = 1 << N
    all_one = (1 << size) - 1
    base = {0: 0, all_one: 0}
    for i in range(N):
        tt = sum(1 << x for x in range(size) if (x >> i) & 1)
        base[tt] = 0
        base[all_one ^ tt] = 0
    return base

def hamming(a, b):
    return bin(a ^ b).count('1')


# ============================================================
# Method 1: Exact BFS for N <= 4
# ============================================================

def exact_min_circuit_bfs(target_int, N, max_gates=20, timeout=120):
    """Exact minimum by enumerating all reachable truth tables."""
    size = 1 << N
    all_one = (1 << size) - 1
    start = time.time()

    cost = make_base_signals(N)
    if target_int in cost:
        return 0

    for s in range(1, max_gates + 1):
        if time.time() - start > timeout:
            return None

        by_cost = defaultdict(list)
        for tt, c in cost.items():
            by_cost[c].append(tt)

        new_fns = {}
        for ca in range(s):
            cb = s - 1 - ca
            if cb < ca: break
            if ca not in by_cost or cb not in by_cost: continue

            la, lb = by_cost[ca], by_cost[cb]
            pairs = ([(la[i], la[j]) for i in range(len(la)) for j in range(i, len(la))]
                     if ca == cb else [(a, b) for a in la for b in lb])

            for a, b in pairs:
                for res in (a & b, a | b, a ^ b):
                    for val in (res, all_one ^ res):
                        if val == target_int:
                            return s
                        if val not in cost and val not in new_fns:
                            new_fns[val] = s
                if time.time() - start > timeout:
                    return None

            if len(new_fns) + len(cost) > 10_000_000:
                return None

        if not new_fns:
            return None
        cost.update(new_fns)
        if target_int in cost:
            return cost[target_int]

    return None


# ============================================================
# Method 2: SAT-based check (verify a given size works)
# ============================================================

def sat_verify_size(target_tt, N, s, timeout=60):
    """Verify that a circuit of s gates CAN compute target_tt (SAT = yes)."""
    try:
        from pysat.solvers import Glucose3
    except ImportError:
        return None

    size = 1 << N
    vc = [0]
    def nv():
        vc[0] += 1; return vc[0]

    t_and = [nv() for _ in range(s)]
    t_or  = [nv() for _ in range(s)]
    t_xor = [nv() for _ in range(s)]

    nb = lambda i: N + i
    sel_a, sel_b, neg_a, neg_b = [], [], [], []
    for i in range(s):
        n = nb(i)
        sel_a.append([nv() for _ in range(n)])
        sel_b.append([nv() for _ in range(n)])
        neg_a.append(nv()); neg_b.append(nv())

    gv = [[nv() for _ in range(size)] for _ in range(s)]
    va = [[nv() for _ in range(size)] for _ in range(s)]
    vb = [[nv() for _ in range(size)] for _ in range(s)]

    cls = []

    for i in range(s):
        cls.append([t_and[i], t_or[i], t_xor[i]])
        cls.append([-t_and[i], -t_or[i]])
        cls.append([-t_and[i], -t_xor[i]])
        cls.append([-t_or[i], -t_xor[i]])

    for i in range(s):
        n = nb(i)
        cls.append(sel_a[i][:]); cls.append(sel_b[i][:])
        for j1 in range(n):
            for j2 in range(j1+1, n):
                cls.append([-sel_a[i][j1], -sel_a[i][j2]])
                cls.append([-sel_b[i][j1], -sel_b[i][j2]])

    bv = [[(x >> j) & 1 for x in range(size)] for j in range(N)]

    for i in range(s):
        n = nb(i)
        for x in range(size):
            for j in range(n):
                if j < N:
                    v = bv[j][x]
                    if v:
                        cls.append([-sel_a[i][j], neg_a[i], va[i][x]])
                        cls.append([-sel_a[i][j], -neg_a[i], -va[i][x]])
                        cls.append([-sel_b[i][j], neg_b[i], vb[i][x]])
                        cls.append([-sel_b[i][j], -neg_b[i], -vb[i][x]])
                    else:
                        cls.append([-sel_a[i][j], neg_a[i], -va[i][x]])
                        cls.append([-sel_a[i][j], -neg_a[i], va[i][x]])
                        cls.append([-sel_b[i][j], neg_b[i], -vb[i][x]])
                        cls.append([-sel_b[i][j], -neg_b[i], vb[i][x]])
                else:
                    gi = j - N
                    gvar = gv[gi][x]
                    cls.append([-sel_a[i][j], neg_a[i], -va[i][x], gvar])
                    cls.append([-sel_a[i][j], neg_a[i], va[i][x], -gvar])
                    cls.append([-sel_a[i][j], -neg_a[i], -va[i][x], -gvar])
                    cls.append([-sel_a[i][j], -neg_a[i], va[i][x], gvar])
                    cls.append([-sel_b[i][j], neg_b[i], -vb[i][x], gvar])
                    cls.append([-sel_b[i][j], neg_b[i], vb[i][x], -gvar])
                    cls.append([-sel_b[i][j], -neg_b[i], -vb[i][x], -gvar])
                    cls.append([-sel_b[i][j], -neg_b[i], vb[i][x], gvar])

    for i in range(s):
        for x in range(size):
            a, b, g = va[i][x], vb[i][x], gv[i][x]
            cls.append([-t_and[i], -a, -b, g])
            cls.append([-t_and[i], a, -g])
            cls.append([-t_and[i], b, -g])
            cls.append([-t_or[i], a, b, -g])
            cls.append([-t_or[i], -a, g])
            cls.append([-t_or[i], -b, g])
            cls.append([-t_xor[i], a, b, -g])
            cls.append([-t_xor[i], -a, -b, -g])
            cls.append([-t_xor[i], -a, b, g])
            cls.append([-t_xor[i], a, -b, g])

    for x in range(size):
        if target_tt[x]: cls.append([gv[s-1][x]])
        else: cls.append([-gv[s-1][x]])

    solver = Glucose3()
    for c in cls: solver.add_clause(c)
    result = solver.solve()
    solver.delete()
    return result


def sat_verify_with_timeout(target_tt, N, s, timeout=15):
    """Run SAT verification in a subprocess with hard timeout."""
    import multiprocessing as mp

    def worker(tt, n, gates, result_queue):
        try:
            r = sat_verify_size(tt, n, gates, timeout=9999)
            result_queue.put(r)
        except Exception as e:
            result_queue.put(None)

    q = mp.Queue()
    p = mp.Process(target=worker, args=(target_tt, N, s, q))
    p.start()
    p.join(timeout=timeout)
    if p.is_alive():
        p.terminate()
        p.join(1)
        if p.is_alive():
            p.kill()
            p.join()
        return None
    try:
        return q.get_nowait()
    except:
        return None


# ============================================================
# Method 3a: Beam-search BFS synthesis
# ============================================================

def beam_search_synthesis(target_int, N, beam_width=2000, max_layers=30, timeout=60):
    """
    Bottom-up BFS with beam search. Keep only the beam_width most
    promising truth tables (by Hamming distance to target).
    """
    size = 1 << N
    all_one = (1 << size) - 1
    start = time.time()

    cost = make_base_signals(N)
    if target_int in cost:
        return 0

    neg_target = all_one ^ target_int

    for layer in range(1, max_layers + 1):
        if time.time() - start > timeout:
            break

        # Select beam: best functions by Hamming distance to target
        scored = sorted(cost.keys(),
                        key=lambda t: min(hamming(t, target_int), hamming(t, neg_target)))
        beam = scored[:beam_width]

        new_fns = {}
        for i in range(len(beam)):
            a = beam[i]
            ca = cost[a]
            for j in range(i, min(i + beam_width // 2, len(beam))):
                b = beam[j]
                cb = cost[b]
                new_cost = ca + cb + 1

                for res in (a & b, a | b, a ^ b):
                    for val in (res, all_one ^ res):
                        if val == target_int:
                            return new_cost
                        if val not in cost:
                            if val not in new_fns or new_fns[val] > new_cost:
                                new_fns[val] = new_cost

            if time.time() - start > timeout:
                break

        if not new_fns:
            break

        # Keep only best new functions
        if len(new_fns) > beam_width * 2:
            best_new = sorted(new_fns.items(),
                              key=lambda kv: (min(hamming(kv[0], target_int),
                                                   hamming(kv[0], neg_target)),
                                               kv[1]))
            new_fns = dict(best_new[:beam_width * 2])

        for tt, c in new_fns.items():
            if tt not in cost or cost[tt] > c:
                cost[tt] = c

        if target_int in cost:
            return cost[target_int]

    return cost.get(target_int, cost.get(neg_target))


# ============================================================
# Method 3b: Shannon decomposition synthesis
# ============================================================

def shannon_synthesis(target_int, N, timeout=60):
    """Top-down Shannon expansion synthesis."""
    size = 1 << N
    all_one = (1 << size) - 1
    start = time.time()

    base = make_base_signals(N)
    cache = dict(base)

    if target_int in cache:
        return 0

    def synth(tt, depth=0):
        if tt in cache: return cache[tt]
        neg = all_one ^ tt
        if neg in cache:
            cache[tt] = cache[neg]; return cache[tt]
        if time.time() - start > timeout or depth > N + 3:
            return None

        best = None
        for var in range(N):
            bit = 1 << var
            f0 = f1 = 0
            for x in range(size):
                if (tt >> x) & 1:
                    if x & bit:
                        f1 |= (1 << (x & ~bit)); f1 |= (1 << (x | bit))
                    else:
                        f0 |= (1 << (x & ~bit)); f0 |= (1 << (x | bit))

            c0 = synth(f0, depth + 1)
            c1 = synth(f1, depth + 1)
            if c0 is None or c1 is None: continue

            if f0 == f1: c = c0
            elif f0 == 0: c = c1 + 1
            elif f0 == all_one: c = c1 + 1
            elif f1 == 0: c = c0 + 1
            elif f1 == all_one: c = c0 + 1
            elif f0 == (all_one ^ f1): c = max(c0, c1) + 1
            else: c = c0 + c1 + 3

            if best is None or c < best: best = c

        if best is not None: cache[tt] = best
        return best

    return synth(target_int)


# ============================================================
# Method 3c: Direct DNF/CNF synthesis (guaranteed, always works)
# ============================================================

def dnf_synthesis(target_tt, N):
    """
    Build circuit from minimized DNF: for each 1 in truth table,
    create an AND of literals, then OR them all together.
    Uses tree reduction for AND/OR chains.
    Returns gate count.
    """
    size = 1 << N
    ones = [x for x in range(size) if target_tt[x]]
    zeros = [x for x in range(size) if not target_tt[x]]

    if not ones: return 0
    if not zeros: return 0
    if len(ones) == 1:
        # Single minterm: AND of N literals = N-1 AND gates
        return max(0, N - 1)
    if len(zeros) == 1:
        # Single maxterm: OR of N literals = N-1 OR gates
        return max(0, N - 1)

    # Use whichever is smaller (DNF or CNF)
    use_dnf = len(ones) <= len(zeros)
    terms = ones if use_dnf else zeros

    # For each term, compute the minterm/maxterm
    # Minterm: AND of all literals matching that assignment
    # Cost per term: N-1 AND gates (using tree)
    # Actually: we can share sub-expressions via tree structure

    # Simple upper bound: each term needs at most N-1 gates, then combine
    # k terms with tree of ORs needs k-1 gates
    # Total: k * (N-1) + (k-1) = k*N - 1

    # Better: use tree reduction
    # Each term is an AND of N literals -> N-1 gates with binary tree
    # But many terms share prefixes. Use a sum-of-products circuit.

    k = len(terms)

    # Shared-tree approach: for each variable position, many terms agree.
    # Simple approach: just compute the bound.
    gates_per_term = max(0, N - 1)  # AND tree
    combine_gates = max(0, k - 1)   # OR/AND tree to combine
    total = k * gates_per_term + combine_gates

    # If using CNF (negating), add 0 gates (output negation is free)
    return total


def lupanov_upper_bound(N):
    """
    Lupanov's theorem: any N-variable Boolean function has a circuit
    of size at most (1 + epsilon) * 2^N / N for large N.
    For small N, return a concrete upper bound.
    """
    return int(math.ceil(2**N / N * 1.5))  # conservative


# ============================================================
# Method 3d: Randomized greedy synthesis
# ============================================================

def random_greedy_synthesis(target_int, N, attempts=10, timeout=30):
    """
    Multiple random greedy attempts. Each attempt builds a circuit
    by repeatedly combining the two functions closest to target.
    """
    size = 1 << N
    all_one = (1 << size) - 1
    neg_target = all_one ^ target_int
    start = time.time()
    best = None

    for attempt in range(attempts):
        if time.time() - start > timeout: break

        lib = dict(make_base_signals(N))
        if target_int in lib: return 0

        for step in range(100):
            if time.time() - start > timeout: break

            # Pick functions to combine: those closest to target
            scored = sorted(lib.keys(),
                            key=lambda t: min(hamming(t, target_int),
                                               hamming(t, neg_target)))
            # Take top candidates with some randomness
            top_k = min(50, len(scored))
            candidates = scored[:top_k]

            # Add some random ones for diversity
            if len(scored) > top_k:
                extra = random.sample(scored[top_k:], min(10, len(scored) - top_k))
                candidates.extend(extra)

            improved = False
            random.shuffle(candidates)

            for i in range(min(len(candidates), 60)):
                for j in range(i, min(len(candidates), 60)):
                    a, b = candidates[i], candidates[j]
                    ca, cb = lib[a], lib[b]
                    nc = ca + cb + 1

                    for res in (a & b, a | b, a ^ b):
                        for val in (res, all_one ^ res):
                            if val == target_int:
                                if best is None or nc < best:
                                    best = nc
                                improved = True
                            elif val not in lib or lib[val] > nc:
                                lib[val] = nc
                                improved = True

            if target_int in lib:
                v = lib[target_int]
                if best is None or v < best:
                    best = v
                break

            if not improved:
                break

    return best


# ============================================================
# Combined: best of all methods
# ============================================================

def find_circuit_size(target_tt, N, timeout=300):
    """Returns (size, method, is_exact)."""
    size_2n = 1 << N
    all_one = (1 << size_2n) - 1
    target_int = tt_to_int(target_tt)
    start = time.time()

    base = make_base_signals(N)
    if target_int in base:
        return 0, "literal", True

    results = {}

    # Phase 1: Shannon decomposition (primary method, give it most time)
    t0 = time.time()
    shannon_timeout = min(timeout * 0.7, 180) if N >= 7 else min(30, timeout/5)
    ub_sh = shannon_synthesis(target_int, N, timeout=shannon_timeout)
    results['shannon'] = ub_sh
    print(f"    [shannon] {'UB='+str(ub_sh) if ub_sh else 'fail'} ({time.time()-t0:.1f}s)", end="", flush=True)

    # Phase 2: Beam search BFS (only useful for N <= 6)
    remaining = timeout - (time.time() - start)
    if N <= 6 and remaining > 5:
        t0 = time.time()
        ub_beam = beam_search_synthesis(target_int, N, beam_width=3000,
                                         timeout=min(remaining * 0.3, 60))
        results['beam'] = ub_beam
        print(f"  [beam] {'UB='+str(ub_beam) if ub_beam else 'fail'} ({time.time()-t0:.1f}s)", end="", flush=True)

    # Phase 3: Random greedy (only useful for N <= 5)
    remaining = timeout - (time.time() - start)
    if N <= 5 and remaining > 5:
        t0 = time.time()
        ub_rg = random_greedy_synthesis(target_int, N, attempts=15,
                                         timeout=min(remaining * 0.2, 30))
        results['rgreedy'] = ub_rg
        print(f"  [rgreedy] {'UB='+str(ub_rg) if ub_rg else 'fail'} ({time.time()-t0:.1f}s)", end="", flush=True)

    # Phase 3.5: DNF synthesis (always works, gives guaranteed bound)
    ub_dnf = dnf_synthesis(target_tt, N)
    results['dnf'] = ub_dnf
    print(f"  [dnf] UB={ub_dnf}", end="", flush=True)

    # Best upper bound
    best_ub = None
    best_method = None
    for name, val in results.items():
        if val is not None and (best_ub is None or val < best_ub):
            best_ub = val
            best_method = name

    # Phase 4: Exact BFS for small N
    exact = None
    remaining = timeout - (time.time() - start)
    if N <= 4 and remaining > 10:
        t0 = time.time()
        exact = exact_min_circuit_bfs(target_int, N,
                                       max_gates=best_ub if best_ub else 20,
                                       timeout=min(remaining * 0.5, 120))
        print(f"  [exact] {'EXACT='+str(exact) if exact is not None else 'fail'} ({time.time()-t0:.1f}s)", end="", flush=True)

    # Phase 5: SAT verification (only for N <= 5, solver blocks so use alarm)
    remaining = timeout - (time.time() - start)
    if exact is None and N <= 5 and best_ub is not None and remaining > 15:
        print(f"  [SAT]", end="", flush=True)
        # Try one step below best_ub
        r = sat_verify_with_timeout(target_tt, N, best_ub - 1, timeout=min(20, remaining * 0.3))
        if r is True:
            best_ub -= 1
            best_method = "SAT"
            print(f" improved={best_ub}", end="", flush=True)
            # Try again
            remaining = timeout - (time.time() - start)
            if remaining > 10:
                r2 = sat_verify_with_timeout(target_tt, N, best_ub - 1, timeout=min(15, remaining * 0.3))
                if r2 is True:
                    best_ub -= 1
                    print(f" improved={best_ub}", end="", flush=True)
        elif r is False:
            exact = best_ub  # confirmed optimal
            print(f" confirmed_exact={best_ub}", end="", flush=True)

    print()

    if exact is not None:
        return exact, "exact", True
    elif best_ub is not None:
        return best_ub, f"UB({best_method})", False
    else:
        return -1, "failed", False


# ============================================================
# Main experiment
# ============================================================

def run_experiment():
    print("=" * 70)
    print("Minimum Boolean Circuit Size for pi(x) mod 2")
    print("Gate basis: {AND, OR, XOR} with free input negation")
    print("=" * 70)

    all_results = {}
    TIMEOUT = 180  # per function

    for N in range(4, 11):
        print(f"\n{'='*60}")
        print(f"N = {N}  (domain 0..{(1<<N)-1}, |TT| = {1<<N})")
        print(f"{'='*60}")

        pi_tt = compute_pi_mod2_table(N)
        pr_tt = compute_isprime_table(N)
        random.seed(42 + N)
        rn_tt = random_balanced_table(N, sum(pi_tt))

        wpi, wpr, wrn = sum(pi_tt), sum(pr_tt), sum(rn_tt)
        print(f"  weights: pi_mod2={wpi}/{1<<N}, is_prime={wpr}/{1<<N}, random={wrn}/{1<<N}")

        for fn, tt in [("pi_mod2", pi_tt), ("is_prime", pr_tt), ("random", rn_tt)]:
            print(f"\n  --- {fn} (N={N}) ---")
            t0 = time.time()
            sz, method, is_exact = find_circuit_size(tt, N, timeout=TIMEOUT)
            elapsed = time.time() - t0
            tag = "EXACT" if is_exact else "UB"
            print(f"    => {sz} gates [{tag}, {method}] ({elapsed:.1f}s)")
            all_results[(fn, N)] = (sz, method, is_exact, elapsed)

    # ============================================================
    # Summary
    # ============================================================
    print("\n\n" + "=" * 70)
    print("RESULTS SUMMARY")
    print("=" * 70)

    print(f"\n{'N':>3} | {'pi(x)mod2':>12} {'type':>6} | {'is_prime':>12} {'type':>6} | {'random':>12} {'type':>6}")
    print("-" * 80)

    pi_s, pr_s, rn_s = {}, {}, {}
    pi_ex = {}

    for N in range(4, 11):
        parts = []
        for fn, dct in [("pi_mod2", pi_s), ("is_prime", pr_s), ("random", rn_s)]:
            s, m, ex, t = all_results.get((fn, N), (-1, "n/a", False, 0))
            tag = "exact" if ex else "UB"
            prefix = "" if ex else "<="
            parts.append(f"{prefix}{s:>10} {tag:>6}")
            if s > 0:
                dct[N] = s
                if fn == "pi_mod2": pi_ex[N] = ex
        print(f"{N:>3} | {' | '.join(parts)}")

    # Growth analysis
    print("\n\nGROWTH RATE ANALYSIS")
    print("=" * 70)

    for label, sd in [("pi(x) mod 2", pi_s), ("is_prime(x)", pr_s), ("random", rn_s)]:
        if len(sd) < 2: continue
        ns = sorted(sd.keys())
        ss = [sd[n] for n in ns]
        print(f"\n{label}: {dict(zip(ns, ss))}")
        for i in range(1, len(ns)):
            if ns[i] == ns[i-1]+1 and ss[i-1] > 0:
                print(f"  s({ns[i]})/s({ns[i-1]}) = {ss[i]/ss[i-1]:.3f}")

        if len(ns) >= 3:
            try:
                import numpy as np
                from scipy.optimize import curve_fit
                x = np.array(ns, dtype=float); y = np.array(ss, dtype=float)
                try:
                    pp, _ = curve_fit(lambda n, a, b: a * n**b, x, y, p0=[1, 2])
                    rp = np.sum((y - pp[0] * x**pp[1])**2)
                    print(f"  Poly fit: {pp[0]:.4f} * N^{pp[1]:.4f}  (SSE={rp:.2f})")
                except: rp = float('inf')
                try:
                    ep, _ = curve_fit(lambda n, a, b: a * b**n, x, y, p0=[0.5, 1.5])
                    re = np.sum((y - ep[0] * ep[1]**x)**2)
                    print(f"  Exp fit:  {ep[0]:.4f} * {ep[1]:.4f}^N  (SSE={re:.2f})")
                except: re = float('inf')
                if rp < float('inf') and re < float('inf'):
                    print(f"  => Better fit: {'POLYNOMIAL' if rp < re else 'EXPONENTIAL'}")
            except ImportError: pass

    print(f"\n\nSHANNON BOUND (random N-var function needs ~2^N/N gates):")
    for N in range(4, 11):
        sh = (1 << N) / N
        s = pi_s.get(N, "?")
        print(f"  N={N}: circuit={s}, Shannon~{sh:.0f}")

    print(f"\nBDD COMPARISON (prior: BDD ~ 2^(0.73*N)):")
    for N in range(4, 11):
        bdd = 2**(0.73*N)
        s = pi_s.get(N, "?")
        print(f"  N={N}: circuit={s}, BDD_est={bdd:.1f}")

    print(f"\nCOMMUNICATION RANK (2^(N/2-1) + 2):")
    for N in range(4, 11):
        cr = 2**(N/2 - 1) + 2
        s = pi_s.get(N, "?")
        print(f"  N={N}: circuit={s}, comm_rank={cr:.0f}")

    print("\nIMPLICATIONS:")
    print("-" * 50)
    if len(pi_s) >= 3:
        ns = sorted(pi_s.keys())
        ss = [pi_s[n] for n in ns]
        ratios = [ss[i]/ss[i-1] for i in range(1, len(ns))
                  if ns[i]==ns[i-1]+1 and ss[i-1]>0]
        if ratios:
            avg = sum(ratios)/len(ratios)
            print(f"  Average growth ratio: {avg:.3f}")
            if avg < 3:
                print("  Growth appears MODERATE -- consistent with polynomial circuits")
            else:
                print("  Growth appears RAPID -- may suggest super-polynomial circuits")
    print("  CAVEAT: synthesis gives UPPER bounds; true minima may be smaller")
    print("  Exact minimization is co-NP-hard for circuits > ~6 gates")

    print("\nDone.")
    return all_results


if __name__ == "__main__":
    results = run_experiment()
