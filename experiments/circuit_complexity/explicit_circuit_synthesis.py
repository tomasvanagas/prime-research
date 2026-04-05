#!/usr/bin/env python3
"""
Explicit Circuit Synthesis for pi(x)
=====================================
For N-bit inputs x (0..2^N-1), compute pi(x) and analyze Boolean circuit complexity.

Measures:
1. BDD sizes with multiple variable orderings per output bit
2. Greedy circuit synthesis (AND/OR/NOT/XOR gates)
3. Exhaustive minimum circuit search for LSB (N <= 7)
4. Threshold gate count estimation
5. Scaling ratios vs N^k

Session 28 experiment.
"""

import sys
import time
import math
import random

# --- Prime counting ---

def sieve_primes(limit):
    if limit < 2:
        return []
    is_prime = bytearray(b'\x01') * (limit + 1)
    is_prime[0] = is_prime[1] = 0
    for i in range(2, int(limit**0.5) + 1):
        if is_prime[i]:
            is_prime[i*i::i] = bytearray(len(is_prime[i*i::i]))
    return [i for i in range(limit + 1) if is_prime[i]]

def compute_pi_table(N):
    limit = (1 << N) - 1
    primes = sieve_primes(limit)
    pi = [0] * (limit + 1)
    idx = 0
    for x in range(limit + 1):
        while idx < len(primes) and primes[idx] <= x:
            idx += 1
        pi[x] = idx
    return pi

def extract_bit(val, bit):
    return (val >> bit) & 1


# --- BDD: use integer bitmask representation ---
# Truth table of n variables = integer of 2^n bits.
# Bit x of the integer = f(x).

def tt_to_int(tt_tuple):
    """Convert tuple truth table to integer bitmask."""
    val = 0
    for i, b in enumerate(tt_tuple):
        if b:
            val |= (1 << i)
    return val

def bdd_size(tt_int, n, var_order):
    """
    Count ROBDD nodes for Boolean function given as integer bitmask,
    with variable ordering var_order (list of bit positions, top-first).

    Uses Shannon decomposition: at each level, split on var_order[level].
    Cofactors are computed by restricting the truth table.
    We memoize on the truth table integer.
    """
    unique_nodes = set()  # set of (level, lo_id, hi_id)
    cache = {}  # tt_int -> node_id
    next_id = [2]  # 0 = False terminal, 1 = True terminal

    def build(tt, level):
        if tt in cache:
            return cache[tt]

        if level == n:
            # Should be 0 or 1 at this point (single entry)
            result = 1 if tt else 0
            cache[tt] = result
            return result

        var = var_order[level]

        # Compute cofactors
        # lo: restrict var=0, hi: restrict var=1
        # lo_tt has bit x set iff tt has bit x set AND bit var of x is 0
        # But we also need to "project out" variable var.
        # Since remaining levels handle remaining variables, we keep full indexing
        # and just rely on memoization — same subfunctions map to same node.

        # Actually: proper cofactor = function with variable fixed.
        # f|_{x_var=0}(x) = f(x with bit var forced to 0)
        # f|_{x_var=1}(x) = f(x with bit var forced to 1)
        # Both are still 2^n-bit integers but effectively depend on n-1 vars.

        lo = 0  # cofactor var=0
        hi = 0  # cofactor var=1
        mask = 1 << var
        for x in range(1 << n):
            if (x & mask) == 0:
                # var=0 position: value goes to lo
                if (tt >> x) & 1:
                    lo |= (1 << x)
                    # Also set the bit with var=1 to same value (project out var)
                    lo |= (1 << (x | mask))
            else:
                # var=1 position: value goes to hi
                if (tt >> x) & 1:
                    hi |= (1 << x)
                    # Also set the bit with var=0 to same value
                    hi |= (1 << (x ^ mask))

        lo_id = build(lo, level + 1)
        hi_id = build(hi, level + 1)

        if lo_id == hi_id:
            cache[tt] = lo_id
            return lo_id

        node = (level, lo_id, hi_id)
        unique_nodes.add(node)
        nid = next_id[0]
        next_id[0] += 1
        cache[tt] = nid
        return nid

    build(tt_int, 0)
    return len(unique_nodes)


def bdd_best_of_orderings(tt_int, n, num_orderings=20):
    """Try multiple variable orderings, return smallest BDD size."""
    best = float('inf')

    orderings = [list(range(n)), list(range(n-1, -1, -1))]

    if n >= 4:
        o = []
        for i in range(n // 2):
            o.append(i)
            o.append(n - 1 - i)
        if n % 2:
            o.append(n // 2)
        orderings.append(o)

    for _ in range(max(0, num_orderings - len(orderings))):
        perm = list(range(n))
        random.shuffle(perm)
        orderings.append(perm)

    for order in orderings:
        t0 = time.time()
        size = bdd_size(tt_int, n, order)
        if size < best:
            best = size
        if time.time() - t0 > 10:
            break  # this ordering class is too slow

    return best


# --- Greedy circuit synthesis ---

def greedy_circuit_size(target_tt, n, max_lib=30000, max_iter=40):
    """
    Build AND/OR/XOR/NOT circuit greedily.
    Returns estimated gate count or -1 if failed.
    Uses integer bitmask representation for speed.
    """
    size = 1 << n
    target = target_tt  # integer bitmask
    ntarget = ((1 << size) - 1) ^ target  # NOT of target

    library = set()
    gate_cost = {}  # tt_int -> gate count

    # Input literals
    for i in range(n):
        tt = 0
        for x in range(size):
            if (x >> i) & 1:
                tt |= (1 << x)
        library.add(tt)
        gate_cost[tt] = 0
        ntt = ((1 << size) - 1) ^ tt
        library.add(ntt)
        gate_cost[ntt] = 1  # 1 NOT gate

    all_zero = 0
    all_one = (1 << size) - 1
    library.add(all_zero)
    library.add(all_one)
    gate_cost[all_zero] = 0
    gate_cost[all_one] = 0

    if target in library:
        return gate_cost[target]
    if ntarget in library:
        return gate_cost[ntarget] + 1

    for iteration in range(max_iter):
        lib_list = list(library)
        if len(lib_list) > 300:
            lib_list = random.sample(lib_list, 300)

        new_fns = {}

        for i in range(len(lib_list)):
            tt1 = lib_list[i]
            c1 = gate_cost.get(tt1, 999)
            for j in range(i, len(lib_list)):
                tt2 = lib_list[j]
                c2 = gate_cost.get(tt2, 999)
                cost = c1 + c2 + 1

                for res in [tt1 & tt2, tt1 | tt2, tt1 ^ tt2]:
                    if res == target:
                        return cost
                    nres = all_one ^ res
                    if nres == target:
                        return cost + 1

                    if res not in library and (res not in new_fns or new_fns[res] > cost):
                        new_fns[res] = cost
                    if nres not in library and (nres not in new_fns or new_fns[nres] > cost + 1):
                        new_fns[nres] = cost + 1

                if len(new_fns) + len(library) > max_lib:
                    break
            if len(new_fns) + len(library) > max_lib:
                break

        if not new_fns:
            break

        for tt, c in new_fns.items():
            library.add(tt)
            if tt not in gate_cost or gate_cost[tt] > c:
                gate_cost[tt] = c

        if target in library:
            return gate_cost[target]
        if ntarget in library:
            return gate_cost[ntarget] + 1

    return gate_cost.get(target, gate_cost.get(ntarget, -1))


# --- Exhaustive minimum circuit search ---

def exhaustive_min_circuit(target_int, n, max_gates=10):
    """
    BFS to find minimum AND/OR/XOR circuit (with free NOT).
    """
    size = 1 << n
    all_one = (1 << size) - 1
    target = target_int
    ntarget = all_one ^ target

    available = set()
    for i in range(n):
        tt = 0
        for x in range(size):
            if (x >> i) & 1:
                tt |= (1 << x)
        available.add(tt)
        available.add(all_one ^ tt)
    available.add(0)
    available.add(all_one)

    if target in available or ntarget in available:
        return 0

    for num_gates in range(1, max_gates + 1):
        new_fns = set()
        all_fns = list(available)

        for i in range(len(all_fns)):
            tt1 = all_fns[i]
            for j in range(i, len(all_fns)):
                tt2 = all_fns[j]
                for res in [tt1 & tt2, tt1 | tt2, tt1 ^ tt2]:
                    if res == target or res == ntarget:
                        return num_gates
                    nres = all_one ^ res
                    if nres == target or nres == ntarget:
                        return num_gates
                    if res not in available:
                        new_fns.add(res)
                    if nres not in available:
                        new_fns.add(nres)

            if len(new_fns) > 500000:
                break

        available.update(new_fns)

        if target in available or ntarget in available:
            return num_gates

        if len(available) > 1000000:
            return -1

    return -1


# --- Threshold gate estimation ---

def estimate_threshold_gates(tt_tuple, n):
    """Estimate threshold gate count from total influence."""
    size = 1 << n
    total_influence = 0.0
    for i in range(n):
        flips = 0
        mask = 1 << i
        for x in range(size):
            x_flip = x ^ mask
            if tt_tuple[x] != tt_tuple[x_flip]:
                flips += 1
        total_influence += flips / size

    # Single linear threshold function: influence <= n
    # Our heuristic: gates ~ max(1, ceil(total_influence / n) * something)
    if total_influence <= n:
        est = 1
    else:
        est = max(1, int(math.ceil(total_influence / 2)))

    return est, total_influence


# --- Main experiment ---

def run_experiment():
    results = {}

    print("=" * 80)
    print("EXPLICIT CIRCUIT SYNTHESIS FOR pi(x)")
    print("=" * 80)
    sys.stdout.flush()

    for N in range(4, 15):
        print(f"\n{'='*60}")
        print(f"N = {N} bits, x in [0, {(1<<N)-1}]")
        print(f"{'='*60}")
        sys.stdout.flush()

        t0 = time.time()

        # 1. Compute pi(x) truth table
        pi_table = compute_pi_table(N)
        max_pi = max(pi_table)
        num_output_bits = max(1, max_pi.bit_length())

        print(f"  pi(2^{N}-1) = {pi_table[-1]}, max = {max_pi}, output bits = {num_output_bits}")

        # Extract per-bit truth tables (both tuple and int forms)
        bit_tts_tuple = []
        bit_tts_int = []
        for b in range(num_output_bits):
            tt_t = tuple(extract_bit(pi_table[x], b) for x in range(1 << N))
            tt_i = tt_to_int(tt_t)
            bit_tts_tuple.append(tt_t)
            bit_tts_int.append(tt_i)

        # 2. BDD analysis
        print(f"\n  --- BDD Analysis ---")
        sys.stdout.flush()
        bdd_sizes = []

        num_ord = 40 if N <= 10 else (15 if N <= 12 else 5)

        skip_remaining = False
        for b in range(num_output_bits):
            if skip_remaining:
                bdd_sizes.append(None)
                continue
            tb = time.time()
            bdd_s = bdd_best_of_orderings(bit_tts_int[b], N, num_ord)
            te = time.time() - tb
            bdd_sizes.append(bdd_s)
            print(f"    Bit {b}: BDD = {bdd_s} nodes ({te:.2f}s)")
            sys.stdout.flush()
            if te > 60:
                print(f"    (remaining bits: timeout)")
                skip_remaining = True

        total_bdd = sum(s for s in bdd_sizes if s is not None)
        print(f"  Total BDD (measured): {total_bdd}")
        sys.stdout.flush()

        # 3. Greedy circuit synthesis
        print(f"\n  --- Greedy Circuit Synthesis ---")
        sys.stdout.flush()
        circuit_sizes = []

        if N <= 11:
            for b in range(num_output_bits):
                ts = time.time()
                cs = greedy_circuit_size(bit_tts_int[b], N)
                te = time.time() - ts
                circuit_sizes.append(cs)
                print(f"    Bit {b}: ~{cs} gates ({te:.2f}s)")
                sys.stdout.flush()
                if te > 60:
                    for b2 in range(b + 1, num_output_bits):
                        circuit_sizes.append(None)
                    break
        else:
            circuit_sizes = [None] * num_output_bits
            print("    (skipped, N > 11)")

        # 4. LSB analysis
        print(f"\n  --- LSB (pi mod 2) ---")
        lsb_tt = bit_tts_tuple[0]
        ones = sum(lsb_tt)
        print(f"    Density: {ones}/{1<<N} = {ones/(1<<N):.4f}")

        lsb_min_circuit = None
        if N <= 7:
            te_start = time.time()
            lsb_min_circuit = exhaustive_min_circuit(bit_tts_int[0], N, max_gates=10)
            te = time.time() - te_start
            print(f"    Exhaustive min circuit: {lsb_min_circuit} gates ({te:.2f}s)")
        elif N == 8:
            te_start = time.time()
            lsb_min_circuit = exhaustive_min_circuit(bit_tts_int[0], N, max_gates=5)
            te = time.time() - te_start
            print(f"    Exhaustive min circuit (depth<=5): {lsb_min_circuit} gates ({te:.2f}s)")
        sys.stdout.flush()

        # 5. Threshold estimation
        print(f"\n  --- Threshold Gates ---")
        threshold_est = []
        for b in range(min(num_output_bits, 4)):
            est, infl = estimate_threshold_gates(bit_tts_tuple[b], N)
            threshold_est.append(est)
            print(f"    Bit {b}: ~{est} threshold gates (influence={infl:.2f})")

        # 6. Scaling
        lsb_bdd = bdd_sizes[0] if bdd_sizes and bdd_sizes[0] is not None else None
        if lsb_bdd is not None and lsb_bdd > 0:
            print(f"\n  --- Scaling ---")
            print(f"    BDD/2^N      = {lsb_bdd / (1<<N):.6f}")
            print(f"    BDD/2^(N/2)  = {lsb_bdd / (2**(N/2)):.4f}")
            for k in [1, 2, 3, 4]:
                print(f"    BDD/N^{k}      = {lsb_bdd / (N**k):.4f}")
            print(f"    log2(BDD)    = {math.log2(lsb_bdd):.2f}")

        elapsed = time.time() - t0
        print(f"\n  Total time: {elapsed:.2f}s")
        sys.stdout.flush()

        results[N] = {
            'max_pi': max_pi,
            'num_output_bits': num_output_bits,
            'bdd_sizes': bdd_sizes,
            'total_bdd': total_bdd,
            'circuit_sizes': circuit_sizes,
            'lsb_min_circuit': lsb_min_circuit,
            'threshold_est': threshold_est,
            'lsb_density': ones / (1 << N),
            'elapsed': elapsed,
        }

        if elapsed > 300:
            print(f"\n  ** Stopping: N={N} took {elapsed:.0f}s **")
            break

    # --- Summary ---
    print_summary(results)
    return results


def print_summary(results):
    print("\n\n" + "=" * 80)
    print("SUMMARY")
    print("=" * 80)

    header = f"{'N':>3} | {'pi_max':>7} | {'bits':>4} | {'LSB_BDD':>8} | {'Tot_BDD':>8} | {'LSB_circ':>9} | {'LSB_min':>7} | {'Thresh':>6} | {'log2BDD':>8} | {'BDD/N^2':>8} | {'BDD/sqrt':>9}"
    print(header)
    print("-" * len(header))

    ns = []
    bdd_vals = []
    for N in sorted(results.keys()):
        r = results[N]
        lb = r['bdd_sizes'][0] if r['bdd_sizes'] and r['bdd_sizes'][0] is not None else None
        lc = r['circuit_sizes'][0] if r['circuit_sizes'] and r['circuit_sizes'][0] is not None else None
        lm = r['lsb_min_circuit']
        th = r['threshold_est'][0] if r['threshold_est'] else None

        if lb is not None and lb > 0:
            ns.append(N)
            bdd_vals.append(lb)
            l2 = f"{math.log2(lb):.2f}"
            rn2 = f"{lb / N**2:.2f}"
            rsq = f"{lb / (2**(N/2)):.4f}"
        else:
            l2 = rn2 = rsq = '-'

        print(f"{N:>3} | {r['max_pi']:>7} | {r['num_output_bits']:>4} | {str(lb or '-'):>8} | {r['total_bdd']:>8} | {str(lc if lc is not None else '-'):>9} | {str(lm if lm is not None else '-'):>7} | {str(th if th is not None else '-'):>6} | {l2:>8} | {rn2:>8} | {rsq:>9}")

    if len(ns) >= 4:
        print("\n\nSCALING FITS (LSB BDD):")

        log_ns = [math.log(n) for n in ns]
        log_bs = [math.log(b) for b in bdd_vals]
        n_pts = len(ns)
        sx = sum(log_ns); sy = sum(log_bs)
        sxy = sum(x*y for x,y in zip(log_ns, log_bs))
        sx2 = sum(x*x for x in log_ns)
        denom = n_pts * sx2 - sx**2
        if denom != 0:
            alpha = (n_pts * sxy - sx * sy) / denom
            c_pow = math.exp((sy - alpha * sx) / n_pts)
            print(f"  Power law: BDD ~ {c_pow:.2f} * N^{alpha:.2f}")

        log2_bs = [math.log2(b) for b in bdd_vals]
        sx_l = sum(ns); sy_l = sum(log2_bs)
        sxy_l = sum(x*y for x,y in zip(ns, log2_bs))
        sx2_l = sum(x*x for x in ns)
        denom_l = n_pts * sx2_l - sx_l**2
        if denom_l != 0:
            beta = (n_pts * sxy_l - sx_l * sy_l) / denom_l
            gamma = (sy_l - beta * sx_l) / n_pts
            print(f"  Exponential: BDD ~ 2^({beta:.4f}*N + {gamma:.2f})")
            print(f"  (sqrt(x) = 2^(N/2) => beta = 0.5)")
            print(f"  (poly(N) => beta ~ 0)")

            if beta > 0.4:
                print(f"  ==> BDD scales ~sqrt(x), consistent with sqrt barrier")
            elif beta > 0.15:
                print(f"  ==> Sub-sqrt, super-poly: 2^({beta:.2f}*N)")
            else:
                print(f"  ==> Polynomial scaling in N!")

            print(f"\n  NOTE: BDD is a restricted model; circuits can be exponentially smaller.")

    write_results(results, ns, bdd_vals)
    print(f"\nDone.")


def write_results(results, ns, bdd_vals):
    lines = []
    lines.append("# Explicit Circuit Synthesis for pi(x) -- Results")
    lines.append("")
    lines.append("## Experiment (Session 28)")
    lines.append("")
    lines.append("For N-bit inputs x in [0, 2^N - 1], compute pi(x) and measure Boolean circuit complexity.")
    lines.append("This is the FIRST explicit circuit synthesis for pi(x) in this project.")
    lines.append("")

    lines.append("### Summary Table")
    lines.append("")
    lines.append("| N | pi(2^N-1) | Bits | LSB BDD | Total BDD | LSB circuit | LSB min exact | Threshold | log2(BDD) | BDD/N^2 | BDD/2^(N/2) |")
    lines.append("|---|-----------|------|---------|-----------|-------------|---------------|-----------|-----------|---------|-------------|")

    for N in sorted(results.keys()):
        r = results[N]
        lb = r['bdd_sizes'][0] if r['bdd_sizes'] and r['bdd_sizes'][0] is not None else '-'
        lc = r['circuit_sizes'][0] if r['circuit_sizes'] and r['circuit_sizes'][0] is not None else '-'
        lm = r['lsb_min_circuit'] if r['lsb_min_circuit'] is not None else '-'
        th = r['threshold_est'][0] if r['threshold_est'] else '-'

        if isinstance(lb, int) and lb > 0:
            l2 = f"{math.log2(lb):.2f}"
            rn2 = f"{lb / N**2:.2f}"
            rsq = f"{lb / 2**(N/2):.4f}"
        else:
            l2 = rn2 = rsq = '-'

        lines.append(f"| {N} | {r['max_pi']} | {r['num_output_bits']} | {lb} | {r['total_bdd']} | {lc} | {lm} | {th} | {l2} | {rn2} | {rsq} |")

    lines.append("")
    lines.append("### Per-bit BDD sizes")
    lines.append("")
    for N in sorted(results.keys()):
        r = results[N]
        sizes = [str(s) if s is not None else '?' for s in r['bdd_sizes']]
        lines.append(f"- N={N}: bits={r['num_output_bits']}, BDD = [{', '.join(sizes)}]")

    lines.append("")
    lines.append("### Scaling Analysis")
    lines.append("")

    if len(ns) >= 4:
        lines.append("| N | LSB BDD | log2(BDD) | BDD/N^2 | BDD/N^3 | BDD/2^(N/2) | BDD/2^(0.8N) |")
        lines.append("|---|---------|-----------|---------|---------|-------------|--------------|")
        for i, N in enumerate(ns):
            b = bdd_vals[i]
            l2 = math.log2(b) if b > 0 else 0
            lines.append(f"| {N} | {b} | {l2:.2f} | {b/N**2:.2f} | {b/N**3:.4f} | {b/2**(N/2):.6f} | {b/2**(0.8*N):.6f} |")

        # Fits
        n_pts = len(ns)
        log_ns = [math.log(n) for n in ns]
        log_bs = [math.log(b) for b in bdd_vals]
        sx = sum(log_ns); sy = sum(log_bs)
        sxy = sum(x*y for x,y in zip(log_ns, log_bs))
        sx2 = sum(x*x for x in log_ns)
        denom = n_pts * sx2 - sx**2

        log2_bs = [math.log2(b) for b in bdd_vals]
        sx_l = sum(ns); sy_l = sum(log2_bs)
        sxy_l = sum(x*y for x,y in zip(ns, log2_bs))
        sx2_l = sum(x*x for x in ns)
        denom_l = n_pts * sx2_l - sx_l**2

        if denom != 0 and denom_l != 0:
            alpha = (n_pts * sxy - sx * sy) / denom
            c_pow = math.exp((sy - alpha * sx) / n_pts)
            beta = (n_pts * sxy_l - sx_l * sy_l) / denom_l
            gamma = (sy_l - beta * sx_l) / n_pts

            lines.append("")
            lines.append("### Fitted Models")
            lines.append("")
            lines.append(f"- **Power law:** BDD ~ {c_pow:.2f} * N^{alpha:.2f}")
            lines.append(f"- **Exponential:** BDD ~ 2^({beta:.4f}*N + {gamma:.2f})")
            lines.append(f"  - sqrt(x) = 2^(N/2) corresponds to beta = 0.5")
            lines.append(f"  - poly(N) corresponds to beta ~ 0")
            lines.append(f"  - **Measured beta = {beta:.4f}**")
            lines.append("")
            lines.append("### Interpretation")
            lines.append("")

            if beta > 0.4:
                lines.append("**The LSB BDD scales approximately as 2^(N/2) = sqrt(x).** This is consistent")
                lines.append("with the universal sqrt barrier observed across all other computational models.")
            elif beta > 0.15:
                lines.append(f"**The LSB BDD scales as 2^({beta:.2f}*N) — sub-sqrt but super-polynomial.**")
            else:
                lines.append("**The LSB BDD appears to scale polynomially in N!**")

    lines.append("")
    lines.append("### Important Caveats")
    lines.append("")
    lines.append("1. **BDD != Circuit.** BDDs (branching programs) are restricted; general circuits can be exponentially more compact.")
    lines.append("2. **BDD upper bound IS a circuit upper bound.** BDD size S => circuit size O(S).")
    lines.append("3. **Variable ordering matters.** Our best-of-k is an upper bound; true min-BDD could be smaller.")
    lines.append("4. **N = 4..14 is small.** Scaling may not extrapolate to N = 30+ (x ~ 10^9).")
    lines.append("")
    lines.append("### Relation to Previous Results")
    lines.append("")
    lines.append("| Model | Scaling | Source |")
    lines.append("|-------|---------|--------|")
    lines.append("| ANF (GF2) | degree N, ~50% nonzero coeffs | Session 13 |")
    lines.append("| OBDD | 2^{0.79*N} | Session 20 |")
    lines.append("| Comm rank | 2^{N/2-1}+2 | Session 17 |")
    if len(ns) >= 4 and denom_l != 0:
        lines.append(f"| BDD (multi-order) | ~2^({beta:.2f}*N) | This experiment |")
    lines.append("")

    outpath = "/apps/aplikacijos/prime-research/experiments/circuit_complexity/explicit_circuit_synthesis_results.md"
    with open(outpath, 'w') as f:
        f.write('\n'.join(lines) + '\n')
    print(f"\nResults written to {outpath}")


if __name__ == '__main__':
    random.seed(42)
    run_experiment()
