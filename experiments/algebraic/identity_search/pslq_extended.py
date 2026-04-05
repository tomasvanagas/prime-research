#!/usr/bin/env python3
"""
Extended PSLQ identity search for f(x) = pi(x) - R(x) at x up to 100000.

Prior work (pslq_identity_results.txt) tested up to x=10000 and found NO valid
relations -- all PSLQ hits were spurious (coeff of f = 0) or point-specific.

This experiment extends the search to x=100000 scale and adds:
  - Functional relations: f(ax) vs f(x) for a=2,3,4
  - Shifted relations: f(x), f(x+1), ..., f(x+10)
  - Full basis with zeta zero oscillations (gamma1, gamma2)

Date: 2026-04-05
"""

import numpy as np
from mpmath import mp, mpf, log, sqrt, li, sin, cos, pslq, cbrt, pi
import time
import sys

mp.dps = 50

# First two zeta zero imaginary parts (high precision)
gamma1 = mpf('14.134725141734693790457251983562470270784257115699')
gamma2 = mpf('21.022039638771554992628479593896902777334340524903')

# Load data
data = np.load('experiments/algebraic/identity_search/fx_data.npz')
x_arr = data['x']      # 2..100000
f_arr = data['f']       # f(x) = pi(x) - R(x)
pi_arr = data['pi']     # pi(x)

def get_f(x_val):
    """Get f(x) from the data array."""
    idx = int(x_val) - 2  # x_arr[0] = 2
    if 0 <= idx < len(f_arr):
        return mpf(str(f_arr[idx]))
    return None

def get_pi(x_val):
    """Get pi(x) from data."""
    idx = int(x_val) - 2
    if 0 <= idx < len(pi_arr):
        return int(pi_arr[idx])
    return None

results = []

def record(section, desc, relation_str, residual, coeff_f, valid):
    """Record a result."""
    results.append({
        'section': section,
        'desc': desc,
        'relation': relation_str,
        'residual': float(residual) if residual is not None else None,
        'coeff_f': coeff_f,
        'valid': valid
    })
    flag = " *** VALID RELATION ***" if valid else ""
    status = "SPURIOUS(coeff_f=0)" if coeff_f == 0 and residual is not None and float(residual) < 1e-10 else ""
    print(f"  [{section}] {desc}: residual={residual}{flag}{' '+status if status else ''}")

# ============================================================
# SECTION A: Standard PSLQ at large x values
# Basis: [f(x), 1, log(x), sqrt(x), x^{1/3}, x^{1/4}, x^{1/5},
#          li(x), li(sqrt(x)), li(x^{1/3}),
#          sin(gamma1*log(x)), cos(gamma1*log(x)),
#          sin(gamma2*log(x)), cos(gamma2*log(x))]
# ============================================================
print("=" * 70)
print("SECTION A: Standard PSLQ at large x values")
print("=" * 70)

test_points_A = [10000, 20000, 50000, 100000]
basis_names = ['f(x)', '1', 'log(x)', 'sqrt(x)', 'x^{1/3}', 'x^{1/4}', 'x^{1/5}',
               'li(x)', 'li(sqrt(x))', 'li(x^{1/3})',
               'sin(g1*log(x))', 'cos(g1*log(x))', 'sin(g2*log(x))', 'cos(g2*log(x))']

for xv in test_points_A:
    fx = get_f(xv)
    if fx is None:
        print(f"  x={xv}: data not available, skipping")
        continue

    xm = mpf(xv)
    logx = log(xm)
    sqrtx = sqrt(xm)
    cbrtx = xm ** (mpf(1)/3)
    x14 = xm ** (mpf(1)/4)
    x15 = xm ** (mpf(1)/5)
    lix = li(xm)
    lisqrt = li(sqrtx)
    licbrt = li(cbrtx)
    sg1 = sin(gamma1 * logx)
    cg1 = cos(gamma1 * logx)
    sg2 = sin(gamma2 * logx)
    cg2 = cos(gamma2 * logx)

    vec = [fx, mpf(1), logx, sqrtx, cbrtx, x14, x15, lix, lisqrt, licbrt, sg1, cg1, sg2, cg2]

    print(f"\n  x={xv}, f(x)={float(fx):.6f}")
    t0 = time.time()
    try:
        rel = pslq(vec, maxcoeff=10**6, maxsteps=5000)
    except Exception as e:
        print(f"    PSLQ failed: {e}")
        record('A-PSLQ', f'x={xv}', 'FAILED', None, None, False)
        continue
    elapsed = time.time() - t0

    if rel is None:
        print(f"    No relation found ({elapsed:.1f}s)")
        record('A-PSLQ', f'x={xv}', 'No relation', None, None, False)
    else:
        # Compute residual
        residual = abs(sum(mpf(c) * v for c, v in zip(rel, vec)))
        coeff_f = rel[0]
        terms = [f"{rel[i]}*{basis_names[i]}" for i in range(len(rel)) if rel[i] != 0]
        rel_str = " + ".join(terms)
        valid = (coeff_f != 0) and (float(residual) < 1e-10)

        print(f"    Relation ({elapsed:.1f}s): {rel_str}")
        print(f"    Residual: {residual}")
        print(f"    coeff(f) = {coeff_f}")
        record('A-PSLQ', f'x={xv}', rel_str, residual, coeff_f, valid)

        # Cross-validate at nearby point if valid
        if valid and xv < 100000:
            xv2 = xv + 1000
            fx2 = get_f(xv2)
            xm2 = mpf(xv2)
            logx2 = log(xm2)
            vec2 = [fx2, mpf(1), logx2, sqrt(xm2), xm2**(mpf(1)/3), xm2**(mpf(1)/4),
                     xm2**(mpf(1)/5), li(xm2), li(sqrt(xm2)), li(xm2**(mpf(1)/3)),
                     sin(gamma1*logx2), cos(gamma1*logx2), sin(gamma2*logx2), cos(gamma2*logx2)]
            res2 = abs(sum(mpf(c)*v for c, v in zip(rel, vec2)))
            print(f"    Cross-check at x={xv2}: residual={res2}")
            if float(res2) > 1e-5:
                print(f"    >> FAILS cross-validation (point-specific)")

# ============================================================
# SECTION B: Functional relations f(ax) vs f(x)
# Basis: [f(x), f(2x), f(3x), f(4x), 1, log(x)]
# ============================================================
print("\n" + "=" * 70)
print("SECTION B: Functional relations f(ax) vs f(x)")
print("=" * 70)

test_points_B = [1000, 2000, 5000, 10000]
func_names = ['f(x)', 'f(2x)', 'f(3x)', 'f(4x)', '1', 'log(x)']

for xv in test_points_B:
    fx = get_f(xv)
    f2x = get_f(2 * xv)
    f3x = get_f(3 * xv)
    f4x = get_f(4 * xv)

    if any(v is None for v in [fx, f2x, f3x, f4x]):
        print(f"  x={xv}: some values not available, skipping")
        continue

    xm = mpf(xv)
    vec = [fx, f2x, f3x, f4x, mpf(1), log(xm)]

    print(f"\n  x={xv}: f(x)={float(fx):.4f}, f(2x)={float(f2x):.4f}, f(3x)={float(f3x):.4f}, f(4x)={float(f4x):.4f}")
    t0 = time.time()
    try:
        rel = pslq(vec, maxcoeff=10**6, maxsteps=5000)
    except Exception as e:
        print(f"    PSLQ failed: {e}")
        record('B-FUNC', f'x={xv}', 'FAILED', None, None, False)
        continue
    elapsed = time.time() - t0

    if rel is None:
        print(f"    No relation found ({elapsed:.1f}s)")
        record('B-FUNC', f'x={xv}', 'No relation', None, None, False)
    else:
        residual = abs(sum(mpf(c) * v for c, v in zip(rel, vec)))
        coeff_f = rel[0]
        terms = [f"{rel[i]}*{func_names[i]}" for i in range(len(rel)) if rel[i] != 0]
        rel_str = " + ".join(terms)
        # For functional: valid if ANY f-coefficient is nonzero
        has_f_coeff = any(rel[i] != 0 for i in range(4))
        valid = has_f_coeff and (float(residual) < 1e-10)

        print(f"    Relation ({elapsed:.1f}s): {rel_str}")
        print(f"    Residual: {residual}")
        record('B-FUNC', f'x={xv}', rel_str, residual, coeff_f, valid)

        # Cross-validate
        if valid:
            xv2 = xv + 500 if xv + 500 * 4 <= 100000 else xv + 100
            fx2 = get_f(xv2)
            f2x2 = get_f(2 * xv2)
            f3x2 = get_f(3 * xv2)
            f4x2 = get_f(4 * xv2)
            if all(v is not None for v in [fx2, f2x2, f3x2, f4x2]):
                vec2 = [fx2, f2x2, f3x2, f4x2, mpf(1), log(mpf(xv2))]
                res2 = abs(sum(mpf(c) * v for c, v in zip(rel, vec2)))
                print(f"    Cross-check at x={xv2}: residual={res2}")
                if float(res2) > 1e-5:
                    print(f"    >> FAILS cross-validation (point-specific)")

# Also test extended functional basis with sqrt terms
print("\n  --- Extended functional basis with sqrt(x) ---")
func_ext_names = ['f(x)', 'f(2x)', 'f(3x)', 'f(4x)', '1', 'log(x)', 'sqrt(x)', '1/log(x)']

for xv in [5000, 10000, 20000]:
    fx = get_f(xv)
    f2x = get_f(2 * xv)
    f3x = get_f(3 * xv)
    f4x = get_f(4 * xv)

    if any(v is None for v in [fx, f2x, f3x, f4x]):
        continue

    xm = mpf(xv)
    vec = [fx, f2x, f3x, f4x, mpf(1), log(xm), sqrt(xm), 1/log(xm)]

    print(f"\n  x={xv} (extended):")
    t0 = time.time()
    try:
        rel = pslq(vec, maxcoeff=10**6, maxsteps=5000)
    except Exception as e:
        print(f"    PSLQ failed: {e}")
        record('B-FUNC-EXT', f'x={xv}', 'FAILED', None, None, False)
        continue
    elapsed = time.time() - t0

    if rel is None:
        print(f"    No relation found ({elapsed:.1f}s)")
        record('B-FUNC-EXT', f'x={xv}', 'No relation', None, None, False)
    else:
        residual = abs(sum(mpf(c) * v for c, v in zip(rel, vec)))
        has_f_coeff = any(rel[i] != 0 for i in range(4))
        valid = has_f_coeff and (float(residual) < 1e-10)
        terms = [f"{rel[i]}*{func_ext_names[i]}" for i in range(len(rel)) if rel[i] != 0]
        rel_str = " + ".join(terms)
        print(f"    Relation ({elapsed:.1f}s): {rel_str}")
        print(f"    Residual: {residual}")
        record('B-FUNC-EXT', f'x={xv}', rel_str, residual, rel[0], valid)

        if valid:
            xv2 = xv + 1000 if (xv + 1000) * 4 <= 100000 else xv + 200
            fvals2 = [get_f(a * xv2) for a in [1, 2, 3, 4]]
            if all(v is not None for v in fvals2):
                xm2 = mpf(xv2)
                vec2 = fvals2 + [mpf(1), log(xm2), sqrt(xm2), 1/log(xm2)]
                res2 = abs(sum(mpf(c) * v for c, v in zip(rel, vec2)))
                print(f"    Cross-check at x={xv2}: residual={res2}")
                if float(res2) > 1e-5:
                    print(f"    >> FAILS cross-validation")

# ============================================================
# SECTION C: Shifted relations f(x), f(x+1), ..., f(x+10), 1
# ============================================================
print("\n" + "=" * 70)
print("SECTION C: Shifted relations f(x), f(x+1), ..., f(x+10)")
print("=" * 70)

test_points_C = [1000, 5000, 10000]
shift_range = 11  # f(x) through f(x+10)

for xv in test_points_C:
    fvals = [get_f(xv + k) for k in range(shift_range)]
    if any(v is None for v in fvals):
        print(f"  x={xv}: data not available, skipping")
        continue

    vec = fvals + [mpf(1)]
    shift_names = [f'f(x+{k})' for k in range(shift_range)] + ['1']

    print(f"\n  x={xv}: f(x)={float(fvals[0]):.4f}, ..., f(x+10)={float(fvals[10]):.4f}")
    t0 = time.time()
    try:
        rel = pslq(vec, maxcoeff=10**6, maxsteps=5000)
    except Exception as e:
        print(f"    PSLQ failed: {e}")
        record('C-SHIFT', f'x={xv}', 'FAILED', None, None, False)
        continue
    elapsed = time.time() - t0

    if rel is None:
        print(f"    No relation found ({elapsed:.1f}s)")
        record('C-SHIFT', f'x={xv}', 'No relation', None, None, False)
    else:
        residual = abs(sum(mpf(c) * v for c, v in zip(rel, vec)))
        has_f_coeff = any(rel[i] != 0 for i in range(shift_range))
        valid = has_f_coeff and (float(residual) < 1e-10)
        terms = [f"{rel[i]}*{shift_names[i]}" for i in range(len(rel)) if rel[i] != 0]
        rel_str = " + ".join(terms)
        print(f"    Relation ({elapsed:.1f}s): {rel_str}")
        print(f"    Residual: {residual}")
        record('C-SHIFT', f'x={xv}', rel_str, residual, rel[0], valid)

        # Cross-validate at different x
        if valid:
            xv2 = xv + 2000 if xv + 2000 + 10 <= 100000 else xv + 500
            fvals2 = [get_f(xv2 + k) for k in range(shift_range)]
            if all(v is not None for v in fvals2):
                vec2 = fvals2 + [mpf(1)]
                res2 = abs(sum(mpf(c) * v for c, v in zip(rel, vec2)))
                print(f"    Cross-check at x={xv2}: residual={res2}")
                if float(res2) > 1e-5:
                    print(f"    >> FAILS cross-validation (point-specific)")

# Also test shorter shift windows for tighter relations
print("\n  --- Short shift window: f(x), f(x+1), f(x+2), f(x+3), 1 ---")
for xv in [1000, 5000, 10000, 50000]:
    fvals = [get_f(xv + k) for k in range(4)]
    if any(v is None for v in fvals):
        continue
    vec = fvals + [mpf(1)]
    names = [f'f(x+{k})' for k in range(4)] + ['1']
    t0 = time.time()
    try:
        rel = pslq(vec, maxcoeff=10**6, maxsteps=5000)
    except Exception as e:
        record('C-SHIFT-SHORT', f'x={xv}', 'FAILED', None, None, False)
        continue
    elapsed = time.time() - t0

    if rel is None:
        print(f"  x={xv}: No relation ({elapsed:.1f}s)")
        record('C-SHIFT-SHORT', f'x={xv}', 'No relation', None, None, False)
    else:
        residual = abs(sum(mpf(c) * v for c, v in zip(rel, vec)))
        has_f = any(rel[i] != 0 for i in range(4))
        valid = has_f and float(residual) < 1e-10
        terms = [f"{rel[i]}*{names[i]}" for i in range(len(rel)) if rel[i] != 0]
        rel_str = " + ".join(terms)
        print(f"  x={xv}: {rel_str} | residual={residual}")
        record('C-SHIFT-SHORT', f'x={xv}', rel_str, residual, rel[0], valid)

        if valid:
            xv2 = xv + 3000 if xv + 3000 + 3 <= 100000 else xv + 300
            fvals2 = [get_f(xv2 + k) for k in range(4)]
            if all(v is not None for v in fvals2):
                vec2 = fvals2 + [mpf(1)]
                res2 = abs(sum(mpf(c) * v for c, v in zip(rel, vec2)))
                print(f"    Cross-check x={xv2}: residual={res2}")
                if float(res2) > 1e-5:
                    print(f"    >> FAILS cross-validation")

# ============================================================
# SUMMARY
# ============================================================
print("\n" + "=" * 70)
print("SUMMARY OF ALL RESULTS")
print("=" * 70)

valid_count = sum(1 for r in results if r['valid'])
total = len(results)
print(f"\nTotal tests: {total}")
print(f"Valid relations (coeff_f != 0, residual < 1e-10): {valid_count}")

print(f"\n{'Rank':<6}{'Section':<16}{'Residual':<20}{'Description'}")
print("-" * 90)

# Sort: valid first, then by residual
sorted_results = sorted(results, key=lambda r: (not r['valid'], r['residual'] if r['residual'] is not None else 1e99))

for i, r in enumerate(sorted_results):
    res_str = f"{r['residual']:.6e}" if r['residual'] is not None else "N/A"
    flag = " <<<VALID>>>" if r['valid'] else ""
    print(f"{i+1:<6}{r['section']:<16}{res_str:<20}{r['desc']}{flag}")

if valid_count == 0:
    print("\n*** NO VALID ALGEBRAIC RELATIONS FOUND ***")
    print("All PSLQ hits are either spurious (coeff_f=0) or point-specific.")
    print("This confirms f(x) = pi(x) - R(x) has no simple algebraic identity")
    print("in the tested bases up to x=100000.")
else:
    print(f"\n*** {valid_count} POTENTIALLY VALID RELATION(S) FOUND -- CHECK CROSS-VALIDATION ***")

print("\nDone.")
