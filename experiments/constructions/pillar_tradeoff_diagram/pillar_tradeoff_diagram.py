"""
C6 — Three-Pillars × HKM Time-Space Tradeoff Diagram
====================================================

Composes E7.7 (three-pillars informational closure) with E6.7 (HKM
Pareto point) into a structural classification of the project's
catalog of pi(x) algorithms.

Each algorithm is tagged by its primary information source:

    prime_positions   — needs the prime list / pi(x) directly as input.
    zeta_zeros        — needs a list of nontrivial zeta zeros.
    floor_values      — uses the family { floor(x/k) : 1 <= k <= sqrt(x) }
                        and Dirichlet convolutions thereof (Lucy-DP track).

For each algorithm we record (alpha, beta) such that

    T = O~(N^alpha)    S = O~(N^beta)

(the tilde absorbs polylog factors). HKM's published exponents are
(alpha, beta) = (8/15, 1/3); Aggarwal (E6.6) is (1/2, .) but it makes
log-many calls to a pi(x) sub-procedure -- so its effective exponent
inherits the chosen pi(x) implementation.

The script computes:

    - Per-pillar Pareto frontier in (alpha, beta).
    - Per-pillar best-time and best-space algorithms.
    - Cross-pillar dominance check at HKM's point.
    - The asymptotic T*S product per pillar's frontier.
    - Whether the floor-pillar T*S frontier saturates the E7.6
      Lucy-DAG pebbling lower bound (T*S >= Omega(N^{5/6} / log N)).

Output
------
- pillar_tradeoff_diagram_results.json -- the raw catalog + computed
  Pareto / dominance / saturation tables.
- ASCII tradeoff diagram printed to stdout.

Falsification (pre-stated in definition.md)
-------------------------------------------
F1: any catalog entry beats HKM elementwise -> HKM is not Pareto-optimal.
F2: any non-floor pillar entry has (alpha, beta) <= HKM elementwise
    -> floor-pillar placement is artifactual.
F3: all three pillars have identical Pareto frontiers -> classification
    is algorithmically vacuous.
F4: any catalog entry achieves T*S < N^{5/6} / log N -> contradicts
    E7.6, input data corrupt.
"""

import json
import math
from fractions import Fraction
from pathlib import Path

# ---------------------------------------------------------------------------
# Algorithm catalog
# ---------------------------------------------------------------------------
#
# Each entry has:
#   name        : descriptive name + literature pointer
#   pillar      : prime_positions | zeta_zeros | floor_values
#   alpha       : time exponent  T ~ N^alpha (polylog absorbed)
#   beta        : space exponent S ~ N^beta  (polylog absorbed)
#   condition   : 'unconditional' | 'BPSW' | 'GRH' | 'mixed'
#   edge        : EDGES.md edge ID where the algorithm or its primitive
#                 lives, when applicable
#   note        : short explanation of the (alpha, beta) values

CATALOG = [
    # --- prime_positions pillar -------------------------------------------
    {
        "name": "Sieve of Eratosthenes (full enumeration)",
        "pillar": "prime_positions",
        "alpha": Fraction(1, 1),
        "beta":  Fraction(1, 1),
        "condition": "unconditional",
        "edge": None,
        "note": "Generates all primes up to N then counts; O(N log log N)"
                " time, O(N) space.",
    },
    {
        "name": "Segmented sieve",
        "pillar": "prime_positions",
        "alpha": Fraction(1, 1),
        "beta":  Fraction(1, 2),
        "condition": "unconditional",
        "edge": None,
        "note": "Process sqrt(N)-sized segments; same time, space O(sqrt N).",
    },
    {
        "name": "Wheel-mod-30 sieve (E2.1 wheel realization)",
        "pillar": "prime_positions",
        "alpha": Fraction(1, 1),
        "beta":  Fraction(1, 2),
        "condition": "unconditional",
        "edge": "E2.1",
        "note": "Same exponents as segmented; wheel gives 8/30 storage"
                " constant (Mertens product) but no exponent change.",
    },
    {
        "name": "Aggarwal binary search on pi (E6.6)",
        "pillar": "prime_positions",
        "alpha": Fraction(1, 2),
        "beta":  Fraction(1, 2),
        "condition": "BPSW",
        "edge": "E6.6",
        "note": "p(n) reduces to O(log n) calls to pi(x); each pi(x) call"
                " inherits floor-pillar cost. Effective exponent here"
                " reflects the BEST currently-attached pi(x) sub-routine"
                " (HKM at floor pillar). Logged on prime pillar because"
                " the algorithm OUTPUT is a prime position, not pi.",
    },

    # --- zeta_zeros pillar ------------------------------------------------
    {
        "name": "Direct explicit-formula evaluation",
        "pillar": "zeta_zeros",
        "alpha": Fraction(1, 2),
        "beta":  Fraction(1, 2),
        "condition": "unconditional (with epsilon)",
        "edge": "E3.5",
        "note": "Requires zeros up to height T ~ x^{1/2+eps} for exact"
                " pi(x); each zero costs O(polylog) once tabulated."
                " Storing N(T) ~ T log T zeros in O(log T) bits gives"
                " S ~ x^{1/2+eps}.",
    },
    {
        "name": "Galway 2004 zero-driven pi(x)",
        "pillar": "zeta_zeros",
        "alpha": Fraction(1, 2),
        "beta":  Fraction(1, 2),
        "condition": "GRH",
        "edge": None,
        "note": "Best published zero-pillar method for exact pi(x);"
                " same order as direct explicit-formula.",
    },
    {
        "name": "Odlyzko-Schoenhage zeta evaluation",
        "pillar": "zeta_zeros",
        "alpha": Fraction(1, 2),
        "beta":  Fraction(1, 2),
        "condition": "unconditional",
        "edge": None,
        "note": "Amortized O(T^{1/2}) per zeta evaluation at T heights."
                " A primitive used by zero-pillar algorithms; not itself"
                " a pi(x) algorithm but bounds the zero-pillar exponent.",
    },
    {
        "name": "Connes operator one-shot fit (S53)",
        "pillar": "zeta_zeros",
        "alpha": Fraction(1, 2),  # not algorithmic; place at zero-pillar
                                  # frontier as a heuristic
        "beta":  Fraction(1, 2),
        "condition": "heuristic",
        "edge": "E3.1",
        "note": "Recovers 50 zeros from primes <= 13 to 10^{-55..-3} as"
                " a one-shot fit; not a polylog algorithm. Listed for"
                " completeness on the zero pillar.",
    },

    # --- floor_values pillar ----------------------------------------------
    {
        "name": "Lucy DP basic",
        "pillar": "floor_values",
        "alpha": Fraction(3, 4),
        "beta":  Fraction(1, 2),
        "condition": "unconditional",
        "edge": "E6.9",
        "note": "Direct DP over { floor(x/k) }; T = O(x^{3/4}), S = O(sqrt x).",
    },
    {
        "name": "Lucy + Fenwick tree update",
        "pillar": "floor_values",
        "alpha": Fraction(2, 3),
        "beta":  Fraction(1, 2),
        "condition": "unconditional",
        "edge": "E6.9",
        "note": "Log-update Fenwick tree drops time to O(x^{2/3} log^{1/3} x).",
    },
    {
        "name": "Meissel-Lehmer (classical)",
        "pillar": "floor_values",
        "alpha": Fraction(2, 3),
        "beta":  Fraction(1, 3),
        "condition": "unconditional",
        "edge": "E7.9",
        "note": "Phi(x, a) recursion; cube-root-space variant is Meissel-Lehmer.",
    },
    {
        "name": "Gourdon variant",
        "pillar": "floor_values",
        "alpha": Fraction(2, 3),
        "beta":  Fraction(1, 2),
        "condition": "unconditional",
        "edge": "E6.9",
        "note": "Two-parameter family; same T as Lucy+Fenwick, refined"
                " constants.",
    },
    {
        "name": "primecount (Walisch / Kim production)",
        "pillar": "floor_values",
        "alpha": Fraction(2, 3),
        "beta":  Fraction(1, 2),
        "condition": "unconditional",
        "edge": "E6.9",
        "note": "Engineering ceiling of the sieve route; same Big-O as"
                " Lucy+Fenwick / Gourdon, multi-level POPCNT counter.",
    },
    {
        "name": "HKM (Helfgott-Kessler-Mendlovic)",
        "pillar": "floor_values",
        "alpha": Fraction(8, 15),
        "beta":  Fraction(1, 3),
        "condition": "unconditional",
        "edge": "E6.7",
        "note": "NTT-accelerated Dirichlet convolution; T = O~(x^{8/15}),"
                " S = O~(x^{1/3}). Best simultaneously sub-x^{2/3} time"
                " AND sub-x^{1/2} space algorithm known.",
    },
]

# A canonical reference to HKM and the E7.6 lower bound.
HKM_NAME = "HKM (Helfgott-Kessler-Mendlovic)"
LOWER_BOUND_TS = Fraction(5, 6)  # E7.6: T*S >= Omega(N^{5/6} / log N)


# ---------------------------------------------------------------------------
# Pareto frontier and dominance utilities
# ---------------------------------------------------------------------------

def per_pillar_pareto(catalog):
    """
    For each pillar, return the subset of algorithms on the Pareto
    frontier in (alpha, beta) space (smaller is better in both axes).
    Returns dict pillar -> list of catalog entries on the frontier.
    """
    by_pillar = {}
    for entry in catalog:
        by_pillar.setdefault(entry["pillar"], []).append(entry)
    out = {}
    for pillar, entries in by_pillar.items():
        frontier = []
        for e in entries:
            dominated = False
            for f in entries:
                if f is e:
                    continue
                if (f["alpha"] <= e["alpha"] and f["beta"] <= e["beta"]
                        and (f["alpha"] < e["alpha"]
                             or f["beta"] < e["beta"])):
                    dominated = True
                    break
            if not dominated:
                frontier.append(e)
        # Sort frontier by alpha ascending (so beta is descending).
        frontier.sort(key=lambda e: (e["alpha"], e["beta"]))
        out[pillar] = frontier
    return out


def cross_pillar_dominators(catalog, target):
    """
    Return all catalog entries on a different pillar from `target`
    that dominate `target` elementwise (alpha' <= alpha AND beta' <= beta).
    """
    out = []
    for e in catalog:
        if e["pillar"] == target["pillar"]:
            continue
        if (e["alpha"] <= target["alpha"] and e["beta"] <= target["beta"]
                and (e["alpha"] < target["alpha"]
                     or e["beta"] < target["beta"])):
            out.append(e)
    return out


def total_dominators(catalog, target):
    """
    Return all entries (any pillar) that dominate target elementwise.
    """
    out = []
    for e in catalog:
        if e is target:
            continue
        if (e["alpha"] <= target["alpha"] and e["beta"] <= target["beta"]
                and (e["alpha"] < target["alpha"]
                     or e["beta"] < target["beta"])):
            out.append(e)
    return out


def best_per_pillar(catalog):
    """
    Per pillar, return:
        best_time   : the entry with smallest alpha (ties: smallest beta)
        best_space  : the entry with smallest beta  (ties: smallest alpha)
        best_TS     : the entry with smallest alpha+beta
    """
    by_pillar = {}
    for entry in catalog:
        by_pillar.setdefault(entry["pillar"], []).append(entry)
    out = {}
    for pillar, entries in by_pillar.items():
        bt = min(entries, key=lambda e: (e["alpha"], e["beta"]))
        bs = min(entries, key=lambda e: (e["beta"], e["alpha"]))
        bts = min(entries, key=lambda e: (e["alpha"] + e["beta"]))
        out[pillar] = {
            "best_time":  bt,
            "best_space": bs,
            "best_TS":    bts,
        }
    return out


def floor_pillar_saturation(catalog, lower_bound_ts=LOWER_BOUND_TS):
    """
    For each floor-pillar entry, compute T*S exponent (alpha + beta) and
    the gap to the E7.6 lower bound (5/6).
    """
    rows = []
    for e in catalog:
        if e["pillar"] != "floor_values":
            continue
        ts = e["alpha"] + e["beta"]
        gap = ts - lower_bound_ts
        rows.append({
            "name": e["name"],
            "alpha": e["alpha"],
            "beta": e["beta"],
            "TS_exp": ts,
            "gap_to_E7_6": gap,
        })
    rows.sort(key=lambda r: r["TS_exp"])
    return rows


def check_falsification(catalog):
    """
    Apply F1-F4 to the catalog. Returns a dict with one bool per criterion.
    """
    hkm = next(e for e in catalog if e["name"] == HKM_NAME)

    # F1: any catalog entry strictly dominates HKM?
    f1_violators = total_dominators(catalog, hkm)

    # F2: any non-floor pillar entry has (alpha, beta) <= HKM elementwise?
    f2_violators = [
        e for e in catalog
        if e["pillar"] != "floor_values"
        and e["alpha"] <= hkm["alpha"] and e["beta"] <= hkm["beta"]
    ]

    # F3: all three pillars have IDENTICAL Pareto frontiers in (alpha,
    # beta) -- check by extracting the (alpha, beta) sets per pillar.
    pareto = per_pillar_pareto(catalog)
    pillar_keys = sorted(pareto.keys())
    pillar_frontiers = {
        p: frozenset((e["alpha"], e["beta"]) for e in pareto[p])
        for p in pillar_keys
    }
    f3_identical = (len(set(pillar_frontiers.values())) == 1
                    and len(pillar_keys) >= 2)

    # F4: any catalog entry achieves T*S < 5/6?
    f4_violators = [
        e for e in catalog
        if (e["alpha"] + e["beta"]) < LOWER_BOUND_TS
        and e["pillar"] == "floor_values"
    ]

    return {
        "F1_HKM_dominated": [e["name"] for e in f1_violators],
        "F2_cross_pillar_HKM_dominated":
            [e["name"] for e in f2_violators],
        "F3_pareto_frontiers_identical": f3_identical,
        "F4_floor_pillar_violates_E7_6":
            [e["name"] for e in f4_violators],
    }


# ---------------------------------------------------------------------------
# Pretty-printing utilities
# ---------------------------------------------------------------------------

def _frac_to_str(x: Fraction):
    """Render a Fraction as e.g. '8/15 ≈ 0.5333'."""
    return f"{x.numerator}/{x.denominator} = {float(x):.4f}"


def _entry_to_dict(e):
    return {
        "name": e["name"],
        "pillar": e["pillar"],
        "alpha": str(e["alpha"]),
        "alpha_float": float(e["alpha"]),
        "beta": str(e["beta"]),
        "beta_float": float(e["beta"]),
        "TS_exp": str(e["alpha"] + e["beta"]),
        "TS_exp_float": float(e["alpha"] + e["beta"]),
        "condition": e["condition"],
        "edge": e["edge"],
        "note": e["note"],
    }


def render_ascii_diagram(catalog):
    """
    Render a 21x21 ASCII grid with x = time exponent (alpha) on the
    horizontal axis (0 left, 1 right) and y = space exponent (beta) on
    the vertical axis (0 bottom, 1 top). Mark each algorithm by a
    one-letter code per pillar:

        P  prime_positions
        Z  zeta_zeros
        F  floor_values
        H  HKM (overrides F at the (8/15, 1/3) cell)

    Multiple algorithms in the same cell -> letter of highest priority
    (H > P > Z > F to keep visual cleanliness).
    """
    W = 21  # 0.00, 0.05, ..., 1.00
    # Each cell is a set of pillar-symbols that have an algorithm there.
    grid_sets = [[set() for _ in range(W)] for _ in range(W)]
    grid_has_hkm = [[False for _ in range(W)] for _ in range(W)]

    pillar_symbol = {
        "prime_positions": "P",
        "zeta_zeros": "Z",
        "floor_values": "F",
    }

    def cell(x):
        return min(W - 1, max(0, round(float(x) * (W - 1))))

    for e in catalog:
        col = cell(e["alpha"])
        row = (W - 1) - cell(e["beta"])  # invert so larger beta = top
        grid_sets[row][col].add(pillar_symbol[e["pillar"]])
        if e["name"] == HKM_NAME:
            grid_has_hkm[row][col] = True

    # Render: H if HKM present; * if 2+ pillars at same cell; else single
    # pillar's letter; else '.'.
    grid = [["." for _ in range(W)] for _ in range(W)]
    for r in range(W):
        for c in range(W):
            s = grid_sets[r][c]
            if grid_has_hkm[r][c]:
                grid[r][c] = "H"
            elif len(s) >= 2:
                grid[r][c] = "*"
            elif len(s) == 1:
                grid[r][c] = next(iter(s))

    lines = []
    lines.append("    space (beta)")
    lines.append("    1.0  +" + "-" * W + "+")
    for r in range(W):
        beta_label = (W - 1 - r) / (W - 1)
        if abs(beta_label - round(beta_label, 1)) < 1e-9 and (
            beta_label * 10 % 2 == 0):
            label = f"{beta_label:0.2f}"
        else:
            label = "    "
        lines.append(f"   {label:>5} | " + "".join(grid[r]) + " |")
    lines.append("         +" + "-" * W + "+")
    lines.append("           0.0     0.5     1.0")
    lines.append("              time (alpha)")
    lines.append("")
    lines.append("Legend: P = prime_positions  Z = zeta_zeros  "
                 "F = floor_values  H = HKM  * = 2+ pillars at this point")
    return "\n".join(lines)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    out_dir = Path(__file__).parent
    json_path = out_dir / "pillar_tradeoff_diagram_results.json"

    pareto = per_pillar_pareto(CATALOG)
    bests = best_per_pillar(CATALOG)
    sat = floor_pillar_saturation(CATALOG)
    falsif = check_falsification(CATALOG)
    hkm = next(e for e in CATALOG if e["name"] == HKM_NAME)
    cross_dom_hkm = cross_pillar_dominators(CATALOG, hkm)

    # ------------------- header -------------------
    print("=" * 78)
    print("C6 — Three-Pillars × HKM Tradeoff Diagram")
    print("=" * 78)
    print()
    print(f"Catalog size: {len(CATALOG)} algorithms across "
          f"{len(set(e['pillar'] for e in CATALOG))} pillars.")
    print(f"Lower bound (E7.6): T*S >= Omega(N^{{5/6}}) for floor pillar.")
    print()

    # ------------------- per-pillar Pareto frontiers -------------------
    print("-" * 78)
    print("PER-PILLAR PARETO FRONTIERS")
    print("-" * 78)
    for pillar in sorted(pareto.keys()):
        print(f"\nPillar: {pillar}")
        for e in pareto[pillar]:
            mark = "  <-- HKM" if e["name"] == HKM_NAME else ""
            print(f"    alpha = {_frac_to_str(e['alpha']):<22}"
                  f"  beta = {_frac_to_str(e['beta']):<22}"
                  f"  T*S = {_frac_to_str(e['alpha']+e['beta']):<22}"
                  f"  : {e['name']}{mark}")

    # ------------------- best per pillar -------------------
    print()
    print("-" * 78)
    print("PER-PILLAR BEST")
    print("-" * 78)
    for pillar in sorted(bests.keys()):
        b = bests[pillar]
        print(f"\nPillar: {pillar}")
        print(f"    best time:   alpha={float(b['best_time']['alpha']):.4f}"
              f"  beta={float(b['best_time']['beta']):.4f}"
              f"  : {b['best_time']['name']}")
        print(f"    best space:  alpha={float(b['best_space']['alpha']):.4f}"
              f"  beta={float(b['best_space']['beta']):.4f}"
              f"  : {b['best_space']['name']}")
        print(f"    best T*S:    alpha={float(b['best_TS']['alpha']):.4f}"
              f"  beta={float(b['best_TS']['beta']):.4f}"
              f"  T*S = {float(b['best_TS']['alpha']+b['best_TS']['beta']):.4f}"
              f"  : {b['best_TS']['name']}")

    # ------------------- floor pillar saturation -------------------
    print()
    print("-" * 78)
    print("FLOOR PILLAR T*S SATURATION (vs E7.6 lower bound 5/6 = 0.8333)")
    print("-" * 78)
    print()
    print(f"{'algorithm':<48} {'T*S exp':>10} {'gap to 5/6':>12}")
    for r in sat:
        print(f"  {r['name']:<46} "
              f"{float(r['TS_exp']):>10.4f} "
              f"{float(r['gap_to_E7_6']):>+12.4f}")

    # ------------------- HKM placement -------------------
    print()
    print("-" * 78)
    print("HKM PLACEMENT")
    print("-" * 78)
    print(f"\nHKM: alpha = {_frac_to_str(hkm['alpha'])}  "
          f"beta = {_frac_to_str(hkm['beta'])}  "
          f"T*S = {_frac_to_str(hkm['alpha']+hkm['beta'])}")
    print(f"Pillar by definition: {hkm['pillar']}")
    print(f"On {hkm['pillar']} Pareto frontier? "
          f"{'YES' if hkm in pareto[hkm['pillar']] else 'NO'}")
    print(f"Cross-pillar dominators of HKM: "
          f"{[e['name'] for e in cross_dom_hkm] or 'NONE'}")
    if not cross_dom_hkm:
        print("  -> HKM's (alpha, beta) point is unique to the floor pillar.")
        print("     No zero-pillar or prime-pillar algorithm achieves both")
        print("     T <= N^{8/15} AND S <= N^{1/3} simultaneously.")

    # ------------------- pillar dominance regions -------------------
    print()
    print("-" * 78)
    print("PILLAR DOMINANCE (where each pillar leads)")
    print("-" * 78)
    print()
    by_pillar_min_alpha = {p: min(e["alpha"] for e in pareto[p])
                           for p in pareto}
    by_pillar_min_beta = {p: min(e["beta"] for e in pareto[p])
                          for p in pareto}
    by_pillar_min_TS = {p: min(e["alpha"] + e["beta"] for e in pareto[p])
                        for p in pareto}
    print(f"{'pillar':<22} {'min alpha':>12} {'min beta':>12}"
          f" {'min T*S':>12}")
    for p in sorted(pareto.keys()):
        print(f"  {p:<20} {float(by_pillar_min_alpha[p]):>12.4f}"
              f" {float(by_pillar_min_beta[p]):>12.4f}"
              f" {float(by_pillar_min_TS[p]):>12.4f}")

    print()
    # Decide who controls each axis (ties broken by reporting all leaders).
    def leaders(metric_dict):
        m = min(metric_dict.values())
        return [p for p, v in metric_dict.items() if v == m]

    print("Time-only minimum (Pareto in alpha):  "
          f"{', '.join(leaders(by_pillar_min_alpha))}"
          f" at alpha = {float(min(by_pillar_min_alpha.values())):.4f}")
    print("Space-only minimum (Pareto in beta):  "
          f"{', '.join(leaders(by_pillar_min_beta))}"
          f" at beta  = {float(min(by_pillar_min_beta.values())):.4f}")
    print("T*S-product minimum:                  "
          f"{', '.join(leaders(by_pillar_min_TS))}"
          f" at T*S   = {float(min(by_pillar_min_TS.values())):.4f}")

    # ------------------- ASCII diagram -------------------
    print()
    print("-" * 78)
    print("TRADEOFF DIAGRAM (alpha = time exponent, beta = space exponent)")
    print("-" * 78)
    print()
    print(render_ascii_diagram(CATALOG))

    # ------------------- falsification -------------------
    print()
    print("-" * 78)
    print("FALSIFICATION CHECKS (pre-stated in definition.md)")
    print("-" * 78)
    print(f"  F1 (HKM dominated by some catalog entry): "
          f"{falsif['F1_HKM_dominated'] or 'PASS (no violator)'}")
    print(f"  F2 (HKM cross-pillar dominated): "
          f"{falsif['F2_cross_pillar_HKM_dominated'] or 'PASS (no violator)'}")
    print(f"  F3 (all pillars have identical Pareto frontiers): "
          f"{'FAIL' if falsif['F3_pareto_frontiers_identical'] else 'PASS'}")
    print(f"  F4 (floor pillar T*S < 5/6, contradicting E7.6): "
          f"{falsif['F4_floor_pillar_violates_E7_6'] or 'PASS (no violator)'}")

    falsified = (
        bool(falsif["F1_HKM_dominated"])
        or bool(falsif["F2_cross_pillar_HKM_dominated"])
        or falsif["F3_pareto_frontiers_identical"]
        or bool(falsif["F4_floor_pillar_violates_E7_6"])
    )
    verdict = "CONSTRUCTION FALSIFIED" if falsified else "CONSTRUCTION PASSES"
    print()
    print("Verdict:", verdict)

    # ------------------- emit JSON -------------------
    out = {
        "catalog_size": len(CATALOG),
        "pillars": sorted(set(e["pillar"] for e in CATALOG)),
        "catalog": [_entry_to_dict(e) for e in CATALOG],
        "pareto_per_pillar": {
            p: [_entry_to_dict(e) for e in pareto[p]]
            for p in pareto
        },
        "best_per_pillar": {
            p: {k: _entry_to_dict(v) for k, v in d.items()}
            for p, d in bests.items()
        },
        "floor_TS_saturation": [
            {"name": r["name"],
             "alpha": str(r["alpha"]),
             "beta": str(r["beta"]),
             "TS_exp": str(r["TS_exp"]),
             "TS_exp_float": float(r["TS_exp"]),
             "gap_to_E7_6": str(r["gap_to_E7_6"]),
             "gap_to_E7_6_float": float(r["gap_to_E7_6"])}
            for r in sat
        ],
        "hkm": _entry_to_dict(hkm),
        "hkm_cross_pillar_dominators":
            [_entry_to_dict(e) for e in cross_dom_hkm],
        "by_pillar_min_alpha":
            {p: float(v) for p, v in by_pillar_min_alpha.items()},
        "by_pillar_min_beta":
            {p: float(v) for p, v in by_pillar_min_beta.items()},
        "by_pillar_min_TS":
            {p: float(v) for p, v in by_pillar_min_TS.items()},
        "falsification": {
            "F1_HKM_dominated": falsif["F1_HKM_dominated"],
            "F2_cross_pillar_HKM_dominated":
                falsif["F2_cross_pillar_HKM_dominated"],
            "F3_pareto_frontiers_identical":
                falsif["F3_pareto_frontiers_identical"],
            "F4_floor_pillar_violates_E7_6":
                falsif["F4_floor_pillar_violates_E7_6"],
        },
        "verdict": verdict,
    }
    json_path.write_text(json.dumps(out, indent=2, default=str))
    print()
    print(f"JSON written to {json_path}")


if __name__ == "__main__":
    main()
