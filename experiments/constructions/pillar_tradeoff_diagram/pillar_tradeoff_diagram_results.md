# C6 — Three-Pillars × HKM Tradeoff Diagram — Results

**Status:** BUILT (S81). All four pre-stated falsification criteria PASS.
**Edges composed:** E7.7 (three-pillars informational closure) + E6.7
(HKM at `(N^{8/15}, N^{1/3})`).
**Verdict:** B-grade structural classification (refinement, not novelty).
No polylog opening. New content is the explicit pillar-by-pillar Pareto
placement, not the individual algorithms.

---

## 1. The catalog

14 algorithms tagged across 3 pillars (full table in
`pillar_tradeoff_diagram_results.json`). Each has `(alpha, beta)` such
that `T = O~(N^alpha)` and `S = O~(N^beta)` (polylog absorbed).

| Pillar           | # algorithms |
|------------------|-------------:|
| prime_positions  | 4 |
| zeta_zeros       | 4 |
| floor_values     | 6 |

The catalog is intentionally a **superset** of the project's EDGES.md
algorithms — naive sieve, segmented sieve, Aggarwal binary search,
direct explicit-formula evaluation, Galway 2004, Odlyzko-Schoenhage,
Connes one-shot fit, Lucy DP, Lucy + Fenwick, Meissel-Lehmer, Gourdon,
primecount, HKM. Each entry's `(alpha, beta)` traces back to the cited
source.

---

## 2. Per-pillar Pareto frontiers

```
Pillar: floor_values
  alpha=8/15=0.5333  beta=1/3=0.3333  T*S=13/15=0.8667  : HKM       <-- HKM

Pillar: prime_positions
  alpha=1/2=0.5000   beta=1/2=0.5000  T*S=1.0000        : Aggarwal binary search

Pillar: zeta_zeros
  alpha=1/2=0.5000   beta=1/2=0.5000  T*S=1.0000        : Direct explicit-formula
  alpha=1/2=0.5000   beta=1/2=0.5000  T*S=1.0000        : Galway 2004
  alpha=1/2=0.5000   beta=1/2=0.5000  T*S=1.0000        : Odlyzko-Schoenhage
  alpha=1/2=0.5000   beta=1/2=0.5000  T*S=1.0000        : Connes one-shot
```

**Floor-pillar frontier collapses to HKM.** All other floor-pillar
algorithms (Lucy DP, Lucy+Fenwick, Meissel-Lehmer, Gourdon, primecount)
are dominated elementwise by HKM. (Meissel-Lehmer at `(2/3, 1/3)` is
dominated by HKM at `(8/15, 1/3)` in time, same space; primecount at
`(2/3, 1/2)` is dominated in both axes.)

**Zeta-pillar frontier is a 4-way tie at `(1/2, 1/2)`.** All
zero-pillar entries land at the same exponent point — they differ only
in polylog factors and conditional assumptions (Galway needs GRH;
direct explicit-formula needs explicit zero density).

**Prime-pillar frontier collapses to Aggarwal at `(1/2, 1/2)`.** Naive
and segmented sieves are dominated. Note: Aggarwal's exponent inherits
from the chosen `pi(x)` sub-routine; the cited `(1/2, 1/2)` reflects the
zero-pillar-coupled use; coupling Aggarwal to HKM would yield
`(8/15, 1/3)` and re-locate Aggarwal-as-meta-algorithm onto the floor
pillar (see §6 for the meta-algorithm note).

---

## 3. Pillar dominance summary

| pillar           | min alpha (time) | min beta (space) | min T*S product |
|------------------|-----------------:|-----------------:|----------------:|
| floor_values     | **0.5333**       | **0.3333**       | **0.8667**      |
| prime_positions  | 0.5000           | 0.5000           | 1.0000          |
| zeta_zeros       | 0.5000           | 0.5000           | 1.0000          |

Three structural claims fall out:

- **Time-only minimum** is achieved on the prime *and* zeta pillars
  (tie at `alpha = 1/2`). The floor pillar is `0.0333` worse on time
  alone (`8/15 ≈ 0.5333`).
- **Space-only minimum** is uniquely achieved on the floor pillar at
  `beta = 1/3`. No zero-pillar or prime-pillar algorithm in the catalog
  achieves `beta < 1/2`.
- **T*S minimum** is uniquely achieved on the floor pillar at
  `T*S = 13/15`. Both other pillars are flat at `T*S = 1`.

---

## 4. HKM placement (the headline structural result)

```
HKM: alpha = 8/15 ≈ 0.5333    beta = 1/3 ≈ 0.3333    T*S = 13/15 ≈ 0.8667
Pillar by definition: floor_values   (NTT on Dirichlet convolution)
On floor_values Pareto frontier? YES (it IS the frontier)
Cross-pillar dominators of HKM:  NONE
```

**Structural claim:**

> HKM's `(alpha, beta) = (8/15, 1/3)` point is achievable on the
> *floor pillar only*. No catalog algorithm on the prime pillar or
> the zeta pillar simultaneously achieves both `T <= N^{8/15}` and
> `S <= N^{1/3}`.

This is the C6-specific composition output: **the cross-pillar
dominance of HKM is empty**, so HKM occupies a region of the
`(T, S)` plane that the other two pillars cannot reach.

---

## 5. E7.6 saturation (floor pillar)

E7.6 (Meissel-Lehmer DAG pebbling) gives the lower bound
`T*S >= Omega(N^{5/6} / log N)`, valid on the floor pillar.

| algorithm                            | T*S exp | gap to 5/6 |
|--------------------------------------|--------:|-----------:|
| HKM                                  | 13/15 ≈ 0.8667 | +0.0333 |
| Meissel-Lehmer (classical)           | 1.0000  | +0.1667    |
| Lucy + Fenwick                       | 7/6 ≈ 1.1667 | +0.3333 |
| Gourdon variant                      | 7/6 ≈ 1.1667 | +0.3333 |
| primecount production                | 7/6 ≈ 1.1667 | +0.3333 |
| Lucy DP basic                        | 5/4 = 1.2500 | +0.4167 |

**HKM saturates E7.6 to within `N^{0.034}`.** The floor pillar is
asymptotically tight: any algorithm beating HKM by a polynomial factor
on this pillar would have to circumvent the Lucy-DAG pebbling lower
bound, which is a black-pebble-game theorem and does not depend on the
specific algorithm structure beyond DAG width × depth.

---

## 6. Tradeoff diagram (ASCII)

```
    space (beta)
    1.0  +---------------------+
    1.00 | ....................P |
         | ..................... |
         | ..................... |
         | ..................... |
    0.80 | ..................... |
         | ..................... |
         | ..................... |
         | ..................... |
    0.60 | ..................... |
         | ..................... |
         | ..........*..F.F....P |
         | ..................... |
    0.40 | ..................... |
         | ...........H.F....... |
         | ..................... |
         | ..................... |
    0.20 | ..................... |
         | ..................... |
         | ..................... |
         | ..................... |
    0.00 | ..................... |
         +---------------------+
           0.0     0.5     1.0
              time (alpha)

Legend:
  P = prime_positions only      Z = zeta_zeros only
  F = floor_values only         H = HKM (floor)
  * = 2+ pillars at this cell   . = no algorithm here
```

Reading the diagram:

- `P` at top-right `(1.0, 1.0)` — naive sieve.
- `P` at right `(1.0, 0.5)` — segmented and wheel sieves (overlapping).
- `*` at center `(0.5, 0.5)` — Aggarwal (prime) overlapping with the
  four zero-pillar algorithms (direct explicit-formula, Galway,
  Odlyzko-Schoenhage, Connes one-shot fit).
- `F` at `(0.667, 0.5)` — Lucy+Fenwick / Gourdon / primecount.
- `F` at `(0.667, 0.333)` — Meissel-Lehmer classical.
- `F` at `(0.75, 0.5)` — Lucy DP basic.
- `H` at `(0.5333, 0.333)` — **HKM, unique to the floor pillar**.

The H cell is the *only* cell in the lower-left quadrant `alpha <= 0.6,
beta <= 0.4`. No other pillar enters that quadrant.

---

## 7. Falsification check (pre-stated)

| ID | criterion                                                    | result |
|----|--------------------------------------------------------------|--------|
| F1 | Some catalog entry beats HKM elementwise                     | PASS — no violator |
| F2 | Non-floor entry has `(alpha, beta) <=` HKM elementwise       | PASS — no violator |
| F3 | All three pillars have identical Pareto frontiers            | PASS — frontiers differ |
| F4 | Floor-pillar entry achieves `T*S < N^{5/6}`, contradicts E7.6 | PASS — no violator |

**Construction passes all four falsifiers.**

---

## 8. Subtleties

### Aggarwal as a meta-algorithm

Aggarwal binary search (E6.6) makes `O(log n)` calls to `pi(x)` plus a
BPSW primality check. Its time exponent inherits from the chosen
`pi(x)` sub-routine. The catalog lists Aggarwal at `(1/2, 1/2)` as a
fixed point representing zero-pillar-coupled use; coupling Aggarwal to
HKM gives Aggarwal-via-HKM at `(8/15, 1/3)`, which migrates the
"effective placement" of Aggarwal into the floor pillar. So Aggarwal's
pillar placement is OUTPUT pillar (prime), while its WORK pillar
follows whichever `pi(x)` it is composed with. This is a natural
property of the Aggarwal reduction (E6.6 is best read as
"`p(n) -> O(log n)` queries to `pi(x)`" rather than as a standalone
algorithm).

### "alpha = 1/2 + epsilon" vs "alpha = 1/2"

The zero-pillar algorithms achieve `O(N^{1/2+epsilon})` time/space
under epsilon-control via known explicit zero-density estimates.
For exponent comparison we round `1/2+epsilon -> 1/2` since the
`epsilon` is a polylog factor. This rounding does NOT affect any of
the F1-F4 conclusions: the gap between zero-pillar `(1/2+eps, 1/2+eps)`
and HKM `(8/15, 1/3)` is `Omega(N^{1/15})` on space, which dwarfs any
sub-polynomial epsilon.

### Conditional assumptions

The catalog mixes unconditional, BPSW-conditional, and GRH-conditional
algorithms. Pareto frontier comparisons are honest only WITHIN a
conditional level (e.g., HKM is unconditional; Galway's `(1/2+eps,
1/2+eps)` requires GRH). Restricting to **unconditional** algorithms
only:

| pillar           | best (T*S) (unconditional) |
|------------------|---------------------------:|
| floor_values     | HKM at 13/15 ≈ 0.8667 |
| prime_positions  | Aggarwal at 1.0 (BPSW required for the primality check, so technically conditional) |
| zeta_zeros       | Direct explicit-formula at 1.0 (with explicit zero-density input) |

Restricting further to **fully unconditional** algorithms with no
primality oracle and no zero-density input drops the prime and zeta
pillars to their next-best (~`(1, 1)` for naive/segmented and
non-existent for zeta-pillar), and HKM still wins on T*S.

---

## 9. What this construction does and does not deliver

### Delivered

- **A pillar-by-pillar Pareto frontier** for the project's catalog of
  pi(x) algorithms, in `(alpha, beta)` exponent coordinates.
- **An explicit cross-pillar dominance check at HKM's point**:
  cross-pillar dominators of HKM are empty.
- **A floor-pillar saturation table** showing HKM is within
  `N^{0.034}` of the E7.6 lower bound.
- **A reproducible JSON / Python catalog** (`pillar_tradeoff_diagram.py`)
  that future agents can extend by adding new entries to `CATALOG`.

### NOT delivered (and why)

- **No polylog opening.** The construction is a structural
  classification, not a new algorithm. Floor pillar is asymptotically
  tight; zero pillar is at `1/2+eps` on both axes; prime pillar
  inherits whichever pillar is used for `pi(x)`. No path in the
  diagram leads to `(polylog, polylog)`.
- **No new lower bound.** E7.6 (Meissel-Lehmer pebbling) is the
  governing lower bound on the floor pillar. The construction does
  not extend E7.6 to other pillars; doing so would require a different
  argument (e.g., zero-counting capacity for the zeta pillar).
- **No proof that the pillar list (E7.7) is exhaustive at the
  algorithmic level.** E7.7 is informational closure; this
  construction only PROBES the algorithmic projection of E7.7 onto
  the (T, S) plane. A fourth pillar with `(T, S)` outside all three
  current frontiers would be a major project event; this construction
  does not preclude it.

---

## 10. Closure outcome

Mode: **DUPLICATE-PLUS** (extends E7.7 + E6.7 with structural placement
data; not a new edge by itself; refines the project's understanding of
how the three pillars project onto the algorithmic Pareto plane).

A new EDGES.md entry **is not warranted** — the structural placement
follows from existing edges. EDGES.md `E6.7` and `E7.7` are annotated
with a back-pointer to this construction.

CLOSED_PATHS row added: `S81-C6 — three-pillars Pareto classification:
HKM uniquely accesses (8/15, 1/3); zero pillar leads on time, floor
pillar leads on space and T*S; floor pillar saturates E7.6 to N^{0.034}`.
