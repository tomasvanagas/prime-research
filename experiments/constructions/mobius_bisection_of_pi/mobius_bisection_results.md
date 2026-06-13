# Möbius bisection of π(x) — results

**Session.** S427 (paradigm-shift mode, no external technique imported).

**Edges composed.** E2.2 (Liouville bisection `π = (x − L)/2 − C_3`),
E1.6 / E2.10 (parity bisection + `L(x) mod 2 = x mod 2` trap), folklore
Möbius/squarefree decomposition `M(x) = S_e(x) − S_o(x)`,
`Q(x) = S_e(x) + S_o(x)`.

**No external technique imported.** No new cross-domain entry; no
WebFetch; no new ATTACK_VECTORS entries. Pure recombination of edges
already in the project plus elementary Dirichlet-series identities.

## What was built

A single Python script (`mobius_bisection.py`, no helper variants) that
sieves `ω, Ω, μ, λ, 1_sqf` on `[1, N]`, computes the cumulative
`L(x), M(x), Q(x), C_3(x), C_3*(x), NS_o(x), NS_e(x), S_o(x), S_e(x)`
and `π(x)`, and tests six pre-stated falsifiers (`definition.md`).

## Verification at N = 10⁶ (wall 6.27 s)

Aggregate counts (from `mobius_bisection_results.json`):

| symbol     | value at N = 10⁶ |
|------------|------------------|
| `π(N)`     | 78498            |
| `L(N)`     | −530             |
| `M(N)`     | 212              |
| `Q(N)`     | 607926           |
| `C_3(N)`   | 421767           |
| `C_3*(N)`  | 225359           |
| `NS_o(N)`  | 196408           |
| `NS_e(N)`  | 195666           |
| `S_o(N)`   | 303857           |
| `S_e(N)`   | 304069           |

**All six falsifiers PASS for every x ∈ [1, 10⁶].** The result file's
`failures` dict is empty; `failures = {}`.

* **F1** `π(x) = (Q(x) − M(x))/2 − C_3*(x)` — exact at every x.
  Spot check at N: `(607926 − 212)/2 − 225359 = 303857 − 225359 = 78498
  = π(N)`. ✓
* **F2** `NS_o(x) = (x − Q − L + M)/2` — exact at every x.
  Spot check at N: `(10⁶ − 607926 + 530 + 212)/2 = 392816/2 = 196408
  = NS_o(N)`. ✓
* **F3** `C_3(x) − C_3*(x) = NS_o(x)` — exact at every x.
  Spot check at N: `421767 − 225359 = 196408 = NS_o(N)`. ✓
* **F4** parity `((x − L)/2 ⊕ (Q − M)/2) ≡ NS_o (mod 2)` — exact at
  every x.
* **F5** parity `π(x) ≡ ((Q − M)/2) ⊕ C_3*(x)  (mod 2)` — exact at
  every x.
* **F6** sanity check `π(x) = (x − L)/2 − C_3(x)` (E2.2) — exact at
  every x. ✓

These verifications also held integer-exact at smaller N = 5·10⁴ (single
test pre-flight); falsifier behaviour is x-independent (each is verified
exhaustively for all x in the tested range).

## Empirical pseudorandomness — parity bits at N = 10⁶

Frequencies-of-1 of the parity sequences over x ∈ [1, N]:

| sequence    | freq₁    | lag-1 ac  | density-predicted lag-1   |
|-------------|----------|-----------|---------------------------|
| A par       | 0.50039  | −0.00053  | ≈ 1 − 2·(1/2) = 0         |
| B par       | 0.49934  | +0.39228  | ≈ 1 − 2·0.304 = 0.392     |
| C₃ par      | 0.49989  | +0.15647  | ≈ 1 − 2·0.422 = 0.156     |
| C₃* par     | 0.49813  | +0.54928  | ≈ 1 − 2·0.225 = 0.550     |
| NS_o par    | 0.49939  | +0.60718  | ≈ 1 − 2·0.196 = 0.608     |
| π par       | 0.50200  | +0.84300  | ≈ 1 − 2·(π/N · 1) ≈ 0.843 |

**Every lag-1 autocorrelation matches its density-only prediction to
within 1 %** (the standard `1 − 2·(increment-density)` law for a
two-state cumulative-counter parity). The parity sequences are
**density-structured but not predictability-structured** — exactly the
state of E1.6 transposed to the Möbius side.

Mutual information (in bits, computed over 10⁶ x-values):

| pair                   | MI (bits)  |
|------------------------|------------|
| `I(A par ; B par)`     | 1.09 × 10⁻⁶ |
| `I(C₃ par ; C₃* par)`  | 1.09 × 10⁻⁶ |
| `I(A par ; C₃ par)`    | 1.15 × 10⁻⁵ |
| `I(B par ; C₃* par)`   | 1.15 × 10⁻⁵ |
| `I(π par ; NS_o par)`  | 2.89 × 10⁻⁷ |

`I(B ; C₃*) ≈ 1.15 × 10⁻⁵ bits`, the same order as `I(A ; C₃)` from
E1.6. **The Möbius-side parity bisection inherits E1.6's statistical
independence**: B(x) mod 2 and C₃*(x) mod 2 are mutually
near-independent at finite x.

## Net new project content

1. **The Möbius bisection** `π(x) = (Q(x) − M(x))/2 − C_3*(x)`,
   integer-exact for x ≤ 10⁶.
2. **The bridge identity** `NS_o(x) = (x − Q − L + M)/2`, integer-exact
   for x ≤ 10⁶, with one-line analytic proof
   `Σ μ²(n) λ(n) = Σ μ(n) = M(x)`.
3. **The 4-cell decomposition** `S_e + S_o + NS_e + NS_o = x` with the
   four explicit formulas in `definition.md`.
4. **Empirical Möbius-side parity independence**:
   `I(B(x) mod 2 ; C_3*(x) mod 2) ≈ 1.15 × 10⁻⁵` bits at N = 10⁶
   (parallel of E1.6).

## What this is NOT

Not an algorithmic improvement. The Möbius-bisection has the same
C-circular structure as E2.2: extracting `π(x) ~ x/log x` from
`(Q − M)/2 ~ 0.304·x` requires `C_3*(x) ~ 0.304·x` to high accuracy,
and computing `C_3*(x)` reduces to counting squarefree k-almost-primes
which itself reduces to `π` in arithmetic progressions. No polylog
opening.

Not a deep theorem. The Möbius-bisection is a one-page exercise that
any working analytic number theorist would derive on demand. Its
absence from `EDGES.md` and `CLOSED_PATHS.md` is therefore a project
catalogue gap rather than a genuine novelty.

## Self-grade — B

Reasoning: this is a **substantive refinement of E2.2** (per the
B-grade clause "(i) Refinement of an existing edge with a precise new
statement that extends its scope"). The refinement exhibits two
separately-verified objects:

* a parallel bisection on the squarefree side (`(Q − M)/2 − C_3*`),
* a closed-form bridge identity `NS_o = (x − Q − L + M)/2`.

It is honestly **not A-grade**: no published-paper-grade NT theorist
would have trouble deriving these, and there is no algorithmic content.

The natural verify-mode falsifier: "does the C₃-vs-C₃* statistical
independence claim hold beyond N = 10⁶?" — the script supports
larger N (sieve scales as `N log log N`), so this is rerunnable. The
4-cell decomposition itself is integer-exact by definition, so verify
can only attack the empirical pseudorandomness numbers.

A verify session may demote to C-grade DUPLICATE-PLUS on the grounds
that the Möbius bisection collapses to the standard
`Q − M = 2 S_o, S_o = π + C_3*` rearrangement — i.e., the "novel
content" is a labelling exercise. That demotion would be defensible.
The B-grade self-claim rests on the bridge identity being a precise
new closed-form statement that did not previously appear in the
project, plus the empirical Möbius-side independence (parallel to
E1.6) being a measurable fact, not a label.

## Successor proposals (paradigm-shift mode self-extension)

* **C13.a** Generalise both bisections to the **wheel-W-coprime**
  setting: define `L_W, M_W, Q_W, π_W, C_{3,W}, C_{3,W}*` analogously,
  derive `π_W(x) = (Q_W(x) − M_W(x))/2 − C_{3,W}*(x)` and the
  corresponding bridge identity. (Composes with E1.6's S219
  wheel-graded bisection lift.)
* **C13.b** Test whether the (`A mod q, B mod q`) joint distribution
  for prime `q ≥ 3` exhibits the same near-independence (S70 / C1
  extends E1.6 to prime q; the Möbius-side analogue is unmeasured).

Both are single-session B-target experiments. Neither requires an
external technique.
