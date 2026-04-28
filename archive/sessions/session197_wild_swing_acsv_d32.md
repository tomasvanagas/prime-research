# Session 197 — Wild swing on §D.D32 (Pemantle-Wilson ACSV diagonal extraction of `χ_P`)

**Mode:** wild_swing (full-session single attempt; permission to fail).
**Target picked:** §D.D32 (ACSV smooth-point diagonal extraction).
**Channelled mathematician:** Pemantle / Flajolet (combinatorial
asymptotics) — but the actual lever is the meta-fact that ACSV requires
D-finite input, plus the rigorous combination of classical complex-
analysis / D-finite theory.
**Self-grade: B** (substantive structural refinement: an unconditional
theorem closing an entire cross-domain technique class with a fresh
edge E7.21 and pre-stated falsifiability conditions; ambitious A-grade
target was a polylog evaluator from ACSV; B-grade outcome is the clean
unconditional barrier).

---

## 0. Why I picked D32

The wild-swing prompt's default preference list (§C1, §A1, §B1, §A3,
§D4, §C2) is exhausted: §C1 closed S71, §A1 closed S84, §B1 closed
S92, §A3 closed E7.12 S79, §D4 closed E7.13 S80, §C2 closed S123. The
prompt says: pick the highest-A-grade-probability OPEN attack the
project hasn't already attempted.

Among the open targets, D32 had three strong properties:

1. The A-grade outcome is a **fundamental closed-form identity for
   `π(n)`** with direct polylog-evaluation consequence (Pemantle-Wilson
   smooth-point asymptotic gives `π(n) ≈ ζ^{*-n} · poly(n)` if `f` is
   the diagonal of a rational `F`). This is the most ambitious A-grade
   outcome of any §D entry.
2. The B-grade fallback was explicitly written as: "prove `Σ χ_P(n) z^n`
   is not the diagonal of any rational `F(z_1, ..., z_d)` for `d ≤ D_*`
   — a NEW negative-shape structural barrier on the algebraic class of
   generating functions admitting `χ_P`-diagonal encodings." This is
   provable using a clean assembly of Pólya-Carlson + Skolem-Mahler-Lech
   + finite-singularity D-finite theorem + Furstenberg-Lipshitz.
3. Cross-domain ACSV was UNUSED in CROSS_DOMAIN_TECHNIQUES §7, and the
   B-grade closure produces a structurally new edge.

The A-grade outcome would have falsified `f` is non-D-finite (a
classical folklore result, but not previously assembled in the project
as an unconditional structural barrier). The B-grade outcome — which is
what materialised — strengthens CLOSED_PATHS row 584 from empirical at
order ≤ 20 to **unconditional at all orders**, which is genuinely new
project content.

---

## 1. What I produced

### 1.1 Theorem (S197) — unconditional, all orders

> **Theorem.** Let `f(z) := Σ_{n≥1} χ_P(n) z^n`. Then:
>
>   - (a) `f` has `|z|=1` as natural boundary.
>   - (b) `f` is not D-finite (not P-recursive).
>   - (c) `f` is not algebraic over `ℚ(z)`.
>   - (d) for every `d ≥ 1`, `f` is not the diagonal of any rational
>     `F(z_1, ..., z_d) ∈ ℚ(z_1, ..., z_d)`.
>   - (e) for every `d ≥ 1`, `f` is not the diagonal of any algebraic
>     multivariate `F`.

### 1.2 Proof assembly

  1. **Cauchy-Hadamard:** ROC of `f` is exactly 1.
  2. **Pólya-Carlson 1916/1921:** integer-coefficient power series
     with ROC = 1 is rational or has `|z|=1` as natural boundary.
  3. **Non-rationality:** bounded-integer LRS is eventually periodic;
     for any candidate period `T` and any prime `p ≥ N_0`, periodicity
     forces `p, p+T, p+2T, ...` all prime, but `p + pT = p(1+T)` is
     composite. Contradiction. Combined with (2): `f` has `|z|=1` as
     natural boundary.
  4. **Stanley *EC2* Thm 6.4.6 / Flajolet-Sedgewick *AC* Thm B.13:**
     D-finite series have at most finitely many singular points.
     Combined with (3): `f` is not D-finite — proves (b).
  5. Algebraic ⊂ D-finite ⇒ (c).
  6. **Furstenberg 1967** (`d=2`) / **Lipshitz 1988** (`d ≥ 2`):
     diagonals of rational/algebraic multivariate are D-finite.
     Combined with (b): no rational/algebraic `F` exists ⇒ (d, e).
     ∎

### 1.3 Empirical companion — `experiments/constructions/acsv_pi_diagonal/`

Five tests in `acsv_pi_diagonal.py`:

  - **T1.** ROC = 1 confirmed at `N = 200000` (sanity check).
  - **T2.** P-recursion search at `d ≤ 30, r ≤ 8` over `χ_P` shows no
    holonomic signal (`χ_P / random` smallest singular value ratio
    `~ 1.1–1.6` in 9/30 cells, `~ 1.5–10` in remaining cells consistent
    with random's coincidental fitting at higher `r`). Extends
    CLOSED_PATHS row 584 (`d ≤ 20`) — but the rigorous theorem above
    obviates further empirical extension.
  - **T3.** Eventual-period search at `T ≤ 200`, `n ∈ [N/2, N]`,
    `N = 200000`. Best agreement = 0.88 at `T = 180` (random density-
    matched coincidence is ≈ 0.85). No period reproducing `χ_P`
    exactly — empirical confirmation of the rigorous step (3).
  - **T4.** **Erdős-Turán equidistribution at work.** Closest root of
    `f_N(t) := Σ_{n ≤ N} χ_P(n) t^n` to `t = 1`:

| `N`  | closest root distance to `t = 1` |
|------|-----------------------------------|
| 64   | 0.0927                            |
| 128  | 0.0476                            |
| 256  | 0.0234                            |
| 512  | 0.0115                            |
| 1024 | 0.0059                            |

    Halves as `N` doubles, consistent with Erdős-Turán 1950 *Ann. Math.*
    51 (integer-coefficient polynomials with bounded `L^∞` have roots
    equidistributing on `|t|=1`). **A stable Pemantle-Wilson critical
    point does not exist in the limit.**
  - **T5.** Algebraic-relation search at degree `≤ 6` on the diagonal
    sequence `π(n)` for `n ≤ 128`: zero kernel dimension at `D ∈
    {2, 3, 4, 5, 6}` — no non-trivial algebraic relation up to degree
    6. (Sanity check; the rigorous proof gives the all-degree result
    via Lemma D.)

### 1.4 New edge E7.21 (added to EDGES.md)

> `f(z) := Σ χ_P(n) z^n` has `|z|=1` as natural boundary (Pólya-
> Carlson); hence is not D-finite, not algebraic, and not the
> diagonal of any rational/algebraic multivariate `F` for any
> `d ≥ 1`. EVS shape.

This is the project's **first complex-analytic structural barrier**
on `χ_P` — distinct from the information-theoretic edges (E1.5),
spectral edges (E7.1–7.20), and archimedean unit-circle quantitative
measurements (E2.20 Mahler / E2.21 Newman / E2.27 wavelet Hölder).

### 1.5 What this rules out, in one citation

  - **Pemantle-Wilson ACSV** smooth-point asymptotics for `χ_P` or
    `π(n)`.
  - **Wilf-Zeilberger creative telescoping** for `χ_P`/`π(n)` sums.
  - **Holonomic gradient methods** (Hibi-Nishiyama-Takayama 2013).
  - **Karr-Schneider Π Σ-symbolic summation** for `χ_P`/`π(n)`.

Each of these requires a D-finite input. By Theorem 1(b), `f` is not
D-finite — so each is closed by E7.21 in one line.

---

## 2. Self-evaluation per CLAUDE.md

> **Q1: What did I produce that was not in the project before this
> session?**

(a) An unconditional theorem: `Σ χ_P(n) z^n` has `|z|=1` as natural
boundary, is not D-finite, is not algebraic, and is not the diagonal
of any rational/algebraic multivariate generating function for any
`d ≥ 1`. This **strictly strengthens CLOSED_PATHS row 584** from
empirical at order ≤ 20 to unconditional all-order.

(b) A new EDGES.md entry **E7.21** — the project's first complex-
analytic structural barrier (as opposed to spectral, information-
theoretic, or archimedean-quantitative).

(c) A new closure of ATTACK_VECTORS §D.D32 — the only known frontier
attack vector that requires diagonal-extraction machinery, now closed
unconditionally.

(d) Cross-domain ACSV §7 promoted PROPOSED → USED I.

(e) Successor proposals D32.a (transcendental-`F` PW) and D32.b
(Schur-class non-rational expansions, NEW cross-domain technique
proposed for §1 of CROSS_DOMAIN_TECHNIQUES).

(f) Empirical companion confirming the theoretical claims at finite
`N`, including the Erdős-Turán equidistribution observation that
`f_N(t)` roots accumulate at `t = 1` as `N → ∞`.

> **Q2: What edges did my work compose or cite?**

Composes: `χ_P` is not eventually periodic (folklore in CLOSED 23) +
Pólya-Carlson 1916/1921 + Stanley/FS finite-singularity of D-finite +
Furstenberg 1967 / Lipshitz 1988. Cites edges E2.20 (Mahler measure on
archimedean unit circle), E2.21 (Newman L^∞), E2.27 (wavelet Hölder
on `D(x)`) as related but structurally distinct.

> **Q3: If my session produced only duplicate closures, why?**

Not applicable — the session produced a new theorem that does not exist
in the project's prior closures. The closest prior result is
CLOSED_PATHS row 584 (empirical D-finite ruling-out at low order); my
theorem strictly strengthens this.

> **Q4: What is the next-action for the next agent?**

Two options recorded in `experiments/constructions/acsv_pi_diagonal/
acsv_pi_diagonal_results.md` §3.4 and in ATTACK_VECTORS.md "Closed
attacks":

  - **D32.a (1 session, same cross-domain technique):** test whether a
    transcendental-but-holonomic-D-module `F(z_1, ..., z_d)` (involving
    `e^{z_i}`, `log z_i`, hypergeometric components) escapes Lipshitz's
    diagonal closure. Pemantle-Wilson saddle-point machinery has
    partial extensions to such `F` but no proven asymptotic theorem.
  - **D32.b (2 sessions, NEW cross-domain technique):** Schur-class
    non-rational generating-function expansions. Schur-class functions
    `S(z) = Σ a_n z^n` with `|S(z)| ≤ 1` on `|z| < 1` admit Adamyan-
    Arov-Krein 1971 / Schur-Pick interpolation machinery, orthogonal
    to ACSV. Does `f / (1 + |f|)` (re-normalised to Schur class) admit
    structurally extractable asymptotics? **NEW cross-domain
    technique** for the project (not currently in
    CROSS_DOMAIN_TECHNIQUES); add under §1.

Neither successor is a wild swing — both are B-grade single- to
two-session probes within new technique classes.

---

## 3. Why this is B-grade and not A or C

**Not A:** A-grade required the *positive* outcome — actually
identifying a rational/algebraic `F` and a smooth critical point, or
achieving partial-positive content from the cross-domain attack. The
Theorem rules this out structurally; my outcome is the negative-shape
B-grade case explicitly anticipated by the D32 attack-vector text.

**Not C:** C-grade is "duplicate closure of CLOSED_PATHS with
non-trivial structural reason." My theorem is strictly stronger than
CLOSED row 584 (unconditional vs empirical at low order), and it
closes a fresh ATTACK_VECTORS frontier that was OPEN at session start.
The structural reasoning is the rigorous combination of four classical
results (Pólya-Carlson, Skolem-Mahler-Lech, finite-singularity D-finite
theorem, Furstenberg-Lipshitz) — not previously assembled in the
project as an explicit ACSV-class barrier.

**B is honest:** the cross-domain ingredient is the meta-fact that
ACSV-class techniques require D-finite input, plus the four-citation
proof. The theorem itself is novel as a project artefact but the
constituent classical results are well-known in their respective
fields. The actual *technique* importation is from the D-finite-theory
literature (Stanley *EC2* §6.4, Flajolet-Sedgewick *AC* App. B), which
has not been used in this project before — but the importation is
classical-result-citation rather than new mathematics. B is correct.

The *honest-failure clause* of CLAUDE.md applies in reverse: this is
not a duplicate closure, it is a substantive structural strengthening
of a previously empirical result. B-grade.

---

## 4. CROSS_DOMAIN_TECHNIQUES update

§7 row "Analytic Combinatorics in Several Variables" promoted
PROPOSED (D32) → USED I (S197, edge E7.21). The cross-domain content
includes the D-finite-theory citations (Stanley *EC2*, Flajolet-
Sedgewick *AC* Appendix B) plus the classical Pólya-Carlson 1916/1921,
Furstenberg 1967, and Lipshitz 1988.

---

## 5. Note on the run.sh state

This session was run-197 in the harness rotation. The .commit_state at
session start had the Connes amortisation arc at slot 5/5 with
next_action "WRAP synthesis". The user's prompt explicitly invoked the
**wild_swing** mode override, so the Connes thread synthesis was NOT
written this session — that remains the next-action for any subsequent
commit-mode slot. The .commit_state should be advanced by the harness
or the next commit slot.

After this session, frontier-saturation pattern stands:
ATTACK_VECTORS.md still has many OPEN entries (§A.A2, §A.A4, §A.A6,
§A.A7 plethysm sub-frame, §B.B3, §B.B4, §B.B5, §C.C3, §C.C6,
§D.D1, §D.D3, §D.D6, §D.D11, §D.D12, §D.D14, §D.D15, §D.D16, §D.D17,
§D.D18, §D.D19, §D.D21, §D.D23, §D.D24, §D.D25, §D.D26, §D.D28,
§D.D33, §D.D34, §D.D35, §D.D36, §D.D37, §D.D38, §D.D39, §D.D40,
§D.D41, §D.D43.b/.c, §D.D45.a-c, §D.D46, §D.D47, §D.D48, plus 2
successor proposals D32.a/.b from this session). No frontier_gen
trigger condition met (open count well above threshold; A-grade
scarcity above threshold but pre-existing — flagged in S192).

---

## 6. Files touched

- **NEW:** `experiments/constructions/acsv_pi_diagonal/acsv_pi_diagonal.py`
- **NEW:** `experiments/constructions/acsv_pi_diagonal/acsv_pi_diagonal_results.md`
- **NEW:** `archive/sessions/session197_wild_swing_acsv_d32.md` (this file)
- `EDGES.md` — added E7.21 between E7.20 and E7.14
- `ATTACK_VECTORS.md` — D32 marked CLOSED inline, full closure entry
  appended to "Closed attacks" section
- `CLOSED_PATHS.md` — appended new D32 closure row at line 819+
- `CROSS_DOMAIN_TECHNIQUES.md` — §7 ACSV row promoted PROPOSED →
  USED I (S197, edge E7.21)
