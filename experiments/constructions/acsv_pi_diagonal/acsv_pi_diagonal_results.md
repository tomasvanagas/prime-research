# D32 / ACSV — Pemantle-Wilson smooth-point diagonal extraction of `χ_P`

**Session:** S197 (wild-swing mode)
**Verdict:** D32 **CLOSED**, mode **I** (negative-shape structural barrier),
**B-grade**.
**Cross-domain technique:** ACSV (Pemantle-Wilson 2013/2024) — promote
PROPOSED → USED (I).
**Edges cited:** E2.20 (Mahler), E2.21 (Newman L^∞), E7.21 (new — see
EDGES.md).
**Falsification:** exhibit a multivariate rational `F(z_1, ..., z_d) ∈
ℚ(z_1, ..., z_d)` and an integer `d ≥ 1` such that
`[z_1^n ⋯ z_d^n] F = π(n)` (or `χ_P(n)`) for all sufficiently large `n`.
The theorem below shows this is impossible.

---

## 0. The wild-swing target

D32 asks: does there exist a rational `F(z_1, ..., z_d)` whose **diagonal
coefficient extraction** equals `π(n)` (or `χ_P(n)`)? If yes, the
Pemantle-Wilson smooth-point method gives a polylog evaluator for `π(n)`
via the asymptotic `[z_1^n ⋯ z_d^n] F = ζ^{*-n} (2πn)^{-(d-1)/2} /
√(H(ζ*))`. **A-grade outcome:** find such `F`. **B-grade outcome:** rule
it out structurally.

This session takes the B-grade outcome — but with an **unconditional
proof** that strengthens CLOSED_PATHS line 584 (which was empirical at
order ≤ 20, polynomial degree ≤ 8) to all orders simultaneously.

---

## 1. Theorem (S197)

> **Theorem 1.** Let `f(z) := Σ_{n≥1} χ_P(n) z^n` be the prime-indicator
> generating function. Then:
>
>   (a) `f` has the unit circle `|z|=1` as a **natural boundary**;
>
>   (b) `f` is **not D-finite** (not P-recursive) — does not satisfy any
>       linear ODE `Σ_{i=0}^k P_i(z) f^{(i)}(z) = 0` with polynomials
>       `P_i ∈ ℂ[z]`;
>
>   (c) `f` is **not algebraic** over `ℚ(z)` (or `ℂ(z)`);
>
>   (d) For every `d ≥ 1`, `f` is **not the diagonal** of any rational
>       `F(z_1, ..., z_d) ∈ ℚ(z_1, ..., z_d)`;
>
>   (e) For every `d ≥ 1`, `f` is **not the diagonal** of any algebraic
>       multivariate `F`.

**Corollary (D32 closure).** The Pemantle-Wilson smooth-point
machinery — applied to any rational or algebraic multivariate `F` —
cannot produce an asymptotic for `χ_P(n)` or `π(n)`.

---

## 2. Proof of Theorem 1

We assemble four classical results.

**Lemma A (radius of convergence).** `f` has radius of convergence
exactly `1`.

*Proof.* `χ_P(n) ∈ {0, 1}`, so `lim sup |χ_P(n)|^{1/n} = 1` (the
sequence equals 1 along the prime subsequence, which is infinite by
Euclid; equals 0 elsewhere). By Cauchy-Hadamard, the radius of
convergence is `1/1 = 1`. ∎ (Verified empirically by T1 in
`acsv_pi_diagonal.py`.)

**Lemma B (Pólya-Carlson 1916/1921).** Let `g(z) = Σ a_n z^n` be a
power series with integer coefficients `a_n ∈ ℤ` and radius of
convergence exactly `1`. Then either `g` is a rational function whose
poles lie on `|z|=1` (and the only pole on `|z|≤1` is at `z=1` if `a_n`
are bounded), **or** `|z|=1` is a natural boundary for `g`.

*Reference.* Pólya 1916 *Sitzungsber. Berlin Math. Ges.*; Carlson 1921
*Math. Z.* 9, 1; standard textbook treatment in Bieberbach,
*Analytische Fortsetzung* (Springer 1955), Ch. III.

**Lemma C (`f` is not rational).** `f(z) = Σ_{n≥1} χ_P(n) z^n` is not
a rational function.

*Proof.* Suppose `f(z) = P(z)/Q(z)` with `P, Q ∈ ℂ[z]`, `Q(0) ≠ 0`.
Then `(χ_P(n))_{n≥0}` satisfies a linear recurrence with constant
coefficients of order `deg Q`. Since `χ_P` is bounded and
integer-valued, the recurrence forces `(χ_P(n))_{n ≥ N_0}` to be
**eventually periodic** with some period `T ≥ 1` (a bounded integer
LRS depends only on `deg Q` previous values, giving finitely many
possible state vectors and hence eventual periodicity by the pigeonhole
principle).

Pick any prime `p ≥ max(N_0, 2)` (possible by Euclid). Periodicity
forces `χ_P(p + kT) = χ_P(p) = 1` for every `k ≥ 0`, i.e., `p + kT` is
prime for every `k ≥ 0`. Take `k = p`: then `p + pT = p(1 + T)` with
`1 + T ≥ 2`, so `p(1 + T)` has both factors `p ≥ 2` and `1 + T ≥ 2`,
hence is composite. Contradiction. ∎

(Note: the Skolem-Mahler-Lech theorem (Skolem 1934, Mahler 1935, Lech
1953) gives a stronger structural assertion — the zero set of any
non-degenerate linear recurrence sequence is a finite union of APs
plus finite — but the bounded-integer special case used here is
elementary via pigeonhole. Verified empirically by T3 in
`acsv_pi_diagonal.py`: best eventual-period agreement at `T ≤ 200` is
≈ 0.88, no period reproducing `χ_P` exactly.)

**Combining A + B + C: `f` has natural boundary on `|z| = 1`.** ∎

**Lemma D (D-finite ⇒ finitely many singularities).** If `f` is
D-finite, satisfying `Σ_{i=0}^k P_i(z) f^{(i)} = 0` with `P_i ∈ ℂ[z]`
and `P_k ≢ 0`, then `f` extends to a holomorphic function on
`ℂ \ Z(P_k)`, where `Z(P_k)` is the (finite) zero set of `P_k`.

*Reference.* Standard linear-ODE theory in the complex domain: see
Ince, *Ordinary Differential Equations* (1956), Ch. XV; Stanley,
*Enumerative Combinatorics II* (CUP 1999), Theorem 6.4.6; Flajolet-
Sedgewick, *Analytic Combinatorics* (CUP 2009), Theorem B.13.

**Conclusion (b).** `f` having `|z|=1` as natural boundary contradicts
Lemma D's finite-singularity assertion. Hence `f` is **not D-finite**. ∎

**Conclusion (c).** Algebraic ⊂ D-finite (Comtet 1964 / standard). So
`f` is not algebraic. ∎

**Lemma E (Furstenberg 1967 / Christol 1990 / Lipshitz 1988).** Let
`F(z_1, ..., z_d) ∈ ℚ[[z_1, ..., z_d]]` be a multivariate power series
with `F` rational (resp. algebraic, resp. D-finite as a multivariate
series). Then the diagonal `Δ F(z) := Σ_n [z_1^n ⋯ z_d^n] F · z^n` is:

  - **algebraic** if `d = 2` and `F` is rational (Furstenberg 1967
    *J. Algebra* 7, 271 = "Algebraic functions over finite fields");
  - **D-finite** for every `d ≥ 2` and `F` rational (Lipshitz 1988
    "The diagonal of a D-finite power series is D-finite", *J. Algebra*
    113, 373);
  - **D-finite** for every `d ≥ 2` and `F` D-finite (Lipshitz 1988,
    same).

**Conclusion (d).** Suppose `f = Δ F` for some rational `F` in `d`
variables. By Lemma E, `f` is D-finite. Contradicts (b). ∎

**Conclusion (e).** Algebraic ⊂ D-finite, so the same argument
applies. ∎

---

## 3. Why this closes D32 / ACSV — and what would not

### 3.1 The smooth-point Pemantle-Wilson machinery

Pemantle-Wilson Theorem 9.5.7 / 10.3.4 (PW 2013/2024) states: if
`F = G/H ∈ ℚ(z_1, ..., z_d)` and `ζ* ∈ ℂ^d` is a **smooth critical
point** of the singular variety `V_F := {H = 0}` minimising
`|ζ_1 ⋯ ζ_d|`, then for fixed direction `r = (1, ..., 1)`:

```
[z_1^n ⋯ z_d^n] F = (ζ*)^{-n·r} · (2πn)^{-(d-1)/2} · C · (1 + O(1/n))
```

where `C = G(ζ*) / √(det H(ζ*))` is computable from `F`. This
asymptotic gives a **polylog-time evaluator** for the diagonal
coefficient sequence — *provided* `F` is rational (or algebraic) and
`ζ*` is a smooth critical point.

### 3.2 What Theorem 1 rules out

Theorem 1(d, e): no `F ∈ ℚ(z_1, ..., z_d)` (rational or algebraic)
gives `f(z)` as a diagonal. So the smooth-point asymptotic above
**cannot be assembled** for `χ_P` or `π(n)`.

### 3.3 Empirical companion (T4)

Even if we IGNORE Theorem 1 and naively try the simplest D32
candidate

```
F_N(x, y) := f_N(xy) / (1 - xy),    where f_N(t) := Σ_{n ≤ N} χ_P(n) t^n,
```

so that `[(xy)^n] F_N = π(min(n, N))` for `n ≤ N`, the smooth-point
machinery still fails: the singular variety
`V_{F_N} = { f_N(xy) = 0 } ∪ { xy = 1 }`
has roots that **do not stabilise** as `N → ∞`. Empirically:

| `N`  | closest root of `f_N(t)` to `t = 1` |
|------|--------------------------------------|
| 64   | 0.0927                               |
| 128  | 0.0476                               |
| 256  | 0.0234                               |
| 512  | 0.0115                               |
| 1024 | 0.0059                               |

The closest root to `t = 1` halves as `N` doubles, consistent with the
**Erdős-Turán equidistribution theorem** (1950 *Ann. Math.*; integer-
coefficient polynomials with `O(1)` `L^∞` norm have roots that
equidistribute on `|t| = 1`). A polylog smooth-point asymptotic
`[(xy)^n] F = ζ^{*-n} · poly(n)` requires a *fixed* `ζ*` — but here
the relevant root moves with `N`, so no `N → ∞` limit asymptotic
exists.

### 3.4 What would NOT close D32

Theorem 1 explicitly leaves two escape doors that ACSV-folklore knows:

  - **Transcendental `F`.** If `F` is allowed to be transcendental in
    one or more variables (e.g., `F(z, e^z)`), Pemantle-Wilson's
    machinery does not directly apply — but neither do its asymptotic
    guarantees. The "ACSV" technique as published is for rational /
    algebraic `F`. (This escape is theoretically real but empty for
    polylog purposes: a transcendental `F` does not produce a
    closed-form polylog evaluator without further machinery.)
  - **Differential-algebraic `F`.** If we allow `F` to be the solution
    of a non-linear differential-algebraic system, Lipshitz's theorem
    does not apply. Same caveat: no published polylog asymptotic
    machinery exists for such `F`.

Neither escape door has produced any computable polylog asymptotic
for any analytic-NT object, anywhere in the literature. The closure is
therefore **structural** within the ACSV domain as published.

---

## 4. Connection to existing closures and edges

### 4.1 Strict strengthening of CLOSED_PATHS line 584

CLOSED_PATHS line 584 (Session 23): "Holonomic (D-finite) recurrence
for `π(n)`. FAIL, mode I. NOT holonomic for any order `d ≤ 20` with
polynomial degree `r ≤ 8`. Test/random ratio ~ 1.0–1.7." This was
**empirical** at finite order. Theorem 1 makes the result
**unconditional at all orders simultaneously**.

(T2 in `acsv_pi_diagonal.py` re-runs the empirical search out to
`d ≤ 30, r ≤ 8` and confirms no holonomic signal — but the
unconditional proof obviates further empirical extension.)

### 4.2 New negative-shape edge E7.21

The non-D-finiteness of `f(z)` is a structural fact that constrains
**any technique that requires a holonomic / algebraic / rational
generating function for `χ_P` or `π(n)`**. This includes, beyond ACSV:

  - **Telescoping / creative-telescoping algorithms** (Zeilberger, Wilf-
    Zeilberger): require D-finite input.
  - **Holonomic gradient methods** (Hibi-Nishiyama-Takayama 2013): same.
  - **Symbolic summation à la Karr-Schneider** in `Π Σ`-fields:
    requires `χ_P` to be in a `Π Σ`-extension.

Each is now closed by E7.21 with a single citation.

### 4.3 Distinct from existing E2.20 / E2.21 (archimedean unit-circle
edges)

E2.20 (Mahler measure of `f_N`) and E2.21 (Newman L^∞ of `f_N`) live
on the **archimedean** unit circle and measure **specific extremal
values** of `f_N` on `|z|=1`. E7.21 (this work) is a **structural
statement about the limit** `f = lim_{N → ∞} f_N`: the natural
boundary fact does not follow from E2.20 / E2.21 individually. It
*does* follow from the prior project rigorous work that `χ_P` is not
linearly recurrent (CLOSED line 23) plus Pólya-Carlson — but the
project had not previously assembled this composition into an explicit
ACSV-class barrier.

### 4.4 Connection to E1.5 (information-theoretic)

E1.5 says `H(π(x) mod m | π(x-1) mod m) = h_2(π(X)/X) + O(1/π(X))`,
independent of `m`. This is an **information-theoretic** measure of
randomness on conditional steps. The new edge E7.21 (non-D-finiteness)
is an **algebraic-analytic** measure of structurelessness on the full
generating function. Both confirm the same intuitive fact (`χ_P` is
"structureless" in distinct mathematical senses); E7.21 is the first
project edge that uses a **complex-analytic** structural argument
rather than a statistical or information-theoretic one.

---

## 5. Falsifiability statement

The closure of D32 by Theorem 1 is falsified by exhibiting any one of
the following:

**F1 (most direct).** A polynomial-coefficient linear recurrence
`Σ_{i=0}^k P_i(n) χ_P(n + i) = 0` with `P_i ∈ ℂ[n]`, `P_k ≢ 0`, valid
for all sufficiently large `n`. (Falsifies (b) — would falsify Theorem
1 outright since it would prove `f` is D-finite.)

**F2 (more lenient).** A polynomial `Q(z, w) ∈ ℂ[z, w]`, `Q ≢ 0`, with
`Q(z, f(z)) = 0` as formal power series. (Falsifies (c).)

**F3 (target form).** A multivariate rational `F(z_1, ..., z_d) ∈
ℚ(z_1, ..., z_d)` and a vector `r ∈ ℤ_{>0}^d` such that
`[z_1^{r_1 n} ⋯ z_d^{r_d n}] F = χ_P(n)` for all `n` sufficiently
large (the "diagonal" can be in any direction `r`, not only `r =
(1, ..., 1)`). (Falsifies (d).)

Either of F1 or F2 contradicts the rigorous combination Pólya-Carlson
+ Skolem-Mahler-Lech + finite-singularity D-finite theorem, so each is
*formally impossible*. F3 contradicts Lipshitz's theorem (Lemma E).
Anyone who finds such an `F` resolves the polylog `π(x)` problem.

---

## 6. Cross-domain ingredient (CROSS_DOMAIN_TECHNIQUES update)

ACSV §7 row promoted: PROPOSED → **USED (I)**. The application is the
**negative-shape barrier** of Theorem 1.

**Channelled mathematician:** Pemantle / Flajolet (combinatorial
asymptotics) — but the actual technique is the *meta*-fact that
ACSV-class techniques require a D-finite / algebraic input. The
classical results assembled (Pólya-Carlson, Skolem-Mahler-Lech,
Furstenberg, Lipshitz) are the genuine cross-domain content.

---

## 7. Summary — what to cite from this construction

If a future agent considers any of the following, cite E7.21:

  - "Maybe `f(z)` admits a holonomic recurrence at some order `d > 30`."
    → No. Theorem 1(b) is unconditional.
  - "Maybe a multivariate rational `F(z_1, ..., z_d)` encodes `π(n)`
    as a diagonal." → No, for any `d`. Theorem 1(d) closes this.
  - "Apply ACSV smooth-point asymptotics to `χ_P`." → Inapplicable.
    Theorem 1(d, e) closes the technique class.
  - "Apply Wilf-Zeilberger creative telescoping / holonomic gradient /
    Karr-Schneider symbolic summation to π(n)." → Each requires
    D-finite input. Theorem 1(b) closes them all.

The cleanest citation: **E7.21 (Pólya-Carlson natural-boundary
barrier; `f` is non-D-finite, hence ACSV-class techniques are
inapplicable to `χ_P`)**.

---

## 8. References

- Pólya 1916 *Sitzungsber. Berlin Math. Ges.* — power series with
  integer coefficients.
- Carlson 1921 "Über Potenzreihen mit ganzzahligen Koeffizienten"
  *Math. Z.* 9, 1.
- Skolem 1934 *Comptes Rendus du 8e Congrès des Mathématiciens
  Scandinaves*; Mahler 1935 *Proc. Akad. Wet. Amsterdam* 38; Lech 1953
  *Ark. Mat.* 2 — the Skolem-Mahler-Lech theorem.
- Furstenberg 1967 "Algebraic functions over finite fields" *J. Algebra*
  7, 271.
- Lipshitz 1988 "The diagonal of a D-finite power series is D-finite"
  *J. Algebra* 113, 373.
- Stanley 1999 *Enumerative Combinatorics II* CUP, Theorem 6.4.6.
- Flajolet-Sedgewick 2009 *Analytic Combinatorics* CUP, Theorem B.13.
- Bieberbach 1955 *Analytische Fortsetzung* Springer, Ch. III.
- Pemantle-Wilson 2013/2024 *Analytic Combinatorics in Several
  Variables* CUP, 2nd ed. — Theorem 9.5.7 / 10.3.4 (smooth-point).
- Erdős-Turán 1950 *Ann. Math.* 51 — equidistribution of polynomial
  roots.

---

## 9. Outcome grade

**Mode:** I (negative-shape structural barrier).
**Grade:** B (substantive refinement: unconditional theorem strengthens
existing empirical CLOSED line 584; closes a fresh ATTACK_VECTORS entry
with a new edge E7.21; ambitious-failure clause: the wild-swing target
A-grade was a polylog evaluator from ACSV; the realised outcome is a
clean unconditional barrier ruling out the entire ACSV technique class
for `χ_P`).
