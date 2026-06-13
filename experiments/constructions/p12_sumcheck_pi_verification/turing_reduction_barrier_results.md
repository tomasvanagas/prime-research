# turing_reduction_barrier.py — results

**Question (open item 3, the COMPUTATIONAL lower-bound face, direction (i)).**
S512 (`sharpP_probe.py`) settled the *parsimonious* (many-one) #P-hardness route:
a reduction `#A → π` with `π(x(w))=#A(w)` realizes the target count `c` by landing
`x ∈ [p_c, p_{c+1}−1]`, so the reduction map `c ↦ x(c)` IS the inverse prime
`p(·)` = the project goal — **circular**. S512 explicitly left open *direction (i):
a Turing (not parsimonious) reduction from a #P-complete count to a π-oracle that
sidesteps the `c↦x(c)` circularity.* This cycle closes direction (i)
**conditionally**, with a fine-grained-complexity argument the parsimonious
analysis cannot see: **π has its own sub-linear algorithm, and that algorithm is
itself the barrier.**

Code: `turing_reduction_barrier.py` (`--selftest`, 52 checks; `--pi-bench`,
`--blowup-nmax`). Self-contained number theory (sieve, n-th prime, a real
sub-linear `pi_lucy`, a brute `count_sat`); no chain imports.

## The barrier (conditional theorem)

Let `α = inf{ a : π(x) ∈ TIME(x^{a+o(1)}) }`.

> **THEOREM.** Assume SETH (Impagliazzo–Paturi). There is **no** poly(N)-time
> Turing reduction from #SAT to a π-**function**-oracle in which every oracle
> query has bit-length `≤ c·n` with `c < 1/α` (n = #SAT variables). Equivalently,
> any such reduction must query π at some `x ≥ 2^{(1/α)·n}`.

*Proof.* If `M^π` computes #SAT in poly(n) time with all queries `≤ c·n` bits and
`c < 1/α`, pick `a'` with `α < a' < 1/c`. Replace each π-query by the time-`x^{a'}`
algorithm (cost `≤ 2^{a'·c·n}` each); `M` makes poly(n) queries + poly(n) work, so
it computes `#SAT(w)` in `2^{a'c·n + o(n)}` with `a'c < 1`. Then SAT (decide
`#SAT>0`) is solved in `2^{(1−δ)n}`, `δ>0`, for every k — refuting SETH. ∎

So the **critical query-blowup is `c* = 1/α`.** The barrier covers *adaptive,
non-parsimonious* Turing reductions; the parsimonious case (S512) is the `c≈1`
slice of it.

## Headline

1. **π is sub-linear ⇒ a barrier `c* > 1` exists for the project's OWN algorithm.**
   `--pi-bench` measures the real `pi_lucy` (Lucy_Hedgehog, `O(x^{3/4})` time,
   `O(√x)` space) over `x = 10⁴…10⁸`: fitted time exponent **α ≈ 0.66–0.70**
   (sub-`3/4` in this range; the `t/x^{0.75}` column DECREASES, confirming the true
   exponent is below 0.75). Any `α < 1` gives `c* = 1/α > 1`.

   | π-algorithm | exponent α | barrier `c* = 1/α` | source |
   |---|---|---|---|
   | **measured** `pi_lucy` | **≈ 0.66–0.70** | **≈ 1.43–1.52** | this script (real timing) |
   | combinatorial LMO / Deléglise–Rivat | 2/3 | 1.5 | proven, cited |
   | analytic Lagarias–Odlyzko | 1/2 | **2.0** | best known, cited |
   | **PROJECT GOAL** (polylog) | 0 | ∞ | hypothetical |

2. **The natural (parsimonious) reduction has blowup `c → 1 < c*` ⇒ forbidden.**
   To encode "count = `C=2^n`" as `π(x)=C`, the smallest `x` is `p_C` (the C-th
   prime), bit-length `≈ n + log n`, so `c = bitlen(p_C)/n → 1` (PNT). Measured
   exactly via sieve: `c` falls `1.50 → 1.20` over `n=4…22`; PNT extrapolation
   `c = 1 + log₂(ln C)/n → 1.00`. Since `1 < c*` for every sub-linear π, **all
   parsimonious-blowup Turing reductions are SETH-forbidden** — independently of,
   and matching, the S512 circularity.

3. **Worked #SAT → π encoding.** A concrete 3-CNF (n=12, 14 clauses, `#SAT=480`):
   `π(x)=480` is realized by `x ∈ [p₄₈₀, p₄₈₁−1] = [3413, 3432]`, smallest
   `x=3413` (verified `π(3413)=480`, `π(3412)=479`); instance blowup `12/12 = 1.0`,
   worst-case-over-n=12-instances (count `~2¹²`) blowup `1.333`. Both `< c*=2`
   (analytic), so **FORBIDDEN under SETH** (would give `#SAT` in `2^{0.667n}` /
   `2^{0.889n}`). Also circular (S512): computing `x=p_C` from the formula is the
   inverse-prime goal.

4. **COROLLARY (the project tie-in): polylog-π XOR #P-hard-π, under P ≠ NP.**
   With `α→0` (the goal), `c* = 1/α = ∞`: a polylog π-oracle answers every query
   in poly(N) time, so a poly-time Turing reduction from #SAT to π would put #SAT
   in poly(n) time ⇒ **P = #P** (⇒ P = NP). Hence, assuming **P ≠ NP**:
   - if π has a polylog algorithm, π is **NOT** #P-hard (poly-time Turing);
   - a #P-hardness **proof** for π (poly-time Turing) proves **no polylog
     algorithm exists** — i.e. proves **this project impossible**.

   The two open questions ("is π #P-hard?" / "does π have a polylog algorithm?")
   are the **two horns of one dichotomy**; the √x algorithm sits exactly between
   them (`c* = 2`).

## What this is / is NOT

- The meta-principle ("a sub-linear target + a tight reduction = a sub-linear
  source") is **folklore fine-grained complexity**; the **contribution** is
  (1) applying it to π with π's own √x algorithm, (2) the quantitative threshold
  `c* = 1/α` and the MEASURED parsimonious blowup `c → 1 < c*`, (3) the
  polylog/#P-hardness mutual-exclusivity corollary tying the barrier to the goal.
- **Conditional** (SETH). SETH suffices because a `2^{(1−ε)n}` #SAT algorithm
  decides every k-SAT in that time, refuting SETH; no separate "#SETH" needed.
- **Scope.** Covers poly-time Turing reductions with query-blowup `< c*`.
  Super-`c*`-blowup reductions (querying π at `x ≥ 2^{(1/α)n}`, *exponentially*
  far from where the answer `p_C ≈ 2^n` lives) are **formally open** — but no
  natural construction realizes them, and there the analytic-π simulation already
  costs `≥ 2^n`, so they yield no #SAT speedup either.
- This is a **negative** result (a barrier to hardness), so the project's
  "circularity check" falsifier (does a claimed reduction secretly compute `p(·)`?)
  is not the relevant discipline here; the relevant honesty check is whether the
  conditional theorem's proof is sound (it is the standard simulation argument)
  and whether the blowup measurement is real (it is — exact via sieve + PNT).

## What would falsify each claim

- **(barrier)** A poly-time Turing reduction `#SAT → π` with all queries `< (1/α)n`
  bits, for the true π exponent α — under SETH this is impossible; exhibiting one
  would refute SETH (a major result), not the logic here.
- **(1)** `pi_lucy` time exponent `≥ 1` (not sub-linear) — not observed
  (measured α ≈ 0.66–0.70 over 10⁴–10⁸; selftest case 16 asserts α ∈ (0.4, 1.0)).
  The cited 2/3 and 1/2 rungs are established algorithms.
- **(2)** Parsimonious blowup not → 1, or `> c*` — not observed (exact c falls
  1.50→1.20 over n=4..22, PNT → 1.00; selftest cases 13–14).
- **(4)** A polylog π algorithm AND a poly-time Turing #P-hardness reduction
  coexisting without collapsing P=#P — impossible by the simulation (selftest
  case 11: `seth_forbidden(α=0, c)` True for every finite c).

## Selftest (52 checks)

`pi_lucy == pi_sieve` over a range incl. boundaries (0,1,2,3,…,222222) and the
known `π(10⁶)=78498`; `nth_prime` exact (1..10, 168→997, 1000→7919); realizing
prime `π(p_C)=C, π(p_C−1)=C−1` for C∈{1,10,100,500,1229}; `count_sat` on trivial
formulas + agreement with an independent enumeration on a random 3-CNF;
`fit_exponent` exact recovery on `x^{1/2}, x^{3/4}, x, 3x^{2/3}`; barrier
arithmetic `c*=1/α` and the SETH-verdict boundary (`c<c*` ⇔ forbidden, incl.
α=0 → every finite c forbidden, and the `c=2.0` exactly-not-`<2` boundary);
exact-vs-PNT parsimonious blowup agreement and monotone decrease to ~1; the
worked example (asserts `π(p_C)=C` internally); measured `pi_lucy` exponent
sub-linear and values == sieve.

## Bottom line for open item 3

Direction (i) — a Turing reduction sidestepping the S512 circularity — is now
**conditionally closed**: under SETH, no poly-time Turing reduction from #SAT to
π has query-blowup `< 1/α`, and the natural (parsimonious) blowup is `c → 1`, so
**every natural Turing reduction is forbidden.** The same √x algorithm that the
verification program (S491–S516) builds to *verify* π succinctly is what *blocks*
π from being #P-hard under natural reductions. And the corollary makes the
project's two big questions one dichotomy: **#P-hardness of π and a polylog
algorithm for π cannot both hold (P ≠ NP).** The remaining genuinely-open lever
for open item 3 is now narrowed to: (a) super-`c*`-blowup Turing reductions
(unnatural, no construction known), or (b) a **non-sieve sub-√x witness** for
`L_π` that also beats the S511 information floor — the other half of the NEXT
ACTION.

Cross-refs: `sharpP_probe.py`/S512 (the parsimonious half), `cert_incompressibility.py`/
S511 (the information floor), CLOSED_PATHS row 175 (corrected S512),
`novel/succinct_verification_of_pi.md` (the √x verification stack this barrier
mirrors). SETH: Impagliazzo–Paturi 2001. Analytic π: Lagarias–Odlyzko 1987.
Combinatorial π: Lagarias–Miller–Odlyzko 1985, Deléglise–Rivat 1996. Lucy_Hedgehog:
folklore (Lucy 2010).
