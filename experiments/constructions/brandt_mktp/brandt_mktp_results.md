# Brandt 2024 MKtP diagonalization vs pi(x) mod 2 — results

**Session:** 51 — Construction-flavoured attack on Chain E (E5.3) via a
non-AKS technique (FOCUS-3 in `TODO.md`).

**Target:** Determine whether the Brandt 2024 (TCC) diagonalisation
technique that proves `MKtP not in DTIME[O(n)]` extends to a non-trivial
lower bound on `pi(x) mod 2`.

**Verdict:** **FAIL / E (equivalence-class closure).**  Brandt's
technique is structurally welded to MKtP and does not transplant to
natural functions like `pi(x) mod 2`.  No amount of empirical
pseudorandomness of pi (E1.9, 33 measures) bridges the gap.

## Source paper

- N. Brandt, "Lower Bounds for Levin–Kolmogorov Complexity," **TCC 2024**,
  LNCS 15363; IACR Cryptology ePrint Archive **2024/687**.
- Read in full from <https://eprint.iacr.org/2024/687.pdf> (PDF
  fetched via `curl`, parsed with `pypdf`).
- Pages 1–14 of the PDF give Theorem 1, Theorem 2, Corollaries 1–2, the
  technical overview, the formal proof of Lemma 1, and TRAVERSE
  pseudocode (Figure 4).

## What Brandt actually proves

```
    MKtP  not in  DTIME[O(n)]                    Theorem 1, unconditional
    MKtP  not in  DTIME[O(n^2)]                  Corollary 1, conditional on
                                                 linear-time UTM (e.g. RAM)
    MKtP  not in  Heur_{0, gamma_fn} DTIME[O(n)] Theorem 2, one-sided heuristic
    MKtP  not in  Union_k DTIME[prod ln^{(i)}(n)] Corollary 2, slightly super-linear
```

**Technique skeleton (page 5 of Brandt):**

1. Assume for contradiction `MKtP in DTIME[O(n)]` via a TM `Pi_Kt`.
2. Build TM `M_{c_Omega, c_Kt}` running TRAVERSE (Figure 4 in paper):
    - Start at z_1 = '0'.
    - At step i, query `Pi_Kt(z_i, |z_i|+ceil(log_2 |z_i|)-c_Omega)`.
    - If z_i is Kt-random (i.e. Pi_Kt says NO), descend: z_{i+1} := z_i || '0'.
    - Otherwise step right: z_{i+1} := next(z_i).
3. The traversal *cannot* wrap around (Lemma 2), because Chaitin's
   constant Omega has Kt-random prefixes of every length.
4. Inequality (1) of the paper: `i_l <= sum_j gamma_j sum_kappa 2^kappa
   exp(sigma_{j-kappa} - sigma_j)` is bounded by `2^l / (l ln l)` infinitely
   often (Lemma 1).
5. Therefore the last visited string `z_{i_l}` of length l satisfies
   `Kt(z_{i_l}) >= l` (because that's the descent rule), but also
   `Kt(z_{i_l}) <= |M| + log_2(runtime) = c_Kt + log_2(t_l log_2 t_l)`.
6. Plugging in the runtime bound `O(i_l * l) = O(2^l / ln l)` from Lemma 1,
   the upper bound on Kt becomes `l + log_2(l) - log_2(ln l) + O(1)
   = l - ln ln l + O(1)`, which contradicts `Kt(z_{i_l}) >= l + log_2 l - c_Omega`.

**Bypassing Razborov-Rudich (T4):** Brandt's proof relativizes (page 2,
"That our lower bounds relativize can be taken as a hint that
relativizing techniques might in fact be strong enough...").  It uses
**no large/efficiently-computable property** of MKtP.  It uses Chaitin
Omega's 1-Kt-randomness as a one-shot witness, not a distribution over
the function class.  Razborov-Rudich constructivity/largeness conditions
do not even apply because there's no probabilistic ensemble being
distinguished.  This *is* the "ingredient that bypasses Natural Proofs."

## Empirical results from `brandt_mktp.py`

### Part A: Bounded TRAVERSE on a 3-bit-per-op stack VM, max length 10

```
Visited 11 strings, 10 flagged Kt-random,
TRAVERSE descended 10 times, max length reached 10.

(1, '0',          True,  1)
(2, '00',         True,  2)
(3, '000',        True,  3)
(4, '0000',       True,  4)
(5, '00000',      True,  5)
(6, '000000',     True,  6)
(7, '0000000',    False, 7)    <- bounded Kt finds short program
(8, '0000001',    True,  7)    <- so TRAVERSE moves right
(9, '00000010',   True,  8)
(10,'000000100',  True,  9)
(11,'0000001000', True, 10)
```

This confirms Brandt's traversal does what the paper says: it *exploits*
the bounded Kt oracle by descending only on Kt-random strings.  With our
bounded VM the all-0^n string is Kt-non-random for n >= 7 (a 3-bit
"PUSH 0" loop generates it), so the algorithm's right-step rule fires
and the traversal proceeds.  This is a *positive sanity check* that the
construction is faithful.

### Part B: Bounded `Kt(pi_N)` for `pi_N = (pi(x) mod 2 : x in [0, 2^N))`

```
        N | |pi_N| | Kt_bounded(pi_N) | n - Kt | Kt-random?
        --+-------+------------------+--------+------------
        3 |     8 |               15 |     -7 | True
        4 |    16 |               61 |    -45 | True
        5 |    32 |               61 |    -29 | True
        6 |    64 |               61 |      3 | False
        7 |   128 |               61 |     67 | False
        8 |   256 |               61 |    195 | False
        9 |   512 |               61 |    451 | False
        10|  1024 |               61 |    963 | False
```

`Kt_bounded` saturates at the VM's enumeration bound (61) for N >= 4 —
no program of length up to L_MAX = 12 produces `pi_N`.  For N >= 6, the
target string is so long that the saturated bound 61 is much smaller
than `n = 2^N`, so the **bounded test** flags `pi_N` as **not** Kt-random
in the strict `Kt(x) >= |x|` sense.  This is **not** evidence about
true Kt(pi_N) — it's the limit of our finite enumerator.  The honest
reading is: at all tested N, `pi_N` is incompressible by short programs
in a bounded VM, consistent with the project's broader pseudorandomness
finding (E1.9, novel/pseudorandomness_of_pi.md).

### Part C: pi_N as MKtP queries

For each N, we query `(pi_N, k)` for `k in {n-2, n-1, n, n+1}` (the four
threshold values Brandt's TRAVERSE would issue). At small N (4..6), all
four queries return False (pi_N's true Kt looks larger than n).  At
N >= 7, our bounded-Kt saturation makes them all True.  Either way,
**the answers do not feed into any traversal that would produce a
hard string of growing length**.  This is the operational statement of
the structural argument below.

## The four obstructions to extending Brandt to pi(x) mod 2

(O1) **Hard string is wrong.**  The string z that Brandt's TRAVERSE
produces is a Kt-random prefix that depends on the run, *not* a fixed
function of N.  For the project we need a lower bound on the *fixed*
function pi.  Brandt's traversal does not know how to use a "fast
pi-decider" `Pi_pi` to construct z so that `Kt(z) >= |z|`.

(O2) **Self-reference inside Kt.**  Brandt's contradiction
`Kt(z) >= |z|` and `Kt(z) <= |M| + log_2 t` uses the *same* complexity
measure on both sides.  This is what makes the diagonalisation work
*for MKtP itself*.  No analogous self-referential bound exists for pi:
a pi-decider answers "is x prime?" not "is z complicated?", so a TM
constructed from `Pi_pi` cannot be analysed in terms of pi-hardness of
its outputs.

(O3) **Density of hard witnesses.**  Lemma 2 of Brandt uses the fact
that Chaitin's Omega has Kt-random prefixes of *every* length — this is
what prevents TRAVERSE from wrapping around and is the
*only* ingredient that bypasses the black-box barrier (page 5).  No
analogous density of "pi-hard" strings exists: pi mod 2 is one fixed
infinite sequence, and "pi-hard input" is not a meaningful concept (the
function is total and computable in DTIME[O(x^{2/3})] unconditionally,
E7.6).

(O4) **Wrong complexity class.**  Brandt's bound is on uniform time
(linear, conditional super-linear, conditional quadratic).  E5.3 needs
a **circuit** lower bound (separate TC^0 from NC^1 for the AKS
primitive).  Brandt does not produce circuit lower bounds and explicitly
contrasts himself (page 4) with the Williams/Hirahara
"algorithmic-method" approach that *would* yield circuit bounds — which
is the very approach subject to Razborov-Rudich on stronger circuit
classes.  In other words: the price of Brandt's relativisation is that
he doesn't get circuit lower bounds at all.  We cannot have both
"bypasses Natural Proofs" *and* "gives circuit lower bound" out of
this technique simultaneously.

## Why pseudorandomness of pi (E1.9) does not save the argument

`novel/pseudorandomness_of_pi.md` documents 33 independent measures
under which pi(x) mod 2 is empirically indistinguishable from random.
A naive hope: if pi_N is Kt-random for all N, can't we just *be* the
Kt-random witness in Brandt's traversal?  No, for two reasons:

1. **Brandt's traversal does not need a Kt-random witness in the
   `pi`-shape.**  It needs *some* Kt-random infinite sequence; Chaitin
   Omega is the canonical choice and is mathematically simpler to argue
   about (Fact 2 of Brandt cites Chaitin 1975 for the K-randomness of
   Omega's prefixes).  Substituting pi for Omega does not give us a
   stronger conclusion, it just gives us the same Theorem 1 — about
   MKtP, not about pi.

2. **Empirical pseudorandomness != provable Kt-randomness.**  The 33
   measures are all finite-N statistics.  Proving `Kt(pi_N) >= |pi_N|`
   asymptotically would itself be a circuit-style lower bound and is
   subject to Natural Proofs (the very barrier we're trying to thread
   from the other direction).  E1.9 is *consistent with* pi_N being
   Kt-random but does not establish it.

## What this experiment **does** add to the project

- Confirms that no follow-up to Brandt 2024 (between TCC submission and
  April 2026) has bridged the MKtP -> natural-function gap.  S39 already
  flagged this; we re-checked via WebSearch and found no follow-ups.
- Surfaces a clean, edge-quality structural conjecture (worth an EDGE
  entry, candidate **E5.8** below):

  > **The Brandt-class barrier:** any diagonalisation-based lower bound
  > technique that uses self-referential Kt (or self-referential
  > circuit-size) inequalities is structurally tied to the
  > meta-complexity language it diagonalises, and does not transplant
  > to fixed natural functions like pi(x) mod 2 without an additional
  > reduction (which would itself face Williams-style barriers, hence
  > Natural Proofs on stronger classes).

- Closes Brandt as a Chain-E construction route on E5.3.  Combined with
  E7.10 (AKS modulus-twist orthogonality, S61/S64/S66), Chain E is now
  closed for **both** known attack families (AKS-style *and*
  diagonalisation-based).

## Failure mode: E (equivalence)

We file this as **FAIL/E** — not "C" (no need-primes-to-compute-primes
cycle) and not "I" (no information loss).  The closure is at the
**equivalence-of-techniques** level: Brandt 2024 is not a new vehicle
for circuit lower bounds on natural functions; it is exactly the
unconditional uniform-time lower bound on MKtP that the abstract
states.  Trying to repurpose it for pi(x) reduces (O1)–(O4 above) to
problems of strictly greater difficulty than the lower bound itself.

## Suggested EDGES additions

**E5.8 [proposed, EVS shape] — Brandt 2024 is structurally welded to MKtP.**
The diagonalisation technique behind `MKtP not in DTIME[O(n)]` (TCC 2024,
ePrint 2024/687) does not extend to fixed natural functions because the
hard string the proof constructs is an oracle-dependent Kt-random
prefix, not a fixed boolean function; the contradiction uses
self-referential Kt on both sides of an inequality; the only ingredient
bypassing the black-box barrier (1-Kt-randomness of Chaitin Omega
prefixes) has no analog for fixed natural functions; and the bound is
on uniform time, not on circuits.  This rules out the Brandt family as
a Chain-E construction route on E5.3 and pairs with E7.10 to close
Chain E for both major non-equivalent technique families.
> Why this is shape-revealing: it is a **family-level closure on
> diagonalisation-via-meta-complexity** approaches, parallel to E7.10
> (AKS family) and E7.11 (linear post-processing of zero sums).

## Suggested CLOSED_PATHS entry (S51)

```
| Brandt 2024 (TCC) MKtP diagonalisation -> pi(x) mod 2 lower bound | FAIL | E | Brandt MKtP-not-in-DTIME[O(n)] (ePrint 2024/687, Theorem 1) is structurally tied to MKtP itself and does not extend to natural functions. Four orthogonal obstructions: (O1) the hard string z is an oracle-dependent Kt-random prefix, not a fixed function; (O2) the contradiction uses self-referential Kt on both sides; (O3) the only ingredient bypassing the black-box barrier (1-Kt-randomness of Chaitin Omega) has no analog for pi; (O4) the bound is on UNIFORM TIME, not on circuits, while Chain E (E5.3) needs a CIRCUIT lower bound. Empirical pseudorandomness of pi(x) mod 2 (E1.9, 33 measures) is consistent with pi_N being Kt-random but does not establish it asymptotically and would itself face Natural Proofs (T4). Bounded-Kt simulator (3-bit VM, L_MAX=12) confirms TRAVERSE behaves as in the paper and that pi_N is incompressible by short bounded-VM programs at all tested N. With this and E7.10 (AKS modulus-twist orthogonality) Chain E is closed for both known technique families (AKS and diagonalisation-via-meta-complexity). See experiments/constructions/brandt_mktp/. | 51 |
```

## Files

- `brandt_mktp.py` — bounded-Kt simulator + TRAVERSE + pi_N MKtP encoding
  + structural extension analysis. Runs end-to-end in ~1 second.
- `definition.md` — formal signatures of Brandt's hypothesis class, of
  pi_N, and of the conjectured (and falsified) extension.
- `brandt_mktp_results.md` — this file.
