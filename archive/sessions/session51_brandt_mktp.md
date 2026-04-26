# Session 51 — FOCUS-3: Brandt 2024 MKtP-diagonalisation vs pi(x) mod 2

**Date:** 2026-04-26
**Status:** CLOSED (FAIL/E)
**Critical-path item:** FOCUS-3 (last untouched construction-flavoured
attack on the only open problem — circuit complexity of pi(x), per
`status/OPEN_PROBLEMS.md`)
**Edges added:** E5.8

## Headline

Brandt 2024 (TCC, IACR ePrint 2024/687) proves `MKtP ∉ DTIME[O(n)]`
unconditionally via a length-monotonic depth-first traversal that uses
1-Kt-randomness of Chaitin's Ω prefixes to bypass the black-box barrier.
The technique relativizes and threads the Natural Proofs barrier (T4) —
this part is real and a genuinely new capability.

**However it does not extend to natural functions like pi(x) mod 2.**
We close FOCUS-3 with FAIL/E (equivalence-class) and a new edge E5.8.

## Source review

- N. Brandt, "Lower Bounds for Levin–Kolmogorov Complexity," TCC 2024,
  LNCS 15363; IACR ePrint 2024/687.
- Read in full from <https://eprint.iacr.org/2024/687.pdf> (PDF fetched
  via curl, parsed with pypdf).
- Theorems used: Theorem 1 (unconditional linear-time bound), Theorem 2
  (conditional), Corollary 1 (RAM quadratic), Corollary 2 (slightly
  super-linear), Lemma 1 (infinitely-often bound on traversal counts),
  Lemma 2 (TRAVERSE never wraps because Ω prefixes are Kt-random).

## Construction (this session)

Built `experiments/constructions/brandt_mktp/`:

- `brandt_mktp.py` — bounded-Kt simulator (3-bit-op stack VM, L_MAX=12,
  T_MAX=4096) that enumerates all programs and computes
  Kt_bounded(target). Implements TRAVERSE on this bounded oracle and
  shows it descends-or-right-steps as the paper says (PART A). Computes
  bounded Kt of pi_N for N=3..10 (PART B) and encodes pi_N as MKtP
  queries (PART C). Embeds the four-clause structural argument
  (PART D). Runs end-to-end in 1.0 second.
- `brandt_mktp_results.md` — full closure with quoted theorems, the
  four obstructions to extending Brandt to pi(x) mod 2, and the
  suggested EDGES + CLOSED_PATHS entries.
- `definition.md` — formal signatures of MKtP, pi_N, the bounded VM,
  and the (falsified) conjectured extension.

## The four obstructions

(O1) The hard string z that Brandt's TRAVERSE produces is an
oracle-dependent Kt-random prefix, not a fixed natural function.

(O2) The contradiction `Kt(z) >= |z|` and `Kt(z) <= |M| + log_2 t`
uses self-referential Kt on both sides; there is no analog of
"f-hardness(z) >= |z|" for fixed Boolean f.

(O3) Lemma 2's 1-Kt-randomness of Ω prefixes is the only ingredient
that bypasses the black-box barrier (page 5). No density of "pi-hard
inputs" exists for a fixed total function.

(O4) Brandt produces uniform-time lower bounds, not circuit lower
bounds. Chain E (E5.3) needs a circuit lower bound. Brandt explicitly
avoids the Williams/Hirahara algorithmic-method route on page 4 —
precisely the route that would yield circuit bounds and is subject to
Natural Proofs on stronger classes. The price of relativisation is no
circuit lower bound out of this technique.

## Pseudorandomness of pi (E1.9) does not save the argument

`novel/pseudorandomness_of_pi.md`'s 33 measures are consistent with
pi_N being Kt-random but cannot establish it asymptotically — proving
asymptotic Kt-randomness IS a circuit-style lower bound and faces T4.

## Cumulative effect

With E7.10 (AKS modulus-twist orthogonality, S61/S64/S66) closing the
AKS family and E5.8 (this session) closing the Brandt /
diagonalisation-via-meta-complexity family, **Chain E is now closed
for both known technique families** on E5.3.

The remaining unconstrained levers on the only open problem are:
- non-AKS TC⁰ primality tests (no candidate currently identified)
- entirely-new circuit lower-bound techniques

## Files written

- `experiments/constructions/brandt_mktp/brandt_mktp.py`
- `experiments/constructions/brandt_mktp/brandt_mktp_results.md`
- `experiments/constructions/brandt_mktp/definition.md`
- `archive/sessions/session51_brandt_mktp.md` (this file)

## Files modified

- `status/CLOSED_PATHS.md` — added the S51 Brandt closure row at end.
- `EDGES.md` — added E5.8; updated footer.
- `TODO.md` — marked FOCUS-3 as CLOSED (S51, E5.8); removed from
  active critical path.
