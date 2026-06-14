# dk_sliding_census.py — results (S523)

The **sliding-window** local pattern census `D_k^slide(x)`, the unaligned
companion of the aligned census (`local_pattern_census.py`,
`dk_first_occurrence.py`). It answers the last named width-spectrum
follow-on (PROGRAM.md SECONDARY LINES, "sliding-window variant
unmeasured").

## Object

- **Aligned** census (prior work): windows `[k·m, k·m+k)`, base `n = k·m
  ≡ 0 (mod k)`. Admissibility of an occupancy pattern `S ⊆ {0..k-1}` has
  two conditions: (i) `q|k ⇒` no offset `o∈S` with `q|o`; (ii) `q∤k ⇒
  {o mod q : o∈S} ≠ ℤ_q`. `A_k^aligned = 13, 106, 3573` for `k=8,16,32`.
- **Sliding** census (this script): windows `[n, n+k)` with `n` over ALL
  residues. Condition (i) disappears (no prime divides every base) and
  admissibility collapses to the **classical** Hardy–Littlewood
  prime-constellation criterion: `S` admissible ⇔ for every prime `q≤k`,
  `{o mod q : o∈S} ≠ ℤ_q`.

Constant-memory **segmented sieve** over the window-start axis (a sliding
k-bit window over the prime indicator; tail handled with a `k−1`
overlap), CLI `--k/--xmax/--seg`, checkpoints at powers of two.

## Headline results

### 1. The parity-doubling identity (structural, exact, selftest-locked)

    A_k^slide = 2·A_k^aligned − 1            (k a power of two)
              = 25, 211, 7145   for k = 8, 16, 32.

The `q=2` clause forces every admissible `S` to be **all-even or
all-odd** (consecutive primes share parity). For `k=2^j`:

- the **odd-offset** admissible family is *exactly* the aligned
  admissible set (aligned condition (i) for `q=2|k` already restricts to
  odd offsets, and the remaining conditions coincide), so `#odd =
  A_k^aligned`;
- the **even-offset** family is its image under the reflection `o ↦
  k−1−o`, an affine bijection mod every prime, hence an *admissibility
  automorphism* (it preserves the residue-coverage count `ω(q)` for
  every `q`). So `#even = #odd`, and the two families share only the
  empty pattern.

Therefore `A_k^slide = #even + #odd + 1(empty) = 2·A_k^aligned − 1`.
Measured even/odd split: `3572 + 3572 + 1 = 7145` at `k=32`. The sliding
variant does **not** reach a denser constellation: `max admissible
weight = ρ*(k) = 9` at `k=32`, identical to aligned and parity-symmetric
(even and odd families both attain it). The H-geometry that caps the
aligned census (`ρ*(k) = W(B(k)) = H()` cross-check, S519) caps the
sliding census at the same weight.

### 2. Full saturation, much earlier than aligned

| k | A_k^slide | E_k^slide (exceptional) | D_k^slide = A+E | saturated by |
|---|---|---|---|---|
| 8 | 25 | 4 | 29 | n < 2^9 (adm), 2^18 verified |
| 16 | 211 | 6 | 217 | n < 2^12 |
| 32 | 7145 | 8 | **7153** | **n < 2^27** (x≈1.3×10⁸) |

`D_k^slide(x) → A_k^slide + E_k^slide`, the sliding analogue of the
aligned `D_k → A_k + 1`. **k=32 saturates at n<2^27**, vs the aligned
census needing x≈1.4×10¹¹ ≈ 2^37 — roughly **1000× earlier in x**. Two
reinforcing reasons: (a) sliding samples `k=32×` more windows at the same
x (every n, not just multiples of k); (b) the densest weight-9
constellations are realized *immediately* at n≈10 by the small-prime
cluster itself (e.g. `n=10`, window `[10,42)`, primes
`{11,13,17,19,23,29,31,37,41}` = an admissible 9-tuple), far ahead of
their HL mean schedule.

### 3. The exceptional set `E_k^slide = q*(k)+1`

Sliding windows over small `n` pass over the small primes
`{2,3,5,7,…}`, realizing **inadmissible** patterns (a small prime `q≡0
(mod q)` sits at an occupied offset and the cluster covers `ℤ_q`). At
`k=32` there are exactly **8**, all at `n=0..7`, covering:

| n | covered primes | weight |
|---|---|---|
| 0,1,2 | mod 2,3,5,7 | 11 |
| 3 | mod 3,5,7 | 10 |
| 4,5 | mod 5,7 | 9 |
| 6,7 | mod 7 | 9 |

The last exceptional is at `n = q*(k)`, the largest prime whose residues
the cluster covers within length `k`: `q* = 3,5,7` for `k=8,16,32` (=
`prime(j−1)` for `k=2^j`). Beyond `n=q*`, the prime `q*` leaves the
window and no residue class can be anchored ⇒ every pattern is
admissible. Hence `E_k^slide = q*(k)+1 = 4,6,8`. (Aligned has `E=1`: only
the `m=0` window is exceptional, because alignment forbids the cluster
elsewhere.)

### 4. Classical singular series validated

The HL count test (`--validate-hl`, occurrence counts vs `N_S(x) =
𝔖(S)·x/(ln x)^w`, classical singular series `𝔖(S) = ∏_q
(1−ω(q)/q)/(1−1/q)^w`, no `/k`, no aligned `q|k` factors) at `k=32`,
`x=2^22`: **`obs/pred = 0.994` at weight 6** — a clean validation of the
classical constant, matching the aligned census's 0.980 at weight 6. Low
weights drift below 1 (HL counts the `≥S` superset, we count *exact*
patterns); the densest (w=8,9) overshoot at small x because the
small-prime cluster realizes them ahead of schedule. The twin constant
falls out as `𝔖({0,2}) = 1.3204 = 2C₂` (selftest).

## Base case (the NEXT ACTION's requested cross-check)

`--aligned` restricts the sliding census to `n ≡ 0 (mod k)`. An aligned
window *is* a sliding window anchored at a multiple of `k`, so this
reproduces the aligned realized-pattern set bit-for-bit. Selftest [5]
asserts `sliding | n≡0 mod k == aligned census` at k=8,16,32.

## What would falsify

- an inadmissible pattern realized at **large n** (`n > q*(k)`):
  impossible by plain divisibility — a bug check; none observed;
- an admissible pattern still missing while its HL mean `N_S(x) ≥ 5`
  (P(0)=e⁻⁵<0.7%) and denser cousins appear: breaks HL ordering. The
  monitor flags only this; it raised nothing (full saturation by 2^27);
- `D_k^slide(x)` plateauing strictly below `A_k^slide + E_k^slide`.

## Selftest (11 groups, ~11 s)

1. `A_slide = {8:25,16:211,32:7145}`; even==odd==`A_aligned−1`;
   odd-offset family == aligned admissible set; `A_slide=2·A_aligned−1`.
2. reflection `o↦k−1−o` is an admissibility automorphism (even↔odd).
3. code↔set round-trip.
4. segmented census == direct full sieve == naive trial-division
   (k=6,8,16,32; 3 segment sizes) — the bit-packing and segmentation.
5. **base case:** sliding | n≡0 mod k == aligned census (k=8,16,32).
6. saturation D_8=29 (2^18), all 211 admissible found at k=16 (2^20).
7. exceptional patterns inadmissible (`𝔖=0`), `E={8:4,16:6,32:8}=q*+1`,
   all at `n=0..q*(3,5,7)`.
8. `D_k^slide(x)` monotone non-decreasing.
9. classical `𝔖`: weight-1==1 (not 2 — no aligned `q|k` factor),
   mod-2/mod-3 cover ⇒ 0, twin `{0,2}=2C₂`, mirror-invariant.
10. HL first-occ finite for all admissible; per-weight median
    non-decreasing, strictly increasing for w≥5 (low weights floor —
    they occur immediately).
11. HL count monotone in x; `ρ*(8,16,32)=(3,5,9)`.

## Honest scope

This is a structural measurement of a finite-x census, not progress on
the polylog `π(x)` goal — it sits in the width-spectrum follow-on track.
The `2·A_aligned−1` identity and the `ρ*(k)`-cap are the genuine new
content (the parity reflection was implicit in the aligned `q=2|k`
restriction but never stated for the unaligned window). Full k=32
saturation is a finite empirical fact consistent with the prime k-tuple
conjecture, not a proof. Run log: `run_slide_k32_2e28.log`.
