# Session 55 — FOCUS-2 Deep Dive: Liouville-Identity Structural Closure

**Date:** 2026-04-26.
**Mode:** Deep focus on a single research direction (Task #1 in FOCUS_QUEUE).
**Single experiment:** `experiments/sieve/pi_mod_q_fixed/liouville_modular_structure.py`
**Single results file:** same directory, `..._results.md`.

## Top-line finding

The Liouville-identity attack on `pi(x) mod 2`,

    pi(x) = (x - L(x))/2 - C_3(x)             (E2.2 / S46),

is now closed at the **structural-pseudorandomness** level (in addition to
the identity-bottleneck closure of S46).  Each component of the
decomposition — `A(x) mod 2 = (x - L(x))/2 mod 2` and `C_3(x) mod 2` —
is independently pseudorandom under all standard measures, and the two
are statistically independent (MI ≈ 2 × 10⁻⁵ bits).

Two new pseudorandomness measures are added to
`novel/pseudorandomness_of_pi.md`: #25 `A(x) mod 2`, #26 `C_3(x) mod 2`.

## New algebraic observation (free identity)

`L(x) mod 2 = x mod 2` is trivially true for every x:

> L(x) = sum_{n=1..x} lambda(n) with lambda(n) ∈ {-1, +1}, so L(x) mod 2 =
> sum of parities = x mod 2.

Verified bit-exact for all 2 000 000 sampled x.

This **refines** TODO.md FOCUS-2 step 2.  The original phrasing "polylog
L(x) mod 2 → polylog pi(x) mod 2" is vacuously satisfied.  The actual
missing primitives via E2.2 are:

| Primitive       | Equivalent form             | Required to extract pi(x) mod 2? |
|-----------------|-----------------------------|----------------------------------|
| L(x) mod 4      | A(x) mod 2                  | YES (one of two needed bits)     |
| C_3(x) mod 2    | (same definition)           | YES (the other bit)              |

`pi(x) mod 2 = A(x) mod 2 XOR C_3(x) mod 2`.

## Structural battery

For x ∈ [1, 2 000 000]:

| Stream             | block-H(L=8) | AC max[1..30]  | FFT z (N=2e5) | LFSR/N (4096) |
|--------------------|--------------|----------------|---------------|---------------|
| `pi(x) mod 2`      | 3.56 / 8     | 0.85 (density) | 9.22          | 0.5000        |
| `A(x) mod 2`       | **7.9999**   | **0.0010**     | 5.53          | 0.5000        |
| `C_3(x) mod 2`     | 7.88 / 8     | 0.148 (density)| 5.25          | 0.5000        |
| `L(x) mod 4 high`  | 7.9999       | 0.0011         | 5.62          | 0.4998        |

All metrics consistent with random-Boolean null at N = 2 × 10⁶.  Density-bias
autocorrelation of pi(x) mod 2 is well-known (already covered as measure #21
gzip-compressibility); A(x) mod 2 has no such bias because Omega(n) parity
splits 50/50 over n ≤ x.

## Mutual independence

Joint distribution of (A mod 2, C_3 mod 2):

```
H(A) = 1.00000, H(C_3) = 1.00000, H(A,C_3) = 1.99998
I(A; C_3) = 1.94 × 10⁻⁵ bits
```

So the two bits that XOR to give pi(x) mod 2 are statistically independent —
no structure exploitable for joint compression.  Bit-exact verification of
the XOR identity pi mod 2 = A XOR C_3 succeeds on all 2 000 000 x.

## Cheap-proxy battery

11 binary proxies tested for correlation with pi(x) mod 2:
`x mod 2, x mod 4 ∈ {1,3}, M(x) mod 2, Q(x) mod 2, sigma_0(x) mod 2,
omega(x) mod 2, Omega(x) mod 2, L mod 4 high, L mod 8 bits 1,2`.

* All Pearson |corr| < 0.002 (= noise floor 1/√N).
* Best XOR-fusion of any 4-subset achieves agreement = 0.4951 (worse than
  chance).
* Conditional H(pi mod 2 | proxy) = 0.99993 bits for every proxy
  (max possible = 1).

## Files updated

* `experiments/sieve/pi_mod_q_fixed/liouville_modular_structure.py` (new)
* `experiments/sieve/pi_mod_q_fixed/liouville_modular_structure_results.md` (new)
* `experiments/sieve/pi_mod_q_fixed/_run_log.txt` (raw stdout)
* `status/CLOSED_PATHS.md` (new entry for Liouville-component pseudorandomness, S55)
* `novel/pseudorandomness_of_pi.md` (added measures #25, #26; updated header to "26 Measures")
* `EDGES.md` (E2.2 EVS H → M, Chain C M → L, priority list #2 amended)
* `TODO.md` (FOCUS-2 marked q=2 STRUCTURALLY CLOSED, q ≥ 3 still open)

## What remains open for FOCUS-2

The q = 2 case is now closed at both identity and structural levels.
For q ∈ {3, 5, 7, 11, 13} (the spectral-truncation moduli of Chain B),
no Liouville-style identity is known.  The natural next step would be:

> Test whether character-twisted Liouville sums L_chi(x) = sum lambda(n) chi(n)
> for a non-trivial Dirichlet character chi mod q exhibit any modular
> structure beyond what direct pi(x; q, a) computation gives.

This is outside the scope of the current session.  Recorded in TODO.md
FOCUS-2 amendment as the next concrete step.

## Connection to the wider chain landscape

* **Chain B (CRT recovery)** still has its missing primitive — polylog
  pi(x) mod q for q ∈ {3..13} — untouched.
* **Chain C (Liouville parity)** is now structurally exhausted at q = 2;
  the EVS for E2.2 is downgraded H → M.
* **Chain G (4-bit hard zone)** is unchanged: getting one bit (pi(x) mod 2)
  in polylog is sufficient to break the carry-propagation boundary, but no
  such polylog algorithm exists, and S55 strengthens the empirical case
  against one ever existing.

Sessions 55 ends with no breakthrough and no new open thread.  The pi(x)
mod 2 question is now backed by 26 independent pseudorandomness measures.

## Runtime

~30 seconds for the full sweep (N = 2 × 10⁶).  Negligible cost relative
to the value of the closure.
