# Cyclotomic CRT Splitting for AKS-style Polynomial Powering

**Script:** `cyclotomic_crt_splitting.py`
**Session:** 47

## Question

Does the CRT decomposition

  Z_n[x] / (x^r - 1)  ≅  ∏_{d | r}  Z_n[x] / Φ_d(x)

reduce the *maximum* matrix-powering dimension below `r`, enough to give a TC^0
shortcut for AKS-style verification of `(x+a)^n = x^n + a (mod x^r - 1, n)`?

This addresses the OPEN entry "AKS in TC^0 via matrix powering"
(`status/CLOSED_PATHS.md` line 232) — the only path remaining at the
TC^0/NC^1 boundary for our problem.

## Method

For each n in a representative sample (primes near 10^k for k=2..6,
Carmichael numbers 561, 1105, 1729, 2465, 2821, 6601, 8911, plus
2^10..2^12, 65537, 524287):

1. Compute the AKS-prescribed `r`: smallest r > 1 with gcd(n,r)=1 and
   ord_r(n) > floor(log₂ n)².
2. Factor x^r - 1 over Q as ∏_{d|r} Φ_d(x). Each factor Φ_d has degree φ(d).
3. The maximum dimension to power in is max_{d|r} φ(d).
4. Compare max_dim / r. For the cyclotomic shortcut to break the polylog-MPOW
   barrier, this ratio would have to tend to 0; for any improvement at all,
   it would have to be substantially below 1.

## Key Findings

```
         n |      r |  log2(n)^2 |  max_dim |  max/r | #factors |  r prime
       101 |     41 |         36 |       40 |  0.976 |        2 |     True
       561 |     89 |         81 |       88 |  0.989 |        2 |     True
      1024 |    227 |        100 |      226 |  0.996 |        2 |     True
      1729 |    107 |        100 |      106 |  0.991 |        2 |     True
     65537 |    271 |        256 |      270 |  0.996 |        2 |     True
    100001 |    289 |        256 |      272 |  0.941 |        3 |    False
   1000003 |    379 |        361 |      378 |  0.997 |        2 |     True
```

(22 samples total; full table in script output.)

- **r is prime in 21 / 22 sampled cases (95.5%).** Whenever r is prime,
  x^r - 1 = (x-1) · Φ_r(x), so the CRT splitting is just
  Z_n × Z_n[x]/Φ_r(x). The non-trivial factor has dimension r-1 — essentially r.
- **Maximum max_dim / r = 0.997, minimum = 0.941, average = 0.990.**
  Even in the one composite-r case (r = 289 = 17²), the max factor still
  has dimension 272 = φ(289) = 17·16, only 5.9% smaller than r.
- **Why r is almost always prime in AKS:** the algorithm picks the *smallest*
  r passing the order test ord_r(n) > log² n. Among small integers, primes
  give large multiplicative order (φ(r) = r-1 is a huge cyclic group), so
  prime r usually wins the search by a wide margin. The asymptotic theorem
  (Lenstra-Pomerance) guarantees r = O(log⁴ n) but does not avoid prime r.

## Implication for TC^0

The Mereghetti–Palano (2000) result places fixed-dimension matrix powering
MPOW_k in TC^0 only for *constant* k. For AKS we have k = r = Θ(log^c n)
**growing** with input size. The cyclotomic CRT decomposition was the most
natural "split into smaller dimensions" attack; this experiment closes it:

  - The split produces ONE factor of dimension ≈ r (the Φ_r factor) plus
    trivial 1-dimensional pieces.
  - max-dim ≈ r, so there is no asymptotic reduction.
  - Constant-factor savings (~6% in worst observed case) do not change
    asymptotic complexity class membership.

This **does NOT close the AKS-in-TC^0 question itself** — that requires
showing growing-dimension MPOW is or is not in TC^0, which is the genuine
TC^0/NC^1 frontier (Krebs-Lange-Reifferscheid 2005 algebraic
characterization). It only closes the cyclotomic-shortcut sub-attack.

## Verdict

**CLOSED** (the cyclotomic-CRT sub-attack on AKS-MPOW)
**Failure Mode:** E (Equivalence) — the CRT decomposition is equivalent in
asymptotic dimension to the original r×r matrix powering.

The parent question "AKS in TC^0 via growing-dimension MPOW" remains OPEN
and unaffected by this result.

## One-Line Summary

Cyclotomic CRT splits Z_n[x]/(x^r - 1) into pieces, but the largest piece
Z_n[x]/Φ_r(x) has dimension r-1 (since AKS-prescribed r is prime in 95%+
of cases), so no dimension reduction below r is achieved.

## CLOSED_PATHS entry to add

```
| Cyclotomic CRT decomposition of AKS ring Z_n[x]/(x^r-1) | FAIL | E |
For AKS-prescribed r (smallest r > 1 with gcd(n,r)=1 and ord_r(n) > log_2(n)^2),
r is empirically prime in 21/22 sampled n (10^2..10^6, including Carmichael
numbers and powers of 2). When r prime, x^r-1 = (x-1)·Phi_r(x), so CRT yields
Z_n × Z_n[x]/Phi_r(x); the non-trivial factor has dimension r-1, no
asymptotic reduction. Worst observed reduction max_dim/r = 0.94 (r=289=17^2).
Cyclotomic shortcut for AKS-MPOW closed; parent OPEN status of "growing-dim
MPOW in TC^0" unaffected. See experiments/circuit_complexity/cyclotomic_crt_splitting.py | 47 |
```
