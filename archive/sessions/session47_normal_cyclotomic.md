# Session 47 — Normal Mode (2026-04-25)

## Summary

Two productive items in steady-state mode: (1) a focused literature scan
covering the 3-week window since the previous full update (2026-04-05 →
2026-04-25), (2) a sharpening experiment on the only OPEN entry in
`status/CLOSED_PATHS.md` — "AKS in TC^0 via matrix powering" (line 232).

## 1. Literature scan (2026-04-05 → 2026-04-25)

Ran 10 topic-targeted WebSearches:
1. Growing-dimension matrix powering, TC^0/NC^1 frontier
2. Iterated matrix product circuit lower bounds
3. PRIMES in TC^0, BPSW, AKS variants
4. New π(x) algorithms (Deléglise-Rivat, Gourdon, primecount)
5. Riemann zeta zero summation, Hiary, Guth-Maynard follow-ups
6. Hardy-Littlewood / Selberg / Helfgott-Thompson sieves
7. Sub-x^{1/2} claims for π(x) or p_n
8. Brandt MKtP diagonalization follow-ups
9. Ono / van Ittersum partition detection extensions
10. Connes / Yakaboylu Hilbert-Pólya updates

Findings:
- **No paper in the window changes the asymptotic barrier.**
- Two minor catalog additions to `literature/state_of_art_2026.md`:
  - Connes arXiv:2602.04022 (Feb 2026 survey + Letter)
  - Yakaboylu arXiv:2408.15135 v10 (March 2026 revision)
- primecount remains at v8.4 (April 2026, the most recent release in the
  prior literature snapshot). No new release in window.
- Aggarwal arXiv:2510.16285 remains the tightest current statement on
  computing p_n.

## 2. Cyclotomic CRT splitting for AKS-MPOW

**Question.** The OPEN entry on AKS in TC^0 hinges on whether r×r matrix
powering in Z_n[x]/(x^r-1), with r = O(polylog n), can be done in TC^0.
The CRT decomposition

    Z_n[x] / (x^r - 1)  ≅  ∏_{d|r}  Z_n[x] / Φ_d(x)

was the most natural "split into smaller dimensions" attack: maybe the
maximum factor dimension is much smaller than r.

**Method.** For each n in a sample of 22 values (primes near 10^k for
k=2..6, Carmichael numbers 561, 1105, 1729, 2465, 2821, 6601, 8911, plus
2^10..2^12, 65537, 524287), compute the AKS-prescribed r and the cyclotomic
factorisation x^r-1 = ∏_{d|r} Φ_d(x). Record max dimension max_{d|r} φ(d).

**Result.** r is **prime in 21/22 cases (95.5%)**. AKS picks the smallest
r passing ord_r(n) > log²n; among small integers, primes give the largest
multiplicative order (φ(r) = r-1 cyclic), so they win the race. When r is
prime, x^r-1 = (x-1)·Φ_r(x) and the non-trivial factor still has dimension
r-1.

| metric             | value |
|---|---|
| Average max_dim/r  | 0.990 |
| Best max_dim/r     | 0.941 (only composite r in sample, 289=17²) |
| Worst max_dim/r    | 0.997 |

**Verdict.** The cyclotomic-CRT shortcut is **CLOSED** as failure mode E
(equivalence to direct r-dim MPOW). No asymptotic dimension reduction.

**Important.** This does **not** close the parent OPEN entry (line 232).
The growing-dimension MPOW question remains open at the TC^0/NC^1 boundary
(Mereghetti-Palano covers fixed k only; Krebs-Lange-Reifferscheid algebraic
characterisation is the relevant frontier). The closure here only rules
out the cyclotomic-decomposition sub-attack on it.

## Files written

- `experiments/circuit_complexity/cyclotomic_crt_splitting.py`
- `experiments/circuit_complexity/cyclotomic_crt_splitting_results.md`
- `status/CLOSED_PATHS.md` (new line 688)
- `status/SESSION_INSIGHTS.md` (new Session 47 entry)
- `literature/state_of_art_2026.md` (Session 47 literature delta paragraph)

## Verdict

Steady-state. Literature delta zero impact. One sub-attack on the only
OPEN frontier closed. The remaining theoretical surface area is small;
sessions of this kind are best used for sharpening rather than discovery.
