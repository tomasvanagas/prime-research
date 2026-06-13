# S425 — Re-verify-closure of E2.28: D_Ω₁ q-reduction identity on H_N

**Mode:** re-verify-closure (adversarial probe of E2.28).
**Target:** S161 closure of E2.28 (Baker-Norine q-reduced form of D_P^N).
**Self-grade:** **B** (case (i) — refinement of E2.28 with a precise new
statement extending the scope of S161's generalised identity to a new
divisor class).

## Adversarial frame

S161 closed the chip-firing route by: (a) the rank invariant `r(D_P)`
matches generic Riemann-Roch for matched-degree random divisors; (b)
"building D_P^N requires knowing primes." The S161 generalised identity
covered divisors **supported on `{1} ∪ {primes}`**. S161 explicitly
tested D_P, D_sqfree, D_μ_pos, D_λ_pos, D_Ω₂; these support classes are
(i) restricted to {1} ∪ {primes} (D_P / D_μ_pos), (ii) cascade through
Dhar (D_λ_pos / D_Ω₂), or (iii) extend to all squarefree (D_sqfree).

**Missed angle: D_Ω₁ (prime-power indicator).** S161 did NOT test the
divisor `D_Ω₁(n) = 1[Ω(n) = 1]` (chips at primes AND at prime powers
p^k with k ≥ 2). This sits between D_P and D_sqfree in support: includes
all prime first-powers (like D_P) PLUS prime higher-powers (which D_P
omits). Question: does the q-reduced form on H_N give a clean identity?

## Probe

`s425_inverse_chipfiring_probe.py` (six divisors × {Γ_N, H_N} × N ∈
{16, 32, 64, 128}, 48 runs) measures sink chips and non-sink support.
`s425_verify_omega1_identity.py` extends to N = 256 and prints full
chip distribution.

## Empirical findings

**Identity S425.1 (refinement of E2.28 statement (1)).** On H_N with
sink q = 1, the q-reduced form of `D_Ω₁` satisfies:

> `D'_Ω₁(1) = π(N)`     (sink chip count)
> `Σ_{v≠1} D'_Ω₁(v) = π(√N) + π(N^{1/3}) + π(N^{1/4}) + ⋯`
> `                  = #{prime powers p^k : k ≥ 2, p^k ≤ N}`     (off-sink total)

**Verified.** N ∈ {16, 32, 64, 128, 256}:

| N   | π(N) | input total | empirical sink | empirical off-sink | predicted off-sink |
|-----|-----:|------------:|---------------:|-------------------:|-------------------:|
|  16 |    6 |          10 |              6 |                  4 |                  4 |
|  32 |   11 |          18 |             11 |                  7 |                  7 |
|  64 |   18 |          27 |             18 |                  9 |                  9 |
| 128 |   31 |          44 |             31 |                 13 |                 13 |
| 256 |   54 |          70 |             54 |                 16 |                 16 |

Sink and off-sink-total identities both hold exactly. **The off-sink
chips are scattered across various non-prime-power vertices** (e.g.,
N=32: chips land at 6, 18, 20, 28 — not just prime powers themselves);
only the TOTAL is invariant.

**Identity S425.2 (refinement of E2.28 statement (2)).** On Γ_N with
sink q = 1:

> `D'_Ω₁(1) = π(N) − π(N/2)`,

identical to D_P's sink value. Verified at N ∈ {16, 32, 64, 128} (both
match exactly).

## Structural mechanism (sharpened closure)

Observe that `D_Ω₁ = D_P + D_higher`, where `D_higher` puts a chip on
each prime power `p^k` with `k ≥ 2`. Q-reduction is linear modulo the
lattice of principal divisors, so

  `D'_Ω₁(1) = D'_P(1) + D'_higher(1)`.

By S161, `D'_P(1) = π(N)`. The empirical identity then forces
`D'_higher(1) = 0`. **Why?** Chips at vertices `p^k` with `k ≥ 2` are
NOT adjacent to `q = 1` in H_N (since H_N's edges are `(a, p·a)` for
`p` prime, so `q ~ v` iff `v` is prime). The q-reduction's
chip-lending phase requires a length-1 path from q to extract chips;
chips at depth ≥ 2 from q are graph-topologically unreachable by the
length-1 lending move and they redistribute among non-sink vertices
without contributing to sink.

**Sharpened closure mechanism (S425):**
- **On H_N (depth-1 collection):** `D'(1) = Σ_{v ~ q} D(v)` for any
  divisor D whose chips at depth ≥ 2 from q are graph-topologically
  isolated (chips at `p^k`, `k ≥ 2`, satisfy this — their only
  down-cover is `p^{k-1}`, blocking single-step lending). Since
  `q`'s depth-1 neighbours in H_N are exactly the primes ≤ N,
  `D'(1)` reads off the **prime-restricted chip total** of D. For
  `D ∈ {D_P, D_Ω₁, D_sqfree, D_μ_pos}` on H_N, sink = π(N), π(N),
  π(N) + 1, π(N) + 1 respectively (the +1 in the sqfree/μ_pos cases
  is `D(1) = 1` since both put a chip on vertex 1 itself).
- **On Γ_N (leaf collection):** sink = #{primes p : N/2 < p ≤ N} for
  any divisor D putting chips on those leaves, with the same
  graph-topological blocking on prime powers `p^k > N/2` having
  multiple proper divisors and hence not being leaves.

**Why this preserves closure.** The depth-1 / leaf graph-topological
extraction mechanism IS what makes π(N) emerge as the sink chip
count. But identifying "depth-1 neighbours of q in H_N" is
**equivalent to enumerating primes ≤ N** (by H_N's construction), and
identifying "leaves of Γ_N" is **equivalent to enumerating primes
in (N/2, N]** (since leaves are exactly degree-1 vertices, which are
primes p with 2p > N). Both reductions are polylog-equivalent to direct
prime enumeration on the relevant set — no algorithmic opening.

## Verdict

**Closure stands**, but with a **precise structural mechanism** that
S161 did not articulate:

- S161 said: "building D_P requires knowing primes" — implicit
  closure on the **input divisor**.
- S425 sharpens to: "the chip-firing extraction mechanism is
  **graph-topological** (depth-1 / leaf collection), and the
  extracted graph-topological quantity (depth-1 neighbour count
  of q on H_N, leaf count on Γ_N) is itself prime-counting" —
  closure at the **graph and extraction-mechanism level**, not just
  the input-divisor level.

D_Ω₁ is a NEW divisor class for which the identity holds without being
"supported on {1} ∪ {primes}" (the S161 generalisation's hypothesis).
The mechanism differs (depth-1 vs leaf), but the closure persists.

## What would falsify this

- Find a polylog-buildable input divisor D on H_N or Γ_N for which
  q-reduction's sink encodes π(N) **without** going through depth-1
  / leaf identification. Such a divisor would have to put chips on
  vertices whose graph-topological role does NOT coincide with
  primes / primes-in-(N/2, N]. **No such divisor is currently known.**
- Find an alternative q-reduction-like operation (multi-sink
  reduction, weighted chip-firing, fractional chip-firing) that
  bypasses the depth-1 / leaf bottleneck.

## Files

- `s425_inverse_chipfiring_probe.py`: 6 divisors × {Γ, H} × 4 N values.
- `s425_verify_omega1_identity.py`: full distribution at N up to 256.
- `s425_inverse_chipfiring_results.json`: raw measurements.

## What this session produced

1. **New identity S425.1** on H_N: `D'_Ω₁(1) = π(N)`,
   off-sink-total = `Σ_{k≥2} π(N^{1/k})`. Not in S161; not in EDGES.md
   prior to this session. Empirically verified on N ∈ {16, 32, 64,
   128, 256}.
2. **New identity S425.2** on Γ_N: `D'_Ω₁(1) = π(N) − π(N/2)`,
   matching D_P. Same verification range.
3. **Sharpened closure formulation** for E2.28:
   depth-1 collection on H_N, leaf collection on Γ_N — the
   chip-firing mechanism is graph-topological and reduces to prime
   detection on the relevant set.
4. **Decomposition lemma** `D'_Ω₁ = D'_P + D'_higher` with
   `D'_higher(q) = 0` (chips at `p^k`, `k ≥ 2`, contribute zero to
   sink).
