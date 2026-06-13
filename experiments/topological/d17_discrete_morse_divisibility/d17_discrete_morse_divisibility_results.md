# D17 — Discrete Morse Complexity of the Divisibility Hasse Diagram

**Session:** 232 (frontier attack on ATTACK_VECTORS §D.D17).
**Channel:** Erdős — combinatorial extremal counting.
**Cross-domain ingredient:** Forman 2002 *Sém. Lothar. Combin.* 48 ("A user's
guide to discrete Morse theory"); Benedetti–Lutz 2014 *Exp. Math.* 23, 66
(random discrete Morse). Status before session: **PROPOSED, never used**
in `CROSS_DOMAIN_TECHNIQUES.md` row "Discrete Morse theory".
**Grade self-claim:** **B** (ambitious failure; new closed-form identity
characterising the structural obstruction).

## Falsification statement (pre-stated)

We test whether the discrete-Morse 1-skeleton complexity of the
divisibility-poset Hasse diagram on `[1, N]` satisfies

  H_A: `m_0(N) = O(polylog N)` (radical compression — A-grade)

**Falsifier 1:** `m_0(N) = Ω(N)` (linear) — closes as no compression.

**Falsifier 2:** `m_0(N)` matches the Erdős–Rényi `G(N, |E_div|)` baseline
within sample noise — closes as 39th pseudorandomness measure (B-grade).

## Setup

Build the Hasse diagram `H_N` of the divisibility poset `(Z_{≥1} ∩ [1, N], |)`:
vertices `V = [1, N]`, edges `(m, mp)` for each prime `p` and `m` with `mp ≤ N`.

Run **Forman's greedy random elementary collapse** (Benedetti–Lutz 2014):

- Pick a uniformly random degree-1 vertex `v` of the current sub-graph (a
  "free face") and pair it with its unique incident edge `e = (v, u)`.
- Discard `v` and `e`, decrement `deg(u)` by 1.
- Repeat until no degree-1 vertex remains.

Output: `m_0` = #(remaining vertices), `m_1` = #(remaining edges).
Forman's Morse–Euler formula `m_0 − m_1 = χ(H_N)` is verified at every run.

**Control:** Erdős–Rényi `G(N, |E_div(N)|)` random graphs (same vertex count,
same edge count, uniform edges). Same greedy collapse algorithm.

20 seeds for `H_N`, 20 seeds for ER baseline, at `N ∈ {64, 128, 256, 512,
1024, 2048, 4096}` (with N=8192 in follow-up).

## Empirical findings

### F1 — Both falsifiers triggered

| N    | `n_v` | `n_e`  | `χ`    | `m_0(div)` | `m_0(div)/N` | `m_0(ER mean)` | `Δ = m_0(div) − m_0(ER)` |
|------|-------|--------|--------|------------|--------------|-----------------|---------------------------|
| 64   | 64    | 102    | −38    | 54         | 0.844        | 55.1            | −1.1                      |
| 128  | 128   | 224    | −96    | 110        | 0.859        | 114.4           | −4.4                      |
| 256  | 256   | 485    | −229   | 229        | 0.895        | 231.1           | −2.1                      |
| 512  | 512   | 1033   | −521   | 464        | 0.906        | 470.6           | −6.6                      |
| 1024 | 1024  | 2181   | −1157  | 942        | 0.920        | 961.5           | −19.5                     |
| 2048 | 2048  | 4568   | −2520  | 1905       | 0.930        | 1935.2          | −30.2                     |
| 4096 | 4096  | 9515   | −5419  | 3831       | 0.935        | 3905.7          | −74.7                     |

`m_0(div)/N → 1` as `N → ∞` (ratio increases monotonically), and
`m_0(div) ≈ m_0(ER)` within ≤ 2 %. Falsifier 1 (linear scaling) and
Falsifier 2 (matches ER baseline) BOTH triggered. **No polylog
compression.**

### F2 — Greedy random Morse output is *deterministic* on `H_N`

For `N ∈ {64, 256, 1024}`, 200 independent random seeds produced
**exactly one** `(m_0, m_1)` value each — the random-greedy collapse is
completely seed-independent on the divisibility Hasse diagram. By
contrast, the ER baseline shows variance (`m_0` range = 5–30 across
seeds) at every `N`. This is a structural *Morse-rigidity* property of
the divisibility lattice that the matched-density random graph lacks.

This rigidity is consistent with Bjorner's classical *EL-shellability*
of the divisibility poset, but the empirical observation here is
cleaner: greedy *elementary* collapse on the 1-skeleton converges to a
unique critical-cell count regardless of order.

### F3 — Sharp closed-form identity for `collapses(N) = N − m_0(N)`

The collapse count decomposes (verified empirically up to `N = 8192`):

  `collapses(N) = (π(N) − π(N/2)) + Π_pow(N) + ε(N)`

where
- `π(N) − π(N/2)` = #(primes in `(N/2, N]`),
- `Π_pow(N) := #{p^k : p prime, k ≥ 2, N/2 < p^k ≤ N}` = #(strict
  prime powers > `N/2`) = `Σ_{k ≥ 2} (π(N^{1/k}) − π((N/2)^{1/k}))`,
- `ε(N) ∈ {0, 1}` for all `N` measured (the chained-collapse residual).

Verification table:

| N    | `collapses` | `π(N) − π(N/2)` | `Π_pow(N)` | sum | `ε(N)` |
|------|-------------|------------------|-------------|-----|--------|
| 64   | 10          | 7                | 2           | 9   | 1      |
| 128  | 18          | 13               | 4           | 17  | 1      |
| 256  | 27          | 23               | 3           | 26  | 1      |
| 512  | 48          | 43               | 4           | 47  | 1      |
| 1024 | 82          | 75               | 6           | 81  | 1      |
| 2048 | 143         | 137              | 5           | 142 | 1      |
| 4096 | 265         | 255              | 9           | 264 | 1      |
| 8192 | 475         | 464              | 10          | 474 | 1      |

**ε(N) is identically 1 across all tested N**, suggesting a constant
chained-collapse term (1 chained collapse always occurs, likely the
final pairing of a level-2 leaf that emerges after the initial wave).

By PNT, `π(N) − π(N/2) ~ N/(2 ln N)` and `Π_pow(N) ~ (1 − 2^{−1/2})
√N / ln √N`, so

  `m_0^{Morse,greedy}(H_N) = N − N/(2 ln N) − Θ(√N / ln N) + O(1).`

### F4 — Why the failure is structural

The greedy-Morse collapse on the divisibility Hasse 1-skeleton can only
peel **prime powers `p^k ∈ (N/2, N]`** as initial leaves (these are the
unique vertices `v ≥ 2` with `ω(v) = 1` and `2v > N`, so degree
`= ω(v) + π(N/v) = 1 + 0 = 1`). Once peeled, the structure of `H_N`
admits at most one further *chained* collapse (the empirically constant
`ε(N) ≡ 1`), because every remaining vertex `v ≤ N/2` still has its
upward cover edges intact. The collapse cascade halts.

**The discrete-Morse approach REDUCES `π(N)` to a primes-in-interval
count `π(N) − π(N/2)`, which is no easier than the original problem.**
This is a *circular* structural reduction — exactly the failure profile
the attack vector predicted ("primes ARE the critical cells of any
reasonable Morse function on `Δ(N)`"), but with a quantitative refinement:
the failure is sharp at the level of `O(√N/log N)` (the prime-power
correction).

## Implications

### Why this is B-grade

- **Falsifier 1 + Falsifier 2 BOTH triggered** — no polylog compression,
  and the divisibility-Hasse Morse complexity is not even
  asymptotically distinguishable from random-graph Morse complexity at
  the leading order (both `Θ(N)`).
- **NEW negative-shape closure**: the failure has a clean *quantitative*
  arithmetic identity (F3) — the discrete-Morse 1-skeleton complexity
  of `H_N` reduces *exactly* to a primes-in-interval count plus a
  prime-power correction.
- **Order-independence (F2)** is a structural rigidity of the
  divisibility lattice not present in matched-density random graphs.
  This is genuinely cross-domain: discrete Morse meets analytic NT.
- **Same outcome predicted by D14** (cellular sheaves on the same
  poset): both are *order-based topological probes* of the divisibility
  lattice; D14 (S103) saw `dim H^k(F) ~ Bernoulli` and we now see
  `m_k(div) ~ ER baseline ± 2%`. Two structurally orthogonal
  topological probes converge on the same conclusion: **divisibility-
  poset topology is too rigid to compress `χ_P` without already knowing
  primes.**

### What would have been A-grade (and why it's blocked)

A-grade required `m_0 = O(polylog N)` for the OPTIMAL Morse matching.
Joswig–Pfetsch 2006 shows optimal-matching is NP-hard in general; for
divisibility-poset 1-skeleton specifically, the Morse inequality
`m_0 ≥ b_0 = 1` does not rule out a small `m_0`, but our greedy upper
bound `m_0 ≈ 0.92 N` is far from this lower bound. A non-greedy
algorithm could in principle produce smaller `m_0`, but:

(i) The **ω(n) ≥ 2 cycles** in `H_N` (every composite `n` has ≥ 2
    incoming covering relations) prevent simple peeling.
(ii) Even an optimal Morse matching satisfying `m_0 = polylog(N)` would
     have to encode `χ_P` in the dual *cohomology basis*, which would
     itself reduce to evaluating `π(N)` — circular by F3.

This generalises the failure: the *logical content* of the
divisibility-poset Morse number is `π(N) − π(N/2)`, and any improvement
on the greedy upper bound would have to come from *outside* the
divisibility relation (e.g., adding extra structure like residue
classes mod q).

### Cross-domain registry update

`CROSS_DOMAIN_TECHNIQUES.md` row "Discrete Morse theory":
**PROPOSED → USED-E (S232)** with edge ID `E2.X` (TBD on EDGES.md update).

### Successor proposal (1 follow-on)

**D17.b** (squarefree restriction + augmented Hasse). Build the
*squarefree-only* divisibility Hasse diagram `H_N^{sqfree}` (vertices
`{n ≤ N : μ(n) ≠ 0}`); this is the order complex of the boolean
lattice on primes ≤ √N together with squarefree composites. The
classical result that the boolean lattice's order complex is shellable
with single critical cell at the top suggests `H_N^{sqfree}` may have
strikingly different Morse behavior. If `m_0(H_N^{sqfree}) = polylog N`,
this would be a 0.5-A-grade result for the *squarefree-prime indicator*
`μ²(n)` even if not for `χ_P`. **Budget: 1 session.**

## Files

- `d17_discrete_morse_divisibility.py` — main scan (greedy Morse on
  divisibility Hasse + ER baseline, 20 seeds each, 7 values of `N`).
- `d17_discrete_morse_divisibility_data.json` — raw scaling data.
- `d17_followup.py` — determinism test (200 seeds × 3 `N` values) +
  collapse-level breakdown verifying the F3 identity.

## Session-end self-evaluation

1. **What did I produce that was not in the project before this session?**
   - The closed-form identity `collapses(N) = π(N) − π(N/2) + Π_pow(N) + 1`
     for the greedy-random discrete Morse collapse of the divisibility
     Hasse 1-skeleton.
   - The empirical structural observation that the greedy-random Morse
     output is *deterministic* on the divisibility lattice (uniformly
     across 200 random seeds at `N ∈ {64, 256, 1024}`).
   - The first computation of any discrete-Morse invariant of the
     divisibility poset.

2. **What edges did my work compose or cite?**
   - Predicted in ATTACK_VECTORS §D.D17 failure profile (E).
   - Companion to D14 (cellular sheaves on the same poset, S103) — both
     are order-based topological probes that find no compression.
   - Distinct from D2/E2.17 (metric PH on prime gaps, S96) — D17 uses
     order/divisibility structure; D2 used metric on a Takens embedding.

3. **If my session produced only duplicate closures, why?**
   Not a duplicate — the F3 closed-form identity is a *new* arithmetic
   characterization of the failure mode, sharper than the qualitative
   "primes are forced critical" prediction. The deterministic-output
   observation (F2) is also new structural content.

4. **Next-action for the next agent.**
   D17 closed at B (mode E). Next pickup options on the topological /
   combinatorial axis:
   - **D17.b** (squarefree-only Hasse) — see successor proposal above.
   - **D14 → richer stalks** (`F(n) = F_2^{Ω(n)}` instead of `F_2^{χ_P(n)}`) —
     sub-attack flagged in D14 mitigation; never executed.
   - **D17.c** (Δ(N) full order complex including chains of length ≥ 2,
     not just 1-skeleton) — significantly heavier: chain count grows like
     `~N · 2^{Ω(N)}`, but at `N = 256` is still tractable (`~10^4` cells).
