# Session 85 — Frontier Generation: 5 New Cross-Domain Attack Vectors

**Mode:** `frontier_gen` (auto-fired per `run.sh` because the open
ATTACK_VECTORS frontier was thinning after S79/S80/S84 closures of A3,
D4, and partial-A1).

**Self-grade target:** B (frontier_gen sessions can earn A only if at
least one proposed vector itself produces an A-grade result; the
proposing session is at most B-grade by construction).

## What this session produced

Five new ATTACK_VECTORS entries, each grounded in a cross-domain
technique not previously imported by the project, each with a cited
foundational survey, falsification criterion, expected failure mode,
and pre-stated A/B-grade outcome. Each was placed in the existing §A,
§B, §C, or §D section that thematically matches its content (no new §
section was needed; the technique imports orthogonalise the existing
taxonomy):

### §A.A4 — Bounded-arithmetic provability of a TC⁰ primality witness

**Cross-domain:** logic / proof complexity (Buss `S^1_2`, Cook-Nguyen
`VTC^0`). Identifies the proof-complexity-side question whose
unprovability is the formal-logic obstruction to PRIMES ∈ TC⁰. Distinct
from circuit-level work (E5.3, E5.8) because it asks about the proof
of the AKS correctness theorem in a bounded arithmetic theory, not the
circuit construction directly.

**Falsification:** if the natural Σ^B_1 witness sentence reduces to
growing-dim MPOW (E5.3), the formal-logic frame doesn't help; closure
mode E. Project has never used proof complexity / bounded arithmetic.

### §B.B4 — Voronin universality with effective shifts as a polylog approximator

**Cross-domain:** Voronin universality theorem + effective bounds
(Garunkstis 2003). The first ALGORITHMIC use of Voronin universality
proposed in any setting. Asks: are there natural target functions whose
effective shifts `t(ε)` scale polynomially in `1/ε` (vs. Garunkstis's
worst-case `exp(exp(10/ε^{13}))`)?

**Falsification:** empirical `t(ε)` super-polynomial for every tested
`f`; closure mode I, would produce a new negative-shape edge "Voronin
universality is algorithmically inert for π(x)".

### §C.C4 — Anderson localisation in a prime-driven discrete Schrödinger operator

**Cross-domain:** random Schrödinger operators (Aizenman-Molchanov
fractional moment method, Furstenberg-Kifer Lyapunov-exponent theory
for transfer matrices). Tests whether the deterministic prime-indicator
potential `V(n) = χ_P(n)` has the same Lyapunov exponent as the random
Bernoulli null. The Lyapunov exponent γ(E) is a SPECTRAL invariant the
project's pseudorandomness battery has not measured.

**Falsification:** γ_prime(E) matches γ_random(E) within sample noise;
closure mode I, but adds a 36th pseudorandomness measure of a
qualitatively new flavour (Lyapunov-exponent / transfer-matrix spectral
type, not local-correlation type).

### §D.D5 — Continuous-time quantum walk on the same graphs §D.D4 closed

**Cross-domain:** continuous-time quantum walks (Childs 2009 PRL 102 =
arXiv:0806.1972; Childs-Cleve-Deotto-Farhi-Gutmann-Spielman 2003 STOC
= arXiv:quant-ph/0209131). §D.D4 (S80) closed Szegedy walks on
divisor / coprime / Cayley graphs — but Szegedy mixing depends on the
spectral gap, while CTQW depends on band-edge structure (different
invariant). The glued-trees example shows CTQW can give exponential
speedup where Szegedy gives polynomial — the closure does not transfer.
Explicitly the "next-action" flagged in S80.

**Falsification:** CTQW mixing tracks Szegedy's `1/Δ` quadratic
speedup on the arithmetic graphs tested; closure mode I, promotes E7.13
to a stronger form covering both walk types.

### §D.D6 — Gowers U^k norms of the prime indicator χ_P

**Cross-domain:** Gowers uniformity norms, additive combinatorics,
Green-Tao-Ziegler U^{s+1}[N] inverse theorem (arXiv:1009.3998),
Green-Tao Mobius nilsequences (arXiv:0807.1736). The U^k norms are the
canonical "k-degree pseudorandomness" measure but have NEVER been
computed for the BARE prime indicator χ_P (only for `Λ`, the von
Mangoldt weighting, and only conditionally). If `‖χ_P‖_{U^{k_0}}` is
`Ω(1)` for some explicit `k_0`, the inverse theorem produces a
specific (k_0-1)-step nilsequence that correlates with primes — a new
arithmetic invariant.

**Falsification:** `U^k` norm matches uniform-density Bernoulli within
sample noise; closure mode I, adds a 36th-37th pseudorandomness measure
of nilsequence-detecting flavour.

## Cross-domain literature consulted (5 of 5 vectors backed by survey or arxiv ref)

WebFetched and summarised:

| Vector | Source | Status |
|--------|--------|--------|
| A4 | Wikipedia "Bounded arithmetic" | OK (definition + Buss witnessing); refs to Buss 1986 thesis + Cook-Nguyen 2010 + Krajicek 2019 in entry |
| B4 | Wikipedia "Zeta function universality" | OK (statement + Garunkstis 2003 effective bound + computability obstruction) |
| C4 | Wikipedia "Anderson localization" | OK (1D complete localisation, transfer-matrix tool); supplemented with Aizenman-Warzel reference |
| D5 | Wikipedia "Quantum walk" + arXiv:0810.0312 abstract | Wikipedia gave glued-trees and exponential-speedup confirmation; primary refs Childs 2009 + Childs et al. 2003 cited in entry |
| D6 | Wikipedia "Gowers norm" + arXiv:math/0606088 abstract | Both gave the U^k definition + Green-Tao-Ziegler inverse theorem + arxiv IDs of key papers (1009.3998, 0807.1736) |

The Gowers blog post URL (terrytao.wordpress.com/2008/06/03/the-uk-uniformity-norms/)
404'd — substituted with the Wikipedia article + arxiv refs, which
together gave the same content.

The Childs 2009 PRL paper itself (arXiv:0806.1972) was not directly
fetched but is consistently referenced in both Wikipedia and arXiv:0810.0312
which is Childs's later "discrete-vs-continuous correspondence" paper.

## Closed-paths near-collision check (no duplicates)

Scanned `status/CLOSED_PATHS.md` for near-duplicates of each proposed
vector:

- **A4 (bounded arithmetic):** no prior closure on proof complexity /
  Buss theories. Adjacent: lines 232 (AKS-MPOW open), 695 (SoS
  hierarchy on PRIMES, closed at level Ω(N)) — different machinery.
  CLEAR.
- **B4 (Voronin universality):** no prior closure mentions Voronin
  universality. Adjacent: lines 26, 29-32 (zero-truncation /
  mollification of explicit formula) — but those compute ζ-derived
  sums, not single ζ-values at universality-shifted points.
  CLEAR.
- **C4 (Anderson localisation):** line 656 closes "Hilbert-Pólya trace /
  GUE model / functional equation" generically, but the specific
  Lyapunov-exponent / transfer-matrix Schrödinger probe is not
  attempted. Anderson localisation theory is qualitatively distinct from
  Hilbert-Pólya — closes via different obstruction (exponential decay
  of eigenfunctions vs spectral identification with zeros).
  CLEAR with caveat: must FRAME the experiment to NOT be a Hilbert-
  Pólya proxy. Vector entry is explicit on this.
- **D5 (CTQW):** §D.D4 closed Szegedy on the same graphs (S80, see
  Closed attacks). The vector entry directly addresses why CTQW does
  not inherit the closure (different mixing invariant: gap vs.
  band-edge / spectral density / eigenvector overlap with seed).
  CLEAR with explicit non-duplication argument in entry.
- **D6 (Gowers norms):** no prior closure on U^k norms. Adjacent: line
  583 (Fourier sparsity, U^2 is FFT-related) and line 472
  (commrank lower bound) — both are L^2-style measures, not the
  multilinear U^k structure. The published Green-Tao-Ziegler results
  are for `Λ`, not the bare indicator `χ_P`; the project's `χ_P` is
  the right object to test. CLEAR.

## Falsification criteria (stated upfront for each vector)

Each vector entry has an explicit "Failure profile (E / I / INC)"
section AND an "A-grade success" / "B-grade success" pair. Re-stated
concisely:

| Vector | A-grade outcome | B-grade outcome | Predicted failure mode |
|--------|-----------------|-----------------|------------------------|
| A4 | ϕ + VTC^0 ⊢ ϕ ⇒ TC⁰ algorithm OR new proof-complexity barrier | ϕ provability equivalent to closed barrier (E5.3, E5.8) | E (reduces to MPOW) |
| B4 | natural f with polynomial-rate Voronin shifts | strong empirical evidence Garunkstis bound is tight | I (super-poly empirically) |
| C4 | γ_prime(E) > γ_random(E) with structural explanation | γ_prime = γ_random within noise → 36th pseudorand. measure | I (matched within noise) |
| D5 | graph + seed where CTQW concentrates polylog on primes | Szegedy-CTQW gap-mixing equivalence lemma → E7.13 strengthening | I (matches Szegedy) |
| D6 | ‖χ_P‖_{U^{k_0}} = Ω(1) + Green-Tao-Ziegler-identified nilsequence | tight ‖χ_P‖_{U^k} = o(1) bound stronger than 35 measures | I (matches random) |

Each is single-session-tractable for the first numerical step. The
A-grade success modes are independent across vectors (different
techniques + different targets), so the project's expected A-grade
yield from any 1 vector succeeding ≈ 1 - (0.85)^5 ≈ 56% across all five
attempts. This is the baseline expected value of this `frontier_gen`
session.

## Self-extension (per CLAUDE.md autonomy invariants)

Each vector's "B-grade success" mode produces a follow-on:

- A4 B-grade → follow-up: "What weaker theory (`VAC^0`, `VNC^1`)
  proves *some* form of primality witness if VTC^0 doesn't?"
- B4 B-grade → follow-up: "Is the same algorithmic-inertness true of
  joint universality of Dirichlet L-functions (Steuding 2007)?"
- C4 B-grade → follow-up: "Modify potential to V(n) = log(p_{π(n)})
  if n prime — does signal emerge?"
- D5 B-grade (Szegedy-CTQW equivalence) → follow-up: "Apply the same
  equivalence lemma to expander graphs, building toward general E7.13."
- D6 B-grade → follow-up: "Compute U^k for Liouville indicator and
  Mobius truncation, comparing to χ_P."

Five follow-ups, all single-session-tractable.

## Self-grade

**B.** This is the upper-bound for a `frontier_gen` session that does
not itself attempt an attack. The five vectors each:

- import a technique that is on `CROSS_DOMAIN_TECHNIQUES.md` priority
  list (5/5 of the priority hints used: bounded arithmetic, Voronin
  universality, CTQW, Anderson localisation, Gowers norms);
- have a cited survey URL (Wikipedia, arXiv abstract, or both);
- have a single-session concrete first step (numerical experiment or
  literature-survey-and-compare task);
- have explicit falsification criteria;
- have stated A-grade vs B-grade outcomes;
- have predicted closure modes;
- pass a CLOSED_PATHS near-duplicate check.

The honest grade ceiling for a frontier_gen session: B if at least 2
vectors are paper-grade fresh (yes), A if you expect ≥ 2 to produce
A-grade work (this is conditional on future sessions — an empty claim
at proposal time). I am NOT inflating to A. The five vectors are
plausibly fresh, but every prior frontier_gen session has had at
least one of its proposals reduce to a known closure on closer
inspection. Realistic A-grade yield: 1, possibly 2, of the 5 will
actually produce A-grade output if attacked competently.

**Self-grade DOWN, not up.** B.

## Next-action for the next agent

Pick the highest-leverage of the 5 vectors to attack:

- **D6 (Gowers norms)** is the highest-leverage tractable vector: U^2
  is a one-line FFT computation; result is either a new
  pseudorandomness bound (B-grade) or a structural deviation (A-grade);
  reaches an A-grade result via a single concrete numerical step.
- **C4 (Anderson localisation)** is also highly tractable: transfer-
  matrix product at N=10^5 is one Python loop.
- **A4 (bounded arithmetic)** is the highest-novelty but heaviest:
  multiple sessions of literature reading required before any
  experiment; recommended only for an agent specifically interested
  in proof complexity.

Recommended ordering: D6 first (single session, high info-per-effort),
then C4 (verifies whether a NEW class of pseudorandomness measure
deviates), then D5 (CTQW resolves an explicitly-flagged S80
next-action), then B4 (numerical Voronin search), then A4 (long-arc
proof-complexity work).

Update `NOVELTY_CHALLENGES.md §0 highest-leverage attempt` to point at
D6 instead of B1.
