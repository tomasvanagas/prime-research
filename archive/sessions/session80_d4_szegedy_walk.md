# Session 80 — Cross-domain D4 frontier attack: Szegedy quantum walks on prime graphs

**Date:** 2026-04-26
**Mode:** Cross-domain attack (frontier ATTACK_VECTORS §D.D4)
**Grade:** B (ambitious cross-domain attempt that failed informatively)
**Status of vector:** CLOSED

---

## What did I produce that was not in the project before this session?

1. **A new experimental artefact:** `experiments/quantum/szegedy_walk_prime_graphs/` —
   four scripts (main, extended sweep, eigenvector inspection, degree-class
   check), all running, with results.md + JSON + log artefacts. Tests Szegedy
   quantum walks (Szegedy 2004 discriminant theorem) on three number-theoretic
   graph families.

2. **A new closed-form analytical observation:** the lazy random walk on the
   coprime graph `C_x = ([1..x], gcd=1)` has stationary prime mass / uniform
   prime mass `→ ζ(2) = π²/6 ≈ 1.6449` (Mertens). Verified empirically at
   x=1000 with deviation −0.022 from the asymptote. Derivation: stationary
   ∝ degree, average degree of vertex n is `x · φ(n)/n + O(2^ω(n))`, so prime
   degree is x while average degree is x · 6/π² (Mertens density), giving
   exact ratio π²/6.

3. **A quantitative quantum-walk extension of E7.12 (S79):** classical Cayley
   `Cay((Z/NZ)*, {2,3,5,inv})` mixing time fits `t_C(N) ~ N^{0.789}` and
   the corresponding Szegedy walk fits `t_Q(N) ~ N^{0.414}`. Quadratic
   speedup empirically realised but neither is polylog. The S79 ω(n)-
   probing barrier extends to quantum walks.

4. **A new EDGES.md entry, E7.13:** "Szegedy walks on arithmetic graphs do
   not yield polylog π(x)" — packages the result above into a single negative-
   shape edge, combining (Cayley: poly mixing) and (coprime/divisor: Ω(1)
   mixing without primality eigenvector) into a joint-impossibility statement.

5. **A new CLOSED_PATHS row** with full data and the joint-incompatibility
   condition stated cleanly.

6. **A new annotation on E7.12** showing that the classical Cayley spectral
   barrier carries to the Szegedy quantisation.

## What edges did my work compose or cite?

- **E5.3** (TC⁰ primality requires growing-dim MPOW)
- **E5.8** (Brandt structural-welding to MKtP)
- **E7.10** (AKS modulus-twist orthogonality)
- **E7.12** (Cayley spectrum probes ω(n), not primality) — directly extended
- **Cross-domain reference:** Szegedy, *Quantum Speed-Up of Markov Chain
  Based Algorithms*, FOCS 2004, arxiv quant-ph/0401053.

## Cross-domain ingredient

**Szegedy 2004 discriminant matrix theorem.** Given a reversible Markov
chain `P` on a graph, define `D(x,y) = sqrt(P(x,y) P(y,x))`. The Szegedy
walk operator `W(P)` (a product of two reflections) has eigenvalues
`e^{±2iθ_k}` where `cos(θ_k)` are the eigenvalues of `D`. Spectral gap
of `W(P)` is `2 arcsin(sqrt(δ))` where `δ` is the gap of `P` —
**quadratic mixing speedup** `O(1/sqrt(δ))` vs classical `O(1/δ)`.

The cross-domain import is genuine: the discriminant matrix construction,
the eigenvalue lifting, and the mixing-time bound are all imported and
used. The empirical confirmation `t_Q ≈ sqrt(t_C)` on the Cayley sweep
(slopes 0.789 vs 0.414, ratio ≈ 0.524 ≈ 0.5) is the real signature
that the cross-domain machinery is doing real work.

## What was the challenge?

ATTACK_VECTORS.md §D.D4: "Define a quantum walk whose stationary
distribution encodes π(x) extraction. Analyse its mixing time."
Concrete first step: divisor graph G_x. Compute Szegedy mixing time
for x ≤ 100. Polylog upper bound?

## Why the attack failed (and the failure is structural)

For a Szegedy walk on an arithmetic graph to give polylog primality
extraction, we need *jointly*:

1. Classical spectral gap `δ = Ω(1/polylog(x))` (so quantum mixing
   `O(sqrt(1/δ)) = O(polylog(x))`).
2. A single discriminant eigenvector with primality-localised mass
   (so projecting onto it after mixing yields a primality signal).

The three test families exhibit exhaustive incompatibility:

- **Cayley:** structure is rich (E7.12 says spectrum encodes ω(n)
  classes, which contain primes) but `δ` shrinks polynomially with N.
  Quantum quadratic speedup helps but doesn't cross the polylog
  threshold.
- **Coprime:** `δ = Ω(1)` (very fast mixing), but the only primality-
  carrying eigenvector is the trivial Perron mode (stationary
  distribution = degree-weighted), which gives a constant `π²/6` bias.
  Constant bias means `O(log x)` samples to hit a prime — no
  amortisation possible.
- **Divisor:** `δ = Ω(1)`, but high-prime-mass eigenvectors localise
  on prime-centered clusters (`{p, 2p, ..., 1}` for individual `p`).
  Each eigenvector picks 1-2 specific primes — same degree-class
  probing as E7.12. Counting all primes via spectrum requires
  diagonalisation = `O(x^ω)`.

The structural lemma: **the joint condition (fast mixing + global
primality eigenvector) is empirically incompatible across these three
canonical arithmetic-graph families.** This is the negative-shape
content that makes E7.13 a real edge.

## What is the next-action for the next agent?

The natural next step would be C-grade refinement: try Szegedy walks
on **non-abelian Cayley graphs** (e.g., `Cay(S_n, transpositions)`),
which by Bourgain-Gamburd-style results can have explicit
`Ω(1/log(n))` spectral gaps. Even there, the discriminant theorem
gives at most `O(sqrt(log) · n)` mixing — no polylog opening, but
the spectral structure is richer. Also worth: random-projection
arithmetic graphs (the N1 follow-on suggested in RESEARCH_AGENDA Arc
4) under Szegedy quantisation.

A more ambitious next step: explore **continuous-time quantum walks**
(`e^{i H t}`) rather than Szegedy walks. The hitting/mixing
phenomenology is different; some graphs exhibit polynomial CTQW
speedups beyond Szegedy.

## Cleanup verification

```
$ find experiments/quantum/szegedy_walk_prime_graphs -name "*.py"
  szegedy_walk_prime_graphs.py
  szegedy_walk_extended.py
  eigenvector_inspection.py
  degree_class_check.py

$ find experiments/quantum/szegedy_walk_prime_graphs -name "*results*.md"
  szegedy_walk_prime_graphs_results.md
```

All Python files have associated artefacts (results.md or .json/.log).
No `__pycache__` left behind. All four scripts type-check and run.
