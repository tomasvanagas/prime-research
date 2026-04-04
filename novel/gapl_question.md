# Novel Framing: Is pi(x) in GapL?

**Status:** Novel observation connecting known concepts. The specific question
"Is pi(x) in GapL?" does not appear to be explicitly studied in the literature.

## The Question

GapL is the class of functions expressible as the difference of two #L functions,
equivalently as determinants of poly(N)-size matrices with logspace-computable entries.
GapL ⊂ NC^2 ⊂ NC.

**"Is pi(x) in GapL?"** asks: does there exist a poly(N) × poly(N) matrix M(x)
such that:
- M(x)[i][j] is computable in O(log x) space
- det(M(x)) = pi(x) (or a simple function of pi(x))

This is STRICTLY SHARPER than "Is pi(x) in NC?" because GapL ⊂ NC^2 ⊂ NC.

## Why This Matters

A positive answer would give a concrete ALGEBRAIC formula for pi(x): the
determinant of a small matrix. This is the most structured form an efficient
algorithm could take.

A negative answer (pi(x) ∉ GapL) would not rule out NC but would eliminate
the most natural algebraic approach.

## Known Related Results

- **Primality testing IS in GapL**: (n-1)! mod n can be computed in GapL
  (Shamir 1991). So the DECISION problem "is n prime?" has a GapL formulation.
  But COUNTING primes requires summing 2^N such decisions.

- **The Redheffer matrix**: det(R_x) = M(x) (Mertens function), but R_x is
  x × x = 2^N × 2^N. No poly(N)-size version is known.

- **Valiant's theorem**: Expressing the permanent as a determinant requires
  exponential-size matrices (unless #P = FP). But pi(x) is not known to be
  #P-hard, so Valiant's barrier may not apply.

## Concrete Search Target

Find a directed acyclic graph G on poly(N) vertices where:
- Edges have integer weights computable in logspace from x
- The signed count of source-to-sink paths equals pi(x)

By the Lindstrom-Gessel-Viennot lemma, this is equivalent to finding a
poly(N)-size matrix whose determinant is pi(x).

The natural candidate (sieve graph) has 2^{Θ(N)} vertices (the floor-value set).
The question is whether a DIFFERENT graph structure exists.

## Connection to Other Open Problems

| Question | Strength | Status |
|----------|----------|--------|
| Is pi(x) in GapL? | Strongest (GapL ⊂ NC^2) | OPEN |
| Is pi(x) in NC^2? | Strong | OPEN |
| Is pi(x) in NC? | Our target equivalence | OPEN |
| Is pi(x) in P (binary input)? | Weakest | OPEN |

Resolving ANY of these positively would be a major breakthrough.
Resolving ANY negatively would settle the impossibility.

Source: Session 13, theoretical analysis by counting-theory agent
Date: 2026-04-04
