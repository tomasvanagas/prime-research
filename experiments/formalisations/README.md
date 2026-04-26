# Formalisations — Lean 4 proofs of edges and novel results

This directory holds **machine-checked proofs** of the project's main
mathematical claims. Each formalisation produces a permanent verifiable
artifact that survives across sessions, models, and any informal-argument
drift.

Created S68+ as part of the project's shift from disciplined-closure
mode to novelty-production mode. See `CLAUDE.md` § "The Novelty Bar"
and `NOVELTY_CHALLENGES.md` §3.

## Layout

Each formalisation lives in its own subdirectory:

```
experiments/formalisations/<edge_id_or_result>/
    <name>.lean            # the Lean 4 source — must type-check
    <name>_notes.md        # informal proof sketch + which edge/result
                           # this formalises + status (statement-only,
                           # skeleton-with-sorries, complete, etc.)
    lakefile.lean          # if needed for the subdirectory
    lean-toolchain         # if pinned
```

## Why Lean

- **Forces precision.** Vague claims fail to type-check.
- **Surfaces hidden structure.** Proving an edge formally often reveals
  ingredients that informal arguments hand-wave over.
- **Permanent artifact.** Lean proofs survive prompt drift, model
  upgrades, and informal-argument rewrites.
- **Connects to AI mathematician tooling.** AlphaProof, mathlib4,
  LeanDojo all consume Lean. A Lean library of π(x) facts is the
  project's natural interface to those tools.

## Discipline

A formalisation is **complete** only when `lake build` (or
`lean --run`) succeeds with no `sorry` and no `axiom` introductions
beyond mathlib4's standard.

A formalisation that has `sorry` placeholders is **in progress** —
update `RESEARCH_AGENDA.md` Arc 2 with the current state and what
the next agent should do.

A formalisation that fails to type-check at any point is a **bug** —
debug or revert before declaring any milestone complete.

## Active targets

See `NOVELTY_CHALLENGES.md` §3 for the current Lean queue:

| Priority | Target | Notes |
|----------|--------|-------|
| L1 | E2.1 MPS bond-dim identity | Cleanest — start here |
| L2 | E1.5 0.537-bit invariant | Number-theoretic |
| L3 | E2.7 communication rank +2 | Combinatorial |
| L4 | E5.1 BPSW ⇒ TC⁰ | Circuit-class library work |
| L5 | E5.8 Brandt 4-obstruction | Most ambitious |

## First-time setup (when this directory is empty)

1. Verify Lean 4 + mathlib4 are installable:
   ```
   curl https://raw.githubusercontent.com/leanprover/elan/master/elan-init.sh -sSf | sh
   ```
2. Initialise a project: `lake new prime_research`
3. Add mathlib4 as a dependency in `lakefile.lean`.
4. Pick L1 from `NOVELTY_CHALLENGES.md`. Start with theorem statement
   only, then proof skeleton with `sorry`, then fill in.
5. Save in-progress state to `RESEARCH_AGENDA.md` Arc 2 milestones.
