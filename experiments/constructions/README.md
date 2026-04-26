# Constructions

This directory holds **purpose-built mathematical objects**: custom
circuits, algebraic structures (rings, modules, groups), transforms,
representations, and algorithms invented in-project rather than imported
from literature.

Created S60 to rebalance the project away from measurement-only work.
Promoted to the project's primary mode in the post-S67 mature-state
rebalance. See `CLAUDE.md` § "The Novelty Bar" for the core criterion
and `NOVELTY_CHALLENGES.md` for active construction targets.

## Layout

Each construction lives in its own subdirectory:

```
experiments/constructions/<descriptive_name>/
    <name>.py              # runnable code that builds + evaluates the object
    <name>_results.md      # what it does, what it doesn't, verdict
                           # (C/E/I or "construction-incoherent")
    definition.md          # signature + intended relationship to π(x);
                           # required for algebraic objects, optional for circuits
```

## Bar for entry

- **MUST run.** A theoretical-only essay is not a construction; turn it
  into code or save it under `archive/ephemeral/`.
- **MUST cite EDGES IDs** that the construction composes, uses, or
  contradicts.
- **MUST include a "what would falsify this" statement** in the results
  file. If your construction is non-falsifiable, sharpen it.
- **Failed constructions are valuable** — file them with the same
  discipline as any other experiment and add a CLOSED_PATHS entry.
- **Genuinely-novel-and-working objects** also get a `novel/<name>.md`
  entry describing the object's significance.

## Existing constructions

- **`brandt_mktp/`** (S51) — bounded-Kt simulator + 4-obstruction
  argument that closed FOCUS-3. First successful use of this directory.
  Resulted in E5.8.

## Active targets

The current top construction targets are listed in
`NOVELTY_CHALLENGES.md` §1 "Composition Challenges":

- **C1** — A⊕C₃ bisection × 0.537-bits invariant (compose E1.6 + E1.5)
- **C2** — MPS bond-dim × free-probability moments (compose E2.1 + free probability)
- **C3** — Brandt obstructions × per-bit difficulty (compose E5.8 + E1.3)
- **C4** — Aggarwal × Dusart × BPSW (compose E6.6 + E6.8 + E5.1)
- **C5** — N/2 universality × non-Boolean function (compose E1.4 + E2.5)
- **C6** — Three-pillars × HKM time-space curve (compose E7.7 + E6.7)

Each of these has a one-paragraph spec including expected effort,
falsification criterion, and target subdirectory name.
