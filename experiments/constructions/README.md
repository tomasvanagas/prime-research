# Constructions

This directory holds purpose-built mathematical objects: custom circuits,
algebraic structures (rings, modules, groups), transforms, representations,
and algorithms invented in-project rather than imported from literature.

Created S60 to rebalance the project away from measurement-only work. See
`CLAUDE.md` § "Construction is encouraged, not just measurement" for the
rules.

## Layout

Each construction lives in its own subdirectory:

```
experiments/constructions/<descriptive_name>/
    <name>.py              # runnable code that builds + evaluates the object
    <name>_results.md      # what it does, what it doesn't, verdict (C/E/I or
                           # construction-incoherent)
    definition.md          # OPTIONAL: signature + intended relationship to pi(x)
                           # required only for algebraic objects, not for circuits
```

## Bar for entry

- **MUST run.** A theoretical-only essay is not a construction; turn it into
  code or save it under `archive/ephemeral/`.
- **Failed constructions ARE valuable** — file them with the same discipline
  as any other experiment and add a CLOSED_PATHS entry.
- **Genuinely-novel-and-working objects** also get a `novel/<name>.md` entry
  describing the object's significance.

## Currently empty

No constructions yet. The three FOCUS-1 sub-attacks in `TODO.md`
(Bernstein 2003 smaller-r AKS, non-cyclotomic ring AKS, Healy-Viola
Frobenius transplant) are the highest-leverage candidates to build first.
