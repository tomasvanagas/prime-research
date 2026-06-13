# D17 follow-up — Determinism test + identity verification

This script is a verification companion to `d17_discrete_morse_divisibility.py`.
Its findings are reported under sections F2 and F3 of
`d17_discrete_morse_divisibility_results.md`.

## What it does

1. **Determinism test** — runs greedy random Morse with 200 distinct
   seeds at `N ∈ {64, 256, 1024}` and counts the number of distinct
   `(m_0, m_1)` outputs. Result: 1 distinct output at every tested `N`.
2. **Chained-collapse breakdown** — runs the collapse and reports
   `collapses(N)`, `π(N)`, `π(N/2)`, `π(N) − π(N/2)`, initial-leaf count,
   primes among initial leaves, and prime powers among initial leaves
   for `N ∈ {64, 128, 256, 512, 1024, 2048, 4096, 8192}`. The output
   verifies the identity `collapses(N) = (π(N) − π(N/2)) + Π_pow(N) + 1`.

## How to run

```
python3 d17_followup.py
```

Output is printed to stdout. No JSON written — the main scaling data
is in `d17_discrete_morse_divisibility_data.json`.

## Cross-references

- All findings from this script are summarised under
  `d17_discrete_morse_divisibility_results.md` §F2 (determinism) and §F3
  (closed-form identity).
- See `archive/sessions/session232_d17_discrete_morse_divisibility.md`
  for the session synthesis and grade.
