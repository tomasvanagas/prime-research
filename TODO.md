# TODO

> **For agents starting a session:**
>
> 1. Read `NOVELTY_CHALLENGES.md` — pick a target.
> 2. Read `RESEARCH_AGENDA.md` — see if any in-flight arc fits.
> 3. Read `EDGES.md` — ground your work in existing facts. Cite IDs.
> 4. Check `status/OPEN_PROBLEMS.md` if the target involves circuit
>    complexity of π(x).
> 5. Then handle the housekeeping below if it's been left undone.

---

# === ACTIVE PRIORITY ===

The project's active research targets now live in two files:

- **`NOVELTY_CHALLENGES.md`** — composition challenges, frame-shift
  questions, Lean formalisation queue, negative-shape conjectures,
  synthesis targets. PICK ONE PER SESSION.
- **`RESEARCH_AGENDA.md`** — multi-session arcs that survive across
  sessions. CONTINUE AN ARC IF ONE IS IN-FLIGHT.

Sessions are evaluated on novelty production (per CLAUDE.md). A session
that produces only DUPLICATE-PLUS closures of fresh-perspective brainstorms
is a **failed session** at this stage, even if every workflow rule was
followed.

---

# === RECURRING ===

## Literature watch (every 2-3 weeks, lightweight)

Last run: S58 (window 2026-04-05 → 2026-04-26, NO-DELTA).
See `archive/sessions/session58_literature_watch.md`.

Update `literature/state_of_art_2026.md` only if the asymptotic landscape
changes.

Sources to scan:
- arXiv math.NT recent submissions for "pi(x)", "p_n", "zeta zero"
- ECCC TR2026 for TC⁰/NC¹ separations and Brandt MKtP follow-ups
- Connes / Yakaboylu / van Ittersum / Ono author streams
- primecount / kim-walisch GitHub for new releases
- **Unconditional level-of-distribution past x^{5/8}** — Pascadi 2025
  (E3.12) reached x^{5/8} = 0.625 unconditionally, within 0.04 of the
  algorithmic Meissel-Lehmer threshold x^{2/3} ≈ 0.667. Any unconditional
  improvement past x^{2/3+ε} would convert several conditional π(x)
  algorithms (Lagarias-Odlyzko under GRH, etc.) into unconditional ones.
  Watch for: Pascadi follow-ups, Maynard streams, Bombieri-Vinogradov
  improvements, Friedlander-Iwaniec extensions.
- **AlphaProof / FunSearch / mathematics-AI tooling** — track tools that
  could attack circuit-construction or mathematical-search subproblems.

Note: most months produce zero algorithmic deltas. Absence of news is
itself information.

---

# === HOUSEKEEPING (non-blocking, human review required) ===

## Flag duplicate scripts

Rule (CLAUDE.md File Placement) says no duplicate scripts. Suspects
flagged S39, awaiting human review (do NOT delete without human approval —
each pair has a companion `_results.md`):

- [ ] `experiments/analytic/weil_optimized.py` vs `weil_optimized_v2.py` vs `weil_optimized_v3.py` (v3 same size as v1, possible revert)
- [ ] `experiments/sieve/ht_transfer_attempt.py` vs `ht_signed_transfer_v2.py`
- [ ] `experiments/circuit_complexity/k_party_nof.py` vs `k_party_nof_v2.py`
- [ ] `experiments/circuit_complexity/approx_degree.py` + `approx_degree_quick.py` + `approx_degree_small.py` + `approx_degree_prime.py` + `approx_degree_counting.py`
- [ ] `experiments/information_theory/information_shortcut.py` vs `information_shortcut_v2.py`
- [ ] `experiments/other/breakthrough_attempt_v2.py` (no v1 exists)
- [ ] `experiments/other/self_correcting_v2.py` (no v1 exists)

## Orphan results.md check

Run before halting any session:
```
find experiments/ -name "*.py" | while read f; do
  r="${f%.py}_results.md"; [ ! -f "$r" ] && echo "MISSING: $r"
done
```

Known orphan as of S57: `experiments/wildcard/reservoir_delta_session58.py`
(closed in critique50 as DUPLICATE of line 683 but no companion results.md).
