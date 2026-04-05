# Focus Queue: Deep Dive Research Tasks

Each focused session works on exactly ONE task from this list.
The task number rotates automatically across sessions.

**NOTE (S37):** All 4 original tasks below are COMPLETED. Do NOT re-run them.

**To add a new task:** Copy the format below (Goal, Experiments, Key question,
Save to) and append it as Task 5, 6, etc. New tasks should target the ONLY
remaining open direction: circuit complexity of pi(x). See `status/OPEN_PROBLEMS.md`.

**If no concrete experiment can be formulated:** Do not invent busywork. Instead,
check `literature/state_of_art_2026.md` for new publications or improve existing
algorithms (see BEST_ALGORITHMS.md comparison table for known gaps vs primecount).

---

## Task 1: Kt Complexity of delta(n) — COMPLETED (Sessions 20, 36)

**Results:** `experiments/information_theory/kt_complexity/SYNTHESIS.md`

Key findings: Kt(delta) ~ 5.58*N (linear, incompressible). Entropy rate 5.8 bits/value.
Transfer entropy n->delta is 0.013 bits (n adds no info). AR(7) captures direct memory.
1/f^1.70 spectrum real but not exploitable. 82% of Fourier modes needed.
**Verdict: delta(n) is incompressible. Attack path CLOSED (S35).**

---

## Task 2: Zeta Zero Structural Patterns — COMPLETED (Session 25)

**Results:** `experiments/analytic/zeta_structure/SESSION_25_SUMMARY.md`

Key findings: No rational relations in pairwise ratios (499,500 tested). No integer
linear relations (PSLQ, 13,000+ tests at 60-digit precision). DFT matches GUE (corr=0.9999).
Equidistributed mod all 10 constants. Sparse matrix model needs O(N) params.
**Verdict: Zeros are GUE-random. All structural approaches CLOSED.**

---

## Task 3: Novel Identity Search — COMPLETED (Session 29)

**Results:** `experiments/algebraic/identity_search/` (5 results files)

Key findings: PSLQ relations fail cross-validation. WZ test fails (differences grow).
LLL finds different minimal polynomials at each x. No ODE, no Volterra kernel.
f(x) = pi(x) - R(x) has no computable identity in any tested basis.
**Verdict: CLOSED. f(x) is algebraically independent of all standard bases.**

---

## Task 4: Conditional Algorithms — COMPLETED (Session 33)

**Results:** `experiments/analytic/conditional/` (4 results files)

Key findings: Under GRH, K_min zeros scales O(sqrt(x)/log(x)) — not polylog. Elliott-
Halberstam shows no info gain from residue classes. GRH-based Miller uses only 2-3 witnesses
but doesn't help counting. Schoenfeld bounds are tighter but sieve cost still O(sqrt(x)*ln(x)).
**Verdict: Even the strongest conjectures don't reduce below O(sqrt(x)).**
