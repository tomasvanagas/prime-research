# Session 69 — F2 closure (information-rate side): pi(x) mod 2^k saturates at h_2(pi(X)/X)

**Date:** 2026-04-26
**Run:** 58 -> 59 (run.sh sequence)
**Mode:** Deep focus, Task #4 slot — but Task #4 was already COMPLETED (S33),
so per CLAUDE.md guidance picked NOVELTY_CHALLENGES.md §2.F2 instead.
**Status:** F2 information-rate side **CLOSED with closed-form refinement of E1.5.**
F2 pseudorandomness-battery side remains open.

---

## 1. Why this session

The /run.sh schedule put me on Task #4 (Conditional Algorithms), but
FOCUS_QUEUE.md explicitly notes "All 4 original tasks below are
COMPLETED. Do NOT re-run them." Per CLAUDE.md, when no concrete
experiment can be formulated for the rotated task, sessions should pick
a novelty target from NOVELTY_CHALLENGES.md.

I picked F2 (mod 2^k saturation) because:

1. Time-budget fit: a single-session empirical measurement.
2. Concrete falsification: numerical curve, predicted closed form.
3. Direct test of a misleading line in EDGES.md — the "0.537 bits
   constant in m" of E1.5 was tested only at x ~ 10^4 and the universal
   constant claim was unverified at other scales. F2 implicitly asks
   what the constant actually is.

## 2. What I produced

### 2.1 Code

`experiments/information_theory/pi_mod_2k_saturation/pi_mod_2k_saturation.py`

Sieve up to X in {10^3, 10^4, 10^5, 10^6, 10^7}, compute pi(x), then for
each modulus m in `2^k for k=1..10` plus cross-check {3,5,7,11,13,30,210}
compute the joint statistic
`P(pi(x) mod m = y, pi(x-1) mod m = s)` over x in [1, X] and derive
`H(Y|X)`. Compare to `h_2(pi(X)/X)`. Total wall-clock 8 s.

### 2.2 Closed-form refinement of E1.5 (the novel artifact)

For every modulus m and every X with `m <= pi(X)/100`,
```
   H( pi(x) mod m | pi(x-1) mod m, x in [1, X] )
       =  h_2( pi(X)/X )  +  O( 1 / pi(X) )
```
where `h_2(p) = -p log_2 p - (1-p) log_2 (1-p)`.

Verified to 7-decimal-place precision at X = 10^7 across moduli 2..1024.

### 2.3 Three structural facts beyond E1.5

1. The "0.537 bits" of E1.5 is **not universal** — it is the X = 10^4
   value of `h_2(pi(X)/X)`. The constant decreases with X.
2. The closed form holds in regime `m << pi(X)`. Outside this regime,
   the conditional entropy collapses (state pi(x-1) mod m starts
   encoding x).
3. Asymptotically, per-step rate vanishes: `h_2(pi(X)/X) -> 0` as
   `X -> infinity`. The CRT-cannot-win argument of E1.5 is
   asymptotically *stronger* than its X = 10^4 reading suggested.

### 2.4 Falsification

Pre-stated falsification clauses (`pi_mod_2k_saturation_results.md` §5):

* F2-A (closed-form prediction in regime `m <= pi(X)/100`): **PASS**
  at all X >= 10^4.
* F2-B (m-saturation in same regime): **PASS** in strict regime; FAIL at
  the boundary m ~ pi(X), as predicted.
* F2-C (linear-in-k scaling, the literal F2 framing): **FALSIFIED.**
  No linear-in-k component above 10^{-4} at X = 10^7.

## 3. What edges this session composes / cites

* **E1.5** — sharpened. Empirical "0.537 bits" replaced with closed form
  `h_2(pi(X)/X)`; saturation regime `m << pi(X)` made explicit.
* **T6** (top-of-EDGES.md theme) — corresponding update made.

## 4. Self-evaluation (CLAUDE.md mandatory questions)

1. **What did I produce that was not in the project before this session?**
   A closed-form derivation of E1.5's empirical "0.537 bits" constant,
   plus a regime characterization. Verified at 5 scales of X across 17
   moduli, falsifiable, all artifacts in
   `experiments/information_theory/pi_mod_2k_saturation/`.
2. **What edges did my work compose or cite?** E1.5 (sharpened); T6
   (top-of-file theme line, updated to reflect the closed form).
3. **If my session produced only duplicate closures, why?**
   Not applicable — produced a novel sharpening, not a duplicate.
4. **Next-action for the next agent:** F2 information-rate side is
   closed. F2 pseudorandomness-battery side (35 measures applied to
   pi(x) mod 4 / 8 / 16 / etc.) remains open if anyone wants to run it.
   Suggested next target: F4 (pi_BD = pi - pi_smooth pseudorandomness)
   or F6 (pi(2^k) parametric family). Both are similarly tractable
   single-session measurements.

## 5. Discipline notes

* Did NOT touch run.sh or FOCUS_QUEUE.md.
* Did NOT create _v2 / _quick / _small variants — single script with
  CLI flag `--quick` for the small run.
* Companion `_results.md` written alongside `.py`.
* No __pycache__ created (verified).
* Updated EDGES.md (E1.5 + T6), NOVELTY_CHALLENGES.md (F2 marked
  partially closed), RESEARCH_AGENDA.md (Arc 5 milestone for F2 ticked).
* Did NOT add to CLOSED_PATHS.md — this is a refinement, not a closure
  of an attack route.
* Did NOT add to novel/ — the result is a sharpening of E1.5, captured
  in EDGES.md and the experiments/ results.md. Per CLAUDE.md
  "Re-deriving E2.1 in a different basis goes in the same edge."
  This is in the same spirit (extension of one edge, not a separate
  paper-grade discovery).

## 6. Open questions raised by this work

These are not session targets — just observations worth recording:

1. The error term `O(1/pi(X))` is empirical; a rigorous bound requires
   a quantitative version of the conditional independence assumption.
   This is a tractable problem in analytic number theory: bound
   `|P[x prime | pi(x-1) mod m = s] - pi(X)/X|` using sieve methods.
2. The boundary `m ~ pi(X)` admits an exact characterization: the
   conditional entropy at m = pi(X) should be exactly `h_2(pi(X)/X) -
   delta(X)` where delta tracks the state-coverage gap. Computing the
   asymptotic of delta is interesting but not obviously useful.
3. Lean formalisation candidate: the closed-form statement is short
   enough to add to the §3 queue as a sibling of L2 (E1.5).
