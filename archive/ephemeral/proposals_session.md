# Session 35: Circuit Complexity Deep-Dive

## Date: 2026-04-05

## Focus
Comprehensive circuit complexity analysis of pi(x) mod 2, exploring the last recommended
direction (MKtP/Brandt framework) and measuring new structural properties.

## Experiments Completed

### 1. Approximate Circuit Complexity — Phase Transition Search
**File:** `experiments/circuit_complexity/approx_circuit_complexity.py`
**Result:** NO phase transition. Rank-k approximation accuracy increases gradually.
Errors are spatially uniform. No "easy core" exists.

### 2. Meta-Complexity / MKtP Framework Analysis
**File:** `experiments/circuit_complexity/meta_complexity_analysis.py`
**Result:** CLOSED. Framework is a reformulation, not a new tool. Kt(T_N) = O(2^N·N)
by sieve regardless of circuit size. Brandt framework too generic.

### 3. GF(2) ANF / SLP Structure Analysis
**File:** `experiments/circuit_complexity/gf2_slp_structure.py`
**Result:** CLOSED. ANF sparsity = 0.50 (exactly random), CSE savings match random (±1%),
SLP length Θ(2^N). No GF(2) algebraic compression possible.

### 4. SAT-Based Minimum Circuit Size (pending agent completion)
**File:** `experiments/circuit_complexity/sat_circuit_minimization.py`

### 5. Threshold Gate Circuit Construction (pending agent completion)
**File:** `experiments/circuit_complexity/threshold_circuit_construction.py`

## Cross-Cutting Insight

pi(x) mod 2 is random-like in EVERY computable structural measure tested (20+).
The function is indistinguishable from random by every measure that doesn't require
solving the original problem.

## New Closed Paths: 7 (total 648+ across 35 sessions)
