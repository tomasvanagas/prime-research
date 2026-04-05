# LGV DAG Construction for pi(x): Results

**Script:** `lgv_dag_construction.py`
**Session:** 14

## What Was Tested
Whether a DAG on poly(N) vertices can be constructed such that the LGV (Lindstrom-Gessel-Viennot) signed count of non-intersecting lattice paths equals pi(x). Tried binary decision DAGs, sieve-based DAGs, recursive (Meissel-Lehmer) DAGs, and hybrids.

## Key Findings
- Sieve-based DAG: has 2^{pi(sqrt(x))} nodes (exponential), one path per subset of primes
- Binary decision DAG: equivalent to a branching program; requires exponential size for random-like functions
- Recursive DAG: the Meissel-Lehmer recursion tree has O(x^{2/3}) nodes -- better than exponential but still not poly(N)
- The LGV lemma requires non-intersecting paths, which imposes additional constraints that prevent compression
- No poly(N)-vertex DAG found

## Verdict
**CLOSED**
**Failure Mode:** Information loss (all DAG constructions require super-polynomial nodes)

## One-Line Summary
No poly(N)-vertex DAG yields pi(x) via LGV non-intersecting path counts; sieve DAGs are exponential, recursive DAGs are O(x^{2/3}).
