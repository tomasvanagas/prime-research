I'll proceed with this proposal session. Let me load the tools I'll need and start brainstorming fresh approaches.
Got enough fresh context. Let me draft the proposals and testable experiments.
Now I have a good sense of what's been tried. Let me design proposals that genuinely differ from the existing ones. I'll save the proposals doc first, then build runnable tests.
Now let me build the three test scripts. I'll start with the MPS / TT one.
Now let me run all three. First a quick smoke test to make sure they don't error before committing to a full run.
Need to speed up `riemann_R_inverse` — at 0.1s/call, N=8192 takes ~15 min. Let me precompute and cache.
Cached. Delta std=42, range [-176, 147] for N=8192. Let me update the three scripts to use the cache, then run them.
Same surgery for the other two scripts:
Now run all three in sequence.
Morton ordering for L=13 (odd) overflows. Let me fix.
Volume-law confirmed for every ordering. Proposal A closed. Now Proposal B.
LASSO collapsed to zero (alpha too aggressive); OLS overfits at K≥256. Let me refine the alpha grid and add a fairer per-K test.
LASSO confirms: zero-mode basis is not sparse for delta. Naïve "always predict 0" beats every sparse fit. Proposal B closed. Now C:
All three closed empirically. Let me verify each `_results.md` exists and has good content, plus check the precompute script has its companion results file.
`session49_precompute_data.py` is missing a `_results.md` companion. Let me write one + the missing-results scan.
Let me also strengthen the auto-written results files with more detail, since the current versions are short.
Now the missing-results scan and pycache cleanup.
All `_results.md` companions present, pycache cleared. Final summary write-up + `.run_state`.
Session 49 done. Three fresh proposals tested, all closed by negative empirical results.

**Proposals attacked:**
- **A. Reordered TT compression of delta(n)** — tested 7 orderings (identity, bit-reversal, Gray, Morton, 2-adic, sort-by-R^{-1}, random). All saturate the volume-law bond ceiling 64 = √8192/2, matching gaussian noise. No ordering helps.
- **B. LASSO compressed-sensing in zero-mode basis** — across alpha sweep {1e-6,…,1e-1} and K up to 1024 zeros (γ ≤ 1447), LASSO selects nnz=0 as best test-RMSE solution at every K. OLS overfits catastrophically at K≥256.
- **C. Learned residues at fixed zeros** — ridge/OLS at K up to 2000 zeros never exceeds the 1.1% naive R^{-1} round-recovery rate. Learned coefficients fail to absorb the truncation tail.

All three close by **failure mode I (information loss)**: delta(n) is empirically indistinguishable from GUE-random in every basis tested.

**Files written:**
- `archive/ephemeral/proposals_session.md` — three fresh proposals
- `experiments/proposals/session49_{precompute_data,mps_delta_orderings,compressed_sensing_zeros,neural_zero_residue}.py` + matching `_results.md`
- `experiments/proposals/session49_data.npz` — cache of p(n), R^{-1}(n), delta(n) for n=1..8192
- `.run_state` set to 30.

**Sources:**
- [Cramér's conjecture / efficient algorithm conditional](https://michaelnielsen.org/polymath/index.php?title=An_efficient_algorithm_exists_if_Cramer's_conjecture_holds)
- [Maynard-Tao sieve weights (2024 / 2025)](https://arxiv.org/abs/2403.19696)
- [Counting primes in Õ(n^{2/3})](https://codeforces.com/blog/entry/91632)
- [Tensor network compression survey 2026](https://arxiv.org/abs/2505.20132)
