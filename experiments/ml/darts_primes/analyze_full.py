"""Combined analysis across PRIMES, control, and calibrated-1-bit conditions."""
import json
import numpy as np
from scipy import stats

R = json.load(open("/apps/aplikacijos/prime-research/experiments/ml/darts_primes/run/results.json"))
C = json.load(open("/apps/aplikacijos/prime-research/experiments/ml/darts_primes/run/calibrated_results.json"))


def grab(d, k):
    return np.array([s["final_loss"] for s in d[k]])


P = grab(R, "primes")
RND = grab(R, "controls")
CAL = grab(C, "calibrated")

print("=== Final-loss summary (BCE, nats) ===")
for label, arr in [("PRIMES (chi_P)", P),
                   ("matched-density random", RND),
                   ("calibrated-1-bit random", CAL)]:
    print(f"  {label:32s}: mean={arr.mean():.4f}  std={arr.std(ddof=1):.4f}  "
          f"  min={arr.min():.4f}  max={arr.max():.4f}")

# Reference baselines
print()
print("Reference predictors at N=12, density 0.1377:")
for nm, val in [("constant predictor (entropy)", 0.4008),
                ("1-bit oddness", 0.2961),
                ("mod-6 wheel", 0.2298),
                ("mod-30 wheel", 0.1905),
                ("mod-210 wheel", 0.1614)]:
    print(f"  {nm:32s}: {val:.4f}")

print()
print("=== Pairwise Welch t-tests on final loss ===")
for (a_lab, a), (b_lab, b) in [
    (("PRIMES", P), ("RANDOM_DENS", RND)),
    (("PRIMES", P), ("CALIBRATED",  CAL)),
    (("RANDOM_DENS", RND), ("CALIBRATED", CAL)),
]:
    t, p = stats.ttest_ind(a, b, equal_var=False)
    print(f"  {a_lab:12s} vs {b_lab:12s}: t={t:+.3f}  p={p:.4f}")

print()
print("=== Mean gaps (in nats) ===")
print(f"  PRIMES  − RANDOM_DENS:  {P.mean() - RND.mean():+.4f}  (oddness reduction predicts {0.2961 - 0.4008:+.4f})")
print(f"  PRIMES  − CALIBRATED :  {P.mean() - CAL.mean():+.4f}  (zero if all PRIMES advantage is oddness)")
print(f"  RANDOM_DENS  − CALIBRATED:  {RND.mean() - CAL.mean():+.4f}")

print()
print("=== Conclusion ===")
oddness_gap = abs(P.mean() - 0.2961)
calib_gap = abs(P.mean() - CAL.mean())
print(f"|PRIMES − oddness baseline| = {oddness_gap:.4f}  (DARTS is at parity floor: {oddness_gap < 0.005})")
print(f"|PRIMES − calibrated|       = {calib_gap:.4f}  (oddness fully explains gap: {calib_gap < 0.005})")
