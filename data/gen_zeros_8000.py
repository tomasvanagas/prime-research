"""Generate zeta zeros 2001..8000 (incremental on top of zeta_zeros_2000.txt)
and write zeta_zeros_8000.txt with all 8000 imaginary parts at 15 digits.

Reduced precision (15 dps vs 30 in earlier files) — sufficient for correlation
analysis to 4 decimal places, ~3x faster than 30 dps.

Per-zero cost ~0.5s on this machine; estimated ETA ~50 minutes for 6000 zeros.
"""
import mpmath
import os
import time

mpmath.mp.dps = 15

src = os.path.join(os.path.dirname(__file__), "zeta_zeros_2000.txt")
out = os.path.join(os.path.dirname(__file__), "zeta_zeros_8000.txt")
chk = os.path.join(os.path.dirname(__file__), "zeta_zeros_8000.checkpoint")

with open(src) as f:
    base = [line.strip() for line in f if line.strip()]
assert len(base) == 2000

new_zeros = []
if os.path.exists(chk):
    with open(chk) as f:
        new_zeros = [line.strip() for line in f if line.strip()]
    print(f"resume from checkpoint with {len(new_zeros)} new zeros", flush=True)

start = 2001 + len(new_zeros)
t0 = time.time()
for k in range(start, 8001):
    z = mpmath.im(mpmath.zetazero(k))
    s = mpmath.nstr(z, 15)
    new_zeros.append(s)
    if k % 50 == 0:
        with open(chk, "w") as f:
            f.write("\n".join(new_zeros) + "\n")
        elapsed = time.time() - t0
        per = elapsed / max(1, k - start + 1)
        remaining = (8000 - k) * per
        print(f"k={k}  elapsed={elapsed:.1f}s  per={per:.2f}s  ETA={remaining:.0f}s", flush=True)

with open(out, "w") as f:
    f.write("\n".join(base + new_zeros) + "\n")
if os.path.exists(chk):
    os.remove(chk)
print(f"wrote {out} with {len(base)+len(new_zeros)} zeros in {time.time()-t0:.0f}s")
