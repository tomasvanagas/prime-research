"""Generate zeta zeros 1001..2000 (incremental on top of zeta_zeros_1000.txt)
and write zeta_zeros_2000.txt with all 2000 imaginary parts at 30 digits."""
import mpmath
import os
import time

mpmath.mp.dps = 30

src = os.path.join(os.path.dirname(__file__), "zeta_zeros_1000.txt")
out = os.path.join(os.path.dirname(__file__), "zeta_zeros_2000.txt")
chk = os.path.join(os.path.dirname(__file__), "zeta_zeros_2000.checkpoint")

with open(src) as f:
    base = [line.strip() for line in f if line.strip()]
assert len(base) == 1000

# resume from checkpoint if present
new_zeros = []
if os.path.exists(chk):
    with open(chk) as f:
        new_zeros = [line.strip() for line in f if line.strip()]
    print(f"resume from checkpoint with {len(new_zeros)} new zeros")

start = 1001 + len(new_zeros)
t0 = time.time()
for k in range(start, 2001):
    z = mpmath.im(mpmath.zetazero(k))
    s = mpmath.nstr(z, 30, strip_zeros=False)
    # ensure full precision string
    s = mpmath.nstr(z, 30)
    new_zeros.append(s)
    if k % 25 == 0:
        with open(chk, "w") as f:
            f.write("\n".join(new_zeros) + "\n")
        elapsed = time.time() - t0
        per = elapsed / max(1, k - start + 1)
        remaining = (2000 - k) * per
        print(f"k={k}  elapsed={elapsed:.1f}s  per={per:.2f}s  ETA={remaining:.0f}s", flush=True)

with open(out, "w") as f:
    f.write("\n".join(base + new_zeros) + "\n")
if os.path.exists(chk):
    os.remove(chk)
print(f"wrote {out} with {len(base)+len(new_zeros)} zeros in {time.time()-t0:.0f}s")
