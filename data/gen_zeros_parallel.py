"""Parallel-worker driver: generate zeta zeros for a given (start, end) range
at 15 dps and write to a chunk file. Designed to be spawned in multiple
processes for disjoint ranges so the main gen_zeros_8000.py process and
these workers cover the (2001..8000) range together.

Usage: python3 gen_zeros_parallel.py <start_k> <end_k> <out_path>
"""
import sys
import time
import mpmath
import os

mpmath.mp.dps = 15

start = int(sys.argv[1])
end = int(sys.argv[2])  # inclusive
out_path = sys.argv[3]

results = []
t0 = time.time()
for k in range(start, end + 1):
    z = mpmath.im(mpmath.zetazero(k))
    s = mpmath.nstr(z, 15)
    results.append(s)
    if (k - start + 1) % 50 == 0:
        elapsed = time.time() - t0
        per = elapsed / (k - start + 1)
        remaining = (end - k) * per
        print(f"[{out_path}] k={k}  done={k-start+1}/{end-start+1} "
              f"elapsed={elapsed:.0f}s ETA={remaining:.0f}s", flush=True)
        # Incremental checkpoint
        with open(out_path + '.tmp', 'w') as f:
            f.write("\n".join(f"{start+i} {z}" for i, z in enumerate(results)) + "\n")

with open(out_path, 'w') as f:
    f.write("\n".join(f"{start+i} {z}" for i, z in enumerate(results)) + "\n")
if os.path.exists(out_path + '.tmp'):
    os.remove(out_path + '.tmp')
print(f"[{out_path}] DONE in {time.time()-t0:.0f}s ({len(results)} zeros)", flush=True)
