"""Combine the partial sequential checkpoint and the parallel chunk files
into a single zeta_zeros_8000.txt file with zeros 1..8000.
"""
import os, sys, glob

ROOT = os.path.dirname(__file__)

# 1. Existing 2000 zeros
with open(os.path.join(ROOT, 'zeta_zeros_2000.txt')) as f:
    base = [line.strip() for line in f if line.strip()]
assert len(base) == 2000
print(f"Loaded {len(base)} zeros from zeta_zeros_2000.txt")

# 2. Sequential checkpoint (covers 2001..end_seq)
seq_chk = os.path.join(ROOT, 'zeta_zeros_8000.checkpoint')
seq = []
if os.path.exists(seq_chk):
    with open(seq_chk) as f:
        seq = [line.strip() for line in f if line.strip()]
    print(f"Loaded {len(seq)} zeros from sequential checkpoint (covers 2001..{2000+len(seq)})")
else:
    seq = []

# 3. Parallel chunks (each line is "k value")
chunks = sorted(glob.glob(os.path.join(ROOT, 'parallel_chunks', 'chunk_*.txt')))
parallel_dict = {}
for chunk_path in chunks:
    if chunk_path.endswith('.tmp') or chunk_path.endswith('.log'):
        continue
    with open(chunk_path) as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) == 2:
                k = int(parts[0])
                val = parts[1]
                parallel_dict[k] = val
print(f"Loaded {len(parallel_dict)} zeros from parallel chunks")

# 4. Build the full list 1..8000
combined = list(base)  # zeros 1..2000
seq_end = 2000 + len(seq)
combined.extend(seq)  # zeros 2001..seq_end
# zeros (seq_end+1)..8000 should come from parallel_dict
missing = []
for k in range(seq_end + 1, 8001):
    if k in parallel_dict:
        combined.append(parallel_dict[k])
    else:
        missing.append(k)

if missing:
    print(f"WARNING: missing {len(missing)} zeros: {missing[:10]}...{missing[-10:]}")
    # Truncate at first missing
    combined = combined[:seq_end + sum(1 for k in range(seq_end + 1, missing[0]) if k in parallel_dict)]
    print(f"Truncated to {len(combined)} zeros (last index: {len(combined)})")

out_path = os.path.join(ROOT, 'zeta_zeros_8000.txt')
with open(out_path, 'w') as f:
    f.write("\n".join(combined) + "\n")
print(f"WROTE {out_path} with {len(combined)} zeros")
