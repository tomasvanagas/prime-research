#!/bin/bash
# K-scan and N-scan for B2 robustness
set -e
cd "$(dirname "$0")"

echo "===== K-scan: N=10000, K_tau ∈ {4,6,8,10,12} ====="
for K in 4 6 8 10 12; do
  python3 automorphic_l_function_basis.py 10000 $K 50 200 10 10 "b2_K${K}_Z50.json" 2>&1 | tail -8
done

echo
echo "===== N-scan: K=8, N ∈ {5000, 10000, 20000} ====="
for N in 5000 10000 20000; do
  python3 automorphic_l_function_basis.py $N 8 50 200 10 10 "b2_N${N}_K8.json" 2>&1 | tail -8
done
