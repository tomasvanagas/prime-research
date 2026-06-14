#!/bin/bash
# Detached reach launcher (S529). Runs the large-x verification reach for
# n=25 then n=26 SEQUENTIALLY, each appending to its own log with a DONE
# marker. Launched with `setsid ... </dev/null &` so it is reparented to
# PID 1 and SURVIVES the research-cycle boundary -- the failure mode that
# killed the S528 reach (harness-tracked background processes are terminated
# when the launching cycle's claude session ends). Ground truth (sympy):
#   pi(2^25-1)=2063689   pi(2^26-1)=3957809
cd "$(dirname "$0")" || exit 1
PY=python3
stamp() { date -u +%Y-%m-%dT%H:%M:%SZ; }

{
  echo "######## n=25 reach (honest + delta_pi cheat) START $(stamp) ########"
} > run_n25.log
$PY large_x_benchmark.py --n 25 --seed 1 >> run_n25.log 2>&1
echo "######## n=25 reach DONE $(stamp) ########" >> run_n25.log

{
  echo "######## n=26 reach (honest, --no-cheat) START $(stamp) ########"
} > run_n26.log
$PY large_x_benchmark.py --n 26 --seed 1 --no-cheat >> run_n26.log 2>&1
echo "######## n=26 reach DONE $(stamp) ########" >> run_n26.log
