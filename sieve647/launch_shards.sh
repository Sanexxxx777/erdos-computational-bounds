#!/bin/bash
# 8 шардов, диапазон m=[1.1e9, 1e12], overlap 1e6. Каждый шард — свой checkpoint/findings.
cd "$(dirname "$0")" || exit 1
mkdir -p logs
VENV="$HOME/Projects/math-prover/.venv/bin/python"
BOUNDS=(1100000000 126100000000 251100000000 376100000000 501100000000 626100000000 751100000000 876100000000 1001100000000)
for i in 0 1 2 3 4 5 6 7; do
  lo=${BOUNDS[$i]}
  hi=${BOUNDS[$((i+1))]}
  nohup nice -n 10 "$VENV" sieve.py --lo "$lo" --hi "$hi" --shard-overlap 1000000 \
      > "logs/shard_${i}.log" 2>&1 &
  echo "shard $i: lo=$lo hi=$hi pid=$!"
done
echo "launched 8 shards"
