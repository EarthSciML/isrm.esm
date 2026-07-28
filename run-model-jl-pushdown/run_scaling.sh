#!/bin/bash
# Sequential ladder, one job at a time. Guard is SIZE-AWARE: a ppl=56 run peaks
# ~1.6 GB, a ppl=441 run much more, so each rung demands its own headroom.
cd /projects/illinois/eng/cee/ctessum/ctessum/code/isrm.esm/run-model-jl-pushdown
run_one() {  # $1=firstn  $2=required_free_gb  $3=heap_hint
  for attempt in 1 2 3 4 5 6 7 8 9 10 11 12; do
    tot=$(ps -eo rss,comm | grep -i julia | awk '{s+=$1} END{printf "%.0f", s/1048576}')
    free=$((40 - ${tot:-40}))
    if [ "$free" -ge "$2" ]; then
      echo "########## FIRSTN=$1  (box ${tot}GB, ${free}GB free, need $2GB)"
      FIRSTN=$1 L3_SR_DIR=scaling_sr_$1 timeout 3600 julia -t 2 --heap-size-hint=$3 \
        --project=. mem_scaling.jl 2>&1 | grep -E "MACHINE-READABLE|took|fastpath|peak RSS DURING"
      echo; return 0
    fi
    echo "### FIRSTN=$1 waiting: ${free}GB free, need $2GB (attempt $attempt)"
    sleep 300
  done
  echo "### FIRSTN=$1 GAVE UP — box stayed busy"
}
run_one 1000  6  4G
run_one 2500  8  6G
run_one 5000 10  8G
run_one 10000 14 12G
echo "LADDER COMPLETE"
