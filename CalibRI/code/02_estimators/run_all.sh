#!/bin/bash
# Full pipeline for the sensitivity analyses (repeated splits; RINet with 1,000 bootstrap resamples). Resumable.
M=code/02_estimators/rinet_models; C=code/02_estimators; O=work/res_b4
mkdir -p $O
python3 $C/rinet_boot.py work/ribundle4 $M $O/rinet_b4.csv 1000 all
python3 $C/rinet_boot.py work/ribundle  $M $O/rinet_sim.csv 1000 sim_small
python3 $C/rinet_boot.py work/ribundle3 $M $O/rinet_b3.csv 1000 all
echo ALL_RINET_DONE
