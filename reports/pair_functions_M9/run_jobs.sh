#!/bin/bash
export THERMINATOR_SHARE=/home/alice/therminator2-hboltz/share
cd /home/alice/therminator2-hboltz/build/t2
for s in 1 2 3 4 5 6 7 8; do
  ../work/ex_pair_functions_mc.exe --harmonics --weak on --events 25000 --seed $s \
    --output m9/jobs/harm_weakon_$s.root > m9/jobs/harm_weakon_$s.log 2>&1 &
done
wait
for s in 101 102 103 104 105 106 107 108; do
  ../work/ex_pair_functions_mc.exe --harmonics --weak off --events 25000 --seed $s \
    --output m9/jobs/harm_weakoff_$s.root > m9/jobs/harm_weakoff_$s.log 2>&1 &
done
wait
hadd -f m9/harm_weakon.root m9/jobs/harm_weakon_?.root > m9/hadd_on.log 2>&1
hadd -f m9/harm_weakoff.root m9/jobs/harm_weakoff_1??.root > m9/hadd_off.log 2>&1
echo done
