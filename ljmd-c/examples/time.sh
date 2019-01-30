#!/bin/bash

#/usr/bin/time ../ljmd-split.x < argon_2916.inp

for N in 2916
do
    /usr/bin/time ../ljmd-split.x < argon_$N.inp 2> tmp.dat
    cat tmp.dat | head -n 1 | cut -d ' ' -f1 1> tmp2.dat
    awk -v N=$N '{print N " Step 1-cell: " $1}' tmp2.dat 1>> ../test/time.txt

rm tmp.dat
rm tmp2.dat
done
