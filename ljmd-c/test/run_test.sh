#!/bin/bash

echo "---------------------------------------------------"
echo "--------------- RUNNING TESTS ---------------------"
echo ""

echo -n "Test force.o ... "
./check_force.x < argon_108.inp
if [[ ! $? -eq 0  ]]; then
	echo " Wrong answer"
	exit 1
fi
echo " Accepted."

echo -n "Test overall ... "
../ljmd-split.x < argon_108.inp >/dev/null
ret=$?
./check_overall_energy.x < argon_108.inp
if [[ (! $? -eq 0 ) || ( ! $ret -eq 0 ) ]]; then
	echo " Wrong answer"
	exit 1
fi
echo " Accepted."

echo -n "Test something else ... "
./check_another_thing.x
if [[ ! $? -eq 0  ]]; then
	echo " Wrong answer"
	exit 1
fi
echo " Accepted."

echo ""
export A=$( { /usr/bin/time -f "%e" echo "1" ; } 2>&1 >/dev/null ) 
t1=$( { /usr/bin/time -f "%e" ../ljmd-serial.x < ../examples/argon_108.inp ; } 2>&1 >/dev/null ) 
t2=$( { /usr/bin/time -f "%e" ../ljmd-split.x < ../examples/argon_108.inp ; } 2>&1 > /dev/null)
sup=$(python -c "print('{0:.2f}'.format($t1/$t2))")

echo "Time for serial: $t1 seconds"
echo "Time for split: $t2 seconds"
echo "Speed up x$sup"

echo ""
echo "--------------- ALL TESTS PASSED ------------------"
echo "---------------------------------------------------"

