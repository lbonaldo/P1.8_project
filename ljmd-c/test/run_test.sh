#!/bin/bash

echo "---------------------------------------------------"
echo "--------------- RUNNING TESTS ---------------------"
echo ""

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
echo "--------------- ALL TESTS PASSED ------------------"
echo "---------------------------------------------------"

