#!/bin/bash

ORDER=4
EXE=brusselator_gsl
export OMP_THREADS=${ORDER}

# convergence study
STEPS=100
echo "running order ${ORDER} ridc with nt = ${STEPS}" >> ${EXE}.log
./${EXE} ${ORDER} ${STEPS} > n1.dat
STEPS=200
echo "running order ${ORDER} ridc with nt = ${STEPS}" >> ${EXE}.log
./${EXE} ${ORDER} ${STEPS} > n2.dat
STEPS=400
echo "running order ${ORDER} ridc with nt = ${STEPS}" >> ${EXE}.log
./${EXE} ${ORDER} ${STEPS} > n3.dat
STEPS=800
echo "running order ${ORDER} ridc with nt = ${STEPS}" >> ${EXE}.log
./${EXE} ${ORDER} ${STEPS} > n4.dat
STEPS=1600
echo "running order ${ORDER} ridc with nt = ${STEPS}" >> ${EXE}.log
./${EXE} ${ORDER} ${STEPS} > n5.dat


rm -f diff.log

# compare with reference solution
diff n1.dat ref/n1.dat >> diff.log
diff n2.dat ref/n2.dat >> diff.log
diff n3.dat ref/n3.dat >> diff.log
diff n4.dat ref/n4.dat >> diff.log
diff n5.dat ref/n5.dat >> diff.log


if [ -s diff.log ]
then
	echo ":test-result: FAIL" >> ${EXE}.trs
else
	echo ":test-result: PASS" >> ${EXE}.trs
fi
