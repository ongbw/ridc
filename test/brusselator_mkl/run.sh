#!/bin/bash

make

ORDER=4
EXE=brusselator.exe
export OMP_THREADS=${ORDER}

# convergence study
STEPS=100
echo "running order ${ORDER} ridc with nt = ${STEPS}"
./${EXE} ${ORDER} ${STEPS} > n1.dat
STEPS=200
echo "running order ${ORDER} ridc with nt = ${STEPS}"
./${EXE} ${ORDER} ${STEPS} > n2.dat
STEPS=400
echo "running order ${ORDER} ridc with nt = ${STEPS}"
./${EXE} ${ORDER} ${STEPS} > n3.dat
STEPS=800
echo "running order ${ORDER} ridc with nt = ${STEPS}"
./${EXE} ${ORDER} ${STEPS} > n4.dat
STEPS=1600
echo "running order ${ORDER} ridc with nt = ${STEPS}"
./${EXE} ${ORDER} ${STEPS} > n5.dat

