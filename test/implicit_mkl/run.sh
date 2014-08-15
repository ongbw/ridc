#!/bin/bash

make

ORDER=4
EXE=implicit_ode_mkl.exe
export OMP_THREADS=${ORDER}

# convergence study
STEPS=10
echo "running order ${ORDER} ridc with nt = ${STEPS}"
./${EXE} ${ORDER} ${STEPS} > n1.dat
STEPS=20
echo "running order ${ORDER} ridc with nt = ${STEPS}"
./${EXE} ${ORDER} ${STEPS} > n2.dat
STEPS=40
echo "running order ${ORDER} ridc with nt = ${STEPS}"
./${EXE} ${ORDER} ${STEPS} > n3.dat
STEPS=80
echo "running order ${ORDER} ridc with nt = ${STEPS}"
./${EXE} ${ORDER} ${STEPS} > n4.dat
STEPS=160
echo "running order ${ORDER} ridc with nt = ${STEPS}"
./${EXE} ${ORDER} ${STEPS} > n5.dat

