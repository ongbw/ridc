#!/bin/bash

# build brusselator.exe
make

ORDER=3
EXE=reaction_diffusion.exe
export OMP_THREADS=${ORDER}

# convergence study
STEPS=10000
echo "running order ${ORDER} ridc with nt = ${STEPS}"
./${EXE} ${ORDER} ${STEPS} > n1.dat
STEPS=20000
echo "running order ${ORDER} ridc with nt = ${STEPS}"
./${EXE} ${ORDER} ${STEPS} > n2.dat
STEPS=30000
echo "running order ${ORDER} ridc with nt = ${STEPS}"
./${EXE} ${ORDER} ${STEPS} > n3.dat
STEPS=40000
echo "running order ${ORDER} ridc with nt = ${STEPS}"
./${EXE} ${ORDER} ${STEPS} > n4.dat
STEPS=80000
echo "running order ${ORDER} ridc with nt = ${STEPS}"
./${EXE} ${ORDER} ${STEPS} > n5.dat

