#!/bin/sh

cd run

# generate input files
./input.sh

# test cases
../parps sinsq64.dat serial.dat 64 1.5 0.001
mpirun -np 4 ../parps sinsq64.dat parallel.dat 64 1.5 0.001

# find gamma values
./rungamma.sh
# run with different numbers of processors
./run.sh

