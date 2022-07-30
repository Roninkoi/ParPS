#!/bin/sh
# determine optimal over relaxation parameter (gamma) in range [gamma0, gamma1]
# for a matrix of certain size with np processes

# run commands depending on environment
run() {
    #srun --mpi=pmix -n $1 ../poisson $2 $3 $4 $5 $6 > out.log
    mpirun --oversubscribe -np $1 ../poisson $2 $3 $4 $5 $6 > out.log
}

if [ $# != 6 ]; then
    echo "Usage: $0 <input file> <size> <np> <gamma0> <gamma1> <steps>"
    exit
fi
infile=$1
size=$2
np=$3
gamma0=$4
gamma1=$5
i=0
n=$6

mingamma=${gamma0}
maxit=1000000
minit=${maxit}

echo -n "" > gamma.log

while [ $i -lt $n ]; do
    gamma=$(echo ${gamma0} ${gamma1} $i $n | awk '{print $1+($2-$1)*$3/($4-1)}')
    run ${np} ${infile} out.dat ${size} ${gamma} 0.0001
    it=$(grep converged out.log | awk '{print $3}')
    if [ ${it} ]; then
	  echo "gamma: ${gamma} iterations: ${it}"
	  echo "gamma, iterations: ${gamma} ${it}" >> gamma.log
    else
	  echo "gamma: ${gamma} diverged"
	  it=${maxit}
    fi
    
    if [ ${it} -lt ${minit} ]; then
	  minit=${it}
	  mingamma=${gamma}
	  cp out.dat mingamma.dat
	  cp out.log mingamma.log
    fi
    
    i=$((i+1))
done

echo "minimum: ${mingamma} ${minit}"
echo "minimum: ${mingamma} ${minit}" >> gamma.log

#gnuplot -p -e "set xlabel 'Gamma'; set ylabel 'Iterations'; plot 'gamma.log' using 1:2 with lines"

