#!/bin/sh
# run with different processor numbers and system sizes

# run command depending on environment
run() {
    #srun --mpi=pmix -n $1 ../poisson $2 $3 $4 $5 $6 > sinsq$4n$1.log
    mpirun --oversubscribe -np $1 ../poisson $2 $3 $4 $5 $6 > sinsq$4n$1.log
}

for np in 1 4 6 9 12 16 20 25 30 36 42 49 56 64; do
    
    echo "n: 64, np: ${np}"
    run ${np} sinsq64.dat sinsq64n${np}.dat 64 1.88 0.001
    
    echo "n: 128, np: ${np}"
    run ${np} sinsq128.dat sinsq128n${np}.dat 128 1.93 0.001
    
    echo "n: 192, np: ${np}"
    run ${np} sinsq192.dat sinsq192n${np}.dat 192 1.95 0.001
    
    echo "n: 256, np: ${np}"
    run ${np} sinsq256.dat sinsq256n${np}.dat 256 1.96 0.001
    
    echo "n: 320, np: ${np}"
    run ${np} sinsq320.dat sinsq320n${np}.dat 320 1.97 0.001
    
    echo "n: 384, np: ${np}"
    run ${np} sinsq384.dat sinsq384n${np}.dat 384 1.97 0.001
    
    echo "n: 448, np: ${np}"
    run ${np} sinsq448.dat sinsq448n${np}.dat 448 1.97 0.001
    
    echo "n: 512, np: ${np}"
    run ${np} sinsq512.dat sinsq512n${np}.dat 512 1.98 0.001
    
    i=$((i+1))
done

