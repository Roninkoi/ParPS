#!/bin/sh
# determine optimal gamma for different system sizes

n=16
echo "n: 64, np: 16"
./findgamma.sh sinsq64.dat 64 ${n} 1.5 1.99 50
cp gamma.log gamma64.log

echo "n: 128, np: 16"
./findgamma.sh sinsq128.dat 128 ${n} 1.5 1.99 50
cp gamma.log gamma128.log

echo "n: 192, np: 16"
./findgamma.sh sinsq192.dat 192 ${n} 1.5 1.99 50
cp gamma.log gamma192.log

echo "n: 256, np: 16"
./findgamma.sh sinsq256.dat 256 ${n} 1.5 1.99 50
cp gamma.log gamma256.log

echo "n: 320, np: 16"
./findgamma.sh sinsq320.dat 320 ${n} 1.5 1.99 50
cp gamma.log gamma320.log

echo "n: 384, np: 16"
./findgamma.sh sinsq384.dat 384 ${n} 1.5 1.99 50
cp gamma.log gamma384.log

echo "n: 448, np: 16"
./findgamma.sh sinsq448.dat 448 ${n} 1.5 1.99 50
cp gamma.log gamma448.log

echo "n: 512, np: 16"
./findgamma.sh sinsq512.dat 512 ${n} 1.5 1.99 50
cp gamma.log gamma512.log

rm gamma.log

