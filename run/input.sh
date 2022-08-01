#!/bin/sh
# generate points in a unit square

sinsq() {
    gawk 'BEGIN {
	srand('$2');
	n = '$1';
	pi = atan2(0, -1);
	for (i = 0; i < n; ++i) {
		for (j = 0; j < n; ++j) {
			printf "%f ", -2.*pi*pi*sin(i/(n-1)*pi)*sin(j/(n-1)*pi);
		}
		printf "\n";
	}
}' > $3
}

SEED=5207364
sinsq 64 ${SEED} sinsq64.dat
sinsq 128 ${SEED} sinsq128.dat
sinsq 192 ${SEED} sinsq192.dat
sinsq 256 ${SEED} sinsq256.dat
sinsq 320 ${SEED} sinsq320.dat
sinsq 384 ${SEED} sinsq384.dat
sinsq 448 ${SEED} sinsq448.dat
sinsq 512 ${SEED} sinsq512.dat

N=100
gawk 'BEGIN {
	srand('${SEED}');
	n = '${N}';
	pi = atan2(0, -1);
	for (i = 0; i < n; ++i) {
		for (j = 0; j < n; ++j) {
			printf "%f ", -2.*pi*pi*sin(i/n*pi)*sin(j/n*pi) + 4.*cos(8.*i/n*pi) + 4.*cos(8.*j/n*pi);
		}
		printf "\n";
	}
}' > sinedge100.dat

