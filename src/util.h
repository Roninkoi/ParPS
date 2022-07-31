/* 
 * Parallel MPI Poisson's equation solver
 * Roni Koitermaa
 */

#ifndef UTIL_H
#define UTIL_H

#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <unistd.h>
#include <string.h>
#include <math.h>
#include <time.h>

#define max(a, b) ((a) > (b) ? (a) : (b))

// print n x m matrix
static void matPrint(double *a, int n, int m)
{
	for (int i = 0; i < n; ++i) {
		for (int j = 0; j < m; ++j) {
			printf("%g ", a[i * m + j]);
		}
		printf("\n");
	}
}

// print n x m matrix to file fd
static void matFilePrint(FILE *fd, double *a, int n, int m)
{
	for (int i = 0; i < n; ++i) {
		for (int j = 0; j < m; ++j) {
			fprintf(fd, "%g ", a[i * m + j]);
		}
		fprintf(fd, "\n");
	}
}

// read n x m matrix from file fd
static void matFileRead(FILE *fd, double *a, int n, int m)
{
	for (int i = 0; i < n; ++i) {
		for (int j = 0; j < m; ++j) {
			fscanf(fd, "%le ", &a[i * m + j]);
		}
	}
}

// set ns x ms submatrix s into n x m matrix a at position i0, j0
static void matSetSub(double *a, int n, int m, int i0, int j0, double *s, int ns, int ms)
{
	for (int i = 0; i < ns; ++i) {
		for (int j = 0; j < ms; ++j) {
			a[(i + i0) * m + j0 + j] = s[i * ms + j];
		}
	}
}

// get ns x ms submatrix s from n x m matrix a at position i0, j0
static void matGetSub(double *a, int n, int m, int i0, int j0, double *s, int ns, int ms)
{
	for (int i = 0; i < ns; ++i) {
		for (int j = 0; j < ms; ++j) {
			s[i * ms + j] = a[(i + i0) * m + j0 + j];
		}
	}
}

// print vector of size n
static void vecPrint(double *v, int n)
{
	for (int i = 0; i < n; ++i) {
		printf("%g ", v[i]);
	}
	printf("\n");
}

#endif

