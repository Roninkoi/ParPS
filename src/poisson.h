/* 
 * Parallel MPI Poisson's equation solver
 * Roni Koitermaa 2022
 */

#ifndef POISSON_H
#define POISSON_H

#include "util.h"

typedef struct {
	double *f; // solution
	double *fs;
	double *g; // RHS
	double *gs;
	
	int it; // iteration
	
	int n; // number of rows
	int ns;
	int m; // number of columns
	int ms;
	
	// Dirichlet BCs
	int bcr; // right
	int bcl; // left
	int bcu; // up
	int bcd; // down
	
	double gamma; // over relaxation parameter
} PoissonSolver;

// solve using successive over-relaxation (serial)
double solveSOR(PoissonSolver *p);
// solve black sublattice
double solveSORB(PoissonSolver *p);
// solve red sublattice
double solveSORR(PoissonSolver *p);

#endif

