/* 
 * Parallel MPI Poisson's equation solver
 * Roni Koitermaa 2022
 */

#include "poisson.h"

// Gauss-Seidel iteration
double iterateGS(double *f, double *g, double gamma, int i, int j, int n, int m)
{
	double fn = 0.;
	 // new iteration values
	if (i > 0)
		fn += f[(i - 1) * m + j];
	if (j > 0)
		fn += f[i * m + j - 1];
	// old values
	if (i < n - 1)
		fn += f[(i + 1) * m + j];
	if (j < m - 1)
		fn += f[i * m + j + 1];
	
	double r = (1. - gamma) * f[i * m + j];
	r += gamma / 4. * fn; // neighbour contribution
	r -= gamma / 4. * g[i * m + j]; // RHS

	return r;
}

double solveSOR(PoissonSolver *p)
{
	double res = 0.; // residual
	
	for (int i = p->bcu; i < p->n - p->bcd; ++i) {
		for (int j = p->bcl; j < p->m - p->bcr; ++j) {
			double fold = p->f[i * p->m + j];

			p->f[i * p->m + j] = iterateGS(p->f, p->g, p->gamma, i, j, p->n, p->m);
			
			res += fabs(p->f[i * p->m + j] - fold); // sum differences
		}
	}
	++p->it;

	return res;
}

/*
 * Black sublattice:
 * O#O#O#O
 * #O#O#O#
 * O#O#O#O
 * #O#O#O#
 */
double solveSORB(PoissonSolver *p)
{
	double res = 0.;
	
	for (int i = p->bcu; i < p->n - p->bcd; i += 1) {
		for (int j = p->bcl + i % 2; j < p->m - p->bcr; j += 2) {
			double fold = p->f[i * p->m + j];

			p->f[i * p->m + j] = iterateGS(p->f, p->g, p->gamma, i, j, p->n, p->m);
			
			res += fabs(p->f[i * p->m + j] - fold);
		}
	}
	++p->it;

	return res;
}

/*
 * Red sublattice:
 * #O#O#O#
 * O#O#O#O
 * #O#O#O#
 * O#O#O#O
 */
double solveSORR(PoissonSolver *p)
{
	double res = 0.;
	
	for (int i = p->bcu; i < p->n - p->bcd; i += 1) {
		for (int j = p->bcl + (i + 1) % 2; j < p->m - p->bcr; j += 2) {
			double fold = p->f[i * p->m + j];

			p->f[i * p->m + j] = iterateGS(p->f, p->g, p->gamma, i, j, p->n, p->m);
			
			res += fabs(p->f[i * p->m + j] - fold);
		}
	}
	++p->it;

	return res;
}

