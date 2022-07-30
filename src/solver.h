/* 
 * Parallel MPI Poisson's equation solver
 * Roni Koitermaa
 */

#ifndef SOLVER_H
#define SOLVER_H

#include "poisson.h"

typedef struct {
	int id, nproc; // process id, number of processes
	MPI_Status status;

	int it; // iteration
	double res; // residual
	
	char *inpath; // in/out files
	char *outpath;
	FILE *infile;
	FILE *outfile;

	// for this process
	int *dims; // process dimensions
	int *coords; // process coordinates

	// arrays for all processes in communicator, indexed by id
	int *coordsi; // coordinates of each process domain
	int *coordsj;
	int *offsetsi;  // offsets of each process domain
	int *offsetsj;
	int *sizen; // sizes of each process domain
	int *sizem;
	int *sizens;
	int *sizems;
	
	MPI_Comm cart_comm; // cartesian communicator
	int cart_id;
	int cart_root; // cartesian id of root process

	int left_id; // neighbour process ids
	int right_id;
	int up_id;
	int down_id;
	
	MPI_Datatype row_t; // datatype for matrix row
	MPI_Datatype col_t; // datatype for matrix column
	MPI_Datatype mat_t; // datatype for matrix

	int n; // number of rows
	int m; // number of columns
	double crit; // convergence criterion

	PoissonSolver psolver; // Poisson's equation solver for this process
} Solver;

void solve(Solver *s);

// initialize solver with arguments, read input
int solverInit(Solver *s, int argc, char *argv[]);
// free solver, close output file
void solverDestroy(Solver *s);

#endif

