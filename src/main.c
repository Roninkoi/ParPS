/* 
 * Parallel MPI Poisson's equation solver
 * Roni Koitermaa 2022
 */

#include "solver.h"

int main(int argc, char *argv[])
{
	Solver solver;

	int status = MPI_Init(&argc, &argv);
	if (status != MPI_SUCCESS) {
		printf("MPI initialization failed\n");
		return status;
	}

	status = solverInit(&solver, argc, argv);
	if (status != 0) {
		printf("Solver initialization failed\n");
		return status;
	}

	solve(&solver);

	solverDestroy(&solver);
	
	MPI_Finalize();
	
	return 0;
}

