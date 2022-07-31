/* 
 * Parallel MPI Poisson's equation solver
 * Roni Koitermaa 2022
 */

#include <signal.h>

#include "solver.h"

// do 2D domain decomposition for nproc processes
int decompose(Solver *s)
{
	double *f = s->psolver.f; // total solution
	double *g = s->psolver.g; // total RHS
	
	const int dim = 2; // 2D
	s->dims = calloc(sizeof(int), dim); // size of each dimension
	s->dims[0] = (int) sqrt((double) s->nproc); // make as square as possible
	//printf("dims: %i %i\n", s->dims[0], s->dims[1]);
	MPI_Dims_create(s->nproc, dim, s->dims); // create 2D dims[0] * dims[1] grid with nproc processes

	// set fixed boundaries
	int periods[dim];
	for (int i = 0; i < dim; ++i)
		periods[i] = 0; // not periodic

	int reorder[dim]; // allow reordering of ids
	for (int i = 0; i < dim; ++i)
		reorder[i] = 1;

	MPI_Cart_create(MPI_COMM_WORLD, dim, s->dims, periods, 1, &s->cart_comm); // create new communicator for grid

	MPI_Comm_rank(s->cart_comm, &s->cart_id); // get process id in new communicator

	s->cart_root = s->cart_id; // root process id in cartesian communicator
	MPI_Bcast(&s->cart_root, 1, MPI_INT, 0, MPI_COMM_WORLD);

	s->coords = calloc(sizeof(int), dim);
	MPI_Cart_coords(s->cart_comm, s->cart_id, dim, s->coords); // get coords of process

	/*if (s->dims[0] != s->dims[1] || s->n % s->dims[1] != 0 || s->m % s->dims[0] != 0) { // not divisible
	  fprintf(stderr, "Error! Make sure processes divide grid!\n");
	  return 11;
	  }*/
	s->psolver.n = s->n / s->dims[1];
	s->psolver.m = s->m / s->dims[0];

	/* if not perfectly divisible, put remaining to left/up processes
	 *  3   2
	 * <->  <>
	 * XXXOOXXOOXXOOXX
	 * XXXOOXXOOXXOOXX
	 * XXXOOXXOOXXOOXX
	 * <------------->
	 *        15
	 * (not divisible into 7 processes)
	 */
	if (s->n % s->dims[1] != 0 && s->coords[1] == 0) {
		s->psolver.n = s->n - s->psolver.n * (s->dims[1] - 1);
	}
	if (s->m % s->dims[0] != 0 && s->coords[0] == 0) {
		s->psolver.m = s->m - s->psolver.m * (s->dims[0] - 1);
	}
	
	MPI_Cart_shift(s->cart_comm, 0, 1, &s->left_id, &s->right_id); // get neighbours, dim 0
	MPI_Cart_shift(s->cart_comm, 1, 1, &s->up_id, &s->down_id); // dim 1

	printf("id %i (%i) neighbours: right = %i, up = %i, left = %i, down = %i\n", s->cart_id, s->id, s->right_id, s->up_id, s->left_id, s->down_id);

	if (s->cart_id == 0) {
		printf("process topology:\n");
		int neigh_coords[dim];
		for (int i = 0; i < s->dims[1]; ++i) { // print table of process ids to show xy coords
			for (int j = 0; j < s->dims[0]; ++j) {
				int neigh_id;
				neigh_coords[1] = i;
				neigh_coords[0] = j;
				MPI_Cart_rank(s->cart_comm, neigh_coords, &neigh_id);
				printf("%i ", neigh_id);
			}
			printf("\n");
		}
	}

	/*
	 * boundaries of process domains
	 * X = unused area
	 * O = boundary condition from initial g
	 * # = boundaries from other process
	 *
	 * outside boundary:
	 * X X X X X
	 * X O O O #
	 * X O . . #
	 * X O . . #
	 * X # # # X
	 *
	 * inside boundary:
	 * X # # # X
	 * # . . . #
	 * # . . . #
	 * # . . . #
	 * X # # # X
	 */
	
	// determine boundary size based on neighbouring processes
	s->psolver.bcr = s->right_id < 0 ? 2 : 1;
	s->psolver.bcl = s->left_id < 0 ? 2 : 1;
	s->psolver.bcu = s->up_id < 0 ? 2 : 1;
	s->psolver.bcd = s->down_id < 0 ? 2 : 1;

	s->psolver.ns = s->psolver.n; // actual solution size
	s->psolver.n += 2; // add boundaries
	
	s->psolver.ms = s->psolver.m; // actual solution size
	s->psolver.m += 2; // add boundaries
	
	// make data types
	MPI_Type_vector(s->psolver.ms, 1, 1, MPI_DOUBLE, &s->row_t);
	MPI_Type_commit(&s->row_t);
	MPI_Type_vector(s->psolver.ns, 1, s->psolver.m, MPI_DOUBLE, &s->col_t);
	MPI_Type_commit(&s->col_t);
	MPI_Type_vector(s->psolver.n * s->psolver.m, 1, 1, MPI_DOUBLE, &s->mat_t);
	MPI_Type_commit(&s->mat_t);

	s->psolver.f = calloc(sizeof(double), s->psolver.n * s->psolver.m);
	s->psolver.g = calloc(sizeof(double), s->psolver.n * s->psolver.m);

	/*
	 * Structure of f/g:
	 * O = fs/gs pointer
	 * # = process boundaries
	 *
	 *     ns
	 *   <--->
	 * # # # # #
	 * # O . . #
	 * # . . . #
	 * # . . . #
	 * # # # # #
	 * <------->
	 *     n
	 */
	s->psolver.fs = s->psolver.f + s->psolver.m + 1; // pointer to start of actual solution
	s->psolver.gs = s->psolver.g + s->psolver.m + 1;

	s->coordsi = calloc(sizeof(int), s->nproc);
	s->coordsj = calloc(sizeof(int), s->nproc);
	s->offsetsi = calloc(sizeof(int), s->nproc);
	s->offsetsj = calloc(sizeof(int), s->nproc);
	s->sizen = calloc(sizeof(int), s->nproc);
	s->sizem = calloc(sizeof(int), s->nproc);
	s->sizens = calloc(sizeof(int), s->nproc);
	s->sizems = calloc(sizeof(int), s->nproc);

	if (s->id == 0) {
		MPI_Gather(&s->coords[0], 1, MPI_INT, s->coordsj, 1, MPI_INT, s->cart_root, s->cart_comm);
		MPI_Gather(&s->coords[1], 1, MPI_INT, s->coordsi, 1, MPI_INT, s->cart_root, s->cart_comm);
		MPI_Gather(&s->psolver.n, 1, MPI_INT, s->sizen, 1, MPI_INT, s->cart_root, s->cart_comm);
		MPI_Gather(&s->psolver.m, 1, MPI_INT, s->sizem, 1, MPI_INT, s->cart_root, s->cart_comm);
		MPI_Gather(&s->psolver.ns, 1, MPI_INT, s->sizens, 1, MPI_INT, s->cart_root, s->cart_comm);
		MPI_Gather(&s->psolver.ms, 1, MPI_INT, s->sizems, 1, MPI_INT, s->cart_root, s->cart_comm);

		int offsi = 0;
		int offsj = 0;
		int sizen = 0;
		int sizem = 0;
		for (int i = 0; i < s->dims[1]; ++i) { // calculate offsets from coordinates and sizes
			offsj = 0;
			for (int j = 0; j < s->dims[0]; ++j) { // go through coordinates
				for (int id = 0; id < s->nproc; ++id) { // find matching process
					if (s->coordsi[id] == i && s->coordsj[id] == j) {
						s->offsetsi[id] = offsi;
						s->offsetsj[id] = offsj;
						sizen = s->sizens[id];
						sizem = s->sizems[id];
						break;
					}
				}
				offsj += sizem;
			}
			offsi += sizen;
		}

		printf("domains:\n");
		for (int i = 0; i < s->nproc; ++i)
			printf("i: %i, j: %i, n: %i, m: %i, ns: %i, ms: %i, offsi: %i, offsj: %i (%s)\n", s->coordsi[i], s->coordsj[i], s->sizen[i], s->sizem[i], s->sizens[i], s->sizems[i], s->offsetsi[i], s->offsetsj[i], (s->offsetsi[i] + s->offsetsj[i]) % 2 == 0 ? "black" : "red");
	}
	else {
		MPI_Gather(&s->coords[0], 1, MPI_INT, NULL, 1, MPI_INT, s->cart_root, s->cart_comm);
		MPI_Gather(&s->coords[1], 1, MPI_INT, NULL, 1, MPI_INT, s->cart_root, s->cart_comm);
		MPI_Gather(&s->psolver.n, 1, MPI_INT, NULL, 1, MPI_INT, s->cart_root, s->cart_comm);
		MPI_Gather(&s->psolver.m, 1, MPI_INT, NULL, 1, MPI_INT, s->cart_root, s->cart_comm);
		MPI_Gather(&s->psolver.ns, 1, MPI_INT, NULL, 1, MPI_INT, s->cart_root, s->cart_comm);
		MPI_Gather(&s->psolver.ms, 1, MPI_INT, NULL, 1, MPI_INT, s->cart_root, s->cart_comm);
	}

	// send coordinates, sizes and offsets to all processes
	MPI_Bcast(s->coordsi, s->nproc, MPI_INT, 0, s->cart_comm);
	MPI_Bcast(s->coordsj, s->nproc, MPI_INT, 0, s->cart_comm);
	MPI_Bcast(s->sizen, s->nproc, MPI_INT, 0, s->cart_comm);
	MPI_Bcast(s->sizem, s->nproc, MPI_INT, 0, s->cart_comm);
	MPI_Bcast(s->sizens, s->nproc, MPI_INT, 0, s->cart_comm);
	MPI_Bcast(s->sizems, s->nproc, MPI_INT, 0, s->cart_comm);
	MPI_Bcast(s->offsetsi, s->nproc, MPI_INT, 0, s->cart_comm);
	MPI_Bcast(s->offsetsj, s->nproc, MPI_INT, 0, s->cart_comm);
	
	// divide f and g to all processes
	if (s->id == 0) {
		double *fgs = calloc(sizeof(double), s->psolver.ns * s->psolver.ms);

		// get process 0 submatrices
		matGetSub(f, s->n, s->m, s->offsetsi[s->cart_id], s->offsetsj[s->cart_id], fgs, s->psolver.ns, s->psolver.ms);
		matSetSub(s->psolver.f, s->psolver.n, s->psolver.m, 1, 1, fgs, s->psolver.ns, s->psolver.ms);
		matGetSub(g, s->n, s->m, s->offsetsi[s->cart_id], s->offsetsj[s->cart_id], fgs, s->psolver.ns, s->psolver.ms);
		matSetSub(s->psolver.g, s->psolver.n, s->psolver.m, 1, 1, fgs, s->psolver.ns, s->psolver.ms);
		free(fgs);

		for (int i = 1; i < s->nproc && s->nproc > 1; ++i) { // get data from other processes
			double *fg = calloc(sizeof(double), s->sizen[i] * s->sizem[i]);
			fgs = calloc(sizeof(double), s->sizens[i] * s->sizems[i]);
			
			matGetSub(f, s->n, s->m, s->offsetsi[i], s->offsetsj[i], fgs, s->sizens[i], s->sizems[i]);
			matSetSub(fg, s->sizen[i], s->sizem[i], 1, 1, fgs, s->sizens[i], s->sizems[i]); // get f submatrix
			MPI_Send(fg, s->sizen[i] * s->sizem[i], MPI_DOUBLE, i, 0, s->cart_comm); // send f
			
			matGetSub(g, s->n, s->m, s->offsetsi[i], s->offsetsj[i], fgs, s->sizens[i], s->sizems[i]);
			matSetSub(fg, s->sizen[i], s->sizem[i], 1, 1, fgs, s->sizens[i], s->sizems[i]); // get g submatrix
			MPI_Send(fg, s->sizen[i] * s->sizem[i], MPI_DOUBLE, i, 0, s->cart_comm); // send g
			
			free(fg);
			free(fgs);
		}
	}
	else {
		MPI_Recv(s->psolver.f, s->psolver.n * s->psolver.m, MPI_DOUBLE, s->cart_root, 0, s->cart_comm, &s->status); // get f from process 0
		MPI_Recv(s->psolver.g, s->psolver.n * s->psolver.m, MPI_DOUBLE, s->cart_root, 0, s->cart_comm, &s->status); // get g from process 0
	}

	free(f);
	free(g);

	return 0;
}

int solverInit(Solver *s, int argc, char *argv[])
{
	int status = 0;

	MPI_Comm_size(MPI_COMM_WORLD, &s->nproc);
	MPI_Comm_rank(MPI_COMM_WORLD, &s->id);

	if (s->id == 0) {
		printf("Parallel MPI Poisson's equation solver, Roni Koitermaa 2022\n");
	
		if (argc < 6) {
			fprintf(stderr, "Usage: %s <infile> <outfile> <n> <gamma> <crit>\n", argv[0]);
			return 2;
		}

		s->inpath = argv[1];
		s->outpath = argv[2];
		s->n = atoi(argv[3]);
		s->psolver.gamma = atof(argv[4]);
		s->crit = atof(argv[5]);

		if (s->n <= 0 || s->psolver.gamma <= 0. || s->crit <= 0. || isnan(s->psolver.gamma) || isnan(s->crit)) {
			fprintf(stderr, "Error! Bad arguments!\n");
			return 3;
		}
		
		for (int i = 1; i < s->nproc && s->nproc > 1; ++i) { // Poisson solver arguments to other processes
			MPI_Send(&s->n, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
			MPI_Send(&s->crit, 1, MPI_DOUBLE, i, 0, MPI_COMM_WORLD);
			MPI_Send(&s->psolver.gamma, 1, MPI_DOUBLE, i, 0, MPI_COMM_WORLD);
		}
	}
	else { // receive arguments
		MPI_Recv(&s->n, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &s->status);
		MPI_Recv(&s->crit, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &s->status);
		MPI_Recv(&s->psolver.gamma, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &s->status);
	}

	s->m = s->n; // assume input is square
	// allocate solution
	s->psolver.f = calloc(sizeof(double), s->n * s->m);
	// allocate RHS
	s->psolver.g = calloc(sizeof(double), s->n * s->m);

	if (s->id == 0) {
		printf("infile: %s, outfile: %s\n", s->inpath, s->outpath);
		s->infile = fopen(s->inpath, "r");
		s->outfile = fopen(s->outpath, "w");

		if (s->infile == NULL || s->outfile == NULL) {
			fprintf(stderr, "Error! Check in/out files!\n");
			return 4;
		}
		
		matFileRead(s->infile, s->psolver.g, s->n, s->m); // read RHS g from file
		
		fclose(s->infile);
		
		// use g as starting solution, get BCs from edges
		memcpy(s->psolver.f, s->psolver.g, sizeof(double) * s->n * s->m);

		for (int i = 0; i < s->n * s->m; ++i) // rescale g to account for system size
			s->psolver.g[i] /= s->n * s->n;
	}

	s->psolver.bcr = 1; // Dirichlet boundary conditions (offset solution edge by one)
	s->psolver.bcl = 1;
	s->psolver.bcu = 1;
	s->psolver.bcd = 1;
	
	s->dims = NULL;
	s->coords = NULL;
	if (s->nproc > 1)
		status = decompose(s);
	else {
		s->psolver.fs = s->psolver.f;
		s->psolver.n = s->n;
		s->psolver.ns = s->n;
		s->psolver.m = s->m;
		s->psolver.ms = s->m;
	}

	printf("id %i/%i, n: %i, m: %i, gamma: %g, crit: %g\n", s->id, s->nproc, s->psolver.ns, s->psolver.ms, s->psolver.gamma, s->crit);

	return status;
}

void solverDestroy(Solver *s)
{
	if (s->id == 0 && s->outfile != NULL) {
		fclose(s->outfile);
	}

	if (s->nproc > 1) {
		free(s->coords);
		free(s->dims);
	
		free(s->coordsi);
		free(s->coordsj);
		free(s->offsetsi);
		free(s->offsetsj);
		free(s->sizen);
		free(s->sizem);
		free(s->sizens);
		free(s->sizems);
	}
	
	free(s->psolver.f);
	free(s->psolver.g);
}

// communicate boundaries between processes
double communicate(Solver *s)
{
	// boundary condition values outside of solution values
	double *left_bc = s->psolver.fs - 1;
	double *right_bc = s->psolver.fs + s->psolver.ms;
	double *up_bc = s->psolver.f + 1;
	double *down_bc = s->psolver.fs + s->psolver.ns * s->psolver.m;

	// outer solution values
	double *left_col = s->psolver.fs;
	double *right_col = s->psolver.fs + s->psolver.ms - 1;
	double *up_row = s->psolver.fs;
	double *down_row = s->psolver.fs + (s->psolver.ns - 1) * s->psolver.m;

	if (s->coords[0] % 2 == 0) { // communicate x boundaries
		/*
		 * O = actual solution values, # = process boundaries
		 *
		 * . . . . #    . O . . .
		 * . . . . #    . O . . .
		 * . . . . # <- . O . . .
		 * . . . . #    . O . . .
		 * . . . . #    . O . . .
		 */
		MPI_Recv(right_bc, 1, s->col_t, s->right_id, 0, s->cart_comm, &s->status);
		
		/*
		 * . . . O .    # . . . .
		 * . . . O .    # . . . .
		 * . . . O . -> # . . . .
		 * . . . O .    # . . . .
		 * . . . O .    # . . . .
		 */
		MPI_Send(right_col, 1, s->col_t, s->right_id, 0, s->cart_comm);
		
		/*
		 * # . . . .    . . . O .
		 * # . . . .    . . . O .
		 * # . . . . <- . . . O .
		 * # . . . .    . . . O .
		 * # . . . .    . . . O .
		 */
		MPI_Recv(left_bc, 1, s->col_t, s->left_id, 0, s->cart_comm, &s->status);
		
		/*
		 * . O . . .    . . . . #
		 * . O . . .    . . . . #
		 * . O . . . -> . . . . #
		 * . O . . .    . . . . #
		 * . O . . .    . . . . #
		 */
		MPI_Send(left_col, 1, s->col_t, s->left_id, 0, s->cart_comm);
	}
	else {
		MPI_Send(left_col, 1, s->col_t, s->left_id, 0, s->cart_comm);
		MPI_Recv(left_bc, 1, s->col_t, s->left_id, 0, s->cart_comm, &s->status);
		
		MPI_Send(right_col, 1, s->col_t, s->right_id, 0, s->cart_comm);
		MPI_Recv(right_bc, 1, s->col_t, s->right_id, 0, s->cart_comm, &s->status);
	}
	if (s->coords[1] % 2 == 0) { // communicate y boundaries
		MPI_Recv(up_bc, 1, s->row_t, s->up_id, 0, s->cart_comm, &s->status);
		MPI_Send(up_row, 1, s->row_t, s->up_id, 0, s->cart_comm);
		
		MPI_Recv(down_bc, 1, s->row_t, s->down_id, 0, s->cart_comm, &s->status);
		MPI_Send(down_row, 1, s->row_t, s->down_id, 0, s->cart_comm);
	}
	else {
		MPI_Send(down_row, 1, s->row_t, s->down_id, 0, s->cart_comm);
		MPI_Recv(down_bc, 1, s->row_t, s->down_id, 0, s->cart_comm, &s->status);
		
		MPI_Send(up_row, 1, s->row_t, s->up_id, 0, s->cart_comm);
		MPI_Recv(up_bc, 1, s->row_t, s->up_id, 0, s->cart_comm, &s->status);
	}

	double res = 0.; // sum of residuals in all processes

	MPI_Reduce(&s->res, &res, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	MPI_Bcast(&res, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

	return res;
}

void solve(Solver *s)
{
	double res = 0.; // relative residual for iteration
	double resold = 0.; // relative residual for last iteration
	int logn = max(1000000/s->n/s->m, 1); // logging period
	double resmax = 1.e6; // maximum residual/iterations to break if exceeded
	int itmax = 1000000;
	
	s->psolver.it = 0; // SOR iteration
	s->it = 0; // solver iteration

	double wtime = MPI_Wtime(); // start wall time
	
	do {
		if (s->nproc == 1) { // serial
			s->res = solveSOR(&s->psolver);
		}
		else { // parallel black red algorithm
			int odd = (s->offsetsi[s->cart_id] + s->offsetsj[s->cart_id]) % 2;
			// solve first sublattice
			if (odd)
				s->res = solveSORR(&s->psolver);
			else
				s->res = solveSORB(&s->psolver);
			
			// communicate process boundaries
			double res1 = communicate(s);

			// solve second sublattice
			if (odd)
				s->res = solveSORB(&s->psolver);
			else
				s->res = solveSORR(&s->psolver);
			
			double res2 = communicate(s);
			
			s->res = res1 + res2; // total residual is sum of both sublattices in all processes
		}

		if (s->id == 0 && logn > 0 && s->psolver.it % logn == 0)
			printf("residual: %g\n", s->res); // relative and (absolute) residual

		if (s->res > resmax || s->it > itmax)
			break;

		++s->it;
	} while (s->res > s->crit); // if residual larger than criterion, solve again

	wtime = MPI_Wtime() - wtime;
	
	if (s->id == 0) {
		if (s->res <= s->crit)
			printf("converged ");
		else
			printf("didn't converge ");
	
		printf("in %i iterations, wall time: %g\n", s->it, wtime);
	}

	if (s->nproc == 1) {
		matFilePrint(s->outfile, s->psolver.f, s->n, s->m);
		
		return;
	}

	MPI_Barrier(MPI_COMM_WORLD);

	// gather all matrices
	if (s->id == 0) {
		double *f = calloc(sizeof(double), s->n * s->m);
		double *fps = calloc(sizeof(double), s->psolver.ns * s->psolver.ms);
		
		matGetSub(s->psolver.f, s->psolver.n, s->psolver.m, 1, 1, fps, s->psolver.ns, s->psolver.ms);
		matSetSub(f, s->n, s->m, s->offsetsi[s->cart_id], s->offsetsj[s->cart_id], fps, s->psolver.ns, s->psolver.ms); // put process 0 matrix in place
		free(fps);
		
		for (int i = 1; i < s->nproc; ++i) { // get data from other processes
			double *fp = calloc(sizeof(double), s->sizen[i] * s->sizem[i]);
			fps = calloc(sizeof(double), s->sizens[i] * s->sizems[i]);
			
			MPI_Recv(fp, s->sizen[i] * s->sizem[i], MPI_DOUBLE, i, 0, s->cart_comm, &s->status); // get submatrix
			matGetSub(fp, s->sizen[i], s->sizem[i], 1, 1, fps, s->sizens[i], s->sizems[i]);
			matSetSub(f, s->n, s->m, s->offsetsi[i], s->offsetsj[i], fps, s->sizens[i], s->sizems[i]); // put received matrix into final matrix
			
			free(fp);
			free(fps);
		}
		
		matFilePrint(s->outfile, f, s->n, s->m); // final assembled matrix to file

		free(f);
	}
	else {
		MPI_Send(s->psolver.f, s->psolver.n * s->psolver.m, MPI_DOUBLE, s->cart_root, 0, s->cart_comm); // send matrix of this process to 0
	}
}

