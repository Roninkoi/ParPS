# Parallel Poisson's equation Solver

![Screenshot from 2022-07-31 02-26-40](https://user-images.githubusercontent.com/12766039/182003413-071048b2-c1a7-4bbd-a82d-74a16c8f306d.png)
![Screenshot from 2022-07-31 02-26-09](https://user-images.githubusercontent.com/12766039/182003432-e4cdee34-c674-4556-b3fd-61c91abbda6c.png)

Solve the 2D Poisson's equation
$$
\frac{\partial^2}{\partial x^2} f(x, y) + \frac{\partial^2}{\partial y^2} f(x, y) = g(x, y)
$$
on a square in parallel using OpenMPI. Code written in C, with Python plotting. See doc.pdf for further documentation.

Usage:
`./parps <infile> <outfile> <n> <gamma> <crit>`

infile = input file for g as an n x n matrix, edges are boundary conditions for f
infile = solution output file for f as an n x n matrix
n = system size
gamma = over-relaxation parameter ]1, 2[
crit = convergence criterion (sum of all differences between iterations)

## Plotting

To plot the example result matrices and compare them:

../plotmat.py serial.dat parallel.dat

To monitor convergence of solutions from log files in real time, a Gnuplot script can be used:

gnuplot residual.gp

