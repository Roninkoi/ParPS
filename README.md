# Parallel Poisson's equation Solver

<p align="center">
<img src="https://user-images.githubusercontent.com/12766039/182003413-071048b2-c1a7-4bbd-a82d-74a16c8f306d.png" width=70% height=70%>
<img src="https://user-images.githubusercontent.com/12766039/182003432-e4cdee34-c674-4556-b3fd-61c91abbda6c.png" width=70% height=70%>
</p>

Solve the 2D Poisson's equation

$$
\frac{\partial^2}{\partial x^2} f(x, y) + \frac{\partial^2}{\partial y^2} f(x, y) = g(x, y)
$$

on a square in parallel using OpenMPI. Successive over-relaxation is used to obtain a solution $f$ from RHS $g$. Solution is parallelized using domain decomposition on the unit square, which allows each CPU to solve a piece off the whole problem. Code written in C, with Python plotting. See doc.pdf for further documentation.

Usage:

`./parps <infile> <outfile> <n> <gamma> <crit>`

infile = input file for g as an n x n matrix, edges are boundary conditions for f

infile = solution output file for f as an n x n matrix

n = system size

gamma = over-relaxation parameter ]1, 2[

crit = convergence criterion (sum of all differences between iterations)

## Plotting

To plot the example result matrices and compare them:

./plotmat.py serial.dat parallel.dat

To monitor convergence of solutions from log files in real time, a Gnuplot script can be used:

gnuplot residual.gp

