# Q-ball Solver
 
This is the source code for a solver of the gauged Q-ball equations of motion in three spatial dimensions with adaptive mesh refinement.

This code was used to generate results for the following publication:
> M. P. Kinach and M. W. Choptuik, [*Dynamics of U(1) Gauged Q-balls in Three Spatial Dimensions*](https://doi.org/10.1103/PhysRevD.110.075033), Phys. Rev. D 110, 075033 (2024), [arXiv:2408.07561](https://doi.org/10.48550/arXiv.2408.07561) [hep-th].

## Requirements
 
* C and Fortran compilers with MPI support
* [RNPL and PAMR libraries](http://laplace.physics.ubc.ca/Group/Software.html)
* [Maple](https://maplesoft.com/) with the [FD Toolkit](https://github.com/rmanak/FD/)

## Running the Code

Assuming the relevant compilers and libraries are properly installed, run using:
```
./build-and-run.sh -m evo -n 4 -v log
```
Simulation parameters can be changed by editing the files in the `run/input/` directory.

<!-- 
## Sample Results

* VIDEO: [video](https://vimeo.com/)
* -->
