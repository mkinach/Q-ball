# Global Q-Ball Solver in 2D

This is the source code for a solver of the global (non-gauged) Q-ball equations of motion in two dimensions. The code utilizes the RNPL and PAMR/AMRD libraries for the finite-difference and adaptive mesh refinement implementations, respectively. 

### Requirements
* C and Fortran MPI compilers
* [RNPL and PAMR](http://laplace.physics.ubc.ca/Group/Software.html)

### Running the Code
Assuming the relevant compilers and libraries are properly installed, run using:
```
make clean
make rnpl
make pamr
cat qball-pamr.fparam qball-pamr.rtparam > qball-pamr.param
mpirun -np 4 qball-pamr qball-pamr.param
```

### Sample Results
* VIDEO: [Relativistic Collisions of Q-balls](https://vimeo.com/660839958)

### References
B. Gutierrez, ["Relativistic Scattering of Solitons in Nonlinear Field Theory"](https://dx.doi.org/10.14288/1.0073874), PhD thesis (The University of British Columbia, 2013)
