# Q-ball

Q-balls are [non-topological solitons](https://en.wikipedia.org/wiki/Non-topological_soliton) that arise in complex scalar field theories with a non-linear attractive potential and a global or gauge U(1) symmetry. They have relevance for various cosmological scenarios such as [baryogenesis](https://en.wikipedia.org/wiki/Baryogenesis), [early-Universe phase transitions](https://en.wikipedia.org/wiki/Cosmological_phase_transition), and the [dark matter problem](https://en.wikipedia.org/wiki/Dark_matter).

This repository contains the source code for an adaptive solver of the gauged Q-ball equations of motion in axisymmetry and in three spatial dimensions. The equations are solved in parallel using the [finite-difference method](https://en.wikipedia.org/wiki/Finite_difference_method) with [adaptive mesh refinement](https://en.wikipedia.org/wiki/Adaptive_mesh_refinement) and [multigrid](https://en.wikipedia.org/wiki/Multigrid_method). Further details about the numerical implementation are provided in the following publications:

> M. P. Kinach and M. W. Choptuik, [*Dynamics of U(1) Gauged Q-balls in Three Spatial Dimensions*](https://doi.org/10.1103/PhysRevD.110.075033), Phys. Rev. D **110**, 075033 (2024), [arXiv:2408.07561](https://doi.org/10.48550/arXiv.2408.07561) [hep-th].

> M. P. Kinach and M. W. Choptuik, [*Relativistic Head-on Collisions of U(1) Gauged Q-balls*](https://doi.org/10.1103/PhysRevD.110.015012), Phys. Rev. D **110**, 015012 (2024), [arXiv:2404.04323](https://doi.org/10.48550/arXiv.2404.04323) [hep-th].

> M. P. Kinach and M. W. Choptuik, [*Dynamical Evolution of U(1) Gauged Q-balls in Axisymmetry*](https://doi.org/10.1103/PhysRevD.107.035022), Phys. Rev. D **107**, 035022 (2023), [arXiv:2210.11198](https://doi.org/10.48550/arXiv.2211.11198) [hep-th].

## Requirements

At a minimum, this software requires a Unix-based operating system with the following software stack:

* C and Fortran compilers with MPI support
* [RNPL](https://fpretori.scholar.princeton.edu/group-resources) compiler and libraries
* [PAMR/AMRD](https://fpretori.scholar.princeton.edu/group-resources) libraries
* [Maple](https://maplesoft.com/) with the [FD Toolkit](https://github.com/rmanak/FD/)

Additionally, one may choose to install the [XVS](https://laplace.physics.ubc.ca/Doc/xvs/) and [DV](https://laplace.physics.ubc.ca/Doc/DV/) visualization utilities. The tarballs for installation of RNPL, PAMR/AMRD, XVS, and DV are available on the [UBC Numerical Relativity FTP server](ftp://laplace.phas.ubc.ca/pub/).

For users of Ubuntu Linux, most of the software stack can be installed by following [these instructions for RNPL, XVS, and DV](https://github.com/mkinach/PhD-tutorials/blob/main/RNPL/rnpletal-installguide.md) and [these instructions for PAMR/AMRD](https://github.com/mkinach/PhD-tutorials/blob/main/PAMR/pamr-installguide.md). Maple and the FD Toolkit must be installed separately because Maple is proprietary software which requires a license.

Please be aware that significant computational resources are needed to achieve accurate physical results. Users should not expect to simulate most scenarios of interest on standalone desktop machines. It is recommended to use a high-performance computing environment which can deliver at least 1TB of memory and 32 to 256 CPU cores.

## Running the Code

Two different versions of the code are included in this repository: `axi/` and `3d/`. These correspond to solvers of the equations of motion in axisymmetry and in three spatial dimensions, respectively. The `axi/` directory uses RNPL to implement a second-order Crank-Nicolson method for the time evolution. The `3d/` directory uses the FD Toolkit to implement a fourth-order classic Runge-Kutta method for the time evolution.

Assuming that the software stack is properly installed, you can build and run the main program using:
```
./build-and-run.sh
```
Pass the `-h` flag to this script to list the available options. General simulation parameters can be modified by editing the files in the respective `run/input/` directories. For the axisymmetric case, some simulation parameters can also be modified in `axi/src/rnpl/init_qball.inc`. To use different gauged Q-ball initial data, you can modify the parameters in `run/input/shoot_log.mpl` and `run/input/shoot_poly.mpl`.

Output will be saved in the respective `run/output/` directories as SDF (Scientific Data Format) files. This is a binary format which is designed to be visualized with XVS and DV. If you wish to convert the SDF files to a plaintext format, you can use the `sdfdump` utility which is included with the RNPL installation.

## Sample Results

* Video: [Relativistic Collisions of Gauged Q-balls](https://vimeo.com/1074858134)
* Video: [Dynamical Stability of Gauged Q-balls](https://vimeo.com/1074847373)
