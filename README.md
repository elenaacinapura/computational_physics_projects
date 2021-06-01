# Computational Physics Projects

Projects for the course of Computational Physics at UniTN. The programs are written in C and plots are done in Gnuplot and Python.

## Requirements
The code in this repository relies on the following libraries of the C language, that must be installed in order to run succesfully:
- GSL (Gnu Scientific Library)
- BLAS (Basic Linear Algebra Subroutine), in particular the OpenBLAS implementation
- LAPACK
  
Furthermore, the programs include several functions from my personal library of numerical methods, which is added as a submodule of this repository. Therefore, make sure to clone this repository together with its submodules, using the `--recurse-submodules` option.

## Compiling
To compile the source code, use the following commands from the uppermost level of the repository:
```
mkdir build
cd build
cmake ..
make -j
```
This compiles all the files in the repository, so it may take a while. You can also skip the last command (`make -j`) and compile-run single files specifically with the instructions provided in the following section.

## Content
Here is a list of the projects of compuational physics contained in this repository.

#### Classical Mechanics
- **Classical scattering** with a Lennard-Jones potential
- Classical scattering in a **dipolar field**
- Motion of a **wave** in a gravitational field

#### Statistical Mechanics
- Computation of the **second virial coefficient**
- **Molecular Dynamics** of a **fluid** of particles interacting with the Lennard-Jones potential
- Calculation of the **structure factor** of a **liquid** of particles interacting with Lennard-Jones potential
- Calculation of the **pair distribution function** of a repulsive liquid using the **Ornstein-Zernicke** equation with the **Percus-Yevick** closure

#### Quantum Mechanics
- **Bound states** of quantum potential
- Solution of the **Schrodinger equation** using the **split-operator** method
- Solution of the **Gross-Pitaevskii** equation for **Bose-Einstein Condensates**
- **Quantum Scattering** with a Lennard-Jones fluid
- Simulation of the **SIR model** for infectious diseases
- **Polytropic model** for the Sun
- **Quantum tunneling** through potential wells

#### Miscellaneous
- Estimating the **age of the universe** with the inflation model
- Implementation of the **variational method** for finding the ground state of a potential