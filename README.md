# Phase-Field Fracture Optimal Control

The provided code refers to the paper 
[Space-time formulation, discretization, and computational performance studies for phase-field fracture optimal control problems: reproduction code Section 5.1.1](https://doi.org/10.1016/j.jcp.2022.111554)
published in 2022. The code we provide is completely novel and offers a concrete implementation of the mathematics developed in the paper, wherein no implementation aspects were covered.

# Summary

The provided codebase is designed to solve phase-field fracture optimal control 
problems, wherein the goal is to achieve a desired fracture in a brittle material
thorugh the application of external forces. This algorithmic framework was developed alongside our recent publication
[1] and enables the accurate and efficient simulation of phase-field optimal control problems in a space-time fashion. 
In this example, the fracture is controlled through Neumann boundary conditions, and the cost functioncal tracks the distance between.
Our code is based on the open source libraries DOpElib [2] and deal.II [3].

# Installation instructions and executing the code

1. Install deal.II 9.5.1 via www.dealii.org \
Download: https://www.dealii.org/download.html \
Installation instructions: https://www.dealii.org/current/readme.html 

2. Install DOpElib via http://www.dopelib.net \
Download: https://github.com/winnifried/dopelib \
Installation instructions: https://winnifried.github.io/dopelib/documentation.html 
in Chapter 2 of the *.pdf manual. 

3. Please put the current code into some new folder on your machine. \
Follow instructions from DOpElib *.pdf manual in Chapter 4 
(i.e, Section 4.4 Creating new examples)
to set up all environment variables (for finding deal.II and DOpElib) correctly. \
Then: build, compile, run as described in Section 4.4 of the dopelib manual.
To this end, building and compiling yields the executable JCP_5_1_1 that can
be run in the (Linux) terminal via ./JCP_5_1_1 
by taking implicitly the parameter file dope.prm into account.

4. The results of this code (see local folder Results/ (this name given in dope.prm (bottom),
which can be thus changed when multiple simulations shall be run simultaneously) ) should then reproduce 
Example 1 (Section 5.1.1) of Khimin et al., JCP, 2022. To compare the terminal output of the current
implementation with own runs, the log file dope_Aug_12_2024.log can be used.

# First steps

As the first steps, we recommend to read [1] first in order to understand the algorithmic background. Next, C++ knowledge in order to handle DOpElib and deal.II is required. Afterwards, the point of departure is [main.cc](main.cc). Therein, all include files from DOpElib and deal.II are taken as well as all the basic functionality is set up. The runtime parameters are adjusted in [dope.prm](dope.prm).
Running the code follows with the instructions previously given.

# Automated testing

# Documentation

[main.cc](main.cc): In this file, first deal.II and DOpElib libraries are included. Then, local files such as [localpde.h](localpde.h), [localfunctional.h](localfunctional.h), [functionals.h](functionals.h), and [my_functions.h](my_functions.h) are included.


[dope.prm](dope.prm): This is the runtime parameter file. Therein, global and local parameters can be adjusted without compiling the code again. These values are model, regularization, material, numerical, optimization parameters. Moreover, the output directory can be chosen.

[localpde.h](localpde.h): In this file, the actual model, i.e., PDE system is implemented. Here, it is a quasi-static phase-field fracture system with two unknowns: vector-valued displacements and a scalar-valued smooth indicator function. The latter is subject to an inequality constraint in time, i.e., crack irreversibility. This is treated with the help of simple penalization.

[localfunctional.h](localfunctional.h): In this file, the optimization objective is implemented in terms of the cost functional. 

[functionals.h](functionals.h): In this file, additional output functionals are implemented, so-called quantities of interest that are specifically interesting to be observed from an physics or engineering point of view.

[my_functions.h](my_functions.h): In this file, problem-specific functions such as non-homogeneous Dirichlet boundary conditions are implemented.



# Contributing guidelines

For contributing to this library, please see [here](CONTRIBUTING.md).

# License

The license is GNU LESSER GENERAL PUBLIC LICENSE (LGPL) Version 2.1. Detailed information can be found [here](LICENSE).

# References

[[1]](https://doi.org/10.1016/j.jcp.2022.111554) Denis Khimin, Marc C. Steinbach, Thomas Wick; Space-time formulation, discretization, and computational performance studies for phase-field fracture optimal control problems,
Journal of Computational Physics (JCP), Vol. 470, 2022.

[[2]](https://doi.org/10.11588/ans.2017.2.11815) Christian Goll, Thomas Wick, Winnifried Wollner; DOpElib: Differential Equations and Optimization Environment; A Goal Oriented Software
Library for Solving PDEs and Optimization Problems with PDEs, Archive of Numerical Software (ANS), Vol. 5(2),  pp. 1-14, 2017.

[[3]](https://doi.org/10.1515/jnma-2023-0089) Daniel Arndt, Wolfgang Bangerth, Maximilian Bergbauer, Marco Feder, Marc Fehling, Johannes Heinz, Timo Heister, Luca Heltai, Martin Kronbichler, Matthias Maier, Peter Munch, Jean-Paul Pelteret, Bruno Turcksin, David Wells, Stefano Zampini. The deal.II Library, Version 9.5, Journal of Numerical Mathematics, vol. 31, no. 3, pages 231-246, 2023. 
