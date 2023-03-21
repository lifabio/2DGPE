# 2DGPE
Gross-Pitajevskii equations sover in 2D through the Crank-Nicholson method with predictor-corrector. Realized in Fortran 95.
It solves the system of two-coupled Gross-Pitajevskii equations describing the temporal dynamics of two interacting Bose-Einstein Condensates (BEC) in two dimensions.

This folder contains different version of the solver featuring the computations of different observables and the possibility to handle different initial conditions.

- BEC2d_2GPE_Vortex_v06 : is capable of handling rotating reference frame, compute expectation of value of Energy, Norm, Momentum and Angular Momentum as well as handle different initial condition like:gaussian wave packet, gaussian vortices, thomas-fermi approximation...etc

- BEC2d_2GPE_AMR_v07 : like the previous one with an Adptive Mesh Refinement.

