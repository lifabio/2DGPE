# 2DGPE Vortex
Gross-Pitajevskii equations sover in 2D through the Crank-Nicholson method with predictor-corrector. Realized in Fortran 95. It solves the system of two-coupled Gross-Pitajevskii equations describing the temporal dynamics of two interacting Bose-Einstein Condensates (BEC) in two dimensions (fixed grid).

This version of the code is capable of handling rotating reference frame, compute expectation of value of Energy, Norm, Momentum and Angular Momentum as well as handle different initial condition like:gaussian wave packet, gaussian vortices, thomas-fermi approximation...etc


- Input Parameter file in the form (Example):
```
 201 201 30000		[Npx Npy Npt]
 100 100 200 50 25	[Npixel_x Npixel_y Nprint printEnerrgy_counter grid_update_counter]
 30. 30. 18.		[Xmax Ymax Tmax]
 1					[flag potential: 1-harmonic...]
 4 3					[flag initial conditon sp1 and sp2]
 0.0					[angular speed of the trap rad/s]
 1.2728 1. 1. 0. 0. 0. 0.		[h_leng sig1 sig2 x01 y01 x02 y02]
 1.0 1.0 1.0 1.0			[wx1 wy1 wx2 wy2]
 1. 1. 10000. 10000. 0.1 0.1 0.0	[m1 m2 N1 N2 g1 g2 g12]
```
- Interactions and Parameters:
	- m1 : atomic mass of species 1
	- N1 : number of particles in BEC of species 1
	- m2 : atomic mass of species 2
	- N2 : number of particles in BEC of species 2
	- g1 : intraspecies interactrion of BEC 1
	- g2 : intraspecies interactrion of BEC 2
	- g12: interspecies interaction between BEC 1 and BEC 2

- Flag initial conditions values:
	- 0 : Psi zero everywhere
	- 2 : Gaussian wave-packet + vortex at the center
	- 3 : Thomas Fermi + vortex at the center
	- 4 : Thomas Fermi approx. for Harmonic Potential
	- 5 : Constant Psi
	- 6 : Constant Psi + Vortex at the center
	- Default : Gaussian wave-packet

- Flag Potential values:
	- 1 : Harmonic potential
	- 2 : Cilindrical potential
	- Default : Harmonic potential  

