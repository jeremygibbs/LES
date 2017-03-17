Email jeremy.gibbs@utah.edu for this data. It exceeds GitHub limits on file size.

This folder contains DNS simulation data for decaying isotropic turbulence based on the 
Compte-Bellot Corrsin (aka CBC) classic experimental data.
(see: Comte-Bellot, G., & Corrsin, S. (1971). Simple Eulerian time correlation of full-and narrow-band velocity signals in grid-generated,“isotropic”turbulence. Journal of Fluid Mechanics, 48(02), 273–337).

The initial condition was generated using an in-house isotropic turbulence generator based
on the cbc spectrum. The data was subsequently used as an initial condition in a finite
volume cfd code (Wasatch). 

Grid Specification:
-------------------
Domain size: L = Lx = Ly = Lz = 2*pi/15 m
Number of grid points: Nx = Ny = Nz = 268
Grid size:   h = dx = dy = dz = L/268
density: 1.0 kg/m3
dynamic viscosity: 2.0e-5 kg/m.s

Files:
------
*cbcspectrum.mat: File containing the original CBC spectrum taken at 3 snapshots in time: 
t = 0, 0.28, and 0.66 s. The first column in the data contains the wave numbers (1/m).

* p,u,v,w_t0.0s.mat: initial condition for pressure and velocity field.

* p,u,v,w_t0.28s.mat: pressure and velocity fields at t = 0.28 s.

* p,u,v,w_t0.66s.mat: pressure and velocity fields at t = 0.66 s.


IMPORTANT:
----------
The velocity field is staggered in space, in the -x, -y, and -z directions.
