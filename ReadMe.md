# Phase-Field Simulation of Multiphase Fluid Flow with Wetting Boundary Conditions

This repository contains a simulation code for multiphase fluid flow using the phase-field method (PFM), with an improved implementation of the wetting boundary condition for off-grid objects. The proposed method, called the immersed-boundary phase-field implementation (IB-PFI), aims to reduce anisotropic errors arising from the use of a rectangular grid.

## Features

- Simulates droplets adhering to circular objects and capillary flow in a parallel-plate channel.
- Implements the wetting boundary condition for off-grid objects using the immersed-boundary formulation of solid–fluid interfaces.
- Suppresses anisotropic errors and improves agreement with theoretical predictions compared to simulations without IB-PFI.
- Extends the applicability of the PFM to simulations of multiphase fluid flows under various geometric conditions.

## How to Build and Run

1. Compile the code using the following command:
   `mpic++ pfmwrite.cpp -lm -O3`
3. Run the simulation with the desired number of CPU cores (`np_num`). Note that `np_num` should satisfy the condition `[XCELL_NUM % (2*np_num) == 0]`. For example, to run with 2 CPU cores:
`mpirun -np 2 ./a.out`
