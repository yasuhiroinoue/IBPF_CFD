# Phase-Field Simulation of Multiphase Fluid Flow with Wetting Boundary Conditions

This repository contains a simulation code for multiphase fluid flow using the phase-field method (PFM), with an improved implementation of the wetting boundary condition for off-grid objects. The proposed method, called the immersed-boundary phase-field implementation (IB-PFI), aims to reduce anisotropic errors arising from the use of a rectangular grid.

## Mathematical Model

- Simulates droplets adhering to circular objects and capillary flow in a parallel-plate channel.
- Implements the wetting boundary condition for off-grid objects using the immersed-boundary formulation of solid–fluid interfaces.
- Suppresses anisotropic errors and improves agreement with theoretical predictions compared to simulations without IB-PFI.
- Extends the applicability of the PFM to simulations of multiphase fluid flows under various geometric conditions.

The details of the mathematical model are described in the following paper:

Inoue, Y., Ishida, K., Takada, N., & Hojo, M. (2015). Reductions in Anisotropic Errors from Implementation of Phase-Field Wetting Boundary Condition for Off-Grid Objects. International Journal of Computational Methods, 12(6), 1550042 (2015). [https://doi.org/10.1142/S0219876215500425](https://doi.org/10.1142/S0219876215500425)

## How to Build and Run

1. Compile the code using the following command:
   `mpic++ pfmwrite.cpp -lm -O3`
3. Run the simulation with the desired number of CPU cores (`np_num`). Note that `np_num` should satisfy the condition `[XCELL_NUM % (2*np_num) == 0]`. For example, to run with 2 CPU cores:
`mpirun -np 2 ./a.out`

## How to Visualize

The simulation data is stored in the ASCII Tecplot (.tec`) format. You can visualize the data using ParaView, which can be obtained from [https://www.paraview.org/](https://www.paraview.org/).

## Acknowledgements

The development of this code was led by Ishida Kazuki, with contributions from Deji Takeji

