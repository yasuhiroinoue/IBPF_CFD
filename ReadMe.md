Computational Fluid Dynamics Simulation for two-phase fluids contacting on solid objects (Liquid-Gas-Solid) using Phase-Field method.
Off-grid strcutre is implemented by Immersed-Boundary method.

**How to build and run**

mpic++ pfmwrite.cpp -lm -O3

mpirun -np 2 ./a.out

**How to visualize**

File format is ASCII Tecplot (tec).
You can visualize the data using ParaView.

**Acknowlege**

Thanks to Ishida and Deji.
