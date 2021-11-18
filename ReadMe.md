Computational Fluid Dynamics Simulation for two-phase fluids contacting on solid objects (Liquid-Gas-Solid) using Phase-Field method.
Off-grid strcutre is implemented by Immersed-Boundary method.

**How to build and run**

mpic++ pfmwrite.cpp -lm -O3

mpirun -np 2 ./a.out

The number of cpu cores (np_num) should statisfy [XCELL_NUM % (2*np_num) == 0]. (e.g. The above expression means np_num = 2)

**How to visualize**

File format is ASCII Tecplot (tec).
You can visualize the data using ParaView.

**Acknowlege**

Thanks to Ishida and Deji.
