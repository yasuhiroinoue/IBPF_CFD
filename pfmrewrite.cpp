/*
//
// For contact information
// Yasuhiro Inoue
// inoue.yasuhiro.4n@kyoto-u.ac.jp
// All comments are translated to English (2024.03.29)
//
// A super minor revision is applied on Nov 17th, 2021
// Extention 'plt' is renamed to 'tec', which indicates ASCII Tecplot format.
// You can visualize the results using ParaView
//
// The main code is written by Ishida Kazuki
// The sweepout is written by Deji Takeji
//
*/

//////////////////////////////////////////////
///////			Ishida Kazuki		//////////
///////	 updated 2010/10/07  		//////////
///////			IB method +PF method//////////
///////			rewrite				//////////
///////			staggered lattice	//////////
///////								//////////
//////////////////////////////////////////////

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <iostream>
#include <math.h>
#include "mpi.h"
#define N_ARRAY 3
#include "sweepout.h"
#undef N_ARRAY


//#define WALLX
//#define WALLY
#define USEMPI //MPI

#define EPS1 1.0e-8
#define EPS2 1.0e-7
#define EPS3 1.0e-8 //1.0e-8 /*Minimum allowable error for SOR*/
#define MAX 1000000 /*Maximum number of SOR iterations*/
#define XMAX 96.0
#define YMAX 120.0
#define XCELL_NUM 96 //Number of cells
#define YCELL_NUM 120 //Number of cells
#define dt 1.25e-2 //Time step width
#define RHOL 1.0 //Density of liquid
#define RHOG 1.25e-3 //Density of gas
#define RHOS 0.5 //Internal density of object
#define MUL 5.12e-1//1.67e-1 /*Viscosity coefficient*/
#define MUG 6.9375e-3//1.3875e-3//2.26e-3 /*Viscosity of gas*/
#define MUS 0.01 //Viscosity inside the object (not used)
#define PHImax 0.405
#define PHImin 0.265
#define PHIG 0.275
#define PHIL 0.380
#define TMAX 1.0
#define gamma 12.0 //Mobility
#define gammas -1.0e-3 //Wettability
#define vdWa 1.0 //Constant of vdW Model
#define vdWb 1.0 //Constant of vdW Model
#define vdWT 0.293 //Temperature of vdW Model
#define kappa1 1.3778//2.94e-2
#define kappa2 0.1
#define gravity 0.0//1.5e-4
#define ALPHA 0.8 //SOR relaxation factor
#define ALPHA2 0.8//1.0
#define ALPHA3 0.8//1.2
#define fiberx 32.0
#define fibery 32.0
#define Diameter 32.0
#define diameter2 144.0
//#define KAIDAN
#define HEIGHT 0 //Used for making a staircase-like shape
#define NORMERR 1.0e-5
#define FLAGERR -10e-5
#define FIBER_NUM 1 //Number of fibers
#define LJE 5.0e-7 //Potential used to avoid fiber collision
#define SIGMAPLUS 3.0 //Added to the fiber diameter to make the collision diameter
#define CUTOFF 3.0 //Cutoff distance from the collision diameter
#define L 32.0
#define H0 100.0
#define DEGREE (180.0-63.0)
#define LRADIUS 20.0
#define XVALUE 48.0
#define YVALUE (30.0+LRADIUS*(1.0+sqrt((1.0-DEGREE/180.0)+(sin(M_PI*DEGREE/180.0))/M_PI)))


struct _cell_type {
    double v_x; // Velocity in x direction
    double v_y; // Velocity in y direction
    double press; // Pressure
    double v_xhypo; // Hypothetical velocity in x direction
    double v_yhypo; // Hypothetical velocity in y direction
    double hx; // Increment of vx by Adams Bashforth method for the previous step
    double hy; // Increment of vy by Adams Bashforth method for the previous step
    double x; // x coordinate
    double y; // y coordinate
    double phi; // Order parameter
    double phi_x; // Value on the boundary in x direction
    double phi_y; // Value on the boundary in y direction
    double hf; // Update of phi by Adams Bashforth method for the previous step
    double rho; // Density
    double rho_x;
    double rho_y;
    double vof; // Volume fraction of the object
    double vof_x;
    double vof_y;
    double ls; // Distance from the surface of the object (Level-set function)
    double ls_x;
    double ls_y;
};


struct _fiber_type{
    double x;
    double y;
    double radius; // Radius
    double v_x;
    double v_y;
    double v_xold; // Velocity of the previous step
    double v_yold;
    double angle; // Angle of the object
    double angvel; // Angular velocity
    double angvelold; // Angular velocity of the previous step
    double hx; // Increment of vx by Adams Bashforth method for the previous step
    double hy;
    double ha; // Increment of angvel by Adams Bashforth method for the previous step
};

struct _flag_type { // Flag indicating the number of cells from the surface of the object
    int c;
    int x;
    int y;
    int fiber; // Indicates the presence of a fiber
    int fiber_x;
    int fiber_y;
};



_flag_type* FLAG;
_fiber_type* FIBER;

struct _mu_type{ // Viscosity
    double c; // Center
    double x; // On the center of the lattice face
    double y;
    double xy; // At the midpoint of the lattice edge
};

struct _vold_type{ // For holding the velocity before updating
    double x;
    double y;
};


double	dx = XMAX / (double) XCELL_NUM;
double	dy = YMAX / (double) YCELL_NUM;
int CELL_SIZE = (XCELL_NUM+6)*(YCELL_NUM+6);

using namespace std;

int n,nout,TCELL_NUM;
double h_inold;

void solidlagpoint(_fiber_type* fiberp,int fibernumber); // Calculation of forces acting on the object, Stage 1
void solidlagpointsecond(_fiber_type* fiberp,int fibernumber); // Calculation of forces acting on the object, Stage 2
void solidmove(void); // Calculation of object movement using Crank Nicolson method
void flagdefine(void); // Determination of volume fraction and flags
void flagbnd(void); // Flag boundary conditions
void vhypobnd(void); // vhypo boundary conditions
void kyokai(void); // Velocity boundary conditions
void vofbnd(void); // Boundary conditions for object volume fraction
void pressbnd(void); // Pressure boundary conditions
void phibnd(void); // phi boundary conditions
void phikyokaibnd(void); // phi boundary conditions on the cell boundary
void phid2bnd(double phid2[(XCELL_NUM+6)*(YCELL_NUM+6)]); // Laplacian phi boundary conditions
void ddrhodefine(void); // Determination of Laplacian rho
void ddrhobnd(void); // Laplacian rho boundary conditions
void Fbnd(double Fx[(XCELL_NUM+6)*(YCELL_NUM+6)],double  Fy[(XCELL_NUM+6)*(YCELL_NUM+6)]); // phi Flux boundary conditions
void nseq(void); // Navier-Stokes equations
void filewrite(void); // Output for Tecplot
void koshin(void); // SOR method for velocity and pressure correction
void rhodefine(void); // Determination of rho
void mudefine(void); // Determination of mu
void CH_eq(void); // Cahn-Hilliard equation
void restartwrite(void); // Creation of restart file
void restartread(void); // Reading of restart file
void indexsetdefine(int i, int j); // Definition of indexes for accessing structures
void indexSORdefine(int i, int j);
void indexnseq(int i, int j);
void indexbndx(int j);
void indexbndy(int i);
double interphi(int cellindex,int phiflag); // Function for extrapolation around the object surface
double ljpower(double sigma,double distance); // Lennard-Jones potential
void discount(void); // Output as necessary
#ifdef USEMPI
void initial(int my_rank); // Initial conditions
void SOR(int chkflag, int my_rank); // SOR method
void hontai(int chkflag, int my_rank);
void projection(int chkflag,int my_rank);
#else
void initial(void);
void SOR(int chkflag);
void hontai(int chkflag);
void projection(int chkflag);
#endif

_cell_type* CELL;
_vold_type* VOLD; // For holding the velocity before updating

_mu_type* TMPMU;
double* wetphi; // -wetting potential/kappa2
double* Fx;
double* Fy;
double* phid2;
double* ddrho;
int fibernum = FIBER_NUM; // Preparation for creating the number of objects, including virtual objects placed in the case of periodic boundaries
double wetphiwall; // Wettability of the wall surface
double* fiberfx;
double* fiberfy;
double* fibern;

// For communication in MPI
double* mail2;
double* mail1;
int process_num;

// Index for accessing structure members and preparation of structures
int index0;
int index_xm;
int index_xp;
int index_yp;
int index_ym;
int index_xpp;
int index_ypp;
int index_xmm;
int index_ymm;
int index_xm_yp;
int index_xp_ym;
int index_x0;
int index_x1;
int index_x2;
int index_x3;
int index_x4;
int index_x5;
int index_xcp0;
int index_xcp1;
int index_xcp2;
int index_xcp3;
int index_xcp4;
int index_xcp5;
int index_x6;
			
_cell_type* cellp_x0;
_cell_type* cellp_x1;
_cell_type* cellp_x2;
_cell_type* cellp_x3;
_cell_type* cellp_x4;
_cell_type* cellp_x5;
_cell_type* cellp_xcp0;
_cell_type* cellp_xcp1;
_cell_type* cellp_xcp2;
_cell_type* cellp_xcp3;
_cell_type* cellp_xcp4;
_cell_type* cellp_xcp5;
_cell_type* cellp_x6;

int index_y0;
int index_y1;
int index_y2;
int index_y3;
int index_y4;
int index_y5;
int index_y6;
int index_ycp0;
int index_ycp1;
int index_ycp2;
int index_ycp3;
int index_ycp4;
int index_ycp5;
			
_cell_type* cellp_y0;
_cell_type* cellp_y1;
_cell_type* cellp_y2;
_cell_type* cellp_y3;
_cell_type* cellp_y4;
_cell_type* cellp_y5;
_cell_type* cellp_y6;
_cell_type* cellp_ycp0;
_cell_type* cellp_ycp1;
_cell_type* cellp_ycp2;
_cell_type* cellp_ycp3;
_cell_type* cellp_ycp4;
_cell_type* cellp_ycp5;



_cell_type* cellp;
_cell_type* cellp_xp;
_cell_type* cellp_yp;
_cell_type* cellp_xm;
_cell_type* cellp_ym;
_cell_type* cellp_xpp;
_cell_type* cellp_ypp;
_cell_type* cellp_xmm;
_cell_type* cellp_ymm;
_cell_type* cellp_xm_yp;
_cell_type* cellp_xp_ym;
_mu_type* mup;
_mu_type* mup_yp;
_mu_type* mup_xm;
_mu_type* mup_xp;
_mu_type* mup_ym;

void indexsetdefine(int i, int j){
	index0     = j  + i*(YCELL_NUM+6);
	index_xm   = j + (i-1)*(YCELL_NUM+6);
	index_xp  =  j     + (i+1)*(YCELL_NUM+6);
	index_yp  =  (j+1) + i*(YCELL_NUM+6);
	index_ym  =  (j-1) + i*(YCELL_NUM+6);
	index_xpp =  j     + (i+2)*(YCELL_NUM+6);
	index_ypp =  (j+2) + i*(YCELL_NUM+6);
	index_xmm =  j     + (i-2)*(YCELL_NUM+6);
	index_ymm =  (j-2) + i*(YCELL_NUM+6);
	
	
	cellp= &CELL[index0];
	cellp_xp = &CELL[index_xp];
	cellp_yp = &CELL[index_yp];
	cellp_xm = &CELL[index_xm];
	cellp_ym = &CELL[index_ym];
	cellp_xpp = &CELL[index_xpp];
	cellp_ypp = &CELL[index_ypp];
	cellp_xmm = &CELL[index_xmm];
	cellp_ymm = &CELL[index_ymm];
}

void indexSORdefine(int i, int j){
	index0   = j     + i*(YCELL_NUM+6);
	index_xm = j     + (i-1)*(YCELL_NUM+6);
	index_xp = j     + (i+1)*(YCELL_NUM+6);
	index_yp = (j+1) + i*(YCELL_NUM+6);
	index_ym = (j-1) + i*(YCELL_NUM+6);
					
	cellp= &CELL[index0];
	cellp_xp = &CELL[index_xp];
	cellp_yp = &CELL[index_yp];
	cellp_xm = &CELL[index_xm];
	cellp_ym = &CELL[index_ym];
}

	
void indexnsdefine(int i, int j){
	index0      = j     + i*(YCELL_NUM+6);
	index_xp    = j     + (i+1)*(YCELL_NUM+6);
	index_xpp   = j     + (i+2)*(YCELL_NUM+6);
	index_xm    = j     + (i-1)*(YCELL_NUM+6);
	index_xmm   = j     + (i-2)*(YCELL_NUM+6);
	index_yp    = (j+1) + i*(YCELL_NUM+6);
	index_ypp   = (j+2) + i*(YCELL_NUM+6);
	index_ym    = (j-1) + i*(YCELL_NUM+6);
	index_ymm   = (j-2) + i*(YCELL_NUM+6);
	index_xm_yp = (j+1) + (i-1)*(YCELL_NUM+6);
	index_xp_ym = (j-1) + (i+1)*(YCELL_NUM+6);
				
	cellp       = &CELL[index0];
	cellp_xp    = &CELL[index_xp];
	cellp_yp    = &CELL[index_yp];
	cellp_xm    = &CELL[index_xm];
	cellp_ym    = &CELL[index_ym];
	cellp_xpp   = &CELL[index_xpp];
	cellp_xmm   = &CELL[index_xmm];
	cellp_ypp   = &CELL[index_ypp];
	cellp_ymm   = &CELL[index_ymm];
	cellp_xm_yp = &CELL[index_xm_yp];
	cellp_xp_ym = &CELL[index_xp_ym];
	mup    = &TMPMU[index0];
	mup_yp = &TMPMU[index_yp];
	mup_xm = &TMPMU[index_xm];
	mup_xp = &TMPMU[index_xp];
	mup_ym = &TMPMU[index_ym];
}

void indexbndx(int j){
#ifdef KAIDAN
	index_x0   = j;
	index_x1   = j + 1*(YCELL_NUM+6);
	index_x2   = j + 2*(YCELL_NUM+6);
	index_x3   = j + 3*(YCELL_NUM+6);
	index_x4   = j + 4*(YCELL_NUM+6);
	index_x5   = j + 5*(YCELL_NUM+6);
	index_x6   = j + 6*(YCELL_NUM+6);
	
		if((j+HEIGHT)<(YCELL_NUM+3)) {
			index_xcp0 = HEIGHT + j + XCELL_NUM*(YCELL_NUM+6);
			index_xcp1 = HEIGHT + j + (XCELL_NUM+1)*(YCELL_NUM+6);
			index_xcp2 = HEIGHT + j + (XCELL_NUM+2)*(YCELL_NUM+6);
			index_xcp3 = HEIGHT + j + (XCELL_NUM+3)*(YCELL_NUM+6);
			index_xcp4 = HEIGHT + j + (XCELL_NUM+4)*(YCELL_NUM+6);
			index_xcp5 = HEIGHT + j + (XCELL_NUM+5)*(YCELL_NUM+6);
		}
		else {
			index_xcp0 = HEIGHT-YCELL_NUM + j + XCELL_NUM*(YCELL_NUM+6);
			index_xcp1 = HEIGHT-YCELL_NUM + j + (XCELL_NUM+1)*(YCELL_NUM+6);
			index_xcp2 = HEIGHT-YCELL_NUM + j + (XCELL_NUM+2)*(YCELL_NUM+6);
			index_xcp3 = HEIGHT-YCELL_NUM + j + (XCELL_NUM+3)*(YCELL_NUM+6);
			index_xcp4 = HEIGHT-YCELL_NUM + j + (XCELL_NUM+4)*(YCELL_NUM+6);
			index_xcp5 = HEIGHT-YCELL_NUM + j + (XCELL_NUM+5)*(YCELL_NUM+6);
//		}
	}
#else
	index_x0   =  j;
	index_x1   =  j + (YCELL_NUM+6);
	index_x2   =  j + 2*(YCELL_NUM+6);
	index_x3   =  j + 3*(YCELL_NUM+6);
	index_x4   =  j + 4*(YCELL_NUM+6);
	index_x5   =  j + 5*(YCELL_NUM+6);
	index_x6   =  j + 6*(YCELL_NUM+6);
	index_xcp0 =  j + XCELL_NUM*(YCELL_NUM+6);
	index_xcp1 =  j + (XCELL_NUM+1)*(YCELL_NUM+6);
	index_xcp2 =  j + (XCELL_NUM+2)*(YCELL_NUM+6);
	index_xcp3 =  j + (XCELL_NUM+3)*(YCELL_NUM+6);
	index_xcp4 =  j + (XCELL_NUM+4)*(YCELL_NUM+6);
	index_xcp5 =  j + (XCELL_NUM+5)*(YCELL_NUM+6);
#endif
	
	cellp_x0   = &CELL[index_x0];
	cellp_x1   = &CELL[index_x1];
	cellp_x2   = &CELL[index_x2];
	cellp_x3   = &CELL[index_x3];
	cellp_x4   = &CELL[index_x4];
	cellp_x5   = &CELL[index_x5];
	cellp_x6   = &CELL[index_x6];
	cellp_xcp0 = &CELL[index_xcp0];
	cellp_xcp1 = &CELL[index_xcp1];
	cellp_xcp2 = &CELL[index_xcp2];
	cellp_xcp3 = &CELL[index_xcp3];
	cellp_xcp4 = &CELL[index_xcp4];
	cellp_xcp5 = &CELL[index_xcp5];
}


void indexbndy(int i){
	index_y0   =                    i*(YCELL_NUM+6);
	index_y1   = 1 + i*(YCELL_NUM+6);
	index_y2   = 2 + i*(YCELL_NUM+6);
	index_y3   = 3 + i*(YCELL_NUM+6);
	index_y4   = 4 + i*(YCELL_NUM+6);
	index_y5   = 5 + i*(YCELL_NUM+6);
	index_y6   = 6 + i*(YCELL_NUM+6);
	index_ycp0 = YCELL_NUM + i*(YCELL_NUM+6);
	index_ycp1 = (YCELL_NUM+1) + i*(YCELL_NUM+6);
	index_ycp2 = (YCELL_NUM+2) + i*(YCELL_NUM+6);
	index_ycp3 = (YCELL_NUM+3) + i*(YCELL_NUM+6);
	index_ycp4 = (YCELL_NUM+4) + i*(YCELL_NUM+6);
	index_ycp5 = (YCELL_NUM+5) + i*(YCELL_NUM+6);
				
	cellp_y0   = &CELL[index_y0];
	cellp_y1   = &CELL[index_y1];
	cellp_y2   = &CELL[index_y2];
	cellp_y3   = &CELL[index_y3];
	cellp_y4   = &CELL[index_y4];
	cellp_y5   = &CELL[index_y5];
	cellp_y6   = &CELL[index_y6];
	cellp_ycp0 = &CELL[index_ycp0];
	cellp_ycp1 = &CELL[index_ycp1];
	cellp_ycp2 = &CELL[index_ycp2];
	cellp_ycp3 = &CELL[index_ycp3];
	cellp_ycp4 = &CELL[index_ycp4];
	cellp_ycp5 = &CELL[index_ycp5];
}


void flagdefine(void){
	for(int f=0; f<FIBER_NUM;f++){
		_fiber_type* fiberp = &FIBER[f];
		
		if(fiberp->x<0.0) fiberp->x+=(double)XCELL_NUM*dx;
		else if(fiberp->x>((double)XCELL_NUM*dx)) fiberp->x-=(double)XCELL_NUM*dx;
		if(fiberp->y<0.0) fiberp->y+=(double)YCELL_NUM*dy;
		else if(fiberp->y>((double)YCELL_NUM*dy)) fiberp->y-=(double)YCELL_NUM*dy;
	}

#ifndef WALLX
	for(int f=0; f<FIBER_NUM;f++){
		FIBER[f+FIBER_NUM]=FIBER[f];
		
		_fiber_type* fiberp = &FIBER[f+FIBER_NUM];
		
		if(fiberp->x<((double)XCELL_NUM*dx*0.5)) fiberp->x+=XCELL_NUM*dx;
		else fiberp->x-=XCELL_NUM*dx;
	}
	#ifndef WALLY
	for(int f=0; f<FIBER_NUM;f++){
		FIBER[f+2*FIBER_NUM]=FIBER[f];
		
		_fiber_type* fiberp = &FIBER[f+2*FIBER_NUM];
		
		if(fiberp->y<((double)YCELL_NUM*dy*0.5)) fiberp->y+=YCELL_NUM*dy;
		else fiberp->y-=YCELL_NUM*dy;
	}
	for(int f=0; f<FIBER_NUM;f++){
		FIBER[f+3*FIBER_NUM]=FIBER[f];
		
		_fiber_type* fiberp = &FIBER[f+3*FIBER_NUM];
		
		if(fiberp->x<((double)XCELL_NUM*dx*0.5)) fiberp->x+=XCELL_NUM*dx;
		else fiberp->x-=XCELL_NUM*dx;
		if(fiberp->y<((double)YCELL_NUM*dy*0.5)) fiberp->y+=YCELL_NUM*dy;
		else fiberp->y-=YCELL_NUM*dy;
	}
	#endif
#else
	#ifndef WALLY
		for(int f=0; f<FIBER_NUM;f++){
			FIBER[f+FIBER_NUM]=FIBER[f];
			
			_fiber_type* fiberp = &FIBER[f+FIBER_NUM];
			
			if(fiberp->y<((double)YCELL_NUM*dy*0.5)) fiberp->y+=YCELL_NUM*dy;
			else fiberp->y-=YCELL_NUM*dy;
		}
	#endif
#endif

	for(int i=0; i<(XCELL_NUM+6); i++){
		for(int j=0; j<(YCELL_NUM+6); j++){
			index0   = j     + i*(YCELL_NUM+6);
			cellp    = &CELL[index0];
			
			int defflag=0;
			int defflag_x=0;
			int defflag_y=0;
			int ftarget,ftarget_x,ftarget_y;
			double upx,downx,upy,downy;
			double distance;
			double mindis=1000.0*(dx+dy);
			double mindis_x=1000.0*(dx+dy);
			double mindis_y=1000.0*(dx+dy);
			_flag_type* flagp    = &FLAG[index0];
			
			for(int f=0; f<fibernum;f++){
				_fiber_type* fiberp = &FIBER[f];
				
				distance=sqrt((cellp->x-fiberp->x)*(cellp->x-fiberp->x)+(cellp->y-fiberp->y)*(cellp->y-fiberp->y))-fiberp->radius;
				if(mindis>distance){
					mindis=distance;
					ftarget=f;
				}
			
				distance=sqrt((cellp->x-0.5*dx-fiberp->x)*(cellp->x-0.5*dx-fiberp->x)+(cellp->y-fiberp->y)*(cellp->y-fiberp->y))-fiberp->radius;
				if(mindis_x>distance){
					mindis_x=distance;
					ftarget_x=f;
				}
				
				distance=sqrt((cellp->x-fiberp->x)*(cellp->x-fiberp->x)+(cellp->y-0.5*dy-fiberp->y)*(cellp->y-0.5*dy-fiberp->y))-fiberp->radius;
				if(mindis_y>distance){
					mindis_y=distance;
					ftarget_y=f;
				}
			}
			
			cellp->ls=mindis;
			cellp->ls_x=mindis_x;
			cellp->ls_y=mindis_y;
			flagp->fiber=ftarget;
			flagp->fiber_x=ftarget_x;
			flagp->fiber_y=ftarget_y;
		}
	}
	
	for(int i=0; i<(XCELL_NUM+6); i++){
		for(int j=0; j<(YCELL_NUM+6); j++){
			index0   = j     + i*(YCELL_NUM+6);
			cellp    = &CELL[index0];
			
			_flag_type* flagp    = &FLAG[index0];
			
			
			if(cellp->ls>0.0) {
				flagp->c=1;
				cellp->vof=0.0;
			}
			else {
				flagp->c=0;
				cellp->vof=1.0;
			}
			
			if(cellp->ls_x>0.0) {
				flagp->x=1;
				cellp->vof_x=0.0;
			}
			else {
				flagp->x=0;
				cellp->vof_x=1.0;
			}
			
			if(cellp->ls_y>0.0) {
				flagp->y=1;
				cellp->vof_y=0.0;
			}
			else {
				flagp->y=0;
				cellp->vof_y=1.0;
			}
		}
	}
	
	for(int i=0; i<(XCELL_NUM+6); i++){
		for(int j=0; j<(YCELL_NUM+6); j++){
			index0   = j     + i*(YCELL_NUM+6);
			cellp    = &CELL[index0];
			
			int defflag=0;
			int defflag_x=0;
			int defflag_y=0;
			int smtarget,smtarget_x,smtarget_y;
			double mindis=dx*dx+dy*dy;
			double mindis_x=dx*dx+dy*dy;
			double mindis_y=dx*dx+dy*dy;
			double upx,downx,upy,downy;
			double distance;
			
			double targetx,targety;
			double nvec_x,nvec_y;
			_flag_type* flagp=&FLAG[index0];
			_fiber_type* fiberp=&FIBER[flagp->fiber];
			
			if((cellp->ls >-1.0*dx)&&(cellp->ls<1.0*dx)){
				nvec_x=(cellp->x-fiberp->x)/sqrt((cellp->x-fiberp->x)*(cellp->x-fiberp->x)+(cellp->y-fiberp->y)*(cellp->y-fiberp->y));
				nvec_y=(cellp->y-fiberp->y)/sqrt((cellp->x-fiberp->x)*(cellp->x-fiberp->x)+(cellp->y-fiberp->y)*(cellp->y-fiberp->y));
				targetx=fiberp->x+fiberp->radius*nvec_x;
				targety=fiberp->y+fiberp->radius*nvec_y;
			}
			
			if((cellp->x < (targetx+0.5*dx))&&(cellp->x > (targetx-0.5*dx))){
				if((cellp->y < (targety+0.5*dy))&&(cellp->y > (targety-0.5*dy))){
					defflag=1;
				}
			}
			
			if(defflag==1){
				if(fabs(nvec_x)>fabs(nvec_y)){
					upx=-nvec_y*(cellp->y+0.5*dy-targety)/nvec_x+targetx;
					downx=-nvec_y*(cellp->y-0.5*dy-targety)/nvec_x+targetx;
					
					if(((upx-cellp->x+0.5*dx)>=0.0)&&((cellp->x+0.5*dx-upx)>=0.0)){
						if((downx-cellp->x+0.5*dx)>=0.0){
							if((cellp->x+0.5*dx-downx)>=0.0){
								cellp->vof=0.5*dy*((upx-cellp->x+0.5*dx)+(downx-cellp->x+0.5*dx))/dx/dy;
							}
							else{
								cellp->vof=(dx*dy-0.5*dy/(downx-upx)*(cellp->x+0.5*dx-upx)*(cellp->x+0.5*dx-upx))/dx/dy;
							}
						}
						else {
							cellp->vof=0.5*dy/(upx-downx)*(upx-cellp->x+0.5*dx)*(upx-cellp->x+0.5*dx)/dx/dy;
						}
					}
					else {
						if((upx-cellp->x+0.5*dx)>=0.0){
							if((cellp->x+0.5*dx-upx)>=0.0){
								cellp->vof=0.5*dy*((downx-cellp->x+0.5*dx)+(upx-cellp->x+0.5*dx))/dx/dy;
							}
							else{
								cellp->vof=(dx*dy-0.5*dy/(upx-downx)*(cellp->x+0.5*dx-downx)*(cellp->x+0.5*dx-downx))/dx/dy;
							}
						}
						else {
							cellp->vof=0.5*dy/(downx-upx)*(downx-cellp->x+0.5*dx)*(downx-cellp->x+0.5*dx)/dx/dy;
						}
					}
					
					if(cellp->vof>1.0)cellp->vof=1.0;
					if(cellp->vof<0.0)cellp->vof=0.0;
					if(nvec_x<0.0) cellp->vof=1.0-cellp->vof;
				}
				else{
					upy=-nvec_x*(cellp->x+0.5*dx-targetx)/nvec_y+targety;
					downy=-nvec_x*(cellp->x-0.5*dx-targetx)/nvec_y+targety;
					
					if(((upy-cellp->y+0.5*dy)>=0.0)&&((cellp->y+0.5*dy-upy)>=0.0)){
						if((downy-cellp->y+0.5*dy)>=0.0){
							if((cellp->y+0.5*dy-downy)>=0.0){
								cellp->vof=0.5*dx*((upy-cellp->y+0.5*dy)+(downy-cellp->y+0.5*dy))/dy/dx;
							}
							else{
								cellp->vof=(dy*dx-0.5*dx/(downy-upy)*(cellp->y+0.5*dy-upy)*(cellp->y+0.5*dy-upy))/dy/dx;
							}
						}
						else {
							cellp->vof=0.5*dx/(upy-downy)*(upy-cellp->y+0.5*dy)*(upy-cellp->y+0.5*dy)/dy/dx;
						}
					}
					else {
						if((upy-cellp->y+0.5*dy)>=0.0){
							if((cellp->y+0.5*dy-upy)>=0.0){
								cellp->vof=0.5*dx*((downy-cellp->y+0.5*dy)+(upy-cellp->y+0.5*dy))/dy/dx;
							}
							else{
								cellp->vof=(dy*dx-0.5*dx/(upy-downy)*(cellp->y+0.5*dy-downy)*(cellp->y+0.5*dy-downy))/dy/dx;
							}
						}
						else {
							cellp->vof=0.5*dx/(downy-upy)*(downy-cellp->y+0.5*dy)*(downy-cellp->y+0.5*dy)/dy/dx;
						}
					}
					
					if(cellp->vof>1.0)cellp->vof=1.0;
					if(cellp->vof<0.0)cellp->vof=0.0;
					if(nvec_y<0.0) cellp->vof=1.0-cellp->vof;
				}
			}
			
			
			if((cellp->ls_x >-1.0*dx)&&(cellp->ls_x<1.0*dx)){
				nvec_x=(cellp->x-0.5*dx-fiberp->x)/sqrt((cellp->x-0.5*dx-fiberp->x)*(cellp->x-0.5*dx-fiberp->x)+(cellp->y-fiberp->y)*(cellp->y-fiberp->y));
				nvec_y=(cellp->y-fiberp->y)/sqrt((cellp->x-0.5*dx-fiberp->x)*(cellp->x-0.5*dx-fiberp->x)+(cellp->y-fiberp->y)*(cellp->y-fiberp->y));
				targetx=fiberp->x+fiberp->radius*nvec_x;
				targety=fiberp->y+fiberp->radius*nvec_y;
			}
			
			
			if((cellp->x < (targetx+dx))&&(cellp->x > targetx)){
				if((cellp->y < (targety+0.5*dy))&&(cellp->y > (targety-0.5*dy))){
					defflag_x=1;
				}
			}
			
			if(defflag_x==1){
				if(fabs(nvec_x)>fabs(nvec_y)){
					upx=-nvec_y*(cellp->y+0.5*dy-targety)/nvec_x+targetx;
					downx=-nvec_y*(cellp->y-0.5*dy-targety)/nvec_x+targetx;
					
					if(((upx-(cellp->x-0.5*dx)+0.5*dx)>=0.0)&&(((cellp->x-0.5*dx)+0.5*dx-upx)>=0.0)){
						if((downx-(cellp->x-0.5*dx)+0.5*dx)>=0.0){
							if(((cellp->x-0.5*dx)+0.5*dx-downx)>=0.0){
								cellp->vof_x=0.5*dy*((upx-(cellp->x-0.5*dx)+0.5*dx)+(downx-(cellp->x-0.5*dx)+0.5*dx))/dx/dy;
							}
							else{
								cellp->vof_x=(dx*dy-0.5*dy/(downx-upx)*((cellp->x-0.5*dx)+0.5*dx-upx)*((cellp->x-0.5*dx)+0.5*dx-upx))/dx/dy;
							}
						}
						else {
							cellp->vof_x=0.5*dy/(upx-downx)*(upx-(cellp->x-0.5*dx)+0.5*dx)*(upx-(cellp->x-0.5*dx)+0.5*dx)/dx/dy;
						}
					}
					else {
						if((upx-(cellp->x-0.5*dx)+0.5*dx)>=0.0){
							if(((cellp->x-0.5*dx)+0.5*dx-upx)>=0.0){
								cellp->vof_x=0.5*dy*((downx-(cellp->x-0.5*dx)+0.5*dx)+(upx-(cellp->x-0.5*dx)+0.5*dx))/dx/dy;
							}
							else{
								cellp->vof_x=(dx*dy-0.5*dy/(upx-downx)*((cellp->x-0.5*dx)+0.5*dx-downx)*((cellp->x-0.5*dx)+0.5*dx-downx))/dx/dy;
							}
						}
						else {
							cellp->vof_x=0.5*dy/(downx-upx)*(downx-(cellp->x-0.5*dx)+0.5*dx)*(downx-(cellp->x-0.5*dx)+0.5*dx)/dx/dy;
						}
					}
					
					if(cellp->vof>1.0)cellp->vof_x=1.0;
					if(cellp->vof<0.0)cellp->vof_x=0.0;
					if(nvec_x<0.0) cellp->vof_x=1.0-cellp->vof_x;
				}
				else{
					upy=-nvec_x*((cellp->x-0.5*dx)+0.5*dx-targetx)/nvec_y+targety;
					downy=-nvec_x*((cellp->x-0.5*dx)-0.5*dx-targetx)/nvec_y+targety;
					
					if(((upy-cellp->y+0.5*dy)>=0.0)&&((cellp->y+0.5*dy-upy)>=0.0)){
						if((downy-cellp->y+0.5*dy)>=0.0){
							if((cellp->y+0.5*dy-downy)>=0.0){
								cellp->vof_x=0.5*dx*((upy-cellp->y+0.5*dy)+(downy-cellp->y+0.5*dy))/dy/dx;
							}
							else{
								cellp->vof_x=(dy*dx-0.5*dx/(downy-upy)*(cellp->y+0.5*dy-upy)*(cellp->y+0.5*dy-upy))/dy/dx;
							}
						}
						else {
							cellp->vof_x=0.5*dx/(upy-downy)*(upy-cellp->y+0.5*dy)*(upy-cellp->y+0.5*dy)/dy/dx;
						}
					}
					else {
						if((upy-cellp->y+0.5*dy)>=0.0){
							if((cellp->y+0.5*dy-upy)>=0.0){
								cellp->vof_x=0.5*dx*((downy-cellp->y+0.5*dy)+(upy-cellp->y+0.5*dy))/dy/dx;
							}
							else{
								cellp->vof_x=(dy*dx-0.5*dx/(upy-downy)*(cellp->y+0.5*dy-downy)*(cellp->y+0.5*dy-downy))/dy/dx;
							}
						}
						else {
							cellp->vof_x=0.5*dx/(downy-upy)*(downy-cellp->y+0.5*dy)*(downy-cellp->y+0.5*dy)/dy/dx;
						}
					}
					
					if(cellp->vof>1.0)cellp->vof_x=1.0;
					if(cellp->vof<0.0)cellp->vof_x=0.0;
					if(nvec_y<0.0) cellp->vof_x=1.0-cellp->vof_x;
				}
			}
			
			
			if((cellp->ls_y >-1.0*dx)&&(cellp->ls_y<1.0*dx)){
				nvec_x=(cellp->x-fiberp->x)/sqrt((cellp->x-fiberp->x)*(cellp->x-fiberp->x)+(cellp->y-0.5*dy-fiberp->y)*(cellp->y-0.5*dy-fiberp->y));
				nvec_y=(cellp->y-0.5*dy-fiberp->y)/sqrt((cellp->x-fiberp->x)*(cellp->x-fiberp->x)+(cellp->y-0.5*dy-fiberp->y)*(cellp->y-0.5*dy-fiberp->y));
				targetx=fiberp->x+fiberp->radius*nvec_x;
				targety=fiberp->y+fiberp->radius*nvec_y;
			}
			
			if((cellp->x < (targetx+0.5*dx))&&(cellp->x > (targetx-0.5*dx))){
				if((cellp->y < (targety+dy))&&(cellp->y > targety)){
					defflag_y=1;
				}
			}
			
			if(defflag_y==1){
				if(fabs(nvec_x)>fabs(nvec_y)){
					upx=-nvec_y*((cellp->y-0.5*dy)+0.5*dy-targety)/nvec_x+targetx;
					downx=-nvec_y*((cellp->y-0.5*dy)-0.5*dy-targety)/nvec_x+targetx;
					
					if(((upx-cellp->x+0.5*dx)>=0.0)&&((cellp->x+0.5*dx-upx)>=0.0)){
						if((downx-cellp->x+0.5*dx)>=0.0){
							if((cellp->x+0.5*dx-downx)>=0.0){
								cellp->vof_y=0.5*dy*((upx-cellp->x+0.5*dx)+(downx-cellp->x+0.5*dx))/dx/dy;
							}
							else{
								cellp->vof_y=(dx*dy-0.5*dy/(downx-upx)*(cellp->x+0.5*dx-upx)*(cellp->x+0.5*dx-upx))/dx/dy;
							}
						}
						else {
							cellp->vof_y=0.5*dy/(upx-downx)*(upx-cellp->x+0.5*dx)*(upx-cellp->x+0.5*dx)/dx/dy;
						}
					}
					else {
						if((upx-cellp->x+0.5*dx)>=0.0){
							if((cellp->x+0.5*dx-upx)>=0.0){
								cellp->vof_y=0.5*dy*((downx-cellp->x+0.5*dx)+(upx-cellp->x+0.5*dx))/dx/dy;
							}
							else{
								cellp->vof_y=(dx*dy-0.5*dy/(upx-downx)*(cellp->x+0.5*dx-downx)*(cellp->x+0.5*dx-downx))/dx/dy;
							}
						}
						else {
							cellp->vof_y=0.5*dy/(downx-upx)*(downx-cellp->x+0.5*dx)*(downx-cellp->x+0.5*dx)/dx/dy;
						}
					}
					
					if(cellp->vof>1.0)cellp->vof_y=1.0;
					if(cellp->vof<0.0)cellp->vof_y=0.0;
					if(nvec_x<0.0) cellp->vof_y=1.0-cellp->vof_y;
				}
				else{
					upy=-nvec_x*(cellp->x+0.5*dx-targetx)/nvec_y+targety;
					downy=-nvec_x*(cellp->x-0.5*dx-targetx)/nvec_y+targety;
					
					if(((upy-(cellp->y-0.5*dy)+0.5*dy)>=0.0)&&(((cellp->y-0.5*dy)+0.5*dy-upy)>=0.0)){
						if((downy-(cellp->y-0.5*dy)+0.5*dy)>=0.0){
							if(((cellp->y-0.5*dy)+0.5*dy-downy)>=0.0){
								cellp->vof_y=0.5*dx*((upy-(cellp->y-0.5*dy)+0.5*dy)+(downy-(cellp->y-0.5*dy)+0.5*dy))/dy/dx;
							}
							else{
								cellp->vof_y=(dy*dx-0.5*dx/(downy-upy)*((cellp->y-0.5*dy)+0.5*dy-upy)*((cellp->y-0.5*dy)+0.5*dy-upy))/dy/dx;
							}
						}
						else {
							cellp->vof_y=0.5*dx/(upy-downy)*(upy-(cellp->y-0.5*dy)+0.5*dy)*(upy-(cellp->y-0.5*dy)+0.5*dy)/dy/dx;
						}
					}
					else {
						if((upy-(cellp->y-0.5*dy)+0.5*dy)>=0.0){
							if(((cellp->y-0.5*dy)+0.5*dy-upy)>=0.0){
								cellp->vof_y=0.5*dx*((downy-(cellp->y-0.5*dy)+0.5*dy)+(upy-(cellp->y-0.5*dy)+0.5*dy))/dy/dx;
							}
							else{
								cellp->vof_y=(dy*dx-0.5*dx/(upy-downy)*((cellp->y-0.5*dy)+0.5*dy-downy)*((cellp->y-0.5*dy)+0.5*dy-downy))/dy/dx;
							}
						}
						else {
							cellp->vof_y=0.5*dx/(downy-upy)*(downy-(cellp->y-0.5*dy)+0.5*dy)*(downy-(cellp->y-0.5*dy)+0.5*dy)/dy/dx;
						}
					}
					
					if(cellp->vof>1.0)cellp->vof_y=1.0;
					if(cellp->vof<0.0)cellp->vof_y=0.0;
					if(nvec_y<0.0) cellp->vof_y=1.0-cellp->vof_y;
				}
			}
		}
	}
	
	
	
	for(int i=1; i<(XCELL_NUM+5); i++){
		for(int j=1; j<(YCELL_NUM+5); j++){
			index0   = j     + i*(YCELL_NUM+6);
			
			cellp    = &CELL[index0];
			
			_flag_type* flagp    = &FLAG[index0];
			
			if(cellp->vof>0.99999) flagp->c=0;
			else if(cellp->vof<0.00001) flagp->c=1;
			else flagp->c=2;
			
			if(cellp->vof_x>0.99999) flagp->x=0;
			else if(cellp->vof_x<0.00001) flagp->x=1;
			else flagp->x=2;
			
			if(cellp->vof_y>0.99999) flagp->y=0;
			else if(cellp->vof_y<0.00001) flagp->y=1;
			else flagp->y=2;
		}
	}
	
	for(int i=3; i<(XCELL_NUM+3); i++){
		for(int j=3; j<(YCELL_NUM+3); j++){
			index0     = j  + i*(YCELL_NUM+6);
			index_xm   = j + (i-1)*(YCELL_NUM+6);
			index_xp  =  j     + (i+1)*(YCELL_NUM+6);
			index_yp  =  (j+1) + i*(YCELL_NUM+6);
			index_ym  =  (j-1) + i*(YCELL_NUM+6);
			index_xpp =  j     + (i+2)*(YCELL_NUM+6);
			index_ypp =  (j+2) + i*(YCELL_NUM+6);
			
			_flag_type* flagp    = &FLAG[index0];
			_flag_type* flagp_xm = &FLAG[index_xm];
			_flag_type* flagp_ym = &FLAG[index_ym];
			_flag_type* flagp_xp = &FLAG[index_xp];
			_flag_type* flagp_yp = &FLAG[index_yp];
			_flag_type* flagp_xpp = &FLAG[index_xpp];
			_flag_type* flagp_ypp = &FLAG[index_ypp];
			
			if(flagp->c==0){
				if((flagp_xp->c==2)||(flagp_xm->c==2)||(flagp_yp->c==2)
					||(flagp_ym->c==2)||(flagp_xm->x==2)||(flagp_ym->y==2)
					||(flagp_xpp->x==2)||(flagp_ypp->y==2)){
						flagp->c=4;
				}
			}
		}
	}
	
/*	for(int i=3; i<(XCELL_NUM+3); i++){
		for(int j=3; j<(YCELL_NUM+3); j++){
			index0     = j  + i*(YCELL_NUM+6);
			index_xm   = j + (i-1)*(YCELL_NUM+6);
			index_xp  =  j     + (i+1)*(YCELL_NUM+6);
			index_yp  =  (j+1) + i*(YCELL_NUM+6);
			index_ym  =  (j-1) + i*(YCELL_NUM+6);
			index_xpp =  j     + (i+2)*(YCELL_NUM+6);
			index_ypp =  (j+2) + i*(YCELL_NUM+6);
			
			_flag_type* flagp    = &FLAG[index0];
			_flag_type* flagp_xm = &FLAG[index_xm];
			_flag_type* flagp_ym = &FLAG[index_ym];
			_flag_type* flagp_xp = &FLAG[index_xp];
			_flag_type* flagp_yp = &FLAG[index_yp];
			_flag_type* flagp_xpp = &FLAG[index_xpp];
			_flag_type* flagp_ypp = &FLAG[index_ypp];
			
			if((flagp->c==2)||(flagp->c==4)){
				if((flagp->x%4)&&(flagp->y%4)&&(flagp_xp->x%4)&&(flagp_yp->y%4)){
						flagp->c=6;
				}
			}
		}
	}*/
}

#ifdef USEMPI
int main(int argc, char* argv[]){
	int restartout;
	int chkflag;
	int my_rank;	/* MPI rank */
	
	//-------- Parallel computation ------//
	
	TCELL_NUM=400000;
	nout=2000;
	restartout=10000;
	chkflag=0;
	n=0;
	

	MPI_Init(&argc, &argv);			// Call this exactly once before calling any other MPI functions
	
	MPI_Comm_size(MPI_COMM_WORLD, &process_num);		// Get the number of processes in the communicator.
	
	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);		// Each process in the communicator obtains its own rank
	initial(my_rank);
	if(my_rank==0) filewrite();
	
	if(my_rank==0){
		for(; n<(TCELL_NUM+1); n++) {
			hontai(chkflag,my_rank);
//			printf("n = %d \n",n);
			if(n==1) filewrite();
			if(n==100) filewrite();
//			discount();
			if(n%nout==0) {
				filewrite();
				
			}
			if(n%restartout==0){
				restartwrite();
			}
		}
	}
	else {
		for(; n<(TCELL_NUM+1); n++){
			SOR(chkflag,my_rank);
		}
	}
	

	MPI_Finalize();
	exit(EXIT_SUCCESS);
	
	delete [] CELL;
	delete [] VOLD;
	delete [] TMPMU;
	delete [] FLAG;
	return 0;
}
#else	// Else for USEMPI
int main(void){
	int restartout;
	int chkflag;
	TCELL_NUM=300000;
	nout=2000;
	restartout=5000;
	chkflag=0;
	n=0;
	
	initial();
	filewrite();
	
	for(; n<(TCELL_NUM+1); n++) {
		hontai(chkflag);
		if(n==1) filewrite();
		if(n%nout==0) {
			filewrite();
//			discount();
		}
		if(n%restartout==0){
			restartwrite();
		}
	}
	
	
	delete [] CELL;
	delete [] VOLD;
	delete [] TMPMU;
	delete [] FLAG;
	delete [] ddrho;
	delete [] fiberfx;
	delete [] fiberfy;
	delete [] fibern;
	return 0;
}
#endif

			//-------------- End of main function ----------------//

			//-------------- Initialization -----------------------//
#ifdef USEMPI
void initial(int my_rank){
#else
void initial(void){
#endif
	
	CELL = new _cell_type[CELL_SIZE];
	memset(CELL,0,CELL_SIZE*sizeof(_cell_type));
	VOLD = new _vold_type[CELL_SIZE];
	memset(VOLD,0,CELL_SIZE*sizeof(_vold_type));
	TMPMU = new _mu_type[CELL_SIZE];
	memset(TMPMU,0,CELL_SIZE*sizeof(_mu_type));
	Fx = new double [CELL_SIZE];
	memset(Fx,0,CELL_SIZE*sizeof(double));
	Fy = new double [CELL_SIZE];
	memset(Fy,0,CELL_SIZE*sizeof(double));
	phid2 = new double [CELL_SIZE];
	memset(phid2,0,CELL_SIZE*sizeof(double));
	ddrho = new double [CELL_SIZE];
	memset(ddrho,0,CELL_SIZE*sizeof(double));
	fiberfx = new double [fibernum];
	memset(fiberfx,0,fibernum*sizeof(double));
	fiberfy = new double [fibernum];
	memset(fiberfy,0,fibernum*sizeof(double));
	fibern = new double [fibernum];
	memset(fibern,0,fibernum*sizeof(double));
	
#ifndef WALLX
		fibernum*=2;
#endif

#ifndef WALLY
		fibernum*=2;
#endif
	
	FIBER = new _fiber_type [fibernum];
	memset(FIBER,0,fibernum*sizeof(_fiber_type));
	wetphi = new double [fibernum];
	memset(wetphi,0,fibernum*sizeof(double));
	
#ifdef USEMPI
	mail2 = new double [YCELL_NUM];
	mail1 = new double [YCELL_NUM];
	memset(mail2,0,YCELL_NUM*sizeof(double));
	memset(mail1,0,YCELL_NUM*sizeof(double));
#endif

	FLAG = new _flag_type[CELL_SIZE];
	memset(FLAG,0,CELL_SIZE*sizeof(_flag_type));

	
	_fiber_type* fiberp=&FIBER[0];
	
	fiberp->x=48.0;
	fiberp->y=30.0;
	fiberp->radius = 20.0;
	
//	fiberp=&FIBER[1];
	
//	fiberp->x=80.0;
//	fiberp->y=48.0;
//	fiberp->radius = 15.0;
	
	wetphiwall=gammas/kappa2;
	
	wetphi[0]=gammas/kappa2;
		
//	wetphi[1]=gammas/kappa2;
	
	for(int f=FIBER_NUM; f<fibernum; f++){
		wetphi[f]=wetphi[f-FIBER_NUM];
	}
	
	for(int i=0; i<(XCELL_NUM+6); i++){
		for(int j=0; j<(YCELL_NUM+6); j++){
			int index = j + i*(YCELL_NUM+6);
			
			_cell_type* cellp= &CELL[index];
			
			cellp->x = ((double)i-2.5)*dx;
			cellp->y = ((double)j-2.5)*dy;
			
//			cellp->v_x = Uinf;				//IB1
		}
	}
	
	
	for(int i=0; i<(XCELL_NUM+6); i++){
		for(int j=0; j<(YCELL_NUM+6); j++){
			int index = j + i*(YCELL_NUM+6);
			
			_cell_type* cellp= &CELL[index];
			
			cellp->x = ((double)i-2.5)*dx;
			cellp->y = ((double)j-2.5)*dy;
			
//			cellp->v_x = Uinf;				//IB1
		}
	}

	flagdefine();
	
	vofbnd();
	vofbnd();
	
/*	for(int i=0; i<(XCELL_NUM+6); i++) {
		for(int j=0; j<(YCELL_NUM+6); j++) {
			index0 = j + i*(YCELL_NUM+6);
			_cell_type* cellp= &CELL[index0];
			_flag_type* flagp= &FLAG[index0];
			
			cellp->phi = PHImax;
			
			if(cellp->vof>0.99) flagp->c=1;
			else flagp->c=1;
			
		}
	}
	
	for(int i=1; i<(XCELL_NUM+6); i++) {
		for(int j=1; j<(YCELL_NUM+6); j++) {
			index0 = j + i*(YCELL_NUM+6);
			index_xm = j + (i-1)*(YCELL_NUM+6);
			index_ym = j-1 + i*(YCELL_NUM+6);
			_cell_type* cellp= &CELL[index0];
			_flag_type* flagp= &FLAG[index0];
			_flag_type* flagp_xm= &FLAG[index_xm];
			_flag_type* flagp_ym= &FLAG[index_ym];
			
			if((flagp->c==1)&&(flagp_xm->c==1)) flagp->x=1;
			else flagp->x=1;
			if((flagp->c==1)&&(flagp_ym->x==1)) flagp->y=1;
			else flagp->y=1;
		}
	}
*/	
	
	flagbnd();
	flagbnd();
	
	for(int i=0; i<(XCELL_NUM+6); i++) {
		for(int j=0; j<7; j++) {
			index0 = j + i*(YCELL_NUM+6);
			_cell_type* cellp= &CELL[index0];
			
//			cellp->phi = PHImax;
			
		}
	}
		
		
	
	for(int i=0; i<(XCELL_NUM+6); i++) {
		for(int j=0; j<(YCELL_NUM+6); j++) {
			double distance;
			index0 = j + i*(YCELL_NUM+6);
			cellp= &CELL[index0];
			
//			if((B*cellp->x-A*cellp->y+D1)>0.0) cellp->phi=PHImax;
//			else cellp->phi=PHImin;
				
			
			distance = sqrt((cellp->x - XVALUE)*(cellp->x - XVALUE)+(cellp->y - YVALUE)*(cellp->y - YVALUE));
			
			if(distance<=(LRADIUS*(sqrt((1.0-DEGREE/180.0)+(sin(M_PI*DEGREE/180.0))/M_PI)))) cellp->phi = PHImax;
			else cellp->phi=PHImin;
			
//			if((cellp->y-0.25*cellp->x-11.0)>0.0){
//				distance = (cellp->x - 48.0)*(cellp->x - 48.0)+(cellp->y - 60.0)*(cellp->y - 60.0);
//				if((distance<diameter2)&&(cellp->vof<0.999)){
//					cellp->phi = PHImax;
//				}
	
//				distance = (cellp->x - ((double)XCELL_NUM*dx+12.0))*(cellp->x - ((double)XCELL_NUM*dx+12.0))+(cellp->y - 51.0)*(cellp->y - 51.0);
	//			distance = (cellp->x - 32.0)*(cellp->x - 32.0)+(cellp->y - 50.0)*(cellp->y - 50.0);
//				if((distance<diameter2)&&(cellp->vof<0.999)){
//					cellp->phi = PHImax;
//				}
	//				else cellp->phi = PHImin;
//			}
		}
	}
	
	if(n!=0) {
		if(my_rank==0) restartread();
	}
	phibnd();
	phibnd();

//	if(my_rank==0) filewrite();
	
	for(int i=3; i<(XCELL_NUM+3); i++) {
		for(int j=3; j<(YCELL_NUM+3); j++) {
			indexsetdefine(i,j);
			_flag_type* flagp=&FLAG[index0];
			
			if((flagp->c==2)||(flagp->c==4)||(flagp->c==6)||((cellp->vof>0.000000001)&&(cellp->vof<0.999999999999))) cellp->phi=interphi(index0,0);
			else if(flagp->c==0) cellp->phi=PHImax;
		}
	}
	
	phibnd();
	phibnd();
	
	for(int i=3; i<(XCELL_NUM+4); i++) {
		for(int j=3; j<(YCELL_NUM+4); j++) {
			indexsetdefine(i,j);

			cellp->phi_x = (9.0*(cellp->phi+cellp_xm->phi)-cellp_xp->phi-cellp_xmm->phi)/16.0/dx;
			cellp->phi_y = (9.0*(cellp->phi+cellp_ym->phi)-cellp_yp->phi-cellp_ymm->phi)/16.0/dy;
		}
	}
	
	phikyokaibnd();
	phikyokaibnd();
	
	rhodefine();

	for(int i=0;i<1000;i++) {
		CH_eq();
	}
}
			//-------------- End of Initialization -----------------//
	
			//-------------- The main of Calculation ---------------------//
#ifdef USEMPI
void hontai(int chkflag, int my_rank){
	for(int i=0; i<(XCELL_NUM+6); i++) {
		for(int j=0; j<(YCELL_NUM+6); j++) {
				index0 = j + i*(YCELL_NUM+6);
				
				_vold_type* voldp = &VOLD[index0];
				cellp = &CELL[index0];
				
				voldp->x = cellp->v_x;
				voldp->y = cellp->v_y;
		}
	}
//	if(my_rank==0)printf("ok");
	projection(chkflag,my_rank);
	
	CH_eq();
	
	flagdefine();
	
	flagbnd();
	flagbnd();
}
	
#else		//Else for USEMPI
	
void hontai(int chkflag){
	for(int i=0; i<(XCELL_NUM+6); i++) {
		for(int j=0; j<(YCELL_NUM+6); j++) {
				index0 = j + i*(YCELL_NUM+6);
				
				_vold_type* voldp = &VOLD[index0];
				cellp = &CELL[index0];
				
				voldp->x = cellp->v_x;							// At this point, cellp->v_x is still the same as before

				voldp->y = cellp->v_y;
			
		}
	}
	
	projection(chkflag);
	
	CH_eq();
	
	flagdefine();
	flagbnd();
	flagbnd();
}	
#endif
			//-------------- End of The main Calculation ---------------//
	
	
			//-------------- projection -------------------//
#ifdef USEMPI
void projection(int chkflag,int my_rank){
	
	nseq();
	
	vhypobnd();
	vhypobnd();


	SOR(chkflag,my_rank);
	
	koshin();
	
	kyokai();
	kyokai();
}
#else //Else for USEMPI
void projection(int chkflag){
	nseq();
	
	vhypobnd();
	vhypobnd();

	SOR(chkflag);
	
	koshin();

	kyokai();
	kyokai();
}
#endif //USEMPI
			//-------------- End of projection -------------//
	
	
			//-------------- Derivation of tentative flow velocity ---------------//
void nseq(void){
	double vx_average,vy_average,vxvxx,vyvxy,vyvyy,vxvyx;
	double dxi, dyi, iryu, kaku,kaimen;
	double muvxdxp,muvxdx,muvydyp,muvydy;
	
	ddrhodefine();
	
	ddrhobnd();
	rhodefine();
	mudefine();

	dxi = 1.0/dx;
	dyi = 1.0/dy;
	
#ifndef WALLX	
	for(int i=3/*4*/; i<(XCELL_NUM+3); i++){					// For periodic boundaries, it's 3

#else
	for(int i=3; i<(XCELL_NUM+3); i++){					
#endif
		for(int j=3; j<(YCELL_NUM+3); j++){
			indexnsdefine(i,j);
				
				vy_average = (cellp_xm_yp->v_y + cellp_yp->v_y + cellp_xm->v_y + cellp->v_y)*0.25;
				
				vxvxx = cellp->v_x * ( -cellp_xpp->v_x + 8.0*(cellp_xp->v_x - cellp_xm->v_x) + cellp_xmm->v_x)/12.0*dxi
						+fabs(cellp->v_x)*(cellp_xpp->v_x - 4.0*cellp_xp->v_x + 6.0*cellp->v_x - 4.0*cellp_xm->v_x + cellp_xmm->v_x)*0.25*dxi;
				vyvxy = vy_average * ( -cellp_ypp->v_x + 8.0*(cellp_yp->v_x - cellp_ym->v_x) + cellp_ymm->v_x)/12.0*dyi
						+fabs(vy_average)*(cellp_ypp->v_x - 4.0*cellp_yp->v_x + 6.0*cellp->v_x - 4.0*cellp_ym->v_x + cellp_ymm->v_x)*0.25*dyi;
				muvxdxp = mup->c    * (cellp_xp->v_x - cellp->v_x   )*dxi;
				muvxdx  = mup_xm->c * (cellp->v_x    - cellp_xm->v_x)*dxi;
				
				iryu=-(vxvxx+vyvxy);
				kaku=(2.0*dxi*(muvxdxp-muvxdx)
					 +dyi*(mup_yp->xy *(dxi*(cellp_yp->v_y - cellp_xm_yp->v_y)+dyi*(cellp_yp->v_x - cellp->v_x))
				-mup->xy    *(dxi*(cellp->v_y    - cellp_xm->v_y)   +dyi*(cellp->v_x    - cellp_ym->v_x)))) / (cellp->rho_x*(1.0-cellp->vof_x)+RHOS*cellp->vof_x);
//					atsu= -(-cellp_xp->press+27.0*(cellp->press-cellp_xm->press)+cellp_xmm->press)*dxi/cellp->rho_x/24.0;
				
					kaimen=(1.0-cellp->vof_x)*kappa1*dxi*(-ddrho[index_xp]+27.0*(ddrho[index0]-ddrho[index_xm])+ddrho[index_xmm])/24.0*cellp->rho_x/(cellp->rho_x*(1.0-cellp->vof_x)+RHOS*cellp->vof_x);

	if(n==0) cellp->v_xhypo = cellp->v_x + dt*(iryu + kaku + kaimen);	//1���
	else cellp->v_xhypo = cellp->v_x + 0.5*dt*(3.0*(iryu + kaku + kaimen)-cellp->hx);
				cellp->hx=iryu + kaku + kaimen;
				//				printf("%f %f %f %f %f %f \n",cellp->x,cellp->y,cellp->v_xhypo,iryu,kaku,kaimen);//debug
		}
	}

		

	for(int i=3; i<(XCELL_NUM+4); i++){
#ifndef WALLY
		for(int j=3/*4*/; j<(YCELL_NUM+3); j++){					// For periodic boundaries, it's 3
#else
		for(int j=3; j<(YCELL_NUM+3); j++){					
#endif
			indexnsdefine(i,j);
			
				vx_average = (cellp_xp_ym->v_x + cellp_ym->v_x + cellp_xp->v_x + cellp->v_x)*0.25;
				
				vyvyy = cellp->v_y *( -cellp_ypp->v_y + 8.0*(cellp_yp->v_y - cellp_ym->v_y) + cellp_ymm->v_y)/12.0*dyi
							+fabs(cellp->v_y)*(cellp_ypp->v_y - 4.0*cellp_yp->v_y + 6.0*cellp->v_y - 4.0*cellp_ym->v_y + cellp_ymm->v_y)*0.25*dyi;
				vxvyx = vx_average *( -cellp_xpp->v_y + 8.0*(cellp_xp->v_y - cellp_xm->v_y) + cellp_xmm->v_y)/12.0*dxi
							+fabs(vx_average)*(cellp_xpp->v_y - 4.0*cellp_xp->v_y + 6.0*cellp->v_y - 4.0*cellp_xm->v_y + cellp_xmm->v_y)*0.25*dxi;
				muvydyp = mup->c    *(cellp_yp->v_y - cellp->v_y)   *dyi;
				muvydy  = mup_ym->c *(cellp->v_y    - cellp_ym->v_y)*dyi;
				
				iryu=-(vyvyy+vxvyx);
				kaku=(2.0*dyi*(muvydyp-muvydy)
					  +dxi*(mup_xp->xy *(dyi*(cellp_xp->v_x - cellp_xp_ym->v_x)+ dxi*(cellp_xp->v_y - cellp->v_y))
						   -mup->xy    *(dyi*(cellp->v_x    - cellp_ym->v_x)   + dxi*(cellp->v_y    - cellp_xm->v_y)))) / (cellp->rho_y*(1.0-cellp->vof_y)+RHOS*cellp->vof_y);

			kaimen=(1.0-cellp->vof_y)*kappa1*dyi*(-ddrho[index_yp]+27.0*(ddrho[index0]-ddrho[index_ym])+ddrho[index_ymm])/24.0*cellp->rho_y/(cellp->rho_y*(1.0-cellp->vof_y)+RHOS*cellp->vof_y);

	if(n==0) cellp->v_yhypo = cellp->v_y + dt*(iryu + kaku + kaimen);	// First time
	else			cellp->v_yhypo = cellp->v_y + dt*0.5*(3.0*(iryu + kaku + kaimen)-cellp->hy);
				cellp->hy=iryu + kaku + kaimen;
		}
	}
		
		
	for(int i=0; i<(XCELL_NUM+6); i++) {
		for(int j=0; j<(YCELL_NUM+6); j++) {
			index0 = j + i*(YCELL_NUM+6);
			
			cellp->rho_x = cellp->rho_x*(1.0-cellp->vof_x)+RHOS*cellp->vof_x;
			cellp->rho_y = cellp->rho_y*(1.0-cellp->vof_y)+RHOS*cellp->vof_y;
		}
	}
}
			//-------------- End of deriving tentative flow velocity ---------//

		
		
		
		
		
		
		
			//-------------- koshin ----------------------//
void koshin(void){
	double pressaverage = 0.0;
	int pressnum=0;
	
/*	for(int i=3; i<(XCELL_NUM+3); i++) {
		for(int j=3; j<(YCELL_NUM+3); j++) {
			index0 = j + i*(YCELL_NUM+6);
			
			cellp = &CELL[index0];
			
			pressaverage+=cellp->press;
			pressnum++;
		}
	}
	
	pressaverage = pressaverage / (double)pressnum;
	
	for(int i=3; i<(XCELL_NUM+3); i++) {
		for(int j=3; j<(YCELL_NUM+3); j++) {
			index0 = j + i*(YCELL_NUM+6);
			
			cellp = &CELL[index0];
			
			cellp->press = cellp->press - pressaverage;
		}
	}
*/	
	pressbnd();
	
	for(int i=0; i<FIBER_NUM; i++){
		_fiber_type* fiberp =& FIBER[i];
		
		fiberfx[i] = 0.0;
		fiberfy[i] = 0.0;
		fibern[i] = 0.0;
		
//		solidlagpoint(fiberp,i);
	}
	
	for(int i=3;i<(XCELL_NUM+3);i++){
		for(int j=3;j<(YCELL_NUM+3);j++){
			index0    = j + i*(YCELL_NUM+6);
			index_xm  = j + (i-1)*(YCELL_NUM+6);
			index_xp  = j + (i+1)*(YCELL_NUM+6);
			index_xmm = j + (i-2)*(YCELL_NUM+6);
			
			cellp= &CELL[index0];
			cellp_xp = &CELL[index_xp];
			cellp_xm = &CELL[index_xm];
			cellp_xmm = &CELL[index_xmm];
			
///				cellp->v_x=cellp->v_xhypo - 1.0*dt/cellp->rho_x*(-cellp_xp->press+27.0*(cellp->press-cellp_xm->press)+cellp_xmm->press)/dx/24.0;
				cellp->v_x=cellp->v_xhypo - 1.0*dt/cellp->rho_x*(cellp->press-cellp_xm->press)/dx;
		}
	}
	for(int i=3;i<(XCELL_NUM+3);i++){
		for(int j=3;j<(YCELL_NUM+3);j++){
			index0    = j     + i*(YCELL_NUM+6);
			index_yp  = (j+1) + i*(YCELL_NUM+6);
			index_ym  = (j-1) + i*(YCELL_NUM+6);
			index_ymm = (j-2) + i*(YCELL_NUM+6);
			
			cellp= &CELL[index0];
			cellp_yp = &CELL[index_yp];
			cellp_ym = &CELL[index_ym];
			cellp_ymm = &CELL[index_ymm];
	
//				cellp->v_y=cellp->v_yhypo - 1.0*dt/cellp->rho_y*(-cellp_yp->press+27.0*(cellp->press-cellp_ym->press)+cellp_ymm->press)/dy/24.0;
				cellp->v_y=cellp->v_yhypo - 1.0*dt/cellp->rho_y*(cellp->press-cellp_ym->press)/dy;
		}
	}
	
	for(int i=3; i<(XCELL_NUM+3); i++) {
		for(int j=3; j<(YCELL_NUM+3); j++) {
			index0 = j + i*(YCELL_NUM+6);
				
			_cell_type* cellp = &CELL[index0];
			_flag_type* flagp=&FLAG[index0];
			_fiber_type* fiberp;
			
			double distance;
			
			fiberp=&FIBER[flagp->fiber_x];
			
				cellp->v_x-=(cellp->v_x-0.0)*cellp->vof_x;
				
			fiberp=&FIBER[flagp->fiber_y];
			
				cellp->v_y-=(cellp->v_y-0.0)*cellp->vof_y;
		}
	}
	
	
	
	for(int i=0; i<FIBER_NUM; i++){
		_fiber_type* fiberp =& FIBER[i];
		
//		solidlagpointsecond(fiberp,i);
	}
	
//	printf("\n fiberfx=%f  fiberfy=%f  fiberfn=%f \n",fiberfx[0],fiberfy[0],fibern[0]);
	
/*	
	for(int i=0; i<(fibernum-1); i++){
		_fiber_type* fiberp=&FIBER[i];
		double disx,disy,dis;
		double sigma;
		
		for(int f=i+1; f<fibernum; f++){
			_fiber_type* fibertargetp=&FIBER[f];
			
			disx=fiberp->x-fibertargetp->x;
			disy=fiberp->y-fibertargetp->y;
			dis=sqrt(disx*disx+disy*disy);
			
//			printf(" dis=%f ",dis);
			sigma=fiberp->radius+fibertargetp->radius+SIGMAPLUS;
			
			if(dis<(sigma*1.122462)){
				fiberfx[i]+=dt*ljpower(sigma,dis)*disx/dis;
				fiberfx[f]-=dt*ljpower(sigma,dis)*disx/dis;
				fiberfy[i]+=dt*ljpower(sigma,dis)*disy/dis;
				fiberfy[f]-=dt*ljpower(sigma,dis)*disy/dis;
			}
		}
	}
*/		
	
	// Consider when an object crosses the boundary under periodic boundary conditions
/*	for(int f=fibernum-1; f>=FIBER_NUM; f--){
		fiberfx[f-FIBER_NUM]+=fiberfx[f];
		fiberfy[f-FIBER_NUM]+=fiberfy[f];
		fibern[f-FIBER_NUM]+=fibern[f];
	}
*/	
	
//	for(int i=0; i<FIBER_NUM; i++){
//		_fiber_type* fiberp=&FIBER[i];
//		
//		fiberp->v_x+=0.5*dt*(3.0*fiberfx[i]-fiberp->hx)/(RHOS*M_PI*fiberp->radius*fiberp->radius);
//		fiberp->hx=fiberfx[i];
//		if(i==0){
//			fiberp->v_y+=0.5*dt*(3.0*(fiberfy[i]+RHOS*gravity*M_PI*fiberp->radius*fiberp->radius)-fiberp->hy)/(RHOS*M_PI*fiberp->radius*fiberp->radius);
//			fiberp->hy=fiberfy[i]+RHOS*gravity*M_PI*fiberp->radius*fiberp->radius;
//		}
//		else{
//			fiberp->v_y+=0.5*dt*(3.0*(fiberfy[i]-RHOS*gravity*M_PI*fiberp->radius*fiberp->radius)-fiberp->hy)/(RHOS*M_PI*fiberp->radius*fiberp->radius);
//			fiberp->hy=fiberfy[i]-RHOS*gravity*M_PI*fiberp->radius*fiberp->radius;
//		}
//		printf("fy=%f ",fiberfy[i]);
//		fiberp->angvel+=0.5*dt*(3.0*fibern[i]-fiberp->ha)*2.0/M_PI/(RHOS*fiberp->radius*fiberp->radius*fiberp->radius*fiberp->radius);
//		printf("angvel=%f ang=%f ",fiberp->angvel,fiberp->angle);
//		fiberp->angvel=0.0;
//		fiberp->ha=fibern[i];
//	}
	
//	solidmove();
	
//	if(n<500){
//		for(int i=0; i<FIBER_NUM; i++){
//			_fiber_type* fiberp=&FIBER[i];
//			
//			fiberp->angvel=0.0;
//			fiberp->angle=0.0;
//			fiberp->x=64.0;
//			fiberp->y=48.0;
//			fiberp->v_x=0.0;
//			fiberp->v_y=0.0;
//		}
//	}
}
			//-------------- End of update ---------------//
			
			
		
		
		
			//-------------- Determination of rho --------------------//
void rhodefine(void) {
	for(int i=0; i<(XCELL_NUM+6); i++) {
		for(int j=0; j<(YCELL_NUM+6); j++) {
			index0 = j + i*(YCELL_NUM+6);
				
			cellp = &CELL[index0];
				if(cellp->phi<=PHIG) cellp->rho=RHOG;
				else if (cellp->phi>=PHIL) cellp->rho=RHOL;
				else cellp->rho=(RHOL+RHOG)*0.5+(RHOL-RHOG)*0.5*sin((cellp->phi-(PHIL+PHIG)*0.5)/(PHIL-PHIG)*M_PI);
			
				if(cellp->phi_x<=PHIG) cellp->rho_x=RHOG;
				else if (cellp->phi_x>=PHIL) cellp->rho_x=RHOL;
				else cellp->rho_x=(RHOL+RHOG)*0.5+(RHOL-RHOG)*0.5*sin((cellp->phi_x-(PHIL+PHIG)*0.5)/(PHIL-PHIG)*M_PI);
				
				if(cellp->phi_y<=PHIG) cellp->rho_y=RHOG;
				else if (cellp->phi_y>=PHIL) cellp->rho_y=RHOL;
				else cellp->rho_y=(RHOL+RHOG)*0.5+(RHOL-RHOG)*0.5*sin((cellp->phi_y-(PHIL+PHIG)*0.5)/(PHIL-PHIG)*M_PI);
		}
	}
	
/*	for(int i=0; i<(XCELL_NUM+6); i++) {
		for(int j=0; j<(YCELL_NUM+6); j++) {
			index0 = j + i*(YCELL_NUM+6);
				
			cellp = &CELL[index0];
			if(cellp->vof>0.99)cellp->rho=RHOS;
			if(cellp->vof_x>0.99)cellp->rho_x=RHOS;
			if(cellp->vof_y>0.99)cellp->rho_y=RHOS;
		}
	}
	*/
}
			//-------------- End of determination of rho --------------//
			
	
void ddrhodefine(void){
	double dx2i = 1.0/dx/dx;
	double dy2i = 1.0/dy/dy;
	double xpphi,xppphi,xmphi,xmmphi,ypphi,yppphi,ymphi,ymmphi;
	double xprho,xpprho,xmrho,xmmrho,yprho,ypprho,ymrho,ymmrho;
	
	for(int i=3; i<(XCELL_NUM+4); i++) {
		for(int j=3; j<(YCELL_NUM+4); j++) {
			indexsetdefine(i,j);
			
				ddrho[index0] = (-cellp_xpp->rho+16.0*(cellp_xp->rho+cellp_xm->rho)-30.0*cellp->rho-cellp_xmm->rho)/12.0*dx2i
								+(-cellp_ypp->rho+16.0*(cellp_yp->rho+cellp_ym->rho)-30.0*cellp->rho-cellp_ymm->rho)/12.0*dy2i;
		}
	}
	
	ddrhobnd();
	
	for(int i=3; i<(XCELL_NUM+3); i++) {
		for(int j=3; j<(YCELL_NUM+3); j++) {
			indexsetdefine(i,j);
			
			_flag_type* flagp=&FLAG[index0];
			
			if((flagp->c==2)||(flagp->c==4)||(flagp->c==6)){
				ddrho[index0]=interphi(index0,3);
			}
		}
	}
}		
		
void ddrhobnd(void){
#ifndef WALLY
	for(int i=0; i<(XCELL_NUM+6); i++){											//Periodic boundary
		indexbndy(i);
		
		ddrho[index_y0]   = ddrho[index_ycp0];
		ddrho[index_y1]   = ddrho[index_ycp1];
		ddrho[index_y2]   = ddrho[index_ycp2];
		ddrho[index_ycp3] = ddrho[index_y3];
		ddrho[index_ycp4] = ddrho[index_y4];
		ddrho[index_ycp5] = ddrho[index_y5];
	}
#endif
#ifndef WALLX
	for(int j=0; j<(YCELL_NUM+6); j++){
		indexbndx(j);
			
		ddrho[index_x0]   = ddrho[index_xcp0];
		ddrho[index_x1]   = ddrho[index_xcp1];
		ddrho[index_x2]   = ddrho[index_xcp2];
		ddrho[index_xcp3] = ddrho[index_x3];
		ddrho[index_xcp4] = ddrho[index_x4];
		ddrho[index_xcp5] = ddrho[index_x5];
	}
#endif
	
#ifdef WALLY
	for(int i=0; i<(XCELL_NUM+6); i++){											//Wall
		indexbndy(i);
			
			ddrho[index_y0]   = ddrho[index_y5];
			ddrho[index_y1]   = ddrho[index_y4];
			ddrho[index_y2]   = ddrho[index_y3];
			ddrho[index_ycp3] = ddrho[index_ycp2];
			ddrho[index_ycp4] = ddrho[index_ycp1];
			ddrho[index_ycp5] = ddrho[index_ycp0];
	}
#endif
#ifdef WALLX
	for(int j=0; j<(YCELL_NUM+6); j++){
		indexbndx(j);
			
		ddrho[index_x0]   = ddrho[index_x5];
		ddrho[index_x1]   = ddrho[index_x4];
		ddrho[index_x2]   = ddrho[index_x3];
		ddrho[index_xcp3] = ddrho[index_xcp2];
		ddrho[index_xcp4] = ddrho[index_xcp1];
		ddrho[index_xcp5] = ddrho[index_xcp0];
	}
#endif

}
		
			
			//----------- Determination of MU ---------------//
void mudefine(void) {
	double phixy;
	double rhoxy;
	
	for (int i=2; i<(XCELL_NUM+5); i++){
		for(int j=2; j<(YCELL_NUM+5); j++){
			int index_xm_ym = (j-1) + (i-1)*(YCELL_NUM+6);
			
			_cell_type* cellp_xm_ym = &CELL[index_xm_ym];
			
			indexsetdefine(i,j);
				
			_mu_type*   mup       = &TMPMU[index0];
			
			
			mup->c= MUG+(MUL-MUG)/(RHOL-RHOG)*(cellp->rho-RHOG);
			mup->x = MUG+(MUL-MUG)/(RHOL-RHOG)*(cellp->rho_x-RHOG);
			mup->y = MUG+(MUL-MUG)/(RHOL-RHOG)*(cellp->rho_y-RHOG);
			
			phixy = (cellp->phi + cellp_xm->phi + cellp_ym ->phi + cellp_xm_ym->phi)*0.25;
			
			if(phixy<=PHIG) rhoxy=RHOG;
			else if (phixy>=PHIL) rhoxy=RHOL;
			else rhoxy=(RHOL+RHOG)*0.5+(RHOL-RHOG)*0.5*sin((phixy-(PHIL+PHIG)*0.5)/(PHIL-PHIG)*M_PI);
			
			
			mup->xy = MUG+(MUL-MUG)/(RHOL-RHOG)*(rhoxy-RHOG);
		}
	}
}
			//----------- End of determination of MU ---------//
		
		
		
			//----------- Cahn-Hilliard eq. ---------------//
void CH_eq(void) {
	double zetax,zetay;
	double phidx3,phidy3;
	double dx2i = 1.0/dx/dx;
	double dy2i = 1.0/dy/dy;
	
	
	for(int i=3; i<(XCELL_NUM+4); i++) {
		for(int j=3; j<(YCELL_NUM+4); j++) {
			indexsetdefine(i,j);
			_flag_type* flagp=&FLAG[index0];
			
			if(flagp->c==1){
				phid2[index0] = (-cellp_xpp->phi+16.0*(cellp_xp->phi+cellp_xm->phi)-30.0*cellp->phi-cellp_xmm->phi)/12.0*dx2i
								+(-cellp_ypp->phi+16.0*(cellp_yp->phi+cellp_ym->phi)-30.0*cellp->phi-cellp_ymm->phi)/12.0*dy2i;
			}
		}
	}
	phid2bnd(phid2);
	
	for(int i=3; i<(XCELL_NUM+3); i++) {
		for(int j=3; j<(YCELL_NUM+3); j++) {
			indexsetdefine(i,j);
			_flag_type* flagp=&FLAG[index0];
			
			if((flagp->c==2)||(flagp->c==4)||(flagp->c==6)) {
				phid2[index0]=interphi(index0,2);
			}
			
		}
	}
	
	
	for(int i=3; i<(XCELL_NUM+3); i++) {
		for(int j=3; j<(YCELL_NUM+3); j++) {
			indexsetdefine(i,j);
			_flag_type* flagp=&FLAG[index0];
			
			if((flagp->c==2)||(flagp->c==4)||(flagp->c==6)) cellp->phi=interphi(index0,1);
			else if(flagp->c==0) cellp->phi=PHImax;
		}
	}
	
	phibnd();
	
	for(int i=2; i<(XCELL_NUM+5); i++) {
		for(int j=2; j<(YCELL_NUM+5); j++) {
			indexsetdefine(i,j);
			
			_flag_type* flagp=&FLAG[index0];
			
			
			_vold_type* voldp = &VOLD[index0];
			
			zetax = vdWT / (1.0-vdWb*cellp->phi_x) / (1.0-vdWb*cellp->phi_x) -2.0*vdWa*cellp->phi_x;
			zetay = vdWT / (1.0-vdWb*cellp->phi_y) / (1.0-vdWb*cellp->phi_y) -2.0*vdWa*cellp->phi_y;
			
//			if((flagp->x==1)||(flagp->x==2)){
				
				phidx3 = (-phid2[index_xp]+27.0*(phid2[index0]-phid2[index_xm])+phid2[index_xmm])/24.0/dx;		//Yi-1/2[i][j][k]
				Fx[index0] = cellp->phi_x*voldp->x+gamma*(-zetax*(-cellp_xp->phi+27.0*(cellp->phi-cellp_xm->phi)+cellp_xmm->phi)/dx/24.0
								+kappa2*cellp->phi_x*phidx3);		//phiSxS	//Yi-1/2[i][j][k]
//			}
			
//			if((flagp->y==1)||(flagp->y==2)){
				phidy3 = (-phid2[index_yp]+27.0*(phid2[index0]-phid2[index_ym])+phid2[index_ymm])/24.0/dy;		//Yj-1/2[i][j][k]
				Fy[index0] = cellp->phi_y*voldp->y+gamma*(-zetay*(-cellp_yp->phi+27.0*(cellp->phi-cellp_ym->phi)+cellp_ymm->phi)/dy/24.0
								+kappa2*cellp->phi_y*phidy3);
//			}
		}
	}
	
	
	Fbnd(Fx,Fy);
	
	
	for(int i=3; i<(XCELL_NUM+3); i++) {
		for(int j=3; j<(YCELL_NUM+3); j++) {
			indexsetdefine(i,j);
			_flag_type* flagp=&FLAG[index0];
			
			if((flagp->c==2)||(flagp->c==4)||(flagp->c==6)) cellp->phi=interphi(index0,0);
			else if(flagp->c==0) cellp->phi=PHImax;
		}
	}
	
	phibnd();
	
	for(int i=3; i<(XCELL_NUM+3); i++) {
		for(int j=3; j<(YCELL_NUM+3); j++) {
			indexsetdefine(i,j);
			_flag_type* flagp=&FLAG[index0];
			
			if(flagp->c==1){
				
				cellp->phi=cellp->phi-0.5*dt*(3.0*((Fx[index_xp]-Fx[index0])/dx
										  +(Fy[index_yp]-Fy[index0])/dy)-cellp->hf);
				cellp->hf= (Fx[index_xp]-Fx[index0])/dx
										  +(Fy[index_yp]-Fy[index0])/dy;
			}
		}
	}
	
	
	phibnd();
	
	
	for(int i=3; i<(XCELL_NUM+3); i++) {
		for(int j=3; j<(YCELL_NUM+3); j++) {
			indexsetdefine(i,j);
			_flag_type* flagp=&FLAG[index0];
			
			if((flagp->c==2)||(flagp->c==4)||(flagp->c==6)) cellp->phi=interphi(index0,0);
			else if((cellp->ls < 0.0)&&(cellp->ls > -3.51)) cellp->phi=interphi(index0,0);
			else if(flagp->c==0) cellp->phi=PHImax;
		}
	}
	
	phibnd();
	
	
	for(int i=2; i<(XCELL_NUM+5); i++) {
		for(int j=2; j<(YCELL_NUM+5); j++) {
			indexsetdefine(i,j);

			cellp->phi_x = (9.0*(cellp->phi+cellp_xm->phi)-cellp_xp->phi-cellp_xmm->phi)/16.0/dx;
			cellp->phi_y = (9.0*(cellp->phi+cellp_ym->phi)-cellp_yp->phi-cellp_ymm->phi)/16.0/dy;
		}
	}
	
	phikyokaibnd();
	phikyokaibnd();
	
	rhodefine();
}
	
			//----------- Cahn-Hilliard eq. --------//


		
void Fbnd(double Fx[(XCELL_NUM+6)*(YCELL_NUM+6)],double Fy[(XCELL_NUM+6)*(YCELL_NUM+6)]){
#ifdef WALLY
	for(int i=0; i<(XCELL_NUM+6); i++){										//Wall
		indexbndy(i);
			
		Fx[index_y0]   = -Fx[index_y5];
		Fx[index_y1]   = -Fx[index_y4];
		Fx[index_y2]   = -Fx[index_y3];
		Fx[index_ycp3] = -Fx[index_ycp2];
		Fx[index_ycp4] = -Fx[index_ycp1];
		Fx[index_ycp5] = -Fx[index_ycp0];
		Fy[index_y0]   = -Fy[index_y6];
		Fy[index_y1]   = -Fy[index_y5];
		Fy[index_y2]   = -Fy[index_y4];
		Fy[index_y3]   = 0.0;
		Fy[index_ycp3] = 0.0;
		Fy[index_ycp4] = -Fy[index_ycp2];
		Fy[index_ycp5] = -Fy[index_ycp1];
		
	}
#endif
#ifdef WALLX
	for (int j=0; j<(YCELL_NUM+6); j++){
		indexbndx(j);
			
		
		Fx[index_x0]   = -Fx[index_x6];
		Fx[index_x1]   = -Fx[index_x5];
		Fx[index_x2]   = -Fx[index_x4];
		Fx[index_x3]   = 0.0;
		Fx[index_xcp3] = 0.0;
		Fx[index_xcp4] = -Fx[index_xcp2];
		Fx[index_xcp5] = -Fx[index_xcp1];
		Fy[index_x0]   = -Fy[index_x5];
		Fy[index_x1]   = -Fy[index_x4];
		Fy[index_x2]   = -Fy[index_x3];
		Fy[index_xcp3] = -Fy[index_xcp2];
		Fy[index_xcp4] = -Fy[index_xcp1];
		Fy[index_xcp5] = -Fy[index_xcp0];
	}
#endif
}
		
			
			//---------- tec -------------//
void filewrite(void){
	double tecx,tecy,tecu,tecv,tecp,tecphi,tecrho;
	double divv;
	double lsflag,theta;
	char filename_tec[30];
	FILE *tecout;

	/* Tecplotpo */
	printf("Storing data for Tecplot\n");
	sprintf(filename_tec,"d-%d.tec",n);
	tecout = fopen(filename_tec,"wb");
	fprintf(tecout,"variables = \"x\",\"y\",\"u\",\"v\",\"pres\",\"divv\",\"phi\",\"rho\",\"flag\",\"flagx\",\"flagy\"\n");
//	fprintf(tecout,"variables = \"x\",\"y\",\"u\",\"v\",\"pres\",\"divv\",\"rho\"\n");
	fprintf(tecout,"zone i=%d j=%d f=point\n",XCELL_NUM+1,YCELL_NUM+1);
//	fe = freeeg();
	for(int j=3;j<=YCELL_NUM+3;j++){
		for(int i=3;i<=XCELL_NUM+3;i++){
			int index_xm_ym = (j-1) + (i-1)*(YCELL_NUM+6);
			
			
			_cell_type* cellp_xm_ym = &CELL[index_xm_ym];
			
		
			index0      = j     + i*(YCELL_NUM+6);
			index_xp    = j     + (i+1)*(YCELL_NUM+6);
			index_yp    = (j+1) + i*(YCELL_NUM+6);
			index_xm    = j     + (i-1)*(YCELL_NUM+6);
			index_ym    = (j-1) + i*(YCELL_NUM+6);
			
			cellp= &CELL[index0];
			cellp_xp = &CELL[index_xp];
			cellp_yp = &CELL[index_yp];
			cellp_xm = &CELL[index_xm];
			cellp_ym = &CELL[index_ym];
			_flag_type* flagp = &FLAG[index0];
			
			tecx = cellp->x-0.5*dx;
			tecy = cellp->y-0.5*dy;
			
			tecu = 0.5*(cellp->v_x + cellp_ym->v_x);
			tecv = 0.5*(cellp->v_y + cellp_xm->v_y);
			
//			tecp = 0.25*(cellp->press + cellp_xm->press + cellp_ym->press + cellp_xm_ym->press);
			tecp=cellp->press;
			divv = (cellp_xp->v_x-cellp->v_x)/dx+(cellp_yp->v_y-cellp->v_y)/dy;
			tecphi = 0.25*(cellp->phi + cellp_xm->phi + cellp_ym->phi + cellp_xm_ym->phi);
//			tecphi=cellp->phi;
//			tecphi = cellp->phi;
			tecrho = 0.25*(cellp->rho + cellp_xm->rho + cellp_ym->rho + cellp_xm_ym->rho);
//			tecrho=cellp->rho;
			lsflag=cellp->ls;
			
//			_fiber_type* fiberp=&FIBER[flagp->fiber];
//			if((cellp->ls>-5.01)&&(cellp->ls<-1.99)) {
//				theta=-fiberp->angle+atan2(cellp->y-fiberp->y,cellp->x-fiberp->x);
//				while((theta<-M_PI)||(theta>(M_PI+10e-6))){
//					if(theta<-M_PI) theta+=2.0*M_PI;
//					else if(theta>(M_PI+10e-6)) theta-=2.0*M_PI;
//				}
//					
//				lsflag=theta-M_PI*0.5;
//			}
			
			
			fprintf(tecout,"%6.4lf %6.4lf %20.16lf %20.16lf %20.16lf %20.16lf %20.16lf %20.16lf %20.16lf %20.16lf %20.16lf \n",tecx,tecy,tecu,tecv,tecp,divv,tecphi,tecrho,cellp->ls,cellp->vof_x,cellp->vof_y);
//			fprintf(tecout,"%6.4lf %6.4lf %20.16lf %20.16lf %20.16lf %20.16lf %20.16lf\n",tecx,tecy,tecu,tecv,tecp,divv,tecrho);
		}
	}
	fclose(tecout);
}


//-------------- Tentative flow velocity boundary conditions ---------------------//
void vhypobnd(void){
	
#ifndef WALLY
	
	for(int i=0; i<(XCELL_NUM+6); i++){										//Periodic boundary
		indexbndy(i);
			
		cellp_y0->v_xhypo   = cellp_ycp0->v_xhypo;
		cellp_y1->v_xhypo   = cellp_ycp1->v_xhypo;
		cellp_y2->v_xhypo   = cellp_ycp2->v_xhypo;
		cellp_ycp3->v_xhypo = cellp_y3->v_xhypo;
		cellp_ycp4->v_xhypo = cellp_y4->v_xhypo;
		cellp_ycp5->v_xhypo = cellp_y5->v_xhypo;
		cellp_y0->v_yhypo   = cellp_ycp0->v_yhypo;
		cellp_y1->v_yhypo   = cellp_ycp1->v_yhypo;
		cellp_y2->v_yhypo   = cellp_ycp2->v_yhypo;
		cellp_ycp3->v_yhypo = cellp_y3->v_yhypo;
		cellp_ycp4->v_yhypo = cellp_y4->v_yhypo;
		cellp_ycp5->v_yhypo = cellp_y5->v_yhypo;
	}
#endif
#ifndef WALLX
	for (int j=0; j<(YCELL_NUM+6); j++){
		indexbndx(j);
		
		
/*		cellp_x0->v_xhypo   = Uinf;		//IB1
		cellp_x1->v_xhypo   = Uinf;
		cellp_x2->v_xhypo   = Uinf;
		cellp_xcp3->v_xhypo = cellp_xcp2->v_xhypo;
		cellp_xcp4->v_xhypo = cellp_xcp2->v_xhypo;
		cellp_xcp5->v_xhypo = cellp_xcp2->v_xhypo;
		cellp_x0->v_yhypo   = 0.0;
		cellp_x1->v_yhypo   = 0.0;
		cellp_x2->v_yhypo   = 0.0;
		cellp_xcp3->v_yhypo = cellp_xcp2->v_yhypo;
		cellp_xcp4->v_yhypo = cellp_xcp2->v_yhypo;
		cellp_xcp5->v_yhypo = cellp_xcp2->v_yhypo;
*/		
		cellp_x0->v_xhypo   = cellp_xcp0->v_xhypo;
		cellp_x1->v_xhypo   = cellp_xcp1->v_xhypo;
		cellp_x2->v_xhypo   = cellp_xcp2->v_xhypo;
		cellp_xcp3->v_xhypo = cellp_x3->v_xhypo;
		cellp_xcp4->v_xhypo = cellp_x4->v_xhypo;
		cellp_xcp5->v_xhypo = cellp_x5->v_xhypo;
		cellp_x0->v_yhypo   = cellp_xcp0->v_yhypo;
		cellp_x1->v_yhypo   = cellp_xcp1->v_yhypo;
		cellp_x2->v_yhypo   = cellp_xcp2->v_yhypo;
		cellp_xcp3->v_yhypo = cellp_x3->v_yhypo;
		cellp_xcp4->v_yhypo = cellp_x4->v_yhypo;
		cellp_xcp5->v_yhypo = cellp_x5->v_yhypo;
	}
#endif

#ifdef WALLY
	for(int i=0; i<(XCELL_NUM+6); i++){										//Wall
		indexbndy(i);
			
/*		cellp_y0->v_xhypo   = -cellp_y5->v_xhypo;			//cavity
		cellp_y1->v_xhypo   = -cellp_y4->v_xhypo;
		cellp_y2->v_xhypo   = -cellp_y3->v_xhypo;
		cellp_ycp3->v_xhypo = -cellp_ycp2->v_xhypo+Uinf;
		cellp_ycp4->v_xhypo = -cellp_ycp1->v_xhypo+Uinf;
		cellp_ycp5->v_xhypo = -cellp_ycp0->v_xhypo+Uinf;
		cellp_y0->v_yhypo   = -cellp_y6->v_yhypo;
		cellp_y1->v_yhypo   = -cellp_y5->v_yhypo;
		cellp_y2->v_yhypo   = -cellp_y4->v_yhypo;
		cellp_y3->v_yhypo   = 0.0;
		cellp_ycp3->v_yhypo = 0.0;
		cellp_ycp4->v_yhypo = -cellp_ycp2->v_yhypo;
		cellp_ycp5->v_yhypo = -cellp_ycp1->v_yhypo;
*/		
		cellp_y0->v_xhypo   = -cellp_y5->v_xhypo;
		cellp_y1->v_xhypo   = -cellp_y4->v_xhypo;
		cellp_y2->v_xhypo   = -cellp_y3->v_xhypo;
		cellp_ycp3->v_xhypo = -cellp_ycp2->v_xhypo;
		cellp_ycp4->v_xhypo = -cellp_ycp1->v_xhypo;
		cellp_ycp5->v_xhypo = -cellp_ycp0->v_xhypo;
		cellp_y0->v_yhypo   = -cellp_y6->v_yhypo;
		cellp_y1->v_yhypo   = -cellp_y5->v_yhypo;
		cellp_y2->v_yhypo   = -cellp_y4->v_yhypo;
		cellp_y3->v_yhypo   = 0.0;
		cellp_ycp3->v_yhypo = 0.0;
		cellp_ycp4->v_yhypo = -cellp_ycp2->v_yhypo;
		cellp_ycp5->v_yhypo = -cellp_ycp1->v_yhypo;
	}
#endif
#ifdef WALLX
	for (int j=0; j<(YCELL_NUM+6); j++){
		indexbndx(j);
			
		cellp_x0->v_xhypo   = -cellp_x6->v_xhypo;
		cellp_x1->v_xhypo   = -cellp_x5->v_xhypo;
		cellp_x2->v_xhypo   = -cellp_x4->v_xhypo;
		cellp_x3->v_xhypo   = 0.0;
		cellp_xcp3->v_xhypo = 0.0;
		cellp_xcp4->v_xhypo = -cellp_xcp2->v_xhypo;
		cellp_xcp5->v_xhypo = -cellp_xcp1->v_xhypo;
		cellp_x0->v_yhypo   = -cellp_x5->v_yhypo;
		cellp_x1->v_yhypo   = -cellp_x4->v_yhypo;
		cellp_x2->v_yhypo   = -cellp_x3->v_yhypo;
		cellp_xcp3->v_yhypo = -cellp_xcp2->v_yhypo;
		cellp_xcp4->v_yhypo = -cellp_xcp1->v_yhypo;
		cellp_xcp5->v_yhypo = -cellp_xcp0->v_yhypo;
	}
#endif
}
			//-------------- End of tentative flow velocity boundary conditions ---------------//

			//-------------- Flow velocity boundary conditions -----------//
void kyokai(void){//"kyokai" is Japanese, it means boundary
#ifndef WALLY
	for(int i=0; i<(XCELL_NUM+6); i++){										//Periodic boundary
		indexbndy(i);
			
		cellp_y0->v_x   = cellp_ycp0->v_x;
		cellp_y1->v_x   = cellp_ycp1->v_x;
		cellp_y2->v_x   = cellp_ycp2->v_x;
		cellp_ycp3->v_x = cellp_y3->v_x;
		cellp_ycp4->v_x = cellp_y4->v_x;
		cellp_ycp5->v_x = cellp_y5->v_x;
		cellp_y0->v_y   = cellp_ycp0->v_y;
		cellp_y1->v_y   = cellp_ycp1->v_y;
		cellp_y2->v_y   = cellp_ycp2->v_y;
		cellp_ycp3->v_y = cellp_y3->v_y;
		cellp_ycp4->v_y = cellp_y4->v_y;
		cellp_ycp5->v_y = cellp_y5->v_y;
	}
#endif
#ifndef WALLX
	for (int j=0; j<(YCELL_NUM+6); j++){
		indexbndx(j);
/*			
		cellp_x0->v_x   = Uinf;		//IB1
		cellp_x1->v_x   = Uinf;
		cellp_x2->v_x   = Uinf;
		cellp_xcp3->v_x = cellp_xcp2->v_x;
		cellp_xcp4->v_x = cellp_xcp2->v_x;
		cellp_xcp5->v_x = cellp_xcp2->v_x;
		cellp_x0->v_y   = 0.0;
		cellp_x1->v_y   = 0.0;
		cellp_x2->v_y   = 0.0;
		cellp_xcp3->v_y = cellp_xcp2->v_y;
		cellp_xcp4->v_y = cellp_xcp2->v_y;
		cellp_xcp5->v_y = cellp_xcp2->v_y;
*/		
		cellp_x0->v_x   = cellp_xcp0->v_x;
		cellp_x1->v_x   = cellp_xcp1->v_x;
		cellp_x2->v_x   = cellp_xcp2->v_x;
		cellp_xcp3->v_x = cellp_x3->v_x;
		cellp_xcp4->v_x = cellp_x4->v_x;
		cellp_xcp5->v_x = cellp_x5->v_x;
		cellp_x0->v_y   = cellp_xcp0->v_y;
		cellp_x1->v_y   = cellp_xcp1->v_y;
		cellp_x2->v_y   = cellp_xcp2->v_y;
		cellp_xcp3->v_y = cellp_x3->v_y;
		cellp_xcp4->v_y = cellp_x4->v_y;
		cellp_xcp5->v_y = cellp_x5->v_y;
	}
#endif
#ifdef WALLY
	for(int i=0; i<(XCELL_NUM+6); i++){										//Wall
		indexbndy(i);
			
/*		cellp_y0->v_x   = -cellp_y5->v_x;	//cavity
		cellp_y1->v_x   = -cellp_y4->v_x;
		cellp_y2->v_x   = -cellp_y3->v_x;
		cellp_ycp3->v_x = -cellp_ycp2->v_x+2.0*Uinf;
		cellp_ycp4->v_x = -cellp_ycp1->v_x+2.0*Uinf;
		cellp_ycp5->v_x = -cellp_ycp0->v_x+2.0*Uinf;
		cellp_y0->v_y   = -cellp_y6->v_y;
		cellp_y1->v_y   = -cellp_y5->v_y;
		cellp_y2->v_y   = -cellp_y4->v_y;
		cellp_y3->v_y   = 0.0;
		cellp_ycp3->v_y = 0.0;
		cellp_ycp4->v_y = -cellp_ycp2->v_y;
		cellp_ycp5->v_y = -cellp_ycp1->v_y;
*/		
		cellp_y0->v_x   = -cellp_y5->v_x;
		cellp_y1->v_x   = -cellp_y4->v_x;
		cellp_y2->v_x   = -cellp_y3->v_x;
		cellp_ycp3->v_x = -cellp_ycp2->v_x;
		cellp_ycp4->v_x = -cellp_ycp1->v_x;
		cellp_ycp5->v_x = -cellp_ycp0->v_x;
		cellp_y0->v_y   = -cellp_y6->v_y;
		cellp_y1->v_y   = -cellp_y5->v_y;
		cellp_y2->v_y   = -cellp_y4->v_y;
		cellp_y3->v_y   = 0.0;
		cellp_ycp3->v_y = 0.0;
		cellp_ycp4->v_y = -cellp_ycp2->v_y;
		cellp_ycp5->v_y = -cellp_ycp1->v_y;
	}
#endif
#ifdef WALLX
	for (int j=1; j<(YCELL_NUM+5); j++){
		indexbndx(j);
			
		cellp_x0->v_x   = -cellp_x6->v_x;
		cellp_x1->v_x   = -cellp_x5->v_x;
		cellp_x2->v_x   = -cellp_x4->v_x;
		cellp_x3->v_x   = 0.0;
		cellp_xcp3->v_x = 0.0;
		cellp_xcp4->v_x = -cellp_xcp2->v_x;
		cellp_xcp5->v_x = -cellp_xcp1->v_x;
		cellp_x0->v_y   = -cellp_x5->v_y;
		cellp_x1->v_y   = -cellp_x4->v_y;
		cellp_x2->v_y   = -cellp_x3->v_y;
		cellp_xcp3->v_y = -cellp_xcp2->v_y;
		cellp_xcp4->v_y = -cellp_xcp1->v_y;
		cellp_xcp5->v_y = -cellp_xcp0->v_y;
	}
#endif
}
			//-------------- End of flow velocity boundary conditions -----//

//-------------- E -----------------//
void pressbnd(void){
/*	for(int i=3; i<(XCELL_NUM+3); i++){
		for(int j=3; j<(YCELL_NUM+3); j++){
			indexSORdefine(i,j);
			_flag_type* flagp=&FLAG[index0];
			
			if((flagp->c==2)||(flagp->c==4)) cellp->press=interpress(index0,4);
//			else if(flagp->c==4) cellp->press=interphi(index0,4);
		}
	}
*/	
#ifndef WALLY
	for(int i=0; i<(XCELL_NUM+6); i++){											//E
		indexbndy(i);
			
		cellp_y0->press   = cellp_ycp0->press;
		cellp_y1->press   = cellp_ycp1->press;
		cellp_y2->press   = cellp_ycp2->press;
		cellp_ycp3->press = cellp_y3->press;
		cellp_ycp4->press = cellp_y4->press;
		cellp_ycp5->press = cellp_y5->press;
	}
#endif
#ifndef WALLX
	for(int j=0; j<(YCELL_NUM+6); j++){
		indexbndx(j);
		
/*		cellp_x0->press   = cellp_x3->press;		//IB1
		cellp_x1->press   = cellp_x3->press;
		cellp_x2->press   = cellp_x3->press;
		cellp_xcp3->press = cellp_xcp2->press;
		cellp_xcp4->press = cellp_xcp2->press;
		cellp_xcp5->press = cellp_xcp2->press;
*/		
		cellp_x0->press   = cellp_xcp0->press;
		cellp_x1->press   = cellp_xcp1->press;
		cellp_x2->press   = cellp_xcp2->press;
		cellp_xcp3->press = cellp_x3->press;
		cellp_xcp4->press = cellp_x4->press;
		cellp_xcp5->press = cellp_x5->press;
	}
#endif

#ifdef WALLY
	for(int i=0; i<(XCELL_NUM+6); i++){
		indexbndy(i);
		
		cellp_y0->press   = cellp_y5->press;
		cellp_y1->press   = cellp_y4->press;
		cellp_y2->press   = cellp_y3->press;
		cellp_ycp3->press = cellp_ycp2->press;
		cellp_ycp4->press = cellp_ycp1->press;
		cellp_ycp5->press = cellp_ycp0->press;
	}
#endif
#ifdef WALLX
	for(int j=0; j<(YCELL_NUM+6); j++){
		indexbndx(j);
		
		cellp_x0->press   = cellp_x5->press;
		cellp_x1->press   = cellp_x4->press;
		cellp_x2->press   = cellp_x3->press;
		cellp_xcp3->press = cellp_xcp2->press;
		cellp_xcp4->press = cellp_xcp1->press;
		cellp_xcp5->press = cellp_xcp0->press;
	}
#endif
}
		
			//---------------- EI ------------//

		
void vofbnd(void){
#ifndef WALLY
	for(int i=0; i<(XCELL_NUM+6); i++){											//E
		indexbndy(i);
			
		cellp_y0->vof   = cellp_ycp0->vof;
		cellp_y1->vof   = cellp_ycp1->vof;
		cellp_y2->vof   = cellp_ycp2->vof;
		cellp_ycp3->vof = cellp_y3->vof;
		cellp_ycp4->vof = cellp_y4->vof;
		cellp_ycp5->vof = cellp_y5->vof;
		cellp_y0->vof_y   = cellp_ycp0->vof_y;
		cellp_y1->vof_y   = cellp_ycp1->vof_y;
		cellp_y2->vof_y   = cellp_ycp2->vof_y;
		cellp_ycp3->vof_y = cellp_y3->vof_y;
		cellp_ycp4->vof_y = cellp_y4->vof_y;
		cellp_ycp5->vof_y = cellp_y5->vof_y;
		cellp_y0->vof_x   = cellp_ycp0->vof_x;
		cellp_y1->vof_x   = cellp_ycp1->vof_x;
		cellp_y2->vof_x   = cellp_ycp2->vof_x;
		cellp_ycp3->vof_x = cellp_y3->vof_x;
		cellp_ycp4->vof_x = cellp_y4->vof_x;
		cellp_ycp5->vof_x = cellp_y5->vof_x;
		
	}
#endif
#ifndef WALLX
	for(int j=0; j<(YCELL_NUM+6); j++){
		indexbndx(j);
		
	
		cellp_x0->vof   = cellp_xcp0->vof;
		cellp_x1->vof   = cellp_xcp1->vof;
		cellp_x2->vof   = cellp_xcp2->vof;
		cellp_xcp3->vof = cellp_x3->vof;
		cellp_xcp4->vof = cellp_x4->vof;
		cellp_xcp5->vof = cellp_x5->vof;
		
		cellp_x0->vof_x   = cellp_xcp0->vof_x;
		cellp_x1->vof_x   = cellp_xcp1->vof_x;
		cellp_x2->vof_x   = cellp_xcp2->vof_x;
		cellp_xcp3->vof_x = cellp_x3->vof_x;
		cellp_xcp4->vof_x = cellp_x4->vof_x;
		cellp_xcp5->vof_x = cellp_x5->vof_x;
		
		cellp_x0->vof_y   = cellp_xcp0->vof_y;
		cellp_x1->vof_y   = cellp_xcp1->vof_y;
		cellp_x2->vof_y   = cellp_xcp2->vof_y;
		cellp_xcp3->vof_y = cellp_x3->vof_y;
		cellp_xcp4->vof_y = cellp_x4->vof_y;
		cellp_xcp5->vof_y = cellp_x5->vof_y;
	}
#endif

#ifdef WALLY
	for(int i=0; i<(XCELL_NUM+6); i++){
		indexbndy(i);
		
		cellp_y0->vof   = cellp_y5->vof;
		cellp_y1->vof   = cellp_y4->vof;
		cellp_y2->vof   = cellp_y3->vof;
		cellp_ycp3->vof = cellp_ycp2->vof;
		cellp_ycp4->vof = cellp_ycp1->vof;
		cellp_ycp5->vof = cellp_ycp0->vof;
		
		cellp_y0->vof_x   = cellp_y5->vof_x;
		cellp_y1->vof_x   = cellp_y4->vof_x;
		cellp_y2->vof_x   = cellp_y3->vof_x;
		cellp_ycp3->vof_x = cellp_ycp2->vof_x;
		cellp_ycp4->vof_x = cellp_ycp1->vof_x;
		cellp_ycp5->vof_x = cellp_ycp0->vof_x;
		
		cellp_y0->vof_y   = cellp_y5->vof_y;
		cellp_y1->vof_y   = cellp_y4->vof_y;
		cellp_y2->vof_y   = cellp_y3->vof_y;
		cellp_ycp3->vof_y = cellp_ycp2->vof_y;
		cellp_ycp4->vof_y = cellp_ycp1->vof_y;
		cellp_ycp5->vof_y = cellp_ycp0->vof_y;
	}
#endif
#ifdef WALLX
	for(int j=0; j<(YCELL_NUM+6); j++){
		indexbndx(j);
		
		cellp_x0->vof   = cellp_x5->vof;
		cellp_x1->vof   = cellp_x4->vof;
		cellp_x2->vof   = cellp_x3->vof;
		cellp_xcp3->vof = cellp_xcp2->vof;
		cellp_xcp4->vof = cellp_xcp1->vof;
		cellp_xcp5->vof = cellp_xcp0->vof;
		
		cellp_x0->vof_x   = cellp_x5->vof_x;
		cellp_x1->vof_x   = cellp_x4->vof_x;
		cellp_x2->vof_x   = cellp_x3->vof_x;
		cellp_xcp3->vof_x = cellp_xcp2->vof_x;
		cellp_xcp4->vof_x = cellp_xcp1->vof_x;
		cellp_xcp5->vof_x = cellp_xcp0->vof_x;
		
		cellp_x0->vof_y   = cellp_x5->vof_y;
		cellp_x1->vof_y   = cellp_x4->vof_y;
		cellp_x2->vof_y   = cellp_x3->vof_y;
		cellp_xcp3->vof_y = cellp_xcp2->vof_y;
		cellp_xcp4->vof_y = cellp_xcp1->vof_y;
		cellp_xcp5->vof_y = cellp_xcp0->vof_y;
	}
#endif
}
		
		
		
void flagbnd(void){
#ifndef WALLY
	for(int i=0; i<(XCELL_NUM+6); i++){											//E
		indexbndy(i);
		
		_flag_type* flagp_y0   = &FLAG[index_y0];
		_flag_type* flagp_y1   = &FLAG[index_y1];
		_flag_type* flagp_y2   = &FLAG[index_y2];
		_flag_type* flagp_y3   = &FLAG[index_y3];
		_flag_type* flagp_y4   = &FLAG[index_y4];
		_flag_type* flagp_y5   = &FLAG[index_y5];
		_flag_type* flagp_y6   = &FLAG[index_y6];
		_flag_type* flagp_ycp0 = &FLAG[index_ycp0];
		_flag_type* flagp_ycp1 = &FLAG[index_ycp1];
		_flag_type* flagp_ycp2 = &FLAG[index_ycp2];
		_flag_type* flagp_ycp3 = &FLAG[index_ycp3];
		_flag_type* flagp_ycp4 = &FLAG[index_ycp4];
		_flag_type* flagp_ycp5 = &FLAG[index_ycp5];
		
		
		flagp_y0->c   = flagp_ycp0->c;
		flagp_y1->c   = flagp_ycp1->c;
		flagp_y2->c   = flagp_ycp2->c;
		flagp_ycp3->c = flagp_y3->c;
		flagp_ycp4->c = flagp_y4->c;
		flagp_ycp5->c = flagp_y5->c;
		flagp_y0->x   = flagp_ycp0->x;
		flagp_y1->x   = flagp_ycp1->x;
		flagp_y2->x   = flagp_ycp2->x;
		flagp_ycp3->x = flagp_y3->x;
		flagp_ycp4->x = flagp_y4->x;
		flagp_ycp5->x = flagp_y5->x;
		flagp_y0->y   = flagp_ycp0->y;
		flagp_y1->y   = flagp_ycp1->y;
		flagp_y2->y   = flagp_ycp2->y;
		flagp_ycp3->y = flagp_y3->y;
		flagp_ycp4->y = flagp_y4->y;
		flagp_ycp5->y = flagp_y5->y;
	}
#endif
#ifndef WALLX
	for(int j=0; j<(YCELL_NUM+6); j++){
		indexbndx(j);
		
		_flag_type* flagp_x0   = &FLAG[index_x0];
		_flag_type* flagp_x1   = &FLAG[index_x1];
		_flag_type* flagp_x2   = &FLAG[index_x2];
		_flag_type* flagp_x3   = &FLAG[index_x3];
		_flag_type* flagp_x4   = &FLAG[index_x4];
		_flag_type* flagp_x5   = &FLAG[index_x5];
		_flag_type* flagp_x6   = &FLAG[index_x6];
		_flag_type* flagp_xcp0 = &FLAG[index_xcp0];
		_flag_type* flagp_xcp1 = &FLAG[index_xcp1];
		_flag_type* flagp_xcp2 = &FLAG[index_xcp2];
		_flag_type* flagp_xcp3 = &FLAG[index_xcp3];
		_flag_type* flagp_xcp4 = &FLAG[index_xcp4];
		_flag_type* flagp_xcp5 = &FLAG[index_xcp5];
		
		flagp_x0->c   = flagp_xcp0->c;
		flagp_x1->c   = flagp_xcp1->c;
		flagp_x2->c   = flagp_xcp2->c;
		flagp_xcp3->c = flagp_x3->c;
		flagp_xcp4->c = flagp_x4->c;
		flagp_xcp5->c = flagp_x5->c;
		flagp_x0->x   = flagp_xcp0->x;
		flagp_x1->x   = flagp_xcp1->x;
		flagp_x2->x   = flagp_xcp2->x;
		flagp_xcp3->x = flagp_x3->x;
		flagp_xcp4->x = flagp_x4->x;
		flagp_xcp5->x = flagp_x5->x;
		flagp_x0->y   = flagp_xcp0->y;
		flagp_x1->y   = flagp_xcp1->y;
		flagp_x2->y   = flagp_xcp2->y;
		flagp_xcp3->y = flagp_x3->y;
		flagp_xcp4->y = flagp_x4->y;
		flagp_xcp5->y = flagp_x5->y;
	}
#endif

#ifdef WALLY
	for(int i=0; i<(XCELL_NUM+6); i++){
		indexbndy(i);
		
		_flag_type* flagp_y0   = &FLAG[index_y0];
		_flag_type* flagp_y1   = &FLAG[index_y1];
		_flag_type* flagp_y2   = &FLAG[index_y2];
		_flag_type* flagp_y3   = &FLAG[index_y3];
		_flag_type* flagp_y4   = &FLAG[index_y4];
		_flag_type* flagp_y5   = &FLAG[index_y5];
		_flag_type* flagp_y6   = &FLAG[index_y6];
		_flag_type* flagp_ycp0 = &FLAG[index_ycp0];
		_flag_type* flagp_ycp1 = &FLAG[index_ycp1];
		_flag_type* flagp_ycp2 = &FLAG[index_ycp2];
		_flag_type* flagp_ycp3 = &FLAG[index_ycp3];
		_flag_type* flagp_ycp4 = &FLAG[index_ycp4];
		_flag_type* flagp_ycp5 = &FLAG[index_ycp5];
		
		flagp_y0->c   = flagp_y5->c;
		flagp_y1->c   = flagp_y4->c;
		flagp_y2->c   = flagp_y3->c;
		flagp_ycp3->c = flagp_ycp2->c;
		flagp_ycp4->c = flagp_ycp1->c;
		flagp_ycp5->c = flagp_ycp0->c;
		flagp_y0->x   = flagp_y5->x;
		flagp_y1->x   = flagp_y4->x;
		flagp_y2->x   = flagp_y3->x;
		flagp_ycp3->x = flagp_ycp2->x;
		flagp_ycp4->x = flagp_ycp1->x;
		flagp_ycp5->x = flagp_ycp0->x;
		flagp_y0->y   = flagp_y5->y;
		flagp_y1->y   = flagp_y4->y;
		flagp_y2->y   = flagp_y3->y;
		flagp_ycp3->y = flagp_ycp2->y;
		flagp_ycp4->y = flagp_ycp1->y;
		flagp_ycp5->y = flagp_ycp0->y;
	}
#endif
#ifdef WALLX
	for(int j=0; j<(YCELL_NUM+6); j++){
		indexbndx(j);
		
		_flag_type* flagp_x0   = &FLAG[index_x0];
		_flag_type* flagp_x1   = &FLAG[index_x1];
		_flag_type* flagp_x2   = &FLAG[index_x2];
		_flag_type* flagp_x3   = &FLAG[index_x3];
		_flag_type* flagp_x4   = &FLAG[index_x4];
		_flag_type* flagp_x5   = &FLAG[index_x5];
		_flag_type* flagp_x6   = &FLAG[index_x6];
		_flag_type* flagp_xcp0 = &FLAG[index_xcp0];
		_flag_type* flagp_xcp1 = &FLAG[index_xcp1];
		_flag_type* flagp_xcp2 = &FLAG[index_xcp2];
		_flag_type* flagp_xcp3 = &FLAG[index_xcp3];
		_flag_type* flagp_xcp4 = &FLAG[index_xcp4];
		_flag_type* flagp_xcp5 = &FLAG[index_xcp5];
		
		flagp_x0->c   = flagp_x5->c;
		flagp_x1->c   = flagp_x4->c;
		flagp_x2->c   = flagp_x3->c;
		flagp_xcp3->c = flagp_xcp2->c;
		flagp_xcp4->c = flagp_xcp1->c;
		flagp_xcp5->c = flagp_xcp0->c;
		flagp_x0->x   = flagp_x5->x;
		flagp_x1->x   = flagp_x4->x;
		flagp_x2->x   = flagp_x3->x;
		flagp_xcp3->x = flagp_xcp2->x;
		flagp_xcp4->x = flagp_xcp1->x;
		flagp_xcp5->x = flagp_xcp0->x;
		flagp_x0->y   = flagp_x5->y;
		flagp_x1->y   = flagp_x4->y;
		flagp_x2->y   = flagp_x3->y;
		flagp_xcp3->y = flagp_xcp2->y;
		flagp_xcp4->y = flagp_xcp1->y;
		flagp_xcp5->y = flagp_xcp0->y;
	}
#endif
}



void phibnd(void){
#ifndef WALLY
	for(int i=0; i<(XCELL_NUM+6); i++){
		indexbndy(i);
		
		cellp_y0->phi   = cellp_ycp0->phi;
		cellp_y1->phi   = cellp_ycp1->phi;
		cellp_y2->phi   = cellp_ycp2->phi;
		cellp_ycp3->phi = cellp_y3->phi;
		cellp_ycp4->phi = cellp_y4->phi;
		cellp_ycp5->phi = cellp_y5->phi;
	}
#endif
#ifndef WALLX
	for(int j=0; j<(YCELL_NUM+6); j++){
		indexbndx(j);
			
		cellp_x0->phi   = cellp_xcp0->phi;
		cellp_x1->phi   = cellp_xcp1->phi;
		cellp_x2->phi   = cellp_xcp2->phi;
		cellp_xcp3->phi = cellp_x3->phi;
		cellp_xcp4->phi = cellp_x4->phi;
		cellp_xcp5->phi = cellp_x5->phi;
	}
#endif
	
#ifdef WALLY
	for(int i=0; i<(XCELL_NUM+6); i++){
		indexbndy(i);
			
		if((cellp_y3->phi>PHIG)&&(cellp_y3->phi<PHIL)){
			cellp_y0->phi   = cellp_y5->phi-wetphiwall*dy*5.0;
			cellp_y1->phi   = cellp_y4->phi-wetphiwall*dy*3.0;
			cellp_y2->phi   = cellp_y3->phi-wetphiwall*dy*1.0;
		}
		else {
			cellp_y0->phi   = cellp_y5->phi;
			cellp_y1->phi   = cellp_y4->phi;
			cellp_y2->phi   = cellp_y3->phi;
		}
		if((cellp_ycp2->phi>PHIG)&&(cellp_ycp2->phi<PHIL)){
			cellp_ycp3->phi = cellp_ycp2->phi-wetphiwall*dy*1.0;
			cellp_ycp4->phi = cellp_ycp1->phi-wetphiwall*dy*3.0;
				cellp_ycp5->phi = cellp_ycp0->phi-wetphiwall*dy*5.0;
		}
		else {
			cellp_ycp3->phi = cellp_ycp2->phi;
			cellp_ycp4->phi = cellp_ycp1->phi;
			cellp_ycp5->phi = cellp_ycp0->phi;
		}
	}
#endif
#ifdef WALLX
	for(int j=0; j<(YCELL_NUM+6); j++){
		indexbndx(j);
			
		if((cellp_x3->phi>PHIG)&&(cellp_x3->phi<PHIL)){
			cellp_x0->phi   = cellp_x5->phi-wetphiwall*dx*5.0;
			cellp_x1->phi   = cellp_x4->phi-wetphiwall*dx*3.0;
			cellp_x2->phi   = cellp_x3->phi-wetphiwall*dx*1.0;
		}
		else{
			cellp_x0->phi   = cellp_x5->phi;
			cellp_x1->phi   = cellp_x4->phi;
			cellp_x2->phi   = cellp_x3->phi;
			}
		if((cellp_xcp2->phi>PHIG)&&(cellp_xcp2->phi<PHIL)){
			cellp_xcp3->phi = cellp_xcp2->phi-wetphiwall*dx*1.0;
			cellp_xcp4->phi = cellp_xcp1->phi-wetphiwall*dx*3.0;
			cellp_xcp5->phi = cellp_xcp0->phi-wetphiwall*dx*5.0;
		}
		else{
			cellp_xcp3->phi = cellp_xcp2->phi;
			cellp_xcp4->phi = cellp_xcp1->phi;
			cellp_xcp5->phi = cellp_xcp0->phi;
		}
	}
#endif
}


void phikyokaibnd(void){
#ifndef WALLY
	for(int i=0; i<(XCELL_NUM+6); i++){
		indexbndy(i);
		
		cellp_y0->phi_x   = cellp_ycp0->phi_x;
		cellp_y1->phi_x   = cellp_ycp1->phi_x;
		cellp_y2->phi_x   = cellp_ycp2->phi_x;
		cellp_ycp3->phi_x = cellp_y3->phi_x;
		cellp_ycp4->phi_x = cellp_y4->phi_x;
		cellp_ycp5->phi_x = cellp_y5->phi_x;
		cellp_y0->phi_y   = cellp_ycp0->phi_y;
		cellp_y1->phi_y   = cellp_ycp1->phi_y;
		cellp_y2->phi_y   = cellp_ycp2->phi_y;
		cellp_ycp3->phi_y = cellp_y3->phi_y;
		cellp_ycp4->phi_y = cellp_y4->phi_y;
		cellp_ycp5->phi_y = cellp_y5->phi_y;
	}
#endif
#ifndef WALLX
	for (int j=0; j<(YCELL_NUM+6); j++){
		indexbndx(j);
			
		cellp_x0->phi_x   = cellp_xcp0->phi_x;
		cellp_x1->phi_x   = cellp_xcp1->phi_x;
		cellp_x2->phi_x   = cellp_xcp2->phi_x;
		cellp_xcp3->phi_x = cellp_x3->phi_x;
		cellp_xcp4->phi_x = cellp_x4->phi_x;
		cellp_xcp5->phi_x = cellp_x5->phi_x;
		cellp_x0->phi_y   = cellp_xcp0->phi_y;
		cellp_x1->phi_y   = cellp_xcp1->phi_y;
		cellp_x2->phi_y   = cellp_xcp2->phi_y;
		cellp_xcp3->phi_y = cellp_x3->phi_y;
		cellp_xcp4->phi_y = cellp_x4->phi_y;
		cellp_xcp5->phi_y = cellp_x5->phi_y;
	}
#endif
#ifdef WALLY
	for(int i=0; i<(XCELL_NUM+6); i++){
		indexbndy(i);
			
		if((cellp_y3->phi_x>PHIG)&&(cellp_y3->phi_x<PHIL)){
			cellp_y0->phi_x   = cellp_y5->phi_x-wetphiwall*dx*5.0;
			cellp_y1->phi_x   = cellp_y4->phi_x-wetphiwall*dx*3.0;
			cellp_y2->phi_x   = cellp_y3->phi_x-wetphiwall*dx*1.0;
		}
		else {
			cellp_y0->phi_x   = cellp_y5->phi_x;
			cellp_y1->phi_x   = cellp_y4->phi_x;
			cellp_y2->phi_x   = cellp_y3->phi_x;
		}
		if((cellp_ycp2->phi_x>PHIG)&&(cellp_ycp2->phi_x<PHIL)){
			cellp_ycp3->phi_x = cellp_ycp2->phi_x-wetphiwall*dx*1.0;
			cellp_ycp4->phi_x = cellp_ycp1->phi_x-wetphiwall*dx*3.0;
			cellp_ycp5->phi_x = cellp_ycp0->phi_x-wetphiwall*dx*5.0;
		}
		else {
			cellp_ycp3->phi_x = cellp_ycp2->phi_x;
			cellp_ycp4->phi_x = cellp_ycp1->phi_x;
			cellp_ycp5->phi_x = cellp_ycp0->phi_x;
		}
		if((cellp_y3->phi_y>PHIG)&&(cellp_y3->phi_y<PHIL)){
			cellp_y0->phi_y   = cellp_y6->phi_y-wetphiwall*dy*6.0;
			cellp_y1->phi_y   = cellp_y5->phi_y-wetphiwall*dy*4.0;
			cellp_y2->phi_y   = cellp_y4->phi_y-wetphiwall*dy*2.0;
		}
		else {
			cellp_y0->phi_y   = cellp_y6->phi_y;
			cellp_y1->phi_y   = cellp_y5->phi_y;
			cellp_y2->phi_y   = cellp_y4->phi_y;
		}
		if((cellp_ycp3->phi_y>PHIG)&&(cellp_ycp3->phi_y<PHIL)){
			cellp_ycp4->phi_y = cellp_ycp2->phi_y-wetphiwall*dy*2.0;
			cellp_ycp5->phi_y = cellp_ycp1->phi_y-wetphiwall*dy*4.0;
		}
		else {
			cellp_ycp4->phi_y = cellp_ycp2->phi_y;
			cellp_ycp5->phi_y = cellp_ycp1->phi_y;
		}
	}
#endif
#ifdef WALLX
	for (int j=0; j<(YCELL_NUM+6); j++){
		indexbndx(j);
			
		if((cellp_x3->phi_x>PHIG)&&(cellp_x3->phi_x<PHIL)){
			cellp_x0->phi_x   = cellp_x6->phi_x-wetphiwall*dx*6.0;
			cellp_x1->phi_x   = cellp_x5->phi_x-wetphiwall*dx*4.0;
			cellp_x2->phi_x   = cellp_x4->phi_x-wetphiwall*dx*2.0;
		}
		else {
			cellp_x0->phi_x   = cellp_x6->phi_x;
			cellp_x1->phi_x   = cellp_x5->phi_x;
			cellp_x2->phi_x   = cellp_x4->phi_x;
		}
		if((cellp_xcp3->phi_x>PHIG)&&(cellp_xcp3->phi_x<PHIL)){
			cellp_xcp4->phi_x = cellp_xcp2->phi_x-wetphiwall*dx*2.0;
			cellp_xcp5->phi_x = cellp_xcp1->phi_x-wetphiwall*dx*4.0;
		}
		else {
			cellp_xcp4->phi_x = cellp_xcp2->phi_x;
			cellp_xcp5->phi_x = cellp_xcp1->phi_x;
		}
		if((cellp_x3->phi_y>PHIG)&&(cellp_x3->phi_y<PHIL)){
			cellp_x0->phi_y   = cellp_x5->phi_y-wetphiwall*dy*5.0;
			cellp_x1->phi_y   = cellp_x4->phi_y-wetphiwall*dy*3.0;
				cellp_x2->phi_y   = cellp_x3->phi_y-wetphiwall*dy*1.0;
		}
		else {
			cellp_x0->phi_y   = cellp_x5->phi_y;
			cellp_x1->phi_y   = cellp_x4->phi_y;
			cellp_x2->phi_y   = cellp_x3->phi_y;
		}
		if((cellp_xcp2->phi_y>PHIG)&&(cellp_xcp2->phi_y<PHIL)){
			cellp_xcp3->phi_y = cellp_xcp2->phi_y-wetphiwall*dy*1.0;
			cellp_xcp4->phi_y = cellp_xcp1->phi_y-wetphiwall*dy*3.0;
			cellp_xcp5->phi_y = cellp_xcp0->phi_y-wetphiwall*dy*5.0;
			}
		else {
			cellp_xcp3->phi_y = cellp_xcp2->phi_y;
			cellp_xcp4->phi_y = cellp_xcp1->phi_y;
			cellp_xcp5->phi_y = cellp_xcp0->phi_y;
		}
	}
#endif
}


void phid2bnd(double phid2[(XCELL_NUM+6)*(YCELL_NUM+6)]){
#ifndef WALLY
	for(int i=0; i<(XCELL_NUM+6); i++){
		indexbndy(i);
			
		phid2[index_y0]   = phid2[index_ycp0];
		phid2[index_y1]   = phid2[index_ycp1];
		phid2[index_y2]   = phid2[index_ycp2];
		phid2[index_ycp3] = phid2[index_y3];
		phid2[index_ycp4] = phid2[index_y4];
		phid2[index_ycp5] = phid2[index_y5];
	}
#endif
#ifndef WALLX	
	for(int j=0; j<(YCELL_NUM+6); j++){
		indexbndx(j);
		
		phid2[index_x0]   = phid2[index_xcp0];
		phid2[index_x1]   = phid2[index_xcp1];
		phid2[index_x2]   = phid2[index_xcp2];
		phid2[index_xcp3] = phid2[index_x3];
		phid2[index_xcp4] = phid2[index_x4];
		phid2[index_xcp5] = phid2[index_x5];
	}
#endif

#ifdef WALLY
	for(int i=0; i<(XCELL_NUM+6); i++){
		indexbndy(i);
		
		phid2[index_y0]   = phid2[index_y5];
		phid2[index_y1]   = phid2[index_y4];
		phid2[index_y2]   = phid2[index_y3];
		phid2[index_ycp3] = phid2[index_ycp2];
		phid2[index_ycp4] = phid2[index_ycp1];
		phid2[index_ycp5] = phid2[index_ycp0];
	}
#endif
#ifdef WALLX
	for(int j=0; j<(YCELL_NUM+6); j++){
		indexbndx(j);
			
		phid2[index_x0]   = phid2[index_x5];
		phid2[index_x1]   = phid2[index_x4];
		phid2[index_x2]   = phid2[index_x3];
		phid2[index_xcp3] = phid2[index_xcp2];
		phid2[index_xcp4] = phid2[index_xcp1];
		phid2[index_xcp5] = phid2[index_xcp0];
	}
#endif
}

		

void restartwrite(void){
	FILE* fp;
	char filename_restart[30];
	
	printf("Storing data for Restart\n");
	sprintf(filename_restart,"restart-%d.dat",n);
	
	if((fp = fopen(filename_restart, "ab"))==NULL){
		printf("file could not be opened\n");
		exit(1);
	}
	
	
	if(fwrite(&n, sizeof(int), 1, fp) != 1){
		printf("file could not be opened");
		exit(1);
	}
	
	for(int i=0; i<(XCELL_NUM+6); i++) {
		for(int j=0; j<(YCELL_NUM+6); j++) {
			index0 = j + i*(YCELL_NUM+6);
			cellp= &CELL[index0];
			
			if(fwrite(&(cellp->v_x), sizeof(double), 1, fp) != 1){
				printf("file could not be opened");
				exit(1);
			}
			if(fwrite(&(cellp->v_y), sizeof(double), 1, fp) != 1){
				printf("file could not be opened");
				exit(1);
			}
			if(fwrite(&(cellp->press), sizeof(double), 1, fp) != 1){
				printf("file could not be opened");
				exit(1);
			}
			if(fwrite(&(cellp->phi), sizeof(double), 1, fp) != 1){
				printf("file could not be opened");
				exit(1);
			}
			if(fwrite(&(cellp->hx), sizeof(double), 1, fp) != 1){
				printf("file could not be opened");
				exit(1);
			}
			if(fwrite(&(cellp->hy), sizeof(double), 1, fp) != 1){
				printf("file could not be opened");
				exit(1);
			}
		}
	}
//	for(int i=0; i<FIBER_NUM; i++){
//		_fiber_type* fiberp = &FIBER[i];
//		
//		if(fwrite(&(fiberp->x), sizeof(double), 1, fp) != 1){
//			printf("file could not be opened");
//			exit(1);
//		}
//		if(fwrite(&(fiberp->y), sizeof(double), 1, fp) != 1){
//			printf("file could not be opened");
//			exit(1);
//		}
//		if(fwrite(&(fiberp->v_x), sizeof(double), 1, fp) != 1){
//			printf("file could not be opened");
//			exit(1);
//		}
//		if(fwrite(&(fiberp->v_y), sizeof(double), 1, fp) != 1){
//			printf("file could not be opened");
//			exit(1);
//		}
//		if(fwrite(&(fiberp->v_xold), sizeof(double), 1, fp) != 1){
//			printf("file could not be opened");
//			exit(1);
//		}
//		if(fwrite(&(fiberp->v_yold), sizeof(double), 1, fp) != 1){
//			printf("file could not be opened");
//			exit(1);
//		}
//		if(fwrite(&(fiberp->angle), sizeof(double), 1, fp) != 1){
//			printf("file could not be opened");
//			exit(1);
//		}
//		if(fwrite(&(fiberp->angvel), sizeof(double), 1, fp) != 1){
//			printf("file could not be opened");
//			exit(1);
//		}
//		if(fwrite(&(fiberp->angvelold), sizeof(double), 1, fp) != 1){
//			printf("file could not be opened");
//			exit(1);
//		}
//		if(fwrite(&(fiberp->hx), sizeof(double), 1, fp) != 1){
//			printf("file could not be opened");
//			exit(1);
//		}
//		if(fwrite(&(fiberp->hy), sizeof(double), 1, fp) != 1){
//			printf("file could not be opened");
//			exit(1);
//		}
//		if(fwrite(&(fiberp->ha), sizeof(double), 1, fp) != 1){
//			printf("file could not be opened");
//			exit(1);
//		}
//	}

	fclose(fp);
}


void restartread(void){
	FILE* fp;
	char filename_restart[30];
	
	
	sprintf(filename_restart,"restart-%d.dat",n);
	
	
	if((fp = fopen(filename_restart, "rb"))==NULL){
		printf("file could not be opened\n");
		exit(1);
	}
	
	if(fread(&n, sizeof(int), 1, fp) != 1){
		printf("file could not be opened\n");
		exit(1);
	}
	n++;
	
	
	for(int i=0; i<(XCELL_NUM+6); i++) {
		for(int j=0; j<(YCELL_NUM+6); j++) {
			index0 = j + i*(YCELL_NUM+6);
			cellp= &CELL[index0];
			
			if(fread(&(cellp->v_x), sizeof(double), 1, fp) != 1){
				printf("file could not be opened");
				exit(1);
			}
			if(fread(&(cellp->v_y), sizeof(double), 1, fp) != 1){
				printf("file could not be opened");
				exit(1);
			}
			if(fread(&(cellp->press), sizeof(double), 1, fp) != 1){
				printf("file could not be opened");
				exit(1);
			}
			if(fread(&(cellp->phi), sizeof(double), 1, fp) != 1){
				printf("file could not be opened");
				exit(1);
			}
			if(fread(&(cellp->hx), sizeof(double), 1, fp) != 1){
				printf("file could not be opened");
				exit(1);
			}
			if(fread(&(cellp->hy), sizeof(double), 1, fp) != 1){
				printf("file could not be opened");
				exit(1);
			}
		}
	}
//	for(int i=0; i<FIBER_NUM; i++){
//		_fiber_type* fiberp = &FIBER[i];
//		
//		if(fread(&(fiberp->x), sizeof(double), 1, fp) != 1){
//			printf("file could not be opened");
//			exit(1);
//		}
//		if(fread(&(fiberp->y), sizeof(double), 1, fp) != 1){
//			printf("file could not be opened");
//			exit(1);
//		}
//		if(fread(&(fiberp->v_x), sizeof(double), 1, fp) != 1){
//			printf("file could not be opened");
//			exit(1);
//		}
//		if(fread(&(fiberp->v_y), sizeof(double), 1, fp) != 1){
//			printf("file could not be opened");
//			exit(1);
//		}
//		if(fread(&(fiberp->v_xold), sizeof(double), 1, fp) != 1){
//			printf("file could not be opened");
//			exit(1);
//		}
//		if(fread(&(fiberp->v_yold), sizeof(double), 1, fp) != 1){
//			printf("file could not be opened");
//			exit(1);
//		}
//		if(fread(&(fiberp->angle), sizeof(double), 1, fp) != 1){
//			printf("file could not be opened");
//			exit(1);
//		}
//		if(fread(&(fiberp->angvel), sizeof(double), 1, fp) != 1){
//			printf("file could not be opened");
//			exit(1);
//		}
//		if(fread(&(fiberp->angvelold), sizeof(double), 1, fp) != 1){
//			printf("file could not be opened");
//			exit(1);
//		}
//		if(fread(&(fiberp->hx), sizeof(double), 1, fp) != 1){
//			printf("file could not be opened");
//			exit(1);
//		}
//		if(fread(&(fiberp->hy), sizeof(double), 1, fp) != 1){
//			printf("file could not be opened");
//			exit(1);
//		}
//		if(fread(&(fiberp->ha), sizeof(double), 1, fp) != 1){
//			printf("file could not be opened");
//			exit(1);
//		}
//	}
	
	printf("Reading data for Restart\n");
	
	fclose(fp);
}


double interphi(int cellindex,int phiflag){
	double ctargetx,ctargety;
	int ctrgt_xp_yp;
	int ctrgt_xp;
	int ctrgt_yp;
	int ctrgt;
	int num=0;
	double array_deji[3][3];
	double b_deji[3];
	double x_deji[3];
	double target = 0.0;
	double targetwall = 0.0;
	double dh1;
	int xflag=0;
	int yflag=0;
	int cellindex_xp=cellindex+YCELL_NUM+6;
	int cellindex_yp=cellindex+1;
	_cell_type* cellp_in = &CELL[cellindex];
	_cell_type* cellp_inxp = &CELL[cellindex_xp];
	_cell_type* cellp_inyp = &CELL[cellindex_yp];
	int cindex,cindex_xp,cindex_yp,cindex_xp_yp;
	_flag_type* flagp_in = &FLAG[cellindex];
	double markx,marky,marknvec_x,marknvec_y;
	
	marknvec_x=(cellp_inxp->ls_x-cellp_in->ls_x)/dx;
	marknvec_y=(cellp_inyp->ls_y-cellp_in->ls_y)/dy;
	markx=cellp_in->x+marknvec_x*fabs(cellp_in->ls);
	marky=cellp_in->y+marknvec_y*fabs(cellp_in->ls);
	
	if(cellp_in->ls>=0.0){
		ctargetx = cellp_in->x;
		ctargety = cellp_in->y;
	}
	else {
		ctargetx = 2.0*markx-cellp_in->x;
		ctargety = 2.0*marky-cellp_in->y;
	}
	
	dh1=sqrt((markx-cellp_in->x)*(markx-cellp_in->x)+(marky-cellp_in->y)*(marky-cellp_in->y));
//		printf("dh1=%f \n",dh1);
	if(dh1<0.1*dx) dh1=0.1*dx;
	for(int hnum=0; hnum<100; hnum++){
		_flag_type* flagp_in = &FLAG[cellindex];
		
		if(cellp_in->ls>=0){
			ctargetx = cellp_in->x + marknvec_x * dh1 * ((double)hnum*0.25);
			ctargety = cellp_in->y + marknvec_y * dh1 * ((double)hnum*0.25);
		}
		else {
			ctargetx = 2.0*markx-cellp_in->x + marknvec_x * dh1 * ((double)hnum*0.25);
			ctargety = 2.0*marky-cellp_in->y + marknvec_y * dh1 * ((double)hnum*0.25);
		}
/*
#ifndef WALLX
		if(ctargetx>(dx*(double)XCELL_NUM)) {
			ctargetx-=dx*(double)XCELL_NUM;
			xflag=1;
		}
		else if(ctargetx<0.0) {
			ctargetx+=dx*(double)XCELL_NUM;
			xflag=2;
		}
		else xflag=0;
#endif
		
#ifndef WALLY
		if(ctargety>(dy*(double)YCELL_NUM)) {
			ctargety-=dy*(double)YCELL_NUM;
			yflag=1;
		}
		else if(ctargety<0.0) {
			ctargety+=dy*(double)YCELL_NUM;
			yflag=2;
		}
		else yflag=0;
#endif
*/		
		ctrgt = (int)((ctargetx+2.5*dx)/dx)*(YCELL_NUM+6)+(int)((ctargety+2.5*dy)/dy);	// The point at the bottom left of the moved point
		ctrgt_xp_yp = ctrgt+YCELL_NUM+6+1;
		ctrgt_xp = ctrgt + YCELL_NUM+6;
		ctrgt_yp = ctrgt + 1;
		
		_flag_type* flagp = &FLAG[ctrgt];
		_flag_type* flagp_xp = &FLAG[ctrgt+YCELL_NUM+6];
		_flag_type* flagp_yp = &FLAG[ctrgt+1];
		_flag_type* flagp_xp_yp = &FLAG[ctrgt+YCELL_NUM+6+1];
		_cell_type* ctrgtp = &CELL[ctrgt];
		_cell_type* ctrgtp_xp = &CELL[ctrgt_xp];
		_cell_type* ctrgtp_yp = &CELL[ctrgt_yp];
		_cell_type* ctrgtp_xp_yp = &CELL[ctrgt_xp_yp];
		
		if(hnum==99) printf("hnum>99 error x=%f y=%f flag=%d\n",cellp_in->x,cellp_in->y,phiflag);
		
		if(marknvec_x>0.92388){	//Vertical
			if((flagp->c%2)&&(flagp_yp->c%2)){
				hnum=100;
				
				array_deji[1][1]=ctrgtp->x;
				array_deji[1][2]=ctrgtp->y;
				array_deji[2][1]=ctrgtp_yp->x;
				array_deji[2][2]=ctrgtp_yp->y;
				if((phiflag==0)||(phiflag==1)){
					b_deji[1]=ctrgtp->phi;
					b_deji[2]=ctrgtp_yp->phi;
				}
				else if(phiflag==2){
					b_deji[1]=phid2[ctrgt];
					b_deji[2]=phid2[ctrgt_yp];
				}
				else if(phiflag==3){
					b_deji[1]=ddrho[ctrgt];
					b_deji[2]=ddrho[ctrgt_yp];
				}
			}
		}
		else if(marknvec_x<-0.92388){
			if((flagp_xp->c%2)&&(flagp_xp_yp->c%2)){
				hnum=100;
				
				array_deji[1][1]=ctrgtp_xp->x;
				array_deji[1][2]=ctrgtp_xp->y;
				array_deji[2][1]=ctrgtp_xp_yp->x;
				array_deji[2][2]=ctrgtp_xp_yp->y;
				if((phiflag==0)||(phiflag==1)){
					b_deji[1]=ctrgtp_xp->phi;
					b_deji[2]=ctrgtp_xp_yp->phi;
				}
				else if(phiflag==2){
					b_deji[1]=phid2[ctrgt_xp];
					b_deji[2]=phid2[ctrgt_xp_yp];
				}
				else if(phiflag==3){
					b_deji[1]=ddrho[ctrgt_xp];
					b_deji[2]=ddrho[ctrgt_xp_yp];
				}
			}
		}
		else if(marknvec_y>0.92388){	// Horizontal
			if((flagp->c%2)&&(flagp_xp->c%2)){
				hnum=100;
				
				array_deji[1][1]=ctrgtp->x;
				array_deji[1][2]=ctrgtp->y;
				array_deji[2][1]=ctrgtp_xp->x;
				array_deji[2][2]=ctrgtp_xp->y;
				if((phiflag==0)||(phiflag==1)){
					b_deji[1]=ctrgtp->phi;
					b_deji[2]=ctrgtp_xp->phi;
				}
				else if(phiflag==2){
					b_deji[1]=phid2[ctrgt];
					b_deji[2]=phid2[ctrgt_xp];
				}
				else if(phiflag==3){
					b_deji[1]=ddrho[ctrgt];
					b_deji[2]=ddrho[ctrgt_xp];
				}
			}
		}
		else if(marknvec_y<-0.92388){	// Horizontal
			if((flagp_yp->c%2)&&(flagp_xp_yp->c%2)){
				hnum=100;
				
				array_deji[1][1]=ctrgtp_yp->x;
				array_deji[1][2]=ctrgtp_yp->y;
				array_deji[2][1]=ctrgtp_xp_yp->x;
				array_deji[2][2]=ctrgtp_xp_yp->y;
				if((phiflag==0)||(phiflag==1)){
					b_deji[1]=ctrgtp_yp->phi;
					b_deji[2]=ctrgtp_xp_yp->phi;
				}
				else if(phiflag==2){
					b_deji[1]=phid2[ctrgt_yp];
					b_deji[2]=phid2[ctrgt_xp_yp];
				}
				else if(phiflag==3){
					b_deji[1]=ddrho[ctrgt_yp];
					b_deji[2]=ddrho[ctrgt_xp_yp];
				}
			}
		}
		else if((marknvec_x*marknvec_y)>0.0){	// Diagonal (negative slope)
			if((flagp_xp->c%2)&&(flagp_yp->c%2)){
				hnum=100;
				
				array_deji[1][1]=ctrgtp_xp->x;
				array_deji[1][2]=ctrgtp_xp->y;
				array_deji[2][1]=ctrgtp_yp->x;
				array_deji[2][2]=ctrgtp_yp->y;
				if((phiflag==0)||(phiflag==1)){
					b_deji[1]=ctrgtp_xp->phi;
					b_deji[2]=ctrgtp_yp->phi;
				}
				else if(phiflag==2){
					b_deji[1]=phid2[ctrgt_xp];
					b_deji[2]=phid2[ctrgt_yp];
				}
				else if(phiflag==3){
					b_deji[1]=ddrho[ctrgt_xp];
					b_deji[2]=ddrho[ctrgt_yp];
				}
			}
		}
		else {	// Diagonal (positive slope)
			if((flagp->c%2)&&(flagp_xp_yp->c%2)){
				hnum=100;
				
				array_deji[1][1]=ctrgtp->x;
				array_deji[1][2]=ctrgtp->y;
				array_deji[2][1]=ctrgtp_xp_yp->x;
				array_deji[2][2]=ctrgtp_xp_yp->y;
				b_deji[1]=ctrgtp->phi;
				b_deji[2]=ctrgtp_xp_yp->phi;
				if((phiflag==0)||(phiflag==1)){
					b_deji[1]=ctrgtp->phi;
					b_deji[2]=ctrgtp_xp_yp->phi;
				}
				else if(phiflag==2){
					b_deji[1]=phid2[ctrgt];
					b_deji[2]=phid2[ctrgt_xp_yp];
				}
				else if(phiflag==3){
					b_deji[1]=ddrho[ctrgt];
					b_deji[2]=ddrho[ctrgt_xp_yp];
				}
			}
		}
	}
	
	if(xflag==1){
		array_deji[1][1]+=dx*(double)XCELL_NUM;
		array_deji[2][1]+=dx*(double)XCELL_NUM;
	}
	else if(xflag==2){
		array_deji[1][1]-=dx*(double)XCELL_NUM;
		array_deji[2][1]-=dx*(double)XCELL_NUM;
	}
	
	if(yflag==1){
		array_deji[1][2]+=dy*(double)YCELL_NUM;
		array_deji[2][2]+=dy*(double)YCELL_NUM;
	}
	else if(yflag==2){
		array_deji[1][2]-=dy*(double)YCELL_NUM;
		array_deji[2][2]-=dy*(double)YCELL_NUM;
	}
	
	array_deji[0][0]=0.0;
	array_deji[1][0]=1.0;
	array_deji[2][0]=1.0;
	array_deji[0][1]=marknvec_x;
	array_deji[0][2]=marknvec_y;
	b_deji[0]=0.0;
	
	lu(array_deji,b_deji,x_deji);	// LU decomposition
	
	
	target=x_deji[0]+x_deji[1]*cellp_in->x+x_deji[2]*cellp_in->y;
	if(phiflag==0){
		if((target<=PHIL)&&(target>=PHIG)){
			b_deji[0]=wetphi[flagp_in->fiber];
			
			lu(array_deji,b_deji,x_deji);
			
			target=x_deji[0]+x_deji[1]*cellp_in->x+x_deji[2]*cellp_in->y;
		}
		
		if(target<=PHImin) target=PHImin;
		if(target>=PHImax) target=PHImax;
	}
	
	return target;
}
		
		
		
void solidmove(void){
	for(int i=0; i<FIBER_NUM; i++){
		_fiber_type* fiberp=&FIBER[i];
		
		fiberp->x += 0.5*dt*(fiberp->v_xold+fiberp->v_x);
		fiberp->y += 0.5*dt*(fiberp->v_yold+fiberp->v_y);
		fiberp->v_xold=fiberp->v_x;
		fiberp->v_yold=fiberp->v_y;
		
//		printf("\n n=%d u=%f v=%f x=%f y=%f\n",n,fiberp->v_x,fiberp->v_y,fiberp->x,fiberp->y);
		
		fiberp->angle += 0.5*dt*(fiberp->angvelold+fiberp->angvel);
		fiberp->angvelold=fiberp->angvel;
	}
}
	

double ljpower(double sigma,double distance){
	double sigma6=sigma*sigma*sigma*sigma*sigma*sigma;
	double dis6=distance*distance*distance*distance*distance*distance;
	
	
	return 4.0*LJE*(12.0*sigma6*sigma6/dis6/distance-6.0*sigma6/dis6/distance);
}
	
	
void solidlagpoint(_fiber_type* fiberp,int fibernumber){
	double fiberlagx,fiberlagy;
	double fiberlagv_x,fiberlagv_y;
	double deltas;
	int lagnum = 5000;
	double deltar = 1.0;
	double *drhodx, *drhody,*dvxdx,*dvxdy,*dvydx,*dvydy;
	double lagv_x,lagv_y,lagdrhodx,lagdrhody,lagddrho,lagmu,lagrho,lagphi;
	double lagdvxdx,lagdvxdy,lagdvydx,lagdvydy,lagpress;
	double lagDxx,lagDxy,lagDyy;
	
	drhodx = new double [CELL_SIZE];
	memset(drhodx,0,CELL_SIZE*sizeof(double));
	drhody = new double [CELL_SIZE];
	memset(drhody,0,CELL_SIZE*sizeof(double));
	dvxdx = new double [CELL_SIZE];
	memset(dvxdx,0,CELL_SIZE*sizeof(double));
	dvxdy = new double [CELL_SIZE];
	memset(dvxdy,0,CELL_SIZE*sizeof(double));
	dvydx = new double [CELL_SIZE];
	memset(dvydx,0,CELL_SIZE*sizeof(double));
	dvydy = new double [CELL_SIZE];
	memset(dvydy,0,CELL_SIZE*sizeof(double));
	
	for(int i=3; i<(XCELL_NUM+3); i++){ // For periodic boundaries, it's 3
		for(int j=3; j<(YCELL_NUM+3); j++){
			indexnsdefine(i,j);
			
			int index_xp_ypp = (j+2) + (i+1)*(YCELL_NUM+6);
			int index_xp_yp  = (j+1) + (i+1)*(YCELL_NUM+6);
			int index_xp_ymm = (j-2) + (i+1)*(YCELL_NUM+6);
			int index_xpp_yp = (j+1) + (i+2)*(YCELL_NUM+6);
			int index_xmm_yp = (j+1) + (i-2)*(YCELL_NUM+6);
			
			_cell_type* cellp_xp_ypp = &CELL[index_xp_ypp];
			_cell_type* cellp_xp_yp  = &CELL[index_xp_yp];
			_cell_type* cellp_xp_ymm = &CELL[index_xp_ymm];
			_cell_type* cellp_xpp_yp = &CELL[index_xpp_yp];
			_cell_type* cellp_xmm_yp = &CELL[index_xmm_yp];
			
			drhodx[index0] = (-cellp_xpp->rho+8.0*(cellp_xp->rho-cellp_xm->rho)+cellp_xmm->rho)/12.0/dx;
			drhody[index0] = (-cellp_ypp->rho+8.0*(cellp_yp->rho-cellp_ym->rho)+cellp_ymm->rho)/12.0/dy;
			dvxdx[index0]  = (-cellp_xpp->v_x+27.0*(cellp_xp->v_x-cellp->v_x)+cellp_xm->v_x)/24.0/dx;
			dvydy[index0]  = (-cellp_ypp->v_y+27.0*(cellp_yp->v_y-cellp->v_y)+cellp_ym->v_y)/24.0/dy;
			dvxdy[index0]  = 0.5*(-cellp_ypp->v_x+8.0*(cellp_yp->v_x-cellp_ym->v_x)+cellp_ymm->v_x
							-cellp_xp_ypp->v_x+8.0*(cellp_xp_yp->v_x-cellp_xp_ym->v_x)+cellp_xp_ymm->v_x)/12.0/dy;
			dvydx[index0]  = 0.5*(-cellp_xpp->v_y+8.0*(cellp_xp->v_y-cellp_xm->v_y)+cellp_xmm->v_y
							-cellp_xpp_yp->v_y+8.0*(cellp_xp_yp->v_y-cellp_xm_yp->v_y)+cellp_xmm_yp->v_y)/12.0/dx;
		}
	}
	
#ifndef WALLY
	for(int i=0; i<(XCELL_NUM+6); i++){											//E
		indexbndy(i);
			
		drhodx[index_y0]   = drhodx[index_ycp0];
		drhodx[index_y1]   = drhodx[index_ycp1];
		drhodx[index_y2]   = drhodx[index_ycp2];
		drhodx[index_ycp3] = drhodx[index_y3];
		drhodx[index_ycp4] = drhodx[index_y4];
		drhodx[index_ycp5] = drhodx[index_y5];
		drhody[index_y0]   = drhody[index_ycp0];
		drhody[index_y1]   = drhody[index_ycp1];
		drhody[index_y2]   = drhody[index_ycp2];
		drhody[index_ycp3] = drhody[index_y3];
		drhody[index_ycp4] = drhody[index_y4];
		drhody[index_ycp5] = drhody[index_y5];
		dvxdx[index_y0]   = dvxdx[index_ycp0];
		dvxdx[index_y1]   = dvxdx[index_ycp1];
		dvxdx[index_y2]   = dvxdx[index_ycp2];
		dvxdx[index_ycp3] = dvxdx[index_y3];
		dvxdx[index_ycp4] = dvxdx[index_y4];
		dvxdx[index_ycp5] = dvxdx[index_y5];
		dvydy[index_y0]   = dvydy[index_ycp0];
		dvydy[index_y1]   = dvydy[index_ycp1];
		dvydy[index_y2]   = dvydy[index_ycp2];
		dvydy[index_ycp3] = dvydy[index_y3];
		dvydy[index_ycp4] = dvydy[index_y4];
		dvydy[index_ycp5] = dvydy[index_y5];
		dvxdy[index_y0]   = dvxdy[index_ycp0];
		dvxdy[index_y1]   = dvxdy[index_ycp1];
		dvxdy[index_y2]   = dvxdy[index_ycp2];
		dvxdy[index_ycp3] = dvxdy[index_y3];
		dvxdy[index_ycp4] = dvxdy[index_y4];
		dvxdy[index_ycp5] = dvxdy[index_y5];
		dvydx[index_y0]   = dvydx[index_ycp0];
		dvydx[index_y1]   = dvydx[index_ycp1];
		dvydx[index_y2]   = dvydx[index_ycp2];
		dvydx[index_ycp3] = dvydx[index_y3];
		dvydx[index_ycp4] = dvydx[index_y4];
		dvydx[index_ycp5] = dvydx[index_y5];
	}
#endif
#ifndef WALLX
	for(int j=0; j<(YCELL_NUM+6); j++){
		indexbndx(j);
		
		drhodx[index_x0]   = drhodx[index_xcp0];
		drhodx[index_x1]   = drhodx[index_xcp1];
		drhodx[index_x2]   = drhodx[index_xcp2];
		drhodx[index_xcp3] = drhodx[index_x3];
		drhodx[index_xcp4] = drhodx[index_x4];
		drhodx[index_xcp5] = drhodx[index_x5];
		drhody[index_x0]   = drhody[index_xcp0];
		drhody[index_x1]   = drhody[index_xcp1];
		drhody[index_x2]   = drhody[index_xcp2];
		drhody[index_xcp3] = drhody[index_x3];
		drhody[index_xcp4] = drhody[index_x4];
		drhody[index_xcp5] = drhody[index_x5];
		dvxdx[index_x0]   = dvxdx[index_xcp0];
		dvxdx[index_x1]   = dvxdx[index_xcp1];
		dvxdx[index_x2]   = dvxdx[index_xcp2];
		dvxdx[index_xcp3] = dvxdx[index_x3];
		dvxdx[index_xcp4] = dvxdx[index_x4];
		dvxdx[index_xcp5] = dvxdx[index_x5];
		dvydy[index_x0]   = dvydy[index_xcp0];
		dvydy[index_x1]   = dvydy[index_xcp1];
		dvydy[index_x2]   = dvydy[index_xcp2];
		dvydy[index_xcp3] = dvydy[index_x3];
		dvydy[index_xcp4] = dvydy[index_x4];
		dvydy[index_xcp5] = dvydy[index_x5];
		dvxdy[index_x0]   = dvxdy[index_xcp0];
		dvxdy[index_x1]   = dvxdy[index_xcp1];
		dvxdy[index_x2]   = dvxdy[index_xcp2];
		dvxdy[index_xcp3] = dvxdy[index_x3];
		dvxdy[index_xcp4] = dvxdy[index_x4];
		dvxdy[index_xcp5] = dvxdy[index_x5];
		dvydx[index_x0]   = dvydx[index_xcp0];
		dvydx[index_x1]   = dvydx[index_xcp1];
		dvydx[index_x2]   = dvydx[index_xcp2];
		dvydx[index_xcp3] = dvydx[index_x3];
		dvydx[index_xcp4] = dvydx[index_x4];
		dvydx[index_xcp5] = dvydx[index_x5];
	}
#endif

#ifdef WALLY
	for(int i=0; i<(XCELL_NUM+6); i++){
		indexbndy(i);
		
		drhodx[index_y0]   = -drhodx[index_y5];
		drhodx[index_y1]   = -drhodx[index_y4];
		drhodx[index_y2]   = -drhodx[index_y3];
		drhodx[index_ycp3] = -drhodx[index_ycp2];
		drhodx[index_ycp4] = -drhodx[index_ycp1];
		drhodx[index_ycp5] = -drhodx[index_ycp0];
		drhody[index_y0]   = -drhody[index_y5];
		drhody[index_y1]   = -drhody[index_y4];
		drhody[index_y2]   = -drhody[index_y3];
		drhody[index_ycp3] = -drhody[index_ycp2];
		drhody[index_ycp4] = -drhody[index_ycp1];
		drhody[index_ycp5] = -drhody[index_ycp0];
		dvxdx[index_y0]   = dvxdx[index_ycp0];
		dvxdx[index_y1]   = dvxdx[index_ycp1];
		dvxdx[index_y2]   = dvxdx[index_ycp2];
		dvxdx[index_ycp3] = dvxdx[index_y3];
		dvxdx[index_ycp4] = dvxdx[index_y4];
		dvxdx[index_ycp5] = dvxdx[index_y5];
		dvydy[index_y0]   = dvydy[index_ycp0];
		dvydy[index_y1]   = dvydy[index_ycp1];
		dvydy[index_y2]   = dvydy[index_ycp2];
		dvydy[index_ycp3] = dvydy[index_y3];
		dvydy[index_ycp4] = dvydy[index_y4];
		dvydy[index_ycp5] = dvydy[index_y5];
		dvxdy[index_y0]   = dvxdy[index_ycp0];
		dvxdy[index_y1]   = dvxdy[index_ycp1];
		dvxdy[index_y2]   = dvxdy[index_ycp2];
		dvxdy[index_ycp3] = dvxdy[index_y3];
		dvxdy[index_ycp4] = dvxdy[index_y4];
		dvxdy[index_ycp5] = dvxdy[index_y5];
		dvydx[index_y0]   = dvydx[index_ycp0];
		dvydx[index_y1]   = dvydx[index_ycp1];
		dvydx[index_y2]   = dvydx[index_ycp2];
		dvydx[index_ycp3] = dvydx[index_y3];
		dvydx[index_ycp4] = dvydx[index_y4];
		dvydx[index_ycp5] = dvydx[index_y5];
	}
#endif
#ifdef WALLX
	for(int j=0; j<(YCELL_NUM+6); j++){
		indexbndx(j);
		
		drhodx[index_x0]   = -drhodx[index_x5];
		drhodx[index_x1]   = -drhodx[index_x4];
		drhodx[index_x2]   = -drhodx[index_x3];
		drhodx[index_xcp3] = -drhodx[index_xcp2];
		drhodx[index_xcp4] = -drhodx[index_xcp1];
		drhodx[index_xcp5] = -drhodx[index_xcp0];
		drhody[index_x0]   = -drhody[index_x5];
		drhody[index_x1]   = -drhody[index_x4];
		drhody[index_x2]   = -drhody[index_x3];
		drhody[index_xcp3] = -drhody[index_xcp2];
		drhody[index_xcp4] = -drhody[index_xcp1];
		drhody[index_xcp5] = -drhody[index_xcp0];
		dvxdx[index_x0]   = dvxdx[index_xcp0];
		dvxdx[index_x1]   = dvxdx[index_xcp1];
		dvxdx[index_x2]   = dvxdx[index_xcp2];
		dvxdx[index_xcp3] = dvxdx[index_x3];
		dvxdx[index_xcp4] = dvxdx[index_x4];
		dvxdx[index_xcp5] = dvxdx[index_x5];
		dvydy[index_x0]   = dvydy[index_xcp0];
		dvydy[index_x1]   = dvydy[index_xcp1];
		dvydy[index_x2]   = dvydy[index_xcp2];
		dvydy[index_xcp3] = dvydy[index_x3];
		dvydy[index_xcp4] = dvydy[index_x4];
		dvydy[index_xcp5] = dvydy[index_x5];
		dvxdy[index_x0]   = dvxdy[index_xcp0];
		dvxdy[index_x1]   = dvxdy[index_xcp1];
		dvxdy[index_x2]   = dvxdy[index_xcp2];
		dvxdy[index_xcp3] = dvxdy[index_x3];
		dvxdy[index_xcp4] = dvxdy[index_x4];
		dvxdy[index_xcp5] = dvxdy[index_x5];
		dvydx[index_x0]   = dvydx[index_xcp0];
		dvydx[index_x1]   = dvydx[index_xcp1];
		dvydx[index_x2]   = dvydx[index_xcp2];
		dvydx[index_xcp3] = dvydx[index_x3];
		dvydx[index_xcp4] = dvydx[index_x4];
		dvydx[index_xcp5] = dvydx[index_x5];
	}
#endif
	
	deltas = 2.0*M_PI*fiberp->radius/(double)lagnum;
	
	for (int lagindex=0; lagindex<lagnum; lagindex++){
		fiberlagx = fiberp->x+(fiberp->radius+deltar)*cos(2.0*M_PI/(double)lagnum*(double)lagindex);
		fiberlagy = fiberp->y+(fiberp->radius+deltar)*sin(2.0*M_PI/(double)lagnum*(double)lagindex);
		fiberlagv_x = fiberp->v_x-(fiberp->radius+deltar)*fiberp->angvel*sin(2.0*M_PI/(double)lagnum*(double)lagindex);
		fiberlagv_y = fiberp->v_y+(fiberp->radius+deltar)*fiberp->angvel*cos(2.0*M_PI/(double)lagnum*(double)lagindex);
		
		double nvec_x = cos(2.0*M_PI/(double)lagnum*(double)lagindex);
		double nvec_y = sin(2.0*M_PI/(double)lagnum*(double)lagindex);
		
		int ctrgt = (int)((fiberlagx+2.5*dx)/dx)*(YCELL_NUM+6)+(int)((fiberlagy+2.5*dy)/dy);	// The point at the bottom left of the moved point
		int ctrgt_xp_yp = ctrgt+YCELL_NUM+6+1;
		int ctrgt_xp = ctrgt + YCELL_NUM+6;
		int ctrgt_yp = ctrgt + 1;
		
		_cell_type* ctrgtp = &CELL[ctrgt];
		_cell_type* ctrgtp_xp = &CELL[ctrgt_xp];
		_cell_type* ctrgtp_yp = &CELL[ctrgt_yp];
		_cell_type* ctrgtp_xp_yp = &CELL[ctrgt_xp_yp];
		
		
//		lagrho = (ctrgtp->rho*(ctrgtp_xp->x-fiberlagx)*(ctrgtp_yp->y-fiberlagy)/dx/dy + ctrgtp_xp->rho*(fiberlagx-ctrgtp->x)*(ctrgtp_yp->y-fiberlagy)/dx/dy
//					+ ctrgtp_yp->rho*(ctrgtp_xp->x-fiberlagx)*(fiberlagy-ctrgtp->y)/dx/dy + ctrgtp_xp_yp->rho*(fiberlagx-ctrgtp->x)*(fiberlagy-ctrgtp->y)/dx/dy);
		
		lagphi = (ctrgtp->phi*(ctrgtp_xp->x-fiberlagx)*(ctrgtp_yp->y-fiberlagy)/dx/dy + ctrgtp_xp->phi*(fiberlagx-ctrgtp->x)*(ctrgtp_yp->y-fiberlagy)/dx/dy
					+ ctrgtp_yp->phi*(ctrgtp_xp->x-fiberlagx)*(fiberlagy-ctrgtp->y)/dx/dy + ctrgtp_xp_yp->phi*(fiberlagx-ctrgtp->x)*(fiberlagy-ctrgtp->y)/dx/dy);
		
		if(lagphi>PHImax) lagphi=PHImax;
		if(lagphi<PHImin) lagphi=PHImin;
		
		if(lagphi<=PHIG) lagrho=RHOG;
		else if (lagphi>=PHIL) lagrho=RHOL;
		else lagrho=(RHOL+RHOG)*0.5+(RHOL-RHOG)*0.5*sin((lagphi-(PHIL+PHIG)*0.5)/(PHIL-PHIG)*M_PI);
		
		lagpress = (ctrgtp->press*(ctrgtp_xp->x-fiberlagx)*(ctrgtp_yp->y-fiberlagy)/dx/dy + ctrgtp_xp->press*(fiberlagx-ctrgtp->x)*(ctrgtp_yp->y-fiberlagy)/dx/dy
					+ ctrgtp_yp->press*(ctrgtp_xp->x-fiberlagx)*(fiberlagy-ctrgtp->y)/dx/dy + ctrgtp_xp_yp->press*(fiberlagx-ctrgtp->x)*(fiberlagy-ctrgtp->y)/dx/dy);
		
		lagddrho = (ddrho[ctrgt]*(ctrgtp_xp->x-fiberlagx)*(ctrgtp_yp->y-fiberlagy)/dx/dy + ddrho[ctrgt_xp]*(fiberlagx-ctrgtp->x)*(ctrgtp_yp->y-fiberlagy)/dx/dy
					+ ddrho[ctrgt_yp]*(ctrgtp_xp->x-fiberlagx)*(fiberlagy-ctrgtp->y)/dx/dy + ddrho[ctrgt_xp_yp]*(fiberlagx-ctrgtp->x)*(fiberlagy-ctrgtp->y)/dx/dy);
		
		lagdrhodx = (drhodx[ctrgt]*(ctrgtp_xp->x-fiberlagx)*(ctrgtp_yp->y-fiberlagy)/dx/dy + drhodx[ctrgt_xp]*(fiberlagx-ctrgtp->x)*(ctrgtp_yp->y-fiberlagy)/dx/dy
					+ drhodx[ctrgt_yp]*(ctrgtp_xp->x-fiberlagx)*(fiberlagy-ctrgtp->y)/dx/dy + drhodx[ctrgt_xp_yp]*(fiberlagx-ctrgtp->x)*(fiberlagy-ctrgtp->y)/dx/dy);
		
		lagdrhody = (drhody[ctrgt]*(ctrgtp_xp->x-fiberlagx)*(ctrgtp_yp->y-fiberlagy)/dx/dy + drhody[ctrgt_xp]*(fiberlagx-ctrgtp->x)*(ctrgtp_yp->y-fiberlagy)/dx/dy
					+ drhody[ctrgt_yp]*(ctrgtp_xp->x-fiberlagx)*(fiberlagy-ctrgtp->y)/dx/dy + drhody[ctrgt_xp_yp]*(fiberlagx-ctrgtp->x)*(fiberlagy-ctrgtp->y)/dx/dy);
		lagdvxdx = (dvxdx[ctrgt]*(ctrgtp_xp->x-fiberlagx)*(ctrgtp_yp->y-fiberlagy)/dx/dy + dvxdx[ctrgt_xp]*(fiberlagx-ctrgtp->x)*(ctrgtp_yp->y-fiberlagy)/dx/dy
					+ dvxdx[ctrgt_yp]*(ctrgtp_xp->x-fiberlagx)*(fiberlagy-ctrgtp->y)/dx/dy + dvxdx[ctrgt_xp_yp]*(fiberlagx-ctrgtp->x)*(fiberlagy-ctrgtp->y)/dx/dy);
		lagdvydy = (dvydy[ctrgt]*(ctrgtp_xp->x-fiberlagx)*(ctrgtp_yp->y-fiberlagy)/dx/dy + dvydy[ctrgt_xp]*(fiberlagx-ctrgtp->x)*(ctrgtp_yp->y-fiberlagy)/dx/dy
					+ dvydy[ctrgt_yp]*(ctrgtp_xp->x-fiberlagx)*(fiberlagy-ctrgtp->y)/dx/dy + dvydy[ctrgt_xp_yp]*(fiberlagx-ctrgtp->x)*(fiberlagy-ctrgtp->y)/dx/dy);
		lagdvxdy = (dvxdy[ctrgt]*(ctrgtp_xp->x-fiberlagx)*(ctrgtp_yp->y-fiberlagy)/dx/dy + dvxdy[ctrgt_xp]*(fiberlagx-ctrgtp->x)*(ctrgtp_yp->y-fiberlagy)/dx/dy
					+ dvxdy[ctrgt_yp]*(ctrgtp_xp->x-fiberlagx)*(fiberlagy-ctrgtp->y)/dx/dy + dvxdy[ctrgt_xp_yp]*(fiberlagx-ctrgtp->x)*(fiberlagy-ctrgtp->y)/dx/dy);
		lagdvydx = (dvydx[ctrgt]*(ctrgtp_xp->x-fiberlagx)*(ctrgtp_yp->y-fiberlagy)/dx/dy + dvydx[ctrgt_xp]*(fiberlagx-ctrgtp->x)*(ctrgtp_yp->y-fiberlagy)/dx/dy
					+ dvydx[ctrgt_yp]*(ctrgtp_xp->x-fiberlagx)*(fiberlagy-ctrgtp->y)/dx/dy + dvydx[ctrgt_xp_yp]*(fiberlagx-ctrgtp->x)*(fiberlagy-ctrgtp->y)/dx/dy);
		
		lagmu = MUG+(MUL-MUG)/(RHOL-RHOG)*(lagrho-RHOG);
		
		ctrgt = (int)((fiberlagx+3.0*dx)/dx)*(YCELL_NUM+6)+(int)((fiberlagy+2.5*dy)/dy);	// The point at the bottom left of the moved point
		ctrgt_xp_yp = ctrgt+YCELL_NUM+6+1;
		ctrgt_xp = ctrgt + YCELL_NUM+6;
		ctrgt_yp = ctrgt + 1;
		
		ctrgtp = &CELL[ctrgt];
		ctrgtp_xp = &CELL[ctrgt_xp];
		ctrgtp_yp = &CELL[ctrgt_yp];
		ctrgtp_xp_yp = &CELL[ctrgt_xp_yp];
		
		lagv_x = (ctrgtp->v_x*(ctrgtp_xp->x-0.5*dx-fiberlagx)*(ctrgtp_yp->y-fiberlagy)/dx/dy + ctrgtp_xp->v_x*(fiberlagx-ctrgtp->x+0.5*dx)*(ctrgtp_yp->y-fiberlagy)/dx/dy
					+ ctrgtp_yp->v_x*(ctrgtp_xp->x-0.5*dx-fiberlagx)*(fiberlagy-ctrgtp->y)/dx/dy + ctrgtp_xp_yp->v_x*(fiberlagx-ctrgtp->x+0.5*dx)*(fiberlagy-ctrgtp->y)/dx/dy);
		
		ctrgt = (int)((fiberlagx+2.5*dx)/dx)*(YCELL_NUM+6)+(int)((fiberlagy+3.0*dy)/dy);	// The point at the bottom left of the moved point
		ctrgt_xp_yp = ctrgt+YCELL_NUM+6+1;
		ctrgt_xp = ctrgt + YCELL_NUM+6;
		ctrgt_yp = ctrgt + 1;
		
		ctrgtp = &CELL[ctrgt];
		ctrgtp_xp = &CELL[ctrgt_xp];
		ctrgtp_yp = &CELL[ctrgt_yp];
		ctrgtp_xp_yp = &CELL[ctrgt_xp_yp];
		
		lagv_y = (ctrgtp->v_y*(ctrgtp_xp->x-fiberlagx)*(ctrgtp_yp->y-0.5*dy-fiberlagy)/dx/dy + ctrgtp_xp->v_y*(fiberlagx-ctrgtp->x)*(ctrgtp_yp->y-0.5*dy-fiberlagy)/dx/dy
					+ ctrgtp_yp->v_y*(ctrgtp_xp->x-fiberlagx)*(fiberlagy-ctrgtp->y+0.5*dy)/dx/dy + ctrgtp_xp_yp->v_y*(fiberlagx-ctrgtp->x)*(fiberlagy-ctrgtp->y+0.5*dy)/dx/dy);
		
		fiberfx[fibernumber] += deltas*(-lagpress*nvec_x+kappa1*((lagrho*lagddrho+0.5*(lagdrhodx*lagdrhodx+lagdrhody*lagdrhody)-lagdrhodx*lagdrhodx)*nvec_x-lagdrhodx*lagdrhody*nvec_y));
		fiberfy[fibernumber] += deltas*(-lagpress*nvec_y+kappa1*((lagrho*lagddrho+0.5*(lagdrhodx*lagdrhodx+lagdrhody*lagdrhody)-lagdrhody*lagdrhody)*nvec_y-lagdrhodx*lagdrhody*nvec_x));
		
		fiberfx[fibernumber] -= deltas*lagrho*lagv_x*((lagv_x-fiberlagv_x)*nvec_x+(lagv_y-fiberlagv_y)*nvec_y);
		fiberfy[fibernumber] -= deltas*lagrho*lagv_y*((lagv_x-fiberlagv_x)*nvec_x+(lagv_y-fiberlagv_y)*nvec_y);
		
		lagDxx = 0.5*(lagdvxdx + lagdvxdx);
		lagDxy = 0.5*(lagdvxdy + lagdvydx);
		lagDyy = 0.5*(lagdvydy + lagdvydy);
		
		fiberfx[fibernumber] += 2.0*deltas*lagmu*(lagDxx*nvec_x+lagDxy*nvec_y);
		fiberfy[fibernumber] += 2.0*deltas*lagmu*(lagDxy*nvec_x+lagDyy*nvec_y);
		
		fibern[fibernumber]-=(deltas*(-lagpress*nvec_x+kappa1*((lagrho*lagddrho+0.5*(lagdrhodx*lagdrhodx+lagdrhody*lagdrhody)-lagdrhodx*lagdrhodx)*nvec_x-lagdrhodx*lagdrhody*nvec_y))
								- deltas*lagrho*lagv_x*((lagv_x-fiberlagv_x)*nvec_x+(lagv_y-fiberlagv_y)*nvec_y)
								+ 2.0*deltas*lagmu*(lagDxx*nvec_x+lagDxy*nvec_y))*(fiberlagy-fiberp->y);
		fibern[fibernumber]+=(deltas*(-lagpress*nvec_y+kappa1*((lagrho*lagddrho+0.5*(lagdrhodx*lagdrhodx+lagdrhody*lagdrhody)-lagdrhody*lagdrhody)*nvec_y-lagdrhodx*lagdrhody*nvec_x))
								- deltas*lagrho*lagv_y*((lagv_x-fiberlagv_x)*nvec_x+(lagv_y-fiberlagv_y)*nvec_y)
								+ 2.0*deltas*lagmu*(lagDxy*nvec_x+lagDyy*nvec_y))*(fiberlagx-fiberp->x);
	}
	
	delete [] drhodx;
	delete [] drhody;
	delete [] dvxdx;
	delete [] dvxdy;
	delete [] dvydx;
	delete [] dvydy;
}

void solidlagpointsecond(_fiber_type* fiberp,int fibernumber){
	double fiberlagx,fiberlagy;
	double deltas;
	int lagnum = 2000;
	double deltar = 1.0;
	double lagv_x,lagv_y,lagrho,lagphi,lagoldv_x,lagoldv_y;
	
	deltas = 2.0*M_PI*fiberp->radius/(double)lagnum;
	
	for (int lagindex=0; lagindex<lagnum; lagindex++){
		fiberlagx = fiberp->x+(fiberp->radius+deltar*0.5)*cos(2.0*M_PI/(double)lagnum*(double)lagindex);
		fiberlagy = fiberp->y+(fiberp->radius+deltar*0.5)*sin(2.0*M_PI/(double)lagnum*(double)lagindex);
		
		int ctrgt = (int)((fiberlagx+2.5*dx)/dx)*(YCELL_NUM+6)+(int)((fiberlagy+2.5*dy)/dy);	// The point at the bottom left of the moved point
		int ctrgt_xp_yp = ctrgt+YCELL_NUM+6+1;
		int ctrgt_xp = ctrgt + YCELL_NUM+6;
		int ctrgt_yp = ctrgt + 1;
		
		_cell_type* ctrgtp = &CELL[ctrgt];
		_cell_type* ctrgtp_xp = &CELL[ctrgt_xp];
		_cell_type* ctrgtp_yp = &CELL[ctrgt_yp];
		_cell_type* ctrgtp_xp_yp = &CELL[ctrgt_xp_yp];
		
		
//		lagrho = (ctrgtp->rho*(ctrgtp_xp->x-fiberlagx)*(ctrgtp_yp->y-fiberlagy)/dx/dy + ctrgtp_xp->rho*(fiberlagx-ctrgtp->x)*(ctrgtp_yp->y-fiberlagy)/dx/dy
//					+ ctrgtp_yp->rho*(ctrgtp_xp->x-fiberlagx)*(fiberlagy-ctrgtp->y)/dx/dy + ctrgtp_xp_yp->rho*(fiberlagx-ctrgtp->x)*(fiberlagy-ctrgtp->y)/dx/dy);
		
		lagphi = (ctrgtp->phi*(ctrgtp_xp->x-fiberlagx)*(ctrgtp_yp->y-fiberlagy)/dx/dy + ctrgtp_xp->phi*(fiberlagx-ctrgtp->x)*(ctrgtp_yp->y-fiberlagy)/dx/dy
					+ ctrgtp_yp->phi*(ctrgtp_xp->x-fiberlagx)*(fiberlagy-ctrgtp->y)/dx/dy + ctrgtp_xp_yp->phi*(fiberlagx-ctrgtp->x)*(fiberlagy-ctrgtp->y)/dx/dy);
		
		if(lagphi>PHImax) lagphi=PHImax;
		if(lagphi<PHImin) lagphi=PHImin;
		
		if(lagphi<=PHIG) lagrho=RHOG;
		else if (lagphi>=PHIL) lagrho=RHOL;
		else lagrho=(RHOL+RHOG)*0.5+(RHOL-RHOG)*0.5*sin((lagphi-(PHIL+PHIG)*0.5)/(PHIL-PHIG)*M_PI);
		
		ctrgt = (int)((fiberlagx+3.0*dx)/dx)*(YCELL_NUM+6)+(int)((fiberlagy+2.5*dy)/dy);	// The point at the bottom left of the moved point
		ctrgt_xp_yp = ctrgt+YCELL_NUM+6+1;
		ctrgt_xp = ctrgt + YCELL_NUM+6;
		ctrgt_yp = ctrgt + 1;
		
		ctrgtp = &CELL[ctrgt];
		ctrgtp_xp = &CELL[ctrgt_xp];
		ctrgtp_yp = &CELL[ctrgt_yp];
		ctrgtp_xp_yp = &CELL[ctrgt_xp_yp];
		_vold_type* voldp = &VOLD[ctrgt];
		_vold_type* voldp_xp = &VOLD[ctrgt_xp];
		_vold_type* voldp_yp = &VOLD[ctrgt_yp];
		_vold_type* voldp_xp_yp = &VOLD[ctrgt_xp_yp];
		
		
		lagv_x = (ctrgtp->v_x*(ctrgtp_xp->x-0.5*dx-fiberlagx)*(ctrgtp_yp->y-fiberlagy)/dx/dy + ctrgtp_xp->v_x*(fiberlagx-ctrgtp->x+0.5*dx)*(ctrgtp_yp->y-fiberlagy)/dx/dy
					+ ctrgtp_yp->v_x*(ctrgtp_xp->x-0.5*dx-fiberlagx)*(fiberlagy-ctrgtp->y)/dx/dy + ctrgtp_xp_yp->v_x*(fiberlagx-ctrgtp->x+0.5*dx)*(fiberlagy-ctrgtp->y)/dx/dy);
		lagoldv_x = (voldp->x*(ctrgtp_xp->x-0.5*dx-fiberlagx)*(ctrgtp_yp->y-fiberlagy)/dx/dy + voldp_xp->x*(fiberlagx-ctrgtp->x+0.5*dx)*(ctrgtp_yp->y-fiberlagy)/dx/dy
					+ voldp_yp->x*(ctrgtp_xp->x-0.5*dx-fiberlagx)*(fiberlagy-ctrgtp->y)/dx/dy + voldp_xp_yp->x*(fiberlagx-ctrgtp->x+0.5*dx)*(fiberlagy-ctrgtp->y)/dx/dy);
		
		
		ctrgt = (int)((fiberlagx+2.5*dx)/dx)*(YCELL_NUM+6)+(int)((fiberlagy+3.0*dy)/dy);	// The point at the bottom left of the moved point
		ctrgt_xp_yp = ctrgt+YCELL_NUM+6+1;
		ctrgt_xp = ctrgt + YCELL_NUM+6;
		ctrgt_yp = ctrgt + 1;
		
		ctrgtp = &CELL[ctrgt];
		ctrgtp_xp = &CELL[ctrgt_xp];
		ctrgtp_yp = &CELL[ctrgt_yp];
		ctrgtp_xp_yp = &CELL[ctrgt_xp_yp];
		voldp = &VOLD[ctrgt];
		voldp_xp = &VOLD[ctrgt_xp];
		voldp_yp = &VOLD[ctrgt_yp];
		voldp_xp_yp = &VOLD[ctrgt_xp_yp];
		
		
		lagv_y = (ctrgtp->v_y*(ctrgtp_xp->x-fiberlagx)*(ctrgtp_yp->y-0.5*dy-fiberlagy)/dx/dy + ctrgtp_xp->v_y*(fiberlagx-ctrgtp->x)*(ctrgtp_yp->y-0.5*dy-fiberlagy)/dx/dy
					+ ctrgtp_yp->v_y*(ctrgtp_xp->x-fiberlagx)*(fiberlagy-ctrgtp->y+0.5*dy)/dx/dy + ctrgtp_xp_yp->v_y*(fiberlagx-ctrgtp->x)*(fiberlagy-ctrgtp->y+0.5*dy)/dx/dy);
		lagoldv_y = (voldp->y*(ctrgtp_xp->x-fiberlagx)*(ctrgtp_yp->y-0.5*dy-fiberlagy)/dx/dy + voldp_xp->y*(fiberlagx-ctrgtp->x)*(ctrgtp_yp->y-0.5*dy-fiberlagy)/dx/dy
					+ voldp_yp->y*(ctrgtp_xp->x-fiberlagx)*(fiberlagy-ctrgtp->y+0.5*dy)/dx/dy + voldp_xp_yp->y*(fiberlagx-ctrgtp->x)*(fiberlagy-ctrgtp->y+0.5*dy)/dx/dy);
		
		
		fiberfx[fibernumber] -= deltar*deltas*lagrho*(lagv_x-lagoldv_x);
		fiberfy[fibernumber] -= deltar*deltas*lagrho*(lagv_y-lagoldv_y);
		
		fibern[fibernumber]+=deltar*deltas*lagrho*(lagv_x-lagoldv_x)*(fiberlagy-fiberp->y);
		fibern[fibernumber]-=deltar*deltas*lagrho*(lagv_y-lagoldv_y)*(fiberlagx-fiberp->x);
	}
}
		
		
		
void discount(void){
/*	char filename_dis[30];
	FILE *disout;
	double phiLG,phiLG2;
	double h_in,a_in;
	double v_average;
	double vx_average,vy_average;
	double a_vx,a_vy;
	double v_average2;
	
	h_in = 0.0;
	a_in = 0.0;
	a_vx = 0.0;
	a_vy = 0.0;
	
	vx_average = 0.0;
	vy_average = 0.0;
	
	if((n%nout==0)||((n+1)%nout==0)){
		phiLG = (PHImax + PHImin) * 0.5;
		phiLG2 = (PHIL + PHIG) * 0.5;
		for(int i=3; i<(XCELL_NUM+3); i++) {
			for(int j=3; j<(YCELL_NUM+3); j++) {
				index0 = j + i*(YCELL_NUM+6);
				cellp= &CELL[index0];
				_flag_type* flagp=&FLAG[index0];
				
				if(flagp->c==1){
					if((cellp->phi<PHIL)&&(cellp->phi>PHIG)){
							h_in += fabs(B*cellp->x-A*cellp->y+D2)/sqrt(A*A+B*B)*(0.5*(PHIL-PHIG)-fabs(cellp->phi-phiLG2))/(PHIL-PHIG)*2.0;
							a_in += (0.5*(PHIL-PHIG)-fabs(cellp->phi-phiLG2))/(PHIL-PHIG)*2.0;
					}
				}
			}
		}
		
		for(int i=8; i<9; i++) {
			for(int j=3; j<(YCELL_NUM+3); j++) {
				index0 = j + i*(YCELL_NUM+6);
				cellp= &CELL[index0];
				_flag_type* flagp=&FLAG[index0];
				
				if(cellp->vof_x<0.99999999999){
					vx_average += cellp->v_x;
					a_vx += 1.0;
				}
				if(cellp->vof_y<0.99999999999) {
					vy_average += cellp->v_y;
					a_vy += 1.0;
				}
			}
		}
		
		h_in /= a_in;
		vx_average /= a_vx;
		vy_average /= a_vy;
		
		if(n==0) h_inold = h_in;
		
		v_average = sqrt(vx_average*vx_average+vy_average*vy_average);
		v_average2 = (h_in - h_inold)/dt;
		
		if(n%nout==0) {
				
			printf("Storing data for Dis\n");
			sprintf(filename_dis,"dis.dat");
			disout = fopen(filename_dis,"a+");
		//	if(n==0)fprintf(disout,"%d %f %f %f \n", n, 0.0, 0.0, 0.0);
		//	else fprintf(disout,"%d %f %f %f \n", n, h_in, h_out, h_in - h_out);
			fprintf(disout,"%d %f %f %f %f %f \n", n, h_in, v_average, v_average2, vx_average, vy_average);
			fclose(disout);
		}
				
		h_inold = h_in;
	}*/
}

		
		
		
		
#ifdef USEMPI
//------------ Sending and receiving boundary surfaces from higher to lower ranks -----------//
void SORmpi1(int my_rank){
	int mailnum = 0;
	int mailsize1;
	MPI_Request req1,req2;
	MPI_Status status;
	int target1,target2;


	target1=my_rank-1;
	target2=my_rank+1;
	mailsize1=YCELL_NUM;
	
	
#ifdef WALLX
	if(my_rank!=0){
		mailnum=0;
		for(int j=3; j<(YCELL_NUM+3); j++){
			int index_m1 = j + (XCELL_NUM/process_num*my_rank+3)*(YCELL_NUM+6);
				
			_cell_type* cellp_m1 = &CELL[index_m1];
				
			mail1[mailnum] = cellp_m1->press;
			mailnum++;
		}
		MPI_Isend(&mail1[0], mailsize1, MPI_DOUBLE,target1,my_rank*10+target1,MPI_COMM_WORLD, &req1);
	}
	if(my_rank!=process_num-1){
		MPI_Irecv(&mail2[0], mailsize1, MPI_DOUBLE,target2,target2*10+my_rank,MPI_COMM_WORLD, &req2);
	}

	if(my_rank!=0) MPI_Wait(&req1, &status);
	if(my_rank!=process_num-1){
		MPI_Wait(&req2, &status);
		mailnum=0;
		for(int j=3; j<(YCELL_NUM+3); j++){
			int index_t2 = j + (XCELL_NUM/process_num*(my_rank+1)+3)*(YCELL_NUM+6);
			
			_cell_type* cellp_t2 = &CELL[index_t2];
			
			cellp_t2->press = mail2[mailnum];
			mailnum++;
		}
	}
#else
	if(target1==-1) target1=process_num-1;
	if(target2==process_num) target2=0;
	
	mailnum=0;
	for(int j=3; j<(YCELL_NUM+3); j++){
		int index_m1 = j + (XCELL_NUM/process_num*my_rank+3)*(YCELL_NUM+6);
			
		_cell_type* cellp_m1 = &CELL[index_m1];
		
		mail1[mailnum] = cellp_m1->press;
		mailnum++;
	}
	MPI_Isend(&mail1[0], mailsize1, MPI_DOUBLE,target1,my_rank*10+target1,MPI_COMM_WORLD, &req1);


	MPI_Irecv(&mail2[0], mailsize1, MPI_DOUBLE,target2,target2*10+my_rank,MPI_COMM_WORLD, &req2);

	MPI_Wait(&req1, &status);

	MPI_Wait(&req2, &status);
	mailnum=0;
	for(int j=3; j<(YCELL_NUM+3); j++){
		int index_t2 = j + (XCELL_NUM/process_num*(my_rank+1)+3)*(YCELL_NUM+6);
		
		if(((XCELL_NUM/process_num*(my_rank+1)+3))==(XCELL_NUM+3)){
			if((j+HEIGHT)<(YCELL_NUM+3)) {
			index_t2 = HEIGHT +j + (XCELL_NUM/process_num*(my_rank+1)+3)*(YCELL_NUM+6);
			}
			else {
				index_t2 = HEIGHT-YCELL_NUM +j + (XCELL_NUM/process_num*(my_rank+1)+3)*(YCELL_NUM+6);
			}
		}
		
		_cell_type* cellp_t2 = &CELL[index_t2];
		
		cellp_t2->press = mail2[mailnum];
		mailnum++;
	}
#endif
}
			//-------------------- End of sending and receiving from higher to lower ranks -----//
#endif


			//-------------- SOR -------------------------//
#ifdef USEMPI
void SOR(int chkflag, int my_rank){
	double dpres;
	double dti;
	double* dxyrho;
	double* poissonS;
	double dpmax;
	double pressaverage;
	double* mailall;
	double* mail0;
	int target1,target2;
	int mailnum = 0;
	MPI_Request req1,req2;
	MPI_Status status;
	int mailsize1, mailsize2;
	double dpmaxall;
	int l;

	
	dxyrho = new double[(XCELL_NUM+6)*(YCELL_NUM+6)];
	poissonS = new double[(XCELL_NUM+6)*(YCELL_NUM+6)];
	mailall = new double[4*(XCELL_NUM+6)*(YCELL_NUM+6)];
	mail0 = new double [XCELL_NUM*YCELL_NUM/process_num];
	
	memset(dxyrho,0,CELL_SIZE*sizeof(double));
	memset(poissonS,0,CELL_SIZE*sizeof(double));
	memset(mailall,0,4*(XCELL_NUM+6)*(YCELL_NUM+6)*sizeof(double));
	memset(mail0,0,XCELL_NUM*YCELL_NUM/process_num*sizeof(double));
	
	pressaverage=0.0;
	
	dti = 1.0/dt;	//1st time
		
	if(my_rank==0){
		for(int i=3; i<(XCELL_NUM+3);i++){
			for(int j=3; j<(YCELL_NUM+3); j++){
				index0   = j     + i*(YCELL_NUM+6);
				index_xp = j     + (i+1)*(YCELL_NUM+6);
				index_yp = (j+1) + i*(YCELL_NUM+6);
				
				cellp    = &CELL[index0];
				cellp_xp = &CELL[index_xp];
				cellp_yp = &CELL[index_yp];
				
				poissonS[index0]=dti*dx*dy*((cellp_xp->v_xhypo - cellp->v_xhypo)*dy+(cellp_yp->v_yhypo - cellp->v_yhypo)*dx);
				
			}
		}

	}
	//----------- MPI_Bcast ------------------//
	if(my_rank==0){

		mailnum=0;
		for(int i=0; i<(XCELL_NUM+6); i++){
			for(int j=0; j<(YCELL_NUM+6); j++){
				index0   = j     + i*(YCELL_NUM+6);
				
				cellp = &CELL[index0];
				
				mailall[mailnum] = cellp->press;
				mailnum++;
			}
		}
		for(int i=0; i<(XCELL_NUM+6); i++){
			for(int j=0; j<(YCELL_NUM+6); j++){
				index0   = j     + i*(YCELL_NUM+6);
					
				cellp = &CELL[index0];
				
				mailall[mailnum] = cellp->rho_x;
				mailnum++;
			}
		}
		for(int i=0; i<(XCELL_NUM+6); i++){
			for(int j=0; j<(YCELL_NUM+6); j++){
				index0   = j     + i*(YCELL_NUM+6);
				
				cellp = &CELL[index0];
				
				mailall[mailnum] = cellp->rho_y;
				mailnum++;
			}
		}
		for(int i=0; i<(XCELL_NUM+6); i++){
			for(int j=0; j<(YCELL_NUM+6); j++){
				index0   = j     + i*(YCELL_NUM+6);
			
				mailall[mailnum] = poissonS[index0];
				mailnum++;
			}
		}
	}



	MPI_Bcast(&mailall[0],4*(XCELL_NUM+6)*(YCELL_NUM+6),MPI_DOUBLE,0,MPI_COMM_WORLD);
	
	
	if(my_rank!=0){
		mailnum=0;
		for(int i=0; i<(XCELL_NUM+6); i++){
			for(int j=0; j<(YCELL_NUM+6); j++){
				index0   = j     + i*(YCELL_NUM+6);
				
				cellp = &CELL[index0];
				
				cellp->press = mailall[mailnum];
				mailnum++;
			}
		}
		for(int i=0; i<(XCELL_NUM+6); i++){
			for(int j=0; j<(YCELL_NUM+6); j++){
				index0   = j     + i*(YCELL_NUM+6);
				
				cellp = &CELL[index0];
				
				cellp->rho_x = mailall[mailnum];
				mailnum++;
			}
		}
		for(int i=0; i<(XCELL_NUM+6); i++){
			for(int j=0; j<(YCELL_NUM+6); j++){
				index0   = j     + i*(YCELL_NUM+6);
				
				cellp = &CELL[index0];
				
				cellp->rho_y = mailall[mailnum];
				mailnum++;
			}
		}
		for(int i=0; i<(XCELL_NUM+6); i++){
			for(int j=0; j<(YCELL_NUM+6); j++){
				index0   = j     + i*(YCELL_NUM+6);
										
				poissonS[index0] = mailall[mailnum];
				mailnum++;
			}
		}
	}
	//-------------- End of MPI_Bcast ----------------//

	mailsize1=YCELL_NUM;
	mailsize2=(XCELL_NUM)*(YCELL_NUM)/process_num;
	target1=my_rank-1;
	target2=my_rank+1;
		
	if(target1==-1) target1=process_num-1;
	if(target2==process_num) target2=0;
	
	for(int i=3;i<(XCELL_NUM+3);i++) {
		for(int j=3;j<(YCELL_NUM+3);j++){
			indexSORdefine(i,j);
			
			dxyrho[index0] = dy*dy*(cellp_xp->rho_x + cellp->rho_x) / cellp_xp->rho_x / cellp->rho_x
								+dx*dx*(cellp_yp->rho_y + cellp->rho_y) / cellp_yp->rho_y / cellp->rho_y;
		}
	}
		
	l=0;
	
	for(;l<MAX;l++) {
		dpmax=0.0;
		
			//------------ Sending and receiving boundary surfaces from lower to higher ranks -----------//
#ifdef WALLX
		if(my_rank!=(process_num-1)){
			mailnum=0;
			for(int j=3; j<(YCELL_NUM+3); j++){
				int index_m2 = j + (XCELL_NUM/process_num*(my_rank+1)+2)*(YCELL_NUM+6);
				
				_cell_type* cellp_m2 = &CELL[index_m2];
				
				mail2[mailnum] = cellp_m2->press;	
				mailnum++;
			}
			MPI_Isend(&mail2[0], mailsize1,MPI_DOUBLE,target2,my_rank*10+target2,MPI_COMM_WORLD,&req2);
		}
		if(my_rank!=0){
			MPI_Irecv(&mail1[0], mailsize1,MPI_DOUBLE,target1,target1*10+my_rank,MPI_COMM_WORLD,&req1);
		}
		
		
		if(my_rank!=(process_num-1)) MPI_Wait(&req2,&status);
		
		if(my_rank!=0) {
			MPI_Wait(&req1,&status);
			mailnum=0;
			for(int j=3; j<(YCELL_NUM+3); j++){
				int index_t1 = j + (XCELL_NUM/process_num*my_rank+2)*(YCELL_NUM+6);
					
				_cell_type* cellp_t1 = &CELL[index_t1];
				
				cellp_t1->press = mail1[mailnum];
				mailnum++;
			}
		}
#else
	mailnum=0;
	for(int j=3; j<(YCELL_NUM+3); j++){
		int index_m2 = j + (XCELL_NUM/process_num*(my_rank+1)+2)*(YCELL_NUM+6);
		
		_cell_type* cellp_m2 = &CELL[index_m2];
		
		mail2[mailnum] = cellp_m2->press;	
		mailnum++;
	}
	MPI_Isend(&mail2[0], mailsize1,MPI_DOUBLE,target2,my_rank*10+target2,MPI_COMM_WORLD,&req2);

	MPI_Irecv(&mail1[0], mailsize1,MPI_DOUBLE,target1,target1*10+my_rank,MPI_COMM_WORLD,&req1);
	
	
	MPI_Wait(&req2,&status);
	
	MPI_Wait(&req1,&status);
	mailnum=0;
		
	for(int j=3; j<(YCELL_NUM+3); j++){
		int index_t1 = j + (XCELL_NUM/process_num*my_rank+2)*(YCELL_NUM+6);
		if(my_rank==0){
			if((j-HEIGHT)>2) {
				index_t1 = -HEIGHT +j + (XCELL_NUM/process_num*my_rank+2)*(YCELL_NUM+6);
			}
			else {
				index_t1 = -HEIGHT+YCELL_NUM +j + (XCELL_NUM/process_num*my_rank+2)*(YCELL_NUM+6);
			}
		}
		
		_cell_type* cellp_t1 = &CELL[index_t1];
		
		cellp_t1->press = mail1[mailnum];
		mailnum++;
	}
		
//	for(int j=3; j<(YCELL_NUM+3); j++){
//		int index_t1 = j + (XCELL_NUM/process_num*my_rank+2)*(YCELL_NUM+6);
//		
//		_cell_type* cellp_t1 = &CELL[index_t1];
//		
//		cellp_t1->press = mail1[mailnum];
//		mailnum++;
//	}
#endif
			//------------ Sending and receiving boundary surfaces from lower to higher ranks -----------//
		
		
		if(!(l%8)){
			//-------------------- From step0 to step1 --------------------//
			for(int i=(XCELL_NUM/process_num*my_rank+3); i<(XCELL_NUM/process_num*(my_rank+1)+3); i+=2){
				for(int j=3; j<(YCELL_NUM+3); j+=2){
					indexSORdefine(i,j);
					
					dpres=(dy*dy*(cellp_xp->press/cellp_xp->rho_x+cellp_xm->press/cellp->rho_x)
							  +dx*dx*(cellp_yp->press/cellp_yp->rho_y+cellp_ym->press/cellp->rho_y)
							  -poissonS[index0])/dxyrho[index0]-cellp->press;
					
					cellp->press+=ALPHA*dpres;
					
					if(dpres<0.0) dpres *= -1.0;
					if(dpres>dpmax) dpmax = dpres;
				}
			}
			for(int i=(XCELL_NUM/process_num*my_rank+3); i<(XCELL_NUM/process_num*(my_rank+1)+3); i+=2){
				for(int j=4; j<(YCELL_NUM+3); j+=2){
					indexSORdefine(i,j);
					
					dpres=(dy*dy*(cellp_xp->press/cellp_xp->rho_x+cellp_xm->press/cellp->rho_x)
							  +dx*dx*(cellp_yp->press/cellp_yp->rho_y+cellp_ym->press/cellp->rho_y)
							  -poissonS[index0])/dxyrho[index0]-cellp->press;
					
					cellp->press+=ALPHA*dpres;
					
					if(dpres<0.0) dpres *= -1.0;
					if(dpres>dpmax) dpmax = dpres;
				}
			}
		
			SORmpi1(my_rank);
		
			//-------------- End of step0 to step1 ---------------//
		
			//-------------- From step2 to step3 ------------------//
			for(int i=(XCELL_NUM/process_num*my_rank+4); i<(XCELL_NUM/process_num*(my_rank+1)+3); i+=2){
				for(int j=3; j<(YCELL_NUM+3); j+=2){
					indexSORdefine(i,j);
					
					dpres=(dy*dy*(cellp_xp->press/cellp_xp->rho_x+cellp_xm->press/cellp->rho_x)
							  +dx*dx*(cellp_yp->press/cellp_yp->rho_y+cellp_ym->press/cellp->rho_y)
							  -poissonS[index0])/dxyrho[index0]-cellp->press;
						
					cellp->press+=ALPHA*dpres;
					
					if(dpres<0.0) dpres *= -1.0;
					if(dpres>dpmax) dpmax = dpres;
				}
			}
			for(int i=(XCELL_NUM/process_num*my_rank+4); i<(XCELL_NUM/process_num*(my_rank+1)+3); i+=2){
				for(int j=4; j<(YCELL_NUM+3); j+=2){
					indexSORdefine(i,j);
					
					dpres=(dy*dy*(cellp_xp->press/cellp_xp->rho_x+cellp_xm->press/cellp->rho_x)
							  +dx*dx*(cellp_yp->press/cellp_yp->rho_y+cellp_ym->press/cellp->rho_y)
							  -poissonS[index0])/dxyrho[index0]-cellp->press;
						
					cellp->press+=ALPHA*dpres;
					
					if(dpres<0.0) dpres *= -1.0;
					if(dpres>dpmax) dpmax = dpres;
				}
			}
		}
		
		if(!(l%8-1)){
			//-------------------- From step0 to step1 --------------------//
			for(int i=(XCELL_NUM/process_num*my_rank+3); i<(XCELL_NUM/process_num*(my_rank+1)+3); i+=2){
				for(int j=(YCELL_NUM+1); j>2; j-=2){
					indexSORdefine(i,j);
					
					dpres=(dy*dy*(cellp_xp->press/cellp_xp->rho_x+cellp_xm->press/cellp->rho_x)
							  +dx*dx*(cellp_yp->press/cellp_yp->rho_y+cellp_ym->press/cellp->rho_y)
							  -poissonS[index0])/dxyrho[index0]-cellp->press;
						
					cellp->press+=ALPHA*dpres;
					
					if(dpres<0.0) dpres *= -1.0;
					if(dpres>dpmax) dpmax = dpres;
				}
			}
			for(int i=(XCELL_NUM/process_num*my_rank+3); i<(XCELL_NUM/process_num*(my_rank+1)+3); i+=2){
				for(int j=(YCELL_NUM+2); j>2; j-=2){
					indexSORdefine(i,j);
					
					dpres=(dy*dy*(cellp_xp->press/cellp_xp->rho_x+cellp_xm->press/cellp->rho_x)
							  +dx*dx*(cellp_yp->press/cellp_yp->rho_y+cellp_ym->press/cellp->rho_y)
							  -poissonS[index0])/dxyrho[index0]-cellp->press;
						
					cellp->press+=ALPHA*dpres;
					
					if(dpres<0.0) dpres *= -1.0;
					if(dpres>dpmax) dpmax = dpres;
				}
			}
		
			SORmpi1(my_rank);
		
			//-------------- End of step0 to step1 ---------------//
		
			//-------------- From step2 to step3 ------------------//
			for(int i=(XCELL_NUM/process_num*my_rank+4); i<(XCELL_NUM/process_num*(my_rank+1)+3); i+=2){
				for(int j=(YCELL_NUM+1); j>2; j-=2){
					indexSORdefine(i,j);
					
					dpres=(dy*dy*(cellp_xp->press/cellp_xp->rho_x+cellp_xm->press/cellp->rho_x)
							  +dx*dx*(cellp_yp->press/cellp_yp->rho_y+cellp_ym->press/cellp->rho_y)
							  -poissonS[index0])/dxyrho[index0]-cellp->press;
						
					cellp->press+=ALPHA*dpres;
					
					if(dpres<0.0) dpres *= -1.0;
					if(dpres>dpmax) dpmax = dpres;
				}
			}
			for(int i=(XCELL_NUM/process_num*my_rank+4); i<(XCELL_NUM/process_num*(my_rank+1)+3); i+=2){
				for(int j=(YCELL_NUM+2); j>2; j-=2){
					indexSORdefine(i,j);
					
					dpres=(dy*dy*(cellp_xp->press/cellp_xp->rho_x+cellp_xm->press/cellp->rho_x)
							  +dx*dx*(cellp_yp->press/cellp_yp->rho_y+cellp_ym->press/cellp->rho_y)
							  -poissonS[index0])/dxyrho[index0]-cellp->press;
						
					cellp->press+=ALPHA*dpres;
					
					if(dpres<0.0) dpres *= -1.0;
					if(dpres>dpmax) dpmax = dpres;
				}
			}
		}
		
		if(!(l%8-2)){
			//-------------------- From step0 to step1 --------------------//
			for(int i=(XCELL_NUM/process_num*(my_rank+1)+1); i>(XCELL_NUM/process_num*(my_rank)+2); i-=2){
				for(int j=(YCELL_NUM+1); j>2; j-=2){
					indexSORdefine(i,j);
					
					dpres=(dy*dy*(cellp_xp->press/cellp_xp->rho_x+cellp_xm->press/cellp->rho_x)
							  +dx*dx*(cellp_yp->press/cellp_yp->rho_y+cellp_ym->press/cellp->rho_y)
							  -poissonS[index0])/dxyrho[index0]-cellp->press;
						
					cellp->press+=ALPHA*dpres;
					
					if(dpres<0.0) dpres *= -1.0;
					if(dpres>dpmax) dpmax = dpres;
				}
			}
			for(int i=(XCELL_NUM/process_num*(my_rank+1)+1); i>(XCELL_NUM/process_num*(my_rank)+2); i-=2){
				for(int j=(YCELL_NUM+2); j>2; j-=2){
					indexSORdefine(i,j);
					
					dpres=(dy*dy*(cellp_xp->press/cellp_xp->rho_x+cellp_xm->press/cellp->rho_x)
							  +dx*dx*(cellp_yp->press/cellp_yp->rho_y+cellp_ym->press/cellp->rho_y)
							  -poissonS[index0])/dxyrho[index0]-cellp->press;
						
					cellp->press+=ALPHA*dpres;
					
					if(dpres<0.0) dpres *= -1.0;
					if(dpres>dpmax) dpmax = dpres;
				}
			}
		
			SORmpi1(my_rank);
		
			//-------------- End of step0 to step1 ---------------//
		
			//-------------- From step2 to step3 ------------------//
			for(int i=(XCELL_NUM/process_num*(my_rank+1)+2); i>(XCELL_NUM/process_num*(my_rank)+2); i-=2){
				for(int j=(YCELL_NUM+1); j>2; j-=2){
					indexSORdefine(i,j);
					
					dpres=(dy*dy*(cellp_xp->press/cellp_xp->rho_x+cellp_xm->press/cellp->rho_x)
							  +dx*dx*(cellp_yp->press/cellp_yp->rho_y+cellp_ym->press/cellp->rho_y)
							  -poissonS[index0])/dxyrho[index0]-cellp->press;
						
					cellp->press+=ALPHA*dpres;
					
					if(dpres<0.0) dpres *= -1.0;
					if(dpres>dpmax) dpmax = dpres;
				}
			}
			for(int i=(XCELL_NUM/process_num*(my_rank+1)+2); i>(XCELL_NUM/process_num*(my_rank)+2); i-=2){
				for(int j=(YCELL_NUM+2); j>2; j-=2){
					indexSORdefine(i,j);
					
					dpres=(dy*dy*(cellp_xp->press/cellp_xp->rho_x+cellp_xm->press/cellp->rho_x)
							  +dx*dx*(cellp_yp->press/cellp_yp->rho_y+cellp_ym->press/cellp->rho_y)
							  -poissonS[index0])/dxyrho[index0]-cellp->press;
						
					cellp->press+=ALPHA*dpres;
					
					if(dpres<0.0) dpres *= -1.0;
					if(dpres>dpmax) dpmax = dpres;
				}
			}
		}
		
		if(!(l%8-3)){
			//-------------------- From step0 to step1 --------------------//
			for(int i=(XCELL_NUM/process_num*(my_rank+1)+1); i>(XCELL_NUM/process_num*(my_rank)+2); i-=2){
				for(int j=3; j<(YCELL_NUM+3); j+=2){
					indexSORdefine(i,j);
					
					dpres=(dy*dy*(cellp_xp->press/cellp_xp->rho_x+cellp_xm->press/cellp->rho_x)
							  +dx*dx*(cellp_yp->press/cellp_yp->rho_y+cellp_ym->press/cellp->rho_y)
							  -poissonS[index0])/dxyrho[index0]-cellp->press;
						
					cellp->press+=ALPHA*dpres;
					
					if(dpres<0.0) dpres *= -1.0;
					if(dpres>dpmax) dpmax = dpres;
				}
			}
			for(int i=(XCELL_NUM/process_num*(my_rank+1)+1); i>(XCELL_NUM/process_num*(my_rank)+2); i-=2){
				for(int j=4; j<(YCELL_NUM+3); j+=2){
					indexSORdefine(i,j);
					
					dpres=(dy*dy*(cellp_xp->press/cellp_xp->rho_x+cellp_xm->press/cellp->rho_x)
							  +dx*dx*(cellp_yp->press/cellp_yp->rho_y+cellp_ym->press/cellp->rho_y)
							  -poissonS[index0])/dxyrho[index0]-cellp->press;
						
					cellp->press+=ALPHA*dpres;
					
					if(dpres<0.0) dpres *= -1.0;
					if(dpres>dpmax) dpmax = dpres;
				}
			}
		
			SORmpi1(my_rank);
		
			//-------------- End of step0 to step1 ---------------//
		
			//-------------- From step2 to step3 ------------------//
			for(int i=(XCELL_NUM/process_num*(my_rank+1)+2); i>(XCELL_NUM/process_num*(my_rank)+2); i-=2){
				for(int j=3; j<(YCELL_NUM+3); j+=2){
					indexSORdefine(i,j);
					
					dpres=(dy*dy*(cellp_xp->press/cellp_xp->rho_x+cellp_xm->press/cellp->rho_x)
							  +dx*dx*(cellp_yp->press/cellp_yp->rho_y+cellp_ym->press/cellp->rho_y)
							  -poissonS[index0])/dxyrho[index0]-cellp->press;
						
					cellp->press+=ALPHA*dpres;
					
					if(dpres<0.0) dpres *= -1.0;
					if(dpres>dpmax) dpmax = dpres;
				}
			}
			for(int i=(XCELL_NUM/process_num*(my_rank+1)+2); i>(XCELL_NUM/process_num*(my_rank)+2); i-=2){
				for(int j=4; j<(YCELL_NUM+3); j+=2){
					indexSORdefine(i,j);
					
					dpres=(dy*dy*(cellp_xp->press/cellp_xp->rho_x+cellp_xm->press/cellp->rho_x)
							  +dx*dx*(cellp_yp->press/cellp_yp->rho_y+cellp_ym->press/cellp->rho_y)
							  -poissonS[index0])/dxyrho[index0]-cellp->press;
						
					cellp->press+=ALPHA*dpres;
					
					if(dpres<0.0) dpres *= -1.0;
					if(dpres>dpmax) dpmax = dpres;
				}
			}
		}
		
		if(!(l%8-4)){
			//-------------------- From step0 to step1 --------------------//
			for(int j=3; j<(YCELL_NUM+3); j+=2){
				for(int i=(XCELL_NUM/process_num*my_rank+3); i<(XCELL_NUM/process_num*(my_rank+1)+3); i+=2){
					indexSORdefine(i,j);
					
					dpres=(dy*dy*(cellp_xp->press/cellp_xp->rho_x+cellp_xm->press/cellp->rho_x)
							  +dx*dx*(cellp_yp->press/cellp_yp->rho_y+cellp_ym->press/cellp->rho_y)
							  -poissonS[index0])/dxyrho[index0]-cellp->press;
						
					cellp->press+=ALPHA*dpres;
					
					if(dpres<0.0) dpres *= -1.0;
					if(dpres>dpmax) dpmax = dpres;
				}
			}
			for(int j=4; j<(YCELL_NUM+3); j+=2){
				for(int i=(XCELL_NUM/process_num*my_rank+3); i<(XCELL_NUM/process_num*(my_rank+1)+3); i+=2){
					indexSORdefine(i,j);
					
					dpres=(dy*dy*(cellp_xp->press/cellp_xp->rho_x+cellp_xm->press/cellp->rho_x)
							  +dx*dx*(cellp_yp->press/cellp_yp->rho_y+cellp_ym->press/cellp->rho_y)
							  -poissonS[index0])/dxyrho[index0]-cellp->press;
						
					cellp->press+=ALPHA*dpres;
					
					if(dpres<0.0) dpres *= -1.0;
					if(dpres>dpmax) dpmax = dpres;
				}
			}
		
			SORmpi1(my_rank);
		
			//-------------- End of step0 to step1 ---------------//
		
			//-------------- From step2 to step3 ------------------//
			for(int j=3; j<(YCELL_NUM+3); j+=2){
				for(int i=(XCELL_NUM/process_num*my_rank+4); i<(XCELL_NUM/process_num*(my_rank+1)+3); i+=2){
					indexSORdefine(i,j);
					
					dpres=(dy*dy*(cellp_xp->press/cellp_xp->rho_x+cellp_xm->press/cellp->rho_x)
							  +dx*dx*(cellp_yp->press/cellp_yp->rho_y+cellp_ym->press/cellp->rho_y)
							  -poissonS[index0])/dxyrho[index0]-cellp->press;
						
					cellp->press+=ALPHA*dpres;
					
					if(dpres<0.0) dpres *= -1.0;
					if(dpres>dpmax) dpmax = dpres;
				}
			}
			for(int j=4; j<(YCELL_NUM+3); j+=2){
				for(int i=(XCELL_NUM/process_num*my_rank+4); i<(XCELL_NUM/process_num*(my_rank+1)+3); i+=2){
					indexSORdefine(i,j);
					
					dpres=(dy*dy*(cellp_xp->press/cellp_xp->rho_x+cellp_xm->press/cellp->rho_x)
							  +dx*dx*(cellp_yp->press/cellp_yp->rho_y+cellp_ym->press/cellp->rho_y)
							  -poissonS[index0])/dxyrho[index0]-cellp->press;
						
					cellp->press+=ALPHA*dpres;
					
					if(dpres<0.0) dpres *= -1.0;
					if(dpres>dpmax) dpmax = dpres;
				}
			}
		}
		
		if(!(l%8-5)){
			//-------------------- From step0 to step1 --------------------//
			for(int j=(YCELL_NUM+1); j>2; j-=2){
				for(int i=(XCELL_NUM/process_num*my_rank+3); i<(XCELL_NUM/process_num*(my_rank+1)+3); i+=2){
					indexSORdefine(i,j);
					
					dpres=(dy*dy*(cellp_xp->press/cellp_xp->rho_x+cellp_xm->press/cellp->rho_x)
							  +dx*dx*(cellp_yp->press/cellp_yp->rho_y+cellp_ym->press/cellp->rho_y)
							  -poissonS[index0])/dxyrho[index0]-cellp->press;
						
					cellp->press+=ALPHA*dpres;
					
					if(dpres<0.0) dpres *= -1.0;
					if(dpres>dpmax) dpmax = dpres;
				}
			}
			for(int j=(YCELL_NUM+2); j>2; j-=2){
				for(int i=(XCELL_NUM/process_num*my_rank+3); i<(XCELL_NUM/process_num*(my_rank+1)+3); i+=2){
					indexSORdefine(i,j);
					
					dpres=(dy*dy*(cellp_xp->press/cellp_xp->rho_x+cellp_xm->press/cellp->rho_x)
							  +dx*dx*(cellp_yp->press/cellp_yp->rho_y+cellp_ym->press/cellp->rho_y)
							  -poissonS[index0])/dxyrho[index0]-cellp->press;
						
					cellp->press+=ALPHA*dpres;
					
					if(dpres<0.0) dpres *= -1.0;
					if(dpres>dpmax) dpmax = dpres;
				}
			}
		
			SORmpi1(my_rank);
		
			//-------------- End of step0 to step1 ---------------//
		
			//-------------- From step2 to step3 ------------------//
			for(int j=(YCELL_NUM+1); j>2; j-=2){
				for(int i=(XCELL_NUM/process_num*my_rank+4); i<(XCELL_NUM/process_num*(my_rank+1)+3); i+=2){
					indexSORdefine(i,j);
					
					dpres=(dy*dy*(cellp_xp->press/cellp_xp->rho_x+cellp_xm->press/cellp->rho_x)
							  +dx*dx*(cellp_yp->press/cellp_yp->rho_y+cellp_ym->press/cellp->rho_y)
							  -poissonS[index0])/dxyrho[index0]-cellp->press;
						
					cellp->press+=ALPHA*dpres;
					
					if(dpres<0.0) dpres *= -1.0;
					if(dpres>dpmax) dpmax = dpres;
				}
			}
			for(int j=(YCELL_NUM+2); j>2; j-=2){
				for(int i=(XCELL_NUM/process_num*my_rank+4); i<(XCELL_NUM/process_num*(my_rank+1)+3); i+=2){
					indexSORdefine(i,j);
					
					dpres=(dy*dy*(cellp_xp->press/cellp_xp->rho_x+cellp_xm->press/cellp->rho_x)
							  +dx*dx*(cellp_yp->press/cellp_yp->rho_y+cellp_ym->press/cellp->rho_y)
							  -poissonS[index0])/dxyrho[index0]-cellp->press;
						
					cellp->press+=ALPHA*dpres;
					
					if(dpres<0.0) dpres *= -1.0;
					if(dpres>dpmax) dpmax = dpres;
				}
			}
		}
		
		if(!(l%8-6)){
			//-------------------- From step0 to step1 --------------------//
			for(int j=(YCELL_NUM+1); j>2; j-=2){
				for(int i=(XCELL_NUM/process_num*(my_rank+1)+1); i>(XCELL_NUM/process_num*(my_rank)+2); i-=2){
					indexSORdefine(i,j);
					
					dpres=(dy*dy*(cellp_xp->press/cellp_xp->rho_x+cellp_xm->press/cellp->rho_x)
							  +dx*dx*(cellp_yp->press/cellp_yp->rho_y+cellp_ym->press/cellp->rho_y)
							  -poissonS[index0])/dxyrho[index0]-cellp->press;
						
					cellp->press+=ALPHA*dpres;
					
					if(dpres<0.0) dpres *= -1.0;
					if(dpres>dpmax) dpmax = dpres;
				}
			}
			for(int j=(YCELL_NUM+2); j>2; j-=2){
				for(int i=(XCELL_NUM/process_num*(my_rank+1)+1); i>(XCELL_NUM/process_num*(my_rank)+2); i-=2){
					indexSORdefine(i,j);
					
					dpres=(dy*dy*(cellp_xp->press/cellp_xp->rho_x+cellp_xm->press/cellp->rho_x)
							  +dx*dx*(cellp_yp->press/cellp_yp->rho_y+cellp_ym->press/cellp->rho_y)
							  -poissonS[index0])/dxyrho[index0]-cellp->press;
						
					cellp->press+=ALPHA*dpres;
					
					if(dpres<0.0) dpres *= -1.0;
					if(dpres>dpmax) dpmax = dpres;
				}
			}
		
			SORmpi1(my_rank);
		
			//-------------- End of step0 to step1 ---------------//
		
			//-------------- From step2 to step3 ------------------//
			for(int j=(YCELL_NUM+1); j>2; j-=2){
				for(int i=(XCELL_NUM/process_num*(my_rank+1)+2); i>(XCELL_NUM/process_num*(my_rank)+2); i-=2){
					indexSORdefine(i,j);
					
					dpres=(dy*dy*(cellp_xp->press/cellp_xp->rho_x+cellp_xm->press/cellp->rho_x)
							  +dx*dx*(cellp_yp->press/cellp_yp->rho_y+cellp_ym->press/cellp->rho_y)
							  -poissonS[index0])/dxyrho[index0]-cellp->press;
						
					cellp->press+=ALPHA*dpres;
					
					if(dpres<0.0) dpres *= -1.0;
					if(dpres>dpmax) dpmax = dpres;
				}
			}
			for(int j=(YCELL_NUM+2); j>2; j-=2){
				for(int i=(XCELL_NUM/process_num*(my_rank+1)+2); i>(XCELL_NUM/process_num*(my_rank)+2); i-=2){
					indexSORdefine(i,j);
					
					dpres=(dy*dy*(cellp_xp->press/cellp_xp->rho_x+cellp_xm->press/cellp->rho_x)
							  +dx*dx*(cellp_yp->press/cellp_yp->rho_y+cellp_ym->press/cellp->rho_y)
							  -poissonS[index0])/dxyrho[index0]-cellp->press;
						
					cellp->press+=ALPHA*dpres;
					
					if(dpres<0.0) dpres *= -1.0;
					if(dpres>dpmax) dpmax = dpres;
				}
			}
		}
		
		if(!(l%8-7)){
			//-------------------- From step0 to step1 --------------------//
			for(int j=3; j<(YCELL_NUM+3); j+=2){
				for(int i=(XCELL_NUM/process_num*(my_rank+1)+1); i>(XCELL_NUM/process_num*(my_rank)+2); i-=2){
					indexSORdefine(i,j);
					
					dpres=(dy*dy*(cellp_xp->press/cellp_xp->rho_x+cellp_xm->press/cellp->rho_x)
							  +dx*dx*(cellp_yp->press/cellp_yp->rho_y+cellp_ym->press/cellp->rho_y)
							  -poissonS[index0])/dxyrho[index0]-cellp->press;
						
					cellp->press+=ALPHA*dpres;
					
					if(dpres<0.0) dpres *= -1.0;
					if(dpres>dpmax) dpmax = dpres;
				}
			}
			for(int j=4; j<(YCELL_NUM+3); j+=2){
				for(int i=(XCELL_NUM/process_num*(my_rank+1)+1); i>(XCELL_NUM/process_num*(my_rank)+2); i-=2){
					indexSORdefine(i,j);
					
					dpres=(dy*dy*(cellp_xp->press/cellp_xp->rho_x+cellp_xm->press/cellp->rho_x)
							  +dx*dx*(cellp_yp->press/cellp_yp->rho_y+cellp_ym->press/cellp->rho_y)
							  -poissonS[index0])/dxyrho[index0]-cellp->press;
						
					cellp->press+=ALPHA*dpres;
					
					if(dpres<0.0) dpres *= -1.0;
					if(dpres>dpmax) dpmax = dpres;
				}
			}
		
			SORmpi1(my_rank);
		
			//-------------- End of step0 to step1 ---------------//
		
		//-------------- From step2 to step3 ------------------//

			for(int j=3; j<(YCELL_NUM+3); j+=2){
				for(int i=(XCELL_NUM/process_num*(my_rank+1)+2); i>(XCELL_NUM/process_num*(my_rank)+2); i-=2){
					indexSORdefine(i,j);
					
					dpres=(dy*dy*(cellp_xp->press/cellp_xp->rho_x+cellp_xm->press/cellp->rho_x)
							  +dx*dx*(cellp_yp->press/cellp_yp->rho_y+cellp_ym->press/cellp->rho_y)
							  -poissonS[index0])/dxyrho[index0]-cellp->press;
						
					cellp->press+=ALPHA*dpres;
					
					if(dpres<0.0) dpres *= -1.0;
					if(dpres>dpmax) dpmax = dpres;
				}
			}
			for(int j=4; j<(YCELL_NUM+3); j+=2){
				for(int i=(XCELL_NUM/process_num*(my_rank+1)+2); i>(XCELL_NUM/process_num*(my_rank)+2); i-=2){
					indexSORdefine(i,j);
					
					dpres=(dy*dy*(cellp_xp->press/cellp_xp->rho_x+cellp_xm->press/cellp->rho_x)
							  +dx*dx*(cellp_yp->press/cellp_yp->rho_y+cellp_ym->press/cellp->rho_y)
							  -poissonS[index0])/dxyrho[index0]-cellp->press;
						
					cellp->press+=ALPHA*dpres;
					
					if(dpres<0.0) dpres *= -1.0;
					if(dpres>dpmax) dpmax = dpres;
				}
			}
		}
		
		
		
		pressbnd();
		
		MPI_Allreduce(&dpmax,&dpmaxall,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
	
//		MPI_Reduce(&dpmax,&dpmaxall,1,MPI_DOUBLE,MPI_MAX,0,MPI_COMM_WORLD);
	
//		MPI_Bcast(&dpmaxall,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
	
//		sprintf(filename_debug,"debug-%d.dat",my_rank);
//		debugout = fopen(filename_debug,"a");
//		fprintf(debugout,"03 OK my_rank= %d \n", my_rank);
//		fclose(debugout);
	
		if(chkflag==1) printf("Poisson Itr : %d  %e\n",l,dpmaxall);
	
		if(dpmaxall<EPS1) break;
	
	}
	
	
	
	
		//-------------- MPI_Gather ----------//
	mailnum=0;
	for(int i=(XCELL_NUM/process_num*my_rank+3); i<(XCELL_NUM/process_num*(my_rank+1))+3; i++){
		for(int j=3; j<(YCELL_NUM+3); j++){
			index0 = j + i*(YCELL_NUM+6);
			
			_cell_type* cellp = &CELL[index0];
			
			mail0[mailnum]=cellp->press;
			mailnum++;
		}
	}
	

	
	
	MPI_Gather(&mail0[0], mailsize2,MPI_DOUBLE,&mailall[0],mailsize2,MPI_DOUBLE,0,MPI_COMM_WORLD);
	
	if(my_rank==0){
		mailnum=0;
		for(int i=3; i<(XCELL_NUM+3); i++){
			for(int j=3; j<(YCELL_NUM+3); j++){
				index0 = j + i*(YCELL_NUM+6);
					
				_cell_type* cellp = &CELL[index0];
				
				cellp->press=mailall[mailnum];
				mailnum++;
			}
		}
	}
	//-------------- End of MPI_Gather ---------------//
	

	for(int i=3; i<(XCELL_NUM+3); i++) {
		for(int j=3; j<(YCELL_NUM+3); j++) {
			index0 = j + i*(YCELL_NUM+6);
				
			_cell_type* cellp = &CELL[index0];
			
			pressaverage += cellp->press;
		}
	}
	
	pressaverage = pressaverage / (double)XCELL_NUM/ (double)YCELL_NUM;
	
/*	for(int i=3; i<(XCELL_NUM+3); i++) {
		for(int j=3; j<(YCELL_NUM+3); j++) {
			index0 = j + i*(YCELL_NUM+6);
			
			_cell_type* cellp = &CELL[index0];
			
			cellp->press = cellp->press - pressaverage;
		}
	}
*/	
	pressbnd();
	
	delete [] dxyrho;
	delete [] poissonS;
	delete [] mailall;
	delete [] mail0;
}
	
#else	// Else for USEMPI	
		
		
void SOR(int chkflag){
	double dpres;
	double dti;
	double* dxyrho;
	double* poissonS;
	double dpmax;
	
	int l;

	
	dxyrho = new double[(XCELL_NUM+6)*(YCELL_NUM+6)];
	poissonS = new double[(XCELL_NUM+6)*(YCELL_NUM+6)];
	
	memset(dxyrho,0,CELL_SIZE*sizeof(double));
	memset(poissonS,0,CELL_SIZE*sizeof(double));
	
	
	dti = 1.0/dt;
		
	
	for(int i=3; i<(XCELL_NUM+3);i++){
		for(int j=3; j<(YCELL_NUM+3); j++){
			index0   = j     + i*(YCELL_NUM+6);
			index_xp = j     + (i+1)*(YCELL_NUM+6);
			index_yp = (j+1) + i*(YCELL_NUM+6);
			
			cellp    = &CELL[index0];
			cellp_xp = &CELL[index_xp];
			cellp_yp = &CELL[index_yp];
			
//			poissonS[index0]=-dx*dy*((cellp_xp->v_xhypo - cellp->v_xhypo)*dy+(cellp_yp->v_yhypo - cellp->v_yhypo)*dx);
			poissonS[index0]=dti*dx*dy*((cellp_xp->v_xhypo - cellp->v_xhypo)*dy+(cellp_yp->v_yhypo - cellp->v_yhypo)*dx);
		}
	}

	
	
	for(int i=3;i<(XCELL_NUM+3);i++) {
		for(int j=3;j<(YCELL_NUM+3);j++){
			indexSORdefine(i,j);
			
			dxyrho[index0] = dy*dy*(cellp_xp->rho_x + cellp->rho_x) / cellp_xp->rho_x / cellp->rho_x
								+dx*dx*(cellp_yp->rho_y + cellp->rho_y) / cellp_yp->rho_y / cellp->rho_y;
		}
	}
		
	l=0;
	
	for(;l<MAX;l++) {
		dpmax=0.0;
		
//		if(!(l%8)){
			//-------------------- From step0 to step1 --------------------//
			for(int i=3; i<(XCELL_NUM+3); i++){
				for(int j=3; j<(YCELL_NUM+3); j++){
					indexSORdefine(i,j);
					_flag_type* flagp=&FLAG[index0];
					_flag_type* flagp_xp=&FLAG[index_xp];
					_flag_type* flagp_yp=&FLAG[index_yp];
					
//					if((flagp->x%4)&&(flagp->y%4)&&(flagp_xp->x%4)&&(flagp_yp->y%4)){
//					if((flagp->c==1)||(flagp->c==6)){
					dpres=(dy*dy*(cellp_xp->press/cellp_xp->rho_x+cellp_xm->press/cellp->rho_x)
							  +dx*dx*(cellp_yp->press/cellp_yp->rho_y+cellp_ym->press/cellp->rho_y)
							  -poissonS[index0])/dxyrho[index0]-cellp->press;
					
					cellp->press+=ALPHA*dpres;
					
//						if((abs(dpres)>8.81e-6)&&(l>15000)) printf("x=%f y=%f \n",cellp->x,cellp->y);
					if(dpres<0.0) dpres *= -1.0;
					if(dpres>dpmax) dpmax = dpres;
//					}
				}
			}
//		}
		
/*		if(!(l%8-1)){
			//-------------------- From step0 to step1 --------------------//
			for(int i=3; i<(XCELL_NUM+3); i++){
				for(int j=(YCELL_NUM+2); j>2; j--){
					indexSORdefine(i,j);
					
					dpres=(dy*dy*(cellp_xp->press/cellp_xp->rho_x+cellp_xm->press/cellp->rho_x)
							  +dx*dx*(cellp_yp->press/cellp_yp->rho_y+cellp_ym->press/cellp->rho_y)
							  -poissonS[index0])/dxyrho[index0]-cellp->press;
						
					cellp->press+=ALPHA*dpres;
					
					if(dpres<0.0) dpres *= -1.0;
					if(dpres>dpmax) dpmax = dpres;
				}
			}
		}
		
		if(!(l%8-2)){
			//-------------------- From step0 to step1 --------------------//
			for(int i=(XCELL_NUM+2); i>2; i--){
				for(int j=(YCELL_NUM+2); j>2; j--){
					indexSORdefine(i,j);
					
					dpres=(dy*dy*(cellp_xp->press/cellp_xp->rho_x+cellp_xm->press/cellp->rho_x)
							  +dx*dx*(cellp_yp->press/cellp_yp->rho_y+cellp_ym->press/cellp->rho_y)
							  -poissonS[index0])/dxyrho[index0]-cellp->press;
						
					cellp->press+=ALPHA*dpres;
					
					if(dpres<0.0) dpres *= -1.0;
					if(dpres>dpmax) dpmax = dpres;
				}
			}
		}
		
		if(!(l%8-3)){
			//-------------------- From step0 to step1 --------------------//
			for(int i=(XCELL_NUM+2); i>2; i--){
				for(int j=3; j<(YCELL_NUM+3); j++){
					indexSORdefine(i,j);
					
					dpres=(dy*dy*(cellp_xp->press/cellp_xp->rho_x+cellp_xm->press/cellp->rho_x)
							  +dx*dx*(cellp_yp->press/cellp_yp->rho_y+cellp_ym->press/cellp->rho_y)
							  -poissonS[index0])/dxyrho[index0]-cellp->press;
						
					cellp->press+=ALPHA*dpres;
					
					if(dpres<0.0) dpres *= -1.0;
					if(dpres>dpmax) dpmax = dpres;
				}
			}
		}
		
		if(!(l%8-4)){
			//-------------------- From step0 to step1 --------------------//
			for(int j=3; j<(YCELL_NUM+3); j++){
				for(int i=3; i<(XCELL_NUM+3); i++){
					indexSORdefine(i,j);
					
					dpres=(dy*dy*(cellp_xp->press/cellp_xp->rho_x+cellp_xm->press/cellp->rho_x)
							  +dx*dx*(cellp_yp->press/cellp_yp->rho_y+cellp_ym->press/cellp->rho_y)
							  -poissonS[index0])/dxyrho[index0]-cellp->press;
						
					cellp->press+=ALPHA*dpres;
					
					if(dpres<0.0) dpres *= -1.0;
					if(dpres>dpmax) dpmax = dpres;
				}
			}
		}
		
		if(!(l%8-5)){
			//-------------------- From step0 to step1 --------------------//
			for(int j=(YCELL_NUM+2); j>2; j--){
				for(int i=3; i<(XCELL_NUM+3); i++){
					indexSORdefine(i,j);
					
					dpres=(dy*dy*(cellp_xp->press/cellp_xp->rho_x+cellp_xm->press/cellp->rho_x)
							  +dx*dx*(cellp_yp->press/cellp_yp->rho_y+cellp_ym->press/cellp->rho_y)
							  -poissonS[index0])/dxyrho[index0]-cellp->press;
						
					cellp->press+=ALPHA*dpres;
					
					if(dpres<0.0) dpres *= -1.0;
					if(dpres>dpmax) dpmax = dpres;
				}
			}
		}
		
		if(!(l%8-6)){
			//-------------------- From step0 to step1 --------------------//
			for(int j=(YCELL_NUM+2); j>2; j--){
				for(int i=(XCELL_NUM+2); i>2; i--){
					indexSORdefine(i,j);
					
					dpres=(dy*dy*(cellp_xp->press/cellp_xp->rho_x+cellp_xm->press/cellp->rho_x)
							  +dx*dx*(cellp_yp->press/cellp_yp->rho_y+cellp_ym->press/cellp->rho_y)
							  -poissonS[index0])/dxyrho[index0]-cellp->press;
						
					cellp->press+=ALPHA*dpres;
					
					if(dpres<0.0) dpres *= -1.0;
					if(dpres>dpmax) dpmax = dpres;
				}
			}
		}
		
		if(!(l%8-7)){
			//-------------------- From step0 to step1 --------------------//
			for(int j=3; j<(YCELL_NUM+3); j++){
				for(int i=(XCELL_NUM+2); i>2; i--){
					indexSORdefine(i,j);
					
					dpres=(dy*dy*(cellp_xp->press/cellp_xp->rho_x+cellp_xm->press/cellp->rho_x)
							  +dx*dx*(cellp_yp->press/cellp_yp->rho_y+cellp_ym->press/cellp->rho_y)
							  -poissonS[index0])/dxyrho[index0]-cellp->press;
						
					cellp->press+=ALPHA*dpres;
					
					if(dpres<0.0) dpres *= -1.0;
					if(dpres>dpmax) dpmax = dpres;
				}
			}
		}
		
		
*/		
		pressbnd();
		
		if(chkflag==1) printf("Poisson Itr : %d  %e\n",l,dpmax);
	
		if(dpmax<EPS1) break;
	
	}
/*	
	for(;l<MAX;l++) {
		dpmax=0.0;
		
		if(!(l%8)){
			//-------------------- From step0 to step1 --------------------//
			for(int i=3; i<(XCELL_NUM+3); i++){
				for(int j=3; j<(YCELL_NUM+3); j++){
					indexSORdefine(i,j);
					
					dpres=(dy*dy*(cellp_xp->press/cellp_xp->rho_x+cellp_xm->press/cellp->rho_x)
							  +dx*dx*(cellp_yp->press/cellp_yp->rho_y+cellp_ym->press/cellp->rho_y)
							  -poissonS[index0])/dxyrho[index0]-cellp->press;
					
					cellp->press+=ALPHA2*dpres;
					
					if(dpres<0.0) dpres *= -1.0;
					if(dpres>dpmax) dpmax = dpres;
				}
			}
		}
		
		if(!(l%8-1)){
			//-------------------- From step0 to step1 --------------------//
			for(int i=3; i<(XCELL_NUM+3); i++){
				for(int j=(YCELL_NUM+2); j>2; j--){
					indexSORdefine(i,j);
					
					dpres=(dy*dy*(cellp_xp->press/cellp_xp->rho_x+cellp_xm->press/cellp->rho_x)
							  +dx*dx*(cellp_yp->press/cellp_yp->rho_y+cellp_ym->press/cellp->rho_y)
							  -poissonS[index0])/dxyrho[index0]-cellp->press;
						
					cellp->press+=ALPHA2*dpres;
					
					if(dpres<0.0) dpres *= -1.0;
					if(dpres>dpmax) dpmax = dpres;
				}
			}
		}
		
		if(!(l%8-2)){
			//-------------------- From step0 to step1 --------------------//
			for(int i=(XCELL_NUM+2); i>2; i--){
				for(int j=(YCELL_NUM+2); j>2; j--){
					indexSORdefine(i,j);
					
					dpres=(dy*dy*(cellp_xp->press/cellp_xp->rho_x+cellp_xm->press/cellp->rho_x)
							  +dx*dx*(cellp_yp->press/cellp_yp->rho_y+cellp_ym->press/cellp->rho_y)
							  -poissonS[index0])/dxyrho[index0]-cellp->press;
						
					cellp->press+=ALPHA2*dpres;
					
					if(dpres<0.0) dpres *= -1.0;
					if(dpres>dpmax) dpmax = dpres;
				}
			}
		}
		
		if(!(l%8-3)){
			//-------------------- From step0 to step1 --------------------//
			for(int i=(XCELL_NUM+2); i>2; i--){
				for(int j=3; j<(YCELL_NUM+3); j++){
					indexSORdefine(i,j);
					
					dpres=(dy*dy*(cellp_xp->press/cellp_xp->rho_x+cellp_xm->press/cellp->rho_x)
							  +dx*dx*(cellp_yp->press/cellp_yp->rho_y+cellp_ym->press/cellp->rho_y)
							  -poissonS[index0])/dxyrho[index0]-cellp->press;
						
					cellp->press+=ALPHA2*dpres;
					
					if(dpres<0.0) dpres *= -1.0;
					if(dpres>dpmax) dpmax = dpres;
				}
			}
		}
		
		if(!(l%8-4)){
			//-------------------- From step0 to step1 --------------------//
			for(int j=3; j<(YCELL_NUM+3); j++){
				for(int i=3; i<(XCELL_NUM+3); i++){
					indexSORdefine(i,j);
					
					dpres=(dy*dy*(cellp_xp->press/cellp_xp->rho_x+cellp_xm->press/cellp->rho_x)
							  +dx*dx*(cellp_yp->press/cellp_yp->rho_y+cellp_ym->press/cellp->rho_y)
							  -poissonS[index0])/dxyrho[index0]-cellp->press;
						
					cellp->press+=ALPHA2*dpres;
					
					if(dpres<0.0) dpres *= -1.0;
					if(dpres>dpmax) dpmax = dpres;
				}
			}
		}
		
		if(!(l%8-5)){
			//-------------------- From step0 to step1 --------------------//
			for(int j=(YCELL_NUM+2); j>2; j--){
				for(int i=3; i<(XCELL_NUM+3); i++){
					indexSORdefine(i,j);
					
					dpres=(dy*dy*(cellp_xp->press/cellp_xp->rho_x+cellp_xm->press/cellp->rho_x)
							  +dx*dx*(cellp_yp->press/cellp_yp->rho_y+cellp_ym->press/cellp->rho_y)
							  -poissonS[index0])/dxyrho[index0]-cellp->press;
						
					cellp->press+=ALPHA2*dpres;
					
					if(dpres<0.0) dpres *= -1.0;
					if(dpres>dpmax) dpmax = dpres;
				}
			}
		}
		
		if(!(l%8-6)){
			//-------------------- From step0 to step1 --------------------//
			for(int j=(YCELL_NUM+2); j>2; j--){
				for(int i=(XCELL_NUM+2); i>2; i--){
					indexSORdefine(i,j);
					
					dpres=(dy*dy*(cellp_xp->press/cellp_xp->rho_x+cellp_xm->press/cellp->rho_x)
							  +dx*dx*(cellp_yp->press/cellp_yp->rho_y+cellp_ym->press/cellp->rho_y)
							  -poissonS[index0])/dxyrho[index0]-cellp->press;
						
					cellp->press+=ALPHA2*dpres;
					
					if(dpres<0.0) dpres *= -1.0;
					if(dpres>dpmax) dpmax = dpres;
				}
			}
		}
		
		if(!(l%8-7)){
			//-------------------- From step0 to step1 --------------------//
			for(int j=3; j<(YCELL_NUM+3); j++){
				for(int i=(XCELL_NUM+2); i>2; i--){
					indexSORdefine(i,j);
					
					dpres=(dy*dy*(cellp_xp->press/cellp_xp->rho_x+cellp_xm->press/cellp->rho_x)
							  +dx*dx*(cellp_yp->press/cellp_yp->rho_y+cellp_ym->press/cellp->rho_y)
							  -poissonS[index0])/dxyrho[index0]-cellp->press;
						
					cellp->press+=ALPHA2*dpres;
					
					if(dpres<0.0) dpres *= -1.0;
					if(dpres>dpmax) dpmax = dpres;
				}
			}
		}
		
		
		
		pressbnd();
		
		if(chkflag==1) printf("Poisson Itr : %d  %e\n",l,dpmax);
	
		if(dpmax<EPS2) break;
	
	}
	
	
	for(;l<MAX;l++) {
		dpmax=0.0;
		
		if(!(l%8)){
			//-------------------- From step0 to step1 --------------------//
			for(int i=3; i<(XCELL_NUM+3); i++){
				for(int j=3; j<(YCELL_NUM+3); j++){
					indexSORdefine(i,j);
					
					dpres=(dy*dy*(cellp_xp->press/cellp_xp->rho_x+cellp_xm->press/cellp->rho_x)
							  +dx*dx*(cellp_yp->press/cellp_yp->rho_y+cellp_ym->press/cellp->rho_y)
							  -poissonS[index0])/dxyrho[index0]-cellp->press;
					
					cellp->press+=ALPHA3*dpres;
					
					if(dpres<0.0) dpres *= -1.0;
					if(dpres>dpmax) dpmax = dpres;
				}
			}
		}
		
		if(!(l%8-1)){
			//-------------------- From step0 to step1 --------------------//
			for(int i=3; i<(XCELL_NUM+3); i++){
				for(int j=(YCELL_NUM+2); j>2; j--){
					indexSORdefine(i,j);
					
					dpres=(dy*dy*(cellp_xp->press/cellp_xp->rho_x+cellp_xm->press/cellp->rho_x)
							  +dx*dx*(cellp_yp->press/cellp_yp->rho_y+cellp_ym->press/cellp->rho_y)
							  -poissonS[index0])/dxyrho[index0]-cellp->press;
						
					cellp->press+=ALPHA3*dpres;
					
					if(dpres<0.0) dpres *= -1.0;
					if(dpres>dpmax) dpmax = dpres;
				}
			}
		}
		
		if(!(l%8-2)){
			//-------------------- From step0 to step1 --------------------//
			for(int i=(XCELL_NUM+2); i>2; i--){
				for(int j=(YCELL_NUM+2); j>2; j--){
					indexSORdefine(i,j);
					
					dpres=(dy*dy*(cellp_xp->press/cellp_xp->rho_x+cellp_xm->press/cellp->rho_x)
							  +dx*dx*(cellp_yp->press/cellp_yp->rho_y+cellp_ym->press/cellp->rho_y)
							  -poissonS[index0])/dxyrho[index0]-cellp->press;
						
					cellp->press+=ALPHA3*dpres;
					
					if(dpres<0.0) dpres *= -1.0;
					if(dpres>dpmax) dpmax = dpres;
				}
			}
		}
		
		if(!(l%8-3)){
			//-------------------- From step0 to step1 --------------------//
			for(int i=(XCELL_NUM+2); i>2; i--){
				for(int j=3; j<(YCELL_NUM+3); j++){
					indexSORdefine(i,j);
					
					dpres=(dy*dy*(cellp_xp->press/cellp_xp->rho_x+cellp_xm->press/cellp->rho_x)
							  +dx*dx*(cellp_yp->press/cellp_yp->rho_y+cellp_ym->press/cellp->rho_y)
							  -poissonS[index0])/dxyrho[index0]-cellp->press;
						
					cellp->press+=ALPHA3*dpres;
					
					if(dpres<0.0) dpres *= -1.0;
					if(dpres>dpmax) dpmax = dpres;
				}
			}
		}
		
		if(!(l%8-4)){
			//-------------------- From step0 to step1 --------------------//
			for(int j=3; j<(YCELL_NUM+3); j++){
				for(int i=3; i<(XCELL_NUM+3); i++){
					indexSORdefine(i,j);
					
					dpres=(dy*dy*(cellp_xp->press/cellp_xp->rho_x+cellp_xm->press/cellp->rho_x)
							  +dx*dx*(cellp_yp->press/cellp_yp->rho_y+cellp_ym->press/cellp->rho_y)
							  -poissonS[index0])/dxyrho[index0]-cellp->press;
						
					cellp->press+=ALPHA3*dpres;
					
					if(dpres<0.0) dpres *= -1.0;
					if(dpres>dpmax) dpmax = dpres;
				}
			}
		}
		
		if(!(l%8-5)){
			//-------------------- From step0 to step1 --------------------//
			for(int j=(YCELL_NUM+2); j>2; j--){
				for(int i=3; i<(XCELL_NUM+3); i++){
					indexSORdefine(i,j);
					
					dpres=(dy*dy*(cellp_xp->press/cellp_xp->rho_x+cellp_xm->press/cellp->rho_x)
							  +dx*dx*(cellp_yp->press/cellp_yp->rho_y+cellp_ym->press/cellp->rho_y)
							  -poissonS[index0])/dxyrho[index0]-cellp->press;
						
					cellp->press+=ALPHA3*dpres;
					
					if(dpres<0.0) dpres *= -1.0;
					if(dpres>dpmax) dpmax = dpres;
				}
			}
		}
		
		if(!(l%8-6)){
			//-------------------- From step0 to step1 --------------------//
			for(int j=(YCELL_NUM+2); j>2; j--){
				for(int i=(XCELL_NUM+2); i>2; i--){
					indexSORdefine(i,j);
					
					dpres=(dy*dy*(cellp_xp->press/cellp_xp->rho_x+cellp_xm->press/cellp->rho_x)
							  +dx*dx*(cellp_yp->press/cellp_yp->rho_y+cellp_ym->press/cellp->rho_y)
							  -poissonS[index0])/dxyrho[index0]-cellp->press;
						
					cellp->press+=ALPHA3*dpres;
					
					if(dpres<0.0) dpres *= -1.0;
					if(dpres>dpmax) dpmax = dpres;
				}
			}
		}
		
		if(!(l%8-7)){
			//-------------------- From step0 to step1 --------------------//
			for(int j=3; j<(YCELL_NUM+3); j++){
				for(int i=(XCELL_NUM+2); i>2; i--){
					indexSORdefine(i,j);
					
					dpres=(dy*dy*(cellp_xp->press/cellp_xp->rho_x+cellp_xm->press/cellp->rho_x)
							  +dx*dx*(cellp_yp->press/cellp_yp->rho_y+cellp_ym->press/cellp->rho_y)
							  -poissonS[index0])/dxyrho[index0]-cellp->press;
						
					cellp->press+=ALPHA3*dpres;
					
					if(dpres<0.0) dpres *= -1.0;
					if(dpres>dpmax) dpmax = dpres;
				}
			}
		}
		
		
		
		pressbnd();
		
		if(chkflag==1) printf("Poisson Itr : %d  %e\n",l,dpmax);
	
		if(dpmax<EPS3) break;
	
	}
*/	
/*	for(int i=3; i<(XCELL_NUM+3); i++) {
		for(int j=3; j<(YCELL_NUM+3); j++) {
			index0 = j + i*(YCELL_NUM+6);
				
			_cell_type* cellp = &CELL[index0];
			
			pressaverage += cellp->press;
		}
	}
	
	pressaverage = pressaverage / (double)XCELL_NUM/ (double)YCELL_NUM;
	
	for(int i=3; i<(XCELL_NUM+3); i++) {
		for(int j=3; j<(YCELL_NUM+3); j++) {
			index0 = j + i*(YCELL_NUM+6);
			
			_cell_type* cellp = &CELL[index0];
			
			cellp->press = cellp->press - pressaverage;
		}
	}
*/	
	pressbnd();
	
	delete [] dxyrho;
	delete [] poissonS;
}
	
#endif	//USEMPI
