#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
//#include <iostream.h>
#include <math.h>

#include "typedef.h"
#include "alloc.h"
#include "vector_algebra.h"
#include "ref_frame_transform.h"

#define Pi  3.1415926535897931
#define DtR Pi/180
#define RtD 180/Pi
#define DBL_EPS 1e-14
#define OUTPUT_PATH "output/"
#define PROGRAM_VERSION "FreeWake2018_Omega"
#define AIRFOIL_PATH "airfoils/"


//global Variables
FILE *test;				//file for test output during debugging


GENERAL info;			//general input information

PANEL *panelPtr;		//pointer holds information on panel geometry

BOUND_VORTEX *elementPtr;//pointer holds information on elementary wings

DVE *surfacePtr;		//pointer to surface Distributed-Vorticity elements
DVE **wakePtr;			//pointer to wake DVE

double Nt_free[2], Nt_ind[2];
						//magnitude of induced and free stream normal
						//forces/density of total wing

double **CN;			//total normal forces, CL,CLi,CY,CYi for each timestep

double *CDi_DVE;		//total ind. drag (Eppler) with DVEs for each timestep
double CDi_finit;		//total induced drag with DVEs after all timesteps

//8-8-07 G.B.
double alpha1,alpha2,alphastep;	//AOA loop, beginning, end, stepsie

