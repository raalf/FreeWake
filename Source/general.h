#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
//#include <iostream.h>
#include <math.h>
#include <sys/stat.h>

#if defined(_WIN32)
#include "direct.h" //required only on windows
//#else 
#endif


#include "typedef.h"
#include "alloc.h"
#include "vector_algebra.h"
#include "ref_frame_transform.h"

#define Pi  3.1415926535897931
#define DtR Pi/180
#define RtD 180/Pi
#define DBL_EPS 1e-14                   //"zero" definition
#define OUTPUT_PATH "output/"           //directory where output is saved
#define PROGRAM_VERSION "FreeWake2020_dev V0.4"
#define AIRFOIL_PATH "airfoils/"        //directory with airfoils
#define CAMBER_PATH "inputs/camber/"  //directory with camber data
#define flagSTARFORCE 1   //flag = 1 skip computing aero laods until end


//global Variables
FILE *test;				//file for test output during debugging


GENERAL info;			//general input information

PANEL *panelPtr;		//pointer holds information on panel geometry

BOUND_VORTEX *elementPtr;//pointer holds information on elementary wings

DVE *surfacePtr;		//pointer to surface Distributed-Vorticity elements
DVE **wakePtr;			//pointer to wake DVE

double Nt_free[3], Nt_ind[3];
						//magnitude of induced and free stream normal
						//forces/density of total wing

double CF[3];           //total forces in wind axis system (CFX,CFY,CFY)
double **Cf;            //section forces in wind axis system (Cfx,Cfy,Cfz)
double *CDi_DVE;		//total ind. drag (Eppler) with DVEs for each timestep
double CDi_finit;		//total induced drag with DVEs after all timesteps
double Cl, Cm, Cn;		//
//8-8-07 G.B.
double alpha1,alpha2,alphastep;	//AOA loop, beginning, end, stepsie

