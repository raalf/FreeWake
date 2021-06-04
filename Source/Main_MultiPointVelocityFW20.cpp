#define _CRT_SECURE_NO_WARNINGS
#include "general.h"
#include "FreeWakeWing.h"
int main()
{
//computes velocity due to time step file generated from Free Wake 2020 at
//multiple points 
//Step 1 read in information about points surveyed
//Step 2 read TDVE file
//Step 3 compute velocities at point P
//Step 4 create output file of velocites at point P
//program checks computation of the induced velocity in P due to a DVE
//
//these program allows to compute the velocity field that is induced by a
//given wing and the corresponding wake
//
//		!!!ATTENTION!!!
//
// The output path is OUTPUT_PATH (info.output), which is defined in general.h!!
//
	struct Points
	{double x;
	double y;
	double z;
	}Points[10000];

	struct Vels
	{double u;
	double v;
	double w;
	}Vels[10000];

	double P[3],w_ind[3],tempA[3];
	int i,j,k;				// loop counter
	int plane,type;
	int imax,jmax,kmax;
	int wing,span,time,n;		//more loop counters
	int numpoints=0;			//number of points
	int timestep,nospan,nosurface;
	double xmin,xmax,ymin,ymax,zmin,zmax;
	double xstep,ystep,zstep;
	double tempS;
	char ch;		//generic character
	

	GENERAL info;			//general info
	DVE *surfaceDVE, **wakeDVE;

	char iofile[125];	//input-output-file
	FILE *fp,*fs;

//===================================================================//
//START  Step 1 
//read in information about point and time step from matlab file
//===================================================================//

	//creates file name for file with point information
	sprintf(iofile,"%s%s",info.output,"pointinfo.txt");
	printf("output file is %s \n",iofile);

	// checks if input file exists
	if ((fp = fopen(iofile, "r"))== NULL)
	{
		printf("File could not be opened, stupid:\n");
		exit(1);
	}
	//opens input file
	fp = fopen(iofile, "r");
	//scan number of points
	fscanf(fp,"%d ", &numpoints);
	
	for (i=0;i<numpoints;i++)
		fscanf(fp,"%d %lf %lf %lf", &n,&Points[i].x,&Points[i].y,&Points[i].z);

	fclose(fp);

//===================================================================//
//END Step 1 
//===================================================================//

//===================================================================//
//START  Step 2
//read in information about point and time step from TDVE file
//===================================================================//

//	1. user input of timestep and intervalls that is to be plotted
	printf("\nThis is the Velocity Field program for %s\n\n",PROGRAM_VERSION);
	printf("The the file with the points needs to be located in the ");
	printf("current directory.\n");
	printf("Please enter full name with paht of _TDVE file ");
	printf("that is to be used: \n");
	scanf("%s",iofile);

//	printf("\nhardwired as %s\n","output/FW_Input/FW_Input_TDVE#1.txt");
//	sprintf(iofile,"%s","output/FW_Input/FW_Input_TDVE#1.txt");1

//	2. reading in data of timestep

	//creates file name timestep##.txt ## is number of timestep
//	sprintf(iofile,"%s%s%d%s",OUTPUT_PATH,"timestep",timestep,".txt");
//	sprintf(iofile,"%s%d%s","timestep",timestep,".txt");

	// checks if input file exists
	if ((fp = fopen(iofile, "r"))== NULL)
	{
		printf("File could not be opened,:\n");
		exit(1);
	}

	//opens input file
	fp = fopen(iofile, "r");

	//find the ':'-sign in input file before program version
	do	ch = fgetc(fp);
	while (ch!=':');

	//find the ':'-sign in input file before input file
	do	ch = fgetc(fp);
	while (ch!=':');

	//find the ':'-sign in input file before flight condition
	do	ch = fgetc(fp);
	while (ch!=':');

	//find the ':'-sign in input file before angle of attack
	do	ch = fgetc(fp);
	while (ch!=':');

	//find the ':'-sign in input file before ref. area
	do	ch = fgetc(fp);
	while (ch!=':');

	//find the ':'-sign in input file before aspect ratio
	do	ch = fgetc(fp);
	while (ch!=':');

	//find the ':'-sign in input file before timestep
	do	ch = fgetc(fp);
	while (ch!=':');
	//reads timestep
	fscanf(fp,"%d", &timestep);

	//find the ':'-sign in input file before number of span-elements
	do	ch = fgetc(fp);
	while (ch!=':');
	//reads number of elements in spanwise direction
	fscanf(fp,"%d", &nospan);

	//find the ':'-sign in input file before number of surface elements
	do	ch = fgetc(fp);
	while (ch!=':');
	fscanf(fp,"%d", &nosurface);

	//find the first '#'-sign at end of header of surface elements
	do	ch = fgetc(fp);
	while (ch!='#');

	//allocates memory for surfaceDVE
	ALLOC1D(&surfaceDVE,nosurface);

	for(n=0;n<nosurface;n++)
	{
		//read ref. point, nu, epsilon, sweeps, local U,span, chord
		fscanf(fp,"%*d %lf %lf %lf \
					%*lf %*lf %*lf %*lf \
				   %lf %lf %lf \
				   %lf %lf \
				   %lf %lf %lf \
				   %lf %lf %lf",\
		   &surfaceDVE[n].xo[0],&surfaceDVE[n].xo[1],&surfaceDVE[n].xo[2],\
		   &surfaceDVE[n].A,&surfaceDVE[n].B,&surfaceDVE[n].C,\
		   &surfaceDVE[n].eta,&surfaceDVE[n].xsi,\
		   &surfaceDVE[n].nu,&surfaceDVE[n].epsilon,&surfaceDVE[n].psi,\
		   &surfaceDVE[n].phiLE,&surfaceDVE[n].phi0,&surfaceDVE[n].phiTE);

			//convert angles to radians
		surfaceDVE[n].nu 		*= DtR;
		surfaceDVE[n].epsilon 	*= DtR;
		surfaceDVE[n].psi		*= DtR;
		surfaceDVE[n].phiLE		*= DtR;
		surfaceDVE[n].phi0		*= DtR;
		surfaceDVE[n].phiTE		*= DtR;

		//find the beginning of next span information
		do	ch = fgetc(fp);
		while (ch!='\n');
	}

	printf(" Done reading in SDVE\n");
	//find the second '#'-sign at end of header of wake elements
	do ch = fgetc(fp);
	while (ch!='#');

    //allocates memory for wakeDVE
	ALLOC2D(&wakeDVE,timestep+1,nospan);

	for(time=0;time<=timestep;time++)
	{
		for(span=0;span<nospan;span++)
		{
  			//read ref. point, nu, epsilon, sweeps, local U,span, chord
			fscanf(fp,"%*d %*d \
				%lf %lf %lf \
				%lf %lf %lf \
				%lf %lf %lf \
				%lf %lf %lf \
				%lf %lf %lf \
				%lf %lf %lf",\
				   &wakeDVE[time][span].xo[0],&wakeDVE[time][span].xo[1],&wakeDVE[time][span].xo[2],\
				   &wakeDVE[time][span].nu,&wakeDVE[time][span].epsilon,&wakeDVE[time][span].psi,\
				   &wakeDVE[time][span].u[0],&wakeDVE[time][span].u[1],&wakeDVE[time][span].u[2],\
				   &wakeDVE[time][span].A,&wakeDVE[time][span].B,&wakeDVE[time][span].C,\
				   &wakeDVE[time][span].eta,&wakeDVE[time][span].xsi,&wakeDVE[time][span].singfct,\
				   &wakeDVE[time][span].phiLE,&wakeDVE[time][span].phi0,&wakeDVE[time][span].phiTE);


			//convert angles to radians
			wakeDVE[time][span].nu 		*= DtR;
			wakeDVE[time][span].epsilon *= DtR;
			wakeDVE[time][span].psi		*= DtR;
			wakeDVE[time][span].phiLE	*= DtR;
			wakeDVE[time][span].phi0	*= DtR;
			wakeDVE[time][span].phiTE	*= DtR;

			//find the beginning of next span information
			do	ch = fgetc(fp);
			while (ch!='\n');
		}
		//find the beginning of next time index
		do	ch = fgetc(fp);
		while (ch!='\n');
	}

	//closes input file
	fclose(fp);
	printf(" Done reading in WDVE\n");
//printf("%d   %lf   %lf  %lf\n", timestep,P[0],P[1],P[2]);

//===================================================================//
					//-- Reading in is DONE --
		//Now the singularity factors are being computed
//===================================================================//

	//singularity factor of lifting surface DVEs
	//computing decaying factor for added singularity at wing tip
	k=0;	//initializing index counter of first DVE of wing
	tempS = surfaceDVE[0].eta;

	//finding shortest halfspan
	for(n=0;n<nosurface;n++)
	{
		if(surfaceDVE[n].eta<tempS)
			tempS = surfaceDVE[n].eta;
	}
	tempS = 0.01*tempS; //singularity factor for surface DVEs

	for(n=0;n<nosurface;n++)
		surfaceDVE[n].singfct = tempS;  //assigning singularity factor
//===================================================================//
		//-- Computation of singularity factors is DONE --
		//Now comes the part where the velocity is computed
//===================================================================//
//===================================================================//
//END Step 2
//===================================================================//
	printf(" Done reading with updating inputs\n");

//===================================================================//
//START Step 3		
//computing velocity in P
//===================================================================//
	info.noelement = nosurface;
	info.nospanelement = nospan;

	for(i=0;i<numpoints;i++)
	{
		P[0] = Points[i].x;
		P[1] = Points[i].y;
		P[2] = Points[i].z;
		//velocity is induced by all surface and wake DVE's in point P
		DVE_Induced_Velocity(info,P,surfaceDVE,wakeDVE,timestep,w_ind);
				 			//subroutine in induced_velocity.cpp

//		printf("%1f  %1f  %1f ", Points[i].x,Points[i].y,Points[i].z);
//		printf("%1f  %1f  %1f\n",w_ind[0],w_ind[1],w_ind[2]);

		Vels[i].u = w_ind[0];
		Vels[i].v = w_ind[1];
		Vels[i].w = w_ind[2];
	}
//===================================================================//
//END Step 3
//===================================================================//
//	START Step 4
//	saving output (coordinates and induced velocity components)
//===================================================================//
//printf("\n Induced Velocities \n");
//printf("%lf  %lf  %lf\n",w_ind[0],w_ind[1],w_ind[2]);


// Create Output file
sprintf(iofile,"%s","velocityinfo.txt");

	// checks if input file exists
	// if ((fp = fopen(iofile, "r"))== NULL)
	//{
	//	printf("File could not be opened, stupid:\n");
	//	exit(1);
	//}
// Open output file and write w_ind	
 fs = fopen(iofile,"w");

 fprintf(fs,"%d\n",numpoints);

for(i=0;i<numpoints;i++)
{
	fprintf(fs,"%1f  %1f  %1f ", Points[i].x,Points[i].y,Points[i].z);
	fprintf(fs,"%1f  %1f  %1f\n", Vels[i].u,Vels[i].v,Vels[i].w);
 }

 fclose(fs);
// End Create Output file

// printf("\n DONE\n");
printf("\nIt's done.  Enter anything, anything: ");
//scanf("%d",&timestep);

return(0);
}
//===================================================================//
		//END of program
//===================================================================//
