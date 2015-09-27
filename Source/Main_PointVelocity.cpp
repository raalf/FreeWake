#include "general.h"
#include "PerfCode.h"
main()
{
//computes velocity due to time step file generated from Free Wake 2007 at
//point P
//Step 1 read in information about point and time step from matlab file
//Step 2 read in time step file for the specified time step
//Step 3 compute velocities at point P
//Step 4 create output file of velocites at point P
//program checks computation of the induced velocity in P due to a DVE
//
//these program allows to compute the velocity field that is induced by a
//given wing and the corresponding wake
//
//		!!!ATTENTION!!!
//
// The output path is OUTPUT_PATH, which is defined in general.h!!
//

	double P[3],w_ind[3],tempA[3];
	int i,j,k;				// loop counter
	int plane,type;
	int imax,jmax,kmax;
	int wing,span,time;		//more loop counters
	double xmin,xmax,ymin,ymax,zmin,zmax;
	double xstep,ystep,zstep;
	double tempS;
	char ch;		//generic character


//	GENERAL info;			//general info
	DVE *surfaceDVE, **wakeDVE;

	int timestep;
	char iofile[125];	//input-output-file
	FILE *fp,*fs;

//===================================================================//
//START  Step 1 
//read in information about point and time step from matlab file
//===================================================================//

	//creates file name for file with point information
	sprintf(iofile,"%s%s",OUTPUT_PATH,"pointinfo.txt");

	// checks if input file exists
	if ((fp = fopen(iofile, "r"))== NULL)
	{
		printf("File could not be opened, stupid:\n");
		exit(1);
	}

	fscanf(fp,"%d %lf %lf %lf", &timestep,&P[0],&P[1],&P[2]);

	//opens input file
	fp = fopen(iofile, "r");

	fclose(fp);


printf("%d   %lf   %lf  %lf\n", timestep,P[0],P[1],P[2]);

//===================================================================//
//END Step 1 
//===================================================================//


//===================================================================//
//	START Step 2	
//	set-up to read in timestep-file
//===================================================================//


	//creates file name timestep##.txt ## is number of timestep
	sprintf(iofile,"%s%s%d%s",OUTPUT_PATH,"timestep",timestep,".txt");

	// checks if input file exists
	if ((fp = fopen(iofile, "r"))== NULL)
	{
		printf("File could not be opened, stupid:\n");
		exit(1);
	}

	//opens input file
	fp = fopen(iofile, "r");

//===================================================================//
//			reads in general info
//===================================================================//
	//find the ':'-sign in input file before program version
	do	ch = fgetc(fp);
	while (ch!=':');

	//find the ':'-sign in input file before reference area
	do	ch = fgetc(fp);
	while (ch!=':');
	fscanf(fp,"%lf", &info.S);
//printf(" %lf\n",info.S);

	//find the ':'-sign in input file before aspect ratio
	do	ch = fgetc(fp);
	while (ch!=':');
	fscanf(fp,"%lf", &info.AR);
//printf(" %lf",info.AR);

	//find the ':'-sign in input file before alpha
	do	ch = fgetc(fp);
	while (ch!=':');
	fscanf(fp,"%lf", &info.alpha);
//	info.alpha *= DtR;

	//find the ':'-sign in input file before beta
	do	ch = fgetc(fp);
	while (ch!=':');
	fscanf(fp,"%lf", &info.beta);
	info.beta *= DtR;

	//find the ':'-sign in input file before symmetry condition
	do	ch = fgetc(fp);
	while (ch!=':');
	fscanf(fp,"%d", &info.sym);

	//find the ':'-sign in input file before timestep
	do	ch = fgetc(fp);
	while (ch!=':');
//	//reads timestep
//	fscanf(fp,"%*d", &timestep);

	//find the ':'-sign in input file before number of elements along span
	do	ch = fgetc(fp);
	while (ch!=':');
	//reads nospanelements
	fscanf(fp,"%d", &info.nospanelement);
//printf("%d\n", info.nospanelement);

	//find the ':'-sign in input file before number of elements along chord
	do	ch = fgetc(fp);
	while (ch!=':');
	//reads nospanelements
	fscanf(fp,"%d", &info.m);
//printf("%d\n", info.m);

	//number of surface DVEs
	info.noelement = info.m*info.nospanelement;

	//find the ':'-sign in input file before number of wings
	do	ch = fgetc(fp);
	while (ch!=':');
	//reads nospanelements
	fscanf(fp,"%d\n", &info.nowing);

	//find the ':'-sign in input file before number of span indices
	do	ch = fgetc(fp);
	while (ch!=':');

	for(wing=0;wing<info.nowing;wing++)
	{
		fscanf(fp,"%d", &info.wing1[wing]);
		fscanf(fp,"%d", &info.wing2[wing]);
	}

	fclose(fp);

//===================================================================//
//			allocates memory
//===================================================================//

	ALLOC1D(&surfaceDVE,info.noelement);
	ALLOC2D(&wakeDVE,timestep+1,info.nospanelement);

//===================================================================//
//		reads info surface info
//===================================================================//

	//read in information on surface and wake DVEs of timestepfile
	Read_Timestep(timestep,surfaceDVE,wakeDVE);
										//Subroutine in read_input.cpp

//===================================================================//
					//-- Reading in is DONE --
		//Now the singularity factors are being computed
//===================================================================//

	//	computes singularity factors for wake
	//loop over timesteps
	for(time=0;time<=timestep;time++)
	//loop over number of wings
	for(wing=0;wing<info.nowing;wing++)
	{
		//is the wing symmetrical or not?
		if(info.sym == 1)	//decay factor is 1% of tip-element half-span
			tempS = 0.01*wakeDVE[time][info.wing2[wing]].eta;
		else//wing has two tips, possibly different in geometry
		{	//in that case, decay factor is 1% of the shorter half-span
			if(  wakeDVE[time][info.wing1[wing]].eta
			   < wakeDVE[time][info.wing2[wing]].eta)
						tempS = 0.01*wakeDVE[time][info.wing1[wing]].eta;
			else 		tempS = 0.01*wakeDVE[time][info.wing2[wing]].eta;
		}

	//loop over wale DVEs of current timestep
		for(span=info.wing1[wing];span<=info.wing2[wing];span++)
			wakeDVE[time][span].singfct = tempS;  //assigning decay factor
	}//next wing

	//singularity factor of lifting surface DVEs
	//computing decaying factor for added singularity at wing tip
	k=0;	//initializing index counter of first DVE of wing
	for(wing=0;wing<info.nowing;wing++)
	{
		//index of last DVE of this wing (located at tip and trail. edge)
		span = k + (info.wing2[wing]-info.wing1[wing]+1)*info.m - 1;

		//is the wing symmetrical or not?
		if(info.sym == 1)	//decay factor is 1% of tip-element half-span
			tempS = 0.01*surfaceDVE[span].eta;
		else//wing has two tips, possibly different in geometry
		{	//in that case, decay factor is 1% of the shorter half-span
			if(  surfaceDVE[k].eta < surfaceDVE[span].eta)
						tempS = 0.01*surfaceDVE[k].eta;
			else 		tempS = 0.01*surfaceDVE[span].eta;
		}

		//loop over surface DVEs of current wing
		for(i=k;i<=span;i++)
			surfaceDVE[i].singfct = tempS;  //assigning decay factor

		k = span+1;	//updating index for next wing
	}//next wing

//===================================================================//
		//-- Computation of singularity factors is DONE --
		//Now comes the part where the velocity is computed
//===================================================================//
//===================================================================//
//END Step 2
//===================================================================//



//===================================================================//
//START Step 3		
//computing velocity in P
//===================================================================//

	//velocity is induced by all surface and wake DVE's in point P
	DVE_Induced_Velocity(info,P,surfaceDVE,wakeDVE,timestep,w_ind);
				 			//subroutine in induced_velocity.cpp


//computing induced velocity of most right point of each wing.
//right wingtip points are stored in xright and the velocity in uright
//for(wing=0;wing<info.nowing;wing++)
//{
//	for(time=0;time<=timestep;time++)
//	{
		//span index of DVE that is the most right on of current wing
//		span = info.wing2[wing];
//				if(time==0) 		type = 3;
//		else 	if(time==timestep)	type = 2;
//		else 						type = 1;
//		DVE_Tip_Induced_Velocity(info,wakeDVE[time][span],P,tempA,type);
		 						//subroutine in induced_velocity.cpp

//		w_ind[0] -= tempA[0];
//		w_ind[1] -= tempA[1];
//		w_ind[2] -= tempA[2];
//	}
//}
//===================================================================//
//END Step 3
//===================================================================//
//	START Step 4
//	saving output (coordinates and induced velocity components)
//===================================================================//
printf("\n Induced Velocities \n");
printf("%lf  %lf  %lf\n",w_ind[0],w_ind[1],w_ind[2]);


// Create Output file
sprintf(iofile,"%s%s",OUTPUT_PATH,"velocityinfo.txt");

	// checks if input file exists
	if ((fp = fopen(iofile, "r"))== NULL)
	{
		printf("File could not be opened, stupid:\n");
		exit(1);
	}
// Open output file and write w_ind	
 fs = fopen(iofile,"w");
 fprintf(fs,"%1f  %1f  %1f", w_ind[0],w_ind[1],w_ind[2]);
 fclose(fs);
// End Create Output file


printf("\n DONE\n");
//scanf("%d",&timestep);

return(0);
}
//===================================================================//
		//END of program
//===================================================================//
