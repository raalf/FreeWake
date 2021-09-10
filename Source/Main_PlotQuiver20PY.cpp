#include "general.h"
#include "PerfCode.h"
int main()
{
//This program creates and runs a python script to plot the surface and wake DVEs 
//	with a vectors defined in the output/quiver.txt file.
// The quiver file must be columns of data as: x,y,z,u,v,w
//		where x,y,z is the point where vector begins and u,v,w are the components
//		of the vector.
//
// ***It is suggested to use the function CreateQuiverFile found in write_output.cpp
// 		to create the quiver file
//
//This requires matplotlib (version 1.5 or newer) and python to be installed
//
//This program is entirely based on the Main_PlotWakeML.cpp but intended to be used 
//by those who do not have a MatLab licence and with the added feature of quiver plots.
//
// COMPILING - Need to put the executable version of this program in the output folder
//			 on linux: g++ Source/Main_PlotQuiverPY.cpp -o output/Quiver
//
// D.F.B. in Braunschweig, Germany, Mar. 2020
	//
	//update: the file plots _TDVE##.txt files.

//General algorithm
//	1. user input of timestep and width intervals that are to be plotted
//	2. reading in data of timestep
//	3. writing to plotting file (with extension .py)
//  4. running the python script in system command 

	// TOGGLE WAKE PLOTTING
	bool plotWAKE = 1; // Turn on and off the wake plotting

	int timestep;					//timestep whose wake is being plotted
	int intervall=1;				//time intervalls being plotted, default is every one
	int tmin,cutoff;				//minum and maximum time intervalls
	int nospan;						//number of elements in span direction
	int nochord;					//number of elements in chord direction of surface
	int nosurface;					//number of elements of surface
	int n,span,time;				// loop counter
	double x1[3],x2[3],x3[3],x4[3];	//corner points
	double tempA[3],tempAA[3];
	char iofile[125];				//input-output-file
	char ch;						//generic character
	FILE *fp;						//output file
	FILE *fp2;						//quiver files
	int j,i;						// generic counters
	double vec[3],pos[3]; 			// xyz and uvw of quiver 

//	1. user input of timestep and intervalls that is to be plotted
	printf("\nThis is the quiver-plotting program for %s\n\n",PROGRAM_VERSION);
	printf("The timestep and quiver files needs to be located in the ");
	printf("current directory.\n");
	printf("The user must have python installed as this used matplotlib.\n\n");
	printf("Please enter full name of _TDVE file ");
	printf("that is to be plotted: ");
	scanf("%s",iofile);

	printf("Plot wake? 1 for yes; 0 for not ");
	scanf("%d",&n);

	plotWAKE = bool(n);

	printf(" moving on %d\n", plotWAKE);

//	printf("Please enter number of timestep ");
//	printf("whose wake needs to be plotted: ");
//	scanf("%d",&timestep);

	//printf("\nWhat time intervalls are desired? ");
//	scanf("%d",&intervall);
	intervall = 1;

	//printf("\nHow far downstream do you want to plot? (0 is first timestep) ");
//	scanf("%d",&cutoff);
	cutoff = 1;

//	printf("\nHow close to the trailing edge? ");
//	printf("(%d is up to trailing edge) ",timestep);
//	scanf("%d",&tmin);
	tmin = timestep;

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
	tmin = timestep;

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

	//allocates memory for surfacePtr
	ALLOC1D(&surfacePtr,nosurface);

	for(n=0;n<nosurface;n++)
	{
		//read ref. point, nu, epsilon, sweeps, local U,span, chord
		fscanf(fp,"%*d %lf %lf %lf \
					%*lf %*lf %*lf %*lf \
				   %*lf %*lf %*lf \
				   %lf %lf \
				   %lf %lf %lf \
				   %lf %lf %lf",\
		   &surfacePtr[n].xo[0],&surfacePtr[n].xo[1],&surfacePtr[n].xo[2],\
		   &surfacePtr[n].eta,&surfacePtr[n].xsi,\
		   &surfacePtr[n].nu,&surfacePtr[n].epsilon,&surfacePtr[n].psi,\
		   &surfacePtr[n].phiLE,&surfacePtr[n].phi0,&surfacePtr[n].phiTE);

			//convert angles to radians
		surfacePtr[n].nu 		*= DtR;
		surfacePtr[n].epsilon 	*= DtR;
		surfacePtr[n].psi		*= DtR;
		surfacePtr[n].phiLE		*= DtR;
		surfacePtr[n].phi0		*= DtR;
		surfacePtr[n].phiTE		*= DtR;

		//find the beginning of next span information
		do	ch = fgetc(fp);
		while (ch!='\n');
	}

	printf(" Done reading in SDVE\n");
	//find the second '#'-sign at end of header of wake elements
	do ch = fgetc(fp);
	while (ch!='#');

    //allocates memory for wakePtr
	ALLOC2D(&wakePtr,timestep+1,nospan);

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
					   %lf %lf\
					   %*lf \
					   %lf %lf %lf",\
				   &wakePtr[time][span].xo[0],&wakePtr[time][span].xo[1],&wakePtr[time][span].xo[2],\
				   &wakePtr[time][span].nu,&wakePtr[time][span].epsilon,&wakePtr[time][span].psi,\
				   &wakePtr[time][span].u[0],&wakePtr[time][span].u[1],&wakePtr[time][span].u[2],\
				   &wakePtr[time][span].A,&wakePtr[time][span].B,&wakePtr[time][span].C,\
				   &wakePtr[time][span].eta,&wakePtr[time][span].xsi,\
				   &wakePtr[time][span].phiLE,&wakePtr[time][span].phi0,&wakePtr[time][span].phiTE);


			//convert angles to radians
			wakePtr[time][span].nu 		*= DtR;
			wakePtr[time][span].epsilon *= DtR;
			wakePtr[time][span].psi		*= DtR;
			wakePtr[time][span].phiLE	*= DtR;
			wakePtr[time][span].phi0	*= DtR;
			wakePtr[time][span].phiTE	*= DtR;

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


//	3. writing to plotting file (with extension .py)

sprintf(iofile,"%s","wakeplot.py");
//sprintf(iofile,"%s","wakeplot.m");
	//opens input file
	fp = fopen(iofile, "w");

// writing Header
	fprintf(fp,"#Plotting wake results that were generated with ");
	fprintf(fp,"#%s\n",PROGRAM_VERSION);
	//comment
	fprintf(fp,"#Plotting wake of timestep %d\n\n",timestep);
	
// importing required librarys
	fprintf(fp,"import matplotlib as mpl\n");
	fprintf(fp,"from mpl_toolkits.mplot3d import Axes3D\n");
	fprintf(fp,"import numpy as np\n");
	fprintf(fp,"import matplotlib.pyplot as plt\n");

// setup an "set_axes_equal" function
	fprintf(fp,"\ndef set_axes_equal(ax):\n");
	fprintf(fp,"\tx_limits = ax.get_xlim3d()\n");
	fprintf(fp,"\ty_limits = ax.get_ylim3d()\n");
	fprintf(fp,"\tz_limits = ax.get_zlim3d()\n");	
	fprintf(fp,"\tx_range = abs(x_limits[1] - x_limits[0])\n");
	fprintf(fp,"\tx_middle = np.mean(x_limits)\n");	
	fprintf(fp,"\ty_range = abs(y_limits[1] - y_limits[0])\n");
	fprintf(fp,"\ty_middle = np.mean(y_limits)\n");
	fprintf(fp,"\tz_range = abs(z_limits[1] - z_limits[0])\n");
	fprintf(fp,"\tz_middle = np.mean(z_limits)\n");
	fprintf(fp,"\tplot_radius = 0.5*max([x_range, y_range, z_range])\n");
	fprintf(fp,"\tax.set_xlim3d([x_middle - plot_radius, x_middle + plot_radius])\n");
	fprintf(fp,"\tax.set_ylim3d([y_middle - plot_radius, y_middle + plot_radius])\n");
	fprintf(fp,"\tax.set_zlim3d([z_middle - plot_radius, z_middle + plot_radius])\n\n");

// setup plots

	fprintf(fp,"mpl.rcParams['legend.fontsize'] = 10\n");
	fprintf(fp,"fig = plt.figure()\n");
	fprintf(fp,"ax = fig.gca(projection='3d')\n");
	fprintf(fp,"ax.set_aspect('equal')\n\n");

	// Plot labels 
	fprintf(fp,"ax.set_xlabel('Global X')\n");
	fprintf(fp,"ax.set_ylabel('Global Y')\n");
	fprintf(fp,"ax.set_zlabel('Global Z')\n\n");

if(plotWAKE){
//Plottine Wake
	//computing the corner points for each DVE and plotting them together.
	for(time=cutoff;time<=tmin;time+=intervall)
	{
		for(span=0;span<nospan;span++)
		{

//printf("  %d %d\n",time,span);	
			//computing left-leading edge point in local ref. frame
			tempA[0] = -wakePtr[time][span].xsi\
					 - wakePtr[time][span].eta*tan(wakePtr[time][span].phiLE);
			tempA[1] = -wakePtr[time][span].eta;
			tempA[2] = 0;

			Star_Glob(tempA,wakePtr[time][span].nu,\
				wakePtr[time][span].epsilon,wakePtr[time][span].psi,tempAA);
			vsum(tempAA,wakePtr[time][span].xo,x1);


			//computing left-trailing edge point in local ref. frame
			tempA[0] = wakePtr[time][span].xsi\
					 - wakePtr[time][span].eta*tan(wakePtr[time][span].phiTE);
			tempA[1] = -wakePtr[time][span].eta;
			tempA[2] = 0;

			Star_Glob(tempA,wakePtr[time][span].nu,\
				wakePtr[time][span].epsilon,wakePtr[time][span].psi,tempAA);
			vsum(tempAA,wakePtr[time][span].xo,x2);


			//computing right-trailing edge point in local ref. frame
			tempA[0] = wakePtr[time][span].xsi\
					 + wakePtr[time][span].eta*tan(wakePtr[time][span].phiTE);
			tempA[1] = wakePtr[time][span].eta;
			tempA[2] = 0;

			Star_Glob(tempA,wakePtr[time][span].nu,\
				wakePtr[time][span].epsilon,wakePtr[time][span].psi,tempAA);
			vsum(tempAA,wakePtr[time][span].xo,x3);


			//computing right-leading edge point in local ref. frame
			tempA[0] = -wakePtr[time][span].xsi\
					 + wakePtr[time][span].eta*tan(wakePtr[time][span].phiLE);
			tempA[1] = wakePtr[time][span].eta;
			tempA[2] = 0;

			Star_Glob(tempA,wakePtr[time][span].nu,\
				wakePtr[time][span].epsilon,wakePtr[time][span].psi,tempAA);
			vsum(tempAA,wakePtr[time][span].xo,x4);


			//print plot coordinates
			fprintf(fp,"ax.plot("); //first plot command
//#			fprintf(fp,"plot("); //first plot command
			fprintf(fp,"[%lf,%lf,%lf,%lf,%lf],",x1[0],x2[0],x3[0],x4[0],x1[0]); //x
			fprintf(fp,"[%lf,%lf,%lf,%lf,%lf],",x1[1],x2[1],x3[1],x4[1],x1[1]); //y
			fprintf(fp,"[%lf,%lf,%lf,%lf,%lf],",x1[2],x2[2],x3[2],x4[2],x1[2]); //z
			fprintf(fp,"'r')\n");//next line
		}
	}
}
//Plotting wing
	for(n=0;n<nosurface;n++)
	{
		//computing left-leading edge point in local ref. frame
		tempA[0] = -surfacePtr[n].xsi\
				 - surfacePtr[n].eta*tan(surfacePtr[n].phiLE);
		tempA[1] = -surfacePtr[n].eta;
		tempA[2] = 0;

		Star_Glob(tempA,surfacePtr[n].nu,\
			surfacePtr[n].epsilon,surfacePtr[n].psi,tempAA);
		vsum(tempAA,surfacePtr[n].xo,x1);


		//computing left-trailing edge point in local ref. frame
		tempA[0] = surfacePtr[n].xsi\
				 - surfacePtr[n].eta*tan(surfacePtr[n].phiTE);
		tempA[1] = -surfacePtr[n].eta;
		tempA[2] = 0;

		Star_Glob(tempA,surfacePtr[n].nu,\
			surfacePtr[n].epsilon,surfacePtr[n].psi,tempAA);
		vsum(tempAA,surfacePtr[n].xo,x2);


		//computing right-trailing edge point in local ref. frame
		tempA[0] = surfacePtr[n].xsi\
				 + surfacePtr[n].eta*tan(surfacePtr[n].phiTE);
		tempA[1] = surfacePtr[n].eta;
		tempA[2] = 0;

		Star_Glob(tempA,surfacePtr[n].nu,\
			surfacePtr[n].epsilon,surfacePtr[n].psi,tempAA);
		vsum(tempAA,surfacePtr[n].xo,x3);


		//computing right-leading edge point in local ref. frame
		tempA[0] = -surfacePtr[n].xsi\
				 + surfacePtr[n].eta*tan(surfacePtr[n].phiLE);
		tempA[1] = surfacePtr[n].eta;
		tempA[2] = 0;

		Star_Glob(tempA,surfacePtr[n].nu,\
			surfacePtr[n].epsilon,surfacePtr[n].psi,tempAA);
		vsum(tempAA,surfacePtr[n].xo,x4);


		//print plot coordinates
		fprintf(fp,"ax.plot("); //first plot command
//#		fprintf(fp,"plot("); //first plot command
		fprintf(fp,"[%lf,%lf,%lf,%lf,%lf],",x1[0],x2[0],x3[0],x4[0],x1[0]); //x
		fprintf(fp,"[%lf,%lf,%lf,%lf,%lf],",x1[1],x2[1],x3[1],x4[1],x1[1]); //y
		fprintf(fp,"[%lf,%lf,%lf,%lf,%lf],",x1[2],x2[2],x3[2],x4[2],x1[2]); //z
		fprintf(fp,"'k')\n");//next line
	}

/***************************** ADDITION FROM Main_PlotWakeML.cpp ************************
	// Add quiver plot
	// checks if input file exists
	if ((fp2 = fopen("quiver1.txt", "r"))== NULL)
	{
		printf("File could not open quiver1.txt:\n");
	//	exit(1);
	}
	else  //of quiver.txt exists, plot normal forces
	{
		j = 0;

		//opens input file
		fp2 = fopen("quiver1.txt", "r");
		do	{
			ch = fgetc(fp2);
			if(ch == '\n'){
			j++;}
			}
		while (ch!=EOF);
		fclose(fp2);
		double temp;
		fp2 = fopen("quiver1.txt", "r");
		for(i=0;i<j;i++){
			fscanf(fp2,"%lf %lf %lf %lf %lf %lf",&pos[0],&pos[1],&pos[2],&vec[0],&vec[1],&vec[2]);
			fprintf(fp,"ax.quiver(%lf,%lf,%lf,%lf,%lf,%lf)\n",pos[0],pos[1],pos[2],vec[0],vec[1],vec[2]);
			do	ch = fgetc(fp2);
			while (ch!='\n');
		}
		fclose(fp2);
	}
	/////////////////////*/
	fprintf(fp,"ax.quiver(%lf,%lf,%lf,%lf,%lf,%lf)\n",pos[0],pos[1],pos[2],vec[0],vec[1],vec[2]);

	// Python lines to show plot and run axis equal
	fprintf(fp,"\nset_axes_equal(ax)\n");
	fprintf(fp,"plt.show()\n");

	fclose(fp);

	//allocates memory for surfacePtr
	FREE1D(&surfacePtr,nosurface);

	//allocates memory for wakePtr
	FREE2D(&wakePtr,timestep+1,nospan);

	// run python script
	system("python wakeplot.py");

	printf("\nIt's all done here folks.\n");


return(0);
}
//===================================================================//
		//END of program
//===================================================================//
