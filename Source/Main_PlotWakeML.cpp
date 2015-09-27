#include "general.h"
#include "PerfCode.h"
main()
{
//This program plots the wake of a wake step.
//	1. user input of timestep and width intervalls that are to be plotted
//	   furthermore, input of how far downstream
//	2. reading in data of timestep
//	3. writing to plotting file (with extension .wrl)
//  4. displays file in VRML

	int timestep;		//timestep whose wake is being plotted
	int intervall=1;		//time intervalls being plotted, default is every one
	int tmin,cutoff;			//minum and maximum time intervalls
	int nospan;			//number of elements in span direction
	int nochord;		//number of elements in chord direction of surface
	int nosurface;		//number of elements of surface
	int n,span,time;		// loop counter
	double x1[3],x2[3],x3[3],x4[3];	//corner points
	double tempA[3],tempAA[3];
	char iofile[125];	//input-output-file
	char ch;			//generic character
	FILE *fp;			//output file

//	1. user input of timestep and intervalls that is to be plotted
	printf("This is the wake-plotting program for %s\n\n",PROGRAM_VERSION);
	printf("The timestep file needs to be located in the ");
	printf("following directory:\n");
	printf("%s (<- loacated in current directory "," ");//OUTPUT_PATH);
	printf("if just a space)\n");
	printf("Please enter number of timestep");
	printf(" whose wake needs to be plotted: ");
	scanf("%d",&timestep);

	printf("\nWhat time intervalls are desired? ");
//	scanf("%d",&intervall);
intervall = 1;

	printf("\nHow far downstream do you want to plot? (0 is first timestep) ");
//	scanf("%d",&cutoff);
cutoff = 1;

	printf("\nHow close to the trailing edge? ");
	printf("(%d is up to trailing edge) ",timestep);
//	scanf("%d",&tmin);
tmin = timestep;

//	2. reading in data of timestep

	//creates file name timestep##.txt ## is number of timestep
//	sprintf(iofile,"%s%s%d%s",OUTPUT_PATH,"timestep",timestep,".txt");
	sprintf(iofile,"%s%d%s","timestep",timestep,".txt");

	// checks if input file exists
	if ((fp = fopen(iofile, "r"))== NULL)
	{
		printf("File could not be opened, stupid:\n");
		exit(1);
	}

	//opens input file
	fp = fopen(iofile, "r");

	//find the ':'-sign in input file before program version
	do	ch = fgetc(fp);
	while (ch!=':');

	//find the ':'-sign in input file before ref. area
	do	ch = fgetc(fp);
	while (ch!=':');

	//find the ':'-sign in input file before aspect ratio
	do	ch = fgetc(fp);
	while (ch!=':');

	//find the ':'-sign in input file before angle of attack
	do	ch = fgetc(fp);
	while (ch!=':');

	//find the ':'-sign in input file before sideslip angle
	do	ch = fgetc(fp);
	while (ch!=':');

	//find the ':'-sign in input file before symmetry condition
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

	//find the ':'-sign in input file before number of span-elements
	do	ch = fgetc(fp);
	while (ch!=':');
	fscanf(fp,"%d", &nochord);

	nosurface = nospan*nochord;  //no. of surface DVEs

	//find the first '#'-sign at end of header of surface elements
	do	ch = fgetc(fp);
	while (ch!='#');


	//allocates memory for surfacePtr
	ALLOC1D(&surfacePtr,nosurface);

	for(n=0;n<nosurface;n++)
	{
		//read ref. point, nu, epsilon, sweeps, local U,span, chord
		fscanf(fp,"%*d %lf %lf %lf %*lf %*lf %*lf \
				   %*lf %*lf %*lf %*lf %lf %lf %lf %lf\
				   %lf %lf %lf %lf",\
		   &surfacePtr[n].xo[0],&surfacePtr[n].xo[1],&surfacePtr[n].xo[2],\
		   &surfacePtr[n].eta,&surfacePtr[n].xsi,
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
			fscanf(fp,"%*d %*d %lf %lf %lf %lf %lf %lf \
					   %lf %lf %lf %lf %lf %lf %lf %lf\
					   %*lf %lf %lf %lf",\
				   &wakePtr[time][span].xo[0],&wakePtr[time][span].xo[1],\
				   &wakePtr[time][span].xo[2],&wakePtr[time][span].nu,\
				   &wakePtr[time][span].epsilon,&wakePtr[time][span].psi,\
				   &wakePtr[time][span].u[0],&wakePtr[time][span].u[1],\
				   &wakePtr[time][span].u[2],&wakePtr[time][span].A,\
				   &wakePtr[time][span].B,&wakePtr[time][span].C,\
				   &wakePtr[time][span].eta,&wakePtr[time][span].xsi,\
				   &wakePtr[time][span].phiLE,&wakePtr[time][span].phi0,\
				   &wakePtr[time][span].phiTE);


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

//	3. writing to plotting file (with extension .m)

sprintf(iofile,"%s","wakeplot.m");
//sprintf(iofile,"%s","wakeplot.m");
	//opens input file
	fp = fopen(iofile, "w");

//writing header
	fprintf(fp,"%%\n"); //madatory header line
	fprintf(fp,"%%Plotting wake results that were generated with ");
	fprintf(fp,"%s\n",PROGRAM_VERSION);
	//comment
	fprintf(fp,"%%plotting wake of timestep %d\n\n",timestep);
	fprintf(fp,"hold off\n");

//Plottine Wake
	//computing the corner points for each DVE and plotting them together.
	for(time=cutoff;time<=tmin;time+=intervall)
	{
		for(span=0;span<nospan;span++)
		{
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
			fprintf(fp,"plot3("); //first plot command
//#			fprintf(fp,"plot("); //first plot command
			fprintf(fp,"[%lf,%lf,%lf,%lf,%lf],",x1[0],x2[0],x3[0],x4[0],x1[0]); //x
			fprintf(fp,"[%lf,%lf,%lf,%lf,%lf],",x1[1],x2[1],x3[1],x4[1],x1[1]); //y
			fprintf(fp,"[%lf,%lf,%lf,%lf,%lf],",x1[2],x2[2],x3[2],x4[2],x1[2]); //z
			fprintf(fp,"'r')\nhold on\n");//next line
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
		fprintf(fp,"plot3("); //first plot command
//#		fprintf(fp,"plot("); //first plot command
		fprintf(fp,"[%lf,%lf,%lf,%lf,%lf],",x1[0],x2[0],x3[0],x4[0],x1[0]); //x
		fprintf(fp,"[%lf,%lf,%lf,%lf,%lf],",x1[1],x2[1],x3[1],x4[1],x1[1]); //y
		fprintf(fp,"[%lf,%lf,%lf,%lf,%lf],",x1[2],x2[2],x3[2],x4[2],x1[2]); //z
		fprintf(fp,"'k')\nhold on\n");//next line
	}


	fclose(fp);


	//allocates memory for surfacePtr
	FREE1D(&surfacePtr,nosurface);

	//allocates memory for wakePtr
	FREE2D(&wakePtr,timestep+1,nospan);

	printf("\nIt's done.  Enter anything, anything: ");
	scanf("%d",&timestep);

return(0);
}
//===================================================================//
		//END of program
//===================================================================//
