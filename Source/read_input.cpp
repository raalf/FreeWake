//reads general information from input file
//void General_Info_from_File(double &,double &,double &,double [3],int &,\
//							double &,double &,double &,double &,int &,int &,\
//							int &, int &, int &,int &,int &);
void General_Info_from_File(GENERAL &,double &,double &,double &);

//reads panel information from input file
void Panel_Info_from_File(PANEL *, const GENERAL);
//reads in V-tail and fuselage, interference drag information
void VT_Fus_Info(GENERAL &,double [5],double [5], int [5],\
				double [20],double &d,int &,double &I);

//This subroutine reads in the data of a particular timestep
void Read_Timestep(int,DVE *,DVE **);

//===================================================================//
		//START OF General_Info_from_File
//===================================================================//
//void General_Info_from_File(double &Uinf,double &alpha, double &beta, \
//							double U[3],int &maxtime,double &deltime,\
//							double &deltae,double &S,double &b,int &steady,\
//							int &linear,int &sym,int &relax,int &nowing,\
//							int &m,int &nopanel)
void General_Info_from_File(GENERAL &info,double &alpha1,double &alpha2,double &alphastep)
{
	//Function 'General_Info_from_File'
	//reads general information from input file:
	//free stream velocity, angle of attack, sideslip angle, ref. area
	//symmetry flag, number of chordwise lifting lines, number of panels

	FILE *fp;		//input file

	char ch;		//generic character
	int tempI;		//temporary integer

	// checks if input file exists
	if ((fp = fopen("input.txt", "r"))== NULL) {
		printf("Input file could not be opened:\n");
		scanf("%c",&ch);
		exit(1);
	}

	//opens input file
	fp = fopen("input.txt", "r");

//=============================================================================

	// read relaxed wake flag
	//(relax =1 -> wake is relaxed, =0 -> it's not)
	//find the '='-sign in input file before sym
	do	
	ch = fgetc(fp);
	while (ch!='=');
	//reads relaxed wake flag
	fscanf(fp,"%d", &info.relax);

	// read steady (=1)/unsteady (=0) aerodynamics flag
	//find the '='-sign in input file before steady flag
	do	ch = fgetc(fp);
	while (ch!='=');
	//reads flag for steady/unsteady aerodynamics
	fscanf(fp,"%d", &info.steady);

	//read viscous (=1)/inviscid (=0) flow flag  added GB 2-18-20
	//find the '='-sign in input file before steady flag
	do	ch = fgetc(fp);
	while (ch!='=');
	//reads flag for viscous aerodynamics
	fscanf(fp,"%d",&tempI);
	info.flagVISCOUS=tempI;

	// read symmetrical geometry flag
	//(sym =1 -> symmetrical conditions, =0 -> asymmetrical)
	//find the '='-sign in input file before sym
	do	ch = fgetc(fp);
	while (ch!='=');
	//reads symmetry flag
	fscanf(fp,"%d", &info.sym);

	// linear theory flag
	//(lin =1 -> linear theory is being applied, =0 -> it's not)
	info.linear = 0;

//=============================================================================
	// read longitudinal trim flag
	//(trim =1 -> aircraft is trimmed m = 1!, trim=0 -> no trim)
	//find the '='-sign in input file before sym
	do	ch = fgetc(fp);
	while (ch!='=');
	//reads symmetry flag
	fscanf(fp,"%d", &info.trim);

	//read max. number of time steps
	//find the '='-sign in input file before maxtime
	do	ch = fgetc(fp);
	while (ch!='=');
	//reads span
	fscanf(fp,"%d", &info.maxtime);

	//read width of time steps
	//find the '='-sign in input file before deltime
	do	ch = fgetc(fp);
	while (ch!='=');
	//reads width of time steps
	fscanf(fp,"%lf", &info.deltime);

	//read deltae
	//find the '='-sign in input file before dele
	do	ch = fgetc(fp);
	while (ch!='=');
	//reads deltae
	fscanf(fp,"%lf", &info.deltae);
	info.deltae = info.deltae*info.deltae; //the square is needed in the do-while loop

//=============================================================================

	// read free stream velocity
	//find the '='-sign in input file before Uinf
	do	ch = fgetc(fp);
	while (ch!='=');
	//reads free stream velocity
	fscanf(fp,"%lf", &info.Uinf);

	// read angle of attack
	//find the '='-sign in input file before alpha
	do	ch = fgetc(fp);
	while (ch!='=');
	//reads angle of attack
//old, changed 8-8-07 G.B.
//	fscanf(fp,"%lf", &alpha);
//	alpha *=DtR;	//changes deg. to radians
	fscanf(fp,"%lf", &alpha1);	
	fscanf(fp,"%lf", &alpha2);
	fscanf(fp,"%lf", &alphastep);
	alpha1 *=DtR;	alpha2 *=DtR;	alphastep *=DtR;	//changes deg. to radians
	
	//read sideslip angle
	//find the '='-sign in input file before beta
	do	ch = fgetc(fp);
	while (ch!='=');
	//reads sideslip
	fscanf(fp,"%lf", &info.beta);
	info.beta *=DtR;	//changes deg. to radians

	//computes free stream velocity vector
	info.U[0]=info.Uinf*cos(alpha1)*cos(info.beta);
	info.U[1]=info.Uinf			*sin(info.beta);
	info.U[2]=info.Uinf*sin(alpha1)*cos(info.beta);
	//computes free stream velocity vector

	//read density
	//find the '='-sign in input file before beta
	do	ch = fgetc(fp);
	while (ch!='=');
	//reads density
	fscanf(fp,"%lf", &info.density);

	//read kinematic viscosity
	//find the '='-sign in input file before beta
	do	ch = fgetc(fp);
	while (ch!='=');
	//reads kinematic viscosity
	fscanf(fp,"%lf", &info.nu);

//=============================================================================

	// read reference area
	//find the '='-sign in input file before S
	do	ch = fgetc(fp);
	while (ch!='=');
	//reads reference area
	fscanf(fp,"%lf", &info.S);

	// read span
	//find the '='-sign in input file before b
	do	ch = fgetc(fp);
	while (ch!='=');
	//reads span
	fscanf(fp,"%lf", &info.b);

	// read mean aerodynamic chord
	//find the '='-sign in input file before b
	do	ch = fgetc(fp);
	while (ch!='=');
	//reads cmac
	fscanf(fp,"%lf", &info.cmac);

	// read aircraft weight
	//find the '='-sign in input file before b
	do	ch = fgetc(fp);
	while (ch!='=');
	//reads aircraft weight
	fscanf(fp,"%lf", &info.W);
	
	// read CG-location
	//find the '='-sign in input file before b
	do	ch = fgetc(fp);
	while (ch!='=');
	//reads CG-location
	fscanf(fp,"%lf", &info.RefPt[0]);
	fscanf(fp,"%lf", &info.RefPt[1]);
	fscanf(fp,"%lf", &info.RefPt[2]);

	// read zero lift moment of wing
	//find the '='-sign in input file before b
	do	ch = fgetc(fp);
	while (ch!='=');
	//reads CMoWing
	fscanf(fp,"%lf", &info.CMoWing);


//=============================================================================

	// read number of wings
	//find the '='-sign in input file before nowing
	do	ch = fgetc(fp);
	while (ch!='=');
	//reads
	fscanf(fp,"%d", &info.nowing);

	// read number of panels
	//find the '='-sign in input file before nopanel
	do	ch = fgetc(fp);
	while (ch!='=');
	//reads
	fscanf(fp,"%d", &info.nopanel);

	// read number of chordwise lifting lines
	//find the '='-sign in input file before m
	do	ch = fgetc(fp);
	while (ch!='=');
	//reads number of chordwise lifting lines
	fscanf(fp,"%d", &info.m);
	if(info.trim==1)   //longitudinal trim routine
	{
		info.m=1;
	  printf("Running Longitudinal Trim Version \n ");
	  printf("\t\t\t\tThat means m=1!!\n");
	  printf("\n\n\t\t\t!!!ONLY TWO WINGSS!!\n");
	  printf("\t\tThat means one main wing and one horizontal tail!!\n");
	}

	
	// read number of airfoils
	//find the '='-sign in input file before m
	do	ch = fgetc(fp);
	while (ch!='=');
	//reads number of airfoils
	fscanf(fp,"%d", &info.noairfoils);

//=============================================================================
	//closes input file
	fclose(fp);
}
//===================================================================//
		//END OF General_Info_from_File
//===================================================================//

//===================================================================//
		//START OF Panel_Info_from_File
//===================================================================//
void Panel_Info_from_File(PANEL *panelPtr, const GENERAL info)
{
	//Function 'Panel_Info_from_File'
	//reads panel information from input file:
	//leading edge coordinates, chords, incident angles, local free stream
	//velocity delta (!= 0 for rotating wing) and adds Uinf
	//	Panel Boundary Conditions:
	//	First Digit: 	0 - undefined circulation strength
	//					1 - zero circulation (free end)
	//					2 - circulation strength equal to neighboring elementary wing
	//	Second Digit: 	0 - undefined circulation slope
	//					1 - zero slope in circulation
	//					2 - circulation slope equal to neighboring elementary wing
	//	Third Digit: 	0 - undefined circulation curvature
	//					1 - zero slope in circulation change
	//					2 - circulation slope equal to neighboring elementary wing

	FILE *fp;		//input file
	char ch;		//generic character
	int i;			//loop counter


	// checks if input file exists
	if ((fp = fopen("input.txt", "r"))== NULL) {
		printf("Input file could not be opened:\n");
		scanf("%c",&ch);
		exit(1);
	}

	//opens input file
	fp = fopen("input.txt", "r");

	//read in loop over number of panels (info.nopanel)
	for (i=0;i<info.nopanel;i++)
	{
		//find the '#'-sign at beginning of panel info
		do	ch = fgetc(fp);
		while (ch!='#');

		//find the '='-sign in first line of panel info
		do	ch = fgetc(fp);
		while (ch!='=');
		//reads number of spanwise elements n
		fscanf(fp,"%d", &(panelPtr[i].n));

//removed 2-18-20 GB
		//find nest '='-sign in first line of panel info
//		do	ch = fgetc(fp);
//		while (ch!='=');
		//reads number of airfoil used for this panel
//		fscanf(fp,"%d", &(panelPtr[i].airfoil));
//			panelPtr[i].airfoil--; //adjust index
		panelPtr[i].airfoil = -1; //default, should cause error
		
		//find the first ':'-sign in second line of panel info
		do	ch = fgetc(fp);
		while (ch!=':');

		//reads number of panel that borders to the left
		fscanf(fp,"%d", &(panelPtr[i].left));

		//find the second ':'-sign in second line of panel info
		do	ch = fgetc(fp);
		while (ch!=':');
		//reads number of panel that borders to the right
		fscanf(fp,"%d", &(panelPtr[i].right));

		//find the beginning of the second line of panel info
		do	ch = fgetc(fp);
		while (ch!='\n');
		//find the beginning of the third line of panel info
		do	ch = fgetc(fp);
		while (ch!='\n');

		//reads info of panel edge 1:
		//l.e. coordinates, chord, incident angle,
		//local free stream vel., boundary condition.
		fscanf(fp,"%lf %lf %lf %lf %lf %d",\
			&(panelPtr[i].x1[0]),&(panelPtr[i].x1[1]),&(panelPtr[i].x1[2]),\
			&(panelPtr[i].c1),&(panelPtr[i].eps1),\
			&(panelPtr[i].BC1));
		panelPtr[i].eps1 *=DtR;	//changes deg. to radians

		//read in airfoil on edge 1 GB 2-14-20
		fscanf(fp," %d",&(panelPtr[i].airfoil1));
			panelPtr[i].airfoil1--; //adjust index

		panelPtr[i].u1[0]=0; panelPtr[i].u1[1]=0; panelPtr[i].u1[2]=0;

		//find the beginning of the fifth line of panel info
		do	ch = fgetc(fp);
		while (ch!='\n');
		//find the beginning of the sixth line of panel info
		do	ch = fgetc(fp);
		while (ch!='\n');

		//reads info of panel edge 2:
		//l.e. coordinates, chord, incident angle,
		//local free stream vel. variation, boundary condition.
		fscanf(fp,"%lf %lf %lf %lf %lf %d",\
			&(panelPtr[i].x2[0]),&(panelPtr[i].x2[1]),&(panelPtr[i].x2[2]),\
			&(panelPtr[i].c2),&(panelPtr[i].eps2),\
			&(panelPtr[i].BC2));
		panelPtr[i].eps2 *=DtR;						//changes deg. to radians

		//read in airfoil on edge 2 GB 2-14-20
		fscanf(fp," %d",&(panelPtr[i].airfoil2));
			panelPtr[i].airfoil2--; //adjust index

			panelPtr[i].u2[0]=0; panelPtr[i].u2[1]=0; panelPtr[i].u2[2]=0;

		vsum(panelPtr[i].u1,info.U,panelPtr[i].u1);	//adds undisturbed free
		vsum(panelPtr[i].u2,info.U,panelPtr[i].u2);	//stream vel.

	}		//END loop over i

	//closes input file
	fclose(fp);

	//debug output
//	for (i=0;i<info.nopanel;i++)
//	{
//		printf("pnale %d x1 y1 z1 (%lf,%lf,%lf)  c1 %lf eps1 %lf, BC1 %d airf1 %d\n",\
//			i,panelPtr[i].x1[0],panelPtr[i].x1[1],panelPtr[i].x1[2],\
//			panelPtr[i].c1,panelPtr[i].eps1,panelPtr[i].BC1,panelPtr[i].airfoil1);
//		printf("pnale %d x2 y2 z2 (%lf,%lf,%lf)  c2 %lf eps2 %lf, BC2 %d airf2 %d\n",\
//			i,panelPtr[i].x2[0],panelPtr[i].x2[1],panelPtr[i].x2[2],\
//			panelPtr[i].c2,panelPtr[i].eps2,panelPtr[i].BC2,panelPtr[i].airfoil2);
//	}


//#############################################################################
	//determining indiceses of elements at the left and right trailing edge
	//corner of the panels.  added 8/12/05 G.B.

	//the first panel
	panelPtr[0].TE1 = panelPtr[0].n*(info.m-1);			//the left DVE @ TE
	panelPtr[0].TE2 = panelPtr[0].TE1 + panelPtr[0].n-1;//the right DVE @ TE

	//loop over number of panels (info.nopanel)
	for (i=1;i<info.nopanel;i++)
	{
		//the left DVE @ TE
		panelPtr[i].TE1 = panelPtr[i-1].TE2 + 1 + panelPtr[i].n*(info.m-1);
		//the right DVE @ TE
		panelPtr[i].TE2 = panelPtr[i].TE1 + panelPtr[i].n-1;
	}
//#############################################################################
}
//===================================================================//
		//END OF Panel_Info_from_File
//===================================================================//

//===================================================================//
		//START OF VT_Fuselage_Info File
//===================================================================//
void VT_Fus_Info(GENERAL &info,\
		double VTchord[5],double VTarea[5], int VTairfoil[5],\
		double FusSectS[20],double &delFus,int &FusLT,double &IFdrag)
{

	//reads in information of vertical tail and fuselage  from input file:
	//
	//	
	// VTchord[5],VTarea[5];//chord and aera of VT panel
	// VTairfoil[5];		//VT airfoil
	// FusSectS[20];		//fuselage section area
	// delFus;				//width of each section
	// FusLT;				//fuselage section where turbulent
	// IF;					//interfernce drag fraction


	FILE *fp;		//input file
	char ch;		//generic character
	int i;			//loop counter


	// checks if input file exists
	if ((fp = fopen("input.txt", "r"))== NULL) {
		printf("Input file could not be opened:\n");
		scanf("%c",&ch);
		exit(1);
	}

	//opens input file
	fp = fopen("input.txt", "r");
	
	//find the '%'-sign at beginning of V-tail info
	do	ch = fgetc(fp);
	while (ch!='%');

	//===========V-TAIL=============//
	//find the '='-sign before no. of V-tail panels
	do	ch = fgetc(fp);
	while (ch!='=');

	//reads number of V-tail panels
	fscanf(fp,"%d", &(info.noVT));
	//checks if max. no. of V-tail panels is exceeded
	if(info.noVT>5) {printf("\nToo many V-tail panels, max. no. is 5!!");
					exit(0);}

	//find the beginning of the second line of V-tail info
	do	ch = fgetc(fp);
	while (ch!='\n');
	//find the beginning of the third line of V-tail info
	do	ch = fgetc(fp);
	while (ch!='\n');

	//reads in chord, area, and airfoil of each V-tail panel
	for (i=0;i<info.noVT;i++)
	{
		fscanf(fp,"%*d %lf %lf  %d",&VTchord[i],&VTarea[i],&VTairfoil[i]);
		  VTairfoil[i]--;//adjust airfoil index
		//find the beginning of next line
		do	ch = fgetc(fp);
		while (ch!='\n');

	}
	//===========V-TAIL=============//

	//==========FUSELAGE============//
	//find the '='-sign before no. of fuselage sections
	do	ch = fgetc(fp);
	while (ch!='=');

	//reads number of fuselage sections
	fscanf(fp,"%d", &(info.noFus));
	//checks if max. no. of fuselage sections is exceeded
	if(info.noFus>20) {printf("\nToo many V-tail panels, max. no. is 5!!");
					exit(0);}

	//find the '='-sign before width of sections
	do	ch = fgetc(fp);
	while (ch!='=');
	//reads width fuselage sections
	fscanf(fp,"%lf", &(delFus));

	//find the '='-sign before info of flow transition
	do	ch = fgetc(fp);
	while (ch!='=');
	//reads in no. of panel where flow transitions from lam to turb
	fscanf(fp,"%d", &(FusLT)); 
		FusLT--;  //adjusting index

	//find the beginning of the next line
	do	ch = fgetc(fp);
	while (ch!='\n');
	//find the beginning of the next line
	do	ch = fgetc(fp);
	while (ch!='\n');

	//reads fuselage section circumferences and computes area
	for (i=0;i<info.noFus;i++)
	{
		fscanf(fp,"%*d %lf",&FusSectS[i]);  //diameter
			FusSectS[i] *= delFus*Pi;
			
		//find the beginning of the next line
		do	ch = fgetc(fp);
		while (ch!='\n');
	}
	//==========FUSELAGE============//

	//========INTERFERENCE==========//
	//find the '='-sign before the interference drag fraction
	do	ch = fgetc(fp);
	while (ch!='=');
	//reads in interference drag fractin
	fscanf(fp,"%lf", &(IFdrag)); 
		IFdrag *=.01;  //adjusting to fraction
	//========INTERFERENCE==========//

	//closes input file
	fclose(fp);
}
//===================================================================//
		//END OFVT_Fuselage_Info File
//===================================================================//

//===================================================================//
		//START OF Read_Timestep
//===================================================================//
void Read_Timestep(const int timestep,DVE *surfaceDVE,DVE **wakeDVE)

{
//This subroutine reads in the data of a particular timestep
//and assigns them to the appropriate variables
//
// input:
//
//	timestep		timestep that is being read in
//	surfaceDVE		holds info of DVEs that model the lifting surface
//	wakeDVE			holds info of DVEs that model the wake
//
// memory has to be allocated for surfaceDVE and wakeDVE in the function
// that calls Read_Timestep
//
	int nosurface;		//number of DVEs on the lifting surface

	int index,span,time;// loop counter
	char iofile[125];	//input-output-file
	char ch;			//generic character
	FILE *fp;			//output file

	//creates file name timestep##.txt ## is number of timestep
	sprintf(iofile,"%s%s%d%s",OUTPUT_PATH,"timestep",timestep,".txt");

	// checks if input file exists
	if ((fp = fopen(iofile, "r"))== NULL)
	{
		printf("Output file could not be opened:\n");
		scanf("%c",&ch);
		exit(1);
	}

	//opens input file
	fp = fopen(iofile, "r");

	//find the ':'-sign in input file before program version
	do	ch = fgetc(fp);
	while (ch!=':');

	//find the ':'-sign in input file before reference area
	do	ch = fgetc(fp);
	while (ch!=':');
	fscanf(fp,"%lf", &info.S);

	//find the ':'-sign in input file before aspect ratio
	do	ch = fgetc(fp);
	while (ch!=':');
	fscanf(fp,"%lf", &info.AR);

	//find the ':'-sign in input file before alpha
	do	ch = fgetc(fp);
	while (ch!=':');
	fscanf(fp,"%lf", &info.alpha);
	info.alpha *= DtR;

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

	//find the ':'-sign in input file before number of elements along chord
	do	ch = fgetc(fp);
	while (ch!=':');
	//reads nospanelements
	fscanf(fp,"%d", &info.m);

	//number of surface DVEs
	nosurface = info.m*info.nospanelement;

	//find the ':'-sign in input file before number of wings
	do	ch = fgetc(fp);
	while (ch!=':');
	//reads nospanelements
	fscanf(fp,"%d", &info.nowing);

	//find the ':'-sign in input file before number of span indices
	do	ch = fgetc(fp);
	while (ch!=':');
	for(index=0;index<info.nowing;index++)
	{
		fscanf(fp,"%d", &info.wing1[index]);
		fscanf(fp,"%d", &info.wing2[index]);
	}


//	a reading in surface-DVE info of timestep

	//find the '#'-sign at end of header of surface data
	do		ch = fgetc(fp);
	while (ch!='#');

	for(index=0;index<nosurface;index++)
	{
		//reading in the data (values in parentheses are ignored):

//index xo yo zo Nlift_tot NLift_ind NY_tot NYi
		//reading: (index)  xo  yo  zo  (Nlift_tot) (NLift_ind) (NY_tot) (NYi)
		fscanf(fp,"%*d %lf %lf %lf %*lf %*lf %*lf %*lf",\
			   &surfaceDVE[index].xo[0],&surfaceDVE[index].xo[1],\
			   &surfaceDVE[index].xo[2]);

//A B C eta xsi
		//reading: A B C eta xsi
		fscanf(fp,"%lf %lf %lf %lf %lf",\
			   &surfaceDVE[index].A,&surfaceDVE[index].B,\
			   &surfaceDVE[index].C,&surfaceDVE[index].eta,\
			   &surfaceDVE[index].xsi);

//nu epsilon psi phiLE phi0 phiTE
		//reading: nu epsilon psi phiLE phi0 phiTE
		fscanf(fp,"%lf %lf %lf %lf %lf %lf",\
			   &surfaceDVE[index].nu,&surfaceDVE[index].epsilon,\
			   &surfaceDVE[index].psi,&surfaceDVE[index].phiLE,\
			   &surfaceDVE[index].phi0,&surfaceDVE[index].phiTE);

		//convert angles to radians
		surfaceDVE[index].nu 		*= DtR;
		surfaceDVE[index].epsilon 	*= DtR;
		surfaceDVE[index].psi		*= DtR;
		surfaceDVE[index].phiLE		*= DtR;
		surfaceDVE[index].phi0		*= DtR;
		surfaceDVE[index].phiTE		*= DtR;

		//find the beginning of next element information
		do	ch = fgetc(fp);
		while (ch!='\n');
	}

//	reading in wake-DVE info of timestep

	//find the '#'-sign at end of header
	do		ch = fgetc(fp);
	while (ch!='#');

	for(time=0;time<=timestep;time++)
	{
		for(span=0;span<info.nospanelement;span++)
		{
			//reading in the data (values in parentheses are ignored):
 //span time xo yo zo nu epsilon psi U V W A B C eta xsi K phiLE phi0 phiTE
			//(span) (time) xo yo zo
			fscanf(fp,"%*d %*d %lf %lf %lf ",&wakeDVE[time][span].xo[0],\
					&wakeDVE[time][span].xo[1],&wakeDVE[time][span].xo[2]);

			//nu epsilon psi
			fscanf(fp,"%lf %lf %lf",&wakeDVE[time][span].nu,\
				   &wakeDVE[time][span].epsilon,&wakeDVE[time][span].psi);

			//U V W
			fscanf(fp,"%lf %lf %lf",&wakeDVE[time][span].u[0],\
				   &wakeDVE[time][span].u[1],&wakeDVE[time][span].u[2]);

			//A B C
			fscanf(fp,"%lf %lf %lf",&wakeDVE[time][span].A,\
				   &wakeDVE[time][span].B,&wakeDVE[time][span].C);

			//eta xsi K
			fscanf(fp,"%lf %lf %lf",&wakeDVE[time][span].eta,\
				   &wakeDVE[time][span].xsi,&wakeDVE[time][span].K);

			//phiLE phi0 phiTE
			fscanf(fp,"%lf %lf %lf",&wakeDVE[time][span].phiLE,\
				   &wakeDVE[time][span].phi0,&wakeDVE[time][span].phiTE);

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

}
//===================================================================//
		//END OF Read_Timestep
//===================================================================//
