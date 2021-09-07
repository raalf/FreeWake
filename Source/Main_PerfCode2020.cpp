#define _CRT_SECURE_NO_WARNINGS
#include "general.h"
#include "PerfCode.h"


int i,ii,a,a2;		//loop counters, max AOA increment
	int k,l,m;			//loop counters
	
	//drag variables
	double cd;				//section profile drag force coefficients
	double cm,CMo;			//section, total zero-lift moment coefficient
	double cd1,cd2,cm1,cm2;	//section values of panel edge 1 and 2 for interpolation; GB 2-14-20
	double CL,CY,CD,CDi=0;	//aircraft lift, total and induced drag coefficients
	double CLi, CYi = 0;	//induced aircraft lift and side force
	double CLinviscid=0;	//inviscid CL without stall correction
	double CDprofile=0;		//wing profile drag coefficient
	double CDfuselage=0;	//fuselage drag coefficeing
	double CDht=0,CDvt=0;	//horizontal and vertical tail drag coefficeing
	double CLtarget = 0;	//for trim lift
	double CLtemp[3];  //for storing temp lift for trim
	double alpha2old; //for lift trim
	double alphatemp1, alphatemp2, alphasteptemp;
	//coefficients with capital letter, e.g. CL, are wrt to wing area
	
	double D=0, Di=0;		//total and induced drag
	double Dprofile=0;		//wing profile 
	double Dfuselage=0;		//fuselage drags
	double Dvt=0,Dht=0;		//vertical tail,horizontal tial drags
	double Dint=0,Dmisc=0;	//interference, miscillanuous drags
	
	double V_inf=0,q_inf=0; //free stream velocity and dynamic pressure
	double Re;				//Reynolds number

	// Camber related variables
	double ***camberPtr; 	// Point for the 3D array of camber data
	int cambRow = 0;		//Counter for rows of camber data
	int cambCol = 0;		//Counter for column of camber data
	// Airfoil related variables
	double ***airfoilPtr; 	// Point for the 3D array of camber data
	int airfoilRow = 0;		//Counter for rows of camber data
	int airfoilCol = 0;		//Counter for column of camber data
	int airfoil=0;			//airfoil index

	//V-tail and fuselage	
	//max of 5 VT panesl!!
	double VTchord[5],VTarea[5];	//chord and aera of VT panel
	int VTairfoil[5];
	//max of 20 fuselage sections
	double FusSectS[20];	//fuselage section area
	double delFus;			//width of each section
	int FusLT;				//fuselage section where turbulent
	double IFdrag;			//interfernce drag fraction


	int n;					//loop counter
	long pos; 				//file position counter
	double tempS;			//temp scalar
	char answer,string[6];	//string is used to update config file
	char filename[160];		//file path and name

	//Input/output files
	FILE *MomSol;			//output file of trim solutions
	FILE *Performance;		//output file for performance results
    FILE *FltConfg;         //output file for flight configuration results




int main(int argc, char *argv[])

{
// This program computes the drag of an aircraft configuration.
//
// Step 0:	read in aircraft configuration, paneling, and airfoil info
//
// Step 1:	compute trim for given angle of attack and determine
//			inviscid induced drag and normal force distribution
//			in LongitudinalTrim.cpp.  Save cn- and cd-span distribution 
//			in extra file
//
// Step 2:	based on inviscid cn-distribution, compute wing and 
//			horizontal tail profile drags
//
// Step 2a:	update wing zero-lift moment for next AOA (works if changes
//			of zero-lift moment are small with change of AOA)
//
// Step 3:	compute vertical tail profile drag
//
// Step 4:	compute fuselage drag
//
// Step 5:	compute cooling and miscillaneuous drags
//
// Step 6:	compute interference drag
//
// Step7:	save results in output file
//
//
//  Wing 1 is the main wing, Wing 2 is the horizontal tail.
//
//ATTENTION: NO SPANWISE VELOCITY VARIATION !!!
//If you don't know what this means, don't worry and move on.
//
//Update December 2013: 
//   - changed profile drag rountine, now m>1 works for multiple panels.
//   - added trim flag for input file.
//	 - fixed pitching moment routine
//
	printf("Starting FreeWake\n");
	printf("===========================================================================\n ");
	printf("\n\tPerformance Code based on FreeWake 2015 (includes stall model)\n\n");
	printf("\t\tRunning Version %s \n ",PROGRAM_VERSION);
	printf("===========================================================================\n ");

	//info.flagVISCOUS = 1; // Turn on/off viscous corrections (O = OFF, 1 = ON)


//===================================================================//
		//START read general and panel info from file 'input.txt'
//===================================================================//
	//reads general information from input file:
	//free stream velocity, angle of attack, sideslip angle, ref. area
	//symmetry flag, number of chordwise lifting lines, number of panels
	//	Uinf		-free stream velocity
	//	alpha		-angle of attack
	//	beta		-sideslip angle
	//  U			-free stream velocity vector
	//	maxtime		-maximal number of time steps
	//	deltime		-time step width
	//	S			-reference area
	// 	b			-reference wing span
	//  steady		-steady(=1)/unsteady (default) aerodynamics flag
	// 	sym			-symmetrical geometry (=1 symmetrical; 0= assymmetrical)
	//	linear		-1 = linear, small approximation flag, 0 = not
	//	relax		-1 = relaxed wake, 0 = fixed wake
	//	nowing		-number of separate wings
	//	m			-number lifting lines/elementary wings in chord direction
	// 	nopanel		-number of panels
	//	part of variable 'info' of type GENERAL
//	General_Info_from_File(info.Uinf,info.alpha,info.beta,info.U,info.maxtime,\
//						   info.deltime,info.deltae,info.S,info.b,info.steady,\
//						   info.linear,info.sym,info.relax,info.nowing,\
//						   info.m,info.nopanel);
//						   //Subroutine in read_input.cpp



// Define   1. input filename, 2. output directory, 3. configuration file
    
	if (argc < 2) //argc is the number of inputs
	{//no input name: assume input.txt in the working dir
		//changed to save default input.txt results to output/inputs/... D.F.B. 5-2021
        sprintf(info.inputfilename, "%s.txt", "input");
        sprintf(info.output, "%s%s/", OUTPUT_PATH,"input");
        sprintf(info.config, "%s%s_cfg.txt", info.output,"input");
	}
	else if(argc ==2)//with an input filename: use it. 
	{
		 //argv[0] will be the .exe, argv[1] will be the filename
		sprintf(info.inputfilename, "%s.txt", argv[1]);
        sprintf(info.output, "%s%s/", OUTPUT_PATH,argv[1]);
        sprintf(info.config, "%s%s_cfg.txt", info.output,argv[1]);



	}
	else
	{
		printf("Incorrect inputs passed to exe. Format should be:\n");
		printf("FreeWake2020.exe input_filename\n");
		scanf("%c", &answer);
		exit(1);
	}

	// mkdir is moved out of if (argc == 2) statement D.F.B. 5-2021 
		int nError = 0; //check which OS we are compiling on
#if defined(_WIN32)
		nError=_mkdir(info.output);//create output directory in output/
#else 
		nError=mkdir(info.output, 0777); // can be used on non-Windows
#endif
		/*if (nError != 0) {
			printf("error creating output directory\n"); 
			scanf("%c", &answer); 
			exit(1);
		}*/

    printf("\nInput Filename : %s\n", info.inputfilename);
    printf("Output directory : %s\n", info.output);
    
    //saves input file to output directory
    Save_Input_File(info.inputfilename,info.output);//in write_output.cpp
    

	General_Info_from_File(info,alpha1,alpha2,alphastep);
						   //Subroutine in read_input.cpp

	info.AR = info.b*info.b/info.S;  //reference aspect ratio
	
//allocates mememory for panel information in 'panelPtr'
	//for 'nopanel'-number panels
	ALLOC1D(&panelPtr,info.nopanel);

	//reads from input file panel information.
	//For each panel in panelPtr:
	//	x2[]		-x,y,z coordinates of leading edge corners
	//	c1,c2		-chord length of panel sides
	//	eps1, eps2	-incident angle of panel sides
	//	u1[], u2[]	-local free stream velocities at panel sides
	//	BC1, BC2	-boundary conditions at panel sides.
	//				 see typedef.h for more info on Boundary Conditions
	//	n 			-number of elements in chord direction
	//  left, right	-neighboring panels
	Panel_Info_from_File(panelPtr, info);	//Subroutine in read_input.cpp

//===================================================================//
		//END read general and panel info from file 'input.txt'
//===================================================================//
    
//===================================================================//
		//Start read V-tial and fuselage info from file 'input.txt'
//===================================================================//
//max of 5 VT panesl!!
//double VTchord[5],VTarea[5];	//chord and aera of VT panel
//int VTairfoil[5];
//max of 20 fuselage sections
//double FusSectS[20];		//fuselage section area
//double delFus;				//width of each section
//int FusLT;					//fuselage section where turbulent
//double IFdrag;					//interfernce drag fraction

	VT_Fus_Info\
	(info,VTchord,VTarea,VTairfoil,FusSectS,delFus,FusLT,IFdrag);
								//Subroutine in read_input.cpp
//===================================================================//
		//END read V-tial and fuselage info from file 'input.txt'
//===================================================================//

//===================================================================//
//*******************************************************************//
//===================================================================//
//
//						HORSTMANN METHOD (fixed wake)
//
//The following part is a redo of Horstmann's multiple lifting line
//method as discussed in his dissertation fron 1986.  It uses multiple
//lifting lines that have a second order spline ciruculation distribution.
//The shed wake is drag free, and consist of fixed vortex sheets with a
//linear vorticity distribution.  The method was expanded to compute the
//the induced drag along the trailing edge, as it is proposed by Eppler
//
//===================================================================//
//*******************************************************************//
//===================================================================//


//===================================================================//
		//START generation of elementary wings
//===================================================================//
	//devides panels into elementary wings and computes some basic
	//elementary wing properties.
	//input:
	// 	For each panel in panelPtr:
	//	x1[], x2[]	-x,y,z coordinates of leading edge corners
	//	c1,c2		-chord length of panel sides
	//	eps1, eps2	-incident angle of panel sides
	//	u1[], u2[]	-local free stream veloocity at panel sides
	//				 (for rotating wings}
	//	n 			-number of elementary wings in span direction
	//
	//ouput:
	//definition of elementary wing properties.
	// 	For each elementary wing in elementPtr:
	//	xo[]		-midspan location of bound vortex in global ref.frame
	//	xA[]		-control point location in global ref. frame
	//	normal[]	-surface normal at control point in global ref. frame
	//	eta			-half span of elementary wing
	//	phi 		-sweep of elementary wing
	//	nu			-dihedral of elementary wing
	//	u[3]		-local free stream velocity at element midspan
	//				 in global ref. frame, varies for rotating wings
	//	BC			-boundary conditions
	//
	// The total number of elementary wings is m*Sum(panelPtr.n)
	info.noelement=0;					//initlizes noelement
    info.nospanelement=0;                    //initlizes nospanelement
	for (i=0;i<info.nopanel;i++)			//loop over number of panels
    {
        info.noelement +=panelPtr[i].n*panelPtr[i].m;//adds no. of elements
        info.nospanelement +=panelPtr[i].n;    //adds no. of spanwise elements
    }
    info.Dsize=3*info.noelement;

//===================================================================//
		//END generation of elementary wings
//===================================================================//

	//allocating memory
	ALLOC1D(&surfacePtr,info.noelement);	//surface DVE
	ALLOC2D(&N_force,info.noelement,9);	//surface DVE normal forces
					//[0]: free stream lift, [1]: induced lift,
					//[2]: free stream side, [3]: induced side force/density
					//[4]: free str. normal, [5]: ind. normal frc/density
                    //[6,7,8]: eN_x, eN_y, eN_z in global ref. frame

    ALLOC1D(&spanPtr,info.nospanelement);//of struct STRIP
    ALLOC2D(&wakePtr,info.maxtime+1,info.nospanelement);	//wake DVE
    ALLOC1D(&CDi_DVE,info.maxtime+1);   //total induced drag (Eppler)
    
//===================================================================//
		//START wing generation
//===================================================================//

	//identifies separate wings
	Wing_Generation(panelPtr,info.nopanel,info.wing1,info.wing2,\
						info.panel1,info.panel2,info.dve1,info.dve2);
									//Subroutine in wing_geometry.cpp
//===================================================================//
		//END wing generation
//===================================================================//


//===================================================================//
		//Read in airfoil and camber data files
//===================================================================//
printf("reading in airfoils\n");

	if(info.flagCAMBER)
	{ //Skip if flagCAMBER is turned off
		// Read in the camber data. Subroutines in read_input.cpp
		Airfoil_or_Camber_Array_Size(info, &cambRow, &cambCol, 2);
		ALLOC3D(&camberPtr,cambRow,cambCol,2);
		Read_Airfoil_or_Camber(info, camberPtr,cambRow, cambCol,2);
	}
printf("Done reading in camber information\n");

	if (info.flagVISCOUS)
	{ // Skip if flagVISCOUS is turned off
		// Read in airfoil data. Subroutines in read_input.cpp
		Airfoil_or_Camber_Array_Size(info, &airfoilRow, &airfoilCol, 1);
		ALLOC3D(&airfoilPtr,airfoilRow,airfoilCol,5);
		Read_Airfoil_or_Camber(info, airfoilPtr,airfoilRow, airfoilCol,1);
	}	
// Moved reading airfoils to functions, D.F.B. 2-14-20	
printf("Done reading in aerodynamic characteristics\n");

//===================================================================//
		//DONE Reading in airfoil data files
//===================================================================//

//===================================================================//
//				Opening output files
//===================================================================//

	//Trim itereation
	//creates file "TrimSol.txt in directory "output"
    sprintf(filename,"%s%s",info.output,"TrimSol.txt");
	//open output file
	MomSol = fopen(filename, "w");

	//write header
	fprintf(MomSol,"cmac= %lf  ",info.cmac);
	fprintf(MomSol,"xcg = %lf  ycg = %lf  zcg = %lf\n",info.RefPt[0],info.RefPt[1],info.RefPt[2]);
	fprintf(MomSol,"alpha      eps      CL        CDi      ");
	fprintf(MomSol,"CMresid      CLht     CMoWing\n");

	//Performance results
	//creates file "Performance.txt in directory "output"
	sprintf(filename,"%s%s",info.output,"Performance.txt");
	//open output file
	Performance = fopen(filename, "w");

	//write header
	fprintf(Performance,"Output file of performance calculations\n");
	fprintf(Performance,"Weight    = %lf  ",info.W);
	fprintf(Performance,"Wing area = %lf  ",info.S);
	fprintf(Performance,"cmac      = %lf  \n",info.cmac);
	fprintf(Performance,"xcg = %lf  ycg = %lf  zcg = %lf\n",\
			info.RefPt[0],info.RefPt[1],info.RefPt[2]);
	fprintf(Performance,"%8s %8s %8s %8s %8s %8s %8s %8s",\
	"alpha","Vinf","CL","CD","Dtotal","L/D","wglide","Preq");
	fprintf(Performance," %8s %8s %8s %8s %8s",\
	"Dind","Dprof","Dht","Dvt","Dfus");
	fprintf(Performance," %8s %8s %8s\n",\
		"Dinter","Dmisc","CMoWing");

	fflush(Performance);


//===================================================================//
       //START setting up configuration file
//===================================================================//
    //four output files:
    //  1. configuration: holds config und summary of each flight condition
    //  2.  ----fltcfg# : holds spanwise information of each wing and case
    //  3.   -----fltcfg#: DVE summary of wings and wakes
        
    //===================================================================//
    //1. configuration file
    // file name was created up around line 148
        
    //need to check if exist - if no --> create and copy input file into it
    //                          yes --> think
    // possible options are CL trim --> append or overwrite
    // alpha sweep --> append

    //saves input file and header to configuration file in output directory
    Save_Config_Head_File(info);                    //in write_output.cpp
 //===================================================================//
        //END setting up configuration file
//===================================================================//
//===================================================================//
//               DONE Opening output files
//===================================================================//

    
//===================================================================//
//iterate CL or Alpha sweep
//===================================================================//
//===================================================================//
	//if trimCL is 1, we will trim alpha here to achieve the correct CL.
	//we will first run alpha 1 and alpha 5, and then interpolate the CL
	//result linearly and replace the first result with the new result. 
	//Then, the second result becomes the old result and we iterate this 
	//process until the condition is met. This should also work for 
	//visc ON, but it will take more iterations to converge. 

	//if trimCL is 0, we will step through the user defined alphas
    //if trimCL is 1, we will trim alpha here to achieve the correct CL.
	if (info.trimCL == 1)
    {
		q_inf = 0.5 * info.density * info.Uinf * info.Uinf;
		CLtarget = info.W / (q_inf * info.S * cos(info.bank));
		printf("\n-----TRIM FOR CL-----Target CL = %.4lf\n", CLtarget);
	//	alpha2 = 5* DtR; //starting alphas to run
	//	alpha2old = alpha2; 
	//	alpha1 = 1* DtR;
		//starting alpha based on trage CL and elliptically loaded wing - 3deg
		alpha1 = CLtarget*(info.AR+2)/info.AR*0.15923-0.0523;
		alpha2 = alpha1 + 0.07; //second alpha2 = alpha1+3
		alpha2old = alpha2; 
			alphastep = alpha2 - alpha1;
		CLtemp[0] = 0;
		CLtemp[1] = 0;
	}
	else CLtarget = 10; //if we aren't trimming for lift we need some logic here
	
    //===================================================================//
    //iterate until CL changes less than this much
    //===================================================================//
    while ((CF[2] - CLtarget) * (CF[2] - CLtarget) > delCLtarget)
    {
		//we could iterate until alpha changes by less than a certain amnt with this:
		//while (sqrt((alpha1-alpha2old) * (alpha1-alpha2old)) > 0.001*DtR)

		a2 = int((alpha2 - alpha1) / alphastep + .5);  //max. number of alpha increments

        //===================================================================//
                //looping over AOA
        //===================================================================//
		for (a = 0; a <= a2; a++)
		{
			//updating AOA info
			//the new AOA
			info.alpha = alpha1 + a * alphastep;

			printf("\nalpha = %.2lf \n", info.alpha * RtD);


			
				alphatemp1 = alpha1;
				alphatemp2 = alpha2; 
				alphasteptemp = alphastep;
				General_Info_from_File(info, alpha1, alpha2, alphastep); //reload all info from input file, since we have rotated it already in panel_rotation
				Panel_Info_from_File(panelPtr, info);
				alpha1 = alphatemp1;
				alpha2 = alphatemp2;
				alphastep = alphasteptemp;

				if (info.gradient == 0) {
					info.gradient = (9.81 * tan(info.bank)) / (info.Uinf); //this Uinf will need to include Ws too!
					if (info.gradient == 0) info.gradient = DBL_EPS;
				}

			//===================================================================//
			//START rotating panels 
			//===================================================================//
			//Rotates panels to account for sideslip, roll and alpha 
			//only applies for turning flight
			if (info.flagCIRC) {
				Panel_Rotation(info, panelPtr);
				//Subroutine in wing_geometry.cpp
			}
			//computes free stream velocity vector. THIS GETS BUILT INSIDE PANEL ROT FUNCTION if CIRC. 
			else
			{
				info.U[0] = info.Uinf * cos(info.alpha) * cos(info.beta);
				info.U[1] = info.Uinf * sin(info.beta);
				info.U[2] = info.Uinf * sin(info.alpha) * cos(info.beta) - info.Ws;
			}

			//new free stream velocities at panel edges
				//ATTENTION: NO SPANWISE VELOCITY VARIATION !!!
			for (i = 0; i < info.nopanel; i++) //THIS GETS REBUILT INSIDE CIRCLING_UINF if CIRC
			{
				panelPtr[i].u1[0] = info.U[0];
				panelPtr[i].u1[1] = info.U[1];
				panelPtr[i].u1[2] = info.U[2];

				panelPtr[i].u2[0] = info.U[0];
				panelPtr[i].u2[1] = info.U[1];
				panelPtr[i].u2[2] = info.U[2];
			}
			//===============================================================//
				//compute induced drag and lift distribution
			//===============================================================//
			LongitudinalTrim(info,panelPtr,surfacePtr,wakePtr,spanPtr, \
								CL,CY,CDi,CLi,CYi,MomSol,camberPtr,N_force);
												//Subroutine in longtrim.cpp
            //===============================================================//
            //DONE compute induced drag and lift distribution
            //===============================================================//

            //===============================================================//
            //computing free stream values
            //===============================================================//
			//q_inf = info.W / (CL * info.S);
			//V_inf = sqrt(2 * q_inf / info.density);
			V_inf = info.Uinf; //will need w added
			//induced drag
			Di = CDi * info.S * q_inf;

			//printf("CL %lf CY %lf CN %lf CDi %lf alpha %.2lf", \
				CL, CY, sqrt(CL * CL + CY * CY), CDi, info.alpha * RtD);
			printf("CL %lf CLi %lf CY %lf CYi %lf CN %lf CDi %lf",\
                   CL,CLi,CY,CYi,sqrt(CL * CL + CY * CY), CDi);
			//printf(" CN %lf CDi %lf",sqrt(CL*CL+CY*CY),CDi_DVE[timestep]);  //###
			printf("\nCFX %lf CFY %lf CFZ %lf\n", \
				CF[0], CF[1], CF[2]);
			printf("Cl %lf Cm %lf Cn %lf\n",Cl, Cm, Cn);//#
			//===============================================================//
				//computing wing/horizontal tail profile drag
			//===============================================================//

			Dprofile = 0; Dht = 0;	CMo = 0;//initialize  variables
			CDprofile = 0;

			tempS = 1 / info.nu;	//inverse of kin. viscosity

			if (info.flagVISCOUS) {
				i = 0;		//intitaliing span index counter
				m = 0;		//index of leading edge DVE
				for (k = 0; k < info.nopanel; k++)  //loop over panels
				{
					for (l = 0; l < panelPtr[k].n; l++)  //loop over span of panel k
					{
						// Calculate Reynolds number.
						if (info.flagCIRC) {//Added for circling flight: D.F.B. 03-2020
							// For circling flight calculate based on the panel left edge velocity.
							Re = norm2(surfacePtr[m].u) * 2 * surfacePtr[m].xsi * tempS * panelPtr[k].m;
						}
						else {
							// For non-circling flight, calculate based on the fixed-lift velocity
							Re = V_inf * 2 * surfacePtr[m].xsi * tempS * panelPtr[k].m;
						}

						//interpolation of airfoil drag between airfoil of panel edge 1 and 2
						//GB 2-14-20
						//computing the normal force coefficient
		 				spanPtr[i].Cd = spanPtr[i].D_force*2\
		 								/(dot(surfacePtr[m].u,surfacePtr[m].u)*spanPtr[i].area);
						spanPtr[i].Cn = sqrt(dot(spanPtr[i].Cf,spanPtr[i].Cf)\
										-spanPtr[i].Cd*spanPtr[i].Cd);
						
						airfoil = surfacePtr[m].airfoil[0];		//airfoil number, panel edge 1
						//computing section drag and moment coefficient based on panel edge 1
						cd1 = SectionDrag(airfoilPtr[airfoil], Re, spanPtr[i].Cn, airfoilCol, cm, m);
						cm1 = cm; //zero-lift moment of panel edge 1 airfoil

										//in drag_force.cpp
						airfoil = surfacePtr[m].airfoil[1];		//airfoil number, panel edge 2
						//computing section drag and moment coefficient based on panel edge 2
						cd2 = SectionDrag(airfoilPtr[airfoil], Re,spanPtr[i].Cn, airfoilCol, cm, m);
						cm2 = cm; //zero-lift moment of panel edge 2 airfoil

						cd = cd1 + surfacePtr[m].ratio * (cd2 - cd1); //weighted cd based on span location 
						cm = cm1 + surfacePtr[m].ratio * (cm2 - cm1); //weighted cmo based on span location 

						Dprofile += cd * surfacePtr[m].S * panelPtr[k].m; //wing drag

						//adding section moment coefficients (*S*chord)
						CMo += cm * surfacePtr[m].S * surfacePtr[m].xsi * 2 * panelPtr[k].m * panelPtr[k].m;

						i++;  //next span index 
						m++;	//index of next leading edge DVE 
					}
					m += panelPtr[k].n * (panelPtr[k].m - 1);  //index of next LE DVE of next panel
				}

				Dprofile *= q_inf; Dht *= q_inf;  //multiplying with dyn. pressure
				if (info.sym == 1)	//symmetrical geometry, double values
				{
					Dprofile *= 2;
					Dht *= 2;
				}

				CMo *= 2 / (info.S * info.cmac);		//non-dimensionalizing 

				//set new wing-zero lift moment to values of previous AOA
				info.CMoWing = CMo;
			}
			//===============================================================//
				//DONE computing wing/horizontal tail profile drag
			//===============================================================//

			//===============================================================//
				//START computing vertical-tail profile drag
			//===============================================================//
			Dvt = 0;	//initialize profile drag variables

			if (info.flagVISCOUS) {
				for (i = 0; i < info.noVT; i++)  //loop over surface DVEs
				{

					Re = V_inf * VTchord[i] / info.nu;	//Reynolds number

					airfoil = VTairfoil[i];		//airfoil number

					//computing the section drag coefficient, cl=0
					tempS = 0;
					cd = SectionDrag(airfoilPtr[airfoil], Re, tempS, airfoilCol, cm, i + info.noelement);
					//subroutine in drag_force.cpp

					Dvt += cd * VTarea[i];	  //h-tail drag
				}
				Dvt *= q_inf;  //multiplying with dyn. pressure
			}

			//===============================================================//
				//END computing vertical-tail profile drag
			//===============================================================//

			//===============================================================//
				//START computing fuselage drag
			//===============================================================//
			Dfuselage = 0;		//initializing

			if (info.flagVISCOUS) {
				tempS = V_inf * delFus / info.nu;	//almost local Re#

				//loop over fuselae sections
				for (i = 0; i < info.noFus; i++)
				{
					//computing local Re#
					Re = (i + 0.5) * tempS;

					if (i < FusLT)		cd = 0.664 / sqrt(Re);	//laminar flow
					else			cd = 0.0576 / pow(Re, 0.2);	//turbulent

					Dfuselage += cd * FusSectS[i];
				}
				Dfuselage *= q_inf * 1.;  //dyn. pressure and correcting pressure drag	
			}
			//===============================================================//
				//END computing fuselage drag
			//===============================================================//

			//===============================================================//
				//START total drag
			//===============================================================//
			D = Di + Dprofile + Dht + Dvt + Dfuselage + Dmisc;
			Dint = D * IFdrag;  //interference drag is a fraction of total
			D += Dint;
			//===============================================================//
				//END computing total drag
			//===============================================================//
		// *
		// *
		// *  Created by Goetz  Bramesfeld on 1/22/11.
		// *
		// *
			//===============================================================//
				//adjusting total CL for stalled sections
			//===============================================================//
			if (info.flagVISCOUS) {
				i = 0;		//intitaliing span index counter
				m = 0;		//index of leading edge DVE
				tempS = 0;        //initializing temporary CL holder
				for (k = 0; k < info.nopanel; k++)  //loop over panels
				{
					for (l = 0; l < panelPtr[k].n; l++)  //loop over span of panel k
					{
						//adding local normal force: lift/roh = cn*area*cos(dihedral)
						tempS += spanPtr[i].Cn * (surfacePtr[m].S * panelPtr[k].m)\
							* cos(surfacePtr[m].nu);

						i++;  //next span index 
						m++;	//index of next leading edge DVE 
					}
					m += panelPtr[k].n * (panelPtr[k].m - 1);  //index of next LE DVE of next panel
				}

				tempS = tempS / info.S * 2; //normalizing force/roh to overall CL
				printf("check output of new, for stall corrected CL = %lf  old was %lf", tempS, CL);
				printf("  both values should be the same (at least very similar) if no stall\n");

				CLinviscid = CL;		//inviscid CL without stall correction
				CL = tempS;        //reassigning CL

			}
			else 
			{
				CLinviscid = CL;	// Creating CLinviscid
				//When viscous corrections are off dont need to re-assign CL
			}
			CLtemp[a] = CF[2]; //storing CL value for iterating 
			//===============================================================//
				//DONE adjusting total CL for stalled sections
			//===============================================================//

			//===============================================================//
				//START write to Performance file
			//===============================================================//
			tempS = CLinviscid / CL;  //correction for stalled CL
			V_inf *= sqrt(tempS);

			CD = D / (q_inf * info.S);

			fprintf(Performance, "%8.1lf %8.2lf %8.2lf %8.5lf", \
				info.alpha * RtD, V_inf, CL, CD);
			fprintf(Performance, " %8.3lf %8.1lf %8.2lf %8.2lf", \
				D, CL / CD, V_inf * CD / CL, D * V_inf);

			fprintf(Performance, " %8.3lf %8.3lf %8.3lf %8.3lf %8.3lf", \
				Di, Dprofile, Dht, Dvt, Dfuselage);
			fprintf(Performance, " %8.3lf %8.3lf %8.3lf", \
				Dint, Dmisc, info.CMoWing);
			fprintf(Performance, " \n");
			fflush(Performance);

			//===============================================================//
				//END write to Performance file
			//===============================================================//

		//		printf(" Dvt %lf Dfus %lf Dint %lf D %lf\n",Dvt,Dfuselage,Dint,D);
            
            //===============================================================//
            //START save to config file if no CL iteration, i.e. only alpha swee
            //===============================================================//
           if (info.trimCL == 0)  //simple alfa sweep -> safe results of flight
            {
                FltConfg = fopen(info.config,"a"); //append file
                
                if (FltConfg == NULL)
                {
                    fclose(FltConfg);
                    printf(" couldn't save to FltConfg file\nPress any key to exit...\n");
                    exit(EXIT_FAILURE);
                }
                fprintf(FltConfg,"%-10d%-10.5lf%-10.3lf%-10.3lf",\
                        a+1,CLtarget,info.alphaset*RtD,info.betaset*RtD);
                fprintf(FltConfg,"%-12lf%-12lf%-12lf%-12lf%-12lf%-12lf%-12lf%-12lf%-12lf",\
                        CL,CY,CDi,CF[0],CF[1],CF[2],Cl,Cm,Cn);
                fprintf(FltConfg,"\n");
                
                fclose(FltConfg);

    	        //save information of spanwise strips
	        	SaveSpanDVEInfo(panelPtr,surfacePtr,wakePtr,spanPtr,N_force,a+1,info.timestep);
        											//in write_output.cpp
            }
            //===============================================================//
            //END save to config file,  no CL iteration
            //===============================================================//

        }//end loop over 'a' angle of attack
        //===================================================================//
        //            END OF LOOP OVER AOA
        //===================================================================//

		if (info.trimCL == 1)
        {
			//calculate new alpha to run for CL iteration
			alpha2 = alpha1 + (CLtarget - CLtemp[0]) / ((CLtemp[1] - CLtemp[0]) / (alpha2old - alpha1));
			alpha2old = alpha1;
			CLtemp[1] = CLtemp[0];

			alpha1 = alpha2;
			//alphastep = alpha2 - alpha1;
		}
		else CLtarget = CF[2]; //set here so we dont iterate again if we aren't trimming
	}
    //===================================================================//
    //END iterate CL
    //===================================================================//
 
    //===============================================================//
    //START save to config file if CL iteration
    //===============================================================//
    if (info.trimCL == 1)  //CL iteration -> safe results of flight config
    {
       int FCno=0; //Flight condition number
       double CLoldtarget,alphaold; //CL and alpha of last rund

        FltConfg = fopen(info.config,"r+"); //append file
        
        if (FltConfg == NULL)  //check if flight config file can be appended
        {
            fclose(FltConfg);
            printf(" couldn't save to FltConfg file\nPress any key to exit...\n");
            exit(EXIT_FAILURE);
        }
        
        //check if there already exist an entry. if not string[5]='-' and FCno set to 0
        fseek(FltConfg, -Linelength-1,SEEK_END);  //move pointer at the beginning of the last in file
  		if(string[5]=='-')  FCno=0; 
 		else
 		{
 			fseek(FltConfg, -Linelength-1,SEEK_END);  //move pointer at the beginning of the last line
     		fscanf(FltConfg,"%d %lf%lf",&FCno,&CLoldtarget,&alphaold);
     	}

         //three cases:
            //1. no entry yet -> append
            //2. same CLtarget -> update last line
            //3. new CLtarget -> append

        if(FCno == 0) //no entry yet
        {
            FCno=1; 
            printf("append empty\n");
            fseek(FltConfg,0,SEEK_END);  //move pointer at the beginning of the last line
                       fprintf(FltConfg,"%-10d%-10.5lf%-10.3lf%-10.3lf",\
                FCno,CLtarget,info.alphaset*RtD,info.betaset*RtD);
            fprintf(FltConfg,"%-12lf%-12lf%-12lf%-12lf%-12lf%-12lf%-12lf%-12lf%-12lf",\
                CL,CY,CDi,CF[0],CF[1],CF[2],Cl,Cm,Cn);
            fprintf(FltConfg,"\n");
        }
        else if((CLtarget-CLoldtarget)*(CLtarget-CLoldtarget) <= delCLtarget*1.01) //update CLtarget 
        {													//if within 1% of convergence criterion
            fseek(FltConfg,0,SEEK_END);  //move pointer at the beginning of the last line
            fseek(FltConfg, -Linelength,SEEK_END);  //move pointer at the beginning of the last line
            printf("CLtarget overwrite\n");
            fprintf(FltConfg,"%-10d%-10.5lf%-10.3lf%-10.3lf",\
                FCno,CLtarget,info.alpha*RtD,info.beta*RtD);
            fprintf(FltConfg,"%-12lf%-12lf%-12lf%-12lf%-12lf%-12lf%-12lf%-12lf%-12lf",\
                CL,CY,CDi,CF[0],CF[1],CF[2],Cl,Cm,Cn);
            fprintf(FltConfg,"\n");
        }
        else //add next target CL
        {
            FCno++; //updating flight condition counter
            fprintf(FltConfg,"%-10d%-10.5lf%-10.3lf%-10.3lf",\
                FCno,CLtarget,info.alpha*RtD,info.beta*RtD);
            fprintf(FltConfg,"%-12lf%-12lf%-12lf%-12lf%-12lf%-12lf%-12lf%-12lf%-12lf",\
                CL,CY,CDi,CF[0],CF[1],CF[2],Cl,Cm,Cn);
            fprintf(FltConfg,"\n");
         printf("New flight configuration added\n");
        }

        //update number of flight condition number in header 
        fseek(FltConfg,0,SEEK_SET);  //move pointer to the beginning of the file
        //find mmarker
        pos = 0; //initializing position counter
        do fscanf(FltConfg,"%s",string);
        while(strcmp(string,"#$&#$&")!=0);

    	//find the '='-sign in input file before "projected ref area"
//		for(n=0;n<3;n++) do	answer=fgetc(FltConfg); while (answer!='=');
//could update with projected area
		//find '=' before Unraveled area
//		do	answer=fgetc(FltConfg); while (answer!='=');
//could update with unraveled area
	
		//find the '='-sign in input file before "number of flight conditions"
		for(n=0;n<8;n++) do	answer=fgetc(FltConfg); while (answer!='=');


 		pos = ftell(FltConfg);
 		fseek(FltConfg,pos,SEEK_SET);  //move pointer to the current position in file
		fprintf(FltConfg," %d",FCno);

        fclose(FltConfg);

        //save information of spanwise strips
        SaveSpanDVEInfo(panelPtr,surfacePtr,wakePtr,spanPtr,N_force,FCno,info.timestep);
        										//in write_output.cpp
   }
    //===============================================================//
    //END save to config file if CL iteration
    //===============================================================//

    //===============================================================//
		// Saving info for VoGen
	//===============================================================//
 	sprintf(filename, "%sVoGen.txt", info.output);
    FltConfg = fopen(filename,"w"); //create file

// 	sprintf(filename, "VoGen.txt" );
//  FltConfg = fopen(filename,"a"); //create file

    if (FltConfg == NULL)  //check if flight config file can be appended
    {
        fclose(FltConfg);
        printf(" couldn't open VoGen summary file\nPress any key to exit...\n");
        exit(EXIT_FAILURE);
    }

    //header of tabulated case data
    fprintf(FltConfg,"%-10s%-10s%-10s","CL","Alpha","Beta");
    fprintf(FltConfg,"%-12s%-12s%-12s","CL","CQ","CDI");
    fprintf(FltConfg,"%-12s%-12s%-12s","CX","CY","CZ");
    fprintf(FltConfg,"%-12s%-12s%-12s\n","CL","CM","CN");
    fprintf(FltConfg,"%-10s%-10s%-10s","target","[deg]","[deg]");
    fprintf(FltConfg,"%-12s%-12s%-12s","(lift)","(side)"," ");
    fprintf(FltConfg,"%-12s%-12s%-12s"," "," "," ");
    fprintf(FltConfg,"%-12s%-12s%-12s\n","(roll)"," "," ");


    fprintf(FltConfg,"%-10.5lf%-10.3lf%-10.3lf",\
        CLtarget,info.alpha*RtD,info.beta*RtD);
    fprintf(FltConfg,"%-12lf%-12lf%-12lf%-12lf%-12lf%-12lf%-12lf%-12lf%-12lf",\
        CL,CY,CDi,CF[0],CF[1],CF[2],Cl,Cm,Cn);
    fprintf(FltConfg,"\n");
    fclose(FltConfg);
	//===============================================================//
		// DONESaving info for VoGen
	//===============================================================//

	//free allocated memory
	FREE1D(&panelPtr,info.nopanel);
	FREE1D(&surfacePtr,info.noelement);
	FREE2D(&N_force,info.noelement,9);	
    FREE1D(&spanPtr,info.nospanelement);
    FREE2D(&wakePtr,info.maxtime+1,info.nospanelement);
    FREE1D(&CDi_DVE,info.maxtime+1);
	if(info.flagCAMBER) 
		FREE3D(&camberPtr,cambRow,cambCol,2);
	if (info.flagVISCOUS)
		FREE3D(&airfoilPtr,airfoilRow,airfoilCol,5);

	fclose(MomSol);//close output file of trim iteration results
	fclose(Performance);//close output file of performance calc's
//printf("done\n");
//printf("push any key and return ",PROGRAM_VERSION);
//scanf("%c",&answer);

	printf("Done FreeWake\n");
	printf("===========================================================================\n");
	printf("===========================================================================\n\n");
	return(0);
}
//===================================================================//
		//END of program
//===================================================================//
