#define _CRT_SECURE_NO_WARNINGS
#include "general.h"
#include "PerfCode.h"

int i,ii,a,a2;		//loop counters, max AOA increment
	int k,l,m;			//loop counters
	int HTpanel;		//index of first HT panel
	
	//drag variables
	double *cn,cd;		//section normal and drag force coefficients
	double cm,CMo;			//section, total zero-lift moment coefficient
	double CL,CD,CDi=0;		//aircraft lift, total and induced drag coefficients
	double CLinviscid=0;	//inviscid CL without stall correction
	double CDprofile=0;		//wing profile drag coefficient
	double CDfuselage=0;	//fuselage drag coefficeing
	double CDht=0,CDvt=0;	//horizontal and vertical tail drag coefficeing
	//coefficients with cap C are wrt to wing area
	
	double D=0, Di=0;		//total and induced drag
	double Dprofile=0;		//wing profile 
	double Dfuselage=0;		//fuselage drags
	double Dvt=0,Dht=0;		//vertical tail,horizontal tial drags
	double Dint=0,Dmisc=0;	//interference, miscillanuous drags
	
	double V_inf=0,q_inf=0; //free stream velocity and dynamic pressure
	double Re;			//Reynolds number
	
	// Airfoil related variables
	double profile[20][2000][5];				//array holding airfoil data 
		//[airfoil#][row in input file][0=alfa,1=cl,2=cd,3=Re,4=cm]
		//array is sorted Re# ascending, CL of each Re# ascending
		//max no. of airfoils = 20; max no. of points per airfoil = 2000 numbers increased G.b. 5/16/11
	int rows[20];				//numbers of rows of airfoil file numbers increased G.b. 5/16/11
	int airfoil=0;			//airfoil index

	//V-tail and fuselage	
	//max of 5 VT panesl!!
	double VTchord[5],VTarea[5];	//chord and aera of VT panel
	int VTairfoil[5];
	//max of 20 fuselage sections
	double FusSectS[20];		//fuselage section area
	double delFus;				//width of each section
	int FusLT;					//fuselage section where turbulent
	double IFdrag;				//interfernce drag fraction


	double tempS;
	char answer ;
	char filename[137];	//file path and name

	//Input/output files
	FILE *AD;			//airfoil data file
	FILE *MomSol;		//output file of trim solutions
	FILE *Performance;	//output file for performance results




int main()

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

printf("===========================================================================\n ");
printf("\n\tPerformance Code based on FreeWake 2015 (includes stall model)\n\n");
printf("\t\tRunning Version %s \n ",PROGRAM_VERSION);
printf("===========================================================================\n ");

	

	//Input/output files
	FILE *AD;			//airfoil data file
	FILE *MomSol;		//output file of trim solutions
	FILE *Performance;	//output file for performance results

//printf("Do you want to start %s? y for YES ",PROGRAM_VERSION);
//scanf("%c",&answer);
//if(answer != 'y' && answer != 'Y')		exit(0);
//Delete_timestep();

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

	General_Info_from_File(info,alpha1,alpha2,alphastep);
						   //Subroutine in read_input.cpp

	info.AR = info.b*info.b/info.S;  //reference aspect ratio

	//allocates mememory for panel information in 'panelPtr'
	//for 'nopanel'-number panels
	ALLOC1D(&panelPtr,info.nopanel);


//printf("\n\n\t\t!!!Only ONE CHORDWISE ROW OF SURFACE DVES!!\n");
//printf("\t\t\t\tThat means m=1!!\n");
//printf("\n\n\t\t\t!!!ONLY TWO WINGSS!!\n");
//printf("\t\tThat means one main wing and one horizontal tail!!\n");
//printf("\t\t\t\tG.B. 8-8-2007\n");
//info.m=1;  //no. of chordwise surface DVEs fixed to one, G.B. 3-5-2007


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
	info.nospanelement=0;					//initlizes noelement
	for (i=0;i<info.nopanel;i++)			//loop over number of panels
		info.nospanelement +=panelPtr[i].n;	//adds no. of spanwise elements
	info.noelement=info.nospanelement*info.m;//the total number of el. wings
	info.Dsize=3*info.noelement;

//===================================================================//
		//END generation of elementary wings
//===================================================================//

//allocating memory
ALLOC1D(&surfacePtr,info.noelement);	//surface DVE
ALLOC1D(&cn,info.nospanelement);	//normal force coeff. of wing section

//===================================================================//
		//START wing generation
//===================================================================//

		//identifies separate wings
		Wing_Generation(panelPtr,info.nopanel,info.wing1,info.wing2,\
						info.panel1,info.panel2);
									//Subroutine in wing_geometry.cpp

//===================================================================//
		//END wing generation
//===================================================================//

//===================================================================//
		//START determining panel index which belong to HT
//===================================================================//

	i = info.wing1[1]/info.m;  //number of elements in spanwise direction

	HTpanel = 0; //starting with the first panel (index 0)

	//subtract no. of spanwise elements until the first HT panel is reached
	do{
		i-=panelPtr[HTpanel].n;
		HTpanel ++;
	}
	while(i>0);

//===================================================================//
		//DONE determining panel index which belong to HT
//===================================================================//

//===================================================================//
		//Read in airfoil data files
//===================================================================//
//	double profile[15][500][5];				//array holding airfoil data 
		//[airfoil#][row in input file][0=alfa,1=cl,2=cd,3=Re,4=cm]
		//array is sorted Re# ascending, CL of each Re# ascending
		//max no. of airfoils = 15; max no. of points per airfoil = 500
//	int rows[15];				//numbers of rows of airfoil file
//	int airfoil=0;			//airfoil index
//	char ch;
	
	//initializing profile
	for(airfoil=0;airfoil<info.noairfoils;airfoil++)
	for(i=0;i<2000;i++)
	for(ii=0;ii<5;ii++)
			profile[airfoil][i][ii]=0;

	//read airfoil data	
	for(airfoil=0;airfoil<info.noairfoils;airfoil++)
	{
		//creates file name airfoil##.dat ## is the number of the airfoil
		sprintf(filename,"%s%s%d%s",AIRFOIL_PATH,"airfoil",airfoil+1,".dat");

		// checks if airfoil file exists
		if ((AD = fopen(filename, "r"))== NULL) {
			printf("Airfoil file %d could not be opened:\n",airfoil+1);
			exit(1);
		}

		//opens airfoil file
		AD = fopen(filename, "r");
		
		//read in number of rows
		do	
		answer = fgetc(AD);
		while (answer!='=');
		//reads relaxed wake flag
		fscanf(AD,"%d", &rows[airfoil]);
		if(rows[airfoil]>2000)
		{
			printf("\n\t\t\t!!!\n");
			printf("Airfoil %d has more than the maximum ",airfoil+1);
			printf(" allowable number of rows (2000 < %d)\n",rows[airfoil]);
			printf("push any key and return ",PROGRAM_VERSION);
			scanf("%c",&answer);
			scanf("%c",&answer);
			exit(0);
		}
		
		//read in airfoil data, row by row
		for(i=0;i<rows[airfoil];i++)
		{
			fscanf(AD,"%lf", &profile[airfoil][i][0]); //angle of attack
				profile[airfoil][i][0] *=DtR;	//conversion to radians
			fscanf(AD,"%lf", &profile[airfoil][i][1]); //cl
			fscanf(AD,"%lf", &profile[airfoil][i][2]); //cd
			fscanf(AD,"%lf", &profile[airfoil][i][3]); //Re
			fscanf(AD,"%lf", &profile[airfoil][i][4]); //cm
		}
		fclose(AD);
	}

//===================================================================//
		//DONE Reading in airfoil data files
//===================================================================//

//===================================================================//
//				Opening output files
//===================================================================//

	//Trim itereation
	//creates file "TrimSol.txt in directory "output"
	sprintf(filename,"%s%s",OUTPUT_PATH,"TrimSol.txt");
	//open output file
	MomSol = fopen(filename, "w");

	//write header
	fprintf(MomSol,"cmac= %lf  ",info.cmac);
	fprintf(MomSol,"xcg = %lf  ycg = %lf  zcg = %lf\n",info.RefPt[0],info.RefPt[1],info.RefPt[2]);
	fprintf(MomSol,"alpha      eps      CL        CDi      ");
	fprintf(MomSol,"CMresid      CLht     CMoWing\n");

	//Performance results
	//creates file "Performance.txt in directory "output"
	sprintf(filename,"%s%s",OUTPUT_PATH,"Performance.txt");
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
//===================================================================//
		//looping over AOA
//===================================================================//
//===================================================================//

	a2 = int((alpha2-alpha1)/alphastep+.5);  //max. number of alpha increments

	for(a=0;a<=a2;a++)
	{
		//updating AOA info
		//the new AOA
		info.alpha = alpha1+a*alphastep;

		printf("\nalpha = %.2lf \n",info.alpha*RtD);

		//computes free stream velocity vector
		info.U[0]=info.Uinf*cos(info.alpha)*cos(info.beta);
		info.U[2]=info.Uinf*sin(info.alpha)*cos(info.beta);

		//new free stream velocities at panel edges
			//ATTENTION: NO SPANWISE VELOCITY VARIATION !!!
		for(i=0;i<info.nopanel;i++)
		{
			panelPtr[i].u1[0]= info.U[0];
			panelPtr[i].u1[2]= info.U[2];

			panelPtr[i].u2[0]= info.U[0];
			panelPtr[i].u2[2]= info.U[2];
		}

	//===============================================================//
		//compute induced drag and lift distribution
	//===============================================================//
		LongitudinalTrim(info,panelPtr,surfacePtr,HTpanel,cn,\
						 CL,CDi,MomSol); //Subroutine in longtrim.cpp
	//===============================================================//
		//DONE compute induced drag and lift distribution
	//===============================================================//

	//===============================================================//
		//computing free stream values
	//===============================================================//
		q_inf = info.W/(CL*info.S);
		V_inf = sqrt(2*q_inf/info.density);

		//induced drag
		Di = CDi*info.S*q_inf;
printf("CL %lf V %lf CDi %lf \n",CL,V_inf,CDi);
	//===============================================================//
		//computing wing/horizontal tail profile drag
	//===============================================================//
		Dprofile=0; Dht=0;	CMo=0;//initialize  variables
		CDprofile=0;

		tempS = 1/info.nu;	//inverse of kin. viscosity

		i=0;		//intitaliing span index counter
		m=0;		//index of leading edge DVE
		for(k=0;k<info.nopanel;k++)  //loop over panels
		{
		  for(l=0;l<panelPtr[k].n;l++)  //loop over span of panel k
		  {
			Re = V_inf*2*surfacePtr[m].xsi*tempS*info.m;	//Reynolds number
			airfoil = surfacePtr[m].airfoil;		//airfoil number

			//computing section drag and moment coefficient
			cd=SectionDrag(profile[airfoil],Re,cn[i],rows[airfoil],cm,m);
							//in drag_force.cpp
			Dprofile += cd*surfacePtr[m].S*info.m; //wing drag

			//adding section moment coefficients (*S*chord)
			CMo += cm*surfacePtr[m].S*surfacePtr[m].xsi*2*info.m*info.m;
			
			i++;  //next span index 
			m++;	//index of next leading edge DVE 
		  }
		  m += panelPtr[k].n*(info.m-1);  //index of next LE DVE of next panel
		}
	
		Dprofile*=q_inf; Dht*=q_inf;  //multiplying with dyn. pressure
		if(info.sym == 1)	//symmetrical geometry, double values
		{
			Dprofile *=2;
			Dht *= 2;
		}	
		
		CMo*=2/(info.S*info.cmac);		//non-dimensionalizing 

		//set new wing-zero lift moment to values of previous AOA
		info.CMoWing = CMo;

	//===============================================================//
		//DONE computing wing/horizontal tail profile drag
	//===============================================================//

	//===============================================================//
		//START computing vertical-tail profile drag
	//===============================================================//
		Dvt=0;	//initialize profile drag variables

		for(i=0;i<info.noVT;i++)  //loop over surface DVEs
		{
			
			Re = V_inf*VTchord[i]/info.nu;	//Reynolds number

			airfoil = VTairfoil[i];		//airfoil number
 
			//computing the section drag coefficient, cl=0
			tempS =0;
			cd=SectionDrag(profile[airfoil],Re,tempS,rows[airfoil],cm,i+info.noelement);
				//subroutine in drag_force.cpp
		
			Dvt += cd*VTarea[i];	  //h-tail drag
		}
		Dvt*=q_inf;  //multiplying with dyn. pressure

	//===============================================================//
		//END computing vertical-tail profile drag
	//===============================================================//

	//===============================================================//
		//START computing fuselage drag
	//===============================================================//
		Dfuselage=0;		//initializing

		tempS = V_inf*delFus/info.nu;	//almost local Re#

		//loop over fuselae sections
		for(i=0;i<info.noFus;i++)
		{
			//computing local Re#
			Re=(i+0.5)*tempS;

			if(i<FusLT)		cd=0.664/sqrt(Re);	//laminar flow
			else			cd=0.0576/pow(Re,0.2);	//turbulent
			
			Dfuselage += cd*FusSectS[i];
		}
		Dfuselage *=q_inf*1.;  //dyn. pressure and correction for pressure drag	
		
	//===============================================================//
		//END computing fuselage drag
	//===============================================================//

	//===============================================================//
		//START total drag
	//===============================================================//
		D = Di+Dprofile+Dht+Dvt+Dfuselage+Dmisc;
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

		i=0;		//intitaliing span index counter
		m=0;		//index of leading edge DVE
        tempS=0;        //initializing temporary CL holder
		for(k=0;k<info.nopanel;k++)  //loop over panels
		{
		  for(l=0;l<panelPtr[k].n;l++)  //loop over span of panel k
		  {
            //adding local normal force: lift/roh = cn*area*cos(dihedral)
            tempS += cn[i]*(surfacePtr[m].S*info.m)\
                    *cos(surfacePtr[m].nu);

 			i++;  //next span index 
			m++;	//index of next leading edge DVE 
		  }
		  m += panelPtr[k].n*(info.m-1);  //index of next LE DVE of next panel
		}

        tempS = tempS/info.S*2; //normalizing force/roh to overall CL
        printf("check output of new, for stall corrected CL = %lf  old was %lf",tempS, CL);
        printf("  both values should be the same (at least very similar) if no stall\n");
        
        CLinviscid=CL;		//inviscid CL without stall correction
		CL=tempS;        //reassigning CL
        
    //===============================================================//
        //DONE adjusting total CL for stalled sections
    //===============================================================//
    
	//===============================================================//
		//START write to Performance file
	//===============================================================//
		tempS = CLinviscid/CL;  //correction for stalled CL
		V_inf *= sqrt(tempS);
		
		CD = D/(q_inf*info.S);

		fprintf(Performance,"%8.1lf %8.2lf %8.2llf %8.5lf",\
					info.alpha*RtD,     V_inf,  CL,CD);
		fprintf(Performance," %8.3lf %8.1lf %8.2lf %8.2lf",\
							D,    CL/CD,  V_inf*CD/CL, D*V_inf);

		fprintf(Performance," %8.3lf %8.3lf %8.3lf %8.3lf %8.3lf",\
		Di,Dprofile,Dht,Dvt,Dfuselage);
		fprintf(Performance," %8.3lf %8.3lf %8.3lf",\
		Dint,Dmisc,info.CMoWing);
		fprintf(Performance," \n");
		fflush(Performance);

	//===============================================================//
		//END write to Performance file
	//===============================================================//

printf(" Dvt %lf Dfus %lf Dint %lf D %lf\n",Dvt,Dfuselage,Dint,D);


	}//end loop over 'a' angle of attack
//===================================================================//
//			END OF LOOP OVER AOA
//===================================================================//

	//free allocated memory
	FREE1D(&panelPtr,info.nopanel);
	FREE1D(&surfacePtr,info.noelement);
	FREE1D(&cn,info.noelement);
	
	fclose(MomSol);//close output file of trim iteration results
	fclose(Performance);//close output file of performance calc's
//printf("done\n");
//printf("push any key and return ",PROGRAM_VERSION);
//scanf("%c",&answer);
//scanf("%c",&answer);
	//return(0);
}
//===================================================================//
		//END of program
//===================================================================//
