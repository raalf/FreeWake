void LongitudinalTrim(GENERAL &,PANEL *,DVE *,DVE **,STRIP *&,double &,\
						double &,double &, double&, double&, FILE *,double ***,double **&);

void LongitudinalTrim(GENERAL &info,PANEL *panelPtr,DVE *surfaceDVEPtr,DVE **wakePtr,STRIP *&spanPtr\
						,double &CL,double &CY,double &CDi, double &CLi, double &CYi,\
						FILE *MomSol,double ***camberPtr,double **&N_force)
{
// This program finds longitudinal trim solutions for twin wing configurations
//
//  Wing 1 is the main wing, Wing 2 is the horizontal tail.
//
//ATTENTION: NO SPANWISE VELOCITY VARIATION !!!
//If you don't know what this means, don't worry and move on.
//
// INPUT
//	info		general information
//	panelPtr	panel information
//
//output
//	surfaceDVEPtr	geometry of surface elements
// 	spanPtr		strip (section) information
//	Cn			normal force coefficient of each span element
//	Cd			drag force coefficient of each span element
// 	CL			total lift
//  CY          total side force
//	CDi			total induced drag
// 	N_force;	surface DVE's normal forces/density
				//[0]: free stream lift, [1]: induced lift,
				//[2]: free stream side, [3]: induced side force/density
				//[4]: free str. normal, [5]: ind. normal frc/density
                //[6,7,8]: eN_x, eN_y, eN_z in global ref. frame
//
//the routine also saves the trim information to TrimSol.txt
//furthermore saved is the spanwise lift distribution and spanwise
//drag distribution in SpanAOA<info.alpha>.txt where <info.alpha>
//holds the numerical value of the current angle of attack 

	int i,j,k,l,m,n=0,p=0;				//loop counter
	double epsilonHT;	//HT incident correction [rad]
	double epsilonAIL; //aileron incident correction [rad]
	double CM_resid=0;	//residual pitch moment
	double Cl_resid=0;	//residual roll moment
	double *mult; //aileron diferential storage
	double eps1;		//the HT incidence angle of one iteration previous
	double eps1ail;		//the aileron incidence angle of one iteration previous
	double CM_old;		//saving pitching moment of previous iteration step
	double Cl_old;		//saving roll moment of previous iteration step
	double CLht,CLhti;	//total and induced lift of HT
	double *cl,*cy;		//section lift and side force coefficients
	double tempS, tempsum, tempR;
    int tempI;

	FILE *spaninfo;			//output file for spanwise information
	char filename[133];	//file path and name for spanwise information
	char answer;

	//allocating memory	
	ALLOC1D(&cl,info.nospanelement);	//section lift coefficient
   	ALLOC1D(&cy,info.nospanelement);	//section side force coefficient
 	ALLOC1D(&mult, 2);	//will store diff aileron info

	//initial HT incident correction
	epsilonHT = 0;//panelPtr[i].eps1;  	//[rad]
	eps1 = 0;	//the HT incidence angle of one iteration previous

	//initial ail incident correction
	epsilonAIL = 0;//panelPtr[i].eps1;  	//[rad]
	eps1ail = 0;	//the HT incidence angle of one iteration previous
	//===========================================//
		//pitching moment computation - Step 0
	//===========================================//

    //computed the residual moment coefficient of wing and tail
	PitchingMoment(info,panelPtr,surfacePtr,wakePtr,info.cmac,\
				info.RefPt,CLht,CLhti,\
				N_force,spanPtr,CL,CY,CDi,camberPtr,CLi,CYi);
                                    //Subroutine in PitchMoment.cpp


	//adding zero lift of wing only if camber is turned off
	if(~info.flagCAMBER){Cm += info.CMoWing;}

	CM_resid = Cm; //Cm and Cl are global
	Cl_resid = Cl;

	if(info.trimPITCH==1 || info.trimROLL==1)   // trim routine for pitch and roll
	{
		if (info.trimPITCH == 1) {
			printf("\n-----TRIM FOR PITCH Cm-----Target Cm = 0.0\n");

			//check if there are any elevators defined on wing 2:
			tempsum = 0;
			for (i = info.panel1[1]; i <= info.panel2[1]; i++) { //for all panels on wing 2
				tempsum = +panelPtr[i].deflect1;
			}
			if (tempsum == 0) {
				printf("No elevator defined on wing 2.\n You must define a deflection on at least one panel on wing 2. \n");
				scanf("%c", &answer);
				exit(1);
			}
		
			if (CM_resid > 0) {

				epsilonHT = 0.035; //deflect HT trailing edge 2deg down
			}
			else {
				epsilonHT = -0.035; //deflect HT trailing edge 2deg up
			}


			for (i = info.panel1[1]; i <= info.panel2[1]; i++) { //for all panels on wing 2
				if (panelPtr[i].deflect1 != 0){ //if the user has defined a deflection for this panel, we will overwrite it with our trim epsilon.				
					panelPtr[i].deflect1 = epsilonHT;
					panelPtr[i].deflect2 = epsilonHT;
				}
			}
		}

		if (info.trimROLL == 1) {
			printf("\n-----TRIM FOR ROLL Cl-----Target Cl = 0.0\n");
			//check if there are any ailerons defined on wing 2:
			tempsum = 0;
			for (i = info.panel1[0]; i <= info.panel2[0]; i++) { //for all panels on wing 1
				tempsum = +panelPtr[i].deflect1;
			}
			if (tempsum == 0) {
				printf("No aileron defined on wing 1.\n You must define a deflection on at least one panel on wing 1. \n");
				scanf("%c", &answer);
				exit(1);
			}

			if (Cl_resid > 0) {

				epsilonAIL = 0.035; //deflect trailing edge 2deg down
			}
			else {
				epsilonAIL = -0.035; //deflect trailing edge 2deg up
			}

			j = 0;
			for (i = info.panel1[0]; i <= info.panel2[0]; i++) { //for all panels on wing 1
				if (panelPtr[i].deflect1 != 0) { 
					mult[j] = sqrt(panelPtr[i].deflect1* panelPtr[i].deflect1) *RtD;//use the input ail deflection as the multiplier for differential ailerons
					j++; 
				}
			}

			j = 0;
			for (i = info.panel1[0]; i <= info.panel2[0]; i++) { //for all panels on wing 1
				if (panelPtr[i].deflect1 != 0) { //if the user has defined a deflection for this panel, we will overwrite it with our trim aileron position. 
					if (j == 0) { //enter deflection for left wing
						panelPtr[i].deflect1 = epsilonAIL;
						panelPtr[i].deflect2 = epsilonAIL;
					}
					else if (j == 1) { //reverse if this is the right wing 
						panelPtr[i].deflect1 = -epsilonAIL;
						panelPtr[i].deflect2 = -epsilonAIL;
					}
					else {
						printf("Maximum of two ailerons.\n");
						scanf("%c", &answer);
						exit(1);
					}

					if (panelPtr[i].deflect1 > 0) { //use first multiplier if deflecting down 
						panelPtr[i].deflect1 *= mult[0];
						panelPtr[i].deflect2 *= mult[0];
					}
					else { //use second multiplier if deflecting up
						panelPtr[i].deflect1 *= mult[1];
						panelPtr[i].deflect2 *= mult[1];
					}
					j++;
				}
			}

			eps1ail = 1*DtR; //first run from earlier had a aileron deflection that was entered in the input
		}
	//========================================//
	//	pitching moment iteration loop
	//========================================//
	  p=0;	//initializing loop counter
	  do
	  {
	//========================================//
	//	pitching moment computation - Step i
	//========================================//
		  if (info.trimPITCH == 1) {
			  printf("\nelevator = %.2lf deg ", epsilonHT * RtD);
		  }
		  if (info.trimROLL == 1) {
			  printf("\naileron = %.8lf deg * multiplier ", epsilonAIL * RtD);
		  }
		//saving pitching moment of previous iteration step
		CM_old = CM_resid;
		Cl_old = Cl_resid; 

		//computed the residual moment coefficient of wing and tail
		PitchingMoment(info,panelPtr,surfacePtr,wakePtr,info.cmac,\
					info.RefPt,CLht,CLhti,\
					N_force,spanPtr,CL,CY,CDi,camberPtr,CLi,CYi);
                                    //Subroutine in PitchMoment.cpp

		//adding zero lift of wing only if camber is turned off
		if (~info.flagCAMBER) { Cm += info.CMoWing; } //changed to account for camber, BBB Apr 2020

		CM_resid = Cm; //Cm and Cl are global
		Cl_resid = Cl;

		if (info.trimPITCH == 1) {
			//computing new HT incidence angle
			tempS = epsilonHT - CM_resid * (epsilonHT - eps1) / (CM_resid - CM_old);
			
			eps1 = epsilonHT;  epsilonHT = tempS;	//reassigning HT angles
			for (i = info.panel1[1]; i <= info.panel2[1]; i++) { //for all panels on wing 2
				if (panelPtr[i].deflect1 != 0) { //if the user has defined a deflection for this panel, we will overwrite it with our trim epsilon. 
					panelPtr[i].deflect1 = epsilonHT;
					panelPtr[i].deflect2 = epsilonHT;
				}
			}
		}
		else {
			CM_resid = 0; //set residual pitch moment to 0 if we aren't trimming pitch
		}
		if (info.trimROLL == 1) {
			//computing new aileron incidence angle
			tempR = epsilonAIL - Cl_resid * (epsilonAIL - eps1ail) / (Cl_resid - Cl_old);
			
			eps1ail = epsilonAIL;  epsilonAIL = tempR;	//reassigning aileron angles
			j = 0;
			for (i = info.panel1[0]; i <= info.panel2[0]; i++) { //for all panels on wing 1
				if (panelPtr[i].deflect1 != 0) { //if the user has defined a deflection for this panel, we will overwrite it with our trim aileron position. 
					if (j == 0) { //enter deflection for left wing
						panelPtr[i].deflect1 = epsilonAIL;
						panelPtr[i].deflect2 = epsilonAIL;
					}
					else if (j == 1) { //reverse if this is the right wing
						panelPtr[i].deflect1 = -epsilonAIL;
						panelPtr[i].deflect2 = -epsilonAIL;
					}
					else {
						printf("Maximum of two ailerons.\n");
						scanf("%c", &answer);
						exit(1);
					}

					if (panelPtr[i].deflect1 > 0) { //use first multiplier if deflecting down 
						panelPtr[i].deflect1 *= mult[0];
						panelPtr[i].deflect2 *= mult[0];
					}
					else { //use second multiplier if deflecting up 
						panelPtr[i].deflect1 *= mult[1];
						panelPtr[i].deflect2 *= mult[1];
					}
					j++;
				}
			}
		}
		else {
			Cl_resid = 0; //set residual roll moment to 0 if we aren't trimming roll
		}
		

		p++;	//incrementing loop counter

	  }while((p<20) && ((CM_resid*CM_resid>.000025) || (Cl_resid * Cl_resid > .000025)));
	//========================================//
	//	pitching moment iteration loop  END
	//========================================*///
	}   //end of longitudinal trim routine

//===================================================================//
//		Compute cn and Cd (spanwise values)
//===================================================================//
/* Removed because different output strucuture

	i=0;		//intitaliing span index counter
	m=0;		//index of leading edge DVE
	for(k=0;k<info.nopanel;k++)  //loop over panels
	{
		for(l=0;l<panelPtr[k].n;l++)  //loop over span of panel k
		{
			//working each spanwise strip
			//initializing
			spanPtr[i].Cn=0;
			cl[i]=0;
			cy[i]=0;
//			S[i] = 0;  //area of current spanwise strip
			
			//adding up chordal values of one spanwise location (indexed i)
			for(m=0;m<panelPtr[k].m;m++)
			{
				j=n+m*panelPtr[k].n;  //counting index along chord
//
                spanPtr[i].Cn += (N_force[j][4]+N_force[j][5]); //adding forces/density along chord
                cl[i] += (N_force[j][0]+N_force[j][1]);
                cy[i] += (N_force[j][2]+N_force[j][3]);
  //              S[i] += surfaceDVEPtr[j].S;             //adding DVE areas along chord
			}
			//Nondimensionalizing values using summed areas and velocity at LE of chordal row of DVEs
            tempI = j-panelPtr[k].n*(panelPtr[k].m-1); // index of DVE at leading edge
			tempS = 2/(spanPtr[i].area*dot(surfaceDVEPtr[tempI].u,surfaceDVEPtr[tempI].u));
            spanPtr[i].Cd = spanPtr[i].D_force*tempS;
			spanPtr[i].Cn *= tempS;
			cl[i] *= tempS;
			cy[i] *= tempS;
  printf("span %d  cn  %lf cl %lf cy %lf\n",i,spanPtr[i].Cn, cl[i],cy[i]);          
			i++;  //next span index 
			n++;	//index of next leading edge DVE 
		}
		n += panelPtr[k].n*(panelPtr[k].m-1);  //index of next LE DVE of next panel
	}
//===================================================================//
//			save spanwise information, lift and drag distribution
//===================================================================//
	//creates file name AOA#.##.txt ## is angle of attack
	sprintf(filename,"%s%s%.2lf%s",info.output,"AOA",info.alpha*RtD,".txt");

	//creates file in subdirectory output
	spaninfo = fopen(filename, "w");

	//write header
	fprintf(spaninfo,"This file contains spanwise information at");
	fprintf(spaninfo," angle of %.2lf\n",info.alpha*RtD);
	fprintf(spaninfo," The total lift (uncorrected for stall) and drag coefficient are ");
	fprintf(spaninfo," CL = %6.6lf (uncorrected)  CDi = %12.8lf  ",CL,CDi);
	fprintf(spaninfo," The horizontal tail is panel %d through %d\n"\
				,info.wing1[1],info.wing2[1]);
	
	//writes header for information on surface elements
	fprintf(spaninfo,"%6s %16s %16s %16s %16s %16s %16s %16s",\
	"index","xo","yo","zo","cn","cl","cy","cd");
	fprintf(spaninfo," %16s %16s %16s %16s %16s %16s %16s %16s %16s %16s\t#\n",\
	"A_TE","B_TE","C_TE","S","span","chord","nu","epsilon","psi","phiLE");

//	for(i=0;i<info.nospanelement;i++)
	i=0;		//intitaliing span index counter
	m=0;		//index of leading edge DVE
	for(k=0;k<info.nopanel;k++)   //loop over panels
	{
	  for(l=0;l<panelPtr[k].n;l++)  //loop over span of panel k
	  {
		if(i==info.wing1[1]) 
					fprintf(spaninfo,"HT\n");  //separates HT data
		j=m+(panelPtr[k].n*(panelPtr[k].m-1));
		//surface element index
		fprintf(spaninfo,"%6d",i);
		//coord. of ref point
		fprintf(spaninfo," %16.12lf %16.12lf %16.12lf",\
				surfacePtr[m].xo[0],surfacePtr[m].xo[1],\
				surfacePtr[m].xo[2]);
		//normal, lift, side, drag force coefficients
		fprintf(spaninfo," %16.12lf %16.12lf %16.12lf %16.12lf",\
				spanPtr[i].Cn,cl[i],cy[i],spanPtr[i].Cd);

		//more info on element
		fprintf(spaninfo," %16.12lf %16.12lf %16.12lf %16.12lf %16.12lf %16.12lf",\
				surfacePtr[j].A,surfacePtr[j].B,surfacePtr[j].C,spanPtr[i].area,\
				surfacePtr[m].eta*2,surfacePtr[m].xsi*2*panelPtr[k].m);
		fprintf(spaninfo," %16.12lf %16.12lf %16.12lf",\
				surfacePtr[m].nu*RtD,surfacePtr[m].epsilon*RtD,\
				surfacePtr[m].psi*RtD);
		fprintf(spaninfo," %16.12lf",\
				(surfacePtr[m].phiLE*RtD));
		fprintf(spaninfo,"\n");

		i++;  //next span index 
		m++;	//index of next leading edge DVE 
	  } 
	  m += panelPtr[k].n*(panelPtr[k].m-1);  //index of next LE DVE of next panel
	}


	fclose(spaninfo);
	*/
//===================================================================//
//		DONE saving spanwise information, lift and drag distribution
//===================================================================//

//===================================================================//
		//write trim results to TrimSol.txt
//===================================================================//
	fprintf(MomSol,"%6.2lf   %6.2lf  ",info.alpha*RtD,epsilonHT*RtD);
	fprintf(MomSol,"%6.3lf %12.8lf  %7.4lf  ",CL,CDi,Cm);
	fprintf(MomSol,"%10.6lf  %8.4lf\n",CLht,info.CMoWing);
	fflush(MomSol);

FREE1D(&cl,info.nospanelement);	//section lift coefficient
FREE1D(&cy,info.nospanelement);	//section side force coefficient
FREE1D(&mult, 2);	//will store diff aileron info

printf("\n");
}
//===================================================================//
		//END of program
//===================================================================//
