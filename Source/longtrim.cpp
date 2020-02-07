void LongitudinalTrim(GENERAL,PANEL *,DVE *,int,double *&,double &,\
						double &,FILE *);

void LongitudinalTrim(GENERAL info,PANEL *panelPtr,DVE *surfaceDVEPtr,int HTpanel,\
						double *&cn,double &CL,double &CDi,\
						FILE *MomSol)
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
//	HTpanel		index of first panel of HT
//
//output
//	surfaceDVEPtr	geometry of surface elements
//	cn			normal force coefficient of each span element
//	cd			drag force coefficient of each span element
// 	CL			total lift
//	CDi			total induced drag
//
//the routine also saves the trim information to TrimSol.txt
//furthermore saved is the spanwise lift distribution and spanwise
//drag distribution in SpanAOA<info.alpha>.txt where <info.alpha> 
//holds the numerical value of the current angle of attack 

	int i,j,k,l,m,n=0;				//loop counter
	double epsilonHT;	//HT incident correction [rad]
	double CM_resid;	//residual moment
	double eps1;		//the HT incidence angle of one iteration previous
	double CM_old;		//saving pitching moment of previous iteration step
	double CLht,CLhti;	//total and induced lift of HT
	double *cl,*cy;		//section lift and side force coefficients
	double *S;			//area of a spanwise strip (sum of areas of DVEs of one span location)
	double tempS;
	double **N_force;		//surface DVE's normal forces/density
					//[0]: free stream lift, [1]: induced lift,
					//[2]: free stream side, [3]: induced side force/density
					//[4]: free str. normal, [5]: ind. normal frc/density
	double *D_force;		//drag forces/density along span
	double *cd;	//section ind. drag coefficient

	FILE *spaninfo;			//output file for spanwise information
	char filename[133];	//file path and name for spanwise information

	//allocating memory	
	ALLOC1D(&cl,info.noelement);	//section lift coefficient
   	ALLOC1D(&cy,info.nospanelement);	//section side force coefficient
   	ALLOC1D(&S,info.nospanelement);	//section area
   	ALLOC1D(&cd,info.nospanelement);	//section ind. drag coefficient
	ALLOC2D(&N_force,info.noelement,6);	//surface DVE normal forces
	ALLOC1D(&D_force,info.nospanelement);//Drag force per span element
	
	//initial HT incident correction
	epsilonHT = 0;//panelPtr[i].eps1;  	//[rad]
	eps1 = 0;	//the HT incidence angle of one iteration previous

	//===========================================//
		//pitching moment computation - Step 0
	//===========================================//

	//computed the residual moment coefficient of wing and tail
	CM_resid = PitchingMoment\
				(info,panelPtr,surfacePtr,info.cmac,epsilonHT,\
				HTpanel,info.RefPt,CLht,CLhti,\
				N_force,D_force,CL,CDi);//Subroutine in PitchMoment.cpp

	//adding zero lift of wing	
	CM_resid += info.CMoWing;

//printf("\neps %.1lf CMresid %.4lf CMold %.4lf",epsilonHT*RtD,CM_resid,CM_old);
	if(info.trim==1)   //longitudinal trim routine
	{
	  if(CM_resid > 0)  	epsilonHT = 0.035; //deflect HT trailing edge 2deg down
	  else 				epsilonHT =-0.035; //deflect HT trailing edge 2deg up

	//========================================//
	//	pitching moment iteration loop
	//========================================//
	  i=0;	//initializing loop counter
	  do
	  {
	//========================================//
	//	pitching moment computation - Step i
	//========================================//
		//saving pitching moment of previous iteration step
		CM_old = CM_resid;
		//computed the residual moment coefficient of wing and tail
		CM_resid = PitchingMoment\
					(info,panelPtr,surfacePtr,info.cmac,epsilonHT,\
					HTpanel,info.RefPt,CLht,CLhti,\
					N_force,D_force,CL,CDi);//Subroutine in PitchMoment.cpp

		//adding zero lift of wing	
		CM_resid += info.CMoWing;

		//computing new HT incidence angle
		tempS = epsilonHT - CM_resid*(epsilonHT-eps1)/(CM_resid-CM_old);
		eps1 = epsilonHT;  epsilonHT = tempS;	//reassigning HT angles
			
		i++;	//incrementing loop counter

	  }while((i<20) && (CM_resid*CM_resid>.000025));
	//========================================//
	//	pitching moment iteration loop  END
	//========================================*///
	}   //end of longitudinal trim routine

//===================================================================//
//		Compute cn and cd (spanwise values)
//===================================================================//

	i=0;		//intitaliing span index counter
	m=0;		//index of leading edge DVE
	for(k=0;k<info.nopanel;k++)  //loop over panels
	{
		for(l=0;l<panelPtr[k].n;l++)  //loop over span of panel k
		{
			//working each spanwise strip
			//initializing
			cn[i]=0;
			cl[i]=0;
			cy[i]=0;
			S[i] = 0;  //area of current spanwise strip
			
			//adding up chordal values of one spanwise location (indexed i)
			for(m=0;m<info.m;m++)
			{
				j=n+m*panelPtr[k].n;
				cn[i] += (N_force[j][4]+N_force[j][5]);
				cl[i] += (N_force[j][0]+N_force[j][1]);
				cy[i] += (N_force[j][2]+N_force[j][3]);
				S[i] += surfaceDVEPtr[j].S;
			}
			//Nondimensionalizing values
			tempS = 2/(info.Uinf*info.Uinf*S[i]);			
			cd[i] = D_force[i]*tempS;		
			cn[i] *= tempS;
			cl[i] *= tempS;
			cy[i] *= tempS;	

			i++;  //next span index 
			n++;	//index of next leading edge DVE 
		}
		n += panelPtr[k].n*(info.m-1);  //index of next LE DVE of next panel
	}


//===================================================================//
//			save spanwise information, lift and drag distribution
//===================================================================//
	//creates file name AOA#.##.txt ## is angle of attack
	sprintf(filename,"%s%s%.2lf%s",OUTPUT_PATH,"AOA",info.alpha*RtD,".txt");

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
	"A","B","C","S","span","chord","nu","epsilon","psi","phiLE");

//	for(i=0;i<info.nospanelement;i++)
	i=0;		//intitaliing span index counter
	m=0;		//index of leading edge DVE
	for(k=0;k<info.nopanel;k++)   //loop over panels
	{
	  for(l=0;l<panelPtr[k].n;l++)  //loop over span of panel k
	  {
		if(i==info.wing1[1]) 
					fprintf(spaninfo,"HT\n");  //separates HT data
		
		//surface element index
		fprintf(spaninfo,"%6d",i);
		//coord. of ref point
		fprintf(spaninfo," %16.12lf %16.12lf %16.12lf",\
				surfacePtr[m].xo[0],surfacePtr[m].xo[1],\
				surfacePtr[m].xo[2]);
		//normal, lift, side, drag force coefficients
		fprintf(spaninfo," %16.12lf %16.12lf %16.12lf %16.12lf",\
				cn[i],cl[i],cy[i],cd[i]);

		//more info on element
		fprintf(spaninfo," %16.12lf %16.12lf %16.12lf %16.12lf %16.12lf %16.12lf",\
				surfacePtr[m].A,surfacePtr[m].B,surfacePtr[m].C,S[i],\
				surfacePtr[m].eta*2,surfacePtr[m].xsi*2*info.m);
		fprintf(spaninfo," %16.12lf %16.12lf %16.12lf",\
				surfacePtr[m].nu*RtD,surfacePtr[m].epsilon*RtD,\
				surfacePtr[m].psi*RtD);
		fprintf(spaninfo," %16.12lf",\
				(surfacePtr[m].phiLE*RtD));
		fprintf(spaninfo,"\n");

		i++;  //next span index 
		m++;	//index of next leading edge DVE 
	  } 
	  m += panelPtr[k].n*(info.m-1);  //index of next LE DVE of next panel
	}


	fclose(spaninfo);
//===================================================================//
//		DONE saving spanwise information, lift and drag distribution
//===================================================================//

//===================================================================//
		//write trim results to TrimSol.txt
//===================================================================//
	fprintf(MomSol,"%6.2lf   %6.2lf  ",info.alpha*RtD,epsilonHT*RtD);
	fprintf(MomSol,"%6.3lf %12.8lf  %7.4lf  ",CL,CDi,CM_resid);
	fprintf(MomSol,"%10.6lf  %8.4lf\n",CLht,info.CMoWing);
	fflush(MomSol);



FREE2D(&N_force,info.noelement,6);
FREE1D(&D_force,info.nospanelement);
FREE1D(&cl,info.nospanelement);	//section lift coefficient
FREE1D(&cy,info.nospanelement);	//section side force coefficient
FREE1D(&S,info.nospanelement);	//section area
FREE1D(&cd,info.nospanelement);	//section side force coefficient

printf("\n");
}
//===================================================================//
		//END of program
//===================================================================//
