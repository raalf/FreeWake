//computes pitching moment
double PitchingMoment(const GENERAL,PANEL *,const double,DVE *&,\
						const double,const int,const double [3],\
						double &,double &,double **&,double *&,\
						double &,double &);

double PitchingMoment(const GENERAL info,PANEL *panelPtr,DVE *&surfacePtr,\
						const double cmac, const double epsilonHT,\
						const int HTpanel,const double xCG[3],\
						double &CLht,double &CLhti,\
						double **&N_force,double *&D_force,\
						double &CL,double &CDi_finit)
{
//this routine computes the pitching moment of a given configuration.
//The moment is computed about the point xCG, which moves with the wing
//The forces are computed using the Horstmann method and the distributed
//vorticity element method.
//
//input
//	info		general information
//	panelPtr	panel information
//  cmac 		mean aerodynamic chord of wing
//	epsilonHT	additional horizontal tail angle [rad]
//				positive - trailing edge down
//	HTpanel		index of first panel of HT
//	xCG			CG location
//
//output
//	surfacePtr	geometry of lifting surface
//	CM_resid	residutal pitching moment coefficient w/o CMo of wing
//	CLht		CL of HT, S is reference area
// 	CLhti		induced lift of HT
//	CL			aircraft CL
//	CDi_finit	induced drag
//
//
//  Wing 1 is the main wing, Wing 2 is the horizontal tail.
//  Thus,
//		Wing1 index = info.wing1[0]..info.wing2[0]
//  	Wing2 index = info.wing1[1]..info.wing2[1]
//
	int i,j,l,timestep=0;	// loop counter
	int saveStep=20;	//number of steps when relaxed wake is saved
	int timestart=0;	//first timestep of relaxed wake scheme
 //GB 2-9-20 	int HTindex=0;		//the first DVE index of the HT
	double CLi,CY,CYi;	//lift and side force coefficients
	double e_old;		//span efficiency of previous time step
	double deltae;		//square of delta_e of curent and previous time step
	double tempS,tempA[3];//a temporary variable for a scalar and vector
	double qc=1/(0.5*info.Uinf*info.Uinf*info.S*cmac);	//1/(dynamic pressure *cmac)
	double delX[3];		//vector from surface DVE ref. pt. to 1/4location
	double MomArm,deltaM;//moment arm of lift forces, moment of single DVE
	double Moment,CM_resid;//residual pitch moment and moment coefficient
	double XCG[3];		//CG location in this routine is moved with wing
				XCG[0] = xCG[0];	XCG[1] = xCG[1];	XCG[2] = xCG[2];

double *R,**D;			//resultant vector and matrix
int *pivot;				//holds information for pivoting D

//===================================================================//
		//Update tail incidence
//===================================================================//

    //add incicent angle increment
    if(info.trim==1)   //only if longitudinal trim routine is ON
//    {	for(i=HTpanel;i<info.nopanel;i++) changed GB 2-10-20
    {    for(i=info.panel1[1];i<=info.panel2[1];i++)
	{
		tempS = panelPtr[i].eps2 - panelPtr[i].eps1; //twist
		panelPtr[i].eps1 = epsilonHT;
		panelPtr[i].eps2 = epsilonHT+tempS;

		//adding number of spanwise DVEs of HT
 //removed GB 2-9-20 		HTindex += panelPtr[i].n;
	}
    }
 
	//computingn the first DVE index of HT
 //removed GB 2-9-20 	HTindex *= -info.m;	HTindex += info.noelement;

//===================================================================//
		//DONE Update tail incidence
//===================================================================//

//===================================================================//
//*******************************************************************//
//===================================================================//
//
//						HORSTMANN METHOD (fixed wake)
//
//The following part is a strongly abreviated part of Horstmann's 
//multiple lifting line method.  It doesn't work anymore, but the 
//information is needed for the relaxed wake part.
//         NO MORE HORSTMANN, ALL REMOVED GB 2-9-20
//===================================================================//
//*******************************************************************//
//===================================================================//

//===================================================================//
		//START generating surface Distributed-Vorticity Elements
//===================================================================//

	Surface_DVE_Generation(info,panelPtr,surfacePtr);
								//Subroutine in wing_geometry.cpp

	//save information on elementary wings to file
//	Save_Surface_DVEs(info,surfacePtr);	//Subroutine in write_output.cpp

//===================================================================//
		//END generating surface Distributed-Vorticity elements
//===================================================================//

	//allocate memory for relaxed wake part
	ALLOC1D(&pivot,info.Dsize);							//pivoting array
	ALLOC2D(&wakePtr,info.maxtime+1,info.nospanelement);	//wake DVE

	ALLOC1D(&CDi_DVE,info.maxtime+1);			//total induced drag (Eppler)
	ALLOC2D(&CN,info.maxtime+1,4);				//total normal forces

	//initalizing
	for(i=0; i<info.maxtime; i++)
	{
		CN[i][0] = 0;	CN[i][1] = 0;	CN[i][2] = 0;	CN[i][3] = 0;
		CDi_DVE[i] = 0;
	}
	for(j=0; j<info.nospanelement; j++)		D_force[j] = 0;

//===================================================================//
		//START generating D matrix
//===================================================================//
    
	//1. Assembyly of upper 2/3 of D-matrix using the boundary conditions
    //between the panels/elements
    //The new kinematic conditions due to the DVE is being recomputed.
	//2. The boundary conditions, as they were computed previously don't
	//change.
    //3. Decompose D-matrix in upper-lower diagonal matrices for faster
    //solving
    
    //allocates mememory for R and D
    ALLOC1D(&R,info.Dsize);
    ALLOC2D(&D,info.Dsize,info.Dsize);

    //initializing D
    for(i=0; i<info.Dsize; i++)
    for(j=0; j<info.Dsize; j++)
        D[i][j]=0.0;
    
    //1. assembly of first part of D, boundary cond. of
    //DVE within each panel
    DVE_BoundaryCond(surfacePtr,panelPtr,info,D);
                        //subroutine in equ_system.cpp

    //2. assemble new lower 1/3 of D-matrix
	DVE_KinCond(surfacePtr,info,panelPtr,D);
										//Subroutine in equ_system.cpp

	//decompose D-matrix into lower/upper matrix,
	//l/u coefficients saved in D -- WARNING: original D-values lost!!
	//pivot holds the pivoting information of D
	//also assignes zero values to appropriate elements of RHS-vector
	LU_Decomposition(D,info.Dsize,pivot);	//Subroutine in gauss.cpp

	//initalizing of R
	for (i=0; i<info.Dsize; i++)
		R[i]=0;
//===================================================================//
		//END generating new kinematic conditions for D matrix
//===================================================================//

//===================================================================//
		//START DVE vorticity distribution
//===================================================================//
//Assembles and solves an equation system that defines the vorticity
//distribution across each surface DV-elementary of the wing.
//First, the resultant vector is being assembled. The elements of
//its upper 2/3 are zeros and the lower 1/3 are the velocity component
//that resulst from the vorticity in the wake and from the free stream.
//only the component that is normal to the surface is considered (in
//accordance with the kinematic condition). The resultant vector
//needs to be assembled under considerations of the pivoting array.
//
//For each surface DVE, the circulation strength of their bound
//vortices is given by:	    	gamma(i) = A + B*etai + C*etai^2
//The vortex sheet inbetween has the strength B+2*etai*C.
//the function returns the coefficients A, B, C for each DVE

	DVE_Vorticity_Distribution\
					(info,panelPtr,surfacePtr,wakePtr,-1,D,R,pivot);
										//Subroutine in equ_system.cpp
//===================================================================//
		//END DVE vorticity distribution
//===================================================================//


	//initial values of timestep and CDiold
	timestep =timestart-1;
	e_old	 = 10; //initial CDi of 'previous' timestep

//	printf("working on timestep: "); //##

//\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
///////////////////////////////////////////////////////////////////////
//===================================================================//
//#######################Loop over time steps########################//
//===================================================================//
///////////////////////////////////////////////////////////////////////

	do
	{
		timestep ++; //advance timestep

		printf("%d ",timestep);
//		fflush(stdout);
//===================================================================//
		//START Move_Wing
//===================================================================//
//the every time step the wing is moved by delx which is deltime*u

		Move_Wing(info,surfacePtr);	//Subroutine in wing_geometry.cpp

		//move CG
		//newXCG -= local U * delta time
		XCG[0] -= surfacePtr[0].u[0] * info.deltime;
		XCG[1] -= surfacePtr[0].u[1] * info.deltime;
		XCG[2] -= surfacePtr[0].u[2] * info.deltime;

//===================================================================//
		//END Move_Wing
//===================================================================//

//===================================================================//
		//START Squirt_out_Wake
//===================================================================//
//creates most recent wake DV element just aft of trailing edge
//its vorticity coefficients, A,B, and C, are found when solving the
//equation system in DVE_Vorticity_Distribution

		Squirt_out_Wake(info,panelPtr,surfacePtr,wakePtr[timestep]);
									//Subroutine in wake_geometry.cpp

//===================================================================//
		//END Squirt_out_Wake
//===================================================================//

//===================================================================//
		//START New_vorticity_coefficients
//===================================================================//
//*		//updates vorticity in wake,
		//computes new vorticity coefficients A, B, and C
		if (info.steady == 1)	//steady airloads,
		{	//location have equally large amount of integrated circulation
			 j=timestep;
//			for(i=0;i<=timestep;i++)
			i=timestep;
			Update_wake_vorticity(info,panelPtr,wakePtr[i],wakePtr[j]);
	 	}								//subroutine in wake_geometry.cpp
		else //unsteady airloads, varying integrated circulation, k
		{ //for(i=0;i<=timestep;i++)
		  {
			i=timestep;
			j=i;
			Update_wake_vorticity(info,panelPtr,wakePtr[i],wakePtr[j]);
		  }									//subroutine in wake_geometry.cpp
		}
//*/
//===================================================================//
		//END New_vorticity_coefficients
//===================================================================//

//===================================================================//
		//START DVE vorticity distribution
//===================================================================//
//Assembles and solves an equation system that defines the vorticity
//distribution across each surface DV-elementary of the wing.
//First, the resultant vector is being assembled. The elements of
//its upper 2/3 are zeros and the lower 1/3 are the velocity component
//that resulst from the vorticity in the wake and from the free stream.
//only the component that is normal to the surface is considered (in
//accordance with the kinematic condition). The resultant vector
//needs to be assembled under considerations of the pivoting array.
//
//For each surface DVE, the circulation strength of their bound
//vortices is given by:	    	gamma(i) = A + B*etai + C*etai^2
//The vortex sheet inbetween has the strength B+2*etai*C.
//the function returns the coefficients A, B, C for each DVE
//*
		DVE_Vorticity_Distribution\
 				(info,panelPtr,surfacePtr,wakePtr,timestep,D,R,pivot);
										//Subroutine in equ_system.cpp

//===================================================================//
		//END DVE vorticity distribution
//===================================================================//

//===================================================================//
		//START New_vorticity_coefficients
//===================================================================//
//*		//updates vorticity in wake,
		//computes new vorticity coefficients A, B, and C
		if (info.steady == 1)	//steady airloads,
		{	//location have equally large amount of integrated circulation
			 j=timestep;
			for(i=0;i<=timestep;i++)
			Update_wake_vorticity(info,panelPtr,wakePtr[i],wakePtr[j]);
		}								//subroutine in wake_geometry.cpp
		else //unsteady airloads, varying integrated circulation, k
		{ for(i=0;i<=timestep;i++)
		  {
			j=i;
			Update_wake_vorticity(info,panelPtr,wakePtr[i],wakePtr[j]);
		  }									//subroutine in wake_geometry.cpp
		}
//*///===================================================================//
		//END New_vorticity_coefficients
//===================================================================//

//===================================================================//
							//START Relax_Wake
//===================================================================//
//relaxing the wake:
//	1. computes local induced velocity at side edges of DVEs
//	2. displaces of side edges of DVEs
//	3. computes new ref. pt of wake DVE
//	4. computes new eta, nu, epsilon, and psi, as well as new xsi
//	5. //DELETED AND REPLACE WITH ROUTINE IN MAIN 5/27/2005 G.B.
//	6. computes new singularity factor for

		//relax only after first two timesteps have been executed
		if(info.relax == 1 && timestep > 1)
			Relax_Wake(info,timestep,surfacePtr,wakePtr);
						 				//Subroutine in wake_geometry.cpp

//===================================================================//
							//END Relax_Wake
//===================================================================//

//===================================================================//
		//START New_vorticity_coefficients
//===================================================================//
//*		//updates vorticity in wake,
		//computes new vorticity coefficients A, B, and C
		if (info.steady == 1)	//steady airloads,
		{	//location have equally large amount of integrated circulation
			 j=timestep;
			for(i=0;i<=timestep;i++)
			Update_wake_vorticity(info,panelPtr,wakePtr[i],wakePtr[j]);
	 	}								//subroutine in wake_geometry.cpp
		else //unsteady airloads, varying integrated circulation, k
		{ for(i=0;i<=timestep;i++)
		  {
			j=i;
			Update_wake_vorticity(info,panelPtr,wakePtr[i],wakePtr[j]);
		  }									//subroutine in wake_geometry.cpp
		}
///*/
//===================================================================//
		//END New_vorticity_coefficients
//===================================================================//

//===================================================================//
		//START DVE lift computation
//===================================================================//
//*/		//computes normal forces/density for each surface DVE

		Surface_DVE_Normal_Forces(info,panelPtr,timestep,wakePtr,\
							  					surfacePtr,N_force);
							  		//Subroutine in lift_force.cpp

		//computes total lift and side force/density, and ascoefficients
		DVE_Wing_Normal_Forces(info,N_force,Nt_free, Nt_ind, CL,CLi,CY,CYi);
						 				//Subroutine in lift_force.cpp

//===================================================================//
		//END DVE lift computation
//===================================================================//

//===================================================================//
		//START Induce_DVE_Drag
//===================================================================//

		//	CDi			- total drag coefficient
		//  D_force 	- local drag force/density along span
		CDi_DVE[timestep] = \
		Induced_DVE_Drag(info,panelPtr,surfacePtr,wakePtr,timestep,D_force);
						 			//Subroutine in drag_force.cpp

//===================================================================//
		//END Induce_DVE_Drag
//===================================================================//*/

		//current span efficiency
		tempS = CL*CL/(Pi*info.AR*CDi_DVE[timestep]);

		//sqare of difference in span efficiencies
		deltae = e_old-tempS;  deltae = deltae*deltae;

		e_old = tempS; //save e of previous time step
//		printf("delta e %lf\n",sqrt(deltae));

	//continue time-stepping loop as long as
	// - e has not converged and
	// - maximum time steps have not been reached
	} while((deltae > info.deltae) && (timestep<info.maxtime));

///////////////////////////////////////////////////////////////////////
//===================================================================//
//#####################END Loop over time steps######################//
//##################### OR convergence of CDi  ######################//
//===================================================================//
///////////////////////////////////////////////////////////////////////
//\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

//===================================================================//
		//compute pitching moment
//===================================================================//

 	Moment = 0;//initializing Moment
 	//loop over no. surface DVE's

//printf("\n");
	for(i=0;i<info.noelement;i++)
	{
		//vector from surfaceDVE ref. pt. to 1/4chord location
		tempA[0]=-surfacePtr[i].xsi;	tempA[1]=0;		tempA[2]=0;
		Star_Glob(tempA,surfacePtr[i].nu,surfacePtr[i].epsilon,\
					surfacePtr[i].psi,delX);

		//computing the moment arm:
		MomArm = (XCG[0]-surfacePtr[i].xo[0]-delX[0])*cos(info.alpha)\
				+(XCG[2]-surfacePtr[i].xo[2]-delX[2])*sin(info.alpha);

		//moment increment is (L+Lind)*MomARm
		deltaM = (N_force[i][0]+N_force[i][1])*MomArm;

		//adding to total pitching moment/density
		Moment += deltaM;
	}

	//the residual-moment coefficient
	CM_resid = Moment*qc;  if(info.sym==1) CM_resid*=2;

//printf("CM resid %lf  \n",CM_resid);
//===================================================================//
		//DONE compute pitching moment
//===================================================================//

	CDi_finit = CDi_DVE[timestep];

	//computing lift forces of HT
	CLht =0; //initializing
	CLhti=0; //initializing
 //GB 2-9-20   for(i=HTindex;i<info.noelement;i++)
    for(i=info.dve1[1];i<=info.dve2[1];i++)
	{
		CLht  += N_force[i][0];
		CLhti += N_force[i][1];
	}

    //non-dimensional
	CLht *= qc*cmac;  CLhti *= qc*cmac;	if(info.sym==1) {CLht*=2; CLhti*=2;}

//===================================================================//
//*******************************************************************//
//===================================================================//
//===================================================================//
				//save final wake shape results to files//
//===================================================================//

		//saves surface and wake information of last timestep to file
		Save_Timestep(info,timestep,wakePtr,surfacePtr,N_force);
											//Subroutine in write_output.cpp
//===================================================================//
						//DONE save results//
//===================================================================//


	//free allocated memory
	FREE1D(&R,info.Dsize);
	FREE2D(&D,info.Dsize,info.Dsize);
	FREE1D(&pivot,info.Dsize);

	FREE1D(&CDi_DVE,info.maxtime+1);
	FREE2D(&CN,info.maxtime+1,4);
	FREE2D(&wakePtr,info.maxtime+1,info.nospanelement);
	
	//returning the residual moment coefficient
	return(CM_resid);
}
//===================================================================//
		//END of program
//===================================================================//
