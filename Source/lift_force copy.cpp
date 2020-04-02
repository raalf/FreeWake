//computes total normal force acting on wing
void DVE_Wing_Normal_Forces(const GENERAL,double **,double [2],double [2],\
							double &,double &,double &,double &);
//computes normal forces of surface DVE
void Surface_DVE_Normal_Forces(const GENERAL,const PANEL *,const int,\
							   DVE **,DVE *,double **);
//computes normal forces of elementary wings, fixed wake model
void Elementary_Wing_Normal_Forces\
					(const PANEL *,const GENERAL, BOUND_VORTEX *);
//computes normal forces of wing, fixed wake model
void Wing_Normal_Forces(const  PANEL* ,const GENERAL,const BOUND_VORTEX*,\
						double [2],double [2],\
						double &,double &,double &,double &);

//===================================================================//
		//START FUNCTION DVE_Wing_Normal_Forces
//===================================================================//
void DVE_Wing_Normal_Forces(const GENERAL info,double **N_force,\
							double Nt_free[2], double Nt_ind[2],\
							double &CL,double &CLi,double &CY,double &CYi)
{
//this routine adds up the DVE's normal forces in order to compute the
//total wing normal forces/density and coefficients based on free stream
//input:
// info		general information on case
// N_force	normal forces/density of each surface DVE, second index is:
//			[0]: free stream lift, [1]: induced lift,
//			[2]: free stream side, [3]: induced side force/density

//
//ouput:
// as part of surfacePtr.:
// N_free	total lift and side forces/density due to free stream flow
// N_ind	total lift and side forces/density due to induced velocities
//
// Nt_free	total lift and side forces/density due to free stream flow
// Nt_ind	total lift and side forces/density due to induced velocities
// CL		total lift coefficient
// CLi		total induced lift coefficient
// CY		total side-force coefficient
// CYi		total induced side-force coefficient

  int l=0;									//counter
  double q=0.5*info.Uinf*info.Uinf*info.S; 	//ref. area* dyn. pressure/density

  Nt_free[0]=0;
  Nt_free[1]=0;
  Nt_ind[0]=0;
  Nt_ind[1]=0;

  //loop over number of panels
  for (l=0;l<info.noelement;l++)
  {
	  //adding the normal forces/density of all elementary wings
	  Nt_free[0]	+=  N_force[l][0];
	  Nt_free[1]	+=  N_force[l][2];

	  Nt_ind[0]	+=  N_force[l][1];
	  Nt_ind[1]	+=  N_force[l][3];;

//#printf("NtX  =%lf\t NtZ  =%lf\n",Nt_free[0],Nt_free[1]);//#
//#printf("NtXi =%lf\t NtZi =%lf\n",Nt_ind[0],Nt_ind[1]);//#
  }

  if (info.sym==1 && info.beta == 0)
  {	//twice the force if symmetric geometry
	  Nt_free[0]*=2;

	  Nt_ind[0]	*=2;
  	  Nt_ind[1]	*=2;
//#printf("NtX  =%lf\t NtZ  =%lf\n",Nt_free[0],Nt_free[1]);//#
  }
  //total lift and side force coefficients
  CL = (Nt_free[0]+Nt_ind[0])/q;
  CY = (Nt_free[1]+Nt_ind[1])/q;

  //total induced lift and side force coefficients
  CLi = Nt_ind[0]/q;
  CYi = Nt_ind[1]/q;

//#printf("CL=%lf\tCLi=%lf\tCY=%lf\tCYi=%lf\n",CL,CLi,CY,CYi);//#

}
//===================================================================//
		//END FUNCTION DVE_Wing_Normal_Forces
//===================================================================//

//===================================================================//
		//FUNCTION Surface_DVE_Normal_Forces
//===================================================================//
void Surface_DVE_Normal_Forces(const GENERAL info,const PANEL *panelPtr,\
							   const int timestep,DVE **wakePtr,\
							   DVE *surfacePtr,double **N_force)
{
//computes lift and side force/density acting on surface DVE's. The local
//lift force is comuted by applying Kutta-Joukowski's theorem to both edges
//and to the center of DVE's bound vortices. These values are integrated,
//using Simpson's rule (with overhang), in order to get the lift force/density
//for the complete surface DVE.  The computational effort is reduced by
//combining the vorticities of two neighboring aft and forwrd bound vortices.
//Furthermore, as of now, any lift force due to the DVE's vortex sheet is
//neglected, since the resulting induced side velocities should be small
//
//
//Note: the velocities are not computed directly at the vortex edges, but
//(1-delta)/2 of the DVE span further towards the center.
//Otherwise, the computed induced velocity would become singular at that
//point due to the next neighboring elementary wing influence.
//
//input:
//	info			general information on case
//	timestep		current time step
// 	wakePtr			wake DVE information
// 	surfacePtr		surface DVE information
//
//ouput:
//	N_force[l][]	l-th surface DVE's normal forces/density
//					[0]: free stream lift, [1]: induced lift,
//					[2]: free stream side, [3]: induced side force/density
//					[4]: free stream normal, [5]: induced normal force/density
//
//
int i,j,k;					//i'th panels
							//j'th chordwise lifting line, j=0..(info.m-1)
							//k'th spanwise elementary wing, k=0..(n-1)
int l=0;					//index of surface DVE whose forces are computed
double A, B, C;				//vorticity distribution coefficient of el. l
double eta, eta8;			//half span of elementary wing l, 90% value of eta
double sbeta=sin(info.beta),\
	   cbeta=cos(info.beta);//sine and cosine of beta
double xoLE[3];				//forward DVE bound vortex mid point
double X[3];				//vector from forward DVE bound vortex mid point
							//to its edge; x1=xoLE-x and x2=xoLE+x
double S[3];				//vector of bound vortex
double UxS;					//abs. value of U x S
double eN[3];				//normal force direction
double eL[3],eS[3];			//lift force and side force direction
double w1[3],wo[3],w2[3];	//ind. vel. at bound vortex sides and center
double gamma1,gammao,gamma2;//vorticity at bound vortex sides and center
double R1[3],Ro[3],R2[3]; 	//resultant ind. force at bound vortex sides and center
double R[3];				//resultant ind. force/density of element l
double N_free;				//magnitude free stream norm. forces/density
double tempA[3],tempAA[3], tempS;
double **U1,**Uo,**U2;		//mid-chord velocities of upstream DVE
double spandir[3];
double Uxspandir;	


//loop over number of panels
for (i=0;i<info.nopanel;i++)
{
  if(panelPtr[i].m>1) //if more than one lifting line
  {    //allocating temporary memory for the induced velocity of upstream DVE
        //needed for averaging velocity induced at lifting lines
        //added GB 2-9-20
        ALLOC2D(&U1,panelPtr[i].n,3);
        ALLOC2D(&Uo,panelPtr[i].n,3);
        ALLOC2D(&U2,panelPtr[i].n,3);
  }

  //loop over number of chordwise elements 'info.m'
  for (j=0;j<panelPtr[i].m;j++)
  {
	 //loop over number of spanwise elements 'n'
	 for (k=0;k<panelPtr[i].n;k++)
	 {
		if(j==0)
		{
			//vector of bound vortex, KHH eq. 66
			//first allong the leading edge
			tempA[0]=tan(surfacePtr[l].phiLE); tempA[1]=1; tempA[2]=0;
			//transforming into local reference frame
			Star_Glob(tempA,surfacePtr[l].nu,surfacePtr[l].epsilon,\
					surfacePtr[l].psi,S); //function in ref_frame_transform.h
		}
		else
		{
			//vector along the bound vortex along LE
			tempA[0]=tan(surfacePtr[l].phi0); tempA[1]=1; tempA[2]=0;
			//transforming into local reference frame
			Star_Glob(tempA,surfacePtr[l].nu,surfacePtr[l].epsilon,\
												surfacePtr[l].psi,S);
								//function in ref_frame_transform.h
		}

		// S calculation works for circling flight D.F.B. 03-2020
		eta  = surfacePtr[l].eta;
		eta8 = eta*0.8;

		//#non-small angles
		//normal force direction
		cross(surfacePtr[l].u,S,tempA);			//#	UxS
		UxS=norm2(tempA);						//	|UxS|
		scalar(tempA,1/UxS,eN);					//	eN=(UxS)/|UxS|
       
        //*** Calculated lift direction based on the freestream velocity direction
		//		Changed by D.F.B. 03-2020
		//vector along the bound vortex along LE
		tempA[0]=0; tempA[1]=1; tempA[2]=0;
		//transforming into local reference frame
		Star_Glob(tempA,0,surfacePtr[l].epsilon,surfacePtr[l].psi,spandir);

		cross(surfacePtr[l].u,spandir,tempA);		//#	UxS
		Uxspandir=norm2(tempA);						//	|UxS|
		scalar(tempA,1/Uxspandir,eL);				//	eL=(UxS)/|UxS|

//        printf("eL\t%f %f %f\n",eL[0],eL[1],eL[2]);
//        printf("spandir\t%f %f %f\n",spandir[0],spandir[1],spandir[2]);
//        printf("DVEu\t%f %f %f\n",surfacePtr[l].u[0],surfacePtr[l].u[1],surfacePtr[l].u[2]);
         

		//***Removed by D.F.B. 03-2020
		//the lift direction  eL=Ux[0,1,0]/|Ux[0,1,0]|
		/*tempS = sqrt(surfacePtr[l].u[0]*surfacePtr[l].u[0]\
		 			+surfacePtr[l].u[2]*surfacePtr[l].u[2]);
		eL[0] = -surfacePtr[l].u[2]/tempS;
		eL[1] =  0;
		eL[2] =  surfacePtr[l].u[0]/tempS;*/
		//printf("eL\t%f %f %f\n",eL[0],eL[1],eL[2]);
		//printf("eN\t%f %f %f\n",eN[0],eN[1],eN[2]);
		//if(i==0 & j==0 & k ==0){CreateQuiverFile(surfacePtr[l].xo, eL,0);}
		//else{CreateQuiverFile(surfacePtr[l].xo, eL,1);}

		//CreateQuiverFile(surfacePtr[l].xo, eL,1);

		//the side force direction eS=UxeL/|UxeL|
		// cross(eL,surfacePtr[l].u,tempA);  \\ Removed by D.F.B. 03-2020
		// Original FW use eL x U... Changed to U x eL :
		cross(surfacePtr[l].u,eL,tempA); 
		tempS=1/norm2(tempA);
		scalar(tempA,tempS,eS);


//#printf(" U = %lf\t%lf\t%lf\n",surfacePtr[l].u[0],surfacePtr[l].u[1],surfacePtr[l].u[2]);//#
//#printf("S = %lf  %lf  %lf  %lf\n",S[0],S[1],S[2],norm2(S));//#
//printf("eN= %lf\t%lf\t%lf\t%lf\n",eN[0],eN[1],eN[2],norm2(eN));//#
//printf("eL= %lf\t%lf\t%lf\t%lf\n",eL[0],eL[1],eL[2],norm2(eL));//#
//printf("eS= %lf\t%lf\t%lf\t%lf\n",eS[0],eS[1],eS[2],norm2(eS));//#
//printf("UXS= %lf\n",(UxS));//#

		if(j==0)
		{	//most forward bound vortex
			A = surfacePtr[l].A;
			B = surfacePtr[l].B;
		 	C = surfacePtr[l].C;
		}
		else
		{	//combined circulation of l-th and (l-n)-th DVE's
			A = surfacePtr[l].A - surfacePtr[l-panelPtr[i].n].A;
			B = surfacePtr[l].B - surfacePtr[l-panelPtr[i].n].B;
			C = surfacePtr[l].C - surfacePtr[l-panelPtr[i].n].C;
		}

//*****************************************************************************
		//computing magnitude of normal force/density due to free stream
//*****************************************************************************
		N_free = (A*2*eta + C/3*2*eta*eta*eta)*UxS;
		// N_free calculation works for circling flight D.F.B. 03-2020
//#printf("N_free =%lf\t L_free =%lf\n",N_free,2*N_free*sqrt(eN[0]*eN[0]+eN[2]*eN[2]));//#

//*****************************************************************************
		 //computing the induced force/density
//*****************************************************************************
		//computing the ind. velocity at left (1) edge of bound vortex
		//vector from mid point of elementary wing bound vortex
		//to edge (1) of bound vortex; x1=xo-x and x2=xo+x
		//due to singular behavior of velocity at element edge
		//velocity is computed 0.1eta away from edge (hence factor 0.8)

		if(j==0)
		{
			//computing the leading edge midpoint
            //first the vector to leading edge
            tempA[0]=-surfacePtr[l].xsi;tempA[1]=0;tempA[2]=0;
            //transforming into local reference frame
            Star_Glob(tempA,surfacePtr[l].nu,surfacePtr[l].epsilon,\
                                            surfacePtr[l].psi,tempAA);
										//function in ref_frame_transform.h
			xoLE[0] = surfacePtr[l].xo[0] + tempAA[0];
			xoLE[1] = surfacePtr[l].xo[1] + tempAA[1];
			xoLE[2] = surfacePtr[l].xo[2] + tempAA[2];

			//computing the ind. velocity at left (1) edge of bound vortex
			scalar(S,-eta8,X);
			vsum(xoLE,X,tempA);
			DVE_Induced_Velocity(info,tempA,surfacePtr,wakePtr,timestep,w1);
				 					//subroutine in induced_velocity.cpp

		 	//computing the ind. velocity at center (0) of bound vortex
			DVE_Induced_Velocity(info,xoLE,surfacePtr,wakePtr,timestep,wo);
					 					//subroutine in induced_velocity.cpp

			//computing the ind. velocity at right (2) edge of bound vortex
			//vector from mid point of elementary wing bound vortex
			//to edge (2) of bound vortex; x1=xo-x and x2=xo+x
			//due to singular behavior of velocity at element edge
			//velocity is computed 0.1eta away from edge (hence factor 0.8)
			scalar(S,eta8,X);
			vsum(xoLE,X,tempA);
			DVE_Induced_Velocity(info,tempA,surfacePtr,wakePtr,timestep,w2);
					 					//subroutine in induced_velocity.cpp

            if(panelPtr[i].m>1)//if more than one lifting line, compute and store
			{			//velocities at half-chord loaction for averaging later
			  xoLE[0] = surfacePtr[l].xo[0];
			  xoLE[1] = surfacePtr[l].xo[1];
			  xoLE[2] = surfacePtr[l].xo[2];

			//###################################
			//added 2/6/15 G.B.
			  //recomputing S along midchord of first row of DVEs in preparation
			  //for averaging with velocities of second row of DVEs ->

			  //vector along the bound vortex along LE
			  tempA[0]=tan(surfacePtr[l].phi0); tempA[1]=1; tempA[2]=0;
			  //transforming into local reference frame  -> S saved in tempAA
			  Star_Glob(tempA,surfacePtr[l].nu,surfacePtr[l].epsilon,\
												surfacePtr[l].psi,tempAA);
								//function in ref_frame_transform.h
			//###################################

			  //computing the ind. velocity at left (1) edge of bound vortex
			  scalar(tempAA,-eta8,tempA);
			  vsum(xoLE,X,tempA);
			  DVE_Induced_Velocity\
			   					(info,tempA,surfacePtr,wakePtr,timestep,U1[k]);
					 					//subroutine in induced_velocity.cpp

		 	  //computing the ind. velocity at center (0) of bound vortex
			  DVE_Induced_Velocity\
			  					(info,xoLE,surfacePtr,wakePtr,timestep,Uo[k]);
					 					//subroutine in induced_velocity.cpp

			  //computing the ind. velocity at right (2) edge of bound vortex
			  //vector from mid point of elementary wing bound vortex
			  //to edge (2) of bound vortex; x1=xo-x and x2=xo+x
			  //due to singular behavior of velocity at element edge
			  //velocity is computed 0.1eta away from edge (hence factor 0.8)
			  scalar(tempAA,eta8,X);
			  vsum(xoLE,X,tempA);
			  DVE_Induced_Velocity\
			  					(info,tempA,surfacePtr,wakePtr,timestep,U2[k]);
					 					//subroutine in induced_velocity.cpp
			}

		}
		else  //i.e. j>0; circulation of l-th and (l-n)-th DVE's combined
		{
		//case of multiple lifting lines along the span
		//the induced velocity at the lifting line is averaged with the
		//velocities at mid chord locations of the DVES upstream and
		//downstream of the bound vortex. Otherwise, the singularity of
		//the bound vortex and the discontinuity of the bound vortex sheet
		//of a wing with twist causes trouble.  G.B. 1/24/06

			xoLE[0] = surfacePtr[l].xo[0];
			xoLE[1] = surfacePtr[l].xo[1];
			xoLE[2] = surfacePtr[l].xo[2];

			//computing the ind. velocity at left (1) edge of bound vortex
			scalar(S,-eta8,X);
			vsum(xoLE,X,tempA);
			DVE_Induced_Velocity(info,tempA,surfacePtr,wakePtr,timestep,tempAA);
					 					//subroutine in induced_velocity.cpp
			//averaging velocity and reassigning velocity to temporary variable
			w1[0] = (U1[k][0]+tempAA[0])*0.5; U1[k][0] = tempAA[0];
			w1[1] = (U1[k][1]+tempAA[1])*0.5; U1[k][1] = tempAA[1];
			w1[2] = (U1[k][2]+tempAA[2])*0.5; U1[k][2] = tempAA[2];

		 	//computing the ind. velocity at center (0) of bound vortex
			DVE_Induced_Velocity(info,xoLE,surfacePtr,wakePtr,timestep,tempAA);
					 					//subroutine in induced_velocity.cpp
			//averaging velocity and reassigning velocity to temporary variable
			wo[0] = (Uo[k][0]+tempAA[0])*0.5; Uo[k][0] = tempAA[0];
			wo[1] = (Uo[k][1]+tempAA[1])*0.5; Uo[k][1] = tempAA[1];
			wo[2] = (Uo[k][2]+tempAA[2])*0.5; Uo[k][2] = tempAA[2];

			//computing the ind. velocity at right (2) edge of bound vortex
			//vector from mid point of elementary wing bound vortex
			//to edge (2) of bound vortex; x1=xo-x and x2=xo+x
			//due to singular behavior of velocity at element edge
			//velocity is computed 0.1eta away from edge (hence factor 0.8)
			scalar(S,eta8,X);
			vsum(xoLE,X,tempA);
			//DVE_Induced_Velocity(info,tempA,surfacePtr,wakePtr,timestep,w2);
			DVE_Induced_Velocity(info,tempA,surfacePtr,wakePtr,timestep,tempAA);
				 					//subroutine in induced_velocity.cpp
			//averaging velocity and reassigning velocity to temporary variable
			w2[0] = (U2[k][0]+tempAA[0])*0.5; U2[k][0] = tempAA[0];
			w2[1] = (U2[k][1]+tempAA[1])*0.5; U2[k][1] = tempAA[1];
			w2[2] = (U2[k][2]+tempAA[2])*0.5; U2[k][2] = tempAA[2];

		}//done computing the velocities that are induced at the bound vortex.

	  //Integration of induced forces with Simpson's Rule
	  //Integration requires overhanging edges!!
	  //See also KHH linees 2953 - 2967, A23SIM

		//Kutta-Joukowski at left (1) edge
		cross(w1,S,tempA);				// w1xS
		gamma1  = A-B*eta8+C*eta8*eta8;//gamma1
		scalar(tempA,gamma1,R1);

		//Kutta-Joukowski at center
		cross(wo,S,tempA);				// woxS
		gammao  =	A;
		scalar(tempA,gammao,Ro);

 		//Kutta-Joukowski at right (2) edge
 		cross(w2,S,tempA);				// w2xS
		gamma2  =	A+B*eta8+C*eta8*eta8;
		scalar(tempA,gamma2,R2);

//#printf("gamma %lf\t%lf\t%lf\n",gamma1,gammao,gamma2);//#
//#printf("R1     %lf\t%lf\t%lf\n",R1[0],R1[1],R1[2]);//#
//#printf("Ro     %lf\t%lf\t%lf\n",Ro[0],Ro[1],Ro[2]);//#
//#printf("R2     %lf\t%lf\t%lf\n",R2[0],R2[1],R2[2]);//#

		//The resultierende induced force of element l is
		//determined by numerically integrating forces acros element
		//using Simpson's Rule with overhaning parts
		R[0]  = (R1[0]+4*Ro[0]+R2[0])*eta8/3;			//Rx
		R[1]  = (R1[1]+4*Ro[1]+R2[1])*eta8/3;			//Ry
		R[2]  = (R1[2]+4*Ro[2]+R2[2])*eta8/3;			//Rz

		//plus overhanging parts
		R[0] += (7*R1[0]-8*Ro[0]+7*R2[0])*(eta-eta8)/3;//Rx
		R[1] += (7*R1[1]-8*Ro[1]+7*R2[1])*(eta-eta8)/3;//Ry
		R[2] += (7*R1[2]-8*Ro[2]+7*R2[2])*(eta-eta8)/3;//Rz
//#printf("R[%d] %lf %lf %lf ",l,R[0],R[1],R[2]);//#
//*****************************************************************************
		 //the NORMAL FORCE/density
//*****************************************************************************
		//free stream lift
		N_force[l][4] = N_free;
		N_force[l][5] = dot(R,eN);			//induced normal force

//*****************************************************************************
		 //the LIFT FORCE/density is the normal force in the x-z plane or
//*****************************************************************************
		//free stream lift
		N_force[l][0] = N_free * sqrt(eN[0]*eN[0]+eN[2]*eN[2]);
		if (eN[2]<0)  N_force[l][0] *= -1;	//neg. if resultend is downward

		N_force[l][1] = dot(R,eL);			//induced lift
//#printf("N_free =%lf\t L_free =%lf\n",N_free,2*N_free*sqrt(eN[0]*eN[0]+eN[2]*eN[2]));//#
		double tempN_FREE[3];
		scalar(eN,N_free/50,tempN_FREE);
		if(i==0 & j==0 & k ==0)
            {CreateQuiverFile(surfacePtr[l].xo, tempN_FREE,0);}
		else{CreateQuiverFile(surfacePtr[l].xo, tempN_FREE,1);}
//*****************************************************************************
	  	 //the SIDE FORCE/density is the force in y-direction or N*eN[y]
//*****************************************************************************
		N_force[l][2] = N_free * eN[1];	//free stream side force
		N_force[l][3] = dot(R,eS);		//inducd side force

		l++; //increment elementary wing index l=0..(noelement-1)
	 }	//End loop over k - loop over span of panel
  }	//End loop over j - loop over chord of panel
    
  if(panelPtr[i].m>1)
  {    //temporary variable of induced velocity of upstream DVE
      // added GB 2-9-20
        FREE2D(&U1,panelPtr[i].n,3);
        FREE2D(&Uo,panelPtr[i].n,3);
        FREE2D(&U2,panelPtr[i].n,3);
  }
} //End loop over i - loop over panels

}
//===================================================================//
		//END FUNCTION Surface_DVE_Normal_Forces
//===================================================================//

//===================================================================//
		//FUNCTION Elementary_Wing_Normal_Forces
		//computes lift and side force/density acting on elem. wings
//===================================================================//
void Elementary_Wing_Normal_Forces(const PANEL* panelPtr, const GENERAL info,\
								   BOUND_VORTEX* elementPtr)
{
//Function computes normal forces/density for each elementary wing
//to do this, Kutta-Joukowski's theorem is applied at both edges and the
//center of the elementary wing in order to get the local lift forces. These
//values are integrated, using Simpson's rule (with overhang), in order to
//get the lift force/density for the complete elementary wing.
//Note: the velocities are not computed directly at the edges, but
//(1-delta)/2 of the elementary span further towards the center.
//Otherwise, the computed induced velocity would become singular at that
//point due to the next neighboring elementary wing influence.
//
//input:
// 	panel		-panel info, especially number of chord and span devisions
//	info		-general information on case
// 	elementPtr	-pointer to elementary-wing information
//	l			-index of elementary wing whose forces are being computed
//				 incremented by one at end of loop 'k' over number of
//				 spanwise elementary wings. Hence, l = 0 .. (noelement-1)
//
//ouput:
// as part of elementPtr.:
//	N_free		-lift and side forces/density due to free stream flow
//	N_ind		-lift and side forces/density due to induced velocities

int i,j,k;					//i'th panels
							//j'th chordwise lifting line, j=0..(info.m-1)
							//k'th spanwise elementary wing, k=0..(n-1)
int l=0;					//index of elementary wing whose forces are computed
							//index is incremented by one at end of loop
							//'k' over number of spanwise elementary wings.
							// Hence, l = 0 .. (noelement-1)
double A, B, C;				//vorticity distribution coefficient of el. l
double eta, eta8;			//half span of elementary wing l, 90% value of eta
double sbeta=sin(info.beta),\
	   cbeta=cos(info.beta);//sine and cosine of beta
double X[3];				//vector from bound vortex of elementary wing mid point
							//to its edge; x1=xo-x and x2=xo+x
double S[3];				//vector of bound vortex
double UxS;					//abs. value of U x S
double eN[3];				//normal force direction
double eL[3],eS[3];			//lift force and side force direction
double w1[3],wo[3],w2[3];	//ind. vel. at bound vortex sides and center
double gamma1,gammao,gamma2;//vorticity at bound vortex sides and center
double R1[3],Ro[3],R2[3]; 	//resultant ind. force at bound vortex sides and center
double R[3];				//resultant ind. force/density of element l
double N_free;				//magnitude free stream norm. forces/density
double tempA[3], tempS;

//loop over number of panels
for (i=0;i<info.nopanel;i++)
{
  //loop over number of chordwise elements 'info.m'
//removed GB 2-9-20  for (j=0;j<info.m;j++)
  for (j=0;j<panelPtr[i].m;j++)
  {
	 eta  =	elementPtr[l].eta;
	 eta8=eta*0.8;
	 //vector of bound vortex, KHH eq. 66
	 S[0] =	tan(elementPtr[l].phi);
	 S[1] = cos(elementPtr[l].nu);
	 S[2] = sin(elementPtr[l].nu);

	 //loop over number of spanwise elements 'n'
	 for (k=0;k<panelPtr[i].n;k++)
	 {
	/*##################################################
	//#############taken out 6/20/03 G.B. ##############
	//##################################################
	///////////////////////////////////////////////////
		 //#non-small angles, non-linear theory
	//////////////////////////////////////////////////
		 if (info.linear==0)
		 {

			 //normal force direction
			 cross(elementPtr[l].u,S,tempA);			//#	UxS
			 UxS=norm2(tempA);							//	|UxS|
			 scalar(tempA,1/UxS,eN);					//	eN=(UxS)/|UxS|

			 //the lift direction  eL=Ux[0,1,0]/|Ux[0,1,0]|
			 tempS = sqrt(elementPtr[l].u[0]*elementPtr[l].u[0]\
			 			 +elementPtr[l].u[2]*elementPtr[l].u[2]);
			 eL[0] = -elementPtr[l].u[2]/tempS;
			 eL[1] =  0;
			 eL[2] =  elementPtr[l].u[0]/tempS;

			 //the side force direction eS=UxeL/|UxeL|
			 cross(eL,elementPtr[l].u,tempA);
			 tempS=1/norm2(tempA);
			 scalar(tempA,tempS,eS);
		 }
		 else
	//////////////////////////////////////////////////
		 //small angle approach, linear theory
	///////////////////////////////////////////////////
	//##################################################//
	//##################################################//
	//##################################################*/
		 {

			 //#small angle approach as done by KHHH
			 tempS=norm2(elementPtr[l].u);
			 tempA[0] =  tempS * sbeta*S[2];			//# UxS
			 tempA[1] = -tempS * cbeta*S[2];			//#
			 tempA[2] =  tempS *(cbeta*S[1]-sbeta*S[0]);//#
			 UxS=norm2(tempA);							//	|UxS|
			 scalar(tempA,1/UxS,eN);					//	eN=(UxS)/|UxS|

			 //lift force direction
			 eL[0] =  0;
			 eL[1] =  0;
			 eL[2] =  1;

			 //side force direction
			 eS[0] =  0;
			 eS[1] =  1;
			 eS[2] =  0;
//#			 printf("small angle approximation for lift computation!\n");
		 }

//#printf("U = %lf\t%lf\t%lf\n",elementPtr[l].u[0],elementPtr[l].u[1],elementPtr[l].u[2]);//#
//#printf("S = %lf\t%lf\t%lf\t%lf\n",S[0],S[1],S[2],norm2(S));//#
//#printf("eN= %lf\t%lf\t%lf\t%lf\n",eN[0],eN[1],eN[2],norm2(eN));//#
//#printf("eL= %lf\t%lf\t%lf\t%lf\n",eL[0],eL[1],eL[2],norm2(eL));//#
//#printf("eS= %lf\t%lf\t%lf\t%lf\n",eS[0],eS[1],eS[2],norm2(eS));//#
//#printf("UxS= %lf\n",(UxS));//#
//#printf("phi= %lf\t%lf\tnu= %lf\n",(elementPtr[l].phi*RtD),elementPtr[l].phi*RtD/cos(elementPtr[l].nu),(elementPtr[l].nu*RtD));//#


		 A	  =	elementPtr[l].A;
		 B	  =	elementPtr[l].B;
		 C	  =	elementPtr[l].C;
//#vorticity at right and left edge, as well as center of elementary wing
//#printf("%lf\t%lf\t%lf\t\n",(A-B*eta+C*eta*eta)/2,A/2,(A+B*eta+C*eta*eta)/2);//#
//#printf("A = %lf\tC = %lf\teta = %lf\n",A,C,eta);//#

//*****************************************************************************
		 //computing magnitude of normal force/density due to free stream
//*****************************************************************************
		 N_free = (A*2*eta + C/3*2*eta*eta*eta)*UxS;
//#printf("N_free =%lf\t L_free =%lf\n",N_free,2*N_free*sqrt(eN[0]*eN[0]+eN[2]*eN[2]));//#

//*****************************************************************************
		 //computing the induced force/density
//*****************************************************************************
		 //computing the ind. velocity at left (1) edge of bound vortex
		 //vector from mid point of elementary wing bound vortex
		 //to edge (1) of bound vortex; x1=xo-x and x2=xo+x
		 //due to singular behavior of velocity at element edge
		 //velocity is computed 0.1eta away from edge (hence factor 0.8)
		 scalar(S,-eta8,X);
		 vsum(elementPtr[l].xo,X,tempA);
		 Induced_Velocity(elementPtr,info,tempA,w1);
		 					//subroutine in induced_velocity.cpp
//printf("x1 %2.5lf %2.5lf %2.5lf\n",tempA[0],tempA[1],tempA[2]);//#

	 	 //computing the ind. velocity at center (0) of bound vortex
		 Induced_Velocity(elementPtr,info,elementPtr[l].xo,wo);
		 					//subroutine in induced_velocity.cpp

		//saving induced velocity at center of bound vortex //added 7-20-05 G.B.
		elementPtr[l].uind[0] = wo[0];
		elementPtr[l].uind[1] = wo[1];
		elementPtr[l].uind[2] = wo[2];

		 //computing the ind. velocity at right (2) edge of bound vortex
		 //vector from mid point of elementary wing bound vortex
		 //to edge (2) of bound vortex; x1=xo-x and x2=xo+x
		 //due to singular behavior of velocity at element edge
		 //velocity is computed 0.1eta away from edge (hence factor 0.8)
		 scalar(S,eta8,X);
		 vsum(elementPtr[l].xo,X,tempA);
		 Induced_Velocity(elementPtr,info,tempA,w2);
		 					//subroutine in induced_velocity.cpp
//printf("x1 %2.5lf %2.5lf %2.5lf\n",tempA[0],tempA[1],tempA[2]);//#

//printf("w1= %lf\t%lf\t%lf\n",w1[0],w1[1],w1[2]);//#
//printf("wo= %lf\t%lf\t%lf\n",wo[0],wo[1],wo[2]);//#
//printf("w2= %lf\t%lf\t%lf\n",w2[0],w2[1],w2[2]);//#
//printf("w1= %.3lf\tw0= %.3lf\tw2= %.3lf\n",w1[2],wo[2],w2[2]);//#

	  //Integration of induced forces with Simpson's Rule
	  //Integration requires overhanging edges!!
	  //See also KHH linees 2953 - 2967, A23SIM

		 //Kutta-Joukowski at left (1) edge
		 cross(w1,S,tempA);				// w1xS
		 gamma1  = A-B*eta8+C*eta8*eta8;//gamma1
		 scalar(tempA,gamma1,R1);

		 //Kutta-Joukowski at center
		 cross(wo,S,tempA);				// woxS
		 gammao  =	A;
		 scalar(tempA,gammao,Ro);

 		 //Kutta-Joukowski at right (2) edge
 		 cross(w2,S,tempA);				// w2xS
		 gamma2  =	A+B*eta8+C*eta8*eta8;
		 scalar(tempA,gamma2,R2);

//#printf("gamma %lf\t%lf\t%lf\n",gamma1,gammao,gamma2);//#
//#printf("R1     %lf\t%lf\t%lf\n",R1[0],R1[1],R1[2]);//#
//#printf("Ro     %lf\t%lf\t%lf\n",Ro[0],Ro[1],Ro[2]);//#
//#printf("R2     %lf\t%lf\t%lf\n",R2[0],R2[1],R2[2]);//#

		 //The resultierende induced force of element l is
		 //determined by numerically integrating forces acros element
		 //using Simpson's Rule with overhaning parts
		 R[0]  = (R1[0]+4*Ro[0]+R2[0])*eta8/3;			//Rx
		 R[1]  = (R1[1]+4*Ro[1]+R2[1])*eta8/3;			//Ry
		 R[2]  = (R1[2]+4*Ro[2]+R2[2])*eta8/3;			//Rz

		 //plus overhanging parts
		 R[0] += (7*R1[0]-8*Ro[0]+7*R2[0])*(eta-eta8)/3;//Rx
		 R[1] += (7*R1[1]-8*Ro[1]+7*R2[1])*(eta-eta8)/3;//Ry
		 R[2] += (7*R1[2]-8*Ro[2]+7*R2[2])*(eta-eta8)/3;//Rz
//#printf("R %lf\t%lf\t%lf\n",R[0],R[1],R[2]);//#
//*****************************************************************************
		 //the LIFT FORCE/density is the normal force in the x-z plane or
//*****************************************************************************
		 elementPtr[l].N_free[0] = N_free * sqrt(eN[0]*eN[0]+eN[2]*eN[2]);
		 //neg. if resultend is downward
		 if (eN[2]<0)  elementPtr[l].N_free[0] *= -1;
		 elementPtr[l].N_ind[0] = dot(R,eL);
//*****************************************************************************
	  	 //the SIDE FORCE/density is the force in y-direction or N*eN[y]
//*****************************************************************************
		 elementPtr[l].N_free[1] = N_free * eN[1];
		 elementPtr[l].N_ind[1]  = dot(R,eS);


//#printf("CLfree =%lf\tCYfree =%lf\n",2*elementPtr[l].N_free[0],elementPtr[l].N_free[1]);//#
//#printf("CLi =%lf\tCYi =%lf\n",2*elementPtr[l].N_ind[0],elementPtr[l].N_ind[1]);//#

		 l++; //increment elementary wing index l=0..(noelement-1)
	 }	//End loop over k
  }	//End loop over j
} //End loop over i
}
//===================================================================//
		//END FUNCTION Elementary_Wing_Normal_Forces
//===================================================================//

//===================================================================//
		//START FUNCTION Wing_Normal_Forces
//===================================================================//
void Wing_Normal_Forces(const  PANEL* panelPtr, const GENERAL info,\
						const BOUND_VORTEX* elementPtr, \
						double Nt_free[2], double Nt_ind[2],\
						double &CL,double &CLi,double &CY,double &CYi)
{
	//input:
	// panel		-panel info, especially number of chord and span devisions
	// info			-general information on case
	// elementPtr	-pointer to elementary-wing information
	//
	//ouput:
	// as part of elementPtr.:
	// N_free	-total lift and side forces/density due to free stream flow
	// N_ind	-total lift and side forces/density due to induced velocities
	//
	// Nt_free	-total lift and side forces/density due to free stream flow
	// Nt_ind	-total lift and side forces/density due to induced velocities
	// CL		-total lift coefficient
	// CLi		-total induced lift coefficient
	// CY		-total side-force coefficient
	// CYi		-total induced side-force coefficient
	// CDi		-total induced drag computed at lifting lines

  int l=0;								 //counter
  double q=1/(0.5*info.Uinf*info.Uinf*info.S); //ref. area* dyn. pressure/density

  Nt_free[0]=0;
  Nt_free[1]=0;
  Nt_ind[0]=0;
  Nt_ind[1]=0;

  //loop over number of panels
  for (l=0;l<info.noelement;l++)
  {
	  //adding the normal forces/density of all elementary wings
	  Nt_free[0]	+=  elementPtr[l].N_free[0];
	  Nt_free[1]	+=  elementPtr[l].N_free[1];

	  Nt_ind[0]	+=  elementPtr[l].N_ind[0];
	  Nt_ind[1]	+=  elementPtr[l].N_ind[1];

//#printf("NtX  =%lf\t NtZ  =%lf\n",Nt_free[0],Nt_free[1]);//#
//#printf("NXi =%lf\t NZi =%lf\n",elementPtr[l].N_ind[0],elementPtr[l].N_ind[1]);//#
//#printf("NtXi =%lf\t NtZi =%lf\n",Nt_ind[0],Nt_ind[1]);//#
  }

  if (info.sym==1 && info.beta == 0)
  {	//twice the lift force if symmetric geometry
	  Nt_free[0]*=2;

	  Nt_ind[0]	*=2;
//#printf("NtX  =%lf\t NtZ  =%lf\n",Nt_free[0],Nt_free[1]);//#
  }
  //total lift and side force coefficients
  CL = (Nt_free[0]+Nt_ind[0])*q;
  CY = (Nt_free[1]+Nt_ind[1])*q;

  //total induced lift and side force coefficients
  CLi = Nt_ind[0]*q;
  CYi = Nt_ind[1]*q;

//#printf("CL=%lf\tCLi=%lf\tCY=%lf\tCYi=%lf\n",CL,CLi,CY,CYi);//#

}
//===================================================================//
		//END FUNCTION Wing_Normal_Forces
//===================================================================//
