//computes total normal force acting on wing
void DVE_Wing_Normal_Forces(const GENERAL,const PANEL *,\
                            const DVE *,const int,double **,\
                            double *, double **,\
                            double [2],double [2],\
							double &,double &,double &,double &);
//computes normal forces of surface DVE
void Surface_DVE_Normal_Forces(const GENERAL,const PANEL *,const int,\
							   DVE **,DVE *,double **);


//===================================================================//
		//START FUNCTION DVE_Wing_Normal_Forces
//===================================================================//
void DVE_Wing_Normal_Forces(const GENERAL info,const PANEL *panelPtr,\
                            const DVE *surfacePtr,const int time,\
                            double **N_force,double *D_force,\
                            double **Span_force,
							double Nt_free[2], double Nt_ind[2],\
							double &CL,double &CLi,double &CY,double &CYi)
{
//this routine adds up the DVE's normal forces in order to compute the
//total wing normal forces/density and coefficients based on free stream
//input:
// info		general information on case
// panelPtr panel informatin
// surfacePtr   surface DVE information
// time     current timestep
// N_force	normal forces/density of each surface DVE, second index is:
//			[0]: free stream lift, [1]: induced lift,
//			[2]: free stream side, [3]: induced side force/density
//          [4]: free stream normal, [5]: induced normal force/density
//          [6,7,8]: eN_x, eN_y, eN_z in global ref. frame
// D_force  drag force/density across span
//
//ouput:
//
// Span_force   force/density vector across span in wind axis system
// Nt_free	total lift, side and normal forces/density due to free stream flow
// Nt_ind	total lift, side and normalforces/density due to induced velocities
// CL		total lift coefficient
// CLi		total induced lift coefficient
// CY		total side-force coefficient
// CYi		total induced side-force coefficient

    int l,m,n,panel,span,index;				//counter
    double eN[3],eD[3];         //normal force and drag directions
    double omega,cosOm,sinOm;               //turn angle
    double tempS,tempA[3];
    double q = 1/(0.5*info.Uinf*info.Uinf*info.S); 	//1/(ref. area* dyn. pressure/density)
    double q_local;         //1/(ref. area dyn. pressure) of DVE
//===================================================================//
                        //span forces
//===================================================================//
     if(info.flagCIRC) //total angle turned if turning flight
     {
         omega=(info.gradient*info.deltime*(time-1)); //turn angle
         cosOm = cos(omega);
         sinOm = sin(omega);
     }

    span = 0;       //initializing span index
    index =0;       //initialzing surface DVE index
    //loop over panels
    for(panel=0;panel<info.nopanel;panel++)
    {
        //loop over panel span (along leading edge indices)
        for(n=panelPtr[panel].LE1;n<=panelPtr[panel].LE2;n++)
        {
            index = n; //setting index to first chordwise DVE of span location
            
            //initializing
            Span_force[span][0]=0; Span_force[span][1]=0; Span_force[span][2]=0;

            //loop over chord of panel
            for(m=0;m<panelPtr[panel].m;m++)
            {
                //normal force direction
                eN[0] = N_force[index][6];
                eN[1] = N_force[index][7];
                eN[2] = N_force[index][8];
 
                //adding normal forces
                tempS = (N_force[index][4]+N_force[index][5]); //total normal force
                Span_force[span][0] += eN[0]*tempS;
                Span_force[span][1] += eN[1]*tempS;
                Span_force[span][2] += eN[2]*tempS;

                index +=panelPtr[panel].n; //surfaceDVE indexing down the chord
            }//next chord element; 'index' should be value of surfaceDVE at trailing edge

            index = index - panelPtr[panel].n; //adjusting surfaceDVE
 
//################################################
            // drag force direction
            if(!info.flagCIRC)
              {
                  //drag force direction
                  eD[0] = surfacePtr[index].U[0];
                  eD[1] = surfacePtr[index].U[1];
                  eD[2] = surfacePtr[index].U[2];
              }
              else
              {
                 // Added by D.F.B. 03-2020 because of circling flight
                 // If there is circling flight, set the DVE drag direction to the
                 // velocitiy direction at the TE
                 tempS = 1/norm2(surfacePtr[span].uTE[0]);
                 eD[0] = surfacePtr[index].uTE[0][0]*tempS;
                 eD[1] = surfacePtr[index].uTE[0][1]*tempS;
                 eD[2] = surfacePtr[index].uTE[0][2]*tempS;
              }

            //adding drag
            Span_force[span][0] += eD[0]*D_force[span];
            Span_force[span][1] += eD[1]*D_force[span];
            Span_force[span][2] += eD[2]*D_force[span];
            
//################################################
            //rotate into wind axis system
            if(info.flagCIRC) //turning flight -> rotate vector to wind-axis frame
            {   //rotating by omega
                tempA[0] = Span_force[span][0]*cosOm - Span_force[span][1]*sinOm;
                tempA[1] = Span_force[span][0]*sinOm + Span_force[span][1]*cosOm;
                tempA[2] = Span_force[span][2];
                //reassigning
                Span_force[span][0]=tempA[0];
                Span_force[span][1]=tempA[1];
                Span_force[span][2]=tempA[2];
             }
    //NOTE! if body-reference frame required, rotation by alpha needed
//################################################
//printf("\n%d Spanforce %lf  %lf  %lf ",\
//       span,tempA[0],tempA[1],tempA[2]);
//printf("\n u   %lf  %lf  %lf ",\
       surfacePtr[index].u[0],surfacePtr[index].u[1],surfacePtr[index].u[2]);

            span++; // increase to next span index
        }//loop over span (n) of panel
    } //next panel
    
//===================================================================//
                      //section loads
//===================================================================//
    //compute section force coefficients
    span=0;
    for(panel=0;panel<info.nopanel;panel++)
    {
       //loop over panel span (along leading edge indices)
       for(n=panelPtr[panel].LE1;n<=panelPtr[panel].LE2;n++)
       {
           //1/(0.5 U^2 S) of chordwise strip
           tempS = 2/(dot(surfacePtr[n].u,surfacePtr[n].u)*\
                      surfacePtr[n].S*panelPtr[panel].m);
           Cf[span][0] = Span_force[span][0]*tempS;
           Cf[span][1] = Span_force[span][1]*tempS;
           Cf[span][2] = Span_force[span][2]*tempS;
           
           span++; //next span section
       }
    }
    
//===================================================================//
                    //total forces
//===================================================================//
    //computing CFX,CFY and CFZ
    CF[0] = 0; CF[1] = 0; CF[2] = 0; //initialzing aircraft coefficients

    for(span=0;span<info.nospanelement;span++)
    {
      CF[0] += Span_force[span][0];
      CF[1] += Span_force[span][1];
      CF[2] += Span_force[span][2];
    }
    scalar(CF,q,CF);

    //initializing
    Nt_free[0]=0;
    Nt_free[1]=0;
    Nt_free[2]=0;
    Nt_ind[0]=0;
    Nt_ind[1]=0;
    Nt_ind[2]=0;

    printf("\n");
    
    //loop over number of surfaceDVEs
    for (l=0;l<info.noelement;l++)
    {
        //adding the normal forces/density of all elementary wings
        Nt_free[0]	+=  N_force[l][0];
        Nt_free[1]	+=  N_force[l][2];
        Nt_free[2]  +=  N_force[l][4];

        Nt_ind[0]	+=  N_force[l][1];
        Nt_ind[1]	+=  N_force[l][3];
        Nt_ind[2]   +=  N_force[l][5];
printf("DVE %d  Nx  =%lf\t Ntxind  =%lf  in lift_force\n",l,N_force[l][0],N_force[l][1]);//#
    //#printf("NtX  =%lf\t NtZ  =%lf\n",Nt_free[0],Nt_free[1]);//#
    //#printf("NtXi =%lf\t NtZi =%lf\n",Nt_ind[0],Nt_ind[1]);//#
    }

    
    if (info.sym==1 && info.beta == 0)
    {	//twice the force if symmetric geometry
        Nt_free[0]*=2;
        Nt_free[2]*=2;

        Nt_ind[0]	*=2;
        Nt_ind[2]	*=2;
    //#printf("NtX  =%lf\t NtZ  =%lf\n",Nt_free[0],Nt_free[1]);//#
    }
    //total lift and side force coefficients
    CL = (Nt_free[0]+Nt_ind[0])*q;
    CY = (Nt_free[1]+Nt_ind[1])*q;

    //total induced lift and side force coefficients
    CLi = Nt_ind[0]*q;
    CYi = Nt_ind[1]*q;

    printf("\nCL=%lf\tCLi=%lf\tCY=%lf\tCYi=%lf \nCFX %lf CFY %lf CFZ %lf  |CF| %lf\n",\
           CL,CLi,CY,CYi,CF[0],CF[1],CF[2],norm2(CF));//#

    //===================================================================//
    //===================================================================//
//################
    
    // N_force    normal forces/density of each surface DVE, second index is:
    //            [0]: free stream lift, [1]: induced lift,
    //            [2]: free stream side, [3]: induced side force/density
    //          [4]: free stream normal, [5]: induced normal force/density
    //          [6,7,8]: eN_x, eN_y, eN_z in global ref. frame
    // D_force  drag force/density across span
    //
    //ouput:
    //
    // Span_force   force/density vector across span in wind axis system
    // Nt_free    total lift, side and normal forces/density due to free stream flow
    // Nt_ind    total lift, side and normalforces/density due to induced velocities
    
    
   // for (l=0;l<info.noelement;l++)
   // {
   //     printf
//}

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
//                  [6,7,8]: eN_x, eN_y, eN_z in global ref. frame
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
			//vector along the mid chord of DVE, for averaging of DVEs aft of LE
			tempA[0]=tan(surfacePtr[l].phi0); tempA[1]=1; tempA[2]=0;
			//transforming into local reference frame
			Star_Glob(tempA,surfacePtr[l].nu,surfacePtr[l].epsilon,\
												surfacePtr[l].psi,S);
								//function in ref_frame_transform.h
		}

		// S calculation works for circling flight D.F.B. 03-2020
		eta  = surfacePtr[l].eta;
		eta8 = eta*0.8;


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
         
        //correction of vector along LE for DVEs aft of LE row
        if(j>0)
        {
             //vector along the bound vortex along LE
             tempA[0]=tan(surfacePtr[l].phiLE); tempA[1]=1; tempA[2]=0;
             //transforming into local reference frame
             Star_Glob(tempA,surfacePtr[l].nu,surfacePtr[l].epsilon,\
                                                 surfacePtr[l].psi,S);
                                 //function in ref_frame_transform.h
        }

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
    //computing normal force vector
//*****************************************************************************

        //#non-small angles
        //normal force direction
        cross(surfacePtr[l].u,S,tempA);            //#    UxS
        UxS=norm2(tempA);                        //    |UxS|
        scalar(tempA,1/UxS,eN);                    //    eN=(UxS)/|UxS|

        //*** Calculated lift direction based on the freestream velocity direction
        //        Changed by D.F.B. 03-2020
        //vector along the bound vortex along LE
        tempA[0]=0; tempA[1]=1; tempA[2]=0;
        //transforming into local reference frame
        Star_Glob(tempA,0,surfacePtr[l].epsilon,surfacePtr[l].psi,spandir);

        cross(surfacePtr[l].u,spandir,tempA);        //#    UxS
        Uxspandir=norm2(tempA);                        //    |UxS|
        scalar(tempA,1/Uxspandir,eL);                //    eL=(UxS)/|UxS|

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
         
/* ************************ Quiver output of lift vector *************************
         if(i==0 & j==0 & k ==0){CreateQuiverFile(surfacePtr[l].xo, eL,0);}
         else{CreateQuiverFile(surfacePtr[l].xo, eL,1);}

         CreateQuiverFile(surfacePtr[l].xo, eL,1);
// ************************* Quiver output of lift vector *************************/

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


//*****************************************************************************
     //computing magnitude of normal force/density due to free stream
//*****************************************************************************
        N_free = (A*2*eta + C/3*2*eta*eta*eta)*UxS;
     // N_free calculation works for circling flight D.F.B. 03-2020

//*****************************************************************************
		 //the NORMAL FORCE/density
//*****************************************************************************
		//free stream lift
		N_force[l][4] = N_free;
		N_force[l][5] = dot(R,eN);			//induced normal force
         
         N_force[l][6] = eN[0];         //saving normal force vector
         N_force[l][7] = eN[1];
         N_force[l][8] = eN[2];

//*****************************************************************************
		 //the LIFT FORCE/density is the normal force in the x-z plane or
//*****************************************************************************
		//free stream lift
		N_force[l][0] = N_free * sqrt(eN[0]*eN[0]+eN[2]*eN[2]);
		if (eN[2]<0)  N_force[l][0] *= -1;	//neg. if resultend is downward

		N_force[l][1] = dot(R,eL);			//induced lift
//#printf("N_free =%lf\t L_free =%lf\n",N_free,2*N_free*sqrt(eN[0]*eN[0]+eN[2]*eN[2]));//#

//* ************************ Quiver output of normal vector *************************
         scalar(eN,N_free/info.Uinf,tempA);
		if(i==0 && j==0 && k ==0)
                CreateQuiverFile(surfacePtr[l].xo, tempA,0);
		else    CreateQuiverFile(surfacePtr[l].xo, tempA,1);
//************************* Quiver output of normal vector *************************/

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
