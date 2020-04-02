//computes induced drag at trailing edge of DVE wing
double Induced_DVE_Drag(const GENERAL,const PANEL*,DVE*,DVE**,\
												const int,double*);
//computes section drag of surface DVE
double SectionDrag(double ***, double,double &,int,double &,const int);

//===================================================================//
	//START Induced_DVE_Drag computation - Drag along trailing edge
//===================================================================//
	//computes induced drag at trailing edge
double Induced_DVE_Drag(const GENERAL info,const PANEL* panelPtr,\
						DVE* surfacePtr,DVE** wakePtr,\
						const int rightnow,double* D_force)
{
//This function is the DVE expansion of the function Induced_Eppler_Drag
//that computes the induced drag at the trailing edge, where
//the spanwise bound vorticity has been collapsed into a single vortex.
//The method is discussed more thoroughly in:
//Eppler and Schmid-Goeller, "A Method to Calculate the Influence of
//Vortex Roll-Up on the Induced Drag of Wings," Finite Approximations
//in Fluid Mechanics II, Hirsche, E. H. (ed.), Notes on Numerical Fluid
//Mechanics, Volume 25, Friedr. Vieweg und Sohn, Braunschweig, 1990
//Kutta-Joukowsky is being applied to the trailing edge at three points
//along each trailing edge element.  Similarly to the lift computation,
//Simpson's Rule is used to compute the total drag of each element.

//Function computes forces/density for each spanwise section of the
//wing by applying Kutta-Joukowski's theorem to either of its edges and
//its center in order to get the local forces. The total spanwise section
//forces are determined by integrating with Simpson's rule (with overhang)
//Note: the velocities are not computed directly at the edges, but
//(1-delta)/2 of the elementary span further towards the center.
//Otherwise, the computed induced velocity would become singular at that
//point due to the next neighboring elementary wing influence.
//
//input:
//	info		- general information
//	panelPtr	- information on panels
//	surfacePtr	- information on surface DVEs
//	wakePtr		- information on wake DVEs
//	rightnow	- current time step
//
//output:
//	CDi			- total drag coefficient
//  D_force 	- local drag force/density along span

int panel,p,span,s,time,k,wing;
int	index;					//span index of surface DVEs along trailing edge
int i;						//span index of wake DVEs
double A, B, C;				//vorticity distribution coefficient along t.e.
double eta, eta8;			//half span of elementary wing l, 90% value of eta
double eD[3],eL[3],eS[3];	//drag, lift, side force direction
double xsiTE,phiTE;			//dist. most aft surf, DVE ref.pt to TE, TE sweep
double S[3];				//trailing edge vector
double X[3][3];				//points along trailing edge element,
							//left X[1], center X[0], right X[2]
double delX[3],Xstar[3];	//delta for corection, corrected X
double w[3],w_ind[3][3];	//delta and total ind. vel. at t.e. edges & center
double gamma1,gammao,gamma2;//vorticity at trailing edge edges and center
double R1[3],Ro[3],R2[3]; 	//res. ind. force at trail. edge edges and center
double R[3];				//resultant ind. force/density of element i
double CDi = 0,CLi=0,CYi=0;	//total induced drag at trailind edge
double tempA[3],tempB[3],tempS;
double tempC[3][3];
double tempProj[3],tempTE[3];
int type;					//type of wake DVE
DVE tempDVE;				//temporary DVE


//#############################################################################
//							FORCE LOOP - START
//#############################################################################
//
//in this loop the induced forces at the trailing edge are computed for three
//points along the trailing egd of the wing in the region of the surface DVE
//with the index 'index'.  The forces are integrated over the surface DVE's
//span in order to get the induced drag contribution of that span location.

	i = 0;  //initializing wake DVE index
	wing =0; //initializing wing index

	//loop over panels
	for (panel=0;panel<info.nopanel;panel++)
	//loop over trailing edge elements of current panel
	for (index = panelPtr[panel].TE1; index <= panelPtr[panel].TE2; index++)
	{
		//increase wing index to next wing
        if(index>info.dve2[wing]) wing++;

        if(!info.flagCIRC)
        {
        	//###
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
			eD[0] = surfacePtr[index].uTE[0][0];
			eD[1] = surfacePtr[index].uTE[0][1];
			eD[2] = surfacePtr[index].uTE[0][2];
		}


	//#########################################################################
		/* The lift and side for calculation have been romoved by D.F.B. 03-2020
		// Calcs removed because the are incorrect but not even used in this function
		// If need in future, see the fcn Surface_DVE_Normal_Forces in lift_force.cpp
		// for the correct calcs.
		//
		//the lift direction  eL={U x [0,1,0]}/|U x [0,1,0]|
		tempS = 1/sqrt(surfacePtr[index].U[0]*surfacePtr[index].U[0]\
		 				+surfacePtr[index].U[2]*surfacePtr[index].U[2]);
		eL[0] = -surfacePtr[index].U[2]*tempS;
		eL[1] =  0;
		eL[2] =  surfacePtr[index].U[0]*tempS;

		//the side force direction eS=UxeL/|UxeL|
		cross(eL,surfacePtr[index].U,tempA);
		tempS=1/norm2(tempA);
		scalar(tempA,tempS,eS);
		*/
	//#########################################################################
		A =	surfacePtr[index].A;
		B =	surfacePtr[index].B;
		C =	surfacePtr[index].C;

	//#########################################################################
		//Computing the three points along the unswept trailing edge,
		//at which Kutta-Joukowsky is applied.

		//The wing-trailing edge is located at the TRAILING EDGE of most aft
		//surface DVE
		xsiTE = surfacePtr[index].xsi;

		//The trailing-edge sweep is phiTE of the DVE
		phiTE = surfacePtr[index].phiTE;

		//the left and right points are 20% of half span away from edge
		//in order to stay away from the singularity along the edge of
		//the DVE
		eta  =	surfacePtr[index].eta;
		eta8=eta*.8;	//0.8 as done for lift computation,

		//X1:
		Edge_Point(surfacePtr[index].xo,surfacePtr[index].nu,\
				   surfacePtr[index].epsilon,surfacePtr[index].psi,\
				   phiTE,-eta8,xsiTE,X[1]);
				   						//Subroutine in wake_geometry.cpp
		//X2:
		Edge_Point(surfacePtr[index].xo,surfacePtr[index].nu,\
				   surfacePtr[index].epsilon,surfacePtr[index].psi,\
				   phiTE,eta8,xsiTE,X[2]);
				   						//Subroutine in wake_geometry.cpp
		//X0 = (X1+X2)/2
		vsum(X[1],X[2],tempA);		scalar(tempA,0.5,X[0]);

		//if(index == panelPtr[panel].TE1){CreateQuiverFile(X[1], eD,0);}
		//else{CreateQuiverFile(X[1], eD,1);}


//printf("\nx0 %2.3lf  %2.3lf  %2.3lf ",X[0][0],X[0][1],X[0][2]);  //#
//printf("x1 %2.3lf  %2.3lf  %2.3lf ",X[1][0],X[1][1],X[1][2]);  //#
//printf("x2 %2.3lf  %2.3lf  %2.3lf ",X[2][0],X[2][1],X[2][2]);  //#

	//#########################################################################
		//computing the normalized vector along the trailing edge
		S[0] = X[2][0] - X[1][0];
		S[1] = X[2][1] - X[1][1];
		S[2] = X[2][2] - X[1][2];
		tempS= 0.5/eta8;//1/norm2(S) //I don't why, but it works G.B. 5/30/05
		scalar(S,tempS,S);

//printf("\nS %lf %lf %lf  TE: %lf  %lf  %lf  ",S[0],S[1],S[2],\
//	surfacePtr[index].TEvc[0],surfacePtr[index].TEvc[1],surfacePtr[index].TEvc[2]);
	//#########################################################################
		//initializing induced velocities of current span location
		w_ind[1][0] = 0;  w_ind[1][1] = 0;  w_ind[1][2] = 0;
		w_ind[0][0] = 0;  w_ind[0][1] = 0;  w_ind[0][2] = 0;
		w_ind[2][0] = 0;  w_ind[2][1] = 0;  w_ind[2][2] = 0;

  //###########################################################################
  //						INDUCED VELOCITY LOOP - START
  //###########################################################################
  //
  //In this loop, the velocity are computed that are induced by the entire flow
  //field at the three points along the wing-trailing edge in the region of the
  //surface DVE 'index'.  This is done for each point at a time.

		//loop over the three points of trailing edge
		for (k=0;k<3;k++)
		{
		  span=0;	//initializing span index for wake DVEs

		  //loop over panels
		  for (p=0;p<info.nopanel;p++)
		  //loop over trailing edge elements of current panel
		  for (s = panelPtr[p].TE1; s <= panelPtr[p].TE2; s++)
		  {
             if(s>=info.dve1[wing] && s<=info.dve2[wing])
			{
			//DVE 's' (the inducer) and DVE 'index' (the induced one) are
			//of the same wing


		  	//################################################################
		  	//New method of moving the TE points of the index (induced DVE) points. 
				//17 Oct 2014. Bill B
			//This method moves the points in the freestream direction into the plane 
				//passing through the TE of the S (inducer) having freestream direction 
				//as the normal
			
				//determine vector on inducer from control point to TE (call it tempB)
				tempB[0] = surfacePtr[s].xsi;
				tempB[1] = 0;
				tempB[2] = 0;

				//make temp B global, call it tempA
				Star_Glob(tempB,surfacePtr[s].nu,surfacePtr[s].epsilon,\
													  surfacePtr[s].psi,tempA);

				//add tempA to location of control point, so we are left with location 
				//of the inducer's (S) TE in global coords. call it tempTE

				tempTE[0]=surfacePtr[s].xo[0]+tempA[0];
				tempTE[1]=surfacePtr[s].xo[1]+tempA[1];
				tempTE[2]=surfacePtr[s].xo[2]+tempA[2];
				
				//vector from TE of S(inducer) to point k on TE of index (induced). call it delX
				delX[0] = X[k][0] - tempTE[0];
				delX[1] = X[k][1] - tempTE[1];
				delX[2] = X[k][2] - tempTE[2];


				if(info.flagCIRC){ //Added by D.F.B. 03-2020 for circling flight
					// NOTE: Move *induced point* in its own velocity direction
					tempS=dot(delX,surfacePtr[index].uTE[k]);
					scalar(surfacePtr[index].uTE[k],tempS,tempB);
				} else{ // if not circling flight, move pts in freestream direction (U)
					//delX projected into the freestream direction (magnitude)
					tempS=dot(delX,surfacePtr[index].U);
					
					//vector from TE of s(inducer) to TE of index (induced) projected into the freestream direction (with direction) to make tempB
					scalar(surfacePtr[index].U,tempS,tempB);
				}

				// (X[k] - tempB) should be global origin to new point of interest
				Xstar[0] = X[k][0] - tempB[0];
				Xstar[1] = X[k][1] - tempB[1];
				Xstar[2] = X[k][2] - tempB[2];

						
				}
			else
			{
			//DVE 's' (the inducer) and DVE 'index' (the induced one) are
			//of different wings

				Xstar[0] = X[k][0];
				Xstar[1] = X[k][1];
				Xstar[2] = X[k][2];
			}

		  //###################################################################
		  //computing the velocity induced by the freshly shed wake in the TE
		  //The very beginning of the shed strip of a vortex sheet belongs to
		  //the first (most recently pooped out) row of wake DVEs.

			//assigning temporary DVE that induces on trailing edge
			//as Schmid-Goeller discusses in his dissertation,
			//it has no sweep and belongs to a spanwise strip of wake
			//elements that starts at the point of interest

			tempDVE.xo[0] 	 = wakePtr[rightnow][span].xo[0];
			tempDVE.xo[1] 	 = wakePtr[rightnow][span].xo[1];
			tempDVE.xo[2] 	 = wakePtr[rightnow][span].xo[2];

			tempDVE.phiLE	 = 0;

			tempDVE.phiTE 	 = 0; //wakePtr[rightnow][span].phiTE;
			tempDVE.nu		 = wakePtr[rightnow][span].nu;
			tempDVE.epsilon  = wakePtr[rightnow][span].epsilon;
			tempDVE.psi		 = wakePtr[rightnow][span].psi;

			tempDVE.eta		 = wakePtr[rightnow][span].eta;
			tempDVE.xsi		 = wakePtr[rightnow][span].xsi;

			tempDVE.A		 = wakePtr[rightnow][span].A;
			tempDVE.B		 = wakePtr[rightnow][span].B;
			tempDVE.C		 = wakePtr[rightnow][span].C;

			tempDVE.singfct	 = 0; //surfacPtr[s].singfct;

	//		type = 4;  //vortex sheet reaches from 0.5xsi to xsi
			type = 1;  //DVE is only a vortex sheet
			
			//computes induced velocity in X[k] due to DVE tempDVE
			Single_DVE_Induced_Velocity(info,tempDVE,Xstar,w,type);
		 						//subroutine in induced_velocity.cpp

			w_ind[k][0] += w[0];  w_ind[k][1] += w[1];  w_ind[k][2] += w[2];

		  //###################################################################
		  //computing the induced velocities of the remaining wake DVEs of the
		  //current spanwise location

			//loop across wake elements
			//for(time=0;time<=rightnow;time++)
			for(time=0;time<rightnow;time++)
			{
			
			tempDVE.xo[0] 	 = wakePtr[time][span].xo[0];
			tempDVE.xo[1] 	 = wakePtr[time][span].xo[1];
			tempDVE.xo[2] 	 = wakePtr[time][span].xo[2];

			tempDVE.phiLE	 = 0;
			tempDVE.phiTE 	 = 0; //wakePtr[rightnow][span].phiTE;
			tempDVE.nu		 = wakePtr[time][span].nu;
			tempDVE.epsilon  = wakePtr[time][span].epsilon;
			tempDVE.psi		 = wakePtr[time][span].psi;

			tempDVE.eta		 = wakePtr[time][span].eta;
			tempDVE.xsi		 = wakePtr[time][span].xsi;

			tempDVE.A		 = wakePtr[time][span].A;
			tempDVE.B		 = wakePtr[time][span].B;
			tempDVE.C		 = wakePtr[time][span].C;

			tempDVE.singfct	 = wakePtr[time][span].singfct;

				if(time!=0) type = 1;  //DVE is only a vortex sheet
				else 		type = 3;//oldest wake is semi-infin. vort. sheet

				//computes induced velocity in X[k] due to remaining wake DVEs
				Single_DVE_Induced_Velocity\
							(info,tempDVE,Xstar,w,type);
		 						//subroutine in induced_velocity.cpp

				w_ind[k][0] += w[0];
				w_ind[k][1] += w[1];
				w_ind[k][2] += w[2];
			}//end loop over time, along a strip in wake

			span ++; //advancing span index of wake DVEs
		  }//end loop over panel span (s), panels (p)
		}//end loop over THE three points (k)
  //###########################################################################
  //						INDUCED VELOCITY LOOP - END
  //###########################################################################

//printf("\nw1 %lf   %lf   %lf  %lf\n",w_ind[1][0],w_ind[1][1],w_ind[1][2],norm2(w_ind[1]));  //#
//printf("w0 %lf   %lf   %lf  %lf\n",w_ind[0][0],w_ind[0][1],w_ind[0][2],norm2(w_ind[0]));  //#
//printf("w2 %lf   %lf   %lf  %lf\n",w_ind[2][0],w_ind[2][1],w_ind[2][2],norm2(w_ind[2]));  //#
  //###########################################################################
  //				AND NOW: THE FORCE INTEGRATION
  //###########################################################################

		//Integration of induced forces with Simpson's Rule
		//Integration requires overhanging edges!!
		//See also KHH linees 2953 - 2967, A23SIM

		//Kutta-Joukowski at left (1) edge
		cross(w_ind[1],S,tempA);			// w1xS
		gamma1  = A-B*eta8+C*eta8*eta8;		//gamma1
		scalar(tempA,gamma1,R1);

		//Kutta-Joukowski at center
		cross(w_ind[0],S,tempA);			// woxS
		gammao  = A;
		scalar(tempA,gammao,Ro);

 		//Kutta-Joukowski at right (2) edge
 		cross(w_ind[2],S,tempA);			// w2xS
		gamma2  = A+B*eta8+C*eta8*eta8;
		scalar(tempA,gamma2,R2);

//#printf("R1 %lf\t%lf\t%lf\n",R1[0],R1[1],R1[2]);//#
//#printf("Ro %lf\t%lf\t%lf\n",Ro[0],Ro[1],Ro[2]);//#
//#printf("R2 %lf\t%lf\t%lf\n",R2[0],R2[1],R2[2]);//#

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
//#printf("Ro %lf\t%lf\t%lf\n",R[0],R[1],R[2]);//#

	//#########################################################################
		//the DRAG FORCE/density is the induce force in eD direction
	//#########################################################################
		D_force[i] = dot(R,eD);

		//add all partial drag/lift/side values [force/density]
		CDi += D_force[i];

//    printf("\tm = %d  eD \t%2.8lf\t%2.8lf\t%2.8lf\n",\
//           i,eD[0],eD[1],eD[2]);  //###
        
		i++;	//advanicing span index of wake DVEs
	}//end loop over trailing edge surface DVEs and over panels

//#############################################################################
//							FORCE LOOP - END
//#############################################################################

	//non-dimensionalize
	tempS = 0.5*info.Uinf*info.Uinf*info.S;
	CDi /= tempS;

	if (info.sym==1)
	{
		CDi*=2;//sym. geometry and flow, twice the drag
	}

	///////////////////////////////////////////////////////////////////////
	//Drag-force contribution to total forces and moments of DVE
	//Moments are with respect to info.RefPt.
	//G.B. 11-24-06

	span = 0; //initializing D_force index
	//loop over panels
	for(panel=0;panel<info.nopanel;panel++)
	{
		//smallest index of panel-1
		index=panelPtr[panel].TE2-panelPtr[panel].n*panelPtr[panel].m;

	  	//loop over trailing edge elements of current panel
	  	for(k=panelPtr[panel].TE1;k<=panelPtr[panel].TE2;k++)
	  	{
			//drag per DVE is saved in tempS
			tempS = D_force[span]*info.density/panelPtr[panel].m;

			//loop over elements of one span location, from TE to LE
			for(i=k;i>index;i-=panelPtr[panel].n)
			{
				//the drag vector of each DVE is stored in R[3]
				R[0] = surfacePtr[i].U[0]*tempS;
				R[1] = surfacePtr[i].U[1]*tempS;
				R[2] = surfacePtr[i].U[2]*tempS;

				//surfacePtr[i].Force was initialized in
				//Surface_DVE_Normal_Forces in lift_force.cpp
				surfacePtr[i].Force[0] += R[0];
				surfacePtr[i].Force[1] += R[1];
				surfacePtr[i].Force[2] += R[2];

			}//end loop i, chordwise elements
			span++;  //increase span index

	  	}//end loop k, spanwise elements
	}//end loop panel, loop over panels
///////////////////////////////////////////////////////////////////////

//#	printf(" L=%lf  Y=%lf",CLi,CYi);
	return CDi;
}
//===================================================================//
	//END Induced_DVE_Drag computation - Drag along trailing edge
//===================================================================//

//================================================================================================================
// SectionDrag
//================================================================================================================
double SectionDrag(double **profiledata, double Re,double &cl,int airfoilCol,\
					double &cm,const int section)
{
	//input:
	//profiledata	airfoil data of airfoil in use
	//Re			Re# of interest
	//cl			cl of interest
	//airfoilCol	max number of rows in airfoil data file
	//section		index of wing section whose profile drag is computed
	//
	//
	//output
	//cm			interpolated section moment coefficient
	//cd			interpolated section drag coefficient

	//This routine finds the section drag coefficient of a particular cl and Re.  For that purpose
	//three interpolations are performend:
	//1. find cd of Re# below Re# of interest
	//2. find cd of Re# above Re# of intersest
	//3. interpolate between Re# results
	//
	//there are several special cases

	//HiRe,LoRe max and min Re# of airfoil data

	// Updated 2-14-20 by D.F.B. in Braunschweig, Germany
	// 	Now uses the new airfoil format and doesn't need number of rows

	int index,index1,index2;				//index of largest CL of Re#<RE, index of CL>cl of Re#<RE
											//index of CL>cl of Re#>RE
	double Re1,Re2;							//Re# above and below Re of interest
	//High and low Re in input file

	// Calcualte how many rows in the airfoil file	
	int rows = 0;
	for(int i = 0; i<airfoilCol;++i){
		if (profiledata[i][3] < DBL_EPS){break;}
		++rows;
	}

	double HiRe=profiledata[rows-1][3];
	double LoRe=profiledata[0][3];
	double m,cd,cd1,cd2;			//interpolation slope, overall and part interpolation results
	double cm1,cm2;				//moment coefficients for interpolation
	double tempS;
//###



//printf("rows %d lowRen %lf HighRe %lf\n",rows,LoRe,HiRe);
		if(LoRe<Re && HiRe>Re)  //Re# of interest falls into interval of available airfoil data
		{
			//find index of highest CL of Re# just below Re
			index=0;
			while(profiledata[index][3] < Re)
				{index++;}
			index--;  //adjust index

		//===================================
		//interpolation at lower Re

			//find index of CL of lower Re# just below cl of interest
			index1 = index;
			Re1 = profiledata[index1][3];	//lower Re#
			while(profiledata[index1][1]>cl && index1>=0 && profiledata[index1][3]==Re1)
				{index1--;}

			//cl falls outside of cl range (stalled or lower cl value) no cl-interpolation
			if(index1==index)
			{
				cd1=profiledata[index1][2]; 			//section stalled
				cm1=profiledata[index1][4];



				printf("1 STALL of section %d @ cl=%.3lf Re=%.0lf!\n",section,cl,Re);
				cl = .825*profiledata[index1][1];            //setting the cl to the 2D cl max of the airfoil data
				printf("1 STALL of section %d; cl from above set to cl max=%.01f\n", section, cl);
			}
			 else if(index1==-1 || profiledata[index1][3]!=Re1)
			{
				cd1=profiledata[index1+1][2];		//section cl less than what is listed
				cm1=profiledata[index1+1][4];		//section cl less than what is listed
				printf("1 Insufficient low cl airfoil data for section ");
				printf("%d cl=%.3lf and Re=%.0lf!\n",section,cl,Re);
			}

			else	//index1 is index of cl value below value of interest, thus interpolation between
					//index1 and index1+1
			{
				//(cl-cl2)/(cl1-cl2)
				tempS=(cl-profiledata[index1][1])\
						/(profiledata[index1][1]-profiledata[index1+1][1]);

				//interpolation of drag; slope (cd1-cd2)/(cl1-cl2)(cl-cl2)
				m =	(profiledata[index1][2]-profiledata[index1+1][2]);
				cd1=m*tempS + profiledata[index1][2];

				//interpolation of moment; slope (cm1-cm2)/(cl1-cl2)
				m =	(profiledata[index1][4]-profiledata[index1+1][4]);
				cm1=m*tempS + profiledata[index1][4];
			}

		//===================================
		//interpolation at higher Re

			//find idnex of CL of higher Re# just above cl of interest
			index2=index+1;
			Re2 = profiledata[index2][3];	//upper Re#
			while(profiledata[index2][1]<cl  && index2<rows && profiledata[index2][3]==Re2)
				{index2++;}

			//cl falls outside cl range no interpolation
			if(index2==rows || profiledata[index2][3]!=Re2)
			{
				cd2=profiledata[index2-1][2]; 			//section stalled
				cm2=profiledata[index2-1][4]; 			//section stalled
				printf("2 STALL of section %d @ cl=%.3lf Re=%.0lf!\n",section,cl,Re);
				cl = .825*profiledata[index2-1][1];            //setting the cl to the 2D cl max of the airfoil data
				printf("2 STALL of section %d; cl from above set to cl max=%.01f\n", section,  cl);
			}
			else if(index2 == (index+1))
			{
				cd2=profiledata[index2][2]; 			//section cl less than what is listed
				cm2=profiledata[index2][4]; 			//section cl less than what is listed
				printf("2 Insufficient low cl airfoil data for section ");
				printf("%d cl=%.3lf and Re=%.0lf!\n",section,cl,Re);
			}
			else	//index2 is index of cl value above value of interest, thus interpolation between
					//index2-1 and index2
			{
				//(cl-cl2)/(cl1-cl2)
				tempS=(cl-profiledata[index2-1][1])/\
					(profiledata[index2][1]-profiledata[index2-1][1]);

				//interpolation slope (cd1-cd2)/(cl1-cl2)
				m =	(profiledata[index2][2]-profiledata[index2-1][2]);
				cd2=m*tempS + profiledata[index2-1][2];

				//interpolation slope (cm1-cm2)/(cl1-cl2)
				m =	(profiledata[index2][4]-profiledata[index2-1][4]);
				cm2=m*tempS + profiledata[index2-1][4];
			}

			//===================================
			//Interpolation between Re#
			tempS=(Re-Re1)/(Re2-Re1);
			cd=(cd2-cd1)*tempS+cd1;  //<=========cd
			//Interpolation between Re#
			cm=(cm2-cm1)*tempS+cm1;  //<=========cm
		}
		//===================================
		//Re# falls outside of interval of available airfoil data
		else
		{	if(Re>HiRe)  //Re# above values of airfoi data file
			{
				printf("Re=%.0lf of section %d exceeding highest available Re#=%.0lf!\n",Re,section,HiRe);

				index=rows-1;

				//find index of CL of lower Re# just below cl of interest
				index1 = index;
				Re1 = profiledata[index1][3];	//lower Re#
				while(profiledata[index1][1]>cl && index1>=0 && profiledata[index1][3]==Re1)
					{index1--;}

				//cl falls outside of cl range (stalled or lower cl value) no cl-interpolation
				if(index1==index)
				{
					cd=profiledata[index1][2]; 			//section stalled
					cm=profiledata[index1][4]; 			//section stalled
					printf("4 STALL of section %d @ cl=%.3lf Re=%.0lf!\n",section,cl,Re);

					cl = .825*profiledata[index1][1];            //setting the cl to the 2D cl max of the airfoil data
				printf("4 STALL of section %d; cl from above set to cl max=%.01f\n", section,  cl);

				}
				else if(index1==-1 || profiledata[index1][3]!=Re1)
				{
					cd=profiledata[index1+1][2];		//section cl less than what is listed
					cm=profiledata[index1+1][4];		//section cl less than what is listed
					printf("4 Insufficient low cl airfoil data for section ");
					printf("%d cl=%.3lf and Re=%.0lf!\n",section,cl,Re);
				}

				else	//index1 is index of cl value below value of interest, thus interpolation between
						//index1 and index1+1
				{
					//interpolation slope (cd1-cd2)/(cl1-cl2)
					tempS= (cl-profiledata[index1][1])/\
						(profiledata[index1][1]-profiledata[index1+1][1]);

					m =	(profiledata[index1][2]-profiledata[index1+1][2]);
					cd=m*tempS + profiledata[index1][2];  //<=========cd

					//interpolation slope (cm1-cm2)/(cl1-cl2)
					m =	(profiledata[index1][4]-profiledata[index1+1][4]);
					cm=m*tempS + profiledata[index1][4];  //<=========cm
				}

			}
			else	//if(Re<LoRe)  //Re# less than values of airfoi data file
			{
				printf("Re=%.0lf of section %d less than available Re#=%.0lf!\n",Re,section,LoRe);

				index = -1;

				//find idnex of CL of higher Re# just above cl of interest
				index2=index+1;
				Re2 = profiledata[index2][3];	//upper Re#
				while(profiledata[index2][1]<cl  && index2<rows && profiledata[index2][3]==Re2)
					{index2++;}

				//cl falls outside cl range no interpolation
				if(index2==rows || profiledata[index2][3]!=Re2)
				{
					cd=profiledata[index2-1][2]; 			//section stalled
					cm=profiledata[index2-1][4]; 			//section stalled
					printf("5 STALL of section %d @ cl=%.3lf Re=%.0lf!\n",section,cl,Re);

					cl = .825*profiledata[index2-1][1];            //setting the cl to the 2D cl max of the airfoil data
					printf("5 STALL of section %d; cl from above set to cl max=%.01f\n", section,  cl);

				}
				else if(index2 == (index+1))
				{
					cd=profiledata[index2][2]; 			//section cl less than what is listed
					cm=profiledata[index2][4]; 			//section cl less than what is listed
					printf("5 Insufficient low cl airfoil data for section ");
					printf("%d cl=%.3lf and Re=%.0lf!\n",section,cl,Re);
						}
				else	//index2 is index of cl value above value of interest, thus interpolation between
						//index2-1 and index2
				{
					//interpolation slope (cd1-cd2)/(cl1-cl2)
					tempS=(cl-profiledata[index2-1][1])/\
						(profiledata[index2][1]-profiledata[index2-1][1]);

					m =	(profiledata[index2][2]-profiledata[index2-1][2]);
					cd=m*tempS + profiledata[index2-1][2];  //<=========cd

					//interpolation slope (cm1-cm2)/(cl1-cl2)
					m =	(profiledata[index2][4]-profiledata[index2-1][4]);
					cm=m*tempS + profiledata[index2-1][4];  //<=========cm
				}
			}
		}
//###
		return(cd);  //<=========cd
}
//================================================================================================================
//End SectionDrag
//================================================================================================================

