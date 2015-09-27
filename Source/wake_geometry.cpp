//creates new wake DVE right aft of trailing edge
void Squirt_out_Wake(const GENERAL,const PANEL *,DVE *,DVE *);
//relaxing the wake
void Relax_Wake(const GENERAL,const int,const DVE*,DVE**);
//computes point along the edge of a DVE
void Edge_Point(const double [3],const double,const double,const double,\
				const double,const double,const double,double [3]);
//computes new circulation&voritcity coefficients of stretched wake DVE
void New_vorticity_coefficients(const GENERAL,const PANEL *,DVE *,const DVE*);

//===================================================================//
		//START FUNCTION Squirt_out_Wake
//===================================================================//
void Squirt_out_Wake(const GENERAL info,const PANEL *panelPtr,\
					 DVE *surfacePtr,DVE *wakePtr)
//routine generates coordinates for the first wake DVE's that are located
//right aft of the trailing edge.   They have the same A,B,C
//coefficients as the surface DVE's that are located along the trailing edge.
//Thus, there influence is part of the D matrix.
//
//input
//	info		general information
//	panelPtr	panel information
//  surfacePtr	geometric information on surface DVE's
//
//output
//	wakePtr		array of wake DVE's that are located along trailing edge
//				at the current time step

{
	int k,j,wing;			//loop counters over panels, span, and wings
	int element = 0;		//index to surface elements along trailing edge
	int span = 0;			//span index for wake DVE's
	double ABS_u;			//absolute velocity
	double xTE[3];			//current surface DVE trailing edge center
	double TEvc[3],xto[3];	//vector along wing t.e., vector from t.e.to xo
	double ABS_TEvc;		//abs value of wing-t.e.
	double tempA[3],tempS;

	//loop over number of panels
	for(k=0;k<info.nopanel;k++)
	{
		//forwards surface-DVE index to trailing egde elements
		element += panelPtr[k].n*(info.m-1);

		//loop over number of spanwise elements of current panel
		for(j=0;j<panelPtr[k].n;j++)
		{
			//Computing current center of surface DVE trailing edge
			//xsi to trailing edge
			tempA[0]=surfacePtr[element].xsi;tempA[1]=0;tempA[2]=0;
			//transformation into global ref. frame
			Star_Glob(tempA,surfacePtr[element].nu,\
					  surfacePtr[element].epsilon,\
					  surfacePtr[element].psi,xTE);
					  				//function in ref_frame_transform.h

			//trailing edge pt.
			vsum(surfacePtr[element].xo,xTE,xTE);


			//local free stream velocity direction defines drag
			//direction and orientation of DVE-vortex sheet
				//motion of trailing edge center
			wakePtr[span].u[0] = surfacePtr[element].xTE[0]-xTE[0];
			wakePtr[span].u[1] = surfacePtr[element].xTE[1]-xTE[1];
			wakePtr[span].u[2] = surfacePtr[element].xTE[2]-xTE[2];
				scalar(wakePtr[span].u,1/info.deltime,wakePtr[span].u);

			ABS_u = norm2(wakePtr[span].u);
			//normalized velocity
			scalar(wakePtr[span].u,1/ABS_u,wakePtr[span].U);

			//delta x wing moved forward
			wakePtr[span].xsi=0.5*ABS_u*info.deltime;

			//computing vector along trailing edge of wing, xte
			//in local reference frame
			tempA[0] = tan(surfacePtr[element].phiTE)*surfacePtr[element].eta;
			tempA[1] = surfacePtr[element].eta;
			tempA[2] = 0;
			//transformation into global ref. frame
			Star_Glob(tempA,surfacePtr[element].nu,\
					  surfacePtr[element].epsilon,\
					  surfacePtr[element].psi,TEvc);
					  				//function in ref_frame_transform.h
			ABS_TEvc = norm2(TEvc); //|xte|

		//*****************************************************************
			//from trailing edge of wing (what is also the leading
			//edge of the first wake DVE) to first, post-wing DVE's ref point
		//*****************************************************************

			//vector from trailing edge of wing to ref. point of wakeDVE
			//this vector is parallel to local vel:
			tempS = ABS_u*0.5*info.deltime;
			xto[0] =  wakePtr[span].U[0]*tempS;
			xto[1] =  wakePtr[span].U[1]*tempS;
			xto[2] =  wakePtr[span].U[2]*tempS;

			//computing new reference point
			wakePtr[span].xo[0] = xTE[0] + xto[0];
			wakePtr[span].xo[1] = xTE[1] + xto[1];
			wakePtr[span].xo[2] = xTE[2] + xto[2];

//*****************************************************************************

			//computing normal of DVE
			cross(xto,TEvc,tempA);
			scalar(tempA,1/norm2(tempA),wakePtr[span].normal);

			//dihedral angle, nu,
			//tan(nu)=-(ny)/(nz)
			if(fabs(wakePtr[span].normal[2]) > DBL_EPS)
			{
				wakePtr[span].nu = \
						-atan(wakePtr[span].normal[1]/wakePtr[span].normal[2]);
				if (wakePtr[span].normal[2] < 0)
						wakePtr[span].nu += Pi;//|nu|>Pi/2
			}
			else  //tan(nu) -> infinity  -> |nu| = Pi/2
			{
				if(wakePtr[span].normal[1]>0)	wakePtr[span].nu = -0.5*Pi;
				else							wakePtr[span].nu =  0.5*Pi;
			}

			//incidence angle, epsilon
			wakePtr[span].epsilon = asin(wakePtr[span].normal[0]);

			//computing yaw angle psi new method G.B. 11-5-06
			//psi is the angle between free stream and xsi-axis

			//transformation of velocity into local ref. frame with psi=0
			Glob_Star(wakePtr[span].U,wakePtr[span].nu,\
					  wakePtr[span].epsilon,0,tempA);
					  				//function in ref_frame_transform.h

			//tan(psi)=(v)/(u)
			if(fabs(tempA[0])>DBL_EPS)
			{					wakePtr[span].psi = atan(tempA[1]/tempA[0]);
				if (tempA[0]<0)	wakePtr[span].psi += Pi;//|psi|>Pi/2
			}
			else  //tan(psi) -> infinity  -> |psi| = Pi/2
			{	if(tempA[1]>0)	wakePtr[span].psi =  0.5*Pi;
				else			wakePtr[span].psi = -0.5*Pi;
			}

			//computing leading edge sweep, phiLE
			//phiLE is the angle between the leading edge, vector xte,
			//and the xsi_star axis, or the vector u
			tempS = dot(TEvc,wakePtr[span].U)/ABS_TEvc;
			wakePtr[span].phiLE = asin(tempS);


			//computing trailing edge phi
			//new method G.B. 11-5-06

			//projection of TE vector of previous timestep
			//on plane of new wakeDVE
			Glob_Star(surfacePtr[element].TEvc,wakePtr[span].nu,\
					  wakePtr[span].epsilon,wakePtr[span].psi,tempA);
					  				//function in ref_frame_transform.h
			//phiTE is the angle between the trailing edge, vector TEvc
			//and the xsi_star axis, or the vector u
			wakePtr[span].phiTE = atan(tempA[0]/tempA[1]);

			//sweep at mid chord
			wakePtr[span].phi0=0.5*(wakePtr[span].phiLE+wakePtr[span].phiTE);

			//DVE half-span, projection of xte onto eta_star
			wakePtr[span].eta = ABS_TEvc*cos(wakePtr[span].phiLE);

			//circulation is shed into wake, is corrected for possibly
			//changed span in wake-relaxation routine
			wakePtr[span].A = surfacePtr[element].A;
			wakePtr[span].B = surfacePtr[element].B;
			wakePtr[span].C = surfacePtr[element].C;

			//total circulation of DVE, stays constant despite stretching
			wakePtr[span].K = surfacePtr[element].A \
							+ surfacePtr[element].eta\
							* surfacePtr[element].eta/3\
							* surfacePtr[element].C;

			//computes point halfway along left edge of DVE
			Edge_Point(wakePtr[span].xo,wakePtr[span].nu,\
					   wakePtr[span].epsilon,wakePtr[span].psi,\
					   wakePtr[span].phi0,-wakePtr[span].eta,\
					   0,wakePtr[span].xleft);
				 					//subroutine in wake_geometry.cpp

			//updatin surface DVE's trailing edge center point
			surfacePtr[element].xTE[0]=xTE[0];
			surfacePtr[element].xTE[1]=xTE[1];
			surfacePtr[element].xTE[2]=xTE[2];
			//updatin surface DVE's trailing edge vector
			surfacePtr[element].TEvc[0]=TEvc[0];
			surfacePtr[element].TEvc[1]=TEvc[1];
			surfacePtr[element].TEvc[2]=TEvc[2];

/*
printf("\nxo %2.3lf %2.3lf %2.3lf\n",\
wakePtr[span].xo[0],wakePtr[span].xo[1],wakePtr[span].xo[2]);
printf("eta %2.3lf xsi %2.3lf",\
wakePtr[span].eta,wakePtr[span].xsi);
printf("K %2.3lf A %2.3lf B %2.3lf C %2.3lf ",\
wakePtr[span].K,wakePtr[span].A,wakePtr[span].B,wakePtr[span].C);
printf("nu %2.3lf eps %2.3lf psi %2.3lf phiTE %2.3lf\n",\
wakePtr[span].nu*RtD,wakePtr[span].epsilon*RtD,wakePtr[span].psi*RtD,wakePtr[span].phiTE*RtD);
printf("xl %2.3lf %2.3lf %2.3lf ",\
wakePtr[span].xleft[0],wakePtr[span].xleft[1],wakePtr[span].xleft[2]);
//*/

//*****************************************************************************

			element ++;		//advance surface-DVE index
			span ++;		//advance wake-DVE index
		}	//end loop over j, DVE's along trailing edge
	}	// end loop k over nopanel

	//##this part has been added 2/9/05 G.B.
	//decaying factor for added singularity at wing tip
	for(wing=0;wing<info.nowing;wing++)
	{
		//is the wing symmetrical or not?
		if(info.sym == 1)//decay factor is 1% of tip-element half-span
			tempS = 0.01*wakePtr[info.wing2[wing]].eta;
		else//wing has two tips, possibly different in geometry
		{	//in that case, decay factor is 1% of the shorter half-span
			if(wakePtr[info.wing1[wing]].eta < wakePtr[info.wing2[wing]].eta)
						tempS = 0.01*wakePtr[info.wing1[wing]].eta;
			else 		tempS = 0.01*wakePtr[info.wing2[wing]].eta;
		}
		//loop over wale DVEs of current timestep
		for(span=info.wing1[wing];span<=info.wing2[wing];span++)
			wakePtr[span].singfct = tempS;  //assigning decay factor
	}//next wing


//DELETED AND REPLACEd WITH ROUTINE IN MAIN 5/27/2005 G.B.
	//adjusting circulation for changed span
//	New_vorticity_coefficients(info,panelPtr,wakePtr,wakePtr);
		 							//subroutine in wake_geometry.cpp
}
//===================================================================//
		//END FUNCTION Squirt_out_Wake
//===================================================================//

//===================================================================//
		//FUNCTION Relax_Wake
		//relaxes wake
//===================================================================//
void Relax_Wake(const GENERAL info,const int rightnow,\
				const DVE* surfacePtr,DVE** wakePtr)
{
	//relaxing the wake:
	//	1. computes local induced velocity at side edges of DVEs
	//	2. displaces of side edges of DVEs
	//	3. computes new ref. pt of wake DVE
	//	4. computes new eta, nu, epsilon, and psi, as well as new xsi
	//  5. DELETED AND REPLACED WITH ROUTINE IN MAIN 5/27/2005 G.B.
	//	6. computes new singularity factor for
	//
	//	Function determines the new ref. point location by moving left and
	//  right edges.
	//The new xo, u, eta, nu are determined with the help of the displaced
	//points at either edge.
	//epsilon and psi depend on the connection the DVE further upstream
	//
	//input:
	//	info		- general information
	//	rightnow	- current time step
	//	surfacePtr	- information on surface DVE's
	//	wakePtr		- information on wake DVE's
	//
	//ouput:
	// 	For each wake DVE in wakePtr updated
	//	xo[]		- reference point in global ref.frame
	//  nu,epsilon  - dihedral, incident of DVE
	//	psi			- yaw angle of DVE
	//  xsi			- half chord length in psi-rotated ref. frame

//displaces point along local streamline
void Displace_Point(const GENERAL,const int,const double [3],\
					const double [3],const double [3],double x[3]);
//computes new ref. pt and vector between pt. at left edge and ref. pt.
void New_xo(const GENERAL,DVE &,double [3],double [3]);
//computes new roll, pitch and yaw angles, half chord, and half span
void New_eta_nu_eps_psi_xsi(DVE &,const DVE);
//computes new ref. pt, roll, pitch, yaw, half chord, xleft of first wake DVE
void New_wakeDVE0(const GENERAL,DVE *,const DVE *);

int time,span,wing,j,k;		//loop counters
int element;				//index of surface elements along trailing edge
double tempS,tempA[3];
double ***xright,***uright;	//location and local velocity on right edge of a wing

//##########################################
//FILE *fp;														//#
//fp = fopen("output\\test.txt", "a");									//#
//fprintf(fp,"timestep %d\n",rightnow);//#
//###########################################################################


//allocating memory
ALLOC3D(&xright,info.nowing,rightnow+1,3);
ALLOC3D(&uright,info.nowing,rightnow+1,3);

//*****************************************************************************
//	1. computes local induced velocity at side edges of DVEs

//&&	for(time=0;time<=rightnow;time++)
	for(time=1;time<=rightnow;time++)	//&&
	{
		//computing the locally induced velocities at left edge points
		for(span=0;span<info.nospanelement;span++)
		{
			DVE_Induced_Velocity(info,wakePtr[time][span].xleft,\
						surfacePtr,wakePtr,rightnow,wakePtr[time][span].u);
				 						//subroutine in induced_velocity.cpp
		}//next span element

		//computing induced velocity of most right point of each wing.
		//right wingtip points are stored in xright and the velocity in uright
		for(wing=0;wing<info.nowing;wing++)
		{
			//span index of DVE that is the most right on of current wing
			span = info.wing2[wing];

			//find point halfway along right edge
			//computes point halfway along left edge of element
			//temporarily stored in .normal of most right DVE
			Edge_Point(wakePtr[time][span].xo,wakePtr[time][span].nu,\
					   wakePtr[time][span].epsilon,wakePtr[time][span].psi,\
					   wakePtr[time][span].phi0,wakePtr[time][span].eta,\
					   0,xright[wing][time]);
					 					//subroutine in wake_geometry.cpp

			//computing velocity induced at right tip of wing
			DVE_Induced_Velocity(info,xright[wing][time],surfacePtr,wakePtr,\
								rightnow,uright[wing][time]);
				 						//subroutine in induced_velocity.cpp

//computing induced velocity of wake DVE [time][span] on point P.
//Single_DVE_Induced_Velocity(info,wakePtr[time][span],xright[wing][time],tempA,1);
				 						//subroutine in induced_velocity.cpp
//uright[wing][time][0] -= tempA[0];
//uright[wing][time][1] -= tempA[1];
//uright[wing][time][2] -= tempA[2];

//##########################################
//fprintf(fp,"xright: %lf %lf %lf ",xright[wing][time][0],xright[wing][time][1],xright[wing][time][2]);//#
//fprintf(fp," %lf %lf %lf ",tempA[0],tempA[1],tempA[2]);//#
//fprintf(fp," %lf %lf %lf\n",uright[wing][time][0],uright[wing][time][1],uright[wing][time][2]);//#
//fprintf(fp,"c3: %lf\t%lf\t%lf\tC %lf\n",c3[0],c3[1],c3[2],DVelement.C);//#
//###########################################################################


		}//next wing
	}//next spanwise strip in wake (next time)

//##########################################
//fclose (fp);														//#
//###########################################################################

//*****************************************************************************
//	2. displaces of side edges of DVEs

//&&	for(time=1;time<rightnow;time++)
	for(time=2;time<rightnow;time++)	//&&
	{
		//displacing left edge points
		for(span=0;span<info.nospanelement;span++)
		{
			Displace_Point(info,rightnow,wakePtr[time-1][span].u,\
							wakePtr[time][span].u,wakePtr[time+1][span].u,\
								   				wakePtr[time][span].xleft);
		}//next span element		//subroutine in wake_geometry.cpp

		//displacing edge points along right wing tip
		for(wing=0;wing<info.nowing;wing++)
		{
			Displace_Point(info,rightnow,uright[wing][time-1],\
							uright[wing][time],uright[wing][time+1],\
							xright[wing][time]);
										//subroutine in wake_geometry.cpp
		}//next wing
	}//next spanwise strip in wake (next time)

	//the zero vector
	tempA[0] = 0; tempA[1] = 0; tempA[2] = 0;

	//displacing left edge points of first and last row of wake DVE's
//&	//(time index = rightnow and 0)
	//(time index = rightnow and 1)		//&
	for(span=0;span<info.nospanelement;span++)
	{
		//post trailing edge
		Displace_Point(info,rightnow,\
						wakePtr[rightnow][span].u,\
						tempA,\
						wakePtr[rightnow-1][span].u,\
						wakePtr[rightnow][span].xleft);
							//subroutine in wake_geometry.cpp

		//second oldest row of DVE's		//&&
		Displace_Point(info,rightnow,wakePtr[1][span].u,wakePtr[1][span].u,\
						wakePtr[1][span].u,wakePtr[1][span].xleft);		//&&
								//subroutine in wake_geometry.cpp
		//very first row of DVE's
//&&		Displace_Point(info,rightnow,wakePtr[0][span].u,wakePtr[0][span].u,\
//&&						wakePtr[0][span].u,wakePtr[0][span].xleft);
	}//next span element		//subroutine in wake_geometry.cpp

	//displacing edge points along right wing tip
	for(wing=0;wing<info.nowing;wing++)
	{
		//post trailing edge
		Displace_Point(info,rightnow,\
						uright[wing][rightnow],\
						tempA,\
						uright[wing][rightnow-1],\
						xright[wing][rightnow]);
									//subroutine in wake_geometry.cpp
		//second oldest row of wake DVEs	//&&
		Displace_Point(info,rightnow,uright[wing][1],uright[wing][1],\
						uright[wing][1],xright[wing][1]);	//&&
		//first row of wake DVEs
//&&		Displace_Point(info,rightnow,uright[wing][0],uright[wing][0],\
//&&						uright[wing][0],xright[wing][0]);
									//subroutine in wake_geometry.cpp
	}//next wing

//*****************************************************************************
//	3. computes new reference point location and vector between the two
//	   displaced points at edge (stored in .normal)
//	   In the case of the first, post-wing row of DVEs, xo is actually the midspan
//     location along the DVEs trailing edges

//&&	for(time=0;time<=rightnow;time++)
	for(time=1;time<=rightnow;time++)		//&&
	//loop over number of wings
	for(wing=0;wing<info.nowing;wing++)
	{
		//loop across span of current wing, up to one before wingtip
		for(span=info.wing1[wing];span<info.wing2[wing];span++)
		{
			New_xo(info,wakePtr[time][span],wakePtr[time][span].xleft,\
											wakePtr[time][span+1].xleft);
	 									//subroutine in wake_geometry.cpp
		}
		//and now the right wing tip
		span=info.wing2[wing];
		New_xo(info,wakePtr[time][span],wakePtr[time][span].xleft,\
				xright[wing][time]);	//subroutine in wake_geometry.cpp
	}

//*****************************************************************************
//  4. computes new eta, nu, epsilon, and psi, as well as new xsi

	//computes new angles for most recently squirted out DVEs
	//Attention! computed along their trailing edge, not mid chord
	element = 0;	//initializing index for surface DVE at trailing edge
	span = 0;		//initializing index for wake DVE along DVE
	//loop over number of panels
	for(k=0;k<info.nopanel;k++)
	{
		//forwards surface-DVE index to trailing egde elements
		element += panelPtr[k].n*(info.m-1);

		//loop over number of spanwise elements of current panel
		for(j=0;j<panelPtr[k].n;j++)
		{
			New_eta_nu_eps_psi_xsi
							(wakePtr[rightnow][span],surfacePtr[element]);
		 								//subroutine in wake_geometry.cpp
			element ++;		//advance surface-DVE index
			span ++;		//advance wake-DVE index
		}	//end loop over j, DVE's along trailing edge
	}	// end loop k over nopanel

	//computes new angles for remaining wake field
	for(span=0;span<info.nospanelement;span++)
//&&	for(time=rightnow;time>0;time--)
	for(time=rightnow;time>1;time--)  //&&
		New_eta_nu_eps_psi_xsi(wakePtr[time-1][span],wakePtr[time][span]);
		 								//subroutine in wake_geometry.cpp

//*****************************************************************************
//	4a. same as 3. and 4, just for the first row of wake DVEs (time index = 0)
//		these are aligned with the free stream velocity vector

	//&&&&
	New_wakeDVE0(info,wakePtr[0],wakePtr[1]);

//*****************************************************************************
//DELETED AND REPLACED WITH ROUTINE IN MAIN 5/27/2005 G.B.
//	5. computes new vorticity coefficients A, B, and C

//	for(time=0;time<=rightnow;time++)
//	{
//		if (info.steady == 1)
//			 j=rightnow;//steady airloads, wake DVEs of each span location have
//						//equal amount integrated circulation
//		else j=time;	//unsteady airloads, varying integrated circulation, k
//
//		New_vorticity_coefficients(info,panelPtr,wakePtr[time],wakePtr[j]);
//		 							//subroutine in wake_geometry.cpp
//	}

//*****************************************************************************
//	6. computes new singularity factor

	//loop over timesteps
	for(time=0;time<=rightnow;time++)
	//loop over number of wings
	for(wing=0;wing<info.nowing;wing++)
	{
		//is the wing symmetrical or not?
		if(info.sym == 1)	//decay factor is 1% of tip-element half-span
			tempS = 0.01*wakePtr[time][info.wing2[wing]].eta;
		else//wing has two tips, possibly different in geometry
		{	//in that case, decay factor is 1% of the shorter half-span
			if(  wakePtr[time][info.wing1[wing]].eta
			   < wakePtr[time][info.wing2[wing]].eta)
						tempS = 0.01*wakePtr[time][info.wing1[wing]].eta;
			else 		tempS = 0.01*wakePtr[time][info.wing2[wing]].eta;
		}
		//loop over wale DVEs of current timestep
		for(span=info.wing1[wing];span<=info.wing2[wing];span++)
			wakePtr[time][span].singfct = tempS;  //assigning decay factor
	}//next wing

//*///****************************************************************************
//freeing allocated memory
FREE3D(&xright,info.nowing,rightnow+1,3);
FREE3D(&uright,info.nowing,rightnow+1,3);
}
//===================================================================//
		//END FUNCTION Relax_Wake
//===================================================================//
//===================================================================//
		//FUNCTION Edge_Point
//===================================================================//
void Edge_Point(const double xo[3],const double nu,const double epsilon,\
				const double psi,const double phi,const double eta,\
				const double xsi,double x[3])
{
	//
	//	function computes the point that is located the y-distance,
	// 	eta, from xo.
	//
	//input:
	//  xo			- reference point
	//	nu			- angle, at which y-distance is measured
	//	epsilon		- incidence of DVE
	//	psi			- yaw angle of DVE
	//	phi			- sweep of mid-chord line of DVE
	//	eta			- y-distance to point of interest
	//	xsi			- x-distance to point of interest
	//
	//ouput:
	// 	For each wake DVE in wakePtr updated
	//	x[]			- point of interest in global ref.frame

	double edge[3];		//edge location in local ref. frame

	edge[0] = xsi+eta*tan(phi); edge[1]=eta; edge[2]=0;

	//transformation of xsi-edge into global reference frame
	Star_Glob(edge,nu,epsilon,psi,x);//function in ref_frame_transform.h

	//with respect to reference point
	x[0] += xo[0];
	x[1] += xo[1];
	x[2] += xo[2];
}
//===================================================================//
		//END FUNCTION Edge_Point
//===================================================================//

//===================================================================//
		//FUNCTION Displace_Point
//===================================================================//
void Displace_Point(const GENERAL info,const int rightnow,\
					const double uk[3],const double ul[3],const double um[3],\
					double x[3])
{
	//
	//	function displaces a point along the local streamline
	//
	//input:
	//	info		- general information
	//	rightnow	- current time step
	//	uk			- local velocity at midchord of DVE downstream
	//	ul			- local velocity at midchord of DVE of interest
	//	um			- local velocity at midchord of DVE upstream
	//  x			- point that is displaced
	//
	//ouput:
	// 	For each wake DVE in wakePtr updated
	//	x[]			- point of interest in global ref.frame

	double u[3];	  		// local induced velocity

	u[0] = 0.25*(uk[0] + 2*ul[0] + um[0]);
	u[1] = 0.25*(uk[1] + 2*ul[1] + um[1]);
	u[2] = 0.25*(uk[2] + 2*ul[2] + um[2]);

//#printf("\nuk %2.3lf %2.3lf %2.3lf",uk[0],uk[1],uk[2]);
//#printf("  ul %2.3lf %2.3lf %2.3lf",ul[0],ul[1],ul[2]);
//#printf("  um %2.3lf %2.3lf %2.3lf",um[0],um[1],um[2]);

	//new xleft = xleft_old + local Uind * delta time
	x[0] += u[0]*info.deltime*1.0;
	x[1] += u[1]*info.deltime*1.0;
	x[2] += u[2]*info.deltime*1.0;
}
//===================================================================//
		//END FUNCTION Displace_Point
//===================================================================//

//===================================================================//
		//START FUNCTION New_xo
//===================================================================//
void New_xo(const GENERAL info,DVE &wakeDVE,double x_left[3],double x_right[3])
{
	//
	//This function computes a new reference point and the vector
	//between left and right edge
	//
	//input:
	//	info		- general information
	//	x_left		- point to the left
	//	x_right		- point to the right
	//	wakeDVE		- information of wake DVE of interest
	//
	//ouput:
	// 	For the wake DVE of interest the following values are being updated
	//	x[]			- reference point in global ref.frame
	//	u[]			- velocity in x[]

	double one_over_time = 1/info.deltime;
	double tempA[3];

	//computing vector between left edge and right edge
	wakeDVE.normal[0] = 0.5*(x_right[0] - x_left[0]);
	wakeDVE.normal[1] = 0.5*(x_right[1] - x_left[1]);
	wakeDVE.normal[2] = 0.5*(x_right[2] - x_left[2]);

	//saving old ref. pt. location temporarily saved
	tempA[0] = wakeDVE.xo[0];
	tempA[1] = wakeDVE.xo[1];
	tempA[2] = wakeDVE.xo[2];

	//determining new reference-point location
	wakeDVE.xo[0]= x_left[0] + wakeDVE.normal[0];
	wakeDVE.xo[1]= x_left[1] + wakeDVE.normal[1];
	wakeDVE.xo[2]= x_left[2] + wakeDVE.normal[2];

	//induced velocity in reference pt is delta x/time
	wakeDVE.u[0] = (wakeDVE.xo[0] - tempA[0])*one_over_time;
	wakeDVE.u[1] = (wakeDVE.xo[1] - tempA[1])*one_over_time;
	wakeDVE.u[2] = (wakeDVE.xo[2] - tempA[2])*one_over_time;
//printf("  u %2.3lf %2.3lf %2.3lf",wakeDVE.u[0],wakeDVE.u[1],wakeDVE.u[2]);

//printf("\nxl %2.3lf %2.3lf %2.3lf xr %2.3lf %2.3lf %2.3lf",x_left[0],x_left[1],x_left[2],x_right[0],x_right[1],x_right[2]);
//printf("\nxo-xl %2.3lf %2.3lf %2.3lf",wakeDVE.normal[0],wakeDVE.normal[1],wakeDVE.normal[2]);
//printf("  xo %2.3lf %2.3lf %2.3lf ",wakeDVE.xo[0],wakeDVE.xo[1],wakeDVE.xo[2]);
}
//===================================================================//
		//END FUNCTION New_xo
//===================================================================//

//===================================================================//
		//START FUNCTION New_eta_nu_eps_psi_xsi
//===================================================================//
void New_eta_nu_eps_psi_xsi(DVE &wakeDVE,const DVE US_DVE)
{
//This function computes the roll, pitch and sweep angles, as well the
//halfspan of a wake DVE after it has been "relaxed".
//nu and epsilon are computed from the normal of the DVE plane.  The normal
//is the crossproduct of the vector delX1 and delX2.  delX1 is the vector
//from the half point of the trailing edge of the next upstream DVE to the new
//reference point. delX2 is the vector between the left side-edge point and the
//reference point.  This vector also is used to determine the change in sweep
//and the new halfspan of the DVE.
//Yaw, psi, is kept constant as is the halfchord xsi.
//
//
//input:
//	wakeDVE	- DVE of interest
//	US_DVE	- DVE upstream
//
//output:
//updated values for wakeDVE
//	nu		- roll angle
//	epsilon - pitch (incidence) angle
//	psi		- yaw angle
//	xsi		- half chord length (in local, psi-rotated system)
//	eta		- half span

	double xteUS[3];	//coord. of trailing edge point of upstream DVE
	double delX1[3];	//vector between ref. pt and xoUS
	double delX2[3];	//vector between ref. pt and left edge point
	double delXSI1[3];	//delX1 in local DVE coordinates
	double delXSI2[3];	//delX2 in local DVE coordinates
	double delPhi;		//change in sweep
	double tempA[3];


//********************************************************************
	//computing trailing edge point of upstream DVE

	//from reference point to trailing edge of upstream DVE:
	//locally, xsi to trailing edge
	tempA[0]=US_DVE.xsi;tempA[1]=0;tempA[2]=0;
	//transformation into global ref. frame
	Star_Glob(tempA,US_DVE.nu,US_DVE.epsilon,US_DVE.psi,xteUS);
					  				//function in ref_frame_transform.h

	//trailing edge pt. temporarily in wakePtr[span].xo
	vsum(xteUS,US_DVE.xo,xteUS);

//********************************************************************
	//computing vector from trailing edge pt.of upstream DVE to
	//ref. pt of wakeDVe
	//P_downstream - _upstream
	delX1[0] = wakeDVE.xo[0] - xteUS[0];
	delX1[1] = wakeDVE.xo[1] - xteUS[1];
	delX1[2] = wakeDVE.xo[2] - xteUS[2];

	//vector from left edge mid-chord location and ctrl. point
	delX2[0] = wakeDVE.normal[0];
	delX2[1] = wakeDVE.normal[1];
	delX2[2] = wakeDVE.normal[2];
//#printf("X2 %2.3lf %2.3lf %2.3lf ",wakeDVE.normal[0],wakeDVE.normal[1],wakeDVE.normal[2]);

	//computing the normal of DVE,
	cross(delX1,delX2,tempA);
	scalar(tempA,1/norm2(tempA),wakeDVE.normal);
//#printf("n %2.3lf %2.3lf %2.3lf %2.2lf ",wakeDVE.normal[0],wakeDVE.normal[1],wakeDVE.normal[2],norm2(wakeDVE.normal));

	//here begins the stuff you been waiting for:

//********************************************************************
	//dihedral angle, nu,
	//tan(nu)=-(ny)/(nz)
	if(fabs(wakeDVE.normal[2]) > DBL_EPS)
	{
		wakeDVE.nu = -atan(wakeDVE.normal[1]/wakeDVE.normal[2]);
		if (wakeDVE.normal[2] < 0)  wakeDVE.nu += Pi;//|nu|>Pi/2
	}
	else  //tan(nu) -> infinity  -> |nu| = Pi/2
	{
		if(wakeDVE.normal[1]>0) 	wakeDVE.nu = -0.5*Pi;
		else 						wakeDVE.nu =  0.5*Pi;
	}
//********************************************************************
	//computing epsilon
	wakeDVE.epsilon = asin(wakeDVE.normal[0]);

//********************************************************************
	//computing yaw angle psi
		//the orientation of the rotational axis of the vortex seet,
		//is parallel to delX1.

	//delX1 in xsi reference frame (psi = 0)
	Glob_Star(delX1,wakeDVE.nu,wakeDVE.epsilon,0,delXSI1);
					  				//function in ref_frame_transform.h
//#printf("\n XSI1 %2.3lf %2.3lf %2.3lf ",delXSI1[0],delXSI1[1],delXSI1[2]);

	//tan(psi)=(eta2)/(xsi1)
	if(delXSI1[0]*delXSI1[0] > DBL_EPS)
	{
		wakeDVE.psi= atan(delXSI1[1]/delXSI1[0]);
			if(delXSI1[0] < 0)  wakeDVE.psi += Pi;//|nu|>Pi/2
	}
	else	//tan(psi) -> infinity  -> |psi| = Pi/2
		wakeDVE.psi = 0.5*Pi*fabs(delXSI1[1])/delXSI1[1];

	//length of projection
	wakeDVE.xsi = sqrt(delXSI1[0]*delXSI1[0] + delXSI1[1]*delXSI1[1]);

//********************************************************************
	//computing sweep and halfspan

	//1. transform delX2 into local frame
	Glob_Star(delX2,wakeDVE.nu,wakeDVE.epsilon,wakeDVE.psi,delXSI2);
					  				//function in ref_frame_transform.h
//#printf(" XSI2 %2.3lf %2.3lf %2.3lf ",delXSI2[0],delXSI2[1],delXSI2[2]);

	//if everything works => delXSI[2] = 0
	if(delXSI2[1] < 0)
	{
		printf(" ohwei! around line 788, wake_geometry.cpp\n %2.2lf %2.16lf\n",delXSI2[1],delXSI2[2]);
				exit(0);
	}

	//computing change in sweep
	delPhi = atan(delXSI2[0]/delXSI2[1]) - wakeDVE.phi0;

	//updating sweeps
	wakeDVE.phiLE += delPhi;
	wakeDVE.phi0  += delPhi;
	wakeDVE.phiTE += delPhi;

	//computing the new half span
	wakeDVE.eta = delXSI2[1];//*cos(wakeDVE.phi0);//fabs(delXSI2[1]);

//********************************************************************
//#printf("\nnu %2.2lf e %2.2lf ph %2.2lf psi %2.2lf eta %2.4lf  ",wakeDVE.nu*RtD,wakeDVE.epsilon*RtD,wakeDVE.phi0*RtD,wakeDVE.psi*RtD,wakeDVE.eta);
}
//===================================================================//
		//END FUNCTION New_eta_nu_eps_psi_xsi
//===================================================================//

//===================================================================//
		//START FUNCTION New_wakeDVE0
//===================================================================//
void New_wakeDVE0(const GENERAL info,DVE *wakeDVE0,const DVE *wakeDVE1)
{
	//
	//This function computes the reference point, angles, span and left
	//edge points of the first row of wake DVEs (time index = 0).
	//These are aligned with the free stream and are attached at their
	//leading edge to the trailing edges of the the upstream DVEs (index 1)
	//added June, 11 2005, G.B.
	//
	//input:
	//	info		- general information
	//	wakeDVE0	- first row of wake DVEs, time index = 0
	//	wakeDVE1	- row of wake DVEs directly upstream of wakeDVE0
	//
	//ouptut:
	//  wakeDVE0	- updated xo, nu, epsilon, psi, eta, xleft

	int span;		//span index
	double xte[3];	//mid location as well as direction of TE of upstream DVE
	double ABS_xte;	//magnitude of xte
	double tempA[3],tempS;

	//the first row of wake DVEs
	for(span=0;span<info.nospanelement;span++)
	{

	//********************************************************************
		//computing the reference pont location

		//computing trailing edge point of upstream DVE
		//from reference point to trailing edge of upstream DVE:
		//locally, xsi to trailing edge
		tempA[0]=wakeDVE1[span].xsi;	tempA[1]=0;	tempA[2]=0;
		//transformation into global ref. frame
		Star_Glob(tempA,wakeDVE1[span].nu,wakeDVE1[span].epsilon,\
						wakeDVE1[span].psi,xte);
					  				//function in ref_frame_transform.h
		xte[0] +=  wakeDVE1[span].xo[0];
		xte[1] +=  wakeDVE1[span].xo[1];
		xte[2] +=  wakeDVE1[span].xo[2];

		//the new ref. point is located along the extension of the free stream
		wakeDVE0[span].xo[0] = xte[0] + wakeDVE0[span].xsi*wakeDVE0[span].U[0];
		wakeDVE0[span].xo[1] = xte[1] + wakeDVE0[span].xsi*wakeDVE0[span].U[1];
		wakeDVE0[span].xo[2] = xte[2] + wakeDVE0[span].xsi*wakeDVE0[span].U[2];
	//********************************************************************
		//determining sweep angles and half span

		//computing vector along trailing edge of wing, xte
		//in local reference frame
		tempA[0] = tan(wakeDVE1[span].phiTE)*wakeDVE1[span].eta;
		tempA[1] = wakeDVE1[span].eta;
		tempA[2] = 0;
		//transformation into global ref. frame
		Star_Glob(tempA,wakeDVE1[span].nu,wakeDVE1[span].epsilon,\
				  wakeDVE1[span].psi,xte);  //function in ref_frame_transform.h

		ABS_xte = norm2(xte); //|xte|

		//computing leading edge sweep, phiLE
		//phiLE is the angle between the leading edge, vector xte,
		//and the xsi_star axis, or the vector u
		tempS = dot(xte,wakeDVE0[span].U)/ABS_xte;
		wakeDVE0[span].phiLE = asin(tempS);
		//the remaining angles are of equal value since
		//they don't really matter for semi-infinite sheet
		wakeDVE0[span].phi0  = wakeDVE0[span].phiLE;
		wakeDVE0[span].phiTE = wakeDVE0[span].phiLE;

		//DVE half-span, projection of xte onto eta_star
		wakeDVE0[span].eta = ABS_xte*cos(wakeDVE0[span].phiLE);

	//********************************************************************
		//determining nu, epsilon, and psi

		//computing the normal of DVE,
		cross(wakeDVE0[span].U,xte,tempA);
		scalar(tempA,1/norm2(tempA),wakeDVE0[span].normal);

		//dihedral angle, nu,
		//tan(nu)=-(ny)/(nz)
		if(fabs(wakeDVE0[span].normal[2]) > DBL_EPS)
		{
			wakeDVE0[span].nu = \
					-atan(wakeDVE0[span].normal[1]/wakeDVE0[span].normal[2]);
			if (wakeDVE0[span].normal[2] < 0)
									wakeDVE0[span].nu += Pi;//|nu|>Pi/2
		}
		else  //tan(nu) -> infinity  -> |nu| = Pi/2
		{
			if(wakeDVE0[span].normal[1]>0)	wakeDVE0[span].nu = -0.5*Pi;
			else							wakeDVE0[span].nu =  0.5*Pi;
		}

		//epsilon
		wakeDVE0[span].epsilon = asin(wakeDVE0[span].normal[0]);

		//yaw angle psi
		//the angle between free stream and xsi-axis
		//or between eta and u minus 90deg
		tempS = (wakeDVE0[span].U[1] * cos(wakeDVE0[span].nu)\
				+wakeDVE0[span].U[2] * sin(wakeDVE0[span].nu));
		wakeDVE0[span].psi = asin(tempS);

	//********************************************************************
//		//updating left edge point
//		Edge_Point(wakeDVE0[span].xo,wakeDVE0[span].nu,\
//				   wakeDVE0[span].epsilon,wakeDVE0[span].psi,\
//				   wakeDVE0[span].phiTE,-wakeDVE0[span].eta,0,\
//				   wakeDVE0[span].xleft);  //subroutine in wake_geometry.cpp

	//********************************************************************
		//assinging a local velocity
		wakeDVE0[span].u[0] = wakeDVE1[span].u[0];
		wakeDVE0[span].u[1] = wakeDVE1[span].u[1];
		wakeDVE0[span].u[2] = wakeDVE1[span].u[2];
	}//next span location
}
//===================================================================//
		//END FUNCTION New_wakeDVE0
//===================================================================//

//===================================================================//
		//START FUNCTION New_vorticity_coefficients
//===================================================================//
void New_vorticity_coefficients(const GENERAL info,const PANEL *panePtr,\
								DVE *wakeDVE,const DVE *newestDVE)
{
//adjusts vorticity distribution in wake strips in order to compensate
//for stretching of wake DVE's.  In an irrotational flwow, the total
//circulation remains constant in the wake.
//Thus 1/2eta_i*Int{(A+eta*B+eta^2*C) dy} = const !!!
//For each timestep strip in the wake of each wing this leads to an system
//of 2n equations for 2n unknowns.  The equations cover the conservation
//of circulation and continuity of vorticity between the DVE's
//
//input:
//	info 		- general information
//	wakeDVE		- wake DVE of one timestep
//	newestDVE	- if steady airloads, the integrated ciruculation, the k-value,
//				  has the uniform value of the first post-trailing edge element
//				  for all DVEs  of one span location.  If unsteady airloads,
//				  k-value is individual for each DVE of a certain span location
//				  DOESN'T WORK, YET!? G.B. 10/30/04
//
//output:
// updated circulation and vorticity coefficients of stretched wakeDVE


	double **D,*R,*BC;				//2nx2n matrix, 2n right hand side, new coeff.
	double twothirds = 2./3.;
	int size=2*info.nospanelement;	//dimension of D matrix and RHS-vector
	int h,i;						//indices of DVEs, h is left i
	int panel,n;					//loop over panel, spanwise elements
	int col,row;					//column and row of D matrix

  	//allocate memory
	ALLOC2D(&D,size,size);
	ALLOC1D(&R,size);
	ALLOC1D(&BC,size);

	//initializing D, R and BC
	for(row=0;row<size;row++)
	{
		for(col=0;col<size;col++)
			D[row][col]=0;
		R[row]  = 0;
		BC[row] = 0;
	}

	///////////////////////
	//Assembling D-matrix//
	///////////////////////

	//initializing counter
	i = 0;

	for(panel=0;panel<info.nopanel;panel++)		//loop over panels
	{
//////////////////////////////////////////////////////////////////
		//finding boundart condition of the first (left) DVE
		//of particular timestep in the wake of a particular panel

		col = i*2;
		row = col;

		switch (panelPtr[panel].BC1)	//edge 1 bound. cond.
		{
			case 100: 	//gamma = 0
				D[row][col]   = wakeDVE[i].eta;
				D[row][col+1] = -twothirds*wakeDVE[i].eta*wakeDVE[i].eta;
				R[row] 		  = newestDVE[i].K;
			break;

			case 110: 	//gamma' = 0
				D[row][col]   = 1.;
				D[row][col+1] = -2*wakeDVE[i].eta;  //##!!
				R[row] 		  = 0.;
			break;

			case 10:	//gamma' = 0
				D[row][col]   = 1.;
				D[row][col+1] = -2*wakeDVE[i].eta;  //##!!
				R[row] 		  = 0.;
			break;

			case 220:	//gamma[i] = gamma[j] AND gamma[i]' = gamma[j]'
				h = i - 1;  //index of element to left
				if (h<0)
				{
					printf("\n WARNING!! Left boundary condition of ");
					printf("panel %d messed up!!!",panel+1);
				}

				D[row][col]   = wakeDVE[i].eta;
				D[row][col+1] = -twothirds*wakeDVE[i].eta*wakeDVE[i].eta;
				D[row][col-2] = wakeDVE[h].eta;
				D[row][col-1] = twothirds*wakeDVE[h].eta*wakeDVE[h].eta;
				R[row] 		  = newestDVE[i].K - newestDVE[h].K;

				row --; //previous row
				D[row][col]   = -1.;
				D[row][col+1] = 2*wakeDVE[i].eta;  //##!!
				D[row][col-2] = 1.;
				D[row][col-1] = 2*wakeDVE[h].eta;  //##!!
				R[row] 		  = 0.;
			break;

			case 22:	//gamma[i]' = gamma[j]' AND gamma[i]" = gamma[j]"
				h = i - 1;  //index of element to left
				if (h<0)
				{
					printf("\n WARNING!! Left boundary condition of ");
					printf("panel %d messed up!!!",panel+1);
				}

				D[row][col]   = 0.;
				D[row][col+1] = -1.;
				D[row][col-2] = 0.;
				D[row][col-1] = 1.;
				R[row] 		  = 0.;

				row --; //next row
				D[row][col]   = -1.;
				D[row][col+1] = 2*wakeDVE[i].eta;  //##!!
				D[row][col-2] = 1.;
				D[row][col-1] = 2*wakeDVE[h].eta;  //##!!
				R[row] 		  = 0.;
			break;

			default:	//default case assumes gamma = 0
				printf("\nWARNING!! \nEdge 1  boundary condition of");
				printf(" of panel %d undefined!!\n",panel+1);
				exit(0); //exit program
			break;
		}		// end switch statement for edge 1 bound. cond.

//////////////////////////////////////////////////////////////////
		//assembling D within in the boundaries of a panel

		i ++; //advancing index

		for(n=1;n<(panelPtr[panel].n);n++) //loop over spanwise elements-1
		{
			col = i*2;
			row = col;
			h = i - 1;  //index of element to left
			if (h<0)
			{
				printf("\n WARNING!! Problem with  ");
				printf("panel %d !!!",panel+1);
				exit(0);//exiting program
			}
//printf("\n i %d  col = %d row= %d h = %d  ",i,col,row,h);//#

			D[row][col]   = wakeDVE[i].eta;
			D[row][col+1] = -twothirds*wakeDVE[i].eta*wakeDVE[i].eta;
			D[row][col-2] = wakeDVE[h].eta;
			D[row][col-1] = twothirds*wakeDVE[h].eta*wakeDVE[h].eta;
			R[row] 		  = newestDVE[i].K - newestDVE[h].K;

			row --; //next row
			D[row][col]   = -1.;
			D[row][col+1] = 2*wakeDVE[i].eta;  //##!!
			D[row][col-2] = 1.;
			D[row][col-1] = 2*wakeDVE[h].eta;  //##!!
			R[row] 		  = 0.;

			i ++; //advancing index
		}//loop over n, next span element

//////////////////////////////////////////////////////////////////
		//right edge of most right DVE in wake of a panel
		i--;	//correcting index counter
		col = i*2;
		row = col+1;

		switch (panelPtr[panel].BC2)	//edge 2 bound. cond.
		{
			case 100: 	//gamma = 0
				D[row][col]   = wakeDVE[i].eta;
				D[row][col+1] = twothirds*wakeDVE[i].eta*wakeDVE[i].eta;
				R[row] 		  = -newestDVE[i].K;
			break;

			case 110:	//gamma' = 0
				D[row][col]   = 1.;
				D[row][col+1] = 2*wakeDVE[i].eta;
				R[row] 		  = 0.;
			break;

			case 10:	//gamma' = 0
				D[row][col]   = 1.;
				D[row][col+1] = 2*wakeDVE[i].eta;
				R[row] 		  = 0.;
			break;

			case 220:	//gamma[i] = gamma[j] AND gamma[i]' = gamma[j]'
			//is being taken care of by next DVE to the right
			break;

			case 22:	//gamma[i]' = gamma[j]' AND gamma[i]" = gamma[j]"
			//is being taken care of by next DVE to the right
			break;

			default:	//default case assumes gamma = 0
				printf("\nWARNING!! \nEdge 2  boundary condition of");
				printf(" of panel %d undefined!!\n",panel+1);
				exit(0); //exit program
			break;
		}		// end switch statement for edge 1 bound. cond.

		i ++; //advancing index
	}//loop over panels


/*//////////////////
char filename[133];	//file path and name
 FILE *fp;

	//creates file name timestep##.txt ## is number of timestep
	sprintf(filename,"%s%s",OUTPUT_PATH,"test.txt");

	 fp = fopen(filename, "a");

	fprintf(fp, "\n");

	for(row=0;row<size;row++)
 	{
		 //row number
		fprintf(fp, "\n%d\t",row);
		for(col=0;col<size;col++)
		{
			fprintf(fp, "%.10lf\t",D[row][col]);
		}
 		//n-th element of R
		fprintf(fp, "\t\t%.10lf \t%.10lf \t%.10lf \t%.10lf",R[row],wakeDVE[row/2].eta,wakeDVE[row/2].K,newestDVE[row/2].K);
 	}
// 	fclose(fp);

//*/////////////////

	//solving equation system with gaussian elimination
 	GaussSolve(D,R,size,BC);	//subroutine in gauss.cpp

	//assign new coefficients
 	for(i=0;i<info.nospanelement;i++)
 	{
		row = 2*i;
		wakeDVE[i].B = BC[row];
		wakeDVE[i].C = BC[row+1];
		wakeDVE[i].A = newestDVE[i].K - 0.5*twothirds*wakeDVE[i].C*\
					   					wakeDVE[i].eta*wakeDVE[i].eta;
	}

/*//////////////////
//char filename[133];	//file path and name
// FILE *fp;

	//creates file name timestep##.txt ## is number of timestep
//	sprintf(filename,"%s%s",OUTPUT_PATH,"test.txt");
//
//	 fp = fopen(filename, "a");

	fprintf(fp, "\n");

	for(i=0;i<info.nospanelement;i++)
 	{
			fprintf(fp, "A %.10f\tB %.10lf\tC %.10lf\n",wakeDVE[i].A,wakeDVE[i].B,wakeDVE[i].C);
	}
 	fclose(fp);

//*/////////////////

	//Free memory
	FREE2D(&D,size,size);
	FREE1D(&R,size);
	FREE1D(&BC,size);

}
//===================================================================//
		//END FUNCTION New_vorticity_coefficients
//===================================================================//

