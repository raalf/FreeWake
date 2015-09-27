//In the code, the terms "KHH" and "Horstmann" refer to equations
//and methods that are described in:
//"Ein Mehrfach-Traglinienverfahren und seine Verwendung fuer Entwurf
//und Nachrenchnung nichtplanarer Fluegelanordnungen",
//by K.-H. Horstmann,published 1987 by DFVLR, Germany, DFVLR-FB 87-51.
//
//References to program lines are with respect to the Fortran code
//"FLUE", original written by K.-H. Horstmann in 1986 as part of his
//dissertation research
//
//decleration of subroutines
//
//computes induced velocity of system of surface and wake DVE's
void DVE_Induced_Velocity(const GENERAL,const double [3],const DVE *,\
						  DVE **,const int,double [3]);
//computes velocity due to surface DVE's
void Surface_DVE_Vel_Induction(const GENERAL,const double [3],\
							   const DVE *,double [3]);
//computes velocity due to wake DVE's
void Wake_DVE_Vel_Induction(const GENERAL,const double [3],\
							DVE **,const int,double [3]);
//computes ind. velocity in a point due to a DVE with possible symmetry
void Single_DVE_Induced_Velocity(const GENERAL,const DVE,const double [3],\
										double [3],const int);
//computes the influence coefficients in P due a DVE and its possible symmetry
void DVE_Influence_Coeff(const DVE,const GENERAL,const double [3],\
						 double [3],double [3],double [3],const int);
//
//fixed wake model routines
//
//computes induced velocity of complete system
void Induced_Velocity(const BOUND_VORTEX *,const GENERAL,\
	 				  double const [3], double [3]);
//computes inducde velocity of fixed wake along trailing edge
void Fixed_Wake_TEinduction(const BOUND_VORTEX *,const GENERAL,\
	 				  	   double const [3], double [3]);
//computes induced velocity of a vortex line with parabolic circ. distr.
void Vortex_Line_Induced_Velocity(const GENERAL,double const [3],\
								  double const [3],double,double,double,\
								  double,double,double,double [3]);
//computes influence coefficients of bound vortex and its trailing wake
void Influence_Coeff(const BOUND_VORTEX,const GENERAL,double const [3],\
	 				 double [3],double [3],double [3]);
//computes influence coefficients of bound vortex j in point i
void BoundVortexInduction(const double [3],const double [3], \
						  const double, const double, const double,\
						  double [3],double [3],double [3]);
//computes influence coefficients of semi-infinite vortex sheet j in point i
void VortexSheetInduction(const double [3],const double [3], \
						  const double, const double, const double,\
						  double [3],double [3],double [3],const double);
//computes influence coefficients of semi-infinite vortex sheet j in point i
//same as VortexSheetInduction, but original from Horstmann
void FixedWakeInduction(const double [3],const double [3], \
						const double, const double, const double,\
						double [3],double [3],double [3]);

//===================================================================//
		//START function DVE_Induced_Velocity
//===================================================================//
void DVE_Induced_Velocity(const GENERAL info,const double P[3],\
						  const DVE *surfacePtr,DVE **wakePtr,\
						  const int timestep,double w_ind[3])
{
//function computes induced velocity at point P due to the vorticity
//of distributed vorticity elements on the lifting surface and in
//the wake
//input:
// info			general element information
// P			point of where induced velocities are being computed
// surfacePtr	surface DVE's information
// wakePtr		wake DVE's information
// timestep		current time step
//
//output:
// w_ind	induced velocity in point P

	double w_surface[3];	//ind. velocity of surface DVE's
	double w_wake[3];		//ind. velocity of wake DVE's

	//computes velocity that is induced by surface DVE's
	Surface_DVE_Vel_Induction(info,P,surfacePtr,w_surface);

	//computes velocity that is induced by wake DVE's
	Wake_DVE_Vel_Induction(info,P,wakePtr,timestep,w_wake);

//#fprintf(test,"w_s: %2.3lf   %2.3lf  %2.3lf ",w_surface[0],w_surface[1],w_surface[2]);//#
//#fprintf(test,"w_w: %2.3lf   %2.8lf  %2.3lf\n",w_wake[0],w_wake[1],w_wake[2]);//#
//printf(" ws: %2.4lf %2.4lf %2.4lf ",w_surface[0],w_surface[1],w_surface[2]);//#
//printf("ww: %2.4lf %2.4lf %2.4lf",w_wake[0],w_wake[1],w_wake[2]);//#

	//adding of velocity components
	w_ind[0] = w_wake[0] + w_surface[0];
	w_ind[1] = w_wake[1] + w_surface[1];
	w_ind[2] = w_wake[2] + w_surface[2];
//#printf("                w_total  : %lf   %lf   %lf\n",w_ind[0],w_ind[1],w_ind[2]);//#
}
//===================================================================//
		//END function DVE_Induced_Velocity
//===================================================================//

//===================================================================//
		//START function Surface_DVE_Vel_Induction
//===================================================================//
void Surface_DVE_Vel_Induction(const GENERAL info,const double P[3],\
							   const DVE *surfacePtr,double w_surface[3])
{
//this function computes the velocity in point P that is indunced by
//the surface DVEs.
//
//input:
// info			general element information
// P			point of where induced velocities are being computed
// surfacePtr	wake DVE information
//
//output:
// w_surface		induced velocity in point P

	int j; 						//counter over number of surface DVE's
	double w_ind[3];			//delta induced velocities

	//initializing w_surface
	w_surface[0]=0;
	w_surface[1]=0;
	w_surface[2]=0;


	//loop over surface DVEs that induce velocities
	for(j=0;j<info.noelement;j++)
	{
			//computing induced velocity of j-th surface DVE on point P.
			Single_DVE_Induced_Velocity(info,surfacePtr[j],P,w_ind,0);

			//adding delta velocities
			w_surface[0] += w_ind[0];
			w_surface[1] += w_ind[1];
			w_surface[2] += w_ind[2];

//#fprintf(test,"P %2.3lf\t%2.3lf\t%2.3lf",P[0],P[1],P[2]);//#
//#fprintf(test,"  j = %d w_ind %lf\t%lf\t%lf\n",j,w_ind[0],w_ind[1],w_ind[2]);//#
	}

//#fprintf(test,"P %2.3lf\t%2.3lf\t%2.3lf",P[0],P[1],P[2]);//#
//#fprintf(test," w_surface %lf\t%lf\t%lf\n",w_surface[0],w_surface[1],w_surface[2]);//#
//#printf("w_surface %lf\t%lf\t%lf\n",w_surface[0],w_surface[1],w_surface[2]);//#
}
//===================================================================//
		//END function Surface_DVE_Vel_Induction
//===================================================================//

//===================================================================//
		//START function Wake_DVE_Vel_Induction
//===================================================================//
void Wake_DVE_Vel_Induction(const GENERAL info,const double P[3],\
							DVE **wakePtr,const int timestep,\
							double w_wake[3])
{
//this function computes the velocity in point P that is indunced by the
//DVEs in the wake.
//
//input:
// info		general element information
// P		point of where induced velocities are being computed
// wakePtr	wake DVE information
// timestep	current timestep
//
//output:
// w_wake	induced velocity in point P

	int time,span; 				//counter up to current time step, over span
	double w_ind[3];			//delta induced velocities

	//initializing w_wake
	w_wake[0]=0;
	w_wake[1]=0;
	w_wake[2]=0;

	for(span=0;span<info.nospanelement;span++)
	{
		//computing induced velocity at P due to most downstream wake DVEs
		//The very first row of wake DVEs consist of semi-infinite vortex
		//sheets.  At the first time step they have a vortex filametn at
		//their leading edge
		if(timestep<1)
			Single_DVE_Induced_Velocity(info,wakePtr[0][span],P,w_ind,-3);
		else		//changed to simple wake DVE, G.B,. 11-5-06
			Single_DVE_Induced_Velocity(info,wakePtr[0][span],P,w_ind,3);
//			Single_DVE_Induced_Velocity(info,wakePtr[0][span],P,w_ind,1);

		//adding "starting" DVE influence
		w_wake[0] += w_ind[0];
		w_wake[1] += w_ind[1];
		w_wake[2] += w_ind[2];

		//computing ind. vel. at P due to wake DVE just aft of trailing edge
		//These DVEs have a vortex filament only at their leading edge
		if(timestep>=1)
		{
			Single_DVE_Induced_Velocity(info,wakePtr[timestep][span],P,w_ind,2);

			//adding "starting" DVE influence
			w_wake[0] += w_ind[0];
			w_wake[1] += w_ind[1];
			w_wake[2] += w_ind[2];
		}
	}

	//loop over wake DVE's, time and span dependent
	for(time=1;time<timestep;time++)
	for(span=0;span<info.nospanelement;span++)
	{
		//computing induced velocity of wake DVE [time][span] on point P.
		Single_DVE_Induced_Velocity(info,wakePtr[time][span],P,w_ind,1);

		//adding delta velocities
		w_wake[0] += w_ind[0];
		w_wake[1] += w_ind[1];
		w_wake[2] += w_ind[2];
	}
//#printf("w_wake %lf\t%lf\t%lf\n",w_wake[0],w_wake[1],w_wake[2]);//#
}
//===================================================================//
		//END function Wake_DVE_Induction
//===================================================================//

//===================================================================//
		//START function Single_DVE_Induced_Velocity
//===================================================================//
void Single_DVE_Induced_Velocity(const GENERAL info,const DVE DVelement,\
								 const double P[3],double w_ind[3],\
								 const int DVE_type)
{
//computes induced velocity in point P due to DVE DVelement.
//also computes induced velocity of symmetry element
//
//input:
//	info		general info
//	DVelement	distributed vorticity element that induces
//	P			point of interest
// 	DVE_type 	type of DVE,
//	 	DVE_type = 0 DVE has vortex filaments at leading and
//		 		 	trailing edge, usually lifting surface DVE
//		DVE_type = 1 DVE has no vortex filaments at leading and
//					trailing edge, usually wake DVE
//		DVE_type = 2 DVE has vortex filament at leading edge,
//					but not at trailing edge
//		DVE_type =-2 DVE has vortex filament at trailing edge,
//					but not at leading edge
//		DVE_type = 3 DVE is a semi infinite vortex sheet without a
//					vortex filaments at leading and trailing edge
//		DVE_type =-3 DVE is a semi infinite vortex sheet with a
//					 vortex filaments at its leading edge
//		DVE_type = 4 DVE is a vortex sheet that is located from 1/2xsi to
//					 xsi aft of the ref. pt. (for CD computation along TE)
//		DVE_type =-4 DVE is a vortex from -xsi to 1/2xsi aft of the ref. pt.
//					 also vortex filament at LE (CD computation along TE)

//
//
//output:
//	w_ind		induced velocity in P due to DVelement

double a3[3],b3[3],c3[3];	//influence coefficient

	//computes influence coefficients, a3,b3, and c3, of
	//DVE on on point P.
	//(symmetry case is considered in subroutine)
	DVE_Influence_Coeff(DVelement,info,P,a3,b3,c3,DVE_type);

//##########################################
//FILE *fp;														//#
//fp = fopen("output\\test.txt", "a");									//#
//fprintf(fp,"a3: %lf\t%lf\t%lf\tA %lf\n",a3[0],a3[1],a3[2],DVelement.A);//#
//fprintf(fp,"b3: %lf\t%lf\t%lf\tB %lf\n",b3[0],b3[1],b3[2],DVelement.B);//#
//fprintf(fp,"c3: %lf\t%lf\t%lf\tC %lf\n",c3[0],c3[1],c3[2],DVelement.C);//#
//fclose (fp);														//#
//###########################################################################
//#fprintf(test,"a3: %lf\t%lf\t%lf ",a3[0],a3[1],a3[2]);//#
//#fprintf(test,"b3: %lf\t%lf\t%lf ",b3[0],b3[1],b3[2]);//#
//#fprintf(test,"c3: %lf\t%lf\t%lf\n",c3[0],c3[1],c3[2]);//#

	//velocities (sort of like KHH eq. 36 and 38)
	w_ind[0] = DVelement.A*a3[0] + DVelement.B*b3[0] + DVelement.C*c3[0];
	w_ind[1] = DVelement.A*a3[1] + DVelement.B*b3[1] + DVelement.C*c3[1];
	w_ind[2] = DVelement.A*a3[2] + DVelement.B*b3[2] + DVelement.C*c3[2];

	scalar(w_ind,-1/(4*Pi),w_ind);

}
//===================================================================//
		//END function Single_DVE_Induced_Velocity
//===================================================================//

//===================================================================//
		//START function DVE_Influence_Coeff
//===================================================================//

//computes the influence coefficient due to the induction of a
//distributed vortex element at point P. Considers symmetrical conditions.

void DVE_Influence_Coeff(const DVE element,const GENERAL info,\
					 	 const double P[3],double a[3],double b[3],\
					 	 double c[3],const int DVE_type)
{
//input:
	// element	DV element that induces
	// P		point in which DVE induces
	// info		general setting info, such as symmetry
	// DVE_type 	type of DVE,
	//	 DVE_type = 0 DVE has vortex filaments at leading and
	//			 	trailing edge, usually lifting surface DVE
	//	 DVE_type = 1 DVE has no vortex filaments at leading and
	//				trailing edge, usually wake DVE
	//	 DVE_type = 2 DVE has vortex filament at leading edge,
	//				but not at trailing edge
	//	 DVE_type =-2 DVE has vortex filament at trailing edge,
	//				but not at leading edge
	//	 DVE_type = 3 DVE is a semi infinite vortex sheet without a
	//				 vortex filaments at leading and trailing edge
	//	DVE_type =-3 DVE is a semi infinite vortex sheet with a
	//				 vortex filaments at its leading edge
	//	DVE_type = 4 DVE is a vortex sheet that is located from 1/2xsi to
	//			 xsi aft of the ref. pt. (for CD computation along TE)
	//	DVE_type =-4 DVE is a vortex from -xsi to 1/2xsi aft of the ref. pt.
	//				also vortex filament at LE (CD computation along TE)
	//
	//
	//output:
	//a, b, c	influence coefficients as described by KHH in Appendix 3

//computes influence coefficients in point i due to a DVE
void DVEInduction(const double [3],const double [3],const double,\
				  const double,const double,const double,const double,\
				  const double,const double,double [3],double [3],double [3],\
				  const int,const double);

	double a3[3],b3[3],c3[3];		//influence coefficients of DVE
	double tempA[3];

//#printf("P=%lf\t%lf\t%lf\n",P[0],P[1],P[2]);//#

	//computes influence coefficients, a3,b3, and c3 of DVE in P
	DVEInduction(P, element.xo,element.nu,element.epsilon,element.phiLE,\
				 element.phiTE,element.psi,element.eta,element.xsi,a3,b3,c3,\
				 DVE_type,element.singfct);
					 					//subroutine in induced_velocity.cpp

//#printf("a3:  %lf\t%lf\t%lf \n",a3[0],a3[1],a3[2]);//#
//#printf("b3:  %lf\t%lf\t%lf \n",b3[0],b3[1],b3[2]);//#
//#printf("c3:  %lf\t%lf\t%lf \n",c3[0],c3[1],c3[2]);//#

	a[0] = a3[0];		b[0] = b3[0];		c[0] = c3[0];
	a[1] = a3[1];		b[1] = b3[1];		c[1] = c3[1];
	a[2] = a3[2];		b[2] = b3[2];		c[2] = c3[2];

	//symmetrical geometry and flight condition
	if ((info.sym == 1)) //# g.b. 9.17.02 not needed condition: && (info.beta == 0))
	{
//#printf("we are in symmetric mode\n");//#
	//the symmetry condition requires the addition of the induced
	//velocities caused by the mirror images
	//the mirror images require:
	//			1. negation of the y-coodinate of the DVE reference point
	//			2. negation of phi, nu, and psi of DVE (the inducing one)
	//			3. negation of the second coefficient of the vorticty fct

		tempA[0] =  element.xo[0];
		tempA[1] = -element.xo[1];	//negates y component
		tempA[2] =  element.xo[2];
//#printf("xo of mirror image: %lf\t %lf\t %lf \n",tempA[0],tempA[1],tempA[2]);//#

		//computes influence coefficients, a3,b3, and c3 of DVE in P
		DVEInduction(P, tempA,-element.nu,element.epsilon,-element.phiLE,\
					 -element.phiTE,-element.psi,element.eta,element.xsi,\
					 a3,b3,c3,DVE_type,element.singfct);
					 					//subroutine in induced_velocity.cpp

//#printf("a3:  %lf\t%lf\t%lf \n",a3[0],a3[1],a3[2]);//#
//#printf("b3:  %lf\t%lf\t%lf \n",b3[0],b3[1],b3[2]);//#
//#printf("c3:  %lf\t%lf\t%lf \n",c3[0],c3[1],c3[2]);//#


	a[0] += a3[0];		b[0] -= b3[0];		c[0] += c3[0];
	a[1] += a3[1];		b[1] -= b3[1];		c[1] += c3[1];
	a[2] += a3[2];		b[2] -= b3[2];		c[2] += c3[2];

//#printf("a=%lf\t%lf\t%lf\n",a[0],a[1],a[2]);//#
//#printf("b=%lf\t%lf\t%lf\n",b[0],b[1],b[2]);//#
//#printf("c=%lf\t%lf\t%lf\n",c[0],c[1],c[2]);//#
	}//end if statement for symmetrical condition
}
//===================================================================//
		//END function DVE_Influence_Coeff
//===================================================================//

//===================================================================//
		//START function DVEInduction
//===================================================================//
void DVEInduction(const double xA[3],double const xo[3],const double nu,\
				  const double eps,const double phiLE,const double phiTE,\
				  const double psi,const double eta,const double xsi,\
				  double a3x[3],double b3x[3],double c3x[3],\
				  const int DVE_type,const double singfct)
{
//this function computes influence coefficients in point i due to
//a distributed vorticity element.
//The element has two vortex lines that are apart by 2*xsi and have a span
//of 2*eta.  A vortex sheet is inbetween the two vortex lines.  If the
//leading vortex line has a circulation of A+B*eta+C*eta^2 the aft line's
//circulation is the negative value.  The sheet's vorticity is the derivative
//B+2*eta*C.
//The influence coefficient of a DVE is the combination of four separate
//influence coefficients computed according to KHH: two of lifting lines
//and two of semi-infinite vortex sheets, one starting at the leading edge
//of the DVE and one of negative magnitude starting at the trailing edge of
//the DVE.
//the return values are in global coordinates
//
//Input:
//	xA		  	point where velocity is induced in global coordinates
//  xo 	  		reference point of DVE that induces
//  nu,eps		DVE roll, pitch angles
//  phiLE,phiTE	DVE leading and trailing edge sweeps
//	psi			DVE yaw angle
//  eta,xsi		DVE half-span and half-chord
//  DVE_type 	type of DVE,
//			DVE_type = 0 DVE has vortex filaments at leading and
//						 trailing edge, usually lifting surface DVE
//			DVE_type = 1 DVE has no vortex filaments at leading and
//						 trailing edge, usually wake DVE
//			DVE_type = 2 DVE has vortex filament at leading edge,
//						 but not at trailing edge
//			DVE_type =-2 DVE has vortex filament at trailing edge,
//						 but not at leading edge
//			DVE_type = 3 DVE is a semi infinite vortex sheet without a
//						 vortex filaments at leading and trailing edge
//			DVE_type =-3 DVE is a semi infinite vortex sheet with a
//						 vortex filaments at its leading edge
//			DVE_type = 4 DVE is a vortex sheet that is located from 1/2xsi to
//						 xsi aft of the ref. pt. (for CD computation along TE)
//			DVE_type =-4 DVE is a vortex from -xsi to 1/2xsi aft of the ref. pt.
//						 also vortex filament at LE (CD computation along TE)
//  singfct rate at which singularity at edge of vortex sheet decays
//			## singfct added 2/8/05 G.B.  This factor is only needed when
//			relaxing the wake, otherwise it is set to zero (0)
//
//Output:
//	a3x[3]	induced coefficient of DVE in global system
//	b3x[3]	induced coefficient of DVE in global system
//	c3x[3]	induced coefficient of DVE in global system

double rA[3];					//vector betwenn point A and DVE ref. point
double xsiA[3];					//point A coordinates in local DVE system
double xsio[3];					//DVE reference point in local frame
double a1le[3],b1le[3],c1le[3];	//influence coefficients of bound vortex
double a1te[3],b1te[3],c1te[3];	//influence coefficients of bound vortex
double a2le[3],b2le[3],c2le[3];	//influence coefficients of (fixed) wake
double a2te[3],b2te[3],c2te[3];	//influence coefficients of (fixed) wake
double a3xi[3],b3xi[3],c3xi[3];	//influence coefficients of DVE, local

	//Expressing point A with respect to the DVE reference point in the
	//local DVE coordinates
	//1. xA-Xo
	rA[0] = xA[0] - xo[0];
	rA[1] = xA[1] - xo[1];
	rA[2] = xA[2] - xo[2];

	//rotation to loca DVE reference frame, first nu, then epsilon, then psi
	Glob_Star(rA,nu,eps,psi,xsiA);		//function in ref_frame_transform

//#printf("rA %lf %lf %lf xsiA %lf %lf %lf\n",rA[0],rA[1],rA[2],xsiA[0],xsiA[1],xsiA[2]);//##
	//******************************************************* //
	//computing the leading edge influence of the DVE //
	//******************************************************* //

	//vortex system reference point midspan of DVE leading edge,
	//unless when needed for CDiEppler computation. added 8/13/05 G.B.
	if(DVE_type == 4)   xsio[0] = 0.5*xsi;
	else 				xsio[0] = -xsi;

	xsio[1] = 0;
	xsio[2] = 0;

	//computes influence coefficients, a1,b1, and c1, of
	//bound leading edge vortex of element on point xsiA
	//nu=0 since rotation already done

	a1le[0] = 0;  a1le[1] = 0; a1le[2] = 0;
	b1le[0] = 0;  b1le[1] = 0; b1le[2] = 0;
	c1le[0] = 0;  c1le[1] = 0; c1le[2] = 0;
	if(DVE_type == 0 || DVE_type == 2 || DVE_type == -3 || DVE_type == -4)
		BoundVortexInduction(xsiA,xsio,0,phiLE,eta,a1le,b1le,c1le);
							//subroutine in induced_velocity.cpp

//#printf("a1le: %lf\t%lf\t%lf\n",a1le[0],a1le[1],a1le[2]); //#
//#printf("b1le: %lf\t%lf\t%lf\n",b1le[0],b1le[1],b1le[2]); //#
//#printf("c1le: %lf\t%lf\t%lf\n",c1le[0],c1le[1],c1le[2]); //#
//#fprintf(test,"1lez: %lf\t%lf\t%lf ",a1le[2],b1le[2],c1le[2]); //#

	//computes influence coefficients, a2,b2, and c2, of
	//semi-infinite vortex sheet starting at leading edge of element
	//on point xsiA
	//nu=0 since rotation already done
	VortexSheetInduction(xsiA,xsio,0,phiLE,eta,a2le,b2le,c2le,singfct);
								//subroutine in induced_velocity.cpp

//#fprintf(test,"2lez: %lf\t%lf\t%lf ",a2le[2],b2le[2],c2le[2]); //#
//#printf("a2le: %lf\t%lf\t%lf\n",a2le[0],a2le[1],a2le[2]); //#
//#printf("b2le: %lf\t%lf\t%lf\n",b2le[0],b2le[1],b2le[2]); //#
//#printf("c2le: %lf\t%lf\t%lf\n",c2le[0],c2le[1],c2le[2]); //#

	a3xi[0] = a1le[0];
	a3xi[1] = a1le[1];
	a3xi[2] = a1le[2];

	b3xi[0] = (b1le[0]+b2le[0]);	c3xi[0] = (c1le[0]+c2le[0]);
	b3xi[1] = (b1le[1]+b2le[1]);	c3xi[1] = (c1le[1]+c2le[1]);
	b3xi[2] = (b1le[2]+b2le[2]);	c3xi[2] = (c1le[2]+c2le[2]);

	//******************************************************** //
	//computing the trailing edge influence of the DVE element //
	//******************************************************** //

	//vortex system reference point midspan of DVE trailing edge
	//unless when neede for the alternative CDiEppler computation
	if(DVE_type == -4)	xsio[0] = 0.5*xsi; //added 8/13/05 G.B.
	else 				xsio[0] = xsi;

	//computes influence coefficients, a1,b1, and c1, of
	//bound trailing edge vortex of element on point xsiA
	//nu=0 since rotation already done
	a1te[0] = 0;  a1te[1] = 0; a1te[2] = 0;
	b1te[0] = 0;  b1te[1] = 0; b1te[2] = 0;
	c1te[0] = 0;  c1te[1] = 0; c1te[2] = 0;
	if(DVE_type == 0 || DVE_type == -2)
		BoundVortexInduction(xsiA,xsio,0,phiTE,eta,a1te,b1te,c1te);
							//subroutine in induced_velocity.cpp

//#fprintf(test,"1tez: %lf\t%lf\t%lf \n",a1te[2],b1te[2],c1te[2]); //#
//#printf("a1te: %lf\t%lf\t%lf\n",-a1te[0],-a1te[1],-a1te[2]); //#
//#printf("b1te: %lf\t%lf\t%lf\n",-b1te[0],-b1te[1],-b1te[2]); //#
//#printf("c1te: %lf\t%lf\t%lf\n",-c1te[0],-c1te[1],-c1te[2]); //#

	//computes influence coefficients, a2,b2, and c2, of
	//semi-infinite vortex sheet starting at trailing edge of element
	//on point xsiA
	//nu=0 since rotation already done
	a2te[0] = 0;  a2te[1] = 0; a2te[2] = 0;
	b2te[0] = 0;  b2te[1] = 0; b2te[2] = 0;
	c2te[0] = 0;  c2te[1] = 0; c2te[2] = 0;

	if(DVE_type != 3 && DVE_type != -3)
		VortexSheetInduction(xsiA,xsio,0,phiTE,eta,a2te,b2te,c2te,singfct);
								//subroutine in induced_velocity.cpp

//#fprintf(test,"2tez: %lf\t%lf\t%lf \n",a2te[2],b2te[2],c2te[2]); //#
//#printf("b2  : %lf\t%lf\t%lf\n",b2le[0]-b2te[0],b2le[1]-b2te[1],b2le[2]-b2te[2]); //#
//#printf("c2te :  %lf\t%lf\t%lf\n",c2le[0]-c2te[0],c2le[1]-c2te[1],c2le[2]-c2te[2]); //#
//#printf("c2te :  %lf\t%lf\t%lf\n",c2te[0],c2te[1],c2te[2]); //#

	a3xi[0] -= a1te[0];
	a3xi[1] -= a1te[1];
	a3xi[2] -= a1te[2];

	b3xi[0] -= (b1te[0]+b2te[0]);	c3xi[0] -= (c1te[0]+c2te[0]);
	b3xi[1] -= (b1te[1]+b2te[1]);	c3xi[1] -= (c1te[1]+c2te[1]);
	b3xi[2] -= (b1te[2]+b2te[2]);	c3xi[2] -= (c1te[2]+c2te[2]);

	//*********************************** //
	//rotating back into global co-system //
	//*********************************** //
	Star_Glob(a3xi,nu,eps,psi,a3x);		//function in ref_frame_transform
	Star_Glob(b3xi,nu,eps,psi,b3x);		//function in ref_frame_transform
	Star_Glob(c3xi,nu,eps,psi,c3x);		//function in ref_frame_transform
//printf("a: %2.3lf %2.3lf %2.3lf",a3xi[0],a3xi[1],a3xi[2]); //#
//printf(" %2.3lf %2.3lf %2.3lf",a3x[0],a3x[1],a3x[2]); //#
//printf(" b: %2.3lf %2.3lf %2.3lf",b3xi[0],b3xi[1],b3xi[2]); //#
//printf(" %2.3lf %2.3lf %2.3lf\n",b3x[0],b3x[1],b3x[2]); //#
//printf("c: %2.3lf%2.3lf%2.3lf",c3xi[0],c3xi[1],c3xi[2]); //#
//printf(" %2.3lf %2.3lf %2.3lf\n",c3x[0],c3x[1],c3x[2]); //#
//##########################################
//FILE *fp;														//#
//fp = fopen("output\\test.txt", "a");									//#
//fprintf(fp,"a: %2.3lf %2.3lf %2.3lf\n",a3xi[0],a3xi[1],a3xi[2]); //#
//fprintf(fp," %2.3lf %2.3lf %2.3lf",a3x[0],a3x[1],a3x[2]); //#
//fprintf(fp," b: %2.3lf %2.3lf %2.3lf\n",b3xi[0],b3xi[1],b3xi[2]); //#
//fprintf(fp,"b2le: %2.3lf%2.3lf%2.3lf ",b2le[0],b2le[1],b2le[2]); //#
//fprintf(fp,"b2te: %2.3lf%2.3lf%2.3lf\n",b2te[0],b2te[1],b2te[2]); //#
//fprintf(fp,"b1le: %2.3lf%2.3lf%2.3lf ",b1le[0],b1le[1],b1le[2]); //#
//fprintf(fp,"b1te: %2.3lf%2.3lf%2.3lf\n",b1te[0],b1te[1],b1te[2]); //#
//fprintf(fp," %2.3lf %2.3lf %2.3lf\n",b3x[0],b3x[1],b3x[2]); //#
//fprintf(fp,"c: %2.3lf%2.3lf%2.3lf\n",c3xi[0],c3xi[1],c3xi[2]); //#
//fprintf(fp," %2.3lf %2.3lf %2.3lf\n",c3x[0],c3x[1],c3x[2]); //#
//fprintf(fp,"c2le: %2.3lf%2.3lf%2.3lf ",c2le[0],c2le[1],c2le[2]); //#
//fprintf(fp,"c2te: %2.3lf%2.3lf%2.3lf\n",c2te[0],c2te[1],c2te[2]); //#
//fprintf(fp,"c1le: %2.3lf%2.3lf%2.3lf ",c1le[0],c1le[1],c1le[2]); //#
//fprintf(fp,"c1te: %2.3lf%2.3lf%2.3lf\n",c1te[0],c1te[1],c1te[2]); //#
//fclose (fp);														//#
//###########################################################################

}
//===================================================================//
		//END function DVEInduction
//===================================================================//

//===================================================================//
		//START function Induced_Velocity
//===================================================================//
void Induced_Velocity(const BOUND_VORTEX *elementPtr,const GENERAL info,\
	 				  double const P[3], double w_ind[3])
{
	//function computes induced velocity at point P due to wake and bound
	//vorticity
	//input:
	//			general element information
	// P		point of where induced velocities are being computed
	//
	//output:
	// w_ind	induced velocity in point P

	int j; 					//counter of elements that induce velocity on P
	double a[3],b[3],c[3];	//influence coefficients
	double A, B, C;			//vorticity coefficients of inducing element

	//initializing w_ind
	w_ind[0]=0;
	w_ind[1]=0;
	w_ind[2]=0;

	for(j=0;j<info.noelement;j++)
	{
		A = elementPtr[j].A;
		B = elementPtr[j].B;
		C = elementPtr[j].C;

		//computes influence coefficients of element j on point P
		//according to KHH Appendix 3, incl. fixed wake influence
		Influence_Coeff(elementPtr[j],info,P,a,b,c);
									//subroutine in induced_velocity.cpp

//#printf("a=%lf\t%lf\t%lf\n",a[0],a[1],a[2]);//#
//#printf("b=%lf\t%lf\t%lf\n",b[0],b[1],b[2]);//#
//#printf("c=%lf\t%lf\t%lf\n",c[0],c[1],c[2]);//#

		//KHH eq. 36 and 38
		w_ind[0] += A*a[0]+B*b[0]+C*c[0];
		w_ind[1] += A*a[1]+B*b[1]+C*c[1];
		w_ind[2] += A*a[2]+B*b[2]+C*c[2];
	}
	scalar(w_ind,-1/(4*Pi),w_ind);
//#printf("w_ind %l\t%lf\t%lf\n",w_ind[0],w_ind[1],w_ind[2]);//#
}
//===================================================================//
		//END function Induced_Velocity
//===================================================================//

//===================================================================//
		//START function Fixed_Wake_TEinduction
//===================================================================//
void Fixed_Wake_TEinduction(const BOUND_VORTEX *trailedgePtr,\
	 					    const GENERAL info,\
	 				        double const P[3], double w_ind[3])
{
	//function computes velocity in point P that is indunced by fixed
	//wake, P is part of trailing edge.  The computed velocity is needed
	//to determine the induced drag according to Eppler along the trailing
	//edge due to a fixed wake.  According to Schmid-Goeller, the ind. drag
	//depends only on the influence of a fixed wake that is without sweep
	//and starts in P[0].
	//
	//input:
	//			general element information
	// P		point of where induced velocities are being computed
	//
	//a, b, c	influence coefficients as described by KHH in Appendix 3
	//
	//output:
	// w_ind	induced velocity in point P

int j; 						//counter of elements that induce velocity on P
double a[3],b[3],c[3];		//influence coefficients
double a2[3],b2[3],c2[3];	//influence coefficients of (fixed) wake
double tempA[3];

//initializing w_ind
w_ind[0]=0;
w_ind[1]=0;
w_ind[2]=0;

for(j=0;j<info.nospanelement;j++)
{
		//computes influence coefficients of element j on point P
		//according to KHH Appendix 3, incl. fixed wake influence

	tempA[0] = P[0];
	tempA[1] = trailedgePtr[j].xo[1];
	tempA[2] = trailedgePtr[j].xo[2];

	//computes influence coefficients, a2,b2, and c2, of
	//trailing vortex sheet (fixed) of element j on point P
VortexSheetInduction(P,tempA,trailedgePtr[j].nu,0,\
					   trailedgePtr[j].eta,a2,b2,c2,0);

//####CHANGED 7/5/04  G.B.
//	FixedWakeInduction(P,tempA,trailedgePtr[j].nu,0,\
//					   trailedgePtr[j].eta,a2,b2,c2);
								//subroutine in induced_velocity.cpp

	a[0] = a2[0];	a[1] = a2[1];	a[2] = a2[2];	//a2 = 0 for wakes!!
	b[0] = b2[0];	b[1] = b2[1];	b[2] = b2[2];
	c[0] = c2[0];	c[1] = c2[1];	c[2] = c2[2];

//#printf("b2:  %lf\t%lf\t%lf \n",b2[0],b2[1],b2[2]);//#
//#printf("c2:  %lf\t%lf\t%lf \n",c2[0],c2[1],c2[2]);//#

	//symmetrical geometry and flight condition
	if ((info.sym == 1)) //# g.b. 9.17.02 not needed condition: && (info.beta == 0))
	{
	//the symmetry condition requires the addition of the induced
	//velocities caused by the mirror images
	//the mirror images require:
	//			1. negation of the y-coodinate of the inducing element
	//			2. negation of the second coefficient of the vorticty fct
	//			3. negation of phi and nu of element j (the inducing one)

		tempA[1] = -tempA[1];	//negates y component

		//computes influence coefficients, a2,b2, and c2, of
		//trailing vortex sheet (fixed) of element j on point P
		VortexSheetInduction(P,tempA,-trailedgePtr[j].nu,0,\
						   trailedgePtr[j].eta,a2,b2,c2,0);
//####CHANGED 7/5/04  G.B.
//		FixedWakeInduction(P,tempA,-trailedgePtr[j].nu,0,\
//						   trailedgePtr[j].eta,a2,b2,c2);
								//subroutine in induced_velocity.cpp

//#		vsum(a,a2,a);			//a2 = 0 for wakes !!!!
		scalar(b2,-1,tempA);  //mirror image!!
		vsum(b,tempA,b);
		vsum(c,c2,c);
	}//end if statement for symmetrical condition


		//KHH eq. 36 and 38
		w_ind[0] += trailedgePtr[j].B*b[0]+trailedgePtr[j].C*c[0];
		w_ind[1] += trailedgePtr[j].B*b[1]+trailedgePtr[j].C*c[1];
		w_ind[2] += trailedgePtr[j].B*b[2]+trailedgePtr[j].C*c[2];
}
scalar(w_ind,-1/(4*Pi),w_ind);
//#printf("w_ind %lf\t%lf\t%lf\n",w_ind[0],w_ind[1],w_ind[2]);//#
}
//===================================================================//
		//END function Fixed_Wake_TEinduction
//===================================================================//

//===================================================================//
		//START function Vortex_Line_Induced_Velocity
//===================================================================//
void Vortex_Line_Induced_Velocity(const GENERAL info,double const P[3],\
								  double const xo[3],double A,double B,\
								  double C,double eta, double phi,double nu,\
								  double w_ind[3])
{
	//function computes the velocity in point P that is indunced by a
	//vortex line with parabolic circulation distribution.The function
	//considers any possible mirrow-vortex line due to the symmetry
	//condition (sym=1)
	//
	//input:
	// info		general information
	// xo		vortex line reference point
	// A,B,C	vortex line coefficient of circulation
	// eta		vortex line half width
	// phi		vortex line sweep
	// nu		vortex line dihedral
	// P		point of where induced velocities are being computed
	//
	//output:
	// w_ind	induced velocity in point P


double a1[3],b1[3],c1[3];	//influence coefficients
double a[3],b[3],c[3];		//influence coefficients
double tempA[3];

			//computes influence coefficients, a1,b1, and c1, of
			BoundVortexInduction(P,xo,nu,phi,eta,a1,b1,c1);
										//subroutine in induced_velocity.cpp

			//symmetrical geometry and flight condition
			if ((info.sym == 1)) //# g.b. 9.17.02 not needed condition: && (info.beta == 0))
			{
			//the symmetry condition requires the addition of the induced
			//velocities caused by the mirror images
			//the mirror images require:
			//		1. negation of the y-coodinate of the inducing element
			//		2. negation of phi and nu of element j (the inducing one)
			//		3. negation of the second coefficient of the vorticty fct

				tempA[0] =  xo[0];
				tempA[1] = -xo[1];	//negates y component
				tempA[2] =  xo[2];

				//computes influence coefficients, a1,b1, and c1, of
				//bound vortex j on point P
				BoundVortexInduction(P,tempA,-nu,-phi,eta,a,b,c);
										//subroutine in induced_velocity.cpp

				//adding mirror image influence coefficients
				a1[0] += a[0];		b1[0] -= b[0];		c1[0] += c[0];
				a1[1] += a[1];		b1[1] -= b[1];		c1[1] += c[1];
				a1[2] += a[2];		b1[2] -= b[2];		c1[2] += c[2];

			}//end if statement for symmetrical condition

			w_ind[0] = A*a1[0] + B*b1[0] + C*c1[0];
			w_ind[1] = A*a1[1] + B*b1[1] + C*c1[1];
			w_ind[2] = A*a1[2] + B*b1[2] + C*c1[2];
//#printf("%lf\t%lf\t%lf \n",norm2(a1),norm2(b1),norm2(c1));//#

			scalar(w_ind,-1/(4*Pi),w_ind);

//#printf("wbv  %lf   %lf   %lf\n",w_ind[0],w_ind[1],w_ind[2]);  //#
}
//===================================================================//
		//END function Vortex_Line_Induced_Velocity
//===================================================================//

//===================================================================//
		//START function Influence_Coefficient
//===================================================================//

//computes the influence coefficient of vortex system j on point P
//A fixed wake is assumed
void Influence_Coeff(const BOUND_VORTEX element_J,const GENERAL info,
					 double const P[3],double a[3],double b[3],double c[3])
{
	//input:
	// 			general element information
	// P		point on which system is inducing
	//

	//output:
	//a, b, c	influence coefficients as described by KHH in Appendix 3

	double a1[3],b1[3],c1[3];		//influence coefficients of bound vortex
	double a2[3],b2[3],c2[3];		//influence coefficients of (fixed) wake
	double tempA[3];

//#printf("P=%lf\t%lf\t%lf\n",P[0],P[1],P[2]);//#

	//computes influence coefficients, a1,b1, and c1, of
	//bound vortex j on point P
	BoundVortexInduction(P, element_J.xo, \
						 element_J.nu,element_J.phi, \
						 element_J.eta,a1,b1,c1);
							//subroutine in induced_velocity.cpp

//#printf("a1:  %lf\t%lf\t%lf \n",a1[0],a1[1],a1[2]);//#
//#printf("b1:  %lf\t%lf\t%lf \n",b1[0],b1[1],b1[2]);//#
//#printf("c1:  %lf\t%lf\t%lf \n",c1[0],c1[1],c1[2]);//#

	//computes influence coefficients, a2,b2, and c2, of
	//trailing vortex sheet (fixed) of element j on point P
		VortexSheetInduction(P,element_J.xo, \
					   element_J.nu,element_J.phi, \
					   element_J.eta,a2,b2,c2,0);
//####CHANGED 7/5/04  G.B.
//	FixedWakeInduction(P, element_J.xo, \
//					   element_J.nu,element_J.phi, \
//					   element_J.eta,a2,b2,c2);
								//subroutine in induced_velocity.cpp

//#printf("a2:  %lf\t%lf\t%lf \n",a2[0],a2[1],a2[2]);//#
//#printf("b2:  %lf\t%lf\t%lf \n",b2[0],b2[1],b2[2]);//#
//#printf("c2:  %lf\t%lf\t%lf \n",c2[0],c2[1],c2[2]);//#

	//Horstmann Eq.A3-1
	vsum(a1,a2,a);
	vsum(b1,b2,b);
	vsum(c1,c2,c);

	//symmetrical geometry and flight condition
	if ((info.sym == 1)) //# g.b. 9.17.02 not needed condition: && (info.beta == 0))
	{
//#printf("we are in symmetric mode\n");//#
	//the symmetry condition requires the addition of the induced
	//velocities caused by the mirror images
	//the mirror images require:
	//			1. negation of the y-coodinate of the inducing element
	//			2. negation of the second coefficient of the vorticty fct
	//			3. negation of phi and nu of element j (the inducing one)

		tempA[0] = element_J.xo[0];
		tempA[1] = -element_J.xo[1];	//negates y component
		tempA[2] = element_J.xo[2];
//#printf("xo of mirror image: %lf\t %lf\t %lf \n",tempA[0],tempA[1],tempA[2]);//#

		//computes influence coefficients, a1,b1, and c1, of
		//bound vortex j on point P
		BoundVortexInduction(P, tempA, \
							 -element_J.nu,-element_J.phi, \
							 element_J.eta,a1,b1,c1);
							//subroutine in induced_velocity.cpp

//#printf("a1:  %lf\t%lf\t%lf \n",a1[0],a1[1],a1[2]);//#
//#printf("b1:  %lf\t%lf\t%lf \n",b1[0],b1[1],b1[2]);//#
//#printf("c1:  %lf\t%lf\t%lf \n",c1[0],c1[1],c1[2]);//#

		//computes influence coefficients, a2,b2, and c2, of
		//trailing vortex sheet (fixed) of element j on point P
		VortexSheetInduction(P, tempA, \
						   -element_J.nu,-element_J.phi, \
						   element_J.eta,a2,b2,c2,0);
//####CHANGED 7/5/04  G.B.
//		FixedWakeInduction(P, tempA, \
//						   -element_J.nu,-element_J.phi, \
//						   element_J.eta,a2,b2,c2);
								//subroutine in induced_velocity.cpp

//#printf("a2:  %lf\t%lf\t%lf \n",a2[0],a2[1],a2[2]);//#
//#printf("b2:  %lf\t%lf\t%lf \n",b2[0],b2[1],b2[2]);//#
//#printf("c2:  %lf\t%lf\t%lf \n",c2[0],c2[1],c2[2]);//#

			//Horstmann Eq.A3-1
		vsum(a1,a2,tempA);
		vsum(a,tempA,a);
		vsum(b1,b2,tempA);
		scalar(tempA,-1,tempA);  //mirror image!!
		vsum(b,tempA,b);
		vsum(c1,c2,tempA);
		vsum(c,tempA,c);

//#printf("a=%lf\t%lf\t%lf\n",a[0],a[1],a[2]);//#
//#printf("b=%lf\t%lf\t%lf\n",b[0],b[1],b[2]);//#
//#printf("c=%lf\t%lf\t%lf\n",c[0],c[1],c[2]);//#

	}//end if statement for symmetrical condition
}
//===================================================================//
		//END function Influence_Coefficient
//===================================================================//

//===================================================================//
		//START function BoundVortexInduction
//===================================================================//
void BoundVortexInduction(const double xAi[3],const double xoj[3], \
						  const double nuj, const double phij, const double etaj,\
						  double a1x[3],double b1x[3],double c1x[3])
{
//this function computes influence coefficients of ctr point i due to
//the bound vortex at elementary wing j.
//See also Appendix 3 of
//"Ein Mehrfach-Traglinienverfahren und seine Verwendung fuer Entwurf
//und  Nachrenchnung nichtplanarer Fluegelanordnungen",
//by K.-H. Horstmann, published 1987 by DFVLR, Germany, DFVLR-FB 87-51.

//the return values are in global coordinates
//
//Input:
//	xAi		control point i in global coordinates
//	x0j		element j midspan point of 1/4c line in global coordinates
//	nuj		element j dihedral
//	phij	element j sweep
//	etaj	element j half span
//
//Output:
//	a1x[3]	induced coefficient of bound vortex in global system
//	b1x[3]	induced coefficient of bound vortex in global system
//	c1x[3]	induced coefficient of bound vortex in global system
const double ZERO=1e-5;				//definition of zero, larger for filaments

double tempA[3],tempAA[3],tempS;	//temporary arrays, scalar
double xsiA, etaA,zetaA,D;			//ctrl pt i coordinates in j
									//element ref. frame
double a1,b1,c1, reta1,reta2;
double G11, G12, G13;
double a1xi[3],b1xi[3],c1xi[3];		//coefficients in local ref. frame

	//transformation of i-control point global coordinates into
	//local j-element coordinates. Horstmann Eq. 9
	scalar(xoj,-1,tempAA);
	vsum(tempAA,xAi,tempA);
	if(nuj*nuj > DBL_EPS)
	{
		rotateX(tempA,nuj,tempAA);
		xsiA	= tempAA[0];
		etaA	= tempAA[1];
		zetaA	= tempAA[2];
	}
	else
	{
		xsiA	= tempA[0];
		etaA	= tempA[1];
		zetaA	= tempA[2];
	}

	D 		= xsiA-etaA*tan(phij);
	//if D=0 then A in plane of bound vortex and zeta axis

//#fprintf(test," \nxsiA %lf  etaA %lf  zetaA %.15lf  D %.15lf ",xsiA,etaA,zetaA,D); //#
//#fprintf(test," etaj %lf phi %lf  nuj  %lf\n",etaj,phij*RtD,nuj*RtD);

	if(fabs(zetaA)<ZERO && fabs(D)<ZERO)
	{
		//1. condition: point A lies in xsiA-etaA-plane
		//	=>only velocities in zeta-direction induced
		//2. condition: point A lies in plane defined by zeta-axis and
		// 	 bound vortex line
		//  =>no velocity in zeta-direction induced
		//or in other words, point lies on line of vortex filament

		//KHH FLUE line 1366 and 2653
		a1xi[0] =	0;
		a1xi[1] =	0;
		a1xi[2] =	0;

		b1xi[0] =	0;
		b1xi[1] =	0;
		b1xi[2] =	0;

		c1xi[0] =	0;
		c1xi[1] =	0;
		c1xi[2] =	0;
	}
	else
	{
		//Horstmann Eq. A3-6
		a1=1+tan(phij)*tan(phij);
		b1=-(etaA+xsiA*tan(phij));
		c1=xsiA*xsiA+etaA*etaA+zetaA*zetaA;
//#printf("a1:%lf\tb1:%lf\tc1:%lf \n ",a1,b1,c1);//#

		reta1=sqrt(etaj*etaj*a1-2*etaj*b1+c1);
		reta2=sqrt(etaj*etaj*a1+2*etaj*b1+c1);
//#printf("reta1:%lf\t reta2:%lf\n ",reta1,reta2);//#

		//Horstmann Eq. A3-3,4,5
		tempS=	((a1*c1-b1*b1)*reta1*reta2);
		G11	 = 	(a1*etaj*(reta1+reta2)+b1*(reta1-reta2))/tempS;
		G12	 = 	(-b1*etaj*(reta1+reta2)-c1*(reta1-reta2))/tempS;
		G13	 = 	((2*b1*b1-a1*c1)*etaj*(reta1+reta2)+\
				b1*c1*(reta1-reta2))/(a1*tempS);
		G13 += 	log((sqrt(a1)*reta2+a1*etaj+b1)/\
				(sqrt(a1)*reta1-a1*etaj+b1))/sqrt(a1*a1*a1);
//#printf("tempS:%lf\n",tempS);//#
//#printf("G11:%lf\tG12:%lf\tG13:%lf \n ",G11,G12,G13);//#

		//Horstmann Eq. A3-2
		if(fabs(zetaA)<ZERO)
		{
		//condition: point A lies in xsiA-etaA-plane
		//	=>only velocities in zeta-direction induced
			a1xi[0] =	0;
			a1xi[1] =	0;

			b1xi[0] =	0;
			b1xi[1] =	0;

			c1xi[0] =	0;
			c1xi[1] =	0;
		}
		else
		{
			a1xi[0] =	-G11*zetaA;
			a1xi[1] =	-a1xi[0]*tan(phij);

			b1xi[0] =	a1xi[0] * G12/G11;
			b1xi[1] =	a1xi[1] * G12/G11;

			c1xi[0] =	a1xi[0] * G13/G11;
			c1xi[1] =	a1xi[1] * G13/G11;
		}

		if(fabs(D)<ZERO)
		{
		//condition: point A lies in plane defined by zeta-axis and
		//bound vortex line
		// =>no velocity in zeta-direction induced
			a1xi[2] =	0;
			b1xi[2] =	0;
			c1xi[2] =	0;
		}
		else
		{
			a1xi[2] =	G11*D;
			b1xi[2] =	a1xi[2] * G12/G11;
			c1xi[2] =	a1xi[2] * G13/G11;
		}
	} // end else line 223 and if line 202

//#printf("a1xi:%lf  %lf  %lf\n",a1xi[0],a1xi[1],a1xi[2]);//#
//#printf("b1xi:%lf  %lf  %lf\n",b1xi[0],b1xi[1],b1xi[2]);//#
//#printf("c1xi:%lf  %lf  %lf\n",c1xi[0],c1xi[1],c1xi[2]);//#

	//transform coefficients a, b, and c back into global co-system
	if(nuj*nuj > DBL_EPS)
	{
		rotateX(a1xi,-nuj,a1x);
		rotateX(b1xi,-nuj,b1x);
		rotateX(c1xi,-nuj,c1x);
	}
	else
	{
		a1x[0] = a1xi[0]; a1x[1] = a1xi[1]; a1x[2] = a1xi[2];
		b1x[0] = b1xi[0]; b1x[1] = b1xi[1]; b1x[2] = b1xi[2];
		c1x[0] = c1xi[0]; c1x[1] = c1xi[1]; c1x[2] = c1xi[2];
	}
}

//===================================================================//
		//END function BoundVortexInduction
//===================================================================//

//===================================================================//
		//START function VortexSheetInduction
//===================================================================//
//*///&&&&&&&&&&&&&&
void VortexSheetInduction(const double xAi[3],const double xoj[3], \
	 					 const double nuj,const double phij,const double etaj,\
						 double a2x[3],double b2x[3],double c2x[3],\
						 const double singfct)
{
//this function computes influence coefficients of ctr point i due to
//a semi-infinite vortex sheet.  The function is based on the function
//FixedWakeInduction that is used in Horstmann's multiple lifting line
//method.  The curren method has been modified to deal with the singularities
//along the edges of the sheet, but is otherwise very similar to the function
//FixedWakeInduction that is listed further below.
//See also Appendix 3 of
//"Ein Mehrfach-Traglinienverfahren und seine Verwendung fuer Entwurf
//und  Nachrenchnung nichtplanarer Fluegelanordnungen",
//by K.-H. Horstmann, 1987, DFVLR, Germany, DFVLR-FB 87-51.
//
//the return values are in global coordinates
//
//Input:
//	xAi		element i control point in global coordinates
//	x0j		element j midspan point of 1/4c line in global coordinates
//	nuj		element j dihedral
//	phij	element j sweep
//	etaj	element j half span
//  singfct rate at which singularity at edge of vortex sheet decays
//			## singfct added 2/8/05 G.B.  This factor is only needed when
//			relaxing the wake, otherwise it is set to zero (0)
//
//Output:
//	a2x[3]	induced coefficient of bound vortex in global system
//	b2x[3]	induced coefficient of bound vortex in global system
//	c2x[3]	induced coefficient of bound vortex in global system

double tempA[3],tempAA[3],tempS, tempSS;
double xsiA,etaA,zetaA,D;	//ctrl pt i coordinates in j element ref. frame
double ABSzetaA,ABSD;			//absolute values
double zetaASQR;				//zetaA*zetaA
double a2,b2,c2,t1,t2,rt1,rt2,b[7],c[7],G2[7];
double G21b,G21c,G25b,G25c;
double beta1,beta2,delta1,delta2,gamma1,gamma2,rho,epsilon;
double mu1t1,mu1t2,mu2t1,mu2t2,mu3t1,mu3t2;
double b2xi[3],c2xi[3];			//coefficients in local ref. frame
double xsij = fabs(xoj[0]);		//half-chord of DVE is used for scaling

//check if singfct has been assigned
if(singfct*singfct >1000)
{
	printf(" hello there might be a problem with singfct");
	printf(" in VortexSheetInduction in induced_velocity.cpp");
	printf(" singfct = %lf  ",singfct);
//	exit(0);
}

  //transformation of i-control point global coordinates into
  //local j-element coordinates. Horstmann Eq. 9
  scalar(xoj,-1,tempAA);
  vsum(tempAA,xAi,tempA);
  if(nuj*nuj > DBL_EPS)
  {
	  rotateX(tempA,nuj,tempAA);
  	  xsiA		= tempAA[0];
	  etaA		= tempAA[1];
	  zetaA		= tempAA[2];
  }
  else
  {
	  xsiA	= tempA[0];
	  etaA	= tempA[1];
	  zetaA	= tempA[2];
  }

  zetaASQR	= zetaA*zetaA;
  if(zetaASQR <= DBL_EPS) ABSzetaA = 0;
  else					  ABSzetaA 	= sqrt(zetaASQR);
  D 		= xsiA-etaA*tan(phij); //if D=0 then A in plane of bound vortex
  ABSD		= fabs(D);			  //and zeta axis


//##########################################
//FILE *fp;														//#
//fp = fopen("output\\test.txt", "a");									//#
//fprintf(fp,"xAi: %lf  %lf  %.15lf ",xAi[0],xAi[1],xAi[2]);		//#
//fprintf(fp,"xoj: %lf  %lf  %.15lf \n",xoj[0],xoj[1],xoj[2]);		//#
//fprintf(fp,"xsiA: %lf  etaA: %lf  zetaA: %.15lf \n",xsiA,etaA,zetaA);	//#
//fprintf(fp,"phij: %.15lf  nuj: %.15lf  etaj: %lf \n",phij*RtD,nuj*RtD,etaj);	//#
//fprintf(fp,"|zeta|: %.15lf  zeta^2: %.15lf  \n",ABSzetaA,zetaASQR);	//#
//fclose (fp);														//#
//###########################################################################

 // if(ABSzetaA<=DBL_EPS && ABSD<=DBL_EPS)
  if(ABSzetaA<=DBL_EPS && ABSD<=DBL_EPS && fabs(phij)<=DBL_EPS)
  { //point A lies on bound vortex line/leading edge of vortexsheet

	b2x[0] = 0;
	b2x[1] = 0;
	b2x[2] = 0;

	c2x[0] = 0;
	c2x[1] = 0;
	c2x[2] = 0;

	if(fabs(phij)<=DBL_EPS)	//at leading edge of an unswept vortex sheet
	{
		b2xi[0] = 0;
		b2xi[1] = 0;

		c2xi[0] = 0;
		c2xi[1] = 0;

		t1 = etaA+etaj;
		t2 = etaA-etaj;

		//Will be singular if t1 or t2 is 0,
		//correction below requires that gamma=0 at tip!!!
		//Thus, old version, changed 11/29/04
		//old tempS = log((t1*t1)/(t2*t2));
		tempS = log((t1*t1+singfct)/(t2*t2+singfct));

		b2xi[2] = tempS*0.5;
		c2xi[2] = (-4*etaj+etaA*tempS);

		//transform coefficients a, b, and c back into global co-system
		if(nuj*nuj > DBL_EPS)
		{	rotateX(b2xi,-nuj,b2x);
			rotateX(c2xi,-nuj,c2x);
		}
		else
		{
			b2x[0] = b2xi[0]; b2x[1] = b2xi[1]; b2x[2] = b2xi[2];
			c2x[0] = c2xi[0]; c2x[1] = c2xi[1]; c2x[2] = c2xi[2];
		}
	}
  }
  else //if(ABSzetaA<=DBL_EPS && ABSD<=DBL_EPS && fabs(phij)<=DBL_EPS)
  {//point of interest does not fall on leading edge of sheet

	//Horstmann Eq. A3-17
	a2		= 1+tan(phij)*tan(phij);
	b2		= D*tan(phij);
	c2		= D*D+zetaASQR;

	t1		= etaA+etaj;
	t2		= etaA-etaj;

	rt1		= sqrt(t1*t1*a2+2*t1*b2+c2);
	rt2		= sqrt(t2*t2*a2+2*t2*b2+c2);

	//Horstmann Eq. A3-18, indices off by 1!  Horstmann 1 is 0 here
	b[0]	= -D;
	b[1]	= zetaASQR*tan(phij);
	b[2]	= 0;
	b[3]	= -tan(phij);
	b[4]	= -1;
	b[5]	= 0;
	b[6]	= 0;

	c[0]	= -2*(b[1]-etaA*b[0]);
	c[2]	= 2*tan(phij);
	c[3]	= 2*(xsiA-etaA*c[2]);
	c[4]	= -2*etaA;
	c[5]	= -2*zetaASQR;
	c[6]	= 2;
	c[1]	= 0.5*c[5]*c[3];

	//Horstmann Eq. A3-10
	epsilon = b[0]*b[0]+b[1]*b[3];
	rho  	= sqrt(epsilon*epsilon+4*zetaASQR*b2*b2);

	tempS 	= 0.5*(rho+epsilon);
	if (tempS<=DBL_EPS) beta1   = 0;
	else 				beta1	= -sqrt(tempS);
	tempS 	= 0.5*(rho-epsilon);
	if (tempS<=DBL_EPS) beta2   = 0;
	else 				beta2	= -sqrt(tempS);
	//see Horstmann program FLU lines 1443 through 1447
 	tempS=zetaA*b2;
	if (fabs(tempS) > DBL_EPS) beta2*=tempS/fabs(tempS);


//#####
	//changed 11/22/04 G.B.
	//old,	gamma1 	= (a2*beta2*zetaA+b2*beta1)/rho;
	//old,	gamma2 	= (a2*beta1*zetaA-b2*beta2)/rho;
	//old,	delta1 	= (b2*beta2*zetaA+c2*beta1)/rho;
	//old,	delta2 	= (b2*beta1*zetaA-c2*beta2)/rho;
	//old,	tempS	= gamma1*t1+delta1-rt1;
	//old,	tempSS	= gamma2*t1+delta2;

	gamma1 	= (a2*beta2*zetaA+b2*beta1);
	gamma2 	= (a2*beta1*zetaA-b2*beta2);

	delta1 	= (b2*beta2*zetaA+c2*beta1);
	delta2 	= (b2*beta1*zetaA-c2*beta2);

	tempS	= gamma1*t1+delta1-rt1*rho;
	tempSS	= gamma2*t1+delta2;

	//changed 4/13/04 G.B. ##
	//old, mu1t1	= (tempS*tempS+tempSS*tempSS)/(t1*t1+zetaA*zetaA);
	mu1t1	= (tempS*tempS+tempSS*tempSS);

	//mu2t1 computed as in Horstmann Eq. A3-10. atan-values then corrected.
	if (t1*t1 <= DBL_EPS)
	mu2t1 = Pi/2*ABSzetaA/zetaA+atan(tempSS/tempS);//|zeta/t1|->infinity
	else
	{
		mu2t1 = atan(zetaA/t1)+atan(tempSS/tempS);
	//Corrections according to Horstmamm program FLU lines 1473 and 1474
		if ((zetaA > 0.0) && (t1 < 0.0))	mu2t1+=Pi;
		if ((zetaA < 0.0) && (t1 < 0.0))	mu2t1-=Pi;
	}
	//Corrections according to Horstmamm program FLU lines 1484 and 1485
	if (tempS < 0.0)						mu2t1+=Pi;
	if ((tempSS < 0.0) && (tempS > 0.0))	mu2t1+=2*Pi;

//#####
	//changed 11/22/04 G.B.
	//old,	tempS	= gamma1*t2+delta1-rt2;
	tempS	= gamma1*t2+delta1-rt2*rho;
	tempSS	= gamma2*t2+delta2;

	//changed 4/13/04 G.B. ##
	//old, mu1t2	= (tempS*tempS+tempSS*tempSS)/(t2*t2+zetaA*zetaA);
	mu1t2	= (tempS*tempS+tempSS*tempSS);

	//mu2t2 computed as in Horstmann Eq. A3-10. atan-values then corrected.
	if (t2*t2 < DBL_EPS)
	mu2t2 = Pi/2*ABSzetaA/zetaA+atan(tempSS/tempS);//|zeta/t2|->infinity
	else
	{
		mu2t2 = atan(zetaA/t2)+atan(tempSS/tempS);
	//Corrections according to Horstmamm program FLU lines 1466 1467
		if ((zetaA > 0.0) && (t2 < 0.0))	mu2t2+=Pi;
		if ((zetaA < 0.0) && (t2 < 0.0))	mu2t2-=Pi;
	}
	//Corrections according to Horstmamm program FLU lines 1479 and 1480
	if (tempS < 0.0)						mu2t2+=Pi;
	if ((tempSS < 0.0) && (tempS > 0.0))	mu2t2+=2*Pi;

//#####
	//Horstmann Eq. A3-13
	//changed 11/23/04 G.B. in order to deal with leading edge singularity
	//old,	mu3t1 	= a2*t1+b2+sqrt(a2)*rt1;
	//old,	mu3t2 	= a2*t2+b2+sqrt(a2)*rt2;
	if(fabs(phij)<=DBL_EPS)
	{
		mu3t1 	= a2*t1+b2+sqrt(a2)*rt1;
		mu3t2 	= a2*t2+b2+sqrt(a2)*rt2;
	}
	else
	{
		mu3t1 	= .0001*xsij + a2*t1+b2+sqrt(a2)*rt1;
		mu3t2 	= .0001*xsij + a2*t2+b2+sqrt(a2)*rt2;
	}


//==================================================================
	//KHH eq. A3-8 through A3-16
	//indices off by 1! e.g. Horstmann G21 is G2[0] here
//==================================================================
	//modified G25
	//accounts for singularity along side edge of sheet
	//used for b and c-influence coef. of zeta velocity
	//added 4/13/04 G.B.
	//added factor of 0.01*etaj in order to blend in faster singularity
	//added 7/5/04 G.B.
	//G25b includes the factor b25 (b[4])
	tempS	= singfct	+ zetaASQR + t1*t1;
	tempSS	= singfct	+ zetaASQR + t2*t2;

	G25b 	= -0.5*log(tempSS/tempS);
	//G25c includes the factor c25 (c[4])		original - + - + (last one on line 1572)
	G25c 	= -etaj*log(tempS*tempSS);
//printf("etaj %lf G25c %lf  ",etaj,G25c);
//printf("tempS %lf tempSS %lf \n",tempS,tempSS);

		if (fabs(t1)>DBL_EPS)	//point of interest not on left edge
				G25c += t1*log(zetaASQR + t1*t1);
//printf("2: %lf  ",t1*log(zetaASQR + t1*t1));
		if (fabs(t2)>DBL_EPS)	//point of interest not on right edge
				G25c -= t2*log(zetaASQR + t2*t2);
//printf("3: %lf %lf \n",t2*log(zetaASQR + t2*t2),G25c);

	if (ABSD<=DBL_EPS)    //implicitly zeta!=0
	{//point A in plane that is spaned by bound vortex and zeta-axis

		G2[0]=0;		//G21
		G21b =0;		//modified G21*b21
		G21c =0;		//modified G21*c21
		G2[1]=0;		//G22

		//G25 Horstmann Eq. A3-14
		//used for b and c-influence coef. of eta velocity
		G2[4]	= 0.5*log((t2*t2+zetaASQR)/(t1*t1+zetaASQR));

		//G26 Horstmann Eq. A3-15, see also KHH program lines 1504 ff
		if(ABSzetaA<=DBL_EPS)
			G2[5]=0;  //if statement and =0 added 11/22/04 G.B.
		else
		{
			tempS=(zetaASQR+t1*t2);
			G2[5]	= atan((t2-t1)*zetaA/tempS);
				if(tempS<0 && (t2/zetaA)>0)		G2[5] += Pi;
				if(tempS<0 && (t2/zetaA)<0)		G2[5] -= Pi;
			G2[5] = G2[5]/zetaA;
		}
	}
	else
	{
		if(ABSzetaA<=DBL_EPS)
		{//point A lies i xsi-eta plane, but NOT on bound vortex line
		 //another possible correction is required further down when
	 	 //computing b2xi[1], c2xi[1],b2xi[2], c2xi[2]

	 	 //KHH line 1459 through 1492, removed 04/13/04 G.B.
		 	G2[0]=0;		//G21

			//modified G21, accounts for singularity along side edge of
			//sheet used for b and c-influence coef. of zeta velocity
			//include factors b21 and c21
			//added 04/13/04 G.B. and modified 02/22.05
			tempS = log(mu1t2/mu1t1);
			G21b = b[0]*beta1*(0.5*tempS  + G25b)/rho;	//G21*b21
			G21c = b[0]*beta1*(etaA*tempS + G25c)/rho;//G21*c21

		 	G2[1]=0;		//G22
			G2[5]=0;		//G26
		}
		else
		{//point A is NEITHER in xsi-eta plane, NOR on bound vortex

			//G25 Horstmann Eq. A3-14
			//used for b and c-influence coef. of eta velocity
			G2[4]	= 0.5*log((t2*t2+zetaASQR)/(t1*t1+zetaASQR));

			tempS	= 0.5*log(mu1t2/mu1t1);
			tempSS	= mu2t2-mu2t1;

			//G21 Horstmann Eq. A3-8
			//used for b and c-influence coef. of eta velocity
			G2[0]	= (beta1*(tempS-G2[4])+beta2*tempSS)/rho;
			G21b=0;G21c=0;

			//G22 Horstmann Eq. A3-9
			G2[1]	= (-beta2*(tempS-G2[4])+beta1*tempSS)/(rho*zetaA);
												//ln(1/a)=-ln(a)!!
//from KHH Flue
//107   G22 = (BETA2/(2.0*ZIA*RO))*ALOG(1.0/ARG1)+(BETA1/(ZIA*RO))*ARTA   14900000

			//G26 Horstmann Eq. A3-15, see also KHH program lines 1504 ff
			tempS=(zetaASQR+t1*t2);
			G2[5]	= atan((t2-t1)*zetaA/tempS);
				if(tempS<0 && (t2/zetaA)>0)		G2[5] += Pi;
				if(tempS<0 && (t2/zetaA)<0)		G2[5] -= Pi;
			G2[5] = G2[5]/zetaA;

		}//end if (ABSzetaA<DBL_EPS)
	}//end if (ABSD<DBL_EPS)

	//G24 Horstmann Eq. A3-12
	if (fabs(mu3t2)<=DBL_EPS || fabs(mu3t1)<=DBL_EPS)
	G2[3]=0;
	else
	G2[3]	= log(mu3t2/mu3t1)/sqrt(a2);

	//G23 Horstmann Eq. A3-11
	G2[2]	= (rt2-rt1-b2*G2[3])/a2;

	//G27 Horstmann Eq. A3-16
	G2[6]	= t2-t1;

//#printf("t1 %lf  t2 %lf zetaASQR  %lf  etaj %lf\n",t1,t2,zetaASQR,etaj);
//#printf("G21*b21 %lf  G21b  %lf  G21*c21 %lf  G21c  %lf\n",G2[0]*b[0],G21b,G2[0]*c[0],G21c);
//#printf("G25*b25 %lf  G25b  %lf  G25*c25 %lf  G25c  %lf\n",G2[4]*b[4],G25b,G2[4]*c[4],G25c);
//==================================================================
	//Horstmann Eq. A3-7
	//Horstmann b2xi here b2xi[0];  b2eta here b2xi[1]
//==================================================================

	if(ABSzetaA<=DBL_EPS)
	{	//point A does lie in xsi-eta plane,
		//wake can only induce zeta velocities

		//Horstmann Eq. A3-7, fist part b2eta, c2eta
		b2xi[1] = 0.0;
		c2xi[1] = 0.0;

		//Horstmann Eq. A3-7, second part bzeta, czeta
		//simplified, 4/14/04 G.B.
		b2xi[2] = G21b + G2[3]*b[3] + G25b;
		c2xi[2] = G21c + G2[2]*c[2] + G2[3]*c[3] + G25c + G2[6]*c[6];

//printf("b2z =  %lf  c2z = %lf\n",b2xi[2],c2xi[2]);
//printf("G21b= %lf G21c=  %lf G25b= %lf G25c= %lf\n",G21b,G21c,G25b,G25c);
//printf("G24 = %lf  c4  %lf  G27 = %lf  c7  %lf\n",G2[3],c[3],G2[6],c[6]);
//#fprintf(fp,"G21b= %lf G21c=  %lf G25b= %lf G25c= %lf ",G21b,G21c,G25b,G25c);

//#fprintf(fp," G21b %lf  G24 %lf  b4 %lf  G25b %lf  G21c %lf ",G21b,G2[3],b[3],G25b,G21c);
//#fprintf(fp," mu1t1 %2.15lf  mu1t2 %2.15lf  ",mu1t1,mu1t2);
//#fprintf(fp," xsiA %2.8lf  rt1 %2.8lf  rt2 %2.8lf  ",xsiA,rt1,rt2);
//#fprintf(fp," a2 %2.8lf  b2 %2.8lf  c2 %2.8lf  ",a2,b2,c2);
//#fprintf(fp," beta1 %2.8lf  beta2 %2.8lf  rho %2.8lf  b21 %2.8lf  ",beta1,beta2,rho,b[0]);
//#fprintf(fp," gamma1 %2.8lf  gamma2 %2.8lf  delta1 %2.8lf  delta1 %2.8lf  ",gamma1,gamma2,delta1,delta2);
//#fprintf(fp," t1 %lf  t2 %lf  ",t1,t2);

	}
	else//	if(ABSzetaA<DBL_EPS)
	{	//point does NOT lie in xsi-eta plane, zetaA NOT= 0
		//THIS IS THE NON-EXCEPTION

		//Horstmann Eq. A3-7, fist part b2eta, c2eta
		b2xi[1] = -zetaA*(G2[0]*b[3]+G2[1]*b[0]+G2[5]*b[4]);
		c2xi[1] = -zetaA*(G2[0]*c[3]+G2[1]*c[0]+G2[3]*c[2]\
				  		 +G2[4]*c[6]+G2[5]*c[4]);

		//Horstmann Eq. A3-7, second part bzeta, czeta
		//simplified, 4/14/04 G.B., modified 02/22/05 G.B.
		b2xi[2] = G2[0]*b[0] + G2[1]*b[1] + G2[3]*b[3] + G2[4]*b[4];
		c2xi[2] = G2[0]*c[0] + G2[1]*c[1] + G2[2]*c[2] + G2[3]*c[3]\
				+ G2[4]*c[4] + G2[5]*c[5] + G2[6]*c[6];

	}//end 	if(ABSzetaA<DBL_EPS)

    //wake only induces velocities in eta and zeta direction
    //CAUTION: that requires small angles or xsi and Uinf aligned
    b2xi[0] =	0;
    c2xi[0] =	0;

    //transform coefficients a, b, and c back into global co-system
    if(nuj*nuj > DBL_EPS)
    {	rotateX(b2xi,-nuj,b2x);
    	rotateX(c2xi,-nuj,c2x);
	}
	else
	{
		b2x[0] = b2xi[0]; b2x[1] = b2xi[1]; b2x[2] = b2xi[2];
		c2x[0] = c2xi[0]; c2x[1] = c2xi[1]; c2x[2] = c2xi[2];
	}

  }// else and if(ABSzetaA<=DBL_EPS && ABSD<=DBL_EPS)


//#fprintf(fp," b %lf  %lf  %lf \n",b2x[0],b2x[1],b2x[2]);

//#fclose(fp);
//#fprintf(fp," b21 %lf  b22 %lf  b23 %lf  b24 %lf  b25 %lf  b26 %lf  b27 %lf",b[0],b[1],b[2],b[3],b[4],b[5],b[6]);
//#fprintf(fp," c21 %lf  c22 %lf  c23 %lf  c24 %lf  c25 %lf  c26 %lf  c27 %lf",c[0],c[1],c[2],c[3],c[4],c[5],c[6]);
//#fprintf(fp,"G21b %lf G21 %lf  G22 %lf  G23 %lf  G24 %lf  G25 %lf  G26 %lf  G27 %lf",G21b,G2[0],G2[1],G2[2],G2[3],G2[4],G2[5],G2[6]);

//a2 coefficients of vortex sheet are zero, since vorticity is B+2*C*eta
a2x[0] =	0;
a2x[1] =	0;
a2x[2] =	0;
//##########################################
//
//FILE *fp;
//fp = fopen("output\\test.txt", "a");
//fprintf(fp,"b2y %lf  c2y %lf  ",b2xi[1],c2xi[1]);
//fprintf(fp,"b2z %lf  c2z %lf  ",b2xi[2],c2xi[2]);
//fprintf(fp,"xsiA: %lf  etaA: %lf  zetaA: %lf ",xsiA,etaA,zetaA);	//#
//fprintf(fp,"G21= %lf  G22= %lf  G23= %lf  G24= %lf  G25= %lf  G26= %lf  G27= %lf ",G2[0],G2[1],G2[2],G2[3],G2[4],G2[5],G2[6]);
//fprintf(fp,"G21*b %lf G21*c %lf G25*b %lf  G25*c %lf  ",G2[0]*b[0],G2[0]*c[0],G2[4]*b[4],G2[4]*c[4]);
//fprintf(fp,"G21b %lf G21c %lf G25b %lf  G25c %lf  ",G21b,G21c,G25b,G25c);
//fprintf(fp,"b1 %lf  b2 %lf  c1 %lf  c2 %lf  ",b[0],b[1],c[0],c[1]);	//#
//fprintf(fp,"mu1t1 %lf mu1t2 %lf ",mu1t1,mu1t2);	//#
//fprintf(fp,"mu2t1 %lf mu2t2 %lf ",mu2t1,mu2t2);	//#
//fprintf(fp,"beta1 %lf  beta2 %lf  rho %lf ",beta1,beta2,rho);	//#
//fprintf(fp,"G21b %lf  G21 %lf  G21c %lf  G21 %lf  ",G21b,G2[0],G21c,G2[0]);	//#
//fprintf(fp,"G25b %lf  G25c %lf G25 %lf t1 %lf  t2 %lf ",G25b,G25c,G2[4],t1,t2);
//fprintf(fp," \n");
//fclose(fp);
//
//#######################################
//fprintf(fp,"b2xsi= %lf  %lf  c2xsi= %lf   %lf\n",b2xi[1],b2xi[2],c2xi[1],c2xi[2]);
//printf("xsiA: %2.2lf %2.2lf %2.2lf ",xsiA,etaA,zetaA);	//#
//printf("b2x= %lf  %lf  c2x= %lf   %lf\n",b2x[1],b2x[2],c2x[1],c2x[2]);
//printf("a2x %2.2lf b2x %2.2lf c2x %2.2lf ",a2x[0],b2x[0],c2x[0]);
//printf("a2x %2.2lf b2x %2.2lf c2x %2.2lf ",a2x[1],b2x[1],c2x[1]);
//printf("a2x %2.2lf b2x %2.2lf c2x %2.2lf\n",a2x[2],b2x[2],c2x[2]);
}
/////&&&&&&&&&&&&&&
//===================================================================//
		//END function VortexSheetInduction
//===================================================================//

//===================================================================//
		//START function FixedWakeInduction
//===================================================================//
/*///&&&&&&&&&&&&&&
void VortexSheetInduction(const double xAi[3],const double xoj[3], \
						const double nuj, const double phij, const double etaj,\
						double a2x[3],double b2x[3],double c2x[3],const double singfct)
//void FixedWakeInduction(const double xAi[3],const double xoj[3], \
//						const double nuj, const double phij, const double etaj,\
//						double a2x[3],double b2x[3],double c2x[3])
{
//this function computes influence coefficients of ctr point i due to
//a fixed wake that stretches from the 1/4c line of elementary wing j
//into infinity.
//See also Appendix 3 of
//"Ein Mehrfach-Traglinienverfahren und seine Verwendung fuer Entwurf
//und  Nachrenchnung nichtplanarer Fluegelanordnungen",
//by K.-H. Horstmann, published 1987 by DFVLR, Germany, DFVLR-FB 87-51.
//
//the return values are in global coordinates
//
//Input:
//	xAi		element i control point in global coordinates
//	x0j		element j midspan point of 1/4c line in global coordinates
//	nuj		element j dihedral
//	phij	element j sweep
//	etaj	element j half span
//
//Output:
//	a2x[3]	induced coefficient of bound vortex in global system
//	b2x[3]	induced coefficient of bound vortex in global system
//	c2x[3]	induced coefficient of bound vortex in global system

double tempA[3],tempAA[3],tempS, tempSS;
double xsiA,etaA,zetaA,D;	//ctrl pt i coordinates in j element ref. frame
double ABSzetaA,ABSD;			//absolute values
double a2,b2,c2,t1,t2,rt1,rt2,b[7],c[7],G2[7];
double beta1,beta2,delta1,delta2,gamma1,gamma2,rho,epsilon;
double mu1t1,mu1t2,mu2t1,mu2t2,mu3t1,mu3t2;
double b2xi[3],c2xi[3];				//coefficients in local ref. frame
int i;								//counter

  //transformation of i-control point global coordinates into
  //local j-element coordinates. Horstmann Eq. 9
  scalar(xoj,-1,tempAA);
  vsum(tempAA,xAi,tempA);
  if(nuj*nuj > DBL_EPS)
  {
	  rotateX(tempA,nuj,tempAA);
	  xsiA	= tempAA[0];
	  etaA	= tempAA[1];
	  zetaA	= tempAA[2];
  }
  else
  {
	  xsiA	= tempA[0];
	  etaA	= tempA[1];
	  zetaA	= tempA[2];
  }

  ABSzetaA 	= fabs(zetaA);
  D 		= xsiA-etaA*tan(phij); //if D=0 then A in plane of bound vortex
  ABSD		= fabs(D);			  //and zeta axis

//##########################################
//FILE *fp;														//#
///##########################################
//#fp = fopen("check.txt", "a");									//#
//#fprintf(fp,"xAi: %lf  %lf  %lf \n",xAi[0],xAi[1],xAi[2]);		//#
//#fprintf(fp,"xoj: %lf  %lf  %lf \n",xoj[0],xoj[1],xoj[2]);		//#
//#fprintf(fp,"xsiA: %lf  etaA: %lf  zetaA: %lf \n",xsiA,etaA,zetaA);	//#
//#fprintf(fp,"phij: %lf  nuj: %lf  etaj: %lf \n",phij*RtD,nuj*RtD,etaj);	//#
//#fclose (fp);														//#
//##########################################
//#printf(" xsiA: %lf  etaA: %lf  zetaA: %lf eta: %lf\n",xsiA,etaA,zetaA,etaj);	//#

  if(ABSzetaA<DBL_EPS && ABSD<DBL_EPS)
  { //point A lies on bound vortex line/leading edge of vortexsheet
    //if phi is not zero than velocity is singular and value
    //is set to zero

	b2x[0] = 0;
	b2x[1] = 0;
	b2x[2] = 0;

	c2x[0] = 0;
	c2x[1] = 0;
	c2x[2] = 0;

	if(fabs(phij)<DBL_EPS)	//at leading edge of an unswept vortex sheet
	{						//page 22-3 of G.B. notes, 8/7/03
							//a problem arises when the point A is located
							//at the edge of an element, i.e. t1 or t2 = 0
							//no solution implemented yet, but a possible
							//solution is to compute the velocities a bit
							//to the right and to the left and average across
		b2xi[0] = 0;
		b2xi[1] = 0;

		c2xi[0] = 0;
		c2xi[1] = 0;

		t1 = etaA+etaj;
		t2 = etaA-etaj;

		tempS = log((t1*t1)/(t2*t2));

		b2xi[2] = tempS*0.5;
		c2xi[2] = (-4*etaj+etaA*tempS);

		//transform coefficients a, b, and c back into global co-system
		if(nuj*nuj > DBL_EPS)
		{	rotateX(b2xi,-nuj,b2x);
			rotateX(c2xi,-nuj,c2x);
		}
		else
		{
			b2x[0] = b2xi[0]; b2x[1] = b2xi[1]; b2x[2] = b2xi[2];
			c2x[0] = c2xi[0]; c2x[1] = c2xi[1]; c2x[2] = c2xi[2];
		}
	}  //end if(fabs(phij)<DBL_EPS)
  }
  else  //if(ABSzetaA<DBL_EPS && ABSD<DBL_EPS)
  {
	//Horstmann Eq. A3-17
	a2		= 1+tan(phij)*tan(phij);
	b2		= D*tan(phij);
	c2		= D*D+zetaA*zetaA;

	t1		= etaA+etaj;	//if t1 or t2 = 0 => see coment a few lines above!
	t2		= etaA-etaj;

	rt1		= sqrt(t1*t1*a2+2*t1*b2+c2);
	rt2		= sqrt(t2*t2*a2+2*t2*b2+c2);

	//Horstmann Eq. A3-18, indices off by 1!  Horstmann 1 is 0 here
	b[0]	= -D;
	b[1]	= zetaA*zetaA*tan(phij);
	b[2]	= 0;
	b[3]	= -tan(phij);
	b[4]	= -1;
	b[5]	= 0;
	b[6]	= 0;

	c[0]	= -2*(b[1]-etaA*b[0]);
	c[2]	= 2*tan(phij);
	c[3]	= 2*(xsiA-etaA*c[2]);
	c[4]	= -2*etaA;
	c[5]	= -2*zetaA*zetaA;
	c[6]	= 2;
	c[1]	= 0.5*c[5]*c[3];

	//Horstmann Eq. A3-10
	epsilon = b[0]*b[0]+b[1]*b[3];
	rho  	= sqrt(epsilon*epsilon+4*zetaA*zetaA*b2*b2);

	tempS 	= 0.5*(rho+epsilon);
	if (tempS<=DBL_EPS) beta1   = 0;
	else 				beta1	= -sqrt(tempS);
	tempS 	= 0.5*(rho-epsilon);
	if (tempS<=DBL_EPS) beta2   = 0;
	else 				beta2	= -sqrt(tempS);
	//see Horstmann program FLU lines 1443 through 1447
 	tempS=zetaA*b2;
	if (fabs(tempS) > DBL_EPS) beta2*=tempS/fabs(tempS);

	gamma1 	= (a2*beta2*zetaA+b2*beta1)/rho;
	gamma2 	= (a2*beta1*zetaA-b2*beta2)/rho;

	delta1 	= (b2*beta2*zetaA+c2*beta1)/rho;
	delta2 	= (b2*beta1*zetaA-c2*beta2)/rho;

	tempS	= gamma1*t1+delta1-rt1;
	tempSS	= gamma2*t1+delta2;
	mu1t1	= (tempS*tempS+tempSS*tempSS)/(t1*t1+zetaA*zetaA);
	//mu2t1 computed as in Horstmann Eq. A3-10. atan-values then corrected.
	if (t1*t1 < DBL_EPS)
	mu2t1 = Pi/2*ABSzetaA/zetaA+atan(tempSS/tempS);//|zeta/t1|->infinity
	else
	{
		mu2t1 = atan(zetaA/t1)+atan(tempSS/tempS);
	//Corrections according to Horstmamm program FLU lines 1473 and 1474
		if ((zetaA > 0.0) && (t1 < 0.0))	mu2t1+=Pi;
		if ((zetaA < 0.0) && (t1 < 0.0))	mu2t1-=Pi;
	}
	//Corrections according to Horstmamm program FLU lines 1484 and 1485
	if (tempS < 0.0)						mu2t1+=Pi;
	if ((tempSS < 0.0) && (tempS > 0.0))	mu2t1+=2*Pi;

	tempS	= gamma1*t2+delta1-rt2;
	tempSS	= gamma2*t2+delta2;
	mu1t2	= (tempS*tempS+tempSS*tempSS)/(t2*t2+zetaA*zetaA);
	//mu2t2 computed as in Horstmann Eq. A3-10. atan-values then corrected.
	if (t2*t2 < DBL_EPS)
	mu2t2 = Pi/2*ABSzetaA/zetaA+atan(tempSS/tempS);//|zeta/t2|->infinity
	else
	{
		mu2t2 = atan(zetaA/t2)+atan(tempSS/tempS);
	//Corrections according to Horstmamm program FLU lines 1466 1467
		if ((zetaA > 0.0) && (t2 < 0.0))	mu2t2+=Pi;
		if ((zetaA < 0.0) && (t2 < 0.0))	mu2t2-=Pi;
	}
	//Corrections according to Horstmamm program FLU lines 1479 and 1480
	if (tempS < 0.0)						mu2t2+=Pi;
	if ((tempSS < 0.0) && (tempS > 0.0))	mu2t2+=2*Pi;

	//Horstmann Eq. A3-13
	mu3t1 	= a2*t1+b2+sqrt(a2)*rt1;
	mu3t2 	= a2*t2+b2+sqrt(a2)*rt2;

//==================================================================
	//KHH eq. A3-8 through A3-16
	//indices off by 1! e.g. Horstmann G21 is G2[0] here
//==================================================================

	//G21 & G26
	if (ABSD<DBL_EPS)
	{//point A in plane that is spaned by bound vortex and zeta-axis
	 //another possible correction is required further down when
	 //computing b2xi[1], c2xi[1],b2xi[2], c2xi[2]

		G2[0]=0;		//G21
		G2[1]=0;		//G22
	}
	else
	{
		if (ABSzetaA<DBL_EPS)
		{//point A lies i xsi-eta plane, but NOT on bound vortex line
		 //another possible correction is required further down when
	 	 //computing b2xi[1], c2xi[1],b2xi[2], c2xi[2]
	 	 //KHH line 1459 through 1492
			G2[0]=(beta1*0.5*log(mu1t2/mu1t1))/rho;	//G21
			G2[1]=0;								//G22
//fp = fopen("test.txt", "a");
//fprintf(fp,"line 1925  G21 %lf  G22 %lf \n",G2[0],G2[1]);
//fclose(fp);

		}
		else
		{//point A is NEITHER in xsi-eta plane, NOR on bound vortex
			tempS	= 0.5*log(mu1t2/mu1t1);
			tempSS	= mu2t2-mu2t1;

			//G21 Horstmann Eq. A3-8
			G2[0]	= (beta1*tempS+beta2*tempSS)/rho;

			//G22 Horstmann Eq. A3-9
			G2[1]	= (-beta2*tempS+beta1*tempSS)/(rho*zetaA);
												//ln(1/a)=-ln(a)!!
//fp = fopen("test.txt", "a");
//fprintf(fp,"line 1941  G21 %lf  G22 %lf \n",G2[0],G2[1]);
//fclose(fp);

//from KHH Flue
//107   G22 = (BETA2/(2.0*ZIA*RO))*ALOG(1.0/ARG1)+(BETA1/(ZIA*RO))*ARTA   14900000

		}//end if (ABSzetaA<DBL_EPS)
	}//end if (ABSD<DBL_EPS)

	//G24 Horstmann Eq. A3-12
	if (fabs(mu3t2)<=DBL_EPS || fabs(mu3t1)<=DBL_EPS)
	G2[3]=0;
	else
	G2[3]	= log(mu3t2/mu3t1)/sqrt(a2);

	//G23 Horstmann Eq. A3-11
	G2[2]	= (rt2-rt1-b2*G2[3])/a2;

	//G25 Horstmann Eq. A3-14
	G2[4]	= 0.5*log((t2*t2+zetaA*zetaA)/(t1*t1+zetaA*zetaA));

	if(ABSzetaA<=DBL_EPS)  		//G26
	{   //point A does lie in xsi-eta plane
		G2[5]=0.0;		// G26
	}
	else
	{	//point A does not lie in xsi-eta plane
		//G26 Horstmann Eq. A3-15, see also KHH program lines 1504 ff
		tempS=(zetaA*zetaA+t1*t2);
		G2[5]	= atan((t2-t1)*zetaA/tempS);
		if(tempS<0 && (t2/zetaA)>0)
			G2[5] += Pi;
		if(tempS<0 && (t2/zetaA)<0)
			G2[5] -= Pi;
		G2[5] = G2[5]/zetaA;
	}//done with G26

	//G27 Horstmann Eq. A3-16
	G2[6]	= t2-t1;

//==================================================================
	//Horstmann Eq. A3-7
	//Horstmann b2xi here b2xi[0];  b2eta here b2xi[1]
//==================================================================

	// Horstmann Eq. A3-7, second part for b2zeta and c3zeta
	b2xi[2] =	0;
	c2xi[2] =	0;
	for(i=0;i<7;i++)
		{
			b2xi[2] += G2[i]*b[i];
			c2xi[2] += G2[i]*c[i];
		}

	if(ABSzetaA<=DBL_EPS)
	{	//point A does lie in xsi-eta plane,
		//wake can only induce zeta velocities

		//Horstmann Eq. A3-7
		b2xi[1] = 0;
		c2xi[1] = 0;

		if (ABSD>=DBL_EPS)
		{	//point A falls into xsi-eta plane, but not on bound vortex line

			tempS	= 0.5*log(mu1t2/mu1t1);
			tempSS	= mu2t2-mu2t1;

 			//G22 KHH line 1536
			G2[1]	= (-beta2*tempS+beta1*tempSS)/(rho		);
												//ln(1/a)=-ln(a)!!
			b2xi[1] = G2[1]*b[0];	//<= SHOULD THIS BE b2xi[2]????????????
			c2xi[1] = G2[1]*c[0];
//fp = fopen("test.txt", "a");
//fprintf(fp,"line 2015  G21 %lf  G22 %lf \n",G2[0],G2[1]);
//fclose(fp);

		}//end if (ABSD>DBL_EPS)
	}
	else//	if(ABSzetaA<DBL_EPS)
	{	//point does NOT lie in xsi-eta plane, zetaA NOT= 0
		//THIS IS THE NON-EXCEPTION

		//Horstmann Eq. A3-7, fist part b2eta, c2eta
		b2xi[1] = -zetaA*(G2[0]*b[3]+G2[1]*b[0]+G2[5]*b[4]);
		c2xi[1] = -zetaA*(G2[0]*c[3]+G2[1]*c[0]+G2[3]*c[2]\
				  +G2[4]*c[6]+G2[5]*c[4]);
	}//end 	if(ABSzetaA<DBL_EPS)

    //wake only induces velocities in eta and zeta direction
    //CAUTION: that requires small angles or xsi and Uinf aligned
    b2xi[0] =	0;
    c2xi[0] =	0;

    //transform coefficients a, b, and c back into global co-system
    if(nuj*nuj > DBL_EPS)
    {	rotateX(b2xi,-nuj,b2x);
    	rotateX(c2xi,-nuj,c2x);
	}
	else
	{
		b2x[0] = b2xi[0]; b2x[1] = b2xi[1]; b2x[2] = b2xi[2];
		c2x[0] = c2xi[0]; c2x[1] = c2xi[1]; c2x[2] = c2xi[2];
	}
  }// else and if(ABSzetaA<DBL_EPS && ABSD<DBL_EPS)

//a2 coefficients of vortex sheet are zero, since vorticity is B+2*C*eta
  a2x[0] =	0;
  a2x[1] =	0;
  a2x[2] =	0;

//##########################################
//FILE *fp;
//fp = fopen("test.txt", "a");
//fprintf(fp,"xsiA: %lf  etaA: %lf  zetaA: %lf ",xsiA,etaA,zetaA);	//#
//fprintf(fp,"b2y= %lf  c2y= %lf ",b2xi[1],c2xi[1]);
//fprintf(fp,"b2z= %lf  c2z= %lf ",b2xi[2],c2xi[2]);
//fprintf(fp,"G21= %lf  G22= %lf  G23= %lf  G24= %lf  G25= %lf  G26= %lf  G27= %lf\n",G2[0],G2[1],G2[2],G2[3],G2[4],G2[5],G2[6]);
//fclose(fp);
//##########################################
//
//FILE *fp;
//fp = fopen("test.txt", "a");
//fprintf(fp,"b2y %lf  c2y %lf  ",b2xi[1],c2xi[1]);
//fprintf(fp,"b2z %lf  c2z %lf  ",b2xi[2],c2xi[2]);
//fprintf(fp,"G21= %lf  G22= %lf  G23= %lf  G24= %lf  G25= %lf  G26= %lf  G27= %lf ",G2[0],G2[1],G2[2],G2[3],G2[4],G2[5],G2[6]);
//fprintf(fp,"G21*b %lf G21*c %lf G22*b %lf  G22*c %lf  ",G2[0]*b[0],G2[0]*c[0],G2[1]*b[1],G2[1]*c[1]);
//fprintf(fp,"G23*b %lf G23*c %lf G24*b %lf  G24*c %lf  ",G2[2]*b[2],G2[2]*c[2],G2[3]*b[3],G2[3]*c[3]);
//fprintf(fp,"G25*b %lf G25*c %lf G26*b %lf  G26*c %lf  ",G2[4]*b[4],G2[4]*c[4],G2[5]*b[5],G2[5]*c[5]);
//fprintf(fp,"G27*b %lf G27*c %lf ",G2[6]*b[6],G2[6]*c[6]);
//fprintf(fp," \n");
//fclose(fp);
//#######################################//
//printf("xsiA: %2.2lf %2.2lf %2.2lf ",xsiA,etaA,zetaA);	//#
//printf("b= %lf  %lf  c= %lf   %lf\n",b2xi[1],b2xi[2],c2xi[1],c2xi[2]);
//#printf("a2x = %lf  b2x = %lf  c2x = %lf\n",a2x[0],b2x[0],c2x[0]);
//#printf("a2x = %lf  b2x = %lf  c2x = %lf\n",a2x[1],b2x[1],c2x[1]);
//#printf("a2x = %lf  b2x = %lf  c2x = %lf\n",a2x[2],b2x[2],c2x[2]);
}
//===================================================================//
		//END function FixedWakeInduction
//===================================================================//
//*///&&&&&&&&&&&&&&
