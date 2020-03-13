// sign convention according to Horstmann' thesis

// definition of type general //
struct GENERAL
{
	char inputfilename[126]; 
	bool flagVISCOUS; // Turn on/off viscous corrections (O = OFF, 1 = ON)
	bool flagCAMBER; // Turn on/off camber (O = OFF, 1 = ON)
	double Uinf;		//free stream velocity
	double alpha, beta;	//angle of attack, sideslip angle
	double S, b,AR;		//reference area, span, aspect ratio
	double U[3];		//free stream velocity vector
	double density;		//density of fluid
	double nu;			//kinematic viscosity
	double RefPt[3];	//global reference point for moments
						//added 1-0-13-2006 G.B. also used for CG
	bool flagHORZ;		// On/off flag for flight in horizontal plane GB 3.9.20
						//if on, geometry is rotated by alpha and Uinf in xy-plane
	
	// Circling flight info D.F.B 02-2020
	bool flagCIRC;		// On/off flag for circling flight
	double bank;		// Bank angle (rad)
	double Ws; 			// Upwind velocity
	double gradient;	// Velocity gradient
	
	double W;			//aircraft weight; added G.B. 8-8-07
	double cmac;		//mean aeodynamic chord; added G.B. 8-8-07
	double CMoWing;		//zero-lift moment coefficient of wing

	int maxtime;		//maximal number of time steps
	double deltime;		//time step width
	double deltae;		//square of deltae that determins convergence criteria
						//of do-while loop in main program.
						//set 0 if only time stepping is desired

	int sym;			//symmetrical geometry flag (=1 sym. 0= assym.)
	int linear;			//linear theory flag (=1 applied, 0= not applied)
	int steady;			//steady (=1)/unsteady (default) aerodynamics flag
	int relax;			//relaxed wake in time stepping scheme (=1 applied, 0=not)
	int trim;			//longitudinal trim flag (=1 for trimming, =0 no trim; =1 m=1!)
	int nowing;			//number of separated wings
	int wing1[5],wing2[5];//span index of edges 1 (left) and 2 (right) of wing
	int panel1[5],panel2[5];//indices of panels at left and right edge of wing
    int dve1[5],dve2[5];//indices of first and last dve of wing
						//G.B. 11-5-06
	int noairfoils;		//no. of airfoils used G.B. 8-9-07
	
	int m; 				//number lifting lines/elementary wings in chord direction
	int nopanel;		//number of panels
	int noelement;		//number of elementary wings
	int nospanelement;	//number of elements in span direction, noelement/m
	int Dsize;			//3*info.noelement, dimension of R and D

	int noVT,noFus;		//number of V-tail and fuselage sections
};

// definition of type panel//
struct PANEL
{
	double x1[3];		//panel side 1 leading edge coordinates
    double c1, eps1;	//panel side 1 chord and incident angle
    double u1[3];		//free stream velocity variation at panel side 1
	int BC1;			//boundary condition at panel side 1
    int airfoil1;       //airfoil at panel side 1  added GB 2-14-20
    double hinge1;		//hinge location in %c on side 1 of panel

    double x2[3];		//panel side 2 leading edge coordinates
    double c2, eps2;	//panel side 2 chord and incident angle
    double u2[3];		//free stream velocity variation at panel side 2
	int BC2;			//boundary condition at panel side 2
    int airfoil2;       //airfoil at panel side 2   added GB 2-14-20
    double hinge2;		//hinge location in %c on side 2 of panel

	int left, right;	//left and right panel neighbors. 0 -> free end
    int n;				//number of spanwise elementary wings
    int m;              //number of chordwise elementary wings GB 2-9-20
	int airfoil;		//airfoil file number G.B. 8-9-07

    int TE1,TE2;		//indices of left and right DVE @ TE of the panel
						//added 8/12/05 G.B.
    int LE1,LE2;		//index of first and right DVE @ LE of panel 
    					//GB 2-20-20

    //	Panel Boundary Conditions:
	//	First Digit: 	0 - undefined circulation strength
	//					1 - zero circulation (free end)
	//					2 - circul. strength equal to neighboring el. wing
	//	Second Digit: 	0 - undefined circulation slope
	//					1 - zero slope in circulation
	//					2 - circul. slope equal to neighboring el. wing
	//	Third Digit: 	0 - undefined circulation curvature
	//					1 - zero circulation curvature
	//					2 - circul. curvature equal to neighboring el. wing
	//
	// BC = 100 - gamma=0
	// BC = 110 - gamma=0 AND gamma'=0 (free tip)
	// BC = 220 - gamma[i] = gamma[j]  AND  gamma[i]' = gamma[j]'
	// BC = 022 - gamma[i]' = gamma[j]'  AND  gamma[i]" = gamma[j]"
	// BC = 010 - gamma' = 0 (center-line case with symmetrical conditions)
};

// definition of type bound vortex //
struct BOUND_VORTEX
{
	double xo[3];		//midspan location of bound vortex in global ref.frame
	double xA[3];		//control point location in global ref. frame
	double normal[3];	//surface normal at control point in global ref. frame

	double eta;			//half span of elementary bound vortex
	double chord;		//mid span chord of elementary wing
	double phi, nu;		//sweep, dihedral of elementary bound vortex

	double u[3];		//local free stream vel. variation in global ref. frame
						//at mid span of elementary bound vortex
	double uind[3];		//induced velocity at this point  //added 7-20-05 G.B.

	double A,B,C;		//coefficients of vort. distribution across elementary
						//bound vortex: Gamma =A+B*y+C*y^2 with y=-eta .. eta

						//elementary bound vortex:
	double N_free[2];	//	- lift and side forces/density due to free stream
	double N_ind[2];	//	- lift and side forces/density due to induced flow

	double CDind;		//elementary bound vortex induced drag coefficient
};

// definition of type Distributed-Vorticity element //
struct DVE
{
	double A,B,C;		//coefficients of circulation distribution
						//across DVE
						//l.e.vortex:   Gamma = A+B*y+C*y^2
						//vortex sheet: gamma =   B  +2*C*y
						//t.e.vortex:   Gamma =-A-B*y-C*y^2
						//with y=-eta .. eta

	//double K;			//total circulation of DVE K = A +1/3*eta^2*C
    
    double A_old, B_eta, Csqeta;    //replaces K-approach to compute new
                        //spanwise circulation distribution
                        //average circulation is maintained, that is
                        //A remains constant and B and C scale with spanwise
                        //stretching of wakeDVEs  GB 2-6-20

	double xo[3];		//DVE referece point in global ref.frame
						//It is located midspan and midchord

	double eta;			//half span of elementary bound vortex
	double xsi;			//half chord of mid span of DVE
	double S;			//area of element G.B. 8-10-07
	double nu,epsilon;	//incident (pitch) and roll angle of DVE
	double psi;			//yaw angle of DVE
	double phiLE,phiTE;	//leading and trailing edge sweep
	double phi0;		//sweep of mid-chordline

	double normal[3];	//DVE surface normal
	double u[3];		//ind. velocity in ref. pt. in global coordinate
	//NOTE! surface DVEs u stores the local free stream vector!!
	double U[3];		//normalized local free stream vel. in ref. pt.
	double xleft[3];	//coordinates of mid-chord of left edge
	double singfct;		//rate at which additional singularity at the edge
						//of the vortex sheet decays
						//## singfct added 2/8/05 G.B.
	double Force[3],Moment[3];//aerodynamic forces and moments of DVE in global
						//reference frame.
						//Moments with respect to global reference point
						//added 10-13-2006 (it's a Friday!) G.B.
	double x1[3],x2[3];	//points at left and right end of leading edge bound vortex
						//needed for displacing surface DVE's
	double u1[3],u2[3];	//freestream velocities in x1 and x2  added 10-29-06 G.B.
	double uTE[3][3];	//freestream at trailing edge. xTE[0] ctr pt, xTE[1] -80%, xTE[2] +80% added D.F.B. 03-2020
	double xTE[3],TEvc[3]; //center point at and vector along trailing edge of previous
						//timestep
	int airfoil[2];		//airfoil file number G.B. 8-9-07; two airfoiils of panel edges 1 and 2 GB 2-14-20
	double ratio;		//interpolation ratio from left (=0) to right side (=1) GB 2-14-20
};

// definition of type wing//
struct Wing
{
	int dve1,dve2;		//first and last DVE index of wing
	int panel1,panel2;	//first and last panel index of wing
	double CL,CDi,e;	//lift, induced drag, span efficiency of wing
};
