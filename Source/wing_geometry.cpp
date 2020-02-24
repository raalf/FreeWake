//generates trailing edge vortex elements of a panel
void Trailing_Edge_Generation\
					(const PANEL, const GENERAL, const BOUND_VORTEX*,\
					 BOUND_VORTEX *, int &);
//identifies edges of separate wings
void Wing_Generation(const PANEL*,int,int[5],int[5],int[5],int[5],\
                     int[5],int[5]);
//generates surface DVE elements
void Surface_DVE_Generation(const GENERAL,const PANEL *,DVE *,double ***);
//moves wing by delta x every time step,
void Move_Wing(const GENERAL, DVE*);
//moves flexible wing by delta x every time step,
void Move_Flex_Wing(const GENERAL, DVE*);
//moves DVEs according do the input camber line
void Apply_Camber(const PANEL*, double[3], double[3] ,double ***, int, int, double,double *,double *,double *,double *);
//===================================================================//
		//FUNCTION Trailing_Edge_Generation
		//generates trailing edge vortex elements of panel 'j'
//===================================================================//

void Trailing_Edge_Generation(const PANEL panel,const GENERAL info,\
							  const BOUND_VORTEX* elementPtr,\
							  BOUND_VORTEX* trailedgePtr, int &l)
{
	//defines trailing edge of panel 'j' and devides it into trailing edge
	//elements for cumulating bound vortex upstream of it
	//
	//input:
	// 	For panel 'j'(in panel):
	//	x1[], x2[]	-x,y,z coordinates of leading edge corners
	//	c1,c2		-chord length of panel sides
	//	eps1, eps2	-incident angle of panel sides
	//	u1[], u2[]	-local free stream velocity variation at panel sides
	//				 (for rotating wings}
	//	n 			-number of elementary wings in span direction
	//
	//  info.m		-number of elementary wings in chord direction
	//  info.U		-free stream velocity
	//
	//	l			-trailing edge element index, incremented by one at end
	//
	//ouput:
	//definition of trailing edge element properties.
	// 	For each trailing edge element in trailedgePtr:
	//	xo[]		-midspan location of trailing edge vortex in global ref.frame
	//	xA[]		-not used
	//	normal[]	-not used
	//  chord 		-not used
	//	eta			-half span of elementary wing
	//	phi 		-sweep of elementary wing
	//	nu			-dihedral of elementary wing
	//	u[3]		-local free stream velocity at mid span
	//				 of trailing edge element in global ref.frame.

int k;								//loop counters, k=0..(panel.n-1)
double phi, nu, eta;  				//panel sweep, dihedral, span
double xte1[3],xte2[3],xte[3];		//trailing edge points of edge 1 and 2
double tempS, tempA[3],tempA1[3];	//temporary scalar, array

	scalar(panel.x1,-1,tempA);	//negates x1
	vsum(panel.x2,tempA,tempA1); //x_1/4. = x1/4_2-x1/4_1
	tempS= sqrt(tempA1[1]*tempA1[1]+tempA1[2]*tempA1[2]); //panel span at 1/4c
	nu = asin(tempA1[2]/tempS);				//panel dihedral at 1/4c

	//computing the left trailing edge point of panel
	xte1[0] = panel.x1[0]+0.75*panel.c1*cos(panel.eps1);
	xte1[1] = panel.x1[1]+0.75*panel.c1*sin(panel.eps1)*sin(nu);
	xte1[2] = panel.x1[2]-0.75*panel.c1*sin(panel.eps1)*cos(nu);

	//computing the right leading edge point of panel
	xte2[0] = panel.x2[0]+0.75*panel.c2*cos(panel.eps2);
	xte2[1] = panel.x2[1]+0.75*panel.c2*sin(panel.eps2)*sin(nu);
	xte2[2] = panel.x2[2]-0.75*panel.c2*sin(panel.eps2)*cos(nu);

	//computing trailing edge vector
	xte[0] =  xte2[0]-xte1[0];
	xte[1] =  xte2[1]-xte1[1];
	xte[2] =  xte2[2]-xte1[2];

//printf("xte\t%2.5lf %2.5lf %2.5lf\n",xte[0]/norm2(xte),xte[1]/norm2(xte),xte[2]/norm2(xte));

	phi = asin(xte[0]/norm2(xte));		//sweep of panel's trailing edge
	tempS= sqrt(xte[1]*xte[1]+xte[2]*xte[2]);	//panel span at TE
	nu = asin(xte[2]/tempS);					//trailing edge dihedral
	eta = tempS/panel.n*0.5;  //halfspan of trailing edge of elementary wings

	//loop over number of spanwise elements 'n'
	for (k=0;k<panel.n;k++)
	{
		tempS= (.5+k)/panel.n;

		//computes mid span point of trailing edge element
		scalar(xte,tempS,tempA);
		vsum(xte1, tempA, trailedgePtr[l].xo);

		trailedgePtr[l].phi = phi;	// elementary wing sweep
		trailedgePtr[l].nu = nu;	// elementary wing dihedral
		trailedgePtr[l].eta = eta;	// elementary wing halfspan

		//computes local free stream velocity at xo location
		//varies along span in case of rotating wings
		scalar(panel.u2,tempS,tempA);		 //u2/n*(.5+k)
		scalar(panel.u1,(1-tempS),tempA1);	 //-u1(1/n*(.5+k)-1)
		vsum(tempA,tempA1,trailedgePtr[l].u);//u2+(u2-u1)/n*(.5+k)

		//KHH small angle approach,
		//ctrl. point lies in shed vortex sheet 3/4chord aft of bound vortex
		//added 6/24/03 G.B.
		if ((info.linear == 1)  && (panel.m == 1))
		{
			trailedgePtr[l].xo[0]=elementPtr[l].xo[0]+elementPtr[l].chord*.75;
			trailedgePtr[l].xo[1]=elementPtr[l].xo[1];
			trailedgePtr[l].xo[2]=elementPtr[l].xo[2];
		}

		l++; //increment elementary wing index l=0..(noelement-1)
	}	//End loop over k

}
//===================================================================//
		//END FUNCTION Trailing_Edge_Generation
//===================================================================//

//===================================================================//
		//FUNCTION Wing_Generation
		//identifies separate wings
//===================================================================//
void Wing_Generation(const PANEL* panelPtr,const int nopanel,\
					 int wing1[5],int wing2[5],int panel1[5],int panel2[5],\
                     int dve1[5], int dve2[5])
{
	//identifies the separate wings and their span indices of their tips
	//
	//input:
	// 	panel		- information on panels
	//  nopanel		- number of panels
	//
	//ouput:
	//
	//	wing1[wing]	- span index of left edge of wing "wing"
	//	wing2[wing]	- span index of right edge of wing "wing"
	//	panel1[wing]- index of left panel of wing "wing"
	//	panel2[wing]- index of right panel of wing "wing"
	//				panel1 and panel2 added G.B. 11-5-06
    //    dve1[wing]- index of first dve of wing "wing"
    //    dve2[wing]- index of last dve of wing "wing" added GB 2-9-20

int k,span=0,wing=0,index=0;				//loop counters, k=0..(panel.n-1)

	//identifies left edges (edge 1) of wings
	for(k=0;k<nopanel;k++)
	{
		//looks at panels left edge (edge 1)
		if(panelPtr[k].left==0)
		{
			wing1[wing]=span;
			panel1[wing]=k;
            dve1[wing]=index;
			wing++;			//advance to next wing
		}
		span += panelPtr[k].n; //move to next panel
        index += panelPtr[k].n*panelPtr[k].m; //index to first dve of next panel
	}

//identifies right edges (edge 2) of wings

	//initialize
	wing = 0;
	span = -1;
    index = -1;

	for(k=0;k<nopanel;k++)
	{
		span += panelPtr[k].n; //move to next panel
        index += panelPtr[k].n*panelPtr[k].m; //index to first dve of next
		//looks at panels left edge (edge 1)
		if(panelPtr[k].right==0)
		{
			wing2[wing]=span;
			panel2[wing]=k;
            dve2[wing]=index;
			wing++;			//advance to next wing
		}
	}
//    for(k=0;k<info.nowing;k++)
//    {
//       printf("  wing[%d] dve1 %d  dve2 %d panel1 %d panel2 %d\n",k,info.dve1[k],info.dve2[k],info.panel1[k],info.panel2[k]);
//    } //###

}
//===================================================================//
		//END FUNCTION Wing_Generation
//===================================================================//

//===================================================================//
		//FUNCTION Surface_DVE_Generation
//===================================================================//
void Surface_DVE_Generation(const GENERAL info,const PANEL* panelPtr,\
							DVE* surfacePtr, double ***camberPtr)
{
//generates surface Distributied-Vorticity elements. The element
//exists of a leading and
//trailing edge vortex with parabolic circulation distributions that
//are gamma=A+B*eta+C*eta^2 and gamma=-A-B*eta-C*eta^2 respectively.
//Inbetween the two vortices is a vortex sheet with a linear vorticity
//distribution.  The elements location is defined with its control point
//that is located at half span and and half chord.
//The element is planar and and is rotated about the x-axis by the
//dihedral or "bank angle" nu and pitched about the new y axis by epsilon.
//in the current version (6/25/03) these angles are determined with the
//surface normal of the elementary wings.
//
//input:
//	l			-elementary wing index, incremented by one at end
//				 of loop 'k' over number of spanwise elementary wings.
//				 Hence, l = 0 .. (noelement-1)
//
//ouput:
//definition of properties of surface DVE.
// 	For each surface DVE in surfacePtr:
//	xo[]		-DVE reference and control point, global ref.frame
//  xsi			-half chord at midspan of DVE
//	eta			-half span of DVE
//	phiLE 		-sweep of leading edge of DVE
//	phiTE		-sweep of trailing edge of DVE
//	nu			-dihedral of DVE
//	epsilon		-incident angle of DVE
//	u[3]		-velocity at xA in global ref.frame.
//  singfct		-rate at which singularity at edge of vortex sheet decays
//				 ## singfct added 2/8/05 G.B.
//	ratio 		-interpolation ratio between left and right edge a panel 14/2/20 G.B.

int i,m,n,wing;			//loop counters
int l=0;			//l => suface DVE counter,
double singfct;		//decay rate of singularity at edge of vortex sheet
double tempS,tempA[3],tempAA[3];

double xquart[3];					//1/4chord line of panel, defined in input
double x1LE[3],x2LE[3],x1[3],x2[3]; //edge points of panel and spanwise rows
double xLE[3],xsiLE[3];	//vector along LE of a spanewise row of surface DVEs
double delchord1,delchord2,delchord; //edge 1&2 chord,chord span-increment,
double deleps,ceps,seps;//increments of epsi; cos/sin of incidence angle
double nu,nu2;		//nus of panel 1/4c line, LE of DVE row
double delTANphi;	//tan(phiLE-phiTE)
double delX[3];		//vector from center of leading edge to control point
double eps1,eps2;
double chord1,chord2;
bool flagCAMBER = 0;

	//loop over number of panels
	for (i=0;i<info.nopanel;i++)
	{
		//xquart: vector between 1/4 chord points of panel edges
		scalar(panelPtr[i].x1,-1,tempA);	//negates x1
		vsum(panelPtr[i].x2,tempA,xquart); //x_l.e. = x1/4_2-x1/4_1

		//panel span
		tempS = sqrt(xquart[1]*xquart[1]+xquart[2]*xquart[2]); 

		//1/4chord line dihedral
		nu = asin(xquart[2]/tempS);

		delchord1 = panelPtr[i].c1/panelPtr[i].m;	//chordwise increment left side
		delchord2 = panelPtr[i].c2/panelPtr[i].m;	//chordwise increment right side
		delchord  = (delchord2-delchord1)/panelPtr[i].n;//chord/span increment

		//tangent of change in sweep angle of each chordwise row
		delTANphi = (delchord2-delchord1)/tempS;

		//spanwise incident increments
		deleps = (panelPtr[i].eps2-panelPtr[i].eps1)/panelPtr[i].n;

		/* Changed 20.02.20 D.F.B. User now defined panel at LE and the wing
		/ twists about the LE of the wing (Instead of the quater chord)
		//computing the left leading edge point of panel
		x1LE[0] = panelPtr[i].x1[0]-0.25*panelPtr[i].c1*cos(panelPtr[i].eps1);
		x1LE[1] = panelPtr[i].x1[1]-0.25*panelPtr[i].c1*sin(panelPtr[i].eps1)\
																	*sin(nu);
		x1LE[2] = panelPtr[i].x1[2]+0.25*panelPtr[i].c1*sin(panelPtr[i].eps1)\
																	*cos(nu);															
		computing the right leading edge point of panel
		x2LE[0] = panelPtr[i].x2[0]-0.25*panelPtr[i].c2*cos(panelPtr[i].eps2);
		x2LE[1] = panelPtr[i].x2[1]-0.25*panelPtr[i].c2*sin(panelPtr[i].eps2)\
																	*sin(nu);
		x2LE[2] = panelPtr[i].x2[2]+0.25*panelPtr[i].c2*sin(panelPtr[i].eps2)\
																	*cos(nu);
		*/

		//computing the left leading edge point of panel
		x1LE[0] = panelPtr[i].x1[0];
		x1LE[1] = panelPtr[i].x1[1];
		x1LE[2] = panelPtr[i].x1[2];	
		//computing the right leading edge point of panel
		x2LE[0] = panelPtr[i].x2[0];
		x2LE[1] = panelPtr[i].x2[1];							
		x2LE[2] = panelPtr[i].x2[2];															
		//loop over number of chordwise elements 'info.m'
		//removed GB 2-9-20		for (m=0;m<info.m;m++)
        for (m=0;m<panelPtr[i].m;m++)
		{
			/* Changed 20.02.20 D.F.B. There was issues with multiple m so now
			// panels are defined in input along LE instead of quater chord
			tempS = (0.25+m); //leading edge location/chord, 1/4c of spanwise row

			//computing left LE locations of current spanwise row of DVEs
			x1[0] = x1LE[0]+delchord1*tempS*cos(panelPtr[i].eps1);
			x1[1] = x1LE[1]+delchord1*tempS*sin(panelPtr[i].eps1)*sin(nu);
			x1[2] = x1LE[2]-delchord1*tempS*sin(panelPtr[i].eps1)*cos(nu);
			

			//computing right LE locations of current spanwise row of DVEs
			x2[0] = x2LE[0]+delchord2*tempS*cos(panelPtr[i].eps2);
			x2[1] = x2LE[1]+delchord2*tempS*sin(panelPtr[i].eps2)*sin(nu);
			x2[2] = x2LE[2]-delchord2*tempS*sin(panelPtr[i].eps2)*cos(nu);
			*/
			//computing left LE locations of current spanwise row of DVEs
			x1[0] = x1LE[0]+delchord1*m*cos(panelPtr[i].eps1);
			x1[1] = x1LE[1]+delchord1*m*sin(panelPtr[i].eps1)*sin(nu);
			x1[2] = x1LE[2]-delchord1*m*sin(panelPtr[i].eps1)*cos(nu);

			//computing right LE locations of current spanwise row of DVEs
			x2[0] = x2LE[0]+delchord2*m*cos(panelPtr[i].eps2);
			x2[1] = x2LE[1]+delchord2*m*sin(panelPtr[i].eps2)*sin(nu);
			x2[2] = x2LE[2]-delchord2*m*sin(panelPtr[i].eps2)*cos(nu);

			if(info.flagCAMBER){
			Apply_Camber(panelPtr,x1,x2,camberPtr, m, i, nu, &eps1, &eps2, &chord1, &chord2);
			}

			//computing vector along LE of current spanwise row of DVEs
			tempS = 1./panelPtr[i].n;
			xLE[0] = (x2[0]-x1[0])*tempS;
			xLE[1] = (x2[1]-x1[1])*tempS;
			xLE[2] = (x2[2]-x1[2])*tempS;

			//computing dihedral of LE of current spanwise row of DVEs
			tempS	= sqrt(xLE[1]*xLE[1]+xLE[2]*xLE[2]);
			nu2 	= asin(xLE[2]/tempS);

			//loop over number of spanwise elements 'n'
			for (n=0;n<panelPtr[i].n;n++)
			{
				tempS = (0.5+n);

				if(info.flagCAMBER){
					//half-chord length at midspan of DVE using camber info
					surfacePtr[l].xsi=0.5*(chord1+tempS*(chord2-chord1)/panelPtr[i].n);
					//temporary incidence angle at half span of DVE using camber info
					surfacePtr[l].epsilon = eps1 + tempS*(eps2-eps1)/panelPtr[i].n;
				} else{
					//half-chord length at midspan of DVE
					surfacePtr[l].xsi=0.5*(delchord1+delchord*tempS); 
					//temporary incidence angle at half span of DVE
					surfacePtr[l].epsilon = panelPtr[i].eps1 + tempS*deleps;
				}

				ceps  = cos(surfacePtr[l].epsilon); //needed for tempA below
				seps  = sin(surfacePtr[l].epsilon);

				//computing the mid-span location of the DVE LE
				//temporarily saved in xleft
				surfacePtr[l].xleft[0] = x1[0]+xLE[0]*tempS;
				surfacePtr[l].xleft[1] = x1[1]+xLE[1]*tempS;
				surfacePtr[l].xleft[2] = x1[2]+xLE[2]*tempS;

				//vector from LE mid-span point to reference point
				delX[0] =  surfacePtr[l].xsi * ceps;
				delX[1] =  surfacePtr[l].xsi * seps * sin(nu2);
				delX[2] = -surfacePtr[l].xsi * seps * cos(nu2);

				//computing the DVE reference point
				surfacePtr[l].xo[0] = surfacePtr[l].xleft[0] + delX[0];
				surfacePtr[l].xo[1] = surfacePtr[l].xleft[1] + delX[1];
				surfacePtr[l].xo[2] = surfacePtr[l].xleft[2] + delX[2];

				//computing the normal of DVE,
				cross(delX,xLE,surfacePtr[l].xleft);
				tempS = 1/norm2(surfacePtr[l].xleft);
				scalar(surfacePtr[l].xleft,tempS,surfacePtr[l].normal);

				//computes local free stream velocity at xo location
				scalar(panelPtr[i].u2,tempS,tempA);		//u2/n*(.5+k)
				scalar(panelPtr[i].u1,(1-tempS),tempAA);//-u1(1/n*(.5+k)-1)
				vsum(tempA,tempAA,surfacePtr[l].u);		//u2+(u2-u1)/n*(.5+k)

				//the normalized free stream direction, G.B. 10/19/04
				tempS=1./norm2(surfacePtr[l].u);
				surfacePtr[l].U[0] = surfacePtr[l].u[0]*tempS;
				surfacePtr[l].U[1] = surfacePtr[l].u[1]*tempS;
				surfacePtr[l].U[2] = surfacePtr[l].u[2]*tempS;

				//*************************************************************
				//dihedral angle, nu, tan(nu)=-(ny)/(nz)
				if(fabs(surfacePtr[l].normal[2]) > DBL_EPS)
				{ tempS=-atan(surfacePtr[l].normal[1]/surfacePtr[l].normal[2]);
				  if (surfacePtr[l].normal[2] < 0)  tempS += Pi;//|nu|>Pi/2
				}
				else  //tan(nu) -> infinity  -> |nu| = Pi/2
				{
				  if(surfacePtr[l].normal[1]>0) tempS = -0.5*Pi;
				  else 						 	tempS =  0.5*Pi;
				}
				surfacePtr[l].nu = tempS;
				//*************************************************************

				//computing epsilon
				surfacePtr[l].epsilon = asin(surfacePtr[l].normal[0]);

				//computing psi, set to zero for surface DVEs
				surfacePtr[l].psi= 0;

				//computing yaw angle psi
				//psi is the angle between vector from LE to Xo and xsi-axis
					//transforming xLE into local reference frame
					Glob_Star(delX,surfacePtr[l].nu,surfacePtr[l].epsilon,\
						surfacePtr[l].psi,tempA);
										//function in ref_frame_transform.h
				surfacePtr[l].psi = atan(tempA[1]/tempA[0]);

				//transforming xLE into local reference frame
				Glob_Star(xLE,surfacePtr[l].nu,surfacePtr[l].epsilon,\
					  	  surfacePtr[l].psi,xsiLE);
					  				//function in ref_frame_transform.h

				//DVE half span, eta
				surfacePtr[l].eta = 0.5*xsiLE[1];

				//GB 2-14-20
				//interpolation ratio for airfoil and camber, =0 at panel1 and =1 at panel2
				//ratio = [span location of X0]/[panel span] reduces to:
				surfacePtr[l].ratio = (n+0.5)/panelPtr[i].n;

				//DVE leading edge sweep
				tempS = xsiLE[0]/xsiLE[1]; //tan(phi)=xsi/eta
				surfacePtr[l].phiLE = atan(tempS);

				//DVE trailing edge sweep
				surfacePtr[l].phiTE = atan(tempS+delTANphi);

				//DVE mid-chord sweep
				surfacePtr[l].phi0 = 0.5 * \
									(surfacePtr[l].phiLE+surfacePtr[l].phiTE);

 				//DVE area 4*eta*xsi  G.B. 8-10-07
				surfacePtr[l].S = 4* surfacePtr[l].eta*surfacePtr[l].xsi;

				//DVE assign airfoil number to element
				//two airfoils to interpolate between panel edge 1 and 2
				//GB 2-14-20
				surfacePtr[l].airfoil[0] = panelPtr[i].airfoil1;
				surfacePtr[l].airfoil[1] = panelPtr[i].airfoil2;
				

//printf("phiLE %2.4lf phiTE %2.4lf phi0 %2.4lf \n",\
//surfacePtr[l].phiLE*RtD,surfacePtr[l].phiTE*RtD,surfacePtr[l].phi0*RtD);

				l++;		//next surface-DVE
			}//END loop over number of spanwise elements 'n'
		}//loop over number of chordwise elements 'panel.m'
	}//END loop over number of panels

	//##this part has been added 2/9/05 G.B.
	//computing decaying factor for added singularity at wing tip
	int k=0;			//k => trailing-edge element counter

	if(info.nowing>1)  	//more than one wing and possible interaction
	{					//between surface and wake
		for(wing=0;wing<info.nowing;wing++)
		{
			//index of last DVE of this wing (located at tip and trail. edge)
			//removed GB 2-9-20 l = k + (info.wing2[wing]-info.wing1[wing]+1)*info.m - 1;
            l=info.dve2[wing];
            
			//is the wing symmetrical or not?
			if(info.sym == 1)	//decay factor is 1% of tip-element half-span
				singfct = 0.01*surfacePtr[l].eta;
			else//wing has two tips, possibly different in geometry
			{	//in that case, decay factor is 1% of the shorter half-span
				if(  surfacePtr[k].eta < surfacePtr[l].eta)
							singfct = 0.01*surfacePtr[k].eta;
				else 		singfct = 0.01*surfacePtr[l].eta;
			}

			//loop over surface DVEs of current wing
			for(i=k;i<=l;i++)
			surfacePtr[i].singfct = singfct;  //assigning decay factor

			k = l+1;	//updating index for next wing
		}//next wing
	}
	else	//if more than one wing and the potential of
	{		//wake-surface interaction exists; 8/16/05 G.B.
		for(i=0;i<info.noelement;i++)	surfacePtr[i].singfct = 0.;
	}

	//computing points that are at right and left edge oof leading edge of element
	//This section was added for FlexWing on 10-29-06 G.B.
	void Edge_Point(const double [3],const double,const double,const double,\
				const double,const double,const double,double [3]);
									//found in wake_geometry.cpp
	for (i=0;i<info.noelement;i++)
	{
		//computes point at left side of ledge edge of DVE
		Edge_Point(surfacePtr[i].xo,surfacePtr[i].nu,\
				   surfacePtr[i].epsilon,surfacePtr[i].psi,\
				   surfacePtr[i].phiLE,-surfacePtr[i].eta,\
				   -surfacePtr[i].xsi,surfacePtr[i].x1);
				 					//subroutine in wake_geometry.cpp

		//computes point at right side of ledge edge of DVE
		Edge_Point(surfacePtr[i].xo,surfacePtr[i].nu,\
				   surfacePtr[i].epsilon,surfacePtr[i].psi,\
				   surfacePtr[i].phiLE,surfacePtr[i].eta,\
				   -surfacePtr[i].xsi,surfacePtr[i].x2);
				 					//subroutine in wake_geometry.cpp

		//initializing velocities in x1 and x2 with the one found in xo
		surfacePtr[i].u1[0]=surfacePtr[i].u[0];
		surfacePtr[i].u1[1]=surfacePtr[i].u[1];
		surfacePtr[i].u1[2]=surfacePtr[i].u[2];

		surfacePtr[i].u2[0]=surfacePtr[i].u[0];
		surfacePtr[i].u2[1]=surfacePtr[i].u[1];
		surfacePtr[i].u2[2]=surfacePtr[i].u[2];

		//computing point at center of trailing edge xTE
		Edge_Point(surfacePtr[i].xo,surfacePtr[i].nu,\
				   surfacePtr[i].epsilon,surfacePtr[i].psi,\
				   surfacePtr[i].phiLE,0,\
				   surfacePtr[i].xsi,surfacePtr[i].xTE);
				 					//subroutine in wake_geometry.cpp

		//computing vector along trailing edge TEvc
			//right edge point of trailinge edge
		Edge_Point(surfacePtr[i].xo,surfacePtr[i].nu,\
				   surfacePtr[i].epsilon,surfacePtr[i].psi,\
				   surfacePtr[i].phiLE,surfacePtr[i].eta,\
				   surfacePtr[i].xsi,tempA);
				 					//subroutine in wake_geometry.cpp
		surfacePtr[i].TEvc[0] = tempA[0]-surfacePtr[i].xTE[0];
		surfacePtr[i].TEvc[1] = tempA[1]-surfacePtr[i].xTE[1];
		surfacePtr[i].TEvc[2] = tempA[2]-surfacePtr[i].xTE[2];
	}

}
//===================================================================//
		//END FUNCTION Surface_DVE_Generation
//===================================================================//

//===================================================================//
		//START FUNCTION Move_Flex_Wing
//===================================================================//
void Move_Flex_Wing(const GENERAL info, DVE* surfacePtr)
{
//new routine added 10-29-06 G.B.
//moves flexible wing by delta x every time step,
//1. moves xo, x1, x2
//2. computes new eta nu, epsilon, psi, xsi of DVE
//3. updates singularity decay factor
//4. updates circulation coefficients as DVE stretches in span
//
//
//	input
//	info			general information
//	surfacePtr		surface DVE's
//

int i,wing,k,l;		//counters
double delx[3];		//increment x1,xo, and x2 move during timestep
double delX1[3];	//vector between center of leading edge and ref. pt
double delX2[3];	//vector between left and right leading edge
double delXSI1[3];	//delX1 in local DVE coordinates
double delXSI2[3];	//delX2 in local DVE coordinates
double delPhi;		//change in sweep
double singfct;		//temporary decay factor
double tempA[3];



//####################################################################
//1. moves xo, x1, x2
	for(i=0;i<info.noelement;i++)
	{
		//moves reference point
		//delta x = local U * delta time
		delx[0] = surfacePtr[i].u[0] * info.deltime;
		delx[1] = surfacePtr[i].u[1] * info.deltime;
		delx[2] = surfacePtr[i].u[2] * info.deltime;

		//move reference point
		surfacePtr[i].xo[0] -= delx[0];
		surfacePtr[i].xo[1] -= delx[1];
		surfacePtr[i].xo[2] -= delx[2];

		//moves point at left edge of DVE-leading edge
		//delta x = local U * delta time
		delx[0] = surfacePtr[i].u1[0] * info.deltime;
		delx[1] = surfacePtr[i].u1[1] * info.deltime;
		delx[2] = surfacePtr[i].u1[2] * info.deltime;

		//move reference point
		surfacePtr[i].x1[0] -= delx[0];
		surfacePtr[i].x1[1] -= delx[1];
		surfacePtr[i].x1[2] -= delx[2];

		//moves point at right edge of DVE-leading edge
		//delta x = local U * delta time
		delx[0] = surfacePtr[i].u2[0] * info.deltime;
		delx[1] = surfacePtr[i].u2[1] * info.deltime;
		delx[2] = surfacePtr[i].u2[2] * info.deltime;

		//move reference point
		surfacePtr[i].x2[0] -= delx[0];
		surfacePtr[i].x2[1] -= delx[1];
		surfacePtr[i].x2[2] -= delx[2];
	}
//####################################################################

//2. computes new eta nu, epsilon, psi, xsi of DVE
	//routine similar to the one found in wake_geometry.cpp for relaxing
	for(i=0;i<info.noelement;i++)
	{
	//********************************************************************
		//computing vector from center of leading edge to reference point
		delX1[0] = surfacePtr[i].xo[0] \
					- 0.5*(surfacePtr[i].x1[0]+surfacePtr[i].x2[0]);
		delX1[1] = surfacePtr[i].xo[1] \
					- 0.5*(surfacePtr[i].x1[1]+surfacePtr[i].x2[1]);
		delX1[2] = surfacePtr[i].xo[2] \
					- 0.5*(surfacePtr[i].x1[2]+surfacePtr[i].x2[2]);

		//vector along leading edge of DVE
		delX2[0] = surfacePtr[i].x2[0] - surfacePtr[i].x1[0];
		delX2[1] = surfacePtr[i].x2[1] - surfacePtr[i].x1[1];
		delX2[2] = surfacePtr[i].x2[2] - surfacePtr[i].x1[2];

		//computing the normal of DVE,
		cross(delX1,delX2,tempA);
		scalar(tempA,1/norm2(tempA),surfacePtr[i].normal);

		//here begins the stuff you been waiting for:

	//********************************************************************
		//dihedral angle, nu,
		//tan(nu)=-(ny)/(nz)
		if(fabs(surfacePtr[i].normal[2]) > DBL_EPS)
		{
			surfacePtr[i].nu=-atan(surfacePtr[i].normal[1]\
									/surfacePtr[i].normal[2]);
			if (surfacePtr[i].normal[2]<0)  surfacePtr[i].nu+= Pi;//|nu|>Pi/2
		}
		else  //tan(nu) -> infinity  -> |nu| = Pi/2
		{
			if(surfacePtr[i].normal[1]>0) 	surfacePtr[i].nu = -0.5*Pi;
			else	 						surfacePtr[i].nu =  0.5*Pi;
		}
	//********************************************************************
		//computing epsilon
		surfacePtr[i].epsilon = asin(surfacePtr[i].normal[0]);

	//********************************************************************
		//computing yaw angle psi
			//the orientation of the rotational axis of the vortex seet,
			//is parallel to delX1.

		//delX1 in xsi reference frame (psi = 0)
		Glob_Star(delX1,surfacePtr[i].nu,surfacePtr[i].epsilon,0,delXSI1);
						  				//function in ref_frame_transform.h

		//tan(psi)=(eta2)/(xsi1)
		if(delXSI1[0]*delXSI1[0] > DBL_EPS)
		{
			surfacePtr[i].psi= atan(delXSI1[1]/delXSI1[0]);
				if(delXSI1[0] < 0)  surfacePtr[i].psi += Pi;//|nu|>Pi/2
		}
		else	//tan(psi) -> infinity  -> |psi| = Pi/2
			surfacePtr[i].psi = 0.5*Pi*fabs(delXSI1[1])/delXSI1[1];

		//length of projection
		surfacePtr[i].xsi = sqrt(delXSI1[0]*delXSI1[0]+delXSI1[1]*delXSI1[1]);

	//********************************************************************
		//computing sweep and halfspan

		//1. transform delX2 into local frame
		Glob_Star(delX2,surfacePtr[i].nu,surfacePtr[i].epsilon,\
											surfacePtr[i].psi,delXSI2);
						  				//function in ref_frame_transform.h

		//if everything works => delXSI[2] = 0
		if(delXSI2[1] < 0)
		{
			printf(" ohwei! around line 765, wing_geometry.cpp\n %2.2lf %2.16lf\n",delXSI2[1],delXSI2[2]);
					exit(0);
		}

		//computing change in sweep
		delPhi = atan(delXSI2[0]/delXSI2[1]) - surfacePtr[i].phiLE;

		//updating sweeps
		surfacePtr[i].phiLE += delPhi;
		surfacePtr[i].phi0  += delPhi;
		surfacePtr[i].phiTE += delPhi;

		//computing the new half span
		surfacePtr[i].eta = 0.5*delXSI2[1];
	//********************************************************************
	}
//####################################################################

//3. updates singularity decay factor
	//computing decaying factor for added singularity at wing tip
	k=0;			//k => trailing-edge element counter

	if(info.nowing>1)  	//more than one wing and possible interaction
	{					//between surface and wake
		for(wing=0;wing<info.nowing;wing++)
		{
			//index of last DVE of this wing (located at tip and trail. edge)
            //removed GB 2-9-20 l = k + (info.wing2[wing]-info.wing1[wing]+1)*info.m - 1;
            l=info.dve2[wing];

			//is the wing symmetrical or not?
			if(info.sym == 1)	//decay factor is 1% of tip-element half-span
				singfct = 0.01*surfacePtr[l].eta;
			else//wing has two tips, possibly different in geometry
			{	//in that case, decay factor is 1% of the shorter half-span
				if(  surfacePtr[k].eta < surfacePtr[l].eta)
							singfct = 0.01*surfacePtr[k].eta;
				else 		singfct = 0.01*surfacePtr[l].eta;
			}

			//loop over surface DVEs of current wing
			for(i=k;i<=l;i++)
			surfacePtr[i].singfct = singfct;  //assigning decay factor

			k = l+1;	//updating index for next wing
		}//next wing
	}
	else	//if more than one wing and the potential of
	{		//wake-surface interaction exists; 8/16/05 G.B.
		for(i=0;i<info.noelement;i++)	surfacePtr[i].singfct = 0.;
	}
//####################################################################

//4. updates circulation coefficients as DVE stretches in span
//not needed since new circulation distribution is computed based on
//newly squirted out wake row
//####################################################################

}
//===================================================================//
		//END FUNCTION Move_Flex_Wing
//===================================================================//

//===================================================================//
		//START FUNCTION Move_Wing
//===================================================================//
void Move_Wing(const GENERAL info, DVE* surfacePtr)
{
//moves wing by delta x every time step,
//function updates xo location of surface DVE's
//
//	input
//	info			general information
//	surfacePtr		surface DVE's
//

int i;
double delx[3];

	for(i=0;i<info.noelement;i++)
	{
		//delta x = local U * delta time
		delx[0] = surfacePtr[i].u[0] * info.deltime;
		delx[1] = surfacePtr[i].u[1] * info.deltime;
		delx[2] = surfacePtr[i].u[2] * info.deltime;

		//move reference point
		surfacePtr[i].xo[0] -= delx[0];
		surfacePtr[i].xo[1] -= delx[1];
		surfacePtr[i].xo[2] -= delx[2];
	}
}
//===================================================================//
		//END FUNCTION Move_Wing
//===================================================================//


//===================================================================//
		//START FUNCTION Apply_Camber 
//===================================================================//
void Apply_Camber(const PANEL* panelPtr, double x1[3], double x2[3], \
	double ***camberPtr, int m, int i, double nu, \
	double *eps1, double *eps2, double *chord1, double *chord2)
{
	// The function Apply_Camber determines changed the LE panel points to follow
	// 		the curvature of the camber line.
	//
	// Function inputs:
	//		PANEL* panelPtr -Panel structure for the panelPtr
	//		x1,x2 			-x,y,z position of the left (1) and right (2) panel LE pts
	//		camberptr 		-camber ptr holding all of the camber data
	//		m 				-chorwise row of interest
	//		i 				-panelPtr of interest
	//		nu 				-dihedral angle of panel
	//		eps1, eps2 		-(see output) 
	//		chord1, chord2 	-(see output)
	//
	// Function outputs:
	//		x1,x2 			-updated left and right panel LE points
	//		eps1, eps2 		-updated left and right edge epsilon angles with camber
	//		chord1,chord2 	-updated left and right edge chord lengths with camber
	//
	// NOTE (1) indicated left and (2) indicated right for the above variables

	// D.F.B. in Braunschweig, Germany, Feb. 2020

	int j; 	//Generic counter
	double tempZ1, tempZ2; 		//Z LE offsets
	double tempZTE1, tempZTE2;	//Z TE offsets
	double tempx1[3],tempx2[3];	//Local reference frame LE pts
	double leftZ, rightZ;		//Delta Z from LE to TE

	// Bring the inputted x1 and x2 into the local DVE ref frame
	Glob_Star(x1,nu,panelPtr[i].eps1,0,tempx1);
	Glob_Star(x2,nu,panelPtr[i].eps2,0,tempx2);

	//First consider LE left side
	// Seach for index where (current m)/(total m) is nearest the camber data
	j = 0;
	do{j++;}
	while(camberPtr[panelPtr[i].airfoil1][j][0]<(double(m)/double(panelPtr[i].m)));

	//Calculate the z/c by linearly interpolate the using the above define index 
	tempZ1 = camberPtr[panelPtr[i].airfoil1][j-1][1] + ((double(m)/double(panelPtr[i].m)-camberPtr[panelPtr[i].airfoil1][j-1][0])*\
			(camberPtr[panelPtr[i].airfoil1][j][1]-camberPtr[panelPtr[i].airfoil1][j-1][1])/\
			(camberPtr[panelPtr[i].airfoil1][j][0]-camberPtr[panelPtr[i].airfoil1][j-1][0]));
	//Scale z/c from camber data to the chordlength of the edge
	tempx1[2] +=(tempZ1*panelPtr[i].c1); 

	//Repeat for TE left side
	j = 0;
	do{j++;}
	while(camberPtr[panelPtr[i].airfoil1][j][0]<(double(m+1)/double(panelPtr[i].m)));

	tempZTE1= camberPtr[panelPtr[i].airfoil1][j-1][1] + ((double(m+1)/double(panelPtr[i].m)-camberPtr[panelPtr[i].airfoil1][j-1][0])*\
			(camberPtr[panelPtr[i].airfoil1][j][1]-camberPtr[panelPtr[i].airfoil1][j-1][1])/\
			(camberPtr[panelPtr[i].airfoil1][j][0]-camberPtr[panelPtr[i].airfoil1][j-1][0]));

	//Repeat for LE right side
	j = 0;
	do{j++;}
	while(camberPtr[panelPtr[i].airfoil2][j][0]<(double(m)/double(panelPtr[i].m)));

	tempZ2 = camberPtr[panelPtr[i].airfoil2][j-1][1] + ((double(m)/double(panelPtr[i].m)-camberPtr[panelPtr[i].airfoil2][j-1][0])*\
			(camberPtr[panelPtr[i].airfoil2][j][1]-camberPtr[panelPtr[i].airfoil2][j-1][1])/\
			(camberPtr[panelPtr[i].airfoil2][j][0]-camberPtr[panelPtr[i].airfoil2][j-1][0]));
	tempx2[2] +=(tempZ2*panelPtr[i].c2); 

	//Repeat for TE right side
	j = 0;
	do{j++;}
	while(camberPtr[panelPtr[i].airfoil2][j][0]<(double(m+1)/double(panelPtr[i].m)));

	tempZTE2 = camberPtr[panelPtr[i].airfoil2][j-1][1] + ((double(m+1)/double(panelPtr[i].m)-camberPtr[panelPtr[i].airfoil2][j-1][0])*\
		(camberPtr[panelPtr[i].airfoil2][j][1]-camberPtr[panelPtr[i].airfoil2][j-1][1])/\
		(camberPtr[panelPtr[i].airfoil2][j][0]-camberPtr[panelPtr[i].airfoil2][j-1][0]));


	//Calculate delta z from LE to TE on left edge and on right edge
	leftZ = (tempZ1-tempZTE1)*panelPtr[i].c1;
	rightZ =(tempZ2-tempZTE2)*panelPtr[i].c2;

	//New section chord on left and right edges
	*chord1 = sqrt(leftZ*leftZ+(panelPtr[i].c1/panelPtr[i].m)*(panelPtr[i].c1/panelPtr[i].m));
	*chord2 = sqrt(rightZ*rightZ+(panelPtr[i].c2/panelPtr[i].m)*(panelPtr[i].c2/panelPtr[i].m));

	// Calculate mid-span eps
	*eps1 = atan(leftZ/(panelPtr[i].c1/panelPtr[i].m));
	*eps2 = atan(rightZ/(panelPtr[i].c2/panelPtr[i].m));


	// Convert new LE points to the global reference frame 
	Star_Glob(tempx1,nu,panelPtr[i].eps1,0,x1);
	Star_Glob(tempx2,nu,panelPtr[i].eps2,0,x2);

}
//===================================================================//
		//END FUNCTION Move_Wing
//===================================================================//