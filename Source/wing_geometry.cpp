//generates trailing edge vortex elements of a panel
void Trailing_Edge_Generation\
					(const PANEL, const GENERAL, const BOUND_VORTEX*,\
					 BOUND_VORTEX *, int &);
//identifies edges of separate wings
void Wing_Generation(const PANEL*,int,int[5],int[5],int[5],int[5],\
                     int[5],int[5]);
//Rotates panels for attitude during turning flight
void Panel_Rotation(GENERAL &,PANEL *);
//generates surface DVE elements
void Surface_DVE_Generation(GENERAL &,PANEL *,DVE *,double ***,\
							const double);
//moves wing by delta x every time step,
void Move_Wing(const GENERAL, DVE*, const double[3],double[3]);
//moves flexible wing by delta x every time step,
void Move_Flex_Wing(const GENERAL, DVE*);
//moves DVEs according do the input camber line
void Apply_Camber(const PANEL*, double[3], double[3] , double[3], double[3],double ***, int, \
	int, double,double *,double *,double *,double *);
//deflects DVEs about the hinge
void DeflectAboutHinge(const PANEL*, const double, double[3], double[3], double[3], double[3], \
	int, int, double, double *,double *m,double[3], double[3]);
//calculates U_inf for circling flight
void Circling_UINF(GENERAL, DVE*,double[3], const double [3]);
// Re-calculated DVE parameters based on DVE LE pts and Control pts
void DVE_LEandCPtoParam(GENERAL, DVE*);

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
 	if(tempA1[1]<0) nu = Pi-nu; 
	//Added to allows for panel to be defined from right to left and dihedral over 90deg 
	//##GB 5/20/21 

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
	if(xte[1]<0) nu = Pi-nu;  
	//Added to allows for panel to be defined from right to left and dihedral over 90deg 
    //##GB 5/20/21 

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

int k,span,wing,index,panel;		//loop counters, k=0..(panel.n-1)
    
    panel1[0]=0;    panel2[0]=0;  //first left panel of wing 1

    for(wing=0;wing<info.nowing;wing++)
    {
        k=panel1[wing]+1;

        //identifies whether panel to left or right as defined in input file is
        //of lesser or equal index as current index and not no-neighober (-1).
        //If so, advance to next panel.
        while(((panelPtr[k].left-1<=panel2[wing] && panelPtr[k].left-1 !=-1) \
          ||  (panelPtr[k].right-1<=panel2[wing] && panelPtr[k].right-1 !=-1)) \
			&& k < info.nopanel)

        {
            panel2[wing]=k;
            panel1[wing+1]=k+1;
            panel2[wing+1]=k+1;
            k++;
        }
    }



	//identify DVE indices
    span=0; wing=0; index=0; //initialize
    
    for(wing=0;wing<info.nowing;wing++) //loop over wings
    {
        wing1[wing] = span;
        dve1[wing] = index;
        //loop over panels of wing
        for(panel=panel1[wing];panel <= panel2[wing];panel++)
        {
            span += panelPtr[panel].n-1;
            index += panelPtr[panel].n*panelPtr[panel].m-1;
            
            wing2[wing] = span;
            dve2[wing] = index;
            
            index++;
            span++;
        }// next panel
    }// next wing

/*  //test output
    for(wing=0;wing<info.nowing;wing++) //loop over wings
    {
        printf("wing %d panel1 %d panel2 %d\n",wing,panel1[wing],panel2[wing]);
    }
// */
    
}
//===================================================================//
		//END FUNCTION Wing_Generation
//===================================================================//

//===================================================================//
		//FUNCTION 101
//===================================================================//
void Panel_Rotation(GENERAL &info,PANEL* panelPtr)
{
//Function important for turning flight simulation.
//Function rotates panels to account for aircraft sideslip, roll angles
//and angle of attack. 
//rotation by alpha on only considered when turning flight
//happens in during level flight.
	//
	//input:
	//	info 		- uses alpha, beta and phi to rotate panels and Ref. Pt
	// 	panel		- information of panels
	//
	//ouput:
    //  I. 	 updated panel geometry (move x1 & x2);
    //  Ia.  update epsilon of panel if horizontal flight
    //  II.     updated CG location (RefPt)
    //  III. adjust beta and if horizontal flight alpha
    //  IV.  adjust freestream vector (info.U), also to include upwind info.Ws
	//
	//1. rotation about x-axis by phi; positive if left wing down
	//2. rotarion about z'-axis; positive nose to left of flow
	//3. rotation about y"; only if in horizontal plane; posit. nose up
	//
	//rotateX, rotateY, ratateZ found in vector_algebra.h

	int panel; //panel counter
	double tempA[3],tempAA[3], tempS; //temporary arrays and scalar
    double eps,psi,nu,wingNu;
    

//    I.      updated panel geometry (move x1 & x2);

    if(info.flagHORZ)   eps = info.alpha;
    else eps = 0;
    psi = info.beta;
    nu = info.bank;
    
	for(panel=0;panel<info.nopanel;panel++)
	{
        //  Ia.  update epsilon of panel if horizontal flight before rotation of X1 and X2
        if(info.flagHORZ)
        {
            {   //finding leading edge dihedral (because if vertical surface, no adjustment needed)
                //code block taken from FUNCTION Surface_DVE_Generation
                //xquart: vector along leading edge of panel
                scalar(panelPtr[panel].x1,-1,tempA);    //negates x1
                vsum(panelPtr[panel].x2,tempA,tempAA); //x_l.e. = x1/4_2-x1/4_1

                //panel span
                tempS = sqrt(tempAA[1]*tempAA[1]+tempAA[2]*tempAA[2]);

                //leading edge dihedral
                wingNu = asin(tempAA[2]/tempS);
                if(tempAA[1]<0) wingNu = Pi-wingNu; 
                //Added to allows for panel to be defined from right to left and dihedral over 90deg 
                //##GB 5/20/21 
            }
            //correcting epsilons for alpha and beta considering dihedral (wingNu)
            panelPtr[panel].eps1 += eps*cos(wingNu) + psi*sin(wingNu);
            panelPtr[panel].eps2 += eps*cos(wingNu) + psi*sin(wingNu);
        } //END if(info.flagHORZ)
        
        //rotate X1
        rotateZ(panelPtr[panel].x1,psi,tempA); //%note! this is beta, NOT yaw for circ. flight. 
		rotateY(tempA, -eps, tempAA);
        rotateX(tempAA,-nu, panelPtr[panel].x1);


		//rotate X2
        rotateZ(panelPtr[panel].x2,psi,tempA);
		rotateY(tempA, -eps, tempAA);
		rotateX(tempAA, -nu, panelPtr[panel].x2);
    }//END loop over panels

//  II.     updated CG location (RefPt)
	//rotate reference point (CG)
	rotateZ(info.RefPt,psi,tempA);
	rotateY(tempA, -eps, tempAA);
	rotateX(tempAA, -nu, info.RefPt);

//  III. adjust beta and if horizontal flight alpha
    info.beta = 0;
    if(info.flagHORZ)   info.alpha=0;

    
//  IV.  adjust freestream vector (info.U), also to include upwind info.Ws
	info.U[0] = info.Uinf*cos(info.alpha)*cos(info.beta); //this is never used.
	info.U[1] = info.Uinf            *sin(info.beta); //this is never used. 
    info.U[2]= info.Uinf*sin(info.alpha)*cos(info.beta)-info.Ws; //this is used to add to the omega velocities in circUINF

		
}
//===================================================================//
		//END FUNCTION Panel_Rotation
//===================================================================//


//===================================================================//
		//FUNCTION Surface_DVE_Generation
//===================================================================//
void Surface_DVE_Generation(GENERAL &info,PANEL* panelPtr,\
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
double tempS, tempSpan,tempA[3],tempAA[3]; //temporary variables. 
double xquart[3];					//1/4chord line of panel, defined in input
double x1LE[3],x2LE[3],x1[3],x2[3]; //edge points of panel and spanwise rows
double xLE[3],xsiLE[3];	//vector along LE of a spanewise row of surface DVEs
double delchord1,delchord2,delchord; //edge 1&2 chord,chord span-increment,
double deleps,ceps,seps;//increments of epsi; cos/sin of incidence angle
double nu,nu2;		//nus of panel 1/4c line, LE of DVE row
double delTANphi;	//tan(phiLE-phiTE)
double delX[3];		//vector from center of leading edge to control point
double eps1,eps2;	//epsilon of chordwise section (needed for camber/hinges)
double epsH1,epsH2;	//epsilon of chordwise section of hinge
double epsC1,epsC2; //epsilon of chordwise section of camber
double chord1,chord2;	//new chord due to camber
double xH1[3],xH2[3];	//hinge location 
double check1, check2;		//Checks to see if hinge is a DVE LE
int closeLL =0;			//closest lifting line to hinge
double adjust1, adjust2;	//distance to shift lifting lines if they are not on a hinge
double adjust1next, adjust2next;	//distance to shift lifting lines if they are not on a hinge
char answer;	//error handling

//init area calculations
info.AREA = 0;
info.projAREA = 0;
info.projSPAN = 0;
info.surfAREA = 0;

	//loop over number of panels
	for (i=0;i<info.nopanel;i++)
	{
		//xquart: vector between leading edge corners
		scalar(panelPtr[i].x1,-1,tempA);	//negates x1
		vsum(panelPtr[i].x2,tempA,xquart); //x_l.e. = x1/4_2-x1/4_1

		//panel span
		tempSpan = sqrt(xquart[1]*xquart[1]+xquart[2]*xquart[2]); 

		//panel area
		panelPtr[i].AREA = tempSpan * (panelPtr[i].c1 + panelPtr[i].c2) / 2;
		info.AREA += panelPtr[i].AREA;

		//panel area projection to xy plane
		panelPtr[i].projAREA = xquart[1] * \
            (panelPtr[i].c1 * cos(panelPtr[i].eps1) + panelPtr[i].c2 * cos(panelPtr[i].eps2)) / 2;
		info.projAREA += panelPtr[i].projAREA;
        info.projSPAN += xquart[1]; //projected span

		//1/4chord line dihedral
		nu = asin(xquart[2]/ tempSpan);
		if(xquart[1]<0) nu = Pi-nu; //Added to allows for panel to be defined from right to left

		//note: delchord is only temporary and may be changed as elements are redistributed
		//below
		delchord1 = panelPtr[i].c1/double(panelPtr[i].m);	//chordwise increment left side
		delchord2 = panelPtr[i].c2/double(panelPtr[i].m);	//chordwise increment right side
		//delchord  = (delchord2-delchord1)/panelPtr[i].n;//chord/span increment

		//tangent of change in sweep angle of each chordwise row
		//removed BB 2020, replaces below for each m
		//delTANphi = (delchord2-delchord1)/tempSpan;

		//spanwise incident increments
		//deleps = (panelPtr[i].eps2-panelPtr[i].eps1)/panelPtr[i].n; //removed BB 2020, replaced with local eps calcs

		/* Changed 20.02.20 D.F.B. User now defined panel at LE and the wing
		/ pitch is adjusted about the LE of the wing (instead of the 1/4 chord)
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


		panelPtr[i].surfAREA = 0; //init panel surf area.

		//should throw an error here if we have a hinge but m = 1
		if ((panelPtr[i].hinge1 != 0.0 || panelPtr[i].hinge2 != 0.0) && panelPtr[i].m == 1) {
			//in this case we have defined a hinge, but m=1, so we will raise an error
			printf("You have defined a hinge, but this is impossible with m=1\n");
			printf("Either increase m or set the hinge location to 0.0\n");
			printf("---Exiting program---\n");
			scanf("%c", &answer);
			exit(0);
		}

		//determine the closest lifting line to the hinge, and how far we have to move the lifting lines
		//to have them line up with the hinge
		adjust1 = 0;
		adjust2 = 0;
		adjust1next = 0;
		adjust2next = 0;
		check1 = (panelPtr[i].m * panelPtr[i].hinge1 - round(panelPtr[i].m * panelPtr[i].hinge1));
		check2 = (panelPtr[i].m * panelPtr[i].hinge2 - round(panelPtr[i].m * panelPtr[i].hinge1));
		if (fabs(check1) > DBL_EPS || fabs(check2) > DBL_EPS) 
		{
			//hinge on either side of panel not at lifting line
			//this will be the lifting line which is the closest to the hinge
			//at the midspan of the panel  ##changed from left edge GB 1-23-2
	//		closeLL = round(panelPtr[i].m * panelPtr[i].hinge1);
			closeLL = round(panelPtr[i].m * (panelPtr[i].hinge1+panelPtr[i].hinge2)*0.5);
			if (closeLL == 0) closeLL = 1; //cannot move LE
			else if (closeLL == panelPtr[i].m) closeLL = panelPtr[i].m - 1; //cannot move TE
			//how far to move the LL to have it be at the hinge:
			//adjust1 = ((double(closeLL) / panelPtr[i].m) - panelPtr[i].hinge1) * -panelPtr[i].c1;
			//adjust2 = ((double(closeLL) / panelPtr[i].m) - panelPtr[i].hinge2) * -panelPtr[i].c2;
		}


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

			//need to move lifting lines if checks failed:
			if (fabs(check1) > DBL_EPS || fabs(check2) > DBL_EPS) {
				if (m < closeLL) {
					//how far to move the LL:
					adjust1 = (((panelPtr[i].hinge1 / double(closeLL)) * double(m)) - (double(m) / double(panelPtr[i].m))) * panelPtr[i].c1;
					adjust2 = (((panelPtr[i].hinge2 / double(closeLL)) * double(m)) - (double(m) / double(panelPtr[i].m))) * panelPtr[i].c2;

					adjust1next = (((panelPtr[i].hinge1 / double(closeLL)) * (double(m) + 1.)) - ((double(m) + 1.) / double(panelPtr[i].m))) * panelPtr[i].c1;
					adjust2next = (((panelPtr[i].hinge2 / double(closeLL)) * (double(m) + 1.)) - ((double(m) + 1.) / double(panelPtr[i].m))) * panelPtr[i].c2;

					//new chord:
					chord1 = (delchord1 * (double(m) + 1.) + adjust1next) - ((delchord1 * double(m)) + adjust1);
					chord2 = (delchord2 * (double(m) + 1.) + adjust2next) - ((delchord2 * double(m)) + adjust2);

				}

				else { //m>closeLL
					adjust1 = ((1.0 - (((1.0 - panelPtr[i].hinge1) / (double(panelPtr[i].m) - double(closeLL))) * (double(panelPtr[i].m) - double(m)))) - (double(m) / double(panelPtr[i].m))) * panelPtr[i].c1;
					adjust2 = ((1.0 - (((1.0 - panelPtr[i].hinge2) / (double(panelPtr[i].m) - double(closeLL))) * (double(panelPtr[i].m) - double(m)))) - (double(m) / double(panelPtr[i].m))) * panelPtr[i].c2;

					adjust1next = ((1. - (((1. - panelPtr[i].hinge1) / (double(panelPtr[i].m) - double(closeLL))) * (double(panelPtr[i].m) - (double(m) + 1.)))) - ((double(m) + 1.) / double(panelPtr[i].m))) * panelPtr[i].c1;
					adjust2next = ((1. - (((1. - panelPtr[i].hinge2) / (double(panelPtr[i].m) - double(closeLL))) * (double(panelPtr[i].m) - (double(m) + 1.)))) - ((double(m) + 1.) / double(panelPtr[i].m))) * panelPtr[i].c2;

					chord1 = (delchord1 * (double(m) + 1.) + adjust1next) - ((delchord1 * double(m)) + adjust1);
					chord2 = (delchord2 * (double(m) + 1.) + adjust2next) - ((delchord2 * double(m)) + adjust2);
				}
			}
			else { //we do not need to move LL
				chord1 = panelPtr[i].c1 / panelPtr[i].m;
				chord2 = panelPtr[i].c2 / panelPtr[i].m;
			}
			
			//computing left LE locations of current spanwise row of DVEs

			x1[0] = x1LE[0]+(delchord1*m + adjust1)*cos(panelPtr[i].eps1);
			x1[1] = x1LE[1]+(delchord1*m + adjust1)*sin(panelPtr[i].eps1)*sin(nu);
			x1[2] = x1LE[2]-(delchord1*m + adjust1)*sin(panelPtr[i].eps1)*cos(nu);

			//computing right LE locations of current spanwise row of DVEs
			x2[0] = x2LE[0] + (delchord2* m + adjust2) * cos(panelPtr[i].eps2);
			x2[1] = x2LE[1] + (delchord2* m + adjust2) * sin(panelPtr[i].eps2) * sin(nu);
			x2[2] = x2LE[2] - (delchord2* m + adjust2) * sin(panelPtr[i].eps2) * cos(nu);



			delTANphi = (chord2 - chord1) / tempSpan; //define this here, before updating chord with camber

				// Assign espilon of m section as the panel edge eps 
			// Needed for camber or trim
			eps1 = panelPtr[i].eps1;
			eps2 = panelPtr[i].eps2;

			// If there is camber, add it to x1 and x2 and compute new
			//chord and epsilon values
			if(info.flagCAMBER){
				Apply_Camber(panelPtr,x1,x2, x1LE,x2LE,camberPtr, m, i, nu, \
					&epsC1, &epsC2, &chord1, &chord2);
				eps1 += epsC1;
				eps2 += epsC2;
			}

			// If there is a hinge deflection, do it now for this element
			if (panelPtr[i].deflect1 != 0) {

				if ((panelPtr[i].deflect1 != panelPtr[i].deflect2)) {
					
					printf("Hinge deflection must be equal for each side of a panel");
					printf("  Error panel %d	\n", i+1);
					printf("---Exiting program---\n");
					scanf("%c", &answer);
					exit(0);
				}

				DeflectAboutHinge(panelPtr, panelPtr[i].deflect1, x1, x2, x1LE,x2LE, m, i, \
					nu, &epsH1, &epsH2, xH1, xH2);
				eps1 += epsH1; //If there is camber, add the epsilon to that value
				eps2 += epsH2;
			}


			// If there is trim adjust the tail at its hinge according to 
			//the deflection of epsilonHT
			//removed and replaced with above block. BB2020
			/*if(info.trimPITCH){ //removed and replaced with new trimming, BB 2020
				DeflectAboutHinge(panelPtr,epsilonHT, x1, x2, m, i, \
								nu, &epsH1, &epsH2, xH1, xH2);
					eps1 +=epsH1; //If there is camber, add the epsilon to that value
					eps2 +=epsH2;
			}*/
			
			//computing vector along LE of current spanwise row of DVEs
			tempS = 1./panelPtr[i].n;
			xLE[0] = (x2[0]-x1[0])*tempS;
			xLE[1] = (x2[1]-x1[1])*tempS;
			xLE[2] = (x2[2]-x1[2])*tempS;

			//computing dihedral of LE of current spanwise row of DVEs
			tempS	= sqrt(xLE[1]*xLE[1]+xLE[2]*xLE[2]);
			nu2 	= asin(xLE[2]/tempS);
			if(xLE[1]<0) nu2 = Pi-nu2;  //Added to allows for panel to be defined from right to left and dihedral over 90deg 

			//loop over number of spanwise elements 'n'
			for (n=0;n<panelPtr[i].n;n++)
			{
				tempS = (0.5+n);

				//half-chord length at midspan of DVE
				surfacePtr[l].xsi = 0.5 * (chord1 + tempS * (chord2 - chord1) / panelPtr[i].n);
				//temporary incidence angle at half span of DVE
				surfacePtr[l].epsilon = eps1 + tempS * (eps2 - eps1) / panelPtr[i].n;

				ceps  = cos(surfacePtr[l].epsilon); //needed for tempA below
				seps  = sin(surfacePtr[l].epsilon);

				//computing the mid-span location of the DVE LE
				//temporarily saved in xleft
				tempA[0] = x1[0]+xLE[0]*tempS;
				tempA[1] = x1[1]+xLE[1]*tempS;
				tempA[2] = x1[2]+xLE[2]*tempS;

				//vector from LE mid-span point to reference point
				delX[0] =  surfacePtr[l].xsi * ceps;
				delX[1] =  surfacePtr[l].xsi * seps * sin(nu2);
				delX[2] = -surfacePtr[l].xsi * seps * cos(nu2);

				//computing the DVE reference point
				surfacePtr[l].xo[0] = tempA[0] + delX[0];
				surfacePtr[l].xo[1] = tempA[1] + delX[1];
				surfacePtr[l].xo[2] = tempA[2] + delX[2];

				//computing the normal of DVE,
				cross(delX,xLE,tempA);
				tempS = 1/norm2(tempA);
				scalar(tempA,tempS,surfacePtr[l].normal);

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

				//add up panel surface area
				panelPtr[i].surfAREA += surfacePtr[l].S;

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

		info.surfAREA += panelPtr[i].surfAREA;
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

	//#################################################################
	//plots geometry for checking
	//uses parts of code from Main_PlotQuiver20PY.cpp

	//General algorithm
	//	DEFAULTS - 1. user input of timestep and width intervals that are to be plotted
	//	DONE - 2. reading in data of timestep  - DONE
	//	3. writing to plotting file (with extension .py)
	//  4. running the python script in system command 

	double x3[3],x4[3];	//corner points
	char iofile[125];				//input-output-file
	FILE *fp;						//output file

//	3. writing to plotting file (with extension .py)
	sprintf(iofile,"%s","geometryCheck.py");
	//opens input file
	fp = fopen(iofile, "w");

// writing Header
	fprintf(fp,"#Plotting wake results that were generated with ");
	fprintf(fp,"#%s\n",PROGRAM_VERSION);
	//comment
	fprintf(fp,"#Plotting wake of timestep %d\n\n",0);
	
// importing required librarys
	fprintf(fp,"import matplotlib as mpl\n");
	fprintf(fp,"from mpl_toolkits.mplot3d import Axes3D\n");
	fprintf(fp,"import numpy as np\n");
	fprintf(fp,"import matplotlib.pyplot as plt\n");

// setup an "set_axes_equal" function
	fprintf(fp,"\ndef set_axes_equal(ax):\n");
	fprintf(fp,"\tx_limits = ax.get_xlim3d()\n");
	fprintf(fp,"\ty_limits = ax.get_ylim3d()\n");
	fprintf(fp,"\tz_limits = ax.get_zlim3d()\n");	
	fprintf(fp,"\tx_range = abs(x_limits[1] - x_limits[0])\n");
	fprintf(fp,"\tx_middle = np.mean(x_limits)\n");	
	fprintf(fp,"\ty_range = abs(y_limits[1] - y_limits[0])\n");
	fprintf(fp,"\ty_middle = np.mean(y_limits)\n");
	fprintf(fp,"\tz_range = abs(z_limits[1] - z_limits[0])\n");
	fprintf(fp,"\tz_middle = np.mean(z_limits)\n");
	fprintf(fp,"\tplot_radius = 0.5*max([x_range, y_range, z_range])\n");
	fprintf(fp,"\tax.set_xlim3d([x_middle - plot_radius, x_middle + plot_radius])\n");
	fprintf(fp,"\tax.set_ylim3d([y_middle - plot_radius, y_middle + plot_radius])\n");
	fprintf(fp,"\tax.set_zlim3d([z_middle - plot_radius, z_middle + plot_radius])\n\n");

// setup plots

	fprintf(fp,"mpl.rcParams['legend.fontsize'] = 10\n");
	fprintf(fp,"fig = plt.figure()\n");
	fprintf(fp,"ax = fig.gca(projection='3d')\n");
	fprintf(fp,"ax.set_aspect('equal')\n\n");

	// Plot labels 
	fprintf(fp,"ax.set_xlabel('Global X')\n");
	fprintf(fp,"ax.set_ylabel('Global Y')\n");
	fprintf(fp,"ax.set_zlabel('Global Z')\n\n");

//Plotting wing
	for(n=0;n<info.noelement;n++)
	{
//		//computing left-leading edge point in local ref. frame
//		tempA[0] = -surfacePtr[n].xsi\
//				 - surfacePtr[n].eta*tan(surfacePtr[n].phiLE);
//		tempA[1] = -surfacePtr[n].eta;
//		tempA[2] = 0;

//		Star_Glob(tempA,surfacePtr[n].nu,\
			surfacePtr[n].epsilon,surfacePtr[n].psi,tempAA);
//		vsum(tempAA,surfacePtr[n].xo,x1);

		x1[0] = surfacePtr[n].x1[0];
		x1[1] = surfacePtr[n].x1[1];
		x1[2] = surfacePtr[n].x1[2];


		//computing left-trailing edge point in local ref. frame
		tempA[0] = surfacePtr[n].xsi\
				 - surfacePtr[n].eta*tan(surfacePtr[n].phiTE);
		tempA[1] = -surfacePtr[n].eta;
		tempA[2] = 0;

		Star_Glob(tempA,surfacePtr[n].nu,\
			surfacePtr[n].epsilon,surfacePtr[n].psi,tempAA);
		vsum(tempAA,surfacePtr[n].xo,x2);


		//computing right-trailing edge point in local ref. frame
		tempA[0] = surfacePtr[n].xsi\
				 + surfacePtr[n].eta*tan(surfacePtr[n].phiTE);
		tempA[1] = surfacePtr[n].eta;
		tempA[2] = 0;

		Star_Glob(tempA,surfacePtr[n].nu,\
			surfacePtr[n].epsilon,surfacePtr[n].psi,tempAA);
		vsum(tempAA,surfacePtr[n].xo,x3);


		//computing right-leading edge point in local ref. frame
//		tempA[0] = -surfacePtr[n].xsi\
//				 + surfacePtr[n].eta*tan(surfacePtr[n].phiLE);
//		tempA[1] = surfacePtr[n].eta;
//		tempA[2] = 0;

//		Star_Glob(tempA,surfacePtr[n].nu,\
//			surfacePtr[n].epsilon,surfacePtr[n].psi,tempAA);
//		vsum(tempAA,surfacePtr[n].xo,x4);
		x4[0] = surfacePtr[n].x2[0];
		x4[1] = surfacePtr[n].x2[1];
		x4[2] = surfacePtr[n].x2[2];



		//print plot coordinates
		fprintf(fp,"ax.plot("); //first plot command
//#		fprintf(fp,"plot("); //first plot command
		fprintf(fp,"[%lf,%lf,%lf,%lf,%lf],",x1[0],x2[0],x3[0],x4[0],x1[0]); //x
		fprintf(fp,"[%lf,%lf,%lf,%lf,%lf],",x1[1],x2[1],x3[1],x4[1],x1[1]); //y
		fprintf(fp,"[%lf,%lf,%lf,%lf,%lf],",x1[2],x2[2],x3[2],x4[2],x1[2]); //z
		fprintf(fp,"'k')\n");//next line
	}

	// Python lines to show plot and run axis equal
	fprintf(fp,"\nset_axes_equal(ax)\n");
	fprintf(fp,"plt.show()\n");

	fclose(fp);

	// run python script 
	//system("python geometryCheck.py");

	printf("\nTo check geometry type: \"python geometryCheck.py\"\n\n");


	//DONE with geometry check 

	//#################################################################

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
void Move_Wing(const GENERAL info, DVE* surfacePtr,const double circCenter[3], double XCG[3])
{
//moves wing by delta x every time step,
//function updates xo location of surface DVE's
//
//	INPUTS:
//		info			general information
//		surfacePtr		surface DVE's
//		circCenter		center point of circling flight Added by D.F.B. 03-2020

int i;
double delx[3],delX1[3],delX2[3];
double tempA[3];
double rotAngle; // How many radian to rotate points


	if(!info.flagCIRC){
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

			surfacePtr[i].x1[0] -= delx[0];
			surfacePtr[i].x1[1] -= delx[1];
			surfacePtr[i].x1[2] -= delx[2];

			surfacePtr[i].x2[0] -= delx[0];
			surfacePtr[i].x2[1] -= delx[1];
			surfacePtr[i].x2[2] -= delx[2];

		}
		// Updating CG location was moved from PitchMoment - D.F.B. 03-2020
		//move CG
		//newXCG -= local U * delta time
		XCG[0] -= surfacePtr[0].u[0] * info.deltime;
		XCG[1] -= surfacePtr[0].u[1] * info.deltime;
		XCG[2] -= surfacePtr[0].u[2] * info.deltime;

	} else{
		for(i=0;i<info.noelement;i++)
		{	
			// -------Move points for circling flight-------
			// D.F.B. in Braunschweig, Germany, Mar. 2020

			// Angle to rotate points
			rotAngle = info.gradient * info.deltime;

			// Apply rotation matrix about the z-axis
			// Calculate vector from rotation center to control point		
			delx[0] = surfacePtr[i].xo[0] - circCenter[0];
			delx[1] = surfacePtr[i].xo[1] - circCenter[1];

			// Move control points
			surfacePtr[i].xo[0] = delx[0]*cos(rotAngle)-delx[1]*sin(rotAngle)+circCenter[0];
			surfacePtr[i].xo[1] = delx[0]*sin(rotAngle)+delx[1]*cos(rotAngle)+circCenter[1];
			surfacePtr[i].xo[2] -= info.U[2] * info.deltime;

			// Move left leading edge point
			delx[0] = surfacePtr[i].x1[0] - circCenter[0];
			delx[1] = surfacePtr[i].x1[1] - circCenter[1];
			surfacePtr[i].x1[0] = delx[0]*cos(rotAngle)-delx[1]*sin(rotAngle)+circCenter[0];
			surfacePtr[i].x1[1] = delx[0]*sin(rotAngle)+delx[1]*cos(rotAngle)+circCenter[1];
			surfacePtr[i].x1[2] -= info.U[2] * info.deltime;

			// Move right leading edge point
			delx[0] = surfacePtr[i].x2[0] - circCenter[0];
			delx[1] = surfacePtr[i].x2[1] - circCenter[1];
			surfacePtr[i].x2[0] = (delx[0]*cos(rotAngle)-delx[1]*sin(rotAngle))+circCenter[0];
			surfacePtr[i].x2[1] = (delx[0]*sin(rotAngle)+delx[1]*cos(rotAngle))+circCenter[1];
			surfacePtr[i].x2[2] -= info.U[2] * info.deltime;				
		}

		// Recompute DVE param based on new LE and Control pts
		DVE_LEandCPtoParam(info, surfacePtr);

		//Move CG with the circling flight
		delx[0] = XCG[0] - circCenter[0];
		delx[1] = XCG[1] - circCenter[1];
		XCG[0] = (delx[0]*cos(rotAngle)-delx[1]*sin(rotAngle))+circCenter[0];
		XCG[1] = (delx[0]*sin(rotAngle)+delx[1]*cos(rotAngle))+circCenter[1];
		XCG[2] -= info.U[2] * info.deltime;

	}
}
//===================================================================//
		//END FUNCTION Move_Wing
//===================================================================//


//===================================================================//
		//START FUNCTION Apply_Camber 
//===================================================================//
void Apply_Camber(const PANEL* panelPtr, double x1[3], double x2[3], \
	double x1LE[3], double x2LE[3],	double ***camberPtr, int m, int i,\
	double nu, double *epsC1, double *epsC2, double *chord1, double *chord2)
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
	double tempA[3], tempAA[3];
	double tempx1LE[3], tempx2LE[3];	//Local reference panel LE pts
	double leftZ, rightZ;		//Delta Z from LE to TE

	// Bring the input x1 and x2 into the local DVE ref frame
	Glob_Star(x1,nu,panelPtr[i].eps1,0,tempx1);
	Glob_Star(x2,nu,panelPtr[i].eps2,0,tempx2);
	Glob_Star(x1LE, nu, panelPtr[i].eps1, 0, tempx1LE);
	Glob_Star(x2LE, nu, panelPtr[i].eps2, 0, tempx2LE);

	//tempA[0] = tempx1[0] - tempx1LE[0];
	//tempA[1] = tempx1[1] - tempx1LE[1];
	//tempA[2] = tempx1[2] - tempx1LE[2];

	//tempAA[0] = tempx2[0] - tempx2LE[0];
	//tempAA[1] = tempx2[1] - tempx2LE[1];
	//tempAA[2] = tempx2[2] - tempx2LE[2];

	//First consider LE left side
	// Seach for index where (current m)/(total m) is nearest the camber data
	j = 0;
	do{j++;}
	//while(camberPtr[panelPtr[i].airfoil1][j][0]<(double(m)/double(panelPtr[i].m)));
	while (camberPtr[panelPtr[i].airfoil1][j][0] < ((tempx1[0] - tempx1LE[0]) / panelPtr[i].c1));

	//Calculate the z/c by linearly interpolate the using the above define index 
	tempZ1 = camberPtr[panelPtr[i].airfoil1][j-1][1] + (((tempx1[0] - tempx1LE[0]) / panelPtr[i].c1 -camberPtr[panelPtr[i].airfoil1][j-1][0])*\
			(camberPtr[panelPtr[i].airfoil1][j][1]-camberPtr[panelPtr[i].airfoil1][j-1][1])/\
			(camberPtr[panelPtr[i].airfoil1][j][0]-camberPtr[panelPtr[i].airfoil1][j-1][0]));
	//tempZ1 = camberPtr[panelPtr[i].airfoil1][j - 1][1] + ((double(m) / double(panelPtr[i].m) - camberPtr[panelPtr[i].airfoil1][j - 1][0]) * \
		(camberPtr[panelPtr[i].airfoil1][j][1] - camberPtr[panelPtr[i].airfoil1][j - 1][1]) / \
		(camberPtr[panelPtr[i].airfoil1][j][0] - camberPtr[panelPtr[i].airfoil1][j - 1][0]));
	//Scale z/c from camber data to the chordlength of the edge
	tempx1[2] +=(tempZ1*panelPtr[i].c1); 

	//Repeat for TE left side
	j = 0;
	tempA[0] = tempx1[0] + *chord1 - tempx1LE[0];
	if (tempA[0] > panelPtr[i].c1) tempA[0] = panelPtr[i].c1; //because of small errors, this can end up off the trailing edge and the interp will throw an error
	do{j++;}
	while(camberPtr[panelPtr[i].airfoil1][j][0]<((tempA[0]) / panelPtr[i].c1) );
	//while (camberPtr[panelPtr[i].airfoil1][j][0] < double(m+1)/double(panelPtr[i].m));
	tempZTE1= camberPtr[panelPtr[i].airfoil1][j-1][1] + (((((tempA[0])) / panelPtr[i].c1) -camberPtr[panelPtr[i].airfoil1][j-1][0])*\
			(camberPtr[panelPtr[i].airfoil1][j][1]-camberPtr[panelPtr[i].airfoil1][j-1][1])/\
			(camberPtr[panelPtr[i].airfoil1][j][0]-camberPtr[panelPtr[i].airfoil1][j-1][0]));
	//tempZTE1 = camberPtr[panelPtr[i].airfoil1][j - 1][1] + ((double(m + 1) / double(panelPtr[i].m) - camberPtr[panelPtr[i].airfoil1][j - 1][0]) * \
		(camberPtr[panelPtr[i].airfoil1][j][1] - camberPtr[panelPtr[i].airfoil1][j - 1][1]) / \
		(camberPtr[panelPtr[i].airfoil1][j][0] - camberPtr[panelPtr[i].airfoil1][j - 1][0]));
	//Repeat for LE right side
	j = 0;
	do{j++;}
	//while(camberPtr[panelPtr[i].airfoil2][j][0]<(double(m)/double(panelPtr[i].m)));
	while (camberPtr[panelPtr[i].airfoil2][j][0] < ((tempx2[0] - tempx2LE[0]) / panelPtr[i].c2));

	//tempZ2 = camberPtr[panelPtr[i].airfoil2][j - 1][1] + ((double(m)/double(panelPtr[i].m) - camberPtr[panelPtr[i].airfoil2][j - 1][0]) * \
		(camberPtr[panelPtr[i].airfoil2][j][1] - camberPtr[panelPtr[i].airfoil2][j - 1][1]) / \
		(camberPtr[panelPtr[i].airfoil2][j][0] - camberPtr[panelPtr[i].airfoil2][j - 1][0]));
	//tempx2[2] += (tempZ2 * panelPtr[i].c2);

	tempZ2 = camberPtr[panelPtr[i].airfoil2][j-1][1] + (((tempx2[0] - tempx2LE[0]) / panelPtr[i].c2 -camberPtr[panelPtr[i].airfoil2][j-1][0])*\
			(camberPtr[panelPtr[i].airfoil2][j][1]-camberPtr[panelPtr[i].airfoil2][j-1][1])/\
			(camberPtr[panelPtr[i].airfoil2][j][0]-camberPtr[panelPtr[i].airfoil2][j-1][0]));
	tempx2[2] +=(tempZ2*panelPtr[i].c2); 

	//Repeat for TE right side
	j = 0;
	tempAA[0] = tempx2[0] + *chord2 - tempx2LE[0];
	if (tempAA[0] > panelPtr[i].c2) tempAA[0] = panelPtr[i].c2; //because of small errors, this can end up off the trailing edge and the interp will throw an error
	do{j++;}
	while(camberPtr[panelPtr[i].airfoil2][j][0]< ((tempAA[0]) / panelPtr[i].c2) );
	//while (camberPtr[panelPtr[i].airfoil2][j][0] < (double(m + 1) / double(panelPtr[i].m)));
	tempZTE2 = camberPtr[panelPtr[i].airfoil2][j-1][1] + (((((tempAA[0]) ) / panelPtr[i].c2) -camberPtr[panelPtr[i].airfoil2][j-1][0])*\
		(camberPtr[panelPtr[i].airfoil2][j][1]-camberPtr[panelPtr[i].airfoil2][j-1][1])/\
		(camberPtr[panelPtr[i].airfoil2][j][0]-camberPtr[panelPtr[i].airfoil2][j-1][0]));
	//tempZTE2 = camberPtr[panelPtr[i].airfoil2][j - 1][1] + ((double(m + 1) / double(panelPtr[i].m) - camberPtr[panelPtr[i].airfoil2][j - 1][0]) * \
		(camberPtr[panelPtr[i].airfoil2][j][1] - camberPtr[panelPtr[i].airfoil2][j - 1][1]) / \
		(camberPtr[panelPtr[i].airfoil2][j][0] - camberPtr[panelPtr[i].airfoil2][j - 1][0]));

	//Calculate delta z from LE to TE on left edge and on right edge
	leftZ = (tempZ1-tempZTE1)*panelPtr[i].c1;
	rightZ =(tempZ2-tempZTE2)*panelPtr[i].c2;

	// Calculate mid-span eps
	//*epsC1 = atan(leftZ / (panelPtr[i].c1 / panelPtr[i].m));
	//*epsC2 = atan(rightZ / (panelPtr[i].c2 / panelPtr[i].m));
	// Calculate mid-span eps
	*epsC1 = atan(leftZ / *chord1);
	*epsC2 = atan(rightZ / *chord2);


	//New section chord on left and right edges
	//*chord1 = sqrt(leftZ*leftZ+(panelPtr[i].c1/panelPtr[i].m)*(panelPtr[i].c1/panelPtr[i].m));
	//*chord2 = sqrt(rightZ*rightZ+(panelPtr[i].c2/panelPtr[i].m)*(panelPtr[i].c2/panelPtr[i].m));
	*chord1 = sqrt(leftZ * leftZ + (*chord1) * (*chord1));
	*chord2 = sqrt(rightZ * rightZ + (*chord2) * (*chord2));

	// Convert new LE points to the global reference frame 
	Star_Glob(tempx1,nu,panelPtr[i].eps1,0,x1);
	Star_Glob(tempx2,nu,panelPtr[i].eps2,0,x2);

}
//===================================================================//
		//END FUNCTION Apply_Camber
//===================================================================//


//===================================================================//
		//START FUNCTION DeflectAboutHinge 
//===================================================================//
void DeflectAboutHinge(const PANEL* panelPtr, const double deflection, \
	double x1[3], double x2[3], double x1LE[3], double x2LE[3],  int m, int i, double nu, \
	double *epsH1, double *epsH2, double xH1[3], double xH2[3])
{
	// The function DeflectAboutHinge deflects the DVEs aft of a hinge by some
	//		deflection angle. This function will work for a cambered wing.
	//
	// Note: There must be a LE at the DVE hinge or the program will exit with warning
	//
	// Function inputs:
	//		PANEL* panelPtr -Panel structure for the panelPtr
	//		deflection 		-Deflection angle (rad)
	//		x1,x2 			-x,y,z position of the left (1) and right (2) panel LE pts
	//		m 				-chorwise row of interest
	//		i 				-panelPtr of interest
	//		nu 				-dihedral angle of panel
	//		epsH1, epsH2 	-(see output) 
	//		xH1,xH2			-(see output) 
	//
	// Function outputs:
	//		x1,x2 			-updated left and right panel LE points
	//		epsH1, epsH2 	-updated left and right edge epsilon angles with deflection
	//		xH1,xH2			-x,y,z position of hinge for left and right panel edge
	//							-this is created one m is add the LE point
	//
	// NOTE (1) indicated left and (2) indicated right for the above variables

	// D.F.B. in Braunschweig, Germany, Feb. 2020

	int j;						//Generic counter
	double check1,check2;		//Checks to make sure that hinge is a DVE LE
	double tempx1[3],tempx2[3];	//Local reference frame LE pts
	double vecX, vecZ;				//x and z of a vector from hinge to LE pt. (in local)
	double tempxH1[3],tempxH2[3];	//Hinge location in local reference frame
	double tempx1LE[3], tempx2LE[3]; //Leading edge of panel in local
	double tempA[3], tempAA[3];
	//removed BB 2020, since we move lifting lines outside and chordwise dist. is no
	//longer uniform, this check will not work. Therefore, we assume that there will be a
	//lifting line at the hinge line by this point.

	// Check if hinge point is at the LE of an element
	//check1 = (panelPtr[i].m*panelPtr[i].hinge1 - round(panelPtr[i].m*panelPtr[i].hinge1));
	//check2 = (panelPtr[i].m*panelPtr[i].hinge2 - round(panelPtr[i].m*panelPtr[i].hinge2));

	// Give a warning message and exit program if left panel hinge isnt at a DVE LE
	//if(fabs(check1)>DBL_EPS){
	//	printf("Hinge line of the left edge of panel %d  is not at a DVE LE.\n",i+1);
	//	printf("Adjust the number of chordwise elements or hinge location.\n");
	//	printf("---Exiting program---\n");
	//	exit(0);
	//}
	// Give a warning message and exit program if right panel hinge isnt at a DVE LE
	//if(fabs(check2)>DBL_EPS){
	//	printf("Hinge line of the right edge of panel %d not at a DVE LE.\n",i+1);
	//	printf("Adjust the number of chordwise elements or hinge location.\n");
	//	printf("---Exiting program---\n");
	//	exit(0);
	//}

	// ================== Deflecting about hinge begins here ==================
	// First deflect left side of panel
	Glob_Star(x1, nu, panelPtr[i].eps1, 0, tempx1);
	Glob_Star(x1LE, nu, panelPtr[i].eps1, 0, tempx1LE);

	tempA[0] = tempx1[0] - tempx1LE[0];
	tempA[2] = tempx1[2] - tempx1LE[2];

	if(tempA[0] /panelPtr[i].c1 < panelPtr[i].hinge1 && (panelPtr[i].hinge1 - (tempA[0] / panelPtr[i].c1)) > (1000*DBL_EPS)) {
		// left point is in front of hinge. If it is, set epsH1 to 0
		*epsH1 = 0;
		
	}else if(((tempA[0] / panelPtr[i].c1) -panelPtr[i].hinge1) < (1000.*DBL_EPS)){
		// left point at hinge
		*epsH1 = deflection;
		
		//save this point as hinge point for future calcs
		for(j=0;j<3;j++){xH1[j]=x1[j];}
	}else{
		// assume element is downstream of hinge		
		// set epsH1 to deflection
		// angle and move the LE point accordingly
		*epsH1 = deflection;

		// Put DVE in local reference frame
		Glob_Star(x1,nu,panelPtr[i].eps1,0,tempx1);
		// Put the hinge location in the local reference frame
		Glob_Star(xH1,nu,panelPtr[i].eps1,0,tempxH1);

		// Create a vector from the hinge point to the LE point of interest
		// Because in local reference frame, only consider the x and z components
		vecX = tempx1[0]-tempxH1[0];
		vecZ = tempx1[2]-tempxH1[2];

		// Multiply the vectors from hinge to LE pt with a rotation matrix
		// Note that down deflection is considered positive thus -ve deflection is used
		tempx1[0] = vecX*cos(-deflection)-vecZ*sin(-deflection)+tempxH1[0];
		tempx1[2] = vecX*sin(-deflection)+vecZ*cos(-deflection)+tempxH1[2];

		// Put new LE point into global reference frame
		Star_Glob(tempx1,nu,panelPtr[i].eps1,0,x1);
	}

	// Repeat everything for right edge (see above for detailed comments)
	Glob_Star(x2, nu, panelPtr[i].eps2, 0, tempx2);
	Glob_Star(x2LE, nu, panelPtr[i].eps2, 0, tempx2LE);

	tempAA[0] = tempx2[0] - tempx2LE[0];
	tempAA[2] = tempx2[2] - tempx2LE[2];

	if(tempAA[0] / panelPtr[i].c2 < panelPtr[i].hinge2 && (panelPtr[i].hinge2-(tempAA[0] / panelPtr[i].c2)) > (1000.*DBL_EPS)){
		*epsH2 = 0;	
	}else if(((tempAA[0] / panelPtr[i].c2 )-panelPtr[i].hinge2)<(1000.*DBL_EPS)){
		*epsH2 = deflection;
		for(j=0;j<3;j++){xH2[j]=x2[j];}
	}else{
		
		*epsH2 = deflection;

		Glob_Star(x2,nu,panelPtr[i].eps2,0,tempx2);
		Glob_Star(xH2,nu,panelPtr[i].eps2,0,tempxH2);

		vecX = tempx2[0]-tempxH2[0];
		vecZ = tempx2[2]-tempxH2[2];
		tempx2[0] = vecX*cos(-deflection)-vecZ*sin(-deflection)+tempxH2[0];
		tempx2[2] = vecX*sin(-deflection)+vecZ*cos(-deflection)+tempxH2[2];
		
		Star_Glob(tempx2,nu,panelPtr[i].eps2,0,x2);
	}
}
//===================================================================//
		//END FUNCTION DeflectAboutHinge
//===================================================================//


//===================================================================//
		//START FUNCTION Circling_UINF
//===================================================================//
void Circling_UINF(GENERAL info, DVE* surfacePtr,const double circCenter[3])
{
	// Calculates the velocity vector at each DVE given the cross product of
	//	Omega x r. Where omega is the velocity gradient given in the input
	//	file and r is a vector from the center of rotation to any given
	//	dve control point. The z velocity is given by info.U[2].
	//
	// Function inputs:
	//		GENERAL info 	- Uses info.noelements and info.U[2]
	//		DVE* surfacePtr	- DVE info
	//		circCenter		- Center point of circling flight
	//
	// Function outputs:
	//		Updated surfacePtr with the following variables changed:
	//			u1,u2,u - Velocity at left edge, right edge and control point
	//			uTE[3][3] - Velocities at the TE. uTE[0] halfspan, uTE[1] -80%, uTE[2] +80%
	//						 This is required for induced drag calculation
	//
	// D.F.B. in Braunschweig, Germany, Mar. 2020

	int i; 					//Generic counter
	double tempA[3];		//Generic temp vector
	double r[3];			//Distance from center of rotation to DVE control pt
	double omega[3];		//Rotational rate [0 0 gradient]
	double X[3][3];			//Trailing edge points of DVE X[0] halfspan, X[1] -80%, X[2] +80%
	double eta8;			// 80% of half span


	// This function first calls calculates the velocities at the control point and LE pts
	// Than calculated the velocities along the TE

	// Need to declare the Edge_Point function in order to use it
	//void CreateQuiverFile(const double[3], const double[3],const int);
	void Edge_Point(const double [3],const double,const double,const double,\
				const double,const double,const double,double [3]);
	// Iterate through number of elements
	for(i=0;i<info.noelement;i++)
	{

		// Create Omega vector using only the component of the rotation in the 
		//XY plane. This assumes that Vk = Vinf*cos(alpha). 
		//This also works with horizontal flight because we force alpha = 0
		//in wing_geometry line 299.
		omega[0] = 0;
		omega[1] = 0;
        omega[2] = -info.gradient; 
		

		// ********************* Calculate u,u1,u2 *********************
		// Calculate r, the vector from the center of rotation to the DVE
		r[0] = surfacePtr[i].xo[0]-circCenter[0];
		r[1] = surfacePtr[i].xo[1]-circCenter[1];
		r[2] = surfacePtr[i].xo[2]-circCenter[2];

		// Calculate the local velocities (Omega cross r)
		cross(omega,r,surfacePtr[i].u);
		surfacePtr[i].u[2] = info.U[2]; // Apply Z velocity based on input

		// Repeat for left edge
		r[0] = surfacePtr[i].x1[0]-circCenter[0];
		r[1] = surfacePtr[i].x1[1]-circCenter[1];
		r[2] = surfacePtr[i].x1[2]-circCenter[2];
		cross(omega,r,surfacePtr[i].u1);
		surfacePtr[i].u1[2] = info.U[2];

		// Repeat for right edge
		r[0] = surfacePtr[i].x2[0]-circCenter[0];
		r[1] = surfacePtr[i].x2[1]-circCenter[1];
		r[2] = surfacePtr[i].x2[2]-circCenter[2];
		cross(omega,r,surfacePtr[i].u2);
		surfacePtr[i].u2[2] = info.U[2];

		// ********************* Calculate uTE  *********************
		// Compuate velocity at TE for induced drag calcs
		// First calculate the location of the TE points: X[3][3]
		// Than use the same method as above to calculate the velocities
		//
		// The following code was based off code found in Induced_DVE_Drag
		
		//the left and right points are 20% of half span away from edge
		//in order to stay away from the singularity along the edge of
		//the DVE
		eta8  =	surfacePtr[i].eta*0.8; //0.8 as done for lift computation,

		//X1:
		Edge_Point(surfacePtr[i].xo,surfacePtr[i].nu,\
				   surfacePtr[i].epsilon,surfacePtr[i].psi,\
				   surfacePtr[i].phiTE,-eta8,surfacePtr[i].xsi,X[1]);
				   						//Subroutine in wake_geometry.cpp
		//X2:
		Edge_Point(surfacePtr[i].xo,surfacePtr[i].nu,\
				   surfacePtr[i].epsilon,surfacePtr[i].psi,\
				   surfacePtr[i].phiTE,eta8,surfacePtr[i].xsi,X[2]);
				   						//Subroutine in wake_geometry.cpp
		//X0 = (X1+X2)/2
		vsum(X[1],X[2],tempA);		scalar(tempA,0.5,X[0]);


		// Calculate the TE velocities
		r[0] = X[0][0]-circCenter[0];
		r[1] = X[0][1]-circCenter[1];
		r[2] = X[0][2]-circCenter[2];
		cross(omega,r,surfacePtr[i].uTE[0]);
		surfacePtr[i].uTE[0][2] = info.U[2];

		r[0] = X[1][0]-circCenter[0];
		r[1] = X[1][1]-circCenter[1];
		r[2] = X[1][2]-circCenter[2];
		cross(omega,r,surfacePtr[i].uTE[1]);
		surfacePtr[i].uTE[1][2] = info.U[2];

		r[0] = X[2][0]-circCenter[0];
		r[1] = X[2][1]-circCenter[1];
		r[2] = X[2][2]-circCenter[2];
		cross(omega,r,surfacePtr[i].uTE[2]);
		surfacePtr[i].uTE[2][2] = info.U[2];
		
		/* Code ussed for plotting quiver
		if(i==0){CreateQuiverFile(X[0], surfacePtr[i].uTE[0],0);}
		else{CreateQuiverFile(X[0], surfacePtr[i].uTE[0],1);}
		CreateQuiverFile(X[1], surfacePtr[i].uTE[1],1);
		CreateQuiverFile(X[1], surfacePtr[i].uTE[1],1);
		CreateQuiverFile(X[2], surfacePtr[i].uTE[2],1);
		CreateQuiverFile(X[2], surfacePtr[i].uTE[2],1);*/

	}
}
//===================================================================//
		//END FUNCTION Circling_UINF
//===================================================================//



//===================================================================//
		//START FUNCTION DVE_LEandCPtoParam
//===================================================================//
void DVE_LEandCPtoParam(const GENERAL info, DVE* surfacePtr)
{
	// Calculated the DVE parameters based on the LE pts and control pts.
	// 	This function will iterate through info.noelements to update the
	//	param for all DVEs.
	//
	// Function inputs:
	//		GENERAL info 	- Uses info.noelements 
	//		DVE* surfacePtr	- DVE pointers. Specifically uses:
	//							surfacePtr.xo - control point location
	//							surfacePtr.x1 - left LE point location
	//							surfacePtr.x2 - right LE point location
	//
	// Function outputs:
	//		surfacePtr with the following variables updated:
	//			.normal, .nu, .epsilon, .psi, .phiLE,phiTE,phi0,.eta
	//
	//	This function is based on code found the function: Move_Flex_Wing
	//
	// D.F.B. in Braunschweig, Germany, Mar. 2020
		

int i;							//Generic counter
double tempA[3];				// temp vector
double delX1[3];				// vector from LE midpoint to reference pt
double delX2[3];				// Vector along DVE LE
double delXSI1[3],delXSI2[3];	//delX1 and delX2 in local DVE coordinates
double delPhi;					//change in sweep.


for(i=0;i<info.noelement;i++)
	{
			// This section of code was taken from move_flex_wing
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
	}
}

//===================================================================//
		//END FUNCTION DVE_LEandCPtoParam
//===================================================================//
