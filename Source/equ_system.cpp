
//declerations of subroutines used in Vorticity_Distribution
void FlexWingVortDist(const GENERAL,const PANEL *,\
					DVE *,DVE **,double **,double *,int);
void DVE_Vorticity_Distribution(const GENERAL,const PANEL *,DVE *,\
								DVE **,const int,\
								double **,double *, const int *);
void DVE_KinCond(DVE *,const GENERAL,const PANEL *,double **);
void Vorticity_Distribution(const GENERAL, const PANEL *, BOUND_VORTEX *,\
							double **, double *);
void Trailing_Edge_Vorticity(const GENERAL, const PANEL *,BOUND_VORTEX *,\
							 BOUND_VORTEX *);

//===================================================================//
		//START FlexWingVortDist
//===================================================================//
void FlexWingVortDist(const GENERAL info,const PANEL *panelPtr,\
					DVE *surfacePtr,DVE **wakePtr,double **D,\
					double *R,int timestep)
{
//Solves equation system Dx=R using a Gaussian algorithm
//input
//info			general information
//panelPtr		panel information
//surfacePtr	surface DVEs
//walePtr		wake DVEs
//D				D matrix 3nx3n, n = no. of surface DVEs
//R				right hand side vector
//timestep		current timestep
//
//
//output
//Surface DVEs vorticity coefficients A, B, and C
void DVE_Resultant(const GENERAL,const PANEL *,const DVE *,DVE **,\
				   const int,double *);

	double *x;
	int m,n;

	//computes resultant vector.
	DVE_Resultant(info,panelPtr,surfacePtr,wakePtr,timestep,R);

	ALLOC1D(&x,(info.Dsize));		//frees memory allocated for A

	GaussSolve(D,R,info.Dsize,x);	//subroutine in gauss.cpp

	//assigns circulation coefficients A, B, and C
	for(n=0;n<info.noelement;n++)
	{
		m=n*3;
		surfacePtr[n].A=x[m];
		surfacePtr[n].B=x[m+1];
		surfacePtr[n].C=x[m+2];
//#printf("A= %lf\tB= %lf\tC= %lf\n",x[m],x[m+1],x[m+2]);//##
	}
	FREE1D(&x,(info.Dsize));		//frees memory allocated for A
}
//===================================================================//
		//END FlexWingVortDist
//===================================================================//

//===================================================================//
		//START function DVE_Vorticity_Distribution
//===================================================================//
//In this subroutine the remaining elements of the resultant vector  are
//assembled and the equation system is solved.  Matrix D is already
//contains the lower and upper surfaces required for an abbreviated Gauss
//algorith.  The returned values are the vorticity coefficients for each
//surface DVE, A, B, and C.
//The equation system has the following form:
//				D A = R
//			D  	lower and upper matrix with (3n)^2 elements,
//				where n is the number of elementary wings
//			A   3n vector with the A, B, and C coefficients of the
//				vorticity distribution across each elementary wing.
//			R	pivoted 3n resultant vector.
//				R and D are pivoted and consist of three major parts:
//					Part 1: Boundary cond. between surface DVE's
//							of the same panel.  Magnitude and slope
//							of the vorticity is preserved across the
//							boundaries
//					Part 2: Boundary cond. between surface DVE's of
//							neighboring panels (including free ends).
//							possible conditions are preservation of
//							vorticity (or 0 at free end), slope, or
//							curvature.
//					Part 3: Kinematic condition at the ctrl. points
//
// input
//	info 			general information
//	surfacePtr		pointer to information on surface DVE's
//	wakePtr			pointer to information on wake DVE's
//	timestep		current time step
//	D				lower/upper matrix
//	R				resultant vector
//	pivot			pivoting array
void DVE_Vorticity_Distribution(const GENERAL info,const PANEL *panelPtr,\
								DVE *surfacePtr,DVE **wakePtr,\
								const int timestep,double **D,double *R,\
								const int *pivot)
{
void DVE_Resultant(const GENERAL,const PANEL *,const DVE *,DVE **,\
				   const int,double *);
int n,m;			//counter

//next two lines removed 2-11-06 G.B. don't remember why they were needed
//int element=0;		//surface DVE index
//for(n=0;n<info.noelement;n++)

	//computes resultant vector.
	DVE_Resultant(info,panelPtr,surfacePtr,wakePtr,timestep,R);

	//Solves L*(U*x) = R equation system, lower/upper matrix contained in D
	//x has vorticity coefficients A, B, and C of each surface DVE
	//The result x is returned in array R
	LU_Solver(D,info.Dsize,pivot,R);
						//subroutine in gauss.cpp

	//assigns circulation coefficients A, B, and C to surface DVE's
	for(n=0;n<info.noelement;n++)
	{
		m=n*3;
		surfacePtr[n].A=R[m];
		surfacePtr[n].B=R[m+1];
		surfacePtr[n].C=R[m+2];
//#printf("A= %lf\tB= %lf\tC= %lf\n",R[m],R[m+1],R[m+2]);//#
	}

//##########################################
//FILE *fp;														//#
//fp = fopen("output\\test.txt", "a");									//#
//	for(n=0;n<info.noelement;n++)
//		fprintf(fp,"time: %d A %lf  B %lf C %lf \n",timestep,surfacePtr[n].A,surfacePtr[n].B,surfacePtr[n].C);		//#
//fclose (fp);														//#
//###########################################################################

//!!!! REMOVED, G.B. MAY 6, 04
	//assigns circulation coefficients A,B, and C to first wake DVE's
	//that are located just aft of trailing edge. The coefficients are
	//the same as the ones of the most aft surface DVE's
/*	for(k=0;k<info.nopanel;k++)		//loop over panels
	{
		//forwards surface-DVE index to trailing egde elements
		element += panelPtr[k].n*(info.m-1);

		//loop over number of spanwise elements of current panel
		for(n=0;n<panelPtr[k].n;n++)
		{
		m=element*3;
		wakePtr[timestep][span].A = R[m];
		wakePtr[timestep][span].B = R[m+1];
		wakePtr[timestep][span].C = R[m+2];

		element ++;		//advance surface-DVE index
		span ++;		//advance wake-DVE index
		}//end loop over spanwise elements of current panel
	}//end loop over panels  */
}
//===================================================================//
		//END function DVE_Vorticity_Distribution
//===================================================================//

//===================================================================//
		//START assembly of DVE resultant vector
//===================================================================//
//assembles remaining parts of resultant vector that are due to the
//velocities at each surface element that are the results of the
//free stream and the vorticity in the wake. The component of that velocity
//that is normal to the surface DVE form the non-zero values of the
//resultant vector.  They are part of the kinematic condition, i.e., no
//flow through the elementary-wing surface at the control point
void DVE_Resultant(const GENERAL info,const PANEL *panelPtr,\
				   const DVE *surfacePtr,DVE **wakePtr,const int timestep,\
				   double *R)
{
  int element;		//loop counter over surface DVE's
  int panel,m,i;	//loops over panels, chord lines, span wise elements
  int imax;			//max. no. of spanwise elements that are not at edge
  int n2=2*info.noelement;//2*number of elements
  double w_wake[3];	//velocity induced by wake in control point
  double w_extern[3];	//free stream and wake induced velocities in ith DVE

  element=0; // initializing element index counter

  for(panel=0;panel<info.nopanel;panel++)	//loop over panels
  {
	imax = panelPtr[panel].n-1;

	//loop over chordwise lift. lines
	for(m=0;m<info.m;m++)
	{
	  //checking if left edge is a free tip
	  if(panelPtr[panel].BC1==110)
	  {	//kinematic condition is not being satisfied at the tips, but
	  	//gamma'= 0, see boundary conditions
			R[element] 					= 0;
			R[element+info.noelement]	= 0;
			R[element+n2]				= 0;

	 		//increase DVE index
	 		element ++;
	  }
	  else
	  {	//compute the velocity induced at suface DVE by wake and free stream
		  if(timestep<0)
		  {
			  w_extern[0]=surfacePtr[element].u[0];
			  w_extern[1]=surfacePtr[element].u[1];
			  w_extern[2]=surfacePtr[element].u[2];
		  }
		  else
		  {
			  Wake_DVE_Vel_Induction\
			  (info,surfacePtr[element].xo,wakePtr,timestep,w_wake);
								//Subroutine in induced_velocity

			  //add wake induced vel. and free stream
			  vsum(w_wake,surfacePtr[element].u,w_extern);
		  }
		  //upper two-thirds are zero
		  R[element] 			  	= 0;
		  R[element+info.noelement] = 0;
		  //computing the external velocity normal component
		  R[element+n2] = 4*Pi*dot(w_extern,surfacePtr[element].normal);

		  //increase DVE index
		  element ++;
	  }//done with left side edge

	  for(i=1;i<imax;i++)	//loop over spanwise elements-1
	  {
			//compute the velocity of ith suface DVE control point that
			//is induced by wake with
			if(timestep<0)
			{
				w_extern[0]=surfacePtr[element].u[0];
				w_extern[1]=surfacePtr[element].u[1];
				w_extern[2]=surfacePtr[element].u[2];
			}
			else
			{
				Wake_DVE_Vel_Induction\
				(info,surfacePtr[element].xo,wakePtr,timestep,w_wake);
									//Subroutine in induced_velocity

				//add wake induced vel. and free stream
				vsum(w_wake,surfacePtr[element].u,w_extern);
			}
			//upper two-thirds are zero
			R[element] 				  = 0;
			R[element+info.noelement] = 0;
			//computing the external velocity normal component
			R[element+n2] = 4*Pi*dot(w_extern,surfacePtr[element].normal);

	    //increase DVE index to next one
	    element ++;
	  }//done with loop over spanwise elements -1

	  if(panelPtr[panel].n > 1)	//panel has more than one spanwise element
		{
			//checking if right edge is a free tip
	  	if(panelPtr[panel].BC2==110)
	  	{	//kinematic condition is not being satisfied at the tips, but
	  		//gamma'= 0, see boundary conditions
				R[element] 					= 0;
				R[element+info.noelement]	= 0;
				R[element+n2]				= 0;

		 		//increase DVE index
		 		element ++;
		 	}
	  	else
	  	{	//compute the velocity ind. at suface DVE by wake and free stream
				if(timestep<0)
				{
					w_extern[0]=surfacePtr[element].u[0];
					w_extern[1]=surfacePtr[element].u[1];
					w_extern[2]=surfacePtr[element].u[2];
				}
				else
				{
					Wake_DVE_Vel_Induction\
					(info,surfacePtr[element].xo,wakePtr,timestep,w_wake);
										//Subroutine in induced_velocity

					//add wake induced vel. and free stream
					vsum(w_wake,surfacePtr[element].u,w_extern);
				}
				//upper two-thirds are zero
				R[element] 				 = 0;
				R[element+info.noelement]= 0;
				//computing the external velocity normal component
				R[element+n2] = 4*Pi*dot(w_extern,surfacePtr[element].normal);
				//used to be dot(w_extern,surfacePtr[i].normal),
				//don't know why,but sometimes doesnt work with i.
				//G.B. 2-11-06

		 		//increase DVE index
		 		element ++;
		 	}//done with right side edge
		}
  }//done loop over m, next chord location
  }//done panel, next panel
//##########################################
//#control output of R
//FILE *fp;														//#
//fp = fopen("test.txt", "a");									//#
//fprintf(fp,"R: timestep %d ",timestep);				//#
//for(i=0;i<info.Dsize;i++)	//#
//	fprintf(fp,"%lf  ",R[i]);		//#
//fprintf(fp,"\n");				//#
//fclose (fp);														//#
//###########################################################################

}
//===================================================================//
		//END assembly of DVE resultant vector
//===================================================================//

/*/===================================================================//
!!!!OLD!!!!		//START assembly of DVE resultant vector !!!!OLD!!!!
//===================================================================//
//assembles remaining parts of resultant vector that are due to the
//velocities at each surface element that are the results of the
//free stream and the vorticity in the wake. The component of that velocity
//that is normal to the surface DVE form the non-zero values of the
//resultant vector.  They are part of the kinematic condition, i.e., no
//flow through the elementary-wing surface at the control point
void DVE_Resultant(const GENERAL info,const DVE *surfacePtr,\
				   DVE **wakePtr,const int timestep,\
				   double *R)
{
	int i;				//loop counter over surface DVE's
	int n2=2*info.noelement;//2*number of elements
	double w_wake[3];	//velocity induced by wake in control point
	double w_extern[3];	//free stream and wake induced velocities in ith DVE

	//loop over surface DVE's
	for(i=0;i<info.noelement;i++)
	{
		//compute the velocity of ith suface DVE control point that
		//is induced by wake with
		if(timestep<0)
		{
			w_extern[0]=surfacePtr[i].u[0];
			w_extern[1]=surfacePtr[i].u[1];
			w_extern[2]=surfacePtr[i].u[2];
		}
		else
		{
			Wake_DVE_Vel_Induction\
			(info,surfacePtr[i].xo,wakePtr,timestep,w_wake);
								//Subroutine in induced_velocity

			//add wake induced vel. and free stream
			vsum(w_wake,surfacePtr[i].u,w_extern);
		}
		//upper two-thirds are zero
		R[i] 				= 0;
		R[i+info.noelement] = 0;
		//computing the external velocity normal component
		R[i+n2] = 4*Pi*dot(w_extern,surfacePtr[i].normal);
	}
//#control output of R
//#for(i=0;i<info.Dsize;i++)	//#
//#printf("%lf  ",R[i]);		//#
//#printf("\n");				//#
}*/
//===================================================================//
//	!!!!OLD!!!!//END assembly of DVE resultant vector !!!!OLD!!!!
//===================================================================//

//===================================================================//
		//START function Vorticity_Distribution
//===================================================================//

	//Assembly of equation system that defines vorticity distribution
	//across each elementary wing.  Equation system has the form :
	//				D A = R
	//			D  	matrix with (3n)^2 elements, where n is the number
	//				of elementary wings
	//			A   3n vector with the A, B, and C coefficients of the
	//				vorticity distribution across each elementary wing.
	//			R	3n resultant vector, equal zero except for Part 3.
	//				R and D consist of three major parts:
	//					Part 1: Boundary cond. between elem. wings due
	//							of the same panel.  Magnitude and slope
	//							of the vorticity is preserved across the
	//							boundaries
	//					Part 2: Boundary cond. between elem. wings of
	//							neighboring panels (including free ends).
	//							possible conditions are preservation of
	//							vorticity (or 0 at free end), slope, or
	//							curvature.
	//					Part 3: Kinematic condition in the ctrl. point.
	//
	//The function also solves the equation system for the vorticity
	//coefficients, A, B, and C, of each elementary wing by applying a
	//Gaussian elimination.
void Vorticity_Distribution(const GENERAL info, const PANEL *panelPtr, \
							BOUND_VORTEX *elementPtr, double **D, double *R)
{
//assembles boundary conditions between panels
void BoundaryCond(const BOUND_VORTEX *,const PANEL *,const GENERAL,double **);
void Resultant(const BOUND_VORTEX *, const GENERAL, double *);
void KinematicCond(const BOUND_VORTEX *,const GENERAL,const PANEL *,double **);

int n, m;					//counter
int Dsize=info.Dsize;		//size of matrix D
double *x;					//temporary vector with Vorticity coefficients A,B,C

	//assmebly of resultant vector R
	Resultant(elementPtr,info,R);
						//subroutine in equ_system.cpp

	//assembly of first part of D, boundary cond. of
	//elementary wings within each panel
	BoundaryCond(elementPtr,panelPtr,info,D);
						//subroutine in equ_system.cpp

	//assembles second part of D that deals with
	//the kinematic condition of the elementary wing
	KinematicCond(elementPtr,info,panelPtr,D);
						//subroutine in equ_system.cpp
	//Solves D x = R equation system with Gaussian elimination
	ALLOC1D(&x,(Dsize)); //allocates memory for x

 /*  ###########################################################
 //save D matrix and resultant vector R in file D_matrix.txt
 FILE *fp;
 fp = fopen("D_matrix.txt", "w");
 //writes header line
 fprintf(fp, "\t");
 for(m=0; m<Dsize; m++)
 fprintf(fp, "%ld\t",m);
 fprintf(fp, "\t\tR");
 for(n=0; n<Dsize; n++)
 {
	 //row number
 	fprintf(fp, "\n%d\t",n);
 	//n-th row of D
 	for(m=0; m<Dsize; m++)
 	fprintf(fp, "%lf\t",D[n][m]);
 	//n-th element of R
	fprintf(fp, "\t\t%lf",R[n]);
 }
 fclose(fp);
 //###########################################################//*/

	GaussSolve(D,R,Dsize,x);
						//subroutine in gauss.cpp

	//assigns circulation coefficients A, B, and C
	for(n=0;n<info.noelement;n++)
	{
		m=n*3;
		elementPtr[n].A=x[m];
		elementPtr[n].B=x[m+1];
		elementPtr[n].C=x[m+2];
//#printf("A= %lf\tB= %lf\tC= %lf\n",x[m],x[m+1],x[m+2]);//##
	}
	FREE1D(&x,(Dsize));		//frees memory allocated for A
}
//===================================================================//
		//END Vorticity_Distribution
//===================================================================//

//===================================================================//
		//START assembly of resultant vector
//===================================================================//
//assembles resultant vector of 1 x 3n (with n being the number of
//elementary wings). The first 2n elements of the vector are zero.
//The remaining n elements are part of the kinematic condition, i.e., no
//flow through the elementary-wing surface at the control point
void Resultant(const BOUND_VORTEX *elementPtr, const GENERAL info, double *R)
{
	int i;							//counter
	int n2=2*info.noelement;		//2*number of elements

	for(i=0; i<n2; i++)				//vector elements 0..2n-1 =0
		R[i]=0;

	for(i=0; i<info.noelement; i++)	//vector elements 2n..3n-1
		if (info.linear == 0) 		//non-linear theory
			R[i+n2]=4*Pi*dot(elementPtr[i].u,elementPtr[i].normal);
		else		//#small angle approximation for beta and alpha
		{			//modified 6/15/03 G.B.
			R[i+n2]=4*Pi*norm2(elementPtr[i].u)*\
					(asin(elementPtr[i].normal[0])+\
 			         sin(elementPtr[i].nu)*info.beta+\
 			   		 cos(elementPtr[i].nu)*info.alpha);
//#			printf("small angle approximation for resultant vector!!\n");
		}
}
//===================================================================//
		//END assembly of resultant vector
//===================================================================//

//===================================================================//
		//START assembly of Part 1&2 of D-matrix
//===================================================================//
//assembles first part of D matrix that is the boundary conditions of
//the elementary wings.  Within a panel, two conditions exists tht are
//continuous vorticity magnitude and slope in spanwise direction.
//Across the panel boundaries the following conditions exist and are
//defined in panelPtr[j].BC1 and panelPtr[j].BC2 with a three digit
//number:
//
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
// BC = 100 - gamma = 0
// BC = 110 - gamma = 0 AND gamma' = 0  (free tip)
// BC = 220 - gamma[i] = gamma[j]  AND  gamma[i]' = gamma[j]'
// BC = 022 - gamma[i]' = gamma[j]'  AND  gamma[i]" = gamma[j]"
// BC = 010 - gamma' = 0
//========================================================================

void BoundaryCond(const BOUND_VORTEX *elementPtr, const PANEL *panelPtr, \
										const GENERAL info, double **D)
{
	int element=0;		//elementary wing index counter
	int panel;			//panel index
	int n,m;			//spanwise element, chordwise lifting line counter
	int col,row;		//column, row index
	int nextelement;	//index of neighboring element of panel to the right
	int k;				//counter

for(panel=0;panel<info.nopanel;panel++)			//loop over panels
	{
		//if it exists, the panel to the "right" is determined
		//then the index difference of the neighboring elementary wings
		//of the two panels is determined.
		if ((panelPtr[panel].right <= panel || \
			 panelPtr[panel].right>info.nopanel)\
			&& panelPtr[panel].right!=0)
			{	printf("INPUT FILE ERROR! Panel numbering messed up!!\n");
				printf("Indices are not ascending\n");
				exit(1);
			}
			//error message if panel numbering is not ascending or out of bound
		else
		{
			nextelement=0;		//initializing
			//adds number of spanwise elements up to right panel
			for (k=panel;k<(panelPtr[panel].right-1);k++)
				nextelement+=panelPtr[k].n;
			nextelement*=(info.m-1);
			nextelement++;
		}

		for(m=0;m<info.m;m++)		//loop over chordwise lifting lines
		{
			col = 3*(element);			//column of D matrix
			switch (panelPtr[panel].BC1)	//edge 1 bound. cond.
			{
				case 100: 	//gamma = 0
					D[element][col]   = -1;
					D[element][col+1] = elementPtr[element].eta;
					D[element][col+2] = -elementPtr[element].eta*\
										elementPtr[element].eta;
//#printf("case 100 element %d  nu %lf\n",element,elementPtr[element].eta);//#
				break;

				//added 5/13/04 G.B.
				case 110: 	//gamma = 0 and gamma' = 0 (free tip)
					D[element][col]   = -1;
					D[element][col+1] = elementPtr[element].eta;
					D[element][col+2] = -elementPtr[element].eta*\
										elementPtr[element].eta;
					row = element+2*info.noelement;
					D[row][col]   = 0;
					D[row][col+1] = 1;
					D[row][col+2] = -2*elementPtr[element].eta;
//#printf("case 110 element %d  nu %lf\n",element,elementPtr[element].eta);//#
				break;

				case 220:	//gamma[i] = gamma[j] AND gamma[i]' = gamma[j]'
					col = 3*(element);			//column of D matrix
					D[element][col]   = -1;
					D[element][col+1] = elementPtr[element].eta;
					D[element][col+2] = -elementPtr[element].eta*\
										elementPtr[element].eta;

					row = element+info.noelement;//+info.m-2;
					D[row][col]   = 0;
					D[row][col+1] = -1;
					D[row][col+2] = 2*elementPtr[element].eta;
//#printf("case 220 element %d  nu %lf\n",element,elementPtr[element].eta);//#
				break;

				case 22:	//gamma[i]' = gamma[j]' AND gamma[i]" = gamma[j]"
					col = 3*(element);			//column of D matrix
					D[element][col]   = 0;
					D[element][col+1] = 0;
					D[element][col+2] = 2;

					row = element+info.noelement;//+info.m-2;
					D[row][col]   = 0;
					D[row][col+1] = -1;
					D[row][col+2] = 2*elementPtr[element].eta;
//#printf("case 22 element %d  nu %lf\n",element,elementPtr[element].eta);//#
				break;

				case 10:	//gamma' = 0
					col = 3*(element);			//column of D matrix
					D[element][col]   = 0;
					D[element][col+1] = -1;
					D[element][col+2] = 2*elementPtr[element].eta;
//#printf("case 10 element %d  nu %lf\n",element,elementPtr[element].eta);//#
				break;

				default:	//default case assumes gamma = 0
					col = 3*(element);			//column of D matrix
					D[element][col]   = -1;
					D[element][col+1] = elementPtr[element].eta;
					D[element][col+2] = -elementPtr[element].eta*\
										elementPtr[element].eta;
					printf("WARNING!! \nEdge 1  boundary condition of");
					printf(" of panel %d undefined!!\n",panel+1);
//#printf("default element %d  nu %lf\n",element,elementPtr[element].eta);//#
				break;
			}		// end switch statement for edge 1 bound. cond.

		for(n=1;n<(panelPtr[panel].n);n++) //loop over spanwise elements-1
			{
				element++;		//increase elementary-wing index by 1
				col = 3*element;
				row = element+info.noelement;//+info.m-2;

				if (panelPtr[panel].left == 0)
				row+=(info.m-m-1); 		//row correction for left-end panels

				//continuous vorticity magnitude
				D[element][col-3] = 1;
				D[element][col-2] = elementPtr[element-1].eta;
				D[element][col-1] = elementPtr[element-1].eta*\
									elementPtr[element-1].eta;
				D[element][col]   = -1;
				D[element][col+1] = elementPtr[element].eta;
				D[element][col+2] = -elementPtr[element].eta\
									*elementPtr[element].eta;

				//continuous vorticity slope
				D[row][col-3] = 0;
				D[row][col-2] = 1;
				D[row][col-1] = 2*elementPtr[element-1].eta;

				D[row][col]   = 0;
				D[row][col+1] = -1;
				D[row][col+2] = 2*elementPtr[element].eta;

			}//end loop over n, number of spanwise elementary wings

		//determine row to put right edge (edge 2) boundary condition
		if (panelPtr[panel].right!=0)
			 row = element+nextelement+\
				  (panelPtr[panelPtr[panel].right-1].n-panelPtr[panel].n)*m;
		else row = element+panelPtr[panel].n*(info.m-(m+1))+m+1;

		switch (panelPtr[panel].BC2)	//edge 2 bound. cond.
			{
				case 100: 	//gamma = 0
					if(panel<info.nopanel-1) row += info.noelement;
					D[row][col]   = 1;
					D[row][col+1] = elementPtr[element].eta;
					D[row][col+2] = elementPtr[element].eta*\
									elementPtr[element].eta;
				break;

				//added 5/13/04 G.B.
				case 110: 	//gamma = 0 and gamma'=0
					if(panel<info.nopanel-1) row += info.noelement;
					D[row][col]   = 1;
					D[row][col+1] = elementPtr[element].eta;
					D[row][col+2] = elementPtr[element].eta*\
									elementPtr[element].eta;
					row = element+2*info.noelement;
					D[row][col]   = 0;
					D[row][col+1] = 1;
					D[row][col+2] = 2*elementPtr[element].eta;
//#printf("case 110 element %d  nu %lf\n",element,elementPtr[element].eta);//#
				break;

				case 220:	//gamma[i] = gamma[j] AND gamma[i]' = gamma[j]'
					D[row][col]   = 1;
					D[row][col+1] = elementPtr[element].eta;
					D[row][col+2] = elementPtr[element].eta*\
										elementPtr[element].eta;
					row += info.noelement;//+info.m-2;
					D[row][col]   = 0;
					D[row][col+1] = 1;
					D[row][col+2] = 2*elementPtr[element].eta;
				break;

				case 22:	//gamma[i]' = gamma[j]' AND gamma[i]" = gamma[j]"
					D[row][col]   = 0;
					D[row][col+1] = 0;
					D[row][col+2] = 1;

					row += info.noelement;//+info.m-2;
					D[row][col]   = 0;
					D[row][col+1] = 1;
					D[row][col+2] = 2*elementPtr[element].eta;
				break;

				case 10:	//gamma' = 0
					D[row][col]   = 0;
					D[row][col+1] = 1;
					D[row][col+2] = 2*elementPtr[element].eta;
				break;

				default:	//default case assumes gamma = 0
					D[row][col]   = 1;
					D[row][col+1] = elementPtr[element].eta;
					D[row][col+2] = elementPtr[element].eta*\
									elementPtr[element].eta;
					printf("WARNING!! \nEdge 2  boundary condition");
					printf(" of panel %d undefined!!\n",panel+1);
				break;
			}		// end switch statement for edge 2 bound. cond.
			element++;	//increment index to next elementary wing
		}//end loop over m, one lifting line at a time
	} //end loop over panel, one panel at a time
}
//===================================================================//
		//END assembly of Part 1&2 of D-matrix
//===================================================================//

//===================================================================//
		//START function DVE_KinCond
//===================================================================//

	//This function computes the third part of the D-matrix for DVE's. It
	//computes the kinematic condition (tangency requirement at surface DVE)
	//at the reference point of a DVE. Specifiaclly, function computes the
	//influence coefficients a,b,c that are due to the DVE systems of the
	//lifting surface and the wake.  The coefficients multiplied with the
	//local bound vortex coefficients, A, B, and C, will yield the locally
	//induced velocity. In order to satisfy the kinematic condition, the
	//component that is normal to the element surface, cancels the normal
	//component of the local velocity that is composed of the free stream
	//and the part induced by the wake. This later part determines the
	//resutant vector.
	//The upper 2/3 of the D-matrix were determined previously in the
	//function BoundaryCond and hasn't changed.

void DVE_KinCond(DVE *surfacePtr,const GENERAL info,\
				 const PANEL *panelPtr,double **D)
{
	//input
	//  surfacePtr			geometric information on surface DVE's
	//	info				general info
	//
	//output
	//	D					D-matrix, third part, boundary conditions remain
	//						the same as in the previous part.

int i,j,m,panel,row,col;//counter
int k,l,span;
int imin,imax,element;
double a[3],b[3],c[3];			//combined influence coefficients
double tempS;
int type;  //DVE type (vortex sheet, LE&TE vortices)

//###!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!//
//###	Removed: 4/29/04, G. Bramesfeld										 //
//###!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!//
//###DVE *tempDVE;					//temporary DVE
//###ALLOC1D(&tempDVE,info.nospanelement);//temporary first post wing wake DVE
//###
//###//computes temporary wake DVE elements that arejust aft of trailing edge
//###Squirt_out_Wake(info,panelPtr,surfacePtr,tempDVE);
//###									//subroutine in wake_geometry.cpp

//loop over surface DVE's at whose control points the induced velocity
//is being computed.  Defines also row of D-matrix (3i)
element=0; // initializing element index counter
for(panel=0;panel<info.nopanel;panel++)	//loop over panels
{
  //checking if left or right side are free tips
  //if so, kinematic condition is not being satisfied at the tips, but
  //gamma'= 0, see boundary conditions
  if(panelPtr[panel].BC1==110)		imin = 1;
  else								imin = 0;
  if(panelPtr[panel].BC2==110)		imax = panelPtr[panel].n-1;
  else								imax = panelPtr[panel].n;

  for(m=0;m<info.m;m++)		//loop over chordwise lift. lines
  {
	//increase index by one if left edge is a free tip
	if(panelPtr[panel].BC1==110)	element ++;

	for(i=imin;i<imax;i++)	//loop over spanwise elements-1
	{
		row=element+2*info.noelement;

//##//loop over surface DVE's at whose control points the induced velocity
//##//is being computed.  Defines also row of D-matrix (3i)
//##	for(i=0;i<info.noelement;i++)
//##	{
//##		row=i+2*info.noelement;

		//loop over surface DVE's that induce velocity on element i
	//#		for(j=0;j<info.noelement;j++)
		j=0; //index of DVE that induces
		for(k=0;k<info.nopanel;k++)	//loop over panels
	  for(l=0;l<info.m;l++)		//loop over chordwise lift. lines
	  for(span=0;span<panelPtr[k].n;span++) //loop over panel span
		{
			//setting singfct of DVE j temporarily to zero; added 8/16/05 GB
			tempS = surfacePtr[j].singfct;
			surfacePtr[j].singfct = 0;

			if(l<info.m-1) type=0; //DVE has TE&LE vortex GB 8/10/21
			else type = 2; //if DVE is at trailing edge, no TE-vortex

 			//computes influence coefficients, a, b, and c, of
			//DVE j on reference point of DVE element
			DVE_Influence_Coeff\
					(surfacePtr[j],info,surfacePtr[element].xo,a,b,c,type);
									//subroutine in induced_velocity.cpp
					
//printf("ind %d PoI %d %lf  %lf %lf\n",j,element,a[0],a[1],a[2]);
//printf("ind %d PoI %d %lf  %lf %lf\n",j,element,b[0],b[1],b[2]);
//printf("ind %d PoI %d %lf  %lf %lf\n\n",j,element,c[0],c[1],c[2]);

			//reassigning singfct of DVE j; added 8/16/05 GB
			surfacePtr[j].singfct = tempS;

			col = 3*j;  //first column

			//Horstmann Eq.43
			D[row][col]		= dot(a,surfacePtr[element].normal);
			D[row][col+1]	= dot(b,surfacePtr[element].normal);
			D[row][col+2]	= dot(c,surfacePtr[element].normal);
//printf("row %d col %d %lf  %lf %lf\n",row,col,D[row][col],D[row][col+1],D[row][col+2]);

			j++; //next DVE that induces
		}//end loop j, element that induces vel. on element i
		element++; //next element index
	}	//end loop i, number of spanwise elements
	//increase index by one if right edge is a free tip
	if(panelPtr[panel].BC2==110)	element ++;
  }//end loop over m, number of chordwise lifting lines
}//end  loop over panel, number of panels

//##	}	//end loop i, element at whose ctrl pt the ind. vel is comp.
//###	FREE1D(&tempDVE,info.nospanelement);

/*  ###########################################################
 //save D matrix and resultant vector R in file D_matrix.txt
 FILE *fp;
 fp = fopen("D_matrix.txt", "w");
 int n;
 //writes header line
 fprintf(fp, "\t");
 for(m=0; m<info.noelement*3; m++)
 fprintf(fp, "%d\t",m);
 fprintf(fp, "\t\tR");
 for(n=0; n<info.noelement*3; n++)
 {
	 //row number
 	fprintf(fp, "\n%d\t",n);
 	//n-th row of D
 	for(m=0; m<info.noelement*3; m++)
 	fprintf(fp, "%lf\t",D[n][m]);
 	//n-th element of R
//	fprintf(fp, "\t\t%lf",R[n]);
 }
 fclose(fp);
 //###########################################################//*/

}
//===================================================================//
		//END function DVE_KinCond
//===================================================================//

//===================================================================//
		//START function KinCond
//===================================================================//

	//computes kinematic condition (induced flow opposite to free-stream flow
	//component through elementary-wing surface) at elementary wing control
	//point. Specifiaclly, function computes for each elementary wing
	//influence coefficients a,b,c, that are due to the bound and trailing
	//vorticity.  The coefficients multiplied with the local bound vortex
	//coefficients, A, B, and C, will yield the locally induced velocity.
	//In order to satisfy the kinematic condition, the component that is
	//normal to the element surface, cancels the normal component of the
	//free stream velocity. This defines the third part of the D-matrix
void KinematicCond(const BOUND_VORTEX *elementPtr,const GENERAL info,\
				   const PANEL *panelPtr,double **D)
{
	int element,i,j,row,panel,m;	//counter
	int imin,imax,col;
	double a[3],b[3],c[3];			//combined influence coefficients
	double normal[3];

//####
//loop over surface DVE's at whose control points the induced velocity
//is being computed.  Defines also row of D-matrix (3i)
element=0; // initializing element index counter
for(panel=0;panel<info.nopanel;panel++)	//loop over panels
{
  //checking if left or right side are free tips
  //if so, kinematic condition is not being satisfied at the tips, but
  //gamma'= 0, see boundary conditions
  if(panelPtr[panel].BC1==110)		imin = 1;
  else								imin = 0;
  if(panelPtr[panel].BC2==110)		imax = panelPtr[panel].n-1;
  else								imax = panelPtr[panel].n;


  for(m=0;m<info.m;m++)		//loop over chordwise lift. lines
  {
	//increase index by one if left edge is a free tip
	if(panelPtr[panel].BC1==110)	element ++;

	for(i=imin;i<imax;i++)	//loop over spanwise elements-1
	{
		row=element+2*info.noelement;

		//loop over elements that induce velocity on element i
		for(j=0;j<info.noelement;j++)
		{
 			//computes influence coefficients, a, b, and c, of
			//bound vortex j and its fixed wake on control point i
			Influence_Coeff(elementPtr[j],info,elementPtr[element].xA,a,b,c);
									//subroutine in induced_velocity.cpp

			//KHH small angle approach, not needed
			//changed 6/15/03 G.B.
			if (info.linear == 1)
			{
				normal[0] = asin(elementPtr[element].normal[0]);
				normal[1] = -sin(elementPtr[element].nu);
				normal[2] = cos(elementPtr[element].nu);
			}
			else	//standard method
			{
				normal[0] = elementPtr[element].normal[0];
				normal[1] = elementPtr[element].normal[1];
				normal[2] = elementPtr[element].normal[2];
			}
			col = 3*j;
			//Horstmann Eq.43
			D[row][col]		= dot(a,normal);
			D[row][col+1]	= dot(b,normal);
			D[row][col+2]	= dot(c,normal);
		}//next element j
		element++; //next element index
	}	//end loop i, number of spanwise elements
	//increase index by one if right edge is a free tip
	if(panelPtr[panel].BC2==110)	element ++;
  }//end loop over m, number of chordwise lifting lines
}//end  loop over panel, number of panels

//####
/*///#########
	//loop over elements at whose control points the induced velocity
	//is being computed.  Defines also row of D-matrix (3i)
	for(i=0;i<info.noelement;i++)
	{
		row=i+2*info.noelement;

		//loop over elements that induce velocity on element i
		for(j=0;j<info.noelement;j++)
		{
 			//computes influence coefficients, a, b, and c, of
			//bound vortex j and its fixed wake on control point i
			Influence_Coeff(elementPtr[j],info,elementPtr[i].xA,a,b,c);
									//subroutine in induced_velocity.cpp

//KHH small angle approach, not needed
//changed 6/15/03 G.B.
if (info.linear == 1)
{
	normal[0] = asin(elementPtr[i].normal[0]);
	normal[1] = -sin(elementPtr[i].nu);
	normal[2] = cos(elementPtr[i].nu);
}
else	//standard method
{
	normal[0] = elementPtr[i].normal[0];
	normal[1] = elementPtr[i].normal[1];
	normal[2] = elementPtr[i].normal[2];
}

			//Horstmann Eq.43
			D[row][3*j]		= dot(a,normal);
			D[row][3*j+1]	= dot(b,normal);
			D[row][3*j+2]	= dot(c,normal);

//#printf("a: %lf  %lf  %lf\n",a[0],a[1],a[2]);
//#printf("b: %lf  %lf  %lf\n",b[0],b[1],b[2]);
//#printf("c: %lf  %lf  %lf\n",c[0],c[1],c[2]);

//#printf("D[%d][%d]: %lf \tD[%d][%d]: %lf\tD[%d][%d]: %lf\n",\
//#row,3*i,D[row][3*i],row,3*i+1,D[row][3*i+1],row,3*i+2,D[row][3*i+2]);//#
		}//end loop j, element that induces vel. on element i
	}	//end loop i, element at whose ctrl pt the ind. vel is comp.
////#########*/

}
//===================================================================//
		//END function KinCond
//===================================================================//

//===================================================================//
		//START function Trailing_Edge_Vorticity_Distribution
//===================================================================//

	//Adds up bound vorticity at trailing edge.  Only vorticity in span
	//direction, i.e. in a plane perpendicular to the free stream, is
	//considered. Just as the bound vorticity, the trailing edge vortex
	//ellements have a parabolic circulation distribution that is
	//described with three coefficients, A, B, and C.

void Trailing_Edge_Vorticity(const GENERAL info,const PANEL *panelPtr,\
					BOUND_VORTEX *elementPtr,BOUND_VORTEX *trailedgePtr)
{
//#FILE *fp;

int k,n, m;		//counter
int l=0;		//index of trailing edge element
int io=0;		//index of first bound vortex element of panel
int i;			//index of bound vortex element

	for(k=0;k<info.nopanel;k++)			//loop over panels
	{
		for(n=0;n<panelPtr[k].n;n++)	//loop over panel span
		{
			//initializing coefficients
			trailedgePtr[l].A=0;
			trailedgePtr[l].B=0;
			trailedgePtr[l].C=0;

			for(m=0;m<info.m;m++)//loop over panel chord
			{
				//index of current bound vortex
				i=io+n+m*panelPtr[k].n;

				//adding up spanwise circulation
				trailedgePtr[l].A+= elementPtr[i].A;
				trailedgePtr[l].B+= elementPtr[i].B;
				trailedgePtr[l].C+= elementPtr[i].C;
			}
			l++;						//next trailing edge element
		}
		io += panelPtr[k].n*info.m;		//next first bound vortex
	}

//#	//save circulation coefficients of trailing edge in file
//#	fp = fopen("ciruculation distr.txt", "w");
//#	//writes header
//#	fprintf(fp, "n  \tA\tB\t\tC\t\tphi\t\tnu\n");
//#
//#	l=0;
//#	io=0;
//#
//#	for(k=0;k<info.nopanel;k++)			//loop over panels
//#	{
//#		for(n=0;n<panelPtr[k].n;n++)	//loop over panel span
//#		{
//#			for(m=0;m<info.m;m++)//loop over panel chord
//#			{
//#				i=io+n+m*panelPtr[k].n;	//index of bound vortex
//#				//bound vortex element index
//#				fprintf(fp, "%d  ",i);
//#				//A,B,C coefficients
//#				fprintf(fp, "%lf\t%lf\t%lf\t%lf\t%lf\n",\
//#				elementPtr[i].A,elementPtr[i].B,elementPtr[i].C,\
//#				elementPtr[i].phi*180/Pi,elementPtr[i].nu*180/Pi);
//#			}
//#			fprintf(fp, "\n\ttrailing edge vortices\n");
//#			//trailing edge vortex element index
//#			fprintf(fp, "%d  ",l);
//#			//A,B,C coefficients
//#			fprintf(fp, "%lf\t%lf\t%lf\t%lf\t%lf\n\n\n",\
//#			trailedgePtr[l].A,trailedgePtr[l].B,trailedgePtr[l].C,\
//#			trailedgePtr[l].phi*180/Pi,trailedgePtr[l].nu*180/Pi);
//#
//#			l++;
//#		}
//#		io += panelPtr[k].n*info.m;		//next first bound vortex
//#	}
//#	fclose(fp);
}
//===================================================================//
		//END function Trailing_Edge_Vorticity_Distribution
//===================================================================//
