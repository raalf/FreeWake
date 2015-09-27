void GaussSolve(double **,double *, const int,double *);
void LU_Decomposition(double **, int, int *);
void LU_Solver(double **,const int,const int *,double *);

//===================================================================//
		//START GaussSolve
//===================================================================//
/*
   Solve a system of n equations in n unknowns using Gaussian Elimination
   Solve an equation in matrix form Ax = b
   The 2D array a is the matrix A with an additional column b.
   This is often written (A:b)

   A0,0    A0,1    A0,2    ....  A0,n-1     b0
   A1,0    A1,1    A1,2    ....  A1,n-1     b1
   A2,0    A2,1    A2,2    ....  A2,n-1     b2
   :       :       :             :          :
   :       :       :             :          :
   An-1,0  An-1,1  An-1,2  ....  An-1,n-1   bn-1

   The result is returned in x, otherwise the function returns FALSE
   if the system of equations is singular.
*/
void GaussSolve(double **A, double *R, const int n,double *x)
{
   int i,j,k,maxrow;
   double tmp, **a;

	//allocates memory for a
	ALLOC2D(&a,(n),(n+1));

   //assigns a=A|R
   for (i=0;i<n;i++)
   {
	   a[i][n]=R[i];
	   for(j=0;j<n;j++)
	   	   a[i][j]=A[i][j];
   }

   for (i=0;i<n;i++)
   {
      /* Find the row with the largest first value */
      maxrow = i;
      for (j=i+1;j<n;j++)
      {
         if (fabs(a[j][i]) > fabs(a[maxrow][i]))
            maxrow = j;
      }

      /* Swap the maxrow and ith row */
      for (k=i;k<n+1;k++)
      {
         tmp = a[i][k];
         a[i][k] = a[maxrow][k];
         a[maxrow][k] = tmp;
      }

      /* Singular matrix? */
      if ((a[i][i]*a[i][i]) < DBL_EPS)
      {
         printf("\nWARNING! Singular matrix - ");
         printf("no unique solution in function GaussSolve\n");
         printf("a[%d][%d] = %lf  a^2[%d][%d] =%.20lf\n",\
         					i,i,a[i][i],i,i,a[i][i]*a[i][i]);
	 }


      /* Eliminate the ith element of the jth row */
      for (j=i+1;j<n;j++)
      {
         for (k=n;k>=i;k--)
         {
            a[j][k] -= a[i][k] * a[j][i] / a[i][i];
         }
      }
   }//end loop over i, rows of a

//#for(i=0;i<n;i++)				//#
//#{	for(j=0;j<n+1;j++)		//#
//#		printf("%lf\t",a[i][j]);//#
//#  	printf("\n");			//#
//#}							//#


   /* Do the back substitution */
   for (j=n-1;j>=0;j--)
   {
      tmp = 0;
      for (k=j+1;k<n;k++)
         tmp += a[j][k] * x[k];
      x[j] = (a[j][n] - tmp) / a[j][j];
   }
   //frees allocated memory
   FREE2D(&a,(n),(n+1));
}
//===================================================================//
		//END GaussSolve
//===================================================================//

//===================================================================//
		//START LU_Decomposition
//===================================================================//
void LU_Decomposition(double **a, int n, int *indx)
{
//
// Decomposes matrix a in an upper and lower matrix using. The lower and
// upper coefficients are saved to a and the original a values are lost!!
//
// a0,0    a0,1    a0,2  ... a0,n-1		   U0,0    U0,1    U0,2  ... U0,n-1
// a1,0    a1,1    a1,2  ... a1,n-1		   L1,0    U1,1    U1,2  ... U1,n-1
// a2,0    a2,1    a2,2  ... a2,n-1    =>  L2,0    L2,1    U2,2  ... U2,n-1
// :       :       :         :	     	   :       :       :         :
// :       :       :         :	    	   :       :       :         :
// an-2,0  an-2,1  an-2,2... an-2,n-1	   Ln-2,0  Ln-2,1  Ln-2,2... Un-2,n-1
// an-1,0  an-1,1  an-1,2... an-1,n-1	   Ln-1,0  Ln-1,1  Ln-1,2... Un-1,n-1

// The result is returned in A. If the system of equations is singular, the
// program exits.
//
// The original method supplied with the book had messed up indices that
// started at value 1 (as it is done for fortran).  I tried to fix that
// as much as possible and it seems to work now. G.B. 11/8/03
//
//input
//  a		matrix nxn
//  indx	array that holds information of pivoting a
//  n		size of matrix
//
//
/* (C) Copr. 1986-92 Numerical Recipes Software 0K.s@%. */

	int i,imax,j,k;
	double big,dum,sum,temp;
	double *vv;


/*/##########################################
 FILE *fp;
 fp = fopen("D_matrix.txt", "a");
 fprintf(fp, "\n\n");
 fprintf(fp, "the new D\n");
 for(i=0; i<n; i++)
 {
	 //row number
	 fprintf(fp, "%d\t",i);
	 //n-th row of D
 	for(j=0; j<n; j++)
 		fprintf(fp, "%lf\t",a[i][j]);
 	fprintf(fp,"\n");
 }
 fclose(fp);
//#######################################*///

	ALLOC1D(&vv,n);		//temp scaling array

	for (i=0;i<n;i++)
	{
		big=0.0;
		for (j=0;j<n;j++)
			if ((temp=fabs(a[i][j])) > big) big=temp;
		if (big == 0.0) printf("\nSingular matrix in routine ludcmp");
		vv[i]=1.0/big;
	}
	for (j=1;j<=n;j++)
	{
		for (i=0;i<j-1;i++)
		{
			sum=a[i][j-1];
			for (k=0;k<i;k++) sum -= a[i][k]*a[k][j-1];
			a[i][j-1]=sum;
		}
		big=0.0;
		for (i=j;i<=n;i++) {
			sum=a[i-1][j-1];
			for (k=1;k<j;k++)
				sum -= a[i-1][k-1]*a[k-1][j-1];
			a[i-1][j-1]=sum;
			if ( (dum=vv[i-1]*fabs(sum)) >= big) {
				big=dum;
				imax=i;
			}
		}
		if (j != imax)
		{
			for (k=1;k<=n;k++)
			{
				dum=a[imax-1][k-1];
				a[imax-1][k-1]=a[j-1][k-1];
				a[j-1][k-1]=dum;
			}
			vv[imax-1]=vv[j-1];
		}
		indx[j-1]=imax;
		if (a[j-1][j-1] == 0.0) a[j-1][j-1]=0.00001;
		if (j != n) {
			dum=1.0/(a[j-1][j-1]);
			for (i=j+1;i<=n;i++) a[i-1][j-1] *= dum;
		}
	}

	FREE1D(&vv,n);			//free memory

/*##########################################
 fp = fopen("D_matrix.txt", "a");
 fprintf(fp, "\n\n");
 fprintf(fp, "the decomposed D\n");
 for(i=0; i<n; i++)
 {
	 //row number
	 fprintf(fp, "%d\t",i);
	 //n-th row of D
 	for(j=0; j<n; j++)
 		fprintf(fp, "%lf\t",a[i][j]);
 	fprintf(fp,"\t%d\n",indx[i]);
 }
 fclose(fp);
//#######################################*///

}
//===================================================================//
		//END LU_Decomposition
//===================================================================//

//===================================================================//
		//START LU_Solver
//===================================================================//
void LU_Solver(double **a,const int n,const int *indx,double *b)
{
//Solves L*(U*x) = b equation system, Lower and upper matrix are
//contained in a.  indx contains the pivoting information for b.
//The result is returned in b, otherwise the function returns a warning
//if the system of equations is singular.
//
// The LU-solver is described in detail in "Numerical Recipes",
// W.H.Press, et al., Cambridge University Press, Cambridge
//
// The original method supplied with the book had messed up indices that
// started at value 1 (as it is done for fortran).  I tried to fix that
// as much as possible and it seems to work now. G.B. 11/8/03
//
/* (C) Copr. 1986-92 Numerical Recipes Software 0K.s@%. */
//
// input:
//  a		matrix nxn
//  n		size of matrix
//  indx	array that holds information of pivoting a
//  b		rhs vector, holds results after execution of this routine
//
	int i,ii=0,ip,j;
	double sum;

	for (i=1;i<=n;i++)
	{
		ip=indx[i-1];
		sum=b[ip-1];
		b[ip-1]=b[i-1];
		if (ii)
			for (j=ii;j<=i-1;j++) sum -= a[i-1][j-1]*b[j-1];
		else if (sum) ii=i;
		b[i-1]=sum;
	}

//#printf("x\n");
//#for (i=0;i<n;i++)
//#	printf(" %lf",b[i]);
//#printf("\n");

	for (i=n-1;i>=0;i--)
	{
		sum=b[i];
		for (j=i+1;j<n;j++) sum -= a[i][j]*b[j];
		b[i]=sum/a[i][i];
	}
}
//===================================================================//
		//END LU_Solver
//===================================================================//
