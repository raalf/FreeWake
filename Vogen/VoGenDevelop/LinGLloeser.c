#include "LinGLloeser.h"
#include "stdio.h"

double dabs(x)
double x;
{
  if(x >= 0.)
    return x;
  return (-x);
}

double det (a, pivot, coln)
int pivot [], coln;
float *a;
{
  int i, sign;
  double d;
  
  sign = 0;
  d = 1.;
  for (i = 0; i < coln; i++)
  {
  if (pivot [i] != i)
    sign++;
    d  *= a INDEX(i,i);
  }
  sign= sign - ((sign >> 1) << 1);
  if (sign)
    d = -d;
  return d;
}

void sscal(n, sa, sx, incx)
int n, incx;
float sa, *sx;
{
  int i, nincx;

  if (n<=0)
    return;
  nincx = incx*n;
  DOV(i, nincx, incx)
    sx[i]=sx[i]*sa;
  return;
}

float sasum (n, sx, incx)
int incx, n;
float *sx;
{
  int i, nincx;
  double stemp;

  stemp=0.0;
  nincx=n*incx;
  if(n<=0)
    return 0.0;
  DOV(i, n, nincx)
    stemp =stemp+abs(sx[i]);
  return (stemp);
}

void saxpy (n, sa, sx, incx, sy, incy)
int n, incx, incy;
float sa, *sx, *sy;
{
  int i, iy, ix;

/*sy =sa * sx +sy */
  if(n<=0)
    return;
  if (sa == 0.)
    return;

  iy = ix =0;
  if (incx<0)
    ix=incx*(1-n);

  if (incy<0)
    iy=incy*(1-n);

  DOFOR(i,n)
  {
    sy [iy] =sy[iy]+sa*sx[ix];
    iy+= incy;
    ix+= incx;
  }
  return;
}

void sswap (n, sx, incx, sy, incy)
int n, incx, incy;
float *sx, *sy;

{
  int ix, iy, i;
  float t;
  
  if (n<=0) return;
  ix = iy =0;
  if (incx <0)
    ix=incx*(1-n);
  if (incy <0)
    iy=incy*(1-n);

  DOFOR(i,n);
  {
    t=sx[ix];
    sx[ix]=sy[iy];
    sy[iy]=t;
    ix +=incx;
    iy +=incy;
  }
  return;
}


double deti (a, pivot, coln, expon)
int pivot [], coln, *expon;
float *a;
{
  int i, sign;
  double d, dabs();

  sign = 0;
  d = 1.;
  *expon = 0;
  for (i = 0; i < coln; i++)
  {
    if (pivot [i] != i)
      sign++;
    d *= a INDEX(i,i);
    if (dabs (d) > 10.)
    {
      (*expon)++;
      d *= .1; 
    }
    else if (dabs (d) < .1)
    {
      (*expon)--;
      d *= 10;
    }
  }
  sign= sign - ((sign >> 1) << 1);
  if (sign)
  d = -d;
  return d;
}

void invm (a, coln, n, pivot, work)
int coln, n, *pivot; 
float *a, *work;
{
  int i, j, k, l, kb, kp1, nm1;
  float t;
  
  nm1 = n - 1;
  /* keine Berechnung der Determinanten */
  /* Invertieren von R */
  DOFOR(k,n)
  {
    a INDEX(k,k) = t = 1. / a INDEX(k,k);
    t = -t;
    sscal (k, t, &(a INDEX(0,k)), coln);
    kp1 = k + 1;
    if (nm1 >= kp1)
    {
      DOBYYY(j,kp1,n)
      {
        t = a INDEX(k,j);
        a INDEX(k,j) = 0.0;
        saxpy (k+1, t, &(a INDEX(0,k)), coln, &(a INDEX(0,j)), coln);
      }
    }
  }
  /* inv(R) * inv(L) */
  if (nm1 >= 1)
  {
    DOFOR(kb,nm1)
    {
      k = nm1 - kb - 1;
      kp1 = k + 1;
      DOBYYY(i,kp1,n)
      {
        work [i] = a INDEX(i,k);
        a INDEX(i,k) = 0.0;
      }
      DOBYYY(j,kp1,n)
      {
        t = work [j];
        saxpy (n, t, &(a INDEX(0,j)), coln, &(a INDEX(0,k)), coln);
      }
      l = pivot [k];
      if (l != k)
        sswap (n, &(a INDEX(0,k)), coln, &(a INDEX(0,l)), coln);
    }
  }
  return; 
}

void lufact (a,coln,n,pivot,info)
float *a;
int coln, n, *pivot, *info;
{
  int i, j, k, l, kp1, nm1, last;
  float t;

  *info=0;
  nm1=n-1;
  if (nm1>=1)
  {
    DOFOR(k,nm1)
    {
      kp1=k+1;
      pivot[k]=l=isamax((n-k), &(a INDEX(k,k)), coln)+k;
      if (a INDEX(l,k)!=0.)
      {
        if (l !=k)
        {
          t=a INDEX(l,k);
          a INDEX(l,k)=a INDEX(k,k);
          a INDEX(k,k)=t;
        }
        t=-1./a INDEX(k,k);
        sscal (nm1-k, t, &(a INDEX(k+1,k)), coln);
        DOBYYY(j,kp1,n)
        {
          t=a INDEX(l,j);
          if(l!=k)
          {
            a INDEX(l,j)=a INDEX(k,j);
            a INDEX(k,j)=t;          
          }
          saxpy (nm1-k, t, &(a INDEX(k+1,k)), coln, &(a INDEX(k+1,j)), coln);
        }
      }
      else
      {
        *info=k;
      }
    }
  }
  pivot  [nm1]=nm1;
  if (a INDEX(nm1, nm1)==0.0)
    *info=nm1;
  return;
}

int isamax(n, sx, incx)
int n, incx;
float *sx;
{
  int maxi, ix, i;
  float temp, smax;

  if  (n<=0)
    return -1;
  if (n==0)
    return 0;
  /*ix=0*/
  maxi=0;
  smax=abs(sx[0]);
  ix=incx; /*ix=ix + incx=incx*/
  DFOR(i, 2, n)
  {
    temp = abs(sx[ix]);
    if (temp>smax)
    {
      smax=temp;
      maxi=i;
    }
    ix+= incx;
  }
  return maxi;
}

float sdot (n, sx, incx, sy, incy)
int n, incx, incy;
float *sx, *sy;
{
  int i, ix, iy;
  float stemp;
  
  if (n<=0)
    return (0.);
  
  iy =ix =0;
  stemp=0.0;
  if (incx<0) ix=incx*(1-n);
  if (incy<0) iy=incy*(1-n);

  DOFOR(i,n)
  {
    stemp +=sy [iy]*sx[ix];
    iy+=incy;
    ix+=incx;
  }
  return stemp;
}

void printm (a, coln, rown, col, row)
int rown, row, col, coln;
float a[];
{
  int i, j, btm, top, count;
  printf("\n");
  btm=top=0;
  while (btm<col)
  {
    top=minLinGLloeser(col,(btm+8));
    printf("Ausgabe der Matrixspalten %d bis %d\n", btm, (top-1));
    DOFOR(j, row)
    {
      for (i=btm; i<top;i++)
      {
        printf(" %e ", a INDEX(j,i));
      }
      printf("\n");
    }
    btm+=8;
  }
  return;
}

void backsub (a, coln, n, pivot, b)
int coln, n, *pivot;
float *a, *b;
{
  int k, l, nm1;
  float t;

  nm1 = n - 1;
  DOFOR(k,n)
  {
    l = pivot [k];
    t = b [l];
    if (l != k)
    {
      b [l] = b [k];
      b [k] = t;
    }
    saxpy (nm1-k, t, &(a INDEX(k+1,k)), coln, &(b[k+1]), 1);
  }
  /* Rx = y lösen */
  DOFOR(l,n)
  {
    k = nm1 - l;
    b [k] = b [k] / a INDEX(k,k);
    t = -b [k]; 
    saxpy (k, t, &(a INDEX(0,k)), coln, b, 1);
  }
  return;
}

void backt (a, coln, n, pivot, b)
int coln, n, *pivot;
float *a, *b;
{
/* vie backsub, außer daß a'x = b gelöst vird */
  int k, l, nm1;
  float t, sdot ();

  nm1 = n - 1;
  /* zuerst R'y = b lösen */
  DOFOR(k,n)
  {
    t = sdot (k, &(a INDEX(0,k)), coln, b, 1); 
    b [k] = (b [k] - t) / a INDEX(k,k);
  }
  /* L^x = y lösen */
  if (nm1 < 1)
    return;
  for (k = nm1 - 1; k >= 0; k--)
  {
    b [k] += sdot (nm1-k, &(a INDEX(k+1,k)), coln, &(b[k+1]), 1);
    l = pivot [k];
    if (l != k)
    {
      t = b[l];
      b [l] = b [k];
      b [k] = t;
    }
  }
  return;
}

void LinGlLoesung (a, b, n)
int n;
float *a, *b;
{
  int coln,pivot[n],info,i;
  double det (), deti (), determ;
  coln = n;
  lufact (a, coln, n, pivot, &info);
  determ = det (a, pivot, coln);
  determ = deti (a, pivot, coln, &i);
  backsub (a, coln, n, pivot, b);
}
