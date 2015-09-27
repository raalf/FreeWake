//***************************************************************************
#ifndef ALLOC_H
#define ALLOC_H
//***************************************************************************
//#include <iostream.h>
//#include <stdlib.h>
//***************************************************************************
static double my_memory_allocated = 0.0;
static double my_memory_deleted = 0.0;
#define GLOBAL_VAR_FOR_ALLOCATION my_memory_allocated
#define GLOBAL_VAR_FOR_DELETION my_memory_deleted
/*
extern double GLOBAL_VAR_FOR_ALLOCATION;
extern double GLOBAL_VAR_FOR_DELETION;
*/
//***************************************************************************
#define IS_ZERO(x) (fabs(x) < 1.0e-8)
void myprintf(char *fmt, ...);
//***************************************************************************
inline double TotalMemoryAllocated()  // in bytes
{
return GLOBAL_VAR_FOR_ALLOCATION;
}
//***************************************************************************
inline double TotalMemoryDeleted()   // in bytes
{
return GLOBAL_VAR_FOR_DELETION;
}
//***************************************************************************
inline double TotalMemoryConsumed()  // in bytes
{
return (GLOBAL_VAR_FOR_ALLOCATION-GLOBAL_VAR_FOR_DELETION);
}
//***************************************************************************
template <class Etype>
inline void MY_MALLOC(Etype **ptr, int n = 1)
{
Etype *a = NULL;
a = new Etype[n];
GLOBAL_VAR_FOR_ALLOCATION += 1.0*n*sizeof(Etype);

if (a == NULL)
	{
	 printf("MY_MALLOC: Insufficient memory, could not allocate %d bytes!\n",sizeof(Etype)*n);
	 exit(-1);
	}

*ptr = a;
}
//***************************************************************************
template <class Etype>
inline void ALLOC1D(Etype **ptr, int m = 1)
{
Etype *a = NULL;

#ifdef ANIDEBUG
printf("Allocating %d x %d = %d bytes\n",m,sizeof(Etype),m*sizeof(Etype));
#endif

a = new Etype[m];
GLOBAL_VAR_FOR_ALLOCATION += 1.0*m*sizeof(Etype);
if (a == NULL)
	{
	 printf("ALLOC1D: Insufficient memory, could not allocate %d bytes!\n",sizeof(Etype)*m);
//	 cout << "ALLOC1D: Insufficient memory, could not allocate " << sizeof(Etype)*m << " bytes!" << endl;
	 exit(-1);
	}

*ptr = a;

#ifdef ANIDEBUG
printf("done\n");
//cout << "done!" << endl << flush;
#endif
}
//***************************************************************************
template <class Etype>
inline void ALLOC2D(Etype ***ptr, int m, int n)
{
Etype **a = NULL;
a = new Etype*[m];
GLOBAL_VAR_FOR_ALLOCATION += 1.0*m*sizeof(Etype*);

if (a == NULL)
	{
	 printf("ALLOC2D: Insufficient memory, could not allocate %d bytes!\n",sizeof(Etype*)*m);
//	 cout << "ALLOC2D: Insufficient memory, could not allocate " << sizeof(Etype*)*m << " bytes!" << endl;
	 exit(-1);
	}

for (int i = 0; i < m; i++)
	ALLOC1D(&a[i], n);

*ptr = a;
}
//***************************************************************************
template <class Etype>
inline void ALLOC3D(Etype ****ptr, int m, int n, int o)
{
Etype ***a = NULL;
a = new Etype**[m];
GLOBAL_VAR_FOR_ALLOCATION += 1.0*m*sizeof(Etype**);

if (a == NULL)
	{
	 printf("ALLOC3*D: Insufficient memory, could not allocate %d bytes!\n",sizeof(Etype**)*m);
//	 cout << "ALLOC3D: Insufficient memory, could not allocate " << sizeof(Etype**)*m << " bytes!" << endl;
	 exit(-1);
	}

for (int i = 0; i < m; i++)
	ALLOC2D(&a[i], n, o);

*ptr = a;
}
//***************************************************************************
template <class Etype>
inline void ALLOC4D(Etype *****ptr, int m, int n, int o, int p)
{
Etype ****a = NULL;
a = new Etype***[m];
GLOBAL_VAR_FOR_ALLOCATION += 1.0*m*sizeof(Etype***);

if (a == NULL)
	{
	 printf("ALLOC4D: Insufficient memory, could not allocate %d bytes!\n",sizeof(Etype***)*m);
//	 cout << "ALLOC4D: Insufficient memory, could not allocate " << sizeof(Etype***)*m << " bytes!" << endl;
	 exit(-1);
	}

for (int i = 0; i < m; i++)
	ALLOC3D(&a[i], n, o, p);

*ptr = a;
}
//***************************************************************************
template <class Etype>
inline void FREE1D(Etype **ptr, int m = 1)
{

Etype *a = *ptr;

if (a == NULL) return;

delete [] a;

GLOBAL_VAR_FOR_DELETION += 1.0*m*sizeof(Etype);

*ptr = NULL;
}
//***************************************************************************
template <class Etype>
inline void FREE2D(Etype ***ptr, int m, int n)
{
Etype **a = *ptr;

if (a == NULL) return;

for (int i = 0; i < m; i++)
	FREE1D(&a[i], n);

delete [] a;
GLOBAL_VAR_FOR_DELETION += 1.0*m*sizeof(Etype*);
*ptr = NULL;
}
//***************************************************************************
template <class Etype>
inline void FREE3D(Etype ****ptr, int m, int n, int o)
{
Etype ***a = *ptr;

if (a == NULL) return;

for (int i = 0; i < m; i++)
	FREE2D(&a[i], n, o);

delete [] a;
GLOBAL_VAR_FOR_DELETION += 1.0*m*sizeof(Etype**);
*ptr = NULL;
}
//***************************************************************************
template <class Etype>
inline void FREE4D(Etype *****ptr, int m, int n, int o, int p)
{
Etype ****a = *ptr;

if (a == NULL) return;

for (int i = 0; i < m; i++)
	FREE3D(&a[i], n, o, p);

delete [] a;
GLOBAL_VAR_FOR_DELETION += 1.0*m*sizeof(Etype***);
*ptr = NULL;
}
//***************************************************************************
#endif // ALLOC_H
//***************************************************************************
