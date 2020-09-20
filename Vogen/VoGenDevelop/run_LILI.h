/******************run_LILI.h********************                  
*                                                *      
*      Vorentwurfs-Generator "VoGen"             *          
*      Version 2.00                              *
*      Date of last modification 22.04.2020      *
*                                                *          
*      Writen by: Jan Himisch                    *      
*             jan.himisch@dlr.de                 *      
*                                                *      
*                                                *
*************************************************/
#ifndef _RUNLILI_H_
#define _RUNLILI_H_

#include "input.h"
//#define MaxZeilenLaenge 400

#define ZR_NULL    {int n; for (n=0;n < MaxZeilenLaenge;n++) {ZEILE_READ[n] = '\0';}}
#define nfgets(x)  {int i; ZR_NULL; for (i=1 ; i<=x; i++) {fgets (ZEILE_READ,MaxZeilenLaenge,fresult);}}

struct GEO_OUT
{
  double X1;
  double Y1;
  double Z1;
  double L1;  
  double X2;
  double Y2;
  double Z2;
  double L2;
  double XA;
  double YA;
  double ZA;
  double LA;      
};

struct CIRC 
{
  double XA;
  double YA;
  double ZA;
  double LA;
  double CNORM01;
  double CNORM02;
  double CNORM0;
  double DNORM01;
  double DNORM02;
  double DNORM0;
  double CNORM;
  double CNORM1;
  double CNORM2;
  double CWI;
};

struct FluegelBeiwerte
{
  double CA; 
  double CQ; 
  double CWI;
  double CL;
  double CM;
  double CN;
};

struct FALL
{
  struct CIRC *beiwerte;
  double Alfa;
  double Beta;
  double CA;  
  double CQ; 
  double CWI;
  double CL;
  double CM;
  double CN;
  struct FluegelBeiwerte *FluegelBw; 
};

struct file14
{
  double LREF_CL;
  double LREF_CM;
  double LREF_CN;
  double BEZUGSFLAECHE; 
  double GRUNDRISSFLAECHE; 
  double ABGEWICKELTE_FLAECHE;  
  double BEZUGSSPANNWEITE;
  double PROJIZIERTE_SPANNWEITE;    
  double BEZUGSSTRECKUNG;    
    
  struct GEO_OUT *GEODAT;
  struct FALL *Fall;
};

struct PltDist
{
 double  X;
 double  Y;
 double  Z;
 double  S;
 double  LOCAL_CIRCULATION;
 double  CFX_FROM_CDI;
 double  CFY;
 double  CFZ;
 double  CFNORMAL;
 double  CMX;
 double  CMX_FROM_CFY;
 double  CMX_FROM_CFZ;
 double  CMY;
 double  CMY_FROM_CDI;
 double  CMY_FROM_CFZ;
 double  CMZ;
 double  CMZ_FROM_CDI;
 double  CMZ_FROM_CFY;
 double  REF_AREA;
 double  REF_SPAN;
 double  REF_LEN_CMX;
 double  REF_LEN_CMY;
 double  REF_LEN_CMZ;
};

struct FLUEGEL
{
  struct PltDist *Distri;
};

struct EinFall
{
  struct FLUEGEL *Fluegel;
};

struct PLTdistribution
{
  int AnzFluegel;
  int AnzFaelle;
  struct EinFall *Fall;
};

void run_Lili(struct inputfile *ParaFile);
void read14 (struct file14 *rfile, struct inputfile *ParaFile);
void write_load_dist(struct file14 *File14, struct inputfile *input);
void readDistributionPLT (struct PLTdistribution *PLTdist, struct inputfile *input);
void run_PROZESS (struct inputfile *input);

#endif
