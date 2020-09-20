/******************run_Polint.h********************                  
*                                                *      
*      Vorentwurfs-Generator "VoGen"             *          
*      Version 1.00                              *          
*      Date of last modification 01.10.2009      *          
*                                                *          
*      Writen by: Jan Himisch                    *      
*             jan.himisch@dlr.de                 *      
*                                                *      
*                                                *
*************************************************/
#ifndef _RUNPOLINT_H_
#define _RUNPOLINT_H_


#include "input.h"

struct PolintDistribution
{
  double *X;
  double *Y;
  double *Z;
  double *S;
  double *LOCAL_CIRCULATION;
  double *CFX;
  double *CFX_FROM_CDI;
  double *CFX_FROM_AIRFOIL;
  double *CFY;
  double *CFZ;
  double *CFNORMAL;
  double *CMX;
  double *CMX_FROM_CFY;
  double *CMX_FROM_CFZ;
  double *CMY;
  double *CMY_FROM_CDI;
  double *CMY_FROM_AIRFOIL;
  double *CMZ;
  double *CMZ_FROM_CDI;
  double *CMZ_FROM_AIRFOIL;
  double *REF_AREA;
  double *REF_SPAN;
  double *REF_LEN_CMX;
  double *REF_LEN_CMY;
  double *REF_LEN_CMZ;
};

struct PolintTotal
{
  double ANGLE_OF_ATTACK;
  double YAW_ANGLE;
  double MACH_NO;
  double CFX;
  double CFX_FROM_CDI;
  double CFX_FROM_AIRFOIL;
  double CFY;
  double CFZ;
  double CMX;
  double CMX_FROM_CFY;
  double CMX_FROM_CFZ;
  double CMY;
  double CMY_FROM_CDI;
  double CMY_FROM_AIRFOIL;
  double CMZ;
  double CMZ_FROM_CDI;
  double CMZ_FROM_AIRFOIL;
  double REF_AREA;
  double REF_SPAN;
  double REF_LEN_CMX;
  double REF_LEN_CMY;
  double REF_LEN_CMZ;
  double MOM_REF_X;
  double MOM_REF_Y;
  double MOM_REF_Z;
};


struct POLINTOUT
{
  int NrOfPanel;
  int NrFaelle;
  struct PolintTotal *total;
  struct PolintDistribution *distribution;
};


void run_Polint(struct inputfile *ParaFile);
void read_Polint_Out (struct inputfile *ParaFile, struct POLINTOUT *PolintOut);

#endif
