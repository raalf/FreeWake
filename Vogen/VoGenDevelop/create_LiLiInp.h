/***************create_LiLiInp.h******************                  
*                                                *      
*      Vorentwurfs-Generator "VoGen"             *          
*      Version 1.07                              *          
*      Date of last modification 22.02.2011      *          
*                                                *          
*      Writen by: Jan Himisch                    *      
*             jan.himisch@dlr.de                 *      
*                                                *      
*                                                *
*************************************************/
#ifndef _CREATE_LILIINP_H_
#define _CREATE_LILIINP_H_

#include "input.h"
#include "errmes.h"
#include <math.h>




struct SectionDistribution 
{
  char name[MAX_LENGTH_PNAME];
  double *etaPos;
  int NrEtaPos;
};

void create_LiLiInp_File (struct inputfile *ParaFile);
struct SectionDistribution *getSecInfo(struct Einfluegel *Fluegel, int *NrSect);     

#endif
