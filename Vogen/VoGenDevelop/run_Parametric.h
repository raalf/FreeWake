/****************run_Parametric.h*****************                  
*                                                *      
*      Vorentwurfs-Generator "VoGen"             *          
*      Version 1.09                              *          
*      Date of last modification 10.04.2013      *          
*                                                *          
*      Writen by: Jan Himisch                    *      
*             jan.himisch@dlr.de                 *      
*                                                *      
*                                                *
*************************************************/
#ifndef _RUNPARAMETRIC_H_
#define _RUNPARAMETRIC_H_

#include "input.h"
#include "create_LiLiInp.h"
#include "run_Polint.h"
#include "errmes.h"
#include "run_LILI.h"

#define deltaS 0.01

struct OutPutFileName
{
  char name[100];
};

void run_parametric (struct inputfile  *input, int laufSet);


#endif
