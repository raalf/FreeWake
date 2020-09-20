/******************run_Polint.h*******************                  
*                                                *      
*      Vorentwurfs-Generator "VoGen"             *          
*      Version 1.15                              *          
*      Date of last modification 16.11.2017      *          
*                                                *          
*      Writen by: Jan Himisch                    *      
*             jan.himisch@dlr.de                 *      
*                                                *      
*                                                *
*************************************************/
#ifndef _RUNPOLAR_H_
#define _RUNPOLAR_H_


#include <stdlib.h>
#include <stdio.h>
#include "errmes.h"
#include "input.h"
#include "create_LiLiInp.h"
#include "run_LILI.h"
#include "run_Polint.h"
#include "integration.h"
#include "methode.h"
#include "make_bsurf.h"
#include "run_STRUCT.h"
#include "trimming.h"
#include "run_Parametric.h"


void run_einzelnt(struct inputfile *ParaFile);

#endif
