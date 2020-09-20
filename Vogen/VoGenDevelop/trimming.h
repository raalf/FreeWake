/********************trimming.h*******************                 
*                                                *      
*      Vorentwurfs-Generator "VoGen"             *          
*      Version 2.00                              *
*      Date of last modification 22.06.2020      *
*                                                *          
*      Writen by: Jan Himisch                    *      
*             jan.himisch@dlr.de                 *      
*                                                *      
*                                                *
*************************************************/
#ifndef _TRIMMING_H_
#define _TRIMMING_H_

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
#include "run_polar.h"
#include "POLAR_INTER.h"
#include "mischer.h"

void checkTrimm(struct inputfile *input);
void Trimm_Schleife(struct inputfile *input);

#endif
