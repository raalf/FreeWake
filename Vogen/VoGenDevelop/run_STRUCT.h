/********************run_STRUCT.h*****************                  
*                                                *      
*      Vorentwurfs-Generator "VoGen"             *          
*      Version 2.00                              *
*      Date of last modification 23.04.2020      *
*                                                *          
*      Writen by: Jan Himisch                    *      
*             jan.himisch@dlr.de                 *      
*                                                *      
*                                                *
*************************************************/
#ifndef _RUNSTRUCT_H_
#define _RUNSTRUCT_H_

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#include <unistd.h>
#include "errmes.h"
#include "run_LILI.h"
#include "input.h"
#include "create_LiLiInp.h"
#include "integration.h"
#include "POLAR_INTER.h"

struct PunktLoesung
{
  double x;
  double y;
  double z;
  double dx;
  double dy;
  double dz;
  double phix;
  double phiy;
  double phiz;  
};

struct STRUCTLOESUNG
{
  int AnzDefoPunkte;
  //struct PunktLoesung *Punkt;
  struct PunktLoesung Punkt[60];
};

void run_stru(struct inputfile *input);

#endif
