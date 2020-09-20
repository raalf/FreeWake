/*****************make_bsurf.h********************                  
*                                                *      
*      Vorentwurfs-Generator "VoGen"             *          
*      Version 1.09                              *          
*      Date of last modification 27.05.2013      *          
*                                                *          
*      Writen by: Jan Himisch                    *      
*             jan.himisch@dlr.de                 *      
*                                                *      
*                                                *
*************************************************/
#ifndef _MAKE_BSURF_H_
#define _MAKE_BSURF_H_

#include "run_LILI.h"
#include "input.h"

struct BSURF_DATEN
{
  double Xpos;
  double Ypos;
  double Zpos;
  double CL;
  double CD;
  double CS;
  double CN;   
};

void make_bsurf_inp( struct inputfile *input);
void read_BSURF( struct inputfile *input,  struct PLTdistribution *PLTdist);
void Comp_Square_Diff(struct inputfile *input, struct PLTdistribution *PLTdist);
#endif
