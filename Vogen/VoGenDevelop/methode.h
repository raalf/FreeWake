/******************methode.h**********************                  
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
#ifndef _METHODE_H_
#define _METHODE_H_

#include"run_LILI.h"
#include"input.h"

struct METHODE
{
  double dcw;
  double dcwAlfa;
  double dcwAerea;
};

void methodeAlt(struct inputfile *input, struct METHODE *method, struct file14 *resultfile);


#endif
