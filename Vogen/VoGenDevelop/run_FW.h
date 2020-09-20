/******************run_FW.h***********************
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
#ifndef _RUNFW_H_
#define _RUNFW_H_

#include "input.h"

struct ProfileSortier
{
    char Airfoil_Str[50];
    int Fluegel;
    int section;
};

void run_FWProzess (struct inputfile *input);

#endif
