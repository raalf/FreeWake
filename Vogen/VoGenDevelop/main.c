/********************main.c***********************                  
*                                                *      
*      Vorentwurfs-Generator "VoGen"             *          
*      Version 2.00                              *
*      Date of last modification 21.04.2020      *
*                                                *          
*      Writen by: Jan Himisch                    *      
*             jan.himisch@dlr.de                 *      
*                                                *      
*                                                *
*************************************************/
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
#include "trimming.h"
#include "run_Parametric.h"

void printHEADER(void)
{
  printf("\n");
  printf("  **************************************************************\n");
  printf("  *              _    __        ______                         *\n");
  printf("  *             | |  / /____   / ____/___   ____               *\n");
  printf("  *             | | / // __ \134 / / __ / _ \134 / __ \134              *\n");
  printf("  *             | |/ // /_/ // /_/ //  __// / / /              *\n");
  printf("  *             |___/ \134____/ \134____/ \134___//_/ /_/               *\n");
  printf("  *                                                            *\n");
  printf("  *                        Version 1.98  (2.0 alpha)           *\n");
  printf("  *                                                            *\n");
  printf("  *     A program for fast generation, parameterisation and    *\n");
  printf("  *       post-processing of Lifting-Line-Input and Output     *\n");
  printf("  *                                                            *\n");
  printf("  *  Developed by: Jan Himisch   jan.himisch@dlr.de            *\n");
  printf("  *                                                            *\n");
  printf("  **************************************************************\n");
  printf("\n");
}

void printBOT(void)
{
  printf("\n");
  printf("  VoGen run completed!!!\n");
  printf("  **************************************************************\n");
  printf("\n");
}

int main (zahl,dateiname)
int zahl;
char *dateiname[];
{
  struct inputfile input;  

  if (zahl != 2)
  {
    errmessage (11);
  }
  printHEADER();
  printf("  Read Input ... ");
  read_input(dateiname[1], &input); 
  printf("finished\n\n"); 

  if (input.AnzahlMischer > 0)
  {
     run_mischer(&input);
  }
  if ((input.runOneByOne == 1) &&((input.numberAlpha + input.numberCL)>1))
  { 
    if (input.OldInput == 0)
    {
        errmessage(69);
    }
    run_einzelnt(&input);
  }else if (input.NumberOfParametSetsinput>0)
  {
    run_parametric(&input, 0);
  }else if (input.Austrimmen > 1){
    Trimm_Schleife(&input);  
  }else if (input.Struckt_Coup_ON == 1){    
    run_stru(&input);
  }else{          
      if (input.LiliStart >0){
          run_PROZESS(&input);
      }else if (input.FwStart >0) {
          run_FWProzess(&input);
      }
  }
  if (input.makeBsurfinp == 1)
  {
    make_bsurf_inp(&input);
  }
  printBOT();  
  return 0;
}
