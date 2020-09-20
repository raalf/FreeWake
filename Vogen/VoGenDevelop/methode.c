/******************methode.c**********************                  
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

#include <stdlib.h>
#include<math.h>
#include"methode.h"


void methodeAlt(struct inputfile *input, struct METHODE *method, struct file14 *resultfile)
{
  int lauf, lauf2;
  for (lauf=0; lauf< (input->numberCL+input->numberAlpha); lauf++)
  {
    method[lauf].dcwAerea=(resultfile->ABGEWICKELTE_FLAECHE- input->UrArea)/resultfile->ABGEWICKELTE_FLAECHE*input->BCw0;    
    method[lauf].dcwAlfa=0;
    for (lauf2=1; lauf2<input->kHOrd;lauf2++)
    {
      method[lauf].dcwAlfa=method[lauf].dcwAlfa+pow((input->BAlafa-resultfile->Fall[lauf].Alfa),lauf2)*input->kH[lauf2];
    }
    method[lauf].dcwAlfa=method[lauf].dcwAlfa*(-1);
    method[lauf].dcw=method[lauf].dcwAerea+method[lauf].dcwAlfa;
  }
}
