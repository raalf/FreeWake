/****************integration.h********************                  
*                                                *      
*      Vorentwurfs-Generator "VoGen"             *          
*      Version 1.03                              *          
*      Date of last modification 04.06.2010      *          
*                                                *          
*      Writen by: Jan Himisch                    *      
*             jan.himisch@dlr.de                 *      
*                                                *      
*                                                *
*************************************************/
#ifndef _INTEGRATION_H_
#define _INTEGRATION_H_


struct KreuzList
{
  int KrNr;
  int AusgangsFluegel;
};

struct DISTRIBUTION
{
  double X;
  double Y;
  double Z;
  double L;
  double breite;
  double Thickness;
  double CN;
  double CBM; //CMX
  double CBT; //CMY
  double CBG; //CMZ
  double FX;
  double FY;
  double FZ;
};

struct FLUEGELint 
{
  struct DISTRIBUTION  *distri;
  double IBM;
  double IBMD;
};

struct BM_FALL
{
  double WRBM;
  double IBM;
  double IBMD;
  struct FLUEGELint *Fluegel;
};


void BM_integration (struct inputfile *input, struct file14  *file14, struct PLTdistribution *PLTdist, struct BM_FALL *BMcase, int Bsurf);
void write_BM_distribution(struct BM_FALL *BMcase, struct file14  *File14, struct inputfile *input, int schalt);

#endif
