/******************run_Polint.c*******************                  
*                                                *      
*      Vorentwurfs-Generator "VoGen"             *          
*      for Version 1.05                          *          
*      Date of last modification 25.08.2010      *          
*                                                *          
*      Writen by: Jan Himisch                    *      
*             jan.himisch@dlr.de                 *      
*                                                *      
*                                                *
*************************************************/

#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include "run_Polint.h"
#include "errmes.h"

void run_Polint(struct inputfile *ParaFile)
{
  char EXECUTRstring[450] = {"\0"};
  int return_value;
  char *ptr; 
  char dummyLName[50];
  
  strcpy(dummyLName, ParaFile->LAUFNAME);
  ptr = strrchr(dummyLName, '.');
  *ptr = '\0';
  strcpy(EXECUTRstring, ParaFile->POLINTEXE);
  strcat(EXECUTRstring, " ");
  strcat(EXECUTRstring, dummyLName);
  strcat(EXECUTRstring, ".lili.");
  strcat(EXECUTRstring,ParaFile->LiLiVersionName); 
  strcat(EXECUTRstring, "/export/");
  strcat(EXECUTRstring, dummyLName);
  strcat(EXECUTRstring, ".xml ");  
  strcat(EXECUTRstring, ParaFile->PROFIPFAD);
  
  return_value = system(EXECUTRstring);
  if (return_value != 0)
  {
    errmessage (25);
  }    
}

void read_Polint_Out(struct inputfile *ParaFile, struct POLINTOUT *PolintOut)
{
  FILE *fopen(),*fresult;
  char ZEILE_READ[MaxZeilenLaenge]; 
  char Resultfile[150] = {"\0"};
  char dummyLName[50];
  char *ptr;
  int lauf,lauf2; 
  char DATNAME[50],EXECNAME[100], dummy[4]; 

  strcpy(dummyLName, ParaFile->LAUFNAME);
  ptr = strrchr(dummyLName, '.');
  *ptr = '\0';
  strcpy(Resultfile, dummyLName);
  strcat(Resultfile, "_total_polint.plt");
  PolintOut->NrFaelle=ParaFile->numberAlpha+ParaFile->numberCL;
  PolintOut->total=(struct PolintTotal *)malloc(PolintOut->NrFaelle*sizeof(struct PolintTotal));
  fresult=fopen(Resultfile,"r");
  fgets (ZEILE_READ,MaxZeilenLaenge,fresult);
  do 
  {
    if (strstr (ZEILE_READ, "ZONE T = \"Configuration\"")!=NULL)
    {
      for (lauf=0; lauf<PolintOut->NrFaelle; lauf++)
      {
	fscanf(fresult," %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n", &PolintOut->total[lauf].ANGLE_OF_ATTACK, &PolintOut->total[lauf].YAW_ANGLE, &PolintOut->total[lauf].MACH_NO, &PolintOut->total[lauf].CFX ,\
	 &PolintOut->total[lauf].CFX_FROM_CDI , &PolintOut->total[lauf].CFX_FROM_AIRFOIL , &PolintOut->total[lauf].CFY , &PolintOut->total[lauf].CFZ , &PolintOut->total[lauf].CMX , &PolintOut->total[lauf].CMX_FROM_CFY , &PolintOut->total[lauf].CMX_FROM_CFZ , &PolintOut->total[lauf].CMY , \
	 &PolintOut->total[lauf].CMY_FROM_CDI , &PolintOut->total[lauf].CMY_FROM_AIRFOIL , &PolintOut->total[lauf].CMZ , &PolintOut->total[lauf].CMZ_FROM_CDI , &PolintOut->total[lauf].CMZ_FROM_AIRFOIL , &PolintOut->total[lauf].REF_AREA , &PolintOut->total[lauf].REF_SPAN ,\
	 &PolintOut->total[lauf].REF_LEN_CMX , &PolintOut->total[lauf].REF_LEN_CMY , &PolintOut->total[lauf].REF_LEN_CMZ , &PolintOut->total[lauf].MOM_REF_X , &PolintOut->total[lauf].MOM_REF_Y , &PolintOut->total[lauf].MOM_REF_Z  );       
      }
    }
    fgets (ZEILE_READ,MaxZeilenLaenge,fresult); 
  }
  while (!feof(fresult));
  fclose(fresult);
  strcpy(Resultfile, dummyLName);
  strcat(Resultfile, "_distribution_polint.plt\0");
  fresult=fopen(Resultfile,"r");
  PolintOut->distribution=(struct PolintDistribution *)malloc(PolintOut->NrFaelle*sizeof(struct PolintDistribution));
  for (lauf=0; lauf < PolintOut->NrFaelle; lauf++)
  {
    PolintOut->distribution[lauf].X=(double *)malloc(ParaFile->AnzahlSpanPanel*sizeof(double));
    PolintOut->distribution[lauf].Y=(double *)malloc(ParaFile->AnzahlSpanPanel*sizeof(double));
    PolintOut->distribution[lauf].Z=(double *)malloc(ParaFile->AnzahlSpanPanel*sizeof(double));
    PolintOut->distribution[lauf].S=(double *)malloc(ParaFile->AnzahlSpanPanel*sizeof(double));
    PolintOut->distribution[lauf].LOCAL_CIRCULATION=(double *)malloc(ParaFile->AnzahlSpanPanel*sizeof(double));
    PolintOut->distribution[lauf].CFX=(double *)malloc(ParaFile->AnzahlSpanPanel*sizeof(double));
    PolintOut->distribution[lauf].CFX_FROM_CDI=(double *)malloc(ParaFile->AnzahlSpanPanel*sizeof(double));
    PolintOut->distribution[lauf].CFX_FROM_AIRFOIL=(double *)malloc(ParaFile->AnzahlSpanPanel*sizeof(double));
    PolintOut->distribution[lauf].CFY=(double *)malloc(ParaFile->AnzahlSpanPanel*sizeof(double));
    PolintOut->distribution[lauf].CFZ=(double *)malloc(ParaFile->AnzahlSpanPanel*sizeof(double));
    PolintOut->distribution[lauf].CFNORMAL=(double *)malloc(ParaFile->AnzahlSpanPanel*sizeof(double));
    PolintOut->distribution[lauf].CMX=(double *)malloc(ParaFile->AnzahlSpanPanel*sizeof(double));
    PolintOut->distribution[lauf].CMX_FROM_CFY=(double *)malloc(ParaFile->AnzahlSpanPanel*sizeof(double));
    PolintOut->distribution[lauf].CMX_FROM_CFZ=(double *)malloc(ParaFile->AnzahlSpanPanel*sizeof(double));
    PolintOut->distribution[lauf].CMY=(double *)malloc(ParaFile->AnzahlSpanPanel*sizeof(double));
    PolintOut->distribution[lauf].CMY_FROM_CDI=(double *)malloc(ParaFile->AnzahlSpanPanel*sizeof(double));
    PolintOut->distribution[lauf].CMY_FROM_AIRFOIL=(double *)malloc(ParaFile->AnzahlSpanPanel*sizeof(double));
    PolintOut->distribution[lauf].CMZ=(double *)malloc(ParaFile->AnzahlSpanPanel*sizeof(double));
    PolintOut->distribution[lauf].CMZ_FROM_CDI =(double *)malloc(ParaFile->AnzahlSpanPanel*sizeof(double));
    PolintOut->distribution[lauf].CMZ_FROM_AIRFOIL=(double *)malloc(ParaFile->AnzahlSpanPanel*sizeof(double));
    PolintOut->distribution[lauf].REF_AREA=(double *)malloc(ParaFile->AnzahlSpanPanel*sizeof(double));
    PolintOut->distribution[lauf].REF_SPAN=(double *)malloc(ParaFile->AnzahlSpanPanel*sizeof(double));
    PolintOut->distribution[lauf].REF_LEN_CMX=(double *)malloc(ParaFile->AnzahlSpanPanel*sizeof(double));
    PolintOut->distribution[lauf].REF_LEN_CMY=(double *)malloc(ParaFile->AnzahlSpanPanel*sizeof(double));
    PolintOut->distribution[lauf].REF_LEN_CMZ=(double *)malloc(ParaFile->AnzahlSpanPanel*sizeof(double));
  }
  do
  {
    fgets (ZEILE_READ,MaxZeilenLaenge,fresult);
  }
  while (strstr (ZEILE_READ, "ZONE")==NULL);
  //printf("\t\tHALLO\n");
  for (lauf=0; lauf < PolintOut->NrFaelle; lauf++)
  {
    for (lauf2=0; lauf2 < ParaFile->AnzahlSpanPanel; lauf2++)
    {
      fgets (ZEILE_READ,MaxZeilenLaenge,fresult);
      sscanf(ZEILE_READ," %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf", &PolintOut->distribution[lauf].X[lauf2] ,&PolintOut->distribution[lauf].Y[lauf2] ,&PolintOut->distribution[lauf].Z[lauf2] ,&PolintOut->distribution[lauf].S[lauf2] ,\
       &PolintOut->distribution[lauf].LOCAL_CIRCULATION[lauf2] , &PolintOut->distribution[lauf].CFX[lauf2], &PolintOut->distribution[lauf].CFX_FROM_CDI[lauf2] , &PolintOut->distribution[lauf].CFX_FROM_AIRFOIL[lauf2] , &PolintOut->distribution[lauf].CFY[lauf2] , \
       &PolintOut->distribution[lauf].CFZ[lauf2] , &PolintOut->distribution[lauf].CFNORMAL[lauf2] , &PolintOut->distribution[lauf].CMX[lauf2] , &PolintOut->distribution[lauf].CMX_FROM_CFY[lauf2] , &PolintOut->distribution[lauf].CMX_FROM_CFZ[lauf2] , \
       &PolintOut->distribution[lauf].CMY[lauf2] , &PolintOut->distribution[lauf].CMY_FROM_CDI[lauf2] , &PolintOut->distribution[lauf].CMY_FROM_AIRFOIL[lauf2] , &PolintOut->distribution[lauf].CMZ[lauf2] , &PolintOut->distribution[lauf].CMZ_FROM_CDI[lauf2] , \
       &PolintOut->distribution[lauf].CMZ_FROM_AIRFOIL[lauf2] , &PolintOut->distribution[lauf].REF_AREA[lauf2] , &PolintOut->distribution[lauf].REF_SPAN[lauf2] , &PolintOut->distribution[lauf].REF_LEN_CMX[lauf2] , &PolintOut->distribution[lauf].REF_LEN_CMY[lauf2] ,\
       &PolintOut->distribution[lauf].REF_LEN_CMZ[lauf2] );           
    }
    fgets (ZEILE_READ,MaxZeilenLaenge,fresult);
    fgets (ZEILE_READ,MaxZeilenLaenge,fresult);
  } 
  fclose (fresult);
  strcpy(DATNAME, dummyLName);
  strcat(DATNAME,"_polars_polint.plt"); 
  strcpy(EXECNAME,"mv ");
  strcat(EXECNAME,DATNAME);
  strcat(EXECNAME," LaufV1.lili.");
  strcat(EXECNAME,ParaFile->LiLiVersionName);  
  system(EXECNAME); 
  strcpy(DATNAME, dummyLName);
  strcat(DATNAME,"_total_polint.plt"); 
  strcpy(EXECNAME,"mv ");
  strcat(EXECNAME,DATNAME);
  strcat(EXECNAME," LaufV1.lili.");
  strcat(EXECNAME,ParaFile->LiLiVersionName);  
  system(EXECNAME); 
  strcpy(DATNAME, dummyLName);
  strcat(DATNAME,"_distribution_polint.plt"); 
  strcpy(EXECNAME,"mv ");
  strcat(EXECNAME,DATNAME);
  strcat(EXECNAME," LaufV1.lili.");
  strcat(EXECNAME,ParaFile->LiLiVersionName);  
  system(EXECNAME); 
  strcpy(DATNAME, dummyLName);
  strcat(DATNAME,"_polint.xml"); 
  strcpy(EXECNAME,"mv ");
  strcat(EXECNAME,DATNAME);
  strcat(EXECNAME," LaufV1.lili.");
  strcat(EXECNAME,ParaFile->LiLiVersionName);  
  system(EXECNAME); 
     
}
