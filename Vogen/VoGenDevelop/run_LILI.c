/******************run_LILI.c********************                  
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
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "run_LILI.h"
#include "errmes.h"
#include "run_Polint.h"
#include "integration.h"
#include "methode.h"
#include "create_LiLiInp.h"
#include "POLAR_INTER.h"
#include "make_bsurf.h"


void run_Lili(struct inputfile *ParaFile)
{
  char EXECUTRstring[150] = {"\0"};
  int return_value;
  
  strcpy(EXECUTRstring, ParaFile->LILIEXE);
  strcat(EXECUTRstring, " -pj:");
  strcat(EXECUTRstring, ParaFile->LAUFNAME);
 
  if (ParaFile->XMLEA == 0)
  {
    strcat(EXECUTRstring, " -xo:0");
  }else{
    strcat(EXECUTRstring, " -xo:1");
  }
  printf("  Run Lifting-Line ... \n");
  return_value = system(EXECUTRstring); 
  if (return_value != 0)
  {
    errmessage (24);
  }
  printf("                   ...finished\n\n"); 
}  

void run_PROZESS(struct inputfile *input)
{
    struct PLTdistribution PLTDIST;
    struct file14 resultfile;
    struct POLINTOUT PolintOut;
    struct BM_FALL *BMDist;
    struct METHODE *method;
    struct GESAMT_POLARint *GPolare;



    if (input->OldInput == 0)
    {
      create_LiLiInp_File (input);
    }
    if (input->LiliStart==1)
    {
      run_Lili(input);
      read14 (&resultfile , input);
      readDistributionPLT (&PLTDIST, input);
      if (input->PolintStart==1)
      {
        run_Polint(input);
        read_Polint_Out (input, &PolintOut);
      }
      BMDist=(struct BM_FALL *)malloc((input->numberCL+input->numberAlpha)*sizeof(struct BM_FALL));
      BM_integration (input, &resultfile, &PLTDIST, BMDist, 0);
      write_BM_distribution(BMDist, &resultfile, input, 0);
      method=(struct METHODE *)malloc((input->numberCL+input->numberAlpha)*sizeof(struct METHODE));
      methodeAlt(input, method, &resultfile);
      write_load_dist(&resultfile, input);
      if (input->InterneInterpolation == 1)
      {
          GPolare=(struct GESAMT_POLARint *)malloc((input->numberCL+input->numberAlpha)*sizeof(struct GESAMT_POLARint));
          PolarInterpolation(input,&PLTDIST,&resultfile,GPolare);
      }
      if (input->TAUbsurfINP == 1)
      {
        read_BSURF(input, &PLTDIST);
        BMDist=(struct BM_FALL *)realloc(BMDist,(input->numberCL+input->numberAlpha+1)*sizeof(struct BM_FALL));
        BM_integration (input, &resultfile, &PLTDIST, BMDist, 1);
        write_BM_distribution(BMDist, &resultfile, input, 1);
        if (input->BsurfSquareDiff == 1)
        {
          Comp_Square_Diff(input,&PLTDIST);
        }
      }
    }
}


void initfile14 (struct file14 *rfile)
{
  rfile->LREF_CL=-1;
  rfile->LREF_CM=-1;
  rfile->LREF_CN=-1;
  rfile->BEZUGSFLAECHE=-1; 
  rfile->GRUNDRISSFLAECHE=-1;
  rfile->ABGEWICKELTE_FLAECHE=-1; 
  rfile->BEZUGSSPANNWEITE=-1;
  rfile->PROJIZIERTE_SPANNWEITE=-1; 
  rfile->BEZUGSSTRECKUNG=-1;
}

void cstring(char string1[], char string2[], int start, int laenge)
{
  int lauf;
  for (lauf = start; lauf < start+laenge ; lauf ++)
  {
    string1[lauf -start] = string2[lauf]; 
  }
  string1[lauf -start] = '\0';
}

void readDistributionPLT (struct PLTdistribution *PLTdist, struct inputfile *input)
{
  int Pointpos,lauf,lauf2, laufWing;
  FILE *fresult;
  char ZEILE_READ[MaxZeilenLaenge]; 
  char Resultfile[150] = {"\0"};
  
  Pointpos=strcspn(input->LAUFNAME,".");
  strncpy(Resultfile,input->LAUFNAME,Pointpos);
  strcat(Resultfile, ".lili.");
  strcat(Resultfile, input->LiLiVersionName);
  strcat(Resultfile, "/export/tecplot/");
  strncat(Resultfile,input->LAUFNAME,Pointpos); 
  strcat(Resultfile, "_distribution.plt");
  fresult=fopen(Resultfile,"r");
  PLTdist->AnzFluegel=input->NbWings;
  PLTdist->AnzFaelle=input->numberCL+input->numberAlpha;
  if (fresult==NULL)
  {
    errmessage(32);
  }
  PLTdist->Fall=(struct EinFall *)malloc(PLTdist->AnzFaelle*sizeof(struct EinFall));
  for (lauf=0; lauf < PLTdist->AnzFaelle; lauf++)
  {
    PLTdist->Fall[lauf].Fluegel=(struct FLUEGEL *)malloc(input->NbWings*sizeof(struct FLUEGEL));
  }
  fgets (ZEILE_READ,MaxZeilenLaenge,fresult);
  lauf=0;
  laufWing=1;
  do 
  {
     if ((strstr (ZEILE_READ, "ZONE T=")!=NULL) && (strstr (ZEILE_READ, "Wing")!=NULL))
     {
       PLTdist->Fall[lauf].Fluegel[laufWing-1].Distri=(struct PltDist *)malloc(input->Fluegel[laufWing-1].AnzahlSpanPanel*sizeof(struct PltDist));
       for (lauf2=0; lauf2 < input->Fluegel[laufWing-1].AnzahlSpanPanel; lauf2++)
       {
         fscanf(fresult," %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n",&PLTdist->Fall[lauf].Fluegel[laufWing-1].Distri[lauf2].X,&PLTdist->Fall[lauf].Fluegel[laufWing-1].Distri[lauf2].Y,&PLTdist->Fall[lauf].Fluegel[laufWing-1].Distri[lauf2].Z,&PLTdist->Fall[lauf].Fluegel[laufWing-1].Distri[lauf2].S
         ,&PLTdist->Fall[lauf].Fluegel[laufWing-1].Distri[lauf2].LOCAL_CIRCULATION,&PLTdist->Fall[lauf].Fluegel[laufWing-1].Distri[lauf2].CFX_FROM_CDI,&PLTdist->Fall[lauf].Fluegel[laufWing-1].Distri[lauf2].CFY,&PLTdist->Fall[lauf].Fluegel[laufWing-1].Distri[lauf2].CFZ,&PLTdist->Fall[lauf].Fluegel[laufWing-1].Distri[lauf2].CFNORMAL,&PLTdist->Fall[lauf].Fluegel[laufWing-1].Distri[lauf2].CMX
         ,&PLTdist->Fall[lauf].Fluegel[laufWing-1].Distri[lauf2].CMX_FROM_CFY,&PLTdist->Fall[lauf].Fluegel[laufWing-1].Distri[lauf2].CMX_FROM_CFZ,&PLTdist->Fall[lauf].Fluegel[laufWing-1].Distri[lauf2].CMY,&PLTdist->Fall[lauf].Fluegel[laufWing-1].Distri[lauf2].CMY_FROM_CDI,&PLTdist->Fall[lauf].Fluegel[laufWing-1].Distri[lauf2].CMY_FROM_CFZ
         ,&PLTdist->Fall[lauf].Fluegel[laufWing-1].Distri[lauf2].CMZ,&PLTdist->Fall[lauf].Fluegel[laufWing-1].Distri[lauf2].CMZ_FROM_CDI,&PLTdist->Fall[lauf].Fluegel[laufWing-1].Distri[lauf2].CMZ_FROM_CFY,&PLTdist->Fall[lauf].Fluegel[laufWing-1].Distri[lauf2].REF_AREA,&PLTdist->Fall[lauf].Fluegel[laufWing-1].Distri[lauf2].REF_SPAN,&PLTdist->Fall[lauf].Fluegel[laufWing-1].Distri[lauf2].REF_LEN_CMX
         ,&PLTdist->Fall[lauf].Fluegel[laufWing-1].Distri[lauf2].REF_LEN_CMY,&PLTdist->Fall[lauf].Fluegel[laufWing-1].Distri[lauf2].REF_LEN_CMZ);
       }
       laufWing++;
       if (lauf > (PLTdist->AnzFaelle*input->NbWings))
       {
         printf("Anzahl Faelle: %d, Lauf:%d\n",PLTdist->AnzFaelle, lauf);
         errmessage(0);
       }
       if (laufWing > input->NbWings)
       {
         laufWing=1;
         lauf++;
       }
     }
     fgets (ZEILE_READ,MaxZeilenLaenge,fresult);
  }                                                                                                                                 
  while (!feof(fresult));
  fclose (fresult);
  
  /*DEBUG
  printf("\n");
  laufWing=1;
  for (lauf=0; lauf < PLTdist->AnzFaelle*input->NbWings; lauf ++)
  {
    printf("\nFall: %d\n",lauf);
    for (lauf2=0; lauf2< input->Fluegel[laufWing-1].AnzahlSpanPanel; lauf2++)
    {       
       printf("B");
       printf(" %lf %lf \n",PLTdist->Fall[lauf].Distri[lauf2].LOCAL_CIRCULATION,PLTdist->Fall[lauf].Distri[lauf2].CMZ);
    }
    laufWing++;
    if (laufWing > input->NbWings)
    {
      laufWing=1;
    }
  }
  DEBUG-END*/                                                                                                         
}                                                                                                                                 
                                                                                                                                 
void read14 (struct file14 *rfile, struct inputfile *ParaFile)                                                                         
{                                                                                                                                 
  FILE *fopen(),*fresult;                                                                                                         
  char Resultfile[150] = {"\0"};                                                                                                 
  int Pointpos,lauf, FallNr, Iteration, check;
  char ZEILE_READ[MaxZeilenLaenge];                                                                                                  
  char dummy1[50],dummy2[50],dummy3[50],dummy4[50], dummy5[5], dummy6[5],dummy7[5];                                                 
  long DatPointer;                                                                                                                 
                                                                                                                                   
  Iteration = 0;                                                                                                                 
  initfile14(rfile);                                                                                                                 
  Pointpos=strcspn(ParaFile->LAUFNAME,".");                                                                                         
  strncpy(Resultfile,ParaFile->LAUFNAME,Pointpos);                                                                                                                                                                     
  strcat(Resultfile, ".lili.");
  strcat(Resultfile, ParaFile->LiLiVersionName);
  strcat(Resultfile, "/results/");  
  strncat(Resultfile,ParaFile->LAUFNAME,Pointpos);                                                                                 
  strcat(Resultfile,".14");
  //printf("\n StrOutput :  %s\n", Resultfile);
  fresult=fopen(Resultfile,"r");
  if (fresult==NULL)
  {
    errmessage(31);
  }
  ZR_NULL;
  fgets (ZEILE_READ,MaxZeilenLaenge,fresult);
  rfile->Fall=(struct FALL *)malloc((ParaFile->numberCL+ParaFile->numberAlpha)*sizeof(struct FALL));
  for (lauf =0; lauf < (ParaFile->numberCL+ParaFile->numberAlpha); lauf++)
  {
    rfile->Fall[lauf].beiwerte=(struct CIRC *)malloc(2*ParaFile->AnzahlPanel*sizeof(struct CIRC));
    rfile->Fall[lauf].FluegelBw=(struct FluegelBeiwerte *)malloc(ParaFile->NbWings*sizeof(struct FluegelBeiwerte));
  }
  do
  {
    if ((strstr (ZEILE_READ, "MOMENTENBEZUGSLAENGEN")!=NULL) && (rfile->LREF_CL==-1))
    {
      Pointpos=strcspn(ParaFile->LAUFNAME,"=");
      sscanf(ZEILE_READ," %s  %s  %s  %lf  %s  %s  %lf %s  %s  %lf ",&dummy1[0], &dummy2[0], &dummy5[0], &rfile->LREF_CL, &dummy3[0], &dummy6[0], &rfile->LREF_CM, &dummy4[0], &dummy7[0], &rfile->LREF_CN);      
      //printf("\n%lf %lf %lf\n",rfile->LREF_CL,rfile->LREF_CM,rfile->LREF_CN);
    }
    if ((strstr (ZEILE_READ, "ABGEWICKELTE FLAECHE")!=NULL) && (rfile->ABGEWICKELTE_FLAECHE==-1))
    {
      Pointpos=strcspn(ParaFile->LAUFNAME,"=");
      sscanf(ZEILE_READ," %s  %s %lf  %s  %s  %lf %s %s %s  %lf ",&dummy1[0], &dummy5[0], &rfile->BEZUGSFLAECHE, &dummy2[0], &dummy6[0], &rfile->GRUNDRISSFLAECHE, &dummy3[0], &dummy4[0], &dummy7[0], &rfile->ABGEWICKELTE_FLAECHE);      
      //printf("\n%lf %lf %lf\n",rfile->BEZUGSFLAECHE,rfile->GRUNDRISSFLAECHE,rfile->ABGEWICKELTE_FLAECHE);
    }
    if ((strstr (ZEILE_READ, "PROJIZIERTE SPANNWEITE")!=NULL) && (rfile->PROJIZIERTE_SPANNWEITE==-1))
    {
      Pointpos=strcspn(ParaFile->LAUFNAME,"=");
      sscanf(ZEILE_READ," %s  %s %lf  %s  %s %s %lf ",&dummy1[0], &dummy5[0], &rfile->BEZUGSSPANNWEITE, &dummy2[0], &dummy3[0], &dummy6[0], &rfile->PROJIZIERTE_SPANNWEITE);           
      //printf("\n%lf %lf\n",rfile->BEZUGSSPANNWEITE,rfile->PROJIZIERTE_SPANNWEITE);
    }
    if ((strstr (ZEILE_READ, "BEZUGSSTRECKUNG")!=NULL) && (rfile->BEZUGSSTRECKUNG==-1))
    {
      Pointpos=strcspn(ParaFile->LAUFNAME,"=");
      sscanf(ZEILE_READ," %s  %s  %lf ",&dummy1[0], &dummy5[0], &rfile->BEZUGSSTRECKUNG);          
      //printf("\n%lf \n",rfile->BEZUGSSTRECKUNG);
    }
    if (strstr (ZEILE_READ, "GEOMETRISCHE DATEN DER FLUEGELELEMENTE")!=NULL) 
    {
      nfgets(3);
      rfile->GEODAT=(struct GEO_OUT *)malloc(2*ParaFile->AnzahlPanel*sizeof(struct GEO_OUT));
      for (lauf =0; lauf < 2*ParaFile->AnzahlPanel; lauf++)
      {
        ZR_NULL;
        fgets (ZEILE_READ,MaxZeilenLaenge,fresult);
        sscanf(ZEILE_READ," %s %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf ", &dummy5[0], &rfile->GEODAT[lauf].X1, &rfile->GEODAT[lauf].Y1, &rfile->GEODAT[lauf].Z1, &rfile->GEODAT[lauf].L1, &rfile->GEODAT[lauf].X2, &rfile->GEODAT[lauf].Y2, &rfile->GEODAT[lauf].Z2, &rfile->GEODAT[lauf].L2, &rfile->GEODAT[lauf].XA, &rfile->GEODAT[lauf].YA, &rfile->GEODAT[lauf].ZA, &rfile->GEODAT[lauf].LA);
      }
    }
    
    if (strstr (ZEILE_READ, "Fall")!=NULL)
    {
      sscanf(ZEILE_READ," %s  %d ",&dummy1[0], &FallNr);
      if (FallNr <= ParaFile->numberAlpha)
      {
        do
        {
          ZR_NULL;
          fgets (ZEILE_READ,MaxZeilenLaenge,fresult);
        }
        while (strstr (ZEILE_READ, "OERTLICHE WERTE VON RESULTIERENDEM AUFTRIEB UND WIDERSTAND")==NULL);
        nfgets(3);
        for (lauf =0; lauf < ParaFile->AnzahlPanel; lauf++)
        {
          ZR_NULL;
          fgets (ZEILE_READ,MaxZeilenLaenge,fresult);
          check=sscanf(ZEILE_READ," %s %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %e", &dummy5[0], &rfile->Fall[FallNr-1].beiwerte[lauf].XA, &rfile->Fall[FallNr-1].beiwerte[lauf].YA, &rfile->Fall[FallNr-1].beiwerte[lauf].ZA, &rfile->Fall[FallNr-1].beiwerte[lauf].LA, &rfile->Fall[FallNr-1].beiwerte[lauf].CNORM01, &rfile->Fall[FallNr-1].beiwerte[lauf].CNORM0, &rfile->Fall[FallNr-1].beiwerte[lauf].CNORM02, &rfile->Fall[FallNr-1].beiwerte[lauf].DNORM01, &rfile->Fall[FallNr-1].beiwerte[lauf].DNORM0, &rfile->Fall[FallNr-1].beiwerte[lauf].DNORM02, &rfile->Fall[FallNr-1].beiwerte[lauf].CNORM1, &rfile->Fall[FallNr-1].beiwerte[lauf].CNORM, &rfile->Fall[FallNr-1].beiwerte[lauf].CNORM2, &rfile->Fall[FallNr-1].beiwerte[lauf].CWI);         
        }
        if (check == EOF)
        {
          errmessage(31);
        }
        nfgets(5);
        cstring(&dummy1[0], &ZEILE_READ[0], 72, 14);
        sscanf(dummy1, "%lf",  &rfile->Fall[FallNr-1].Alfa);
        cstring(&dummy1[0], &ZEILE_READ[0], 109, 10);
        sscanf(dummy1, "%lf",  &rfile->Fall[FallNr-1].Beta);
        nfgets(5);
        cstring(&dummy1[0], &ZEILE_READ[0], 71, 14);
        sscanf(dummy1, "%lf",  &rfile->Fall[FallNr-1].CA);
        ZR_NULL;
        fgets (ZEILE_READ,MaxZeilenLaenge,fresult);
        cstring(&dummy1[0], &ZEILE_READ[0], 71, 14);
        sscanf(dummy1, "%lf",  &rfile->Fall[FallNr-1].CQ);
        ZR_NULL;
        fgets (ZEILE_READ,MaxZeilenLaenge,fresult);
        cstring(&dummy1[0], &ZEILE_READ[0], 71, 14);
        sscanf(dummy1, "%lf",  &rfile->Fall[FallNr-1].CWI);
        ZR_NULL;
        fgets (ZEILE_READ,MaxZeilenLaenge,fresult);
        cstring(&dummy1[0], &ZEILE_READ[0], 71, 14);
        sscanf(dummy1, "%lf",  &rfile->Fall[FallNr-1].CL);
        ZR_NULL;
        fgets (ZEILE_READ,MaxZeilenLaenge,fresult);
        cstring(&dummy1[0], &ZEILE_READ[0], 71, 14);
        sscanf(dummy1, "%lf",  &rfile->Fall[FallNr-1].CM);
        ZR_NULL;
        fgets (ZEILE_READ,MaxZeilenLaenge,fresult);
        cstring(&dummy1[0], &ZEILE_READ[0], 71, 14);
        sscanf(dummy1, "%lf",  &rfile->Fall[FallNr-1].CN);
        nfgets(11);
        for (lauf =1; lauf <= ParaFile->NbWings; lauf++)
        {
           cstring(&dummy1[0], &ZEILE_READ[0], 71, 14);
           sscanf(dummy1, "%lf",  &rfile->Fall[FallNr-1].FluegelBw[lauf-1].CA);
           ZR_NULL;
           fgets (ZEILE_READ,MaxZeilenLaenge,fresult);
           cstring(&dummy1[0], &ZEILE_READ[0], 71, 14);
           sscanf(dummy1, "%lf",  &rfile->Fall[FallNr-1].FluegelBw[lauf-1].CQ);
           ZR_NULL;
           fgets (ZEILE_READ,MaxZeilenLaenge,fresult);
           cstring(&dummy1[0], &ZEILE_READ[0], 71, 14);
           sscanf(dummy1, "%lf",  &rfile->Fall[FallNr-1].FluegelBw[lauf-1].CWI);  
           ZR_NULL;
           fgets (ZEILE_READ,MaxZeilenLaenge,fresult);
           cstring(&dummy1[0], &ZEILE_READ[0], 71, 14);
           sscanf(dummy1, "%lf",  &rfile->Fall[FallNr-1].FluegelBw[lauf-1].CL);
           ZR_NULL;
           fgets (ZEILE_READ,MaxZeilenLaenge,fresult);
           cstring(&dummy1[0], &ZEILE_READ[0], 71, 14);
           sscanf(dummy1, "%lf",  &rfile->Fall[FallNr-1].FluegelBw[lauf-1].CM);
           ZR_NULL;
           fgets (ZEILE_READ,MaxZeilenLaenge,fresult);
           cstring(&dummy1[0], &ZEILE_READ[0], 71, 14);
           sscanf(dummy1, "%lf",  &rfile->Fall[FallNr-1].FluegelBw[lauf-1].CN);  
           nfgets(7);           
        }
      }else{
        DatPointer=ftell(fresult);
        Iteration=FallNr;
        do
        {
          ZR_NULL;
          fgets (ZEILE_READ,MaxZeilenLaenge,fresult);
          if (strstr (ZEILE_READ, "Fall") != NULL)
          {
            sscanf(ZEILE_READ," %s  %d ",&dummy1[0], &FallNr);
            if (FallNr==Iteration)
            {
              DatPointer=ftell(fresult);
            }else{             
                  break;
            }
          }
        }
        while ((!feof(fresult)) && (Iteration==FallNr));  
        if (!feof(fresult))
        {
          FallNr--;
        }
        fseek(fresult, DatPointer, SEEK_SET);
        do
        {
          ZR_NULL;
          fgets (ZEILE_READ,MaxZeilenLaenge,fresult);
        }
        while (strstr (ZEILE_READ, "OERTLICHE WERTE VON RESULTIERENDEM AUFTRIEB UND WIDERSTAND")==NULL);
        nfgets(3);
        for (lauf =0; lauf < ParaFile->AnzahlPanel; lauf++)
        {
          ZR_NULL;
          fgets (ZEILE_READ,MaxZeilenLaenge,fresult);
          check=sscanf(ZEILE_READ," %s %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %e", &dummy5[0], &rfile->Fall[FallNr-1].beiwerte[lauf].XA, &rfile->Fall[FallNr-1].beiwerte[lauf].YA, &rfile->Fall[FallNr-1].beiwerte[lauf].ZA, &rfile->Fall[FallNr-1].beiwerte[lauf].LA, &rfile->Fall[FallNr-1].beiwerte[lauf].CNORM01, &rfile->Fall[FallNr-1].beiwerte[lauf].CNORM0, &rfile->Fall[FallNr-1].beiwerte[lauf].CNORM02, &rfile->Fall[FallNr-1].beiwerte[lauf].DNORM01, &rfile->Fall[FallNr-1].beiwerte[lauf].DNORM0, &rfile->Fall[FallNr-1].beiwerte[lauf].DNORM02, &rfile->Fall[FallNr-1].beiwerte[lauf].CNORM1, &rfile->Fall[FallNr-1].beiwerte[lauf].CNORM, &rfile->Fall[FallNr-1].beiwerte[lauf].CNORM2, &rfile->Fall[FallNr-1].beiwerte[lauf].CWI);         
        }
        if (check == EOF)
        {
          errmessage(31);
        }
        nfgets(5);
        cstring(&dummy1[0], &ZEILE_READ[0], 71, 14);
        sscanf(dummy1, "%lf",  &rfile->Fall[FallNr-1].Alfa);
        cstring(&dummy1[0], &ZEILE_READ[0], 109, 10);
        sscanf(dummy1, "%lf",  &rfile->Fall[FallNr-1].Beta);
        nfgets(6);
        cstring(&dummy1[0], &ZEILE_READ[0], 71, 14);
        sscanf(dummy1, "%lf",  &rfile->Fall[FallNr-1].CA);
        ZR_NULL;
        fgets (ZEILE_READ,MaxZeilenLaenge,fresult);
        cstring(&dummy1[0], &ZEILE_READ[0], 71, 14);
        sscanf(dummy1, "%lf",  &rfile->Fall[FallNr-1].CQ);
        ZR_NULL;
        fgets (ZEILE_READ,MaxZeilenLaenge,fresult);
        cstring(&dummy1[0], &ZEILE_READ[0], 71, 14);
        sscanf(dummy1, "%lf",  &rfile->Fall[FallNr-1].CWI);
        ZR_NULL;
        fgets (ZEILE_READ,MaxZeilenLaenge,fresult);
        cstring(&dummy1[0], &ZEILE_READ[0], 71, 14);
        sscanf(dummy1, "%lf",  &rfile->Fall[FallNr-1].CL);
        ZR_NULL;
        fgets (ZEILE_READ,MaxZeilenLaenge,fresult);
        cstring(&dummy1[0], &ZEILE_READ[0], 71, 14);
        sscanf(dummy1, "%lf",  &rfile->Fall[FallNr-1].CM);
        ZR_NULL;
        fgets (ZEILE_READ,MaxZeilenLaenge,fresult);
        cstring(&dummy1[0], &ZEILE_READ[0], 71, 14);
        sscanf(dummy1, "%lf",  &rfile->Fall[FallNr-1].CN);
        nfgets(11);
        for (lauf =1; lauf <= ParaFile->NbWings; lauf++)
        {
          cstring(&dummy1[0], &ZEILE_READ[0], 71, 14);
          sscanf(dummy1, "%lf",  &rfile->Fall[FallNr-1].FluegelBw[lauf-1].CA);
          ZR_NULL;
          fgets (ZEILE_READ,MaxZeilenLaenge,fresult);
          cstring(&dummy1[0], &ZEILE_READ[0], 71, 14);
          sscanf(dummy1, "%lf",  &rfile->Fall[FallNr-1].FluegelBw[lauf-1].CQ);
          ZR_NULL;
          fgets (ZEILE_READ,MaxZeilenLaenge,fresult);
          cstring(&dummy1[0], &ZEILE_READ[0], 71, 14);
          sscanf(dummy1, "%lf",  &rfile->Fall[FallNr-1].FluegelBw[lauf-1].CWI);  
          ZR_NULL;
          fgets (ZEILE_READ,MaxZeilenLaenge,fresult);
          cstring(&dummy1[0], &ZEILE_READ[0], 71, 14);
          sscanf(dummy1, "%lf",  &rfile->Fall[FallNr-1].FluegelBw[lauf-1].CL);
          ZR_NULL;
          fgets (ZEILE_READ,MaxZeilenLaenge,fresult);
          cstring(&dummy1[0], &ZEILE_READ[0], 71, 14);
          sscanf(dummy1, "%lf",  &rfile->Fall[FallNr-1].FluegelBw[lauf-1].CM);
          ZR_NULL;
          fgets (ZEILE_READ,MaxZeilenLaenge,fresult);
          cstring(&dummy1[0], &ZEILE_READ[0], 71, 14);
          sscanf(dummy1, "%lf",  &rfile->Fall[FallNr-1].FluegelBw[lauf-1].CN);  
          nfgets(7);               
        }
      }
    }
    ZR_NULL;
    fgets (ZEILE_READ,MaxZeilenLaenge,fresult);
  }
  while (!feof(fresult)); 
  fclose(fresult);
  
 /*DEBUG 
  for (lauf=1; lauf<= (ParaFile->numberAlpha+ParaFile->numberCL); lauf++)
  {
    for (zahl=1; zahl <= ParaFile->NbWings; zahl ++)
    {
      printf("Fluegel %d  CA:  %lf          CM: %lf\n",zahl ,rfile->Fall[lauf-1].FluegelBw[zahl-1].CA,rfile->Fall[lauf-1].FluegelBw[zahl-1].CM );
    }
    printf("Fluegel ges  Alfa: %lf     CA:  %lf          CN: %lf\n",rfile->Fall[lauf-1].Alfa,rfile->Fall[lauf-1].CA,rfile->Fall[lauf-1].CN );
    printf("XA2: %lf      CNORM2: %lf\n\n", rfile->Fall[lauf-1].beiwerte[2].XA, rfile->Fall[lauf-1].beiwerte[2].CNORM);
  }*/   
}

void write_load_dist(struct file14 *File14, struct inputfile *input)
{
  int lauf, lauf2, lauf3, lauf4, lauf5;
  double CN, Vwert,Xwert,Lwert,CQ,CL;
  char ExecName[150] = {"\0"}; 
  FILE *fresult;
  fresult = fopen("Load_Distribution.plt","w+");
  fprintf(fresult,"TITLE = \"Coefficient-Distribution of the configuration\" \n");
  fprintf(fresult,"  VARIABLES = X, Y, Z, L, V, CN, CN*l, CL, CL*l, CQ, CQ*l \n");
  for (lauf=0; lauf< (input->numberCL+input->numberAlpha); lauf++)
  {
    for (lauf2 = 0; lauf2 < input->NbWings; lauf2++)
    {
      fprintf(fresult,"ZONE T=\"FALL %d , Fluegel %d\" \n", (lauf+1), (lauf2+1));
      if (lauf2==0)
      {
        lauf4=0;
      }
      for (lauf3=0; lauf3< input->Fluegel[lauf2].AnzahlSpanPanel; lauf3++)
      {
        CN=0;
        for (lauf5=0; lauf5< input->Fluegel[lauf2].AnzahlTiefe;lauf5++)
        {
          CN=CN+File14->Fall[lauf].beiwerte[lauf4+lauf5*input->Fluegel[lauf2].AnzahlSpanPanel].CNORM;
          //printf("LAUFGesamt:%d CN:%lf\n",lauf4+lauf5*input->Fluegel[lauf2].AnzahlSpanPanel,CN);
            }
        CN=CN/input->Fluegel[lauf2].AnzahlTiefe;
        Lwert=input->Fluegel[lauf2].AnzahlTiefe*File14->Fall[lauf].beiwerte[lauf4].LA;
        Xwert=File14->Fall[lauf].beiwerte[lauf4].XA-File14->Fall[lauf].beiwerte[lauf4].LA/2+Lwert/4;
        //printf("CN: %lf Lwert: %lf Xwert:%lf\n",CN, Lwert, Xwert );
        CL=0;
        CQ=0;
        if ((File14->GEODAT[lauf4].Y1-File14->GEODAT[lauf4].Y2)==0)
        {
          CL=0;          
          if (File14->GEODAT[lauf4].Z1>File14->GEODAT[lauf4].Z2)
          {
            CQ=CN*(-1);
            Vwert=90;
          }else{
            Vwert=-90;
            CQ=CN;
          }
        }else if ((File14->GEODAT[lauf4].Z1-File14->GEODAT[lauf4].Z2)==0) {
          CQ=0;          
          if (File14->GEODAT[lauf4].Y1>File14->GEODAT[lauf4].Y2)
          {
            Vwert=0;
            CL=CN; 
          }else{
            CL=-CN;
            Vwert=180;
          }        
        }else{
          if ((File14->GEODAT[lauf4].Z1>File14->GEODAT[lauf4].Z2)&&(File14->GEODAT[lauf4].Y1>File14->GEODAT[lauf4].Y2))
          {
            Vwert=atan((File14->GEODAT[lauf4].Z1-File14->GEODAT[lauf4].Z2)/(File14->GEODAT[lauf4].Y1-File14->GEODAT[lauf4].Y2));
            CL=cos(Vwert)*CN;
            CQ=sin(Vwert)*CN*(-1);
            Vwert=Vwert*180/3.14159265359;
          } else if ((File14->GEODAT[lauf4].Z1<File14->GEODAT[lauf4].Z2)&&(File14->GEODAT[lauf4].Y1>File14->GEODAT[lauf4].Y2)){
            Vwert=atan((File14->GEODAT[lauf4].Z2-File14->GEODAT[lauf4].Z1)/(File14->GEODAT[lauf4].Y1-File14->GEODAT[lauf4].Y2));
            CL=cos(Vwert)*CN;
            CQ=sin(Vwert)*CN;
            Vwert=Vwert*180/3.14159265359*(-1);
          } else if ((File14->GEODAT[lauf4].Z1>File14->GEODAT[lauf4].Z2)&&(File14->GEODAT[lauf4].Y1<File14->GEODAT[lauf4].Y2)){
            Vwert=atan((File14->GEODAT[lauf4].Y2-File14->GEODAT[lauf4].Y1)/(File14->GEODAT[lauf4].Z1-File14->GEODAT[lauf4].Z2));
            CL=sin(Vwert)*CN*(-1);
            CQ=cos(Vwert)*CN*(-1);
            Vwert=Vwert*180/3.14159265359+90;
          } else{
            Vwert=atan((File14->GEODAT[lauf4].Y2-File14->GEODAT[lauf4].Y1)/(File14->GEODAT[lauf4].Z2-File14->GEODAT[lauf4].Z1));
            CL=sin(Vwert)*CN*(-1);
            CQ=cos(Vwert)*CN;
            Vwert=(Vwert*180/3.14159265359+90)*(-1);
          }
        }
        if ((CL==0)&&(CQ==0) &&(CN!=0))
        {
          errmessage(36);
        }                            
        fprintf(fresult," %lf %lf %lf %lf %.2lf %e %e %e %e %e %e\n",Xwert,File14->Fall[lauf].beiwerte[lauf4].YA, File14->Fall[lauf].beiwerte[lauf4].ZA, Lwert, Vwert, CN, CN*Lwert, CL, CL*Lwert,  CQ, CQ*Lwert);
        lauf4++;
      }
      lauf4=lauf4+input->Fluegel[lauf2].AnzahlSpanPanel*(input->Fluegel[lauf2].AnzahlTiefe-1);
    }
  }
  fclose(fresult);  
  strcpy(ExecName,"mv Load_Distribution.plt LaufV1.lili.");
  strcat(ExecName, input->LiLiVersionName);
  system(ExecName);
}
