/****************run_Parametric.c*****************                  
*                                                *      
*      Vorentwurfs-Generator "VoGen"             *          
*      for Version 2.00                          *
*      Date of last modification 11.05.2020      *
*                                                *          
*      Writen by: Jan Himisch                    *      
*             jan.himisch@dlr.de                 *      
*                                                *      
*                                                *
*************************************************/

#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <math.h>
#include "run_Parametric.h"
#include "basis.h"
#include "subsplin.h"
#include "splintab.h"
#include "integration.h"
#include "methode.h"
#include "run_STRUCT.h"
#include "POLAR_INTER.h"
#include "trimming.h"

void getLaufS(struct inputfile  *input, int lauf, double *LaufS, int laufP)
{
  int lauf2;
  int Fluegel, StartS, EndS;
  LaufS[0]=0;
  Fluegel = input->ParaVar[laufP].SplineVarPunkt[lauf].Fluegel-1;
  StartS=input->ParaVar[laufP].SplineVarPunkt[lauf].StartSchnitt;
  EndS=input->ParaVar[laufP].SplineVarPunkt[lauf].EndSchnitt;
  
  for (lauf2 = 0; lauf2 < (EndS-StartS); lauf2++)
  {
    LaufS[lauf2+1]=LaufS[lauf2] + sqrt(pow((input->Fluegel[Fluegel].Schnitte[StartS+lauf2-1].posx-input->Fluegel[Fluegel].Schnitte[StartS+lauf2].posx),2)+pow((input->Fluegel[Fluegel].Schnitte[StartS+lauf2-1].posy-input->Fluegel[Fluegel].Schnitte[StartS+lauf2].posy),2)+pow((input->Fluegel[Fluegel].Schnitte[StartS+lauf2-1].posz-input->Fluegel[Fluegel].Schnitte[StartS+lauf2].posz),2));
  }
}

int getPos(double KUReta, double *eta ,int NrStP)
{
  int lauf;
  for (lauf=1; lauf < NrStP; lauf++)
  {
    if (eta[lauf]  > KUReta)
    {
      break;
    }
  }
  return (lauf-1);
} 

void run_parametric (struct inputfile  *input, int laufSet)
{  
  struct POLINTOUT PolintOut;
  struct OutPutFileName *Out;
  struct file14 resultfile;
  struct PLTdistribution PLTDIST;
  struct BM_FALL *BMDist;
  struct METHODE *method;
  struct GESAMT_POLARint *GPolare;
  int laufP,lauf,lauf2,ANZCHAR, check, checkOptFour;
  double ZwiParaWert;
  FILE *fParaInp, *fParaOut;
  char LINE[600];
  long DatPointer;
  
  checkOptFour=0;
  GPolare=(struct GESAMT_POLARint *)malloc((input->numberCL+input->numberAlpha)*sizeof(struct GESAMT_POLARint));
  
  for (laufP=0;laufP<input->NumberOfParametSetsinput;laufP++)
  {
    printf("Link Numb %d\n",input->ParaVar[laufP].link);
    if (input->ParaVar[laufP].link>0)
    {
      int linknum;
      if (input->ParaVar[laufP].parametric==1)
      {
        for (lauf=0; lauf<input->ParaVar[laufP].anzVar; lauf++)
        {
          printf("Hier\n");
          linknum=input->ParaVar[laufP].link-1;
          printf("Hier Linknum %d laufP %d Fluegel:%d Schnitt:%d\n", linknum, laufP,input->ParaVar[laufP].VarParameter[lauf].Fluegel-1,input->ParaVar[laufP].VarParameter[lauf].Schnitt-1);
          input->Fluegel[input->ParaVar[laufP].VarParameter[lauf].Fluegel-1].Schnitte[input->ParaVar[laufP].VarParameter[lauf].Schnitt-1].posx=input->Fluegel[input->ParaVar[linknum].VarParameter[lauf].Fluegel-1].Schnitte[input->ParaVar[linknum].VarParameter[lauf].Schnitt-1].posx;
          input->Fluegel[input->ParaVar[laufP].VarParameter[lauf].Fluegel-1].Schnitte[input->ParaVar[laufP].VarParameter[lauf].Schnitt-1].posy=input->Fluegel[input->ParaVar[linknum].VarParameter[lauf].Fluegel-1].Schnitte[input->ParaVar[linknum].VarParameter[lauf].Schnitt-1].posy*(-1);
          input->Fluegel[input->ParaVar[laufP].VarParameter[lauf].Fluegel-1].Schnitte[input->ParaVar[laufP].VarParameter[lauf].Schnitt-1].posz=input->Fluegel[input->ParaVar[linknum].VarParameter[lauf].Fluegel-1].Schnitte[input->ParaVar[linknum].VarParameter[lauf].Schnitt-1].posz;
          input->Fluegel[input->ParaVar[laufP].VarParameter[lauf].Fluegel-1].Schnitte[input->ParaVar[laufP].VarParameter[lauf].Schnitt-1].AlphaPlus=input->Fluegel[input->ParaVar[linknum].VarParameter[lauf].Fluegel-1].Schnitte[input->ParaVar[linknum].VarParameter[lauf].Schnitt-1].AlphaPlus;
          input->Fluegel[input->ParaVar[laufP].VarParameter[lauf].Fluegel-1].Schnitte[input->ParaVar[laufP].VarParameter[lauf].Schnitt-1].tiefe=input->Fluegel[input->ParaVar[linknum].VarParameter[lauf].Fluegel-1].Schnitte[input->ParaVar[linknum].VarParameter[lauf].Schnitt-1].tiefe;
          printf("Hier2\n");
        }
      }
      if (input->ParaVar[laufP].parametric==2)
      {
        int zahli,StartS;
        for (lauf=0; lauf<input->ParaVar[laufP].anzSpVar; lauf++)
        {
          for (lauf2=0; lauf2<=fabs(input->ParaVar[laufP].SplineVarPunkt[lauf].EndSchnitt-input->ParaVar[laufP].SplineVarPunkt[lauf].StartSchnitt); lauf2++)
          {
            if (input->ParaVar[laufP].SplineVarPunkt[lauf].EndSchnitt < input->ParaVar[laufP].SplineVarPunkt[lauf].StartSchnitt)
            {
              zahli=lauf2*(-1);
            }else{
              zahli=lauf2;            
            }
            StartS=input->ParaVar[laufP].SplineVarPunkt[lauf].StartSchnitt;
            linknum=input->ParaVar[laufP].link-1;
            printf("Lauf2: %d   zahli  %d  CurSchnitt %d\n", lauf2, zahli,input->ParaVar[laufP].SplineVarPunkt[lauf].StartSchnitt-1+zahli);
            input->Fluegel[input->ParaVar[laufP].SplineVarPunkt[lauf].Fluegel-1].Schnitte[input->ParaVar[laufP].SplineVarPunkt[lauf].StartSchnitt-1+zahli].posx=input->Fluegel[input->ParaVar[linknum].SplineVarPunkt[lauf].Fluegel-1].Schnitte[input->ParaVar[linknum].SplineVarPunkt[lauf].StartSchnitt-1+lauf2].posx;
            input->Fluegel[input->ParaVar[laufP].SplineVarPunkt[lauf].Fluegel-1].Schnitte[input->ParaVar[laufP].SplineVarPunkt[lauf].StartSchnitt-1+zahli].posy=input->Fluegel[input->ParaVar[linknum].SplineVarPunkt[lauf].Fluegel-1].Schnitte[input->ParaVar[linknum].SplineVarPunkt[lauf].StartSchnitt-1+lauf2].posy*(-1);
            input->Fluegel[input->ParaVar[laufP].SplineVarPunkt[lauf].Fluegel-1].Schnitte[input->ParaVar[laufP].SplineVarPunkt[lauf].StartSchnitt-1+zahli].posz=input->Fluegel[input->ParaVar[linknum].SplineVarPunkt[lauf].Fluegel-1].Schnitte[input->ParaVar[linknum].SplineVarPunkt[lauf].StartSchnitt-1+lauf2].posz;
            input->Fluegel[input->ParaVar[laufP].SplineVarPunkt[lauf].Fluegel-1].Schnitte[input->ParaVar[laufP].SplineVarPunkt[lauf].StartSchnitt-1+zahli].AlphaPlus=input->Fluegel[input->ParaVar[linknum].SplineVarPunkt[lauf].Fluegel-1].Schnitte[input->ParaVar[linknum].SplineVarPunkt[lauf].StartSchnitt-1+lauf2].AlphaPlus;
            input->Fluegel[input->ParaVar[laufP].SplineVarPunkt[lauf].Fluegel-1].Schnitte[input->ParaVar[laufP].SplineVarPunkt[lauf].StartSchnitt-1+zahli].tiefe=input->Fluegel[input->ParaVar[linknum].SplineVarPunkt[lauf].Fluegel-1].Schnitte[input->ParaVar[linknum].SplineVarPunkt[lauf].StartSchnitt-1+lauf2].tiefe;
          }
        }
      }
      continue;
    }
    fParaInp = fopen(input->ParaVar[laufP].PARAINP,"r");  
    if (fParaInp ==NULL)
    {
      errmessage(27);
    } 
    
    //Folgender Abschnitt ist auskommentiert, da er einen Fehler verursacht. (Der erste Wert darf nicht negativ sein) 
    /*DatPointer=ftell(fParaInp);
    fgets (LINE,500,fParaInp); 
    
    // tested ob es dich bei der Zeile um eine Komentar-,Leer- oder Parameter-Zeile handelt 
    for (lauf=0; lauf < 50; lauf++)
    {
      if (isspace(LINE[lauf])==0)
      {
        if (isdigit(LINE[lauf])!=0)
        {
          lauf=50;
        }else{
          DatPointer=ftell(fParaInp);
          fgets (LINE,600,fParaInp); 
          lauf=-1;
        }
      }
      if (lauf == 49)
      {
        DatPointer=ftell(fParaInp);
        fgets (LINE,600,fParaInp); 
        lauf=-1;
      }
    }
    // Liest Parameter ein
    fseek(fParaInp, DatPointer, SEEK_SET);*/
    if ((laufSet>1) && (input->ParaVar[laufP].Paraset4eachCLA==1))
    {      
      for (lauf=1;lauf<laufSet;lauf++)
      {
        fgets (LINE,600,fParaInp);
        if (feof(fParaInp))
        {
          errmessage(27);
        }

      }
    }
    printf("ParaDatei: %s \n",input->ParaVar[laufP].PARAINP);
    if (input->ParaVar[laufP].parametric==1)
    {
       for (lauf=0; lauf < input->ParaVar[laufP].anzVar; lauf ++)
       {  
         if (input->ParaVar[laufP].VarParameter[lauf].X==1)
         {
           printf("lese X\n");
           check=fscanf(fParaInp," %lf ",&ZwiParaWert);
           if (check == EOF)
           {
             errmessage(26);
           }
           if (input->ParaVar[laufP].AbsOrDiffer==0)
           {
             input->Fluegel[input->ParaVar[laufP].VarParameter[lauf].Fluegel-1].Schnitte[input->ParaVar[laufP].VarParameter[lauf].Schnitt-1].posx=ZwiParaWert*input->scalfac;
           }else{
             input->Fluegel[input->ParaVar[laufP].VarParameter[lauf].Fluegel-1].Schnitte[input->ParaVar[laufP].VarParameter[lauf].Schnitt-1].posx=input->Fluegel[input->ParaVar[laufP].VarParameter[lauf].Fluegel-1].Schnitte[input->ParaVar[laufP].VarParameter[lauf].Schnitt-1].posx+ZwiParaWert*input->scalfac;
           }
         }
         if (input->ParaVar[laufP].VarParameter[lauf].Y==1)
         {
           printf("lese Y\n");
           check=fscanf(fParaInp," %lf ",&ZwiParaWert);
           if (check == EOF)
           {
             errmessage(26);
           }
           if (input->ParaVar[laufP].AbsOrDiffer==0)
           {
             input->Fluegel[input->ParaVar[laufP].VarParameter[lauf].Fluegel-1].Schnitte[input->ParaVar[laufP].VarParameter[lauf].Schnitt-1].posy=ZwiParaWert*input->scalfac;
           }else{
             input->Fluegel[input->ParaVar[laufP].VarParameter[lauf].Fluegel-1].Schnitte[input->ParaVar[laufP].VarParameter[lauf].Schnitt-1].posy=input->Fluegel[input->ParaVar[laufP].VarParameter[lauf].Fluegel-1].Schnitte[input->ParaVar[laufP].VarParameter[lauf].Schnitt-1].posy+ZwiParaWert*input->scalfac;
           }
         }
         if (input->ParaVar[laufP].VarParameter[lauf].Z==1)
         {
           printf("lese Z\n");
           check=fscanf(fParaInp," %lf ",&ZwiParaWert);
           if (check == EOF)
           {
             errmessage(26);
           }
           if (input->ParaVar[laufP].AbsOrDiffer==0)
           {
             input->Fluegel[input->ParaVar[laufP].VarParameter[lauf].Fluegel-1].Schnitte[input->ParaVar[laufP].VarParameter[lauf].Schnitt-1].posz=ZwiParaWert*input->scalfac;
           }else{
             input->Fluegel[input->ParaVar[laufP].VarParameter[lauf].Fluegel-1].Schnitte[input->ParaVar[laufP].VarParameter[lauf].Schnitt-1].posz=input->Fluegel[input->ParaVar[laufP].VarParameter[lauf].Fluegel-1].Schnitte[input->ParaVar[laufP].VarParameter[lauf].Schnitt-1].posz+ZwiParaWert*input->scalfac;
           }
         }
         if (input->ParaVar[laufP].VarParameter[lauf].tiefe==1)
         {
           printf("lese Tiefe\n");
           check=fscanf(fParaInp," %lf ",&ZwiParaWert);
           if (check == EOF)
           {
             errmessage(26);
           }
           if (input->ParaVar[laufP].AbsOrDiffer==0)
           {
             input->Fluegel[input->ParaVar[laufP].VarParameter[lauf].Fluegel-1].Schnitte[input->ParaVar[laufP].VarParameter[lauf].Schnitt-1].tiefe=ZwiParaWert*input->scalfac;
           }else{
             input->Fluegel[input->ParaVar[laufP].VarParameter[lauf].Fluegel-1].Schnitte[input->ParaVar[laufP].VarParameter[lauf].Schnitt-1].tiefe=input->Fluegel[input->ParaVar[laufP].VarParameter[lauf].Fluegel-1].Schnitte[input->ParaVar[laufP].VarParameter[lauf].Schnitt-1].tiefe+ZwiParaWert*input->scalfac;
           }
         }
         if (input->ParaVar[laufP].VarParameter[lauf].twist==1)
         {
           printf("lese twist\n");
           check=fscanf(fParaInp," %lf ",&ZwiParaWert);
           if (check == EOF)
           {
             errmessage(26);
           }
           if (input->ParaVar[laufP].AbsOrDiffer==0)
           {
             input->Fluegel[input->ParaVar[laufP].VarParameter[lauf].Fluegel-1].Schnitte[input->ParaVar[laufP].VarParameter[lauf].Schnitt-1].AlphaPlus=ZwiParaWert;
           }else{
             input->Fluegel[input->ParaVar[laufP].VarParameter[lauf].Fluegel-1].Schnitte[input->ParaVar[laufP].VarParameter[lauf].Schnitt-1].AlphaPlus=input->Fluegel[input->ParaVar[laufP].VarParameter[lauf].Fluegel-1].Schnitte[input->ParaVar[laufP].VarParameter[lauf].Schnitt-1].AlphaPlus+ZwiParaWert;
           }
         }
         if (input->ParaVar[laufP].VarParameter[lauf].Vloc==1)
         {
           printf("lese Vloc\n");
           check=fscanf(fParaInp," %lf ",&ZwiParaWert);
           if (check == EOF)
           {
             errmessage(26);
           }
           if (input->ParaVar[laufP].AbsOrDiffer==0)
           {
             input->Fluegel[input->ParaVar[laufP].VarParameter[lauf].Fluegel-1].Schnitte[input->ParaVar[laufP].VarParameter[lauf].Schnitt-1].Vloc=ZwiParaWert;
           }else{
             input->Fluegel[input->ParaVar[laufP].VarParameter[lauf].Fluegel-1].Schnitte[input->ParaVar[laufP].VarParameter[lauf].Schnitt-1].Vloc=input->Fluegel[input->ParaVar[laufP].VarParameter[lauf].Fluegel-1].Schnitte[input->ParaVar[laufP].VarParameter[lauf].Schnitt-1].Vloc+ZwiParaWert;
           }
         }       
       }
    }
  
    if (input->ParaVar[laufP].parametric==2)
    {
      int  lauf3, NrStP,lentab;
      double *eta, *value, etaAkt;
      double *LaufS,*bi,*ci,*di;
      double xstep,*xtab,*ytab;
      char DATNAME[50],EXECNAME[100], dummy[4];
      
      FILE *SplOut;
      
      for (lauf =0; lauf < input->ParaVar[laufP].anzSpVar; lauf++)
      {
        //printf("%d Start\n",(input->ParaVar[laufP].SplineVarPunkt[lauf].EndSchnitt-input->ParaVar[laufP].SplineVarPunkt[lauf].StartSchnitt+1));
        LaufS=(double *)malloc((input->ParaVar[laufP].SplineVarPunkt[lauf].EndSchnitt-input->ParaVar[laufP].SplineVarPunkt[lauf].StartSchnitt+1)*sizeof(double));
        if (!LaufS) printf("%lu byte allocation failed\n",(input->ParaVar[laufP].SplineVarPunkt[lauf].EndSchnitt-input->ParaVar[laufP].SplineVarPunkt[lauf].StartSchnitt+1)*sizeof(double));
        getLaufS(input, lauf, LaufS, laufP);
        eta=(double *)malloc((input->ParaVar[laufP].SplineVarPunkt[lauf].AnzahlStuetzpunkte*2)*sizeof(double));
        if (!eta) printf("%lu byte allocation failed\n",(input->ParaVar[laufP].SplineVarPunkt[lauf].AnzahlStuetzpunkte*2)*sizeof(double));
  
        value=(double *)malloc((input->ParaVar[laufP].SplineVarPunkt[lauf].AnzahlStuetzpunkte*2)*sizeof(double));
        for (lauf2=0; lauf2 < input->ParaVar[laufP].SplineVarPunkt[lauf].AnzahlStuetzpunkte; lauf2++)
        {
          check=fscanf(fParaInp," %lf ",&eta[lauf2*2]);
          if (check == EOF)
          {
            errmessage(28);
          }
          /*if ((eta[lauf2*2]>1) || (eta[lauf2*2]<0))
          {
            errmessage(29);
          }*/
          check=fscanf(fParaInp," %lf ",&value[lauf2*2]);
          if (check == EOF)
          {
            errmessage(28);
          }        
        }
        for (lauf2=0; lauf2 < input->ParaVar[laufP].SplineVarPunkt[lauf].AnzahlStuetzpunkte; lauf2++)
        {
          if (lauf2==0)
          {
            eta[1]=eta[0]+deltaS;
            value[1]=(value[2]-value[0])/((eta[2]-deltaS)-eta[0])*(eta[1]-eta[0])+value[0];
          }else if (lauf2==(input->ParaVar[laufP].SplineVarPunkt[lauf].AnzahlStuetzpunkte-1)){
            eta[lauf2*2+1]=eta[lauf2*2];
            value[lauf2*2+1]=value[lauf2*2];
            eta[lauf2*2]=eta[lauf2*2]-deltaS;
            value[lauf2*2]=(value[lauf2*2+1]-value[lauf2*2-1])/(eta[lauf2*2+1]-eta[lauf2*2-1])*(eta[lauf2*2]-eta[lauf2*2-1])+value[lauf2*2-1];
          }else{
            eta[lauf2*2+1]=eta[lauf2*2]+deltaS;
            eta[lauf2*2]=eta[lauf2*2]-deltaS;
            value[lauf2*2+1]=value[lauf2*2];
          }
        }
        /*for (lauf3=0; lauf3 < input->ParaVar[laufP].SplineVarPunkt[lauf].AnzahlStuetzpunkte*2; lauf3++)
        {
          printf("%d : %lf \n ",lauf3, eta[lauf3]);
        } */     
        
        bi=(double *)malloc((input->ParaVar[laufP].SplineVarPunkt[lauf].AnzahlStuetzpunkte*2)*sizeof(double));
        ci=(double *)malloc((input->ParaVar[laufP].SplineVarPunkt[lauf].AnzahlStuetzpunkte*2)*sizeof(double));
        di=(double *)malloc((input->ParaVar[laufP].SplineVarPunkt[lauf].AnzahlStuetzpunkte*2)*sizeof(double));
        NrStP=input->ParaVar[laufP].SplineVarPunkt[lauf].AnzahlStuetzpunkte*2-1;
        check=akima(&NrStP,NrStP,eta,value,0,-1,bi,ci,di);
        if (check != 0)
        {
          printf("Fehler: %d ",check);
          errmessage(30);
        } 
        xstep  = eta[NrStP] - eta[0];
        xstep  = xstep/ (100 - 4 - NrStP);        
        xtab = (REAL *)malloc(100*sizeof(double));        
        ytab = (REAL *)malloc(100*sizeof(double));
        check=sptab(NrStP, eta[0], eta[NrStP-1], xstep, 100 - 1, eta, value, bi, ci, di, xtab, ytab, &lentab);
        strcpy(DATNAME,"akimaSpline_");
        sprintf(dummy,"%d", lauf+1);
        strcat(DATNAME, dummy); 
        strcat(DATNAME, ".dat"); 
        SplOut=fopen(DATNAME,"w+");
        fprintf(SplOut," VARIABLES = \"xtab\" \"ytab\" \n" );
        for (lauf3=0; lauf3< lentab; lauf3++)
        {
          fprintf(SplOut,"%lf  %lf\n ",xtab[lauf3],ytab[lauf3]);
        }
        fclose(SplOut);
        strcpy(DATNAME,"NeueSchnittWerte_");
        sprintf(dummy,"%d", lauf+1);
        strcat(DATNAME, dummy); 
        strcat(DATNAME, ".dat"); 
        SplOut=fopen(DATNAME,"w+");
        fprintf(SplOut," VARIABLES = \"eta\" \"S\" \"X\" \"Y\" \"Z\" \"value\"\n" );
        for (lauf3=0; lauf3<=input->ParaVar[laufP].SplineVarPunkt[lauf].EndSchnitt-input->ParaVar[laufP].SplineVarPunkt[lauf].StartSchnitt; lauf3++)
        {
          int i;  
          double newvalue;
          etaAkt=LaufS[lauf3]/LaufS[input->ParaVar[laufP].SplineVarPunkt[lauf].EndSchnitt-input->ParaVar[laufP].SplineVarPunkt[lauf].StartSchnitt];
          i=getPos(etaAkt, eta,NrStP); 
          newvalue=value[i]+bi[i]*(etaAkt-eta[i])+ci[i]*pow((etaAkt-eta[i]),2)+di[i]*pow((etaAkt-eta[i]),3);
          //printf("i:%d newvalue:%lf etakt:%lf etai:%lf etai+1:%lf\n", i, newvalue, etaAkt, eta[i],eta[i+1]);
          if (input->ParaVar[laufP].SplineVarPunkt[lauf].Variable.X==ON) 
          {            
            if (input->ParaVar[laufP].AbsOrDiffer==0)
            {
              input->Fluegel[input->ParaVar[laufP].SplineVarPunkt[lauf].Fluegel-1].Schnitte[input->ParaVar[laufP].SplineVarPunkt[lauf].StartSchnitt-1+lauf3].posx=newvalue*input->scalfac;
            }else{
              input->Fluegel[input->ParaVar[laufP].SplineVarPunkt[lauf].Fluegel-1].Schnitte[input->ParaVar[laufP].SplineVarPunkt[lauf].StartSchnitt-1+lauf3].posx=input->Fluegel[input->ParaVar[laufP].SplineVarPunkt[lauf].Fluegel-1].Schnitte[input->ParaVar[laufP].SplineVarPunkt[lauf].StartSchnitt-1+lauf3].posx+newvalue*input->scalfac;
            }
          }
          if (input->ParaVar[laufP].SplineVarPunkt[lauf].Variable.Y==ON) 
          {
            if (input->ParaVar[laufP].AbsOrDiffer==0)
            {
              input->Fluegel[input->ParaVar[laufP].SplineVarPunkt[lauf].Fluegel-1].Schnitte[input->ParaVar[laufP].SplineVarPunkt[lauf].StartSchnitt-1+lauf3].posy=newvalue*input->scalfac;
            }else{
              input->Fluegel[input->ParaVar[laufP].SplineVarPunkt[lauf].Fluegel-1].Schnitte[input->ParaVar[laufP].SplineVarPunkt[lauf].StartSchnitt-1+lauf3].posy=input->Fluegel[input->ParaVar[laufP].SplineVarPunkt[lauf].Fluegel-1].Schnitte[input->ParaVar[laufP].SplineVarPunkt[lauf].StartSchnitt-1+lauf3].posy+newvalue*input->scalfac;
            }
          }
          if (input->ParaVar[laufP].SplineVarPunkt[lauf].Variable.Z==ON) 
          {
            if (input->ParaVar[laufP].AbsOrDiffer==0)
            {
              input->Fluegel[input->ParaVar[laufP].SplineVarPunkt[lauf].Fluegel-1].Schnitte[input->ParaVar[laufP].SplineVarPunkt[lauf].StartSchnitt-1+lauf3].posz=newvalue*input->scalfac;
            }else{
              input->Fluegel[input->ParaVar[laufP].SplineVarPunkt[lauf].Fluegel-1].Schnitte[input->ParaVar[laufP].SplineVarPunkt[lauf].StartSchnitt-1+lauf3].posz=input->Fluegel[input->ParaVar[laufP].SplineVarPunkt[lauf].Fluegel-1].Schnitte[input->ParaVar[laufP].SplineVarPunkt[lauf].StartSchnitt-1+lauf3].posz+newvalue*input->scalfac;
            }            
          }
          if (input->ParaVar[laufP].SplineVarPunkt[lauf].Variable.tiefe==ON) 
          {
            if (input->ParaVar[laufP].AbsOrDiffer==0)
            {
              input->Fluegel[input->ParaVar[laufP].SplineVarPunkt[lauf].Fluegel-1].Schnitte[input->ParaVar[laufP].SplineVarPunkt[lauf].StartSchnitt-1+lauf3].tiefe=newvalue*input->scalfac;
            }else{
              input->Fluegel[input->ParaVar[laufP].SplineVarPunkt[lauf].Fluegel-1].Schnitte[input->ParaVar[laufP].SplineVarPunkt[lauf].StartSchnitt-1+lauf3].tiefe=input->Fluegel[input->ParaVar[laufP].SplineVarPunkt[lauf].Fluegel-1].Schnitte[input->ParaVar[laufP].SplineVarPunkt[lauf].StartSchnitt-1+lauf3].tiefe+newvalue*input->scalfac;
            }            
          }
          if (input->ParaVar[laufP].SplineVarPunkt[lauf].Variable.twist==ON) 
          {
            //input->Fluegel[input->ParaVar[laufP].SplineVarPunkt[lauf].Fluegel-1].Schnitte[input->ParaVar[laufP].SplineVarPunkt[lauf].StartSchnitt-1+lauf3].AlphaPlus=newvalue+input->Fluegel[input->ParaVar[laufP].SplineVarPunkt[lauf].Fluegel-1].Schnitte[input->ParaVar[laufP].SplineVarPunkt[lauf].StartSchnitt-1+lauf3].OrgTwist;
            //Mir ist hier nicht mehr ganz klar warum ich eigentlich den Org-Twist speziel zwischenspeicher mustte!!!
            if (input->ParaVar[laufP].AbsOrDiffer==0)
            {
              input->Fluegel[input->ParaVar[laufP].SplineVarPunkt[lauf].Fluegel-1].Schnitte[input->ParaVar[laufP].SplineVarPunkt[lauf].StartSchnitt-1+lauf3].AlphaPlus=newvalue;
            }else{
              input->Fluegel[input->ParaVar[laufP].SplineVarPunkt[lauf].Fluegel-1].Schnitte[input->ParaVar[laufP].SplineVarPunkt[lauf].StartSchnitt-1+lauf3].AlphaPlus=input->Fluegel[input->ParaVar[laufP].SplineVarPunkt[lauf].Fluegel-1].Schnitte[input->ParaVar[laufP].SplineVarPunkt[lauf].StartSchnitt-1+lauf3].AlphaPlus+newvalue;
            }            
          }
          fprintf(SplOut,"%lf  %lf %lf %lf  %lf %lf\n ",etaAkt,LaufS[lauf3], input->Fluegel[input->ParaVar[laufP].SplineVarPunkt[lauf].Fluegel-1].Schnitte[input->ParaVar[laufP].SplineVarPunkt[lauf].StartSchnitt-1+lauf3].posx ,input->Fluegel[input->ParaVar[laufP].SplineVarPunkt[lauf].Fluegel-1].Schnitte[input->ParaVar[laufP].SplineVarPunkt[lauf].StartSchnitt-1+lauf3].posy 
          ,input->Fluegel[input->ParaVar[laufP].SplineVarPunkt[lauf].Fluegel-1].Schnitte[input->ParaVar[laufP].SplineVarPunkt[lauf].StartSchnitt-1+lauf3].posz ,newvalue);
        }
        fclose(SplOut);
        
        free(LaufS);free(eta);free(value); 
        free(bi);free(ci);free(di);
        free(xtab);free(ytab);
      }
    }  
    if (input->ParaVar[laufP].parametric==3)
    {
      int lauf4;
      for (lauf4=0; lauf4<input->AnzahlKlappen;lauf4++)
      {
        check=fscanf(fParaInp," %lf ",&ZwiParaWert);
        if (check == EOF)
        {
          errmessage(26);
        }
        if (input->ParaVar[laufP].AbsOrDiffer==0)
        {
          input->Klappen[lauf4].winkel=ZwiParaWert;
        }else{
          input->Klappen[lauf4].winkel=input->Klappen[lauf4].winkel+ZwiParaWert;
        }            
      }
    }
    if (input->ParaVar[laufP].parametric==4)
    {
      int lauf4;
      checkOptFour=1;
      for (lauf4=0; lauf4<input->AnzahlMixVar;lauf4++)
      {
        check=fscanf(fParaInp," %lf ",&ZwiParaWert);
        if (check == EOF)
        {
          errmessage(26);
        }
        input->MixVar[lauf4]=ZwiParaWert;
      }
    }
    fclose (fParaInp);
  }  
  
  if (checkOptFour==1)
  {
    run_mischer(input);
  }  
  
  // Lifting-Line-Prozesskette Starten  
  if (input->LiliStart==1)
  {
    if (input->Austrimmen > 1)
    {
      Trimm_Schleife(input);
      read14 (&resultfile , input);
      readDistributionPLT (&PLTDIST, input); 
      BMDist=(struct BM_FALL *)malloc((input->numberCL+input->numberAlpha)*sizeof(struct BM_FALL));
      if (!BMDist) printf("%lu byte allocation failed\n",(input->numberCL+input->numberAlpha)*sizeof(struct BM_FALL));
      BM_integration (input, &resultfile, &PLTDIST, BMDist,0);
      method=(struct METHODE *)malloc((input->numberCL+input->numberAlpha)*sizeof(struct METHODE));
      if (!method) printf("%lu byte allocation failed\n",(input->numberCL+input->numberAlpha)*sizeof(struct METHODE));
      methodeAlt(input, method, &resultfile);
      if (input->InterneInterpolation == 1)                      
      {
        PolarInterpolation(input,&PLTDIST,&resultfile,GPolare); 
      }
    }else if (input->Struckt_Coup_ON==1){
      /* if (input->Log_Iter_ON==1)
      {
        printf("!!!Warnung!!!: Log Itteration abgeschaltet!\n");
        input->Log_Iter_ON=0;
      }*/
      run_stru(input);
      read14 (&resultfile , input);
      readDistributionPLT (&PLTDIST, input); 
      BMDist=(struct BM_FALL *)malloc((input->numberCL+input->numberAlpha)*sizeof(struct BM_FALL));
      if (!BMDist) printf("%lu byte allocation failed\n",(input->numberCL+input->numberAlpha)*sizeof(struct BM_FALL));
      BM_integration (input, &resultfile, &PLTDIST, BMDist,0);
      method=(struct METHODE *)malloc((input->numberCL+input->numberAlpha)*sizeof(struct METHODE));
      if (!method) printf("%lu byte allocation failed\n",(input->numberCL+input->numberAlpha)*sizeof(struct METHODE));
      methodeAlt(input, method, &resultfile);     
    }else{  
        run_PROZESS(input);
    }
    if (input->ParaVar[laufP].parametric==2)
    {
      for (lauf =0; lauf < input->ParaVar[laufP].anzSpVar; lauf++)
      {
        char DATNAME[50],EXECNAME[100], dummy[4];
        strcpy(DATNAME,"akimaSpline_");
        sprintf(dummy,"%d", lauf+1);
        strcat(DATNAME, dummy); 
        strcat(DATNAME, ".dat"); 
        strcpy(EXECNAME,"mv ");
        strcat(EXECNAME,DATNAME);
        strcat(EXECNAME," LaufV1.lili.");
        strcat(EXECNAME,input->LiLiVersionName); 
        
        printf ("%s\n", EXECNAME);
        system(EXECNAME);    
    
        strcpy(DATNAME,"NeueSchnittWerte_");
        sprintf(dummy,"%d", lauf+1);
        strcat(DATNAME, dummy); 
        strcat(DATNAME, ".dat");     strcpy(EXECNAME,"mv ");
        strcat(EXECNAME,DATNAME);
        strcat(EXECNAME," LaufV1.lili.");
        strcat(EXECNAME,input->LiLiVersionName); 
        printf ("%s\n", EXECNAME);
        system(EXECNAME);      
      }
    }


    /*write parmeter-out*/
    lauf=0;
    fParaOut=fopen("paramet_out.dat","r");
    if (fParaOut != NULL)
    {
      fseek(fParaOut, -1L, SEEK_END);
      do 
      {
        fseek(fParaOut, -2L, SEEK_CUR);
      }while (fgetc(fParaOut) != 10 ); 
      fscanf(fParaOut,"%d", &lauf);
              
    }else{
      fParaOut = fopen("paramet_out.dat","w");
      fprintf(fParaOut,"TITLE = \" Ergebniss Parameter-Variation\" \n");  
      fprintf(fParaOut,"VARIABLES = \"Lauf\" \"ALPHA\" \"CL\" \"CDI\" ");
      if (input->InterneInterpolation == 1)
      {
        fprintf(fParaOut," \"CDP\" \"CD\"");
      }    
      if ((input->numberCL+input->numberAlpha) > 1)
      {
        fprintf(fParaOut," \"dCL/dalfa\"");
      }
      if (input->BCwi != 0)
      {
        fprintf(fParaOut," \"dCDi [DC]\" ");
      }
      if ((input->BCwi != 0) && (input->BCw0 != 0) && (input->BAlafa != 0) && (input->UrArea != 0) && (input->kHOrd != 0))
      {
        fprintf(fParaOut," \"dCD [DC]\" \"dCDs [DC]\"  \"dCDalfa [DC]\" \"dCDi+dCDs [DC]\" \"unRolledArea [m2]\"");
      }
      if (input->PolintStart==1)
      {
        fprintf(fParaOut," \"CD_POL\" ");
        if (input->BASIS_CD_POL!=0)
        {
          fprintf(fParaOut," \"dCD_POL [DC]\"");	  
        }
      }    
      fprintf(fParaOut," \"WRBM\" ");
      if (input->BWRBM != 0)
      {
        fprintf(fParaOut," \"dWRBM [%]\"");
      }
      fprintf(fParaOut," \"IBM\" ");
      if (input->BasisIntBiegeMoment!=0)
      {
        fprintf(fParaOut," \"dIBM [%]\" ");
      }
      fprintf(fParaOut," \"IDBM\" ");
      if (input->BasisIntDickenBiegeMoment!=0)
      {
        fprintf(fParaOut," \"dIDBM [%]\" ");
      }
      fprintf(fParaOut,"\n");                 
    }
    fclose(fParaOut);
    fParaOut = fopen("paramet_out.dat","a+");
    lauf++;
    for (lauf2=0;lauf2 <(input->numberCL+input->numberAlpha); lauf2++ )
    {
      fprintf(fParaOut," %d %.3lf %.4lf %.8lf",lauf, resultfile.Fall[lauf2].Alfa,  resultfile.Fall[lauf2].CA, resultfile.Fall[lauf2].CWI);
      if (input->InterneInterpolation == 1)
      {
        fprintf(fParaOut," %.8lf %.8lf",GPolare[lauf2].CDP,GPolare[lauf2].CDP+resultfile.Fall[lauf2].CWI+input->BASIS_CD_POL);
      }    
      if ((input->numberCL+input->numberAlpha) > 1)
      {      
        fprintf(fParaOut," %.8lf", (resultfile.Fall[0].CA-resultfile.Fall[1].CA)/(resultfile.Fall[0].Alfa-resultfile.Fall[1].Alfa));
      }        
      if (input->BCwi != 0)
      {
         fprintf(fParaOut," %.8lf", (resultfile.Fall[lauf2].CWI-input->BCwi)*10000);   
      }
      if ((input->BCwi != 0) && (input->BCw0 != 0) && (input->BAlafa != 0) && (input->UrArea != 0) && (input->kHOrd != 0))
      {
        fprintf(fParaOut," %.8lf %.8lf %.8lf %.8lf %.8lf", (method[lauf2].dcw + resultfile.Fall[lauf2].CWI -input->BCwi)*10000, method[lauf2].dcwAerea*10000, method[lauf2].dcwAlfa*10000, (method[lauf2].dcwAerea+resultfile.Fall[lauf2].CWI-input->BCwi)*10000 , resultfile.ABGEWICKELTE_FLAECHE);   
      }    
      if (input->PolintStart==1)
      {
        fprintf(fParaOut," %.8lf",PolintOut.total[lauf2].CFX );
        if (input->BASIS_CD_POL!=0)
        {
          fprintf(fParaOut," %.8lf", (PolintOut.total[lauf2].CFX - input->BASIS_CD_POL)*10000);      
        }
      }    
      fprintf(fParaOut," %e", BMDist[lauf2].WRBM );
      if (input->BWRBM != 0)
      {
        fprintf(fParaOut," %e",(BMDist[lauf2].WRBM-input->BWRBM)/input->BWRBM*100 );
      }
      fprintf(fParaOut," %e", BMDist[lauf2].IBM );
      if (input->BasisIntBiegeMoment!=0)
      {
        fprintf(fParaOut," %e", (BMDist[lauf2].IBM-input->BasisIntBiegeMoment)/input->BasisIntBiegeMoment*100);
      }
      fprintf(fParaOut," %e", BMDist[lauf2].IBMD );
      if (input->BasisIntDickenBiegeMoment!=0)
      {
        fprintf(fParaOut," %e", (BMDist[lauf2].IBMD-input->BasisIntDickenBiegeMoment)/input->BasisIntDickenBiegeMoment*100);
      }
      fprintf(fParaOut,"\n");
    }
    fclose(fParaOut);
  }
}
