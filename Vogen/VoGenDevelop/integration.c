/****************integration.c********************                  
*                                                *      
*      Vorentwurfs-Generator "VoGen"             *          
*      Version 1.15                              *          
*      Date of last modification 16.11.2017      *          
*                                                *          
*      Writen by: Jan Himisch                    *      
*             jan.himisch@dlr.de                 *      
*                                                *      
*                                                *
*************************************************/
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include "input.h"
#include "run_LILI.h"
#include "integration.h"
#include "errmes.h"

struct KreuzList *Klist={NULL};
int *Klist2={NULL};
int AKlist2=0;

int isSchnitt(int SpanPos, int Wing,struct inputfile *input)
{
  int lauf, NRSpanPan,returnval;
  SpanPos++;
  NRSpanPan=0;
  for (lauf = (input->Fluegel[Wing].anzahlSchnitte-1); lauf >=0; lauf--)
  {
    if (SpanPos==NRSpanPan)
    {
      returnval= lauf;
      lauf=-1;
    } else if (NRSpanPan>SpanPos) {
      returnval= -1;
      lauf=-1;
    } else {
      NRSpanPan=NRSpanPan+input->Fluegel[Wing].Schnitte[lauf].AnzahlPan;
    }
  }
  return returnval;
}

int checkKopplung (struct inputfile *input, int Wing, int SpanPos)
{
  int lauf,iist ,lauf2;
  int check=0, SchnittNR;
  iist=0;
  SchnittNR=isSchnitt(SpanPos, Wing, input);


  if (input->AnzKreuz==0)
  {
    return 0;
  }else{
    
    for (lauf=0; lauf < input->AnzKreuz; lauf++)
    { 
      if (Klist==NULL)
      {
        int i;
        Klist=(struct KreuzList *)malloc(input->AnzKreuz*sizeof(struct KreuzList));
        for (i=0; i < input->AnzKreuz; i++)
        {
          Klist[i].KrNr=0;
          Klist[i].AusgangsFluegel=0;
        }
      }
      if (((input->Kreuz[lauf].KreuzFl1==(Wing+1)) && (input->Kreuz[lauf].KreuzPosFl1==(SchnittNR+1)) && (input->Kreuz[lauf].KopelArt!=3)) 
      || ((input->Kreuz[lauf].KreuzFl2==(Wing+1)) &&  (input->Kreuz[lauf].KreuzPosFl2==(SchnittNR+1)) && (input->Kreuz[lauf].KopelArt!=3)))    
      {
        if ((Wing+1)==input->Kreuz[lauf].KreuzFl1)
        {
          if (iist>0)
          {
            for (lauf2=0; lauf2<iist; lauf2++)
            {
              if (Klist[lauf2].AusgangsFluegel==input->Kreuz[lauf].KreuzFl2)
              {
                errmessage(34);
              }
              if (Klist[lauf2].KrNr==lauf)
              {
                check=1;
              }
            }
          }
          if (check==0)
          {
            Klist[iist].KrNr=lauf;
            Klist[iist].AusgangsFluegel=(Wing+1);
            return input->Kreuz[lauf].KreuzFl2;
            check=1;
          }
        }
        if ((Wing+1)==input->Kreuz[lauf].KreuzFl2)
        {
          if (iist>0)
          {
            for (lauf2=0; lauf2<iist; lauf2++)
            {
              if (Klist[lauf2].AusgangsFluegel==input->Kreuz[lauf].KreuzFl1)
              {
                errmessage(34);
              }
              if (Klist[lauf2].KrNr==lauf)
              {
                check=1;
              }
            }
          }
          if (check==0)
          {
            Klist[iist].KrNr=lauf;
            Klist[iist].AusgangsFluegel=(Wing+1);
            return input->Kreuz[lauf].KreuzFl1;
            check=1;
          }
        }        
      }
    }
    if (check==0)
    {
      return 0;
    }
  } 
}

double get_CBM(struct PLTdistribution *PLTdist, struct inputfile *input, int CASE, int FLUEGEL, int SpanPanel,double yakt,double zakt)
{
  int lauf,IsKop;
  double returnValue;
  returnValue=0;
  for (lauf = SpanPanel; lauf >= 0; lauf--)
  {    
    returnValue=returnValue+((PLTdist->Fall[CASE].Fluegel[FLUEGEL].Distri[lauf].Y - yakt)*PLTdist->Fall[CASE].Fluegel[FLUEGEL].Distri[lauf].CFZ*PLTdist->Fall[CASE].Fluegel[FLUEGEL].Distri[lauf].REF_AREA 
    - (PLTdist->Fall[CASE].Fluegel[FLUEGEL].Distri[lauf].Z - zakt)*PLTdist->Fall[CASE].Fluegel[FLUEGEL].Distri[lauf].CFY*PLTdist->Fall[CASE].Fluegel[FLUEGEL].Distri[lauf].REF_AREA
    + PLTdist->Fall[CASE].Fluegel[FLUEGEL].Distri[lauf].CMX*PLTdist->Fall[CASE].Fluegel[FLUEGEL].Distri[lauf].REF_AREA*PLTdist->Fall[CASE].Fluegel[FLUEGEL].Distri[lauf].REF_LEN_CMX)
    / (input->refAerea*input->refSpan); 

    IsKop=checkKopplung(input, FLUEGEL, lauf);
        
    /*Achtung Fluegel 1 kan z.B. kein BM auf Fluegel 2 verursachen. Nur Umgekehrt!!*/
    if ((IsKop!=0) && (IsKop> (FLUEGEL+1))) 
    {
      //printf ("IsKop:%d  SpanPan:%d\n", IsKop-1,input->Fluegel[IsKop-1].AnzahlSpanPanel-1);
      returnValue=returnValue+get_CBM(PLTdist, input, CASE, IsKop-1, input->Fluegel[IsKop-1].AnzahlSpanPanel-1, yakt , zakt);
    }
  }
  if (Klist != NULL)
  {  
    free(Klist);
    Klist=NULL;
  }
  return returnValue;
}

double get_CBT(struct PLTdistribution *PLTdist, struct inputfile *input, int CASE, int FLUEGEL, int SpanPanel,double xakt,double zakt)
{
  int lauf,IsKop;
  double returnValue;
  returnValue=0;
  for (lauf = SpanPanel; lauf >= 0; lauf--)
  {
    returnValue=returnValue+((PLTdist->Fall[CASE].Fluegel[FLUEGEL].Distri[lauf].X - xakt)*PLTdist->Fall[CASE].Fluegel[FLUEGEL].Distri[lauf].CFZ*PLTdist->Fall[CASE].Fluegel[FLUEGEL].Distri[lauf].REF_AREA 
    - (PLTdist->Fall[CASE].Fluegel[FLUEGEL].Distri[lauf].Z - zakt)*PLTdist->Fall[CASE].Fluegel[FLUEGEL].Distri[lauf].CFX_FROM_CDI*PLTdist->Fall[CASE].Fluegel[FLUEGEL].Distri[lauf].REF_AREA
    + PLTdist->Fall[CASE].Fluegel[FLUEGEL].Distri[lauf].CMY*PLTdist->Fall[CASE].Fluegel[FLUEGEL].Distri[lauf].REF_AREA*PLTdist->Fall[CASE].Fluegel[FLUEGEL].Distri[lauf].REF_LEN_CMY)
    / (input->refAerea*input->refSpan); 
    
    IsKop=checkKopplung(input, FLUEGEL, lauf);
        
    /*Achtung Fluegel 1 kan z.B. kein BM auf Fluegel 2 verursachen. Nur Umgekehrt!!*/
    if ((IsKop!=0) && (IsKop> (FLUEGEL+1))) 
    {
      //printf ("IsKop:%d  SpanPan:%d\n", IsKop-1,input->Fluegel[IsKop-1].AnzahlSpanPanel-1);
      returnValue=returnValue+get_CBT(PLTdist, input, CASE, IsKop-1, input->Fluegel[IsKop-1].AnzahlSpanPanel-1, xakt , zakt);
    }
  }
  if (Klist != NULL)
  {  
    free(Klist);
    Klist=NULL;
  }
  return returnValue;
}

double get_CBG(struct PLTdistribution *PLTdist, struct inputfile *input, int CASE, int FLUEGEL, int SpanPanel,double xakt,double yakt)
{
  int lauf,IsKop;
  double returnValue;
  returnValue=0;
  for (lauf = SpanPanel; lauf >= 0; lauf--)
  {

    returnValue=returnValue+((PLTdist->Fall[CASE].Fluegel[FLUEGEL].Distri[lauf].X - xakt)*PLTdist->Fall[CASE].Fluegel[FLUEGEL].Distri[lauf].CFX_FROM_CDI*PLTdist->Fall[CASE].Fluegel[FLUEGEL].Distri[lauf].REF_AREA 
    - (PLTdist->Fall[CASE].Fluegel[FLUEGEL].Distri[lauf].Y - yakt)*PLTdist->Fall[CASE].Fluegel[FLUEGEL].Distri[lauf].CFY*PLTdist->Fall[CASE].Fluegel[FLUEGEL].Distri[lauf].REF_AREA
    + PLTdist->Fall[CASE].Fluegel[FLUEGEL].Distri[lauf].CMZ*PLTdist->Fall[CASE].Fluegel[FLUEGEL].Distri[lauf].REF_AREA*PLTdist->Fall[CASE].Fluegel[FLUEGEL].Distri[lauf].REF_LEN_CMZ)
    / (input->refAerea*input->refSpan); 
    
    IsKop=checkKopplung(input, FLUEGEL, lauf);
        
    /*Achtung Fluegel 1 kan z.B. kein BM auf Fluegel 2 verursachen. Nur Umgekehrt!!*/
    if ((IsKop!=0) && (IsKop> (FLUEGEL+1))) 
    {
      //printf ("IsKop:%d  SpanPan:%d\n", IsKop-1,input->Fluegel[IsKop-1].AnzahlSpanPanel-1);
      returnValue=returnValue+get_CBG(PLTdist, input, CASE, IsKop-1, input->Fluegel[IsKop-1].AnzahlSpanPanel-1, xakt , yakt);
    }
  }
  if (Klist != NULL)
  {  
    free(Klist);
    Klist=NULL;
  }
  return returnValue;
}

double get_FXYZ(struct PLTdistribution *PLTdist, struct inputfile *input, int CASE, int FLUEGEL, int SpanPanel, int XYZ)
{
  int lauf,IsKop;
  double returnValue;
  returnValue=0;
  for (lauf = SpanPanel; lauf >= 0; lauf--)
  {
    if (XYZ==1)
    {
      returnValue=returnValue+PLTdist->Fall[CASE].Fluegel[FLUEGEL].Distri[lauf].CFX_FROM_CDI*PLTdist->Fall[CASE].Fluegel[FLUEGEL].Distri[lauf].REF_AREA/input->refAerea;
    }
    if (XYZ==2)
    {
      returnValue=returnValue+PLTdist->Fall[CASE].Fluegel[FLUEGEL].Distri[lauf].CFY*PLTdist->Fall[CASE].Fluegel[FLUEGEL].Distri[lauf].REF_AREA/input->refAerea;
    }
    if (XYZ==3)
    {
      returnValue=returnValue+PLTdist->Fall[CASE].Fluegel[FLUEGEL].Distri[lauf].CFZ*PLTdist->Fall[CASE].Fluegel[FLUEGEL].Distri[lauf].REF_AREA/input->refAerea;
    }

    IsKop=checkKopplung(input, FLUEGEL, lauf);
        
    /*Achtung Fluegel 1 kan z.B. kein BM auf Fluegel 2 verursachen. Nur Umgekehrt!!*/
    if ((IsKop!=0) && (IsKop> (FLUEGEL+1))) 
    {
      //printf ("IsKop:%d  SpanPan:%d\n", IsKop-1,input->Fluegel[IsKop-1].AnzahlSpanPanel-1);
      returnValue=returnValue+get_FXYZ(PLTdist, input, CASE, IsKop-1, input->Fluegel[IsKop-1].AnzahlSpanPanel-1, XYZ);
    }
  }
  if (Klist != NULL)
  {  
    free(Klist);
    Klist=NULL;
  }
  return returnValue;
}


double getIBMCASE(struct inputfile *input,struct BM_FALL FALL,int fluegel,int cont)
{
  int lauf, lauf2, check;
  double returnvalue;
  if (cont==0)
  {
    returnvalue=FALL.Fluegel[fluegel-1].IBM;
  }else{
    returnvalue=FALL.Fluegel[fluegel-1].IBMD;  
  }
  if (input->AnzKreuz!=0)
  {
    check=0;
    if (Klist2==NULL)
    {
      Klist2=(int *)malloc((input->AnzKreuz+1)*sizeof(int));
      AKlist2=1;
      Klist2[0]=fluegel;
    }    
    for (lauf =0; lauf < input->AnzKreuz; lauf++)
    {
      if ((input->Kreuz[lauf].KreuzFl1==fluegel) && (input->Kreuz[lauf].KreuzFl2>fluegel))
      {
        for (lauf2=0;lauf2<AKlist2; lauf2++)
        {
          if (Klist2[lauf2]==input->Kreuz[lauf].KreuzFl2) {check=1;}
        }
        if ((check==0) && (input->Kreuz[lauf].KopelArt!=3))
        {
          Klist2[AKlist2]=input->Kreuz[lauf].KreuzFl2;
          AKlist2++;
          if (cont==0)
          {
            returnvalue=returnvalue+getIBMCASE(input,FALL, input->Kreuz[lauf].KreuzFl2,0);
          }else{
            returnvalue=returnvalue+getIBMCASE(input,FALL, input->Kreuz[lauf].KreuzFl2,1);          
          }
        }
      }
      if ((input->Kreuz[lauf].KreuzFl2==fluegel) && (input->Kreuz[lauf].KreuzFl1>fluegel))
      {
        for (lauf2=0;lauf2<AKlist2; lauf2++)
        {
          if (Klist2[lauf2]==input->Kreuz[lauf].KreuzFl1) {check=1;}
        }
        if ((check==0) && (input->Kreuz[lauf].KopelArt!=3))
        {
          Klist2[AKlist2]=input->Kreuz[lauf].KreuzFl1;
          AKlist2++;
          if (cont==0)
          {
            returnvalue=returnvalue+getIBMCASE(input,FALL, input->Kreuz[lauf].KreuzFl1,0);
          }else{
            returnvalue=returnvalue+getIBMCASE(input,FALL, input->Kreuz[lauf].KreuzFl1,1);          
          }
        }
      }
    }
    if (Klist2!=NULL)
    {
      free(Klist2);
      Klist2=NULL;
      AKlist2=0;
    }
  }
  return returnvalue;
}

void BM_integration (struct inputfile *input, struct file14  *File14, struct PLTdistribution *PLTdist, struct BM_FALL *BMcase, int Bsurf)
{
  int lauf, lauf2, lauf3, lauf4,lauf5, SNR, NrFaelle;
  
  NrFaelle=input->numberCL+input->numberAlpha;
  if (Bsurf==1)
  {
    NrFaelle++;
  }
  /*Geometrische Eckdaten ermitteln*/
  for (lauf = 0; lauf < NrFaelle; lauf++) 
  {  
    if (Bsurf==1)
    {
      lauf=NrFaelle-1;
    }
    lauf5=0;
    BMcase[lauf].Fluegel=(struct FLUEGELint *)malloc(input->NbWings*sizeof(struct FLUEGELint));
    for (lauf2 = 0; lauf2 < input->NbWings; lauf2++)
    {
      BMcase[lauf].Fluegel[lauf2].distri=(struct DISTRIBUTION *)malloc((input->Fluegel[lauf2].AnzahlSpanPanel+1)*sizeof(struct DISTRIBUTION));
      for (lauf3=0; lauf3< input->Fluegel[lauf2].AnzahlSpanPanel; lauf3++)
      {
        BMcase[lauf].Fluegel[lauf2].distri[lauf3].X=(File14->GEODAT[lauf5].L2*input->Fluegel[lauf2].AnzahlTiefe/4)+(File14->GEODAT[lauf5].X2-File14->GEODAT[lauf5].L2/4);
        BMcase[lauf].Fluegel[lauf2].distri[lauf3].Y=File14->GEODAT[lauf5].Y2;
        BMcase[lauf].Fluegel[lauf2].distri[lauf3].Z=File14->GEODAT[lauf5].Z2;
        BMcase[lauf].Fluegel[lauf2].distri[lauf3].L=File14->GEODAT[lauf5].L2*input->Fluegel[lauf2].AnzahlTiefe;
        BMcase[lauf].Fluegel[lauf2].distri[lauf3].breite=sqrt(pow(File14->GEODAT[lauf5].X1-File14->GEODAT[lauf5].X2,2)+pow(File14->GEODAT[lauf5].Y1-File14->GEODAT[lauf5].Y2,2)+pow(File14->GEODAT[lauf5].Z1-File14->GEODAT[lauf5].Z2,2));
        lauf5++;
        if (input->readRelDick==1)
        {
          SNR=-1;
          lauf4=lauf3;
          SNR=isSchnitt(lauf3, lauf2, input);
          //printf("\n\n SNR:%d\n",SNR);
          if (SNR==-1)
          {
            do 
            {              
              SNR=isSchnitt(lauf4, lauf2, input);
              lauf4--;
              //printf("Lauf3:%d  Lauf4:%d  SNR:%d \n",lauf3,lauf4, SNR);
            }
            while ((SNR==-1)&&(lauf4>0));
          }
          //printf("Lauf3:%d  Lauf4:%d  SNR:%d \n",lauf3,lauf4, SNR);
          if (((lauf3==lauf4) && (lauf4>0))||(SNR==0))
          {
            BMcase[lauf].Fluegel[lauf2].distri[lauf3].Thickness=input->Fluegel[lauf2].Schnitte[SNR].tiefe*input->Fluegel[lauf2].Schnitte[SNR].relDicke;
            //printf("Tiefe:%lf  relDicke:%lf\n",input->Fluegel[lauf2].Schnitte[SNR].tiefe,input->Fluegel[lauf2].Schnitte[SNR].relDicke);
          }else{
            double weg1;
            double weg2;
            if (lauf4<=0)
            {
              SNR=input->Fluegel[lauf2].anzahlSchnitte-1;
            }
            weg1=sqrt(pow((input->Fluegel[lauf2].Schnitte[SNR].posy-BMcase[lauf].Fluegel[lauf2].distri[lauf3].Y),2)+pow((input->Fluegel[lauf2].Schnitte[SNR].posz-BMcase[lauf].Fluegel[lauf2].distri[lauf3].Z),2));
            weg2=sqrt(pow((input->Fluegel[lauf2].Schnitte[SNR-1].posy-BMcase[lauf].Fluegel[lauf2].distri[lauf3].Y),2)+pow((input->Fluegel[lauf2].Schnitte[SNR-1].posz-BMcase[lauf].Fluegel[lauf2].distri[lauf3].Z),2));
            BMcase[lauf].Fluegel[lauf2].distri[lauf3].Thickness=(input->Fluegel[lauf2].Schnitte[SNR-1].tiefe*input->Fluegel[lauf2].Schnitte[SNR-1].relDicke - input->Fluegel[lauf2].Schnitte[SNR].tiefe*input->Fluegel[lauf2].Schnitte[SNR].relDicke)
            /(weg1+weg2)*weg1+(input->Fluegel[lauf2].Schnitte[SNR].tiefe*input->Fluegel[lauf2].Schnitte[SNR].relDicke);
            //printf("weg1:%lf weg2:%lf\n", weg1,weg2);
          }
          //printf("X:%lf Y:%lf Z:%lf TH:%lf\n",BMcase[lauf].Fluegel[lauf2].distri[lauf3].X,BMcase[lauf].Fluegel[lauf2].distri[lauf3].Y,BMcase[lauf].Fluegel[lauf2].distri[lauf3].Z,BMcase[lauf].Fluegel[lauf2].distri[lauf3].Thickness);
        }
      }
      lauf5=lauf5+input->Fluegel[lauf2].AnzahlSpanPanel*(input->Fluegel[lauf2].AnzahlTiefe-1);
    }
  }
  //printf ("Geometriepostionen ermittelt!\n");
  //printf ("\n");
  /*Biegemomentverteilung ermitteln*/
  for (lauf = 0; lauf < NrFaelle; lauf++) 
  {  
    if (Bsurf==1)
    {
      lauf=NrFaelle-1;
    }
    for (lauf2 = 0; lauf2 < input->NbWings; lauf2++)
    {
      for (lauf3=0; lauf3< input->Fluegel[lauf2].AnzahlSpanPanel; lauf3++)
      {
        BMcase[lauf].Fluegel[lauf2].distri[lauf3].CBM=get_CBM(PLTdist, input, lauf, lauf2, lauf3,BMcase[lauf].Fluegel[lauf2].distri[lauf3].Y,BMcase[lauf].Fluegel[lauf2].distri[lauf3].Z);
        // Hier Stand BMcase[lauf].Fluegel[lauf2].distri[lauf3+1].CBM=get_CBM(PLTdist, input, lauf, lauf2, lauf3,BMcase[lauf].Fluegel[lauf2].distri[lauf3].Y,BMcase[lauf].Fluegel[lauf2].distri[lauf3].Z);
        // Ich kann nicht nachvollziehen ob hinter lauf3+1 nicht doch ein Sinn stand. Noch mal durch rechenen!
        BMcase[lauf].Fluegel[lauf2].distri[lauf3].CBT=get_CBT(PLTdist, input, lauf, lauf2, lauf3,BMcase[lauf].Fluegel[lauf2].distri[lauf3].X,BMcase[lauf].Fluegel[lauf2].distri[lauf3].Z);
        BMcase[lauf].Fluegel[lauf2].distri[lauf3].CBG=get_CBG(PLTdist, input, lauf, lauf2, lauf3,BMcase[lauf].Fluegel[lauf2].distri[lauf3].X,BMcase[lauf].Fluegel[lauf2].distri[lauf3].Y);

        BMcase[lauf].Fluegel[lauf2].distri[lauf3].FX=get_FXYZ(PLTdist, input, lauf, lauf2, lauf3,1);
        BMcase[lauf].Fluegel[lauf2].distri[lauf3].FY=get_FXYZ(PLTdist, input, lauf, lauf2, lauf3,2);
        BMcase[lauf].Fluegel[lauf2].distri[lauf3].FZ=get_FXYZ(PLTdist, input, lauf, lauf2, lauf3,3);
      }
    }
  }  
  //printf("CBM ermittelt!\n");

  /*IBM ermitteln pro Flügel ermitteln*/
  for (lauf = 0;  lauf < NrFaelle; lauf++) 
  {
    if (Bsurf==1)
    {
      lauf=NrFaelle-1;
    }
    for (lauf2 = 0; lauf2 < input->NbWings; lauf2++)
    {
      int check;
      check=-1;
      for (lauf3=0; lauf3< input->Fluegel[lauf2].AnzahlSpanPanel; lauf3++)
      {
        SNR=-1;
        if (input->IBM_BD_WING==(lauf2+1))
        {
          SNR=isSchnitt(lauf3, lauf2, input);
          if (SNR!=-1)
          {
            check=SNR;
          }
        }
        //printf("LAUF3input->numberCL+input->numberAlpha;:%d SNR:%d CHECK:%d\n",lauf3,SNR, check);
        if (lauf3==0)
        {
          BMcase[lauf].Fluegel[lauf2].IBM=0;
        }
        if (((lauf2+1)>=input->IBM_BD_WING) && ((check==-1) || ((check+1)>= input->IBM_BD_SEC)))
        {
          BMcase[lauf].Fluegel[lauf2].IBM=BMcase[lauf].Fluegel[lauf2].IBM+BMcase[lauf].Fluegel[lauf2].distri[lauf3].CBM*BMcase[lauf].Fluegel[lauf2].distri[lauf3].breite;
          if ((check+1)==input->IBM_BD_SEC)
          {
            break;
          }
        }
      }
    }
  }
  //printf("IBM pro Fluegel ermittelt!\n");
  
  /*IBMD ermitteln pro Flügel ermitteln*/
  if (input->readRelDick==1)
  {
    for (lauf = 0;  lauf < NrFaelle; lauf++) 
    {
      if (Bsurf==1)
      {
        lauf=NrFaelle-1;
      }
      for (lauf2 = 0; lauf2 < input->NbWings; lauf2++)
      {
        int check;
        check=-1;
        for (lauf3=0; lauf3< input->Fluegel[lauf2].AnzahlSpanPanel; lauf3++)
        {
          SNR=-1;
          if (input->IBM_BD_WING==(lauf2+1))
          { 
            SNR=isSchnitt(lauf3, lauf2, input);
            if (SNR!=-1)
            {
              check=SNR;
            }
          }
          if (lauf3==0)
          {
            BMcase[lauf].Fluegel[lauf2].IBMD=0;
          }
          if (((lauf2+1)>=input->IBM_BD_WING) && ((check==-1) || ((check+1)>= input->IBM_BD_SEC)))
          {
            BMcase[lauf].Fluegel[lauf2].IBMD=BMcase[lauf].Fluegel[lauf2].IBMD+BMcase[lauf].Fluegel[lauf2].distri[lauf3].CBM*BMcase[lauf].Fluegel[lauf2].distri[lauf3].breite/BMcase[lauf].Fluegel[lauf2].distri[lauf3].Thickness;
            if ((check+1)==input->IBM_BD_SEC)
            {
              break;
            }
          }            
        }
      }
    }
    //printf("IBMD pro Fluegel ermittelt!\n");
  }
  
  /*Gesamt IBM ermitteln Flügel 1 ist die Referenz*/ 
  for (lauf = 0;  lauf < NrFaelle ; lauf++) 
  {
    if (Bsurf==1)
    {
      lauf=NrFaelle-1;
    }
    BMcase[lauf].IBM=getIBMCASE(input, BMcase[lauf],1,0);
  }
  
  /*Gesamt IBMD ermitteln Flügel 1 ist die Referenz*/ 
  if (input->readRelDick==1)
  {
    for (lauf = 0;  lauf < NrFaelle; lauf++) 
    {
      if (Bsurf==1)
      {
        lauf=NrFaelle-1;
      }
      BMcase[lauf].IBMD=getIBMCASE(input, BMcase[lauf],1,1);
    }
  }
  
  /*WRBM ermitteln*/
  if (input->WRBM_BD_WING>0)
  {
    for (lauf=0;lauf< NrFaelle; lauf++) 
    {
      if (Bsurf==1)
      {
        lauf=NrFaelle-1;
      }
      if(input->WRBM_BD_SEC>0)
      {
        for (lauf3=0; lauf3< input->Fluegel[input->WRBM_BD_WING-1].AnzahlSpanPanel; lauf3++)
	{
	  SNR=isSchnitt(lauf3, input->WRBM_BD_WING-1, input);
	  if (SNR==(input->WRBM_BD_SEC-1))
	  {
	    BMcase[lauf].WRBM=BMcase[lauf].Fluegel[input->WRBM_BD_WING-1].distri[lauf3].CBM;
	    break;
	  }
	}
      }else{
        BMcase[lauf].WRBM=BMcase[lauf].Fluegel[input->WRBM_BD_WING-1].distri[input->Fluegel[input->WRBM_BD_WING-1].AnzahlSpanPanel-1].CBM;
      }
    }
  }else{
    for (lauf=0;lauf< NrFaelle; lauf++) 
    {
      if (Bsurf==1)
      {
        lauf=NrFaelle-1;
      }
      BMcase[lauf].WRBM=BMcase[lauf].Fluegel[0].distri[input->Fluegel[0].AnzahlSpanPanel-1].CBM;
    }  
  }
}

void write_BM_distribution(struct BM_FALL *BMcase, struct file14  *File14, struct inputfile *input, int schalt)
{
  int lauf, lauf2, lauf3, lauf4, lauf5,NrFaelle;
  double  CN;
  char ExecName[150] = {"\0"}; 
  FILE *fout, *f2out;
   
  fout=fopen("BM_Distribution.plt","w+");
  fprintf(fout,"TITLE = \"Coefficient-Distribution of the configuration\" \n");
  if (schalt==1)
  {
    fprintf(fout,"  VARIABLES = X, Y, Z, CMB \n"); 
  }else{
    fprintf(fout,"  VARIABLES = X, Y, Z, CN, CMB \n");  
  }
  f2out=fopen("BM_Coefficents.plt","w+");
  fprintf(f2out,"TITLE = \"Coefficient-BM\" \n");
  if (schalt==1)
  {
    fprintf(f2out,"VARIABLES =");
  }else{
    fprintf(f2out,"VARIABLES = CL, CWI");  
  }
  for (lauf2 = 0; lauf2 < input->NbWings; lauf2++)
  {
    if (lauf2==0 && schalt==1)
    {
      fprintf(f2out," IBM_W%d",lauf2+1);    
    }else{
      fprintf(f2out,", IBM_W%d",lauf2+1);
    }
    if (input->readRelDick==1)
    {
      fprintf(f2out,", IBMD_W%d",lauf2+1);
    }
  }     
  fprintf(f2out,", WRBM, IBM");  
  if (input->readRelDick==1)
  {
    fprintf(f2out,", IBMD");
  }
  fprintf(f2out,"\n");  
  NrFaelle=input->numberCL+input->numberAlpha;
  if (schalt==1)
  {
    NrFaelle++;
  }
  for (lauf=0; lauf< NrFaelle; lauf++)
  {
    if (schalt==1)
    {
      lauf=NrFaelle-1;
    }
    if (schalt!=1)
    {
      fprintf(f2out," %e \t%e", File14->Fall[lauf].CA, File14->Fall[lauf].CWI);
    }
    for (lauf2 = 0; lauf2 < input->NbWings; lauf2++)
    {
      fprintf(fout,"ZONE T=\"FALL %d , Fluegel %d\" \n", (lauf+1), (lauf2+1));
      if (lauf2==0)
      {
        lauf4=0;
      }
      for (lauf3=0; lauf3< input->Fluegel[lauf2].AnzahlSpanPanel; lauf3++)
      {
        CN=0;
        if (schalt==0)
        {
          for (lauf5=0; lauf5< input->Fluegel[lauf2].AnzahlTiefe;lauf5++)
          {
            CN=CN+File14->Fall[lauf].beiwerte[lauf4+lauf5*input->Fluegel[lauf2].AnzahlSpanPanel].CNORM2;
          }
          CN=CN/input->Fluegel[lauf2].AnzahlTiefe;
          fprintf(fout," %lf %lf %lf %e %e\n",BMcase[lauf].Fluegel[lauf2].distri[lauf3].X,BMcase[lauf].Fluegel[lauf2].distri[lauf3].Y,BMcase[lauf].Fluegel[lauf2].distri[lauf3].Z ,CN, BMcase[lauf].Fluegel[lauf2].distri[lauf3].CBM);
        }else{
          fprintf(fout," %lf %lf %lf %e\n",BMcase[lauf].Fluegel[lauf2].distri[lauf3].X,BMcase[lauf].Fluegel[lauf2].distri[lauf3].Y,BMcase[lauf].Fluegel[lauf2].distri[lauf3].Z , BMcase[lauf].Fluegel[lauf2].distri[lauf3].CBM);        
        }
        lauf4++;
      }

      lauf4=lauf4+input->Fluegel[lauf2].AnzahlSpanPanel*(input->Fluegel[lauf2].AnzahlTiefe-1);
      printf("  IBM wing: %d  case %d  = %e\n", lauf2+1, lauf+1, BMcase[lauf].Fluegel[lauf2].IBM);
      fprintf(f2out," \t%e", BMcase[lauf].Fluegel[lauf2].IBM);
      if (input->readRelDick==1)
      {
        printf("  IBMD wing: %d  case %d  = %e\n", lauf2+1, lauf+1, BMcase[lauf].Fluegel[lauf2].IBMD);
	fprintf(f2out," \t%e", BMcase[lauf].Fluegel[lauf2].IBMD);
      }
    }
    printf("  IBM case %d  = %e\n", lauf+1, BMcase[lauf].IBM);
    printf("  WRBM case %d  = %e\n", lauf+1, BMcase[lauf].WRBM);
    fprintf(f2out," \t%e \t%e", BMcase[lauf].WRBM, BMcase[lauf].IBM);
    if (input->readRelDick==1)
    {
      printf("  IBMD case %d  = %e\n", lauf+1, BMcase[lauf].IBMD);
      fprintf(f2out," \t%e", BMcase[lauf].IBMD);
    }
    fprintf(f2out," \n");
  }
  fclose(fout);
  fclose(f2out);
  if (schalt == 0)
  {
    strcpy(ExecName,"mv BM_Distribution.plt LaufV1.lili.");
    strcat(ExecName, input->LiLiVersionName);
    system(ExecName);
    strcpy(ExecName,"mv BM_Coefficents.plt LaufV1.lili.");
    strcat(ExecName, input->LiLiVersionName);
    system(ExecName);

  } else {
    system("mv BM_Distribution.plt BM_TAU_Distribution.plt");
    system("mv BM_Coefficents.plt BM_TAU_Coefficents.plt");  
  }
}
