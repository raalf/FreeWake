/***************create_LiLiInp.c******************                  
*                                                *      
*      Vorentwurfs-Generator "VoGen"             *          
*      Version 1.12                              *          
*      Date of last modification 02.01.2014      *          
*                                                *          
*      Writen by: Jan Himisch                    *      
*             jan.himisch@dlr.de                 *      
*                                                *      
*                                                *
*************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "create_LiLiInp.h"

const int ZusAufPunkte = 0; 
const double GeoScale = 1;

const char HEAD1[]= "|===============================================================================================================================|\n";
const char HEAD2[]= " LIFTING_LINE INPUTFILE\n";
const char HEAD3[]= "|-------------------------------------------------------------------------------------------------------------------------------|\n";
//const char HEAD4[]= " VERSION V2.3beta \n";
const char HEAD5[]= "|---------------------------------------------------------NAME_OF_DATASET-------------------------------------------------------|\n";
//const char HEAD6[]= " LIFTING_LINE V2.3beta , EXAMPLE: SIMPLE\n";
const char HEAD7[]= "|SYMM|COMPRESS_MA-NO|-GEO_SCALING--|N_PART_WINGS|N_ANG_OF_ATTACK|N_ADD_LINES|RESULT_CHANNEL|XMLDAT|-------------------------------|\n";
const char HEAD7_KreisFlug[]= "|SYMM|COMPRESS_MA-NO|-GEO_SCALING--|N_PART_WINGS|N_ANG_OF_ATTACK|CIRCLING_MODE|N_ADD_LINES|RESULT_CHANNEL|XMLDAT|T_DISTRIB|-------|\n";
const char HEAD8[]= "|---REF_AREA---|---REF_SPAN---|-REF_LEN_CMX--|-REF_LEN_CMY--|-REF_LEN_CMZ--|--MOM_REF_X---|--MOM_REF_Y---|--MOM_REF_Z---|---------|\n";
const char HEAD9[]= "|I_PW|N_P-|----X_PW_1----|----Y_PW_1----|----Z_PW_1----|--CHORD_PW_1--|----X_PW_2----|----Y_PW_2----|----Z_PW_2----|--CHORD_PW_2--|\n";
const char HEAD10[]= "|COUPLING_CONDITIONS_1|COUPLING_CONDITIONS_2|---TWIST_1----|---TWIST_2----|FACT_CL-ALPHA-|FUSELAGE_INF|P_DISTRIB|ELLIPT_CHORD|----|\n";
const char HEAD11[]= "|-----BETA-----|-----ALPHA-----|\n";
const char HEAD12[]= "|N_TARGET-CFZ|--TARGET-CFZ--|--TARGET-CFZ--|--TARGET-CFZ--|--TARGET-CFZ--|--TARGET-CFZ--|--TARGET-CFZ--|--TARGET-CFZ--|--TARGET_...\n";
const char HEAD13[]= "|N_WING|--------------------------------------------------------------------------------------------------------------------------|\n";
const char HEAD14[]= "|N_PW_SPAN|N_PW_CHORD|------------------------------------------------------------------------------------------------------------|\n";
const char HEAD15[]= "|-----XPA------|-----YPA------|-----ZPA------|-----XPE------|-----YPE------|-----ZPE------|-NPU-|---------------------------------|\n";
const char HEAD16[]= "|KINEMAT_VISCOS|-VELOCITY_INF-|-MACH-NO_INF--|------------------------------------------------------------------------------------|\n";
const char HEAD17[]= "|N_AIRFOILS|----------------------------------------------------------------------------------------------------------------------|\n";
const char HEAD18[]= "|AIRFOILNAME|---------------------------------------------------------------------------------------------------------------------|\n";
const char HEAD19[]= "|N_ETA_COORD|ETA_COORDINATE|ETA_COORDINATE|ETA_COORDINATE|ETA_COORDINATE|ETA_COORDINATE|ETA_COORDINATE|ETA_COORDINATE|ETA_COORD...\n";
const char HEAD20[]= "|==================================================END OF LIFTING_LINE INPUTFILE==================================================|\n";
const char KreisflugPara[]="|-----Vk[m/s]----|-----PHI[°]-----|-----Ws[m/s]----|--gradient[1/s]-|----------------|------------CIRCLING PARAMETERS-------------|\n";
const char ARAND[]= "T F F F F F ";
const char IRAND[]= "F T F F F F ";
const char ZRAND[]= "F F T T F F ";
const char ZRAND2[]= "F F T F F T ";
const char ZRAND3[]= "F F T F F F ";
const char KRAND[]= "F F T F F T ";
const char K2RAND[]= "F F T F F F ";

int *Kreuznr=NULL;

double getL (struct inputfile *ParaFile,int laufFL,int laufTiefe,int laufSchnitte)
{
  double LPANEL;
  if (ParaFile->Fluegel[laufFL].Schnitte[laufSchnitte].anzPanelKlappe>0)
  {
     if ((laufTiefe) < (ParaFile->Fluegel[laufFL].AnzahlTiefe-ParaFile->Fluegel[laufFL].Schnitte[laufSchnitte].anzPanelKlappe))
     {
       LPANEL=(1.0-ParaFile->Fluegel[laufFL].Schnitte[laufSchnitte].etaSchnitt)*ParaFile->Fluegel[laufFL].Schnitte[laufSchnitte].tiefe/(ParaFile->Fluegel[laufFL].AnzahlTiefe-ParaFile->Fluegel[laufFL].Schnitte[laufSchnitte].anzPanelKlappe);
     }else{
       LPANEL=ParaFile->Fluegel[laufFL].Schnitte[laufSchnitte].etaSchnitt*ParaFile->Fluegel[laufFL].Schnitte[laufSchnitte].tiefe/ParaFile->Fluegel[laufFL].Schnitte[laufSchnitte].anzPanelKlappe;     
     }
  }else{
    LPANEL=ParaFile->Fluegel[laufFL].Schnitte[laufSchnitte].tiefe/ParaFile->Fluegel[laufFL].AnzahlTiefe;
  }
  return LPANEL;
}

void XYZtransformation(struct inputfile *ParaFile,int laufFL,int laufTiefe,int laufSchnitt, double *XWERT,double *YWERT,double *ZWERT)
{
  int laufT;
  double SchnittTiefe;	
  SchnittTiefe=ParaFile->Fluegel[laufFL].Schnitte[laufSchnitt].tiefe/ParaFile->Fluegel[laufFL].AnzahlTiefe;
  *XWERT=0.0;
  for (laufT=0; laufT < laufTiefe; laufT++)
  {
    *XWERT=*XWERT+getL(ParaFile, laufFL, laufT, laufSchnitt);
  }
  *XWERT=*XWERT+getL(ParaFile, laufFL, laufT, laufSchnitt)*0.25;
  /*Wert+NasenPunkt +Drehen um Target Alfa oder Ziel Alfa */
  *XWERT=*XWERT+ ParaFile->Fluegel[laufFL].Schnitte[laufSchnitt].posx * cos(ParaFile->Alfa_Rot * 3.141592653/180*-(1))+ ParaFile->Fluegel[laufFL].Schnitte[laufSchnitt].posz * sin(ParaFile->Alfa_Rot * 3.141592653/180*-(1));
  *YWERT= ParaFile->Fluegel[laufFL].Schnitte[laufSchnitt].posy;
  *ZWERT= ParaFile->Fluegel[laufFL].Schnitte[laufSchnitt].posx * sin(ParaFile->Alfa_Rot * 3.141592653/180*-(1))*(-1)+ ParaFile->Fluegel[laufFL].Schnitte[laufSchnitt].posz * cos(ParaFile->Alfa_Rot * 3.141592653/180*-(1));
}	

int IsKoppel(int lauf1,int lauf2,int lauf3,int AKreuz,struct Kreuzung *Kreuz, struct Einfluegel *Fluegel)
{
  int zahl, zahl2, returnflag, zwiNr1, zwiNr2;
  returnflag=0;
  if (AKreuz==0)
  {
    return 0;
  }else{
    if (Kreuznr == NULL)
    {
      Kreuznr = (int *)malloc(AKreuz*sizeof (int)); 
      for (zahl=0; zahl<AKreuz; zahl++)
      {
        zwiNr1=Fluegel[Kreuz[zahl].KreuzFl1].AnzahlTiefe-(Kreuz[zahl].KreuzTiefFl1-1);
        zwiNr2=Fluegel[Kreuz[zahl].KreuzFl2].AnzahlTiefe-(Kreuz[zahl].KreuzTiefFl2-1);	
        if (zwiNr1>zwiNr2)
        {
          Kreuznr[zahl]=zwiNr2;      
	}else{
          Kreuznr[zahl]=zwiNr1; 
	}
      }
    }
    for (zahl=0; zahl<AKreuz; zahl++)
    {
      //printf("KreuzSchnitt:%d  LAuf1: %d Lauf2: %d Lauf3: %d\n",Kreuz[zahl].KreuzPosFl1, lauf1, lauf2, lauf3);
      if ((lauf1==(Kreuz[zahl].KreuzFl1-1)) && (lauf3==Kreuz[zahl].KreuzPosFl1) && ((Kreuz[zahl].KopelArt==1)||(Kreuz[zahl].KopelArt==3)))
      {
	if ((lauf2 >= Kreuz[zahl].KreuzTiefFl1)  && ((lauf2 - Kreuz[zahl].KreuzTiefFl1)<= (Fluegel[Kreuz[zahl].KreuzFl2-1].AnzahlTiefe - Kreuz[zahl].KreuzTiefFl2)))
        {
	  returnflag=lauf2+1-Kreuz[zahl].KreuzTiefFl1;
	  for (zahl2=0;zahl2<zahl;zahl2++)
          {
            returnflag=returnflag+Kreuznr[zahl2];
	  }
	}
      }   
      if ((lauf1==(Kreuz[zahl].KreuzFl2-1)) && (lauf3==Kreuz[zahl].KreuzPosFl2) && ((Kreuz[zahl].KopelArt==1)||(Kreuz[zahl].KopelArt==3)))
      {
        if ((lauf2 >= Kreuz[zahl].KreuzTiefFl2)  && (lauf2 - Kreuz[zahl].KreuzTiefFl2<= Fluegel[Kreuz[zahl].KreuzFl1-1].AnzahlTiefe - Kreuz[zahl].KreuzTiefFl1))
        {
          returnflag=lauf2+1-Kreuz[zahl].KreuzTiefFl2;
	  for (zahl2=0;zahl2<zahl;zahl2++)
          {
            returnflag=returnflag+Kreuznr[zahl2];
	  }
	}
      }   
    }
  }
  //printf ("Returnflag:%d\n",returnflag);
  return returnflag;  
}

int whichKreuzNR(int lauf1,int lauf2,int lauf3,int AKreuz,struct Kreuzung *Kreuz,  struct Einfluegel *Fluegel)
{
  int zahl, returnflag;
  returnflag=0;
  for (zahl=0; zahl<AKreuz; zahl++)
  {
    if ((lauf1==(Kreuz[zahl].KreuzFl1-1)) && (lauf3==Kreuz[zahl].KreuzPosFl1))
    {
      if ((lauf2 >= Kreuz[zahl].KreuzTiefFl1)  && ((lauf2 - Kreuz[zahl].KreuzTiefFl1)<= (Fluegel[Kreuz[zahl].KreuzFl2-1].AnzahlTiefe - Kreuz[zahl].KreuzTiefFl2)))
      {
        returnflag=zahl;
      }    
    }
    if ((lauf1==(Kreuz[zahl].KreuzFl2-1)) && (lauf3==Kreuz[zahl].KreuzPosFl2))
    {
      if ((lauf2 >= Kreuz[zahl].KreuzTiefFl2)  && (lauf2 - Kreuz[zahl].KreuzTiefFl2<= Fluegel[Kreuz[zahl].KreuzFl1-1].AnzahlTiefe - Kreuz[zahl].KreuzTiefFl1))
      {
        returnflag=zahl;
      }
    }  
  }
  return returnflag;
}

/*Zur Verbesserung der Genauigkeit Camber Bestimmung als Durchschnitswert des betrachteten Abschnittes*/
double Get_Camber2(struct EinSchnitt *Schnitt, int tiefe, int maxtiefe,struct EineKlappen *Klappen, int IArand, double startCamb, double endCamb)
{
  int lauf, XYpos1, XYpos2;
  double searchval, searchval2, Zval1, Zval2, Xval1, Xval2, Winkel;
  
  XYpos1=-1;
  XYpos2=-1;
  if (Schnitt->anzPanelKlappe>0)
  {
    if ((tiefe+1)<=(maxtiefe-Schnitt->anzPanelKlappe))
    { 
      searchval = (double)(1-Schnitt->etaSchnitt)*(tiefe+startCamb) * 1/(maxtiefe-Schnitt->anzPanelKlappe);  
    }else{
      searchval = (double)(1-Schnitt->etaSchnitt)+Schnitt->etaSchnitt/Schnitt->anzPanelKlappe*(tiefe-(maxtiefe-Schnitt->anzPanelKlappe)+startCamb);      
    }
  }else{
    searchval = (double)(tiefe+0.5) * 1/maxtiefe;  
  }
  for (lauf = 0; lauf <= Schnitt->AnzahlXY; lauf++)
  {
    if ((XYpos1==-1) && (Schnitt->xcamb[lauf]> searchval))
    {
      XYpos1=lauf;
      break;
    }
  }
  if (Schnitt->anzPanelKlappe>0)
  {
    if ((tiefe+1)<=(maxtiefe-Schnitt->anzPanelKlappe))
    {
      searchval2 = (double)(1-Schnitt->etaSchnitt)*(tiefe+endCamb) * 1/(maxtiefe-Schnitt->anzPanelKlappe);  
    }else{
      searchval2 = (double)(1-Schnitt->etaSchnitt)+Schnitt->etaSchnitt/Schnitt->anzPanelKlappe*(tiefe-(maxtiefe-Schnitt->anzPanelKlappe)+endCamb);      
    }
  }else{
    searchval2 = (double)(tiefe+1) * 1/maxtiefe;  
  }
  
  if (searchval2 !=1.0) 
  {
    for (lauf = 0; lauf <= Schnitt->AnzahlXY; lauf++)
    {
      if ((XYpos2==-1) && (Schnitt->xcamb[lauf]> searchval2))
      {
        XYpos2=lauf;
	break;
      }
    }
  } else {
    XYpos2=Schnitt->AnzahlXY;
  }      
  if (XYpos1!=0) 
  {
    if (fabs(searchval-Schnitt->xcamb[XYpos1])>0.001)
    {
      Zval1=(searchval-Schnitt->xcamb[XYpos1-1])*(Schnitt->zcamb[XYpos1]-Schnitt->zcamb[XYpos1-1])/(Schnitt->xcamb[XYpos1]-Schnitt->xcamb[XYpos1-1])+Schnitt->zcamb[XYpos1-1];
      Xval1=searchval;
    }else{
      Zval1=Schnitt->zcamb[XYpos1];
      Xval1=Schnitt->xcamb[XYpos1];
    }
  } else {
    if (fabs(searchval-Schnitt->xcamb[XYpos1])>0.001)
    {
      Zval1=(Schnitt->xcamb[XYpos1]-searchval)*(Schnitt->zcamb[XYpos1+1]-Schnitt->zcamb[XYpos1])/(Schnitt->xcamb[XYpos1+1]-Schnitt->xcamb[XYpos1])+Schnitt->zcamb[XYpos1];
      Xval1=searchval;
    }else{
      Zval1=Schnitt->zcamb[XYpos1];
      Xval1=Schnitt->xcamb[XYpos1];
    }
  }
  if (XYpos2!=1) 
  {
    if (fabs(searchval2-Schnitt->xcamb[XYpos2])>0.001)
    { 
      Zval2=(Schnitt->xcamb[XYpos2]-searchval2)*(Schnitt->zcamb[XYpos2+1]-Schnitt->zcamb[XYpos2])/(Schnitt->xcamb[XYpos2+1]-Schnitt->xcamb[XYpos2])+Schnitt->zcamb[XYpos2];
      Xval2=searchval2;
    }else{
      Zval2=Schnitt->zcamb[XYpos2];
      Xval2=Schnitt->xcamb[XYpos2];
    }
  } else {
    if (fabs(searchval2-Schnitt->xcamb[XYpos2])>0.001)
    {     
      Zval2=(searchval2-Schnitt->xcamb[XYpos2-1])*(Schnitt->zcamb[XYpos2]-Schnitt->zcamb[XYpos2-1])/(Schnitt->xcamb[XYpos2]-Schnitt->xcamb[XYpos2-1])+Schnitt->zcamb[XYpos2-1];    
      Xval2=searchval2;
    }else{
      Zval2=Schnitt->zcamb[XYpos2];
      Xval2=Schnitt->xcamb[XYpos2];
    }
  }
  
  if ((XYpos1+1)!=XYpos2)
  {
    int laufW;
    double TeilerLaenge;
    TeilerLaenge=0.0;
    Winkel=0.0;
    for (laufW=XYpos1;laufW<XYpos2;laufW++)
    {
      Winkel=Winkel+atan((Schnitt->zcamb[laufW+1]-Schnitt->zcamb[laufW])/(Schnitt->xcamb[laufW+1]-Schnitt->xcamb[laufW]))*180/PI*(-1)*(Schnitt->xcamb[laufW+1]-Schnitt->xcamb[laufW]);
      TeilerLaenge=TeilerLaenge+(Schnitt->xcamb[laufW+1]-Schnitt->xcamb[laufW]);
    }
    if (Schnitt->zcamb[XYpos2]!=Zval2)
    {
      Winkel=Winkel+atan((Zval2-Schnitt->zcamb[XYpos2])/(Xval2-Schnitt->xcamb[XYpos2]))*180/PI*(-1)*(Xval2-Schnitt->xcamb[XYpos2]);
      TeilerLaenge=TeilerLaenge+(Xval2-Schnitt->xcamb[XYpos2]);      
    }
    Winkel=Winkel/TeilerLaenge;
  }else{
    Winkel=atan((Zval2-Zval1)/(Xval2-Xval1))*180/PI*(-1);
  }
  Winkel=Winkel+Schnitt->AlphaPlus;
  //printf("Serchval 1&2 Xval1 Zval1 Xval2 Zval2 Winkel: %lf  %lf %lf  %lf %lf  %lf %lf \n",searchval,searchval2, Xval1,Zval1,Xval2,Zval2,Winkel);  
  if ( (Schnitt->KlappeNrI > 0) && (IArand == 1) && ((tiefe+1)>(maxtiefe-Schnitt->anzPanelKlappe)))
  {
    Winkel=Winkel+Klappen[Schnitt->KlappeNrI-1].winkel;
  }
  if ( (Schnitt->KlappeNrA > 0) && (IArand == 2) && ((tiefe+1)>(maxtiefe-Schnitt->anzPanelKlappe)))
  {
    Winkel=Winkel+Klappen[Schnitt->KlappeNrA-1].winkel;
  }  
  return (Winkel);
}


/*Cameber Bestimmung wie sie in der Lifting Line Documentation empfohlen wird*/
double Get_Camber(struct EinSchnitt *Schnitt, int tiefe, int maxtiefe,struct EineKlappen *Klappen, int IArand, double startCamb, double endCamb)
{
  int lauf, XYpos1, XYpos2;
  double searchval, searchval2, Zval1, Zval2, Xval1, Xval2, Winkel;
  
  XYpos1=-1;
  XYpos2=-1;
  if (Schnitt->anzPanelKlappe>0)
  {
    if ((tiefe+1)<=(maxtiefe-Schnitt->anzPanelKlappe))
    { 
      searchval = (double)(1-Schnitt->etaSchnitt)*(tiefe+startCamb) * 1/(maxtiefe-Schnitt->anzPanelKlappe);  
    }else{
      searchval = (double)(1-Schnitt->etaSchnitt)+Schnitt->etaSchnitt/Schnitt->anzPanelKlappe*(tiefe-(maxtiefe-Schnitt->anzPanelKlappe)+startCamb);      
    }
  }else{
    searchval = (double)(tiefe+startCamb) * 1/maxtiefe;  
  }
  for (lauf = 0; lauf <= Schnitt->AnzahlXY; lauf++)
  {
    if ((XYpos1==-1) && (Schnitt->xcamb[lauf]> searchval))
    {
      XYpos1=lauf;
      break;
    }
  }
  if (Schnitt->anzPanelKlappe>0)
  {
    if ((tiefe+1)<=(maxtiefe-Schnitt->anzPanelKlappe))
    {
      searchval2 = (double)(1-Schnitt->etaSchnitt)*(tiefe+endCamb) * 1/(maxtiefe-Schnitt->anzPanelKlappe);  
    }else{
      searchval2 = (double)(1-Schnitt->etaSchnitt)+Schnitt->etaSchnitt/Schnitt->anzPanelKlappe*(tiefe-(maxtiefe-Schnitt->anzPanelKlappe)+endCamb);      
    }
  }else{
    searchval2 = (double)(tiefe+endCamb) * 1/maxtiefe;  
  }
  
  if (searchval2 !=1.0) 
  {
    for (lauf = 0; lauf <= Schnitt->AnzahlXY; lauf++)
    {
      if ((XYpos2==-1) && (Schnitt->xcamb[lauf]> searchval2))
      {
        XYpos2=lauf;
	break;
      }
    }
  } else {
    XYpos2=Schnitt->AnzahlXY;
  }    
  if (XYpos1!=0) 
  {
    if (fabs(searchval-Schnitt->xcamb[XYpos1])>0.001)
    {
      Zval1=(searchval-Schnitt->xcamb[XYpos1-1])*(Schnitt->zcamb[XYpos1]-Schnitt->zcamb[XYpos1-1])/(Schnitt->xcamb[XYpos1]-Schnitt->xcamb[XYpos1-1])+Schnitt->zcamb[XYpos1-1];
      Xval1=searchval;
    }else{
      Zval1=Schnitt->zcamb[XYpos1];
      Xval1=Schnitt->xcamb[XYpos1];
    }
  } else {
    if (fabs(searchval-Schnitt->xcamb[XYpos1])>0.001)
    {
      Zval1=(Schnitt->xcamb[XYpos1]-searchval)*(Schnitt->zcamb[XYpos1+1]-Schnitt->zcamb[XYpos1])/(Schnitt->xcamb[XYpos1+1]-Schnitt->xcamb[XYpos1])+Schnitt->zcamb[XYpos1];
      Xval1=searchval;
    }else{
      Zval1=Schnitt->zcamb[XYpos1];
      Xval1=Schnitt->xcamb[XYpos1];
    }
  }
  if (XYpos2!=1) 
  {
    if (fabs(searchval2-Schnitt->xcamb[XYpos2])>0.001)
    { 
      Zval2=(Schnitt->xcamb[XYpos2]-searchval2)*(Schnitt->zcamb[XYpos2+1]-Schnitt->zcamb[XYpos2])/(Schnitt->xcamb[XYpos2+1]-Schnitt->xcamb[XYpos2])+Schnitt->zcamb[XYpos2];
      Xval2=searchval2;
    }else{
      Zval2=Schnitt->zcamb[XYpos2];
      Xval2=Schnitt->xcamb[XYpos2];
    }
  } else {
    if (fabs(searchval2-Schnitt->xcamb[XYpos2])>0.001)
    {     
      Zval2=(searchval2-Schnitt->xcamb[XYpos2-1])*(Schnitt->zcamb[XYpos2]-Schnitt->zcamb[XYpos2-1])/(Schnitt->xcamb[XYpos2]-Schnitt->xcamb[XYpos2-1])+Schnitt->zcamb[XYpos2-1];    
      Xval2=searchval2;
    }else{
      Zval2=Schnitt->zcamb[XYpos2];
      Xval2=Schnitt->xcamb[XYpos2];
    }
  }
  Winkel=atan((Zval2-Zval1)/(Xval2-Xval1))*180/PI*(-1)+Schnitt->AlphaPlus;
  //printf("Serchval 1&2 Xval1 Zval1 Xval2 Zval2 Winkel: %lf  %lf %lf  %lf %lf  %lf %lf \n",searchval,searchval2, Xval1,Zval1,Xval2,Zval2,Winkel);
  if ( (Schnitt->KlappeNrI > 0) && (IArand == 1) && ((tiefe+1)>(maxtiefe-Schnitt->anzPanelKlappe)))
  {
    Winkel=Winkel+Klappen[Schnitt->KlappeNrI-1].winkel;
  }
  if ( (Schnitt->KlappeNrA > 0) && (IArand == 2) && ((tiefe+1)>(maxtiefe-Schnitt->anzPanelKlappe)))
  {
    Winkel=Winkel+Klappen[Schnitt->KlappeNrA-1].winkel;
  }
  return (Winkel);
}


struct SectionDistribution *getSecInfo(struct Einfluegel *Fluegel, int *NrSect)
{
  int lauf,lauf2,cont;
  double laufS=0;
  struct SectionDistribution *SectDist;

  *NrSect=1;
  SectDist = (struct SectionDistribution *)malloc(sizeof (struct SectionDistribution));
  strcpy(SectDist[0].name, Fluegel->Schnitte[0].ProfilName);
  SectDist[0].NrEtaPos=1;
  SectDist[0].etaPos = (double  *)malloc(sizeof (double));
  SectDist[0].etaPos = 0;
  if (Fluegel->SymExt==1)
  {
    laufS=Fluegel->Schnitte[0].posy;
    SectDist->NrEtaPos=2;
    SectDist->etaPos = (double *)realloc(SectDist->etaPos,SectDist->NrEtaPos*sizeof (double));
    SectDist->etaPos[1] = laufS;
  }

  for (lauf=1;lauf<Fluegel->anzahlSchnitte; lauf++)
  {
    laufS=laufS+sqrt(pow(Fluegel->Schnitte[lauf].posx-Fluegel->Schnitte[lauf-1].posx, 2)+pow(Fluegel->Schnitte[lauf].posy-Fluegel->Schnitte[lauf-1].posy, 2)+pow(Fluegel->Schnitte[lauf].posz-Fluegel->Schnitte[lauf-1].posz, 2));
    cont=1;
    for (lauf2=0; lauf2<*NrSect; lauf2++)
    {
      if (strcmp(Fluegel->Schnitte[lauf].ProfilName, SectDist[lauf2].name)==0)
      {
        cont=0;
        SectDist[lauf2].NrEtaPos=SectDist[lauf2].NrEtaPos+1;
        SectDist[lauf2].etaPos = (double *)realloc(SectDist[lauf2].etaPos,SectDist[lauf2].NrEtaPos*sizeof (double));
        SectDist[lauf2].etaPos[SectDist[lauf2].NrEtaPos-1] = laufS;
        //printf("NrEtas:%d  Cur.SecName: %s    CurEtAPos: %lf\n", SectDist[lauf2].NrEtaPos, Fluegel->Schnitte[lauf].ProfilName,SectDist[lauf2].etaPos[SectDist->NrEtaPos-1]);
      }
    }
    if (cont==1)
    {
      *NrSect=(*NrSect)+1;
      SectDist = (struct SectionDistribution *)realloc(SectDist,(*NrSect)*sizeof (struct SectionDistribution));
      strcpy(SectDist[*NrSect-1].name, Fluegel->Schnitte[lauf].ProfilName);
      SectDist[*NrSect-1].NrEtaPos=1;
      SectDist[*NrSect-1].etaPos = (double *)malloc(sizeof (double));
      *SectDist[*NrSect-1].etaPos = laufS;
      //printf("NrSect:%d  Cur.SecName: %s  CurEtAPos: %lf\n", *NrSect, Fluegel->Schnitte[lauf].ProfilName,*SectDist[*NrSect-1].etaPos);
    }
  }
  for (lauf2=0; lauf2<*NrSect; lauf2++)
  { 
    for (lauf=0; lauf<SectDist[lauf2].NrEtaPos; lauf++)
    {
      SectDist[lauf2].etaPos[lauf]=SectDist[lauf2].etaPos[lauf]/laufS;
    }
  }
  return SectDist;
}


void create_LiLiInp_File (struct inputfile *ParaFile)
{
  FILE *fLinp;
  int lauf1, lauf2, lauf3, lauf4, anzPanel,anzLiliDefPunkte,Kopel,Kopel2;
  double XWert, YWert ,ZWert, Lschnitt, Alfa;
  double DynVisc;
  char HEAD4[25],HEAD6[100];
  
  strcpy(HEAD4," VERSION ");
  strcat(HEAD4,ParaFile->LiLiVersionName);
  strcat(HEAD4," \n");
  strcpy(HEAD6," LIFTING_LINE ");
  strcat(HEAD6,ParaFile->LiLiVersionName);
  strcat(HEAD6," , EXAMPLE: SIMPLE\n");
  
  if (ParaFile->MASSE_SPEED==1)
  {
    if ((ParaFile->POLINT_MASSE<=0)||(ParaFile->RefDensity<=0)||(ParaFile->refAerea<=0))
    {
      errmessage (48);
    }
    if (ParaFile->numberAlpha>0)
    {
      errmessage (49);
    }
    ParaFile->RefSpeed=sqrt((ParaFile->POLINT_MASSE*9.81*2)/(ParaFile->RefDensity*ParaFile->refAerea*ParaFile->targetCl[0]));
  }
  fLinp = fopen("LaufV1.inp","w+");  
  if (fLinp == NULL)
  {
    errmessage (19);
  }
  fputs(HEAD1, fLinp);
  fputs(HEAD2, fLinp);
  fputs(HEAD3, fLinp);
  fputs(HEAD4, fLinp);
  fputs(HEAD5, fLinp);
  fputs(HEAD6, fLinp);
  if (strcmp(ParaFile->LiLiVersion,"V2p3b")==0)
  {
    fputs(HEAD7, fLinp);  
  }else if (strcmp(ParaFile->LiLiVersion,"V2p3bCirc")==0)    
  {
    fputs(HEAD7_KreisFlug, fLinp);
  }else{
    errmessage (63);
  }  
  anzPanel=0;
  anzLiliDefPunkte=0;
  ParaFile->AnzahlSpanPanel=0;
  for (lauf1=0; lauf1< ParaFile->NbWings; lauf1++)
  {
    for(lauf2=0; lauf2< ParaFile->Fluegel[lauf1].anzahlSchnitte; lauf2++)
    {
      if (lauf2==0)
      {
        if (ParaFile->Fluegel[lauf1].SymExt == 1)
	{
	  anzPanel=anzPanel + ParaFile->Fluegel[lauf1].Schnitte[lauf2].AnzahlPan * ParaFile->Fluegel[lauf1].AnzahlTiefe;
	  anzLiliDefPunkte = anzLiliDefPunkte + ParaFile->Fluegel[lauf1].AnzahlTiefe;
	  ParaFile->AnzahlSpanPanel = ParaFile->AnzahlSpanPanel + ParaFile->Fluegel[lauf1].Schnitte[lauf2].AnzahlPan;
	  ParaFile->Fluegel[lauf1].AnzahlSpanPanel = ParaFile->Fluegel[lauf1].Schnitte[lauf2].AnzahlPan;
	} else {
	  ParaFile->Fluegel[lauf1].AnzahlSpanPanel = 0;
	}	
      } else 
      {
        anzPanel=anzPanel + ParaFile->Fluegel[lauf1].Schnitte[lauf2].AnzahlPan * ParaFile->Fluegel[lauf1].AnzahlTiefe;
	anzLiliDefPunkte = anzLiliDefPunkte + ParaFile->Fluegel[lauf1].AnzahlTiefe;
	ParaFile->AnzahlSpanPanel = ParaFile->AnzahlSpanPanel + ParaFile->Fluegel[lauf1].Schnitte[lauf2].AnzahlPan;
	ParaFile->Fluegel[lauf1].AnzahlSpanPanel = ParaFile->Fluegel[lauf1].AnzahlSpanPanel + ParaFile->Fluegel[lauf1].Schnitte[lauf2].AnzahlPan;
      }
    }
  }
  ParaFile->AnzahlPanel=anzPanel; 
  if (strcmp(ParaFile->LiLiVersion,"V2p3b")==0)
  {
    fprintf(fLinp, "    %d       %lf       %lf           %d               %d             %d           14",ParaFile->GlobSym ,ParaFile->MachNumber, GeoScale ,anzLiliDefPunkte, ParaFile->numberAlpha, ZusAufPunkte);    
  }else if (strcmp(ParaFile->LiLiVersion,"V2p3bCirc")==0)
  {
    fprintf(fLinp, "    %d       %lf       %lf           %d               %d             %d           %d             14",ParaFile->GlobSym ,ParaFile->MachNumber, GeoScale ,anzLiliDefPunkte, ParaFile->numberAlpha, ParaFile->KreisflugModus, ZusAufPunkte);        
  }
  if (ParaFile->XMLKOP == 0)
  {
    fprintf(fLinp, "      0");
  } else 
  {
    fprintf(fLinp, "      1");
  }
  if (strcmp(ParaFile->LiLiVersion,"V2p3bCirc")==0)
  {
    fprintf(fLinp, "      0\n");  
  } else {
    fprintf(fLinp, "\n");  
  }   
  fputs(HEAD8, fLinp);
  
  fprintf(fLinp, " %lf  %lf  %lf   %lf  %lf  %lf   %lf  %lf", ParaFile->refAerea, ParaFile->refSpan , ParaFile->refSpan , ParaFile->refChordlength ,ParaFile->refChordlength , ParaFile->origXYZ[0] , ParaFile->origXYZ[1] , ParaFile->origXYZ[2] ); 
  
  if (ParaFile->Alfa_Rot != 0)
  {
    fprintf(fLinp, "  %lf   %lf  %lf\n", ParaFile->origXYZ[0] * cos(ParaFile->Alfa_Rot * 3.141592653/180*-(1))+ ParaFile->origXYZ[2] * sin(ParaFile->Alfa_Rot * 3.141592653/180*-(1)) , ParaFile->origXYZ[1] , ParaFile->origXYZ[0] * sin(ParaFile->Alfa_Rot * 3.141592653/180*-(1))*(-1)+ ParaFile->origXYZ[2] * cos(ParaFile->Alfa_Rot * 3.141592653/180*-(1)) );
  }else{
    fprintf(fLinp, "  %lf   %lf  %lf\n", ParaFile->origXYZ[0] , ParaFile->origXYZ[1] , ParaFile->origXYZ[2] );
  }
  lauf4 = 0; /* lauf4 z�hlt die Anzhal der gesetzten definitionsschitte*/
  /* Schleife �ber die ANzahl der Fl�gel*/
  for (lauf1=0; lauf1<ParaFile->NbWings; lauf1++)
  {
    for (lauf2=0;lauf2<ParaFile->Fluegel[lauf1].AnzahlTiefe; lauf2++)
    {
      for (lauf3=0; lauf3<ParaFile->Fluegel[lauf1].anzahlSchnitte; lauf3++) 
      {
        lauf4++;
        fputs(HEAD9, fLinp);
        Lschnitt=getL(ParaFile,lauf1,lauf2,ParaFile->Fluegel[lauf1].anzahlSchnitte-lauf3-1);
	//Lschnitt=ParaFile->Fluegel[lauf1].Schnitte[ParaFile->Fluegel[lauf1].anzahlSchnitte-lauf3-1].tiefe/ParaFile->Fluegel[lauf1].AnzahlTiefe;
	XYZtransformation(ParaFile,lauf1,lauf2,ParaFile->Fluegel[lauf1].anzahlSchnitte-lauf3-1, &XWert, &YWert,&ZWert);
	fprintf(fLinp, " %d  %d  %lf  %lf  %lf  %lf ", lauf4, ParaFile->Fluegel[lauf1].Schnitte[ParaFile->Fluegel[lauf1].anzahlSchnitte-lauf3-1].AnzahlPan, XWert , YWert , ZWert ,Lschnitt );            
       
        if (lauf3==ParaFile->Fluegel[lauf1].anzahlSchnitte-1)
	{
          XYZtransformation(ParaFile,lauf1,lauf2,ParaFile->Fluegel[lauf1].anzahlSchnitte-lauf3-1, &XWert, &YWert,&ZWert);
	  YWert=0;	
	} else	
        {
          //Lschnitt=ParaFile->Fluegel[lauf1].Schnitte[ParaFile->Fluegel[lauf1].anzahlSchnitte-lauf3-2].tiefe/ParaFile->Fluegel[lauf1].AnzahlTiefe;
          Lschnitt=getL(ParaFile,lauf1,lauf2,ParaFile->Fluegel[lauf1].anzahlSchnitte-lauf3-2);
	  XYZtransformation(ParaFile,lauf1,lauf2,ParaFile->Fluegel[lauf1].anzahlSchnitte-lauf3-2, &XWert, &YWert,&ZWert);
	}
	fprintf(fLinp, " %lf  %lf  %lf  %lf \n",  XWert , YWert , ZWert ,Lschnitt );                    
	fputs(HEAD10, fLinp);
	//Panel-Randbedingung Rand1
	Kopel=IsKoppel(lauf1,lauf2+1,ParaFile->Fluegel[lauf1].anzahlSchnitte-lauf3, ParaFile->AnzKreuz,ParaFile->Kreuz, ParaFile->Fluegel);
        if ((Kopel==0) && (lauf3==0))
        {
          fprintf(fLinp," %s %d 0 ", ARAND , Kopel);
	} else if ((Kopel > 0) ) {
          Kopel2=whichKreuzNR(lauf1,lauf2+1,ParaFile->Fluegel[lauf1].anzahlSchnitte-lauf3, ParaFile->AnzKreuz,ParaFile->Kreuz, ParaFile->Fluegel);
	  if (((ParaFile->Kreuz[Kopel2].KreuzFl1==(lauf1+1)) && (ParaFile->Kreuz[Kopel2].GabBedFl1==2)) || ((ParaFile->Kreuz[Kopel2].KreuzFl2==(lauf1+1)) && (ParaFile->Kreuz[Kopel2].GabBedFl2==2)))
	  {
	    fprintf(fLinp," %s %d 0 ", K2RAND , Kopel);
          }else{
	    fprintf(fLinp," %s %d 0 ", KRAND , Kopel);
	  }
	} else {
          fprintf(fLinp," %s %d 0 ", ZRAND , Kopel); 
	}
        //Panel-Randbedingung Rand2
	Kopel=IsKoppel(lauf1,lauf2+1,ParaFile->Fluegel[lauf1].anzahlSchnitte-lauf3-1, ParaFile->AnzKreuz,ParaFile->Kreuz, ParaFile->Fluegel);
        if ((Kopel==0) && (lauf3 < ParaFile->Fluegel[lauf1].anzahlSchnitte-2))
        {
          fprintf(fLinp," %s %d 0 ", ZRAND , Kopel);
	} else if ((Kopel==0) && (lauf3 == ParaFile->Fluegel[lauf1].anzahlSchnitte-2) && (ParaFile->Fluegel[lauf1].SymExt == 1)) {
	  fprintf(fLinp," %s %d 0 ", ZRAND , Kopel); 	
	} else if ((Kopel>0) && (lauf3 < ParaFile->Fluegel[lauf1].anzahlSchnitte-2)) {
          Kopel2=whichKreuzNR(lauf1,lauf2+1,ParaFile->Fluegel[lauf1].anzahlSchnitte-lauf3, ParaFile->AnzKreuz,ParaFile->Kreuz, ParaFile->Fluegel);
	  if (((ParaFile->Kreuz[Kopel2].KreuzFl1==(lauf1+1)) && (ParaFile->Kreuz[Kopel2].GabBedFl1==1)) || ((ParaFile->Kreuz[Kopel2].KreuzFl2==(lauf1+1)) && (ParaFile->Kreuz[Kopel2].GabBedFl2==1)))
	  {
	    fprintf(fLinp," %s %d 0 ", KRAND , Kopel);
          }else{
	    fprintf(fLinp," %s %d 0 ", K2RAND , Kopel);
	  }
	} else if ((Kopel>0) && (lauf3 == ParaFile->Fluegel[lauf1].anzahlSchnitte-2) && (ParaFile->Fluegel[lauf1].SymExt == 0) && (ParaFile->Fluegel[lauf1].SymRandBed == 0)) {
          Kopel2=whichKreuzNR(lauf1,lauf2+1,ParaFile->Fluegel[lauf1].anzahlSchnitte-lauf3, ParaFile->AnzKreuz,ParaFile->Kreuz, ParaFile->Fluegel);
	  if (((ParaFile->Kreuz[Kopel2].KreuzFl1==(lauf1+1)) && (ParaFile->Kreuz[Kopel2].GabBedFl1==1)) || ((ParaFile->Kreuz[Kopel2].KreuzFl2==(lauf1+1)) && (ParaFile->Kreuz[Kopel2].GabBedFl2==1)))
	  {
	    fprintf(fLinp," %s %d 0 ", KRAND , Kopel);
          }else{
	    fprintf(fLinp," %s %d 0 ", K2RAND , Kopel);
	  }
	} else if ((Kopel>0) && (lauf3 == ParaFile->Fluegel[lauf1].anzahlSchnitte-2) && (ParaFile->Fluegel[lauf1].SymExt == 0) && (ParaFile->Fluegel[lauf1].SymRandBed == 1)) {
          errmessage(20);
        } else if ((Kopel>0) && (lauf3 == ParaFile->Fluegel[lauf1].anzahlSchnitte-2) && (ParaFile->Fluegel[lauf1].SymExt == 0))  {
          Kopel2=whichKreuzNR(lauf1,lauf2+1,ParaFile->Fluegel[lauf1].anzahlSchnitte-lauf3, ParaFile->AnzKreuz,ParaFile->Kreuz, ParaFile->Fluegel);
	  if (((ParaFile->Kreuz[Kopel2].KreuzFl1==(lauf1+1)) && (ParaFile->Kreuz[Kopel2].GabBedFl1==1)) || ((ParaFile->Kreuz[Kopel2].KreuzFl2==(lauf1+1)) && (ParaFile->Kreuz[Kopel2].GabBedFl2==1)))
	  {
	    fprintf(fLinp," %s %d 0 ", KRAND , Kopel);
          }else{
	    fprintf(fLinp," %s %d 0 ", K2RAND , Kopel);
	  }
	} else if ((Kopel==0) && (lauf3 == ParaFile->Fluegel[lauf1].anzahlSchnitte-2) && (ParaFile->Fluegel[lauf1].SymExt == 0))  {
          if (ParaFile->Fluegel[lauf1].SymRandBed == 1)
          {
            fprintf(fLinp," %s %d 0 ", IRAND , Kopel);
	  } else {
            fprintf(fLinp," %s %d 0 ", ARAND , Kopel);
	  }   
	} else if ((lauf3 == ParaFile->Fluegel[lauf1].anzahlSchnitte-1))  {
          if (ParaFile->Fluegel[lauf1].SymRandBed == 1)
          {
            fprintf(fLinp," %s %d 0 ", IRAND , Kopel);
	  } else {
            fprintf(fLinp," %s %d 0 ", ARAND , Kopel);
	  }   
	}
        if (ParaFile->CamberMethod==1)
        {
          Alfa=Get_Camber(&ParaFile->Fluegel[lauf1].Schnitte[ParaFile->Fluegel[lauf1].anzahlSchnitte-lauf3-1],lauf2, ParaFile->Fluegel[lauf1].AnzahlTiefe,ParaFile->Klappen ,1,ParaFile->CamberStart,ParaFile->CamberEnd); // Die 1 zum Schluss heisst InnenRand vom Schnitt (gilt von Symmetrieebene)
	}else{
          Alfa=Get_Camber2(&ParaFile->Fluegel[lauf1].Schnitte[ParaFile->Fluegel[lauf1].anzahlSchnitte-lauf3-1],lauf2, ParaFile->Fluegel[lauf1].AnzahlTiefe,ParaFile->Klappen ,1,ParaFile->CamberStart,ParaFile->CamberEnd); // Die 1 zum Schluss heisst InnenRand vom Schnitt (gilt von Symmetrieebene)        
        }
        fprintf(fLinp," %lf ", Alfa);
        if ((lauf3 == ParaFile->Fluegel[lauf1].anzahlSchnitte-1))
        {
          fprintf(fLinp," %lf ", Alfa+ParaFile->AlfaPlusRumpf);
        }else {
          if (ParaFile->CamberMethod==1)
          {
            Alfa=Get_Camber(&ParaFile->Fluegel[lauf1].Schnitte[ParaFile->Fluegel[lauf1].anzahlSchnitte-lauf3-2],lauf2, ParaFile->Fluegel[lauf1].AnzahlTiefe,ParaFile->Klappen,2,ParaFile->CamberStart,ParaFile->CamberEnd); // Die 2 zum Schluss heisst AussenRand vom Schnitt (gilt von Symmetrieebene)
          }else{
            Alfa=Get_Camber2(&ParaFile->Fluegel[lauf1].Schnitte[ParaFile->Fluegel[lauf1].anzahlSchnitte-lauf3-2],lauf2, ParaFile->Fluegel[lauf1].AnzahlTiefe,ParaFile->Klappen,2,ParaFile->CamberStart,ParaFile->CamberEnd); // Die 2 zum Schluss heisst AussenRand vom Schnitt (gilt von Symmetrieebene)          
          }
          fprintf(fLinp," %lf ", Alfa);
	}
	
	// Im folgenden Block k�nnten noch funktionen hinsichtlich der Allgemeing�ltigkeit erg�nzt werden (Auftribsanstieg, Fl�gel die ein Panel Breit sind)	
        if (lauf3==0) 
        {
          for (lauf3=0; lauf3<ParaFile->Fluegel[lauf1].AnzahlTiefe; lauf3++) 
	  {
	    Kopel=IsKoppel(lauf1,lauf3+1,ParaFile->Fluegel[lauf1].anzahlSchnitte, ParaFile->AnzKreuz,ParaFile->Kreuz, ParaFile->Fluegel);
	    if (Kopel!=0)
	    {
	      break;
	    }
          }
	  if (Kopel!=0)
	  {
	    fprintf(fLinp," 1 0 0 0 \n");
	  }else{
	    fprintf(fLinp," 1 0 1 0 \n");
	  }
	  lauf3=0;
        } else {
          if (lauf3==(ParaFile->Fluegel[lauf1].anzahlSchnitte-2) && ParaFile->GlobSym==0 )
          {
            int zahl;
            for (zahl=0; zahl<ParaFile->Fluegel[lauf1].AnzahlTiefe; zahl++) 
	    {
	      Kopel=IsKoppel(lauf1,zahl+1,(lauf3+1), ParaFile->AnzKreuz,ParaFile->Kreuz, ParaFile->Fluegel);
	      if (Kopel!=0)
	      {
	        break;
	      }
            }
            if (Kopel!=0)
            {
              fprintf(fLinp," 1 0 0 0 \n");
            }else{
              fprintf(fLinp," 1 0 3 0 \n");
            }
          }else{
            //printf("%d %d\n",lauf3,ParaFile->Fluegel[lauf1].anzahlSchnitte );
            fprintf(fLinp," 1 0 0 0 \n");          
          }
        }
	
	if ((lauf3==ParaFile->Fluegel[lauf1].anzahlSchnitte-2) && (ParaFile->Fluegel[lauf1].SymExt == 0))
        {
          lauf3++;
	}
      }
    }
  }
  fputs(HEAD11, fLinp);
  fprintf(fLinp, " %lf ", ParaFile->beta);
  if (ParaFile->numberAlpha > 0)
  {
    for (lauf1=0; lauf1 < ParaFile->numberAlpha ; lauf1++)
    {
      fprintf(fLinp, " %lf ", ParaFile->targetAlpha[lauf1]);
    } 
  }
  fprintf(fLinp, " \n");
  fputs(HEAD12, fLinp);
  fprintf(fLinp, " %d ", ParaFile->numberCL);
  if (ParaFile->numberCL > 0)  
  {
    for (lauf1=0; lauf1 < ParaFile->numberCL ; lauf1++)
    {
      fprintf(fLinp, " %lf ", ParaFile->targetCl[lauf1]);
    } 
  }
  fprintf(fLinp, " \n");
  fputs(HEAD13, fLinp);
  fprintf(fLinp, " %d \n", ParaFile->NbWings);
  for (lauf1=0; lauf1<ParaFile->NbWings; lauf1++)
  {
    fputs(HEAD14, fLinp);
    if (ParaFile->Fluegel[lauf1].SymExt==1)
    {
      fprintf(fLinp, " %d  %d\n", ParaFile->Fluegel[lauf1].anzahlSchnitte, ParaFile->Fluegel[lauf1].AnzahlTiefe );
    } else {
      fprintf(fLinp, " %d  %d\n", ParaFile->Fluegel[lauf1].anzahlSchnitte-1, ParaFile->Fluegel[lauf1].AnzahlTiefe );
    }
  }
  if (ParaFile->XMLKOP==1)
  {
    int NrOfSect;
    struct SectionDistribution *SectDist;
    
    fputs(HEAD16,fLinp);
    if ((ParaFile->RefSpeed==-1) || (ParaFile->RefDensity==-1) || (ParaFile->RefTemp==-1))
    {
      errmessage(23);
    }
    DynVisc=(1.458*0.00001*sqrt(pow(ParaFile->RefTemp,3))/(ParaFile->RefTemp+110.4));
    fprintf(fLinp, "  %lf  %lf   %lf\n",MUE ,ParaFile->RefSpeed ,ParaFile->MachNumber); 
    for (lauf1=1;lauf1<=ParaFile->NbWings;lauf1++)
    {
      SectDist=getSecInfo(&ParaFile->Fluegel[lauf1-1], &NrOfSect);
      fputs(HEAD17,fLinp);
      fprintf(fLinp, " %d \n",NrOfSect);
      for (lauf2=0; lauf2 < NrOfSect; lauf2++)
      {
        fputs(HEAD18,fLinp);
        fprintf(fLinp,"%s\n",SectDist[lauf2].name );
        fputs(HEAD19,fLinp);
        fprintf(fLinp," %d  ",SectDist[lauf2].NrEtaPos);
        for (lauf3=0;lauf3<SectDist[lauf2].NrEtaPos;lauf3++)
        {
          fprintf(fLinp," %lf  ",SectDist[lauf2].etaPos[lauf3]);
	}
        fprintf(fLinp,"\n");
      }
    }
  }
  if (strcmp(ParaFile->LiLiVersion,"V2p3bCirc")==0)
  {
    fputs(KreisflugPara,fLinp);
    ParaFile->RefSpeed=sqrt(pow(ParaFile->RefSpeed/sqrt(cos(ParaFile->HaengeWinkel* 3.141592653/180)),2)+ pow(ParaFile->Steigen,2));
    ParaFile->PhiKreis=ParaFile->RefSpeed/(pow(ParaFile->RefSpeed,2)/(9.81*tan(ParaFile->HaengeWinkel* 3.141592653/180)))*cos(ParaFile->HaengeWinkel* 3.141592653/180)/ParaFile->RefSpeed;    
    fprintf(fLinp,"        %lf         %lf         %lf         %lf                                 0.0000\n",ParaFile->RefSpeed,ParaFile->HaengeWinkel,ParaFile->Steigen,ParaFile->PhiKreis);
  }
  fputs(HEAD20,fLinp);
  fclose(fLinp);
}
  
