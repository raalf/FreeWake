/********************run_STRUCT.c*****************                  
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

#include "run_STRUCT.h"


void run_stru(struct inputfile *input)
{
  int lauf, lauf2 ,lauf3, lauf4, lauf5,lauf6,laufdum,StopSchnitt, AnzP;
  struct file14 resultfile;
  struct PLTdistribution PLTDIST;
  struct BM_FALL *BMDist=NULL;
  double *xOrg,*yOrg,*zOrg,*Syz,*Slaeng;
  double REFSTROE;
  char EXECUTRstring[550] = {"\0"};
  int return_value;
  char ZEILE_READ[500];
  char dummy[100];
  struct STRUCTLOESUNG JIG2FLIGHT;
  struct STRUCTLOESUNG JIG2CURRENT;
  struct Einfluegel *OrgFluegel;
  //struct GESAMT_POLARint *GPolare;
  
  FILE *beamFILE;
  
  get_Vloc(input);
  
  OrgFluegel=(struct Einfluegel *)malloc(input->NbWings*sizeof(struct Einfluegel));
  for (lauf=0; lauf< input->NbWings; lauf++)
  {
    OrgFluegel[lauf].anzahlSchnitte=input->Fluegel[lauf].anzahlSchnitte;
    OrgFluegel[lauf].Schnitte=(struct EinSchnitt  *)malloc(OrgFluegel[lauf].anzahlSchnitte*sizeof(struct EinSchnitt)); 
    for (lauf2=0; lauf2<OrgFluegel[lauf].anzahlSchnitte; lauf2++)
    {
       OrgFluegel[lauf].Schnitte[lauf2].posx=input->Fluegel[lauf].Schnitte[lauf2].posx;
       OrgFluegel[lauf].Schnitte[lauf2].posy=input->Fluegel[lauf].Schnitte[lauf2].posy;
       OrgFluegel[lauf].Schnitte[lauf2].posz=input->Fluegel[lauf].Schnitte[lauf2].posz;
       OrgFluegel[lauf].Schnitte[lauf2].AlphaPlus=input->Fluegel[lauf].Schnitte[lauf2].AlphaPlus;
       OrgFluegel[lauf].Schnitte[lauf2].Vloc=input->Fluegel[lauf].Schnitte[lauf2].Vloc;
       OrgFluegel[lauf].Schnitte[lauf2].tiefe=input->Fluegel[lauf].Schnitte[lauf2].tiefe;
    }
  }
  
  if (input->NumStructDih > 0)
  {
    input->StrDihDist[0].Syz=0;
    for (lauf=1; lauf<input->NumStructDih; lauf++)
    {      
      input->StrDihDist[lauf].Syz=input->StrDihDist[lauf-1].Syz+(input->StrDihDist[lauf].etaPos-input->StrDihDist[lauf-1].etaPos)*input->StructModBezSpan
      /cos(input->StrDihDist[lauf-1].etaV*PI/180);
    }
  }
  
  for (lauf=-1; lauf < input->Numb_Iter; lauf++)
  {
    if (lauf==-1)
    {
      lauf=0;
    }
    if ((input->numberCL+input->numberAlpha) >1)
    {
      errmessage(41);
    }
    run_PROZESS(input);
/*    if (BMDist==NULL)
    {    */
      BMDist=(struct BM_FALL *)realloc(BMDist,(input->numberCL+input->numberAlpha)*sizeof(struct BM_FALL));
/*    }*/
    BM_integration (input, &resultfile, &PLTDIST, BMDist,0);
    write_BM_distribution(BMDist, &resultfile, input, 0); 
    beamFILE=fopen("fxyz4beam.txt","w+");
    if (beamFILE==NULL)
    {
      errmessage(42);
    }
    fprintf(beamFILE,"!    No    x[mm]               y[mm]               z[mm]                 Fx[N]               Fy[N]               Fz[N]                Mx[Nm]              My[Nm]              Mz[Nm]           \n");
    fprintf(beamFILE,"!---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------\n");	
    if (input->ProjektionsFluegel>input->NbWings)
    {
      errmessage(43);
    }
    if (input->ProjektionsFluegel<=0)
    {
      input->ProjektionsFluegel=input->NbWings;
    }    
    lauf4=0; // lauf4 z�hlt die Anzal der St�tzstellen �ber alle Fl�gel mit
    AnzP=0;
    for (lauf2=0; lauf2<input->ProjektionsFluegel; lauf2++)
    {
      if (((lauf2+1)==input->ProjektionsFluegel) && (input->ProjektionsSchnitt< input->Fluegel[lauf2].anzahlSchnitte) && (input->ProjektionsSchnitt>0)) 
      {
        StopSchnitt=input->ProjektionsSchnitt-1;
      }else{
        StopSchnitt=input->Fluegel[lauf2].anzahlSchnitte-1;
      }
      if (lauf==0)
      {
        if (lauf2==0)
        {
          xOrg=(double *)malloc(input->Fluegel[lauf2].AnzahlSpanPanel*sizeof(double));
          yOrg=(double *)malloc(input->Fluegel[lauf2].AnzahlSpanPanel*sizeof(double));
          zOrg=(double *)malloc(input->Fluegel[lauf2].AnzahlSpanPanel*sizeof(double));
          Syz=(double *)malloc(input->Fluegel[lauf2].AnzahlSpanPanel*sizeof(double));
          Slaeng=(double *)malloc(input->Fluegel[lauf2].AnzahlSpanPanel*sizeof(double));
        }else{
          xOrg=(double *)realloc(xOrg,(input->Fluegel[lauf2].AnzahlSpanPanel+lauf4)*sizeof(double));
          yOrg=(double *)realloc(yOrg,(input->Fluegel[lauf2].AnzahlSpanPanel+lauf4)*sizeof(double));
          zOrg=(double *)realloc(zOrg,(input->Fluegel[lauf2].AnzahlSpanPanel+lauf4)*sizeof(double));        
          Syz=(double *)realloc(yOrg,(input->Fluegel[lauf2].AnzahlSpanPanel+lauf4)*sizeof(double));
          Slaeng=(double *)realloc(yOrg,(input->Fluegel[lauf2].AnzahlSpanPanel+lauf4)*sizeof(double));
        }
      }
      
      for (lauf3=0; lauf3<=StopSchnitt; lauf3++)
      {
        if ((lauf3==0) && (input->Fluegel[lauf2].SymExt==0))
        {
          lauf6=0; // lauf6 z�hlt die Anzahl der St�tzstellen im aktuellem Fl�gel mit 
          continue;
        }
        if (lauf3==0)
        {
          lauf6=0; // lauf6 z�hlt die Anzahl der St�tzstellen im aktuellem Fl�gel mit 
        }
        for (lauf5=0; lauf5<input->Fluegel[lauf2].Schnitte[lauf3].AnzahlPan; lauf5++)
        {
          int zwi, laufA;
          double xNEW;
          double yNEW;
          double zNEW;
          double curDhi;
          double Scur;
          zwi=input->Fluegel[lauf2].AnzahlSpanPanel-1;
          if (lauf==0)
          {
            if ((lauf4==0) || (input->NumStructDih==0))
            {
              Syz[lauf4]=0;
              Slaeng[lauf4]=0;
              xNEW=PLTDIST.Fall[0].Fluegel[lauf2].Distri[zwi-lauf6].X+(input->deltaStructX*input->scalfac);
              yNEW=PLTDIST.Fall[0].Fluegel[lauf2].Distri[zwi-lauf6].Y;
              zNEW=PLTDIST.Fall[0].Fluegel[lauf2].Distri[zwi-lauf6].Z+(input->deltaStructZ*input->scalfac);
            } else {
              Syz[lauf4]=Syz[lauf4-1]+sqrt(
              pow(PLTDIST.Fall[0].Fluegel[lauf2].Distri[zwi-lauf6].Y-PLTDIST.Fall[0].Fluegel[lauf2].Distri[zwi-lauf6+1].Y,2)+
              pow(PLTDIST.Fall[0].Fluegel[lauf2].Distri[zwi-lauf6].Z-PLTDIST.Fall[0].Fluegel[lauf2].Distri[zwi-lauf6+1].Z,2));
              Slaeng[lauf4]=Slaeng[lauf4-1]+sqrt(
              pow(PLTDIST.Fall[0].Fluegel[lauf2].Distri[zwi-lauf6].X-PLTDIST.Fall[0].Fluegel[lauf2].Distri[zwi-lauf6+1].X,2)+
              pow(PLTDIST.Fall[0].Fluegel[lauf2].Distri[zwi-lauf6].Y-PLTDIST.Fall[0].Fluegel[lauf2].Distri[zwi-lauf6+1].Y,2)+
              pow(PLTDIST.Fall[0].Fluegel[lauf2].Distri[zwi-lauf6].Z-PLTDIST.Fall[0].Fluegel[lauf2].Distri[zwi-lauf6+1].Z,2));
              for (laufA=0; laufA<input->NumStructDih; laufA++)
              {
                //printf("input->StrDihDist.Syz:%lf  SYZcur:%lf\n",input->StrDihDist[laufA].Syz,Syz[lauf4-1]);
                if (input->StrDihDist[laufA].Syz<=Syz[lauf4-1]*1000)
                {                  
                  curDhi=input->StrDihDist[laufA].etaV;
                }
              }
              Scur=Syz[lauf4]-Syz[lauf4-1];
              //printf("churDhi:%lf   laufA:%d input->NumStructDih:%d\n",curDhi,laufA,input->NumStructDih);
              //printf("deltaZ: %lf  SCUR:%lf zOrg[lauf4-1]:%lf\n\n",sin(curDhi*PI/180)*Scur, Scur,zOrg[lauf4-1]);
              yNEW=yOrg[lauf4-1]+cos(curDhi*PI/180)*Scur;
              zNEW=zOrg[lauf4-1]+sin(curDhi*PI/180)*Scur;
              xNEW=xOrg[lauf4-1]+cos(curDhi*PI/180)*Scur
              *(PLTDIST.Fall[0].Fluegel[lauf2].Distri[zwi-lauf6].X-PLTDIST.Fall[0].Fluegel[lauf2].Distri[zwi-lauf6+1].X)
              /(PLTDIST.Fall[0].Fluegel[lauf2].Distri[zwi-lauf6].Y-PLTDIST.Fall[0].Fluegel[lauf2].Distri[zwi-lauf6+1].Y);
            }
            xOrg[lauf4]=xNEW;
            yOrg[lauf4]=yNEW;
            zOrg[lauf4]=zNEW;
            
          }
          REFSTROE=(input->RefDensity*input->RefSpeed*input->RefSpeed*PLTDIST.Fall[0].Fluegel[lauf2].Distri[zwi-lauf6].REF_AREA)/2;
          /*printf("Refstore: %1.12e  RefDensity: %1.12e  RefSpeed:%1.12e Aerea:%1.12e\n",REFSTROE,input->RefDensity ,input->RefSpeed, PLTDIST.Fall[0].Fluegel[lauf2].Distri[zwi-lauf6].REF_AREA );*/
          fprintf(beamFILE,"   %d  %1.12e  %1.12e  %1.12e  %1.12e  %1.12e  %1.12e  %1.12e  %1.12e  %1.12e\n",lauf4+1 , xOrg[lauf4]/input->scalfac*input->scalefactStruct,  yOrg[lauf4]/input->scalfac*input->scalefactStruct, zOrg[lauf4]/input->scalfac*input->scalefactStruct, 
          PLTDIST.Fall[0].Fluegel[lauf2].Distri[zwi-lauf6].CFX_FROM_CDI*REFSTROE, PLTDIST.Fall[0].Fluegel[lauf2].Distri[zwi-lauf6].CFY*REFSTROE, 
          PLTDIST.Fall[0].Fluegel[lauf2].Distri[zwi-lauf6].CFZ*REFSTROE, 
          PLTDIST.Fall[0].Fluegel[lauf2].Distri[zwi-lauf6].CMX*REFSTROE*PLTDIST.Fall[0].Fluegel[lauf2].Distri[zwi-lauf6].REF_LEN_CMX, 
          PLTDIST.Fall[0].Fluegel[lauf2].Distri[zwi-lauf6].CMY*REFSTROE*PLTDIST.Fall[0].Fluegel[lauf2].Distri[zwi-lauf6].REF_LEN_CMY, 
          PLTDIST.Fall[0].Fluegel[lauf2].Distri[zwi-lauf6].CMZ*REFSTROE*PLTDIST.Fall[0].Fluegel[lauf2].Distri[zwi-lauf6].REF_LEN_CMZ);
          /*printf("   %d %1.12e %1.12e  %1.12e  %1.12e\n  %1.12e  %1.12e  %1.12e \n  %1.12e  %1.12e  %1.12e  \n%1.12e  %1.12e  %1.12e\n\n"
          ,lauf4+1 ,REFSTROE, xOrg[lauf4]/input->scalfac,  yOrg[lauf4]/input->scalfac, zOrg[lauf4]/input->scalfac, 
          PLTDIST.Fall[0].Fluegel[lauf2].Distri[zwi-lauf6].CFX_FROM_CDI, PLTDIST.Fall[0].Fluegel[lauf2].Distri[zwi-lauf6].CFY, 
          PLTDIST.Fall[0].Fluegel[lauf2].Distri[zwi-lauf6].CFZ, 
          PLTDIST.Fall[0].Fluegel[lauf2].Distri[zwi-lauf6].CMX,
          PLTDIST.Fall[0].Fluegel[lauf2].Distri[zwi-lauf6].CMY,
          PLTDIST.Fall[0].Fluegel[lauf2].Distri[zwi-lauf6].CMZ,
          PLTDIST.Fall[0].Fluegel[lauf2].Distri[zwi-lauf6].REF_LEN_CMX, 
          PLTDIST.Fall[0].Fluegel[lauf2].Distri[zwi-lauf6].REF_LEN_CMY, 
          PLTDIST.Fall[0].Fluegel[lauf2].Distri[zwi-lauf6].REF_LEN_CMZ);*/
          lauf4++;
          lauf6++;

          if ((lauf2+1)==input->ProjektionsFluegel)
          {
            AnzP++;
          }
        }
      }              
    }    
    lauf5=input->Fluegel[input->ProjektionsFluegel-1].AnzahlSpanPanel-1;
    if (AnzP<lauf5)
    {
      int laufA;
      double Scur;
      double curDhi;
      double xNEW;
      double yNEW;
      double zNEW;
      if ((lauf==0) && (input->NumStructDih==0))
      {   
        xOrg[lauf4]=BMDist->Fluegel[input->ProjektionsFluegel-1].distri[lauf5-AnzP].X+input->deltaStructX/1000;
        yOrg[lauf4]=BMDist->Fluegel[input->ProjektionsFluegel-1].distri[lauf5-AnzP].Y;
        zOrg[lauf4]=BMDist->Fluegel[input->ProjektionsFluegel-1].distri[lauf5-AnzP].Z+input->deltaStructZ/1000; 
      }else if (lauf==0){
        Syz[lauf4]=Syz[lauf4-1]+sqrt(
        pow(BMDist->Fluegel[input->ProjektionsFluegel-1].distri[lauf5-AnzP].Y-PLTDIST.Fall[0].Fluegel[input->ProjektionsFluegel-1].Distri[lauf5-lauf6+1].Y,2)+
        pow(BMDist->Fluegel[input->ProjektionsFluegel-1].distri[lauf5-AnzP].Z-PLTDIST.Fall[0].Fluegel[input->ProjektionsFluegel-1].Distri[lauf5-lauf6+1].Z,2));
        Slaeng[lauf4]=Slaeng[lauf4-1]+sqrt(
        pow(BMDist->Fluegel[input->ProjektionsFluegel-1].distri[lauf5-AnzP].X-PLTDIST.Fall[0].Fluegel[input->ProjektionsFluegel-1].Distri[lauf5-lauf6+1].X,2)+
        pow(BMDist->Fluegel[input->ProjektionsFluegel-1].distri[lauf5-AnzP].Y-PLTDIST.Fall[0].Fluegel[input->ProjektionsFluegel-1].Distri[lauf5-lauf6+1].Y,2)+
        pow(BMDist->Fluegel[input->ProjektionsFluegel-1].distri[lauf5-AnzP].Z-PLTDIST.Fall[0].Fluegel[input->ProjektionsFluegel-1].Distri[lauf5-lauf6+1].Z,2));
        for (laufA=0; laufA<input->NumStructDih; laufA++)
        {
          if (input->StrDihDist[laufA].Syz<=Syz[lauf4-1])
          {
            curDhi=input->StrDihDist[laufA].etaV;
          }
        }
        Scur=Syz[lauf4]-Syz[lauf4-1];
        //printf("churDhi:%lf\n",curDhi);
        yNEW=yOrg[lauf4-1]+cos(curDhi*PI/180)*Scur;
        zNEW=zOrg[lauf4-1]+sin(curDhi*PI/180)*Scur;
        xNEW=xOrg[lauf4-1]+cos(curDhi*PI/180)*Scur
        *(BMDist->Fluegel[input->ProjektionsFluegel-1].distri[lauf5-AnzP].X-PLTDIST.Fall[0].Fluegel[input->ProjektionsFluegel-1].Distri[lauf5-lauf6+1].X)
        /(BMDist->Fluegel[input->ProjektionsFluegel-1].distri[lauf5-AnzP].Y-PLTDIST.Fall[0].Fluegel[input->ProjektionsFluegel-1].Distri[lauf5-lauf6+1].Y);        
        xOrg[lauf4]=xNEW;
        yOrg[lauf4]=yNEW;
        zOrg[lauf4]=zNEW;
      }
      
      REFSTROE=input->RefDensity/2*input->RefSpeed*input->RefSpeed*input->refAerea;
      fprintf(beamFILE,"   %d  %1.12e  %1.12e  %1.12e  %1.12e  %1.12e  %1.12e  %1.12e  %1.12e  %1.12e\n",lauf4+1 , xOrg[lauf4]* 1000 *input->scalefactStruct,  yOrg[lauf4]* 1000*input->scalefactStruct, zOrg[lauf4]* 1000*input->scalefactStruct, 
      BMDist->Fluegel[input->ProjektionsFluegel-1].distri[lauf5-AnzP].FX*REFSTROE, BMDist->Fluegel[input->ProjektionsFluegel-1].distri[lauf5-AnzP].FY, 
      BMDist->Fluegel[input->ProjektionsFluegel-1].distri[lauf5-AnzP].FZ*REFSTROE, 
      BMDist->Fluegel[input->ProjektionsFluegel-1].distri[lauf5-AnzP].CBM*REFSTROE*input->refSpan, 
      BMDist->Fluegel[input->ProjektionsFluegel-1].distri[lauf5-AnzP].CBT*REFSTROE*input->refSpan, 
      BMDist->Fluegel[input->ProjektionsFluegel-1].distri[lauf5-AnzP].CBG*REFSTROE*input->refSpan);
    }
    fclose(beamFILE); 
    strcpy(EXECUTRstring, "cp fxyz4beam.txt \0");
    strcat(EXECUTRstring, input->StructProg_EXEdir);
    return_value = system(EXECUTRstring); 
    strcpy(EXECUTRstring, "cd \0");
    strcat(EXECUTRstring, input->StructProg_EXEdir);
    return_value = system(EXECUTRstring); 
    strcpy(EXECUTRstring, input->StructProg_EXEdir);
    strcat(EXECUTRstring, input->StructProg_EXEFile);
    return_value = system(EXECUTRstring);
    strcpy(EXECUTRstring, "cp \0");
    strcat(EXECUTRstring, input->StructProg_EXEdir);
    strcat(EXECUTRstring, input->StructProg_ResultSubDir);
    strcat(EXECUTRstring, input->StructProg_ResultFile);
    strcat(EXECUTRstring, " \0");
    strcat(EXECUTRstring, input->LILIDIR);
    strcat(EXECUTRstring, "/struct_beam.dat \0");
    return_value = system(EXECUTRstring);
    strcpy(EXECUTRstring, "cd \0");
    strcat(EXECUTRstring, input->LILIDIR);
    return_value = system(EXECUTRstring);     
    
    
    beamFILE=fopen("struct_beam.dat","r");
    if (beamFILE == NULL)
    {
      printf("\n Konnte Datei struct_beam.dat nicht �ffnen\n");
      errmessage (12);
    }
    fgets (ZEILE_READ,499,beamFILE);
    lauf2=0;
    laufdum=0;
    do
    {
      if (isdigit(ZEILE_READ[lauf2])!=0)
      {
        dummy[laufdum]=ZEILE_READ[lauf2];
        laufdum++;
      }
      lauf2++;      
    } while (ZEILE_READ[lauf2]!='\0');
    sscanf(dummy,"%d",&JIG2CURRENT.AnzDefoPunkte);
    //if (JIG2CURRENT.Punkt==NULL)
    //{
    //  JIG2CURRENT.Punkt=(struct PunktLoesung *)malloc(JIG2CURRENT.AnzDefoPunkte*sizeof(struct PunktLoesung));
    //}
    //printf("JIG2CURRENT.AnzDefoPunkte:%d\n",JIG2CURRENT.AnzDefoPunkte);
    for (lauf2=0; lauf2 < JIG2CURRENT.AnzDefoPunkte; lauf2++)
    {
      fscanf(beamFILE,"%s  %lf %lf %lf",dummy, &JIG2CURRENT.Punkt[lauf2].x, &JIG2CURRENT.Punkt[lauf2].y, &JIG2CURRENT.Punkt[lauf2].z );
      JIG2CURRENT.Punkt[lauf2].x=JIG2CURRENT.Punkt[lauf2].x/input->scalefactStruct;
      JIG2CURRENT.Punkt[lauf2].y=JIG2CURRENT.Punkt[lauf2].y/input->scalefactStruct;
      JIG2CURRENT.Punkt[lauf2].z=JIG2CURRENT.Punkt[lauf2].z/input->scalefactStruct;
      fgets (ZEILE_READ,499,beamFILE);
    }
    for (lauf2=0; lauf2 < JIG2CURRENT.AnzDefoPunkte; lauf2++)
    {
      fscanf(beamFILE,"%s  %lf %lf %lf",dummy, &JIG2CURRENT.Punkt[lauf2].dx, &JIG2CURRENT.Punkt[lauf2].dy, &JIG2CURRENT.Punkt[lauf2].dz );
      JIG2CURRENT.Punkt[lauf2].dx=JIG2CURRENT.Punkt[lauf2].dx/input->scalefactStruct;
      JIG2CURRENT.Punkt[lauf2].dy=JIG2CURRENT.Punkt[lauf2].dy/input->scalefactStruct;
      JIG2CURRENT.Punkt[lauf2].dz=JIG2CURRENT.Punkt[lauf2].dz/input->scalefactStruct;
      fgets (ZEILE_READ,499,beamFILE);
      //printf("%d %s  %lf %lf %lf\n",lauf2,dummy, JIG2CURRENT.Punkt[lauf2].dx, JIG2CURRENT.Punkt[lauf2].dy, JIG2CURRENT.Punkt[lauf2].dz );
    }
    for (lauf2=0; lauf2 < JIG2CURRENT.AnzDefoPunkte; lauf2++)
    {
      fscanf(beamFILE,"%s  %lf %lf %lf",dummy, &JIG2CURRENT.Punkt[lauf2].phix, &JIG2CURRENT.Punkt[lauf2].phiy, &JIG2CURRENT.Punkt[lauf2].phiz );
      fgets (ZEILE_READ,499,beamFILE);
    }
    fclose(beamFILE);
    
    strcpy(EXECUTRstring, "mv \0");
    strcat(EXECUTRstring, input->LILIDIR);
    strcat(EXECUTRstring, "/struct_beam.dat \0");
    strcat(EXECUTRstring, input->LILIDIR);
    strcat(EXECUTRstring, "/LaufV1.lili.");  
    strcat(EXECUTRstring, input->LiLiVersionName);
    strcat(EXECUTRstring, "/ \0");
    return_value = system(EXECUTRstring);
    strcpy(EXECUTRstring, "mv \0");
    strcat(EXECUTRstring, input->LILIDIR);
    strcat(EXECUTRstring, "/fxyz4beam.txt \0");
    strcat(EXECUTRstring, input->LILIDIR);
    strcat(EXECUTRstring, "/LaufV1.lili.");  
    strcat(EXECUTRstring, input->LiLiVersionName);
    strcat(EXECUTRstring, "/ \0");
    return_value = system(EXECUTRstring);
    
    //printf(" lauf:%d input->useJIG2FLIGHT:%d\n\n",lauf,input->useJIG2FLIGHT);
    if ((lauf==0) && (input->useJIG2FLIGHT==1))
    {
      printf("Bin Im Jig 2 Flight_lesen\n");
      beamFILE=fopen(input->StructProg_Jig2FlightFile,"r");
      if (beamFILE==NULL)
      {
        errmessage (45);
      }      
      fgets (ZEILE_READ,499,beamFILE);
      lauf2=0;
      laufdum=0;
      do
      {
        if (isdigit(ZEILE_READ[lauf2])!=0)
        {
          dummy[laufdum]=ZEILE_READ[lauf2];
          laufdum++;
        }
        lauf2++;      
      } while (ZEILE_READ[lauf2]!='\0');
      sscanf(dummy,"%d",&JIG2FLIGHT.AnzDefoPunkte);
      //JIG2FLIGHT.Punkt=(struct PunktLoesung *)malloc(JIG2FLIGHT.AnzDefoPunkte*sizeof(struct PunktLoesung));
      for (lauf2=0; lauf2 < JIG2FLIGHT.AnzDefoPunkte; lauf2++)
      {        
        fscanf(beamFILE,"%s  %lf %lf %lf",dummy, &JIG2FLIGHT.Punkt[lauf2].x, &JIG2FLIGHT.Punkt[lauf2].y, &JIG2FLIGHT.Punkt[lauf2].z );
        fgets (ZEILE_READ,499,beamFILE);
      }
      for (lauf2=0; lauf2 < JIG2FLIGHT.AnzDefoPunkte; lauf2++)
      {
        fscanf(beamFILE,"%s  %lf %lf %lf",dummy, &JIG2FLIGHT.Punkt[lauf2].dx, &JIG2FLIGHT.Punkt[lauf2].dy, &JIG2FLIGHT.Punkt[lauf2].dz );
        fgets (ZEILE_READ,499,beamFILE);
      }
      for (lauf2=0; lauf2 < JIG2FLIGHT.AnzDefoPunkte; lauf2++)
      {
        fscanf(beamFILE,"%s  %lf %lf %lf",dummy, &JIG2FLIGHT.Punkt[lauf2].phix, &JIG2FLIGHT.Punkt[lauf2].phiy, &JIG2FLIGHT.Punkt[lauf2].phiz );
        fgets (ZEILE_READ,499,beamFILE);
      }
      fclose(beamFILE);
      if (JIG2FLIGHT.AnzDefoPunkte!=JIG2CURRENT.AnzDefoPunkte)
      {
        printf("Fehler1 JIG2FLIGHT.AnzDefoPunkte:%d  JIG2CURRENT.AnzDefoPunkte:%d\n",JIG2FLIGHT.AnzDefoPunkte,JIG2CURRENT.AnzDefoPunkte);
        errmessage(44);       
      }
      for (lauf2=0; lauf2 < JIG2FLIGHT.AnzDefoPunkte; lauf2++)
      {
        if ((JIG2FLIGHT.Punkt[lauf2].x != JIG2CURRENT.Punkt[lauf2].x) || (JIG2FLIGHT.Punkt[lauf2].y != JIG2CURRENT.Punkt[lauf2].y) ||(JIG2FLIGHT.Punkt[lauf2].z != JIG2CURRENT.Punkt[lauf2].z))
        {
          printf("Fehler2 JIG2FLIGHT.Punkt[lauf2].x: %lf  JIG2CURRENT.Punkt[lauf2].x: %lf \n JIG2FLIGHT.Punkt[lauf2].y: %lf  JIG2CURRENT.Punkt[lauf2].y: %lf\n JIG2FLIGHT.Punkt[lauf2].z: %lf  IG2CURRENT.Punkt[lauf2].z: %lf \n"
          ,JIG2FLIGHT.Punkt[lauf2].x ,JIG2CURRENT.Punkt[lauf2].x,JIG2FLIGHT.Punkt[lauf2].y ,JIG2CURRENT.Punkt[lauf2].y,
          JIG2FLIGHT.Punkt[lauf2].z ,JIG2CURRENT.Punkt[lauf2].z);
          errmessage(44);
        }
      }
    }
    if ((lauf==0) && (input->useJIG2FLIGHT==0))
    {
      JIG2FLIGHT.AnzDefoPunkte=JIG2CURRENT.AnzDefoPunkte;
      //JIG2FLIGHT.Punkt=(struct PunktLoesung *)malloc(JIG2FLIGHT.AnzDefoPunkte*sizeof(struct PunktLoesung));
      for (lauf2=0; lauf2 < JIG2FLIGHT.AnzDefoPunkte; lauf2++)
      {
        JIG2FLIGHT.Punkt[lauf2].x=0;
        JIG2FLIGHT.Punkt[lauf2].y=0;
        JIG2FLIGHT.Punkt[lauf2].z=0;
        JIG2FLIGHT.Punkt[lauf2].dx=0;
        JIG2FLIGHT.Punkt[lauf2].dy=0;
        JIG2FLIGHT.Punkt[lauf2].dz=0;
        JIG2FLIGHT.Punkt[lauf2].phix=0;
        JIG2FLIGHT.Punkt[lauf2].phiy=0;
        JIG2FLIGHT.Punkt[lauf2].phiz=0;
      }
    }
    printf ("input->Log_Iter_ON:%d\n",input->Log_Iter_ON);
    if (input->Log_Iter_ON == 1)
    {
      double dzahl;
      char laufNR[4];
      strcpy(EXECUTRstring, "mv LaufV1.lili."); 
      strcat(EXECUTRstring, input->LiLiVersionName);
      strcat(EXECUTRstring, " LaufV1.lili.");
      strcat(EXECUTRstring, input->LiLiVersionName);
      strcat(EXECUTRstring, "_\0");
      
      if (lauf>=100)
      {
        errmessage(46);
      }
      if (lauf>=10)
      {      
        double *zahl;
        *zahl=10;
        dzahl=modf(lauf, zahl);
        laufNR[0]=48+*zahl;
        laufNR[1]= 48 + fmod(lauf, *zahl*10);
        laufNR[2]='\0';
      } else {
        laufNR[0]=48+lauf;
        laufNR[1]='\0';
      }     
      strcat(EXECUTRstring, laufNR);
      return_value = system(EXECUTRstring); 
      printf ("RETURN_VALUE : %d ",return_value);
    }
    for (lauf2=0; lauf2< input->NbWings; lauf2++)
    {
      double deltaX, deltaY, deltaZ, deltaPHIX, deltaPHIY, deltaPHIZ,abstX, abstY, abstZ;
      double deltaX2, deltaY2, deltaZ2, deltaPHIX2, deltaPHIY2, deltaPHIZ2,abstX2, abstY2, abstZ2;
      double abstS1 ,abstS2;
      int lauf11, lauf12,dPzaehler;   

      for (lauf3=0; lauf3< input->Fluegel[lauf2].anzahlSchnitte; lauf3++)
      {
        //printf("MAL schauen wie oft ich hier ankomme %d!!\n", lauf3);
        if ((lauf2 >= input->ProjektionsFluegel) || (((lauf2+1) == input->ProjektionsFluegel) && ((lauf3+1)>=input->ProjektionsSchnitt)))
        {          
          //printf("AUSSENfluegel\n");
          //deltaX=JIG2CURRENT.Punkt[JIG2FLIGHT.AnzDefoPunkte-1].dx-JIG2FLIGHT.Punkt[JIG2FLIGHT.AnzDefoPunkte-1].dx;
          deltaX=JIG2CURRENT.Punkt[JIG2FLIGHT.AnzDefoPunkte-1].dx;
          deltaY=JIG2CURRENT.Punkt[JIG2FLIGHT.AnzDefoPunkte-1].dy-JIG2FLIGHT.Punkt[JIG2FLIGHT.AnzDefoPunkte-1].dy;
          deltaZ=JIG2CURRENT.Punkt[JIG2FLIGHT.AnzDefoPunkte-1].dz-JIG2FLIGHT.Punkt[JIG2FLIGHT.AnzDefoPunkte-1].dz;
          deltaPHIX=JIG2CURRENT.Punkt[JIG2FLIGHT.AnzDefoPunkte-1].phix-JIG2FLIGHT.Punkt[JIG2FLIGHT.AnzDefoPunkte-1].phix;
          deltaPHIY=JIG2CURRENT.Punkt[JIG2FLIGHT.AnzDefoPunkte-1].phiy-JIG2FLIGHT.Punkt[JIG2FLIGHT.AnzDefoPunkte-1].phiy;
          deltaPHIZ=JIG2CURRENT.Punkt[JIG2FLIGHT.AnzDefoPunkte-1].phiz-JIG2FLIGHT.Punkt[JIG2FLIGHT.AnzDefoPunkte-1].phiz;
          
          abstX=OrgFluegel[lauf2].Schnitte[lauf3].posx-JIG2CURRENT.Punkt[JIG2FLIGHT.AnzDefoPunkte-1].x;
          abstY=OrgFluegel[lauf2].Schnitte[lauf3].posy-JIG2CURRENT.Punkt[JIG2FLIGHT.AnzDefoPunkte-1].y;
          abstZ=OrgFluegel[lauf2].Schnitte[lauf3].posz-JIG2CURRENT.Punkt[JIG2FLIGHT.AnzDefoPunkte-1].z;
          
          //Die Postions�nderung durch den AbstandX um den Drehpunkt muss zu null gesetzt werden, um bei Mehrfachfl�gel
          //, welche nicht den gleichen Nasenpunkt haben, keine L�cken durch Drehung zu erzeugen.
          abstX=0;
          
          input->Fluegel[lauf2].Schnitte[lauf3].posx=OrgFluegel[lauf2].Schnitte[lauf3].posx+deltaX
          +(abstX*cos(deltaPHIY)+abstZ*sin(deltaPHIY)-abstX) + (abstX*cos(deltaPHIZ)-abstY*sin(deltaPHIZ)-abstX);          
          
          input->Fluegel[lauf2].Schnitte[lauf3].posy=OrgFluegel[lauf2].Schnitte[lauf3].posy+deltaY
          +(abstY*cos(deltaPHIX)-abstZ*sin(deltaPHIX)-abstY)+(abstY*cos(deltaPHIZ)+abstX*sin(deltaPHIZ)-abstY);        
          
          input->Fluegel[lauf2].Schnitte[lauf3].posz=OrgFluegel[lauf2].Schnitte[lauf3].posz+deltaZ
          +(abstZ*cos(deltaPHIX)+abstY*sin(deltaPHIX)-abstZ)+(abstZ*cos(deltaPHIY)-abstX*sin(deltaPHIY)-abstZ);
           
          input->Fluegel[lauf2].Schnitte[lauf3].AlphaPlus=OrgFluegel[lauf2].Schnitte[lauf3].AlphaPlus+cos(OrgFluegel[lauf2].Schnitte[lauf3].Vloc/180*PI)*deltaPHIY*180/PI
          +sin(OrgFluegel[lauf2].Schnitte[lauf3].Vloc/180*PI)*deltaPHIZ*180/PI;
          
          input->Fluegel[lauf2].Schnitte[lauf3].Vloc=OrgFluegel[lauf2].Schnitte[lauf3].Vloc+deltaPHIX*180/PI;
        }else{
          //printf("Innenfluegel\n");
          dPzaehler=0;       
          for (lauf11=0;lauf11<=lauf2; lauf11++)
          {
            for (lauf12=0;lauf12<=lauf3; lauf12++)
            {
              if (lauf12==0)
              {
                if (input->Fluegel[lauf11].SymExt==1)
                {
                  dPzaehler=dPzaehler+input->Fluegel[lauf11].Schnitte[lauf12].AnzahlPan;
                }
              }else{
                dPzaehler=dPzaehler+input->Fluegel[lauf11].Schnitte[lauf12].AnzahlPan;
              }
            }
          }
          deltaX=JIG2CURRENT.Punkt[dPzaehler-1].dx-JIG2FLIGHT.Punkt[dPzaehler-1].dx;
          deltaY=JIG2CURRENT.Punkt[dPzaehler-1].dy-JIG2FLIGHT.Punkt[dPzaehler-1].dy;
          deltaZ=JIG2CURRENT.Punkt[dPzaehler-1].dz-JIG2FLIGHT.Punkt[dPzaehler-1].dz;
          deltaPHIX=JIG2CURRENT.Punkt[dPzaehler-1].phix-JIG2FLIGHT.Punkt[dPzaehler-1].phix;
          deltaPHIY=JIG2CURRENT.Punkt[dPzaehler-1].phiy-JIG2FLIGHT.Punkt[dPzaehler-1].phiy;
          deltaPHIZ=JIG2CURRENT.Punkt[dPzaehler-1].phiz-JIG2FLIGHT.Punkt[dPzaehler-1].phiz;
          abstX=OrgFluegel[lauf2].Schnitte[lauf3].posx-JIG2CURRENT.Punkt[dPzaehler-1].x;
          abstY=OrgFluegel[lauf2].Schnitte[lauf3].posy-JIG2CURRENT.Punkt[dPzaehler-1].y;
          abstZ=OrgFluegel[lauf2].Schnitte[lauf3].posz-JIG2CURRENT.Punkt[dPzaehler-1].z;
          abstX2=abs(OrgFluegel[lauf2].Schnitte[lauf3].posx-JIG2CURRENT.Punkt[dPzaehler].x);
          abstY2=abs(OrgFluegel[lauf2].Schnitte[lauf3].posy-JIG2CURRENT.Punkt[dPzaehler].y);
          abstZ2=abs(OrgFluegel[lauf2].Schnitte[lauf3].posz-JIG2CURRENT.Punkt[dPzaehler].z);
          deltaX2=JIG2CURRENT.Punkt[dPzaehler].dx-JIG2FLIGHT.Punkt[dPzaehler].dx;
          deltaY2=JIG2CURRENT.Punkt[dPzaehler].dy-JIG2FLIGHT.Punkt[dPzaehler].dy;
          deltaZ2=JIG2CURRENT.Punkt[dPzaehler].dz-JIG2FLIGHT.Punkt[dPzaehler].dz;
          deltaPHIX2=JIG2CURRENT.Punkt[dPzaehler].phix-JIG2FLIGHT.Punkt[dPzaehler].phix;
          deltaPHIY2=JIG2CURRENT.Punkt[dPzaehler].phiy-JIG2FLIGHT.Punkt[dPzaehler].phiy;
          deltaPHIZ2=JIG2CURRENT.Punkt[dPzaehler].phiz-JIG2FLIGHT.Punkt[dPzaehler].phiz;
          /*printf("deltaX:%lf  deltaY:%lf   deltaZ:%lf\n",deltaX,deltaY,deltaZ);
          printf("deltaPHIX:%lf  deltaPHIY:%lf   deltaPHIZ:%lf\n",deltaPHIX,deltaPHIY,deltaPHIZ);*/
          abstS1=sqrt(abstX*abstX+abstY*abstY+abstZ*abstZ);
          abstS2=sqrt(abstX2*abstX2+abstY2*abstY2+abstZ2*abstZ2);
          deltaX=(deltaX*abstS2+deltaX2*abstS1)/(abstS1+abstS2);
          deltaY=(deltaY*abstS2+deltaY2*abstS1)/(abstS1+abstS2);
          deltaZ=(deltaZ*abstS2+deltaZ2*abstS1)/(abstS1+abstS2);
          deltaPHIX=(deltaPHIX*abstS2+deltaPHIX2*abstS1)/(abstS1+abstS2);         //+(abstY*cos(deltaPHIX)-abstZ*sin(deltaPHIX)-abstY)+(abstY*cos(deltaPHIZ)+0*sin(deltaPHIZ)-abstY);
          deltaPHIY=(deltaPHIY*abstS2+deltaPHIY2*abstS1)/(abstS1+abstS2);
          deltaPHIZ=(deltaPHIZ*abstS2+deltaPHIZ2*abstS1)/(abstS1+abstS2);
          /*printf("deltaX:%lf  deltaY:%lf   deltaZ:%lf\n",deltaX,deltaY,deltaZ);
          printf("deltaPHIX:%lf  deltaPHIY:%lf   deltaPHIZ:%lf\n",deltaPHIX,deltaPHIY,deltaPHIZ);*/
          input->Fluegel[lauf2].Schnitte[lauf3].posx=OrgFluegel[lauf2].Schnitte[lauf3].posx+deltaX;
          input->Fluegel[lauf2].Schnitte[lauf3].posy=OrgFluegel[lauf2].Schnitte[lauf3].posy+deltaY;
          input->Fluegel[lauf2].Schnitte[lauf3].posz=OrgFluegel[lauf2].Schnitte[lauf3].posz+deltaZ;
          input->Fluegel[lauf2].Schnitte[lauf3].AlphaPlus=OrgFluegel[lauf2].Schnitte[lauf3].AlphaPlus+cos(OrgFluegel[lauf2].Schnitte[lauf3].Vloc*PI/180)*deltaPHIY*180/PI
          +sin(OrgFluegel[lauf2].Schnitte[lauf3].Vloc*PI/180)*deltaPHIZ*180/PI;
          input->Fluegel[lauf2].Schnitte[lauf3].Vloc=OrgFluegel[lauf2].Schnitte[lauf3].Vloc+deltaPHIX*180/PI;

        }
      }
    } 
  }
  if (input->Numb_Iter > 0)
  {
    if ((input->numberCL+input->numberAlpha) >1)
    {
      errmessage(41);
    }
    run_PROZESS(input);
    if (input->Log_Iter_ON == 1)
    {
      strcpy(EXECUTRstring, "mv LaufV1.lili."); 
      strcat(EXECUTRstring, input->LiLiVersionName);
      strcat(EXECUTRstring, "_* LaufV1.lili.");
      strcat(EXECUTRstring, input->LiLiVersionName);
      strcat(EXECUTRstring, "\0");     
      return_value = system(EXECUTRstring); 
      printf ("RETURN_VALUE : %d ",return_value);
    }
  }
}
