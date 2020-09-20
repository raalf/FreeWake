/*****************make_bsurf.c********************                  
*                                                *      
*      Vorentwurfs-Generator "VoGen"             *          
*      Version 1.16                              *          
*      Date of last modification 21.04.2018      *          
*                                                *          
*      Writen by: Jan Himisch                    *      
*             jan.himisch@dlr.de                 *      
*                                                *      
*                                                *
*************************************************/
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "make_bsurf.h"
#include "errmes.h"

void make_bsurf_inp( struct inputfile *input)
{
  int lauf, lauf2, lauf3, check;
  double LY, LZ, LG;
  char DATNAME[50], OUTNAME[20], dummy[4];
  FILE *fBsurf;

  strcpy(DATNAME,"BsurInpCoord.dat");
  fBsurf=fopen(DATNAME,"w+");
  
  fprintf(fBsurf,"alpha_angle !!Achtung Hier noch Winkel aus TAU-Rechnung eintragen!!\n");
  fprintf(fBsurf,"set Span_half  %lf\n", input->refSpan);
  fprintf(fBsurf,"set Cref %lf\n", input->refChordlength);
  fprintf(fBsurf,"set aref %lf\n", input->refAerea);
  
  for (lauf=0; lauf< input->NbWings; lauf++)
  {    
    
    strcpy(OUTNAME,"OutDistFluegel_");
    sprintf(dummy,"%d", lauf+1);
    strcat(OUTNAME, dummy);  
     
    
    for (lauf2=1; lauf2<=input->NUM_OF_ZONES; lauf2++)
    {
      check=0;
      for (lauf3=1; lauf3<=input->BSURF_FLUEGEL[lauf].Zone[0];lauf3++)
      {
        if (input->BSURF_FLUEGEL[lauf].Zone[lauf3]==lauf2)
	{
	  check=1;
	}
      }
      if (lauf>0)
      {
	for (lauf3=1; lauf3<=input->BSURF_FLUEGEL[lauf-1].Zone[0];lauf3++)
	{
	  if ((input->BSURF_FLUEGEL[lauf-1].Zone[lauf3]==lauf2) && (check==1))
	  {
	    check=2;
	  }
	  if ((input->BSURF_FLUEGEL[lauf-1].Zone[lauf3]==lauf2) && (check==0))   
	  {
	    check=3;
	  }
	}
	if (check==1)
	{
	  fprintf(fBsurf,"visable_index %d\n", lauf2);
	}
      }
      if (((lauf==0) && (check==0)) || (check==3))
      {
        fprintf(fBsurf,"disable_index %d\n", lauf2);
      }
    }
    fprintf(fBsurf,"\n\n");
    
    for (lauf2=0; lauf2<(input->Fluegel[lauf].anzahlSchnitte-1); lauf2++)
    {
      fprintf(fBsurf,"eta_x_min %lf \n", input->Fluegel[lauf].Schnitte[lauf2].posx+input->Fluegel[lauf].Schnitte[lauf2].tiefe/4);
      fprintf(fBsurf,"eta_y_min %lf \n", input->Fluegel[lauf].Schnitte[lauf2].posy);
      fprintf(fBsurf,"eta_z_min %lf \n", input->Fluegel[lauf].Schnitte[lauf2].posz);
      fprintf(fBsurf,"eta_x_max %lf \n", input->Fluegel[lauf].Schnitte[lauf2+1].posx+input->Fluegel[lauf].Schnitte[lauf2].tiefe/4);
      fprintf(fBsurf,"eta_y_max %lf \n", input->Fluegel[lauf].Schnitte[lauf2+1].posy);
      fprintf(fBsurf,"eta_z_max %lf \n", input->Fluegel[lauf].Schnitte[lauf2+1].posz);
      LY=input->Fluegel[lauf].Schnitte[lauf2+1].posy-input->Fluegel[lauf].Schnitte[lauf2].posy;
      LZ=input->Fluegel[lauf].Schnitte[lauf2+1].posz-input->Fluegel[lauf].Schnitte[lauf2].posz;
      LG=sqrt(pow(LY,2)+pow(LZ,2));
      fprintf(fBsurf,"cutplane normal 0 %lf %lf\n", LY/LG, LZ/LG);
      fprintf(fBsurf,"cutplane normal2 0 %lf %lf\n\n", LY/LG, LZ/LG);
      for (lauf3=0; lauf3 < (input->Fluegel[lauf].Schnitte[lauf2+1].AnzahlPan*2); lauf3++)
      {
        LG=(double)1/(input->Fluegel[lauf].Schnitte[lauf2+1].AnzahlPan*2);
        LG=LG*(lauf3+1);
        fprintf(fBsurf,"cutplane eta_dist %f  %s\n", LG,OUTNAME); 
      }
      fprintf(fBsurf,"\n\n");   
    }
    strcat(OUTNAME, ".tec"); 
    fprintf(fBsurf,"add_zones \"%s\" \n\n", OUTNAME);    
  }
  fprintf(fBsurf,"\nexit\n");  
  fclose(fBsurf);
}
void get_dist(double *dist1,double *dist2,int FLNR, int laufPLT,int lauf,int  lauf2, struct PLTdistribution *PLTdist, struct BSURF_DATEN *BSURFout)
{
  *dist1=sqrt(pow((PLTdist->Fall[0].Fluegel[FLNR].Distri[laufPLT].Y-BSURFout[lauf].Ypos),2) + pow((PLTdist->Fall[0].Fluegel[FLNR].Distri[laufPLT].Z-BSURFout[lauf].Zpos),2));
  *dist2=sqrt(pow((PLTdist->Fall[0].Fluegel[FLNR].Distri[laufPLT].Y-BSURFout[lauf2].Ypos),2) + pow((PLTdist->Fall[0].Fluegel[FLNR].Distri[laufPLT].Z-BSURFout[lauf2].Zpos),2));
}

void read_BSURF(struct inputfile *input, struct PLTdistribution *PLTdist)
{
  int lauf, lauf2, lauf3, lauf4, check,PLTpos;
  /*lauf: läuft über anzahl der Flügel;  
  lauf2: läuft über Anzahl der Bsurf-Linien: Nach einlesen Anzahl der Bsurf Linien
  lauf3: läuft über anzahl der Lifting-LinePanel
  lauf4: sucht die Bsurf Position die am nächsten zurm aktuellen LiLi-Panel-liegt und speichert diese in check*/
  double dummy1 , dummy2, dummy3, dummy4, dummy5, dummy6, dummy7, dummy8, dummy9, dummy10, dummy11, dummy12, dummy13, dummy14, dummy15, dummy16, dummy17, dummy18, dummy19, dummy20;
  double dummy21, dummy22, dummy23, dummy24, dummy25, dummy26, dummy27, dummy28, dummy29, dummy30, dummy31, dummy32, dummy33, dummy34, dummy35, dummy36, dummy37, dummy38, dummy39, dummy40;
  double dist1, dist2;
  char DATNAME[50], dummy[4],ZEILE_READ[MaxZeilenLaenge]; ;
  FILE *fBsurf;
  struct BSURF_DATEN *BSURFout;
  
  if (PLTdist->Fall==NULL)
  {
    PLTdist->Fall=(struct EinFall *)malloc(1*sizeof(struct EinFall));
    PLTpos=0;
  }else{
    PLTdist->Fall=(struct EinFall *)realloc(PLTdist->Fall,(PLTdist->AnzFaelle+1)*sizeof(struct EinFall));
    PLTpos=PLTdist->AnzFaelle;    
  }
  PLTdist->Fall[PLTpos].Fluegel=(struct FLUEGEL *)malloc(input->NbWings*sizeof(struct FLUEGEL));
  for (lauf=0; lauf< input->NbWings; lauf++)
  {    
    PLTdist->Fall[PLTpos].Fluegel[lauf].Distri=(struct PltDist *)malloc(input->Fluegel[lauf].AnzahlSpanPanel*sizeof(struct PltDist));
    strcpy(DATNAME,"OutDistFluegel_");
    sprintf(dummy,"%d", lauf+1);
    strcat(DATNAME, dummy);  
    strcat(DATNAME, ".load");  
    fBsurf=fopen(DATNAME,"r");
    if (fBsurf ==NULL)
    {
      errmessage(47);
    } 
    for (lauf2=0; lauf2<= 16; lauf2++)
    {
      fgets (ZEILE_READ,MaxZeilenLaenge,fBsurf);
    }
    lauf2=1;
    BSURFout=(struct BSURF_DATEN *)malloc(lauf2*sizeof(struct BSURF_DATEN));
    do
    {
      fgets (ZEILE_READ,MaxZeilenLaenge,fBsurf); 
      check=sscanf(ZEILE_READ," %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf ", &dummy1
      , &dummy2, &dummy3, &dummy4, &dummy5, &dummy6, &dummy7, &dummy8, &dummy9, &dummy10, &dummy11, &dummy12, &dummy13, &dummy14, &dummy15, &dummy16, &dummy17, &dummy18, &dummy19, &dummy20, &dummy21,
      &dummy22, &dummy23, &dummy24, &dummy25, &dummy26, &dummy27, &dummy28, &dummy29, &dummy30, &dummy31, &dummy32, &dummy33, &dummy34, &dummy35, &dummy36, &dummy37, &dummy38, &dummy39, &dummy40);  
      if (check == EOF)
      {
        errmessage(31);
      }
      BSURFout[lauf2-1].Xpos=dummy19;
      BSURFout[lauf2-1].Ypos=dummy20;
      BSURFout[lauf2-1].Zpos=dummy21;
      BSURFout[lauf2-1].CL=dummy3;
      BSURFout[lauf2-1].CD=dummy4;
      BSURFout[lauf2-1].CS=dummy6;
      BSURFout[lauf2-1].CN=dummy35;
      if (!feof(fBsurf))
      {
        lauf2++;
	BSURFout=(struct BSURF_DATEN *)realloc(BSURFout, lauf2*sizeof(struct BSURF_DATEN));
      }
    }
    while  (!feof(fBsurf)); 
    fclose (fBsurf);    
    for (lauf3=0; lauf3< input->Fluegel[lauf].AnzahlSpanPanel; lauf3++)
    {
      int VorZ, VorY;
      dummy1=9999;
      /*Init PLT*/
      PLTdist->Fall[PLTpos].Fluegel[lauf].Distri[lauf3].X=PLTdist->Fall[0].Fluegel[lauf].Distri[lauf3].X;
      PLTdist->Fall[PLTpos].Fluegel[lauf].Distri[lauf3].Y=PLTdist->Fall[0].Fluegel[lauf].Distri[lauf3].Y;
      PLTdist->Fall[PLTpos].Fluegel[lauf].Distri[lauf3].Z=PLTdist->Fall[0].Fluegel[lauf].Distri[lauf3].Z;
      PLTdist->Fall[PLTpos].Fluegel[lauf].Distri[lauf3].S=PLTdist->Fall[0].Fluegel[lauf].Distri[lauf3].S;
      PLTdist->Fall[PLTpos].Fluegel[lauf].Distri[lauf3].LOCAL_CIRCULATION=PLTdist->Fall[0].Fluegel[lauf].Distri[lauf3].LOCAL_CIRCULATION;
      PLTdist->Fall[PLTpos].Fluegel[lauf].Distri[lauf3].CFX_FROM_CDI=PLTdist->Fall[0].Fluegel[lauf].Distri[lauf3].CFX_FROM_CDI;
      PLTdist->Fall[PLTpos].Fluegel[lauf].Distri[lauf3].CFY=PLTdist->Fall[0].Fluegel[lauf].Distri[lauf3].CFY;
      PLTdist->Fall[PLTpos].Fluegel[lauf].Distri[lauf3].CFZ=PLTdist->Fall[0].Fluegel[lauf].Distri[lauf3].CFZ;
      PLTdist->Fall[PLTpos].Fluegel[lauf].Distri[lauf3].CFNORMAL=PLTdist->Fall[0].Fluegel[lauf].Distri[lauf3].CFNORMAL;
      PLTdist->Fall[PLTpos].Fluegel[lauf].Distri[lauf3].CMX=PLTdist->Fall[0].Fluegel[lauf].Distri[lauf3].CMX;
      PLTdist->Fall[PLTpos].Fluegel[lauf].Distri[lauf3].CMX_FROM_CFY=PLTdist->Fall[0].Fluegel[lauf].Distri[lauf3].CMX_FROM_CFY;
      PLTdist->Fall[PLTpos].Fluegel[lauf].Distri[lauf3].CMX_FROM_CFZ=PLTdist->Fall[0].Fluegel[lauf].Distri[lauf3].CMX_FROM_CFZ;
      PLTdist->Fall[PLTpos].Fluegel[lauf].Distri[lauf3].CMY=PLTdist->Fall[0].Fluegel[lauf].Distri[lauf3].CMY;
      PLTdist->Fall[PLTpos].Fluegel[lauf].Distri[lauf3].CMY_FROM_CDI=PLTdist->Fall[0].Fluegel[lauf].Distri[lauf3].CMY_FROM_CDI;
      PLTdist->Fall[PLTpos].Fluegel[lauf].Distri[lauf3].CMY_FROM_CFZ=PLTdist->Fall[0].Fluegel[lauf].Distri[lauf3].CMY_FROM_CFZ;
      PLTdist->Fall[PLTpos].Fluegel[lauf].Distri[lauf3].CMZ=PLTdist->Fall[0].Fluegel[lauf].Distri[lauf3].CMZ;
      PLTdist->Fall[PLTpos].Fluegel[lauf].Distri[lauf3].CMZ_FROM_CDI=PLTdist->Fall[0].Fluegel[lauf].Distri[lauf3].CMZ_FROM_CDI;
      PLTdist->Fall[PLTpos].Fluegel[lauf].Distri[lauf3].CMZ_FROM_CFY=PLTdist->Fall[0].Fluegel[lauf].Distri[lauf3].CMZ_FROM_CFY;
      PLTdist->Fall[PLTpos].Fluegel[lauf].Distri[lauf3].REF_AREA=PLTdist->Fall[0].Fluegel[lauf].Distri[lauf3].REF_AREA;
      PLTdist->Fall[PLTpos].Fluegel[lauf].Distri[lauf3].REF_SPAN=PLTdist->Fall[0].Fluegel[lauf].Distri[lauf3].REF_SPAN;
      PLTdist->Fall[PLTpos].Fluegel[lauf].Distri[lauf3].REF_LEN_CMX=PLTdist->Fall[0].Fluegel[lauf].Distri[lauf3].REF_LEN_CMX;
      PLTdist->Fall[PLTpos].Fluegel[lauf].Distri[lauf3].REF_LEN_CMY=PLTdist->Fall[0].Fluegel[lauf].Distri[lauf3].REF_LEN_CMY;
      PLTdist->Fall[PLTpos].Fluegel[lauf].Distri[lauf3].REF_LEN_CMZ=PLTdist->Fall[0].Fluegel[lauf].Distri[lauf3].REF_LEN_CMZ;
      
      for (lauf4=0; lauf4<lauf2; lauf4++)
      {
        if (dummy1> sqrt(pow((PLTdist->Fall[PLTpos].Fluegel[lauf].Distri[lauf3].Y-BSURFout[lauf4].Ypos),2)+pow((PLTdist->Fall[PLTpos].Fluegel[lauf].Distri[lauf3].Z-BSURFout[lauf4].Zpos),2)))
        {
          dummy1=sqrt(pow((PLTdist->Fall[PLTpos].Fluegel[lauf].Distri[lauf3].Y-BSURFout[lauf4].Ypos),2)+pow((PLTdist->Fall[PLTpos].Fluegel[lauf].Distri[lauf3].Z-BSURFout[lauf4].Zpos),2));
          check=lauf4;
        }
      }
      VorZ=0;
      VorY=0;
      if (check == 0)
      {
        if ((PLTdist->Fall[PLTpos].Fluegel[lauf].Distri[lauf3].Y-BSURFout[0].Ypos) > 0 && (PLTdist->Fall[PLTpos].Fluegel[lauf].Distri[lauf3].Y-BSURFout[1].Ypos) < 0)
	{
	  VorY=1;
	}
        if ((PLTdist->Fall[PLTpos].Fluegel[lauf].Distri[lauf3].Y-BSURFout[0].Ypos) < 0 && (PLTdist->Fall[PLTpos].Fluegel[lauf].Distri[lauf3].Y-BSURFout[1].Ypos) > 0)
	{
	  VorY=1;
	}
        if ((PLTdist->Fall[PLTpos].Fluegel[lauf].Distri[lauf3].Z-BSURFout[0].Zpos) > 0 && (PLTdist->Fall[PLTpos].Fluegel[lauf].Distri[lauf3].Z-BSURFout[1].Zpos) < 0)
	{
	  VorZ=1;
	}
        if ((PLTdist->Fall[PLTpos].Fluegel[lauf].Distri[lauf3].Y-BSURFout[0].Zpos) < 0 && (PLTdist->Fall[PLTpos].Fluegel[lauf].Distri[lauf3].Y-BSURFout[1].Zpos) > 0)
	{
	  VorZ=1;
	}
        if ((VorZ==0) && (VorY==0))
	{
	  PLTdist->Fall[PLTpos].Fluegel[lauf].Distri[lauf3].CFX_FROM_CDI=BSURFout[0].CD;
	  PLTdist->Fall[PLTpos].Fluegel[lauf].Distri[lauf3].CFZ=BSURFout[0].CL;
	  PLTdist->Fall[PLTpos].Fluegel[lauf].Distri[lauf3].CFY=BSURFout[0].CS;
          PLTdist->Fall[PLTpos].Fluegel[lauf].Distri[lauf3].CFNORMAL=BSURFout[0].CN;
	} else {
	  get_dist(&dist1, &dist2, lauf, lauf3 ,0, 1, PLTdist ,BSURFout);
	  PLTdist->Fall[PLTpos].Fluegel[lauf].Distri[lauf3].CFX_FROM_CDI=(BSURFout[0].CD*dist2+BSURFout[1].CD*dist1)/(dist1+dist2);
	  PLTdist->Fall[PLTpos].Fluegel[lauf].Distri[lauf3].CFZ=(BSURFout[0].CL*dist2+BSURFout[1].CL*dist1)/(dist1+dist2);
	  PLTdist->Fall[PLTpos].Fluegel[lauf].Distri[lauf3].CFY=(BSURFout[0].CS*dist2+BSURFout[1].CS*dist1)/(dist1+dist2);
          PLTdist->Fall[PLTpos].Fluegel[lauf].Distri[lauf3].CFNORMAL=(BSURFout[0].CN*dist2+BSURFout[1].CN*dist1)/(dist1+dist2);
	}
      } else if (check == (lauf2-1)) {
        if ((PLTdist->Fall[PLTpos].Fluegel[lauf].Distri[lauf3].Y-BSURFout[check].Ypos) > 0 && (PLTdist->Fall[PLTpos].Fluegel[lauf].Distri[lauf3].Y-BSURFout[check-1].Ypos) < 0)
	{
	  VorY=1;
	}
        if ((PLTdist->Fall[PLTpos].Fluegel[lauf].Distri[lauf3].Y-BSURFout[check].Ypos) < 0 && (PLTdist->Fall[PLTpos].Fluegel[lauf].Distri[lauf3].Y-BSURFout[check-1].Ypos) > 0)
	{
	  VorY=1;
	}
        if ((PLTdist->Fall[PLTpos].Fluegel[lauf].Distri[lauf3].Z-BSURFout[check].Zpos) > 0 && (PLTdist->Fall[PLTpos].Fluegel[lauf].Distri[lauf3].Z-BSURFout[check-1].Zpos) < 0)
	{
	  VorZ=1;
	}
        if ((PLTdist->Fall[PLTpos].Fluegel[lauf].Distri[lauf3].Y-BSURFout[check].Zpos) < 0 && (PLTdist->Fall[PLTpos].Fluegel[lauf].Distri[lauf3].Y-BSURFout[check-1].Zpos) > 0)
	{
	  VorZ=1;
	}
        if ((VorZ==0) && (VorY==0))
	{
	  PLTdist->Fall[PLTpos].Fluegel[lauf].Distri[lauf3].CFX_FROM_CDI=BSURFout[check].CD;
	  PLTdist->Fall[PLTpos].Fluegel[lauf].Distri[lauf3].CFZ=BSURFout[check].CL;
	  PLTdist->Fall[PLTpos].Fluegel[lauf].Distri[lauf3].CFY=BSURFout[check].CS;
          PLTdist->Fall[PLTpos].Fluegel[lauf].Distri[lauf3].CFNORMAL=BSURFout[check].CN;
	} else {
	  get_dist(&dist1, &dist2, lauf, lauf3 ,check, check-1, PLTdist ,BSURFout);
	  PLTdist->Fall[PLTpos].Fluegel[lauf].Distri[lauf3].CFX_FROM_CDI=(BSURFout[check].CD*dist2+BSURFout[check-1].CD*dist1)/(dist1+dist2);
	  PLTdist->Fall[PLTpos].Fluegel[lauf].Distri[lauf3].CFZ=(BSURFout[check].CL*dist2+BSURFout[check-1].CL*dist1)/(dist1+dist2);
	  PLTdist->Fall[PLTpos].Fluegel[lauf].Distri[lauf3].CFY=(BSURFout[check].CS*dist2+BSURFout[check-1].CS*dist1)/(dist1+dist2);
          PLTdist->Fall[PLTpos].Fluegel[lauf].Distri[lauf3].CFNORMAL=(BSURFout[check].CN*dist2+BSURFout[check-1].CN*dist1)/(dist1+dist2);
	}      
      } else {
        VorY=0;
	VorZ=0;
	if ((PLTdist->Fall[PLTpos].Fluegel[lauf].Distri[lauf3].Y-BSURFout[check].Ypos) > 0 && (PLTdist->Fall[PLTpos].Fluegel[lauf].Distri[lauf3].Y-BSURFout[check-1].Ypos) < 0)
	{
	  VorY--;
	}
        if ((PLTdist->Fall[PLTpos].Fluegel[lauf].Distri[lauf3].Y-BSURFout[check].Ypos) < 0 && (PLTdist->Fall[PLTpos].Fluegel[lauf].Distri[lauf3].Y-BSURFout[check-1].Ypos) > 0)
	{
	  VorY--;
	}
        if ((PLTdist->Fall[PLTpos].Fluegel[lauf].Distri[lauf3].Z-BSURFout[check].Zpos) > 0 && (PLTdist->Fall[PLTpos].Fluegel[lauf].Distri[lauf3].Z-BSURFout[check-1].Zpos) < 0)
	{
	  VorZ--; 
	}
        if ((PLTdist->Fall[PLTpos].Fluegel[lauf].Distri[lauf3].Z-BSURFout[check].Zpos) < 0 && (PLTdist->Fall[PLTpos].Fluegel[lauf].Distri[lauf3].Z-BSURFout[check-1].Zpos) > 0)
	{
	  VorZ--;
	}

	if ((PLTdist->Fall[PLTpos].Fluegel[lauf].Distri[lauf3].Y-BSURFout[check].Ypos) > 0 && (PLTdist->Fall[PLTpos].Fluegel[lauf].Distri[lauf3].Y-BSURFout[check+1].Ypos) < 0)
	{
	  VorY++;
	}
        if ((PLTdist->Fall[PLTpos].Fluegel[lauf].Distri[lauf3].Y-BSURFout[check].Ypos) < 0 && (PLTdist->Fall[PLTpos].Fluegel[lauf].Distri[lauf3].Y-BSURFout[check+1].Ypos) > 0)
	{
	  VorY++;
	}
        if ((PLTdist->Fall[PLTpos].Fluegel[lauf].Distri[lauf3].Z-BSURFout[check].Zpos) > 0 && (PLTdist->Fall[PLTpos].Fluegel[lauf].Distri[lauf3].Z-BSURFout[check+1].Zpos) < 0)
	{
	  VorZ++;
	}
        if ((PLTdist->Fall[PLTpos].Fluegel[lauf].Distri[lauf3].Z-BSURFout[check].Zpos) < 0 && (PLTdist->Fall[PLTpos].Fluegel[lauf].Distri[lauf3].Z-BSURFout[check+1].Zpos) > 0)
	{
	  VorZ++;
	}
        /*if ((VorZ == 0) && (VorY == 0))
	{
	  printf ("check: %d; \n",check);
          printf ("BsurfYpos: %lf; \n",BSURFout[check].Ypos);
          printf ("PLTYpos: %lf; \n",PLTdist->Fall[PLTpos].Fluegel[lauf].Distri[lauf3].Y);
          errmessage(40);
	}*/
	//if ((VorY >= 0) && (VorZ >= 0)) Ok hier erwarte ich wohl noch eine rein spannweitige Änderung
        if ((VorY > 0) ) 
	{
	  get_dist(&dist1, &dist2, lauf, lauf3 ,check, check+1, PLTdist ,BSURFout);
	  PLTdist->Fall[PLTpos].Fluegel[lauf].Distri[lauf3].CFX_FROM_CDI=(BSURFout[check].CD*dist2+BSURFout[check+1].CD*dist1)/(dist1+dist2);
	  PLTdist->Fall[PLTpos].Fluegel[lauf].Distri[lauf3].CFZ=(BSURFout[check].CL*dist2+BSURFout[check+1].CL*dist1)/(dist1+dist2);
	  PLTdist->Fall[PLTpos].Fluegel[lauf].Distri[lauf3].CFY=(BSURFout[check].CS*dist2+BSURFout[check+1].CS*dist1)/(dist1+dist2);
          PLTdist->Fall[PLTpos].Fluegel[lauf].Distri[lauf3].CFNORMAL=(BSURFout[check].CN*dist2+BSURFout[check+1].CN*dist1)/(dist1+dist2);
	}     
	//if ((VorY <= 0) && (VorZ <= 0)) Ok hier erwarte ich wohl noch eine rein spannweitige Änderung
        if ((VorY < 0))
	{
	  get_dist(&dist1, &dist2, lauf, lauf3 ,check, check-1, PLTdist ,BSURFout);
	  PLTdist->Fall[PLTpos].Fluegel[lauf].Distri[lauf3].CFX_FROM_CDI=(BSURFout[check].CD*dist2+BSURFout[check-1].CD*dist1)/(dist1+dist2);
	  PLTdist->Fall[PLTpos].Fluegel[lauf].Distri[lauf3].CFZ=(BSURFout[check].CL*dist2+BSURFout[check-1].CL*dist1)/(dist1+dist2);
	  PLTdist->Fall[PLTpos].Fluegel[lauf].Distri[lauf3].CFY=(BSURFout[check].CS*dist2+BSURFout[check-1].CS*dist1)/(dist1+dist2);
          PLTdist->Fall[PLTpos].Fluegel[lauf].Distri[lauf3].CFNORMAL=(BSURFout[check].CN*dist2+BSURFout[check-1].CN*dist1)/(dist1+dist2);
	}
        if ((VorY == 0))
	{
	  get_dist(&dist1, &dist2, lauf, lauf3 ,check, check-1, PLTdist ,BSURFout);
	  PLTdist->Fall[PLTpos].Fluegel[lauf].Distri[lauf3].CFX_FROM_CDI=BSURFout[check].CD;
	  PLTdist->Fall[PLTpos].Fluegel[lauf].Distri[lauf3].CFZ=BSURFout[check].CL;
	  PLTdist->Fall[PLTpos].Fluegel[lauf].Distri[lauf3].CFY=BSURFout[check].CS;
          PLTdist->Fall[PLTpos].Fluegel[lauf].Distri[lauf3].CFNORMAL=BSURFout[check].CN;
	}
   
      }  
    }
    free(BSURFout);
  }
}
void Comp_Square_Diff(struct inputfile *input, struct PLTdistribution *PLTdist)
{
  int lauf1 , lauf2;
  double difference,AEREA,chord,CN;
  char ExecName[150] = {"\0"}; 
  FILE *fSquare, *fplt;
  if (PLTdist->AnzFaelle > 1)
  {
    errmessage (58);
  }
  
  fplt=fopen("LoadCompare.plt","w");
  fprintf(fplt,"TITLE = \"Coefficient-Distribution for Load-Comparrison\" \n");
  fprintf(fplt,"  VARIABLES = X, Y, Z, CN, CN*l, CN*S \n");
  fprintf(fplt,"ZONE T=\"Reference Distribution\" \n");
  for (lauf1=0; lauf1 < input->NbWings; lauf1++)
  {
    for (lauf2=0; lauf2 < input->Fluegel[lauf1].AnzahlSpanPanel; lauf2++)
    {
      AEREA=PLTdist->Fall[0].Fluegel[lauf1].Distri[lauf2].REF_AREA;
      chord=PLTdist->Fall[0].Fluegel[lauf1].Distri[lauf2].REF_LEN_CMX;
      CN=PLTdist->Fall[0].Fluegel[lauf1].Distri[lauf2].CFNORMAL;
      fprintf(fplt,"%e %e %e %e %e %e %e\n",PLTdist->Fall[0].Fluegel[lauf1].Distri[lauf2].X,PLTdist->Fall[0].Fluegel[lauf1].Distri[lauf2].Y,PLTdist->Fall[0].Fluegel[lauf1].Distri[lauf2].Z,PLTdist->Fall[0].Fluegel[lauf1].Distri[lauf2].S,CN,CN*chord,CN*AEREA);
    }
  }
  fprintf(fplt,"ZONE T=\"Lifting-Line Distribution\" \n");
  for (lauf1=0; lauf1 < input->NbWings; lauf1++)
  {
    for (lauf2=0; lauf2 < input->Fluegel[lauf1].AnzahlSpanPanel; lauf2++)
    {
      AEREA=PLTdist->Fall[0].Fluegel[lauf1].Distri[lauf2].REF_AREA;
      chord=PLTdist->Fall[0].Fluegel[lauf1].Distri[lauf2].REF_LEN_CMX;
      CN=PLTdist->Fall[1].Fluegel[lauf1].Distri[lauf2].CFNORMAL;
      fprintf(fplt,"%e %e %e %e %e %e %e\n",PLTdist->Fall[1].Fluegel[lauf1].Distri[lauf2].X,PLTdist->Fall[1].Fluegel[lauf1].Distri[lauf2].Y,PLTdist->Fall[1].Fluegel[lauf1].Distri[lauf2].Z,PLTdist->Fall[1].Fluegel[lauf1].Distri[lauf2].S,CN,CN*chord,CN*AEREA);
    }
  }
  fclose(fplt);
  strcpy(ExecName,"mv LoadCompare.plt LaufV1.lili.");
  strcat(ExecName, input->LiLiVersionName);
  system(ExecName);
  
  
  difference=0.0;
  for (lauf1=0; lauf1 < input->NbWings; lauf1++)
  {
    for (lauf2=0; lauf2 < input->Fluegel[lauf1].AnzahlSpanPanel; lauf2++)
    {
      difference=difference+pow(PLTdist->Fall[0].Fluegel[lauf1].Distri[lauf2].CFNORMAL-PLTdist->Fall[1].Fluegel[lauf1].Distri[lauf2].CFNORMAL,2);
    }
  }
  printf("  square-diffrence_CN: %e\n",difference);
  fSquare=fopen("SquareDiff_TAU-LiLi_CN.dat","w");
  fprintf(fSquare,"%e",difference);
  fclose(fSquare);  

  difference=0.0;
  for (lauf1=0; lauf1 < input->NbWings; lauf1++)
  {
    for (lauf2=0; lauf2 < input->Fluegel[lauf1].AnzahlSpanPanel; lauf2++)
    {
      AEREA=PLTdist->Fall[0].Fluegel[lauf1].Distri[lauf2].REF_AREA;
      difference=difference+pow(PLTdist->Fall[0].Fluegel[lauf1].Distri[lauf2].CFNORMAL*AEREA-PLTdist->Fall[1].Fluegel[lauf1].Distri[lauf2].CFNORMAL*AEREA,2);
    }
  }
  printf("  square-diffrence_CNArea: %e\n",difference);
  fSquare=fopen("SquareDiff_TAU-LiLi_CNArea.dat","w");
  fprintf(fSquare,"%e",difference);
  fclose(fSquare);  

  difference=0.0;
  for (lauf1=0; lauf1 < input->NbWings; lauf1++)
  {
    for (lauf2=0; lauf2 < input->Fluegel[lauf1].AnzahlSpanPanel; lauf2++)
    {
      chord=PLTdist->Fall[0].Fluegel[lauf1].Distri[lauf2].REF_LEN_CMX;
      difference=difference+pow(PLTdist->Fall[0].Fluegel[lauf1].Distri[lauf2].CFNORMAL*chord-PLTdist->Fall[1].Fluegel[lauf1].Distri[lauf2].CFNORMAL*chord,2);
    }
  }
  printf("  square-diffrence_CNchord: %e\n",difference);
  fSquare=fopen("SquareDiff_TAU-LiLi_CNchord.dat","w");
  fprintf(fSquare,"%e",difference);
  fclose(fSquare);  

}
