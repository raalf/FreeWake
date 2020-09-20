/******************POALR_INTER.h******************                  
*                                                *      
*      Vorentwurfs-Generator "VoGen"             *          
*      Version 1.15                              *          
*      Date of last modification 17.11.2017      *          
*                                                *          
*      Writen by: Jan Himisch                    *      
*             jan.himisch@dlr.de                 *      
*                                                *      
*                                                *
*************************************************/
#ifndef _POLARINTER_H_
#define _POLARINTER_H_

#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "run_LILI.h"
#include "input.h"
#include "errmes.h"
#include "LinGLloeser.h"

struct EinWert
{
  double Alfa;
  double Ca;
  double Cd;
  double CDP;
  double CDPW;   
  double CDPP;   
  double CDPF;   
  double Cm;
};

struct EinePolare
{
  double Eta_NUM;
  int anzahlWerte;
  struct EinWert *Werte;
  struct EinePolare *next;
};

struct PolarenRe
{
  double RE_NUM;
  int anzahlEtaNum;
  struct EinePolare *EtaPolare;
  struct PolarenRe *next;
};

struct PolarenMa
{
  double MA_NUM;
  int anzahlReNum;
  struct PolarenRe *RePolaren;
  struct PolarenMa *next;
};

struct EinProfiel
{
  int AnzMaNummern;
  char Profiel_Name[50];
  struct PolarenMa *MaPolare;  
  struct EinProfiel *next;
};

struct Profiele
{
  int Anzahl_InputProfile;
  struct EinProfiel *Profiel;
};

struct VERTEILUNG_POLARint
{
  double X;
  double Y;
  double Z;
  double S;
  double Alfa;
  double CA;
  double CD;
  double CDI;
  double CDPW;
  double CDPF;
  double CDPP;
  double CDP;
  double CM;
  double tiefe;
  double RE_LOC;
  double Vloc;
  double areaPan;
  char SEC_Name[50];
  double SEC_Percent;
  char SEC2_Name[50];
  double Flap_Set;
};

struct FLUEGELP
{
  struct VERTEILUNG_POLARint *Distri;
};

struct EinFallP
{
  struct FLUEGELP *Fluegel;
};

struct Polar_Integration_Distribution
{
  int AnzFluegel;
  int AnzFaelle;
  struct EinFallP *Fall;
};

struct GESAMT_POLARint
{
  double Alfa;
  double CA;
  double CD;
  double CDI;
  double CDPW;
  double CDPF;
  double CDPP;
  double CDP;
  double CM;
  double AERA;
  double Vinf;
  double RE_Flug;
  double Ma_Flug;
};

void PolarInterpolation(struct inputfile *input, struct PLTdistribution *PLTdist, struct file14 *rfile, struct GESAMT_POLARint *GPolare);
#endif
