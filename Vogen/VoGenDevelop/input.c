/********************input.c**********************                  
*                                                *      
*      Vorentwurfs-Generator "VoGen"             *          
*      Version 2.00                              *
*      Date of last modification 21.04.2019      *
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
#include <ctype.h>

#include "input.h"
#include "errmes.h"
#include "trimming.h"
#include "mischer.h"

void get_Vloc(struct inputfile *input)
{
  int lauf, lauf2;
  double V1, V2;
  for (lauf=0; lauf < input->NbWings; lauf++)
  {
    for (lauf2=0; lauf2 < input->Fluegel[lauf].anzahlSchnitte; lauf2++)
    {
      if ((input->Fluegel[lauf].SymExt==1) && (lauf2==0))
      {
        V1=180/PI*atan((input->Fluegel[lauf].Schnitte[lauf2+1].posz-input->Fluegel[lauf].Schnitte[lauf2].posz)/(input->Fluegel[lauf].Schnitte[lauf2+1].posy-input->Fluegel[lauf].Schnitte[lauf2].posy));
        input->Fluegel[lauf].Schnitte[lauf2].Vloc=V1/2;
      }else if ((input->Fluegel[lauf].SymExt==0) && (lauf2==0))
      {
        V1=180/PI*atan((input->Fluegel[lauf].Schnitte[lauf2+1].posz-input->Fluegel[lauf].Schnitte[lauf2].posz)/(input->Fluegel[lauf].Schnitte[lauf2+1].posy-input->Fluegel[lauf].Schnitte[lauf2].posy));
        input->Fluegel[lauf].Schnitte[lauf2].Vloc=V1;
      }else if (lauf2==(input->Fluegel[lauf].anzahlSchnitte-1))
      {
        V2=180/PI*atan((input->Fluegel[lauf].Schnitte[lauf2].posz-input->Fluegel[lauf].Schnitte[lauf2-1].posz)/(input->Fluegel[lauf].Schnitte[lauf2].posy-input->Fluegel[lauf].Schnitte[lauf2-1].posy));
        input->Fluegel[lauf].Schnitte[lauf2].Vloc=V2;
      } else {
        V1=180/PI*atan((input->Fluegel[lauf].Schnitte[lauf2+1].posz-input->Fluegel[lauf].Schnitte[lauf2].posz)/(input->Fluegel[lauf].Schnitte[lauf2+1].posy-input->Fluegel[lauf].Schnitte[lauf2].posy));
        V2=180/PI*atan((input->Fluegel[lauf].Schnitte[lauf2].posz-input->Fluegel[lauf].Schnitte[lauf2-1].posz)/(input->Fluegel[lauf].Schnitte[lauf2].posy-input->Fluegel[lauf].Schnitte[lauf2-1].posy));      
        input->Fluegel[lauf].Schnitte[lauf2].Vloc=(V1+V2)/2;      
      }
    }
  }  
}

int roundfloat(float x)
{
  if(x>0) return (int)(x + 0.5);

  return (int)(x - 0.5);
}

void init_Einfluegel(struct Einfluegel *Fluegel)
{
  Fluegel->SymExt=1;
  Fluegel->SymRandBed=-5;
  Fluegel->Spiegeln=0;
}

struct EinSchnitt copy_Schnitt(struct EinSchnitt *ESchnitt)
{
   struct EinSchnitt NSchnitt;
   int lauf;
   
   NSchnitt.posx=ESchnitt->posx;
   NSchnitt.posy=ESchnitt->posy;
   NSchnitt.posz=ESchnitt->posz; 
   NSchnitt.OrgTwist=ESchnitt->OrgTwist;
   NSchnitt.AlphaPlus=ESchnitt->AlphaPlus;
   NSchnitt.AnzahlPan=ESchnitt->AnzahlPan;
   NSchnitt.tiefe=ESchnitt->tiefe;
   //if (readRelDick == 1)
   //{
     NSchnitt.relDicke=ESchnitt->relDicke;
   //}
   NSchnitt.AnzahlXY=ESchnitt->AnzahlXY;
   for (lauf=0;lauf<=NSchnitt.AnzahlXY;lauf++)
   {
     NSchnitt.xcamb[lauf]=ESchnitt->xcamb[lauf];
     NSchnitt.zcamb[lauf]=ESchnitt->zcamb[lauf]; 
   }
   sscanf(ESchnitt->ProfilName,"%s",NSchnitt.ProfilName);
   NSchnitt.KlappeNrI=0;
   NSchnitt.KlappeNrA=0;
   NSchnitt.etaSchnitt=0.0;
   NSchnitt.anzPanelKlappe=0;
   return NSchnitt;     
}

struct EinSchnitt read_Schnitt(char Datei[], double scalfac, int readRelDick, int XMLKOP)
{
  int zahl,zahl2, check;
  char ZEILE_READ[MaxZeilenLaenge];
  struct EinSchnitt Schnitt;
  
  FILE *fopen(),*ffile; 
  strcpy(Schnitt.ProfilName,"\0");
  ffile = fopen(Datei,"r");
  if (ffile == NULL)
  {
    printf("\n Could not open file %s\n" , Datei);
    errmessage (12);
  }
  fgets (ZEILE_READ,MaxZeilenLaenge,ffile);
  fgets (ZEILE_READ,MaxZeilenLaenge,ffile);
  sscanf(ZEILE_READ," %lf %lf %lf %lf \n", &Schnitt.posx, &Schnitt.posy, &Schnitt.posz, &Schnitt.AlphaPlus);
  Schnitt.posx=Schnitt.posx*scalfac;
  Schnitt.posy=Schnitt.posy*scalfac;
  Schnitt.posz=Schnitt.posz*scalfac;
  Schnitt.OrgTwist=Schnitt.AlphaPlus;
  fgets (ZEILE_READ,MaxZeilenLaenge,ffile);
  fgets (ZEILE_READ,MaxZeilenLaenge,ffile);
  if (readRelDick == 1)
  {
    check=sscanf(ZEILE_READ," %d %lf %lf\n",  &Schnitt.AnzahlPan, &Schnitt.tiefe, &Schnitt.relDicke);
    if (check<3)
    {
      errmessage(35);
    }
  }else{
    check=sscanf(ZEILE_READ," %d %lf \n", &Schnitt.AnzahlPan, &Schnitt.tiefe);
  }
  //printf("Tiefe:%lf\n",Schnitt.tiefe);
  Schnitt.tiefe=Schnitt.tiefe*scalfac;
  fgets (ZEILE_READ,MaxZeilenLaenge,ffile);
  if (strstr(ZEILE_READ, "|--name-section--|")!=NULL)
  {
    //Schnitt.ProfilName=(char *)malloc(MAX_LENGTH_PNAME*sizeof(char));
    fgets (ZEILE_READ,MaxZeilenLaenge,ffile);
    sscanf(ZEILE_READ," %s \n",Schnitt.ProfilName);
    fgets (ZEILE_READ,MaxZeilenLaenge,ffile);
    //printf(" %s \n",Schnitt.ProfilName);
  }else{
    if (XMLKOP==1)
    {
      errmessage(22);
    }
  }
  zahl=0;
  zahl2=0;
  if (strstr(ZEILE_READ, "|--section-file--|")!=NULL)
  {
    char Profildatei[256];   
    fgets (ZEILE_READ,MaxZeilenLaenge,ffile);
    fclose(ffile);
    sscanf(ZEILE_READ," %s \n",Profildatei);
    ffile = fopen(Profildatei,"r");
    if (ffile == NULL)
    {
      printf("\n Could not open file %s\n" , Profildatei);
      errmessage (12);
    }
  }
  fgets (ZEILE_READ,MaxZeilenLaenge,ffile);
  do 
  {
    check=sscanf(ZEILE_READ," %lf %lf\n",&Schnitt.xcamb[zahl], &Schnitt.zcamb[zahl]);
    //printf("zCAMB: %lf   CHECK:%d\n",Schnitt.zcamb[zahl], check);
    if (check<2)
    {
      break;
    }
    zahl++;    
    fgets (ZEILE_READ,MaxZeilenLaenge,ffile);     
  } 
  while (!feof(ffile)); 
  //zahl--;
  Schnitt.AnzahlXY=zahl-1;
  //printf("SchnittName: %s     AnzahlXY:%d\n\n",Datei ,zahl-1);
  if  (zahl-1 < 1)
  {
    errmessage(21);
  }
  if (Schnitt.xcamb[0]!=0)
  {
    for (zahl2=1; zahl2<zahl ; zahl2++)
    {
      Schnitt.xcamb[zahl2] = Schnitt.xcamb[zahl2] - Schnitt.xcamb[0];
    }
    Schnitt.xcamb[0]=0;
  }
  if (Schnitt.xcamb[zahl-1]!=1)
  {
    for (zahl2=0; zahl2<zahl ; zahl2++)
    {
      Schnitt.xcamb[zahl2] = Schnitt.xcamb[zahl2] / Schnitt.xcamb[zahl-1];
    }
  }
  fclose(ffile);
  /*printf("Posx %lf\n",Schnitt.posx); */
  // setzen der Klappen inits
  Schnitt.KlappeNrI=0;
  Schnitt.KlappeNrA=0;
  Schnitt.etaSchnitt=0.0;
  Schnitt.anzPanelKlappe=0;
  return Schnitt;  
}

void init_input(struct inputfile *input)
{
  
  input->debuglevel=0;
  
  input->MachNumber = 0; 
  input->beta = 0;
  input->refSpan = 0;
  input->refAerea = 0;                       
  input->origXYZ[0] = 0.0;
  input->origXYZ[1] = 0.0; 
  input->origXYZ[2] = 0.0;                             
                                                                       
  input->Alfa_Rot = 0;                          
                                                            
  input->numberCL = 0;                         
  input->targetCl = NULL;                        
  input->numberAlpha = 0;                          
  input->targetAlpha = NULL;                       
  input->LiliStart = 1;
  //ZEIGER =  input.LILIEXE[0];                        
  strcpy( input->LiLiVersion, "V2p3b");
  strcpy( input->LILIEXE, "./LIFTING_LINE_LINUX_32BIT.EXE");          
  strcpy( input->LAUFNAME, "LaufV1.inp");
  input->runOneByOne = 0;
  input->XMLEA = 0;                               
  input->OldInput = 0;

  strcpy( input->FwVersion, "V0p1");
  strcpy( input->FwEXE, "./fw");
  input->FwStart = 0;
  input->FwSteadyUnstedy = 0;
  input->FwRelaxWake = 0;
  input->FwnumbTimeStep = 10;
  input->FwtimeStepLength = 0.0025;
  input->FwConvergDelta = 0.00;
                                             
  input->GlobSym = 1;
  input->scalfac = 1; 
  input->NbWings = 1;                          
  input->CamberMethod = 1;
  input->CamberStart = 0.5;
  input->CamberEnd = 1.0;    
  input->AnzKreuz = 0;                             
  input->readRelDick = 0;                            
  input->AlfaPlusRumpf = 0;
  input->Fluegel = NULL;
  input->Kreuz = NULL;
  input->AnzahlKlappen = 0;
  input->Klappen = NULL;


  input->AnzahlMischer = 0;

  input->NumberOfParametSetsinput = 0;
  //input->parametric = 0;                         
  //input->NrVariations = 1;
  //input->anzVar = 0;
  //input->anzSpVar = 0;                       
  //input->VarParameter = NULL;
                                         
  input->Austrimmen = 0;                         
  input->MaxCmIter = 3;                                                 
  input->ZielCm = 0;
  input->CmSteuerFl = 0;
  input->CmSteuerKlappe = 0;
  input->ArotAdjust = 0;
  input->InitStepSize = 1.0;
  strcpy(input->EpsAusgabe, "CmTrimAngles.dat");

  //input->zusAufP=0; 
  input->IBM_BD_WING = 1;
  input->IBM_BD_SEC = 1;
  input->WRBM_BD_WING = 1;
  input->WRBM_BD_SEC = 1;
  input->RefSpeed = -1;
  input->RefDensity = -1;
  input->RefTemp = -1;
  input->BasisIntBiegeMoment = 0;
  input->BasisIntDickenBiegeMoment = 0;
  
  input->makeBsurfinp = 0;
  input->TAUbsurfINP = 0;
  input->NUM_OF_ZONES = 0;
  input->BSURF_FLUEGEL=NULL;
  input->BsurfSquareDiff = 0; 	   
  
  input->kHOrd = 0;         
  input->kH = NULL;
  input->BCwi = 0;             
  input->BWRBM = 0;            
  input->BCw0 = 0;
  input->UrArea = 0; 
  input->BAlafa = 0;
  
  input->Struckt_Coup_ON = 0;
  input->Load_Case_ON = 0;
  input->Numb_Iter = 0;
  input->Log_Iter_ON = 0;
  input->useJIG2FLIGHT = 0;
  input->Target_CL_LoadCase = 0;
  input->Ref_Density_LoadCase = 0;
  input->Ref_Speed_LoadCase = 0;
  input->ProjektionsFluegel = 1;
  input->ProjektionsSchnitt = -1;
  input->deltaStructX = 0;
  input->deltaStructZ = 0;
  input->NumStructDih = 0;
  input->StructModBezSpan = 0;
  input->NumStructDih = 0;
  input->scalefactStruct = 1.0;
  strcpy( input->LILIDIR, "\0");
  strcpy( input->StructProg_ResultSubDir, "\0");
  strcpy( input->StructProg_EXEdir, "\0");
  
  input->KreisflugModus = 0;
  input->HaengeWinkel = 0;
  input->Steigen = 0;
  
  input->XMLKOP=0;
  input->MASSE_SPEED=0;
  input->PolintStart=0;
  input->InterneInterpolation=0;
  input->BASIS_CD_POL=0;
  input->POLINT_MASSE=0;
  input->useCdpw=0;
  input->useCdpp=0;
  input->useCdpf=0;
  input->useCm=1;
  input->ClMinMaxHandling=0;
  input->PENALTYfact=1.1;
}

void read_input(char datei[],struct inputfile *input)
{
  int zahl, zaehler2, lauf;
  int IstVar,NewVar;
  double aEbene;
  char ZEILE_READ[MaxZeilenLaenge]; 
  char ZEILE_CLEAN[MaxZeilenLaenge];
  char *ZEIGER, *ZEIGER2;
  char laufNR[3];
  FILE *fopen(),*fpara;
  
  
  IstVar=0;
  NewVar=0;
  init_input(input);
  fpara = fopen(datei,"r");  
  if (fpara == NULL)
  {
    printf("\n Could not read specified parameter file %s \n",datei);
    errmessage (12);
  }else {
    printf("reading parameter file\n");
  }
  
  fgets (ZEILE_READ,MaxZeilenLaenge,fpara);
  do
  {
    zahl =0;
    zaehler2=0;
    do  
    {    
      if (isspace(ZEILE_READ[zahl])==0)
      { 
        ZEILE_CLEAN[zaehler2]  = ZEILE_READ[zahl];
        zahl++;
        zaehler2 = zaehler2 + 1;
      }else
      {
        zahl++;
      }
    } 
    while (ZEILE_READ[zahl] != '\0');
    for (zahl=0;zahl < MaxZeilenLaenge;zahl++)
    {
      ZEILE_READ[zahl] = '\0';
    }
    
    zahl = 0;
    zaehler2 = 0; 
    do 
    {
      ZEILE_READ[zahl]  = ZEILE_CLEAN[zahl];
      zahl++;
    }
    while ((ZEILE_CLEAN[zahl] != ':') && (ZEILE_CLEAN[zahl] != '\0'));
    ZEIGER=&ZEILE_CLEAN[zahl+1];
    //printf("Aktuelle Zeile: %s\n", ZEILE_CLEAN);
    if (strcmp( ZEILE_READ, "debuglevel") == 0 ) 
    {
      sscanf(ZEIGER, "%d", &input->debuglevel);
    } 
    //HEADER Lifting-Line
    if (strcmp( ZEILE_READ, "machnumber") == 0 ) 
    {
      sscanf(ZEIGER, "%lf", &input->MachNumber);
    } 
    if (strcmp( ZEILE_READ, "sideslipangle") == 0 ) 
    {
      sscanf(ZEIGER, "%lf", &input->beta);
    } 
    if ( strcmp( ZEILE_READ, "referencespan[m]") == 0 ) 
    {
      sscanf(ZEIGER, "%lf", &input->refSpan);
    }
    if ( strcmp( ZEILE_READ, "referenceaerea[m2]") == 0 ) 
    {
      sscanf(ZEIGER, "%lf", &input->refAerea);
    }   
    if ( strcmp( ZEILE_READ, "referencechordlength[m]") == 0 ) 
    {
      sscanf(ZEIGER, "%lf", &input->refChordlength);
    }
    if ( strcmp( ZEILE_READ, "originX[m]") == 0 ) 
    {
      sscanf(ZEIGER, "%lf", &input->origXYZ[0]);
    }
    if ( strcmp( ZEILE_READ, "originY[m]") == 0 ) 
    {
      sscanf(ZEIGER, "%lf", &input->origXYZ[1]);
    }
    if ( strcmp( ZEILE_READ, "originZ[m]") == 0 ) 
    {
      sscanf(ZEIGER, "%lf", &input->origXYZ[2]);
    }
    if ( strcmp( ZEILE_READ, "targetalfarot") == 0)
    {
      sscanf(ZEIGER,"%lf", &input->Alfa_Rot);
    }

    //Lifting-Line-Steuerung  
    if ( strcmp( ZEILE_READ, "numberoftargetCL") == 0 ) 
    {
      if (input->targetCl == NULL)
      {
        sscanf(ZEIGER, "%d", &input->numberCL);
      } else {
        errmessage(13);
      }
    }
    if ( strncmp( ZEILE_READ , "targetCL", 8) == 0)
    {
      for (lauf=1; lauf <= input->numberCL ; lauf++)
      {   
        if (input->targetCl == NULL)
        {
          input->targetCl = (double *)malloc(input->numberCL*sizeof(double));
        } 
	if (input->numberCL != 0)
        {
          char DUMMY[11 ]= {"\0"};    
          strcpy(DUMMY,"targetCL");
	  if (lauf < 10)
          {
            laufNR[0]=48+lauf;
            laufNR[1]='\0';
	    laufNR[2]='\0';
          }else
          {  
            laufNR[0]=48+(int)(lauf/10);
            laufNR[1]=48+lauf-((int)(lauf/10)*10);
            laufNR[2]='\0';
          }      
          strcat(DUMMY,laufNR); 
          if ( strcmp( ZEILE_READ, DUMMY) == 0 ) 
          {
            sscanf(ZEIGER, "%lf", &input->targetCl[lauf-1]);
          }
        }
      }
    }
    if ( strcmp( ZEILE_READ, "numberoftargetalfa") == 0 ) 
    {
      if (input->targetAlpha == NULL)
    {
        sscanf(ZEIGER, "%d", &input->numberAlpha);
      } else {
        errmessage(14);
    }
    }
    if ( strncmp( ZEILE_READ , "targetalfa", 9) == 0)
    {
      for (lauf=1; lauf<=input->numberAlpha; lauf++)
      {   
        if (input->targetAlpha  == NULL)
        {
          input->targetAlpha = (double *)malloc(input->numberAlpha*sizeof(double));
        } 
	if (input->numberAlpha != 0)
        {
          char DUMMY[14]="targetalfa";    
          if (lauf < 10)
          {
            laufNR[0]=48+lauf;
            laufNR[1]='\0';
          }else
          {
            laufNR[0]=48+(int)(lauf/10);
            laufNR[1]=48+lauf-((int)(lauf/10)*10);
            laufNR[2]='\0';
          }      
          strcat(DUMMY,laufNR);     
          if ( strcmp( ZEILE_READ, DUMMY) == 0 ) 
          {
            sscanf(ZEIGER, "%lf", &input->targetAlpha[lauf-1]);;
          }
        }
      }
    }
    if ( strcmp( ZEILE_READ, "startLiftingLine(0=off/1=input/2=run/3=+PP)") == 0 )
    {
      sscanf(ZEIGER, "%d", &input->LiliStart);
    }    
    if ( strcmp( ZEILE_READ, "LIFTINGLINEVERSION") == 0 ) 
    {
      sscanf(ZEIGER, "%s", &input->LiLiVersion[0]);
    }        
    if ( strcmp( ZEILE_READ, "LIFTINGLINEEXECUTECOMAND") == 0 ) 
    {
      sscanf(ZEIGER, "%s", &input->LILIEXE[0]);
    }
    if ( strcmp( ZEILE_READ, "runonebyone(on=1/off=0)") == 0 ) 
    {
      sscanf(ZEIGER, "%d", &input->runOneByOne);
    }
    if ( strcmp( ZEILE_READ, "XMLoutput(on=1/off=0)") == 0 ) 
    {
      sscanf(ZEIGER, "%d", &input->XMLEA);
    }
    if ( strcmp( ZEILE_READ, "useoldinputfile(on=1/off=0)") == 0 ) 
    {
      sscanf(ZEIGER, "%d", &input->OldInput);
    }
    
    //FreeWake
    if ( strcmp( ZEILE_READ, "startFreeWake(0=off/1=input/2=run/3=+PP)") == 0 )
    {
      sscanf(ZEIGER, "%d", &input->FwStart);
    }
    if ( strcmp( ZEILE_READ, "FreeWakeVersion") == 0 )
    {
      sscanf(ZEIGER, "%s", &input->FwVersion[0]);
    }
    if ( strcmp( ZEILE_READ, "FreeWakeEXECUTECOMAND") == 0 )
    {
      sscanf(ZEIGER, "%s", &input->FwEXE[0]);
    }
    if ( strcmp( ZEILE_READ, "RelaxedWake") == 0 )
    {
      sscanf(ZEIGER, "%d", &input->FwRelaxWake);
    }
    if ( strcmp( ZEILE_READ, "Unsteady(on=1/off=0)") == 0 )
    {
      sscanf(ZEIGER, "%d", &input->FwSteadyUnstedy);
    }
    if ( strcmp( ZEILE_READ, "Max.numberoftimesteps") == 0 )
    {
      sscanf(ZEIGER, "%d", &input->FwnumbTimeStep);
    }
    if ( strcmp( ZEILE_READ, "Widthofeachtimestep(sec)") == 0 )
    {
      sscanf(ZEIGER, "%lf", &input->FwtimeStepLength);
    }
    if ( strcmp( ZEILE_READ, "Convergencedelta-spaneffic.") == 0 )
    {
      sscanf(ZEIGER, "%lf", &input->FwConvergDelta);
    }


    //Geometrische Fl�gelbedingungen    
    
    if ( strcmp( ZEILE_READ, "globalsymmetry(on=1/off=0)") == 0 ) 
    {
      sscanf(ZEIGER, "%d", &input->GlobSym);
    }      
    if ( strcmp( ZEILE_READ, "scalefactor") == 0 ) 
    {
      sscanf(ZEIGER, "%lf", &input->scalfac);
    }  
    if ( strcmp( ZEILE_READ, "camberanalysismethod(classic=1/integral=2)") == 0 ) 
    {
      sscanf(ZEIGER, "%d", &input->CamberMethod);    
    }  
    if ( strcmp( ZEILE_READ, "camberstartvalue(typical=0.5)") == 0 ) 
    {
      sscanf(ZEIGER, "%lf", &input->CamberStart);    
    }  
    if ( strcmp( ZEILE_READ, "camberstopvalue(typical=1.0)") == 0 ) 
    {
      sscanf(ZEIGER, "%lf", &input->CamberEnd);    
    }  
    if ( strcmp( ZEILE_READ, "numberofwings") == 0 ) 
    {
      if (input->Fluegel == NULL)
      {
        sscanf(ZEIGER, "%d", &input->NbWings);
        input->Fluegel = (struct Einfluegel *)malloc(input->NbWings*sizeof(struct Einfluegel));
      } else {
        errmessage(15);
      }
    }
    if ( strcmp( ZEILE_READ, "numberofjunctions") == 0 ) 
    {
      int i;
      if (input->Kreuz == NULL)
      {
    	sscanf(ZEIGER, "%d", &input->AnzKreuz);
    	input->Kreuz = (struct Kreuzung *)malloc(input->AnzKreuz*sizeof(struct Kreuzung));
      } else {
    	errmessage(16);
      }
      for (i=0; i < input->AnzKreuz;i++)
      {
    	input->Kreuz[i].KopelArt=1;
      }
    }       
    if ( strcmp( ZEILE_READ, "readrelativethicknesses(on=1/off=0)") == 0 )  
    {
      sscanf(ZEIGER, "%d", &input->readRelDick);
    }
    if ( strcmp( ZEILE_READ, "fuselagechambertwist") == 0 ) 
    {
      sscanf(ZEIGER, "%lf", &input->AlfaPlusRumpf);
    }
    if ( strncmp( ZEILE_READ , "wing", 4) == 0)
    {
      lauf = 0;
      zahl = 0;
      do 
      {
        if (zahl==0)
        {
          zahl = ZEILE_READ[lauf+4] - 48; 
        } else {
          zahl=zahl*10;
          zahl = zahl + ZEILE_READ[lauf+4] - 48; 
        }
        lauf++;
      }  
      while ((ZEILE_READ[lauf+4] >= 48) && (ZEILE_READ[lauf+4] <= 57));
      if (zahl <= input->NbWings)
      {
        ZEIGER2=&ZEILE_READ[lauf+4];
        if ( strcmp( ZEIGER2, "numberofsections") == 0 )
        {
          sscanf(ZEIGER, "%d", &input->Fluegel[zahl-1].anzahlSchnitte); 
        }
        if ( strcmp( ZEIGER2, "numberofpanelchord") == 0 )
        {
          sscanf(ZEIGER, "%d", &input->Fluegel[zahl-1].AnzahlTiefe);
        }
        if ( strcmp( ZEIGER2, "extendtosymmetry(on=1/off=0)") == 0 )
        {
          sscanf(ZEIGER, "%d", &input->Fluegel[zahl-1].SymExt);
        }
        if ( strcmp( ZEIGER2, "sym.boundarycondition(on=1/off=0)") == 0 )
        {
          sscanf(ZEIGER, "%d", &input->Fluegel[zahl-1].SymRandBed);
        }
        if ( strcmp( ZEIGER2, "createsymmetry(off=0/on=1/newwing=2)") == 0 )
        {
          sscanf(ZEIGER, "%d", &input->Fluegel[zahl-1].Spiegeln);        
        }
        if ( strcmp( ZEIGER2, "FileSurfix") == 0 )
        {
          sscanf(ZEIGER, "%s", &input->Fluegel[zahl-1].Surfix[0]);
        }
        if ( strcmp( ZEIGER2, "FilePrefix") == 0 )
        {
          sscanf(ZEIGER, "%s", &input->Fluegel[zahl-1].Prefix[0]);
        }
        if ( strcmp( ZEIGER2, "zones") == 0 )
	{
          char *zeig, *zeig2;
          int lauf5,lauf52;
	  
	  lauf5 = 1;
	  zeig = NULL;
	  zeig = strchr(ZEIGER, (int)'{');
	  if (zeig==NULL)
	  {
	     errmessage(37);
	  }
	  if (strchr(zeig, (int)'}')==NULL) 
	  {
	     errmessage(37);
	  }
	  zeig2 = strtok((zeig), "_} ");
	  //printf("ZEIG2: %s\n", (zeig2));
	  if (input->BSURF_FLUEGEL==NULL)
	  {
            input->BSURF_FLUEGEL=(struct BSURF_LISTZONE *)malloc(input->NbWings*sizeof(struct BSURF_LISTZONE));   
          }
	  input->BSURF_FLUEGEL[zahl-1].Zone=(int *)malloc(2*sizeof(int));
          do 
	  {
            input->BSURF_FLUEGEL[zahl-1].Zone[lauf5]=0;
	    lauf52=1;
	    do
            {  
	      if (lauf52 == 1)
	      {
	        if (lauf5==1)
		{
		  input->BSURF_FLUEGEL[zahl-1].Zone[lauf5]=input->BSURF_FLUEGEL[zahl-1].Zone[lauf5]+zeig2[1]-48;
		  lauf52++; 
		}else{
		  input->BSURF_FLUEGEL[zahl-1].Zone[lauf5]=input->BSURF_FLUEGEL[zahl-1].Zone[lauf5]+zeig2[0]-48;           
	        }
	      }else{
	        input->BSURF_FLUEGEL[zahl-1].Zone[lauf5]=input->BSURF_FLUEGEL[zahl-1].Zone[lauf5]*10;
		input->BSURF_FLUEGEL[zahl-1].Zone[lauf5]=input->BSURF_FLUEGEL[zahl-1].Zone[lauf5]+zeig2[lauf52-1]-48; 
	      }
              lauf52++;
	      //printf("NUMB: %d Strlen: %d\n",input->BSURF_FLUEGEL[zahl-1].Zone[lauf5] ,strlen(zeig2));
	    } while ((int)strlen(zeig2) >= lauf52);
	    zeig2 = strtok(NULL, "_} ");
	    //printf("ZEIG2: %s\n", (zeig2));
	    if (zeig2!=NULL)
	    {
	      lauf5++;
	      input->BSURF_FLUEGEL[zahl-1].Zone=(int *)realloc(input->BSURF_FLUEGEL[zahl-1].Zone, (lauf5+2)*sizeof(int));
	    }
          } while  (zeig2!=NULL);
	  input->BSURF_FLUEGEL[zahl-1].Zone[0]=lauf5;
        }
      }
    }
    if ( strncmp( ZEILE_READ , "junction", 8) == 0)
    {
      lauf = 0;
      zahl = 0;
      do 
      {  
        if (zahl ==0)
        {
          zahl = ZEILE_READ[lauf+8] - 48; 
        } else {
          zahl=zahl*10;
          zahl = zahl + ZEILE_READ[lauf+8] - 48; 
        }
        lauf++;
      }  
      while ((ZEILE_READ[lauf+8] >= 48) && (ZEILE_READ[lauf+8] <= 57));
      /*printf("Zah Kreuz einlsen: %d \n",zahl);*/
      if (zahl <= input->AnzKreuz)
      {
        ZEIGER2=&ZEILE_READ[lauf+8];
        if ( strcmp( ZEIGER2, "wing1") == 0 )
        {
          sscanf(ZEIGER, "%d", &input->Kreuz[zahl-1].KreuzFl1); 
        }
        if ( strcmp( ZEIGER2, "wing2") == 0 )
        {
          sscanf(ZEIGER, "%d", &input->Kreuz[zahl-1].KreuzFl2); 
        }
        if ( strcmp( ZEIGER2, "wing1element") == 0 )
        {
          sscanf(ZEIGER, "%d", &input->Kreuz[zahl-1].KreuzPosFl1); 
        }
        if ( strcmp( ZEIGER2, "wing2element") == 0 )
        {
          sscanf(ZEIGER, "%d", &input->Kreuz[zahl-1].KreuzPosFl2); 
        }
        if ( strcmp( ZEIGER2, "wing1chordstartpos.") == 0 )
        {
          sscanf(ZEIGER, "%d", &input->Kreuz[zahl-1].KreuzTiefFl1); 
        }
        if ( strcmp( ZEIGER2, "wing2chordstartpos.") == 0 )
        {
          sscanf(ZEIGER, "%d", &input->Kreuz[zahl-1].KreuzTiefFl2); 
        }
        if ( strcmp( ZEIGER2, "wing1junctionconstraint") == 0 )
        {
          sscanf(ZEIGER, "%d", &input->Kreuz[zahl-1].GabBedFl1); 
        }
        if ( strcmp( ZEIGER2, "wing2junctionconstraint") == 0 )
        {
          sscanf(ZEIGER, "%d", &input->Kreuz[zahl-1].GabBedFl2); 
        }
        if ( strcmp( ZEIGER2, "typofjunctionconstrain") == 0 )
        {
          sscanf(ZEIGER, "%d", &input->Kreuz[zahl-1].KopelArt); 
        }	
      }
    }      
    //Klappen
    if ( strcmp( ZEILE_READ, "numberofflaps") == 0 ) 
    {
      int a;
      if (input->Klappen == NULL)
      {
        sscanf(ZEIGER, "%d", &input->AnzahlKlappen);
        input->Klappen = (struct EineKlappen *)malloc(input->AnzahlKlappen*sizeof(struct EineKlappen));
        for (a=0; a<input->AnzahlKlappen; a++)
        {
          input->Klappen[a].AnzPanel=0;
          input->Klappen[a].Spiegeln=0;
        }
      } else {
        errmessage(50);
      }
    }     
    if ( strncmp( ZEILE_READ , "flap", 4) == 0)
    {
      lauf = 0;
      zahl = 0;
      do 
      {  
        if (zahl ==0)
        {
          zahl = ZEILE_READ[lauf+4] - 48; 
        } else {
          zahl=zahl*10;
          zahl = zahl + ZEILE_READ[lauf+4] - 48; 
        }
        lauf++;
      }  
      while ((ZEILE_READ[lauf+4] >= 48) && (ZEILE_READ[lauf+4] <= 57));
      //printf("Zahl Klappe einlsen: %d \n",zahl);
      if (zahl <= input->AnzahlKlappen)
      {
        ZEIGER2=&ZEILE_READ[lauf+4];
        if ( strcmp( ZEIGER2, "wing") == 0 )
        {
          sscanf(ZEIGER, "%d", &input->Klappen[zahl-1].Fluegel); 
        }
        if ( strcmp( ZEIGER2, "section1") == 0 )
        {
          sscanf(ZEIGER, "%d", &input->Klappen[zahl-1].Schnitt1); 
        }
        if ( strcmp( ZEIGER2, "section2") == 0 )
        {
          sscanf(ZEIGER, "%d", &input->Klappen[zahl-1].Schnitt2); 
        }
        if ( strcmp( ZEIGER2, "section1relativflapdepth") == 0 )
        {
          sscanf(ZEIGER, "%lf", &input->Klappen[zahl-1].eta1); 
        }
        if ( strcmp( ZEIGER2, "section2relativflapdepth") == 0 )
        {
          sscanf(ZEIGER, "%lf", &input->Klappen[zahl-1].eta2); 
        }
        if ( strcmp( ZEIGER2, "deflectionangle[deg]") == 0 )
        {
          sscanf(ZEIGER, "%lf", &input->Klappen[zahl-1].winkel); 
        }
        if ( strcmp( ZEIGER2, "numberofpanels") == 0 )
        {
          sscanf(ZEIGER, "%d", &input->Klappen[zahl-1].AnzPanel); 
        }      
      }
    }      
    //Mischer
    
    if ( strcmp( ZEILE_READ, "devvariable") == 0 ) 
    {
      char DUMMY[11]= {"\0"};
      zahl = 0;
      zaehler2 = 0; 
      do 
      {
        DUMMY[zahl]  = ZEIGER[zahl];
        zahl++;
      }
      while ((ZEIGER[zahl] != '=') && (ZEIGER[zahl] != '\0'));
      ZEIGER=&ZEIGER[zahl+1];
      
      if ( strncmp( DUMMY , "x", 1) != 0)
      {
        errmessage(65);
      }
      ZEIGER2=&DUMMY[1];
      sscanf(ZEIGER2, "%d", &NewVar);
      if (NewVar>IstVar)
      {
        if (IstVar==0)
        {
          input->MixVar = (double *)malloc(NewVar*sizeof(double));
          for(lauf=0;lauf<NewVar;lauf++)
          {
             input->MixVar[NewVar-1]=0;
          }
          IstVar=NewVar;
        }else{
          input->MixVar = (double *)realloc(input->MixVar, NewVar*sizeof(double));
          for(lauf=IstVar-1;lauf<NewVar;lauf++)
          {
             input->MixVar[NewVar-1]=0.0;
          }
          IstVar=NewVar;
        }
      }
      sscanf( ZEIGER, "%lf", &input->MixVar[NewVar-1]);      
      input->AnzahlMixVar=IstVar;
    }
        
    if ( strcmp( ZEILE_READ, "numberofmixersrules") == 0 ) 
    {
       sscanf(ZEIGER, "%d", &input->AnzahlMischer);
       input->Mischer = (struct EinMischer *)malloc(input->AnzahlMischer*sizeof(struct EinMischer));
    }
    if ( strncmp( ZEILE_READ , "mixer", 5) == 0)
    {
      lauf = 0;
      zahl = 0;
      do 
      {  
        if (zahl ==0)
        {
          zahl = ZEILE_READ[lauf+5] - 48; 
        } else {
          zahl=zahl*10;
          zahl = zahl + ZEILE_READ[lauf+5] - 48; 
        }
        lauf++;
      }  
      while ((ZEILE_READ[lauf+5] >= 48) && (ZEILE_READ[lauf+5] <= 57));
      if ((input->debuglevel==3) || (input->debuglevel>99))
      {
	printf("\nRegelLese:%s   Zahl:%d\n",ZEILE_READ,zahl);
      }
      if (zahl <= input->AnzahlMischer)
      {
        ZEIGER2=&ZEILE_READ[lauf+5];
        if ((input->debuglevel==3) || (input->debuglevel>99))
        {
          printf("RegelLese:%s\n",ZEIGER);
        }
        if ( strcmp( ZEIGER2, "rule") == 0 )
        {
          sscanf(ZEIGER, "%s", &input->Mischer[zahl-1].Regel[0]);
        }
      }  
    }
    //Parametric-Variations  
    if ( strcmp( ZEILE_READ, "parametervariation(off=0/directvalues=1/spline=2/flap=3/flapfunc=4)") == 0 ) 
    {
      input->NumberOfParametSetsinput++;
      if (input->NumberOfParametSetsinput==1)
      { 
        input->ParaVar=(struct PARAvariation *)malloc(input->NumberOfParametSetsinput*sizeof(struct PARAvariation));
      }else{
        input->ParaVar=(struct PARAvariation *)realloc(input->ParaVar, input->NumberOfParametSetsinput*sizeof(struct PARAvariation));      
      }
      input->ParaVar[input->NumberOfParametSetsinput-1].Paraset4eachCLA=0;
      input->ParaVar[input->NumberOfParametSetsinput-1].AbsOrDiffer=0;
      input->ParaVar[input->NumberOfParametSetsinput-1].link=0;
      sscanf(ZEIGER, "%d", &input->ParaVar[input->NumberOfParametSetsinput-1].parametric);
    }  
    if ( strcmp( ZEILE_READ, "numberofparametricvariations") == 0 ) 
    {
      sscanf(ZEIGER, "%d", &input->ParaVar[input->NumberOfParametSetsinput-1].NrVariations);
    }  
    if ( strcmp( ZEILE_READ, "nameofparametricinputfile") == 0 ) 
    {
      sscanf(ZEIGER, "%s", &input->ParaVar[input->NumberOfParametSetsinput-1].PARAINP[0]);
    }   
    if ( strcmp( ZEILE_READ, "seperatesetofparameterforeachcaoralfa(0=off/1=on)") == 0 ) 
    {
      sscanf(ZEIGER, "%d", &input->ParaVar[input->NumberOfParametSetsinput-1].Paraset4eachCLA);
    }           
    if ( strcmp( ZEILE_READ, "valueabsolutordiffrence(absolut=0/diffrence=1)") == 0 ) 
    {
      sscanf(ZEIGER, "%d", &input->ParaVar[input->NumberOfParametSetsinput-1].AbsOrDiffer);
    }             
    if ( strcmp( ZEILE_READ, "parametersforoptimization") == 0 ) 
    {    
      char *zeig, *zeig2;
      struct tstring { char teilstring[50]; };
      struct tstring *pointstring;
      int posClose, posPoint, raus;
      input->ParaVar[input->NumberOfParametSetsinput-1].anzVar=0;
      raus=0;
      zeig=NULL;
      do 
      {
        zeig = strchr(ZEIGER, (int)'{');
        if (zeig == NULL)
        {
          lauf=0;
          for (lauf=0;lauf <= MaxZeilenLaenge;lauf++)
          {
            ZEILE_READ[lauf] = '\0';
          }
          fgets (ZEILE_READ,MaxZeilenLaenge,fpara);
          ZEIGER=&ZEILE_READ[0];
        }
      }
      while ((zeig == NULL) && (!feof(fpara)));
      if (feof(fpara)) 
      {
        errmessage(17);
      }
      posClose = strcspn( zeig, "}" );
      posPoint = strcspn( zeig, "." );
      //printf("posClose:%d ; posPoint:%d\n",posClose,posPoint);
      if (posClose == posPoint)
      {
        do
        {
          fgets (ZEILE_READ,MaxZeilenLaenge,fpara);
          posClose = strcspn( ZEILE_READ, "}" );
          posPoint = strcspn( ZEILE_READ, "." );
        }while ((posClose!=posPoint) && (!feof(fpara))); 
      }
      if ((posClose < posPoint) )
      {
        printf ("NoDataInput\n");
        raus=1;
      } else if  (feof(fpara)) {
        errmessage(17);
      } else if ((posClose > posPoint))
      {
        zeig=(zeig+1);
        do
        {
          zeig2 = strtok((zeig), "_ \n\t");
          do
          {
            posClose = strcspn( zeig2, "}" );
            posPoint = strcspn( zeig2, "." );
            if ((*(zeig2)==(int)'}') || (posClose < posPoint))
            {
              raus=1;
            }else 
            {
              posClose = strcspn( zeig2, "}" );
              posPoint = strcspn( zeig2, "\0" );
              if (posClose < posPoint)
              {
                 raus=1;
                 zeig2 = strtok((zeig2), "}");
              }
              
              if (input->ParaVar[input->NumberOfParametSetsinput-1].anzVar==0)
              {
                pointstring=(struct tstring *)malloc(sizeof(struct tstring));           
              }else {
                pointstring=(struct tstring *)realloc(pointstring, (input->ParaVar[input->NumberOfParametSetsinput-1].anzVar+1)*sizeof(struct tstring));                                 
              }              
              strcpy(pointstring[input->ParaVar[input->NumberOfParametSetsinput-1].anzVar].teilstring, zeig2);
              input->ParaVar[input->NumberOfParametSetsinput-1].anzVar++;                             
            }                               
            if (raus==0){                         
              zeig2 = strtok(NULL, "_ \n\t");                
            }                               
          }                                 
          while ((zeig2!= NULL) && (raus==0));               
          lauf=0;                              
          for (lauf=0;lauf <= MaxZeilenLaenge;lauf++)            
          {                                 
            ZEILE_READ[lauf] = '\0';                   
          } 
          fgets (ZEILE_READ,MaxZeilenLaenge,fpara);        
          zeig=&ZEILE_READ[0];                  
        }
        while ((raus!=1) && (!feof(fpara)));
      }
      if (input->ParaVar[input->NumberOfParametSetsinput-1].anzVar >0)
      {
        input->ParaVar[input->NumberOfParametSetsinput-1].VarParameter=(struct VariationsPunkt *)malloc(input->ParaVar[input->NumberOfParametSetsinput-1].anzVar*sizeof(struct VariationsPunkt));
      }
      for (lauf=0; lauf < input->ParaVar[input->NumberOfParametSetsinput-1].anzVar; lauf++)
      {      
        input->ParaVar[input->NumberOfParametSetsinput-1].VarParameter[lauf].twist=OFF;
        input->ParaVar[input->NumberOfParametSetsinput-1].VarParameter[lauf].X=OFF;
        input->ParaVar[input->NumberOfParametSetsinput-1].VarParameter[lauf].Y=OFF;
        input->ParaVar[input->NumberOfParametSetsinput-1].VarParameter[lauf].Z=OFF;
        input->ParaVar[input->NumberOfParametSetsinput-1].VarParameter[lauf].Vloc=OFF;
        zeig=pointstring[lauf].teilstring;
        zeig2 = strtok(zeig, ".");
        sscanf(zeig2,"%d", &input->ParaVar[input->NumberOfParametSetsinput-1].VarParameter[lauf].Fluegel);
        zeig2 = strtok(NULL, ".");
        sscanf(zeig2,"%d", &input->ParaVar[input->NumberOfParametSetsinput-1].VarParameter[lauf].Schnitt);
        input->ParaVar[input->NumberOfParametSetsinput-1].AnzahlParameter=0;
        do 
        {
          if ((zeig2 = strtok(NULL, ".")) != NULL)
          {         
            if (strcmp(zeig2,"tw")==0){ input->ParaVar[input->NumberOfParametSetsinput-1].VarParameter[lauf].twist=ON; input->ParaVar[input->NumberOfParametSetsinput-1].AnzahlParameter++;}   
            if (strcmp(zeig2,"x")==0) { input->ParaVar[input->NumberOfParametSetsinput-1].VarParameter[lauf].X=ON; input->ParaVar[input->NumberOfParametSetsinput-1].AnzahlParameter++;}     
            if (strcmp(zeig2,"y")==0) { input->ParaVar[input->NumberOfParametSetsinput-1].VarParameter[lauf].Y=ON; input->ParaVar[input->NumberOfParametSetsinput-1].AnzahlParameter++;}     
            if (strcmp(zeig2,"z")==0) { input->ParaVar[input->NumberOfParametSetsinput-1].VarParameter[lauf].Z=ON; input->ParaVar[input->NumberOfParametSetsinput-1].AnzahlParameter++;}     
            if (strcmp(zeig2,"depth")==0) {input->ParaVar[input->NumberOfParametSetsinput-1].VarParameter[lauf].tiefe=ON; input->ParaVar[input->NumberOfParametSetsinput-1].AnzahlParameter++;}
          }
        }while(zeig2!= NULL);
      }
    }    
    
    if ( strcmp( ZEILE_READ, "parametersforsplineoptimization") == 0 )
    {
      char *zeig, *zeig2;
      struct tstring { char teilstring[50]; };
      struct tstring *pointstring;
      int posClose, posPoint, raus;
      
      input->ParaVar[input->NumberOfParametSetsinput-1].anzSpVar=0;
      /*printf("Jetzt aber\n");*/
      raus=0;
      zeig=NULL;
      do 
      {
        zeig = strchr(ZEIGER, (int)'{');
        if (zeig == NULL)
        {
          lauf=0;
          for (lauf=0;lauf <= MaxZeilenLaenge;lauf++)
          {
            ZEILE_READ[lauf] = '\0';
          }
          fgets (ZEILE_READ,MaxZeilenLaenge,fpara);
          ZEIGER=&ZEILE_READ[0];
        }
      }
      while ((zeig == NULL) && (!feof(fpara)));
      if (feof(fpara)) 
      {
        errmessage(17);
      }
      posClose = strcspn( zeig, "}" );
      posPoint = strcspn( zeig, "." );
      //printf("posClose:%d ; posPoint:%d\n",posClose,posPoint);
      if (posClose == posPoint)
      {
        do
        {
          fgets (ZEILE_READ,MaxZeilenLaenge,fpara);
          posClose = strcspn( ZEILE_READ, "}" );
          posPoint = strcspn( ZEILE_READ, "." );
        }while ((posClose!=posPoint) && (!feof(fpara))); 
      }
      if ((posClose < posPoint) )
      {
        printf ("NoDataInput\n");
        raus=1;
      } else if  (feof(fpara)) {
        errmessage(17);
      } else if ((posClose > posPoint))
      {
        zeig=(zeig+1);
        do
        {
          zeig2 = strtok((zeig), "_ \n\t");
          do
          {
            posClose = strcspn( zeig2, "}" );
            posPoint = strcspn( zeig2, "." );
            if ((*(zeig2)==(int)'}') || (posClose < posPoint))
            {
              raus=1;
            }else 
            {
              posClose = strcspn( zeig2, "}" );
              posPoint = strcspn( zeig2, "\0" );
              if (posClose < posPoint)
              {
                 raus=1;
                 zeig2 = strtok((zeig2), "}");
              }              
              if (input->ParaVar[input->NumberOfParametSetsinput-1].anzSpVar==0)
              {
                pointstring=(struct tstring *)malloc(sizeof(struct tstring));
              }else {
                pointstring=(struct tstring *)realloc(pointstring, (input->ParaVar[input->NumberOfParametSetsinput-1].anzSpVar+1)*sizeof(struct tstring));
              }
              strcpy(pointstring[input->ParaVar[input->NumberOfParametSetsinput-1].anzSpVar].teilstring, zeig2);   
              input->ParaVar[input->NumberOfParametSetsinput-1].anzSpVar++;                             
            }                               
            if (raus==0){                         
              zeig2 = strtok(NULL, "_ \n\t");                
            }                               
          }                                 
          while ((zeig2!= NULL) && (raus==0));               
          lauf=0;                              
          for (lauf=0;lauf <= MaxZeilenLaenge;lauf++)            
          {                                 
            ZEILE_READ[lauf] = '\0';                   
          } 
          fgets (ZEILE_READ,MaxZeilenLaenge,fpara);        
          zeig=&ZEILE_READ[0];                  
        }
        while ((raus!=1) && (!feof(fpara)));
      } 
      if (input->ParaVar[input->NumberOfParametSetsinput-1].anzSpVar >0)
      {
        input->ParaVar[input->NumberOfParametSetsinput-1].SplineVarPunkt=(struct SplineVariationsPunkt *)malloc(input->ParaVar[input->NumberOfParametSetsinput-1].anzSpVar*sizeof(struct SplineVariationsPunkt));
      }
      for (lauf=0; lauf < input->ParaVar[input->NumberOfParametSetsinput-1].anzSpVar; lauf++)
      {      
        input->ParaVar[input->NumberOfParametSetsinput-1].SplineVarPunkt[lauf].Variable.twist=OFF;
        input->ParaVar[input->NumberOfParametSetsinput-1].SplineVarPunkt[lauf].Variable.X=OFF; 
        input->ParaVar[input->NumberOfParametSetsinput-1].SplineVarPunkt[lauf].Variable.Y=OFF;
        input->ParaVar[input->NumberOfParametSetsinput-1].SplineVarPunkt[lauf].Variable.Z=OFF; 
        input->ParaVar[input->NumberOfParametSetsinput-1].SplineVarPunkt[lauf].Variable.tiefe=OFF;
        input->ParaVar[input->NumberOfParametSetsinput-1].SplineVarPunkt[lauf].Variable.Vloc=OFF;
        
        //printf("Teilstring:%s\n",pointstring[lauf].teilstring);
        zeig=pointstring[lauf].teilstring;
        zeig2 = strtok(zeig, ".");
        sscanf(zeig2,"%d", &input->ParaVar[input->NumberOfParametSetsinput-1].SplineVarPunkt[lauf].Fluegel);
        zeig2 = strtok(NULL, ".");
        sscanf(zeig2,"%d", &input->ParaVar[input->NumberOfParametSetsinput-1].SplineVarPunkt[lauf].StartSchnitt);
        zeig2 = strtok(NULL, ".");
        sscanf(zeig2,"%d", &input->ParaVar[input->NumberOfParametSetsinput-1].SplineVarPunkt[lauf].EndSchnitt);
        zeig2 = strtok(NULL, ".");
        if (strcmp(zeig2,"tw")==0){ input->ParaVar[input->NumberOfParametSetsinput-1].SplineVarPunkt[lauf].Variable.twist=ON; }   
        if (strcmp(zeig2,"x")==0) { input->ParaVar[input->NumberOfParametSetsinput-1].SplineVarPunkt[lauf].Variable.X=ON; }	  
        if (strcmp(zeig2,"y")==0) { input->ParaVar[input->NumberOfParametSetsinput-1].SplineVarPunkt[lauf].Variable.Y=ON;}	 
        if (strcmp(zeig2,"z")==0) { input->ParaVar[input->NumberOfParametSetsinput-1].SplineVarPunkt[lauf].Variable.Z=ON; }	  
        if (strcmp(zeig2,"tiefe")==0) {input->ParaVar[input->NumberOfParametSetsinput-1].SplineVarPunkt[lauf].Variable.tiefe=ON; }
        zeig2 = strtok(NULL, ".");
        sscanf(zeig2,"%d", &input->ParaVar[input->NumberOfParametSetsinput-1].SplineVarPunkt[lauf].AnzahlStuetzpunkte);
        /*printf("anzSpVar: %d , Fl:%d, StS: %d EnS:%d AnzSt: %d\n ", lauf,input->ParaVar[input->NumberOfParametSetsinput-1].SplineVarPunkt[lauf].Fluegel,input->ParaVar[input->NumberOfParametSetsinput-1].SplineVarPunkt[lauf].StartSchnitt
        ,input->ParaVar[input->NumberOfParametSetsinput-1].SplineVarPunkt[lauf].EndSchnitt,input->ParaVar[input->NumberOfParametSetsinput-1].SplineVarPunkt[lauf].AnzahlStuetzpunkte);*/
      }
      //printf("anzSpVar: %d \n ", input->ParaVar[input->NumberOfParametSetsinput-1].anzSpVar);
    }

    //Trimmen
    if ( strcmp( ZEILE_READ, "trimming(on=1/off=0)") == 0 ) 
    {
      sscanf(ZEIGER, "%d", &input->Austrimmen);
    }  
    if ( strcmp( ZEILE_READ, "maxnumberofiterrationen") == 0 ) 
    {
      sscanf(ZEIGER, "%d", &input->MaxCmIter);
    }  
    if ( strcmp( ZEILE_READ, "trimsurfacewingnumber") == 0 ) 
    {
      sscanf(ZEIGER, "%d", &input->CmSteuerFl);
    }  
    if ( strcmp( ZEILE_READ, "trimsurfaceflapnumber") == 0 ) 
    {
      sscanf(ZEIGER, "%d", &input->CmSteuerKlappe);
    }  
    if ( strcmp( ZEILE_READ, "Epsilonlogfile") == 0 ) 
    {
      sscanf(ZEIGER, "%s", &input->EpsAusgabe[0]);
    }
    if ( strcmp( ZEILE_READ, "targetCm") == 0 ) 
    {
      sscanf(ZEIGER, "%lf", &input->ZielCm);
    }  
    if ( strcmp( ZEILE_READ, "alfarotadjustment") == 0 ) 
    {
      sscanf(ZEIGER, "%d", &input->ArotAdjust);
    }  
    if ( strcmp( ZEILE_READ, "InitTrimStepSize") == 0 )
    {
      sscanf(ZEIGER, "%lf", &input->InitStepSize);
    }
    

    //Post-Processing      
    if ( strcmp( ZEILE_READ, "IBMreferencesection") == 0 ) 		     
    {									     
      char *zeig, *zeig2, dummy[19];							   
      sscanf(ZEIGER, "%s", &dummy[0]);			 
      zeig = strtok(dummy, ".");
      zeig2 = strtok(NULL, ".");					     	  
      sscanf(zeig,"%d", &input->IBM_BD_WING);
      sscanf(zeig2,"%d", &input->IBM_BD_SEC);				     	  
    }  
    if ( strcmp( ZEILE_READ, "WRBMdefinitionpoint") == 0 ) 		     
    {									     
      char *zeig, *zeig2, dummy[19];							   
      sscanf(ZEIGER, "%s", &dummy[0]);			 
      zeig = strtok(dummy, ".");
      zeig2 = strtok(NULL, ".");					     	  
      sscanf(zeig,"%d", &input->WRBM_BD_WING);
      sscanf(zeig2,"%d", &input->WRBM_BD_SEC);				     	  
    }  

    if ( strcmp( ZEILE_READ, "referencespeed[m/s]") == 0 ) 
    {
      sscanf(ZEIGER, "%lf", &input->RefSpeed);
    }  
    if ( strcmp( ZEILE_READ, "referencedensity[kg/m3]") == 0 ) 
    {
      sscanf(ZEIGER, "%lf", &input->RefDensity);
    }
    if ( strcmp( ZEILE_READ, "referencetemperature(Kelvin)") == 0 ) 
    {
      sscanf(ZEIGER, "%lf", &input->RefTemp);
    }  
    if ( strcmp( ZEILE_READ, "IBMreferencevalue") == 0 ) 
    {
      sscanf(ZEIGER, "%lf", &input->BasisIntBiegeMoment);
    }  
    if ( strcmp( ZEILE_READ, "IDBMreferencevalue") == 0 ) 
    {
      sscanf(ZEIGER, "%lf", &input->BasisIntDickenBiegeMoment);
    } 
    if ( strcmp( ZEILE_READ, "makeBSURFinput(on=1/off=0)") == 0 ) 
    {
      sscanf(ZEIGER, "%d", &input->makeBsurfinp);
    }	
    if ( strcmp( ZEILE_READ, "readBSURFloaddist.(on=1/off=0)") == 0 ) 
    {
      sscanf(ZEIGER, "%d", &input->TAUbsurfINP);
    }
    if ( strcmp( ZEILE_READ, "numberofzones") == 0 ) 
    {
      sscanf(ZEIGER, "%d", &input->NUM_OF_ZONES);
    }	
    if ( strcmp( ZEILE_READ, "computesquarediffrence(on=1/off=0)") == 0 ) 
    {
      sscanf(ZEIGER, "%d", &input->BsurfSquareDiff);
    }	

    //Die Methode    
    if ( strcmp( ZEILE_READ, "orderofalfainfluence(max.9)") == 0 ) 
    {
      sscanf(ZEIGER, "%d", &input->kHOrd);
      if (input->kHOrd > 9)
      {
        errmessage(19);
      }
      if (input->kHOrd > 0)
      {
        input->kH = (double *)malloc(input->kHOrd*sizeof(double));
      }
    }  
    if (( strcmp( ZEILE_READ, "kH1") == 0 ) && (input->kHOrd != 0)) 
    {
      sscanf(ZEIGER, "%lf", &input->kH[0]);
    }  
    if (( strcmp( ZEILE_READ, "kH2") == 0 ) && (input->kHOrd != 0)) 
    {
      sscanf(ZEIGER, "%lf", &input->kH[1]);
    }  
    if (( strcmp( ZEILE_READ, "kH3") == 0 ) && (input->kHOrd != 0)) 
    {
      sscanf(ZEIGER, "%lf", &input->kH[2]);
    }  
    if (( strcmp( ZEILE_READ, "kH4") == 0 ) && (input->kHOrd != 0)) 
    {
      sscanf(ZEIGER, "%lf", &input->kH[3]);
    }  
    if (( strcmp( ZEILE_READ, "kH5") == 0 ) && (input->kHOrd != 0)) 
    {
      sscanf(ZEIGER, "%lf", &input->kH[4]);
    }  
    if (( strcmp( ZEILE_READ, "kH6") == 0 ) && (input->kHOrd != 0)) 
    {
      sscanf(ZEIGER, "%lf", &input->kH[5]);
    }  
    if (( strcmp( ZEILE_READ, "kH7") == 0 ) && (input->kHOrd != 0)) 
    { 
      sscanf(ZEIGER, "%lf", &input->kH[6]);
    }  
    if (( strcmp( ZEILE_READ, "kH8") == 0 ) && (input->kHOrd != 0)) 
    {
      sscanf(ZEIGER, "%lf", &input->kH[7]);
    }  
    if (( strcmp( ZEILE_READ, "kH9") == 0 ) && (input->kHOrd != 0)) 
    {
      sscanf(ZEIGER, "%lf", &input->kH[8]);
    }  
    if ( strcmp( ZEILE_READ, "referenceCDi") == 0 )
    {
      sscanf(ZEIGER, "%lf", &input->BCwi);
    }  
    if ( strcmp( ZEILE_READ, "referenceCl") == 0 )
    {
      sscanf(ZEIGER, "%lf", &input->BWRBM);
    }  
    if ( strcmp( ZEILE_READ, "referenceCD0") == 0 )
    {
      sscanf(ZEIGER, "%lf", &input->BCw0);
    }  
    if ( strcmp( ZEILE_READ, "referenceunrolledarea") == 0 )
    {
      sscanf(ZEIGER, "%lf", &input->UrArea);
    }  
    if ( strcmp( ZEILE_READ, "referencealfa") == 0 )
    {
      sscanf(ZEIGER, "%lf", &input->BAlafa);
    }  

    //Struct-Kopplung
    if ( strcmp( ZEILE_READ, "struct.coupling(on=1/off=0)") == 0 ) 
    {
      sscanf(ZEIGER, "%d", &input->Struckt_Coup_ON);
    }  
    if ( strcmp( ZEILE_READ, "struct.loadcase(on=1/off=0)") == 0 ) 
    {
      sscanf(ZEIGER, "%d", &input->Load_Case_ON);
    }  
    if ( strcmp( ZEILE_READ, "struct.numberofitterations") == 0 ) 
    {
      sscanf(ZEIGER, "%d", &input->Numb_Iter);
    }  
    if ( strcmp( ZEILE_READ, "logiterations(on=1/off=0)") == 0 ) 
    {
      sscanf(ZEIGER, "%d", &input->Log_Iter_ON);
    }  
    if ( strcmp( ZEILE_READ, "loadreferencesection(wing.section)") == 0 ) 
    {
      char *zeig, *zeig2;
      sscanf(ZEIGER, "%s", &input->ProjektsionsSchnitt[0]);
      zeig = strtok(input->ProjektsionsSchnitt, ".");
      zeig2 = strtok(NULL, ".");
      sscanf(zeig, "%d", &input->ProjektionsFluegel);
      sscanf(zeig2, "%d", &input->ProjektionsSchnitt);
    }  
    if ( strcmp( ZEILE_READ, "structprogrammworkdirectory") == 0 ) 
    {
      sscanf(ZEIGER, "%s", &input->StructProg_EXEdir[0]);
    }  
    if ( strcmp( ZEILE_READ, "structprogrammexecutebal") == 0 ) 
    {
      sscanf(ZEIGER, "%s", &input->StructProg_EXEFile[0]);
    }  
    if ( strcmp( ZEILE_READ, "structprogrammexecutebal2") == 0 ) 
    {
      sscanf(ZEIGER, "%s", &input->StructProg_EXEFile2[0]);
    }  
    if ( strcmp( ZEILE_READ, "resultsubdir") == 0 ) 
    {
      sscanf(ZEIGER, "%s", &input->StructProg_ResultSubDir[0]);
    }  
    if ( strcmp( ZEILE_READ, "LIFTINGLINEworkdirectory") == 0 ) 
    {
      sscanf(ZEIGER, "%s", &input->LILIDIR[0]);
    }  
    if ( strcmp( ZEILE_READ, "structBEAMTXT") == 0 ) 
    {
      sscanf(ZEIGER, "%s", &input->StructProg_ResultFile[0]);
    }  
    if ( strcmp( ZEILE_READ, "structBEAMJIG2FLIGHTTXT") == 0 ) 
    {
      sscanf(ZEIGER, "%s", &input->StructProg_Jig2FlightFile[0]);
    }                                                                                                         
    if ( strcmp( ZEILE_READ, "additionaltargetCL") == 0 )                                                     
    {
      sscanf(ZEIGER, "%lf", &input->Target_CL_LoadCase);
    }                                                                                                 
    if ( strcmp( ZEILE_READ, "additionalreferencespeed[m/s]") == 0 )                             
    {
      sscanf(ZEIGER, "%lf", &input->Ref_Density_LoadCase);
    }  
    if ( strcmp( ZEILE_READ, "additionalreferencedensity[kg/m3]") == 0 )
    {
      sscanf(ZEIGER, "%lf", &input->Ref_Speed_LoadCase);
    }  
    if ( strcmp( ZEILE_READ, "offsetstruct.modeltoLiftinglinegeometrydeltaX[m]") == 0 )   
    {
      sscanf(ZEIGER, "%lf", &input->deltaStructX);
    }  
    if ( strcmp( ZEILE_READ, "offsetstruct.modeltoLiftinglinegeometrydeltaZ[m]") == 0 )
    {
      sscanf(ZEIGER, "%lf", &input->deltaStructZ);
    }  
    if ( strcmp( ZEILE_READ, "structmodelreferencespan[m]") == 0 )
    {
      sscanf(ZEIGER, "%lf", &input->StructModBezSpan);
    }  
    if ( strcmp( ZEILE_READ, "useJIG2FLIGHTTXT(on=1/off=0)") == 0 )
    {
      sscanf(ZEIGER, "%d", &input->useJIG2FLIGHT);
    }
    if ( strcmp( ZEILE_READ, "scalefactortostructmodel") == 0 )
    {
      sscanf(ZEIGER, "%lf", &input->scalefactStruct);
    }
    if ( strcmp( ZEILE_READ, "structmodeldihedraldistribution") == 0 ) 
    {    
      char *zeig, *zeig2;
      struct tstring { char teilstring[50]; };
      struct tstring *pointstring;
      int posClose, posPoint, raus;
      raus=0;
      zeig=NULL;
      do 
      {
        zeig = strchr(ZEIGER, (int)'{');
        if (zeig == NULL)
        {
          lauf=0;
          for (lauf=0;lauf <= MaxZeilenLaenge;lauf++)
          {
            ZEILE_READ[lauf] = '\0';
          }
          fgets (ZEILE_READ,MaxZeilenLaenge,fpara);
          ZEIGER=&ZEILE_READ[0];
        }
      }
      while ((zeig == NULL) && (!feof(fpara)));
      if (feof(fpara)) 
      {
        errmessage(17);
      }
      posClose = strcspn( zeig, "}" );
      posPoint = strcspn( zeig, "." );
      //printf("posClose:%d ; posPoint:%d\n",posClose,posPoint);
      if (posClose == posPoint)
      {
        do
        {
          fgets (ZEILE_READ,MaxZeilenLaenge,fpara);
          posClose = strcspn( ZEILE_READ, "}" );
          posPoint = strcspn( ZEILE_READ, "." );
        }while ((posClose!=posPoint) && (!feof(fpara))); 
      }
      if ((posClose < posPoint) )
      {
        printf ("NoDataInput\n");
        raus=1;
      } else if  (feof(fpara)) {
        errmessage(17);
      } else if ((posClose > posPoint))
      {
        zeig=(zeig+1);
        do
        {
          zeig2 = strtok((zeig), "_ \n\t");
          do
          {
            posClose = strcspn( zeig2, "}" );
            posPoint = strcspn( zeig2, "." );
            if ((*(zeig2)==(int)'}') || (posClose < posPoint))
            {
              raus=1;
            }else 
            {
              posClose = strcspn( zeig2, "}" );
              posPoint = strcspn( zeig2, "\0" );
              if (posClose < posPoint)
              {
                 raus=1;
                 zeig2 = strtok((zeig2), "}");
              }
              input->NumStructDih++;
              if (input->NumStructDih==1)
              {
                pointstring=(struct tstring *)malloc(sizeof(struct tstring));
                strcpy(pointstring[0].teilstring, zeig2);              
              }else {
                pointstring=(struct tstring *)realloc(pointstring, input->NumStructDih*sizeof(struct tstring));
                strcpy(pointstring[input->NumStructDih-1].teilstring, zeig2);                 
              }                             
            }                               
            if (raus==0){                         
              zeig2 = strtok(NULL, "_ \n\t");                
            }                               
          }                                 
          while ((zeig2!= NULL) && (raus==0));               
          lauf=0;                              
          for (lauf=0;lauf <= MaxZeilenLaenge;lauf++)            
          {                                 
            ZEILE_READ[lauf] = '\0';                   
          } 
          fgets (ZEILE_READ,MaxZeilenLaenge,fpara);        
          zeig=&ZEILE_READ[0];                  
        }
        while ((raus!=1) && (!feof(fpara)));
      }
      
      if (input->NumStructDih >0)
      {
        input->StrDihDist=(struct Struct_Dihedral_Dist *)malloc(input->NumStructDih*sizeof(struct Struct_Dihedral_Dist));
      }
      for (lauf=0; lauf < input->NumStructDih; lauf++)
      {      
        zeig=pointstring[lauf].teilstring;
        zeig2 = strtok(zeig, ",");
        sscanf(zeig2,"%lf", &input->StrDihDist[lauf].etaPos);
        zeig2 = strtok(NULL, ",\0 ");
        sscanf(zeig2,"%lf", &input->StrDihDist[lauf].etaV);
      }
      /*for (lauf=0; lauf < input->NumStructDih; lauf++)
      {      
        printf("etatPos: %lf    etaV: %lf\n",input->StrDihDist[lauf].etaPos, input->StrDihDist[lauf].etaV);
      }*/    
    }    
    //Kreisflug-Optionen
    if ( strcmp( ZEILE_READ, "circlingmodus(on=1/off=0)") == 0 )
    {
      sscanf(ZEIGER, "%d", &input->KreisflugModus);
    }
    if ( strcmp( ZEILE_READ, "bankangle[deg]") == 0 )
    {
      sscanf(ZEIGER, "%lf", &input->HaengeWinkel);
    }        
    if ( strcmp( ZEILE_READ, "uplift[m/s]") == 0 )
    {
      sscanf(ZEIGER, "%lf", &input->Steigen);
    }        
    //XML-Kopplung
    if ( strcmp( ZEILE_READ, "XML-POLINT-outputactive(on=1/off=0)") == 0 )
    {
      sscanf(ZEIGER, "%d", &input->XMLKOP);
    }    
    if ( strcmp( ZEILE_READ, "startPolint(on=1/off=0)") == 0 )
    {
      sscanf(ZEIGER, "%d", &input->PolintStart);
    }    
    if ( strcmp( ZEILE_READ, "useinternalpolarinterpolation(on=1/off=0)") == 0 )
    {
      sscanf(ZEIGER, "%d", &input->InterneInterpolation);
    }        
    if ( strcmp( ZEILE_READ, "usemassorspeed(0=speed/1=mass)") == 0 )
    {
      sscanf(ZEIGER, "%d", &input->MASSE_SPEED);
    }        
    if ( strcmp( ZEILE_READ, "Polintexecutebal") == 0 ) 
    {
      sscanf(ZEIGER, "%s", &input->POLINTEXE[0]);
    }  
    if ( strcmp( ZEILE_READ, "directorytosectionpolars") == 0 ) 
    {
      sscanf(ZEIGER, "%s", &input->PROFIPFAD[0]);
    } 
    if ( strcmp( ZEILE_READ, "additionalCD") == 0 ) 
    {
      int pos,pos2;
      pos=strcspn(ZEIGER, "AaBbCcDdEeFfGgHhIiJjKkLlMmNnOoPpQqRrSsTtUuVvWwXxYyZz");
      pos2=strcspn(ZEIGER, "0123456789");
      if (pos<pos2)
      {
        sscanf(ZEIGER, "%s", &input->BASIS_CD_POL_File[0]);
        input->BASIS_CD_POL=-1;
      }else{ 
        sscanf(ZEIGER, "%lf", &input->BASIS_CD_POL);
      }      
    }  
    if ( strcmp( ZEILE_READ, "referencemass[kg]") == 0 ) 
    {
      sscanf(ZEIGER, "%lf", &input->POLINT_MASSE);
    }          
    if ( strcmp( ZEILE_READ, "useCm(on=1/off=0)") == 0 )          
    {                                                             
      sscanf(ZEIGER, "%d", &input->useCm);                        
    }                                                             
    if ( strcmp( ZEILE_READ, "useCdpp(on=1/off=0)") == 0 ) 
    {
      sscanf(ZEIGER, "%d", &input->useCdpp);
    }      
    if ( strcmp( ZEILE_READ, "useCdpw(on=1/off=0)") == 0 ) 
    {
      sscanf(ZEIGER, "%d", &input->useCdpw);
    }      
    if ( strcmp( ZEILE_READ, "useCdpf(on=1/off=0)") == 0 ) 
    {
      sscanf(ZEIGER, "%d", &input->useCdpf);
    }         
    if ( strcmp( ZEILE_READ, "polarboundarybehavior(0-5)") == 0 ) 
    {
      sscanf(ZEIGER, "%d", &input->ClMinMaxHandling);
    }      
    if ( strcmp( ZEILE_READ, "penealtyfactor(default=1.1)") == 0 ) 
    {
      sscanf(ZEIGER, "%lf", &input->PENALTYfact);
    }         
    for (zahl=0;zahl < MaxZeilenLaenge;zahl++)
    {
      ZEILE_CLEAN[zahl] = '\0';
      ZEILE_READ[zahl] = '\0';
    }
    fgets (ZEILE_READ,MaxZeilenLaenge,fpara);
  }
  while (!feof(fpara)); 
  fclose(fpara); 
  //printf("Hallo!!\n\n");
  if (input->OldInput == 0)
  {
    int lauf, lauf2, lauf3, NewNumberOfWings;
    NewNumberOfWings=input->NbWings;
    for (lauf=0; lauf <= (input->NbWings-1)  ;lauf++ )
    {
      input->Fluegel[lauf].Schnitte = (struct EinSchnitt *)malloc(input->Fluegel[lauf].anzahlSchnitte*sizeof(struct EinSchnitt));
      for (lauf2=0; lauf2<input->Fluegel[lauf].anzahlSchnitte; lauf2++)
      {
         char dummy[MaxZeilenLaenge] = {"\0"};
         char dummy2[20] = {"\0"};

         sscanf(input->Fluegel[lauf].Surfix, "%s", &dummy[0]);
         //n=strlen(dummy);
         sprintf(dummy2,"%d", lauf2+1);
         strcat(dummy,dummy2);
         strcat(dummy,input->Fluegel[lauf].Prefix);
         input->Fluegel[lauf].Schnitte[lauf2] = read_Schnitt(dummy, input->scalfac, input->readRelDick,input->XMLKOP); 
         //printf("X:%lf  Y:%lf  Z: %lf  T:%lf\n", input->Fluegel[lauf].Schnitte[lauf2].posx, input->Fluegel[lauf].Schnitte[lauf2].posy, input->Fluegel[lauf].Schnitte[lauf2].posz, input->Fluegel[lauf].Schnitte[lauf2].tiefe);
      }
      //An dieser Stelle m�ssen nicht nur die Schnitte gespiegelt werden sondern auch die Kreuzungen und Klappen 
      //Wenn keine Symmettrierandbedingung f�r den Fl�gel definiert ist, muss der Fl�gel nicht in dem Fl�gel gespiegelt weden, sondern ein neuer Fl�gel gespiegelter Fl�gel eingef�gt werden  
      if (input->Fluegel[lauf].Spiegeln==1)
      {
        int newSize;
        printf("\n");
        /*Begin Debug
        for (lauf2=0; lauf2<input->Fluegel[lauf].anzahlSchnitte;lauf2++)
        {
          printf("Fluegel: %d \tSchnitt: %d \tPosY: %lf\n",lauf+1,lauf2+1,input->Fluegel[lauf].Schnitte[lauf2].posy);
        }  End Debug*/
        if (input->Fluegel[lauf].SymExt==1)
        {
          newSize=input->Fluegel[lauf].anzahlSchnitte*2+1;          
        }else{
          newSize=input->Fluegel[lauf].anzahlSchnitte*2-1;        
        }
        input->Fluegel[lauf].Schnitte =(struct EinSchnitt *)realloc(input->Fluegel[lauf].Schnitte, newSize*sizeof(struct EinSchnitt));
        for (lauf2=0; lauf2<input->Fluegel[lauf].anzahlSchnitte; lauf2++)
        {
          input->Fluegel[lauf].Schnitte[newSize-1-lauf2] = copy_Schnitt(&input->Fluegel[lauf].Schnitte[input->Fluegel[lauf].anzahlSchnitte-1-lauf2]);
          /*if ((input->Fluegel[lauf].SymExt==0) &&(lauf2==(input->Fluegel[lauf].anzahlSchnitte-2)))
          {
            lauf2++;
          }*/
        }
        if (input->Fluegel[lauf].SymExt==1)
        {
          input->Fluegel[lauf].Schnitte[input->Fluegel[lauf].anzahlSchnitte] = copy_Schnitt(&input->Fluegel[lauf].Schnitte[0]);
          input->Fluegel[lauf].Schnitte[input->Fluegel[lauf].anzahlSchnitte].posy=0.0;
        }
        for (lauf2=0; lauf2<input->Fluegel[lauf].anzahlSchnitte; lauf2++)
        {
          input->Fluegel[lauf].Schnitte[lauf2] = copy_Schnitt(&input->Fluegel[lauf].Schnitte[newSize-1-lauf2]);
          if (lauf2 > 0)
          {
            input->Fluegel[lauf].Schnitte[lauf2].AnzahlPan=input->Fluegel[lauf].Schnitte[newSize-1-(lauf2-1)].AnzahlPan;
          }
          input->Fluegel[lauf].Schnitte[lauf2].posy=input->Fluegel[lauf].Schnitte[lauf2].posy*(-1);
        }               
        if (input->AnzKreuz>0)
        {
          for (lauf2=0;lauf2<input->AnzKreuz;lauf2++)
          {
            if (input->Kreuz[lauf2].KreuzFl1==(lauf+1))
            {
              input->Kreuz[lauf2].KreuzPosFl1=input->Kreuz[lauf2].KreuzPosFl1+newSize-input->Fluegel[lauf].anzahlSchnitte;
            }
            if (input->Kreuz[lauf2].KreuzFl2==(lauf+1))
            {
              input->Kreuz[lauf2].KreuzPosFl2=input->Kreuz[lauf2].KreuzPosFl2+newSize-input->Fluegel[lauf].anzahlSchnitte;
            }

          }          
        }
        input->Fluegel[lauf].anzahlSchnitte=newSize; 
        if (input->AnzahlKlappen > 0)
        {
          for (lauf2=0;lauf2<input->AnzahlKlappen;lauf2++)
          {
            if (input->Klappen[lauf2].Fluegel==lauf+1)
            {
              input->Klappen[lauf2].Spiegeln=1;
              if (input->Fluegel[lauf].SymExt==1)
              {
                input->Klappen[lauf2].Schnitt1=input->Klappen[lauf2].Schnitt1+(int)(floor(newSize/2))+1;
                input->Klappen[lauf2].Schnitt2=input->Klappen[lauf2].Schnitt2+(int)(floor(newSize/2))+1;
              }else{
                input->Klappen[lauf2].Schnitt1=input->Klappen[lauf2].Schnitt1+(int)(floor(newSize/2));
                input->Klappen[lauf2].Schnitt2=input->Klappen[lauf2].Schnitt2+(int)(floor(newSize/2));
              }
            }
          }
        }
        if (input->NumberOfParametSetsinput > 0)
        {
          for (lauf2=0;lauf2<input->NumberOfParametSetsinput; lauf2++)
          {
            if (input->ParaVar[lauf2].parametric == 1)
            {
              for (lauf3=0;lauf3 < input->ParaVar[lauf2].anzVar; lauf3++)
              {
                if (input->ParaVar[lauf2].VarParameter[lauf3].Fluegel==(lauf+1))
                {
                  if (input->Fluegel[lauf].SymExt==1)
                  {
                    input->ParaVar[lauf2].VarParameter[lauf3].Schnitt=input->ParaVar[lauf2].VarParameter[lauf3].Schnitt+(int)(floor(newSize/2))+1;
                  }else{
                    input->ParaVar[lauf2].VarParameter[lauf3].Schnitt=input->ParaVar[lauf2].VarParameter[lauf3].Schnitt+(int)(floor(newSize/2));
                  }
                  input->ParaVar[lauf2].link=-1;
                }                
              }
            }
            
            if (input->ParaVar[lauf2].parametric == 2)
            {
              for (lauf3=0;lauf3 < input->ParaVar[lauf2].anzSpVar; lauf3++)
              {
                if (input->ParaVar[lauf2].SplineVarPunkt[lauf3].Fluegel==(lauf+1))
                {
                  if (input->Fluegel[lauf].SymExt==1)
                  {
                    input->ParaVar[lauf2].SplineVarPunkt[lauf3].StartSchnitt=input->ParaVar[lauf2].SplineVarPunkt[lauf3].StartSchnitt+(int)(floor(newSize/2))+1;
                    input->ParaVar[lauf2].SplineVarPunkt[lauf3].EndSchnitt=input->ParaVar[lauf2].SplineVarPunkt[lauf3].EndSchnitt+(int)(floor(newSize/2))+1;

                  }else{
                    input->ParaVar[lauf2].SplineVarPunkt[lauf3].StartSchnitt=input->ParaVar[lauf2].SplineVarPunkt[lauf3].StartSchnitt+(int)(floor(newSize/2));
                    input->ParaVar[lauf2].SplineVarPunkt[lauf3].EndSchnitt=input->ParaVar[lauf2].SplineVarPunkt[lauf3].EndSchnitt+(int)(floor(newSize/2));
                  }
                  input->ParaVar[lauf2].link=-1;
                }
              }            
            }
            
            //input->ParaVar[lauf2].anzSpVar
            //input->ParaVar[lauf2].anzVar
          }
        }
        input->Fluegel[lauf].SymExt=0;
        input->Fluegel[lauf].SymRandBed=0;

        /*Begin Debug
        printf("\n");
        for (lauf2=0; lauf2<input->Fluegel[lauf].anzahlSchnitte;lauf2++)
        {
          printf("Fluegel: %d \tSchnitt: %d \tPosY: %lf\tProfilName: %s\n",lauf+1,lauf2+1,input->Fluegel[lauf].Schnitte[lauf2].posy,input->Fluegel[lauf].Schnitte[lauf2].ProfilName);
        }
        End Debug*/        
      }
      if (input->Fluegel[lauf].Spiegeln==2)
      {
        int laufS;
        NewNumberOfWings++;
        input->Fluegel = (struct Einfluegel *)realloc(input->Fluegel,NewNumberOfWings*sizeof(struct Einfluegel));
        input->Fluegel[NewNumberOfWings-1].anzahlSchnitte=input->Fluegel[lauf].anzahlSchnitte;
        input->Fluegel[NewNumberOfWings-1].AnzahlTiefe=input->Fluegel[lauf].AnzahlTiefe;
        input->Fluegel[NewNumberOfWings-1].SymExt=input->Fluegel[lauf].SymExt;
        input->Fluegel[NewNumberOfWings-1].SymRandBed=input->Fluegel[lauf].SymRandBed;
        input->Fluegel[NewNumberOfWings-1].AnzahlSpanPanel=input->Fluegel[lauf].AnzahlSpanPanel;
        input->Fluegel[NewNumberOfWings-1].Spiegeln=0;
        input->Fluegel[NewNumberOfWings-1].Schnitte = (struct EinSchnitt *)malloc(input->Fluegel[NewNumberOfWings-1].anzahlSchnitte*sizeof(struct EinSchnitt));
        for (laufS=0; laufS<input->Fluegel[lauf].anzahlSchnitte; laufS++)
        {
          input->Fluegel[NewNumberOfWings-1].Schnitte[laufS]=copy_Schnitt(&input->Fluegel[lauf].Schnitte[laufS]); 
          input->Fluegel[NewNumberOfWings-1].Schnitte[laufS].posy=input->Fluegel[NewNumberOfWings-1].Schnitte[laufS].posy * (-1);
        }
        input->Fluegel[lauf].Spiegeln=NewNumberOfWings;
        if (input->AnzahlKlappen > 0)
        {
          for (lauf2=0;lauf2<input->AnzahlKlappen;lauf2++)
          {
            if (input->Klappen[lauf2].Fluegel==(lauf+1))
            {
              input->Klappen[lauf2].Spiegeln=NewNumberOfWings;
            }
          }
        }
        
        
        if (input->NumberOfParametSetsinput > 0)
        {
          for (lauf2=0;lauf2<input->NumberOfParametSetsinput; lauf2++)
          {
            if (input->ParaVar[lauf2].parametric == 1)
            {
              for (lauf3=0;lauf3 < input->ParaVar[lauf2].anzVar; lauf3++)
              {
                if (input->ParaVar[lauf2].VarParameter[lauf3].Fluegel==(lauf+1))
                {
                  input->ParaVar[lauf2].link=NewNumberOfWings*(-1);
                }                
              }
            }
            
            if (input->ParaVar[lauf2].parametric == 2)
            {
              for (lauf3=0;lauf3 < input->ParaVar[lauf2].anzSpVar; lauf3++)
              {
                if (input->ParaVar[lauf2].SplineVarPunkt[lauf3].Fluegel==(lauf+1))
                {
                  input->ParaVar[lauf2].link=NewNumberOfWings*(-1);
                }
              }            
            }
            
            //input->ParaVar[lauf2].anzSpVar
            //input->ParaVar[lauf2].anzVar
          }
        }
        //Neuen Fluegel erstellen
      
      
      }      
    }
    printf("Alle Schnitte eingelesen\n");
    input->NbWings=NewNumberOfWings;
  }
  //Restliches Scalieren
  if (input->scalfac != 1)
  {
    input->refSpan=input->refSpan*input->scalfac;
    //input->refAerea=input->refAerea*input->scalfac*input->scalfac;
    input->refAerea=input->refAerea*input->scalfac;
    input->refChordlength=input->refChordlength*input->scalfac;
    input->origXYZ[0]=input->origXYZ[0]*input->scalfac;
    input->origXYZ[1]=input->origXYZ[1]*input->scalfac;
    input->origXYZ[2]=input->origXYZ[2]*input->scalfac;
  }
  //Check Klappen-Input
  if (input->AnzahlKlappen > 0)
  {
    FILE *fKout;
    for (lauf=0;lauf<input->AnzahlKlappen;lauf++)
    {
      
      if (input->Klappen[lauf].Spiegeln==1)
      {
        int ZWIschnitt;
        input->AnzahlKlappen= input->AnzahlKlappen +1;
        input->Klappen = (struct EineKlappen *)realloc(input->Klappen,input->AnzahlKlappen*sizeof(struct EineKlappen));
        input->Klappen[input->AnzahlKlappen-1].Fluegel=input->Klappen[lauf].Fluegel;
        ZWIschnitt=input->Fluegel[input->Klappen[lauf].Fluegel-1].anzahlSchnitte-(input->Klappen[lauf].Schnitt2-1);
        input->Klappen[input->AnzahlKlappen-1].Schnitt2=input->Fluegel[input->Klappen[lauf].Fluegel-1].anzahlSchnitte-(input->Klappen[lauf].Schnitt1-1);
        input->Klappen[input->AnzahlKlappen-1].Schnitt1=ZWIschnitt;        
        input->Klappen[input->AnzahlKlappen-1].AnzPanel=input->Klappen[lauf].AnzPanel;
        input->Klappen[input->AnzahlKlappen-1].Spiegeln=0;
        input->Klappen[input->AnzahlKlappen-1].eta1=input->Klappen[lauf].eta2;
        input->Klappen[input->AnzahlKlappen-1].eta2=input->Klappen[lauf].eta1;
        input->Klappen[input->AnzahlKlappen-1].winkel=input->Klappen[lauf].winkel;
        input->Klappen[lauf].Spiegeln=input->AnzahlKlappen-1;
      }else if (input->Klappen[lauf].Spiegeln>1)
      {
        input->AnzahlKlappen= input->AnzahlKlappen +1;
        input->Klappen = (struct EineKlappen *)realloc(input->Klappen,input->AnzahlKlappen*sizeof(struct EineKlappen));
        input->Klappen[input->AnzahlKlappen-1].Fluegel=input->Klappen[lauf].Spiegeln;
        input->Klappen[input->AnzahlKlappen-1].Schnitt1=input->Klappen[lauf].Schnitt1;
        input->Klappen[input->AnzahlKlappen-1].Schnitt2=input->Klappen[lauf].Schnitt2;
        input->Klappen[input->AnzahlKlappen-1].AnzPanel=input->Klappen[lauf].AnzPanel;
        input->Klappen[input->AnzahlKlappen-1].Spiegeln=0;
        input->Klappen[input->AnzahlKlappen-1].eta1=input->Klappen[lauf].eta1;
        input->Klappen[input->AnzahlKlappen-1].eta2=input->Klappen[lauf].eta2;
        input->Klappen[input->AnzahlKlappen-1].winkel=input->Klappen[lauf].winkel;
        input->Klappen[lauf].Spiegeln=input->AnzahlKlappen-1;
      }           
      // Folgende 2 Abfragen k�nnen hoffentlich nach Lifting-Line-Update wegfallen
      /*if (input->Klappen[lauf].eta1 != input->Klappen[lauf].eta2)
      {
        printf("\nZur Zeit m�ssen eta1 und eta2 einer Klappe gleich sein!\n");
        errmessage(1);
      }*/
      /*if (modf(input->Klappen[lauf].eta1/(1.0/input->Fluegel[input->Klappen[lauf].Fluegel-1].AnzahlTiefe),&dummyZahl2)!=0.0)
      {
        printf("\nZur Zeit muss eta der Klappe ein ein- oder vielfaches von 1/AnzahlFl�gelPanel sein!\n");
        errmessage(1);
      } */    
      // Ende des Blocks der nach Lifting-Line update hoffentlich wegf�llt
      
      
      if (input->Klappen[lauf].Fluegel>input->NbWings)
      {
        errmessage(51);
      }else if (input->Klappen[lauf].Schnitt1>=input->Klappen[lauf].Schnitt2)
      {
        printf("Klappe: %d\tFluegel: %d\tKlappenschnitt 1: %d\t Klappenschnitt 2: %d\n",lauf+1,input->Klappen[lauf].Fluegel,input->Klappen[lauf].Schnitt1,input->Klappen[lauf].Schnitt2);
        errmessage(53);
      }else if (input->Fluegel[input->Klappen[lauf].Fluegel-1].anzahlSchnitte<input->Klappen[lauf].Schnitt2)
      {
        printf("Fluegel: %d\tAnzahl Schnitte auf Fluegel: %d\t Klappenschnitt 2: %d\n",input->Klappen[lauf].Fluegel,input->Fluegel[input->Klappen[lauf].Fluegel-1].anzahlSchnitte,input->Klappen[lauf].Schnitt2);
        errmessage(54);     
      }else if (input->Fluegel[input->Klappen[lauf].Fluegel-1].Schnitte[input->Klappen[lauf].Schnitt1-1].KlappeNrI != 0 && input->Fluegel[input->Klappen[lauf].Fluegel-1].Schnitte[input->Klappen[lauf].Schnitt1-1].KlappeNrI==input->Fluegel[input->Klappen[lauf].Fluegel-1].Schnitte[input->Klappen[lauf].Schnitt1-1].KlappeNrA)
      {
	errmessage(52);
      }else{
        input->Fluegel[input->Klappen[lauf].Fluegel-1].Schnitte[input->Klappen[lauf].Schnitt1-1].KlappeNrA=lauf+1;
	input->Fluegel[input->Klappen[lauf].Fluegel-1].Schnitte[input->Klappen[lauf].Schnitt2-1].KlappeNrI=lauf+1;
        input->Fluegel[input->Klappen[lauf].Fluegel-1].Schnitte[input->Klappen[lauf].Schnitt1-1].etaSchnitt=input->Klappen[lauf].eta1;
        input->Fluegel[input->Klappen[lauf].Fluegel-1].Schnitte[input->Klappen[lauf].Schnitt2-1].etaSchnitt=input->Klappen[lauf].eta2;
	if (input->Klappen[lauf].AnzPanel==0)
        {
          input->Klappen[lauf].AnzPanel=roundfloat(((input->Klappen[lauf].eta1+input->Klappen[lauf].eta2)/2)/(1/(float)input->Fluegel[input->Klappen[lauf].Fluegel-1].AnzahlTiefe));
        }
        //printf("%lf \t %d \t%lf\n",(input->Klappen[lauf].eta1+input->Klappen[lauf].eta2)/2 ,input->Klappen[lauf].Fluegel, 1/(float)input->Fluegel[input->Klappen[lauf].Fluegel-1].AnzahlTiefe);
	if (input->Klappen[lauf].AnzPanel < 1)
	{
	  input->Klappen[lauf].AnzPanel=1;
	}
	if (input->Fluegel[input->Klappen[lauf].Fluegel-1].Schnitte[input->Klappen[lauf].Schnitt1-1].anzPanelKlappe<input->Klappen[lauf].AnzPanel)
	{
	  input->Fluegel[input->Klappen[lauf].Fluegel-1].Schnitte[input->Klappen[lauf].Schnitt1-1].anzPanelKlappe=input->Klappen[lauf].AnzPanel;
	}
	if (input->Fluegel[input->Klappen[lauf].Fluegel-1].Schnitte[input->Klappen[lauf].Schnitt2-1].anzPanelKlappe<input->Klappen[lauf].AnzPanel)
	{
          input->Fluegel[input->Klappen[lauf].Fluegel-1].Schnitte[input->Klappen[lauf].Schnitt2-1].anzPanelKlappe=input->Klappen[lauf].AnzPanel;
        }
	if (input->Klappen[lauf].Schnitt1+1 < input->Klappen[lauf].Schnitt2)
	{
	  int lauf2;
	  double eta1X, eta1Y, eta1Z,eta2X, eta2Y, eta2Z, alfa,Leta2XYZ, Abstand;
	  for (lauf2=input->Klappen[lauf].Schnitt1+1; lauf2 < input->Klappen[lauf].Schnitt2; lauf2++)
	  {	    
	    input->Fluegel[input->Klappen[lauf].Fluegel-1].Schnitte[lauf2-1].KlappeNrI=lauf+1;
	    input->Fluegel[input->Klappen[lauf].Fluegel-1].Schnitte[lauf2-1].KlappeNrA=lauf+1;
	    input->Fluegel[input->Klappen[lauf].Fluegel-1].Schnitte[lauf2-1].anzPanelKlappe=input->Klappen[lauf].AnzPanel;
	    eta1X=input->Fluegel[input->Klappen[lauf].Fluegel-1].Schnitte[input->Klappen[lauf].Schnitt1-1].posx+input->Fluegel[input->Klappen[lauf].Fluegel-1].Schnitte[input->Klappen[lauf].Schnitt1-1].tiefe*(1.0-input->Fluegel[input->Klappen[lauf].Fluegel-1].Schnitte[input->Klappen[lauf].Schnitt1-1].etaSchnitt);
	    eta1Y=input->Fluegel[input->Klappen[lauf].Fluegel-1].Schnitte[input->Klappen[lauf].Schnitt1-1].posy;
	    eta1Z=input->Fluegel[input->Klappen[lauf].Fluegel-1].Schnitte[input->Klappen[lauf].Schnitt1-1].posz;
	    eta2X=input->Fluegel[input->Klappen[lauf].Fluegel-1].Schnitte[input->Klappen[lauf].Schnitt2-1].posx+input->Fluegel[input->Klappen[lauf].Fluegel-1].Schnitte[input->Klappen[lauf].Schnitt2-1].tiefe*(1.0-input->Fluegel[input->Klappen[lauf].Fluegel-1].Schnitte[input->Klappen[lauf].Schnitt2-1].etaSchnitt);
	    eta2Y=input->Fluegel[input->Klappen[lauf].Fluegel-1].Schnitte[input->Klappen[lauf].Schnitt2-1].posy;
	    eta2Z=input->Fluegel[input->Klappen[lauf].Fluegel-1].Schnitte[input->Klappen[lauf].Schnitt2-1].posz;
	    eta2X=eta2X-eta1X;
	    eta2Y=eta2Y-eta1Y;
	    eta2Z=eta2Z-eta1Z;
	    Leta2XYZ=sqrt(pow(eta2X,2)+pow(eta2Y,2)+pow(eta2Z,2));
	    //printf("eta1: %lf %lf %lf eta2Vec:  %lf %lf %lf  Leta: %lf\n",eta1X,eta1Y,eta1Z,eta2X,eta2Y,eta2Z,Leta2XYZ);
	    eta2X=eta2X/Leta2XYZ; //eta2X,Y,Z ist jetzt Normalenvector von der Ebene durch den Schnitt!
	    eta2Y=eta2Y/Leta2XYZ;
	    eta2Z=eta2Z/Leta2XYZ;
	    
            //Neu �ber Winkel
            alfa=asin(eta2X);
            eta1Y=input->Fluegel[input->Klappen[lauf].Fluegel-1].Schnitte[input->Klappen[lauf].Schnitt1-1].posy;
            eta2Y=input->Fluegel[input->Klappen[lauf].Fluegel-1].Schnitte[lauf2-1].posy;
            eta1Z=input->Fluegel[input->Klappen[lauf].Fluegel-1].Schnitte[input->Klappen[lauf].Schnitt1-1].posz;
            eta2Z=input->Fluegel[input->Klappen[lauf].Fluegel-1].Schnitte[lauf2-1].posz;
            Abstand=sqrt(pow(eta1Y-eta2Y,2)+pow(eta1Z-eta2Z,2));
            eta2X=eta1X+tan(alfa)*Abstand;
            input->Fluegel[input->Klappen[lauf].Fluegel-1].Schnitte[lauf2-1].etaSchnitt=1.0-((eta2X-input->Fluegel[input->Klappen[lauf].Fluegel-1].Schnitte[lauf2-1].posx)/input->Fluegel[input->Klappen[lauf].Fluegel-1].Schnitte[lauf2-1].tiefe);
            //printf("eta1X:%lf\teta2X:%lf\tposx:%lf\ttiefe:%lf\n",eta1X,eta2X,input->Fluegel[input->Klappen[lauf].Fluegel-1].Schnitte[lauf2-1].posx,input->Fluegel[input->Klappen[lauf].Fluegel-1].Schnitte[lauf2-1].tiefe);
            if ((input->Fluegel[input->Klappen[lauf].Fluegel-1].Schnitte[lauf2-1].etaSchnitt >1) ||(input->Fluegel[input->Klappen[lauf].Fluegel-1].Schnitte[lauf2-1].etaSchnitt <0))
            {
              printf("\nError: A straight hingeline cant be realized for flap %d\n",lauf+1);
              errmessage(1);
            }
            
            /* Alt falsch mit ebene??
            aEbene=eta2X*input->Fluegel[input->Klappen[lauf].Fluegel-1].Schnitte[lauf2-1].posx+eta2Y*input->Fluegel[input->Klappen[lauf].Fluegel-1].Schnitte[lauf2-1].posy+eta2Z*input->Fluegel[input->Klappen[lauf].Fluegel-1].Schnitte[lauf2-1].posz;
	    Abstand=fabs(eta1X*eta2X+eta1Y*eta2Y+eta1Z*eta2Z-aEbene);
	    //printf("Abstand: %lf  posXSchnitt: %lf\n",Abstand,(Abstand*eta2X)+eta1X);
	    input->Fluegel[input->Klappen[lauf].Fluegel-1].Schnitte[lauf2-1].etaSchnitt=1.0-((((Abstand*eta2X)+eta1X)-input->Fluegel[input->Klappen[lauf].Fluegel-1].Schnitte[lauf2-1].posx)/input->Fluegel[input->Klappen[lauf].Fluegel-1].Schnitte[lauf2-1].tiefe);
            */
          }
	}
      }
  
    }
    fKout=fopen("Flap_list.dat","w");
    fprintf(fKout,"FlapNr\twing\tSpanpos1\tSpanpos2\teta1\t\teta2\t\tFlapdeflection\n");
    for (lauf=0;lauf<input->AnzahlKlappen;lauf++)
    {
      fprintf(fKout,"%d\t",lauf+1);
      fprintf(fKout,"%d\t",input->Klappen[lauf].Fluegel);
      fprintf(fKout,"%d\t\t",input->Klappen[lauf].Schnitt1);
      fprintf(fKout,"%d\t\t",input->Klappen[lauf].Schnitt2);
      fprintf(fKout,"%e\t",input->Klappen[lauf].eta1);
      fprintf(fKout,"%e\t",input->Klappen[lauf].eta2);
      fprintf(fKout,"%f\n",input->Klappen[lauf].winkel);
    }
    fclose(fKout);     
    
    //Debug Klappen begin
    if (input->debuglevel==1 || input->debuglevel>99)
    {
      for (lauf=0;lauf<input->AnzahlKlappen;lauf++)
      {
        printf("FlapNr\twing\tSpanpos1\tSpanpos2\teta1\t\teta2\t\tFlapdeflection\n");
        printf("%d\t",lauf+1);
        printf("%d\t",input->Klappen[lauf].Fluegel);
        printf("%d\t\t",input->Klappen[lauf].Schnitt1);
        printf("%d\t\t",input->Klappen[lauf].Schnitt2);
        printf("%e\t",input->Klappen[lauf].eta1);
        printf("%e\t",input->Klappen[lauf].eta2);
        printf("%f\n",input->Klappen[lauf].winkel);
      }
      printf("\n");
      for (lauf=0; lauf<input->NbWings; lauf++)
      {
         int lauf2;
         for (lauf2=0; lauf2<input->Fluegel[lauf].anzahlSchnitte; lauf2++)
         {
           printf("Fluegel: %d \t Schnitt: %d \t KlappeI: %d \t KlappeA: %d\tanzPanel: %d \t etaKlappe: %lf\t Profil-Name: %s\n",lauf+1,lauf2+1,input->Fluegel[lauf].Schnitte[lauf2].KlappeNrI,input->Fluegel[lauf].Schnitte[lauf2].KlappeNrA,input->Fluegel[lauf].Schnitte[lauf2].anzPanelKlappe,input->Fluegel[lauf].Schnitte[lauf2].etaSchnitt,input->Fluegel[lauf].Schnitte[lauf2].ProfilName );
           /*for (lauf3=0;lauf3<input->Fluegel[lauf].Schnitte[lauf2].AnzahlXY;lauf3++)
           {
             printf("xcamb: %lf \tzcamb: %lf\n",input->Fluegel[lauf].Schnitte[lauf2].xcamb[lauf3],input->Fluegel[lauf].Schnitte[lauf2].zcamb[lauf3]);  
           }*/
         }
      }
      printf("\n");
      //exit (EXIT_FAILURE);
    }
    //Debug Klappen end
  }
  //check Mischer
  /*if (input->debuglevel==3 || input->debuglevel>99)
  {
    if ( input->AnzahlMischer > 0 )
    {
      check_Mischer(input);
    }
  }*/
  //Check optimization input
  //printf ("Number of Parameter Sets in Input:%d\n",input->NumberOfParametSetsinput);
  if (input->AnzKreuz>0)
  {
    int newNumbKreuz;
    newNumbKreuz=input->AnzKreuz;
    for (lauf=0;lauf< input->AnzKreuz; lauf++)
    {
      //printf ("\nCheck Kreuzspiegelung: nFl1:%d\tnFl2s:%d\tFl1:%d\tFl2:%d\n",input->Fluegel[input->Kreuz[lauf].KreuzFl1-1].Spiegeln, input->Fluegel[input->Kreuz[lauf].KreuzFl2-1].Spiegeln, input->Kreuz[lauf].KreuzFl1,input->Kreuz[lauf].KreuzFl2);
      if(input->Fluegel[input->Kreuz[lauf].KreuzFl1-1].Spiegeln == 1 && input->Fluegel[input->Kreuz[lauf].KreuzFl2-1].Spiegeln > 1 && input->Kreuz[lauf].KreuzFl2 != input->Kreuz[lauf].KreuzFl1)
      {
        newNumbKreuz++;
        input->Kreuz = (struct Kreuzung *)realloc(input->Kreuz,newNumbKreuz*sizeof(struct Kreuzung));
        input->Kreuz[newNumbKreuz-1].KreuzFl1=input->Kreuz[lauf].KreuzFl1;
        input->Kreuz[newNumbKreuz-1].KreuzFl2=input->Fluegel[input->Kreuz[lauf].KreuzFl2-1].Spiegeln;        
        input->Kreuz[newNumbKreuz-1].KreuzPosFl1=input->Fluegel[input->Kreuz[lauf].KreuzFl1-1].anzahlSchnitte-(input->Kreuz[lauf].KreuzPosFl1-1);
        input->Kreuz[newNumbKreuz-1].KreuzPosFl2=input->Kreuz[lauf].KreuzPosFl2;
        input->Kreuz[newNumbKreuz-1].KreuzTiefFl1=input->Kreuz[lauf].KreuzTiefFl1;
        input->Kreuz[newNumbKreuz-1].KreuzTiefFl2=input->Kreuz[lauf].KreuzTiefFl2;
        input->Kreuz[newNumbKreuz-1].GabBedFl1=input->Kreuz[lauf].GabBedFl1;
        input->Kreuz[newNumbKreuz-1].GabBedFl2=input->Kreuz[lauf].GabBedFl2;
        input->Kreuz[newNumbKreuz-1].KopelArt=input->Kreuz[lauf].KopelArt;
      }      
      if(input->Fluegel[input->Kreuz[lauf].KreuzFl1-1].Spiegeln > 1 && input->Fluegel[input->Kreuz[lauf].KreuzFl2-1].Spiegeln == 1 && input->Kreuz[lauf].KreuzFl2 != input->Kreuz[lauf].KreuzFl1)
      {
        newNumbKreuz++;
        input->Kreuz = (struct Kreuzung *)realloc(input->Kreuz,newNumbKreuz*sizeof(struct Kreuzung));
        input->Kreuz[newNumbKreuz-1].KreuzFl1=input->Fluegel[input->Kreuz[lauf].KreuzFl1-1].Spiegeln;    
        input->Kreuz[newNumbKreuz-1].KreuzFl2=input->Kreuz[lauf].KreuzFl2;    
        input->Kreuz[newNumbKreuz-1].KreuzPosFl1=input->Kreuz[lauf].KreuzPosFl1;
        input->Kreuz[newNumbKreuz-1].KreuzPosFl2=input->Fluegel[input->Kreuz[lauf].KreuzFl2-1].anzahlSchnitte-(input->Kreuz[lauf].KreuzPosFl2-1);
        input->Kreuz[newNumbKreuz-1].KreuzTiefFl1=input->Kreuz[lauf].KreuzTiefFl1;
        input->Kreuz[newNumbKreuz-1].KreuzTiefFl2=input->Kreuz[lauf].KreuzTiefFl2;
        input->Kreuz[newNumbKreuz-1].GabBedFl1=input->Kreuz[lauf].GabBedFl1;
        input->Kreuz[newNumbKreuz-1].GabBedFl2=input->Kreuz[lauf].GabBedFl2;
        input->Kreuz[newNumbKreuz-1].KopelArt=input->Kreuz[lauf].KopelArt;
      } 
      /*if(input->Fluegel[input->Kreuz[lauf].KreuzFl1-1].Spiegeln == 1 && input->Fluegel[input->Kreuz[lauf].KreuzFl2-1].Spiegeln == 0) 
      {
        input->Kreuz[lauf].KreuzPosFl1=input->Fluegel[input->Kreuz[lauf].KreuzFl1-1].anzahlSchnitte-(input->Kreuz[lauf].KreuzPosFl1-1);
      }
      if(input->Fluegel[input->Kreuz[lauf].KreuzFl1-1].Spiegeln == 0 && input->Fluegel[input->Kreuz[lauf].KreuzFl2-1].Spiegeln == 1) 
      {
        input->Kreuz[lauf].KreuzPosFl2=input->Fluegel[input->Kreuz[lauf].KreuzFl2-1].anzahlSchnitte-(input->Kreuz[lauf].KreuzPosFl2-1);
      } Doppelt gemoppelt h�lt in diesem Fall nicht besser! Wird bereits oben im Fl�gelspiegeln umgesetzt!*/
      if(input->Fluegel[input->Kreuz[lauf].KreuzFl1-1].Spiegeln > 1 && input->Fluegel[input->Kreuz[lauf].KreuzFl2-1].Spiegeln == 0 && input->Kreuz[lauf].KreuzFl2 != input->Kreuz[lauf].KreuzFl1)
      {
        errmessage(59);
      }
      if(input->Fluegel[input->Kreuz[lauf].KreuzFl1-1].Spiegeln == 0 && input->Fluegel[input->Kreuz[lauf].KreuzFl2-1].Spiegeln > 1 && input->Kreuz[lauf].KreuzFl2 != input->Kreuz[lauf].KreuzFl1)
      {
        errmessage(59);      
      }
      if(input->Fluegel[input->Kreuz[lauf].KreuzFl1-1].Spiegeln > 1 && input->Fluegel[input->Kreuz[lauf].KreuzFl2-1].Spiegeln > 1 && input->Kreuz[lauf].KreuzFl2 != input->Kreuz[lauf].KreuzFl1) 
      {
        newNumbKreuz++;
        input->Kreuz = (struct Kreuzung *)realloc(input->Kreuz,newNumbKreuz*sizeof(struct Kreuzung));
        input->Kreuz[newNumbKreuz-1].KreuzFl1=input->Fluegel[input->Kreuz[lauf].KreuzFl1-1].Spiegeln;    
        input->Kreuz[newNumbKreuz-1].KreuzFl2=input->Fluegel[input->Kreuz[lauf].KreuzFl2-1].Spiegeln;  
        input->Kreuz[newNumbKreuz-1].KreuzPosFl1=input->Kreuz[lauf].KreuzPosFl1;
        input->Kreuz[newNumbKreuz-1].KreuzPosFl2=input->Kreuz[lauf].KreuzPosFl2;
        input->Kreuz[newNumbKreuz-1].KreuzTiefFl1=input->Kreuz[lauf].KreuzTiefFl1;
        input->Kreuz[newNumbKreuz-1].KreuzTiefFl2=input->Kreuz[lauf].KreuzTiefFl2;
        input->Kreuz[newNumbKreuz-1].GabBedFl1=input->Kreuz[lauf].GabBedFl1;
        input->Kreuz[newNumbKreuz-1].GabBedFl2=input->Kreuz[lauf].GabBedFl2;
        input->Kreuz[newNumbKreuz-1].KopelArt=input->Kreuz[lauf].KopelArt;
      }
      if(input->Fluegel[input->Kreuz[lauf].KreuzFl1-1].Spiegeln == 1 && input->Fluegel[input->Kreuz[lauf].KreuzFl2-1].Spiegeln == 1 && input->Kreuz[lauf].KreuzFl2 != input->Kreuz[lauf].KreuzFl1) 
      {
        newNumbKreuz++;      
        input->Kreuz = (struct Kreuzung *)realloc(input->Kreuz,newNumbKreuz*sizeof(struct Kreuzung));
        input->Kreuz[newNumbKreuz-1].KreuzFl1=input->Kreuz[lauf].KreuzFl1;  
        input->Kreuz[newNumbKreuz-1].KreuzFl2=input->Kreuz[lauf].KreuzFl2;    
        input->Kreuz[newNumbKreuz-1].KreuzPosFl1=input->Fluegel[input->Kreuz[lauf].KreuzFl1-1].anzahlSchnitte-(input->Kreuz[lauf].KreuzPosFl1-1);
        input->Kreuz[newNumbKreuz-1].KreuzPosFl2=input->Fluegel[input->Kreuz[lauf].KreuzFl2-1].anzahlSchnitte-(input->Kreuz[lauf].KreuzPosFl2-1);
        input->Kreuz[newNumbKreuz-1].KreuzTiefFl1=input->Kreuz[lauf].KreuzTiefFl1;
        input->Kreuz[newNumbKreuz-1].KreuzTiefFl2=input->Kreuz[lauf].KreuzTiefFl2;
        input->Kreuz[newNumbKreuz-1].GabBedFl1=input->Kreuz[lauf].GabBedFl1;
        input->Kreuz[newNumbKreuz-1].GabBedFl2=input->Kreuz[lauf].GabBedFl2;
        input->Kreuz[newNumbKreuz-1].KopelArt=input->Kreuz[lauf].KopelArt;
      }      
      if(input->Fluegel[input->Kreuz[lauf].KreuzFl1-1].Spiegeln > 0 && input->Kreuz[lauf].KreuzFl2 == input->Kreuz[lauf].KreuzFl1  && input->Fluegel[input->Kreuz[lauf].KreuzFl1-1].Spiegeln == input->Fluegel[input->Kreuz[lauf].KreuzFl2-1].Spiegeln ) 
      {
        newNumbKreuz++;     
        input->Kreuz = (struct Kreuzung *)realloc(input->Kreuz,newNumbKreuz*sizeof(struct Kreuzung));
        if(input->Fluegel[input->Kreuz[lauf].KreuzFl1-1].Spiegeln ==1)
        {  
          input->Kreuz[newNumbKreuz-1].KreuzFl1=input->Kreuz[lauf].KreuzFl1;  
          input->Kreuz[newNumbKreuz-1].KreuzFl2=input->Kreuz[lauf].KreuzFl2;            
          input->Kreuz[newNumbKreuz-1].KreuzPosFl1=input->Fluegel[input->Kreuz[lauf].KreuzFl1-1].anzahlSchnitte-(input->Kreuz[lauf].KreuzPosFl1-1);
          input->Kreuz[newNumbKreuz-1].KreuzPosFl2=input->Fluegel[input->Kreuz[lauf].KreuzFl2-1].anzahlSchnitte-(input->Kreuz[lauf].KreuzPosFl2-1);
        }else{
          input->Kreuz[newNumbKreuz-1].KreuzFl1=input->Fluegel[input->Kreuz[lauf].KreuzFl1-1].Spiegeln;    
          input->Kreuz[newNumbKreuz-1].KreuzFl2=input->Fluegel[input->Kreuz[lauf].KreuzFl2-1].Spiegeln;  
          input->Kreuz[newNumbKreuz-1].KreuzPosFl1=input->Kreuz[lauf].KreuzPosFl1;
          input->Kreuz[newNumbKreuz-1].KreuzPosFl2=input->Kreuz[lauf].KreuzPosFl2;
        }
        input->Kreuz[newNumbKreuz-1].KreuzTiefFl1=input->Kreuz[lauf].KreuzTiefFl1;
        input->Kreuz[newNumbKreuz-1].KreuzTiefFl2=input->Kreuz[lauf].KreuzTiefFl2;
        input->Kreuz[newNumbKreuz-1].GabBedFl1=input->Kreuz[lauf].GabBedFl1;
        input->Kreuz[newNumbKreuz-1].GabBedFl2=input->Kreuz[lauf].GabBedFl2;
        input->Kreuz[newNumbKreuz-1].KopelArt=input->Kreuz[lauf].KopelArt;      
      }



      /*if(input->Fluegel[input->Kreuz[lauf].KreuzFl2-1].Spiegeln > 1)
      if(input->Fluegel[input->Kreuz[lauf].KreuzFl1-1].Spiegeln == 1)  
      if(input->Fluegel[input->Kreuz[lauf].KreuzFl2-1].Spiegeln == 1)*/
    }
    input->AnzKreuz=newNumbKreuz;
  }
  if (input->Austrimmen > 0)
  {
    checkTrimm(input);
  }
  
  if (input->NumberOfParametSetsinput>=1)
  {
    struct PARAvariation *ParaVarNeu;
    int NewNumberOfParametSetsinput, lauf2,ZwiZahl;
    NewNumberOfParametSetsinput=0;
    for (lauf=0; lauf<input->NumberOfParametSetsinput; lauf++)
    {
      if (input->ParaVar[lauf].parametric>0)
      {
        if ((input->ParaVar[lauf].Paraset4eachCLA==1)&&(input->runOneByOne==0)&&((input->numberAlpha+input->numberCL)>1))
        {
          input->runOneByOne=1;
          printf("Warning: run one by one (on=1/off=0) has been switched on to meet requierments of parametrization!\n");
        }
        NewNumberOfParametSetsinput++;
        if (NewNumberOfParametSetsinput==1)
        {
          ParaVarNeu=(struct PARAvariation *)malloc(NewNumberOfParametSetsinput*sizeof(struct PARAvariation));              
        }else{
          ParaVarNeu=(struct PARAvariation *)realloc(ParaVarNeu, NewNumberOfParametSetsinput*sizeof(struct PARAvariation));
        }
        ParaVarNeu[NewNumberOfParametSetsinput-1]=input->ParaVar[lauf];
      }        
    }
    input->NumberOfParametSetsinput=NewNumberOfParametSetsinput;
    input->ParaVar=ParaVarNeu;
    for (lauf=0; lauf<input->NumberOfParametSetsinput; lauf++)
    {
      if (NewNumberOfParametSetsinput==1)
      {
        for(lauf2=0; lauf2<input->ParaVar[lauf].anzVar; lauf2++)
        {
          if (lauf2>0)
          {
            if (ZwiZahl!=input->ParaVar[lauf].VarParameter[lauf2].Fluegel)
            {
              errmessage(60);
            } 
          }
          ZwiZahl=input->ParaVar[lauf].VarParameter[lauf2].Fluegel;
        }
      }
      if (NewNumberOfParametSetsinput==2)
      {
        for(lauf2=0; lauf2<input->ParaVar[lauf].anzSpVar; lauf2++)
        {
          if (lauf2>0)
          {
            if (ZwiZahl!=input->ParaVar[lauf].SplineVarPunkt[lauf2].Fluegel)
            {
              errmessage(60);
            } 
          }
          ZwiZahl=input->ParaVar[lauf].SplineVarPunkt[lauf2].Fluegel;
        }
      }           
    }
  }    
  //printf ("Final Number of Parameter Sets in Input:%d\n",input->NumberOfParametSetsinput);
  if (input->NumberOfParametSetsinput>=2)
  {
    int lauf2;
    for (lauf=0; lauf<(input->NumberOfParametSetsinput-1); lauf++)
    {
      for (lauf2=(lauf+1); lauf2<input->NumberOfParametSetsinput; lauf2++)
      {
        if (!strcmp(input->ParaVar[lauf].PARAINP,input->ParaVar[lauf2].PARAINP))
        {
          errmessage(56);
        }
      }
    }
  }
  zahl=input->NumberOfParametSetsinput;
  for (lauf=0;lauf<input->NumberOfParametSetsinput; lauf++)
  {
    if (input->ParaVar[lauf].link<=-1)
    {
      int lauf2;
      zahl++;
      input->ParaVar=(struct PARAvariation *)realloc(input->ParaVar, zahl*sizeof(struct PARAvariation));      
      input->ParaVar[zahl-1].parametric=input->ParaVar[lauf].parametric;      
      input->ParaVar[zahl-1].NrVariations=input->ParaVar[lauf].NrVariations;    
      input->ParaVar[zahl-1].anzVar=input->ParaVar[lauf].anzVar;
      input->ParaVar[zahl-1].Paraset4eachCLA=input->ParaVar[lauf].Paraset4eachCLA;
      input->ParaVar[zahl-1].AbsOrDiffer=input->ParaVar[lauf].AbsOrDiffer;
      input->ParaVar[zahl-1].anzSpVar=input->ParaVar[lauf].anzSpVar;
      input->ParaVar[zahl-1].AnzahlParameter=input->ParaVar[lauf].AnzahlParameter;      
      if (input->ParaVar[lauf].parametric==1)
      {
        input->ParaVar[zahl-1].VarParameter=(struct VariationsPunkt *)malloc(input->ParaVar[lauf].anzVar*sizeof(struct VariationsPunkt));
        for(lauf2=0;lauf2<input->ParaVar[lauf].anzVar; lauf2++)
        {
          input->ParaVar[zahl-1].VarParameter[lauf2].Fluegel=input->ParaVar[lauf].VarParameter[lauf2].Fluegel;
          if (input->ParaVar[lauf].link==-1) 
          {
            input->ParaVar[zahl-1].VarParameter[lauf2].Schnitt=input->Fluegel[input->ParaVar[lauf].VarParameter[lauf2].Fluegel-1].anzahlSchnitte-(input->ParaVar[lauf].VarParameter[lauf2].Schnitt-1);
            input->ParaVar[zahl-1].VarParameter[lauf2].Fluegel=input->ParaVar[lauf].VarParameter[lauf2].Fluegel;
          }else{            
            input->ParaVar[zahl-1].VarParameter[lauf2].Schnitt=input->ParaVar[lauf].VarParameter[lauf2].Schnitt;
            input->ParaVar[zahl-1].VarParameter[lauf2].Fluegel=input->ParaVar[lauf].link*(-1);          
          }
           input->ParaVar[zahl-1].VarParameter[lauf2].twist=input->ParaVar[lauf].VarParameter[lauf2].twist;
          input->ParaVar[zahl-1].VarParameter[lauf2].X=input->ParaVar[lauf].VarParameter[lauf2].X;
          input->ParaVar[zahl-1].VarParameter[lauf2].Y=input->ParaVar[lauf].VarParameter[lauf2].Y;
          input->ParaVar[zahl-1].VarParameter[lauf2].Z=input->ParaVar[lauf].VarParameter[lauf2].Z;
          input->ParaVar[zahl-1].VarParameter[lauf2].tiefe=input->ParaVar[lauf].VarParameter[lauf2].tiefe;
          input->ParaVar[zahl-1].VarParameter[lauf2].Vloc=input->ParaVar[lauf].VarParameter[lauf2].Vloc;
        }
      }
      if (input->ParaVar[lauf].parametric==2)
      {
        for(lauf2=0;lauf2<input->ParaVar[lauf].anzSpVar; lauf2++)
        {
          input->ParaVar[zahl-1].SplineVarPunkt=(struct SplineVariationsPunkt *)malloc(input->ParaVar[lauf].anzSpVar*sizeof(struct SplineVariationsPunkt));
          input->ParaVar[zahl-1].SplineVarPunkt[lauf2].AnzahlStuetzpunkte=input->ParaVar[lauf].SplineVarPunkt[lauf2].AnzahlStuetzpunkte;          
          if (input->ParaVar[lauf].link==-1)
          {
            input->ParaVar[zahl-1].SplineVarPunkt[lauf2].StartSchnitt=input->Fluegel[input->ParaVar[lauf].SplineVarPunkt[lauf2].Fluegel-1].anzahlSchnitte-(input->ParaVar[lauf].SplineVarPunkt[lauf2].StartSchnitt-1);
            input->ParaVar[zahl-1].SplineVarPunkt[lauf2].EndSchnitt=input->Fluegel[input->ParaVar[lauf].SplineVarPunkt[lauf2].Fluegel-1].anzahlSchnitte-(input->ParaVar[lauf].SplineVarPunkt[lauf2].EndSchnitt-1);           
            input->ParaVar[zahl-1].SplineVarPunkt[lauf2].Fluegel=input->ParaVar[lauf].SplineVarPunkt[lauf2].Fluegel;
          }else{ 
            input->ParaVar[zahl-1].SplineVarPunkt[lauf2].StartSchnitt=input->ParaVar[lauf].SplineVarPunkt[lauf2].StartSchnitt;
            input->ParaVar[zahl-1].SplineVarPunkt[lauf2].EndSchnitt=input->ParaVar[lauf].SplineVarPunkt[lauf2].EndSchnitt; 
            input->ParaVar[zahl-1].SplineVarPunkt[lauf2].Fluegel=input->ParaVar[lauf].link*(-1);          
          }
          input->ParaVar[zahl-1].SplineVarPunkt[lauf2].Variable.twist=input->ParaVar[lauf].SplineVarPunkt[lauf2].Variable.twist;
          input->ParaVar[zahl-1].SplineVarPunkt[lauf2].Variable.X=input->ParaVar[lauf].SplineVarPunkt[lauf2].Variable.X;
          input->ParaVar[zahl-1].SplineVarPunkt[lauf2].Variable.Y=input->ParaVar[lauf].SplineVarPunkt[lauf2].Variable.Y;
          input->ParaVar[zahl-1].SplineVarPunkt[lauf2].Variable.Z=input->ParaVar[lauf].SplineVarPunkt[lauf2].Variable.Z;
          input->ParaVar[zahl-1].SplineVarPunkt[lauf2].Variable.tiefe=input->ParaVar[lauf].SplineVarPunkt[lauf2].Variable.tiefe;
          input->ParaVar[zahl-1].SplineVarPunkt[lauf2].Variable.Vloc=input->ParaVar[lauf].SplineVarPunkt[lauf2].Variable.Vloc;
        }      
      }
      input->ParaVar[zahl-1].link=lauf+1;
      input->ParaVar[lauf].link=(zahl-1)*(-1);

    }
  }
  input->NumberOfParametSetsinput=zahl;
  //printf("Numb Para %d\n",input->NumberOfParametSetsinput);
  if ((input->BASIS_CD_POL==-1) && (input->PolintStart==1 || input->InterneInterpolation==1))
  {
    int check;
    FILE *fCDInp;
    fCDInp = fopen(input->BASIS_CD_POL_File,"r");
    if (fCDInp ==NULL)
    {
      errmessage(61);      
    }
    check=fscanf(fCDInp," %lf ",&input->BASIS_CD_POL); 
    if (check == EOF)
    {
      errmessage(62);
    }
  }
  //Setzten der des Liftin-Line Verzeichniss Namen 
  if (strcmp(input->LiLiVersion,"V2p3bCirc")==0)
  {
    strcpy(input->LiLiVersionName,"V2.3");
  }else if (strcmp(input->LiLiVersion,"V2p3b")==0){
    strcpy(input->LiLiVersionName,"V2.3beta");
  }
  
  //Einige generelle Checks
  if (input->KreisflugModus==1 && strcmp(input->LiLiVersion,"V2p3bCirc")!=0)
  {
    errmessage(64);
  }
  if (input->LiliStart>0 && input->FwStart>0)
  {
    errmessage(67);
  }
  if (input->LiliStart<=0 && input->FwStart<=0)
  {
    errmessage(68);
  }
  
  
  
/*  if (input->MASSE_SPEED!=0)
  {
    input->runOneByOne =1;
  }*/
  /*Debug Check*/
  /*for (zahl =0; zahl<input->anzSpVar; zahl++)
  {
    printf ("VariableCheck: Fluegel:%d   StartSchnitt: %d   EndSchnitt: %d \n",input->SplineVarPunkt[zahl].Fluegel, input->SplineVarPunkt[zahl].StartSchnitt, input->SplineVarPunkt[zahl].EndSchnitt);
  }  */
  //printf ("VariableCheck: %lf\n", input->Alfa_Rot);
}

  
  
