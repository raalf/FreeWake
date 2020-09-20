/******************POALR_INTER.c******************                  
*                                                *      
*      Vorentwurfs-Generator "VoGen"             *          
*      Version 1.17                              *
*      Date of last modification 04.01.2019      *
*                                                *          
*      Writen by: Jan Himisch                    *      
*             jan.himisch@dlr.de                 *      
*                                                *      
*                                                *
*************************************************/

#include "POLAR_INTER.h"

int CheckWarn=0;

void WarnungInterpolation(char Profilname[], double SollValue, double GrenzValue, int MaReEtaCl,struct inputfile *input)    
{  
  char FileName[150] = {"\0"};
  FILE *fWarn; 
       
  strcpy(FileName,"LaufV1.lili.");
  strcat(FileName,input->LiLiVersionName); 
  strcat(FileName,"/VoGen_PolarInterpolatin_Warnings.out");
  fWarn=fopen(FileName,"a+");
  if (MaReEtaCl==0)
  {
    fprintf(fWarn,"\n Profile: %s \t Mach-value requested: %lf \t  Mach-value limit: %lf \n",Profilname, SollValue, GrenzValue);
  }
  if (MaReEtaCl==1)
  {
    fprintf(fWarn,"\n Profile: %s \t Reynoldsnumber requested: %lf \t  Reynoldsnumber limit: %lf \n",Profilname, SollValue, GrenzValue);
  }
  if (MaReEtaCl==2)
  {
    fprintf(fWarn,"\n Profile: %s \t Flapdeflection requested: %lf \t  Flapdeflection limit: %lf \n",Profilname, SollValue, GrenzValue);
  }
  if (MaReEtaCl==3)
  {
    fprintf(fWarn,"\n Profile: %s \t CL-value requested: %lf \t  CL-value limit: %lf \n",Profilname, SollValue, GrenzValue);
  } 
  if (MaReEtaCl==4)
  {
    fprintf(fWarn,"\n Profile: %s \t CL-value requested: %lf \t  CL-value limit: %lf \nCD extrapolated!",Profilname, SollValue, GrenzValue);
  } 
  if (MaReEtaCl==5)
  {
    fprintf(fWarn,"\n Profile: %s \t CL-value requested: %lf \t  CL-value limit: %lf \nCD extrapolated! Penalty added",Profilname, SollValue, GrenzValue);
  } 
  if (MaReEtaCl==6)
  {
    fprintf(fWarn,"\n Profile: %s \t CL-value requested: %lf \t  CL-value limit: %lf \n",Profilname, SollValue, GrenzValue);
  } 

  if (CheckWarn==0)
  {
    CheckWarn=1;
    printf("\n Warnings occured during profile-drag integration!\n\tCheck VoGen_PolarInterpolatin_Warnings.out!\n");
  } 
  if (input->ClMinMaxHandling==1)
  {
    CheckWarn=2;
  }
  fclose(fWarn);
}

/*Die Sotierfunktionen habe ich aus "C von A bis Z von J�rgen Wolf
Das umfassende Handbuch f�r Linux, Unix und Windows
\u2013 2., aktualisierte und erweiterte Auflage 2006" 
und entsprechend angepasst und korigiert*/

void selectionSortPolar(struct EinWert *array, int elemente) 
{
   int index,index_klein,wert; 
   struct EinWert wert_klein;
   /* Schleife durchl�uft von links nach rechts */
   for(index = 0; index < elemente; index++) 
   {
      /* Aktuelle Position */
      wert=index;
      /* Schleife durchl�uft bis kleineres Element als
       * aktuelle Pos. gefunden wurde oder bis zum Ende,
       * was bedeutet, die aktuelle Position ist schon
       * das kleinste */
      for(index_klein = index+1; index_klein <= elemente;index_klein++) 
      { /* Ein kleineres Element gefunden? */
         if(array[index_klein].Alfa < array[wert].Alfa)
            /* Neues kleinstes Element */
            wert=index_klein;
      }
      /* Kleinstes Element an die aktuelle
       * Position falls n�tig */
      if(wert != index) 
      {
         wert_klein.Alfa =array[wert].Alfa;
         wert_klein.Ca   =array[wert].Ca;
         wert_klein.Cd   =array[wert].Cd;
         wert_klein.CDP  =array[wert].CDP;
         wert_klein.CDPW =array[wert].CDPW;
         wert_klein.CDPP =array[wert].CDPP;
         wert_klein.CDPF =array[wert].CDPF;
         wert_klein.Cm   =array[wert].Cm;         
         
         array[wert].Alfa  =array[index].Alfa;
         array[wert].Ca    =array[index].Ca;
         array[wert].Cd    =array[index].Cd;
         array[wert].CDP   =array[index].CDP;
         array[wert].CDPW  =array[index].CDPW;
         array[wert].CDPP  =array[index].CDPP;
         array[wert].CDPF  =array[index].CDPF;
         array[wert].Cm    =array[index].Cm;
                           
         array[index].Alfa  =wert_klein.Alfa;
         array[index].Ca    =wert_klein.Ca;  
         array[index].Cd    =wert_klein.Cd;  
         array[index].CDP   =wert_klein.CDP; 
         array[index].CDPW  =wert_klein.CDPW;
         array[index].CDPP  =wert_klein.CDPP;
         array[index].CDPF  =wert_klein.CDPF;
         array[index].Cm    =wert_klein.Cm;
         index --;  
      }
   }
}
void slectionSortMa(struct PolarenMa **MaStart,int elemente)
{
  int index,index_klein,wert;
  struct PolarenMa *MaRunner;
  struct PolarenMa *MaRunnerStart;
  struct PolarenMa *MaRunnerZwi;
  struct PolarenMa **MaIndexArray;
  
  MaIndexArray=malloc((elemente+1)*sizeof(struct PolarenMa*));
  MaRunner=*MaStart;
  MaIndexArray[0]=MaRunner;
  for (index = 1; index <= elemente; index++) 
  {
    MaRunner=MaRunner->next;
    MaIndexArray[index]=MaRunner;
  }
  for (index = 0; index <= elemente; index++) 
  {
    MaRunner=MaIndexArray[index];
  }
  for(index = 0; index < elemente; index++) 
  {
    wert=index;
    MaRunnerStart=MaIndexArray[index];
    for(index_klein = index+1; index_klein <= elemente;index_klein++) 
    { 
      MaRunner=MaIndexArray[index_klein];
      if(MaRunner->MA_NUM < MaRunnerStart->MA_NUM)
      {
        wert=index_klein;
      }
    }    
    if(wert != index) 
    {
      MaRunnerZwi=MaIndexArray[wert];
      MaIndexArray[wert]=MaIndexArray[index];
      MaIndexArray[index]=MaRunnerZwi;
      
      *MaStart=MaIndexArray[0];
      MaRunner=MaIndexArray[0];
      MaRunner->next=MaIndexArray[1];
      MaRunner=MaIndexArray[1];
      for (wert=1; wert  <= elemente; wert++)
      {
        if (wert==elemente)
        {
          MaRunner->next=NULL;
        }else{
          MaRunner->next=MaIndexArray[wert+1];
          MaRunner=MaIndexArray[wert+1];          
        }
      }
      index--;
    }
  }
}

void slectionSortRe(struct PolarenRe **ReStart,int elemente)
{
  int index,index_klein,wert;
  struct PolarenRe *ReRunner;
  struct PolarenRe *ReRunnerStart;
  struct PolarenRe *ReRunnerZwi;
  struct PolarenRe **ReIndexArray;
  
  ReIndexArray=malloc((elemente+1)*sizeof(struct PolarenRe*));
  ReRunner=*ReStart;
  ReIndexArray[0]=ReRunner;
  for (index = 1; index <= elemente; index++) 
  {
    ReRunner=ReRunner->next;
    ReIndexArray[index]=ReRunner;
  }
  for (index = 0; index <= elemente; index++) 
  {
    ReRunner=ReIndexArray[index];
  }
  for(index = 0; index < elemente; index++) 
  {
    wert=index;
    ReRunnerStart=ReIndexArray[index];
    for(index_klein = index+1; index_klein <= elemente;index_klein++) 
    { 
      ReRunner=ReIndexArray[index_klein];
      if(ReRunner->RE_NUM < ReRunnerStart->RE_NUM)
      {
        wert=index_klein;
      }
    }    
    if(wert != index) 
    {
      ReRunnerZwi=ReIndexArray[wert];
      ReIndexArray[wert]=ReIndexArray[index];
      ReIndexArray[index]=ReRunnerZwi;
      
      *ReStart=ReIndexArray[0];
      ReRunner=ReIndexArray[0];
      ReRunner->next=ReIndexArray[1];
      ReRunner=ReIndexArray[1];
      for (wert=1; wert  <= elemente; wert++)
      {
        if (wert==elemente)
        {
          ReRunner->next=NULL;
        }else{
          ReRunner->next=ReIndexArray[wert+1];
          ReRunner=ReIndexArray[wert+1];          
        }
      }
      index--;
    }
  }
}


void slectionSortEta(struct EinePolare **EtaStart,int elemente)
{
  int index,index_klein,wert;
  struct EinePolare *EtaRunner;
  struct EinePolare *EtaRunnerStart;
  struct EinePolare *EtaRunnerZwi;
  struct EinePolare **EtaIndexArray;
  
  EtaIndexArray=malloc((elemente+1)*sizeof(struct EinePolare*));
  EtaRunner=*EtaStart;
  EtaIndexArray[0]=EtaRunner;
  for (index = 1; index <= elemente; index++) 
  {
    EtaRunner=EtaRunner->next;
    EtaIndexArray[index]=EtaRunner;
  }
  for (index = 0; index <= elemente; index++) 
  {
    EtaRunner=EtaIndexArray[index];
  }
  for(index = 0; index < elemente; index++) 
  {
    wert=index;
    EtaRunnerStart=EtaIndexArray[index];
    for(index_klein = index+1; index_klein <= elemente;index_klein++) 
    { 
      EtaRunner=EtaIndexArray[index_klein];
      if(EtaRunner->Eta_NUM < EtaRunnerStart->Eta_NUM)
      {
        wert=index_klein;
      }
    }    
    if(wert != index) 
    {
      EtaRunnerZwi=EtaIndexArray[wert];
      EtaIndexArray[wert]=EtaIndexArray[index];
      EtaIndexArray[index]=EtaRunnerZwi;
      
      *EtaStart=EtaIndexArray[0];
      EtaRunner=EtaIndexArray[0];
      EtaRunner->next=EtaIndexArray[1];
      EtaRunner=EtaIndexArray[1];
      for (wert=1; wert  <= elemente; wert++)
      {
        if (wert==elemente)
        {
          EtaRunner->next=NULL;
        }else{
          EtaRunner->next=EtaIndexArray[wert+1];
          EtaRunner=EtaIndexArray[wert+1];          
        }
      }
      index--;
    }
  }
}

struct EinWert *getETA(double Klappe, int I,struct PolarenRe *CurRe)
{
  struct EinePolare *EtaRunner;
  int schalter,control;
  
  if (CurRe->anzahlEtaNum==0 && CurRe->EtaPolare==NULL)
  {
    CurRe->EtaPolare=(struct EinePolare *)malloc(sizeof(struct EinePolare));
    EtaRunner=CurRe->EtaPolare;
    EtaRunner->next=NULL;
    EtaRunner->anzahlWerte=I;
    EtaRunner->Eta_NUM=Klappe;    
    EtaRunner->Werte=(struct EinWert *)malloc(I*sizeof(struct EinWert));;    
    CurRe->anzahlEtaNum++;
  }else{
    EtaRunner=CurRe->EtaPolare;
    schalter=0;
    control=0;
    do
    {
      if (schalter==1)
      {
            EtaRunner=EtaRunner->next;
      }
      if (EtaRunner->Eta_NUM==Klappe)
      {
       printf("ReNum:%lf \t Flap:%lf�\n",CurRe->RE_NUM, EtaRunner->Eta_NUM);
       errmessage(55);
       control=1;
      }
      schalter=1;
    }while (EtaRunner->next!=NULL);
    if (control==0)
    {
      EtaRunner->next=(struct EinePolare *)malloc(sizeof(struct EinePolare));
      EtaRunner=EtaRunner->next;
      EtaRunner->next=NULL;
      EtaRunner->anzahlWerte=I;
      EtaRunner->Eta_NUM=Klappe;    
      EtaRunner->Werte=(struct EinWert *)malloc(I*sizeof(struct EinWert));;    
      CurRe->anzahlEtaNum++;
    }
  }
  if (CurRe->anzahlEtaNum > 1)
  {
    slectionSortEta(&CurRe->EtaPolare,CurRe->anzahlEtaNum-1);
  }
  return (EtaRunner->Werte);  
}

struct EinWert *getRE(double RE,double Klappe, int I,struct PolarenMa *CurMa)
{
   struct PolarenRe *ReRunner;
   struct EinWert *WritePolare;     
   int schalter,control;
   if (CurMa->anzahlReNum==0 && CurMa->RePolaren==NULL)
   {
     CurMa->RePolaren=(struct PolarenRe *)malloc(sizeof(struct PolarenRe));
     ReRunner=CurMa->RePolaren;
     ReRunner->next=NULL;
     ReRunner->anzahlEtaNum=0;
     ReRunner->RE_NUM=RE;
     ReRunner->EtaPolare=NULL;
     WritePolare=getETA(Klappe,I,ReRunner);
     CurMa->anzahlReNum++;  
   }else{
     ReRunner=CurMa->RePolaren;     
     schalter=0;
     control=0;
     do
     {
       if (schalter==1)
       {
              ReRunner=ReRunner->next;
       }
       if (ReRunner->RE_NUM==RE)
       {
              WritePolare=getETA(Klappe,I,ReRunner);
              control=1;
       }
       schalter=1;
     }while (ReRunner->next!=NULL);
     if (control==0)
     {
       ReRunner->next=(struct PolarenRe *)malloc(sizeof(struct PolarenRe));       
       ReRunner=ReRunner->next;
       ReRunner->next=NULL;
       ReRunner->RE_NUM=RE;
       ReRunner->anzahlEtaNum=0;
       ReRunner->EtaPolare=NULL;
       WritePolare=getETA(Klappe,I,ReRunner);
       CurMa->anzahlReNum++;  
     }             
   }
   if (CurMa->anzahlReNum>1)
   {
     slectionSortRe(&CurMa->RePolaren,CurMa->anzahlReNum-1);
   }
   return (WritePolare);
}


void readProfil (struct EinSchnitt *Schnitt, char PROFIPFAD[],struct EinProfiel *Cur_Profiel,struct inputfile *input)
{
  FILE *fopen(),*ffile;  
  char dummy[MaxZeilenLaenge] = {"\0"};
  char dummy2[100] = {"\0"};
  char ZEILE_READ[MaxZeilenLaenge]; 
  double Mach, Re, Eta;  
  int Inumber,schalter,control, lauf, lauf2;
  struct PolarenMa *MaRunner;
  struct EinWert *WritePolare; 
  char *ptr;
  
  //printf("Hallo %s\n\n",Schnitt->ProfilName);
  strcpy(Cur_Profiel->Profiel_Name,Schnitt->ProfilName);
  Cur_Profiel->MaPolare=NULL;
  Cur_Profiel->AnzMaNummern=0;
  //printf("dummy: %s\n\n",dummy);
  sscanf(PROFIPFAD, "%s", &dummy[0]);
  sscanf(Schnitt->ProfilName, "%s", &dummy2[0]);
  strcat(dummy,dummy2);
  strcat(dummy,".plt");
  //printf("Oefne Polare: %s\n",dummy);
  ffile = fopen(dummy,"r");
  if (ffile == NULL)
  {
    printf("\nError: Konnte Datei: %s fuer die Polareniterpolation nicht oeffnen!\n", dummy);
    errmessage(1);
  }
  do
  {
    fgets (ZEILE_READ,MaxZeilenLaenge,ffile);
    if( strstr(ZEILE_READ, "ZONE") != NULL)
    {   
      Mach=0.0; Re=0.0; Eta=0.0;
      ptr = strtok(ZEILE_READ, ";\"= ");
      while(ptr != NULL) 
      {
         
        ptr = strtok(NULL, ";\"= ");
        //printf("ptr: %s\n",ptr);
        if( strstr(ptr, "Mach") != NULL)
        {
          ptr = strtok(NULL, ";\"= ");
          sscanf(ptr, "%lf", &Mach);        
        }
        if( strstr(ptr, "Reynoldszahl") != NULL)
        {
          ptr = strtok(NULL, ";\"= ");
          sscanf(ptr, "%lf", &Re);
        }        
        if( strstr(ptr, "Klappe") != NULL)
        {
          ptr = strtok(NULL, ";\"= ");
          sscanf(ptr, "%lf", &Eta);
        }        
        if( strstr(ptr, "I") != NULL)
        {
          ptr = strtok(NULL, ";\"= ");
          sscanf(ptr, "%d", &Inumber);
          break;
        }                        
      }
      if (Cur_Profiel->MaPolare==NULL && Cur_Profiel->AnzMaNummern==0)
      {
        Cur_Profiel->MaPolare=(struct PolarenMa *)malloc(sizeof(struct PolarenMa));
        MaRunner=Cur_Profiel->MaPolare;
        MaRunner->next=NULL;
        MaRunner->MA_NUM=Mach;
        MaRunner->anzahlReNum=0;
        MaRunner->RePolaren=NULL;
        WritePolare=getRE(Re,Eta,Inumber,MaRunner);
        Cur_Profiel->AnzMaNummern++;
      }else{
        MaRunner=Cur_Profiel->MaPolare;
        schalter=0;
        control=0;
        do
        {
          if (schalter==1)
          {
            MaRunner=MaRunner->next;
          }
          if (MaRunner->MA_NUM==Mach)
          {
            WritePolare=getRE(Re,Eta,Inumber,MaRunner);
            control=1;
          }
          schalter=1;
        }while (MaRunner->next!=NULL);
        if (control==0)
        {
          MaRunner->next=(struct PolarenMa *)malloc(sizeof(struct PolarenMa));
          MaRunner=MaRunner->next;
          MaRunner->next=NULL;
          MaRunner->MA_NUM=Mach;
          MaRunner->anzahlReNum=0;
          MaRunner->RePolaren=NULL;
          WritePolare=getRE(Re,Eta,Inumber,MaRunner);
          Cur_Profiel->AnzMaNummern++;
        }              
      }
      //printf ("Mach: %lf\tRe: %lf\tKlappe: %lf\tI: %d\n", Mach, Re, Eta, Inumber);          
      for (lauf=0; lauf < Inumber; lauf++)
      {
        int readCM,readCDPW,readCDPP,readCDPF, check;
        check=0;
        if (feof(ffile))
        {
          check=1;
        }else{
          fgets (ZEILE_READ,MaxZeilenLaenge,ffile);
        }
        if( strstr(ZEILE_READ, "ZONE") != NULL  || check==1)
        {
           printf("Im Profil \"%s\" Mach-Zahl:%lf  Re:%lf  Klappe:%lf  enthaelt weniger Polarenpunkte als erwartet\n",Cur_Profiel->Profiel_Name, Mach, Re, Eta);
           errmessage(1);
        }
        do
        {
          ptr = strtok(ZEILE_READ, " \t");
        }
        while ((ptr != NULL) && (strlen(ptr)==0));
        lauf2=0;      
        if (input->useCm) {readCM=3;}else{readCM=2;}
        if (input->useCdpw) {readCDPW=readCM+1;}else{readCDPW=readCM;}
        if (input->useCdpp) {readCDPP=readCDPW+1;}else{readCDPP=readCDPW;}
        if (input->useCdpf) {readCDPF=readCDPP+1;}else{readCDPF=readCDPP;}
        check=0;
        while(ptr != NULL) 
        {
          if (lauf2==0)        {sscanf(ptr, "%lf", &WritePolare[lauf].Alfa);}
          if (lauf2==1)        {sscanf(ptr, "%lf", &WritePolare[lauf].Cd);}
          if (lauf2==2)        {sscanf(ptr, "%lf", &WritePolare[lauf].Ca);}          
          if (input->useCm) {if (lauf2==readCM) {sscanf(ptr, "%lf", &WritePolare[lauf].Cm);}}
          if (input->useCdpw) {if (lauf2==readCDPW) {sscanf(ptr, "%lf", &WritePolare[lauf].CDPW);}}
          if (input->useCdpp) {if (lauf2==readCDPP) {sscanf(ptr, "%lf", &WritePolare[lauf].CDPP);}}          
          if (input->useCdpf) {if (lauf2==readCDPF) {sscanf(ptr, "%lf", &WritePolare[lauf].CDPF);}}          

          do
          {
            ptr = strtok(NULL, " \t");
          }
          while ((ptr != NULL) && (strlen(ptr)==0));
          lauf2++;
          check=0;
        }
      }
      selectionSortPolar(WritePolare, Inumber-1); 
    }        
  } while (!feof(ffile)); 
  fclose(ffile); 
  /*Sortieren*/
  if (Cur_Profiel->AnzMaNummern > 1)
  {
    slectionSortMa(&Cur_Profiel->MaPolare,Cur_Profiel->AnzMaNummern-1);
  } 
}

void print_screan_debug_Polaren(struct Profiele AlleProfile, struct inputfile *input)
{
  int lauf;
  struct EinProfiel *EinP_runner;
  struct PolarenMa *Ma_runner;
  struct PolarenRe *Re_runner;
  struct EinePolare *Eta_runner;
  
  EinP_runner=AlleProfile.Profiel;  
  do
  {
    printf("\nProfile-Name: %s\n", EinP_runner->Profiel_Name);
    Ma_runner=EinP_runner->MaPolare;
    do
    {
      printf("MachZahl: %lf\n", Ma_runner->MA_NUM);
      Re_runner=Ma_runner->RePolaren;
      do
      {
        printf("ReZahl: %lf\n", Re_runner->RE_NUM);
        Eta_runner=Re_runner->EtaPolare;
        do
        {
          printf("Klappen-Winkel: %lf\n", Eta_runner->Eta_NUM);
          for (lauf=0;lauf<Eta_runner->anzahlWerte;lauf++)
          {
            printf("%lf\t%lf\t%lf",Eta_runner->Werte[lauf].Alfa,Eta_runner->Werte[lauf].Ca,Eta_runner->Werte[lauf].Cd);
            if (input->useCm) {printf("\t%lf",Eta_runner->Werte[lauf].Cm);}
            if (input->useCdpw) {printf("\t%lf",Eta_runner->Werte[lauf].CDPW);}                           
            if (input->useCdpp) {printf("\t%lf",Eta_runner->Werte[lauf].CDPP);}                             
            if (input->useCdpf) {printf("\t%lf",Eta_runner->Werte[lauf].CDPF);}                             
            printf("\n");                                                                                 
          }
          printf("\n");
          Eta_runner=Eta_runner->next;
        }while(Eta_runner!=NULL);         
        Re_runner=Re_runner->next;
      }while(Re_runner!=NULL);     
      Ma_runner=Ma_runner->next;
    }while(Ma_runner!=NULL);    
    EinP_runner=EinP_runner->next;
  }while(EinP_runner!=NULL);
}



struct EinWert getProfilBeiwertePOLARE(struct EinePolare *ETA_runner,double CA, struct inputfile *input,char Profilname[])
{
  int lauf, s1;
  struct EinWert *WERT1,WERTR;
  
  //Debug
  //printf("Rein:getProfilBeiwertePOLARE\n");
    
  WERT1=ETA_runner->Werte;  
  
  //Debug
  //printf("CA:%lf ProfName:%s ETA_Prof:%lf\n",CA,Profilname,ETA_runner->Eta_NUM);
  
  s1=0;
  for(lauf=0;lauf<ETA_runner->anzahlWerte;lauf++)
  {
    //Debug
    //printf ("CA_such:%lf CA_ist:%lf CW_ist:%lf lauf:%d Antal Werte:%d\n",CA,WERT1[lauf].Ca,WERT1[lauf].Cd,lauf,ETA_runner->anzahlWerte);
    if (((lauf+1)==ETA_runner->anzahlWerte && s1==0) || (WERT1[lauf].Ca>=CA && s1==0))
    {
      s1=1;
      if (input->useCm) {WERTR.Cm=WERT1[lauf].Cm;}
      if (input->useCdpw) {WERTR.CDPW=WERT1[lauf].CDPW;}  
      if (input->useCdpp) {WERTR.CDPP=WERT1[lauf].CDPP;}  
      if (input->useCdpf) {WERTR.CDPF=WERT1[lauf].CDPF;}  
      WERTR.Alfa=WERT1[lauf].Alfa;
      if (WERT1[lauf].Ca!=CA)
      {
        
        float a[3][3],b[3];
        switch (input->ClMinMaxHandling)
        {
        case 0: 
            WarnungInterpolation(Profilname, CA, WERT1[lauf].Ca,3,input);
            break;
        case 1:  
            WarnungInterpolation(Profilname, CA, WERT1[lauf].Ca,3,input);
            break;
        case 2: 
            WarnungInterpolation(Profilname, CA, WERT1[lauf].Ca,4,input);
            if (lauf==0)
            {
              WERTR.CDP=(WERT1[0].Cd-WERT1[1].Cd)/(WERT1[0].Ca-WERT1[1].Ca)*(CA-WERT1[0].Ca)+WERT1[0].Cd;
            }else{
              WERTR.CDP=(WERT1[lauf].Cd-WERT1[lauf-1].Cd)/(WERT1[lauf].Ca-WERT1[lauf-1].Ca)*(CA-WERT1[lauf].Ca)+WERT1[lauf].Cd;            
            }
            break;
        case 3: 
            WarnungInterpolation(Profilname, CA, WERT1[lauf].Ca,5,input);
            if (lauf==0)
            {
              WERTR.CDP=(WERT1[0].Cd-WERT1[1].Cd)/(WERT1[0].Ca-WERT1[1].Ca)*(CA-WERT1[0].Ca)*input->PENALTYfact+WERT1[0].Cd;
            }else{
              WERTR.CDP=(WERT1[lauf].Cd-WERT1[lauf-1].Cd)/(WERT1[lauf].Ca-WERT1[lauf-1].Ca)*(CA-WERT1[lauf].Ca)*input->PENALTYfact+WERT1[lauf].Cd;            
            }
            break;
        case 4: 
            WarnungInterpolation(Profilname, CA, WERT1[lauf].Ca,4,input);
            if (lauf==0)
            {
              a[0][0]=WERT1[0].Ca*WERT1[0].Ca;
              a[0][1]=WERT1[0].Ca;
              a[0][2]=1.0;
              b[0]=WERT1[0].Cd;
              
              a[1][0]=WERT1[1].Ca*WERT1[1].Ca;
              a[1][1]=WERT1[1].Ca;
              a[1][2]=1.0;
              b[1]=WERT1[1].Cd;

              a[2][0]=WERT1[2].Ca*WERT1[2].Ca;
              a[2][1]=WERT1[2].Ca;
              a[2][2]=1.0;
              b[2]=WERT1[2].Cd;
              
              LinGlLoesung(&a[0][0],&b[0],3);                           
              
              WERTR.CDP=b[0]*CA*CA+b[1]*CA+b[2];
            }else{
              a[0][0]=WERT1[lauf-2].Ca*WERT1[lauf-2].Ca;
              a[0][1]=WERT1[lauf-2].Ca;
              a[0][2]=1.0;
              b[0]=WERT1[lauf-2].Cd;
              
              a[1][0]=WERT1[lauf-1].Ca*WERT1[lauf-1].Ca;
              a[1][1]=WERT1[lauf-1].Ca;
              a[1][2]=1.0;
              b[1]=WERT1[lauf-1].Cd;

              a[2][0]=WERT1[lauf].Ca*WERT1[lauf].Ca;
              a[2][1]=WERT1[lauf].Ca;
              a[2][2]=1.0;
              b[2]=WERT1[lauf].Cd;
              
              LinGlLoesung(&a[0][0],&b[0],3);                           
               
              WERTR.CDP=b[0]*CA*CA+b[1]*CA+b[2];
            }
            break;
        case 5: 
            WarnungInterpolation(Profilname, CA, WERT1[lauf].Ca,5,input);
            if (lauf==0)
            {
              a[0][0]=WERT1[0].Ca*WERT1[0].Ca;
              a[0][1]=WERT1[0].Ca;
              a[0][2]=1.0;
              b[0]=WERT1[0].Cd;
              
              a[1][0]=WERT1[1].Ca*WERT1[1].Ca;
              a[1][1]=WERT1[1].Ca;
              a[1][2]=1.0;
              b[1]=WERT1[1].Cd;

              a[2][0]=WERT1[2].Ca*WERT1[2].Ca;
              a[2][1]=WERT1[2].Ca;
              a[2][2]=1.0;
              b[2]=WERT1[2].Cd;
              
              LinGlLoesung(&a[0][0],&b[0],3);                           
              
              WERTR.CDP=b[0]*CA*CA+b[1]*CA+b[2]+fabs(WERT1[0].Cd-(b[0]*CA*CA+b[1]*CA+b[2]))*input->PENALTYfact;
            }else{
              a[0][0]=WERT1[lauf-2].Ca*WERT1[lauf-2].Ca;
              a[0][1]=WERT1[lauf-2].Ca;
              a[0][2]=1.0;
              b[0]=WERT1[lauf-2].Cd;
              
              a[1][0]=WERT1[lauf-1].Ca*WERT1[lauf-1].Ca;
              a[1][1]=WERT1[lauf-1].Ca;
              a[1][2]=1.0;
              b[1]=WERT1[lauf-1].Cd;

              a[2][0]=WERT1[lauf].Ca*WERT1[lauf].Ca;
              a[2][1]=WERT1[lauf].Ca;
              a[2][2]=1.0;
              b[2]=WERT1[lauf].Cd;
              
              LinGlLoesung(&a[0][0],&b[0],3);                           
              
              WERTR.CDP=b[0]*CA*CA+b[1]*CA+b[2]+fabs(WERT1[lauf].Cd-b[0]*CA*CA+b[1]*CA+b[2])*input->PENALTYfact;
            }    
            break;
        }       
      }else{
        WERTR.CDP=WERT1[lauf].Cd;
      }
    }else{
      if (s1==0)
      {
        if ((WERT1[lauf+1].Ca>=CA) && (WERT1[lauf].Ca < CA))
        {
          s1=1;
          WERTR.CDP=(WERT1[lauf+1].Cd-WERT1[lauf].Cd)/(WERT1[lauf+1].Ca-WERT1[lauf].Ca)*(CA-WERT1[lauf].Ca)+WERT1[lauf].Cd;
          if (input->useCm) {WERTR.Cm=(WERT1[lauf+1].Cm-WERT1[lauf].Cm)/(WERT1[lauf+1].Ca-WERT1[lauf].Ca)*(CA-WERT1[lauf].Ca)+WERT1[lauf].Cm;}
          if (input->useCdpw) {WERTR.CDPW=(WERT1[lauf+1].CDPW-WERT1[lauf].CDPW)/(WERT1[lauf+1].Ca-WERT1[lauf].Ca)*(CA-WERT1[lauf].Ca)+WERT1[lauf].CDPW;}
          if (input->useCdpp) {WERTR.CDPP=(WERT1[lauf+1].CDPP-WERT1[lauf].CDPP)/(WERT1[lauf+1].Ca-WERT1[lauf].Ca)*(CA-WERT1[lauf].Ca)+WERT1[lauf].CDPP;}
          if (input->useCdpf) {WERTR.CDPF=(WERT1[lauf+1].CDPF-WERT1[lauf].CDPF)/(WERT1[lauf+1].Ca-WERT1[lauf].Ca)*(CA-WERT1[lauf].Ca)+WERT1[lauf].CDPF;}
          WERTR.Alfa=(WERT1[lauf+1].Alfa-WERT1[lauf].Alfa)/(WERT1[lauf+1].Ca-WERT1[lauf].Ca)*(CA-WERT1[lauf].Ca)+WERT1[lauf].Alfa;;
        }
      }
    }
  }
  //Debug
  //printf("POLAREsearch Cdp CM Alfa: %lf %lf %lf\n",WERTR.CDP, WERTR.Cm, WERTR.Alfa);
  
  return(WERTR);
  
  //Debug
  //printf("Raus:getProfilBeiwertePOLARE\n");
}


struct EinWert getProfilBeiwerteETA(struct PolarenRe *Re_runner,double CA, double Flap_Set,struct inputfile *input,char Profilname[])
{
  struct EinePolare *ETA_runner;
  int s1;
  struct EinWert WERT1,WERT2,WERTR;
  
  //Debug
  //printf("Rein:getProfilBeiwerteETA\n");
  
  ETA_runner=Re_runner->EtaPolare;
  
  //DEBUG
  //printf("RePolare:%lf\n",Re_runner->RE_NUM);
  
  s1=0;
  do
  {
    if ((ETA_runner->next==NULL && s1==0) || (ETA_runner->Eta_NUM>=Flap_Set && s1==0))
    {
      s1=1;
      WERT1=getProfilBeiwertePOLARE(ETA_runner, CA, input,Profilname);
      WERTR.CDP=WERT1.CDP;
      if (input->useCm) {WERTR.Cm=WERT1.Cm;}
      if (input->useCdpw) {WERTR.CDPW=WERT1.CDPW;}  
      if (input->useCdpp) {WERTR.CDPP=WERT1.CDPP;}  
      if (input->useCdpf) {WERTR.CDPF=WERT1.CDPF;}  
      WERTR.Alfa=WERT1.Alfa;
      if (ETA_runner->Eta_NUM!=Flap_Set)
      {
        WarnungInterpolation(Profilname, Flap_Set, ETA_runner->Eta_NUM,2,input); 
      }
    }else{
      if (s1==0)
      {
        if ((ETA_runner->next->Eta_NUM>=Flap_Set) && (ETA_runner->Eta_NUM < Flap_Set))
        {
          s1=1;
          WERT1=getProfilBeiwertePOLARE(ETA_runner, CA, input,Profilname);
          WERT2=getProfilBeiwertePOLARE(ETA_runner->next, CA, input,Profilname);
          WERTR.CDP=(WERT2.CDP-WERT1.CDP)/(ETA_runner->next->Eta_NUM-ETA_runner->Eta_NUM)*(Flap_Set-ETA_runner->Eta_NUM)+WERT1.CDP;         
          if (input->useCm) {WERTR.Cm=(WERT2.Cm-WERT1.Cm)/(ETA_runner->next->Eta_NUM-ETA_runner->Eta_NUM)*(Flap_Set-ETA_runner->Eta_NUM)+WERT1.Cm;}
          if (input->useCdpw) {WERTR.CDPW=(WERT2.CDPW-WERT1.CDPW)/(ETA_runner->next->Eta_NUM-ETA_runner->Eta_NUM)*(Flap_Set-ETA_runner->Eta_NUM)+WERT1.CDPW;}
          if (input->useCdpp) {WERTR.CDPP=(WERT2.CDPP-WERT1.CDPP)/(ETA_runner->next->Eta_NUM-ETA_runner->Eta_NUM)*(Flap_Set-ETA_runner->Eta_NUM)+WERT1.CDPP;}
          if (input->useCdpf) {WERTR.CDPF=(WERT2.CDPF-WERT1.CDPF)/(ETA_runner->next->Eta_NUM-ETA_runner->Eta_NUM)*(Flap_Set-ETA_runner->Eta_NUM)+WERT1.CDPF;}
          WERTR.Alfa=(WERT2.Alfa-WERT1.Alfa)/(ETA_runner->next->Eta_NUM-ETA_runner->Eta_NUM)*(Flap_Set-ETA_runner->Eta_NUM)+WERT1.Alfa;;
        }
      }
    }
    ETA_runner=ETA_runner->next;
  }while(ETA_runner!=NULL);
  //printf("ETAsearch Cdp CM Alfa: %lf %lf %lf\n",WERTR.CDP, WERTR.Cm, WERTR.Alfa);
  return(WERTR);
  //Debug
  //printf("Raus:getProfilBeiwerteETA\n");
}


struct EinWert getProfilBeiwerteRE(struct PolarenMa *Ma_runner,double CA, double RE_LOC, double Flap_Set,struct inputfile *input,char Profilname[])
{
  struct PolarenRe *Re_runner;
  int s1;
  struct EinWert WERT1,WERT2,WERTR;
  
  //Debug
  //printf("Rein:getProfilBeiwerteRE\n");
  //printf("CA:%lf ProfName:%s Eta:%lf RE_LOC:%lf\n",CA,Profilname,Flap_Set,RE_LOC);
  
  Re_runner=Ma_runner->RePolaren;
  s1=0;
  do
  {
    if ((Re_runner->next==NULL && s1==0) || (Re_runner->RE_NUM>=RE_LOC && s1==0))
    {
      s1=1;
      WERT1=getProfilBeiwerteETA(Re_runner, CA, Flap_Set, input,Profilname);
      WERTR.CDP=WERT1.CDP;
      if (input->useCm) {WERTR.Cm=WERT1.Cm;}
      if (input->useCdpw) {WERTR.CDPW=WERT1.CDPW;} 
      if (input->useCdpp) {WERTR.CDPP=WERT1.CDPP;} 
      if (input->useCdpf) {WERTR.CDPF=WERT1.CDPF;} 
      WERTR.Alfa=WERT1.Alfa;
      if (Re_runner->RE_NUM>=RE_LOC)
      {
        WarnungInterpolation(Profilname, RE_LOC, Re_runner->RE_NUM,1,input);   
      }
    }else{
      if (s1==0)
      {
        if ((Re_runner->next->RE_NUM>=RE_LOC) && (Re_runner->RE_NUM < RE_LOC))
        {
          s1=1;
          WERT1=getProfilBeiwerteETA(Re_runner, CA, Flap_Set, input,Profilname);
          WERT2=getProfilBeiwerteETA(Re_runner->next, CA, Flap_Set, input,Profilname);
          WERTR.CDP=(WERT2.CDP-WERT1.CDP)/(Re_runner->next->RE_NUM-Re_runner->RE_NUM)*(RE_LOC-Re_runner->RE_NUM)+WERT1.CDP;
          if (input->useCm) {WERTR.Cm=(WERT2.Cm-WERT1.Cm)/(Re_runner->next->RE_NUM-Re_runner->RE_NUM)*(RE_LOC-Re_runner->RE_NUM)+WERT1.Cm;}
          if (input->useCdpw) {WERTR.CDPW=(WERT2.CDPW-WERT1.CDPW)/(Re_runner->next->RE_NUM-Re_runner->RE_NUM)*(RE_LOC-Re_runner->RE_NUM)+WERT1.CDPW;}
          if (input->useCdpp) {WERTR.CDPP=(WERT2.CDPP-WERT1.CDPP)/(Re_runner->next->RE_NUM-Re_runner->RE_NUM)*(RE_LOC-Re_runner->RE_NUM)+WERT1.CDPP;}
          if (input->useCdpf) {WERTR.CDPF=(WERT2.CDPF-WERT1.CDPF)/(Re_runner->next->RE_NUM-Re_runner->RE_NUM)*(RE_LOC-Re_runner->RE_NUM)+WERT1.CDPF;}
          WERTR.Alfa=(WERT2.Alfa-WERT1.Alfa)/(Re_runner->next->RE_NUM-Re_runner->RE_NUM)*(RE_LOC-Re_runner->RE_NUM)+WERT1.Alfa;;
        }
      }
    }
    Re_runner=Re_runner->next;
  }while(Re_runner!=NULL);
  return(WERTR);
  //Debug
  //printf("Raus:getProfilBeiwerteRE\n");
}

struct EinWert getProfilBeiwerteMA( struct EinProfiel *EinP_runner, double Mach, double CA, double RE_LOC, double Flap_Set,struct inputfile *input,char Profilname[])
{
  int s1;
  struct EinWert WERT1,WERT2,WERTR;
  struct PolarenMa *Ma_runner;

  Ma_runner=EinP_runner->MaPolare;
  s1=0;
  do
  {
    if ((Ma_runner->next==NULL && s1==0) || (Ma_runner->MA_NUM>=Mach && s1==0))
    {
      s1=1;
      WERT1=getProfilBeiwerteRE(Ma_runner, CA, RE_LOC, Flap_Set, input,Profilname);
      WERTR.CDP=WERT1.CDP;
      if (input->useCm) {WERTR.Cm=WERT1.Cm;}
      if (input->useCdpw) {WERTR.CDPW=WERT1.CDPW;}
      if (input->useCdpp) {WERTR.CDPP=WERT1.CDPP;}
      if (input->useCdpf) {WERTR.CDPF=WERT1.CDPF;}
      WERTR.Alfa=WERT1.Alfa;
      if (Ma_runner->MA_NUM!=Mach)
      {
        WarnungInterpolation(EinP_runner->Profiel_Name, Mach, Ma_runner->MA_NUM,0,input);
      }
    }else{
      if (s1==0)
      {
        if ((Ma_runner->next->MA_NUM > Mach) && (Ma_runner->MA_NUM < Mach))
        {
          s1=1;
          WERT1=getProfilBeiwerteRE(Ma_runner, CA, RE_LOC, Flap_Set, input,Profilname);
          WERT2=getProfilBeiwerteRE(Ma_runner->next, CA, RE_LOC, Flap_Set, input,Profilname);
          WERTR.CDP=(WERT2.CDP-WERT1.CDP)/(Ma_runner->next->MA_NUM-Ma_runner->MA_NUM)*(Mach-Ma_runner->MA_NUM)+WERT1.CDP;
          if (input->useCm) {WERTR.Cm=(WERT2.Cm-WERT1.Cm)/(Ma_runner->next->MA_NUM-Ma_runner->MA_NUM)*(Mach-Ma_runner->MA_NUM)+WERT1.Cm;}
          if (input->useCdpw) {WERTR.CDPW=(WERT2.CDPW-WERT1.CDPW)/(Ma_runner->next->MA_NUM-Ma_runner->MA_NUM)*(Mach-Ma_runner->MA_NUM)+WERT1.CDPW;}
          if (input->useCdpp) {WERTR.CDPP=(WERT2.CDPP-WERT1.CDPP)/(Ma_runner->next->MA_NUM-Ma_runner->MA_NUM)*(Mach-Ma_runner->MA_NUM)+WERT1.CDPP;}
          if (input->useCdpf) {WERTR.CDPF=(WERT2.CDPF-WERT1.CDPF)/(Ma_runner->next->MA_NUM-Ma_runner->MA_NUM)*(Mach-Ma_runner->MA_NUM)+WERT1.CDPF;}
          WERTR.Alfa=(WERT2.Alfa-WERT1.Alfa)/(Ma_runner->next->MA_NUM-Ma_runner->MA_NUM)*(Mach-Ma_runner->MA_NUM)+WERT1.Alfa;
        }
      }
    }
    Ma_runner=Ma_runner->next;
  }while(Ma_runner!=NULL);
  return(WERTR);
}

void getProfilBeiwerte(struct VERTEILUNG_POLARint *Datensatz,struct Profiele AlleProfile ,double Mach,struct inputfile *input)
{
  struct EinWert WERT1,WERT2;
  struct EinProfiel *EinP_runner;

  //Debug
  //printf("Rein:getProfilBeiwerte\n");
  
  EinP_runner=AlleProfile.Profiel;
  do
  {
    if (strcmp (EinP_runner->Profiel_Name, Datensatz->SEC_Name)==0)
    {
      /*Ma_runner=EinP_runner->MaPolare;
      s1=0;
      do
      {
        if ((Ma_runner->next==NULL && s1==0) || (Ma_runner->MA_NUM>=Mach && s1==0))
        {
          s1=1;
          WERT1=getProfilBeiwerteRE(Ma_runner, Datensatz->CA, Datensatz->RE_LOC, Datensatz->Flap_Set, input,EinP_runner->Profiel_Name);
          Datensatz->CDP=WERT1.CDP;
          if (input->useCm) {Datensatz->CM=WERT1.Cm;}
          if (input->useCdpw) {Datensatz->CDPW=WERT1.CDPW;}
          if (input->useCdpp) {Datensatz->CDPP=WERT1.CDPP;}
          if (input->useCdpf) {Datensatz->CDPF=WERT1.CDPF;}
          Datensatz->Alfa=WERT1.Alfa;
          Datensatz->CD=WERT1.CDP+Datensatz->CDI; 
          if (Ma_runner->MA_NUM!=Mach)
          {
            WarnungInterpolation(EinP_runner->Profiel_Name, Mach, Ma_runner->MA_NUM,0,input);       
          }
        }else{
          if (s1==0)
          {
            if ((Ma_runner->next->MA_NUM > Mach) && (Ma_runner->MA_NUM < Mach))
            {
              s1=1;
              WERT1=getProfilBeiwerteRE(Ma_runner, Datensatz->CA, Datensatz->RE_LOC, Datensatz->Flap_Set, input,EinP_runner->Profiel_Name);
              WERT2=getProfilBeiwerteRE(Ma_runner->next, Datensatz->CA, Datensatz->RE_LOC, Datensatz->Flap_Set, input,EinP_runner->Profiel_Name);
              Datensatz->CDP=(WERT2.CDP-WERT1.CDP)/(Ma_runner->next->MA_NUM-Ma_runner->MA_NUM)*(Mach-Ma_runner->MA_NUM)+WERT1.CDP;
              if (input->useCm) {Datensatz->CM=(WERT2.Cm-WERT1.Cm)/(Ma_runner->next->MA_NUM-Ma_runner->MA_NUM)*(Mach-Ma_runner->MA_NUM)+WERT1.Cm;}
              if (input->useCdpw) {Datensatz->CDPW=(WERT2.CDPW-WERT1.CDPW)/(Ma_runner->next->MA_NUM-Ma_runner->MA_NUM)*(Mach-Ma_runner->MA_NUM)+WERT1.CDPW;}
              if (input->useCdpp) {Datensatz->CDPP=(WERT2.CDPP-WERT1.CDPP)/(Ma_runner->next->MA_NUM-Ma_runner->MA_NUM)*(Mach-Ma_runner->MA_NUM)+WERT1.CDPP;}
              if (input->useCdpf) {Datensatz->CDPF=(WERT2.CDPF-WERT1.CDPF)/(Ma_runner->next->MA_NUM-Ma_runner->MA_NUM)*(Mach-Ma_runner->MA_NUM)+WERT1.CDPF;}
              Datensatz->Alfa=(WERT2.Alfa-WERT1.Alfa)/(Ma_runner->next->MA_NUM-Ma_runner->MA_NUM)*(Mach-Ma_runner->MA_NUM)+WERT1.Alfa;
              Datensatz->CD=Datensatz->CDP+Datensatz->CDI;
            }
          }
        }
        Ma_runner=Ma_runner->next;
      }while(Ma_runner!=NULL);*/
	WERT1=getProfilBeiwerteMA(EinP_runner, Mach, Datensatz->CA, Datensatz->RE_LOC, Datensatz->Flap_Set, input,EinP_runner->Profiel_Name);
    }
    if (strcmp (EinP_runner->Profiel_Name, Datensatz->SEC2_Name)==0)
    {
	WERT2=getProfilBeiwerteMA(EinP_runner, Mach, Datensatz->CA, Datensatz->RE_LOC, Datensatz->Flap_Set, input,EinP_runner->Profiel_Name);
    }
    EinP_runner=EinP_runner->next;
  }while(EinP_runner!=NULL);
  Datensatz->CDP=(WERT1.CDP*Datensatz->SEC_Percent/100+WERT2.CDP*(100-Datensatz->SEC_Percent)/100);
  if (input->useCm) {Datensatz->CM=(WERT1.Cm*Datensatz->SEC_Percent/100+WERT2.Cm*(100-Datensatz->SEC_Percent)/100);}
  if (input->useCdpw) {Datensatz->CDPW=(WERT1.CDPW*Datensatz->SEC_Percent/100+WERT2.CDPW*(100-Datensatz->SEC_Percent)/100);}
  if (input->useCdpp) {Datensatz->CDPP=(WERT1.CDPP*Datensatz->SEC_Percent/100+WERT2.CDPP*(100-Datensatz->SEC_Percent)/100);}
  if (input->useCdpf) {Datensatz->CDPF=(WERT1.CDPF*Datensatz->SEC_Percent/100+WERT2.CDPF*(100-Datensatz->SEC_Percent)/100);}
  Datensatz->Alfa=(WERT1.Alfa*Datensatz->SEC_Percent/100+WERT2.Alfa*(100-Datensatz->SEC_Percent)/100);
  Datensatz->CD=Datensatz->CDP+Datensatz->CDI;
  //Debug
  //printf("Raus:getProfilBeiwerte\n");
}

void PolarInterpolation(struct inputfile *input, struct PLTdistribution *PLTdist, struct file14 *rfile, struct GESAMT_POLARint *GPolare)
{
  int lauf,lauf1, lauf2,lauf3,lauf4, schalter;
  struct Profiele AlleProfile;
  struct EinProfiel *EinP_runner;
  struct EinProfiel *TEST_runner;
  struct Polar_Integration_Distribution Pdistri;
  char ExecName[250] = {"\0"};
  double Speed,ws;
  FILE *fout;
  
  AlleProfile.Anzahl_InputProfile=0;
  for (lauf1=0;lauf1<input->NbWings;lauf1++)
  {    
    for (lauf2=0; lauf2<input->Fluegel[lauf1].anzahlSchnitte; lauf2++)
    {
      schalter=0;
      if (AlleProfile.Anzahl_InputProfile==0)
      {
        AlleProfile.Profiel=(struct EinProfiel *)malloc(sizeof(struct EinProfiel));
        EinP_runner=AlleProfile.Profiel;
        EinP_runner->next=NULL;        
        readProfil (&input->Fluegel[lauf1].Schnitte[lauf2], input->PROFIPFAD,EinP_runner,input);
        AlleProfile.Anzahl_InputProfile=1;
      }else{
        EinP_runner=AlleProfile.Profiel;
        do
        {
          TEST_runner=EinP_runner;
          //printf ("%s %s \n",EinP_runner->Profiel_Name,input->Fluegel[lauf1].Schnitte[lauf2].ProfilName);          
          if (strcmp(EinP_runner->Profiel_Name,input->Fluegel[lauf1].Schnitte[lauf2].ProfilName)==0)
          {
            //printf("Schalt1\n");
            schalter=1;
          }
          if (EinP_runner->next!=NULL)
          {
            EinP_runner=EinP_runner->next;
          }
        }while ( TEST_runner->next != NULL );
        //printf ("Schalter: %d\n", schalter);
        if (schalter==0 && EinP_runner->next == NULL)
        {
          EinP_runner->next=(struct EinProfiel *)malloc(sizeof(struct EinProfiel));
          if (EinP_runner->next==NULL)
          {
            printf("Error emory alloc\n\n");
          }
          EinP_runner=EinP_runner->next;
          EinP_runner->next=NULL;
          AlleProfile.Anzahl_InputProfile++;
          readProfil (&input->Fluegel[lauf1].Schnitte[lauf2], input->PROFIPFAD,EinP_runner,input);
        }      
      }      
    }   
  } 
  printf("AnzProfile: %d\n",AlleProfile.Anzahl_InputProfile);
  Pdistri.AnzFaelle=PLTdist->AnzFaelle;
  Pdistri.AnzFluegel=input->NbWings;
  Pdistri.Fall=(struct EinFallP *)malloc(PLTdist->AnzFaelle*sizeof(struct EinFallP));
  for (lauf=0;lauf<PLTdist->AnzFaelle;lauf++)
  {
    Pdistri.Fall[lauf].Fluegel=(struct FLUEGELP *)malloc(Pdistri.AnzFluegel*sizeof(struct FLUEGELP));
    for (lauf2=0;lauf2<Pdistri.AnzFluegel;lauf2++)
    {
      //DEBUG
      //printf("\nFluegelNr.:%d\n",lauf2+1);
      Pdistri.Fall[lauf].Fluegel[lauf2].Distri=(struct VERTEILUNG_POLARint *)malloc(input->Fluegel[lauf2].AnzahlSpanPanel*sizeof(struct VERTEILUNG_POLARint));
      for (lauf3=0;lauf3<input->Fluegel[lauf2].AnzahlSpanPanel;lauf3++)
      {
        //DEBUG
        //printf("\nSpanPanelNr.:%d\n",lauf3+1);
        Pdistri.Fall[lauf].Fluegel[lauf2].Distri[lauf3].X=PLTdist->Fall[lauf].Fluegel[lauf2].Distri[lauf3].X;
        Pdistri.Fall[lauf].Fluegel[lauf2].Distri[lauf3].Y=PLTdist->Fall[lauf].Fluegel[lauf2].Distri[lauf3].Y;
        Pdistri.Fall[lauf].Fluegel[lauf2].Distri[lauf3].Z=PLTdist->Fall[lauf].Fluegel[lauf2].Distri[lauf3].Z;
        Pdistri.Fall[lauf].Fluegel[lauf2].Distri[lauf3].S=PLTdist->Fall[lauf].Fluegel[lauf2].Distri[lauf3].S;        
        Pdistri.Fall[lauf].Fluegel[lauf2].Distri[lauf3].CA=PLTdist->Fall[lauf].Fluegel[lauf2].Distri[lauf3].CFNORMAL;
        Pdistri.Fall[lauf].Fluegel[lauf2].Distri[lauf3].CDI=PLTdist->Fall[lauf].Fluegel[lauf2].Distri[lauf3].CFX_FROM_CDI;
        Pdistri.Fall[lauf].Fluegel[lauf2].Distri[lauf3].tiefe=PLTdist->Fall[lauf].Fluegel[lauf2].Distri[lauf3].REF_LEN_CMY;
        Pdistri.Fall[lauf].Fluegel[lauf2].Distri[lauf3].areaPan=PLTdist->Fall[lauf].Fluegel[lauf2].Distri[lauf3].REF_AREA;
        //DEBUG
        //printf("\nPosX:%lf PosY:%lf PosZ:%lf PosS:%lf Ca:%lf Cd:%lf cord:%lf area:%lf\n",Pdistri.Fall[lauf].Fluegel[lauf2].Distri[lauf3].X,Pdistri.Fall[lauf].Fluegel[lauf2].Distri[lauf3].Y,Pdistri.Fall[lauf].Fluegel[lauf2].Distri[lauf3].Z,Pdistri.Fall[lauf].Fluegel[lauf2].Distri[lauf3].S,Pdistri.Fall[lauf].Fluegel[lauf2].Distri[lauf3].CA,Pdistri.Fall[lauf].Fluegel[lauf2].Distri[lauf3].CDI,Pdistri.Fall[lauf].Fluegel[lauf2].Distri[lauf3].tiefe,Pdistri.Fall[lauf].Fluegel[lauf2].Distri[lauf3].areaPan);
        if (input->MASSE_SPEED==1)
        {
          if ((input->POLINT_MASSE<=0)||(input->RefDensity<=0)||(input->refAerea<=0))
          {
            errmessage (48);
          }
          if (input->numberAlpha>0)
          {
            errmessage (49);
          }         
          Speed=sqrt((input->POLINT_MASSE*9.81*2)/(input->RefDensity*input->refAerea*input->targetCl[lauf]));
          if (input->KreisflugModus==1)
          {
              Speed=sqrt(pow(Speed/sqrt(cos(input->HaengeWinkel* 3.141592653/180)),2)+ pow(input->Steigen,2));
              Speed=Speed*(1-Pdistri.Fall[lauf].Fluegel[lauf2].Distri[lauf3].Y*input->PhiKreis)*cos(input->beta* 3.141592653/180);
          }
        }else{
          //Speed=sqrt((input->POLINT_MASSE*9.81*2)/(input->RefDensity*input->refAerea*input->targetCl[lauf]));
          Speed=input->RefSpeed;
          if (input->KreisflugModus==1)
          {  
              Speed=sqrt(pow(Speed/sqrt(cos(input->HaengeWinkel* 3.141592653/180)),2)+ pow(input->Steigen,2));
              Speed=Speed*(1-Pdistri.Fall[lauf].Fluegel[lauf2].Distri[lauf3].Y*input->PhiKreis)*cos(input->beta* 3.141592653/180);
          }
        }
        if (Speed!=input->RefSpeed)
        {
          Pdistri.Fall[lauf].Fluegel[lauf2].Distri[lauf3].CA=Pdistri.Fall[lauf].Fluegel[lauf2].Distri[lauf3].CA*pow(input->RefSpeed,2)/pow(Speed,2);
        }
        
        //Debug
        //printf("Y-Pos: %lf  Speed: %lf  PhiKreis: %lf RefSpeed: %lf\n",Pdistri.Fall[lauf].Fluegel[lauf2].Distri[lauf3].Y,Speed,input->PhiKreis,input->RefSpeed);
        Pdistri.Fall[lauf].Fluegel[lauf2].Distri[lauf3].Vloc=Speed;
        Pdistri.Fall[lauf].Fluegel[lauf2].Distri[lauf3].RE_LOC=Speed*PLTdist->Fall[lauf].Fluegel[lauf2].Distri[lauf3].REF_LEN_CMY/MUE;
        schalter=0;
        lauf4=1;
        do
        {
          schalter=schalter+input->Fluegel[lauf2].Schnitte[input->Fluegel[lauf2].anzahlSchnitte-lauf4].AnzahlPan;
          if (schalter<(lauf3+1))
          {
            lauf4++;
          }   
        }while (schalter<(lauf3+1));
        strcpy(Pdistri.Fall[lauf].Fluegel[lauf2].Distri[lauf3].SEC_Name,input->Fluegel[lauf2].Schnitte[input->Fluegel[lauf2].anzahlSchnitte-lauf4].ProfilName);
        Pdistri.Fall[lauf].Fluegel[lauf2].Distri[lauf3].SEC_Percent=(((double)schalter+1)-((double)lauf3+1.0))/((double)input->Fluegel[lauf2].Schnitte[input->Fluegel[lauf2].anzahlSchnitte-lauf4].AnzahlPan+1.0)*100;
        strcpy(Pdistri.Fall[lauf].Fluegel[lauf2].Distri[lauf3].SEC2_Name,input->Fluegel[lauf2].Schnitte[input->Fluegel[lauf2].anzahlSchnitte-lauf4-1].ProfilName);
        if (input->Fluegel[lauf2].Schnitte[input->Fluegel[lauf2].anzahlSchnitte-lauf4].KlappeNrI>0)
        {
          Pdistri.Fall[lauf].Fluegel[lauf2].Distri[lauf3].Flap_Set=input->Klappen[input->Fluegel[lauf2].Schnitte[input->Fluegel[lauf2].anzahlSchnitte-lauf4].KlappeNrI-1].winkel;
        }else{
          Pdistri.Fall[lauf].Fluegel[lauf2].Distri[lauf3].Flap_Set=0.0;
        }
        getProfilBeiwerte(&Pdistri.Fall[lauf].Fluegel[lauf2].Distri[lauf3],AlleProfile,0.0, input);
      }
    }
  }
  if (CheckWarn==2)
  {
    printf("Programm stopt aufgrund  einer Interpolation-Warning!\n");
    errmessage(1);
  }

  printf ("\nInterpolation durchgefuehrt!\n");
    
  fout=fopen("PolarInterpolation_Distribution.plt","w+");
  fprintf(fout,"TITLE = \"Polar-Interpolation-Coefficient-Distribution\" \n");
  fprintf(fout,"VARIABLES = X, Y, Z, S, Alfa, CL, CD, CDI, CDP");
  if (input->useCm) {fprintf(fout,", CM");}
  if (input->useCdpw) {fprintf(fout,", CDPW");}
  if (input->useCdpp) {fprintf(fout,", CDPP");}
  if (input->useCdpf) {fprintf(fout,", CDPF");}
  fprintf(fout,", tiefe, RE_LOC, V_LOC, Flap-Set\n");
  for (lauf=0;lauf<PLTdist->AnzFaelle;lauf++)
  {
    GPolare[lauf].CDP=0;
    GPolare[lauf].AERA=0;
    for (lauf2=0;lauf2<Pdistri.AnzFluegel;lauf2++)
    {
      fprintf(fout,"ZONE T=\"FALL %d , Fluegel %d\" \n", (lauf+1), (lauf2+1));
      for (lauf3=0;lauf3<input->Fluegel[lauf2].AnzahlSpanPanel;lauf3++)
      {
        //Polare verteilung plotten        
        fprintf(fout, "%lf %lf %lf %lf %lf %lf %lf %lf %lf", Pdistri.Fall[lauf].Fluegel[lauf2].Distri[lauf3].X, Pdistri.Fall[lauf].Fluegel[lauf2].Distri[lauf3].Y, Pdistri.Fall[lauf].Fluegel[lauf2].Distri[lauf3].Z, Pdistri.Fall[lauf].Fluegel[lauf2].Distri[lauf3].S, Pdistri.Fall[lauf].Fluegel[lauf2].Distri[lauf3].Alfa, Pdistri.Fall[lauf].Fluegel[lauf2].Distri[lauf3].CA, Pdistri.Fall[lauf].Fluegel[lauf2].Distri[lauf3].CD, Pdistri.Fall[lauf].Fluegel[lauf2].Distri[lauf3].CDI,Pdistri.Fall[lauf].Fluegel[lauf2].Distri[lauf3].CDP);
        if (input->useCm) {fprintf(fout, " %lf", Pdistri.Fall[lauf].Fluegel[lauf2].Distri[lauf3].CM);}                 
        if (input->useCdpw) {fprintf(fout, " %lf", Pdistri.Fall[lauf].Fluegel[lauf2].Distri[lauf3].CDPW);}
        if (input->useCdpp) {fprintf(fout, " %lf", Pdistri.Fall[lauf].Fluegel[lauf2].Distri[lauf3].CDPP);}
        if (input->useCdpf) {fprintf(fout, " %lf", Pdistri.Fall[lauf].Fluegel[lauf2].Distri[lauf3].CDPF);}        
        fprintf(fout, " %lf %lf %lf %lf",Pdistri.Fall[lauf].Fluegel[lauf2].Distri[lauf3].tiefe ,Pdistri.Fall[lauf].Fluegel[lauf2].Distri[lauf3].RE_LOC,Pdistri.Fall[lauf].Fluegel[lauf2].Distri[lauf3].Vloc,Pdistri.Fall[lauf].Fluegel[lauf2].Distri[lauf3].Flap_Set);
        
        //DEBUG
        if (input->debuglevel==2 || input->debuglevel>99)
        {
          fprintf(fout, " %s %s %lf\n",Pdistri.Fall[lauf].Fluegel[lauf2].Distri[lauf3].SEC_Name,Pdistri.Fall[lauf].Fluegel[lauf2].Distri[lauf3].SEC2_Name,Pdistri.Fall[lauf].Fluegel[lauf2].Distri[lauf3].SEC_Percent);
        }
        fprintf(fout, "\n");
        
        // Rueckrechnung der lokalen Beiwerte fuer die globale Integration
        if (Speed!=input->RefSpeed)
        {
          Pdistri.Fall[lauf].Fluegel[lauf2].Distri[lauf3].CA=Pdistri.Fall[lauf].Fluegel[lauf2].Distri[lauf3].CA*pow(Speed,2)/pow(input->RefSpeed,2);
          Pdistri.Fall[lauf].Fluegel[lauf2].Distri[lauf3].CD=Pdistri.Fall[lauf].Fluegel[lauf2].Distri[lauf3].CD*pow(Speed,2)/pow(input->RefSpeed,2);
          if (input->useCm) {Pdistri.Fall[lauf].Fluegel[lauf2].Distri[lauf3].CM=Pdistri.Fall[lauf].Fluegel[lauf2].Distri[lauf3].CM*pow(Speed,2)/pow(input->RefSpeed,2);}
          if (input->useCdpw) {Pdistri.Fall[lauf].Fluegel[lauf2].Distri[lauf3].CDPW=Pdistri.Fall[lauf].Fluegel[lauf2].Distri[lauf3].CDPW*pow(Speed,2)/pow(input->RefSpeed,2);}
          if (input->useCdpp) {Pdistri.Fall[lauf].Fluegel[lauf2].Distri[lauf3].CDPP=Pdistri.Fall[lauf].Fluegel[lauf2].Distri[lauf3].CDPP*pow(Speed,2)/pow(input->RefSpeed,2);}
          if (input->useCdpf) {Pdistri.Fall[lauf].Fluegel[lauf2].Distri[lauf3].CDPF=Pdistri.Fall[lauf].Fluegel[lauf2].Distri[lauf3].CDPF*pow(Speed,2)/pow(input->RefSpeed,2);}
        }
        
        //Integration
        GPolare[lauf].CDP=GPolare[lauf].CDP+Pdistri.Fall[lauf].Fluegel[lauf2].Distri[lauf3].CDP*Pdistri.Fall[lauf].Fluegel[lauf2].Distri[lauf3].areaPan;
        GPolare[lauf].AERA=GPolare[lauf].AERA+Pdistri.Fall[lauf].Fluegel[lauf2].Distri[lauf3].areaPan;
        if (input->useCdpw) {GPolare[lauf].CDPW=GPolare[lauf].CDPW+Pdistri.Fall[lauf].Fluegel[lauf2].Distri[lauf3].CDPW*Pdistri.Fall[lauf].Fluegel[lauf2].Distri[lauf3].areaPan;}
        if (input->useCdpp) {GPolare[lauf].CDPP=GPolare[lauf].CDPP+Pdistri.Fall[lauf].Fluegel[lauf2].Distri[lauf3].CDPP*Pdistri.Fall[lauf].Fluegel[lauf2].Distri[lauf3].areaPan;}
        if (input->useCdpf) {GPolare[lauf].CDPF=GPolare[lauf].CDPF+Pdistri.Fall[lauf].Fluegel[lauf2].Distri[lauf3].CDPF*Pdistri.Fall[lauf].Fluegel[lauf2].Distri[lauf3].areaPan;}
      }
    }
    if (input->GlobSym==1)
    {
      GPolare[lauf].CDP=GPolare[lauf].CDP/(input->refAerea/2);
      if (input->useCdpw) {GPolare[lauf].CDPW=GPolare[lauf].CDPW/(input->refAerea/2);}
      if (input->useCdpp) {GPolare[lauf].CDPP=GPolare[lauf].CDPP/(input->refAerea/2);}
      if (input->useCdpf) {GPolare[lauf].CDPF=GPolare[lauf].CDPF/(input->refAerea/2);}
    }else{
      GPolare[lauf].CDP=GPolare[lauf].CDP/input->refAerea;
      if (input->useCdpw) {GPolare[lauf].CDPW=GPolare[lauf].CDPW/input->refAerea;}
      if (input->useCdpp) {GPolare[lauf].CDPP=GPolare[lauf].CDPP/input->refAerea;}
      if (input->useCdpf) {GPolare[lauf].CDPF=GPolare[lauf].CDPF/input->refAerea;}
    }
    printf ("Area Integriert & Input: %lf %lf\n",GPolare[lauf].AERA, input->refAerea );
  }
  
  fclose(fout);
  strcpy(ExecName,"mv PolarInterpolation_Distribution.plt LaufV1.lili.");
  strcat(ExecName,input->LiLiVersionName); 
  system(ExecName);

  fout=fopen("PolarInterpolation_Coefficents.plt","w+");
  fprintf(fout,"TITLE = \"Polar-Interpolation-Coefficient\" \n");
  fprintf(fout,"VARIABLES = Alfa, CL, CD, CDI, CDP");
  if (input->useCdpw) {fprintf(fout,", CDPW");}
  if (input->useCdpp) {fprintf(fout,", CDPP");}
  if (input->useCdpf) {fprintf(fout,", CDPF");}
  if (input->MASSE_SPEED) {fprintf(fout,", V, E, ws");}
  fprintf(fout,"\n");
  fprintf(fout,"ZONE T=\"polare\" \n");
  for (lauf=0;lauf<PLTdist->AnzFaelle;lauf++)
  {
    fprintf(fout, "%lf %lf %lf %lf %lf", rfile->Fall[lauf].Alfa, rfile->Fall[lauf].CA, rfile->Fall[lauf].CWI+GPolare[lauf].CDP+input->BASIS_CD_POL, rfile->Fall[lauf].CWI, GPolare[lauf].CDP);
    if (input->useCdpw) {fprintf(fout, " %lf", GPolare[lauf].CDPW);}
    if (input->useCdpp) {fprintf(fout, " %lf", GPolare[lauf].CDPP);}
    if (input->useCdpf) {fprintf(fout, " %lf", GPolare[lauf].CDPF);}
    if (input->MASSE_SPEED) 
    {
      Speed=sqrt((input->POLINT_MASSE*9.81*2)/(input->RefDensity*input->refAerea*rfile->Fall[lauf].CA));
      ws=Speed*(rfile->Fall[lauf].CWI+GPolare[lauf].CDP+input->BASIS_CD_POL)/rfile->Fall[lauf].CA;
      if (input->KreisflugModus==1)
      {
          Speed=sqrt(pow(Speed/sqrt(cos(input->HaengeWinkel* 3.141592653/180)),2)+ pow(input->Steigen,2));
          ws=sqrt((input->POLINT_MASSE*9.81*2)/(input->RefDensity*input->refAerea))*(rfile->Fall[lauf].CWI+GPolare[lauf].CDP+input->BASIS_CD_POL)/(pow(cos(input->HaengeWinkel* 3.141592653/180),1.5)*pow(rfile->Fall[lauf].CA,1.5));             
      }
      fprintf(fout," %lf %lf %lf",Speed*3.6,rfile->Fall[lauf].CA/(rfile->Fall[lauf].CWI+GPolare[lauf].CDP+input->BASIS_CD_POL),ws );
    }
    fprintf(fout, "\n");
  }
  fclose(fout);

  strcpy(ExecName,"mv PolarInterpolation_Coefficents.plt LaufV1.lili.");
  strcat(ExecName,input->LiLiVersionName); 
  system(ExecName);

  //print_screan_debug_Polaren(AlleProfile, input);
  /*EinP_runner=AlleProfile.Profiel;
  do
  {
    printf("ProfileName: %s\n",EinP_runner->Profiel_Name);
    EinP_runner=EinP_runner->next;
  }while (EinP_runner!=NULL);*/
}
