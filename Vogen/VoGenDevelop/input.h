/********************input.h**********************                  
*                                                *      
*      Vorentwurfs-Generator "VoGen"             *          
*      Version 2.00                              *
*      Date of last modification 20.04.2020      *
*                                                *          
*      Writen by: Jan Himisch                    *      
*             jan.himisch@dlr.de                 *      
*                                                *      
*                                                *
*************************************************/
#ifndef _INPUT_H_
#define _INPUT_H_


#define MAXCAMB 1000
#define MaxSurfixLength 250
#define MaxPrefixLength 10
#define MaxZeilenLaenge 800
#define MAX_LENGTH_PNAME 150
#define T_BODEN 288,15
#define P_BODEN 101325
#define R_GASConst 287.05287
#define GAMMA -6.5
#define PI 3.141592653
#define MUE 0.0000149
#define MaxRegelLength 250

enum{ OFF, ON };

struct EinSchnitt 
{
  double posx;
  double posy;
  double posz;
  double AlphaPlus;
  double OrgTwist; // Dieser Wert wird nur gebraucht um bei parametrischer Ausf�hrung den inizialen Twist zu behalten!
  double Vloc;
  double tiefe;
  double relDicke;
  int AnzahlPan;
  int AnzahlXY;
  char ProfilName[50];
  double xcamb[MAXCAMB];
  double zcamb[MAXCAMB];
  int KlappeNrI;
  int KlappeNrA;
  double etaSchnitt;
  int anzPanelKlappe;
};


struct Kreuzung
{
  int KreuzFl1;
  int KreuzFl2;
  int KreuzPosFl1;
  int KreuzPosFl2;
  int KreuzTiefFl1;
  int KreuzTiefFl2;
  int GabBedFl1;
  int GabBedFl2;
  int KopelArt;
};

struct Einfluegel
{
  int anzahlSchnitte;
  int AnzahlTiefe;
  int SymExt;
  int SymRandBed;
  int AnzahlSpanPanel;
  int Spiegeln;
  char Surfix[MaxSurfixLength];
  char Prefix[MaxPrefixLength];
  struct EinSchnitt *Schnitte;
};

struct VariationsPunkt 
{
  int Fluegel;
  int Schnitt;
  unsigned twist:1;
  unsigned X:1;
  unsigned Y:1;
  unsigned Z:1;
  unsigned tiefe:1;
  unsigned Vloc:1;
};

struct SplineVariationsPunkt 
{
  int AnzahlStuetzpunkte;
  int Fluegel;
  int StartSchnitt;
  int EndSchnitt;
  struct VariationsPunkt Variable;
};

struct BSURF_LISTZONE
{
  int *Zone;
};

struct Struct_Dihedral_Dist
{
  double etaPos;
  double etaV;
  double Syz;
};

struct EineKlappen
{
  int Fluegel;
  int Schnitt1;
  int Schnitt2;
  int AnzPanel;
  int Spiegeln;
  double eta1;
  double eta2;
  double winkel;
};

struct EinMischer
{
  char Regel[MaxRegelLength];
};

struct PARAvariation
{
  int parametric;                          
  int NrVariations;                         
  char PARAINP[100];                        
  int anzVar;
  int Paraset4eachCLA;
  int AbsOrDiffer;
  int anzSpVar;
  int link;                      
  struct VariationsPunkt *VarParameter;
  struct SplineVariationsPunkt *SplineVarPunkt;
  //info
  int AnzahlParameter;
};

struct inputfile 
{ 
  //debug
  int debuglevel;
  //HEADER Lifting-Line
  double MachNumber; 
  double beta; 
  double refSpan;
  double refAerea;
  double refChordlength;                        
  double origXYZ[3];                           
                                                              
  double Alfa_Rot;
                          
  //Lifting-Line-Steuerung                      
  int numberCL;                         
  double *targetCl;                         
  int numberAlpha;                          
  double *targetAlpha;                        
  int LiliStart;
  int runOneByOne; 
  char LiLiVersion[25]; //wird vom Nutzer gesetzt
  char LiLiVersionName[25]; //wird durchs Programm in Abhängigkeit von LiLiVersion gesetzt                       
  char LILIEXE[50];  
  char LAUFNAME[50];          
  int XMLEA;                               
  int OldInput;

  //FreeWake Steuerung
  char FwVersion[25];
  char FwEXE[50];
  int FwStart;
  int FwSteadyUnstedy;
  int FwRelaxWake;
  int FwnumbTimeStep;
  double FwtimeStepLength;
  double FwConvergDelta;
                                    
  //Geometrische Fl�gelbedingungen                  
  int GlobSym;
  double scalfac; 
  int NbWings;                          
  int AnzKreuz;                             
  int readRelDick;
  int CamberMethod;
  double CamberStart;
  double CamberEnd;                            
  double AlfaPlusRumpf;
  struct Einfluegel *Fluegel;
  struct Kreuzung *Kreuz;     
  
  //KlappenSettings                  
  int AnzahlKlappen;
  struct EineKlappen *Klappen;

  //Mischer                 
  int AnzahlMischer;
  int AnzahlMixVar;
  double *MixVar;
  struct EinMischer *Mischer;
                                    
  //Parametric-Variations 
  int NumberOfParametSetsinput;                   
  struct PARAvariation *ParaVar;                               
 
  //Trimmen                            
  int Austrimmen;
  int ArotAdjust;                         
  int MaxCmIter;                         
  int CmSteuerFl;
  int CmSteuerKlappe;                             
  char  EpsAusgabe[100];
  double ZielCm;
  double InitStepSize;
  
  //Zus�tzliche Aufpunkte?? 
//  int zusAufP; 
  
  //Post-Processing
  int IBM_BD_WING;
  int IBM_BD_SEC;
  int WRBM_BD_WING;
  int WRBM_BD_SEC;
  double RefSpeed;
  double RefDensity;
  double RefTemp;
  double BasisIntBiegeMoment;
  double BasisIntDickenBiegeMoment;

  //BSURF - Input - File 
  int makeBsurfinp;
  int NUM_OF_ZONES;
  int TAUbsurfINP;
  struct BSURF_LISTZONE *BSURF_FLUEGEL;
  int BsurfSquareDiff;
  
  //Die Methode
  int kHOrd;
  double *kH;             
  double BCwi;             
  double BWRBM;             
  double BCw0;
  double UrArea; 
  double BAlafa;
  
  //Struct-Kopplung
  int Struckt_Coup_ON;
  int Load_Case_ON;
  int Numb_Iter;
  int Log_Iter_ON;
  int useJIG2FLIGHT;
  /*int Span_Cor_ON; Sollte mit der neusten Version von Fame nicht mehr n�tig sein*/
  char ProjektsionsSchnitt[6];
  int ProjektionsFluegel;
  int ProjektionsSchnitt;
  char StructProg_EXEdir[256];
  char StructProg_EXEFile[50];
  char StructProg_EXEFile2[50];
  char StructProg_ResultSubDir[256];
  char LILIDIR[256];
  char StructProg_ResultFile[50];
  char StructProg_Jig2FlightFile[50];
  double Target_CL_LoadCase;
  double Ref_Density_LoadCase;
  double Ref_Speed_LoadCase;
  double deltaStructX;  
  double deltaStructZ;
  double StructModBezSpan;
  double scalefactStruct;
  int NumStructDih;
  struct Struct_Dihedral_Dist *StrDihDist; 
  
  //Kreisflug-Optionen
  int KreisflugModus;
  double HaengeWinkel;
  double Steigen;
  double PhiKreis;
   
  //XML-KOPLUNG
  int XMLKOP;
  int PolintStart;
  int InterneInterpolation; 
  int MASSE_SPEED;
  int useCdpw;
  int useCdpp;
  int useCdpf;
  int useCm;
  int ClMinMaxHandling;
  double PENALTYfact;
  char POLINTEXE[MAX_LENGTH_PNAME];
  char PROFIPFAD[200];
  double BASIS_CD_POL;
  char BASIS_CD_POL_File[100];
  double POLINT_MASSE;
  
    
  //INFO                                            
  int AnzahlPanel;
  int AnzahlSpanPanel;
};

void read_input(char *datei,struct inputfile *input);
void init_input(struct inputfile *input);
void init_Einfluegel(struct Einfluegel *Fluegel);
void get_Vloc(struct inputfile *input);

#endif
