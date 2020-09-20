/******************run_FW.c***********************
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
#include "run_FW.h"
#include "errmes.h"

const char HEADFW1[]= "FreeWake-Input created by VoGen\n";
const char HEADFW2[]= "xleft    yleft   zleft   chord   epsilon   Bound.Cond.   Airfoil   Hinge   Deflection\n";
const char HEADFW3[]= "xright   yright  zright  chord   epsilon   Bound.Cond.   Airfoil   Hinge   Deflection\n";
const char HEADFW5[]= "%<- special identifier\n";
const char HEADFW6[]= "Vertical tail information:\n";
const char HEADFW7[]= "Number of panels (max 5) = 0\n";
const char HEADFW8[]= "no.        chord   area    airfoil\n";
const char HEADFW9[]= "Fuselage information:\n";
const char HEADFW10[]= "Number of sections (max 20) =     0\n";
const char HEADFW11[]= "Width of each section =           0.000\n";
const char HEADFW12[]= "Panel where transition occurs =   0\n";
const char HEADFW13[]= "No.       Diamter\n";
const char HEADFW14[]= "Interference drag = 0.0%\n";


int IsKoppFW(int FLNr,int SchnNR,struct inputfile *input)
{
   int lauf, lauf2,lauf3, returnVal;
   returnVal=0;
   if (input->AnzKreuz!=0)
   {
       for (lauf=0; input->AnzKreuz; lauf++)
       {
           if (input->Kreuz[lauf].KreuzFl1==FLNr && input->Kreuz[lauf].KreuzPosFl1==SchnNR)
           {
               returnVal++;
               for (lauf2=1;lauf2<=FLNr;lauf2++)
               {
                   for(lauf3=1;lauf3<=input->Fluegel[lauf2-1].anzahlSchnitte;lauf3++)
                   {
                       if (lauf2==FLNr)
                       {
                           if (lauf2<input->Kreuz[lauf].KreuzPosFl2)
                           {
                               returnVal++;
                           }
                       }else{
                           returnVal++;
                       }
                   }
               }
           }
           if (input->Kreuz[lauf].KreuzFl2==FLNr && input->Kreuz[lauf].KreuzPosFl2==SchnNR)
           {
               returnVal++;
               for (lauf2=1;lauf2<=FLNr;lauf2++)
               {
                   for(lauf3=1;lauf3<=input->Fluegel[lauf2-1].anzahlSchnitte;lauf3++)
                   {
                        if (lauf2==FLNr)
                        {
                           if (lauf2<input->Kreuz[lauf].KreuzPosFl1)
                           {
                              returnVal++;
                           }
                        }else{
                           returnVal++;
                        }
                   }
               }
           }
       }
   }

   return returnVal;
}

int whichAirfoil(struct Einfluegel *Fluegel,int laufSchnitt, struct ProfileSortier *Airfoils,int anzAirfoils)
{
   int lauf, returnval;
   returnval=0;
   for (lauf=0; lauf<anzAirfoils; lauf++)
   {
       if (strcmp(Fluegel->Schnitte[laufSchnitt].ProfilName, Airfoils[lauf].Airfoil_Str)==0)
       {
           returnval=lauf+1;
           break;
       }
   }
   return returnval;
}

void isFlap(double *eta, double *FlapAng, struct EinSchnitt Schnitt,int SchnittNr,struct EineKlappen *Klappen, int IA)
{
   int FlapNr;
   if (IA==1)
   {
       FlapNr=Schnitt.KlappeNrA;
   }else{
       FlapNr=Schnitt.KlappeNrI;
   }
   *eta=0.0;
   *FlapAng=0.0;
   if (FlapNr!=0)
   {

       if (Klappen[FlapNr-1].Schnitt1==(SchnittNr+1))
       {
           *eta=Klappen[FlapNr-1].eta1;
           *FlapAng=Klappen[FlapNr-1].winkel;
       }else if (Klappen[FlapNr-1].Schnitt2==(SchnittNr+1))
       {
           *eta=Klappen[FlapNr-1].eta2;
           *FlapAng=Klappen[FlapNr-1].winkel;
       }else if ((SchnittNr+1)>=Klappen[FlapNr-1].Schnitt1 && (SchnittNr+1)<=Klappen[FlapNr-1].Schnitt2 )
       {
           *eta=Schnitt.etaSchnitt;
           *FlapAng=Klappen[FlapNr-1].winkel;
       }else{
           printf ("Der Sinn des Lebens ist 42! Und Tschüss\n");
           printf ("SchnitNr:%d Schnitt1:%d Schnitt2:%d\n",SchnittNr+1,Klappen[FlapNr-1].Schnitt1,Klappen[FlapNr-1].Schnitt2);
           errmessage(1);
       }
   }

}
/*void checkDir()
{

}*/

void copyCamber(struct inputfile *input, struct ProfileSortier *Airfoils, int Anzahl)
{
    int lauf,lauf2,lauf3;
    char ExecName[150] = {"\0"};
    char FileName[50] = {"\0"};
    char laufNR[4];

    FILE *fchamber;
    strcpy(ExecName,"mkdir camber");
    system(ExecName);

    for (lauf=1; lauf<=Anzahl; lauf++)
    {
        strcpy(FileName,"camber");
        if (lauf>=10)
        {
           laufNR[0]=48+(int)(lauf/10);
           laufNR[1]=48+lauf-((int)(lauf/10)*10);
           laufNR[2]='\0';
        } else {
           laufNR[0]=48+lauf;
           laufNR[1]='\0';
        }
        strcat(FileName, laufNR);
        strcat(FileName,".camb");
        fchamber = fopen(FileName,"w+");
        fprintf(fchamber,"%s generated by VoGen\n",FileName);
        for (lauf2=0;lauf2<=input->Fluegel[Airfoils[lauf].Fluegel].Schnitte[Airfoils[lauf].section].AnzahlXY;lauf2++)
        {
            fprintf(fchamber," %lf  %lf\n",input->Fluegel[Airfoils[lauf].Fluegel].Schnitte[Airfoils[lauf].section].xcamb[lauf2],input->Fluegel[Airfoils[lauf].Fluegel].Schnitte[Airfoils[lauf].section].zcamb[lauf2]);
        }

        fclose(fchamber);
        strcpy(ExecName,"mv camber");
        strcat(ExecName, laufNR);
        strcat(ExecName,".camb camber");
        system(ExecName);
        for (lauf3=0; lauf3<150; lauf3++)
        {
            ExecName[lauf3]="\0";
        }
        for (lauf3=0; lauf3<50; lauf3++)
        {
            FileName[lauf3]="\0";
        }
    }

}



void create_FW_File (struct inputfile *input)
{
   FILE *fFWinp;
   int lauf1,lauf2,lauf3,check1,anzPanel,anzAirfoils,KoopelPanel;
   double FlapEta, FlapAngle;
   struct ProfileSortier Airfoils[15];

   fFWinp = fopen("FW_Input.txt","w+");
   if (fFWinp == NULL)
   {
       errmessage (70);
   }

   /* Header Part */
   fputs(HEADFW1, fFWinp);
   fprintf(fFWinp,"Relaxed wake (yes 1, no 0):                           relax =    %d\n",input->FwRelaxWake);
   fprintf(fFWinp,"Steady (1) or unsteady (2):                    aerodynamics =    %d\n",input->FwSteadyUnstedy);
   fprintf(fFWinp,"Viscous solutions (1) or inviscid (0)  viscous              =     0\n");
   fprintf(fFWinp,"Using camber (1) or no camber (0)  camber                   =     1\n");
   fprintf(fFWinp,"Symmetrical geometry (yes 1, no 0):                     sym =    %d\n",input->GlobSym);
   fprintf(fFWinp,"Pitch trim (yes 1, no 0):                         trimPITCH =     0\n");
   fprintf(fFWinp,"Roll trim (yes 1, no 0):                           trimROLL =     0\n");
   fprintf(fFWinp,"Lift trim (yes 1, no 0):                             trimCL =     0\n");
   //fprintf(fFWinp,"Longitudinal trim (yes 1, no 0):                       trim =     0\n\n");
   fprintf(fFWinp,"Max. number of time steps:                          maxtime =    %d\n",input->FwnumbTimeStep);
   fprintf(fFWinp,"Width of each time step (sec):                      deltime =    %lf\n",input->FwtimeStepLength);
   fprintf(fFWinp,"Convergence delta-span effic.:                       deltae =    %lf     (0 if only timestepping)\n\n",input->FwConvergDelta);
   if (input->RefSpeed != -1)
   {
        fprintf(fFWinp,"Freestream velocity (leave value 1):                   Uinf =    %lf\n",input->RefSpeed);
   }else{
        fprintf(fFWinp,"Freestream velocity (leave value 1):                   Uinf =    -1\n");
   }
   fprintf(fFWinp,"AOA beginning, end, step size [deg]:                  alpha =    %lf %lf 1.0\n", input->targetAlpha[0], input->targetAlpha[0]);
   fprintf(fFWinp,"Sideslip angle [deg]:                                  beta =    %lf\n\n", input->beta);
   fprintf(fFWinp,"Circling flight information\n");
   fprintf(fFWinp,"Circling flight on (1) off (0):                    circling =    %d\n",input->KreisflugModus);
   fprintf(fFWinp,"Turning in horizontal plane (1)                  horizontal =    0    if (0), a/c descends\n");
   fprintf(fFWinp,"Bank angle (deg)                                        phi =    %lf\n",input->PhiKreis);
   fprintf(fFWinp,"Upwind velocity (m/s)                                    Ws =    %lf\n",input->Steigen);
   fprintf(fFWinp,"Velocity gradient (1/s)                            gradient =    0 .00000000001  .3864  .03141592654  \n\n");
   fprintf(fFWinp,"Density:                                            density =    %lf\n",input->RefDensity);
   fprintf(fFWinp,"Kinematic viscosity:                                     nu =    1.420000e-05\n\n");
   fprintf(fFWinp,"Reference area:                                           S =    %lf\n", input->refAerea);
   fprintf(fFWinp,"Reference span:                                           b =    %lf\n\n", input->refSpan);
   fprintf(fFWinp,"Mean aerodynamic chord:                                cmac =    %lf\n",input->refChordlength);
   fprintf(fFWinp,"Aircraft weight (N):                                      W =    %lf\n",input->POLINT_MASSE);
   fprintf(fFWinp,"CG location (x y z):                                     cg =    %lf %lf %lf\n", input->origXYZ[0] , input->origXYZ[1] , input->origXYZ[2] );
   fprintf(fFWinp,"CMo of wing:                                            CMo =    0.000000\n\n");
   fprintf(fFWinp,"No. of wings (max. 5):                                wings =    %d\n",input->NbWings);
   anzPanel=0;
   for (lauf1=0; lauf1< input->NbWings; lauf1++)
   {
       for(lauf2=0; lauf2< input->Fluegel[lauf1].anzahlSchnitte; lauf2++)
       {
           if (lauf2==0)
           {
              if (input->Fluegel[lauf1].SymExt == 1)
              {
                  anzPanel++;
              }
           }else{
               anzPanel++;
           }
       }
   }
   fprintf(fFWinp,"No. of panels:                                       panels =    %d\n",anzPanel);
   anzAirfoils=0;

   for (lauf1=0; lauf1< input->NbWings; lauf1++)
   {
       for(lauf2=0; lauf2< input->Fluegel[lauf1].anzahlSchnitte; lauf2++)
       {
           if (lauf2==0 && lauf1==0)
           {
               anzAirfoils++;
               strcpy(&Airfoils[0].Airfoil_Str[0],input->Fluegel[lauf1].Schnitte[lauf2].ProfilName);
           }else{
               check1=0;
               for (lauf3=0;lauf3<anzAirfoils;lauf3++)
               {
                   if (strcmp(input->Fluegel[lauf1].Schnitte[lauf2].ProfilName,&Airfoils[lauf3].Airfoil_Str[0])==0)
                   {
                       check1=1;
                   }
               }
               if (check1==0)
               {
                   strcpy(&Airfoils[anzAirfoils].Airfoil_Str[0],input->Fluegel[lauf1].Schnitte[lauf2].ProfilName);
                   anzAirfoils++;
               }
           }
       }
    }
   //printf ("Anzahl Airfoils: %d\n",anzAirfoils);
   if (anzAirfoils>15)
   {
       printf("No. of airfoils (max. 15):                  airfoils = %d\n",anzAirfoils);
       errmessage (71);
   }
   fprintf(fFWinp,"No. of airfoils (max. 15):                         airfoils =    %d\n\n\n",anzAirfoils);
   fprintf(fFWinp,"Panel boundary conditions:\n");
   fprintf(fFWinp,"        Symmetry line -         10\n");
   fprintf(fFWinp,"        Between panels -        220\n");
   fprintf(fFWinp,"        Free end -              100\n\n");
   fprintf(fFWinp,"Hinge is in percentage chord:\n");
   fprintf(fFWinp,"        Midchord -      0.5\n");
   fprintf(fFWinp,"        Leading edge -  0\n");
   fprintf(fFWinp,"        NO HINGE -      1\n");
   fprintf(fFWinp,"Defines leading edge of wing, all measured in metres:\n\n\n");


   anzPanel=1;
   for (lauf1=0; lauf1< input->NbWings; lauf1++)
   {
       for(lauf2=0; lauf2< input->Fluegel[lauf1].anzahlSchnitte; lauf2++)
       {
           if (lauf2==0)
           {
              if (input->Fluegel[lauf1].SymExt == 1)
              {
                  fprintf(fFWinp,"Panel #: %d. Number of elements across span n = %d and chord m = %d\n",anzPanel, input->Fluegel[lauf1].Schnitte[lauf2].AnzahlPan, input->Fluegel[lauf1].AnzahlTiefe);
                  fprintf(fFWinp,"Neighbouring panels (0 for none) left: 0 right: %d\n",anzPanel+1);
                  anzPanel++;
                  fputs(HEADFW2, fFWinp);
                  fprintf(fFWinp,"%.3lf   0.000   %.3lf   %.3lf   0.000    ",input->Fluegel[lauf1].Schnitte[lauf2].posx,input->Fluegel[lauf1].Schnitte[lauf2].posz,input->Fluegel[lauf1].Schnitte[lauf2].tiefe);
                  if (input->GlobSym == 1)
                  {
                      fprintf(fFWinp," 010    ");
                  }else{
                      fprintf(fFWinp," 100    ");
                  }
                  fprintf(fFWinp,"     %d    0    0\n",whichAirfoil(&input->Fluegel[lauf1],lauf2,Airfoils,anzAirfoils));
                  fputs(HEADFW3, fFWinp);
                  fprintf(fFWinp,"%.3lf   %.3lf   %.3lf   %.3lf   0.000    ",input->Fluegel[lauf1].Schnitte[lauf2].posx,input->Fluegel[lauf1].Schnitte[lauf2].posy,input->Fluegel[lauf1].Schnitte[lauf2].posz,input->Fluegel[lauf1].Schnitte[lauf2].tiefe);
                  fprintf(fFWinp," 220    ");
                  fprintf(fFWinp,"     %d    0    0\n\n",whichAirfoil(&input->Fluegel[lauf1],lauf2,Airfoils,anzAirfoils));
              }
           }else if (lauf2==input->Fluegel[lauf1].anzahlSchnitte-1){
               fprintf(fFWinp,"Panel #: %d. Number of elements across span n = %d and chord m = %d\n",anzPanel, input->Fluegel[lauf1].Schnitte[lauf2].AnzahlPan, input->Fluegel[lauf1].AnzahlTiefe);
               fprintf(fFWinp,"Neighbouring panels (0 for none) left: %d right: ",anzPanel-1);

               KoopelPanel=IsKoppFW(lauf1+1,lauf2+1,input);
               fprintf(fFWinp,"%d\n",KoopelPanel);

               anzPanel++;
               fputs(HEADFW2, fFWinp);
               fprintf(fFWinp,"%.3lf   %.3lf   %.3lf   %.3lf   0.000    ",input->Fluegel[lauf1].Schnitte[lauf2-1].posx,input->Fluegel[lauf1].Schnitte[lauf2-1].posy ,input->Fluegel[lauf1].Schnitte[lauf2-1].posz,input->Fluegel[lauf1].Schnitte[lauf2-1].tiefe);

               fprintf(fFWinp," 220    ");

               fprintf(fFWinp,"    %d    ",whichAirfoil(&input->Fluegel[lauf1],lauf2,Airfoils,anzAirfoils));
               if (input->Fluegel[lauf1].Schnitte[lauf2-1].posy >=0.0)
               {
                   isFlap(&FlapEta, &FlapAngle, input->Fluegel[lauf1].Schnitte[lauf2-1],lauf2-1, input->Klappen,0);
               }else{
                   isFlap(&FlapEta, &FlapAngle, input->Fluegel[lauf1].Schnitte[lauf2-1],lauf2-1, input->Klappen,1);
               }
               fprintf(fFWinp,"%lf    %lf\n",FlapEta,FlapAngle);
               fputs(HEADFW3, fFWinp);
               fprintf(fFWinp,"%.3lf   %.3lf   %.3lf   %.3lf   0.000    ",input->Fluegel[lauf1].Schnitte[lauf2].posx,input->Fluegel[lauf1].Schnitte[lauf2].posy,input->Fluegel[lauf1].Schnitte[lauf2].posz,input->Fluegel[lauf1].Schnitte[lauf2].tiefe);
               if (KoopelPanel!=0)
               {
                   fprintf(fFWinp," 220    ");
               }else{
                   fprintf(fFWinp," 100    ");
               }
               fprintf(fFWinp,"     %d   ",whichAirfoil(&input->Fluegel[lauf1],lauf2,Airfoils,anzAirfoils));
               if (input->Fluegel[lauf1].Schnitte[lauf2].posy >=0.0)
               {
                   isFlap(&FlapEta, &FlapAngle, input->Fluegel[lauf1].Schnitte[lauf2],lauf2 ,input->Klappen,1);
               }else{
                   isFlap(&FlapEta, &FlapAngle, input->Fluegel[lauf1].Schnitte[lauf2],lauf2 ,input->Klappen,0);
               }
               fprintf(fFWinp,"%lf    %lf\n\n",FlapEta,FlapAngle);
           }else{
               fprintf(fFWinp,"Panel #: %d. Number of elements across span n = %d and chord m = %d\n",anzPanel, input->Fluegel[lauf1].Schnitte[lauf2].AnzahlPan, input->Fluegel[lauf1].AnzahlTiefe);
               fprintf(fFWinp,"Neighbouring panels (0 for none) left: %d right: %d\n",anzPanel-1,anzPanel+1);
               anzPanel++;
               fputs(HEADFW2, fFWinp);
               fprintf(fFWinp,"%.3lf   %.3lf   %.3lf   %.3lf   0.000    ",input->Fluegel[lauf1].Schnitte[lauf2-1].posx,input->Fluegel[lauf1].Schnitte[lauf2-1].posy ,input->Fluegel[lauf1].Schnitte[lauf2-1].posz,input->Fluegel[lauf1].Schnitte[lauf2-1].tiefe);
               if ((lauf2-1)==0)
               {
                   fprintf(fFWinp," 100    ");
               }else {
                   fprintf(fFWinp," 220    ");
               }
               fprintf(fFWinp,"     %d    ",whichAirfoil(&input->Fluegel[lauf1],lauf2,Airfoils,anzAirfoils));
               if (input->Fluegel[lauf1].Schnitte[lauf2-1].posy >=0.0)
               {
                  isFlap(&FlapEta, &FlapAngle, input->Fluegel[lauf1].Schnitte[lauf2-1],lauf2-1, input->Klappen,0);
               }else{
                  isFlap(&FlapEta, &FlapAngle, input->Fluegel[lauf1].Schnitte[lauf2-1],lauf2-1, input->Klappen,1);
               }
               fprintf(fFWinp,"%lf    %lf\n",FlapEta,FlapAngle);
               fputs(HEADFW3, fFWinp);
               fprintf(fFWinp,"%.3lf   %.3lf   %.3lf   %.3lf   0.000    ",input->Fluegel[lauf1].Schnitte[lauf2].posx,input->Fluegel[lauf1].Schnitte[lauf2].posy,input->Fluegel[lauf1].Schnitte[lauf2].posz,input->Fluegel[lauf1].Schnitte[lauf2].tiefe);
               fprintf(fFWinp," 220    ");
               fprintf(fFWinp,"     %d   ",whichAirfoil(&input->Fluegel[lauf1],lauf2,Airfoils,anzAirfoils));
               if (input->Fluegel[lauf1].Schnitte[lauf2].posy >=0.0)
               {
                  isFlap(&FlapEta, &FlapAngle, input->Fluegel[lauf1].Schnitte[lauf2],lauf2, input->Klappen,1);
               }else{
                  isFlap(&FlapEta, &FlapAngle, input->Fluegel[lauf1].Schnitte[lauf2],lauf2, input->Klappen,0);
               }
               fprintf(fFWinp,"%lf    %lf\n\n",FlapEta,FlapAngle);
           }
       }
   }

   fprintf(fFWinp,"\n\n\n");
   fputs(HEADFW5, fFWinp);
   fputs(HEADFW6, fFWinp);
   fputs(HEADFW7, fFWinp);
   fputs(HEADFW8, fFWinp);
   fprintf(fFWinp,"\n");
   fputs(HEADFW9, fFWinp);
   fputs(HEADFW10, fFWinp);
   fputs(HEADFW11, fFWinp);
   fputs(HEADFW12, fFWinp);
   fputs(HEADFW13, fFWinp);
   fprintf(fFWinp,"\n");
   fputs(HEADFW14, fFWinp);
   fprintf(fFWinp,"##############\n");

  /* for (lauf1=0; lauf1 < anzPanel; lauf1++)
   {
       fprintf(fFWinp,"Panel #: %d. Number of elements across span n =
   }*/
   //checkDir();
   copyCamber(input, Airfoils,anzAirfoils);

}

void run_FWProzess (struct inputfile *input)
{
    create_FW_File(input);
}

