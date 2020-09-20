/********************mischer.c********************                  
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
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "mischer.h"

double getXvalue(double *MixVar, char Xvar[20])
{
  int xWert;
  sscanf(&Xvar[1],"%d",&xWert);
      //printf("Xvar:%s   XWert:%d\n",Xvar,xWert);
  return(MixVar[xWert-1]);
}

void setregelNull(int Start,char *regelZeig)
{
  int lauf;
  for (lauf=Start;lauf<MaxRegelLength;lauf++)
  {
    if (regelZeig[lauf]=='\0')
    {
      break;
    }
    regelZeig[lauf]='\0';
  }
}

void Minus2PlusMinus(char *regelZeig)
{
  int lauf,lauf2;
  char regelzwi[MaxRegelLength]={"\0"};  
  char regelzwi2[MaxRegelLength]={"\0"};
  if (strcspn(regelZeig,"-") < (strlen(regelZeig)-1))
  {
    lauf2=0;
    for(lauf=0;lauf<strlen(regelZeig);lauf++)
    {
      if ((regelZeig[lauf]=='-') && (lauf!=0))
      {
        if (isdigit(regelZeig[lauf-1]))
        {
          regelzwi[lauf2]='+';
          lauf2++;
          regelzwi[lauf2]=regelZeig[lauf];
          lauf2++;
        }else{
          regelzwi[lauf2]=regelZeig[lauf];
          lauf2++;
        }        
      }else{
        regelzwi[lauf2]=regelZeig[lauf];
        lauf2++;
      }      
    }
    for(lauf=0;lauf<MaxRegelLength;lauf++)
    {
      if (regelzwi[lauf]=='\0')
      {
        break;
      }
      regelZeig[lauf]=regelzwi[lauf];
    }    
   
    /*lauf2=0;
    for(lauf=0;lauf<(strlen(regelZeig));lauf++)
    {
      if (regelZeig[lauf]=='-')
      {
        if (regelZeig[lauf+1]=='-')
        {
          regelzwi2[lauf2]='+';
          lauf++;
          lauf2++;                
        }else{
          regelzwi2[lauf2]=regelZeig[lauf];
          lauf2++;
        }
        
      }else{
        regelzwi2[lauf2]=regelZeig[lauf];
        lauf2++;
      }
    }
    for(lauf=0;lauf<MaxRegelLength;lauf++)
    {
      regelZeig[lauf]=regelzwi2[lauf];
    } */
    
           
  }
}

void StringPlusMinus(char *regelZeig, double *MixVar)
{
  int lauf,lauf2;
  char a[20]= {"\0"}, b[20]= {"\0"},zwi[20]= {"\0"},zwi2[MaxRegelLength]= {"\0"};
  double aWert,bWert,abPlusMinus;
  int Lstop, Lstart;  
  Minus2PlusMinus(regelZeig);
  //printf("RegelForPlusMinus: %s\n",regelZeig);
  if (strcspn(regelZeig,"+") < (strlen(regelZeig)-1))
  {
    lauf=1;
    do 
    {
      b[lauf-1]=regelZeig[strcspn(regelZeig,"+")+lauf];
      lauf++;
    }
    while ((regelZeig[lauf] != '\0') && ( isdigit(regelZeig[strcspn(regelZeig,"+")+lauf]) || (regelZeig[strcspn(regelZeig,"+")+lauf]=='x') || (regelZeig[strcspn(regelZeig,"+")+lauf]=='.') || ((regelZeig[strcspn(regelZeig,"+")+lauf]=='-')&& (lauf==2))));  
    Lstop=strcspn(regelZeig,"+")+lauf;
    lauf=1;
    do 
    {
      zwi[lauf-1]=regelZeig[strcspn(regelZeig,"+")-lauf];
      lauf++;
    }
    while ((regelZeig[lauf] != '\0') && ( isdigit(regelZeig[strcspn(regelZeig,"+")-lauf]) || (regelZeig[strcspn(regelZeig,"+")-lauf]=='x') || (regelZeig[strcspn(regelZeig,"+")-lauf]=='.') || (regelZeig[strcspn(regelZeig,"+")-lauf]=='-') ));
    Lstart=strcspn(regelZeig,"+")-lauf+1;
    for(lauf2=0;lauf2<(lauf-1);lauf2++)
    {
      a[lauf2]=zwi[lauf-2-lauf2];
    }
    if (a[0]!='x')
    {
      sscanf(a,"%lf",&aWert);
    }else{
      aWert=getXvalue(MixVar,a);
    }
    if (b[0]!='x')
    {
      sscanf(b,"%lf",&bWert);
    }else{
      bWert=getXvalue(MixVar,b);     
    }
    //printf("aPM:%lf  bPM:%lf   a:%s   b:%s\n",aWert,bWert,a,b);
    if (regelZeig[strcspn(regelZeig,"+")]=='+')
    {
      abPlusMinus=aWert+bWert;
    }else{
      abPlusMinus=aWert-bWert;
    }
    
    sprintf(zwi,"%lf",abPlusMinus);
    for (lauf=Lstop;lauf<strlen(regelZeig); lauf++)
    {
      zwi2[lauf-Lstop]=regelZeig[lauf];
    }
    strcpy(&regelZeig[Lstart],zwi);
    strcpy(&regelZeig[Lstart+strlen(zwi)],zwi2); 
    setregelNull((Lstart+strlen(zwi)+strlen(zwi2)),regelZeig);
    //printf("PlusMinusRegel:%s\n",regelZeig);
    if (strcspn(regelZeig,"+") < (strlen(regelZeig)-1))
    {
      StringPlusMinus(regelZeig,MixVar);
    }
  }
}

void StringMalGeteilt(char *regelZeig, double *MixVar)
{
  int lauf,lauf2;
  char a[20]= {"\0"}, b[20]= {"\0"},zwi[20]= {"\0"},zwi2[MaxRegelLength]= {"\0"};
  double aWert,bWert,abPlusMinus;
  int Lstop, Lstart;  
  
  Minus2PlusMinus(regelZeig);
  //printf("RegelForPlusMinus: %s\n",regelZeig);
  if (strcspn(regelZeig,"*/") < (strlen(regelZeig)-1))
  {
    lauf=1;
    do 
    {
      b[lauf-1]=regelZeig[strcspn(regelZeig,"*/")+lauf];
      lauf++;
    }
    while ((regelZeig[lauf] != '\0') && ( isdigit(regelZeig[strcspn(regelZeig,"*/")+lauf]) || (regelZeig[strcspn(regelZeig,"*/")+lauf]=='x') || (regelZeig[strcspn(regelZeig,"*/")+lauf]=='.') || ((regelZeig[strcspn(regelZeig,"*/")+lauf]=='-')&& (lauf==2))));  
    Lstop=strcspn(regelZeig,"*/")+lauf;
    lauf=1;
    do 
    {
      zwi[lauf-1]=regelZeig[strcspn(regelZeig,"*/")-lauf];
      if (zwi[lauf-1]=='-')
      {
        lauf++;
        break;
      }
      lauf++;
    }
    while ((regelZeig[lauf] != '\0') && ( isdigit(regelZeig[strcspn(regelZeig,"*/")-lauf]) || (regelZeig[strcspn(regelZeig,"*/")-lauf]=='x') || (regelZeig[strcspn(regelZeig,"*/")-lauf]=='.') || (regelZeig[strcspn(regelZeig,"*/")-lauf]=='-') ));
    Lstart=strcspn(regelZeig,"*/")-lauf+1;
    for(lauf2=0;lauf2<(lauf-1);lauf2++)
    {
      a[lauf2]=zwi[lauf-2-lauf2];
    }
    if (a[0]!='x')
    {
      sscanf(a,"%lf",&aWert);
    }else{
      aWert=getXvalue(MixVar,a);
    }
    if (b[0]!='x')
    {
      sscanf(b,"%lf",&bWert);
    }else{
      bWert=getXvalue(MixVar,b);     
    }
    //printf("aPM:%lf  bPM:%lf   a:%s   b:%s\n",aWert,bWert,a,b);
    if (regelZeig[strcspn(regelZeig,"*/")]=='*')
    {
      abPlusMinus=aWert*bWert;
    }else{
      abPlusMinus=aWert/bWert;
    }
    
    sprintf(zwi,"%lf",abPlusMinus);
    for (lauf=Lstop;lauf<strlen(regelZeig); lauf++)
    {
      zwi2[lauf-Lstop]=regelZeig[lauf];
    }
    strcpy(&regelZeig[Lstart],zwi);
    strcpy(&regelZeig[Lstart+strlen(zwi)],zwi2); 
    setregelNull((Lstart+strlen(zwi)+strlen(zwi2)),regelZeig);
    //printf("PlusMinusRegel:%s\n",regelZeig);
    if (strcspn(regelZeig,"*/") < (strlen(regelZeig)-1))
    {
      StringMalGeteilt(regelZeig,MixVar);
    }
  }
}

void StringPow(char *regelZeig, double *MixVar)
{
  int lauf,lauf2;
  if (strcspn(regelZeig,"*") < (strlen(regelZeig)-1))  //Achtung hier muss ich ran, da nur das erste * gecheckt wird!
  {
    if (regelZeig[strcspn(regelZeig,"*")+1] == '*')
    {
      char a[20]= {"\0"}, b[20]= {"\0"},zwi[20]= {"\0"},zwi2[100]= {"\0"};
      double aWert,bWert,abPow;
      int Lstop, Lstart;
      lauf=1;
      do 
      {
        b[lauf-1]=regelZeig[strcspn(regelZeig,"*")+lauf+1];
        lauf++;
      }
      while ((regelZeig[lauf] != '\0') && ( isdigit(regelZeig[strcspn(regelZeig,"*")+lauf+1]) || (regelZeig[strcspn(regelZeig,"*")+lauf+1]=='x') || (regelZeig[strcspn(regelZeig,"*")+lauf+1]=='.') ));  
      Lstop=strcspn(regelZeig,"*")+lauf;
      lauf=1;
      do 
      {
        zwi[lauf-1]=regelZeig[strcspn(regelZeig,"*")-lauf];
        lauf++;
      }
      while ((regelZeig[lauf] != '\0') && ( isdigit(regelZeig[strcspn(regelZeig,"*")-lauf]) || (regelZeig[strcspn(regelZeig,"*")-lauf]=='x') || (regelZeig[strcspn(regelZeig,"*")-lauf]=='.'))); 
      Lstart=strcspn(regelZeig,"*")-lauf+1;
      for(lauf2=0;lauf2<(lauf-1);lauf2++)
      {
        a[lauf2]=zwi[lauf-2-lauf2];
      }
      if (a[0]!='x')
      {
        sscanf(a,"%lf",&aWert);
      }else{
        aWert=getXvalue(MixVar,a);
      }
      if (b[0]!='x')
      {
        sscanf(b,"%lf",&bWert);
      }else{
        bWert=getXvalue(MixVar,b);      
      }
      abPow=pow(aWert,bWert);
      sprintf(zwi,"%lf",abPow);   
      for (lauf=Lstop;lauf<strlen(regelZeig); lauf++)
      {
        zwi2[lauf-Lstop]=regelZeig[lauf];
      }
      strcpy(&regelZeig[Lstart],zwi);
      strcpy(&regelZeig[Lstart+strlen(zwi)],zwi2); 
      setregelNull((Lstart+strlen(zwi)+strlen(zwi2)),regelZeig);
      StringPow(regelZeig, MixVar);
    }else{
      StringPow(&regelZeig[strcspn(regelZeig,"*")+1],MixVar);
    }
  }
}

void StringSin(char *regelZeig, double *MixVar)
{
  int lauf,lauf2;
  double aWert,Ergebnis;
  char Such[] = {"sin"};
  char *RegelZeiger;
  char zwi[20]= {"\0"},zwi3[20]= {"\0"}, zwi2[MaxRegelLength]= {"\0"};
   
  if (strstr(regelZeig,Such)!=NULL )
  {
    RegelZeiger=strstr(regelZeig,Such);
    lauf=0;
    lauf2=0;
    do 
    {
      zwi[lauf]=RegelZeiger[lauf+3];
      lauf++;
    }
    while((RegelZeiger[lauf+3]=='x') || (RegelZeiger[lauf+3]=='.') || (isdigit(regelZeig[lauf])));
    do 
    {
      zwi2[lauf2]=RegelZeiger[lauf+3];
      lauf++;
      lauf2++;
    }
    while(RegelZeiger[lauf+3]!='\0');
    setregelNull(0,RegelZeiger);
    if (zwi[0]!='x')
    {
      sscanf(zwi,"%lf",&aWert);
    }else{
      aWert=getXvalue(MixVar,zwi);
    }
    Ergebnis=sin(aWert/180*PI);
    sprintf(zwi3,"%lf",Ergebnis);
    strcpy(&RegelZeiger[0],zwi3);
    strcpy(&RegelZeiger[strlen(zwi3)],zwi2);
    setregelNull((strlen(zwi3)+strlen(zwi2)),RegelZeiger);
    if (strstr(regelZeig,Such)!=NULL )
    {
      StringSin(regelZeig,MixVar);
    }
  }  
}

void StringCos(char *regelZeig, double *MixVar)
{
  int lauf,lauf2;
  double aWert,Ergebnis;
  char Such[] = {"cos"};
  char *RegelZeiger;
  char zwi[20]= {"\0"},zwi3[20]= {"\0"}, zwi2[MaxRegelLength]= {"\0"};
   
  if (strstr(regelZeig,Such)!=NULL )
  {
    RegelZeiger=strstr(regelZeig,Such);
    lauf=0;
    lauf2=0;
    do 
    {
      zwi[lauf]=RegelZeiger[lauf+3];
      lauf++;
    }
    while((RegelZeiger[lauf+3]=='x') || (RegelZeiger[lauf+3]=='.') || (isdigit(regelZeig[lauf])));
    do 
    {
      zwi2[lauf2]=RegelZeiger[lauf+3];
      lauf++;
      lauf2++;
    }
    while(RegelZeiger[lauf+3]!='\0');
    setregelNull(0,RegelZeiger);
    if (zwi[0]!='x')
    {
      sscanf(zwi,"%lf",&aWert);
    }else{
      aWert=getXvalue(MixVar,zwi);
    }
    Ergebnis=cos(aWert/180*PI);
    sprintf(zwi3,"%lf",Ergebnis);
    strcpy(&RegelZeiger[0],zwi3);
    strcpy(&RegelZeiger[strlen(zwi3)],zwi2);
    setregelNull((strlen(zwi3)+strlen(zwi2)),RegelZeiger);
    if (strstr(regelZeig,Such)!=NULL )
    {
      StringSin(regelZeig,MixVar);
    }
  }  
}

void StringTan(char *regelZeig, double *MixVar)
{
  int lauf,lauf2;
  double aWert,Ergebnis;
  char Such[] = {"tan"};
  char *RegelZeiger;
  char zwi[20]= {"\0"},zwi3[20]= {"\0"}, zwi2[MaxRegelLength]= {"\0"};
   
  if (strstr(regelZeig,Such)!=NULL )
  {
    RegelZeiger=strstr(regelZeig,Such);
    lauf=0;
    lauf2=0;
    do 
    {
      zwi[lauf]=RegelZeiger[lauf+3];
      lauf++;
    }
    while((RegelZeiger[lauf+3]=='x') || (RegelZeiger[lauf+3]=='.') || (isdigit(regelZeig[lauf])));
    do 
    {
      zwi2[lauf2]=RegelZeiger[lauf+3];
      lauf++;
      lauf2++;
    }
    while(RegelZeiger[lauf+3]!='\0');
    setregelNull(0,RegelZeiger);
    if (zwi[0]!='x')
    {
      sscanf(zwi,"%lf",&aWert);
    }else{
      aWert=getXvalue(MixVar,zwi);
    }
    Ergebnis=tan(aWert/180*PI);
    sprintf(zwi3,"%lf",Ergebnis);
    strcpy(&RegelZeiger[0],zwi3);
    strcpy(&RegelZeiger[strlen(zwi3)],zwi2);
    setregelNull((strlen(zwi3)+strlen(zwi2)),RegelZeiger);
    if (strstr(regelZeig,Such)!=NULL )
    {
      StringSin(regelZeig,MixVar);
    }
  }  
}



double funktionsAnlyse(double *MixVar, char regel[MaxRegelLength], int dbLevel)
{
  char zwiRegel[MaxRegelLength]={"\0"};
  char zwiRegel2[MaxRegelLength]={"\0"};
  char zwiChar[50]={"\0"};
  int lauf,lauf2,lauf3,lauf4,KLauf,KLcheckauf;
  double Ergebnis,Return;
  if (dbLevel==3 || dbLevel>99)
  {
    printf("Haupt-Funktion-Start:%s \n",regel);
  }
  KLauf=0;
  KLcheckauf=-1;
  //Klammern
  lauf=0;
  do
  {
    if ( regel[lauf] == '(' )
    {
      KLauf++;
      if (KLauf==1)
      {
        KLcheckauf=lauf;
      }
    }
    if ( regel[lauf]== ')' )
    {
       KLauf--;
       if ( (KLauf == 0) && (KLcheckauf > -1) ) //Check ob die aktuell äusserste Klammer geschlossen wurde
       {
         for (lauf2=KLcheckauf+1;lauf2<lauf;lauf2++)
         {
            zwiRegel[lauf2-KLcheckauf-1]=regel[lauf2];
            zwiRegel[lauf2-KLcheckauf]='\0';
         }
         lauf2=0;
         lauf3=0;
         for (lauf4=0;lauf4<MaxRegelLength;lauf4++)
         {
           if (zwiRegel2[lauf4]=='\0')
           {
             break;
           }
           zwiRegel2[lauf4]='\0';
         }
         lauf4=0;
         //printf("zwiRegel2 befor read:%s\n",zwiRegel2);
         do 
         {
           if (lauf2 < KLcheckauf) 
           {
              regel[lauf3]=regel[lauf2];
              lauf3++;
           }
           if (lauf2 == KLcheckauf)
           {                          
             Ergebnis=funktionsAnlyse(MixVar,zwiRegel,dbLevel);
             if (dbLevel==3 || dbLevel>99)
             {
               printf("\nErgebnis Funktionsanalyse nach zwischen Aufruf:%s\n", Ergebnis);
               printf("Regel Funktionsanalyse nach zwischen Aufruf:%s\n\n", regel);
             }
             sprintf(zwiChar,"%lf",Ergebnis);
             lauf3=lauf3+strlen(zwiChar);
           }
           
           if (lauf2>lauf)
           {
              zwiRegel2[lauf4]=regel[lauf2];
              lauf4++; 
              lauf3++;          
           }
           lauf2++;
         }
         while((regel[lauf2] != '\0') && (lauf2<MaxRegelLength) );
         //printf("zwiChar:%s\t",zwiChar);
         //printf("zwiRegel2:%s\n",zwiRegel2);
         strcpy(&regel[KLcheckauf],zwiChar);
         if (lauf4>0)
         {
           strcpy(&regel[KLcheckauf+strlen(zwiChar)],zwiRegel2);       
           setregelNull((KLcheckauf+strlen(zwiChar)+strlen(zwiRegel2)),regel);
         }else{
           setregelNull((KLcheckauf+strlen(zwiChar)),regel);    
         }
         KLcheckauf=-1; //??? 
       }  //Ende äusserste Klammer
    } //Klammer Zu Check Ende
    lauf++;
  }
  while((regel[lauf] != '\0') && (lauf<MaxRegelLength) );
  //printf("ZwiFunktion nach Klammer:%s\n",regel);

  if (strcspn(regel,"(") < (strlen(regel)-1))
  {
    Ergebnis=funktionsAnlyse(MixVar,regel,dbLevel);
  }  
    
  //SinCosTan
  StringSin(regel,MixVar);
  StringCos(regel,MixVar);
  StringTan(regel,MixVar);
  
  //Hochzahlen
  StringPow(regel,MixVar);
  //printf("ZwiFunktion nach Hoch:%s\n",regel);
    
  //Malgeteilt
  StringMalGeteilt(regel,MixVar);
  //printf("ZwiFunktion nach MalGeteilt:%s\n",regel);

  //Plus Minus
  StringPlusMinus(regel,MixVar);
    
    
  if (regel[0]!='x')
  {
    sscanf(regel,"%lf",&Return);
  }else{
    Return=getXvalue(MixVar,regel);
  }      
  //printf("ZwiFunktionEnd:%lf\n",Return);
  return(Return);
}

void run_mischer(struct inputfile *input)
{
  int lauf,lauf2,zahl,FuncTarget;
  double FunctionsWert;
  char zwiTarget[5];
  char orgFunktion [MaxRegelLength];
  char *ZEIGER,*ZEIGER2;
  
  if (input->debuglevel==3 || input->debuglevel>99)
  {
    for(lauf=0;lauf<input->AnzahlMixVar;lauf++)
    {
      printf("X%d=%lf\n",lauf+1,input->MixVar[lauf]);
    }
    printf("\n");
    for(lauf=0;lauf<input->AnzahlMischer;lauf++)
    {
      printf("MixerRule:%d :   %s\n",lauf+1,input->Mischer[lauf].Regel);
    }
  }

  for(lauf=0;lauf<input->AnzahlMischer;lauf++)
  {
    for (lauf2=0;lauf2<MaxRegelLength;lauf2++)
    {
        orgFunktion[lauf2]=input->Mischer[lauf].Regel[lauf2];
    }
    zahl=0;
    do
    {
      zwiTarget[zahl]  = input->Mischer[lauf].Regel[zahl];
      zahl++;
    }
    while ((input->Mischer[lauf].Regel[zahl] != '=') && (input->Mischer[lauf].Regel[zahl] != '\0'));
    //printf("\nRegel vor Aufruf:%s %s\n",zwiTarget,input->Mischer[lauf].Regel);
    ZEIGER=&input->Mischer[lauf].Regel[zahl+1];
    //printf("\nRegel vor Aufruf:%s\n",ZEIGER);
    if ( strncmp( zwiTarget , "f", 1) != 0)
    {
      errmessage(66);
    }
    ZEIGER2=&zwiTarget[1];
    sscanf(ZEIGER2, "%d", &FuncTarget);
    if (input->debuglevel==3 || input->debuglevel>99)
    {
      printf("\nFunctionsZiel:%d\n",FuncTarget);
      printf("Regel vor Aufruf:%s\n",ZEIGER);
    }
    FunctionsWert=funktionsAnlyse(input->MixVar,ZEIGER,input->debuglevel);
    if (input->debuglevel==3 || input->debuglevel>99)
    {
       printf("Functions Ergebnis:%lf\n",FunctionsWert);
       printf("\n\n\n");
    }
    input->Klappen[FuncTarget-1].winkel=FunctionsWert;
    for (lauf2=0;lauf2<MaxRegelLength;lauf2++)
    {
        input->Mischer[lauf].Regel[lauf2]=orgFunktion[lauf2];
    }
  }
}


