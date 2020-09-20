/******************run_Polint.h*******************                  
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


#include "run_polar.h"

void run_einzelnt(struct inputfile *input)
{
  int Numb_CL, Numb_Alfa,lauf;
  double *tCl;    
  double *tAlpha; 
  struct file14 resultfile;
  struct POLINTOUT PolintOut;
  struct PLTdistribution PLTDIST;
  struct BM_FALL *BMDist;
  struct METHODE *method;
  char EXECUTRstring[550] = {"\0"};
  char laufNR[4];
  int laufD;
  int return_value;
  struct GESAMT_POLARint *GPolare;

  GPolare=(struct GESAMT_POLARint *)malloc((input->numberCL+input->numberAlpha)*sizeof(struct GESAMT_POLARint));
  Numb_CL=input->numberCL;
  Numb_Alfa=input->numberAlpha;
  if (Numb_CL>0)
  {
    tCl = (double *)malloc(Numb_CL*sizeof(double));
    input->numberAlpha=0;
    input->numberCL=1;
    for (lauf=0; lauf<Numb_CL;lauf++)
    { 
      tCl[lauf]=input->targetCl[lauf];
    }
    for (lauf=0; lauf<Numb_CL;lauf++)
    { 
      input->targetCl[0]=tCl[lauf];
      if (input->NumberOfParametSetsinput>0)                                                                           
      {
        run_parametric(input,  lauf+1);
      }else if (input->Austrimmen > 1){
        Trimm_Schleife(input);  
      }else if (input->Struckt_Coup_ON == 1){    
        run_stru(input);
      }else{      
          run_PROZESS(input);
      }
      strcpy(EXECUTRstring, "mv LaufV1.lili.");
      strcat(EXECUTRstring,input->LiLiVersionName); 
      strcat(EXECUTRstring, " LaufV1.lili.");
      strcat(EXECUTRstring,input->LiLiVersionName); 
      strcat(EXECUTRstring, "_ClA\0");
      
      laufD=lauf+1;
      if (laufD>=100)
      {
        errmessage(46);
      }
      if (laufD>=10)
      {      
    
        laufNR[0]=48+(int)(laufD/10);
        laufNR[1]=48+laufD-((int)(laufD/10)*10);
        laufNR[2]='\0';
      
        /*double *zahl;
        *zahl=10;
        dzahl=modf(laufD, zahl);
        laufNR[0]=48+*zahl;
        laufNR[1]= 48 + fmod(laufD, *zahl*10);
        laufNR[2]='\0';*/
      } else {
        laufNR[0]=48+laufD;
        laufNR[1]='\0';
      }     
      strcat(EXECUTRstring, laufNR);
      return_value = system(EXECUTRstring); 
    }
  }


  if (Numb_Alfa>0)
  {
    tAlpha = (double *)malloc(Numb_Alfa*sizeof(double));
    input->numberAlpha=1;
    input->numberCL=0;
    for (lauf=0; lauf<Numb_Alfa;lauf++)
    { 
      tAlpha[lauf]=input->targetAlpha[lauf];
    }
    for (lauf=0; lauf<Numb_Alfa;lauf++)
    { 
      input->targetAlpha[0]=tAlpha[lauf];
      if (input->NumberOfParametSetsinput>0)                                                                           
      {
        run_parametric(input, lauf+1+Numb_CL);
      }else if (input->Austrimmen > 1){
        Trimm_Schleife(input);  
      }else if (input->Struckt_Coup_ON == 1){    
        run_stru(input);
      }else{      
        create_LiLiInp_File (input);    
        if (input->LiliStart==1)
        {
          run_Lili(input);
          read14 (&resultfile , input);
          readDistributionPLT (&PLTDIST, input); 
          if (input->PolintStart==1)
          {
            run_Polint(input);
            read_Polint_Out (input, &PolintOut);
          }
          BMDist=(struct BM_FALL *)malloc((input->numberCL+input->numberAlpha)*sizeof(struct BM_FALL));
          BM_integration (input, &resultfile, &PLTDIST, BMDist,0);
          write_BM_distribution(BMDist, &resultfile, input, 0);    
          method=(struct METHODE *)malloc((input->numberCL+input->numberAlpha)*sizeof(struct METHODE));
          methodeAlt(input, method, &resultfile); 
          write_load_dist(&resultfile, input);   
        }
      }
    }
    strcpy(EXECUTRstring, "mv LaufV1.lili.");
    strcat(EXECUTRstring,input->LiLiVersionName); 
    strcat(EXECUTRstring, " LaufV1.lili.");
    strcat(EXECUTRstring,input->LiLiVersionName); 
    strcat(EXECUTRstring, "_ClA\0");
    
    laufD=lauf+1;
    if (laufD>=100)
    {
      errmessage(46);
    }
    if (laufD>=10)
    {      
      laufNR[0]=48+(int)(laufD/10);
      laufNR[1]=48+laufD-((int)(laufD/10)*10);
      laufNR[2]='\0';
      
      /*double *zahl;
      *zahl=10;
      dzahl=modf(laufD, zahl);
      laufNR[0]=48+*zahl;
      laufNR[1]= 48 + fmod(laufD, *zahl*10);
      laufNR[2]='\0';*/
    } else {
      laufNR[0]=48+laufD;
      laufNR[1]='\0';
    }     
    strcat(EXECUTRstring, laufNR);
    return_value = system(EXECUTRstring); 
  }
  strcpy(EXECUTRstring, "mkdir LaufV1.lili.");
  strcat(EXECUTRstring,input->LiLiVersionName); 
  strcat(EXECUTRstring, "\0");  
  return_value = system(EXECUTRstring); 
  
  strcpy(EXECUTRstring, "mv  LaufV1.lili."); 
  strcat(EXECUTRstring,input->LiLiVersionName); 
  strcat(EXECUTRstring,"_* LaufV1.lili.");
  strcat(EXECUTRstring,input->LiLiVersionName);
  strcat(EXECUTRstring, "\0");  
  return_value = system(EXECUTRstring); 
}
