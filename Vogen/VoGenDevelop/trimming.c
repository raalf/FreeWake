/********************trimming.c*******************                 
*                                                *      
*      Vorentwurfs-Generator "VoGen"             *          
*      Version 2.00                              *
*      Date of last modification 22.06.2020      *
*                                                *          
*      Writen by: Jan Himisch                    *      
*             jan.himisch@dlr.de                 *      
*                                                *      
*                                                *
*************************************************/
#include "trimming.h"

void checkTrimm(struct inputfile *input)
{
  if ((input->CmSteuerFl > 0) && (input->CmSteuerKlappe > 0))
  {
    printf("Either a wing or a flap has to be defined for trimming\nNot both!\n");
    errmessage(1);
  }  
  if ((input->CmSteuerFl == 0) && (input->CmSteuerKlappe == 0) && (input->Austrimmen < 2 ))
  {
    printf("Either a wing or a flap has to be defined for trimming\n");
    errmessage(1);
  }
  if ((input->Austrimmen == 2 ) && (input->AnzahlMixVar < 3))
  {
    printf("At least three mixer-variables and functions are to be defined if Trimming mode 2 is used!\n");
    errmessage(1);
  }
  if (input->MaxCmIter < 3)
  {
    printf("With a lower number than 3 iterations trimming with VoGen doesn't work\n");
    errmessage(1);  
  }
  if (input->CmSteuerKlappe > 0)
  {
    if (input->CmSteuerKlappe > input->AnzahlKlappen)
    {
      printf("The spezified flap for trimming isn't definend\n");
      errmessage(1);  
    }
  }
  if (input->CmSteuerFl > 0)
  {
    if (input->CmSteuerFl < input->NbWings)
    {
      printf("The spezified wing for trimming isn't definend\n");
      errmessage(1);  
    }
  }
  if ((input->runOneByOne==0)&&((input->numberAlpha+input->numberCL)>1))
  {
    input->runOneByOne=1;
    printf("Warning: Run einzelnt (on=1/off=0)  has been switched on to meet requierments of\n trimming if more than one target Cl + target alfa is defined!\n");
  }
}


void Trimm_Schleife(struct inputfile *input)
{
  struct file14 resultfile;
  struct PLTdistribution PLTDIST;
  struct BM_FALL *BMDist=NULL;
  struct POLINTOUT PolintOut;
  struct METHODE *method;
  int lauf,lauf2,lauf3,checkIter;
  double dCMzuDeta;
  double *CL, *CM, *CN, *EPS, *EPS_2, *EPS_3, *OrgAlfasWing;
  float TrimMatrix[3][3];
  float cmln[3];
  FILE *fopen(),*fpara;  
  struct GESAMT_POLARint *GPolare;
  char EXECUTRstring[550] = {"\0"};
  int return_value;
  
  GPolare=(struct GESAMT_POLARint *)malloc((input->numberCL+input->numberAlpha)*sizeof(struct GESAMT_POLARint));
  fpara = fopen(input->EpsAusgabe,"w+");
  if (input->Austrimmen == 1)
  {
    fprintf(fpara," VARIABLES = \"Iter\" \"Cmx\" \"Cmy\" \"Cmz\" \"Eps\" \"DCm/DEps\"\n" );
  }
  if (input->Austrimmen == 2)
  {
    fprintf(fpara," VARIABLES = \"Iter\" \"Cmx\" \"Cmy\" \"Cmz\" \"X1\"  \"X2\"  \"X3\" \n" );
  }

  for (lauf=1; lauf<=input->MaxCmIter; lauf++)
  {    
    if (input->AnzahlMischer > 0)
    {
       run_mischer(input);
    }

    if (input->Struckt_Coup_ON == 1)
    {
      run_stru(input);
    }else{    
      if (input->OldInput == 0)
      {
        create_LiLiInp_File (input);    
      }
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
        BM_integration (input, &resultfile, &PLTDIST, BMDist, 0);
        write_BM_distribution(BMDist, &resultfile, input, 0);    
        method=(struct METHODE *)malloc((input->numberCL+input->numberAlpha)*sizeof(struct METHODE));
        methodeAlt(input, method, &resultfile); 
        write_load_dist(&resultfile, input);
        if (input->InterneInterpolation == 1)
        {
          PolarInterpolation(input,&PLTDIST,&resultfile, GPolare); 
        }
      }  
    }  

    if (input->Austrimmen == 1)
    {
      if (lauf==1)
      {
	CM=(double *)malloc(input->MaxCmIter*sizeof(double));
	EPS=(double *)malloc(input->MaxCmIter*sizeof(double));
	CM[lauf-1]=resultfile.Fall[0].CM;
	EPS[lauf-1]=0.0;
	EPS[lauf]=1.0;
	if (fabs(input->ZielCm-CM[lauf-1]) < 0.00001)
	{
	  lauf=input->MaxCmIter;
	  continue;
	}
	if (input->CmSteuerFl>0)
	{
	  OrgAlfasWing=(double *)malloc(input->Fluegel[input->CmSteuerFl-1].anzahlSchnitte*sizeof(double));
	  for (lauf2=0; lauf2 < input->Fluegel[input->CmSteuerFl-1].anzahlSchnitte; lauf2++)
	  {
	    OrgAlfasWing[lauf2]=input->Fluegel[input->CmSteuerFl-1].Schnitte[lauf2].AlphaPlus;
	    input->Fluegel[input->CmSteuerFl-1].Schnitte[lauf2].AlphaPlus=OrgAlfasWing[lauf2]+EPS[lauf];
	  }
	}
	if (input->CmSteuerKlappe>0)
	{
	  EPS[lauf-1]=input->Klappen[input->CmSteuerKlappe-1].winkel;
	  EPS[lauf]=input->Klappen[input->CmSteuerKlappe-1].winkel+1.0;
	  input->Klappen[input->CmSteuerKlappe-1].winkel=EPS[lauf];
	}
	fprintf(fpara,"%d %lf %lf 0.0\n",lauf,CM[lauf-1],EPS[lauf-1]);
      }else{
	CM[lauf-1]=resultfile.Fall[0].CM;
	dCMzuDeta=((CM[lauf-1]-CM[lauf-2])/(EPS[lauf-1]-EPS[lauf-2]));
	//dCMzuDeta=((CM[lauf-1]-CM[0])/(EPS[lauf-1]-EPS[0]));
	fprintf(fpara,"%d %lf %lf %lf\n",lauf,CM[lauf-1],EPS[lauf-1],dCMzuDeta);
	fflush(fpara);
	if (fabs(input->ZielCm-CM[lauf-1]) < 0.00001)
	{

	  lauf=input->MaxCmIter;
	  continue;
	}
	if (input->CmSteuerFl>0)
	{
	  if (lauf<input->MaxCmIter)
	  {
	    EPS[lauf]=(input->ZielCm-CM[lauf-1])/dCMzuDeta+EPS[lauf-1];
	    for (lauf2=0; lauf2 < input->Fluegel[input->CmSteuerFl-1].anzahlSchnitte; lauf2++)
	    {
	      input->Fluegel[input->CmSteuerFl-1].Schnitte[lauf2].AlphaPlus=OrgAlfasWing[lauf2]+EPS[lauf];
	    }
	  }else{
	    for (lauf2=0; lauf2 < input->Fluegel[input->CmSteuerFl-1].anzahlSchnitte; lauf2++)
	    {
	      input->Fluegel[input->CmSteuerFl-1].Schnitte[lauf2].AlphaPlus=OrgAlfasWing[lauf2];
	    }
	  }
	}
	if (input->CmSteuerKlappe>0)
	{
	  if (lauf<input->MaxCmIter)
	  {
	    EPS[lauf]=(input->ZielCm-CM[lauf-1])/dCMzuDeta+EPS[lauf-1];
	    input->Klappen[input->CmSteuerKlappe-1].winkel=EPS[lauf];
	  }else{
	    input->Klappen[input->CmSteuerKlappe-1].winkel=EPS[0];
	  }
	  if (input->Klappen[input->CmSteuerKlappe-1].Spiegeln > 0)
	  {
	    input->Klappen[input->Klappen[input->CmSteuerKlappe-1].Spiegeln].winkel=input->Klappen[input->CmSteuerKlappe-1].winkel;
	  }
	}
      }
    }else{
	if (lauf==1)
	{
	    CM=(double *)malloc(input->MaxCmIter*sizeof(double));
	    CN=(double *)malloc(input->MaxCmIter*sizeof(double));
	    CL=(double *)malloc(input->MaxCmIter*sizeof(double));
	    EPS=(double *)malloc(input->MaxCmIter*sizeof(double));
	    EPS_2=(double *)malloc(input->MaxCmIter*sizeof(double));
	    EPS_3=(double *)malloc(input->MaxCmIter*sizeof(double));
	    checkIter=1;
	    EPS[0]=input->MixVar[0];
	    EPS_2[0]=input->MixVar[1];
	    EPS_3[0]=input->MixVar[2];
	    lauf3=1;
	}
	CM[lauf-1]=resultfile.Fall[0].CM;
	CN[lauf-1]=resultfile.Fall[0].CN;
	CL[lauf-1]=resultfile.Fall[0].CL;

	if (checkIter==1)
	{
            /*if (  (fabs(CM[lauf-1]) < 0.00001) && (fabs(CN[lauf-1]) < 0.00001) && (fabs(CL[lauf-1]) < 0.00001) )
            {
                lauf=input->MaxCmIter;
                continue;
            }*/
            if (lauf==1)
            {
                input->MixVar[0]=input->MixVar[0]+input->InitStepSize;
                EPS[lauf3]=input->MixVar[0];
            }else{
                EPS[lauf3]=EPS[lauf3-1]+EPS[lauf3-1]-EPS[lauf3-2];
                input->MixVar[0]=EPS[lauf3];
            }
	}
        if (checkIter==2)
        {
            /*if (  (fabs(CM[lauf-1]) < 0.00001) && (fabs(CN[lauf-1]) < 0.00001) && (fabs(CL[lauf-1]) < 0.00001) )
            {
                lauf=input->MaxCmIter;
                continue;
            }*/
            TrimMatrix[0][0] = (CM[lauf-1]-CM[lauf-2])/(EPS[lauf3]-EPS[lauf3-1]);
            TrimMatrix[1][0] = (CL[lauf-1]-CL[lauf-2])/(EPS[lauf3]-EPS[lauf3-1]);
            TrimMatrix[2][0] = (CN[lauf-1]-CN[lauf-2])/(EPS[lauf3]-EPS[lauf3-1]);
            if (lauf==2)
            {
                input->MixVar[1]=input->MixVar[1]+input->InitStepSize;
                EPS_2[lauf3]=input->MixVar[1];
            }else{
                EPS_2[lauf3]=EPS_2[lauf3-1]+EPS_2[lauf3-1]-EPS_2[lauf3-2];
                input->MixVar[1]=EPS_2[lauf3];
            }
        }
        if (checkIter==3)
        {
            /*if (  (fabs(CM[lauf-1]) < 0.00001) && (fabs(CN[lauf-1]) < 0.00001) && (fabs(CL[lauf-1]) < 0.00001) )
            {
                lauf=input->MaxCmIter;
                continue;
            }*/
            TrimMatrix[0][1] = (CM[lauf-1]-CM[lauf-2])/(EPS_2[lauf3]-EPS_2[lauf3-1]);
            TrimMatrix[1][1] = (CL[lauf-1]-CL[lauf-2])/(EPS_2[lauf3]-EPS_2[lauf3-1]);
            TrimMatrix[2][1] = (CN[lauf-1]-CN[lauf-2])/(EPS_2[lauf3]-EPS_2[lauf3-1]);
            if (lauf==3)
            {
                input->MixVar[2]=input->MixVar[2]+input->InitStepSize;
                EPS_3[lauf3]=input->MixVar[2];
            }else{
                EPS_3[lauf3]=EPS_3[lauf3-1]+EPS_3[lauf3-1]-EPS_3[lauf3-2];
                input->MixVar[2]=EPS_3[lauf3];
            }
        }
        if (checkIter==4)
        {
            /*if (  (fabs(CM[lauf-1]) < 0.00001) && (fabs(CN[lauf-1]) < 0.00001) && (fabs(CL[lauf-1]) < 0.00001) )
            {
                lauf=input->MaxCmIter;
                continue;
            }*/
            TrimMatrix[0][2] = (CM[lauf-1]-CM[lauf-2])/(EPS_3[lauf3]-EPS_3[lauf3-1]);
            TrimMatrix[1][2] = (CL[lauf-1]-CL[lauf-2])/(EPS_3[lauf3]-EPS_3[lauf3-1]);
            TrimMatrix[2][2] = (CN[lauf-1]-CN[lauf-2])/(EPS_3[lauf3]-EPS_3[lauf3-1]);
            cmln[0]=CM[lauf-1]*(-1.0);
            cmln[1]=CL[lauf-1]*(-1.0);
            cmln[2]=CN[lauf-1]*(-1.0);
            LinGlLoesung(&TrimMatrix[0][0],&cmln[0],3);
            lauf3++;
            EPS[lauf3]=input->MixVar[0]+cmln[0];
            EPS_2[lauf3]=input->MixVar[1]+cmln[1];
            EPS_3[lauf3]=input->MixVar[2]+cmln[2];
            input->MixVar[0]=EPS[lauf3];
            input->MixVar[1]=EPS_2[lauf3];
            input->MixVar[2]=EPS_3[lauf3];
        }
        fprintf(fpara,"%d %lf %lf %lf %lf %lf %lf\n",lauf,CL[lauf-1],CM[lauf-1],CN[lauf-1],EPS[lauf3-1],EPS_2[lauf3-1],EPS_3[lauf3-1]);
        if (checkIter==4)
        {
            checkIter=1;
            lauf3++;
        }else{
            checkIter++;
        }
    }
    if (input->ArotAdjust==1)
    {
      input->Alfa_Rot=resultfile.Fall[0].Alfa;
    }
  }  
  fclose(fpara);
  strcpy(EXECUTRstring, "mv \0");
  strcat(EXECUTRstring, input->EpsAusgabe);
  strcat(EXECUTRstring, " LaufV1.lili.");
  strcat(EXECUTRstring,input->LiLiVersionName); 
  strcat(EXECUTRstring, " \0");
    
  return_value = system(EXECUTRstring);
  free(CM);
  free(EPS); 
  if (input->CmSteuerFl>0)
  {
    free(OrgAlfasWing);
  }

}
