//this file includes all the subroutine that handle writing to files.
//the path to the directory in which the files are stored is defined in
//OUTPUT_PATH and stored in info.output, which is defined in general.h

//deletes previous timestep files
void Delete_timestep();
//saves information on elementary wings to file
void Save_Elementary_Wings(const GENERAL,const BOUND_VORTEX* );
//saves information of trailing edge to file
void Save_Trailing_Edge(const GENERAL info,const BOUND_VORTEX* trailedgePtr);
//saves results of Horstmann's method to file
void Horstmann_Results(const GENERAL info,const BOUND_VORTEX* ,const double,\
					   const double,const double,const double,\
					   const double,const double,const double);
void Header(const GENERAL info,const BOUND_VORTEX* ,const double,\
					   const double,const double,const double,\
					   const double,const double,const double);
//saves results from time-stepping method to file
void Time_Stepping_Results(const GENERAL,double **,double *,const double,\
						   const double,const double,const double,\
						   const double);
//saves final result from time-stepping method to file
void Time_Stepping_End_Results(const GENERAL,const int,const int,double **,\
							double *,const double,const double,\
						   const double,const double,const double);
//save information of surface DVEs to file
void Save_Surface_DVEs(const GENERAL info,const DVE *surfacePtr);
//save input file to output directory
void Save_Input_File(const char [126],const char [126]);
//save input file and header to config directory
void Save_Config_Head_File(GENERAL);
//saves results of current timestep to file
void Save_Timestep(const GENERAL,const int,DVE **,const DVE *,double **);
//saves forces and moments of surface DVEs
void Save_SurfaceDVE_Loads(const GENERAL,const int,const DVE *);

void CreateQuiverFile(const double[3], const double[3],const int,const int);
//saves flight condition files
void SaveSpanDVEInfo(PANEL *,DVE *,DVE **,STRIP *,double **,const int,const int );

//===================================================================//
		//START OF File_Initializing
//===================================================================//
void Delete_timestep()
{
	//deletes previous timestep files

	char comand[160];	//system command to delete previous timestep files

	//deletes previous timestep files
	sprintf(comand,"%s%s%s","del ",info.output,"timestep*");
	system(comand);
}
//===================================================================//
		//END OF File_Initializing
//===================================================================//
//===================================================================//
		//START OF Save_Elementary_Wings
//===================================================================//
void Save_Elementary_Wings(const GENERAL info,const BOUND_VORTEX* elementPtr)
{
//save information on elementary wings to file

	int l;			//loop counter
	FILE *fp;		//output file
	char filename[160];	//file path and name

	//creates file "Elementary_Wings.txt in directory "output"
	sprintf(filename,"%s%s",info.output,"Elementary_Wings.txt");
	fp = fopen(filename, "w");

	//writes header
	fprintf(fp,"\n\n\nProgram Version: %s\n",PROGRAM_VERSION);
	fprintf(fp,"\tlifting line ctr\t\tcontrol point\thalf-span");
	fprintf(fp,"\tsweep\tdihedral\n");
	fprintf(fp,"n\txo\t\tyo\t\tzo\t\txA\t\tyA\tzA\t\teta\t\tchord\t");
	fprintf(fp,"\tphi\t\tnu\t\t\tsurface normal\t\t\t\tlocal vel\n");

	for(l=0; l<info.noelement; l++)
		{
			//elementary wing number
			fprintf(fp, "%d\t",l);
			//coord. of lifting line center
			fprintf(fp, "%lf\t%lf\t%lf\t",\
			elementPtr[l].xo[0],elementPtr[l].xo[1],elementPtr[l].xo[2]);
			//coord. of control point
			fprintf(fp, "%lf\t%lf\t%lf\t",\
			elementPtr[l].xA[0],elementPtr[l].xA[1],elementPtr[l].xA[2]);
			//half span, chord sweep, dihedral of elemenatry wing
			fprintf(fp, "%lf\t%lf\t",\
			elementPtr[l].eta,elementPtr[l].chord);
			fprintf(fp, "%lf\t%lf\t",\
			elementPtr[l].phi*RtD,elementPtr[l].nu*RtD);
			//normal of elemenatry wing
			fprintf(fp, "%lf\t%lf\t%lf\t",\
			elementPtr[l].normal[0],elementPtr[l].normal[1],elementPtr[l].normal[2]);
			//local free stream velocity at center
			fprintf(fp, "%lf\t%lf\t%lf\n",\
			elementPtr[l].u[0],elementPtr[l].u[1],elementPtr[l].u[2]);
		}
	fclose(fp);
}
//===================================================================//
		//END OF Save_Elementary_Wings
//===================================================================//
//===================================================================//
		//START OF Save_Trailing_Edge
//===================================================================//
void Save_Trailing_Edge(const GENERAL info,const BOUND_VORTEX* trailedgePtr)
{
	//save information of trailing edge to file "Elementary_Wings.txt
	//in directory "output"

	int l;			//loop counter
	FILE *fp;		//output file
	char filename[160];//file path and name

	//opens file for appending
	sprintf(filename,"%s%s",info.output,"Elementary_Wings.txt");
	fp = fopen(filename, "a");

	//writes header
	fprintf(fp,"\n\n\nProgram Version: %s\n",PROGRAM_VERSION);
	fprintf(fp, "\n\t\ttrailing edge elements\n");
	fprintf(fp, "\ttrailing edge ctr\t\t\thalf-span\tsweep\t\tdihedral\n");
	fprintf(fp, "n  \txo\tyo\t\tzo\t\teta\t\tphi\t\tnu\t\tlocal vel\n");

	for(l=0; l<info.nospanelement; l++)
		{
			//trailing edge element index
			fprintf(fp, "%d  ",l);
			//coord. of trailing edge element center
			fprintf(fp, "%lf\t%lf\t%lf\t",\
			trailedgePtr[l].xo[0],trailedgePtr[l].xo[1],trailedgePtr[l].xo[2]);
			//half span, sweep, dihedral of elemenatry wing
			fprintf(fp, "%lf\t%lf\t%lf\t",\
			trailedgePtr[l].eta,trailedgePtr[l].phi*RtD,trailedgePtr[l].nu*RtD);
			//local free stream velocity at center
			fprintf(fp, "%lf\t%lf\t%lf\n",\
			trailedgePtr[l].u[0],trailedgePtr[l].u[1],trailedgePtr[l].u[2]);
		}
	fclose(fp);
}
//===================================================================//
		//END OF Save_Trailing_Edge
//===================================================================//
//===================================================================//
		//START of Horstmann_Results
//===================================================================//
void Horstmann_Results(const GENERAL info,const BOUND_VORTEX* elementPtr,\
					   const double CL,const double CLi,\
					   const double CY,const double CYi,\
					   const double CDi_ellipt,const double CDi_Trefftz,\
					   const double CDi_Eppler)
{
//saves results of Horstmann's method to file

	int i;				//loop counter
	double tempS;
	double cl,cy,cn;  //temporary stores for cl, cy, and cn
	FILE *fp;			//output file
	char filename[160];//file path and name

	//creates file "output\results.txt"
	sprintf(filename,"%s%s",info.output,"results.txt");
	fp = fopen(filename, "w"); 			//###/

	fprintf(fp,"\n\n\nProgram Version: %s\n",PROGRAM_VERSION);
	fprintf(fp,"Results with Horstmann's method, fixed wake\n");

	fprintf(fp,"Ref. Area   : %lf\nAspect Ratio: %lf\n",info.S,info.AR);
	fprintf(fp,"Alpha       : %lf\nBeta        : %lf\n",\
											info.alpha*RtD,info.beta*RtD);
	if (info.sym==1)		fprintf(fp,"symmetrical geometry: 1\n");
	else					fprintf(fp,"asymmetrical geometry: 0\n");

	fprintf(fp,"no. of spanwise element:  %d\n",info.nospanelement);
	fprintf(fp,"no. of chordwise element: %d\n",info.m);

	fprintf(fp,"no. of wings: %d  first and last span indices of each wing:",info.nowing);
	for(i=0;i<info.nowing;i++)
		fprintf(fp,"  %d  %d",info.wing1[i],info.wing2[i]);
	fprintf(fp,"\n");


	fprintf(fp,"\n");
	fprintf(fp,"wing elements:\n");
//	fprintf(fp,"No\txo\t\tyo\t\tzo\t\tchord\t\teta\t\tphi\t\tnu\t\t");
//	fprintf(fp,"A\t\tB\t\tC\t\tcl\t\tcy\t\tcd\t");
//	fprintf(fp,"u\t\tv\t\tw\t\tu_ind\t\tv_ind\t\tw_ind\t#\n");

	fprintf(fp,"%2s %14s %12s %12s %12s %12s %12s %12s",\
				"No","xo","yo","zo","chord","eta","phi","nu");
	fprintf(fp," %12s %12s %12s %12s %12s %12s %12s",\
				"A","B","C","cl","cy","cN","cd");
	fprintf(fp," %12s %12s %12s %12s %12s %12s",\
				"u","v","w","u_ind","v_ind","w_ind");
	fprintf(fp,"  # \n");

	for(i=0;i<info.noelement;i++)
	{
		//local dyn. pressure times element area
		tempS  = 1/dot(elementPtr[i].u,elementPtr[i].u);  //U(y)^2
		tempS /= (elementPtr[i].eta*elementPtr[i].chord);

		fprintf(fp,"%3d %13.5lf %12.5lf %12.5lf %12.5lf %12.5lf %12.5lf %12.5lf",i,\
		elementPtr[i].xo[0],elementPtr[i].xo[1],elementPtr[i].xo[2],\
		elementPtr[i].chord,elementPtr[i].eta,elementPtr[i].phi*RtD,\
		elementPtr[i].nu*RtD);
		
		cl=(elementPtr[i].N_free[0]+elementPtr[i].N_ind[0])*tempS;
		cy=(elementPtr[i].N_free[1]+elementPtr[i].N_ind[1])*tempS;
		cn=sqrt(cl*cl+cy*cy);
			if(cl<0) cn*=-1;

		fprintf(fp," %12.5lf %12.5lf %12.5lf %12.5lf %12.5lf %12.5lf %12.5lf",\
		elementPtr[i].A,elementPtr[i].B,elementPtr[i].C,\
		cl,cy,cn,elementPtr[i].CDind);

		fprintf(fp," %12.5lf %12.5lf %12.5lf %12.5lf %12.5lf %12.5lf",\
		elementPtr[i].u[0],elementPtr[i].u[1],elementPtr[i].u[2],\
		elementPtr[i].uind[0],elementPtr[i].uind[1],elementPtr[i].uind[2]);


		fprintf(fp,"\n");

	}


	fprintf(fp,"\n\n\nTotal Wing loads:\n\n");

	fprintf(fp,"CL=%lf\tCLi=%lf\tCLfree=%lf\nCY=%lf\tCYi=%lf\tCYfree=%lf\n\n",\
		CL,CLi,CL-CLi,CY,CYi,CY-CYi);//#

	fprintf(fp,"CDi_ellipt   = %lf\n",CDi_ellipt);
	fprintf(fp,"CDi_Trefftz  = %lf\n",CDi_Trefftz);
	fprintf(fp,"CDi_Eppler   = %lf\n\n",CDi_Eppler);

	fprintf(fp,"e_Trefftz  = %lf  k_Trefftz  = %lf\n",
						CDi_ellipt/CDi_Trefftz,CDi_Trefftz/CDi_ellipt);
	fprintf(fp,"e_Eppler   = %lf  k_Eppler   = %lf\n",\
						CDi_ellipt/CDi_Eppler,CDi_Eppler/CDi_ellipt);

	//results of the time-stepping method to file
	//header
	fprintf(fp,"\n\n\nProgram Version: %s\n",PROGRAM_VERSION);
	fprintf(fp,"Results with time-stepping method and with a ");
	if(info.relax == 1) 		fprintf(fp,"relaxed wake\n");
	else						fprintf(fp,"fixed wake\n");
	if(info.steady ==1)			fprintf(fp,"steady aerodynamics\n\n");
	else						fprintf(fp,"unsteady aerodynamics\n\n");

	fprintf(fp,"max. no. of timesteps 	 : %d\n",info.maxtime);
	fprintf(fp,"time increment        	 : %lf\n",info.deltime);
	fprintf(fp,"step size (based on Uinf): %lf\n\n",info.deltime*info.Uinf);

	fprintf(fp,"%4s %16s %16s %16s %16s %16s %16s %16s %16s\t#\n",\
	"step","time","CD_ind","CL_tot","CL_ind","CY_tot","CY_ind","CD_ellipt","e   ");



fclose(fp);			//###/
}
//===================================================================//
		//END of Horstmann_Results
//===================================================================//

//===================================================================//
		//START of Horstmann_Results
//===================================================================//
void Header(const GENERAL info,const BOUND_VORTEX* elementPtr,\
					   const double CL,const double CLi,\
					   const double CY,const double CYi,\
					   const double CDi_ellipt,const double CDi_Trefftz,\
					   const double CDi_Eppler)
{
//saves results of Horstmann's method to file

	int i;				//loop counter
//	double tempS;
	FILE *fp;			//output file
	char filename[160];//file path and name

	//creates file "output\results.txt"
	sprintf(filename,"%s%s",info.output,"results.txt");
	fp = fopen(filename, "w"); 			//###/

	fprintf(fp,"\n\n\nProgram Version: %s\n",PROGRAM_VERSION);

	fprintf(fp,"Ref. Area   : %lf\nAspect Ratio: %lf\n",info.S,info.AR);
	fprintf(fp,"Alpha       : %lf\nBeta        : %lf\n",\
											info.alpha*RtD,info.beta*RtD);
	fprintf(fp,"density:  %lf\n",info.density);

	if (info.sym==1)		fprintf(fp,"symmetrical geometry: 1\n");
	else					fprintf(fp,"asymmetrical geometry: 0\n");

	fprintf(fp,"no. of spanwise element:  %d\n",info.nospanelement);
	fprintf(fp,"no. of chordwise element: %d\n",info.m);

	fprintf(fp,"no. of wings: %d  first and last span indices of each wing:",info.nowing);
	for(i=0;i<info.nowing;i++)
		fprintf(fp,"  %d  %d",info.wing1[i],info.wing2[i]);
	fprintf(fp,"\n");


	//results of the time-stepping method to file
	//header
	fprintf(fp,"\n\n\nProgram Version: %s\n",PROGRAM_VERSION);
	fprintf(fp,"Results with time-stepping method and with a ");
	if(info.relax == 1) 		fprintf(fp,"relaxed wake\n");
	else						fprintf(fp,"fixed wake\n");
	if(info.steady ==1)			fprintf(fp,"steady aerodynamics\n\n");
	else						fprintf(fp,"unsteady aerodynamics\n\n");

	fprintf(fp,"max. no. of timesteps 	 : %d\n",info.maxtime);
	fprintf(fp,"time increment        	 : %lf\n",info.deltime);
	fprintf(fp,"step size (based on Uinf): %lf\n\n",info.deltime*info.Uinf);

	fprintf(fp,"%4s %16s %16s %16s %16s %16s %16s %16s %16s",\
	"step","time","CD_ind","CL_tot","CL_ind","CY_tot","CY_ind","CD_ellipt","e   ");

	fprintf(fp,"%16s %16s %16s %16s %16s %16s",\
	"  Fx","  Fy","  Fz","  Mx","  My","  Mz");

	fprintf(fp,"\t#\n");

fclose(fp);			//###/
}
//===================================================================//
		//END of Header
//===================================================================//

//===================================================================//
		//START of Time_Stepping_Results
//===================================================================//
void Time_Stepping_Results(const GENERAL info,int const first, int const last,\
						   double **CN,double *CDi,\
						   const double CL_finit,const double CLi_finit,\
						   const double CY_finit,const double CYi_finit,\
						   const double CDi_finit)
{
//saves results of the time-stepping method to file

	int i;				//loop counter
	double time,CDi_ellipt,e;
	FILE *fp;			//output file
	char filename[160];//filepath and name

	//creates file "output\results.txt"
	sprintf(filename,"%s%s",info.output,"results.txt");
	fp = fopen(filename, "a"); 			//###/

/*	//header
	fprintf(fp,"\n\n\nProgram Version: %s\n",PROGRAM_VERSION);
	fprintf(fp,"Results with time-stepping method and with a ");
	if(info.relax == 1) 		fprintf(fp,"relaxed wake\n");
	else						fprintf(fp,"fixed wake\n");
	if(info.steady ==1)			fprintf(fp,"steady aerodynamics\n\n");
	else						fprintf(fp,"unsteady aerodynamics\n\n");

	fprintf(fp,"max. no. of timesteps 	 : %d\n",info.maxtime);
	fprintf(fp,"time increment        	 : %lf\n",info.deltime);
	fprintf(fp,"step size (based on Uinf): %lf\n\n",info.deltime*info.Uinf);

	fprintf(fp,"%4s%16s%16s%16s%16s%16s%16s%16s%16s\n",\
	"step","time","CD_ind","CL_tot","CL_ind","CY_tot","CY_ind","CD_ellipt","e");
*/
	for (i=first; i<=last; i++)
	{
	  time   		= info.deltime*(1+i);
	  CDi_ellipt 	= (CN[i][0]*CN[i][0]+CN[i][2]*CN[i][2])/(info.AR*Pi);
	  e 	   		= CDi_ellipt/CDi[i];
	fprintf(fp,"%4d %16.4lf %16.8lf %16.8lf %16.8lf %16.8lf %16.8lf %16.8lf %16.8lf",\
	  		  i,time,CDi[i],CN[i][0],CN[i][1],CN[i][2],CN[i][3],CDi_ellipt,e);

//	fprintf(fp,"%16.8lf",CN[i][0]*CN[i][0]/(info.AR*Pi*CDi[i]));
	fprintf(fp,"\n");
	}

fclose(fp);			//###/
}
//===================================================================//
		//END of Time_Stepping_Results
//===================================================================//
//===================================================================//
		//START of Time_Stepping_End_Results
//===================================================================//
void Time_Stepping_End_Results(const GENERAL info,const int steps,\
							const int timestep,double **CN,double *CDi,\
						   const double CL_finit,const double CLi_finit,\
						   const double CY_finit,const double CYi_finit,\
						   const double CDi_finit)
{
//saves results of the time-stepping method to file

	int i,first;				//loop counter
	double time,CDi_ellipt,e;
	FILE *fp;			//output file
	char filename[160];//file path and name

	//creates file "output\results.txt"
	sprintf(filename,"%s%s",info.output,"results.txt");
	fp = fopen(filename, "a"); 			//###/

	first=int(timestep/steps+0.5);

	for(i=first;i<=timestep;i++)
	{
		time   		= info.deltime*(1+i);
		CDi_ellipt 	= (CN[i][0]*CN[i][0]+CN[i][2]*CN[i][2])/(info.AR*Pi);
		e 	   		= CDi_ellipt/CDi[i];
		fprintf(fp,"%4d %16.4lf %16.8lf %16.8lf %16.8lf",\
					i,time,CDi[i],CN[i][0],CN[i][1]);
		fprintf(fp," %16.8lf %16.8lf %16.8lf %16.8lf",\
	  		  		CN[i][2],CN[i][3],CDi_ellipt,e);
		fprintf(fp,"\n");
	}

	fprintf(fp,"\nAnd that's it!!\n");

fclose(fp);			//###/
}
//===================================================================//
		//END of Time_Stepping_End_Results
//===================================================================//

//===================================================================//
		//START of Save_Surface_DVEs
//===================================================================//
void Save_Surface_DVEs(const GENERAL info,const DVE *surfacePtr)
{
//save information of surface DVEs to file

	int l;			//loop counter
	FILE *fp;		//output file
	char filename[160];//file path and name

	//creates file "output\Surface_DVE.txt"
	sprintf(filename,"%s%s",info.output,"Surface_DVE.txt");
	fp = fopen(filename, "w");

	//writes header
	fprintf(fp,"\n\n\nProgram Version: %s\n",PROGRAM_VERSION);
	fprintf(fp, "\n\nsurface DVE\n");
	fprintf(fp, "\tcontrol point\thalf-span  sweep  dihedral\n");
	fprintf(fp, "n\txo\t\tyo\t\tzo\t\teta\t\txsi\t\tnu");
	fprintf(fp, "\t\tphiLE\t\tphiTE\t\tepsilon\t\tnormal\n");

	for(l=0; l<info.noelement; l++)
		{
			//elementary wing number
			fprintf(fp, "%d\t",l);
			//coord. of DVE center
			fprintf(fp, "%lf\t%lf\t%lf\t",\
			surfacePtr[l].xo[0],surfacePtr[l].xo[1],surfacePtr[l].xo[2]);
			//half span, sweep, dihedral of elemenatry wing
			fprintf(fp, "%lf\t%lf\t%lf\t",\
			surfacePtr[l].eta,surfacePtr[l].xsi,surfacePtr[l].nu*RtD);
			//leading and trailing edge sweeps, pitch
			fprintf(fp, "%lf\t%lf\t%lf\t",\
			surfacePtr[l].phiLE*RtD,surfacePtr[l].phiTE*RtD,surfacePtr[l].epsilon*RtD);
			//surface normal
			fprintf(fp, "%lf\t%lf\t%lf\t",\
			surfacePtr[l].normal[0],surfacePtr[l].normal[1],surfacePtr[l].normal[2]);

			fprintf(fp, "\n");
		}
	fclose(fp);
}
//===================================================================//
		//END of Save_Surface_DVEs
//===================================================================//

//===================================================================//
        //START of Save_Input_File
//===================================================================//
//save input file to output directory
void Save_Input_File(const char sourcefile[126],const char targetdir[126])
{
    //copies the input file from input directory to output directory
    //adapted from https://www.programmingsimplified.com/c-program-copy-file
    //sourcefile  - input file for FreeWake
    //targetdrive - output directory
    
    FILE *source,*target;
    char ch,targetfile[126];

    source = fopen(sourcefile, "r");
    printf("source = %s\n", sourcefile);

    if (source == NULL)
    {
    	printf("ERROR IS IN LINE 494 [write_output.cpp]\n\n");
       	printf(" couldn't copy input file\nPress any key to exit...\n");
		exit(EXIT_FAILURE);
    }

    // WB - FEB 13, 2025 - added the backslashes
    char path_separator = '/';
    #ifdef _WIN32
    	path_separator = '\\';
    #endif
    sprintf(targetfile, "%s%c%s",targetdir, path_separator, sourcefile);
    
    //sprintf(targetfile, "%s%s",targetdir,sourcefile);
    //printf("Attempting to write to: %s\n", targetfile);
	//printf("targetdir = %s\n", targetdir);
	//printf("targetfile = %s\n", targetfile);

    target = fopen(targetfile, "w");

    if (target == NULL)
    {
       fclose(source);
       printf("ERROR IS IN LINE 517 [write_output.cpp]\n\n");
       printf(" couldn't copy input file\nPress any key to exit...\n");
       
       exit(EXIT_FAILURE);
    }

    while ((ch = fgetc(source)) != EOF)
       fputc(ch, target);

    fclose(source);
    fclose(target);
}
//===================================================================//
        //END of Save_Input_File
//===================================================================//

//===================================================================//
        //START of Save_Config_Head_File
//===================================================================//
//saves input file and header to configuration file in output directory
void Save_Config_Head_File(GENERAL info)
{
    //copies the input file from input directory to config file located in
    //output directory
    //adapted from https://www.programmingsimplified.com/c-program-copy-file
    //sourcefile  - input file for FreeWake
    //targetfile - configuration file (with directory path)
    
    FILE *source,*target;
    char ch,str[100];

    source = fopen(info.inputfilename, "r");

    if (source == NULL)
    {
       printf(" couldn't copy input file\nPress any key to exit...\n");
       exit(EXIT_FAILURE);
    }

    //If configuration file does not exist yet, it is being created
    if ((target = fopen(info.config, "r"))== NULL)
    {
        printf("Configuration file %s is created\n",info.output);
  
        target = fopen(info.config, "w"); //create file to write to

        if (target == NULL)
        {
            fclose(target);
            printf(" couldn't copy target file\nPress any key to exit...\n");
            exit(EXIT_FAILURE);
        }

        while ((ch = fgetc(source)) != EOF)
            fputc(ch, target);

        fclose(source);

        fprintf(target,"\n\nThe next line is a special marker\n#$&#$&\n\n"); //sets a marker
        
        //writing general information heading
        fprintf(target,"Ref. length for rolling mom.  = %lf\n",info.b);
        fprintf(target,"Ref. lenght for pitching mom. = %lf\n",info.cmac);
        fprintf(target,"Ref. length for yawing mom.   = %lf\n",info.b);
        fprintf(target,"Reference area                = %lf\n",info.S);
     //   fprintf(target,"Projected planform area       = %lf\n",info.projAREA);
     //   fprintf(target,"Unraveled area                = %lf\n",info.AREA);
        fprintf(target,"Referece span                 = %lf\n",info.b);
     //   fprintf(target,"Projected span                = %lf\n",info.projSPAN);
        fprintf(target,"Reference aspect ratio        = %lf\n",info.AR);
        fprintf(target,"number of wings               = %d\n",info.nowing);
        fprintf(target,"saved no. of flight conditns  = ");
        if (info.trimCL == 0)	fprintf(target,"\t\t<-- alpha sweep performed, value without meaning\n\n");
        else 					fprintf(target,"\n\n");
        
        //header of tabulated case data
        fprintf(target,"%-10s%-10s%-10s%-10s","Flight","CLtarget","Alpha","Beta");
        fprintf(target,"%-12s%-12s%-12s","CL","CQ","CDI");
        fprintf(target,"%-12s%-12s%-12s","CX","CY","CZ");
        fprintf(target,"%-12s%-12s%-12s\n","CL","CM","CN");
        fprintf(target,"%-10s%-10s%-10s%-10s","Condtn","deflt=10","[deg]","[deg]");
        fprintf(target,"%-12s%-12s%-12s","(lift)","(side)"," ");
        fprintf(target,"%-12s%-12s%-12s"," "," "," ");
        fprintf(target,"%-12s%-12s%-12s\n","(roll)"," "," ");
        fprintf(target,"-------------------------------------------------------------");
        fprintf(target,"--------------------------");
        fprintf(target,"------------------------------------------------------------#\n");

        fclose(target);
    }
    else  //configuration file already exists, determine flight config. number
    {
        printf("Configuration file %s already exists\n",info.config);
    }

}
//===================================================================//
        //END of Save_Config_Head_File
//===================================================================//


//===================================================================//
		//START of Save_Timestep
//===================================================================//
void Save_Timestep(const GENERAL info,const int timestep,DVE **wakePtr,\
								const DVE *surfacPtr,double **N_force)
{
//saves results of current timestep to file

	int time,span;		//loop counter
	FILE *fp;			//output file
	char filename[160];//file path and name

	//creates file name timestep##.txt ## is number of timestep
	sprintf(filename,"%s%s%d%s",info.output,"timestep",timestep,".txt");

	//creates file in subdirectory output
	fp = fopen(filename, "w");

	//writes header
	fprintf(fp,"\n\n\nProgram Version: %s\n",PROGRAM_VERSION);
	if(info.relax == 1) 		fprintf(fp,"relaxed wake\n");
	else						fprintf(fp,"fixed wake\n");
	if(info.steady ==1)			fprintf(fp,"steady aerodynamics\n");
	else						fprintf(fp,"unsteady aerodynamics\n");

	fprintf(fp,"Ref. Area   : %lf\nAspect Ratio: %lf\n",info.S,info.AR);
	fprintf(fp,"Alpha       : %lf\nBeta        : %lf\n",\
											info.alpha*RtD,info.beta*RtD);

	if (info.sym==1)			fprintf(fp,"symmetrical geometry: 1\n");
	else						fprintf(fp,"asymmetrical geometry: 0\n");

	fprintf(fp,"lifting surface after timestep: %d\n",timestep);
	fprintf(fp,"elements in span direction: %d\n",info.nospanelement);
	fprintf(fp,"number of surface elements: %d\n",info.noelement);  //updated GB 2-25-20

	fprintf(fp,"no. of wings: %d  first and last span indices of each wing:",info.nowing);
	for(span=0;span<info.nowing;span++)
		fprintf(fp,"  %d  %d",info.wing1[span],info.wing2[span]);
	fprintf(fp,"\n");
	fprintf(fp,"\n");

	//writes header for information on surface elements
	fprintf(fp,"%6s %16s %16s %16s %16s %16s %16s %16s",\
	"index","xo","yo","zo","Nlift_tot","NLift_ind","NY_tot","NYi");
	fprintf(fp," %16s %16s %16s %16s %16s %16s %16s %16s %16s %16s %16s\t#\n",\
	"A","B","C","eta","xsi","nu","epsilon","psi","phiLE","phi0","phiTE");

	for(span=0;span<info.noelement;span++)
	{
		//surface element index
		fprintf(fp,"%6d",span);
		//coord. of ref point
		fprintf(fp," %16.12lf %16.12lf %16.12lf",\
				surfacePtr[span].xo[0],surfacePtr[span].xo[1],\
				surfacePtr[span].xo[2]);
		//normal forces per density
		fprintf(fp," %16.12lf %16.12lf %16.12lf %16.12lf",\
				N_force[span][0]+N_force[span][1],N_force[span][1],\
				N_force[span][2]+N_force[span][3],N_force[span][3]);
		//more info on element
		fprintf(fp," %16.12lf %16.12lf %16.12lf %16.12lf %16.12lf",\
				surfacePtr[span].A,surfacePtr[span].B,surfacePtr[span].C,\
				surfacePtr[span].eta,surfacePtr[span].xsi);
		fprintf(fp," %16.12lf %16.12lf %16.12lf",\
				surfacePtr[span].nu*RtD,surfacePtr[span].epsilon*RtD,\
				surfacePtr[span].psi*RtD);
		fprintf(fp," %16.12lf %16.12lf %16.12lf",surfacePtr[span].phiLE*RtD,\
				(surfacePtr[span].phiLE+surfacePtr[span].phiTE)*0.5*RtD,\
				surfacePtr[span].phiTE*RtD);
		fprintf(fp,"\n");
	}


	//writes header for wake information
	fprintf(fp,"\n\nwake shape after timestep: %d\n",timestep);
	fprintf(fp,"%5s%5s%16s %16s %16s %16s %16s %16s",\
				"span","time","xo","yo","zo","nu","epsilon","psi");
	fprintf(fp," %16s %16s %16s %16s %16s %16s %16s %16s %16s",\
				"U","V","W","A","B","C","eta","xsi","singfct");
	fprintf(fp," %16s %16s %16s\t#\n","phiLE","phi0","phiTE");

	for (time=0;time<=timestep;time++)
	{
		//loop across wake elements of one time/downstream location
		for(span=0;span<info.nospanelement;span++)
		{
			//trailing edge element index
			fprintf(fp,"%5d%5d",span,time);
			//coord. of ref point
			fprintf(fp," %16.12lf %16.12lf %16.12lf",\
						wakePtr[time][span].xo[0],wakePtr[time][span].xo[1],\
						wakePtr[time][span].xo[2]);
			//nu,epsilon, sweep, dihedral of elemenatry wing
			fprintf(fp," %16.12lf %16.12lf %16.12lf",\
						wakePtr[time][span].nu*RtD,\
						wakePtr[time][span].epsilon*RtD,\
						wakePtr[time][span].psi*RtD);
			//local free stream velocity at center
			fprintf(fp," %16.12lf %16.12lf %16.12lf",\
						wakePtr[time][span].u[0],wakePtr[time][span].u[1],\
						wakePtr[time][span].u[2]);
			fprintf(fp," %16.12lf %16.12lf %16.12lf",\
						wakePtr[time][span].A,wakePtr[time][span].B,\
						wakePtr[time][span].C);
			//element half span and half chord
			fprintf(fp," %16.12lf %16.12lf",\
						wakePtr[time][span].eta,wakePtr[time][span].xsi);
			//element half span and half chord
            fprintf(fp," %16.12lf",wakePtr[time][span].singfct);//    999.999);
	//					wakePtr[time][span].K);//changed GB 2/6/20
			//leading-edge, mid-chord, and trailing edge sweep2
			fprintf(fp," %16.12lf %16.12lf %16.12lf",\
				wakePtr[time][span].phiLE*RtD,\
				wakePtr[time][span].phi0*RtD,wakePtr[time][span].phiTE*RtD);

//			fprintf(fp," %16.12lf",wakePtr[time][span].singfct);

//			//left edge half chord and factor for treating the tip singularity
//			fprintf(fp," %16.12lf %16.12lf %16.12lf %16.12lf",\
//					wakePtr[time][span].xleft[0],wakePtr[time][span].xleft[1],\
//					wakePtr[time][span].xleft[2],wakePtr[time][span].singfct);

			fprintf(fp,"\n");
		}
		fprintf(fp,"\n");
	}
fclose(fp);

//	for (i=0; i<info.maxtime; i++)
//	{
//	  time   		= info.deltime*(1+i);
//	  CDi_ellipt 	= (CN[i][0]*CN[i][0]+CN[i][2]*CN[i][2])/(info.AR*Pi);
//	  e 	   		= CDi_ellipt/CDi[i];
//	fprintf(fp,"%4d%8.4lf%16.5lf%16.5lf%16.5lf%16.5lf%16.5lf%16.5lf%16.5lf\n",\
//	  		  i,time,CDi[i],CN[i][0],CN[i][1],CN[i][2],CN[i][3],CDi_ellipt,e);
//	}

}
//===================================================================//
		//Ende of Save_Timestep
//===================================================================//

//===================================================================//
		//START of Save_SurfaceDVE_Loads
//===================================================================//
void Save_SurfaceDVE_Loads(const GENERAL info,const int timestep,\
						const DVE *surfacPtr)
{
//saves current timestep forces and moments of surface DVE's

	int time,span;		//loop counter
	FILE *fp;			//output file
	char filename[160];//file path and name

	//creates file name timestep##.txt ## is number of timestep
	sprintf(filename,"%s%s%d%s",info.output,"SDVE_loads",timestep,".txt");

	//creates file in subdirectory output
	fp = fopen(filename, "w");

	//writes header
	fprintf(fp,"\n\n\nProgram Version: %s\n",PROGRAM_VERSION);
	if(info.relax == 1) 		fprintf(fp,"relaxed wake\n");
	else						fprintf(fp,"fixed wake\n");
	if(info.steady ==1)			fprintf(fp,"steady aerodynamics\n");
	else						fprintf(fp,"unsteady aerodynamics\n");

	fprintf(fp,"Ref. Area   : %lf\nAspect Ratio: %lf\n",info.S,info.AR);
	fprintf(fp,"Alpha       : %lf\nBeta        : %lf\n",\
											info.alpha*RtD,info.beta*RtD);

	if (info.sym==1)			fprintf(fp,"symmetrical geometry: 1\n");
	else						fprintf(fp,"asymmetrical geometry: 0\n");

	fprintf(fp,"lifting surface after timestep: %d\n",timestep);
	fprintf(fp,"elements in span direction: %d\n",info.nospanelement);
	fprintf(fp,"elements in chord direction: %d\n",info.m);

	fprintf(fp,"no. of wings: %d  first and last span indices of each wing:",info.nowing);
	for(span=0;span<info.nowing;span++)
		fprintf(fp,"  %d  %d",info.wing1[span],info.wing2[span]);
	fprintf(fp,"\n");
	fprintf(fp,"\n");

	//writes header for information on surface elements
	fprintf(fp,"%6s %16s %16s %16s %16s %16s %16s",\
			"index","xo","yo","zo","xoLE","yoLE","zoLE");
	fprintf(fp," %16s %16s %16s %16s %16s %16s",\
				"Fx","Fy","Fz","Mx","My","Mz");
	fprintf(fp," %16s %16s %16s %16s %16s %16s\t#\n",\
				"xsi","eta","nu","epsilon","psi","phi0");

	for(span=0;span<info.noelement;span++)
	{
		//surface element index
		fprintf(fp,"%6d",span);
		//coord. of ref point
		fprintf(fp," %16.12lf %16.12lf %16.12lf",\
				surfacePtr[span].xo[0],surfacePtr[span].xo[1],\
				surfacePtr[span].xo[2]);
		//coord. of LE center
		fprintf(fp," %16.12lf %16.12lf %16.12lf",\
				0.5*(surfacePtr[span].x1[0]+surfacePtr[span].x2[0]),\
				0.5*(surfacePtr[span].x1[1]+surfacePtr[span].x2[1]),\
				0.5*(surfacePtr[span].x1[2]+surfacePtr[span].x2[2]));
		//Force
//		fprintf(fp," %16.12lf %16.12lf %16.12lf",\
//				surfacePtr[span].Force[0],\
//				surfacePtr[span].Force[1],\
//				surfacePtr[span].Force[2]);
		//Moment
//		fprintf(fp," %16.12lf %16.12lf %16.12lf",\
 //                surfacePtr[span].Moment[0],\
//				surfacePtr[span].Moment[1],\
//				surfacePtr[span].Moment[2]);
		//more info on element
		fprintf(fp," %16.12lf %16.12lf",\
				surfacePtr[span].xsi,surfacePtr[span].eta);
		fprintf(fp," %16.12lf %16.12lf %16.12lf %16.12lf",\
				surfacePtr[span].nu*RtD,surfacePtr[span].epsilon*RtD,\
				surfacePtr[span].psi*RtD,surfacePtr[span].phi0*RtD);
		fprintf(fp,"\n");
	}


fprintf(fp,"\n");
fclose(fp);


}
//===================================================================//
		//Ende of Save_SurfaceDVE_Loads
//===================================================================//

//===================================================================//
		//START of Test
//===================================================================//
//void Test(const GENERAL,double **);
void Test(const GENERAL info,double **D,const double *R)
{
 //  ###########################################################
 //save D matrix and resultant vector R in file D_matrix.txt
 int m,n;
char filename[160];//file path and name
 FILE *fp;

	//creates file name timestep##.txt ## is number of timestep
	sprintf(filename,"%s%s",info.output,"test.txt");

	 fp = fopen(filename, "a");
	 //writes header line
	 fprintf(fp, "\t");
	 for(m=0; m<info.Dsize; m++)
		 fprintf(fp, "%d\t",m);
	 fprintf(fp, "\t\tR");

	for(n=0; n<info.Dsize; n++)
 	{
		 //row number
 		fprintf(fp, "\n%d\t",n);
 		//n-th row of D
 		for(m=0; m<info.Dsize; m++)
 		fprintf(fp, "%lf\t",D[n][m]);
 		//n-th element of R
		fprintf(fp, "\t\t%lf",R[n]);
 	}
 	fclose(fp);

}
 //###########################################################//*/
//===================================================================//
		//Ende of Test
//===================================================================//


//===================================================================//
		//START of CreateQuiverFile
//===================================================================//
void CreateQuiverFile(const double pos[3], const double vec[3], const int idx, const int filenum)
{
	// This function create the quiver file used for Main_PlotQuiverPY.cpp
	// 	The vector (vec) will be added to the quiver.txt
	//
	// Function inputs:
	//		pos[3]		- Positional values of x,y,z where vector should begin
	//		vec[3]		- Vector values of u,v,w (ie component of the vector to be plotted)
	//		idx 		- Determine if adding or starting the file
	//						idx == 0 opens file and begins new file
	//						idx != 0 adds to file  
	//		filenum		- Number to add to filename (ex. quiver4.txt). Zero will not add any value
	// 
	// Example of use to plot eN
	//	for(i=0;i<info.noelements;i++){
	//	if(i==0){CreateQuiverFile(surfacePtr[l].xo, eN,0,0);}
	//	else{CreateQuiverFile(surfacePtr[l].xo, eN,1,0);}
	//	}
	//
	//
	// D.F.B. in Braunschweig, Germany, Mar. 2020


	FILE *fp;		//output file
	char filename[160];//file path and name 

	if (filenum == 0) 
	{
		// Either open file or add to file
		if (idx == 0) { fp = fopen("output/quiver.txt", "w"); }
		else { fp = fopen("output/quiver.txt", "a"); }
	}
	else 
	{
		sprintf(filename, "%s%d%s","output/quiver", filenum, ".txt");
		if (idx == 0) { fp = fopen(filename, "w"); }
		else { fp = fopen(filename, "a"); }
	}

if (idx == 2) fprintf(fp, "\n");
// Write vector position and vector components to file
fprintf(fp,"%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t\n",pos[0],pos[1],pos[2],vec[0],vec[1],vec[2]);

//Close file
fclose(fp);

}
//===================================================================//
		//CreateQuiverFile
//===================================================================//


//===================================================================//
		//START of SaveSpanInfo
//===================================================================//
void SaveSpanDVEInfo(PANEL *panelPtr,DVE *surfacePtr,DVE **wakePtr,STRIP *spanPtr,\
					double **N_force,const int FCno,const int timestep)
{
	// This function saves for a particular flight condition:
	//			A. Saving information of spanwise strips 
	//			B. Saves wing and wake of last timestep (similar to timestep)
	//			C. Save surface DVEs
	// (i.e. anlge of attack or CL) to a file. 
	// The file name is <input>span<FCno>.txt
	// GB 11-24-20
	//
	// Function inputs:
	//		panelPtr	panel information
	//		surfacePtr 	surface DVE information
	//		wakePtr		wake DVE information
	//		spanPtr		strip information	
	// 		N_force;	surface DVE's normal forces/density
					//[0]: free stream lift, [1]: induced lift,
					//[2]: free stream side, [3]: induced side force/density
					//[4]: free str. normal, [5]: ind. normal frc/density
               		//[6,7,8]: eN_x, eN_y, eN_z in global ref. frame
	//		
	//		FCno		flight condition number (identifies CLtarget or alpha)
	//		timestep 	number of timesteps that were computed
	//
	// Before saving data to file, several values still have to be determined.
	// These values are part of the STRIP structure [where value is computed]:
	//						
	// xref[3],x1[3],x2[3]; //reference point of strip, left and right points [here]
	// Cf[3];       	//strip force coefficients Cfx,Cfy,Cfz (includes ind. drag) [here]
	// Cn;				//strip normal force coef. (w/o drag), used for cd_profie [longtrim.cpp]
	// Cm[3];			//strip moments coefficients Cl, Cm, Cn [lift_force.cpp]
	// area,chord,span; //reference area, chord and span of strip [here]
	// momarm[3];		//reference lengths for Cl, Cm and Cn [here]
	// chord1,chord2;	//strip chord length at left and right edge [here]
	// A,B,C;			//strip circulation coefficients [here]
	// Gamma1,Gamma0,Gamma2;//circulation at left, center and right of strip [here]
	// cn1,cn0,cn2;		//section force coeff. (VxGamma) left, center and right of strip [here]

	// Span_force[3]; 	//x,y,z aerodynamic force/density (includes drag) [lift_force.cpp]
	// Moment[3];		//x,y,z moment about strip ref. point [lift_force.cpp]
	// D_force;			//induced drag forces/density [lift_force.cpp]
	// Cd;				//ind. drag coefficient of strip, [longtrim.cpp]
	// cd_profile;		//strip profile drag coefficient [MainPerfCode]

	//
	//===================================================================//

    int span,index,panel,n,m,time;	//loop counter
    unsigned long length=strlen(info.config); //lenght of configuration file name with path
    double tempS;
    double eps, psi, nu,tempA[3],tempAA[3]; //required to rotate geometry back into level
    double delC,delSpan;	//change in chord across a panel
	char *FCfile,*Ftstep,*SDVE;	//file path and name for flight condition and last timestep 
    FILE *fp;		//output file

	//===================================================================//
    //create file names 
	ALLOC1D(&FCfile,length+11);  //allocate memory for file name, max. 999 flight conditions
	ALLOC1D(&Ftstep,length+11);  //allocate memory for file name, max. 999 flight conditions
	ALLOC1D(&SDVE,length+11);	 //allocate memory for file name, max. 999 flight conditions

	//extracting core of filenames, i.e. path + input filename without .txt
	n=0;
	m=int(length)-7;
	while(n<m)
	{
		SDVE[n]=info.config[n];
		n++;
	}
	SDVE[n]='\0';

	//filename is outputpath/input+"DVE#"+<FCno>+".txt"
	sprintf(Ftstep,"%s%s%d%s",SDVE,"TDVE#",FCno,".txt");
	//filename is outputpath/input+"FC#"+<FCno>+".txt"
	sprintf(FCfile,"%s%s%d%s",SDVE,"FC#",FCno,".txt");
	//filename is outputpath/input+"SDVE#"+<FCno>+".txt"
	sprintf(SDVE,"%s%s%d%s",SDVE,"SDVE#",FCno,".txt");
	
    //===================================================================//
		//START A. Saving information of spanwise strips 
    //===================================================================//
	//===================================================================//
    //assemble missing data for each strip

	//determining missing information of spanwise strips
	span=0;  //index along wingspan
	index=0; //index of DVE at TE

	if(info.flagHORZ)   eps = -info.alphaset;
    else eps = 0;
    psi = -info.betaset;
    nu = -info.bank;

    //loop over panels
    for(panel=0;panel<info.nopanel;panel++)
    {
    	//computing chords of strips, change in chord across panel span
    	delC = (panelPtr[panel].c2-panelPtr[panel].c1)/panelPtr[panel].n;

    	//computing the span of the strips, which is the panel span/n
    		tempA [0]=panelPtr[panel].x1[0]-panelPtr[panel].x2[0];
    		tempA [1]=panelPtr[panel].x1[1]-panelPtr[panel].x2[1];
    		tempA [2]=panelPtr[panel].x1[2]-panelPtr[panel].x2[2];
    	delSpan = norm2(tempA)/panelPtr[panel].n;

		for (n = panelPtr[panel].LE1; n <= panelPtr[panel].LE2; n++)
		{
			index = n+panelPtr[panel].n*(panelPtr[panel].m-1); //index of DVE at TE			
			
			//chords left, center and right of strip
    		spanPtr[span].chord1 = panelPtr[panel].c1+delC*(n-panelPtr[panel].LE1);
    		spanPtr[span].chord = spanPtr[span].chord1+delC*0.5;
			spanPtr[span].chord2 = spanPtr[span].chord1+delC;

			//computing the span and area of each strip, 
			spanPtr[span].span = delSpan;
			spanPtr[span].area = delSpan * spanPtr[span].chord;

			//edge points of strips, x1 and x2, were saved in PitchmMoment.cpp

			//rotation x1 and x2 back to wings-level reference frame, 
			//reversed process of Panel_Rotation in wing_geometry.cpp
			//rotate X1
			rotateX(spanPtr[span].x1,-nu,tempA);
			rotateY(tempA, -eps,tempAA);
			rotateZ(tempAA,psi,spanPtr[span].x1); //%note! this is beta, NOT yaw for circ. flight. 
			//rotate X2
			rotateX(spanPtr[span].x2,-nu,tempA);
			rotateY(tempA,-eps, tempAA);
			rotateZ(tempAA,psi,spanPtr[span].x2); //%note! this is beta, NOT yaw for circ. flight. 

			//reference point is at the center of the leading edge of each strip
			spanPtr[span].xref[0] = (spanPtr[span].x1[0] + spanPtr[span].x2[0]) / 2;
			spanPtr[span].xref[1] = (spanPtr[span].x1[1] + spanPtr[span].x2[1]) / 2;
			spanPtr[span].xref[2] = (spanPtr[span].x1[2] + spanPtr[span].x2[2]) / 2;

			//circulation values of strip are based on TE values
			spanPtr[span].A = surfacePtr[index].A;
			spanPtr[span].B = surfacePtr[index].B;
			spanPtr[span].C = surfacePtr[index].C;

			spanPtr[span].Gamma1 = surfacePtr[index].A - surfacePtr[index].eta*surfacePtr[index].B \
								 + surfacePtr[index].eta*surfacePtr[index].eta*surfacePtr[index].C;
			spanPtr[span].Gamma0 = surfacePtr[index].A;
			spanPtr[span].Gamma2 = surfacePtr[index].A + surfacePtr[index].eta*surfacePtr[index].B \
								 + surfacePtr[index].eta*surfacePtr[index].eta*surfacePtr[index].C;
			
			//moment arm for computing strip moment coefficients about ref. point.
			spanPtr[span].momarm[0] = spanPtr[span].span;
			spanPtr[span].momarm[1] = spanPtr[span].chord;
			spanPtr[span].momarm[2] = spanPtr[span].span;

			//compute section force coefficients by doing VXGamma at 1,0, and 2
			//
			//since the section normal force coefficients are mainly used for structural anlaysis
			//a simplified approach is used, that is the section force is based on Kutta-Joukowski
			//theorem.
			tempS = 2/norm2(surfacePtr[index].u);

			spanPtr[span].cn1 = spanPtr[span].Gamma1*tempS/spanPtr[span].chord1;
			spanPtr[span].cn0 = spanPtr[span].Gamma0*tempS/spanPtr[span].chord;
			spanPtr[span].cn2 = spanPtr[span].Gamma2*tempS/spanPtr[span].chord2;

	    	span++; // increase to next span index
  	    }//loop over span (n) of panel
 	} //next panel

	//===================================================================//
    //section load coefficients
    span=0;
    for(panel=0;panel<info.nopanel;panel++)
    {
       //loop over panel span (along leading edge indices)
       for(n=panelPtr[panel].LE1;n<=panelPtr[panel].LE2;n++)
       {
       		//1/(0.5 U^2 S) of chordwise strip
         	tempS = 2/(dot(surfacePtr[n].u,surfacePtr[n].u)*\
                      spanPtr[span].area); //this area is found above as total strip area
         	spanPtr[span].Cf[0] = spanPtr[span].Span_force[0]*tempS;
         	spanPtr[span].Cf[1] = spanPtr[span].Span_force[1]*tempS;
         	spanPtr[span].Cf[2] = spanPtr[span].Span_force[2]*tempS;
 
         	spanPtr[span].Cd = spanPtr[span].D_force*tempS;
           
		 	spanPtr[span].Cm[0]= spanPtr[span].Moment[0]*tempS/spanPtr[span].momarm[0];
		 	spanPtr[span].Cm[1]= spanPtr[span].Moment[1]*tempS/spanPtr[span].momarm[1];
		 	spanPtr[span].Cm[2]= spanPtr[span].Moment[2]*tempS/spanPtr[span].momarm[2];

		 	//strip drag coefficient
		 	spanPtr[span].Cd = spanPtr[span].D_force*tempS;
		 	//srip normal force coefficient is Cf-Cd Warning no sign information
		 	spanPtr[span].Cn = sqrt(dot(spanPtr[span].Cf,spanPtr[span].Cf)\
		 						-spanPtr[span].Cd*spanPtr[span].Cd);
		 	if(spanPtr[span].A*12+spanPtr[span].span*spanPtr[span].span*spanPtr[span].C<0) 
		 		spanPtr[span].Cn *= -1; //if integrierte circulation of strip < 0 ==> Cn<0

        	span++; //next span section
       }
    }
   	//===================================================================//
	//START quiver plot file of spanwise distribution of force coeff.
	//===================================================================//
	//for(span=0;span<info.nospanelement;span++)
	//	CreateQuiverFile(spanPtr[span].xref,spanPtr[span].Cf,1,0);
	//Cf includes ind. drag coeff., which is also added to file 
   	//===================================================================//
	//END quiver plot file of spanwise distribution of force coeff.
	//===================================================================//


   	//===================================================================//
		//END determining missing information of spanwise strips
	//===================================================================//



	fp = fopen(FCfile, "w");

    if (fp == NULL)
    {
       printf(" couldn't create flight condition file\nPress any key to exit...\n");
       exit(EXIT_FAILURE);
    }

	//===================================================================//
    //make header
    fprintf(fp,"Flight configuration defined in input file = %s\n",info.inputfilename);
    if (info.trimCL == 1) fprintf(fp,"CL trim final alpha = %.3lf\n",info.alpha*RtD);
    else fprintf(fp,"Angle of attack sweep alpha = %.3lf\n",info.alpha*RtD);
 	fprintf(fp,"Flight Condition number = %d\n",FCno);
 	fprintf(fp,"number of timesteps = %d\n",timestep);
 	fprintf(fp,"number of spanwise strips = %d\n\n",info.nospanelement);

    //header of tabulated data
    // stip number || leading edge of strip 
    fprintf(fp,"%-10s%-14s%-14s%-14s","strip","x_ref","y_ref","z_ref");
    // circulation at ref. point || force coefficients || normal force coeff (w/o ind. drag)||ind. drag coeff
    fprintf(fp,"%-14s%-14s%-14s%-14s%-14s%-14s","Gamma0","Cfx","Cfy","Cfz","Cfn (no drag)","Cd (inviscid)");
    // moment coefficints of strip
    fprintf(fp,"%-14s%-14s%-14s","Cl (roll)","Cm (pitch)","Cn (yaw)");
    //ref area || ref span || ref. length
    fprintf(fp,"%-14s%-14s%-14s","Ref. area","Ref. span","Ref. lngth");
    //reference length for moment coefficients
    fprintf(fp,"%-14s%-14s%-14s","Ref-l cl","Ref-l cm","Ref-l cn");
    // leading edge point at left side of strip || left edge chord 
    fprintf(fp,"%-14s%-14s%-14s%-14s","x1","y1","z1","chord1");
    // leading edge point at right side of strip || right edge chord 
    fprintf(fp,"%-14s%-14s%-14s%-14s","x2","y2","z2","chord2");
    // circulation on left edge || circulation on right edge
    fprintf(fp,"%-14s%-14s","Gamma1","Gamma2");
    // section force coefficeints at left, center and right (VxGamma)
    fprintf(fp,"%-14s%-14s%-14s","cn1","cn0","cn2");
	// circulation coefficients of strip (given by trailing edge)
    fprintf(fp,"%-14s%-14s%-14s","A_te","B_te","C_te");
    fprintf(fp,"\n");

    fprintf(fp,"------------------------------------------------------------"); 
	fprintf(fp,"------------------------------------------------------------");
	fprintf(fp,"------------------------------------------------------------");
	fprintf(fp,"------------------------------------------------------------");
	fprintf(fp,"------------------------------------------------------------");
	fprintf(fp,"------------------------------------------------------------");
	fprintf(fp,"------------------------------------------------------------");
 	fprintf(fp,"---------------------");
	fprintf(fp,"----------------------------------------------#\n");
 

	//===================================================================//
    //write data to file

    for(index=0;index<info.nospanelement;index++)
    {
	     fprintf(fp,"%-10d%-14lf%-14lf%-14lf",index,\
	     	spanPtr[index].xref[0],spanPtr[index].xref[1],spanPtr[index].xref[2]);
	    // circulation at ref. point || force coefficients  
	    fprintf(fp,"%-14lf%-14lf%-14lf%-14lf",spanPtr[index].Gamma0,\
	    	spanPtr[index].Cf[0],spanPtr[index].Cf[1],spanPtr[index].Cf[2]);
	    // normal force coeff (w/o ind. drag) 
	    fprintf(fp,"%-14lf%-14lf",spanPtr[index].Cn,spanPtr[index].Cd);
	    // moment coefficints of strip
	    fprintf(fp,"%-14lf%-14lf%-14lf",\
	    	spanPtr[index].Cm[0],spanPtr[index].Cm[1],spanPtr[index].Cm[2]);
	    //ref area || ref span || ref. length
	    fprintf(fp,"%-14lf%-14lf%-14lf",\
	    	spanPtr[index].area,spanPtr[index].span,spanPtr[index].chord);
	    //reference length for moment coefficients
	    fprintf(fp,"%-14lf%-14lf%-14lf",\
	    	spanPtr[index].momarm[0],spanPtr[index].momarm[1],spanPtr[index].momarm[2]);
	    // leading edge point at left side of strip || left edge chord 
	    fprintf(fp,"%-14lf%-14lf%-14lf%-14lf",\
	    	spanPtr[index].x1[0],spanPtr[index].x1[1],spanPtr[index].x1[2],spanPtr[index].chord1);
	    // leading edge point at right side of strip || right edge chord 
	    fprintf(fp,"%-14lf%-14lf%-14lf%-14lf",\
	    	spanPtr[index].x2[0],spanPtr[index].x2[1],spanPtr[index].x2[2],spanPtr[index].chord2);
	    // circulation on left edge || circulation on right edge
	    fprintf(fp,"%-14lf%-14lf",spanPtr[index].Gamma1,spanPtr[index].Gamma2);
	    // section force coefficeints at left, center and right (VxGamma)
	    fprintf(fp,"%-14lf%-14lf%-14lf",spanPtr[index].cn1,spanPtr[index].cn0,spanPtr[index].cn2);
		// circulation coefficients of strip (given by trailing edge)
	    fprintf(fp,"%-14lf%-14lf%-14lf",spanPtr[index].A,spanPtr[index].B,spanPtr[index].C);
 
		fprintf(fp,"\n");
	}
	fclose(fp);
	
 	//===================================================================//
		//END A. Saving information of spanwise strips 
 	//===================================================================//

	//===================================================================//
		//START B. Saves wing and wake of last timestep
	//===================================================================//
	
	
    //open output file for DVE information
    fp = fopen(Ftstep, "w");

    if (fp == NULL)
    {
       printf(" couldn't create flight condition file\nPress any key to exit...\n");
       exit(EXIT_FAILURE);
    }


	//===================================================================//
	//writes header
	fprintf(fp,"\n\n\nProgram Version: %s\n",PROGRAM_VERSION);
	fprintf(fp,"Flight configuration defined in input file : %s\n",info.inputfilename);
	fprintf(fp,"Flight Condition number : %d\n",FCno);
 	fprintf(fp,"This file holds surface and wake DVE informtion of last timestep (old timestep##.txt)\n");
    if (info.trimCL == 1) fprintf(fp,"CL trim final alpha : %.3lf\n",info.alpha*RtD);
    else fprintf(fp,"Angle of attack sweep alpha : %.3lf\n",info.alpha*RtD);

	fprintf(fp,"Ref. Area   : %lf\nAspect Ratio: %lf\n",info.S,info.AR);

	fprintf(fp,"surface and wake DVE after number of timesteps: %d\n",timestep);
	fprintf(fp,"elements in span direction: %d\n",info.nospanelement);
	fprintf(fp,"number of surface elements: %d\n",info.noelement);  //updated GB 2-25-20

	fprintf(fp,"no. of wings: %d  first and last span indices of each wing:",info.nowing);
	for(span=0;span<info.nowing;span++)
		fprintf(fp,"  %d  %d",info.wing1[span],info.wing2[span]);
	fprintf(fp,"\n");
	fprintf(fp,"\n");

	//writes header for information on surface elements
	fprintf(fp,"%6s %16s %16s %16s %16s %16s %16s %16s",\
	"index","xo","yo","zo","Nlift_tot","NLift_ind","NY_tot","NYi");
	fprintf(fp," %16s %16s %16s %16s %16s %16s %16s %16s %16s %16s %16s\t#\n",\
	"A","B","C","eta","xsi","nu","epsilon","psi","phiLE","phi0","phiTE");

	for(span=0;span<info.noelement;span++)
	{
		//surface element index
		fprintf(fp,"%6d",span);
		//coord. of ref point
		fprintf(fp," %16.12lf %16.12lf %16.12lf",\
				surfacePtr[span].xo[0],surfacePtr[span].xo[1],\
				surfacePtr[span].xo[2]);
		//normal forces per density
		fprintf(fp," %16.12lf %16.12lf %16.12lf %16.12lf",\
				N_force[span][0]+N_force[span][1],N_force[span][1],\
				N_force[span][2]+N_force[span][3],N_force[span][3]);
		//more info on element
		fprintf(fp," %16.12lf %16.12lf %16.12lf %16.12lf %16.12lf",\
				surfacePtr[span].A,surfacePtr[span].B,surfacePtr[span].C,\
				surfacePtr[span].eta,surfacePtr[span].xsi);
		fprintf(fp," %16.12lf %16.12lf %16.12lf",\
				surfacePtr[span].nu*RtD,surfacePtr[span].epsilon*RtD,\
				surfacePtr[span].psi*RtD);
		fprintf(fp," %16.12lf %16.12lf %16.12lf",surfacePtr[span].phiLE*RtD,\
				(surfacePtr[span].phiLE+surfacePtr[span].phiTE)*0.5*RtD,\
				surfacePtr[span].phiTE*RtD);
		fprintf(fp,"\n");
	}

	//writes header for wake information
	fprintf(fp,"\n\nwake shape after timestep: %d\n",timestep);
	fprintf(fp,"%5s%5s%16s %16s %16s %16s %16s %16s",\
				"span","time","xo","yo","zo","nu","epsilon","psi");
	fprintf(fp," %16s %16s %16s %16s %16s %16s %16s %16s %16s",\
				"U","V","W","A","B","C","eta","xsi","singfct");
	fprintf(fp," %16s %16s %16s\t#\n","phiLE","phi0","phiTE");


	for (time=0;time<=timestep;time++)
	{
		//loop across wake elements of one time/downstream location
		for(span=0;span<info.nospanelement;span++)
		{
			//trailing edge element index
			fprintf(fp,"%5d%5d",span,time);
			//coord. of ref point
			fprintf(fp," %16.12lf %16.12lf %16.12lf",\
						wakePtr[time][span].xo[0],wakePtr[time][span].xo[1],\
						wakePtr[time][span].xo[2]);
			//nu,epsilon, sweep, dihedral of elemenatry wing
			fprintf(fp," %16.12lf %16.12lf %16.12lf",\
						wakePtr[time][span].nu*RtD,\
						wakePtr[time][span].epsilon*RtD,\
						wakePtr[time][span].psi*RtD);
			//local free stream velocity at center
			fprintf(fp," %16.12lf %16.12lf %16.12lf",\
						wakePtr[time][span].u[0],wakePtr[time][span].u[1],\
						wakePtr[time][span].u[2]);
			fprintf(fp," %16.12lf %16.12lf %16.12lf",\
						wakePtr[time][span].A,wakePtr[time][span].B,\
						wakePtr[time][span].C);
			//element half span and half chord
			fprintf(fp," %16.12lf %16.12lf",\
						wakePtr[time][span].eta,wakePtr[time][span].xsi);
			//element half span and half chord
            fprintf(fp," %16.12lf",wakePtr[time][span].singfct);//    999.999);
			//leading-edge, mid-chord, and trailing edge sweep2
			fprintf(fp," %16.12lf %16.12lf %16.12lf",\
				wakePtr[time][span].phiLE*RtD,\
				wakePtr[time][span].phi0*RtD,wakePtr[time][span].phiTE*RtD);

			fprintf(fp,"\n");
		}
		fprintf(fp,"\n");
	}
	fclose(fp);

	//===================================================================//
		//END B. Saves wing and wake of last timestep
	//===================================================================//

	//===================================================================//
		//START C. Saves Surface DVE of last timestep
	//===================================================================//
	
    //open output file for DVE information
    fp = fopen(SDVE,"w");

    if (fp == NULL)
    {
       printf(" couldn't create flight condition file\nPress any key to exit...\n");
       exit(EXIT_FAILURE);
    }

	//===================================================================//
	//writes header
	fprintf(fp,"Program Version = %s\n",PROGRAM_VERSION);
	fprintf(fp,"Flight configuration defined in input file = %s\n",info.inputfilename);
	fprintf(fp,"Flight Condition number = %d\n",FCno);
 	fprintf(fp,"This file holds surface DVE informtion of last timestep\n");
    if (info.trimCL == 1) fprintf(fp,"CL trim final alpha = %.3lf\n",info.alpha*RtD);
    else fprintf(fp,"Angle of attack sweep alpha = %.3lf\n",info.alpha*RtD);

	fprintf(fp,"surface and wake DVE after number of timesteps = %d\n",timestep);
	fprintf(fp,"elements in span direction = %d\n",info.nospanelement);
	fprintf(fp,"number of surface elements = %d\n",info.noelement);

	fprintf(fp,"Coefficients are based on DVE area and DVE freestreaem velocity\n");
		// 		N_force;	surface DVE's normal forces/density
					//[0]: free stream lift, [1]: induced lift,
					//[2]: free stream side, [3]: induced side force/density
					//[4]: free str. normal, [5]: ind. normal frc/density
               		//[6,7,8]: eN_x, eN_y, eN_z in global ref. frame


	//writes header for information on surface elements
	fprintf(fp,"%6s %16s %16s %16s %16s ",\
	"index","xo","yo","zo","cfn");
	fprintf(fp,"%16s %16s %16s ",\
	"normal-x","normal-y","normal-z");
	fprintf(fp,"%16s %16s %16s %16s %16s %16s ",\
	"u","v","w","area","eta","xsi");
	fprintf(fp,"%16s %16s %16s %16s %16s %16s ",\
	"nu","epsilon","psi","phiLE","Phi0","phiTE");
	fprintf(fp,"%16s %16s %16s %16s %16s %16s ",\
	"x1","y1","z1","x2","y2","z2");
	fprintf(fp,"%16s %16s %16s \t#\n","A","B","C");
	
	for(index=0;index<info.noelement;index++)
	{
		//surface element index
		fprintf(fp,"%6d ",index);
		//coord. of ref point
		fprintf(fp,"%16lf %16lf %16lf ",\
				surfacePtr[index].xo[0],surfacePtr[index].xo[1],\
				surfacePtr[index].xo[2]);
		//normal forces coefficient
		fprintf(fp,"%16lf ",(N_force[index][4]+N_force[index][5])*0.5/\
				(dot(surfacePtr[index].u,surfacePtr[index].u)*surfacePtr[index].S));
		//DVE normal vector
		fprintf(fp,"%16lf %16lf %16lf ",\
				surfacePtr[index].normal[0],surfacePtr[index].normal[1],\
				surfacePtr[index].normal[2]);
		//free stream velocity in control point
		fprintf(fp,"%16lf %16lf %16lf ",\
				surfacePtr[index].u[0],surfacePtr[index].u[1],surfacePtr[index].u[2]);
		// area, eta, xsi
		fprintf(fp,"%16lf %16lf %16lf ",\
				surfacePtr[index].S,surfacePtr[index].eta,surfacePtr[index].xsi);
		//nu, epsilon, psi
		fprintf(fp,"%16lf %16lf %16lf ",\
				surfacePtr[index].nu*RtD,surfacePtr[index].epsilon*RtD,\
				surfacePtr[index].psi*RtD);
		//phi leading and trailing edge
		fprintf(fp,"%16lf %16lf %16lf ",surfacePtr[index].phiLE*RtD,\
				(surfacePtr[index].phiLE+surfacePtr[index].phiTE)*0.5*RtD,\
				surfacePtr[index].phiTE*RtD);
		//x1 left leading edge point
		fprintf(fp,"%16lf %16lf %16lf ",\
				surfacePtr[index].x1[0],surfacePtr[index].x1[1],surfacePtr[index].x1[2]);
		//x2 right leading edge point
		fprintf(fp,"%16lf %16lf %16lf ",\
				surfacePtr[index].x2[0],surfacePtr[index].x2[1],surfacePtr[index].x2[2]);
		//circulation coefficients
		fprintf(fp,"%16lf %16lf %16lf ",\
				surfacePtr[index].A,surfacePtr[index].B,surfacePtr[index].C);
		fprintf(fp,"\n");
	}

	fclose(fp);

	//===================================================================//
		//END C. Saves Surface DVE of last timestep
	//===================================================================//
	FREE1D(&FCfile,length+11);
	FREE1D(&Ftstep,length+11);
	FREE1D(&SDVE,length+11);

}
//===================================================================//
		//END SaveSpanInfo
//===================================================================//
