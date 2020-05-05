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
//void Save_Input_File(const GENERAL);
//saves results of current timestep to file
void Save_Timestep(const GENERAL,const int,DVE **,const DVE *,double **);
//saves forces and moments of surface DVEs
void Save_SurfaceDVE_Loads(const GENERAL,const int,const DVE *);

void CreateQuiverFile(const double[3], const double[3],const int,const int);

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
	char filename[137];	//file path and name

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
	char filename[137];	//file path and name

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
	char filename[126];	//file path and name

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
	char filename[126];	//file path and name

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
	char filename[126];	//file path and name

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
	char filename[126];	//file path and name

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
	char filename[132];	//file path and name

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

    if (source == NULL)
    {
       printf(" couldn't copy input file\nPress any key to exit...\n");
       exit(EXIT_FAILURE);
    }

    sprintf(targetfile, "%s%s",targetdir,sourcefile);

    target = fopen(targetfile, "w");

    if (target == NULL)
    {
       fclose(source);
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
		//START of Save_Timestep
//===================================================================//
void Save_Timestep(const GENERAL info,const int timestep,DVE **wakePtr,\
								const DVE *surfacPtr,double **N_force)
{
//saves results of current timestep to file

	int time,span;		//loop counter
	FILE *fp;			//output file
	char filename[133];	//file path and name

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
	char filename[133];	//file path and name

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
char filename[133];	//file path and name
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
	//		filenum		- Number to add to filename (ex. quiver4.txt) Zero will not add any value
	// 
	// Example of use to plot eN
	//	for(i=0;i<info.noelements;i++){
	//	if(i==0){CreateQuiverFile(surfacePtr[l].xo, eN,0);}
	//	else{CreateQuiverFile(surfacePtr[l].xo, eN,1);}
	//	}
	//
	//
	// D.F.B. in Braunschweig, Germany, Mar. 2020


	FILE *fp;		//output file
	char filename[126];	//file path and name 

	if (filenum == 0) {
		// Either open file or add to file
		sprintf(filename, "%s%s", info.output,"quiver.txt");
		if (idx == 0) { fp = fopen(filename, "w"); }
		else { fp = fopen(filename, "a"); }
	}
	else {
		sprintf(filename, "%s%s%d%s",info.output,"quiver", filenum, ".txt");
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
