/*******************************************************************************
***********						tAo main source file					********
********************************************************************************
	For compilation and installation check the file ../tao/README.
	Main author since 1995 Daniel Garcia-Castellanos, danielgc@ictja.csic.es. 
	Copyright details and other information in ../tao/doc/ 
********************************************************************************

Memory debugging with: 
valgrind --dsymutil=yes --track-origins=yes --tool=memcheck --leak-check=full `which tisc` linear_range -tf0

	COMMENTS (programmer's agenda)
	-Track rock particles in thrusts and sediments. 
	-Implement grain size in transport. This to calculate grain size distribution in sediment blocks, and as a first step for the next point once transitory flow is implemented. 
	-Implement sediment load effect on erosion (Sklar). Intersting for acceleration of erosion during lake overtopping.
	-Implement transitory water flow.
	-DONE. Implement sediment compaction (easy in calculate_topo and when writting the pfl file).
*/


#include "tao.h"
#include "taolib.c"
#include "taoio.c"



int main(int argc, char **argv)
{
	/*get input parameters and files*/
	inputs (argc, argv) ;

	fprintf(stdout, "\nT= %.4f My", Time/Matosec);

	if (switch_dt_output) {calculate_topo(topo); Write_Ouput();}

	/*MAIN LOOP: In this loop time increases from Timeini to Timefinal*/
	do { 
		/*Remember topography before tectonics and flexure*/
		calculate_topo(topo);

		/*Calculate tectonic deformation and tectonic load*/
		tectload();

		/*Sea level variations*/
		calculate_sea_level();

		/*Define & solve elastic equation*/
		Elastoplastic_Deflection();

		/*Define & solve viscoelastic flexure equation*/
		Viscous_Relaxation();

		/*Calculates surface water flow and sediment-erosion*/
		surface_processes();

		Time += dt;
		fprintf(stdout, "\nT= %.4f My", Time/Matosec);

		if (switch_dt_output) Write_Ouput();
	} while (Time < Timefinal-dt/10);

	The_End();

	return(1);
}



/**************************************************/
/****	  SUBROUTINES  IN  RANDOM  ORDER	 *****/
/**************************************************/



int inputs(int argc, char **argv)
{
	int 	reformat=0 ;
	char 	resume_filename[MAXLENFILE], 
		command[MAXLENLINE], 
		load_file_name[MAXLENLINE];
	BOOL	success_def_prm=NO, switch_initial_geom=NO;

	run_type = 0;
	nmax_input_points = 5000;
	switch_strs_history = YES;
	setbuf(stdout, NULL);

	putenv("tao_dir=" TAODIR); 
	
	/*Version of tAo will be matched against the parameters file *.PRM*/
	/*¡¡ UPDATE template.PRM !!*/
	strcpy(version, "tAo_2018-08-07");

	/*Default parameter values are read from ./tao/doc/template.PRM*/
	sprintf(projectname, "%s/doc/template", TAODIR);
	success_def_prm = read_file_parameters(0, 0);
	sprintf(projectname, "");

	for (int iarg=1; iarg<argc; iarg++) {
		if (argv[iarg][0] == '-') {
			float 	value;
			char 	prm[MAXLENLINE];
			for (int ilet=2; ilet < strlen(argv[iarg])+2; ilet++) 
				prm[ilet-2] = argv[iarg][ilet];
			value=atof(prm);
			switch (argv[iarg][1]) {
				case 'f':
					reformat=1;
					if (argv[iarg][2]) reformat = value;
					break;
				case 'F':
					run_type=2;
					if (strlen(prm)>0) strcpy (resume_filename, prm);
					else		   sprintf(resume_filename, "%s"".all", projectname);
					break;
				case 'h':
					switch (argv[iarg][2]) {
						case 'p':
							fprintf(stderr, "\nFile ./tao/doc/template.PRM (sample parameters file) follows in stdout:\n") ;
						sprintf(command, "cat %s/doc/template.PRM", TAODIR);
						system(command) ;
						break;
						case 'c':
							fprintf(stderr, "\nFile ./tao/doc/template.PRM (sample parameters file) follows in stdout:\n") ;
						sprintf(command, "cat %s/doc/template.PRM | %s/script/cleanPRM", TAODIR, TAODIR);
						system(command) ;
						break;
						case 'u':
							fprintf(stderr, "\nFile ./tao/doc/template.UNIT (sample unit file) follows in stdout:\n") ;
						sprintf(command, "cat %s/doc/template.UNIT", TAODIR);
						system(command) ;
						break;
						default:
						fprintf(stderr, "\nFile ./tao/doc/tao.info.txt follows:\n") ;
						sprintf(command, "more %s/doc/tao.info.txt", TAODIR);
						AUTHORSHIP;
						system(command) ;
						break;
					}
					fprintf(stderr, "\n") ;
					exit(0);
				case 'Q':
					run_type=1;
					strcpy(load_file_name, prm);
					break;
				case 'V':
					verbose_level = 1;
					if (argv[iarg][2]) verbose_level = value;
					break;
			}
		}
		else {
			if (run_type != 2) run_type=10;
			if (strlen(projectname)<1) strcpy(projectname, argv[iarg]);
		}
	}

	if (verbose_level>=1) {
	  fprintf(stdout, 	
	  	"\n**** tAo: 2D LITHOSPHERIC FLEXURE, EROSION/SEDIMENTATION, AND FORELAND BASIN MODELING ****"
	  	"\nVersion: %s", version);
	  fflush(stdout);
	}

	if (!run_type) {
		syntax();
		fprintf(stdout, "\n\nType %s -h for further information.\n", argv[0]);
		exit(0);
	}

	nloads=0; n_image=0; nlakes=0;
	numBlocks=0; i_first_Block_load=0; i_Block_insert=0;
	nwrotenfiles=0; switch_topoest=NO;

	switch (run_type)
	{
		case 0:
			fprintf(stdout, "\n\n*** END of run *** \n\n");
			exit(0);
		case 1:
			interpr_command_line_opts(argc, argv);
			Direct_mode(load_file_name);
			exit(0);
		case 2:
			read_file_resume(resume_filename);
			interpr_command_line_opts(argc, argv);
			if (verbose_level>=1) fprintf(stdout, "\nResuming project '%s'. Timefinal=%.1f My", projectname, Timefinal/Matosec);
			if (switch_dt_output) n_image--; /*Don't produce 2 jpg's of the same stage of restart*/
			return(1);
		case 10:
			if (!read_file_parameters(verbose_level>=1, 0)) {
				syntax();
				fprintf(stderr, "\n\tAvailable parameter files in this directory:\n");	
				system("ls *.PRM");
				if (!success_def_prm) {
					PRINT_ERROR("\t\aDefault parameters file './tao/doc/template.PRM' could not be read.\n"); 
				}
				exit(0);
			}
			if (reformat) {
				sprintf(projectname, "%s/doc/template", TAODIR);
				read_file_parameters(0, reformat);
				exit(0);
			}
			interpr_command_line_opts(argc, argv);
			Allocate_Memory();
			break;
	}

	{
		char filename[MAXLENLINE]; FILE *file;
		sprintf(filename, "%s.out", projectname);
		if (switch_file_out) {
			if ((file = fopen(filename, "w")) == NULL) {
				PRINT_ERROR("Cannot open standard output file %s.\n", filename);
			}
			else {
				if (verbose_level>=1) fprintf(stdout, "\nInfo: standard output redirected to %s.\n", filename);
			stdout=file;
			}
		}
		else {
		remove(filename);
		}
	}

	dx = (xf-x0) / (Nx-1);
	dz = (100000) / (Nz-1);
	dt *= Matosec;
	tau *= Matosec; 
	dt_record *= Matosec;
	Timefinal *= Matosec;
	Timeini *= Matosec;
	Time = Timeini;
	last_time_file_time = Timeini - 9999*dt_record;	/*very old*/
	dt_eros *= Matosec;
	Kerosdif /= secsperyr;
	Keroseol /= Matosec;
	Ksedim   /= Matosec;
	rain /= secsperyr*1e3;
	if (hydro_model==1) Krain *= 1e6/Matosec/1e3/1e3;
	evaporation_ct /= secsperyr*1e3;
	lost_rate *= 1e-2 * 1e-3;
	temp_sea_level += TEMP_FREEZE_WATER;
	switch_write_file_Blocks = 1;


	if (verbose_level>=3) {
		fprintf(stdout, "\nPlate model: %d \tBC=%d\tNx=%d from %.2f to %.2f (dx=%f m)", isost_model, boundary_conds, Nx, x0, xf, dx);
		fprintf(stdout, "\nDensities: asthenosphere =%f; mantle=%f; crust=%.2f", densasthen, densmantle, denscrust);
		fprintf(stdout, "\nTiming: from %.2f to %.2f every %.2f My\n", Timeini/Matosec, Timefinal/Matosec, dt/Matosec);
	}


	/*Test of incompatibilities between parameters*/
	if (densenv && water_load)		{ water_load=NO; 	PRINT_WARNING("Sea not possible when densenv<>0. Sea switch turned off.") ; }
	if (!erosed_model && (Ksedim || Kerosdif || Keroseol))	{ Ksedim=Kerosdif=Keroseol=0 ; if (verbose_level>=1) PRINT_WARNING("Erosion-sedimentation unswitched. Parameters Ksedim, Kerosdif & Keroseol have no effect.") ; }
	if (!water_load && Ksedim)		{ Ksedim=0;			PRINT_WARNING("Warning: Sea presence isn't switched, Ksedim has no effect.");}
	if (!erosed_model)  			{ Ksedim=Kerosdif=Keroseol=0 ; }
	if (!hydro_model && erosed_model>1)  	{ erosed_model=1 ; PRINT_WARNING("Warning: hydro_model is 0; erosed_model is switched to 1 accordingly.") ; }
	if (switch_ps && !switch_write_file)	{ if (verbose_level>=3)	PRINT_WARNING("Warning: switch_write_file needed to make a postscript. Postscript may not be done.") ; }
	if (isost_model!=2)  			{ tau=0; }
	if (tau<=0 && isost_model==2)  		{ isost_model=1; }
	if (boundary_conds == 0)		{ appmoment = 0; }
	if (boundary_conds == 2)		{ vert_force=0; appmoment = 0; }


	sprintf(command, "rm -f %s.temp0 %s.mtrz %s.grv_mod", projectname, projectname, projectname);
	system(command);

	read_file_sea_level(); calculate_sea_level();
	read_file_horiz_record_time();
	read_file_Te();
	read_file_Crust_Thick(crust_thick_default);
	read_file_Upper_Crust_Thick(upper_crust_thick_default);	
	read_file_YSE();
	read_file_Temperature() ;
	Init_Stress();

	switch_initial_geom = read_file_initial_deflection(w) + read_file_initial_topo(topo) ;
	for (int i=0; i<Nx; i++)  { 
		topo[i] += random_topo * ((((float) rand()) / ((float) RAND_MAX)) -.5);
		topo[i] += zini;
		Blocks_base[i] = topo[i];
		topo[i] -= w[i];
		h_water[i] = MAX_2(sea_level-topo[i], 0);
	}
	if (switch_initial_geom && verbose_level>=1) {
		float	altmax=-1e9, altmin=1e9;
		for (int i=0; i<Nx; i++) {
			if (altmax<topo[i])	altmax=topo[i];
			if (altmin>topo[i])	altmin=topo[i];
		}
		fprintf(stdout, "\n  alt.Init.:  max = %9.1f m	 min = %9.1f m   ", altmax, altmin);
	}

	return(1);
}



int interpr_command_line_opts(int argc, char **argv) 
{
	/*Interpretates the command line options described in tao.info.txt*/

	PRINT_INFO("Enetering command line interpretation.");
	for (int iarg=1; iarg<argc; iarg++) {
		if (argv[iarg][0] == '-') {
			float 	value, value2;
			char 	prm[MAXLENLINE], prm2[MAXLENLINE], *ptr;
			for (int ilet=2; ilet < strlen(argv[iarg])+2; ilet++) 
				prm[ilet-2] = argv[iarg][ilet];
			for (int ilet=3; ilet < strlen(argv[iarg])+2; ilet++) 
				prm2[ilet-3] = argv[iarg][ilet];
			value  = atof(prm);
			value2 = atof(prm2);
			PRINT_DEBUG("\aArgument: %s", argv[iarg]);
			switch (argv[iarg][1]) {
				case 'A':
					grav_anom_type = value;
					break;
				case 'B':
					boundary_conds = value;
					break;
				case 'D':
					if (run_type!=2) {
						x0 = atof(strtok(prm, "/"));
						xf = atof(strtok(NULL, "/"));
					}
					if (xmin<x0) xmin=x0;
					if (xmax>xf) xmax=xf;
					else fprintf(stdout, "\nWarning: Impossible to change the domain when resuming a model.");
					break;
				case 'd':
					if (run_type!=2) Nx = (int) (xf-x0)/value + 1;
					else fprintf(stdout, "\nWarning: Imposible to change dx or Nx value when restarting a model.");
					break;
				case 'L': /*OLD*/
					if (run_type!=2) x0 = value;
					else fprintf(stdout, "\nWarning: Imposible to change x0 value when restarting a model.");
					if (xmin<x0) xmin=x0;
					break;
				case 'M':
					isost_model = value;
					if (prm[1] == 't') switch_strs_history = NO;
					break;
				case 'm':
					appmoment = value;
					break;
				case 'N':
					Nx = value;
					break;
				case 'o':
					switch_file_out=YES;
					break;
				case 'P':
					switch_ps=YES;
					switch_write_file_Blocks=YES;
					if (argv[iarg][2] == 'c') {
						switch_dt_output=YES;
						strcpy(gif_geom, "");
						if (strlen(prm2)>0) strcpy(gif_geom, prm2);
					}
					break;
				case 'p':
					horz_force = value;
					break;
				case 'q':
					ptr = strtok(prm, "=");
					match_parameter(ptr, strtok(NULL, "/"), 1, 0, prm);
					break;
				case 'R': /*OLD*/
					if (run_type!=2) xf = value;
					else fprintf(stdout, "\nWarning: Imposible to change xf value when restarting a model.");
					if (xmax>xf) xmax=xf;
					break;
				case 'r':
					switch (argv[iarg][2]) {
						case 'e':	densenv	= value2; 	break;
						case 'c':	denscrust  = value2; 	break;
						case 'i':	densinfill = value2; 	break;
						case 'm':	densmantle = value2; 	break;
						case 'a':	densasthen = value2; 	break;
					}
					break;
				case 'S':
					{
						int iblock, nblocks;
						struct BLOCK_1D	Block_aux;
						iblock = atoi(strtok(prm, "/"));
						nblocks = atoi(strtok(NULL, "/"));
						PRINT_INFO("Block %d will be moved by %d positions", iblock, nblocks);
						Block_aux=Blocks[iblock];
						if (nblocks>0) {
							for (int iu=iblock; iu<iblock+nblocks; iu++) {
								Blocks[iu]=Blocks[iu+1];
								PRINT_INFO("%d = %d", iu, iu+1);
							}
						}
						else {
							for (int iu=iblock; iu>iblock+nblocks; iu--) {
								Blocks[iu]=Blocks[iu-1];
								PRINT_INFO("%d = %d", iu, iu-1);
							}
						}
						Blocks[iblock+nblocks]=Block_aux;
						PRINT_INFO("%d = %d", iblock+nblocks, numBlocks-1);
					}
					break;
				case 's':
					vert_force = value;
					break;
				case 'T':
					Te_default = value;
					break;
				case 't':
					switch (argv[iarg][2]) {
						case 'i':	Timeini	= value2; 	if (run_type==2) Timeini *= Matosec;	break;
						case 'f':	Timefinal  = value2;if (run_type==2) Timefinal *= Matosec;	break;
						case 'd':	dt		 = value2; 	if (run_type==2) dt *= Matosec;			break;
						case 'e':	dt_eros	= value2; 	if (run_type==2) dt_eros *= Matosec;	break;
						case 'v':	tau		= value2; 	if (run_type==2) tau *= Matosec; 		break;
						case 'r':	dt_record = value2; if (run_type==2) dt_record *= Matosec;	break;
					}
					break;
				case 'V':
					verbose_level = 1;
					if (argv[iarg][2]) verbose_level = value;
					break;
				case 'v':
					{
						float density, velocity;
						density = atof(strtok(prm, "/"));
						velocity = atof(strtok(NULL, "/"));
						for (int iu=0; iu<numBlocks; iu++) {
							if (Blocks[iu].density==-density || iu==density) {
								Blocks[iu].vel=velocity*1e3/Matosec;
								Blocks[iu].last_vel_time=Time-dt;/*!!*/
								Blocks[iu].last_shift=0;
							}
						}
					}
					break;
			}
		}
	}
	return(1);
}



int Direct_mode(char *load_file_name)
{
	int 	i;
	FILE 	*file;

	/*Solves flexure problem in direct mode: taking a single external load file and 
	writing deflection in standar output. There are no other input files neither 
	output files*/

	PRINT_INFO("Entering direct mode. x0,xf,Nx= %.1f,%.1f,%d\n", x0,xf, Nx);
	dx = (xf-x0) / (Nx-1) ;
	Allocate_Memory();
	if (strcmp(load_file_name, "")) {
		if ((file = fopen(load_file_name, "rt")) == NULL) {
			fprintf(stderr, "\nLoad file '%s' not found.\n", load_file_name);
			exit(0);
		}
		readinterplin(file, h_last_unit, Nx, x0, xf) ;
		fclose(file);
	}
	for (i=0; i<Nx; i++) {
		D[i] = ET2RIG(Te_default); 
		Dq[i] = h_last_unit[i];
	}
	Elastoplastic_Deflection();
	fprintf(stdout, "\n\n#x[km]\t\tw[m]\t\tpressure[Pa]\n"); 
	for (i=0; i<Nx; i++) 
		fprintf(stdout, "%8.1f\t%8.1f\t%8.1f\n", 
			(x0+i*dx)/1e3, w[i], h_last_unit[i]);
	fprintf(stdout, "\n"); 
}




int tectload()
{
	/*
	CALCULATES NEW LOAD INCREMENT FROM UNIT FILES, Returns 1 if elastic
	flexure must be done (i.e, if changes in load  occurred), 0 otherwise.
	*/

	PRINT_ARRAY_INFO(topo, "topogr.", "m", "m2") 

	/*Reads external load from file(s)*/
	while (read_file_unit());

	/*Moves Blocks*/
	move_Blocks();

	/*Interpolates loads through time*/
	gradual_Block();

	Repare_Blocks();

	return (1);
}




int Elastoplastic_Deflection()
{
	int 	i;
	BOOL	load_changes=NO;

	/*
	CALLS SUBROUTINES TO SOLVE ELASTIC OR ELASTIC-PLASTIC FLEXURE WITH 
	THE NEW LOAD INCREMENT AND ADDS NEW DEFLECTION INCREMENT TO THE 
	TOTAL DEFLECTION.
	*/

	for (i=0; i<Nx; i++) if (Dq[i]) load_changes = YES;
	if (isost_model>0 && (load_changes || (Time==Timeini && (horz_force || vert_force || appmoment)))) {
	  fprintf(stdout, " e");	  fflush(stdout);
	  if (isost_model<3) {
			if (!Te_default) {
			/*LOCAL ISOSTASY*/
			float Krest;
			for (i=0; i<Nx; i++) {
				GET_KREST(Krest, q, i)
				Dw[i] = Dq[i] / Krest;
			}
			}
			else {
				/*REGIONAL ISOSTASY*/
		/*Pure Elastic Flexure without rheological inputs*/
		float 	*moment;
		double	**A,		/*Linear System Matrix (diagonal terms)*/
			*b;		/*Independent Column*/
		int 	NDs=3, NDi=3;
		A = alloc_matrix_dbl (Nx, NDi+1+NDs);
		b = (double *) calloc (Nx , sizeof(double));
		moment = (float *) calloc (Nx , sizeof(float));
		LES_matrix(A, b, D, q, Dq, w, NO) ;
		solveLES(A, b, Nx, NDs, NDi, Dw) ;
		for (i=0; i<Nx; i++) free(A[i]);
		free(A); free(b); 

		for (i=0;i<Nx;i++) {
			w[i] += Dw[i] ;
			if (i != 0 && i != Nx-1) 
				moment[i] += 
					-D[i] * (Dw[i-1] - 2*Dw[i] + Dw[i+1]) 
					/ pow(dx,2) ;
		}
		if (isost_model!=2) flexural_stats(moment);
		free(moment);

		}
	  }
	  else {
		/*Elasto-plastic flexure with EET rheological calculation*/
		Rheo_Flex_Iter();
	  }

	  if (switch_topoest) {
	  	  /*Defines the thickness of last infill Block*/
	  	  for (i=0; i<Nx; i++)
	  		  Blocks[i_first_Block_load-1].thick[i] +=  MAX_2(Dw[i], 0) ;
	  }
	}

	PRINT_ARRAY_INFO(Dq, "load_incrm", "N/m2", "N/m") 

	/*Adds all this load to the total load array*/
	for (int i=0; i<Nx; i++)  {
		q[i] += Dq[i];
		if (Dq[i]) load_changes = YES;
	}
	if (load_changes) fprintf(stdout, " d");
	/*Resets deflection and load arrays*/
	for (i=0;i<Nx;i++)  Dq[i]=Dw[i]=0;

	fflush(stdout);
	return(1);
}



int surface_processes()
{
	/* 
	CALCULATES EROSION AND SEDIMENTATION:
	*/
	float 	eros_level;
	BOOL	switch_horiz_record=NO;

	PRINT_DEBUG("erosed_model=%d hydro_model=%d", erosed_model, hydro_model);

	total_sed_mass=total_bedrock_eros_mass=0;

	if (!erosed_model && !hydro_model) return (0);
	switch_topoest=NO;
	

	/*Creates a new sediment Block if necessary*/
	if (erosed_model) {
		int i;
		float TimelastBlock=-9999*Matosec;
		for (int i=0; i<Nx; i++) eros_now[i]=0;
		for (int i=0; i<numBlocks; i++) 
		if (Blocks[i].age > TimelastBlock && Blocks[i].density==denssedim) TimelastBlock = Blocks[i].age;
		for (int i=0; i<n_record_times; i++) 
		if (Time>horiz_record_time[i]-dt/2 && Time<=horiz_record_time[i]+dt/2) 
			switch_horiz_record=YES;
		if (Time == Timeini 
		  || ((Time-TimelastBlock)>(dt_record-.001*dt) && dt_record && !n_record_times)
		  || switch_horiz_record) {
			insert_new_Block(numBlocks);
			Blocks[numBlocks-1].type = 'S' ;
			Blocks[numBlocks-1].density = denssedim ;
			Blocks[numBlocks-1].erodibility = erodibility_sed ;
			Blocks[numBlocks-1].detr_ratio = calloc(Nx, sizeof(float));
			Blocks[numBlocks-1].detr_grsize = calloc(Nx, sizeof(float));
		}
	}
	/*Calculates height at every point*/
	calculate_topo(topo);

	/*Diffusive Erosion*/
	Diffusive_Eros_1D (topo, Kerosdif, dt, dt_eros/5);

	/*Adds background erosion and sea sedimentation*/
	constant_rate_eros (topo, Keroseol, Ksedim, sea_level, water_load, dt, dt_eros, Time, n_eros_level_input_points, var_eros_level, &eros_level);

	/*Adds fluvial transport*/
	Surface_Transport (topo, dt, dt_eros, erosed_model);


	/*Calculates water column load*/
	calculate_water_load();



	if (n_eros_level_input_points)
		fprintf(stdout,
			"\teros_lvl.: %8.1f m", eros_level);

	if (verbose_level>=1) {
		if (erosed_model>1) {
			int i, i_biggest_nosea=0, n_biggest_nosea=0;
			double error_water;
			for (i=1; i<=nlakes; i++) {
				if (Lake[i].n > n_biggest_nosea) {
					if (Lake[i].n_sd) {
					if (topo[Lake[i].sd[0]] > sea_level)
						n_biggest_nosea = Lake[i].n; i_biggest_nosea = i;
					}
					else	n_biggest_nosea = Lake[i].n; i_biggest_nosea = i;
				}
			}
			for (i=1; i<=nlakes; i++) {
				float vol=0;
				if (i==i_biggest_nosea || (Lake[i].n>ceil((double) Nx/100) && verbose_level>=3) || verbose_level>=4) {
				PRINT_SUMLINE("lake %d/%d: %6.2f km2 %6.1f km%5.0f m ", 
					i, nlakes, Lake[i].vol/1e6, Lake[i].n*dx/1e3, Lake[i].alt);
				/*write the position of the lowest node*/
				if (Lake[i].n) fprintf(stdout, "bott@ %.0f km", (Lake[i].cell[0]*dx+x0)/1e3);
				if (Lake[i].n_sd) {
					if (topo[Lake[i].sd[0]]>sea_level) {
					 fprintf(stdout, " %d out@ %.0f km %8.2e m3/s", Lake[i].n_sd, (Lake[i].sd[0]*dx+x0)/1e3, drainage[Lake[i].sd[0]].discharge/*Lake_Input_Discharge(i)*/);
					}
					else fprintf(stdout, " Sea	   %8.2e m3/s", Lake_Input_Discharge(i));
				}
				else	 fprintf(stdout, " Endorheic %8.2e m3/s", Lake_Input_Discharge(i));
				}
			}
		}
	}
	if (verbose_level>=1) {
		int i, imaxerosion=SIGNAL, imaxsediment=SIGNAL, imaxaltsediment=SIGNAL;
		float top_Block, totalerosion=0, volsediment=0, maxerosion=0, maxsediment=0, max_altit_seds=-9999;
		for (i=(xmin-x0)/dx; i<Nx-(xf-xmax)/dx; i++) {
			totalerosion  -= MASS2SEDTHICK_1D(total_erosion[i])*dx;
			if (maxerosion<total_erosion[i])
				{maxerosion = total_erosion[i]; imaxerosion=i;}
			if (maxsediment<-total_erosion[i])
				{maxsediment = -total_erosion[i]; imaxsediment=i;}
			top_Block = Blocks_base[i]-w[i];
			for (int j=0; j<numBlocks; j++) {
				top_Block+=Blocks[j].thick[i];
				if (Blocks[j].density==denssedim) volsediment+=Blocks[j].thick[i]*dx;
				if (Blocks[j].density==denssedim && top_Block>max_altit_seds && Blocks[j].thick[i]>1) {
					max_altit_seds = top_Block; imaxaltsediment=i;
				}
			}
		}
		PRINT_SUMLINE("total_cumul_eros-sed= %8.1f km2\tvol.sediment= %8.1f km2",
			totalerosion/1e6, volsediment/1e6
		);
		if (imaxerosion!=SIGNAL) fprintf(stdout,
			"\n  erosion_max.   = %8.1f m  \t@x= %5.1f km", 
			MASS2SEDTHICK_1D(maxerosion),  (x0+imaxerosion*dx)/1000
		);
		if (imaxsediment!=SIGNAL) fprintf(stdout,
			"\n  sediment.max.  = %8.1f m  \t@x= %5.1f km" 
			"\n  max.altit.sed. = %8.1f m  \t@x= %5.1f km", 
			MASS2SEDTHICK_1D(maxsediment), (x0+imaxsediment*dx)/1000, 
			max_altit_seds, (x0+imaxaltsediment*dx)/1000
		);
	}

	fflush(stdout);
	return(1);
}



int flexural_stats (float *moment) {
	if (verbose_level>=1) {
		/*prints flexural statistics*/
		int 	i, iwmindt=SIGNAL, iwmaxdt=SIGNAL, idwmindt=SIGNAL, idwmaxdt=SIGNAL, ihmaxdt=SIGNAL;
		float	shear=0, shearmax=-1e19, xshearmax=-1e19, shearmin=+1e19, xshearmin=+1e19, 
			momentmax=-1e19, xmomentmax=-1e19, momentmin=+1e19, xmomentmin=+1e19, 
			xfirstnodo=0, x, 
			wmaxdt=-1e19, wmindt=+1e19, dwmaxdt=-1e19, dwmindt=+1e19;
		float 	Warea=0, Dmean=0;

		for (i=1; i<Nx-1; i++) {
			if (w[i]*w[i+1] <= 0  &&  xfirstnodo == 0 )
				xfirstnodo = i*dx+x0 ;
			if ( i < Nx-2 )
				shear = - D[i] * (w[i+2] - 3*w[i+1] + 3*w[i] -w[i-1]) / pow(dx,3);
			if (shearmin > shear)
				{ shearmin = shear;	xshearmin=x0+dx*(i+.5); }
			if (shearmax < shear)
				{ shearmax = shear;	xshearmax=x0+dx*(i+.5); }
			if (shearmax < shear)
				{ shearmax = shear;	xshearmax=x0+dx*(i+.5); }
			if (momentmax < moment[i] && x0+dx*i>xmin && x0+dx*i<xmax)
				{ momentmax = moment[i];	xmomentmax=x0+dx*i; }
			if (momentmin > moment[i] && x0+dx*i>xmin && x0+dx*i<xmax)
				{ momentmin = moment[i];	xmomentmin=x0+dx*i; }
		}
		for (i=1; i<Nx-1; i++) {
			if (x0+i*dx>=xmin && x0+i*dx<=xmax) {
			if (wmindt>w[i]) {		wmindt=w[i]; 	iwmindt=i; }
			if (wmaxdt<w[i]) {		wmaxdt=w[i]; 	iwmaxdt=i; }
			if (dwmindt>Dw[i]) {		dwmindt=Dw[i]; 	idwmindt=i; }
			if (dwmaxdt<Dw[i]) {		dwmaxdt=Dw[i]; 	idwmaxdt=i; }
			Warea += w[i]*dx;
			Dmean += D[i]/Nx;
			}
		}
		PRINT_SUMLINE("moment_max.	= %10.3e N  \t@x= %5.1f km", momentmax, xmomentmax/1000);
		PRINT_SUMLINE("moment_min.	= %10.3e N  \t@x= %5.1f km", momentmin, xmomentmin/1000);
		PRINT_SUMLINE("shear_max.	 = %10.3e N/m\t@x= %5.1f km", shearmax, xshearmax/1000);
		PRINT_SUMLINE("shear_min.	 = %10.3e N/m\t@x= %5.1f km", shearmin, xshearmin/1000);
		PRINT_SUMLINE("deflection_max.= %8.1f m   \t@x= %5.1f km", wmaxdt, (x0+iwmaxdt*dx)/1000);
		PRINT_SUMLINE("deflection_min.= %8.1f m   \t@x= %5.1f km", wmindt, (x0+iwmindt*dx)/1000);
		PRINT_SUMLINE("defl.vel.max.  = %8.1f m/My\t@x= %5.1f km", dwmaxdt/(dt/Matosec), (x0+idwmaxdt*dx)/1000);
		PRINT_SUMLINE("defl.vel.min.  = %8.1f m/My\t@x= %5.1f km", dwmindt/(dt/Matosec), (x0+idwmindt*dx)/1000);
		PRINT_SUMLINE("first zero @x= %8.1f km ", xfirstnodo/1000);
		PRINT_SUMLINE("deflectn.area  = %8.1f km2 \tmean_rigid.= %.2e N m/m", Warea/1e6, Dmean);
	}
}



int gravity_anomaly()
{
	/*
	USES gravanompolyg() TO CALCULATE GRAVITY ANOMALY ALONG HORIZONS PROFILE 
			(results in m/s2) 
	USES geoidanompolyg() TO CALCULATE GEOID ANOMALY 
			(results in m) 
		z-values relative to sea level. This means that zero 
		deflection is at an elevation of 'zini'.
	*/

	register int   	i, ix, k_Block, i_Block, np_aux_Block, ix_min, ix_max;
	float 	z_max_grav_model = -zini+crust_thick_default+40000, 
		*Block_aux_x, *Block_aux_z, 
		*upper_hori_aux, *lower_hori_aux, 
		*alt_measurement,
		aux;
	float	*gravanom, 		/*Gravity anomaly	[mgal]*/
		*geoidanom;		/*Geoid anomaly		[m]   */
	FILE 	*file;
	char 	filename[MAXLENFILE];

	sprintf(filename,"%s.xg", projectname);
	remove(filename);

	if (!grav_anom_type) return(0);

	fprintf(stdout, " a");	fflush(stdout);

	geoidanom =  (float *) calloc(Nx, sizeof(float));
	gravanom =   (float *) calloc(Nx, sizeof(float));
	Block_aux_x = (float *) calloc (Nx*2+5, sizeof(float));
	Block_aux_z = (float *) calloc (Nx*2+5, sizeof(float));
	alt_measurement = (float *) calloc (Nx, sizeof(float));
	lower_hori_aux = (float *) calloc (Nx, sizeof(float));
	upper_hori_aux = (float *) calloc (Nx, sizeof(float));
	for (ix=0; ix<Nx; ix++) {gravanom[ix] = 0;  geoidanom[ix] = 0;}

	ix_min = MAX_2((xmin-x0-.1*dx)/dx, 0) ;	ix_max = MIN_2(floor((xmax-x0+.1*dx)/dx) + 2, Nx);

	Repare_Blocks();

	sprintf(filename,"%s.grv_mod", projectname);
	remove(filename);
	if ((file = fopen(filename, "wt")) == NULL) {
		if (verbose_level>=3) fprintf(stderr, "Warning: Cannot open output file '%s'.\n", filename);
	}
	fprintf(file, "#Time: %.2fMy\tGravity model bodies: x(km)-z(m)", Time/Matosec);

	/*Calculates topography at every point*/
	calculate_topo(topo);
	for (ix=0; ix<Nx; ix++)  alt_measurement[ix] = MAX_2(topo[ix], 0) + 10;

	/*Calculates gravity atraction of the Blocks*/
	np_aux_Block=Nx*2+4;
	for (i_Block=0; i_Block<numBlocks; i_Block++) {
		/*Constructs the auxiliar Block*/
		for (ix=0; ix<Nx; ix++) {
			upper_hori_aux[ix] = w[ix]-Blocks_base[ix];
			for (k_Block=0; k_Block<=i_Block; k_Block++) upper_hori_aux[ix] -= Blocks[k_Block].thick[ix];
			lower_hori_aux[ix] = w[ix]-Blocks_base[ix];
			for (k_Block=0; k_Block<i_Block; k_Block++)  lower_hori_aux[ix] -= Blocks[k_Block].thick[ix];
		}
		make_gravi_body (upper_hori_aux, lower_hori_aux, Block_aux_x, Block_aux_z);
		fprintf(file, "\n>");
		for (i=0; i<np_aux_Block; i++) fprintf(file, "\n%9.2f %9.1f", Block_aux_x[i]/1000, -Block_aux_z[i]);

		/*Calculates anomaly for this Block*/
		for (ix=ix_min; ix<ix_max; ix++) {
			gravanom[ix]  += gravanompolyg (Block_aux_x, Block_aux_z, np_aux_Block, x0+ix*dx, -alt_measurement[ix], Blocks[i_Block].density);
			geoidanom[ix] += geoidanompolyg(Block_aux_x, Block_aux_z, np_aux_Block, x0+ix*dx, -10e3, Blocks[i_Block].density);
		}
	}

	/*Calculates gravity atraction of sea water*/
	if (water_load) {
		for (ix=0; ix<Nx; ix++) {
			upper_hori_aux[ix] = 0;
			lower_hori_aux[ix] = MAX_2(-topo[ix], 0);
		}
		make_gravi_body (upper_hori_aux, lower_hori_aux, Block_aux_x, Block_aux_z);
		fprintf(file, "\n>");
		for (i=0; i<np_aux_Block; i++) fprintf(file, "\n%9.2f %9.1f", Block_aux_x[i]/1000, -Block_aux_z[i]);

		for (ix=ix_min; ix<ix_max; ix++) {
			gravanom[ix]  += gravanompolyg (Block_aux_x, Block_aux_z, np_aux_Block, x0+ix*dx, -alt_measurement[ix], denswater);
			geoidanom[ix] += geoidanompolyg(Block_aux_x, Block_aux_z, np_aux_Block, x0+ix*dx, -10e3, denswater);
		}
	}

	/*Calculates gravity atraction of crust*/
	for (ix=0; ix<Nx; ix++) {
		upper_hori_aux[ix] = w[ix]-Blocks_base[ix];
		lower_hori_aux[ix] = w[ix]+crust_thick[ix];
	}
	make_gravi_body (upper_hori_aux, lower_hori_aux, Block_aux_x, Block_aux_z);
	fprintf(file, "\n>");
	for (i=0; i<np_aux_Block; i++) fprintf(file, "\n%9.2f %9.1f", Block_aux_x[i]/1000, -Block_aux_z[i]);

	for (ix=ix_min; ix<ix_max; ix++) {
		gravanom[ix]  += gravanompolyg (Block_aux_x, Block_aux_z, np_aux_Block, x0+ix*dx, -alt_measurement[ix], denscrust);
		aux=geoidanompolyg(Block_aux_x, Block_aux_z, np_aux_Block, x0+ix*dx, -10e3, denscrust);
		geoidanom[ix] += aux;
	}

	/*Calculates gravity atraction of mantle*/
	for (ix=0; ix<Nx; ix++) {
		upper_hori_aux[ix] = w[ix]+crust_thick[ix] - zini;
		lower_hori_aux[ix] = z_max_grav_model;
	}
	make_gravi_body (upper_hori_aux, lower_hori_aux, Block_aux_x, Block_aux_z);
	fprintf(file, "\n>");
	for (i=0; i<np_aux_Block; i++) fprintf(file, "\n%9.2f %9.1f", Block_aux_x[i]/1000, -Block_aux_z[i]);

	for (ix=ix_min; ix<ix_max; ix++) {
		aux=gravanompolyg (Block_aux_x, Block_aux_z, np_aux_Block, x0+ix*dx, -alt_measurement[ix], densmantle);
		gravanom[ix]  += aux;
		geoidanom[ix] += geoidanompolyg(Block_aux_x, Block_aux_z, np_aux_Block, x0+ix*dx, -10e3, densmantle);
	}

	/*Substracts the atraction of normal water layer*/
	if (water_load && zini<0) {
		Block_aux_x[0] = x0 - 1e8;	Block_aux_z[0] = - zini;
		Block_aux_x[1] = x0 - 1e8;	Block_aux_z[1] = 0;
		Block_aux_x[2] = xf + 1e8;	Block_aux_z[2] = 0;
		Block_aux_x[3] = xf + 1e8;	Block_aux_z[3] = - zini;
		np_aux_Block = 4;
		fprintf(file, "\n>");
		for (i=0; i<np_aux_Block; i++) fprintf(file, "\n%9.2f %9.1f", Block_aux_x[i]/1000, -Block_aux_z[i]);

		for (ix=ix_min; ix<ix_max; ix++) {
			gravanom[ix]  -= gravanompolyg (Block_aux_x, Block_aux_z, np_aux_Block, x0+ix*dx, -alt_measurement[ix], denswater);
			geoidanom[ix] -= geoidanompolyg(Block_aux_x, Block_aux_z, np_aux_Block, x0+ix*dx, -10e3, denswater);
		}
	}

	/*Substracts the atraction of normal crust*/
	Block_aux_x[0] = x0 - 1e8;	Block_aux_z[0] = - zini + crust_thick_default;
	Block_aux_x[1] = x0 - 1e8;	Block_aux_z[1] = - zini;
	Block_aux_x[2] = xf + 1e8;	Block_aux_z[2] = - zini;
	Block_aux_x[3] = xf + 1e8;	Block_aux_z[3] = - zini + crust_thick_default;
	np_aux_Block = 4;
	fprintf(file, "\n>");
	for (i=0; i<np_aux_Block; i++) {fprintf(file, "\n%9.2f %9.1f", Block_aux_x[i]/1000, -Block_aux_z[i]);}

	for (ix=ix_min; ix<ix_max; ix++) {
		aux = gravanompolyg (Block_aux_x, Block_aux_z, np_aux_Block, x0+ix*dx, -alt_measurement[ix], denscrust);
		gravanom[ix]  -= aux;
		geoidanom[ix] -= geoidanompolyg(Block_aux_x, Block_aux_z, np_aux_Block, x0+ix*dx, -10e3, denscrust);
	}

	/*Substracts the atraction of normal mantle*/
	Block_aux_x[0] = x0 - 1e8;	Block_aux_z[0] = z_max_grav_model;
	Block_aux_x[1] = x0 - 1e8;	Block_aux_z[1] = - zini + crust_thick_default;
	Block_aux_x[2] = xf + 1e8;	Block_aux_z[2] = - zini + crust_thick_default;
	Block_aux_x[3] = xf + 1e8;	Block_aux_z[3] = z_max_grav_model;
	np_aux_Block = 4;
	fprintf(file, "\n>");
	for (i=0; i<np_aux_Block; i++) fprintf(file, "\n%9.2f %9.1f", Block_aux_x[i]/1000, -Block_aux_z[i]);

	for (ix=ix_min; ix<ix_max; ix++) {
		aux = gravanompolyg (Block_aux_x, Block_aux_z, np_aux_Block, x0+ix*dx, -alt_measurement[ix], densmantle);
		gravanom[ix]  -= aux;
		geoidanom[ix] -= geoidanompolyg(Block_aux_x, Block_aux_z, np_aux_Block, x0+ix*dx, -10e3, densmantle);
	}


	/* Now we add Bouguer correction:
		Free-air correction is not applied since in a 2D model there are no
		vertical variations of the regional gravity field related to the 
		distance to the center of earth (Turcotte & Schubert, 1982).
	*/
	if (grav_anom_type==1) {
		for (ix=0; ix<Nx; ix++) if (x0+ix*dx >= xmin && x0+ix*dx <= xmax) {
			/*Assumed standar 2.67 g/cm3*/
			if (topo[ix]>sea_level) 
				gravanom[ix] -= 2*pi * 2670 * CGU * topo[ix];
			if (topo[ix]<sea_level) 
				gravanom[ix] -= 2*pi*(2670 - ((water_load)? 1000:0))* CGU * topo[ix];
		}
	}

	write_file_grav_anom(gravanom, geoidanom);

	free(Block_aux_x);
	free(Block_aux_z);
	free(alt_measurement);
	free(upper_hori_aux);
	free(lower_hori_aux);
	fprintf(file, "\n>");
	fclose(file);

	return(1);
}




int move_Blocks()
{
	int	*nshift;
	float	*new_thick;

	/*
	  Moves the Blocks and calculates the isostatic load and thickness change
	  Deforms the sediment.
	*/

	new_thick = calloc(Nx, sizeof(float));
	nshift = calloc(numBlocks, sizeof(int));

	/*CRAWL UP along all Blocks to find the relevant moving ones first*/
	for (int iu=0; iu<numBlocks; iu++) {
		PRINT_DEBUG("move_Block:  Block %d; nshift=%d", iu, nshift[iu]);
		if (Blocks[iu].density == denssedim) {
		/*DEFORM SEDIMENT UNITS*/
		for (int i=0; i<Nx; i++) new_thick[i] = Blocks[iu].thick[i];
		for (int i=0; i<Nx; i++) {
			float sedthick;
			/*CRAWL DOWN Blocks to find the uppermost moving Block below this point*/
			for (int ju=iu-1, sedthick=Blocks[iu].thick[i]; ju>=0; ju--) {
			/*Calculate the thickness of sediments between the top of this sed. Block and the moving Block*/
			if (Blocks[ju].density == denssedim) {
				sedthick += Blocks[ju].thick[i];
			}
			else {
   			  /*Amount of cells to propagate the deformation: ~20 deg assumed.*/
			  int nprop = SIGN(nshift[ju]) * (int) ceil(sedthick*5/dx);
			  int i_unprop = i-nprop;
			  DOMAIN_LIMIT_1D(i_unprop);
			  if (Blocks[ju].thick[i_unprop]>.1) {
				if (!nshift[ju]) 
					break;
				else {
					int i_shift =   i+nshift[ju];
					/*If block ju is moving below [i] then shift seds.*/
					if (deform_sed && IN_DOMAIN_1D(i_shift)) 
					new_thick[i_shift] += Blocks[iu].thick[i];
					if (deform_sed && IN_DOMAIN_1D(i))  
					new_thick[i]	   -= Blocks[iu].thick[i];
					break;
				}
			  }
				}
			}
		}
		}
		else 
		{
		/*MOVE BLOCK UNITS and define nshift[]*/
		float theor_shift = Blocks[iu].vel * (Time-Blocks[iu].last_vel_time);
		nshift[iu] = floor((theor_shift - Blocks[iu].last_shift) /dx +.5);
		if (Time > Blocks[iu].time_stop + .1*dt) {nshift[iu]=0;}
		Blocks[iu].shift += nshift[iu]*dx;
		Blocks[iu].last_shift += nshift[iu]*dx;
		for (int i=0; i<Nx; i++) {
			int i_unshifted = i-nshift[iu];	
			DOMAIN_LIMIT_1D(i_unshifted);
			new_thick[i] = Blocks[iu].thick[i_unshifted];
		}
		}
		for (int i=0; i<Nx; i++) {
			if (new_thick[i]<0) PRINT_ERROR("\aBlock %d has a negative thickness: %.2f m", iu, new_thick[i]);
				Dq[i] += g * (new_thick[i] - Blocks[iu].thick[i]) * Blocks[iu].density;
				Blocks[iu].thick[i] = new_thick[i];
		}
	}
	fflush(stdout);
	free(new_thick);
	free(nshift);
	return(1);
}




int read_file_unit()
{
	/*
	  READS UNIT FILE NAMED 'projectnameNUM.UNIT' WHERE 'NUM' IS 1 FOR THE
	  FIRST UNIT, 2 FOR THE SECOND, ETC. Interpolates this input.
	  Creates new Block to store its properties and cuts sediment
	  Blocks when file contains fault depth rather than a thickness itself.
	*/

	int 	cut_Block, nparams=0;
	float	time_stop=9999/*My*/, time_unit, 
		erodibility_aux=NO_DATA, fill_up_to=NO_DATA, 
		vel=0, density=NO_DATA;
	BOOL 	insert, cut_seds, cut_Blocks, cut_all, top, fault, switch_move, 
		ride, hidden, z_absol;
	FILE 	*file;
	char 	filename[MAXLENFILE];

	/*Read the next unit age*/
	sprintf(filename, "%s%d.UNIT", projectname, nloads+1);
	if ((file = fopen(filename, "rt")) == NULL) {
		PRINT_WARNING("Cannot read unit file '%s'.", filename);
		return (0);
	}
	time_unit = Timeini/Matosec;
	{
		int nlines=0, nread, show, replace=0;
		char str1[MAXLENLINE], str2[MAXLENLINE], line[MAXLENLINE+200], *lineptr;
		show=(verbose_level>=3)? 1 : 0;
		rewind(file);
		while ((lineptr=fgets(line, MAXLENLINE+200-1, file)) != NULL && nlines<NMAXHEADERLINES) {
		nlines++; nread=sscanf(lineptr, "%s %s", str1, str2);
		if (nread == 2) {
			Match_Param_Replace_flt ( "time",  	time_unit, 0 )
			/*Old versions:*/
		}
		}
	}
	time_unit *= Matosec;
	/*Return if it isn't time yet to read the new unit file*/
	if (time_unit>Time+.1*dt || time_unit<Timeini) return(0);

	PRINT_INFO("Reading '%s'", filename);
	switch_move = fault = switch_gradual = 
		insert = hidden = cut_seds = cut_Blocks = cut_all = 
		top = ride = z_absol = NO;
	i_Block_insert = numBlocks;
	cut_Block = 0;

	/*READS AND INTERPOLATES UNIT/LOAD FILE*/
	{
		int nlines=0, nread, show, replace=0;
		char str1[MAXLENLINE], str2[MAXLENLINE], line[MAXLENLINE+200], *lineptr;
		rewind(file); 
		while ((lineptr=fgets(line, MAXLENLINE+200-1, file)) != NULL && nlines<NMAXHEADERLINES) {
				nlines++; nread=sscanf(lineptr, "%s %s", str1, str2);
				if (nread == 2) {
			Match_Param_Replace_flt ( "vel",		vel, 	0 )
			Match_Param_Replace_flt ( "time_stop",  	time_stop, 	0 )
			Match_Param_Replace_flt ( "density", 	 	density, 	0 )
			Match_Param_Replace_flt ( "erodibility",	erodibility_aux, 	0 )
			Match_Param_Replace_int ( "gradual",		switch_gradual, 	0 )
			Match_Param_Replace_int ( "hidden",		hidden, 	0 )
			Match_Param_Replace_int ( "ride",		ride, 	0 )
			Match_Param_Replace_int ( "insert",		insert, 	0 )
			Match_Param_Replace_int ( "top",		top, 	0 )
			Match_Param_Replace_int ( "move",  		switch_move, 	0 )
			Match_Param_Replace_int ( "fault",	  	fault, 	0 )
			Match_Param_Replace_int ( "z_absol",	  	z_absol,   	0 )
			Match_Param_Replace_int ( "cut_seds",  		cut_seds, 	0 )
			Match_Param_Replace_int ( "cut_Block",  		cut_Block, 	0 )
			Match_Param_Replace_int ( "cut_Blocks",  	cut_Blocks, 	0 )
			Match_Param_Replace_int ( "cut_all",  		cut_all, 	0 )
			Match_Param_Replace_int ( "topoest",		switch_topoest, 	0 )
			Match_Param_Replace_flt ( "fill_up_to", 	fill_up_to, 	0 )
			/*Old versions:*/
			Match_Param_Replace_int ( "cut_unit",  		cut_Block, 	0 )
			Match_Param_Replace_int ( "cut_units",  	cut_Blocks, 	0 )
			Match_Param_Replace_flt ( "erodability",	erodibility_aux, 	1 )
				}
				if (strcmp(str1, "thickness_distribution")==0) break;
		}
		rewind(file); 
	}
	if (fill_up_to == NO_DATA) 
		readinterplin(file, h_last_unit, Nx, x0, xf) ;
	else {
		for (int i=0; i<Nx; i++) h_last_unit[i] = MAX_2(0, fill_up_to-topo[i]);
	}
	fclose(file);

	nloads++;

	/*Check incompatibilities between unit file signals*/
	if (switch_gradual && switch_move) {
		PRINT_WARNING("Gradual+moving units not implemented. It won't be gradual.");
		switch_gradual = NO;
	}

	vel *= 1e3/Matosec;
	time_stop *= Matosec;

	/*ACT ACCORDING TO THE SIGNALS*/
	if (fault) {
		switch_move = YES;
	}
	/*Creates a Block of infill if switch_topoest; it will be filled during Deflection*/
	if (switch_topoest) {
		insert_new_Block(i_first_Block_load);
		Blocks[i_first_Block_load].type = 'I'; /*stands for Infill*/
		Blocks[i_first_Block_load].density = densinfill;
		if (densinfill<2550) Blocks[i_first_Block_load].erodibility = erodibility_sed;
		i_first_Block_load++; 	i_Block_insert++;
	}
	if (insert) {
		i_Block_insert = 0;
	}
	if (top) {
		for (int k=numBlocks-1; k>=0; k--) {
			if (Blocks[k].density != denssedim) {
				i_Block_insert = k+1;
				break;
			}
		}
	}
	if (cut_all) {
		cut_Blocks = YES;
	}

	if (fault && !top && !cut_all) i_Block_insert = 0;

	if (fault && !cut_Block) {
		int numBlocks0=numBlocks;
		/*Make a copy of all Blocks*/
		PRINT_DEBUG("Cutting Blocks: numBlocks= %d", numBlocks);
		for (int k=0; k<numBlocks0; k++) {
			float *thick_aux;
			insert_new_Block(numBlocks);
			thick_aux = Blocks[numBlocks-1].thick;
			Blocks[numBlocks-1] = Blocks[k];
			Blocks[numBlocks-1].thick = thick_aux;
			Blocks[numBlocks-1].vel = vel;
			Blocks[numBlocks-1].last_vel_time = Time-dt;/*!!*/
			Blocks[numBlocks-1].last_shift = 0;
			Blocks[numBlocks-1].time_stop = time_stop;
			if (Blocks[numBlocks-1].type == 'S') {
				Blocks[numBlocks-1].detr_ratio  = calloc(Nx, sizeof(float));
				Blocks[numBlocks-1].detr_grsize = calloc(Nx, sizeof(float));
			}
			if (density		 != NO_DATA && Blocks[numBlocks-1].type != 'S') Blocks[numBlocks-1].density	 = density;
			if (erodibility_aux != NO_DATA && Blocks[numBlocks-1].type != 'S') Blocks[numBlocks-1].erodibility = erodibility_aux;
		}
		/*Modify Blocks_base and cut above the fault*/
		for (int i=0; i<Nx; i++) {
			float z_fault = -h_last_unit[i], base_of_Block = Blocks_base[i];
			if (z_absol) base_of_Block -= w[i];
			h_last_unit[i] = MAX_2(0, Blocks_base[i] - z_fault);
			Blocks_base[i]  = MIN_2(Blocks_base[i], z_fault);
			if (cut_Blocks) {
				for (int k=0; k<numBlocks0; k++) {
						float top_of_Block=base_of_Block+Blocks[k].thick[i];
					if (Blocks[k].density == denssedim && !cut_all) {
						break;
					}
					if (z_fault <= base_of_Block) {
						Blocks[k+numBlocks0].thick[i] += Blocks[k].thick[i];
						Blocks[k].thick[i]			= 0;
					}
					else {
						Blocks[k+numBlocks0].thick[i] += MAX_2 (0, top_of_Block-z_fault);
						Blocks[k].thick[i]		   -= MAX_2 (0, top_of_Block-z_fault);
					}
						base_of_Block = top_of_Block;
				}
			}
		}
	}
	if (density		==NO_DATA) density		 = denscrust;
	if (erodibility_aux==NO_DATA) erodibility_aux = erodibility;

	/*Create a new Block for this file*/
	insert_new_Block(i_Block_insert);

	/*Add the thickness in file to the new Block; Thin Blocks and basement if the thickness is negative*/
	if (!switch_gradual && !hidden && !cut_Block) {
		for (int i=0; i<Nx; i++) {
		if (h_last_unit[i]>0) 
			Blocks[i_Block_insert].thick[i] = h_last_unit[i];
		else {
			float   h_load_aux, h_load_aux2;
			int	 j;
			h_load_aux = fabs(h_last_unit[i]);
			for (j=i_Block_insert-1; h_load_aux>0 && j>=0; j--) {
					h_load_aux2 = MIN_2(Blocks[j].thick[i], h_load_aux);
					h_load_aux -= h_load_aux2;
					Blocks[j].thick[i] -= h_load_aux2;
			}
			/*j is the deepest eroded Block in i*/
			if (j==-1) {
					Blocks_base[i] -= h_load_aux;
			}
		}
		}
	}
	if (hidden) Blocks[i_Block_insert].type = 'H';
	Blocks[i_Block_insert].density = density;
	Blocks[i_Block_insert].erodibility = erodibility_aux;
	Blocks[i_Block_insert].vel = vel;
	Blocks[i_Block_insert].time_stop = time_stop;

	if (ride) {
		PRINT_DEBUG("ride %d %d", i_Block_insert, numBlocks);
		for (int i_Block=i_Block_insert+1; i_Block<numBlocks; i_Block++) {
			Blocks[i_Block].vel		   = Blocks[i_Block_insert].vel; 
			Blocks[i_Block].last_shift	= 0; 
			Blocks[i_Block].last_vel_time = Time-dt;/*!!*/ 
			Blocks[i_Block].time_stop	 = Blocks[i_Block_insert].time_stop; 
		}
	}

	if (cut_Block) {
		int numBlocks0=numBlocks;
		float *thick_aux;
		if (cut_Block<0) {
			/*Search the biggest Block with density close to -cut_Block*/
			float vol, volmax=0, Blockvolmax=-1;
			for (int k=0; k<numBlocks; k++) {
				for (int i=vol=0; i<Nx; i++)  vol += Blocks[k].thick[i];
				if (vol>volmax && fabs(Blocks[k].density-fabs(cut_Block))<.5) {volmax=vol; Blockvolmax=k;}
			}
			cut_Block = Blockvolmax+1;
		}
		PRINT_INFO("Cutting Block %d", cut_Block);
		insert_new_Block(cut_Block);
		thick_aux = Blocks[cut_Block].thick;
		Blocks[cut_Block] = Blocks[cut_Block-1];
		Blocks[cut_Block].thick = thick_aux;
		Blocks[cut_Block].vel = vel;
		Blocks[cut_Block].density = density;
		Blocks[cut_Block].last_vel_time = Time-dt;/*!!*/
		Blocks[cut_Block].last_shift = 0;
		Blocks[cut_Block].time_stop = time_stop;
		/*Cut the Block*/
		for (int i=0; i<Nx; i++) {
			float z_fault = -h_last_unit[i], base_of_Block = Blocks_base[i], top_of_Block;
			if (z_absol) base_of_Block -= w[i];
			for (int k=0; k<cut_Block-1; k++) base_of_Block += Blocks[k].thick[i];
			top_of_Block = base_of_Block + Blocks[cut_Block-1].thick[i];
			if (z_fault <= base_of_Block) {
				Blocks[cut_Block].thick[i]   += Blocks[cut_Block-1].thick[i];
				Blocks[cut_Block-1].thick[i]  = 0;
			}
			else {
				Blocks[cut_Block].thick[i]   += MAX_2 (0, top_of_Block-z_fault);
				Blocks[cut_Block-1].thick[i] -= MAX_2 (0, top_of_Block-z_fault);
			}
		}
	}

	/*Don't Repare_Blocks() in case of: 
		Gradual load, because then h_last_unit[] will be empty until tectload()
		Topoest load, because the infill Block will be filled upon deflection.
	*/
	if (!switch_gradual && !switch_topoest) Repare_Blocks();

	/*Increment the isostatic load for this time interval*/
	if (!switch_gradual && !fault && !cut_Block) /*!cut_Block condition needed for Alice's setup (redefining densities without 'fault' option*/
		for (int i=0; i<Nx; i++) Dq[i] += (density-densenv)*g*h_last_unit[i];

	PRINT_INFO("Unit read from '%s'. ", filename);
	PRINT_DEBUG("%d params; dens=%.0f kg/m3; erodibility=%.1e; ", nparams, density, erodibility_aux);
	if (switch_gradual) PRINT_INFO("Will be gradually loaded until %.2fMy.", time_stop/Matosec);
	if (switch_move) PRINT_INFO("Vel= %.2f km/My till T=%.1f My", vel*Matosec/1000, time_stop/Matosec);

	fprintf(stdout, " l") ;
	if (fault)  		fprintf(stdout, "F") ;
	if (hidden)		fprintf(stdout, "H") ;
	if (insert) 		fprintf(stdout, "I") ;
	if (top) 		fprintf(stdout, "P") ;
	if (switch_gradual)	fprintf(stdout, "G") ;
	if (switch_move) 	fprintf(stdout, "M") ;
	if (ride) 		fprintf(stdout, "R") ;
	if (switch_topoest)	fprintf(stdout, "T") ;

	return(1);
}




int syntax()
{
	/*
		Displays the command line syntax of the program
	*/
	char 	filename[MAXLENFILE], line[MAXLENLINE], *lineptr;
	FILE 	*file;
	BOOL 	print=NO;

	sprintf(filename, "%s/doc/tao.info.txt", TAODIR);
	if ((file = fopen(filename, "rt")) == NULL) {
		if (verbose_level>=3) fprintf(stderr, "\nWarning: Cannot find file %s. You shouldn't have moved tAo directory after compilation...", filename);
		return(0);
	}
	fprintf(stderr, "\nSyntax:\n");
	while (1) {
		fgets(line, MAXLENLINE-1, file);
		if (lineptr=strstr(line, "Signals")) break;;
		if (print==YES) fprintf(stderr, "%s", line);
		if (strstr(line, "SYNTAX")) print=YES;
	}
	fclose(file);
	return (1);
}





int The_End()
{
	int	i, j;
	char	command[MAXLENLINE];	
	float	total_load=0, total_restitutive_force=0, total_Blocks_mass=0, total_sed_mass=0, total_sed_grain_mass=0, 
		x, xleft, xright, Krest;

	Krest =  (switch_topoest) ?  (densasthen-densinfill)*g : (densasthen-densenv)*g ;

	gravity_anomaly();
	write_file_Temperature();

	if (verbose_level>=1) {
		fprintf(stdout, "\n\nFinal statistics:");
		for (i=0; i<Nx; i++)  {total_load += q[i]; total_restitutive_force += (Krest*w[i]);}
		total_load *= (dx*Nx/(Nx+1));
		total_restitutive_force *= (dx*Nx/(Nx+1));

		for (i=0; i<numBlocks; i++)  for (j=0; j<Nx; j++) {
			total_Blocks_mass += Blocks[i].thick[j] * Blocks[i].density;
			if (Blocks[i].density == denssedim) {
				total_sed_mass += Blocks[i].thick[j] * Blocks[i].density;
				total_sed_grain_mass += THICK2SEDMASS_1D(Blocks[i].thick[j]);
			}
		}
		total_Blocks_mass *= (dx*Nx/(Nx+1));
		total_sed_mass  *= (dx*Nx/(Nx+1));
		total_sed_grain_mass  *= ((float) Nx/ (float) (Nx+1));
		fprintf(stdout, "\n\tTotal_load_weight = %10.3e N/m", total_load);
		fprintf(stdout, "\n\tTotal_rest._force = %10.3e N/m", total_restitutive_force);
		fprintf(stdout, "\n\tTotal_Blocks_mass = %10.3e kg/m", total_Blocks_mass);
		fprintf(stdout, "\n\tTotal_sedim_mass  = %10.3e kg/m (%.3g kg/m of grain)", total_sed_mass, total_sed_grain_mass);
	}

	fprintf(stdout, "\n\n%d Blocks:", numBlocks);
	fprintf(stdout, "\nNo.\tDensity\tAge \tVolume \tVel   \tShift \tFrom x\tTo x \tStop\ttype");
	if (erosed_model) fprintf(stdout, "\tErodabil.");
	fprintf(stdout, "\n   \t(kg/m3)\t(My)\t(km2) \t(km/My)\t(km)  \t(km)  \t(km) \t (My)\t(-)");
	if (erosed_model) fprintf(stdout, "\t		 ");
	for (i=numBlocks-1; i>=0; i--) {
		float volume;
		for (xleft=x0,j=0; j<Nx; j++) {
			x=x0+j*dx; 
			if (Blocks[i].thick[j]>1) break;
			xleft=x;
		}
		for (xright=xf,j=Nx-1; j>=0; j--) {
			x=x0+j*dx; 
			if (Blocks[i].thick[j]>1) break;
			xright=x;
		}
		for (volume=j=0; j<Nx; j++) {
			volume += Blocks[i].thick[j];
		}
		volume *= dx*Nx/(Nx+1);
		fprintf(stdout, "\n%2d:\t%.0f \t%.1f\t%.1f\t%.2f\t%.1f\t%.1f\t%.1f\t%.1f\t%c ", 
			i, Blocks[i].density, Blocks[i].age/Matosec, volume/1e6, 
			Blocks[i].vel/1e3*Matosec, Blocks[i].shift/1e3, xleft/1e3, xright/1e3, Blocks[i].time_stop/Matosec, Blocks[i].type); 
		if (erosed_model) fprintf(stdout, "\t%.2e", Blocks[i].erodibility);
	}
	fprintf(stdout, "\n -:\t%.0f\t%.1f\t-\t0\t-\t-\t-\t-\tbasement\n", denscrust, Timeini/Matosec);

	write_file_resume();
	if (!switch_dt_output) Write_Ouput();

	sprintf(command, "rm -f %s*.tao.tmp", projectname);
	system(command);
	{
		char filename[MAXLENFILE];
		FILE *file;
		sprintf(filename, "%s/.taodefaults", getenv ("HOME"));
		if ((file = fopen(filename, "rt")) == NULL) {
			sprintf(command, "echo First use of %s by `whoami` at `hostname` > %s; date >> %s", version, filename, filename);
			if (verbose_level>=3) fprintf(stdout, "\n%s", command);
			system(command);
			sprintf(command, "mail danielgc@ictja.csic.es < %s", filename);
			if (verbose_level>=3) fprintf(stdout, "\n%s", command);
			system(command);
			fclose (file);
		}
	}
	if (verbose_level>=3) AUTHORSHIP;
	fprintf(stdout, "\n\n");

	if (run_type==10 || run_type==2) exit(0);
	return(1);
}




int Viscous_Relaxation()
{
	int 	i, NDs=3, NDi=3;
	double	**A, *b;
	float	*dwdt, *moment;

	/*
	CALL SUBROUTINES TO CALCULATE VISCOELASTIC RELAXATION FROM PRESENT 
	DEFLECTION AND TOTAL LOAD.   ADD NEW DEFLECTION INCREMENT TO THE 
	PREVIOUS DEFLECTION
	*/

	if (isost_model!=2 || !Te_default)  return(0);

	A = alloc_matrix_dbl (Nx, NDi+1+NDs);
	b = (double *) calloc (Nx , sizeof(double));
	moment = (float *) calloc (Nx , sizeof(float));
	dwdt   = (float *) calloc (Nx , sizeof(float));

	fprintf(stdout, " v");
	LES_matrix(A, b, D, q, Dq, w, YES) ;
	solveLES(A, b, Nx, NDs, NDi, dwdt) ;
	for (i=0; i<Nx; i++) {
		w[i] += dwdt[i]*dt;	Dw[i] = dwdt[i]*dt;
	}
	for (i=1; i<Nx-1; i++) {
		moment[i] = (- D[i]*(1-nu*nu)/(1-.25)  * tau * 
				(dwdt[i+1] - 2*dwdt[i] + dwdt[i-1]) / 
				pow(dx, 2) + moment[i]*tau/dt ) / 
				(1 + tau/dt) ;
		/*If the deflection was made constant then:*/
		/*  	/= exp(dt/tau) ;   */
	}
	moment[0]=moment[1];
	moment[Nx-1]=moment[Nx-2];

	if (switch_topoest) {
		/*Defines the thickness of last infill Block*/
		for (i=0; i<Nx; i++) 
			Blocks[i_first_Block_load-1].thick[i] +=  MAX_2(Dw[i], 0) ;
	}

	flexural_stats(moment);

	free(moment);
	for (i=0; i<Nx; i++) free(A[i]);
	free(A);
	free(b);
	free(dwdt);
	return(1);
}



int Write_Ouput()
{
	write_file_time(w, topo);

	write_file_Blocks();
	write_file_erosed(total_erosion);
	write_file_Te();
	write_file_stress();
	write_file_maxmompoint();

	/*Make GMT Postscript*/
	if (switch_ps) {
		char 	command[200];
		sprintf(command, "tao.gmt.job %s %.2f %.2f %.2f %.2f %d %d", 
			projectname, xmin/1000, xmax/1000, zmin, zmax, water_load, (isost_model<3)? 0:1);
		if (verbose_level>=3) 
			fprintf(stdout, "\nPostscript file '%s.ps' is going to be produced with command:", projectname) ;
		if (verbose_level>=3) 
			fprintf(stdout, "\n%s\n", command) ;
		system(command);
		if (switch_dt_output) {
			/*crop by default to the border*/
			if (strlen(gif_geom)<2) sprintf(gif_geom, "-trim -background Khaki -label 'tAo software: %s' -gravity South -append", projectname);
			sprintf(command, "convert -density 200 %s.ps %s -interlace NONE  %s%03d.jpg", /*-fill \"#ffffff\" -draw \"rectangle 70,10 130,25\" -fill \"#000000\" -font helvetica -draw \"text 74,22 t_%+3.2f_My \" */
				projectname, gif_geom, projectname, n_image);
			if (verbose_level>=3)
				fprintf(stdout, "\n%s\n", command) ;
			system(command);
			n_image++;
		}
	}

	return (1);
}


