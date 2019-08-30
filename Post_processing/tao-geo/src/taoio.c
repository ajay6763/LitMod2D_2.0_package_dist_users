/*
INPUT/OUTPUT  FUNCTIONS  FOR  tao.c
*/

int read_file_horiz_record_time()
{
	/*
		Reads file containig times of horizons recording named 'projectname.REC'
	*/

	int 	i;
	FILE 	*file;
	char 	filename[MAXLENFILE];
	float	*aux1;

	Read_Open_Filename_Return(".REC", "rt", "Horizon recording times")

	n_record_times=0;
	aux1 = calloc(nmax_input_points, sizeof(float));
	for (;;) {
		TAKE_LINE_1(aux1[n_record_times]);
		n_record_times++;
		if (n_record_times>=nmax_input_points-1 ) {
			PRINT_ERROR("Too many points (>%d) in horizon recording times file.", n_record_times-1);
			break;
		}
	}
	fclose(file);
	if (verbose_level>=1) fprintf(stdout, "\nHorizon recording times at '%s'. %d times were read.", filename, n_record_times);

	horiz_record_time = calloc(n_record_times, sizeof(float));
	for (i=0; i<n_record_times; i++) {
		horiz_record_time[i] = aux1[i]*Matosec;
	}
	
	/*desactivate automatic sediment Blocks generation*/
	if (n_record_times) dt_record=0;
	
	free(aux1); 
	return(1);
}




int read_file_initial_topo(float *h)
{
	int 	i ;
	FILE 	*file ;

	/*READS INITIAL TOPOGRAPHY IF EXISTS FILE WITH NAME 'projectname.ZINI'*/

	Read_Open_Filename_Return(".ZINI", "rt", "Initial topography")

	readinterplin(file, h, Nx, x0, xf) ;

	fclose(file);
	return (1);
}



int read_file_initial_deflection(float *w)
{
	int 	i ;
	FILE 	*file ;

	/*READS INITIAL DEFLECTION IF DOES EXIST FILE WITH NAME 'projectname.WINI'*/

	Read_Open_Filename_Return(".WINI", "rt", "Initial deflection")

	readinterplin(file, w, Nx, x0, xf) ;

	fclose(file);
	return (1);
}




int read_file_Crust_Thick(float crust_thick_default)
{
	/*CARGA EL FICHERO DE ESPESOR CORTICAL, DE NOMBRE  'projectname.CRUST'
			Interpola entre los puntos del fichero.*/


	int 	i;
	FILE 	*file;

	for (i=0;i<Nx;i++) {crust_thick[i] = crust_thick_default;}

	Read_Open_Filename_Return(".CRUST", "rt", "Crustal thickness")

	readinterplin(file, crust_thick, Nx, x0, xf) ;

	fclose(file);
	return(1);
}



int read_file_Upper_Crust_Thick(float crust_thick_supdefault)
{
	/*CARGA EL FICHERO DE ESPESOR DE CORTEZA SUPERIOR, DE NOMBRE  'projectname.UCRUST'
			Interpola entre los puntos del fileero.*/


	int 	i;
	FILE 	*file;

	for (i=0;i<Nx;i++) {upper_crust_thick[i] = upper_crust_thick_default;}

	Read_Open_Filename_Return(".UCRUST", "rt", "Upper crustal thickness")

	readinterplin(file, upper_crust_thick, Nx, x0, xf) ;

	for (i=0; i<Nx; i++) upper_crust_thick[i] = MIN_2(upper_crust_thick[i], crust_thick[i]) ;

	fclose(file);
	return(1);
}



int read_file_parameters (int show, int reformat) 
{
	int	nread, nparams=0, nline=0, verbose_level_ant=verbose_level;
	char 	*lineptr, str1[MAXLENLINE], str2[MAXLENLINE], 
		line[MAXLENLINE+200], PRMfilename[MAXLENFILE];
	FILE 	*file;
	BOOL	switch_matched_vers=NO;

	/*
	READ THE PARAMETERS FILE NAMED  'projectname.PRM'
	You have an explanation of these parameters in the example project 
	parameters file (doc/template.PRM), where you have also the name of 
	the variable related to each parameter. The use of most variables is 
	described in the include files '*.h'.
	*/

	sprintf(PRMfilename, "%s.PRM", projectname);
	if (show && verbose_level>=3) fprintf(stdout, "\nCurrent tAo project: %s", projectname);
	if ((file = fopen(PRMfilename, "rt")) == NULL) {
		PRINT_ERROR("Can't open parameters file '%s'.\n", PRMfilename);
		return(0);
	}

	x0=NO_DATA; xf=NO_DATA;
	if (show) fprintf(stdout, "\nParameters at '%s'.", PRMfilename);
	while ((lineptr=fgets(line, MAXLENLINE+200-1, file)) != NULL) {
		int status;
		status=0;
		nread=sscanf(lineptr, "%s %s", str1, str2);
		if (nread >= 2) {
			if (!strcmp(str1, "version")) {
				strcpy(version_input, str2);
				if (!strcmp(version_input, version)) {
					switch_matched_vers = YES;
					nparams++;
				}
				if ((show && verbose_level>=3) || (reformat && verbose_level>=3)) PRINT_INFO("Input file version: %s\n", version_input);
				if (reformat) {fprintf(stdout, "version\t%s\n", version); }
			}
			status=match_parameter(str1, str2, show, reformat, line);
			nparams += status; 
			if (!strcmp(str1, "version")) status=1;
		}
		/*If no parameter matched then just reproduce the entire line*/
		if (reformat==1 && !status) fprintf(stdout, "%s", line);
		nline++;
	}
	if (verbose_level_ant>=2) fprintf(stdout, " (%d parameters matched)", nparams);
	if (!switch_matched_vers) {
		if (verbose_level_ant>=2) 
		    fprintf(stderr,
			"\nInfo: Present version '%s' not matched in PRM file. See 'tisc/doc/template.PRM'.", version);
		if (nparams<2) {
		    PRINT_ERROR("\aWrong format in parameters file '%s'. Only %d parameters were matched. Current version is '%s'. "
			"\nSee example file 'tisc/doc/template.PRM'.\nEND.\n", 
			PRMfilename, nparams, version);
		    exit(0); 
		}
	}

	if (x0==NO_DATA) x0 = xmin;
	if (xf==NO_DATA) xf = xmax;
	if (Kerosdif>1e5) PRINT_WARNING("Kerosdif expected in Blocks [m2/yr] !") ;

	fclose(file);
	return(1);
}




int read_file_resume(char *filename)
{
	int 	i, numBlocks_aux, i_Block_insert_aux, run_type_aux, end_check;
	FILE 	*file;
	char  	version_aux[LENGTHVERS], projectname_aux[MAXLENFILE];

	/*
	  READS A BINARY FILE WITH ALL THE INFORMATION 
	  REQUIRED TO RESTART THE PROGRAM
	*/

	if ((file = fopen(filename,"r")) == NULL) {PRINT_ERROR("Cannot read ""Resume"" input file '%s'.\n", filename); exit(0);}\
	if (verbose_level>=1) fprintf(stdout, "\n""Model-run data at '%s'", filename);

	/*Defined in universal.h:*/
	fread(&Nx, 		sizeof(int),		1, 	file);
	fread(&Nz, 		sizeof(int),		1, 	file);
	fread(&verbose_level, 	sizeof(int),		1, 	file);

	fread(&dx, 		sizeof(float),		1, 	file);
	fread(&dz, 		sizeof(float),		1, 	file);

	fread(version_aux,	sizeof(char),		LENGTHVERS, 	file);
	if (strcmp(version, version_aux)) PRINT_WARNING("restart file '%s' does not match present tAo version '%s'.", filename, version);
	fread(version_input,	sizeof(char),		LENGTHVERS, 	file);

	fread(&switch_geograph_coor, sizeof(BOOL),	1, 	file);
	fread(&switch_ps, 	sizeof(BOOL),		1, 	file);
	fread(&switch_write_file, sizeof(BOOL),		1, 	file);


	/*Defined in geomodel.h:*/
	fread(&grav_anom_type, 	sizeof(int),		1, 	file);
	fread(&isost_model, 	sizeof(int),		1, 	file);

	fread(&Te_default, 	sizeof(float),		1, 	file);
	fread(&crust_thick_default, sizeof(float),	1, 	file);
	fread(&upper_crust_thick_default, sizeof(float),1, 	file);
	fread(&densasthen, 	sizeof(float),		1, 	file);
	fread(&densmantle, 	sizeof(float),		1, 	file);
	fread(&denssedim, 	sizeof(float),		1, 	file);
	fread(&denscrust, 	sizeof(float),		1, 	file);
	fread(&densenv, 	sizeof(float),		1, 	file);
	fread(&densinfill, 	sizeof(float),		1, 	file);
	fread(&sea_level, 	sizeof(float),		1, 	file);
	fread(&temp_sea_level, 	sizeof(float),		1, 	file);
	fread(&Time, 		sizeof(float),		1, 	file);
	fread(&Timefinal, 	sizeof(float),		1, 	file);
	fread(&Timeini, 	sizeof(float),		1, 	file);
	fread(&dt, 		sizeof(float),		1, 	file);
	fread(&dt_eros, 		sizeof(float),		1, 	file);
	fread(&tau, 		sizeof(float),		1, 	file);

	fread(projectname_aux, 	sizeof(char),	MAXLENFILE, 	file);
	if (strcmp(projectname, projectname_aux)) fprintf(stdout, "\nERROR: restart file '%s' does not match present tAo project name '%s'.", filename, version);
	fread(gif_geom, 	sizeof(char),	MAXLENLINE, 	file);

	fread(&water_load, 	sizeof(BOOL),		1, 	file);


	/*Defined in tao+tisc.h:*/
	fread(&boundary_conds, 	sizeof(int),		1, 	file);
	fread(&nloads, 		sizeof(int),		1, 	file);
	fread(&nmax_input_points, sizeof(int),		1, 	file);
	fread(&n_sea_level_input_points, sizeof(int),	1, 	file);
	fread(&n_eros_level_input_points, sizeof(int),	1, 	file);
	fread(&n_record_times, 	sizeof(int),		1, 	file);
	fread(&i_first_Block_load, sizeof(int),	1, 	file);
	fread(&i_Block_insert_aux, sizeof(int),		1, 	file);
	fread(&numBlocks_aux, 	sizeof(int),		1, 	file);
	fread(&nwrotenfiles, 	sizeof(int),		1, 	file);
	fread(&run_type_aux, 	sizeof(int),		1, 	file);

	fread(&zini, 		sizeof(float),		1, 	file);
	fread(&dt_record, 	sizeof(float),		1, 	file);
	fread(&sed_porosity, 	sizeof(float),		1, 	file);
	fread(&compact_depth, 	sizeof(float),		1, 	file);
	fread(&Kerosdif, 	sizeof(float),		1, 	file);
	fread(&last_time_file_time, 	sizeof(float),		1, 	file);
	fread(&random_topo, 	sizeof(float),		1, 	file);

	fread(&switch_dt_output, 	sizeof(BOOL),		1, 	file);
	fread(&switch_file_out, 	sizeof(BOOL),		1, 	file);
	fread(&switch_gradual, 	sizeof(BOOL),		1, 	file);
	fread(&switch_insert_load, 	sizeof(BOOL),		1, 	file);
	fread(&switch_topoest, 		sizeof(BOOL),		1, 	file);
	fread(&switch_write_file_Blocks, sizeof(BOOL),		1, 	file);
	fread(&deform_sed, sizeof(BOOL),		1, 	file);

	/*Defined in tao.h:*/
	fread(&imomentmax, 	sizeof(int),		1, 	file);
	fread(&nx_temp_input, 	sizeof(int),		1, 	file);
	fread(&nbasins, 	sizeof(int),		1, 	file);
	fread(&nlakes,	 	sizeof(int),		1, 	file);
	fread(&n_image, 	sizeof(int),		1, 	file);
	fread(&hydro_model,	sizeof(int),		1, 	file);
	fread(&erosed_model,	sizeof(int),		1, 	file);
	fread(&eros_bound_cond,	sizeof(char),		2, 	file);

	fread(&x0, 		sizeof(float),		1, 	file);
	fread(&xf, 		sizeof(float),		1, 	file);
	fread(&xmin, 		sizeof(float),		1, 	file);
	fread(&xmax, 		sizeof(float),		1, 	file);
	fread(&zmin, 		sizeof(float),		1, 	file);
	fread(&zmax, 		sizeof(float),		1, 	file);
	fread(&horz_force, 	sizeof(float),		1, 	file);
	fread(&vert_force, 	sizeof(float),		1, 	file);
	fread(&appmoment, 	sizeof(float),		1, 	file);
	fread(&Keroseol, 	sizeof(float),		1, 	file);
	fread(&Ksedim, 		sizeof(float),		1, 	file);
	fread(&critical_slope, 	sizeof(float),		1, 	file);
	fread(&K_river_cap, 	sizeof(float),		1, 	file);
	fread(&erodibility, 	sizeof(float),		1, 	file);
	fread(&erodibility_sed, sizeof(float),		1, 	file);
	fread(&l_fluv_sedim, 	sizeof(float),		1, 	file);
	fread(&lost_rate, 	sizeof(float),		1, 	file);
	fread(&evaporation_ct, 	sizeof(float),		1, 	file);
	fread(&riverbasinwidth, 	sizeof(float),		1, 	file);
	fread(&rain, 		sizeof(float),		1, 	file);
	fread(&Krain,	 	sizeof(float),		1, 	file);
	fread(&CXrain, 		sizeof(float),		1, 	file);

	fread(&switch_strs_history, 	sizeof(BOOL),		1, 	file);
	fread(&switch_YSE_file, 	sizeof(BOOL),		1, 	file);


	/*Arrays:*/
	Allocate_Memory();
	fread(w, 		sizeof(float),	Nx, 	file);
	fread(D, 		sizeof(float),	Nx, 	file);
	fread(q, 		sizeof(float),	Nx, 	file);
	fread(Dw, 		sizeof(float),	Nx, 	file);
	fread(Dq, 		sizeof(float),	Nx, 	file);
	fread(h_water, 		sizeof(float),	Nx, 	file);
	fread(h_last_unit, 	sizeof(float),	Nx, 	file);
	fread(Te, 		sizeof(float),	Nx, 	file);
	fread(crust_thick, 	sizeof(float),	Nx, 	file);
	fread(upper_crust_thick,sizeof(float),	Nx, 	file);
	fread(topo, 		sizeof(float),	Nx, 	file);
	fread(Blocks_base, 	sizeof(float),	Nx, 	file);
	fread(yieldcompres, 	sizeof(float),	Nz, 	file);
	fread(yieldextens, 	sizeof(float),	Nz, 	file);

	horiz_record_time = calloc(n_record_times, sizeof(float));
	fread(horiz_record_time, sizeof(float), n_record_times, file);

	var_sea_level = calloc(n_sea_level_input_points, sizeof(float *));
	for (i=0; i<n_sea_level_input_points; i++) {
		var_sea_level[i] = calloc(2, sizeof(float));
		fread(var_sea_level[i], sizeof(float), 2, file);
	}
	var_eros_level = calloc(n_eros_level_input_points, sizeof(float *));
	for (i=0; i<n_eros_level_input_points; i++) {
		var_eros_level[i] = calloc(2, sizeof(float));
		fread(var_eros_level[i], sizeof(float), 2, file);
	}

	if (isost_model>=3 && !switch_YSE_file) {
		Temperature = alloc_matrix(Nx, Nz);
		for (i=0; i<Nx; i++) 
			fread(Temperature[i], sizeof(float), Nz, file);
	}
	if (isost_model>=3) {
		stress = alloc_matrix(Nx, Nz);
		for (i=0; i<Nx; i++) 
			fread(stress[i], sizeof(float), Nz, file);
	}


	for (i=0; i<numBlocks_aux; i++) {
		float *ptr1, *ptr2, *ptr3;
		insert_new_Block(numBlocks);
		ptr1 = Blocks[numBlocks-1].thick; 
		fread(&Blocks[numBlocks-1], sizeof(struct BLOCK_1D), 1, file);
		Blocks[numBlocks-1].thick = ptr1;
	}
	if (numBlocks_aux != numBlocks) PRINT_ERROR("in '%s', %d Blocks?!", filename, numBlocks_aux);
	i_Block_insert=i_Block_insert_aux;
	for (i=0; i<numBlocks; i++) {
		fread(Blocks[i].thick, 		sizeof(float),	Nx, 	file);
	}
	for (i=0; i<numBlocks; i++) {
	    if (Blocks[i].type == 'S') {
		Blocks[i].detr_ratio = calloc(Nx, sizeof(float));
		Blocks[i].detr_grsize = calloc(Nx, sizeof(float));
		fread(Blocks[i].detr_ratio, 	sizeof(float),	Nx, 	file);
		fread(Blocks[i].detr_grsize, 	sizeof(float),	Nx, 	file);
	    }
	}

	if (erosed_model) {
		fread(eros_now, 	sizeof(float),	Nx, 	file);
		fread(total_erosion, 	sizeof(float),	Nx, 	file);
	}
	if (hydro_model) {
		int j;
		fread(precipitation, 	sizeof(float),	Nx, 	file);
		fread(evaporation, 	sizeof(float),	Nx, 	file);
		Lake = calloc (nlakes+1, sizeof(struct LAKE_INFO_1D));
		fread(Lake, sizeof(struct LAKE_INFO_1D), nlakes+1, file);
		for (j=1; j<=nlakes; j++) {
			Lake[j].cell = calloc (Lake[j].n, sizeof(int));
			fread(Lake[j].cell, sizeof(int), Lake[j].n, file);
			Lake[j].sd = calloc (Lake[j].n_sd, sizeof(int));
			fread(Lake[j].sd, sizeof(int), Lake[j].n_sd, file);
		}
	}

	fread(&end_check,	sizeof(int),		1, 	file);
	if (end_check != 12345) {
		PRINT_ERROR("\achecking the end of resume file (%d).\n", end_check);
		exit(0);
	}
	else {
		PRINT_INFO("Check of resume file '%s' is ok.", filename);
	}

	if (switch_file_out){
		char filename[MAXLENLINE];
		FILE *file2;
		sprintf(filename, "%s.out", projectname);
		if ((file2 = fopen(filename, "a")) == NULL) {
			PRINT_ERROR("Cannot open standard output file %s.\n", filename);
		}
		else {
			PRINT_INFO("standard output redirected to %s.\n", filename);
		}
		stdout=file2;
	}

	fclose(file);
	return(1);
}



int read_file_sea_level()
{
	/*
		Reads file with sea level along time named 'projectname.SLV'
	*/

	int 	i, j;
	FILE 	*file;
	float	*aux1, *aux2, *aux3;

	sea_level = 0;

	Read_Open_Filename_Return(".SLV", "rt", "Sea level")

	n_sea_level_input_points=n_eros_level_input_points=0;
	aux1 = calloc(nmax_input_points, sizeof(float));
	aux2 = calloc(nmax_input_points, sizeof(float));
	aux3 = calloc(nmax_input_points, sizeof(float));
	for (i=0; i<nmax_input_points; i++) aux3[i]=NO_DATA;
	
	for (;;) {
		/*TAKE_LINE_2(aux1[n_sea_level_input_points], aux2[n_sea_level_input_points]);*/
		{char auxstr[MAXLENLINE], *lin; int nfields=0; while (nfields<2) {lin=fgets(auxstr, MAXLENLINE-1, file); if (lin==NULL) break; nfields=sscanf(lin, "%f %f %f", &aux1[n_sea_level_input_points], &aux2[n_sea_level_input_points], &aux3[n_sea_level_input_points]);}; if (lin==NULL) break;}
		n_sea_level_input_points++;
		if (aux3[n_sea_level_input_points-1]!=NO_DATA) 
			n_eros_level_input_points++;
		if (n_sea_level_input_points>=nmax_input_points-1 ) {
			PRINT_ERROR("Too many points (%d) in sea level file.", n_sea_level_input_points);
			break;
		}
	}
	fclose(file); 
	if (verbose_level>=1) fprintf(stdout, "\nSea level variations file contains %d points.", n_sea_level_input_points);
	var_sea_level = calloc(n_sea_level_input_points, sizeof(float *));
	var_eros_level = calloc(n_eros_level_input_points, sizeof(float *));
	for (i=0, j=0; i<n_sea_level_input_points; i++) {
		var_sea_level[i] = calloc(2, sizeof(float));
		var_sea_level[i][0] = aux1[i]*Matosec;
		var_sea_level[i][1] = aux2[i];
		if (aux3[i]!=NO_DATA) {
			var_eros_level[j] = calloc(2, sizeof(float));
			var_eros_level[j][0] = aux1[i]*Matosec;
			var_eros_level[j][1] = aux3[i];
			j++;
		}
	}
	free(aux1); free(aux2); free(aux3);
	return(1);
}



int read_file_Te()
{
	/*
		Reads elastic thickness file named 'projectname.EET'
		Linearly interpolates between given points.
	*/

	int 	i;
	FILE 	*file;
	float 	Dref=ET2RIG(Te_default) ;
	char 	filename[MAXLENFILE];

	sprintf(filename, "%s.EET", projectname);
	if ((file = fopen(filename, "rt")) == NULL) {
		PRINT_WARNING("Cannot read EET file '%s'. Constant rigidity = %3.4e", filename, Dref);
		if (isost_model>=3 && Te_default==-1) {
			/*Read the eet file as initial value for Te in Rheoiter*/
			sprintf(filename, "%s.eeth", projectname);
			if ((file = fopen(filename, "rt")) == NULL) {
				PRINT_WARNING("EET file '%s' not found. Constant rigidity = %3.4e", filename, Dref);
			}
			else {
				/*Note that x in this file is in km!*/
				if (verbose_level>=1) fprintf(stdout, "\nInitial elastic thickness read from '%s'.", filename);
				readinterplin(file, Te, Nx, x0/1000, xf/1000);
				fclose(file);
				for (i=0;i<Nx;i++) {D[i]=ET2RIG(Te[i]);}
				return(1);
			}
		}
		for (i=0;i<Nx;i++) {D[i]=Dref; Te[i]=Te_default; }
		return(0);
	}
	if (verbose_level>=1) fprintf(stdout, "\nElastic thickness at '%s'.", filename);
	readinterplin(file, Te, Nx, x0, xf);
	for (i=0;i<Nx;i++) {D[i]=ET2RIG(Te[i]);}
	fclose(file);
	return(1);
}



int read_file_Temperature()
{
	/*READS TEMPERATURE FILE 'projectname.TMP'*/

	int 	i, n_fld, ix, iz, i_x_temp, i_zinp, 
		nz_temp_input, n_max_input_temp=100;
	FILE 	*file;
	BOOL	last_line=NO;
	char 	filename[MAXLENFILE], *lin, 
		linea[MAXLENLINE];
	float	a, b, x_first, x_last, 
		*temp, *z_temp, x, z,
		**temp_input,		/*Temperature at input points of .TMP file [ºC]*/
		*x_temp_input; 		/*positions of the given geotherms of .TMP file*/

	if (isost_model<3 || switch_YSE_file) return(0);

	Temperature = alloc_matrix(Nx, Nz);

	sprintf(filename, "%s.TMP", projectname) ;
	if ((file = fopen(filename, "rt")) == NULL) {
		PRINT_WARNING("Temperature file '%s' not found. Plasticity unswitched.", filename);
		isost_model = 1 ;
		return(0);
	}
	else if (verbose_level>=1) fprintf(stdout, "\nTemperature at '%s'.", filename);

	x_temp_input =	(float *) calloc(n_max_input_temp, sizeof(float));
	temp_input =	(float **)calloc(n_max_input_temp, sizeof(float *));
	z_temp = (float *) calloc(500, sizeof(float));
	temp =   (float *) calloc(500, sizeof(float));
	nz_temp_input = 0; 
	nx_temp_input = 0; 
	do {
		lin=fgets(linea, MAXLENLINE-1, file);
		if (!lin) last_line=YES;
		n_fld = sscanf(linea, "%f %f", &a, &b);
		if ((n_fld==2 && nx_temp_input==0)) { 
			nx_temp_input=1;
			x_temp_input[0]=xmin;
			temp_input[nx_temp_input-1] =   (float *) calloc(Nz, sizeof(float));
			if (verbose_level>=1) fprintf(stdout, " Temperature laterally constant.");
		}
		if ((n_fld==1 && nx_temp_input>0) || !lin) {
			/*Interpolates vertically the geotherm (temperature-z)*/
			for (iz=0; iz<Nz; iz++) {
				z = iz*dz;
				for (i_zinp=0; i_zinp<nz_temp_input-1; i_zinp++)
					if (z<=z_temp[i_zinp+1]) break;
				temp_input[nx_temp_input-1][iz] = temp[i_zinp] + (z-z_temp[i_zinp]) * 
					(temp[i_zinp+1]-temp[i_zinp]) / (z_temp[i_zinp+1]-z_temp[i_zinp]);
				if (z <= z_temp[0]) 	temp_input[nx_temp_input-1][iz]=temp[0]; 
				if (z >= z_temp[nz_temp_input-1]) temp_input[nx_temp_input-1][iz]=temp[nz_temp_input-1]; 
			}
		}
		if (n_fld==1) { 
			nx_temp_input++; 
			nz_temp_input = 0; 
			x_temp_input[nx_temp_input-1] = a;
			temp_input[nx_temp_input-1] =   (float *) calloc(Nz, sizeof(float));
		}
		if (n_fld==2) { 
			nz_temp_input++; 
			z_temp[nz_temp_input-1] = a; 
			temp[nz_temp_input-1]   = b;
		}
	} while(!last_line);
	fclose(file);
	
	x_first = x_temp_input[0];
	x_last = x_temp_input[nx_temp_input-1];
	if (verbose_level>=1 && nx_temp_input>1) fprintf(stdout, " %d geotherms read from x=%.1f to %.1f km", nx_temp_input, x_first/1000, x_last/1000);

	/*Interpolates horizontally the geotherm (temperature-z)*/
	for (ix=0; ix<Nx; ix++) {
		x=ix*dx+x0;
		for (i_x_temp=0; i_x_temp<nx_temp_input-1; i_x_temp++)
			 if (x<=x_temp_input[i_x_temp+1]) break;
		for (iz=0; iz<Nz; iz++) {
			if (nx_temp_input>1) {
			    if (x>x_temp_input[0] && x<x_temp_input[nx_temp_input-1]) {
				Temperature[ix][iz] = temp_input[i_x_temp][iz] +
					(x-x_temp_input[i_x_temp]) *
					(temp_input[i_x_temp+1][iz]-temp_input[i_x_temp][iz]) /
					(x_temp_input[i_x_temp+1]-x_temp_input[i_x_temp]);
			    }
			    if (x<=x_temp_input[0]) 		
				Temperature[ix][iz] = temp_input[0][iz];
			    if (x>=x_temp_input[nx_temp_input-1]) 	
				Temperature[ix][iz] = temp_input[nx_temp_input-1][iz];
			}
			else 	Temperature[ix][iz]=temp_input[0][iz];
		}
	}
	
	write_file_Temperature_initial();
	
	free(z_temp);
	free(temp);
	free(x_temp_input);
	free(temp_input);
	return(1);
}





/******************************************************************************/






int write_file_erosed (float *total_erosion)
{
	int 	i, ix_min, ix_max;
	float	x ;
	FILE 	*file;

	/*PRINTS A FILE with surface transport results*/

	Write_Open_Filename_Return (".eros", "wt", !erosed_model);

	fprintf(file, "#Time: %.2fMy\n#x[km]\terosion[m]\teros_rate[m/My]\ttopo[m]", Time/Matosec);
		if (erosed_model>=2) 
			fprintf(file, "\tdischarge[m3/s]\tsedload[kg/s]\ttype\tprecipt[l/m2/yr]\tevap[l/m2/yr]");
	ix_min = MAX_2((xmin-x0-.1*dx)/dx, 0) ;	ix_max = MIN_2(floor((xmax-x0+.1*dx)/dx) + 2, Nx);
	for(i=ix_min; i<ix_max; i++) {
		float drainage_aux;
		x=x0+dx*i;
		fprintf(file, "\n%7.2f\t%8.1f\t%8.2f\t%7.1f", 
			x/1000, total_erosion[i] / dx/riverbasinwidth/denscrust, 
			eros_now[i]/(dt/Matosec) / dx/riverbasinwidth/denscrust, topo[i] );
		if (erosed_model>=2) {
			drainage_aux = drainage[i].discharge;
			if (drainage[i].lake) if (Lake[drainage[i].lake].n_sd) 
				drainage_aux = drainage[Lake[drainage[i].lake].sd[0]].discharge;
			fprintf(file, "\t%8.6f\t%8.6f\t%c\t%8.1f\t%8.1f", 
				drainage_aux, 
				drainage[i].masstr, drainage[i].type, 
				precipitation[i]*secsperyr*1e3, evaporation[i]*secsperyr*1e3);
		}
	}
	fprintf(file, "\n");
	fclose(file);
	return(1);
}


int write_file_grav_anom (float *gravanom, float *geoidanom)
{
	int 	i, ix_min, ix_max;
	float	x ;
	FILE 	*file;

	/*PRINTS A FILE WITH GRAVITY ANOMALY IN mGal*/

	Write_Open_Filename_Return (".xg", "wt", !grav_anom_type);

	fprintf(file, "#Time: %.2fMy\n#x[km]\tdg[mGal]  N[m]", Time/Matosec);
	ix_min = MAX_2((xmin-x0-.1*dx)/dx, 0) ;	ix_max = MIN_2(floor((xmax-x0+.1*dx)/dx) + 2, Nx);
	for(i=ix_min; i<ix_max; i++) {
		x=x0+dx*i; 
		fprintf(file, "\n%6.2f\t%6.1f\t%6.2f", 
			x/1000, (gravanom[i]-gravanom[ix_max-2]) * 1e5, geoidanom[i]-geoidanom[ix_max-2]);
	}
	fprintf(file, "\n");
	fclose(file);
	return(1);
}


int write_file_maxmompoint ()
{
	int	i ;
	FILE 	*file ;

	Write_Open_Filename_Return (".ysen", "wt", isost_model<3);

	fprintf(file, 	"#MMP: \tx=%.2f km\n"
			"#z(m)  Compress(MPa) Extens(MPa)  Temp(C)  Stress [MPa]\n",
			(imomentmax*dx+x0)/1e3);
	for (i=0; i<Nz; i++) {
		fprintf(file, "%.2f\t%.2f\t%.2f\t%.2f\t%.2f\n", 
			i*dz/1e3, 
			yieldcompres[i]/1e6, yieldextens[i]/1e6, 
			(switch_YSE_file)? NO_DATA:Temperature[imomentmax][i], 
			stress[imomentmax][i]/1e6);
	}
	fclose(file);
	return(1);
}



int write_file_stress ()
{
	int	ix, iz, ix_min, ix_max ;
	FILE 	*file ;

	Write_Open_Filename_Return (".strs", "wt", isost_model<3);

	fprintf(file, "#2D Stress grid (x-z) distribution.\n") ;
	fprintf(file, "#x(km)\tz(km)\tstress(MPa)\ttemper.(C)\n");
	/*Write stresses in the grid file*/
	ix_min = MAX_2((xmin-x0-.1*dx)/dx, 0) ;	ix_max = MIN_2(floor((xmax-x0+.1*dx)/dx) + 2, Nx);
	for (ix=ix_min; ix<ix_max; ix++) for (iz=0; iz<Nz; iz++) {
		fprintf (file, "%.2f\t%.2f\t%.2f\t%.2f\n", 
			(ix*dx+x0)/1e3, (iz*dz+w[ix]-zini)/1e3, stress[ix][iz]/1e6, (switch_YSE_file)? NO_DATA:Temperature[ix][iz]);
	}
	fclose(file);
	return (1);
}



int write_file_Te () 
{
	int	i;
	float	x;
	FILE 	*file ;

	Write_Open_Filename_Return (".eeth", "wt", isost_model<3);

	/*Write a file with the calculated EET*/
	fprintf(file, "# %s: Calculated EET.\n", projectname) ;
	for (i=0; i<Nx; i++) {
		x=x0+dx*i;
		if (x>xmin-dx && x<xmax+dx) {
 			fprintf(file, "%.2f\t%.2f\n", (i*dx+x0)/1e3, Te[i]);
		}
	}
	fclose(file);
	return(1);
}



int write_file_Temperature ()
{
	int	ix, iz, ix_min, ix_max ;
	FILE 	*file ;

	Write_Open_Filename_Return (".temp", "wt", isost_model<3 || switch_YSE_file);

	fprintf(file, "# 2D Temperature grid (x-z) distribution.\n") ;
	fprintf(file, "# x (km)\t z (km)\t temperature (ºC)\n") ;
	/*Write temperatures in the grid file*/
	ix_min = MAX_2((xmin-x0-.1*dx)/dx, 0) ;	ix_max = MIN_2(floor((xmax-x0+.1*dx)/dx) + 2, Nx);
	for (ix=ix_min; ix<ix_max; ix++) for (iz=0; iz<Nz; iz++) {
		fprintf (file, "%.2f\t%.2f\t%.2f\n", 
			(ix*dx+x0)/1e3, (topo[ix]-iz*dz)/1e3, Temperature[ix][iz]);
	}
	fclose(file);
	return (1);
}



int write_file_Temperature_initial ()
{
	int	ix, iz, ix_min, ix_max ;
	FILE 	*file ;

	Write_Open_Filename_Return (".tempini", "wt", isost_model<3 || verbose_level<1);

	calculate_topo(topo);
	ix_min = MAX_2((xmin-x0-.1*dx)/dx, 0) ;	ix_max = MIN_2(floor((xmax-x0+.1*dx)/dx) + 2, Nx);
	fprintf(file, "# 2D Temperature grid (x-z) distribution.\n") ;
	fprintf(file, "# x (km)\t z (km)\t temperature (ºC)\n") ;
	/*Write temperatures in the grid file*/
	for (ix=ix_min; ix<ix_max; ix++) for (iz=0; iz<Nz; iz++) {
		fprintf (file, "%.2f\t%.2f\t%.2f\n", 
			(ix*dx+x0)/1e3, (topo[ix]-iz*dz)/1e3, Temperature[ix][iz]);
	}
	fclose(file);
	return (1);
}



int write_file_time (float *w, float *topo)
{
	/*
	  WRITES deflection and elevation along time file
	*/

	int 	i;
	FILE 	*file;
	char 	filename[MAXLENLINE], filename1[MAXLENLINE], filename2[MAXLENLINE], 
		command[MAXLENLINE];
	BOOL	return_cond;
	float	youngest_age=-1e16;

	for (i=0; i<numBlocks; i++) youngest_age = MAX_2(Blocks[i].age, youngest_age);
	return_cond = (!switch_write_file && (Timefinal-Time) >= dt) 
		|| !isost_model 
		|| (((Time-last_time_file_time) < dt_record*.9999 || (!dt_record && youngest_age!=Time)) && (Timefinal-Time) >= dt) 
		|| (Time-last_time_file_time) == 0;

	calculate_topo(topo);
	if (!switch_write_file) nwrotenfiles=0;
	if (nwrotenfiles==0) {
		Write_Open_Filename_Return (".xzt", "wt", return_cond);
		fprintf(file, "#Time: \t%.3f My\n#x(km)\tw(m)\th(m)", Time/Matosec);
	}
	else {
		if (return_cond) return (0);
		sprintf(filename, "%s.xzt", projectname);
		if (verbose_level>=3) fprintf(stdout, "\nInfo: Writing file '%s' (%d times).", filename, nwrotenfiles+1);
		sprintf(filename1, "%s.aux1.xzt.tao.tmp", projectname);
		sprintf(filename2, "%s.aux2.xzt.tao.tmp", projectname);
		if ((file = fopen(filename1, "wt")) == NULL) {
			PRINT_ERROR("Cannot open auxiliar output file %s.\n", filename1);
			return (0);
		}
		fprintf(file,      "%.3f My\nw(m)\th(m)", Time/Matosec);
	}

	for(i=0; i<Nx; i++) {
		float x;
		x=x0+dx*i; 
		if (x > xmin-dx && x < (xmax+dx)) {
			if (nwrotenfiles==0)	fprintf(file, "\n%1.2f\t%1.3f\t%1.3f",
							x/1000, w[i], topo[i] );
			else			fprintf(file, "\n%1.3f\t%1.3f", 
							w[i], topo[i] );
		}
	}
	fprintf(file, "\n");
	fclose(file);
	if (nwrotenfiles > 0) {
		/*Paste the new columns into the .xzt file.*/
		sprintf(command, 
			"paste %s %s > %s", filename, filename1, filename2);
		system(command);
		rename(filename2, filename); remove(filename1);
	}

	switch_write_file_Blocks=YES;
	nwrotenfiles++;
	last_time_file_time = Time;

	return(1);
}


int write_file_Blocks()
{
	FILE 	*file ;

	/*PRINTS A FILE WITH HORIZON ALTITUDES IN COLUMNS*/

	Write_Open_Filename_Return (".pfl", "wt", !switch_write_file_Blocks);

	fprintf(file, "#tAo output file of project '%s'. t= %.2f My", projectname, Time/Matosec);
	fprintf(file, "\n#Densities:\t%8.0f", denscrust);
	for (int i=0; i<numBlocks; i++) fprintf(file, "\t%8.0f", Blocks[i].density);
	if (erosed_model>=2) fprintf(file, "\t%8.0f", denswater);
	fprintf(file, "\n#x(km),Ages->\t%8.2f", Timeini/Matosec);
	for (int i=0; i<numBlocks; i++) fprintf(file, "\t%8.2f", Blocks[i].age/Matosec);
	if (erosed_model>=2) fprintf(file, "\t   water");
	for (int i=0; i<Nx; i++) {
		float x;
		x=x0+dx*i; 
		if (x > xmin-dx && x < xmax+dx) {
			float thickness_above=0, top_block;
			fprintf(file, "\n%8.2f", x/1000);
			for (int i_Block=0; i_Block<numBlocks; i_Block++) 
				thickness_above += Blocks[i_Block].thick[i];
			top_block = Blocks_base[i]-w[i];
			fprintf(file, "\t%8.1f",  top_block);
			for (int i_Block=0; i_Block<numBlocks; i_Block++) {
				thickness_above -= Blocks[i_Block].thick[i];
				top_block += Blocks[i_Block].thick[i];
				if (Blocks[i_Block].density==denssedim) top_block -= compaction(sed_porosity, compact_depth, thickness_above, thickness_above+Blocks[i_Block].thick[i]);
				fprintf(file, "\t%8.1f",  top_block);
			}
			if (erosed_model>=2) {
				float top_water=top_block;
				if (drainage[i].lake) top_water = Lake[drainage[i].lake].alt;
				fprintf(file, "\t%8.1f", top_water);
			}
		}
	}
	fclose(file);
	return(1);
}



int write_file_resume()
{
	int 	i, j, end_check=12345;
	FILE 	*file ;

	/*
	  PRINTS A BINARY FILE WITH ALL THE INFORMATION 
	  REQUIRED TO RESTART THE PROGRAM
	*/

	Write_Open_Filename_Return (".all", "wt", !switch_write_file_Blocks);

	/*Defined in universal.h:*/
	fwrite(&Nx, 		sizeof(int),		1, 	file);
	fwrite(&Nz, 		sizeof(int),		1, 	file);
	fwrite(&verbose_level, 	sizeof(int),		1, 	file);

	fwrite(&dx, 		sizeof(float),		1, 	file);
	fwrite(&dz, 		sizeof(float),		1, 	file);

	fwrite(version,		sizeof(char),		LENGTHVERS, 	file);
	fwrite(version_input,	sizeof(char),		LENGTHVERS, 	file);

	fwrite(&switch_geograph_coor, sizeof(BOOL),	1, 	file);
	fwrite(&switch_ps, 	sizeof(BOOL),		1, 	file);
	fwrite(&switch_write_file, sizeof(BOOL),		1, 	file);


	/*Defined in geomodel.h:*/
	fwrite(&grav_anom_type, 	sizeof(int),		1, 	file);
	fwrite(&isost_model, 	sizeof(int),		1, 	file);

	fwrite(&Te_default, 	sizeof(float),		1, 	file);
	fwrite(&crust_thick_default, sizeof(float),	1, 	file);
	fwrite(&upper_crust_thick_default, sizeof(float),1, 	file);
	fwrite(&densasthen, 	sizeof(float),		1, 	file);
	fwrite(&densmantle, 	sizeof(float),		1, 	file);
	fwrite(&denssedim, 	sizeof(float),		1, 	file);
	fwrite(&denscrust, 	sizeof(float),		1, 	file);
	fwrite(&densenv, 	sizeof(float),		1, 	file);
	fwrite(&densinfill, 	sizeof(float),		1, 	file);
	fwrite(&sea_level, 	sizeof(float),		1, 	file);
	fwrite(&temp_sea_level, 	sizeof(float),		1, 	file);
	fwrite(&Time, 		sizeof(float),		1, 	file);
	fwrite(&Timefinal, 	sizeof(float),		1, 	file);
	fwrite(&Timeini, 	sizeof(float),		1, 	file);
	fwrite(&dt, 		sizeof(float),		1, 	file);
	fwrite(&dt_eros, 		sizeof(float),		1, 	file);
	fwrite(&tau, 		sizeof(float),		1, 	file);

	fwrite(projectname, 	sizeof(char),	MAXLENFILE, 	file);
	fwrite(gif_geom, 	sizeof(char),	MAXLENLINE, 	file);

	fwrite(&water_load, 	sizeof(BOOL),		1, 	file);


	/*Defined in tao+tisc.h:*/
	fwrite(&boundary_conds, 	sizeof(int),		1, 	file);
	fwrite(&nloads, 		sizeof(int),		1, 	file);
	fwrite(&nmax_input_points, sizeof(int),		1, 	file);
	fwrite(&n_sea_level_input_points, sizeof(int),	1, 	file);
	fwrite(&n_eros_level_input_points, sizeof(int),	1, 	file);
	fwrite(&n_record_times, 	sizeof(int),		1, 	file);
	fwrite(&i_first_Block_load, sizeof(int),	1, 	file);
	fwrite(&i_Block_insert, sizeof(int),		1, 	file);
	fwrite(&numBlocks, 	sizeof(int),		1, 	file);
	fwrite(&nwrotenfiles, 	sizeof(int),		1, 	file);
	fwrite(&run_type, 	sizeof(int),		1, 	file);

	fwrite(&zini, 		sizeof(float),		1, 	file);
	fwrite(&dt_record, 	sizeof(float),		1, 	file);
	fwrite(&sed_porosity, 	sizeof(float),		1, 	file);
	fwrite(&compact_depth, 	sizeof(float),		1, 	file);
	fwrite(&Kerosdif, 	sizeof(float),		1, 	file);
	fwrite(&last_time_file_time, 	sizeof(float),		1, 	file);
	fwrite(&random_topo, 	sizeof(float),		1, 	file);

	fwrite(&switch_dt_output, 	sizeof(BOOL),		1, 	file);
	fwrite(&switch_file_out, 	sizeof(BOOL),		1, 	file);
	fwrite(&switch_gradual, 	sizeof(BOOL),		1, 	file);
	fwrite(&switch_insert_load, 	sizeof(BOOL),		1, 	file);
	fwrite(&switch_topoest, 		sizeof(BOOL),		1, 	file);
	fwrite(&switch_write_file_Blocks, sizeof(BOOL),		1, 	file);
	fwrite(&deform_sed, sizeof(BOOL),		1, 	file);

	/*Defined in tao.h:*/
	fwrite(&imomentmax, 	sizeof(int),		1, 	file);
	fwrite(&nx_temp_input, 	sizeof(int),		1, 	file);
	fwrite(&nbasins, 	sizeof(int),		1, 	file);
	fwrite(&nlakes,	 	sizeof(int),		1, 	file);
	fwrite(&n_image, 	sizeof(int),		1, 	file);
	fwrite(&hydro_model,	sizeof(int),		1, 	file);
	fwrite(&erosed_model,	sizeof(int),		1, 	file);
	fwrite(&eros_bound_cond,	sizeof(char),		2, 	file);

	fwrite(&x0, 		sizeof(float),		1, 	file);
	fwrite(&xf, 		sizeof(float),		1, 	file);
	fwrite(&xmin, 		sizeof(float),		1, 	file);
	fwrite(&xmax, 		sizeof(float),		1, 	file);
	fwrite(&zmin, 		sizeof(float),		1, 	file);
	fwrite(&zmax, 		sizeof(float),		1, 	file);
	fwrite(&horz_force, 	sizeof(float),		1, 	file);
	fwrite(&vert_force, 	sizeof(float),		1, 	file);
	fwrite(&appmoment, 	sizeof(float),		1, 	file);
	fwrite(&Keroseol, 	sizeof(float),		1, 	file);
	fwrite(&Ksedim, 		sizeof(float),		1, 	file);
	fwrite(&critical_slope, 	sizeof(float),		1, 	file);
	fwrite(&K_river_cap, 	sizeof(float),		1, 	file);
	fwrite(&erodibility, 	sizeof(float),		1, 	file);
	fwrite(&erodibility_sed, sizeof(float),		1, 	file);
	fwrite(&l_fluv_sedim, 	sizeof(float),		1, 	file);
	fwrite(&lost_rate, 	sizeof(float),		1, 	file);
	fwrite(&evaporation_ct, 	sizeof(float),		1, 	file);
	fwrite(&riverbasinwidth, 	sizeof(float),		1, 	file);
	fwrite(&rain, 		sizeof(float),		1, 	file);
	fwrite(&Krain,	 	sizeof(float),		1, 	file);
	fwrite(&CXrain, 		sizeof(float),		1, 	file);

	fwrite(&switch_strs_history, 	sizeof(BOOL),		1, 	file);
	fwrite(&switch_YSE_file, 	sizeof(BOOL),		1, 	file);


	/*Arrays:*/
	fwrite(w, 		sizeof(float),	Nx, 	file);
	fwrite(D, 		sizeof(float),	Nx, 	file);
	fwrite(q, 		sizeof(float),	Nx, 	file);
	fwrite(Dw, 		sizeof(float),	Nx, 	file);
	fwrite(Dq, 		sizeof(float),	Nx, 	file);
	fwrite(h_water, 		sizeof(float),	Nx, 	file);
	fwrite(h_last_unit, 	sizeof(float),	Nx, 	file);
	fwrite(Te, 		sizeof(float),	Nx, 	file);
	fwrite(crust_thick, 	sizeof(float),	Nx, 	file);
	fwrite(upper_crust_thick,sizeof(float),	Nx, 	file);
	fwrite(topo, 		sizeof(float),	Nx, 	file);
	fwrite(Blocks_base, 	sizeof(float),	Nx, 	file);
	fwrite(yieldcompres, 	sizeof(float),	Nz, 	file);
	fwrite(yieldextens, 	sizeof(float),	Nz, 	file);
	fwrite(horiz_record_time, sizeof(float), n_record_times, file);

	for (i=0; i<n_sea_level_input_points; i++) 
		fwrite(var_sea_level[i], sizeof(float), 2, file);
	for (i=0; i<n_eros_level_input_points; i++) 
		fwrite(var_eros_level[i], sizeof(float), 2, file);

	if (isost_model>=3 && !switch_YSE_file) {
		for (i=0; i<Nx; i++) 
			fwrite(Temperature[i], sizeof(float), Nz, file);
	}
	if (isost_model>=3) {
		for (i=0; i<Nx; i++) 
			fwrite(stress[i], sizeof(float), Nz, file);
	}

	fwrite(Blocks, 		sizeof(struct BLOCK_1D),	numBlocks, file);
	for (i=0; i<numBlocks; i++) {
		fwrite(Blocks[i].thick, 		sizeof(float),	Nx, 	file);
	}
	for (i=0; i<numBlocks; i++) {
	    if (Blocks[i].type == 'S') {
		fwrite(Blocks[i].detr_ratio, 	sizeof(float),	Nx, 	file);
		fwrite(Blocks[i].detr_grsize, 	sizeof(float),	Nx, 	file);
	    }
	}

	if (erosed_model) {
		fwrite(eros_now, 	sizeof(float),	Nx, 	file);
		fwrite(total_erosion, 	sizeof(float),	Nx, 	file);
	}
	if (hydro_model) {
		fwrite(precipitation, 	sizeof(float),	Nx, 	file);
		fwrite(evaporation, 	sizeof(float),	Nx, 	file);
		fwrite(Lake, sizeof(struct LAKE_INFO_1D), nlakes+1, file);
		for (j=1; j<=nlakes; j++) {
			fwrite(Lake[j].cell, sizeof(int), Lake[j].n, file);
			fwrite(Lake[j].sd, sizeof(int), Lake[j].n_sd, file);
		}
	}

	fwrite(&end_check,	sizeof(int),		1, 	file);

	fclose(file);
	return(1);
}

