/*
	GENERAL FUNCTION LIBRARY FOR tao.c
	More info in tao/doc/ 
*/




int Allocate_Memory()
{
	/* Allocates dynamic memory for the arrays and initializes them to zero*/

	w = 		(float *) calloc(Nx, sizeof(float));
	D = 		(float *) calloc(Nx, sizeof(float));
	q = 		(float *) calloc(Nx, sizeof(float));
	Dw = 		(float *) calloc(Nx, sizeof(float));
	Dq = 		(float *) calloc(Nx, sizeof(float));
	h_water = 	(float *) calloc(Nx, sizeof(float));
	h_last_unit =	(float *) calloc(Nx, sizeof(float));
	Te =	 	(float *) calloc(Nx, sizeof(float));
	crust_thick = 	(float *) calloc(Nx, sizeof(float));
	upper_crust_thick = (float *) calloc(Nx, sizeof(float));
	topo =  	(float *) calloc(Nx, sizeof(float));
	Blocks_base = 	(float *) calloc(Nx, sizeof(float));
	yieldcompres =	(float *) calloc(Nz, sizeof(float));
	yieldextens = 	(float *) calloc(Nz, sizeof(float));

	Blocks =		(struct BLOCK_1D *) calloc(NmaxBlocks, sizeof(struct BLOCK_1D));

	if (hydro_model) {
		int i; 
		precipitation =	(float *) calloc(Nx, sizeof(float));
		evaporation =	(float *) calloc(Nx, sizeof(float));
		sortcell = (int *) calloc(Nx, sizeof(int));
		for (i=0; i<Nx; i++) {sortcell[i]=i;}
		drainage = (struct DRAINAGE_1D *) calloc(Nx, sizeof(struct DRAINAGE));
	}
	if (erosed_model) {
		eros_now =	(float *) calloc(Nx, sizeof(float));
		total_erosion =	(float *) calloc(Nx, sizeof(float));
	}

	/*lake 0 is not used; thus create it now with no content*/
	Lake = (struct LAKE_INFO_1D *) calloc(1, sizeof(struct LAKE_INFO_1D));

	return(1);
}





float calculate_topo(float *topo_new)
{
	float mean=0;
	/*Calculates current topography based on Blocks, Blocks_base and deflection*/
	PRINT_DEBUG("Entering");
	for (int i=0; i<Nx; i++) {
		float thickness_above=0;
		for (int i_Block=0; i_Block<numBlocks; i_Block++) 
			thickness_above += Blocks[i_Block].thick[i];
		topo_new[i] = Blocks_base[i]-w[i];
		for (int i_Block=0; i_Block<numBlocks; i_Block++) {
			thickness_above -= Blocks[i_Block].thick[i];
			topo_new[i] += Blocks[i_Block].thick[i];
			if (Blocks[i_Block].density==denssedim) topo_new[i] -= compaction(sed_porosity, compact_depth, thickness_above, thickness_above+Blocks[i_Block].thick[i]);
		}
		mean += topo_new[i];
	}
	mean /= Nx;
	return (mean);
}



int Delete_Block(int i_Block)
{
	int  	k;

	/*Deallocates one Block*/

	PRINT_DEBUG("Block being deleted: %d ; numBlocks= %d ; i_first_Block_load = %d ; i_Block_insert = %d", i_Block, numBlocks, i_first_Block_load, i_Block_insert);
	free(Blocks[i_Block].thick);
	if (Blocks[i_Block].type == 'S') {
		free(Blocks[i_Block].detr_ratio);
		free(Blocks[i_Block].detr_grsize);
	}
	for (k=i_Block; k<numBlocks-1; k++) Blocks[k] = Blocks[k+1];
	numBlocks--;
	return(1);
}



int gradual_Block()
{
	float 	Dhl;

	/*Non-instantaneous loading of a file load (distributed along time).*/

	/*interpolation can only last until next load because then the original load shape is then lost*/
	if (!switch_gradual || Time>Blocks[i_Block_insert].time_stop-dt*.1) return(0);

	for (int i=0; i<Nx; i++) { 
		Dhl = h_last_unit[i]*dt/(Blocks[i_Block_insert].time_stop-Blocks[i_Block_insert].age);
		/*Increments the load for this time interval*/
		if (Blocks[i_Block_insert].type != 'H') {
			if (Dhl>=0) {
				Blocks[i_Block_insert].thick[i] += Dhl;
			}
			else {
				float 	h_load_aux, h_load_aux2;
				int	k;
				h_load_aux = fabs((double) Dhl);
				for (k=i_Block_insert-1; h_load_aux>0 && k>=0; k--) {
					h_load_aux2 = MIN_2(Blocks[k].thick[i], h_load_aux);
					h_load_aux -= h_load_aux2;
					Blocks[k].thick[i] -= h_load_aux2;
				}
				/*k is the deepest eroded Block*/
				if (k==-1) {
					Blocks_base[i] -= h_load_aux;
				}
			}
		}
		Dq[i] += (Blocks[i_Block_insert].density-densenv)*g*Dhl;
	}

	fflush(stdout);
	return(1);
}




int Init_Stress()
{
	int 	ix, iz, iter, numiter=50;
	float 	iterforce, criterio, mechanical_thickness, refstress;


	if (isost_model<3) return (0);
	
	stress = alloc_matrix(Nx, Nz);

	/*Distributes the tectonic force along the strong lithosphere.*/
	if (switch_strs_history && horz_force) {
		for (ix=0; ix<Nx; ix++) {
			/*Finds the YSE from temperature & geometry*/
			if (!switch_YSE_file) yield_stress_envelope (
				Temperature[ix], Nz, dz, 0, 
				upper_crust_thick[ix], crust_thick[ix],
				isost_model,
				yieldcompres, yieldextens,
				&mechanical_thickness
			);

			for (iter=0, refstress=-horz_force/Nz/dz ;  iter<numiter ; iter++) {
				for (iz=iterforce=0; iz<Nz; iz++) {
					stress[ix][iz] = ((refstress>0) ? MIN_2(refstress, yieldextens[iz]) : MAX_2(refstress, yieldcompres[iz]) ) ;
					iterforce += stress[ix][iz] * dz ;
				}
				/*horz_force tiene el criterio de signos contrario a stress.*/
				criterio = (-horz_force - iterforce) ;
				if ( fabs(criterio) < 1e9 ) break ;
				refstress += criterio/fabs(criterio) * 1000e6 / pow(2,iter);
			}
			if (iter>=numiter && refstress<1000e6) {fprintf(stderr, "ERROR: horz_force is bigger than lithospheric strength!"); exit(0);}
			if (verbose_level>=4) fprintf (stderr, "\nRefstress=%.2e\tIterforce:%.2e\tCriterio:%.2e\tIters:%d", refstress, iterforce, criterio, iter);
		}
	}

	return(1);
}



int insert_new_Block(int num_new_Block)
{
	struct BLOCK_1D	Block_aux;

	/*Creates a new Block and increments numBlocks by 1.
		num_new_Block ranges from 0 to numBlocks.
		Blocks number num_new_Block and above are shifted upwards.
		If num_new_Block == numBlocks then a new block is created on top of all. 
	*/

	if (verbose_level>=2) fprintf(stdout, "  b"); fflush(stdout);
	PRINT_DEBUG("New Block being created: %d ; numBlocks= %d ; i_first_Block_load = %d ; i_Block_insert = %d", num_new_Block, numBlocks, i_first_Block_load, i_Block_insert);
	if (numBlocks>NmaxBlocks-5) PRINT_WARNING("Lots of Blocks! "); 

	Blocks[numBlocks].thick = 	(float *) calloc(Nx, sizeof(float));

	//shift blocks upwards
	Block_aux = Blocks[numBlocks];
	for (int j_Block=numBlocks; j_Block>num_new_Block; j_Block--) 
		Blocks[j_Block] = Blocks[j_Block-1];
	Blocks[num_new_Block] = Block_aux;

	/*Default properties*/
	Blocks[num_new_Block].type = '-';
	Blocks[num_new_Block].age = Time;
	Blocks[num_new_Block].density = 0;
	Blocks[num_new_Block].erodibility = erodibility;
	Blocks[num_new_Block].vel = 0;
	Blocks[num_new_Block].last_vel_time = Time;
	Blocks[num_new_Block].time_stop = 9999*Matosec;
	Blocks[num_new_Block].shift = 0;
	Blocks[num_new_Block].last_shift = 0;
	
	numBlocks++;

	return(1);
}



int match_parameter(char *str1, char *str2, int show, int replace, char *line)
{
	BOOL switch_debug=NO;
	int nparams=0;

	Match_Param_Replace_int ( "Nx", 	Nx, 	0 )
	Match_Param_Replace_int ( "Nz", 	Nz, 	0 )
	Match_Param_Replace_flt ( "x0", 	x0, 	0 )
	Match_Param_Replace_flt ( "xf", 	xf, 	0 )
	Match_Param_Replace_flt ( "xmin", 	xmin, 	0 )
	Match_Param_Replace_flt ( "xmax", 	xmax, 	0 )
	Match_Param_Replace_flt ( "zmin", 	zmin, 	0 )
	Match_Param_Replace_flt ( "zmax", 	zmax, 	0 )
	Match_Param_Replace_flt ( "Te",		Te_default, 	0 )
	Match_Param_Replace_flt ( "crust",	crust_thick_default, 	0 )
	Match_Param_Replace_flt ( "ucrust",	upper_crust_thick_default, 	0 )
	Match_Param_Replace_flt ( "zini",	zini, 	0 )
	Match_Param_Replace_flt ( "random_topo",	random_topo, 	0 )
	Match_Param_Replace_flt ( "densasthen",	densasthen, 	0 )
	Match_Param_Replace_flt ( "densmantle",	densmantle, 	0 )
	Match_Param_Replace_flt ( "denscrust",	denscrust, 	0 )
	Match_Param_Replace_flt ( "densinfill",	densinfill, 	0 )
	Match_Param_Replace_flt ( "denssedim",	denssedim, 	0 )
	Match_Param_Replace_flt ( "densenv",	densenv, 	0 )
	Match_Param_Replace_flt ( "sed_porosity",	sed_porosity, 	0 )
	Match_Param_Replace_flt ( "compact_depth",	compact_depth, 	0 )
	Match_Param_Replace_int ( "erosed_model",	erosed_model, 	0 )
	Match_Param_Replace_flt ( "Kerosdif",	Kerosdif, 	0 )
	Match_Param_Replace_flt ( "Keroseol",	Keroseol, 	0 )
	Match_Param_Replace_flt ( "Ksedim",	Ksedim, 	0 )
	Match_Param_Replace_int ( "hydro_model",	hydro_model, 	0 )
	Match_Param_Replace_flt ( "rain",	rain, 	0 )
	Match_Param_Replace_flt ( "Krain",	Krain, 	0 )
	Match_Param_Replace_flt ( "relhumid",	relative_humidity, 	0 )
	Match_Param_Replace_flt ( "CXrain",	CXrain, 	0 )
	Match_Param_Replace_flt ( "evaporation",	evaporation_ct, 	0 )
	Match_Param_Replace_flt ( "riverbasinwidth",	riverbasinwidth, 	0 )
	Match_Param_Replace_flt ( "lost_rate",	lost_rate, 	0 )
	Match_Param_Replace_flt ( "K_river_cap",	K_river_cap, 	0 )
	Match_Param_Replace_flt ( "erodibility",	erodibility, 	0 )
	Match_Param_Replace_flt ( "erodibility_sed",	erodibility_sed, 	0 )
	Match_Param_Replace_flt ( "l_fluv_sedim",	l_fluv_sedim, 	0 )
	Match_Param_Replace_flt ( "temp_sea_level",	temp_sea_level, 	0 )
	Match_Param_Replace_int ( "deform_sed",  	deform_sed, 	0 )
	Match_Param_Replace_chr ( "eros_bound_cond",	eros_bound_cond,  	0 )
	Match_Param_Replace_flt ( "Timeini",	Timeini, 	0 )
	Match_Param_Replace_flt ( "Timefinal",	Timefinal, 	0 )
	Match_Param_Replace_flt ( "tau",	tau, 	0 )
	Match_Param_Replace_flt ( "dt",		dt, 	0 )
	Match_Param_Replace_flt ( "dt_eros",	dt_eros, 	0 )
	Match_Param_Replace_flt ( "dt_record",	dt_record, 	0 )
	Match_Param_Replace_int ( "boundary_conds",	boundary_conds, 	0 )
	Match_Param_Replace_flt ( "horz_force",	horz_force, 	0 )
	Match_Param_Replace_flt ( "vert_force",	vert_force, 	0 )
	Match_Param_Replace_flt ( "moment",	appmoment, 	0 )
	Match_Param_Replace_int ( "isost_model",	isost_model, 	0 )
	Match_Param_Replace_int ( "water_load", 	water_load, 	0 )
	Match_Param_Replace_int ( "switch_topoest", switch_topoest, 	0 )
	Match_Param_Replace_int ( "grav_anom",	grav_anom_type, 	0 )
	Match_Param_Replace_int ( "switch_files",	switch_write_file, 	0 )
	Match_Param_Replace_int ( "switch_ps",  	switch_ps, 	0 )
	Match_Param_Replace_int ( "verbose_level",	verbose_level, 	0 )

	/*old version*/
	Match_Param_Replace_flt ( "erodability",	erodibility, 	1 )
	Match_Param_Replace_flt ( "erodability_sed",	erodibility_sed, 	1 )
	Match_Param_Replace_int ( "switch_verbose",	verbose_level, 	1 )
	Match_Param_Replace_int ( "switch_debug",	switch_debug, 	1 )
	if (switch_debug) verbose_level=3;
	Match_Param_Replace_flt ( "alt0",	zini, 	1 )
	Match_Param_Replace_int ( "lith_type",	isost_model, 	1 )
	Match_Param_Replace_int ( "erosed_type",	erosed_model, 	1 )
	Match_Param_Replace_int ( "switch_erosed",	erosed_model, 	1 )
	Match_Param_Replace_int ( "switch_sea", 	water_load, 	1 )
	Match_Param_Replace_flt ( "ymin", 	zmin, 	1 )
	Match_Param_Replace_flt ( "ymax", 	zmax, 	1 )
	Match_Param_Replace_flt ( "dtmemounit",	dt_record, 	1 )

	return (nparams);
}



int make_gravi_body(float *upper_hor, float *lower_hor, float *body_x, float *body_z)
{
	/*Makes a gravity body from two horizons to be used
	with subrutine anompolig() to calculate its gravitational effect
	*/
	int 	ix, np_body=0;

	for (ix=0; ix<Nx; ix++) {
		body_x[np_body] = x0+ix*dx;	body_z[np_body] = upper_hor[ix];
		np_body++ ;
	}
	body_x[np_body] = body_x[np_body-1]+1e8;	body_z[np_body] = upper_hor[Nx-1];
	np_body++ ;
	body_x[np_body] = body_x[np_body-1];		body_z[np_body] = lower_hor[Nx-1];
	np_body++ ;
	for (ix=Nx-1; ix>=0; ix--) {
		body_x[np_body] = x0+ix*dx;	body_z[np_body] = lower_hor[ix];
		np_body++ ;
	}
	body_x[np_body] = body_x[np_body-1]-1e8;	body_z[np_body] = lower_hor[0];
	np_body++ ;
	body_x[np_body] = body_x[np_body-1];		body_z[np_body] = upper_hor[0];
	np_body++ ;

	return(np_body);
}




int LES_matrix (double **A,	/* Matrix of the lineal equation system. */
		double *b, 	/* Independent term. */
		float *D, 	/* Rigidity [N·m] */
		float *q,	/* Total load (not including restitutive force).*/
		float *Dq,	/* Load increment [N] which effect (deflection) is to be calculated (not including restitutive force).*/
		float *w, 
		BOOL doing_visco)	/* Uses present deflection for certain terms. */
{
	register int 	i, j, NDi=3, NDs=3;
	double		Krest, dx2, dx3, dx4;
	char		filename[MAXLENFILE];

	/* 
	  DEFINE LA MATRIZ CASI DIAGONAL (en banda) DEL SISTEMA DE ECUACIONES 
	  RESULTANTE DE LA ECUACION DIFERENCIAL ELASTICA O VISCOELASTICA  
	  A·w = b DISCRETIZADA POR EL METODO DE LAS DIFERENCIAS FINITAS
	*/

	dx2 = dx*dx;
	dx3 = dx2*dx;
	dx4 = dx2*dx2;
	for (i=0; i<Nx; i++) {
		for (j=0; j<NDi+1+NDs; j++) A[i][j] = 0.;
		b[i]= 0.;
	}

	/*In viscoelastic model nu=.5*/
	if (doing_visco)  {  for (i=0; i<Nx; i++) D[i] *= (1-nu*nu)/(1-.25) ;	}	


	for (i=2; i<Nx-2; i++)
	{
		/*Krest is the bouyancy or restitution force constant*/
		GET_KREST(Krest, Dq, i)

		/*Matrix definition. A[i][3] are the elements of the diagonal.*/
		/*DON'T take out the (double) conversions!*/
		A[i][1] = (double) +1*D[i-1];
		A[i][2] = (double) -2*D[i]   -2*D[i-1]		+1*horz_force*dx2;
		A[i][3] = (double) +1*D[i+1] +4*D[i]  +1*D[i-1]	-2*horz_force*dx2 + Krest*dx4; 
		A[i][4] = (double) -2*D[i]   -2*D[i+1]		+1*horz_force*dx2;
		A[i][5] = (double) +1*D[i+1];
		
		switch (doing_visco) {
			case 0:	/* ELASTIC EQUATION */
				/*Term of horz_force is justified in Karner, 1992.*/
				b[i] = (double) dx4 * ((Dq[i] - horz_force*(w[i+1]-2*w[i]+w[i-1])/dx2)); 
				break;

			case 1:	/* VISCOELASTIC EQUATION */
				/*Viscoelastic equation, see Nadai (1963); Lambeck (1983)*/
				b[i] = (double) dx4 * ((q[i] - w[i]*Krest - horz_force*(w[i+1]-2*w[i]+w[i-1])/dx2) / tau /*termino que deberia estar:+Dq[i]/dt*/);
				break;
		}
	}

	/* BOUNDARY CONDITIONS */
	/* The boundary conditions, in the finite differences method, are not given exactly 
		in the boundary but in the next (2nd derivative) or second next point (3rd
		derivative). So, in order to minimize the effect of this, when 2nd or 3rd 
		derivatives are going to be imposed (boundary_conds !=2), high values of Nx 
		(number of discretizing points) are required.*/

	/*Boundary Conditions at LEFT.*/
	switch (boundary_conds) {
	  case 0:	/* B.C. simetric continuous half-plate (the model only considers half of a simetric plate). */
		A[0][3]	= 1 ;
		A[0][5]	= -1;
		b[0] 	= 0 ;		/* 1ª Derivada 0 en el extremo izda.*/
		A[1][2]	= (double) -1;
		A[1][3]	= (double) +3 - horz_force/D[1]*dx2;
		A[1][4]	= (double) -3 + horz_force/D[1]*dx2;	/* Derivada 3ª en el extremo igual a */
		A[1][5]	= (double) +1;				/* la applied vertical shear force.  */
		b[1] 	= (double) -vert_force/2 /D[1]*dx3 / ((doing_visco)? tau : 1) ;
		break;
	  case 1:	/* B.C. broken plate at left, free at right.*/
	  case 3:	/* B.C. free ends. */
		/* 2nd derivate at point #2 = applied moment */
		A[0][3]	= +1;
		A[0][4]	= -2;
		A[0][5]	= +1;
		b[0] 	= (Time==Timeini || doing_visco || (isost_model>=3 && !switch_strs_history)) ? -appmoment*dx2/D[1] : 0 ;
		/* 3rd derivate (point #2.5) = applied vertical force*/
		A[1][2]	= (double) -1;
		A[1][3]	= (double) +3 - horz_force/D[1]*dx2;
		A[1][4]	= (double) -3 + horz_force/D[1]*dx2;
		A[1][5]	= (double) +1;
		b[1] 	= (double) (Time==Timeini || doing_visco || (isost_model>=3 && !switch_strs_history)) ? -vert_force/D[1]*dx3 / ((doing_visco)? tau : 1) : 0 ;
		break;
	  case 2:	/* B.C. continuous plate*/
		A[0][4]	= 1 ;
		b[0] 	= 0 ;
		A[1][2]	= 1 ;
		A[1][4]	= -1;
		b[1] 	= 0 ;
		break;
	  case 4:	/* B.C. hunging plate*/
	  case 5:	/* B.C. Broken plate at right, free at left.*/
		A[0][3]	= 1 ;
		A[0][4]	= -2;		/* Derivada 2ª en el pto. 1 nula */
		A[0][5]	= 1 ;		/*	  (no moment).			 */
		b[0] 	= 0 ;
		A[1][3]	= 1 ;		/* Deflexi¢n 0 en el extremo.	*/
		b[1] 	= 0 ;
		break;
	}

	/*B.C. at RIGHT.*/
	switch (boundary_conds) {
	  case 0:	/* B.C. Simetric continuous half-plate*/
	  case 2:	/* B.C. placa continua*/
		A[Nx-2][2] = 1 ;
		A[Nx-2][4] = -1;		/* 1ª Derivada 0 en el extremo.*/
		b[Nx-2]	   = 0 ;
		A[Nx-1][2] = 1 ;		/* Flexi¢n 0 en el extremo.*/
		b[Nx-1]	   = 0 ;
		break;
	  case 1:	/* B.C. Broken at left, free at right*/
	  case 3:	/* B.C. Free ends*/
		A[Nx-1][3] = 1 ;
		A[Nx-1][2] = -2;		/* Derivada 2ª en el pto. 1 nula*/
		A[Nx-1][1] = 1 ;		/*	  (no moment).*/
		b[Nx-1]	= 0 ;
		A[Nx-2][1] = -1;
		A[Nx-2][2] = +3;
		A[Nx-2][3] = -3;		/*	 Derivada 3ª en el 2 nula*/
		A[Nx-2][4] = 1 ;		/*  (no vertical shear force).*/
		b[Nx-2]	= 0 ;
		break;
	  case 4:	/* B.C. hunging plate*/
		A[Nx-1][3] = 1 ;
		A[Nx-1][2] = -2;		/* Derivada 2ª en el pto. 1 nula*/
		A[Nx-1][1] = 1 ;		/*	  (no hay momento).*/
		b[Nx-1]	= 0 ;
		A[Nx-2][3] = 1 ;		/* Null deflection 0 at the end.*/
		b[Nx-2]	   = 0 ;
		break;
	  case 5:	/* B.C. Broken plate at right, free at left.*/
		A[Nx-1][3] = 1 ;
		A[Nx-1][2] = -2;		/* Derivada 2ª en el pto. 1 nula*/
		A[Nx-1][1] = 1 ;		/*	  (moment).*/
		b[Nx-1]	= (double) (Time==Timeini || doing_visco || (isost_model>=3 && !switch_strs_history)) ? -appmoment*dx2/D[Nx-2]/2 : 0 ;
		A[Nx-2][1] = -1;
		A[Nx-2][2] = +3;
		A[Nx-2][3] = -3;		/* Derivada 3ª en el 2 nula*/
		A[Nx-2][4] = 1 ;		/*   (shear force).*/
		b[Nx-2]	= (double) (Time==Timeini || doing_visco || (isost_model>=3 && !switch_strs_history)) ? -vert_force*dx3/D[Nx-2]/ ((doing_visco)? tau : 1) : 0 ;
		break;
	}
	/* Restores the original value of rigidity */
	if (doing_visco)  { for (i=0; i<Nx; i++) D[i] /= (1-nu*nu)/(1-.25) ;}

	sprintf(filename, "%s.mtrx", projectname); 
	remove(filename);
	if (verbose_level>=4 && switch_write_file && Nx<=1001) { 
		WriteAlmostDiagonalMatrix(A, b, Nx, filename, NDs, NDi);
	}
	return(1);
}




float moment_calculator (float 	d2wdx2, 
			float 	horz_force, 
			float 	*yieldcompres, 
			float 	*yieldextens, 
			float 	*stress, 
			int 	isost_model, 		/*3:oceanic; 4:cont always decoupled; 5:cont. never deoupled; 6:cont. auto deoupled*/
			float 	decoupl_depth, 		/*In m. Only used when isost_model==4*/
			int 	Nz, 
			float 	dz, 
			float 	*refstressdir, 
			int 	*ncapas) 		/*Number of decoupled layers*/
{
	/* RETURNS THE CALCULATED MOMENT IN THE PLATE WHITH THE YIELD STRESS 
	ENVELOPE TAKING INTO ACCOUNT DECOUPLING IN SUBLAYERS */

	int	i, iter, j, itop, ifloor, layer=0, numlayers, numiter=50, k;
	float	stress_distrib_slope,			/*Slope of the linear part.*/
			total_moment=0, momentlayer, 	/*Total and layer moments.*/
			decoupl_stress_limit=50e6, 		/*Default decoupling yield stress.*/
			yield_stress_minim=10e6, 		/*Minimum yield stress. This defines mechanical thickness. Ranalli, 1994.*/
			z, 	/*Depth.*/
			ztoplayer[10], zfloorlayer[10], /*Top & base of each layer.*/
			linearstress, z_null_strs, 
			pressurelayer, pressure, 
			abspressure, mecanthick=0, 
			criterio, 
			Dsigma, 
			backpressure, refstress=0, 
			*refstressv; 					/*Stress due to tect. force.*/
	BOOL	switch_saturatedlayer;			/*YES if that decoupled layer is entirely at the yield stress due to horz_force*/

	refstressv = (float *) calloc(Nz, sizeof(float));

	if (!d2wdx2) d2wdx2 = 1e-12; /*very little*/
	stress_distrib_slope = - d2wdx2 * E / (1-nu*nu);

	/*Distributes the tectonic force along the strong lithosphere.*/
	if (horz_force) {
		refstress = -horz_force/Nz/dz ;
		for (iter=0; iter<numiter; iter++) {
			backpressure=0 ;
			for (i=0; i<Nz; i++) {
				refstressv[i] = ((refstress>0) ? MIN_2(refstress, yieldextens[i]) : MAX_2(refstress, yieldcompres[i]) )  ;
				backpressure += refstressv[i] * dz ;
			}
			/*horz_force has opposite sign criterion than stress.*/
			criterio = (-horz_force - backpressure) ;
			if ( fabs(criterio) < 1e9 ) break ;
			refstress += criterio/fabs(criterio) * 1000e6 / pow(2,iter)  ;
		}
		if (iter>=numiter && refstress<1000e6) fprintf(stderr, "_!_ ") ;
	}


	/*Define decoupled layers.*/
	for (i=0; i<Nz; i++) {
		z = i*dz ;
		if (yieldextens[i]-yieldcompres[i] > 2*yield_stress_minim /*decoupl_stress_limit*/) {
			ztoplayer[layer] = z;
			for (; i<Nz; i++) {
				z = i*dz;
				Dsigma = yieldextens[i]-yieldcompres[i];
				if (
					(Dsigma < decoupl_stress_limit && isost_model == 6)
					|| (decoupl_depth>i*dz && decoupl_depth<=(i+1)*dz && isost_model == 4)
					|| i==Nz-1
				   ) {
					int iz;
					/*Base of layer is controlled by yield_stress_minim*/
					for (	iz=i; 
						iz>=0 && (fabs(yieldextens[iz])<yield_stress_minim 
						|| fabs(yieldcompres[iz])<yield_stress_minim); 
						iz--
					) ;
					zfloorlayer[layer] = iz*dz ;
					mecanthick += zfloorlayer[layer] - ztoplayer[layer];
					layer++;
					break;
				}
			}
		}
	}
	numlayers=layer;
	/*printf("\nd2wdx2=%.2e, horz_force=%.2e, refstrs=%.2e, Nz=%d, dz=%.2e, numlayers=%d",
		d2wdx2,horz_force,refstress,Nz,dz, numlayers);*/


	for (i=0 ; i<Nz; i++) stress[i]=0 ;
	pressure = 0 ;


	/*Distribute bending stresses along each decoupled layer.*/
	for (layer=0; layer<numlayers ; layer++) {
		/*We need to find the depth where flexural bending stress is zero, crossing from positive to negative. This depth must accomplish that the integrated force equals horz_force*/
		z_null_strs = (zfloorlayer[layer]+ztoplayer[layer])/2 ;
		itop=ztoplayer[layer]/dz ; ifloor=zfloorlayer[layer]/dz ;
		for (j=itop, switch_saturatedlayer=YES; j<=ifloor; j++) {
			if ((refstressv[j]>0 && refstressv[j]<yieldextens[j]) || (refstressv[j]<=0 && refstressv[j]>yieldcompres[j])) {
				switch_saturatedlayer=NO ;
			}
		}
		if (switch_saturatedlayer==NO) {
			for (i=1 ; i<=numiter; i++) {
				/*Iterates until null stress point position reach convergence.*/
				pressurelayer = abspressure = momentlayer = 0 ;
				for (j=itop; j<=ifloor; j++) {
					z = j*dz ;
					linearstress = refstress+ stress_distrib_slope * (z-z_null_strs) ;
					if (linearstress<yieldextens[j] && linearstress>yieldcompres[j]) 
											stress[j] = linearstress;
					else {
						if (linearstress>yieldextens[j])	stress[j] = yieldextens[j];
						else					stress[j] = yieldcompres[j];
					}
					pressurelayer += (stress[j] - refstressv[j])* dz ;
					abspressure += fabs(stress[j]- refstressv[j]) * dz ;
					momentlayer += (stress[j] - refstressv[j]) * (z-z_null_strs) * dz;
				}
				criterio = pressurelayer ;
				if (fabs(criterio)/abspressure < .001) break ;
				z_null_strs -= (criterio*d2wdx2)/fabs(criterio*d2wdx2) * dz*(Nz-1) /pow(2, i) ;
			}
			if (i>=numiter && horz_force == 0) fprintf(stdout, "**%d", layer) ;
		}
		else {
			momentlayer = 0 ;
			for (i=itop ; i<ifloor; i++) stress[i]=refstressv[i] ;
		}
		total_moment += momentlayer ; 
		pressure += pressurelayer ;
	}

	free(refstressv);
	*refstressdir = refstress ;
	*ncapas = numlayers ;
	return (total_moment) ;
}



float moment_calculator_hist (
			float 	d2wdx2, 
			float 	horz_force, 		/*not used ! */
			float 	*yieldcompres, 
			float 	*yieldextens, 
			float 	*stress, 
			int 	isost_model, 		/*3:oceanic; 4:cont always decoupled; 5:cont. never deoupled; 6:cont. auto deoupled*/
			float 	decoupl_depth, 		/*In m. Only used when isost_model==4*/
			int 	Nz, 
			float 	dz, 
			float 	*totalmoment,  		/*Total (cumulative) moment at this point*/
			int 	*nlayers) 		/*Number of decoupled layers*/
{
	/* RETURNS THE CALCULATED MOMENT IN THE PLATE WITH THE YIELD STRESS 
	ENVELOPE TAKING INTO ACCOUNT DECOUPLING IN SUBLAYERS */

	int	i, ix, iz, iter, j, k, itop, ibot, layer=0, 
		numlayers, numiter=20;
	float	stress_distrib_slope,		/*Slope of the linear part.*/
		cumulmoment=0, cumulmomentlayer,/*Total and layer cumulative moments.*/
		incremoment=0, incremomentlayer,/*Total and layer increment moments.*/
		decoupl_stress_limit=50e6, 	/*Default decoupling yield stress.*/
		yield_stress_minim=10e6, 	/*Minimum yield stress. This defines mechanical thickness. Ranalli, 1994.*/
		z, 				/*Depth.*/
		ztoplayer[10], zbotlayer[10],	/*Top & base of each layer.*/
		layerforceincre, 
		totalforceincre, 
		layerabsforceincre, 
		mecanthick=0, 
		Dsigma, 
		*newstress;
	BOOL	switch_saturatedlayer;

	/*Define the n decoupled layers (top and bottom).*/
	for (i=0; i<Nz; i++) {
		z = i*dz ;
		if (yieldextens[i]-yieldcompres[i] > 2*yield_stress_minim /*decoupl_stress_limit*/) {
			ztoplayer[layer] = z;
			/*fprintf(stdout, "\n\tTop of layer: %.2f", ztoplayer[layer]);*/
			for (; i<Nz; i++) {
				z = i*dz;
				Dsigma = yieldextens[i]-yieldcompres[i];
				if (
					(Dsigma < decoupl_stress_limit && isost_model == 6)
					|| (decoupl_depth>i*dz && decoupl_depth<=(i+1)*dz && isost_model == 4)
					|| i==Nz-1
				   ) {
					/*Base of layer is controlled by yield_stress_minim*/
					for (	iz=i; 
						iz>=0 && 
						(fabs(yieldextens[iz])<yield_stress_minim 
						|| fabs(yieldcompres[iz])<yield_stress_minim); 
						iz--
					) ;
					zbotlayer[layer] = iz*dz ;
					/*fprintf(stdout, "\tBase of layer: %.2f", zbotlayer[layer]);*/
					mecanthick += zbotlayer[layer] - ztoplayer[layer];
					layer++;
					break;
				}
			}
		}
	}
	numlayers=layer;

	totalforceincre = 0;
	if (!d2wdx2) d2wdx2 = 1e-12; /*very little*/
	stress_distrib_slope = - d2wdx2 * E / (1-nu*nu);
	newstress = calloc(Nz, sizeof(float));

	/*Distribute bending stresses along each decoupled layer.*/
	for (layer=0; layer<numlayers; layer++) {
		float 	z_null_strs = (zbotlayer[layer]+ztoplayer[layer])/2,
			layer_thick = zbotlayer[layer]-ztoplayer[layer];

		itop=ztoplayer[layer]/dz; ibot=zbotlayer[layer]/dz;
		for (j=itop, switch_saturatedlayer=YES; j<=ibot; j++) {
			if ((stress[j]>0 && stress[j]<yieldextens[j]) || (stress[j]<=0 && stress[j]>yieldcompres[j])) {
				switch_saturatedlayer=NO ;
			}
		}
		if (switch_saturatedlayer==NO) {
			/*Iterates until null stress point position reach convergence.*/
			for (i=1; i<=numiter; i++) {
				layerforceincre = layerabsforceincre = cumulmomentlayer = incremomentlayer = 0;
				for (iz=itop; iz<=ibot; iz++) {
					float linearstress;
					z = iz*dz;
					linearstress = stress_distrib_slope * (z-z_null_strs);
					newstress[iz] = stress[iz] + linearstress;
					if (newstress[iz]>yieldextens[iz])	newstress[iz] = yieldextens[iz];
					if (newstress[iz]<yieldcompres[iz])	newstress[iz] = yieldcompres[iz];
					layerforceincre += (newstress[iz]-stress[iz]) * dz;
					layerabsforceincre += fabs(newstress[iz]-stress[iz]) * dz;
					cumulmomentlayer += newstress[iz] * (z-z_null_strs) * dz;
					incremomentlayer += (newstress[iz]-stress[iz]) * (z-z_null_strs) * dz;
				}
				/*fprintf(stdout, "\n\tIter.%d; Pres.lay. %.2e; Abspres. %.2e;\tNullpoint:%.2f", i, layerforceincre, layerabsforceincre, z_null_strs);*/
				/*Exits if convergence ocurred; modifies null point depth otherwise:modifies*/
				if (fabs(layerforceincre) < .001e12) break;   /*if (fabs(layerforceincre)/layerabsforceincre < .001) break;*/
				z_null_strs -= (layerforceincre*d2wdx2)/fabs(layerforceincre*d2wdx2) * layer_thick/pow(2,i);
			}
			if (i>=numiter) fprintf(stdout, "*%d", layer);
			for (iz=itop; iz<=ibot; iz++)  {
				z = iz*dz;
				stress[iz] = newstress[iz];
			}
		}
		else {
			incremomentlayer = 0;
			cumulmomentlayer = 0;
		}

		cumulmoment += cumulmomentlayer;
		incremoment += incremomentlayer;
		totalforceincre += layerforceincre;
		/*if (verbose_level>=3) fprintf(stderr, "\n\t%d iters. en la capa %d, que tiene momento %.2e. Fuerza acumulada=%.2e", i, layer+1, momentlayer, totalforceincre);*/
	}
	free(newstress);

	*totalmoment = cumulmoment;
	*nlayers = numlayers ;
	return (incremoment) ;
}



int Repare_Blocks()
{
	int  	i_Block_max_arrange;

	/*Avoids Blocks to have negative thickness or zero area*/

	i_Block_max_arrange = (switch_topoest)? i_first_Block_load : numBlocks;
	for (int i_Block=1; i_Block < i_Block_max_arrange; i_Block++) {
		for (int i=0; i<Nx; i++) {
			Blocks[i_Block-1].thick[i] = MAX_2(Blocks[i_Block-1].thick[i], 0);
		}
	}
	/*Delete empty Blocks except if hidden*/
	for (int i_Block=0; i_Block<numBlocks; i_Block++) {
		float Block_area=0;
		for (int i=0; i<Nx; i++)  Block_area += Blocks[i_Block].thick[i];
		Block_area *= dx;
		if (Block_area<1e2 && Blocks[i_Block].type!='H' && Blocks[i_Block].type != 'I') {Delete_Block(i_Block); i_Block--;}
	}
	return(1);
}





int read_file_YSE()
{
	/*READS THE YIELD STRESS ENVELOPE FILE 'projectname.YSE'*/

	int 	iz, nz_input;
	FILE 	*file;
	char 	filename[MAXLENFILE];
	float	*yse_comp, *yse_extn, *z_yse, z;

	if (isost_model<3) return(0);
	sprintf(filename, "%s.YSE", projectname);
	if ((file = fopen(filename, "rt")) == NULL) {
		if (verbose_level>=3) fprintf(stderr, "\nYSE file '%s' not found. ", filename);
		switch_YSE_file = NO;
		return(0);
	}
	else if (verbose_level>=1) fprintf(stdout, "\nYSE at '%s'.", filename);

	yse_comp = calloc(10000, sizeof(float));
	yse_extn = calloc(10000, sizeof(float));
	z_yse =	calloc(10000, sizeof(float));
	nz_input = 0;
	for (;;) {
		TAKE_LINE_3(z_yse[nz_input], yse_comp[nz_input], yse_extn[nz_input]); 
		nz_input++;
	}
	fclose (file);
	switch_YSE_file = YES;

	/*Interpolates vertically the YSE*/
	for (iz=0; iz<Nz; iz++) {
		z = iz*dz;
		yieldcompres[iz] = interpol_in_xy_data(z_yse, yse_comp, nz_input, z);
		yieldextens[iz]  = interpol_in_xy_data(z_yse, yse_extn, nz_input, z);		
	}

	free(z_yse); free(yse_comp); free(yse_extn);
	return(1);
}




int Rheo_Flex_Iter () {
	int 	i, ix, iz, i_x_temp, j, rheoiter, 
		ncapas, NDs=3, NDi=3;
	float 	d2wdx2=SIGNAL, criterioconv, 
		Te_ant, 
		refstress, x, incremoment, 
		mechanical_thickness;	/*[m] */
	FILE 	*file, *file_temp ;
	char  	fileout[MAXLENFILE], gmtcommand[MAXLENLINE] ;
	float	*want, *x_stress, momentmax, *moment;
	double  **A, *b;

	/*
		SOLVES THE ELASTIC-PLASTIC FLEXURE WITH AN ITERATION THAT BRINGS
		TO AN EET AND DEFLECTION BOTH CONSISTENT WITH THE YSE. 
		See McNutt (1988) and BUROV & DIAMENT (1992, 1995).
	*/

	fprintf(stdout, " Rh  ");

	A = alloc_matrix_dbl (Nx, NDi+1+NDs);
	b = (double *) calloc (Nx , sizeof(double));
	moment = (float *) calloc (Nx , sizeof(float));
		if (!switch_strs_history) {
		/*
		  Note that in this case (no stress history) deflection is 
		  calculated with all the present load. Therefore, the 
		  resulting deflection is not added to the previous one 
		  and there is no possibility for viscous relaxation to 
		  be taken into account.
		*/
			float max_Te_var;
		want =   (float *) calloc (Nx , sizeof(float));

			for (i=0;i<Nx;i++) want[i] = w[i];

			/*Calculates a pure elastic Initial deflection*/
			LES_matrix(A, b, D, q, Dq, w, NO);
			solveLES(A, b, Nx, 3, 3, w);

			/*
			  Uses this w deflection as a first value for an iteration that
			  finds succesive EET and deflection distributions
			  until convergence is reached.
			*/
			for (rheoiter=0; rheoiter<NMAXRHEOITERS; rheoiter++) {
				max_Te_var=0;
		fprintf(stdout, "\b\b%2d", rheoiter);  fflush(stdout);
		/*For each x position calculates stress distribution & EET:*/
		for (ix=0, criterioconv=momentmax=0; ix<Nx; ix++) {
			x=x0+dx*ix;
			/*Curvature from previous deflection (positive at forebulge)*/
			if (ix>0 && ix<Nx-1)	d2wdx2 = (w[ix+1] -2*w[ix] + w[ix-1]) /dx/dx ;
			if (ix==0)			d2wdx2 = (Dw[2] -2*Dw[1] + Dw[0]) /dx/dx;
			if (ix==Nx-1) 		d2wdx2 = (Dw[Nx] -2*Dw[Nx-1] + Dw[Nx-2]) /dx/dx;
			if (fabs(d2wdx2)<1e-12) d2wdx2 = -1e-12;

			/*Finds the yield stress envelope from temperature
				and geometry:*/
			if (!switch_YSE_file) yield_stress_envelope (
				Temperature[ix], Nz, dz, 0, 
				upper_crust_thick[ix], crust_thick[ix], 
				isost_model, 
				yieldcompres, yieldextens, 
				&mechanical_thickness
			);

			/*Distributes stresses along vertical profile 
				depending on preliminar plate curvature and 
				returns the moment at this x position.*/
			moment[ix] = moment_calculator (
				d2wdx2, horz_force, 
				yieldcompres, yieldextens, 
				stress[ix], isost_model, crust_thick[ix], 
				Nz, dz, 
				&refstress, &ncapas
			);

			Te_ant = Te[ix] ;
			D[ix] = -moment[ix] / d2wdx2 ;	
			Te[ix] = RIG2ET(D[ix]) ;
			criterioconv += fabs(Te[ix] - Te_ant) ;
			max_Te_var=MAX_2(max_Te_var, fabs(Te[ix] - Te_ant));
			if (fabs(momentmax)<fabs(moment[ix])) {
				momentmax=moment[ix]; imomentmax=ix;
			}
		}

		/*Calculates new deflection with present EET*/
		LES_matrix(A, b, D, q, Dq, w, NO) ;
		solveLES(A, b, Nx, 3, 3, w) ;

		/*Checks convergence*/
		if ((criterioconv*dx < MAXETERR && max_Te_var<MAX_Te_LOC_VAR) || rheoiter>=NMAXRHEOITERS-1) break;
			}
			fprintf(stdout, "\b\b") ;
			if (rheoiter>=NMAXRHEOITERS-1) {
				fprintf(stdout, "! \b");
				if (verbose_level>=3)
					fprintf(stderr, "\nERROR: Lack of convergence in EET!. EET error area = %.2f km2", criterioconv*dx/1e6);
			}

			for (i=0;i<Nx;i++) Dw[i] = w[i] - want[i];
			free(want);
	}


	else {
			BOOL 	switch_last_repeat=NO;
			float 	point_moment, max_Te_var;

			x_stress = calloc (Nz , sizeof(float));

			/*Calculates a pure elastic initial deflection*/
			LES_matrix(A, b, D, q, Dq, w, NO);
			solveLES(A, b, Nx, 3, 3, Dw);

			/*
			  Uses this Dw deflection increment as a first value for an iteration
			  that finds succesive EET and deflection increments until
			  convergence is reached.
			*/
			for (rheoiter=0; rheoiter<NMAXRHEOITERS; rheoiter++) {
				max_Te_var=0;
		fprintf(stdout, "\b\b%2d", rheoiter);  fflush(stdout);
		/*For each x position calculates stress distribution & EET:*/
		for (ix=0, criterioconv=momentmax=0; ix<Nx; ix++) {
			x=x0+dx*ix;
			/*Curvature increment from previous deflection (positive at forebulge)*/
			if (ix>0 && ix<Nx-1) 	d2wdx2 = (Dw[ix+1] -2*Dw[ix] + Dw[ix-1]) /dx/dx;
			if (ix==0)	  	d2wdx2 = (Dw[2] -2*Dw[1] + Dw[0]) /dx/dx;
			if (ix==Nx-1)   	d2wdx2 = (Dw[Nx-1] -2*Dw[Nx-2] + Dw[Nx-3]) /dx/dx;
			if (fabs(d2wdx2)<1e-12) d2wdx2 = -1e-12;

			/*Finds the yield stress envelope from temperature
				and geometry:*/
			if (!switch_YSE_file) yield_stress_envelope (
				Temperature[ix], Nz, dz, 0, 
				upper_crust_thick[ix], crust_thick[ix], 
				isost_model, 
				yieldcompres, yieldextens, 
				&mechanical_thickness
			);

			/*
			  Distributes stresses along vertical profile 
			  depending on preliminar plate curvature and 
			  returns the moment at this x position.
			*/
			for (iz=0; iz<Nz; iz++) x_stress[iz]=stress[ix][iz];
			incremoment = moment_calculator_hist (
				d2wdx2, horz_force, 
				yieldcompres, yieldextens, 
				x_stress, isost_model, crust_thick[ix], 
				Nz, dz, 
				&point_moment, &ncapas
			);
			if (switch_last_repeat) {
				moment[ix] += incremoment;
				if (fabs(momentmax)<fabs(moment[ix])) {
					momentmax=moment[ix]; imomentmax=ix;
				}
				for (iz=0; iz<Nz; iz++) stress[ix][iz]=x_stress[iz];
			}

			Te_ant = Te[ix];
			D[ix] = -incremoment / d2wdx2;
			Te[ix] = RIG2ET(D[ix]);
			criterioconv += fabs(Te[ix] - Te_ant);
			max_Te_var=MAX_2(max_Te_var, fabs(Te[ix] - Te_ant));
		}

		/*Calculates new deflection with present EET*/
		LES_matrix(A, b, D, q, Dq, w, NO) ;
		solveLES(A, b, Nx, 3, 3, Dw) ;
		
		/*Checks convergence*/
		if ((criterioconv*dx < MAXETERR && max_Te_var<MAX_Te_LOC_VAR) || switch_last_repeat) {
			if (switch_last_repeat) break;
			else {switch_last_repeat=YES; rheoiter--;}
		}
			}
			fprintf(stdout, "\b\b") ;
			if (rheoiter>=NMAXRHEOITERS-1) {
				fprintf(stdout, "!(%s) \b", (criterioconv*dx<MAXETERR)? "M" : ((max_Te_var<MAX_Te_LOC_VAR)? "m" : "B"));
				if (verbose_level>=3)
					fprintf(stderr, "\nERROR: Lack of convergence in EET!. EET error area = %.2f km2  ", criterioconv*dx/1e6);
			}

			for (i=0;i<Nx;i++) w[i] += Dw[i];
			free(x_stress);
	}

	free (moment);
	return(1);
}





int solveLES (double **A, double *b, int Nx, int nds, int ndi, float *x)
{
	/*THIS ROUTINE CALLS THE ONES WHICH SOLVE THE LINEAR EQUATION SYSTEM */

	if (TriangularizeAlmostDiagonalEquationSystem(A, b, Nx, nds, ndi)) {
		fprintf(stdout, "\nUNDETERMINED SYSTEM !!!");
		getchar(); return 0;
	}
	SolveAlmostDiagonalTriangularEquationSystem(A, b, Nx, nds, ndi, x);
	return 1;
}




int yield_stress_envelope_semibrittle (
			float *Temperatura, 	/*Temperature array*/
			int   Nz, 		/*Number of depth finite differences points*/
			float dz, 		/*Depth gridding interval*/
			float z0, 		/*Equivalent meters of crust to determine the pressure and strength at the top of the plate*/
			float Uppercrustbase, 
			float Lowercrustbase, 
			int   isost_model, 	/*3: oceanic; >=4: Continental*/
			float *yieldcompres, 
			float *yieldextens, 
			float *espmecanptr )
{
	/*
		CALCULATES THE YIELD DIFFERENTIAL STRESS ENVELOPE FOR 
		COMPRESSIVE AND EXTENSIONAL STRESSES FROM GIVEN TEMPERATURE 
		DISTRIBUTION AND LITHOSPHERE VERTICAL GEOMETRY FOR A CERTAIN 
		X POSITION.
	*/

	int 	i, 
		numlayers=0;
		
	float	alpha_exten=.75, 
		alpha_compr=3, 
		alpha_strik=1.2, 
		lambda=.4, 
		Quc = 140e3, 			/*Power flow activation energies*/
		Qlc = 250e3, 
		Qml = 520e3, 
		epsuc = 2.5e-26, 		/*Reference power law strain rates in upper, lower crust and mantle in [s-1 Pa-3] */
		epslc = 3.2e-21, 
		epsml = 7e-14,
		n_exp = 3, 			/*Power law exponent*/
		R = 8.317, 			/*in J/mol/K */
		mecanlimit = 10e6, 		/*Mechanical thickness criteria [Pa]*/
		z, Q = SIGNAL, T, 
		brittle_exten, 
		brittle_compr, 
		brittle_strik, 
		espmecan=0;
	double	strainrate = 1e-16,		/*This is the assumed strain rate, which corresponds with a strain of .1 in 31.6 Ma */
		strainrateref = SIGNAL, 
		ductil_power_law, 
		ductil_Dorn_law, 
		ductil, 
		aux1,aux2,aux3;
	BOOL	switch_competente=NO;


	/*Oceans: Bodine et al. 1981 */
	if (isost_model==3) {
		Lowercrustbase = 3000;
		Qml = 522.5e3; 
	}

	for (i=0; i<Nz; i++) {
		z = i*dz;
		T = Temperatura[i]+273.3;
		if (z<Uppercrustbase) 				{Q=Quc; strainrateref=epsuc;}
		if (z>=Uppercrustbase && z<Lowercrustbase) 	{Q=Qlc; strainrateref=epslc;}
		if (z>=Lowercrustbase || isost_model==3) 	{Q=Qml; strainrateref=epsml;}

		/*Brittle failure: (Ranalli, p. 248)*/
		brittle_exten =  alpha_exten * g * denscrust * (z + z0) * (1-lambda);
		brittle_compr = -alpha_compr * g * denscrust * (z + z0) * (1-lambda);
		brittle_strik =  alpha_strik * g * denscrust * (z + z0) * (1-lambda);

		/*Ductil flow ('plastic'):*/
		/*Power law flow (Bodine et al. 1981; Lynch & Morgan, 1987):*/
		aux1 = strainrate / (double) strainrateref;
		ductil_power_law = pow(aux1, 1/n_exp) * (double) exp(Q/(n_exp*R*T));


		/*Limits stress to less than 200 MPa (semibrittle failure: Renshaw, 2004, JGR; Ranalli et al., 2006, Journal of Geodynamics)*/
		ductil = MIN_2(200e6, ductil_power_law);
			
		/*fprintf(stderr, "z=%.2f\tbc=%.2e\tbe=%.2e\tdf=%.2e\tdp=%.2e\tdd=%.2e\tT=%.2f\n", z, brittle_compr, brittle_exten, ductil, ductil_power_law, T);*/

		yieldcompres[i] = MAX_2(brittle_compr, -ductil) ;
		yieldextens[i]  = MIN_2(brittle_exten,  ductil) ;
		if (yieldextens[i] > mecanlimit || yieldcompres[i] < -mecanlimit || z<=Lowercrustbase) { 
			espmecan += dz ;
			if (!switch_competente) numlayers++ ; 
			switch_competente = YES;
		}
		else	switch_competente = NO;
	}
	*espmecanptr = espmecan;

	return (numlayers);
}



int yield_stress_envelope (
			float *Temperatura, 	/*Temperature array*/
			int   Nz, 		/*Number of depth finite differences points*/
			float dz, 		/*Depth gridding interval*/
			float z0, 		/*Equivalent meters of crust to determine the pressure and strength at the top of the plate*/
			float Uppercrustbase, 
			float Lowercrustbase, 
			int   isost_model, 	/*3: oceanic; >=4: Continental*/
			float *yieldcompres, 
			float *yieldextens, 
			float *espmecanptr )
{
	/*
		CALCULATES THE YIELD DIFFERENTIAL STRESS ENVELOPE FOR 
		COMPRESSIVE AND EXTENSIONAL STRESSES FROM GIVEN TEMPERATURE 
		DISTRIBUTION AND LITHOSPHERE VERTICAL GEOMETRY FOR A CERTAIN 
		X POSITION.
	Por cambiar aun: anyadir el efecto del peso de la capa de agua.
	*/

	int 	i, 
		numlayers=0;
		
	float	betaexten = 16e3, 		/*Brittle failure coefs.*/
		betacomp = -40e3,
		Quc = 140e3, 			/*Power flow activation energies*/
		Qlc = 250e3, 
		Qml = 520e3, 
		epsuc = 2.5e-26, 		/*Reference power law strain rates in upper, lower crust and mantle in [s-1 Pa-3] */
		epslc = 3.2e-21, 
		epsml = 7e-14,
		Dorn_Qml = 545e3, 		/*Dorn law activation energy; Sonder & England (1985) say 540e3: it's changed to fit with power law*/
		Dorn_strnrateref = 5.7e11, 
		Dorn_stressref = 8.5e9, 	/*[Pa] Bodine et al, 1981; Sonder & England, 1985*/
		n_exp = 3, 			/*Power law exponent*/
		R = 8.317, 			/*in J/mol/K */
		mecanlimit = 10e6, 		/*Mechanical thickness criteria [Pa]*/
		z, Q = SIGNAL, T, 
		brittleext, 
		brittlecomp, 
		espmecan=0;
	double	strainrate = 1e-16,		/*This is the assumed strain rate, which corresponds with a strain of .1 in 31.6 Ma */
		strainrateref = SIGNAL, 
		ductil_power_law, 
		ductil_Dorn_law, 
		ductil, 
		aux1,aux2,aux3;
	BOOL	switch_competente=NO;


	/*Oceans: Bodine et al. 1981 */
	if (isost_model==3) {
		Lowercrustbase = 3000;
		Qml = 522.5e3;			/*Power law*/
		Dorn_Qml = 535e3, 		/*Dorn law activation energy*/
		Dorn_stressref = 1.7019e10; 	/*[Pa] */
	}
	/*Continents: */
	else {
	}

	for (i=0; i<Nz; i++) {
		z = i*dz;
		T = Temperatura[i]+273.3;
		if (z<Uppercrustbase) 				{Q=Quc; strainrateref=epsuc;}
		if (z>=Uppercrustbase && z<Lowercrustbase) 	{Q=Qlc; strainrateref=epslc;}
		if (z>=Lowercrustbase || isost_model==3) 	{Q=Qml; strainrateref=epsml;}

		/*Brittle failure:*/
		brittleext = betaexten * (z + z0);
		brittlecomp = betacomp * (z + z0);

		/*Ductil flow:*/
		/*Power law flow (Bodine et al. 1981; Lynch & Morgan, 1987):*/
		aux1 = strainrate / (double) strainrateref;
		ductil_power_law = pow(aux1, 1/n_exp) * (double) exp(Q/(n_exp*R*T));

		/*Dorn flow when stress is higher than 200 MPa (Goetze & Evans, 1979; Bodine et al. 1981):*/
		ductil_Dorn_law = Dorn_stressref * (1 - sqrt(R*T/Dorn_Qml*log(Dorn_strnrateref/strainrate)));	if (ductil_Dorn_law<0) ductil_Dorn_law = 0;


		/*Power and Dorn laws have to coincide at T=983.26 K, i.e., 200 MPa*/
		if (ductil_power_law>200e6 && ductil_Dorn_law>0 && z>=Lowercrustbase) 
			ductil = ductil_Dorn_law;
		else
			ductil = ductil_power_law;
			
		/*
		fprintf(stderr, "z=%.2f\tbc=%.2e\tbe=%.2e\tdf=%.2e\tdp=%.2e\tdd=%.2e\tT=%.2f\n", z, brittlecomp, brittleext, ductil, ductil_power_law, ductil_Dorn_law, T);
		*/

		yieldcompres[i] = (brittlecomp > -ductil)? 	brittlecomp : -ductil ;
		yieldextens[i]  = (brittleext < ductil)?	brittleext  :  ductil ;
		if (yieldextens[i] > mecanlimit || yieldcompres[i] < -mecanlimit || z<=Lowercrustbase) { 
			espmecan += dz ;
			if (!switch_competente) numlayers++ ; 
			switch_competente = YES;
		}
		else	switch_competente = NO;
	}
	*espmecanptr = espmecan;

	return (numlayers);
}




float calculate_sea_level()
{
	/*
	  Calculates the sea level and the load related to changes 
	  in the water column (sea and lakes).
	*/

	if (!water_load) return (0);

	/*Calculates sea level*/
	if (n_sea_level_input_points) {
		int i;
		for (i=0; i<n_sea_level_input_points; i++) {
		if (var_sea_level[i][0]>=Time) break;
			}
			if (i!=0 && i!=n_sea_level_input_points) {
		sea_level = 	( (Time-var_sea_level[i-1][0])*var_sea_level[i][1] + 
				  (var_sea_level[i][0]-Time)*var_sea_level[i-1][1] ) 
					 / (var_sea_level[i][0]-var_sea_level[i-1][0]);
		}
 		else {
		if (i==0)			 sea_level = var_sea_level[0][1];
		if (i==n_sea_level_input_points) sea_level = var_sea_level[n_sea_level_input_points-1][1];
		}
	}
	else sea_level = 0;
	return(sea_level);
}



int calculate_water_load()
{
	float	water_volume=0;

	/*
	  Calculates sea/lake water load due to the water column
	*/

	calculate_topo(topo);

	if (!water_load) return(0);

	for (int i=0; i<Nx; i++) {
		int il;
		float Dq_water, h_water_now=0;
		if (hydro_model) {
			if (il=drainage[i].lake) {
				/*sea lake already has its proper level defined*/
			h_water_now = MAX_2(0, Lake[il].alt-topo[i]);
			}
		}
		else {
			h_water_now = MAX_2(0, sea_level-topo[i]);
		}
		Dq_water = (h_water_now-h_water[i]) * g * (denswater-densenv);
//PRINT_ERROR("%.2f, %.2f %.2f %d %.2f %.2f %.2e", i*dx+x0, h_water_now, h_water[i], il, sea_level, topo[i], Dq_water);
		h_water[i] = h_water_now;
		/*Don't load the initial water column*/
		if (Time>Timeini) Dq[i] += Dq_water;
		water_volume += h_water[i];		
	}

	PRINT_ARRAY_INFO(h_water, "water", "m", "m2") 
	if (n_sea_level_input_points) {
		PRINT_SUMLINE("sea_level: %8.1f m   sea_volume = %.1f km2", sea_level, water_volume*dx*Nx/(Nx-1)/1e6);
	}
	fflush(stdout);
	return(1);
}




