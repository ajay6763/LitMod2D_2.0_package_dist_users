/*
INCLUDE FILE FOR tao.c
*/

#include "geomodel.h"			/*General definitions and types for geophysical models*/


#define	MAXETERR	10e6		/*Maximum error area allowed as Te convergence criteria in the stress distrib. algorithm (m2)*/
#define	MAX_Te_LOC_VAR	200		/*Maximum error in Te at any point. Also a convergence criteria in the stress distrib. algorithm (m2)*/
#define	NMAXRHEOITERS 	100		/*Maximum number of iterations to find Te convergence in the stress distrib. algorithm*/
#define	NmaxBlocks	400			/*Maximum number of Blocks to be recorded*/
#define Read_Param(x, y)  {char line[MAXLENLINE+200]; while ((sscanf(lineptr=fgets(line, MAXLENLINE+200, file), x, y )) < 1) \
	if (lineptr==NULL) {fprintf(stderr, "\nERROR reading parameters file: parameter not found.\n");	break;}}

#define PRINT_ARRAY_INFO(array, name, units, unitsintegr) {\
	    int i,imin=SIGNAL,imax=SIGNAL; float max=-1e19, min=1e19, vol=0;\
	    for (i=(xmin-x0)/dx; i<Nx-(xf-xmax)/dx; i++) {\
	    	    vol += array[i]*dx;\
	    	    if (max<array[i])  {max=array[i]; imax=i;}\
	    	    if (min>array[i])  {min=array[i]; imin=i;}}\
	    PRINT_SUMLINE("%s:\tmax= %.2e %s @x= %.1f km\tmin= %.2e %s @x= %.1f km\tintegr= %.2e %s", name, max, units, (x0+imax*dx)/1e3, min, units, (x0+imin*dx)/1e3, vol, unitsintegr);}

/*convert between sediment thickness and sediment grain mass*/
#define MASS2SEDTHICK_1D(mass)	((mass) /(denssedim-sed_porosity*denswater)/dx/riverbasinwidth)	/*converts sediment mass into sediment thickness*/
#define THICK2SEDMASS_1D(thick)	((thick)*(denssedim-sed_porosity*denswater)*dx*riverbasinwidth)	/*converts sediment thickness into sediment mass*/
#define GET_KREST(Krest, q, i) { /*Define the restoring force value.*/\
	    if (switch_topoest) {\
		/*If the current i knot is below the load then the compensation density is densinfill.*/\
		if (q[i])   Krest = (densasthen-densinfill)*g;\
		/*Otherwise use the sediment density.*/\
		else	    Krest = (densasthen-denssedim)*g;\
	    } else    Krest = (densasthen-densenv)*g;}

int *sortcell;
struct DRAINAGE_1D  *drainage;
struct LAKE_INFO_1D *Lake;	/*Lake[0] does not exist; Lake[1] is the sea or the first normal lake.*/
struct BLOCK_1D	*Blocks ; 	/*Blocks array: Blocks[numBlocks] first Block is 0, last is numBlocks-1;*/

int 	imomentmax, 
	nmax_input_points,	/*Maximum number of input points in moving loads*/
	nx_temp_input;		/*Number of x points with defined geotherm in .TMP file*/

char	eros_bound_cond[2];
int	erosed_model, 
	hydro_model, 
	mode_interp, 		
	nbasins, 
	nlakes, 		/*number of lakes >= 0 */
	boundary_conds,
	n_image;

float 	evaporation_ct, 		/*[m3/s/m2].*/
	riverbasinwidth, 		/*[m]*/
	total_evap_water,
	total_lost_water,
	total_rain, 
	x0, xf, 			/*X horizontal direction logical limits of the model [m]*/
	xmin, xmax, zmin, zmax, 	/*Screen limits for the graphic output [m]*/
	horz_force, vert_force, appmoment; 	/*Horizontal tectonic force and Vertical boundary applied force (in N/m); Boundary applied moment (in N m/m)*/


/*Arrays:*/
float	*crust_thick,			/*Crust thickness array*/
	*D, 				/*Rigidity		[N m] */
	*Dq, 				/*Load increment	[N/m/m] */
	*Dw, 				/*Deflection increment	[m] */
	*h_water, 				/*Water load thickness	[m] */
	*h_last_unit,		 	/*Height at last unit file [m] */
	*q, 				/*load 			[N/m/m] */
	*Te,				/*E.E.T. array*/
	*eros_now,
	*precipitation,
	*evaporation,
	*total_erosion, 		/*Cumulated erosion at each point*/
	*topo,  			/*Topography over sea level*/
	*Blocks_base,			/*Base of Blocks[0] (the lowest) measured from the current position of the original datum surface (now deflected by isostasy).*/
	*upper_crust_thick,		/*Upper crust thickness array*/
	*w, 				/*Deflection		[m] */
	**Temperature, 			/*[ºC]*/
	**stress, 			/*[Pa]*/
	*yieldcompres,			/*yield stress envelope -- compresive  [Pa]*/
	*yieldextens;			/*yield stress envelope -- extensive  [Pa]*/


char 	gif_geom[MAXLENLINE];


/*Boolean type variables:*/
BOOL	switch_insert_load=NO,		/*Si to insert the load height beneath the first load Block*/
	switch_strs_history, 
	switch_YSE_file=NO;			/*SI if yield stress is given directly from a .YSE file*/




/*Function prototypes in libreria.c:*/
int yield_stress_envelope (
    float *Temperatura,     /*Temperature array*/
    int   Nz,		    /*Number of depth finite differences points*/
    float dz,		    /*Depth gridding interval*/
    float z0, 
    float Uppercrustbase,
    float Lowercrustbase,
    int   lith_type,	    /*2: oceanic; >=3: Continental*/
    float *yieldcompres,
    float *yieldextens,
    float *espmecanptr 
);

/*Function prototypes in libreria.c:*/
float **alloc_matrix  	(int num_fil, int num_col);
double **alloc_matrix_dbl 	(int num_fil, int num_col);
float geoidanompolig 		(float *, float *, int, float, float, float);
float gravanompolig		(float *, float *, int, float, float, float);
int	readinterplin		(FILE *, float *, int, float, float);
float interpol_in_xy_data 	(float *x_array, float *y_array, int n_x, float x);

/*Declaration of functions in taosp.c*/
int Surface_Transport (float *topo, float dt, float dt_eros, int erosed_model);
int Fluvial_Transport (float *topo, float dt, float dt_eros, int erosed_model);

int Add_Node_To_Lake (int row, int i_lake);
int Add_Saddle_To_Lake (int row_sd, int row_tr, int i_lake);
int Attempt_Delete_Node_From_Lake (int row);
int constant_rate_eros (float *topo, float Keroseol, float Ksedim, float sea_level, BOOL switch_sea, 	float dt, float dt_eros, float Time, 	int n_eros_level_input_points, float **var_eros_level, 	float *eros_level);
int Deallocate_Lake (int i_lake);
int Define_Lake (int i_lake);
int Delete_Node_From_Lake (int row);
int Diffusive_Eros_1D (float *topo, float Kerosdif, float dt, float dt_eros); 
int Divide_Lake (int row);
int Erode (double d_mass, int row);
float Lake_Input_Discharge (int ilake);
int	Lake_Node_Number(int row);
int	Lake_Saddle_Number (int row);
int move_Blocks();
int	New_Lake ();
int	Repare_Blocks();
int	Sediment (double d_mass, int row);
int	tectload();
int	Unify_Lakes (int i_lake, int i_lake_to_delete);
int	Landslide_Transport (float critical_slope, float dt, float dt_eros);



#include "tao+tisc.h"
