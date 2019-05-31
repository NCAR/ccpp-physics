#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include "fregrid_util.h"
#include "mpp.h"
#include "mpp_io.h"
#include "tool_util.h"
#include "mosaic_util.h"
#include "read_mosaic.h"
#include "gradient_c2l.h"
#include "globals.h"
#include "interp.h"

#define D2R (M_PI/180)
#define R2D (180/M_PI)
#define EPSLN10 (1.e-10)
#define REL_COEF ( 0.9 )
#define MAX_ITER 4000

void init_halo(double *var, int nx, int ny, int nz, int halo);
void update_halo(int nx, int ny, int nz, double *data, Bound_config *bound, Data_holder *dHold);
void setup_boundary(const char *mosaic_file, int ntiles, Grid_config *grid, Bound_config *bound, int halo, int position);
void delete_bound_memory(int ntiles, Bound_config *bound);
void copy_var_config(const Var_config *var_in, Var_config *var_out);
void init_var_config(Var_config *var, int interp_method);
void fill_boundaries(int ni, int nj, double *data, int is_cyclic);
int parse_string(const char *str1, const char *str2, char *strOut, char *errmsg);
void do_extrapolate (int ni, int nj, int nk, const double *lon, const double *lat, const double *data_in,
		     double *data_out, int is_cyclic, double missing_value, double stop_crit);
/*******************************************************************************
  void setup_tile_data_file(Mosaic_config mosaic, const char *filename)
  This routine will setup the data file name for each tile.
*******************************************************************************/

void set_mosaic_data_file(int ntiles, const char *mosaic_file, const char *dir, File_config *file,
			  const char *filename)
{
  char   str1[STRING]="", str2[STRING]="", tilename[STRING]="";
  int    i, n, len, fid, vid;
  size_t start[4], nread[4];

  len = strlen(filename); 
  if( strcmp(filename+len-3, ".nc") ==0 ) 
    strncpy(str1, filename, len-3);
  else
    strcpy(str1, filename);
  if(dir) {
    if(strlen(dir)+strlen(str1) >= STRING)mpp_error("set_mosaic_data_file(fregrid_util): length of str1 + "
						    "length of dir should be no greater than STRING");
    sprintf(str2, "%s/%s", dir, str1);
  }
  else
    strcpy(str2, str1);
  
  for(i=0; i<4; i++) {
    start[i] = 0; nread[i] = 1;
  }
  nread[1] = STRING;
  if(ntiles > 1) {
    if(!mosaic_file) mpp_error("fregrid_util: when ntiles is greater than 1, mosaic_file should be defined");
    fid = mpp_open(mosaic_file, MPP_READ);
    vid = mpp_get_varid(fid, "gridtiles");
  }
  for(i = 0; i < ntiles; i++) {
    start[0] = i;
    if(ntiles > 1) {
      mpp_get_var_value_block(fid, vid, start, nread, tilename);
      if(strlen(str2) + strlen(tilename) > STRING -5) mpp_error("set_mosaic_data_file(fregrid_util): length of str2 + "
								"length of tilename should be no greater than STRING-5");
      sprintf(file[i].name, "%s.%s.nc", str2, tilename);
    }
    else 
      sprintf(file[i].name, "%s.nc", str2);
  }

}; /* setup_data_file */

/*******************************************************************************
  void set_scalar_var()
*******************************************************************************/
void set_field_struct(int ntiles, Field_config *field, int nvar, char * varname, File_config *file)
{
  int  n, i;
  
  if(nvar == 0) return;
 
  for(n=0; n<ntiles; n++) {
    field[n].area = NULL;
    field[n].file = file[n].name;
    field[n].fid = &(file[n].fid);
    field[n].nvar = nvar;
    field[n].var = (Var_config *)malloc(nvar*sizeof(Var_config)); 
    for(i=0; i<nvar; i++) 
      strcpy(field[n].var[i].name, varname+i*STRING);
  }

}; /* set_field_var */


/*******************************************************************************
  void set_weight_inf(int ntiles, Grid_config *grid, const char *weight_file, const char *weight_field);
  read the weight information
*******************************************************************************/
void set_weight_inf(int ntiles, Grid_config *grid, const char *weight_file, const char *weight_field,
		    int has_cell_measure_att)
{
  int n, fid, vid;
  char file[512];
  
  for(n=0; n<ntiles; n++) grid[n].weight_exist = 0;
  if( !weight_field || has_cell_measure_att ) return;
  
  for(n=0; n<ntiles; n++) {
    if(ntiles==1)
      strcpy( file, weight_file );
    else
      sprintf(file, "%s.tile%d.nc", weight_file, n+1);

    grid[n].weight_exist = 1;
    grid[n].weight = (double *)malloc(grid[n].nx*grid[n].ny*sizeof(double));
    fid = mpp_open(file, MPP_READ);
    vid = mpp_get_varid(fid, weight_field);
    mpp_get_var_value(fid, vid, grid[n].weight);
    mpp_close(fid);
  }
    

};/* set_weight_inf */


/*******************************************************************************
  void get_input_grid()
  To decrease memory usage, use grid_out latitude to decide the input grid to be stored.
*******************************************************************************/
void get_input_grid(int ntiles, Grid_config *grid, Bound_config *bound_T, const char *mosaic_file,
		    unsigned int opcode, int *great_circle_algorithm, int save_weight_only)
{
  int         n, m1, m2, i, j, l, ind1, ind2, nlon, nlat;
  int         ts, tw, tn, te, halo, nbound;
  int         m_fid, g_fid, vid;
  int         *nx, *ny;
  double      *x, *y;
  char         grid_file[256], filename[256], dir[256];
  size_t        start[4], nread[4];
  Data_holder *dHold;
  Bound_config *bound_C;
  int          read_tgrid;
  
  halo = 0;
  if(opcode & BILINEAR) halo = 1;
  for(n=0; n<4; n++) {
    start[n] = 0; nread[n] = 1;
  }

  read_tgrid = 1;
  if(save_weight_only || (opcode & CONSERVE_ORDER1) ) read_tgrid = 0;
  if(opcode & EXTRAPOLATE ) read_tgrid = 1;
  bound_C = (Bound_config *)malloc(ntiles*sizeof(Bound_config));
  nx = (int *)malloc(ntiles * sizeof(int) );
  ny = (int *)malloc(ntiles * sizeof(int) );
  
  *great_circle_algorithm = 0;
  m_fid = mpp_open(mosaic_file, MPP_READ);
  get_file_path(mosaic_file, dir);
  for(n=0; n<ntiles; n++) {
    start[0] = n; start[1] = 0; nread[0] = 1; nread[1] = STRING;
    vid = mpp_get_varid(m_fid, "gridfiles");
    mpp_get_var_value_block(m_fid, vid, start, nread, filename);
    sprintf(grid_file, "%s/%s", dir, filename);
    g_fid = mpp_open(grid_file, MPP_READ);

    if(n==0) *great_circle_algorithm = get_great_circle_algorithm(g_fid);
    nx[n] = mpp_get_dimlen(g_fid, "nx");
    ny[n] = mpp_get_dimlen(g_fid, "ny");
    if(nx[n]%2) mpp_error("fregrid_util(get_input_grid): the size of dimension nx should be even (on supergrid)");
    if(ny[n]%2) mpp_error("fregrid_util(get_input_grid): the size of dimension ny should be even (on supergrid)");
    nx[n] /= 2;
    ny[n] /= 2;
    grid[n].halo = halo;
    grid[n].nx   = nx[n];
    grid[n].ny   = ny[n];
    grid[n].nxc  = nx[n];
    grid[n].nyc  = ny[n];
    grid[n].is_cyclic = 0;    
    /* get supergrid */
    x = (double *)malloc((2*nx[n]+1)*(2*ny[n]+1)*sizeof(double));
    y = (double *)malloc((2*nx[n]+1)*(2*ny[n]+1)*sizeof(double));

    vid = mpp_get_varid(g_fid, "x");
    mpp_get_var_value(g_fid, vid, x);
    vid = mpp_get_varid(g_fid, "y");
    mpp_get_var_value(g_fid, vid, y);

    grid[n].lonc = (double *) malloc((nx[n]+1+2*halo)*(ny[n]+1+2*halo)*sizeof(double));
    grid[n].latc = (double *) malloc((nx[n]+1+2*halo)*(ny[n]+1+2*halo)*sizeof(double));
    grid[n].lont1D = (double *) malloc(nx[n]*sizeof(double));
    grid[n].latt1D = (double *) malloc(ny[n]*sizeof(double));
    for(i=0; i<nx[n]; i++) grid[n].lont1D[i] = x[2*nx[n]+1+2*i+1] * D2R;
    for(j=0; j<ny[n]; j++) grid[n].latt1D[j] = y[(2*j+1)*(2*nx[n]+1) + 1] * D2R;
    
    if(halo>0) {
      init_halo(grid[n].lonc, nx[n]+1, ny[n]+1, 1, halo);
      init_halo(grid[n].latc, nx[n]+1, ny[n]+1, 1, halo);
    }
    for(j=0; j<=ny[n]; j++) for(i=0; i<=nx[n]; i++) {
      ind1 = (j+halo)*(nx[n]+1+2*halo)+i+halo;
      ind2 = 2*j*(2*nx[n]+1)+2*i;
      grid[n].lonc[ind1] = x[ind2]*D2R;
      grid[n].latc[ind1] = y[ind2]*D2R;
    }
    if(read_tgrid) {
      grid[n].lont = (double *) malloc((nx[n]+2)*(ny[n]+2)*sizeof(double));
      grid[n].latt = (double *) malloc((nx[n]+2)*(ny[n]+2)*sizeof(double));
      for(j=0; j<ny[n]; j++) for(i=0; i<nx[n]; i++) {
	ind1 = (j+1)*(nx[n]+2)+i+1;
	ind2 = (2*j+1)*(2*nx[n]+1)+2*i+1;
	grid[n].lont[ind1] = x[ind2]*D2R;
	grid[n].latt[ind1] = y[ind2]*D2R;
      }
    
      init_halo(grid[n].lont, nx[n], ny[n], 1, 1);
      init_halo(grid[n].latt, nx[n], ny[n], 1, 1);
    }
    
    /* if vector, need to get rotation angle */
    /* we assume the grid is orthogonal */
    if( opcode & VECTOR ) {
      if( opcode & AGRID) {
	double *angle;
      	angle          = (double *) malloc((2*nx[n]+1)*(2*ny[n]+1)*sizeof(double));
	grid[n].cosrot = (double *) malloc(nx[n]*ny[n]*sizeof(double));	 
	grid[n].sinrot = (double *) malloc(nx[n]*ny[n]*sizeof(double));
	vid = mpp_get_varid(g_fid, "angle_dx");
	mpp_get_var_value(g_fid, vid, angle);
	grid[n].rotate = 0;
	for(j=0; j<ny[n]; j++) for(i=0; i<nx[n]; i++) {
          m1 = j*nx[n]+i;
	  m2 = (2*j+2)*(2*nx[n]+2)+2*i+2;
	  grid[n].cosrot[m1] = cos(angle[m2]*D2R);
	  grid[n].sinrot[m1] = sin(angle[m2]*D2R);
	  if(fabs(grid[n].sinrot[m1]) > EPSLN10) grid[n].rotate = 1;
	}
	free(angle);
      }
    }
    free(x);
    free(y);
    mpp_close(g_fid);
  }

  mpp_close(m_fid);

  /* get the boundary condition */
  setup_boundary(mosaic_file, ntiles, grid, bound_T, 1, CENTER);
  if(read_tgrid) {
    for(n=0; n<ntiles; n++) {
      nlon = grid[n].nx;
      nlat = grid[n].ny;
      nbound = bound_T[n].nbound;
      if(nbound > 0 ) {
	dHold = (Data_holder *)malloc(nbound*sizeof(Data_holder));
	for(l=0; l<nbound; l++) {
	  dHold[l].data = grid[bound_T[n].tile2[l]].lont;
	  dHold[l].nx = grid[bound_T[n].tile2[l]].nx+2;
	  dHold[l].ny = grid[bound_T[n].tile2[l]].ny+2;
	}
	update_halo(nlon+2, nlat+2, 1, grid[n].lont, &(bound_T[n]), dHold );
	for(l=0; l<nbound; l++) dHold[l].data = grid[bound_T[n].tile2[l]].latt;
	update_halo(nlon+2, nlat+2, 1, grid[n].latt, &(bound_T[n]), dHold );
	
	for(l=0; l<nbound; l++) dHold[l].data = NULL;
	free(dHold);
      }
    }
  }
  
  /* for bilinear interpolation, need to get cell-center grid */
  if(opcode & BILINEAR) {
    setup_boundary(mosaic_file, ntiles, grid, bound_C, 1, CORNER);
    /*--- fill the halo of corner cell */
    for(n=0; n<ntiles; n++) {
      nlon = grid[n].nx;
      nlat = grid[n].ny;
      grid[n].vlon_t = (double *) malloc(3*(nlon+2)*(nlat+2)*sizeof(double));
      grid[n].vlat_t = (double *) malloc(3*(nlon+2)*(nlat+2)*sizeof(double));
      grid[n].xt     = (double *) malloc(  (nlon+2)*(nlat+2)*sizeof(double));
      grid[n].yt     = (double *) malloc(  (nlon+2)*(nlat+2)*sizeof(double));
      grid[n].zt     = (double *) malloc(  (nlon+2)*(nlat+2)*sizeof(double));
      latlon2xyz((nlon+2)*(nlat+2), grid[n].lont, grid[n].latt, grid[n].xt, grid[n].yt, grid[n].zt);
      unit_vect_latlon((nlon+2)*(nlat+2), grid[n].lont, grid[n].latt, grid[n].vlon_t, grid[n].vlat_t);
      nbound = bound_C[n].nbound;
      if(nbound > 0) {
	dHold = (Data_holder *)malloc(nbound*sizeof(Data_holder));
	for(l=0; l<nbound; l++) {
	  dHold[l].data = grid[bound_C[n].tile2[l]].lonc;
	  dHold[l].nx   = grid[bound_C[n].tile2[l]].nx + 1 + 2*halo;
	  dHold[l].ny   = grid[bound_C[n].tile2[l]].ny + 1 + 2*halo;
	}
	update_halo(nlon+1+2*halo, nlat+1+2*halo, 1, grid[n].lonc, &(bound_C[n]), dHold);
	for(l=0; l<nbound; l++) dHold[l].data = grid[bound_C[n].tile2[l]].latc;
	
	update_halo(nlon+1+2*halo, nlat+1+2*halo, 1, grid[n].latc, &(bound_C[n]), dHold);
	for(l=0; l<nbound; l++) dHold[l].data = NULL;
	free(dHold);
      }
    }
  }
  else if(opcode & CONSERVE_ORDER2 && (!save_weight_only) ) {
    double p1[3], p2[3], p3[3], p4[3];
    
    for(n=0; n<ntiles; n++) {
      int is_true = 1;
      nlon = grid[n].nx;
      nlat = grid[n].ny;
      /* calculate dx, dy, area */
      grid[n].dx     = (double *)malloc(nlon    *(nlat+1)*sizeof(double));
      grid[n].dy     = (double *)malloc((nlon+1)*nlat    *sizeof(double));
      grid[n].area   = (double *)malloc(nlon    *nlat    *sizeof(double));
      grid[n].edge_w = (double *)malloc(         (nlat+1)*sizeof(double));
      grid[n].edge_e = (double *)malloc(         (nlat+1)*sizeof(double));
      grid[n].edge_s = (double *)malloc((nlon+1)         *sizeof(double));
      grid[n].edge_n = (double *)malloc((nlon+1)         *sizeof(double));
      grid[n].en_n   = (double *)malloc(3*nlon  *(nlat+1)*sizeof(double));
      grid[n].en_e   = (double *)malloc(3*(nlon+1)*nlat  *sizeof(double));
      grid[n].vlon_t = (double *) malloc(3*nlon*nlat     *sizeof(double));
      grid[n].vlat_t = (double *) malloc(3*nlon*nlat     *sizeof(double));
      calc_c2l_grid_info(&nlon, &nlat, grid[n].lont, grid[n].latt, grid[n].lonc, grid[n].latc,
			 grid[n].dx, grid[n].dy, grid[n].area, grid[n].edge_w, grid[n].edge_e,
			 grid[n].edge_s, grid[n].edge_n, grid[n].en_n, grid[n].en_e,
			 grid[n].vlon_t, grid[n].vlat_t, &is_true, &is_true, &is_true, &is_true);
    }
  }

  if(opcode & BILINEAR) delete_bound_memory(ntiles, bound_C);
  free(bound_C);
  free(nx);
  free(ny);
}; /* get_input_grid */ 


void get_input_output_cell_area(int ntiles_in, Grid_config *grid_in, int ntiles_out, Grid_config *grid_out, unsigned int opcode)
{
  double *lonc=NULL, *latc=NULL;
  int n, nx, ny, ind1, ind2, i, j, halo;

  halo = grid_in->halo;
  for(n=0; n<ntiles_in; n++) {
    nx = grid_in[n].nx;
    ny = grid_in[n].ny;

    grid_in[n].cell_area = (double *)malloc(nx*ny*sizeof(double));
    lonc = (double *)malloc((nx+1)*(ny+1)*sizeof(double));
    latc = (double *)malloc((nx+1)*(ny+1)*sizeof(double));
      
    for(j=0; j<=ny; j++) for(i=0; i<=nx; i++) {
      ind1 = j*(nx+1)+i;
      ind2 = (j+halo)*(nx+1+2*halo)+i+halo;
      lonc[ind1] = grid_in[n].lonc[ind2];
      latc[ind1] = grid_in[n].latc[ind2];
    }

    /*calculate grid_in cell area */
    if( opcode & GREAT_CIRCLE ) 
      get_grid_great_circle_area(&nx, &ny, lonc, latc, grid_in[n].cell_area);
    else 
      get_grid_area(&nx, &ny, lonc, latc, grid_in[n].cell_area);  

    free(lonc);
    free(latc);
  }

  for(n=0; n<ntiles_out; n++) {
    nx = grid_out[n].nxc;
    ny = grid_out[n].nyc;

    grid_out[n].cell_area = (double *)malloc(nx*ny*sizeof(double));
      
    /*calculate grid_in cell area */
    if( opcode & GREAT_CIRCLE ) 
      get_grid_great_circle_area(&nx, &ny, grid_out[n].lonc, grid_out[n].latc, grid_out[n].cell_area);
    else 
      get_grid_area(&nx, &ny, grid_out[n].lonc, grid_out[n].latc, grid_out[n].cell_area);  

  }
  
}

/*******************************************************************************
  void get_output_grid_from_mosaic(Mosaic_config *mosaic)

*******************************************************************************/
void get_output_grid_from_mosaic(int ntiles, Grid_config *grid, const char *mosaic_file, unsigned int opcode,
				 int *great_circle_algorithm)
{
  int         n, i, j, ii, jj, npes, layout[2];
  int         m_fid, g_fid, vid, ind;
  int         *nx, *ny;
  double      *x, *y;
  size_t        start[4], nread[4];
  char         grid_file[256], filename[256], dir[256];
  
  if(opcode & BILINEAR) mpp_error("fregrid_util: for bilinear interpolation, output grid can not got from mosaic file");
  
  npes   = mpp_npes();
  nx = (int *)malloc(ntiles * sizeof(int) );
  ny = (int *)malloc(ntiles * sizeof(int) );

  for(n=0; n<4; n++) {
    start[n] = 0;
    nread[n] = 1;
  }

  *great_circle_algorithm = 0;

  /* check if the grid is tripolar grid or not */
  for(n=0; n<ntiles; n++) grid[n].is_tripolar = 0;
  if(ntiles == 1 && mpp_field_exist(mosaic_file, "contacts") ) {
    int ncontacts, m;
    char errmsg[512];
    int tile1[2], tile2[2];
    int istart1[2], iend1[2], jstart1[2], jend1[2];
    int istart2[2], iend2[2], jstart2[2], jend2[2];

    ncontacts = read_mosaic_ncontacts(mosaic_file);
    if(ncontacts < 1) {
      sprintf(errmsg, "fregrid_util.c: number of contacts should be larger than 0 when field contacts exist in file %s",
	      mosaic_file );
      mpp_error(errmsg);
    }
    if(ncontacts > 2) {
      sprintf(errmsg, "fregrid_util.c: "
	      "number of contacts should be no larger than 2 in file %s",mosaic_file );
      mpp_error(errmsg);
    }
    read_mosaic_contact( mosaic_file, tile1, tile2, istart1, iend1, jstart1, jend1,
			 istart2, iend2, jstart2, jend2 );
      
    for(m=0; m<ncontacts; m++) {
      if( jstart1[m] == jend1[m] ) {  /* y-direction contact, cyclic or folded-north */
	if(jstart2[m] != jend2[m] )
	  mpp_error("fregrid_util.c: only cyclic/folded-north condition is allowed for y-boundary");
	if( jstart1[m] == jstart2[m] )  /* folded north */
	  grid->is_tripolar = 1;
      }
    }
  }
  
  m_fid = mpp_open(mosaic_file, MPP_READ);
  get_file_path(mosaic_file, dir);
  
  for(n=0; n<ntiles; n++) {
    start[0] = n; start[1] = 0; nread[0] = 1; nread[1] = STRING;
    vid = mpp_get_varid(m_fid, "gridfiles");
    mpp_get_var_value_block(m_fid, vid, start, nread, filename);
    sprintf(grid_file, "%s/%s", dir, filename);
    g_fid = mpp_open(grid_file, MPP_READ);
    if(n==0) *great_circle_algorithm = get_great_circle_algorithm(g_fid);
    nx[n] = mpp_get_dimlen(g_fid, "nx");
    ny[n] = mpp_get_dimlen(g_fid, "ny");
    if(nx[n]%2) mpp_error("fregrid_util(get_output_grid_from_mosaic): the size of dimension nx should be even (on supergrid)");
    if(ny[n]%2) mpp_error("fregrid_util(get_output_grid_from_mosaic): the size of dimension ny should be even (on supergrid)");
    nx[n] /= 2;
    ny[n] /= 2;
    grid[n].nx = nx[n];
    grid[n].ny = ny[n];
    /* to be able to reprocessor count, layout need to be set as follwoing */
    layout[0] = 1;
    layout[1] = npes;
    mpp_define_domain2d(grid[n].nx, grid[n].ny, layout, 0, 0, &(grid[n].domain));
    mpp_get_compute_domain2d(grid[n].domain, &(grid[n].isc), &(grid[n].iec), &(grid[n].jsc), &(grid[n].jec));
    grid[n].nxc = grid[n].iec - grid[n].isc + 1;
    grid[n].nyc = grid[n].jec - grid[n].jsc + 1;

    grid[n].lonc  = (double *) malloc((grid[n].nxc+1)*(grid[n].nyc+1)*sizeof(double));
    grid[n].latc  = (double *) malloc((grid[n].nxc+1)*(grid[n].nyc+1)*sizeof(double));
    
    grid[n].lont1D = (double *) malloc(nx[n]*sizeof(double));
    grid[n].latt1D = (double *) malloc(ny[n]*sizeof(double));
    grid[n].lonc1D = (double *) malloc((nx[n]+1)*sizeof(double));
    grid[n].latc1D = (double *) malloc((ny[n]+1)*sizeof(double));
    x = (double *) malloc((2*nx[n]+1)*(2*ny[n]+1)*sizeof(double));
    y = (double *) malloc((2*nx[n]+1)*(2*ny[n]+1)*sizeof(double));
    vid = mpp_get_varid(g_fid, "x");
    mpp_get_var_value(g_fid, vid, x);
    vid = mpp_get_varid(g_fid, "y");
    mpp_get_var_value(g_fid, vid, y);
    for(j=0; j<=grid[n].nyc; j++) for(i=0; i<=grid[n].nxc; i++) {
      jj = 2*(j + grid[n].jsc);
      ii = 2*(i + grid[n].isc);
      grid[n].lonc[j*(grid[n].nxc+1)+i] = x[jj*(2*nx[n]+1)+ii] * D2R;
      grid[n].latc[j*(grid[n].nxc+1)+i] = y[jj*(2*nx[n]+1)+ii] * D2R;
    }
    
    if(grid[n].is_tripolar)
      ind = nx[n]/4;
    else
      ind = 0;
    for(i=0; i<nx[n]+1; i++) grid[n].lonc1D[i] = x[2*i] * D2R;
    for(j=0; j<ny[n]+1; j++) grid[n].latc1D[j] = y[2*j*(2*nx[n]+1)+2*ind] * D2R;
    if(!grid[n].is_tripolar) ind = 1;
    
    for(i=0; i<nx[n]; i++) grid[n].lont1D[i] = x[2*nx[n]+1+2*i+1] * D2R;
    for(j=0; j<ny[n]; j++) grid[n].latt1D[j] = y[(2*j+1)*(2*nx[n]+1) + 2*ind] * D2R;
    free(x);
    free(y);
      
    /* if vector, need to get rotation angle */
    /* we assume the grid is orthogonal */
    if(opcode & VECTOR) {
      if(opcode & AGRID) {
	double *angle; 
	angle          = (double *) malloc((2*nx[n]+1)*(2*ny[n]+1)*sizeof(double));
	grid[n].cosrot = (double *) malloc(nx[n]*ny[n]*sizeof(double));	 
	grid[n].sinrot = (double *) malloc(nx[n]*ny[n]*sizeof(double));
	vid = mpp_get_varid(g_fid, "angle_dx");
	mpp_get_var_value(g_fid, vid, angle);
	grid[n].rotate = 0;
	for(j=0; j<grid[n].nyc; j++) for(i=0; i<grid[n].nxc; i++) {
	  jj = 2*(j + grid[n].jsc);
	  ii = 2*(i + grid[n].isc);
	  grid[n].cosrot[j*grid[n].nxc+i] = cos(angle[jj*(2*nx[n]+1)+ii]*D2R);
	  grid[n].sinrot[j*grid[n].nxc+i] = sin(angle[jj*(2*nx[n]+1)+ii]*D2R);
	  if(fabs(grid[n].sinrot[j*grid[n].nxc+i]) > EPSLN10) grid[n].rotate = 1;
	}
	free(angle);
      }
    }
    mpp_close(g_fid);
  }

  mpp_close(m_fid);
  
  free(nx);
  free(ny);
}; /* get_output_grid_from_mosaic*/ 

/*******************************************************************************
  void get_output_grid_by_size(Mosaic_config *mosaic, int nlon, int nlat, int finer_steps, unsigned int opcode)
  calculate output grid based on nlon, nlat and finer steps.

*******************************************************************************/
void get_output_grid_by_size(int ntiles, Grid_config *grid, double lonbegin, double lonend, double latbegin, double latend,
			     int nlon, int nlat, int finer_steps, int center_y, unsigned int opcode)
{
  double      dlon, dlat, lon_fine, lat_fine, lon_range, lat_range;
  int         nx_fine, ny_fine, i, j, layout[2];
  int nxc, nyc, ii, jj;
  
  if(ntiles !=1) mpp_error("fregrid_utils: ntiles of output mosaic must be 1 for bilinear interpolation");
  if(finer_steps && !(opcode&BILINEAR)) mpp_error("fregrid_util: finer_steps must be 0 when interp_method is not bilinear");
  
  grid->nx      = nlon;
  grid->ny      = nlat;
  grid->nx_fine = pow(2,finer_steps)*nlon;
  grid->ny_fine = pow(2,finer_steps)*(nlat-1)+1;
  nx_fine       = grid->nx_fine;
  ny_fine       = grid->ny_fine;
  lon_range     = lonend - lonbegin;
  lat_range     = latend - latbegin;
  grid->is_tripolar = 0;
  grid->lont1D = (double *)malloc(nlon*sizeof(double));
  grid->latt1D = (double *)malloc(nlat*sizeof(double));
  grid->lonc1D = (double *)malloc((nlon+1)*sizeof(double));
  grid->latc1D = (double *)malloc((nlat+1)*sizeof(double));

  dlon=lon_range/nlon;
  for(i=0; i<nlon; i++) grid->lont1D[i]  = (lonbegin + (i + 0.5)*dlon)*D2R;
  for(i=0; i<=nlon; i++) grid->lonc1D[i] = (lonbegin + i*dlon)*D2R;

  layout[0] = 1;
  layout[1] = mpp_npes();
  mpp_define_domain2d(grid->nx, grid->ny, layout, 0, 0, &(grid->domain));
  mpp_get_compute_domain2d(grid->domain, &(grid->isc), &(grid->iec), &(grid->jsc), &(grid->jec));
  grid->nxc = grid->iec - grid->isc + 1;
  grid->nyc = grid->jec - grid->jsc + 1;
  nxc       = grid->nxc;
  nyc       = grid->nyc;
  if(center_y) {
    dlat=lat_range/nlat;
    for(j=0; j<nlat; j++) grid->latt1D[j] = (latbegin+(j+0.5)*dlat)*D2R;
    for(j=0; j<=nlat; j++) grid->latc1D[j] = (latbegin+j*dlat)*D2R;
  }
  else {
    dlat=lat_range/(nlat-1);
    for(j=0; j<nlat; j++) grid->latt1D[j] = (latbegin+j*dlat)*D2R;
    for(j=0; j<=nlat; j++) grid->latc1D[j] = (latbegin+(j-0.5)*dlat)*D2R;
  }
  
  if(opcode & BILINEAR) {
    grid->latt1D_fine = (double *)malloc(ny_fine*sizeof(double));
    grid->lont   = (double *)malloc(nx_fine*ny_fine*sizeof(double));
    grid->latt   = (double *)malloc(nx_fine*ny_fine*sizeof(double));
    grid->xt     = (double *)malloc(nx_fine*ny_fine*sizeof(double));
    grid->yt     = (double *)malloc(nx_fine*ny_fine*sizeof(double));  
    grid->zt     = (double *)malloc(nx_fine*ny_fine*sizeof(double));  
    grid->vlon_t = (double *)malloc(3*nx_fine*ny_fine*sizeof(double));
    grid->vlat_t = (double *)malloc(3*nx_fine*ny_fine*sizeof(double));  

    dlon = lon_range/nx_fine;
    for(i=0; i<nx_fine; i++) {
      lon_fine = (lonbegin + (i + 0.5)*dlon)*D2R;
      for(j=0; j<ny_fine; j++) grid->lont[j*nx_fine+i] = lon_fine;
    }
    if(center_y) {
      dlat=lat_range/ny_fine;
      for(j=0; j<ny_fine; j++)grid->latt1D_fine[j] = (latbegin+(j+0.5)*dlat)*D2R;
    }
    else {
      dlat = lat_range/(ny_fine-1);
      for(j=0; j<ny_fine; j++) grid->latt1D_fine[j] = (latbegin+j*dlat)*D2R;
     
    }
    for(j=0; j<ny_fine; j++) for(i=0; i<nx_fine; i++) {
      grid->latt[j*nx_fine+i] = grid->latt1D_fine[j];
    }
    /* get the cartesian coordinates */  
    latlon2xyz(nx_fine*ny_fine, grid->lont, grid->latt, grid->xt, grid->yt, grid->zt);

    unit_vect_latlon(nx_fine*ny_fine, grid->lont, grid->latt, grid->vlon_t, grid->vlat_t);
  
  }

    grid->lonc  = (double *) malloc((nxc+1)*(nyc+1)*sizeof(double));
    grid->latc  = (double *) malloc((nxc+1)*(nyc+1)*sizeof(double));
    for(j=0; j<=nyc; j++) {
      jj = j + grid->jsc;
      for(i=0; i<=nxc; i++) {
	ii = i + grid->isc;
	grid->lonc[j*(nxc+1)+i] = grid->lonc1D[ii];
	grid->latc[j*(nxc+1)+i] = grid->latc1D[jj];
      }
    }
    if(opcode & VECTOR) { /* no rotation is needed for regular lat-lon grid. */
      grid->rotate = 0;
    }
  
}; /* get_output_grid_by_size */


/*******************************************************************************
void init_var_config(Var_config *var)
*******************************************************************************/
void init_var_config(Var_config *var, int interp_method)
{
  var->nz            = 1;
  var->nn            = 1;
  var->has_naxis     = 0;
  var->has_zaxis     = 0;
  var->has_taxis     = 0;  
  var->kstart        = 0;
  var->kend          = 0;
  var->lstart        = 0;
  var->lend          = 0;
  var->ndim          = 0;
  var->interp_method = interp_method;
  
}

/*******************************************************************************
void copy_var_config(Var_config *var)
*******************************************************************************/
void copy_var_config(const Var_config *var_in, Var_config *var_out)
{
  int i;
  
  var_out->nz            = var_in->nz;
  var_out->nn            = var_in->nn;
  var_out->has_naxis     = var_in->has_naxis;
  var_out->has_zaxis     = var_in->has_zaxis;
  var_out->has_taxis     = var_in->has_taxis;
  var_out->kstart        = var_in->kstart;
  var_out->kend          = var_in->kend ;
  var_out->lstart        = var_in->lstart;
  var_out->lend          = var_in->lend;
  var_out->ndim          = var_in->ndim;
  var_out->interp_method = var_in->interp_method;
  for(i=0; i<var_in->ndim; i++) var_out->index[i] = var_in->index[i];
}

/* We assume all the tiles have the same vgrid */
void get_output_vgrid( VGrid_config *vgrid, const char *vgrid_file )
{
  int fid, nz, vid, k;
  double *z=NULL;
  
  /* first get number of levels */
  fid = mpp_open(vgrid_file, MPP_READ);
  nz = mpp_get_dimlen(fid, "nzv");
  if((nz-1)%2) mpp_error("fregrid_util: size of dimension nzv should be 2*nlev+1");
  z = (double *)malloc(nz*sizeof(double));
  vid = mpp_get_varid(fid, "zeta");
  mpp_get_var_value(fid, vid, z);
  mpp_close(fid);
  
  nz = (nz-1)/2;
  vgrid->nz = nz;
  vgrid->z = (double *)malloc(nz*sizeof(double));
  vgrid->zb = (double *)malloc((nz+1)*sizeof(double));
  for(k=0; k<nz; k++) vgrid->z[k] = z[2*k+1];
  for(k=0; k<=nz; k++) vgrid->zb[k] = z[2*k];
  free(z);
}

void get_input_vgrid( VGrid_config *vgrid, const char *vgrid_file, const char *field )
{
  int fid, vid, vid2, ndim, i, nz;
  char dimname[32];
  char cart;
  
  /* first get number of levels */
  fid = mpp_open(vgrid_file, MPP_READ);
  
  vid = mpp_get_varid(fid, field);
  ndim = mpp_get_var_ndim(fid, vid);
  nz = 0;
  for(i=0; i<ndim; i++) {
    mpp_get_var_dimname(fid, vid, i, dimname);
    vid2 = mpp_get_varid(fid, dimname);
    cart = mpp_get_var_cart(fid, vid2);
    if(cart == 'Z') {
      nz = mpp_get_dimlen(fid, dimname);
      vgrid->nz = nz;
      vgrid->z = (double *)malloc(nz*sizeof(double));
      mpp_get_var_value(fid, vid2, vgrid->z);
    }
  }

  mpp_close(fid);
  
  if(nz == 0) mpp_error("fregrid_util: no vertical levels found in the input file");
  
}

void setup_vertical_interp(VGrid_config *vgrid_in, VGrid_config *vgrid_out)
{
  int nk1, nk2, kstart, kend, k;
  
  nk1 = vgrid_in->nz;
  nk2 = vgrid_out->nz;


  for(kstart=0; kstart<nk2; kstart++) {
    if(vgrid_out->z[kstart] >= vgrid_in->z[0]) break;
  }
  for(kend=nk2-1; kend>=0; kend--) {
    if(vgrid_out->z[kend] <= vgrid_in->z[nk1-1]) break;
  }


  if(kstart >0 && mpp_pe()==mpp_root_pe()) {
    printf("NOTE from fregrid_util: the value from level 0 to level %d will be set to the value at the shallowest source levle.\n", kstart-1);
  }
  if(kend <nk2-1 && mpp_pe()==mpp_root_pe()) {
    printf("NOTE from fregrid_util: the value from level %d to level %d will be set to the value at the deepest source levle.\n", kend+1, nk2-1);
  }  
  vgrid_out->kstart = kstart;
  vgrid_out->kend   = kend;
  vgrid_out->need_interp = 1;
  if(nk1 == nk2 ){
    for(k=0; k<nk1; k++) {
      if(fabs(vgrid_out->z[k]-vgrid_in->z[k]) > EPSLN10 ) break;
    }
    if(k==nk1) vgrid_out->need_interp = 0;
  }
}

void do_vertical_interp(VGrid_config *vgrid_in, VGrid_config *vgrid_out, Grid_config *grid_out, Field_config *field, int varid)
{
  int nk1, nk2, nx, ny, kstart, kend, i, k;
  double *tmp;
  
  if(vgrid_out->need_interp && field->var[varid].has_zaxis ) {
    nk1 = vgrid_in->nz;
    nk2 = vgrid_out->nz;
    nx  = grid_out->nx;
    ny  = grid_out->ny;
    tmp = (double *)malloc(nx*ny*nk1*sizeof(double));
    for(i=0; i<nx*ny*nk1; i++) tmp[i] = field->data[i];
    if(nk1 != nk2 ) {
      free(field->data);
      field->data =  (double *)malloc(nx*ny*nk2*sizeof(double));
    }
    
    kstart = vgrid_out->kstart;
    kend   = vgrid_out->kend;
    for(k=0; k<kstart; k++) {
      for(i=0; i<nx*ny; i++) field->data[k*nx*ny+i] = tmp[i];
    }
    for(k=kend; k<nk2; k++) {
      for(i=0; i<nx*ny; i++) field->data[k*nx*ny+i] = tmp[(nk1-1)*nx*ny+i];
    }
    nk2 = kend - kstart + 1;
    linear_vertical_interp(nx, ny, nk1, nk2, vgrid_in->z, vgrid_out->z+kstart, tmp, field->data+kstart*nx*ny); 
    free(tmp);
  }
  
}

/*******************************************************************************
  void get_input_metadata(Mosaic_config, *mosaic)
*******************************************************************************/
void get_input_metadata(int ntiles, int nfiles, File_config *file1, File_config *file2,
		        Field_config *scalar, Field_config *u_comp, Field_config *v_comp,
			const Grid_config *grid, int kbegin, int kend, int lbegin, int lend,
			unsigned int opcode, char *associated_file_dir)
{
  int     n, m, i, l, ll, nscalar, nvector, nfield;
  int     ndim, dimsize[5], nz;
  nc_type type[5];
  char    cart[5];
  char    dimname[5][STRING], bndname[5][STRING], errmsg[STRING];
  File_config  *file  = NULL;
  Field_config *field = NULL;
  size_t start[4], nread[4];
  int interp_method, use_bilinear, use_conserve;
  int len, found;
  
  
  use_bilinear = 0;
  use_conserve = 0;
  if(opcode & CONSERVE_ORDER1) {
    use_conserve = 1;
    interp_method = CONSERVE_ORDER1;
  }
  else if(opcode & CONSERVE_ORDER2) {
    use_conserve = 1;
    interp_method = CONSERVE_ORDER2;
  }
  else if(opcode & BILINEAR) {
    use_bilinear = 1;
    interp_method = BILINEAR;
  }
    
  /* First find out how many fields in file and file2. */
  nscalar = 0;
  nvector = 0;
  if( scalar) nscalar = scalar->nvar;
  if( u_comp) nvector = u_comp->nvar;

  for(n=0; n<4; n++) {
    start[n] = 0; nread[n] = 1;
  }
  
  for(m=0; m<nfiles; m++) {
    file = m==0? file1:file2;  
    for(n=0; n<ntiles; n++) {
      file[n].nt        = 1;
      file[n].axis      = (Axis_config *)malloc(MAXDIM*sizeof(Axis_config));
      file[n].ndim      = 0;
      file[n].has_tavg_info = 0;
    }
  }
             
  for(l=0; l<nscalar; l++) for(n=0; n<ntiles; n++) {
    init_var_config(scalar[n].var+l, interp_method);
  }
    
  for(l=0; l<nvector; l++) for(n=0; n<ntiles; n++) {
    init_var_config(u_comp[n].var+l, interp_method);
    init_var_config(v_comp[n].var+l, interp_method);
  }
  
  nfield = (nfiles == 1)? (nscalar+2*nvector):nvector;  /* when nfiles = 2, no scalar */
  
  for(m=0; m<nfiles; m++) {
    file = m==0? file1:file2;
    for(n=0; n<ntiles; n++) {
      const int MAXLIST = 10;
      char file_list[MAXLIST][STRING];
      char file_id[MAXLIST];
      int nfile_list;

      for(i=0; i<MAXLIST; i++) {
	file_list[i][0] = '\0';
	file_id[i] = 0;
      }
      nfile_list = 0;
      
      file[n].has_cell_measure_att = 0;
      for(l=0; l<nfield; l++) {
	if(nfiles == 1) {
	  if(l<nscalar) {
	    field = scalar;
	    ll = l;
	  }
	  else if(l<nscalar + nvector) {
	    field = u_comp;
	    ll = l - nscalar;
	  }
	  else {
	    field = v_comp;
	    ll = l - nscalar - nvector;
	  }
	}
	else {
	  ll = l;
	  field = (m==0)?u_comp:v_comp;
	}

      	field[n].var[ll].vid = mpp_get_varid(file[n].fid, field[n].var[ll].name);
	field[n].var[ll].type  = mpp_get_var_type(file[n].fid, field[n].var[ll].vid);
	if(field[n].var[ll].type != NC_SHORT && field[n].var[ll].type != NC_INT &&
	   field[n].var[ll].type != NC_FLOAT && field[n].var[ll].type != NC_DOUBLE ) {
	  sprintf(errmsg, "fregrid_util(get_input_metadata): field %s in file %s has an invalid type, "
		  "the type should be NC_DOUBLE, NC_FLOAT, NC_INT or NC_SHORT", field[n].var[ll].name, file[n].name );
	  mpp_error(errmsg);
	}
        /* check if time_avg_info attribute existed in any variables */
	if(!file[n].has_tavg_info) {
	  if(mpp_var_att_exist(file[n].fid, field[n].var[ll].vid, "time_avg_info") ) file[n].has_tavg_info = 1;
	}

	/* check the cell_methods */
        field[n].var[ll].cell_methods=CELL_METHODS_MEAN;	
	if(mpp_var_att_exist(file[n].fid, field[n].var[ll].vid, "cell_methods")) {
	  char attval[STRING] = "";
          char areaval[STRING] = "";
          int status;
          char errout[STRING], errmsg[STRING];
          mpp_get_var_att(file[n].fid, field[n].var[ll].vid, "cell_methods", attval);
          status = parse_string(attval, "area:", areaval, errout);
          if(status==-1) {
              sprintf(errmsg, "fregrid_util(get_input_metadata): %s for cell_methods "
                      "attribute of field %s in file %s", errout, field[n].var[ll].name, file[n].name );
              mpp_error(errmsg);
          }
          else if(status ==1) {
	    if(strcmp(areaval,"mean")==0)
	      field[n].var[ll].cell_methods=CELL_METHODS_MEAN;
	    else if(strcmp(areaval,"sum")==0)
	      field[n].var[ll].cell_methods=CELL_METHODS_SUM;
    	    else {
	      sprintf(errmsg, "fregrid_util(get_input_metadata): field %s in file %s attribute cell_methods "
	    	       "should have value 'mean' or 'sum' after area: ", field[n].var[ll].name, file[n].name );
	      mpp_error(errmsg);
	    }
          }
	}
	/* check if exist attribute "cell_measures" and get the name of area */
        if(mpp_var_att_exist(file[n].fid, field[n].var[ll].vid, "cell_measures")) {
          char attval[STRING] = "";
          char *str2=NULL;
	  char errout[STRING];
          int status;
 
	  file[n].has_cell_measure_att = 1;  
          /* --source_grid must be set when cell_measures attribute exists */
	  mpp_get_var_att(file[n].fid, field[n].var[ll].vid, "cell_measures", attval);
          status = parse_string(attval, "area:", field[n].var[ll].area_name, errout);
	  if(status==-1) {
	      sprintf(errmsg, "fregrid_util(get_input_metadata): %s for cell_measure "
		      "attribute of field %s in file %s", errout, field[n].var[ll].name, file[n].name );
	      mpp_error(errmsg);
	  }
          else if(status ==1)
	    field[n].var[ll].cell_measures=1;
	  else
	    field[n].var[ll].cell_measures=0;
	  if(field[n].var[ll].cell_measures==1) {
            char associated_file[512], dimname[32];
	    int vid2;
	    char cart;
	    /* make sure field[n].var[ll].area_name exist in the current file or associated file */
	    if( mpp_var_exist(file[n].fid, field[n].var[ll].area_name) ) {
              strcpy(associated_file, file[n].name );
              field[n].var[ll].area_fid = file[n].fid;
            }
            else {
	      char globalatt[1024], file1[512], str1[STRING], name[STRING], file2[STRING];

              /* check if the variable is in associated_files or not */
              if( !mpp_global_att_exist(file[n].fid, "associated_files") ) {
                 sprintf(errmsg, "fregrid_util(get_input_metadata):  field %s does not exist in file %s, "
                                 "and the file also does not have global attribute 'associated_files' ",
                                  field[n].var[ll].area_name, file[n].name );
                 mpp_error(errmsg);
              }           
              mpp_get_global_att(file[n].fid, "associated_files", globalatt);
	      sprintf(name, "%s:", field[n].var[ll].area_name);
              status = parse_string(globalatt, name, file2, errout);
	      if(status==0) {
		sprintf(errmsg, "fregrid_util(get_input_metadata): global sttribute associated_files "
			"does not contains string %s in file %s",
			name, file[n].name );
		mpp_error(errmsg);
	      }
	      else if(status==-1) {
		sprintf(errmsg, "fregrid_util(get_input_metadata): %s for associated_files "
			"global attribute in file %s", errout, file[n].name );
		mpp_error(errmsg);
	      }
	      if(associated_file_dir)
		sprintf(file1, "%s/%s", associated_file_dir, file2);
	      else
		strcpy(file1, file2);
	      
              /* check if the file exist or not, if not add tile# */
              if(mpp_file_exist(file1)) 
                 strcpy(associated_file, file1);
              else {    /* add tile number if there is more than one tile */
                 if(ntiles==1) {
                    sprintf(errmsg, "fregrid_util(get_input_metadata): ntiles==1 and file %s does not exist", file1);
                    mpp_error(errmsg);
                 }

		 len = strlen(file1);
		 if( strcmp(file1+len-3, ".nc") ==0 ) {
		    strncpy(str1, file1, len-3);
		    str1[len-3] = '\0';
		 }
		 else
		    strcpy(str1, file1);
	         sprintf(associated_file, "%s.tile%d.nc", str1, n+1);
	      
                 if( ! mpp_file_exist(associated_file) ) {
                    sprintf(errmsg, "fregrid_util(get_input_metadata): both %s and %s do not exist", file1, associated_file);
                    mpp_error(errmsg);
                 }
              }
              {
		/*search through the list */
		int found_file;
		
		found_file = 0;
		for(i=0; i<nfile_list; i++) {
		  if( !strcmp(associated_file, file_list[i]) ) {
		    field[n].var[ll].area_fid = file_id[i];
		    found_file = 1;
		    break;
		  }
		}
		if(!found_file) {
                  field[n].var[ll].area_fid = mpp_open(associated_file, MPP_READ);
		  if(nfile_list >= MAXLIST)mpp_error("fregrid_util(get_input_metadata): nfile_list > MAXLIST");
		  strcpy(file_list[i], associated_file);
		  file_id[i] = field[n].var[ll].area_fid;
		  nfile_list++;
		}
              }
	    }
	    /* get the vid of area_name and check if area_name is time dependent or not */
	    field[n].var[ll].area_vid = mpp_get_varid(field[n].var[ll].area_fid, field[n].var[ll].area_name);
	    ndim = mpp_get_var_ndim(field[n].var[ll].area_fid, field[n].var[ll].area_vid);
	    /* check if it has T-axis */
	    mpp_get_var_dimname(field[n].var[ll].area_fid, field[n].var[ll].area_vid, 0, dimname);
	    vid2 = mpp_get_varid(field[n].var[ll].area_fid, dimname);
	    cart = mpp_get_var_cart(field[n].var[ll].area_fid, vid2);
    
	    if( cart == 'T' ) {
	      if(ndim <=2) {
		sprintf(errmsg, "fregrid_util(get_input_metadata):  number of dimension of field %s in file %s has t-axis and <3" ,
		      field[n].var[ll].area_name, associated_file );
		mpp_error(errmsg);
	      }
	      field[n].var[ll].area_has_taxis = 1;
	    }
	    else {
	      field[n].var[ll].area_has_taxis = 0;
	    }
	    /* check if has n-axis (diurnal data). */
	    field[n].var[ll].area_has_naxis = 0;
	    if(ndim>2) {
	      if(field[n].var[ll].area_has_taxis)
		mpp_get_var_dimname(field[n].var[ll].area_fid, field[n].var[ll].area_vid, 1, dimname);
	      else
		mpp_get_var_dimname(field[n].var[ll].area_fid, field[n].var[ll].area_vid, 0, dimname);
	      vid2 = mpp_get_varid(field[n].var[ll].area_fid, dimname);
	      cart = mpp_get_var_cart(field[n].var[ll].area_fid, vid2);
	      if(cart == 'N') {
		field[n].var[ll].area_has_naxis = 1;
		if(ndim==3 && field[n].var[ll].area_has_taxis) {
		  sprintf(errmsg, "fregrid_util(get_input_metadata): ndim=3, has_taxis=T and hax_naxis=T for field %s in file %s",
			  field[n].var[ll].area_name, associated_file );
		}
		else if(ndim==4 && !field[n].var[ll].area_has_taxis) {
		  sprintf(errmsg, "fregrid_util(get_input_metadata): ndim=4, has_taxis=F for field %s in file %s",
			  field[n].var[ll].area_name, associated_file );
		}  
	      }
            }

	    if(ndim>4) {
	      sprintf(errmsg, "fregrid_util(get_input_metadata):  number of dimension of field %s in file %s > 4" ,
		      field[n].var[ll].area_name, associated_file );
	      mpp_error(errmsg);
	    }
	    if( mpp_var_att_exist(field[n].var[ll].area_fid, field[n].var[ll].area_vid, "missing_value") ) {
	      mpp_get_var_att_double(field[n].var[ll].area_fid, field[n].var[ll].area_vid, "missing_value", &(field[n].var[ll].area_missing));
	    }
	    else if( mpp_var_att_exist(field[n].var[ll].area_fid, field[n].var[ll].area_vid, "_FillValue") ) {
	      mpp_get_var_att_double(field[n].var[ll].area_fid, field[n].var[ll].area_vid, "_FillValue", &(field[n].var[ll].area_missing));
	    }
	    else
	      field[n].var[ll].area_missing = 0;	
	  }
        }
	else {
	    field[n].var[ll].cell_measures=0;
	    field[n].var[ll].area_missing = 0;
	    field[n].var[ll].area_has_taxis = 0;
	}
	 
	/* get the interp_method from the field attribute if existing 
           when interp_method is not conserve_order2_monotonic
        */
        if( !(opcode & MONOTONIC) ) {
	  if(mpp_var_att_exist(file[n].fid, field[n].var[ll].vid, "interp_method")) {
            char    remap_method[STRING] = "";
	    mpp_get_var_att(file[n].fid, field[n].var[ll].vid, "interp_method", remap_method);
            if(!strcmp(remap_method, "conserve_order1") ) {
  	      use_conserve = 1;
  	      field[n].var[ll].interp_method = CONSERVE_ORDER1;
	    }
  	    else if(!strcmp(remap_method, "conserve_order2") ) {
	      use_conserve = 1;
	      field[n].var[ll].interp_method = CONSERVE_ORDER2;
	    }
	    else if(!strcmp(remap_method, "bilinear") ) {
	      use_bilinear = 1;
	      field[n].var[ll].interp_method = BILINEAR;
	    }
	    else {
	      sprintf(errmsg, "get_input_metadata(fregrid_util.c): in file %s, attribute interp_method of field %s has value = %s"
 		      "is not suitable, it should be conserve_order1, conserve_order2 or bilinear", file[n].name,
		      field[n].var[ll].name, remap_method);
      	      mpp_error(errmsg);
	    }
	  }
	}
	ndim = mpp_get_var_ndim(file[n].fid, field[n].var[ll].vid);
	if(ndim <2 || ndim>5) mpp_error("get_input_metadata(fregrid_util.c): ndim should be no less than 2 and no larger than 5");	
	for(i=0; i<ndim; i++) {
	  int vid;
          mpp_get_var_dimname(file[n].fid, field[n].var[ll].vid, i, dimname[i]);
	  dimsize[i] = mpp_get_dimlen(file[n].fid, dimname[i]);
	  vid = mpp_get_varid(file[n].fid, dimname[i]);
	  cart[i] = mpp_get_var_cart(file[n].fid, vid);
	  type[i] = mpp_get_var_type(file[n].fid, vid);
	  mpp_get_var_bndname(file[n].fid, vid, bndname[i]);
      	}
	field[n].var[ll].ndim = ndim;
	if(cart[ndim-1] != 'X') mpp_error("get_input_metadata(fregrid_util.c): the last dimension cartesian should be 'X'");
	if(cart[ndim-2] != 'Y') mpp_error("get_input_metadata(fregrid_util.c): the second last dimension cartesian should be 'Y'");
	if(dimsize[ndim-1] != grid[n].nx) mpp_error("get_input_metadata(fregrid_util.c): x-size in grid file in not the same as in data file");
	if(dimsize[ndim-2] != grid[n].ny) mpp_error("get_input_metadata(fregrid_util.c): y-size in grid file in not the same as in data file");
        if(ndim > 2) {
	  if(cart[ndim-3] == 'Z') {
	    field[n].var[ll].has_zaxis = 1;
	    field[n].var[ll].nz        = dimsize[ndim-3];
	    if(kend > field[n].var[ll].nz) {
	      sprintf(errmsg, "get_input_metadata(fregrid_util.c): KlevelEnd should be no larger than "
		      "number of vertical levels of field %s in file %s.", field[n].var[ll].name, file[n].name);
	      mpp_error(errmsg);
	    }
	    if(kbegin>0) {
	      field[n].var[ll].kstart = kbegin - 1;
	      field[n].var[ll].kend   = kend - 1;
	      field[n].var[ll].nz     = kend - kbegin + 1;
	    }
	    else {
	      field[n].var[ll].kstart = 0;
	      field[n].var[ll].kend   = field[n].var[ll].nz - 1;
	    }	    
	  }
	  else if(cart[ndim-3] == 'N') {
	    field[n].var[ll].has_naxis = 1;
	    field[n].var[ll].nn        = dimsize[ndim-3];
	  }
	}
	if(ndim > 3) {
	  if(cart[ndim-4] == 'Z') {
	    mpp_error("get_input_metadata(fregrid_util.c): the Z-axis must be the third dimension");
	  }
	  if(cart[ndim-4] == 'N') {
	    field[n].var[ll].has_naxis = 1;
	    field[n].var[ll].nn        = dimsize[ndim-4];
	  }
	}
	  
	if(cart[0] == 'T') {
	  field[n].var[ll].has_taxis = 1;
	  if(lend > dimsize[0]) {
	    sprintf(errmsg, "get_input_metadata(fregrid_util.c): LstepEnd should be no larger than "
		    "number of time levels of field %s in file %s.", field[n].var[ll].name, file[n].name);
	    mpp_error(errmsg);
	  }
	  if(lbegin>0) {
	    field[n].var[ll].lstart = lbegin - 1;
	    field[n].var[ll].lend   = lend - 1;
	    file[n].nt              = lend - lbegin + 1;
	  }
	  else {
	    field[n].var[ll].lstart = 0;
	    field[n].var[ll].lend   = dimsize[0] - 1;
	    file[n].nt              = dimsize[0];
	  }
	    
	}
	for(i=0; i<ndim; i++) {
	  /* loop through all the file dimensions to see if the dimension already exist or not */
	  int found, j; 
	  found = 0;
	  for(j=0; j<file[n].ndim; j++) {
	    if(!strcmp(dimname[i], file[n].axis[j].name) ) {
	      found = 1;
	      field[n].var[ll].index[i] = j;
	      break;
	    }
	  }
	  if(!found) {
	    j                         = file[n].ndim;
	    field[n].var[ll].index[i] = file[n].ndim;
	    file[n].ndim++;
            if(	file[n].ndim > MAXDIM) mpp_error("get_input_metadata(fregrid_util.c):ndim is greater than MAXDIM");

	    file[n].axis[j].cart = cart[i];
	    file[n].axis[j].type = type[i];
	    strcpy(file[n].axis[j].name, dimname[i]);
	    strcpy(file[n].axis[j].bndname, bndname[i]);
	    file[n].axis[j].vid = mpp_get_varid(file[n].fid, dimname[i]);
	    if(cart[i] == 'T') {
	      start[0] = field[n].var[ll].lstart;
	      file[n].axis[j].size = file[n].nt;
	    }
	    else if(cart[i] == 'Z') {
	      start[0] = field[n].var[ll].kstart;
	      file[n].axis[j].size = field[n].var[ll].nz;
	    }
	    else {
	      start[0] = 0;
               file[n].axis[j].size = dimsize[i];
	    }
	    file[n].axis[j].data = (double *)malloc(file[n].axis[j].size*sizeof(double));
	    nread[0] = file[n].axis[j].size;
	    mpp_get_var_value_block(file[n].fid, file[n].axis[j].vid, start, nread, file[n].axis[j].data);
	    file[n].axis[j].bndtype = 0;
	    if(strcmp(bndname[i], "none") ) {
	      file[n].axis[j].bndid = mpp_get_varid(file[n].fid, bndname[i]);
	      if(mpp_get_var_ndim(file[n].fid,file[n].axis[j].bndid) == 1) {
		file[n].axis[j].bndtype = 1;
		file[n].axis[j].bnddata = (double *)malloc((file[n].axis[j].size+1)*sizeof(double));
		nread[0] = file[n].axis[j].size+1;
	      }
	      else {
	        file[n].axis[j].bndtype = 2;
		file[n].axis[j].bnddata = (double *)malloc(2*file[n].axis[j].size*sizeof(double));
		nread[0] = file[n].axis[j].size; nread[1] = 2;
	      }
	      mpp_get_var_value_block(file[n].fid, file[n].axis[j].bndid, start, nread, file[n].axis[j].bnddata);
	    }
	    else if( cart[i] == 'X' || cart[i] == 'Y' ) {
	      sprintf(file[n].axis[j].bndname, "%s_bnds", file[n].axis[j].name);
	    }
	      
	  }
	} /*ndim*/
      }  /* nvar */
    } /* ntile */
    /* make sure the consistency between tiles */
    for(n=1; n<ntiles; n++) {
      if(file[n].has_tavg_info != file[0].has_tavg_info)
	mpp_error("get_input_metadata(fregrid_util.c): mismatch between tiles for field attribute has_tavg_info");
      
      if(file[n].ndim != file[0].ndim)
	mpp_error("get_input_metadata(fregrid_util.c): mismatch between tiles for file ndim");
      for(l=0; l<file[n].ndim; l++) {
	if(strcmp(file[n].axis[l].name, file[0].axis[l].name) )
	   mpp_error("get_input_metadata(fregrid_util.c): mismatch between tiles for file axis name");
      }
      for(l=0; l<nscalar+2*nvector; l++) {
	if(l<nscalar) {
	  field = scalar;
	  ll = l;
	}
	else if(l<nscalar + nvector) {
	  field = u_comp;
	  ll = l - nscalar;
	}
	else {
	  field = v_comp;
	  ll = l - nscalar - nvector;
	}      

        if(field[n].var[ll].cell_measures != field[0].var[ll].cell_measures) {
	    sprintf(errmsg, "get_input_metadata(fregrid_util.c): mismatch of attribute `cell_measures` between tiles "
		    "for field %s in file %s", field[n].var[ll].name, file[n].name);
	    mpp_error(errmsg);
	}
        if(field[n].var[ll].area_has_taxis != field[0].var[ll].area_has_taxis) {
	    sprintf(errmsg, "get_input_metadata(fregrid_util.c): mismatch of has_taxis between tiles "
		    "for field %s in file %s", field[n].var[ll].area_name, file[n].name);
	    mpp_error(errmsg);
	}	
	if(field[n].var[ll].ndim != field[0].var[ll].ndim)
	  mpp_error("get_input_metadata(fregrid_util.c): mismatch between tiles for var ndim");
	if(field[n].var[ll].interp_method != field[0].var[ll].interp_method)
	  mpp_error("get_input_metadata(fregrid_util.c): mismatch between tiles for interp_method");
	for(i=0; i<field[n].var[ll].ndim; i++) {  
	  if(field[n].var[ll].index[i] != field[0].var[ll].index[i])
	    mpp_error("get_input_metadata(fregrid_util.c): mismatch between tiles for var dimindex");
	}
      }
    }
    /* close the file and get the tavg_info */
    for(n=0; n<ntiles; n++) {
      if(file[n].has_tavg_info) {

	file[n].t1    = (double *)malloc(file[n].nt*sizeof(double));
	file[n].t2    = (double *)malloc(file[n].nt*sizeof(double));	
	file[n].dt    = (double *)malloc(file[n].nt*sizeof(double));
	file[n].id_t1 = mpp_get_varid(file[n].fid, "average_T1");
	file[n].id_t2 = mpp_get_varid(file[n].fid, "average_T2");
	file[n].id_dt = mpp_get_varid(file[n].fid, "average_DT");
	if(lbegin > 0) 
	  start[0] = lbegin-1;
	else
	  start[0] = 0;
	nread[0] = file[n].nt; nread[1] = 1; 
	mpp_get_var_value_block(file[n].fid, file[n].id_t1, start, nread, file[n].t1);
	mpp_get_var_value_block(file[n].fid, file[n].id_t2, start, nread, file[n].t2);	
	mpp_get_var_value_block(file[n].fid, file[n].id_dt, start, nread, file[n].dt);
      }
    }

  } /*nfile*/

  /* make sure bilinear and conservative interpolation do not co-exist. */
  if(use_bilinear && use_conserve) mpp_error("get_input_metadata(fregrid_util.c): bilinear interpolation and conservative "
					     "interpolation can not co-exist, check you option interp_method in command "
					     "line and field attribute interp_method in source file");
  
}; /* get_input_metadata */
  
/* get the string after str2 in str1 and save it into strOut
   return 1 if the string is found, return 0 if not, return -1 if error found
*/
int parse_string(const char *str1, const char *str2, char *strOut, char *errmsg)
{

  char *str=NULL;
  int len2, len, istart, attlen, i;

  len2 = strlen(str2);
  str = strstr(str1, str2);
  if( str ) { /* str2 is found */
    str = str+len2;
    len = strlen(str);
	  
    /* find the start position */
    istart = len;
    for(i=0; i<len; i++) {
      if(str[i] != ' ') {
	istart = i;
	break;
      }
    }
    if(istart == len) {
      sprintf(errmsg, "empty after %s", str2);
      return -1;
    }
	  
    /* find the length of the value */
    attlen=0;
    for(i=istart; i<len; i++) {
      if(str[i] == ' ') break;
      attlen++;
    }
	
    strncpy(strOut, str+istart, attlen);
    strOut[attlen] = '\0';
    return 1;
  }

  return 0;

}



/*******************************************************************************
void set_output_metadata ( Mosaic_config *mosaic)
*******************************************************************************/

void set_output_metadata (int ntiles_in, int nfiles, const File_config *file1_in, const File_config *file2_in,
			  const Field_config *scalar_in, const Field_config *u_in, const Field_config *v_in,
			  int ntiles_out, File_config *file1_out, File_config *file2_out, Field_config *scalar_out,
			  Field_config *u_out, Field_config *v_out, const Grid_config *grid_out, const VGrid_config *vgrid_out,
			  const char *history, const char *tagname, unsigned int opcode)
{
  int j;
  int m, n, ndim, i, l, dims[5];
  int dim_bnds, dim_time;
  int nscalar, nvector;
  const File_config *file_in = NULL;
  File_config *file_out = NULL;
  int dst_is_latlon;
  int standard_dimension;

  standard_dimension = opcode & STANDARD_DIMENSION;

  dst_is_latlon = 1;
  for(n=0; n<ntiles_out; n++) {
    for(j=0; j<= grid_out[n].nyc; j++) {
      for(i=1; i<=grid_out[n].nxc; i++) {
	if(grid_out[n].latc[j*(grid_out[n].nx+1)+i] != grid_out[n].latc[j*(grid_out[n].nx+1)]) {
	  dst_is_latlon = 0;
          goto LOCATION1;
	}
      }
    }
    for(i=0; i<= grid_out[n].nxc; i++) {
      for(j=1; j<= grid_out[n].nyc; j++) {
	if(grid_out[n].lonc[j*(grid_out[n].nx+1)+i] != grid_out[n].lonc[i]) {
	  dst_is_latlon = 0;
          goto LOCATION1;
	}
      }
    }
  }
  
  LOCATION1:  nscalar = 0;
  mpp_sum_int(1, &dst_is_latlon);
  if(dst_is_latlon == mpp_npes())
    dst_is_latlon = 1;
  else
    dst_is_latlon = 0;
  
  nvector = 0;
  if( scalar_in) nscalar = scalar_in->nvar;
  if( u_in)      nvector = u_in->nvar;
  for(n=0; n<ntiles_out; n++) {
    for(l=0; l<nscalar; l++) {
      copy_var_config(scalar_in[0].var+l, scalar_out[n].var+l);
      if( vgrid_out->nz > 0) scalar_out[n].var[l].nz = vgrid_out->nz; 
    }
    
    for(l=0; l<nvector; l++) {
      copy_var_config(u_in[0].var+l, u_out[n].var+l);
      copy_var_config(v_in[0].var+l, v_out[n].var+l);
    }
  }
  
  for(m=0; m<nfiles; m++) {
    file_in  = m==0? file1_in:file2_in;
    file_out = m==0? file1_out:file2_out;
    for(n=0; n<ntiles_out; n++) {
      char tilename[STRING];
      file_out[n].nt        = file_in[0].nt;
      ndim                  = file_in[0].ndim;
      file_out[n].ndim      = ndim;
      file_out[n].axis      = (Axis_config *)malloc(ndim*sizeof(Axis_config));
      file_out[n].has_tavg_info = file_in[0].has_tavg_info;
      
      for(i=0; i<ndim; i++) {
	file_out[n].axis[i].cart = file_in[0].axis[i].cart;
	if( file_in[0].axis[i].cart == 'Z' && vgrid_out->nz >0 )
	  file_out[n].axis[i].size = vgrid_out->nz;
	else
	  file_out[n].axis[i].size = file_in[0].axis[i].size;
	if(file_out[n].axis[i].cart == 'X') file_out[n].axis[i].size = grid_out[n].nx;
	if(file_out[n].axis[i].cart == 'Y') file_out[n].axis[i].size = grid_out[n].ny;	
	file_out[n].axis[i].type    = file_in[0].axis[i].type;
	file_out[n].axis[i].bndtype = file_in[0].axis[i].bndtype;
	if(standard_dimension &&  (file_out[n].axis[i].cart == 'X' || file_out[n].axis[i].cart == 'Y') )
	  file_out[n].axis[i].bndtype = 3;
	if(file_out[n].axis[i].bndtype ==0 && (file_out[n].axis[i].cart == 'X' || file_out[n].axis[i].cart == 'Y')
	   && dst_is_latlon) file_out[n].axis[i].bndtype = 3;
	strcpy(file_out[n].axis[i].name, file_in[0].axis[i].name);
	strcpy(file_out[n].axis[i].bndname, file_in[0].axis[i].bndname);
      }
      for(i=0; i<ndim; i++) {
	file_out[n].axis[i].data = (double *)malloc(file_out[n].axis[i].size*sizeof(double));
	if( file_out[n].axis[i].cart == 'X' ) {   /* x-axis */
	  for(l=0; l<file_out[n].axis[i].size; l++) 
	    file_out[n].axis[i].data[l] = grid_out[n].lont1D[l]*R2D; /* T-cell center */
	}
	else if ( file_out[n].axis[i].cart == 'Y') { /* y-axis */
	  for(l=0; l<file_out[n].axis[i].size; l++) 
	    file_out[n].axis[i].data[l] = grid_out[n].latt1D[l]*R2D; /* T-cell center */
	}	  
	else if( file_out[n].axis[i].cart == 'Z' && vgrid_out->nz > 0 ) { /* z-axis */
	  for(l=0; l<file_out[n].axis[i].size; l++) 
	    file_out[n].axis[i].data[l] = vgrid_out->z[l];
	}
	else {
	  for(l=0; l<file_out[n].axis[i].size; l++) 
	    file_out[n].axis[i].data[l] = file_in[0].axis[i].data[l];
	}
	switch( file_out[n].axis[i].bndtype ) {
      	case 1:
	  file_out[n].axis[i].bnddata = (double *)malloc((file_out[n].axis[i].size+1)*sizeof(double));
	  if( file_out[n].axis[i].cart == 'X' ) {   /* x-axis */
	    for(l=0; l<=file_out[n].axis[i].size; l++) file_out[n].axis[i].bnddata[l] = grid_out[n].lonc1D[l]*R2D;
	  }
	  else if(file_out[n].axis[i].cart == 'Y') {
	    for(l=0; l<=file_out[n].axis[i].size; l++) file_out[n].axis[i].bnddata[l  ] = grid_out[n].latc1D[l]*R2D;
	  }
	  else if(file_out[n].axis[i].cart == 'Z' && vgrid_out->nz > 0) {
	    for(l=0; l<=file_out[n].axis[i].size; l++) file_out[n].axis[i].bnddata[l  ] = vgrid_out->zb[l];
	  }
	  else{
	    for(l=0; l<=file_out[n].axis[i].size; l++) file_out[n].axis[i].bnddata[l] = file_in[0].axis[i].bnddata[l];
	  }
	  break;
     	case 2:
	  file_out[n].axis[i].bnddata = (double *)malloc(2*file_out[n].axis[i].size*sizeof(double));
	  if( file_out[n].axis[i].cart == 'X' ) {   /* x-axis */
	    for(l=0; l<file_out[n].axis[i].size; l++) {
	      file_out[n].axis[i].bnddata[2*l  ] = grid_out[n].lonc1D[l]*R2D;
	      file_out[n].axis[i].bnddata[2*l+1] = grid_out[n].lonc1D[l+1]*R2D;
	    }
	  }
	  else if(file_out[n].axis[i].cart == 'Y') {
	    for(l=0; l<file_out[n].axis[i].size; l++) {
	      file_out[n].axis[i].bnddata[2*l  ] = grid_out[n].latc1D[l]*R2D;
	      file_out[n].axis[i].bnddata[2*l+1] = grid_out[n].latc1D[l+1]*R2D;
	    }
	  }
	  else if(file_out[n].axis[i].cart == 'Z' && vgrid_out->nz > 0) {
	    for(l=0; l<file_out[n].axis[i].size; l++) {
	      file_out[n].axis[i].bnddata[2*l  ] = vgrid_out->zb[l];
	      file_out[n].axis[i].bnddata[2*l+1] = vgrid_out->zb[l+1];
	    }
	  }
	  else {
	    for(l=0; l<file_out[n].axis[i].size; l++) {
	      file_out[n].axis[i].bnddata[2*l  ] = file_in[0].axis[i].bnddata[2*l];
	      file_out[n].axis[i].bnddata[2*l+1] = file_in[0].axis[i].bnddata[2*l+1];
	    }
	  }
	  break;
       	case 3:
	  file_out[n].axis[i].bnddata = (double *)malloc(2*file_out[n].axis[i].size*sizeof(double));
	  if( file_out[n].axis[i].cart == 'X' ) {   /* x-axis */
	    for(l=0; l<file_out[n].axis[i].size; l++) {
	      file_out[n].axis[i].bnddata[2*l  ] = grid_out[n].lonc1D[l]*R2D;
	      file_out[n].axis[i].bnddata[2*l+1] = grid_out[n].lonc1D[l+1]*R2D;
	    }
	  }
	  else if(file_out[n].axis[i].cart == 'Y') {
	    for(l=0; l<file_out[n].axis[i].size; l++) {
	      file_out[n].axis[i].bnddata[2*l  ] = grid_out[n].latc1D[l]*R2D;
	      file_out[n].axis[i].bnddata[2*l+1] = grid_out[n].latc1D[l+1]*R2D;
	    }
	  }
	  break;	  	  
	}
	
      }
	for(l=0; l<nscalar; l++)scalar_out[n].var[l].type = scalar_in[0].var[l].type;
      
      if(mpp_pe() == mpp_root_pe()) {
	file_out[n].fid = mpp_open(file_out[n].name, MPP_WRITE);
	mpp_copy_global_att(file_in[0].fid, file_out[n].fid);
	mpp_def_global_att(file_out[n].fid, "history", history);
	mpp_def_global_att(file_out[n].fid, "code_version", tagname);
	/* define dim_bnds if bnds axis exist */
	for(i=0; i<ndim; i++) {
	  if(file_out[n].axis[i].bndtype == 2 || file_out[n].axis[i].bndtype == 3 ||
	     (file_out[n].axis[i].bndtype == 1 && standard_dimension) ) {
	    dim_bnds = mpp_def_dim(file_out[n].fid, "bnds", 2);
	    break;
	  }
	}
	
	for(i=0; i<ndim; i++) {
	  if(file_out[n].axis[i].cart=='T') {
	    file_out[n].axis[i].dimid = mpp_def_dim(file_out[n].fid, file_out[n].axis[i].name, NC_UNLIMITED);
	    dim_time = file_out[n].axis[i].dimid;
	  }
	  else {
	    if((file_out[n].axis[i].type == NC_INT || standard_dimension) && file_out[n].axis[i].cart == 'X' ) 
	      file_out[n].axis[i].dimid = mpp_def_dim(file_out[n].fid, "lon", file_out[n].axis[i].size);
	    else if((file_out[n].axis[i].type == NC_INT || standard_dimension) && file_out[n].axis[i].cart == 'Y' ) 
	      file_out[n].axis[i].dimid = mpp_def_dim(file_out[n].fid, "lat", file_out[n].axis[i].size);
	    else 
	      file_out[n].axis[i].dimid = mpp_def_dim(file_out[n].fid, file_out[n].axis[i].name, file_out[n].axis[i].size);
	  }
	}
	for(i=0; i<ndim; i++) {
	  if(file_out[n].axis[i].type == NC_INT && file_out[n].axis[i].cart == 'X' ) {
	    file_out[n].axis[i].vid = mpp_def_var(file_out[n].fid, "lon", NC_DOUBLE, 1,
						  &(file_out[n].axis[i].dimid), 0);
	    mpp_def_var_att(file_out[n].fid, file_out[n].axis[i].vid, "units", "degrees");
	    mpp_def_var_att(file_out[n].fid, file_out[n].axis[i].vid, "axis", "X");
	    mpp_def_var_att(file_out[n].fid, file_out[n].axis[i].vid, "standard_name", "grid_longitude");
	    if( file_out[n].axis[i].bndtype != 3 ) mpp_error("fregrid_util.c: axis bndtype must be 3 when axis vartype is NC_INT");
	    mpp_def_var_att(file_out[n].fid, file_out[n].axis[i].vid, "bounds", "lon_bnds");
	    /* define bounds variable */
	    dims[0] = file_out[n].axis[i].dimid;
	    dims[1] = dim_bnds;
	    file_out[n].axis[i].bndid = mpp_def_var(file_out[n].fid, "lon_bnds", NC_DOUBLE, 2, dims , 0);
	    mpp_def_var_att(file_out[n].fid, file_out[n].axis[i].bndid, "units", "degrees");
	    mpp_def_var_att(file_out[n].fid, file_out[n].axis[i].bndid, "standard_name", "grid_longitude_bounds");
	  }
	  else if(file_out[n].axis[i].type == NC_INT && file_out[n].axis[i].cart == 'Y' ) {
	    file_out[n].axis[i].vid = mpp_def_var(file_out[n].fid, "lat", NC_DOUBLE, 1,
						  &(file_out[n].axis[i].dimid), 0);
	    mpp_def_var_att(file_out[n].fid, file_out[n].axis[i].vid, "units", "degrees");
	    mpp_def_var_att(file_out[n].fid, file_out[n].axis[i].vid, "axis", "Y");
	    mpp_def_var_att(file_out[n].fid, file_out[n].axis[i].vid, "standard_name", "grid_latitude");
	    if( file_out[n].axis[i].bndtype != 3 ) mpp_error("fregrid_util.c: axis bndtype must be 3 when axis vartype is NC_INT");
	    mpp_def_var_att(file_out[n].fid, file_out[n].axis[i].vid, "bounds", "lat_bnds");
	    /* define bounds variable */
	    dims[0] = file_out[n].axis[i].dimid;
	    dims[1] = dim_bnds;
	    file_out[n].axis[i].bndid = mpp_def_var(file_out[n].fid, "lat_bnds", NC_DOUBLE, 2, dims , 0);
	    mpp_def_var_att(file_out[n].fid, file_out[n].axis[i].bndid, "units", "degrees");
	    mpp_def_var_att(file_out[n].fid, file_out[n].axis[i].bndid, "standard_name", "grid_latitude_bounds");
	  }
	  else if( standard_dimension && file_out[n].axis[i].cart == 'X' ) {
	    file_out[n].axis[i].vid = mpp_def_var(file_out[n].fid, "lon", NC_DOUBLE, 1,
						  &(file_out[n].axis[i].dimid), 4, "long_name", "longitude",
						  "units", "degrees_E", "cartesian_axis", "X", "bounds", "lon_bnds" );
	    if( file_out[n].axis[i].bndtype != 3 ) mpp_error("fregrid_util.c: x axis bndtype must be 3 when standard_dimension is true");
	    /* define bounds variable */
	    dims[0] = file_out[n].axis[i].dimid;
	    dims[1] = dim_bnds;
	    file_out[n].axis[i].bndid = mpp_def_var(file_out[n].fid, "lon_bnds", NC_DOUBLE, 2, dims , 3,
						    "long_name", "longitude bounds", "units", "degrees_E", "cartesian_axis", "X");
	  }
	  else if( standard_dimension && file_out[n].axis[i].cart == 'Y' ) {
	    file_out[n].axis[i].vid = mpp_def_var(file_out[n].fid, "lat", NC_DOUBLE, 1,
						  &(file_out[n].axis[i].dimid), 4, "long_name", "latitude",
						  "units", "degrees_N", "cartesian_axis", "Y", "bounds", "lat_bnds" );
	    if( file_out[n].axis[i].bndtype != 3 ) mpp_error("fregrid_util.c: y axis bndtype must be 3 when standard_dimension is true");
	    /* define bounds variable */
	    dims[0] = file_out[n].axis[i].dimid;
	    dims[1] = dim_bnds;
	    file_out[n].axis[i].bndid = mpp_def_var(file_out[n].fid, "lat_bnds", NC_DOUBLE, 2, dims , 3,
						    "long_name", "latitude bounds", "units", "degrees_N", "cartesian_axis", "Y");
	  }	  	  
	  else {
	    file_out[n].axis[i].vid = mpp_def_var(file_out[n].fid, file_out[n].axis[i].name, file_out[n].axis[i].type, 1,
						  &(file_out[n].axis[i].dimid), 0);
	    mpp_copy_var_att(file_in[0].fid, file_in[0].axis[i].vid, file_out[n].fid, file_out[n].axis[i].vid);
	    if(file_out[n].axis[i].bndtype == 3) mpp_def_var_att(file_out[n].fid, file_out[n].axis[i].vid, "bounds", file_out[n].axis[i].bndname);
	  
	    switch( file_out[n].axis[i].bndtype ) {
	    case 3:
	      dims[0] = file_out[n].axis[i].dimid;
	      dims[1] = dim_bnds;
	      file_out[n].axis[i].bndid = mpp_def_var(file_out[n].fid, file_out[n].axis[i].bndname, file_out[n].axis[i].type, 2, dims , 0);
	      mpp_copy_var_att(file_in[0].fid, file_in[0].axis[i].vid, file_out[n].fid, file_out[n].axis[i].bndid);
	      break;
	    case 1:
	      /* the bndname axis may already defined. */
	      for(j=0; j<ndim; j++) {
		if(strcmp(file_out[n].axis[i].bndname, file_out[n].axis[j].name)==0 ) {
		  dims[0] = file_out[n].axis[j].dimid;
		  file_out[n].axis[i].bndtype = 0; 
		  break;
		}
	      }
	      if(j==ndim){ /* not found */
		dims[0] = mpp_def_dim(file_out[n].fid, file_out[n].axis[i].bndname, file_out[n].axis[i].size+1);
		file_out[n].axis[i].bndid = mpp_def_var(file_out[n].fid, file_out[n].axis[i].bndname, file_out[n].axis[i].type, 1, dims , 0);
		mpp_copy_var_att(file_in[0].fid, file_in[0].axis[i].bndid, file_out[n].fid, file_out[n].axis[i].bndid);
	      }
	      break;
	    case 2:
	      dims[0] = file_out[n].axis[i].dimid;
	      dims[1] = dim_bnds;
	      file_out[n].axis[i].bndid = mpp_def_var(file_out[n].fid, file_out[n].axis[i].bndname, file_out[n].axis[i].type, 2, dims , 0);
	      mpp_copy_var_att(file_in[0].fid, file_in[0].axis[i].bndid, file_out[n].fid, file_out[n].axis[i].bndid);
	      break;
	    }
	  }
	}

	/* define the field meta data */
	for(l=0; l<nscalar; l++) {
	  int natts, i;
	  for(i=0; i<scalar_out[n].var[l].ndim; i++) dims[i] = file_out[n].axis[scalar_out[n].var[l].index[i]].dimid;
	  scalar_out[n].var[l].vid = mpp_def_var(file_out[n].fid, scalar_out[n].var[l].name, scalar_out[n].var[l].type,
						   scalar_out[n].var[l].ndim, dims, 0);
	  /* when standard_dimension is true and output grid is lat-lon grid, do not output attribute coordinates */
	  natts = mpp_get_var_natts(file_in[0].fid, scalar_in[0].var[l].vid);
          for(i=0; i<natts; i++) {
	    char name[256];
	    mpp_get_var_attname(file_in[0].fid, scalar_in[0].var[l].vid, i, name);
	    if( !standard_dimension || !dst_is_latlon || strcmp(name, "coordinates") )
	      mpp_copy_att_by_name(file_in[0].fid, scalar_in[0].var[l].vid, file_out[n].fid, scalar_out[n].var[l].vid,name);
	  }
		
	  if(scalar_out[n].var[l].interp_method == CONSERVE_ORDER1)
             mpp_def_var_att(file_out[n].fid, scalar_out[n].var[l].vid, "interp_method", "conserve_order1");
	  else if(scalar_out[n].var[l].interp_method == CONSERVE_ORDER2)
             mpp_def_var_att(file_out[n].fid, scalar_out[n].var[l].vid, "interp_method", "conserve_order2");
	  else if(scalar_out[n].var[l].interp_method == BILINEAR)
             mpp_def_var_att(file_out[n].fid, scalar_out[n].var[l].vid, "interp_method", "bilinear");
	}
	for(l=0; l<nvector; l++) {
	  if(m==0) {
	    u_out[n].var[l].type = u_in[0].var[l].type;
	    for(i=0; i<u_out[n].var[l].ndim; i++) dims[i] = file_out[n].axis[u_out[n].var[l].index[i]].dimid;
	    u_out[n].var[l].vid = mpp_def_var(file_out[n].fid, u_out[n].var[l].name, u_out[n].var[l].type,
			        		u_out[n].var[l].ndim, dims, 0); 
	    mpp_copy_var_att(file_in[0].fid, u_in[0].var[l].vid, file_out[n].fid, u_out[n].var[l].vid);
	    if(u_out[n].var[l].interp_method == CONSERVE_ORDER1)
	      mpp_def_var_att(file_out[n].fid, u_out[n].var[l].vid, "remapping_method", "conserve_order1");
	    else if(u_out[n].var[l].interp_method == CONSERVE_ORDER2)
	      mpp_def_var_att(file_out[n].fid, u_out[n].var[l].vid, "remapping_method", "conserve_order2");
	    else if(u_out[n].var[l].interp_method == BILINEAR)
	      mpp_def_var_att(file_out[n].fid, u_out[n].var[l].vid, "remapping_method", "bilinear");
	    
	  }
	  if(m==1 || nfiles == 1) {
	    v_out[n].var[l].type = v_in[0].var[l].type;
	    for(i=0; i<v_out[n].var[l].ndim; i++) dims[i] = file_out[n].axis[v_out[n].var[l].index[i]].dimid;
	    v_out[n].var[l].vid = mpp_def_var(file_out[n].fid, v_out[n].var[l].name, v_out[n].var[l].type,
						v_out[n].var[l].ndim, dims, 0); 
	    mpp_copy_var_att(file_in[0].fid, v_in[0].var[l].vid, file_out[n].fid, v_out[n].var[l].vid);
	    if(v_out[n].var[l].interp_method == CONSERVE_ORDER1)
	      mpp_def_var_att(file_out[n].fid, v_out[n].var[l].vid, "remapping_method", "conserve_order1");
	    else if(v_out[n].var[l].interp_method == CONSERVE_ORDER2)
	      mpp_def_var_att(file_out[n].fid, v_out[n].var[l].vid, "remapping_method", "conserve_order2");
	    else if(v_out[n].var[l].interp_method == BILINEAR)
	      mpp_def_var_att(file_out[n].fid, v_out[n].var[l].vid, "remapping_method", "bilinear");
	  }
	}
	/* define time avg info variables */
	if(file_out[n].has_tavg_info) {
	  int ll;
	  file_out[n].id_t1 = mpp_def_var(file_out[n].fid, "average_T1", NC_DOUBLE, 1, &dim_time, 0);
	  mpp_copy_var_att(file_in[0].fid, file_in[0].id_t1, file_out[n].fid, file_out[n].id_t1);
	  file_out[n].id_t2 = mpp_def_var(file_out[n].fid, "average_T2", NC_DOUBLE, 1, &dim_time, 0);
	  mpp_copy_var_att(file_in[0].fid, file_in[0].id_t2, file_out[n].fid, file_out[n].id_t2);	  
	  file_out[n].id_dt = mpp_def_var(file_out[n].fid, "average_DT", NC_DOUBLE, 1, &dim_time, 0);
	  mpp_copy_var_att(file_in[0].fid, file_in[0].id_dt, file_out[n].fid, file_out[n].id_dt);
	  file_out[n].t1 = (double *)malloc(file_out[n].nt*sizeof(double));
	  file_out[n].t2 = (double *)malloc(file_out[n].nt*sizeof(double));
	  file_out[n].dt = (double *)malloc(file_out[n].nt*sizeof(double));
	  for(ll=0; ll<file_out[n].nt; ll++) {
	    file_out[n].t1[ll] = file_in[0].t1[ll];
	    file_out[n].t2[ll] = file_in[0].t2[ll];
	    file_out[n].dt[ll] = file_in[0].dt[ll];
	  }
	}
	mpp_end_def(file_out[n].fid);
	for(i=0; i<ndim; i++) {
	  if(file_out[n].axis[i].cart == 'T') continue;
	  mpp_put_var_value(file_out[n].fid, file_out[n].axis[i].vid, file_out[n].axis[i].data);
	  if( file_out[n].axis[i].bndtype > 0 )
	    mpp_put_var_value(file_out[n].fid, file_out[n].axis[i].bndid, file_out[n].axis[i].bnddata);
	}
      }
    }
  }
  
}; /* set_output_metadata */

/*******************************************************************************
   void get_field_attribute( )
   *******************************************************************************/
void get_field_attribute( int ntiles, Field_config *field)
{
  int n, l, nfield;
  char str[128];
  
  nfield = field->nvar;

  for(l=0; l<nfield; l++) {
    for(n=0; n<ntiles; n++) {
      field[n].var[l].missing = 0;
      field[n].var[l].scale   = 0;
      field[n].var[l].offset  = 0;
      field[n].var[l].vid = mpp_get_varid(*(field[n].fid), field[n].var[l].name);
      if( field[n].var[l].has_missing = mpp_var_att_exist(*(field[n].fid), field[n].var[l].vid, "missing_value") ) {
	mpp_get_var_att_double(*(field[n].fid), field[n].var[l].vid, "missing_value", &(field[n].var[l].missing));
      }
      /* also check for _FillValue */
      if( !field[n].var[l].has_missing ) {
         if( field[n].var[l].has_missing = mpp_var_att_exist(*(field[n].fid), field[n].var[l].vid, "_FillValue") ) {
           mpp_get_var_att_double(*(field[n].fid), field[n].var[l].vid, "_FillValue", &(field[n].var[l].missing));
        }
      }
      if(mpp_var_att_exist(*(field[n].fid), field[n].var[l].vid, "scale_factor")) {
	mpp_get_var_att_double(*(field[n].fid), field[n].var[l].vid, "scale_factor", &(field[n].var[l].scale) );
      }
      if(mpp_var_att_exist(*(field[n].fid), field[n].var[l].vid, "add_offset")) {
	mpp_get_var_att_double(*(field[n].fid), field[n].var[l].vid, "add_offset", &(field[n].var[l].offset) );
      }      
    }
  }

}; /* get_field_attribute */

/*******************************************************************************
void copy_field_attribute()
*******************************************************************************/
void copy_field_attribute( int ntiles_out, Field_config *field_in, Field_config *field_out)
{

  int n, l;
  
  for(n=0; n<ntiles_out; n++) 
    for(l=0; l<field_out->nvar; l++) {
      field_out[n].var[l].missing = field_in->var[l].missing;
      field_out[n].var[l].scale   = field_in->var[l].scale;
      field_out[n].var[l].offset  = field_in->var[l].offset;
    }
  

}
  

/*******************************************************************************
void set_remap_file( )
*******************************************************************************/
void set_remap_file( int ntiles, const char *mosaic_file, const char *remap_file, Interp_config *interp,
		     unsigned int *opcode, int save_weight_only)
{
  int    i, len, m, fid, vid;
  size_t start[4], nread[4];
  char str1[STRING], tilename[STRING];
  int file_exist;
  
  if(!remap_file) return;
  
  for(i=0; i<4; i++) {
    start[i] = 0; nread[i] = 1;
  }
  nread[1] = STRING;
  
  len = strlen(remap_file);
  if(len >= STRING) mpp_error("setoutput_remap_file(fregrid_util): length of remap_file should be less than STRING");  
  if( strcmp(remap_file+len-3, ".nc")==0 ) {
    strncpy(str1, remap_file, len-3);
    str1[len-3] = 0;
  }
  else
    strcpy(str1, remap_file);

  (*opcode) |= WRITE;

  if(ntiles>1) {
     fid = mpp_open(mosaic_file, MPP_READ);
     vid   = mpp_get_varid(fid, "gridtiles");
  }
  
  for(m=0; m<ntiles; m++) {
    interp[m].file_exist = 0;
    if(ntiles > 1) {
      start[0] = m;
      mpp_get_var_value_block(fid, vid, start, nread, tilename);
      if(strlen(str1) + strlen(tilename) > STRING -5) mpp_error("set_output_remap_file(fregrid_util): length of str1 + "
								"length of tilename should be no greater than STRING-5");
      sprintf(interp[m].remap_file, "%s.%s.nc", str1, tilename);
    }
    else
      sprintf(interp[m].remap_file, "%s.nc", str1);
    /* check xgrid file to be read (=1) or write ( = 2) */
    if(!save_weight_only && mpp_file_exist(interp[m].remap_file)) {
      (*opcode) |= READ;
      interp[m].file_exist = 1;
    }
      
  }

  if(ntiles>1) mpp_close(fid);
  
};/* set_remap_file */


/*----------------------------------------------------------------------
  void write_output_axis_data( )
  write out time axis data of the output data file
  --------------------------------------------------------------------*/
void write_output_time(int ntiles, File_config *file, int level)
{
  int         i, n;
  size_t      start[4], nwrite[4];

  for(i=0; i<4; i++) {
    start[i] = 0; nwrite[i] = 1;
  }
  start[0] = level;   
  if( mpp_pe() == mpp_root_pe()) {
    for(n=0; n<ntiles; n++) {
      for(i=0; i<file[n].ndim; i++) {
	if(file[n].axis[i].cart == 'T') {
	  nwrite[1] = 1;
	  mpp_put_var_value_block(file[n].fid, file[n].axis[i].vid, start,
				 nwrite, &(file[n].axis[i].data[level]));
	  if(strcmp(file[n].axis[i].bndname, "none") ) {
	    nwrite[1] = 2;
	    mpp_put_var_value_block(file[n].fid, file[n].axis[i].bndid, start,
				   nwrite, &(file[n].axis[i].bnddata[level*2]));
	  }
	}
      }
      /* write out time_avg_info if exist */
      if(file[n].has_tavg_info) {
	nwrite[1] = 1;
	mpp_put_var_value_block(file[n].fid, file[n].id_t1, start, nwrite, &(file[n].t1[level]) );
	mpp_put_var_value_block(file[n].fid, file[n].id_t2, start, nwrite, &(file[n].t2[level]) );
	mpp_put_var_value_block(file[n].fid, file[n].id_dt, start, nwrite, &(file[n].dt[level]) );
      }
    }
  }
}; /* write_output_time */

/*---------------------------------------------------------------------------
  void get_input_data(Mosaic_config *input, int l)
  get the input data for the number l variable.
  -------------------------------------------------------------------------*/
void get_input_data(int ntiles, Field_config *field, Grid_config *grid, Bound_config *bound,
		    int varid, int level_z, int level_n, int level_t, int extrapolate, double stop_crit)
{
  int         halo, i, j, k, i1, i2, n, p;
  int         memsize, nx, ny, ndim, nbound, l, pos, nz;
  size_t      *start, *nread;
  double      *data;
  short       *data_i2;
  int         *data_i4;
  Data_holder *dHold;
  int         interp_method;
  double      missing_value;
  size_t      start2[4], nread2[4];
  
  missing_value = field->var[varid].missing;
  interp_method = field->var[varid].interp_method;
  if(interp_method == CONSERVE_ORDER1)
    halo = 0;
  else
    halo = 1;

  nz = 1;
  if( level_z < 0 ) nz = field->var[varid].nz;
  ndim = field->var[varid].ndim;
  if(ndim < 2) mpp_error("fregrid_util(get_input_data): ndim must be no less than 2");
  nread  = (size_t *)malloc(ndim*sizeof(size_t));
  start  = (size_t *)malloc(ndim*sizeof(size_t));
  for(i=0; i<ndim; i++) {
    start[i] = 0; nread[i] = 1;
  }
  pos = 0;
  if(field->var[varid].has_taxis) start[pos++] = level_t;
  if(field->var[varid].has_naxis) start[pos++] = level_n;
  if(field->var[varid].has_zaxis) {
    if( level_z <0 ) {
      nread[pos]   = field->var[varid].nz;
      start[pos++] = field->var[varid].kstart;
    }
    else
      start[pos++] = level_z; 
  }
  if(ndim != pos + 2) mpp_error("fregrid_util(get_input_data): mimstch between ndim and has_taxis/has_zaxis/has_naxis");
  
  /* first read input data for each tile */
  for(n=0; n<ntiles; n++) {
    nx = grid[n].nx;
    ny = grid[n].ny;
    memsize = 1;
    nread[pos]   = ny;
    nread[pos+1] = nx;
    memsize = (nx+2*halo)*(ny+2*halo)*nz;
    
    field[n].data = (double *)malloc(memsize*sizeof(double));
    
    if(halo ==0 )
      data = field[n].data;
    else {
      data = (double *)malloc(nx*ny*nz*sizeof(double));
      init_halo(field[n].data, nx, ny, nz, 1);
    }
    
    switch(field[n].var[varid].type) {
    case NC_SHORT:
      data_i2 = (short *)malloc(memsize*sizeof(short));
      mpp_get_var_value_block(*(field[n].fid), field[n].var[varid].vid, start, nread, data_i2);
      for(i=0; i<nx*ny*nz; i++) data[i] = data_i2[i];
      free(data_i2);
      break;
    case NC_INT:
      data_i4 = (int *)malloc(memsize*sizeof(int));
      mpp_get_var_value_block(*(field[n].fid), field[n].var[varid].vid, start, nread, data_i4);
      for(i=0; i<nx*ny*nz; i++) data[i] = data_i4[i];
      free(data_i4);
      break;      
    case NC_DOUBLE:case NC_FLOAT:
      mpp_get_var_value_block(*(field[n].fid), field[n].var[varid].vid, start, nread, data);
      break;
    default:
      mpp_error("fregrid_util(get_input_data): field type should be NC_INT, NC_SHORT, NC_FLOAT or NC_DOUBLE");
    }

    if(field[n].var[varid].scale != 0) {
      for(i=0; i<nx*ny*nz; i++)
	if(data[i] != missing_value) data[i] *= field[n].var[varid].scale;
    }
    if(field[n].var[varid].offset != 0) {
      for(i=0; i<nx*ny*nz; i++)
	if(data[i] != missing_value) data[i] += field[n].var[varid].offset;
    }      

    /* extrapolate data if needed */
    if(extrapolate) {
      double *tmp = NULL;

      tmp = (double *)malloc(nx*ny*nz*sizeof(double));
      for(i=0; i<nx*ny*nz; i++) tmp[i] = data[i];
      do_extrapolate(nx, ny, nz, grid[n].lont1D, grid[n].latt1D, tmp, data, grid[n].is_cyclic,
		     field[n].var[varid].missing, stop_crit );
      field[n].var[varid].has_missing = 0;
      free(tmp);
    }      
    if(halo != 0) {
      /* copy the data onto compute domain */
      for(k=0; k<nz; k++) for(j=0; j<ny; j++) for(i=0; i<nx; i++) {
	i1 = k*(nx+2*halo)*(ny+2*halo)+(j+halo)*(nx+2*halo)+i+halo;
	i2 = k*nx*ny + j*nx+i;
	field[n].data[i1] = data[i2];
     }
      free(data);
    }

    /* currently we are assuming area is in double or float type */
    if(field[n].var[varid].cell_measures) {
      int i;

      if( !field[n].area ) {
	field[n].area = (double *) malloc(nx*ny*sizeof(double));
      }
      for(i=0; i<4; i++) { start2[i] = 0; nread2[i] = 1; }
      i = 0;
      if(field[n].var[varid].area_has_taxis) start2[i++] = level_t;
      if(field[n].var[varid].area_has_naxis) start2[i++] = level_n;
      nread2[i] = ny;
      nread2[i+1] = nx;
      mpp_get_var_value_block(field[n].var[varid].area_fid, field[n].var[varid].area_vid, start2, nread2, field[n].area );

    }
  } 

  /* update halo when halo > 0 */
  if(halo > 0) {
    for(n=0; n<ntiles; n++) {
      nbound = bound[n].nbound;
      if(nbound > 0) {
	dHold = (Data_holder *)malloc(nbound*sizeof(Data_holder));
	for(l=0; l<nbound; l++) {
	  dHold[l].data = field[bound[n].tile2[l]].data;
	  dHold[l].nx = grid[bound[n].tile2[l]].nx+2;
	  dHold[l].ny = grid[bound[n].tile2[l]].ny+2;
	}
	update_halo(grid[n].nx+2, grid[n].ny+2, nz, field[n].data, &(bound[n]), dHold );
	for(l=0; l<nbound; l++) dHold[l].data = NULL;
	free(dHold);
      }
    }
    /* second order conservative interpolation, gradient need to be calculated */
    if(interp_method == CONSERVE_ORDER2) {
      for(n=0; n<ntiles; n++) {
	int is_true = 1;
	nx = grid[n].nx;
	ny = grid[n].ny;
	field[n].grad_x = (double *)malloc(nx*ny*nz*sizeof(double));
	field[n].grad_y = (double *)malloc(nx*ny*nz*sizeof(double));
	field[n].grad_mask = (int *)malloc(nx*ny*nz*sizeof(int));
	for(k=0; k<nz; k++) for(j=0; j<ny; j++) for(i=0; i<nx; i++) {
	  field[n].grad_mask[k*nx*ny+j*nx+i] = 0;
	}
	for(k=0; k<nz; k++) {
	  p = k*(nx+2)*(ny+2);
	  grad_c2l(&(grid[n].nx), &(grid[n].ny), field[n].data+p, grid[n].dx, grid[n].dy, grid[n].area,
		   grid[n].edge_w, grid[n].edge_e, grid[n].edge_s, grid[n].edge_n,
		   grid[n].en_n, grid[n].en_e, grid[n].vlon_t, grid[n].vlat_t, 
		   field[n].grad_x, field[n].grad_y, &is_true, &is_true, &is_true, &is_true);
	}
	/* where there is missing and using second order conservative interpolation, need to calculate mask for gradient */
	if( field[n].var[varid].has_missing ) {
	  int ip1, im1, jp1, jm1,kk,ii,jj;
	  double missing;
	  missing = field[n].var[varid].missing;
	  for(k=0; k<nz; k++) for(j=0; j<ny; j++) for(i=0; i<nx; i++) {
	    ii=i+1; ip1 = ii+1; im1 = ii-1; jj = j+1; jp1 = jj+1; jm1 = jj-1;
	    l = k*(nx+2)*(ny+2);
	    if(field[n].data[l+jm1*(nx+2)+im1] == missing || field[n].data[l+jm1*(nx+2)+ii ] == missing ||
	       field[n].data[l+jm1*(nx+2)+ip1] == missing || field[n].data[l+jj *(nx+2)+im1] == missing ||
	       field[n].data[l+jj *(nx+2)+ip1] == missing || field[n].data[l+jp1*(nx+2)+im1] == missing ||
	       field[n].data[l+jp1*(nx+2)+ii ] == missing || field[n].data[l+jp1*(nx+2)+ip1] == missing  )
	      field[n].grad_mask[k*nx*ny+j*nx+i] = 1;
	  }
	}
      }
    }
  }


  
}; /* get_input_data */

/*---------------------------------------------------------------------------
  void get_input_data(Mosaic_config *input, int l)
  get the input data for the number l variable.
  -------------------------------------------------------------------------*/
void get_test_input_data(char *test_case, double test_param, int ntiles, Field_config *field,
			 Grid_config *grid, Bound_config *bound, unsigned int opcode)
{
  int         halo, i, j, k, ii, n, nx, ny, l, nbound;
  double      *data;
  Data_holder *dHold;
  char input_file[128];
  int  fid, vid, dim[2]; 
  
  if(opcode & CONSERVE_ORDER1)
    halo = 0;
  else
    halo = 1;
  
  for(n=0; n<ntiles; n++) {
    nx = grid[n].nx;
    ny = grid[n].ny;
    field[n].data = (double *)malloc((nx+2*halo)*(ny+2*halo)*sizeof(double));
    data          = (double *)malloc(nx*ny*sizeof(double));
    if(!strcmp(test_case,"tanh_cosphi_costheta") ) {
      for(j=0; j<ny; j++) for(i=0; i<nx; i++) {
	data[j*ny+i] = tanh(test_param*cos(grid[n].lont[(j+1)*(nx+2)+i+1])*cos(grid[n].latt[(j+1)*(nx+2)+i+1]));
      }
    }
    if(!strcmp(test_case,"tanh_sinphi_sintheta") ) {
      for(j=0; j<ny; j++) for(i=0; i<nx; i++) {
	data[j*ny+i] = tanh(test_param*sin(grid[n].lont[(j+1)*(nx+2)+i+1])*sin(grid[n].latt[(j+1)*(nx+2)+i+1]));
      }
    }    
    else if(!strcmp(test_case,"cosphi_costheta") ) {
      for(j=0; j<ny; j++) for(i=0; i<nx; i++) {
	data[j*ny+i] = cos(grid[n].lont[(j+1)*(nx+2)+i+1])*cos(grid[n].latt[(j+1)*(nx+2)+i+1]);
      }
    }
    else if(!strcmp(test_case,"sinphi_costheta") ) {
      for(j=0; j<ny; j++) for(i=0; i<nx; i++) {
	data[j*ny+i] = sin(grid[n].lont[(j+1)*(nx+2)+i+1])*cos(grid[n].latt[(j+1)*(nx+2)+i+1]);
      }
    }
    else
      mpp_error("fregrid_util: invalid choice of test_case");

    /* write out input data */
    sprintf(input_file, "%s.input.tile%d.nc", test_case, n+1);
    fid = mpp_open(input_file, MPP_WRITE);
    dim[0] = mpp_def_dim(fid, "grid_y", ny);
    dim[1] = mpp_def_dim(fid, "grid_x", nx);
    vid    = mpp_def_var(fid, "data", NC_DOUBLE, 2, dim, 0);
    mpp_end_def(fid);
    mpp_put_var_value(fid, vid, data);
    mpp_close(fid);
    
    for(j=0; j<ny; j++) for(i=0; i<nx; i++) {
      ii = (j+halo)*(nx+2*halo)+i+halo;
      field[n].data[ii] = data[j*nx+i];
    }
    free(data);
  }
    
  /* update halo when halo > 0 */
  if(halo > 0) {
    for(n=0; n<ntiles; n++) {
      nbound = bound[n].nbound;
      if(nbound > 0) {
	dHold = (Data_holder *)malloc(nbound*sizeof(Data_holder));
	for(l=0; l<nbound; l++) {
	  dHold[l].data = field[bound[n].tile2[l]].data;
	  dHold[l].nx = grid[bound[n].tile2[l]].nx+2;
	  dHold[l].ny = grid[bound[n].tile2[l]].ny+2;
	}
	update_halo(grid[n].nx+2, grid[n].ny+2, 1, field[n].data, &(bound[n]), dHold );
	for(l=0; l<nbound; l++) dHold[l].data = NULL;
	free(dHold);
      }
    }
    /* second order conservative interpolation, gradient need to be calculated */
    if(opcode & CONSERVE_ORDER2) {
      for(n=0; n<ntiles; n++) {
	int is_true = 1;
	field[n].grad_x = (double *)malloc((grid[n].nx+2)*(grid[n].ny+2)*sizeof(double));
	field[n].grad_y = (double *)malloc((grid[n].nx+2)*(grid[n].ny+2)*sizeof(double));
	grad_c2l(&(grid[n].nx), &(grid[n].ny), field[n].data, grid[n].dx, grid[n].dy, grid[n].area,
		 grid[n].edge_w, grid[n].edge_e, grid[n].edge_s, grid[n].edge_n,
		 grid[n].en_n, grid[n].en_e, grid[n].vlon_t, grid[n].vlat_t,
		 field[n].grad_x, field[n].grad_y, &is_true, &is_true, &is_true, &is_true);
      }
    }
  }
  
}; /* get_test_input_data */


void allocate_field_data(int ntiles, Field_config *field, Grid_config *grid, int nz)
{
  int n, i;
  size_t memsize;
  
  for(n=0; n<ntiles; n++) {
    memsize = grid[n].nx*grid[n].ny*nz;
    field[n].data = (double *)malloc(memsize*sizeof(double));
  }

  
}; /* allocate_field_data */


/*-------------------------------------------------------------------------
  write_field_data(Mosaic_config *output)
  write data to output file
  -----------------------------------------------------------------------*/
void write_field_data(int ntiles, Field_config *field, Grid_config *grid, int varid, int level_z, int level_n, int level_t)
{
  double *gdata;
  double missing_value;
  short  *data_i2;
  int    *data_i4;
  int    nx, ny, nz, n, ndim, i, j, data_size, pos;
  size_t *nwrite, *start;

  ndim = field->var[varid].ndim;
  if(ndim < 2) mpp_error("fregrid_util(write_field_data): ndim must be no less than 2");

  nwrite = (size_t *)malloc(ndim*sizeof(size_t));
  start  = (size_t *)malloc(ndim*sizeof(size_t));
  nz = 1;
  if(level_z<0) nz = field->var[varid].nz;
  for(i=0; i<ndim; i++) {
    start[i] = 0; nwrite[i] = 1;
  }

  missing_value = field->var[varid].missing;
  pos = 0;
  if(field->var[varid].has_taxis) start[pos++] = level_t;
  if(field->var[varid].has_naxis) start[pos++] = level_n;
  if(field->var[varid].has_zaxis) {
    if(level_z < 0)
      nwrite[pos++] = nz;
    else
      start[pos++] = level_z;
  }
  if(ndim != pos + 2) mpp_error("fregrid_util(write_field_data): mimstch between ndim and has_taxis/has_zaxis/has_naxis");
  
  for(n=0; n<ntiles; n++) {
    /* global data onto root pe */
    nx = grid[n].nx;
    ny = grid[n].ny;
    nwrite[pos]   = ny;
    nwrite[pos+1] = nx;    
    data_size = nx*ny*nz;

    if(field[n].var[varid].offset != 0) {
      for(i=0; i<nx*ny*nz; i++)
	if(field[n].data[i] != missing_value) field[n].data[i] -= field[n].var[varid].offset;
    }
    if(field[n].var[varid].scale != 0) {
      for(i=0; i<nx*ny*nz; i++)
	if(field[n].data[i] != missing_value) field[n].data[i] /= field[n].var[varid].scale;
    }

    if(mpp_npes() == 1) {
      gdata = field[n].data;
    }
    else {
      gdata = (double *)malloc(nx*ny*nz*sizeof(double));
      mpp_global_field_double_3D(grid[n].domain, grid[n].nxc, grid[n].nyc, nz,
			      field[n].data, gdata);
    }
    
    switch(field[n].var[varid].type) {
    case NC_SHORT:
      data_i2 = (short *)malloc(data_size*sizeof(short));
      for(i=0; i<data_size; i++) data_i2[i] = (short)gdata[i];
      if(mpp_pe() == mpp_root_pe()) mpp_put_var_value_block(*(field[n].fid), field[n].var[varid].vid, start, nwrite, data_i2);
      free(data_i2);
      break;
    case NC_INT:
      data_i4 = (int *)malloc(data_size*sizeof(int));
      for(i=0; i<data_size; i++) data_i4[i] = (int)gdata[i];
      if(mpp_pe() == mpp_root_pe()) mpp_put_var_value_block(*(field[n].fid), field[n].var[varid].vid, start, nwrite, data_i4);
      free(data_i4);
      break;
    case NC_DOUBLE:case NC_FLOAT:
      if(mpp_pe() == mpp_root_pe())mpp_put_var_value_block(*(field[n].fid), field[n].var[varid].vid, start, nwrite, gdata);
      break;
    default:
      mpp_error("fregrid_util(write_field_data): field type should be NC_SHORT, NC_FLOAT or NC_DOUBLE");
    }
    if(mpp_npes() != 1) free(gdata);
  }

};/* write_output_data */

void get_contact_direction(int ncontact, const int *tile, const int *istart, const int *iend,
			   const int *jstart, const int *jend, int *dir)
{
  int n;

  for(n=0; n<ncontact; n++) {
    if(istart[n] == iend[n] && jstart[n] == jend[n])
      mpp_error("fregrid_util: istart = iend and jstart = jend can not be both true for one contact");
    if(istart[n] != iend[n] && jstart[n] != jend[n])
      mpp_error("fregrid_util: either istart = iend or jstart = jend need to be true");
    if(istart[n] == iend[n]) {
      if(istart[n] == 0)
	dir[n] = WEST;
      else
	dir[n] = EAST;
    }
    else {
      if(jstart[n] == 0)
	dir[n] = SOUTH;
      else
	dir[n] = NORTH;
    }
  }

}

void setup_boundary(const char *mosaic_file, int ntiles, Grid_config *grid, Bound_config *bound, int halo, int position)
{
  int ncontacts, shift, n, nbound, l, l2, nb, nx, ny;
  int *tile, *dir;
  int *istart, *iend, *jstart, *jend;

  ncontacts = read_mosaic_ncontacts(mosaic_file);
  if(ncontacts == 0) {
    for(n=0; n<ntiles; n++) bound[n].nbound = 0;
    return;
  }
  if(ntiles ==1) grid[0].is_cyclic = 1;
  
  if(position == CENTER)
    shift = 0;
  else if(position == CORNER)
    shift = 1;
  else
    mpp_error("fregrid_util: position should be CENTER or CORNER");
  
  tile   = (int *)malloc(2*ncontacts*sizeof(int));
  istart = (int *)malloc(2*ncontacts*sizeof(int));
  iend   = (int *)malloc(2*ncontacts*sizeof(int));
  jstart = (int *)malloc(2*ncontacts*sizeof(int));
  jend   = (int *)malloc(2*ncontacts*sizeof(int));
  dir    = (int *)malloc(2*ncontacts*sizeof(int));
  read_mosaic_contact(mosaic_file, tile, tile+ncontacts, istart, iend, jstart, jend,
		      istart+ncontacts, iend+ncontacts, jstart+ncontacts, jend+ncontacts);
  for(n=0; n<2*ncontacts; n++) --tile[n];
  get_contact_direction(2*ncontacts, tile, istart, iend, jstart, jend, dir);
  
  /* First find number of boundary for each tile */
  for(n=0; n<ntiles; n++) {
    nbound = 0;
    nx = grid[n].nx;
    ny = grid[n].ny;
    for(l=0; l<2*ncontacts; l++) {
      if(tile[l] == n) nbound++;
    }
    bound[n].nbound = nbound;
    if(nbound > 0) {
      bound[n].is1    = (int *)malloc(nbound*sizeof(int));
      bound[n].ie1    = (int *)malloc(nbound*sizeof(int));
      bound[n].js1    = (int *)malloc(nbound*sizeof(int));
      bound[n].je1    = (int *)malloc(nbound*sizeof(int));
      bound[n].is2    = (int *)malloc(nbound*sizeof(int));
      bound[n].ie2    = (int *)malloc(nbound*sizeof(int));
      bound[n].js2    = (int *)malloc(nbound*sizeof(int));
      bound[n].je2    = (int *)malloc(nbound*sizeof(int));      
      bound[n].rotate = (int *)malloc(nbound*sizeof(int));
      bound[n].tile2  = (int *)malloc(nbound*sizeof(int));
      nb = 0;
      for(l=0; l<2*ncontacts; l++) {
	if(tile[l] != n) continue;
	switch(dir[l]) {
	case WEST:
	  bound[n].is1[nb]  = 0;
	  bound[n].ie1[nb]  = halo-1;
	  bound[n].js1[nb]  = min(jstart[l],jend[l])+halo;
	  bound[n].je1[nb]  = max(jstart[l],jend[l])+halo+shift;
	  break;
	case EAST:
	  bound[n].is1[nb]  = nx+shift+halo;
	  bound[n].ie1[nb]  = nx+shift+halo+halo-1;
	  bound[n].js1[nb]  = min(jstart[l],jend[l])+halo;
	  bound[n].je1[nb]  = max(jstart[l],jend[l])+halo+shift;
	  break;
	case SOUTH:
	  bound[n].is1[nb]  = min(istart[l],iend[l])+halo;
	  bound[n].ie1[nb]  = max(istart[l],iend[l])+halo+shift;
	  bound[n].js1[nb]  = 0;
	  bound[n].je1[nb]  = halo-1;
	  break;
	case NORTH:
	  bound[n].is1[nb]  = min(istart[l],iend[l])+halo;
	  bound[n].ie1[nb]  = max(istart[l],iend[l])+halo+shift;
	  bound[n].js1[nb]  = ny+shift+halo;
	  bound[n].je1[nb]  = ny+shift+halo+halo-1;
	  break;	    
	}
	l2 = (l+ncontacts)%(2*ncontacts);
	bound[n].tile2[nb] = tile[l2];
	switch(dir[l2]) {
	case WEST:
	  bound[n].is2[nb]  = halo+shift;
	  bound[n].ie2[nb]  = halo+shift+halo-1;
	  bound[n].js2[nb]  = min(jstart[l2],jend[l2])+halo;
	  bound[n].je2[nb]  = max(jstart[l2],jend[l2])+halo+shift;
	  break;
	case EAST:
	  bound[n].is2[nb]  = nx-halo+1;
	  bound[n].ie2[nb]  = nx;
	  bound[n].js2[nb]  = min(jstart[l2],jend[l2])+halo;
	  bound[n].je2[nb]  = max(jstart[l2],jend[l2])+halo+shift;
	  break;
	case SOUTH:
	  bound[n].is2[nb]  = min(istart[l2],iend[l2])+halo;
	  bound[n].ie2[nb]  = max(istart[l2],iend[l2])+halo+shift;
	  bound[n].js2[nb]  = halo+shift;
	  bound[n].je2[nb]  = halo+shift+halo-1;
	  break;
	case NORTH:
	  bound[n].is2[nb]  = min(istart[l2],iend[l2])+halo;
	  bound[n].ie2[nb]  = max(istart[l2],iend[l2])+halo+shift;
	  bound[n].js2[nb]  = ny-halo+1;
	  bound[n].je2[nb]  = ny;
	  break;	    
	}
	bound[n].rotate[nb] = ZERO;
	if(dir[l] == WEST && dir[l2] == NORTH) bound[n].rotate[nb] = NINETY;
	if(dir[l] == EAST && dir[l2] == SOUTH) bound[n].rotate[nb]= NINETY;
	if(dir[l] == SOUTH && dir[l2] == EAST) bound[n].rotate[nb] = MINUS_NINETY;
	if(dir[l] == NORTH && dir[l2] == WEST) bound[n].rotate[nb] = MINUS_NINETY;      
	if(dir[l] == NORTH && dir[l2] == NORTH) bound[n].rotate[nb] = ONE_HUNDRED_EIGHTY;
	/* make sure the size match at the boundary */
	if( (bound[n].ie2[nb]-bound[n].is2[nb]+1)*(bound[n].je2[nb]-bound[n].js2[nb]+1) !=
	    (bound[n].ie1[nb]-bound[n].is1[nb]+1)*(bound[n].je1[nb]-bound[n].js1[nb]+1) )
	  mpp_error("fregrid_util: size mismatch between the boundary");
	nb++;
      }      
    }
  }
}; /* setup_boundary */

void delete_bound_memory(int ntiles, Bound_config *bound)
{
  int n;
  
  for(n=0; n<ntiles; n++) {
    if(bound[n].nbound > 0) {
      free(bound[n].is1);
      free(bound[n].ie1);
      free(bound[n].js1);
      free(bound[n].je1);
      free(bound[n].is2);
      free(bound[n].ie2);
      free(bound[n].js2);
      free(bound[n].je2);
      free(bound[n].tile2);
      free(bound[n].rotate);
    }
  }
}
  

/*-----------------------------------------------------------------------------
  void init_halo(double *var, int nx, int ny, int nz, int halo)
  initialze the halo data to be zero.
  ---------------------------------------------------------------------------*/
void init_halo(double *var, int nx, int ny, int nz, int halo)
{
  int i, j, k;
  int nxd, nyd, nall;

  nxd = nx+2*halo;
  nyd = ny+2*halo;
  nall = nxd*nyd;
  
  for(k=0; k<nz; k++) {
    for(j=0; j<nyd; j++) for(i=0; i<halo; i++) var[k*nall+j*nxd+i] = 0; /* west halo */
    for(j=0; j<nyd; j++) for(i=nx+halo; i<nxd; i++) var[k*nall+j*nxd+i] = 0; /* east halo */
    for(j=0; j<halo; j++) for(i=0; i<nxd; i++) var[k*nall+j*nxd+i] = 0; /* south halo */
    for(j=ny+halo; j<nyd; j++) for(i=0; i<nxd; i++) var[k*nall+j*nxd+i] = 0; /* north halo */
  }

};/* init_halo */

 
void update_halo(int nx, int ny, int nz, double *data, Bound_config *bound, Data_holder *dHold)
{
  int nbound, n, i, j, k, l, size1, size2, nx2, ny2;
  int is1, ie1, js1, je1, is2, ie2, js2, je2, bufsize;
  double *buffer;
  
  nbound = bound->nbound;
  size1  = nx*ny;

  for(n=0; n<nbound; n++) {
    is1 = bound->is1[n];
    ie1 = bound->ie1[n];
    js1 = bound->js1[n];
    je1 = bound->je1[n];
    is2 = bound->is2[n];
    ie2 = bound->ie2[n];
    js2 = bound->js2[n];
    je2 = bound->je2[n];
    nx2 = dHold[n].nx;
    ny2 = dHold[n].ny;
    size2 = nx2*ny2;
    bufsize = nz*(ie2-is2+1)*(je2-js2+1);
    buffer = (double *)malloc(bufsize*sizeof(double));
    /* fill the buffer */
    l = 0;
    switch(bound->rotate[n]) {
    case ZERO:
      for(k=0; k<nz; k++) for(j=js2; j<=je2; j++) for(i=is2; i<=ie2; i++) buffer[l++] = dHold[n].data[k*size2+j*nx2+i];
      break;
    case NINETY:
      for(k=0; k<nz; k++) for(i=ie2; i>=is2; i--) for(j=js2; j<=je2; j++) buffer[l++] = dHold[n].data[k*size2+j*nx2+i];
      break;
    case MINUS_NINETY:
      for(k=0; k<nz; k++) for(i=is2; i<=ie2; i++) for(j=je2; j>=js2; j--) buffer[l++] = dHold[n].data[k*size2+j*nx2+i];
      break;
    case ONE_HUNDRED_EIGHTY:
      for(k=0; k<nz; k++) for(j=je2; j>=js2; j--) for(i=ie2; i>=is2; i--) buffer[l++] = dHold[n].data[k*size2+j*nx2+i];
      break;
    }
    l = 0;
    for(k=0; k<nz; k++) for(j=js1; j<=je1; j++) for(i=is1; i<=ie1; i++) data[k*size1+j*nx+i] = buffer[l++];
    free(buffer);
  }  

}

/* do_extrapolate assume the input data is on a lat-lon grid */

void do_extrapolate (int ni, int nj, int nk, const double *lon, const double *lat, const double *data_in,
		     double *data_out, int is_cyclic, double missing_value, double stop_crit)
{
  int i, j, k, n, n1, n2, n3, n4, n5, n6;
  double initial_guess = 0.0;
  double latp, latm, cfc;
  double cstr, csm, csj;
  double resmax;
  double *dyu=NULL, *dyt=NULL;
  double *dxu=NULL, *dxt=NULL;
  double *cfn=NULL, *cfe=NULL;
  double *cfs=NULL, *cfw=NULL;
  double *tmp=NULL, *sor=NULL, *res=NULL;
  
  /* construct grid factors for a sphere */
  dxu = (double *)malloc(ni*sizeof(double));
  dxt = (double *)malloc(ni*sizeof(double));
  dyu = (double *)malloc(nj*sizeof(double));
  dyt = (double *)malloc(nj*sizeof(double));
  
  for(j=0; j<nj-1; j++) dyu[j] = lat[j+1] - lat[j];
  dyu[nj-1] = dyu[nj-2];
  for(j=1; j<nj; j++) dyt[j] = 0.5*(dyu[j] + dyu[j-1]);
  dyt[0] = dyt[1];
  for(i=0; i<ni-1; i++) dxu[i] = lon[i+1] - lon[i];
  dxu[ni-1] = dxu[ni-2];
  for(i=1; i<ni; i++) dxt[i] = 0.5*(dxu[i] + dxu[i-1]);
  dxt[0] = dxt[1];  

  cfn = (double *)malloc(ni*nj*sizeof(double));
  cfe = (double *)malloc(ni*nj*sizeof(double));
  cfs = (double *)malloc(ni*nj*sizeof(double));
  cfw = (double *)malloc(ni*nj*sizeof(double));
  for(j=0; j<nj; j++) {
    if (j == nj-1)
      latp = lat[j] + 0.5*(lat[j] - lat[j-1]);
    else
      latp = 0.5*(lat[j] + lat[j+1]);

    if (j == 0) 
      latm = lat[j] - 0.5*(lat[j+1] - lat[j]);
    else
      latm = 0.5*(lat[j] + lat[j-1]);
    csj  = cos(latp);
    csm  = cos(latm);
    cstr = 1.0/cos(lat[j]);
    for(i=0; i<ni; i++) {
      n = j*ni+i;
      cfn[n] = csj*cstr/(dyt[j]*dyu[j]);
      cfs[n] = csm*cstr/(dyt[j]*dyu[max(j-1,0)]);
      cfe[n] = cstr*cstr/(dxu[i]*dxt[i]);
      cfw[n] = cstr*cstr/(dxu[max(i-1,0)]*dxt[i]);
      cfc    = 1.0/(cfn[n]+cfs[n]+cfe[n]+cfw[n]);
      cfn[n] = cfn[n]*cfc;
      cfs[n] = cfs[n]*cfc;
      cfe[n] = cfe[n]*cfc;
      cfw[n] = cfw[n]*cfc;        
    }
  }

  tmp = (double *)malloc((ni+2)*(nj+2)*sizeof(double));
  sor = (double *)malloc( ni   *nj    *sizeof(double));
  res = (double *)malloc( ni   *nj    *sizeof(double));
  
  for(n=0; n<(ni+2)*(nj+2); n++) tmp[n] = 0.0;

  for(j=0; j<nj; j++) for(i=0; i<ni; i++) {
    n = j*ni+i;
    n1 = (j+1)*(ni+2)+i+1;
    if( fabs(data_in[n] - missing_value) <= EPSLN10 ) tmp[n1] = initial_guess;
  }

  for(k=0; k<nk; k++) {
    for(j=0; j<nj; j++) for(i=0; i<ni; i++) {
      n = k*ni*nj + j*ni + i;
      n1 = j*ni + i;
      n2 = (j+1)*(ni+2) + i+1;
      if(fabs(data_in[n] - missing_value) <= EPSLN10 )
	sor[n1] = REL_COEF;
      else {
	tmp[n2] = data_in[n];
	sor[n1] = 0.0;
      }
    }

    fill_boundaries(ni, nj, tmp, is_cyclic);

    /* iterate */
    for(n=0; n<MAX_ITER; n++) {
      resmax=0.0;
      for(j=0; j<nj; j++) for(i=0; i<ni; i++) {
	n1 = j*ni + i;
	n2 = (j+1)*(ni+2) + i+1;
	n3 = (j+1)*(ni+2) + i  ;
	n4 = (j+1)*(ni+2) + i+2;
	n5 = (j  )*(ni+2) + i+1;
	n6 = (j+2)*(ni+2) + i+1;
	res[n1] = cfw[n1]*tmp[n3] + cfe[n1]*tmp[n4] + cfs[n1]*tmp[n5] + cfn[n1]*tmp[n6] - tmp[n2];
      }
      for(j=0; j<nj; j++) for(i=0; i<ni; i++) {
	n1 = j*ni + i;
	n2 = (j+1)*(ni+2) + i+1;
	res[n1] *= sor[n1];
	tmp[n2] += res[n1];
	resmax = max(fabs(res[n1]), resmax);
      }

      if(resmax <= stop_crit || n == MAX_ITER-1) {
	for(j=0; j<nj; j++) for(i=0; i<ni; i++) {
	  n1 = k*ni*nj + j*ni + i;
	  n2 = (j+1)*(ni+2) + i+1;
	  data_out[n1] = tmp[n2];
	}
	break;
      }

      /* update boundaries */

      fill_boundaries(ni, nj, tmp, is_cyclic);
    }

    if(mpp_pe() == mpp_root_pe() ) printf("Stopped after %d iterations, maxres = %g\n", n, resmax);
  }

  /* release memory */
  free(dxt);
  free(dxu);
  free(dyt);
  free(dyu);
  free(cfn);
  free(cfe);
  free(cfs);
  free(cfw);
  free(tmp);
  free(sor);
  free(res);
  
} /* do_extrapolate */


void fill_boundaries(int ni, int nj, double *data, int is_cyclic)
{
  int i,j;

  if(is_cyclic) {
    for(j=1; j<=nj; j++) {
      data[j*(ni+2)]      = data[j*(ni+2)+ni];
      data[j*(ni+2)+ni+1] = data[j*(ni+2)+1];
    }
  }
}
 
