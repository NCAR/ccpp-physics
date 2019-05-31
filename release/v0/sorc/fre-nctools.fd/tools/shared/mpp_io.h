/****************************************************************
                    mpp_io.h
   This headers defines interface to read and write netcdf file. All the data
will be written out from root pe. 

   contact: Zhi.Liang@noaa.gov

****************************************************************/
#ifndef MPP_IO_H_
#define MPP_IO_H_
#include <netcdf.h>

#define MPP_WRITE 100
#define MPP_READ  200
#define MPP_APPEND  300
#define MPP_INT NC_INT
#define MPP_DOUBLE NC_DOUBLE
#define MPP_CHAR NC_CHAR
#define HEADER_BUFFER_VALUE (16384)
#define MPP_FILL_INT NC_FILL_INT
#define MPP_FILL_DOUBLE NC_FILL_DOUBLE

int mpp_open(const char *file, int action);
void mpp_close(int ncid);
int mpp_get_nvars(int fid);
void mpp_get_varname(int fid, int varid, char *name);
int mpp_get_varid(int fid, const char *varname);
int mpp_get_dimid(int fid, const char *dimname);
int mpp_get_dimlen(int fid, const char *name);
int mpp_get_var_natts(int fid, int vid);
void mpp_get_var_attname(int fid, int vid, int i, char *name);
void mpp_copy_att_by_name(int fid_in, int vid_in, int fid_out, int vid_out, const char *name);
void mpp_get_var_value(int fid, int vid, void *data);
void mpp_get_var_value_block(int fid, int vid, const size_t *start, const size_t *nread, void *data);
void mpp_get_var_att(int fid, int vid, const char *name, void *val);
void mpp_get_var_att_double(int fid, int vid, const char *name, double *val);
void mpp_get_global_att(int fid, const char *name, void *val);
int mpp_get_var_ndim(int fid, int vid);
nc_type mpp_get_var_type(int fid, int vid);
char mpp_get_var_cart(int fid, int vid);
void mpp_get_var_dimname(int fid, int vid, int ind, char *name);
char mpp_get_dim_cart(int fid, const char *name);
void mpp_get_var_bndname(int fid, int vid, char *bndname);
int mpp_var_att_exist(int fid, int vid, const char *att);
int mpp_global_att_exist(int fid, const char *att);
int mpp_def_dim(int fid, const char* name, int size);
int mpp_def_var(int fid, const char* name, nc_type type, int ndim, const int *dims, int natts, ...);
void mpp_def_global_att(int fid, const char *name, const char *val);
void mpp_def_global_att_double(int fid, const char *name, size_t len, const double *val);
void mpp_def_var_att(int fid, int vid, const char *attname, const char *attval);
void mpp_def_var_att_double(int fid, int vid, const char *attname, double attval);
void mpp_copy_var_att(int fid_in, int vid_in, int fid_out, int vid_out);
void mpp_copy_global_att(int fid_in, int fid_out);
void mpp_put_var_value(int fid, int vid, const void* data);
void mpp_put_var_value_block(int fid, int vid, const size_t *start, const size_t *nread, const void *data);
void mpp_end_def(int fid);
void mpp_redef(int fid);
int mpp_file_exist(const char *file);
int mpp_field_exist(const char *file, const char *field);
int mpp_var_exist(int fid, const char *field);
int mpp_dim_exist(int fid, const char *dimname);
int get_great_circle_algorithm(int fid);
#endif
