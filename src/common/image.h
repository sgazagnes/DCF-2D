#ifndef IMAGE_H
#define IMAGE_H

#include <FreeImage.h>
#include "cmdline.h"
#include "types.h"

/* Read */
value *ReadBasic(const char *fname, size_t *dims, int *bitpix);
value *ReadHDF5(const char *fname, const char *dataset, int *n_dims, size_t *dims, size_t *dims_T,
		size_t *grid);
value *ReadFITS(char *fname,  int *n_dims, size_t *dims, size_t *dims_t,
		int *bitpix, size_t* grid);
value *ReadInput( struct gengetopt_args_info args, int *n_dims, size_t *dims,
		  size_t *dims_t, int *bitpix);
/* Write */
void WriteOutput(struct gengetopt_args_info args, value *img, int bitpix, int n_dims,
		 size_t *dims, size_t *dims_T);
void WriteDAP( struct gengetopt_args_info args, value *out, value *outDH, value *outScale,
	       int bitpix, int n_dims, size_t *dims_T, size_t *dims);
void WriteBasic(const char *fname_out, value *im, int bitpix, size_t *dims);
void WriteHDF5(const char* fname_in, const char* fname_out,  const char *dataset_out, value *out,
	       int n_dims, size_t *dims_T, size_t *grid);
void WriteFITS(char *fnamein, char *fnameout,  value *out,  int bitpix, int n_dims,
	       size_t *dims_T, size_t *grid);
/* Annexes */
FIBITMAP* GenericLoader(const char* lpszPathName, int flag);
value *UpsideDown(value *img, size_t width, size_t height);
#endif
