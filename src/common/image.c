#define _GNU_SOURCE
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <FreeImage.h>
#include <string.h>
#include <hdf5.h>
#include <fitsio.h>

#include "image.h"
#include "checks.h"
#include "logc.h"
#include "types.h"
#include "mpihelper.h"


/*+++++++++++++++++++++++++++++++++++++++*/
/*     					 */
/*	     Functions Read     	 */
/*                                       */
/*+++++++++++++++++++++++++++++++++++++++*/


value *ReadInput( struct gengetopt_args_info args, int *n_dims, size_t *dims,
		   size_t *dims_T, int *bitpix){
  /* Main Read Function */
  value 	*img;
  char 		*fname;
  int		myrank = rank();
  size_t 	grid[2] = {args.grid_arg[0], args.grid_arg[1]};

  if (!strcmp(args.intype_arg, "h5")){
    /* HDF5 files */	
    info("Reading HDF5 file");
    if (asprintf(&fname, "%s.%s", args.inprefix_arg, args.intype_arg) < 0){
      error("Could not allocate memory for input filename");
      MPI_Abort(MPI_COMM_WORLD, 001);
    }
    img = ReadHDF5(fname, args.dataset_arg, n_dims, dims, dims_T, grid);
    if (args.bitsperpixel_orig == NULL)
      warn("Bits per pixel not specified, 8 bits is assumed for the HDF5 dataset.");    
    *bitpix = args.bitsperpixel_arg;
    
  } else if (!strcmp(args.intype_arg, "fits") || !strcmp(args.intype_arg, "fit")){
     /* FITS files */
    info("Reading FITS Image (not adapted for FITS Tables !!!!)");
    if (asprintf(&fname, "%s.%s", args.inprefix_arg, args.intype_arg) < 0) {
      error("Could not allocate memory for input filename");
      MPI_Abort(MPI_COMM_WORLD, 001);
    }
    img = ReadFITS(fname, n_dims, dims, dims_T, bitpix, grid);
    if (args.bitsperpixel_orig != NULL){
      /* If the user specificed a bit depth different from the FITS BITPIX keyword */
      if (*bitpix != args.bitsperpixel_arg)
	warn("The specified bit depth of the image (%d) does not match the image header (%d)",
	     args.bitsperpixel_arg, *bitpix);
      *bitpix = args.bitsperpixel_arg;
    }
    
  } else{
    /* Assuming pgm, png, bmp, jpeg, tif ... */
    if (asprintf(&fname, "%s-%d.%s", args.inprefix_arg, myrank, args.intype_arg) < 0) {
      error("Could not allocate memory for input filename");
      MPI_Abort(MPI_COMM_WORLD, 001);
    }
    info("Reading image %s with FreeImage", fname);
    img = ReadBasic(fname,  dims, bitpix);
    if (args.bitsperpixel_orig != NULL){
      /* If the user specificed a bit depth different from the FITS BITPIX keyword */
      if (*bitpix != args.bitsperpixel_arg)
	warn("The specified bit depth of the image (%d) does not match the image header (%d)",
	     args.bitsperpixel_arg, *bitpix);
      *bitpix = args.bitsperpixel_arg;
    }
    *n_dims = 2;
    dims_T[0] = grid[0] * dims[0];
    dims_T[1] = grid[1] * dims[1];
  }
  free(fname);
  return img;
} /* ReadImage */



value *ReadFITS(char *fname,  int *n_dims, size_t *dims, size_t *dims_T,
		 int *bitpix, size_t* grid){
  /* Read Fits images */
  fitsfile 	*fptr;   		/* FITS file pointer */
  value 	*img;			/* Pointer to data */
  int 		status = 0;   		/* CFITSIO status value MUST be initialized to zero! */
  int 		myrank = rank();	/* Rank of process */
  long 		inc[2] = {1,1};		/* Increment for fits read functiion */
  long 		naxes[2] = {1,1};	/* Dimensions (2D max) */
  long 		counts[2] = {1,1};	/* Pixels to read in fits function */
  long 		offsets[2] = {1,1};	/* Pixel index to start in fist function */
  
  if (!fits_open_file(&fptr, fname, READONLY, &status)) {	/* Open file */
    if (!fits_get_img_param(fptr, 2,  bitpix,  n_dims,
			    naxes, &status)){			/* Get image parameter */
      if (*n_dims > 2 || *n_dims == 0){	   			/* 2D Only */
	error("only 1D or 2D images are supported");
	MPI_Abort(MPI_COMM_WORLD, 002);
      } else if (*bitpix < 1 || *bitpix > 16){ 			/* 16 bits max */
	error("only <= 16 bits  supported");
	MPI_Abort(MPI_COMM_WORLD, 003);
      } else {
	counts[0]  = naxes[0] / grid[0]; 			/* Number of pixels to read (width) */
	counts[1]  = naxes[1] / grid[1]; 			/* Number of pixels to read (width) */
	offsets[0] = (myrank % grid[0]) * counts[0] + 1;   	/* First pixel to read (width) */   
	offsets[1] = (myrank / grid[0]) * counts[1] + 1;	/* First pixel to read (height) */

	if((myrank%grid[0]) < (naxes[0] % grid[0])){
	  counts[0]++;
	  offsets[0] += (myrank%grid[0]);
	}else
	  offsets[0] += (naxes[0] % grid[0]);
    
	if((myrank/grid[0]) < (naxes[1] % grid[1])){
	  counts[1]++;
	  offsets[1] += (myrank/grid[0]);
	}else
	  offsets[1] += (naxes[1] % grid[1]);
		
	dims_T[0]  = naxes[0];					/* Dimensions of the image in dims_T */
	dims_T[1]  = naxes[1];
	dims[0] = counts[0];			       		/* Tile dimensions stored in dims */
	dims[1] = counts[1];
	img     =  malloc(counts[1]*counts[0] * sizeof(value));	/* Allocating for the tile size */
	check_alloc(img, 001);
	counts[0] += offsets[0]-1;				/* Last pixel to read (width) */
	counts[1] += offsets[1]-1;				/* Last pixel to read (height) */
	    
	info("Tile dim: %ld by %ld , offsets %ld %ld, counts %ld and %ld ",
	     dims[0],dims[1],offsets[0],offsets[1],counts[0],counts[1]);
	fits_read_subset(fptr, TUSHORT, offsets, counts, inc, NULL, img,
			 NULL,&status); 			/* Reading the image subset */
      }
    }
    fits_close_file(fptr, &status);
  } else {
    error("Cannot open the file");
    MPI_Abort(MPI_COMM_WORLD, 004);
  }    
  info("Image of %ld x %ld  pixels. Bits per pixel = %d.",
       dims_T[0], dims_T[1], *bitpix);
  return img;
}

value *ReadHDF5(const char *fname, const char *dataset, int *n_dims, size_t *dims, size_t *dims_T,
		 size_t *grid){
  /* Read HDF5 files */
  hid_t       	file_id, dataset_id, dataspace, type;  		/* identifiers */
  int 		myrank = rank();
  hsize_t 	*hdims, hn_dims;
  hsize_t 	counts[2] = {1,1};
  hsize_t 	offsets[2] = {1,1};
  
  file_id = H5Fopen(fname, H5F_ACC_RDONLY, H5P_DEFAULT);	/* Open HDF5 file */
  if (file_id < 0){
    error("Could not open file ");
    MPI_Abort(MPI_COMM_WORLD, 004);
  }
  
  dataset_id = H5Dopen2(file_id,dataset, H5P_DEFAULT);		/* Open dataset */
  if (dataset_id < 0) {
    error("Could not open dataset ");
    MPI_Abort(MPI_COMM_WORLD, 005);
  }
  
  dataspace = H5Dget_space(dataset_id);				/* Copy dataset */
  hn_dims = H5Sget_simple_extent_ndims(dataspace);
 /* Get number of dimensions */
  if(hn_dims > 2){
    error("Only handle 2D data");
    MPI_Abort(MPI_COMM_WORLD, 002);
  }
  
  hdims = calloc(hn_dims,sizeof(hsize_t));			/* Allocate dimension array */
  H5Sget_simple_extent_dims(dataspace, hdims, NULL);		/* Get dimensions */
  info("Image dimensions: %ld by %ld \n", hdims[0], hdims[1]);	/*WARN:  0 is height,  1 is width */
  dims_T[1] = hdims[0];
  dims_T[0] = hdims[1];
  free(hdims);
  
  counts[0]  = dims_T[1] / grid[1];				/* Number of pixels to read (height) */
  offsets[0] = (myrank / grid[0]) * counts[0];			/* First pixel index (width) */
  counts[1]  = dims_T[0] / grid[0];				/* Number of pixels to read (width) */
  offsets[1] = (myrank % grid[0]) * counts[1];			/* First pixel index (height) */

  if((myrank/grid[0]) < (dims_T[1] % grid[1])){
    counts[0]++;
    offsets[0] += (myrank/grid[0]);
  }else
    offsets[0] += (dims_T[1] % grid[1]);

  if((myrank%grid[0]) < (dims_T[0] % grid[0])){
    counts[1]++;
    offsets[1] += (myrank%grid[0]);
  }else
    offsets[1] += (dims_T[0] % grid[0]);
    
  info("Tile dimension: height %ld, width %ld.\n First pixel offset: height %ld, width %ld",
       counts[0], counts[1], offsets[0], offsets[1]);
  
  value *img = malloc(counts[0]*counts[1] * sizeof(value));
  check_alloc(img,001);
  hid_t memory_window = H5Screate_simple(hn_dims, counts, NULL); /* Allocate memory to read the data */
  H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, offsets, NULL,
		      counts, NULL);				/* Select the data subset */
  if(sizeof(value) == 1)					/* Define type of data to read */
    type =  H5T_STD_U8LE;
  else if(sizeof(value) == 2)
    type =  H5T_STD_U16LE;
  else if(sizeof(value) == 4)
    type =  H5T_STD_U32LE;
  else{
    error("Issue with format");
    MPI_Abort(MPI_COMM_WORLD, 006);
  }

  herr_t err = H5Dread(dataset_id, type, memory_window,
		       dataspace, H5P_DEFAULT, img); 		/* Read the data */
  if (err < 0){ 
    error("Could not read data from the dataset");
    MPI_Abort(MPI_COMM_WORLD, 005);
  }
  
  dims[0] = counts[1];
  dims[1] = counts[0];
  *n_dims =  hn_dims;
  
  H5Dclose(dataset_id);
  H5Sclose(dataspace);
  H5Fclose(file_id);
  return img;
}

value *ReadBasic(const char *fname, size_t *dims, int *bitpix){
  /* Reading classic types of images */
  FIBITMAP 	*dib = GenericLoader(fname,0);
  size_t 	i = 0;/* Generic Loader from FreeImage */
  if (dib == NULL){
    error("An issue happened while reading the image");
    MPI_Abort(MPI_COMM_WORLD, 004);
  }
  
  //*bitpix = FreeImage_GetBPP(dib);				/* Get dynamic range */
  dims[0] = FreeImage_GetWidth(dib);				/* Get width */
  dims[1] = FreeImage_GetHeight(dib);				/* Get height */

  FREE_IMAGE_COLOR_TYPE nc = FreeImage_GetColorType(dib);

  /*  if(nc != 1){
      FIBITMAP* hImage= dib;
      dib = FreeImage_ConvertTo32Bits( hImage );
      hImage= dib;
      dib = FreeImage_ConvertTo8Bits( hImage );
      FreeImage_Unload( hImage );
      }*/
  if(nc != 1)
    dib =  FreeImage_ConvertToGreyscale(dib);
  *bitpix = FreeImage_GetBPP(dib);		
  value *img = malloc(dims[0]*dims[1] * sizeof(value));		
  check_alloc(img,001);
   
  if (*bitpix <= 8){
    for(size_t y = 0; y < dims[1]; y++) {
      BYTE *bits = (BYTE *)FreeImage_GetScanLine(dib, dims[1] - y -1);
      for(size_t x = 0; x < dims[0]; x++,i++)
	img[i] = bits[x];
    }
    FreeImage_Unload(dib);
    return img;
    
  } else if (*bitpix <= 16){
    for(size_t y = 0; y < dims[1]; y++) {
      unsigned short *bits = (unsigned short *)FreeImage_GetScanLine(dib, y);
      for(size_t x = 0; x < dims[0]; x++,i++) 
	img[i] = bits[x];
    }
    FreeImage_Unload(dib);
    return img;  
  }else{
    error("> 16 bits are not supported ");
    MPI_Abort(MPI_COMM_WORLD, 003);
    FreeImage_Unload(dib);
    return NULL;
  }
}

/*+++++++++++++++++++++++++++++++++++++++*/
/*     					 */
/*	     Functions Write      	 */
/*                                       */
/*+++++++++++++++++++++++++++++++++++++++*/

void WriteOutput(struct gengetopt_args_info args, value *img, int bitpix, int n_dims,
		 size_t *dims, size_t *dims_T){
  /*Main Write Function */
  char 		*fname_out, *fname_in, *dataset_out;
  int 		myrank = rank();
  size_t 	grid[2] = {args.grid_arg[0], args.grid_arg[1]};
  
  if(args.outtype_arg == NULL)
     args.outtype_arg = args.intype_arg;
  
  if(!strcmp(args.outtype_arg, "h5")){
    info("Writing HDF5 file");
    if (asprintf(&fname_in, "%s.%s", args.inprefix_arg, args.intype_arg) < 0) {
      error("Could not allocate memory for input filename");
      MPI_Abort(MPI_COMM_WORLD, 006);
    }
    if (asprintf(&fname_out, "%s.%s", args.outprefix_arg, args.outtype_arg) < 0) {
      error("Could not allocate memory for output filename");
      MPI_Abort(MPI_COMM_WORLD, 006);
    }
    asprintf(&dataset_out, "%s%ld", args.filter_arg,args.lambda_arg);
    WriteHDF5(fname_in, fname_out,dataset_out, img, n_dims, dims_T, grid);
  }else if(!strcmp(args.outtype_arg, "fits") || !strcmp(args.outtype_arg, "fit")){
    info("Writing FITS file");
    if(!strcmp(args.intype_arg, "fits")){
      if (asprintf(&fname_in, "%s.%s", args.inprefix_arg, args.intype_arg) < 0) {
	error("Could not allocate memory for input filename");
	MPI_Abort(MPI_COMM_WORLD, 006);
      }
    }
    if (asprintf(&fname_out, "%s.%s", args.outprefix_arg, args.outtype_arg) < 0) {
      error("Could not allocate memory for output filename");
      MPI_Abort(MPI_COMM_WORLD, 006);
    }
    WriteFITS(fname_in,fname_out, img, bitpix, n_dims, dims_T, grid);
  }else{
    info("Writing using FreeImage");
    if(!strcmp(args.intype_arg, "fits")){
      img = UpsideDown(img, dims[0], dims[1]);
      myrank = (grid[1]-1-myrank/grid[0])*grid[0] + myrank%grid[0];
    }
    if (asprintf(&fname_out, "%s-%d.%s", args.outprefix_arg, myrank, args.outtype_arg) < 0) {
      error("Could not allocate memory for output filename");
      MPI_Abort(MPI_COMM_WORLD, 006);
    }
    WriteBasic(fname_out, img, bitpix, dims);
  }

} /* WritePGM */


void WriteDAP( struct gengetopt_args_info args, value *out, value *outDH, value *outScale,
	       int bitpix, int n_dims, size_t *dims_T, size_t *dims) {
  /* Specific DAP function */
  char		*fname_in = NULL, *fname_out;
  int		myrank = rank();
  size_t 	grid[2] = {args.grid_arg[0], args.grid_arg[1]};

  if(args.outtype_arg == NULL)
    args.outtype_arg = args.intype_arg;
   if(!strcmp(args.outtype_arg, "h5")){
    info("Writing HDF5 file");
    if (asprintf(&fname_in, "%s.%s", args.inprefix_arg, args.intype_arg) < 0) {
      error("Could not allocate memory for input filename");
      MPI_Abort(MPI_COMM_WORLD, 006);
    }
    if (asprintf(&fname_out, "%s.%s", args.outprefix_arg, args.outtype_arg) < 0) {
      error("Could not allocate memory for output filename");
      MPI_Abort(MPI_COMM_WORLD, 006);
    }
    info("Writing file %s", fname_out);
    WriteHDF5(fname_in,fname_out,"Contrast", out, n_dims, dims_T, grid);
    MPI_Barrier(MPI_COMM_WORLD);
    WriteHDF5(fname_out,fname_out,"Luminance", outDH, n_dims, dims_T, grid);
    MPI_Barrier(MPI_COMM_WORLD);
    WriteHDF5(fname_out,fname_out,"Scale", outScale, n_dims, dims_T, grid);

   }else if(!strcmp(args.outtype_arg, "fits") || !strcmp(args.outtype_arg, "fit")){
     info("Writing FITS file");
     if(!strcmp(args.intype_arg, "fits")){
       if (asprintf(&fname_in, "%s.%s", args.inprefix_arg, args.intype_arg) < 0) {
	 error("Could not allocate memory for input filename");
	 MPI_Abort(MPI_COMM_WORLD, 006);
       }
     }
     if (asprintf(&fname_out, "%sLuminance.%s", args.outprefix_arg, args.outtype_arg) < 0) {
       error("Could not allocate memory for output filename");
       MPI_Abort(MPI_COMM_WORLD, 006);
     }
     WriteFITS(fname_in,fname_out, out, bitpix, n_dims, dims_T, grid);
     free(fname_out);
     MPI_Barrier(MPI_COMM_WORLD);

     if (asprintf(&fname_out, "%sContrast.%s", args.outprefix_arg, args.outtype_arg) < 0) {
       error("Could not allocate memory for output filename");
       MPI_Abort(MPI_COMM_WORLD, 006);
     }
     WriteFITS(fname_in,fname_out, outDH, bitpix, n_dims, dims_T, grid);
     free(fname_out);
     MPI_Barrier(MPI_COMM_WORLD);

     if (asprintf(&fname_out, "%sScale.%s", args.outprefix_arg, args.outtype_arg) < 0) {
       error("Could not allocate memory for output filename");
       MPI_Abort(MPI_COMM_WORLD, 006);
     }
     WriteFITS(fname_in,fname_out, outScale, bitpix, n_dims, dims_T, grid);
     free(fname_out);
   } else{
      if(!strcmp(args.intype_arg, "fits")){
       out = UpsideDown(out, dims[0], dims[1]);
       outDH = UpsideDown(outDH, dims[0], dims[1]);
       outScale = UpsideDown(outScale, dims[0], dims[1]);
       myrank = (grid[1]-1-myrank/grid[0])*grid[0] + myrank%grid[0];
      }
     FREE_IMAGE_FORMAT fif = FIF_UNKNOWN;
     if (asprintf(&fname_out, "%sL-%d.%s", args.outprefix_arg, myrank, args.outtype_arg) < 0) {
       error("could not allocate memory for output filename");
       MPI_Abort(MPI_COMM_WORLD, 006);
     }
  
     fif = FreeImage_GetFileType(fname_out, 0);
     if(fif == FIF_UNKNOWN) 
       fif = FreeImage_GetFIFFromFilename(fname_out);
     if((fif != FIF_UNKNOWN) && FreeImage_FIFSupportsReading(fif)) {
       WriteBasic(fname_out, out, bitpix, dims);
       info("Wrote image Luminance");   
       free(fname_out);
    
       if (asprintf(&fname_out, "%sC-%d.%s", args.outprefix_arg, myrank, args.outtype_arg) < 0) {
	 error("could not allocate memory for output filename");
	 MPI_Abort(MPI_COMM_WORLD, 006);
       }
       WriteBasic(fname_out, outDH, bitpix, dims);
       info("Wrote image Contrast");
       free(fname_out);

       if (asprintf(&fname_out, "%sS-%d.%s", args.outprefix_arg, myrank, args.outtype_arg) < 0) {
	 error("could not allocate memory for output filename");
	 MPI_Abort(MPI_COMM_WORLD, 006);
       }
       WriteBasic(fname_out, outScale, bitpix, dims);
       info("Wrote image Scale");
       free(fname_out);
     }
     else{
       error("FreeImage couldn't create the new images");
       MPI_Abort(MPI_COMM_WORLD, 006);
     }
   }
} /* writeImage */

 
void WriteBasic(const char *fname_out, value *im, int bitpix, size_t *dims){
  /*Write png, pgm, bmp, tif ... */
  FIBITMAP 	*outmap;
  int 		j,y;
  size_t 	width = dims[0];
  size_t 	height = dims[1];
  
  FREE_IMAGE_FORMAT fif = FreeImage_GetFIFFromFilename(fname_out);
  if (fif == FIF_UNKNOWN) 
    fif = FreeImage_GetFIFFromFilename(fname_out);
  if ((fif != FIF_UNKNOWN) && FreeImage_FIFSupportsWriting(fif)) {
    if (bitpix <= 8){
      ubyte *imagebuf = malloc(width *height* sizeof(ubyte));;
    
      outmap = FreeImage_AllocateT(FIT_BITMAP,width,height,8,0xFF,0xFF,0xFF);
      for ( j=height-1,  y=0; j>=0; j--, y++){
	imagebuf = FreeImage_GetScanLine(outmap,j);      
	for (size_t x=0;x<width;x++)
	  imagebuf[x]=im[width*y + x];
      }	 
    } else if (bitpix <= 16){
      unsigned short *imagebuf;
      outmap = FreeImage_AllocateT(FIT_UINT16,width,height,16,0xFFFF,0xFFFF,0xFFFF);
      for (j=height-1, y=0; j>=0; j--, y++){      
	imagebuf = (unsigned short *)FreeImage_GetScanLine(outmap,j);
	for (size_t x=0; x<width ;x++)
	  imagebuf[x]=im[(width)*j + x];	
      }
    } else{
      error("not handling more than 16 bits");
      MPI_Abort(MPI_COMM_WORLD, 003);
    }
       
    FreeImage_Save(fif,outmap,fname_out,0); 
    FreeImage_Unload(outmap);
  }
  else{
    error("FreeImage couldn't write the output \n");
    MPI_Abort(MPI_COMM_WORLD, 007);
  }
}


void WriteHDF5(const char* fname_in, const char* fname_out,  const char *dataset_out, value *out,
	       int n_dims, size_t *dims_T, size_t *grid){
  /* Writing HDF5 files */
  hid_t       	file_id, dataset_id, dataspace, type;  /* identifiers */

  hsize_t	hdims[2] = {dims_T[1], dims_T[0]};
  hsize_t 	offsets[2] = {1,1};
  hsize_t 	counts[2] = {hdims[0] / grid[1], hdims[1] / grid[0]};
  int         	message;
  int 		myrank = rank();

  if(sizeof(value) == 1)
    type =  H5T_STD_U8LE;
  else if(sizeof(value) == 2)
    type =  H5T_STD_U16LE;
  else if(sizeof(value) == 4)
    type =  H5T_STD_U32LE;
  else{
    error("Non supported value type");
    MPI_Abort(MPI_COMM_WORLD, 007);
  }
  
  /*Writing new HDF5 file */
  if(myrank ==0){
    if(!strcmp(fname_in, fname_out)){
      file_id = H5Fopen(fname_out, H5F_ACC_RDWR, H5P_DEFAULT);
      info("Writing in the same input file");
    }else{
      file_id = H5Fcreate(fname_out, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
      info("File created ");
    }
    hid_t global_memory_space = H5Screate_simple(n_dims, hdims, NULL);
    H5Eset_auto (H5E_DEFAULT, NULL, NULL);
    dataset_id = H5Dopen(file_id,dataset_out, H5P_DEFAULT);		/* Open dataset*/
    if (dataset_id < 0) {
      info("Creating dataset");
      dataset_id = H5Dcreate(file_id,dataset_out, type, global_memory_space,
			     H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);	/* Create dataset */
    }
    dataspace = H5Dget_space(dataset_id);			/* Copy dataset */
  }else{
    MPI_Recv(&message, 1, MPI_INT, myrank - 1, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    file_id = H5Fopen(fname_out, H5F_ACC_RDWR, H5P_DEFAULT);
    dataset_id = H5Dopen(file_id,dataset_out, H5P_DEFAULT);
    dataspace = H5Dget_space(dataset_id);
  }

  offsets[0] = (myrank / grid[0]) * counts[0];			/* First pixel index (height) */
  offsets[1] = (myrank % grid[0]) * counts[1];			/* First pixel index (width) */
  
  if((myrank/grid[0]) < (hdims[0] % grid[1])){
    counts[0]++;
    offsets[0] += (myrank/grid[0]);
  }else
    offsets[0] += (hdims[0] % grid[1]);

  if((myrank%grid[0]) < (hdims[1] % grid[0])){
    counts[1]++;
    offsets[1] += (myrank%grid[0]);
  }else
    offsets[1] += (hdims[1] % grid[0]);

  info("Tile dimension: height %ld, width %ld.\n First pixel offset: height %ld, width %ld",
       counts[0], counts[1], offsets[0], offsets[1]);
  hid_t memory_window = H5Screate_simple((hsize_t)n_dims, counts, NULL);
  H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, offsets, NULL, counts, NULL);
  H5Dwrite(dataset_id, type, memory_window, dataspace, H5P_DEFAULT, out);

  H5Dclose(dataset_id);
  H5Sclose(dataspace);
  H5Fclose(file_id);
  info("HDF5 closed");
  if(myrank < np()-1)
    MPI_Send(&myrank, 1, MPI_INT, myrank+1, 1, MPI_COMM_WORLD);
}

void WriteFITS(char *fnamein, char *fnameout,  value *out,  int bitpix, int n_dims,
	       size_t *dims_T, size_t *grid){
  /* Writing FITS images */
  fitsfile 	*outfptr, *infptr; 	/* FITS file pointers */
  int 		status = 0; 		/* CFITSIO status value MUST be initialized to zero! */
  int 		message;
  int  		myrank = rank(); 
  long 		naxes[2] = {dims_T[0],dims_T[1]};      
  long 		counts[2] = {naxes[0]/grid[0],naxes[1]/grid[1]};
  long 		offsets[2] = {1,1};
  char 		str1[100];
  strcpy(str1, "!");   		/* '!' symbol makes the output image, if existing, to be overwritten */

  offsets[0] = (myrank % grid[0]) * counts[0] + 1;   	/* First pixel to read (width) */   
  offsets[1] = (myrank / grid[0]) * counts[1] + 1;	/* First pixel to read (height) */

  if((myrank%grid[0]) < (naxes[0] % grid[0])){
    counts[0]++;
    offsets[0] += (myrank%grid[0]);
  }else
    offsets[0] += (naxes[0] % grid[0]);
    
  if((myrank/grid[0]) < (naxes[1] % grid[1])){
    counts[1]++;
    offsets[1] += (myrank/grid[0]);
  }else
    offsets[1] += (naxes[1] % grid[1]);

  counts[0] += offsets[0]-1;				/* Last pixel to read (width) */
  counts[1] += offsets[1]-1;
      
  if(myrank == 0){
    if (!fits_create_file(&outfptr, strcat(str1, fnameout), &status)){
      // copy all the header keywords from original image to new output file
      info("Image created");
      if(fnamein != NULL){
	fits_open_file(&infptr, fnamein, READONLY, &status); // open input images
	fits_copy_header(infptr, outfptr, &status);
	fits_get_img_param(infptr, 2, &bitpix, &n_dims, naxes, &status);
	fits_close_file(infptr, &status);

      }else{
	info("Creating");
	if(bitpix <=8)
	  fits_create_img(outfptr, BYTE_IMG, n_dims, naxes, &status);
	else if(bitpix <= 16)
	  fits_create_img(outfptr, USHORT_IMG, n_dims, naxes, &status);
	else
	  fits_create_img(outfptr, ULONG_IMG, n_dims, naxes, &status);
      }

      // info("Tile dimension: width %ld, height %ld.\n First pixel offset: width %ld, height %ld",
      //	   counts[0], counts[1], offsets[0], offsets[1]);

      fits_write_subset(outfptr, TUSHORT, offsets, counts, out, &status);
    }else {
      error("Cannot create the file");
      MPI_Abort(MPI_COMM_WORLD, 007);
    }
  }else{
    MPI_Recv(&message, 1, MPI_INT, myrank - 1, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    fits_open_file(&outfptr, fnameout, READWRITE, &status); // open input images	

    // info("Tile dimension: width %ld, height %ld.\n First pixel offset: width %ld, height %ld",
    //	counts[0], counts[1], offsets[0], offsets[1]);
    
    fits_write_subset(outfptr, TUSHORT, offsets, counts, out, &status);
  }
  fits_close_file(outfptr, &status);
  if(rank() < np()-1){
    MPI_Send(&status, 1, MPI_INT, myrank+1, 1, MPI_COMM_WORLD);
  }
} /* WriteFITS */


/*+++++++++++++++++++++++++++++++++++++++*/
/*     					 */
/*	     Functions annexes  	 */
/*                                       */
/*+++++++++++++++++++++++++++++++++++++++*/


value *UpsideDown(value *img, size_t width, size_t height){
  value *imgB = calloc(width*height,sizeof(value));
  for(size_t i = 0; i<height*width ;i++)
    imgB[i] = img[i];
  
  for(size_t i = 0; i <height; i++)
    for(size_t j = 0; j<width; j++)
      img[i*width +j] = imgB[(height-i-1)*width + j];
  free(imgB);
  return img;
} /* UpsideDown */


FIBITMAP* GenericLoader(const char* lpszPathName, int flag) {
  FREE_IMAGE_FORMAT fif = FIF_UNKNOWN;
  //  check the file signature and deduce its format
  fif = FreeImage_GetFileType(lpszPathName, 0);
  if(fif == FIF_UNKNOWN) {
    // no signature ? try to guess the file format from the file extension
    fif = FreeImage_GetFIFFromFilename(lpszPathName);
  }
  // check that the plugin has reading capabilities ...
  if((fif != FIF_UNKNOWN) && FreeImage_FIFSupportsReading(fif)) {
    FIBITMAP *dib = FreeImage_Load(fif, lpszPathName, flag);
    // unless a bad file format, we are done !
    return dib;
  }
  else{
    error("FreeImage couldn't read the input \n");
    MPI_Abort(MPI_COMM_WORLD, 132);
  }
  return NULL;
} /* GenericLoader */
