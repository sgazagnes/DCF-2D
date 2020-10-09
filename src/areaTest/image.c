#define _GNU_SOURCE
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <FreeImage.h>
#include <string.h>

#include "image.h"
#include "checks.h"
#include "logc.h"
#include "types.h"
#include "mpihelper.h"

/** Generic image loader
@param lpszPathName Pointer to the full file name
@param flag Optional load flag constant
@return Returns the loaded dib if successful, returns NULL otherwise
*/

/*value FindMax(value *im,size_t* width, size_t *height){
  value max = im[0];
  ulong i,j, size=width ;
  for (i = 1; i<=size;i++){
    max = (im->Pixmap[i]>max) ? im->Pixmap[i] : max;
  }
  return max;
  }*/



value *ReadImage( struct gengetopt_args_info args, size_t* width, size_t *height, short *bitsperpixel)
{   
  value *img;
  char *imgfilename;
  if (asprintf(&imgfilename, "%s%d.%s", args.inprefix_arg, rank(), args.intype_arg) < 0) {
    fprintf(stderr, "%d: could not allocate memory for input filename\n", rank());
    MPI_Abort(MPI_COMM_WORLD, 132);
  }
  img = ReadBasic(imgfilename,width,height, bitsperpixel);

  return img;
} /* ReadImage */


value *ReadBasic(char *fnm, size_t *width, size_t *height, short *bitsperpixel){
   FIBITMAP *dib = GenericLoader(fnm,0);

   if (dib == NULL) return NULL;
     
   *bitsperpixel =  FreeImage_GetBPP(dib);
   *height = FreeImage_GetHeight(dib);
   *width = FreeImage_GetWidth(dib);
   
   size_t imsize = (*width)*(*height);
   value *im = calloc(imsize, sizeof(value));
   check_alloc(im,010);
   
   if(*bitsperpixel == 8){
     size_t i=0;
     for(size_t y = 0; y < *height; y++) {
       BYTE *bits = (BYTE *)FreeImage_GetScanLine(dib, *height - y -1);
       for(size_t x = 0; x < *width; x++,i++) {
	 im[i] = bits[x];
       }
     }
     FreeImage_Unload(dib);
     return im;
   }else if(*bitsperpixel == 16){
     size_t i=0;
     for(size_t y = 0; y < *height; y++) {
       unsigned short *bits = (unsigned short *)FreeImage_GetScanLine(dib, y);
       for(size_t x = 0; x < *width; x++,i++) {
	 im[i] = bits[x];
       }
     }
     FreeImage_Unload(dib);
     return im;  
   }else{ 
     FreeImage_Unload(dib);
     return NULL;
   }
}

void WriteImages( struct gengetopt_args_info args, value *out, value *outDH, value *outScale, short bitsperpixel, size_t width, size_t height) {
  char *outfilename;
  if(!strcmp(args.intype_arg, "pgm") || !strcmp(args.intype_arg, "tif") || !strcmp(args.intype_arg, "png") || !strcmp(args.intype_arg, "bmp") || !strcmp(args.intype_arg, "jpg")){
   
    if (asprintf(&outfilename, "%sL%d.%s", args.outprefix_arg, rank(), args.outtype_arg) < 0) {
      error("could not allocate memory for output filename\n");
      MPI_Abort(MPI_COMM_WORLD, 131);
    }

    info("Writing image to %s, %zux%zu pixels", outfilename, width, height);
    WriteBasic(outfilename, out, bitsperpixel, width, height);
    info("Wrote image");   
    free(outfilename);
    
    if (asprintf(&outfilename, "%sC%d.%s", args.outprefix_arg, rank(), args.outtype_arg) < 0) {
      error("could not allocate memory for output filename\n");
      MPI_Abort(MPI_COMM_WORLD, 131);
    }

    info("Writing image to %s, %zux%zu pixels", outfilename, width, height);
    WriteBasic(outfilename, outDH, bitsperpixel,  width, height);
    info("Wrote image");
    free(outfilename);
    
    if (asprintf(&outfilename, "%sS%d.%s", args.outprefix_arg, rank(), args.outtype_arg) < 0) {
      error("could not allocate memory for output filename\n");
      MPI_Abort(MPI_COMM_WORLD, 131);
    }

    info("Writing image to %s, %zux%zu pixels", outfilename, width, height);
    WriteBasic(outfilename, outScale, bitsperpixel,  width, height);
    info("Wrote image");
  }else{
    error("Code does not handle other files than png, bmp, jpg, pgm of tif");
    MPI_Abort(MPI_COMM_WORLD, 132);
  }
    free(outfilename);
} /* writeImage */


void WriteImage(struct gengetopt_args_info args, value *img, short bitsperpixel, size_t width, size_t height)
{
  char *fname;
  if (asprintf(&fname, "%s%d.%s", args.outprefix_arg, rank(), args.outtype_arg) < 0) {
    fprintf(stderr, "%d: could not allocate memory for input filename\n", rank());
    MPI_Abort(MPI_COMM_WORLD, 132);
  }
   WriteBasic(fname, img, bitsperpixel, width, height);
} /* WritePGM */




 
void WriteBasic( char *outfilename, value *im, short bitspp, size_t width, size_t height){
  FIBITMAP *outmap;

  int j,y;
  FREE_IMAGE_FORMAT fif = FreeImage_GetFIFFromFilename(outfilename);
  if(fif == FIF_UNKNOWN) {
    // no signature ?
    // try to guess the file format from the file extension
    fif = FreeImage_GetFIFFromFilename(outfilename);
  }
  // check that the plugin has reading capabilities ...
  if((fif != FIF_UNKNOWN) && FreeImage_FIFSupportsWriting(fif)) {

     
    if (bitspp == 8){
      ubyte *imagebuf = malloc(width *height* sizeof(ubyte));;
    
      outmap = FreeImage_AllocateT(FIT_BITMAP,width,height,bitspp,0xFF,0xFF,0xFF);
      for ( j=height-1,  y=0; j>=0; j--, y++){
	imagebuf = FreeImage_GetScanLine(outmap,j);      
	for (size_t x=0;x<width;x++)
	  imagebuf[x]=im[width*y + x];
      
      }	 
    } else if  (bitspp == 16){
      unsigned short *imagebuf;
      outmap = FreeImage_AllocateT(FIT_UINT16,width,height,16,0xFFFF,0xFFFF,0xFFFF);
      for (j=height-1, y=0; j>=0; j--, y++){      
	imagebuf = (unsigned short *)FreeImage_GetScanLine(outmap,j);
	for (size_t x=0; x<width ;x++)
	  imagebuf[x]=im[(width)*j + x];
	
      }
    } else{
      error("not handling more than 16 bits");
      MPI_Abort(MPI_COMM_WORLD, 132);
    }
       
    FreeImage_Save(fif,outmap,outfilename,0); 
    FreeImage_Unload(outmap);
  }
  else{
    error("FreeImage couldn't write the output \n");
    MPI_Abort(MPI_COMM_WORLD, 132);
  }
 

}

FIBITMAP* GenericLoader(const char* lpszPathName, int flag) {
  FREE_IMAGE_FORMAT fif = FIF_UNKNOWN;
//  check the file signature and deduce its format
// (the second argument is currently not used by FreeImage)
  fif = FreeImage_GetFileType(lpszPathName, 0);
  if(fif == FIF_UNKNOWN) {
    // no signature ?
    // try to guess the file format from the file extension
    fif = FreeImage_GetFIFFromFilename(lpszPathName);
  }
  // check that the plugin has reading capabilities ...
  if((fif != FIF_UNKNOWN) && FreeImage_FIFSupportsReading(fif)) {
    // ok, let's load the file
    FIBITMAP *dib = FreeImage_Load(fif, lpszPathName, flag);
    // unless a bad file format, we are done !
    return dib;
  }
  else{
    error("FreeImage couldn't read the input \n");
    MPI_Abort(MPI_COMM_WORLD, 132);
  }
  return NULL;
}


