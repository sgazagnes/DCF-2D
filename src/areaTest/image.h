#ifndef IMAGE_H
#define IMAGE_H

#include <FreeImage.h>
#include "cmdline.h"
#include "types.h"

FIBITMAP* GenericLoader(const char* lpszPathName, int flag);
value *ReadImage( struct gengetopt_args_info args, size_t* width, size_t *height, short *bitsperpixel);
value *ReadBasic(char *fnm, size_t *width, size_t *height, short *bitsperpixel);
void WriteImage(struct gengetopt_args_info args, value *img, short bitsperpixel, size_t width, size_t height);
void WriteImages( struct gengetopt_args_info args, value *out, value *outDH, value *outScale, short bitsperpixel, size_t width, size_t height);
void WriteBasic( char *outfilename, value *im, short bitspp, size_t width, size_t height);
#endif
