#ifndef CHECKST_H
#define CHECKST_H

#include "types.h"

void check_alloc(void *array, int code);
void check_not_null(void *ptr, int code);
void check_mpi_error(int errorval, int code);
void check_file_close(int errorval, const char *filename);
void check_area_size(size_t size, ulong area);
void check_boundary(Boundary *b);
void printfrank(const char *format, int my_rank, ...);
void PrintParents(MaxNode *maxtree, size_t width, size_t height, int myrank);
void PrintBound(Boundary *b, int myrank );

#endif
