#include <stdlib.h>
#include <stdio.h>
#include <mpi.h>
#include <stdarg.h>

#include "mpihelper.h"
#include "checks.h"
#include "logc.h"
#include "types.h"

#define MAX_MSG_SIZE 256

void check_alloc(void *array, int code) {
  if (NULL == array) {
    error("Allocation failed, out of memory? code: %d", code);
    MPI_Abort(MPI_COMM_WORLD, code);
  }
} /* check_alloc */

void check_not_null(void *ptr, int code) {
  if (NULL == ptr) {
    error("Pointer is null. code: %d", code);
    MPI_Abort(MPI_COMM_WORLD, code);
  }
} /* check_not_null */

void check_mpi_error(int errorval, int code) {
  if (MPI_SUCCESS != errorval) {
    error("MPI error %d! %d", errorval, code);
    MPI_Abort(MPI_COMM_WORLD, code);
  }
} /* check_mpi_error */

void check_file_close(int errorval, const char *filename) {
  if (0 != errorval) {
    error("Error closing file %s : errorval", filename);
    MPI_Abort(MPI_COMM_WORLD, errorval);
  }
}

void check_area_size(size_t size, ulong area) {
  if (area != size) {
    error("Area does not equal size after flooding, aborting.");
    MPI_Abort(MPI_COMM_WORLD, 133);
  }
}


void check_boundary(Boundary *b) {
  info("Boundary exists of array %p offset %zu %zu %zu %zu %zu, size %zu \n",
	 (void *)  b->array, b->offset[0], b->offset[1], b->offset[2], b->offset[3], b->offset[4], b->size);
  /* for (size_t a = 0; a < b->size; a++) {
    printf("%zu: index %ld parent %ld area %lu gval %hi filter %hi borderindex %ld borderparent %ld \n",
          a, (b->array[a]).index, (b->array[a]).parent, (b->array[a]).area, (b->array[a]).gval, (b->array[a]).filter,
          (b->array[a]).borderindex, b->borderparent[a].i);
	  }*/
} /* check_boundary */


void printfrank(const char *format, int my_rank, ...)
{

  char buffer[MAX_MSG_SIZE] = "";  

    va_list args;

    va_start(args,my_rank);
    vsnprintf(buffer, MAX_MSG_SIZE, format, args);
    va_end(args);   

    if (rank() == my_rank) printf("IN PROC %d, %s\n", my_rank, buffer);
}

void PrintParents(MaxNode *maxtree, size_t width, size_t height, int myrank)
{
  size_t i,j,p;
  printf("MaxTree in proc %d\n\n", myrank);
  for(i=0; i < height; i++)
    {
      for(j = 0; j<width; j++)
	{
	  p =  i*width +j;
	  printf("i:%ld gval:%hi par:%ld bind:%ld \t",p,maxtree[p].gval,maxtree[p].parent, maxtree[p].borderindex);
	}
        printf("\n");
    }
    printf("\n\n");
}


void PrintBound(Boundary *b, int myrank )
{

  if(b != NULL){
    for(size_t i=0; i < b->size; i++)
      {
	printfrank("[BOUNDARY %d], node %ld (node %ld gval %hi area %lu in maxtree %d) has parent %ld in boundary %d, othborderlr is %ld is in boundary %d",myrank,b,i,b->array[i].index,b->array[i].gval,b->array[i].area,b->array[i].process,b->borderparent[i].i, b->borderparent[i].b,b->othborderlr[i].i, b->othborderlr[i].b);// has parent %d in boundary %d (proc %d)
      }
    printf("\n\n");
  }
}
