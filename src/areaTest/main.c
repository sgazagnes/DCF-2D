/*
 * Compile using "make area" Makefile, C99 is assumed, needs openmpi or another MPI 2.0 compliant library, and freeimage
 *
 * Author: Simon Gazagnes, using the code previously implemented by Jan Kazemier
 * 
 * This program converts an image file into a max tree in order to filter using connected operators
 * The application supports distributed memory parallel machines, through message passing (MPI)
 * In this example a filtering is performed using the area of a flat zone as an attribute for filtering
 *
 */


#define _GNU_SOURCE
#include <assert.h>
#include <limits.h>
#include <math.h>
#include <stddef.h>
#include <stdio.h>
#include <errno.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/times.h>
#include <mpi.h>

#include "main.h"
#include "arguments.h"
#include "boundary.h"
#include "checks.h"
#include "constants.h"
#include "image.h"
#include "logc.h"
#include "maxtree.h"
#include "mpihelper.h"
#include "queue.h"
#include "types.h"
#include "writefile.h"

extern MPI_Datatype mpi_max_node_type;
extern int verbosity;

struct tms tstruct;

int main(int argc, char** argv) {
  /* Initialization */
  init_mpi();
  create_mpi_maxtree_type();
  
  Arguments args;
  parse_args(argc, argv, &args);
  print_args(args);
  size_t width, height;
  short bitsperpixel;
  value *gvals_buf = ReadImage(args, &width, &height, &bitsperpixel);
  MAXGREYVAL = pow(2,bitsperpixel);
  size_t size = width * height;
  size_t tree_size = size;
  clock_t start; 
  
  if (rank()==0)
    start = times(&tstruct);

  /*Build local component tree */
  MaxNode *maxtree = create_empty_maxtree(size);
  ulong *histogram =  calloc(MAXGREYVAL, sizeof(ulong));
  idx min_idx = set_gvals(maxtree, size, gvals_buf, histogram);
  value mingval = gvals_buf[min_idx];
  free(gvals_buf);
  Queue *q = create_queue(size, MAXGREYVAL);
  set_queue_offsets(q, histogram, MAXGREYVAL);
  free(histogram);
  queue_add(q, mingval, min_idx);
  bool *reached = calloc(size, sizeof(bool));
  reached[min_idx] = true;
  idx levelroots[MAXGREYVAL];
  for (size_t i = 0; i < MAXGREYVAL; i++)
    levelroots[i] = BOTTOM;
  levelroots[mingval] = min_idx;
  ulong area = 0;
  flood_tree(width, height, maxtree, q, levelroots, reached, mingval, &area);
  check_area_size(area, size);
  free(reached);
  free_queue(q);
  
  /*Correct the local component trees with the boundary trees */  
  maxtree = correct_borders(maxtree, width, height, args, &tree_size);

  /* Filter the tree */
  value *out = calloc(size , sizeof(value));
  check_alloc(out, 010);
  filter_on_area(size, maxtree, out, args.lambda_arg);

  /* Finalize */
  MPI_Barrier(MPI_COMM_WORLD);
  if (rank()==0)
    printf("Distributed memory: wallclock time = %f\n",(float)(times(&tstruct) - start)/(float)sysconf(_SC_CLK_TCK));
  WriteImage(args, out , bitsperpixel, width, height);

  /* clean up */
  free(maxtree);
  cmdline_parser_free(&args);
  free(out);
  finalize_mpi();
  return 0;
} /* main */
