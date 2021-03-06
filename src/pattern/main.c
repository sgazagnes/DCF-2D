
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
#include "lambdavec.h"
#include "maxtree.h"
#include "mpihelper.h"
#include "queue.h"
#include "types.h"
#include "writefile.h"

extern MPI_Datatype mpi_max_node_type;
extern int verbosity;

struct tms tstruct;

int main(int argc, char** argv) {
  init_mpi();
  create_mpi_maxtree_type();
  Arguments 	args;
  parse_args(argc, argv, &args);
  print_args(args);

  bool 		flood_func;
  int	 	bitpix, n_dims;
  size_t 	size, tree_size, width, height;
  size_t 	*dims   = calloc(MAXDIM,sizeof(size_t));
  size_t 	*dims_T = calloc(MAXDIM,sizeof(size_t));
  value 	*gvals_buf = ReadInput(args, &n_dims, dims, dims_T, &bitpix); /* Reading the data */
  clock_t 	start;

  width = dims[0];
  height = dims[1];
  size = tree_size = width * height;
  MAXGREYVAL = pow(2,bitpix);
  info("Max grey value %ld", MAXGREYVAL);

  if(args.flood_orig == NULL){
    if(bitpix >= 16)
      flood_func = 1;
    else
      flood_func = 0;
  }else
    flood_func = args.flood_arg;
   
  if (rank()==0)    start = times(&tstruct);
  LambdaVec *lvec = LambdaVectorRead(args.lvec_arg, 1.0);

 /*Build local component tree */
  MaxNode 	*maxtree   = calloc(size, sizeof(MaxNode));
  check_alloc(maxtree, 001);

  if(flood_func == 0){
    ulong 	*histogram =  calloc(MAXGREYVAL, sizeof(ulong));
    ulong 	min_idx = init_tree_hist(maxtree, size, gvals_buf, histogram);
    Queue 	*queue = create_queue(size, MAXGREYVAL);
    set_queue_offsets(queue, histogram, MAXGREYVAL);
    queue_add(queue, gvals_buf[min_idx], min_idx);
    bool 	*reached = calloc(size, sizeof(bool));
    reached[min_idx] = true;
    idx 	levelroots[MAXGREYVAL];
    for (size_t i = 0; i < MAXGREYVAL; i++)    levelroots[i] = BOTTOM;
    levelroots[gvals_buf[min_idx]] = min_idx;
    ulong 	area = 0;
    flood_tree_s(width, height, maxtree, queue, levelroots, reached, gvals_buf[min_idx], &area);
    check_area_size(area, size);
    free(reached);
    free(histogram);
    free_queue(queue);
  }else{
    ulong 	min_idx = init_tree(maxtree, size, gvals_buf);
    pQueue 	*queue = pQueueCreate(size);
    pStack 	*stack = pStackCreate(size, queue->array);
    flood_tree_w(queue, stack, gvals_buf,  width,  height, maxtree, min_idx);
    pQueueDelete(queue);
    pStackDelete(stack);
  }
  
  if(np() > 1) maxtree = correct_borders(maxtree, width, height, args, &tree_size);

  ulong *threshold = calloc(lvec->NumLambdas, sizeof(ulong));
  ulong *spectrum = calloc(lvec->NumLambdas, sizeof(ulong));
  for (int i=0;i<lvec->NumLambdas-1;i++)
    threshold[i]=(ulong) lvec->Lambda[i+1];
  // info("Before pattern %f",(float)(times(&tstruct) - start)/(float)sysconf(_SC_CLK_TCK));
  areaPatternSpectrum(size, maxtree, threshold, lvec->NumLambdas-1, spectrum);
  //info("adding specy %f",(float)(times(&tstruct) - start)/(float)sysconf(_SC_CLK_TCK));
   if(np() > 1) add_spec(args, spectrum,  lvec->NumLambdas-1);
   // info("free mtree %f",(float)(times(&tstruct) - start)/(float)sysconf(_SC_CLK_TCK));
  free(maxtree);

  MPI_Barrier(MPI_COMM_WORLD);
  //info("waiting");
  if (rank()==0)
    printf("wallclock time = %f\n",(float)(times(&tstruct) - start)/(float)sysconf(_SC_CLK_TCK));

  if(rank() == 0){
    FILE *f = fopen("pattern.txt", "w");
    if (f == NULL) {
      printf("Error opening file!\n");
      exit(1);
    }
    for(size_t i = 0; i< (size_t) lvec->NumLambdas; i++){
      fprintf(f, "Spec[%ld]= %ld \n", i, spectrum[i]);
    }
    fclose(f);
  }

  free(spectrum);
  LambdaVectorDelete(lvec);

  /* clean up */
  cmdline_parser_free(&args);
  free(threshold);

  finalize_mpi();
  return 0;
} /* main */
