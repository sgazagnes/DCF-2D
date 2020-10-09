/*
 * Compile using "make area" Makefile, C99 is assumed, needs openmpi, or another MPI 2.0 compliant library
 *
 * Author: Jan J. Kazemier
 * Written as a part of my master thesis
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
  size_t 	combine = 0;
  int	 	bitpix, n_dims;
  size_t 	size, tree_size, width, height;
  size_t 	*dims   = calloc(MAXDIM,sizeof(size_t));
  size_t 	*dims_T = calloc(MAXDIM,sizeof(size_t));
  value 	*gvals_buf = ReadInput(args, &n_dims, dims, dims_T, &bitpix); /* Reading the data */
  clock_t 	start;
  LambdaVec *lvec = LambdaVectorRead(args.lvec_arg, 1.0);

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


  value *outLum = calloc(tree_size, sizeof(value));
  check_alloc(outLum, 021);
  value *outCont = calloc(tree_size, sizeof(value));
  check_alloc(outCont, 022);
  value *outScale = calloc(tree_size, sizeof(value));
  check_alloc(outScale, 023);
  
  ulong *threshold = calloc(lvec->NumLambdas, sizeof(ulong));
  check_alloc(outScale, 024);

  for (int i=0;i<lvec->NumLambdas-1;i++)
    threshold[i]=(ulong) lvec->Lambda[i+1];

  MaxTreeSetDAPseg( size, maxtree, outCont, outScale, outLum, threshold, lvec->NumLambdas -1  );
  free(maxtree);
  
  if(combine){
     info("Combining");
    value *inverted = calloc(size, sizeof(value));
    check_alloc(inverted, 030);
    for(size_t i = 0; i < size ; i++)
      inverted[i] = MAXGREYVAL - 1 - gvals_buf[i];
  
    maxtree   = calloc(size, sizeof(MaxNode));
    check_alloc(maxtree, 001);

    if(flood_func == 0){
      ulong 	*histogram =  calloc(MAXGREYVAL, sizeof(ulong));
      ulong 	min_idx = init_tree_hist(maxtree, size, inverted, histogram);
      Queue 	*queue = create_queue(size, MAXGREYVAL);
      set_queue_offsets(queue, histogram, MAXGREYVAL);
      queue_add(queue, inverted[min_idx], min_idx);
      bool 	*reached = calloc(size, sizeof(bool));
      reached[min_idx] = true;
      idx 	levelroots[MAXGREYVAL];
      for (size_t i = 0; i < MAXGREYVAL; i++)    levelroots[i] = BOTTOM;
      levelroots[gvals_buf[min_idx]] = min_idx;
      ulong 	area = 0;
      flood_tree_s(width, height, maxtree, queue, levelroots, reached, inverted[min_idx], &area);
      check_area_size(area, size);
      free(reached);
      free(histogram);
      free_queue(queue);
    }else{
      ulong 	min_idx = init_tree(maxtree, size, inverted);
      pQueue 	*queue = pQueueCreate(size);
      pStack 	*stack = pStackCreate(size, queue->array);
      flood_tree_w(queue, stack, inverted,  width,  height, maxtree, min_idx);
      pQueueDelete(queue);
      pStackDelete(stack);
    }

    free(inverted);
  
    value *temp;
    size_t tree_size2 = size;
    maxtree = correct_borders(maxtree, width, height, args, &tree_size2);
    value *outLum2 = calloc(tree_size2, sizeof(value));
    check_alloc(outLum2, 024);
    value *outCont2 = calloc(tree_size2, sizeof(value));
    check_alloc(outCont2, 025);
    value *outScale2 = calloc(tree_size2, sizeof(value));
    check_alloc(outScale2, 026);

    temp = outLum;
    outLum = outLum2;
    outLum2 = temp;
    temp =outScale;
    outScale = outScale2;
    outScale2 = temp;
    temp =outCont;
    outCont = outCont2;
    outCont2 = temp;
    MaxTreeSetDAPseg(size, maxtree, outCont, outScale, outLum, threshold,lvec->NumLambdas -1 );

    combineResults(size, maxtree, outCont, outCont2, outScale, outScale2, outLum, outLum2, lvec->NumLambdas -1);

    free(maxtree);
    free(outLum2);
    free(outCont2);
    free(outScale2);
  }

  free(gvals_buf);

  MPI_Barrier(MPI_COMM_WORLD);

  if (rank()==0)
    printf("wallclock time = %f\n",(float)(times(&tstruct) - start)/(float)sysconf(_SC_CLK_TCK));
 
   WriteDAP(  args, outLum, outCont, outScale, bitpix, n_dims, dims_T, dims);
  LambdaVectorDelete(lvec);
  MPI_Barrier(MPI_COMM_WORLD);
     
  /* clean up */
  cmdline_parser_free(&args);
  free(outLum);
  free(outCont);
  free(outScale);
  free(threshold);

  finalize_mpi();
  return 0;
} /* main */
