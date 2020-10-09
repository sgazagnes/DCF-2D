/*
 * Compile using "make area" Makefile, C99 is assumed, needs openmpi or another MPI 2.0 compliant library, and freeimage
 *
 * Author: Simon Gazagnes, using the code previously implemented by Jan Kazemier
 * 
 * This program converts an image file into a max tree in order to filter using connected operators
 * The application supports distributed memoryy parallel machines, through message passing (MPI)
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

long memoryy;
long maxmemoryy;
extern MPI_Datatype 	mpi_max_node_type;
extern int 		verbosity;
struct tms 		tstruct;

int main(int argc, char** argv) {
  /* Initialization */
  init_mpi();
  create_mpi_maxtree_type();
  memoryy =0;
  maxmemoryy = 0;
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
  info("Max grey value: %lu", MAXGREYVAL);
  memoryy += size * sizeof(value);
  maxmemoryy = MAX(memoryy, maxmemoryy);

  
  if(args.flood_orig == NULL){
    if(bitpix >= 16)
      flood_func = 1;
    else
      flood_func = 0;
  }else
    flood_func = args.flood_arg;

  
  if (rank()==0)    start = times(&tstruct);

  /*Build local component tree */
  MaxNode 	*maxtree   = calloc(size, sizeof(MaxNode));
  check_alloc(maxtree, 001);
  memoryy += size * sizeof(MaxNode);
  maxmemoryy = MAX(memoryy, maxmemoryy);
  // printfrank("Init size of mtree: %ld", 0, size*sizeof(MaxNode));
  if(flood_func == 0){
    ulong 	*histogram =  calloc(MAXGREYVAL, sizeof(ulong));
    memoryy += MAXGREYVAL * sizeof(ulong);
    maxmemoryy = MAX(memoryy, maxmemoryy);
    ulong 	min_idx = init_tree_hist(maxtree, size, gvals_buf, histogram);
    Queue 	*queue = create_queue(size, MAXGREYVAL);
    
    memoryy += size * sizeof(ulong) + MAXGREYVAL * sizeof(Queue);
    maxmemoryy = MAX(memoryy, maxmemoryy);
    
    set_queue_offsets(queue, histogram, MAXGREYVAL);
    queue_add(queue, gvals_buf[min_idx], min_idx);
    bool 	*reached = calloc(size, sizeof(bool));
    
    memoryy += size * sizeof(bool);
    maxmemoryy = MAX(memoryy, maxmemoryy);
    reached[min_idx] = true;
    idx 	*levelroots = calloc(MAXGREYVAL, sizeof(idx));
    memoryy += MAXGREYVAL * sizeof(idx);
    maxmemoryy = MAX(memoryy, maxmemoryy);
    
    for (size_t i = 0; i < MAXGREYVAL; i++)    levelroots[i] = BOTTOM;
    levelroots[gvals_buf[min_idx]] = min_idx;
    ulong 	area = 0;
    flood_tree_s(width, height, maxtree, queue, levelroots, reached, gvals_buf[min_idx], &area);
    check_area_size(area, size);
    
    free(reached);
    free(levelroots);
    free(histogram);
    free_queue(queue);
    
  }else{
    ulong 	min_idx = init_tree(maxtree, size, gvals_buf);
    pQueue 	*queue = pQueueCreate(size);
    memoryy += sizeof(pQueue) + (size+1)*sizeof(long);
    maxmemoryy = MAX(memoryy, maxmemoryy);
    
    pStack 	*stack = pStackCreate(size, queue->array);
    memoryy += sizeof(pStack);
    maxmemoryy = MAX(memoryy, maxmemoryy);
   
    flood_tree_w(queue, stack, gvals_buf,  width,  height, maxtree, min_idx);
    pQueueDelete(queue);
    pStackDelete(stack);
  }
     
  free(gvals_buf);
  // printfrank("After tree building, memory used: %ld", 0, maxmemoryy);
  memoryy = size * sizeof(MaxNode);
  maxmemoryy = MAX(memoryy, maxmemoryy);

  /*Correct the local component trees with the boundary trees */  
  if(np() > 1) maxtree = correct_borders(maxtree, width, height, args, &tree_size);

  /* Filter the tree */
  value 	*out = calloc(size , sizeof(value));
  check_alloc(out, 010);
  
  memoryy += size * sizeof(value);
  maxmemoryy = MAX(memoryy, maxmemoryy);
  
  if(!strcmp(args.filter_arg, "area")){
    info("Filtering with following rule: area");
    filter_on_area(size, maxtree, out, args.lambda_arg);
  }
  
  /* Finalize */
  MPI_Barrier(MPI_COMM_WORLD);
  if (rank()==0)
    printf("Distributed memory: wallclock time = %f \n",
	   (float)(times(&tstruct) - start)/(float)sysconf(_SC_CLK_TCK));
  //printfrank("Max memory used %ld \n", 0,maxmemoryy);

  
   WriteOutput(args, out, bitpix, n_dims, dims, dims_T);
  
  /* clean up */
  free(maxtree);
  cmdline_parser_free(&args);
  free(out);
  finalize_mpi();
  return 0;
} /* main */
