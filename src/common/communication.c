
#include <math.h>
#include <stddef.h>
#include <stdio.h>
#include <stdbool.h>
#include <mpi.h>
#include <stdarg.h>
#include "communication.h"
#include "boundary.h"

#include "mpihelper.h"
#include "writefile.h"
#include "cmdline.h"
#include "types.h"
#include "constants.h"
#include "logc.h"
#include "checks.h"

long memoryy;
long maxmemoryy;
extern MPI_Datatype mpi_max_node_type;


void send_boundary(Boundary *b, int dest) {

  /* Prepare parent index*/
  idx *borderparent_indices_only = malloc(b->size * sizeof(idx));
  for (size_t i = 0; i < b->size; i++) 
    borderparent_indices_only[i] = b->borderparent[i].i;
  
  /*Send boundary through MPI */
  int tag = 1; // MPI Tag is ignored 
  int err;
  err = MPI_Send(b->array, b->size, mpi_max_node_type, dest, tag, MPI_COMM_WORLD);
  check_mpi_error(err, 5311);
  err = MPI_Send(b->offset, 5, MPI_LONG, dest, tag, MPI_COMM_WORLD);
  check_mpi_error(err, 5312);
  err = MPI_Send(borderparent_indices_only, b->size, MPI_LONG, dest, tag, MPI_COMM_WORLD);
  check_mpi_error(err, 5313);

  free(borderparent_indices_only);
} /* send_boundary */


Boundary *receive_boundary(int src) {  
  /* Initialize MPI */
  int tag = 1; // MPI Tag is ignored 
  MPI_Status status;

  int recv_size;
  MPI_Probe(src, tag, MPI_COMM_WORLD, &status);
  MPI_Get_count(&status, mpi_max_node_type, &recv_size);

  /* Initialize b */
  Boundary *b = calloc(1, sizeof(Boundary));
  check_alloc(b, 5321);
  memoryy +=  sizeof(Boundary);
  maxmemoryy = MAX(memoryy, maxmemoryy);
 
  b->array = malloc(recv_size * sizeof(MaxNode));
  check_alloc(b->array, 5322);
  memoryy += recv_size * sizeof(MaxNode);
  
  b->offset = calloc(5 , sizeof(size_t));
  check_alloc(b->offset, 5323);
  memoryy += 5*sizeof(size_t);
  
  b->reached = calloc(recv_size, sizeof(bool));
  check_alloc( b->reached, 5325);
  memoryy += recv_size * sizeof(bool);

  b->borderparent = malloc(recv_size * sizeof(BorderIndex));
  check_alloc(b->borderparent, 5326);
  memoryy+= recv_size * sizeof(BorderIndex);

  b->othborderlr = malloc(recv_size * sizeof(BorderIndex));
  check_alloc(b->othborderlr, 5328);
  memoryy+= recv_size * sizeof(BorderIndex);
  maxmemoryy = MAX(memoryy, maxmemoryy);
  
  for (size_t i = 0; i < (size_t) recv_size; ++i) 
    b->othborderlr[i] = (BorderIndex) {.b = b, .i = BOTTOM};
  b->size = recv_size;
  b->initsize = recv_size;
  b->allocsize = recv_size;
  
  idx *borderparent_indices_only = malloc(recv_size * sizeof(idx));
  check_alloc(borderparent_indices_only,5327);
  memoryy+= recv_size * sizeof(idx);
  maxmemoryy = MAX(memoryy, maxmemoryy);
  
  /* Receive from other processes with MPI */
  int err;
  err = MPI_Recv(b->array, recv_size, mpi_max_node_type, src, tag, MPI_COMM_WORLD, &status);
  check_mpi_error(err, 5341);
  err = MPI_Recv(b->offset, 5, MPI_LONG, src, tag, MPI_COMM_WORLD, &status);
  check_mpi_error(err, 5342);
  err = MPI_Recv(borderparent_indices_only, recv_size, MPI_LONG, src, tag, MPI_COMM_WORLD, &status);
  check_mpi_error(err, 5343);

  /* Assign borderparents and allocate otherborderlr */
  for (size_t i = 0; i < b->size; i++) 
    b->borderparent[i] =  (BorderIndex) {.b = b, .i = borderparent_indices_only[i]};
  free(borderparent_indices_only);
  memoryy -= recv_size * sizeof(idx);

  return b;
} /* receive_boundary */


void send_updated_boundary(Boundary *b, int dest) {
  /* Prepare parent index*/
  idx *borderparent_indices_only = malloc(b->size * sizeof(idx));
  check_alloc(borderparent_indices_only, 553);

  ulong *area_only = malloc(b->initsize*sizeof(ulong));
  check_alloc(area_only, 554);

  for (size_t i = 0; i < b->size; i++){ 
    borderparent_indices_only[i] = b->borderparent[i].i;
    if(i<b->initsize)
      area_only[i] = b->array[i].area;
  }
  
  /*Send information through MPI*/
  int tag = 1; // MPI Tag is ignored 
  int err;

  err = MPI_Send(borderparent_indices_only, b->size, MPI_LONG, dest, tag, MPI_COMM_WORLD);
  check_mpi_error(err, 543);
  err = MPI_Send(area_only, b->initsize, MPI_LONG, dest, tag, MPI_COMM_WORLD);
  check_mpi_error(err, 542);


  if(b->size > b->initsize){
    MaxNode *changes = &b->array[b->initsize] ;
    err = MPI_Send(changes, b->size - b->initsize, mpi_max_node_type, dest, tag, MPI_COMM_WORLD);
    check_mpi_error(err, 541);
  }

  free(borderparent_indices_only);
  free(area_only);

} /* send_updated_bobundary */

Boundary *receive_updated_boundary( Boundary *b, int src) {
  /* Initialize MPI */
  int tag = 1; // MPI Tag is ignored 
  MPI_Status status;
  int err;
  int recv_size;
  MPI_Probe(src, tag, MPI_COMM_WORLD, &status);
  MPI_Get_count(&status, MPI_LONG, &recv_size);

  idx *borderparent_indices_only = malloc(recv_size * sizeof(idx));
  check_alloc(borderparent_indices_only, 553);
  
  memoryy += recv_size * sizeof(idx);
   
  ulong *area_only = malloc(b->initsize * sizeof(ulong));
  check_alloc(area_only, 554);
  memoryy += b->initsize * sizeof(ulong);
  
  err = MPI_Recv(borderparent_indices_only, recv_size, MPI_LONG, src, tag, MPI_COMM_WORLD, &status);
  check_mpi_error(err, 556);
  err = MPI_Recv(area_only, b->size, MPI_LONG, src, tag, MPI_COMM_WORLD, &status);
  check_mpi_error(err, 556);
  
  if((size_t) recv_size > b->size){
    /*Reallocation of b to receive the info */
    b->array = realloc(b->array, recv_size * sizeof(MaxNode));
    check_alloc(b->array, 551);
    memoryy -= b->size * sizeof(MaxNode);
    memoryy += recv_size * sizeof(MaxNode);
    maxmemoryy = MAX(memoryy, maxmemoryy);
  
    
    b->borderorigin = realloc(b->borderorigin,recv_size * sizeof(BorderIndex));
    check_alloc(b->borderorigin, 552);
    memoryy -= b->size * sizeof(BorderIndex);
    memoryy += recv_size * sizeof(BorderIndex);
   
    for (size_t i = b->size; i < (size_t) recv_size; ++i) 
      b->borderorigin[i] =   (BorderIndex) {.b = b, .i = BOTTOM};
    b->borderparent = realloc(b->borderparent, recv_size * sizeof(BorderIndex));
    check_alloc(b->borderparent, 554);
    memoryy -= b->size * sizeof(BorderIndex);
    memoryy += recv_size * sizeof(BorderIndex);
    
    b->othborderlr = realloc(b->othborderlr,recv_size * sizeof(BorderIndex));
    check_alloc(b->othborderlr, 552);
    memoryy -= b->size * sizeof(BorderIndex);
    memoryy += recv_size * sizeof(BorderIndex);
    maxmemoryy = MAX(memoryy, maxmemoryy);

    for (size_t i = b->size; i < (size_t) recv_size; ++i) 
      b->othborderlr[i] =   (BorderIndex) {.b = b, .i = BOTTOM};
    b->allocsize = recv_size;

    /* Received supp nodes */
    err = MPI_Recv(b->array+b->size, recv_size-b->size, mpi_max_node_type, src, tag, MPI_COMM_WORLD, &status);
    check_mpi_error(err, 555);

    b->size = recv_size;
  }

  /* Assigning area and parents */

  for (size_t i = 0; i < b->size; i++){
    b->borderparent[i] = (BorderIndex) {.b = b, .i = borderparent_indices_only[i]};
    if(i<b->initsize)
      b->array[i].area = area_only[i] ;
  }
  
  free(borderparent_indices_only);
  memoryy -= recv_size * sizeof(idx);

  free(area_only);
  memoryy -= b->initsize * sizeof(ulong);

  return b;
} /* receive_boundary */
