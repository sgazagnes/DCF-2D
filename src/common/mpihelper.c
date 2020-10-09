#include <math.h>
#include <stddef.h>
#include <stdio.h>
#include <stdbool.h>
#include <mpi.h>

#include "mpihelper.h"
#include "types.h"

MPI_Datatype mpi_max_node_type;

int rank(void) {
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  return rank;
} /* rank */

int np(void) {
  int numprocs;
  MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
  return numprocs;
} /* np */


void init_mpi(void) {
  int initialized;

  MPI_Initialized(&initialized);
  if (!initialized) {
    MPI_Init(NULL, NULL);
  }
} /* init_mpi */

void finalize_mpi(void) {
  int finalized;
  MPI_Finalized(&finalized);
  if (!finalized) {
    MPI_Finalize();
  }
} /* finalize_mpi */


void create_mpi_maxtree_type(void) {
  /* Create MPI_Type for MaxTree Nodes */
  const int nitems = 7;
  int blocklengths[7] = {1, 1, 1, 1, 1, 1, 1};
  MPI_Aint offsets[7];
  offsets[0] = offsetof(MaxNode, index);
  offsets[1] = offsetof(MaxNode, parent);
  offsets[2] = offsetof(MaxNode, area);
  offsets[3] = offsetof(MaxNode, gval);
  offsets[4] = offsetof(MaxNode, filter);
  // offsets[5] = offsetof(MaxNode, reached);
  offsets[5] = offsetof(MaxNode, borderindex);
  offsets[6] = offsetof(MaxNode, process);

  /* skipped global index here */

  MPI_Datatype types[8] = {
    MPI_LONG,
    MPI_LONG,
    MPI_UNSIGNED_LONG,
    MPI_INT,
    MPI_INT,
    //  MPI_BOOL,
    MPI_LONG,
    MPI_INT
  };

  MPI_Type_create_struct(nitems, blocklengths, offsets, types, &mpi_max_node_type);
  MPI_Type_commit(&mpi_max_node_type);
} /* create_mpi_maxtree_type */
