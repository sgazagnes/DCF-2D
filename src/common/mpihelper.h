#ifndef MPIHELPER_H
#define MPIHELPER_H

int rank(void);
int np(void);
void init_mpi(void);
void finalize_mpi(void);
void create_mpi_maxtree_type(void);
  
#endif
