#include <math.h>
#include <stddef.h>
#include <stdio.h>
#include <stdbool.h>
#include <mpi.h>

#include "types.h"
#include "arguments.h"
#include "mpihelper.h"
#include "cmdline.h"
#include "logc.h"

Arguments *parse_args(int argc, char** argv, Arguments *args) {
  if (cmdline_parser (argc, argv, args) != 0) {
    error("Error while parsing command-line arguments: %s", *argv);
    MPI_Abort(MPI_COMM_WORLD, 110);
  }

  if(rank()==0 && args->grid_arg[0]*args->grid_arg[1] != np()){
    error("This program assumes an equal amount of processes to the amount of grid tiles");
    MPI_Abort(MPI_COMM_WORLD, 111);
  }
  set_verbosity(args->verbosity_arg);

  return args;
} /* parse_args */


void print_args(Arguments args) {
  info("input prefix: %s",  args.inprefix_arg);
  info("output prefix: %s", args.outprefix_arg);
  info("intype: %s",        args.intype_arg);
  info("outtype: %s",       args.outtype_arg);
  info("grid: %d x %d",     args.grid_arg[0], args.grid_arg[1]);
  info("filter: %s",        args.filter_arg);
  info("lambda: %ld",       args.lambda_arg);
  //info("alloc: %ld",       args.allocation_arg);
  info("verbosity: %s",     args.verbosity_arg);
} /* print_args */
