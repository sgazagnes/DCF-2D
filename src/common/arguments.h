#ifndef ARGUMENTS_H
#define ARGUMENTS_H
#include "types.h"

Arguments *parse_args(int argc, char** argv, Arguments *args);
void print_args(Arguments args);

#endif