#ifndef WRITEFILE_H
#define WRITEFILE_H

#include <stdlib.h>
#include "types.h"

void write_boundary_file_ascii(Boundary *b, const char *filename);
void write_boundary_file_binary(Boundary *b, const char *filename);
void write_maxtree_file_ascii(MaxNode *b, size_t size, const char *filename);

void queues_to_file(Queue *hq, value maxvalue, const char *filename);

#endif