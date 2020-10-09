#ifndef TYPES_H
#define TYPES_H
#include <stdlib.h>
#include <stdbool.h>

#define MIN(a,b)  ((a<=b) ? (a) : (b))
#define MAX(a,b)  ((a>=b) ? (a) : (b))

typedef struct gengetopt_args_info Arguments;
typedef unsigned long  ulong;
typedef long idx;
typedef unsigned short value;
typedef unsigned char ubyte;

typedef struct _MaxNode {
  idx 	index;
  idx 	parent;
  ulong area;
  value gval;
  short filter; /* gray level after filtering */
  idx 	borderindex; /* index in boundary */
  int 	process; /* which rank it originally belongs to */
} MaxNode;

typedef struct _BorderIndex BorderIndex;

typedef struct _Boundary {
  MaxNode 	*array;
  size_t 	*offset; /* offsets for north, east, south, west and ancestors respectively; */
  size_t 	size;
  size_t 	initsize;
  size_t 	allocsize;
  BorderIndex 	*borderparent; /* Parent in the boundary tree */
  BorderIndex 	*borderorigin; /* Keep track of the original boundary tree node used in the combined boundary tree during the merge/combining stage */
  BorderIndex 	*othborderlr; /* Keep track of the levelroot in each boundary tree during the merge stage */
  bool 		*reached;
} Boundary;

typedef struct _Queue {
  ulong *pixels;
  ulong head;
  ulong tail; /* First free place in queue, empty if head=tail */
} Queue;

struct _BorderIndex {
  Boundary 	*b;
  idx 		i;
};

typedef struct {
  ulong  curpos;
  ulong  maxsize;
  ulong *array;  
} pStack;


typedef struct _pQueue {
  ulong size;
  ulong maxsize;
  ulong *array;
} pQueue;


typedef enum {
  HORIZONTAL = 0,
  VERTICAL = 1
} Direction;


#endif
