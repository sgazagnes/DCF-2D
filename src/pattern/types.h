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
  idx index;
  idx parent;
  ulong area;
  value gval;
  short filter;
  ulong privateArea;
  value gval_par;
  idx borderindex; /* index in boundary */
  int process; /* which rank it originally belongs to */
} MaxNode;

typedef struct _BorderIndex BorderIndex;

typedef struct _Boundary {
  MaxNode *array;
  size_t *offset; /* offsets for north, east, south, west and ancestors respectively; north is at always 0, just for convenience */
  size_t size;
  size_t initsize;
  size_t allocsize;
  BorderIndex *borderparent;
  BorderIndex *borderorigin;
  BorderIndex *othborderlr;/* pointer to boundary and index */
  bool *reached; 
} Boundary;

typedef struct _Queue {
  unsigned long *pixels;
  unsigned long head;
  unsigned long tail; /* First free place in queue, empty if head=tail */
} Queue;

struct _BorderIndex {
  Boundary *b;
  /* Boundary *other; */ /* unused, as we overwrite whole BorderIndex */
  idx i;
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
