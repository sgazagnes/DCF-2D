#ifndef MAXTREE_H
#define MAXTREE_H

#include "types.h"
#include "queue.h"

MaxNode *create_empty_maxtree(size_t size);
void set_to_zero(MaxNode *tree, size_t size);
idx set_gvals(MaxNode *tree, size_t size, value *gvals, ulong *hist);
void filter_on_area(size_t size, MaxNode *maxtree, value *out, ulong lambda);
size_t get_neighbors (idx *neighbors, idx p, ulong x, ulong y, ulong width, ulong height);
int flood_tree(size_t width, size_t height, MaxNode *maxtree, Queue *q, idx *levelroot, bool *reached, int level, ulong *thisarea);
idx get_levelroot(MaxNode *maxtree, idx x);
idx levelroot(MaxNode *maxtree, MaxNode this);
idx get_parent(MaxNode *maxtree, idx x);

#endif
