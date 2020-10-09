#ifndef MAXTREE_H
#define MAXTREE_H

#include "types.h"
#include "queue.h"

ulong init_tree(MaxNode *tree, size_t size, value *gvals);
ulong init_tree_hist(MaxNode *tree, size_t size, value *gvals, ulong *hist);
void filter_on_area(size_t size, MaxNode *maxtree, value *out, ulong lambda);
uint get_neighbors (ulong *neighbors, ulong p, ulong x, ulong y, size_t width, size_t height);
int flood_tree_s(size_t width, size_t height, MaxNode *maxtree, Queue *q, idx *levelroot, bool *reached, value level, ulong *thisarea);
void flood_tree_w(pQueue *queue, pStack *stack, value *gval, size_t width, size_t height, MaxNode *node, ulong min_idx);
idx get_levelroot(MaxNode *maxtree, idx x);
idx levelroot(MaxNode *maxtree, MaxNode this);
idx get_parent(MaxNode *maxtree, idx x);

#endif
