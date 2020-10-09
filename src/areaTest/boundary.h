#ifndef BOUNDARY_H
#define BOUNDARY_H

#include "types.h"

MaxNode *correct_borders(MaxNode *maxtree, size_t width, size_t height, Arguments args, size_t *treesize);
Boundary *update( Boundary *b, Boundary *c, double b_alloc);
Boundary *update_branch(BorderIndex x, BorderIndex z, BorderIndex s);
void update_node_attribute(BorderIndex x, BorderIndex s);
void adding_node(BorderIndex x, BorderIndex s);
Boundary *combine(Boundary *a, Boundary *b, Direction d);
void b_add(Boundary *c,  Boundary *b, size_t s, size_t i);
void reset_borderindex(Boundary *b);
MaxNode *apply_changes(MaxNode *tree, Boundary *b, size_t tree_size);
void merge(Boundary *a, Boundary *b, Direction d);
BorderIndex idx_i(Boundary *b, size_t c, Direction d, size_t length);
BorderIndex idx_j(Boundary *b, size_t c, Direction d);
Boundary *create_boundary(MaxNode *maxtree, size_t width, size_t height, double b_alloc);
void add_side(Boundary *b, MaxNode *maxtree, int side, size_t length, size_t offset, long increment);
Boundary *add_ancestors(Boundary *b, MaxNode *maxtree);
Boundary *realloc_b(Boundary *b, size_t newsize);
void free_boundary(Boundary *b);

BorderIndex b_levelroot(BorderIndex bi);
BorderIndex b_levelrootinbound(BorderIndex bi);
bool is_bottom(BorderIndex bi);
bool bi_equal(BorderIndex ai, BorderIndex bi);
BorderIndex b_parent_lr(BorderIndex bi);
BorderIndex b_parent(BorderIndex bi);
value b_gval(BorderIndex bi);
MaxNode b_node(BorderIndex bi);

#endif
