#include <math.h>
#include <stddef.h>
#include <stdio.h>
#include <stdbool.h>
#include <mpi.h>

#include "maxtree.h"
#include "arguments.h"
#include "cmdline.h"
#include "types.h"
#include "constants.h"
#include "logc.h"
#include "checks.h"
#include "queue.h"
#include "mpihelper.h"


MaxNode *create_empty_maxtree(size_t size) {
  MaxNode *tree = calloc(size, sizeof(MaxNode));
  check_alloc(tree, 310);
  set_to_zero(tree, size);
  return tree;
} /* create_empty_maxtree */

void set_to_zero(MaxNode *tree, size_t size) {
  for (size_t x = 0; x < size; x++) {
    tree[x].index = x;
    tree[x].parent = BOTTOM;
    tree[x].area = 0;
    tree[x].gval = 0;
    tree[x].filter = BOTTOM;
    tree[x].borderindex = BOTTOM;
    tree[x].process = rank();
  }
} /* set_to_zero */

idx set_gvals(MaxNode *tree, size_t size, value *gvals, ulong *hist) {
  /* sets gvals in tree, and also returns index of lowest intensity: min_gval_idx */

  value min_gval = SHRT_MAX;
  idx min_gval_idx = LONG_MAX;

  for (size_t x = 0; x < size; x++) {
    value gval = gvals[x]; /* get value from buffer */
    tree[x].gval = gval; /* set value in tree */
    hist[gval]++; /* increment histogram */
    if (gval < min_gval) { /* check for minimum */
      min_gval = gval; /* new minimum */
      min_gval_idx = x; /* new index */
    }
  }
  return min_gval_idx;
} /* set_gvals */

void filter_on_area(size_t size, MaxNode *maxtree, value *out, ulong lambda) {
  int val;
  for (size_t v = 0; v < size; v++) {
    if (maxtree[v].filter == BOTTOM) { /* not filtered yet */
      idx w = (idx)v;
      idx parent = maxtree[w].parent;
      /* repeat while we're not at the bottom, th */
      while ((parent != BOTTOM) && (maxtree[w].filter == BOTTOM) && ((maxtree[w].gval == maxtree[parent].gval ) || (maxtree[w].area < lambda))) {
        w = parent;
        parent = maxtree[w].parent;
      }

      if (maxtree[w].filter != BOTTOM) {
        /* criterion satisfied at level maxtree[w].filter */
        val = maxtree[w].filter;
      } else if (maxtree[w].area >= lambda) {
        /* w satisfies criterion */
        val = maxtree[w].gval; // image->gval[w];
      } else {
        /* criterion cannot be satisfied */
        val = 0;
      }

      /* set filt along par-path from v to w */
      idx u = v;
      while (u != w) {
        /* if(0<=u && u<size){ is always true */
        maxtree[u].filter = val;
        u = maxtree[u].parent;

      }
      /* if(0<=w && w<size){ is always true */
      maxtree[w].filter = val;
    }
    out[v] = maxtree[v].filter;
  }
} /* filter_on_area */

size_t get_neighbors (idx *neighbors, idx p, ulong x, ulong y, ulong width, ulong height) {
  int n = 0;

  if (x < (width - 1))     neighbors[n++] = p + 1;
  if (y > 0)               neighbors[n++] = p - width;
  if (x > 0)               neighbors[n++] = p - 1;
  if (y < (height - 1))    neighbors[n++] = p + width;

  return n;
} /* get_neighbors */

int flood_tree(size_t width, size_t height, MaxNode *maxtree, Queue *q, idx *levelroot, bool *reached, int level, ulong *thisarea) {
  /* flood gray level level */
  idx neighbors[CONNECTIVITY]; /* CONNECTIVITY is set in constants.h */
  ulong area = *thisarea;

  /* propagation */
  while (queue_is_not_empty(q, level)) {
    area++;
    idx p = queue_first(q, level); /* current pixel */

    idx x = p % width; /* x coordinate of p */
    idx y = p / width; /* y coordinate of p */

    /* 3D: y = (p % size) / image->width; */
    /* 3D: z = p / size; */

    size_t numneighbors = get_neighbors(neighbors, p, x, y, width, height);

    for (size_t i = 0; i < numneighbors; i++) {
      idx c = neighbors[i]; /* current neighbor */
      /* trace("%zu i: %zu reached: %d", p, i, (maxtree[c].flags & REACHED)); */
      if (!reached[c]) {
        reached[c] = true;
        value fc = maxtree[c].gval; // image->gval[c]; /* grey value of current neighbour */
        /*trace("%zu i: %zu fc: %d", p, i, fc);*/
        if (levelroot[fc] == (idx)BOTTOM) { /* current grey value still has bottom as levelroot */
          levelroot[fc] = c; /* now set to current pixel */
        } else {
          maxtree[c].parent = levelroot[fc]; /* already has level root, set parent to it */
        }

        queue_add(q, fc, c);

        if (fc > level) { /* current value is larger than current level */
          ulong childarea = 0;
          do {
            fc = flood_tree(width, height, maxtree, q, levelroot, reached, fc, &childarea);
            if (fc >= MAXGREYVAL) { /* TODO: reinstate MAXLEVELS? */
              return fc;
            }
          } while (fc != level);
          area += childarea;
        }
      }
    }
  }
  /* let m := the highest level that has a level root */
  int m = level - 1;
  while (m > 0 && levelroot[m] == BOTTOM) m--; /* find parent */
  if (m >= 0) maxtree[levelroot[level]].parent = levelroot[m]; /* connect its parent */
  maxtree[levelroot[level]].area = area;
  levelroot[level] = BOTTOM;
  *thisarea = area;
  return m;
} /* flood_tree */

idx get_levelroot(MaxNode *maxtree, idx x) {
  /* index based, check whether this node is bottom */
  idx r = x;
  if (r == BOTTOM) return BOTTOM;
  trace("get_levelroot %ld %p", x, (void*)maxtree);
  value gv = maxtree[x].gval; // image->gval[x];
  while ((maxtree[r].parent != BOTTOM) && (gv == maxtree[maxtree[r].parent].gval)) { // image->gval[maxtree[r].parent]
    trace("going up r:%ld gv: %d", r, gv);
    r = maxtree[r].parent;
  }
  /* tree compression */
  while (x != r) {
    trace("compressing");
    idx y = maxtree[x].parent;
    maxtree[x].parent = r;
    x = y;
  }
  return r;
} /* get_levelroot */

idx levelroot(MaxNode *maxtree, MaxNode this) {
  return get_levelroot(maxtree, this.index);
} /* levelroot */

idx get_parent(MaxNode *maxtree, idx x) {
  /* returns index of level root of the parent */
  return get_levelroot(maxtree, maxtree[x].parent);
} /* get_parent */
