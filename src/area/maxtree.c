#include <math.h>
#include <stddef.h>
#include <stdio.h>
#include <stdbool.h>
#include <mpi.h>
#include <assert.h>

#include "maxtree.h"
#include "arguments.h"
#include "cmdline.h"
#include "types.h"
#include "constants.h"
#include "logc.h"
#include "checks.h"
#include "queue.h"
#include "mpihelper.h"


ulong init_tree(MaxNode *tree, size_t size, value *gvals) {
  ulong min_idx = 0;
  for (size_t x = 0; x < size; x++) {
    tree[x].index = x;
    tree[x].parent = BOTTOM;
    tree[x].area = 0;
    tree[x].gval = gvals[x];
    tree[x].filter = BOTTOM;
    tree[x].borderindex = BOTTOM;
    tree[x].process = rank();
    if(gvals[x] < gvals[min_idx]) min_idx = x;
  }
  return min_idx;
} /* set_to_zero */

ulong init_tree_hist(MaxNode *tree, size_t size, value *gvals, ulong *hist) {
  ulong min_idx = 0;  
  for (size_t x = 0; x < size; x++) {
    tree[x].index = x;
    tree[x].parent = BOTTOM;
    tree[x].area = 0;
    tree[x].gval = gvals[x];
    tree[x].filter = BOTTOM;
    tree[x].borderindex = BOTTOM;
    tree[x].process = rank();
    hist[gvals[x]]++;
    if(gvals[x] < gvals[min_idx]) min_idx = x;
  }
    return min_idx;
} /* set_to_zero */


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

uint get_neighbors (ulong *neighbors, ulong p, ulong x, ulong y, size_t width, size_t height) {
  uint n = 0;

  if (x < (width - 1))     neighbors[n++] = p + 1;
  if (y > 0)               neighbors[n++] = p - width;
  if (x > 0)               neighbors[n++] = p - 1;
  if (y < (height - 1))    neighbors[n++] = p + width;

  return n;
} /* get_neighbors */

int flood_tree_s(size_t width, size_t height, MaxNode *maxtree, Queue *q, idx *levelroot, bool *reached, value level, ulong *thisarea) {
  /* flood gray level level */
  ulong neighbors[CONNECTIVITY]; /* CONNECTIVITY is set in constants.h */
  ulong area = *thisarea;
  ulong p,x,y;
  /* propagation */
  while (queue_is_not_empty(q, level)) {
    area++;
    p = queue_first(q, level); /* current pixel */
    x = p % width; /* x coordinate of p */
    y = p / width; /* y coordinate of p */


    uint numneighbors = get_neighbors(neighbors, p, x, y, width, height);

    for (size_t i = 0; i < numneighbors; i++) {
      idx c = neighbors[i]; /* current neighbor */
      /* trace("%zu i: %zu reached: %d", p, i, (maxtree[c].flags & REACHED)); */
      if (!reached[c]) {
        reached[c] = true;
        value fc = maxtree[c].gval; // image->gval[c]; /* grey value of current neighbour */
        /*trace("%zu i: %zu fc: %d", p, i, fc);*/
        if (levelroot[fc] == BOTTOM) { /* current grey value still has bottom as levelroot */
          levelroot[fc] = c; /* now set to current pixel */
        } else {
          maxtree[c].parent = levelroot[fc]; /* already has level root, set parent to it */
        }

        queue_add(q, fc, c);

        if (fc > level) { /* current value is larger than current level */
          ulong childarea = 0;
          do {
            fc = flood_tree_s(width, height, maxtree, q, levelroot, reached, fc, &childarea);
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


void flood_tree_w(pQueue *queue, pStack *stack, value *gval, size_t width, size_t height, MaxNode *node, ulong min_idx){
  ulong neighbors[CONNECTIVITY];
  ulong i, nextpix, p, q, oldtop, oldpix;
  uint numneighbors;
  ulong x,y;

  node[min_idx].parent = min_idx;
  pStackPush(stack, min_idx);
  node[min_idx].area = 1;
 
  nextpix = min_idx;

  do{ 
    p = nextpix;
    assert( ( p == pStackTop(stack) ) || ( p == pQueueFront(queue) ) );
    oldpix = p;
    x = p % width; 
    y = p / width;

    numneighbors = get_neighbors(neighbors, p, x, y, width, height);    

    for (i=0; i<numneighbors; i++) {
      q = neighbors[i];
      if (node[q].parent == BOTTOM){
	node[q].parent=q;
	node[q].area=1;
        if (gval[q]>gval[p]){
	  pStackPush(stack,q);
	  nextpix = q;

          assert (nextpix == pStackTop(stack) );
          assert (nextpix != oldpix );
	  break;
	}
	pQueuePush(queue,gval,q);
      }
    }

    if (nextpix==p){      /* No break occurred */
      if (p != pStackTop(stack)){      /* p is in queue and 
				          processing is finished  */
	assert(gval[p] == gval[pStackTop(stack)]);

	p =  pQueuePop(queue,gval);     /* remove from queue */
	node[p].parent = pStackTop(stack);
	node[pStackTop(stack)].area++;

	if (!IsEmpty(queue)){
	  nextpix = pQueueFront(queue);

	  if (gval[nextpix] < gval[p]){   
	    /* moving down, but first process top of stack */	    
	    nextpix = pStackTop(stack);
	  }
	} else {
	    nextpix = pStackTop(stack);
	}
   
	assert( ( nextpix != oldpix ) &&
	        (( nextpix == pStackTop(stack) ) || 
		 ( nextpix == pQueueFront(queue) ) ));
        	
      } else {                          /* p is top of stack */

	if (!IsEmpty(queue)){	  
	  nextpix = pQueueFront(queue);    /* candidate next pixel */
	  assert( oldpix != nextpix );

	  if (gval[nextpix] < gval[p]){ /* moving down, top of stack done */
	    oldtop = pStackPop(stack);
	    p = pStackTop(stack);
	    assert( oldpix != p);
 
	    if (gval[nextpix] > gval[p]){
	      nextpix = pQueuePop(queue, gval); /* move nextpix from 					     queue to stack     */
	      pStackPush(stack, nextpix);
              p = nextpix; 
	      assert( oldpix != nextpix );
        
	    } else if (gval[nextpix] < gval[p]){
              nextpix = p;              /* flood from current stack top */
	    }
	    node[oldtop].parent = p;           /* == stack top */
	    node[p].area += node[oldtop].area; /* update area */
           
	  }
	  assert( oldpix != nextpix );
	} else {

	    oldtop = pStackPop(stack);
            if (!IsEmptyStack(stack)){
	      p = pStackTop(stack);
	      node[oldtop].parent = p;           /* == stack top */
	      node[p].area += node[oldtop].area; /* update area */
	      nextpix = p;
	      assert( oldpix != nextpix );
	    }
	}
	
      }

    }

  } while (!IsEmpty(queue) || (!IsEmptyStack(stack)));   
  node[min_idx].parent=BOTTOM;
}



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
