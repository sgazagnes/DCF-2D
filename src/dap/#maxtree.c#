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

#define ISROOT(a) ((a) == BOTTOM) 
#define ISLEVELROOT(node, v)  ((ISROOT(node[v].parent)) || (node[v].gval!=node[node[v].parent].gval)) 

ulong init_tree(MaxNode *tree, size_t size, value *gvals) {
  ulong min_idx = 0;
  for (size_t x = 0; x < size; x++) {
    tree[x].index = x;
    tree[x].parent = BOTTOM;
    tree[x].area = 0;
    tree[x].curDH = 0;
    tree[x].scale = 0;
    tree[x].valid = false;
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
    tree[x].curDH = 0;
    tree[x].scale = 0;
    tree[x].valid = false;
    tree[x].gval = gvals[x];
    tree[x].filter = BOTTOM;
    tree[x].borderindex = BOTTOM;
    tree[x].process = rank();
    hist[gvals[x]]++;
    if(gvals[x] < gvals[min_idx]) min_idx = x;
  }
    return min_idx;
} /* set_to_zero */



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
        if (levelroot[fc] == (idx)BOTTOM) { /* current grey value still has bottom as levelroot */
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

void combineResults(size_t size, MaxNode *maxtree, value *outDH, value *outDH2, value *outScale, value *outScale2, value* outOrig, value* outOrig2, int numscales) {
  
  for (size_t p = 0; p<size; p++){
    if (outDH2[p]>outDH[p]) { /* maxDH of opening > maxDH of closing instance */
      outScale[p] = numscales + 1 + outScale2[p];
      outDH[p] = outDH2[p];
      outOrig[p] = outOrig2[p];
    } else if (outDH2[p]<outDH[p]){
      outOrig[p] = MAXGREYVAL -1  -  outOrig[p];
      outScale[p]++;            /* to match convention of Ispra code */
    } else { /* equality */      
      outScale[p]=0;
      outOrig[p] = (outDH[p]!=0) * (MAXGREYVAL-1 - maxtree[p].gval);
    }    
  }
}/* combineResults */

int findScale(ulong area, ulong *lambda, int numscales ){
  int upper = numscales-1, lower = 0, mid;
  
  //if (area <lambda[lower])
  //return lower;
  if (area >= lambda[upper])
    return upper+1;

  mid = (upper + lower) / 2;
  while(mid!=lower){
    if (area >= lambda[mid])
      lower = mid;
    else
      upper = mid;
    mid = (upper + lower) / 2;
  }
  return lower;
} /* findScale */

void NodeSetDAPseg( MaxNode * node, ulong current, 
                    value *curDH,
		    value *maxDH, value *maxScale, 
		    value *maxOrig, value *curScale,  
		    value *outDH, value *outScale, 
		    value *outOrig, ulong *lambda, 
		    int numscales      ){
  /* pre: node[current].valid == false */
  ulong parent; 
  value scale, DH;

  if (ISLEVELROOT(node, current)){
    scale = findScale(node[current].area,lambda,numscales);
    DH = ((!ISROOT(node[current].parent)) && (scale < numscales)) ?
      node[current].gval - node[node[current].parent].gval : 0;
  }
  if ( ISLEVELROOT(node, current) && 
       ( ( scale == numscales) || 
         ( ISROOT(node[current].parent) )       )  ){

    if ( !ISROOT(node[current].parent) ) {
      *maxScale = numscales;
      *maxOrig = 0;
    } else {
      *maxScale = scale;
      *maxOrig = node[current].gval;
    }

    *maxDH = 0;
    *curDH = 0;
    *curScale = *maxScale;

  } else {

    parent = node[current].parent;
    
    if (!node[parent].valid) { 
      /* go into recursion to set parent values correctly */

      NodeSetDAPseg( node, parent, 
		     curDH, maxDH, maxScale, 
		     maxOrig, curScale,  
		     outDH, outScale, 
		     outOrig, lambda,numscales);

    } else { /* if the parent is valid, copy relevant values */

      *maxScale = outScale[parent];
      *maxDH = outDH[parent];
      *maxOrig = outOrig[parent];
      *curScale = node[parent].scale;
      *curDH = node[parent].curDH;

    } 
    
    if ( ISLEVELROOT(node, current) ){  
      /* if I have a level root, some things might change */
      
      if ( scale == *curScale) {
        /* parent's area is in same scale class, 
	   add current pixel's curDH */
	*curDH += DH;
      } else {
	/* at scale class change, update current scale and DH */
	*curDH = DH;
	*curScale = scale;
      }
      
      if ( *curDH >= *maxDH ) {
	/* If updated curDH is higher than or equal to the maximum DH found
	   update maxDH, maxScale, and outOrig */
	*maxDH = *curDH;
	*maxScale = *curScale;
	*maxOrig = node[current].gval;
      } 
      
    } 
 
  }

  outScale[current] = *maxScale;
  outDH[current] = *maxDH;
  outOrig[current] = *maxOrig;      
  node[current].scale = *curScale;
  node[current].curDH = *curDH;
  node[current].valid = true;
  
  return;
} /*NodesetDAPseg */

void MaxTreeSetDAPseg( size_t size, MaxNode * node, value *outDH, value *outScale, 
		       value *outOrig, ulong *lambda, 
		       int numscales     ){
  uint v;
  value curScale, maxScale;
  value maxDH,  curDH, maxOrig;

  for (v=0; v<size; v++) {
    if (!node[v].valid) {
      NodeSetDAPseg(node, v,
		    &curDH, &maxDH, 
		    &maxScale, &maxOrig, 
		    &curScale,  outDH, 
		    outScale, outOrig,lambda,numscales);
    }
  }
} /* MaxTreeSetDAPseg */


