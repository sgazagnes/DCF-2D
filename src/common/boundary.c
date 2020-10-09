#include <math.h>
#include <stddef.h>
#include <stdio.h>
#include <stdbool.h>
#include <mpi.h>

#include "boundary.h"
#include "communication.h"
#include "mpihelper.h"
#include "writefile.h"
#include "cmdline.h"
#include "types.h"
#include "constants.h"
#include "logc.h"
#include "checks.h"

extern MPI_Datatype mpi_max_node_type;

MaxNode *correct_borders(MaxNode *maxtree, size_t width, size_t height, Arguments args, size_t* tree_size) {
  
  /* Initialize variables */
  int merged = 0;
  int base_w = 1;
  int base_h = 1;
  int g_w = args.grid_arg[0];
  int g_h = args.grid_arg[1];
  int grid_w = g_w;
  int grid_h = g_h;
  int my_rank = rank();
  int my_row = my_rank / grid_w;
  int connect = 0;
  int connect_nb = 0;
  

  /* Find the size of the boundary tree array for this process */
  if(my_rank % 2)
    connect = 0;
  else if(my_row % 2)
    connect = 1;
  else{
    base_w = g_w - my_rank % g_w;
    base_h = g_h - my_row % g_h;
    while(!(base_w % 2)){
      connect++;
      base_w /= 2;
    }
    while(!(base_h %2)){
      connect++;
      base_h /= 2;
    }
    base_w = 1; 
    base_h = 1;
  }  
  Boundary **b=  malloc((1+2*connect)*sizeof(Boundary*));
  long memoryybef = memoryy;
  /* Initialize the boundary tree of the local component tree */
  b[0] = create_boundary(maxtree, width, height, grid_w, grid_h);      
  long memoryybound = memoryy - memoryybef;
  MPI_Allreduce(&memoryybound, &memoryybound, 1, MPI_LONG, MPI_MAX, MPI_COMM_WORLD);
  // printfrank("Size of boundtree %ld", 0, memoryybound);

  /* Start Merging the boundary trees */
  while (grid_w > 1 || grid_h > 1) {   
    if (grid_w > 1 && grid_w >= grid_h) {  // MERGE HORIZONTALLY       
      if (!merged) {
        if (my_rank % (2 * base_w) == 0) {
          int other_rank = my_rank + base_w;
          b[connect_nb+1] = receive_boundary(other_rank);
	  merge(b[connect_nb], b[connect_nb+1], HORIZONTAL);
	  b[connect_nb+2] = combine(b[connect_nb], b[connect_nb+1], HORIZONTAL);
	  connect_nb += 2;
	}
	else {
	  int other_rank = my_rank - base_w;
	  send_boundary(b[connect_nb], other_rank);
	  merged++;
	}
      }
      else {
	merged++;
      }
      base_w *= 2;
      grid_w = grid_w / 2;
    }

    if (grid_h > 1 && grid_h >= grid_w) {  // MERGE VERTICALLY
      if (!merged) {
        if (my_row % (2 * base_h) == 0) {
          int other_rank = my_rank + base_h * g_w;
          b[connect_nb+1] = receive_boundary(other_rank);
          merge(b[connect_nb], b[connect_nb+1], VERTICAL);
          b[connect_nb+2] = combine(b[connect_nb], b[connect_nb+1], VERTICAL);
	  connect_nb += 2;
	} else {
	  int other_rank = my_rank - base_h * g_w;	  
          send_boundary(b[connect_nb], other_rank);
          merged++;
        }
      } else {
	merged++;
      }
      base_h *= 2;
      grid_h = grid_h / 2;    }    
  }
  base_h /= 2;
  base_w /= 2;
  merged--;

  /* Updating the boundary trees */
  while (grid_h != g_h || grid_w != g_w) {
    if (base_h >= base_w) {   // Update vertically
      if (merged<=0) {
	if (my_row % (2 * base_h) == 0) {
	  int other_rank = my_rank + base_h * g_w;
	  b[connect_nb-2] = update(b[connect_nb-2],  b[connect_nb]);
	  if( b[connect_nb-2]->size > b[connect_nb-2]->initsize){
	    b[connect_nb-2]->borderorigin = realloc( b[connect_nb-2]->borderorigin, b[connect_nb-2]->size * sizeof(BorderIndex));
	    memoryy -= b[connect_nb-2]->initsize * sizeof(BorderIndex);
	    memoryy += b[connect_nb-2]->size * sizeof(BorderIndex);
	    maxmemoryy = MAX(memoryy, maxmemoryy);
  
	    check_alloc(b[connect_nb-2]->borderorigin, 4610);
	    for(size_t i =b[connect_nb-2]->initsize; i<b[connect_nb-2]->size; i++)
	      b[connect_nb-2]->borderorigin[i] = (BorderIndex) {.b = b[connect_nb-2], .i = BOTTOM};
	  }
	  b[connect_nb-1] = update(b[connect_nb-1],  b[connect_nb]);
	  free_boundary(b[connect_nb]);
	  send_updated_boundary(b[connect_nb-1], other_rank);
	  free_boundary(b[connect_nb-1]);
      	  connect_nb -= 2;
	} else {
	  int other_rank = my_rank - base_h * g_w;
	  b[connect_nb] = receive_updated_boundary(b[connect_nb], other_rank);
	  merged--;
	}
      }
      else {
	merged--;
      }
      base_h /= 2;
      grid_h *= 2;
    }     
    if (base_w >= base_h && base_w != 0) { // Update horizontally
      if (merged<=0) {
	if (my_rank % (2 * base_w) == 0) {
	  int other_rank = my_rank + base_w;
	  b[connect_nb-2] = update(b[connect_nb-2],  b[connect_nb]);
	  if( b[connect_nb-2]->size > b[connect_nb-2]->initsize){
	    b[connect_nb-2]->borderorigin = realloc( b[connect_nb-2]->borderorigin, b[connect_nb-2]->size * sizeof(BorderIndex));
	    check_alloc(b[connect_nb-2]->borderorigin, 4610);
	    for(size_t i =b[connect_nb-2]->initsize; i<b[connect_nb-2]->size; i++)
	      b[connect_nb-2]->borderorigin[i] = (BorderIndex) {.b = b[connect_nb-2], .i = BOTTOM};
	  }
	  memoryy -= b[connect_nb-2]->initsize * sizeof(BorderIndex);
	  memoryy += b[connect_nb-2]->size * sizeof(BorderIndex);
	  maxmemoryy = MAX(memoryy, maxmemoryy);
  
	  b[connect_nb-1] = update(b[connect_nb-1],  b[connect_nb]);
	  free_boundary(b[connect_nb]);
	  send_updated_boundary(b[connect_nb-1], other_rank);
	  free_boundary(b[connect_nb-1]);
	  connect_nb -=  2;
	} else {
	  int other_rank = my_rank - base_w;
	  b[connect_nb] = receive_updated_boundary(b[connect_nb], other_rank);
	  merged--;
	}
      }
      else {
	merged--;
      }      
      base_w /= 2;
      grid_w = grid_w*2; 
    }
  }

  /*Applying the changes to the local component tree */
  // long node_mtree = *tree_size;
  maxtree = apply_changes(maxtree, b[0], tree_size);
  // long node_init_btree = b[0]->initsize;
  // long node_added = b[0]->size - b[0]->initsize;
  //  MPI_Allreduce(&node_init_btree, &node_init_btree, 1, MPI_LONG, MPI_MAX, MPI_COMM_WORLD);
  //  MPI_Allreduce(&node_added, &node_added, 1, MPI_LONG, MPI_MAX, MPI_COMM_WORLD);
  //*tree_size = *tree_size + (b[0]->size - b[0]->initsize);
  //printfrank("nodes in maxtree %ld  nodes in boundtree %ld, nodes added %ld, size added %ld (%ld)", 0, node_mtree,node_init_btree, node_added, node_added * (sizeof(MaxNode)+24), sizeof(MaxNode));
  free_boundary(b[0]);
  return maxtree;
} /* correct_borders */


MaxNode *apply_changes(MaxNode *tree, Boundary *b, size_t *tree_size) {
  /*Reallocating the local tree to include the new nodes */
  size_t newsize = (*tree_size + ( b->size - b->initsize));
  tree = realloc(tree, newsize * sizeof(MaxNode));  
  check_alloc(tree, 310); 
  memoryy -= *tree_size * sizeof(MaxNode);
  memoryy += newsize * sizeof(MaxNode);
  maxmemoryy = MAX(memoryy, maxmemoryy);
  

  
  /*Adding the new nodes */
  for(size_t i = b->initsize; i < b->size; i++){
    b->array[i].index = *tree_size;
    tree[(*tree_size)++] = b->array[i];
  }

  /* Updating the nodes's parents*/
  for(size_t i =0; i < b->size; i++){
    if(i<b->initsize)
      tree[b->array[i].index].area = b->array[i].area;    
    idx b_parent = b->borderparent[i].i;
    if(b_parent == BOTTOM)
      tree[b->array[i].index].parent = BOTTOM;    
    else
      tree[b->array[i].index].parent = b->array[b_parent].index;
  }

  return tree;
} /* apply_changes */


Boundary *update( Boundary *b, Boundary *c) {

  /* reallocate memorym */
  // b = realloc_b(b, 1.5*b->size);

  for(size_t i = 0; i< b->size; i++){
    BorderIndex x = (BorderIndex) {.b = b, .i = i};
    BorderIndex z = b_levelroot(x);
    BorderIndex s;
    x = b_levelrootinbound(x);
    if(bi_equal(z,x))
      s = (BorderIndex) {.b = c, .i = z.b->array[z.i].borderindex};
    else
      s = (BorderIndex) {.b = c, .i = x.b->othborderlr[x.i].b->array[x.b->othborderlr[x.i].i].borderindex};    
    if(!(x.b->reached[x.i]))
      b = update_branch(x,z,s);
  }    
  printfrank("Mem update: %ld",0, memoryy);
  if(b->allocsize>b->size)
    b = realloc_b(b, b->size);

  return b;

}

Boundary *update_branch(BorderIndex x, BorderIndex z, BorderIndex s){

  /* If the node is not in the combined boundary (s bottom) we update it using the info from the two local btres */
  while(is_bottom(s)  && !(x.b->reached[x.i])){  
    update_node_attribute(x,z);
    z = b_parent_lr(z);
    if(x.b == z.b)
      x.b->borderparent[x.i] = z;
    else{
      if(!is_bottom(z.b->othborderlr[z.i]))
	x.b->borderparent[x.i] = z.b->othborderlr[z.i];
      else{
	if(x.b->size == x.b->allocsize)
	  x.b = realloc_b(x.b, 1.5*x.b->size);
	
	adding_node(x,z);
	x.b->array[x.b->size-1].borderindex = z.b->array[z.i].borderindex;
	x.b->othborderlr[x.b->borderparent[x.i].i] = z;
	z.b->othborderlr[z.i] = x.b->borderparent[x.i];
      }
    }
    x.b->reached[x.i] = true;
    x = x.b->borderparent[x.i];
    s.i = z.b->array[z.i].borderindex;
  }
  BorderIndex origin;
  s = b_levelroot(s);

   while(!(x.b->reached[x.i]) && !is_bottom(s)){
    update_node_attribute(x,s);
    s = b_parent_lr(s);
    
    if(!is_bottom(s)){
      origin = s.b->borderorigin[s.i];
      if(x.b == origin.b){
	x.b->borderparent[x.i] = origin;
      } else if(!is_bottom(origin) && !is_bottom(origin.b->othborderlr[origin.i])){
	x.b->borderparent[x.i] = origin.b->othborderlr[origin.i];
      } else{
	if(x.b->size == x.b->allocsize)
	  x.b = realloc_b(x.b, 1.5*x.b->size);       
	adding_node(x,s);
	x.b->array[x.b->size-1].borderindex = s.i;	 
	if(!is_bottom(origin)){
	  x.b->othborderlr[x.b->size-1] = origin;
	  origin.b->othborderlr[origin.i] = (BorderIndex) {.b = x.b, .i = x.b->size-1};
	}
	else
	  s.b->borderorigin[s.i] =  (BorderIndex) {.b = x.b, .i = x.b->size-1};
      }
       x.b->reached[x.i] = true;
      x = x.b->borderparent[x.i];
    }
    else{
      x.b->borderparent[x.i].i = BOTTOM;
      x.b->reached[x.i] = true;
    }   	  
   }
   return x.b;
}

void update_node_attribute(BorderIndex x, BorderIndex s){
  x.b->array[x.i].area = s.b->array[s.i].area;
}

void adding_node(BorderIndex x, BorderIndex s){
  x.b->array[x.b->size] = s.b->array[s.i];
  x.b->reached[x.b->size] =  false;
  x.b->borderparent[x.i].b = x.b;
  x.b->borderparent[x.i].i = x.b->size;
  x.b->size++;
}
  
Boundary *combine(Boundary *a, Boundary *b, Direction d) {

  /* Given two boundary trees a and b, return the combined tree c */
  size_t upper_bound = a->size + b->size;
    
  Boundary *c = calloc(1, sizeof(Boundary));
  check_alloc(c, 441);
  memoryy += sizeof(Boundary);
  maxmemoryy = MAX(memoryy, maxmemoryy);
  
  c->offset = calloc(5, sizeof(long));
  check_alloc(c->offset, 442);
  memoryy += 5 * sizeof(long);
  maxmemoryy = MAX(memoryy, maxmemoryy);
  
  c->array = calloc(upper_bound, sizeof(MaxNode));
  check_alloc(c->array, 443);
  memoryy += upper_bound * sizeof(MaxNode);
  maxmemoryy = MAX(memoryy, maxmemoryy);
  
  
  c->borderparent = calloc(upper_bound, sizeof(BorderIndex));
  check_alloc(c->borderparent, 444);
  memoryy += upper_bound * sizeof(BorderIndex);
  maxmemoryy = MAX(memoryy, maxmemoryy);
  
  c->borderorigin = calloc(upper_bound,sizeof(BorderIndex));
  check_alloc(c->borderorigin, 446);
  memoryy += upper_bound * sizeof(BorderIndex);
  maxmemoryy = MAX(memoryy, maxmemoryy);
  
 
  reset_borderindex(a);
  reset_borderindex(b);

  size_t s = 0; /* points to entry point in c */
  size_t i; /* index local to either a or b */

  if (d == HORIZONTAL) {
    /* new top = a,b,c */
    /* a,b */
    c->offset[0] = s; /* superfluous, but for clarity */
    for (i = a->offset[0]; i <= a->offset[1]; i++, s++) { /* '<=' for b */
      b_add(c, a, s, i);
    }
    /* c */
    for (i = b->offset[0]; i < b->offset[1]; i++, s++) {
      b_add(c,  b, s, i);
    }

    /* new right side = d*/
    c->offset[1] = s;
    for (i = b->offset[1]; i < b->offset[2]; i++, s++) {
      b_add(c,  b, s, i);
    }

    /* new bottom-side = e,f,g */
    /* e */
    c->offset[2] = s;
    for (i = a->offset[2]; i < a->offset[3]; i++, s++) {
      b_add(c, a, s, i);
    }

    /* f */
    i = b->offset[4] - 1;
    b_add(c, b, s, i);
    s++;

    /* g */
    for (i = b->offset[2]; i < b->offset[3]; i++, s++) {
      b_add(c, b, s, i);
    }

    /* new left side = h */
    c->offset[3] = s;
    for (i = a->offset[3]; i < a->offset[4]; i++, s++) {
      b_add(c, a, s, i);
    }

  } else { /* d == VERTICAL */
    /* new top = a */
    /* a */
    c->offset[0] = s; /* superfluous, but for clarity */
    for (i = a->offset[0]; i < a->offset[1]; i++, s++) { /* '<=' for b */
      b_add(c, a, s, i);
    }

    /* new right side = b,c,d */
    /* b */
    c->offset[1] = s;
    for (i = a->offset[1]; i < a->offset[2]; i++, s++) { /* '<=' for b */
      b_add(c, a, s, i);
    }
    /* c */
    i = a->offset[3] - 1;
    b_add(c,  a, s, i);
    s++;

    /* d */
    for (i = b->offset[1]; i < b->offset[2]; i++, s++) {
      b_add(c, b, s, i);
    }

    /* new bottom-side = e */
    c->offset[2] = s;
    for (i = b->offset[2]; i < b->offset[3]; i++, s++) {
      b_add(c,  b, s, i);
    }

    /* new left side = f,g,h */
    c->offset[3] = s;
    /* f */
    for (i = a->offset[3]; i < a->offset[4]; i++, s++) {
      b_add(c,  a, s, i);
    }
    /* g */
    i = b->offset[0];
    b_add(c,  b, s, i);
    s++;
    /* h */
    for (i = b->offset[3]; i < b->offset[4]; i++, s++) {
      b_add(c,  b, s, i);
    }
  }

  // debug("New border contains %ld nodes (excluding ancestors)", s);

  /* add ancestors */

  c->offset[4] = s;
  size_t ex = s;
  
  for (i = 0; i < ex; i++) {
    size_t curr = i;
    BorderIndex origin = c->borderorigin[curr];
    BorderIndex parent = b_levelroot(b_parent(origin));
    
    while (true) {
      if (is_bottom(parent)) {
        /* case 1: parent is bottom */
        c->borderparent[curr] = (BorderIndex) {.b = c, .i = BOTTOM};
        break;
      }

      idx c_idx = b_node(parent).borderindex; /* index in c, if not BOTTOM */
      if (c_idx != BOTTOM) {
        /* case 2: parent is in c, avoid duplicates */
        c->borderparent[curr] = (BorderIndex) {.b = c, .i = c_idx};
        break;
      }
      /* case 3: add parent to c */
      (parent.b)->array[parent.i].borderindex = s; /* to avoid duplicates in case 2 */
      c->borderorigin[s] = parent;
      c->array[s] = b_node(parent);
      c->borderparent[curr] = (BorderIndex) {.b = c, .i = s};
      curr = s;
      s++;
      parent =b_levelroot(b_parent(parent));;
    }
  }

  c->size = s;
  
  // info("Added %ld ancestors for %ld nodes", s - ex, s);

  c->array = realloc(c->array, s * sizeof(MaxNode));
  check_alloc(c->array, 447);
  memoryy -= upper_bound * sizeof(MaxNode);
  memoryy += s * sizeof(MaxNode);
  maxmemoryy = MAX(memoryy, maxmemoryy);
    
  c->borderparent = realloc(c->borderparent, s * sizeof(BorderIndex));
  check_alloc(c->borderparent, 448);
  memoryy -= upper_bound * sizeof(BorderIndex);
  memoryy += s * sizeof(BorderIndex);
  maxmemoryy = MAX(memoryy, maxmemoryy);
 

  c->reached = calloc(s, sizeof(bool));
  check_alloc(c->reached, 449);
  
  memoryy += s * sizeof(bool);
  maxmemoryy = MAX(memoryy, maxmemoryy);
 
  
  c->borderorigin = realloc(c->borderorigin, s * sizeof(BorderIndex));
  check_alloc(c->borderorigin, 4491);
  memoryy -= upper_bound * sizeof(BorderIndex);
  memoryy += s * sizeof(BorderIndex);
  maxmemoryy = MAX(memoryy, maxmemoryy);
 

  c->othborderlr = malloc(s * sizeof(BorderIndex));
  check_alloc(c->othborderlr, 4492);
  memoryy += s * sizeof(BorderIndex);
  maxmemoryy = MAX(memoryy, maxmemoryy);
 
  for (i = 0; i < c->size; ++i) 
    c->othborderlr[i] = (BorderIndex)  {.b = c, .i = BOTTOM};
  c->allocsize = s;
  c->initsize = s;
  return c;
}

void b_add(Boundary *c, Boundary *b, size_t s, size_t i) {
  c->borderorigin[s]  = (BorderIndex) {.b = b, .i = i};
  b->array[i].borderindex = s;
  c->array[s] = b->array[i];
}

void merge(Boundary *a, Boundary *b, Direction d) {
  /* there are two cases:
     HORIZONTALLY: merge the right side of a with the left side of b 
     VERTICALLY: merge the bottom of a with top of b, */
  
  size_t length;
  if(d==HORIZONTAL){
    length = a->offset[2] - a->offset[1]; /* here we forget 1, fixed below */
  } else {
    length = a->offset[3] - a->offset[2]; /* here we forget 1, fixed below */
  }

  /* traverse border */
  for (size_t c = 0; c < length + 1; c++) { /* +1 because the forgotten one */
        
    ulong area = 0, area1 = 0;
    BorderIndex x = idx_i(a, c, d, length); /* get node i border from border a */
    BorderIndex y = idx_j(b, c, d); /* get node i index from border b */
    BorderIndex z,h;
    
    x = b_levelroot(x);
    y = b_levelroot(y);
    
    if (b_gval(x) < b_gval(y)) {
      h = x; x = y; y = h;
    }

    while (!bi_equal(x, y) && !is_bottom(y)) {
      z = b_parent_lr(x);
       
      if (!is_bottom(z) && (b_gval(z) >= b_gval(y))){
	x.b->array[x.i].area += area;
	x = z;
      } else{
	if(b_gval(x) == b_gval(y)){
	  if(!is_bottom(x.b->othborderlr[x.i]) && !is_bottom(y.b->othborderlr[y.i])){ // If a link between the levelroots in each boundary has already been done 
	    if(x.b == y.b){ 
	      x.b->borderparent[x.i] = y;
	      (x.b->othborderlr[x.i].b)->borderparent[x.b->othborderlr[x.i].i] = y.b->othborderlr[y.i];
	    }else {
	      x.b->borderparent[x.i] = y.b->othborderlr[y.i];
	      x.b->othborderlr[x.i].b->borderparent[x.b->othborderlr[x.i].i] = y;
	    }
	    x.b->othborderlr[x.i].b->othborderlr[x.b->othborderlr[x.i].i].i = BOTTOM;
	    x.b->othborderlr[x.i].i = BOTTOM;
	  } else if(!is_bottom(x.b->othborderlr[x.i])){
	    area = y.b->array[y.i].area;
	    h = b_parent_lr(y);
	    if(x.b == y.b)
	      y.b->borderparent[y.i] = x;
	    else
	      y.b->borderparent[y.i] = x.b->othborderlr[x.i];
	    y = h;
	    continue;
	  } else if( !is_bottom(y.b->othborderlr[y.i])){
	    if(x.b == y.b)
	      x.b->borderparent[x.i] = y;
	    else
	      x.b->borderparent[x.i] = y.b->othborderlr[y.i];
	  }else {
	    if(x.b != y.b){
	      y.b->othborderlr[y.i] = x;
	      x.b->othborderlr[x.i] = y;
	    }
	    x.b->borderparent[x.i] = y;
	  }
	}else
	  x.b->borderparent[x.i] = y;

	area1 = x.b->array[x.i].area + area;
	area = x.b->array[x.i].area;
	x.b->array[x.i].area = area1;

	x = y;
	y = z;
      }	
    }
	 
    if (is_bottom(y)) {
      while (!is_bottom(x)) {
	x.b->array[x.i].area += area;
	x = b_parent_lr(x);
      }
    }
  }
} //* merge */


BorderIndex idx_i(Boundary *a, size_t c, Direction d, size_t length) {
  BorderIndex bi;
  bi.b = a;

  if (d == VERTICAL) { /* return indices from bottom of a */
    if (c == 0) {
      bi.i = (idx) (a->offset[4] - 1); /* bottomleft */
    } else {
      bi.i = (idx) (a->offset[2] + c - 1); /* rightward from bottomleft+1 */
    }
  } else {
    /* HORIZONTAL */
    /* return indexes from right side of a */
    if (c != length) {
      bi.i = (idx) (a->offset[1] + c); /* downward from topright+1 */
    } else { /* last pixel */
      bi.i = (idx) (a->offset[3] - 1); /* bottomright */
    }
  }
  return bi;
}

BorderIndex idx_j(Boundary *b, size_t c, Direction d) {
  BorderIndex bi;
  bi.b = b;

  if (d == VERTICAL) { /* return index from top of b */
    bi.i = (idx) c;
  } else {
    /* HORIZONTAL */
    /* return indexes from left side of b */
    if (c == 0) {
      bi.i = (idx) 0; /* topleft */
    } else {
      bi.i = (idx) (b->offset[3] + c - 1); /* downward from topleft+1 */
    }
  }
  return bi;
}
   
Boundary *create_boundary(MaxNode *maxtree, size_t width, size_t height, size_t grid_w, size_t grid_h) {
  
  /*  We define a boundary of a as an array of MaxNodes, an array of offsets and the size of the array
      The array exists of {north,east,south,west,ancestors} concatenated in one piece of contiguous memoryy
      The ancestors are all nodes in the tree from every node in the border to the bottom following the parent relation
      The array of offset is {0,+width-1, +height-1, +width-1,}; which are the starting points of the respective borders in the array
  */
  size_t myrow = rank() /grid_w;
  size_t mycol = rank() % grid_w;
  
  Boundary *b = calloc(1, sizeof(Boundary));
  check_alloc(b, 4200);
  memoryy +=  sizeof(Boundary);
  maxmemoryy = MAX(memoryy, maxmemoryy);
 
  b->offset = calloc(5, sizeof(long));
  check_alloc(b->offset, 4201);
  memoryy += 5 * sizeof(long);
  maxmemoryy = MAX(memoryy, maxmemoryy);
 

  size_t b_limit = (2*width + 2*height -4);
  b->allocsize = b_limit;

  b->array = malloc( b_limit * sizeof(MaxNode));
  check_alloc(b->array, 4202);
  memoryy += b_limit * sizeof(MaxNode);
  maxmemoryy = MAX(memoryy, maxmemoryy);
 

  b->borderparent = malloc( b_limit * sizeof(BorderIndex));
  check_alloc(b->borderparent, 4202);
  memoryy += b_limit * sizeof(BorderIndex);
  maxmemoryy = MAX(memoryy, maxmemoryy);
 
  for (size_t i = 0; i < b_limit; ++i) 
    b->borderparent[i] = (BorderIndex) {.b = b, .i = BOTTOM};

  /*
    Order:
    0: topleft to topright-1
    1: topright to bottomright-1
    2: bottomleft+1 to bottomright
    3: topleft+1 to bottomleft
    4: ancestors
  */
  if(myrow > 0)
    add_side(b, maxtree, 0, width - 1, 0, 1);
  else if(mycol > 0)
    add_side(b, maxtree, 0, 1, 0, 1);
  else
    add_side(b, maxtree, 0, 0, 0, 1);


  if(mycol != grid_w -1)
    add_side(b, maxtree, 1, height - 1, width - 1, width);
  else if(myrow > 0)
    add_side(b, maxtree, 1, 1, width - 1, width);
  else
    add_side(b, maxtree, 1, 0, width - 1, width);

  if(myrow != grid_h -1)
    add_side(b, maxtree, 2, width - 1, (height - 1)*width + 1, 1);
  else if(mycol != grid_w-1)
    add_side(b, maxtree, 2, 1, (height - 1)*width + width - 1, 1);
  else
    add_side(b, maxtree, 2, 0, (height - 1)*width + 1, 1);

  if(mycol > 0)
    add_side(b, maxtree, 3, height - 1, width, width);
  else if(myrow != grid_h -1)
    add_side(b, maxtree, 3, 1, (height - 1)*width, width);
  else
    add_side(b, maxtree, 3, 0, width, width);

  b = add_ancestors(b, maxtree);
  check_boundary(b);

   /* shrink */
  b->array = realloc(b->array, b->size * sizeof(MaxNode));
  check_alloc(b->array, 4223);
  memoryy -= b_limit * sizeof(MaxNode);
  memoryy += b->size * sizeof(MaxNode);
  maxmemoryy = MAX(memoryy, maxmemoryy);
 

  
  b->borderparent = realloc(b->borderparent, b->size * sizeof(BorderIndex));
  check_alloc(b->borderparent, 4224);
  memoryy -= b_limit * sizeof(BorderIndex);
  memoryy += b->size * sizeof(BorderIndex);
  maxmemoryy = MAX(memoryy, maxmemoryy);
 

  // Allocating the remaining variables
  b->reached = calloc(b->size, sizeof(bool));
  check_alloc(b->reached, 4204);
  memoryy += b->size * sizeof(bool);
  maxmemoryy = MAX(memoryy, maxmemoryy);
   
  b->othborderlr = malloc(b->size * sizeof(BorderIndex));
  check_alloc(b->othborderlr, 4206);
  memoryy += b->size * sizeof(BorderIndex);
  maxmemoryy = MAX(memoryy, maxmemoryy);
 

  for (size_t i = 0; i < b->size; ++i) 
    b->othborderlr[i] = (BorderIndex) {.b = b, .i = BOTTOM};
  b->allocsize = b->size;
  b->initsize = b->size;
  // printfrank("Btree size %d, balloc size: %d",0, b->size, b->allocsize);

  return b;
} /* create_boundary */

void add_side(Boundary *b, MaxNode *maxtree, int side, size_t length, size_t offset, long increment) {
  // debug("Adding subborder from node %p length %zu offset %zu increment %ld", (void *)maxtree, length, offset, increment);
  b->offset[side] = b->size;
  b->size += length;
  
  for (size_t i = 0; i < length; ++i) {
    size_t j = offset + (i * increment);
    size_t indx = b->offset[side] + i;
    maxtree[j].borderindex = indx;
    b->array[indx] = maxtree[j];
  }
} /* add_side */

Boundary *add_ancestors(Boundary *b, MaxNode *maxtree) {
  /* We have to keep track of all nodes that are ancestors which are not in the border */
  /* Note that in order to visit every node only once, we need to iterate over the border, and only add ancestors that are not in the border nor already in the map. Administration is done here. */
  b->offset[4] = b->size;
  size_t origsize = b->size;

  /*
  The theoretical upper limit of the boundary including ancestors is:
  (G-1) * (L/2)
  where G is the number of gray levels, and L is the length of the border
  We therefore allocate the upper limit and shrink afterwards.
  To save memoryy, one could grow the array slowly, at the expense of speed.

  As the size of the ancestors map is dynamic and depending on the input image, we then have to make an assumption for a good initial size of the map. We here chose to make it as big as the border itself, and will double when needed, which is O(n).

    if(b->size >= *buffersize){
      *buffersize = 2*(*buffersize);
      debug("Doubling ancestor size to %zu", *buffersize);
      ancestors = realloc(ancestors, (*buffersize)*sizeof(MaxNode));
      check_alloc(ancestors, 117);
    }
  */
  
  // size_t ulimit = (size_t) ((MAXGREYVAL - 1) * (origsize / 2l) + 1);
  /* another limit is the number of nodes left */

  /* allocate memoryy */
  // b->array = realloc(b->array, ulimit * sizeof(MaxNode));
  //  check_alloc(b->array, 4221);
  // b->borderparent = malloc(ulimit * sizeof(BorderIndex));
  // check_alloc(b->borderparent, 4222);

  /* initialize */
  /* for (size_t i = 0; i < ulimit; ++i) {
    BorderIndex bi = {.b = b, .i = BOTTOM};
    b->borderparent[i] = bi;
    }*/
  
  /* b->size is the next free index of the boundary tree */
  for (size_t i = 0; i < origsize; ++i) {
    idx curr = b->array[i].index; /* index in maxtree */

    while (true) {
      idx parent = maxtree[curr].parent; /* index in maxtree */
      if (parent == BOTTOM) {
        /* borderparent stays BOTTOM */
        break; /* next! */
      } else if (maxtree[parent].borderindex != BOTTOM) {
        /* parent is already in the border */
        /* point there */
        idx bx = maxtree[curr].borderindex; /* index in border */
        b->borderparent[bx] = (BorderIndex) {.b=b, .i = maxtree[parent].borderindex}; /* equals b->size */

        /* b->borderparent[bx] is unchanged (will be set after receive ) */
        break; /* next! */
      } else {
        /* parent is not BOTTOM, parent is not in boundary tree yet */
	if(b->size == b->allocsize){
	  b->array = realloc(b->array, 1.5*b->size * sizeof(MaxNode));
	  check_alloc(b->array, 4601);
	  b->borderparent = realloc(b->borderparent,1.5*b->size * sizeof(BorderIndex));
	  check_alloc(b->borderparent, 4602);
	  b->allocsize = 1.5*b->size;
	  for (size_t j = b->size; j < b->allocsize; ++j) 
	    b->borderparent[j] = (BorderIndex) {.b = b, .i = BOTTOM};
	}
        /* add parent to border */
        maxtree[parent].borderindex = b->size;
        b->array[b->size] = maxtree[parent]; /* copy */

        /* set borderparent */
        idx bx = maxtree[curr].borderindex;
        b->borderparent[bx] = (BorderIndex) {.b=b, .i = maxtree[parent].borderindex}; /* equals b->size */

        b->size++; /* move to next empty index */
        curr = parent; /* next! */
      }
    }
  }
  return b;
} /* add_ancestors */

Boundary *realloc_b(Boundary *b, size_t newsize){
  b->array = realloc(b->array, newsize * sizeof(MaxNode));
  check_alloc(b->array, 4601);
  memoryy -= b->size * sizeof(MaxNode);
  memoryy += newsize * sizeof(MaxNode);
  maxmemoryy = MAX(memoryy, maxmemoryy);
 
  b->borderparent = realloc(b->borderparent,newsize * sizeof(BorderIndex));
  check_alloc(b->borderparent, 4602);
  memoryy -= b->size * sizeof(BorderIndex);
  memoryy += newsize * sizeof(BorderIndex);
  maxmemoryy = MAX(memoryy, maxmemoryy);
 

  b->reached = realloc(b->reached, newsize * sizeof(bool));
  check_alloc(b->reached, 4603);
  memoryy -= b->size * sizeof(bool);
  memoryy += newsize * sizeof(bool);
  maxmemoryy = MAX(memoryy, maxmemoryy);
 
  b->othborderlr = realloc( b->othborderlr, newsize * sizeof(BorderIndex));
  check_alloc(b->othborderlr, 4604);
  
  memoryy -= b->size * sizeof(BorderIndex);
  memoryy += newsize * sizeof(BorderIndex);
  maxmemoryy = MAX(memoryy, maxmemoryy);
 

  if(newsize > b->size){
    for(size_t i = b->size; i<newsize; i++)
      b->othborderlr[i] = (BorderIndex) {.b = b, .i = BOTTOM};
  }
  //b->borderorigin = realloc( b->borderorigin, newsize * sizeof(BorderIndex));
  //check_alloc(b->borderorigin, 4604);
  b->allocsize = newsize;
  return b;
}


void reset_borderindex(Boundary *b) {
  for (size_t i = 0; i < b->size; i++) 
    b->array[i].borderindex = BOTTOM;
  
}

BorderIndex b_levelroot(BorderIndex bi) {
  // trace("Getting levelroot of %ld in %p", bi.i, (void *)bi.b);

  BorderIndex ri = bi;
  if (is_bottom(ri)) {
    return bi; /* still BOTTOM */
  }
  else{
  value gv = b_gval(bi);

  while ( !is_bottom(b_parent(ri)) && (gv == b_gval(b_parent(ri))) ) {
    ri = b_parent(ri);
    if(bi_equal(ri, b_parent(ri))){
	error("[CONNEXION ERROR] Next levelroot %ld in %d");
	MPI_Abort(MPI_COMM_WORLD, 118);
      }
  }

  return ri;
  }
} /* b_levelroot */


 BorderIndex b_levelrootinbound(BorderIndex bi) {
   //trace("Searching boundary levelroot of %ld in %p", bi.i, (void *)bi.b);

  if (is_bottom(bi)) {
    return bi; /* still BOTTOM */
  }
  
  BorderIndex ri=bi, pi=bi;  
  value gv = b_gval(bi);

  while ( !is_bottom(b_parent(pi)) && gv == b_gval(b_parent(pi)) && b_parent(pi).b == bi.b ) {
    pi = b_parent(pi);
    ri = pi;
  }

  return ri;
} /* b_levelroot_bound */


bool is_bottom(BorderIndex bi) {
  return bi.i == BOTTOM;
} /* is_bottom */
 
bool bi_equal(BorderIndex ai, BorderIndex bi) {
  return (ai.b == bi.b) && (ai.i == bi.i);
} /* bi_equal */

BorderIndex b_parent_lr(BorderIndex bi) {
  /* get parent index */
  BorderIndex pi = b_parent(bi);
  /* get levelroot of parent */
  BorderIndex ri = b_levelroot(pi);
  // trace("b_parent_lr pi: %p %ld ri: %p %ld", pi.i, (void *)pi.b, ri.i, (void *)ri.b);
  return ri;
} /* b_parent_lr */

BorderIndex b_parent(BorderIndex bi) {
  BorderIndex pi = bi.b->borderparent[bi.i];
  return pi;
} /* b_parent */

MaxNode b_node(BorderIndex bi) {
  if (bi.b == NULL) {
    error("asking node of non-initilized border %ld", bi.i);
  }
  return (bi.b)->array[bi.i];
}

value b_gval(BorderIndex bi) {
  if (bi.b == NULL) {
    error("asking gval of non-initilized border");
    MPI_Abort(MPI_COMM_WORLD, 118);
  }
  if ((size_t)bi.i > ((bi.b)->size)) {
    warn("asking gval of index (%ld) higher than size (%ld) of border", bi.i, (bi.b)->size);
    /* MPI_Abort(MPI_COMM_WORLD, 117); */
    return BOTTOM;
  }

  return ((bi.b)->array[bi.i]).gval;
} /* b_gval */


void free_boundary(Boundary *b) {
  free(b->array);
  memoryy -= b->size * sizeof(MaxNode);
  free(b->offset);
  memoryy -= 5 * sizeof(long);
  free(b->borderparent);
  memoryy -= b->size * sizeof(BorderIndex);
  free(b->othborderlr);
  memoryy -= b->size * sizeof(BorderIndex);
  free(b->borderorigin);
  memoryy -= b->size * sizeof(BorderIndex);
  memoryy -= b->size * sizeof(bool);
  free(b->reached);
  free(b);
  memoryy -= sizeof(Boundary);
} /* free_boundary */

