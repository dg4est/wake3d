//
//  alloc.c
//  cdriver
//
//  Created by Michael Brazell on 1/5/17.
//  Copyright © 2017 Michael J. Brazell. All rights reserved.
//

/* some of this taken from p4est sc */


#include "alloc.h"
#include <stdlib.h>

static int alloc_count = 0;

void * my_malloc (size_t size){

  void *ret;

//  if(size == 0) printf("oh no allocating a length zero array\n");

  ret = malloc (size);

  if (size > 0) {
    ++alloc_count;
  }
  else {
    alloc_count += ((ret == NULL) ? 0 : 1);
  }

  return ret;
}

void * my_calloc (size_t n, size_t size){

  void *ret;

  //  if(size == 0) printf("oh no allocating a length zero array\n");

  ret = calloc (n,size);

  if (n*size > 0) {
    ++alloc_count;
  }
  else {
    alloc_count += ((ret == NULL) ? 0 : 1);
  }

  return ret;
}

void my_free (void *ptr) {
  if (ptr != NULL) {
    --alloc_count;
  }
  free (ptr);
}

int get_alloc_count(){
  return alloc_count;
}
