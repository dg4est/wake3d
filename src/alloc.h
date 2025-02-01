//
//  alloc.h
//  cdriver
//
//  Created by Michael Brazell on 1/5/17.
//  Copyright © 2017 Michael J. Brazell. All rights reserved.
//

#ifndef alloc_h
#define alloc_h

#include <stdio.h>

#define my_alloc(t,n)      (t *) my_malloc ((n) * sizeof(t))
#define my_alloc_zero(t,n) (t *) my_calloc ((size_t) (n), sizeof(t))

void *my_malloc(size_t size);
void *my_calloc(size_t n, size_t size);
void my_free(void *ptr);
int get_alloc_count();

#endif /* alloc_h */
