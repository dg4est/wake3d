//
//  input.h
//  cdriver
//
//  Created by Michael Brazell on 1/6/17.
//  Copyright © 2017 Michael J. Brazell. All rights reserved.
//

#ifndef input_h
#define input_h

int find_keyword_integer(char* filename, char *keyword, int *integer, char mandatory);
int find_keyword_two_integers(char* filename, char *keyword, int *integer1, int *integer2, char mandatory);
int find_keyword_string(char *filename, char *keyword, char *string, char mandatory);
int find_keyword_three_doubles(char* filename, char *keyword, double *dbl1, double *dbl2, double *dbl3, char mandatory);

#endif /* input_h */
