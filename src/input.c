//
//  input.c
//  cdriver
//
//  Created by Michael Brazell on 1/6/17.
//  Copyright © 2017 Michael J. Brazell. All rights reserved.
//

#include <string.h>
#include <stdlib.h>
#include <stdio.h>

#include "input.h"
#include "defs.h"

int find_keyword_three_doubles(char* filename, char *keyword, double *dbl1, double *dbl2, double *dbl3, char mandatory){
  FILE *fp;
  char buff[buff_size];
  fp = fopen(filename,"r");

  unsigned long length = strlen(keyword);

  while(fp!=NULL && fgets(buff,sizeof(buff),fp) !=NULL){
    if(strstr(buff,keyword)){
      sscanf(&buff[length],"%le %le %le",dbl1,dbl2,dbl3);
      sscanf(&buff[length],"%le, %le, %le",dbl1,dbl2,dbl3);
      sscanf(&buff[length],"%le %le, %le",dbl1,dbl2,dbl3);
      sscanf(&buff[length],"%le, %le %le",dbl1,dbl2,dbl3);

      fclose(fp);
      return 0;
    }
  }

  if(mandatory) printf("\033[1;31m[wake3d] Mandatory Keyword [%s] not found!\033[0m\n",keyword);
  fclose(fp);
  return 1;
}

int find_keyword_integer(char* filename, char *keyword, int *integer, char mandatory){
  FILE *fp;
  char buff[buff_size];
  fp = fopen(filename,"r");

  unsigned long length = strlen(keyword);

  while(fp!=NULL && fgets(buff,sizeof(buff),fp) !=NULL){
    if(strstr(buff,keyword)){
      *integer = atoi(&buff[length]);
      fclose(fp);
      return 0;
    }
  }

  if(mandatory) printf("\033[1;31m[wake3d] Mandatory Keyword [%s] not found!\033[0m\n",keyword);
  fclose(fp);
  return 1;
}

int find_keyword_two_integers(char* filename, char *keyword, int *integer1, int *integer2, char mandatory){
  FILE *fp;
  char buff[buff_size];
  int ret;
  fp = fopen(filename,"r");

  unsigned long length = strlen(keyword);

  while(fp!=NULL && fgets(buff,sizeof(buff),fp) !=NULL){
    if(strstr(buff,keyword)){
      ret = atoi(&buff[length]);
      sscanf(&buff[length],"%d %d",integer1,integer2);
      fclose(fp);
      return 0;
    }
  }

  if(mandatory) printf("\033[1;31m[wake3d] Mandatory Keyword [%s] not found!\033[0m\n",keyword);
  fclose(fp);
  return 1;
}


int find_keyword_string(char *filename, char *keyword, char *string, char mandatory){
  FILE *fp;
  char buff[buff_size];
  unsigned long i,length;

  fp = fopen(filename,"r");
  length = strlen(keyword);

  while(fp!=NULL && fgets(buff,sizeof(buff),fp) !=NULL){
    if(strstr(buff,keyword)){

      int ind=0;
      for(i=length;i<buff_size;i++){
        if(buff[i] != ' '){
          string[ind] = buff[i];
          ind = ind+1;
        }
      }

      string[strcspn(string, "\n")] = 0;

      if(strlen(string) >= buff_size){
        printf("[wake3d] Warning: string length greater than buff_size set in driver.h : %lu out of %d\n",strlen(string),buff_size);
      }

      fclose(fp);
      return 0;
    }
  }

  if(mandatory) printf("\033[1;31m[wake3d] Mandatory Keyword [%s] not found!\033[0m\n",keyword);
  fclose(fp);
  return 1;
}


