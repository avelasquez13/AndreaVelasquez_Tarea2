#include "struct.h"
#include <stdio.h>

#define fname "sedov.dat"
#define ft "tiempo.dat"

void print_L(physics_grid *G){
  fprintf(stdout, "L_x = %f\n", G->L_x);  
  fprintf(stdout, "L_y = %f\n", G->L_y);
  fprintf(stdout, "L_z = %f\n", G->L_z); 
}

void print_list(double *list, int length){
  int i;
  FILE *fp;
  fp = fopen(fname, "a");
  
  for (i=0;i<length;i++){
    fprintf(fp,"%f ",list[i]);
  }
  fprintf(fp,"\n");
  fclose(fp);
}

void print_time(double tiempo){
  FILE *fp;
  fp = fopen(ft, "a");
  
  fprintf(fp,"%f\n",tiempo);
  
  fclose(fp);
}
