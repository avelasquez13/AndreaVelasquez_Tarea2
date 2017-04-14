#include "struct.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

void init_to_zero(FLOAT *p, int n_points){
  int i;
  for(i=0;i<n_points;i++){
    p[i] = 0.0;
  }
}


physics_grid * create_physics_grid(void){
  physics_grid *G;
  if(!(G = malloc(sizeof(physics_grid)))){
    fprintf(stderr, "Problem with data allocation\n");fflush(stdout);
    exit(0);
  } 
  
  G->L_x=0.0;
  G->L_y=0.0;
  G->L_z=0.0;
  G->delta_x=0.0;
  G->delta_y=0.0;
  G->delta_z=0.0;
  G->N_x=0;
  G->N_y=0;
  G->N_z=0;
  G->N_cells=0;
  G->P=NULL;
  return G;
}


U_grid * create_U_grid(void){
  U_grid *G;
  if(!(G = malloc(sizeof(U_grid)))){
    fprintf(stderr, "Problem with data allocation\n");fflush(stdout);
    exit(0);
  } 
  G->N_x=0.0;
  G->N_y=0.0;
  G->N_z=0.0;
  G->N_cells=0.0;
  G->U=NULL;
  return G;
}

F_grid * create_F_grid(void){
  F_grid *G;
  if(!(G = malloc(sizeof(F_grid)))){
    fprintf(stderr, "Problem with data allocation\n");fflush(stdout);
    exit(0);
  } 
  G->N_x=0.0;
  G->N_y=0.0;
  G->N_z=0.0;
  G->N_cells=0.0;
  G->F=NULL;
  return G;
}


void init_problem(physics_grid *P, U_grid *U, F_grid *F_p, F_grid *F_m){
  
  P->L_x = 256.0;
  P->L_y = 256.0;
  P->L_z = 256.0;    
  P->delta_x = 2.0;
  P->delta_y = 2.0;
  P->delta_z = 2.0;
  P->N_x = (int)(P->L_x/P->delta_x);
  P->N_y = (int)(P->L_y/P->delta_y);
  P->N_z = (int)(P->L_z/P->delta_z);
  P->N_cells = P->N_x * P->N_y * P->N_z;
  
  U->N_x = P->N_x;
  U->N_y = P->N_y;
  U->N_z = P->N_z;
  U->N_cells = P->N_cells;
  
  F_p->N_x = P->N_x;
  F_p->N_y = P->N_y;
  F_p->N_z = P->N_z;
  F_p->N_cells = P->N_cells;

  F_m->N_x = P->N_x;
  F_m->N_y = P->N_y;
  F_m->N_z = P->N_z;
  F_m->N_cells = P->N_cells;

  if(!(P->P = malloc(P->N_cells * (NDIM +2) * sizeof(FLOAT)))){
    fprintf(stderr, "Problem with pressure allocation");
    exit(1);
  }
  init_to_zero(P->P, P->N_cells * (NDIM +2));
  
  if(!(U->U = malloc(U->N_cells * (NDIM +2) * sizeof(FLOAT)))){
    fprintf(stderr, "Problem with U_1 allocation");
    exit(1);
  }
  init_to_zero(U->U, U->N_cells * (NDIM +2));
  
  if(!(F_p->F = malloc(F_p->N_cells * (NDIM) * (NDIM + 2) * sizeof(FLOAT)))){
    fprintf(stderr, "Problem with F allocation");
    exit(1);
  }
  init_to_zero(F_p->F, F_p->N_cells * NDIM * (NDIM +2));

  if(!(F_m->F = malloc(F_m->N_cells * (NDIM) * (NDIM + 2) * sizeof(FLOAT)))){
    fprintf(stderr, "Problem with F allocation");
    exit(1);
  }
  init_to_zero(F_m->F, F_m->N_cells * NDIM * (NDIM +2));
}

void init_radios(physics_grid *P, double *radios, int *posiciones){
  int x, y, z, pos;
  double rad_cuadrado;
  if(!(radios = malloc(P->N_cells*sizeof(FLOAT)))){
    fprintf(stderr, "Problem with F allocation");
    exit(1);
  }

  if(!(posiciones = malloc(P->N_cells*sizeof(int)))){
    fprintf(stderr, "Problem with F allocation");
    exit(1);
  }

  for (z=0;z<P.N_z;z++){ // Guarda radio para cada posicion
    for (y=0;y<P.N_y;y++){
      for (x=0;x<P.N_x;x++){
	pos = x + P.N_x*y + P.N_x*P.N_y*z;
	posiciones[pos] = pos;
	rad_cuadrados = ((x+0.5)*P.delta_x - 0.5*P.L_x)*((x+0.5)*P.delta_x - 0.5*P.L_x) + ((y+0.5)*P.delta_y - 0.5*P.L_y)*((y+0.5)*P.delta_y - 0.5*P.L_y) + ((z+0.5)*P.delta_z - 0.5*P.L_z)*((z+0.5)*P.delta_z - 0.5*P.L_z);
	radios[pos] = sqrt(rad_cuadrados);
      }
    }
  }

  // TODO: algoritmo para ordenar radios y con eso ordenar ordenar correspondientemente posiciones
}
