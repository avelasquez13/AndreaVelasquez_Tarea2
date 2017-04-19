#include "struct.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "init.h"
#include "space.h"

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

double* create_listNdoubles(int npoints){
	double *G;
	  if(!(G = malloc(npoints*sizeof(FLOAT)))){
	    fprintf(stderr, "Problem with F allocation");
	    exit(1);
	  }
	  init_to_zero(G, npoints);
	  return G;
}
int* create_listNints(int npoints){
	int *G;
	  if(!(G = malloc(npoints*sizeof(int)))){
	    fprintf(stderr, "Problem with F allocation");
	    exit(1);
	  }

	  int i;
	  for(i=0;i<npoints;i++){
	    G[i] = 0;
	  }
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

/**
 * Inicializa la lista de radios posibles para las celdas. estos son los radios que sean raiz de todos los enteros entre 0 y Nx^2+Ny^2+Nz^2
 */
void init_radios(physics_grid *P, double *radios, double *rho_avg, int *contador, int length){
	int i;
	for (i = 0; i < length; ++i) {
		radios[i]=sqrt(i);
	}
}

/**
 * Calcula la longitud del arreglo de radios como el cuadrado maximo de la distancia (ej 64^2+64^2+64^2)
 */
int length_radios(physics_grid *P){
	return pow((P->N_x/2),2)+pow((P->N_y/2),2)+pow((P->N_z/2),2)+1;
}
/**
 * Inicializa las condiciones iniciales de la explosion
 */
void init_conditions(U_grid *U, physics_grid *P){
  // TODO falta definir las condiciones iniciales
	int i;
	int ncells=U->N_cells;
	for (i = 0; i < ncells; ++i) {
		U->U[0*ncells+i]=1.177;
		U->U[4*ncells+i]=101325/0.4;
	}
	int pos=posi(P->N_x/2, P->N_y/2, P->N_z/2, P->N_x, P->N_y);
	U->U[4*ncells+pos]=1.177*pow(10,10);
}

