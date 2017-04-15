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

double* create_listNdoubles(physics_grid *P){
	double *G;
	  if(!(G = malloc(P->N_cells*sizeof(FLOAT)))){
	    fprintf(stderr, "Problem with F allocation");
	    exit(1);
	  }
	  init_to_zero(G, P->N_cells);
	  return G;
}
int* create_listNints(physics_grid *P){
	int *G;
	  if(!(G = malloc(P->N_cells*sizeof(int)))){
	    fprintf(stderr, "Problem with F allocation");
	    exit(1);
	  }

	  int i;
	  for(i=0;i<P->N_cells;i++){
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
 * Inicializa la lista de radios de las celdas y la lista de posiciones ordenada por radio ascendente
 */
int init_radios(physics_grid *P, double *radios, double *dist, double *rho, int *posiciones){
  int i, x, y, z, pos,length;
  double rad_cuadrados;

  for (z=0;z<P->N_z;z++){ // Guarda radio para cada posicion
    for (y=0;y<P->N_y;y++){
      for (x=0;x<P->N_x;x++){
    	  pos = posi(x,y,z, P->N_x,P->N_y);
    	  posiciones[pos] = pos;
    	  rad_cuadrados = pow(P->delta_x*(x - P->N_x/2),2) + pow(P->delta_y*(y - P->N_y/2),2) + pow(P->delta_z*(z - P->N_z/2),2);
    	 radios[pos] = sqrt(rad_cuadrados);
      }
    }
  }
  ordenarPorRadios(radios, posiciones, P->N_cells);

  pos = 1;
  dist [0] = radios[0];
  length = 1;
  for (i=1;i<P->N_cells/8;i++){ // Recorre dist para llenarlo
    while (dist[i-1] == radios[pos-1] && pos<P->N_cells){ // Busca cambio en radios
      pos++;
    }
    if (pos<P->N_cells){ // Si sigue en una celda valida asigna dist
      dist[i] = radios[pos-1];
      length = i+1;
    }
    else{ // Ya recorrio radios, los valores que quedan de la lista dist sobran
      dist[i] = -1;
    }
  }
  return length;
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

/**
 * Ordena los radios y las posiciones por radio ascendente. Utiliza el algoritmo bubble sort.
 */
void ordenarPorRadios(double *radios, int *posiciones, int length){
	int i,j;
	double rad,pos;
	for (i = 1; i < length; ++i) {
		for (j = 0; j < length-i; ++j) {
		    printf("Pos ord: %d,%d\n",i,j);
			rad=radios[j];
			pos=posiciones[j];
			if(rad>radios[j+1]){
						radios[j]=radios[j+1];
						posiciones[j]=posiciones[j+1];
						radios[j+1]=rad;
						posiciones[j+1]=pos;
				//TODO arreglar esto
					}
		}
	}
}
