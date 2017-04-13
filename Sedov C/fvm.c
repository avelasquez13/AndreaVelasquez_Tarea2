/*
 * fvm.c
 *
 *  Created on: Apr 11, 2017
 *      Author: felipe
 */
#include <stdio.h>
#include "struct.h"
#include "fvm.h"
#include "space.h"

/**
 * Calcula la evolución de la onda de choque hasta un r_final que entra por parametro
 *
 */
void evolve(physics_grid *P, U_grid *U, F_grid *F, double r_final){
	//TODO
}

/**
 *Calcula una iteración del problema (avanza un dt)
 */
void step(physics_grid *P, U_grid *U, F_grid *F){
	//TODO
}

/**
 * Actualiza los parametros del physics_grid con base en los valores del U_grid
 */
void actualizarP(physics_grid *P, U_grid *U){
	//TODO
        int val=0;
	int x;
	int y;
	int z;
	int pos;
	for (pos=0;pos<P.N_cells;pos++){ // Densidad
	  P.P[val*P.N_cells+pos] = U.U[pos];
	}

	val++;
	for (z=0;z<P.N_z;z++){ // Presion
	  for (y=0;y<P.N_y;y++){
	    for (x=0:x<P.N_x;x++){
	      pos = pos(x,y,z,P.N_x,P.N_y);
	      P.P[val*P.N_cells+pos] = (GAMMA - 1)*(U.U[4*P.N_cells+pos] - 0.5*(U.U[1*P.N_cells+pos]*U.U[1*P.N_cells+pos] + U.U[2*P.N_cells+pos]*U.U[2*P.N_cells+pos] + U.U[3*P.N_cells+pos]*U.U[3*P.N_cells+pos])/U.U[pos]);
	    }
	  }
	}

	val++;
	for (z=0;z<P.N_z;z++){ // Vel x
	  for (y=0;y<P.N_y;y++){
	    for (x=0:x<P.N_x;x++){
	      pos = pos(x,y,z,P.N_x,P.N_y);
	      P.P[val*P.N_cells+pos] = U.U[(val-1)*P.N_cells+pos]/U.U[pos];
	    }
	  }
	}

	val++;
	for (z=0;z<P.N_z;z++){ // Vel y
	  for (y=0;y<P.N_y;y++){
	    for (x=0:x<P.N_x;x++){
	      pos = pos(x,y,z,P.N_x,P.N_y);
	      P.P[val*P.N_cells+pos] = U.U[(val-1)*P.N_cells+pos]/U.U[pos];
	    }
	  }
	}

	val++;
	for (z=0;z<P.N_z;z++){ // Vel z
	  for (y=0;y<P.N_y;y++){
	    for (x=0:x<P.N_x;x++){
	      pos = pos(x,y,z,P.N_x,P.N_y);
	      P.P[val*P.N_cells+pos] = U.U[(val-1)*P.N_cells+pos]/U.U[pos];
	    }
	  }
	}
}

/**
 * Actualiza los valores del F_grid con base en los valores del U_grid
 */
void actualizarF(U_grid *U, F_grid *Fp, F_grid *Fm){
	//TODO
}

/**
 * Encuentra la posición de la onda de choque
 */
double radioChoque(physics_grid *P){
	//TODO
	double r;

	return r;
}

/**
 * Actualiza los valores de entalpía
 */
void h(physics_grid *P, double* h){
	//TODO

}

/**
 * Actualiza los valores de cs
 */
void cs(physics_grid *P, double* cs){
	//TODO

}

/**
 * Encuentra el valor de (u, v, ó w)+cs máximo
 */
double vmax(physics_grid *P, double* cs){
	double vmax;
	//TODO
	return vmax;
}

/**
 * Encuentra el dt apropiado para las condiciones actuales
 */
double dt(physics_grid *P, double* cs){
	double dt;
	//TODO
	return dt;
}

/**
 * Devuelve un vector con la propiedad deseada del physics_grid
 * @param int propiedad toma los valores de las constnates:
 * RHO= 0
 * PRESSURE= 1
 * VX= 2
 * VY= 3
 * VZ= 4
 */
double* propiedad(physics_grid *P, int propiedad){
	//TODO
}
