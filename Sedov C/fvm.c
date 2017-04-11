/*
 * fvm.c
 *
 *  Created on: Apr 11, 2017
 *      Author: felipe
 */
#include <stdio.h>
#include "struct.h"
#include "fvm.h"

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
