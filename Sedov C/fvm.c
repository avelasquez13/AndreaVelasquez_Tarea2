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
 * Calcula la evoluci�n de la onda de choque hasta un r_final que entra por parametro
 *
 */
void evolve(physics_grid *P, U_grid *U, F_grid *F, double r_final){
	//TODO
}

/**
 *Calcula una iteraci�n del problema (avanza un dt)
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
void actualizarF(int eje, U_grid *U, F_grid *F){
	//TODO
}

/**
 * Encuentra la posici�n de la onda de choque
 */
double radioChoque(U_grid *U){
	//TODO
	double r;

	return r;
}

/**
 * Actualiza los valores de entalp�a
 */
void h(U_grid *U, double* h){
	//TODO

}

/**
 * Actualiza los valores de cs
 */
void cs(U_grid *U, double* cs){
	//TODO

}

/**
 * Encuentra el valor de (u, v, � w)+cs m�ximo
 */
double vmax(U_grid *U, double* cs){
	double vmax;

	return vmax;
}

/**
 * Encuentra el dt apropiado para las condiciones actuales
 */
double dt(U_grid *U, double* cs){
	double dt;

	return dt;
}
