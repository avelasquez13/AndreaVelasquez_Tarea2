/*
 * fvm.h
 *
 *  Created on: Apr 11, 2017
 *      Author: felip
 */

#ifndef FVM_H_
#define FVM_H_

void evolve(physics_grid *P, U_grid *U, F_grid *F, double r_final);
void step(physics_grid *P, U_grid *U, F_grid *F);
void actualizarP(physics_grid *P, U_grid *U);
void actualizarF(int eje, U_grid *U, F_grid *F);
double radioChoque(physics_grid *P);
void h(physics_grid *P, double* h);
void cs(physics_grid *P, double* cs);
double vmax(physics_grid *P, double* cs);
double dt(physics_grid *P, double* cs);

#endif /* FVM_H_ */
