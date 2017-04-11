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
double radioChoque(U_grid *U);
void h(U_grid *U, double* h);
void cs(U_grid *U, double* cs);
double vmax(U_grid *U, double* cs);
double dt(U_grid *U, double* cs);

#endif /* FVM_H_ */
