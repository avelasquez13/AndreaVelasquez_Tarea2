/*
 * fvm.h
 *
 *  Created on: Apr 11, 2017
 *      Author: felip
 */

#ifndef FVM_H_
#define FVM_H_

void evolve(physics_grid *P, U_grid *U, F_grid *F, double r_final);
void step(physics_grid *P, U_grid *U, F_grid *Fp);
void actualizarP(physics_grid *P, U_grid *U);
void actualizarF(U_grid *U, F_grid *Fp, F_grid *Fm);
double radioChoque(physics_grid *P);
void h(U_grid *U, double* h);
void cs(U_grid *U, double* cs);
double vmax(physics_grid *P, U_grid *U);
double dt(physics_grid *P, U_grid *U);
void propiedad(physics_grid *P, int propiedad, double *prop);
double presion(double *u_cell);

#endif /* FVM_H_ */