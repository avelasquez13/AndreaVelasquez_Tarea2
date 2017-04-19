/*
 * space.h
 *
 *  Created on: Apr 11, 2017
 *      Author: felip
 */

#ifndef SPACE_H_
#define SPACE_H_

#define ejX 0
#define ejY 1
#define ejZ 2


int posi(int x, int y, int z, int Nx, int Ny);
void coord(int *coord, int pos, int Nx, int Ny);
int posF(int x, int y, int z, int eje, int valor, int Nx, int Ny, int Nz);
int radioSq(physics_grid *P, int x);

#endif /* SPACE_H_ */
