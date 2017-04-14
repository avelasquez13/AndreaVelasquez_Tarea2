/*
 * space.c
 *
 *  Created on: Apr 11, 2017
 *      Author: felip
 */
#include <stdio.h>
#include "struct.h"
#include "space.h"

/**
 * Retorna la posición en el arreglo del que corresponde a las coordenadas (x,y,z)
 * Primero se recorre x, luego y, luego z
 */
int pos(int x, int y, int z, int Nx, int Ny){
	int pos;
	pos = x + Nx*y + Nx*Ny*z;
	return pos;
}

/**
 * Retorna las coordenadas que corresponden a la posicion  del arreglo pos.
 */
void coord(int *coord, int pos, int Nx, int Ny){
	coord[2] = pos/(Nx*Ny);
	coord[1] = (pos - Nx*Ny*coord[2])/Nx;
	coord[0] = (pos - Nx*Ny*coord[2] - Nx*coord[1]);
}

/**
 * Retorna la posición en el arreglo de F que corresponde a las coordenadas (x,y,z) y el eje que entra por parametro
 * @param eje=0 para Fx, eje=1 para Fy, eje=2 para Fz
 */
int posF(int x, int y, int z, int eje, int valor, int Nx, int Ny, int Nz){
	int pos;
	pos = eje*(NDIM+2)*Nx*Ny*Nz + valor*Nx*Ny*Nz + pos(x,y,z,Nx,Ny);
	return pos;
}
