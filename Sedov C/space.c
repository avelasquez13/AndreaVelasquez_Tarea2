/*
 * space.c
 *
 *  Created on: Apr 11, 2017
 *      Author: felipe
 */
#include <stdio.h>
#include <stdlib.h>
#include "struct.h"
#include "space.h"
#include <math.h>
/**
 * Retorna la posición en el arreglo del que corresponde a las coordenadas (x,y,z)
 * Primero se recorre x, luego y, luego z
 */
int posi(int x, int y, int z, int Nx, int Ny){
	int pos;
	pos = x + Nx*y + Nx*Ny*z;
	return pos;
}

/**
 * Retorna las coordenadas que corresponden a la posicion  del arreglo pos.
 */
void coord(int *coord, int pos, int Nx, int Ny){
	coord[2] = pos/(Nx*Ny);
	int rem=pos%(Nx*Ny);
	coord[1] = rem/Nx;
	rem= rem%Nx;
	coord[0] = rem;

	//Tambien se pueden calcular las coordenadas asi
	/*coord[2] = pos/(Nx*Ny);
	coord[1] = (pos - Nx*Ny*coord[2])/Nx;
	coord[0] = (pos - Nx*Ny*coord[2] - Nx*coord[1]);*/


	printf("coord ind %d: %d,%d,%d\n",pos,coord[0],coord[1],coord[2]);
}

/**
 * Retorna la posición en el arreglo de F que corresponde a las coordenadas (x,y,z) y el eje que entra por parametro
 * @param eje=0 para Fx, eje=1 para Fy, eje=2 para Fz
 */
int posF(int x, int y, int z, int eje, int valor, int Nx, int Ny, int Nz){
	int p;
	p = eje*(NDIM+2)*Nx*Ny*Nz + valor*Nx*Ny*Nz + posi(x,y,z,Nx,Ny);
	return p;
}

/**
 * Encuentra en que posicion de la lista de radios se encuentra la celda x
 */
int radioSq(physics_grid *P, int x){
	int c[3];
	coord(c,x,P->N_x,P->N_x);
	int radsq=pow((c[0]-P->N_x/2),2)+pow((c[1]-P->N_y/2),2)+pow((c[2]-P->N_z/2),2);
	return radsq;
}
