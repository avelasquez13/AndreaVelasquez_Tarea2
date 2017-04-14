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
        int val=0;
	int x;
	int y;
	int z;
	int pos;
	int i;
	double p;
	double u[5];
	for (pos=0;pos<P.N_cells;pos++){ // Densidad
	  P.P[val*P.N_cells+pos] = U.U[pos];
	}

	val++;
	for (z=0;z<P.N_z;z++){ // Presion
	  for (y=0;y<P.N_y;y++){
	    for (x=0:x<P.N_x;x++){
	      pos = pos(x,y,z,P.N_x,P.N_y);
	      for (i=0;i<5;i++){
		u[i] = U.U[i*P.N_cells+pos];
	      }
	      p = presion(u);
	      P.P[val*P.N_cells+pos] = p;
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
        int x, y, z, eje, val, i, pos, posf;
	double u_celda[5], u_sig[5], u_ant[5], p;
	for (z=1;z<U.N_z-1;z++){
	  for (y=1;y<U.N_y-1;y++){
	    for (x=1:x<U.N_x-1;x++){
	      pos = pos(x,y,z,U.N_x,U.N_y);
	      u_celda = {U.U[0*U.N_cells+pos], U.U[1*U.N_cells+pos], U.U[2*U.N_cells+pos], U.U[3*U.N_cells+pos], U.U[4*U.N_cells+pos]};
	      
	      eje = 0; // Fx
	      pos = pos(x+1,y,z,U.N_x,U.N_y); // Avanza en x
	      u_sig = {U.U[0*U.N_cells+pos], U.U[1*U.N_cells+pos], U.U[2*U.N_cells+pos], U.U[3*U.N_cells+pos], U.U[4*U.N_cells+pos]};

	      pos = pos(x-1,y,z,U.N_x,U.N_y); // Retrocede en x
	      u_ant = {U.U[0*U.N_cells+pos], U.U[1*U.N_cells+pos], U.U[2*U.N_cells+pos], U.U[3*U.N_cells+pos], U.U[4*U.N_cells+pos]};

	      for (i=0;i<5;i++){ // Promedio backward y forward
		u_ant[i] = 0.5*(u_celda[i] + u_ant[i]);
		u_sig[i] = 0.5*(u_celda[i] + u_sig[i]);
	      }

	      val = 0; // Fpmx 0
	      posf = posF(x,y,z,eje,val,U.N_x,U.N_y,U.N_z);
	      Fp[posf] = u_sig[1];
	      Fm[posf] = u_ant[1];

	      val = 1; // Fpx 1
	      p = presion(u_sig);
	      posf = posF(x,y,z,eje,val,U.N_x,U.N_y,U.N_z);
	      Fp[posf] = u_sig[1]*u_sig[1]/u_sig[0] + p;
	      val = 4; // Fpx 4
	      posf = posF(x,y,z,eje,val,U.N_x,U.N_y,U.N_z);
	      Fp[posf] = u_sig[1]*u_sig[4]/u_sig[0] + p*u_sig[1]/u_sig[0];

	      val = 1; // Fmx 1
	      p = presion(u_ant);
	      posf = posF(x,y,z,eje,val,U.N_x,U.N_y,U.N_z);
	      Fm[posf] = u_ant[1]*u_ant[1]/u_ant[0] + p;
	      val = 4; // Fmx 4
	      posf = posF(x,y,z,eje,val,U.N_x,U.N_y,U.N_z);
	      Fm[posf] = u_ant[1]*u_ant[4]/u_ant[0] + p*u_ant[1]/u_ant[0];

	      val = 2; // Fpmx 2
	      posf = posF(x,y,z,eje,val,U.N_x,U.N_y,U.N_z);
	      Fp[posf] = u_sig[1]*u_sig[2]/u_sig[0];
	      Fm[posf] = u_ant[1]*u_ant[2]/u_ant[0];
	      
	      val = 3; // Fpmx 3
	      posf = posF(x,y,z,eje,val,U.N_x,U.N_y,U.N_z);
	      Fp[posf] = u_sig[1]*u_sig[3]/u_sig[0];
	      Fm[posf] = u_ant[1]*u_ant[3]/u_ant[0];

	      eje = 1; // Fy
	      pos = pos(x,y+1,z,U.N_x,U.N_y); // Avanza en y
	      u_sig = {U.U[0*U.N_cells+pos], U.U[1*U.N_cells+pos], U.U[2*U.N_cells+pos], U.U[3*U.N_cells+pos], U.U[4*U.N_cells+pos]};

	      pos = pos(x,y-1,z,U.N_x,U.N_y); // Retrocede en y
	      u_ant = {U.U[0*U.N_cells+pos], U.U[1*U.N_cells+pos], U.U[2*U.N_cells+pos], U.U[3*U.N_cells+pos], U.U[4*U.N_cells+pos]};

	      for (i=0;i<5;i++){ // Promedio backward y forward
		u_ant[i] = 0.5*(u_celda[i] + u_ant[i]);
		u_sig[i] = 0.5*(u_celda[i] + u_sig[i]);
	      }

	      val = 0; // Fpmy 0
	      posf = posF(x,y,z,eje,val,U.N_x,U.N_y,U.N_z);
	      Fp[posf] = u_sig[2];
	      Fm[posf] = u_ant[2];

	      val = 1; // Fpmy 1
	      posf = posF(x,y,z,eje,val,U.N_x,U.N_y,U.N_z);
	      Fp[posf] = u_sig[1]*u_sig[2]/u_sig[0];
	      Fm[posf] = u_ant[1]*u_ant[2]/u_ant[0];

	      val = 2; // Fpy 2
	      p = presion(u_sig);
	      posf = posF(x,y,z,eje,val,U.N_x,U.N_y,U.N_z);
	      Fp[posf] = u_sig[2]*u_sig[2]/u_sig[0] + p;
	      val = 4; // Fpy 4
	      posf = posF(x,y,z,eje,val,U.N_x,U.N_y,U.N_z);
	      Fp[posf] = u_sig[2]*u_sig[4]/u_sig[0] + p*u_sig[2]/u_sig[0];

	      val = 2; // Fmy 2
	      p = presion(u_ant);
	      posf = posF(x,y,z,eje,val,U.N_x,U.N_y,U.N_z);
	      Fm[posf] = u_ant[2]*u_ant[2]/u_ant[0] + p;
	      val = 4; // Fmy 4
	      posf = posF(x,y,z,eje,val,U.N_x,U.N_y,U.N_z);
	      Fm[posf] = u_ant[2]*u_ant[4]/u_ant[0] + p*u_ant[2]/u_ant[0];

	      val = 3; // Fpmy 3
	      posf = posF(x,y,z,eje,val,U.N_x,U.N_y,U.N_z);
	      Fp[posf] = u_sig[2]*u_sig[3]/u_sig[0];
	      Fm[posf] = u_ant[2]*u_ant[3]/u_ant[0];

	      eje = 2; // Fz
	      pos = pos(x,y,z+1,U.N_x,U.N_y); // Avanza en z
	      u_sig = {U.U[0*U.N_cells+pos], U.U[1*U.N_cells+pos], U.U[2*U.N_cells+pos], U.U[3*U.N_cells+pos], U.U[4*U.N_cells+pos]};

	      pos = pos(x,y,z-1,U.N_x,U.N_y); // Retrocede en z
	      u_ant = {U.U[0*U.N_cells+pos], U.U[1*U.N_cells+pos], U.U[2*U.N_cells+pos], U.U[3*U.N_cells+pos], U.U[4*U.N_cells+pos]};

	      for (i=0;i<5;i++){ // Promedio backward y forward
		u_ant[i] = 0.5*(u_celda[i] + u_ant[i]);
		u_sig[i] = 0.5*(u_celda[i] + u_sig[i]);
	      }

	      val = 0; // Fpmz 0
	      posf = posF(x,y,z,eje,val,U.N_x,U.N_y,U.N_z);
	      Fp[posf] = u_sig[3];
	      Fm[posf] = u_ant[3];

	      val = 1; // Fpmz 1
	      posf = posF(x,y,z,eje,val,U.N_x,U.N_y,U.N_z);
	      Fp[posf] = u_sig[1]*u_sig[3]/u_sig[0];
	      Fm[posf] = u_ant[1]*u_ant[3]/u_ant[0];

	      val = 2; // Fpmz 2
	      posf = posF(x,y,z,eje,val,U.N_x,U.N_y,U.N_z);
	      Fp[posf] = u_sig[2]*u_sig[3]/u_sig[0];
	      Fm[posf] = u_ant[2]*u_ant[3]/u_ant[0];

	      val = 3; // Fpz 3
	      p = presion(u_sig);
	      posf = posF(x,y,z,eje,val,U.N_x,U.N_y,U.N_z);
	      Fp[posf] = u_sig[3]*u_sig[3]/u_sig[0] + p;
	      val = 4; // Fpz 4
	      posf = posF(x,y,z,eje,val,U.N_x,U.N_y,U.N_z);
	      Fp[posf] = u_sig[3]*u_sig[4]/u_sig[0] + p*u_sig[3]/u_sig[0];

	      val = 3; // Fmz 3
	      p = presion(u_ant);
	      posf = posF(x,y,z,eje,val,U.N_x,U.N_y,U.N_z);
	      Fm[posf] = u_ant[3]*u_ant[3]/u_ant[0] + p;
	      val = 4; // Fmz 4
	      posf = posF(x,y,z,eje,val,U.N_x,U.N_y,U.N_z);
	      Fm[posf] = u_ant[3]*u_ant[4]/u_ant[0] + p*u_ant[3]/u_ant[0];
	    }
	  }
	}
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

/**
 * Devuelve la presion en una celda dado como parametro u en esa celda
 */
double presion(double *u_cell){
        double p;
	p = (GAMMA - 1)*(u_cell[4] - 0.5*(u_cell[1]*u_cell[1] + u_cell[2]*u_cell[2] + u_cell[3]*u_cell[3])/u_cell[0]);
	return p;
}
