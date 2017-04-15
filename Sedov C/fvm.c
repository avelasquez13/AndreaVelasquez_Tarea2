/*
 * fvm.c
 *
 *  Created on: Apr 11, 2017
 *      Author: felipe
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "struct.h"
#include "fvm.h"
#include "space.h"

/**
 * Calcula la evolución de la onda de choque hasta un r_final que entra por parametro
 *
 */
double evolve(physics_grid *P, U_grid *U, F_grid *Fp, F_grid *Fm, double r_final, double *radios, double *rho, double *dist, int *posiciones, int length){
	double radio=0;
	double tiempo=0;

	actualizarP(P,U);
	radio = radioChoque(P,radios,rho, dist,posiciones,length);
	while(radio<r_final){
	  fprintf(stdout, "Tiempo = %f\n", tiempo);
	  tiempo += step(P,U,Fp,Fm);
	  radio = radioChoque(P,radios,rho, dist,posiciones,length);
	}
	return tiempo;
}

/**
 *Calcula una iteración del problema (avanza un dt), devuelve el tiempo de avance
 */
double step(physics_grid *P, U_grid *U, F_grid *Fp, F_grid *Fm){
	int x, y, z, eje, val, pos, posf;
	double delta_t;
	actualizarF(U,Fp,Fm);
	delta_t = dt(P,U);
	for (z=1;z<U->N_z-1;z++){
	  for (y=1;y<U->N_y-1;y++){
	    for (x=1;x<U->N_x-1;x++){
	      pos = posi(x,y,z,U->N_x,U->N_y);
	      for (val=0;val<(NDIM+2);val++){
	    	//U->U[val*U.N_cells+pos] = U->U[val*U.N_cells+pos]; //Que quiere hacer aca?

			eje = 0;
			posf = posF(x,y,z,eje,val,U->N_x,U->N_y,U->N_z);
			U->U[val*U->N_cells+pos] += (Fm->F[posf] - Fp->F[posf])*delta_t/P->delta_x;

			eje = 1;
			posf = posF(x,y,z,eje,val,U->N_x,U->N_y,U->N_z);
			U->U[val*U->N_cells+pos] += (Fm->F[posf] - Fp->F[posf])*delta_t/P->delta_y;

			eje = 2;
			posf = posF(x,y,z,eje,val,U->N_x,U->N_y,U->N_z);
			U->U[val*U->N_cells+pos] += (Fm->F[posf] - Fp->F[posf])*delta_t/P->delta_z;
	      }
	    }
	  }
	}
	actualizarP(P,U);
	return delta_t;
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
	for (pos=0;pos<P->N_cells;pos++){ // Densidad
	  P->P[val*P->N_cells+pos] = U->U[pos];
	}

	val++;
	for (z=0;z<P->N_z;z++){ // Presion
	  for (y=0;y<P->N_y;y++){
	    for (x=0;x<P->N_x;x++){
	      pos = posi(x,y,z,P->N_x,P->N_y);
	      Ucelda(pos,U,u);
	      p = presion(u);
	      P->P[val*P->N_cells+pos] = p;
	    }
	  }
	}

	val++;
	for (z=0;z<P->N_z;z++){ // Vel x
	  for (y=0;y<P->N_y;y++){
	    for (x=0;x<P->N_x;x++){
	      pos = posi(x,y,z,P->N_x,P->N_y);
	      P->P[val*P->N_cells+pos] = U->U[(val-1)*P->N_cells+pos]/U->U[pos];
	    }
	  }
	}

	val++;
	for (z=0;z<P->N_z;z++){ // Vel y
	  for (y=0;y<P->N_y;y++){
	    for (x=0;x<P->N_x;x++){
	      pos = posi(x,y,z,P->N_x,P->N_y);
	      P->P[val*P->N_cells+pos] = U->U[(val-1)*P->N_cells+pos]/U->U[pos];
	    }
	  }
	}

	val++;
	for (z=0;z<P->N_z;z++){ // Vel z
	  for (y=0;y<P->N_y;y++){
	    for (x=0;x<P->N_x;x++){
	      pos = posi(x,y,z,P->N_x,P->N_y);
	      P->P[val*P->N_cells+pos] = U->U[(val-1)*P->N_cells+pos]/U->U[pos];
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
	for (z=1;z<U->N_z-1;z++){
	  for (y=1;y<U->N_y-1;y++){
	    for (x=1;x<U->N_x-1;x++){

	      pos = posi(x,y,z,U->N_x,U->N_y);
	      Ucelda(pos, U,u_celda);
	      
	      eje = 0; // Fx
	      pos = posi(x+1,y,z,U->N_x,U->N_y); // Avanza en x
	      Ucelda(pos, U,u_sig);

	      pos = posi(x-1,y,z,U->N_x,U->N_y); // Retrocede en x
	      Ucelda(pos, U, u_ant);

	      for (i=0;i<5;i++){ // Promedio backward y forward
	    	  u_ant[i] = 0.5*(u_celda[i] + u_ant[i]);
	    	  u_sig[i] = 0.5*(u_celda[i] + u_sig[i]);
	      }

	      val = 0; // Fpmx 0
	      posf = posF(x,y,z,eje,val,U->N_x,U->N_y,U->N_z);
	      Fp->F[posf] = u_sig[1];
	      Fm->F[posf] = u_ant[1];

	      val = 1; // Fpx 1
	      p = presion(u_sig);
	      posf = posF(x,y,z,eje,val,U->N_x,U->N_y,U->N_z);
	      Fp->F[posf] = u_sig[1]*u_sig[1]/u_sig[0] + p;
	      val = 4; // Fpx 4
	      posf = posF(x,y,z,eje,val,U->N_x,U->N_y,U->N_z);
	      Fp->F[posf] = u_sig[1]*u_sig[4]/u_sig[0] + p*u_sig[1]/u_sig[0];

	      val = 1; // Fmx 1
	      p = presion(u_ant);
	      posf = posF(x,y,z,eje,val,U->N_x,U->N_y,U->N_z);
	      Fm->F[posf] = u_ant[1]*u_ant[1]/u_ant[0] + p;
	      val = 4; // Fmx 4
	      posf = posF(x,y,z,eje,val,U->N_x,U->N_y,U->N_z);
	      Fm->F[posf] = u_ant[1]*u_ant[4]/u_ant[0] + p*u_ant[1]/u_ant[0];

	      val = 2; // Fpmx 2
	      posf = posF(x,y,z,eje,val,U->N_x,U->N_y,U->N_z);
	      Fp->F[posf] = u_sig[1]*u_sig[2]/u_sig[0];
	      Fm->F[posf] = u_ant[1]*u_ant[2]/u_ant[0];
	      
	      val = 3; // Fpmx 3
	      posf = posF(x,y,z,eje,val,U->N_x,U->N_y,U->N_z);
	      Fp->F[posf] = u_sig[1]*u_sig[3]/u_sig[0];
	      Fm->F[posf] = u_ant[1]*u_ant[3]/u_ant[0];

	      eje = 1; // Fy
	      pos = posi(x,y+1,z,U->N_x,U->N_y); // Avanza en y
	      Ucelda(pos, U,u_sig);

	      pos = posi(x,y-1,z,U->N_x,U->N_y); // Retrocede en y
	      Ucelda(pos, U,u_ant);

	      for (i=0;i<5;i++){ // Promedio backward y forward
	    	  u_ant[i] = 0.5*(u_celda[i] + u_ant[i]);
	    	  u_sig[i] = 0.5*(u_celda[i] + u_sig[i]);
	      }

	      val = 0; // Fpmy 0
	      posf = posF(x,y,z,eje,val,U->N_x,U->N_y,U->N_z);
	      Fp->F[posf] = u_sig[2];
	      Fm->F[posf] = u_ant[2];

	      val = 1; // Fpmy 1
	      posf = posF(x,y,z,eje,val,U->N_x,U->N_y,U->N_z);
	      Fp->F[posf] = u_sig[1]*u_sig[2]/u_sig[0];
	      Fm->F[posf] = u_ant[1]*u_ant[2]/u_ant[0];

	      val = 2; // Fpy 2
	      p = presion(u_sig);
	      posf = posF(x,y,z,eje,val,U->N_x,U->N_y,U->N_z);
	      Fp->F[posf] = u_sig[2]*u_sig[2]/u_sig[0] + p;
	      val = 4; // Fpy 4
	      posf = posF(x,y,z,eje,val,U->N_x,U->N_y,U->N_z);
	      Fp->F[posf] = u_sig[2]*u_sig[4]/u_sig[0] + p*u_sig[2]/u_sig[0];

	      val = 2; // Fmy 2
	      p = presion(u_ant);
	      posf = posF(x,y,z,eje,val,U->N_x,U->N_y,U->N_z);
	      Fm->F[posf] = u_ant[2]*u_ant[2]/u_ant[0] + p;
	      val = 4; // Fmy 4
	      posf = posF(x,y,z,eje,val,U->N_x,U->N_y,U->N_z);
	      Fm->F[posf] = u_ant[2]*u_ant[4]/u_ant[0] + p*u_ant[2]/u_ant[0];

	      val = 3; // Fpmy 3
	      posf = posF(x,y,z,eje,val,U->N_x,U->N_y,U->N_z);
	      Fp->F[posf] = u_sig[2]*u_sig[3]/u_sig[0];
	      Fm->F[posf] = u_ant[2]*u_ant[3]/u_ant[0];

	      eje = 2; // Fz
	      pos = posi(x,y,z+1,U->N_x,U->N_y); // Avanza en z
	      Ucelda(pos, U,u_sig);

	      pos = posi(x,y,z-1,U->N_x,U->N_y); // Retrocede en z
	      Ucelda(pos, U, u_ant);

	      for (i=0;i<5;i++){ // Promedio backward y forward
	    	  u_ant[i] = 0.5*(u_celda[i] + u_ant[i]);
	    	  u_sig[i] = 0.5*(u_celda[i] + u_sig[i]);
	      }

	      val = 0; // Fpmz 0
	      posf = posF(x,y,z,eje,val,U->N_x,U->N_y,U->N_z);
	      Fp->F[posf] = u_sig[3];
	      Fm->F[posf] = u_ant[3];

	      val = 1; // Fpmz 1
	      posf = posF(x,y,z,eje,val,U->N_x,U->N_y,U->N_z);
	      Fp->F[posf] = u_sig[1]*u_sig[3]/u_sig[0];
	      Fm->F[posf] = u_ant[1]*u_ant[3]/u_ant[0];

	      val = 2; // Fpmz 2
	      posf = posF(x,y,z,eje,val,U->N_x,U->N_y,U->N_z);
	      Fp->F[posf] = u_sig[2]*u_sig[3]/u_sig[0];
	      Fm->F[posf] = u_ant[2]*u_ant[3]/u_ant[0];

	      val = 3; // Fpz 3
	      p = presion(u_sig);
	      posf = posF(x,y,z,eje,val,U->N_x,U->N_y,U->N_z);
	      Fp->F[posf] = u_sig[3]*u_sig[3]/u_sig[0] + p;
	      val = 4; // Fpz 4
	      posf = posF(x,y,z,eje,val,U->N_x,U->N_y,U->N_z);
	      Fp->F[posf] = u_sig[3]*u_sig[4]/u_sig[0] + p*u_sig[3]/u_sig[0];

	      val = 3; // Fmz 3
	      p = presion(u_ant);
	      posf = posF(x,y,z,eje,val,U->N_x,U->N_y,U->N_z);
	      Fm->F[posf] = u_ant[3]*u_ant[3]/u_ant[0] + p;
	      val = 4; // Fmz 4
	      posf = posF(x,y,z,eje,val,U->N_x,U->N_y,U->N_z);
	      Fm->F[posf] = u_ant[3]*u_ant[4]/u_ant[0] + p*u_ant[3]/u_ant[0];
	    }
	  }
	}
}

/**
 * Encuentra la posición de la onda de choque
 */
double radioChoque(physics_grid *P, double *radios, double *rho, double *dist, int *posiciones, int length){
	double r, der, dermax, *pres;
	int i;
	int max=0;

	if(!(pres = malloc((P->N_cells/8)*sizeof(double)))){
	  fprintf(stderr, "Problem with data allocation\n");fflush(stdout);
	  exit(0);
	}

	perfilRadial(P,radios,posiciones,length,rho,pres);

	dermax = 0;
	der = 0;
	for (i=1;i<length;i++){
	  der = (pres[i] - pres[i-1])/(dist[i] - dist[i-1]);
	  if (dermax<der){
	    dermax = der;
	    max = i;
	  }
	}
	free(pres);

	r = dist[max];
	return r;
}

/**
 * Actualiza los valores de entalpía
 */
void h(U_grid *U, double* h){
        int pos;
	double p;
	double u_celda[5];
	for (pos=0;pos<U->N_cells;pos++){
	  Ucelda(pos, U, u_celda);
	  p = presion(u_celda);
	  h[pos] = (u_celda[4] + p)/u_celda[0];
	}
}

/**
 * Actualiza los valores de cs
 */
void cs(U_grid *U, double* cs){
        double *h_act;
	int pos;

	if(!(h_act = malloc(U->N_cells*sizeof(double)))){
	  fprintf(stderr, "Problem with data allocation\n");fflush(stdout);
	  exit(0);
	}
	h(U,h_act);
	for (pos=0;pos<U->N_cells;pos++){
	  cs[pos] = sqrt((GAMMA + 1)*h_act[pos]);
	}
	free(h_act);
}

/**
 * Encuentra el valor de (u, v, ó w)+cs máximo
 */
double vmax(physics_grid *P, U_grid *U){
        double vmax, vel, *cs_act;
	int pos,eje;

	if(!(cs_act = malloc(U->N_cells*sizeof(double)))){
	  fprintf(stderr, "Problem with data allocation\n");fflush(stdout);
	  exit(0);
	}
	cs(U,cs_act);

	vmax = P->P[2*P->N_cells] + cs_act[0];
	for (pos=0;pos<P->N_cells;pos++){
	  for (eje=0;eje<3;eje++){
	    vel = P->P[(2+eje)*P->N_cells+pos] + cs_act[pos];
	    if (vel > vmax){
	      vmax = vel;
	    }
	  }
	}
	free(cs_act);
	return vmax;
}

/**
 * Encuentra el dt apropiado para las condiciones actuales
 */
double dt(physics_grid *P, U_grid *U){
	double dt, vm;
	vm = vmax(P,U);
	dt = 0.5*(P->delta_x/vm);
	return dt;
}

/**
 * Devuelve un vector con la propiedad deseada del physics_grid
 * @param int propiedad toma los valores de las constantes:
 * RHO= 0
 * PRESSURE= 1
 * VX= 2
 * VY= 3
 * VZ= 4
 */
void propiedad(physics_grid *P, int propiedad, double *prop){
	int pos;
	for (pos=0;pos<P->N_cells;pos++){
	  prop[pos] = P->P[propiedad*P->N_cells+pos];
	}
}

/**
 * Devuelve la presion en una celda dado como parametro el vector U en esa celda
 */
double presion(double *u_cell){
        double p;
	p = (GAMMA - 1)*(u_cell[4] - 0.5*(u_cell[1]*u_cell[1] + u_cell[2]*u_cell[2] + u_cell[3]*u_cell[3])/u_cell[0]);
	return p;
}

/**
 * Evalua y produce una lista de promedio radial de densidad y presion
 */
void perfilRadial(physics_grid *P, double *radios, int *posiciones, int length, double *dens, double *pres){
	int i, j, pos;
	double dens_m, pres_m, *dens_ord, *pres_ord;

	int nCells=P->N_cells;

	if(!(dens_ord = malloc(nCells*sizeof(double)))){
	  fprintf(stderr, "Problem with data allocation\n");fflush(stdout);
	  exit(0);
	}
	if(!(pres_ord = malloc(nCells*sizeof(double)))){
	  fprintf(stderr, "Problem with data allocation\n");fflush(stdout);
	  exit(0);
	}
	
	for (i=0;i<nCells;i++){ // Lista de dens y pres ordenadas por radio
	  dens_ord[i] = P->P[0*nCells+posiciones[i]];
	  pres_ord[i] = P->P[1*nCells+posiciones[i]];
	}

	pos = 0;
	for (j=0;j<length;j++){ // Listas de promedios radiales
	  dens_m = pres_m = 0;
	  if (pos<nCells){
	    dens_m += dens_ord[pos];
	    pres_m += pres_ord[pos];
	  }

	  i = 1;
	  while ((pos + 1)<nCells && radios[pos] == radios[pos+1]){
	    pos++;
	    i++;
	    dens_m += dens_ord[pos];
	    pres_m += pres_ord[pos];
	  }

	  dens[j] = dens_m/i;
	  pres[j] = pres_m/i;
	  pos++;
	}
	free(dens_ord);
	free(pres_ord);
}

/**
 * Retorna el vector U para una celda especifica como una lista de 5 doubles
 */
void Ucelda(int pos, U_grid* U_act, double* u)
{
	int ncells=U_act->N_cells;
	u[0]= U_act->U[0*ncells+pos];
	u[1]= U_act->U[1*ncells+pos];
	u[2]= U_act->U[2*ncells+pos];
	u[3]= U_act->U[3*ncells+pos];
	u[4]= U_act->U[4*ncells+pos];
}
