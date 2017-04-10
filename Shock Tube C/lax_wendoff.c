/*
 * lax_wendoff.c
 *
 *  Created on: Apr 9, 2017
 *      Author: felip
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "estructuras.h"
#include "lax_wendoff.h"

int I;
int i;

u* lax_wendoff(double tmax, int pI){
	I=pI;
	double dt=0.01;
	double t=0;

	u* U_n= init();
	u* U_n1;
	f* F_n=F(U_n);

	while(t<=tmax){
		t=t+dt;
		printf("Tiempo= %f\n",t);

		U_n1=step(U_n, F_n, dt);
		F_n=F(U_n1);
		U_n=U_n1;

		if (u_max(U_n)!=0){
			dt=0.5/u_max(U_n)*dx;
		}


	}

	return U_n;

}

u* init(void){

	u* U_0=u_malloc();

	for (i = 0; i < I/2; ++i) {
		U_0[i].U1=1;
		U_0[i].U2=0;
		U_0[i].U3=1/0.4;
	}
	for (i = I/2; i < I; ++i) {
		U_0[i].U1=0.125;
		U_0[i].U2=0;
		U_0[i].U3=0.1/0.4;
	}
	return U_0;
}

f* F(u* U_n){
	f* F_n;

	if(!(F_n = malloc(I*sizeof(f)))){
	    fprintf(stderr, "Problem with data allocation\n");fflush(stdout);
	    exit(0);
	}

	for (i = 0; i < I; ++i) {
		double P= 0.4*(U_n[i].U3-1/2*(U_n[i].U2*U_n[i].U2)/U_n[i].U1);
		F_n[i].F1= U_n[i].U2;
		F_n[i].F2= (U_n[i].U2*U_n[i].U2)/U_n[i].U1+P;
		F_n[i].F3= U_n[i].U2/U_n[i].U1*(U_n[i].U3+P);
	}

	return F_n;
}

u* step(u* U_n, f* F_n, double dt){
	u* U_n1=u_malloc();

	//Lax-Wendoff para todos los puntos intermedios
	for (i = 1; i < I-1; ++i) {
		U_n1[i].U1=	U_n[i].U1-dt/(2*dx)*(F_n[i+1].F1-F_n[i-1].F1)+1/4*(dt*dt)/(dx*dx)*((U_n[i+1].U1+U_n[i].U1)*(F_n[i+1].F1-F_n[i].F1)-(U_n[i].U1+U_n[i-1].U1)*(F_n[i].F1-F_n[i-1].F1));
		U_n1[i].U2=	U_n[i].U2-dt/(2*dx)*(F_n[i+1].F2-F_n[i-1].F2)+1/4*(dt*dt)/(dx*dx)*((U_n[i+1].U2+U_n[i].U2)*(F_n[i+1].F2-F_n[i].F2)-(U_n[i].U2+U_n[i-1].U2)*(F_n[i].F2-F_n[i-1].F2));
		U_n1[i].U3=	U_n[i].U3-dt/(2*dx)*(F_n[i+1].F3-F_n[i-1].F3)+1/4*(dt*dt)/(dx*dx)*((U_n[i+1].U3+U_n[i].U3)*(F_n[i+1].F3-F_n[i].F3)-(U_n[i].U3+U_n[i-1].U3)*(F_n[i].F3-F_n[i-1].F3));
	}

	//Reflexion en las fronteras
	U_n1[0].U1=U_n1[1].U1;
	U_n1[0].U2=-U_n1[1].U2;
	U_n1[0].U3=U_n1[1].U3;

	U_n1[I-1].U1=U_n1[I-2].U1;
	U_n1[I-1].U2=-U_n1[I-2].U2;
	U_n1[I-1].U3=U_n1[I-2].U3;

	return U_n1;
}

double u_max(u* U){
	double max=0;
	for (i = 1; i < I-1; ++i) {
		double vel= fabs(U[i].U2/U[i].U1);
		if(vel>max){
			max=vel;
		}
	}
	return max;
}

u* u_malloc(void){
	u* U;
		if(!(U = malloc(I*sizeof(u)))){
		    fprintf(stderr, "Problem with data allocation\n");fflush(stdout);
		    exit(0);
		}
		else return U;
}
