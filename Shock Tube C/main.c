/*
 * main.c
 *
 *  Created on: Apr 9, 2017
 *      Author: felipe
 */
#include <stdio.h>
#include <stdlib.h>
#include "estructuras.h"
#include "lax_wendoff.h"
#include "io.h"

int I;

int main(){

	I=1/dx+2;
	u* U;
	if(!(U = malloc(I*sizeof(u)))){
			    fprintf(stderr, "Problem with data allocation\n");fflush(stdout);
			    exit(0);
			}
	U=lax_wendoff(0.22, I);
	save(U, I);

	return 0;
}
