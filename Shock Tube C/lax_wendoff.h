/*
 * lax_wendoff.h
 *
 *  Created on: Apr 9, 2017
 *      Author: felip
 */

#ifndef LAX_WENDOFF_H_
#define LAX_WENDOFF_H_

u* init(void);
f* F(u* U_n);
u* step(u* U_n, f* F_n, double dt);
double u_max(u* U);
u* lax_wendoff(double tmax, int pI);
u* u_malloc(void);

#endif /* LAX_WENDOFF_H_ */
