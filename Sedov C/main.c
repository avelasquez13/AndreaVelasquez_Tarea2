#include <stdio.h>
#include "struct.h"
#include "init.h"
#include "io.h"
#include "fvm.h"


int main(int argc, char **argv){
  physics_grid * P_state;
  U_grid * U_state;
  F_grid  *Fp, *Fm;
  double *radios; // radios[i] = radio de la celda en posiciones[i]
  double *dist, *rho; // Guardan valores de distancia (sin repetir) y promedio de rho a esa distancia
  double tiempo;
  int *posiciones; // Guarda las posiciones ordenadas ascendentemente por radio
  int length; // Tamaño del SET de radios
  int i;

  P_state = create_physics_grid();
  U_state = create_U_grid();
  Fp = create_F_grid();
  Fm = create_F_grid();
   
  init_problem(P_state, U_state, Fp, Fm);
  length = init_radios(P_state, radios, dist, rho, posiciones);

  init_conditions(U_state);

  print_list(dist,length);//Imprime radios?

  //Imprime los valores de rho para las posiciones de la onda de choque 10m,60m, y 120m
  double pos[3]={10,60,120};
  for (i=0;i<3;i++){
	  tiempo = evolve(P_state, U_state, Fp, Fm, pos[i], radios, rho, dist, posiciones, length);
	  print_list(rho,length);
  }


  print_L(P_state); // Supongo que esta linea sobra
  
  return 0;
}
