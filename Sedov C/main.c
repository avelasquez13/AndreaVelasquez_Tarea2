#include <stdio.h>
#include "struct.h"
#include "init.h"
#include "io.h"
#include "fvm.h"


int main(int argc, char **argv){
  fprintf(stdout, "---Sedov--- \n");
  physics_grid * P_state;
  U_grid * U_state;
  F_grid  *Fp, *Fm;
  double *radios; // radios[i] = radio de la celda en posiciones[i]
  double *rho; // Guarda valores de rho promedio a un radio
  int *contador;// numero de celdas con ese radio
  double tiempo = 0;
  int length; // Tamaño del SET de radios
  int i;

  P_state = create_physics_grid();
  U_state = create_U_grid();
  Fp = create_F_grid();
  Fm = create_F_grid();
  init_problem(P_state, U_state, Fp, Fm);
  printf("Size arrays: %d\n",P_state->N_cells);

	length = length_radios(P_state);
	radios=create_listNdoubles(length);
	rho=create_listNdoubles(length);
	contador=create_listNints(length);

  init_radios(P_state, radios, rho, contador, length);
  printf("Size radios: %d\n",length);

  init_conditions(U_state,P_state);
  fprintf(stdout, "Inicializo \n");
  print_list(radios,length);//Imprime radios

  //Imprime los valores de rho para las posiciones de la onda de choque 10m,60m, y 120m
  double pos[3]={10,60,120};
  for (i=0;i<3;i++){
	  tiempo += evolve(P_state, U_state, Fp, Fm, pos[i], radios, rho, contador, length);
	  print_list(rho,length);
	  print_time(tiempo);
	  fprintf(stdout, "Evoluciono %fm \n", pos[i]);
  }


  print_L(P_state); // Supongo que esta linea sobra

  return 0;
}
