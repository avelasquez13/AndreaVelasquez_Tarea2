#include <stdio.h>
#include "struct.h"
#include "init.h"
#include "io.h"


int main(int argc, char **argv){
  physics_grid * P_state;
  U_grid * U_state;
  F_grid  *Fp, *Fm;

  P_state = create_physics_grid();
  U_state = create_U_grid();
  Fp = create_F_grid();
  Fm = create_F_grid();
   
  init_problem(P_state, U_state, Fp, Fm);

  print_L(P_state);
  
  return 0;
}
