void init_problem(physics_grid *P, U_grid *U, F_grid *F_p, F_grid *F_m);
physics_grid * create_physics_grid(void);
U_grid * create_U_grid(void);
F_grid * create_F_grid(void);
double* create_listNdoubles(physics_grid *P);
int* create_listNints(physics_grid *P);
int init_radios(physics_grid *P, double *radios, double *dist, double *rho, int *posiciones);
void init_conditions(U_grid *U, physics_grid *P);
void ordenarPorRadios(double *radios, int *posiciones, int length);
