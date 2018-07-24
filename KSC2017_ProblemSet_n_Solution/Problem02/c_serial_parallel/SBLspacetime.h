#ifndef SBLspacetime_h
#define SBLspacetime_h

/* Schwarzschild spacetime in Boyer-Lindquist coordinate */
/* All quantities are unitless by scaling with M, a mass of black hole, in the geometrized unit. */

#include <stdbool.h>
#include <math.h>

/* Quantities           */
enum _FUNCTIONS {X_r, X_theta, X_phi, U_t, U_r, U_theta, U_phi, C_NUMBER_OF_FUNCTIONS};

/* Derivative functions */
extern double (* const derivative_function[C_NUMBER_OF_FUNCTIONS])(double*);

/* unit basis vectors */
double U_r_unit      (double *var);
double U_theta_unit  (double *var);
double U_phi_unit    (double *var);

/* constraint */
double U_t_constraint(double *var);


double get_x(double *var);
double get_y(double *var);
double get_z(double *var);

bool bool_near_horizon(double *var);

#endif /* SBLspacetime_h */
