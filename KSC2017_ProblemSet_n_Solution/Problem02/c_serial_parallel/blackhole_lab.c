#include "blackhole_lab.h"

const double RK_factor1[] = {1. / 2.                           };
const double RK_factor2[] = {     0., 1. / 2.                  };
const double RK_factor3[] = {     0.,      0.,      1.         };
const double RK_factor4[] = {1. / 6., 1. / 3., 1. / 3., 1. / 6.};
const double * const RK_factor[C_RK_ORDER] = {RK_factor1, RK_factor2, RK_factor3, RK_factor4};

#define X get_x(var)
#define Y get_y(var)
#define Z get_z(var)

bool bool_outside_boundary(double *var)
{
    return X > C_L1 || X < -C_L2 || fabs(Y) > C_l || fabs(Z) > C_l;
}

bool bool_cross_picture   (double *var)
{
    return X < - C_L2 && fabs(Y) < C_l / 2. && fabs(Z) < C_l / 2.;
}

bool bool_stop_condition  (double *var)
{
    return
        bool_outside_boundary(var) ||
        bool_cross_picture   (var) ||
        bool_near_horizon    (var);
}

#define ARGS_OF_INITIAL_VALUE int w, int h, double *var

static double initial_value_X_r    (ARGS_OF_INITIAL_VALUE)
{
    return C_L1;
}

static double initial_value_X_theta(ARGS_OF_INITIAL_VALUE)
{
    return C_PI / 2.;
}

static double initial_value_X_phi  (ARGS_OF_INITIAL_VALUE)
{
    return 0.;
}

static double initial_value_U_r    (ARGS_OF_INITIAL_VALUE)
{
    return C_d * U_r_unit(var);
}

static double initial_value_U_theta(ARGS_OF_INITIAL_VALUE)
{
    return - (-0.5 + (double)h / ((double)(C_RESOLUTION_HEIGHT - 1))) * U_theta_unit(var);
}

static double initial_value_U_phi  (ARGS_OF_INITIAL_VALUE)
{
    return - (-0.5 + (double)w / ((double)(C_RESOLUTION_WIDTH  - 1))) * U_phi_unit(var);
}

static double initial_value_U_t    (ARGS_OF_INITIAL_VALUE)
{
    return U_t_constraint(var);
}

double (* const initial_value_function[C_NUMBER_OF_FUNCTIONS])(int, int, double*) = {
    initial_value_X_r,
    initial_value_X_theta,
    initial_value_X_phi,
    initial_value_U_t,
    initial_value_U_r,
    initial_value_U_theta,
    initial_value_U_phi,
};

int init_order[C_NUMBER_OF_FUNCTIONS] = {X_r, X_theta, X_phi, U_r, U_theta, U_phi, U_t};

