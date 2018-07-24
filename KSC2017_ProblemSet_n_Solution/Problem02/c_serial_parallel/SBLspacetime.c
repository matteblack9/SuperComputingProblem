#include "SBLspacetime.h"
#include <math.h>

#define R       var[X_r]
#define R2      (R*R)
#define R3      (R*R*R)
#define THETA   var[X_theta]
#define PHI     var[X_phi]
#define S       sin(THETA)
#define S2      (S*S)
#define C       cos(THETA)
#define U(a)    var[U_t + (a)]

typedef enum _TWO_INDEX {tt, tr, th, tp, rt, rr, rh, rp, ht, hr, hh, hp, pt, pr, ph, pp} TWO_INDEX;
#define I2(a, b) (4*(a) + b)

static double g(TWO_INDEX i, double *var)
{
    switch(i) {
        case tt:
            return - (1. - 2. / R);
        case rr:
            return - 1. / g(tt, var);
        case hh:
            return R*R;
        case pp:
            return R*R*S*S;
        default:
            return 0.;
    }
}

#define G(a,b)  g(I2(a, b), var)

static double U_unit(int a, double *var)
{
    return 1. / sqrt(G(a, a));
}

double U_r_unit    (double *var)
{
    return U_unit(1, var);
}

double U_theta_unit(double *var)
{
    return U_unit(2, var);
}

double U_phi_unit  (double *var)
{
    return U_unit(3, var);
}

static double u_t(double *var)
{
    int i;
    double ret = 0.;
    for (i = 1; i <= 3; i++)
        ret += G(0, i) * U(i);
    return ret;
}

static double u_square(double *var)
{
    int i, j;
    double ret = 0.;
    for (i = 1; i <= 3; i++)
        for (j = 1; j <= 3; j++)
            ret += G(i, j) * U(i) * U(i);
    return ret;
}

double U_t_constraint(double *var)
{
    return -(
        + u_t(var)
        + sqrt(
            + u_t(var) * u_t(var)
            - G(0, 0) * u_square(var)
        )
    ) / G(0, 0);
}

static double derivative_X(int a, double *var)
{
    return U(a) / U(0);
}

static double derivative_X_r    (double *var)
{
    return derivative_X(1, var);
}

static double derivative_X_theta(double *var)
{
    return derivative_X(2, var);
}

static double derivative_X_phi  (double *var)
{
    return derivative_X(3, var);
}

typedef enum _THREE_INDEX {
    ttt, ttr, tth, ttp, trt, trr, trh, trp, tht, thr, thh, thp, tpt, tpr, tph, tpp,
    rtt, rtr, rth, rtp, rrt, rrr, rrh, rrp, rht, rhr, rhh, rhp, rpt, rpr, rph, rpp,
    htt, htr, hth, htp, hrt, hrr, hrh, hrp, hht, hhr, hhh, hhp, hpt, hpr, hph, hpp,
    ptt, ptr, pth, ptp, prt, prr, prh, prp, pht, phr, phh, php, ppt, ppr, pph, ppp
} THREE_INDEX;
#define I3(a, b, c) (16*(a) + 4*(b) + c)

static double Gamma(THREE_INDEX a, double *var)
{
    switch (a) {
        case ttr:
        case trt:
            return 1. / (R2 - 2.*R);
        case rtt:
            return (R - 2.) / R3;
        case rrr:
            return 1. / (2.*R - R2);
        case rhh:
            return 2. - R;
        case rpp:
            return (2. - R)*S2;
        case hrh:
        case hhr:
            return 1. / R;
        case hpp:
            return -C*S;
        case prp:
        case ppr:
            return 1. / R;
        case php:
        case pph:
            return C / S;
        default:
            return 0.;
    }
}

static double derivative_U(int a, double *var)
{
    int i, j;
    double ret = 0.;
    for (i = 0; i <= 3; i++)
        for (j = 0; j <= 3; j++)
            ret += Gamma(I3(a, i, j), var) * U(i) * U(j);
    return - ret / U(0);
}

static double derivative_U_t    (double *var)
{
    return derivative_U(0, var);
}

static double derivative_U_r    (double *var)
{
    return derivative_U(1, var);
}

static double derivative_U_theta(double *var)
{
    return derivative_U(2, var);
}

static double derivative_U_phi  (double *var)
{
    return derivative_U(3, var);
}

double (* const derivative_function[C_NUMBER_OF_FUNCTIONS])(double*) = {
    derivative_X_r,
    derivative_X_theta,
    derivative_X_phi,
    derivative_U_t,
    derivative_U_r,
    derivative_U_theta,
    derivative_U_phi,
};

double get_x(double *var)
{
    return R * S * cos(PHI);
}

double get_y(double *var)
{
    return R * S * sin(PHI);
}

double get_z(double *var)
{
    return R * C;
}

bool bool_near_horizon(double *var)
{
    return R <= 3.;
}
