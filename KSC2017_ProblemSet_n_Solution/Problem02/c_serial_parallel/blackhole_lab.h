#ifndef blackhole_lab_h
#define blackhole_lab_h

#include "SBLspacetime.h"

#include <math.h>

/* Parameters for lab setting         */
#define C_PI            acos(-1.)
#define C_L1            10.                     /* distance between observer and black hole */
#define C_L2            10.                     /* distance between black hole and picture  */
#define C_l             30.                     /* length of side of picture                */
#define C_VISUAL_ANGLE  120.                    /* visual angle of vision                   */
#define C_d             (1. / (sqrt(2.) * tan(C_VISUAL_ANGLE / 2. / 180. * C_PI)))
                                                /* distance between observer and vision     */
/* Parameters for resolution of vision */
#define C_OUTPUT_FILENAME         "data"
#define C_RESOLUTION_WIDTH        100
#define C_RESOLUTION_HEIGHT       100
#define C_RESOLUTION_TOTAL_PIXELS (C_RESOLUTION_WIDTH * C_RESOLUTION_HEIGHT)

/* Parameters for time integration     */
#define C_DELTA_TIME                -0.05
#define C_TOTAL_STEP                20000
#define C_RK_ORDER                  4
extern const double * const RK_factor[C_RK_ORDER];

/* For classification of rays          */
enum _RAY_STATUS {HIT, FALL, OUTSIDE, YET};
bool bool_outside_boundary(double *var);
bool bool_cross_picture   (double *var);
bool bool_stop_condition  (double *var);

/* For initialization of quantities    */
extern double (* const initial_value_function[C_NUMBER_OF_FUNCTIONS])(int, int, double*);
extern int init_order[C_NUMBER_OF_FUNCTIONS];

#endif /* blackhole_lab_h */
