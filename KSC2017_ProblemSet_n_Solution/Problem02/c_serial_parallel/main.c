#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <time.h>
#include "blackhole_lab.h"

void ray_trace(int w, int h, int *status, double *y, double *z)
{
    double value[C_NUMBER_OF_FUNCTIONS];
    
    // set initial value
    for (int i = 0; i < C_NUMBER_OF_FUNCTIONS; i++) {
        value[init_order[i]] = initial_value_function[init_order[i]](w, h, value);
    }
    // time integration
    for (int i = 0; i < C_TOTAL_STEP && !bool_stop_condition(value); i++) {
        double value_temp[C_NUMBER_OF_FUNCTIONS];
        double derivative[C_RK_ORDER][C_NUMBER_OF_FUNCTIONS];
        for (int k = 0; k <= C_RK_ORDER; k++) {
            for (int j = 0; j < C_NUMBER_OF_FUNCTIONS; j++) {
                value_temp[j] = value[j];
            }
            for (int l = 0; l < k; l++) {
                for (int j = 0; j < C_NUMBER_OF_FUNCTIONS; j++) {
                    value_temp[j] += RK_factor[k - 1][l] * derivative[l][j] * C_DELTA_TIME;
                }
            }
            if (C_RK_ORDER == k) {
                for (int j = 0 ;j < C_NUMBER_OF_FUNCTIONS; j++) {
                    value[j] = value_temp[j];
                }
            } else {
                for (int j = 0; j < C_NUMBER_OF_FUNCTIONS; j++) {
                    derivative[k][j] = derivative_function[j](value_temp);
                }
            }
        }
    }
    // information return
    if (bool_cross_picture(value)) {
        (*status) = HIT;
        (*y) = get_y(value) / C_l;
        (*z) = get_z(value) / C_l;
    } else if (bool_near_horizon(value)) {
        (*status) = FALL;
    } else if (bool_outside_boundary(value)) {
        (*status) = OUTSIDE;
    } else {
        (*status) = YET;
    }
}

void run(int w_start, int w_end, int h_start, int h_end, int *status, double *y, double *z)
{
    for (int h = h_start; h <= h_end; h++) {
        for (int w = w_start; w <= w_end; w++) {
            int k = h * C_RESOLUTION_WIDTH + w;
            ray_trace(w, h, status + k, y + k, z + k);
        }
    }
}

void export(int *status, double *y, double *z)
{
    FILE *data_fp = fopen(C_OUTPUT_FILENAME, "wb");
    
    int resolution[2] = {C_RESOLUTION_WIDTH, C_RESOLUTION_HEIGHT};
    fwrite(resolution, sizeof(int), 2, data_fp);
    
    fwrite(status, sizeof(int   ), C_RESOLUTION_TOTAL_PIXELS, data_fp);
    fwrite(     y, sizeof(double), C_RESOLUTION_TOTAL_PIXELS, data_fp);
    fwrite(     z, sizeof(double), C_RESOLUTION_TOTAL_PIXELS, data_fp);
    
    fclose(data_fp);
}

int main()
{
    int    *status = malloc(sizeof(int   )*C_RESOLUTION_TOTAL_PIXELS);
    double      *y = malloc(sizeof(double)*C_RESOLUTION_TOTAL_PIXELS);
    double      *z = malloc(sizeof(double)*C_RESOLUTION_TOTAL_PIXELS);

    int w_start = 0;
    int w_end   = C_RESOLUTION_WIDTH  - 1;
    int h_start = 0;
    int h_end   = C_RESOLUTION_HEIGHT - 1;
    
    clock_t start = clock();
    run(w_start, w_end, h_start, h_end, status, y, z);
    printf("Elapsed time: %fs\n", (clock() - start) / ((double)CLOCKS_PER_SEC));
    
    export(status, y, z);
    
    free(status);
    free(y);
    free(z);
    
    return 0;
}
