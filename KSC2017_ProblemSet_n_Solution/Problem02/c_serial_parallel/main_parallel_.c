#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <time.h>
#include "blackhole_lab.h"
#include <mpi.h>

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

#define CEIL_DIV(x, y) (((x) + (y) - 1) / (y))

#define CHUNK_SIZE 5
#define NUMBER_OF_JOBS_WIDTH    CEIL_DIV(C_RESOLUTION_WIDTH , CHUNK_SIZE)
#define NUMBER_OF_JOBS_HEIGHT   CEIL_DIV(C_RESOLUTION_HEIGHT, CHUNK_SIZE)
#define NUMBER_OF_JOBS          (NUMBER_OF_JOBS_WIDTH * NUMBER_OF_JOBS_HEIGHT)

enum _MESSAGE_TAG {TAG_NEW_JOB, TAG_JOB_COMPLETED, TAG_DATA};
enum _JOB_STATUS {STAT_READY, STAT_START, STAT_COMPLETED};

int main(int argc, char **argv)
{
    int *status = malloc(sizeof(int   )*C_RESOLUTION_TOTAL_PIXELS);
    double   *y = malloc(sizeof(double)*C_RESOLUTION_TOTAL_PIXELS);
    double   *z = malloc(sizeof(double)*C_RESOLUTION_TOTAL_PIXELS);

    int w_start[NUMBER_OF_JOBS];
    int w_end[NUMBER_OF_JOBS];
    int h_start[NUMBER_OF_JOBS];
    int h_end[NUMBER_OF_JOBS];
    
    clock_t start = clock();

    for (int j = 0; j < NUMBER_OF_JOBS; j++) {
        int j_h = j / NUMBER_OF_JOBS_WIDTH;
        int j_w = j % NUMBER_OF_JOBS_WIDTH;
        w_start[j] = CHUNK_SIZE * j_w;
        w_end  [j] = (j_w == NUMBER_OF_JOBS_WIDTH  - 1) ? C_RESOLUTION_WIDTH  - 1 : CHUNK_SIZE * (j_w + 1) - 1;
        h_start[j] = CHUNK_SIZE * j_h;
        h_end  [j] = (j_h == NUMBER_OF_JOBS_HEIGHT - 1) ? C_RESOLUTION_HEIGHT - 1 : CHUNK_SIZE * (j_h + 1) - 1;
    }
    
    // mpi initialize
    int my_rank, num_procs;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    int num_working_procs = ((num_procs - 1) < NUMBER_OF_JOBS) ? num_procs - 1: NUMBER_OF_JOBS;

    if (my_rank > num_working_procs) {
        MPI_Finalize();
        return 0;
    }
    
    // data type
    MPI_Datatype chunk_int_type[NUMBER_OF_JOBS];
    MPI_Datatype chunk_double_type[NUMBER_OF_JOBS];
    for (int j = 0; j < NUMBER_OF_JOBS; j++) {
        int size_array[2] = {C_RESOLUTION_HEIGHT, C_RESOLUTION_WIDTH};
        int subsize_array[2] = {h_end[j] - h_start[j] + 1, w_end[j] - w_start[j] + 1};
        int start_array[2] = {h_start[j], w_start[j]};

        MPI_Type_create_subarray(2, size_array, subsize_array, start_array, MPI_ORDER_C, MPI_INT, chunk_int_type + j);
        MPI_Type_create_subarray(2, size_array, subsize_array, start_array, MPI_ORDER_C, MPI_DOUBLE, chunk_double_type + j);
        MPI_Type_commit(chunk_int_type + j);
        MPI_Type_commit(chunk_double_type + j);
    }

    // communicator
    MPI_Group world_group, group;
    MPI_Comm_group(MPI_COMM_WORLD, &world_group);
    int ranges[1][3] = {{0, num_working_procs, 1}};
    MPI_Group_range_incl(world_group, 1, ranges, &group);
    MPI_Comm comm;
    MPI_Comm_create_group(MPI_COMM_WORLD, group, 0, &comm);
    MPI_Group_free(&world_group);
    MPI_Group_free(&group);

    if (0 == my_rank) {
        enum _JOB_STATUS job_status[NUMBER_OF_JOBS];
        for (int j = 0; j < NUMBER_OF_JOBS; j++) {
            job_status[j] = STAT_READY;
        }
        int job_number[num_working_procs];
        MPI_Request request_send[num_working_procs], request_recv[num_working_procs];
        int initial_count = 0;
        while (1) {
            int i;
            int rank;
            if (initial_count < num_working_procs) {
                i = initial_count++;
                rank = i + 1;
            } else {
                MPI_Status status;
                MPI_Waitany(num_working_procs, request_recv, &i, &status);
                if (MPI_UNDEFINED == i) break;
                rank = status.MPI_SOURCE;
                job_status[job_number[i]] = STAT_COMPLETED;
                MPI_Request_free(request_send + i);
            }
            job_number[i] = -1;
            for (int j = 0; j < NUMBER_OF_JOBS; j++) {
                if (STAT_READY == job_status[j]) {
                    job_number[i] = j;
                    job_status[j] = STAT_START;
                    break;
                }
            }
            MPI_Isend(job_number + i, 1, MPI_INT, rank, TAG_NEW_JOB, comm, request_send + i);
            if (job_number[i] >= 0) {
                MPI_Irecv(job_number + i, 1, MPI_INT, rank, TAG_JOB_COMPLETED, comm, request_recv + i);
            }
        }
        while(1) {
            int i;
            MPI_Waitany(num_working_procs, request_send, &i, MPI_STATUS_IGNORE);
            if (MPI_UNDEFINED == i) break;
        }
        
        // receive data from slaves
        for (int j = 0; j < NUMBER_OF_JOBS; j++) {
            MPI_Recv(status, 1,    chunk_int_type[j], MPI_ANY_SOURCE, TAG_DATA + j, comm, MPI_STATUS_IGNORE);
            MPI_Recv(     y, 1, chunk_double_type[j], MPI_ANY_SOURCE, TAG_DATA + j, comm, MPI_STATUS_IGNORE);
            MPI_Recv(     z, 1, chunk_double_type[j], MPI_ANY_SOURCE, TAG_DATA + j, comm, MPI_STATUS_IGNORE);
        }
    } else {
        bool my_job[NUMBER_OF_JOBS] = {false};
        while (1) {
            // receive a new job message from the master
            int job_number;
            MPI_Recv(&job_number, 1, MPI_INT, 0, TAG_NEW_JOB, comm, MPI_STATUS_IGNORE);
            if (job_number < 0) break;
            
            // get a job domain from the job_number and run
            run(w_start[job_number], w_end[job_number], h_start[job_number], h_end[job_number], status, y, z);
            my_job[job_number] = true;
            
            // send a job completed message to the master
            MPI_Send(&job_number, 1, MPI_INT, 0, TAG_JOB_COMPLETED, comm);
        }
        for (int j = 0; j < NUMBER_OF_JOBS; j++) {
            if (my_job[j]) {
                // send data to the master
                MPI_Send(status, 1,    chunk_int_type[j], 0, TAG_DATA + j, comm);
                MPI_Send(y     , 1, chunk_double_type[j], 0, TAG_DATA + j, comm);
                MPI_Send(z     , 1, chunk_double_type[j], 0, TAG_DATA + j, comm);
            }
        }
    }

    MPI_Comm_free(&comm);
    for (int j = 0; j < NUMBER_OF_JOBS; j++) {
        MPI_Type_free(chunk_double_type + j);
    }
    MPI_Finalize();

    if (0 == my_rank) {
        printf("Elapsed time: %fs\n", (clock() - start) / ((double)CLOCKS_PER_SEC));
        
        export(status, y, z);
    }
    
    free(status);
    free(y);
    free(z);

    return 0;
}
