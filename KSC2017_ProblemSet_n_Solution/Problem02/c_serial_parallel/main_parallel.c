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

#define CEIL_DIV(x, y) (((x) + (y) - 1) / (y))
#define CHUNK_SIZE 20

#define NUMBER_OF_JOBS_WIDTH    CEIL_DIV(C_RESOLUTION_WIDTH , CHUNK_SIZE)
#define NUMBER_OF_JOBS_HEIGHT   CEIL_DIV(C_RESOLUTION_HEIGHT, CHUNK_SIZE)
#define NUMBER_OF_JOBS          (NUMBER_OF_JOBS_WIDTH * NUMBER_OF_JOBS_HEIGHT)

void get_job_domain(int job_number, int *w_start, int *w_end, int *h_start, int *h_end)
{
    int i = job_number / NUMBER_OF_JOBS_WIDTH;
    int j = job_number % NUMBER_OF_JOBS_WIDTH;
    
    *w_start = CHUNK_SIZE * j;
    *w_end   = (j == NUMBER_OF_JOBS_WIDTH  - 1) ? C_RESOLUTION_WIDTH  - 1 : CHUNK_SIZE * (j + 1) - 1;
    *h_start = CHUNK_SIZE * i;
    *h_end   = (i == NUMBER_OF_JOBS_HEIGHT - 1) ? C_RESOLUTION_HEIGHT - 1 : CHUNK_SIZE * (i + 1) - 1;
}

void export(int my_rank, bool *performed, int *status, double *y, double *z)
{
    MPI_File fh;
    MPI_File_open(MPI_COMM_WORLD, C_OUTPUT_FILENAME, MPI_MODE_WRONLY | MPI_MODE_CREATE, MPI_INFO_NULL, &fh);
    MPI_File_set_size(fh, sizeof(int)*2 + (sizeof(int) + 2*sizeof(double))*C_RESOLUTION_TOTAL_PIXELS);
    
    if (0 == my_rank) {
        int resolution[2] = {C_RESOLUTION_WIDTH, C_RESOLUTION_HEIGHT};
        MPI_File_write_at(fh, 0, resolution, 2, MPI_INT, MPI_STATUS_IGNORE);
    } else {
        for (int i = 0; i < NUMBER_OF_JOBS; i++) {
            if (performed[i]) {
                int w_start;
                int w_end;
                int h_start;
                int h_end;
                get_job_domain(i, &w_start, &w_end, &h_start, &h_end);
                for (int h = h_start; h <= h_end; h++) {
                    int index = h * C_RESOLUTION_WIDTH + w_start;
                    MPI_Offset offset = sizeof(int) * 2;
                    MPI_File_write_at(fh, offset + sizeof(int   ) * index, status + index, w_end - w_start + 1, MPI_INT   , MPI_STATUS_IGNORE);
                    offset += sizeof(int   ) * C_RESOLUTION_TOTAL_PIXELS;
                    MPI_File_write_at(fh, offset + sizeof(double) * index,      y + index, w_end - w_start + 1, MPI_DOUBLE, MPI_STATUS_IGNORE);
                    offset += sizeof(double) * C_RESOLUTION_TOTAL_PIXELS;
                    MPI_File_write_at(fh, offset + sizeof(double) * index,      z + index, w_end - w_start + 1, MPI_DOUBLE, MPI_STATUS_IGNORE);
                }
            }
        }
    }
    
    MPI_File_close(&fh);
}

enum _MESSAGE_TAG {TAG_NEW_JOB, TAG_JOB_COMPLETED};
enum _JOB_STATUS {STAT_READY, STAT_START, STAT_COMPLETED};

int main(int argc, char **argv)
{
    int my_rank, num_procs;
    int my_rank, num_procs;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

    // variables for the master
    clock_t start = 0;
    int num_working_procs = ((num_procs - 1) < NUMBER_OF_JOBS) ? num_procs - 1: NUMBER_OF_JOBS;
    MPI_Request request[num_working_procs + 1];

    // variables for slaves
    bool performed[NUMBER_OF_JOBS] = {false};
    int *status = NULL;
    double   *y = NULL;
    double   *z = NULL;

    if (0 == my_rank) {
        start = clock();
        enum _JOB_STATUS job_status[NUMBER_OF_JOBS];
        for (int i = 0; i < NUMBER_OF_JOBS; i++) {
            job_status[i] = STAT_READY;
        }
        
        // isend new job messages to all slaves
        for (int i = 0; i < num_working_procs; i++) {
            int rank = i + 1;
            MPI_Isend(&i, 1, MPI_INT, rank, TAG_NEW_JOB, MPI_COMM_WORLD, request + rank);
            //send jobs to all processor and request that is MPI_WAIT or MPI TEST etc...
            job_status[i] = STAT_START;
        }
        
        int number_of_remaining_jobs = NUMBER_OF_JOBS - num_working_procs;
        while (number_of_remaining_jobs > 0) {
            int job_number;
            
            // receive a job completed message from any slave
            MPI_Status mpi_status;
            MPI_Recv(&job_number, 1, MPI_INT, MPI_ANY_SOURCE, TAG_JOB_COMPLETED, MPI_COMM_WORLD, &mpi_status);
            int rank = mpi_status.MPI_SOURCE;
            MPI_Request_free(request + rank);
            job_status[job_number] = STAT_COMPLETED;

            // isend a new job message to the above slave
            for (int i = 0; i < NUMBER_OF_JOBS; i++) {
                if (STAT_READY == job_status[i]) {
                    job_number = i;
                    break;
                }
            }
            job_status[job_number] = STAT_START;
            number_of_remaining_jobs--;
            MPI_Isend(&job_number, 1, MPI_INT, rank, TAG_NEW_JOB, MPI_COMM_WORLD, request + rank);
        }

        for (int i = 1; i <= num_working_procs; i++) {
            int job_number;
            // receive a last job completed message from all slaves
            MPI_Status mpi_status;
            MPI_Recv(&job_number, 1, MPI_INT, MPI_ANY_SOURCE, TAG_JOB_COMPLETED, MPI_COMM_WORLD, &mpi_status);
            int rank = mpi_status.MPI_SOURCE;
            MPI_Request_free(request + rank);
            job_status[job_number] = STAT_COMPLETED;
            
            // isend a no job message to the above slave
            job_number = -1;
            MPI_Isend(&job_number, 1, MPI_INT, i, TAG_NEW_JOB, MPI_COMM_WORLD, request + rank);
        }
    } else if (my_rank <= num_working_procs) {
        status = malloc(sizeof(int   )*C_RESOLUTION_TOTAL_PIXELS);
        y      = malloc(sizeof(double)*C_RESOLUTION_TOTAL_PIXELS);
        z      = malloc(sizeof(double)*C_RESOLUTION_TOTAL_PIXELS);

        while (1) {
            // receive a new job message from the master
            int job_number;
            MPI_Recv(&job_number, 1, MPI_INT, 0, TAG_NEW_JOB, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            if (job_number < 0) break;
            
            // get a job domain from the job_number and run
            int w_start;
            int w_end;
            int h_start;
            int h_end;
            get_job_domain(job_number, &w_start, &w_end, &h_start, &h_end);
            run(w_start, w_end, h_start, h_end, status, y, z);
            performed[job_number] = true;
            
            // send a job completed message to the master
            MPI_Send(&job_number, 1, MPI_INT, 0, TAG_JOB_COMPLETED, MPI_COMM_WORLD);
        }
    }
    
    MPI_Barrier(MPI_COMM_WORLD);

    if (0 == my_rank) {
        for (int rank = 1; rank <= num_working_procs; rank++) {
            MPI_Request_free(request + rank);
        }
    }

    export(my_rank, performed, status, y, z);

    MPI_Finalize();

    if (0 == my_rank) {
        printf("Elapsed time: %fs\n", (clock() - start) / ((double)CLOCKS_PER_SEC));
    } else if (my_rank <= num_working_procs) {
        free(status);
        free(y);
        free(z);
    }
    
    return 0;
}
