#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>
#include <omp.h>

#define MASTER      0                  /* taskid of first process */
#define NONE        0                  /* indicates no neighbor */
#define BEGIN       1                  /* message tag */
#define LTAG        2                  /* message tag */
#define RTAG        3                  /* message tag */
#define TERMINATION    4                    /* message tag */
#define CONVERGE    5
#define OFFSET        6
#define ROW            7
#define RESULT        8                          /* message tag */
#define MAXITER     10              /* maximum iterations to go through */


double **alloc_2d_double(int rows, int cols);

int dowork(double **data, double **local_data, int offset, int rows, int cols,int rank, double epsilon);
void print_2d_array(double **A, int rows, int cols);
void write_output(int size, double **result, char *output_file);
void master(double **data, int size, int numworkers, MPI_Datatype rowtype, char *output_file);
void worker(double **data, int size, double epsilon, int numworkers, int rank,
            MPI_Datatype rowtype);

int main(int argc, char *argv[]) {
    int world_size, rank, numworkers;


    MPI_Datatype rowtype;
    int size = 0;
    int i, j;
    double **data;
    double epsilon;
    double start = omp_get_wtime();
    char *output_file;

    // Initializes the MPI execution environment
    int provided;
    MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    numworkers = world_size - 1;
    MPI_Request recv_request;

    char str[10];
    if (argc != 4) {
        printf("Incorrect arguments. Should have one argument for size, one for epsilon and one for output.");
        exit(EXIT_FAILURE);
    }
    size = atoi(argv[1]);
    epsilon = atof(argv[2]);
    output_file = argv[3];

    // Initialize the array to store the original data				
    data = alloc_2d_double(size, size);

    // Initialize the input
    for (i = 0; i < size; i++) {
        if (i < size / 5) {
            for (j = 0; j < size; j++) {
                data[i][j] = 100;
            }
        } else {
            for (j = 0; j < size; j++) {
                if (j > size / 5 * 2 && j < size / 5 * 3) {
                    data[i][j] = 100;
                } else {
                    data[i][j] = 0;
                }
            }
        }
    }

    MPI_Type_vector(1, size, 1, MPI_DOUBLE, &rowtype);
    MPI_Type_commit(&rowtype);


    /************************* master code *******************************/
    if (rank == MASTER) {
        master(data, size, numworkers, rowtype, output_file);
    }
    /************************* workers code **********************************/
    if (rank != MASTER) {
        worker(data, size, epsilon, numworkers, rank, rowtype);
    }

    MPI_Finalize();
    printf("Time: %f\n", omp_get_wtime() - start);

    return 0;
}

void master(double **data, int size, int numworkers, MPI_Datatype rowtype, char *output_file) {
    int source, msgtype, msg, rows, offset, i;
#pragma omp parallel private(rows,offset,source,msgtype,msg) shared(data)
    {
        /* wait for results from all worker tasks */
#pragma omp for
        for (i=1; i<=numworkers; i++)
        {
            source = i;
            MPI_Recv(&offset, 1, MPI_INT, source, OFFSET, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Recv(&rows, 1, MPI_INT, source, ROW, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Recv(&data[offset][0], rows, rowtype, source,
                     RESULT, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
    }
    /* All threads join master thread */
    /* Write final output */
    write_output(size, data, output_file);
}

void worker(double **data, int size, double epsilon, int numworkers, int rank,
            MPI_Datatype rowtype) {
    int status,left_status,right_status;    // status 0: not converged
                                            // status 1: converged, but termination condition not met
                                            // status 2: termination condition met
    int averow, extra, rows, offset, left, right, i, j;
    status = 0;
    averow = size / numworkers;
    extra = size % numworkers;
    rows = (rank <= extra) ? averow + 1 : averow;
    if (rank <= extra) {
        rows = averow + 1;
        offset = (rank - 1) * rows;
        //end_row = start_row+rows-1;
    } else {
        rows = averow;
        offset = extra * (averow + 1) + (rank - extra) * averow - rows;
    }

    /* Find out the left, right neighbours*/
    if (rank == 1)
        left = NONE;
    else
        left = rank - 1;
    if (rank == numworkers)
        right = NONE;
    else
        right = rank + 1;

    left_status = left ? 0 : 2;
    right_status = right? 0 : 2;

    /* Allocate memory for local data */
    double **local_data;
    /* my partition, plus two neighbouring rows */
    local_data = alloc_2d_double(rows + 2, size);
    /* initialise the two neighbouring rows to zero */
    for (i = 0; i < rows + 2; i++) {
        if (i == 0 || i == rows + 1) {
            for (j = 0; j < size; j++) {
                local_data[i][j] = 0;
            }
        } else {
            for (j = 0; j < size; j++) {
                local_data[i][j] = data[offset + i - 1][j];
            }
        }
    }

    printf("rank: %d, offset: %d,rows:%d\n", rank, offset, rows);
    int num_iter = 0;
    int signal = 0;

    MPI_Buffer_attach(malloc(4 * size * sizeof(double)), 4 * size * sizeof(double));
    while (1) {
        num_iter++;

        // For the nodes with even rank, first communicate with left then right
        if(rank % 2==0){
            if (left != NONE) {
                if (left_status!=2) {
                    MPI_Bsend(&local_data[1][0], 1, rowtype, left, RTAG, MPI_COMM_WORLD);
                    MPI_Recv(&local_data[0][0], 1, rowtype, left, LTAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                }
            }
            if (right != NONE) {
                if (right_status!=2) {
                    MPI_Bsend(&local_data[rows][0], 1, rowtype, right, LTAG, MPI_COMM_WORLD);
                    MPI_Recv(&local_data[rows + 1][0], 1, rowtype, right, RTAG, MPI_COMM_WORLD,
                             MPI_STATUS_IGNORE);
                }
            }
        } else {
            // For the nodes with odd rank, first communicate with right then left
            if (right != NONE) {
                if (right_status!=2) {
                    MPI_Bsend(&local_data[rows][0], 1, rowtype, right, LTAG, MPI_COMM_WORLD);
                    MPI_Recv(&local_data[rows + 1][0], 1, rowtype, right, RTAG, MPI_COMM_WORLD,
                             MPI_STATUS_IGNORE);
                }
            }
            if (left != NONE) {
                if (left_status!=2) {
                    MPI_Bsend(&local_data[1][0], 1, rowtype, left, RTAG, MPI_COMM_WORLD);
                    MPI_Recv(&local_data[0][0], 1, rowtype, left, LTAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                }
            }
        }


        //only update the first to the second-last row in the local data
        if(dowork(data, local_data, offset, rows, size, rank, epsilon) && !status ){
            status = 1;
        }

        if(status && left_status && right_status){
            status = 2;
        }

        // Send to neighbour when converged
        if(rank%2==0){
            if (left != NONE && left_status!=2){
                MPI_Send(&status, 1, MPI_INT, left, CONVERGE, MPI_COMM_WORLD);
                MPI_Recv(&left_status, 1, MPI_INT, left, CONVERGE, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
            if (right != NONE && right_status!=2){
                MPI_Send(&status, 1, MPI_INT, right, CONVERGE, MPI_COMM_WORLD);
                MPI_Recv(&right_status, 1, MPI_INT, right, CONVERGE, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
        }else{
            if (right != NONE && right_status!=2){
                MPI_Send(&status, 1, MPI_INT, right, CONVERGE, MPI_COMM_WORLD);
                MPI_Recv(&right_status, 1, MPI_INT, right, CONVERGE, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
            if (left != NONE && left_status!=2){
                MPI_Send(&status, 1, MPI_INT, left, CONVERGE, MPI_COMM_WORLD);
                MPI_Recv(&left_status, 1, MPI_INT, left, CONVERGE, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
        }

        if (status==2) {
            break;
        }

    }

    /* Finally, send my portion of final results back to master */
    MPI_Send(&offset, 1, MPI_INT, MASTER, OFFSET, MPI_COMM_WORLD);
    MPI_Send(&rows, 1, MPI_INT, MASTER, ROW, MPI_COMM_WORLD);
    MPI_Send(&local_data[1][0], rows, rowtype, MASTER, RESULT, MPI_COMM_WORLD);

}

double **alloc_2d_double(int rows, int cols) {
    double *data = (double *) malloc(rows * cols * sizeof(double));
    double **array = (double **) malloc(rows * sizeof(double *));
    int i;
    for (i = 0; i < rows; i++)
        array[i] = &(data[cols * i]);

    return array;
}

void print_2d_array(double **A, int rows, int cols) {
    int i, j;
    for (i = 0; i < rows; i++) {
        for (j = 0; j < cols; j++) {
            printf("%f ", A[i][j]);
        }
        printf("\n");
    }
}

int dowork(double **data, double **local_data, int offset, int rows, int cols, int rank, double epsilon) {
    double max_change, max_black, max_red;
    max_change = 1;
    max_black = 0;
    max_red = 0;
    int i, j;


#pragma omp parallel for shared (data,local_data) private (i, j) reduction(max: max_red)
    for (i = 1; i <= rows; i++) {
        // Calculate all the red points
        j = i % 2 == 0 ? 0 : 1;
        for (; j < cols; j += 2) {
            // Only calculate the point without an initial value
            if (data[offset + i - 1][j] == 0) {
                // Set the point next to the current point to be 0 if the point is on boundary
                double left, right, up, down, diff, prev;
                up = local_data[i - 1][j];
                down = local_data[i + 1][j];
                left = j - 1 < 0 ? 0 : local_data[i][j - 1];
                right = j + 1 >= cols ? 0 : local_data[i][j + 1];
                prev = local_data[i][j];
                local_data[i][j] = 0.25 * (left + right + up + down);
                diff = fabs(prev - local_data[i][j]);
                if (diff > max_red) {
                    max_red = diff;
                }
            }
        }
    }

#pragma omp parallel for shared (data,local_data) private (i, j) reduction(max: max_black)
    for (i = 1; i <= rows; i++) {
        // Calculate all the black points
        j = i % 2 == 0 ? 1 : 0;
        for (; j < cols; j += 2) {
            // Only calculate the point without an initial value
            if (data[offset + i - 1][j] == 0) {
                // Set the point next to the current point to be 0 if the point is on boundary
                double left, right, up, down, diff, prev;
                up = local_data[i - 1][j];
                down = local_data[i + 1][j];
                left = j - 1 < 0 ? 0 : local_data[i][j - 1];
                right = j + 1 >= cols ? 0 : local_data[i][j + 1];
                prev = local_data[i][j];
                local_data[i][j] = 0.25 * (left + right + up + down);
                diff = fabs(prev - local_data[i][j]);
                if (diff > max_black) {
                    max_black = diff;
                }
            }
        }
    }

    max_change = max_black > max_red ? max_black : max_red;
//    printf("rank: %d, max_change: %f\n",rank,max_change);
    if (max_change <= epsilon) {
        return 1;
    } else {
        return 0;
    }

}

void write_output(int size, double **result, char *output_file) {
    int i, j;
    FILE *output = fopen(output_file, "w");
    fprintf(output, "%d\n", size);
    for (i = 0; i < size; i++) {
        for (j = 0; j < size; j++) {
            fprintf(output, "%f ", result[i][j]);
        }
        fputs("\n", output);
    }
    fclose(output);
}