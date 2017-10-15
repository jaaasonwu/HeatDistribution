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
#define EPSILON    1e-2


double **alloc_2d_double(int rows, int cols);

int dowork(double **data, double **local_data, int offset, int rows, int cols,int rank);

void print_2d_array(double **A, int rows, int cols);

void write_output(int size, double **result);

int main(int argc, char *argv[]) {
    int world_size, rank,
            numworkers,                    /* number of worker processes */
    averow, rows, extra,
            offset,                        /* for sending rows of data */
            dest, source,               /* to - from for message send-receive */
            left, right,                /* neighbor tasks */
            msgtype,                  /* for message types */
    msg;

    MPI_Datatype rowtype;
    int size = 0;
    int i, j;
    double **data;

    // Initializes the MPI execution environment
    int provided;
    //printf("BEGIN\n");
    MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided);
    //printf("INIT_THEAD_DONE\n");
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    numworkers = world_size - 1;
    MPI_Status Stat;
    MPI_Request recv_request;

    char str[10];
    FILE *input = fopen("input.dat", "r");
    if (!input) {
        printf("File read error");
        exit(EXIT_FAILURE);
    }
    fgets(str, sizeof(str), input);
    size = atoi(str);

    // Initialize the array to store the original data				
    data = alloc_2d_double(size, size);

    // Read the file
    for (i = 0; i < size; i++) {
        for (j = 0; j < size; j++) {
            fscanf(input, "%lf ", &(data[i][j]));
        }
    }
    fclose(input);

    averow = size / numworkers;
    extra = size % numworkers;
    MPI_Type_vector(1, size, 1, MPI_DOUBLE, &rowtype);
    MPI_Type_commit(&rowtype);

    //MPI_Bcast(&size,1, MPI_INT ,0, MPI_COMM_WORLD);
    //printf("rank: %d, size: %d\n",rank,size);


    /************************* master code *******************************/
    if (rank == MASTER) {
        printf("Starting mpi_heat2D with %d worker tasks.\n", numworkers);

		#pragma omp parallel private(rows,offset,source,msgtype,msg,Stat)
        {
            #pragma omp for
            for (i = 1; i<=numworkers; i++)
            {
                source = i;
                msgtype = TERMINATION;
                MPI_Recv(&msg, 1, MPI_INT, source, msgtype, MPI_COMM_WORLD, &Stat);
                printf("%d\n", i);
            }

			/* wait for results from all worker tasks */
			#pragma omp for
            for (i=1; i<=numworkers; i++)
            {
                source = i;
                MPI_Recv(&offset, 1, MPI_INT, source, OFFSET, MPI_COMM_WORLD, &Stat);
                //printf("A\n");
                MPI_Recv(&rows, 1, MPI_INT, source, ROW, MPI_COMM_WORLD, &Stat);
                //printf("B\n");
                MPI_Recv(&data[offset][0], rows, rowtype, source,
                  RESULT, MPI_COMM_WORLD, &Stat);
                //printf("%d done\n",i);
            }
        }
		/* All threads join master thread */
		/* Write final output */
        write_output(size, data);

    }
    /************************* workers code **********************************/
    if (rank != MASTER) {
        int status,left_status,right_status;    // status 0: not converged
                                                // status 1: converged, but termination condition not met
                                                // status 2: termination condition met
        status = 0;
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
        // if(rank==1){
        // 	print_2d_array(local_data,rows+2,size);
        // }

        // if (left != NONE)
        // {
        // 		// send the first row of the worker's partition to the left neighbour
        // 	MPI_Bsend(&local_data[1][0], 1,rowtype, left, RTAG, MPI_COMM_WORLD);
        // 	source = left;
        // 	msgtype = LTAG;
        // 		// receive from left neighbour
        // 	MPI_Recv(&local_data[0][0], 1, rowtype, source,msgtype, MPI_COMM_WORLD, &Stat);
        // }
        // if (right != NONE)
        // {
        // 	MPI_Bsend(&local_data[rows][0], 1,rowtype, right, LTAG, MPI_COMM_WORLD);
        // 	source = right;
        // 	msgtype = RTAG;
        // 	MPI_Recv(&local_data[rows+1][0], 1, rowtype, source,msgtype, MPI_COMM_WORLD, &Stat);
        // }

        // if (rank==1)
        // {
        // 	print_2d_array(local_data,rows+2,size);
        // }


        int num_iter = 0;
        int signal = 0;

        MPI_Buffer_attach(malloc(4 * size * sizeof(double)), 4 * size * sizeof(double));
        while (1) 
        {
            num_iter++;
            printf("rank: %d, iter: %d\n", rank,num_iter);

            if(rank%2==0){
                if (left != NONE) {                   
                    // if(rank==2){
                    //     printf("rank 1,2 %d %d,printed by %d\n",left_converged,converged,rank);
                    // }
                    if (left_status!=2) {
                        MPI_Bsend(&local_data[1][0], 1, rowtype, left, RTAG, MPI_COMM_WORLD);                       
                        MPI_Recv(&local_data[0][0], 1, rowtype, left, LTAG, MPI_COMM_WORLD, &Stat); 
                    }
                }
                if (right != NONE) {                                 
                    if (right_status!=2) {
                        MPI_Bsend(&local_data[rows][0], 1, rowtype, right, LTAG, MPI_COMM_WORLD); 
                        MPI_Recv(&local_data[rows + 1][0], 1, rowtype, right, RTAG, MPI_COMM_WORLD,
                         &Stat);        
                    }                      
                }
            }else
            {
                if (right != NONE) {                                 
                    // if(rank==1){
                    //     printf("rank 1,2 %d %d,printed by %d\n",converged,right_converged,rank);
                    // }
                    if (right_status!=2) {
                        MPI_Bsend(&local_data[rows][0], 1, rowtype, right, LTAG, MPI_COMM_WORLD); 
                        MPI_Recv(&local_data[rows + 1][0], 1, rowtype, right, RTAG, MPI_COMM_WORLD,
                         &Stat);        
                    }                      
                }
                if (left != NONE) {                   
                    if (left_status!=2) {
                        MPI_Bsend(&local_data[1][0], 1, rowtype, left, RTAG, MPI_COMM_WORLD);                       
                        MPI_Recv(&local_data[0][0], 1, rowtype, left, LTAG, MPI_COMM_WORLD, &Stat); 
                    }
                }
            }

            // if (left != NONE) {
            //     // send the first row of the worker's partition to the left neighbour
            //     MPI_Bsend(&local_data[1][0], 1, rowtype, left, RTAG, MPI_COMM_WORLD);
            //     source = left;
            //     msgtype = LTAG;
            //     MPI_Bsend(&converged, 1, MPI_INT, source, CONVERGE, MPI_COMM_WORLD);
            //     if (!(left_converged == 1 && converged == 1)) {
            //         MPI_Recv(&local_data[0][0], 1, rowtype, source, msgtype, MPI_COMM_WORLD, &Stat);
            //         MPI_Recv(&left_converged, 1, MPI_INT, source, CONVERGE, MPI_COMM_WORLD, &Stat);
            //     }
            // }
            // if (right != NONE) {
            //     MPI_Bsend(&local_data[rows][0], 1, rowtype, right, LTAG, MPI_COMM_WORLD);
            //     source = right;
            //     msgtype = RTAG;
            //     MPI_Bsend(&converged, 1, MPI_INT, source, CONVERGE, MPI_COMM_WORLD);
            //     if (!(right_converged == 1 && converged == 1)) {
            //         MPI_Recv(&local_data[rows + 1][0], 1, rowtype, source, msgtype, MPI_COMM_WORLD,
            //          &Stat);
            //         MPI_Recv(&right_converged, 1, MPI_INT, source, CONVERGE, MPI_COMM_WORLD, &Stat);
            //     }
            // }


            //only update the first to the second-last row in the local data            
            if(dowork(data, local_data, offset, rows, size,rank) && !status ){
                status = 1;
                //signal = 1;
            }

            if(status && left_status && right_status){
                status = 2;
            }

            if(rank%2==0){
                if (left != NONE && left_status!=2){
                    MPI_Send(&status, 1, MPI_INT, left, CONVERGE, MPI_COMM_WORLD);
                    MPI_Recv(&left_status, 1, MPI_INT, left, CONVERGE, MPI_COMM_WORLD, &Stat);
                }
                if (right != NONE && right_status!=2){
                    MPI_Send(&status, 1, MPI_INT, right, CONVERGE, MPI_COMM_WORLD);
                    MPI_Recv(&right_status, 1, MPI_INT, right, CONVERGE, MPI_COMM_WORLD, &Stat);
                }
            }else{
                if (right != NONE && right_status!=2){
                    MPI_Send(&status, 1, MPI_INT, right, CONVERGE, MPI_COMM_WORLD);
                    MPI_Recv(&right_status, 1, MPI_INT, right, CONVERGE, MPI_COMM_WORLD, &Stat);
                }
                if (left != NONE && left_status!=2){
                    MPI_Send(&status, 1, MPI_INT, left, CONVERGE, MPI_COMM_WORLD);
                    MPI_Recv(&left_status, 1, MPI_INT, left, CONVERGE, MPI_COMM_WORLD, &Stat);
                }
            }

            if (status==2) {
                //printf("to send TERMINATION %d\n",rank );
                MPI_Ssend(&status, 1, MPI_INT, MASTER, TERMINATION, MPI_COMM_WORLD);
                break;
            } 
            
        }
        //printf("%d break\n",rank );

        /* Finally, send my portion of final results back to master */
        MPI_Send(&offset, 1, MPI_INT, MASTER, OFFSET, MPI_COMM_WORLD);

        MPI_Send(&rows, 1, MPI_INT, MASTER, ROW, MPI_COMM_WORLD);

        MPI_Send(&local_data[1][0], rows, rowtype, MASTER, RESULT, MPI_COMM_WORLD);
        

    }

    MPI_Finalize();

    return 0;
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

int dowork(double **data, double **local_data, int offset, int rows, int cols,int rank) {
    double max_change, max_black, max_red;
    max_change = 1;
    max_black = 0;
    max_red = 0;
    int i, j;


#pragma omp parallel for shared (data,local_data) private (i, j) reduction(max: max_red)
    for (i = 1; i <= rows; i++) {
        //printf("printed by thread %d\n",omp_get_thread_num());
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
    //printf("max_red: %f\n",max_red);

#pragma omp parallel for shared (data,local_data) private (i, j) reduction(max: max_black)
    for (i = 1; i <= rows; i++) {
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
    printf("rank: %d, max_change: %f\n",rank,max_change);
    if (max_change <= EPSILON) {
//		int msg = 1;
//		MPI_Bsend(&msg, 1, MPI_INT, MASTER, TERMINATION, MPI_COMM_WORLD);
        return 1;
    } else {
        return 0;
    }

}

void write_output(int size, double **result) {
    int i, j;
    FILE *output = fopen("output.dat", "w");
    fprintf(output, "%d\n", size);
    for (i = 0; i < size; i++) {
        for (j = 0; j < size; j++) {
            fprintf(output, "%f ", result[i][j]);
        }
        fputs("\n", output);
    }
    fclose(output);
}