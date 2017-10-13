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
#define DONE        4                  /* message tag */
#define TERMINATION	5					/* message tag */
#define MAXITER     10              /* maximum iterations to go through */
#define EPSILON 	5e-2


double **alloc_2d_double(int rows, int cols);
void dowork(double **data,double **local_data, int offset, int rows, int cols);
void print_2d_array(double **A,int rows,int cols);
void write_output(int size, double **result);

int main(int argc, char *argv[])
{
	int world_size, rank,
		numworkers,					/* number of worker processes */
		averow,rows,extra,
	 	offset,						/* for sending rows of data */
	 	dest, source,               /* to - from for message send-receive */
	 	left,right,       			/* neighbor tasks */
		msgtype,                  /* for message types */
		msg;

	MPI_Datatype rowtype;
	int size = 0;
	int i, j;
	double **data;

	// Initializes the MPI execution environment
	int provided;
	//printf("BEGIN\n");
	MPI_Init_thread(&argc,&argv,MPI_THREAD_MULTIPLE,&provided);
	//printf("INIT_THEAD_DONE\n");
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &world_size);
	numworkers = world_size-1;
	MPI_Status Stat;

	char str[10];
	FILE *input = fopen("input.dat", "r");
	if (!input) {
		printf("File read error");
		exit(EXIT_FAILURE);
	}
	fgets(str, sizeof(str), input);
	size = atoi(str);

    // Initialize the array to store the original data				
	data = alloc_2d_double(size,size);

    // Read the file
	for (i = 0; i < size; i++) {			
		for (j = 0; j < size; j++) {
			fscanf(input, "%lf ", &(data[i][j]));
		}
	}
	fclose(input);

	averow = size/numworkers;
	extra = size%numworkers;
	MPI_Type_vector(1, size, 1, MPI_DOUBLE, &rowtype);
	MPI_Type_commit(&rowtype);

	//MPI_Bcast(&size,1, MPI_INT ,0, MPI_COMM_WORLD);
	//printf("rank: %d, size: %d\n",rank,size);


	/************************* master code *******************************/
	if(rank == MASTER){
		printf ("Starting mpi_heat2D with %d worker tasks.\n", numworkers);

		#pragma omp parallel private(rows,offset,source,msgtype,msg,Stat)
		{
			#pragma omp for
			for (i = 0; i<=numworkers; i++)
			{
				source = i;
				msgtype = TERMINATION;
				MPI_Recv(&msg, 1, MPI_INT, source, msgtype, MPI_COMM_WORLD, &Stat);
			}

			/* wait for results from all worker tasks */
			#pragma omp for 
			for (i=1; i<=numworkers; i++)
			{
				source = i;
				msgtype = DONE;
				// MPI_Recv(&offset, 1, MPI_INT, source, msgtype, MPI_COMM_WORLD, &Stat);
				// MPI_Recv(&rows, 1, MPI_INT, source, msgtype, MPI_COMM_WORLD, &Stat);
				// MPI_Recv(&data[offset][0], rows, rowtype, source,
				// 	msgtype, MPI_COMM_WORLD, &Stat);
			}
		}
		/* All threads join master thread */
		/* Write final output */
		//write_output(size, data);

	}
	/************************* workers code **********************************/
	if (rank != MASTER)
	{
		int converged = 0;
		rows = (rank <= extra) ? averow+1 : averow; 
		if(rank<=extra){
			rows = averow+1;
			offset = (rank-1)*rows;
			//end_row = start_row+rows-1;
		}else{
			rows = averow;
			offset = extra*(averow+1)+(rank-extra)*averow - rows;
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

		/* Allocate memory for local data */
		double **local_data;
		/* my partition, plus two neighbouring rows */
		local_data = alloc_2d_double(rows+2,size);
		/* initialise the two neighbouring rows to zero */
		for(i=0;i<rows+2;i++){
			if(i==0 || i==rows+1)
			{
				for(j=0;j<size;j++){
					local_data[i][j]=0;
				}		
			}else
			{
				for(j=0;j<size;j++){
					local_data[i][j]=data[offset+i-1][j];
				}		
			}
		}

		printf("rank: %d, offset: %d,rows:%d\n",rank,offset,rows);
		// if(rank==1){
		// 	print_2d_array(local_data,rows+2,size);
		// }
		
		// if (left != NONE)
		// {
		// 		// send the first row of the worker's partition to the left neighbour
		// 	MPI_Send(&local_data[1][0], 1,rowtype, left, RTAG, MPI_COMM_WORLD);
		// 	source = left;
		// 	msgtype = LTAG;
		// 		// receive from left neighbour
		// 	MPI_Recv(&local_data[0][0], 1, rowtype, source,msgtype, MPI_COMM_WORLD, &Stat);
		// }
		// if (right != NONE)
		// {
		// 	MPI_Send(&local_data[rows][0], 1,rowtype, right, LTAG, MPI_COMM_WORLD);
		// 	source = right;
		// 	msgtype = RTAG;
		// 	MPI_Recv(&local_data[rows+1][0], 1, rowtype, source,msgtype, MPI_COMM_WORLD, &Stat);
		// }

		// if (rank==1)
		// {
		// 	print_2d_array(local_data,rows+2,size);
		// }


		int num_iter = 0;
		
		//while(!converged)
		while(num_iter<=10)
		{
			num_iter++;
			printf("iter: %d\n", num_iter);
			if (left != NONE)
			{
				// send the first row of the worker's partition to the left neighbour
				MPI_Send(&local_data[1][0], 1,rowtype, left, RTAG, MPI_COMM_WORLD);
				source = left;
				msgtype = LTAG;
				// receive from left neighbour
				MPI_Recv(&local_data[0][0], 1, rowtype, source,msgtype, MPI_COMM_WORLD, &Stat);
			}
			if (right != NONE)
			{
				MPI_Send(&local_data[rows][0], 1,rowtype, right, LTAG, MPI_COMM_WORLD);
				source = right;
				msgtype = RTAG;
				MPI_Recv(&local_data[rows+1][0], 1, rowtype, source,msgtype, MPI_COMM_WORLD, &Stat);
			}
			//only update the first to the second-last row in the local data
			dowork(data,local_data,offset,rows,size);
		}
		
		 /* Finally, send my portion of final results back to master */
		// MPI_Send(&offset, 1, MPI_INT, MASTER, DONE, MPI_COMM_WORLD);
		// MPI_Send(&rows, 1, MPI_INT, MASTER, DONE, MPI_COMM_WORLD);
		// MPI_Send(&local_data[1][0], rows, rowtype, MASTER, DONE, MPI_COMM_WORLD);

	}

	MPI_Finalize();
	
	return 0;
}

double **alloc_2d_double(int rows, int cols) {
	double *data = (double *)malloc(rows*cols*sizeof(double));
	double **array= (double **)malloc(rows*sizeof(double*));
	int i;
	for (i=0; i<rows; i++)
		array[i] = &(data[cols*i]);

	return array;
}

void print_2d_array(double **A,int rows,int cols)
{
	int i,j;
	for (i = 0; i < rows; i++)
	{
		for(j=0;j<cols;j++)
		{
			printf("%f ",A[i][j]);
		}
		printf("\n");
	}
}

void dowork(double **data,double **local_data, int offset, int rows, int cols)
{
	double max_change, max_black, max_red;
	max_change = 1;
	max_black = 0;
	max_red = 0;
	int i,j;

	
	#pragma omp parallel for shared (data,local_data) private (i, j) reduction(max: max_red)
	for (i = 1; i <= rows; i++) {
		//printf("printed by thread %d\n",omp_get_thread_num());
		j = i % 2 == 0 ? 0 : 1;
		for (; j < cols; j += 2) {
            // Only calculate the point without an initial value
			if (data[offset+i-1][j] == 0) {
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
			if (data[offset+i-1][j] == 0) {
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
	printf("max_change: %f\n",max_change);
	if(max_change<=EPSILON){
		int msg = 1;
		MPI_Send(&msg, 1, MPI_INT, MASTER, TERMINATION, MPI_COMM_WORLD);
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