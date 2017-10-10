#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>

#define EPSILON 5e-2


void write_output(int num_lines, double **result) {
    int i, j;
    FILE *output = fopen("output2.dat", "w");
    fprintf(output, "%d\n", num_lines);
    for (i = 0; i < num_lines; i++) {
        for (j = 0; j < num_lines; j++) {
            fprintf(output, "%f ", result[i][j]);
        }
        fputs("\n", output);
    }
    fclose(output);
}

int main(int argc, char **argv) {
    char str[10];
    int num_lines, i, j;
    double start, end, time, diff, my_diff;
    start = omp_get_wtime();
    FILE *input = fopen("input.dat", "r");
    if (!input) {
        printf("File read error");
        exit(EXIT_FAILURE);
    }
    fgets(str, sizeof(str), input);
    num_lines = atoi(str);

    // Initialize the array to store the original data
    double **data;
    data = malloc(sizeof(double *) * num_lines);

    // Read the file
    for (i = 0; i < num_lines; i++) {
        data[i] = malloc(sizeof(double) * num_lines);
        for (j = 0; j < num_lines; j++) {
            fscanf(input, "%lf ", &(data[i][j]));
        }
    }
    fclose(input);


    // Initialize the array to store the final result
    double **result;
    result = malloc(sizeof(double *) * num_lines);
    for (i = 0; i < num_lines; i++) {
        result[i] = malloc(sizeof(double) * num_lines);
        for (j = 0; j < num_lines; j++) {
            result[i][j] = data[i][j];
        }
    }

    double **prev;
    prev = malloc(sizeof(double *) * num_lines);
    for (i = 0; i < num_lines; i++) {
        prev[i] = malloc(sizeof(double) * num_lines);
    }

    double max_change = 1;
    int num_iter = 0;
//    omp_set_num_threads(2);
    // Stop when the change of all points is smaller than the EPSILON
    while (max_change > EPSILON) {
        num_iter++;
        max_change = 0;

#pragma omp parallel shared (result, prev) private ( i, j )
        {
#pragma omp for
            for (i = 0; i < num_lines; i++) {
                for (j = 0; j < num_lines; j++) {
                    prev[i][j] = result[i][j];
                }
            }
#pragma omp for
            for (i = 0; i < num_lines; i++) {
                for (j = 0; j < num_lines; j++) {
                    // Only calculate the point without an initial value
                    if (data[i][j] == 0) {
                        // Set the point next to the current point to be 0 if the point is on boundary
                        double left, right, up, down, diff;
                        left = i - 1 < 0 ? 0 : prev[i - 1][j];
                        right = i + 1 >= num_lines ? 0 : prev[i + 1][j];
                        up = j - 1 < 0 ? 0 : prev[i][j - 1];
                        down = j + 1 >= num_lines ? 0 : prev[i][j + 1];

                        result[i][j] = 0.25 * (left + right + up + down);
                    }
                }
            }
        }


# pragma omp parallel shared ( diff, result, prev ) private ( i, j, my_diff )
        {
#pragma omp for
            for ( i = 0; i < num_lines; i++ )
            {
                for ( j = 0; j < num_lines; j++ )
                {
                    if (my_diff < fabs(prev[i][j] - result[i][j]))
                    {
                        my_diff = fabs ( prev[i][j] - result[i][j] );
                    }
                }
            }
# pragma omp critical
            {
                if ( max_change < my_diff )
                {
                    max_change = my_diff;
                }
            }
        }

        printf("%d %f\n", num_iter, max_change);
    }
    write_output(num_lines, result);
    end = omp_get_wtime();
    time = end - start;
    printf("%lf", time);
}