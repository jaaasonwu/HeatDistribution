#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>

#define EPSILON 5e-2


void write_output(int num_lines, double **result) {
    int i, j;
    FILE *output = fopen("omp_output.dat", "w");
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
    double start, end, time;
    start = omp_get_wtime();
    if (argc == 1) {
        printf("Missing input file");
        exit(EXIT_FAILURE);
    }
    FILE *input = fopen(argv[1], "r");
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


    double max_change = 1, max_black, max_red;
    int num_iter = 0;
    // Stop when the change of all points is smaller than the EPSILON

    while (max_change > EPSILON) {
        num_iter++;
        max_change = 0;
        max_black = 0;
        max_red = 0;
#pragma omp parallel for shared (result) private (i, j) reduction(max: max_red)
        // Calculate all the red dots
        for (i = 0; i < num_lines; i++) {
            j = i % 2 == 0 ? 0 : 1;
            for (; j < num_lines; j += 2) {
                // Only calculate the point without an initial value
                if (data[i][j] == 0) {
                    // Set the point next to the current point to be 0 if the point is on boundary
                    double left, right, up, down, diff, prev;
                    left = i - 1 < 0 ? 0 : result[i - 1][j];
                    right = i + 1 >= num_lines ? 0 : result[i + 1][j];
                    up = j - 1 < 0 ? 0 : result[i][j - 1];
                    down = j + 1 >= num_lines ? 0 : result[i][j + 1];

                    prev = result[i][j];
                    result[i][j] = 0.25 * (left + right + up + down);

                    diff = fabs(prev - result[i][j]);
                    if (diff > max_red) {
                        max_red = diff;
                    }
                }
            }
        }

#pragma omp parallel for shared(result) private(i, j) reduction(max:max_black)
        // calculate all the black dots
        for (i = 0; i < num_lines; i++) {
            j = i % 2 == 0 ? 1 : 0;
            for (; j < num_lines; j += 2) {
                // Only calculate the point without an initial value
                if (data[i][j] == 0) {
                    // Set the point next to the current point to be 0 if the point is on boundary
                    double left, right, up, down, diff, prev;
                    left = i - 1 < 0 ? 0 : result[i - 1][j];
                    right = i + 1 >= num_lines ? 0 : result[i + 1][j];
                    up = j - 1 < 0 ? 0 : result[i][j - 1];
                    down = j + 1 >= num_lines ? 0 : result[i][j + 1];

                    prev = result[i][j];
                    result[i][j] = 0.25 * (left + right + up + down);

                    diff = fabs(prev - result[i][j]);
                    if (diff > max_black) {
                        max_black = diff;
                    }
                }
            }
        }

        max_change = max_black > max_red ? max_black : max_red;
//        printf("%d %f\n", num_iter, max_change);
    }
    write_output(num_lines, result);
    end = omp_get_wtime();
    time = end - start;
    printf("%lf", time);
}