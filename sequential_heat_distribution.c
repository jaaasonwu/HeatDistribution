#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define EPSILON 1e-3


void write_output(int num_lines, double **result) {
    int i, j;
    FILE *output = fopen("output.dat", "w");
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

    double max_change = 1;
    int num_iter = 0;
    // Stop when the change of all points is smaller than the EPSILON
    while (max_change > EPSILON) {
        num_iter++;
        max_change = 0;

        for (i = 0; i < num_lines; i++) {
            for (j = 0; j < num_lines; j++) {
                // Only calculate the point without an initial value
                if (data[i][j] == 0) {
                    // Set the point next to the current point to be 0 if the point is on boundary
                    double left, right, up, down, prev, diff;
                    left = i - 1 < 0 ? 0 : result[i - 1][j];
                    right = i + 1 >= num_lines ? 0 : result[i + 1][j];
                    up = j - 1 < 0 ? 0 : result[i][j - 1];
                    down = j + 1 >= num_lines ? 0 : result[i][j + 1];

                    prev = result[i][j];
                    result[i][j] = 0.25 * (left + right + up + down);
                    // Get the changes between two iterations and update the max_change
                    diff = fabs(prev - result[i][j]);
                    max_change = diff > max_change ? diff : max_change;
                }
            }
        }
        printf("%d %f\n", num_iter, max_change);
    }
    write_output(num_lines, result);
}