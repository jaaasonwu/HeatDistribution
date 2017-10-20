#!/bin/bash

mpicc -o mpi_heat_distribution mpi_heat_distribution.c -fopenmp -O3
