Run the default configuration: (size 1000, epsilon 0.001):
Compile using: bash compile.sh
Run using: sbatch run.sh

Change the size of input or epsilon:
Go to run.sh. The line that runs the program is in the format of:
srun mpi_heat_distribution {size} {epsilon} {output_filename}
Just change the corresponding fields.

Visualize the output:
On a computer with python, numpy and matplotlib installed:
python visualize.py {filename}