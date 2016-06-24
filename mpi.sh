export OMP_NUM_THREADS=40; mpirun -n 1 ./mpi
export OMP_NUM_THREADS=10; mpirun -n 4 ./mpi
export OMP_NUM_THREADS=2; mpirun -n 16 ./mpi
