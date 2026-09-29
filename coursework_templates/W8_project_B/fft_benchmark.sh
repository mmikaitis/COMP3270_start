#!/bin/bash
#SBATCH --partition=teachingnodes
#SBATCH --account=2026COMP3270
#SBATCH --job-name=fft_mpi_benchmark
#SBATCH --output=fft_%j.out
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=32
#SBATCH --cpus-per-task=1
#SBATCH --mem=35G
#SBATCH --time=00:15:00
set -euo pipefail
module load openmpi
export OMP_NUM_THREADS=1

RESULTS=fft_mpi_benchmark_results.txt
echo "procs repeat runtime error" > "$RESULTS"
for p in 8 16 32 64; do
  mpirun -n "$p" --map-by "ppr:$((p / 2)):node" --bind-to core \
    ./fft_mpi 16384 3 >> "$RESULTS"
done
