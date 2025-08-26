#!/bin/bash

#SBATCH --job-name=foamVTK
#SBATCH --partition=grace
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
# MPI ranks per node (hybrid: ranks * threads = total cores)
#SBATCH --cpus-per-task=1
# OpenMP threads per rank
#SBATCH --mem=16G
#SBATCH --time=1-00:00:00
#SBATCH --output=test_case_%j.out
#SBATCH --error=test_case_%j.err

# Load required modules
module purge
module load OpenFOAM/v2312-foss-2023a  # Adjust version as needed

# Set OpenFOAM environment
source $FOAM_BASH

# Change to case directory
cd /Isambaseball/jobs/results/v085_s2200_m015_a030

# Print job information
echo "Job started at: $(date)"
echo "Running on node: $(hostname)"
echo "Working directory: $(pwd)"
echo "Number of MPI tasks: $SLURM_NTASKS"

# Check if case directory exists and contains necessary files
if [ ! -f "system/controlDict" ]; then
    echo "Error: controlDict not found in system/ directory"
    exit 1
fi

# Run foamToVTK
echo "Starting foamToVTK conversion..."
foamToVTK -parallel

# Alternative: If you want to run without parallel flag
# foamToVTK

# Check if conversion was successful
if [ $? -eq 0 ]; then
    echo "foamToVTK completed successfully at: $(date)"
    
    # Optional: Print VTK directory contents
    if [ -d "VTK" ]; then
        echo "VTK files created:"
        ls -la VTK/
    fi
else
    echo "foamToVTK failed with exit code: $?"
    exit 1
fi

echo "Job completed at: $(date)"
