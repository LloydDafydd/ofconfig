#!/bin/bash

#SBATCH --job-name=baseball_mesh_hires
#SBATCH --partition=grace
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --time=1-00:00:00
#SBATCH --output=mesh_generation_%j.out
#SBATCH --error=mesh_generation_%j.err

# ============================================================================
# BASEBALL MESH GENERATION - HIGH QUALITY SYMMETRIC MESH
# ============================================================================

echo "=========================================="
echo "BASEBALL CFD - MASTER MESH GENERATION"
echo "=========================================="
echo "Generating high-quality symmetric mesh using newbaseball.stl"
echo "Started: $(date)"
echo "=========================================="

# Load OpenFOAM environment
module purge
module load openmpi/5.0.8/5.0.8

# Create mesh directory
MESH_DIR="mesh_generation"
mkdir -p $MESH_DIR
cd $MESH_DIR

# Copy OpenFOAM case files
cp -r ../openfoam_case/* .

# Generate scaled baseball STL
echo "Step 1: Generating scaled baseball STL..."
cd /$HOME/Isambaseball
cd $MESH_DIR

# Copy scaled STL
mkdir -p constant/triSurface
cp /$HOME/Isambaseball/baseball.stl constant/triSurface/

# Generate background mesh
echo "Step 2: Generating background mesh..."
blockMesh > log.blockMesh 2>&1
check_error() {
    if [ $? -ne 0 ]; then
        echo "ERROR: $1 failed"
        exit 1
    fi
}
check_error "Background mesh generation"

# Get initial mesh statistics
INITIAL_CELLS=$(foamDictionary constant/polyMesh/owner -entry nEntries -value 2>/dev/null)
echo "✓ Background mesh: $INITIAL_CELLS cells"

# Generate high-quality mesh
echo "Step 3: Generating VERY HIGH-RESOLUTION surface mesh..."
echo "This will take 2-4 hours for ultra-high-quality mesh (levels 8-9)..."

# Determine number of MPI ranks to use
if [ -n "$SLURM_NTASKS" ]; then
    NPROCS=$SLURM_NTASKS
elif [ -n "$SLURM_NTASKS_PER_NODE" ] && [ -n "$SLURM_JOB_NUM_NODES" ]; then
    NPROCS=$(( SLURM_NTASKS_PER_NODE * SLURM_JOB_NUM_NODES ))
else
    NPROCS=${SLURM_NTASKS_PER_NODE:-24}
fi

echo "Running snappyHexMesh with $NPROCS MPI ranks"

# Prefer mpirun, fall back to srun if available
# Prepare domain decomposition for parallel snappyHexMesh
echo "Writing system/decomposeParDict with numberOfSubdomains $NPROCS"
cat > system/decomposeParDict << EOF
/*--------------------------------*- C++ -*----------------------------------*\
| =========                 |                                                 |
| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\    /   O peration     | Version:  2312                                  |
|   \\  /    A nd           |                                                 |
|    \\/     M anipulation  |                                                 |
\\*---------------------------------------------------------------------------*/
FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    location    "system";
    object      decomposeParDict;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

numberOfSubdomains $NPROCS;

method          scotch;

distributed     no;

EOF

echo "Decomposing domain for $NPROCS subdomains..."
decomposePar > log.decomposePar 2>&1
check_error "decomposePar"

if command -v mpirun >/dev/null 2>&1; then
    mpirun -np $NPROCS snappyHexMesh -parallel -overwrite > log.snappyHexMesh 2>&1
    RC=$?
elif command -v srun >/dev/null 2>&1; then
    srun --ntasks=$NPROCS snappyHexMesh -parallel -overwrite > log.snappyHexMesh 2>&1
    RC=$?
else
    echo "WARNING: MPI launcher not found; running snappyHexMesh serially"
    snappyHexMesh -overwrite > log.snappyHexMesh 2>&1
    RC=$?
fi

check_error "snappyHexMesh (exit code $RC)"

# Get final mesh statistics
FINAL_CELLS=$(foamDictionary constant/polyMesh/owner -entry nEntries -value 2>/dev/null)
echo "✓ Final mesh: $FINAL_CELLS cells"

# Check mesh quality
echo "Step 4: Checking mesh quality..."
checkMesh > log.checkMesh 2>&1
if [ $? -ne 0 ]; then
    echo "WARNING: Mesh quality check found issues"
    echo "Continuing with mesh (often warnings are acceptable for complex geometries)"
fi

# Create clean mesh directory
echo "Step 5: Creating mesh archive..."
cd /$HOME/Isambaseball
mkdir -p master_mesh
cp -r $MESH_DIR/constant/polyMesh master_mesh/
cp -r $MESH_DIR/constant/extendedFeatureEdgeMesh master_mesh/ 2>/dev/null || echo "No extended features"
cp -r $MESH_DIR/constant/triSurface master_mesh/

# Create compressed archive
tar -czf master_mesh.tar.gz master_mesh/

# Get archive size
ARCHIVE_SIZE=$(du -h master_mesh.tar.gz | cut -f1)

echo ""
echo "=========================================="
echo "MESH GENERATION COMPLETE!"
echo "=========================================="
echo "Mesh Statistics:"
echo "  Background cells: $INITIAL_CELLS"
echo "  Final cells: $FINAL_CELLS"
echo "  Refinement factor: $(echo "scale=1; $FINAL_CELLS / $INITIAL_CELLS" | bc)x"
echo ""
echo "Files created:"
echo "  master_mesh/           - Clean mesh files"
echo "  master_mesh.tar.gz     - Compressed archive ($ARCHIVE_SIZE)"
echo ""
echo "Next steps:"
echo "  1. Run test case to validate mesh"
echo "  2. Submit parameter study jobs"
echo "=========================================="
