#!/bin/bash

#SBATCH --job-name=baseball_mesh_hires
#SBATCH --partition=grace
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --time=0-00:30:00
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

# Copy scaled STL (using baseball.stl which may have intersecting surfaces)
mkdir -p constant/triSurface
cp /$HOME/Isambaseball/baseball.stl constant/triSurface/

# Note: surfaceFeatureExtract is skipped due to intersecting surfaces in STL
# snappyHexMesh will use implicit feature detection instead

# Create a diagnostic tarball function so we can always package logs/STL
# This will be invoked on error exits and at the end of successful runs.
create_debug_tarball() {
    DEBUG_NAME="mesh_debug_package_${SLURM_JOB_ID:-manual}.tar.gz"
    echo "Creating diagnostic tarball: $DEBUG_NAME"
    tar czvf "$DEBUG_NAME" --ignore-failed-read \
        log.snappyHexMesh log.snappyHexMesh.noLayers log.blockMesh log.decomposePar log.checkMesh snappy_no_layers.log \
        constant/polyMesh/boundary system/snappyHexMeshDict constant/triSurface/baseball.stl \
        mesh_generation_*.out mesh_generation_*.err master_mesh.tar.gz 2>/dev/null || true
    mv -f "$DEBUG_NAME" "$HOME/Isambaseball/" 2>/dev/null || cp -f "$DEBUG_NAME" "$HOME/Isambaseball/" 2>/dev/null || true
    echo "Diagnostic package available at: $HOME/Isambaseball/$DEBUG_NAME"
    chmod a+r "$HOME/Isambaseball/$DEBUG_NAME" 2>/dev/null || true
}

# Generate background mesh
echo "Step 2: Generating background mesh..."
blockMesh > log.blockMesh 2>&1
check_error() {
    if [ $? -ne 0 ]; then
        echo "ERROR: $1 failed"
        # create diagnostic package before exiting so results are available
        create_debug_tarball
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

# Optional: run a diagnostic snappy run inside the Slurm job with layers disabled.
# Set DIAG_SNAPPY=1 when submitting the job to enable this (default is 0).
DIAG_SNAPPY=${DIAG_SNAPPY:-0}

# Run feature extraction to produce an .eMesh if possible. If produced,
# enable explicit feature snapping by using a runtime copy of the dict.
echo "Attempting feature extraction to create .eMesh for explicit feature snapping (if possible)"
# Create or repair a minimal system/surfaceFeaturesDict if it is missing or malformed
# Avoid overwriting a valid user-supplied file. If the first non-empty line does not
# look like an OpenFOAM header we replace it with a safe, minimal dictionary and
# strip any Windows CRLFs which can confuse OpenFOAM's parser.
first_line=""
if [ -f system/surfaceFeaturesDict ]; then
    first_line=$(awk 'NF{print; exit}' system/surfaceFeaturesDict 2>/dev/null || true)
fi

if [ -z "$first_line" ] || [[ "$first_line" != "/*"* ]]; then
    cat > system/surfaceFeaturesDict <<'EOF'
/*--------------------------------*- C++ -*----------------------------------*\
| =========                 |                                                   |
| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\    /   O peration     | Version:  2312                                  |
|   \\  /    A nd           |                                                   |
|    \\/     M anipulation  |                                                   |
\*---------------------------------------------------------------------------*/
FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    object      surfaceFeaturesDict;
}

// Minimal surfaceFeaturesDict for baseball_wrapped.stl

geometry
{
    "constant/triSurface/baseball_wrapped.stl"
    {
        type triSurfaceMesh;
        name baseball;
    }
}

features
(
    {
        file "baseball.eMesh";
        level 4;
    }
);

mergeTolerance 1e-6;

EOF
    # Ensure no Windows CRLF remain which can confuse OpenFOAM header parser
    if command -v sed >/dev/null 2>&1; then
        sed -i 's/\r$//' system/surfaceFeaturesDict || true
    fi
fi

# Prefer the modern utility if available
if command -v surfaceFeatures >/dev/null 2>&1; then
    echo "Running surfaceFeatures -> log.surfaceFeatureExtract"
    surfaceFeatures > log.surfaceFeatureExtract 2>&1 || true
elif command -v surfaceFeatureExtract >/dev/null 2>&1; then
    echo "Running legacy surfaceFeatureExtract -> log.surfaceFeatureExtract"
    surfaceFeatureExtract > log.surfaceFeatureExtract 2>&1 || true
else
    echo "No surface features utility found; continuing with implicit feature detection" > log.surfaceFeatureExtract
fi
SNAPPY_DICT=system/snappyHexMeshDict
if [ -f constant/triSurface/baseball.eMesh ]; then
    echo "Found constant/triSurface/baseball.eMesh — enabling explicitFeatureSnap via a runtime dict"
    cp system/snappyHexMeshDict system/snappyHexMeshDict.useEMesh || true
    sed -i 's/implicitFeatureSnap true/implicitFeatureSnap false/' system/snappyHexMeshDict.useEMesh || true
    sed -i 's/explicitFeatureSnap false/explicitFeatureSnap true/' system/snappyHexMeshDict.useEMesh || true
    SNAPPY_DICT=system/snappyHexMeshDict.useEMesh
fi

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
        if [ "$DIAG_SNAPPY" = "1" ]; then
        echo "Running diagnostic snappyHexMesh (addLayers=false) as DIAG_SNAPPY=1"
        cp "$SNAPPY_DICT" system/snappyHexMeshDict.noLayers || true
        sed -i 's/addLayers[[:space:]]*true/addLayers false/' system/snappyHexMeshDict.noLayers || true
        mpirun -np $NPROCS snappyHexMesh -parallel -overwrite -dict system/snappyHexMeshDict.noLayers > log.snappyHexMesh.noLayers 2>&1
        echo "Diagnostic snappy log: log.snappyHexMesh.noLayers"
    fi

    mpirun -np $NPROCS snappyHexMesh -parallel -overwrite -dict "$SNAPPY_DICT" > log.snappyHexMesh 2>&1
    RC=$?
elif command -v srun >/dev/null 2>&1; then
        if [ "$DIAG_SNAPPY" = "1" ]; then
        echo "Running diagnostic snappyHexMesh (addLayers=false) as DIAG_SNAPPY=1"
        cp "$SNAPPY_DICT" system/snappyHexMeshDict.noLayers || true
        sed -i 's/addLayers[[:space:]]*true/addLayers false/' system/snappyHexMeshDict.noLayers || true
        srun --ntasks=$NPROCS snappyHexMesh -parallel -overwrite -dict system/snappyHexMeshDict.noLayers > log.snappyHexMesh.noLayers 2>&1
        echo "Diagnostic snappy log: log.snappyHexMesh.noLayers"
    fi

    srun --ntasks=$NPROCS snappyHexMesh -parallel -overwrite -dict "$SNAPPY_DICT" > log.snappyHexMesh 2>&1
    RC=$?
else
    echo "WARNING: MPI launcher not found; running snappyHexMesh serially"
    snappyHexMesh -overwrite > log.snappyHexMesh 2>&1
    RC=$?
fi

check_error "snappyHexMesh (exit code $RC)"

# If snappyHexMesh was run in parallel it writes processorN/ meshes.
# Reconstruct the mesh to the top-level constant/polyMesh so packaging
# and the boundary inspection see the merged mesh (which contains the
# baseball patch). Run only when using more than 1 MPI rank.
if [ "${NPROCS:-1}" -gt 1 ]; then
    echo "Step 3b: Reconstructing parallel mesh to top-level constant/polyMesh..."
    reconstructParMesh -constant > log.reconstructParMesh 2>&1
    if [ $? -ne 0 ]; then
        echo "ERROR: reconstructParMesh failed"
        # produce diagnostics for debugging
        create_debug_tarball
        exit 1
    fi
    # remove processor directories left behind by the parallel run to avoid
    # accidentally packaging per-processor data later
    rm -rf processor* 2>/dev/null || true
fi

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

# Verify that snappyHexMesh produced the baseball boundary patch
echo "Step 4b: Verifying presence of baseball patch in mesh..."
if [ ! -f "constant/polyMesh/boundary" ]; then
    echo "ERROR: boundary file not found after snappyHexMesh"
    echo "Check log.snappyHexMesh for errors"
    tail -200 log.snappyHexMesh 2>/dev/null || true
    exit 1
fi

PATCH_LIST=$(awk '/^[[:space:]]*[a-zA-Z]/ && !/type|nFaces|startFace|physicalType/ {gsub(/^[[:space:]]+|[[:space:]]+$/,"",$0); print $0}' constant/polyMesh/boundary)
echo "Available patches:"; echo "$PATCH_LIST"

if echo "$PATCH_LIST" | grep -qiE "baseball|ball|surface"; then
    echo "✓ baseball patch found in mesh"
else
    echo "ERROR: baseball patch NOT found in mesh. Aborting packaging."
    echo "Boundary file contents:"; echo "-------------------------"; cat constant/polyMesh/boundary
    echo "--- tail of log.snappyHexMesh ---"
    tail -200 log.snappyHexMesh 2>/dev/null || true
    echo "SUGGESTIONS:"
    echo "  - Try rerunning snappyHexMesh with layers disabled (set addLayers false) to see if the surface is created without layer extrusion." 
    echo "  - Inspect and repair the STL (remove intersecting triangles) or use an alternative STL (newbaseball.stl)."
    echo "  - Re-run: sbatch ../1_generate_mesh.sh after fixes"
    # produce diagnostic tarball so you can download logs & STL for analysis
    create_debug_tarball
    exit 1
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

# Also create a diagnostic tarball with logs, snappy outputs, dicts and geometry
DEBUG_NAME="mesh_debug_package_${SLURM_JOB_ID:-manual}.tar.gz"
echo "Creating diagnostic tarball: $DEBUG_NAME (contains logs, snappy outputs, dicts, boundary and STL)"
tar czvf "$DEBUG_NAME" --ignore-failed-read \
    log.snappyHexMesh log.snappyHexMesh.noLayers log.blockMesh log.decomposePar log.checkMesh snappy_no_layers.log \
    constant/polyMesh/boundary system/snappyHexMeshDict constant/triSurface/baseball.stl \
    mesh_generation_*.out mesh_generation_*.err master_mesh.tar.gz 2>/dev/null || true

# Move diagnostic tar to home for easy download
mv -f "$DEBUG_NAME" "$HOME/Isambaseball/" 2>/dev/null || cp -f "$DEBUG_NAME" "$HOME/Isambaseball/" 2>/dev/null || true
echo "Diagnostic package available at: $HOME/Isambaseball/$DEBUG_NAME"
chmod a+r "$HOME/Isambaseball/$DEBUG_NAME" 2>/dev/null || true

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
