#!/bin/bash

#SBATCH --job-name=baseball_test
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --time=1-00:00:00
#SBATCH --partition=grace
#SBATCH --output=test_case_%j.out
#SBATCH --error=test_case_%j.err

# ============================================================================
# BASEBALL TEST CASE - VALIDATION SIMULATION
# ============================================================================
# Test parameters: 85 mph, 2000 RPM, 0° seam, 0° spin angle
# Should produce near-zero average horizontal force for validation
# ============================================================================

echo "=========================================="
echo "BASEBALL CFD - TEST CASE SIMULATION"
echo "=========================================="
echo "Test conditions:"
echo "  Velocity: 85 mph (38.0 m/s)"
echo "  Spin rate: 2000 RPM"
echo "  Seam orientation: 0°"
echo "  Spin angle: 0°"
echo "  Expected: Near-zero lateral force"
echo "Started: $(date)"
echo "=========================================="

# Load OpenFOAM environment
module purge
module load openmpi/5.0.8/5.0.8

# Create working directory
WORK_DIR="test_case_v85_s2000"
mkdir -p $WORK_DIR
cd $WORK_DIR

# Copy base case
cp -r /$HOME/Isambaseball/openfoam_case/* .

# Install pre-generated mesh
if [ -f "/$HOME/Isambaseball/master_mesh.tar.gz" ]; then
    echo "Installing pre-generated mesh..."
    tar -xzf /$HOME/Isambaseball/master_mesh.tar.gz
    cp -r master_mesh/* constant/
    rm -rf master_mesh
    echo "✓ Mesh installed"
else
    echo "ERROR: Master mesh not found! Run mesh generation first."
    exit 1
fi

# Configure test parameters
VELOCITY=38.0      # 85 mph in m/s
SPIN_RATE=2000     # RPM
SEAM_ANGLE=0       # degrees
SPIN_AXIS_ANGLE=0  # degrees (pure backspin/topspin)

# Calculate spin rate in rad/s
OMEGA=$(echo "scale=6; $SPIN_RATE * 2 * 3.14159 / 60" | bc)

echo "Configuring boundary conditions..."
echo "  Velocity: $VELOCITY m/s"
echo "  Angular velocity: $OMEGA rad/s"

# Update velocity in boundary conditions
sed -i "s/uniform (38.0 0 0)/uniform ($VELOCITY 0 0)/" 0/U

# Update rotational velocity for spinning baseball
cat > 0/omega << EOF
/*--------------------------------*- C++ -*----------------------------------*\
FoamFile
{
    format      ascii;
    class       volVectorField;
    object      omega;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

dimensions      [0 0 -1 0 0 0 0];

internalField   uniform (0 $OMEGA 0);  // Pure backspin

boundaryField
{
    baseball
    {
        type            fixedValue;
        value           uniform (0 $OMEGA 0);
    }
    
    inlet
    {
        type            zeroGradient;
    }
    
    outlet  
    {
        type            zeroGradient;
    }
    
    sides
    {
        type            zeroGradient;
    }
    
    top
    {
        type            zeroGradient;
    }
    
    bottom
    {
        type            zeroGradient;
    }
}
EOF

# Domain decomposition for parallel execution
echo "Decomposing domain for $SLURM_NTASKS cores..."
decomposePar > log.decomposePar 2>&1

if [ $? -ne 0 ]; then
    echo "ERROR: Domain decomposition failed"
    exit 1
fi

# Verify decomposition created processor directories and that each has the 0/ files
echo "Verifying decomposition and presence of 0/ files in processor directories..."
NUM_PROCS=${SLURM_NTASKS:-1}
MISSING_PROC=0
for ((i=0;i<NUM_PROCS;i++)); do
    PDIR="processor${i}"
    if [ ! -d "$PDIR/0" ]; then
        echo "ERROR: Expected directory $PDIR/0 does not exist"
        MISSING_PROC=1
    else
        # If specific runtime fields (e.g. omega) are missing, copy the main 0/ files as a fallback
        if [ ! -f "$PDIR/0/omega" ]; then
            echo "Notice: $PDIR/0/omega missing — copying 0/* -> $PDIR/0/"
            cp -r 0/* "$PDIR/0/" || {
                echo "ERROR: Failed to copy 0/* to $PDIR/0/"
                exit 1
            }
        fi
    fi
done

if [ "$MISSING_PROC" -ne 0 ]; then
    echo "ERROR: One or more processor directories are missing; check log.decomposePar"
    ls -l processor* || true
    tail -n 80 log.decomposePar || true
    exit 1
fi

# Run simulation
echo "Starting CFD simulation..."
echo "This will take 2-4 hours..."

mpirun -np $SLURM_NTASKS pimpleFoam -parallel > log.pimpleFoam 2>&1

if [ $? -ne 0 ]; then
    echo "ERROR: Simulation failed"
    tail -50 log.pimpleFoam
    exit 1
fi

# Reconstruct parallel results
echo "Reconstructing results..."
reconstructPar > log.reconstructPar 2>&1

# Extract force coefficients
echo "Extracting force data..."
grep -E "^[0-9]" postProcessing/forces/0/forceCoeffs.dat > forces_clean.dat 2>/dev/null || echo "Warning: Force coefficients not found"

# Calculate average forces
if [ -f "forces_clean.dat" ]; then
    echo "Calculating force statistics..."
    
    # Skip initial transient (first 0.1 seconds)
    awk '$1 > 0.1 {print $2}' forces_clean.dat > cd_data.tmp
    awk '$1 > 0.1 {print $3}' forces_clean.dat > cl_data.tmp
    awk '$1 > 0.1 {print $4}' forces_clean.dat > cm_data.tmp
    
    # Calculate averages using awk
    CD_AVG=$(awk '{sum+=$1; n++} END {if(n>0) print sum/n; else print "N/A"}' cd_data.tmp)
    CL_AVG=$(awk '{sum+=$1; n++} END {if(n>0) print sum/n; else print "N/A"}' cl_data.tmp)
    CM_AVG=$(awk '{sum+=$1; n++} END {if(n>0) print sum/n; else print "N/A"}' cm_data.tmp)
    
    # Calculate standard deviations
    CD_STD=$(awk -v avg=$CD_AVG '{sum+=($1-avg)^2; n++} END {if(n>1) print sqrt(sum/(n-1)); else print "N/A"}' cd_data.tmp)
    CL_STD=$(awk -v avg=$CL_AVG '{sum+=($1-avg)^2; n++} END {if(n>1) print sqrt(sum/(n-1)); else print "N/A"}' cl_data.tmp)
    CM_STD=$(awk -v avg=$CM_AVG '{sum+=($1-avg)^2; n++} END {if(n>1) print sqrt(sum/(n-1)); else print "N/A"}' cm_data.tmp)
    
    rm -f *.tmp
    
    echo ""
    echo "=========================================="
    echo "FORCE COEFFICIENT ANALYSIS"
    echo "=========================================="
    echo "Average force coefficients (excluding 0.1s transient):"
    echo "  Drag coefficient (Cd):    $CD_AVG ± $CD_STD"
    echo "  Lift coefficient (Cl):    $CL_AVG ± $CL_STD"  
    echo "  Moment coefficient (Cm):  $CM_AVG ± $CM_STD"
    echo ""
    
    # Validation check
    CL_ABS=$(echo "$CL_AVG" | awk '{print ($1<0)?-$1:$1}')
    if awk "BEGIN {exit !($CL_ABS < 0.05)}"; then
        echo "✓ VALIDATION PASSED: Lateral force is acceptably small"
        echo "  |Cl| = $CL_ABS < 0.05 (expected for symmetric case)"
    else
        echo "⚠ VALIDATION WARNING: Lateral force higher than expected"
        echo "  |Cl| = $CL_ABS >= 0.05 (check mesh quality and setup)"
    fi
else
    echo "ERROR: Force data not found. Check simulation logs."
fi

echo ""
echo "=========================================="
echo "TEST CASE COMPLETE"
echo "=========================================="
echo "Completed: $(date)"
echo "Check results in: $WORK_DIR"
echo "Next step: Run parameter study if validation passed"
echo "=========================================="
