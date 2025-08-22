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

# Update magUInf in force coefficients to match velocity
sed -i "s/magUInf.*38.0;/magUInf         $VELOCITY;/" system/controlDict

# Create omega field before decomposition
echo "Creating omega field..."
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

# Check if baseball patch exists in the mesh
echo "Checking mesh patches..."

# First check if boundary file exists
if [ ! -f "constant/polyMesh/boundary" ]; then
    echo "ERROR: Mesh boundary file not found! Check mesh installation."
    exit 1
fi

# Extract patch names from boundary file
ACTUAL_PATCHES=$(awk '/^[[:space:]]*[a-zA-Z]/ && !/type|nFaces|startFace|physicalType/ {gsub(/^[[:space:]]+|[[:space:]]+$/,"",$0); print $0}' constant/polyMesh/boundary)
echo "Available patches:"
echo "$ACTUAL_PATCHES"

# Find the baseball patch (or closest match)
BASEBALL_PATCH=""
for patch in $ACTUAL_PATCHES; do
    if echo "$patch" | grep -q -E "(baseball|ball|wall|surface)"; then
        BASEBALL_PATCH="$patch"
        echo "Found potential baseball patch: $patch"
        break
    fi
done

if [ -z "$BASEBALL_PATCH" ]; then
    echo "ERROR: No baseball/wall patch found in mesh!"
    echo "Available patches: $ACTUAL_PATCHES"
    echo ""
    echo "Mesh boundary file contents:"
    cat constant/polyMesh/boundary
    exit 1
fi

echo "Using baseball patch: $BASEBALL_PATCH"

# Update controlDict with correct patch name
sed -i "s/patches.*baseball.*/patches         ($BASEBALL_PATCH);/" system/controlDict

# Domain decomposition for parallel execution
echo "Decomposing domain for $SLURM_NTASKS cores..."
decomposePar > log.decomposePar 2>&1

if [ $? -ne 0 ]; then
    echo "ERROR: Domain decomposition failed"
    cat log.decomposePar
    exit 1
fi

# Distribute omega field to processor directories
echo "Distributing omega field to processor directories..."
redistributePar -overwrite > log.redistributePar 2>&1

if [ $? -ne 0 ]; then
    echo "ERROR: Field redistribution failed"
    tail -50 log.redistributePar
    exit 1
fi

# Run simulation using modern OpenFOAM syntax
echo "Starting CFD simulation..."
echo "This will take 2-4 hours..."

# Use srun if available (Slurm native), otherwise mpirun
if command -v srun > /dev/null 2>&1; then
    echo "Launching with srun (Slurm) on $SLURM_NTASKS tasks..."
    srun --ntasks=$SLURM_NTASKS foamRun -solver incompressibleFluid -parallel > log.foamRun 2>&1
else
    echo "Launching with mpirun on $SLURM_NTASKS ranks..."
    mpirun -np $SLURM_NTASKS foamRun -solver incompressibleFluid -parallel > log.foamRun 2>&1
fi

if [ $? -ne 0 ]; then
    echo "ERROR: Simulation failed"
    tail -50 log.foamRun
    exit 1
fi

# Reconstruct parallel results
echo "Reconstructing results..."
reconstructPar > log.reconstructPar 2>&1

if [ $? -ne 0 ]; then
    echo "ERROR: Reconstruction failed"
    tail -20 log.reconstructPar
    exit 1
fi

# Post-process force data
echo "Post-processing force data..."

# Check if force data exists
if [ ! -f "postProcessing/forces/0/force.dat" ]; then
    echo "ERROR: Force data not found. Running postProcess..."
    
    # Try to generate force data using postProcess
    foamPostProcess -func forces -latestTime > log.postProcess 2>&1
    foamPostProcess -func forceCoeffs -latestTime >> log.postProcess 2>&1
fi

# Extract and calculate force coefficients
if [ -f "postProcessing/forces/0/force.dat" ]; then
    echo "Processing force data from force.dat..."
    
    # Constants for coefficient calculation
    RHO=1.225
    U=$VELOCITY
    LREF=0.0732
    AREF=0.004218
    
    # Calculate coefficients from force.dat
    # Columns: time Fx Fy Fz Mx My Mz ...
    awk -v rho="$RHO" -v U="$U" -v A="$AREF" -v L="$LREF" '
        BEGIN {
            q = 0.5 * rho * U * U
            print "# time Cd Cl Cm"
        }
        NF >= 7 && $1 > 0.1 {
            Cd = $2 / (q * A)        # Drag from Fx
            Cl = $4 / (q * A)        # Lift from Fz
            Cm = $6 / (q * A * L)    # Moment from My
            printf "%.6f %.8f %.8f %.8f\n", $1, Cd, Cl, Cm
        }
    ' postProcessing/forces/0/force.dat > forces_coeffs.dat
    
    # Calculate statistics
    awk '
        NR > 1 {
            n++
            sum_cd += $2; sum_cl += $3; sum_cm += $4
            sum_cd2 += $2*$2; sum_cl2 += $3*$3; sum_cm2 += $4*$4
        }
        END {
            if (n > 0) {
                cd_avg = sum_cd / n
                cl_avg = sum_cl / n  
                cm_avg = sum_cm / n
                
                if (n > 1) {
                    cd_std = sqrt((sum_cd2 - n*cd_avg*cd_avg)/(n-1))
                    cl_std = sqrt((sum_cl2 - n*cl_avg*cl_avg)/(n-1))
                    cm_std = sqrt((sum_cm2 - n*cm_avg*cm_avg)/(n-1))
                } else {
                    cd_std = 0; cl_std = 0; cm_std = 0
                }
                
                printf "CD_AVG=%.6f\n", cd_avg
                printf "CL_AVG=%.6f\n", cl_avg
                printf "CM_AVG=%.6f\n", cm_avg
                printf "CD_STD=%.6f\n", cd_std
                printf "CL_STD=%.6f\n", cl_std
                printf "CM_STD=%.6f\n", cm_std
            }
        }
    ' forces_coeffs.dat > force_stats.txt
    
    # Load statistics
    source force_stats.txt
    
elif [ -f "postProcessing/forceCoeffs/0/coefficient.dat" ]; then
    echo "Processing force data from coefficient.dat..."
    
    # Extract coefficients directly
    awk '$1 > 0.1 {print $1, $2, $3, $4}' postProcessing/forceCoeffs/0/coefficient.dat > forces_coeffs.dat
    
    # Calculate statistics  
    awk '
        {
            n++
            sum_cd += $2; sum_cl += $3; sum_cm += $4
            sum_cd2 += $2*$2; sum_cl2 += $3*$3; sum_cm2 += $4*$4
        }
        END {
            if (n > 0) {
                cd_avg = sum_cd / n
                cl_avg = sum_cl / n
                cm_avg = sum_cm / n
                
                if (n > 1) {
                    cd_std = sqrt((sum_cd2 - n*cd_avg*cd_avg)/(n-1))
                    cl_std = sqrt((sum_cl2 - n*cl_avg*cl_avg)/(n-1))
                    cm_std = sqrt((sum_cm2 - n*cm_avg*cm_avg)/(n-1))
                } else {
                    cd_std = 0; cl_std = 0; cm_std = 0
                }
                
                printf "CD_AVG=%.6f\n", cd_avg
                printf "CL_AVG=%.6f\n", cl_avg
                printf "CM_AVG=%.6f\n", cm_avg
                printf "CD_STD=%.6f\n", cd_std
                printf "CL_STD=%.6f\n", cl_std
                printf "CM_STD=%.6f\n", cm_std
            }
        }
    ' forces_coeffs.dat > force_stats.txt
    
    # Load statistics
    source force_stats.txt
    
else
    echo "ERROR: No force data found in postProcessing directory"
    echo "Available postProcessing contents:"
    find postProcessing -name "*.dat" 2>/dev/null || echo "No .dat files found"
    
    CD_AVG="N/A"; CL_AVG="N/A"; CM_AVG="N/A"
    CD_STD="N/A"; CL_STD="N/A"; CM_STD="N/A"
fi

# Display results
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
if [ "$CL_AVG" != "N/A" ]; then
    CL_ABS=$(echo "$CL_AVG" | awk '{print ($1<0)?-$1:$1}')
    if awk "BEGIN {exit !($CL_ABS < 0.05)}"; then
        echo "✓ VALIDATION PASSED: Lateral force is acceptably small"
        echo "  |Cl| = $CL_ABS < 0.05 (expected for symmetric case)"
    else
        echo "⚠ VALIDATION WARNING: Lateral force higher than expected"
        echo "  |Cl| = $CL_ABS >= 0.05 (check mesh quality and setup)"
    fi
fi

echo ""
echo "=========================================="
echo "TEST CASE COMPLETE"
echo "=========================================="
echo "Completed: $(date)"
echo "Check results in: $WORK_DIR"
echo "Force coefficients saved to: forces_coeffs.dat"
echo "Next step: Run parameter study if validation passed"
echo "========================================"
