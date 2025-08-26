#!/bin/bash

#SBATCH --job-name=baseball_test
#SBATCH --partition=grace
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24   # MPI ranks per node (hybrid: ranks * threads = total cores)
#SBATCH --cpus-per-task=4     # OpenMP threads per rank
#SBATCH --mem=128G
#SBATCH --time=1-00:00:00
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

# Configure OpenMP thread count to avoid oversubscription when running many MPI ranks.
# Use SLURM's cpus-per-task when available; fallback to 1.
if [ -n "$SLURM_CPUS_PER_TASK" ]; then
    export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
else
    export OMP_NUM_THREADS=1
fi
echo "Set OMP_NUM_THREADS=$OMP_NUM_THREADS"

# Create working directory
WORK_DIR="test_case_v85_s2000"
mkdir -p $WORK_DIR
cd $WORK_DIR

# Hybrid MPI+OpenMP configuration: ranks-per-node * cpus-per-task should match cores per node
# RANKS_PER_NODE is used to instruct mpirun how to map ranks across sockets
RANKS_PER_NODE=24

# Copy base case
cp -r $HOME/Isambaseball/openfoam_case/* .

# Ensure functions/forceCoeffs (contains pitchAxis) is present in working case
if [ -f "$HOME/Isambaseball/openfoam_case/system/functions/forceCoeffs" ]; then
    mkdir -p system/functions
    cp -f "$HOME/Isambaseball/openfoam_case/system/functions/forceCoeffs" system/functions/forceCoeffs
    echo "Copied system/functions/forceCoeffs into case"
fi

# ControlDict is edited directly in openfoam_case/system/controlDict
# The job script will no longer modify controlDict at runtime so the case
# uses the authoritative settings in the case directory.

# ...existing code...

# Install pre-generated mesh
if [ -f "$HOME/Isambaseball/master_mesh.tar.gz" ]; then
    echo "Installing pre-generated mesh..."
    tar -xzf $HOME/Isambaseball/master_mesh.tar.gz
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
/*--------------------------------*- C++ -*----------------------------------*/
FoamFile
{
    format      ascii;
    class       volScalarField;
    location    "0";
    object      omega;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

dimensions      [0 0 -1 0 0 0 0];

// scalar omega (rad/s)
internalField   uniform $OMEGA;

boundaryField
{
    inlet
    {
        type            fixedValue;
        value           uniform $OMEGA;
    }
    
    outlet
    {
        type            zeroGradient;
    }
    
    sides
    {
        type            slip;
    }
    
    top
    {
        type            slip;
    }
    
    bottom
    {
        type            slip;
    }
    
    baseball
    {
        type            omegaWallFunction;
        value           uniform $OMEGA;
    }
}
EOF
# Normalize line endings (remove CR from CRLF) in case this script was edited on Windows
if command -v dos2unix >/dev/null 2>&1; then
    dos2unix 0/omega >/dev/null 2>&1 || true
else
    # portable fallback: remove \r characters
    sed -i 's/\r$//' 0/omega || true
fi
# Remove UTF-8 BOM if present (some editors add BOM which breaks OpenFOAM headers)
awk 'NR==1{sub(/\xef\xbb\xbf/,"")};{print}' 0/omega > 0/omega.nobom || true
mv -f 0/omega.nobom 0/omega || true
# Ensure file ends with a newline
tail -c1 0/omega | od -An -t u1 | grep -q . || echo >> 0/omega || true

# Sanitize initial field files in 0/ to avoid malformed headers (CRLF, BOM, missing location)
echo "Sanitizing 0/ field files to avoid header/encoding issues..."
sanitize_field() {
    f="$1"
    [ -f "$f" ] || return
    # remove UTF-8 BOM if present
    awk 'NR==1{sub(/\xef\xbb\xbf/,"")};{print}' "$f" > "$f.nobom" || true
    mv -f "$f.nobom" "$f" || true
    # remove CR from CRLF
    if command -v dos2unix >/dev/null 2>&1; then
        dos2unix "$f" >/dev/null 2>&1 || true
    else
        sed -i 's/\r$//' "$f" || true
    fi
    # ensure file ends with newline
    tail -c1 "$f" | od -An -t u1 | grep -q . || echo >> "$f" || true

    # If file contains a FoamFile but no location, insert location "0" after class line
    if grep -q "FoamFile" "$f" 2>/dev/null; then
        if ! grep -q 'location[[:space:]]*"0"' "$f" 2>/dev/null; then
            awk 'BEGIN{inF=0;done=0} /FoamFile/{print; inF=1; next} inF && /class[[:space:]]+/ && !done {print; print "    location    \"0\";"; done=1; next} {print}' "$f" > "$f.tmp" && mv -f "$f.tmp" "$f" || true
        fi
    fi
}

for fld in 0/*; do
    sanitize_field "$fld"
done

# Domain decomposition for parallel execution
echo "Decomposing domain for $SLURM_NTASKS cores..."
# Use -copyZero to copy the entire 0/ directory into processor*/0/ instead of
# attempting to pass a fields list (decomposePar's -fields flag does not take
# an argument and passing one caused a FOAM error). -copyZero is the safest
# way to ensure custom initial fields like omega are present on each processor.
# Generate a decomposeParDict for this run so numberOfSubdomains matches SLURM
cat > system/decomposeParDict << DOC
/* Decomposition dictionary auto-generated by 2_test_case.sh */
FoamFile
{
    format      ascii;
    class       dictionary;
    location    "system";
    object      decomposeParDict;
}

numberOfSubdomains $SLURM_NTASKS;

method          scotch;

scotchCoeffs
{
    // Let scotch choose defaults; processorWeights omitted to avoid mismatch
}

distributed     no;
roots           ( );

// End of auto-generated file
DOC

decomposePar > log.decomposePar 2>&1

if [ $? -ne 0 ]; then
    echo "ERROR: Domain decomposition failed"
    exit 1
fi

# Ensure custom fields (like omega) are present in each processor directory
# DecomposePar sometimes does not create or copy custom files; copy 0/omega
# into each processor*/0/omega so parallel runs can read them.
echo "Distributing custom initial fields to processor folders..."
for pd in processor[0-9]*; do
    if [ -d "$pd" ]; then
        mkdir -p "$pd/0"
        if [ -f "0/omega" ]; then
            if [ -f "$pd/0/omega" ]; then
                echo "  updating $pd/0/omega (preserving procBoundary entries)"
                # Replace internalField uniform value
                sed -i -E 's/^(\s*internalField\s+uniform\s+)[^;]+;/\1'$OMEGA';/' "$pd/0/omega" || true
                # Replace value uniform lines (inlet, baseball, etc.) with scalar omega
                sed -i -E 's/^(\s*value\s+uniform\s+)[^;]+;/\1'$OMEGA';/' "$pd/0/omega" || true
                # Normalize and remove BOM if any
                awk 'NR==1{sub(/\xef\xbb\xbf/,"")};{print}' "$pd/0/omega" > "$pd/0/omega.nobom" || true
                mv -f "$pd/0/omega.nobom" "$pd/0/omega" || true
                if command -v dos2unix >/dev/null 2>&1; then
                    dos2unix "$pd/0/omega" >/dev/null 2>&1 || true
                else
                    sed -i 's/\r$//' "$pd/0/omega" || true
                fi
                tail -c1 "$pd/0/omega" | od -An -t u1 | grep -q . || echo >> "$pd/0/omega" || true
            else
                cp -f "0/omega" "$pd/0/omega"
                echo "  copied 0/omega -> $pd/0/omega"
                # Remove BOM from processor copy and normalize line endings
                awk 'NR==1{sub(/\xef\xbb\xbf/,"")};{print}' "$pd/0/omega" > "$pd/0/omega.nobom" || true
                mv -f "$pd/0/omega.nobom" "$pd/0/omega" || true
                if command -v dos2unix >/dev/null 2>&1; then
                    dos2unix "$pd/0/omega" >/dev/null 2>&1 || true
                else
                    sed -i 's/\r$//' "$pd/0/omega" || true
                fi
                tail -c1 "$pd/0/omega" | od -An -t u1 | grep -q . || echo >> "$pd/0/omega" || true
            fi
        fi
    fi
done

# Verification: ensure processor directories contain expected initial fields
echo "Verifying processor initial fields and boundary entries..."
missing=0
for pd in processor[0-9]*; do
    if [ -d "$pd" ]; then
        # Check omega presence
        if [ ! -f "$pd/0/omega" ]; then
            echo "MISSING: $pd/0/omega"
            missing=1
        fi

        # Quick check for procBoundary entries in p boundaryField (common cause of IO errors)
        if [ -f "$pd/0/p" ]; then
            if ! grep -q "procBoundary" "$pd/0/p"; then
                echo "WARNING: $pd/0/p does not contain procBoundary entries; decomposition may have been done incorrectly"
                missing=1
            fi
        else
            echo "MISSING: $pd/0/p"
            missing=1
        fi
    fi
done

if [ $missing -ne 0 ]; then
    echo "ERROR: Missing initial fields or processor boundary entries detected. Inspect log.decomposePar and processor*/constant/polyMesh/boundary"

    print_diagnostics() {
        echo
        echo "=========================================="
        echo "CASE DIAGNOSTICS"
        echo "=========================================="

        echo "--- last 200 lines of log.decomposePar ---"
        if [ -f log.decomposePar ]; then
            tail -n 200 log.decomposePar
        else
            echo "(no log.decomposePar found)"
        fi

        echo
        echo "--- processor directories (summary) ---"
        ls -ld processor* 2>/dev/null || echo "(no processor* directories found)"

        echo
        echo "--- Search for procBoundary in processor*/0/p ---"
        for pd in processor[0-9]*; do
            if [ -d "$pd" ]; then
                echo ">>> $pd/0/p <<<"
                if [ -f "$pd/0/p" ]; then
                    echo "(first 200 lines)"
                    sed -n '1,200p' "$pd/0/p" || true
                    echo "(lines containing procBoundary)"
                    grep -n "procBoundary" "$pd/0/p" || echo "(none)"
                else
                    echo "(missing)"
                fi
                echo
            fi
        done

        echo "--- processor constant/polyMesh/boundary files ---"
        for pd in processor[0-9]*; do
            if [ -d "$pd" ]; then
                b="$pd/constant/polyMesh/boundary"
                echo ">>> $b <<<"
                if [ -f "$b" ]; then
                    sed -n '1,200p' "$b" || true
                else
                    echo "(missing)"
                fi
                echo
            fi
        done

        echo "--- Original 0/ field snippets ---"
        for f in omega U p; do
            echo ">>> 0/$f <<<"
            if [ -f "0/$f" ]; then
                sed -n '1,120p' "0/$f" || true
            else
                echo "(missing)"
            fi
            echo
        done

        echo "--- Quick grep for procBoundary across processors ---"
        grep -R --line-number "procBoundary" processor* 2>/dev/null || echo "(no procBoundary found in any processor*/0/p)"

        echo "=========================================="
    }

    print_diagnostics
    exit 1
fi


# Run simulation
echo "Starting CFD simulation..."
echo "This will take 2-4 hours..."

echo "Running pimpleFoam with hybrid MPI+OpenMP: OMP_NUM_THREADS=$OMP_NUM_THREADS"
# Compute total ranks (if SLURM_NTASKS is not set, use NNODES*RANKS_PER_NODE)
if [ -n "$SLURM_NTASKS" ]; then
    TOTAL_RANKS=$SLURM_NTASKS
elif [ -n "$SLURM_NNODES" ]; then
    TOTAL_RANKS=$(( SLURM_NNODES * RANKS_PER_NODE ))
else
    # fallback to single node RANKS_PER_NODE
    TOTAL_RANKS=$RANKS_PER_NODE
fi

# Use ppr (processes-per-resource) mapping to distribute ranks across sockets and
# bind ranks to cores; export OMP_NUM_THREADS to processes.
mpirun --report-bindings --bind-to core --map-by ppr:${RANKS_PER_NODE}:socket -x OMP_NUM_THREADS -np $TOTAL_RANKS pimpleFoam -parallel > log.pimpleFoam 2>&1

if [ $? -ne 0 ]; then
    echo "ERROR: Simulation failed"
    tail -50 log.pimpleFoam
    exit 1
fi

# Reconstruct parallel results
echo "Reconstructing results..."
reconstructPar > log.reconstructPar 2>&1

# Ensure postProcessing function-object outputs exist for the whole time series.
# Running foamPostProcess without -latestTime will process all time directories
# and (re)generate postProcessing/forceCoeffs/0/coefficient.dat and
# postProcessing/forces/0/force.dat from the reconstructed data.
echo "Running foamPostProcess to (re)generate postProcessing force outputs..."
if command -v foamPostProcess >/dev/null 2>&1; then
    foamPostProcess -func forces > log.foamPostProcess 2>&1 || true
    foamPostProcess -func forceCoeffs >> log.foamPostProcess 2>&1 || true
    echo "foamPostProcess complete (see log.foamPostProcess)"
else
    echo "Warning: foamPostProcess not found in PATH; skipping post-processing step"
fi

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
