#!/bin/bash

#SBATCH --job-name=baseball_minimal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --time=1-00:00:00
#SBATCH --partition=grace
#SBATCH --output=test_case_minimal_%j.out
#SBATCH --error=test_case_minimal_%j.err

# Minimal OpenFOAM parallel run script
# - Writes a canonical ASCII volVectorField 0/omega
# - Runs decomposePar, redistributePar -overwrite
# - Launches the parallel solver via srun (falls back to mpirun)
# This intentionally avoids defensive workarounds; it assumes the case and
# cluster are configured correctly for OpenFOAM parallel runs.

set -euo pipefail

WORK_DIR="test_case_v85_s2000_minimal"
mkdir -p "$WORK_DIR"
cd "$WORK_DIR"

# copy base case
cp -r /$HOME/Isambaseball/openfoam_case/* .

# install mesh if provided
if [ -f "/$HOME/Isambaseball/master_mesh.tar.gz" ]; then
    tar -xzf /$HOME/Isambaseball/master_mesh.tar.gz
    cp -r master_mesh/* constant/
    rm -rf master_mesh
fi

# Parameters
VELOCITY=38.0
SPIN_RATE=2000
# compute omega with bc -l and format to fixed 6 decimals to avoid unexpected tokens
OMEGA=$(echo "$SPIN_RATE * 2 * 3.14159 / 60" | bc -l)
OMEGA=$(printf '%.6f' "$OMEGA")

# Update simple U and magUInf if present
if [ -f 0/U ]; then
    sed -i "s/uniform (38.0 0 0)/uniform ($VELOCITY 0 0)/" 0/U || true
fi
if [ -f system/controlDict ]; then
    sed -i "s/magUInf.*38.0;/magUInf         $VELOCITY;/" system/controlDict || true
fi

# Write canonical ASCII volScalarField 0/omega (turbulence field)
mkdir -p 0
cat > 0/omega <<EOF
/*--------------------------------*- C++ -*----------------------------------*\\
| =========                 |                                                 |
| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\    /   O peration     | Version:  1.0                                   |
|   \\  /    A nd           | Website:  https://openfoam.org                  |
|    \/     M anipulation  |                                                 |
\\*---------------------------------------------------------------------------*/
FoamFile
{
    version     2.0;
    format      ascii;
    class       volScalarField;
    object      omega;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

dimensions      [0 0 -1 0 0 0 0];

# scalar omega used by turbulence models
internalField   uniform $OMEGA;

boundaryField
{
    inlet { type fixedValue; value uniform $OMEGA; }
    outlet { type zeroGradient; }
    sides { type zeroGradient; }
    top { type zeroGradient; }
    bottom { type zeroGradient; }
    baseball
    {
        type            omegaWallFunction;
        value           uniform $OMEGA;
    }
}

EOF

# Also update the rotatingWallVelocity omega in 0/U so the boundary matches the requested spin
if [ -f 0/U ]; then
    # Replace a numeric omega value inside the U file (keeps formatting/comments)
    sed -i -E "s/(^\s*omega\s+)[0-9]+(\.[0-9]+)?;/\1$OMEGA;/" 0/U || true
fi

# Ensure ascii line endings
if command -v dos2unix >/dev/null 2>&1; then
    dos2unix 0/omega || true
else
    sed -i 's/\r$//' 0/omega || true
fi

# Basic mesh check
if [ ! -f constant/polyMesh/boundary ]; then
    echo "ERROR: Mesh constant/polyMesh/boundary not found. Aborting."
    exit 1
fi

# decompose
# show the generated 0/omega for debugging (before decomposition)
echo "--- generated 0/omega ---"
sed -n '1,40p' 0/omega || true
echo "------------------------"

decomposePar > log.decomposePar 2>&1

# redistribute (ensure fields are placed in processor directories)
redistributePar -overwrite > log.redistributePar 2>&1

# show a short preview
echo "--- preview processor0/0/omega (if present) ---"
if [ -f processor0/0/omega ]; then
    head -n 40 processor0/0/omega || true
else
    echo "processor0/0/omega not present after redistributePar"
fi
echo "-----------------------------------------------"

# Launch parallel solver (use srun if available)
NTASKS=${SLURM_NTASKS:-24}
LAUNCH_RC=0
if command -v srun >/dev/null 2>&1; then
    echo "Launching with srun on $NTASKS tasks..."
    srun --ntasks=$NTASKS foamRun -solver incompressibleFluid -parallel > log.foamRun.srun 2>&1 || LAUNCH_RC=$?
    if [ $LAUNCH_RC -ne 0 ]; then
        echo "srun launch failed (code $LAUNCH_RC); falling back to mpirun"
        mpirun -np $NTASKS foamRun -solver incompressibleFluid -parallel > log.foamRun.mpirun 2>&1 || LAUNCH_RC=$?
    fi
else
    echo "srun not found; launching with mpirun on $NTASKS ranks..."
    mpirun -np $NTASKS foamRun -solver incompressibleFluid -parallel > log.foamRun.mpirun 2>&1 || LAUNCH_RC=$?
fi

if [ $LAUNCH_RC -ne 0 ]; then
    echo "ERROR: Parallel solver failed (exit $LAUNCH_RC). Check log.foamRun.*"
    tail -n 200 log.foamRun.srun 2>/dev/null || true
    tail -n 200 log.foamRun.mpirun 2>/dev/null || true
    exit 1
fi

# reconstruct
reconstructPar > log.reconstructPar 2>&1

# post-process forces (attempt)
foamPostProcess -func forces -latestTime > log.postProcess 2>&1 || true
foamPostProcess -func forceCoeffs -latestTime >> log.postProcess 2>&1 || true

echo "Minimal parallel run script completed. Check logs and postProcessing for results." 
