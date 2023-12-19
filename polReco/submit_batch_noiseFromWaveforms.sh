#!/bin/bash

#SBATCH --job-name=getPol
#SBATCH --nodes=1 --ntasks-per-node=1
#SBATCH --account=PAS0654
#SBATCH --mail-type=FAIL
#SBATCH --time=2:00:00
#SBATCH --output=run_Pol.log   # Standard output and error log

eval 'source /users/PAS0654/jflaherty13/.bashrc' #source yours
eval 'source /users/PAS0654/jflaherty13/.bash_profile'
eval 'cvmfs'
export XDG_RUNTIME_DIR=/users/PAS0654/jflaherty13/araAnalysis/thesisWork/polReco/temp/
export RUNLEVEL=3
export QT_QPA_PLATFORM='offscreen'
cd /users/PAS0654/jflaherty13/araAnalysis/thesisWork/polReco #go to wherever you have the code

run=${SLURM_ARRAY_TASK_ID}

# python noiseFromWaveforms.py $run ~/source/AraSim/outputs/updatedIdlPulseModelHpolGain/AraOut.setup_variablePsi.txt.run${run}.root noiseFromWaveforms/

# Get noise from simulated pulser waveforms
python noiseFromWaveforms.py $run /users/PAS0654/jflaherty13/source/AraSim/outputs/20230710_pulserSims_1000m/withNoise/deconvolution/deconvolvedWaveforms_run_$run.root /users/PAS0654/jflaherty13/source/AraSim/outputs/20230710_pulserSims_1000m/withNoise/deconvolution/noiseFromWaveforms/