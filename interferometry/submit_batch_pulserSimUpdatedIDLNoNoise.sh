#!/bin/bash

#SBATCH --job-name=getPol
#SBATCH --nodes=1 --ntasks-per-node=1
#SBATCH --account=PAS0654
#SBATCH --mail-type=FAIL
#SBATCH --time=20:00:00
#SBATCH --output=run_Pol.log   # Standard output and error log

eval 'source /users/PAS0654/jflaherty13/.bashrc' #source yours
eval 'source /users/PAS0654/jflaherty13/.bash_profile'
eval 'cvmfs'
export XDG_RUNTIME_DIR=/users/PAS0654/jflaherty13/araAnalysis/araRecoAndSourceSearch/ARA_Reconstruction/brianInterferometry/temp/
export RUNLEVEL=3
export QT_QPA_PLATFORM='offscreen'
cd /users/PAS0654/jflaherty13/araAnalysis/araRecoAndSourceSearch/ARA_Reconstruction/brianInterferometry #go to wherever you have the code

index=${SLURM_ARRAY_TASK_ID}

# ./calculatearrivaldir_interf 2 ${index} /users/PAS0654/jflaherty13/source/AraSim/outputs/updatedIdlPulseModel/AraOut.setup_variablePsi.txt.run$index.root /users/PAS0654/jflaherty13/source/AraSim/outputs/updatedIdlPulseModel/interferometry

# ./calculatearrivaldir_interf 2 ${index} /users/PAS0654/jflaherty13/source/AraSim/outputs/updatedIdlPulseModelNoNoise/AraOut.setup_variablePsi.txt.run$index.root /users/PAS0654/jflaherty13/source/AraSim/outputs/updatedIdlPulseModelNoNoise/interferometry

# ./calculatearrivaldir_interf 2 ${index} /users/PAS0654/jflaherty13/source/AraSim/outputs/updatedIdlPulseModelHpolGain/AraOut.setup_variablePsi.txt.run$index.root /users/PAS0654/jflaherty13/source/AraSim/outputs/updatedIdlPulseModelHpolGain/interferometry

# ./calculatearrivaldir_interf 2 ${index} /users/PAS0654/jflaherty13/source/AraSim/outputs/alisaIdlPulseModelNoNoise/AraOut.setup_variablePsi.txt.run$index.root /users/PAS0654/jflaherty13/source/AraSim/outputs/alisaIdlPulseModelNoNoise/interferometry

# ./calculatearrivaldir_interf 2 ${index} /users/PAS0654/jflaherty13/source/AraSim/outputs/alisaIdlPulseModel/AraOut.setup_variablePsi.txt.run$index.root /users/PAS0654/jflaherty13/source/AraSim/outputs/alisaIdlPulseModel/interferometry

# ./calculatearrivaldir_interf 2 ${index} /users/PAS0654/jflaherty13/source/AraSim/outputs/updatedIdlPulseModelHpolGainNoNoise/AraOut.setup_variablePsi\=${index}.txt.run$index.root /users/PAS0654/jflaherty13/source/AraSim/outputs/updatedIdlPulseModelHpolGainNoNoise/interferometry

# ./calculatearrivaldir_interf 2 ${index} /users/PAS0654/jflaherty13/source/AraSim/outputs/updatedIdlPulseModelHpolGain/AraOut.setup_variablePsi.txt.run$index.root test/

# ./calculatearrivaldir_interf 2 ${index} /users/PAS0654/jflaherty13/source/AraSim/outputs/updatedIdlPulseModelHpolGain/AraOut.setup_variablePsi.txt.run$index.root /users/PAS0654/jflaherty13/source/AraSim/outputs/updatedIdlPulseModelHpolGainNoNoiseNoResponse/interferometry

# ./calculatearrivaldir_interf 2 ${index} /users/PAS0654/jflaherty13/source/AraSim/outputs/test/AraOut.setup_variablePsi\=$index.txt.run$index.root /users/PAS0654/jflaherty13/source/AraSim/outputs/test/interferometry


# ##Testing deconvolution
# ./calculatearrivaldir_interf 2 ${index} /users/PAS0654/jflaherty13/source/AraSim/outputs/testNeutrino/AraOut.setup_variablePsi\=$index.txt.run$index.root /users/PAS0654/jflaherty13/source/AraSim/outputs/testNeutrino/interferometry


# ##20230710 Pulser Sims
# ./calculatearrivaldir_interf 2 ${index} /users/PAS0654/jflaherty13/source/AraSim/outputs/20230710_pulserSims_1000m/withNoise/AraOut.setup_variablePsi\=$index.txt.run$index.root /users/PAS0654/jflaherty13/source/AraSim/outputs/20230710_pulserSims_1000m/interferometry


##20230710 Pulser Sims but with updated SNR calculation in the vertex reco, so we're saving it to another location for testing.
./calculatearrivaldir_interf 2 ${index} /users/PAS0654/jflaherty13/source/AraSim/outputs/20230710_pulserSims_1000m/withNoise/AraOut.setup_variablePsi\=$index.txt.run$index.root /users/PAS0654/jflaherty13/source/AraSim/outputs/20230710_pulserSims_1000m/interferometryNewSnr