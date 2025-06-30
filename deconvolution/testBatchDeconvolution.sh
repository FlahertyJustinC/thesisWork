#!/bin/bash

#SBATCH --job-name=decon
#SBATCH --nodes=1 --ntasks-per-node=1
#SBATCH --account=PAS0654
#SBATCH --mail-type=FAIL
#SBATCH --time=2:00:00
#SBATCH --array=42-160
#SBATCH --output=run_Pol.log   # Standard output and error log

eval 'source /users/PAS0654/jflaherty13/.bashrc' #source yours
eval 'source /users/PAS0654/jflaherty13/.bash_profile'
eval 'cvmfs'
export XDG_RUNTIME_DIR=/users/PAS0654/jflaherty13/source/AraSim/temp/
export RUNLEVEL=3
export QT_QPA_PLATFORM='offscreen'
cd /users/PAS0654/jflaherty13/source/AraSim #go to wherever you have the code

index=$((10*$SLURM_ARRAY_TASK_ID))

rawDirectory="/users/PAS0654/jflaherty13/source/AraSim/outputs/20240321_AlanPulserSims_A24/AraSimOutput/"
recoDirectory="/users/PAS0654/jflaherty13/source/AraSim/outputs/20240321_AlanPulserSims_A24/interferometry/"
outDirectory="/users/PAS0654/jflaherty13/araAnalysis/thesisWork/deconvolution/testBatch"

#A2 config 6 without and with birefringence
station=2
config=6

# if ! test -f /users/PAS0654/jflaherty13/source/AraSim/outputs/20240321_AlanPulserSims_A24/deconvolution/A${station}/birefringence/deconvolvedWaveforms_run_${index}.root; then

./deconvolveWaveform $station $config $index ${rawDirectory}SPICE_birefringence_forJustin/AraOut.setup_birefringence_A$station\=$index.txt.run$index.root ${recoDirectory}/A$station/birefringence/recangle_reco_out_run_$index.root $outDirectory/A$station/birefringence ${rawDirectory}/../setup_birefringence_A$station.txt
    
# fi

# if ! test -f /users/PAS0654/jflaherty13/source/AraSim/outputs/20240321_AlanPulserSims_A24/interferometry/A${station}/noBirefringence/deconvolvedWaveforms_run_${index}.root; then

./deconvolveWaveform $station $config $index ${rawDirectory}SPICE_birefringence_forJustin/AraOut.setup_noBirefringence_A$station\=$index.txt.run$index.root ${recoDirectory}/A$station/noBirefringence/recangle_reco_out_run_$index.root $outDirectory/A$station/noBirefringence ${rawDirectory}/../setup_noBirefringence_A$station.txt

# fi



#A4 config 4 with and without birefringence
station=4
config=4
# if ! test -f /users/PAS0654/jflaherty13/source/AraSim/outputs/20240321_AlanPulserSims_A24/deconvolution/A${station}/birefringence/deconvolvedWaveforms_run_${index}.root; then

./deconvolveWaveform $station $config $index ${rawDirectory}SPICE_birefringence_forJustin/AraOut.setup_birefringence_A$station\=$index.txt.run$index.root ${recoDirectory}/A$station/birefringence/recangle_reco_out_run_$index.root $outDirectory/A$station/birefringence ${rawDirectory}/../setup_birefringence_A$station.txt
    
# fi

# if ! test -f /users/PAS0654/jflaherty13/source/AraSim/outputs/20240321_AlanPulserSims_A24/interferometry/A${station}/noBirefringence/deconvolvedWaveforms_run_${index}.root; then

./deconvolveWaveform $station $config $index ${rawDirectory}SPICE_birefringence_forJustin/AraOut.setup_noBirefringence_A$station\=$index.txt.run$index.root ${recoDirectory}/A$station/noBirefringence/recangle_reco_out_run_$index.root $outDirectory/A$station/noBirefringence ${rawDirectory}/../setup_noBirefringence_A$station.txt

# fi

