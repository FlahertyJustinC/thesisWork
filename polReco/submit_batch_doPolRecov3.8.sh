#!/bin/bash

#SBATCH --job-name=getPol
#SBATCH --nodes=1 --ntasks-per-node=1
#SBATCH --account=PAS0654
#SBATCH --mail-type=FAIL
#SBATCH --time=8:00:00
#SBATCH --output=run_Pol.log   # Standard output and error log

eval 'source /users/PAS0654/jflaherty13/.bashrc' #source yours
eval 'source /users/PAS0654/jflaherty13/.bash_profile'
eval 'cvmfs'
export XDG_RUNTIME_DIR=/users/PAS0654/jflaherty13/araAnalysis/araRecoAndSourceSearch/SPIceCoreReco/temp/
export RUNLEVEL=3
export QT_QPA_PLATFORM='offscreen'
cd /users/PAS0654/jflaherty13/araAnalysis/thesisWork/polReco/ #go to wherever you have the code

station=2
psi=${SLURM_ARRAY_TASK_ID}

# python doPolReco_v3.8.py $station $psi 0 /users/PAS0654/jflaherty13/source/AraSim/outputs/alisaIdlPulseModel/interferometry/recangle_reco_out_run_${psi}.root /users/PAS0654/jflaherty13/source/AraSim/outputs/alisaIdlPulseModel/AraOut.setup_variablePsi.txt.run${psi}.root output/A${station}/ 1 0 0 0 1

# python doPolReco_v3.8.py $station $psi 0 /users/PAS0654/jflaherty13/source/AraSim/outputs/alisaIdlPulseModelNoNoise/interferometry/recangle_reco_out_run_${psi}.root /users/PAS0654/jflaherty13/source/AraSim/outputs/alisaIdlPulseModelNoNoise/AraOut.setup_variablePsi.txt.run${psi}.root output/A${station}/NoNoise/ 1 0 0 0 1

# python doPolReco_v3.8.py $station $psi 0 /users/PAS0654/jflaherty13/source/AraSim/outputs/updatedIdlPulseModel/interferometry/recangle_reco_out_run_${psi}.root /users/PAS0654/jflaherty13/source/AraSim/outputs/updatedIdlPulseModel/AraOut.setup_variablePsi.txt.run${psi}.root output/A${station}/ 1 0 0 0 1

# python doPolReco_v3.8.py $station $psi 0 /users/PAS0654/jflaherty13/source/AraSim/outputs/updatedIdlPulseModelNoNoise/interferometry/recangle_reco_out_run_${psi}.root /users/PAS0654/jflaherty13/source/AraSim/outputs/updatedIdlPulseModelNoNoise/AraOut.setup_variablePsi.txt.run${psi}.root output/A${station}/NoNoise/ 1 0 0 0 1

# python doPolReco_v3.8.py $station $psi 0 /users/PAS0654/jflaherty13/source/AraSim/outputs/updatedIdlPulseModelHpolGain/interferometry/recangle_reco_out_run_${psi}.root /users/PAS0654/jflaherty13/source/AraSim/outputs/updatedIdlPulseModelHpolGain/AraOut.setup_variablePsi.txt.run${psi}.root output/A${station}/HpolGain/ 1 0 0 0 1

# python doPolReco_v3.8.py $station $psi 0 /users/PAS0654/jflaherty13/source/AraSim/outputs/updatedIdlPulseModelHpolGain/interferometry/recangle_reco_out_run_${psi}.root /users/PAS0654/jflaherty13/source/AraSim/outputs/updatedIdlPulseModelHpolGain/AraOut.setup_variablePsi.txt.run${psi}.root output/A${station}/HpolGain/ 1 0 1 1 1

# python doPolReco_v3.8.py $station $psi 0 /users/PAS0654/jflaherty13/source/AraSim/outputs/updatedIdlPulseModelHpolGainNoNoise/interferometry/recangle_reco_out_run_${psi}.root /users/PAS0654/jflaherty13/source/AraSim/outputs/updatedIdlPulseModelHpolGainNoNoise/AraOut.setup_variablePsi.txt.run${psi}.root output/A${station}/HpolGainNoNoise/ 1 0 1 1 1


# # Test polReco code on new deconvolution data with no Noise
# python doPolReco_v3.8.py $station $psi 0 /users/PAS0654/jflaherty13/source/AraSim/outputs/20230710_pulserSims_1000m/interferometry/recangle_reco_out_run_${psi}.root /users/PAS0654/jflaherty13/source/AraSim/outputs/20230710_pulserSims_1000m/deconvolution/deconvolvedWaveforms_run_${psi}.root output/A${station}/testDeconvolvedWaveforms/ 1 0 0 0 0

# Test polReco code on new deconvolution data with noise but no noise subtraction
# python doPolReco_v3.8.py $station $psi 0 /users/PAS0654/jflaherty13/source/AraSim/outputs/20230710_pulserSims_1000m/interferometry/recangle_reco_out_run_${psi}.root /users/PAS0654/jflaherty13/source/AraSim/outputs/20230710_pulserSims_1000m/deconvolution/deconvolvedWaveforms_run_${psi}.root output/A${station}/testDeconvolvedWaveforms/ 1 0 0 0 0

# Test deconvolution script on SpiceCore data Note: psi=subset number in this case.
# echo "Opening python file"
# python doPolReco_v3.8.py $station 12559 $psi /users/PAS0654/jflaherty13/araAnalysis/araRecoAndSourceSearch/ARA_Reconstruction/brianInterferometry/spicecoreInterfAllRuns/A2/run_012559/recangle_reco_out_run_12559_${psi}.root /users/PAS0654/jflaherty13/source/AraSim/outputs/20230713_SpiceReco/deconvolvedWaveforms_run_12559_${psi}.root output/A${station}/testDeconvolvedWaveformsSpice/ 2 0 0 0 0

#Doing polReco on simulations with NFOUR=8192
# python doPolReco_v3.8.py $station $psi 0 /users/PAS0654/jflaherty13/source/AraSim/outputs/20230710_pulserSims_1000m/interferometry/recangle_reco_out_run_${psi}.root /users/PAS0654/jflaherty13/source/AraSim/outputs/20230716_pulserSims_1000m_NFOUR_8192/deconvolution/deconvolvedWaveforms_run_${psi}.root output/A${station}/20230716_pulserSims_1000m_NFOUR_8192/ 1 0 0 0 0


#20230710 pulser simulations using 45mV noise in vertex Reco
# python doPolReco_v3.8.py $station $psi 0 /users/PAS0654/jflaherty13/source/AraSim/outputs/20230710_pulserSims_1000m/interferometry/recangle_reco_out_run_${psi}.root /users/PAS0654/jflaherty13/source/AraSim/outputs/20230710_pulserSims_1000m/noNoise/deconvolution/deconvolvedWaveforms_run_${psi}.root /users/PAS0654/jflaherty13/source/AraSim/outputs/20230710_pulserSims_1000m/noNoise/polReco/ 1 0 0 0 0

# python doPolReco_v3.8.py $station $psi 0 /users/PAS0654/jflaherty13/source/AraSim/outputs/20230710_pulserSims_1000m/interferometry/recangle_reco_out_run_${psi}.root /users/PAS0654/jflaherty13/source/AraSim/outputs/20230710_pulserSims_1000m/withNoise/deconvolution/deconvolvedWaveforms_run_${psi}.root /users/PAS0654/jflaherty13/source/AraSim/outputs/20230710_pulserSims_1000m/withNoise/polReco/ 1 0 0 0 0

# 20230710 pulser simulations using per-channel noise in vertex Reco
# python doPolReco_v3.8.py $station $psi 0 /users/PAS0654/jflaherty13/source/AraSim/outputs/20230710_pulserSims_1000m/interferometryNewSnr/recangle_reco_out_run_${psi}.root /users/PAS0654/jflaherty13/source/AraSim/outputs/20230710_pulserSims_1000m/noNoise/deconvolutiondeconvolvedWaveforms_run_${psi}.root /users/PAS0654/jflaherty13/source/AraSim/outputs/20230710_pulserSims_1000m/noNoise/polRecoNewSnr/ 1 0 0 0 0


#For ICRC
# python doPolReco_v3.8.py $station $psi 0 /users/PAS0654/jflaherty13/source/AraSim/outputs/20230710_pulserSims_1000m/interferometryNewSnr/recangle_reco_out_run_${psi}.root /users/PAS0654/jflaherty13/source/AraSim/outputs/20230710_pulserSims_1000m/withNoise/deconvolution/deconvolvedWaveforms_run_${psi}.root /users/PAS0654/jflaherty13/source/AraSim/outputs/20230710_pulserSims_1000m/withNoise/polRecoNewSnr/ 1 0 0 0 0

#Testing NFOUR = 4096 in deconvolution
# python doPolReco_v3.8.py $station $psi 0 /users/PAS0654/jflaherty13/source/AraSim/outputs/20230710_pulserSims_1000m/interferometryNewSnr/recangle_reco_out_run_${psi}.root /users/PAS0654/jflaherty13/source/AraSim/outputs/20230710_pulserSims_1000m/withNoise/deconvolution_NFOUR\=4096/deconvolvedWaveforms_run_${psi}.root /users/PAS0654/jflaherty13/source/AraSim/outputs/20230710_pulserSims_1000m/withNoise/polRecoNewSnr_NFOUR\=4096/ 1 0 0 0 0

#Pol Reco on SpiceCore data after new deconvolution 10/11/2023
# python doPolReco_v3.8.py $station 12559 $psi /users/PAS0654/jflaherty13/araAnalysis/araRecoAndSourceSearch/ARA_Reconstruction/brianInterferometry/spicecoreInterfAllRuns/A2/run_012559/recangle_reco_out_run_12559_${psi}.root /users/PAS0654/jflaherty13/source/AraSim/outputs/20230713_SpiceReco/deconvolvedWaveforms_run_12559_${psi}.root output/A${station}/testDeconvolvedWaveformsSpice_20231011_withTimeshift/ 1 0 0 0 0

#Pol Reco on SpiceCore data without deconvolution as a sanity check 10/15/2023
# python doPolReco_v3.8.py $station 12559 $psi /users/PAS0654/jflaherty13/araAnalysis/araRecoAndSourceSearch/ARA_Reconstruction/brianInterferometry/spicecoreInterfAllRuns/A2/run_012559/recangle_reco_out_run_12559_${psi}.root /users/PAS0654/jflaherty13/araAnalysis/thesisWork/data/A2/run_012559/split/event012559__${psi}.root output/A${station}/polRecoRaw12559_20231015/ 0 0 0 0 0


# #Pol Reco on Spicecore data using the same parameters as 10/11/2023, but with proper tagging of calpulsers and soft triggers.
# python doPolReco_v3.8.py $station 12559 $psi /users/PAS0654/jflaherty13/araAnalysis/araRecoAndSourceSearch/ARA_Reconstruction/brianInterferometry/spicecoreInterfAllRuns/A2/run_012559/recangle_reco_out_run_12559_${psi}.root /users/PAS0654/jflaherty13/source/AraSim/outputs/20231018_SpiceReco/deconvolvedWaveforms_run_12559_${psi}.root output/A${station}/testDeconvolvedWaveformsSpice_20231018_updatedEventTags/ 1 0 0 0 0

#Pol Reco on Spice data with narrower band pass filters at Amy's request - 11/2/2023
# freqMin=150
# freqMax=200
# python doPolReco_v3.8.py $station 12559 $psi /users/PAS0654/jflaherty13/araAnalysis/araRecoAndSourceSearch/ARA_Reconstruction/brianInterferometry/spicecoreInterfAllRuns/A2/run_012559/recangle_reco_out_run_12559_${psi}.root /users/PAS0654/jflaherty13/source/AraSim/outputs/20231102_SpiceReco/${freqMin}to${freqMax}//deconvolvedWaveforms_run_12559_${psi}.root output/A${station}/20231102_SpiceReco/${freqMin}to${freqMax}/ 1 0 0 0 0

# freqMin=200
# freqMax=250
# python doPolReco_v3.8.py $station 12559 $psi /users/PAS0654/jflaherty13/araAnalysis/araRecoAndSourceSearch/ARA_Reconstruction/brianInterferometry/spicecoreInterfAllRuns/A2/run_012559/recangle_reco_out_run_12559_${psi}.root /users/PAS0654/jflaherty13/source/AraSim/outputs/20231102_SpiceReco/${freqMin}to${freqMax}//deconvolvedWaveforms_run_12559_${psi}.root output/A${station}/20231102_SpiceReco/${freqMin}to${freqMax}/ 1 0 0 0 0

# freqMin=250
# freqMax=300
# python doPolReco_v3.8.py $station 12559 $psi /users/PAS0654/jflaherty13/araAnalysis/araRecoAndSourceSearch/ARA_Reconstruction/brianInterferometry/spicecoreInterfAllRuns/A2/run_012559/recangle_reco_out_run_12559_${psi}.root /users/PAS0654/jflaherty13/source/AraSim/outputs/20231102_SpiceReco/${freqMin}to${freqMax}//deconvolvedWaveforms_run_12559_${psi}.root output/A${station}/20231102_SpiceReco/${freqMin}to${freqMax}/ 1 0 0 0 0

# freqMin=150
# freqMax=300
# python doPolReco_v3.8.py $station 12559 $psi /users/PAS0654/jflaherty13/araAnalysis/araRecoAndSourceSearch/ARA_Reconstruction/brianInterferometry/spicecoreInterfAllRuns/A2/run_012559/recangle_reco_out_run_12559_${psi}.root /users/PAS0654/jflaherty13/source/AraSim/outputs/20231102_SpiceReco/${freqMin}to${freqMax}//deconvolvedWaveforms_run_12559_${psi}.root output/A${station}/20231102_SpiceReco/${freqMin}to${freqMax}/ 1 0 0 0 0

#Pol reco of Alan's birefringence simulations
freqMin=150
freqMax=300
python doPolReco_v3.8.py $station ${psi}0 0 /fs/scratch/PAS0654/SPICE_birefringence_A2_interferometry/recangle_reco_out_run_${psi}0.root /fs/scratch/PAS0654/SPICE_birefringence_A2_deconvolution/${freqMin}to${freqMax}/deconvolvedWaveforms_run_${psi}0.root /fs/scratch/PAS0654/SPICE_birefringence_A2_polReco/ 1 0 0 0 0