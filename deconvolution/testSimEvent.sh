#!/bin/bash

# index=45
# # setupDir=/users/PAS0654/jflaherty13/source/AraSim/outputs/20240505_coherentSumOnICRC
# # setupDir=/users/PAS0654/jflaherty13/source/AraSim/outputs/20240510_simplePulserSims_raisedSaturationVoltage
# # setupDir=/users/PAS0654/jflaherty13/source/AraSim/outputs/20240513_simplePulserSims_noNoiseWithCorrectArrivalTimes
# # setupDir=/users/PAS0654/jflaherty13/source/AraSim/outputs/20240513_simplePulserSims_withCorrectArrivalTimes

# setupDir=/users/PAS0654/jflaherty13/source/AraSim/outputs/20240516_simplePulserSims_withCorrectedArrivalTimes_weinerDecon
# # setupDir=/users/PAS0654/jflaherty13/source/AraSim/outputs/20240516_simplePulserSims_noNoiseWithCorrectedArrivalTimes_weinerDecon

# # ./deconvolveWaveform $index 225 300 $setupDir/setup_variablePsi.txt $setupDir/AraSimOutputs/AraOut.setup_variablePsi\=$index.txt.run$index.root $setupDir/interferometry/recangle_reco_out_run_$index.root .

# # ./deconvolveWaveform $index 150 300 $setupDir/setup_variablePsi.txt $setupDir/AraSimOutputs/AraOut.setup_variablePsi\=$index.txt.run$index.root ../interferometry/recangle_reco_out_run_$index.root .

# ./deconvolveWaveform $index 200 300 $setupDir/setup_variablePsi.txt $setupDir/AraSimOutputs/AraOut.setup_variablePsi\=$index.txt.run$index.root $setupDir/interferometry/recangle_reco_out_run_$index.root .  




setupDir=/users/PAS0654/jflaherty13/source/AraSim/outputs/20240516_A24_AlanPulserSims

index=800

rawDirectory=$setupDir/AraSimOutput/
recoDirectory=$setupDir/interferometry/
outDirectory=$setupDir/deconvolution/

#A2 config 6 without and with birefringence
station=4
config=4

# if ! test -f /users/PAS0654/jflaherty13/source/AraSim/outputs/20240321_AlanPulserSims_A24/deconvolution/A${station}/birefringence/deconvolvedWaveforms_run_${index}.root; then

./deconvolveWaveform $index 200 300 ${setupDir}/setup_birefringence_A$station.txt ${rawDirectory}SPICE_birefringence_forJustin/AraOut.setup_birefringence_A$station\=$index.txt.run$index.root ${recoDirectory}/A$station/birefringence/recangle_reco_out_run_$index.root . 