#!/bin/bash

setupDir=/users/PAS0654/jflaherty13/source/AraSim/outputs/20240516_A24_AlanPulserSims

index=800

rawDirectory=$setupDir/AraSimOutput/
recoDirectory=$setupDir/interferometry/
outDirectory=$setupDir/deconvolution/

#A2 config 6 without and with birefringence
station=4
config=4


./deconvolveWaveform $index 200 300 ${setupDir}/setup_birefringence_A$station.txt ${rawDirectory}SPICE_birefringence_forJustin/AraOut.setup_birefringence_A$station\=$index.txt.run$index.root ${recoDirectory}/A$station/birefringence/recangle_reco_out_run_$index.root . 
