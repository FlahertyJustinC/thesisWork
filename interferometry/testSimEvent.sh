#!/bin/bash

index=45
# setupDir=/users/PAS0654/jflaherty13/source/AraSim/outputs/20240505_coherentSumOnICRC
# setupDir=/users/PAS0654/jflaherty13/source/AraSim/outputs/20240510_simplePulserSims_raisedSaturationVoltage
setupDir=/users/PAS0654/jflaherty13/source/AraSim/outputs/20240513_simplePulserSims_noNoiseWithCorrectArrivalTimes

./calculatearrivaldir_interf 2 $index $setupDir/AraSimOutputs/AraOut.setup_variablePsi\=$index.txt.run$index.root . $setupDir/setup_variablePsi.txt 