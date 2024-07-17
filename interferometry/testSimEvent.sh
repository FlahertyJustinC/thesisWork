#!/bin/bash

index=0
# setupDir=/users/PAS0654/jflaherty13/source/AraSim/outputs/20240505_coherentSumOnICRC
# setupDir=/users/PAS0654/jflaherty13/source/AraSim/outputs/20240510_simplePulserSims_raisedSaturationVoltage
# setupDir=/users/PAS0654/jflaherty13/source/AraSim/outputs/20240513_simplePulserSims_noNoiseWithCorrectArrivalTimes
# setupDir=/users/PAS0654/jflaherty13/source/AraSim/outputs/20240513_simplePulserSims_withCorrectArrivalTimes
setupDir=/users/PAS0654/jflaherty13/source/AraSim/outputs/20240707_testingMcTruthLaunchAngles

./calculatearrivaldir_interf 4 $index $setupDir/AraOut.setup_noBirefringence_A4.txt.run0.root . $setupDir/setup_noBirefringence_A4.txt debug