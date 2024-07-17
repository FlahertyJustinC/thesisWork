#!/bin/bash

# ./deconvolveWaveform 12559 150 300 ~/source/AraSim/outputs/20240417_updatedVertexRecoAnalysis/realDataSpice/setup_KU_antenna_A2.txt ../data/A2/run_012559/event012559.root ../interferometry/recangle_reco_out_run_12559.root . 

# ./deconvolveWaveform 6128 150 300 ~/source/AraSim/outputs/20240417_updatedVertexRecoAnalysis/realDataSpice/setup_KU_antenna_A4.txt ../data/A4/run_006128/event006128.root ../interferometry/recangle_reco_out_run_6128.root . 

run=6146
index=0
station=4
dataDirectory=/users/PAS0654/jflaherty13/source/AraSim/outputs/20240516_A24_SpiceReco/
rawDirectory=$dataDirectory

./deconvolveWaveform ${run}_$index 200 300  ${rawDirectory}/setup_KU_antenna_A$station.txt ${rawDirectory}A$station/run_00$run/splitPurified/purifiedSample_run_${run}_$index.root ${rawDirectory}A$station/run_00$run/interferometryPurified/recangle_reco_out_run_${run}_$index.root . debug 3

# ./deconvolveWaveform $index 200 300  ${rawDirectory}/setup_KU_antenna_A$station.txt ${rawDirectory}A$station/run_00$run/splitPurified/purifiedSample_run_${run}_$index.root ../interferometry/recangle_reco_out_run_${run}_$index.root . debug 3

# run=12576
# index=185
# station=2
# dataDirectory=/users/PAS0654/jflaherty13/source/AraSim/outputs/20240516_A24_SpiceReco/
# rawDirectory=$dataDirectory

# ./deconvolveWaveform $index 200 300  ${rawDirectory}/setup_KU_antenna_A$station.txt ${rawDirectory}A$station/run_0$run/splitPurified/purifiedSample_run_${run}_$index.root ${rawDirectory}A$station/run_0$run/interferometryPurified/recangle_reco_out_run_${run}_$index.root . debug 3

# ./deconvolveWaveform $index 200 300  ${rawDirectory}/setup_KU_antenna_A$station.txt ${rawDirectory}A$station/run_0$run/splitPurified/purifiedSample_run_${run}_$index.root ../interferometry/recangle_reco_out_run_${run}_$index.root . debug 0

