#!/bin/bash

run=6146
index=0
station=4
dataDirectory=/users/PAS0654/jflaherty13/source/AraSim/outputs/20240516_A24_SpiceReco/
rawDirectory=$dataDirectory

./deconvolveWaveform ${run}_$index 200 300  ${rawDirectory}/setup_KU_antenna_A$station.txt ${rawDirectory}A$station/run_00$run/splitPurified/purifiedSample_run_${run}_$index.root ${rawDirectory}A$station/run_00$run/interferometryPurified/recangle_reco_out_run_${run}_$index.root . debug 0

