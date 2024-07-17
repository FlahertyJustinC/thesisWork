#!/bin/bash

# ./calculatearrivaldir_interf 4 6119 ../data/A4/run_006119/event006119.root . ~/source/AraSim/outputs/20240417_updatedVertexRecoAnalysis/realDataSpice/setup_KU_antenna_A4.txt
# ./calculatearrivaldir_interf 4 6128 ../data/A4/run_006128/event006128.root . ~/source/AraSim/outputs/20240417_updatedVertexRecoAnalysis/realDataSpice/setup_KU_antenna_A4.txt

# ./calculatearrivaldir_interf 2 12559 ../data/A2/run_012559/event012559.root . ~/source/AraSim/outputs/20240417_updatedVertexRecoAnalysis/realDataSpice/setup_KU_antenna_A2.txt

# ./calculatearrivaldir_interf 2 12559_1 ../data/A2/run_012559/splitPurified/purifiedSample_run_12559_1.root . ~/source/AraSim/outputs/20240417_updatedVertexRecoAnalysis/realDataSpice/setup_KU_antenna_A2.txt debug

run=6137
index=10
station=4
dataDirectory=/users/PAS0654/jflaherty13/source/AraSim/outputs/20240516_A24_SpiceReco/
rawDirectory=$dataDirectory

./calculatearrivaldir_interf $station ${run}_$index ${rawDirectory}A$station/run_00$run/splitPurified/purifiedSample_run_${run}_${index}.root . ${rawDirectory}/setup_KU_antenna_A4.txt debug 0


# run=12576
# index=500
# station=2
# dataDirectory=/users/PAS0654/jflaherty13/source/AraSim/outputs/20240516_A24_SpiceReco/
# rawDirectory=$dataDirectory

# ./calculatearrivaldir_interf $station ${run}_$index ${rawDirectory}A$station/run_0$run/splitPurified/purifiedSample_run_${run}_${index}.root . ${rawDirectory}/setup_KU_antenna_A4.txt debug 3