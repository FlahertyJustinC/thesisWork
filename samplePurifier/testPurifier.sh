#!/bin/bash

rawDirectory=/users/PAS0654/jflaherty13/source/AraSim/outputs/20240417_updatedVertexRecoAnalysis/realDataSpice/
# station=4
# run=6128
# index=723
echo $rawDirectory

#Uses old argument syntax
# ./samplePurifier $station $run ${rawDirectory}A$station/run_00$run/split/event00${run}__$index.root . ${rawDirectory}/setup_KU_antenna_A$station.txt 

# ./samplePurifier 4 6128 ../data/A4/run_006128/event006128.root . ~/source/AraSim/outputs/20240417_updatedVertexRecoAnalysis/realDataSpice/setup_KU_antenna_A4.txt

# ./samplePurifier 2 12559 ../data/A2/run_012559/event012559.root . ~/source/AraSim/outputs/20240417_updatedVertexRecoAnalysis/realDataSpice/setup_KU_antenna_A2.txt

# #Uses new syntax 
# station=4
# run=6128
# index=723
# numChannelsPassed=3
# ./samplePurifier ${run}__$index $numChannelsPassed ${rawDirectory}/setup_KU_antenna_A$station.txt ${rawDirectory}A$station/run_00$run/split/event00${run}__$index.root .

#Test on a full run
station=2
run=12557
numChannelsPassed=3
./samplePurifier ${run} $numChannelsPassed ${rawDirectory}/setup_KU_antenna_A$station.txt ${rawDirectory}A$station/run_0$run/event0${run}.root .