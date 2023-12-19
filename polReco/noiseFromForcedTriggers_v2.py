#Import necessary packages
from ROOT import TCanvas, TGraph
from ROOT import gROOT
import ROOT
import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from ROOT import gInterpreter, gSystem
from ROOT import TChain, TSelector, TTree
import scipy
import sys
import polReco_util as util
import pyrex
import warnings
import itertools
warnings.filterwarnings("ignore")
#import argparse

#add headers from AraSim. Not sure if all of them are needed, and I'm lazy to check that. MAK SURE to change the location of the headers
gInterpreter.ProcessLine('#include "/cvmfs/ara.opensciencegrid.org/trunk/centos7/source/AraSim/Position.h"')
gInterpreter.ProcessLine('#include "/cvmfs/ara.opensciencegrid.org/trunk/centos7/source/AraSim/Report.h"')
gInterpreter.ProcessLine('#include "/cvmfs/ara.opensciencegrid.org/trunk/centos7/source/AraSim/Detector.h"')
gInterpreter.ProcessLine('#include "/cvmfs/ara.opensciencegrid.org/trunk/centos7/source/AraSim/Settings.h"')
gInterpreter.ProcessLine('#include "/cvmfs/ara.opensciencegrid.org/trunk/centos7/source/libRootFftwWrapper/include/FFTtools.h"')

gSystem.Load('/cvmfs/ara.opensciencegrid.org/trunk/centos7/source/AraSim/libAra.so') #load the simulation event library. You might get an error asking for the eventSim dictionry. To solve that, go to where you compiled AraSim, find that file, and copy it to where you set LD_LIBRARY_PATH.
gSystem.Load('/cvmfs/ara.opensciencegrid.org/trunk/centos7/ara_build/lib/libAraEvent.so')
gSystem.Load("/cvmfs/ara.opensciencegrid.org/trunk/centos7/source/libRootFftwWrapper/build/libRootFftwWrapper.so")

#Read console input
if len( sys.argv ) > 3:
    runNumber = int(sys.argv[1])
    subsetNumber = int(sys.argv[2])
    rawSrcFolder = str(sys.argv[3])
    outputFolder = str(sys.argv[4])
else:
    print("Run number is required as well at data subset index \n run as $python doPolReco_interf.py <run-number> <subset-number>")
    sys.exit()

# runNumber = 12559
# subsetNumber = 10
# recoSrcFolder = "../ARA_Reconstruction/brianInterferometry/spicecoreInterfFirstPeakRAbove1050/"
# rawSrcFolder = "../ARA_Reconstruction/data/run_012559/split/"
# outputFolder = "forcedTriggerNoiseV2"

print("runNumber = " + str(runNumber))
print("Sourcing raw data from " + str(rawSrcFolder))

#Reduced event file
rawFilename = rawSrcFolder + "/event0"+str(runNumber)+"__"+str(subsetNumber)+".root"
# recoFilename = recoSrcFolder + "/recangle_reco_out_run_"+str(runNumber)+"_"+str(subsetNumber)+".root"

rawInputFile = ROOT.TFile.Open(rawFilename)
# recoInputFile = ROOT.TFile.Open(recoFilename)

eventTree = rawInputFile.Get("eventTree")
# vertexReco = recoInputFile.Get("vertexReco")

rawEvent = ROOT.RawAtriStationEvent()
eventTree.SetBranchAddress("event",ROOT.AddressOf(rawEvent))

totalRawEvents = eventTree.GetEntries()
print('total raw events:', totalRawEvents)

softTriggerEventList = []
#Import Waveforms and list event numbers by event type (RF, soft trigger, calpulser)
for evt in range(totalRawEvents):
    eventTree.GetEntry(evt)
    usefulEvent = ROOT.UsefulAtriStationEvent(rawEvent,ROOT.AraCalType.kLatestCalib)
    isSoftTrigger = rawEvent.isSoftwareTrigger()
    #Set soft trigger events
    if (isSoftTrigger):
        softTriggerEventList.append(evt)
        
eventList = np.array(softTriggerEventList)
numEvents = len(eventList)

#Initialize arrays that will be saved to .pkl
evt_num = np.zeros(numEvents)
runNumberOut = np.zeros(numEvents)
runSubsetNumberOut = np.zeros(numEvents)
unixtime = np.zeros(numEvents)
timeStamp = np.zeros(numEvents)

sampleWindow = 400 #Gets sampling of 400 datapoints from soft trigger hilbert envelopes
hilbertEnvelopeOut = np.zeros((numEvents,16,sampleWindow))
    
# #Calculate noise envelope and RMS
# noiseEnvelope, noiseRms = util.findMeanNoise(softTriggerEventList, eventTree, rawEvent, ROOT)    
    
for index in range(numEvents):
    evt = eventList[index]
    eventTree.GetEntry(evt)
    usefulEvent = ROOT.UsefulAtriStationEvent(rawEvent,ROOT.AraCalType.kLatestCalib)
    
    unixtime[index] = usefulEvent.unixTime
    evt_num[index] = usefulEvent.eventNumber
    runNumberOut[index] = runNumber
    runSubsetNumberOut[index] = subsetNumber
    
    timeStamp[index] = usefulEvent.timeStamp
    
    hilbertEnvelopeOut[index] = util.getHilbertEnvelopesAllChannels(usefulEvent, sampleWindow=sampleWindow)
    
print("Events processed!")

original_df = pd.DataFrame({"runNumber":runNumberOut, "eventNumber":evt_num, "runSubsetNumber":runSubsetNumberOut.tolist(),
                           "timeStamp":timeStamp.tolist(), "hilbertEnvelope":hilbertEnvelopeOut.tolist(), "unixTime":unixtime.tolist()})

isExist = os.path.exists(outputFolder)
if not isExist:
  # Create a new directory because it does not exist 
  os.makedirs(outputFolder)


original_df.to_pickle(outputFolder + "/forcedTriggerNoise_run%i_%i.pkl"%(runNumber,subsetNumber))
print("Output saved to " + outputFolder + "/forcedTriggerNoise_run%i_%i.pkl"%(runNumber,subsetNumber))
