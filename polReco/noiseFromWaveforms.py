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
if len( sys.argv ) > 2:
    runNumber = int(sys.argv[1])
    rawFilename = str(sys.argv[2])
    outputFolder = str(sys.argv[3])
else:
    print("Run number is required \n run as $python noiseFromWaveforms.py <run-number> <rawFileName> <outputFolder>")
    sys.exit()
    
print("runNumber = " + str(runNumber))
print("Sourcing raw data from " + str(rawFilename))

rawInputFile = ROOT.TFile.Open(rawFilename)

file_list = []
file_list.append(rawFilename)

eventTree = TChain("eventTree") #Define chain and tree that needs to be read. "VTree" in this case.
# SimTree = TChain("AraTree2")

for line in file_list:
    eventTree.AddFile(line)
    # SimTree.AddFile(line)

# reportPtr = ROOT.Report()#report pointer
# eventPtr = ROOT.Event()

usefulEvent = ROOT.UsefulAtriStationEvent()
rawEvent = ROOT.RawAtriStationEvent()
eventTree.SetBranchAddress("UsefulAtriStationEvent",ROOT.AddressOf(usefulEvent))
eventTree.SetBranchAddress("RawAtriStationEvent",ROOT.AddressOf(rawEvent))
# SimTree.SetBranchAddress("report",ROOT.AddressOf(reportPtr))
# SimTree.SetBranchAddress("event", ROOT.AddressOf(eventPtr))

# eventTree = rawInputFile.Get("eventTree")
# vertexReco = recoInputFile.Get("vertexReco")

# rawEvent = ROOT.RawAtriStationEvent()
# eventTree.SetBranchAddress("event",ROOT.AddressOf(rawEvent))

totalRawEvents = eventTree.GetEntries()
print('total raw events:', totalRawEvents)

numEvents = totalRawEvents

#Initialize arrays that will be saved to .pkl
evt_num = np.zeros(numEvents)
runNumberOut = np.zeros(numEvents)
runSubsetNumberOut = np.zeros(numEvents)
unixtime = np.zeros(numEvents)
timeStamp = np.zeros(numEvents)

sampleWindow = 100 #Gets sampling of 100 from beginning of waveform for noise.
hilbertEnvelopeOut = np.zeros((numEvents,16,sampleWindow))

for evt in range(numEvents):
    eventTree.GetEntry(evt)
    
    unixtime[evt] = usefulEvent.unixTime
    evt_num[evt] = usefulEvent.eventNumber
    runNumberOut[evt] = runNumber
    
    timeStamp[evt] = usefulEvent.timeStamp
    
    hilbertEnvelopeOut[evt] = util.getHilbertEnvelopesAllChannels(usefulEvent, sampleWindow=sampleWindow)
    
print("Events processed!")

#Create Histogram plotting functions for all quantities that also implement an SNR cut
def makeHistogramByChannelPair(array, addMask=0, snrCutoff=5, pol='', xLims=0, numBins=0, binWidth=1, xLabel='', titleLabel='', fileLabel='', depthCut=0, fit=None, degreeMinCut=0, degreeMaxCut=0):
    # print("Array length pre degree cut: " + str(len(array[mask])))
    #sns.set(rc = {'figure.figsize':(18,18)},font_scale = 2)
    font = {'family' : 'normal',
            'weight' : 'bold',
            'size'   : 16}

    plt.rc('font', **font)
    fig, ax = plt.subplots(nrows=2,ncols=4, figsize=(18,9),sharex=True, sharey=True)
    ax = ax.ravel()
    
    extraFileLabel = ''
    extraTitleLabel = ''
    # if (degreeMinCut != 0):
    #     #extraFileLabel += '-degreeMin=' + str(degreeMinCut)
    #     #extraTitleLabel +=  ', Angle $>$ ' +  str(degreeMinCut)
    #     mask *= array>degreeMinCut
    # if (degreeMaxCut != 0):
    #     #extraFileLabel += '-degreeMax=' + str(degreeMaxCut)
    #     print("using degree max cut.")
    #     mask *= array<degreeMaxCut
    # print("Array length post cut: " + str(len(array[mask])))
    maxBinOut = np.empty(8)
    for pair in range(8):
        yData, xData, figure = ax[pair].hist(array[:,pair], bins = numBins, density=True)
        maxBin = np.median(array[:,pair])
        #ax[pair].axvline(maxBin, label=str(maxBin), color='black', linestyle='--')
        #ax[pair].axvline(maxBin, label="Median = %5.g"%(maxBin), color='black', linestyle='--')
        ax[pair].axvline(maxBin, label="Median = {:.2e}".format(maxBin), color='black', linestyle='--')
        if pol in ['vpol', 'VPol', 'Vpol', 'v', 'V']:
            ax[pair].set_title("Channel " + str(pair))# + '/' + str(pair+8))
        elif pol in ['hpol', 'HPol', 'Hpol', 'h', 'H']:
            ax[pair].set_title("Channel " + str(pair+8))# + '/' + str(pair+8))
        else:
            ax[pair].set_title("Channel " + str(pair) + '/' + str(pair+8))
        ax[pair].legend()
        maxBinOut[pair] = maxBin
    fig.suptitle(titleLabel + extraTitleLabel + ", N = " + str(int(len(array)/400)))
    fig.supylabel("Number of Events")
    fig.supxlabel(xLabel)
    #plt.xlabel(xLabel)
    if (xLims != 0):
        ax[0].set_xlim(xLims)
    print(maxBinOut)
    return maxBinOut
    #plt.ylabel("Number of Events")
    #plt.legend()

    #plt.savefig('./'+plotFolder+'/snrCutoff='+str(snrCutoff)+'/'+fileLabel+'Hist-snrCutoff='+str(snrCutoff)+extraFileLabel+'-PerChannel.png')
    
    
#Flatten hilbert envelopes along channels

hilbertEnvelope = hilbertEnvelopeOut

totalEvents = numEvents
hilbertFlattenedVPol = np.zeros((totalEvents*sampleWindow,8))
hilbertFlattenedHPol = np.zeros((totalEvents*sampleWindow,8))
for ch in range(8):
    hilbertFlattenedVPol[:,ch] = hilbertEnvelope[:,ch,:].flatten()
    hilbertFlattenedHPol[:,ch-8] = hilbertEnvelope[:,ch+8,:].flatten()
    
    
powerVNoiseMedians = makeHistogramByChannelPair(hilbertFlattenedVPol, xLims=[0,200], numBins=100, pol='v', titleLabel="Hilbert Envelopes of VPol Noise of Soft Trigger Events (400 Values Per Event, Covering 200 ns)", xLabel="Hilbert Envelope of Waveform [mV?]")

powerHNoiseMedians = makeHistogramByChannelPair(hilbertFlattenedHPol, xLims=[0,200], numBins=100, pol='h', titleLabel="Hilbert Envelope Values of HPol Noise of Soft Trigger Events (400 Values Per Event, Covering 200 ns)", xLabel="Hilbert Envelope of Waveform [mV?]")

powerNoiseMedians = np.concatenate((powerVNoiseMedians,powerHNoiseMedians))
    


original_df = pd.DataFrame({"runNumber":runNumberOut, "eventNumber":evt_num,
                           "timeStamp":timeStamp.tolist(), "hilbertEnvelope":hilbertEnvelopeOut.tolist(), "unixTime":unixtime.tolist()})



isExist = os.path.exists(outputFolder)
if not isExist:
  # Create a new directory because it does not exist 
  os.makedirs(outputFolder)


original_df.to_pickle(outputFolder + "/noiseFromWaveforms_run%i.pkl"%(runNumber))
print("Output saved to " + outputFolder + "/noiseFromWaveforms_run%i.pkl"%(runNumber))

# original_df = pd.DataFrame({"powerNoiseConstant":powerNoiseMedians})
original_df.to_pickle(outputFolder + "/mediansOfWaveformNoise_run%i.pkl"%(runNumber))

print("Output saved to " + outputFolder + "/mediansOfWaveformNoise_run%i.pkl"%(runNumber))