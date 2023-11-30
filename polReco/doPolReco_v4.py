"""
#####doPolReco.py#####

This script is part of a series of steps in my polarization reconstruction of ARA signals. The steps are:
    1. Vertex reconstruction via interferometry/calculatearrivaldir_interf.cpp
    2. Deconvolution of signal via deconvolvution/deconvolveWaveform.cpp
    3. Reconstruction of polarization direction via polReco/doPolReco.py
    4. Plotting of reconstruction resolution via polReco/plotRecoAngles.ipynb
    
Script is still under construction, with a planned outline below:
    
Doing an overhaul of my polReco code now that deconvolution is separate.  Here's a list of things I need when running it, so I can set up proper keyword arguments:  

- station number (for looking up geometry? It should just be able ot grab it from the input file)
- dataType (whether it's real or simulated data. Maybe I can have that encoded in the deconvolution file? Like if it has a AraTree?  Yeah I should do that, as raw data will have that feature so I can just check for the AraTree.)
- run number (I should fold in the subset number for when I split data, maybe just in the filepath)
- filepath to raw or deconvolved data
- filepath to vertex reco data
- filepath for output
- keyword argument for choosing channels?  Allows for doing reco on just top or bottom antennas.
- keyword argument for deltaT tolerance between H and V peaks?  Currently defaults to 10 ns.  I'll need to use this for keyword arguments:  https://docs.python.org/dev/library/argparse.html

Things I should definitely remove from my polReco script:

- Deconvolution using pyrex
- noise subtraction in time-domain (doing frequency domain in the deconvolution now)
- gain balancing (handled by deconvolution now).

Things I should add:
- Coherent summing of waveforms. 
- Should save the input parameters, like tolerance, to the output file for easy reference. Kinda like how AraSim outputs save the values of the setup file.

Author: Justin Flaherty
Date: November 26, 2023
"""

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
from datetime import datetime
import argparse
warnings.filterwarnings("ignore")

gInterpreter.ProcessLine('#include "/users/PAS0654/jflaherty13/source/AraSim/Position.h"')
gInterpreter.ProcessLine('#include "/users/PAS0654/jflaherty13/source/AraSim/Report.h"')
gInterpreter.ProcessLine('#include "/users/PAS0654/jflaherty13/source/AraSim/Detector.h"')
gInterpreter.ProcessLine('#include "/users/PAS0654/jflaherty13/source/AraSim/Settings.h"')
gInterpreter.ProcessLine('#include "/cvmfs/ara.opensciencegrid.org/trunk/centos7/source/libRootFftwWrapper/include/FFTtools.h"')

gSystem.Load('/cvmfs/ara.opensciencegrid.org/trunk/centos7/ara_build/lib/libAraEvent.so')
gSystem.Load("/cvmfs/ara.opensciencegrid.org/trunk/centos7/source/libRootFftwWrapper/build/libRootFftwWrapper.so")

#Read console input
if len( sys.argv ) > 3:
    runNumber = int(sys.argv[1])
    recoSrcFolder = str(sys.argv[2])
    rawSrcFolder = str(sys.argv[3])
    outputFolder = str(sys.argv[4])

else:
    print("ERROR:  Missing necessary flags. \nRun as:  python doPolReco_interf.py <run-number> <subset-number> <reco source folder> <raw source folder> <real or MC flag> <RF, soft-trigger, or calpulser flag> <noise type flag> <gainBalance flag> <deconvolution flag>")
    sys.exit()

#Initialize flag values before argparser
subset = ''


print("runNumber = " + str(runNumber))
print("Sourcing vertex reconstruction from " + str(recoSrcFolder))
print("Sourcing raw data from " + str(rawSrcFolder))

rawFilename = rawSrcFolder
recoFilename = recoSrcFolder


rawInputFile = ROOT.TFile.Open(rawFilename)
recoInputFile = ROOT.TFile.Open(recoFilename)

try:

    eventTree = rawInputFile.Get("eventTree")
    vertexReco = recoInputFile.Get("vertexReco")

    rawEvent = ROOT.RawAtriStationEvent()
    eventTree.SetBranchAddress("event",ROOT.AddressOf(rawEvent))

    totalRawEvents = eventTree.GetEntries()
    print("Using dataLike condition")
    print('total raw events:', totalRawEvents)

    totalRecoEvents = vertexReco.GetEntries()
    print('total reco events:', totalRecoEvents)
    dataLike = True

except TypeError:
    
    print("Using simLike condition")
    
    file_list = []
    file_list.append(rawFilename)
    
    eventTree = TChain("eventTree") #Define chain and tree that needs to be read. "VTree" in this case.
    SimTree = TChain("AraTree2")
    vertexReco = recoInputFile.Get("vertexReco")

    for line in file_list:
        eventTree.AddFile(line)
        SimTree.AddFile(line)
        
    reportPtr = ROOT.Report()#report pointer
    eventPtr = ROOT.Event()

    usefulEvent = ROOT.UsefulAtriStationEvent()
    rawEvent = ROOT.RawAtriStationEvent()
    eventTree.SetBranchAddress("UsefulAtriStationEvent",ROOT.AddressOf(usefulEvent))
    eventTree.SetBranchAddress("RawAtriStationEvent",ROOT.AddressOf(rawEvent))
    # SimTree.SetBranchAddress("report",ROOT.AddressOf(reportPtr))
    # SimTree.SetBranchAddress("event", ROOT.AddressOf(eventPtr))

    totalRawEvents = eventTree.GetEntries()
    print('total raw events:', totalRawEvents) 
    
    totalRecoEvents = vertexReco.GetEntries()
    print('total reco events:', totalRecoEvents)
    dataLike = False
    
    
if (totalRecoEvents != totalRawEvents):
    print("ERROR: Raw event file and reconstructed event file have different number of events! Ensure that you are using compatible files!")
    sys.exit
    
#Obtain list of RF events, soft-trigger events, and calpulser events.
rfEventList = []
softTriggerEventList = []
calpulserEventList = []

#Import Waveforms and list event numbers by event type (RF, soft trigger, calpulser)
if (dataLike):
    for evt in range(totalRawEvents):
        vertexReco.GetEntry(evt)
        eventTree.GetEntry(evt)
        usefulEvent = ROOT.UsefulAtriStationEvent(rawEvent,ROOT.AraCalType.kLatestCalib)
        if (vertexReco.eventNumber != usefulEvent.eventNumber):
            print("ERROR:  Mismatch of event numbers!  Check your input files and ensure that they are compatible!")
            sys.exit()
        else:
            isCalpulser = rawEvent.isCalpulserEvent()
            isSoftTrigger = rawEvent.isSoftwareTrigger()
            timeStampTemp = usefulEvent.timeStamp
            unixTimeTemp = usefulEvent.unixTime
            # fittedTimeStamp = timeStampFit(unixTimeTemp)
            #Set soft trigger events
            if (isSoftTrigger):
                softTriggerEventList.append(evt)
            #Set calpulser events
            elif (isCalpulser):
                calpulserEventList.append(evt)
            else:
                rfEventList.append(evt)   
            
if (not dataLike):
    for evt in range(totalRawEvents):
        # vertexReco.GetEntry(evt)
        eventTree.GetEntry(evt)
        isCalpulser = usefulEvent.isCalpulserEvent()
        isSoftTrigger = usefulEvent.isSoftwareTrigger()
        # Set soft trigger events
        if (isSoftTrigger):
            print("Event " + str(evt) + " is a soft trigger")
            softTriggerEventList.append(evt)
        #Set calpulser events
        elif (isCalpulser):
            print("Event " + str(evt) + " is a calpulser")
            calpulserEventList.append(evt)
        else:
            print("Event " + str(evt) + " is a RF event.")
            rfEventList.append(evt)
        #Some AraSim results set ALL events as soft triggers, so this covers that condition.
        if (len(rfEventList) == 0):
            rfEventList = softTriggerEventList      
        
#Convert event lists into arrays    
rfEventList = np.array(rfEventList)
softTriggerEventList = np.array(softTriggerEventList)
calpulserEventList = np.array(calpulserEventList)

#Find number of each event type
numRfEvents = len(rfEventList)
numSoftTriggerEvents = len(softTriggerEventList)
numCalpulserEvents = len(calpulserEventList)

#Using event type flag set by user, set which event type will be analyzed (RF, soft trigger, or calpulser)
eventList = rfEventList[:]
solution = 'direct'
numEvents = len(eventList)  

#Initialize arrays that will be saved to .pkl
evt_num = np.zeros(numEvents)
runNumberOut = np.zeros(numEvents)
runSubsetNumberOut = np.zeros(numEvents)
thetaRecoOut = np.zeros((numEvents, 16))
thetaVertexOut = np.zeros(numEvents)
phiRecoOut = np.zeros((numEvents, 16))
phiVertexOut = np.zeros(numEvents)
whichSol = np.zeros(numEvents)
recoROut = np.zeros((numEvents, 8))
psiRecoOut = np.zeros((numEvents, 8))
unixtime = np.zeros(numEvents)
timeStamp = np.zeros(numEvents)
pulserDepth = np.zeros(numEvents)
    
#Define arrays
snrsOut = np.zeros((numEvents,16))
vSnrOut = np.zeros(numEvents)
hSnrOut = np.zeros(numEvents)
hilbertPeakOut = np.zeros((numEvents,16))
peakTimeOut = np.zeros((numEvents,16))

#Initialize array for timing cut values.
runEventNumber = np.empty(numEvents).astype(str)   

noiseEnvelope, noiseRms = util.findMeanNoiseFromWaveform(eventList, eventTree, usefulEvent, ROOT)

print("Setting noise envelope to None.  Fix this later Justin.")
noiseEnvelope = None

#Import unixtime of first event and figure out the date, then use that date to calculate the pulser depth
evt = eventList[0]
vertexReco.GetEntry(evt)
eventTree.GetEntry(evt)
try:
    unixtimeTemp = usefulEvent.unixTime
    date = datetime.utcfromtimestamp(unixtimeTemp)
    #Convert date object from UTC to NZ local time
    import pytz
    from datetime import datetime, timezone
    nz = pytz.timezone('NZ')
    utc = pytz.timezone('UTC')
    def utc_to_local(utc_dt):
        return utc_dt.replace(tzinfo=timezone.utc).astimezone(tz=nz)

    date = utc_to_local(date)
    year = date.year
    month = date.month
    day = date.day

    pulserDepthFile = "spicePulser_"+str(month)+str(day)+"Depth.txt"

    depthFile = pd.read_csv(pulserDepthFile)
    print("Sourcing pulser depth from " + "./"+outputFolder+"/"+pulserDepthFile) #Debugging JCF 9/14/2022
    time = pd.to_datetime(depthFile.NZ_Time)
    time.head()
    newTime = time.apply(lambda dt: dt.replace(day=day, month = month, year = year))
    # newTime#Still in NZ local time. Need to translate to UTC
    df = pd.DataFrame(1, index=newTime, columns=['X'])
    import pytz
    nz = pytz.timezone('NZ')
    utc = pytz.timezone('UTC')
    df.index = df.index.tz_localize(nz).tz_convert(utc)
    unixTimeDepth = (df.index - pd.Timestamp("1970-01-01").tz_localize(utc)) // pd.Timedelta('1s')#This is unix time 
    f = scipy.interpolate.interp1d(unixTimeDepth, depthFile.depth,bounds_error=False, fill_value=0.)
    print("Min Pulser Unix Time = " + str(unixTimeDepth.min()))
    print("Max Pulser Unix Time = " + str(unixTimeDepth.max()))
    pulserDepthCalculation = True
except OverflowError:
    print("Unix time of event undefined.  Bypassing pulser depth calculation.")
    pulserDepthCalculation = False
    
    
# print("eee")
for index in range(numEvents):
    evt = eventList[index]
    print("evt = " + str(evt))
    vertexReco.GetEntry(evt)
    eventTree.GetEntry(evt)
    
    evt_num[index] = usefulEvent.eventNumber
    runNumberOut[index] = runNumber
    runEventNumber[index] = str(int(runNumber)) + "_" + str(int(evt_num[index]))
    thetaVertexOut[index] = 90 - vertexReco.bestTheta_out
    phiVertexOut[index] = vertexReco.bestPhi_out % 360
    
    thetaRecoOut[index] = np.degrees(np.array(vertexReco.reco_arrivalThetas_out)) % 180
    phiRecoOut[index] = np.degrees(np.array(vertexReco.reco_arrivalPhis_out)) % 360
    
    timeStamp[index] = usefulEvent.timeStamp
    
    hilbertPeakOut[index], peakTimeOut[index], snrsOut[index] = util.peakHilbert(usefulEvent, vertexReco, noiseEnvelope, noiseRms, solution="single", timeShift=14.1)
    recoROut[index], psiRecoOut[index] = util.calculatePsiAndR(hilbertPeakOut[index]**2)
    
    ##Adding line to export snrs from vertexReco
    snrsOut[index] = np.array(vertexReco.snrs_out)
    vSnrOut[index] = vertexReco.v_snr_out
    hSnrOut[index] = vertexReco.h_snr_out
    
    if (pulserDepthCalculation):
        pulserDepth[index] = f(usefulEvent.unixTime)
    else:
        # print("Pulser depth calculation bypassed.  Setting run number to pulser depth.")
        pulserDepth[index] = runNumber
    unixtime[index] = usefulEvent.unixTime 
    
    
deltaTPeakOut = peakTimeOut[:,8:] - peakTimeOut[:,:8] 

print("Events processed!")  

#Save data to pandas file
original_df = pd.DataFrame({"runNumber":runNumberOut, 
                            "eventNumber":evt_num.tolist(),  
                            "timeStamp":timeStamp.tolist(), 
                            "runEventNumber":runEventNumber.tolist(),
                            "thetaReco":thetaRecoOut.tolist(), 
                            "phiReco":phiRecoOut.tolist(), 
                            "thetaVertex":thetaVertexOut.tolist(), 
                            "phiVertex":phiVertexOut.tolist(), 
                            "SNR":snrsOut.tolist(), 
                            "vSNR":vSnrOut.tolist(), 
                            "hSNR":hSnrOut.tolist(), 
                            "whichSol":whichSol.tolist(), 
                            "unixtime":unixtime.tolist(), 
                            "psiReco":psiRecoOut.tolist(), 
                            "pulserDepth":pulserDepth.tolist(),
                            "hilbertPeaks":hilbertPeakOut.tolist(), 
                            "peakTimes":peakTimeOut.tolist(), 
                            "deltaTPeaks":deltaTPeakOut.tolist()
                           })
 

isExist = os.path.exists(outputFolder)
if not isExist:
  # Create a new directory because it does not exist 
  os.makedirs(outputFolder)

outfilePath = outputFolder + "/polReco_run%i.pkl"%(runNumber)  

original_df.to_pickle(outfilePath)
print("Output saved to " + outfilePath)