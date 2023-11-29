# #####doPolReco.py#####

# This script is part of a series of steps in my polarization reconstruction of ARA signals. The steps are:
#     1. Vertex reconstruction via interferometry/calculatearrivaldir_interf.cpp
#     2. Deconvolution of signal via deconvolvution/deconvolveWaveform.cpp
#     3. Reconstruction of polarization direction via polReco/doPolReco.py
#     4. Plotting of reconstruction resolution via polReco/plotRecoAngles.ipynb
    
# Script is still under construction, with a planned outline below:
    
# Doing an overhaul of my polReco code now that deconvolution is separate.  Here's a list of things I need when running it, so I can set up proper keyword arguments:  

# - station number (for looking up geometry? It should just be able ot grab it from the input file)
# - dataType (whether it's real or simulated data. Maybe I can have that encoded in the deconvolution file? Like if it has a AraTree?  Yeah I should do that, as raw data will have that feature so I can just check for the AraTree.)
# - run number (I should fold in the subset number for when I split data, maybe just in the filepath)
# - filepath to raw or deconvolved data
# - filepath to vertex reco data
# - filepath for output
# - keyword argument for choosing channels?  Allows for doing reco on just top or bottom antennas.
# - keyword argument for deltaT tolerance between H and V peaks?  Currently defaults to 10 ns.  I'll need to use this for keyword arguments:  https://docs.python.org/dev/library/argparse.html

# Things I should definitely remove from my polReco script:

# - Deconvolution using pyrex
# - noise subtraction in time-domain (doing frequency domain in the deconvolution now)
# - gain balancing (handled by deconvolution now).

# Things I should add:
# - Coherent summing of waveforms. 
# - Should save the input parameters, like tolerance, to the output file for easy reference. Kinda like how AraSim outputs save the values of the setup file.

# Author: Justin Flaherty
# Date: November 26, 2023

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
warnings.filterwarnings("ignore")

gInterpreter.ProcessLine('#include "/users/PAS0654/jflaherty13/source/AraSim/Position.h"')
gInterpreter.ProcessLine('#include "/users/PAS0654/jflaherty13/source/AraSim/Report.h"')
gInterpreter.ProcessLine('#include "/users/PAS0654/jflaherty13/source/AraSim/Detector.h"')
gInterpreter.ProcessLine('#include "/users/PAS0654/jflaherty13/source/AraSim/Settings.h"')
gInterpreter.ProcessLine('#include "/cvmfs/ara.opensciencegrid.org/trunk/centos7/source/libRootFftwWrapper/include/FFTtools.h"')

gSystem.Load('/cvmfs/ara.opensciencegrid.org/trunk/centos7/ara_build/lib/libAraEvent.so')
gSystem.Load("/cvmfs/ara.opensciencegrid.org/trunk/centos7/source/libRootFftwWrapper/build/libRootFftwWrapper.so")

#Read console input
if len( sys.argv ) > 10:
    stationID = int(sys.argv[1])
    runNumber = int(sys.argv[2])
    subsetNumber = int(sys.argv[3])
    recoSrcFolder = str(sys.argv[4])
    rawSrcFolder = str(sys.argv[5])
    outputFolder = str(sys.argv[6])
    dataTypeFlag = int(sys.argv[7])
    signalTypeFlag = int(sys.argv[8])
    noiseTypeFlag = int(sys.argv[9])
    gainBalance = int(sys.argv[10])
    deconvolution = int(sys.argv[11])
    # softTriggerHpol = int(sys.argv[10])

else:
    print("ERROR:  Missing necessary flags. \nRun as:  python doPolReco_interf.py <run-number> <subset-number> <reco source folder> <raw source folder> <real or MC flag> <RF, soft-trigger, or calpulser flag> <noise type flag> <gainBalance flag> <deconvolution flag>")
    sys.exit()

tolerance=10
dt=0
testChannels = [0,2,4]

print("runNumber = " + str(runNumber))
print("Sourcing vertex reconstruction from " + str(recoSrcFolder))
print("Sourcing raw data from " + str(rawSrcFolder))

rawFilename = rawSrcFolder
recoFilename = recoSrcFolder

if (dataTypeFlag == 0):
    rawInputFile = ROOT.TFile.Open(rawFilename)
    recoInputFile = ROOT.TFile.Open(recoFilename)

    eventTree = rawInputFile.Get("eventTree")
    vertexReco = recoInputFile.Get("vertexReco")

    rawEvent = ROOT.RawAtriStationEvent()
    eventTree.SetBranchAddress("event",ROOT.AddressOf(rawEvent))

    totalRawEvents = eventTree.GetEntries()
    print('total raw events:', totalRawEvents)

    totalRecoEvents = vertexReco.GetEntries()
    print('total reco events:', totalRecoEvents)
    
elif ((dataTypeFlag==1) or (dataTypeFlag==2)):
    rawInputFile = ROOT.TFile.Open(rawFilename)
    recoInputFile = ROOT.TFile.Open(recoFilename)
    
    # eventTree = rawInputFile.Get("eventTree")
    # simTree = rawInputFile.Get("AraTree2")
    
    file_list = []
    file_list.append(rawFilename)
    
    eventTree = TChain("eventTree") #Define chain and tree that needs to be read. "VTree" in this case.
    # SimTree = TChain("AraTree2")
    vertexReco = recoInputFile.Get("vertexReco")

    for line in file_list:
        eventTree.AddFile(line)
        # SimTree.AddFile(line)
        
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
    
if (totalRecoEvents != totalRawEvents):
    print("ERROR: Raw event file and reconstructed event file have different number of events! Ensure that you are using compatible files!")
    sys.exit
    
#TODO:  Need to phase out this section of the code, as it's no longer relevant.
#Import noise
if (noiseTypeFlag == 0):
    print("Using no-noise model")
    powerNoiseConstant = None

elif (noiseTypeFlag == 1):
    noisePath = "./noiseFromWaveforms/mediansOfWaveformNoise_run%i.pkl"%(runNumber)
    print("Using waveform noise model.")
    print("Sourcing noise from " + noisePath)
    name = os.path.join(noisePath)
    data = pd.read_pickle(name)
    powerNoiseConstant = np.empty(16)
    for i,entry in zip(range(0,len(data)), data.itertuples()):
        powerNoiseConstant[i] = entry.powerNoiseConstant
else:
    print("No noise model selected.  Setting to None.")
    powerNoiseConstant = None
    
    
    
#Obtain list of RF events, soft-trigger events, and calpulser events.
rfEventList = []
softTriggerEventList = []
calpulserEventList = []

#Import Waveforms and list event numbers by event type (RF, soft trigger, calpulser)
if (dataTypeFlag == 0):
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
            
if (dataTypeFlag == 1):  #Swapping this with condition below because things are breaking.
    #TODO: Make commented code below work with simulated events.
    for evt in range(totalRawEvents):
        # vertexReco.GetEntry(evt)
        eventTree.GetEntry(evt)
        isCalpulser = usefulEvent.isCalpulserEvent()
        isSoftTrigger = usefulEvent.isSoftwareTrigger()
        print("isSoftTrigger = " +str(isSoftTrigger))
        rfEventList.append(evt) #All AraSim events are softtrigger events, so just populate with all events for now.
        # # Set soft trigger events
        # if (isSoftTrigger):
        #     print("Event " + str(evt) + " is a soft trigger")
        #     softTriggerEventList.append(evt)
        # #Set calpulser events
        # elif (isCalpulser):
        #     print("Event " + str(evt) + " is a calpulser")
        #     calpulserEventList.append(evt)
        # else:
        #     print("Event " + str(evt) + " is a RF event.")
        #     rfEventList.append(evt)
        # if (len(rfEventList) == 0):
        #     rfEventList = softTriggerEventList
        
if (dataTypeFlag == 2):
    #TODO: Make commented code below work with simulated events.
    for evt in range(totalRawEvents):
        # vertexReco.GetEntry(evt)
        eventTree.GetEntry(evt)
        isCalpulser = usefulEvent.isCalpulserEvent()
        isSoftTrigger = usefulEvent.isSoftwareTrigger()
        # rfEventList.append(evt) #All AraSim events are softtrigger events, so just populate with all events for now.
        #Set soft trigger events
        if (isSoftTrigger):
            softTriggerEventList.append(evt)
        #Set calpulser events
        elif (isCalpulser):
            calpulserEventList.append(evt)
        else:
            rfEventList.append(evt)           
        
#Convert event lists into arrays    
rfEventList = np.array(rfEventList)
softTriggerEventList = np.array(softTriggerEventList)
calpulserEventList = np.array(calpulserEventList)

#Find number of each event type
numRfEvents = len(rfEventList)
numSoftTriggerEvents = len(softTriggerEventList)
numCalpulserEvents = len(calpulserEventList)

#Using event type flag set by user, set which event type will be analyzed (RF, soft trigger, or calpulser)
if (signalTypeFlag == 0):
    eventList = rfEventList[:]
    solution = 'direct'
    numEvents = len(eventList)  
    # if (len(rfEventList) > len(softTriggerEventList)):
    #     numEvents = len(softTriggerEventList)
    # else:
    #     numEvents = len(eventList)    
elif (signalTypeFlag == 1):
    solution = "noise"
    eventList = softTriggerEventList[:]
    numEvents = len(eventList)
elif (signalTypeFlag == 2):
    solution = 'calpulser'
    eventList = calpulserEventList[:]
    numEvents = len(eventList)
else:
    print("Unknown event type flag.  Please choose 0, 1, or 2 for RF Events, Soft-Trigger Events, or Calpulser Events.")
    sys.exit()        
    
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
recoRSoftTriggerHpolOut = np.zeros((numEvents, 8))
psiRecoSoftTriggerHpolOut = np.zeros((numEvents, 8))
recoRadiusOut = np.zeros(numEvents)
unixtime = np.zeros(numEvents)
timeStamp = np.zeros(numEvents)
pulserDepth = np.zeros(numEvents)  #TODO:  Need to include boolean for getting pulserDepth.
    
#Define arrays to house the powerV and powerH values.
powerOut = np.zeros((numEvents,16))
powerSoftTriggerHpolOut = np.zeros((numEvents,16))
powerNoiseFromSoftTriggerOut = np.zeros((numEvents,16))
snrsOut = np.zeros((numEvents,16))
vSnrOut = np.zeros(numEvents)
hSnrOut = np.zeros(numEvents)
hilbertPeakOut = np.zeros((numEvents,16))
hilbertPeakSoftTriggerHpolOut = np.zeros((numEvents,16))
peakTimeOut = np.zeros((numEvents,16))
# deltaTPeakOut = np.zeros((numEvents,8))

#Initialize array for timing cut values.
timingChannels = np.array([0,2,4])
timingChannelsOut = np.zeros((numEvents,3))
timingChannelsOut[:] = timingChannels
timingCutPass = np.zeros(numEvents)
timingDifferenceMeasured = np.zeros((numEvents,3))
timingDifferenceExpected = np.zeros((numEvents,3))
runEventNumber = np.empty(numEvents).astype(str)   

#Commenting this out because it's no longer used.
#Obtain antenna coordinates
geomTool = ROOT.AraGeomTool.Instance()
# antennaPositions = np.zeros((16,3))
# for ch in range(16):
#     for coord in range(3):
#         antennaPositions[ch, coord] = geomTool.getStationInfo(stationID).getAntennaInfo(ch).antLocation[coord]

#Calculate noise envelope and RMS
#TODO: Add boolean flag for getting noise from soft trigers or event waveforms
noiseEnvelope, noiseRms = util.findMeanNoise(softTriggerEventList, eventTree, rawEvent, ROOT, waveformSampleNs = 50)

#Hardcoding in noise RMS from distribution of deconvolved pulser waveforms - JCF 7/27/2023
print("USING HARDCODED NOISE TO CALCULATE DECONVOLVED SNR.  FIX THIS FOR FUTURE")
noiseRms = np.array([0.097771  , 0.05745579, 0.06394957, 0.08120066, 0.06068223,
       0.11035589, 0.06460599, 0.07275382, 0.05654482, 0.09287774,
       0.0919257 , 0.08902603, 0.07838613, 0.12511595, 0.08134452,
       0.03727756])
# noiseEnvelope, noiseRms = util.findMeanNoiseFromWaveform(eventList, eventTree, usefulEvent, ROOT)

#Debugging and trying to use the median noise for gainCorrection as the noise envelope. - JCF 6/16/2023
noiseEnvelope = powerNoiseConstant

if ((dataTypeFlag == 0) or (dataTypeFlag == 1)):
# if (True):  #Forcing this to be true for the spiceCore reconstruction using my new deconvolution script.
    #Obtain Spicecore Coordinates
    a2EastNorth = np.array([35481,45369])*0.3048 #Convert feet to meters.
    spiceEastNorth = np.array([42359,48974])*0.3048 #Convert feet to meters

    #Calculate Spice position relative to A2.
    spiceInA2Frame = spiceEastNorth - a2EastNorth




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
    except OverflowError:
        print("Unix time of event undefined.  Bypassing pulser depth calculation.")
    
    
for index in range(numEvents):
    evt = eventList[index]
    print("evt = " + str(evt))
    vertexReco.GetEntry(evt)
    eventTree.GetEntry(evt)
    powerOut[index], snrsOut[index] = util.powerFromWaveformSubtractHilbertNoise(rawEvent, usefulEvent, vertexReco, ROOT, noiseEnvelope, noiseRms, gainBalance = gainBalance, gainCorrection = powerNoiseConstant)
    
    evt_num[index] = usefulEvent.eventNumber
    runNumberOut[index] = runNumber
    # runSubsetNumberOut[index] = subsetNumber
    runEventNumber[index] = str(int(runNumber)) + "_" + str(int(evt_num[index]))
    thetaVertexOut[index] = 90 - vertexReco.bestTheta_out
    phiVertexOut[index] = vertexReco.bestPhi_out % 360
    
    thetaRecoOut[index] = np.degrees(np.array(vertexReco.reco_arrivalThetas_out)) % 180
    phiRecoOut[index] = np.degrees(np.array(vertexReco.reco_arrivalPhis_out)) % 360
    
    timeStamp[index] = usefulEvent.timeStamp
    
    hilbertPeakOut[index], peakTimeOut[index], snrsOut[index] = util.peakHilbert(usefulEvent, vertexReco, noiseEnvelope, noiseRms, gainBalance=gainBalance, gainCorrection=powerNoiseConstant, deconvolution=deconvolution, tolerance=tolerance, solution="single", timeShift=14.1)
    recoROut[index], psiRecoOut[index] = util.calculatePsiAndR(hilbertPeakOut[index]**2)
    
    ##Adding line to export snrs from vertexReco
    # snrsOut[index] = np.array(vertexReco.snrs_out)
    # vSnrOut[index] = vertexReco.v_snr_out
    # hSnrOut[index] = vertexReco.h_snr_out
    
    ##Try SNR from Hilbert envelope:
    ##Nope can't do this with the deconvolved files, as the noise is basically zero after deconvolution.  This is giving me infinite SNR.
    vSnrOut[index] = np.median(snrsOut[index,:8])
    hSnrOut[index] = np.median(snrsOut[index,8:])
    
    try:
        pulserDepth[index] = f(usefulEvent.unixTime)
    except NameError:
        print("Pulser depth calculation bypassed.  Setting run number to pulser depth.")
        pulserDepth[index] = runNumber
    unixtime[index] = usefulEvent.unixTime
    # print(f(usefulEvent.unixTime))
    print(pulserDepth[index])
    
    
    
deltaTPeakOut = peakTimeOut[:,8:] - peakTimeOut[:,:8] 

print("Events processed!")

#Write calculations for condition where we swap HPol RF with HPol soft-trigger.
powerSoftTriggerHpolOut[:,:8] = powerOut[:,:8]
hilbertPeakSoftTriggerHpolOut[:,:8] = hilbertPeakOut[:,:8]
softTriggerPeakTime = np.zeros((numEvents,16))

#Import MCTruth values for AraSim results
if (dataTypeFlag == 1):
    #Initialize arrays that will be saved to .pkl
    thetaTrueOut = np.zeros((numEvents, 16))
    thetaVertexTrueOut = np.zeros(numEvents)
    phiTrueOut = np.zeros((numEvents, 16))
    phiVertexTrueOut = np.zeros(numEvents)
    whichSolTrue = np.zeros(numEvents)
    trueROut = np.zeros((numEvents, 8))
    psiTrueOut = np.zeros((numEvents, 8))
    trueRadiusOut = np.zeros(numEvents)
    
    

#Save data to pandas file
original_df = pd.DataFrame({"runNumber":runNumberOut, "runSubsetNumber":runSubsetNumberOut, "eventNumber":evt_num.tolist(),   "timeStamp":timeStamp.tolist(), "runEventNumber":runEventNumber.tolist(),
                            "thetaReco":thetaRecoOut.tolist(), "phiReco":phiRecoOut.tolist(), "thetaVertex":thetaVertexOut.tolist(), "phiVertex":phiVertexOut.tolist(), 
                            "SNR":snrsOut.tolist(), "vSNR":vSnrOut.tolist(), "hSNR":hSnrOut.tolist(), "whichSol":whichSol.tolist(), "unixtime":unixtime.tolist(), 
                            "power":powerOut.tolist(), "recoR":recoROut.tolist(), "psiReco":psiRecoOut.tolist(),
                            "powerNoiseFromSoftTrigger":powerNoiseFromSoftTriggerOut.tolist(), "pulserDepth":pulserDepth.tolist(),
                           "hilbertPeaks":hilbertPeakOut.tolist(), "peakTimes":peakTimeOut.tolist(), "deltaTPeaks":deltaTPeakOut.tolist(), "hilbertPeaksSoftTriggerHpol":hilbertPeakSoftTriggerHpolOut.tolist(),
                           "timingChannels":timingChannelsOut.tolist(), "timingCutPass":timingCutPass.tolist(), "timingDifferenceMeasured":timingDifferenceMeasured.tolist(), "timingDifferenceExpected":timingDifferenceExpected.tolist()})
 

#Create output directory if it doesn't already exist
# outputFolder += "/run_0"+str(runNumber)+"/

# if (signalTypeFlag == 0):
#     outputFolder += "/rfEvents/"
# elif (signalTypeFlag == 1):
#     outputFolder += "/softTriggerEvents/"
# elif (signalTypeFlag == 2):
#     outputFolder += "/calpulserEvents/"
    
# if (deconvolution):
#     outputFolder += "/deconvolution/"
# else:
#     outputFolder += "/noDeconvolution/"
    
# if gainBalance:
#     outputFolder += "/gainCorrected/"
isExist = os.path.exists(outputFolder)
if not isExist:
  # Create a new directory because it does not exist 
  os.makedirs(outputFolder)

# Commenting this out on 20231011 because I need to save the subset info for spiceCore results.
# if (dataTypeFlag == 2):
#     outfilePath = outputFolder + "/polReco_run%i_%i.pkl"%(runNumber,subsetNumber)
# else:
#     outfilePath = outputFolder + "/polReco_run%i.pkl"%(runNumber)  

# outfilePath = outputFolder + "/polReco_run%i_%i.pkl"%(runNumber,subsetNumber)

outfilePath = outputFolder + "/polReco_run%i.pkl"%(runNumber)  

original_df.to_pickle(outfilePath)
print("Output saved to " + outfilePath)