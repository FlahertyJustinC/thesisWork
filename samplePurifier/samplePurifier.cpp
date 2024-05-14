#include <iostream>

// ROOT Includes
#include "TFile.h"
#include "TTree.h"
#include "TGraph.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TLine.h"
#include "TMultiGraph.h"

// ARA Includes
#include "AraGeomTool.h"
#include "RayTraceCorrelator.h"
#include "UsefulAtriStationEvent.h"
#include "/cvmfs/ara.opensciencegrid.org/trunk/centos7/source/libRootFftwWrapper/include/FFTtools.h"

UsefulAtriStationEvent *usefulAtriEvPtr;
UsefulAtriStationEvent *usefulAtriEvPtrOut;

// AraSim includes
#include "Detector.h"
#include "Event.h"
#include "Position.h"
#include "RaySolver.h"
#include "IceModel.h"

#include "tools.h"

bool debugMode = false;
// int passedChannelThreshold = 3;

double calculateNoiseRMS(TGraph *gr, int sampleNumber=100) {
    int totalSampleNumber = abs(sampleNumber);
    int waveformLength = gr->GetN();
    double voltageSubset[totalSampleNumber];
    if (sampleNumber<0){
        //Loop over beginning of waveform for noise
        for(int j=0; j<totalSampleNumber; j++){
            voltageSubset[j] = gr->GetY()[j];
        }
    }
    else {
        //Loop over end of waveform for noise
        for(int j=0; j<totalSampleNumber; j++){
            voltageSubset[j] = gr->GetY()[waveformLength-1-j];
        }        
    }
    double noiseRms = TMath::RMS(totalSampleNumber,voltageSubset);
    
    if (noiseRms == 0) {
        noiseRms=1;
    }
    
    return noiseRms;
}

int main(int argc, char **argv)
{
    if(argc<5) {
        std::cout << "Usage\n" << argv[0] << " <runnum> <number of channel pairs with double peaks> <setup_file> <input file> <output_dir> \n";
        std::cout << "e.g.\n" << argv[0] << " 2 4 setup.txt event012559.root output/\n";
        return 0;
    }

    //Todo: check if these values are the same across stations and configs, or if I need to make it dynamic 4/9/2024
    double interpV = 0.4;
    double interpH = 0.625;

    // int station = atoi(argv[1]);  //This will be made redundant when we include the setupfile. 4/6/2024
    char* runNum = argv[1];
    int passedChannelThreshold = atoi(argv[2]);
    char* setupfile = argv[3];
    char* inputFile = argv[4];
    char* outputDir = argv[5];

    char outfile_name[400];
    sprintf(outfile_name, "%s/purifiedSample_run_%s.root", outputDir, runNum);
    std::cout<<"Output name is "<<outfile_name<<std::endl;
    TFile *fpOut = TFile::Open(outfile_name, "RECREATE");
    if(!fpOut){ std::cerr<<"Cannot open output file "<<fpOut<<std::endl; return -1; }

    TTree *outTree = new TTree("eventTree", "eventTree");
    outTree->Branch("UsefulAtriStationEvent", &usefulAtriEvPtrOut);  
    
    int passedEventCounter = 0;

    printf("------------------\n");
    printf("Output File Setup Complete.\n");
    printf("------------------\n");

    printf("Opening file...\n");
    TFile *fp = TFile::Open(inputFile);
    if(!fp) { std::cerr << "Can't open file\n"; return -1; }
    printf("File opened!\n");
    
    //TODO: Implement the code below that differentiates between real and simulated data. 4/6/2024
    
    //Import eventTree
    TTree *eventTree = (TTree*) fp->Get("eventTree");
    printf("Event tree opened!\n");
    if(!eventTree) { std::cerr << "Can't find eventTree\n"; exit(0); return -1; }
    Long64_t numEntries=eventTree->GetEntries();
    cout << "eventTree has " << numEntries << " entries." << endl;
    RawAtriStationEvent *rawAtriEvPtr=0;
    
    //Try importing paramters using the setup file
    Settings *settings1 = new Settings();
    settings1->ReadFile(setupfile);
    IceModel *icemodel = new IceModel(settings1->ICE_MODEL + settings1->NOFZ*10,settings1->CONSTANTICETHICKNESS * 1000 + settings1->CONSTANTCRUST * 100 + settings1->FIXEDELEVATION * 10 + 0,settings1->MOOREBAY);// creates Antarctica ice model
    Detector *detector = new Detector(settings1, icemodel, setupfile);
    
    //Use the getTrigMasking function to use the same channels that triggering used for the reconstruction
    std::vector<int> excludedChannels;
    for (int i=0; i<16; i++){
        // cout << "detector->GetTrigMasking(i) = " << detector->GetTrigMasking(i) << endl;
        if (not detector->GetTrigMasking(i)){
            cout << "Excluding channel " << i << endl;
            excludedChannels.push_back(i);
        }
    }  
    
    // Check if sim or real data file by checking for existence of AraTree
    TTree *simSettingsTree;
    simSettingsTree=(TTree*) fp->Get("AraTree");
    bool dataLike;
    Event *eventPtr = 0; // it is apparently incredibly important that this be initialized to zero...
    //data like
    if(!simSettingsTree) { 
        dataLike = true;            
        std::cerr << "Can't find AraTree.  Importing as real data.\n";
        eventTree->SetBranchAddress("event",&rawAtriEvPtr);
    }
    // sim like
    else {
        dataLike = false;
        std::cerr << "AraTree exists.  Importing as simulated data.\n";
        eventTree->SetBranchAddress("UsefulAtriStationEvent", &usefulAtriEvPtr);
    }

    printf("------------------\n");
    printf("Input files loaded. Begin looping events.\n");
    printf("------------------\n");
    
    for(Long64_t event=0;event<numEntries;event++) {
    // for(Long64_t event=0;event<500;event++) {    
        //Initialize counter for channels with double-peak
        int doublePeakCounter=0;
        
        fp->cd();
        eventTree->GetEntry(event);
        
        // if (rawAtriEvPtr->eventNumber != 53535){  //Test for A4 Run 6128
        //     continue;
        // }
        
        // if (rawAtriEvPtr->eventNumber != 53625){  //Test for A4 Run 6128
        //     continue;
        // }        

        // if (rawAtriEvPtr->eventNumber != 387){  //Test for A4 Run 6119
        //     continue;
        // }
        
        // if (rawAtriEvPtr->eventNumber != 679){  //Test for A2 Run 12559
        //     continue;
        // }        
        
        
        if (dataLike) {
            // cout << "Triggering data-like condition." << endl;
            delete usefulAtriEvPtr;  //Need to delete the initialized pointer in order to create a new one.
            usefulAtriEvPtr = new UsefulAtriStationEvent(rawAtriEvPtr, AraCalType::kLatestCalib);
        }
        // cout << "*";
        if (event%10 == 0) {
            std::cout<<"Event number: \t"<<event<<std::endl;        
        }
        
        std::map<int, TGraph*> interpolatedWaveforms;
        std::vector<double> noiseRms(16);
        
        //TODO: Implement independent loop for double peak finder that checks parter VPol and HPol channels to find the cutoff times.
        for (int i=0; i<8; i++){
            
            int primaryChannel;  //Initializing variable to count primary channel used in the peak finding.
            TGraph *grV = usefulAtriEvPtr->getGraphFromRFChan(i);
            TGraph *grH = usefulAtriEvPtr->getGraphFromRFChan(i+8);
            TGraph *grIntV = FFTtools::getInterpolatedGraph(grV, 0.5);  //Real data interpolation
            TGraph *grIntH = FFTtools::getInterpolatedGraph(grH, 0.5);  //Real data interpolation
            

 
            vector<double> vvHitTimes; // vector of hit times
            vector<double> vvPeakIntPowers; // vector of peak values   
            vector<double> hhHitTimes; // vector of hit times
            vector<double> hhPeakIntPowers; // vector of peak values
            vector<double> primaryHitTimes; // vector of hit times
            vector<double> primaryPeakIntPowers; // vector of peak values            
            
            int numSearchPeaks=2; //only look for two peaks  //TODO:  Maybe make this a console argument?  Though the cropping is based on having two peaks 4/10/2024
            double peakSeparation=50.0;  //Minimum separation between peaks.  Closest seperation is expected to be ~100 ns.
            getAbsMaximum_N(grIntV, numSearchPeaks, peakSeparation, vvHitTimes, vvPeakIntPowers);
            getAbsMaximum_N(grIntH, numSearchPeaks, peakSeparation, hhHitTimes, hhPeakIntPowers);
            
            double noiseRmsChV = calculateNoiseRMS(grIntV);
            double noiseRmsChH = calculateNoiseRMS(grIntH);          
            
            //Compare Vpol and Hpol peaks, choosing the larger peaked channel to serve as the primary peak.
            double peakThresholdV = noiseRmsChV*4;
            double peakThresholdH = noiseRmsChH*4;
            double peakThreshold;
            // cout << "peakThreshold = " << peakThreshold << endl;
            double firstPeak;
            double secondPeak;
            
            //Adding condition where if the primary peaks are less than their peak thresholds, we grab noise from the end of the waveform and calculate a new threshold.
            if (vvPeakIntPowers[0] < peakThresholdV) {
                noiseRmsChV = calculateNoiseRMS(grIntV, -100);
                peakThresholdV = noiseRmsChV*4;
            }
            if (hhPeakIntPowers[0] < peakThresholdH) {
                noiseRmsChH = calculateNoiseRMS(grIntH, -100);
                peakThresholdH = noiseRmsChH*4;
            }            
            
            
            if (vvPeakIntPowers[0]/noiseRmsChV > hhPeakIntPowers[0]/noiseRmsChH) {
                primaryChannel = i;
                primaryPeakIntPowers = vvPeakIntPowers;
                primaryHitTimes = vvHitTimes;
                peakThreshold = peakThresholdV;
            }
            else {
                primaryChannel = i+8;
                primaryPeakIntPowers = hhPeakIntPowers;
                primaryHitTimes = hhHitTimes;
                peakThreshold = peakThresholdH;
            }
                
            noiseRms[i] = noiseRmsChV;
            noiseRms[i+8]=noiseRmsChH;
            
            
            if (primaryPeakIntPowers[1] > peakThreshold) {
                // cout << "Assuming double-peak signal." << endl;
                if(primaryHitTimes[1]>primaryHitTimes[0]) {
                    firstPeak = primaryHitTimes[0];
                    secondPeak = primaryHitTimes[1];
                } 
                else {
                    firstPeak = primaryHitTimes[1];
                    secondPeak = primaryHitTimes[0];
                }
                
                //Add condition where if the primary channel isn't in the excluded channel list, it adds to the double-peak channel counter
                bool channelNotExcluded = std::find(excludedChannels.begin(), excludedChannels.end(), primaryChannel) == excludedChannels.end();
                if (channelNotExcluded) {
                    doublePeakCounter++;
                }
            }
            interpolatedWaveforms[i] = grIntV;          
            interpolatedWaveforms[i+8] = grIntH;       
            delete grV, grH, grIntV, grIntH;
            
        }
    
        std::map<int, double> snrs; // map of channels to SNRs
        for(int i=0; i<16; i++){
            
            double peak_max = TMath::MaxElement(interpolatedWaveforms[i]->GetN(), interpolatedWaveforms[i]->GetY());
            double peak_min = abs(TMath::MinElement(interpolatedWaveforms[i]->GetN(), interpolatedWaveforms[i]->GetY()));
            if(peak_min > peak_max){
                snrs[i] = peak_min/noiseRms[i];
            }
            else{
                snrs[i] = peak_max/noiseRms[i];
            }
        }
        std::vector<double> snrs_v;
        std::vector<double> snrs_h;
        for(int i=0; i<8; i++) snrs_v.push_back(snrs[i]);
        for(int i=8; i<16; i++) snrs_h.push_back(snrs[i]);
        //Grabs median-ish snr (third largest of eight) and stores as the VSNR and HSNR. 
        sort(snrs_v.begin(), snrs_v.end(), greater<double>()); // sort largest to smallest
        sort(snrs_h.begin(), snrs_h.end(), greater<double>()); // sort largest to smallest
        double the_snr_v = snrs_v[2];
        double the_snr_h = snrs_h[2];
        
        //Adding condition where if the SNR is less than 8, we bypass the event.          
        if (usefulAtriEvPtr->isCalpulserEvent()) {
            // cout << "Event is calpulser.  Bypassing event." << endl;
            continue;
        }
        if (usefulAtriEvPtr->isSoftwareTrigger()) {
            // cout << "Event is Software Trigger.  Bypassing event." << endl;
            continue;
        }   
        if (the_snr_v < 8 and the_snr_h < 8) {
            // cout << "Event below SNR threshold.  Bypassing event." << endl;
            // outTree->Fill();
            continue;
        }
        //End SNR cut
        if (doublePeakCounter < passedChannelThreshold) {
            // cout << "Event has too few channels with double peaks.  Bypassing event." << endl;
            // outTree->Fill();
            continue;            
        }  
        
        usefulAtriEvPtrOut=usefulAtriEvPtr;
        
        passedEventCounter++;

        outTree->Fill();
        
        // delete usefulAtriEvPtrOut;
        // delete usefulAtriEvPtr;
        // delete interpolatedWaveforms;

    
    }

    // write output
    fpOut->Write();
    fpOut->Close();
    delete fpOut;

    cout << "*************************" << endl;
    cout << passedEventCounter << " of " << numEntries << " passed." << endl;
    cout<<"Output written to:\t " <<outfile_name<<endl;
    
    fp->Close();
    delete fp;
}