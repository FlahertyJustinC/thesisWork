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
    
    double bestTheta_out;
    double bestPhi_out;
    double bestCorr_out;
    double bestR_out;
    int bestSol_out;
    double reco_arrivalThetas_out[16];
    double reco_arrivalPhis_out[16];

    double trueTheta_out;
    double truePhi_out;
    double trueR_out;
    int likelySol_out;
    double true_arrivalThetas_out[16];
    double true_arrivalPhis_out[16];
    double cutoffTime[16];    
    
    int eventNumber;
    int unixTime;
    int unixTimeUs;
    int timeStamp;    

    double weight_out;
    double weight;
    double snrs_out[16];
    double v_snr_out;
    double h_snr_out;
    
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
    if(!eventTree) { std::cerr << "Can't find eventTree\n"; return -1; }
    Long64_t numEntries=eventTree->GetEntries();
    cout << "eventTree has " << numEntries << " entries." << endl;
    RawAtriStationEvent *rawAtriEvPtr=0;
    
    //Try importing paramters using the setup file
    Settings *settings1 = new Settings();
    settings1->ReadFile(setupfile);
    IceModel *icemodel = new IceModel(settings1->ICE_MODEL + settings1->NOFZ*10,settings1->CONSTANTICETHICKNESS * 1000 + settings1->CONSTANTCRUST * 100 + settings1->FIXEDELEVATION * 10 + 0,settings1->MOOREBAY);// creates Antarctica ice model
    Detector *detector = new Detector(settings1, icemodel, setupfile);
    // Report *report = new Report(detector, settings1);
    // RaySolver *raySolver = new RaySolver;
    //End attempt at importing setup parameters    
    

    // printf("------------------\n");
    // printf("Begin looping events.\n");
    // printf("------------------\n");    
    
    
    //Use the getTrigMasking function to use the same channels that triggering used for the reconstruction
    std::vector<int> excludedChannels;
    // if (settings1->DETECTOR_STATION==4){
    //     cout << "Using station 4.  Excluding channels 0,4,8 in addition to the masking." << endl;
    //     excludedChannels.push_back(0);  //Testing removing a channel for A4 as only the R pulse makes it into the waveform.
    //     excludedChannels.push_back(4); 
    //     excludedChannels.push_back(8);
    // }
    for (int i=0; i<16; i++){
        // cout << "detector->GetTrigMasking(i) = " << detector->GetTrigMasking(i) << endl;
        if (not detector->GetTrigMasking(i)){
            cout << "Excluding channel " << i << endl;
            excludedChannels.push_back(i);
        }
    }  
    
    // Check if sim or real data file by checking for existence of AraTree
    TTree *simSettingsTree;
    TTree *simTree;
    simSettingsTree=(TTree*) fp->Get("AraTree");
    bool dataLike;
    Event *eventPtr = 0; // it is apparently incredibly important that this be initialized to zero...
    // Report *reportPtr = 0; 
    //data like
    if(!simSettingsTree) { 
        dataLike = true;            
        std::cerr << "Can't find AraTree.  Importing as real data.\n";
        eventTree->SetBranchAddress("event",&rawAtriEvPtr);
        weight = 1;
    }
    // sim like
    else {
        dataLike = false;
        std::cerr << "AraTree exists.  Importing as simulated data.\n";
        eventTree->SetBranchAddress("UsefulAtriStationEvent", &usefulAtriEvPtr);
        weight;
        eventTree->SetBranchAddress("weight", &weight);
        // Detector *detector = 0;
        simSettingsTree=(TTree*) fp->Get("AraTree");
        if(!simSettingsTree) { std::cerr << "Can't find AraTree\n"; return -1; }
        simSettingsTree->SetBranchAddress("detector", &detector);
        simSettingsTree->GetEntry(0);

        simTree=(TTree*) fp->Get("AraTree2");
        if(!simTree) { std::cerr << "Can't find AraTree2\n"; return -1; }
        simTree->SetBranchAddress("event", &eventPtr);
        // simTree->SetBranchAddress("report", &reportPtr);
    }

    printf("------------------\n");
    printf("Input files loaded. Begin looping events.\n");
    printf("------------------\n");
    
    for(Long64_t event=0;event<numEntries;event++) {
    // for(Long64_t event=0;event<100;event++) {        
    // for(Long64_t event=53534;event<53634;event++) {        
    // for(Long64_t event=0;event<10;event++) {  //Debugging and running on first few events.
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
        else {
            // cout << "Triggering sim-like condition." << endl;
            simTree->GetEntry(event);
        }
        // cout << "*";
        std::cout<<"Event number: \t"<<event<<std::endl;        
        
        //TODO: The double peak finder and cropping isn't playing nice with a simulated event with only one peak.
        std::map<int, TGraph*> interpolatedWaveforms;
        std::vector<double> noiseRms(16);
        
        TCanvas *c = new TCanvas("","", 1600, 1600);
        c->Divide(4,4);
        
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
                cutoffTime[i] = firstPeak + (secondPeak-firstPeak)/2;
                cutoffTime[i+8] = firstPeak + (secondPeak-firstPeak)/2;
                
                //Add condition where if the primary channel isn't in the excluded channel list, it adds to the double-peak channel counter
                bool channelNotExcluded = std::find(excludedChannels.begin(), excludedChannels.end(), primaryChannel) == excludedChannels.end();
                // cout << "PrimaryChannel not excluded = " << channelNotExcluded << endl;
                if (channelNotExcluded) {
                    doublePeakCounter++;
                    // cout << "doublePeakCounter = " << doublePeakCounter << endl;
                }
            }
            else {
                // cout << "Assuming single-peak signal." << endl;                
                cutoffTime[i] = grIntV->GetX()[grIntV->GetN() - 1];
                cutoffTime[i+8] = grIntH->GetX()[grIntH->GetN() - 1];
            }
            if (debugMode){
                //Draw Vpol
                c->cd(i+1); gPad->SetGrid(1,1);
                grV->Draw();
                char vTitle[500];
                sprintf(vTitle,"Ch. %.2d", i);
                grV->SetTitle(vTitle);
                TLine *l1v = new TLine(vvHitTimes[0], -1500, vvHitTimes[0], 1500);
                l1v->SetLineColorAlpha(kBlue, 1);
                l1v->Draw();    
                TLine *l2v = new TLine(vvHitTimes[1], -1500, vvHitTimes[1], 1500);
                l2v->SetLineColorAlpha(kRed, 1);
                l2v->Draw();    
                TLine *l3v = new TLine(-5000, vvPeakIntPowers[0], 5000, vvPeakIntPowers[0]);
                l3v->SetLineColorAlpha(kBlue, 1);
                l3v->Draw();    
                TLine *l4v = new TLine(-5000, vvPeakIntPowers[1], 5000, vvPeakIntPowers[1]);
                l4v->SetLineColorAlpha(kRed, 1);
                l4v->Draw();  
                TLine *l5v = new TLine(-5000, peakThresholdV, 5000, peakThresholdV);
                l5v->SetLineColorAlpha(kGreen, 1);
                l5v->Draw(); 
                TLine *l6v = new TLine(cutoffTime[i], -1500, cutoffTime[i], 1500);
                // l6->SetLineColorAlpha(kPink, 1);
                l6v->SetLineStyle(kDashed);
                l6v->Draw();

                //Draw Hpol
                c->cd(i+1+8); gPad->SetGrid(1,1);
                grH->Draw();
                char hTitle[500];
                sprintf(hTitle,"Ch. %.2d", i+8);
                grH->SetTitle(hTitle);
                TLine *l1h = new TLine(hhHitTimes[0], -1500, hhHitTimes[0], 1500);
                l1h->SetLineColorAlpha(kBlue, 1);
                l1h->Draw();    
                TLine *l2h = new TLine(hhHitTimes[1], -1500, hhHitTimes[1], 1500);
                l2h->SetLineColorAlpha(kRed, 1);
                l2h->Draw();    
                TLine *l3h = new TLine(-5000, hhPeakIntPowers[0], 5000, hhPeakIntPowers[0]);
                l3h->SetLineColorAlpha(kBlue, 1);
                l3h->Draw();    
                TLine *l4h = new TLine(-5000, hhPeakIntPowers[1], 5000, hhPeakIntPowers[1]);
                l4h->SetLineColorAlpha(kRed, 1);
                l4h->Draw();  
                TLine *l5h = new TLine(-5000, peakThresholdH, 5000, peakThresholdH);
                l5h->SetLineColorAlpha(kGreen, 1);
                l5h->Draw(); 
                TLine *l6h = new TLine(cutoffTime[i+8], -1500, cutoffTime[i+8], 1500);
                // l6->SetLineColorAlpha(kPink, 1);
                l6h->SetLineStyle(kDashed);
                l6h->Draw();
            }
            
            if (debugMode){
                cout << "**********************" << endl;           
                cout << "Ch = " << i << endl;            
                cout << "initial waveform length = " << grV->GetN() << endl;            
                cout << "interpolated waveform length = " << grIntV->GetN() << endl;            
                cout << "cutoffTime at " << cutoffTime[i] << endl; 
            }
            //Crop waveform to first peak
            grIntV = FFTtools::cropWave(grIntV, grIntV->GetX()[0], cutoffTime[i]);
            if (debugMode){
                cout << "cropped waveform length = " << grIntV->GetN() << endl;            
            }
            //Either pad or crop waveform to fit NFOUR/2
            if (grIntV->GetN() < settings1->NFOUR/2) {
                grIntV = FFTtools::padWaveToLength(grIntV, settings1->NFOUR/2);
            }
            else if (grIntV->GetN() > settings1->NFOUR/2) {
                grIntV = FFTtools::cropWave(grIntV, grIntV->GetX()[0], grIntV->GetX()[settings1->NFOUR/2-1]);
            }
            if (debugMode){
                cout << "padded waveform length = " << grIntV->GetN() << endl;
            }
            interpolatedWaveforms[i] = grIntV;          
            
            if (debugMode){
                cout << "**********************" << endl;           
                cout << "Ch = " << i+8 << endl;            
                cout << "initial waveform length = " << grH->GetN() << endl;            
                cout << "interpolated waveform length = " << grIntH->GetN() << endl;            
                cout << "cutoffTime at " << cutoffTime[i+8] << endl;  
            }
            //Crop waveform to first peak
            grIntH = FFTtools::cropWave(grIntH, grIntH->GetX()[0], cutoffTime[i+8]);
            if (debugMode){
                cout << "cropped waveform length = " << grIntH->GetN() << endl;    
            }
            //Either pad or crop waveform to fit NFOUR/2
            if (grIntH->GetN() < settings1->NFOUR/2) {
                grIntH = FFTtools::padWaveToLength(grIntH, settings1->NFOUR/2);
            }
            else if (grIntH->GetN() > settings1->NFOUR/2) {
                grIntH = FFTtools::cropWave(grIntH, grIntH->GetX()[0], grIntH->GetX()[settings1->NFOUR/2-1]);
            }
            if (debugMode){
                cout << "padded waveform length = " << grIntV->GetN() << endl;
            }
            interpolatedWaveforms[i+8] = grIntH;              
            
            
        }

        if (debugMode){
            char title[500];
            sprintf(title, "waveform.png");
            c->Print(title);
        }
        
        delete c;
    
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
        // cout << "the_snr_v = " << the_snr_v << endl;
        // cout << "the_snr_h = " << the_snr_h << endl;
        v_snr_out = the_snr_v;
        h_snr_out = the_snr_h;           
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
        
        delete usefulAtriEvPtrOut;
        delete usefulAtriEvPtr;

    
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