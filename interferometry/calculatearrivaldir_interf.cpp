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

// AraSim includes
#include "Detector.h"
#include "Event.h"
#include "Position.h"
#include "RaySolver.h"
#include "IceModel.h"

#include "helper.h"
#include "tools.h"

void getCorrMapPeak( TH2D *theCorrMap_input, double &peakTheta, double &peakPhi, double &peakCorr) {

    int _peakZ, _peakTheta, _peakPhi;
    theCorrMap_input->GetMaximumBin(_peakPhi, _peakTheta, _peakZ);

    double maxCorr = theCorrMap_input->GetMaximum();
    double maxPhi = theCorrMap_input->GetXaxis()->GetBinCenter(_peakPhi);
    double maxTheta = theCorrMap_input->GetYaxis()->GetBinCenter(_peakTheta);

    peakCorr = maxCorr;
    peakPhi = maxPhi;
    peakTheta = maxTheta;
} 

//TODO: make radii and numScanned a little more dynamic so that it doesn't have to be hardcoded.  Maybe had it scan the rt_tables folder to have ti compile a list of radii that match the station, then defining the radii array using those values?  numScanned would simply just be the length of radii[]. Should also have it scan AFTER the station number is declared.  4/6/2024
double radii[] = {
      0,  150,  300,  450,  600,  750,  900, 1050, 1200, 1350, 1500, 1650, 1800, 1950,
 2100, 2250, 2400, 2550, 2700, 2850, 3000, 3150, 3300, 3450, 3600, 3750, 3900, 4050,
 4200, 4350, 4500, 4650, 4800, 4950
};
const int numScanned = 34;

int main(int argc, char **argv)
{
    if(argc<6) {
        std::cout << "Usage\n" << argv[0] << " <station> <runnum> <input file> <output_dir> <setup_file> \n";
        std::cout << "e.g.\n" << argv[0] << " 2 http://www.hep.ucl.ac.uk/uhen/ara/monitor/root/run1841/event1841.root setup.txt\n";
        return 0;
    }

    //Todo: check if these values are the same across stations and configs, or if I need to make it dynamic 4/9/2024
    double interpV = 0.4;
    double interpH = 0.625;

    int station = atoi(argv[1]);  //This will be made redundant when we include the setupfile. 4/6/2024
    int runNum = atoi(argv[2]);
    char* inputFile = argv[3];
    char* outputDir = argv[4];
    char* setupfile = argv[5];
    
    

    char outfile_name[400];
    sprintf(outfile_name, "%s/recangle_reco_out_run_%d.root", outputDir, runNum);
    std::cout<<"Output name is "<<outfile_name<<std::endl;
    TFile *fpOut = TFile::Open(outfile_name, "RECREATE");
    if(!fpOut){ std::cerr<<"Cannot open output file "<<fpOut<<std::endl; return -1; }

    TTree *outTree = new TTree("vertexReco", "vertexReco");

    double peakCorrs_out[numScanned];
    double peakThetas_out[numScanned];
    double peakPhis_out[numScanned];
    int peakPol_out[numScanned];
    int peakSol_out[numScanned];
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
    outTree->Branch("peakCorrs", peakCorrs_out, TString::Format("peakCoors_out[%d]/D",numScanned));
    outTree->Branch("peakThetas", peakThetas_out, TString::Format("peakThetas_out[%d]/D",numScanned));
    outTree->Branch("peakPhis", peakPhis_out, TString::Format("peakPhis_out[%d]/D",numScanned));
    outTree->Branch("peakPols", peakPol_out, TString::Format("peakPol_out[%d]/I",numScanned));
    outTree->Branch("peakSols", peakSol_out, TString::Format("peakSol_out[%d]/I",numScanned));
    outTree->Branch("bestTheta", &bestTheta_out, "bestTheta_out/D");
    outTree->Branch("bestPhi", &bestPhi_out, "bestPhi_out/D");
    outTree->Branch("bestCorr", &bestCorr_out, "bestCorr_out/D");
    outTree->Branch("bestR", &bestR_out, "bestR_out/D");
    outTree->Branch("bestSol", &bestSol_out, "bestSol_out/I");
    outTree->Branch("reco_arrivalThetas", reco_arrivalThetas_out, "reco_arrivalThetas_out[16]/D");
    outTree->Branch("reco_arrivalPhis", reco_arrivalPhis_out, "reco_arrivalPhis_out[16]/D");

    //These are simulation specific, and should be set to zero or delete if using a real data reconstruction
    outTree->Branch("trueTheta", &trueTheta_out, "trueTheta_out/D");
    outTree->Branch("truePhi", &truePhi_out, "truePhi_out/D");
    outTree->Branch("trueR", &trueR_out, "trueR_out/D");
    outTree->Branch("trueSol", &likelySol_out, "likelySol_out/D");
    outTree->Branch("true_arrivalThetas", true_arrivalThetas_out, "true_arrivalThetas_out[16]/D");
    outTree->Branch("true_arrivalPhis", true_arrivalPhis_out, "true_arrivalPhis_out[16]/D");
    
    outTree->Branch("weight", &weight_out, "weight_out/D");
    outTree->Branch("snrs", snrs_out, "snrs_out[16]/D");
    outTree->Branch("v_snr", &v_snr_out, "v_snr_out/D");
    outTree->Branch("h_snr", &h_snr_out, "h_snr_out/D");
    
    //These are to accomodate a real data reconstruction
    outTree->Branch("eventNumber", &eventNumber, "eventNumber/I");
    outTree->Branch("unixTime", &unixTime, "unixTime/I");
    outTree->Branch("unixTimeUs", &unixTimeUs, "unixTimeUs/I");
    outTree->Branch("timeStamp", &timeStamp, "timeStamp/I");
    
    //This should work for any double-peak event, as simulated events can show D and R solution in a large enough trigger window.
    outTree->Branch("cutoffTime", &cutoffTime, "cutoffTime[16]/D");    

    printf("------------------\n");
    printf("Output File Setup Complete. Begin correlator setup.\n");
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
    IceModel *icemodel=new IceModel(settings1->ICE_MODEL + settings1->NOFZ*10,settings1->CONSTANTICETHICKNESS * 1000 + settings1->CONSTANTCRUST * 100 + settings1->FIXEDELEVATION * 10 + 0,settings1->MOOREBAY);// creates Antarctica ice model
    Detector *detector = new Detector(settings1, icemodel, setupfile);  
    Report *report = new Report(detector, settings1);
    RaySolver *raySolver = new RaySolver;
    //End attempt at importing setup parameters
    
    // set up correlator and pairs
    RayTraceCorrelator *theCorrelators[numScanned];
    for(int r=0; r<numScanned; r++){
        
        // setup the paths to our ray tracing tables
        double radius = radii[r];
        printf("Loading Radius %.2f\n",radius);
        double angular_size = 1.;
        int icemodel = 0;  //This should be replaced with whatever ice model is used in the setup file. 4/6/2024
        // cout << "Using icemodel = " << icemodel << endl;
        char dirPath[500];
        char refPath[500];
        //std::string topDir = "/mnt/home/baclark/ara/rt_tables";
        std::string topDir = "rt_tables";
        sprintf(dirPath, "%s/arrivaltimes_station_%d_icemodel_%d_radius_%.2f_angle_%.2f_solution_0.root",
            topDir.c_str(), station, icemodel, radius, angular_size
        );
        sprintf(refPath, "%s/arrivaltimes_station_%d_icemodel_%d_radius_%.2f_angle_%.2f_solution_1.root",
            topDir.c_str(), station, icemodel, radius, angular_size
        );

        int numAntennas = 16; //While this should be fine for all stations, it should be dynamic to the settings file as well. 4/6/2024
        theCorrelators[r] = new RayTraceCorrelator(station, numAntennas, radius, angular_size, dirPath, refPath);
        theCorrelators[r]->LoadTables();        
    }
    
    AraGeomTool *geomTool = AraGeomTool::Instance();

    printf("------------------\n");
    printf("Correlator Setup Complete. Begin looping events.\n");
    printf("------------------\n");    
    
    
    //Use the getTrigMasking function to use the same channels that triggering used for the reconstruction
    std::vector<int> excludedChannels;
    for (int i=0; i<16; i++){
        // cout << "detector->GetTrigMasking(i) = " << detector->GetTrigMasking(i) << endl;
        if (not detector->GetTrigMasking(i)){
            cout << "Excluding channel " << i << endl;
            excludedChannels.push_back(i);
        }
    }
    std::map< int, std::vector<int> > pairs_V = theCorrelators[0]->SetupPairs(station, geomTool, AraAntPol::kVertical, excludedChannels);
    std::map< int, std::vector<int> > pairs_H = theCorrelators[0]->SetupPairs(station, geomTool, AraAntPol::kHorizontal, excludedChannels);    
    
    // Check if sim or real data file by checking for existence of AraTree
    TTree *simSettingsTree;
    TTree *simTree;
    simSettingsTree=(TTree*) fp->Get("AraTree");
    bool dataLike = false;
    Event *eventPtr = 0; // it is apparently incredibly important that this be initialized to zero...
    Report *reportPtr = 0; 
    std::vector<double> average_position;
    //data like
    if(!simSettingsTree) { 
        dataLike = true;            
        std::cerr << "Can't find AraTree.  Importing as real data.\n";
        eventTree->SetBranchAddress("event",&rawAtriEvPtr);
        weight = 1;
    }
    // sim like
    else {
        std::cerr << "AraTree exists.  Importing as simulated data.\n";
        eventTree->SetBranchAddress("UsefulAtriStationEvent", &usefulAtriEvPtr);
        weight;
        eventTree->SetBranchAddress("weight", &weight);
        // Detector *detector = 0;
        simSettingsTree=(TTree*) fp->Get("AraTree");
        if(!simSettingsTree) { std::cerr << "Can't find AraTree\n"; return -1; }
        simSettingsTree->SetBranchAddress("detector", &detector);
        simSettingsTree->GetEntry(0);
        average_position = get_detector_cog(detector);
        printf("Detector Center %.2f, %.2f, %.2f \n", average_position[0], average_position[1], average_position[2]);

        simTree=(TTree*) fp->Get("AraTree2");
        if(!simTree) { std::cerr << "Can't find AraTree2\n"; return -1; }
        simTree->SetBranchAddress("event", &eventPtr);
        simTree->SetBranchAddress("report", &reportPtr);
    }

    printf("------------------\n");
    printf("Input files loaded. Set up ray tracing business\n");
    printf("------------------\n");
    
    for(Long64_t event=0;event<numEntries;event++) {
        fp->cd();
        eventTree->GetEntry(event);
        if (dataLike) {
            UsefulAtriStationEvent *usefulAtriEvPtr = new UsefulAtriStationEvent(rawAtriEvPtr, AraCalType::kLatestCalib);
        }
        else {
            simTree->GetEntry(event);
        }
        Position diff_true;
        if (not dataLike){
            std::vector<double> event_position;
            event_position.push_back(eventPtr->Nu_Interaction[0].posnu.GetX());
            event_position.push_back(eventPtr->Nu_Interaction[0].posnu.GetY());
            event_position.push_back(eventPtr->Nu_Interaction[0].posnu.GetZ());
            // printf("Posnu %.2f, %.2f, %.2f \n", event_position[0], event_position[1], event_position[2]);

            std::vector<double> difference;
            difference.push_back(event_position[0] - average_position[0]);
            difference.push_back(event_position[1] - average_position[1]);
            difference.push_back(event_position[2] - average_position[2]);
            // printf("Difference %.2f, %.2f, %.2f \n", difference[0], difference[1], difference[2]);

            
            diff_true.SetXYZ(difference[0], difference[1], difference[2]);
            printf("  Vertex Theta %.2f, Phi %.2f, R %.2f\n", diff_true.Theta()*TMath::RadToDeg(), diff_true.Phi()*TMath::RadToDeg(), diff_true.R());            
        }

        std::cout<<"Looking at event number "<<event<<std::endl;        
        
        //TODO: The double peak finder and cropping isn't playing nice with a simulated event with only one peak.
        std::map<int, TGraph*> interpolatedWaveforms;
        std::vector<double> noiseRms(16);
        
        TCanvas *c = new TCanvas("","", 1600, 1600);
        c->Divide(4,4);
        
        //TODO: Implement independent loop for double peak finder that checks parter VPol and HPol channels to find the cutoff times.
        for (int i=0; i<8; i++){
            TGraph *grV = usefulAtriEvPtr->getGraphFromRFChan(i);
            TGraph *grH = usefulAtriEvPtr->getGraphFromRFChan(i+8);
            TGraph *grIntV = FFTtools::getInterpolatedGraph(grV, 0.5);  //Real data interpolation
            TGraph *grIntH = FFTtools::getInterpolatedGraph(grH, 0.5);  //Real data interpolation
            
            int numNoiseSamples = 100;
            double vPolvoltageSubset[numNoiseSamples];
            double hPolvoltageSubset[numNoiseSamples];
            for(int j=0; j<100; j++){
                vPolvoltageSubset[j] = grIntV->GetY()[j];
                hPolvoltageSubset[j] = grIntH->GetY()[j];
            }
            double noiseRmsChV = TMath::RMS(numNoiseSamples,vPolvoltageSubset);
            double noiseRmsChH = TMath::RMS(numNoiseSamples,hPolvoltageSubset);
            if (noiseRmsChV == 0) {
                noiseRms[i]=1;
            }
            else {
                noiseRms[i]=noiseRmsChV;
            }            
            if (noiseRmsChH == 0) {
                noiseRms[i+8]=1;
            }
            else {
                noiseRms[i+8]=noiseRmsChH;
            }                    
 
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
            
            //Compare Vpol and Hpol peaks, choosing the larger peaked channel to serve as the primary peak.
            double peakThresholdV = noiseRmsChV*4;
            double peakThresholdH = noiseRmsChH*4;
            double peakThreshold;
            // cout << "peakThreshold = " << peakThreshold << endl;
            double firstPeak;
            double secondPeak;
            
            if (vvPeakIntPowers[0]/noiseRmsChV > hhPeakIntPowers[0]/noiseRmsChH) {
                primaryPeakIntPowers = vvPeakIntPowers;
                primaryHitTimes = vvHitTimes;
                peakThreshold = peakThresholdV;
            }
            else {
                primaryPeakIntPowers = hhPeakIntPowers;
                primaryHitTimes = hhHitTimes;
                peakThreshold = peakThresholdH;
            }            
            
            
            if (primaryPeakIntPowers[1] > peakThreshold) {
                cout << "Assuming double-peak signal." << endl;
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
                // cutoffTime[i] = grInt->GetX()[grInt->GetN() - 1];
            }
            else {
                cout << "Assuming single-peak signal." << endl;
                // cout << "grInt->GetN() = " << grInt->GetN() << endl;
                // cout << "grInt->GetX()[grInt->GetN() - 1] = " << grInt->GetX()[grIntV->GetN() - 1] << endl;
                
                cutoffTime[i] = grIntV->GetX()[grIntV->GetN() - 1];
                cutoffTime[i+8] = grIntH->GetX()[grIntH->GetN() - 1];
                // cutoffTime[i] = grInt->GetX()[-1];
            }      
            
            //Draw Vpol
            c->cd(i+1); gPad->SetGrid(1,1);
            grV->Draw();
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
        
        
        for(int i=0; i<16; i++){
            TGraph *gr = usefulAtriEvPtr->getGraphFromRFChan(i);
            // c->cd(i+1); gPad->SetGrid(1,1);
            // gr->Draw();            

            
            cout << "**********************" << endl;
            
            cout << "Ch = " << i << endl;
            
            cout << "initial waveform length = " << gr->GetN() << endl;

            //Put Double peak finder stuff here
            TGraph *grInt = FFTtools::getInterpolatedGraph(gr, 0.5);  //Real data interpolation
            
            cout << "cutoffTime at " << cutoffTime[i] << endl;
            
            //Crop waveform to first peak
            grInt = FFTtools::cropWave(grInt, grInt->GetX()[0], cutoffTime[i]);
            cout << "cropped waveform length = " << grInt->GetN() << endl;
            
            //Either pad or crop waveform to fit NFOUR/2
            if (grInt->GetN() < settings1->NFOUR/2) {
                grInt = FFTtools::padWaveToLength(grInt, settings1->NFOUR/2);
            }
            else if (grInt->GetN() > settings1->NFOUR/2) {
                grInt = FFTtools::cropWave(grInt, grInt->GetX()[0], grInt->GetX()[settings1->NFOUR/2-1]);
            }
            
            cout << "padded waveform length = " << grInt->GetN() << endl;
            interpolatedWaveforms[i] = grInt;                              
            
        }
        char title[500];
        sprintf(title, "waveform.png");   
        c->Print(title);
    
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

        std::map<int, double> weights_V; // V weights
        double tot_weight_v = 0.;
        for(auto iter = pairs_V.begin(); iter!= pairs_V.end(); ++iter){
            int pairNum = iter->first;
            int ant1 = iter->second[0];
            int ant2 = iter->second[1];

            auto gr1_iter = snrs.find(ant1);
            auto gr2_iter = snrs.find(ant2);
            double snr1 = gr1_iter->second;
            double snr2 = gr2_iter->second;
            double snr_prod = snr1 * snr2;
            tot_weight_v+=snr_prod;
            weights_V[pairNum] = snr_prod;
        }
        for(int i=0; i<weights_V.size(); i++){ weights_V[i]/=tot_weight_v; }

        std::map<int, double> weights_H; // H weights
        double tot_weight_h = 0.;
        for(auto iter = pairs_H.begin(); iter!= pairs_H.end(); ++iter){
            int pairNum = iter->first;
            int ant1 = iter->second[0];
            int ant2 = iter->second[1];

            auto gr1_iter = snrs.find(ant1);
            auto gr2_iter = snrs.find(ant2);
            double snr1 = gr1_iter->second;
            double snr2 = gr2_iter->second;
            double snr_prod = snr1 * snr2;
            tot_weight_h+=snr_prod;
            weights_H[pairNum] = snr_prod;
        }
        for(int i=0; i<weights_H.size(); i++){ weights_H[i]/=tot_weight_h; }

        std::vector<TGraph*> corrFunctions_V = theCorrelators[0]->GetCorrFunctions(pairs_V, interpolatedWaveforms, true); 
        std::vector<TGraph*> corrFunctions_H = theCorrelators[0]->GetCorrFunctions(pairs_H, interpolatedWaveforms, true); 

        std::vector<double> peakCorrs;
        std::vector<double> peakThetas;
        std::vector<double> peakPhis;
        std::vector<int> peakPol;
        std::vector<int> peakSol;

        for(int r=0; r<numScanned; r++){
            std::vector<TH2D*> maps;

            maps.push_back(theCorrelators[r]->GetInterferometricMap(pairs_V, corrFunctions_V, 0, weights_V)); // direct solution)
            maps.push_back(theCorrelators[r]->GetInterferometricMap(pairs_V, corrFunctions_V, 1, weights_V)); // reflected solution
            maps.push_back(theCorrelators[r]->GetInterferometricMap(pairs_H, corrFunctions_H, 0, weights_H)); // direct solution
            maps.push_back(theCorrelators[r]->GetInterferometricMap(pairs_H, corrFunctions_H, 1, weights_H)); // reflected solution
            
            // Debugging - JCF 6/7/2023
            TCanvas *c = new TCanvas("","", 1200, 950);            
            maps[0]->Draw("colz");
            maps[0]->GetXaxis()->SetTitle("Phi [deg]");
            maps[0]->GetYaxis()->SetTitle("Theta [deg]");
            maps[0]->GetYaxis()->SetTitleSize(0.05);
            maps[0]->GetYaxis()->SetLabelSize(0.03);
            maps[0]->GetYaxis()->SetTitleOffset(0.6);
            maps[0]->GetXaxis()->SetTitleSize(0.05);
            maps[0]->GetXaxis()->SetLabelSize(0.03);
            maps[0]->GetXaxis()->SetTitleOffset(0.6);
            gStyle->SetOptStat(0);
            maps[0]->GetXaxis()->CenterTitle();
            maps[0]->GetYaxis()->CenterTitle();
            gPad->SetRightMargin(0.15);
            char title[500];
            sprintf(title,"%.2f.png", radii[r]);
            c->SaveAs(title);
            delete c;
            // End debugging

            std::vector<double> bestOne;
            for(int i=0; i<4; i++){
                double peakCorr, peakTheta, peakPhi;
                getCorrMapPeak(maps[i], peakTheta, peakPhi, peakCorr);
                // printf("      i %d, peakCorr %e\n", i, peakCorr);
                bestOne.push_back(peakCorr);
            }
            auto it = max_element(std::begin(bestOne), std::end(bestOne));
            int element = distance(bestOne.begin(), it);
            // printf("    Best Option is %d with corr %e \n",element, *it);
            double peakCorr, peakTheta, peakPhi;
            getCorrMapPeak(maps[element], peakTheta, peakPhi, peakCorr);
            peakCorrs.push_back(peakCorr);

            peakThetas.push_back(peakTheta); 
            peakPhis.push_back(peakPhi);
            
            if(element < 2){ peakPol.push_back(0); }
            else{peakPol.push_back(1);}
            
            if(element==0 || element ==2){ peakSol.push_back(0);}
            else if(element==1 || element ==3 ){ peakSol.push_back(1);}

            printf("    Correlated radius %.2f, Corr %.4f \n",radii[r], peakCorr);

            for(int i=0; i<4; i++){
                delete maps[i];
            }
        }
        auto it = max_element(std::begin(peakCorrs), std::end(peakCorrs));
        int element = distance(peakCorrs.begin(), it);

        // stash the output
        fpOut->cd();
        // reconstructed quantities first
        for(int r=0 ;r<numScanned; r++){
            peakCorrs_out[r] = peakCorrs[r];
            peakThetas_out[r] = peakThetas[r];
            peakPhis_out[r] = peakPhis[r];
            peakPol_out[r] = peakPol[r];
            peakSol_out[r] = peakSol_out[r];
        }
        bestTheta_out = peakThetas[element];
        bestPhi_out = peakPhis[element];
        bestCorr_out = peakCorrs[element];
        bestR_out = radii[element];
        bestSol_out = peakSol[element];
        for(int i=0; i<16; i++){

            int this_binTheta, this_binPhi;
            theCorrelators[element]->ConvertAngleToBins(peakThetas[element], peakPhis[element], 
                this_binTheta, this_binPhi);

            double this_arrivalTheta, this_arrivalPhi;
            theCorrelators[element]->LookupArrivalAngles(i, peakSol[element], 
                this_binTheta, this_binPhi,
                this_arrivalTheta, this_arrivalPhi
            );
            reco_arrivalThetas_out[i] = this_arrivalTheta*180/PI; //Previous saved in radians, Converting to degrees - JCF 4/11/2024
            reco_arrivalPhis_out[i] = this_arrivalPhi*180/PI; 
        }
        
        // then the true quantities (if applicable)
        if (not dataLike) {
            // int likely_sol = guess_triggering_solution(eventPtr, reportPtr);
            int likely_sol = 0;  //Forcing solution to zero since we set up the double peak finder to look for the D pulse.
            std::map<int, double> thetas_truth = get_value_from_mc_truth("theta", likely_sol, reportPtr);
            std::map<int, double> phis_truth = get_value_from_mc_truth("phi", likely_sol, reportPtr);
        
            for(int i=0; i<16; i++){
                double this_true_theta = thetas_truth.find(i)->second;
                double this_true_phi = phis_truth.find(i)->second;
                // printf("  Ant %d, True Arrival Theta %.2f, Reco Arrival Phi %.2f \n",
                //     i, this_true_theta, this_true_phi
                // );
                true_arrivalThetas_out[i] = this_true_theta*180/PI;
                true_arrivalPhis_out[i] = this_true_phi*180/PI; 
            }
            trueTheta_out = diff_true.Theta() * TMath::RadToDeg();  //Converted to match the zenith in the reconstruction calculation.
            truePhi_out = diff_true.Phi() * TMath::RadToDeg();
            trueR_out = diff_true.R();
            likelySol_out = likely_sol; //The output for this is a large negative number, even when likely_sol is hardcoded to zero.
            
        }
        
        weight_out = weight;
        for(int i=0; i<16; i++){
            snrs_out[i] = snrs[i];
        }
        v_snr_out = the_snr_v;
        h_snr_out = the_snr_h;

        outTree->Fill();

    
    }

    // write output
    fpOut->Write();
    fpOut->Close();
    delete fpOut;


    fp->Close();
    delete fp;
}