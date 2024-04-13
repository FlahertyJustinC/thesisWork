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
    if(argc<5) {
        std::cout << "Usage\n" << argv[0] << " <station> <runnum> <input file> <output_dir> \n";
        std::cout << "e.g.\n" << argv[0] << " 2 http://www.hep.ucl.ac.uk/uhen/ara/monitor/root/run1841/event1841.root\n";
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

    // set up correlator and pairs
    RayTraceCorrelator *theCorrelators[numScanned];
    for(int r=0; r<numScanned; r++){
        
        // setup the paths to our ray tracing tables
        double radius = radii[r];
        printf("Loading Radius %.2f\n",radius);
        double angular_size = 1.;
        int iceModel = 0;  //This should be replaced with whatever ice model is used in the setup file. 4/6/2024
        char dirPath[500];
        char refPath[500];
        //std::string topDir = "/mnt/home/baclark/ara/rt_tables";
        std::string topDir = "rt_tables";
        sprintf(dirPath, "%s/arrivaltimes_station_%d_icemodel_%d_radius_%.2f_angle_%.2f_solution_0.root",
            topDir.c_str(), station, iceModel, radius, angular_size
        );
        sprintf(refPath, "%s/arrivaltimes_station_%d_icemodel_%d_radius_%.2f_angle_%.2f_solution_1.root",
            topDir.c_str(), station, iceModel, radius, angular_size
        );

        int numAntennas = 16; //While this should be fine for all stations, it should be dynamic to the settings file as well. 4/6/2024
        theCorrelators[r] = new RayTraceCorrelator(station, numAntennas, radius, angular_size, dirPath, refPath);
        theCorrelators[r]->LoadTables();        
    }
    
    AraGeomTool *geomTool = AraGeomTool::Instance();

    printf("------------------\n");
    printf("Correlator Setup Complete. Make Output Files\n");
    printf("------------------\n");

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
    printf("Output File Setup Complete. Begin looping events\n");
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
        // simSettingsTree->SetBranchAddress("detector", &detector);
        // simSettingsTree->GetEntry(0);
        // std::vector<double> average_position = get_detector_cog(detector);
        // printf("Detector Center %.2f, %.2f, %.2f \n", average_position[0], average_position[1], average_position[2]);


        // // TTree *simTree;
        // Event *eventPtr = 0; // it is apparently incredibly important that this be initialized to zero...
        // Report *reportPtr = 0;
        // simTree=(TTree*) fp->Get("AraTree2");
        // if(!simTree) { std::cerr << "Can't find AraTree2\n"; return -1; }
        // simTree->SetBranchAddress("event", &eventPtr);
        // simTree->SetBranchAddress("report", &reportPtr);
    }

    printf("------------------\n");
    printf("Input files loaded. Set up ray tracing business\n");
    printf("------------------\n");
    
    //Old setup method
    // RaySolver *raySolver = new RaySolver;
    // IceModel *iceModel = new IceModel(0 + 1*10, 0, 0); //TODO: Should be dictated by setup file. 4/6/2024
    // Settings *settings = new Settings();
    // // configure settings    TODO: All of these should be dictated by the setup file. 4/6/2024
    // settings->Z_THIS_TOLERANCE = 1; // higher tolerance
    // settings->Z_TOLERANCE = 0.05;
    // settings->NOFZ=1; // make sure n(z) is turned on
    // settings->RAY_TRACE_ICE_MODEL_PARAMS = 0; // set the ice model as user requested
    //End old setup method.
    
    // //Try importing paramters using the setup file
    // Settings *settings1 = new Settings();
    // settings1->ReadFile(setupfile); 
    // IceModel *icemodel=new IceModel(settings1->ICE_MODEL + settings1->NOFZ*10,settings1->CONSTANTICETHICKNESS * 1000 + settings1->CONSTANTCRUST * 100 + settings1->FIXEDELEVATION * 10 + 0,settings1->MOOREBAY);// creates Antarctica ice model
    // Detector *detector = new Detector(settings1, icemodel, setupfile);  
    // Report *report = new Report(detector, settings1);
    // RaySolver *raySolver = new RaySolver;
    // //End attempt at importing setup parameters
    
    for(Long64_t event=0;event<numEntries;event++) {
        fp->cd();
        eventTree->GetEntry(event);
        if (dataLike) {
            UsefulAtriStationEvent *usefulAtriEvPtr = new UsefulAtriStationEvent(rawAtriEvPtr, AraCalType::kLatestCalib);
        }
        // else {
        //     simTree->GetEntry(event);
        // }

        std::cout<<"Looking at event number "<<event<<std::endl;
        
        // std::map<int, TGraph*> interpolatedWaveforms;
        // for(int i=0; i<16; i++){
        //     TGraph *gr = usefulAtriEvPtr->getGraphFromRFChan(i);
        //     // TGraph *grInt = FFTtools::getInterpolatedGraph(gr,i<8?interpV:interpH);
        //     TGraph *grInt = FFTtools::getInterpolatedGraph(gr,0.5);
        //     cout << "waveform length = " << grInt->GetN() << endl;
        //     interpolatedWaveforms[i] = grInt;
        //     delete gr;
        // }        
        
        
        //TODO: The double peak finder and cropping isn't playing nice with a simulated event with only one peak.
        std::map<int, TGraph*> interpolatedWaveforms;
        
        //Debugging, making plots of waveforms
        TGraph *g[16];
        TCanvas *c = new TCanvas("","", 1200, 1200);
        TMultiGraph *mg = new TMultiGraph();        
        for(int i=0; i<16; i++){
            TGraph *gr = usefulAtriEvPtr->getGraphFromRFChan(i);
            
            cout << "initial waveform length = " << gr->GetN() << endl;
            // cout << "time[0] = " << gr->GetX()[0] << endl;
            // cout << "time[1] = " << gr->GetX()[1] << endl;
            // cout << "time[2] = " << gr->GetX()[2] << endl;

            //Put Double peak finder stuff here
            TGraph *grInt = FFTtools::getInterpolatedGraph(gr, 0.5);  //Real data interpolation
            
            vector<double> vvHitTimes; // vector of hit times
            vector<double> vvPeakIntPowers; // vector of peak values
            int numSearchPeaks=2; //only look for two peaks  //TODO:  Maybe make this a console argument?  Though the cropping is based on having two peaks 4/10/2024
            double peakSeparation=50.0;  //Minimum separation between peaks.  Closest seperation is expected to be ~100 ns.
            getAbsMaximum_N(grInt, numSearchPeaks, peakSeparation, vvHitTimes, vvPeakIntPowers);
            double firstPeak;
            double secondPeak;
            if(vvHitTimes[1]>vvHitTimes[0]) {
                firstPeak = vvHitTimes[0];
                secondPeak = vvHitTimes[1];
            } 
            else {
                firstPeak = vvHitTimes[1];
                secondPeak = vvHitTimes[0];
            }
            cutoffTime[i] = firstPeak + (secondPeak-firstPeak)/2;
            //End double peak finder stuff
            
            //Crop waveform to first peak
            // grInt = FFTtools::cropWave(grInt, grInt->GetX()[0], cutoffTime[i]);
            
            //Test cropping waveform to 1024
            grInt = FFTtools::cropWave(grInt, grInt->GetX()[0], grInt->GetX()[1024-1]);  
            
            //Testing padding waveform to 2048
            // grInt = FFTtools::padWaveToLength(grInt, 2048);

            
            cout << "interpolated waveform length = " << grInt->GetN() << endl;
            // TGraph *grInt;
            // if (grCrop->GetN() < 2048) {
            //     grInt = FFTtools::padWaveToLength(grCrop, 2048);
            // }
            // else {
            //     grInt = FFTtools::padWaveToLength(grCrop, 4096);
            // }
            
            //Adding dynamic waveform padding to get length for factor of 2 for FFT purposes. - JCF 4/11/2024
            // int n=1; //Initial power of two
            // do {
            //     // cout << "Comparing waveform against " << TMath::Power(2,n) << endl;
            //     if (grInt->GetN() == TMath::Power(2,n)) {
            //         // cout << "Waveform already a power of two in length." << endl;
            //         n++;
            //         continue;
            //     }
            //     else if (grInt->GetN() < TMath::Power(2,n)) {
            //         cout << "Padding waveform to length of " << TMath::Power(2,n) << endl;
            //         grInt = FFTtools::padWaveToLength(grInt, TMath::Power(2,n));
            //     }
            //     else {
            //         n++;
            //         // cout << "Increasing iterator to " << TMath::Power(2,n) << endl;
            //     }
            // }
            // while (grInt->GetN() > TMath::Power(2,n-1));
            
            cout << "padded waveform length = " << grInt->GetN() << endl;
            interpolatedWaveforms[i] = grInt;
            
            //Debugging - plotting waveforms
            if (i == 0) {
                g[i] = grInt;
                mg->Add(g[i]);
                mg->Draw("AL");
                char title[500];
                sprintf(title, "waveform%d.png", i);
                c->Print(title);
            }
            // TLine *l1 = new TLine(firstPeak, -1000, firstPeak, 1000);
            // l1->Draw();            
            
            // Debugging drawing the waveform after interpolation - JCF 6/7/2023
            // TCanvas *c = new TCanvas("","", 1200, 950);            
            // grInt->Draw("colz");
            // grInt->GetXaxis()->SetTitle("Voltage [mV]");
            // grInt->GetYaxis()->SetTitle("Time [ns]");
            // grInt->GetYaxis()->SetTitleSize(0.05);
            // grInt->GetYaxis()->SetLabelSize(0.03);
            // grInt->GetYaxis()->SetTitleOffset(0.6);
            // grInt->GetXaxis()->SetTitleSize(0.05);
            // grInt->GetXaxis()->SetLabelSize(0.03);
            // grInt->GetXaxis()->SetTitleOffset(0.6);
            // gStyle->SetOptStat(0);
            // grInt->GetXaxis()->CenterTitle();
            // grInt->GetYaxis()->CenterTitle();
            // gPad->SetRightMargin(0.15);
            // char title[500];
            // sprintf(title,"maps_ev%d_channel%d.png", event, i);
            // c->SaveAs(title);
            // delete c;
            // End debugging            
            // delete gr, grCrop, grInt;
            delete gr, grInt;
        }
        
        //Debugging - plotting waveforms
        //c->SaveAs("testplot.png");
        // mg->Draw("AL");
        // c->Print("waveform.png");
        //mg->SaveAs("testplot.png");        

        //TODO: Use the data-driven noise model that AraSim has implemented 4/9/2024
        //Adding noise calculation on a per channel basis rather than hard-coding a noise value - JCF 7/24/2023
        // double noise = 46.; // the noise is basically 45 mV (independent of channel) in MC
        double noise; // Initializing noise for calculation from waveform below:
        // double noise = 0.; // Debugging and testing for noiseless waveforms. - JCF 6/25/2023
        std::map<int, double> snrs; // map of channels to SNRs
        for(int i=0; i<16; i++){
            //Calculate noise from rms of first 50 ns of waveform
            double voltageSubset[100];
            for(int j=0; j<100; j++){
                voltageSubset[j] = (interpolatedWaveforms[i]->GetY())[j];
            }
            noise = TMath::RMS(100,voltageSubset);
            
            double peak_max = TMath::MaxElement(interpolatedWaveforms[i]->GetN(), interpolatedWaveforms[i]->GetY());
            // cout << "peak_max = " << peak_max << endl;
            double peak_min = abs(TMath::MinElement(interpolatedWaveforms[i]->GetN(), interpolatedWaveforms[i]->GetY()));
            // cout << "peak_min = " << peak_min << endl;
            if(peak_min > peak_max){
                snrs[i] = peak_min/noise;
            }
            else{
                snrs[i] = peak_max/noise;
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

        // int likely_sol = guess_triggering_solution(eventPtr, reportPtr);
        // int likely_sol = 0; //Forcing direct solutions only for debugging - JCF 7/5/2023

        // stash the output
        fpOut->cd();
        // reconstructed quantities first
        for(int r=0 ;r<numScanned; r++){
            peakCorrs_out[r] = peakCorrs[r];
            peakThetas_out[r] = peakThetas[r];//*-1 + 90;  //In degrees, but with zero at the horizontal.  Changing to match zero at vertical. JCF 4/11/2024
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

        // likelySol_out = likely_sol;
        
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