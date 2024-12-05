/*
deconvolveWaveform.cpp
Justin Flaherty
4/6/2024

A simple deconvoultion script that inverts the antenna and electronics response applied in AraSim to convert a waveform in voltage into the electric field.

Requirements:
An installation of AraRoot and AraSim
Symlinking of the AraSim/data directory into the directory of this script.  Wherever you have this script hosted, do the following:

ln -s /path/to/AraSim/data .


*/



#include "/cvmfs/ara.opensciencegrid.org/trunk/centos7/source/libRootFftwWrapper/include/FFTtools.h"

// ROOT includes
#include "TFile.h"
#include "TRandom3.h" 
#include "TTree.h"
#include "TLatex.h"
#include "TLine.h"

// AraSim includes
//vector and position must be first
#include "Vector.h"
#include "Position.h"

#include "AraGeomTool.h"
#include "Constants.h"
#include "counting.hh"
#include "Detector.h"
#include "EarthModel.h"
#include "Efficiencies.h"
#include "Event.h"
#include "IceModel.h"
#include "Primaries.h"
#include "Ray.h"
#include "Report.h"
#include "RaySolver.h"
#include "secondaries.hh"
#include "Settings.h"
#include "signal.hh"
#include "Spectra.h"
#include "Tools.h"
#include "Trigger.h"



using namespace std;

// #ifdef ARA_UTIL_EXISTS
#include "UsefulIcrrStationEvent.h"
ClassImp(UsefulIcrrStationEvent);
#include "UsefulAtriStationEvent.h"
ClassImp(UsefulAtriStationEvent);
// #endif

//Including my tools file.  Should update name to avoid confusion with Tools.h above.
#include "polRecoTools.h"

//TODO: Have outgoing pointer equal the incoming pointer by using a filler function to copy the original information, then I replace the voltage info with my deconvolved voltage.
UsefulAtriStationEvent *usefulAtriEvPtr;
UsefulAtriStationEvent *usefulAtriEvPtrOut;
UsefulAtriStationEvent *usefulAtriCswPtrOut;

//Initializing some keyword arguments
bool debugMode = false;
bool separateCsw = true;
bool useMCTruth = false;
bool forceSinglePeak = false;
bool weinerCorrection = true;
bool spectralSubtract = false;
bool findPolarity = false;
int debugEvent=0;
double minCorr = 0;
double sampleNs=80;  //Sample window for noise and signal in the PSD SNR calculation.
double toleranceNs=30;
double snrThreshold=4;
// double dt=0.1; //Sample rate for interpolation

int main(int argc, char **argv)
{
    if(argc<8) {
        std::cout << "Usage\n" << argv[0] << " <runnum> <band-pass minimum frequency (MHz)> <band-pass maximum frequency (MHz)> <setup file> <input root file> <input reco file> <output directory>\n";
        std::cout << "e.g.\n" << argv[0] << " 1000 150 300 setup.txt AraOut.root recangle_out_run1000.root output/\n";        
        return 0;
    }
    
    // if (argc>8){
    //     std::string action(argv[8]);
    //     if(action == "debug") {
    //         debugMode = true;
    //     }
    // }     
    if (argc>8){
        int i = 8;
        while(argv[i] != NULL) {
            std::string action(argv[i]);
            if (action == "forceSinglePeak") {
                cout << "Enforcing single peak condition." << endl;
                forceSinglePeak = true;
            }         
            if (action == "findPolarity") {
                cout << "Enforcing find polarity condition." << endl;
                findPolarity = true;
            }              
            if (action == "debug") {
                cout << "Entering debug mode." << endl;
                debugMode = true;
                std::string event(argv[i+1]);
                if (atoi(argv[i+1])) {
                    debugEvent = atoi(argv[i+1]);
                }
            }
               
            i++;
        }
    }     
    
    //Import argument parameters
    char* runNumber = argv[1];
    double freqMin = atof(argv[2])*1e6;
    double freqMax = atof(argv[3])*1e6;
    char* setupfile = argv[4];
    char* araFile = argv[5];
    char* recoFile = argv[6];
    char* outputDir = argv[7];
    
    //Import AraRoot file
    printf("Opening root file...\n");
    TFile *fp = TFile::Open(araFile);
    if(!fp) { std::cerr << "Can't open file\n"; return -1; }
    printf("Root File opened!\n");
    
    //Import eventTree
    TTree *eventTree = (TTree*) fp->Get("eventTree");
    printf("Event tree opened!\n");
    if(!eventTree) { std::cerr << "Can't find eventTree\n"; return -1; }
    Long64_t numEntries=eventTree->GetEntries();
    cout << "eventTree has " << numEntries << " entries." << endl;
    RawAtriStationEvent *rawAtriEvPtr=0;
    
    // Check if sim or real data file by checking for existence of AraTree
    TTree *simSettingsTree;
    simSettingsTree=(TTree*) fp->Get("AraTree");
    TTree *AraTree2;
    AraTree2 = (TTree*) fp->Get("AraTree2");
    bool dataLike = false;
    bool calibrated;
    
    AraGeomTool *geomTool = AraGeomTool::Instance();
    Event *eventPtr = new Event();
    // Report *report = new Report();
    double weight;

    
    //Trying condition where it checks for usefulAtriStation branch and imports according to that.
    if(!simSettingsTree) { 
        dataLike = true;            
        std::cerr << "Importing as real data.\n";
        // TTree* atriExists=(TTree*) fp->Get("UsefulAtriStationEvent");
        // if (!atriExists) {
        //Checks if usefulAtri branch exists.  If not, fata gets imported as uncalibrated.
        if (not eventTree->GetBranch("UsefulAtriStationEvent")) {
            calibrated = false;
            cout << "Importing as uncalibrated." << endl;
            eventTree->SetBranchAddress("event",&rawAtriEvPtr);
        }
        else {
            calibrated = true;
            eventTree->SetBranchAddress("UsefulAtriStationEvent", &usefulAtriEvPtr);
        }
        weight = 1;
    }
    // sim like
    else {
        dataLike = false;
        calibrated = true;
        std::cerr << "AraTree exists.  Importing as simulated data.\n";
        eventTree->SetBranchAddress("UsefulAtriStationEvent", &usefulAtriEvPtr);
        // double weight;
        eventTree->SetBranchAddress("weight", &weight);   
        AraTree2->SetBranchAddress("event", &eventPtr);
        // AraTree2->SetBranchAddress("report", &report);
    }
    //End new import method    
        
    //Import vertex reco file
    printf("Opening reco file...\n");
    TFile *fp2 = TFile::Open(recoFile);
    if(!fp2) { std::cerr << "Can't open file\n"; return -1; }
    printf("Reco File opened!\n");
    TTree *vertexReco = (TTree*) fp2->Get("vertexReco");
    double reco_arrivalThetas[16];
    double reco_arrivalPhis[16];
    double cutoffTime[16];
    double cutoffTimeCsw[16];
    double arrivalTimes[16];
    double directPeakTimes[16];
    double refPeakTimes[16];   
    double psiReco[8];  
    double psiRecoMean;
    double psiRecoCsw[8];  
    double psiRecoCsw2[8];  
    double psiRecoCsw3[8];  
    double psiRecoCsw4[8];  
    double psiTrue[8];
    double psiTrue2[8];
    double psiTrue3[8];
    double psiTrue4[8];
    double bestCorr;
    double hilbertPeakOut[16];
    double hilbertPeakCswOut[16];
    double peakTimeOut[16];
    double peakTimeCswOut[16];
    double launch_theta[16];
    double launch_phi[16];
    double theta_nutraject[8];
    double theta_nutrajectCsw[8];
    double phi_nutraject[8];
    double phi_nutrajectCsw[8];
    double vertexRadius;
    double vertexTheta;
    double vertexPhi;
    int bestSol;
    int bestSol_out;
    double weight_out;
    
    double true_arrivalThetas[16];
    double true_arrivalPhis[16];
    double true_launch_theta[16];
    double true_launch_phi[16];   
    double trueVertexRadius;
    double trueVertexTheta;
    double trueVertexPhi;
    double true_theta_nutraject[8];
    double true_phi_nutraject[8];    
    double true_theta_nutraject2[8];
    double true_phi_nutraject2[8];
    double true_theta_nutraject3[8];
    double true_phi_nutraject3[8];
    double true_theta_nutraject4[8];
    double true_phi_nutraject4[8]; 
    
    double mean_theta_nutraject;
    double mean_phi_nutraject;

    double v_snr;
    double h_snr;
    double snr[16];
    
    double true_viewingAngle[8];
    double cherenkovAngle;
    
    // double weight;
    
    // Testing using the true rf angles
    if (useMCTruth){
        vertexReco->SetBranchAddress("true_arrivalThetas", reco_arrivalThetas);
        vertexReco->SetBranchAddress("true_arrivalPhis", reco_arrivalPhis);  
        vertexReco->SetBranchAddress("true_launchThetas", launch_theta);
        vertexReco->SetBranchAddress("true_launchPhis", launch_phi);
        vertexReco->SetBranchAddress("trueR", &vertexRadius);
        vertexReco->SetBranchAddress("trueTheta", &vertexTheta);
        vertexReco->SetBranchAddress("truePhi", &vertexPhi);           
    }
    // end testing
    else {
        vertexReco->SetBranchAddress("reco_arrivalThetas", reco_arrivalThetas);
        vertexReco->SetBranchAddress("reco_arrivalPhis", reco_arrivalPhis);
        vertexReco->SetBranchAddress("reco_launchThetas", launch_theta);
        vertexReco->SetBranchAddress("reco_launchPhis", launch_phi);
        vertexReco->SetBranchAddress("bestR", &vertexRadius);
        vertexReco->SetBranchAddress("bestTheta", &vertexTheta);
        vertexReco->SetBranchAddress("bestPhi", &vertexPhi);        
    }
    vertexReco->SetBranchAddress("cutoffTime", cutoffTime);
    vertexReco->SetBranchAddress("arrivalTimes", arrivalTimes);
    // vertexReco->SetBranchAddress("excludedChannels", &excludedChannels);    
    vertexReco->SetBranchAddress("v_snr", &v_snr);
    vertexReco->SetBranchAddress("h_snr", &h_snr);
    vertexReco->SetBranchAddress("snrs", &snr);
    vertexReco->SetBranchAddress("directPeakTimes", &directPeakTimes);  
    vertexReco->SetBranchAddress("refPeakTimes", &refPeakTimes);   
    vertexReco->SetBranchAddress("bestCorr", &bestCorr);
    vertexReco->SetBranchAddress("bestSol", &bestSol);
    // vertexReco->SetBranchAddress("weight", &weight);
    
    if (not dataLike) {
        vertexReco->SetBranchAddress("true_arrivalThetas", true_arrivalThetas);
        vertexReco->SetBranchAddress("true_arrivalPhis", true_arrivalPhis);  
        vertexReco->SetBranchAddress("true_launchThetas", true_launch_theta);
        vertexReco->SetBranchAddress("true_launchPhis", true_launch_phi);
        vertexReco->SetBranchAddress("trueR", &trueVertexRadius);
        vertexReco->SetBranchAddress("trueTheta", &trueVertexTheta);
        vertexReco->SetBranchAddress("truePhi", &trueVertexPhi);          
    }

    Long64_t numEntriesVertex=vertexReco->GetEntries();
    cout << "Vertex Reco tree opened! Has " << numEntriesVertex << " entries!" << endl;;
    
    printf("------------------\n");
    printf("Input files loaded.  Setting up detector stuff.\n");
    printf("------------------\n");
    
    // string setupfile;
    // setupfile = argv[2];
    Settings *settings1 = new Settings();
    settings1->ReadFile(setupfile); 
    
    AraEventCalibrator *cal = AraEventCalibrator::Instance();
    if (dataLike) {
        setPedestalFile(cal, settings1->DETECTOR_STATION, runNumber);
    }
    
    IceModel *icemodel=new IceModel(settings1->ICE_MODEL + settings1->NOFZ*10,settings1->CONSTANTICETHICKNESS * 1000 + settings1->CONSTANTCRUST * 100 + settings1->FIXEDELEVATION * 10 + 0,settings1->MOOREBAY);// creates Antarctica ice model
    Detector *detector = new Detector(settings1, icemodel, setupfile);  
    // Report *report = new Report(detector, settings1);
    Report *report = new Report(detector, settings1);
    // Event *eventPtr = new Event();
    
    //Use the getTrigMasking function to use the same channels that triggering used for the reconstruction
    std::vector<int> excludedChannels; 
    getExcludedChannels(excludedChannels, settings1, detector);

    //dt interval for waveform interpolation.  Matching to the timestep specified in the setup file.
    double dt = settings1->TIMESTEP*1e9;
    
    printf("------------------\n");
    printf("Make Output Files\n");
    printf("------------------\n");

    char outfile_name[400];
    sprintf(outfile_name, "%s/deconvolvedWaveforms_run_%s.root", outputDir, runNumber);    
    
    std::cout<<"Output name is "<<outfile_name<<std::endl;
    TFile *fpOut = TFile::Open(outfile_name, "RECREATE");
    if(!fpOut){ std::cerr<<"Cannot open output file "<<fpOut<<std::endl; return -1; }    
    
    // TTree *outTree = new TTree("eventTree", "eventTree");
    TTree *outTree = eventTree->CloneTree(0);
    outTree->Branch("UsefulAtriStationEvent", &usefulAtriEvPtrOut);
    outTree->Branch("bestCorr", &bestCorr, "bestCorr/D");
    outTree->Branch("psiReco", &psiReco, "psiReco[8]/D");
    outTree->Branch("psiTrue", &psiTrue, "psiTrue[8]/D");
    outTree->Branch("theta_nutraject", &theta_nutraject, "theta_nutraject[8]/D");
    outTree->Branch("phi_nutraject", &phi_nutraject, "phi_nutraject[8]/D");
    outTree->Branch("bestSol", &bestSol_out, "bestSol/I");
    outTree->Branch("weight", &weight_out, "weight/D");
    // outTree->Branch("pulserDepth", &pulserDepth, "pulserDepth/D")
    
    TTree *outTreeCsw = new TTree("coherentSum", "coherentSum");
    outTreeCsw->Branch("UsefulAtriStationEvent", &usefulAtriCswPtrOut); 
    outTreeCsw->Branch("bestCorr", &bestCorr, "bestCorr/D");
    outTreeCsw->Branch("psiReco", &psiRecoCsw, "psiReco[8]/D");
    outTreeCsw->Branch("psiReco2", &psiRecoCsw2, "psiReco2[8]/D");
    outTreeCsw->Branch("psiReco3", &psiRecoCsw3, "psiReco3[8]/D");
    outTreeCsw->Branch("psiReco4", &psiRecoCsw4, "psiReco4[8]/D");
    outTreeCsw->Branch("theta_nutraject", &theta_nutrajectCsw, "theta_nutraject[8]/D");
    outTreeCsw->Branch("phi_nutraject", &phi_nutrajectCsw, "phi_nutraject[8]/D");
    outTreeCsw->Branch("mean_theta_nutraject", &mean_theta_nutraject, "mean_theta_nutraject/D");
    outTreeCsw->Branch("mean_phi_nutraject", &mean_phi_nutraject, "mean_phi_nutraject/D");
    outTreeCsw->Branch("cherenkovAngle", &cherenkovAngle, "cherenkovAngle/D");
    // outTreeCsw->Branch("pulserDepth", &pulserDepth, "pulserDepth/D")
    
    // TTree *outTreeVertexReco = new TTree("vertexReco", "vertexReco");
    TTree *outTreeVertexReco = vertexReco->CloneTree();
    TTree *outTreeMCTruth;
    if (not dataLike) {
        // outTreeMCTruth = new TTree("AraTree2", "AraTree2");   
        outTreeMCTruth = AraTree2->CloneTree();
        
        outTreeCsw->Branch("true_theta_nutraject", &true_theta_nutraject, "true_theta_nutraject[8]/D");
        outTreeCsw->Branch("true_theta_nutraject2", &true_theta_nutraject2, "true_theta_nutraject2[8]/D");
        outTreeCsw->Branch("true_theta_nutraject3", &true_theta_nutraject3, "true_theta_nutraject3[8]/D");
        outTreeCsw->Branch("true_theta_nutraject4", &true_theta_nutraject4, "true_theta_nutraject4[8]/D");   
        
        outTreeCsw->Branch("true_phi_nutraject", &true_phi_nutraject, "true_phi_nutraject[8]/D");
        outTreeCsw->Branch("true_phi_nutraject2", &true_phi_nutraject2, "true_phi_nutraject2[8]/D");
        outTreeCsw->Branch("true_phi_nutraject3", &true_phi_nutraject3, "true_phi_nutraject3[8]/D");
        outTreeCsw->Branch("true_phi_nutraject4", &true_phi_nutraject4, "true_phi_nutraject4[8]/D");  
        outTreeCsw->Branch("true_viewingAngle", &true_viewingAngle, "true_viewingAngle[8]/D");
    }
    
    //Need to grab lengths of voltage and time arrays from eventTree to initialize the branches in the outfile.
    Int_t fNumChannels; ///< The number of channels
    std::map< Int_t, std::vector <Double_t> > fTimesOut; ///< The times of samples
    std::map< Int_t, std::vector <Double_t> > fVoltsOut; ///< The voltages of samples   
    
    //Initialize voltage and time arrays for the coherent sum.
    Int_t fNumChannelsCsw = 2; ///< The number of channels (one for Vpol, one for Hpol)
    std::map< Int_t, std::vector <Double_t> > fTimesCswOut; ///< The times of samples
    std::map< Int_t, std::vector <Double_t> > fVoltsCswOut; ///< The voltages of samples      
     
    
    //Loop over events
    int loopStart;
    int loopEnd;
    if (debugMode) {
        loopStart=debugEvent;
        loopEnd=loopStart+1;
    }
    else {
        loopStart=0;
        loopEnd=numEntries;
    }
    
    for(Long64_t event=loopStart;event<loopEnd;event++) {
        double pulserDepth;
        std::cout<<"Looking at event number "<<event<<std::endl;
        fp->cd();

        eventTree->GetEntry(event);
        if (not dataLike) {
            AraTree2->GetEntry(event);
        }
        
        cout << "Importing vertex reco info" << endl;
        vertexReco->GetEntry(event);
        bestSol_out = bestSol;
        weight_out = weight;
        
        //Adding peak correlation cut
        if (bestCorr < minCorr) {
            cout << "Event correlation of " << bestCorr << " is below threshold of " << minCorr << ". Bypassing event." << endl;
            continue;
        }
        
        if (debugMode){cout<<"bestCorr = " << bestCorr << endl;}
        
        if (not calibrated) {
            cout << "Triggering datalike condition." << endl;
            delete usefulAtriEvPtr;         
            usefulAtriEvPtr = new UsefulAtriStationEvent(rawAtriEvPtr, AraCalType::kLatestCalib);
            
            
        }
        
        // outTree = eventTree->CloneTree();  
        // outTree->Fill(); 
        
        // outTreeVertexReco = vertexReco->CloneTree();
        // outTreeVertexReco->Fill();
        
        try {
            pulserDepth = getPulserDepth(usefulAtriEvPtr->unixTime);
            // pulserDepth = getPulserDepth(usefulAtriEvPtrOut->unixTime);
            cout << "Pulser depth at " << pulserDepth << endl;
        }
        catch(std::out_of_range) {
            cout << "Pulser depth at " << runNumber << endl;
        }         
        
        
        //Need this for importing the vertexReco from simulated data, as simulated data gets stored by string, then antenna.
        int vertexRecoElectToRFChan[] = {14,2,6,10,  //VPols
                                         12,0,4,8,   //VPols
                                         15,3,7,11,  //HPols
                                         13,1,5,9};  //HPols

        for(int i=0; i<16; i++){
            // usefulAtriEvPtr->stationId = settings1->DETECTOR_STATION_ARAROOT;
            if (debugMode) {
                cout << "########################################" << endl;
                cout << "Channel = " << i << endl;
                cout << "usefulAtriEvPtr->stationId = " << usefulAtriEvPtr->stationId << endl;
                cout << "settings1->DETECTOR_STATION_ARAROOT = " << settings1->DETECTOR_STATION_ARAROOT << endl;
                cout << "usefulAtriEvPtr->getStationId() = " << usefulAtriEvPtr->getStationId() << endl;
            }

            TGraph *gr = usefulAtriEvPtr->getGraphFromRFChan(i);
            TGraph *grOriginal = gr;

            //Save initial and final time for truncating the padded arrays before output.
            double timeStart = gr->GetX()[0];
            double timeEnd = gr->GetX()[gr->GetN()-1];
            
            if (debugMode) {
                cout << "timeStart = " << timeStart << endl;
                cout << "timeEnd = " << timeEnd << endl;
            }
            
            //Grabbing noise sample of waveform to try spectral subtraction
            //Check if peakTime is in the sample region at beginning of waveform
            double tNoiseMin;
            double tNoiseMax;  
            if (forceSinglePeak) {
                refPeakTimes[i] = directPeakTimes[i];
            }
            if (directPeakTimes[i] < gr->GetX()[int(sampleNs/dt)]) {
                
                tNoiseMax = timeEnd-dt;
                tNoiseMin = tNoiseMax - sampleNs;                
            }
            else {              
                tNoiseMin = timeStart+dt;
                tNoiseMax = tNoiseMin + sampleNs;                

            }   
            
            //Apply time shift to center waveform
            double timeShift = (timeStart+timeEnd)/2;
            double initialWaveformLength = gr->GetN();
            for(int k=0; k<initialWaveformLength; k++){
                gr->GetX()[k] = gr->GetX()[k] - timeShift;
            }               
            
            // if (debugMode) {cout << "rrr" << endl;}
            //Interpolate graph to 0.5 ns resolution (or whatever TIMESTEP is specified in the setup file)
            gr = FFTtools::getInterpolatedGraph(gr,dt);
            grOriginal = FFTtools::getInterpolatedGraph(grOriginal,dt);
            if (debugMode) {
                cout << "gr->GetN() = " << gr->GetN() << endl;
                cout << "gr->GetX()[0] = " << gr->GetX()[0] << endl;
            }
            TGraph *grNoise = FFTtools::cropWave(gr, tNoiseMin-timeShift, tNoiseMax-timeShift);

            double noiseSampleLength = grNoise->GetN();

            //Pad waveform to a factor of two as specified in the setup file. - JCF 9/27/2023
            gr = resizeForFFT(gr, settings1->NFOUR);

            if (debugMode) {
                cout << "gr->GetN() = " << gr->GetN() << endl;
                cout << "gr->GetX()[0] = " << gr->GetX()[0] << endl;
            }            
            grNoise = resizeForFFT(grNoise, settings1->NFOUR);

            // Get length of padded waveform
            int waveform_bin = gr->GetN();            
            
            
            double heff_lastbin;
            double freq_lastbin;
            double time[waveform_bin];
            double voltage[waveform_bin];
            double timeNoise[waveform_bin];
            double voltageNoise[waveform_bin];            
            double volts_forint[settings1->NFOUR / 2];
            double T_forint[settings1->NFOUR / 2];
            
            //TODO: This init_T isn't dynamic to the imported data.  Should make have it defined based on the input waveform.
            double init_T = settings1->TIMESTEP *-1.e9 *((double) settings1->NFOUR / 4);
            for(int k=0; k<waveform_bin; k++){
                time[k] = gr->GetX()[k];
                voltage[k] = gr->GetY()[k];
            }            
            
            for(int k=0; k<waveform_bin; k++){
                timeNoise[k] = grNoise->GetX()[k];
                voltageNoise[k] = grNoise->GetY()[k];
            }              
            // delete gr;
            // cout << "time[0] = " << time[0] << endl;
            
            for (int m = 0; m < settings1->NFOUR / 2; m++)
            {
                T_forint[m] = init_T + m*dt;   // in ns
            }
            
            //Importing the cutoff time between spicecore peaks
            double cutoffTimeChannel;
               
            if (!cutoffTime or forceSinglePeak) {
                cout << "Using single-peak condition" << endl;
                cutoffTime[i] = time[waveform_bin-1];             
            }        
    
            double freq_tmp, heff, antenna_theta, antenna_phi;  // values needed for apply antenna gain factor and prepare fft, trigger
            
            
            //TODO:  Verify this is grabbing the correct channel in simulated versus real data (both calibrated and un calibrated)
            if (dataLike) {     
                //Import RF angles use RF channel mapping
                antenna_theta = reco_arrivalThetas[i];//*180/PI;
                antenna_phi = reco_arrivalPhis[i];//*180/PI;
            }
            else {
                //Import RF angles using electric channel mapping
                antenna_theta = reco_arrivalThetas[vertexRecoElectToRFChan[i]];//*180/PI;
                antenna_phi = reco_arrivalPhis[vertexRecoElectToRFChan[i]];//*180/PI;              
            }
            
            if (debugMode) {
                cout << "antenna_theta = " << antenna_theta << endl;
                cout << "antenna_phi = " << antenna_phi << endl;
            }
            
            //Calculate polarization vector that inverts the polarization factor (makes dot products equal to one)
            double newPol_vectorX = -sin(antenna_phi*PI/180);
            double newPol_vectorY = cos(antenna_phi*PI/180);
            double newPol_vectorZ = -1/sin(antenna_theta*PI/180);       

            Vector Pol_vector = Vector(newPol_vectorX, newPol_vectorY, newPol_vectorZ);

            double dF_Nnew;

            double nice = 1.79;  //TODO: Have this use the n(z) in the icemodel. 4/30/2023


            int pol_ant;
            int gain_ch_no = i;          
            double Pol_factor;
            
            if (i < 8) {
                pol_ant=0;
            } else {
                pol_ant=1;
            }
            
            double dT_forfft = time[1] - time[0]; 
            
            int Nnew = createFFTsize(dT_forfft, settings1);
            
            //Initialize and populate the voltage and time arrays for the fft.
            double V_forfft[Nnew];
            double T_forfft[Nnew];
            
            createFFTarrays(voltage, time, settings1, V_forfft, T_forfft, waveform_bin, Nnew);
            
            //Testing taking FFT of noise sample and subtracting from signal.
            double Vnoise_forfft[Nnew];
            double Tnoise_forfft[Nnew];           
            createFFTarrays(voltageNoise, timeNoise, settings1, Vnoise_forfft, Tnoise_forfft, waveform_bin, Nnew);            
            
            //Testing get PSD SNR function
            
            double snrPsd[int(Nnew/2)];
            double freqPsd[int(Nnew/2)];
            getPowerSpectrumSNR(grOriginal, tNoiseMin-timeShift, tNoiseMax-timeShift, directPeakTimes[i], refPeakTimes[i], waveform_bin, Nnew, settings1, snrPsd, freqPsd, dt, sampleNs, debugMode);

            delete gr;
            
            // get spectrum with zero padded WF
            Tools::realft(V_forfft, 1, Nnew); 
            Tools::realft(Vnoise_forfft, 1, Nnew); 
            
            if (debugMode) {
                cout << "*****************************" << endl;
                cout << "V_forfft (before factors) = " << endl;
                // for (int n=0; n<Nnew/2; n++) {
                for (int n=0; n<10; n++) {
                    cout << V_forfft[2*n] << " + i " << V_forfft[2*n+1] << ", ";
                }
                cout << endl;         
            }                     
            
            dF_Nnew = 1. / ((double)(Nnew) *(dT_forfft) *1.e-9);    // in Hz

            freq_tmp = dF_Nnew *((double) Nnew / 2. + 0.5); // in Hz 0.5 to place the middle of the bin and avoid zero freq 

            freq_lastbin = freq_tmp;
            
            heff_lastbin = report->GaintoHeight(detector->GetGain_1D_OutZero(freq_tmp *1.E-6,   // to MHz
                                                antenna_theta, antenna_phi, pol_ant),
                                                freq_tmp, nice);                  
            if (debugMode) {cout << "aaa" << endl;}
            for (int n = 0; n < Nnew / 2; n++)
            // for (int n = 0; n < settings1->NFOUR / 2; n++)            
            {
                //Generate correction factors for the Weiner deconvolution.
                double realWeinerCorr = 1;
                double imagWeinerCorr = 1;               
                
                freq_tmp = dF_Nnew *((double) n + 0.5); // in Hz 0.5 to place the middle of the bin and avoid zero freq          

                heff = report->GaintoHeight(detector->GetGain_1D_OutZero(freq_tmp *1.E-6,   // to MHz
                                            antenna_theta, antenna_phi, pol_ant),
                                            freq_tmp, nice);

                if (spectralSubtract) {
                    V_forfft[2*n] -= Vnoise_forfft[2*n];                    
                }
                // invert entire elect chain gain, phase
                if (n > 0)
                {  
                if (debugMode and n == 500) {
                    cout << "Before Elec " << endl;                    
                    cout <<  "V_forfft[2*n] = " <<  V_forfft[2*n] << endl;
                    cout <<  "V_forfft[2*n+1] = " <<  V_forfft[2*n+1] << endl;      
                }                           
                    report->InvertElect_Tdomain(
                        freq_tmp *1.e-6, 
                        detector, 
                        V_forfft[2 *n], 
                        V_forfft[2 *n + 1], 
                        gain_ch_no,
                        settings1);
                    
                    //Weiner Deconvolution Correction
                    report->InvertElect_Tdomain(
                        freq_tmp *1.e-6, 
                        detector, 
                        realWeinerCorr, 
                        imagWeinerCorr, 
                        gain_ch_no,
                        settings1);                    
                    
                }
                // if (debugMode) {cout << "ccc" << endl;}
                else
                {
                    report->InvertElect_Tdomain_FirstTwo(
                        freq_tmp *1.e-6, 
                        freq_lastbin *1.e-6, 
                        detector, 
                        V_forfft[2 *n], 
                        V_forfft[2 *n + 1], 
                        gain_ch_no,
                        settings1);
                    
                    //Weiner Deconvolution Correction
                    report->InvertElect_Tdomain_FirstTwo(
                        freq_tmp *1.e-6, 
                        freq_lastbin *1.e-6, 
                        detector, 
                        realWeinerCorr, 
                        imagWeinerCorr, 
                        gain_ch_no,
                        settings1);                    
                }
                if (debugMode and n == 500) {
                    cout << "After Elec, before Ant "<< endl;                    
                    cout <<  "V_forfft[2*n] = " <<  V_forfft[2*n] << endl;
                    cout <<  "V_forfft[2*n+1] = " <<  V_forfft[2*n+1] << endl;      
                }    
                // if (debugMode) {cout << "ddd" << endl;}
                if (n > 0)
                {
                    report->InvertAntFactors_Tdomain(
                        detector->GetAntPhase_1D(freq_tmp *1.e-6, antenna_theta, antenna_phi, pol_ant),
                        heff, 
                        Pol_vector, 
                        pol_ant, 
                        Pol_factor, 
                        V_forfft[2 *n], 
                        V_forfft[2 *n + 1],
                        settings1,
                        antenna_theta, 
                        antenna_phi,
                        freq_tmp);   
                    
                    //Weiner Deconvolution Correction
                    report->InvertAntFactors_Tdomain(
                        detector->GetAntPhase_1D(freq_tmp *1.e-6, antenna_theta, antenna_phi, pol_ant),
                        heff, 
                        Pol_vector, 
                        pol_ant, 
                        Pol_factor, 
                        realWeinerCorr, 
                        imagWeinerCorr,
                        settings1,
                        antenna_theta, 
                        antenna_phi,
                        freq_tmp);                       
                }
                else
                {
                    report->InvertAntFactors_Tdomain_FirstTwo(
                        heff, 
                        heff_lastbin, 
                        Pol_vector, 
                        pol_ant, 
                        Pol_factor, 
                        V_forfft[2 *n], 
                        V_forfft[2 *n + 1],
                        antenna_theta, 
                        antenna_phi,
                        freq_tmp);
                    
                    //Weiner Correction
                    report->InvertAntFactors_Tdomain_FirstTwo(
                        heff, 
                        heff_lastbin, 
                        Pol_vector, 
                        pol_ant, 
                        Pol_factor, 
                        realWeinerCorr, 
                        imagWeinerCorr,
                        antenna_theta, 
                        antenna_phi,
                        freq_tmp);                    

                }
                
                if (debugMode and n == 500) {
                    cout << "After ant " << endl;
                    cout << "realWeinerCorr = " << realWeinerCorr << endl;
                    cout << "imagWeinerCorr = " << imagWeinerCorr << endl;                    
                    cout <<  "V_forfft[2*n] = " <<  V_forfft[2*n] << endl;
                    cout <<  "V_forfft[2*n+1] = " <<  V_forfft[2*n+1] << endl;
                    cout << "snrPsd[n] = " << snrPsd[n] << endl;
                }                
                
                if (weinerCorrection) {
                    //Correction for weiner deconvolution
                    wienerDeconvolution(V_forfft[2 *n], V_forfft[2 *n + 1], realWeinerCorr, imagWeinerCorr, snrPsd[n]);
                }
                
                else {
                    if (freq_tmp > 850.*1.e6 or freq_tmp < 100.*1.e6) {              
                        V_forfft[2*n] = 0;
                        V_forfft[2*n+1] = 0;
                    }  
                }
                
                if (debugMode and n == 500) {
                    cout << "After wiener " << endl;                  
                    cout <<  "V_forfft[2*n] = " <<  V_forfft[2*n] << endl;
                    cout <<  "V_forfft[2*n+1] = " <<  V_forfft[2*n+1] << endl;      
                }     
                
                applyBandpassBin(V_forfft[2*n], V_forfft[2*n+1], freq_tmp, freqMin, freqMax);          
                //End Butterworth filter 
            }   // end for freq bin
                            
            // now get time domain waveform back by inv fft              
            Tools::realft(V_forfft, -1, Nnew);
           
            Tools::SincInterpolation(Nnew, T_forfft, V_forfft, settings1->NFOUR / 2, T_forint, volts_forint);  

            //Now write deconvolved voltage data to file.
            for (int n = 0; n < settings1->NFOUR / 2; n++)
            {
                int elecChan;
                if (settings1->DETECTOR_STATION == 1) {
                    elecChan = geomTool->getElecChanFromRFChan(i, 100);
                }
                else {
                    elecChan = geomTool->getElecChanFromRFChan(i, settings1->DETECTOR_STATION_ARAROOT);
                }

                //Testing cropping the waveform padding
                if ((T_forint[n] > timeStart-timeShift) and (T_forint[n] < timeEnd-timeShift)) { 
                    usefulAtriEvPtrOut->fVolts[elecChan].push_back(settings1->ARBITRARY_EVENT_ATTENUATION *volts_forint[n] *2. / (double)(Nnew));  // 2/N for inverse FFT normalization factor
                    usefulAtriEvPtrOut->fTimes[elecChan].push_back(T_forint[n]+timeShift);    
                }
            }
            usefulAtriEvPtrOut->stationId = usefulAtriEvPtr->stationId;
            usefulAtriCswPtrOut->stationId = usefulAtriEvPtr->stationId;            
            
        } //channel loop
        usefulAtriEvPtrOut->eventNumber = usefulAtriEvPtr->eventNumber;
        usefulAtriEvPtrOut->unixTime = usefulAtriEvPtr->unixTime;
        // cout << "eee" << endl;
        
        //Calculate hilbert peaks and polarization angles.  
        peakHilbert(usefulAtriEvPtrOut, hilbertPeakOut, peakTimeOut, cutoffTime, dt, toleranceNs=toleranceNs, findPolarity=findPolarity, debugMode=debugMode);
        calculatePsi(hilbertPeakOut, psiReco, findPolarity);       
        double nofz_atVertex = getNofzAtVertex(geomTool, usefulAtriEvPtrOut, icemodel, vertexRadius, vertexTheta);
        cherenkovAngle = calculateCherenkovAngle(nofz_atVertex);
        double nofz_atTrueVertex;
        
        if (not dataLike) { 
            nofz_atTrueVertex = getNofzAtVertex(geomTool, usefulAtriEvPtrOut, icemodel, trueVertexRadius, trueVertexTheta);
        }
        
        // if (not dataLike) {
        //     cout << "asdasdasdasdasdsad " << calculateTruePsi(eventPtr, report, launch_theta[0]*PI/180, launch_phi[0]*PI/180) << endl;
        // }

        for (int i=0; i<8; i++) {
            calculateNuTrajectory(psiReco[i]*PI/180, launch_theta[i]*PI/180, launch_phi[i]*PI/180, cherenkovAngle, theta_nutraject[i], phi_nutraject[i]);
            // if (not dataLike) {
            //     psiTrue[i] = calculateTruePsi(eventPtr, report, launch_theta[i]*PI/180, launch_phi[i]*PI/180);
            //     // calculateNuTrajectory(psiTrue[i]*PI/180, launch_theta[i]*PI/180, launch_phi[i]*PI/180, nofz_atVertex, theta_nutraject[i], phi_nutraject[i]);
            //     cout << "theta_nutraject_true = " << eventPtr->Nu_Interaction[0].nnu.Theta()*180/PI << endl;
            //     cout << "phi_nutraject_true = " << eventPtr->Nu_Interaction[0].nnu.Phi()*180/PI << endl;
            // }       
            
        }
        
        //Assign timestamp values to help identify calpulsers
        usefulAtriEvPtrOut->timeStamp = usefulAtriEvPtr->timeStamp;
        // Assign triggerInfo values to identify RF and software triggers.
        for (int bit = 0; bit < 4; bit++) {
            usefulAtriEvPtrOut->triggerInfo[bit] = usefulAtriEvPtr->triggerInfo[bit];
        } 
        
        
        usefulAtriCswPtrOut->eventNumber = usefulAtriEvPtr->eventNumber;
        usefulAtriCswPtrOut->unixTime = usefulAtriEvPtr->unixTime;

        //Assign timestamp values to help identify calpulsers
        usefulAtriCswPtrOut->timeStamp = usefulAtriEvPtr->timeStamp;
        // Assign triggerInfo values to identify RF and software triggers.
        for (int bit = 0; bit < 4; bit++) {
            usefulAtriCswPtrOut->triggerInfo[bit] = usefulAtriEvPtr->triggerInfo[bit];
        }   

        //Plotting deconvolved Waveform before CSW
        
        //Create array of channel pairs to exclude in the coherent sum
        std::vector<int> cswExcludedChannelPairs;
        for(int i=0; i<8; i++){
            bool checkVpol = std::find(excludedChannels.begin(), excludedChannels.end(), i) != excludedChannels.end();
            bool checkHpol = std::find(excludedChannels.begin(), excludedChannels.end(), i+8) != excludedChannels.end();
            if (checkVpol or checkHpol or (snr[i] < snrThreshold) or (snr[i+8] < snrThreshold)) {
                cout << "Excluding channel pair " << i << endl;
                cout << "snr[i] = " << snr[i] << endl;
                cout << "snr[i+8] = " << snr[i+8] << endl;
                cswExcludedChannelPairs.push_back(i);
            }

        }        
        //TCanvas for the coherently summed waveform separated by VPol and HPol
        if (debugMode) {
            TCanvas *cTimeshift = new TCanvas("","", 1600, 1600);
            cTimeshift->Divide(4,4);          

//             //Create array of channel pairs to exclude in the coherent sum
//             std::vector<int> cswExcludedChannelPairs;
//             for(int i=0; i<8; i++){
//                 bool checkVpol = std::find(excludedChannels.begin(), excludedChannels.end(), i) != excludedChannels.end();
//                 bool checkHpol = std::find(excludedChannels.begin(), excludedChannels.end(), i+8) != excludedChannels.end();
//                 if (checkVpol or checkHpol) {
//                     cout << "Excluding channel pair " << i << endl;
//                     cswExcludedChannelPairs.push_back(i);
//                 }

//             }
            
            // cout << "bbb" << endl;
            
            double arrivalTimeMax = *max_element(arrivalTimes, arrivalTimes + 16);
            double arrivalTimeMin = *min_element(arrivalTimes, arrivalTimes + 16);        

            for(int i=0; i<16; i++){

                //Check if channel is in the excluded channel list
                bool checkExcluded = std::find(cswExcludedChannelPairs.begin(), cswExcludedChannelPairs.end(), i%8) != cswExcludedChannelPairs.end();

                TGraph *gr = usefulAtriEvPtrOut->getGraphFromRFChan(i);

                for(int k=0; k<gr->GetN(); k++){
                    gr->GetX()[k] = gr->GetX()[k];
                }

                gr = FFTtools::getInterpolatedGraph(gr,dt);

                cTimeshift->cd(i+1); gPad->SetGrid(1,1);
                gr->Draw();
                if (checkExcluded) {
                    gr->SetLineColor(2);
                }
                TLine *lPeak = new TLine(peakTimeOut[i], -1500, peakTimeOut[i], 1500);
                lPeak->SetLineColorAlpha(kBlue, 1);
                lPeak->Draw();     
                TLine *l6v = new TLine(cutoffTime[i], -1500, cutoffTime[i], 1500);
                l6v->SetLineStyle(kDashed);
                l6v->Draw();
                // TGraph *hilbert = FFTtools::getHilbertEnvelope(gr);
                // hilbert->Draw();       
                // hilbert->GetXaxis()->SetLimits(250, 350);   
                char vTitle[500];
                // sprintf(vTitle, "A%d Run %s Event %d Ch. %.2d %0.2f", settings1->DETECTOR_STATION_ARAROOT, runNumber, usefulAtriEvPtr->eventNumber, i, pulserDepth);
                // sprintf(vTitle, "%0.2f", peakTimeOut[i]);
                sprintf(vTitle, "(%0.4f, %0.4f)", peakTimeOut[i], hilbertPeakOut[i]);
                // sprintf(vTitle, "A%d Run %s Event %d Ch. %.2d %0.2f", settings1->DETECTOR_STATION_ARAROOT, runNumber, usefulAtriEvPtr->eventNumber, i, TMath::MaxElement(gr->GetN(),gr->GetY()));                
                gr->SetTitle(vTitle);  
                // gr->GetXaxis()->SetLimits(250, 350);   


            }


            char title[500];
            // sprintf(title, "%s/waveformDeconvolved_%s.png", outputDir, runNumber);
            sprintf(title, "%s/waveformDeconvolved.png", outputDir);
            cTimeshift->Print(title);  
        }
        
        //Plotting CSW deconvolved waveform
        //Create coherent sum logic
        TGraph *cswVpol = new TGraph();
        TGraph *cswHpol = new TGraph();
        
        
        calculateCSW(usefulAtriEvPtrOut, excludedChannels, cutoffTime, arrivalTimes, cswVpol, cswHpol, runNumber, sampleNs, dt, debugMode);  //Disabling this because of a seg fault - JCF 8/28/2024
        // calculateCSW(usefulAtriEvPtrOut, excludedChannels, cutoffTime, directPeakTimes, cswVpol, cswHpol, dt, debugMode);

        if (debugMode) {
            cout << "Sizes of csw TGraphs: "<< endl;
            cout << "cswVpol->GetN() = " << cswVpol->GetN() << endl;
            cout << "cswHpol->GetN() = " << cswHpol->GetN() << endl;

            cout << "Creating TCanvas for CSW" << endl;
            TCanvas *cCsw = new TCanvas("","", 1600, 1600);
            cCsw->Divide(1,2); 

            cCsw->cd(1); gPad->SetGrid(1,1);
            cswVpol->Draw();
            char vCswTitle[500];
            sprintf(vCswTitle,"A%d Run %s Event %d VPol CSW", settings1->DETECTOR_STATION_ARAROOT, runNumber, usefulAtriEvPtr->eventNumber);
            cswVpol->SetTitle(vCswTitle);

            cCsw->cd(2); gPad->SetGrid(1,1);
            cswHpol->Draw();
            char hCswTitle[500];
            sprintf(hCswTitle,"A%d Run %s Event %d HPol CSW", settings1->DETECTOR_STATION_ARAROOT, runNumber, usefulAtriEvPtr->eventNumber);
            cswHpol->SetTitle(hCswTitle);

            char cswTitle[500];
            // sprintf(cswTitle, "%s/waveformCSW_%s.png", outputDir, runNumber);
            sprintf(cswTitle, "%s/waveformCSW.png", outputDir);
            cCsw->Print(cswTitle);   
        }
        
        for (int ch=0; ch<16; ch++) {
            for (int i=0; i<cswVpol->GetN(); i++){
                int elecChan = geomTool->getElecChanFromRFChan(ch, settings1->DETECTOR_STATION_ARAROOT);
                //Vpol
                if (ch<8) {
                    usefulAtriCswPtrOut->fVolts[elecChan].push_back(cswVpol->GetY()[i]);
                    usefulAtriCswPtrOut->fTimes[elecChan].push_back(cswVpol->GetX()[i]);
                }
                //HPol
                else {
                    usefulAtriCswPtrOut->fVolts[elecChan].push_back(cswHpol->GetY()[i]);
                    usefulAtriCswPtrOut->fTimes[elecChan].push_back(cswHpol->GetX()[i]);
                }
            }
        }        
        
        if (not separateCsw) {
            usefulAtriEvPtrOut->fVolts.clear();
            usefulAtriEvPtrOut->fTimes.clear();             
            for (int ch=0; ch<16; ch++) {
                for (int i=0; i<cswVpol->GetN(); i++){
                    int elecChan = geomTool->getElecChanFromRFChan(ch, settings1->DETECTOR_STATION_ARAROOT);
                    //Vpol
                    if (ch<8) {
                        usefulAtriEvPtrOut->fVolts[elecChan].push_back(cswVpol->GetY()[i]);
                        usefulAtriEvPtrOut->fTimes[elecChan].push_back(cswVpol->GetX()[i]);
                    }
                    //HPol
                    else {
                        usefulAtriEvPtrOut->fVolts[elecChan].push_back(cswHpol->GetY()[i]);
                        usefulAtriEvPtrOut->fTimes[elecChan].push_back(cswHpol->GetX()[i]);
                    }
                }
            }       
        }
        for (int i = 0; i<16; i++) {
            cutoffTimeCsw[i] = 1e100;
        }
        double arrivalTimeMax = *max_element(arrivalTimes, arrivalTimes + 16);
        double arrivalTimeMin = *min_element(arrivalTimes, arrivalTimes + 16);
        double arrivalTimeAvg = 0;
        for (int i=0; i<16; i++) {
            arrivalTimeAvg += arrivalTimes[i]/16;
        } 
        //Set condition for empty events or those that resolve in the shadow zone.
        if (arrivalTimeAvg != -1000 and cswVpol->GetN() != 0 and cswHpol->GetN() != 0) {
            peakHilbert(usefulAtriCswPtrOut, hilbertPeakCswOut, peakTimeCswOut, cutoffTimeCsw, dt, toleranceNs=toleranceNs, findPolarity=findPolarity, debugMode=debugMode);
            cout << "$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$" << endl;
            for (int i=0; i<16; i++) {
                
                cout << hilbertPeakCswOut[i] << endl;
            }
            calculatePsi(hilbertPeakCswOut, psiRecoCsw, findPolarity);
            calculatePsiSolutions(psiRecoCsw, psiRecoCsw2, psiRecoCsw3, psiRecoCsw4);
        }
        
        for (int i=0; i<8; i++) {
            
            calculateNuTrajectory(psiRecoCsw[i]*PI/180, launch_theta[i]*PI/180, launch_phi[i]*PI/180, cherenkovAngle, theta_nutrajectCsw[i], phi_nutrajectCsw[i]);
            
            if (not dataLike) {
                 // psiTrue[i] = calculateTruePsi(eventPtr, report, launch_theta[i]*PI/180, launch_phi[i]*PI/180);
                psiTrue[i] = calculateTruePsi(eventPtr, report, true_launch_theta[i]*PI/180, true_launch_phi[i]*PI/180);
                Vector launch_vector = calculateLaunchVector(true_launch_theta[i]*PI/180, true_launch_phi[i]*PI/180);
                double viewingAngle = calculateViewingAngle(eventPtr->Nu_Interaction[0].nnu, launch_vector);
                true_viewingAngle[i] = viewingAngle;
                // calculateNuTrajectory(psiRecoCsw[i]*PI/180, launch_theta[i]*PI/180, launch_phi[i]*PI/180, viewingAngle, theta_nutrajectCsw[i], phi_nutrajectCsw[i]);
                // calculateNuTrajectory(psiTrue[i]*PI/180, launch_theta[i]*PI/180, launch_phi[i]*PI/180, cherenkovAngle, theta_nutrajectCsw[i], phi_nutrajectCsw[i]);
            // }        
            // calculatePsiSolutions(psiTrue, psiTrue2, psiTrue3, psiTrue4);
            // if (not dataLike) {
                // psiTrue[i] = calculateTruePsi(eventPtr, report, true_launch_theta[i]*PI/180, true_launch_phi[i]*PI/180);
                calculateNuTrajectory(psiTrue[i]*PI/180, true_launch_theta[i]*PI/180, true_launch_phi[i]*PI/180, viewingAngle, true_theta_nutraject[i], true_phi_nutraject[i]);
                calculateNuTrajectory(psiTrue[i]*PI/180, true_launch_theta[i]*PI/180, true_launch_phi[i]*PI/180, cherenkovAngle, true_theta_nutraject2[i], true_phi_nutraject2[i]);
                calculateNuTrajectory(psiRecoCsw[i]*PI/180, true_launch_theta[i]*PI/180, true_launch_phi[i]*PI/180, viewingAngle, true_theta_nutraject3[i], true_phi_nutraject3[i]);
                calculateNuTrajectory(psiTrue[i]*PI/180, launch_theta[i]*PI/180, launch_phi[i]*PI/180, viewingAngle, true_theta_nutraject4[i], true_phi_nutraject4[i]);
                // calculateNuTrajectory(psiTrue[i]*PI/180, true_launch_theta[i]*PI/180, true_launch_phi[i]*PI/180, nofz_atTrueVertex, theta_nutrajectCsw[i], phi_nutrajectCsw[i]);  //TODO: Remove this, as it's only for debugging.
                // cout << "theta_nutraject_true = " << eventPtr->Nu_Interaction[0].nnu.Theta()*180/PI << endl;
                // cout << "phi_nutraject_true = " << eventPtr->Nu_Interaction[0].nnu.Phi()*180/PI << endl;
            }              
        }
        
        //Average the trajectories calculated in the channels that passed the CSW
        int numChannels=0;
        mean_theta_nutraject = 0;
        mean_phi_nutraject = 0;
        cout << "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@" << endl;
        for (int i=0; i<8; i++) {
            cout << "i = " << i << endl;
            bool checkExcluded = std::find(cswExcludedChannelPairs.begin(), cswExcludedChannelPairs.end(), i) != cswExcludedChannelPairs.end();
            if (not checkExcluded) {
                cout << "Including i = " << i << endl;
                cout << "theta_nutrajectCsw[i] = " << theta_nutrajectCsw[i] << endl;
                cout << "phi_nutrajectCsw[i] = " << phi_nutrajectCsw[i] << endl;
                mean_theta_nutraject += theta_nutrajectCsw[i];
                mean_phi_nutraject += phi_nutrajectCsw[i];
                numChannels++;                
            }
            else {
                cout << "Excluding i = " << i << endl;
            }
            cout << "mean_theta_nutraject = " << mean_theta_nutraject << endl;
            cout << "mean_phi_nutraject = " << mean_phi_nutraject << endl;            
            
        }
        mean_theta_nutraject=mean_theta_nutraject/numChannels;
        mean_phi_nutraject=mean_phi_nutraject/numChannels;
        cout << "mean_theta_nutraject = " << mean_theta_nutraject << endl;
        cout << "mean_phi_nutraject = " << mean_phi_nutraject << endl;          
        
        
        fpOut->cd();
        outTree->Fill();
        outTreeCsw->Fill();
        // outTreeVertexReco->Fill();
        // if (not dataLike) {
        //     // outTreeMCTruth = AraTree2->CloneTree();
        //     outTreeMCTruth->Fill();
        // }

        if (debugMode) {
            TCanvas *c2 = new TCanvas("","", 1600, 1600);
            c2->Divide(4,4);          

//             //Create array of channel pairs to exclude in the coherent sum
//             std::vector<int> cswExcludedChannelPairs;
//             for(int i=0; i<8; i++){
//                 bool checkVpol = std::find(excludedChannels.begin(), excludedChannels.end(), i) != excludedChannels.end();
//                 bool checkHpol = std::find(excludedChannels.begin(), excludedChannels.end(), i+8) != excludedChannels.end();
//                 if (checkVpol or checkHpol) {
//                     cout << "Excluding channel pair " << i << endl;
//                     cswExcludedChannelPairs.push_back(i);
//                 }

//             }

            // double arrivalTimeMax = *max_element(arrivalTimes, arrivalTimes + 16);
            // double arrivalTimeMin = *min_element(arrivalTimes, arrivalTimes + 16);
            // double arrivalTimeAvg = 0;
            // for (int i=0; i<16; i++) {
            //     arrivalTimeAvg += arrivalTimes[i]/16;
            // }

            for(int i=0; i<16; i++){

                //Check if channel is in the excluded channel list
                bool checkExcluded = std::find(cswExcludedChannelPairs.begin(), cswExcludedChannelPairs.end(), i%8) != cswExcludedChannelPairs.end();

                TGraph *gr = usefulAtriEvPtrOut->getGraphFromRFChan(i);

                for(int k=0; k<gr->GetN(); k++){
                    gr->GetX()[k] = gr->GetX()[k] - arrivalTimes[i] + arrivalTimeAvg;
                    // gr->GetX()[k] = gr->GetX()[k] - directPeakTimes[i];
                }

                gr = FFTtools::getInterpolatedGraph(gr,dt);
                // gr->GetXaxis()->SetLimits(-30, 130);

                c2->cd(i+1); gPad->SetGrid(1,1);
                gr->Draw();
                if (checkExcluded) {
                    gr->SetLineColor(2);
                }
                char vTitle[500];
                sprintf(vTitle,"A%d Run %s Event %d Ch. %.2d", settings1->DETECTOR_STATION_ARAROOT, runNumber, usefulAtriEvPtr->eventNumber, i);
                gr->SetTitle(vTitle);            

            }


            char title[500];
            // sprintf(title, "%s/waveformDeconvolvedTimeshifted_%s.png", outputDir, runNumber);
            sprintf(title, "%s/waveformDeconvolvedTimeshifted.png", outputDir);
            c2->Print(title); 
        }
        
        if (debugMode) {
            //Testing peakHilbert function
            // peakHilbert(usefulAtriEvPtrOut, hilbertPeakOut, peakTimeOut, cutoffTime, dt, debugMode);
            cout << "hilbertPeakOut = " << endl;
            for (int k=0; k<16; k++) {
                cout << hilbertPeakOut[k] << ", ";
            }
            cout << endl;
            cout << "peakTimeOut = " << endl;
            for (int k=0; k<16; k++) {
                cout << peakTimeOut[k] << ", ";
            }
            cout << endl;
            
            calculatePsi(hilbertPeakOut, psiReco, findPolarity);
            
            cout << "psi = " << endl;
            for (int k=0; k<8; k++) {
                cout << psiReco[k] << ", ";
            }
            cout << endl; 
            cout << "psiRecoCsw = " << endl;
            for (int k=0; k<8; k++) {
                cout << psiRecoCsw[k] << ", ";
            }
            cout << endl;              
            cout << "psiTrue = " << endl;
            for (int k=0; k<8; k++) {
                cout << psiTrue[k] << ", ";
            }
            cout << endl;               
            
            
        }
        
        usefulAtriCswPtrOut->fVolts.clear();
        usefulAtriCswPtrOut->fTimes.clear();        
        usefulAtriEvPtrOut->fVolts.clear();
        usefulAtriEvPtrOut->fTimes.clear();   
        if (not dataLike and debugMode) {        
            for (int i=0; i<16; i++) {
                cout << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << endl;
                cout << "i = " << i << endl;
                cout << "arrivalTimes[i] = " << arrivalTimes[i] << endl;
                cout << "psiReco[i%8] = " << psiReco[i%8] << endl;
                cout << "psiRecoCsw[i%8] = " << psiRecoCsw[i%8] << endl;
                // cout << "psiRecoCsw2[i%8] = " << psiRecoCsw2[i%8] << endl;
                // cout << "psiRecoCsw3[i%8] = " << psiRecoCsw3[i%8] << endl;
                // cout << "psiRecoCsw4[i%8] = " << psiRecoCsw4[i%8] << endl;
                cout << "psiTrue[i%8] = " << psiTrue[i%8] << endl;
                // cout << "psiTrue2[i%8] = " << psiTrue2[i%8] << endl;
                // cout << "psiTrue3[i%8] = " << psiTrue3[i%8] << endl;
                // cout << "psiTrue4[i%8] = " << psiTrue4[i%8] << endl;
                cout << "launch_theta[i] = " << launch_theta[i] << endl;
                cout << "true_launch_theta[i] = " << true_launch_theta[i] << endl;
                cout << "launch_phi[i] = " << launch_phi[i] << endl;  
                cout << "true_launch_phi[i] = " << true_launch_phi[i] << endl;
                cout << "theta_nutraject[i%8] = " << theta_nutraject[i%8] << endl;
                cout << "theta_nutrajectCsw[i%8] = " << theta_nutrajectCsw[i%8] << endl;
                cout << "true_theta_nutraject[i%8] = " << true_theta_nutraject[i%8] << endl;
                // cout << "true_theta_nutraject2[i%8] = " << true_theta_nutraject2[i%8] << endl;
                // cout << "true_theta_nutraject3[i%8] = " << true_theta_nutraject3[i%8] << endl;
                // cout << "true_theta_nutraject4[i%8] = " << true_theta_nutraject4[i%8] << endl;
                cout << "eventPtr->Nu_Interaction[0].nnu.Theta()*180/PI = " << eventPtr->Nu_Interaction[0].nnu.Theta()*180/PI << endl;
                cout << "phi_nutraject[i%8] = " << phi_nutraject[i%8] << endl;  
                cout << "phi_nutrajectCsw[i%8] = " << phi_nutrajectCsw[i%8] << endl; 
                cout << "true_phi_nutraject[i%8] = " << true_phi_nutraject[i%8] << endl;
                // cout << "true_phi_nutraject2[i%8] = " << true_phi_nutraject2[i%8] << endl;
                // cout << "true_phi_nutraject3[i%8] = " << true_phi_nutraject3[i%8] << endl;
                // cout << "true_phi_nutraject4[i%8] = " << true_phi_nutraject4[i%8] << endl;
                cout << "eventPtr->Nu_Interaction[0].nnu.Phi()*180/PI = " << eventPtr->Nu_Interaction[0].nnu.Phi()*180/PI << endl;            
                testTrajectoryCalculation(eventPtr, report, true_launch_theta[i]*PI/180, true_launch_phi[i]*PI/180, eventPtr->Nu_Interaction[0].indexN);

                }
//             cout << "mean_theta_nutraject = " << mean_theta_nutraject << endl;
//             cout << "mean_phi_nutraject = " << mean_phi_nutraject << endl;
            
//             cout << "event->Nu_Interaction[0].indexN = " << eventPtr->Nu_Interaction[0].indexN << endl;
//             cout << "event->Nu_Interaction[0].changle*180/PI = " << eventPtr->Nu_Interaction[0].changle*180/PI << endl;    
//             cout << "nofz_atVertex = " << nofz_atVertex << endl;
//             cout << "acos(1/nofz_atVertex)*180/PI = " << acos(1/nofz_atVertex)*180/PI << endl; 

        }
    } //event loop

    fpOut->Write();
    fpOut->Close();
    fp->Close();
    fp2->Close();
    delete fpOut;
    delete fp;
    delete fp2;
    
    cout << "**************************************************" << endl;
    cout << "Output written to:\t " <<outfile_name<<endl;    
    cout << "**************************************************" << endl;
    
} //main    