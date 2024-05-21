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
#include "tools.h"

//TODO: Have outgoing pointer equal the incoming pointer by using a filler function to copy the original information, then I replace the voltage info with my deconvolved voltage.
UsefulAtriStationEvent *usefulAtriEvPtr;
UsefulAtriStationEvent *usefulAtriEvPtrOut;
UsefulAtriStationEvent *usefulAtriCswPtrOut;

bool debugMode = false;
bool separateCsw = true;
bool useMCTruth = true;
bool forceSinglePeak = true;
bool weinerCorrection = true;
double sampleNs=80;  //Sample window for noise and signal in the PSD SNR calculation.
// double dt=0.1; //Sample rate for interpolation

int main(int argc, char **argv)
{
    // if(argc<7) {
    if(argc<8) {
        // std::cout << "Usage\n" << argv[0] << " <station> <config> <runnum> <input root file> <input reco file> <output_dir> <setup_file>\n";
        // std::cout << "e.g.\n" << argv[0] << " 2 6 AraOut.root recangle_out_run<runnum>.root output setup.txt /\n";
        std::cout << "Usage\n" << argv[0] << " <runnum> <band-pass minimum frequency (MHz)> <band-pass maximum frequency (MHz)> <setup file> <input root file> <input reco file> <output directory>\n";
        std::cout << "e.g.\n" << argv[0] << " 1000 150 300 setup.txt AraOut.root recangle_out_run1000.root output/\n";        
        return 0;
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
    bool dataLike = false;
    bool calibrated;

    
    //Trying condition where it checks for usefulAtriStation branch and imports according to that.
    if(!simSettingsTree) { 
        dataLike = true;            
        std::cerr << "Can't find AraTree.  Importing as real data.\n";
        // TTree* atriExists=(TTree*) fp->Get("UsefulAtriStationEvent");
        // if (!atriExists) {
        //Checks if usefulAtri branch exists.  If not, fata gets imported as uncalibrated.
        if (not eventTree->GetBranch("UsefulAtriStationEvent")) {
            calibrated = false;
            cout << "Can't find UsefulAtriStationEvent Tree.  Importing as uncalibrated." << endl;
            eventTree->SetBranchAddress("event",&rawAtriEvPtr);
        }
        else {
            calibrated = true;
            eventTree->SetBranchAddress("UsefulAtriStationEvent", &usefulAtriEvPtr);
        }
        double weight = 1;
    }
    // sim like
    else {
        dataLike = false;
        calibrated = true;
        std::cerr << "AraTree exists.  Importing as simulated data.\n";
        eventTree->SetBranchAddress("UsefulAtriStationEvent", &usefulAtriEvPtr);
        double weight;
        eventTree->SetBranchAddress("weight", &weight);   
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
    double arrivalTimes[16];
    double directPeakTimes[16];
    double refPeakTimes[16];    

    double v_snr;
    double h_snr;
    
    // Testing using the true rf angles
    if (useMCTruth){
        vertexReco->SetBranchAddress("true_arrivalThetas", reco_arrivalThetas);
        vertexReco->SetBranchAddress("true_arrivalPhis", reco_arrivalPhis);  
    }
    // end testing
    else {
        vertexReco->SetBranchAddress("reco_arrivalThetas", reco_arrivalThetas);
        vertexReco->SetBranchAddress("reco_arrivalPhis", reco_arrivalPhis);
    }
    vertexReco->SetBranchAddress("cutoffTime", cutoffTime);
    vertexReco->SetBranchAddress("arrivalTimes", arrivalTimes);
    // vertexReco->SetBranchAddress("excludedChannels", &excludedChannels);    
    vertexReco->SetBranchAddress("v_snr", &v_snr);
    vertexReco->SetBranchAddress("h_snr", &h_snr);   
    vertexReco->SetBranchAddress("directPeakTimes", &directPeakTimes);  
    vertexReco->SetBranchAddress("refPeakTimes", &refPeakTimes);      
    Long64_t numEntriesVertex=vertexReco->GetEntries();
    cout << "Vertex Reco tree opened! Has " << numEntriesVertex << " entries!" << endl;;
    
    printf("------------------\n");
    printf("Input files loaded.  Setting up detector stuff.\n");
    printf("------------------\n");
    
    // string setupfile;
    // setupfile = argv[2];
    Settings *settings1 = new Settings();
    settings1->ReadFile(setupfile); 
    IceModel *icemodel=new IceModel(settings1->ICE_MODEL + settings1->NOFZ*10,settings1->CONSTANTICETHICKNESS * 1000 + settings1->CONSTANTCRUST * 100 + settings1->FIXEDELEVATION * 10 + 0,settings1->MOOREBAY);// creates Antarctica ice model
    Detector *detector = new Detector(settings1, icemodel, setupfile);  
    Report *report = new Report(detector, settings1);
    
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
    
    TTree *outTree = new TTree("eventTree", "eventTree");
    TTree *outTreeCsw = new TTree("coherentSum", "coherentSum");
    outTree->Branch("UsefulAtriStationEvent", &usefulAtriEvPtrOut);
    outTreeCsw->Branch("UsefulAtriStationEvent", &usefulAtriCswPtrOut);
    
    //Need to grab lengths of voltage and time arrays from eventTree to initialize the branches in the outfile.
    Int_t fNumChannels; ///< The number of channels
    std::map< Int_t, std::vector <Double_t> > fTimesOut; ///< The times of samples
    std::map< Int_t, std::vector <Double_t> > fVoltsOut; ///< The voltages of samples   
    
    //Initialize voltage and time arrays for the coherent sum.
    Int_t fNumChannelsCsw = 2; ///< The number of channels (one for Vpol, one for Hpol)
    std::map< Int_t, std::vector <Double_t> > fTimesCswOut; ///< The times of samples
    std::map< Int_t, std::vector <Double_t> > fVoltsCswOut; ///< The voltages of samples      
     
    
    //Loop over events
    // int loopedEntries = numEntries;
    int loopStart;
    int loopEnd;
    if (debugMode) {
        loopStart=0;
        loopEnd=1;
    }
    else {
        loopStart=0;
        loopEnd=numEntries;
    }
    
    // for(Long64_t event=0;event<loopedEntries;event++) {
    for(Long64_t event=loopStart;event<loopEnd;event++) {
    // for(Long64_t event=679;event<680; event++) {  //Debugging and running over events enar desired event to save time in loop.      
    // for(Long64_t event=53530;event<53539;event++) {  //Debugging and running over events near desired event to save time in loop.    
        std::cout<<"Looking at event number "<<event<<std::endl;
        fp->cd();

        eventTree->GetEntry(event);
        
        // if (rawAtriEvPtr->eventNumber != 5572){  //Test for A4 Run 6136
        //     continue;
        // }  

        // if (rawAtriEvPtr->eventNumber != 53535){  //Test for A4 Run 6128
        //     continue;
        // }  

        // if (rawAtriEvPtr->eventNumber != 679){  //Test for A2 Run 12559
        //     continue;
        // }
        cout << "Importing vertex reco info" << endl;
        vertexReco->GetEntry(event);
        // vertexReco->GetEntry(0);  //For debugging purposes.  TODO: Change back!
        
        if (not calibrated) {
            cout << "Triggering datalike condition." << endl;
            delete usefulAtriEvPtr;         
            usefulAtriEvPtr = new UsefulAtriStationEvent(rawAtriEvPtr, AraCalType::kLatestCalib);
            
            
        }
        
        //Need this for importing the vertexReco from simulated data, as simulated data gets stored by string, then antenna.
        int vertexRecoElectToRFChan[] = {14,2,6,10,  //VPols
                                         12,0,4,8,   //VPols
                                         15,3,7,11,  //HPols
                                         13,1,5,9};  //HPols

        for(int i=0; i<16; i++){
            TGraph *gr = usefulAtriEvPtr->getGraphFromRFChan(i);
            //Save initial and final time for truncating the padded arrays before output.
            double timeStart = gr->GetX()[0];
            double timeEnd = gr->GetX()[gr->GetN()-1];
            //Apply time shift to center waveform
            double timeShift = (timeStart+timeEnd)/2;
            for(int k=0; k<gr->GetN(); k++){
                gr->GetX()[k] = gr->GetX()[k] - timeShift;
            }               
            
            //Interpolate graph to 0.5 ns resolution (or whatever TIMESTEP is specified in the setup file)
            gr = FFTtools::getInterpolatedGraph(gr,dt);
            //Pad waveform to a factor of two as specified in the setup file. - JCF 9/27/2023
            if (gr->GetN() < settings1->NFOUR/2) {
                gr = FFTtools::padWaveToLength(gr, settings1->NFOUR/2);
            }
            // Get length of padded waveform
            int waveform_bin = gr->GetN();
            
            
            double heff_lastbin;
            double freq_lastbin;
            double time[waveform_bin];
            double voltage[waveform_bin];
            double volts_forint[settings1->NFOUR / 2];
            double T_forint[settings1->NFOUR / 2];
            
            //TODO: This init_T isn't dynamic to the imported data.  Should make have it defined based on the input waveform.
            double init_T = settings1->TIMESTEP *-1.e9 *((double) settings1->NFOUR / 4);
            for(int k=0; k<waveform_bin; k++){
                time[k] = gr->GetX()[k];
                voltage[k] = gr->GetY()[k];
            }            
            // delete gr;
            // cout << "time[0] = " << time[0] << endl;
            
            for (int m = 0; m < settings1->NFOUR / 2; m++)
            {
                T_forint[m] = init_T + m*dt;   // in ns
            }
            
            //Importing the cutoff time between spicecore peaks
            double cutoffTimeChannel;
            // if (forceSinglePeak) {
            //     cutoffTime[i] = time[waveform_bin-1]; //TODO: Adding this for debugging.  Be sure to remove.
            // }
               
            if (!cutoffTime or forceSinglePeak) {
                cout << "Using single-peak condition" << endl;
                // cutoffTimeChannel = time[waveform_bin-1];
                // cout << "cutoffTimeChannel = " << cutoffTimeChannel << endl;
                // cout << "waveform_bin = " << waveform_bin << endl;
                cutoffTime[i] = time[waveform_bin-1];
                cout << "cutoffTime[i] = " << cutoffTime[i] << endl;
                cout << "waveform_bin = " << waveform_bin << endl;                
            } 
            // else {
            //     cutoffTimeChannel = cutoffTime[i];
            // }            
    
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
        
            //Testing using the actual polarization vector
            // double psi = argv[3]*PI/180;
            // double newPol_vectorX = -cos(psi)*cos(antenna_theta*PI/180)*cos(antenna_phi*PI/180) + sin(psi)*sin(antenna_phi*PI/180);
            // double newPol_vectorY = -cos(psi)*cos(antenna_theta*PI/180)*sin(antenna_phi*PI/180) - sin(psi)*cos(antenna_phi*PI/180);
            // double newPol_vectorZ = cos(psi)*sin(antenna_theta*PI/180);            

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
            
            
        
            int Ntmp = settings1->TIMESTEP *1.e9 / dT_forfft;
            
            int Nnew = 1;         
            while (Ntmp > 1)
            {
                Ntmp = Ntmp / 2;
                Nnew = Nnew *2;             
            }
            Nnew = Nnew * settings1->NFOUR / 2;  
            
            // cout << "Nnew = " << Nnew << endl;
            
            //Initialize and populate the voltage and time arrays for the fft.
            double V_forfft[Nnew];
            double T_forfft[Nnew];
            
            createFFTarrays(voltage, time, settings1, V_forfft, T_forfft, waveform_bin, Nnew);
            
            //Testing get PSD SNR function
            
            double snrPsd[int(Nnew/2)];
            double freqPsd[int(Nnew/2)];
            // cout << "directPeakTimes[" << i << "] = " << directPeakTimes[i] << endl; 
            getPowerSpectrumSNR(gr, directPeakTimes[i], waveform_bin, Nnew, settings1, snrPsd, freqPsd, dt, sampleNs, debugMode);
            
//             if (debugMode) {
//                 cout << "snrPsd = " << endl;
//                 for (int k=0; k<Nnew/2; k++) {
//                     cout << snrPsd[k] << ", ";
//                 }
//                 cout << endl;   
                
//                 cout << "freqPsd = " << endl;
//                 for (int i=0; i<Nnew/2; i++) {
//                     cout << freqPsd[i] << ", ";
//                 }
//                 cout << endl;                
//             }
            

            
            delete gr;
            
            // get spectrum with zero padded WF
            Tools::realft(V_forfft, 1, Nnew); 
            
            // if (debugMode) {
            //     cout << "*****************************" << endl;
            //     cout << "V_forfft (before factors) = " << endl;
            //     for (int n=0; n<Nnew/2; n++) {
            //         cout << V_forfft[2*n] << " + i " << V_forfft[2*n+1] << ", ";
            //     }
            //     cout << endl;         
            // }
            
            dF_Nnew = 1. / ((double)(Nnew) *(dT_forfft) *1.e-9);    // in Hz
            
            // cout << "dF_Nnew = " << dF_Nnew << endl;

            freq_tmp = dF_Nnew *((double) Nnew / 2. + 0.5); // in Hz 0.5 to place the middle of the bin and avoid zero freq
            
            // if (i == 9) {
            //     cout << "*****************************" << endl;
            //     cout << "freq = " << endl;
            //     for (int k=0; k<Nnew/2; k++) {
            //         cout << dF_Nnew *((double) k + 0.5)*1e-6 << ", ";
            //     }
            //     cout << endl;
            // }            

            freq_lastbin = freq_tmp;
            
            heff_lastbin = report->GaintoHeight(detector->GetGain_1D_OutZero(freq_tmp *1.E-6,   // to MHz
                                                antenna_theta, antenna_phi, pol_ant),
                                                freq_tmp, nice);                  
         
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
                // invert entire elect chain gain, phase
                if (n > 0)
                {  
                // if (debugMode) {
                //     cout << "Before Elec " << endl;                    
                //     cout <<  "V_forfft[2*n] = " <<  V_forfft[2*n] << endl;
                //     cout <<  "V_forfft[2*n+1] = " <<  V_forfft[2*n+1] << endl;      
                // }                           
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
                // if (debugMode) {
                //     cout << "After Elec, before Ant "<< endl;                    
                //     cout <<  "V_forfft[2*n] = " <<  V_forfft[2*n] << endl;
                //     cout <<  "V_forfft[2*n+1] = " <<  V_forfft[2*n+1] << endl;      
                // }                                
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
                
                // if (debugMode) {
                //     cout << "After ant " << endl;
                //     cout << "realWeinerCorr = " << realWeinerCorr << endl;
                //     cout << "imagWeinerCorr = " << imagWeinerCorr << endl;                    
                //     cout <<  "V_forfft[2*n] = " <<  V_forfft[2*n] << endl;
                //     cout <<  "V_forfft[2*n+1] = " <<  V_forfft[2*n+1] << endl;      
                // }                
                
                if (weinerCorrection) {
                    //Correction for weiner deconvolution
                    // cout << "******************************************" << endl;
                    // cout << "freq_tmp = " << freq_tmp*1e-6 << endl;
                    wienerDeconvolution(V_forfft[2 *n], V_forfft[2 *n + 1], realWeinerCorr, imagWeinerCorr, snrPsd[n]);
                }
                
                else {
                    if (freq_tmp > 850.*1.e6 or freq_tmp < 100.*1.e6) {
                    // if (freq_tmp > 600.*1.e6 or freq_tmp < 200.*1.e6) {                
                        V_forfft[2*n] = 0;
                        V_forfft[2*n+1] = 0;
                    }  
                }
                

  
                
                double weight = 1;  // Setting initial weight to one, then applying bandpass.  Weight is then multiplied by signal in this bin.
                int order = 8;
                weight /= sqrt(1 + TMath::Power(freqMin/freq_tmp, 4*order));
                weight /= sqrt(1 + TMath::Power(freq_tmp/freqMax, 4*order));
                V_forfft[2*n] *= weight;
                V_forfft[2*n+1] *= weight; 
                //End Butterworth filter 
            }   // end for freq bin
                            
            // now get time domain waveform back by inv fft  
            // if (debugMode) {
            //     cout << "*****************************" << endl;
            //     cout << "V_forfft (before InvFFT) = " << endl;
            //     for (int n=0; n<Nnew/2; n++) {
            //         cout << V_forfft[2*n] << " + i " << V_forfft[2*n+1] << ", ";
            //     }
            //     cout << endl;         
            // }
            
            Tools::realft(V_forfft, -1, Nnew);
            
            // if (debugMode) {
            //     cout << "*****************************" << endl;
            //     cout << "V_forfft (after InvFFT) = " << endl;
            //     for (int n=0; n<Nnew/2; n++) {
            //         cout << V_forfft[2*n] << " + i " << V_forfft[2*n+1] << ", ";
            //     }
            //     cout << endl;         
            // }
            
            //TODO: Make this 160 more data-driven.  Shouldn't be hard-coded, but constants that can be changed by the user.
            // if (antenna_theta > 160 or antenna_theta < 90) {
            //     cout << "Event outside of theta range.  Setting voltage to zero and moving to next event." << endl;
            //     for (int i = 0; i < sizeof(V_forfft) / sizeof(V_forfft[0]); i++) {
            //       V_forfft[i] = 1;
            //     }                
            //     // continue;
            // }                 
            Tools::SincInterpolation(Nnew, T_forfft, V_forfft, settings1->NFOUR / 2, T_forint, volts_forint);   
                           
            
            //Now write deconvolved voltage data to file.
            for (int n = 0; n < settings1->NFOUR / 2; n++)
            // for (int n = 0; n < waveform_bin; n++)
            {
                int elecChan = AraGeomTool::Instance()->getElecChanFromRFChan(i, settings1->DETECTOR_STATION);
                // not pure noise mode (we need signal)
                usefulAtriEvPtrOut->fVolts[elecChan].push_back(settings1->ARBITRARY_EVENT_ATTENUATION *volts_forint[n] *2. / (double)(Nnew));  // 2/N for inverse FFT normalization factor
                usefulAtriEvPtrOut->fTimes[elecChan].push_back(T_forint[n]+timeShift);    
                // cout << "T_forint[" << n << "] = " << T_forint[n] << endl;;
            }
            usefulAtriEvPtrOut->stationId = settings1->DETECTOR_STATION;
            usefulAtriCswPtrOut->stationId = settings1->DETECTOR_STATION;
            
        } //channel loop
        usefulAtriEvPtrOut->eventNumber = usefulAtriEvPtr->eventNumber;
        usefulAtriEvPtrOut->unixTime = usefulAtriEvPtr->unixTime;

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
        //TCanvas for the coherently summed waveform separated by VPol and HPol
        if (debugMode) {
            TCanvas *cTimeshift = new TCanvas("","", 1600, 1600);
            cTimeshift->Divide(4,4);          

            //Create array of channel pairs to exclude in the coherent sum
            std::vector<int> cswExcludedChannelPairs;
            for(int i=0; i<8; i++){
                bool checkVpol = std::find(excludedChannels.begin(), excludedChannels.end(), i) != excludedChannels.end();
                bool checkHpol = std::find(excludedChannels.begin(), excludedChannels.end(), i+8) != excludedChannels.end();
                if (checkVpol or checkHpol) {
                    cout << "Excluding channel pair " << i << endl;
                    cswExcludedChannelPairs.push_back(i);
                }

            }
            
            // cout << "bbb" << endl;
            
            double arrivalTimeMax = *max_element(arrivalTimes, arrivalTimes + 16);
            double arrivalTimeMin = *min_element(arrivalTimes, arrivalTimes + 16);        

            for(int i=0; i<16; i++){

                //Check if channel is in the excluded channel list
                bool checkExcluded = std::find(cswExcludedChannelPairs.begin(), cswExcludedChannelPairs.end(), i%8) != cswExcludedChannelPairs.end();

                TGraph *gr = usefulAtriEvPtrOut->getGraphFromRFChan(i);
                // gr = FFTtools::cropWave(gr, gr->GetX()[0], cutoffTime[i]);

                for(int k=0; k<gr->GetN(); k++){
                    gr->GetX()[k] = gr->GetX()[k];// - arrivalTimes[i] + arrivalTimeMax;
                }

                // cout << "Ch " << i << " time[0] = " << gr->GetX()[0] << endl;

                gr = FFTtools::getInterpolatedGraph(gr,dt);
                // gr = FFTtools::getInterpolatedGraphFreqDom(gr,0.1);

                cTimeshift->cd(i+1); gPad->SetGrid(1,1);
                gr->Draw();
                if (checkExcluded) {
                    gr->SetLineColor(2);
                }
                char vTitle[500];
                sprintf(vTitle,"Ch. %.2d", i);
                gr->SetTitle(vTitle);            

            }


            char title[500];
            sprintf(title, "waveformDeconvolved_%s.png", runNumber);
            cTimeshift->Print(title);  
        }
        
        //Plotting CSW deconvolved waveform
        //Create coherent sum logic
        TGraph *cswVpol = new TGraph();
        TGraph *cswHpol = new TGraph();
        
        calculateCSW(usefulAtriEvPtrOut, excludedChannels, cutoffTime, arrivalTimes, cswVpol, cswHpol, dt, debugMode);

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
            sprintf(vCswTitle,"VPol");
            cswVpol->SetTitle(vCswTitle);

            cCsw->cd(2); gPad->SetGrid(1,1);
            cswHpol->Draw();
            char hCswTitle[500];
            sprintf(hCswTitle,"HPol");
            cswHpol->SetTitle(hCswTitle);

            char cswTitle[500];
            sprintf(cswTitle, "waveformCSW_%s.png", runNumber);
            cCsw->Print(cswTitle);   
        }
        
        for (int ch=0; ch<16; ch++) {
            for (int i=0; i<cswVpol->GetN(); i++){
                int elecChan = AraGeomTool::Instance()->getElecChanFromRFChan(ch, settings1->DETECTOR_STATION);
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
                    int elecChan = AraGeomTool::Instance()->getElecChanFromRFChan(ch, settings1->DETECTOR_STATION);
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
        
        fpOut->cd();
        outTree->Fill();
        outTreeCsw->Fill();

        if (debugMode) {
            TCanvas *c2 = new TCanvas("","", 1600, 1600);
            c2->Divide(4,4);          

            //Create array of channel pairs to exclude in the coherent sum
            std::vector<int> cswExcludedChannelPairs;
            for(int i=0; i<8; i++){
                bool checkVpol = std::find(excludedChannels.begin(), excludedChannels.end(), i) != excludedChannels.end();
                bool checkHpol = std::find(excludedChannels.begin(), excludedChannels.end(), i+8) != excludedChannels.end();
                if (checkVpol or checkHpol) {
                    cout << "Excluding channel pair " << i << endl;
                    cswExcludedChannelPairs.push_back(i);
                }

            }

            double arrivalTimeMax = *max_element(arrivalTimes, arrivalTimes + 16);
            double arrivalTimeMin = *min_element(arrivalTimes, arrivalTimes + 16);
            double arrivalTimeAvg = 0;
            for (int i=0; i<16; i++) {
                arrivalTimeAvg += arrivalTimes[i]/16;
            }

            for(int i=0; i<16; i++){

                //Check if channel is in the excluded channel list
                bool checkExcluded = std::find(cswExcludedChannelPairs.begin(), cswExcludedChannelPairs.end(), i%8) != cswExcludedChannelPairs.end();

                TGraph *gr = usefulAtriEvPtrOut->getGraphFromRFChan(i);

                for(int k=0; k<gr->GetN(); k++){
                    gr->GetX()[k] = gr->GetX()[k] - arrivalTimes[i] + arrivalTimeAvg;
                }

                gr = FFTtools::getInterpolatedGraph(gr,dt);
                // gr->GetXaxis()->SetLimits(16,22);

                c2->cd(i+1); gPad->SetGrid(1,1);
                gr->Draw();
                if (checkExcluded) {
                    gr->SetLineColor(2);
                }
                char vTitle[500];
                sprintf(vTitle,"Ch. %.2d", i);
                gr->SetTitle(vTitle);            

            }


            char title[500];
            sprintf(title, "waveformDeconvolvedTimeshifted_%s.png", runNumber);
            c2->Print(title); 
        }
        
        usefulAtriCswPtrOut->fVolts.clear();
        usefulAtriCswPtrOut->fTimes.clear();        
        usefulAtriEvPtrOut->fVolts.clear();
        usefulAtriEvPtrOut->fTimes.clear();        
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