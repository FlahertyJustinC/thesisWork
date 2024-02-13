#include <iostream>

// ROOT Includes
#include "TFile.h"
#include "TTree.h"
#include "TGraph.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TLine.h"

// ARA Includes
#include "/users/PAS0654/jflaherty13/source/AraRoot/AraRoot_build/include/AraGeomTool.h"
#include "/users/PAS0654/jflaherty13/source/AraRoot/AraRoot_build/include/RayTraceCorrelator.h"
#include "/users/PAS0654/jflaherty13/source/AraRoot/AraRoot_build/include/UsefulAtriStationEvent.h"
#include "/cvmfs/ara.opensciencegrid.org/trunk/centos7/source/libRootFftwWrapper/include/FFTtools.h"


// AraSim includes
// #include "/cvmfs/ara.opensciencegrid.org/trunk/centos7/source/AraSim/Detector.h"
// #include "/cvmfs/ara.opensciencegrid.org/trunk/centos7/source/AraSim/Event.h"
// #include "/cvmfs/ara.opensciencegrid.org/trunk/centos7/source/AraSim/Position.h"
// #include "/cvmfs/ara.opensciencegrid.org/trunk/centos7/source/AraSim/RaySolver.h"
// #include "/cvmfs/ara.opensciencegrid.org/trunk/centos7/source/AraSim/IceModel.h"

#include "/users/PAS0654/jflaherty13/source/AraSim/Detector.h"
#include "/users/PAS0654/jflaherty13/source/AraSim/Event.h"
#include "/users/PAS0654/jflaherty13/source/AraSim/Position.h"
#include "/users/PAS0654/jflaherty13/source/AraSim/RaySolver.h"
#include "/users/PAS0654/jflaherty13/source/AraSim/IceModel.h"
#include "/users/PAS0654/jflaherty13/source/AraSim/Report.h"
#include "/users/PAS0654/jflaherty13/source/AraSim/Settings.h"

using namespace std;

#ifdef ARA_UTIL_EXISTS
    #include "UsefulIcrrStationEvent.h"
    ClassImp(UsefulIcrrStationEvent);
    #include "UsefulAtriStationEvent.h"
    ClassImp(UsefulAtriStationEvent);
#endif

UsefulAtriStationEvent *usefulAtriEvPtr;

int main(int argc, char **argv)
{
    if(argc<6) {
        std::cout << "Usage\n" << argv[0] << " <station> <config> <runnum> <input file> <output_dir> \n";
        std::cout << "e.g.\n" << argv[0] << " 2 6 AraOut.root output/\n";
        return 0;
    }
    
    double interpV = 0.4;
    double interpH = 0.625;
    
    printf("Opening file...\n");
    TFile *fp = TFile::Open(argv[4]);
    if(!fp) { std::cerr << "Can't open file\n"; return -1; }
    printf("File opened!\n");

    // data like
    TTree *eventTree = (TTree*) fp->Get("eventTree");
    printf("Event tree opened!\n");
    if(!eventTree) { std::cerr << "Can't find eventTree\n"; return -1; }
    eventTree->SetBranchAddress("UsefulAtriStationEvent", &usefulAtriEvPtr);
    double weight;
    eventTree->SetBranchAddress("weight", &weight);
    Long64_t numEntries=eventTree->GetEntries();
    
    printf("------------------\n");
    printf("Input files loaded.  Setting up detector stuff.\n");
    printf("------------------\n");
    
    string setupfile;
    setupfile = "SETUP/setup_variablePsi.txt";
    IceModel *iceModel = new IceModel(0 + 1*10, 0, 0);
    Settings *settings = new Settings();
    Detector *detector = new Detector(settings, iceModel, setupfile);
    Report *report = new Report(detector, settings);
    // configure settings    
    settings->Z_THIS_TOLERANCE = 1; // higher tolerance
    settings->Z_TOLERANCE = 0.05;
    settings->NOFZ=1; // make sure n(z) is turned on
    settings->RAY_TRACE_ICE_MODEL_PARAMS = 0; // set the ice model as user requested
    
    for(Long64_t event=0;event<numEntries;event++) {
        fp->cd();
        eventTree->GetEntry(event);
    
        
        
        if (event != 0){
               continue;
        }
        
        std::cout<<"Looking at event number "<<event<<std::endl;
        
        std::map<int, TGraph*> interpolatedWaveforms;
        for(int i=0; i<16; i++){
            TGraph *gr = usefulAtriEvPtr->getGraphFromRFChan(i);
            // TGraph *grInt = FFTtools::getInterpolatedGraph(gr,i<8?interpV:interpH);
            // interpolatedWaveforms[i] = grInt;
            int lenGraph = gr->GetN();
            delete gr;
        
    
            double freq_tmp, heff, antenna_theta, antenna_phi;  // values needed for apply antenna gain factor and prepare fft, trigger

            double dF_Nnew;

            int NFOUR = 1024;
            double nice = 1.79;

            double V_forfft[NFOUR];
            double T_forfft[NFOUR];
            int pol_ant;
            int gain_ch_no;
            Position Pol_vector;
            double Pol_factor;

            double heff_lastbin;
            double freq_lastbin;
            double time[lenGraph];
            double voltage[lenGraph];

            for(int k=0; k<lenGraph; k++){
                time[k] = gr->GetX()[k];
                voltage[k] = gr->GetY()[k];
            }
        


        //Stealing antenna and electronic response steps from AraSim, but applying the inverse functions instead.

            for (int n = 0; n < NFOUR / 2; n++)
            {

                freq_tmp = dF_Nnew *((double) n + 0.5); // in Hz 0.5 to place the middle of the bin and avoid zero freq

                heff = report->GaintoHeight(detector->GetGain_1D_OutZero(freq_tmp *1.E-6,   // to MHz
                        antenna_theta, antenna_phi, pol_ant),
                    freq_tmp, nice);

                // stations[i].strings[j].antennas[k].Heff[ray_sol_cnt].push_back(heff);

                if (n > 0)
                {

                    report->InvertAntFactors_Tdomain(detector->GetAntPhase_1D(freq_tmp *1.e-6, antenna_theta, antenna_phi, pol_ant),
                       heff, Pol_vector, pol_ant, Pol_factor, V_forfft[2 *n], V_forfft[2 *n + 1], antenna_theta, antenna_phi);
                }
                else
                {
                    report->InvertAntFactors_Tdomain_FirstTwo(heff, heff_lastbin, Pol_vector, pol_ant, Pol_factor, V_forfft[2 *n], V_forfft[2 *n + 1], antenna_theta, antenna_phi);

                }

                //
                // apply entire elect chain gain, phase
                //
                if (n > 0)
                {                                             
                    report->InvertElect_Tdomain(freq_tmp *1.e-6, detector, V_forfft[2 *n], V_forfft[2 *n + 1], gain_ch_no);                                                
                }
                else
                {
                    report->InvertElect_Tdomain_FirstTwo(freq_tmp *1.e-6, freq_lastbin *1.e-6, detector, V_forfft[2 *n], V_forfft[2 *n + 1], gain_ch_no);
                }
            }   // end for freq bin
        }
    }
    // now get time domain waveform back by inv fft
    // Tools::realft(V_forfft, -1, stations[i].strings[j].antennas[k].Nnew[ray_sol_cnt]);    
    
    
    
    
    
    
    
    
       
    
//     printf("------------------\n");
//     printf("Make Output Files\n");
//     printf("------------------\n");

//     char outfile_name[400];
//     sprintf(outfile_name, "%s/deconvolvedWaveforms_run_%d.root", argv[5], runNum);    
    
//     std::cout<<"Output name is "<<outfile_name<<std::endl;
//     TFile *fpOut = TFile::Open(outfile_name, "RECREATE");
//     if(!fpOut){ std::cerr<<"Cannot open output file "<<fpOut<<std::endl; return -1; }    
    
//     TTree *outTree = new TTree("deconvolvedWaveform", "deconvolvedWaveform");
    
//     double 
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
} //main    