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
UsefulAtriStationEvent *usefulAtriEvPtr;

// AraSim includes
#include "/cvmfs/ara.opensciencegrid.org/trunk/centos7/source/AraSim/Detector.h"
#include "/cvmfs/ara.opensciencegrid.org/trunk/centos7/source/AraSim/Event.h"
#include "/cvmfs/ara.opensciencegrid.org/trunk/centos7/source/AraSim/Position.h"
#include "/cvmfs/ara.opensciencegrid.org/trunk/centos7/source/AraSim/RaySolver.h"
#include "/cvmfs/ara.opensciencegrid.org/trunk/centos7/source/AraSim/IceModel.h"

// #include "/users/PAS0654/jflaherty13/source/AraSim/Detector.h"
// #include "/users/PAS0654/jflaherty13/source/AraSim/Event.h"
// #include "/users/PAS0654/jflaherty13/source/AraSim/Position.h"
// #include "/users/PAS0654/jflaherty13/source/AraSim/RaySolver.h"
// #include "/users/PAS0654/jflaherty13/source/AraSim/IceModel.h"

#include "helper.h"

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

double radii[] = {
      0,  150,  300,  450,  600,  750,  900, 1050, 1200, 1350, 1500, 1650, 1800, 1950,
 2100, 2250, 2400, 2550, 2700, 2850, 3000, 3150, 3300, 3450, 3600, 3750, 3900, 4050,
 4200, 4350, 4500, 4650, 4800, 4950
};
const int numScanned = 34;

// double radii[] = {
//     1500, 3300
// };
// const int numScanned = 2;

int main(int argc, char **argv)
{
    if(argc<5) {
        std::cout << "Usage\n" << argv[0] << " <station> <runnum> <input file> <output_dir> \n";
        std::cout << "e.g.\n" << argv[0] << " 2 http://www.hep.ucl.ac.uk/uhen/ara/monitor/root/run1841/event1841.root\n";
        return 0;
    }

    double interpV = 0.4;
    double interpH = 0.625;

    int station = atoi(argv[1]);
    int runNum = atoi(argv[2]);

    // set up correlator and pairs
    RayTraceCorrelator *theCorrelators[numScanned];
    for(int r=0; r<numScanned; r++){
        
        // setup the paths to our ray tracing tables
        double radius = radii[r];
        printf("Loading Radius %.2f\n",radius);
        double angular_size = 1.;
        int iceModel = 0;
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

        int numAntennas = 16;
        theCorrelators[r] = new RayTraceCorrelator(station, numAntennas, radius, angular_size, dirPath, refPath);
        theCorrelators[r]->LoadTables();        
    }
    
    AraGeomTool *geomTool = AraGeomTool::Instance();
    std::vector<int> excludedChannels = {15};
    std::map< int, std::vector<int> > pairs_V = theCorrelators[0]->SetupPairs(station, geomTool, AraAntPol::kVertical, excludedChannels);
    std::map< int, std::vector<int> > pairs_H = theCorrelators[0]->SetupPairs(station, geomTool, AraAntPol::kHorizontal, excludedChannels);

    printf("------------------\n");
    printf("Correlator Setup Complete. Make Output Files\n");
    printf("------------------\n");

    char outfile_name[400];
    sprintf(outfile_name, "%s/recangle_reco_out_run_%d.root", argv[4], runNum);
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

    double weight_out;
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

    printf("------------------\n");
    printf("Output File Setup Complete. Begin looping events\n");
    printf("------------------\n");

    printf("Opening file...\n");
    TFile *fp = TFile::Open(argv[3]);
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
    
    // simulation
    TTree *simSettingsTree;
    Detector *detector = 0;
    simSettingsTree=(TTree*) fp->Get("AraTree");
    if(!simSettingsTree) { std::cerr << "Can't find AraTree\n"; return -1; }
    simSettingsTree->SetBranchAddress("detector", &detector);
    simSettingsTree->GetEntry(0);
    std::vector<double> average_position = get_detector_cog(detector);
    printf("Detector Center %.2f, %.2f, %.2f \n", average_position[0], average_position[1], average_position[2]);

    
    TTree *simTree;
    Event *eventPtr = 0; // it is apparently incredibly important that this be initialized to zero...
    Report *reportPtr = 0;
    simTree=(TTree*) fp->Get("AraTree2");
    if(!simTree) { std::cerr << "Can't find AraTree2\n"; return -1; }
    simTree->SetBranchAddress("event", &eventPtr);
    simTree->SetBranchAddress("report", &reportPtr);

    printf("------------------\n");
    printf("Input files loaded. Set up ray tracing business\n");
    printf("------------------\n");
    
    RaySolver *raySolver = new RaySolver;
    IceModel *iceModel = new IceModel(0 + 1*10, 0, 0);
    Settings *settings = new Settings();
    // configure settings    
    settings->Z_THIS_TOLERANCE = 1; // higher tolerance
    settings->Z_TOLERANCE = 0.05;
    settings->NOFZ=1; // make sure n(z) is turned on
    settings->RAY_TRACE_ICE_MODEL_PARAMS = 0; // set the ice model as user requested
    
    for(Long64_t event=0;event<numEntries;event++) {
        fp->cd();
        eventTree->GetEntry(event);
        simTree->GetEntry(event);

        std::cout<<"Looking at event number "<<event<<std::endl;
        
        // Debugging - JCF 6/7/2023
        // if (event!=0){
        //     continue;
        // }
        
        // if (rawAtriEvPtr->eventNumber != 0){
        //     continue;
        // }
        // End debugging - JCF 6/7/2023

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

        Position diff_true;
        diff_true.SetXYZ(difference[0], difference[1], difference[2]);
        printf("  Vertex Theta %.2f, Phi %.2f, R %.2f\n", diff_true.Theta()*TMath::RadToDeg(), diff_true.Phi()*TMath::RadToDeg(), diff_true.R());


        /*
            First, we need to do correlation maps to get directions
        */

        std::map<int, TGraph*> interpolatedWaveforms;
        for(int i=0; i<16; i++){
            TGraph *gr = usefulAtriEvPtr->getGraphFromRFChan(i);
            TGraph *grInt = FFTtools::getInterpolatedGraph(gr,i<8?interpV:interpH);
            interpolatedWaveforms[i] = grInt;
            delete gr;
        }

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
            double peak_min = abs(TMath::MinElement(interpolatedWaveforms[i]->GetN(), interpolatedWaveforms[i]->GetY()));
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
            // TCanvas *c = new TCanvas("","", 1200, 950);            
            // maps[0]->Draw("colz");
            // maps[0]->GetXaxis()->SetTitle("Phi [deg]");
            // maps[0]->GetYaxis()->SetTitle("Theta [deg]");
            // maps[0]->GetYaxis()->SetTitleSize(0.05);
            // maps[0]->GetYaxis()->SetLabelSize(0.03);
            // maps[0]->GetYaxis()->SetTitleOffset(0.6);
            // maps[0]->GetXaxis()->SetTitleSize(0.05);
            // maps[0]->GetXaxis()->SetLabelSize(0.03);
            // maps[0]->GetXaxis()->SetTitleOffset(0.6);
            // gStyle->SetOptStat(0);
            // maps[0]->GetXaxis()->CenterTitle();
            // maps[0]->GetYaxis()->CenterTitle();
            // gPad->SetRightMargin(0.15);
            // char title[500];
            // sprintf(title,"debuggingPlots/maps_ev%d_rad%.2f.png", event, radii[r]);
            // c->SaveAs(title);
            // delete c;
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
            
            // // fix the theta range
            // if(peakTheta < 0){
            //     peakTheta = abs(peakTheta) + 90;
            // }
            // else if(peakTheta>0){
            //     peakTheta = 90 - peakTheta;
            // }
            // else{
            //     peakTheta+=90;
            // }

            // if(peakPhi < 0) peakPhi =360.0 - abs(peakPhi);

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
        // printf("  Best Radius is %.2f, with Corr %.4f, Theta %.2f, Phi %.2f \n", 
        //     radii[element], peakCorrs[element],
        //     peakThetas[element], peakPhis[element]
        //     );
        
        // double temp_bestTheta = peakThetas[element];
        // double temp_bestPhi = peakPhis[element];
        // int temp_binTheta, temp_binPhi;
        // theCorrelators[element]->ConvertAngleToBins(temp_bestTheta, temp_bestPhi, temp_binTheta, temp_binPhi);
        // double temp_arrivalTheta, temp_arrivalPhi;
        // int temp_bestSol = peakSol[element];
        // theCorrelators[element]->LookupArrivalAngles(0, temp_bestSol, 
        //     temp_binTheta, temp_binPhi,
        //     temp_arrivalTheta, temp_arrivalPhi
        // );

        // rerange_theta_phi(temp_bestTheta, temp_bestPhi);

        // printf("  Best Vertex Theta/Phi: %.2f/%.2f. Corresponding Arrival Theta/Phi %.2f/%.2f, and Bin %d/%d \n",
        //     temp_bestTheta, temp_bestPhi,
        //     TMath::RadToDeg()*temp_arrivalTheta, TMath::RadToDeg()*temp_arrivalPhi,
        //     temp_binTheta, temp_binPhi
        // );

        // int likely_sol = guess_triggering_solution(eventPtr, reportPtr);
        int likely_sol = 0; //Forcing direct solutions only for debugging - JCF 7/5/2023
        std::map<int, double> thetas_truth = get_value_from_mc_truth("theta", likely_sol, reportPtr);
        std::map<int, double> phis_truth = get_value_from_mc_truth("phi", likely_sol, reportPtr);
        // for(int i=0; i<16; i++){
        //     printf("    I %d, True Theta/Phi is %.2f, %.2f\n",i, 
        //         TMath::RadToDeg()*thetas_truth.find(i)->second,
        //         TMath::RadToDeg()*phis_truth.find(i)->second)
        //         ;
        // }

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
            // printf("  Ant %d, Reco Arrival Theta %.2f, Reco Arrival Phi %.2f \n",
            //     i, this_arrivalTheta, this_arrivalPhi
            // );
            reco_arrivalThetas_out[i] = this_arrivalTheta;
            reco_arrivalPhis_out[i] = this_arrivalPhi; 
        }

        // then the true quantities
        for(int i=0; i<16; i++){
            double this_true_theta = thetas_truth.find(i)->second;
            double this_true_phi = phis_truth.find(i)->second;
            // printf("  Ant %d, True Arrival Theta %.2f, Reco Arrival Phi %.2f \n",
            //     i, this_true_theta, this_true_phi
            // );
            true_arrivalThetas_out[i] = this_true_theta;
            true_arrivalPhis_out[i] = this_true_phi; 
        }
        trueTheta_out = diff_true.Theta() * TMath::RadToDeg();
        truePhi_out = diff_true.Phi() * TMath::RadToDeg();
        trueR_out = diff_true.R();
        likelySol_out = likely_sol;
        
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