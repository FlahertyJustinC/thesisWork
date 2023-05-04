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
    if(argc<2) {
        std::cout << "Usage\n" << argv[0] << " <input file>  \n";
        return 0;
    }

    TFile *fp = TFile::Open(argv[1]);
    if(!fp) { std::cerr << "Can't open file\n"; return -1; }
    
    // data like
    TTree *eventTree = (TTree*) fp->Get("eventTree");
    if(!eventTree) { std::cerr << "Can't find eventTree\n"; return -1; }
    
    RawAtriStationEvent *rawAtriEvPtr=0;
    eventTree->SetBranchAddress("event",&rawAtriEvPtr);    
    UsefulAtriStationEvent *usefulAtriEvPtr = new UsefulAtriStationEvent(rawAtriEvPtr);
    double weight;
    eventTree->SetBranchAddress("weight", &weight);
    Long64_t numEntries=eventTree->GetEntries();
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    

}