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

int main (int argc, char **argv) {
    
    if(argc<3) {
        std::cout << "Usage\n" << argv[0] << " <input file> <event number> \n";
        std::cout << "e.g.\n" << argv[0] << " AraOut.root 0\n";
        return 0;
    }   
    
    char* inputFile = argv[1];
    int eventNumber = atoi(argv[2]);  

    printf("Opening file...\n");
    TFile *fp = TFile::Open(inputFile);
    if(!fp) { std::cerr << "Can't open file\n"; return -1; }
    printf("File opened!\n");    
    
    //Import eventTree
    TTree *eventTree = (TTree*) fp->Get("eventTree");
    printf("Event tree opened!\n");
    if(!eventTree) { std::cerr << "Can't find eventTree\n"; return -1; }
    Long64_t numEntries=eventTree->GetEntries();
    cout << "eventTree has " << numEntries << " entries." << endl;
    RawAtriStationEvent *rawAtriEvPtr=0; 
    
    
    Settings *settings1 = new Settings();
    // settings1->ReadFile(setupfile); 
    IceModel *icemodel=new IceModel(settings1->ICE_MODEL + settings1->NOFZ*10,settings1->CONSTANTICETHICKNESS * 1000 + settings1->CONSTANTCRUST * 100 + settings1->FIXEDELEVATION * 10 + 0,settings1->MOOREBAY);// creates Antarctica ice model
    Detector *detector = new Detector();  
    Report *report = new Report(detector, settings1);
    RaySolver *raySolver = new RaySolver;    
    
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
        // weight = 1;
    }
    // sim like
    else {
        std::cerr << "AraTree exists.  Importing as simulated data.\n";
        eventTree->SetBranchAddress("UsefulAtriStationEvent", &usefulAtriEvPtr);
        // weight;
        // eventTree->SetBranchAddress("weight", &weight);
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
    
 
    for(Long64_t event=0;event<numEntries;event++) {
        fp->cd();
        eventTree->GetEntry(event);
        if (dataLike) {
            UsefulAtriStationEvent *usefulAtriEvPtr = new UsefulAtriStationEvent(rawAtriEvPtr, AraCalType::kLatestCalib);
        }
        else {
            simTree->GetEntry(event);
        }    
    
        TGraph *g[16];
        auto *c = new TCanvas();
        c->Divide(4,4);
        //Testing waveform plotting
        for(int i=0; i<16; i++){
            TGraph *gr = usefulAtriEvPtr->getGraphFromRFChan(i);                  
            c->cd(i+1); gPad->SetGrid(1,1);
            gr->Draw();
        }
        char title[500];
        sprintf(title, "waveform.png");   
        c->Print(title);
    
    }
    
}