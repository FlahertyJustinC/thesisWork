////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////	analysis.cxx
////
////	January 2019,  clark.2668@osu.edu
////	This is an example of how you analyze ARA data
////	We will learn how to get a waveform
////	How to make a spectrum
////	We will run this as *compiled* code, NOT run as a ROOT macro
////	This code executes over raw (L0) data
////
////	We will also learn how to use the ARA geom tool
////////////////////////////////////////////////////////////////////////////////

// C/C++ Includes
#include <iostream>
#include <cstdlib>
#include <sstream>
#include <fstream>
#include <stdio.h>
#include <vector>
#include <cmath>
#include <algorithm>

//AraRoot Includes
#include "RawAtriStationEvent.h"
#include "UsefulAtriStationEvent.h"
#include "FFTtools.h"
#include "AraGeomTool.h"
#include "AraQualCuts.h"

//ROOT Includes
#include "TTree.h"
#include "TFile.h"
#include "TGraph.h"
#include "TCanvas.h"


using namespace std;

int main(int argc, char **argv)
{
    if(argc<2) {  // Check to make sure there are enough arguments to do something meaningful
		std::cout << "Usage requires you to provide input parameter of the form " << basename(argv[0]) << " <input data file>" << std::endl;
		return -1;
	}
	
	TFile *fpIn = new TFile(argv[1], "OLD"); //we're going to open the data file
	if(!fpIn){
		std::cerr<< "Can not open the old file: " <<argv[1]<<endl;
		return -1;
	} //throw a warning if you can't open it
	
	fpIn->cd(); //go into that file
	TTree *eventTree = (TTree*) fpIn->Get("eventTree"); //load in the event free for this file
	if(!eventTree){
		std::cerr<<"Can't find eventTree in file" <<argv[1]<<" with filename " <<argv[1]<<endl;
		return -1;
	} //throw a warning if you can't open it
	 //set the tree address to access our raw data type
	RawAtriStationEvent *rawAtriEvPtr=0;
	eventTree->SetBranchAddress("event",&rawAtriEvPtr);
	int runNum;
	//we can also get the run number
	eventTree->SetBranchAddress("run",&runNum);
	
	double numEntries = eventTree -> GetEntries(); //get the number of entries in this file

	AraQualCuts *qual = AraQualCuts::Instance();
	numEntries=2;
	for(int event=0; event<numEntries; event++){ //loop over those entries
		
		eventTree->GetEntry(event); //get the event

		//we can see if it has cal pulser or software trigger timing
		bool isCalpulser = rawAtriEvPtr->isCalpulserEvent();
		bool isSoftTrigger = rawAtriEvPtr->isSoftwareTrigger();

		//and also at what unixtime and unixtime microsecond the event was read out by the SBC
		int unixTime=(int)rawAtriEvPtr->unixTime;
		int unixTimeUs=(int)rawAtriEvPtr->unixTimeUs;

		//we can also get the (real) firmware event number the event was stamped with
		int eventNumber=(int)rawAtriEvPtr->eventNumber;
		printf("Converting raw event to Atri event...\n");
		//make a *useful* event out of the *raw* event, which functionally just calibrates it
		UsefulAtriStationEvent *realAtriEvPtr = new UsefulAtriStationEvent(rawAtriEvPtr, AraCalType::kLatestCalib);
        printf("Event converted!\n");
        event->Show()
		vector<string> conditioning = realAtriEvPtr->fConditioningList;
		for(int i=0; i<conditioning.size(); i++){
			cout<<"conditioning value "<<i<<" is "<<conditioning[i]<<endl;
		}

		bool this_qual = qual->isGoodEvent(realAtriEvPtr); //get the quality
		if(!this_qual) continue;
		cout<<"Qual is "<<this_qual<<endl;
	
		//now, we'll get the waveform from channel 2
		TGraph *waveform = realAtriEvPtr->getGraphFromRFChan(0);
		TGraph *interpolated_waveform = FFTtools::getInterpolatedGraph(waveform, 0.5); //get an interpolated waveform with 0.5 ns interpolation
		TGraph *padded_waveform = FFTtools::padWaveToLength(interpolated_waveform,2048);
		TGraph *spectrum = FFTtools::makePowerSpectrumMilliVoltsNanoSecondsdB(padded_waveform); //now make a spectrum

		TCanvas *c = new TCanvas("c","c",2*1000,1000);
		c->Divide(2,2);
		c->cd(1);
			waveform->Draw("ALP");
			waveform->SetTitle("Ch 0 Raw Waveform; Time (ns); Voltage (mV)");
			waveform->SetLineColor(kBlack);
		c->cd(2);
			interpolated_waveform->Draw("ALP");
			interpolated_waveform->SetTitle("Ch 0 Interpolated Waveform; Time (ns); Voltage (mV)");
			interpolated_waveform->SetLineColor(kRed);
		c->cd(3);
			padded_waveform->Draw("ALP");
			padded_waveform->SetTitle("Ch 0 Padded & Interpolated Waveform; Time (ns); Voltage (mV)");
			padded_waveform->SetLineColor(kBlue);
		c->cd(4);
			spectrum->Draw("ALP");
			spectrum->SetTitle("Ch 0 Spectrum; Frequency (MHz); Power (dB)");
			spectrum->SetLineColor(kGreen);
		char title[150];
		sprintf(title,"run%d_event%d.png",runNum,eventNumber);
		c->SaveAs(title);

		//now do some cleanup
		delete c;
		delete spectrum;
		delete interpolated_waveform;
		delete padded_waveform;
		delete waveform;
		delete realAtriEvPtr;
	}


	/*
	Now we can also see how to use the geom tool to get measurement antenna information
	*/
/*
	int station=2;
	AraGeomTool *araGeom = AraGeomTool::Instance();

	for(int i=0; i<16; i++){ //loop over antennas
		int pol = (int) araGeom->getStationInfo(station)->getAntennaInfo(i)->polType; //polarization
		double X = araGeom->getStationInfo(station)->getAntennaInfo(i)->antLocation[0]; //antenna X location
		double Y = araGeom->getStationInfo(station)->getAntennaInfo(i)->antLocation[1]; //antenna Y location
		double Z = araGeom->getStationInfo(station)->getAntennaInfo(i)->antLocation[2]; //antenna Z location
		double delay = araGeom->getStationInfo(station)->getCableDelay(i); //the associated cable delay
	}

	/*
	Now we can also see how to use the geom tool to get cal-pulser antenna information
	*/
/*
	for(int i=0; i<araGeom->getStationInfo(station)->getNumCalAnts(); i++){ //loop over number of cal antennas
		double X = araGeom->getStationInfo(station)->getCalAntennaInfo(i)->antLocation[0];
		double Y = araGeom->getStationInfo(station)->getCalAntennaInfo(i)->antLocation[1];
		double Z = araGeom->getStationInfo(station)->getCalAntennaInfo(i)->antLocation[2];
		string locName(&araGeom->getStationInfo(station)->getCalAntennaInfo(i)->locationName[0]);
		string antName(araGeom->getStationInfo(station)->getCalAntennaInfo(i)->getCalAntName());
	}
*/
printf("End of code reached.\n");
}//close the main program