#include <fstream>
#include <iostream>
#include <stdio.h>
#include <cstdlib>
#include <sstream>
#include <functional>
#include <numeric>
#include <deque>
using namespace std;

//ROOT Includes                                                                                                                 
#include "TROOT.h"
#include "TCanvas.h"
#include "TTree.h"
#include "TFile.h"
#include "TTree.h"
#include "math.h"
#include "TF1.h"
#include "TH2F.h"

#include "TGraph.h"
#include "TMath.h"
#include "TFitResult.h"
#include "TF1.h"
#include "Math/Interpolator.h"
#include "Math/InterpolationTypes.h"

#include "FFTtools.h"
#include "Tools.h"
#include "helper.h"

double integrateBinPower( TGraph *plot, int numBinsToIntegrate, vector<double> &integratedBins)
{
  int nPoints = plot->GetN();
  if (nPoints < numBinsToIntegrate){
    return 0;
  }

  //  Double_t *xVals = plot->GetX();                                                                                                          
  Double_t *yVals = plot->GetY();
  std::deque<double> integrator;
  double sum = 0.;
  for (int i = 0; i < numBinsToIntegrate; i++){
    integrator.push_back(pow(yVals[i], 2));
    sum = accumulate(integrator.begin(), integrator.end(), 0);
  }
  double max = 0.;
  integratedBins.push_back(sum);

  for (int i = 0+numBinsToIntegrate; i < nPoints; i++){

    sum = sum - integrator[0];
    integrator.pop_front();
    integrator.push_back(pow(yVals[i], 2));
    sum = sum + integrator[numBinsToIntegrate-1];

    integratedBins.push_back(sum);

    if (sum > max){
      max = sum;
    }
  }

  return max;
}

vector<TGraph*> makeIntegratedBinPowerGraphs(vector<TGraph*> graphsInput, int numBinsToIntegrate, string xlabel, string ylabel, vector<string> titles){

  //  cout << graphsInput.size() << endl;
  vector<TGraph*> graphsOutput;

  for (int i = 0; i < graphsInput.size(); i++){
    vector<double> integratedBins;
    double maxIntPower = integrateBinPower(graphsInput[i], numBinsToIntegrate, integratedBins);
    double *volts = &integratedBins[0];
    TGraph* gIntPower = new TGraph(integratedBins.size(), graphsInput[i]->GetX(), volts);
    gIntPower->GetXaxis()->SetTitle(xlabel.c_str());
    gIntPower->GetYaxis()->SetTitle(ylabel.c_str());
    gIntPower->SetTitle(titles[i].c_str());
    graphsOutput.push_back(gIntPower);
    integratedBins.clear();
    //    delete gIntPower;
  }

  return graphsOutput;
}

void getAbsMaximum_N(TGraph* plot, int nPeaks, double timeApart, vector<double> &xs, vector<double> &ys)
{
  xs.clear();
  ys.clear();

  int nPoints = plot->GetN();
  if (nPoints < nPeaks){
    cerr << "Number of points in waveform is fewer than the number of peaks requested. Decreasing number of peaks requested to match number points." << endl;
    nPeaks = nPoints;
  } 

  double x_temp, y_temp;
  double y_good, x_good;
  int test;
  double y_upper;

  for (int iPeak = 0; iPeak < nPeaks; iPeak++){
    y_good = -9.0E99;
    y_upper = 9.0E99;
    if (iPeak > 0){
      y_upper = ys[iPeak-1];
    }
    for (int i = 0; i < nPoints; i++){
      test = plot->GetPoint(i, x_temp, y_temp);
      if (iPeak == 0){
        if (y_temp > y_good){
          x_good = x_temp;
          y_good = y_temp;
        }
      } 
      else {
        for (int iTimeTest = 0; iTimeTest < xs.size(); iTimeTest++){
          if (y_temp > y_good && y_temp < y_upper && abs(x_temp - xs[iTimeTest]) > timeApart){
            x_good = x_temp;
            y_good = y_temp;
          }
        }
      }
    }
    xs.push_back(x_good);
    ys.push_back(y_good);
    //cout << iPeak << " : " << xs[iPeak] << " : " << ys[iPeak] << endl;
  }
  return;
}

void getAbsMaximum_N(vector<TGraph*> graphs, int nPeaks, double timeApart, vector<vector<double> > &xs, vector<vector<double> > &ys){

  //  vector<double> xs;
  // vector<double> ys;

  xs.clear();
  ys.clear();

  vector<double> xs_temp;
  vector<double> ys_temp;

  for (int i = 0; i < graphs.size(); i++){
    getAbsMaximum_N(graphs[i], nPeaks, timeApart, xs_temp, ys_temp);
    xs.push_back(xs_temp);
    ys.push_back(ys_temp);
  }
  
  //    cout << ys.size() << endl;
  //  cout << ys[0].size() << endl;
}



void getExcludedChannels(std::vector<int> &excludedChannels, Settings *settings1, Detector *detector) {
    if (settings1->DETECTOR_STATION==2){
        // cout << "Using station 2.  Excluding channels 1 and 3 in addition to the masking." << endl;
        // excludedChannels.push_back(1);  //Removing extra A2 channels as the D-pulse gets clipped.
        // excludedChannels.push_back(3); 
        // excludedChannels.push_back(5); 
        // excludedChannels.push_back(6); 
    }    
    if (settings1->DETECTOR_STATION==3){
        // cout << "Using station 3.  Excluding channel 1, 3 in addition to the masking." << endl;
        // // excludedChannels.push_back(0);
        // excludedChannels.push_back(1);
        // excludedChannels.push_back(3);
    }    
    if (settings1->DETECTOR_STATION==4){
        // cout << "Using station 4.  Excluding channels 0 in addition to the masking." << endl;
        // excludedChannels.push_back(0);
    }
    if (settings1->DETECTOR_STATION==5){
        // cout << "Using station 5.  Excluding channels 2, 6, 10, 14 in addition to the masking." << endl;
        // excludedChannels.push_back(2);
        // excludedChannels.push_back(6); 
        // excludedChannels.push_back(10);
        // excludedChannels.push_back(14);
    }    
    for (int i=0; i<16; i++){
        // cout << "detector->GetTrigMasking(i) = " << detector->GetTrigMasking(i) << endl;
        if (not detector->GetTrigMasking(i)){
            // cout << "Excluding channel " << i << endl;
            excludedChannels.push_back(i);
        }
    }
    cout << "Excluded Channels: ";
    for (int i: excludedChannels) {
        cout << i << ", ";
    }
    cout << endl;
}

double calculateNoiseRMS(TGraph *gr, int sampleNumber=100) {
    int totalSampleNumber = abs(sampleNumber);
    int waveformLength = gr->GetN();
    double voltageSubset[totalSampleNumber];
    if (sampleNumber>0){
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
    
    // cout << "noiseRms = " << noiseRms << endl;
    
    if (noiseRms == 0) {
        noiseRms=1;
    }
    
    return noiseRms;
}

//Creating CSW function
void calculateCSW(UsefulAtriStationEvent *usefulAtriEvPtr, std::vector<int> excludedChannels, double cropTimes[16], double arrivalTimes[16], TGraph *cswVpolOut, TGraph *cswHpolOut, char* run, double sampleNs=80, double dt=0.5, bool debugMode=false) {
    //Harc-coding interpolation interval for now.
    // const double dt = 0.5;

    //Loop over excluded channels to form excluded pairs list
    std::vector<int> cswExcludedChannelPairs;
    if (not excludedChannels.empty()) {
        for(int i=0; i<8; i++){
            bool checkVpol = std::find(excludedChannels.begin(), excludedChannels.end(), i) != excludedChannels.end();
            bool checkHpol = std::find(excludedChannels.begin(), excludedChannels.end(), i+8) != excludedChannels.end();
            if (checkVpol or checkHpol) { 
                cswExcludedChannelPairs.push_back(i);
            }
            //Add check if arrival time is -1000, so that channel gets excluded
            else if (arrivalTimes[i] == -1000 or arrivalTimes[i+8] == -1000) {
                cswExcludedChannelPairs.push_back(i);
            }

        }
    }
    else {cout << "Excluded Channel array empty!" << endl;}


    //Make map of kept Vpol and kept HPol channels
    map<int, TGraph*> mapWaveforms;
    map<int, double> mins; //Map of the min time for each waveform.
    map<int, double> maxes; // Map of the max time sfor ecah waveform.
    double tMin=-1e100;  //Will represent the LARGEST min time across channels.
    double tMax=1e100;  //Will represent the SMALLEST max time across channels.

    double arrivalTimeMax = *max_element(arrivalTimes, arrivalTimes + 16);
    double arrivalTimeMin = *min_element(arrivalTimes, arrivalTimes + 16);
    double arrivalTimeAvg = 0;
    for (int i=0; i<16; i++) {
        if (arrivalTimes[i] == -1000) {
            cout << "Non-physical arrival times.  Bypassing event." << endl;
            return;
        }        
        arrivalTimeAvg += arrivalTimes[i]/16;
    }   
    // if (arrivalTimeAvg == -1000) {
    //     cout << "Non-physical arrival times.  Bypassing event." << endl;
    //     return;
    // }

    
    //Import waveforms and Apply time delays
    TGraph *gr;
    for (int i=0; i<16; i++) {
        
        bool checkExcluded = std::find(cswExcludedChannelPairs.begin(), cswExcludedChannelPairs.end(), i%8) != cswExcludedChannelPairs.end();
        if (not checkExcluded) {        
            //Import waveform as Tgraph.
            gr = usefulAtriEvPtr->getGraphFromRFChan(i);


            if (cropTimes[i] < gr->GetX()[0]+sampleNs) {
                // if (debugMode) {cout << "Bypassing channel" << endl;}
                cswExcludedChannelPairs.push_back(i%8);
                continue;
            }
            //Crop to D-pulse only.
            gr = FFTtools::cropWave(gr, gr->GetX()[0], cropTimes[i]);

            //Interpolate to 0.5 ns
            gr = FFTtools::getInterpolatedGraph(gr,dt);        

            for(int k=0; k<gr->GetN(); k++){
                gr->GetX()[k] = gr->GetX()[k] - arrivalTimes[i] + arrivalTimeAvg;
            }
            mapWaveforms[i] = gr;
            mins[i] = gr->GetX()[0];
            maxes[i] = gr->GetX()[gr->GetN()-1];

            // bool checkExcluded = std::find(cswExcludedChannelPairs.begin(), cswExcludedChannelPairs.end(), i%8) != cswExcludedChannelPairs.end();
            // if (not checkExcluded) {
            if (tMin < mins[i] and mins[i] < tMax) {tMin = mins[i];}
            if (tMax > maxes[i] and maxes[i] > tMin) {tMax = maxes[i];}
        }
        


        // else {
        //     cout << "Excluding channel " << i << endl;
        // }
    }
    

    
    //Extra padding of tMin and tMin to get it away from the edges
    tMin += 2*dt;
    tMax -= 2*dt;    
    
    
    if (tMin > tMax) {
        cout << "WARNING: tMin > tMax! Run, Event = " << run << ", " << usefulAtriEvPtr->eventNumber << endl;
        return;
    }
    
    if ((tMax-tMin) < 80) {
        cout << "Waveform too small for CSW.  Bypassing event. " << endl;
        return;
    }    

    
    //Initialize time and voltage arrays for the coherent sum.
    const int nT = int((tMax-tMin)/dt);
    double timeCsw[nT];
    double voltageCswVpol[nT];    
    double voltageCswHpol[nT];      
    

    //Loop over waveforms and crop to the global tMin and tMax
    for (int i=0; i<16; i++) {

        bool checkExcluded = std::find(cswExcludedChannelPairs.begin(), cswExcludedChannelPairs.end(), i%8) != cswExcludedChannelPairs.end();

        if (not checkExcluded) {
            mapWaveforms[i] = FFTtools::cropWave(mapWaveforms[i], tMin, tMax);  //Seg fault is here!
            
            if (mapWaveforms[i]->GetN() == 0) {
                if (debugMode) {cout << "Bypassing empty channel." << endl;}
                cswExcludedChannelPairs.push_back(i);
                continue;
            }


            //Interpolate to same time domain.
            std::vector<double> tVec;
            std::vector<double> vVec;

            Int_t numIn=mapWaveforms[i]->GetN();
            // cout << "numIn = " << numIn << endl;
            Double_t tIn,vIn;       

            for (int samp=0;samp<numIn;samp++) {
                // mapWaveforms[i]->GetPoint(samp,tIn,vIn);
                tVec.push_back(mapWaveforms[i]->GetX()[samp]);
                vVec.push_back(mapWaveforms[i]->GetY()[samp]);
            }

            // cout << "tVec[0] = " << tVec[0] << endl;
            // cout << "vVec[0] = " << vVec[0] << endl;

            ROOT::Math::Interpolator chanInterp(tVec,vVec,ROOT::Math::Interpolation::kAKIMA);

            Int_t roughPoints=Int_t((tMax-tMin)/dt);

            Double_t *newTimes = new Double_t[roughPoints]; //Will change this at some point, but for now
            Double_t *newVolts = new Double_t[roughPoints]; //Will change this at some point, but for now
            Int_t numPoints=0;
            for(Double_t time=tMin+dt;time<tMax-dt;time+=dt) {
              newTimes[numPoints]=time;
              newVolts[numPoints]=chanInterp.Eval(time);  
              numPoints++;
            }
            // cout << "newTimes[0] = " << newTimes[0] << endl;
            // cout << "newVolts[0] = " << newVolts[0] << endl;
            mapWaveforms[i] = new TGraph(numPoints,newTimes,newVolts);
            // cout << "mapWaveforms[i]->GetN() = " << mapWaveforms[i]->GetN() << endl;
            // delete [] newTimes;
            // delete [] newVolts;  
            
        }
    }
    
    //Add condition to check if all channel pairs are excluded
    if (cswExcludedChannelPairs.size() < 8) {
    
        //Now loop over the time bins to sum the waveforms
        for (int j=0; j<nT-1; j++) {
            //Initilize time bin
            timeCsw[j] = tMin + j*dt;
            // //Initialize voltage bin
            double thisVpolVoltage=0;
            double thisHpolVoltage=0;
            for (int i=0; i<8; i++) {
                bool checkExcluded = std::find(cswExcludedChannelPairs.begin(), cswExcludedChannelPairs.end(), i%8) != cswExcludedChannelPairs.end();
                if (not checkExcluded) {
                    thisVpolVoltage+=mapWaveforms[i]->GetY()[j];
                    thisHpolVoltage+=mapWaveforms[i+8]->GetY()[j];
                }
            }
            cswVpolOut->SetPoint(j, timeCsw[j], thisVpolVoltage);
            cswHpolOut->SetPoint(j, timeCsw[j], thisHpolVoltage);
        }
    }
    else {
        if (debugMode){cout << "No passing channels for CSW.  Return empty array" << endl;}
        
    }
    
}

void createFFTarrays(double* voltage, double* time, Settings *settings1, double* V_forfft, double* T_forfft, int waveform_bin, int Nnew){
    double dT_forfft = time[1] - time[0];

    for (int n = 0; n < Nnew; n++)
    {
        T_forfft[n] = time[waveform_bin / 2] - (dT_forfft *(double)(Nnew / 2 - n));

        if ((n >= Nnew / 2 - waveform_bin / 2) &&
            (n < Nnew / 2 + waveform_bin / 2))
        {
            V_forfft[n] = voltage[n - (Nnew / 2 - waveform_bin / 2)];
        }
        else
            V_forfft[n] = 0.;
    }         
}

void wienerDeconvolution(double &vm_real, double &vm_imag, double corr_real, double corr_imag, double psdSnr) {
    double corr_magnitudeSquared = corr_real*corr_real + corr_imag*corr_imag;
    
    double secondTerm = (1 / (1 + corr_magnitudeSquared/psdSnr));
    
    // cout << "psdSnr = " << psdSnr << endl;
    // cout << "corr_real = " << corr_real << endl;
    // cout << "corr_imag = " << corr_imag << endl;
    // cout << "corr_magnitudeSquared = " << corr_magnitudeSquared << endl;
    // cout << "secondTerm = " << secondTerm << endl;
    // cout << "corr_real*secondTerm = " << corr_real*secondTerm << endl;
    // cout << "corr_imag*secondTerm = " << corr_imag*secondTerm << endl;
    
    // double realOut = (vm_real*corr_real - vm_imag*corr_imag)*secondTerm;
    // double imagOut = (vm_imag*corr_real + vm_real*corr_imag)*secondTerm;
    
    double realOut = vm_real*secondTerm;
    double imagOut = vm_imag*secondTerm;
    
    if (std::isinf(corr_magnitudeSquared) or std::isnan(corr_magnitudeSquared)) {
        vm_real = 0;
        vm_imag = 0;        
    }
    else {
        vm_real = realOut;
        vm_imag = imagOut;
    }
}

void wienerDeconvolution2(double &vm_real, double &vm_imag, double corr_real, double corr_imag, double psdSnr) {
    double corr_magnitudeSquared = corr_real*corr_real + corr_imag*corr_imag;
    
    double secondTerm = (1 / (1 + 1/(corr_magnitudeSquared*psdSnr)));
    
    // cout << "corr_real = " << corr_real << endl;
    // cout << "corr_imag = " << corr_imag << endl;
    // cout << "secondTerm = " << secondTerm << endl;
    
    double realOut = (vm_real*corr_real + vm_imag*corr_imag)/corr_magnitudeSquared*secondTerm;
    double imagOut = (vm_imag*corr_real - vm_real*corr_imag)/corr_magnitudeSquared*secondTerm;
    
    if (std::isinf(1/corr_magnitudeSquared) or std::isnan(1/corr_magnitudeSquared) or std::isnan(realOut*realOut+imagOut*imagOut)) {
        vm_real = 0;
        vm_imag = 0;        
    }
    else {
        vm_real = realOut;
        vm_imag = imagOut;
    }
}

TGraph* resizeForFFT(TGraph *gr, int NFOUR) {
    if (gr->GetN() < NFOUR/2) {
        gr = FFTtools::padWaveToLength(gr, NFOUR/2);
    }
    else if (gr->GetN() > NFOUR/2) {
        gr = FFTtools::cropWave(gr, gr->GetX()[0], gr->GetX()[NFOUR/2-1]);
    }
}

void getPowerSpectrumSNR(TGraph *gr, double tNoiseMin, double tNoiseMax, double directPeakTime, double refPeakTime, int waveform_bin, int Nnew, Settings *settings1, double* snrPsd, double* freqPsd, double dt=0.5, double sampleNs=80, bool debugMode=false) {
    //Initialize TGraphs for noise and signal
    TGraph *grNoise;
    TGraph *grSignal;
    
    // double tNoiseMin;
    // double tNoiseMax;
    double tSignalMin = directPeakTime-0.25*sampleNs;
    double tSignalMax = directPeakTime+0.75*sampleNs;
    
//     if (directPeakTime < gr->GetX()[int(sampleNs/dt)]) {
//         // cout << "aaa" << endl;
//         // grNoise = FFTtools::cropWave(gr, gr->GetX()[gr->GetN()-1] - sampleNs, gr->GetX()[gr->GetN()-1]);
//         tNoiseMax = gr->GetX()[gr->GetN()-1]-dt;
//         tNoiseMin = tNoiseMax - sampleNs;
//         // tNoiseMin = refPeakTime + sampleNs;
//     }
//     else {
//         // cout << "bbb" << endl;
//         // grNoise = FFTtools::cropWave(gr, gr->GetX()[0], gr->GetX()[int(sampleNs/dt)-1]);
//         tNoiseMin = gr->GetX()[0]+dt;
//         tNoiseMax = tNoiseMin + sampleNs;
//         // tNoiseMax = directPeakTime - sampleNs;

//     }        
    
    //Check if peakTime is too close to beginning or end of waveform.
    if (tSignalMin < gr->GetX()[0]) {
        tSignalMin = gr->GetX()[0]+dt;
        tSignalMax = tSignalMin+sampleNs;
            
    }
    else if (tSignalMax > gr->GetX()[gr->GetN()-1]) {
        tSignalMax = gr->GetX()[gr->GetN()-1]-dt;
        tSignalMin = tSignalMax-sampleNs;
        
    }    

    grNoise = FFTtools::cropWave(gr, tNoiseMin, tNoiseMax);
    grSignal = FFTtools::cropWave(gr, tSignalMin, tSignalMax);

    
    //Pad waveforms to match the frequency binning used in the deconvolution
    grSignal = resizeForFFT(grSignal, settings1->NFOUR);
    grNoise = resizeForFFT(grNoise, settings1->NFOUR);
    // grSignal = FFTtools::padWaveToLength(grSignal, settings1->NFOUR/2);
    // grNoise = FFTtools::padWaveToLength(grNoise, settings1->NFOUR/2);
        
    //Make PSD for signal and noise sample
    TGraph *grPsdSignal = FFTtools::makePowerSpectrumMilliVoltsNanoSeconds(grSignal);
    TGraph *grPsdNoise = FFTtools::makePowerSpectrumMilliVoltsNanoSeconds(grNoise);
    
    //Push PSD to vectors for interpolation into the frequency binning used in AraSim
    std::vector<double> psdSignal;
    std::vector<double> psdNoise;
    
    std::vector<double> freqSignal;
    std::vector<double> freqNoise;    
    
    std::vector<double> snrVec;
    std::vector<double> freqVec;  
    
    //Loop over Tgraphs and push PSD's to vectors
    for (int i=0; i<grPsdSignal->GetN(); i++){
        
        double this_freqSignal = grPsdSignal->GetX()[i];        
        double this_psdSignal = grPsdSignal->GetY()[i];
        
        double this_freqNoise = grPsdNoise->GetX()[i]; 
        double this_psdNoise = grPsdNoise->GetY()[i]; 
        
        psdSignal.push_back(this_psdSignal);
        freqSignal.push_back(this_freqSignal);
        
        psdNoise.push_back(this_psdNoise);
        freqNoise.push_back(this_freqNoise);          
    }
    
    //Create interpolators of the noise and signal vectors
    ROOT::Math::Interpolator signalInterp(freqSignal,psdSignal,ROOT::Math::Interpolation::kAKIMA);
    ROOT::Math::Interpolator noiseInterp(freqNoise,psdNoise,ROOT::Math::Interpolation::kAKIMA);
    
    //Get frequency interval for interpolation
    double df = freqSignal[1]-freqSignal[0];
    
    //Interpolate signal and noise and calculate SNR.
    for (int k=0; k<Nnew/2; k++) {
        freqPsd[k] = freqSignal[k] + df/2;
        double this_interpSignal = signalInterp.Eval(freqPsd[k]);
        double this_interpNoise = noiseInterp.Eval(freqPsd[k]);
        if (this_interpNoise == 0) {
            snrPsd[k] = this_interpSignal/1;
        }
        else {
            snrPsd[k] = this_interpSignal/this_interpNoise;
        }
    }  
}

int signum(double value) {
    return (0 < value) - (value < 0);
}

double peakHilbert(UsefulAtriStationEvent *usefulAtriEvPtr, double* hilbertPeakOut, double* peakTimeOut, double* cutoffTime, double dt, double toleranceNs=20, bool findPolarity=false, bool debugMode=false) {
    //Loop over channel pairs and import VPol and HPol waveforms.
    for (int channelPair=0; channelPair<8; channelPair++) {
        //Import waveforms
        TGraph *grV = usefulAtriEvPtr->getGraphFromRFChan(channelPair);
        TGraph *grH = usefulAtriEvPtr->getGraphFromRFChan(channelPair+8);
        //Interpolate to desired dt
        grV = FFTtools::getInterpolatedGraph(grV,dt);
        grH = FFTtools::getInterpolatedGraph(grH,dt);
        //Crop waves to their cutoff time to maintain Direct solution
        //Add condition where waveform gets cropped only if cutoff time is greater than the zero time of both waveforms.
        if (cutoffTime[channelPair] > grV->GetX()[0] and cutoffTime[channelPair+8] > grH->GetX()[0]) {
            grV = FFTtools::cropWave(grV, grV->GetX()[0], cutoffTime[channelPair]);
            grH = FFTtools::cropWave(grH, grH->GetX()[0], cutoffTime[channelPair+8]);
        }
        //Get Hilbert envelope of each waveform
        TGraph *hilbertV = FFTtools::getHilbertEnvelope(grV);
        TGraph *hilbertH = FFTtools::getHilbertEnvelope(grH);
        
        //Get peaks of Hilbert Envelopes
        double peakV = TMath::MaxElement(hilbertV->GetN(), hilbertV->GetY());
        double peakH = TMath::MaxElement(hilbertH->GetN(), hilbertH->GetY());
        
        //Identify primary channel that has a larger Hilbert peak and add conditional statement to set primary and secondary waveforms
        TGraph *primaryHilbert;
        TGraph *primaryGr;
        TGraph *secondaryHilbert;
        TGraph *secondaryGr;
        double *peakPrimary;
        double *peakTimePrimary;
        double *peakSecondary;
        double *peakTimeSecondary;
        if (peakV >= peakH) {
            primaryHilbert = hilbertV;
            primaryGr = grV;
            secondaryHilbert = hilbertH;
            secondaryGr = grH;
            peakPrimary = &hilbertPeakOut[channelPair];
            peakTimePrimary = &peakTimeOut[channelPair];
            peakSecondary = &hilbertPeakOut[channelPair+8];
            peakTimeSecondary = &peakTimeOut[channelPair+8];          
            
        }
        else {
            primaryHilbert = hilbertH;
            primaryGr = grH;
            secondaryHilbert = hilbertV; 
            secondaryGr = grV;
            
            peakPrimary = &hilbertPeakOut[channelPair+8];
            peakTimePrimary = &peakTimeOut[channelPair+8];
            peakSecondary = &hilbertPeakOut[channelPair];
            peakTimeSecondary = &peakTimeOut[channelPair];               
        }
        
        //Find peak of primary channel and get its time
        int peakPrimaryIndex = TMath::LocMax(primaryHilbert->GetN(), primaryHilbert->GetY());
        *peakPrimary = primaryHilbert->GetY()[peakPrimaryIndex];
        int primaryPolarity;
        if (TMath::LocMax(primaryGr->GetN(), primaryGr->GetY()) < TMath::LocMin(primaryGr->GetN(), primaryGr->GetY())) {        
            primaryPolarity = 1;
        }
        else {
            primaryPolarity = -1;
        }
        if (findPolarity) {*peakPrimary *= primaryPolarity;} 
        *peakTimePrimary = primaryHilbert->GetX()[peakPrimaryIndex];
        
        //Use peak location of primary channel to crop secondary channel to time window tolerance and find secondary peak
        //Crop secondary Hilbert to be centered around the primary peak with the tolerance window
        
        //Adding condition where crop region is outside of the domain of the waveform.
        if (*peakTimePrimary+toleranceNs < secondaryHilbert->GetX()[0]) {
            secondaryHilbert = FFTtools::cropWave(secondaryHilbert, secondaryHilbert->GetX()[0], secondaryHilbert->GetX()[0]+2*toleranceNs);
        }
        else if (*peakTimePrimary-toleranceNs > secondaryHilbert->GetX()[secondaryHilbert->GetN()-1]) {
            secondaryHilbert = FFTtools::cropWave(secondaryHilbert, secondaryHilbert->GetX()[secondaryHilbert->GetN()-1]-2*toleranceNs, secondaryHilbert->GetX()[secondaryHilbert->GetN()-1]);
        }
        else {
            secondaryHilbert = FFTtools::cropWave(secondaryHilbert, *peakTimePrimary-toleranceNs, *peakTimePrimary+toleranceNs);
        }
        int peakSecondaryIndex = TMath::LocMax(secondaryHilbert->GetN(), secondaryHilbert->GetY());
        *peakSecondary = secondaryHilbert->GetY()[peakSecondaryIndex];
        int secondaryPolarity;
        if (TMath::LocMax(secondaryGr->GetN(), secondaryGr->GetY()) < TMath::LocMin(secondaryGr->GetN(), secondaryGr->GetY())) {
            secondaryPolarity = 1;
        }
        else {
            secondaryPolarity = -1;
        }
        if (findPolarity) {*peakSecondary *= secondaryPolarity;}         
        if (findPolarity) {*peakSecondary *= signum(secondaryGr->GetY()[peakSecondaryIndex]);}
        *peakTimeSecondary = secondaryHilbert->GetX()[peakSecondaryIndex];        
        
    }
    
}

Vector calculatePolVector(double psi, double theta, double phi) {
    Vector pol_Vector{-cos(psi)*cos(theta)*cos(phi)+sin(psi)*sin(phi),
                      -cos(psi)*cos(theta)*sin(phi)-sin(psi)*cos(phi),
                       cos(psi)*sin(theta)};
    return pol_Vector;
}

Vector calculateLaunchVector(double theta, double phi) {
    Vector launch_vector{cos(phi)*sin(theta),
                         sin(phi)*sin(theta),
                         cos(theta)};
    return launch_vector;
}

Vector calculateThetaHat(double theta, double phi) {
    Vector theta_hat{cos(theta)*cos(phi), 
                     cos(theta)*sin(phi), 
                     -sin(theta)};
    return theta_hat;
}

Vector calculatePhiHat(double phi) {
    Vector phi_hat{-sin(phi), 
                   cos(phi), 
                   0};
    return phi_hat;
}

void calculatePsi(double* hilbertPeaks, double* psiOut, bool findPolarity=false) {
    // double psiOut[8];
    double tanY;
    double tanX;
    
    for (int channelPair=0; channelPair<8; channelPair++) {
        tanY = hilbertPeaks[channelPair+8];
        tanX = hilbertPeaks[channelPair];        
        if (findPolarity) {
            psiOut[channelPair] = TMath::ATan2(-tanY,-abs(tanX))*180/PI;
        }
        else {
            psiOut[channelPair] = TMath::ATan2(tanY,tanX)*180/PI;
        }
    }
}

double calculateTruePsi(Event *event, Report *report, double launch_theta, double launch_phi) {
    
    Vector launch_vector = calculateLaunchVector(launch_theta, launch_phi);
    
    
    Vector Pol_vector = report->GetPolarization(event->Nu_Interaction[0].nnu, launch_vector);
    
    Vector theta_hat = calculateThetaHat(launch_theta, launch_phi);
    Vector phi_hat = calculatePhiHat(launch_phi);
    
    double psi = atan2(-Pol_vector.Dot(phi_hat),-Pol_vector.Dot(theta_hat))*180/PI;
    
    return psi;
    
}

void calculatePsiSolutions(double* psiIn, double* psi2, double* psi3, double* psi4) {
    for (int channelPair=0; channelPair<8; channelPair++) {
        psi2[channelPair] = -psiIn[channelPair];
        psi3[channelPair] = 180 - abs(psiIn[channelPair]);
        psi4[channelPair] = -1*psi3[channelPair];
    }
}


double calculateCherenkovAngle(double nofz) {
        return acos(1/nofz);
}

double calculateViewingAngle(Vector nnu, Vector launch_vector) {
    return nnu.Angle(launch_vector);
}

void calculateNuTrajectory(double psi, double launch_theta, double launch_phi, double viewing_angle, double &theta_nutraject, double &phi_nutraject) {
    
    //Calculate polarization and launch vectors from the input angles    
    Vector launch_vector = calculateLaunchVector(launch_theta, launch_phi);    
    Vector polarization_vector = calculatePolVector(psi, launch_theta, launch_phi);       
    Vector nutraject_vector = cos(viewing_angle)*launch_vector - sin(viewing_angle)*polarization_vector;    
    theta_nutraject = nutraject_vector.Theta()*180/PI;
    phi_nutraject = nutraject_vector.Phi()*180/PI;     
}

void testTrajectoryCalculation(Event *event, Report *report, double launch_theta, double launch_phi, double nofz) {
    
    Vector launch_vector = calculateLaunchVector(launch_theta, launch_phi);    
    Vector Pol_vector = report->GetPolarization(event->Nu_Interaction[0].nnu, launch_vector);  
    double changle = event->Nu_Interaction[0].changle;
    double angleNnuLaunch = event->Nu_Interaction[0].nnu.Angle(launch_vector);
    Vector nutraject_vector = cos(angleNnuLaunch)*launch_vector - sin(angleNnuLaunch)*Pol_vector; 
    double theta_nutraject = nutraject_vector.Theta()*180/PI;
    double phi_nutraject = nutraject_vector.Phi()*180/PI;  
}



double getStationDepth(AraGeomTool *geomTool, int stationId) {
    double station_z = 0;
    for (int i=0; i<16; i++) {
        station_z += geomTool->getStationInfo(stationId)->getAntennaInfo(i)->antLocation[2]/16;
    }
    return station_z;
    
}

double getNofzAtVertex(AraGeomTool *geomTool, UsefulAtriStationEvent *usefulAtri, IceModel *icemodel, double vertexRadius, double vertexTheta) {
    //Calculate station center wrt surface
    double station_z = 0;
    station_z = getStationDepth(geomTool, usefulAtri->stationId);
    cout << "station_z = " << station_z << endl;
    double vertex_z = station_z + vertexRadius*sin(vertexTheta);
    cout << "vertex_z = " << vertex_z << endl;
    //Create position object for station location
    double nofz_atVertex;
    if (vertex_z < 0) {
        nofz_atVertex = icemodel->GetN(vertex_z);
    }
    else {  //Case for event resolving in the air.
        // cout << "Vertex resolves in the air." << endl;
        // nofz_atVertex = 1;
        nofz_atVertex = icemodel->GetN(-vertex_z);
    }
    cout << "nofz_atVertex = " << nofz_atVertex << endl;
    return nofz_atVertex;
}

void setPedestalFile(AraEventCalibrator *cal, int station, char* runNum) {
    const char* pedestalLocation = "/fs/project/PAS0654/ARA_DATA/pedestalFiles/"; //Hardcoding pedestal path
    char* pedestalFile;
    fstream fileStream;
    char test1[400];
    char test2[400];
    char test3[400];
    
    // Find the position of the underscore character
    const char* underscorePos = strchr(runNum, '_');
    
    // Create a buffer to hold the extracted part
    char extractedNum[256]; // Adjust size as needed

    // If an underscore is found, copy the part before the underscore
    if (underscorePos) {
        // Calculate the length of the substring before the underscore
        size_t length = underscorePos - runNum;
        
        // Copy the substring to the buffer and null-terminate it
        strncpy(extractedNum, runNum, length);
        extractedNum[length] = '\0'; // Null-terminate the string
    } else {
        // If no underscore is found, copy the entire original string
        strncpy(extractedNum, runNum, sizeof(extractedNum) - 1);
        extractedNum[sizeof(extractedNum) - 1] = '\0'; // Null-terminate the string
    }    
    
    
    sprintf(test1, "%s/A%d/ped_full_values_A%d_run0%s.dat", pedestalLocation, station, station, extractedNum);
    sprintf(test2, "%s/A%d/ped_full_values_A%d_run00%s.dat", pedestalLocation, station, station, extractedNum);
    sprintf(test3, "%s/A%d/ped_full_values_A%d_R%s.dat", pedestalLocation, station, station, extractedNum);
    fileStream.open(test1);
    cout << "Checking " << test1 << endl;
    if (fileStream) {
         pedestalFile = test1;
        // sprintf(pedestalFile, test1);
    } 
    else {
        fileStream.open(test2);
        cout << "Checking " << test2 << endl;
        if (fileStream) {
             pedestalFile = test2;
            // sprintf(pedestalFile, test2);
        }         
        else {
            fileStream.open(test3);
            cout << "Checking " << test3 << endl;
            if (fileStream) {
                 pedestalFile = test3;
                // sprintf(pedestalFile, test3);
            }    
            else {
                cout << "Pedestal file not found!" << endl;
                // exit(0);
            }
        }
    }     
    
    // AraEventCalibrator *cal = AraEventCalibrator::Instance();
    cout << "Using pedestal file: " << pedestalFile << endl;
    if (station == 100){
        cal->setAtriPedFile(pedestalFile, 1);
    }
    else {
        cal->setAtriPedFile(pedestalFile, station);    
    }
}

//Functions to read in the pulser depth file and interpolate it so we can provide a unixtime and get the pulser depth.
double getPulserDepth(double queryTime) {
    std::string filename = "/users/PAS0654/jflaherty13/thesisWork/tools/spicePulser_2018Depth_unixTime.txt";
    std::vector<double> times;
    std::vector<double> positions;

    // Read the file and populate times and positions vectors
    std::ifstream file(filename);
    if (!file.is_open()) {
        throw std::runtime_error("Could not open file");
    }

    std::string line;
    while (std::getline(file, line)) {
        std::stringstream ss(line);
        std::string timeStr, positionStr;
        
        std::getline(ss, timeStr, ',');
        std::getline(ss, positionStr, ',');

        try {
            double time = std::stod(timeStr);
            double position = std::stod(positionStr);
            times.push_back(time);
            positions.push_back(position);
        } catch (const std::invalid_argument&) {
            // Skip lines that cannot be parsed
            continue;
        }
    }
    file.close();

    // Perform linear interpolation
    if (queryTime < times.front() || queryTime > times.back()) {
        throw std::out_of_range("Query time is out of range");
    }

    for (size_t i = 0; i < times.size() - 1; ++i) {
        if (queryTime >= times[i] && queryTime <= times[i + 1]) {
            double t1 = times[i];
            double t2 = times[i + 1];
            double p1 = positions[i];
            double p2 = positions[i + 1];

            // Linear interpolation formula
            return p1 + (p2 - p1) * (queryTime - t1) / (t2 - t1);
        }
    }

    throw std::runtime_error("Interpolation error: data not found");
    
}

void getPulserVertex(Detector *detector, Settings *settings1, AraGeomTool *geomTool, double pulserDepth, double& vertexTheta, double& vertexPhi, double& vertexRadius, bool debugMode=false) {
    
    //Hard code latitude and longitude of SpiceCore
    double sourceLatitude = -89.97953;
    double sourceLongitude = -100.78595;
    
    //Calculate array coordinates of SpiceCore
    Int_t stationId = settings1->DETECTOR_STATION_ARAROOT;
    double sourceEasting = AraGeomTool::getArrayEastingFromLatLong(sourceLatitude, sourceLongitude);
    double sourceNorthing = AraGeomTool::getArrayNorthingFromLatLong(sourceLatitude, sourceLongitude);    
    
    //Construct vector of source in array coordinates
    TVector3 sourceArrayVector;
    sourceArrayVector[0] = sourceEasting;
    sourceArrayVector[1] = sourceNorthing;
    sourceArrayVector[2] = -1*pulserDepth - getStationDepth(geomTool, settings1->DETECTOR_STATION);
    
    //Convert array vector into station-centric coordinates
    TVector3 sourceStationVector = geomTool->convertArrayToStationCoords(stationId, sourceArrayVector); 
    
    //Calculate vertex in station-centric coordinates
    vertexRadius = sqrt(pow(sourceStationVector[0],2) + pow(sourceStationVector[1],2) + pow(sourceStationVector[2],2));
    vertexPhi = (360+(atan2(sourceStationVector[1], sourceStationVector[0]))*180/PI)*PI/180;
    vertexTheta = acos((sourceStationVector[2])/vertexRadius);      
    
    if (debugMode) {
        cout << "Pulser latitude and longitude: (" << sourceLatitude << ", " << sourceLongitude << ")" << endl;
        cout << "Pulser array Easting and Northing: (" << sourceEasting << ", " << sourceNorthing << ")" << endl;
        cout << "Pulser station-centric Easting and Northing: (" << sourceStationVector[0] << ", " << sourceStationVector[1] << ")" << endl;
    }
}

//Write function that takes inputs of calculated vertex radius, theta, and phi to find the nearest bins for the ray tracing tabel.  Theta bins are centered around the 0.5 degree mark, so we need to take that into consideration.
void getNearestVertexBin(double radiusIn, double thetaIn, double phiIn, double angleResolution, double radiusResolution, double& radiusOut, double& thetaOut, double& phiOut) {
    radiusOut = std::round(radiusIn/radiusResolution)*radiusResolution;
    thetaOut = std::round((thetaIn-angleResolution/2)/angleResolution)*angleResolution+angleResolution/2;
    phiOut = std::round((phiIn-angleResolution/2)/angleResolution)*angleResolution+angleResolution/2;
    
}

void applyBandpassBin(double &vm_real, double &vm_img, double freqBin, double freqMin, double freqMax, int order=8) {
        double weight = 1;  // Setting initial weight to one, then applying bandpass.  Weight is then multiplied by signal in this bin.
        // int order = 8;
        weight /= sqrt(1 + TMath::Power(freqMin/freqBin, 4*order));
        weight /= sqrt(1 + TMath::Power(freqBin/freqMax, 4*order));
        vm_real *= weight;
        vm_img *= weight;
    
}

int createFFTsize(double dT_forfft, Settings *settings1) {
    // double dT_forfft = time[1] - time[0];
    int Ntmp = settings1->TIMESTEP *1.e9 / dT_forfft;
    int Nnew = 1;         
    while (Ntmp > 1)
    {
        Ntmp = Ntmp / 2;
        Nnew = Nnew *2;             
    }
    Nnew = Nnew * settings1->NFOUR / 2;  
    return Nnew;
}


TGraph *butterworthFilter(TGraph *grWave, double_t minFreq, double_t maxFreq, int order=8) {
    double *oldY = grWave->GetY();
    double *oldX = grWave->GetX();
    double deltaT=oldX[1]-oldX[0];
    int length=grWave->GetN();
    FFTWComplex *theFFT=FFTtools::doFFT(length,oldY);

    int newLength=(length/2)+1;

    double deltaF=1/(deltaT*length); //Hz
    deltaF*=1e3; //MHz

    double tempF=0;
    for(int i=0;i<newLength;i++) {
        double weight = 1;  // Setting initial weight to one, then applying bandpass.  Weight is then multiplied by signal in this bin.
        weight /= sqrt(1 + TMath::Power(minFreq/tempF, 4*order));
        weight /= sqrt(1 + TMath::Power(tempF/maxFreq, 4*order));  
        theFFT[i].re *= weight;
        theFFT[i].im *= weight;
        tempF+=deltaF;
    }

    double *filteredVals = FFTtools::doInvFFT(length,theFFT);


    TGraph *grFiltered = new TGraph(length,oldX,filteredVals);
    delete [] theFFT;
    delete [] filteredVals;
    return grFiltered;        
}
