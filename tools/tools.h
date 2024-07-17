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
        cout << "Using station 2.  Excluding channels 1 and 3 in addition to the masking." << endl;
        excludedChannels.push_back(1);  //Removing extra A2 channels as the D-pulse gets clipped.
        excludedChannels.push_back(3); 
        excludedChannels.push_back(5); 
        excludedChannels.push_back(6); 
        //Testing keeping flat channels in CSW histogram
        // excludedChannels.push_back(0);
        // excludedChannels.push_back(2);
    }    
    if (settings1->DETECTOR_STATION==4){
        cout << "Using station 4.  Excluding channels 0 in addition to the masking." << endl;
        excludedChannels.push_back(0);  //Testing removing a channel for A4 as only the R pulse makes it into the waveform.
        // excludedChannels.push_back(4); 
        // excludedChannels.push_back(8);
    }
    for (int i=0; i<16; i++){
        // cout << "detector->GetTrigMasking(i) = " << detector->GetTrigMasking(i) << endl;
        if (not detector->GetTrigMasking(i)){
            cout << "Excluding channel " << i << endl;
            excludedChannels.push_back(i);
        }
    }       
}

double calculateNoiseRMS(TGraph *gr, int sampleNumber=100) {
    int totalSampleNumber = abs(sampleNumber);
    int waveformLength = gr->GetN();
    double voltageSubset[totalSampleNumber];
    if (sampleNumber<0){
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
        arrivalTimeAvg += arrivalTimes[i]/16;
    }    
    //Import waveforms and Apply time delays
    TGraph *gr;
    for (int i=0; i<16; i++) {
        bool checkExcluded = std::find(cswExcludedChannelPairs.begin(), cswExcludedChannelPairs.end(), i%8) != cswExcludedChannelPairs.end();
        if (not checkExcluded) {        
            //Import waveform as Tgraph.
            gr = usefulAtriEvPtr->getGraphFromRFChan(i);
            if (debugMode){
                cout << "i = " << i << endl;
                cout << "gr->GetN() Before crop = " << gr->GetN() << endl;
                cout << "gr->GetX()[0] = " << gr->GetX()[0] << endl;
                cout << "gr->GetX()[gr->GetN()-1] = " << gr->GetX()[gr->GetN()-1] << endl;
                cout << "cropTimes[i] = " << cropTimes[i]<< endl; 
                cout << "arrivalTimes[i] = " << arrivalTimes[i] << endl;
                cout << "arrivalTimeAvg = " << arrivalTimeAvg << endl;
                
            }
            if (cropTimes[i] < gr->GetX()[0]+sampleNs) {
                cswExcludedChannelPairs.push_back(i);
                continue;
            }
            //Crop to D-pulse only.
            gr = FFTtools::cropWave(gr, gr->GetX()[0], cropTimes[i]);
            if (debugMode) {
                cout << "gr->GetN() After crop = " << gr->GetN() << endl; 
            }

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
            if (debugMode){
                cout << "i = " << i << endl;
                cout << "mins[i] = " << mins[i] << endl;
                cout << "maxes[i] = " << maxes[i] << endl;      
            }
            if (tMin < mins[i] and mins[i] < tMax) {tMin = mins[i];}
            if (tMax > maxes[i] and maxes[i] > tMin) {tMax = maxes[i];}
            if (debugMode) {
                cout << "tMin = " << tMin << endl;
                cout << "tMax = " << tMax << endl; 
                cout << "*****************************" << endl;
            }
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
    
    //Initialize time and voltage arrays for the coherent sum.
    const int nT = int((tMax-tMin)/dt);
    double timeCsw[nT];
    double voltageCswVpol[nT];    
    double voltageCswHpol[nT];      
    
    //Loop over waveforms and crop to the global tMin and tMax
    for (int i=0; i<16; i++) {
        bool checkExcluded = std::find(cswExcludedChannelPairs.begin(), cswExcludedChannelPairs.end(), i%8) != cswExcludedChannelPairs.end();
        if (not checkExcluded) {        
            mapWaveforms[i] = FFTtools::cropWave(mapWaveforms[i], tMin, tMax);
            
            if (mapWaveforms[i]->GetN() == 0) {
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
    
    if (debugMode) {
        cout << "Initial gr = " << gr->GetN() << endl;
        cout << "gr->GetX()[0] = " << gr->GetX()[0] << endl;
        cout << "gr->GetX()[gr->GetN()-1] = " << gr->GetX()[gr->GetN()-1] << endl;
    }
    
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
    if (debugMode) {
        cout << "tNoiseMin = " << tNoiseMin << endl;
        cout << "tNoiseMax = " << tNoiseMax << endl;
        cout << "tSignalMin = " << tSignalMin << endl;
        cout << "tSignalMax = " << tSignalMax << endl;
        cout << "Initial grSignal = " << grSignal->GetN() << endl;
        cout << "Initial grNoise = " << grNoise->GetN() << endl;            
    }

    
    //Pad waveforms to match the frequency binning used in the deconvolution
    grSignal = resizeForFFT(grSignal, settings1->NFOUR);
    grNoise = resizeForFFT(grNoise, settings1->NFOUR);
    // grSignal = FFTtools::padWaveToLength(grSignal, settings1->NFOUR/2);
    // grNoise = FFTtools::padWaveToLength(grNoise, settings1->NFOUR/2);
    if (debugMode) {
        cout << "Cropped grSignal = " << grSignal->GetN() << endl;
        cout << "Cropped grNoise = " << grNoise->GetN() << endl;
    }
        
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