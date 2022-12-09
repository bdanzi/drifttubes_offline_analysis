//contenitore di forme d'onda
#ifndef DATAUTIL_H
#define DATAUTIL_H

#include <vector>
#include <map>

#include "TH1F.h"
#include "TH2F.h"
#include "TMath.h"
#include <ostream>
#include <fstream>
static int skipFstBin[5] = {5,5,5,5,5}; //525 for El Cal   //100 dati proto
static int skipLstBin[5] = {10,10,10,10,10};  //475 for El Cal   //50	dati proto
// if (isNov2021TestBeam){
// static int ChannelDiameter[13] = {-1,-1,-1,-1,10,15,20,20,25,25,20,25,40}; // Old test beam Nov 2022
// static float ChannelCellSize[13] = {-1.0,-1.0,-1.0,-1.0,1.0,1.0,1.0,1.0,1.0,1.0,2.0,2.0,2.0}; // Old test beam Nov 2022
// }
// else if(!isNov2021TestBeam){
static int ChannelDiameter[16] = {20,20,15,15,25,20,20,15,20,25,25,15,-1,-1,-1,-1}; // Old test beam Nov 2022
static float ChannelCellSize[16] = {1.5,1.0,1.0,1.0,1.5,1.0,1.0,1.5,1.0,1.0,1.0,1.5,-1.0,-1.0,-1.0,-1.0}; // Old test beam Nov 2022
// }
//static int fbin_rms[5]={30,30,30,30,30};
//static int fbin_bsl[5]={30,30,30,30,30};
static float binsize = 0.0;
static int fbin_rms = 0.0; //30 ns
static int fbin_bsl = 0.0; //30 ns
// fbin_bsl = (int) (30 ns)/binsize
// skipFst = (int) fixed to be changable

static int isl = 0;
//static float invfbin=1.0/((float)fbin_rms);
static float invfbin;
static float _gsample = 0.0;

struct wave {
  std::vector<float> Y; //vettore dei dati.
  //variabili di appoggio
  float bsln; //baseline
  float rms;
  float integ;
  float max;
  float integInR;
  float maxInR;
  float min;
  float minInR;
  int   maxPos;
  int   maxInRPos;
  int   minPos;
  int   minInRPos;
  float sumX;
  float sumY;
  float sumXY;
  float sumX2;
  int   nPtReg;
  float regA;
  float regB;
  std::vector<float> deriv;
  std::vector<float> sderiv;
  float charge_integInR;
  float charge_integ;
  
wave() : bsln(0.0), rms(0.0), integ(0.0), max(-9999), integInR(0.0), maxInR(-9999), min(9999), minInR(9999)
    , maxPos(-1), maxInRPos(-1), minPos(-1), minInRPos(-1)
    , sumX(0.0), sumY(0.0), sumXY(0.0), sumX2(0.0), nPtReg(0), regA(0.0), regB(0.0), charge_integInR(0.0),  charge_integ(0.0) {
  Y.clear();
  deriv.clear();
  sderiv.clear();
}/////inizializzazione
  
  void clear() {
    Y.clear();
    bsln=0.0;
    rms=0.0;
    integ=0.0;
    charge_integInR = 0.0;
    charge_integ = 0.0;
    max=-9999;
    integInR=0.0;
    maxInR=-9999;
    min=9999;
    minInR=9999;
    maxPos=-1;
    maxInRPos=-1;
    minPos=-1;
    minInRPos=-1;
    sumX=0.0;
    sumY=0.0;
    sumXY=0.0;
    sumX2=0.0;
    regA=0.0;
    regB=0.0;
    deriv.clear();
    sderiv.clear();
  }
  
  int nPt() { return Y.size(); }
  int nPtInR() { return (Y.size()-(skipFstBin[isl]+skipLstBin[isl])); }
  float nInteg() {return (integ-bsln*((float)nPt()));}//funzione float che deve ritornare l'integrale  a cui togliamo la baseline. ////
  float nIntegInR() { return (integInR); }//funzione float che deve ritornare l'integrale  a cui togliamo la baseline. ////
  float nMax() { return (max-bsln); }
  float nMaxInR() { return (maxInR-bsln); }
  float nMin() {return (min-bsln);}
  float nMinInR(){return (minInR-bsln);}
  
  float n1Integ() { return (integ-norm1(0.0,(float)(nPt()-1))); }//funzione float che deve ritornare l'integrale  a cui togliamo la baseline. ////
  float n1IntegInR() { return (integInR-norm1((float)skipFstBin[isl],(float)(nPt()-1-skipLstBin[isl]))); }//funzione float che deve ritornare l'integrale  a cui togliamo la baseline. ////
  float n1Max() { return (max-(RegA()+RegB()*((float)maxPos))); }
  float n1MaxInR() { return (maxInR-(RegA()+RegB()*((float)maxInRPos))); }
  float n1Min() {return (min-(RegA()+RegB()*((float)minPos)));}
  float n1MinInR(){return (minInR-(RegA()+RegB()*((float)minInRPos)));}
  float norm1(float x1, float x2) { return (RegA()*(x2-x1)+0.5*RegB()*(x2*x2-x1*x1)); }
  
  float n1At(int i) { return (Y[i]-(RegA()+RegB()*((float)i))); }
  float nAt(int i) { return (Y[i]-bsln); }
  
  float nnAt(int i, float cut=0.05) { //0.05 proto_cosmic
    if (nMaxInR()>cut) { return nAt(i); }
    else { return n1At(i); }
  }
  float nnInteg(float cut=0.05) {//0.05 proto_cosmic
    if (nMaxInR()>cut) { return nInteg(); }
    else { return n1Integ(); }
  }
  float nnIntegInR(float cut=0.05) {//0.05 proto_cosmic
    if (nMaxInR()>cut) { return nIntegInR(); }
    else { return n1IntegInR(); }
  }
  
  float RegA() {
    if (regA==0.0) {
      float N=nPtReg;//skipFstBin+skipLstBin;
      regA=(sumY*sumX2-sumX*sumXY)/(N*sumX2-sumX*sumX);
    }
    return regA;
  }
  float RegB() {
    if (regB==0.0) {
      float N=nPtReg;//skipFstBin+skipLstBin;
      regB=(N*sumXY-sumX*sumY)/(N*sumX2-sumX*sumX);
    }
    return regB;
  }
  //metodi per riempire la waveform
  //void fillWave(int nPt, __uint16_t *arr) {
  void fillWave(std::vector<__uint16_t> &arr, int maxDim=-1, float _gsample= 1.0, float bslnTimeInterval = 25.0 ) {
    float tmpval;
    float binsize = (float) (1/(_gsample)); // units of nanosecond
    //std::cout << "Binsize first function: "<< binsize << std::endl; 
    int fbin_rms = (int) ((bslnTimeInterval)/binsize); // integer
    //std::cout << "fbin_rms" << fbin_rms << std::endl;
    int fbin_bsl = (int) ((bslnTimeInterval)/binsize);
    float integ_partial = 0.0;
    float invfbin =  1.0/((float) fbin_rms);
    //std::cout << "fbin_rms first function: "<< fbin_rms << std::endl; 
    //std::cout << "Here is the sampling rate "<<_gsample; 
    int nPt = arr.size();   //dimensione del vettore arr.
    if (maxDim!=-1 && maxDim<nPt) { nPt=maxDim; }
    clear();               //richiama la funzione sopra definita. 
    //riempimento della variabili di appoggio.
    for (int ipt=skipFstBin[isl]; ipt<fbin_rms+skipFstBin[isl]; ++ipt) { // calcolo integ per rms
      tmpval=((1.0/65536.0)*arr[ipt]-0.5); //Volt from DRS
      tmpval*=-1;
      integ_partial+=(tmpval); //Volt from DRS
    }
    bsln = (float) integ_partial/fbin_rms; // Volt/ number of bins
    for (int ipt=skipFstBin[isl]; ipt<fbin_rms+skipFstBin[isl]; ++ipt) {
      tmpval=((1.0/65536.0)*arr[ipt]-0.5);//Volt from DRS
      tmpval*=-1;//Volt from DRS
      rms+=(tmpval - bsln)*(tmpval - bsln);
    }
    rms = sqrt((rms)/(fbin_rms)); // Volt / number of bins
    //std::cout << "Rms first function: " << rms <<"\n";
    for (int ipt=0; ipt<nPt; ++ipt) {
      tmpval=((1.0/65536.0)*arr[ipt]-0.5);
      tmpval*=-1;
      //std::cout << "i-th point "<<ipt<< "tmpval: " << tmpval <<"\n";
      Y.push_back(tmpval-bsln);
      integ+=(Y.back());   // Volt        //l'integrale è in sostanza una somma, integ � inizializzata a zero e poi viene incrementata.
      //Rms somma Y nei primi 100 bin
      // if (ipt==fbin_rms[isl]-1) {
      //std::cout << "i-th point "<<ipt<<"Bsl: " << bsln <<"\n";
      //if (ipt==fbin_rms[isl]-1/*100*/) {
      if (Y.back()>max) { max=Y.back(); maxPos=ipt;} //std::cout << "i-th point "<<ipt<<"max: " << max <<"\n";}
      if (Y.back()<min) { min=Y.back(); minPos=ipt; }
      if (ipt>=skipFstBin[isl] && ipt<(nPt-skipLstBin[isl])) {
	integInR+=(Y.back());
	if (Y.back()>maxInR) { maxInR=Y.back(); maxInRPos=ipt; }
	if (Y.back()<minInR) { minInR=Y.back(); minInRPos=ipt; }
      } else {
	sumX+=ipt;
	sumX2+=ipt*ipt;
	sumY+=Y.back();
	sumXY+=((float)ipt)*Y.back();
	++nPtReg;
      }
    }//fine ciclo sui punti.
    integ = (integ) * binsize; //V * ns
    //std::cout << "Bsl: " << bsln <<"\n";
    integInR = (integInR) * binsize; // Volt * ns
    charge_integ = (float)(integ/(50.0) * 1.0e+3); //integral in pC 
    charge_integInR = (float)(integInR/(50.0) * 1.0e+3);
    for(int i=0;i<nPt;++i){
      // if you arrive at the end 
      deriv.push_back((Y[(i+1)>(nPt-1)?(nPt-1):(i+1)]-Y[(i-1)<0?0:(i-1)])/(2*binsize));
    }
    for(int i=0;i<nPt;++i){
      sderiv.push_back((deriv[(i+1)>(nPt-1)?(nPt-1):(i+1)]-deriv[(i-1)<0?0:(i-1)])/(2*binsize));
    }
    
  }//fine funzione riempimento.
  
  //void fillWave(int nPt, float *arr) {
  void fillWave(std::vector<float> &arr, int maxDim=-1, float _gsample= 1.0, float bslnTimeInterval = 25.0 ) { // to fill Wave_signal_1
    int nPt=arr.size();
    float tmpval;
    //std::cout << "GigaSample second function: "<< _gsample << std::endl;
    float binsize = (float) (1/(_gsample)); // units of nanosecond
    int fbin_rms = (int) ((bslnTimeInterval)/binsize); // trigger
    int fbin_bsl = (int) ((bslnTimeInterval)/binsize);
    float integ_partial = 0.0;
    float charge_integ = 0.0;
    float charge_integInR = 0.0;
    rms = 0;
    float invfbin =  1.0/((float) fbin_rms);
    //std::cout << "fbin_rms second function: "<< fbin_rms << std::endl; 
    if (maxDim!=-1 && maxDim<nPt) { nPt=maxDim; }
    clear();
    // Rms computation
    for (int ipt=skipFstBin[isl]; ipt<fbin_rms+skipFstBin[isl]; ++ipt) { // calcolo integ per rms
      tmpval=((1.0/65536.0)*arr[ipt]-0.5);
      tmpval*=-1;
      integ_partial+=(tmpval); 
    }
    bsln = (float) integ_partial/fbin_rms; // V * ns / #bins
    for (int ipt=skipFstBin[isl]; ipt<fbin_rms+skipFstBin[isl]; ++ipt) {
      tmpval=((1.0/65536.0)*arr[ipt]-0.5);
      tmpval*=-1;
      rms+=(tmpval - bsln)*(tmpval - bsln);
    }
    rms = sqrt((rms)/(fbin_rms));
    //std::cout << "Rms second function: " << rms <<"\n";
    //definizione delle variabili d'appoggio
    for (int ipt=0; ipt<nPt; ++ipt) {
      Y.push_back(arr[ipt]);        
      integ+=(Y.back());
      if (Y.back()>max) { max=Y.back(); maxPos=ipt; }
      if (Y.back()<min) { min=Y.back(); minPos=ipt; }
      if (ipt>=skipFstBin[isl] && ipt<(nPt-skipLstBin[isl])) {
	integInR+=Y.back();
	if (Y.back()>maxInR) { maxInR=Y.back(); maxInRPos=ipt; }
	if(Y.back()<minInR)  { minInR=Y.back(); minInRPos=ipt; }
      } else {
	sumX+=ipt;
	sumX2+=ipt*ipt;
	sumY+=Y.back();
	sumXY+=((float)ipt)*Y.back();
	      ++nPtReg;
      }
      //if (print) std::cout<<"ipt "<<ipt<<" pnt "<<Y.back()<<" integ "<<integ<<std::endl;
    }
    integ = (integ)*binsize;
    integInR = (integInR) * binsize; // Volt * ns
    charge_integ = (float)(integ/(50.0) * 1.0e+3); //integral in pC 
    charge_integInR = (float)(integInR/(50.0) * 1.0e+3);
    for(int i=skipFstBin[isl];i<nPt-skipLstBin[isl];++i){
      deriv.push_back((Y[(i+1)>(nPt-1)?(nPt-1):(i+1)]-Y[(i-1)<0?0:(i-1)])/(2*binsize));
    }
    for(int i=skipFstBin[isl];i<nPt-skipLstBin[isl];++i){
      sderiv.push_back((deriv[(i+1)>(nPt-1)?(nPt-1):(i+1)]-deriv[(i-1)<0?0:(i-1)])/(2*binsize));
    }
  }
  
  /* void addPnt(float &val, bool InRng=true) {
     Y.push_back(val);
     int ipt=nPt()-1;
     integ+=val;
     
     if (ipt<fbin_rms) {
     rms+=val*val;
     }
     if (ipt==fbin_rms-1) {
     bsln=integ*invfbin;
     rms=sqrt(rms*invfbin-bsln*bsln);
     }
     if (val>max) { max=val; maxPos=ipt; }
     if (val<min) { min=val; minPos=ipt; }
     if (ipt>=skipFstBin[isl] && InRng) {
     integInR+=Y.back();
     if (Y.back()>maxInR) { maxInR=Y.back(); maxInRPos=ipt; }
     if (Y.back()<minInR) { minInR=Y.back(); minInRPos=ipt; }
     } else {
     sumX+=ipt;
     sumX2+=ipt*ipt;
     sumY+=Y.back();
     sumXY+=((float)ipt)*Y.back();
     ++nPtReg;
     }
     } */
}; //fine della struct wave


typedef std::map<int,wave> WvCont;                    //typedef assegna un altro nome al tipo specificato.
typedef std::map<std::pair<int,int>,wave> diffWvCont; //serve per fare la differenza



//struttura per costruire istogrammi
struct hstPerCh { //istogrammi per tutti i canali dell'oscilloscopio.
  //isto che verranno riempiti con wave filtrate
  TH1F *hTimeDifference;
  TH1F *hTimeDifference_clust;
  TH1F *hHPeaks;
  TH2F *hHNPeaks;
  TH2F *hNPeakFPeak;
  TH2F *hNClusterFCluster;
  TH1F *hTPeaks;
  TH1F *hTFstPeaks;
  TH1F *hNPeaks_1;
  TH1F *hBsl;
  TH1F *hIntegN;
  TH1F *hNeventSignals;
  TH1F *hNPeaks;
  TH1F *hNPeaks_clust;
  TH1F *hRms;
  TH1F *hMaxVInR;
  //Derivative study
  TH1F *hFirstDeriv;	
  TH1F *hSecDeriv;
  TH1F *hNElectrons_per_cluster;

  
  hstPerCh(int Ch=0, int isChannel_1cm=0, int isChannel_2cm=0, int isChannel_1p5cm = 0, float alpha = 0., float _gsample= 0., int isRuns_80_20 = 0, int isRuns_90_10 = 0, int isRuns_85_15= 0) {
    //gRootDir->cd();
    //TDirectory fldch(Form("H-Ch%d",Ch),Form("folder for ch%d",Ch));
    //fldch.cd();
    //senza smooth
    hTimeDifference = new TH1F (Form("hTimeDifference_ch%d",Ch),Form("Time Difference between Two Consecutive Electrons - Ch %d - Sense Wire Diameter %d um - Cell Size %0.1f cm - Track Angle %0.1f - %0.1f GSa/s - Gas Mixture 80/20 %d - 90/10 %d - 85/15 %d" ,Ch, ChannelDiameter[Ch],ChannelCellSize[Ch], alpha, _gsample, isRuns_80_20, isRuns_90_10, isRuns_85_15),40,0.,40.); 
    hTimeDifference_clust = new TH1F (Form("hTimeDifference_clust_ch%d",Ch),Form("Time Difference between Two Consecutive Clusters - Ch %d - Sense Wire Diameter %d um - Cell Size %0.1f cm - Track Angle %0.1f - %0.1f GSa/s - Gas Mixture 80/20 %d - 90/10 %d - 85/15 %d" ,Ch, ChannelDiameter[Ch],ChannelCellSize[Ch], alpha, _gsample, isRuns_80_20, isRuns_90_10, isRuns_85_15),100,0.,100.); 
    //hP0=new TH1F (Form("hP0_ch%d",Ch),Form("HP0 - Ch %d - Sense Wire Diameter %d um - Cell Size %0.1f cm - Track Angle %0.1f - %0.1f GSa/s - Gas Mixture 80/20 %d - 90/10 %d - 85/15 %d" ,Ch, ChannelDiameter[Ch],ChannelCellSize[Ch], alpha, _gsample, isRuns_80_20, isRuns_90_10, isRuns_85_15),1000,-1,1);
    //hder = new TH1F (Form("hder_ch%d",Ch),Form("Hder - Ch %d - Sense Wire Diameter %d um - Cell Size %0.1f cm - Track Angle %0.1f - %0.1f GSa/s - Gas Mixture 80/20 %d - 90/10 %d - 85/15 %d" ,Ch, ChannelDiameter[Ch],ChannelCellSize[Ch], alpha, _gsample, isRuns_80_20, isRuns_90_10, isRuns_85_15),1000,-1e+7,1e+7);
    //hchiFT = new TH1F (Form("hchiFT_ch%d",Ch),Form("HchiFT - Ch %d - Sense Wire Diameter %d um - Cell Size %0.1f cm - Track Angle %0.1f - %0.1f GSa/s - Gas Mixture 80/20 %d - 90/10 %d - 85/15 %d" ,Ch, ChannelDiameter[Ch],ChannelCellSize[Ch], alpha, _gsample, isRuns_80_20, isRuns_90_10, isRuns_85_15),1000,-10,10);
    if(isChannel_1cm){
      hNPeaks = new TH1F (Form("hNPeaks_ch%d",Ch),Form("N Electron Peaks found - Ch %d - Sense Wire Diameter %d um - Cell Size %0.1f cm - Track Angle %0.1f - %0.1f GSa/s - Gas Mixture 80/20 %d - 90/10 %d - 85/15 %d" ,Ch, ChannelDiameter[Ch],ChannelCellSize[Ch], alpha, _gsample, isRuns_80_20, isRuns_90_10, isRuns_85_15),100,0.,100);
    }
    if(isChannel_2cm || isChannel_1p5cm){
      hNPeaks = new TH1F (Form("hNPeaks_ch%d",Ch),Form("N Electron Peaks found - Ch %d - Sense Wire Diameter %d um - Cell Size %0.1f cm - Track Angle %0.1f - %0.1f GSa/s - Gas Mixture 80/20 %d - 90/10 %d - 85/15 %d" ,Ch, ChannelDiameter[Ch],ChannelCellSize[Ch], alpha, _gsample, isRuns_80_20, isRuns_90_10, isRuns_85_15),200,0.,200);
    }
    hNPeaks_clust = new TH1F (Form("hNPeaks_clust_ch%d",Ch),Form("N Cluster Peaks found - Ch %d - Sense Wire Diameter %d um - Cell Size %0.1f cm - Track Angle %0.1f - %0.1f GSa/s - Gas Mixture 80/20 %d - 90/10 %d - 85/15 %d" ,Ch, ChannelDiameter[Ch],ChannelCellSize[Ch], alpha, _gsample, isRuns_80_20, isRuns_90_10, isRuns_85_15),100,0.,100);
    hNeventSignals = new TH1F (Form("hNeventSignals_ch%d",Ch),Form("N event Signals - Ch %d - Sense Wire Diameter %d um - Cell Size %0.1f cm - Track Angle %0.1f - %0.1f GSa/s - Gas Mixture 80/20 %d - 90/10 %d - 85/15 %d" ,Ch, ChannelDiameter[Ch],ChannelCellSize[Ch], alpha, _gsample, isRuns_80_20, isRuns_90_10, isRuns_85_15),2,-0.5,1.5); 
    hNElectrons_per_cluster = new TH1F (Form("hNElectrons_per_cluster_ch%d",Ch),Form("N Electrons per Each Cluster found - Ch %d - Sense Wire Diameter %d um - Cell Size %0.1f cm - Track Angle %0.1f - %0.1f GSa/s - Gas Mixture 80/20 %d - 90/10 %d - 85/15 %d" ,Ch, ChannelDiameter[Ch],ChannelCellSize[Ch], alpha, _gsample, isRuns_80_20, isRuns_90_10, isRuns_85_15),10,0.,10);
    hHPeaks = new TH1F (Form("hHPeaks_ch%d",Ch),Form("Height of Electron Peaks found - Ch %d - Sense Wire Diameter %d um - Cell Size %0.1f cm - Track Angle %0.1f - %0.1f GSa/s - Gas Mixture 80/20 %d - 90/10 %d - 85/15 %d" ,Ch, ChannelDiameter[Ch],ChannelCellSize[Ch], alpha, _gsample, isRuns_80_20, isRuns_90_10, isRuns_85_15),250,0,0.250);
    hHNPeaks = new TH2F (Form("hHNPeaks_ch%d",Ch),Form("Height vs N of Electron Peaks found - Ch %d - Sense Wire Diameter %d um - Cell Size %0.1f cm - Track Angle %0.1f - %0.1f GSa/s - Gas Mixture 80/20 %d - 90/10 %d - 85/15 %d" ,Ch, ChannelDiameter[Ch],ChannelCellSize[Ch], alpha, _gsample, isRuns_80_20, isRuns_90_10, isRuns_85_15),100,0,100,250,0,0.5);
    if(isChannel_1cm){
      hNPeakFPeak = new TH2F (Form("hNPeakFPeak_ch%d",Ch),Form("N of Electrons found vs T of First Electron found - Ch %d - Sense Wire Diameter %d um - Cell Size %0.1f cm - Track Angle %0.1f - %0.1f GSa/s - Gas Mixture 80/20 %d - 90/10 %d - 85/15 %d" ,Ch, ChannelDiameter[Ch],ChannelCellSize[Ch], alpha, _gsample, isRuns_80_20, isRuns_90_10, isRuns_85_15),150,0,300,250,0,250);
    }
    if(isChannel_2cm || isChannel_1p5cm){
      hNPeakFPeak = new TH2F (Form("hNPeakFPeak_ch%d",Ch),Form("N of Electrons found vs T of First Electron found - Ch %d - Sense Wire Diameter %d um - Cell Size %0.1f cm - Track Angle %0.1f - %0.1f GSa/s - Gas Mixture 80/20 %d - 90/10 %d - 85/15 %d" ,Ch, ChannelDiameter[Ch],ChannelCellSize[Ch], alpha, _gsample, isRuns_80_20, isRuns_90_10, isRuns_85_15),150,0,600,250,0,250);
    }
    if(isChannel_1cm){
	    hNClusterFCluster = new TH2F (Form("hNClusterFCluster_ch%d",Ch),Form("N of Clusters found vs T of First Cluster found - Ch %d - Sense Wire Diameter %d um - Cell Size %0.1f cm - Track Angle %0.1f - %0.1f GSa/s - Gas Mixture 80/20 %d - 90/10 %d - 85/15 %d" ,Ch, ChannelDiameter[Ch],ChannelCellSize[Ch], alpha, _gsample, isRuns_80_20, isRuns_90_10, isRuns_85_15),150,0,300,100,0,100);
    }
    if(isChannel_2cm || isChannel_1p5cm){
	    hNClusterFCluster = new TH2F (Form("hNClusterFCluster_ch%d",Ch),Form("N of Clusters found vs T of First Cluster found - Ch %d - Sense Wire Diameter %d um - Cell Size %0.1f cm - Track Angle %0.1f - %0.1f GSa/s - Gas Mixture 80/20 %d - 90/10 %d - 85/15 %d" ,Ch, ChannelDiameter[Ch],ChannelCellSize[Ch], alpha, _gsample, isRuns_80_20, isRuns_90_10, isRuns_85_15),150,0,600,100,0,100);
    }
    hTPeaks = new TH1F (Form("hTPeaks_ch%d",Ch),Form("Time of Electron Peaks found - Ch %d - Sense Wire Diameter %d um - Cell Size %0.1f cm - Track Angle %0.1f - %0.1f GSa/s - Gas Mixture 80/20 %d - 90/10 %d - 85/15 %d" ,Ch, ChannelDiameter[Ch],ChannelCellSize[Ch], alpha, _gsample, isRuns_80_20, isRuns_90_10, isRuns_85_15),400,0,800);
    if(isChannel_1cm){
	    hTFstPeaks = new TH1F (Form("hTFstPeaks_ch%d",Ch),Form("Time of First Electron Peak found - Ch %d - Sense Wire Diameter %d um - Cell Size %0.1f cm - Track Angle %0.1f - %0.1f GSa/s - Gas Mixture 80/20 %d - 90/10 %d - 85/15 %d" ,Ch, ChannelDiameter[Ch],ChannelCellSize[Ch], alpha, _gsample, isRuns_80_20, isRuns_90_10, isRuns_85_15),200,0,400);
    }
    if(isChannel_2cm || isChannel_1p5cm){
      hTFstPeaks = new TH1F (Form("hTFstPeaks_ch%d",Ch),Form("Time of First Electron Peak found - Ch %d - Sense Wire Diameter %d um - Cell Size %0.1f cm - Track Angle %0.1f - %0.1f GSa/s - Gas Mixture 80/20 %d - 90/10 %d - 85/15 %d" ,Ch, ChannelDiameter[Ch],ChannelCellSize[Ch], alpha, _gsample, isRuns_80_20, isRuns_90_10, isRuns_85_15),200,0,800);
    }
    if(isChannel_1cm){
	    hNPeaks_1 = new TH1F (Form("hNPeaks_1_ch%d",Ch),Form("Time of Last Electron Peak found - Ch %d - Sense Wire Diameter %d um - Cell Size %0.1f cm - Track Angle %0.1f - %0.1f GSa/s - Gas Mixture 80/20 %d - 90/10 %d - 85/15 %d" ,Ch, ChannelDiameter[Ch],ChannelCellSize[Ch], alpha, _gsample, isRuns_80_20, isRuns_90_10, isRuns_85_15),100,0,500);
    }
    if(isChannel_2cm || isChannel_1p5cm){
      hNPeaks_1 = new TH1F (Form("hNPeaks_1_ch%d",Ch),Form("Time of Last Electron Peak found - Ch %d - Sense Wire Diameter %d um - Cell Size %0.1f cm - Track Angle %0.1f - %0.1f GSa/s - Gas Mixture 80/20 %d - 90/10 %d - 85/15 %d" ,Ch, ChannelDiameter[Ch],ChannelCellSize[Ch], alpha, _gsample, isRuns_80_20, isRuns_90_10, isRuns_85_15),200,0,800);
    }
      
      //hFirstDeriv= new TH1F (Form("hFirstDeriv_ch%d",Ch),Form("First derivative- Ch %d - Sense Wire Diameter %d um - Cell Size %0.1f cm - Track Angle %0.1f - %0.1f GSa/s - Gas Mixture 80/20 %d - 90/10 %d - 85/15 %d" ,Ch, ChannelDiameter[Ch],ChannelCellSize[Ch], alpha, _gsample, isRuns_80_20, isRuns_90_10, isRuns_85_15),10000,-0.01,0.01);
      //hSecDeriv= new TH1F (Form("hSecDeriv_ch%d",Ch),Form("Second derivative- Ch %d - Sense Wire Diameter %d um - Cell Size %0.1f cm - Track Angle %0.1f - %0.1f GSa/s - Gas Mixture 80/20 %d - 90/10 %d - 85/15 %d" ,Ch, ChannelDiameter[Ch],ChannelCellSize[Ch], alpha, _gsample, isRuns_80_20, isRuns_90_10, isRuns_85_15),10000,-0.01,0.01);
      
      
      
    hBsl = new TH1F (Form("hBsl_ch%d",Ch),Form("Base line - Ch %d - Sense Wire Diameter %d um - Cell Size %0.1f cm - Track Angle %0.1f - %0.1f GSa/s - Gas Mixture 80/20 %d - 90/10 %d - 85/15 %d" ,Ch, ChannelDiameter[Ch],ChannelCellSize[Ch], alpha, _gsample, isRuns_80_20, isRuns_90_10, isRuns_85_15),1000,-0.8,0.); 
    if(isChannel_1cm){
	    hIntegN = new TH1F (Form("hIntegN_ch%d",Ch),Form("1 cm Integral Charge in mV * ns /Ohm (pC) - Ch %d - Sense Wire Diameter %d um - Cell Size %0.1f cm - Track Angle %0.1f - %0.1f GSa/s - Gas Mixture 80/20 %d - 90/10 %d - 85/15 %d" ,Ch, ChannelDiameter[Ch],ChannelCellSize[Ch], alpha, _gsample, isRuns_80_20, isRuns_90_10, isRuns_85_15),100,0.,250.);
    }
    if(isChannel_2cm || isChannel_1p5cm){
      hIntegN = new TH1F (Form("hIntegN_ch%d",Ch),Form("2 cm Integral Charge in mV * ns /Ohm (pC) - Ch %d - Sense Wire Diameter %d um - Cell Size %0.1f cm - Track Angle %0.1f - %0.1f GSa/s - Gas Mixture 80/20 %d - 90/10 %d - 85/15 %d" ,Ch, ChannelDiameter[Ch],ChannelCellSize[Ch], alpha, _gsample, isRuns_80_20, isRuns_90_10, isRuns_85_15),120,0.,600.);
    }
    //if(Ch<=9){
    //hIntegInR = new TH1F (Form("hIntegInR_ch%d",Ch),Form("Integral - Ch %d - Sense Wire Diameter %d um - Cell Size %0.1f cm - Track Angle %0.1f - %0.1f GSa/s - Gas Mixture 80/20 %d - 90/10 %d - 85/15 %d" ,Ch, ChannelDiameter[Ch],ChannelCellSize[Ch], alpha, _gsample, isRuns_80_20, isRuns_90_10, isRuns_85_15),400,0.,20.);
    //}
    //else {
    //  hIntegInR = new TH1F (Form("hIntegInR_ch%d",Ch),Form("Integral - Ch %d - Sense Wire Diameter %d um - Cell Size %0.1f cm - Track Angle %0.1f - %0.1f GSa/s - Gas Mixture 80/20 %d - 90/10 %d - 85/15 %d" ,Ch, ChannelDiameter[Ch],ChannelCellSize[Ch], alpha, _gsample, isRuns_80_20, isRuns_90_10, isRuns_85_15),250,0.,25.);
    //}
    //hIntegNInR = new TH1F (Form("hIntegNInR_ch%d",Ch),Form("Integral minius PDS - Ch %d - Sense Wire Diameter %d um - Cell Size %0.1f cm - Track Angle %0.1f - %0.1f GSa/s - Gas Mixture 80/20 %d - 90/10 %d - 85/15 %d" ,Ch, ChannelDiameter[Ch],ChannelCellSize[Ch], alpha, _gsample, isRuns_80_20, isRuns_90_10, isRuns_85_15),200,-10.,10.);
    //hIntegNInRC1 = new TH1F (Form("hIntegNInRC1_ch%d",Ch),Form("Integral minius PDS Norm. on NPeak - Ch %d - Sense Wire Diameter %d um - Cell Size %0.1f cm - Track Angle %0.1f - %0.1f GSa/s - Gas Mixture 80/20 %d - 90/10 %d - 85/15 %d" ,Ch, ChannelDiameter[Ch],ChannelCellSize[Ch], alpha, _gsample, isRuns_80_20, isRuns_90_10, isRuns_85_15),200,-10.,10.);
    //hIntegNInRC2 = new TH1F (Form("hIntegNInRC2_ch%d",Ch),Form("Integral minius PDS Norm. on NPeak and loss - Ch %d - Sense Wire Diameter %d um - Cell Size %0.1f cm - Track Angle %0.1f - %0.1f GSa/s - Gas Mixture 80/20 %d - 90/10 %d - 85/15 %d" ,Ch, ChannelDiameter[Ch],ChannelCellSize[Ch], alpha, _gsample, isRuns_80_20, isRuns_90_10, isRuns_85_15),200,-10.,10.);
    //hIntegNInRC3 = new TH1F (Form("hIntegNInRC3_ch%d",Ch),Form("Integral minius PDS whitout norm - Ch %d - Sense Wire Diameter %d um - Cell Size %0.1f cm - Track Angle %0.1f - %0.1f GSa/s - Gas Mixture 80/20 %d - 90/10 %d - 85/15 %d" ,Ch, ChannelDiameter[Ch],ChannelCellSize[Ch], alpha, _gsample, isRuns_80_20, isRuns_90_10, isRuns_85_15),200,-10.,20.);
    //hIntegNInRC4 = new TH1F (Form("hIntegNInRC4_ch%d",Ch),Form("Integral minius PDS Norm. on loss - Ch %d - Sense Wire Diameter %d um - Cell Size %0.1f cm - Track Angle %0.1f - %0.1f GSa/s - Gas Mixture 80/20 %d - 90/10 %d - 85/15 %d" ,Ch, ChannelDiameter[Ch],ChannelCellSize[Ch], alpha, _gsample, isRuns_80_20, isRuns_90_10, isRuns_85_15),200,-10.,20.);
    hRms = new TH1F (Form("hRms_ch%d",Ch),Form("noise RMS - Ch %d - Sense Wire Diameter %d um - Cell Size %0.1f cm - Track Angle %0.1f - %0.1f GSa/s - Gas Mixture 80/20 %d - 90/10 %d - 85/15 %d" ,Ch, ChannelDiameter[Ch],ChannelCellSize[Ch], alpha, _gsample, isRuns_80_20, isRuns_90_10, isRuns_85_15),125,0.,5.);
    //distribution of derivative for meg FE
    //hMaxV = new TH1F (Form("hMaxV_ch%d",Ch),Form("Max val - Ch %d - Sense Wire Diameter %d um - Cell Size %0.1f cm - Track Angle %0.1f - %0.1f GSa/s - Gas Mixture 80/20 %d - 90/10 %d - 85/15 %d" ,Ch, ChannelDiameter[Ch],ChannelCellSize[Ch], alpha, _gsample, isRuns_80_20, isRuns_90_10, isRuns_85_15),150,-0.02,0.3);
    //hMaxVN = new TH1F (Form("hMaxVN_ch%d",Ch),Form("Max val over base line - Ch %d - Sense Wire Diameter %d um - Cell Size %0.1f cm - Track Angle %0.1f - %0.1f GSa/s - Gas Mixture 80/20 %d - 90/10 %d - 85/15 %d" ,Ch, ChannelDiameter[Ch],ChannelCellSize[Ch], alpha, _gsample, isRuns_80_20, isRuns_90_10, isRuns_85_15),150,-0.02,0.3);
    //hMaxVNSmooth= new TH1F (Form("hMaxVNSmooth_ch%d",Ch),Form("Max val Smooth over base line - Ch %d - Sense Wire Diameter %d um - Cell Size %0.1f cm - Track Angle %0.1f - %0.1f GSa/s - Gas Mixture 80/20 %d - 90/10 %d - 85/15 %d" ,Ch, ChannelDiameter[Ch],ChannelCellSize[Ch], alpha, _gsample, isRuns_80_20, isRuns_90_10, isRuns_85_15),150,-0.02,0.3);
    hMaxVInR = new TH1F (Form("hMaxVInR_ch%d",Ch),Form("Max val - Ch %d - Sense Wire Diameter %d um - Cell Size %0.1f cm - Track Angle %0.1f - %0.1f GSa/s - Gas Mixture 80/20 %d - 90/10 %d - 85/15 %d" ,Ch, ChannelDiameter[Ch],ChannelCellSize[Ch], alpha, _gsample, isRuns_80_20, isRuns_90_10, isRuns_85_15),100,0.,0.5);
    //hMaxVNInR = new TH1F (Form("hMaxVNInR_ch%d",Ch),Form("Max val over base line - Ch %d - Sense Wire Diameter %d um - Cell Size %0.1f cm - Track Angle %0.1f - %0.1f GSa/s - Gas Mixture 80/20 %d - 90/10 %d - 85/15 %d" ,Ch, ChannelDiameter[Ch],ChannelCellSize[Ch], alpha, _gsample, isRuns_80_20, isRuns_90_10, isRuns_85_15),150,-0.02,0.3);
    
    //////////////////////////////////////
    hNPeaks->GetYaxis()->SetTitle("Entries");
    hNElectrons_per_cluster->GetYaxis()->SetTitle("Entries");
    hNElectrons_per_cluster->GetXaxis()->SetTitle("Number of Electrons per Each Cluster found");
    hNPeaks_clust->GetYaxis()->SetTitle("Entries");
    hNeventSignals->GetYaxis()->SetTitle("Entries");
    hNPeaks_clust->GetXaxis()->SetTitle("Number of Clusters found");
    hHPeaks->GetYaxis()->SetTitle("Entries");
    hHPeaks->GetXaxis()->SetTitle("Height [V]");
    hNPeakFPeak->GetYaxis()->SetTitle("Number of Electrons found");
    hNPeakFPeak->GetXaxis()->SetTitle("Time of First Electron found [ns]");
    hNClusterFCluster->GetYaxis()->SetTitle("Number of Clusters found");
    hNClusterFCluster->GetXaxis()->SetTitle("Time of First Cluster found [ns]");
    hHNPeaks->GetYaxis()->SetTitle("Height of Electron Peaks found");
    hHNPeaks->GetXaxis()->SetTitle("Number of Electron Peaks found");
    hTPeaks->GetYaxis()->SetTitle("Entries");
    hTFstPeaks->GetYaxis()->SetTitle("Entries");
    hTFstPeaks->GetXaxis()->SetTitle("Time [ns]");
    hNPeaks_1->GetYaxis()->SetTitle("Entries");
    hBsl->GetYaxis()->SetTitle("Entries");
    hBsl->GetXaxis()->SetTitle("Voltage [V]");
    //hInteg->GetYaxis()->SetTitle("Entries");
    hIntegN->GetYaxis()->SetTitle("Entries");
    hIntegN->GetXaxis()->SetTitle("Charge [pC]");
    //hIntegInR->GetYaxis()->SetTitle("Entries");
    //hIntegNInR->GetYaxis()->SetTitle("Entries");
    //hIntegNInRC1->GetYaxis()->SetTitle("Entries");
    //hIntegNInRC2->GetYaxis()->SetTitle("Entries");
    //hIntegNInRC3->GetYaxis()->SetTitle("Entries");
    //hIntegNInRC4->GetYaxis()->SetTitle("Entries");
    //hMaxV->GetYaxis()->SetTitle("Entries");
    //hMaxVN->GetYaxis()->SetTitle("Entries");
    hMaxVInR->GetYaxis()->SetTitle("Entries");
    hMaxVInR->GetXaxis()->SetTitle("Voltage [V]");
    //hMaxVNInR->GetYaxis()->SetTitle("Entries");
    //hIntegInRoriginalW->GetYaxis()->SetTitle("Entries");
    //hIntegNInRoriginalW->GetYaxis()->SetTitle("Entries");
    
    //hder->GetYaxis()->SetTitle("Entries");
    //hchiFT->GetYaxis()->SetTitle("Entries");
    //hP0->GetYaxis()->SetTitle("Entries");
    //hRmsOriginalW->GetYaxis()->SetTitle("Entries");
    //hMaxVNSmooth->GetYaxis()->SetTitle("Entries");
    hTimeDifference->GetYaxis()->SetTitle("Entries");
    hTimeDifference_clust->GetYaxis()->SetTitle("Entries");
      
    /////////////////////////////////////////////////////////
    hTimeDifference->GetXaxis()->SetTitle("Difference in time [ns]");
    hTimeDifference_clust->GetXaxis()->SetTitle("Difference in time [ns]");
    //hMaxV->GetXaxis()->SetTitle("Volt [V]");
    hNPeaks->GetXaxis()->SetTitle("N Peaks");
    hNPeaks->GetYaxis()->SetTitle("Entries");
    hBsl->GetXaxis()->SetTitle("Baseline [V]");
    hBsl->GetYaxis()->SetTitle("Entries");
    //hIntegNInR-> GetXaxis()->SetTitle("Integral [V]");
    //hIntegNInR-> GetYaxis()->SetTitle("Entries");
    //hIntegNInRC1-> GetXaxis()->SetTitle("Integral [V]");
    //hIntegNInRC1-> GetYaxis()->SetTitle("Entries");
    //hIntegNInRC2-> GetXaxis()->SetTitle("Integral [V]");
    //hIntegNInRC2-> GetYaxis()->SetTitle("Entries");
    hRms->GetXaxis()->SetTitle("Rms [mV]");
    hRms->GetYaxis()->SetTitle("Entries");
    //hMaxVNInR->GetXaxis()->SetTitle("Max_value [V]");
    //hMaxVNInR->GetYaxis()->SetTitle("Entries");
    //hIntegNInRoriginalW-> GetXaxis()->SetTitle("Integral [mA]");
    //hIntegNInRoriginalW-> GetYaxis()->SetTitle("Entries");
    hNPeaks_1->GetXaxis()->SetTitle("Time [ns]");
    
      
      
    
  }
  
  ~hstPerCh() {
    //		delete hBsl;
    //		delete hInteg;
    //		delete hIntegN;
    //		delete hRms;
    //		delete hMaxV;
    //		delete hMaxVN;
  }
};

struct hstDiffCh {//non serve
  
  TH1F *hBsl;
  TH1F *hInteg;
  TH1F *hIntegN;
  TH1F *hRms;
  TH1F *hMaxV;
  TH1F *hMaxVN;
  
  hstDiffCh(int Ch1=0, int Ch2=1, bool ped=false) {
    //gRootDir->cd();
    //TDirectory fldch(Form("H-Ch%d",Ch),Form("folder for ch%d",Ch));
    //fldch.cd();
    std::string pfix="";
    if (ped) { pfix="_ped"; }
    hBsl = new TH1F (Form("hBsl_ch%d-ch%d%s",Ch1,Ch2,pfix.c_str()),Form("Base line - Ch %d",Ch1,Ch2),1000,-0.5,0.5);
    hInteg = new TH1F (Form("hInteg_ch%d-ch%d%s",Ch1,Ch2,pfix.c_str()),Form("Integral - Ch %d ",Ch1,Ch2),600,20.,80.);
    hIntegN = new TH1F (Form("hIntegN_ch%d-ch%d%s",Ch1,Ch2,pfix.c_str()),Form("Integral minius PDS - Ch %d ",Ch1,Ch2),1000,-10.,10.);
    hRms = new TH1F (Form("hRms_ch%d-ch%d%s",Ch1,Ch2,pfix.c_str()),Form("noise RMS - Ch %d ",Ch1,Ch2),500,0,0.05);
    hMaxV = new TH1F (Form("hMaxV_ch%d-ch%d%s",Ch1,Ch2,pfix.c_str()),Form("max val - Ch %d ",Ch1,Ch2),600,-0.02,0.1);
    hMaxVN = new TH1F (Form("hMaxVN_ch%d-ch%d%s",Ch1,Ch2,pfix.c_str()),Form("max val over base line - Ch %d ",Ch1,Ch2),300,-0.02,0.1);
  }
  
  ~hstDiffCh() {
    //		delete hBsl;
    //		delete hInteg;
    //		delete hIntegN;
    //		delete hRms;
    //		delete hMaxV;
    //		delete hMaxVN;
  }
  
};//fine della struct per gli istogrammi.

typedef std::map<std::pair<int,int>, hstDiffCh> hDiffCont;
//struttura per fare la trasformata di Fourier
struct wavefft {
  int nPt;
  double *real;
  double *img;
  double *ampl;
  double *phi;
  double *omega;
  
wavefft(int npt=0) : nPt(npt) {   
  if (npt>0) {
    real = new double[npt];      //new alloca zone di memoria quando non si sa quanta memoria ci serve. Il suo complemento � delete
    img  = new double[npt];
    ampl = new double[npt];
    phi  = new double[npt];
    omega = new double[npt];
  }
  else { real = img = ampl = phi = omega =0x0; }//il puntatore punta a zero:0x0
}
  
  ~wavefft() {
    if (real != 0x0 ) { delete real; }
    if (img  != 0x0 ) { delete img; }
    if (ampl != 0x0 ) { delete ampl; }
    if (phi  != 0x0 ) { delete phi; }
    if (omega!= 0x0 ) { delete omega; }
  }
  
  void setSize(int npt=0) {
    if (npt>0 && npt!=nPt) {
      if (real != 0x0 ) { delete real; }
      if (img  != 0x0 ) { delete img; }
      if (ampl != 0x0 ) { delete ampl; }
      if (phi  != 0x0 ) { delete phi; }
      if (omega!= 0x0 ) { delete omega; }
      nPt=npt;
      real = new double[npt];
      img  = new double[npt];
      ampl = new double[npt];
      phi  = new double[npt];
      omega = new double[npt];
    }
  }
  
  void clear() {
    for (int iPt=0; iPt<nPt; ++iPt) {
      real[iPt]=0.0;
      img[iPt]=0.0;
      ampl[iPt]=0.0;
      phi[iPt]=0.0;
      omega[iPt]=0.0;
    }
  }
  
  void fillRM(int npt, double *Real, double *Img, double *Omega) {
    if (nPt!=npt) {setSize(npt);}
    for (int iPt=0; iPt<nPt; ++iPt) {
      real[iPt]=Real[iPt];
      img[iPt]=Img[iPt];
      ampl[iPt]=TMath::Sqrt(Real[iPt]*Real[iPt]+Img[iPt]*Img[iPt]);
      phi[iPt]=(Img[iPt]||Real[iPt])?TMath::ATan2(Img[iPt],Real[iPt]):0.0;
      omega[iPt]=Omega[iPt];
      //if (Omega!=0x0) { omega[iPt]=Omega[iPt];}
      //else { omega[iPt]=TMath::TwoPi()*((Double_t) (iPt+1))/tmax; }
    }
  }
  
  void fillAP(int npt, double *Ampl, double *Phi, double *Omega) {
    if (nPt!=npt) {setSize(npt);}
    for (int iPt=0; iPt<nPt; ++iPt) {
      ampl[iPt]=Ampl[iPt];
      phi[iPt]=Phi[iPt];
      omega[iPt]=Omega[iPt];
      real[iPt]=ampl[iPt]*TMath::Cos(phi[iPt]);
      img[iPt]=ampl[iPt]*TMath::Sin(phi[iPt]);
    }
  }
  
};//fine struttura wavefft.
typedef std::map<int,wavefft> fftCont; //mappa per le forme d'onda filtrate.


#endif
