/**  
 *
 *  Authors: B. D'Anzi - University and INFN Bari
 * 			F. Cuna - University and INFN Lecce
 *
 **/
#define read_data_cxx
#include "read_data.h"
#include "funcUtils.h"
#include "FindPeak-algo.C"
#include "Clusterization.C"
#include "TPaveText.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TLegend.h>
#include <TLine.h> 
#include<TVector.h>
#include "TF1.h"
#include "TDirectory.h"    
//#include "TVirtualFFT.h"

#include "TMath.h"
#include "TFitResult.h"
#include "TMinuit.h"
#include "TGraphErrors.h"
//#include <map>
#include <fstream>
#include <string>
#include "TSpectrum.h"
#include "TAxis.h"
#include "TVirtualFitter.h"
#include "TMarker.h"


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// 														Define Data Containers															////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
vector<float>X;//definisce l'asse x

// Temporary container for tmp Graphs and Canvas
std::vector<TGraph *> tmpWaves;       //grafico delle forme d'onda
std::vector<TCanvas *> tmpCv;         //canvas forme d'onda
std::vector<TGraph *> tmpFFTAmpl;     //grafico per l'ampiezza
std::vector<TGraph *> tmpFFTFltAmpl;  //grafico per ampiezza filtrata
std::vector<TGraph *> tmpFFTphi;      //grafico per la fase
std::vector<TCanvas *> tmpCvFFT;      //canvas per la fft
std::vector<TGraph *> tmpFltWaves;    //grafici per funzioni d'onda filtrate (primo filtro:retta lineare di taglio)
std::vector<TCanvas *> tmpCvFlt;
//canvas forme d'onda filtrate(primo filtro:retta lineare di taglio)
std::vector<TH1F *>  tmpFFTAmplPeak;  //grafico per lo spettro in frequenza dell'ampiezza con filtro di smooth
std::vector<TCanvas *> tmpCvPeak_1;
std::vector<TH1F *>  tmpFFTAmplPeak_1;  //grafico per lo spettro in frequenza dell'ampiezza con filtro di smooth (0-50MHz)
//std::vector<TGraph *> tmpnoise_1;   //grafici per il noise individuato dal filtro RC.(singolo filtro)
//std::vector<TCanvas *>tmpCvnoise_1; //canvas per il noise individuato dal filtro RC(singolo filtro).
std::vector<TGraph *> tmpflt_batt;    //grafici per il segnale con noise di battimento+fitt
std::vector<TCanvas *> tmpCvflt_batt; //canvas per il segnale con noise di battimento+fitt
std::vector<TCanvas *> tmpCvPeak;
std::vector<TGraph *>tmpflt_batt_2; //secondo battimento
std::vector<TCanvas *>tmpCvflt_batt_2;
std::vector<TGraph *> tmpsignal_1;    //grafici per il segnale finale dopo i filtri
std::vector<TCanvas *> tmpCvsignal_1; //canvas per il segnale finale dopo i filtri
std::vector<TGraph *> tmpNoise_total;  //grafici per il noise totale
std::vector<TCanvas *>tmpCvNoise_total; //canvas per il noise totale
std::vector<float> X_negative; //vettore X per la parte negativa della waveform
std::vector<float> Y_negative; //vettore Y per la parte negativa della waveform
std::vector<TGraph *> neg_wavefbat; //grafico per waveform negative di battimento
std::vector< TGraphErrors *> tmpNegBatt; //grafico per le wave di battimento negative con errore
std::vector<TCanvas *> tmpCvNegBatt;
std::vector<TCanvas *> tmpCvMidNotch;
std::vector<TGraph *> tmpMid_Notch;
std::vector<float> N_signalevents;
std::vector<float>::iterator it; 



////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////														MAIN LOOP function																										////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void read_data::Loop(Char_t *output, Int_t MidEv,Int_t eventn,  Float_t _gsample, Float_t N_1, Float_t N_2, Float_t N_3, Float_t N_4, Float_t bslnTimeInterval, Int_t _dim, Float_t _scale_cut)
{
  isNov2021TestBeam = false;
  isJuly2022TestBeam = true;
  dim = _dim; // to be changed
  float scale_cut = _scale_cut;
  //_tmax = 853.3e-9;//1.0e-6;//853.3e-9;//1.0e-6;//853.3e-9;//1.0e-6;
  _tmax = (float) 1/(_gsample)*dim;
  if(isNov2021TestBeam){
    nMaxCh = 12;
  }
  else if(isJuly2022TestBeam){
    nMaxCh = 15;
  }
  tmax=_tmax;  //to fix compilation problem
  double timeRes;
  int isl = 0;
#ifndef _OSC
  if(_gsample<=1.001){
    isl = 0;
  }else if(_gsample<=1.2001)
    {
      isl = 1; 
    }
  else if(_gsample<=1.5001)
    {
      isl = 2;
    }
  else if(_gsample<=2.001)
    {
      isl = 3;
    }
  invfbin=1.0/((float)fbin_rms);
#else
  // skipFstBin=250; //275;//set Bins to Skip
  // skipLstBin=300;//375
#endif
  
  /*********** Histogram creation***********/
  
  char *basename(char *path);
  cout << "Basename " << basename(output) << endl;
  Char_t outChar[500];
  Char_t outChar_complete[500];
  TString out = basename(output);
  out.ReplaceAll(".root","");
#ifndef _OSC
  sprintf(outChar,"histosTB_%s.root",out.Data()); 
  sprintf(outChar_complete,"histosTB_%s_N1_%0.1f_N2_%0.1f_N3_%0.1f_N4_%0.1f_cut_scale_%0.2f_sampling_%0.1f.root",out.Data(),N_1,N_2,N_3,N_4,scale_cut,_gsample);   //Data sono i dati convertiti nei file.root
#else
  sprintf(outChar,"histosOSC_%s.root",out.Data());  //scrive il titolo della rootupla
#endif
  cout << "The output file is " << outChar_complete << endl;
  TFile * theFile = new TFile(outChar_complete,"RECREATE");  //crea un nuovo file.root e se ne esiste già uno lo sovrascrive
  
  std::map<int, hstPerCh *> HstPerCh;               //mappa per istogrammi sui segnali filtrati.
  int count =0;
  
  
  /***********Directory creation***********/
  theFile->cd("/");
  TDirectory *waveDir = theFile->mkdir("Waves");
  theFile->cd("/");
  TDirectory *waveFFTDir = theFile->mkdir("WavesFFT");
  theFile->cd("/");
  TDirectory *waveFltDir = theFile->mkdir("WavesFlt");//reham 
  theFile->cd("/");
  TDirectory *signal = theFile->mkdir("signal_Afterflt"); //directory per la waveform di segnale dopo tutti i filtri
  
  
  
  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  
  //                                               //  Loop on entries //                                                  //
  
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  
  std::vector<int> trigCh; //vettore di interi per i canali di trigger.
  trigCh.clear();
#ifndef _OSC // For test beam usage
  if(isNov2021TestBeam){
    trigCh.push_back(0);
    trigCh.push_back(1);
    trigCh.push_back(2);
    trigCh.push_back(3);
  }
  else if(isJuly2022TestBeam){
    trigCh.push_back(12);
    trigCh.push_back(13);
  trigCh.push_back(14);
  trigCh.push_back(15);
  }
#else //for muon data - Oscilloscope in the lab
  trigCh.push_back(7);//osc1
  trigCh.push_back(8);//osc2
#endif
  
  
  int nSig=3; // numero di sigma, not used
  std::map<int, int> nSelEv;
  nSelEv.clear();
  
  if (fChain == 0) return; //TTree pointer
  
  Long64_t nentries = fChain->GetEntriesFast(); //Return the number of entries as of the last check.
  cout << "Total Number of entries: " << nentries << endl;
  //  Long64_t nbytes = 0, nb = 0;
  
  Int_t firstEv=0;                              //inizializzo il primo evento a zero
  Int_t lastEv=nentries;                        //l'ultimo evento coincide con le nentries, che si prende con la funzione get...sopra
  //Int_t MidEv;
  if (eventn>=0&&eventn<lastEv) {
    firstEv=eventn;
    lastEv=eventn+1;
  }else if(MidEv>=0&&eventn<-1 ){
    firstEv=MidEv;
    lastEv=-eventn+firstEv;
  }else if (eventn<-1 && eventn>-nentries){
    lastEv=-eventn;}
  //cout<< " MidEv "<< MidEv<< "lastEv "<<lastEv<<endl;
  
  Float_t alpha= 0.;
  Float_t cos_alpha = 0.;
  Float_t expected_electrons =0.;
  Float_t expected_cluster =0.;
  Float_t cluster_per_cm_mip = 18.;
  Float_t drift_size =0.;
  Float_t relativistic_rise = 1.3;
  Float_t cluster_population = 1.6;
  Int_t NPeak;
  string name_file_compact = outChar;
  // if(isNov2021TestBeam){
  // string Runs_80_20[] = {"histosTB_run_127.root", "histosTB_run_117.root"}; // 2021 Nov Test Beam
  // string Runs_90_10[] = {"histosTB_run_72.root", "histosTB_run_73.root", "histosTB_run_74.root", "histosTB_run_86.root", "histosTB_run_87.root", "histosTB_run_88.root", "histosTB_run_89.root", "histosTB_run_90.root", "histosTB_run_91.root", "histosTB_run_92.root", "histosTB_run_93.root", "histosTB_run_94.root", "histosTB_run_95.root", "histosTB_run_96.root", "histosTB_run_97.root" ,"histosTB_run_98.root", "histosTB_run_99.root", "histosTB_run_100.root", "histosTB_run_101.root"};
  // string Runs_85_15[] = {"histosTB_run_0.root"}; // No data for Nov 2021, 2022 Test Beam to be inserted
  // string Runs_alpha_0[] = {"histosTB_run_99.root", "histosTB_run_117.root", "histosTB_run_86.root", "histosTB_run_100.root", "histosTB_run_72.root", "histosTB_run_73.root", "histosTB_run_74.root"}; // 2021 Test Beam
  // string Runs_alpha_15[] = {"histosTB_run_98.root", "histosTB_run_87.root", "histosTB_run_97.root"}; // 2021 Nov Test Beam
  // string Runs_alpha_30[] = {"histosTB_run_96.root", "histosTB_run_88.root", "histosTB_run_95.root"}; // 2021 Nov Test Beam
  // string Runs_alpha_45[] = {"histosTB_run_94.root", "histosTB_run_89.root", "histosTB_run_93.root"}; // 2021 Nov Test Beam
  // string Runs_alpha_60[] = {"histosTB_run_91.root", "histosTB_run_127.root", "histosTB_run_90.root" ,"histosTB_run_92.root"}; // 2021 Nov Test Beam
  // }
  // else if(isJuly2022TestBeam){
  //string Runs_80_20[] = {"histosTB_run_39.root","histosTB_run_40.root", "histosTB_run_41.root", "histosTB_run_42.root","histosTB_run_43.root","histosTB_run_44.root","histosTB_run_45.root","histosTB_run_46.root"}; // 2022 July Test Beam
  //Scan Study: HV
  //string Runs_80_20[] = {"histosTB_run_40.root","histosTB_run_41.root","histosTB_run_42.root","histosTB_run_43.root"};
  //Scan Study: Angle
  //string Runs_80_20[] = {"histosTB_run_44.root", "histosTB_run_45.root","histosTB_run_41.root","histosTB_run_46.root"};
  //Scan Study: gas mixture
  string Runs_80_20[] = {"histosTB_run_41.root"};
  // Scan Study: Sampling rate
  string Runs_90_10[] = {"histosTB_run_9.root","histosTB_run_10.root"}; // 2022 July Test Beam
  string Runs_85_15[] = {"histosTB_run_63.root"}; // 2022 July Test Beam
  string Runs_alpha_0[] = {"histosTB_run_44.root"}; // 2022 July Test Beam
  string Runs_alpha_15[] = {"histosTB_run_0.root"}; // We didn't take data at 15°, 2022 July Test Beam
  string Runs_alpha_30[] = {"histosTB_run_45.root"}; // 2022 July Test Beam
  string Runs_alpha_45[] = {"histosTB_run_39.root","histosTB_run_40.root", "histosTB_run_41.root", "histosTB_run_42.root","histosTB_run_43.root","histosTB_run_9.root","histosTB_run_10.root","histosTB_run_63.root"}; // 2022 July Test Beam
  string Runs_alpha_60[] = {"histosTB_run_46.root"}; // 2022 July Test Beam
  //   }
  int Runs_80_20_size = sizeof Runs_80_20 / sizeof Runs_80_20[0];
  int Runs_90_10_size = sizeof Runs_90_10 / sizeof Runs_90_10[0];
  int Runs_85_15_size = sizeof Runs_85_15 / sizeof Runs_85_15[0];
  int Runs_alpha_0_size = sizeof Runs_alpha_0 / sizeof Runs_alpha_0[0];
  int Runs_alpha_15_size =sizeof Runs_alpha_15 / sizeof Runs_alpha_15[0];
  int Runs_alpha_30_size = sizeof Runs_alpha_30 / sizeof Runs_alpha_30[0];
  int Runs_alpha_45_size = sizeof Runs_alpha_45 / sizeof Runs_alpha_45[0];
  int Runs_alpha_60_size = sizeof Runs_alpha_60 / sizeof Runs_alpha_60[0];
  int isRuns_80_20 = 0;
  int isRuns_90_10 = 0;
  int isRuns_85_15 = 0;
  int isRuns_alpha_0 = 0;
  int isRuns_alpha_15 = 0;
  int isRuns_alpha_30 = 0;
  int isRuns_alpha_45 = 0;
  int isRuns_alpha_60 = 0;
  for (int i = 0; i < Runs_80_20_size; i++) {
    if (Runs_80_20[i] == name_file_compact) {
      isRuns_80_20 = 1;
      cluster_per_cm_mip = 18.;
      printf("Gas mixture changed to 80/20 with Number of Cluster/cm (MIP) to be %0.1f!\n",cluster_per_cm_mip);
      break;
    }
  }
  for (int i = 0; i < Runs_90_10_size; i++) {
    if (Runs_90_10[i] == name_file_compact) {
      isRuns_90_10 = 1;
      cluster_per_cm_mip = 12.;
      printf("Gas mixture changed to 90/10 with Number of Cluster/cm (MIP) to be %0.1f!\n",cluster_per_cm_mip);
      break;
    }
  }
  for (int i = 0; i < Runs_85_15_size; i++) {
    if (Runs_85_15[i] == name_file_compact) {
      isRuns_85_15 = 1;
      cluster_per_cm_mip = 15.;
      printf("Gas mixture changed to 85/15 with Number of Cluster/cm (MIP) to be %0.1f!\n",cluster_per_cm_mip);
      break;
    }
  }
  for (int i = 0; i < Runs_alpha_0_size; i++) {
    if (Runs_alpha_0[i] == name_file_compact) {
      isRuns_alpha_0 = 1;
      alpha = 0.;
      printf("Track Angle Changed to be %0.1f!\n",alpha);
      break;
    }
  }
  for (int i = 0; i < Runs_alpha_15_size; i++) {
    if (Runs_alpha_15[i] == name_file_compact) {
      isRuns_alpha_15 = 1;
      alpha = 15.;
      printf("Track Angle Changed to be %0.1f!\n",alpha);
      break;
    }
  }
  for (int i = 0; i < Runs_alpha_30_size; i++) {
    if (Runs_alpha_30[i] == name_file_compact) {
      isRuns_alpha_30 = 1;
      alpha = 30.;
      printf("Track Angle Changed to be %0.1f!\n",alpha);
      break;
    }
  }
  for (int i = 0; i < Runs_alpha_45_size; i++) {
    if (Runs_alpha_45[i] == name_file_compact) {
      isRuns_alpha_45 = 1;
      alpha = 45.;
      printf("Track Angle Changed to be %0.1f!\n",alpha);
      break;
    }
  }
  for (int i = 0; i < Runs_alpha_60_size; i++) {
    if (Runs_alpha_60[i] == name_file_compact) {
      isRuns_alpha_60 = 1;
      alpha = 60.;
      printf("Track Angle Changed to be %0.1f!\n",alpha);
      break;
    }
  }
  
  cos_alpha = TMath::Cos(alpha*TMath::DegToRad());
  
  
  for (Long64_t jentry=firstEv; jentry<lastEv;jentry++) {
    
    
    N_signalevents.assign(nMaxCh+1,0.0);
    
    
    
    //cout << jentry << "\n" << endl;
    //cout << "Initial channel hits: \n"; 
    // for (int i = 0; i < N_signalevents.size(); i++) 
    //   cout << "Ch: " << i << " Hits: " << N_signalevents[i] << "\n"; 
    // cout << "\n";
    
    
    Long64_t ientry = LoadTree(jentry); //The function finds the corresponding Tree and returns the entry number in this tree.  
    if (ientry < 0) break;
    //nb = fChain->GetEntry(jentry);   nbytes += nb;
    fChain->GetEntry(jentry);
    WvCont Waves;
    WvCont FltWaves;
    fftCont Wffts;
    WvCont tmpWSG_signal;
    WvCont SmtSGWaves;
    WvCont Waves_signal_1;
    
    FltWaves.clear();
    Wffts.clear();
    Waves.clear();
    Waves.clear();
    tmpWSG_signal.clear();
    SmtSGWaves.clear();
    
    //.root files reading
    
    if(jentry==firstEv) {
      // Units are mV and ns !!!
      // 1024 channel of ADC = 0.833*1024 ns // => every channel is 0.853 nano sec 
      // _tmax = (float) 1/(_gsample*1.0e+9)*dim;  
      // timeRes = (_tmax*1.0e+9)/((float)dim);//in nano sec = 0.833s // bin
      timeRes = (float) 1/(_gsample);
      //cout<< "timeRes "<<timeRes;
      //cout<< "_tmax "<<_tmax;
      X.clear();
      //cout<< "dim "<<dim<<" time "<<timeRes<<endl;
      for (int i = 0; i < dim; ++i) { X.push_back(timeRes *(i+1));  }
    }
    
    cout << "Event analyzed: #" << jentry << "\r";
    fflush(stdout);
    
    bool firstEntering=true;
    bool firstEntering_filter=true;
    int counting=0;
    int counting_filter=0;
    int nTriggerChannels=0;
    int counter=0;
    int counter_1cm=0;
    int counter_1p5cm=0;
    int counter_2cm=0;
    // if(isNov2021TestBeam){
    // int Channel_1cm[] = {4,5,6,7,8,9}; // Nov 2021 Test Beam
    // int Channel_2cm[] = {10,11,12}; // Nov 2021 Test Beam
    // int Channel_1p5cm[] = {0}; // NO Nov 2021 Test Beam
    // }
    // else if(isJuly2022TestBeam){
    int Channel_1cm[] = {1,2,3,5,6,8,9,10}; // July 2022 Test Beam
    int Channel_2cm[] = {15}; // NO in July 2022 Test Beam
    int Channel_1p5cm[] = {0,4,7,11}; // New Test Beam
    // }
    int Channel_1cm_size = sizeof Channel_1cm / sizeof Channel_1cm[0];
    int Channel_2cm_size = sizeof Channel_2cm / sizeof Channel_2cm[0];
    int Channel_1p5cm_size = sizeof Channel_1p5cm / sizeof Channel_1p5cm[0];
    int isChannel_1cm = 0;
    int isChannel_2cm = 0;
    int isChannel_1p5cm = 0;
    // if(isNov2021TestBeam){
    // int top_first_driftTubes_line[] = {4,5,6,7,8,9}; // 2021 Nov Test Beam
    // int top_second_driftTubes_line[] = {13,14}; // 2021 Nov Test Beam
    // int top_third_driftTubes_line[] = {10,11,12}; // 2021 Nov Test Beam
    // }
    // else if(isJuly2022TestBeam){
    int top_first_driftTubes_line[] = {8,6,5,1}; // 2022 July Test Beam all 1.0 cm upstream downstream
    int top_second_driftTubes_line[] = {7,4,11}; // 2022 July Test Beam 1.5 cm
    //int top_third_driftTubes_line[] = {0}; // 2022 July Test Beam 1.5 cm
    int top_third_driftTubes_line[] = {3,2,10,9}; // 2022 July Test Beam 1.0 cm
    // }
    int top_first_driftTubes_line_size = sizeof top_first_driftTubes_line / sizeof top_first_driftTubes_line[0];
    int top_second_driftTubes_line_size = sizeof top_second_driftTubes_line / sizeof top_second_driftTubes_line[0];
    int top_third_driftTubes_line_size = sizeof top_third_driftTubes_line / sizeof top_third_driftTubes_line[0];
    int istop_first_driftTubes_line = 0;
    int istop_second_driftTubes_line = 0;
    int istop_third_driftTubes_line = 0;
    int counter_top_first_driftTubes_line_size = 0;
    int counter_top_second_driftTubes_line_size = 0;
    int counter_top_third_driftTubes_line_size = 0;
    // Loop to fill the waveforms  - Customization for Jul 2022 and Nov 2021 Beam Tests
    for (auto point : wd->getX742Data()) {
      
      for (int channel=0; channel<=nMaxCh; channel++){
	
	for(auto trgCh : trigCh) { if (channel==trgCh && firstEntering) {nTriggerChannels++;}}
	
      } //get the number of TriggerChannels
      
      int channel = point.first;
      
      
      bool isTrg=false;
      for(auto trgCh : trigCh) { if (channel==trgCh) { isTrg=true;}}
      
      if (!isTrg && channel<=nMaxCh) { 
	
	Waves[channel].fillWave(point.second,dim, _gsample, bslnTimeInterval);
	//cout << "Rms after normalization(V):" << Waves[channel].rms << endl;
	//cout << "Integ from charge_integInR after normalization(pC):" << Waves[channel].charge_integInR << endl;
	//cout << "Integ from charge_integ before normalization(pC):" << Waves[channel].charge_integ << endl;
	//cout << "Bsln before normalization(V):" << Waves[channel].bsln << endl;
	//cout << "Gsample before normalization:" << _gsample << endl;
	//for(int i=0; i<Waves[channel].nPt();i++){
	//  Waves[channel].Y[i]=Waves[channel].Y[i]-Waves[channel].bsln;
	//  //cout << "Here is the update waveform: "<< Waves[channel].Y[i] << endl; 
	//}
	isChannel_2cm = 0;
	isChannel_1cm = 0;
	isChannel_1p5cm = 0;
	istop_first_driftTubes_line = 0;
	istop_second_driftTubes_line = 0;
	istop_third_driftTubes_line = 0;
	int channel = point.first;
	
	for (int i = 0; i < Channel_1cm_size; i++) {
	  if (Channel_1cm[i] == channel) {
            isChannel_1cm = 1;
            break;
	  }
	}
	for (int i = 0; i < Channel_2cm_size; i++) {
	  if (Channel_2cm[i] == channel) {
            isChannel_2cm = 1;
            break;
	  }
	}
	for (int i = 0; i < Channel_1p5cm_size; i++) {
	  if (Channel_1p5cm[i] == channel) {
            isChannel_1p5cm = 1;
            break;
	  }
	}
	for (int i = 0; i < top_first_driftTubes_line_size; i++) {
	  if (top_first_driftTubes_line[i] == channel) {
            istop_first_driftTubes_line = 1;
            break;
	  }
	}
	for (int i = 0; i < top_second_driftTubes_line_size; i++) {
	  if (top_second_driftTubes_line[i] == channel) {
            istop_second_driftTubes_line = 1;
            break;
	  }
	}
	for (int i = 0; i < top_third_driftTubes_line_size; i++) {
	  if (top_third_driftTubes_line[i] == channel) {
            istop_third_driftTubes_line = 1;
            break;
	  }
	}
	
	// Cout chekcs
	//cout << "nPt in Waves:" << Waves[channel].nPt() << endl;
	//cout << "Rms in Waves:" << Waves[channel].rms << endl;
	//cout << "integ in Waves:" << Waves[channel].integ << endl;
	//cout << "Gigasample:" << _gsample << endl;
	//Waves[channel].fillWave((Waves[channel].Y),Waves[channel].nPt(), _gsample);
	//cout << "Rms after normalization(V):" << Waves[channel].rms << endl;
	//cout << "Integ after normalization(pC):" << Waves[channel].charge_integInR << endl;
	//cout << "Bsln after normalization(V):" << Waves[channel].bsln << endl;
	if (!isTrg && channel<=nMaxCh && Waves[channel].max>10*(Waves[channel].rms) && firstEntering && ((wave)Waves[channel]).charge_integInR > 20) {
	  if(isChannel_1cm){
	    counter_1cm=counter_1cm+1;
	  }
	  
	  if(istop_first_driftTubes_line){
	    counter_top_first_driftTubes_line_size = counter_top_first_driftTubes_line_size+1;
	  }
	  
	  if(istop_second_driftTubes_line){
	    counter_top_second_driftTubes_line_size = counter_top_second_driftTubes_line_size+1;
	  }
	  
	  if(istop_third_driftTubes_line){
	    counter_top_third_driftTubes_line_size = counter_top_third_driftTubes_line_size+1;
	  }

	  if(isChannel_2cm){
	    counter_2cm=counter_2cm+1;
	  }
	  
	  if(isChannel_1p5cm){
	    counter_1p5cm=counter_1p5cm+1;
	  }
	} //get the number of Signal 1 cm and 1.5 cm channels
	
      }
      //Waves[channel].fillWave(FltWaves[channel].Y,FltWaves[channel].nPt());
      //cout << "Event: " << jentry << " Number of 1 cm channels hit: " << counter_1cm <<"\n";
      //cout << "Event: " << jentry << " Number of 1.5 cm channels hit: " << counter_1.5cm <<"\n";
    }
    
    // End of fill of Waveforms
    
    for (auto point : wd->getX742Data()) {
      counter=counter+1;
      //cout << "Counter:"<< counter <<"\n"<<endl;
      //cout << "Channel analysed from the map:" << point.first <<endl;
      
      for (int channel=0; channel<=nMaxCh; channel++){
	
	for(auto trgCh : trigCh) { if (channel==trgCh && firstEntering) {nTriggerChannels++;}}
	
      } //get the number of TriggerChannels
      isChannel_2cm = 0;
      isChannel_1cm = 0;
      isChannel_1p5cm = 0;
      int channel = point.first;
      
      for (int i = 0; i < Channel_1cm_size; i++) {
        if (Channel_1cm[i] == channel) {
	  isChannel_1cm = 1;
	  break;
        }
      }
      for (int i = 0; i < Channel_2cm_size; i++) {
        if (Channel_2cm[i] == channel) {
	  isChannel_2cm = 1;
	  break;
        }
      }
      
      for (int i = 0; i < Channel_1p5cm_size; i++) {
        if (Channel_1p5cm[i] == channel) {
	  isChannel_1p5cm = 1;
	  break;
        }
      }
      

      bool isTrg=false;
      for(auto trgCh : trigCh) { if (channel==trgCh) { isTrg=true;}}
      
      if (!isTrg && channel<=nMaxCh) { 
	if ( HstPerCh.find(channel)==HstPerCh.end() ) {
	  TDirectory *chDir = theFile->mkdir(Form("H-Ch%d_signal",channel));
	  chDir->cd();
	  HstPerCh.insert(make_pair(channel,new hstPerCh(channel,isChannel_1cm,isChannel_2cm, isChannel_1p5cm, alpha,_gsample, isRuns_80_20, isRuns_90_10, isRuns_85_15)));
	  theFile->cd("/");
	}
	//cout << "Channel is = " << channel << endl;
	Waves[channel].fillWave(point.second,dim,_gsample, bslnTimeInterval);
	/****-------------------------------TRASFORMATA DI FOURIER---------------------------------*****/
	
	//FFT(Waves[channel],Wffts[channel]);//trasformata di Fourier not used!!!
	//double *realFltFFT = new double[Waves[channel].nPt()];
	//double *imgFltFFT = new double[Waves[channel].nPt()];
	//
	//filterWaveBsl(Wffts[channel],realFltFFT,imgFltFFT);//Baseline filter.  not used!                                   ///
	//InverseFFT(realFltFFT,imgFltFFT,Waves[channel].nPt(),FltWaves[channel]);//trasformata inversa
	//for(int i=0; i<Waves[channel].nPt();i++){
	//  Waves[channel].Y[i]=Waves[channel].Y[i]-Waves[channel].bsln;
	//}
	//Waves[channel].fillWave((Waves[channel].Y),Waves[channel].nPt(),_gsample);
	bool saveEvents = false; // waveFltDir is not used!!
	if (saveEvents) {
	  waveFltDir->cd();
	  tmpCvFlt.push_back( new TCanvas(Form("CvFlt-Ch%d_ev%d",channel,jentry),Form("tmpFltWave-Ch%d_ev%d",channel,jentry)) );
	  tmpCvFlt.back()->cd();
	  tmpFltWaves.push_back( new TGraph ( dim, &X[0], &FltWaves[channel].Y[0]) );
	  tmpFltWaves.back()->GetXaxis()->SetTitle("time [ns]");
	  tmpFltWaves.back()->SetTitle(Form("tmpFltWave-Ch%d_ev%d",channel,jentry));
	  tmpFltWaves.back()->GetYaxis()->SetTitleOffset(1.4);
	  tmpFltWaves.back()->GetYaxis()->SetTitle("Volt");
	  tmpFltWaves.back()->GetYaxis()->SetRangeUser(-0.1,0.4);
	  tmpFltWaves.back()->Draw("AP");
	  tmpCvFlt.back()->Write();
	  theFile->cd("/");
	  
	} 
	bool saveWave=false; //non salvo, quindi non vediamo l'out di questo pezzo di codice.
	saveEvents=false;
	
	if (saveWave) {
	  waveDir->cd(); // waveDir is not used!!
	  tmpCv.push_back(new TCanvas(Form("Cv-Ch%d_ev%d",channel,jentry),Form("tmpWave-Ch%d_ev%d",channel,jentry)));
	  tmpCv.back()->cd();
	  tmpWaves.push_back(new TGraph (dim, &X[0], &Waves[channel].Y[0]));
	  tmpWaves.back()->GetXaxis()->SetTitle("time [ns]");
	  tmpWaves.back()->SetTitle(Form("tmpWave-Ch%d_ev%d",channel,jentry));
	  tmpWaves.back()->GetYaxis()->SetTitleOffset(1.4);
	  tmpWaves.back()->GetYaxis()->SetTitle("Voltage [V]");
	  tmpWaves.back()->SetMarkerSize(1);
	  tmpWaves.back()->SetMarkerStyle(2);
	  tmpWaves.back()->Draw("APL");
	  tmpCv.back()->Write();
	  theFile->cd("/");
	}
	
	if (saveEvents) {
	  waveDir->cd(); // waveDir is not used!!
	  if(firstEntering){
	    tmpCv.push_back( new TCanvas(Form("Cv-ev%d",jentry),Form("tmpWave-ev%d",jentry)) );
	    tmpCv.back()->Divide(3,4);
	    firstEntering=false;
	  }
	  tmpCv.back()->cd(channel-nTriggerChannels+1);
	  //cout<< "nTriggerchannels"<<nTriggerChannels<<" event "<<jentry<<endl;
	  //cout<< "channel"<<channel<<" event "<<jentry<<endl;
	  tmpWaves.push_back( new TGraph ( dim, &X[0], &Waves[channel].Y[0]) );
	  tmpWaves.back()->GetXaxis()->SetTitle("time [ns]");
	  tmpWaves.back()->SetTitle(Form("tmpWave-Ch%d_ev%d",channel,jentry));
	  tmpWaves.back()->GetYaxis()->SetTitleOffset(1.4);
	  tmpWaves.back()->GetYaxis()->SetTitle("Voltage [V]");
	  tmpWaves.back()->Draw("AP");
	  counting++;
	  if(counting==(nMaxCh-nTriggerChannels+1)){
	    tmpCv.back()->Write();
	  }
	  theFile->cd("/");
	}
	
	
      } //if !nottrigger
      
      theFile->cd("/");
      
      std::vector< std::pair<int,int> > chToDiff;  //vettore di pair: canali da sottrarre
      chToDiff.push_back( make_pair(1,2) );
      
      /////New /////
      
      Int_t pkPos[250];
      Float_t pkHgt[250];
      Int_t nElectrons_per_cluster[250];
      Float_t cut_cluster_ns = 3.5;
      
      Int_t pkPos_1[250];
      Float_t pkHgt_1[250];
      //Clusterization variables//
      Int_t NPeak_clust = 0;
      Int_t pkPos_clust[250];
      Int_t average_pkPos_clust[250];
      Float_t average_pkHgt_clust[250];
      Float_t pkHgt_clust[250];
      
      //maps for full waves
      vector<pair <int, int > > waveFull;
      waveFull.clear();
      pair<int,int> Full;
      NPeak=0;
      
      /* Adding SG filter smoothing
	 int m,k;
	 m=13; //number of bin interested by the SG smoothing 
	 k=3; //order of the polinomial used
	 std::vector<float> tmpWSG_signal_23=smoothSG(Waves[channel].Y,m,k);
	 tmpWSG_signal[channel].fillWave(tmpWSG_signal_23);
      */
      
    //   float scaleInt=1.0;
      
      //if(!isTrg && Waves[channel].max>10*(Waves[channel].rms) && channel<=nMaxCh && ((channel<=9 && counter_1cm>=4) || (channel>=10 && counter_2cm>=3))&&((wave)Waves[channel]).nnIntegInR()>0.1){ //nPtInR == Y.size - first,lastBin; search peak when max amplitude > 5 mV
      //if(!isTrg && Waves[channel].max>10 * Waves[channel].rms){ //nPtInR == Y.size - first,lastBin; search peak when max amplitude > 5 mV
      //if(!isTrg && Waves[channel].max> 10 * Waves[channel].rms ){
      //cout << "Waveform max:" << Waves[channel].max << endl; 
		//cout << "Waveform rms x 10:" << 10 * Waves[channel].rms << endl; 
      //if(!isTrg && Waves[channel].max > 0.05 && ((isChannel_1cm && Waves[channel].charge_integInR > 20 && Waves[channel].charge_integInR < 160) || ((isChannel_2cm || isChannel_1p5cm) && Waves[channel].charge_integInR > 60 && Waves[channel].charge_integInR < 200)) && (counter_top_first_driftTubes_line_size>=2 || counter_top_third_driftTubes_line_size>=2 || counter_top_second_driftTubes_line_size>=2)){  // we are dealing with Volts , charge in pC
      if(!isTrg && Waves[channel].max > 0.05 && ((isChannel_1cm && Waves[channel].charge_integInR > 20 ) || ((isChannel_2cm || isChannel_1p5cm) && Waves[channel].charge_integInR > 60)) && (counter_top_first_driftTubes_line_size>=2 || counter_top_third_driftTubes_line_size>=2 || counter_top_second_driftTubes_line_size>=2)){  // we are dealing with Volts , charge in pC
	
	N_signalevents[channel]= 1.0;
	//cout << channel << endl; 
	//cout <<"\n";  
	((hstPerCh*)HstPerCh[channel])->hNeventSignals->Fill(1.0);
	NPeak = FindPeaks(jentry,skipFstBin[isl],channel,((wave)Waves[channel]).nPt(),&((wave)Waves[channel]).Y[0],((wave)Waves[channel]).rms,&((wave)Waves[channel]).deriv[0],&((wave)Waves[channel]).sderiv[0],pkPos,pkHgt, timeRes, N_1, N_2, N_3, N_4, bslnTimeInterval, isChannel_1cm, isChannel_2cm, isChannel_1p5cm);
	//NPeak = FindPeaks(((wave)FltWaves[channel]).nPtInR(),&((wave)FltWaves[channel]).Y[skipFstBin],1.2e-3/*0.625*((wave)Waves[channel]).rms*/,6,3,pkPos,pkHgt);
	//npt, Float_t *amplitude, Float_t sig, Int_t nrise,Int_t checkUpTo, Int_t *pkPos, Float_t *pkHgt) {
	//cout<<"rms "<<((wave)Waves[channel]).rms<<endl;
	//0.625*rms= 2 sigma
	//cout <<"Is signal"<<endl;
	//cout <<"NPeak "<<NPeak <<endl;
	
	NPeak_clust = ClusterizationFindPeaks(cut_cluster_ns,nElectrons_per_cluster,jentry,skipFstBin[isl],channel,((wave)Waves[channel]).nPt(),((wave)Waves[channel]).rms,pkPos_clust,pkHgt_clust,pkPos,pkHgt,NPeak,timeRes, isChannel_1cm, isChannel_2cm, isChannel_1p5cm, scale_cut, isRuns_90_10, isRuns_85_15, isRuns_80_20);
	// cout <<"NPeak_clust" << NPeak_clust <<endl;
	if(NPeak_clust > 2 && ((float) X[pkPos_clust[0]] > 30. && (float) X[pkPos_clust[NPeak_clust - 1]]< 800.) && ((isChannel_1cm && (float) X[pkPos_clust[NPeak_clust-1]] - (float) X[pkPos_clust[0]]< 270.) || ((isChannel_2cm) && (float) X[pkPos_clust[NPeak_clust-1]] - (float) X[pkPos_clust[0]]< 580.) || ((isChannel_1p5cm) && (float) X[pkPos_clust[NPeak_clust-1]] - (float) X[pkPos_clust[0]]< 800.))){
	  
	  for(int m=0;m < NPeak_clust;m++){
	    ((hstPerCh*)HstPerCh[channel])->hNElectrons_per_cluster->Fill((float)nElectrons_per_cluster[m]);
	  }
	  for(int k=0;k<NPeak;k++){
	    if(k<NPeak-1){
	      ((hstPerCh*)HstPerCh[channel])->hTimeDifference->Fill((float)X[pkPos[k+1]] - (float)X[pkPos[k]]);
	    }
	  }
	  for (int ipk=0; ipk <NPeak; ++ipk){
	    ((hstPerCh*)HstPerCh[channel])->hHPeaks->Fill(pkHgt[ipk]);
	    ((hstPerCh*)HstPerCh[channel])->hHNPeaks->Fill(ipk+1,pkHgt[ipk]);
	  }
	  
	  ((hstPerCh*)HstPerCh[channel])->hNClusterFCluster->Fill((float)X[pkPos_clust[0]],NPeak_clust);
	  ((hstPerCh*)HstPerCh[channel])->hNPeakFPeak->Fill((float)X[pkPos[0]],NPeak);
	  
	  ((hstPerCh*)HstPerCh[channel])->hNPeaks->Fill((float)NPeak);
	  if(((isChannel_2cm) && (float) X[pkPos_clust[0]]< 200.) || isChannel_1cm || isChannel_1p5cm){
	    ((hstPerCh*)HstPerCh[channel])->hNPeaks_clust->Fill((float)NPeak_clust);
	  }
	  ((hstPerCh*)HstPerCh[channel])->hNPeaks_1->Fill(X[pkPos[NPeak-1]]); 
	  ((hstPerCh*)HstPerCh[channel])->hTFstPeaks->Fill(X[pkPos[0]]);
	  
	  for (int ipk=0; ipk < NPeak; ++ipk){
	    ((hstPerCh*)HstPerCh[channel])->hTPeaks->Fill(X[pkPos[ipk]]);
	    
	  }
	  ((hstPerCh*)HstPerCh[channel])->hRms->Fill(((wave)Waves[channel]).rms*1000);      
	  ((hstPerCh*)HstPerCh[channel])->hMaxVInR->Fill(((wave)Waves[channel]).maxInR);
	  ((hstPerCh*)HstPerCh[channel])->hBsl->Fill(((wave)Waves[channel]).bsln);
	  ((hstPerCh*)HstPerCh[channel])->hIntegN->Fill(((wave)Waves[channel]).charge_integInR);
	}
	//((hstPerCh*)HstPerCh[channel])->hMaxVNSmooth->Fill(((wave)tmpWSG_signal[channel]).nMax());
	//((hstPerCh*)HstPerCh[channel])->hInteg->Fill(((wave)Waves[channel]).integ);
	//((hstPerCh*)HstPerCh[channel])->hIntegInR->Fill(((wave)Waves[channel]).integInR);
	//((hstPerCh*)HstPerCh[channel])->hIntegNInRC1->Fill(((wave)Waves[channel]).nnIntegInR()/((float)NPeak));
	//istogrammi rms wave non filtrate      
	//((hstPerCh*)HstPerCh[channel])->hRmsOriginalW->Fill(((wave)Waves[channel]).rms);
	//((hstPerCh*)HstPerCh[channel])->hIntegNInRC2->Fill( ( ((wave)Waves[channel]).nnIntegInR()/((float)NPeak) )/scaleInt );   
	
	//if ((NPeak<10 || ((X[pkPos[0]+skipFstBin[isl]])< 20. || (X[pkPos[NPeak-1]+skipFstBin[isl]])> 300. || (X[pkPos[0]+skipFstBin[isl]])>350.)) && channel <= 10 && channel != 4 && channel != 0 && channel !=7) {
	//    
	//  cout << "Event 1cm tube having low NPeak || Wrong Peak position: " << jentry << " Ch: " << channel << " NPeaks: " <<NPeak<< " First Peak Position:"<< X[pkPos[0]+skipFstBin[isl]] <<" Last peak position:"<<X[pkPos[NPeak-1]+skipFstBin[isl]]<<endl; 
	//}
	//
	//if ((NPeak<20 || ((X[pkPos[0]+skipFstBin[isl]])< 20. || (X[pkPos[NPeak-1]+skipFstBin[isl]])> 650. || (X[pkPos[0]+skipFstBin[isl]])>600.)) && channel == 0 || channel == 4 || channel == 7 || channel == 11) {
	//    
	//  cout << "Event 1.5cm tube having low NPeak || Wrong Peak position: " << jentry << " Ch: " << channel << " NPeaks: " <<NPeak<< " First Peak Position:"<< X[pkPos[0]+skipFstBin[isl]] <<" Last peak position:"<<X[pkPos[NPeak-1]+skipFstBin[isl]]<<endl; 
	//}
	N_signalevents[channel]= 1.0;
	
      } //if for finding peaks
      
      
      
      //if(!isTrg && Waves[channel].max<=10*(Waves[channel].rms) && channel<=nMaxCh ){
      //  //cout <<"Is NOT a signal channel"<<endl;
      //  N_signalevents[channel]= 0.0; 
      //  //cout << channel << endl; 
      //  //cout <<"\n";   
      //  ((hstPerCh*)HstPerCh[channel])->hNeventSignals->Fill(0.0);//(Double_t) N_signalevents[channel]);;
      //}
      
      // cout << "\nAfter the threshold voltage is set, the channel hits are: \n"; 
      //for (int i = 0; i < N_signalevents.size(); i++) 
      //cout << "Event: " << jentry << " Ch: " << i << " Hits: " << N_signalevents[i] << "\n"; 
      //cout << "\n"; 
      
      bool savesignal_1 = true;
      bool uniqueCanvas = false;
      //Waves.clear();
      if (savesignal_1 && !isTrg && channel <=nMaxCh) { 
	//if (savesignal_1 && !isTrg && point.first == channel ) { //new graphs with arrows on the found peaks
	//Waves[channel].fillWave(point.second,dim,_gsample);
	signal->cd();
	
	if(uniqueCanvas){
	  if (firstEntering_filter){
	    tmpCvsignal_1.push_back( new TCanvas(Form("CvSignal_1_ev%d",jentry),Form("tmpSignal_1_ev%d",jentry)) );
	    tmpCvsignal_1.back()->Divide(3,4);
	    firstEntering_filter=false;
	  }
	  
	  tmpCvsignal_1.back()->cd(channel-nTriggerChannels+1);
	  tmpsignal_1.push_back( new TGraph ( dim, &X[0], &Waves[channel].Y[0]) );
	  tmpsignal_1.back()->GetXaxis()->SetTitle("time [ns]");
	  tmpsignal_1.back()->SetTitle(Form("tmpSignal_afterFlt_Ch%d_ev%d",channel,jentry));
	  tmpsignal_1.back()->SetTitle("");
	  tmpsignal_1.back()->GetYaxis()->SetTitleOffset(1.4);
	  tmpsignal_1.back()->GetYaxis()->SetTitle("Volt [V]");
	  //tmpsignal_1.back()->GetXaxis()->SetRangeUser(0.,400.);
	  tmpsignal_1.back()->GetYaxis()->SetRangeUser(-0.1,0.4);
	  tmpsignal_1.back()->SetMarkerStyle(21);
	  tmpsignal_1.back()->Draw("APL");
	  TLegend *leg= new TLegend(0.5,0.75,0.85,0.85); 
	  leg->AddEntry(tmpsignal_1.back(),"Waveform"); 
	  tmpsignal_1.push_back( new TGraph ( dim, &X[0], &Waves[channel].deriv[0]));
	  tmpsignal_1.back()->SetLineColor(kBlue);
	  tmpsignal_1.back()->Draw("Lsame");
	  leg->AddEntry(tmpsignal_1.back(),"First Derivative (Bin method)");
	  tmpsignal_1.push_back( new TGraph ( dim, &X[0], &Waves[channel].sderiv[0]));
	  tmpsignal_1.back()->SetLineColor(kRed);
	  tmpsignal_1.back()->Draw("Lsame");
	  
	  leg->AddEntry(tmpsignal_1.back(),"Second Derivative (Bin method)");
	  TLine *line = new TLine(X[0+30], (Waves[channel].rms)/(sqrt(2)), X[Waves[channel].nPt() - 350],(Waves[channel].rms)/(sqrt(2)));
	  line->SetLineColor(kOrange);
	  line->SetLineWidth(2);
	  line->Draw("same");
	  leg->AddEntry(line,"Sigma of the First Derivative (Bin method)");
	  TLine *line_1 = new TLine(X[0+30], (Waves[channel].rms)/(2), X[Waves[channel].nPt() - 350],(Waves[channel].rms)/2);
	  line_1->SetLineColor(kPink);
	  line_1->SetLineWidth(2);
	  line_1->Draw("same");
	  leg->AddEntry(line_1,"Sigma of the Second Derivative (Bin method)");
	  TLine *line_2 = new TLine(X[0+30], (Waves[channel].rms), X[Waves[channel].nPt() - 350],Waves[channel].rms);
	  line_2->SetLineColor(kViolet);
	  line_2->SetLineWidth(2);
	  line_2->Draw("same");
	  leg->AddEntry(line_2,"Rms (Bin method)");
	  leg->Draw("same");
	  counting_filter++;
	}
	//if(!uniqueCanvas && jentry<=200 && NPeak>0 && ((wave)Waves[channel]).nnIntegInR()>0.1 && ((channel <= 10 && channel != 4 && channel != 0 && channel !=7 && counter_1cm>=4) || (channel == 0 || channel == 4 || channel == 7 || channel == 11 && counter_1p5cm>=3))){
	if(!uniqueCanvas && NPeak>0 && NPeak_clust>0 && jentry<=1000 && Waves[channel].max > 0.05  && ((isChannel_1cm && Waves[channel].charge_integInR > 20 ) || ((isChannel_2cm || isChannel_1p5cm) && Waves[channel].charge_integInR > 60)) && (counter_top_first_driftTubes_line_size>=2 || counter_top_third_driftTubes_line_size>=2 || counter_top_second_driftTubes_line_size>=2)){
	  
	  bool Studies = true;
	  // Nov2021
	  //   int ChannelDiameter[13] = {-1,-1,-1,-1,10,15,20,20,25,25,20,25,40}; // Old test beam correspondance
	  //   float ChannelCellSize[13] = {-1.0,-1.0,-1.0,-1.0,1.0,1.0,1.0,1.0,1.0,1.0,2.0,2.0,2.0};
	  static int ChannelDiameter[16] = {20,20,15,15,25,20,20,15,20,25,25,15,-1,-1,-1,-1}; // Old test beam Nov 2022
	  static float ChannelCellSize[16] = {1.5,1.0,1.0,1.0,1.5,1.0,1.0,1.5,1.0,1.0,1.0,1.5,-1.0,-1.0,-1.0,-1.0}; // Old test beam Nov 2022
	  //else if(!uniqueCanvas && (X[pkPos[NPeak-1]< 20 || X[pkPos[NPeak-1]> 250 || X[pkPos[0]>350)){
	  tmpCvsignal_1.push_back( new TCanvas(Form("CvSignal_1_Ch%d_ev%d",channel,jentry),Form("tmpSignal_1_Ch%d_ev%d",channel,jentry)) );
	  tmpCvsignal_1.back()->cd();
	  tmpsignal_1.push_back( new TGraph ( dim, &X[0], &Waves[channel].Y[0]) );
	  tmpsignal_1.back()->GetXaxis()->SetTitle("time [ns]");
	  tmpsignal_1.back()->SetTitle(Form("Waveform signal Ch%d - Event %d - Sense Wire Diameter %d um - Cell Size %0.1f cm - Track Angle %.1f - %s - %0.1f GSa/s - Gas Mixture 80/20 %d - 90/10 %d - 85/15 %d",channel,jentry,ChannelDiameter[channel], ChannelCellSize[channel], alpha, out.Data(),_gsample, isRuns_80_20, isRuns_90_10, isRuns_85_15));
	  tmpsignal_1.back()->GetYaxis()->SetTitleOffset(1.4);
	  tmpsignal_1.back()->GetYaxis()->SetTitle("Volt [V]");
	  //tmpsignal_1.back()->GetXaxis()->SetRangeUser(0.,550);
	  //tmpsignal_1.back()->GetYaxis()->SetRangeUser(-0.05,0.250);
	  tmpsignal_1.back()->SetMarkerSize(1);
	  tmpsignal_1.back()->SetMarkerStyle(2);
	  tmpsignal_1.back()->Draw("APL");
	  TLegend *leg= new TLegend(0.5,0.65,0.85,0.85); 
	  leg->AddEntry(tmpsignal_1.back(),"Waveform"); 
	  tmpsignal_1.push_back( new TGraph ( dim, &X[0], &Waves[channel].deriv[0]) );
	  tmpsignal_1.back()->SetLineColor(kBlue);
	  if(Studies){
	    tmpsignal_1.back()->Draw("Lsame");
	    leg->AddEntry(tmpsignal_1.back(),"First Derivative (Bin method)");
	  }
	  tmpsignal_1.push_back( new TGraph ( dim, &X[0], &Waves[channel].sderiv[0]) );
	  tmpsignal_1.back()->SetLineColor(kRed);
	  if(Studies){
	    tmpsignal_1.back()->Draw("Lsame");
	    leg->AddEntry(tmpsignal_1.back(),"Second Derivative (Bin method)");
	  }
	  TLine *line = new TLine(X[0], (Waves[channel].rms)/(sqrt(2)), X[Waves[channel].nPt() - 1],(Waves[channel].rms)/(sqrt(2)));
	  line->SetLineColor(kOrange);
	  line->SetLineWidth(2);
	  if(Studies){
	    line->Draw("same");
	    leg->AddEntry(line,"Sigma of the First Derivative (Bin method)");
	  }
	  TLine *line_1 = new TLine(X[0], (Waves[channel].rms)/(2), X[Waves[channel].nPt() -1],(Waves[channel].rms)/(2));
	  line_1->SetLineColor(kPink);
	  line_1->SetLineWidth(2);
	  if(Studies){
	    line_1->Draw("same");
	    leg->AddEntry(line_1,"Sigma of the Second Derivative (Bin method)");
	  }
	  TLine *line_2 = new TLine(X[0], (Waves[channel].rms), X[Waves[channel].nPt() - 1],Waves[channel].rms);
	  line_2->SetLineColor(kViolet);
	  line_2->SetLineWidth(2);
	  if(Studies){
	    line_2->Draw("same");
	    leg->AddEntry(line_2,"Rms (Bin method)");
	  }
	  leg->Draw("same");
	  
	  
	  for (int ipk=0; ipk<NPeak; ipk++){
	    TMarker *tm = new  TMarker(X[pkPos[ipk]], pkHgt[ipk], 23);
	    tm->SetMarkerSize(1.2);
	    tm->SetMarkerColor(kRed);
	    tm->Draw("same");
	    if (ipk == NPeak -1){
	      leg->AddEntry(tm,"Electron Peaks from Primary Ionization Clusters");
	    }
	  }
	  
	  TPaveText *results_found = new TPaveText(0.5,0.48,0.8,0.63, "NDC");
	  results_found->SetTextSize(0.03);
	  results_found->SetTextColor(kRed);
	  results_found->SetFillColor(0);
	  results_found->SetTextAlign(12);
	  if(isChannel_1cm){
	    drift_size = 0.8;
	  }  
	  else if(isChannel_2cm){
	    drift_size = 1.8;
	  }
	  else if(isChannel_1p5cm){
	    drift_size = 1.2;
	  }
	  expected_cluster = cluster_per_cm_mip * drift_size * relativistic_rise * 1/cos_alpha;
	  // 1.5 mm per each wall -> 1.5 cm becomes 1.2 cm
	  //δ cluster/cm (M.I.P.) * drift tube size [cm] * 1.3 (relativisticrise) * 1.6 electrons/cluster * 1/cos(α)
	  expected_electrons = cluster_per_cm_mip * drift_size * relativistic_rise * cluster_population * 1/cos_alpha;
	  results_found->AddText(Form("Electrons found: %d",NPeak));
	  results_found->AddText(Form("Expected Electrons: %.1f",expected_electrons));
	  results_found->AddText(Form("Clusters found: %d",NPeak_clust));
	  results_found->AddText(Form("Expected Clusters: %.1f",expected_cluster));
	  results_found->AddText(Form("Track Angle: %.1f",alpha));
	  results_found->Draw("same");
	  
	  int electron_index = 0;
	  
	  for (int ipk=0; ipk < NPeak_clust; ipk++){
	    average_pkPos_clust[ipk] = 0;
	    average_pkHgt_clust[ipk] = 0;
	    for(int l=0; l < nElectrons_per_cluster[ipk];l++){
	      average_pkPos_clust[ipk] = average_pkPos_clust[ipk] + pkPos[electron_index];
	      //Difference in time of consecutive electrons belonging to the same cluster
	      //if(l<nElectrons_per_cluster[ipk]-1){
	      //((hstPerCh*)HstPerCh[channel])->hTimeDifference->Fill((float)X[pkPos[electron_index+1]+skipFstBin] - (float)X[pkPos[electron_index]+skipFstBin]);
	      //}
	      //cout << (float)X[pkPos[ipk+1]+skipFstBin] - (float)X[pkPos[ipk]+skipFstBin] <<"\n"<<endl;
	      electron_index = electron_index + 1 ;
	    }
	    average_pkPos_clust[ipk] = (int) (average_pkPos_clust[ipk]/nElectrons_per_cluster[ipk]);
	    average_pkHgt_clust[ipk] = (float) ((float) pkHgt_clust[ipk]/(float)nElectrons_per_cluster[ipk]);
	    //cout << "Average Position cluster "<< average_pkPos_clust[ipk] << " Position cluster "<< pkPos_clust[ipk] <<"\n";
	    //cout << "Average height cluster "<< average_pkHgt_clust[ipk] << " Height cluster "<< pkHgt_clust[ipk] <<"\n";
	    TMarker *tm_clust = new  TMarker(X[pkPos_clust[ipk]], pkHgt_clust[ipk], 23);
	    tm_clust->SetMarkerSize(1.2);
	    tm_clust->SetMarkerColor(kBlue);
	    tm_clust->Draw("same");
	    if (ipk == NPeak_clust -1){
	      leg->AddEntry(tm_clust,"Primary Ionization Clusters");
	    }
	  }
	  
	  for(int ipk=0; ipk<NPeak_clust;ipk++){
	    if(ipk<NPeak_clust-1){
	      ((hstPerCh*)HstPerCh[channel])->hTimeDifference_clust->Fill((float)X[pkPos_clust[ipk+1]] - (float)X[pkPos_clust[ipk]]);
	    }
	  }
	} 
	if(counting_filter==(nMaxCh-nTriggerChannels+1) && uniqueCanvas){
	  tmpCvsignal_1.back()->Write();
	}
	//else if(!uniqueCanvas && jentry<=200 && NPeak>0 && ((wave)Waves[channel]).nnIntegInR()>0.1 && ((channel <= 10 && channel != 4 && channel != 0 && channel !=7 && counter_1cm>=4) || (channel == 0 || channel == 4 || channel == 7 || channel == 11 && counter_1p5cm>=3))){
	else if(!uniqueCanvas && NPeak>0 && NPeak_clust>0 && jentry<=1000 && Waves[channel].max > 0.05  && ((isChannel_1cm && Waves[channel].charge_integInR > 20) || ((isChannel_2cm || isChannel_1p5cm) && Waves[channel].charge_integInR > 60)) && (counter_top_first_driftTubes_line_size>=2 || counter_top_third_driftTubes_line_size>=2 || counter_top_second_driftTubes_line_size>=2)){
	  //else if(!uniqueCanvas && NPeak>0 && NPeak_clust>0 && jentry<=200 && Waves[channel].max > 0.05 && ((isChannel_1cm && Waves[channel].charge_integInR > 20 ) || ((isChannel_2cm || isChannel_1p5cm) && Waves[channel].charge_integInR > 60)) && (counter_top_first_driftTubes_line_size>=3 || counter_top_third_driftTubes_line_size>=2 || counter_top_second_driftTubes_line_size>=2)){
	  tmpCvsignal_1.back()->Write();          
	}
	theFile->cd("/");
      } //if on representing the found peaks 
      
      
      if (!isTrg && channel<=nMaxCh && NPeak>0) { 
	
	if(false){
	  string nn_file = "";
	  nn_file = Form("nn_ch%d_alpha%.1f",channel,alpha);
	  nn_file = nn_file + out.Data() + ".txt";
	  ofstream myfile_nn (nn_file,ios::app);
	  if (myfile_nn.is_open())
	    {
	      drift_size = 0.8;
	      
	      //δ cluster/cm (M.I.P.) * drift tube size [cm] * 1.3 (relativisticrise) * 1.6 electrons/cluster * 1/cos(α)
	      expected_electrons = cluster_per_cm_mip * drift_size * relativistic_rise * cluster_population * 1/cos_alpha;
	      if(Waves[channel].max<=10*(Waves[channel].rms)&& channel<=nMaxCh ){
		expected_electrons = 0;
		myfile_nn << expected_electrons;
		myfile_nn << " ";
		//myfile_nn << NPeak;
		//myfile_nn << " ";
	      }
	      if(Waves[channel].max>10*(Waves[channel].rms)&& channel<=nMaxCh){
		myfile_nn << expected_electrons;
		myfile_nn << " ";
		//myfile_nn << NPeak;
		//myfile_nn << " ";
	      }
	      
	      if(channel <= 10 && channel != 4 && channel != 0 && channel !=7){
		//for(int i=0; i<(Waves[channel].nPt()-424);i++){
		for(int i=29; i<(Waves[channel].nPt()-64);i++){
		  if(i==29 || (i==Waves[channel].nPt()-65)){
		    myfile_nn<<(Waves[channel].Y[i]);
		    myfile_nn << " ";
		  }
		  if(i<Waves[channel].nPt()-65){
		    myfile_nn<<(Waves[channel].Y[i+1]-Waves[channel].Y[i]);
		    myfile_nn << " ";
		  }
		}
	      }
	      if(channel == 0 || channel == 4 || channel == 7 || channel == 11 ){
		for(int i=29; i<(Waves[channel].nPt()-64);i++){
		  //for(int i=0; i<(Waves[channel].nPt()-64);i++){
		  if(i==29 || (i==Waves[channel].nPt()-65)){
		    myfile_nn<<(Waves[channel].Y[i]);
		    myfile_nn << " ";
		  }
		  if(i<Waves[channel].nPt()-65){
		    myfile_nn<<(Waves[channel].Y[i+1]-Waves[channel].Y[i]);
		    myfile_nn << " ";
		  }
		}
	      }
	      
	      myfile_nn << "\n";
	      myfile_nn.close();
	    }
	  
	  else cout << "Unable to open file"; 
	}
	
	
      }
      
    } //getX742Data() loop
    
  } //entries loop
  
  
  theFile->cd();
  theFile->Write();
  theFile->Close();
  cout << "\n WELL DONE YOU HAVE FINISHED! \n"; 
  
  
  
}

