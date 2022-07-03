#define read_data_cxx
#include "read_data.h"
#include "funcUtils.h"
#include "FindPeak-algo.C"
#include "Clusterization.C"

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
//////////////Define Data Containers//////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
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
std::vector<TGraph *> tmpsignal_1;    	//grafici per il segnale finale dopo i filtri
std::vector<TCanvas *> tmpCvsignal_1; 	//canvas per il segnale finale dopo i filtri
std::vector<TGraph *> tmpNoise_total;  	//grafici per il noise totale
std::vector<TCanvas *>tmpCvNoise_total; //canvas per il noise totale
std::vector<float> X_negative; 			//vettore X per la parte negativa della waveform
std::vector<float> Y_negative; 			//vettore Y per la parte negativa della waveform
std::vector<TGraph *> neg_wavefbat; 	//grafico per waveform negative di battimento
std::vector< TGraphErrors *> tmpNegBatt; //grafico per le wave di battimento negative con errore
std::vector<TCanvas *> tmpCvNegBatt;
std::vector<TCanvas *> tmpCvMidNotch;
std::vector<TGraph *> tmpMid_Notch;
std::vector<float> N_signalevents;
std::vector<float>::iterator it; 



////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////MAIN LOOP function//////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void read_data::Loop(Char_t *output, Int_t MidEv,Int_t eventn,  Bool_t evalWaveCut,TString fOutName)
{
  tmax=_tmax;  //to fix compilation problem
  double timeRes;
#ifndef _OSC
  skipFstBin=0; //275;//set Bins to Skip
  skipLstBin=10;//375
  fbin_rms=30;
  invfbin=1.0/((float)fbin_rms);
#else
  skipFstBin=250; //275;//set Bins to Skip
  skipLstBin=300;//375
#endif
  
  /*********** creazione degli istogrammi***********/
  
  char *basename(char *path);
  cout << "Basename " << basename(output) << endl;
  Char_t outChar[500];
  TString out = basename(output);
  out.ReplaceAll(".root","");
#ifndef _OSC
  sprintf(outChar,"histosTB_%s.root",out.Data());   //Data sono i dati convertiti nei file.root
#else
  sprintf(outChar,"histosOSC_%s.root",out.Data());  //scrive il titolo della rootupla
#endif
  cout << "The output file is " << outChar << endl;
  TFile * theFile = new TFile(outChar,"RECREATE");  //crea un nuovo file.root e se ne esiste già uno lo sovrascrive
  
  std::map<int, hstPerCh *> HstPerCh;               //mappa per istogrammi sui segnali filtrati.
  int count =0;
  
  
  /***********creazione delle directory***********/
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
#ifndef _OSC
  trigCh.push_back(0);
  trigCh.push_back(1);
  trigCh.push_back(2);
  trigCh.push_back(3);
#else //for muon data
  trigCh.push_back(7);//osc1
  trigCh.push_back(8);//osc2
#endif
  
  
  int nSig=3; // numero di sigma
  std::map<int, int> nSelEv;
  nSelEv.clear();
  
  if (fChain == 0) return; //è il puntatore TTree
  
  Long64_t nentries = fChain->GetEntriesFast(); //Return the number of entries as of the last check.
  cout << "Number of entries in the tree for real data is= " << nentries << endl;
  //  Long64_t nbytes = 0, nb = 0;
  
  Int_t firstEv=0;                               //inizializzo il primo evento a zero
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
  
  
  for (Long64_t jentry=firstEv; jentry<lastEv;jentry++) {
    
    
    N_signalevents.assign(nMaxCh+1,0.0);
    cout << "New event is analyzed: "; 
    cout << jentry << "\n" << endl;
    cout << "Initial channel hits: \n"; 
    for (int i = 0; i < N_signalevents.size(); i++) 
      cout << "Ch: " << i << " Hits: " << N_signalevents[i] << "\n"; 
    cout << "\n";
    
    
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
    Waves_signal_1.clear();
    tmpWSG_signal.clear();
    SmtSGWaves.clear();
    
    //lettura dei file di dati .root
    
    if(jentry==firstEv) {
      // 1024 channel of ADC = 0.833*1024 ns // => every channel is 0.853 nano sec   
      timeRes = (_tmax*1.0e+9)/((float)dim);//in nano sec = 0.833s
      //timeRes = 1;    
      X.clear();
      cout<< "dim "<<dim<<" time "<<timeRes<<endl;
      for (int i = 0; i < dim; ++i) { X.push_back(timeRes *(i+1));  }
    }
    
    bool firstEntering=true;
    bool firstEntering_filter=true;
    int counting=0;
    int counting_filter=0;
    int nTriggerChannels=0;
    int counter=0;
    int counter_1cm=0;
    int counter_2cm=0;
    
    for (auto point : wd->getX742Data()) {
      
      for (int channel=0; channel<=nMaxCh; channel++){
	
	for(auto trgCh : trigCh) { if (channel==trgCh && firstEntering) {nTriggerChannels++;}}
	
      } //get the number of TriggerChannels
      
      int channel = point.first;
      //for (int channel=0; channel<=nMaxCh; channel++){
      
      bool isTrg=false;
      for(auto trgCh : trigCh) { if (channel==trgCh) { isTrg=true;}}
      //if (point.first == channel && !isTrg) {
      if (!isTrg && channel<=nMaxCh) { 
	
	Waves[channel].fillWave(point.second,dim);
	for(int i=0; i<Waves[channel].nPt();i++){
	  Waves[channel].Y[i]=Waves[channel].Y[i]-Waves[channel].bsln;
	}
	Waves_signal_1[channel].fillWave((Waves[channel].Y),Waves[channel].nPt());
	if (!isTrg && channel<=nMaxCh && Waves_signal_1[channel].max>10*(Waves_signal_1[channel].rms) && firstEntering && ((wave)Waves_signal_1[channel]).nnIntegInR()>0.1) {
	  if(channel<=9){
	    counter_1cm=counter_1cm+1;
	    
	  }
	  if(channel>=10){
	    counter_2cm=counter_2cm+1;
	  }
	  
	} //get the number of Signal 1cm and 2 cm channels
	
      }
    }
    
    for (auto point : wd->getX742Data()) {
      counter=counter+1;
      //cout << "Counter:"<< counter <<"\n"<<endl;
      cout << "Channel analysed from the map:" << point.first <<"\n"<<endl;
      
      for (int channel=0; channel<=nMaxCh; channel++){
	
	for(auto trgCh : trigCh) { if (channel==trgCh && firstEntering) {nTriggerChannels++;}}
	
      } //get the number of TriggerChannels
      
      int channel = point.first;
      //for (int channel=0; channel<=nMaxCh; channel++){
      
      bool isTrg=false;
      for(auto trgCh : trigCh) { if (channel==trgCh) { isTrg=true;}}
      //if (point.first == channel && !isTrg) {
      if (!isTrg && channel<=nMaxCh) { 
	if ( HstPerCh.find(channel)==HstPerCh.end() ) {
	  TDirectory *chDir = theFile->mkdir(Form("H-Ch%d_signal",channel));
	  chDir->cd();
	  HstPerCh.insert(make_pair(channel,new hstPerCh(channel)));
	  theFile->cd("/");
	}
	//cout << "Channel is= " << channel << endl;
	Waves[channel].fillWave(point.second,dim);
	/****-------------------------------TRASFORMATA DI FOURIER---------------------------------*****/
	
	FFT(Waves[channel],Wffts[channel]);//trasformata di Fourier
	double *realFltFFT = new double[Waves[channel].nPt()];
	double *imgFltFFT = new double[Waves[channel].nPt()];
	
	filterWaveBsl(Wffts[channel],realFltFFT,imgFltFFT);//filtro sulla baseline.                                        ///
	InverseFFT(realFltFFT,imgFltFFT,Waves[channel].nPt(),FltWaves[channel]);//trasformata inversa
	for(int i=0; i<Waves[channel].nPt();i++){
	  Waves[channel].Y[i]=Waves[channel].Y[i]-Waves[channel].bsln;
	}
	Waves_signal_1[channel].fillWave((Waves[channel].Y),Waves[channel].nPt());
	
	//Waves_signal_1[channel].fillWave(FltWaves[channel].Y,FltWaves[channel].nPt());
	cout << "Event: " << jentry << " Number of 1 cm channels hit: " << counter_1cm <<"\n";
	cout << "Event: " << jentry << " Number of 2 cm channels hit: " << counter_2cm <<"\n";
	
	
	bool saveEvents=false;
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
	  waveDir->cd();
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
	  waveDir->cd();
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
      Int_t NPeak;
      Int_t pkPos[250];
      Float_t pkHgt[250];
      Int_t nElectrons_per_cluster[250];
      Int_t NPeak_1;
      Int_t cut_cluster = 3;
      
      Int_t pkPos_1[250];
      Float_t pkHgt_1[250];
      //Clusterization variables//
      Int_t NPeak_clust=0;
      Int_t pkPos_clust[250];
      Int_t average_pkPos_clust[250];
      Float_t average_pkHgt_clust[250];
      Float_t pkHgt_clust[250];
      
      //maps for full waves
      vector<pair <int, int > > waveFull;
      waveFull.clear();
      pair<int,int> Full;
      
      NPeak=0;								
      NPeak_1=0;
      
      /* Adding SG filter smoothing
	 int m,k;
	 m=13; //number of bin interested by the SG smoothing 
	 k=3; //order of the polinomial used
	 std::vector<float> tmpWSG_signal_23=smoothSG(Waves[channel].Y,m,k);
	 tmpWSG_signal[channel].fillWave(tmpWSG_signal_23);
      */
      
      float scaleInt=1.0;	
      float sumamplitude=0.0;
      if(!isTrg && channel<=nMaxCh){
        for(int i=0;i<Waves_signal_1[channel].nPt();i++){
          sumamplitude =  sumamplitude + (Waves_signal_1[channel]).Y[i];
          //cout <<"Waveform values:"<< "i-eth: " <<i <<" channel :"<<channel<<(Waves_signal_1[channel]).Y[i] << "\n"<<endl;
        }   
        //if(jentry==1 && channel==9){ 
        //cout << sumamplitude << "\n"<<endl;//}
  	((hstPerCh*)HstPerCh[channel])->hSum->Fill(sumamplitude);	 
  	
      }	
      
      if(!isTrg && Waves_signal_1[channel].max>10*(Waves_signal_1[channel].rms) && channel<=nMaxCh && ((channel<=9 && counter_1cm>=4) || (channel>=10 && counter_2cm>=3))&&((wave)Waves_signal_1[channel]).nnIntegInR()>0.1){ //nPtInR == Y.size - first,lastBin; search peak when max amplitude > 5 mV
	
        cout <<"Is signal \n"<<endl;
        N_signalevents[channel]= 1.0;
        //cout << channel << endl; 
        //cout <<"\n";   
        ((hstPerCh*)HstPerCh[channel])->hNeventSignals->Fill(1.0);
	NPeak = FindPeaks(jentry,skipFstBin,channel,((wave)Waves_signal_1[channel]).nPtInR(),&((wave)Waves_signal_1[channel]).Y[skipFstBin],((wave)Waves_signal_1[channel]).rms,&((wave)Waves_signal_1[channel]).deriv[skipFstBin],&((wave)Waves_signal_1[channel]).sderiv[skipFstBin],pkPos,pkHgt);
	//NPeak = FindPeaks(((wave)FltWaves[channel]).nPtInR(),&((wave)FltWaves[channel]).Y[skipFstBin],1.2e-3/*0.625*((wave)Waves[channel]).rms*/,6,3,pkPos,pkHgt);
	//npt, Float_t *amplitude, Float_t sig, Int_t nrise,Int_t checkUpTo, Int_t *pkPos, Float_t *pkHgt) {
	//cout<<"rms "<<((wave)Waves[channel]).rms<<endl;
	//0.625*rms= 2 sigma
	NPeak_clust = NPeak;
	NPeak_clust = ClusterizationFindPeaks(cut_cluster,nElectrons_per_cluster,jentry,skipFstBin,channel,((wave)Waves_signal_1[channel]).nPtInR(),((wave)Waves_signal_1[channel]).rms,pkPos_clust,pkHgt_clust,pkPos,pkHgt,NPeak_clust);
	
	
	
	if ((NPeak>1 && channel<=9 && NPeak<100)|| (NPeak>9 && NPeak<150 && channel>=10 && channel<=12)) {
	  
	  for(int m=0;m < NPeak_clust;m++){
	    ((hstPerCh*)HstPerCh[channel])->hNElectrons_per_cluster->Fill((float)nElectrons_per_cluster[m]);
	  }
	  for(int k=0;k<NPeak;k++){
	    if(k<NPeak-1){
	      ((hstPerCh*)HstPerCh[channel])->hTimeDifference->Fill((float)X[pkPos[k+1]+skipFstBin] - (float)X[pkPos[k]+skipFstBin]);
	    }
	  }
	  ((hstPerCh*)HstPerCh[channel])->hRms->Fill(((wave)Waves_signal_1[channel]).rms*1000);	      
	  ((hstPerCh*)HstPerCh[channel])->hMaxV->Fill(((wave)Waves_signal_1[channel]).max);
	  ((hstPerCh*)HstPerCh[channel])->hMaxVN->Fill(((wave)Waves_signal_1[channel]).nMax());
	  ((hstPerCh*)HstPerCh[channel])->hMaxVInR->Fill(((wave)Waves_signal_1[channel]).maxInR);
	  ((hstPerCh*)HstPerCh[channel])->hMaxVNInR->Fill(((wave)Waves_signal_1[channel]).nMaxInR());
	  
	  ((hstPerCh*)HstPerCh[channel])->hMinV->Fill(((wave)Waves_signal_1[channel]).min);
	  ((hstPerCh*)HstPerCh[channel])->hMinVN->Fill(((wave)Waves_signal_1[channel]).nMin());
	  ((hstPerCh*)HstPerCh[channel])->hMinVInR->Fill(((wave)Waves_signal_1[channel]).minInR);
	  ((hstPerCh*)HstPerCh[channel])->hMinVNInR->Fill(((wave)Waves_signal_1[channel]).nMinInR());
	  
	  ((hstPerCh*)HstPerCh[channel])->hIntegNInRoriginalW->Fill(((wave)Waves_signal_1[channel]).nIntegInR());
	  
	  ((hstPerCh*)HstPerCh[channel])->hBsl->Fill(((wave)Waves[channel]).bsln);
	  //((hstPerCh*)HstPerCh[channel])->hMaxVNSmooth->Fill(((wave)tmpWSG_signal[channel]).nMax());
	  ((hstPerCh*)HstPerCh[channel])->hInteg->Fill(((wave)Waves_signal_1[channel]).integ);
	  ((hstPerCh*)HstPerCh[channel])->hIntegN->Fill(((wave)Waves_signal_1[channel]).nnInteg());
	  ((hstPerCh*)HstPerCh[channel])->hIntegNInR->Fill(((wave)Waves_signal_1[channel]).nnIntegInR());
	  ((hstPerCh*)HstPerCh[channel])->hIntegInR->Fill(((wave)Waves_signal_1[channel]).integInR);
	  ((hstPerCh*)HstPerCh[channel])->hIntegNInRC1->Fill(((wave)Waves_signal_1[channel]).nnIntegInR()/((float)NPeak));
	  
	  /*istogrammi sul minimo delle wave non filtrate*/
	  ((hstPerCh*)HstPerCh[channel])->hMinVoriginalW->Fill(((wave)Waves_signal_1[channel]).min);
	  ((hstPerCh*)HstPerCh[channel])->hMinVNoriginalW->Fill(((wave)Waves_signal_1[channel]).nMin());
	  ((hstPerCh*)HstPerCh[channel])->hMinVInRoriginalW->Fill(((wave)Waves_signal_1[channel]).minInR);
	  ((hstPerCh*)HstPerCh[channel])->hMinVNInRoriginalW->Fill(((wave)Waves_signal_1[channel]).nMinInR());
	  
	  /*istogrammi sul massimo delle wave non filtrate*/
	  ((hstPerCh*)HstPerCh[channel])->hMaxVoriginalW->Fill(((wave)Waves_signal_1[channel]).max);
	  ((hstPerCh*)HstPerCh[channel])->hMaxVNoriginalW->Fill(((wave)Waves_signal_1[channel]).nMax());
	  ((hstPerCh*)HstPerCh[channel])->hMaxVInRoriginalW->Fill(((wave)Waves_signal_1[channel]).maxInR);
	  ((hstPerCh*)HstPerCh[channel])->hMaxVNInRoriginalW->Fill(((wave)Waves_signal_1[channel]).nMaxInR());
	  
	  //istogrammi rms wave non filtrate	      
	  ((hstPerCh*)HstPerCh[channel])->hRmsOriginalW->Fill(((wave)Waves_signal_1[channel]).rms);
	  
	  for (int ipk=0; ipk <NPeak; ++ipk){
	    ((hstPerCh*)HstPerCh[channel])->hHPeaks->Fill(pkHgt[ipk]);
	    ((hstPerCh*)HstPerCh[channel])->hHNPeaks->Fill(ipk+1,pkHgt[ipk]);
	  }
	  
	  ((hstPerCh*)HstPerCh[channel])->hNPeaks->Fill((float)NPeak);
	  ((hstPerCh*)HstPerCh[channel])->hNPeaks_clust->Fill((float)NPeak_clust);
	  ((hstPerCh*)HstPerCh[channel])->hNPeaks_1->Fill(X[pkPos[NPeak-1]+skipFstBin]); 
	  ((hstPerCh*)HstPerCh[channel])->hTFstPeaks->Fill(X[pkPos[0]+skipFstBin]);
	  
	  for (int ipk=0; ipk < NPeak; ++ipk){
	    ((hstPerCh*)HstPerCh[channel])->hTPeaks->Fill(X[pkPos[ipk]+skipFstBin]);
	    
	  }
	  
	  float minDist=1e+20;
	  float tmpDist;
	  int clstFreq=0;
	  int numOfPeaks=0;
	  
	  ((hstPerCh*)HstPerCh[channel])->hIntegNInRC2->Fill( ( ((wave)Waves_signal_1[channel]).nnIntegInR()/((float)NPeak) )/scaleInt );   
	  
	}
	if ((NPeak<10 || ((X[pkPos[0]+skipFstBin])< 20. || (X[pkPos[NPeak-1]+skipFstBin])> 300. || (X[pkPos[0]+skipFstBin])>350.)) && channel<=9) {
	  
	  cout << "Event 1cm tube having low NPeak || Wrong Peak position: " << jentry << " Ch: " << channel << " NPeaks: " <<NPeak<< " First Peak Position:"<< X[pkPos[0]+skipFstBin] <<" Last peak position:"<<X[pkPos[NPeak-1]+skipFstBin]<<"\n"; 
	}
	
	if ((NPeak<20 || ((X[pkPos[0]+skipFstBin])< 20. || (X[pkPos[NPeak-1]+skipFstBin])> 650. || (X[pkPos[0]+skipFstBin])>600.)) && channel>=10 && channel<=12) {
	  
	  cout << "Event 2cm tube having low NPeak || Wrong Peak position: " << jentry << " Ch: " << channel << " NPeaks: " <<NPeak<< " First Peak Position:"<< X[pkPos[0]+skipFstBin] <<" Last peak position:"<<X[pkPos[NPeak-1]+skipFstBin]<<"\n"; 
	}
        N_signalevents[channel]= 1.0;
	
	
	
	
      } //if for finding peaks
      
      
      
      if(!isTrg && Waves_signal_1[channel].max<=10*(Waves_signal_1[channel].rms) && channel<=nMaxCh && ((channel<=9 && counter_1cm<4) || (channel>=10 && counter_2cm<2))){
	cout <<"Is NOT signal \n"<<endl;
	N_signalevents[channel]= 0.0; 
	//cout << channel << endl; 
	cout <<"\n";   
	((hstPerCh*)HstPerCh[channel])->hNeventSignals->Fill(0.0);//(Double_t) N_signalevents[channel]);;
      }
      
      // cout << "\nAfter the threshold voltage is set, the channel hits are: \n"; 
      //for (int i = 0; i < N_signalevents.size(); i++) 
      //cout << "Event: " << jentry << " Ch: " << i << " Hits: " << N_signalevents[i] << "\n"; 
      cout << "\n"; 
      
      bool savesignal_1 = true;
      bool uniqueCanvas = false;
      Waves.clear();
      if (savesignal_1 && !isTrg && channel <=nMaxCh) { 
	//if (savesignal_1 && !isTrg && point.first == channel ) { //new graphs with arrows on the found peaks
	Waves[channel].fillWave(point.second,dim);
	signal->cd();
	
	if(uniqueCanvas){
	  if (firstEntering_filter){
	    tmpCvsignal_1.push_back( new TCanvas(Form("CvSignal_1_ev%d",jentry),Form("tmpSignal_1_ev%d",jentry)) );
	    tmpCvsignal_1.back()->Divide(3,4);
	    firstEntering_filter=false;
	  }
	  
	  tmpCvsignal_1.back()->cd(channel-nTriggerChannels+1);
	  tmpsignal_1.push_back( new TGraph ( dim, &X[0], &Waves_signal_1[channel].Y[0]) );
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
	  tmpsignal_1.push_back( new TGraph ( dim, &X[0], &Waves_signal_1[channel].deriv[0]));
	  tmpsignal_1.back()->SetLineColor(kBlue);
	  tmpsignal_1.back()->Draw("Lsame");
	  leg->AddEntry(tmpsignal_1.back(),"First Derivative (Bin method)");
	  tmpsignal_1.push_back( new TGraph ( dim, &X[0], &Waves_signal_1[channel].sderiv[0]));
	  tmpsignal_1.back()->SetLineColor(kRed);
	  tmpsignal_1.back()->Draw("Lsame");
	  
	  leg->AddEntry(tmpsignal_1.back(),"Second Derivative (Bin method)");
	  TLine *line = new TLine(X[0+30], (Waves_signal_1[channel].rms)/(sqrt(2)), X[Waves_signal_1[channel].nPt() - 350],(Waves_signal_1[channel].rms)/(sqrt(2)));
	  line->SetLineColor(kOrange);
	  line->SetLineWidth(2);
	  line->Draw("same");
	  leg->AddEntry(line,"Sigma of the First Derivative (Bin method)");
	  TLine *line_1 = new TLine(X[0+30], (Waves_signal_1[channel].rms)/(2), X[Waves_signal_1[channel].nPt() - 350],(Waves_signal_1[channel].rms)/2);
	  line_1->SetLineColor(kPink);
	  line_1->SetLineWidth(2);
	  line_1->Draw("same");
	  leg->AddEntry(line_1,"Sigma of the Second Derivative (Bin method)");
	  TLine *line_2 = new TLine(X[0+30], (Waves_signal_1[channel].rms), X[Waves_signal_1[channel].nPt() - 350],Waves_signal_1[channel].rms);
	  line_2->SetLineColor(kViolet);
	  line_2->SetLineWidth(2);
	  line_2->Draw("same");
	  leg->AddEntry(line_2,"Rms (Bin method)");
	  leg->Draw("same");
	  counting_filter++;
	}
	if(!uniqueCanvas && jentry<=200 && NPeak>0 && ((wave)Waves_signal_1[channel]).nnIntegInR()>0.1 && ((channel<=9 && counter_1cm>=4) || (channel>=10 && counter_2cm>=3))){
	  
	  for(int i=0; i<Waves[channel].nPt();i++){
	    Waves_signal_1[channel].deriv[i]=(Waves_signal_1[channel].deriv[i])*10.;
	    Waves_signal_1[channel].sderiv[i]=(Waves_signal_1[channel].sderiv[i])*10.;
	  }
	  //else if(!uniqueCanvas && (X[pkPos[NPeak-1]+skipFstBin]< 20 || X[pkPos[NPeak-1]+skipFstBin]> 250 || X[pkPos[0]+skipFstBin]>350)){
	  tmpCvsignal_1.push_back( new TCanvas(Form("CvSignal_1_Ch%d_ev%d",channel,jentry),Form("tmpSignal_1_Ch%d_ev%d",channel,jentry)) );
	  tmpCvsignal_1.back()->cd();
	  tmpsignal_1.push_back( new TGraph ( dim, &X[0], &Waves_signal_1[channel].Y[0]) );
	  tmpsignal_1.back()->GetXaxis()->SetTitle("time [ns]");
	  tmpsignal_1.back()->SetTitle(Form("tmpSignal_afterFlt_Ch%d_ev%d_%s",channel,jentry,out.Data()));
	  tmpsignal_1.back()->GetYaxis()->SetTitleOffset(1.4);
	  tmpsignal_1.back()->GetYaxis()->SetTitle("Volt [V]");
	  //tmpsignal_1.back()->GetXaxis()->SetRangeUser(0.,400);
	  tmpsignal_1.back()->GetYaxis()->SetRangeUser(-0.1,0.250);
	  tmpsignal_1.back()->SetMarkerSize(1);
	  tmpsignal_1.back()->SetMarkerStyle(2);
	  tmpsignal_1.back()->Draw("APL");
	  TLegend *leg= new TLegend(0.5,0.65,0.85,0.85); 
	  leg->AddEntry(tmpsignal_1.back(),"Waveform"); 
	  tmpsignal_1.push_back( new TGraph ( dim, &X[0], &Waves_signal_1[channel].deriv[0]) );
	  tmpsignal_1.back()->SetLineColor(kBlue);
	  //tmpsignal_1.back()->Draw("Lsame");
	  leg->AddEntry(tmpsignal_1.back(),"First Derivative (Bin method) x10");
	  tmpsignal_1.push_back( new TGraph ( dim, &X[0], &Waves_signal_1[channel].sderiv[0]) );
	  tmpsignal_1.back()->SetLineColor(kRed);
	  //tmpsignal_1.back()->Draw("Lsame");
	  leg->AddEntry(tmpsignal_1.back(),"Second Derivative (Bin method) x10");
	  TLine *line = new TLine(X[0], 20*(Waves_signal_1[channel].rms)/(sqrt(2)), X[Waves_signal_1[channel].nPt() - 1],20*(Waves_signal_1[channel].rms)/(sqrt(2)));
	  line->SetLineColor(kOrange);
	  line->SetLineWidth(2);
	  //line->Draw("same");
	  leg->AddEntry(line,"Sigma of the First Derivative (Bin method) x20");
	  TLine *line_1 = new TLine(X[0], 20*(Waves_signal_1[channel].rms)/(2), X[Waves_signal_1[channel].nPt() -1],20*(Waves_signal_1[channel].rms)/(2));
	  line_1->SetLineColor(kPink);
	  line_1->SetLineWidth(2);
	  //line_1->Draw("same");
	  leg->AddEntry(line_1,"Sigma of the Second Derivative (Bin method) x20");
	  TLine *line_2 = new TLine(X[0], 20*(Waves_signal_1[channel].rms), X[Waves_signal_1[channel].nPt() - 1],20*Waves_signal_1[channel].rms);
	  line_2->SetLineColor(kViolet);
	  line_2->SetLineWidth(2);
	  //line_2->Draw("same");
	  leg->AddEntry(line_2,"Rms (Bin method) x20");
	  //leg->Draw("same");
        
	}
	for (int ipk=0; ipk<NPeak; ipk++){
	  TMarker *tm = new  TMarker(X[pkPos[ipk]+skipFstBin]+0.5*0.833333, pkHgt[ipk], 23);
	  tm->SetMarkerSize(1.5);
	  tm->SetMarkerColor(2);
	  tm->Draw("same");
	}
	
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
	  
	  
	  cout << "Average Position cluster "<< average_pkPos_clust[ipk] << " Position cluster "<< pkPos_clust[ipk] <<"\n";
	  cout << "Average height cluster "<< average_pkHgt_clust[ipk] << " Height cluster "<< pkHgt_clust[ipk] <<"\n";
	  
	  TMarker *tm_clust = new  TMarker(X[pkPos_clust[ipk]+skipFstBin]+0.5*0.833333, 0., 23);
	  tm_clust->SetMarkerSize(1.5);
	  tm_clust->SetMarkerColor(kBlue);
	  tm_clust->Draw("same");
	}
	for(int ipk=0; ipk<NPeak_clust;ipk++){
	  if(ipk<NPeak_clust-1){
	    ((hstPerCh*)HstPerCh[channel])->hTimeDifference_clust->Fill((float)X[pkPos_clust[ipk+1]+skipFstBin] - (float)X[pkPos_clust[ipk]+skipFstBin]);
	  }
	}
	
	if(counting_filter==(nMaxCh-nTriggerChannels+1) && uniqueCanvas){					
	  tmpCvsignal_1.back()->Write();
	}
        else if(!uniqueCanvas && jentry<=200 && NPeak>0 && ((wave)Waves_signal_1[channel]).nnIntegInR()>0.1 && ((channel<=9 && counter_1cm>=4) || (channel>=10 && counter_2cm>=3))){
	  tmpCvsignal_1.back()->Write();          
        }
	theFile->cd("/");
      } //if on representing the found peaks 
      
      
      //} //channel loop
      
    } //getX742Data() loop
    
  }	 //entries loop
  
  
  /************************************Isto su rms e max*************************************/
  std::vector<float> rms_m;
  std::vector<float> max_m;
  std::vector<float> maxCut_m;
  std::vector<float> x_channel;
  for (int channel=0; channel<=nMaxCh; channel++){
    bool isTrg=false;
    for(auto trgCh : trigCh) { if (channel==trgCh) { isTrg=true; } }
    if (!isTrg) {
      if ( HstPerCh.find(channel)!=HstPerCh.end() ) {
        x_channel.push_back(channel);
        /*rms_m[channel] =*/  rms_m.push_back( ((hstPerCh*)HstPerCh[channel])->hRms->GetMean() );
	max_m.push_back( ((hstPerCh*)HstPerCh[channel])->hMaxVNInR->GetMean() );
	((hstPerCh*)HstPerCh[channel])->hMaxVNInR->GetXaxis()->SetRangeUser(0.008,0.1);
	/*max_m [channel] =*/ maxCut_m.push_back( ((hstPerCh*)HstPerCh[channel])->hMaxVNInR->GetMean() );
      }
    }
  }  
  
  TCanvas *cvStatSum = new TCanvas("cvStatSum","ChStatSum");
  cvStatSum->Divide(1,3);
  TGraph* rms_d= new TGraph(x_channel.size(),&x_channel[0],&rms_m[0]);
  TGraph* max_d= new TGraph(x_channel.size(),&x_channel[0],&max_m[0]);
  TGraph* maxCut_d= new TGraph(x_channel.size(),&x_channel[0],&maxCut_m[0]);
  cvStatSum->cd(1);
  rms_d->SetMarkerStyle(8);
  rms_d->SetTitle("rms per channel");
  rms_d->GetXaxis()->SetTitle("channel");
  rms_d->GetXaxis()->SetTitleSize(0.06);
  rms_d->GetXaxis()->SetTitleOffset(0.6);
  rms_d->GetYaxis()->SetTitle("rms");
  rms_d->GetYaxis()->SetTitleSize(0.06);
  rms_d->GetYaxis()->SetTitleOffset(0.6);
  rms_d->Draw("AP");
  cvStatSum->cd(2);
  max_d->SetMarkerStyle(8);
  max_d->SetTitle("max per channel");
  max_d->GetXaxis()->SetTitle("channel");
  max_d->GetXaxis()->SetTitleOffset(0.5);
  max_d->GetXaxis()->SetTitleSize(0.06);
  max_d->GetYaxis()->SetTitle("max");
  max_d->GetYaxis()->SetTitleOffset(0.6);
  max_d->GetYaxis()->SetTitleSize(0.06);
  max_d->Draw("AP");
  cvStatSum->cd(3);
  maxCut_d->SetMarkerStyle(8);
  maxCut_d->SetTitle("max (>0.01) per channel");
  maxCut_d->GetXaxis()->SetTitle("channel");
  maxCut_d->GetXaxis()->SetTitleOffset(0.5);
  maxCut_d->GetXaxis()->SetTitleSize(0.06);
  maxCut_d->GetYaxis()->SetTitle("max");
  maxCut_d->GetYaxis()->SetTitleOffset(0.6);
  maxCut_d->GetYaxis()->SetTitleSize(0.06);
  maxCut_d->Draw("AP");
  cvStatSum->Write();
  
  TH1F *hSnr = new TH1F("hSnr","peak val over noise val;ratio;Entries",500,0,50);
  TH1F *hSnr1 = new TH1F("hSnr1","peak val over noise val;ratio;Entries",500,0,50);
  for (int ipt=0; ipt<x_channel.size(); ++ipt) {
    hSnr->Fill(max_m[ipt]/rms_m[ipt]);
    hSnr1->Fill(maxCut_m[ipt]/rms_m[ipt]);
  }
  
  ///////mi dice quanti eventi ho per ogni canale.
  bool prtNSelEv=true;
  if (prtNSelEv) {
    for (int channel=0; channel<=nMaxCh; channel++){
      if (nSelEv.find(channel)!=nSelEv.end()) {      
      }
    }
  }
  
  
  theFile->cd();
  theFile->Write();
  theFile->Close();
  cout << "\n WELL DONE YOU HAVE FINISHED! \n"; 
}

