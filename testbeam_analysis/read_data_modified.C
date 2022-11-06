#define read_data_cxx
#include "read_data.h"
#include "funcUtils.h"
#include "FindPeak-algo.C"

#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TGraph.h>
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
  skipFstBin=1; //275;//set Bins to Skip
  skipLstBin=10;//375
  fbin_rms=50;
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
    cout << "The vector contains: "; 
    cout << N_signalevents.size();
    for (int i = 0; i < N_signalevents.size(); i++) 
        cout << N_signalevents[i] << " "; 
        cout << "\n";
        

    Long64_t ientry = LoadTree(jentry); //The function finds the corresponding Tree and returns the entry number in this tree.  
    if (ientry < 0) break;
    //nb = fChain->GetEntry(jentry);   nbytes += nb;
    fChain->GetEntry(jentry);
    WvCont Waves;
    WvCont tmpWSG_signal;
    WvCont SmtSGWaves;
    Waves.clear();
    tmpWSG_signal.clear();
    SmtSGWaves.clear();

    //lettura dei file di dati .root

    if(jentry==firstEv) {
      // 1024 channel of ADC = 1 micro sec => every channel is 0.9766 nano sec   
      timeRes = (_tmax*1.0e+9)/((float)dim);//in nano sec      
      X.clear();
      cout<< "dim "<<dim<<" time "<<timeRes<<endl;
      for (int i = 0; i < dim; ++i) { X.push_back(timeRes *(i+1));  }
    }
    
    bool firstEntering=true;
    bool firstEntering_filter=true;
    int counting=0;
    int counting_filter=0;
    int nTriggerChannels=0;
    for (auto point : wd->getX742Data()) {
    
    for (int channel=0; channel<=nMaxCh; channel++){
        
	  for(auto trgCh : trigCh) { if (channel==trgCh && firstEntering) {nTriggerChannels++;}}} //get the number of TriggerChannels
    
    for (int channel=0; channel<=nMaxCh; channel++){
        
	bool isTrg=false;
	for(auto trgCh : trigCh) { if (channel==trgCh) { isTrg=true;}}
	if (point.first == channel && !isTrg) {
          
	  if ( HstPerCh.find(channel)==HstPerCh.end() ) {
	    TDirectory *chDir = theFile->mkdir(Form("H-Ch%d_signal",channel));
	    chDir->cd();
	    HstPerCh.insert(make_pair(channel,new hstPerCh(channel)));
	    theFile->cd("/");
	  }
	  //cout << "Channel is= " << channel << endl;
	  Waves[channel].fillWave(point.second,dim);
	  
	  bool saveWave=false; //non salvo, quindi non vediamo l'out di questo pezzo di codice.
	  //						bool saveWave=true;
	  bool saveEvents=true;
	  
	  
	  if (saveWave) {
	    waveDir->cd();
	    tmpCv.push_back(new TCanvas(Form("Cv-Ch%d_ev%d",channel,jentry),Form("tmpWave-Ch%d_ev%d",channel,jentry)));
	    tmpCv.back()->cd();
	    tmpWaves.push_back(new TGraph (dim, &X[0], &Waves[channel].Y[0]));
	    tmpWaves.back()->GetXaxis()->SetTitle("time [ns]");
	    tmpWaves.back()->SetTitle(Form("tmpWave-Ch%d_ev%d",channel,jentry));
	    tmpWaves.back()->GetYaxis()->SetTitleOffset(1.4);
	    tmpWaves.back()->GetYaxis()->SetTitle("Voltage [V]");
	    tmpWaves.back()->Draw("AL");
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
	    tmpWaves.back()->Draw("AL");
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
    Int_t pkPos[100];
    Float_t pkHgt[100];
    Int_t NPeak_1;
    Int_t pkPos_1[100];
    Float_t pkHgt_1[100];
    
    //maps for full waves
    vector<pair <int, int > > waveFull;
    waveFull.clear();
    pair<int,int> Full;
    
    NPeak=0;								
 	  NPeak_1=0;
    
    int m,k;
	  m=13; //number of bin interested by the SG smoothing 
	  k=3; //order of the polinomial used
	  std::vector<float> tmpWSG_signal_23=smoothSG(Waves[channel].Y,m,k);
	  tmpWSG_signal[channel].fillWave(tmpWSG_signal_23);

   
	  float scaleInt=1.0;		
	    if(!isTrg&&Waves[channel].max>0.005){ //nPtInR == Y.size - first,lastBin; search peak when max amplitude > 5 mV
        cout <<"Is signal \n"<<endl;
            N_signalevents[channel]= 1.0;
            cout << channel << endl; 
           cout <<"\n";   
           ((hstPerCh*)HstPerCh[channel])->hNeventSignals->Fill(1.0);
        NPeak = FindPeaks(((wave)Waves[channel]).nPtInR(),&((wave)Waves[channel]).Y[skipFstBin],1.2e-3/*0.625*((wave)Waves[channel]).rms*/,6,3,pkPos,pkHgt);
	      //npt, Float_t *amplitude, Float_t sig, Int_t nrise,Int_t checkUpTo, Int_t *pkPos, Float_t *pkHgt) {
	      //cout<<"rms "<<((wave)Waves[channel]).rms<<endl;
	      //0.625*rms= 2 sigma
	      float sig=15e-3;//*2.5*((wave)Waves_signal_1[channel]).rms;//3.0*/*1.414**//*1.0e-3;/*/((wave)Waves_signal_1[channel]).rms;				
	      float sig2=1.*sig;//1.0*sig;
	      float sig3=1.414*sig;
	      float meanValueLocal;

	      for (int ipk=0; ipk <NPeak; ++ipk) {
                bool skip=false;
                bool doCheck=false;
                int iPkBin=skipFstBin+pkPos[ipk]; //Bin related to the bin found
                if ( ((wave)Waves[channel]).nnAt(skipFstBin+pkPos[ipk]) > sig2 ) {
                  if ( ( ((wave)Waves[channel]).nMaxInR()>0.01 && iPkBin>((wave)Waves[channel]).maxInRPos )) {
                    doCheck=true;
		  }
		} else {
		  if ( ((wave)Waves[channel]).nnAt(skipFstBin+pkPos[ipk]) < 3*sig ) {
                    skip=true;
                  }else {
                    doCheck=true;
		  }
		}

		if (doCheck) {
                  for (int ick=1; ick <2; ++ick) {
                    if ( fabs(((wave)Waves[channel]).Y[iPkBin]-((wave)Waves[channel]).Y[iPkBin-ick])>sig3
		      || fabs(((wave)Waves[channel]).Y[iPkBin]-((wave)Waves[channel]).Y[iPkBin+ick])>sig3 ) {
		      skip=true;
		      break;
		    }
		  }
		}
		if (!skip) {
		  pkPos_1[NPeak_1]=pkPos[ipk];
		  pkHgt_1[NPeak_1]=pkHgt[ipk];
		  ++NPeak_1;
		}
              } 

	      ((hstPerCh*)HstPerCh[channel])->hBsl->Fill(((wave)Waves[channel]).bsln);
        ((hstPerCh*)HstPerCh[channel])->hMaxVNSmooth->Fill(((wave)tmpWSG_signal[channel]).nMax());
	      ((hstPerCh*)HstPerCh[channel])->hInteg->Fill(((wave)Waves[channel]).integ);
	      ((hstPerCh*)HstPerCh[channel])->hIntegN->Fill(((wave)Waves[channel]).nnInteg());
	      ((hstPerCh*)HstPerCh[channel])->hIntegInR->Fill(((wave)Waves[channel]).integInR);
	      ((hstPerCh*)HstPerCh[channel])->hIntegNInR->Fill(((wave)Waves[channel]).nnIntegInR());
	      
	      if (NPeak_1>2/*0*/) {
                ((hstPerCh*)HstPerCh[channel])->hNPeaks->Fill(NPeak_1);
		for (int ipk=0; ipk <NPeak_1; ++ipk){
		  ((hstPerCh*)HstPerCh[channel])->hHPeaks->Fill(pkHgt_1[ipk]);
		  ((hstPerCh*)HstPerCh[channel])->hHNPeaks->Fill(ipk+1,pkHgt_1[ipk]);
		}

		((hstPerCh*)HstPerCh[channel])->hIntegNInRC1->Fill(((wave)Waves[channel]).nnIntegInR()/((float)NPeak_1));
		if (((wave)Waves[channel]).nnIntegInR()>0.1/*0.2*/) {
		  ((hstPerCh*)HstPerCh[channel])->hNPeaks_1->Fill(NPeak_1);
		  ((hstPerCh*)HstPerCh[channel])->hTFstPeaks->Fill(X[pkPos_1[0]+skipFstBin]);
		  for (int ipk=0; ipk <NPeak_1; ++ipk){
		    ((hstPerCh*)HstPerCh[channel])->hTPeaks->Fill(X[pkPos_1[ipk]+skipFstBin]);
		  }
												
		  float minDist=1e+20;
		  float tmpDist;
		  int clstFreq=0;
		  int numOfPeaks=0;

		  ((hstPerCh*)HstPerCh[channel])->hIntegNInRC2->Fill( ( ((wave)Waves[channel]).nnIntegInR()/((float)NPeak_1) )/scaleInt );   
                }
              }
        
	      ((hstPerCh*)HstPerCh[channel])->hRms->Fill(((wave)Waves[channel]).rms);	      
	      ((hstPerCh*)HstPerCh[channel])->hMaxV->Fill(((wave)Waves[channel]).max);
	      ((hstPerCh*)HstPerCh[channel])->hMaxVN->Fill(((wave)Waves[channel]).nMax());
	      ((hstPerCh*)HstPerCh[channel])->hMaxVInR->Fill(((wave)Waves[channel]).maxInR);
	      ((hstPerCh*)HstPerCh[channel])->hMaxVNInR->Fill(((wave)Waves[channel]).nMaxInR());
	      
	      ((hstPerCh*)HstPerCh[channel])->hMinV->Fill(((wave)Waves[channel]).min);
	      ((hstPerCh*)HstPerCh[channel])->hMinVN->Fill(((wave)Waves[channel]).nMin());
	      ((hstPerCh*)HstPerCh[channel])->hMinVInR->Fill(((wave)Waves[channel]).minInR);
	      ((hstPerCh*)HstPerCh[channel])->hMinVNInR->Fill(((wave)Waves[channel]).nMinInR());
	      
	      ((hstPerCh*)HstPerCh[channel])->hIntegNInRoriginalW->Fill(((wave)Waves[channel]).nIntegInR());
	      
	      /*istogrammi sul minimo delle wave non filtrate*/
	      ((hstPerCh*)HstPerCh[channel])->hMinVoriginalW->Fill(((wave)Waves[channel]).min);
	      ((hstPerCh*)HstPerCh[channel])->hMinVNoriginalW->Fill(((wave)Waves[channel]).nMin());
	      ((hstPerCh*)HstPerCh[channel])->hMinVInRoriginalW->Fill(((wave)Waves[channel]).minInR);
	      ((hstPerCh*)HstPerCh[channel])->hMinVNInRoriginalW->Fill(((wave)Waves[channel]).nMinInR());
	      
	      /*istogrammi sul massimo delle wave non filtrate*/
	      ((hstPerCh*)HstPerCh[channel])->hMaxVoriginalW->Fill(((wave)Waves[channel]).max);
	      ((hstPerCh*)HstPerCh[channel])->hMaxVNoriginalW->Fill(((wave)Waves[channel]).nMax());
	      ((hstPerCh*)HstPerCh[channel])->hMaxVInRoriginalW->Fill(((wave)Waves[channel]).maxInR);
	      ((hstPerCh*)HstPerCh[channel])->hMaxVNInRoriginalW->Fill(((wave)Waves[channel]).nMaxInR());

	      //istogrammi rms wave non filtrate	      
	      ((hstPerCh*)HstPerCh[channel])->hRmsOriginalW->Fill(((wave)Waves[channel]).rms);

            } //if for finding peaks

        //  cout << "The vector contains: "; 
        //for (int i = 0; i < N_signalevents.size(); i++) 
        //cout << N_signalevents[i] << " "; 

       // cout <<"\n";  
        
          

          if(!isTrg&&Waves[channel].max<=0.005){
            cout <<"If not signal \n"<<endl;
            N_signalevents[channel]= 0.0; 
            cout << channel << endl; 
           cout <<"\n";   
           //((hstPerCh*)HstPerCh[channel])->hNeventSignals->Fill(0.0);//(Double_t) N_signalevents[channel]);;
          }

   //   cout << "The vector contains second: ";
   //   cout << N_signalevents.size(); 
   //    for (int i = 0; i < N_signalevents.size(); i++) 
   //     cout << N_signalevents[i] << " ";
   //     cout <<"\n";   
		bool savesignal_1=true;
		Waves.clear();
		
			if (savesignal_1 && !isTrg && point.first == channel ) { //new graphs with arrows on the found peaks
			Waves[channel].fillWave(point.second,dim);
			signal->cd();
			if (firstEntering_filter){
			tmpCvsignal_1.push_back( new TCanvas(Form("CvSignal_1_ev%d",jentry),Form("tmpSignal_1_ev%d",jentry)) );
			tmpCvsignal_1.back()->Divide(3,4);
			firstEntering_filter=false;
			}

			tmpCvsignal_1.back()->cd(channel-nTriggerChannels+1);
			tmpsignal_1.push_back( new TGraph ( dim, &X[0], &Waves[channel].Y[0]) );
			tmpsignal_1.back()->GetXaxis()->SetTitle("time [ns]");
			tmpsignal_1.back()->SetTitle(Form("tmpSignal_afterFlt-Ch%d_ev%d",channel,jentry));
			tmpsignal_1.back()->GetYaxis()->SetTitleOffset(1.4);
			tmpsignal_1.back()->GetYaxis()->SetTitle("Voltage [V]");
      
			//tmpsignal_1.back()->GetYaxis()->SetRangeUser(-0.1,0.4);
			tmpsignal_1.back()->Draw("AL");
			counting_filter++;
			
			for (int ipk=0; ipk<NPeak_1; ipk++){
				TMarker *tm = new  TMarker(X[pkPos_1[ipk]+skipFstBin], pkHgt_1[ipk], 23);
				tm->SetMarkerSize(1.5);
				tm->SetMarkerColor(2);
				tm->Draw();
			}
	    	if(counting_filter==(nMaxCh-nTriggerChannels+1)){					
				tmpCvsignal_1.back()->Write();
				}
				theFile->cd("/");
			} //if on representing the found peaks 
	    

		} //channel loop
        
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

}

