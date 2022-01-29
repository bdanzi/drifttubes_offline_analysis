#define read_data_cxx
#include "read_data.h"
#include "funcUtils.h"

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


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////Define Data Containers//////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
vector<float>X;//definisce l'asse x


// Temporary container for tmp Graphs and Canvas
std::vector<TGraph *> tmpWaves;          //grafico delle forme d'onda
std::vector<TCanvas *> tmpCv;            //canvas forme d'onda
std::vector<TGraph *> tmpFFTAmpl;        //grafico per l'ampiezza
std::vector<TGraph *> tmpFFTFltAmpl;     //grafico per ampiezza filtrata
std::vector<TGraph *> tmpFFTphi;         //grafico per la fase
std::vector<TCanvas *> tmpCvFFT;         //canvas per la fft
std::vector<TGraph *> tmpFltWaves;     //grafici per funzioni d'onda filtrate (primo filtro:retta lineare di taglio)
std::vector<TCanvas *> tmpCvFlt;       //canvas forme d'onda filtrate(primo filtro:retta lineare di taglio)
std::vector<TH1F *>  tmpFFTAmplPeak;   //grafico per lo spettro in frequenza dell'ampiezza con filtro di smooth
//std::vector<TGraph *> tmpnoise_1;      //grafici per il noise individuato dal filtro RC.(singolo filtro)
//std::vector<TCanvas *> tmpCvnoise_1;   //canvas per il noise individuato dal filtro RC(singolo filtro).
std::vector<TCanvas *> tmpCvPeak;
std::vector<TGraph *> tmpsignal_1;       //grafici per il segnale finale dopo i filtri 
std::vector<TCanvas *> tmpCvsignal_1;    //canvas per il segnale finale dopo i filtri 
std::vector<TGraph *> tmpNoise_total;    //grafici per il noise totale 
std::vector<TCanvas *> tmpCvNoise_total; //canvas per il noise totale 

//void doChDiff(int evN, std::pair<int,int> &chDiff, std::map<int, bool> &isFull, diffWvCont &diffWaves, WvCont &Waves, WvCont &Smt4Waves, int SF, hDiffCont &hDiff, hDiffCont &hDiffPed, TFile *theFile);
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////MAIN LOOP function//////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void read_data::Loop(Char_t *output, Int_t MidEv,Int_t eventn, Bool_t evalWaveCut)
{

	tmax=_tmax;  //to fix compilation problem

	skipFstBin=580;//525; //set Bins to Skip for El Cal.
	skipLstBin=530;//475;

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
	TFile * theFile = new TFile(outChar,"RECREATE");  //crea un nuovo file.root e se ne esiste gi� uno lo sovrascrive

	std::map<int, hstPerCh *> HstPerCh;               //mappa per istogrammi sui segnali filtrati.
	//std::map<int, hstPerCh *> HstSF10PerCh;
   	//   theFile->cd("/");
   	//   for (int channel=0; channel<=33; channel++){
   	//	   TDirectory *chDir = theFile->mkdir(Form("H-Ch%d",channel));
   	//	   chDir->cd();
   	//	   HstPerCh.insert(make_pair(channel,new hstPerCh(channel)));
   	//	   theFile->cd("/");
   	//   }
   
	TH1F *hIntRatio = new TH1F("hIntRatio","Ratio of the integrals between chs",1000,0,100);
	TH1F *hMaxRatio = new TH1F("hMaxRatio","Ratio of the maximums between chs",1000,0,100);

   	//   //Fitting functions
   	//   TF1* f1 = new TF1("f1", "gaus",  0.0004, 0.004);
   	//   TF1* f2 = new TF1("f2", "gaus",  0.0004, 0.004);
   	//
   	//   TF1* f3 = new TF1("f3", "landau", 0.2,5.);
   	//   TF1* f4 = new TF1("f4", "gaus", -5.0,0.2);
   	//   TF1* total = new TF1("total", "f4 + f3", -5.0, 5.0);
   
   	//Define counter
	  int count =0; 
 
/***********creazione delle directory***********/
	theFile->cd("/");
	TDirectory *waveDir = theFile->mkdir("Waves");
	theFile->cd("/");
//	TDirectory *waveDiffDir = theFile->mkdir("WavesDiff");
//	theFile->cd("/");
	TDirectory *waveFFTDir = theFile->mkdir("WavesFFT");
	theFile->cd("/");
	//TDirectory *waveFltDir = theFile->mkdir("WavesFlt");
	//theFile->cd("/");
  TDirectory *waveFltPeakDir = theFile->mkdir("Amplfft");      //directory per lo spettro in frequenza dell'ampiezza con smooth
  theFile->cd("/");
  //TDirectory *noise_RC = theFile->mkdir("noise_RC");           //directory per la waveform di noise individuata dal filtro RC(noise_1,singolo filtro)
  //theFile->cd("/");
  TDirectory *signal = theFile->mkdir("signal_Afterflt"); //directory per la waveform di segnale individuata dopo tutti i filtri
  theFile->cd("/");
  TDirectory *noise = theFile->mkdir("noise_total");      //directory per il noise totale
  theFile->cd("/");
  
  
//TDirectory *FilterblsDir = theFile->mkdir("Filterbls");//directory per la differenza bin per bin
//theFile->cd("/");
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	//                                               //  Loop on entries //                                                  //

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	std::vector<int> trigCh; //vettore di interi per i canali di trigger.
	trigCh.clear();
#ifndef _OSC
	trigCh.push_back(32);
	trigCh.push_back(33);
#else
	trigCh.push_back(7);//osc1
//	trigCh.push_back(8);//osc2
#endif
//mappa per i tagli: definire cosa � segnale e cosa � rumore
	std::map<int, std::pair<float, float> > waveCut;
	waveCut.clear();
	int nSig=3;// numero di sigma.(era 3)
	std::map<int, int> nSelEv;
	nSelEv.clear();

	if (fChain == 0) return;                      //� il puntatore TTree

	Long64_t nentries = fChain->GetEntriesFast(); //Return the number of entries as of the last check.
	cout << "Number of entries in the tree for real data is= " << nentries << endl;
	//   Long64_t nbytes = 0, nb = 0;

	Int_t frstEv=0;                               //inizializzo il primo evento a zero
	Int_t lastEv=nentries;                        //l'ultimo evento coincide con le nentries, che si prende con la funzione get...sopra
	if (eventn>=0&&eventn<lastEv) {
		frstEv=eventn;
		lastEv=eventn+1;
	} else if (eventn<-1 && eventn>-nentries){
		lastEv=-eventn;
	}

	//   bool evalWaveCut=false;
	TString wCutFName = Form("wave-Cut-%s.txt",out.Data());

	if (evalWaveCut) {

		std::ofstream cutf;
		cutf.open(wCutFName.Data(), std::ofstream::out);//stampa su file. 

		int SF=10;

		//	   std::map<int, TF1*> tmpFitRes;

		for (Long64_t jentry=frstEv; jentry<lastEv;jentry++) {
			Long64_t ientry = LoadTree(jentry); //The function finds the corresponding Tree and returns the entry number in this tree.  
			if (ientry < 0) break;
			//     nb = fChain->GetEntry(jentry);   nbytes += nb;
			fChain->GetEntry(jentry);
			WvCont Waves;
			Waves.clear();
			WvCont Smt10Waves;
			Smt10Waves.clear();
#ifndef _OSC
			if(jentry==frstEv) {
				// 1024 channel of ADC = 1 micro sec => every channel is 0.9766 nano sec
				double time = (_tmax*1.0e+9)/((float)dim);//in nano sec
				X.clear();
				for (int i = 0; i < dim; ++i) { X.push_back(time *(i+1)); }
			}
			for (auto point : wd->getX742Data()) {
				for (int channel=0; channel<=nMaxCh+2; channel++){
#else //lettura dei file di dati .root 
			if(jentry==frstEv) {
				nMaxCh=wd->getRunHeader()->getNumberOfDevicesInRun()*8-1;//i primi canali sono dell'oscilloscopio 1, gli altri sono dell'oscilloscopio 2.
				dim=wd->getRunHeader()->getNumberOfPointsPerChannel();
				_tmax=(dim+1)*wd->getRunHeader()->getTimeResolOfChannel();
				double time = (_tmax*1.0e+9)/((float)dim);//in nano sec
				X.clear();
				for (int i = 0; i < dim; ++i) { X.push_back(time *(i+1)); }
			}
			for (auto point : wd->getXOSCData()) {
				for (int channel=0; channel<=nMaxCh; channel++){
#endif
					bool isTrg=false;
					for(auto trgCh : trigCh) { if (channel==trgCh) { isTrg=true; } }
					if (point.first == channel && !isTrg) {
						if ( HstPerCh.find(channel)==HstPerCh.end() ) {
							TDirectory *chDir = theFile->mkdir(Form("H-Ch%d_signal",channel));
							chDir->cd();
							HstPerCh.insert(make_pair(channel,new hstPerCh(channel)));
							//HstSF10PerCh.insert(make_pair(channel,new hstPerCh(channel,SF)));
							theFile->cd("/");
						}
						//cout << "Channel is= " << channel << endl;
						//Waves[channel].fillWave(point.second.size(),point.second.data());
						Waves[channel].fillWave(point.second);
						std::vector<float> tmpWSF10=smooth(Waves[channel].Y,SF);
           	Smt10Waves[channel].fillWave( tmpWSF10 );

						/*((hstPerCh*)HstSF10PerCh[channel])->hBsl->Fill(((wave)Smt10Waves[channel]).bsln);
						((hstPerCh*)HstSF10PerCh[channel])->hInteg->Fill(((wave)Smt10Waves[channel]).integ);
						((hstPerCh*)HstSF10PerCh[channel])->hIntegN->Fill(((wave)Smt10Waves[channel]).nInteg());
						((hstPerCh*)HstSF10PerCh[channel])->hRms->Fill(((wave)Smt10Waves[channel]).rms);
						((hstPerCh*)HstSF10PerCh[channel])->hMaxV->Fill(((wave)Smt10Waves[channel]).max);
						((hstPerCh*)HstSF10PerCh[channel])->hMaxVN->Fill(((wave)Smt10Waves[channel]).nMax());*/
            
						bool saveWave=false; //non salvo, quindi non vediamo l'out di questo pezzo di codice.
						if (saveWave) {
							//if (channel==0 && (Smt10Waves[channel].nMax()>0.0035 && Smt10Waves[channel].nMax()<0.005) ) {
							if (channel==11 && (Smt10Waves[channel].n1Max()>0.0125/*(0.003501+5*0.0009636)*/) ) {
								waveDir->cd();
								tmpCv.push_back( new TCanvas(Form("Cv-Ch%d_ev%d",channel,jentry),Form("tmpWave-Ch%d_ev%d",channel,jentry)) );
								tmpCv.back()->cd();
								tmpWaves.push_back( new TGraph ( dim, &X[0], &Waves[channel].Y[0]) );
								tmpWaves.back()->GetXaxis()->SetTitle("time [ns]");
								tmpWaves.back()->SetTitle(Form("tmpWave-Ch%d_ev%d",channel,jentry));
								tmpWaves.back()->GetYaxis()->SetTitleOffset(1.4);
								tmpWaves.back()->GetYaxis()->SetTitle("Volt");
								tmpWaves.back()->Draw("AL");
								tmpWaves.push_back( new TGraph ( dim-2*SF, &X[SF], &Smt10Waves[channel].Y[0]) );
								tmpWaves.back()->SetLineColor(kBlue);
								tmpWaves.back()->Draw("Lsame");
								tmpCv.back()->Write();
								theFile->cd("/");
							}
						}
					}
				}
			}
		}
////-------------------FIT GAUSS---------------
	/*	for (int channel=0; channel<=nMaxCh; channel++){
			if (HstSF10PerCh.find(channel)!=HstSF10PerCh.end()){
				((hstPerCh*)HstSF10PerCh[channel])->hMaxVN->Fit("gaus","","",0,0.005);
				//tmpFitRes.insert( make_pair(channel, ((hstPerCh*)HstSF10PerCh[channel])->hMaxVN->GetFunction("gaus") ) );
				TF1 *tmpFit = ((hstPerCh*)HstSF10PerCh[channel])->hMaxVN->GetFunction("gaus") ;
				float tmpMean = tmpFit->GetParameter(1);
				float tmpSigma = tmpFit->GetParameter(2);
				waveCut.insert( make_pair(channel, make_pair(tmpMean,tmpSigma) ) );                   //riempiamo la mapa waveCut definita all'inizio
			}
		}*/
		std::cout<<"Fit of Intergeal SWF 10 fit results:"<<std::endl;
		for (int channel=0; channel<=nMaxCh; channel++){
			if (waveCut.find(channel)!=waveCut.end()) {                                             //fintanto che i canali non sono finiti scrivi su file.
				std::cout<<channel<<"\t"<<waveCut[channel].first<<"\t"<<waveCut[channel].second<<std::endl; 
        //scrittura, first � la media, second � la sigma, per ogni canale.
				cutf<<channel<<"\t"<<waveCut[channel].first<<"\t"<<waveCut[channel].second<<std::endl;//cutf , scrive le info sul file.
			}
		}

		cutf.close();

	} else {
		std::ifstream cutf(wCutFName.Data());
		int tmpCh;
		float tmpMean, tmpSgm;
		while (!cutf.eof()) {
			cutf>>tmpCh>>tmpMean>>tmpSgm;
			waveCut.insert( make_pair(tmpCh, make_pair(tmpMean,tmpSgm) ) );
		}
		//	   std::cout<<"Fit of Intergeal SWF 10 fit results:"<<std::endl;
		//	   for (int channel=0; channel<=nMaxCh; channel++){
		//		   if (waveCut.find(channel)!=waveCut.end()) {
		//			   std::cout<<channel<<"\t"<<waveCut[channel].first<<"\t"<<waveCut[channel].second<<std::endl;
		//		   }
		//	   }
		cutf.close();
	}

	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	//Try second loop

	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	bool testSelCut=true;
	int SF = 10;

	theFile->cd("/");

	std::vector< std::pair<int,int> > chToDiff;  //vettore di pair: canali da sottrarre
#ifndef _OSC
	chToDiff.push_back( make_pair(1,2) );
#else

#endif

/*contenitori per la differenza tra canali*/
	//hDiffCont hDiff;
	//hDiffCont hDiffPed;

/*	for ( auto chsDiff : chToDiff ) {
		TDirectory *diffDir = theFile->mkdir(Form("DiffCh%d%d",chsDiff.first,chsDiff.second));
		diffDir->cd();
		hDiff.insert( make_pair(chsDiff, hstDiffCh(chsDiff.first,chsDiff.second) ) );
		hDiffPed.insert( make_pair(chsDiff, hstDiffCh(chsDiff.first,chsDiff.second,false) ) );
		theFile->cd("/");
	}*/
 
/*Ricerca dei picchi di noise sulle FFT e plot dei risultati*/
  TSpectrum *spectrum = new TSpectrum(10,1);
 // istogramma distribuzione dei picchi
 // TH1F *peak_HF= new TH1F("peak_HF","High frequency peak distribution", 100,3.2e+9,3.7e+9);
  /*peak_HF->GetXaxis()->SetTitle("Peak position (#omega [Hz])");
   peak_HF->GetYaxis()->SetTitle("Entries");*/

  TH1F *peak_LF= new TH1F("peak_LF","Low frequency peak distribution",100, 2.0e+7,9.0e+8);
  peak_LF->GetXaxis()->SetTitle("Peak position (#omega [Hz])");
  peak_LF->GetYaxis()->SetTitle("Entries");

/*Istogramma distribuzione di min value e min pos */
  TH1F *hmin_value= new TH1F("hmin_value","Min value distribution", 100,0,-0.01e-4);
  hmin_value->GetXaxis()->SetTitle("Min_value");
  hmin_value->GetYaxis()->SetTitle("Entries");
  TH1F *hmin_pos= new TH1F("hmin_pos","Position min value distribution", 1000,0,1000);
  hmin_pos->GetXaxis()->SetTitle("Position min value");
  hmin_pos->GetYaxis()->SetTitle("Entries");
 /**array di errore*/
std::vector<float> Xerror;
std::vector<float> Yerror;


  //***********************mappa per la differenza bin per bin di forme d'onda e forme d'onda dopo il filtro.
 // std::map<int,TH1F*> bin_bin;
  /*for (int channel=0; channel<=nMaxCh+2; channel++){
   bin_bin.insert( make_pair(channel, new TH1F(Form("bin_bin_ch%d",channel),"",1000,-0.5,0.5)) );
  }*/
  
/*Contenitori di waveform*/  
  WvCont Waves;
  WvCont Waves_diff;      
  WvCont FltWaves;
  WvCont FltWavesSM;
  WvCont noise_1;
  WvCont Waves_signal_1; 
  WvCont Waves_noise_tot;
  WvCont Smt10Waves;
  WvCont Smt4Waves;
  
	//diffWvCont diffWaves;
	fftCont Wffts; 
  fftCont CR;
  fftCont RC;
  fftCont Notch_1;
  fftCont Notch_2;
  fftCont Notch_3;
  fftCont Notch_4;
  fftCont Notch_5;
//  fftCont ContforNotch;
 //-------------------------------------------creazione mappa per smooth di SG------------------------------------------------------
  WvCont SmtSGWaves;
  WvCont SmtSGWaves_7;
  WvCont SmtSGWaves_9;
  WvCont SmtSGWaves_13;
  WvCont SmtSGWaves_23;
  WvCont SmtSGWaves_4_9;
  WvCont SmtSGWaves_4_7;
  WvCont smtSG;
  
//---------------------------------------------------------------SECONDO loop, da qui comincia la vera analisi.----------------------------------------
	for (Long64_t jentry=frstEv; jentry<lastEv;jentry++) {

		Long64_t ientry = LoadTree(jentry);
		if (ientry < 0) break;

		//     nb = fChain->GetEntry(jentry);   nbytes += nb;
		fChain->GetEntry(jentry);
		// if (Cut(ientry) < 0) continue;

		std::map<int, bool> isFull; //mappa che contiene indice di canale e mi dice se questo � pieno o vuoto.
		for (int ich=0; ich<=nMaxCh; ++ich ) { isFull[ich]=false; }

		Waves.clear();
		Smt10Waves.clear();
		Smt4Waves.clear();
		Wffts.clear();
		CR.clear();
		RC.clear();
		Notch_1.clear();
		Notch_2.clear();
		Notch_3.clear();
		Notch_4.clear();
		Notch_5.clear();
		FltWaves.clear();
    FltWavesSM.clear();
		noise_1.clear();
		SmtSGWaves.clear();
		SmtSGWaves_7.clear();
		SmtSGWaves_9.clear();
		SmtSGWaves_13.clear();
		SmtSGWaves_4_9.clear();
		SmtSGWaves_4_7.clear();
		SmtSGWaves_23.clear();
		Waves_signal_1.clear();
		Waves_noise_tot.clear();
   
		//for ( auto chsDiff : chToDiff ) { diffWaves[chsDiff].clear(); }

#ifndef _OSC
		if(jentry==frstEv) {
			// 1024 channel of ADC = 1 micro sec => every channel is 0.9766 nano sec
			double time = (_tmax*1.0e+9)/((float)dim);//in nano sec
			X.clear();
			for (int i = 0; i < dim; ++i) { X.push_back(time *(i+1)); }
		}
		for (auto point : wd->getX742Data()) {
			//cout<<"Discrimentor "<<endl;
			//cout<<"**********************"<<endl;

			for (int channel=0; channel<=nMaxCh+2; channel++){
#else
		if(jentry==frstEv) {
			nMaxCh=wd->getRunHeader()->getNumberOfDevicesInRun()*8-1;
			dim=wd->getRunHeader()->getNumberOfPointsPerChannel();
			_tmax=(dim+1)*wd->getRunHeader()->getTimeResolOfChannel();
			tmax=_tmax;
			//std::cout<<"MaxCh "<<nMaxCh<<" dim "<<dim<<" tmax "<<_tmax<<std::endl;
			double time = (_tmax*1.0e+9)/((float)dim);//in nano sec
			X.clear();
			for (int i = 0; i < dim; ++i) { X.push_back(time *(i+1)); }
			std::cout<<"Entries "<<jentry<<" Event n. "<<wd->getEventNumber()<<std::endl;
		}
		for (auto point : wd->getXOSCData()) {
			for (int channel=0; channel<=nMaxCh; channel++){
      
#endif
       
/*Esclusione dei canali di trigger.*/
				bool isTrg=false; //serve a controllare che il canale su cui facciamo l'analisi non sia quello di trigger.
				for(auto trgCh : trigCh) { if (channel==trgCh) { isTrg=true; } }
				if (point.first == channel  ) {
					//cout << "Channel is= " << channel << endl;
					//		 Waves[channel].fillWave(point.second.size(),point.second.data());
//					std::cout<<"Channel "<<channel<<" original wave"<<std::endl;
					Waves[channel].fillWave(point.second);
//					std::cout<<"original wave end"<<std::endl;
					//std::cout<<"nPt "<<Waves[channel].nPt()<<std::endl;
                            
  /****----------------------------TRASFORMATA DI FOURIER---------------------------------*****/                          
					FFT(Waves[channel],Wffts[channel]);//trasformata di Fourier
					double *realFltFFT = new double[Waves[channel].nPt()];
					double *imgFltFFT = new double[Waves[channel].nPt()];
                 
 /****-----------------------------FUNZIONI DI FILTRO-------------------------------------*****/
					//filterWave(Wffts[channel],realFltFFT,imgFltFFT,3e+6,5e+6);  //funzione filtro(taglio)
					// filterWaveSM(Wffts[channel],realFltFFT,imgFltFFT,9,2);     //smooth su SG
					//InverseFFT(realFltFFT,imgFltFFT,Waves[channel].nPt(),FltWavesSM[channel]);//trasformata inversa, con singolo filtro.
 					//double *amplFltFFT = new double[Waves[channel].nPt()];
 					//for (int ip=0; ip<Waves[channel].nPt()/2; ++ip) {
						//std::cout<<"omega "<<Wffts[channel].omega[ip]<<" Flt: real "<<realFltFFT[ip]<<" img "<<imgFltFFT[ip]<<std::endl;
					//amplFltFFT[ip]=TMath::Sqrt(realFltFFT[ip]*realFltFFT[ip]+imgFltFFT[ip]*imgFltFFT[ip]);
 				//	}
 				//	for (int ip=Waves[channel].nPt()/2; ip<Waves[channel].nPt(); ++ip) {
 				//		amplFltFFT[ip]=0.0;
 				//	} 
           
//**----------- RIEMPIMENTO grafici PER L'AMPIEZZA CON SMOOTH(cerchiamo i picchi di noise sulle forme d'onda con smooth di sg)---------*/

/*--------  CONVERSIONE A ISTOGRAMMA-----------------*/
               Double_t *peaks=0x0 ;
               int numOfPeaks=0;
					     bool saveAmplfft = false;//true
//               if(/*saveAmplfft &&*/ ! isTrg){
//               waveFltPeakDir->cd();
//
//               //************RIEMPIMENTO CANVAS*********
//               tmpCvPeak.push_back( new TCanvas(Form("CvPeak-Ch%d_ev%d",channel,jentry),Form("tmpAmplfft-Ch%d_ev%d",channel,jentry)) );
//               //tmpCvPeak.back()->cd()->SetLogy();
//
//               //************RIEMPIMENTO GRAFICI*****
//               tmpFFTAmplPeak.push_back( new TH1F(Form("FFTampl-Ch%d_ev%d",channel,jentry),";#omega [Hz];Amplitude",
//                                           dim/2,Wffts[channel].omega[0],Wffts[channel].omega[dim/2-1]));
//
//               std::vector<float> fftAmplFrNtch;
//               fftAmplFrNtch.clear();
//               std::vector<int> fftIdxFrNtch;
//               fftIdxFrNtch.clear();
//               for(int i=0; i < (dim/2); ++i) {
//               if(Wffts[channel].omega[i]>2.e+8&&Wffts[channel].omega[i]<5.5e+8){
//                   fftAmplFrNtch.push_back(Wffts[channel].ampl[i]);
//                   fftIdxFrNtch.push_back(i);
//                 }
//               }
//
//               std::vector<float> fftAmplFrNtch_smt=smoothSG(fftAmplFrNtch,7,2);
//
//               for (int i=3; i<fftIdxFrNtch.size()-3; ++i) { //2=(int)5/2
//                 if (fftAmplFrNtch_smt[i]>0.0 && fftAmplFrNtch_smt[i]<100.0) tmpFFTAmplPeak.back()->SetBinContent(fftIdxFrNtch[i]+1,fftAmplFrNtch_smt[i]);
//               }
//               tmpFFTAmplPeak.back()->Draw();
//
//               //**********Ricerca dei picchi nello spettro di frequenza******
//
//               Int_t npeaks =  spectrum -> Search(tmpFFTAmplPeak.back(),4.0,"nobackground",0.7);
//               cout << "Peaks found " << npeaks << " evento "<< jentry << " canale " << channel <<endl ; ;
//               peaks = spectrum->GetPositionX() ;
//               numOfPeaks=spectrum->GetNPeaks();
//               for (int j=0 ; j<npeaks ; j++){
////                  cout << "Peak at = " << peaks[j] << " evento "<< jentry << " canale " << channel <<endl ;
//                  /*if(peaks[j]>3.e+9){
//                    peak_HF->Fill(peaks[j]);
//                   }*/
//                 // if( peaks[j]>2.5e+8 && peaks[j]<5.e+8 ){
//                    peak_LF->Fill(peaks[j]);
//                  // }
//                  }
//               //scrittura
//               tmpCvPeak.back()->Write();
//               theFile->cd("/");
//            }
                                                                                 
  /****---------------------------------------------------------------------------------------*****/           

          if (true/*peaks!=0x0*/){
            wavefft tmpWfft;
            tmpWfft.fillAP(Wffts[channel].nPt,Wffts[channel].ampl,Wffts[channel].phi,Wffts[channel].omega);         
            for(int i=0; i<numOfPeaks; ++i){
              cout << "Peak at = " << peaks[i] << " evento "<< jentry << " canale " << channel <<endl ;
  					  filterNotch(tmpWfft,realFltFFT,imgFltFFT,peaks[i]/*(TMath::TwoPi()*66e+6)*/,TMath::TwoPi()*20e+6,2);////
              tmpWfft.clear();
              tmpWfft.fillRM(Wffts[channel].nPt,realFltFFT,imgFltFFT,Wffts[channel].omega);
            }
           	Notch_1[channel].fillAP(tmpWfft.nPt,tmpWfft.ampl,tmpWfft.phi,tmpWfft.omega);
          }
				//	filterNotch(Notch_1[channel],realFltFFT,imgFltFFT,(TMath::TwoPi()*46e+6),TMath::TwoPi()*10e+6,2);//46 MHz 10 MHz///
				//	Notch_2[channel].fillRM(Wffts[channel].nPt,realFltFFT,imgFltFFT,Wffts[channel].omega);                          ///
//				  filterNotch(Notch_1[channel],realFltFFT,imgFltFFT,(TMath::TwoPi()*538e+6),TMath::TwoPi()*20e+6,2);              ///
//					Notch_2[channel].fillRM(Wffts[channel].nPt,realFltFFT,imgFltFFT,Wffts[channel].omega);                          ///
					//filterNotch(Notch_3[channel],realFltFFT,imgFltFFT,(TMath::TwoPi()*5e+6),(TMath::TwoPi()*50e+6),3);//passa alto
					//filterWaveCR(Notch_3[channel],realFltFFT,imgFltFFT,1.0/(TMath::TwoPi()*5e+6));
					//Notch_4[channel].fillRM(Wffts[channel].nPt,realFltFFT,imgFltFFT,Wffts[channel].omega);
					//filterNotch(Notch_4[channel],realFltFFT,imgFltFFT,(TMath::TwoPi()*700e+6),TMath::TwoPi()*10e+6,0);//passa basso
				//	filterWaveRC(Notch_2[channel],realFltFFT,imgFltFFT,1.0/(TMath::TwoPi()*500e+6));
					//Notch_3[channel].fillRM(Wffts[channel].nPt,realFltFFT,imgFltFFT,Wffts[channel].omega);
          //filterpole_zero(Notch_4, realFltFFT,imgFltFFT, 50e+6,1,1);
    	    //Notch_5[channel].fillRM(Wffts[channel].nPt,realFltFFT,imgFltFFT,Wffts[channel].omega);                           
					filterWaveBsl(Notch_1[channel],realFltFFT,imgFltFFT);//filtro sulla baseline.                                        ///
					InverseFFT(realFltFFT,imgFltFFT,Waves[channel].nPt(),FltWaves[channel]);//trasformata inversa
					Waves_signal_1[channel].fillWave(FltWaves[channel].Y);
					double *amplFltFFT_2 = new double[Waves[channel].nPt()];
					for (int ip=0; ip<Waves[channel].nPt()/2; ++ip) {
						amplFltFFT_2[ip]=TMath::Sqrt(realFltFFT[ip]*realFltFFT[ip]+imgFltFFT[ip]*imgFltFFT[ip]);
					}
					for (int ip=Waves[channel].nPt()/2; ip<Waves[channel].nPt(); ++ip) {
						amplFltFFT_2[ip]=0.0;
					}

					//CREAZIONE DIRECTORY E RIEMPIMENTO ISTOGRAMMI, SUI SEGNALI FILTRATI!!!
					if(!isTrg){

						if ( HstPerCh.find(channel)==HstPerCh.end() ){
							TDirectory *chDir = theFile->mkdir(Form("H-Ch_%d_signal",channel));
							chDir->cd();
							HstPerCh.insert(make_pair(channel,new hstPerCh(channel)));
							//	if (!evalWaveCut)  HstSF10PerCh.insert(make_pair(channel,new hstPerCh(channel,SF))); //make_pair funzione per creare una coppia.
							theFile->cd("/");
						}

						((hstPerCh*)HstPerCh[channel])->hBsl->Fill(((wave)Waves_signal_1[channel]).bsln);                               
            ((hstPerCh*)HstPerCh[channel])->hBsl->GetXaxis()->SetTitle("Baseline [V]");         
            ((hstPerCh*)HstPerCh[channel])->hBsl->GetYaxis()->SetTitle("Entries");         
						((hstPerCh*)HstPerCh[channel])->hInteg->Fill(((wave)Waves_signal_1[channel]).integ);
						std::cout<<"ch "<<channel<<" integ "<<((wave)Waves_signal_1[channel]).integ<<std::endl;
						((hstPerCh*)HstPerCh[channel])->hIntegN->Fill(((wave)Waves_signal_1[channel]).n1Integ());
						((hstPerCh*)HstPerCh[channel])->hIntegInR->Fill(((wave)Waves_signal_1[channel]).integInR);                                                
						((hstPerCh*)HstPerCh[channel])->hIntegNInR->Fill(((wave)Waves_signal_1[channel]).n1IntegInR());
            ((hstPerCh*)HstPerCh[channel])->hIntegNInR-> GetXaxis()->SetTitle("Integral [V]");
            ((hstPerCh*)HstPerCh[channel])->hIntegNInR-> GetYaxis()->SetTitle("Entries");
						((hstPerCh*)HstPerCh[channel])->hRms->Fill(((wave)Waves_signal_1[channel]).rms);
						((hstPerCh*)HstPerCh[channel])->hRms->GetXaxis()->SetTitle("Rms [V]");
						((hstPerCh*)HstPerCh[channel])->hRms->GetYaxis()->SetTitle("Entries");   
                                                                                                          
						((hstPerCh*)HstPerCh[channel])->hMaxV->Fill(((wave)Waves_signal_1[channel]).max);
						((hstPerCh*)HstPerCh[channel])->hMaxVN->Fill(((wave)Waves_signal_1[channel]).n1Max());
						((hstPerCh*)HstPerCh[channel])->hMaxVInR->Fill(((wave)Waves_signal_1[channel]).maxInR);
						((hstPerCh*)HstPerCh[channel])->hMaxVNInR->Fill(((wave)Waves_signal_1[channel]).n1MaxInR());
						((hstPerCh*)HstPerCh[channel])->hMaxVNInR->GetXaxis()->SetTitle("Max_value [V]");
						((hstPerCh*)HstPerCh[channel])->hMaxVNInR->GetYaxis()->SetTitle("Entries");
                                                                                                                      
						((hstPerCh*)HstPerCh[channel])->hMinV->Fill(((wave)Waves_signal_1[channel]).min);
						((hstPerCh*)HstPerCh[channel])->hMinVN->Fill(((wave)Waves_signal_1[channel]).n1Min());
						((hstPerCh*)HstPerCh[channel])->hMinVInR->Fill(((wave)Waves_signal_1[channel]).minInR);
						((hstPerCh*)HstPerCh[channel])->hMinVNInR->Fill(((wave)Waves_signal_1[channel]).n1MinInR());
						((hstPerCh*)HstPerCh[channel])->hMinVNInR->GetYaxis()->SetTitle("Entries");
						((hstPerCh*)HstPerCh[channel])->hMinVNInR->GetXaxis()->SetTitle("Min value [V]"); 
      
             /*istogrammi sul minimo delle wave non filtrate*/  
            ((hstPerCh*)HstPerCh[channel])->hMinVoriginalW->Fill(((wave)Waves[channel]).min);
						((hstPerCh*)HstPerCh[channel])->hMinVNoriginalW->Fill(((wave)Waves[channel]).n1Min());
						((hstPerCh*)HstPerCh[channel])->hMinVInRoriginalW->Fill(((wave)Waves[channel]).minInR);
						((hstPerCh*)HstPerCh[channel])->hMinVNInRoriginalW->Fill(((wave)Waves[channel]).n1MinInR());
						((hstPerCh*)HstPerCh[channel])->hMinVNInRoriginalW->GetXaxis()->SetTitle("Min value [V]");
						((hstPerCh*)HstPerCh[channel])->hMinVNInRoriginalW->GetYaxis()->SetTitle("Entries");
      
            /*istogrammi sul massimo delle wave non filtrate*/    
            ((hstPerCh*)HstPerCh[channel])->hMaxVoriginalW->Fill(((wave)Waves[channel]).max);
						((hstPerCh*)HstPerCh[channel])->hMaxVNoriginalW->Fill(((wave)Waves[channel]).n1Max());
						((hstPerCh*)HstPerCh[channel])->hMaxVInRoriginalW->Fill(((wave)Waves[channel]).maxInR);
						((hstPerCh*)HstPerCh[channel])->hMaxVNInRoriginalW->Fill(((wave)Waves[channel]).n1MaxInR());
						((hstPerCh*)HstPerCh[channel])->hMaxVNInRoriginalW->GetXaxis()->SetTitle("Min_value [V]");
						((hstPerCh*)HstPerCh[channel])->hMaxVNInRoriginalW->GetYaxis()->SetTitle("Entries");
						if (/*Waves_signal_1[channel].n1IntegInR()>0.0 &&*/ Waves_signal_1[channel].n1IntegInR()< 0.03){
							cout<< " canale "<< channel << " evento del cavolo "<< jentry << " valore integrale " << Waves_signal_1[channel].n1IntegInR() << " baseline " << Waves_signal_1[channel].bsln << endl;
							std::cout<<"reg A "<<((wave)Waves_signal_1[channel]).RegA()<<" reg B "<<((wave)Waves_signal_1[channel]).RegB()<<std::endl;
						}
						if (channel==6) {
							if ( Waves_signal_1[channel].n1IntegInR()> 0.16 && Waves_signal_1[channel].n1IntegInR()< 0.18) {
								cout<< " canale "<< channel << " evento picco 1 "<< jentry << " valore integrale " << Waves_signal_1[channel].n1IntegInR() << " baseline " << Waves_signal_1[channel].bsln << endl;
								std::cout<<"reg A "<<((wave)Waves_signal_1[channel]).RegA()<<" reg B "<<((wave)Waves_signal_1[channel]).RegB()<<std::endl;
							} else if ( Waves_signal_1[channel].n1IntegInR()> 0.22 && Waves_signal_1[channel].n1IntegInR()< 0.24) {
								cout<< " canale "<< channel << " evento picco 2 "<< jentry << " valore integrale " << Waves_signal_1[channel].n1IntegInR() << " baseline " << Waves_signal_1[channel].bsln << endl;
								std::cout<<"reg A "<<((wave)Waves_signal_1[channel]).RegA()<<" reg B "<<((wave)Waves_signal_1[channel]).RegB()<<std::endl;
							}

						}
                                                                                                      
					}
             
					std::vector<float> tmpWSF10=smooth(Waves[channel].Y,SF);
					Smt10Waves[channel].fillWave( tmpWSF10 );

					/*	if (!evalWaveCut) {
						((hstPerCh*)HstSF10PerCh[channel])->hBsl->Fill(((wave)Smt10Waves[channel]).bsln);
						((hstPerCh*)HstSF10PerCh[channel])->hInteg->Fill(((wave)Smt10Waves[channel]).integ);
						((hstPerCh*)HstSF10PerCh[channel])->hIntegN->Fill(((wave)Smt10Waves[channel]).nInteg());
						((hstPerCh*)HstSF10PerCh[channel])->hRms->Fill(((wave)Smt10Waves[channel]).rms);
						((hstPerCh*)HstSF10PerCh[channel])->hMaxV->Fill(((wave)Smt10Waves[channel]).max);
						((hstPerCh*)HstSF10PerCh[channel])->hMaxVN->Fill(((wave)Smt10Waves[channel]).nMax());
					}*/

					//-------------------------------------Smooth sulle forme d'onda------------------------------------------------
					//media mobile
					SF=4;
					std::vector<float> tmpWSF4=smooth(Waves[channel].Y,SF); //richiama la funzione di smooth in funcUtils.h
					Smt4Waves[channel].fillWave( tmpWSF4 );//richiama la funzione per riempire le forme d'onda.

					//smooth SG

					/*int m=5;
					int k=2;
					std::vector<float> tmpWSG=smoothSG(Waves[channel].Y,m,k);
					//cout<<"vettore di SMOOTH"<<tmpWSG.at(0)<<"\n";
					SmtSGWaves[channel].fillWave(tmpWSG);*/

					/*int m=7;
					int k=2;
					std::vector<float> tmpWSG_7=smoothSG(Waves[channel].Y,m,k);
					SmtSGWaves_7[channel].fillWave(tmpWSG_7);*/


				  int	m=23;
					int	k=2;
					std::vector<float> tmpWSG_23=smoothSG(Waves[channel].Y,m,k);
					SmtSGWaves_23[channel].fillWave(tmpWSG_23);
                                                     
				  /*m=9;
					k=2;
					std::vector<float> tmpWSG_9=smoothSG(Waves[channel].Y,m,k);
					SmtSGWaves_9[channel].fillWave(tmpWSG_9);*/
                                                     
					/*m=13;
					k=2;
					std::vector<float> tmpWSG_13=smoothSG(Waves[channel].Y,m,k);
					SmtSGWaves_13[channel].fillWave(tmpWSG_13);*/

					/*m=9;
					k=4;
					std::vector<float> tmpWSG_4_9=smoothSG(Waves[channel].Y,m,k);
					SmtSGWaves_4_9[channel].fillWave(tmpWSG_4_9);*/

				/*	m=7;
					k=4;
					std::vector<float> tmpWSG_4_7=smoothSG(Waves[channel].Y,m,k);
					SmtSGWaves_4_7[channel].fillWave(tmpWSG_4_7);*/

					//---------------------------------------disegno dei grafici delle forme d'onda sperimentali con smooth---------------------------------
					bool saveWave=true; //true
					if (saveWave) {
						waveDir->cd();
						tmpCv.push_back( new TCanvas(Form("Cv-Ch%d_ev%d",channel,jentry),Form("tmpWave-Ch%d_ev%d",channel,jentry)) );
						tmpCv.back()->cd();
						tmpWaves.push_back( new TGraph ( dim, &X[0], &Waves[channel].Y[0]) );
						tmpWaves.back()->GetXaxis()->SetTitle("time [ns]");
						tmpWaves.back()->SetTitle(Form("tmpWave-Ch%d_ev%d",channel,jentry));
						tmpWaves.back()->GetYaxis()->SetTitleOffset(1.4);
						tmpWaves.back()->GetYaxis()->SetTitle("Volt");
						tmpWaves.back()->Draw("AL");

				/*	SF=10;
              tmpWaves.push_back( new TGraph ( dim-2*SF, &X[SF], &Smt10Waves[channel].Y[0]) );
   						tmpWaves.back()->SetLineColor(kBlue);
   						tmpWaves.back()->Draw("Lsame");*/

			/*		SF=4;
  						tmpWaves.push_back( new TGraph ( dim-2*SF, &X[SF], &Smt4Waves[channel].Y[0]) );
  						tmpWaves.back()->SetLineColor(kRed);
              tmpWaves.back()->Draw("Lsame");*/

						/*int m=5;
                        int k=2;
                        tmpWaves.push_back( new TGraph ( dim-2*m, &X[m], &SmtSGWaves[channel].Y[0]) );
                        tmpWaves.back()->SetLineColor(kGreen);
  						tmpWaves.back()->Draw("Lsame");*/ 

					/*	m=7;
						k=2;
						tmpWaves.push_back(new TGraph(dim-2*m,&X[m],&SmtSGWaves_7[channel].Y[0]));
						tmpWaves.back()->SetLineColor(kOrange);
						tmpWaves.back()->Draw("Lsame");  */

					/*	m=13;
						k=2;
						tmpWaves.push_back(new TGraph(dim-2*m,&X[m],&SmtSGWaves_13[channel].Y[0]));
						tmpWaves.back()->SetLineColor(kMagenta);
						tmpWaves.back()->Draw("Lsame"); */

            m=23;
						k=2;
						tmpWaves.push_back(new TGraph(dim-m,&X[(m-1)/2],&SmtSGWaves_23[channel].Y[0]));
						tmpWaves.back()->SetLineColor(kGreen);
						tmpWaves.back()->Draw("Lsame");

			/*			m=9;
						k=2;
						tmpWaves.push_back(new TGraph(dim-2*m,&X[m],&SmtSGWaves_9[channel].Y[0]));
						tmpWaves.back()->SetLineColor(kPink);
						tmpWaves.back()->Draw("Lsame");   */

				/*		m=7;
						k=4;
						tmpWaves.push_back(new TGraph(dim-2*m,&X[m],&SmtSGWaves_4_7[channel].Y[0]));
						tmpWaves.back()->SetLineColor(kAzure-1);
						tmpWaves.back()->Draw("Lsame");*/

						tmpCv.back()->Write();
						theFile->cd("/");
					}
//--------------------------------------Disegno dei grafici della trasformata-------------------------------
					bool saveWaveFFT=true; //true
					if (saveWaveFFT) {

						waveFFTDir->cd();

						tmpCvFFT.push_back( new TCanvas(Form("CvFFT-Ch%d_ev%d",channel,jentry),Form("tmpWaveFFT-Ch%d_ev%d",channel,jentry)) );
						tmpCvFFT.back()->Divide(1,3);
						tmpCvFFT.back()->SetLogx();
						tmpCvFFT.back()->cd(1)->SetLogy();//primo riquadro : ampiezza
						tmpFFTAmpl.push_back( new TGraph ( dim/2, &Wffts[channel].omega[0], &Wffts[channel].ampl[0]) );
						tmpFFTAmpl.back()->GetXaxis()->SetTitle("#omega [Hz]");
						tmpFFTAmpl.back()->SetTitle(Form("tmpWaveFFT-Ampl-Ch%d_ev%d",channel,jentry));
						tmpFFTAmpl.back()->GetYaxis()->SetTitleOffset(1.4);
						tmpFFTAmpl.back()->GetYaxis()->SetTitle("Amplitude");
						tmpFFTAmpl.back()->Draw("AL");
						tmpCvFFT.back()->cd(2);//secondo riquadro :fase
						tmpFFTphi.push_back( new TGraph ( dim/2, &Wffts[channel].omega[0], &Wffts[channel].phi[0]) );
						tmpFFTphi.back()->GetXaxis()->SetTitle("#omega [Hz]");
						tmpFFTphi.back()->SetTitle(Form("tmpWaveFFT-Phi-Ch%d_ev%d",channel,jentry));
						tmpFFTphi.back()->GetYaxis()->SetTitleOffset(1.4);
						tmpFFTphi.back()->GetYaxis()->SetTitle("#Phi [rad]");
						tmpFFTphi.back()->Draw("AL");
						tmpCvFFT.back()->cd(3)->SetLogy();//terzo riquadro: taglio.
						tmpFFTAmpl.push_back( new TGraph ( dim/2, &Wffts[channel].omega[0],amplFltFFT_2) );
						tmpFFTAmpl.back()->GetXaxis()->SetTitle("#omega [Hz]");
						tmpFFTAmpl.back()->SetTitle(Form("tmpWaveFlt-Ampl-Ch%d_ev%d",channel,jentry));
						tmpFFTAmpl.back()->GetYaxis()->SetTitleOffset(1.4);
						tmpFFTAmpl.back()->GetYaxis()->SetTitle("Amplitude");  
						tmpFFTAmpl.back()->Draw("AL");

						tmpCvFFT.back()->Write();
						theFile->cd("/");

					}
           
          
              
					//********************SOTTRAZIONE BIN PER BIN********************
                                            
					//std::map<int,TH1F*> bin_bin;
					/*   if(!isTrg){
               if(bin_bin.find(channel)==bin_bin.end()){
               bin_bin.insert( make_pair(channel, new TH1F(Form("bin_bin_ch%d",channel),"",1000,-0.5,0.5)) );
                }
               float diff;
               for(int ip=0; ip<Waves[channel].nPt(); ++ip){
               diff= Waves[channel].Y[ip]-FltWaves[channel].Y[ip];
               bin_bin[channel]->Fill(diff); 
               }
               }*/ 
                             
//-------------------------Disegno delle forme d'onda filtrate, tramite funzione filtro con taglio-----------------------
					/*	bool saveFltWave=true; 
					if (saveFltWave && !isTrg) {

						waveFltDir->cd();

						tmpCvFlt.push_back( new TCanvas(Form("CvFlt-Ch%d_ev%d",channel,jentry),Form("tmpFltWave-Ch%d_ev%d",channel,jentry)) );
						tmpCvFlt.back()->cd();
						tmpFltWaves.push_back( new TGraph ( dim, &X[0], &FltWaves[channel].Y[0]) );                                
						tmpFltWaves.back()->GetXaxis()->SetTitle("time [ns]");
						tmpFltWaves.back()->SetTitle(Form("tmpFltWave-Ch%d_ev%d",channel,jentry));
						tmpFltWaves.back()->GetYaxis()->SetTitleOffset(1.4);
						tmpFltWaves.back()->GetYaxis()->SetTitle("Volt");
						tmpFltWaves.back()->Draw("AL");
						tmpCvFlt.back()->Write();
						theFile->cd("/");
					}*/
				
//------------------------Disegno di una sola componente del filtro------------------------
				/*	bool savenoise_1=true; //waveform di noise individuate dall'RC
					  if (savenoise_1 && !isTrg) {

						noise_RC->cd();

						tmpCvnoise_1.push_back( new TCanvas(Form("Cvnoise_1-Ch%d_ev%d",channel,jentry),Form("tmpNoise_1-Ch%d_ev%d",channel,jentry)) );
						tmpCvnoise_1.back()->cd();
						tmpnoise_1.push_back( new TGraph ( dim, &X[0], &noise_1[channel].Y[0]) );                                
						tmpnoise_1.back()->GetXaxis()->SetTitle("time [ns]");
						tmpnoise_1.back()->SetTitle(Form("tmpNoise_RC-Ch%d_ev%d",channel,jentry));
						tmpnoise_1.back()->GetYaxis()->SetTitleOffset(1.4);
						tmpnoise_1.back()->GetYaxis()->SetTitle("Volt");
						tmpnoise_1.back()->Draw("AL");
						tmpCvnoise_1.back()->Write();
						theFile->cd("/");
					} */
          

//--------------------Disegno della wave di segnale: signal_1--------------------------------------------------------
          //creazione della mappa per la sovrapposizione dello smooth di SG
    	    m=23;
					k=2;
					std::vector<float> tmpWSG_23_signal=smoothSG(Waves_signal_1[channel].Y,m,k);
					SmtSGWaves_23[channel].fillWave(tmpWSG_23_signal);
          //////////////////////////////////////////////////////////////////////////////                                                  
					bool savesignal_1=true; //true
					if (savesignal_1 && !isTrg) {

						signal->cd();

						tmpCvsignal_1.push_back( new TCanvas(Form("CvSignal_1-Ch%d_ev%d",channel,jentry),Form("tmpSignal_1-Ch%d_ev%d",channel,jentry)) );
						tmpCvsignal_1.back()->cd();
						tmpsignal_1.push_back( new TGraph ( dim, &X[0], &Waves_signal_1[channel].Y[0]) );                                
						tmpsignal_1.back()->GetXaxis()->SetTitle("time [ns]");
						tmpsignal_1.back()->SetTitle(Form("tmpSignal_afterFlt-Ch%d_ev%d",channel,jentry));
						tmpsignal_1.back()->GetYaxis()->SetTitleOffset(1.4);
						tmpsignal_1.back()->GetYaxis()->SetTitle("Volt");
						tmpsignal_1.back()->Draw("AL");
           	m=23;
						k=2; 
						tmpsignal_1.push_back(new TGraph(dim-m,&X[(m-1)/2],&SmtSGWaves_23[channel].Y[0]));
						tmpsignal_1.back()->SetLineColor(kBlue);
						tmpsignal_1.back()->Draw("Lsame");                                                
						tmpCvsignal_1.back()->Write(); 
						theFile->cd("/");
					}
           
//------------------Disegno della wave di noise totale: noise_total---------------------------------------------
					bool savenoise_total=false;
					if (savenoise_total && !isTrg) {

						noise->cd();

						tmpCvNoise_total.push_back( new TCanvas(Form("CvNoise_total-Ch%d_ev%d",channel,jentry),Form("tmpNoise_total-Ch%d_ev%d",channel,jentry)) );
						tmpCvNoise_total.back()->cd();
						tmpNoise_total.push_back( new TGraph ( dim, &X[0], &Waves_noise_tot[channel].Y[0]) );                                
						tmpNoise_total.back()->GetXaxis()->SetTitle("time [ns]");
						tmpNoise_total.back()->SetTitle(Form("tmpNoise_TotalNoise-Ch%d_ev%d",channel,jentry));
						tmpNoise_total.back()->GetYaxis()->SetTitleOffset(1.4);
						tmpNoise_total.back()->GetYaxis()->SetTitle("Volt");
						tmpNoise_total.back()->Draw("AL");
						tmpCvNoise_total.back()->Write();
						theFile->cd("/");
					}

					//--------------------------------ALGORITMO PER TROVARE FORME D'ONDA PIENE E FORME D'ONDA VUOTE.
				/*	if ( Smt10Waves[channel].nMax()>( waveCut[channel].first + nSig*waveCut[channel].second ) ) {
						isFull[channel]=true;
						std::cout<<"Ev "<<jentry<<" Ch "<<channel<<" Max "<<Smt10Waves[channel].nMax()<<" cut "<<( waveCut[channel].first + nSig*waveCut[channel].second )<<" bsl "<<Smt10Waves[channel].bsln<<" max "<<Smt10Waves[channel].max<<std::endl;
						if (nSelEv.find(channel)==nSelEv.end()) {
							nSelEv[channel]=1;
						} else { ++nSelEv[channel]; }
					}*/
//------------------------Pulizia delle variabili--------------------

					delete [] realFltFFT;
					delete [] imgFltFFT;
				   //delete [] amplFltFFT;
				  //delete [] amplFltFFT_1;
					delete [] amplFltFFT_2;
       //delete [] amplFltFFT_3;

				}// if point

			} // end of channel loop
       
		}// auto point
    
		/*for ( auto chsDiff : chToDiff ) {
			doChDiff(jentry,chsDiff,isFull,diffWaves,Waves,Smt4Waves,4,hDiff,hDiffPed,theFile);
		}*/

		if (((wave)Waves_signal_1[3]).n1IntegInR()<=0) {
			std::cout<<"Strange Intergal ch 3, ev "<<jentry<<std::endl;
		}
		std::cout<<"ch 3, InNInR "<<((wave)Waves_signal_1[3]).n1IntegInR()<<" IntInR "<<((wave)Waves_signal_1[3]).integInR<<" bsl "<<((wave)Waves_signal_1[3]).bsln<<std::endl;
		std::cout<<"InN "<<((wave)Waves_signal_1[3]).n1Integ()<<" Int "<<((wave)Waves_signal_1[3]).integ<<std::endl;
		std::cout<<"reg A "<<((wave)Waves_signal_1[3]).RegA()<<" reg B "<<((wave)Waves_signal_1[3]).RegB()<<std::endl;

		hIntRatio->Fill(((wave)Waves_signal_1[3]).n1IntegInR()/((wave)Waves_signal_1[6]).n1IntegInR()); //FIX cambiare la selezione dei canali!
		hMaxRatio->Fill(((wave)Waves_signal_1[3]).n1MaxInR()/((wave)Waves_signal_1[6]).n1MaxInR()); //FIX cambiare la selezione dei canali!
   
	} // end of loop on entry (second loop)

	//bin_bin->Write();

	/************************************Isto su rms e max*************************************/
//	float *rms_m = new float[nMaxCh];
//	float *max_m = new float[nMaxCh];
//	float *x_channel=new float[nMaxCh];
	std::vector<float> rms_m;
	std::vector<float> max_m;
	std::vector<float> maxCut_m;
	std::vector<float> x_channel;
	for (int channel=0; channel<=nMaxCh; channel++){
//		x_channel[channel]= channel;
      
		bool isTrg=false;
		for(auto trgCh : trigCh) { if (channel==trgCh) { isTrg=true; } }
		if (!isTrg) {
			if ( HstPerCh.find(channel)!=HstPerCh.end() ) {
				x_channel.push_back(channel);
				/*rms_m[channel] =*/  rms_m.push_back( ((hstPerCh*)HstPerCh[channel])->hRms->GetMean() );
				max_m.push_back( ((hstPerCh*)HstPerCh[channel])->hMaxVNInR->GetMean() );
				((hstPerCh*)HstPerCh[channel])->hMaxVNInR->GetXaxis()->SetRangeUser(0.008,0.1);
				/*max_m [channel] =*/ maxCut_m.push_back( ((hstPerCh*)HstPerCh[channel])->hMaxVNInR->GetMean() );
			} /*else{
//				rms_m[channel]=0;
//				max_m [channel]=0;
			}*/
		} /*else{
//			rms_m[channel]=0;
//			max_m [channel]=0;
		}*/
	}

	TCanvas *cvStatSum = new TCanvas("cvStatSum","ChStatSum");
	cvStatSum->Divide(1,3);
//	TGraph* rms_d= new TGraph(nMaxCh,x_channel,rms_m);
//	TGraph* max_d= new TGraph(nMaxCh,x_channel,max_m);
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
				std::cout<<channel<<"\t"<<nSelEv[channel]<<std::endl;
			}
		}
	}

	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	//                                               //  Write Histograms //                                                    //

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	theFile->cd();
	theFile->Write();
	theFile->Close();

	//   for (int channel=0; channel<=33; channel++){
	//	   delete HstPerCh[channel];
	//   }

}


/*void read_data::doChDiff(int evN, std::pair<int,int> &chDiff, std::map<int, bool> &isFull, diffWvCont &diffWaves, WvCont &Waves, WvCont &Smt4Waves, int SF, hDiffCont &hDiff, hDiffCont &hDiffPed, TFile *theFile){

	bool diffEv=false;
	int frstch=chDiff.first;
	int scndCh=chDiff.second;
	if (isFull[frstch] && !isFull[scndCh]) {
		diffEv=true;
		diffWaves[chDiff].clear();
		for (int ipt=0; ipt<(dim-2*SF); ++ipt) {
			float tmpVal = Waves[frstch].Y[ipt+SF]-Smt4Waves[scndCh].Y[ipt];
			diffWaves[chDiff].addPnt(tmpVal);
		}
	}
	else if (!isFull[frstch] && isFull[scndCh]) {
		diffEv=true;
		diffWaves[chDiff].clear();
		for (int ipt=0; ipt<(dim-2*SF); ++ipt) {
			float tmpVal = Waves[scndCh].Y[ipt+SF]-Smt4Waves[frstch].Y[ipt];
			diffWaves[chDiff].addPnt(tmpVal);
		}
	}

	if (diffEv) {
		hDiff[chDiff].hBsl->Fill(diffWaves[chDiff].bsln);
		hDiff[chDiff].hInteg->Fill(diffWaves[chDiff].integ);
		hDiff[chDiff].hIntegN->Fill(diffWaves[chDiff].nInteg());
		hDiff[chDiff].hRms->Fill(diffWaves[chDiff].rms);
		hDiff[chDiff].hMaxV->Fill(diffWaves[chDiff].max);
		hDiff[chDiff].hMaxVN->Fill(diffWaves[chDiff].nMax());

		bool saveWave=false;
		if (saveWave) {

			//			 waveDiffDir->cd();
			theFile->cd("WavesDiff");

			tmpCv.push_back( new TCanvas(Form("Cv-Ch%d-Ch%d_ev%d",frstch,scndCh,evN),Form("tmpWave-Ch%d-Ch%d_ev%d",frstch,scndCh,evN)) );
			tmpCv.back()->cd();
			tmpWaves.push_back( new TGraph ( dim-2*SF, &X[SF], &diffWaves[chDiff].Y[0]) );
			tmpWaves.back()->GetXaxis()->SetTitle("time [ns]");
			tmpWaves.back()->SetTitle(Form("tmpWave-Ch%d-Ch%d_ev%d",frstch,scndCh,evN));
			tmpWaves.back()->GetYaxis()->SetTitleOffset(1.4);
			tmpWaves.back()->GetYaxis()->SetTitle("Volt");
			tmpWaves.back()->Draw("AL");
			//tmpWaves.push_back( new TGraph ( dim-2*SF, &X[SF], &Smt10Waves[channel].Y[0]) );
			//tmpWaves.back()->SetLineColor(kBlue);
			//tmpWaves.back()->Draw("Lsame");
			tmpCv.back()->Write();
			theFile->cd("/");
		}

	}

	if (!isFull[frstch] && !isFull[scndCh]) {
		diffWaves[chDiff].clear();
		for (int ipt=0; ipt<(dim-2*SF); ++ipt) {
			float tmpVal = Waves[frstch].Y[ipt+SF]-Smt4Waves[scndCh].Y[ipt];
			diffWaves[chDiff].addPnt(tmpVal);
		}
		hDiffPed[chDiff].hBsl->Fill(diffWaves[chDiff].bsln);
		hDiffPed[chDiff].hInteg->Fill(diffWaves[chDiff].integ);
		hDiffPed[chDiff].hIntegN->Fill(diffWaves[chDiff].nInteg());
		hDiffPed[chDiff].hRms->Fill(diffWaves[chDiff].rms);
		hDiffPed[chDiff].hMaxV->Fill(diffWaves[chDiff].max);
		hDiffPed[chDiff].hMaxVN->Fill(diffWaves[chDiff].nMax());
	}

}*/
/*----------------sottrazione waveform-waveform filtrate---------------------------------*/
/*void read_data::doChDiff_waves(int evN, int channel, WvCont & Waves_diff, WvCont &Waves, WvCont FltWaves, TFile *theFile){
	bool diffEv=true;
		diffEv=true;
		Waves_diff[channel].clear();
		for (int ipt=0; ipt<(dim); ++ipt) {
			float tmpVal = Waves[channel].Y[ipt]-FltWaves[channel].Y[ipt];
			Waves_diff[channel].addPnt(tmpVal);
		}
		bool saveWave=true;
		if (saveWave) {
			theFile->cd("WavesDiff");
			tmpCv.push_back( new TCanvas(Form("Cv-Ch%d-Ch%d_ev%d",channel,channel,evN),Form("tmpWave-Ch%d-Ch%d_ev%d",channel,channel,evN)) );
			tmpCv.back()->cd();
			tmpWaves.push_back( new TGraph ( dim, &X[0], &Waves_diff[channel].Y[0]) );
			tmpWaves.back()->GetXaxis()->SetTitle("time [ns]");
			tmpWaves.back()->SetTitle(Form("tmpWave-Ch%d_waves-Ch%d_FltWaves_ev%d",channel,channel,evN));
			tmpWaves.back()->GetYaxis()->SetTitleOffset(1.4);
			tmpWaves.back()->GetYaxis()->SetTitle("Volt");
			tmpWaves.back()->Draw("AL");
			tmpCv.back()->Write();
			theFile->cd("/");
		}
  }*/


