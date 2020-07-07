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

Float_t omegaRefs[4]={296e+6, 383e+6, 443e+6, 485e+6};
Float_t integLoss[4]={0.993954, 0.996336, 0.997252, 0.997618};

//void doChDiff(int evN, std::pair<int,int> &chDiff, std::map<int, bool> &isFull, diffWvCont &diffWaves, WvCont &Waves, WvCont &Smt4Waves, int SF, hDiffCont &hDiff, hDiffCont &hDiffPed, TFile *theFile);
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////MAIN LOOP function//////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void read_data::Loop(Char_t *output, Int_t MidEv,Int_t eventn,  Bool_t evalWaveCut,TString fOutName)
{




	tmax=_tmax;  //to fix compilation problem
	double timeRes;
#ifndef _OSC
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
	//std::map<int, hstPerCh *> HstSF10PerCh;
	//   theFile->cd("/");
	//   for (int channel=0; channel<=33; channel++){
	//	   TDirectory *chDir = theFile->mkdir(Form("H-Ch%d",channel));
	//	   chDir->cd();
	//	   HstPerCh.insert(make_pair(channel,new hstPerCh(channel)));
	//	   theFile->cd("/");
	//   }

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
	TDirectory *waveFltDir = theFile->mkdir("WavesFlt");//reham 
	theFile->cd("/");//reham
	TDirectory *waveFltPeakDir = theFile->mkdir("Amplfft");    //directory per lo spettro in frequenza dell'ampiezza con smooth
	theFile->cd("/");
	TDirectory *waveFltPeakDir_1 = theFile->mkdir("Amplfft-LOWFreq");    //directory per lo spettro in frequenza dell'ampiezza con smooth(0-50MHz)
	theFile->cd("/");
	//	TDirectory *noise_RC = theFile->mkdir("noise_RC");         //directory per la waveform di noise individuata dal filtro RC(noise_1,singolo filtro)
	//	theFile->cd("/");
	TDirectory *flt_battimento = theFile->mkdir("Flt_battimento"); //directory per la waveform con noise di battimento + fit
	theFile->cd("/");
	TDirectory *flt_battimento2 = theFile->mkdir("Flt_battimento2"); //directory per la waveform con noise di battimento: seconda frequenza di noise
	theFile->cd("/");
	TDirectory *signal = theFile->mkdir("signal_Afterflt"); //directory per la waveform di segnale dopo tutti i filtri
	theFile->cd("/");
	TDirectory *noise = theFile->mkdir("noise_total");      //directory per il noise totale
	theFile->cd("/");
	TDirectory *NegWaveDir = theFile->mkdir("Wave_negative_batt"); //directory per la wave negativa
	theFile->cd("/");
	TDirectory *After_MidNotch = theFile->mkdir("Wave after Notch at middle frequency "); //directory per la wave negativa
	theFile->cd("/");

	//	TDirectory *FilterblsDir = theFile->mkdir("Filterbls");//directory per la differenza bin per bin
	//	theFile->cd("/");

	//	sovrapposizione onda di battimento:funzione di fit
	TF1 *pbatt=new TF1("pbatt","-2*[0]*TMath::Cos(0.5*TMath::TwoPi()*([2]-[1])*(x-[3])*1e-9)*TMath::Sin(0.5*TMath::TwoPi()*([2]+[1])*(x-[3])*1e-9)",0,1000);
	pbatt->SetParLimits(0,0,0.1);
	//  pbatt->SetParLimits(1,10e+6,16e+6);
	pbatt->SetParLimits(1,10e+6,16e+6);
	pbatt->SetParLimits(2,10e+6,16e+6);
	pbatt->SetParLimits(3,0,1000);
	pbatt->SetNpx(1000);

	//	istogramma battimenti
	TH1F *f_batt= new TH1F("f_batt","beat frequency", 1200,-3e+6,3e+6);
	f_batt->GetYaxis()->SetTitle("Entries");
	f_batt->GetXaxis()->SetTitle("Beat frequency [Hz]");
	TH1F *f_mean= new TH1F("f_mean","mean frequency", 1000,1e+6,20e+6);
	TH1F *t_0= new TH1F("t_0","phase", 1000,0,1000);
	TSpectrum *spectrum = new TSpectrum(10,1);
	//istogramma distribuzione dei picchi
	//TH1F *peak_HF= new TH1F("peak_HF","High frequency peak distribution", 100,3.2e+9,3.7e+9);
	/*peak_HF->GetXaxis()->SetTitle("Peak position (#omega [Hz])");
  	   	  peak_HF->GetYaxis()->SetTitle("Entries");*/

	TH1F *peak_LF= new TH1F("peak_LF","Low frequency peak distribution",100, 2.0e+7,9.0e+8);
	peak_LF->GetXaxis()->SetTitle("Peak position (#omega [Hz])");
	peak_LF->GetYaxis()->SetTitle("Entries");
	TH1F *nPeak_LF= new TH1F("nPeak_LF","Low frequency n. of peaks",10, 0,10);
	nPeak_LF->GetXaxis()->SetTitle("n. Peaks");
	nPeak_LF->GetYaxis()->SetTitle("Entries");

	TH1F *peak_VLF= new TH1F("peak_VLF","Very low frequency peak distribution(<100MHz)",100, 2.0e+7,9.0e+8);
	peak_VLF->GetXaxis()->SetTitle("Peak position (#omega [Hz])");
	peak_VLF->GetYaxis()->SetTitle("Entries");
	TH1F *nPeak_VLF= new TH1F("nPeak_VLF","Very low frequency n. of peaks (<100MHz)",10, 0,10);
	nPeak_VLF->GetXaxis()->SetTitle("n. Peaks");
	nPeak_VLF->GetYaxis()->SetTitle("Entries");


	/*Istogramma di chi quadro e ampiezza del fit del battimento*/
	TH1F *chi2= new TH1F("chi2","Chi2 distribution", 100000,-1e-5,1e-2);
	TH1F *ampl_battimento= new TH1F("ampl_battimento","Ampl distribution", 1000,0,1e-2);
	TH1F *ampl_battimento2= new TH1F("ampl_battimento2","Ampl distribution", 1000,0,1e-2);
	TH1F *f_battcut= new TH1F("f_battcut","beat frequency after cut on min value", 1000,-2e+6,4e+6);
	chi2->GetXaxis()->SetTitle("chi quadro ");
	chi2->GetYaxis()->SetTitle("Entries");
	ampl_battimento->GetXaxis()->SetTitle("ampiezza");
	ampl_battimento->GetYaxis()->SetTitle("Entries");
	ampl_battimento2->GetXaxis()->SetTitle("ampiezza");
	ampl_battimento2->GetYaxis()->SetTitle("Entries");
	f_battcut->GetXaxis()->SetTitle("beat ampl");
	f_battcut->GetYaxis()->SetTitle("Entries");
	TMinuit *gMinuit= new TMinuit();
	//TVirtualFitter::SetDefaultFitter("Migrad");

	/*Istogramma distribuzione di min value e min pos */
	TH1F *hmin_value= new TH1F("hmin_value","Min value distribution", 100,0,-0.01e-4);
	hmin_value->GetXaxis()->SetTitle("Min_value");
	hmin_value->GetYaxis()->SetTitle("Entries");
	TH1F *hmin_pos= new TH1F("hmin_pos","Position min value distribution", 1000,0,1000);
	hmin_pos->GetXaxis()->SetTitle("Position min value");
	hmin_pos->GetYaxis()->SetTitle("Entries");


	//TString fOutName="Full_Wave.txt";
	ofstream myfile;
	if (!fOutName.IsNull()) {
		myfile.open (fOutName.Data());
	}


	  //}

	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	//                                               //  Loop on entries //                                                  //

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	std::vector<int> trigCh; //vettore di interi per i canali di trigger.
	trigCh.clear();
#ifndef _OSC
	trigCh.push_back(32);
	trigCh.push_back(33);
#else //for muon data
	trigCh.push_back(7);//osc1
	trigCh.push_back(8);//osc2
#endif

	std::map<int, std::pair<float, float> > waveCut;
	waveCut.clear();
	int nSig=3; // numero di sigma
	std::map<int, int> nSelEv;
	nSelEv.clear();

	if (fChain == 0) return; //è il puntatore TTree

	Long64_t nentries = fChain->GetEntriesFast(); //Return the number of entries as of the last check.
	//cout << "Number of entries in the tree for real data is= " << nentries << endl;
	//  Long64_t nbytes = 0, nb = 0;

	Int_t frstEv=0;                               //inizializzo il primo evento a zero
	Int_t lastEv=nentries;                        //l'ultimo evento coincide con le nentries, che si prende con la funzione get...sopra
	//Int_t MidEv;
	if (eventn>=0&&eventn<lastEv) {
		frstEv=eventn;
		lastEv=eventn+1;
	}else if(MidEv>=0&&eventn<-1 ){
		frstEv=MidEv;
		lastEv=-eventn+frstEv;
	}else if (eventn<-1 && eventn>-nentries){
		lastEv=-eventn;}
	//cout<< " MidEv "<< MidEv<< "lastEv "<<lastEv<<endl;



	//bool evalWaveCut=false;
	TString wCutFName = Form("wave-Cut-%s.txt",out.Data());
	std::ofstream cutf;
	if (evalWaveCut) {

		std::ofstream cutf;
		cutf.open(wCutFName.Data(), std::ofstream::out);//stampa su file. 

		int SF=10;

		//	   std::map<int, TF1*> tmpFitRes;

		for (Long64_t jentry=frstEv; jentry<lastEv;jentry++) {
			Long64_t ientry = LoadTree(jentry); //The function finds the corresponding Tree and returns the entry number in this tree.  
			if (ientry < 0) break;
			//	        nb = fChain->GetEntry(jentry);   nbytes += nb;
			fChain->GetEntry(jentry);
			WvCont Waves;
			Waves.clear();
			WvCont Smt10Waves;
			Smt10Waves.clear();
			//lettura dei file di dati .root
#ifndef _OSC
			if(jentry==frstEv) {
				// 1024 channel of ADC = 1 micro sec => every channel is 0.9766 nano sec
				timeRes = (_tmax*1.0e+9)/((float)dim);//in nano sec
				X.clear();
				for (int i = 0; i < dim; ++i) { X.push_back(timeRes *(i+1)); }
			}
			for (auto point : wd->getX742Data()) {
				for (int channel=0; channel<=nMaxCh+2; channel++){
#else
					if(jentry==frstEv) {
						nMaxCh=wd->getRunHeader()->getNumberOfDevicesInRun()*8-1;//i primi canali sono dell'oscilloscopio 1, gli altri sono dell'oscilloscopio 2.
						dim=wd->getRunHeader()->getNumberOfPointsPerChannel()-1;//-1 fix the bug of the last point //needed only for run < 41
						_tmax=(dim+1)*wd->getRunHeader()->getTimeResolOfChannel();
						timeRes = (_tmax*1.0e+9)/((float)dim);//in nano sec
						X.clear();
						//std::cout<<"dim "<<dim<<" tmax "<<_tmax<<std::endl;
						for (int i = 0; i < dim; ++i) { X.push_back(timeRes *(i+1)); }
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
								Waves[channel].fillWave(point.second,dim);
								std::vector<float> tmpWSF10=smooth(Waves[channel].Y,SF);
								Smt10Waves[channel].fillWave( tmpWSF10 );

								/*((hstPerCh*)HstSF10PerCh[channel])->hBsl->Fill(((wave)Smt10Waves[channel]).bsln);
						((hstPerCh*)HstSF10PerCh[channel])->hInteg->Fill(((wave)Smt10Waves[channel]).integ);
						((hstPerCh*)HstSF10PerCh[channel])->hIntegN->Fill(((wave)Smt10Waves[channel]).nInteg());
						((hstPerCh*)HstSF10PerCh[channel])->hRms->Fill(((wave)Smt10Waves[channel]).rms);
						((hstPerCh*)HstSF10PerCh[channel])->hMaxV->Fill(((wave)Smt10Waves[channel]).max);
						((hstPerCh*)HstSF10PerCh[channel])->hMaxVN->Fill(((wave)Smt10Waves[channel]).nMax());*/

								bool saveWave=false; //non salvo, quindi non vediamo l'out di questo pezzo di codice.
								//						bool saveWave=true;
								if (saveWave) {
									//if (channel==0 && (Smt10Waves[channel].nMax()>0.0035 && Smt10Waves[channel].nMax()<0.005) ) {
									if (channel==11 && (Smt10Waves[channel].nMax()>0.0125/*(0.003501+5*0.0009636)*/) ) {
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
				waveCut.insert( make_pair(channel, make_pair(tmpMean,tmpSigma) ) );   //riempiamo la mapa waveCut definita all'inizio
			}
		}*/
				//std::cout<<"Fit of Intergeal SWF 10 fit results:"<<std::endl;
				for (int channel=0; channel<=nMaxCh; channel++){
					if (waveCut.find(channel)!=waveCut.end()) {              //fintanto che i canali non sono finiti scrivi su file.
						//std::cout<<channel<<"\t"<<waveCut[channel].first<<"\t"<<waveCut[channel].second<<std::endl;
						//scrittura, first è la media, second è la sigma, per ogni canale.
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

			/************Ricerca dei picchi di noise sulle FFT e plot dei risultati*************/


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
			WvCont FltMidNotch;
			WvCont FltWavesSM;
			WvCont noise_1;
			WvCont flt_batt;
			WvCont flt_batt2;
			WvCont Waves_signal_1;
			WvCont Waves_noise_tot;
			WvCont Smt10Waves;
			WvCont Smt4Waves;
			WvCont tmpWSG_signal;
//			diffWvCont diffWaves;
			fftCont Wffts;
			fftCont CR;
			fftCont RC;
			fftCont Notch_1;
			fftCont Notch_2;
			fftCont Notch_3;
			fftCont Notch_4;
			fftCont Notch_5;
//			fftCont ContforNotch;
//-------------------------------------------creazione mappa per smooth di SG------------------------------/
			WvCont SmtSGWaves;
			WvCont SmtSGWaves_7;
			WvCont SmtSGWaves_9;
			WvCont SmtSGWaves_13;
			WvCont SmtSGWaves_23;
			WvCont SmtSGWaves_4_9;
			WvCont SmtSGWaves_4_7;
			WvCont smtSG;
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



//---------------------------------------------------------------SECONDO loop, da qui comincia la vera analisi.----------------------------------------
			for (Long64_t jentry=frstEv; jentry<lastEv;jentry++) {

				Long64_t ientry = LoadTree(jentry);
				if (ientry < 0) break;

//				nb = fChain->GetEntry(jentry);   nbytes += nb;
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
				FltMidNotch.clear();
				noise_1.clear();
				flt_batt.clear();
				flt_batt2.clear();
				SmtSGWaves.clear();
				SmtSGWaves_7.clear();
				SmtSGWaves_9.clear();
				SmtSGWaves_13.clear();
				SmtSGWaves_4_9.clear();
				SmtSGWaves_4_7.clear();
				SmtSGWaves_23.clear();
				Waves_signal_1.clear();
				Waves_noise_tot.clear();
				tmpWSG_signal.clear();

//				for ( auto chsDiff : chToDiff ) { diffWaves[chsDiff].clear(); }

#ifndef _OSC
				if(jentry==frstEv) {
					// 1024 channel of ADC = 1 micro sec => every channel is 0.9766 nano sec
					timeRes = (_tmax*1.0e+9)/((float)dim);//in nano sec
					X.clear();
					for (int i = 0; i < dim; ++i) { X.push_back(timeRes *(i+1)); }
				}
				for (auto point : wd->getX742Data()) {
					//cout<<"Discrimentor "<<endl;
					//cout<<"**********************"<<endl;

					for (int channel=0; channel<=nMaxCh+2; channel++){
#else
						if(jentry==frstEv) {
							nMaxCh=wd->getRunHeader()->getNumberOfDevicesInRun()*8-1;
							dim=wd->getRunHeader()->getNumberOfPointsPerChannel()-1;//-1 fix the bug of the last point //needed only for run < 41
							_tmax=(dim+1)*wd->getRunHeader()->getTimeResolOfChannel();
							tmax=_tmax;
							//std::cout<<"MaxCh "<<nMaxCh<<" dim "<<dim<<" tmax "<<_tmax<<std::endl;
							timeRes = (_tmax*1.0e+9)/((float)dim);//in nano sec
							X.clear();
							for (int i = 0; i < dim; ++i) { X.push_back(timeRes *(i+1)); }
						}
						for (auto point : wd->getXOSCData()) {
							for (int channel=0; channel<=nMaxCh; channel++){

#endif

								bool isTrg=false; //serve a controllare che il canale su cui facciamo l'analisi non sia quello di trigger.
								for(auto trgCh : trigCh) { if (channel==trgCh) { isTrg=true; } }
								if (point.first == channel  ) {
									//cout << "Channel is= " << channel << endl;
//									Waves[channel].fillWave(point.second.size(),point.second.data());
									Waves[channel].fillWave(point.second,dim);
									//std::cout<<"nPt "<<Waves[channel].nPt()<<std::endl;

/****-------------------------------TRASFORMATA DI FOURIER---------------------------------*****/
									FFT(Waves[channel],Wffts[channel]);//trasformata di Fourier
									double *realFltFFT = new double[Waves[channel].nPt()];
									double *imgFltFFT = new double[Waves[channel].nPt()];
/****-------------------------------FUNZIONI DI FILTRO-------------------------------------*****/
//									filterWave(Wffts[channel],realFltFFT,imgFltFFT,3e+6,5e+6);  //funzione filtro(taglio)
//									 filterWaveSM(Wffts[channel],realFltFFT,imgFltFFT,9,2);     //smooth su SG
//									InverseFFT(realFltFFT,imgFltFFT,Waves[channel].nPt(),FltWavesSM[channel]);//trasformata inversa, con singolo filtro.
//									double *amplFltFFT = new double[Waves[channel].nPt()];
//									for (int ip=0; ip<Waves[channel].nPt()/2; ++ip) {
//									std::cout<<"omega "<<Wffts[channel].omega[ip]<<" Flt: real "<<realFltFFT[ip]<<" img "<<imgFltFFT[ip]<<std::endl;
//									amplFltFFT[ip]=TMath::Sqrt(realFltFFT[ip]*realFltFFT[ip]+imgFltFFT[ip]*imgFltFFT[ip]);
//										}
//										for (int ip=Waves[channel].nPt()/2; ip<Waves[channel].nPt(); ++ip) {
//											amplFltFFT[ip]=0.0;
//										}

/***------------------------------- RIEMPIMENTO grafici PER L'AMPIEZZA CON SMOOTH(cerchiamo i picchi di noise sulle forme d'onda con smooth di sg)---*****/


									//ricerca picchi a freq <100MHz
									bool saveAmpl_1=false;
									Double_t *peaks_VLF=0x0;//VLF=very low freq
									int numOfPeaks_VLF=0;
									if(! isTrg){
										waveFltPeakDir_1->cd();
										if(saveAmpl_1){
											//************RIEMPIMENTO CANVAS*********
											tmpCvPeak_1.push_back( new TCanvas(Form("CvPeak_1-Ch%d_ev%d",channel,jentry),Form("tmpAmplfft_1-Ch%d_ev%d",channel,jentry)) );
											tmpFFTAmplPeak_1.push_back( new TH1F(Form("FFTampl-Ch%d_ev%d",channel,jentry),";#omega [Hz];Amplitude",
													dim/2,Wffts[channel].omega[0],Wffts[channel].omega[dim/2-1]));
											// search peaks///
											std::vector<float> fftAmplFrNtch_1;
											fftAmplFrNtch_1.clear();
											std::vector<int> fftIdxFrNtch_1;
											fftIdxFrNtch_1.clear();

											for(int i=0; i < (dim/2); ++i) {
												if(Wffts[channel].omega[i]> 0.0&&Wffts[channel].omega[i]<70e+6){
													fftAmplFrNtch_1.push_back(Wffts[channel].ampl[i]);
													fftIdxFrNtch_1.push_back(i);
												}
											}

											//std::vector<float> fftAmplFrNtch_smt_1=smoothSG(fftAmplFrNtch_1,5,2);//for fedrica
											for (int i=1; i<fftIdxFrNtch_1.size(); ++i) { //2=(int)5/2 i=2 if smoth 5 and i =3 if smooth 7 (3,central,3)
												/*if (fftAmplFrNtch_smt_1[i]>0.0 && fftAmplFrNtch_smt_1[i]<100.0) */tmpFFTAmplPeak_1.back()->SetBinContent(fftIdxFrNtch_1[i],fftAmplFrNtch_1[i]);
											}
											tmpFFTAmplPeak_1.back()->Draw();



										numOfPeaks_VLF =  spectrum -> Search(tmpFFTAmplPeak_1.back(),1,"nobackground",0.05);
										//cout << "Peaks found for VLF" << numOfPeaks_VLF << " evento "<< jentry << " canale " << channel <<endl ; ;
										peaks_VLF = spectrum->GetPositionX() ;
										if (numOfPeaks_VLF>0) { nPeak_VLF->Fill(numOfPeaks_VLF); }
										for (int j=0 ; j<numOfPeaks_VLF ; j++){

											//cout<<"peaks found <100MHz "<<peaks_VLF[j]<< " evento "<< jentry << " canale " << channel <<endl ;
											peak_VLF->Fill(peaks_VLF[j]);
										}
									}///lalala
										//*****************************************************************************************************//
										if (peaks_VLF!=0x0){ //notch picchi a a bassa frequenza(<100MHz)
											wavefft tmpWfft2;
											tmpWfft2.fillAP(Wffts[channel].nPt,Wffts[channel].ampl,Wffts[channel].phi,Wffts[channel].omega);
											for(int i=0; i<numOfPeaks_VLF; ++i){
												//cout << "Peak at very low freq = " << peaks_VLF[i] << " evento "<< jentry << " canale " << channel <<endl ;
												filterNotch(Wffts[channel],realFltFFT,imgFltFFT,peaks_VLF[i],TMath::TwoPi()*5e+6,2);//10
												tmpWfft2.clear();
												tmpWfft2.fillRM(Wffts[channel].nPt,realFltFFT,imgFltFFT,Wffts[channel].omega);
											}
											Notch_1[channel].fillAP(tmpWfft2.nPt,tmpWfft2.ampl,tmpWfft2.phi,tmpWfft2.omega);
										}
									}

//									/*--------------Ricerca dei picchi a "medie frequenze!" (>200MHz)*/
									Double_t *peaks=0x0 ;
									int numOfPeaks=0;
//									//									double avg_frst=0.0;
//									//									double avg_lst=0.0;
									bool saveAmplfft = true;
									//									float maxLwFnFreq=0;
									if(! isTrg){
										waveFltPeakDir->cd();
										if(saveAmplfft){
											//************RIEMPIMENTO CANVAS*********
											tmpCvPeak.push_back( new TCanvas(Form("CvPeak-Ch%d_ev%d",channel,jentry),Form("tmpAmplfft-Ch%d_ev%d",channel,jentry)) );
											//tmpCvPeak.back()->cd()->SetLogy();
											//} si chiudeva qui il saveamplfft
											//************RIEMPIMENTO GRAFICI*****
											tmpFFTAmplPeak.push_back( new TH1F(Form("FFTampl-Ch%d_ev%d",channel,jentry),";#omega [Hz];Amplitude",
													dim/2,Wffts[channel].omega[0],Wffts[channel].omega[dim/2-1]));
											// search peaks///
											std::vector<float> fftAmplFrNtch;
											fftAmplFrNtch.clear();
											std::vector<int> fftIdxFrNtch;
											fftIdxFrNtch.clear();
											int nPntAvgLst=0;
											for(int i=1; i < (dim/2); ++i) {
//												if(Wffts[channel].omega[i]>0.0e+6&&Wffts[channel].omega[i]<80e+6){//for me(50e+6)
//													if(i<3) {
//														avg_frst= avg_frst + Wffts[channel].ampl[i];
//														if (Wffts[channel].ampl[i]>maxLwFnFreq){maxLwFnFreq=Wffts[channel].ampl[i];}
//													}
//													else {
//														avg_lst = avg_lst + Wffts[channel].ampl[i];
//														++nPntAvgLst;
//													}
//												}

											if(Wffts[channel].omega[i]>2.5e+8&&Wffts[channel].omega[i]<5e+8){
												fftAmplFrNtch.push_back(Wffts[channel].ampl[i]);
												fftIdxFrNtch.push_back(i);
											}

										}
										//										avg_frst/=2.0;
										//										avg_lst/=(float)nPntAvgLst;

										std::vector<float> fftAmplFrNtch_smt=smoothSG(fftAmplFrNtch,5,2);//for fedrica
										//										std::vector<float> fftAmplFrNtch_smt=smooth(fftAmplFrNtch.data(),fftAmplFrNtch.size(),1);
										//										cout<<"fftAmplFrNtch_smt size = "<<fftAmplFrNtch_smt.size()<<endl;
										//										for (int ii=0; ii<fftAmplFrNtch_smt.size();ii++){cout<<"fftAmplFrNtch_smt = "<<fftAmplFrNtch_smt.at(ii)<<endl;}
										for (int i=0; i<fftIdxFrNtch.size(); ++i) { //2=(int)5/2 i=2 if smoth 5 and i =3 if smooth 7 (3,central,3)
											 tmpFFTAmplPeak.back()->SetBinContent(fftIdxFrNtch[i]+1,fftAmplFrNtch[i]);//fftAmplFrNtch_smt
										}
										tmpFFTAmplPeak.back()->Draw();
									}
									//**********************************Ricerca dei picchi nello spettro di frequenza: media freq*********************//

									numOfPeaks =  spectrum -> Search(tmpFFTAmplPeak.back(),1,"nobackground",0.95);//0.5
									//cout << "Peaks found " << numOfPeaks << " evento "<< jentry << " canale " << channel <<endl ; ;
									peaks = spectrum->GetPositionX() ;
									//										numOfPeaks=spectrum->GetNPeaks();
									if (numOfPeaks>0) { nPeak_LF->Fill(numOfPeaks); }
									for (int j=0 ; j<numOfPeaks ; j++){
										//std::cout << "Peak at = " << peaks[j] << " evento "<< jentry << " canale " << channel <<std::endl ;
										peak_LF->Fill(peaks[j]);
									}
									if(saveAmplfft){
										tmpCvPeak.back()->Write();
									}
									theFile->cd("/");
								}

//***********************************Filtraggio noise, tramite Notch //notch picchi a media frequenza(dopo 200 Mhz)**************************************//
									if (peaks!=0x0){
										wavefft tmpWfft;
										tmpWfft.fillAP(Wffts[channel].nPt,Wffts[channel].ampl,Wffts[channel].phi,Wffts[channel].omega);
										for(int i=0; i<numOfPeaks; ++i){
//											cout << "Peak at = " << peaks[i] << " evento "<< jentry << " canale " << channel <<endl ;
											filterNotch(tmpWfft,realFltFFT,imgFltFFT,peaks[i],TMath::TwoPi()*20e+6,2);
											//cout<<"peaks to notch  = "<<peaks[i]<< " evento "<< jentry << " canale " << channel <<std::endl ;
											tmpWfft.clear();
											tmpWfft.fillRM(Wffts[channel].nPt,realFltFFT,imgFltFFT,Wffts[channel].omega);
										}
										Notch_2[channel].fillAP(tmpWfft.nPt,tmpWfft.ampl,tmpWfft.phi,tmpWfft.omega);
									}
									double *amplFltFFT_3 = new double[Waves[channel].nPt()];
									for (int ip=0; ip<Waves[channel].nPt()/2; ++ip) {
										amplFltFFT_3[ip]=TMath::Sqrt(realFltFFT[ip]*realFltFFT[ip]+imgFltFFT[ip]*imgFltFFT[ip]);
									}
									for (int ip=Waves[channel].nPt()/2; ip<Waves[channel].nPt(); ++ip) {
										amplFltFFT_3[ip]=0.0;
									}


//***********************************REHAM:  notch a bassissime frequenze,cn taglio pesato sull'ampiezza del picco di noise++++
//									wavefft tmpWfft2; //filtro a bassa frequenza(<100MHz)
//									tmpWfft2.fillAP(Wffts[channel].nPt,Wffts[channel].ampl,Wffts[channel].phi,Wffts[channel].omega);
//									//condition for noise WF
////									if ((avg_frst/avg_lst)>2){ //filterNotch(Notch_1[channel],realFltFFT,imgFltFFT,(TMath::TwoPi()*3.5e+6),TMath::TwoPi()*5e+6,2);//notch frequency 3 MHZ
////										cout<<" maxLwFnFreq "<<maxLwFnFreq<<endl;
////										filterNotch(Notch_1[channel],realFltFFT,imgFltFFT,maxLwFnFreq,TMath::TwoPi()*10e+6,2);
//										filterNotch(Wffts[channel],realFltFFT,imgFltFFT,maxLwFnFreq,TMath::TwoPi()*10e+6,2);
//										tmpWfft2.clear();
//										tmpWfft2.fillRM(Wffts[channel].nPt,realFltFFT,imgFltFFT,Wffts[channel].omega);
////									}
//									Notch_2[channel].fillAP(tmpWfft2.nPt,tmpWfft2.ampl,tmpWfft2.phi,tmpWfft2.omega);

//*****************************************************************************************************************//


//									filterWaveBsl(Notch_2[channel],realFltFFT,imgFltFFT);//filtro sulla baseline.                                        ///
//									InverseFFT(realFltFFT,imgFltFFT,Waves[channel].nPt(),FltMidNotch[channel]);//trasformata inversa, dopo notch a medie frequenze
									filterWaveBsl(Notch_2[channel],realFltFFT,imgFltFFT);//filtro sulla baseline.                                        ///
									InverseFFT(realFltFFT,imgFltFFT,Waves[channel].nPt(),FltWaves[channel]);//trasformata inversa
									// band pass for noise frequence around 13 MHz
									filterNotch(Notch_2[channel],realFltFFT,imgFltFFT,(TMath::TwoPi()*14e+6),TMath::TwoPi()*20e+6,1);//1 passa banda.    ////
									InverseFFT(realFltFFT,imgFltFFT,Waves[channel].nPt(),flt_batt[channel]);//trasformata inversa necessaria per il battimento
//second trying
									filterNotch(Notch_2[channel],realFltFFT,imgFltFFT,(TMath::TwoPi()*45e+6),TMath::TwoPi()*10e+6,1);//1 passa banda.    ////
									InverseFFT(realFltFFT,imgFltFFT,Waves[channel].nPt(),flt_batt2[channel]);//trasformata inversa necessaria per il battimento

									//the waveform after just notch filter
									bool saveWave=false;//true qui
									if (saveWave) {
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
									}



//*************************************Riempimento dei vettori per la waveform di battimento negativa+ricerca del range********************
									if(!isTrg){
										Y_negative.clear();
										X_negative.clear();
										//float fitRange;
										float min_value=0.,min_pos=0.;
										float peak_temp=0;
										float delta = 1.25e+9/13e+6;//dipende dal digitalizzatore
										float thr = -1.15e-3-3*0.64e-3;
										float thr1 = -1.3e-3;//soglia messa a "mano" osservando i segnali
										int ipmin=0;

										for(int ip=0;ip<flt_batt[channel].nPt();++ip){
											if(flt_batt[channel].Y[ip]<0){
												Y_negative.push_back(flt_batt[channel].Y[ip]);
												X_negative.push_back(X[ip]);
												if(min_value>Y_negative.back()) {
													min_value=Y_negative.back();
													min_pos=X[ip]; ipmin=ip;
												}

											}
										}
										//cout<<"jentry "<< jentry<<" channel "<<channel<<" min value "<<min_value<<" min_pos " <<min_pos<<endl;
										int ipRmin=0;
										int ipRmax=flt_batt[channel].nPt()-1;
										float min_range,max_range;
										if(min_value<thr/*&&min_value>-0.01*/){//Reham
											int nPl=1;
											while (true) {
												int tip=ipmin-nPl*delta;
												if(tip<0) { ipRmin=0; break;}
												bool found=false;
												int ip1=tip-50;
												int ip2=tip+50;
												if(ip1<0) { ip1=0; }
												if(ip2>=(flt_batt[channel].nPt()-1)) { ip2=flt_batt[channel].nPt()-1; }
												for(int ip=ip1; ip<ip2; ++ip){
													if (flt_batt[channel].Y[ip]<thr1) {found=true; break;}
												}
												if (found) {++nPl; ipRmin=ip1;}
												else { break; }
											}
											nPl=1;
											while (true) {
												int tip=ipmin+nPl*delta;
												if(tip>(flt_batt[channel].nPt()-1)) { ipRmax=flt_batt[channel].nPt()-1; break;}
												bool found=false;
												int ip1=tip-50;
												int ip2=tip+50;
												if(ip1<0) { ip1=0; }
												if(ip2>=(flt_batt[channel].nPt()-1)) { ip2=flt_batt[channel].nPt()-1; }
												for(int ip=ip1;ip<ip2;++ip){
													if (flt_batt[channel].Y[ip]<thr1) {found=true; break;}
												}
												if (found) {++nPl; ipRmax=ip2;}
												else { break; }
											}
										}
										min_range=X[ipRmin];
										max_range=X[ipRmax];
										//cout << min_range << " min range " << max_range << " max range "<<channel << " ch " << jentry <<" entry "<<endl;
//										 fitRange = max_range-min_range;
										hmin_value->Fill(min_value);
										hmin_pos->Fill(min_pos);


										//*******************************************************************************************************

										/*estrapolazione del segnale */
										float middle = 0;//(max_range+min_range)*0.5;
										float denom = 0;
										for (int ip=0; ip<dim; ++ip) {
											if ( (X[ip]>min_range && X[ip]<max_range) && flt_batt[channel].Y[ip]<0.0) {
												middle+=flt_batt[channel].Y[ip]*X[ip];
												denom+=flt_batt[channel].Y[ip];
											}
										}
										middle/=denom;
										//cout <<  " middle " << middle << channel << " canale " << jentry << " evento " << endl;
										Xerror.clear();
										Yerror.clear();
										/*riempimento vettori di errori per x e per y*/
										for (int i =0; i<dim;++i){
											Xerror.push_back(0.0);
										}
										for(int i=0;i<flt_batt[channel].nPt();++i){
//											Yerror.push_back(0.55e-3);
											Yerror.push_back(0.);
										//cout<< "value Yneg"<<Y_negative[i]<<endl;
										}
//										for(int i =0;i< Y_negative.size();++i){
//											cout<<" Yvalue ="<<Y_negative[i]<<endl;
//										}

										tmpflt_batt.push_back( new TGraph ( dim, &X[0], &flt_batt[channel].Y[0]) );
										//tmpflt_batt.back()->GetYaxis()->SetRangeUser(-0.01,0.01);
//										tmpflt_batt.back()->Fit(pbatt,"","",X[0]+range*(0.1),X[0]+range*(0.9));
//										tmpflt_batt.back()->Fit(pbatt,"","",X[0]+range*(0.1),X[0]+range*(0.9));

										tmpNegBatt.push_back(new TGraphErrors(X_negative.size(),&X_negative[0],&Y_negative[0],&Xerror[0],&Yerror[0]));
										tmpNegBatt.back()->GetYaxis()->SetRangeUser(-0.008,0.0);

//										tmpNegBatt.back()->GetXaxis()->SetRangeUser(0,1000);//me
										//cout<< "MIN RANGE prima del fit = "<< min_range << " MAX RANGE "<< max_range <<" evento "<< jentry << " canale " << channel <<endl;

										float chi,F_batt=0.;
										bool fitCut=false;
										if(min_value<thr /*&& min_value>-0.01*/){//////////me
//											pbatt->SetParameters(0.002,11e+6,13.5e+6,middle);//Fed
											//cout<< "MIN RANGE nel fit = "<< min_range << " MAX RANGE "<< max_range <<" evento "<< jentry << " canale " << channel <<endl;
											pbatt->SetParameters(0.008,12e+6,15e+6,middle);//me
//											gMinuit->mnsimp();
											tmpNegBatt.back()->Fit("pbatt","W","",min_range,max_range);
//											if (fabs(pbatt->GetParameter(1)-1.0e+7)<1e3 || fabs(pbatt->GetParameter(1)-1.6e+7)<1e3) {pbatt->SetParameter(1,1.11e+7);}
//											@//if (fabs(pbatt->GetParameter(1)-1.0e+7)<1e3 || fabs(pbatt->GetParameter(1)-1.6e+7)<1e3) {pbatt->SetParameter(1,1.0e+7);}
//											if (fabs(pbatt->GetParameter(2)-1.0e+7)<1e3 || fabs(pbatt->GetParameter(2)-1.6e+7)<1e3) {pbatt->SetParameter(2,1.35e+7);}
											tmpNegBatt.back()->Fit("pbatt","M","",min_range,max_range);
											tmpNegBatt.back()->Fit("pbatt","G","",min_range,max_range);
//											if (fabs(pbatt->GetParameter(1)-1.0e+7)<1e3 || fabs(pbatt->GetParameter(1)-1.6e+7)<1e3) {pbatt->SetParameter(1,1.11e+7);}
//											if (fabs(pbatt->GetParameter(2)-1.0e+7)<1e3 || fabs(pbatt->GetParameter(2)-1.6e+7)<1e3) {pbatt->SetParameter(2,1.35e+7);}
//											gMinuit->Migrad();
											TFitResultPtr result = tmpNegBatt.back()->Fit(pbatt,"S","",min_range,max_range);
											f_mean->Fill(0.5*(result->Parameter(1)+result->Parameter(2)));
											F_batt=0.5*(result->Parameter(2)-result->Parameter(1));
											f_batt->Fill(F_batt);
											t_0->Fill(result->Parameter(3));
											int ndf= result->Ndf();
											bool fit_stat=result->IsValid();
//											cout << " ndf " << ndf << endl;
//											cout<< " fit stat " << result->IsValid() << " " << result->CovMatrixStatus()<< endl;
											//float chi=(pbatt->GetChisquare())/(pbatt->GetNDF());//chi quadro normalizzato
											chi=(result->Chi2())/ndf;
											//riempimento istogrammi per ampiezza e chi
											chi2->Fill(chi);
											if(chi<0){cout << chi << " chi_quadro " << " jentry "<<jentry<<" channel "<< channel<<endl;}
											float ampl_b=result->Parameter(0);
											ampl_battimento->Fill(ampl_b);
											//if ((F_batt<386.3-3*2693||F_batt>386.3+3*2693)){ ampl_battimento2->Fill(ampl_b);}
											//cout << ampl_b << " p0 " << endl;
//											if(min_value>-0.007&& min_value<thr && ampl_b>=(3.6e-4+3*1.2e-4)){
//												f_battcut->Fill(F_batt);
//												//if(F_batt>2990e+3){cout << " canale " << channel << " evento sospetto " << jentry << endl;}
//											}
											fitCut=fit_stat && (F_batt<386.3-3*2693||F_batt>386.3+3*2693)&&((TMath::Abs(2*F_batt))<6e+6)&&ampl_b>=(3.6e-4+3*1.2e-4)/*&&chi<(1.91e-7+3*0.44e-7)*/;

											//if(!fitCut){cout << " evento non sottratto " << jentry << " canale " << channel << endl;}
										}
										for (int ip=0; ip<dim; ++ip) {
											float signal_1=FltWaves[channel].Y[ip];
											float noise_tot=Waves[channel].Y[ip]-FltWaves[channel].Y[ip]/*noise_1[channel].Y[ip]*/+pbatt->Eval(X[ip]);
											//float signal_1=FltWaves[channel].Y[ip]-pbatt->Eval(X[ip]);
											//if ((F_batt<294.3-3*3861||F_batt>294.3+3*3861)&&ampl_b>=(1.7e-4+3*0.6e-4)&&chi<(2e-7+3*0.49e-7)) {signal_1+=-pbatt->Eval(X[ip]);}
											//subtract beat fn

											if (min_value<thr&&(X[ip]>min_range && X[ip]<max_range)&& (chi<9.75e-7+5*4.04e-7|| F_batt<539.1-3*4250||F_batt>539.1+3*4250)){
												signal_1=signal_1-(pbatt->Eval(X[ip]));
											}

											Waves_signal_1[channel].addPnt(signal_1);//, (ip<(dim-skipLstBin)) );
											Waves_noise_tot[channel].addPnt(noise_tot);

										}

										//stampa i risultati del fit

										/*	tmpflt_batt_2.push_back( new TGraph ( dim, &X[0], &flt_batt_2[channel].Y[0]) );
											tmpflt_batt_2.back()->GetYaxis()->SetRangeUser(-0.005,0.01);
											float range_2 = X[dim-1]-X[0];
											pbatt->SetParameters(1,28e+6,28e+6,(range_2)*0.5);
											tmpflt_batt_2.back()->Fit(pbatt);//,"","",X[0]+range*(0.1),X[0]+range*(0.9));
											tmpflt_batt_2.back()->Fit(pbatt);//,"","",X[0]+range*(0.1),X[0]+range*(0.9));
											tmpflt_batt_2.back()->Fit(pbatt);
											TFitResultPtr result_2 = tmpflt_batt.back()->Fit(pbatt,"S");  //stampa i risultati del fit
											float F_batt_2=0.5*(result->Parameter(2)-result->Parameter(1));*/
									}

									NPeak=0;
									NPeak_1=0;
									//smooth sulla forma d'onda di segnale
									int m,k;
									m=23;
									k=2;
									std::vector<float> tmpWSG_signal_23=smoothSG(Waves_signal_1[channel].Y,m,k);
									tmpWSG_signal[channel].fillWave(tmpWSG_signal_23);


//*************************************** RIEMPIMENTO ISTOGRAMMI SUI SEGNALI FILTRATI!!!*************
									float scaleInt=1.0;

									if(!isTrg){

										NPeak = FindPeaks(((wave)Waves_signal_1[channel]).nPtInR(),&((wave)Waves_signal_1[channel]).Y[skipFstBin],1.2e-3/*0.67*((wave)Waves_signal_1[channel]).rms*/,pkPos,pkHgt);
										float sig=2.5*/*1.414**//*1.0e-3;/*/((wave)Waves_signal_1[channel]).rms;
										float sig2=2.0*sig;
										float sig3=1.414*sig;
										//std::cout<<"ev "<<jentry<<" ch "<<channel<<" rms "<<((wave)Waves_signal_1[channel]).rms<<" sig "<<sig<<std::endl;
										//std::cout<<"max "<<((wave)Waves_signal_1[channel]).nMaxInR()<<" pos "<<((wave)Waves_signal_1[channel]).maxInRPos<<std::endl;
										for (int ipk=0; ipk <NPeak; ++ipk) {
											bool skip=false;
											bool doCheck=false;
											int iPkBin=skipFstBin+pkPos[ipk];
											//std::cout<<"Check pk "<<ipk <<" pos "<<iPkBin<<" tpos "<<timeRes*(iPkBin+1)<<" val "<<((wave)Waves_signal_1[channel]).nnAt(skipFstBin+pkPos[ipk])<<std::endl;
											if ( ((wave)Waves_signal_1[channel]).nnAt(skipFstBin+pkPos[ipk]) > sig2 ) {
												if ( ( ((wave)Waves_signal_1[channel]).nMaxInR()>0.03 && iPkBin>((wave)Waves_signal_1[channel]).maxInRPos )
														//									|| ((wave)Waves_signal_1[channel]).nnAt(skipFstBin+pkPos[ipk]) >sig
												) {
													doCheck=true;
												}
											} else {
												if ( ((wave)Waves_signal_1[channel]).nnAt(skipFstBin+pkPos[ipk]) <1.2*sig ) { skip=true; }
												else { doCheck=true; }
											}
											if (doCheck) {
												for (int ick=1; ick <4; ++ick) {
													//std::cout<<"ick "<<ick<<" diff "<<fabs(((wave)Waves_signal_1[channel]).Y[iPkBin]-((wave)Waves_signal_1[channel]).Y[iPkBin-ick])<<" diff1 "<<fabs(((wave)Waves_signal_1[channel]).Y[iPkBin]-((wave)Waves_signal_1[channel]).Y[iPkBin+ick]);
													//std::cout<<" vpk "<<((wave)Waves_signal_1[channel]).Y[iPkBin]<<" vpk-ick "<<((wave)Waves_signal_1[channel]).Y[iPkBin-ick]<<" vpk+ick "<<((wave)Waves_signal_1[channel]).Y[iPkBin+ick]<<std::endl;
													if ( fabs(((wave)Waves_signal_1[channel]).Y[iPkBin]-((wave)Waves_signal_1[channel]).Y[iPkBin-ick])>sig3
															|| fabs(((wave)Waves_signal_1[channel]).Y[iPkBin]-((wave)Waves_signal_1[channel]).Y[iPkBin+ick])>sig3	) {
														skip=true;
														//std::cout<<"Peak reg - iPk "<<ipk <<" tpos "<<timeRes*(iPkBin+1)<<std::endl;
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





										if ( HstPerCh.find(channel)==HstPerCh.end() ){
											TDirectory *chDir = theFile->mkdir(Form("H-Ch_%d_signal",channel));
											chDir->cd();
											HstPerCh.insert(make_pair(channel,new hstPerCh(channel)));
											//	if (!evalWaveCut)  HstSF10PerCh.insert(make_pair(channel,new hstPerCh(channel,SF))); //make_pair funzione per creare una coppia.
											theFile->cd("/");

											((hstPerCh*)HstPerCh[channel])->hNPeaks->GetXaxis()->SetTitle("N Peaks #");
											((hstPerCh*)HstPerCh[channel])->hNPeaks->GetYaxis()->SetTitle("Entries");
											((hstPerCh*)HstPerCh[channel])->hBsl->GetXaxis()->SetTitle("Baseline [V]");
											((hstPerCh*)HstPerCh[channel])->hBsl->GetYaxis()->SetTitle("Entries");
											((hstPerCh*)HstPerCh[channel])->hIntegNInR-> GetXaxis()->SetTitle("Integral [V]");
											((hstPerCh*)HstPerCh[channel])->hIntegNInR-> GetYaxis()->SetTitle("Entries");
											((hstPerCh*)HstPerCh[channel])->hIntegNInRC1-> GetXaxis()->SetTitle("Integral [V]");
											((hstPerCh*)HstPerCh[channel])->hIntegNInRC1-> GetYaxis()->SetTitle("Entries");
											((hstPerCh*)HstPerCh[channel])->hIntegNInRC2-> GetXaxis()->SetTitle("Integral [V]");
											((hstPerCh*)HstPerCh[channel])->hIntegNInRC2-> GetYaxis()->SetTitle("Entries");
											((hstPerCh*)HstPerCh[channel])->hRms->GetXaxis()->SetTitle("Rms [V]");
											((hstPerCh*)HstPerCh[channel])->hRms->GetYaxis()->SetTitle("Entries");

											((hstPerCh*)HstPerCh[channel])->hMaxVNInR->GetXaxis()->SetTitle("Max_value [V]");
											((hstPerCh*)HstPerCh[channel])->hMaxVNInR->GetYaxis()->SetTitle("Entries");
											((hstPerCh*)HstPerCh[channel])->hMinVNInR->GetYaxis()->SetTitle("Entries");
											((hstPerCh*)HstPerCh[channel])->hMinVNInR->GetXaxis()->SetTitle("Min value [V]");

											((hstPerCh*)HstPerCh[channel])->hIntegNInRoriginalW-> GetXaxis()->SetTitle("Integral [mA]");
											((hstPerCh*)HstPerCh[channel])->hIntegNInRoriginalW-> GetYaxis()->SetTitle("Entries");

											((hstPerCh*)HstPerCh[channel])->hMinVNInRoriginalW->GetXaxis()->SetTitle("Min value [V]");
											((hstPerCh*)HstPerCh[channel])->hMinVNInRoriginalW->GetYaxis()->SetTitle("Entries");
											((hstPerCh*)HstPerCh[channel])->hMaxVNInRoriginalW->GetXaxis()->SetTitle("Min_value [V]");
											((hstPerCh*)HstPerCh[channel])->hMaxVNInRoriginalW->GetYaxis()->SetTitle("Entries");


										}

										((hstPerCh*)HstPerCh[channel])->hBsl->Fill(((wave)Waves_signal_1[channel]).bsln);
										((hstPerCh*)HstPerCh[channel])->hInteg->Fill(((wave)Waves_signal_1[channel]).integ);
										((hstPerCh*)HstPerCh[channel])->hIntegN->Fill(((wave)Waves_signal_1[channel]).nnInteg());
										((hstPerCh*)HstPerCh[channel])->hIntegInR->Fill(((wave)Waves_signal_1[channel]).integInR);
										((hstPerCh*)HstPerCh[channel])->hIntegNInR->Fill(((wave)Waves_signal_1[channel]).nnIntegInR());

										if (NPeak_1>0) {
											((hstPerCh*)HstPerCh[channel])->hNPeaks->Fill(NPeak_1);
											for (int ipk=0; ipk <NPeak_1; ++ipk){
												((hstPerCh*)HstPerCh[channel])->hHPeaks->Fill(pkHgt_1[ipk]);
												((hstPerCh*)HstPerCh[channel])->hHNPeaks->Fill(ipk+1,pkHgt_1[ipk]);
											}

											((hstPerCh*)HstPerCh[channel])->hIntegNInRC1->Fill(((wave)Waves_signal_1[channel]).nnIntegInR()/((float)NPeak_1));

											if (((wave)Waves_signal_1[channel]).nnIntegInR()>0.2) {
												((hstPerCh*)HstPerCh[channel])->hNPeaks_1->Fill(NPeak_1);
												((hstPerCh*)HstPerCh[channel])->hTFstPeaks->Fill(X[pkPos_1[0]+skipFstBin]);
												for (int ipk=0; ipk <NPeak_1; ++ipk){
													((hstPerCh*)HstPerCh[channel])->hTPeaks->Fill(X[pkPos_1[ipk]+skipFstBin]);
												}

//												float scaleInt=1.0;
												float minDist=1e+20;
												float tmpDist;
												int clstFreq=0;
												if (numOfPeaks>0) {
													for (int j=0 ; j<numOfPeaks ; j++){
														for (int iFpk=0; iFpk<4; ++iFpk) {
															tmpDist=fabs(peaks[j]-omegaRefs[iFpk]);
															if (tmpDist<minDist) {
																minDist=tmpDist;
																clstFreq=iFpk;
															}
														}
														scaleInt*=integLoss[clstFreq];
													}
												}
//												cout<< " scaleint "<< scaleInt << endl;
												((hstPerCh*)HstPerCh[channel])->hIntegNInRC2->Fill( ( ((wave)Waves_signal_1[channel]).nnIntegInR()/((float)NPeak_1) )/scaleInt );

											}

										}

										((hstPerCh*)HstPerCh[channel])->hRms->Fill(((wave)Waves_signal_1[channel]).rms);

										((hstPerCh*)HstPerCh[channel])->hMaxV->Fill(((wave)Waves_signal_1[channel]).max);
										((hstPerCh*)HstPerCh[channel])->hMaxVN->Fill(((wave)Waves_signal_1[channel]).nMax());
										((hstPerCh*)HstPerCh[channel])->hMaxVNSmooth->Fill(((wave)tmpWSG_signal[channel]).nMax());
										((hstPerCh*)HstPerCh[channel])->hMaxVInR->Fill(((wave)Waves_signal_1[channel]).maxInR);
										((hstPerCh*)HstPerCh[channel])->hMaxVNInR->Fill(((wave)Waves_signal_1[channel]).nMaxInR());

										((hstPerCh*)HstPerCh[channel])->hMinV->Fill(((wave)Waves_signal_1[channel]).min);
										((hstPerCh*)HstPerCh[channel])->hMinVN->Fill(((wave)Waves_signal_1[channel]).nMin());
										((hstPerCh*)HstPerCh[channel])->hMinVInR->Fill(((wave)Waves_signal_1[channel]).minInR);
										((hstPerCh*)HstPerCh[channel])->hMinVNInR->Fill(((wave)Waves_signal_1[channel]).nMinInR());

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
										//if (Waves_signal_1[channel].nnIntegInR()>0 && Waves_signal_1[channel].nnIntegInR()<= 0.02){cout<< " canale "<< channel << " evento del cavolo "<< jentry << " valore integrale " << Waves_signal_1[channel].nnIntegInR() << " baseline " << Waves_signal_1[channel].bsln << endl;}
										//istogrammi rms wave non filtrate

										((hstPerCh*)HstPerCh[channel])->hRmsOriginalW->Fill(((wave)Waves[channel]).rms);

	//_________________________________________SEARCHFOR FULL CHANNEL ___________________________________________//

										float max_ch=0.;

										for(int i=0;i<dim;++i){
											if(max_ch<tmpWSG_signal[channel].Y[i]){
												max_ch=tmpWSG_signal[channel].Y[i];
											}
										}
//__________________________________________search for wave without smooth__________________________________//
//										for(int i=0;i<dim;++i){
//											if(max_ch<Waves_signal_1[channel].Y[i]){
//												max_ch=Waves_signal_1[channel].Y[i];
//											}
//										}
//__________________________________________________________________________________________________________//

//_________________________________ channel by channel control_______________________________________//
										if(channel==23){
//											if(max_ch>=0.0054+3*0.00033){//no sm
											if(max_ch>=0.0026+3*0.0006){
												//												((hstPerCh*)HstPerCh[channel])->hIntegNInRFullW->Fill(( ((wave)Waves_signal_1[channel]).nnIntegInR()/((float)NPeak_1) )/scaleInt);
												Full.first=jentry;
												Full.second=channel;
												waveFull.push_back(Full);
											}
										}

										if(channel==8){
//											if(max_ch>=0.0075+3*0.0014){//no sm
											if(max_ch>=0.004+3*0.0012){
												Full.first=jentry;
												Full.second=channel;
												waveFull.push_back(Full);
											}
										}

										if(channel==5){
//											if(max_ch>=0.0071+3*0.0009){//nosm
											if(max_ch>=0.003+3*0.0007){
											Full.first=jentry;
												Full.second=channel;
												waveFull.push_back(Full);
											}
										}

										if(channel==3){
//											if(max_ch>=0.0075+3*0.0012){//nosm
											if(max_ch>=0.0038+3*0.0011){
												Full.first=jentry;
												Full.second=channel;
												waveFull.push_back(Full);
											}
										}

										if(channel==2){
//											if(max_ch>=0.0068+3*0.0009){//nosm
											if(max_ch>=0.0034+3*0.0009){
												Full.first=jentry;
												Full.second=channel;
												waveFull.push_back(Full);
											}
										}

										if(channel==16){
//											if(max_ch>=0.0066+3*0.0009){//nosm
											if(max_ch>=0.0029+3*0.0007){
												Full.first=jentry;
												Full.second=channel;
												waveFull.push_back(Full);
											}
										}

										if(channel==0){
//											if(max_ch>=0.0072+3*0.0012){//nosm
											if(max_ch>=0.0034+3*0.0009){
												Full.first=jentry;
												Full.second=channel;
												waveFull.push_back(Full);
											}
										}

										if(channel==14){
//											if(max_ch>=0.0066+3*0.0008){//nosm
											if(max_ch>=0.0028+3*0.0005){
												Full.first=jentry;
												Full.second=channel;
												waveFull.push_back(Full);
											}
										}

										if(channel==21){
//											if(max_ch>=0.0065+3*0.0008){//nosm
											if(max_ch>=0.0028+3*0.0005){
												Full.first=jentry;
												Full.second=channel;
												waveFull.push_back(Full);
											}
										}

										if(channel==18){
//											if(max_ch>=0.0065+3*0.0008){//nosm
											if(max_ch>=0.0026+3*0.0005){
												Full.first=jentry;
												Full.second=channel;
												waveFull.push_back(Full);
											}
										}

										if(channel==1){
//											if(max_ch>=0.0070+3*0.0011){//nosm
											if(max_ch>=0.0031+3*0.0007){
												Full.first=jentry;
												Full.second=channel;
												waveFull.push_back(Full);
											}
										}

										if(channel==13){
//											if(max_ch>=0.0070+3*0.0009){//nosm
											if(max_ch>=0.0033+3*0.0007){
												Full.first=jentry;
												Full.second=channel;
												waveFull.push_back(Full);
											}
										}

										if(channel==12){
//											if(max_ch>=0.0070+3*0.0011){//nosm
											if(max_ch>=0.0031+3*0.0007){
												Full.first=jentry;
												Full.second=channel;
												waveFull.push_back(Full);
											}
										}

										if(channel==11){
//											if(max_ch>=0.0079+5*0.0012){//nosm
											if(max_ch>=0.0042+3*0.0010){
												Full.first=jentry;
												Full.second=channel;
												waveFull.push_back(Full);
											}
										}

										if(channel==10){
//											if(max_ch>=0.0076+5*0.0011){//nosm
											if(max_ch>=0.0042+3*0.0010){
												Full.first=jentry;
												Full.second=channel;
												waveFull.push_back(Full);
											}
										}

										if(channel==9){
//											if(max_ch>=0.0077+5*0.0010){//nosm
											if(max_ch>=0.0045+3*0.0011){
												Full.first=jentry;
												Full.second=channel;
												waveFull.push_back(Full);
											}
										}
									}//fine di if is trg


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


									m=23;
									k=2;
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

									/*m=7;
					                  k=4;
					                  std::vector<float> tmpWSG_4_7=smoothSG(Waves[channel].Y,m,k);
					                  SmtSGWaves_4_7[channel].fillWave(tmpWSG_4_7);*/



//**********---------------------------------------disegno dei grafici delle forme d'onda sperimentali con smooth---------------------------------
									 saveWave=true; //true qui
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

									  /*SF=10;
              	  	  	  	  	  	  	tmpWaves.push_back( new TGraph ( dim-2*SF, &X[SF], &Smt10Waves[channel].Y[0]) );
   									    tmpWaves.back()->SetLineColor(kBlue);
   										tmpWaves.back()->Draw("Lsame");*/

									 /*SF=4;
  						               tmpWaves.push_back( new TGraph ( dim-2*SF, &X[SF], &Smt4Waves[channel].Y[0]) );
  						               tmpWaves.back()->SetLineColor(kRed);
                                       tmpWaves.back()->Draw("Lsame");*/

									/*int m=5;
                        			  int k=2;
                                      tmpWaves.push_back( new TGraph ( dim-2*m, &X[m], &SmtSGWaves[channel].Y[0]) );
                                      tmpWaves.back()->SetLineColor(kGreen);
  						              tmpWaves.back()->Draw("Lsame");*/

							    	/*m=7;
					        	      k=2;
					      	          tmpWaves.push_back(new TGraph(dim-2*m,&X[m],&SmtSGWaves_7[channel].Y[0]));
						              tmpWaves.back()->SetLineColor(kOrange);
						              tmpWaves.back()->Draw("Lsame");  */

								    /*m=13;
						              k=2;
						              tmpWaves.push_back(new TGraph(dim-2*m,&X[m],&SmtSGWaves_13[channel].Y[0]));
						              tmpWaves.back()->SetLineColor(kMagenta);
						              tmpWaves.back()->Draw("Lsame"); */

//									  m=23;
//									  k=2;
//									  tmpWaves.push_back(new TGraph(dim-m,&X[(m-1)/2],&SmtSGWaves_23[channel].Y[0]));
//									  tmpWaves.back()->SetLineColor(kGreen);
//									  tmpWaves.back()->Draw("Lsame");

									/*m=9;
						              k=2;
						              tmpWaves.push_back(new TGraph(dim-2*m,&X[m],&SmtSGWaves_9[channel].Y[0]));
						              tmpWaves.back()->SetLineColor(kPink);
						              tmpWaves.back()->Draw("Lsame");   */

									/*m=7;
						              k=4;
						              tmpWaves.push_back(new TGraph(dim-2*m,&X[m],&SmtSGWaves_4_7[channel].Y[0]));
						              tmpWaves.back()->SetLineColor(kAzure-1);
						              tmpWaves.back()->Draw("Lsame");*/

									  tmpCv.back()->Write();
									  theFile->cd("/");
									}
//***--------------------------------------Disegno dei grafici della wave di battimento negativa--------------------------------------
									bool saveNegative=true;
//									bool saveNegative=false;
									if(saveNegative && !isTrg){

										NegWaveDir->cd();

										tmpCvNegBatt.push_back(new TCanvas(Form("CvNegWave-Ch%d_ev%d",channel,jentry),Form("tmpNegWave-Ch%d_ev%d",channel,jentry)));
										tmpNegBatt.back()->GetXaxis()->SetTitle("Time (ns)");
										tmpNegBatt.back()->GetYaxis()->SetTitle("Volt");
										tmpNegBatt.back()->GetYaxis()->SetTitleOffset(0.5);
										tmpNegBatt.back()->SetTitle(Form("tmpWave-Negative-Ch%d_ev%d",channel,jentry));
										tmpNegBatt.back()->Draw("AL");

										tmpCvNegBatt.back()->Write();
										theFile->cd("/");
									}
//****--------------------------------------Disegno dei grafici della trasformata-------------------------------
									bool saveWaveFFT=false;//true qui
									//bool saveWaveFFT=false;
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
										tmpFFTAmpl.push_back( new TGraph ( dim/2, &Wffts[channel].omega[0],amplFltFFT_3) );
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

	//****************-------------------------Disegno della flt_batt per visualizzare il fit del battimento1----------------------------

									bool saveflt_batt=false;//true qui
//									bool saveflt_batt=false;
									if (saveflt_batt && !isTrg) {

										flt_battimento->cd();

										tmpCvflt_batt.push_back( new TCanvas(Form("Cvflt_batt-Ch%d_ev%d",channel,jentry),Form("tmpflt_batt-Ch%d_ev%d",channel,jentry)) );
										tmpCvflt_batt.back()->cd();
										tmpflt_batt.back()->GetXaxis()->SetTitle("time [ns]");
										tmpflt_batt.back()->SetTitle(Form("tmpflt_batt-Ch%d_ev%d",channel,jentry));
										tmpflt_batt.back()->GetYaxis()->SetTitleOffset(1.4);
										tmpflt_batt.back()->GetYaxis()->SetTitle("Volt");
										tmpflt_batt.back()->Draw("AL");
										tmpCvflt_batt.back()->Write();
										theFile->cd("/");
									}

											//---------------------Disegno della flt_batt per visualizzare il fit del battimento2----------------------------
									 bool saveflt_batt2=false;
										if (saveflt_batt2 && !isTrg) {

											flt_battimento2->cd();

											tmpCvflt_batt_2.push_back( new TCanvas(Form("Cvflt_batt_2-Ch%d_ev%d",channel,jentry),Form("tmpflt_batt_2-Ch%d_ev%d",channel,jentry)) );
											tmpCvflt_batt_2.back()->cd();
											tmpflt_batt_2.push_back( new TGraph ( dim, &X[0], &flt_batt2[channel].Y[0]) );
											tmpflt_batt_2.back()->GetXaxis()->SetTitle("time [ns]");
											tmpflt_batt_2.back()->SetTitle(Form("tmpflt_batt_2-Ch%d_ev%d",channel,jentry));
											tmpflt_batt_2.back()->GetYaxis()->SetTitleOffset(1.4);
											tmpflt_batt_2.back()->GetYaxis()->SetTitle("Volt");
											tmpflt_batt_2.back()->Draw("AL");
											tmpCvflt_batt_2.back()->Write();
											theFile->cd("/");
										}


//									bool save_AfterMiddleN=true;
//									if (save_AfterMiddleN && !isTrg) {
//										 After_MidNotch->cd();
//											tmpCvMidNotch.push_back( new TCanvas(Form("Cvflt_MidNotch-Ch%d_ev%d",channel,jentry),Form("tmpflt_MidNotch-Ch%d_ev%d",channel,jentry)) );
//											tmpCvMidNotch.back()->cd();
//											tmpMid_Notch.push_back( new TGraph ( dim, &X[0], &FltMidNotch[channel].Y[0]) );
//											tmpMid_Notch.back()->GetXaxis()->SetTitle("time [ns]");
//											tmpMid_Notch.back()->SetTitle(Form("tmpWave_afterMidNotch-Ch%d_ev%d",channel,jentry));
//											tmpMid_Notch.back()->GetYaxis()->SetTitleOffset(1.4);
//											tmpMid_Notch.back()->GetYaxis()->SetTitle("Volt");
//											tmpMid_Notch.back()->Draw("AL");
//											tmpCvMidNotch.back()->Write();
//											theFile->cd("/");
//
//									}




									//--------------------Disegno della wave di segnale FINALE: signal_1--------------------------------------------------------
									//creazione della mappa per la sovrapposizione dello smooth di SG
//									m=23;
//									k=2;
//									std::vector<float> tmpWSG_23_signal=smoothSG(Waves_signal_1[channel].Y,m,k);
//									SmtSGWaves_23[channel].fillWave(tmpWSG_23_signal);


//									if (!isTrg) {
//										NPeak = FindPeaks(((wave)SmtSGWaves_23[channel]).nPtInR(),&((wave)SmtSGWaves_23[channel]).Y[skipFstBin-(m-1)/2],1.0e-3/*0.67*((wave)Waves_signal_1[channel]).rms*/,pkPos,pkHgt);
//										//						float sig=3.0*1.414*/*1.0e-3;/*/((wave)Waves_signal_1[channel]).rms;
//										//						std::cout<<"ev "<<jentry<<" ch "<<channel<<" rms "<<((wave)Waves_signal_1[channel]).rms<<" sig "<<sig<<std::endl;
//										for (int ipk=0; ipk <NPeak; ++ipk) {
//											bool skip=false;
//											int iPkBin=skipFstBin+pkPos[ipk];
//											std::cout<<"iPk "<<ipk <<" tpos "<<timeRes*(iPkBin+1)<<std::endl;
//											if (((wave)Waves_signal_1[channel]).Y[iPkBin]<0.02) {
//												for (int ick=1; ick <3; ++ick) {
//													std::cout<<"ick "<<ick<<" diff "<<fabs(((wave)Waves_signal_1[channel]).Y[iPkBin]-((wave)Waves_signal_1[channel]).Y[iPkBin-ick])<<" diff1 "<<fabs(((wave)Waves_signal_1[channel]).Y[iPkBin]-((wave)Waves_signal_1[channel]).Y[iPkBin+ick]);
//													std::cout<<" vpk "<<((wave)Waves_signal_1[channel]).Y[iPkBin]<<" vpk-ick "<<((wave)Waves_signal_1[channel]).Y[iPkBin-ick]<<" vpk+ick "<<((wave)Waves_signal_1[channel]).Y[iPkBin+ick]<<std::endl;
//													if ( fabs(((wave)Waves_signal_1[channel]).Y[iPkBin]-((wave)Waves_signal_1[channel]).Y[iPkBin-ick])>sig
//															|| fabs(((wave)Waves_signal_1[channel]).Y[iPkBin]-((wave)Waves_signal_1[channel]).Y[iPkBin+ick])>sig	) { skip=true; break; }
//												}
//											}
//											if (!skip /*((wave)Waves_signal_1[channel]).nnAt(skipFstBin+pkPos[ipk]) > 2.5*((wave)Waves_signal_1[channel]).rms*/ ) {
//												pkPos_1[NPeak_1]=pkPos[ipk];
//												pkHgt_1[NPeak_1]=pkHgt[ipk];
//												++NPeak_1;
//											}
//										}
//
//									}

//									bool savesignal_1=false;
									bool savesignal_1=true;
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

										for (int ipk=0; ipk<NPeak_1; ipk++){
											TMarker *tm = new  TMarker(X[pkPos_1[ipk]+skipFstBin/*+(m-1)/2*/], /*fWaveSign**/pkHgt_1[ipk]/*+PeakBaseShift_1*/, 23);
											tm->SetMarkerSize(1.5);
											tm->SetMarkerColor(2);
											tm->Draw();
//											TMarker *tmb = new  TMarker(timeRes*PeakBasePos_1[ipk], fWaveSign*(PeakHei_1[ipk]-PeakRise_1[ipk])+PeakBaseShift_1, 2);
//											tmb->SetMarkerSize(1.5);
//											tmb->SetMarkerColor(50);
//											tmb->Draw();
										}


										m=23;
										k=2;
										tmpsignal_1.push_back(new TGraph(dim-m,&X[(m-1)/2],&tmpWSG_signal[channel].Y[0]));
										tmpsignal_1.back()->SetLineColor(kBlue);
										tmpsignal_1.back()->Draw("Lsame");
										tmpCvsignal_1.back()->Write();
										theFile->cd("/");
									}



									//------------------Disegno della wave di noise totale: noise_total---------------------------------------------
//									bool savenoise_total=true;
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
//									delete [] amplFltFFT;
//									delete [] amplFltFFT_1;
//									delete [] amplFltFFT_2;
									delete [] amplFltFFT_3;


								}// if point



							} // end of channel loop

						}// auto point

					} // end of loop on entry (second loop)
					//					bin_bin->Write();



					/************************************Isto su rms e max*************************************/
//					float *rms_m = new float[nMaxCh];
//					float *max_m = new float[nMaxCh];
//					float *x_channel=new float[nMaxCh];
					std::vector<float> rms_m;
					std::vector<float> max_m;
					std::vector<float> maxCut_m;
					std::vector<float> x_channel;
					for (int channel=0; channel<=nMaxCh; channel++){
//						x_channel[channel]= channel;
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
					//	TGraph* rms_d= new TGraph(nMa<<xCh,x_channel,rms_m);
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
					if (myfile.is_open()){
						for(int i=0;i<waveFull.size();++i){
						myfile << "full_event "<< waveFull[i].first<<" channel "<< waveFull[i].second<< endl;
						}
					}
					if (myfile.is_open()) { myfile.close(); }
				}





//				void read_data::doChDiff(int evN, std::pair<int,int> &chDiff, std::map<int, bool> &isFull, diffWvCont &diffWaves, WvCont &Waves, WvCont &Smt4Waves, int SF, hDiffCont &hDiff, hDiffCont &hDiffPed, TFile *theFile){
//
//					bool diffEv=false;
//					int frstch=chDiff.first;
//					int scndCh=chDiff.second;
//					if (isFull[frstch] && !isFull[scndCh]) {
//						diffEv=true;
//						diffWaves[chDiff].clear();
//						for (int ipt=0; ipt<(dim-2*SF); ++ipt) {
//							float tmpVal = Waves[frstch].Y[ipt+SF]-Smt4Waves[scndCh].Y[ipt];
//							diffWaves[chDiff].addPnt(tmpVal);
//						}
//					}
//					else if (!isFull[frstch] && isFull[scndCh]) {
//						diffEv=true;
//						diffWaves[chDiff].clear();
//						for (int ipt=0; ipt<(dim-2*SF); ++ipt) {
//							float tmpVal = Waves[scndCh].Y[ipt+SF]-Smt4Waves[frstch].Y[ipt];
//							diffWaves[chDiff].addPnt(tmpVal);
//						}
//					}
//
//					if (diffEv) {
//						hDiff[chDiff].hBsl->Fill(diffWaves[chDiff].bsln);
//						hDiff[chDiff].hInteg->Fill(diffWaves[chDiff].integ);
//						hDiff[chDiff].hIntegN->Fill(diffWaves[chDiff].nInteg());
//						hDiff[chDiff].hRms->Fill(diffWaves[chDiff].rms);
//						hDiff[chDiff].hMaxV->Fill(diffWaves[chDiff].max);
//						hDiff[chDiff].hMaxVN->Fill(diffWaves[chDiff].nMax());
//
//						bool saveWave=false;
//						if (saveWave) {
//
//							//			 waveDiffDir->cd();
//							theFile->cd("WavesDiff");
//
//							tmpCv.push_back( new TCanvas(Form("Cv-Ch%d-Ch%d_ev%d",frstch,scndCh,evN),Form("tmpWave-Ch%d-Ch%d_ev%d",frstch,scndCh,evN)) );
//							tmpCv.back()->cd();
//							tmpWaves.push_back( new TGraph ( dim-2*SF, &X[SF], &diffWaves[chDiff].Y[0]) );
//							tmpWaves.back()->GetXaxis()->SetTitle("time [ns]");
//							tmpWaves.back()->SetTitle(Form("tmpWave-Ch%d-Ch%d_ev%d",frstch,scndCh,evN));
//							tmpWaves.back()->GetYaxis()->SetTitleOffset(1.4);
//							tmpWaves.back()->GetYaxis()->SetTitle("Volt");
//							tmpWaves.back()->Draw("AL");
//							//tmpWaves.push_back( new TGraph ( dim-2*SF, &X[SF], &Smt10Waves[channel].Y[0]) );
//							//tmpWaves.back()->SetLineColor(kBlue);
//							//tmpWaves.back()->Draw("Lsame");
//							tmpCv.back()->Write();
//							theFile->cd("/");
//						}
//
//					}
//
//					if (!isFull[frstch] && !isFull[scndCh]) {
//						diffWaves[chDiff].clear();
//						for (int ipt=0; ipt<(dim-2*SF); ++ipt) {
//							float tmpVal = Waves[frstch].Y[ipt+SF]-Smt4Waves[scndCh].Y[ipt];
//							diffWaves[chDiff].addPnt(tmpVal);
//						}
//						hDiffPed[chDiff].hBsl->Fill(diffWaves[chDiff].bsln);
//						hDiffPed[chDiff].hInteg->Fill(diffWaves[chDiff].integ);
//						hDiffPed[chDiff].hIntegN->Fill(diffWaves[chDiff].nInteg());
//						hDiffPed[chDiff].hRms->Fill(diffWaves[chDiff].rms);
//						hDiffPed[chDiff].hMaxV->Fill(diffWaves[chDiff].max);
//						hDiffPed[chDiff].hMaxVN->Fill(diffWaves[chDiff].nMax());
//					}
//
//				}

				/*----------------sottrazione waveform-waveform filtrate---------------------------------*/
//				void read_data::doChDiff_waves(int evN, int channel, WvCont & Waves_diff, WvCont &Waves, WvCont FltWaves, TFile *theFile){
//					bool diffEv=true;
//					diffEv=true;
//					Waves_diff[channel].clear();
//					for (int ipt=0; ipt<(dim); ++ipt) {
//						float tmpVal = Waves[channel].Y[ipt]-FltWaves[channel].Y[ipt];
//						Waves_diff[channel].addPnt(tmpVal);
//					}
//					bool saveWave=true;
//					if (saveWave) {
//						theFile->cd("WavesDiff");
//						tmpCv.push_back( new TCanvas(Form("Cv-Ch%d-Ch%d_ev%d",channel,channel,evN),Form("tmpWave-Ch%d-Ch%d_ev%d",channel,channel,evN)) );
//						tmpCv.back()->cd();
//						tmpWaves.push_back( new TGraph ( dim, &X[0], &Waves_diff[channel].Y[0]) );
//						tmpWaves.back()->GetXaxis()->SetTitle("time [ns]");
//						tmpWaves.back()->SetTitle(Form("tmpWave-Ch%d_waves-Ch%d_FltWaves_ev%d",channel,channel,evN));
//						tmpWaves.back()->GetYaxis()->SetTitleOffset(1.4);
//						tmpWaves.back()->GetYaxis()->SetTitle("Volt");
//						tmpWaves.back()->Draw("AL");
//						tmpCv.back()->Write();
//						theFile->cd("/");
//					}
//				}
//
//
