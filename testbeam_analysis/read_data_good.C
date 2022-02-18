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
	trigCh.push_back(0);
	trigCh.push_back(1);
	trigCh.push_back(2);
	trigCh.push_back(3);
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
	cout << "Number of entries in the tree for real data is= " << nentries << endl;
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
				cout<< "dim "<<dim<<" time "<<timeRes<<endl;
				for (int i = 0; i < dim; ++i) { X.push_back(timeRes *(i+1));  }
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
						std::cout<<"OSCILLOSCOPE dim "<<dim<<" tmax "<<_tmax<<std::endl;
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
			WvCont FltWaves2;
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
//				cout<< nMaxCh<<endl;
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
				FltWaves2.clear();
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
//					cout<<nMaxCh<<endl;
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

/****-------------------------------TRASFORMATA DI FOURIER---------------------------------*****/
									FFT(Waves[channel],Wffts[channel]);//trasformata di Fourier
									double *realFltFFT = new double[Waves[channel].nPt()];
									double *imgFltFFT = new double[Waves[channel].nPt()];

									filterWaveBsl(Wffts[channel],realFltFFT,imgFltFFT);//filtro sulla baseline.                                        ///
									InverseFFT(realFltFFT,imgFltFFT,Waves[channel].nPt(),FltWaves[channel]);//trasformata inversa
									Waves_signal_1[channel].fillWave(Waves/*FltWaves*/[channel].Y);

									bool saveWave=false;
									if (saveWave) {
										waveFltDir->cd();
										tmpCvFlt.push_back( new TCanvas(Form("CvFlt-Ch%d_ev%d",channel,jentry),Form("tmpFltWave-Ch%d_ev%d",channel,jentry)) );
										tmpCvFlt.back()->cd();
										tmpFltWaves.push_back( new TGraph ( dim, &X[0], &FltWaves[channel].Y[0]) );
										tmpFltWaves.back()->GetXaxis()->SetTitle("time [ns]");
										tmpFltWaves.back()->SetTitle(Form("tmpFltWave-Ch%d_ev%d",channel,jentry));
										tmpFltWaves.back()->GetYaxis()->SetTitleOffset(1.4);
										tmpFltWaves.back()->GetYaxis()->SetTitle("Volt");
										tmpFltWaves.back()->GetYaxis()->SetRangeUser(-0.1,0.4);
										tmpFltWaves.back()->Draw("AL");

										tmpCvFlt.back()->Write();
										theFile->cd("/");

									}



//*************************************Riempimento dei vettori per la waveform di battimento negativa+ricerca del range********************


									NPeak=0;
									NPeak_1=0;
									//smooth sulla forma d'onda di segnale
									int m,k;
									m=23; //number of bin interested by the SG smoothing 
									k=2; //order of the polinomial used
									std::vector<float> tmpWSG_signal_23=smoothSG(Waves_signal_1[channel].Y,m,k);
									tmpWSG_signal[channel].fillWave(tmpWSG_signal_23);


//*************************************** RIEMPIMENTO ISTOGRAMMI SUI SEGNALI FILTRATI!!!*************
									float scaleInt=1.0;

									if(!isTrg&&Waves_signal_1[channel].max>0.005/*0.01*/){
										NPeak = FindPeaks(((wave)Waves_signal_1[channel]).nPtInR(),&((wave)Waves_signal_1[channel]).Y[skipFstBin],1.2e-3/*0.625*((wave)Waves_signal_1[channel]).rms*/,6,3,pkPos,pkHgt);
//										if(jentry==1&&channel==9)cout<<"rms "<<((wave)Waves_signal_1[channel]).rms<<endl;
										//0.625*rms= 2 sigma
										float sig=15e-3;//*2.5*((wave)Waves_signal_1[channel]).rms;//3.0*/*1.414**//*1.0e-3;/*/((wave)Waves_signal_1[channel]).rms;
//										float sig=((wave)Waves_signal_1[channel]).rms/2.5;
										float sig2=1.*sig;//1.0*sig;
										float sig3=1.414*sig;;
										//std::cout<<"ev "<<jentry<<" ch "<<channel<<" rms "<<((wave)Waves_signal_1[channel]).rms<<" sig "<<sig<<std::endl;
										//std::cout<<"max "<<((wave)Waves_signal_1[channel]).nMaxInR()<<" pos "<<((wave)Waves_signal_1[channel]).maxInRPos<<std::endl;
										float meanValueLocal;
										for (int ipk=0; ipk <NPeak; ++ipk) {
											bool skip=false;
											bool doCheck=false;
//											meanValueLocal=(((wave)Waves_signal_1[channel]).Y[skipFstBin+pkPos[ipk]-2]+((wave)Waves_signal_1[channel]).Y[skipFstBin+pkPos[ipk]+2])/2;
//											if ( (((wave)Waves_signal_1[channel]).Y[skipFstBin+pkPos[ipk]]-meanValueLocal) <1.2*sig ) { skip=true; }
											int iPkBin=skipFstBin+pkPos[ipk];
//											std::cout<<"Check pk "<<ipk <<" pos "<<iPkBin<<" tpos "<<timeRes*(iPkBin+1)<<" val "<<((wave)Waves_signal_1[channel]).nnAt(skipFstBin+pkPos[ipk])<<std::endl;
											if ( ((wave)Waves_signal_1[channel]).nnAt(skipFstBin+pkPos[ipk]) > sig2 ) {
//												if(channel==11&&jentry==9)cout<<"test 1 "<<jentry<<" channel "<<channel<<endl;
												if ( ( ((wave)Waves_signal_1[channel]).nMaxInR()>0.01 && iPkBin>((wave)Waves_signal_1[channel]).maxInRPos )
														//									|| ((wave)Waves_signal_1[channel]).nnAt(skipFstBin+pkPos[ipk]) >sig
												) {
//													if(channel==11&&jentry==9)cout<<"test 2 "<<jentry<<" channel "<<channel<<endl;

													doCheck=true;
												}
											} else {
//												if(channel==11&&jentry==9)cout<<"test 3 "<<jentry<<" channel "<<channel<<endl;

												if ( ((wave)Waves_signal_1[channel]).nnAt(skipFstBin+pkPos[ipk]) <3*sig ) {
//													if(channel==11&&jentry==9)cout<<"test 4 "<<jentry<<" channel "<<channel<<endl;
													skip=true;
												}else {
													doCheck=true;

//													if(channel==11&&jentry==9)cout<<"test 5 "<<jentry<<" channel "<<channel<<endl;
												}
											}
											if (doCheck) {
//												if(channel==11&&jentry==9)cout<<"test 6 "<<jentry<<" channel "<<channel<<endl;
												for (int ick=1; ick <2; ++ick) {
//													std::cout<<"ick "<<ick<<" diff "<<fabs(((wave)Waves_signal_1[channel]).Y[iPkBin]-((wave)Waves_signal_1[channel]).Y[iPkBin-ick])<<" diff1 "<<fabs(((wave)Waves_signal_1[channel]).Y[iPkBin]-((wave)Waves_signal_1[channel]).Y[iPkBin+ick]);
//													std::cout<<" vpk "<<((wave)Waves_signal_1[channel]).Y[iPkBin]<<" vpk-ick "<<((wave)Waves_signal_1[channel]).Y[iPkBin-ick]<<" vpk+ick "<<((wave)Waves_signal_1[channel]).Y[iPkBin+ick]<<std::endl;
													if ( fabs(((wave)Waves_signal_1[channel]).Y[iPkBin]-((wave)Waves_signal_1[channel]).Y[iPkBin-ick])>sig3
															|| fabs(((wave)Waves_signal_1[channel]).Y[iPkBin]-((wave)Waves_signal_1[channel]).Y[iPkBin+ick])>sig3	) {
														skip=true;
//														std::cout<<"Peak reg - iPk "<<ipk <<" tpos "<<timeRes*(iPkBin+1)<<std::endl;
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

										((hstPerCh*)HstPerCh[channel])->hBsl->Fill(((wave)Waves_signal_1[channel]).bsln);
										((hstPerCh*)HstPerCh[channel])->hInteg->Fill(((wave)Waves_signal_1[channel]).integ);
										((hstPerCh*)HstPerCh[channel])->hIntegN->Fill(((wave)Waves_signal_1[channel]).nnInteg());
										((hstPerCh*)HstPerCh[channel])->hIntegInR->Fill(((wave)Waves_signal_1[channel]).integInR);
										((hstPerCh*)HstPerCh[channel])->hIntegNInR->Fill(((wave)Waves_signal_1[channel]).nnIntegInR());

										if (NPeak_1>2/*0*/) {
//											if(NPeak_1==1||NPeak_1==2) {
//												cout<<"Ev. "<<jentry<<" channel "<<channel<<" N.of Peaks "<<NPeak_1<<endl;
//											}
											((hstPerCh*)HstPerCh[channel])->hNPeaks->Fill(NPeak_1);
											for (int ipk=0; ipk <NPeak_1; ++ipk){
												((hstPerCh*)HstPerCh[channel])->hHPeaks->Fill(pkHgt_1[ipk]);
												((hstPerCh*)HstPerCh[channel])->hHNPeaks->Fill(ipk+1,pkHgt_1[ipk]);
											}

											((hstPerCh*)HstPerCh[channel])->hIntegNInRC1->Fill(((wave)Waves_signal_1[channel]).nnIntegInR()/((float)NPeak_1));

											if (((wave)Waves_signal_1[channel]).nnIntegInR()>0.1/*0.2*/) {
												((hstPerCh*)HstPerCh[channel])->hNPeaks_1->Fill(NPeak_1);
												((hstPerCh*)HstPerCh[channel])->hTFstPeaks->Fill(X[pkPos_1[0]+skipFstBin]);
												for (int ipk=0; ipk <NPeak_1; ++ipk){
													((hstPerCh*)HstPerCh[channel])->hTPeaks->Fill(X[pkPos_1[ipk]+skipFstBin]);
												}

												//												float scaleInt=1.0;
												float minDist=1e+20;
												float tmpDist;
												int clstFreq=0;
												int numOfPeaks=0;

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

									}//fine di if is trg


									std::vector<float> tmpWSF10=smooth(Waves[channel].Y,SF);
									Smt10Waves[channel].fillWave( tmpWSF10 );

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
									 saveWave=false;//true; //true qui
									if (saveWave) {
										waveDir->cd();
										tmpCv.push_back( new TCanvas(Form("Cv-Ch%d_ev%d",channel,jentry),Form("tmpWave-Ch%d_ev%d",channel,jentry)) );
										tmpCv.back()->cd();
										tmpWaves.push_back( new TGraph ( dim, &X[0], &Waves[channel].Y[0]) );
										tmpWaves.back()->GetXaxis()->SetTitle("time [ns]");
										tmpWaves.back()->SetTitle(Form("tmpWave-Ch%d_ev%d",channel,jentry));
										tmpWaves.back()->GetYaxis()->SetTitleOffset(1.4);
										tmpWaves.back()->GetYaxis()->SetTitle("Volt");
										tmpWaves.back()->GetYaxis()->SetRangeUser(-0.1,0.4);
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

//****--------------------------------------Disegno dei grafici della trasformata-------------------------------
									bool saveWaveFFT=false;//true qui
									//bool saveWaveFFT=false;
//									if (saveWaveFFT) {
//
//										waveFFTDir->cd();
//
//										tmpCvFFT.push_back( new TCanvas(Form("CvFFT-Ch%d_ev%d",channel,jentry),Form("tmpWaveFFT-Ch%d_ev%d",channel,jentry)) );
//										tmpCvFFT.back()->Divide(1,3);
//										tmpCvFFT.back()->SetLogx();
//										tmpCvFFT.back()->cd(1)->SetLogy();//primo riquadro : ampiezza
//										tmpFFTAmpl.push_back( new TGraph ( dim/2, &Wffts[channel].omega[0], &Wffts[channel].ampl[0]) );
//										tmpFFTAmpl.back()->GetXaxis()->SetTitle("#omega [Hz]");
//										tmpFFTAmpl.back()->SetTitle(Form("tmpWaveFFT-Ampl-Ch%d_ev%d",channel,jentry));
//										tmpFFTAmpl.back()->GetYaxis()->SetTitleOffset(1.4);
//										tmpFFTAmpl.back()->GetYaxis()->SetTitle("Amplitude");
//										tmpFFTAmpl.back()->Draw("AL");
//										tmpCvFFT.back()->cd(2);//secondo riquadro :fase
//										tmpFFTphi.push_back( new TGraph ( dim/2, &Wffts[channel].omega[0], &Wffts[channel].phi[0]) );
//										tmpFFTphi.back()->GetXaxis()->SetTitle("#omega [Hz]");
//										tmpFFTphi.back()->SetTitle(Form("tmpWaveFFT-Phi-Ch%d_ev%d",channel,jentry));
//										tmpFFTphi.back()->GetYaxis()->SetTitleOffset(1.4);
//										tmpFFTphi.back()->GetYaxis()->SetTitle("#Phi [rad]");
//										tmpFFTphi.back()->Draw("AL");
//										tmpCvFFT.back()->cd(3)->SetLogy();//terzo riquadro: taglio.
//										tmpFFTAmpl.push_back( new TGraph ( dim/2, &Wffts[channel].omega[0],amplFltFFT_3) );
//										tmpFFTAmpl.back()->GetXaxis()->SetTitle("#omega [Hz]");
//										tmpFFTAmpl.back()->SetTitle(Form("tmpWaveFlt-Ampl-Ch%d_ev%d",channel,jentry));
//										tmpFFTAmpl.back()->GetYaxis()->SetTitleOffset(1.4);
//										tmpFFTAmpl.back()->GetYaxis()->SetTitle("Amplitude");
//										tmpFFTAmpl.back()->Draw("AL");
//
//										tmpCvFFT.back()->Write();
//										theFile->cd("/");
//
//									}

//									bool savesignal_1=false;
									bool savesignal_1=false;//true;
									if (savesignal_1 && !isTrg) {

										signal->cd();

										tmpCvsignal_1.push_back( new TCanvas(Form("CvSignal_1-Ch%d_ev%d",channel,jentry),Form("tmpSignal_1-Ch%d_ev%d",channel,jentry)) );
										tmpCvsignal_1.back()->cd();
										tmpsignal_1.push_back( new TGraph ( dim, &X[0], &Waves_signal_1[channel].Y[0]) );
										tmpsignal_1.back()->GetXaxis()->SetTitle("time [ns]");
										tmpsignal_1.back()->SetTitle(Form("tmpSignal_afterFlt-Ch%d_ev%d",channel,jentry));
										tmpsignal_1.back()->GetYaxis()->SetTitleOffset(1.4);
										tmpsignal_1.back()->GetYaxis()->SetTitle("Volt");
										tmpsignal_1.back()->GetYaxis()->SetRangeUser(-0.1,0.4);
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

//
//										m=23;
//										k=2;
//										tmpsignal_1.push_back(new TGraph(dim-m,&X[(m-1)/2],&tmpWSG_signal[channel].Y[0]));
//										tmpsignal_1.back()->SetLineColor(kBlue);
//										tmpsignal_1.back()->Draw("Lsame");
										tmpCvsignal_1.back()->Write();
										theFile->cd("/");
									}



									//------------------------Pulizia delle variabili--------------------

									delete [] realFltFFT;
									delete [] imgFltFFT;
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
//								std::cout<<channel<<"\t"<<nSelEv[channel]<<std::endl;
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
//					if (myfile.is_open()){
//						for(int i=0;i<waveFull.size();++i){
//						myfile << "full_event "<< waveFull[i].first<<" channel "<< waveFull[i].second<< endl;
//						}
//					}
//					if (myfile.is_open()) { myfile.close(); }
				}


