//#include "read_data.h"
#include "funcUtils.h"

#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TGraph.h>
//#include<TVector.h>
#include "TF1.h"
#include "TDirectory.h"    
//#include "TVirtualFFT.h"
#include "TMath.h"
#include "TFitResult.h"
#include "TMinuit.h"
#include "TGraphErrors.h"
#include "TSpectrum.h"
#include "TAxis.h"
#include "TVirtualFitter.h"
#include "TRandom.h"

#include "IsignalFromCharges.C"

//#include <map>
#include <fstream>
#include <string>
//#include <vector>
#include <iostream>

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////Define Data Containers//////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
std::vector<float>X;//definisce l'asse x
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
std::vector<TGraph *> tmpsignal_1;       //grafici per il segnale finale dopo i filtri 
std::vector<TCanvas *> tmpCvsignal_1;    //canvas per il segnale finale dopo i filtri 
std::vector<TGraph *> tmpNoise_total;    //grafici per il noise totale 
std::vector<TCanvas *> tmpCvNoise_total; //canvas per il noise totale 

Float_t pulsPerCh[5]={0, 200e+6, 383e+6, 443e+6, 485e+6};

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////MAIN LOOP function//////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void DoSim(Char_t *output, Int_t eventn, Int_t nChannels)
{

	int dim = 1024;
	float _tmax = 1.0e-6;
	float timeRes= (_tmax*1.0e+9)/((float)dim);//in nano sec

	/*********** creazione degli istogrammi***********/

	//	char *basename(char *path);
	//	cout << "Basename " << basename(output) << endl;
	Char_t outChar[500];
	TString out = basename(output);
	//	out.ReplaceAll(".root","");
#ifndef _OSC
	sprintf(outChar,"hSimTB_%s.root",out.Data());   //Data sono i dati convertiti nei file.root
#else
	sprintf(outChar,"hSimOSC_%s.root",out.Data());  //scrive il titolo della rootupla
	dim = 1246;
	timeRes=0.8;
	_tmax = dim*timeRes*1e-9;
	skipFstBin=250; //set Bins to Skip
	skipLstBin=300;
#endif

	tmax=_tmax;  //to fix compilation problem

	std::cout << "The output file is " << outChar << std::endl;
	TFile * theFile = new TFile(outChar,"RECREATE");  //crea un nuovo file.root e se ne esiste gi� uno lo sovrascrive
	/***********creazione delle directory***********/
	theFile->cd("/");
	TDirectory *waveDir = theFile->mkdir("Waves");
	theFile->cd("/");
	TDirectory *waveFFTDir = theFile->mkdir("WavesFFT");
	theFile->cd("/");
	TDirectory *signal = theFile->mkdir("signal_Afterflt"); //directory per la waveform di segnale individuata dopo tutti i filtri
	theFile->cd("/");
	TDirectory *noise = theFile->mkdir("noise_total");      //directory per il noise totale
	theFile->cd("/");
	/*Contenitori di waveform*/
	WvCont Waves;
	WvCont Waves_diff;      
	WvCont FltWaves;
	WvCont FltWavesSM;
	WvCont noise_1;
	WvCont flt_batt;
	WvCont Waves_signal_1; 

	//diffWvCont diffWaves;
	fftCont Wffts; 
	fftCont CR;
	fftCont RC;
	fftCont Notch_1;
	fftCont Notch_2;
	fftCont Notch_3;
	fftCont Notch_4;
	fftCont Notch_5;
	//-------------------------------------------creazione mappa per smooth di SG------------------------------------------------------
	WvCont SmtSGWaves;
	WvCont SmtSGWaves_7;
	WvCont SmtSGWaves_9;
	WvCont SmtSGWaves_13;
	WvCont SmtSGWaves_23;
	WvCont SmtSGWaves_4_9;
	WvCont SmtSGWaves_4_7;
	WvCont smtSG;
	std::map<int, hstPerCh *> HstPerCh;

	//std::vector<float> Y; //vettore da riempire con il segnale simulato.

	Float_t ChargesSequence[1];
	Float_t timeSequence[1];
	ChargesSequence[0]=5e+5;
	Int_t outCharge;

	X.clear();
	for (int i = 0; i < dim; ++i) { X.push_back(timeRes *(i+1)); }
//	int channel = 1;

	for (int jentry=0; jentry<eventn; ++jentry) {

		Waves.clear();
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
		flt_batt.clear();
		//flt_batt_2.clear();
		SmtSGWaves.clear();
		SmtSGWaves_7.clear();
		SmtSGWaves_9.clear();
		SmtSGWaves_13.clear();
		SmtSGWaves_4_9.clear();
		SmtSGWaves_4_7.clear();
		SmtSGWaves_23.clear();
		Waves_signal_1.clear();


		timeSequence[0]=( skipFstBin+50+(dim-(skipFstBin+skipLstBin+100))*gRandom->Uniform() )*timeRes*1e-3; //in us
		Float_t *data;
		outCharge=0;
		data=IsignalFromCharges( ChargesSequence,timeSequence, 1, outCharge, _tmax, 0.0010, 0.5, 1500.0);//, 10.4e-6);

		for (int channel=0; channel<nChannels; ++channel) {

			for(int i =0; i<dim; ++i){
				Waves[channel].addPnt(data[i], (i<(dim-skipLstBin)));
	//			std::cout << Waves[channel].Y[i] << " wave " << std::endl;
			}

			/****----------------------------TRASFORMATA DI FOURIER---------------------------------*****/
			FFT(Waves[channel],Wffts[channel]);//trasformata di Fourier
			double *realFltFFT = new double[Waves[channel].nPt()];
			double *imgFltFFT = new double[Waves[channel].nPt()];
			/****---------------------------------------------------------------------------------------*****/
			//****Filtraggio noise, tramite Notch*****
			if (channel==0) {
//				filterNotch(Wffts[channel],realFltFFT,imgFltFFT,(TMath::TwoPi()*538e+6),TMath::TwoPi()*20e+6,2);
				filterNotch(Wffts[channel],realFltFFT,imgFltFFT,(TMath::TwoPi()*200e+6),TMath::TwoPi()*200e+6,0);
			}
			else {
				filterNotch(Wffts[channel],realFltFFT,imgFltFFT,pulsPerCh[channel]/*(TMath::TwoPi()*46e+6)*/,TMath::TwoPi()*200e+6,0);
				Notch_2[channel].fillRM(Wffts[channel].nPt,realFltFFT,imgFltFFT,Wffts[channel].omega);
				//filterNotch(Notch_2[channel],realFltFFT,imgFltFFT,(TMath::TwoPi()*538e+6),TMath::TwoPi()*20e+6,2);
			}
			InverseFFT(realFltFFT,imgFltFFT,Waves[channel].nPt(),FltWaves[channel]);//trasformata inversa
			double *amplFltFFT_2 = new double[Waves[channel].nPt()];
			for (int ip=0; ip<Waves[channel].nPt()/2; ++ip) {
				amplFltFFT_2[ip]=TMath::Sqrt(realFltFFT[ip]*realFltFFT[ip]+imgFltFFT[ip]*imgFltFFT[ip]);
			}
			for (int ip=Waves[channel].nPt()/2; ip<Waves[channel].nPt(); ++ip) {
				amplFltFFT_2[ip]=0.0;
			}

			for (int ip=0; ip<dim; ++ip) {
				float signal_1=FltWaves[channel].Y[ip];
				float noise_tot=Waves[channel].Y[ip]-FltWaves[channel].Y[ip];
				Waves_signal_1[channel].addPnt(signal_1, (ip<(dim-skipLstBin)) ); //SEGNALE FILTRATO!!!
			}



			//CREAZIONE DIRECTORY E RIEMPIMENTO ISTOGRAMMI, SUI SEGNALI FILTRATI!!!
			if ( HstPerCh.find(channel)==HstPerCh.end() ){
				TDirectory *chDir = theFile->mkdir(Form("H-Ch_%d_signal",channel));
				chDir->cd();
				HstPerCh.insert(std::make_pair(channel,new hstPerCh(channel)));
				//	if (!evalWaveCut)  HstSF10PerCh.insert(make_pair(channel,new hstPerCh(channel,SF))); //make_pair funzione per creare una coppia.
				theFile->cd("/");

				((hstPerCh*)HstPerCh[channel])->hIntegInR->SetBins(1000,0,1);
				((hstPerCh*)HstPerCh[channel])->hIntegNInR->SetBins(1000,0,1);
				((hstPerCh*)HstPerCh[channel])->hIntegInRoriginalW->SetBins(1000,0,1);
				((hstPerCh*)HstPerCh[channel])->hIntegNInRoriginalW->SetBins(1000,0,1);

				((hstPerCh*)HstPerCh[channel])->hBsl->GetXaxis()->SetTitle("Baseline [mA]");
				((hstPerCh*)HstPerCh[channel])->hBsl->GetYaxis()->SetTitle("Entries");
				((hstPerCh*)HstPerCh[channel])->hIntegNInR-> GetXaxis()->SetTitle("Integral [mA]");
				((hstPerCh*)HstPerCh[channel])->hIntegNInR-> GetYaxis()->SetTitle("Entries");
				((hstPerCh*)HstPerCh[channel])->hRms->GetXaxis()->SetTitle("Rms [mA]");
				((hstPerCh*)HstPerCh[channel])->hRms->GetYaxis()->SetTitle("Entries");
				((hstPerCh*)HstPerCh[channel])->hMaxVNInR->GetXaxis()->SetTitle("Max_value [mA]");
				((hstPerCh*)HstPerCh[channel])->hMaxVNInR->GetYaxis()->SetTitle("Entries");
				((hstPerCh*)HstPerCh[channel])->hMinVNInR->GetYaxis()->SetTitle("Entries");
				((hstPerCh*)HstPerCh[channel])->hMinVNInR->GetXaxis()->SetTitle("Min value [mA]");
				((hstPerCh*)HstPerCh[channel])->hIntegInRoriginalW-> GetXaxis()->SetTitle("Integral [mA]");
				((hstPerCh*)HstPerCh[channel])->hIntegInRoriginalW-> GetYaxis()->SetTitle("Entries");
				((hstPerCh*)HstPerCh[channel])->hIntegNInRoriginalW-> GetXaxis()->SetTitle("Integral [mA]");
				((hstPerCh*)HstPerCh[channel])->hIntegNInRoriginalW-> GetYaxis()->SetTitle("Entries");
				((hstPerCh*)HstPerCh[channel])->hMinVNInRoriginalW->GetXaxis()->SetTitle("Min value [mA]");
				((hstPerCh*)HstPerCh[channel])->hMinVNInRoriginalW->GetYaxis()->SetTitle("Entries");
				((hstPerCh*)HstPerCh[channel])->hMaxVNInRoriginalW->GetXaxis()->SetTitle("Min_value [mA]");
				((hstPerCh*)HstPerCh[channel])->hMaxVNInRoriginalW->GetYaxis()->SetTitle("Entries");
			}


			((hstPerCh*)HstPerCh[channel])->hBsl->Fill(((wave)Waves_signal_1[channel]).bsln);
			((hstPerCh*)HstPerCh[channel])->hInteg->Fill(((wave)Waves_signal_1[channel]).integ);
			((hstPerCh*)HstPerCh[channel])->hIntegN->Fill(((wave)Waves_signal_1[channel]).nInteg());
			((hstPerCh*)HstPerCh[channel])->hIntegInR->Fill(((wave)Waves_signal_1[channel]).integInR);
			((hstPerCh*)HstPerCh[channel])->hIntegNInR->Fill(((wave)Waves_signal_1[channel]).nIntegInR());
			((hstPerCh*)HstPerCh[channel])->hRms->Fill(((wave)Waves_signal_1[channel]).rms);

			((hstPerCh*)HstPerCh[channel])->hMaxV->Fill(((wave)Waves_signal_1[channel]).max);
			((hstPerCh*)HstPerCh[channel])->hMaxVN->Fill(((wave)Waves_signal_1[channel]).nMax());
			((hstPerCh*)HstPerCh[channel])->hMaxVInR->Fill(((wave)Waves_signal_1[channel]).maxInR);
			((hstPerCh*)HstPerCh[channel])->hMaxVNInR->Fill(((wave)Waves_signal_1[channel]).nMaxInR());

			((hstPerCh*)HstPerCh[channel])->hMinV->Fill(((wave)Waves_signal_1[channel]).min);
			((hstPerCh*)HstPerCh[channel])->hMinVN->Fill(((wave)Waves_signal_1[channel]).nMin());
			((hstPerCh*)HstPerCh[channel])->hMinVInR->Fill(((wave)Waves_signal_1[channel]).minInR);
			((hstPerCh*)HstPerCh[channel])->hMinVNInR->Fill(((wave)Waves_signal_1[channel]).nMinInR());

			((hstPerCh*)HstPerCh[channel])->hIntegInRoriginalW->Fill(((wave)Waves[channel]).integInR);
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

			//		((hstPerCh*)HstPerCh[channel])->hRmsOriginal->Fill(((wave)Waves[channel]).rms);
			//		((hstPerCh*)HstPerCh[channel])->hRmsOriginal->GetXaxis()->SetTitle("Rms [mA]");
			//		((hstPerCh*)HstPerCh[channel])->hRmsOriginal->GetYaxis()->SetTitle("Entries");

			//---------------------------------------disegno dei grafici delle forme d'onda sperimentali--------------------------------
			bool saveWave=true;
//			bool saveWave=false;
			if (saveWave) {
				waveDir->cd();
				tmpCv.push_back( new TCanvas(Form("Cv-Ch%d_ev%d",channel,jentry),Form("tmpWave-Ch%d_ev%d",channel,jentry)) );
				tmpCv.back()->cd();
				tmpWaves.push_back( new TGraph ( dim, &X[0], &Waves[channel].Y[0]) );
				tmpWaves.back()->GetXaxis()->SetTitle("time [ns]");
				tmpWaves.back()->SetTitle(Form("tmpWave-Ch%d_ev%d",channel,jentry));
				tmpWaves.back()->GetYaxis()->SetTitleOffset(1.4);
				tmpWaves.back()->GetYaxis()->SetTitle("mA");
				tmpWaves.back()->Draw("AL");
				tmpCv.back()->Write();
				theFile->cd("/");
			}
			//--------------------------------------Disegno dei grafici della trasformata-------------------------------
//			bool saveWaveFFT=true;
			bool saveWaveFFT=false;
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
			//creazione della mappa
			bool savesignal_1=true;
//			bool savesignal_1=false;
			if (savesignal_1 ) {

				signal->cd();
				tmpCvsignal_1.push_back( new TCanvas(Form("CvSignal_1-Ch%d_ev%d",channel,jentry),Form("tmpSignal_1-Ch%d_ev%d",channel,jentry)) );
				tmpCvsignal_1.back()->cd();
				tmpsignal_1.push_back( new TGraph ( dim, &X[0], &Waves_signal_1[channel].Y[0]) );
				tmpsignal_1.back()->GetXaxis()->SetTitle("time [ns]");
				tmpsignal_1.back()->SetTitle(Form("tmpSignal_afterFlt-Ch%d_ev%d",channel,jentry));
				tmpsignal_1.back()->GetYaxis()->SetTitleOffset(1.4);
				tmpsignal_1.back()->GetYaxis()->SetTitle("mA");
				tmpsignal_1.back()->Draw("AL");
				tmpCvsignal_1.back()->Write();
			}

			//------------------------Pulizia delle variabili--------------------

			delete [] realFltFFT;
			delete [] imgFltFFT;
			delete [] amplFltFFT_2;
			//delete [] amplFltFFT_3;
		}  // end loop on channels
	} // end loop on events

	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//                                               //  Write Histograms //                                                    //
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	theFile->cd();
	theFile->Write();
	theFile->Close();

}







