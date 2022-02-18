f#define read_data_cxx
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
#include "TTree.h"

#include "HFitInterface.h"
#include "Fit/Fitter.h"
#include "Math/MinimizerOptions.h"
//////////////Define Data Containers//////////////

vector<float>X;//definisce l'asse x
// Temporary container for tmp Graphs and Canvas
std::vector<TGraph *> tmpWaves;          //grafico delle forme d'onda
std::vector<TCanvas *> tmpCv;            //canvas forme d'onda
std::vector<TGraph *> tmpFFTAmpl;        //grafico per l'ampiezza
std::vector<TGraph *> tmpFFTFltAmpl;     //grafico per ampiezza filtrata
std::vector<TGraph *> tmpFFTphi;        //grafico per la fase
std::vector<TCanvas *> tmpCvFFT;        //canvas per la fft
std::vector<TGraph *> tmpFltWaves;     	//grafici per funzioni d'onda filtrate
std::vector<TCanvas *> tmpCvFlt;       //canvas forme d'onda filtrate(primo filtro:retta lineare di taglio)
std::vector<TH1F *>  tmpFFTAmplPeak;   //grafico per lo spettro in frequenza dell'ampiezza con filtro di smooth
//std::vector<TGraph *> tmpnoise_1;      //grafici per il noise individuato dal filtro RC.(singolo filtro)
//std::vector<TCanvas *> tmpCvnoise_1;   //canvas per il noise individuato dal filtro RC(singolo filtro).
std::vector<TCanvas *> tmpCvPeak;
std::vector<TGraph *> tmpsignal_1;       //grafici per il segnale finale dopo i filtri
std::vector<TCanvas *> tmpCvsignal_1;    //canvas per il segnale finale dopo i filtri
std::vector<TGraph *> tmpNoise_total;    //grafici per il noise totale
std::vector<TCanvas *> tmpCvNoise_total; //canvas per il noise totale
std::vector<TCanvas *> tmpCvsignal2;
//std::vector<TH1F *>  hDrvDistribution;
//void doChDiff(int evN, std::pair<int,int> &chDiff, std::map<int, bool> &isFull, diffWvCont &diffWaves, WvCont &Waves, WvCont &Smt4Waves, int SF, hDiffCont &hDiff, hDiffCont &hDiffPed, TFile *theFile);
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////MAIN LOOP function//////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void read_data::Loop(Char_t *output, Int_t MidEv,Int_t eventn, Bool_t evalWaveCut)
{

	tmax=_tmax;  //to fix compilation problem

	/*********** creazione degli istogrammi***********/

	char *basename(char *path);
	//cout << "Basename " << basename(output) << endl;
	Char_t outChar[500];
	TString out = basename(output);
	out.ReplaceAll(".root","");
#ifndef _OSC
	sprintf(outChar,"histosTB_%s.root",out.Data());   //Data sono i dati convertiti nei file.root
#else
	sprintf(outChar,"PlotFT_%s.root",out.Data());  //scrive il titolo della rootupla
#endif
	//	cout << "The output file is " << outChar << endl;
	TFile * theFile = new TFile(outChar,"RECREATE");  //crea un nuovo file.root e se ne esiste gi� uno lo sovrascrive

	std::map<int, hstPerCh *> HstPerCh;               //mappa per istogrammi sui segnali filtrati.
	//std::map<int, hstPerCh *> HstSF10PerCh;
	int runId=atoi(out(4,5).Data());
	cout<< " run ID " << runId<< endl;
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
	TDirectory *waveFltDir = theFile->mkdir("WavesFlt");
	theFile->cd("/");
	TDirectory *waveFltPeakDir = theFile->mkdir("Amplfft");      //directory per lo spettro in frequenza dell'ampiezza con smooth
	theFile->cd("/");
	//TDirectory *noise_RC = theFile->mkdir("noise_RC");           //directory per la waveform di noise individuata dal filtro RC(noise_1,singolo filtro)
	//theFile->cd("/");
	TDirectory *signal = theFile->mkdir("signal_Afterflt"); //directory per la waveform di segnale individuata dopo tutti i filtri
	theFile->cd("/");
	TDirectory *noise = theFile->mkdir("noise_total");      //directory per il noise totale
	theFile->cd("/");
	TDirectory *WavesCh2 = theFile->mkdir("WavesCh2");      //directory per il noise totale
	theFile->cd("/");

	//fit sinusoidale
	TF1 *sinFit=new TF1("sinFit","[0]*sin(TMath::TwoPi()*([1]*[3]*x-[2]))+[4]",0,10);
	sinFit->SetParLimits(0,0,10);
	sinFit->SetParLimits(1,0,1500);
	sinFit->SetParLimits(2,0,1);
	sinFit->FixParameter(3,1e+6);
	sinFit->SetParLimits(4,-1,1);
	sinFit->SetNpx(1000);

	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	//                                               //  Loop on entries //                                                  //

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	std::vector<int> trigCh; //vettore di interi per i canali di trigger.
	trigCh.clear();
#ifndef _OSC
	trigCh.push_back(32);
	trigCh.push_back(33);
#else
	trigCh.push_back(3);//canale in ingresso
	//	trigCh.push_back(8);//osc2
#endif
	//mappa per i tagli: definire cosa � segnale e cosa � rumore
	std::map<int, std::pair<float, float> > waveCut;
	waveCut.clear();
	int nSig=3;// numero di sigma.(era 3)
	std::map<int, int> nSelEv;
	nSelEv.clear();

	if (fChain == 0) return;                      // puntatore TTree

	Long64_t nentries = fChain->GetEntriesFast(); //Return the number of entries as of the last check.
	//cout << "Number of entries in the tree for real data is= " << nentries << endl;
	//   Long64_t nbytes = 0, nb = 0;

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
	cout<< " MidEv "<< MidEv<< "lastEv "<<lastEv<<endl;

	//   bool evalWaveCut=false;
	TString wCutFName = Form("wave-Cut-%s.txt",out.Data());

	if (evalWaveCut) {

		std::ofstream cutf;
		//cutf.open(wCutFName.Data(), std::ofstream::out);//stampa su file.

		int SF=10;


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
				double time = (_tmax)/((float)dim);//in nano sec
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
						double time = (_tmax)/((float)dim);//
						X.clear();
						for (int i = 0; i < dim; ++i) { X.push_back(time *(i+1)); }
					}
					for (auto point : wd->getXOSCData()) {
						for (int channel=0; channel<=nMaxCh; channel++){
#endif
							bool isTrg=true;
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
				//std::cout<<"Fit of Intergeal SWF 10 fit results:"<<std::endl;
				for (int channel=0; channel<=nMaxCh; channel++){
					if (waveCut.find(channel)!=waveCut.end()) {                                             //fintanto che i canali non sono finiti scrivi su file.
						//std::cout<<channel<<"\t"<<waveCut[channel].first<<"\t"<<waveCut[channel].second<<std::endl;
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


			/*Ricerca dei picchi di noise sulle FFT e plot dei risultati*/
			TSpectrum *spectrum = new TSpectrum(10,1);

			/**array di errore*/
			std::vector<float> Xerror;
			std::vector<float> Yerror;


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
			//TVirtualFitter::SetDefaultFitter("Minuit2");

			float maxFreq[9]={1e+6, 10e+6, 300e+6, 500e+6, 1e+9, 2e+9, 1.5e+9, 700e+6, 1e+9};
			int maxNP=0.0;
			float ampl[]={0.0};
			float phase[]={0.0};
			float freq[]={0};
			TMinuit *gMinuit= new TMinuit();
			//dichiarazione struttura per i branches
			struct data0{
				float_t ampl0=0.0;
				float_t Eampl0=0.0;
				float_t phase0=0.0;
				float_t Ephase0=0.0;
				float_t freq0=0.0;
				float_t Efreq0=0.0;
			}spettro0;

			struct data1{
				float_t ampl1=0.0;
				float_t Eampl1=0.0;
				float_t phase1=0.0;
				float_t Ephase1=0.0;
				float_t freq1=0.0;
				float_t Efreq1=0.0;
			}spettro1;

			struct data2{
				float_t ampl2=0.0;
				float_t Eampl2=0.0;
				float_t phase2=0.0;
				float_t Ephase2=0.0;
				float_t freq2=0.0;
				float_t Efreq2=0.0;
			}spettro2;
			//Creazione file .root
			float dimX[10]={1,0.75,0.15,0.008,0.007,0.0016,1,0.0031,1,1};
			TFile *MyFile = new TFile(Form("TransFctData-%s.root",out.Data()),"RECREATE");
			if ( MyFile->IsOpen() ){  cout <<"file opened"<<endl;}
			TTree *tree =new TTree("tree","spectrum analysis");
			tree->Branch("branch0",&spettro0.ampl0,"ampl0/F:Eampl0/F:phase0/F:Ephase0/F:freq0/F:Efreq0/F");
			tree->Branch("branch1",&spettro1.ampl1,"ampl1/F:Eampl1/F:phase1/F:Ephase1/F:freq1/F:Efreq1/F");
			tree->Branch("branch2",&spettro2.ampl2,"ampl2/F:Eampl2/F:phase2/F:Ephase2/F:freq2/F:Efreq2/F");

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
					double time = (_tmax)/((float)dim);//
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
							std::cout<<"MaxCh "<<nMaxCh<<" dim "<<dim<<endl;
							double time = (_tmax)/((float)dim);
							cout<< " time is "<< time << " and tmax is "<< tmax<< endl;
							X.clear();
							for (int i = 0; i < dim; ++i) { X.push_back(time *(i+1)); }
							std::cout<<"Entries "<<jentry<<" Event n. "<<wd->getEventNumber()<<std::endl;
							maxNP=1.0/(maxFreq[wd->getRunHeader()->getRunNumber()-1]*1.1*time*2);
							cout<< " maxNP "<< maxNP<< endl;
						}
						int Ch_sr=2;
						float p0=0.;
						float p1=0.;
						float period=0.;
						float pos_zero=0.;
						int pos_zero_1=0.;
						int pos_zero1_2=0.;
						float pos_zero2=0.;
						//float zero_cut=1e-3;//run2
						//float zero_cut=1e-4;//run1
						//float zero_cut=8e-4;//run3
						float zero_cut=0.;
						TF1 *pol0=new TF1 ("pol0","x+[0]",0,1000);
						pol0->SetParameter(0,0);
						pol0->SetParLimits(0,0.0,1);
						pol0->SetNpx(1000);
						for (auto point : wd->getXOSCData()) {
							if (point.first == Ch_sr ) {
								Waves[Ch_sr].fillWave(point.second);

								//draw of signal ch2
								bool SaveCh2=false;
								if(SaveCh2){
									WavesCh2->cd();
									tmpCvsignal2.push_back( new TCanvas(Form("CvSignal2-Ch2_ev%d",jentry),Form("tmpSignal_1-Ch2_ev%d",jentry)) );
									tmpCvsignal2.back()->cd();
									TGraph* tmpsignalC2= new TGraph ( dim, &X[0], &Waves[Ch_sr].Y[0]);
									tmpsignalC2->GetXaxis()->SetTitle("time [s]");
									tmpsignalC2->SetTitle(Form("tmpSignal_afterFlt-Ch2_ev%d",jentry));
									tmpsignalC2->GetYaxis()->SetTitleOffset(1.4);
									tmpsignalC2->GetXaxis()->SetRangeUser(0,X[dim*dimX[runId]]);
									tmpsignalC2->GetYaxis()->SetTitle("Volt");
									tmpsignalC2->Draw("AL");
									tmpsignalC2->Fit(pol0,"","",0,X[dim*dimX[runId]*0.1]);
									tmpsignalC2->Fit(pol0,"","",0,X[dim*dimX[runId]]);
									tmpCvsignal2.back()->Write();
									theFile->cd("/");
									delete tmpsignalC2;
								}
								//zero_cut=TMath::Abs(pol0->GetParameter(0));
								//cout << "zero_cut "<<zero_cut<<endl;
								//ricerca del periodo //
								for(int ip=0;ip<Waves[Ch_sr].nPt();++ip){
									//cout<< " asse Y = "<< Waves[Ch_sr].Y[ip]<< endl;
									//cout<< " asse X = "<< X[ip]<<endl;
//									if( X[ip]<X[dim*dimX[runId]]){//run 1 if( X[ip]<0.05e-3,run2 X[ip]<4e-6)
//										if(pos_zero_1 ==0 && TMath::Abs(Waves[Ch_sr].Y[ip])<zero_cut){//<
//											pos_zero=X[ip];
//											pos_zero_1=ip;
//										}
//										if(pos_zero_1!=0 && ip>pos_zero_1 && pos_zero1_2==0 && TMath::Abs(Waves[Ch_sr].Y[ip])<zero_cut){
//											//							if((ip-pos_zero_1)>500){
//											if((ip-pos_zero_1)>maxNP){
//												pos_zero2=X[ip];
//												pos_zero1_2=ip;
//											}
//										}
//									}
									if((pos_zero_1!=0 && ip>pos_zero_1 && pos_zero1_2==0 )||(pos_zero_1 ==0) ){
										if(Waves[Ch_sr].Y[0]>0){
											if(pos_zero_1 ==0){
												if(Waves[Ch_sr].Y[ip]<0){
													pos_zero=X[ip];
													pos_zero_1=ip;
												}
											} else {
												if(Waves[Ch_sr].Y[ip]>0){
													if((ip-pos_zero_1)>maxNP){
														pos_zero2=X[ip];
														pos_zero1_2=ip;
													}
												}
											}
										}
										else if(Waves[Ch_sr].Y[0]<0){
											if(pos_zero_1 ==0){
													if(Waves[Ch_sr].Y[ip]>0){
														pos_zero=X[ip];
														pos_zero_1=ip;
													}
											}else {
												if(Waves[Ch_sr].Y[ip]<0){
													if((ip-pos_zero_1)>maxNP){
														pos_zero2=X[ip];
														pos_zero1_2=ip;
													}
												}
											}
										}
									} else break;
								}
							}
						}
						period=(pos_zero2-pos_zero)*2;
						//cout << " PERIODO " <<period<<endl;
						cout<< " canale " << Ch_sr<<" evento "<< jentry <<" primo_zero(t) "<< pos_zero<<" secondo _zero(t) "<<pos_zero2<<" posizione primo zero "<< pos_zero_1 << " posizione secondo zero "<<pos_zero1_2<< endl;
						Waves.clear();
						for (auto point : wd->getXOSCData()){
							for (int channel=0; channel<=nMaxCh; channel++){
#endif

								/*Esclusione dei canali di trigger.*/
								bool isTrg=true; //serve a controllare che il canale su cui facciamo l'analisi non sia quello di trigger.
								float max_value=0.0;
								float max_pos=0.;
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


									filterWaveBsl(Wffts[channel],realFltFFT,imgFltFFT);//filtro sulla baseline.                                        ///
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
									//if(!isTrg){

										if ( HstPerCh.find(channel)==HstPerCh.end() ){
											TDirectory *chDir = theFile->mkdir(Form("H-Ch_%d_signal",channel));
											chDir->cd();
											HstPerCh.insert(make_pair(channel,new hstPerCh(channel)));
											//	if (!evalWaveCut)  HstSF10PerCh.insert(make_pair(channel,new hstPerCh(channel,SF))); //make_pair funzione per creare una coppia.
											theFile->cd("/");
										}


									//}

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


									int	m=23;
									int	k=2;
									std::vector<float> tmpWSG_23=smoothSG(Waves[channel].Y,m,k);
									SmtSGWaves_23[channel].fillWave(tmpWSG_23);
									//---------------------------------------disegno dei grafici delle forme d'onda sperimentali con smooth---------------------------------
									bool saveWave=false; //true
									if (saveWave) {
										waveDir->cd();
										tmpCv.push_back( new TCanvas(Form("Cv-Ch%d_ev%d",channel,jentry),Form("tmpWave-Ch%d_ev%d",channel,jentry)) );
										tmpCv.back()->cd();
										tmpWaves.push_back( new TGraph ( dim, &X[0], &Waves[channel].Y[0]) );
										tmpWaves.back()->GetXaxis()->SetTitle("time [s]");
										tmpWaves.back()->SetTitle(Form("tmpWave-Ch%d_ev%d",channel,jentry));
										tmpWaves.back()->GetYaxis()->SetTitleOffset(1.4);
										tmpWaves.back()->GetYaxis()->SetTitle("Volt");
										tmpWaves.back()->Draw("AL");
										tmpCv.back()->Write();
										theFile->cd("/");
									}
									//--------------------------------------Disegno dei grafici della trasformata-------------------------------
									bool saveWaveFFT=false; //true
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




									//--------------------Disegno della wave di segnale: signal_1--------------------------------------------------------

									//ricerca max e min in wave_signal_1

									for(int ip=0;ip<Waves_signal_1[channel].nPt();++ip){
										if( X[ip]<X[dim*dimX[runId]]){//4e-6run2
											if(Waves_signal_1[channel].Y[ip]>max_value){
												max_value=Waves_signal_1[channel].Y[ip];
												max_pos=X[ip];
											}
										}
									}
									//cout<<" channel "<< channel<<" entry "<<jentry<< " valore max "<< max_value<< " max_pos "<<max_pos<<endl;
									p0=max_value;
									//derivata waveform
									float der=0.0;
									float mean1=0.0;
									float mean2=0.0;
									int mnOnPt=3;
									int rngDer=15;
									for(int ip=0;ip<mnOnPt;++ip){
										mean1+=Waves_signal_1[channel].Y[ip];
										mean2+=Waves_signal_1[channel].Y[ip+rngDer];
									}
									mean1/=(float)mnOnPt;
									mean2/=(float)mnOnPt;

									der=(mean2-mean1)/(X[rngDer]-X[0]);
									//cout<< " canale " << channel<<" evento "<< jentry<<" derivata "<< der<<endl;
									p1=1.0/period*1e-6;
									float p2=0.;
									//p2=(TMath::ASin(Waves_signal_1[channel].Y[0]/p0));

									//if (TMath::Abs(der)>1e5) {
										if(der>0){
											p2=-TMath::ASin(Waves_signal_1[channel].Y[0]/p0);

											if(Waves_signal_1[channel].Y[0]>0){
												p2+=TMath::TwoPi();
											}
											p2/=TMath::TwoPi();
										}
										if(der<0){
											p2=((TMath::ASin(Waves_signal_1[channel].Y[0]/p0))+(TMath::Pi()))/TMath::TwoPi();
										}
									//}
									cout<< " canale " << channel<<" evento "<< jentry<<" derivata "<< der<<" p0 "<<p0<<" p1 "<< p1 << " p2 "<<p2<<endl;
									bool savesignal_1=true;
									if (savesignal_1 ) {

										signal->cd();

										tmpCvsignal_1.push_back( new TCanvas(Form("CvSignal_1-Ch%d_ev%d",channel,jentry),Form("tmpSignal_1-Ch%d_ev%d",channel,jentry)) );
										tmpCvsignal_1.back()->cd();
										//tmpsignal_1.push_back( new TGraph ( dim, &X[0], &Waves_signal_1[channel].Y[0]) );
										TGraph* tmpsignal= new TGraph ( dim, &X[0], &Waves_signal_1[channel].Y[0]);
										/*tmpsignal_1.back()*/tmpsignal->GetXaxis()->SetTitle("time [s]");
										/*tmpsignal_1.back()*/tmpsignal->SetTitle(Form("tmpSignal_afterFlt-Ch%d_ev%d",channel,jentry));
										/*tmpsignal_1.back()*/tmpsignal->GetYaxis()->SetTitleOffset(1.4);
										/*tmpsignal_1.back()*/tmpsignal->GetXaxis()->SetRangeUser(0,X[dim*dimX[runId]]);
										/*tmpsignal_1.back()*/tmpsignal->GetYaxis()->SetTitle("Volt");
										/*tmpsignal_1.back()*/tmpsignal->Draw("AL");
										sinFit->SetParameter(0,p0);
										sinFit->SetParameter(1,p1);
										sinFit->SetParameter(2,p2);

										//int statusFit=gMinuit->GetStatus();
										//cout<< "status fit "<< statusFit<< endl;
										//per il run 1 X[dim*0.75]
										//ROOT::Math::MinimizerOptions::SetDefaultMinimizer("Minuit2","Migrad");
										//gMinuit->mnsimp();//Minimization using the simplex method of Nelder and Mead.
										gMinuit->Migrad();
										//gMinuit->SetMaxIterations(1000);

										tmpsignal->Fit(sinFit,"M","",0,X[dim*dimX[runId]*0.1]);
										tmpsignal->Fit(sinFit,"M","",0,X[dim*dimX[runId]*0.2]);
										tmpsignal->Fit(sinFit,"M","",0,X[dim*dimX[runId]*0.5]);
										//tmpsignal->Fit(sinFit,"M","",0,X[dim*dimX[runId]*0.5]);
//										tmpsignal->Fit(sinFit,"M","",0,X[dim*dimX[runId]*0.8]);
										tmpsignal->Fit(sinFit,"M","",0,X[dim*dimX[runId]]);//QUI
										//TFitResultPtr result=/*tmpsignal_1.back()*/tmpsignal->Fit(sinFit,"S","",0,X[dim*0.15]);
										//bool fit_stat=result->IsValid();
										//cout<< " fit stat " << result->IsValid() << " matrix " << result->CovMatrixStatus()<< endl;

										/*Double_t fmin;
										Double_t fedm;
										Double_t errdef;
										Int_t  npari;
										Int_t  nparx;
										Int_t istat ;
										gMinuit->mnstat(fmin,fedm,errdef,npari,nparx,istat);*/
										//cout<<"fmin" <<fmin<< " "<<fedm<< " "<<errdef<< " "<<npari<< " "<<nparx<< " "<<istat<< " "<<endl;
//										if (fit_stat==0){
//
//											sinFit->SetParameter(0,p0);
//											sinFit->SetParameter(1,p1);
//											sinFit->SetParameter(2,p2);
//											/*tmpsignal_1.back()*/tmpsignal->Fit(sinFit,"M","",0,X[dim*0.15]);//X[dim*0.75]);
//											///*tmpsignal_1.back()*/tmpsignal->Fit(sinFit,"M","",0,X[dim*0.18]);//X[dim*0.75]);//r1
//											///*tmpsignal_1.back()*/tmpsignal->Fit(sinFit,"M","",0,X[dim*0.18]);//X[dim*0.75]);//r1
//
//										} // da if(fit stat==0)fino alla fine della parentesi è per il run 1.

										if(sinFit->GetParameter(0)<0.01||sinFit->GetParameter(1)<2.||sinFit->GetChisquare()>1){
												sinFit->SetParameter(0,p0);
												sinFit->SetParameter(1,p1);
												sinFit->SetParameter(2,p2);

												tmpsignal->Fit(sinFit,"","",0,X[dim*dimX[runId]*0.01]);
												tmpsignal->Fit(sinFit,"","",0,X[dim*dimX[runId]*0.02]);
												tmpsignal->Fit(sinFit,"M","",0,X[dim*dimX[runId]*0.02]);
												tmpsignal->Fit(sinFit,"M","",0,X[dim*dimX[runId]*0.03]);
												tmpsignal->Fit(sinFit,"M","",0,X[dim*dimX[runId]*0.04]);
												tmpsignal->Fit(sinFit,"M","",0,X[dim*dimX[runId]*0.05]);
//												tmpsignal->Fit(sinFit,"M","",0,X[dim*dimX[runId]*0.07]);
//												tmpsignal->Fit(sinFit,"M","",0,X[dim*dimX[runId]*0.09]);
												tmpsignal->Fit(sinFit,"M","",0,X[dim*dimX[runId]*0.1]);
												tmpsignal->Fit(sinFit,"M","",0,X[dim*dimX[runId]*0.15]);
												tmpsignal->Fit(sinFit,"M","",0,X[dim*dimX[runId]*0.2]);
//												tmpsignal->Fit(sinFit,"M","",0,X[dim*dimX[runId]*0.3]);
//												tmpsignal->Fit(sinFit,"M","",0,X[dim*dimX[runId]*0.4]);
												tmpsignal->Fit(sinFit,"M","",0,X[dim*dimX[runId]*0.5]);
												tmpsignal->Fit(sinFit,"M","",0,X[dim*dimX[runId]*0.8]);
//												tmpsignal->Fit(sinFit,"M","",0,X[dim*dimX[runId]*0.9]);
//												tmpsignal->Fit(sinFit,"M","",0,X[dim*dimX[runId]*0.95]);
												tmpsignal->Fit(sinFit,"M","",0,X[dim*dimX[runId]]);
												cout<< " controllo "<< channel << " "<< jentry <<endl;

										}

										tmpCvsignal_1.back()->Write();
										theFile->cd("/");
										delete tmpsignal;
									}



									//cout<< " canale "<< channel<< " evento "<< jentry << " fase "<< phase[channel]<< " ampl "<< ampl[channel]<< " freq "<< freq[channel]<< endl;
									float ndf=sinFit->GetNDF();
									float chi=sinFit->GetChisquare();
									if(/*chi>1||*/sinFit->GetParameter(0)<0.01){cout << " entry not fit "<< jentry<<" canale "<<channel<<endl;}
									//cout<< " ch "<<channel<<" entry "<<jentry<<" chi "<<chi<< " ndf "<<ndf<<endl;
									//if(chi/ndf<0.02){//1.5 run 1

									if(chi/ndf<0.02){
										if (channel ==0){
											spettro0.ampl0=sinFit->GetParameter(0);
											spettro0.Eampl0=sinFit->GetParError(0);
											//cout<< " ampiezza "<<spettro0.ampl0<<"errore_ampl"<< spettro0.Eampl0<<endl;
											spettro0.phase0=sinFit->GetParameter(2);
											spettro0.Ephase0=sinFit->GetParError(2);
											spettro0.freq0=sinFit->GetParameter(1);
											spettro0.Efreq0=sinFit->GetParError(1);

										}

										if (channel ==1){
											//if(sinFit->GetParameter(0)<0.01){cout << " entry not fit "<< jentry<<" canale "<<channel<<endl;}

											spettro1.ampl1=sinFit->GetParameter(0);
											spettro1.Eampl1=sinFit->GetParError(0);
											spettro1.phase1=sinFit->GetParameter(2);
											spettro1.Ephase1=sinFit->GetParError(2);
											spettro1.freq1=sinFit->GetParameter(1);
											spettro1.Efreq1=sinFit->GetParError(1);
										}

										if (channel ==2){
											spettro2.ampl2=sinFit->GetParameter(0);
											spettro2.Eampl2=sinFit->GetParError(0);
											spettro2.phase2=sinFit->GetParameter(2);
											spettro2.Ephase2=sinFit->GetParError(2);
											spettro2.freq2=sinFit->GetParameter(1);
											spettro2.Efreq2=sinFit->GetParError(1);
											if(spettro2.Efreq2>1){cout << "evento scemo"<<jentry<<endl;}

										}
				float tmpEG=(1/spettro1.ampl1)*TMath::Sqrt((spettro0.ampl0*spettro0.ampl0/spettro1.ampl1*spettro1.ampl1)*(spettro1.Eampl1*spettro1.Eampl1)+(spettro0.Eampl0*spettro0.Eampl0));
				if(tmpEG>1.){cout << " err "<<jentry<<" "<< channel<<endl; }
				float tmpG=	spettro0.ampl0/spettro1.ampl1;
				//if(tmpG>6){cout << " not correct gain"<< jentry<<" "<< channel<<endl; }
									}



									//((hstPerCh*)HstPerCh[channel])->hder->Fill(der);
									((hstPerCh*)HstPerCh[channel])->hchiFT->Fill(chi/ndf);
									((hstPerCh*)HstPerCh[channel])->hP0->Fill(sinFit->GetParameter(0));

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

		tree->Fill();

	} // end of loop on entry (second loop)

					//hDrvDistribution.back()->Draw();

					theFile->cd();
					theFile->Write();
					theFile->Close();

					MyFile->cd();
					MyFile->Write();
					/*tree->Scan("ampl2:phase2:freq2");*/
					//tree->Scan("ampl0:Eampl0:phase0:Ephase0:freq0:Efreq0");
					//tree->Scan("ampl2:phase2:freq2:Efreq2");
					//tree->Print();
					MyFile->Close();

}
