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

//#include <map>
#include <fstream>
#include <string>

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//////////////Define Data Containers//////////////

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

vector<float>X;


// Temporary container for tmp Graphs and Canvas
std::vector<TGraph *> tmpWaves;
std::vector<TCanvas *> tmpCv;
std::vector<TGraph *> tmpFFTAmpl;
std::vector<TGraph *> tmpFFTFltAmpl;
std::vector<TGraph *> tmpFFTphi;
std::vector<TCanvas *> tmpCvFFT;
std::vector<TGraph *> tmpFltWaves;
std::vector<TCanvas *> tmpCvFlt;

//void doChDiff(int evN, std::pair<int,int> &chDiff, std::map<int, bool> &isFull, diffWvCont &diffWaves, WvCont &Waves, WvCont &Smt4Waves, int SF, hDiffCont &hDiff, hDiffCont &hDiffPed, TFile *theFile);

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//////////////MAIN LOOP function//////////////

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void read_data::Loop(Char_t *output, Int_t eventn, Bool_t evalWaveCut)
{

	tmax=_tmax;  //to fix compilation problem

	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	//                                                              // Define Histograms //                                       //

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	char *basename(char *path);
	cout << "Basename " << basename(output) << endl;

	Char_t outChar[500];
	TString out = basename(output);
	out.ReplaceAll(".root","");
#ifndef _OSC
	sprintf(outChar,"histosTB_%s.root",out.Data());
#else
	sprintf(outChar,"histosOSC_%s.root",out.Data());
#endif
	cout << "The output file is " << outChar << endl;

	TFile * theFile = new TFile(outChar,"RECREATE");

	std::map<int, hstPerCh *> HstPerCh;
	std::map<int, hstPerCh *> HstSF10PerCh;

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

	theFile->cd("/");
	TDirectory *waveDir = theFile->mkdir("Waves");
	theFile->cd("/");
	TDirectory *waveDiffDir = theFile->mkdir("WavesDiff");
	theFile->cd("/");
	TDirectory *waveFFTDir = theFile->mkdir("WavesFFT");
	theFile->cd("/");
	TDirectory *waveFltDir = theFile->mkdir("WavesFlt");
	theFile->cd("/");

	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	//                                               //  Loop on entries //                                                  //

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	std::vector<int> trigCh;
	trigCh.clear();
#ifndef _OSC
	trigCh.push_back(32);
	trigCh.push_back(33);
#else
	trigCh.push_back(7);
	trigCh.push_back(8);
#endif

	std::map<int, std::pair<float, float> > waveCut;
	waveCut.clear();
	int nSig=3;
	std::map<int, int> nSelEv;
	nSelEv.clear();

	if (fChain == 0) return;

	Long64_t nentries = fChain->GetEntriesFast();
	cout << "Number of entries in the tree for real data is= " << nentries << endl;
	//   Long64_t nbytes = 0, nb = 0;

	Int_t frstEv=0;
	Int_t lastEv=nentries;
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
		cutf.open(wCutFName.Data(), std::ofstream::out);

		int SF=10;

		//	   std::map<int, TF1*> tmpFitRes;

		for (Long64_t jentry=frstEv; jentry<lastEv;jentry++) {

			Long64_t ientry = LoadTree(jentry);
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
#else
			if(jentry==frstEv) {
				nMaxCh=wd->getRunHeader()->getNumberOfDevicesInRun()*8-1;
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
							TDirectory *chDir = theFile->mkdir(Form("H-Ch%d",channel));
							chDir->cd();
							HstPerCh.insert(make_pair(channel,new hstPerCh(channel)));
							HstSF10PerCh.insert(make_pair(channel,new hstPerCh(channel,SF)));
							theFile->cd("/");
						}

						//cout << "Channel is= " << channel << endl;
						//Waves[channel].fillWave(point.second.size(),point.second.data());
						Waves[channel].fillWave(point.second);
						std::vector<float> tmpWSF10=smooth(Waves[channel].Y,SF);
						Smt10Waves[channel].fillWave( tmpWSF10 );

						((hstPerCh*)HstSF10PerCh[channel])->hBsl->Fill(((wave)Smt10Waves[channel]).bsln);
						((hstPerCh*)HstSF10PerCh[channel])->hInteg->Fill(((wave)Smt10Waves[channel]).integ);
						((hstPerCh*)HstSF10PerCh[channel])->hIntegN->Fill(((wave)Smt10Waves[channel]).nInteg());
						((hstPerCh*)HstSF10PerCh[channel])->hRms->Fill(((wave)Smt10Waves[channel]).rms);
						((hstPerCh*)HstSF10PerCh[channel])->hMaxV->Fill(((wave)Smt10Waves[channel]).max);
						((hstPerCh*)HstSF10PerCh[channel])->hMaxVN->Fill(((wave)Smt10Waves[channel]).nMax());

						bool saveWave=false;
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

		for (int channel=0; channel<=nMaxCh; channel++){
			if (HstSF10PerCh.find(channel)!=HstSF10PerCh.end()){
				((hstPerCh*)HstSF10PerCh[channel])->hMaxVN->Fit("gaus","","",0,0.005);
				//tmpFitRes.insert( make_pair(channel, ((hstPerCh*)HstSF10PerCh[channel])->hMaxVN->GetFunction("gaus") ) );
				TF1 *tmpFit = ((hstPerCh*)HstSF10PerCh[channel])->hMaxVN->GetFunction("gaus") ;
				float tmpMean = tmpFit->GetParameter(1);
				float tmpSigma = tmpFit->GetParameter(2);
				waveCut.insert( make_pair(channel, make_pair(tmpMean,tmpSigma) ) );
			}
		}
		std::cout<<"Fit of Intergeal SWF 10 fit results:"<<std::endl;
		for (int channel=0; channel<=nMaxCh; channel++){
			if (waveCut.find(channel)!=waveCut.end()) {
				std::cout<<channel<<"\t"<<waveCut[channel].first<<"\t"<<waveCut[channel].second<<std::endl;
				cutf<<channel<<"\t"<<waveCut[channel].first<<"\t"<<waveCut[channel].second<<std::endl;
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

	std::vector< std::pair<int,int> > chToDiff;
#ifndef _OSC
	chToDiff.push_back( make_pair(1,2) );
#else
#endif

	hDiffCont hDiff;
	hDiffCont hDiffPed;

	for ( auto chsDiff : chToDiff ) {
		TDirectory *diffDir = theFile->mkdir(Form("DiffCh%d%d",chsDiff.first,chsDiff.second));
		diffDir->cd();
		hDiff.insert( make_pair(chsDiff, hstDiffCh(chsDiff.first,chsDiff.second) ) );
		hDiffPed.insert( make_pair(chsDiff, hstDiffCh(chsDiff.first,chsDiff.second,true) ) );
		theFile->cd("/");
	}

	WvCont Waves;
	WvCont Smt10Waves;
	WvCont Smt4Waves;
	diffWvCont diffWaves;
	fftCont Wffts;
	WvCont FltWaves;

	for (Long64_t jentry=frstEv; jentry<lastEv;jentry++) {
		if ((jentry-frstEv)%200==0) { std::cout<<"Processing event: "<<jentry<<std::endl; }
		Long64_t ientry = LoadTree(jentry);
		if (ientry < 0) break;

		//     nb = fChain->GetEntry(jentry);   nbytes += nb;
		fChain->GetEntry(jentry);
		// if (Cut(ientry) < 0) continue;

		std::map<int, bool> isFull;
		for (int ich=0; ich<=nMaxCh; ++ich ) { isFull[ich]=false; }

		Waves.clear();
		Smt10Waves.clear();
		Smt4Waves.clear();
		Wffts.clear();
		FltWaves.clear();
		for ( auto chsDiff : chToDiff ) { diffWaves[chsDiff].clear(); }

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
		}
		for (auto point : wd->getXOSCData()) {
			for (int channel=0; channel<=nMaxCh; channel++){
#endif
				bool isTrg=false;
				for(auto trgCh : trigCh) { if (channel==trgCh) { isTrg=true; } }
				if (point.first == channel) {
					//cout << "Channel is= " << channel << endl;

					//		 Waves[channel].fillWave(point.second.size(),point.second.data());
					Waves[channel].fillWave(point.second);
//					std::cout<<"nPt "<<Waves[channel].nPt()<<std::endl;
					FFT(Waves[channel],Wffts[channel]);
					double *realFltFFT = new double[Waves[channel].nPt()];
					double *imgFltFFT = new double[Waves[channel].nPt()];
#ifndef _OSC
					filterWave(Wffts[channel],realFltFFT,imgFltFFT);
#else
					filterWave(Wffts[channel],realFltFFT,imgFltFFT,50e+6,200e+6);
#endif
					InverseFFT(realFltFFT,imgFltFFT,Waves[channel].nPt(),FltWaves[channel]);
					double *amplFltFFT = new double[Waves[channel].nPt()];
					for (int ip=0; ip<Waves[channel].nPt()/2; ++ip) {
//						std::cout<<"omega "<<Wffts[channel].omega[ip]<<" Flt: real "<<realFltFFT[ip]<<" img "<<imgFltFFT[ip]<<std::endl;
						amplFltFFT[ip]=TMath::Sqrt(realFltFFT[ip]*realFltFFT[ip]+imgFltFFT[ip]*imgFltFFT[ip]);
					}
					for (int ip=Waves[channel].nPt()/2; ip<Waves[channel].nPt(); ++ip) {
						amplFltFFT[ip]=0.0;
					}

					if ( HstPerCh.find(channel)==HstPerCh.end() ) {
						TDirectory *chDir = theFile->mkdir(Form("H-Ch%d",channel));
						chDir->cd();
						HstPerCh.insert(make_pair(channel,new hstPerCh(channel)));
						if (!evalWaveCut)  HstSF10PerCh.insert(make_pair(channel,new hstPerCh(channel,SF)));
						theFile->cd("/");
					}

					((hstPerCh*)HstPerCh[channel])->hBsl->Fill(((wave)Waves[channel]).bsln);
					((hstPerCh*)HstPerCh[channel])->hInteg->Fill(((wave)Waves[channel]).integ);
					((hstPerCh*)HstPerCh[channel])->hIntegN->Fill(((wave)Waves[channel]).nInteg());
					((hstPerCh*)HstPerCh[channel])->hRms->Fill(((wave)Waves[channel]).rms);
					((hstPerCh*)HstPerCh[channel])->hMaxV->Fill(((wave)Waves[channel]).max);
					((hstPerCh*)HstPerCh[channel])->hMaxVN->Fill(((wave)Waves[channel]).nMax());

					std::vector<float> tmpWSF10=smooth(Waves[channel].Y,SF);
					Smt10Waves[channel].fillWave( tmpWSF10 );

					if (!evalWaveCut) {
						((hstPerCh*)HstSF10PerCh[channel])->hBsl->Fill(((wave)Smt10Waves[channel]).bsln);
						((hstPerCh*)HstSF10PerCh[channel])->hInteg->Fill(((wave)Smt10Waves[channel]).integ);
						((hstPerCh*)HstSF10PerCh[channel])->hIntegN->Fill(((wave)Smt10Waves[channel]).nInteg());
						((hstPerCh*)HstSF10PerCh[channel])->hRms->Fill(((wave)Smt10Waves[channel]).rms);
						((hstPerCh*)HstSF10PerCh[channel])->hMaxV->Fill(((wave)Smt10Waves[channel]).max);
						((hstPerCh*)HstSF10PerCh[channel])->hMaxVN->Fill(((wave)Smt10Waves[channel]).nMax());
					}

					SF=4;
					std::vector<float> tmpWSF4=smooth(Waves[channel].Y,SF);
					Smt4Waves[channel].fillWave( tmpWSF4 );

					bool saveWave=true;
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
						SF=10;
						tmpWaves.push_back( new TGraph ( dim-2*SF, &X[SF], &Smt10Waves[channel].Y[0]) );
						tmpWaves.back()->SetLineColor(kBlue);
						tmpWaves.back()->Draw("Lsame");
						SF=4;
						tmpWaves.push_back( new TGraph ( dim-2*SF, &X[SF], &Smt4Waves[channel].Y[0]) );
						tmpWaves.back()->SetLineColor(kRed);
						tmpWaves.back()->Draw("Lsame");
						tmpCv.back()->Write();
						theFile->cd("/");
					}
					bool saveWaveFFT=true;
					if (saveWaveFFT) {

						waveFFTDir->cd();

						tmpCvFFT.push_back( new TCanvas(Form("CvFFT-Ch%d_ev%d",channel,jentry),Form("tmpWaveFFT-Ch%d_ev%d",channel,jentry)) );
						tmpCvFFT.back()->Divide(1,3);
						tmpCvFFT.back()->SetLogx();
						tmpCvFFT.back()->cd(1)->SetLogy();
						tmpFFTAmpl.push_back( new TGraph ( dim/2, &Wffts[channel].omega[0], &Wffts[channel].ampl[0]) );
						tmpFFTAmpl.back()->GetXaxis()->SetTitle("#omega [Hz]");
						tmpFFTAmpl.back()->SetTitle(Form("tmpWaveFFT-Ampl-Ch%d_ev%d",channel,jentry));
						tmpFFTAmpl.back()->GetYaxis()->SetTitleOffset(1.4);
						tmpFFTAmpl.back()->GetYaxis()->SetTitle("Amplitude");
						tmpFFTAmpl.back()->Draw("AL");
						tmpCvFFT.back()->cd(2);
						tmpFFTphi.push_back( new TGraph ( dim/2, &Wffts[channel].omega[0], &Wffts[channel].phi[0]) );
						tmpFFTphi.back()->GetXaxis()->SetTitle("#omega [Hz]");
						tmpFFTphi.back()->SetTitle(Form("tmpWaveFFT-Phi-Ch%d_ev%d",channel,jentry));
						tmpFFTphi.back()->GetYaxis()->SetTitleOffset(1.4);
						tmpFFTphi.back()->GetYaxis()->SetTitle("#Phi [rad]");
						tmpFFTphi.back()->Draw("AL");
						tmpCvFFT.back()->cd(3)->SetLogy();
						tmpFFTAmpl.push_back( new TGraph ( dim/2, &Wffts[channel].omega[0], amplFltFFT) );
						tmpFFTAmpl.back()->GetXaxis()->SetTitle("#omega [Hz]");
						tmpFFTAmpl.back()->SetTitle(Form("tmpWaveFltFFT-Ampl-Ch%d_ev%d",channel,jentry));
						tmpFFTAmpl.back()->GetYaxis()->SetTitleOffset(1.4);
						tmpFFTAmpl.back()->GetYaxis()->SetTitle("Amplitude");
						tmpFFTAmpl.back()->Draw("AL");
						tmpCvFFT.back()->Write();
						theFile->cd("/");
					}
					bool saveFltWave=true;
					if (saveFltWave) {

						waveFltDir->cd();

						tmpCvFlt.push_back( new TCanvas(Form("CvFlt-Ch%d_ev%d",channel,jentry),Form("tmpFltWave-Ch%d_ev%d",channel,jentry)) );
						tmpCvFlt.back()->cd();
						tmpFltWaves.push_back( new TGraph ( dim, &X[0], &FltWaves[channel].Y[0]) );
						tmpFltWaves.back()->GetXaxis()->SetTitle("time [ns]");
						tmpFltWaves.back()->SetTitle(Form("tmpFltWave-Ch%d_ev%d",channel,jentry));
						tmpFltWaves.back()->GetYaxis()->SetTitleOffset(1.4);
						tmpFltWaves.back()->GetYaxis()->SetTitle("Volt");
						tmpFltWaves.back()->Draw("AL");
						//    					 SF=10;
						//    					 tmpWaves.push_back( new TGraph ( dim-2*SF, &X[SF], &Smt10Waves[channel].Y[0]) );
						//    					 tmpWaves.back()->SetLineColor(kBlue);
						//    					 tmpWaves.back()->Draw("Lsame");
						//    					 SF=4;
						//    					 tmpWaves.push_back( new TGraph ( dim-2*SF, &X[SF], &Smt4Waves[channel].Y[0]) );
						//    					 tmpWaves.back()->SetLineColor(kRed);
						//    					 tmpWaves.back()->Draw("Lsame");
						tmpCvFlt.back()->Write();
						theFile->cd("/");
					}

					if ( Smt10Waves[channel].nMax()>( waveCut[channel].first + nSig*waveCut[channel].second ) ) {
						isFull[channel]=true;
						std::cout<<"Ev "<<jentry<<" Ch "<<channel<<" Max "<<Smt10Waves[channel].nMax()<<" cut "<<( waveCut[channel].first + nSig*waveCut[channel].second )<<" bsl "<<Smt10Waves[channel].bsln<<" max "<<Smt10Waves[channel].max<<std::endl;
						if (nSelEv.find(channel)==nSelEv.end()) {
							nSelEv[channel]=1;
						} else { ++nSelEv[channel]; }
					}

					delete [] realFltFFT;
					delete [] imgFltFFT;
					delete [] amplFltFFT;

				}// if point

			} // end of channel loop

		}// auto point

		for ( auto chsDiff : chToDiff ) {
			doChDiff(jentry,chsDiff,isFull,diffWaves,Waves,Smt4Waves,4,hDiff,hDiffPed,theFile);
		}

	} // end of loop on entry (second loop)

	bool prtNSelEv=false;
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

void read_data::doChDiff(int evN, std::pair<int,int> &chDiff, std::map<int, bool> &isFull, diffWvCont &diffWaves, WvCont &Waves, WvCont &Smt4Waves, int SF, hDiffCont &hDiff, hDiffCont &hDiffPed, TFile *theFile){

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

		bool saveWave=true;
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

}
