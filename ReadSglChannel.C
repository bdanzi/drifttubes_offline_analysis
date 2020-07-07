#include <TStyle.h>
#include <TCanvas.h>
#include <TGraph.h>
#include<TVector.h>
#include "TF1.h"
#include "TFitResult.h"
#include "TMinuit.h"
#include "TGraphErrors.h"
#include <fstream>
#include <string>
#include "TSpectrum.h"
#include "TAxis.h"

#include "TMath.h"
#include "TH1F.h"
#include "iostream"
#include "TFile.h"
#include "TLegend.h"

using namespace std;

void ReadSglChannel(TString fname="", int channel=23, int ev=10){


	TFile *file = new TFile(fname.Data(),"read");
	TCanvas *max=new TCanvas("max","max",1000,1000);
	TCanvas *bsl=new TCanvas("bsl","bsl",1000,1000);
	for(int i =0; i<=ev; ++i){
		TGraph *h1=(TGraph*)file->Get(Form("signal_Afterflt/CvSignal_1-Ch%i_ev%i",channel,i));
		if (h1==0x0) { continue; }
		h1->Draw();
	}
	for(int i =0; i<=23; ++i){
		TH1F *h2=(TH1F*)file->Get(Form("H-Ch%i_signal/hMaxVN_ch%i",channel,i));
		if (h2==0x0) { continue; }
		max->cd();
		h2->Draw();
		TH1F *h3=(TH1F*)file->Get(Form("H-Ch%i_signal/hBsl_ch%i",channel,i));
		if (h3==0x0) { continue; }
		bsl->cd();
		h3->Draw();
	}



}
