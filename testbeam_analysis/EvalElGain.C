#include <stdio.h>
#include <dirent.h>
#include <Riostream.h>
#include <map>
#include <vector>
#include <iterator>
#include "TGraphErrors.h"
#include "TGraph.h"
#include "TCanvas.h"
#include <fstream>      // std::fstream


void EvalElGain(Int_t rNum=4, TString fOutName="") {

	TFile *fIn= TFile::Open(Form("histosTB_run_%d.root",rNum));
	TH1F *hMax = (TH1F*) fIn->Get("H-Ch6_signal/hMaxVNInR_ch6");
	TH1F *hIntRt = (TH1F*) fIn->Get("hIntRatio");
	TH1F *hMaxRt = (TH1F*) fIn->Get("hMaxRatio");
	//hIntRt->Rebin(10);
	//hMaxRt->Rebin(6);
	gStyle->SetOptFit(1);
	TF1 *fg= new TF1("fg","gaus");
	TCanvas *cvmax = new TCanvas();
	hMax->Draw();
	cvmax->Print(Form("/lustrehome/bdanzi/offline_analysis/testbeam_analysis/Plots/histosTB_run_4.root/cvmax_Ch%d.png",6), "png");
	TCanvas *cvIntRt = new TCanvas();
	//hIntRt->Draw();
	//hIntRt->Fit(fg,"","",0,13);
	std::ofstream fOut;
	if (!fOutName.IsNull()) { fOut.open(fOutName.Data(), std::ofstream::out | std::ofstream::app);}
	if (fOut.is_open()) {
		fOut<<rNum<<"\t"<<hMax->GetMean()<<"\t"<<hMax->GetMeanError()<<"\t"<<fg->GetParameter(1)<<"\t"<<fg->GetParError(1)<<"\t"<<fg->GetParameter(2)<<"\t"<<fg->GetParError(2);
	}
	TCanvas *cvMaxRt = new TCanvas();
	//hMaxRt->Draw();
	//hMaxRt->Fit(fg);
	if (fOut.is_open()) {
		fOut<<"\t"<<fg->GetParameter(1)<<"\t"<<fg->GetParError(1)<<"\t"<<fg->GetParameter(2)<<"\t"<<fg->GetParError(2)<<std::endl;
	}
	fOut.close();
}
