//fit integrali e pezzo di codice per creare txt
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
#include "TFitResult.h"
#include <iostream>

using namespace std;

void integ_fit(TString fname="/lustrehome/bdanzi/offline_analysis/testbeam_analysis/histosTB_run_4.root", TString sname = "/lustrehome/bdanzi/offline_analysis/testbeam_analysis/histosTB_run_4.root", TString fOutName=""){ //"/mnt/c/file_txt/ly2-3-4_90-10/1650.txt"

  ofstream myfile;
  if (!fOutName.IsNull()) {
    myfile.open (fOutName.Data());
  }
  
	gStyle->SetOptFit(1);
	gStyle->SetStatX(1);

	TFile *file = new TFile(fname.Data(),"read"); //NOISE
	TFile *file2= new TFile(sname.Data(),"read"); //SIGNAL

	TF1 *fGaus = new TF1("fGaus", "gaus");
	fGaus->SetNpx(1000);
	TF1 *fSig = new TF1("fSig", "gaus(0)+landau(3)");
	fSig->SetNpx(1000);
	fSig->SetParLimits(0,0,10000);
	fSig->SetParLimits(1,-2,2);
	fSig->SetParLimits(2,0,1);
	fSig->SetParLimits(3,0,1000);
	fSig->SetParLimits(4,0,5);
	fSig->SetParLimits(5,0,10);

	TF1 *fSig2 = new TF1("fSig2", "gaus(0)+gaus(3)+landau(6)");
	fSig2->SetNpx(1000);
	fSig2->SetParLimits(0,0,10000);
	fSig2->SetParLimits(1,-2,2);
	fSig2->SetParLimits(2,0,1);
	fSig2->SetParLimits(3,0,1000);
	fSig2->SetParLimits(4,-2,3);
	fSig2->SetParLimits(5,0,4);
	fSig2->SetParLimits(6,0,1000);
	fSig2->SetParLimits(7,0,5);
	fSig2->SetParLimits(8,0,10);
	fSig2->SetLineColor(kRed);
	
	/* TF1  *fLandau=new TF1 ("fLandau", "landau"); */
	/* fLandau->SetNpx(1000); */
	//TF1 *total=new TF1("total,"fGaus_2+fLandau");

	TCanvas *c1 = new TCanvas ("c1","Integral for noise");
	TCanvas *c2 = new TCanvas ("c2","Integral for signal");
	
	Int_t nPresCh=0;
	Int_t nPresCH_1=0;
	Float_t chn[12]={0,0,0,0,0,0,0,0,0,0,0,0};
	Float_t meanNoise[12]={0,0,0,0,0,0,0,0,0,0,0,0};
	Float_t sgmNoise[12]={0,0,0,0,0,0,0,0,0,0,0,0};
	/* Float_t meanSignal[12]={0,0,0,0,0,0,0,0,0,0,0,0}; */
	/* Float_t sgmSignal[12]={0,0,0,0,0,0,0,0,0,0,0,0}; */
	Float_t meanL[12]={0,0,0,0,0,0,0,0,0,0,0,0};
	Float_t sigmaL[12]={0,0,0,0,0,0,0,0,0,0,0,0};
	Float_t errMean[12]={0,0,0,0,0,0,0,0,0,0,0,0};
	Float_t errSigma[12]={0,0,0,0,0,0,0,0,0,0,0,0};
		

	c1->cd()->SetLogy();
	for(int i =4; i<=14; ++i){
     TH1F *h1=(TH1F*)file->Get(Form("H-Ch%i_signal/hIntegNInR_ch%i",i,i));// NOISE	    
		if (h1==0x0) { continue; }
		h1->Rebin(2);
		h1->SetLineColor(kBlue);
		h1->GetXaxis()->SetTitle("Integral[V]");
		h1->GetYaxis()->SetTitle("Entries");
		h1->Draw();
		h1->Fit(fGaus,"","",-8,8);
		gPad->Update();
		
		TPaveStats *stat1 = (TPaveStats*)h1->GetListOfFunctions()->FindObject("stats");
		stat1->SetY1NDC(0.5);    //set new y start position
		stat1->SetY2NDC(0.9);    //set new y end position
		stat1->SetX1NDC(0.1);    //set new x start position
		stat1->SetX2NDC(0.4);    //set new x end position
		stat1->SetTextColor(kBlue);
		
	       	/* TLegend *leg = new TLegend(0.45, 0.9, 0.6, 0.8); */
		/* leg->SetFillColor(0); */
		/* //leg->SetHeader("Histogram"); */
		/* leg->AddEntry(h1, "Noise Integral", "l"); */
		/* leg->Draw(); */
		c1->Print(Form("/lustrehome/bdanzi/offline_analysis/testbeam_analysis/Plots/histosTB_run_4.root/H%i_Integral_noise.png",i), "png");
		chn[nPresCh]=i;
		meanNoise[nPresCh]=fGaus->GetParameter(1); //noise
		sgmNoise[nPresCh]=fGaus->GetParameter(2);
		/* ratioSN_noise[nPresCh]=meanMinFlt[nPresCh]/sgmMinFlt[nPresCh]; */
		/* ratioSN_signal[nPresCh]=meanMinOrg[nPresCh]/sgmMinOrg[nPresCh]; */
		/* ratioMean[nPresCh]=meanSignal[nPresCh]/meanNoise[nPresCh]; */
		/* ratioSigma[nPresCh]=sgmSignal[nPresCh]/sgmNoise[nPresCh]; */
		++nPresCh;

	}
 
   c2->cd()->SetLogy();
	nPresCh=0;
	for(int i =4; i<=14; ++i){
  		TH1F *h2=(TH1F*)file2->Get(Form("H-Ch%i_signal/hIntegNInRC2_ch%i",i,i));//SIGNAL
		if (h2==0x0) { continue; }
		h2->Rebin(2);
		h2->SetLineColor(kGreen+2);
		h2->Draw();
//       	fSig->SetParameters(100,0,0.5,70,2.5,0.6);
       	fSig2->SetParameters(1000,0,0.5,20,0.5,1,70,2.5,0.6);
		fSig2->FixParameter(1,meanNoise[nPresCh]);
		fSig2->FixParameter(2,sgmNoise[nPresCh]);
//		fSig2->FixParameter(3,0);
//		fSig2->FixParameter(4,0);
//		fSig2->FixParameter(5,1);
       	h2->Fit(fSig2,"L");
       	//fSig2->SetParameters(1000,0,0.2,100,2.5,0.5);
		fSig2->ReleaseParameter(1);
		fSig2->ReleaseParameter(2);
//		fSig2->ReleaseParameter(3);
//		fSig2->ReleaseParameter(4);
//		fSig2->ReleaseParameter(5);
//		fSig2->SetParameter(3,500);
		if (fSig->GetParameter(4)<1e-5) { fSig->SetParameter(4,1); }
		h2->Fit(fSig2,"");
		//h2->Fit(fSig2,"");
		gPad->Draw();
		TPaveStats *stat2 = (TPaveStats*)h2->GetListOfFunctions()->FindObject("stats");
		stat2->SetY1NDC(0.3);    //set new y start position
		stat2->SetY2NDC(0.9);    //set new y end position
		stat2->SetX1NDC(0.1);    //set new x start position
		stat2->SetX2NDC(0.4);    //set new x end position
		stat2->SetTextColor(kGreen+2);
		/* TLegend *leg = new TLegend(0.45, 0.9, 0.6, 0.8); */
		/* leg->SetFillColor(0); */
		/* leg->AddEntry(h2, "Signal integral", "l"); */
		/* leg->Draw(); */
		c2->Print(Form("/lustrehome/bdanzi/offline_analysis/testbeam_analysis/Plots/histosTB_run_4.root/H%i_Integral_signal.png",i), "png");
		chn[nPresCh]=i;
		meanL[nPresCh]=fSig2->GetParameter(7/*4*/);//segnale
		sigmaL[nPresCh]=fSig2->GetParameter(8/*5*/);
		errMean[nPresCh]=fSig2->GetParError(7/*4*/);
		errSigma[nPresCh]=fSig2->GetParError(8/*5*/);
		/* cout << meanL[nPresCh] << " media landau "  <<endl; */
		/* cout << sigmaL[nPresCh] << " sigma landau "  <<endl; */
		if (myfile.is_open()) {
	    myfile << chn[nPresCh]<< "\t" << meanL[nPresCh] << "\t" <<errMean[nPresCh]<< "\t" << sigmaL[nPresCh] << "\t"<< errSigma[nPresCh]<<endl;
		  
		}
		++nPresCh;
	}
	if (myfile.is_open()) { myfile.close(); }
        
}
