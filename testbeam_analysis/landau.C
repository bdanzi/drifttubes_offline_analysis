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
#include <cstdlib>
using namespace std;

void landau(/*TString fname="/mnt/c/dati_root/histosOSC_run-00010.root",*/ TString sname = "/lustrehome/bdanzi/offline_analysis/testbeam_analysis/histosTB_run_4.root", TString fOutName=""){ //"/mnt/c/file_txt/ly2-3-4_90-10/1650.txt"

  ofstream myfile;
  if (!fOutName.IsNull()) {
    myfile.open (fOutName.Data());
  }
	gStyle->SetOptFit(1);
	gStyle->SetStatX(1);

	TFile *file2= new TFile(sname.Data(),"read"); //SIGNAL


	TF1 *fSig = new TF1("fSig", "landau");
	fSig->SetNpx(1000);
	fSig->SetParLimits(0,0,10000);
	fSig->SetParLimits(1,-2,2);
	fSig->SetParLimits(2,0,1);	
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
 
        c2->cd();
		//c2->cd()->SetLogy();
		nPresCh=0;
		for(int i =0; i<=11; ++i){

  		TH1F *h2=(TH1F*)file2->Get(Form("H-Ch%i_signal/hIntegNInRC2_ch%i",i,i));//SIGNAL
		if (h2==0x0) { continue; }
		h2->Rebin(4);
		h2->SetLineColor(kGreen+2);
		h2->Draw();	
	    fSig->SetParameters(500,0.5,0.04);
	 	h2->Fit(fSig,"L");
		h2->Fit(fSig,"");
		gPad->Draw();
		TPaveStats *stat2 = (TPaveStats*)h2->GetListOfFunctions()->FindObject("stats");
		stat2->SetY1NDC(0.3);    //set new y start position
		stat2->SetY2NDC(0.7);    //set new y end position
		stat2->SetX1NDC(0.1);    //set new x start position
		stat2->SetX2NDC(0.4);    //set new x end position
		stat2->SetTextColor(kGreen+2);
		/* TLegend *leg = new TLegend(0.45, 0.9, 0.6, 0.8); */
		/* leg->SetFillColor(0); */
		/* leg->AddEntry(h2, "Signal integral", "l"); */
		/* leg->Draw(); */

		c2->Print(Form("/lustrehome/bdanzi/offline_analysis/testbeam_analysis/Plots/histosTB_run_4.root/H%i_Integral_signal.png",i), "png");


		chn[nPresCh]=i;
		meanL[nPresCh]=fSig->GetParameter(1);//segnale
		sigmaL[nPresCh]=fSig->GetParameter(2);
		errMean[nPresCh]=fSig->GetParError(1);
		errSigma[nPresCh]=fSig->GetParError(2);
		/* cout << meanL[nPresCh] << " media landau "  <<endl; */
		/* cout << sigmaL[nPresCh] << " sigma landau "  <<endl; */
		if (myfile.is_open()) {
	    myfile << chn[nPresCh]<< "\t" << meanL[nPresCh] << "\t" <<errMean[nPresCh]<< "\t" << sigmaL[nPresCh] << "\t"<< errSigma[nPresCh]<<endl;
		  
		}
		++nPresCh;
	}
	if (myfile.is_open()) { myfile.close(); }

        
}
