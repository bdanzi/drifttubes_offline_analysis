/*
 * lognormal.C
 *
 *  Created on: 15 nov 2019
 *      Author: federica
 */
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

void lognormal(/*TString fname="/mnt/c/dati_root/histosOSC_run-00010.root",*/ TString sname = "/lustrehome/bdanzi/offline_analysis/testbeam_analysis/histosTB_run_4.root", TString fOutName=""){ //"/mnt/c/file_txt/ly2-3-4_90-10/1650.txt"
  ofstream myfile;
  if (!fOutName.IsNull()) {
    myfile.open (fOutName.Data());
  }

	gStyle->SetOptFit(1);
	gStyle->SetStatX(1);

	TFile *file2= new TFile(sname.Data(),"read");


	TF1 *fSig = new TF1("fSig", "[0]*ROOT::Math::lognormal_pdf(x,[1],[2])",-10,100);
	fSig->SetNpx(1000);
	TCanvas *c2 = new TCanvas ("c2","Integral for signal");
	Int_t nPresCh=0;
	Int_t nPresCH_1=0;
	Float_t chn[12]={0,0,0,0,0,0,0,0,0,0,0,0};
	Float_t meanLog[12]={0,0,0,0,0,0,0,0,0,0,0,0};
	Float_t sigmaLog[12]={0,0,0,0,0,0,0,0,0,0,0,0};
	Float_t errMeanLog[12]={0,0,0,0,0,0,0,0,0,0,0,0};
	Float_t errSigmaLog[12]={0,0,0,0,0,0,0,0,0,0,0,0};

       c2->cd();
       //c2->cd()->SetLogy();
       nPresCh=0;
       double p[3];
       for(int i =0; i<=11; ++i){
    	   TH1F *h2=(TH1F*)file2->Get(Form("H-Ch%i_signal/hIntegNInRC3_ch%i",i,i));//SIGNAL without norm
    	   //TH1F *h2=(TH1F*)file2->Get(Form("H-Ch%i_signal/hIntegNInRC4_ch%i",i,i));//SIGNAL with norm loss only
    	   if (h2==0x0) { continue; }
    	   h2->Rebin(4);
    	   h2->SetLineColor(kGreen+2);
    	   h2->Draw();

    	   h2->GetXaxis()->SetRangeUser(0.5,20);
    	   p[0] = h2->GetEntries()*h2->GetXaxis()->GetBinWidth(1);
		// find median of histogram
    	   double prob[] = {0.5};
    	   double q[1];
    	   h2->GetQuantiles(1,q,prob);
    	   double median = q[0];
		// find mode of histogram

    	   double  mode = h2->GetBinCenter( h2->GetMaximumBin());
    	   cout << "histogram mode is " << mode  << " median is " << median << std::endl;
    	   h2->GetXaxis()->UnZoom();
		// m is log(median)
    	   p[1] = std::log(median);
		// s2 is  log(median) - log(mode)
    	   if(mode==0){ p[2]=0.5;}
    	   else {
    		   double tmplog=std::log(median/mode);
    		   if(tmplog<0){p[2]=0.5;}
    		   else{
    			   p[2] = std::sqrt( std::log(median/mode) );
    		   }
    	   }
    	   cout << "p0"<<p[0]<<" p1 "<< p[1] << " p2 "<< p[2]<<endl;


    	   fSig->SetParameters(p);

    	   fSig->SetParName(0,"A");
    	   fSig->SetParName(1,"m");
    	   fSig->SetParName(2,"s");

    	   h2->Fit(fSig,"G","",mode*0.8,20.0);
    	   h2->Fit(fSig,"G","",mode*0.8,20.0);
    	   //h2->Fit(fSig,"L","",mode*0.55,15.0);
    	   //h2->Fit(fSig,"","");
    	   gPad->Draw();
    	   TPaveStats *stat2 = (TPaveStats*)h2->GetListOfFunctions()->FindObject("stats");
    	   stat2->SetY1NDC(0.4);    //set new y start position
    	   stat2->SetY2NDC(0.7);    //set new y end position
    	   stat2->SetX1NDC(0.5);    //set new x start position
    	   stat2->SetX2NDC(0.9);    //set new x end position
    	   stat2->SetTextColor(kGreen+2);
    	   c2->Print(Form("/lustrehome/bdanzi/offline_analysis/testbeam_analysis/Plots/histosTB_run_4.root/H%i_Integral_control.png",i), "png");
		  //c2->Print(Form("/home/federica/eclipse-workspace/test_beam/Plot/Integral_control/norm_loss/H%i_Integral_control.png",i), "png");
    	   chn[nPresCh]=i;
    	   meanLog[nPresCh]=fSig->GetParameter(1);//segnale
    	   sigmaLog[nPresCh]=fSig->GetParameter(2);
    	   errMeanLog[nPresCh]=fSig->GetParError(1);
    	   errSigmaLog[nPresCh]=fSig->GetParError(2);
    	   cout <<  " media log " <<meanLog[nPresCh] <<endl;
    	   cout << " sigma log " << sigmaLog[nPresCh] <<endl;
    	   if (myfile.is_open()) {
    		   myfile << chn[nPresCh]<< "\t" << meanLog[nPresCh] << "\t" <<errMeanLog[nPresCh]<< "\t" << sigmaLog[nPresCh] << "\t"<< errSigmaLog[nPresCh]<<endl;
    	   }
    	   ++nPresCh;
	}
	if (myfile.is_open()) { myfile.close(); }

}



