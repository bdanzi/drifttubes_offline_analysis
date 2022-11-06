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


void integral(/*char file[300]="H-CH0_signal/hMinVNInR_ch0"*/){
   
TF1 *fit_function = new TF1("fit_function", "gaus");
fit_function->SetNpx(1000);

TFile *_file0 = new TFile("/lustrehome/bdanzi/offline_analysis/testbeam_analysis/histosTB_run_4.root","read");
TFile *_file1 = new TFile("/lustrehome/bdanzi/offline_analysis/testbeam_analysis/histosTB_run_4.root","read");

//_file0->ls();
 _file0->cd("H-Ch10_signal");                                         
 TH1F *h1=(TH1F*)gDirectory->Get("hIntegNInR_ch10");                 
h1->SetLineColor(kBlue);


// _file1->ls();
 _file1->cd("H-Ch10_signal");                                        
 TH1F *h2=(TH1F*)gDirectory->Get("hIntegNInR_ch10");                  
h2->SetLineColor(kGreen);
 
TCanvas *c1=new TCanvas("c1","c1");
c1->cd();
c1->SetLogy();
h1->GetXaxis()->SetTitle("Integral[V]");
h1->GetYaxis()->SetTitle("Entries");
//h1->GetYaxis()->SetRangeUser(0,45); 
gStyle->SetOptFit(1);
gStyle->SetStatX(1);
 
h1->Draw();
h1->Fit("fit_function","","",-2,2);
gPad->Update();
TPaveStats *stat = (TPaveStats*)h1->GetListOfFunctions()->FindObject("stats");
stat->SetY1NDC(0.7);    //set new y start position
stat->SetY2NDC(0.9);    //set new y end position
stat->SetX1NDC(0.1);    //set new y start position
stat->SetX2NDC(0.4);//set new y end position
stat->SetTextColor(kBlue);

h2->Draw("sames");
h2->Fit("fit_function","","",-2,2);
gPad->Draw();
TPaveStats *stat2 = (TPaveStats*)h2->GetListOfFunctions()->FindObject("stats");
stat2->SetY1NDC(0.5);    //set new y start position
stat2->SetY2NDC(0.7);    //set new y end position
stat2->SetX1NDC(0.1);    //set new x start position
stat2->SetX2NDC(0.4);    //set new x end position
stat2->SetTextColor(kGreen);

TLegend *leg = new TLegend(0.6, 0.9, 0.75, 0.8);
leg->SetFillColor(0);
//leg->SetHeader("Histogram");
leg->AddEntry(h1, "Filtered signal", "l");
leg->AddEntry(h2, "Filtered noise", "l"); 
leg->Draw();
 c1->Print("/lustrehome/bdanzi/offline_analysis/testbeam_analysis/Plots/histosTB_run_4.root/H10_IntegNInRangeLog.png", "png");                 
}

