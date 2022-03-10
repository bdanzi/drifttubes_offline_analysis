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
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include "TMath.h"
#include "TH1F.h"
#include "iostream"
#include "TFile.h"
#include "TLegend.h"

using namespace std;

void ReadSglChannel_test(){
  gStyle->SetOptFit(1011);

  int channel=0;
  int ev=0;
  struct stat st = {0};
  char filename[100];
  TString fname("");
  bool isInteractive = false;
  bool isdoubleCanvas = false;
  Float_t cos_alpha = TMath::Cos(60*TMath::DegToRad());
  Float_t expected_electrons =0.;
  Float_t expected_cluster =0.;
  Float_t cluster_per_cm_mip = 12.;
  Float_t drift_size =0.;
  Float_t relativistic_rise = 1.3;
  Float_t cluster_population = 1.6;
	 
  if(isInteractive){
    
    printf("Enter *.root name:");
    scanf("%s",fname.Data());
    printf("Enter the number of channels to be shown:");
    scanf("%d",&channel);
    printf("Enter the number of events (waveForms) to be shown:");
    scanf("%d",&ev);
    if (stat(Form("/lustrehome/bdanzi/offline_analysis/testbeam_analysis/Plots/%s/",fname.Data()), &st) == -1) {
      mkdir(Form("/lustrehome/bdanzi/offline_analysis/testbeam_analysis/Plots/%s/",fname.Data()), 0700);
    }
    if (stat(Form("/lustrehome/bdanzi/offline_analysis/testbeam_analysis/Plots/%s/Waves/",fname.Data()), &st) == -1) {
      mkdir(Form("/lustrehome/bdanzi/offline_analysis/testbeam_analysis/Plots/%s/Waves/",fname.Data()), 0700);
    }
  }
  else{
    char inputfile[100] ="plots.txt";
    FILE *fp = fopen(inputfile, "r");
    if (fp == NULL)
      {
        printf("Error: could not open file %s", inputfile);
       	return 1;
      }
    
    while ((fscanf(fp, "%s", fname.Data() ) ) != -1){
      
      printf("Opening the file %s\n",fname.Data());
      fscanf(fp,"%d",&channel); //number of channels for which I show plots (4-14)
      fscanf(fp,"%d",&ev); //number of events for which I show plots (0-...)
      
      
      
      if (stat(Form("/lustrehome/bdanzi/offline_analysis/testbeam_analysis/Plots/%s/",fname.Data()), &st) == -1) {
	mkdir(Form("/lustrehome/bdanzi/offline_analysis/testbeam_analysis/Plots/%s/",fname.Data()), 0700);
      }
      if (stat(Form("/lustrehome/bdanzi/offline_analysis/testbeam_analysis/Plots/%s/Waves/",fname.Data()), &st) == -1) {
	mkdir(Form("/lustrehome/bdanzi/offline_analysis/testbeam_analysis/Plots/%s/Waves/",fname.Data()), 0700);
      }
      
      TFile *file = new TFile(fname.Data(),"read");
      TCanvas *max_1cm=new TCanvas("max_1","max",3500,1500);
      TCanvas *max_2cm=new TCanvas("max_2","max",3500,1500);
      TCanvas *min=new TCanvas("min","min",3500,1500);
      TCanvas *bsl_1cm=new TCanvas("bsl_1","bsl",3500,1500);
      TCanvas *bsl_2cm=new TCanvas("bsl_2","bsl",3500,1500);
      TCanvas *npeaks_1cm= new TCanvas("npeaks_1","npeaks",3500,1500);
      TCanvas *npeaks_2cm= new TCanvas("npeaks_2","npeaks",3500,1500);
	  TCanvas *npeaks_clustser_1cm= new TCanvas("npeaks_1_cluster","npeaks",3500,1500);
      TCanvas *npeaks_clustser_2cm= new TCanvas("npeaks_2_cluster","npeaks",3500,1500);
      TCanvas *hpeaks_1cm= new TCanvas("hpeaks_1","hpeaks",3500,1500);
      TCanvas *hpeaks_2cm= new TCanvas("hpeaks_2","hpeaks",3500,1500);
      TCanvas *hnpeaks_1cm= new TCanvas("hnpeaks_1","hnpeaks",3500,1500);
      TCanvas *hnpeaks_2cm= new TCanvas("hnpeaks_2","hnpeaks",3500,1500);
      TCanvas *tpeaks_1cm= new TCanvas("tpeaks_1","tpeaks",3500,1500);
      TCanvas *tpeaks_2cm= new TCanvas("tpeaks_2","tpeaks",3500,1500);
      TCanvas *tfpeaks_1cm= new TCanvas("tfpeaks_1","tfpeaks",3500,1500);
      TCanvas *tfpeaks_2cm= new TCanvas("tfpeaks_2","tfpeaks",3500,1500);
      TCanvas *tlpeaks_1cm= new TCanvas("tlpeaks_1","tlpeaks",3500,1500);
      TCanvas *tlpeaks_2cm= new TCanvas("tlpeaks_2","tlpeaks",3500,1500);
      TCanvas *integ_1cm= new TCanvas("integ_1","integ",3500,1500);
      TCanvas *integ_2cm= new TCanvas("integ_2","integ",3500,1500);
      TCanvas *rms_1cm= new TCanvas("rms_1","rms",3500,1500);
      TCanvas *rms_2cm= new TCanvas("rms_2","rms",3500,1500);
	  TCanvas *cluster_population_canvas_1cm= new TCanvas("cluster_population_1cm","cluster",3500,1500);
	  TCanvas *cluster_population_canvas_2cm= new TCanvas("cluster_population_2cm","cluster",3500,1500);
	  
      cluster_population_canvas_1cm->Divide(2,3);
	  cluster_population_canvas_2cm->Divide(2,2);
      npeaks_clustser_1cm->Divide(2,3);
	  npeaks_clustser_2cm->Divide(2,2);
      //integ->Divide(4,2);
      //max->Divide(4,2);
      min->Divide(4,1);
      bsl_1cm->Divide(2,3);
      bsl_2cm->Divide(2,2);
      integ_1cm->Divide(2,3);
      integ_2cm->Divide(2,2);
      rms_1cm->Divide(2,3);
      rms_2cm->Divide(2,2);
      max_1cm->Divide(2,3);
      max_2cm->Divide(2,2);
      npeaks_1cm->Divide(2,3);
      npeaks_2cm->Divide(2,2);
      hpeaks_1cm->Divide(2,3);
      hpeaks_2cm->Divide(2,2);
      hnpeaks_1cm->Divide(2,3);
      hnpeaks_2cm->Divide(2,2);
      tpeaks_1cm->Divide(2,3);
      tpeaks_2cm->Divide(2,2);
      tfpeaks_1cm->Divide(2,3);
      tfpeaks_2cm->Divide(2,2);
      tlpeaks_1cm->Divide(2,3);
      tlpeaks_2cm->Divide(2,2);
      
      
      if(isdoubleCanvas){
	for(int i =0; i<ev; ++i){
	  TGraph *h1=(TGraph*)file->Get(Form("signal_Afterflt/CvSignal_1_ev%i",i));
	  if (h1==0x0) { continue; }
	  h1->SaveAs(Form("/lustrehome/bdanzi/offline_analysis/testbeam_analysis/Plots/%s/Waves/waves_ev%i.pdf",fname.Data(),i));
	}
      }
      else if(!isdoubleCanvas){
	
	for(int i = 0; i<ev; ++i){
	  for(int j = 4; j<=channel; ++j){
	    TGraph *h1=(TGraph*)file->Get(Form("signal_Afterflt/CvSignal_1_Ch%i_ev%i",j,i));
	    if (h1==0x0) { continue; }
	    //h1->GetYaxis()->SetRangeUser(-0.1,0.7);
	    h1->SaveAs(Form("/lustrehome/bdanzi/offline_analysis/testbeam_analysis/Plots/%s/Waves/waves_ev%i_Ch%i.pdf",fname.Data(),i,j));
	  }
	}
      }
      
      for(int i = 4; i<=9 ;++i){ //looping over all channels
	//for(int i =channel; i<=channel; ++i){ //looping over one channel
	drift_size = 0.8;
	//δ cluster/cm (M.I.P.) *drift tube size [cm] *1.3 (relativisticrise)*1.6 electrons/cluster*1/cos(α)
	expected_electrons = cluster_per_cm_mip * drift_size*relativistic_rise * cluster_population * 1/cos_alpha;
	expected_cluster = cluster_per_cm_mip * drift_size*relativistic_rise * 1/cos_alpha;
	
	
	//TH1F *h2=(TH1F*)file->Get(Form("H-Ch%i_signal/hMaxVN_ch%i",i,i));
	//if (h2==0x0) { continue; }
	//max->cd(1);
	//h2->Draw();
	
	TH1F *h17=(TH1F*)file->Get(Form("H-Ch%i_signal/hMaxVInR_ch%i",i,i));
	if (h17==0x0) { continue; }
	max_1cm->cd(i-3);
	//max->cd(2);
	gPad->SetLogy(1);
	gPad->SetLogx(0);
	h17->Draw( "same");
	
	//TH1F *h18=(TH1F*)file->Get(Form("H-Ch%i_signal/hMaxV_ch%i",i,i));
	//if (h18==0x0) { continue; }
	//max->cd(3);
	//h18->Draw( );
	
	//TH1F *h19=(TH1F*)file->Get(Form("H-Ch%i_signal/hMaxVNInR_ch%i",i,i));
	//if (h19==0x0) { continue; }
	//max->cd(4);
	//h19->Draw( );
	
	//TH1F *h30=(TH1F*)file->Get(Form("H-Ch%i_signal/hMaxVoriginalW_ch%i",i,i));
	//if (h30==0x0) { continue; }
	//max->cd(5);
	//h30->Draw( );
	
	//TH1F *h31=(TH1F*)file->Get(Form("H-Ch%i_signal/hMaxVNoriginalW_ch%i",i,i));
	//if (h31==0x0) { continue; }
	//max->cd(6);
	//h31->Draw( );
	
	//TH1F *h32=(TH1F*)file->Get(Form("H-Ch%i_signal/hMaxVInRoriginalW_ch%i",i,i));
	//if (h32==0x0) { continue; }
	//max->cd(7);
	//h32->Draw( );
	
	//TH1F *h33=(TH1F*)file->Get(Form("H-Ch%i_signal/hMaxVNInRoriginalW_ch%i",i,i));
	//if (h33==0x0) { continue; }
	//max->cd(8);
	//h33->Draw( );
	
	
	TH1F *h3=(TH1F*)file->Get(Form("H-Ch%i_signal/hBsl_ch%i",i,i));
	if (h3==0x0) { continue; }
	bsl_1cm->cd(i-3);
	gPad->SetLogy(0);
	gPad->SetLogx(0);
	h3->GetXaxis()->SetRangeUser(-0.5,-0.4);
	h3->Draw("same");
	
	
	TH1F *h4=(TH1F*)file->Get(Form("H-Ch%i_signal/hNPeaks_ch%i",i,i));
	if (h4==0x0) { continue; }
	npeaks_1cm->cd(i-3);
	h4->GetXaxis()->SetRangeUser(0.,90.);
	h4->Fit("gaus");
	gPad->SetLogy(1);
	h4->Draw("same");
	TPaveText *pt_1cm = new TPaveText(0.72,0.2,0.8,0.35,"NDC");
	gPad->SetLogy(1);
	pt_1cm->SetTextSize(0.05);
	pt_1cm->SetTextColor(kRed);
	pt_1cm->SetFillColor(0);
	pt_1cm->SetTextAlign(12);
	pt_1cm->AddText(Form("Expected elecrons: %.1f",expected_electrons));
	pt_1cm->Draw("same");

	TH1F *h20=(TH1F*)file->Get(Form("H-Ch%i_signal/hNPeaks_clust_ch%i",i,i));
	if (h20==0x0) { continue; }
	npeaks_clustser_1cm->cd(i-3);
	h20->GetXaxis()->SetRangeUser(0.,90.);
	h20->Fit("gaus");
	gPad->SetLogy(1);
	h20->Draw("same");
	TPaveText *pt_1cm_cluster = new TPaveText(0.72,0.2,0.8,0.35,"NDC");
	gPad->SetLogy(1);
	pt_1cm_cluster->SetTextSize(0.04);
	pt_1cm_cluster->SetTextColor(kRed);
	pt_1cm_cluster->SetTextAlign(12);
	pt_1cm_cluster->SetFillColor(0);
	pt_1cm_cluster->AddText(Form("Expected Clusters: %.1f",expected_cluster));
	gPad->SetLogy(1);
	pt_1cm_cluster->Draw("same");

	TH1F *h31=(TH1F*)file->Get(Form("H-Ch%i_signal/hNElectrons_per_cluster_ch%i",i,i));
	if (h31==0x0) { continue; }
	cluster_population_canvas_1cm->cd(i-3);
	//h31->GetXaxis()->SetRangeUser(0.,5.);
	//h31->Fit("expo");
	h31->Draw("same");
	TPaveText *cluster_population_1cm = new TPaveText(0.72,0.2,0.8,0.35,"NDC");
	gPad->SetLogy(1);
	cluster_population_1cm->SetTextSize(0.04);
	cluster_population_1cm->SetTextColor(kRed);
	cluster_population_1cm->SetTextAlign(12);
	cluster_population_1cm->SetFillColor(0);
	cluster_population_1cm->AddText(Form("Expected Electrons per Cluster: %.1f",cluster_population));
	gPad->SetLogy(1);
	cluster_population_1cm->Draw("same");
	
	TH1F *h5=(TH1F*)file->Get(Form("H-Ch%i_signal/hHPeaks_ch%i",i,i));
	if (h5==0x0) { continue; }
	hpeaks_1cm->cd(i-3);
	h5->Draw("same");
	
	TH1F *h6=(TH1F*)file->Get(Form("H-Ch%i_signal/hHNPeaks_ch%i",i,i));
	if (h6==0x0) { continue; }
	hnpeaks_1cm->cd(i-3);
	h6->GetYaxis()->SetTitle("Height of Peaks found [V]");
	h6->GetXaxis()->SetRangeUser(0.,90.);
	h6->Draw("colz");
	gPad->Update();
	TPaveStats *st_1cm = (TPaveStats*)h6->FindObject("stats");
	st_1cm->SetX1NDC(0.75); //new x start position
	st_1cm->SetX2NDC(0.85); //new x end position
	st_1cm->SetY1NDC(0.65); //new x start position
	st_1cm->SetY2NDC(0.85); //new x end position
	
	TH1F *h7=(TH1F*)file->Get(Form("H-Ch%i_signal/hTPeaks_ch%i",i,i));
	if (h7==0x0) { continue; }
	tpeaks_1cm->cd(i-3);
	h7->Draw("same");
	
	TH1F *h8=(TH1F*)file->Get(Form("H-Ch%i_signal/hTFstPeaks_ch%i",i,i));
	if (h8==0x0) { continue; }
	tfpeaks_1cm->cd(i-3);
	h8->Draw("same");
	
	TH1F *h9=(TH1F*)file->Get(Form("H-Ch%i_signal/hNPeaks_1_ch%i",i,i));
	if (h9==0x0) { continue; }
	tlpeaks_1cm->cd(i-3);
	h9->Draw("same");
	
	//TH1F *h10=(TH1F*)file->Get(Form("H-Ch%i_signal/hInteg_ch%i",i,i));
	//if (h10==0x0) { continue; }
	//integ->cd(1);
	//h10->Draw();
	
	//TH1F *h11=(TH1F*)file->Get(Form("H-Ch%i_signal/hIntegN_ch%i",i,i));
	//if (h11==0x0) { continue; }
	//integ->cd(2);
	//h11->Draw( );
	
	TH1F *h12=(TH1F*)file->Get(Form("H-Ch%i_signal/hIntegInR_ch%i",i,i));
	if (h12==0x0) { continue; }
	integ_1cm->cd(i-3);
	//integ->cd(3);
	h12->Draw("same");
	
	//TH1F *h13=(TH1F*)file->Get(Form("H-Ch%i_signal/hIntegNInR_ch%i",i,i));
	//if (h13==0x0) { continue; }
	//integ->cd(4);
	//h13->Draw( );
	
	//TH1F *h14=(TH1F*)file->Get(Form("H-Ch%i_signal/hIntegNInRC1_ch%i",i,i));
	//if (h14==0x0) { continue; }
	//integ->cd(5);
	//h14->Draw( );
	
	//TH1F *h15=(TH1F*)file->Get(Form("H-Ch%i_signal/hIntegNInRC2_ch%i",i,i));
	//if (h15==0x0) { continue; }
	//integ->cd(6);
	//h15->Draw( );
	
	//TH1F *h25=(TH1F*)file->Get(Form("H-Ch%i_signal/hIntegNInRoriginalW_ch%i",i,i));
	//if (h25==0x0) { continue; }
	//integ->cd(7);
	//h25->Draw( );
	
	//TH1F *h26=(TH1F*)file->Get(Form("H-Ch%i_signal/hIntegNInRoriginalW_ch%i",i,i));
	//if (h26==0x0) { continue; }
	//integ->cd(8);
	//h26->Draw( );
	
	TH1F *h16=(TH1F*)file->Get(Form("H-Ch%i_signal/hRms_ch%i",i,i));
	if (h16==0x0) { continue; }
	rms_1cm->cd(i-3);
	h16->Draw("same");
	
	/*
	  TH1F *h20=(TH1F*)file->Get(Form("H-Ch%i_signal/hMinV_ch%i",i,i));
	  if (h20==0x0) { continue; }
	  min->cd(1);
	  h20->Draw();
	  
	  TH1F *h21=(TH1F*)file->Get(Form("H-Ch%i_signal/hMinVN_ch%i",i,i));
	  min->cd(2);
	  h21->Draw( );
	  
	  TH1F *h22=(TH1F*)file->Get(Form("H-Ch%i_signal/hMinVInR_ch%i",i,i));
	  if (h22==0x0) { continue; }
	  min->cd(3);
	  h22->Draw( );
	  
	  TH1F *h23=(TH1F*)file->Get(Form("H-Ch%i_signal/hMinVNInR_ch%i",i,i));
	  if (h23==0x0) { continue; }
	  min->cd(4);
	  h23->Draw( );
	*/
	/*TH1F *h27=(TH1F*)file->Get(Form("H-Ch%i_signal/hMinVoriginalW_ch%i",i,i));
	  if (h27==0x0) { continue; }
	  min->cd(5);
	  h27->Draw( );
	  
	  TH1F *h28=(TH1F*)file->Get(Form("H-Ch%i_signal/hMinVNoriginalW_ch%i",i,i));
	  if (h28==0x0) { continue; }
	  min->cd(6);
	  h28->Draw( );
	  
	  TH1F *h29=(TH1F*)file->Get(Form("H-Ch%i_signal/hMinVInRoriginalW_ch%i",i,i));
	  if (h29==0x0) { continue; }
	  min->cd(7);
	  h29->Draw( );*/
	
	
	 
	
      }

	bool savePlots = true;
	
	if(savePlots){
	  
	  npeaks_1cm->SaveAs(Form("/lustrehome/bdanzi/offline_analysis/testbeam_analysis/Plots/%s/npeaks_1cm.pdf",fname.Data()));
	  tpeaks_1cm->SaveAs(Form("/lustrehome/bdanzi/offline_analysis/testbeam_analysis/Plots/%s/tpeaks_1cm.pdf",fname.Data()));
	  tfpeaks_1cm->SaveAs(Form("/lustrehome/bdanzi/offline_analysis/testbeam_analysis/Plots/%s/tfpeaks_1cm.pdf",fname.Data()));
	  tlpeaks_1cm->SaveAs(Form("/lustrehome/bdanzi/offline_analysis/testbeam_analysis/Plots/%s/tlpeaks_1cm.pdf",fname.Data()));
	  hnpeaks_1cm->SaveAs(Form("/lustrehome/bdanzi/offline_analysis/testbeam_analysis/Plots/%s/hnpeaks_1cm.pdf",fname.Data()));
	  hpeaks_1cm->SaveAs(Form("/lustrehome/bdanzi/offline_analysis/testbeam_analysis/Plots/%s/hpeaks_1cm.pdf",fname.Data()));
	  bsl_1cm->SaveAs(Form("/lustrehome/bdanzi/offline_analysis/testbeam_analysis/Plots/%s/bsl_1cm.pdf",fname.Data()));
	  max_1cm->SaveAs(Form("/lustrehome/bdanzi/offline_analysis/testbeam_analysis/Plots/%s/max_1cm.pdf",fname.Data()));
	  integ_1cm->SaveAs(Form("/lustrehome/bdanzi/offline_analysis/testbeam_analysis/Plots/%s/integ_1cm.pdf",fname.Data()));
	  rms_1cm->SaveAs(Form("/lustrehome/bdanzi/offline_analysis/testbeam_analysis/Plots/%s/rms_1cm.pdf",fname.Data()));
	  //min->SaveAs(Form("/lustrehome/bdanzi/offline_analysis/testbeam_analysis/Plots/%s/min_ch%i.pdf",fname.Data(),i));
	  npeaks_clustser_1cm->SaveAs(Form("/lustrehome/bdanzi/offline_analysis/testbeam_analysis/Plots/%s/hnpeaks_cluster_1cm.pdf",fname.Data()));
	  cluster_population_canvas_1cm->SaveAs(Form("/lustrehome/bdanzi/offline_analysis/testbeam_analysis/Plots/%s/cluster_population_1cm.pdf",fname.Data()));
	}
      
      for(int i = 10; i<=channel; ++i){ //looping over all channels
	//for(int i =channel; i<=channel; ++i){ //looping over one channel
	drift_size = 1.8;
	expected_electrons = cluster_per_cm_mip * drift_size*relativistic_rise * cluster_population * 1/cos_alpha;
	expected_cluster = cluster_per_cm_mip * drift_size*relativistic_rise * 1/cos_alpha;
	
	
	
	//TH1F *h2=(TH1F*)file->Get(Form("H-Ch%i_signal/hMaxVN_ch%i",i,i));
	//if (h2==0x0) { continue; }
	//max->cd(1);
	//h2->Draw();
	
	TH1F *h17=(TH1F*)file->Get(Form("H-Ch%i_signal/hMaxVInR_ch%i",i,i));
	if (h17==0x0) { continue; }
	max_2cm->cd(i-9);
	//max->cd(2);
	gPad->SetLogy(1);
	gPad->SetLogx(0);
	h17->Draw("same" );
	
	//TH1F *h18=(TH1F*)file->Get(Form("H-Ch%i_signal/hMaxV_ch%i",i,i));
	//if (h18==0x0) { continue; }
	//max->cd(3);
	//h18->Draw( );
	
	//TH1F *h19=(TH1F*)file->Get(Form("H-Ch%i_signal/hMaxVNInR_ch%i",i,i));
	//if (h19==0x0) { continue; }
	//max->cd(4);
	//h19->Draw( );
	
	//TH1F *h30=(TH1F*)file->Get(Form("H-Ch%i_signal/hMaxVoriginalW_ch%i",i,i));
	//if (h30==0x0) { continue; }
	//max->cd(5);
	//h30->Draw( );
	
	//TH1F *h31=(TH1F*)file->Get(Form("H-Ch%i_signal/hMaxVNoriginalW_ch%i",i,i));
	//if (h31==0x0) { continue; }
	//max->cd(6);
	//h31->Draw( );
	
	//TH1F *h32=(TH1F*)file->Get(Form("H-Ch%i_signal/hMaxVInRoriginalW_ch%i",i,i));
	//if (h32==0x0) { continue; }
	//max->cd(7);
	//h32->Draw( );
	
	//TH1F *h33=(TH1F*)file->Get(Form("H-Ch%i_signal/hMaxVNInRoriginalW_ch%i",i,i));
	//if (h33==0x0) { continue; }
	//max->cd(8);
	//h33->Draw( );
	
	
	TH1F *h3=(TH1F*)file->Get(Form("H-Ch%i_signal/hBsl_ch%i",i,i));
	if (h3==0x0) { continue; }
	bsl_2cm->cd(i-9);
	gPad->SetLogy(0);
	gPad->SetLogx(0);
	h3->GetXaxis()->SetRangeUser(-0.5,-0.4);
	h3->Draw("same");
	
	
	
	TH1F *h4=(TH1F*)file->Get(Form("H-Ch%i_signal/hNPeaks_ch%i",i,i));
	if (h4==0x0) { continue; }
	npeaks_2cm->cd(i-9);
	h4->GetXaxis()->SetRangeUser(0.,200.);
	h4->Fit("gaus");
	gPad->SetLogy(1);
	h4->Draw("same");
	TPaveText *pt_2cm = new TPaveText(0.72,0.2,0.8,0.35,"NDC");
	pt_2cm->SetTextSize(0.04);
	pt_2cm->SetTextColor(kRed);
	pt_2cm->SetFillColor(0);
	pt_2cm->SetTextAlign(12);
	pt_2cm->AddText(Form("Expected elecrons: %.1f",expected_electrons));
	pt_2cm->Draw("same");
	
	TH1F *h5=(TH1F*)file->Get(Form("H-Ch%i_signal/hHPeaks_ch%i",i,i));
	if (h5==0x0) { continue; }
	hpeaks_2cm->cd(i-9);
	h5->Draw("same");
	
	TH1F *h6=(TH1F*)file->Get(Form("H-Ch%i_signal/hHNPeaks_ch%i",i,i));
	if (h6==0x0) { continue; }
	hnpeaks_2cm->cd(i-9);
	h6->GetYaxis()->SetTitle("Height of Peaks found [V]");
	h6->GetXaxis()->SetRangeUser(0.,200.);
	h6->Draw( "colz");
	gPad->Update();
	TPaveStats *st_2cm = (TPaveStats*)h6->FindObject("stats");
	st_2cm->SetX1NDC(0.75); //new x start position
	st_2cm->SetX2NDC(0.85); //new x end position
	st_2cm->SetY1NDC(0.65); //new x start position
	st_2cm->SetY2NDC(0.85); //new x end position
	
	
	TH1F *h7=(TH1F*)file->Get(Form("H-Ch%i_signal/hTPeaks_ch%i",i,i));
	if (h7==0x0) { continue; }
	tpeaks_2cm->cd(i-9);
	h7->Draw( "same");
	
	TH1F *h8=(TH1F*)file->Get(Form("H-Ch%i_signal/hTFstPeaks_ch%i",i,i));
	if (h8==0x0) { continue; }
	tfpeaks_2cm->cd(i-9);
	h8->Draw("same" );
	
	TH1F *h9=(TH1F*)file->Get(Form("H-Ch%i_signal/hNPeaks_1_ch%i",i,i));
	if (h9==0x0) { continue; }
	tlpeaks_2cm->cd(i-9);
	h9->Draw( "same");
	//TH1F *h10=(TH1F*)file->Get(Form("H-Ch%i_signal/hInteg_ch%i",i,i));
	//if (h10==0x0) { continue; }
	//integ->cd(1);
	//h10->Draw();
	
	//TH1F *h11=(TH1F*)file->Get(Form("H-Ch%i_signal/hIntegN_ch%i",i,i));
	//if (h11==0x0) { continue; }
	//integ->cd(2);
	//h11->Draw( );
	
	TH1F *h12=(TH1F*)file->Get(Form("H-Ch%i_signal/hIntegInR_ch%i",i,i));
	if (h12==0x0) { continue; }
	integ_2cm->cd(i-9);
	//integ->cd(3);
	h12->Draw("same" );
	
	//TH1F *h13=(TH1F*)file->Get(Form("H-Ch%i_signal/hIntegNInR_ch%i",i,i));
	//if (h13==0x0) { continue; }
	//integ->cd(4);
	//h13->Draw( );
	
	//TH1F *h14=(TH1F*)file->Get(Form("H-Ch%i_signal/hIntegNInRC1_ch%i",i,i));
	//if (h14==0x0) { continue; }
	//integ->cd(5);
	//h14->Draw( );
	
	//TH1F *h15=(TH1F*)file->Get(Form("H-Ch%i_signal/hIntegNInRC2_ch%i",i,i));
	//if (h15==0x0) { continue; }
	//integ->cd(6);
	//h15->Draw( );
	
	//TH1F *h25=(TH1F*)file->Get(Form("H-Ch%i_signal/hIntegNInRoriginalW_ch%i",i,i));
	//if (h25==0x0) { continue; }
	//integ->cd(7);
	//h25->Draw( );
	
	//TH1F *h26=(TH1F*)file->Get(Form("H-Ch%i_signal/hIntegNInRoriginalW_ch%i",i,i));
	//if (h26==0x0) { continue; }
	//integ->cd(8);
	//h26->Draw( );
	
	TH1F *h16=(TH1F*)file->Get(Form("H-Ch%i_signal/hRms_ch%i",i,i));
	if (h16==0x0) { continue; }
	rms_2cm->cd(i-9);
	h16->Draw("same");
			
	/*
	  TH1F *h20=(TH1F*)file->Get(Form("H-Ch%i_signal/hMinV_ch%i",i,i));
	  if (h20==0x0) { continue; }
	  min->cd(1);
	  h20->Draw();
	  
	  TH1F *h21=(TH1F*)file->Get(Form("H-Ch%i_signal/hMinVN_ch%i",i,i));
	  min->cd(2);
	  h21->Draw( );
	  
	  TH1F *h22=(TH1F*)file->Get(Form("H-Ch%i_signal/hMinVInR_ch%i",i,i));
	  if (h22==0x0) { continue; }
	  min->cd(3);
	  h22->Draw( );
	  
	  TH1F *h23=(TH1F*)file->Get(Form("H-Ch%i_signal/hMinVNInR_ch%i",i,i));
	  if (h23==0x0) { continue; }
	  min->cd(4);
	  h23->Draw( );
	*/
	/*TH1F *h27=(TH1F*)file->Get(Form("H-Ch%i_signal/hMinVoriginalW_ch%i",i,i));
	  if (h27==0x0) { continue; }
	  min->cd(5);
	  h27->Draw( );
	  
	  TH1F *h28=(TH1F*)file->Get(Form("H-Ch%i_signal/hMinVNoriginalW_ch%i",i,i));
	  if (h28==0x0) { continue; }
	  min->cd(6);
	  h28->Draw( );
	  
	  TH1F *h29=(TH1F*)file->Get(Form("H-Ch%i_signal/hMinVInRoriginalW_ch%i",i,i));
	  if (h29==0x0) { continue; }
	  min->cd(7);
	  h29->Draw( );*/
		TH1F *h20=(TH1F*)file->Get(Form("H-Ch%i_signal/hNPeaks_clust_ch%i",i,i));
		if (h20==0x0) { continue; }
		npeaks_clustser_2cm->cd(i-9);
		h20->GetXaxis()->SetRangeUser(0.,90.);
		h20->Fit("gaus");
		gPad->SetLogy(1);
		h20->Draw("same");
		TPaveText *pt_2cm_cluster = new TPaveText(0.72,0.2,0.8,0.35,"NDC");
		pt_2cm_cluster->SetTextSize(0.04);
		pt_2cm_cluster->SetTextColor(kRed);
		pt_2cm_cluster->SetFillColor(0);
		pt_2cm_cluster->SetTextAlign(12);
		pt_2cm_cluster->AddText(Form("Expected Clusters: %.1f",expected_cluster));
		pt_2cm_cluster->Draw("same");

		TH1F *h31=(TH1F*)file->Get(Form("H-Ch%i_signal/hNElectrons_per_cluster_ch%i",i,i));
		if (h31==0x0) { continue; }
		cluster_population_canvas_2cm->cd(i-9);
		//h31->GetXaxis()->SetRangeUser(0.,5.);
		//h31->Fit("expo");
		h31->Draw("same");
		TPaveText *cluster_population_2cm = new TPaveText(0.72,0.2,0.8,0.35,"NDC");
		gPad->SetLogy(1);
		cluster_population_2cm->SetTextSize(0.04);
		cluster_population_2cm->SetTextColor(kRed);
		cluster_population_2cm->SetTextAlign(12);
		cluster_population_2cm->SetFillColor(0);
		cluster_population_2cm->AddText(Form("Expected Electrons per Cluster: %.1f",cluster_population));
		gPad->SetLogy(1);
		cluster_population_2cm->Draw("same");
	

	
	
      }

		
	  if(savePlots){
	  
	  npeaks_2cm->SaveAs(Form("/lustrehome/bdanzi/offline_analysis/testbeam_analysis/Plots/%s/npeaks_2cm.pdf",fname.Data()));
	  tpeaks_2cm->SaveAs(Form("/lustrehome/bdanzi/offline_analysis/testbeam_analysis/Plots/%s/tpeaks_2cm.pdf",fname.Data()));
	  tfpeaks_2cm->SaveAs(Form("/lustrehome/bdanzi/offline_analysis/testbeam_analysis/Plots/%s/tfpeaks_2cm.pdf",fname.Data()));
	  tlpeaks_2cm->SaveAs(Form("/lustrehome/bdanzi/offline_analysis/testbeam_analysis/Plots/%s/tlpeaks_2cm.pdf",fname.Data()));
	  hnpeaks_2cm->SaveAs(Form("/lustrehome/bdanzi/offline_analysis/testbeam_analysis/Plots/%s/hnpeaks_2cm.pdf",fname.Data()));
	  npeaks_clustser_2cm->SaveAs(Form("/lustrehome/bdanzi/offline_analysis/testbeam_analysis/Plots/%s/hnpeaks_cluster_2cm.pdf",fname.Data()));
	  hpeaks_2cm->SaveAs(Form("/lustrehome/bdanzi/offline_analysis/testbeam_analysis/Plots/%s/hpeaks_2cm.pdf",fname.Data()));
	  bsl_2cm->SaveAs(Form("/lustrehome/bdanzi/offline_analysis/testbeam_analysis/Plots/%s/bsl_2cm.pdf",fname.Data()));
	  max_2cm->SaveAs(Form("/lustrehome/bdanzi/offline_analysis/testbeam_analysis/Plots/%s/max_2cm.pdf",fname.Data()));
	  integ_2cm->SaveAs(Form("/lustrehome/bdanzi/offline_analysis/testbeam_analysis/Plots/%s/integ_2cm.pdf",fname.Data()));
	  rms_2cm->SaveAs(Form("/lustrehome/bdanzi/offline_analysis/testbeam_analysis/Plots/%s/rms_2cm.pdf",fname.Data()));
	  cluster_population_canvas_2cm->SaveAs(Form("/lustrehome/bdanzi/offline_analysis/testbeam_analysis/Plots/%s/cluster_population_2cm.pdf",fname.Data()));
	  //min->SaveAs(Form("/lustrehome/bdanzi/offline_analysis/testbeam_analysis/Plots/%s/min_ch%i.pdf",fname.Data(),i));
	  
	  
	}
      
    }
    fclose(fp);
  }
  
  
}
