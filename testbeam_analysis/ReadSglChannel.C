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

void ReadSglChannel(){
  
  int channel=0;
  int ev=0;
  struct stat st = {0};
  char filename[100];
  TString fname("");
  bool isInteractive = false;
  bool isdoubleCanvas = false;
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
      TCanvas *max=new TCanvas("max","max",900,1200);
      TCanvas *min=new TCanvas("min","min",900,1200);
      TCanvas *bsl=new TCanvas("bsl","bsl",800,800);
      TCanvas *peaks= new TCanvas("peaks","peaks",1500,1500);
      TCanvas *integ= new TCanvas("integ","integ",900,1200);
      TCanvas *rms= new TCanvas("rms","rms",900,1200);
      
      peaks->Divide(2,3);
      //integ->Divide(4,2);
      //max->Divide(4,2);
      min->Divide(4,1);
      
      
      if(isdoubleCanvas){
	for(int i =0; i<ev; ++i){
	  TGraph *h1=(TGraph*)file->Get(Form("signal_Afterflt/CvSignal_1_ev%i",i));
	  if (h1==0x0) { continue; }
	  h1->SaveAs(Form("/lustrehome/bdanzi/offline_analysis/testbeam_analysis/Plots/%s/Waves/waves_ev%i.pdf",fname.Data(),i));
	}
      }
      else if(!isdoubleCanvas){
	
	for(int i = 0; i<100; ++i){
	  for(int j = 4; j<=channel; ++j){
	    TGraph *h1=(TGraph*)file->Get(Form("signal_Afterflt/CvSignal_1_Ch%i_ev%i",j,i));
	    if (h1==0x0) { continue; }
	    h1->SaveAs(Form("/lustrehome/bdanzi/offline_analysis/testbeam_analysis/Plots/%s/Waves/waves_ev%i_Ch%i.pdf",fname.Data(),i,j));
	  }
	}
      }
      
      for(int i = 4; i<=channel; ++i){ //looping over all channels
	//for(int i =channel; i<=channel; ++i){ //looping over one channel
	
	
	
	//TH1F *h2=(TH1F*)file->Get(Form("H-Ch%i_signal/hMaxVN_ch%i",i,i));
	//if (h2==0x0) { continue; }
	//max->cd(1);
	//h2->Draw();
	
	TH1F *h17=(TH1F*)file->Get(Form("H-Ch%i_signal/hMaxV_ch%i",i,i));
	if (h17==0x0) { continue; }
	max->cd();
	//max->cd(2);
	h17->Draw( );
	
	//TH1F *h18=(TH1F*)file->Get(Form("H-Ch%i_signal/hMaxVInR_ch%i",i,i));
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
	bsl->cd();
	h3->Draw();
	
	TH1F *h4=(TH1F*)file->Get(Form("H-Ch%i_signal/hNPeaks_ch%i",i,i));
	if (h4==0x0) { continue; }
	peaks->cd(1);
	h4->Draw();
	
	TH1F *h5=(TH1F*)file->Get(Form("H-Ch%i_signal/hHPeaks_ch%i",i,i));
	if (h5==0x0) { continue; }
	peaks->cd(2);
	h5->Draw( );
	
	TH1F *h6=(TH1F*)file->Get(Form("H-Ch%i_signal/hHNPeaks_ch%i",i,i));
	if (h6==0x0) { continue; }
	peaks->cd(3);
	h6->Draw( "colz");
	
	TH1F *h7=(TH1F*)file->Get(Form("H-Ch%i_signal/hTPeaks_ch%i",i,i));
	if (h7==0x0) { continue; }
	peaks->cd(4);
	h7->Draw( );
	
	TH1F *h8=(TH1F*)file->Get(Form("H-Ch%i_signal/hTFstPeaks_ch%i",i,i));
	if (h8==0x0) { continue; }
	peaks->cd(5);
	h8->Draw( );
	
	TH1F *h9=(TH1F*)file->Get(Form("H-Ch%i_signal/hNPeaks_1_ch%i",i,i));
	if (h9==0x0) { continue; }
	peaks->cd(6);
	h9->Draw( );
	
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
	integ->cd();
	//integ->cd(3);
	h12->Draw( );
	
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
	rms->cd();
	h16->Draw();
	
	
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
	
	
	
	bool savePlots = true;
	
	if(savePlots){
	  
	  peaks->SaveAs(Form("/lustrehome/bdanzi/offline_analysis/testbeam_analysis/Plots/%s/peaks_ch%i.pdf",fname.Data(),i));
	  bsl->SaveAs(Form("/lustrehome/bdanzi/offline_analysis/testbeam_analysis/Plots/%s/bsl_ch%i.pdf",fname.Data(),i));
	  max->SaveAs(Form("/lustrehome/bdanzi/offline_analysis/testbeam_analysis/Plots/%s/max_ch%i.pdf",fname.Data(),i));
	  integ->SaveAs(Form("/lustrehome/bdanzi/offline_analysis/testbeam_analysis/Plots/%s/integ_ch%i.pdf",fname.Data(),i));
	  rms->SaveAs(Form("/lustrehome/bdanzi/offline_analysis/testbeam_analysis/Plots/%s/rms_ch%i.pdf",fname.Data(),i));
	  min->SaveAs(Form("/lustrehome/bdanzi/offline_analysis/testbeam_analysis/Plots/%s/min_ch%i.pdf",fname.Data(),i));
	  
	  
	}
	
      }
    }
    fclose(fp);
  }
  
  
}
