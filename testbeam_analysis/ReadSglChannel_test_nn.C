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
  
  gStyle->SetOptFit(1111);
  
  int channel=0;
  int ev=0;
  struct stat st = {0};
  char filename[100];
  TString fname("");
  bool isInteractive = false;
  bool isdoubleCanvas = false;
  Float_t alpha= 0.;
  Float_t cos_alpha = 0.;
  Float_t expected_electrons =0.;
  Float_t expected_cluster =0.;
  Float_t cluster_per_cm_mip = 18.;
  Float_t drift_size =0.;
  Float_t relativistic_rise = 1.3;
  Float_t cluster_population = 1.6;
  Float_t mean_electrons_1cm = 0.;
  Float_t mean_clusters_1cm = 0.;
  Float_t rms_clusters_1cm = 0.;
  Float_t rms_electrons_1cm = 0.;
  Float_t mean_electrons_2cm = 0.;
  Float_t mean_clusters_2cm = 0.;
  Float_t rms_clusters_2cm = 0.;
  Float_t rms_electrons_2cm = 0.;
  Float_t maximum_1cm = 0.;
  Float_t maximum_2cm = 0.;
  Float_t rms_maximum_1cm = 0.;
  Float_t rms_maximum_2cm = 0.;
  Float_t integral_1cm = 0.;
  Float_t integral_2cm = 0.;
  Float_t epc_1cm = 0.;
  Float_t epc_2cm = 0.;
  Float_t rms_epc_1cm = 0.;
  Float_t rms_epc_2cm = 0.;
  Float_t bsl_1cm_var = 0.;
  Float_t bsl_2cm_var = 0.;
  Float_t rms_bsl_1cm = 0.;
  Float_t rms_bsl_2cm = 0.;
  Float_t aveph_1cm = 0.;
  Float_t aveph_2cm = 0.;
  Float_t rms_aveph_1cm = 0.;
  Float_t rms_aveph_2cm = 0.;
  Float_t rms_1cm_var = 0.;
  Float_t rms_2cm_var = 0.;
  Float_t rms_rms_1cm = 0.;
  Float_t rms_rms_2cm = 0.;
  Float_t rms_integral_1cm = 0.;
  Float_t rms_integral_2cm = 0.;
  Int_t counter_filling_electrons_1cm = 0;
  Int_t counter_filling_clusters_1cm = 0;
  Int_t counter_filling_electrons_2cm = 0;
  Int_t counter_filling_clusters_2cm = 0;

  
  TCanvas *aveph_summary_1cm = new TCanvas("aveph_summary_1cm","Aveph 1 cm cell size Drift Tubes",200,10,500,300);
  TCanvas *aveph_summary_2cm = new TCanvas("aveph_summary_2cm","Aveph 2 cm cell size Drift Tubes",200,10,500,300);
  TGraphErrors* gr_aveph_summary_1cm = new TGraphErrors(5);
  TGraphErrors* gr_aveph_summary_2cm = new TGraphErrors(3);
  gr_aveph_summary_2cm->SetTitle("Aveph 2 cm cell size Drift Tubes");
  gr_aveph_summary_2cm->SetMarkerColor(kBlue);
  gr_aveph_summary_2cm->GetXaxis()->SetTitle("Track angle (deg)");
  gr_aveph_summary_2cm->GetYaxis()->SetTitle("Average Pulse height (V)");
  gr_aveph_summary_1cm->SetTitle("Aveph 1 cm cell size Drift Tubes");
  gr_aveph_summary_1cm->SetMarkerColor(kBlue);
  gr_aveph_summary_1cm->GetYaxis()->SetTitle("Average Pulse Height (V)");
  gr_aveph_summary_1cm->GetXaxis()->SetTitle("Track angle (deg)");
  
  TCanvas *rms_summary_1cm = new TCanvas("rms_summary_1cm","Rms 1 cm cell size Drift Tubes",200,10,500,300);
  TCanvas *rms_summary_2cm = new TCanvas("rms_summary_2cm","Rms 2 cm cell size Drift Tubes",200,10,500,300);
  TGraphErrors* gr_rms_summary_1cm = new TGraphErrors(5);
  TGraphErrors* gr_rms_summary_2cm = new TGraphErrors(3);
  gr_rms_summary_2cm->SetTitle("Rms 2 cm cell size Drift Tubes");
  gr_rms_summary_2cm->SetMarkerColor(kBlue);
  gr_rms_summary_2cm->GetXaxis()->SetTitle("Track angle (deg)");
  gr_rms_summary_2cm->GetYaxis()->SetTitle("Rms (mV)");
  gr_rms_summary_1cm->SetTitle("Rms 1 cm cell size Drift Tubes");
  gr_rms_summary_1cm->SetMarkerColor(kBlue);
  gr_rms_summary_1cm->GetYaxis()->SetTitle("Rms (mV)");
  gr_rms_summary_1cm->GetXaxis()->SetTitle("Track angle (deg)");
  

  TCanvas *bsl_summary_1cm = new TCanvas("bsl_summary_1cm","Bsl 1 cm cell size Drift Tubes",200,10,500,300);
  TCanvas *bsl_summary_2cm = new TCanvas("bsl_summary_2cm","Bsl 2 cm cell size Drift Tubes",200,10,500,300);
  TGraphErrors* gr_bsl_summary_1cm = new TGraphErrors(5);
  TGraphErrors* gr_bsl_summary_2cm = new TGraphErrors(3);
  gr_bsl_summary_2cm->SetTitle("Bsl 2 cm cell size Drift Tubes");
  gr_bsl_summary_2cm->SetMarkerColor(kBlue);
  gr_bsl_summary_2cm->GetXaxis()->SetTitle("Track angle (deg)");
  gr_bsl_summary_2cm->GetYaxis()->SetTitle("Baseline (V)");
  gr_bsl_summary_1cm->SetTitle("Bsl 1 cm cell size Drift Tubes");
  gr_bsl_summary_1cm->SetMarkerColor(kBlue);
  gr_bsl_summary_1cm->GetYaxis()->SetTitle("Baseline (V)");
  gr_bsl_summary_1cm->GetXaxis()->SetTitle("Track angle (deg)");
  

  TCanvas *epc_summary_1cm = new TCanvas("epc_summary_1cm","Epc 1 cm cell size Drift Tubes",200,10,500,300);
  TCanvas *epc_summary_2cm = new TCanvas("epc_summary_2cm","Epc 2 cm cell size Drift Tubes",200,10,500,300);
  TGraphErrors* gr_epc_summary_1cm = new TGraphErrors(5);
  TGraphErrors* gr_epc_summary_2cm = new TGraphErrors(3);
  gr_epc_summary_2cm->SetTitle("Epc 2 cm cell size Drift Tubes");
  gr_epc_summary_2cm->SetMarkerColor(kBlue);
  gr_epc_summary_2cm->GetXaxis()->SetTitle("Track angle (deg)");
  gr_epc_summary_2cm->GetYaxis()->SetTitle("Electrons per Cluster");
  gr_epc_summary_1cm->SetTitle("Epc 1 cm cell size Drift Tubes");
  gr_epc_summary_1cm->SetMarkerColor(kBlue);
  gr_epc_summary_1cm->GetYaxis()->SetTitle("Electrons per Cluster");
  gr_epc_summary_1cm->GetXaxis()->SetTitle("Track angle (deg)");
  
  TCanvas *integral_summary_1cm = new TCanvas("integral_summary_1cm","Integral 1 cm cell size Drift Tubes",200,10,500,300);
  TCanvas *integral_summary_2cm = new TCanvas("integral_summary_2cm","Integral 2 cm cell size Drift Tubes",200,10,500,300);
  TGraphErrors* gr_integral_summary_1cm = new TGraphErrors(5);
  TGraphErrors* gr_integral_summary_2cm = new TGraphErrors(3);
  gr_integral_summary_2cm->SetTitle("Integral 2 cm cell size Drift Tubes");
  gr_integral_summary_2cm->SetMarkerColor(kBlue);
  gr_integral_summary_2cm->GetYaxis()->SetTitle("Charge (pC)");
  gr_integral_summary_2cm->GetXaxis()->SetTitle("Track angle (deg)");
  gr_integral_summary_1cm->SetTitle("Integral 1 cm cell size Drift Tubes");
  gr_integral_summary_1cm->SetMarkerColor(kBlue);
  gr_integral_summary_1cm->GetXaxis()->SetTitle("Track angle (deg)");
  gr_integral_summary_1cm->GetYaxis()->SetTitle("Charge (pC)");
  
  TCanvas *maximum_summary_1cm = new TCanvas("maximum_summary_1cm","Maximum 1 cm cell size Drift Tubes",200,10,500,300);
  TCanvas *maximum_summary_2cm = new TCanvas("maximum_summary_2cm","Maximum 2 cm cell size Drift Tubes",200,10,500,300);
  TGraphErrors* gr_maximum_summary_1cm = new TGraphErrors(5);
  TGraphErrors* gr_maximum_summary_2cm = new TGraphErrors(3);
  gr_maximum_summary_2cm->SetTitle("Maximum 2 cm cell size Drift Tubes");
  gr_maximum_summary_2cm->SetMarkerColor(kBlue);
  gr_maximum_summary_2cm->GetYaxis()->SetTitle("Voltage (V)");
  gr_maximum_summary_2cm->GetXaxis()->SetTitle("Track angle (deg)");
  gr_maximum_summary_1cm->SetTitle("Maximum 1 cm cell size Drift Tubes");
  gr_maximum_summary_1cm->SetMarkerColor(kBlue);
  gr_maximum_summary_1cm->GetXaxis()->SetTitle("Track angle (deg)");
  gr_maximum_summary_1cm->GetYaxis()->SetTitle("Voltage (V)");
  

  TCanvas *efficiency_electrons_1cm = new TCanvas("efficiency_electrons_1cm","Electron Finding Efficiency 1 cm cell size Drift Tubes",200,10,500,300);
  TCanvas *efficiency_electrons_2cm = new TCanvas("efficiency_electrons_2cm","Electron Finding Efficiency 2 cm cell size Drift Tubes",200,10,500,300);
  TGraphErrors* gr_efficiency_electrons_1cm = new TGraphErrors(5);
  TGraphErrors* gr_efficiency_electrons_2cm = new TGraphErrors(3);
  
  TCanvas *efficiency_clusters_1cm = new TCanvas("efficiency_clusters_1cm","Cluster Finding Efficiency 1 cm cell size Drift Tubes",200,10,500,300);
  TCanvas *efficiency_clusters_2cm = new TCanvas("efficiency_clusters_2cm","Cluster Finding Efficiency 2 cm cell size Drift Tubes",200,10,500,300);
  TGraphErrors* gr_efficiency_clusters_2cm = new TGraphErrors(3);
  TGraphErrors* gr_efficiency_clusters_1cm = new TGraphErrors(5);
  
  gr_efficiency_electrons_2cm->SetTitle("Electrons Finding Efficiency 2 cm cell size Drift Tubes");
  gr_efficiency_clusters_2cm->SetTitle("Clusters Finding Efficiency 2 cm cell size Drift Tubes");
  gr_efficiency_electrons_1cm->SetTitle("Electrons Finding Efficiency 1 cm cell size Drift Tubes");
  gr_efficiency_clusters_1cm->SetTitle("Clusters Finding Efficiency 1 cm cell size Drift Tubes");
  
  gr_efficiency_electrons_2cm->SetMarkerColor(kBlue);
  gr_efficiency_clusters_2cm->SetMarkerColor(kBlue);
  gr_efficiency_electrons_1cm->SetMarkerColor(kBlue);
  gr_efficiency_clusters_1cm->SetMarkerColor(kBlue);
  
  gr_efficiency_electrons_2cm->GetXaxis()->SetTitle("Track Angle (deg)");
  gr_efficiency_clusters_2cm->GetXaxis()->SetTitle("Track Angle (deg)");
  gr_efficiency_electrons_1cm->GetXaxis()->SetTitle("Track Angle (deg)");
  gr_efficiency_clusters_1cm->GetXaxis()->SetTitle("Track Angle (deg)");
  
  gr_efficiency_electrons_2cm->GetYaxis()->SetTitle("Measured Average Number of Electrons / Expected Number of Electrons");
  gr_efficiency_clusters_2cm->GetYaxis()->SetTitle("Measured Average Number of Clusters / Expected Number of Clusters");
  gr_efficiency_electrons_1cm->GetYaxis()->SetTitle("Measured Average Number of Electrons / Expected Number of Electrons");
  gr_efficiency_clusters_1cm->GetYaxis()->SetTitle("Measured Average Number of Clusters / Expected Number of Clusters");
  
  
  
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
      string name_file = fname.Data();
      fscanf(fp,"%d",&channel); //number of channels for which I show plots (4-14)
      fscanf(fp,"%d",&ev); //number of events for which I show plots (0-...)
      //runinputs10000events=("89" "90" "91" "92")
      //runinputs5000events=("86" "87" "88" "93" "94" "95" "96" "97" "98" "99" "127" "117")
      if((name_file =="histosTB_run_127.root") || (name_file == "histosTB_run_117.root")){
	cluster_per_cm_mip = 18.;
	printf("CHANGED CLUSTER PER CM PHYSICAL QUANTITY to %f!\n",cluster_per_cm_mip);
      }
      else{
	cluster_per_cm_mip = 12.;
	printf("CHANGED CLUSTER PER CM PHYSICAL QUANTITY to %f!\n",cluster_per_cm_mip);
      }
      
      
      if((name_file =="histosTB_run_99.root") || (name_file == "histosTB_run_117.root") || (name_file == "histosTB_run_86.root")|| (name_file == "histosTB_run_100.root")|| (name_file == "histosTB_run_101.root")|| (name_file == "histosTB_run_72.root")|| (name_file == "histosTB_run_73.root")|| (name_file == "histosTB_run_74.root")){
	alpha = 0;
	printf("CHANGED angle Alpha PHYSICAL QUANTITY to %f!\n",alpha);
	
      }
      else if((name_file =="histosTB_run_98.root")|| (name_file == "histosTB_run_87.root")|| (name_file == "histosTB_run_97.root")){
	alpha = 15;
	printf("CHANGED angle Alpha PHYSICAL QUANTITY to %f!\n",alpha);
	
	
      }
      else if((name_file =="histosTB_run_96.root")|| (name_file == "histosTB_run_88.root")|| (name_file == "histosTB_run_95.root")|| (name_file == "histosTB_run_96.root")){
	alpha = 30;
	printf("CHANGED angle Alpha PHYSICAL QUANTITY to %f!\n",alpha);
	
	
      }
      else if((name_file =="histosTB_run_94.root") || (name_file == "histosTB_run_89.root")|| (name_file == "histosTB_run_93.root")|| (name_file == "histosTB_run_94.root")){
	alpha = 45;
	printf("CHANGED angle Alpha PHYSICAL QUANTITY to %f!\n",alpha);
	
	
      }
      else if((name_file =="histosTB_run_91.root") || (name_file =="histosTB_run_127.root")|| (name_file == "histosTB_run_90.root")|| (name_file == "histosTB_run_92.root")){
	alpha = 60;
	printf("CHANGED angle Alpha PHYSICAL QUANTITY to %f!\n",alpha);
	
      }
      
      cos_alpha = TMath::Cos(alpha*TMath::DegToRad());

	mean_electrons_1cm = 0.;
  	mean_clusters_1cm = 0.;
  	rms_clusters_1cm = 0.;
  	rms_electrons_1cm = 0.;
  	mean_electrons_2cm = 0.;
  	mean_clusters_2cm = 0.;
  	rms_clusters_2cm = 0.;
  	rms_electrons_2cm = 0.;
  	maximum_1cm = 0.;
  	maximum_2cm = 0.;
  	rms_maximum_1cm = 0.;
  	rms_maximum_2cm = 0.;
  	integral_1cm = 0.;
  	integral_2cm = 0.;
  	epc_1cm = 0.;
  	epc_2cm = 0.;
  	rms_epc_1cm = 0.;
  	rms_epc_2cm = 0.;
  	bsl_1cm_var = 0.;
  	bsl_2cm_var = 0.;
  	rms_bsl_1cm = 0.;
  	rms_bsl_2cm = 0.;
  	aveph_1cm = 0.;
  	aveph_2cm = 0.;
  	rms_aveph_1cm = 0.;
  	rms_aveph_2cm = 0.;
  	rms_1cm_var = 0.;
  	rms_2cm_var = 0.;
  	rms_rms_1cm = 0.;
  	rms_rms_2cm = 0.;
  	rms_integral_1cm = 0.;
  	rms_integral_2cm = 0.;
      
      if (stat(Form("/lustrehome/bdanzi/offline_analysis/testbeam_analysis/Plots/%s/",fname.Data()), &st) == -1) {
	mkdir(Form("/lustrehome/bdanzi/offline_analysis/testbeam_analysis/Plots/%s/",fname.Data()), 0700);
      }
      if (stat(Form("/lustrehome/bdanzi/offline_analysis/testbeam_analysis/Plots/%s/Waves/",fname.Data()), &st) == -1) {
	mkdir(Form("/lustrehome/bdanzi/offline_analysis/testbeam_analysis/Plots/%s/Waves/",fname.Data()), 0700);
      }
      
      TFile *file = new TFile(fname.Data(),"read");
      TCanvas *max_1cm=new TCanvas("max_1","max",3500,1500);
      TCanvas *waveform=new TCanvas("waveform","waveform",3500,1500);
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
      TCanvas *timediff_1cm= new TCanvas("timediff_1","timediff",3500,1500);
      TCanvas *timediff_2cm= new TCanvas("timediff_2","timediff",3500,1500);
      TCanvas *timediff_clust_1cm= new TCanvas("timediff_clust_1","timediff_clust",3500,1500);
      TCanvas *timediff_clust_2cm= new TCanvas("timediff_clust_2","timediff_clust",3500,1500);
      
	
      cluster_population_canvas_1cm->Divide(2,3);
      cluster_population_canvas_2cm->Divide(2,2);
      timediff_clust_1cm->Divide(2,3);
      timediff_clust_2cm->Divide(2,2);
      timediff_1cm->Divide(2,3);
      timediff_2cm->Divide(2,2);
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
	   //h1->SaveAs(Form("/lustrehome/bdanzi/offline_analysis/testbeam_analysis/Plots/%s/Waves/waves_ev%i.png",fname.Data(),i));
	}
      }
      else if(!isdoubleCanvas){
	
	for(int i = 0; i<ev; ++i){
	  for(int j = 4; j<=channel; ++j){

		string nn_file = "";
		nn_file = Form("nn_ch%d_alpha%f",j,alpha) + fname.Data() + ".txt";
		ofstream myfile_nn (nn_file);
 		if (myfile_nn.is_open())
  		{
  		  
		  drift_size = 0.8;
      	  //δ cluster/cm (M.I.P.) * drift tube size [cm] * 1.3 (relativisticrise) * 1.6 electrons/cluster * 1/cos(α)
      	  expected_electrons = cluster_per_cm_mip * drift_size * relativistic_rise * cluster_population * 1/cos_alpha;
      	  myfile_nn << expected_electrons;
		  
    	  myfile_nn << "This is another line.\n";
		  TGraph *h1=(TGraph*)file->Get(Form("signal_Afterflt/CvSignal_1_Ch%i_ev%i",j,i));
	      if (h1==0x0) { continue; }
		  
	      h1->SetTitle(Form("Channel %d event %d Alpha %.1f Run %s",j,i,alpha, fname.Data()));
	      //h1->GetYaxis()->SetRangeUser(-0.1,0.7);
	      h1->SaveAs(Form("/lustrehome/bdanzi/offline_analysis/testbeam_analysis/Plots/%s/Waves/waves_ev%i_Ch%i.pdf",fname.Data(),i,j));
	      // h1->SaveAs(Form("/lustrehome/bdanzi/offline_analysis/testbeam_analysis/Plots/%s/Waves/waves_ev%i_Ch%i.png",fname.Data(),i,j));
		  myfile_nn.close();
  		}

  		else cout << "Unable to open file"; 
	    
	  }
	}
      }
      
      drift_size = 0.8;
      
      //δ cluster/cm (M.I.P.) * drift tube size [cm] * 1.3 (relativisticrise) * 1.6 electrons/cluster * 1/cos(α)
      expected_electrons = cluster_per_cm_mip * drift_size * relativistic_rise * cluster_population * 1/cos_alpha;
      expected_cluster = cluster_per_cm_mip * drift_size * relativistic_rise * 1/cos_alpha;
      
      for(int i = 4; i<=9; ++i){ //looping over all channels
	//for(int i =channel; i<=channel; ++i){ //looping over one channel
	
	//2cm/1cm 1,8/0,8 = 2,25
	
	//TH1F *h2=(TH1F*)file->Get(Form("H-Ch%i_signal/hMaxVN_ch%i",i,i));
	//if (h2==0x0) { continue; }
	//max->cd(1);
	//h2->Draw();
	
	TH1F *h17=(TH1F*)file->Get(Form("H-Ch%i_signal/hMaxVInR_ch%i",i,i));
	if (h17==0x0) { continue; }
	max_1cm->cd(i-3);
	if((name_file =="histosTB_run_99.root") || (name_file == "histosTB_run_98.root") || (name_file == "histosTB_run_96.root")|| (name_file == "histosTB_run_94.root")|| (name_file == "histosTB_run_91.root")){
	  maximum_1cm = maximum_1cm+ h17->GetMean();
	  rms_maximum_1cm = rms_maximum_1cm + h17->GetRMS();
	}
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
	if((name_file =="histosTB_run_99.root") || (name_file == "histosTB_run_98.root") || (name_file == "histosTB_run_96.root")|| (name_file == "histosTB_run_94.root")|| (name_file == "histosTB_run_91.root")){
	  bsl_1cm_var = bsl_1cm_var + h3->GetMean();
	  rms_bsl_1cm = rms_bsl_1cm + h3->GetRMS();
	}
	gPad->SetLogy(0);
	gPad->SetLogx(0);
	h3->GetXaxis()->SetRangeUser(-0.5,-0.4);
	h3->Draw("same");
	
	
	TH1F *h4=(TH1F*)file->Get(Form("H-Ch%i_signal/hNPeaks_ch%i",i,i));
	if (h4==0x0) { continue; }
	npeaks_1cm->cd(i-3);
	h4->GetXaxis()->SetRangeUser(0.,90.);
	if((name_file =="histosTB_run_99.root") || (name_file == "histosTB_run_98.root") || (name_file == "histosTB_run_96.root")|| (name_file == "histosTB_run_94.root")|| (name_file == "histosTB_run_91.root")){
	  mean_electrons_1cm = mean_electrons_1cm + h4->GetMean();
	  rms_electrons_1cm = rms_electrons_1cm + h4->GetRMS();
	}
	//h4->Fit("landau");
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
	TPaveText *pt_1cm_alpha = new TPaveText(0.72,0.13,0.8,0.17,"NDC");
	pt_1cm_alpha->SetTextSize(0.04);
	pt_1cm_alpha->SetTextColor(kRed);
	pt_1cm_alpha->SetFillColor(0);
	pt_1cm_alpha->SetTextAlign(12);
	pt_1cm_alpha->AddText(Form("Alpha angle (deg): %.1f",alpha));
	pt_1cm_alpha->Draw("same");
	
	TH1F *h20=(TH1F*)file->Get(Form("H-Ch%i_signal/hNPeaks_clust_ch%i",i,i));
	if (h20==0x0) { continue; }
	npeaks_clustser_1cm->cd(i-3);
	h20->GetXaxis()->SetRangeUser(0.,90.);
	if((name_file =="histosTB_run_99.root") || (name_file == "histosTB_run_98.root") || (name_file == "histosTB_run_96.root")|| (name_file == "histosTB_run_94.root")|| (name_file == "histosTB_run_91.root")){
	  mean_clusters_1cm = mean_clusters_1cm + h20->GetMean();
	  rms_clusters_1cm = rms_clusters_1cm + h20->GetRMS();
	}
	if(expected_cluster>0){
	  h20->Fit("gaus");
	  //printf("Hello");
	}
	else{
	  TF1 *f11=new TF1("f11","[0]*TMath::Poisson(x,[1])",0,60.);                                                          
	  f11->SetParName(0,"Normalisation");
	  f11->SetParName(1,"#mu");
	  f11->SetParameters(0,3);
	  f11->SetParameters(0,h20->GetMean());
	  h20->Fit("f11","R");
	}
	h20->Draw("same");
	TPaveText *pt_1cm_cluster = new TPaveText(0.72,0.2,0.8,0.35,"NDC");
	pt_1cm_cluster->SetTextSize(0.04);
	pt_1cm_cluster->SetTextColor(kRed);
	pt_1cm_cluster->SetTextAlign(12);
	pt_1cm_cluster->SetFillColor(0);
	pt_1cm_cluster->AddText(Form("Expected Clusters: %.1f",expected_cluster));
	gPad->SetLogy(1);
	pt_1cm_cluster->Draw("same");
	pt_1cm_alpha->Draw("same");
	
	TH1F *h31=(TH1F*)file->Get(Form("H-Ch%i_signal/hNElectrons_per_cluster_ch%i",i,i));
	if (h31==0x0) { continue; }
	cluster_population_canvas_1cm->cd(i-3);
	if((name_file =="histosTB_run_99.root") || (name_file == "histosTB_run_98.root") || (name_file == "histosTB_run_96.root")|| (name_file == "histosTB_run_94.root")|| (name_file == "histosTB_run_91.root")){
	  epc_1cm = epc_1cm+ h31->GetMean();
	  rms_epc_1cm = rms_epc_1cm + h31->GetRMS();
	}
	//h31->GetXaxis()->SetRangeUser(0.,5.);
	//h31->Fit("expo");
	h31->Draw("same");
	TPaveText *cluster_population_1cm = new TPaveText(0.72,0.2,0.8,0.35,"NDC");
	cluster_population_1cm->SetTextSize(0.04);
	cluster_population_1cm->SetTextColor(kRed);
	cluster_population_1cm->SetTextAlign(12);
	cluster_population_1cm->SetFillColor(0);
	cluster_population_1cm->AddText(Form("Expected Electrons per Cluster: %.1f",cluster_population));
	gPad->SetLogy(1);
	pt_1cm_alpha->Draw("same");
	cluster_population_1cm->Draw("same");
	
	TH1F *h32=(TH1F*)file->Get(Form("H-Ch%i_signal/hTimeDifference_ch%i",i,i));
	if (h32==0x0) { continue; }
	timediff_1cm->cd(i-3);
	TF1  *f4 = new TF1("f4","[0]*exp(-x/[1])",0,10);
	//TF1  *f4 = new TF1("f4","[0]*exp(-x/[1])",0,12);
	f4->SetParameters(0,4000.);
	f4->SetParameters(1,3);
	h32->Fit("f4","R");
	gPad->SetLogy(1);
	h32->Draw("same");
	//pt_1cm_alpha->Draw("same");
	//h31->GetXaxis()->SetRangeUser(0.,5.);
	//h31->Fit("expo");
	pt_1cm_alpha->Draw("same");
	
	
	TH1F *h33=(TH1F*)file->Get(Form("H-Ch%i_signal/hTimeDifference_clust_ch%i",i,i));
	if (h33==0x0) { continue; }
	timediff_clust_1cm->cd(i-3);
	//pt_1cm_alpha->Draw("same");
	TF1  *f1 = new TF1("f1","[0]*exp(-x/[1])",10,40);
	f1->SetParameters(0,1000);
	f1->SetParameters(1,10);
	h33->Fit("f1","R");
	gPad->SetLogy(1);
	h33->Draw("same");
	pt_1cm_alpha->Draw("same");
	
	TH1F *h5=(TH1F*)file->Get(Form("H-Ch%i_signal/hHPeaks_ch%i",i,i));
	if (h5==0x0) { continue; }
	hpeaks_1cm->cd(i-3);
	if((name_file =="histosTB_run_99.root") || (name_file == "histosTB_run_98.root") || (name_file == "histosTB_run_96.root")|| (name_file == "histosTB_run_94.root")|| (name_file == "histosTB_run_91.root")){
	  aveph_1cm = aveph_1cm+ h5->GetMean();
	  rms_aveph_1cm = rms_aveph_1cm + h5->GetRMS();
	}
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
	if((name_file =="histosTB_run_99.root") || (name_file == "histosTB_run_98.root") || (name_file == "histosTB_run_96.root")|| (name_file == "histosTB_run_94.root")|| (name_file == "histosTB_run_91.root")){
	  integral_1cm = integral_1cm+ h12->GetMean();
	  rms_integral_1cm = rms_integral_1cm + h12->GetRMS();
	}
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
	if((name_file =="histosTB_run_99.root") || (name_file == "histosTB_run_98.root") || (name_file == "histosTB_run_96.root")|| (name_file == "histosTB_run_94.root")|| (name_file == "histosTB_run_91.root")){
	  rms_1cm_var = rms_1cm_var + h16->GetMean();
	  rms_rms_1cm = rms_rms_1cm + h16->GetRMS();
	}
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
	
	
	
	
      } //loop on 1 cm tubes
      if((name_file =="histosTB_run_99.root") || (name_file == "histosTB_run_98.root") || (name_file == "histosTB_run_96.root")|| (name_file == "histosTB_run_94.root")|| (name_file == "histosTB_run_91.root")){
	aveph_1cm = aveph_1cm/6.;
	rms_aveph_1cm = rms_aveph_1cm/6.;
	
	rms_1cm_var = rms_1cm_var/6.;
	rms_rms_1cm = rms_rms_1cm/6.;
    
	epc_1cm = epc_1cm/6.;
	rms_epc_1cm = rms_epc_1cm/6.;
	
	integral_1cm = (integral_1cm/6.) * (1000/500);
	rms_integral_1cm = (rms_integral_1cm/6.) * (1000/500);

	bsl_1cm_var = bsl_1cm_var/6.;
	rms_bsl_1cm = rms_bsl_1cm/6.;
	
	maximum_1cm = maximum_1cm/6.;
	rms_maximum_1cm = rms_maximum_1cm/6.;
	
	mean_clusters_1cm = mean_clusters_1cm/6.;
	rms_clusters_1cm = rms_clusters_1cm/6.;
	
	mean_electrons_1cm = mean_electrons_1cm/6.;
	rms_electrons_1cm = rms_electrons_1cm/6.;
	
	//Normalization
	
	mean_clusters_1cm = mean_clusters_1cm/expected_cluster;
	rms_clusters_1cm = rms_clusters_1cm/expected_cluster;
	
	mean_electrons_1cm = mean_electrons_1cm/expected_electrons;
	rms_electrons_1cm = rms_electrons_1cm/expected_electrons;
	
	aveph_summary_1cm->cd();   
	gr_aveph_summary_1cm->SetPoint(counter_filling_electrons_1cm,alpha,aveph_1cm);
	gr_aveph_summary_1cm->SetPointError(counter_filling_electrons_1cm,0.,rms_aveph_1cm);
	gr_aveph_summary_1cm->Draw("sameAC*");
	
 
	rms_summary_1cm->cd();   
	gr_rms_summary_1cm->SetPoint(counter_filling_electrons_1cm,alpha,rms_1cm_var);
	gr_rms_summary_1cm->SetPointError(counter_filling_electrons_1cm,0.,rms_rms_1cm);
	gr_rms_summary_1cm->Draw("sameAC*");
	
	epc_summary_1cm->cd();   
	gr_epc_summary_1cm->SetPoint(counter_filling_electrons_1cm,alpha,epc_1cm);
	gr_epc_summary_1cm->SetPointError(counter_filling_electrons_1cm,0.,rms_epc_1cm);
	gr_epc_summary_1cm->Draw("sameAC*");
	
	bsl_summary_1cm->cd();   
	gr_bsl_summary_1cm->SetPoint(counter_filling_electrons_1cm,alpha,bsl_1cm_var);
	gr_bsl_summary_1cm->SetPointError(counter_filling_electrons_1cm,0.,rms_bsl_1cm);
	gr_bsl_summary_1cm->Draw("sameAC*");
	
	maximum_summary_1cm->cd();   
	gr_maximum_summary_1cm->SetPoint(counter_filling_electrons_1cm,alpha,maximum_1cm);
	gr_maximum_summary_1cm->SetPointError(counter_filling_electrons_1cm,0.,rms_maximum_1cm);
	gr_maximum_summary_1cm->Draw("sameAC*");
	
	integral_summary_1cm->cd();
	gr_integral_summary_1cm->SetPoint(counter_filling_electrons_1cm,alpha,integral_1cm);
	gr_integral_summary_1cm->SetPointError(counter_filling_electrons_1cm,0.,rms_integral_1cm);
	gr_integral_summary_1cm->Draw("sameAC*");
	
	efficiency_electrons_1cm->cd();	 
	gr_efficiency_electrons_1cm->SetPoint(counter_filling_electrons_1cm,alpha,mean_electrons_1cm);
	gr_efficiency_electrons_1cm->SetPointError(counter_filling_electrons_1cm,0.,rms_electrons_1cm);
	gr_efficiency_electrons_1cm->Draw("sameAC*");
	
	efficiency_clusters_1cm->cd();
	gr_efficiency_clusters_1cm->SetPoint(counter_filling_clusters_1cm,alpha,mean_clusters_1cm);
	gr_efficiency_clusters_1cm->SetPointError(counter_filling_clusters_1cm,0.,rms_clusters_1cm);
	
	counter_filling_clusters_1cm = counter_filling_clusters_1cm +1;
	counter_filling_electrons_1cm = counter_filling_electrons_1cm +1;
	gr_efficiency_clusters_1cm->Draw("sameAC*");
	
	if(name_file =="histosTB_run_99.root"){
	  efficiency_clusters_1cm->SaveAs("/lustrehome/bdanzi/offline_analysis/testbeam_analysis/Plots/efficiency_clusters_1cm.pdf",fname.Data());
	  efficiency_electrons_1cm->SaveAs("/lustrehome/bdanzi/offline_analysis/testbeam_analysis/Plots/efficiency_electrons_1cm.pdf",fname.Data());
	  epc_summary_1cm->SaveAs("/lustrehome/bdanzi/offline_analysis/testbeam_analysis/Plots/epc_summary_1cm.pdf",fname.Data());
	  integral_summary_1cm->SaveAs("/lustrehome/bdanzi/offline_analysis/testbeam_analysis/Plots/integral_summary_1cm.pdf",fname.Data());
	  maximum_summary_1cm->SaveAs("/lustrehome/bdanzi/offline_analysis/testbeam_analysis/Plots/maximum_summary_1cm.pdf",fname.Data());
	  bsl_summary_1cm->SaveAs("/lustrehome/bdanzi/offline_analysis/testbeam_analysis/Plots/bsl_summary_1cm.pdf",fname.Data());
	  rms_summary_1cm->SaveAs("/lustrehome/bdanzi/offline_analysis/testbeam_analysis/Plots/rms_summary_1cm.pdf",fname.Data());
	  aveph_summary_1cm->SaveAs("/lustrehome/bdanzi/offline_analysis/testbeam_analysis/Plots/aveph_summary_1cm.pdf",fname.Data());

	  efficiency_clusters_1cm->SaveAs("/lustrehome/bdanzi/offline_analysis/testbeam_analysis/Plots/efficiency_clusters_1cm.png",fname.Data());
      efficiency_electrons_1cm->SaveAs("/lustrehome/bdanzi/offline_analysis/testbeam_analysis/Plots/efficiency_electrons_1cm.png",fname.Data());
      epc_summary_1cm->SaveAs("/lustrehome/bdanzi/offline_analysis/testbeam_analysis/Plots/epc_summary_1cm.png",fname.Data());
      integral_summary_1cm->SaveAs("/lustrehome/bdanzi/offline_analysis/testbeam_analysis/Plots/integral_summary_1cm.png",fname.Data());
      maximum_summary_1cm->SaveAs("/lustrehome/bdanzi/offline_analysis/testbeam_analysis/Plots/maximum_summary_1cm.png",fname.Data());
      bsl_summary_1cm->SaveAs("/lustrehome/bdanzi/offline_analysis/testbeam_analysis/Plots/bsl_summary_1cm.png",fname.Data());
      rms_summary_1cm->SaveAs("/lustrehome/bdanzi/offline_analysis/testbeam_analysis/Plots/rms_summary_1cm.png",fname.Data());
      aveph_summary_1cm->SaveAs("/lustrehome/bdanzi/offline_analysis/testbeam_analysis/Plots/aveph_summary_1cm.png",fname.Data());
      
	}
      }
      
      
      bool savePlots = true;
      
      if(savePlots){
	timediff_clust_1cm->SaveAs(Form("/lustrehome/bdanzi/offline_analysis/testbeam_analysis/Plots/%s/clust_difference_1cm.pdf",fname.Data()));
	timediff_1cm->SaveAs(Form("/lustrehome/bdanzi/offline_analysis/testbeam_analysis/Plots/%s/electrons_difference_1cm.pdf",fname.Data()));
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
	
	timediff_clust_1cm->SaveAs(Form("/lustrehome/bdanzi/offline_analysis/testbeam_analysis/Plots/%s/clust_difference_1cm.png",fname.Data()));
    timediff_1cm->SaveAs(Form("/lustrehome/bdanzi/offline_analysis/testbeam_analysis/Plots/%s/electrons_difference_1cm.png",fname.Data()));
    npeaks_1cm->SaveAs(Form("/lustrehome/bdanzi/offline_analysis/testbeam_analysis/Plots/%s/npeaks_1cm.png",fname.Data()));
    tpeaks_1cm->SaveAs(Form("/lustrehome/bdanzi/offline_analysis/testbeam_analysis/Plots/%s/tpeaks_1cm.png",fname.Data()));
    tfpeaks_1cm->SaveAs(Form("/lustrehome/bdanzi/offline_analysis/testbeam_analysis/Plots/%s/tfpeaks_1cm.png",fname.Data()));
    tlpeaks_1cm->SaveAs(Form("/lustrehome/bdanzi/offline_analysis/testbeam_analysis/Plots/%s/tlpeaks_1cm.png",fname.Data()));
    hnpeaks_1cm->SaveAs(Form("/lustrehome/bdanzi/offline_analysis/testbeam_analysis/Plots/%s/hnpeaks_1cm.png",fname.Data()));
    hpeaks_1cm->SaveAs(Form("/lustrehome/bdanzi/offline_analysis/testbeam_analysis/Plots/%s/hpeaks_1cm.png",fname.Data()));
    bsl_1cm->SaveAs(Form("/lustrehome/bdanzi/offline_analysis/testbeam_analysis/Plots/%s/bsl_1cm.png",fname.Data()));
    max_1cm->SaveAs(Form("/lustrehome/bdanzi/offline_analysis/testbeam_analysis/Plots/%s/max_1cm.png",fname.Data()));
    integ_1cm->SaveAs(Form("/lustrehome/bdanzi/offline_analysis/testbeam_analysis/Plots/%s/integ_1cm.png",fname.Data()));
    rms_1cm->SaveAs(Form("/lustrehome/bdanzi/offline_analysis/testbeam_analysis/Plots/%s/rms_1cm.png",fname.Data()));
    //min->SaveAs(Form("/lustrehome/bdanzi/offline_analysis/testbeam_analysis/Plots/%s/min_ch%i.png",fname.Data(),i));
    npeaks_clustser_1cm->SaveAs(Form("/lustrehome/bdanzi/offline_analysis/testbeam_analysis/Plots/%s/hnpeaks_cluster_1cm.png",fname.Data()));
    cluster_population_canvas_1cm->SaveAs(Form("/lustrehome/bdanzi/offline_analysis/testbeam_analysis/Plots/%s/cluster_population_1cm.png",fname.Data()));
    
      }
      
      
      drift_size = 1.8;
      expected_electrons = cluster_per_cm_mip * drift_size*relativistic_rise * cluster_population * 1/cos_alpha;
      expected_cluster = cluster_per_cm_mip * drift_size*relativistic_rise * 1/cos_alpha;
      
      for(int i = 10; i<=channel; ++i){ //looping over all channels
	//for(int i =channel; i<=channel; ++i){ //looping over one channel
	
	
	//TH1F *h2=(TH1F*)file->Get(Form("H-Ch%i_signal/hMaxVN_ch%i",i,i));
	//if (h2==0x0) { continue; }
	//max->cd(1);
	//h2->Draw();
	
	
	TH1F *h17=(TH1F*)file->Get(Form("H-Ch%i_signal/hMaxVInR_ch%i",i,i));
	if (h17==0x0) { continue; }
	max_2cm->cd(i-9);
	if((name_file =="histosTB_run_99.root") || (name_file == "histosTB_run_98.root") || (name_file == "histosTB_run_96.root")|| (name_file == "histosTB_run_94.root")|| (name_file == "histosTB_run_91.root")){
	  maximum_2cm = maximum_2cm + h17->GetMean();
	  rms_maximum_2cm = rms_maximum_2cm + h17->GetRMS();
	}
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
	if((name_file =="histosTB_run_99.root") || (name_file == "histosTB_run_98.root") || (name_file == "histosTB_run_96.root")|| (name_file == "histosTB_run_94.root")|| (name_file == "histosTB_run_91.root")){
	  bsl_2cm_var = bsl_2cm_var+ h3->GetMean();
	  rms_bsl_2cm = rms_bsl_2cm + h3->GetRMS();
	}
	gPad->SetLogy(0);
	gPad->SetLogx(0);
	h3->GetXaxis()->SetRangeUser(-0.5,-0.4);
	h3->Draw("same");
	
	
	
	TH1F *h4=(TH1F*)file->Get(Form("H-Ch%i_signal/hNPeaks_ch%i",i,i));
	if (h4==0x0) { continue; }
	npeaks_2cm->cd(i-9);
	h4->GetXaxis()->SetRangeUser(0.,200.);
	if((name_file =="histosTB_run_99.root") || (name_file == "histosTB_run_98.root") || (name_file == "histosTB_run_96.root")|| (name_file == "histosTB_run_94.root")|| (name_file == "histosTB_run_91.root")){
	  mean_electrons_2cm = mean_electrons_2cm + h4->GetMean();
	  rms_electrons_2cm = rms_electrons_2cm + h4->GetRMS();
	}
	//h4->Fit("landau");
	gPad->SetLogy(1);
	h4->Draw("same");
	TPaveText *pt_2cm = new TPaveText(0.72,0.2,0.8,0.35,"NDC");
	pt_2cm->SetTextSize(0.04);
	pt_2cm->SetTextColor(kRed);
	pt_2cm->SetFillColor(0);
	pt_2cm->SetTextAlign(12);
	pt_2cm->AddText(Form("Expected elecrons: %.1f",expected_electrons));
	pt_2cm->Draw("same");
	TPaveText *pt_2cm_alpha = new TPaveText(0.72,0.13,0.8,0.17,"NDC");
	pt_2cm_alpha->SetTextSize(0.04);
	pt_2cm_alpha->SetTextColor(kRed);
	pt_2cm_alpha->SetFillColor(0);
	pt_2cm_alpha->SetTextAlign(12);
	pt_2cm_alpha->AddText(Form("Alpha angle (deg): %.1f",alpha));
	pt_2cm_alpha->Draw("same");
	
	TH1F *h5=(TH1F*)file->Get(Form("H-Ch%i_signal/hHPeaks_ch%i",i,i));
	if (h5==0x0) { continue; }
	hpeaks_2cm->cd(i-9);
	if((name_file =="histosTB_run_99.root") || (name_file == "histosTB_run_98.root") || (name_file == "histosTB_run_96.root")|| (name_file == "histosTB_run_94.root")|| (name_file == "histosTB_run_91.root")){
	  aveph_2cm = aveph_2cm+ h5->GetMean();
	  rms_aveph_2cm = rms_aveph_2cm + h5->GetRMS();
	}
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
	
	TH1F *h33=(TH1F*)file->Get(Form("H-Ch%i_signal/hTimeDifference_clust_ch%i",i,i));
	if (h33==0x0) { continue; }
	timediff_clust_2cm->cd(i-9);
	TF1  *f2 = new TF1("f2","[0]*exp(-x/[1])",10,40);
	f2->SetParameters(0,4000);
	f2->SetParameters(1,10);
	h33->Fit("f2","R");
	gPad->SetLogy(1);
	//pt_2cm_alpha->Draw("same");
	//h33->Fit("expo");
	h33->Draw("same");
	pt_2cm_alpha->Draw("same");
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
	if((name_file =="histosTB_run_99.root") || (name_file == "histosTB_run_98.root") || (name_file == "histosTB_run_96.root")|| (name_file == "histosTB_run_94.root")|| (name_file == "histosTB_run_91.root")){
	  integral_2cm = integral_2cm + h12->GetMean();
	  rms_integral_2cm = rms_integral_2cm + h12->GetRMS();
	}
	
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
	rms_2cm->cd(i-9);
	if((name_file =="histosTB_run_99.root") || (name_file == "histosTB_run_98.root") || (name_file == "histosTB_run_96.root")|| (name_file == "histosTB_run_94.root")|| (name_file == "histosTB_run_91.root")){
	  rms_2cm_var = rms_2cm_var + h16->GetMean();
	  rms_rms_2cm = rms_rms_2cm + h16->GetRMS();
	}
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
	if((name_file =="histosTB_run_99.root") || (name_file == "histosTB_run_98.root") || (name_file == "histosTB_run_96.root")|| (name_file == "histosTB_run_94.root")|| (name_file == "histosTB_run_91.root")){
	  mean_clusters_2cm = mean_clusters_2cm + h20->GetMean();
	  rms_clusters_2cm = rms_clusters_2cm + h20->GetRMS();
	}
	if(expected_cluster>0){
	  h20->Fit("gaus");
	  //printf("Hello");
	}
	else{
	  TF1 *f10=new TF1("f10","[0]*TMath::Poisson(x,[1])",0,90.);                                                          
	  f10->SetParName(0,"Normalisation");
	  f10->SetParName(1,"#mu");
	  f10->SetParameters(0,1000);
	  f10->SetParameters(1,h20->GetMean());
	  h20->Fit("f10","R");
	}
	gPad->SetLogy(1);
	h20->Draw("same");
	TPaveText *pt_2cm_cluster = new TPaveText(0.72,0.2,0.8,0.35,"NDC");
	pt_2cm_cluster->SetTextSize(0.04);
	pt_2cm_cluster->SetTextColor(kRed);
	pt_2cm_cluster->SetFillColor(0);
	pt_2cm_cluster->SetTextAlign(12);
	pt_2cm_cluster->AddText(Form("Expected Clusters: %.1f",expected_cluster));
	pt_2cm_cluster->Draw("same");
	pt_2cm_alpha->Draw("same");
	
	TH1F *h31=(TH1F*)file->Get(Form("H-Ch%i_signal/hNElectrons_per_cluster_ch%i",i,i));
	if (h31==0x0) { continue; }
	cluster_population_canvas_2cm->cd(i-9);
	if((name_file =="histosTB_run_99.root") || (name_file == "histosTB_run_98.root") || (name_file == "histosTB_run_96.root")|| (name_file == "histosTB_run_94.root")|| (name_file == "histosTB_run_91.root")){
	  epc_2cm = epc_2cm+ h31->GetMean();
	  rms_epc_2cm = rms_epc_2cm + h31->GetRMS();
	}
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
	pt_2cm_alpha->Draw("same");
	cluster_population_2cm->Draw("same");
	
	
	TH1F *h32=(TH1F*)file->Get(Form("H-Ch%i_signal/hTimeDifference_ch%i",i,i));
	if (h32==0x0) { continue; }
	timediff_2cm->cd(i-9);
	TF1  *f3 = new TF1("f3","[p0]*exp(-x/[p1])",0,10);
	//TF1  *f3 = new TF1("f3","[0]*exp(-x/[1])",0,12);
	f3->SetParameters(0,4000.);
	f3->SetParameters(1,3);
	h32->Fit("f3","R");
	h32->Draw("same");
	gPad->SetLogy(1);
	pt_2cm_alpha->Draw("same");
	//h31->GetXaxis()->SetRangeUser(0.,5.);
	//h31->Fit("expo");
	
      }
      if((name_file =="histosTB_run_99.root") || (name_file == "histosTB_run_98.root") || (name_file == "histosTB_run_96.root")|| (name_file == "histosTB_run_94.root")|| (name_file == "histosTB_run_91.root")){
	aveph_2cm = aveph_2cm/3.;
	rms_aveph_2cm = rms_aveph_2cm/3.;
	
	epc_2cm = epc_2cm/3.;
	rms_epc_2cm = rms_epc_2cm/3.;
	
	integral_2cm = (integral_2cm/3.)* (1000/500);
	rms_integral_2cm = (rms_integral_2cm/3.)* (1000/500);
	
	bsl_2cm_var = bsl_2cm_var/3.;
	rms_bsl_2cm = rms_bsl_2cm/3.;

	rms_2cm_var = rms_2cm_var/3.;
	rms_rms_2cm = rms_rms_2cm/3.;
	
	mean_clusters_2cm = mean_clusters_2cm/3.;
	rms_clusters_2cm = rms_clusters_2cm/3.;
	
	mean_electrons_2cm = mean_electrons_2cm/3.;
	rms_electrons_2cm = rms_electrons_2cm/3.;

	maximum_2cm = maximum_2cm/3.;
	rms_maximum_2cm = rms_maximum_2cm/3.;
	//Normalization

	mean_clusters_2cm = mean_clusters_2cm/expected_cluster;
	mean_electrons_2cm = mean_electrons_2cm/expected_electrons;
	rms_clusters_2cm = rms_clusters_2cm/expected_cluster;
	rms_electrons_2cm = rms_electrons_2cm/expected_electrons;
	aveph_summary_2cm->cd();  
	gr_aveph_summary_2cm->SetPoint(counter_filling_electrons_2cm,alpha,aveph_2cm);
	gr_aveph_summary_2cm->SetPointError(counter_filling_electrons_2cm,0.,rms_aveph_2cm);
	gr_aveph_summary_2cm->Draw("sameAC*");

	
	rms_summary_2cm->cd();  
	gr_rms_summary_2cm->SetPoint(counter_filling_electrons_2cm,alpha,rms_2cm_var);
	gr_rms_summary_2cm->SetPointError(counter_filling_electrons_2cm,0.,rms_rms_2cm);
	gr_rms_summary_2cm->Draw("sameAC*");
	
	
	epc_summary_2cm->cd();  
	gr_epc_summary_2cm->SetPoint(counter_filling_electrons_2cm,alpha,epc_2cm);
	gr_epc_summary_2cm->SetPointError(counter_filling_electrons_2cm,0.,rms_epc_2cm);
	gr_epc_summary_2cm->Draw("sameAC*");
	
	integral_summary_2cm->cd();
	//gPad->SetLogy(1); 
	gr_integral_summary_2cm->SetPoint(counter_filling_electrons_2cm,alpha,integral_2cm);
	gr_integral_summary_2cm->SetPointError(counter_filling_electrons_2cm,0.,rms_integral_2cm);
	gr_integral_summary_2cm->Draw("sameAC*");
	
	maximum_summary_2cm->cd();	 
	gr_maximum_summary_2cm->SetPoint(counter_filling_electrons_2cm,alpha,maximum_2cm);
	gr_maximum_summary_2cm->SetPointError(counter_filling_electrons_2cm,0.,rms_maximum_2cm);
	gr_maximum_summary_2cm->Draw("sameAC*");
	
	bsl_summary_2cm->cd();  
	gr_bsl_summary_2cm->SetPoint(counter_filling_electrons_2cm,alpha,bsl_2cm_var);
	gr_bsl_summary_2cm->SetPointError(counter_filling_electrons_2cm,0.,rms_bsl_2cm);
	gr_bsl_summary_2cm->Draw("sameAC*");
	
	
	efficiency_electrons_2cm->cd();  
	gr_efficiency_electrons_2cm->SetPoint(counter_filling_electrons_2cm,alpha,mean_electrons_2cm);
	gr_efficiency_electrons_2cm->SetPointError(counter_filling_electrons_2cm,0.,rms_electrons_2cm);
	gr_efficiency_electrons_2cm->Draw("sameAC*");
	
	efficiency_clusters_2cm->cd();
	gr_efficiency_clusters_2cm->SetPoint(counter_filling_clusters_2cm,alpha,mean_clusters_2cm);
	gr_efficiency_clusters_2cm->SetPointError(counter_filling_clusters_2cm,0.,rms_clusters_2cm);
	counter_filling_clusters_2cm = counter_filling_clusters_2cm +1;
	counter_filling_electrons_2cm = counter_filling_electrons_2cm +1;
	gr_efficiency_clusters_2cm->Draw("sameAC*");
	
	if(name_file =="histosTB_run_99.root"){
	  efficiency_clusters_2cm->SaveAs("/lustrehome/bdanzi/offline_analysis/testbeam_analysis/Plots/efficiency_clusters_2cm.pdf",fname.Data());
	  efficiency_electrons_2cm->SaveAs("/lustrehome/bdanzi/offline_analysis/testbeam_analysis/Plots/efficiency_electrons_2cm.pdf",fname.Data()); 
	  integral_summary_2cm->SaveAs("/lustrehome/bdanzi/offline_analysis/testbeam_analysis/Plots/integral_summary_2cm.pdf",fname.Data());
	  maximum_summary_2cm->SaveAs("/lustrehome/bdanzi/offline_analysis/testbeam_analysis/Plots/maximum_summary_2cm.pdf",fname.Data());
	  epc_summary_2cm->SaveAs("/lustrehome/bdanzi/offline_analysis/testbeam_analysis/Plots/epc_summary_2cm.pdf",fname.Data());
	  bsl_summary_2cm->SaveAs("/lustrehome/bdanzi/offline_analysis/testbeam_analysis/Plots/bsl_summary_2cm.pdf",fname.Data());
	  rms_summary_2cm->SaveAs("/lustrehome/bdanzi/offline_analysis/testbeam_analysis/Plots/rms_summary_2cm.pdf",fname.Data());
	  aveph_summary_2cm->SaveAs("/lustrehome/bdanzi/offline_analysis/testbeam_analysis/Plots/aveph_summary_2cm.pdf",fname.Data());
	  
	  efficiency_clusters_2cm->SaveAs("/lustrehome/bdanzi/offline_analysis/testbeam_analysis/Plots/efficiency_clusters_2cm.png",fname.Data());
      efficiency_electrons_2cm->SaveAs("/lustrehome/bdanzi/offline_analysis/testbeam_analysis/Plots/efficiency_electrons_2cm.png",fname.Data()); 
      integral_summary_2cm->SaveAs("/lustrehome/bdanzi/offline_analysis/testbeam_analysis/Plots/integral_summary_2cm.png",fname.Data());
      maximum_summary_2cm->SaveAs("/lustrehome/bdanzi/offline_analysis/testbeam_analysis/Plots/maximum_summary_2cm.png",fname.Data());
      epc_summary_2cm->SaveAs("/lustrehome/bdanzi/offline_analysis/testbeam_analysis/Plots/epc_summary_2cm.png",fname.Data());
      bsl_summary_2cm->SaveAs("/lustrehome/bdanzi/offline_analysis/testbeam_analysis/Plots/bsl_summary_2cm.png",fname.Data());
      rms_summary_2cm->SaveAs("/lustrehome/bdanzi/offline_analysis/testbeam_analysis/Plots/rms_summary_2cm.png",fname.Data());
      aveph_summary_2cm->SaveAs("/lustrehome/bdanzi/offline_analysis/testbeam_analysis/Plots/aveph_summary_2cm.png",fname.Data());
      


	}
      }
      
      
      
      
      
      if(savePlots){
	timediff_clust_2cm->SaveAs(Form("/lustrehome/bdanzi/offline_analysis/testbeam_analysis/Plots/%s/clust_difference_2cm.pdf",fname.Data()));
	timediff_2cm->SaveAs(Form("/lustrehome/bdanzi/offline_analysis/testbeam_analysis/Plots/%s/electrons_difference_2cm.pdf",fname.Data()));
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
	
	timediff_clust_2cm->SaveAs(Form("/lustrehome/bdanzi/offline_analysis/testbeam_analysis/Plots/%s/clust_difference_2cm.png",fname.Data()));
    timediff_2cm->SaveAs(Form("/lustrehome/bdanzi/offline_analysis/testbeam_analysis/Plots/%s/electrons_difference_2cm.png",fname.Data()));
    npeaks_2cm->SaveAs(Form("/lustrehome/bdanzi/offline_analysis/testbeam_analysis/Plots/%s/npeaks_2cm.png",fname.Data()));
    tpeaks_2cm->SaveAs(Form("/lustrehome/bdanzi/offline_analysis/testbeam_analysis/Plots/%s/tpeaks_2cm.png",fname.Data()));
    tfpeaks_2cm->SaveAs(Form("/lustrehome/bdanzi/offline_analysis/testbeam_analysis/Plots/%s/tfpeaks_2cm.png",fname.Data()));
    tlpeaks_2cm->SaveAs(Form("/lustrehome/bdanzi/offline_analysis/testbeam_analysis/Plots/%s/tlpeaks_2cm.png",fname.Data()));
    hnpeaks_2cm->SaveAs(Form("/lustrehome/bdanzi/offline_analysis/testbeam_analysis/Plots/%s/hnpeaks_2cm.png",fname.Data()));
    npeaks_clustser_2cm->SaveAs(Form("/lustrehome/bdanzi/offline_analysis/testbeam_analysis/Plots/%s/hnpeaks_cluster_2cm.png",fname.Data()));
    hpeaks_2cm->SaveAs(Form("/lustrehome/bdanzi/offline_analysis/testbeam_analysis/Plots/%s/hpeaks_2cm.png",fname.Data()));
    bsl_2cm->SaveAs(Form("/lustrehome/bdanzi/offline_analysis/testbeam_analysis/Plots/%s/bsl_2cm.png",fname.Data()));
    max_2cm->SaveAs(Form("/lustrehome/bdanzi/offline_analysis/testbeam_analysis/Plots/%s/max_2cm.png",fname.Data()));
    integ_2cm->SaveAs(Form("/lustrehome/bdanzi/offline_analysis/testbeam_analysis/Plots/%s/integ_2cm.png",fname.Data()));
    rms_2cm->SaveAs(Form("/lustrehome/bdanzi/offline_analysis/testbeam_analysis/Plots/%s/rms_2cm.png",fname.Data()));
    cluster_population_canvas_2cm->SaveAs(Form("/lustrehome/bdanzi/offline_analysis/testbeam_analysis/Plots/%s/cluster_population_2cm.png",fname.Data()));
    //min->SaveAs(Form("/lustrehome/bdanzi/offline_analysis/testbeam_analysis/Plots/%s/min_ch%i.png",fname.Data(),i));
    


      }
      
    }
    
    
 
    
    
    fclose(fp);
  }
  
  
}
