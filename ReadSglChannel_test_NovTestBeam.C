/**  
 *
 *  Authors: B. D'Anzi - University and INFN Bari
 * 			F. Cuna - University and INFN Lecce
 *
 **/

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

void ReadSglChannel_test_NovTestBeam(){
  
  gStyle->SetOptFit(1111);
  // For the axis titles:
  gStyle->SetTitleColor(1, "XYZ");
  gStyle->SetTitleFont(42, "XYZ");
  //gStyle->SetTitleSize(0.02, "XYZ");
  gStyle->SetTitleXSize(0.03); // Another way to set the size?
  gStyle->SetTitleYSize(0.03);
  gStyle->SetTitleXOffset(1.6);//0.9);
  gStyle->SetTitleYOffset(1.8); // => 1.15 if exponents
  // theStyle->SetTitleOffset(1.1, "Y"); // Another way to set the Offset  
  // For the axis labels:
  gStyle->SetLabelColor(1, "XYZ");
  gStyle->SetLabelFont(42, "XYZ");
  gStyle->SetLabelOffset(0.007, "XYZ");
  gStyle->SetLabelSize(0.03, "XYZ");  
  // For the axis:
  gStyle->SetAxisColor(1, "XYZ");
  gStyle->SetStripDecimals(kTRUE);
  gStyle->SetTickLength(0.03, "XYZ");
  gStyle->SetNdivisions(510, "XYZ");
  gStyle->SetPadTickX(1);  // To get tick marks on the opposite side of the
  // For the Global title:
  gStyle->SetTitleFont(42);
  gStyle->SetTitleColor(1);
  gStyle->SetTitleTextColor(1);
  gStyle->SetTitleFillColor(10);
  gStyle->SetTitleFontSize(0.05);
  gStyle->SetHistLineColor(kBlack);
  gStyle->SetHistLineStyle(0);
  gStyle->SetHistLineWidth(1);
  // For the frame:
  gStyle->SetFrameBorderMode(0);
  gStyle->SetFrameBorderSize(1);
  gStyle->SetFrameFillColor(0);
  gStyle->SetFrameFillStyle(0);
  gStyle->SetFrameLineColor(1);
  gStyle->SetFrameLineStyle(1);
  gStyle->SetFrameLineWidth(1);
  
  // For the statistics box:
  gStyle->SetOptFile(0);
  //theStyle->SetOptStat(0); // To display the mean and RMS:   SetOptStat("mr");
  gStyle->SetOptStat("emr");
  gStyle->SetStatColor(kWhite);
  gStyle->SetStatFont(42);
  // theStyle->SetStatFontSize(0.05);
  gStyle->SetStatTextColor(1);
  gStyle->SetStatFormat("6.4g");
  gStyle->SetStatBorderSize(1);
  gStyle->SetPadBorderMode(0);
  gStyle->SetFrameBorderMode(0);
  gStyle->SetPadBottomMargin(0.12);
  gStyle->SetPadLeftMargin(0.14);
  gStyle->cd();  
  // For the canvas:
  gStyle->SetCanvasBorderMode(0);
  gStyle->SetCanvasColor(kWhite);
  gStyle->SetCanvasDefH(600); //Height of canvas
  gStyle->SetCanvasDefW(900); //Width of canvas
  gStyle->SetCanvasDefX(0);   //POsition on screen
  gStyle->SetCanvasDefY(0);
  
  int channel=0;
  int ev=0;
  struct stat st = {0};
  char filename[100];
  TString fname("");
  bool isInteractive = false;
  bool isdoubleCanvas = false;
  bool isLogarithm = true;
  bool isNov2021TestBeam = false;
  bool isJuly2022TestBeam = true;
  bool isScanAngleStudy = false;
  bool isScanHVStudy = false;
  bool isGasMixtureStudy = false;
  bool isScanSamplingStudy = true;
  
  //   if(isNov2021TestBeam){tubes.Data() = "2";}
  //else if(isJuly2022TestBeam)
  //{
  TString tubes("1.5");
  //}
  
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
  Int_t number_of_trials = 0;
  Int_t counter_filling_electrons_1cm = 0;
  Int_t counter_filling_clusters_1cm = 0;
  Int_t counter_filling_electrons_2cm = 0;
  Int_t counter_filling_clusters_2cm = 0;
  Int_t counter_plots_1cm = 1;
  Int_t counter_plots_2cm = 1;
  Int_t counter_plots_1p5cm = 1;
  // if(isNov2021TestBeam){
  // int Channel_1cm[] = {4,5,6,7,8,9}; // Nov 2021 Test Beam
  // int Channel_2cm[] = {10,11,12}; // Nov 2021 Test Beam
  // int Channel_1p5cm[] = {0}; // NO Nov 2021 Test Beam
  // }
  // else if(isJuly2022TestBeam){
  int Channel_1cm[] = {1,2,3,5,6,8,9,10}; // July 2022 Test Beam
  int Channel_2cm[] = {15}; // NO in July 2022 Test Beam
  int Channel_1p5cm[] = {0,4,7,11}; // New Test Beam Channel 0 is weird, I would skip it!
  // }
  int Channel_1cm_size = sizeof Channel_1cm / sizeof Channel_1cm[0];
  int Channel_2cm_size = sizeof Channel_2cm / sizeof Channel_2cm[0];
  int Channel_1p5cm_size = sizeof Channel_1p5cm / sizeof Channel_1p5cm[0];
  int isChannel_1cm = 0;
  int isChannel_2cm = 0;
  int isChannel_1p5cm = 0;
  int number_of_graphs_summary = 0;
  
  if (isScanAngleStudy && isNov2021TestBeam){ number_of_graphs_summary = 5;}
  if ((isGasMixtureStudy) && isNov2021TestBeam){ number_of_graphs_summary = 3;}
  if ((isScanSamplingStudy) && isNov2021TestBeam){ number_of_graphs_summary = 2;}
  if ((isScanHVStudy) && isNov2021TestBeam){ number_of_graphs_summary = 2;}
  if ((isScanAngleStudy || isScanHVStudy) && isJuly2022TestBeam){ number_of_graphs_summary = 4;}
  if ((isGasMixtureStudy) && isJuly2022TestBeam){ number_of_graphs_summary = 3;}
  if ((isScanSamplingStudy) && isJuly2022TestBeam){ number_of_graphs_summary = 2;}
  TCanvas *aveph_summary_1cm = new TCanvas("aveph_summary_1cm","Aveph 1 cm cell size Drift Tubes",600,600);
  TCanvas *aveph_summary_2cm = new TCanvas("aveph_summary_2cm",Form("Aveph %s cm cell size Drift Tubes",tubes.Data()),600,600);
  TGraphErrors* gr_aveph_summary_1cm = new TGraphErrors(number_of_graphs_summary);
  TGraphErrors* gr_aveph_summary_2cm = new TGraphErrors(number_of_graphs_summary);
  gr_aveph_summary_2cm->SetTitle(Form("Aveph %s cm cell size Drift Tubes",tubes.Data()));
  gr_aveph_summary_2cm->SetMarkerColor(kBlue);
  gr_aveph_summary_2cm->SetLineColor(kBlue);
  if(isScanAngleStudy){
    gr_aveph_summary_2cm->GetXaxis()->SetTitle("Track angle (deg)");
  }
  else if(isScanHVStudy){
    gr_aveph_summary_2cm->GetXaxis()->SetTitle("HV Configuration");
  }
  else if(isScanSamplingStudy){
    gr_aveph_summary_2cm->GetXaxis()->SetTitle("Sampling Rate (GSa/s)");
  }
  else if(isGasMixtureStudy){
    gr_aveph_summary_2cm->GetXaxis()->SetTitle("Gas Mixtures");
  }
  gr_aveph_summary_2cm->GetYaxis()->SetTitle("Average Pulse height (V)");
  
  gr_aveph_summary_1cm->SetTitle("Aveph 1 cm cell size Drift Tubes");
  gr_aveph_summary_1cm->GetYaxis()->SetTitle("Average Pulse Height (V)");
  gr_aveph_summary_1cm->SetMarkerColor(kBlue);
  gr_aveph_summary_1cm->SetLineColor(kBlue);
  if(isScanAngleStudy){
    gr_aveph_summary_1cm->GetXaxis()->SetTitle("Track angle (deg)");
  }
  else if(isScanHVStudy){
    gr_aveph_summary_1cm->GetXaxis()->SetTitle("HV Configuration");
  }
  else if(isScanSamplingStudy){
    gr_aveph_summary_1cm->GetXaxis()->SetTitle("Sampling Rate (GSa/s)");
  }
  else if(isGasMixtureStudy){
    gr_aveph_summary_1cm->GetXaxis()->SetTitle("Gas Mixtures");
  }
  
  TCanvas *rms_summary_1cm = new TCanvas("rms_summary_1cm","Rms 1 cm cell size Drift Tubes",600,600);
  TCanvas *rms_summary_2cm = new TCanvas("rms_summary_2cm",Form("Rms %s cm cell size Drift Tubes",tubes.Data()),600,600);
  TGraphErrors* gr_rms_summary_1cm = new TGraphErrors(number_of_graphs_summary);
  TGraphErrors* gr_rms_summary_2cm = new TGraphErrors(number_of_graphs_summary);
  gr_rms_summary_2cm->SetTitle(Form("Rms %s cm cell size Drift Tubes",tubes.Data()));
  gr_rms_summary_2cm->SetMarkerColor(kBlue);
  gr_rms_summary_2cm->SetLineColor(kBlue);
  if(isScanAngleStudy){
    gr_rms_summary_2cm->GetXaxis()->SetTitle("Track angle (deg)");
  }
  else if(isScanHVStudy){
    gr_rms_summary_2cm->GetXaxis()->SetTitle("HV Configuration");
  }
  else if(isScanSamplingStudy){
    gr_rms_summary_2cm->GetXaxis()->SetTitle("Sampling Rate (GSa/s)");
  }
  else if(isGasMixtureStudy){
    gr_rms_summary_2cm->GetXaxis()->SetTitle("Gas Mixtures");
  }
  gr_rms_summary_2cm->GetYaxis()->SetTitle("Rms (mV)");
  gr_rms_summary_1cm->SetTitle("Rms 1 cm cell size Drift Tubes");
  gr_rms_summary_1cm->SetMarkerColor(kBlue);
  gr_rms_summary_1cm->SetLineColor(kBlue);
  gr_rms_summary_1cm->GetYaxis()->SetTitle("Rms (mV)");
  if(isScanAngleStudy){
    gr_rms_summary_1cm->GetXaxis()->SetTitle("Track angle (deg)");
  }
  else if(isScanHVStudy){
    gr_rms_summary_1cm->GetXaxis()->SetTitle("HV Configuration");
  }
  else if(isScanSamplingStudy){
  gr_rms_summary_1cm->GetXaxis()->SetTitle("Sampling Rate (GSa/s)");
  }
  else if(isGasMixtureStudy){
    gr_rms_summary_1cm->GetXaxis()->SetTitle("Gas Mixtures");
  }
  
  TCanvas *bsl_summary_1cm = new TCanvas("bsl_summary_1cm","Bsl 1 cm cell size Drift Tubes",600,600);
  TCanvas *bsl_summary_2cm = new TCanvas("bsl_summary_2cm",Form("Bsl %s cm cell size Drift Tubes",tubes.Data()),600,600);
  TGraphErrors* gr_bsl_summary_1cm = new TGraphErrors(number_of_graphs_summary);
  TGraphErrors* gr_bsl_summary_2cm = new TGraphErrors(number_of_graphs_summary);
  gr_bsl_summary_2cm->SetTitle(Form("Bsl %s cm cell size Drift Tubes",tubes.Data()));
  gr_bsl_summary_2cm->SetMarkerColor(kBlue);
  gr_bsl_summary_2cm->SetLineColor(kBlue);
  if(isScanAngleStudy){
    gr_bsl_summary_2cm->GetXaxis()->SetTitle("Track angle (deg)");
  }
  else if(isScanHVStudy){
    gr_bsl_summary_2cm->GetXaxis()->SetTitle("HV Configuration");
  }
  else if(isScanSamplingStudy){
    gr_bsl_summary_2cm->GetXaxis()->SetTitle("Sampling Rate (GSa/s)");
  }
  else if(isGasMixtureStudy){
    gr_bsl_summary_2cm->GetXaxis()->SetTitle("Gas Mixtures");
  }
  gr_bsl_summary_2cm->GetYaxis()->SetTitle("Baseline (V)");
  gr_bsl_summary_1cm->SetTitle("Bsl 1 cm cell size Drift Tubes");
  gr_bsl_summary_1cm->SetMarkerColor(kBlue);
  gr_bsl_summary_1cm->SetLineColor(kBlue);
  gr_bsl_summary_1cm->GetYaxis()->SetTitle("Baseline (V)");
  if(isScanAngleStudy){
    gr_bsl_summary_1cm->GetXaxis()->SetTitle("Track angle (deg)");
  }
  else if(isScanHVStudy){
    gr_bsl_summary_1cm->GetXaxis()->SetTitle("HV Configuration");
  }
  else if(isScanSamplingStudy){
    gr_bsl_summary_1cm->GetXaxis()->SetTitle("Sampling Rate (GSa/s)");
  }
  else if(isGasMixtureStudy){
    gr_bsl_summary_1cm->GetXaxis()->SetTitle("Gas Mixtures");
  }
  
  
  TCanvas *epc_summary_1cm = new TCanvas("epc_summary_1cm","Epc 1 cm cell size Drift Tubes",600,600);
  TCanvas *epc_summary_2cm = new TCanvas("epc_summary_2cm",Form("Epc %s cm cell size Drift Tubes",tubes.Data()),600,600);
  TGraphErrors* gr_epc_summary_1cm = new TGraphErrors(number_of_graphs_summary);
  TGraphErrors* gr_epc_summary_2cm = new TGraphErrors(number_of_graphs_summary);
  TLegend *leg_epc_summary_1cm= new TLegend(0.5,0.75,0.85,0.85); 
  gr_epc_summary_2cm->SetTitle(Form("Epc %s cm cell size Drift Tubes",tubes.Data()));
  gr_epc_summary_2cm->SetMarkerColor(kBlue);
  gr_epc_summary_2cm->SetLineColor(kBlue);
  if(isScanAngleStudy){
    gr_epc_summary_2cm->GetXaxis()->SetTitle("Track angle (deg)");
  }
  else if(isScanHVStudy){
    gr_epc_summary_2cm->GetXaxis()->SetTitle("HV Configuration");
  }
  else if(isScanSamplingStudy){
    gr_epc_summary_2cm->GetXaxis()->SetTitle("Sampling Rate (GSa/s)");
  }
  else if(isGasMixtureStudy){
    gr_epc_summary_2cm->GetXaxis()->SetTitle("Gas Mixtures");
  }
  gr_epc_summary_2cm->GetYaxis()->SetTitle("Electrons per Cluster");
  gr_epc_summary_1cm->SetTitle("Epc 1 cm cell size Drift Tubes");
  gr_epc_summary_1cm->SetMarkerColor(kBlue);
  gr_epc_summary_1cm->SetLineColor(kBlue);
  gr_epc_summary_1cm->GetYaxis()->SetTitle("Electrons per Cluster");
  if(isScanAngleStudy){
    gr_epc_summary_1cm->GetXaxis()->SetTitle("Track angle (deg)");
  }
  else if(isScanHVStudy){
    gr_epc_summary_1cm->GetXaxis()->SetTitle("HV Configuration");
  }
  else if(isScanSamplingStudy){
    gr_epc_summary_1cm->GetXaxis()->SetTitle("Sampling Rate (GSa/s)");
  }
  else if(isGasMixtureStudy){
    gr_epc_summary_1cm->GetXaxis()->SetTitle("Gas Mixtures");
  }
  
  TCanvas *integral_summary_1cm = new TCanvas("integral_summary_1cm","Integral 1 cm cell size Drift Tubes",600,600);
  TCanvas *integral_summary_2cm = new TCanvas("integral_summary_2cm",Form("Integral %s cm cell size Drift Tubes",tubes.Data()),600,600);
  TGraphErrors* gr_integral_summary_1cm = new TGraphErrors(number_of_graphs_summary);
  TGraphErrors* gr_integral_summary_2cm = new TGraphErrors(number_of_graphs_summary);
  gr_integral_summary_2cm->SetTitle(Form("Integral %s cm cell size Drift Tubes",tubes.Data()));
  gr_integral_summary_2cm->SetMarkerColor(kBlue);
  gr_integral_summary_2cm->SetLineColor(kBlue);
  gr_integral_summary_2cm->GetYaxis()->SetTitle("Charge (pC)");
  if(isScanAngleStudy){
    gr_integral_summary_2cm->GetXaxis()->SetTitle("Track angle (deg)");
  }
  else if(isScanHVStudy){
    gr_integral_summary_2cm->GetXaxis()->SetTitle("HV Configuration");
  }
  else if(isScanSamplingStudy){
    gr_integral_summary_2cm->GetXaxis()->SetTitle("Sampling Rate (GSa/s)");
  }
  else if(isGasMixtureStudy){
    gr_integral_summary_2cm->GetXaxis()->SetTitle("Gas Mixtures");
  }
  gr_integral_summary_1cm->SetTitle("Integral 1 cm cell size Drift Tubes");
  gr_integral_summary_1cm->SetMarkerColor(kBlue);
  gr_integral_summary_1cm->SetLineColor(kBlue);
  if(isScanAngleStudy){
    gr_integral_summary_1cm->GetXaxis()->SetTitle("Track angle (deg)");
  }
  else if(isScanHVStudy){
    gr_integral_summary_1cm->GetXaxis()->SetTitle("HV Configuration");
  }
  else if(isScanSamplingStudy){
    gr_integral_summary_1cm->GetXaxis()->SetTitle("Sampling Rate (GSa/s)");
  }
  else if(isGasMixtureStudy){
    gr_integral_summary_1cm->GetXaxis()->SetTitle("Gas Mixtures");
  }
  gr_integral_summary_1cm->GetYaxis()->SetTitle("Charge (pC)");
  
  TCanvas *maximum_summary_1cm = new TCanvas("maximum_summary_1cm","Maximum 1 cm cell size Drift Tubes",600,600);
  TCanvas *maximum_summary_2cm = new TCanvas("maximum_summary_2cm",Form("Maximum %s cm cell size Drift Tubes",tubes.Data()),600,600);
  TGraphErrors* gr_maximum_summary_1cm = new TGraphErrors(number_of_graphs_summary);
  TGraphErrors* gr_maximum_summary_2cm = new TGraphErrors(number_of_graphs_summary);
  gr_maximum_summary_2cm->SetTitle(Form("Maximum %s cm cell size Drift Tubes",tubes.Data()));
  gr_maximum_summary_2cm->SetMarkerColor(kBlue);
  gr_maximum_summary_2cm->SetLineColor(kBlue);
  gr_maximum_summary_2cm->GetYaxis()->SetTitle("Voltage (V)");
  if(isScanAngleStudy){
    gr_maximum_summary_2cm->GetXaxis()->SetTitle("Track angle (deg)");
  }
  else if(isScanHVStudy){
    gr_maximum_summary_2cm->GetXaxis()->SetTitle("HV Configuration");
  }
  else if(isScanSamplingStudy){
    gr_maximum_summary_2cm->GetXaxis()->SetTitle("Sampling Rate (GSa/s)");
  }
  else if(isGasMixtureStudy){
    gr_maximum_summary_2cm->GetXaxis()->SetTitle("Gas Mixtures");
  }
  gr_maximum_summary_1cm->SetTitle("Maximum 1 cm cell size Drift Tubes");
  gr_maximum_summary_1cm->SetMarkerColor(kBlue);
  gr_maximum_summary_1cm->SetLineColor(kBlue);
  if(isScanAngleStudy){
    gr_maximum_summary_1cm->GetXaxis()->SetTitle("Track angle (deg)");
  }
  else if(isScanHVStudy){
    gr_maximum_summary_1cm->GetXaxis()->SetTitle("HV Configuration");
  }
  else if(isScanSamplingStudy){
    gr_maximum_summary_1cm->GetXaxis()->SetTitle("Sampling Rate (GSa/s)");
  }
  else if(isGasMixtureStudy){
    gr_maximum_summary_1cm->GetXaxis()->SetTitle("Gas Mixtures");
  }
  gr_maximum_summary_1cm->GetYaxis()->SetTitle("Voltage (V)");
  
  
  TCanvas *efficiency_electrons_1cm = new TCanvas("efficiency_electrons_1cm","Electron Finding Efficiency 1 cm cell size Drift Tubes",600,600);
  TCanvas *efficiency_electrons_2cm = new TCanvas("efficiency_electrons_2cm",Form("Electron Finding Efficiency %s cm cell size Drift Tubes",tubes.Data()),600,600);
  TGraphErrors* gr_efficiency_electrons_1cm = new TGraphErrors(number_of_graphs_summary);
  TGraphErrors* gr_efficiency_electrons_2cm = new TGraphErrors(number_of_graphs_summary);
  
  TCanvas *efficiency_clusters_1cm = new TCanvas("efficiency_clusters_1cm","Cluster Finding Efficiency 1 cm cell size Drift Tubes",600,600);
  TCanvas *efficiency_clusters_2cm = new TCanvas("efficiency_clusters_2cm",Form("Cluster Finding Efficiency %s cm cell size Drift Tubes",tubes.Data()),600,600);
  TGraphErrors* gr_efficiency_clusters_2cm = new TGraphErrors(number_of_graphs_summary);
  TGraphErrors* gr_efficiency_clusters_1cm = new TGraphErrors(number_of_graphs_summary);
  
  gr_efficiency_electrons_2cm->SetTitle(Form("Electrons Finding Efficiency %s cm cell size Drift Tubes",tubes.Data()));
  gr_efficiency_clusters_2cm->SetTitle(Form("Clusters Finding Efficiency %s cm cell size Drift Tubes",tubes.Data()));
  gr_efficiency_electrons_1cm->SetTitle("Electrons Finding Efficiency 1 cm cell size Drift Tubes");
  gr_efficiency_clusters_1cm->SetTitle("Clusters Finding Efficiency 1 cm cell size Drift Tubes");
  
  gr_efficiency_electrons_2cm->SetMarkerColor(kBlue);
  gr_efficiency_clusters_2cm->SetMarkerColor(kBlue);
  gr_efficiency_electrons_2cm->SetLineColor(kBlue);
  gr_efficiency_clusters_2cm->SetLineColor(kBlue);
  
  if(isScanAngleStudy){
    gr_efficiency_electrons_2cm->GetXaxis()->SetTitle("Track Angle (deg)");
    gr_efficiency_clusters_2cm->GetXaxis()->SetTitle("Track Angle (deg)");
    gr_efficiency_electrons_1cm->GetXaxis()->SetTitle("Track Angle (deg)");
    gr_efficiency_clusters_1cm->GetXaxis()->SetTitle("Track Angle (deg)");
  }
  else if(isScanHVStudy){
    gr_efficiency_electrons_2cm->GetXaxis()->SetTitle("HV Configuration");
    gr_efficiency_clusters_2cm->GetXaxis()->SetTitle("HV Configuration");
    gr_efficiency_electrons_1cm->GetXaxis()->SetTitle("HV Configuration");
    gr_efficiency_clusters_1cm->GetXaxis()->SetTitle("HV Configuration");
  }
  else if(isScanSamplingStudy){
    gr_efficiency_electrons_2cm->GetXaxis()->SetTitle("Sampling Rate (GSa/s)");
    gr_efficiency_clusters_2cm->GetXaxis()->SetTitle("Sampling Rate (GSa/s)");
    gr_efficiency_electrons_1cm->GetXaxis()->SetTitle("Sampling Rate (GSa/s)");
    gr_efficiency_clusters_1cm->GetXaxis()->SetTitle("Sampling Rate (GSa/s)");
  }
  else if(isGasMixtureStudy){
    gr_efficiency_electrons_2cm->GetXaxis()->SetTitle("Gas Mixtures");
    gr_efficiency_clusters_2cm->GetXaxis()->SetTitle("Gas Mixtures");
    gr_efficiency_electrons_1cm->GetXaxis()->SetTitle("Gas Mixtures");
    gr_efficiency_clusters_1cm->GetXaxis()->SetTitle("Gas Mixtures");
  }
  gr_efficiency_electrons_2cm->GetYaxis()->SetTitle("Measured Average Number of Electrons / Expected Number of Electrons");
  gr_efficiency_clusters_2cm->GetYaxis()->SetTitle("Measured Average Number of Clusters / Expected Number of Clusters");
  gr_efficiency_electrons_1cm->GetYaxis()->SetTitle("Measured Average Number of Electrons / Expected Number of Electrons");
  gr_efficiency_clusters_1cm->GetYaxis()->SetTitle("Measured Average Number of Clusters / Expected Number of Clusters");
  
  
  if(isInteractive){ // Deprecated option, to be updated
    
    printf("Enter *.root name:");
    scanf("%s",fname.Data());
    printf("Enter the number of channels to be shown:");
    scanf("%d",&channel);
    printf("Enter the number of events (waveForms) to be shown:");
    scanf("%d",&ev);
    if (stat(Form("/lustrehome/bdanzi/TestBeam2022/analysistestbeam2022/analysis_2021/drifttubes_offline_analysis/testbeam_analysis/Plots/%s/",fname.Data()), &st) == -1) {
      mkdir(Form("/lustrehome/bdanzi/TestBeam2022/analysistestbeam2022/analysis_2021/drifttubes_offline_analysis/testbeam_analysis/Plots/%s/",fname.Data()), 0700);
    }
    if (stat(Form("/lustrehome/bdanzi/TestBeam2022/analysistestbeam2022/analysis_2021/drifttubes_offline_analysis/testbeam_analysis/Plots/%s/Waves/",fname.Data()), &st) == -1) {
      mkdir(Form("/lustrehome/bdanzi/TestBeam2022/analysistestbeam2022/analysis_2021/drifttubes_offline_analysis/testbeam_analysis/Plots/%s/Waves/",fname.Data()), 0700);
    }
  }
  else{
    char inputfile[100] ="plots_oldTestBeam.txt";
    FILE *fp = fopen(inputfile, "r");
    if (fp == NULL)
      {
        printf("Error: could not open file %s", inputfile);
       	return 1;
      }
    
    while ((fscanf(fp, "%s", fname.Data() ) ) != -1){
      
      counter_plots_1cm = 1;
      counter_plots_2cm = 1;
      counter_plots_1p5cm = 1;
      printf("Opening the file %s\n",fname.Data());
      string name_file = fname.Data();
      float _gsample = 0.0;
      TString run("");
      TString N1(""), N2(""), N3(""), N4(""), scale_cut(""),gsample("");
      fscanf(fp,"%s",N1.Data());
      fscanf(fp,"%s",N2.Data());
      fscanf(fp,"%s",N3.Data());
      fscanf(fp,"%s",N4.Data());
      fscanf(fp,"%s",scale_cut.Data());
      fscanf(fp,"%s",gsample.Data());
      _gsample = atof(gsample.Data());
      fscanf(fp,"%s",run.Data());
      printf("Opening the file with the reduced name %s, GSa/s %0.1f \n",run.Data(), _gsample);
      string name_file_compact = run.Data();
      fscanf(fp,"%d",&channel); //number of channels for which I show plots (4-14)
      fscanf(fp,"%d",&ev); //number of events for which I show plots (0-...)
      //   if(isNov2021TestBeam){
      // string Runs_80_20[] = {"histosTB_run_117.root", "histosTB_run_127.root"}; // 2021 Nov Test Beam
      // //string Runs_90_10[] = {"histosTB_run_72.root", "histosTB_run_73.root", "histosTB_run_74.root", "histosTB_run_86.root", "histosTB_run_87.root", "histosTB_run_88.root", "histosTB_run_89.root", "histosTB_run_90.root", "histosTB_run_91.root", "histosTB_run_92.root", "histosTB_run_93.root", "histosTB_run_94.root", "histosTB_run_95.root", "histosTB_run_96.root", "histosTB_run_97.root" ,"histosTB_run_98.root", "histosTB_run_99.root", "histosTB_run_100.root", "histosTB_run_101.root"};
      // // Scan Angle Study
      // string Runs_90_10[] = {"histosTB_run_99.root","histosTB_run_98.root", "histosTB_run_96.root", "histosTB_run_94.root","histosTB_run_91.root"};
      // string Runs_85_15[] = {"histosTB_run_0.root"}; // No data for Nov 2021, 2022 Test Beam to be inserted
      // string Runs_alpha_0[] = {"histosTB_run_99.root", "histosTB_run_117.root", "histosTB_run_86.root", "histosTB_run_100.root", "histosTB_run_72.root", "histosTB_run_73.root", "histosTB_run_74.root"}; // 2021 Test Beam
      // string Runs_alpha_15[] = {"histosTB_run_98.root", "histosTB_run_87.root", "histosTB_run_97.root"}; // 2021 Nov Test Beam
      // string Runs_alpha_30[] = {"histosTB_run_96.root", "histosTB_run_88.root", "histosTB_run_95.root"}; // 2021 Nov Test Beam
      // string Runs_alpha_45[] = {"histosTB_run_94.root", "histosTB_run_89.root", "histosTB_run_93.root"}; // 2021 Nov Test Beam
      // string Runs_alpha_60[] = {"histosTB_run_91.root", "histosTB_run_127.root", "histosTB_run_90.root" ,"histosTB_run_92.root"}; // 2021 Nov Test Beam
      //   }
      //   else if(isJuly2022TestBeam){
      //string Runs_80_20[] = {"histosTB_run_39.root","histosTB_run_40.root", "histosTB_run_41.root", "histosTB_run_42.root","histosTB_run_43.root","histosTB_run_44.root","histosTB_run_45.root","histosTB_run_46.root"}; // 2022 July Test Beam
      ////Scan Study: HV
      //string Runs_80_20[] = {"histosTB_run_40.root","histosTB_run_41.root","histosTB_run_42.root","histosTB_run_43.root"};
      //// Scan Study: Angle
      //string Runs_80_20[] = {"histosTB_run_44.root", "histosTB_run_45.root","histosTB_run_41.root","histosTB_run_46.root"};
      ////Scan Study: Gas mixture
      string Runs_80_20[] = {"histosTB_run_41.root"};
      // Scan Study: Sampling rate
      string Runs_90_10[] = {"histosTB_run_9.root","histosTB_run_10.root"}; // 2022 July Test Beam
      string Runs_85_15[] = {"histosTB_run_63.root"}; // 2022 July Test Beam
      string Runs_alpha_0[] = {"histosTB_run_44.root"}; // 2022 July Test Beam
      string Runs_alpha_15[] = {"histosTB_run_0.root"}; // We didn't take data at 15°, 2022 July Test Beam
      string Runs_alpha_30[] = {"histosTB_run_45.root"}; // 2022 July Test Beam
      string Runs_alpha_45[] = {"histosTB_run_39.root","histosTB_run_40.root", "histosTB_run_41.root", "histosTB_run_42.root","histosTB_run_43.root","histosTB_run_9.root","histosTB_run_10.root","histosTB_run_63.root"}; // 2022 July Test Beam
      string Runs_alpha_60[] = {"histosTB_run_46.root"}; // 2022 July Test Beam
      //   }
      int Runs_80_20_size = sizeof Runs_80_20 / sizeof Runs_80_20[0];
      int Runs_90_10_size = sizeof Runs_90_10 / sizeof Runs_90_10[0];
      int Runs_85_15_size = sizeof Runs_85_15 / sizeof Runs_85_15[0];
      int Runs_alpha_0_size = sizeof Runs_alpha_0 / sizeof Runs_alpha_0[0];
      int Runs_alpha_15_size =sizeof Runs_alpha_15 / sizeof Runs_alpha_15[0];
      int Runs_alpha_30_size = sizeof Runs_alpha_30 / sizeof Runs_alpha_30[0];
      int Runs_alpha_45_size = sizeof Runs_alpha_45 / sizeof Runs_alpha_45[0];
      int Runs_alpha_60_size = sizeof Runs_alpha_60 / sizeof Runs_alpha_60[0];
      int isRuns_80_20 = 0;
      int isRuns_90_10 = 0;
      int isRuns_85_15 = 0;
      int isRuns_alpha_0 = 0;
      int isRuns_alpha_15 = 0;
      int isRuns_alpha_30 = 0;
      int isRuns_alpha_45 = 0;
      int isRuns_alpha_60 = 0;
      TString gasMixture("");
      for (int i = 0; i < Runs_80_20_size; i++) {
	if (Runs_80_20[i] == name_file_compact) {
	  gasMixture = gasMixture.Append("Gas Mixture He:IsoB 80/20");
	  cout << gasMixture.Data() << endl;
	  isRuns_80_20 = 1;
	  cluster_per_cm_mip = 18.;
	  printf("Gas mixture changed to 80/20 with Number of Cluster/cm (MIP) to be %0.1f!\n",cluster_per_cm_mip);
	  break;
	}
      }
      for (int i = 0; i < Runs_90_10_size; i++) {
	if (Runs_90_10[i] == name_file_compact) {
	  gasMixture = gasMixture.Append("Gas Mixture He:IsoB 90/10");
	  cout << gasMixture.Data() << endl;
	  isRuns_90_10 = 1;
	  cluster_per_cm_mip = 12.;
	  printf("Gas mixture changed to 90/10 with Number of Cluster/cm (MIP) to be %0.1f!\n",cluster_per_cm_mip);
	  break;
	}
      }
      for (int i = 0; i < Runs_85_15_size; i++) {
	if (Runs_85_15[i] == name_file_compact) {
	  gasMixture = gasMixture.Append("Gas Mixture He:IsoB 85/15");
	  cout << gasMixture.Data() << endl;
	  isRuns_85_15 = 1;
	  cluster_per_cm_mip = 15.;
	  printf("Gas mixture changed to 85/15 with Number of Cluster/cm (MIP) to be %0.1f!\n",cluster_per_cm_mip);
	  break;
	}
      }
      for (int i = 0; i < Runs_alpha_0_size; i++) {
	if (Runs_alpha_0[i] == name_file_compact) {
	  isRuns_alpha_0 = 1;
	  alpha = 0.;
	  printf("Track Angle Changed to be %0.1f!\n",alpha);
	  break;
	}
      }
      for (int i = 0; i < Runs_alpha_15_size; i++) {
	if (Runs_alpha_15[i] == name_file_compact) {
	  isRuns_alpha_15 = 1;
	  alpha = 15.;
	  printf("Track Angle Changed to be %0.1f!\n",alpha);
	  break;
	}
      }
      for (int i = 0; i < Runs_alpha_30_size; i++) {
	if (Runs_alpha_30[i] == name_file_compact) {
	  isRuns_alpha_30 = 1;
	  alpha = 30.;
	  printf("Track Angle Changed to be %0.1f!\n",alpha);
	  break;
	}
      }
      for (int i = 0; i < Runs_alpha_45_size; i++) {
	if (Runs_alpha_45[i] == name_file_compact) {
	  isRuns_alpha_45 = 1;
	  alpha = 45.;
	  printf("Track Angle Changed to be %0.1f!\n",alpha);
	  break;
	}
      }
      for (int i = 0; i < Runs_alpha_60_size; i++) {
	if (Runs_alpha_60[i] == name_file_compact) {
	  isRuns_alpha_60 = 1;
	  alpha = 60.;
	  printf("Track Angle Changed to be %0.1f!\n",alpha);
	  break;
	}
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
      
      
      if (stat(Form("/lustrehome/bdanzi/TestBeam2022/analysistestbeam2022/analysis_2021/drifttubes_offline_analysis/testbeam_analysis/Plots/%s/",fname.Data()), &st) == -1) {
	mkdir(Form("/lustrehome/bdanzi/TestBeam2022/analysistestbeam2022/analysis_2021/drifttubes_offline_analysis/testbeam_analysis/Plots/%s/",fname.Data()), 0700);
      }
      if (stat(Form("/lustrehome/bdanzi/TestBeam2022/analysistestbeam2022/analysis_2021/drifttubes_offline_analysis/testbeam_analysis/Plots/%s/Waves/",fname.Data()), &st) == -1) {
	mkdir(Form("/lustrehome/bdanzi/TestBeam2022/analysistestbeam2022/analysis_2021/drifttubes_offline_analysis/testbeam_analysis/Plots/%s/Waves/",fname.Data()), 0700);
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
      TCanvas *hncluster_1cm= new TCanvas("hncluster_1","hncluster",3500,1500);
      TCanvas *hncluster_2cm= new TCanvas("hncluster_2","hncluster",3500,1500);
      TCanvas *hnelectron_1cm= new TCanvas("hnelectron_1","hnelectron",3500,1500);
      TCanvas *hnelectron_2cm= new TCanvas("hnelectron_2","hnelectron",3500,1500);
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
      
	 // if Nov2021 Test Beam
      if(isNov2021TestBeam){
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
	hncluster_1cm->Divide(2,3);
	hncluster_2cm->Divide(2,2);
	hnelectron_1cm->Divide(2,3);
	hnelectron_2cm->Divide(2,2);
	tpeaks_1cm->Divide(2,3);
	tpeaks_2cm->Divide(2,2);
	tfpeaks_1cm->Divide(2,3);
	tfpeaks_2cm->Divide(2,2);
	tlpeaks_1cm->Divide(2,3);
	tlpeaks_2cm->Divide(2,2);
      }
      else if(isJuly2022TestBeam){
	cluster_population_canvas_1cm->Divide(2,4);
	cluster_population_canvas_2cm->Divide(2,2);
	timediff_clust_1cm->Divide(2,4);
	timediff_clust_2cm->Divide(2,2);
	timediff_1cm->Divide(2,4);
	timediff_2cm->Divide(2,2);
	npeaks_clustser_1cm->Divide(2,4);
	npeaks_clustser_2cm->Divide(2,2);
	//integ->Divide(4,2);
	//max->Divide(4,2);
	min->Divide(4,1);
	bsl_1cm->Divide(2,4);
	bsl_2cm->Divide(2,2);
	integ_1cm->Divide(2,4);
	integ_2cm->Divide(2,2);
	rms_1cm->Divide(2,4);
	rms_2cm->Divide(2,2);
	max_1cm->Divide(2,4);
	max_2cm->Divide(2,2);
	npeaks_1cm->Divide(2,4);
	npeaks_2cm->Divide(2,2);
	hpeaks_1cm->Divide(2,4);
	hpeaks_2cm->Divide(2,2);
	hnpeaks_1cm->Divide(2,4);
	hnpeaks_2cm->Divide(2,2);
	hncluster_1cm->Divide(2,4);
	hncluster_2cm->Divide(2,2);
	hnelectron_1cm->Divide(2,4);
	hnelectron_2cm->Divide(2,2);
	tpeaks_1cm->Divide(2,4);
	tpeaks_2cm->Divide(2,2);
	tfpeaks_1cm->Divide(2,4);
	tfpeaks_2cm->Divide(2,2);
	tlpeaks_1cm->Divide(2,4);
	tlpeaks_2cm->Divide(2,2);
      }
      if(isdoubleCanvas){
	for(int i =0; i<ev; ++i){
	  TGraph *h1=(TGraph*)file->Get(Form("signal_Afterflt/CvSignal_1_ev%i",i));
	  if (h1==0x0) { continue; }
	  //h1->SaveAs(Form("/lustrehome/bdanzi/TestBeam2022/analysistestbeam2022/analysis_2021/drifttubes_offline_analysis/testbeam_analysis/Plots/%s/Waves/waves_ev%i.pdf",fname.Data(),i));
	  //h1->SaveAs(Form("/lustrehome/bdanzi/offline_analysis/testbeam_analysis/Plots/%s/Waves/waves_ev%i.png",fname.Data(),i));
	}
      }
      else if(!isdoubleCanvas){
	
	for(int i = 0; i<ev; ++i){
	  for(int j = 0; j<=channel; ++j){
	    TGraph *h1=(TGraph*)file->Get(Form("signal_Afterflt/CvSignal_1_Ch%i_ev%i",j,i));
	    if (h1==0x0) { continue; }
	    //h1->SaveAs(Form("/lustrehome/bdanzi/TestBeam2022/analysistestbeam2022/analysis_2021/drifttubes_offline_analysis/testbeam_analysis/Plots/%s/Waves/waves_ev%i_Ch%i.pdf",fname.Data(),i,j));
	    // h1->SaveAs(Form("/lustrehome/bdanzi/offline_analysis/testbeam_analysis/Plots/%s/Waves/waves_ev%i_Ch%i.png",fname.Data(),i,j));
	    
	  }
	}
      }
      
      
      for(int i = 0; i<=channel; ++i){ //looping over all channels
	isChannel_1cm = 0;
	
	for (int j = 0; j < Channel_1cm_size; j++) {
	  if (Channel_1cm[j] == i) {
            isChannel_1cm = 1;
            break;
	  }
	}
	
	if(isChannel_1cm){
	  //2cm/1cm 1,8/0,8 = 2,25
	  drift_size = 0.8;
	  expected_electrons = cluster_per_cm_mip * drift_size*relativistic_rise * cluster_population * 1/cos_alpha;
	  expected_cluster = cluster_per_cm_mip * drift_size*relativistic_rise * 1/cos_alpha;
	  
	  
	  TH1F *h17=(TH1F*)file->Get(Form("H-Ch%i_signal/hMaxVInR_ch%i",i,i));
	  if (h17==0x0) { continue; }
	  if(isNov2021TestBeam) {max_1cm->cd(i-3);}
	  else if(isJuly2022TestBeam) {max_1cm->cd(counter_plots_1cm);}
	  if(isRuns_90_10 || isRuns_80_20 || isRuns_85_15){
	    maximum_1cm = maximum_1cm+ h17->GetMean();
	    rms_maximum_1cm = rms_maximum_1cm + h17->GetRMS();
	  }
	  //max->cd(2);
	  h17->Fit("landau");
	  gPad->SetLogy(0);
	  gPad->SetLogx(0);
	  h17->Draw("same");
	  
	  
	  
	  TH1F *h3=(TH1F*)file->Get(Form("H-Ch%i_signal/hBsl_ch%i",i,i));
	  if (h3==0x0) { continue; }
	  if(isNov2021TestBeam) {bsl_1cm->cd(i-3);}
	  else if(isJuly2022TestBeam) {bsl_1cm->cd(counter_plots_1cm);}
	  if(isRuns_90_10 || isRuns_80_20 || isRuns_85_15){
	    bsl_1cm_var = bsl_1cm_var + h3->GetMean();
	    rms_bsl_1cm = rms_bsl_1cm + h3->GetRMS();
	  }
	  gPad->SetLogy(0);
	  gPad->SetLogx(0);
	  h3->GetXaxis()->SetRangeUser(-0.5,-0.4);
	  h3->Draw("same");
	  
	  
	  TH1F *h4=(TH1F*)file->Get(Form("H-Ch%i_signal/hNPeaks_ch%i",i,i));
	  if (h4==0x0) { continue; }
	  if(isNov2021TestBeam) {npeaks_1cm->cd(i-3);}
	  else if(isJuly2022TestBeam) {npeaks_1cm->cd(counter_plots_1cm);}
	  if(isRuns_90_10 || isRuns_80_20 || isRuns_85_15){
	    mean_electrons_1cm = mean_electrons_1cm + h4->GetMean();
	    rms_electrons_1cm = rms_electrons_1cm + h4->GetRMS();
	  }
	  //h4->Fit("landau");
	  if(isLogarithm){
	    gPad->SetLogy(1);
	  }
	  h4->Draw("same");
	  TPaveText *pt_1cm = new TPaveText(0.1,0.84,0.6,0.89,"NDC");
	  pt_1cm->SetTextSize(0.06);
	  pt_1cm->SetTextColor(kRed);
	  pt_1cm->SetFillColor(kYellow);
	  pt_1cm->SetTextAlign(12);
	  pt_1cm->AddText(Form("Expected Electron Peaks: %.1f - Track angle (deg) %0.1f",expected_electrons, alpha));
	  pt_1cm->Draw("same");
	  TPaveText *pt_1cm_alpha = new TPaveText(0.2,0.84,0.62,0.88,"NDC");
	  pt_1cm_alpha->SetTextSize(0.05);
	  pt_1cm_alpha->SetTextColor(kRed);
	  pt_1cm_alpha->SetFillColor(0);
	  pt_1cm_alpha->SetTextAlign(12);
	  pt_1cm_alpha->AddText(Form("Track angle (deg): %.1f",alpha));
	  //pt_1cm_alpha->Draw("same");
	  
	  TH1F *h20=(TH1F*)file->Get(Form("H-Ch%i_signal/hNPeaks_clust_ch%i",i,i));
	  if (h20==0x0) { continue; }
	  if(isNov2021TestBeam) {npeaks_clustser_1cm->cd(i-3);}
	  else if(isJuly2022TestBeam) {npeaks_clustser_1cm->cd(counter_plots_1cm);}
	  h20->GetXaxis()->SetRangeUser(0.,90.);
	  if(isRuns_90_10 || isRuns_80_20 || isRuns_85_15){
	    mean_clusters_1cm = mean_clusters_1cm + h20->GetMean();
	    rms_clusters_1cm = rms_clusters_1cm + h20->GetRMS();
	  }
	  if(expected_cluster > 0.){
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
	  if(isLogarithm){
	    gPad->SetLogy(1);
	  }
	  h20->Draw("same");
	  TPaveText *pt_1cm_cluster =  new TPaveText(0.2,0.84,0.6,0.88,"NDC");
	  pt_1cm_cluster->SetTextSize(0.07);
	  pt_1cm_cluster->SetTextColor(kRed);
	  pt_1cm_cluster->SetFillColor(kYellow);
	  pt_1cm_cluster->SetTextAlign(12);
	  pt_1cm_cluster->AddText(Form("Expected Clusters: %.1f - Track Angle (deg) %0.1f",expected_cluster,alpha));
	  pt_1cm_cluster->Draw("same");
	  
	  //pt_1cm_alpha->Draw("same");
	  
	  TH1F *h31=(TH1F*)file->Get(Form("H-Ch%i_signal/hNElectrons_per_cluster_ch%i",i,i));
	  if (h31==0x0) { continue; }
	  if(isNov2021TestBeam) {cluster_population_canvas_1cm->cd(i-3);}
	  else if(isJuly2022TestBeam) {cluster_population_canvas_1cm->cd(counter_plots_1cm);}
	  if(isRuns_90_10 || isRuns_80_20 || isRuns_85_15){
	    epc_1cm = epc_1cm+ h31->GetMean();
	    rms_epc_1cm = rms_epc_1cm + h31->GetRMS();
	  }
	  //h31->GetXaxis()->SetRangeUser(0.,5.);
	  //h31->Fit("expo");
	  h31->Draw("same");
	  //   if(isLogarithm){
	  //   	gPad->SetLogy(1);
	  //   }
	  TPaveText *cluster_population_1cm = new TPaveText(0.4,0.85,0.7,0.88,"NDC");
	  cluster_population_1cm->SetTextSize(0.055);
	  cluster_population_1cm->SetTextColor(kRed);
	  cluster_population_1cm->SetTextAlign(12);
	  cluster_population_1cm->SetFillColor(kYellow);
	  cluster_population_1cm->AddText(Form("Expected Electrons per Cluster: %.1f - Track angle (deg) %0.1f",cluster_population,alpha));
	  //pt_1cm_alpha->Draw("same");
	  cluster_population_1cm->Draw("same");
	  
	  TH1F *h32=(TH1F*)file->Get(Form("H-Ch%i_signal/hTimeDifference_ch%i",i,i));
	  if (h32==0x0) { continue; }
	  if(isNov2021TestBeam) {timediff_1cm->cd(i-3);}
	  else if(isJuly2022TestBeam) {timediff_1cm->cd(counter_plots_1cm);}
	  TF1  *f4 = new TF1("f4","[0]*exp(-x/[1])",0,8);
	  //TF1  *f4 = new TF1("f4","[0]*exp(-x/[1])",0,12);
	  f4->SetParameters(0,4000.);
	  f4->SetParameters(1,3);
	  //   h32->Fit("f4","R");
	  //   if(isLogarithm){
	  //   	gPad->SetLogy(1);
	//   }
	  h32->Draw("same");
	  //pt_1cm_alpha->Draw("same");
	  //h31->GetXaxis()->SetRangeUser(0.,5.);
	  //h31->Fit("expo");
	  pt_1cm_alpha->Draw("same");
	  
	  
	  TH1F *h33=(TH1F*)file->Get(Form("H-Ch%i_signal/hTimeDifference_clust_ch%i",i,i));
	  if (h33==0x0) { continue; }
	  if(isNov2021TestBeam) {timediff_clust_1cm->cd(i-3);}
	  else if(isJuly2022TestBeam) {timediff_clust_1cm->cd(counter_plots_1cm);}
	  //pt_1cm_alpha->Draw("same");
	  TF1  *f1 = new TF1("f1","[0]*exp(-x/[1])",10,40);
	  f1->SetParameters(0,1000);
	  f1->SetParameters(1,10);
	//   h33->Fit("f1","R");
	//   if(isLogarithm){
	//   	gPad->SetLogy(1);
	//   }
	  h33->Draw("same");
	  
	  TH1F *h5=(TH1F*)file->Get(Form("H-Ch%i_signal/hHPeaks_ch%i",i,i));
	  if (h5==0x0) { continue; }
	  if(isNov2021TestBeam) {hpeaks_1cm->cd(i-3);}
	  else if(isJuly2022TestBeam) {hpeaks_1cm->cd(counter_plots_1cm);}
	  if(isRuns_90_10 || isRuns_80_20 || isRuns_85_15){
	    aveph_1cm = aveph_1cm + h5->GetMean();
	    rms_aveph_1cm = rms_aveph_1cm + h5->GetRMS();
	  }
	  h5->Draw("same");
	  
	  TH2F *h6=(TH2F*)file->Get(Form("H-Ch%i_signal/hHNPeaks_ch%i",i,i));
	  if (h6==0x0) {continue; }
	  if(isNov2021TestBeam) {hnpeaks_1cm->cd(i-3);}
	  else if(isJuly2022TestBeam) {hnpeaks_1cm->cd(counter_plots_1cm);}
	  h6->GetYaxis()->SetTitle("Height of Peaks found [V]");
	  h6->GetXaxis()->SetRangeUser(0.,90.);
	  h6->Draw("colz");
	  gPad->Update();
	  TPaveStats *st_1cm = (TPaveStats*)h6->FindObject("stats");
	  st_1cm->SetX1NDC(0.75); //new x start position
	  st_1cm->SetX2NDC(0.85); //new x end position
	  st_1cm->SetY1NDC(0.65); //new x start position
	  st_1cm->SetY2NDC(0.85); //new x end position
	  
	  TH2F *h40=(TH2F*)file->Get(Form("H-Ch%i_signal/hNClusterFCluster_ch%i",i,i));
	  if (h40==0x0) { continue; }
	  if(isNov2021TestBeam) {hncluster_1cm->cd(i-3);}
	  else if(isJuly2022TestBeam) {hncluster_1cm->cd(counter_plots_1cm);}
	  auto prof_px1_1 = h40->ProfileX();
	  prof_px1_1->GetYaxis()->SetTitle("Number of Clusters found");
	  prof_px1_1->Draw("same");
	  TF1 *fit_1cm=new TF1("fit_1cm","pol1",20,200.);     
	  prof_px1_1->Fit("fit_1cm","R");
	  //h40->Draw("colz");
	  gPad->Update();
	  
	  
	  TH2F *h41=(TH2F*)file->Get(Form("H-Ch%i_signal/hNPeakFPeak_ch%i",i,i));
	  if (h41==0x0) { continue; }
	  if(isNov2021TestBeam) {hnelectron_1cm->cd(i-3);}
	  else if(isJuly2022TestBeam) {hnelectron_1cm->cd(counter_plots_1cm);}
	  auto prof_px = h41->ProfileX();
	  prof_px->GetYaxis()->SetTitle("Number of Electron Peaks found");
	  prof_px->Draw("same");
	  //h41->Draw("colz");
	  TF1 *fit1_1cm=new TF1("fit1_1cm","pol1",20,200.);     
	  prof_px->Fit("fit1_1cm","R");
	  gPad->Update();
	  
	  TH1F *h7=(TH1F*)file->Get(Form("H-Ch%i_signal/hTPeaks_ch%i",i,i));
	  if (h7==0x0) { continue; }
	  if(isNov2021TestBeam) {tpeaks_1cm->cd(i-3);}
	  else if(isJuly2022TestBeam) {tpeaks_1cm->cd(counter_plots_1cm);}
	  h7->Draw("same");
	  
	  TH1F *h8=(TH1F*)file->Get(Form("H-Ch%i_signal/hTFstPeaks_ch%i",i,i));
	  if (h8==0x0) { continue; }
	  if(isNov2021TestBeam) {tfpeaks_1cm->cd(i-3);}
	  else if(isJuly2022TestBeam) {tfpeaks_1cm->cd(counter_plots_1cm);}
	  h8->Draw("same");
	  
	  TH1F *h9=(TH1F*)file->Get(Form("H-Ch%i_signal/hNPeaks_1_ch%i",i,i));
	  if (h9==0x0) { continue; }
	  if(isNov2021TestBeam) {tlpeaks_1cm->cd(i-3);}
	  else if(isJuly2022TestBeam) {tlpeaks_1cm->cd(counter_plots_1cm);}
	  h9->Draw("same");
	  
	  
	  TH1F *h12=(TH1F*)file->Get(Form("H-Ch%i_signal/hIntegN_ch%i",i,i));
	  if (h12==0x0) { continue; }
	  if(isRuns_90_10 || isRuns_80_20 || isRuns_85_15){
	    integral_1cm = integral_1cm+ h12->GetMean();
	    rms_integral_1cm = rms_integral_1cm + h12->GetRMS();
	  }
	  if(isNov2021TestBeam) {integ_1cm->cd(i-3);}
	  else if(isJuly2022TestBeam) {integ_1cm->cd(counter_plots_1cm);}
	  h12->Draw("same");
	  
	  
	  TH1F *h16=(TH1F*)file->Get(Form("H-Ch%i_signal/hRms_ch%i",i,i));
	  if (h16==0x0) { continue; }
	  if(isNov2021TestBeam) {rms_1cm->cd(i-3);}
	  else if(isJuly2022TestBeam) {rms_1cm->cd(counter_plots_1cm);}
	  if(isRuns_90_10 || isRuns_80_20 || isRuns_85_15){
	    rms_1cm_var = rms_1cm_var + h16->GetMean();
	    rms_rms_1cm = rms_rms_1cm + h16->GetRMS();
	  }
	  h16->Draw("same");
	  
	  counter_plots_1cm = counter_plots_1cm + 1;
	} //if condition for 1 cm tubes
	
      } //loop on all tubes
      
	  	
      
      if(isRuns_90_10 || isRuns_80_20 || isRuns_85_15){
	aveph_1cm = aveph_1cm/ (float) Channel_1cm_size;
	rms_aveph_1cm = rms_aveph_1cm/ (float) Channel_1cm_size;
	
	rms_1cm_var = rms_1cm_var/ (float) Channel_1cm_size;
	rms_rms_1cm = rms_rms_1cm/ (float) Channel_1cm_size;
	
	epc_1cm = epc_1cm/ (float) Channel_1cm_size;
	rms_epc_1cm = rms_epc_1cm/ (float) Channel_1cm_size;
	
	integral_1cm = (integral_1cm/ (float) Channel_1cm_size) ;
	rms_integral_1cm = (rms_integral_1cm/ (float) Channel_1cm_size);
	
	bsl_1cm_var = bsl_1cm_var/ (float) Channel_1cm_size;
	rms_bsl_1cm = rms_bsl_1cm/ (float) Channel_1cm_size;
	
	maximum_1cm = maximum_1cm/ (float) Channel_1cm_size;
	rms_maximum_1cm = rms_maximum_1cm/ (float) Channel_1cm_size;
	
	mean_clusters_1cm = mean_clusters_1cm/ (float) Channel_1cm_size;
	rms_clusters_1cm = rms_clusters_1cm/ (float) Channel_1cm_size;
	
	mean_electrons_1cm = mean_electrons_1cm/ (float) Channel_1cm_size;
	rms_electrons_1cm = rms_electrons_1cm/ (float) Channel_1cm_size;
	
	//Normalization
	
	mean_clusters_1cm = mean_clusters_1cm/expected_cluster;
	rms_clusters_1cm = rms_clusters_1cm/expected_cluster;
	
	mean_electrons_1cm = mean_electrons_1cm/expected_electrons;
	rms_electrons_1cm = rms_electrons_1cm/expected_electrons;
	
	aveph_summary_1cm->cd();
	aveph_summary_1cm->SetGridx();
  	aveph_summary_1cm->SetGridy();
	if(isScanAngleStudy){
	  gr_aveph_summary_1cm->SetPoint(counter_filling_electrons_1cm,alpha,aveph_1cm);
	}
	else if(isScanHVStudy || isGasMixtureStudy){
	  gr_aveph_summary_1cm->SetPoint(counter_filling_electrons_1cm,counter_filling_electrons_1cm,aveph_1cm);
	}
	else if(isScanSamplingStudy){
	  gr_aveph_summary_1cm->SetPoint(counter_filling_electrons_1cm,_gsample,aveph_1cm);
	}
	gr_aveph_summary_1cm->SetPointError(counter_filling_electrons_1cm,0.,rms_aveph_1cm);
	gr_aveph_summary_1cm->SetMarkerStyle(21);
	gr_aveph_summary_1cm->SetMarkerSize(0.5);
	gr_aveph_summary_1cm->SetMarkerColor(kBlue);
	gr_aveph_summary_1cm->SetLineColor(kBlue);
	TLatex *latex_aveph_summary_1cm = new TLatex(gr_aveph_summary_1cm->GetX()[number_of_trials], gr_aveph_summary_1cm->GetY()[number_of_trials],Form("%s DERIV N1 %s N2 %s N3 %s N4 %s scale_cut %s",gasMixture.Data(),N1.Data(),N2.Data(),N3.Data(),N4.Data(),scale_cut.Data())); 
	latex_aveph_summary_1cm->SetTextSize(0.01); 
	latex_aveph_summary_1cm->SetTextColor(kRed);
	gr_aveph_summary_1cm->GetListOfFunctions()->Add(latex_aveph_summary_1cm);
	gr_aveph_summary_1cm->Draw("sameAP");
	
	
	rms_summary_1cm->cd(); 
	rms_summary_1cm->SetGridx();
  	rms_summary_1cm->SetGridy();
	if(isScanAngleStudy){
	gr_rms_summary_1cm->SetPoint(counter_filling_electrons_1cm,alpha,rms_1cm_var);
	}
	else if(isScanHVStudy || isGasMixtureStudy){
	  gr_rms_summary_1cm->SetPoint(counter_filling_electrons_1cm,counter_filling_electrons_1cm,rms_1cm_var);
	}
	else if(isScanSamplingStudy){
	gr_rms_summary_1cm->SetPoint(counter_filling_electrons_1cm,_gsample,rms_1cm_var);
	} 
	gr_rms_summary_1cm->SetPointError(counter_filling_electrons_1cm,0.,rms_rms_1cm);
	gr_rms_summary_1cm->SetMarkerStyle(21);
	gr_rms_summary_1cm->SetMarkerSize(0.5);
	gr_rms_summary_1cm->SetMarkerColor(kBlue);
	gr_rms_summary_1cm->SetLineColor(kBlue);
	TLatex *latex_rms_summary_1cm = new TLatex(gr_rms_summary_1cm->GetX()[number_of_trials], gr_rms_summary_1cm->GetY()[number_of_trials],Form("%s DERIV N1 %s N2 %s N3 %s N4 %s scale_cut %s",gasMixture.Data(),N1.Data(),N2.Data(),N3.Data(),N4.Data(),scale_cut.Data())); 
	latex_rms_summary_1cm->SetTextSize(0.01); 
	latex_rms_summary_1cm->SetTextColor(kRed);
	gr_rms_summary_1cm->GetListOfFunctions()->Add(latex_rms_summary_1cm);
	gr_rms_summary_1cm->Draw("sameAP");
	
	epc_summary_1cm->cd();
	epc_summary_1cm->SetGridx();
  	epc_summary_1cm->SetGridy();
	if(isScanAngleStudy){
	  gr_epc_summary_1cm->SetPoint(counter_filling_electrons_1cm,alpha,epc_1cm);
	}
	else if(isScanHVStudy || isGasMixtureStudy){
	  gr_epc_summary_1cm->SetPoint(counter_filling_electrons_1cm,counter_filling_electrons_1cm,epc_1cm);
	}
	else if(isScanSamplingStudy){
	  gr_epc_summary_1cm->SetPoint(counter_filling_electrons_1cm,_gsample,epc_1cm);
	}
	gr_epc_summary_1cm->SetPointError(counter_filling_electrons_1cm,0.,rms_epc_1cm);
	gr_epc_summary_1cm->SetMarkerStyle(21);
	gr_epc_summary_1cm->SetMarkerSize(0.5);
	gr_epc_summary_1cm->SetMarkerColor(kBlue);
	gr_epc_summary_1cm->SetLineColor(kBlue);
	TLatex *latex_epc_summary_1cm = new TLatex(gr_epc_summary_1cm->GetX()[number_of_trials], gr_epc_summary_1cm->GetY()[number_of_trials],Form("%s DERIV N1 %s N2 %s N3 %s N4 %s scale_cut %s",gasMixture.Data(),N1.Data(),N2.Data(),N3.Data(),N4.Data(),scale_cut.Data())); 
	latex_epc_summary_1cm->SetTextSize(0.01); 
	latex_epc_summary_1cm->SetTextColor(kRed);
	gr_epc_summary_1cm->GetListOfFunctions()->Add(latex_epc_summary_1cm);
	gr_epc_summary_1cm->Draw("sameAP");
	
	// TLegend *leg= new TLegend(0.5,0.75,0.85,0.85); 
	// leg->AddEntry(tmpsignal_1.back(),"Waveform"); 
	
	bsl_summary_1cm->cd();
	bsl_summary_1cm->SetGridx();
  	bsl_summary_1cm->SetGridy();
	if(isScanAngleStudy){
	  gr_bsl_summary_1cm->SetPoint(counter_filling_electrons_1cm,alpha,bsl_1cm_var);
	}
	else if(isScanHVStudy || isGasMixtureStudy){
	  gr_bsl_summary_1cm->SetPoint(counter_filling_electrons_1cm,counter_filling_electrons_1cm,bsl_1cm_var);
	}
	else if(isScanSamplingStudy){
	  gr_bsl_summary_1cm->SetPoint(counter_filling_electrons_1cm,_gsample,bsl_1cm_var);
	}
	gr_bsl_summary_1cm->SetPointError(counter_filling_electrons_1cm,0.,rms_bsl_1cm);
	gr_bsl_summary_1cm->SetMarkerStyle(21);
	gr_bsl_summary_1cm->SetMarkerSize(0.5);
	gr_bsl_summary_1cm->SetMarkerColor(kBlue);
	gr_bsl_summary_1cm->SetLineColor(kBlue);
	TLatex *latex_bsl_summary_1cm = new TLatex(gr_bsl_summary_1cm->GetX()[number_of_trials], gr_bsl_summary_1cm->GetY()[number_of_trials],Form("%s DERIV N1 %s N2 %s N3 %s N4 %s scale_cut %s",gasMixture.Data(),N1.Data(),N2.Data(),N3.Data(),N4.Data(),scale_cut.Data())); 
	latex_bsl_summary_1cm->SetTextSize(0.01); 
	latex_bsl_summary_1cm->SetTextColor(kRed);
	gr_bsl_summary_1cm->GetListOfFunctions()->Add(latex_bsl_summary_1cm);
	gr_bsl_summary_1cm->Draw("sameAP");
	
	maximum_summary_1cm->cd(); 
	maximum_summary_1cm->SetGridx();
  	maximum_summary_1cm->SetGridy();
	if(isScanAngleStudy){  
	  gr_maximum_summary_1cm->SetPoint(counter_filling_electrons_1cm,alpha,maximum_1cm);
	}
	else if(isScanHVStudy || isGasMixtureStudy){
	  gr_maximum_summary_1cm->SetPoint(counter_filling_electrons_1cm,counter_filling_electrons_1cm,maximum_1cm);
	}
	else if(isScanSamplingStudy){
	  gr_maximum_summary_1cm->SetPoint(counter_filling_electrons_1cm,_gsample,maximum_1cm);
	}
	gr_maximum_summary_1cm->SetPointError(counter_filling_electrons_1cm,0.,rms_maximum_1cm);
	gr_maximum_summary_1cm->SetMarkerStyle(21);
	gr_maximum_summary_1cm->SetMarkerSize(0.5);
	gr_maximum_summary_1cm->SetMarkerColor(kBlue);
	gr_maximum_summary_1cm->SetLineColor(kBlue);
	TLatex *latex_maximum_summary_1cm = new TLatex(gr_maximum_summary_1cm->GetX()[number_of_trials], gr_maximum_summary_1cm->GetY()[number_of_trials],Form("%s DERIV N1 %s N2 %s N3 %s N4 %s scale_cut %s",gasMixture.Data(),N1.Data(),N2.Data(),N3.Data(),N4.Data(),scale_cut.Data())); 
	latex_maximum_summary_1cm->SetTextSize(0.01); 
	latex_maximum_summary_1cm->SetTextColor(kRed);
	gr_maximum_summary_1cm->GetListOfFunctions()->Add(latex_maximum_summary_1cm);
	gr_maximum_summary_1cm->Draw("sameAP");
	
	integral_summary_1cm->cd();
	integral_summary_1cm->SetGridx();
  	integral_summary_1cm->SetGridy();
	if(isScanAngleStudy){
	  gr_integral_summary_1cm->SetPoint(counter_filling_electrons_1cm,alpha,integral_1cm);
	}
	else if(isScanHVStudy || isGasMixtureStudy){
	gr_integral_summary_1cm->SetPoint(counter_filling_electrons_1cm,counter_filling_electrons_1cm,integral_1cm);
	}
	else if(isScanSamplingStudy){
	  gr_integral_summary_1cm->SetPoint(counter_filling_electrons_1cm,_gsample,integral_1cm);
	}
	gr_integral_summary_1cm->SetPointError(counter_filling_electrons_1cm,0.,rms_integral_1cm);
	gr_integral_summary_1cm->SetMarkerStyle(21);
	gr_integral_summary_1cm->SetMarkerSize(0.5);
	gr_integral_summary_1cm->SetMarkerColor(kBlue);
	gr_integral_summary_1cm->SetLineColor(kBlue);
	TLatex *latex_integral_summary_1cm = new TLatex(gr_integral_summary_1cm->GetX()[number_of_trials], gr_integral_summary_1cm->GetY()[number_of_trials],Form("%s DERIV N1 %s N2 %s N3 %s N4 %s scale_cut %s",gasMixture.Data(),N1.Data(),N2.Data(),N3.Data(),N4.Data(),scale_cut.Data())); 
	latex_integral_summary_1cm->SetTextSize(0.01); 
	latex_integral_summary_1cm->SetTextColor(kRed);
	gr_integral_summary_1cm->GetListOfFunctions()->Add(latex_integral_summary_1cm);
	gr_integral_summary_1cm->Draw("sameAP");
	
	efficiency_electrons_1cm->cd();
	efficiency_electrons_1cm->SetGridx();
  	efficiency_electrons_1cm->SetGridy();
	if(isScanAngleStudy){
	  gr_efficiency_electrons_1cm->SetPoint(counter_filling_electrons_1cm,alpha,mean_electrons_1cm);
	}
	else if(isScanHVStudy || isGasMixtureStudy){
	  gr_efficiency_electrons_1cm->SetPoint(counter_filling_electrons_1cm,counter_filling_electrons_1cm,mean_electrons_1cm);
	}
	else if(isScanSamplingStudy){
	  gr_efficiency_electrons_1cm->SetPoint(counter_filling_electrons_1cm,_gsample,mean_electrons_1cm);
	}
	gr_efficiency_electrons_1cm->SetPointError(counter_filling_electrons_1cm,0.,rms_electrons_1cm);
	gr_efficiency_electrons_1cm->SetMarkerStyle(21);
	gr_efficiency_electrons_1cm->SetMarkerSize(0.5);
	gr_efficiency_electrons_1cm->SetMarkerColor(kBlue);
	gr_efficiency_electrons_1cm->SetLineColor(kBlue);
	TLatex *latex_efficiency_electrons_1cm = new TLatex(gr_efficiency_electrons_1cm->GetX()[number_of_trials], gr_efficiency_electrons_1cm->GetY()[number_of_trials],Form("%s DERIV N1 %s N2 %s N3 %s N4 %s scale_cut %s",gasMixture.Data(),N1.Data(),N2.Data(),N3.Data(),N4.Data(),scale_cut.Data())); 
	latex_efficiency_electrons_1cm->SetTextSize(0.01); 
	latex_efficiency_electrons_1cm->SetTextColor(kRed);
	gr_efficiency_electrons_1cm->GetListOfFunctions()->Add(latex_efficiency_electrons_1cm);
	gr_efficiency_electrons_1cm->Draw("sameAP");
	
	efficiency_clusters_1cm->cd();
	efficiency_clusters_1cm->SetGridx();
  	efficiency_clusters_1cm->SetGridy();
	if(isScanAngleStudy){
	  gr_efficiency_clusters_1cm->SetPoint(counter_filling_clusters_1cm,alpha,mean_clusters_1cm);
	}
	else if(isScanHVStudy || isGasMixtureStudy){
	  gr_efficiency_clusters_1cm->SetPoint(counter_filling_clusters_1cm,counter_filling_clusters_1cm,mean_clusters_1cm);
	}
	else if(isScanSamplingStudy){
	  gr_efficiency_clusters_1cm->SetPoint(counter_filling_clusters_1cm,_gsample,mean_clusters_1cm);
	}
	gr_efficiency_clusters_1cm->SetPointError(counter_filling_clusters_1cm,0.,rms_clusters_1cm);
	gr_efficiency_clusters_1cm->SetMarkerStyle(21);
	gr_efficiency_clusters_1cm->SetMarkerSize(0.5);
	gr_efficiency_clusters_1cm->SetMarkerColor(kBlue);
	gr_efficiency_clusters_1cm->SetLineColor(kBlue);
	TLatex *latex_efficiency_clusters_1cm = new TLatex(gr_efficiency_clusters_1cm->GetX()[number_of_trials], gr_efficiency_clusters_1cm->GetY()[number_of_trials],Form("%s DERIV N1 %s N2 %s N3 %s N4 %s scale_cut %s",gasMixture.Data(),N1.Data(),N2.Data(),N3.Data(),N4.Data(),scale_cut.Data())); 
	latex_efficiency_clusters_1cm->SetTextSize(0.01); 
	latex_efficiency_clusters_1cm->SetTextColor(kRed);
	gr_efficiency_clusters_1cm->GetListOfFunctions()->Add(latex_efficiency_clusters_1cm);
	gr_efficiency_clusters_1cm->Draw("sameAP");
	
	counter_filling_clusters_1cm = counter_filling_clusters_1cm +1;
	counter_filling_electrons_1cm = counter_filling_electrons_1cm +1;
	
	
	if((name_file_compact == Runs_90_10[Runs_90_10_size - 1]) || (name_file_compact == Runs_80_20[Runs_80_20_size - 1])){
	  efficiency_clusters_1cm->SaveAs("/lustrehome/bdanzi/TestBeam2022/analysistestbeam2022/analysis_2021/drifttubes_offline_analysis/testbeam_analysis/Plots/efficiency_clusters_1cm.pdf");
	  efficiency_electrons_1cm->SaveAs("/lustrehome/bdanzi/TestBeam2022/analysistestbeam2022/analysis_2021/drifttubes_offline_analysis/testbeam_analysis/Plots/efficiency_electrons_1cm.pdf");
	  epc_summary_1cm->SaveAs("/lustrehome/bdanzi/TestBeam2022/analysistestbeam2022/analysis_2021/drifttubes_offline_analysis/testbeam_analysis/Plots/epc_summary_1cm.pdf");
	  integral_summary_1cm->SaveAs("/lustrehome/bdanzi/TestBeam2022/analysistestbeam2022/analysis_2021/drifttubes_offline_analysis/testbeam_analysis/Plots/integral_summary_1cm.pdf");
	  maximum_summary_1cm->SaveAs("/lustrehome/bdanzi/TestBeam2022/analysistestbeam2022/analysis_2021/drifttubes_offline_analysis/testbeam_analysis/Plots/maximum_summary_1cm.pdf");
	  bsl_summary_1cm->SaveAs("/lustrehome/bdanzi/TestBeam2022/analysistestbeam2022/analysis_2021/drifttubes_offline_analysis/testbeam_analysis/Plots/bsl_summary_1cm.pdf");
	  rms_summary_1cm->SaveAs("/lustrehome/bdanzi/TestBeam2022/analysistestbeam2022/analysis_2021/drifttubes_offline_analysis/testbeam_analysis/Plots/rms_summary_1cm.pdf");
	  aveph_summary_1cm->SaveAs("/lustrehome/bdanzi/TestBeam2022/analysistestbeam2022/analysis_2021/drifttubes_offline_analysis/testbeam_analysis/Plots/aveph_summary_1cm.pdf");
	  
	  efficiency_clusters_1cm->SaveAs("/lustrehome/bdanzi/TestBeam2022/analysistestbeam2022/analysis_2021/drifttubes_offline_analysis/testbeam_analysis/Plots/efficiency_clusters_1cm.png");
	  efficiency_electrons_1cm->SaveAs("/lustrehome/bdanzi/TestBeam2022/analysistestbeam2022/analysis_2021/drifttubes_offline_analysis/testbeam_analysis/Plots/efficiency_electrons_1cm.png");
	  epc_summary_1cm->SaveAs("/lustrehome/bdanzi/TestBeam2022/analysistestbeam2022/analysis_2021/drifttubes_offline_analysis/testbeam_analysis/Plots/epc_summary_1cm.png");
	  integral_summary_1cm->SaveAs("/lustrehome/bdanzi/TestBeam2022/analysistestbeam2022/analysis_2021/drifttubes_offline_analysis/testbeam_analysis/Plots/integral_summary_1cm.png");
	  maximum_summary_1cm->SaveAs("/lustrehome/bdanzi/TestBeam2022/analysistestbeam2022/analysis_2021/drifttubes_offline_analysis/testbeam_analysis/Plots/maximum_summary_1cm.png");
	  bsl_summary_1cm->SaveAs("/lustrehome/bdanzi/TestBeam2022/analysistestbeam2022/analysis_2021/drifttubes_offline_analysis/testbeam_analysis/Plots/bsl_summary_1cm.png");
	  rms_summary_1cm->SaveAs("/lustrehome/bdanzi/TestBeam2022/analysistestbeam2022/analysis_2021/drifttubes_offline_analysis/testbeam_analysis/Plots/rms_summary_1cm.png");
	  aveph_summary_1cm->SaveAs("/lustrehome/bdanzi/TestBeam2022/analysistestbeam2022/analysis_2021/drifttubes_offline_analysis/testbeam_analysis/Plots/aveph_summary_1cm.png");
	  
	}
      }
      
      
      bool savePlots = true;
      
      if(savePlots){
	timediff_clust_1cm->SaveAs(Form("/lustrehome/bdanzi/TestBeam2022/analysistestbeam2022/analysis_2021/drifttubes_offline_analysis/testbeam_analysis/Plots/%s/clust_difference_1cm.pdf",fname.Data()));
	timediff_1cm->SaveAs(Form("/lustrehome/bdanzi/TestBeam2022/analysistestbeam2022/analysis_2021/drifttubes_offline_analysis/testbeam_analysis/Plots/%s/electrons_difference_1cm.pdf",fname.Data()));
	npeaks_1cm->SaveAs(Form("/lustrehome/bdanzi/TestBeam2022/analysistestbeam2022/analysis_2021/drifttubes_offline_analysis/testbeam_analysis/Plots/%s/nelectrons_1cm.pdf",fname.Data()));
	//tpeaks_1cm->SaveAs(Form("/lustrehome/bdanzi/TestBeam2022/analysistestbeam2022/analysis_2021/drifttubes_offline_analysis/testbeam_analysis/Plots/%s/tpeaks_1cm.pdf",fname.Data()));
	tfpeaks_1cm->SaveAs(Form("/lustrehome/bdanzi/TestBeam2022/analysistestbeam2022/analysis_2021/drifttubes_offline_analysis/testbeam_analysis/Plots/%s/tfpeaks_1cm.pdf",fname.Data()));
	tlpeaks_1cm->SaveAs(Form("/lustrehome/bdanzi/TestBeam2022/analysistestbeam2022/analysis_2021/drifttubes_offline_analysis/testbeam_analysis/Plots/%s/tlpeaks_1cm.pdf",fname.Data()));
	//hnpeaks_1cm->SaveAs(Form("/lustrehome/bdanzi/TestBeam2022/analysistestbeam2022/analysis_2021/drifttubes_offline_analysis/testbeam_analysis/Plots/%s/hnpeaks_1cm.pdf",fname.Data()));
	hncluster_1cm->SaveAs(Form("/lustrehome/bdanzi/TestBeam2022/analysistestbeam2022/analysis_2021/drifttubes_offline_analysis/testbeam_analysis/Plots/%s/hncluster_1cm.pdf",fname.Data()));
	hpeaks_1cm->SaveAs(Form("/lustrehome/bdanzi/TestBeam2022/analysistestbeam2022/analysis_2021/drifttubes_offline_analysis/testbeam_analysis/Plots/%s/hpeaks_1cm.pdf",fname.Data()));
	bsl_1cm->SaveAs(Form("/lustrehome/bdanzi/TestBeam2022/analysistestbeam2022/analysis_2021/drifttubes_offline_analysis/testbeam_analysis/Plots/%s/bsl_1cm.pdf",fname.Data()));
	max_1cm->SaveAs(Form("/lustrehome/bdanzi/TestBeam2022/analysistestbeam2022/analysis_2021/drifttubes_offline_analysis/testbeam_analysis/Plots/%s/max_1cm.pdf",fname.Data()));
	integ_1cm->SaveAs(Form("/lustrehome/bdanzi/TestBeam2022/analysistestbeam2022/analysis_2021/drifttubes_offline_analysis/testbeam_analysis/Plots/%s/integ_1cm.pdf",fname.Data()));
	rms_1cm->SaveAs(Form("/lustrehome/bdanzi/TestBeam2022/analysistestbeam2022/analysis_2021/drifttubes_offline_analysis/testbeam_analysis/Plots/%s/rms_1cm.pdf",fname.Data()));
	//min->SaveAs(Form("/lustrehome/bdanzi/TestBeam2022/analysistestbeam2022/analysis_2021/drifttubes_offline_analysis/testbeam_analysis/Plots/%s/min_ch%i.pdf",fname.Data(),i));
	npeaks_clustser_1cm->SaveAs(Form("/lustrehome/bdanzi/TestBeam2022/analysistestbeam2022/analysis_2021/drifttubes_offline_analysis/testbeam_analysis/Plots/%s/nclusters_1cm.pdf",fname.Data()));
	cluster_population_canvas_1cm->SaveAs(Form("/lustrehome/bdanzi/TestBeam2022/analysistestbeam2022/analysis_2021/drifttubes_offline_analysis/testbeam_analysis/Plots/%s/cluster_population_1cm.pdf",fname.Data()));
	hnelectron_1cm->SaveAs(Form("/lustrehome/bdanzi/TestBeam2022/analysistestbeam2022/analysis_2021/drifttubes_offline_analysis/testbeam_analysis/Plots/%s/hnelectron_1cm.pdf",fname.Data()));
	
	timediff_clust_1cm->SaveAs(Form("/lustrehome/bdanzi/TestBeam2022/analysistestbeam2022/analysis_2021/drifttubes_offline_analysis/testbeam_analysis/Plots/%s/clust_difference_1cm.png",fname.Data()));
	timediff_1cm->SaveAs(Form("/lustrehome/bdanzi/TestBeam2022/analysistestbeam2022/analysis_2021/drifttubes_offline_analysis/testbeam_analysis/Plots/%s/electrons_difference_1cm.png",fname.Data()));
	npeaks_1cm->SaveAs(Form("/lustrehome/bdanzi/TestBeam2022/analysistestbeam2022/analysis_2021/drifttubes_offline_analysis/testbeam_analysis/Plots/%s/nelectrons_1cm.png",fname.Data()));
	//tpeaks_1cm->SaveAs(Form("/lustrehome/bdanzi/TestBeam2022/analysistestbeam2022/analysis_2021/drifttubes_offline_analysis/testbeam_analysis/Plots/%s/tpeaks_1cm.pdf",fname.Data()));
	tfpeaks_1cm->SaveAs(Form("/lustrehome/bdanzi/TestBeam2022/analysistestbeam2022/analysis_2021/drifttubes_offline_analysis/testbeam_analysis/Plots/%s/tfpeaks_1cm.png",fname.Data()));
	tlpeaks_1cm->SaveAs(Form("/lustrehome/bdanzi/TestBeam2022/analysistestbeam2022/analysis_2021/drifttubes_offline_analysis/testbeam_analysis/Plots/%s/tlpeaks_1cm.png",fname.Data()));
	//hnpeaks_1cm->SaveAs(Form("/lustrehome/bdanzi/TestBeam2022/analysistestbeam2022/analysis_2021/drifttubes_offline_analysis/testbeam_analysis/Plots/%s/hnpeaks_1cm.pdf",fname.Data()));
	hncluster_1cm->SaveAs(Form("/lustrehome/bdanzi/TestBeam2022/analysistestbeam2022/analysis_2021/drifttubes_offline_analysis/testbeam_analysis/Plots/%s/hncluster_1cm.png",fname.Data()));
	hpeaks_1cm->SaveAs(Form("/lustrehome/bdanzi/TestBeam2022/analysistestbeam2022/analysis_2021/drifttubes_offline_analysis/testbeam_analysis/Plots/%s/hpeaks_1cm.png",fname.Data()));
	bsl_1cm->SaveAs(Form("/lustrehome/bdanzi/TestBeam2022/analysistestbeam2022/analysis_2021/drifttubes_offline_analysis/testbeam_analysis/Plots/%s/bsl_1cm.png",fname.Data()));
	max_1cm->SaveAs(Form("/lustrehome/bdanzi/TestBeam2022/analysistestbeam2022/analysis_2021/drifttubes_offline_analysis/testbeam_analysis/Plots/%s/max_1cm.png",fname.Data()));
	integ_1cm->SaveAs(Form("/lustrehome/bdanzi/TestBeam2022/analysistestbeam2022/analysis_2021/drifttubes_offline_analysis/testbeam_analysis/Plots/%s/integ_1cm.png",fname.Data()));
	rms_1cm->SaveAs(Form("/lustrehome/bdanzi/TestBeam2022/analysistestbeam2022/analysis_2021/drifttubes_offline_analysis/testbeam_analysis/Plots/%s/rms_1cm.png",fname.Data()));
	//min->SaveAs(Form("/lustrehome/bdanzi/TestBeam2022/analysistestbeam2022/analysis_2021/drifttubes_offline_analysis/testbeam_analysis/Plots/%s/min_ch%i.pdf",fname.Data(),i));
	npeaks_clustser_1cm->SaveAs(Form("/lustrehome/bdanzi/TestBeam2022/analysistestbeam2022/analysis_2021/drifttubes_offline_analysis/testbeam_analysis/Plots/%s/nclusters_1cm.png",fname.Data()));
	cluster_population_canvas_1cm->SaveAs(Form("/lustrehome/bdanzi/TestBeam2022/analysistestbeam2022/analysis_2021/drifttubes_offline_analysis/testbeam_analysis/Plots/%s/cluster_population_1cm.png",fname.Data()));
	hnelectron_1cm->SaveAs(Form("/lustrehome/bdanzi/TestBeam2022/analysistestbeam2022/analysis_2021/drifttubes_offline_analysis/testbeam_analysis/Plots/%s/hnelectron_1cm.png",fname.Data()));
	
      }
      
      
      
      for(int i = 0; i<=channel; ++i){ //looping over all channels
	isChannel_2cm = 0;
	isChannel_1p5cm = 0;
	
	for (int j = 0; j < Channel_2cm_size; j++) {
	  if (Channel_2cm[j] == i) {
            isChannel_2cm = 1;
            break;
	  }
	}
	
	for (int j = 0; j < Channel_1p5cm_size; j++) {
	  if (Channel_1p5cm[j] == i) {
            isChannel_1p5cm = 1;
            break;
	  }
	}
	
	
	if((isChannel_2cm && isNov2021TestBeam) || (isChannel_1p5cm && isJuly2022TestBeam)){
	  //2cm/1cm 1,8/0,8 = 2,25
	  if(isChannel_2cm){
	    drift_size = 1.8;
	  }
	  else if(isChannel_1p5cm){
	    drift_size = 1.2;
	  }
	  expected_electrons = cluster_per_cm_mip * drift_size*relativistic_rise * cluster_population * 1/cos_alpha;
	  expected_cluster = cluster_per_cm_mip * drift_size*relativistic_rise * 1/cos_alpha;
	  
	  TH1F *h17=(TH1F*)file->Get(Form("H-Ch%i_signal/hMaxVInR_ch%i",i,i));
	  if (h17==0x0) { continue; }
	  if(isNov2021TestBeam) {max_2cm->cd(i-9);}
	  else if(isJuly2022TestBeam) {max_2cm->cd(counter_plots_1p5cm);}
	  if(isRuns_90_10 || isRuns_80_20 || isRuns_85_15){
	    maximum_2cm = maximum_2cm + h17->GetMean();
	    rms_maximum_2cm = rms_maximum_2cm + h17->GetRMS();
	  }
	  //max->cd(2);
	  h17->Fit("landau");
	  gPad->SetLogy(0);
	  gPad->SetLogx(0);
	  h17->Draw("same");
	  
	  TH1F *h3=(TH1F*)file->Get(Form("H-Ch%i_signal/hBsl_ch%i",i,i));
	  if (h3==0x0) { continue; }
	  if(isNov2021TestBeam) {bsl_2cm->cd(i-9);}
	  else if(isJuly2022TestBeam) {bsl_2cm->cd(counter_plots_1p5cm);}
	  if(isRuns_90_10 || isRuns_80_20 || isRuns_85_15){
	    bsl_2cm_var = bsl_2cm_var+ h3->GetMean();
	    rms_bsl_2cm = rms_bsl_2cm + h3->GetRMS();
	  }
	  gPad->SetLogy(0);
	  gPad->SetLogx(0);
	  h3->GetXaxis()->SetRangeUser(-0.5,-0.4);
	  h3->Draw("same");
	  
	  
	  
	  TH1F *h4=(TH1F*)file->Get(Form("H-Ch%i_signal/hNPeaks_ch%i",i,i));
	  if (h4==0x0) { continue; }
	  if(isNov2021TestBeam) {npeaks_2cm->cd(i-9);}
	  else if(isJuly2022TestBeam) {npeaks_2cm->cd(counter_plots_1p5cm);}
	  if(isRuns_90_10 || isRuns_80_20 || isRuns_85_15){
	    mean_electrons_2cm = mean_electrons_2cm + h4->GetMean();
	    rms_electrons_2cm = rms_electrons_2cm + h4->GetRMS();
	  }
	  //h4->Fit("landau");
	  if(isLogarithm){
	    gPad->SetLogy(1);
	  }
	  //   h4->Fit("landau");
	  h4->Draw("same");
	  TPaveText *pt_2cm = new TPaveText(0.11,0.84,0.88,0.88,"NDC");
	  pt_2cm->SetTextSize(0.06);
	  pt_2cm->SetTextColor(kRed);
	  pt_2cm->SetFillColor(kYellow);
	  pt_2cm->SetTextAlign(12);
	  pt_2cm->AddText(Form("Expected Electron Peaks: %.1f - Track angle (deg) %0.1f",expected_electrons,alpha));
	  pt_2cm->Draw("same");
	  TPaveText *pt_2cm_alpha = new TPaveText(0.2,0.84,0.62,0.88,"NDC");
	  pt_2cm_alpha->SetTextSize(0.05);
	  pt_2cm_alpha->SetTextColor(kRed);
	  pt_2cm_alpha->SetFillColor(0);
	  pt_2cm_alpha->SetTextAlign(12);
	  pt_2cm_alpha->AddText(Form("Track angle (deg): %.1f",alpha));
	  //pt_2cm_alpha->Draw("same");
	  
	  TH1F *h5=(TH1F*)file->Get(Form("H-Ch%i_signal/hHPeaks_ch%i",i,i));
	  if (h5==0x0) { continue; }
	  if(isNov2021TestBeam) {hpeaks_2cm->cd(i-9);}
	  else if(isJuly2022TestBeam) {hpeaks_2cm->cd(counter_plots_1p5cm);}
	  if(isRuns_90_10 || isRuns_80_20 || isRuns_85_15){
	    aveph_2cm = aveph_2cm+ h5->GetMean();
	    rms_aveph_2cm = rms_aveph_2cm + h5->GetRMS();
	  }
	  h5->Draw("same");
	  
	  TH2F *h6=(TH2F*)file->Get(Form("H-Ch%i_signal/hHNPeaks_ch%i",i,i));
	  if (h6==0x0) { continue; }
	  if(isNov2021TestBeam) {hnpeaks_2cm->cd(i-9);}
	  else if(isJuly2022TestBeam) {hnpeaks_2cm->cd(counter_plots_1p5cm);}
	  h6->GetYaxis()->SetTitle("Height of Electron Peaks found [V]");
	  h6->GetXaxis()->SetRangeUser(0.,200.);
	  h6->Draw( "colz");
	  gPad->Update();
	  TPaveStats *st_2cm = (TPaveStats*)h6->FindObject("stats");
	  st_2cm->SetX1NDC(0.75); //new x start position
	  st_2cm->SetX2NDC(0.85); //new x end position
	  st_2cm->SetY1NDC(0.65); //new x start position
	  st_2cm->SetY2NDC(0.85); //new x end position
	  
	  TH2F *h40=(TH2F*)file->Get(Form("H-Ch%i_signal/hNClusterFCluster_ch%i",i,i));
	  if (h40==0x0) { continue; }
	  if(isNov2021TestBeam) {hncluster_2cm->cd(i-9);}
	  else if(isJuly2022TestBeam) {hncluster_2cm->cd(counter_plots_1p5cm);}
	  auto prof_px2_2 = h40->ProfileX();
	  prof_px2_2->GetYaxis()->SetTitle("Number of Clusters found");
	  prof_px2_2->Draw("same");
	  TF1 *fit_2cm=new TF1("fit_2cm","pol1",20,450.);       
	  prof_px2_2->Fit("fit_2cm","R");
	  //h40->Draw("colz");
	  gPad->Update();
	  //TPaveStats *st_cluster_2cm = (TPaveStats*)h40->FindObject("stats");
	  //st_cluster_2cm->SetX1NDC(0.75); //new x start position
	  //st_cluster_2cm->SetX2NDC(0.85); //new x end position
	  //st_cluster_2cm->SetY1NDC(0.65); //new x start position
	  //st_cluster_2cm->SetY2NDC(0.85); //new x end position
	  
	  TH2F *h41=(TH2F*)file->Get(Form("H-Ch%i_signal/hNPeakFPeak_ch%i",i,i));
	  if (h41==0x0) { continue; }
	  if(isNov2021TestBeam) {hnelectron_2cm->cd(i-9);}
	  else if(isJuly2022TestBeam) {hnelectron_2cm->cd(counter_plots_1p5cm);}
	  auto prof_px_2 = h41->ProfileX();
	  prof_px_2->GetYaxis()->SetTitle("Number of Electron Peaks found");
	  prof_px_2->Draw("same");
	  TF1 *fit2_2cm=new TF1("fit2_2cm","pol1",20,450.);       
	  prof_px_2->Fit("fit2_2cm","R");
	  //h41->Draw("colz");
	  gPad->Update();
	  //TPaveStats *st_electron_2cm = (TPaveStats*)h41->FindObject("stats");
	  //st_electron_2cm->SetX1NDC(0.75); //new x start position
	  //st_electron_2cm->SetX2NDC(0.85); //new x end position
	  //st_electron_2cm->SetY1NDC(0.65); //new x start position
	  //st_electron_2cm->SetY2NDC(0.85); //new x end position
	  
	  TH1F *h7=(TH1F*)file->Get(Form("H-Ch%i_signal/hTPeaks_ch%i",i,i));
	  if (h7==0x0) { continue; }
	  if(isNov2021TestBeam) {tpeaks_2cm->cd(i-9);}
	  else if(isJuly2022TestBeam) {tpeaks_2cm->cd(counter_plots_1p5cm);}
	  h7->Draw( "same");
	  
	  TH1F *h8=(TH1F*)file->Get(Form("H-Ch%i_signal/hTFstPeaks_ch%i",i,i));
	  if (h8==0x0) { continue; }
	  if(isNov2021TestBeam) {tfpeaks_2cm->cd(i-9);}
	  else if(isJuly2022TestBeam) {tfpeaks_2cm->cd(counter_plots_1p5cm);}
	  h8->Draw("same" );
	  
	  TH1F *h9=(TH1F*)file->Get(Form("H-Ch%i_signal/hNPeaks_1_ch%i",i,i));
	  if (h9==0x0) { continue; }
	  if(isNov2021TestBeam) {tlpeaks_2cm->cd(i-9);}
	  else if(isJuly2022TestBeam) {tlpeaks_2cm->cd(counter_plots_1p5cm);}
	  h9->Draw( "same");
	  
	  TH1F *h33=(TH1F*)file->Get(Form("H-Ch%i_signal/hTimeDifference_clust_ch%i",i,i));
	  if (h33==0x0) { continue; }
	  if(isNov2021TestBeam) {timediff_clust_2cm->cd(i-9);}
	  else if(isJuly2022TestBeam) {timediff_clust_2cm->cd(counter_plots_1p5cm);}
	  TF1  *f2 = new TF1("f2","[0]*exp(-x/[1])",10,40);
	  f2->SetParameters(0,4000);
	  f2->SetParameters(1,10);
	  // h33->Fit("f2","R");
	  //   if(isLogarithm){
	  //   	gPad->SetLogy(1);
	  //   }
	  //pt_2cm_alpha->Draw("same");
	  //h33->Fit("expo");
	  h33->Draw("same");
	  
	  
	  TH1F *h12=(TH1F*)file->Get(Form("H-Ch%i_signal/hIntegN_ch%i",i,i));
	  if (h12==0x0) { continue; }
	  if(isNov2021TestBeam) {integ_2cm->cd(i-9);}
	  else if(isJuly2022TestBeam) {integ_2cm->cd(counter_plots_1p5cm);}
	  if(isRuns_90_10 || isRuns_80_20 || isRuns_85_15){
	    integral_2cm = integral_2cm + h12->GetMean();
	    rms_integral_2cm = rms_integral_2cm + h12->GetRMS();
	}
	  //integ->cd(3);
	  h12->Draw("same");
	  
	  TH1F *h16=(TH1F*)file->Get(Form("H-Ch%i_signal/hRms_ch%i",i,i));
	  if (h16==0x0) { continue; }
	  if(isNov2021TestBeam) {rms_2cm->cd(i-9);}
	  else if(isJuly2022TestBeam) {rms_2cm->cd(counter_plots_1p5cm);}
	  if(isRuns_90_10 || isRuns_80_20 || isRuns_85_15){
	    rms_2cm_var = rms_2cm_var + h16->GetMean();
	    rms_rms_2cm = rms_rms_2cm + h16->GetRMS();
	  }
	  h16->Draw("same");
	  
	  
	  TH1F *h20=(TH1F*)file->Get(Form("H-Ch%i_signal/hNPeaks_clust_ch%i",i,i));
	  if (h20==0x0) { continue; }
	  if(isNov2021TestBeam) {npeaks_clustser_2cm->cd(i-9);}
	  else if(isJuly2022TestBeam) {npeaks_clustser_2cm->cd(counter_plots_1p5cm);}
	  h20->GetXaxis()->SetRangeUser(0.,90.);
	  
	  if(isRuns_90_10 || isRuns_80_20 || isRuns_85_15){
	    mean_clusters_2cm = mean_clusters_2cm + h20->GetMean();
	    rms_clusters_2cm = rms_clusters_2cm + h20->GetRMS();
	  }
	  if(expected_cluster > 0.){
	    h20->Fit("gaus");
	  }
	  else{
	    TF1 *f10=new TF1("f10","[0]*TMath::Poisson(x,[1])",0,90.);                                                          
	    f10->SetParName(0,"Normalisation");
	    f10->SetParName(1,"#mu");
	    f10->SetParameters(0,1000);
	    f10->SetParameters(1,h20->GetMean());
	    h20->Fit("f10","R");
	  }
	  if(isLogarithm){
	    gPad->SetLogy(1);
	  }
	  h20->Draw("same");
	  TPaveText *pt_2cm_cluster =  new TPaveText(0.12,0.85,0.6,0.88,"NDC");
	  pt_2cm_cluster->SetTextSize(0.055);
	  pt_2cm_cluster->SetTextColor(kRed);
	  pt_2cm_cluster->SetFillColor(kYellow);
	  pt_2cm_cluster->SetTextAlign(12);
	  pt_2cm_cluster->AddText(Form("Expected Clusters: %.1f - Track Angle (deg) %0.1f",expected_cluster,alpha));
	  pt_2cm_cluster->Draw("same");
	  //pt_2cm_alpha->Draw("same");
	  
	  TH1F *h31=(TH1F*)file->Get(Form("H-Ch%i_signal/hNElectrons_per_cluster_ch%i",i,i));
	  if (h31==0x0) { continue; }
	  if(isNov2021TestBeam) {cluster_population_canvas_2cm->cd(i-9);}
	  else if(isJuly2022TestBeam) {cluster_population_canvas_2cm->cd(counter_plots_1p5cm);}
	  if(isRuns_90_10 || isRuns_80_20 || isRuns_85_15){
	    epc_2cm = epc_2cm+ h31->GetMean();
	    rms_epc_2cm = rms_epc_2cm + h31->GetRMS();
	  }
	  //h31->GetXaxis()->SetRangeUser(0.,5.);
	  //h31->Fit("expo");
	  h31->Draw("same");
	  TPaveText *cluster_population_2cm = new TPaveText(0.13,0.85,0.6,0.88,"NDC");
	  gPad->SetLogy(0);
	  cluster_population_2cm->SetTextSize(0.055);
	  cluster_population_2cm->SetTextColor(kRed);
	  cluster_population_2cm->SetTextAlign(12);
	  cluster_population_2cm->SetFillColor(kYellow);
	  cluster_population_2cm->AddText(Form("Expected Electrons per Cluster: %.1f - Track Angle (deg) %0.1f",cluster_population,alpha));
	  //pt_2cm_alpha->Draw("same");
	  cluster_population_2cm->Draw("same");
	  
	  
	  TH1F *h32=(TH1F*)file->Get(Form("H-Ch%i_signal/hTimeDifference_ch%i",i,i));
	  if (h32==0x0) { continue; }
	  if(isNov2021TestBeam) {timediff_2cm->cd(i-9);}
	  else if(isJuly2022TestBeam) {timediff_2cm->cd(counter_plots_1p5cm);}
	  TF1  *f3 = new TF1("f3","[p0]*exp(-x/[p1])",0,8);
	  //TF1  *f3 = new TF1("f3","[0]*exp(-x/[1])",0,12);
	  f3->SetParameters(0,4000.);
	  f3->SetParameters(1,3);
	  // h32->Fit("f3","R");
	  h32->Draw("same");
	  // if(isLogarithm){
	  //   	gPad->SetLogy(1);
	  // }
	  pt_2cm_alpha->Draw("same");
	  //h31->GetXaxis()->SetRangeUser(0.,5.);
	  //h31->Fit("expo");
	  counter_plots_1p5cm = counter_plots_1p5cm + 1;
	} // if on 2 cm or 1.5 cm tubes
      } // loop on all tubes 
      
      if(isRuns_90_10 || isRuns_80_20 || isRuns_85_15){ 

	aveph_2cm = aveph_2cm/(float) ((isNov2021TestBeam) > 0 ? Channel_2cm_size : Channel_1p5cm_size);
	rms_aveph_2cm = rms_aveph_2cm/(float) ((isNov2021TestBeam) > 0 ? Channel_2cm_size : Channel_1p5cm_size);
	
	epc_2cm = epc_2cm/(float) ((isNov2021TestBeam) > 0 ? Channel_2cm_size : Channel_1p5cm_size);
	rms_epc_2cm = rms_epc_2cm/(float) ((isNov2021TestBeam) > 0 ? Channel_2cm_size : Channel_1p5cm_size);
	
	integral_2cm = (integral_2cm/(float) ((isNov2021TestBeam) > 0 ? Channel_2cm_size : Channel_1p5cm_size));
	rms_integral_2cm = (rms_integral_2cm/(float) ((isNov2021TestBeam) > 0 ? Channel_2cm_size : Channel_1p5cm_size));
	
	bsl_2cm_var = bsl_2cm_var/(float) ((isNov2021TestBeam) > 0 ? Channel_2cm_size : Channel_1p5cm_size);
	rms_bsl_2cm = rms_bsl_2cm/(float) ((isNov2021TestBeam) > 0 ? Channel_2cm_size : Channel_1p5cm_size);
	
	rms_2cm_var = rms_2cm_var/(float) ((isNov2021TestBeam) > 0 ? Channel_2cm_size : Channel_1p5cm_size);
	rms_rms_2cm = rms_rms_2cm/(float) ((isNov2021TestBeam) > 0 ? Channel_2cm_size : Channel_1p5cm_size);
	
	mean_clusters_2cm = mean_clusters_2cm/(float) ((isNov2021TestBeam) > 0 ? Channel_2cm_size : Channel_1p5cm_size);
	rms_clusters_2cm = rms_clusters_2cm/(float) ((isNov2021TestBeam) > 0 ? Channel_2cm_size : Channel_1p5cm_size);
	
	cout <<"Mean" << mean_clusters_2cm <<"Track Angle"<< alpha << endl;
	mean_electrons_2cm = mean_electrons_2cm/(float) ((isNov2021TestBeam) > 0 ? Channel_2cm_size : Channel_1p5cm_size);
	rms_electrons_2cm = rms_electrons_2cm/(float) ((isNov2021TestBeam) > 0 ? Channel_2cm_size : Channel_1p5cm_size);
	
	maximum_2cm = maximum_2cm/(float) ((isNov2021TestBeam) > 0 ? Channel_2cm_size : Channel_1p5cm_size);
	rms_maximum_2cm = rms_maximum_2cm/(float) ((isNov2021TestBeam) > 0 ? Channel_2cm_size : Channel_1p5cm_size);

	//Normalization
	//cout <<"Expected Cluster" << expected_cluster <<"Track Angle"<< alpha << endl;
	mean_clusters_2cm = mean_clusters_2cm/expected_cluster;
	//cout <<"After division" << mean_clusters_2cm << endl <<"Track Angle"<< alpha << endl;
	mean_electrons_2cm = mean_electrons_2cm/expected_electrons;
	rms_clusters_2cm = rms_clusters_2cm/expected_cluster;
	rms_electrons_2cm = rms_electrons_2cm/expected_electrons;
	
	aveph_summary_2cm->cd();
	aveph_summary_2cm->SetGridx();
  	aveph_summary_2cm->SetGridy();
	if(isScanAngleStudy){
	  gr_aveph_summary_2cm->SetPoint(counter_filling_electrons_2cm,alpha,aveph_2cm);
	}
	else if(isScanHVStudy || isGasMixtureStudy){
	  gr_aveph_summary_2cm->SetPoint(counter_filling_electrons_2cm,counter_filling_electrons_2cm,aveph_2cm);
	}
	else if(isScanSamplingStudy){
	  gr_aveph_summary_2cm->SetPoint(counter_filling_electrons_2cm,_gsample,aveph_2cm);
	}
	gr_aveph_summary_2cm->SetPointError(counter_filling_electrons_2cm,0.,rms_aveph_2cm);
	gr_aveph_summary_2cm->SetMarkerStyle(21);
	gr_aveph_summary_2cm->SetMarkerSize(0.5);
	gr_aveph_summary_2cm->SetMarkerColor(kBlue);
	gr_aveph_summary_2cm->SetLineColor(kBlue);
	TLatex *latex_aveph_summary_2cm = new TLatex(gr_aveph_summary_2cm->GetX()[number_of_trials], gr_aveph_summary_2cm->GetY()[number_of_trials],Form("%s DERIV N1 %s N2 %s N3 %s N4 %s scale_cut %s",gasMixture.Data(),N1.Data(),N2.Data(),N3.Data(),N4.Data(),scale_cut.Data())); 
	latex_aveph_summary_2cm->SetTextSize(0.01); 
	latex_aveph_summary_2cm->SetTextColor(kRed);
	gr_aveph_summary_2cm->GetListOfFunctions()->Add(latex_aveph_summary_2cm);
	gr_aveph_summary_2cm->Draw("sameAP");
	
	rms_summary_2cm->cd();
	rms_summary_2cm->SetGridx();
  	rms_summary_2cm->SetGridy();
	if(isScanAngleStudy){
	  gr_rms_summary_2cm->SetPoint(counter_filling_electrons_2cm,alpha,rms_2cm_var);
	}
	else if(isScanHVStudy || isGasMixtureStudy){
	  gr_rms_summary_2cm->SetPoint(counter_filling_electrons_2cm,counter_filling_electrons_2cm,rms_2cm_var);
	}
	else if(isScanSamplingStudy){
	  gr_rms_summary_2cm->SetPoint(counter_filling_electrons_2cm,_gsample,rms_2cm_var);
	}
	gr_rms_summary_2cm->SetPointError(counter_filling_electrons_2cm,0.,rms_rms_2cm);
	gr_rms_summary_2cm->SetMarkerStyle(21);
	gr_rms_summary_2cm->SetMarkerSize(0.5);
	gr_rms_summary_2cm->SetMarkerColor(kBlue);
	gr_rms_summary_2cm->SetLineColor(kBlue);
	TLatex *latex_rms_summary_2cm = new TLatex(gr_rms_summary_2cm->GetX()[number_of_trials], gr_rms_summary_2cm->GetY()[number_of_trials],Form("%s DERIV N1 %s N2 %s N3 %s N4 %s scale_cut %s",gasMixture.Data(),N1.Data(),N2.Data(),N3.Data(),N4.Data(),scale_cut.Data())); 
	latex_rms_summary_2cm->SetTextSize(0.01); 
	latex_rms_summary_2cm->SetTextColor(kRed);
	gr_rms_summary_2cm->GetListOfFunctions()->Add(latex_rms_summary_2cm);
	gr_rms_summary_2cm->Draw("sameAP");
	
	epc_summary_2cm->cd();
	epc_summary_2cm->SetGridx();
  	epc_summary_2cm->SetGridy();
	if(isScanAngleStudy){
	  gr_epc_summary_2cm->SetPoint(counter_filling_electrons_2cm,alpha,epc_2cm);
	}
	else if(isScanHVStudy || isGasMixtureStudy){
	  gr_epc_summary_2cm->SetPoint(counter_filling_electrons_2cm,counter_filling_electrons_2cm,epc_2cm);
	}
	else if(isScanSamplingStudy){
	  gr_epc_summary_2cm->SetPoint(counter_filling_electrons_2cm,_gsample,epc_2cm);
	}
	gr_epc_summary_2cm->SetPointError(counter_filling_electrons_2cm,0.,rms_epc_2cm);
	gr_epc_summary_2cm->SetMarkerStyle(21);
	gr_epc_summary_2cm->SetMarkerSize(0.5);
	gr_epc_summary_2cm->SetMarkerColor(kBlue);
	gr_epc_summary_2cm->SetLineColor(kBlue);
	TLatex *latex_epc_summary_2cm = new TLatex(gr_epc_summary_2cm->GetX()[number_of_trials], gr_epc_summary_2cm->GetY()[number_of_trials],Form("%s DERIV N1 %s N2 %s N3 %s N4 %s scale_cut %s",gasMixture.Data(),N1.Data(),N2.Data(),N3.Data(),N4.Data(),scale_cut.Data())); 
	latex_epc_summary_2cm->SetTextSize(0.01); 
	latex_epc_summary_2cm->SetTextColor(kRed);
	gr_epc_summary_2cm->GetListOfFunctions()->Add(latex_epc_summary_2cm);
	gr_epc_summary_2cm->Draw("sameAP");
	
	// TLegend *leg= new TLegend(0.5,0.75,0.85,0.85); 
	// leg->AddEntry(tmpsignal_1.back(),"Waveform");
	
	bsl_summary_2cm->cd();
	bsl_summary_2cm->SetGridx();
  	bsl_summary_2cm->SetGridy();
	if(isScanAngleStudy){
	  gr_bsl_summary_2cm->SetPoint(counter_filling_electrons_2cm,alpha,bsl_2cm_var);
	}
	else if(isScanHVStudy || isGasMixtureStudy){
	  gr_bsl_summary_2cm->SetPoint(counter_filling_electrons_2cm,counter_filling_electrons_2cm,bsl_2cm_var);
	}
	else if(isScanSamplingStudy){
	  gr_bsl_summary_2cm->SetPoint(counter_filling_electrons_2cm,_gsample,bsl_2cm_var);
	}
	gr_bsl_summary_2cm->SetPointError(counter_filling_electrons_2cm,0.,rms_bsl_2cm);
	gr_bsl_summary_2cm->SetMarkerStyle(21);
	gr_bsl_summary_2cm->SetMarkerSize(0.5);
	gr_bsl_summary_2cm->SetMarkerColor(kBlue);
	gr_bsl_summary_2cm->SetLineColor(kBlue);
	TLatex *latex_bsl_summary_2cm = new TLatex(gr_bsl_summary_2cm->GetX()[number_of_trials], gr_bsl_summary_2cm->GetY()[number_of_trials],Form("%s DERIV N1 %s N2 %s N3 %s N4 %s scale_cut %s",gasMixture.Data(),N1.Data(),N2.Data(),N3.Data(),N4.Data(),scale_cut.Data())); 
	latex_bsl_summary_2cm->SetTextSize(0.01); 
	latex_bsl_summary_2cm->SetTextColor(kRed);
	gr_bsl_summary_2cm->GetListOfFunctions()->Add(latex_bsl_summary_2cm);
	gr_bsl_summary_2cm->Draw("sameAP");
	
	maximum_summary_2cm->cd();
	maximum_summary_2cm->SetGridx();
  	maximum_summary_2cm->SetGridy();
	if(isScanAngleStudy){
	  gr_maximum_summary_2cm->SetPoint(counter_filling_electrons_2cm,alpha,maximum_2cm);
	}
	else if(isScanHVStudy || isGasMixtureStudy){
	  gr_maximum_summary_2cm->SetPoint(counter_filling_electrons_2cm,counter_filling_electrons_2cm,maximum_2cm);
	}
	else if(isScanSamplingStudy){
	  gr_maximum_summary_2cm->SetPoint(counter_filling_electrons_2cm,_gsample,maximum_2cm);
	}
	gr_maximum_summary_2cm->SetPointError(counter_filling_electrons_2cm,0.,rms_maximum_2cm);
	gr_maximum_summary_2cm->SetMarkerStyle(21);
	gr_maximum_summary_2cm->SetMarkerSize(0.5);
	gr_maximum_summary_2cm->SetMarkerColor(kBlue);
	gr_maximum_summary_2cm->SetLineColor(kBlue);
	TLatex *latex_maximum_summary_2cm = new TLatex(gr_maximum_summary_2cm->GetX()[number_of_trials], gr_maximum_summary_2cm->GetY()[number_of_trials],Form("%s DERIV N1 %s N2 %s N3 %s N4 %s scale_cut %s",gasMixture.Data(),N1.Data(),N2.Data(),N3.Data(),N4.Data(),scale_cut.Data())); 
	latex_maximum_summary_2cm->SetTextSize(0.01); 
	latex_maximum_summary_2cm->SetTextColor(kRed);
	gr_maximum_summary_2cm->GetListOfFunctions()->Add(latex_maximum_summary_2cm);
	gr_maximum_summary_2cm->Draw("sameAP");
	
	integral_summary_2cm->cd();
	integral_summary_2cm->SetGridx();
  	integral_summary_2cm->SetGridy();
	if(isScanAngleStudy){
	  gr_integral_summary_2cm->SetPoint(counter_filling_electrons_2cm,alpha,integral_2cm);
	}
	else if(isScanHVStudy || isGasMixtureStudy){
	  gr_integral_summary_2cm->SetPoint(counter_filling_electrons_2cm,counter_filling_electrons_2cm,integral_2cm);
	}
	else if(isScanSamplingStudy){
	  gr_integral_summary_2cm->SetPoint(counter_filling_electrons_2cm,_gsample,integral_2cm);
	}
	gr_integral_summary_2cm->SetPointError(counter_filling_electrons_2cm,0.,rms_integral_2cm);
	gr_integral_summary_2cm->SetMarkerStyle(21);
	gr_integral_summary_2cm->SetMarkerSize(0.5);
	gr_integral_summary_2cm->SetMarkerColor(kBlue);
	gr_integral_summary_2cm->SetLineColor(kBlue);
	TLatex *latex_integral_summary_2cm = new TLatex(gr_integral_summary_2cm->GetX()[number_of_trials], gr_integral_summary_2cm->GetY()[number_of_trials],Form("%s DERIV N1 %s N2 %s N3 %s N4 %s scale_cut %s",gasMixture.Data(),N1.Data(),N2.Data(),N3.Data(),N4.Data(),scale_cut.Data())); 
	latex_integral_summary_2cm->SetTextSize(0.01); 
	latex_integral_summary_2cm->SetTextColor(kRed);
	gr_integral_summary_2cm->GetListOfFunctions()->Add(latex_integral_summary_2cm);
	gr_integral_summary_2cm->Draw("sameAP");
	
	efficiency_electrons_2cm->cd();
	efficiency_electrons_2cm->SetGridx();
  	efficiency_electrons_2cm->SetGridy();
	if(isScanAngleStudy){
	  gr_efficiency_electrons_2cm->SetPoint(counter_filling_electrons_2cm,alpha,mean_electrons_2cm);
	}
	else if(isScanHVStudy || isGasMixtureStudy){
	  gr_efficiency_electrons_2cm->SetPoint(counter_filling_electrons_2cm,counter_filling_electrons_2cm,mean_electrons_2cm);
	}
	else if(isScanSamplingStudy){
	  gr_efficiency_electrons_2cm->SetPoint(counter_filling_electrons_2cm,_gsample,mean_electrons_2cm);
	}    
	gr_efficiency_electrons_2cm->SetPointError(counter_filling_electrons_2cm,0.,rms_electrons_2cm);
	gr_efficiency_electrons_2cm->SetMarkerStyle(21);
	gr_efficiency_electrons_2cm->SetMarkerSize(0.5);
	gr_efficiency_electrons_2cm->SetMarkerColor(kBlue);
	gr_efficiency_electrons_2cm->SetLineColor(kBlue);
	TLatex *latex_efficiency_electrons_2cm = new TLatex(gr_efficiency_electrons_2cm->GetX()[number_of_trials], gr_efficiency_electrons_2cm->GetY()[number_of_trials],Form("%s DERIV N1 %s N2 %s N3 %s N4 %s scale_cut %s",gasMixture.Data(),N1.Data(),N2.Data(),N3.Data(),N4.Data(),scale_cut.Data())); 
	latex_efficiency_electrons_2cm->SetTextSize(0.01); 
	latex_efficiency_electrons_2cm->SetTextColor(kRed);
	gr_efficiency_electrons_2cm->GetListOfFunctions()->Add(latex_efficiency_electrons_2cm);
	gr_efficiency_electrons_2cm->Draw("sameAP");
	
	efficiency_clusters_2cm->cd();
	efficiency_clusters_2cm->SetGridx();
  	efficiency_clusters_2cm->SetGridy();
	if(isScanAngleStudy){
	  gr_efficiency_clusters_2cm->SetPoint(counter_filling_clusters_2cm,alpha,mean_clusters_2cm);
	}
	else if(isScanHVStudy || isGasMixtureStudy){
	  gr_efficiency_clusters_2cm->SetPoint(counter_filling_clusters_2cm,counter_filling_clusters_2cm,mean_clusters_2cm);
	}
	else if(isScanSamplingStudy){
	  gr_efficiency_clusters_2cm->SetPoint(counter_filling_clusters_2cm,_gsample,mean_clusters_2cm);
	}    
	gr_efficiency_clusters_2cm->SetPointError(counter_filling_clusters_2cm,0.,rms_clusters_2cm);
	gr_efficiency_clusters_2cm->SetMarkerStyle(21);
	gr_efficiency_clusters_2cm->SetMarkerSize(0.5);
	gr_efficiency_clusters_2cm->SetMarkerColor(kBlue);
	gr_efficiency_clusters_2cm->SetLineColor(kBlue);
	TLatex *latex_efficiency_clusters_2cm = new TLatex(gr_efficiency_clusters_2cm->GetX()[number_of_trials], gr_efficiency_clusters_2cm->GetY()[number_of_trials],Form("%s DERIV N1 %s N2 %s N3 %s N4 %s scale_cut %s",gasMixture.Data(),N1.Data(),N2.Data(),N3.Data(),N4.Data(),scale_cut.Data())); 
	latex_efficiency_clusters_2cm->SetTextSize(0.01); 
	latex_efficiency_clusters_2cm->SetTextColor(kRed);
	gr_efficiency_clusters_2cm->GetListOfFunctions()->Add(latex_efficiency_clusters_2cm);
	gr_efficiency_clusters_2cm->Draw("sameAP");
	
	

	counter_filling_clusters_2cm = counter_filling_clusters_2cm +1;
	counter_filling_electrons_2cm = counter_filling_electrons_2cm +1;
	
	
	if((name_file_compact == Runs_90_10[Runs_90_10_size - 1]) || (name_file_compact == Runs_80_20[Runs_80_20_size - 1])){
	  
	  efficiency_clusters_2cm->SaveAs(Form("/lustrehome/bdanzi/TestBeam2022/analysistestbeam2022/analysis_2021/drifttubes_offline_analysis/testbeam_analysis/Plots/efficiency_clusters_%scm.pdf",tubes.Data()));
	  efficiency_electrons_2cm->SaveAs(Form("/lustrehome/bdanzi/TestBeam2022/analysistestbeam2022/analysis_2021/drifttubes_offline_analysis/testbeam_analysis/Plots/efficiency_electrons_%scm.pdf",tubes.Data())); 
	  integral_summary_2cm->SaveAs(Form("/lustrehome/bdanzi/TestBeam2022/analysistestbeam2022/analysis_2021/drifttubes_offline_analysis/testbeam_analysis/Plots/integral_summary_%scm.pdf",tubes.Data()));
	  maximum_summary_2cm->SaveAs(Form("/lustrehome/bdanzi/TestBeam2022/analysistestbeam2022/analysis_2021/drifttubes_offline_analysis/testbeam_analysis/Plots/maximum_summary_%scm.pdf",tubes.Data()));
	  epc_summary_2cm->SaveAs(Form("/lustrehome/bdanzi/TestBeam2022/analysistestbeam2022/analysis_2021/drifttubes_offline_analysis/testbeam_analysis/Plots/epc_summary_%scm.pdf",tubes.Data()));
	  bsl_summary_2cm->SaveAs(Form("/lustrehome/bdanzi/TestBeam2022/analysistestbeam2022/analysis_2021/drifttubes_offline_analysis/testbeam_analysis/Plots/bsl_summary_%scm.pdf",tubes.Data()));
	  rms_summary_2cm->SaveAs(Form("/lustrehome/bdanzi/TestBeam2022/analysistestbeam2022/analysis_2021/drifttubes_offline_analysis/testbeam_analysis/Plots/rms_summary_%scm.pdf",tubes.Data()));
	  aveph_summary_2cm->SaveAs(Form("/lustrehome/bdanzi/TestBeam2022/analysistestbeam2022/analysis_2021/drifttubes_offline_analysis/testbeam_analysis/Plots/aveph_summary_%scm.pdf",tubes.Data()));
	  
	  efficiency_clusters_2cm->SaveAs(Form("/lustrehome/bdanzi/TestBeam2022/analysistestbeam2022/analysis_2021/drifttubes_offline_analysis/testbeam_analysis/Plots/efficiency_clusters_%scm.png",tubes.Data()));
	  efficiency_electrons_2cm->SaveAs(Form("/lustrehome/bdanzi/TestBeam2022/analysistestbeam2022/analysis_2021/drifttubes_offline_analysis/testbeam_analysis/Plots/efficiency_electrons_%scm.png",tubes.Data())); 
	  integral_summary_2cm->SaveAs(Form("/lustrehome/bdanzi/TestBeam2022/analysistestbeam2022/analysis_2021/drifttubes_offline_analysis/testbeam_analysis/Plots/integral_summary_%scm.png",tubes.Data()));
	  maximum_summary_2cm->SaveAs(Form("/lustrehome/bdanzi/TestBeam2022/analysistestbeam2022/analysis_2021/drifttubes_offline_analysis/testbeam_analysis/Plots/maximum_summary_%scm.png",tubes.Data()));
	  epc_summary_2cm->SaveAs(Form("/lustrehome/bdanzi/TestBeam2022/analysistestbeam2022/analysis_2021/drifttubes_offline_analysis/testbeam_analysis/Plots/epc_summary_%scm.png",tubes.Data()));
	  bsl_summary_2cm->SaveAs(Form("/lustrehome/bdanzi/TestBeam2022/analysistestbeam2022/analysis_2021/drifttubes_offline_analysis/testbeam_analysis/Plots/bsl_summary_%scm.png",tubes.Data()));
	  rms_summary_2cm->SaveAs(Form("/lustrehome/bdanzi/TestBeam2022/analysistestbeam2022/analysis_2021/drifttubes_offline_analysis/testbeam_analysis/Plots/rms_summary_%scm.png",tubes.Data()));
	  aveph_summary_2cm->SaveAs(Form("/lustrehome/bdanzi/TestBeam2022/analysistestbeam2022/analysis_2021/drifttubes_offline_analysis/testbeam_analysis/Plots/aveph_summary_%scm.png",tubes.Data()));
	  
	  
	  
	}
      }
      
      
      
      
      
      if(savePlots){
	timediff_clust_2cm->SaveAs(Form("/lustrehome/bdanzi/TestBeam2022/analysistestbeam2022/analysis_2021/drifttubes_offline_analysis/testbeam_analysis/Plots/%s/clust_difference_%scm.pdf",fname.Data(),tubes.Data()));
	timediff_2cm->SaveAs(Form("/lustrehome/bdanzi/TestBeam2022/analysistestbeam2022/analysis_2021/drifttubes_offline_analysis/testbeam_analysis/Plots/%s/electrons_difference_%scm.pdf",fname.Data(),tubes.Data()));
	npeaks_2cm->SaveAs(Form("/lustrehome/bdanzi/TestBeam2022/analysistestbeam2022/analysis_2021/drifttubes_offline_analysis/testbeam_analysis/Plots/%s/nelectrons_%scm.pdf",fname.Data(),tubes.Data()));
	//tpeaks_2cm->SaveAs(Form("/lustrehome/bdanzi/TestBeam2022/analysistestbeam2022/analysis_2021/drifttubes_offline_analysis/testbeam_analysis/Plots/%s/tpeaks_%scm.pdf",fname.Data(),tubes.Data()));
	tfpeaks_2cm->SaveAs(Form("/lustrehome/bdanzi/TestBeam2022/analysistestbeam2022/analysis_2021/drifttubes_offline_analysis/testbeam_analysis/Plots/%s/tfpeaks_%scm.pdf",fname.Data(),tubes.Data()));
	tlpeaks_2cm->SaveAs(Form("/lustrehome/bdanzi/TestBeam2022/analysistestbeam2022/analysis_2021/drifttubes_offline_analysis/testbeam_analysis/Plots/%s/tlpeaks_%scm.pdf",fname.Data(),tubes.Data()));
	//hnpeaks_2cm->SaveAs(Form("/lustrehome/bdanzi/TestBeam2022/analysistestbeam2022/analysis_2021/drifttubes_offline_analysis/testbeam_analysis/Plots/%s/hnpeaks_%scm.pdf",fname.Data()));
	hncluster_2cm->SaveAs(Form("/lustrehome/bdanzi/TestBeam2022/analysistestbeam2022/analysis_2021/drifttubes_offline_analysis/testbeam_analysis/Plots/%s/hncluster_%scm.pdf",fname.Data(),tubes.Data()));
	npeaks_clustser_2cm->SaveAs(Form("/lustrehome/bdanzi/TestBeam2022/analysistestbeam2022/analysis_2021/drifttubes_offline_analysis/testbeam_analysis/Plots/%s/nclusters_%scm.pdf",fname.Data(),tubes.Data()));
	hpeaks_2cm->SaveAs(Form("/lustrehome/bdanzi/TestBeam2022/analysistestbeam2022/analysis_2021/drifttubes_offline_analysis/testbeam_analysis/Plots/%s/hpeaks_%scm.pdf",fname.Data(),tubes.Data()));
	bsl_2cm->SaveAs(Form("/lustrehome/bdanzi/TestBeam2022/analysistestbeam2022/analysis_2021/drifttubes_offline_analysis/testbeam_analysis/Plots/%s/bsl_%scm.pdf",fname.Data(),tubes.Data()));
	max_2cm->SaveAs(Form("/lustrehome/bdanzi/TestBeam2022/analysistestbeam2022/analysis_2021/drifttubes_offline_analysis/testbeam_analysis/Plots/%s/max_%scm.pdf",fname.Data(),tubes.Data()));
	integ_2cm->SaveAs(Form("/lustrehome/bdanzi/TestBeam2022/analysistestbeam2022/analysis_2021/drifttubes_offline_analysis/testbeam_analysis/Plots/%s/integ_%scm.pdf",fname.Data(),tubes.Data()));
	rms_2cm->SaveAs(Form("/lustrehome/bdanzi/TestBeam2022/analysistestbeam2022/analysis_2021/drifttubes_offline_analysis/testbeam_analysis/Plots/%s/rms_%scm.pdf",fname.Data(),tubes.Data()));
	cluster_population_canvas_2cm->SaveAs(Form("/lustrehome/bdanzi/TestBeam2022/analysistestbeam2022/analysis_2021/drifttubes_offline_analysis/testbeam_analysis/Plots/%s/cluster_population_%scm.pdf",fname.Data(),tubes.Data()));
	//min->SaveAs(Form("/lustrehome/bdanzi/offline_analysis/testbeam_analysis/Plots/%s/min_ch%i.pdf",fname.Data(),i));
	hnelectron_2cm->SaveAs(Form("/lustrehome/bdanzi/TestBeam2022/analysistestbeam2022/analysis_2021/drifttubes_offline_analysis/testbeam_analysis/Plots/%s/hnelectron_%scm.pdf",fname.Data(),tubes.Data()));
	
	timediff_clust_2cm->SaveAs(Form("/lustrehome/bdanzi/TestBeam2022/analysistestbeam2022/analysis_2021/drifttubes_offline_analysis/testbeam_analysis/Plots/%s/clust_difference_%scm.png",fname.Data(),tubes.Data()));
	timediff_2cm->SaveAs(Form("/lustrehome/bdanzi/TestBeam2022/analysistestbeam2022/analysis_2021/drifttubes_offline_analysis/testbeam_analysis/Plots/%s/electrons_difference_%scm.png",fname.Data(),tubes.Data()));
	npeaks_2cm->SaveAs(Form("/lustrehome/bdanzi/TestBeam2022/analysistestbeam2022/analysis_2021/drifttubes_offline_analysis/testbeam_analysis/Plots/%s/nelectrons_%scm.png",fname.Data(),tubes.Data()));
	//tpeaks_2cm->SaveAs(Form("/lustrehome/bdanzi/TestBeam2022/analysistestbeam2022/analysis_2021/drifttubes_offline_analysis/testbeam_analysis/Plots/%s/tpeaks_%scm.pdf",fname.Data(),tubes.Data()));
	tfpeaks_2cm->SaveAs(Form("/lustrehome/bdanzi/TestBeam2022/analysistestbeam2022/analysis_2021/drifttubes_offline_analysis/testbeam_analysis/Plots/%s/tfpeaks_%scm.png",fname.Data(),tubes.Data()));
	tlpeaks_2cm->SaveAs(Form("/lustrehome/bdanzi/TestBeam2022/analysistestbeam2022/analysis_2021/drifttubes_offline_analysis/testbeam_analysis/Plots/%s/tlpeaks_%scm.png",fname.Data(),tubes.Data()));
	//hnpeaks_2cm->SaveAs(Form("/lustrehome/bdanzi/TestBeam2022/analysistestbeam2022/analysis_2021/drifttubes_offline_analysis/testbeam_analysis/Plots/%s/hnpeaks_%scm.pdf",fname.Data()));
	hncluster_2cm->SaveAs(Form("/lustrehome/bdanzi/TestBeam2022/analysistestbeam2022/analysis_2021/drifttubes_offline_analysis/testbeam_analysis/Plots/%s/hncluster_%scm.png",fname.Data(),tubes.Data()));
	npeaks_clustser_2cm->SaveAs(Form("/lustrehome/bdanzi/TestBeam2022/analysistestbeam2022/analysis_2021/drifttubes_offline_analysis/testbeam_analysis/Plots/%s/nclusters_%scm.png",fname.Data(),tubes.Data()));
	hpeaks_2cm->SaveAs(Form("/lustrehome/bdanzi/TestBeam2022/analysistestbeam2022/analysis_2021/drifttubes_offline_analysis/testbeam_analysis/Plots/%s/hpeaks_%scm.png",fname.Data(),tubes.Data()));
	bsl_2cm->SaveAs(Form("/lustrehome/bdanzi/TestBeam2022/analysistestbeam2022/analysis_2021/drifttubes_offline_analysis/testbeam_analysis/Plots/%s/bsl_%scm.png",fname.Data(),tubes.Data()));
	max_2cm->SaveAs(Form("/lustrehome/bdanzi/TestBeam2022/analysistestbeam2022/analysis_2021/drifttubes_offline_analysis/testbeam_analysis/Plots/%s/max_%scm.png",fname.Data(),tubes.Data()));
	integ_2cm->SaveAs(Form("/lustrehome/bdanzi/TestBeam2022/analysistestbeam2022/analysis_2021/drifttubes_offline_analysis/testbeam_analysis/Plots/%s/integ_%scm.png",fname.Data(),tubes.Data()));
	rms_2cm->SaveAs(Form("/lustrehome/bdanzi/TestBeam2022/analysistestbeam2022/analysis_2021/drifttubes_offline_analysis/testbeam_analysis/Plots/%s/rms_%scm.png",fname.Data(),tubes.Data()));
	cluster_population_canvas_2cm->SaveAs(Form("/lustrehome/bdanzi/TestBeam2022/analysistestbeam2022/analysis_2021/drifttubes_offline_analysis/testbeam_analysis/Plots/%s/cluster_population_%scm.png",fname.Data(),tubes.Data()));
	//min->SaveAs(Form("/lustrehome/bdanzi/offline_analysis/testbeam_analysis/Plots/%s/min_ch%i.pdf",fname.Data(),i));
	hnelectron_2cm->SaveAs(Form("/lustrehome/bdanzi/TestBeam2022/analysistestbeam2022/analysis_2021/drifttubes_offline_analysis/testbeam_analysis/Plots/%s/hnelectron_%scm.png",fname.Data(),tubes.Data()));
      }
      number_of_trials = number_of_trials + 1;  
    }
    
    
    
    
    
    fclose(fp);
  }
  
  
}
