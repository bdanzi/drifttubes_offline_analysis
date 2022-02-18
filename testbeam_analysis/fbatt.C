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


void fbatt(TString fname=""){

	TFile *file = new TFile(fname.Data(),"read");
	TH1F *h1=(TH1F*)file->Get(Form("f_batt"));
	h1->Draw();
	Double_t par[6];
	TF1 *g1    = new TF1("g1","gaus",-50e+3,0);
//	g1->SetParameter(0,1e+3);
//	g1->SetParameter(1,1e+4);
//	g1->SetParameter(1,4e+4);

	TF1 *g2    = new TF1("g2","gaus",0,50e+3);
//	g2->SetParameter(0,1e+3);
//	g2->SetParameter(1,478);
//	g2->SetParameter(2,4e+4);

	TF1 *total = new TF1("total","gaus(0)+gaus(3)",-50e+3,50e+3);
	total->SetLineColor(2);
	h1->Fit(g1,"R","0");
	h1->Fit(g2,"R+","0");
	g1->GetParameters(&par[0]);
	g2->GetParameters(&par[3]);
	total->SetParameters(par);

	h1->Fit(total,"R+");
	h1->Fit(total,"","",-100e+3,100e+3);
	cout<<" mean "<<total->GetParameter(1)<< " sigma "<< total->GetParameter(2)<<endl;


}
