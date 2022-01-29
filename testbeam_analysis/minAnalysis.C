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


void minAnalysis(TString fname=""){

	TFile *file = new TFile(fname.Data(),"read");
	TH1F *h1=(TH1F*)file->Get(Form("hmin_value"));
	h1->Draw();
	TF1 *gausFit = new TF1("gausFit", "gaus",-0.1,0);
	gausFit->SetNpx(1000);

	int ItgTot,itg,binx2,binx1;
	ItgTot=h1->Integral();
	binx2=h1->GetXaxis()->FindBin(0.);
	int bintot=100;
	float meanGauss,sgmGauss,thrMin=0.;
	for(int i=0;i<bintot;++i){
		if(h1->Integral(i,binx2)>=ItgTot*0.9){
			//cout<<" valore i="<<i<<endl;
			binx1=i;
		}
	}
	cout<< " bin 1 "<<binx1 << " bin 2 "<< binx2<<" x value of bin 1 = "<<h1->GetXaxis()->GetBinCenter(binx1)<<endl;
	h1->Fit("gausFit","","",h1->GetXaxis()->GetBinCenter(binx1),0.0);
	meanGauss=gausFit->GetParameter(1);
	sgmGauss=gausFit->GetParameter(2);
	thrMin=meanGauss-3*sgmGauss;
	cout<< "MEDIA GAUSS = "<< meanGauss << " SIGMA GAUSS = "<< sgmGauss<<" Soglia "<< thrMin <<endl;

}
