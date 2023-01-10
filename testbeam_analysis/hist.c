/*macro per la sovrapposizione di istogrammi*/
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

void hist(TString fname="/mnt/c/dati_root/histosOSC_run-00033.root"){

	gStyle->SetOptFit(1);
	gStyle->SetStatX(1);

	TFile *file = new TFile(fname.Data(),"read");
	//file->ls();
	/* TF1 *fGaus = new TF1("fGaus", "gaus"); */
	/* fGaus->SetNpx(1000); */
	/* TF1 *fGaus_2 = new TF1("fGaus_2", "gaus",-0.02,-0.002); */
	/* fGaus_2->SetNpx(1000); */

	/* //crystal */
	TF1 *f1  = new TF1("f1","crystalball",-0.03,0);
	//       	f1->SetParameters(1, -0.007, 0.001, 1e-2, 1);
	f1->SetParLimits(0,0,1000);
	f1->SetParLimits(1,-0.03,0);//1,-0.03,0
	f1->SetParLimits(2,0,1);//2,0,1
	f1->SetParLimits(3,0,10);//3,0,10
	f1->SetParLimits(4,0,1000);//4,0,1000
	f1->SetLineColor(kRed);
	f1->SetNpx(1000);
	TF1 *f2  = new TF1("f2","crystalball",-0.03,0);
	//	f2->SetParameters(1, -0.008, 0.002, 1e-2, 1);
	f2->SetParLimits(0,0,1000);
	f2->SetParLimits(1,-0.03,0);
	f2->SetParLimits(2,0,1);
	f2->SetParLimits(3,0,10);
	f2->SetParLimits(4,0,100);
	f2->SetLineColor(kYellow);
	f2->SetNpx(1000);

	/* TF1 *f2  = new TF1("f2","fGaus_2(0)+f1(3)",-0.02,-0.06); */
	/* f2->SetLineColor(kMagenta); */


	TCanvas *c1=new TCanvas("c1","c1");
	c1->cd();

	Int_t nPresCh=0;
	Float_t chn[20]={0};
	Float_t meanMinFlt[20]={0};
	Float_t sgmMinFlt[20]={0};
	Float_t meanMinOrg[20]={0};
	Float_t sgmMinOrg[20]={0};
	Float_t ratioFlt[20]={0};
	Float_t ratioOrg[20]={0};
	Float_t ratioMean[20]={0};
	Float_t ratioSigma[20]={0};

	for(int i =0; i<=23; ++i){
		//file->cd("H-Ch0_signal");
		TH1F *h1=(TH1F*)file->Get(Form("H-Ch%i_signal/hMinVNInR_ch%i",i,i));
		if (h1==0x0) { continue; }
		TH1F *h2=(TH1F*)file->Get(Form("H-Ch%i_signal/hMinVNInRoriginalW_ch%i",i,i));
		if (h2==0x0) { continue; }
		h1->SetLineColor(kBlue);
		//	h1->Rebin();
		h2->SetLineColor(kGreen+2);
		//	h2->Rebin();
		h1->GetXaxis()->SetTitle("Min Value[V]");
		h1->GetYaxis()->SetTitle("Entries");
		h1->Draw();
		f1->SetParameters(30,h1->GetMean(),h1->GetStdDev(),1,10);

		h1->Fit(f1,"","",-0.03,-0.001);//,-0.012,-0.004);	//fit con cb sulla forma d'onda originale
		gPad->Update();
		TPaveStats *stat = (TPaveStats*)h1->GetListOfFunctions()->FindObject("stats");
		stat->SetY1NDC(0.7);    //set new y start position
		stat->SetY2NDC(0.9);    //set new y end position
		stat->SetX1NDC(0.1);    //set new y start position
		stat->SetX2NDC(0.4);//set new y end position
		stat->SetTextColor(kBlue);
		h2->Draw("sames");
		// h1->SetStats(111);
		//	h2->Fit(fGaus_2,"R","",-0.018,-0.006);
		std::cout << " sigma_h2 " << h2->GetStdDev() << std::endl;
		std::cout << " mean_h2 " << h2->GetMean() << std::endl;

		f2->SetParameters(30,h2->GetMean(),h2->GetStdDev(),1,10);
		h2->Fit(f2,"","",-0.03,-0.001);//,-0.018,-0.006);
		h2->Fit(f2,"","",-0.03,-0.001);//,-0.018,-0.006);
		gPad->Draw();
		TPaveStats *stat2 = (TPaveStats*)h2->GetListOfFunctions()->FindObject("stats");
		stat2->SetY1NDC(0.5);    //set new y start position
		stat2->SetY2NDC(0.7);    //set new y end position
		stat2->SetX1NDC(0.1);    //set new x start position
		stat2->SetX2NDC(0.4);    //set new x end position
		stat2->SetTextColor(kGreen+2);


		TLegend *leg = new TLegend(0.45, 0.9, 0.6, 0.8);
		leg->SetFillColor(0);
		//leg->SetHeader("Histogram");
		leg->AddEntry(h1, "Filtered signal", "l");
		leg->AddEntry(h2, "Original signal", "l");
		/* leg->SetLegendTextSize(0.); */
		/* leg->SetLegendBorderSize(1); */
		/* leg->SetLegendFillColor(0); */
		/* leg->SetLegendFont(42); */
		leg->Draw();
		c1->Print(Form("/home/federica/eclipse-workspace/test_beam/file_root/image_minvalue/H%i_minvalue.png",i), "png");

		chn[nPresCh]=i;
		meanMinFlt[nPresCh]=f1/*fGaus*/->GetParameter(1);
		cout << " media non  filtrata " << meanMinFlt[nPresCh] << endl;
		sgmMinFlt[nPresCh]=f1/*fGaus*/->GetParameter(2);
		cout << " sigma non filtrata " << sgmMinFlt[nPresCh] << endl;
		meanMinOrg[nPresCh]=f2/*fGaus_2*/->GetParameter(1);
		sgmMinOrg[nPresCh]=f2/*fGaus_2*/->GetParameter(2);
		ratioFlt[nPresCh]=meanMinFlt[nPresCh]/sgmMinFlt[nPresCh];
		ratioOrg[nPresCh]=meanMinOrg[nPresCh]/sgmMinOrg[nPresCh];
		ratioMean[nPresCh]=meanMinFlt[nPresCh]/meanMinOrg[nPresCh];
		ratioSigma[nPresCh]=sgmMinFlt[nPresCh]/sgmMinOrg[nPresCh];
		++nPresCh;

	}



	//stima dell'rms.
	TCanvas *cRms=new TCanvas("cRms","cRms");
	cRms->cd();
	int ch=0;
	Float_t rms[23]={0};
	Float_t SgmRms[23]={0};
	Float_t err[23]={0};
	Float_t errNoFlt[23]={0};
	for(int j =0; j<=23; ++j){
		TH1F *h3=(TH1F*)file->Get(Form("H-Ch%i_signal/hRms_ch%i",j,j));
		if(h3==0x0){continue;}
		h3->Draw();
		chn[ch]=j;
		rms[ch]=h3->GetMean();
		SgmRms[ch]=h3->GetStdDev();
		cout<< rms[ch]<<endl;
		err[ch]=TMath::Abs(meanMinFlt[ch]/rms[ch]);
		errNoFlt[ch]=TMath::Abs(meanMinOrg[ch]/rms[ch]);
		++ch;
	}


	TCanvas *c2= new TCanvas ("c2","c2");
	c2->Divide(1,2);
	TCanvas *ca= new TCanvas ("ca","c2.0");
	ca->Divide(1,2);

	TCanvas *c3= new TCanvas("c3","c3");
	c3->Divide(1,2);
	TCanvas *c4 =new TCanvas("c4","c4");
	c4->Divide(1,2);

	c2->cd(1);
	TGraph *grMeanMinFlt= new TGraph(nPresCh,chn,meanMinFlt);
	grMeanMinFlt->SetTitle("Mean MinValue for Filter Wave");
	grMeanMinFlt->GetXaxis()->SetTitle("Channel");
	grMeanMinFlt->GetXaxis()->SetTitleSize(0.05);
	grMeanMinFlt->GetXaxis()->SetTitleOffset(1.);
	grMeanMinFlt->GetXaxis()->SetLabelSize(0.05); 	
	grMeanMinFlt->GetYaxis()->SetTitle("Mean_Min_Flt [V]");
	grMeanMinFlt->GetYaxis()->SetTitleOffset(1.);
	grMeanMinFlt->GetYaxis()->SetTitleSize(0.05);
	grMeanMinFlt->GetYaxis()->SetLabelSize(0.05); 
	grMeanMinFlt->Draw("AP");
	grMeanMinFlt->SetMarkerStyle(8);
	c2->cd(2);
	TGraph *grsgmMinFlt= new TGraph(nPresCh,chn,sgmMinFlt);
	grsgmMinFlt->SetTitle("Sigma MinValue for Filter Wave");
	grsgmMinFlt->GetXaxis()->SetTitle("Channel");
	grsgmMinFlt->GetXaxis()->SetTitleSize(0.05);
	grsgmMinFlt->GetXaxis()->SetTitleOffset(1.);
	grsgmMinFlt->GetXaxis()->SetLabelSize(0.05); 	
	grsgmMinFlt->GetYaxis()->SetTitle("Sigma_Min_Flt [V]");
	grsgmMinFlt->GetYaxis()->SetTitleOffset(1.); 
	grsgmMinFlt->GetYaxis()->SetTitleSize(0.05);
	grsgmMinFlt->GetYaxis()->SetLabelSize(0.05);
	grsgmMinFlt->Draw("AP");
	grsgmMinFlt->SetMarkerStyle(8);

	ca->cd(1);
	TGraph *grMeanMinOrg= new TGraph(nPresCh,chn,meanMinOrg);
	grMeanMinOrg->SetTitle("Mean MinValue for Original Wave");
	grMeanMinOrg->GetXaxis()->SetTitle("Channel");
	grMeanMinOrg->GetXaxis()->SetTitleOffset(1.);
	grMeanMinOrg->GetXaxis()->SetTitleSize(0.05);
	grMeanMinOrg->GetXaxis()->SetLabelSize(0.05);	
	grMeanMinOrg->GetYaxis()->SetTitle("Mean_Min_Org [V]");
	grMeanMinOrg->GetYaxis()->SetTitleOffset(1.);
	grMeanMinOrg->GetYaxis()->SetTitleSize(0.05);
	grMeanMinOrg->GetYaxis()->SetLabelSize(0.05);	
	grMeanMinOrg->Draw("AP*");
	grMeanMinOrg->SetMarkerStyle(8);
	ca->cd(2);
	TGraph *grsgmMinOrg= new TGraph(nPresCh,chn,sgmMinOrg);
	grsgmMinOrg->SetTitle("Sigma MinValue for Original Wave");
	grsgmMinOrg->GetXaxis()->SetTitle("Channel");
	grsgmMinOrg->GetXaxis()->SetTitleSize(0.05);
	grsgmMinOrg->GetXaxis()->SetTitleOffset(1.);
	grsgmMinOrg->GetXaxis()->SetLabelSize(0.05);	
	grsgmMinOrg->GetYaxis()->SetTitle("Sigma_Min_Org [V]");
	grsgmMinOrg->GetYaxis()->SetTitle("Channel");
	grsgmMinOrg->GetYaxis()->SetTitleSize(0.05);
	grsgmMinOrg->GetYaxis()->SetTitleOffset(1.);
	grsgmMinOrg->GetYaxis()->SetLabelSize(0.05);	
	grsgmMinOrg->Draw("AP*");
	grsgmMinOrg->SetMarkerStyle(8);
	c3->cd(2);
	TGraph *grRatioFlt= new TGraph(nPresCh,chn,ratioFlt);
	grRatioFlt->SetTitle("Mean/Sigma for Filter Wave");
	grRatioFlt->GetXaxis()->SetTitle("Channel");
	grRatioFlt->GetXaxis()->SetTitleSize(0.05);
	grRatioFlt->GetYaxis()->SetTitle("Mean/Sigma");
	grRatioFlt->Draw("AP*");
	grRatioFlt->SetMarkerStyle(8);
	grRatioFlt->GetYaxis()->SetTitleOffset(0.7);
	grRatioFlt->GetYaxis()->SetTitleSize(0.05);
	c3->cd(1);
	TGraph *grRatioOrg= new TGraph(nPresCh,chn,ratioOrg);
	grRatioOrg->SetTitle("Mean/Sigma for Original Wave");
	grRatioOrg->GetXaxis()->SetTitle("Channel");
	grRatioOrg->GetXaxis()->SetTitleSize(0.05);
	grRatioOrg->GetYaxis()->SetTitle("Mean/Sigma");
	grRatioOrg->GetYaxis()->SetTitleOffset(0.7);
	grRatioOrg->GetYaxis()->SetTitleSize(0.05);
	grRatioOrg->Draw("AP*");
	grRatioOrg->SetMarkerStyle(8);
	c2->Print("/home/federica/eclipse-workspace/test_beam/file_root/image_minvalue/MeanSigmaFlt.png", "png");
	ca->Print("/home/federica/eclipse-workspace/test_beam/file_root/image_minvalue/MeanSigmaOrg.png", "png");

	c3->Print("/home/federica/eclipse-workspace/test_beam/file_root/image_minvalue/RatioMeanvsSigma.png", "png");
	//grafici rapporto tra media(sigma) filtrata e media(sigma) non filtrata,	
	c4->cd(1);
	TGraph *grRatioMean= new TGraph(nPresCh,chn,ratioMean);
	grRatioMean->SetTitle("Mean Ratio FltWave/OriginalWave");
	grRatioMean->GetXaxis()->SetTitle("Channel");
	grRatioMean->GetXaxis()->SetTitleSize(0.05);
	grRatioMean->GetXaxis()->SetLabelSize(0.05);
	grRatioMean->GetYaxis()->SetTitle("Mean Ratio");
	grRatioMean->GetYaxis()->SetTitleOffset(0.7);
	grRatioMean->GetYaxis()->SetTitleSize(0.05);
	grRatioMean->GetYaxis()->SetLabelSize(0.05);
	grRatioMean->Draw("AP*");
	grRatioMean->SetMarkerStyle(8);
	c4->cd(2);
	TGraph *grRatioSigma= new TGraph(nPresCh,chn,ratioSigma);
	grRatioSigma->SetTitle("Sigma Ratio FltWave/OriginalWave");
	grRatioSigma->GetXaxis()->SetTitle("Channel");
	grRatioSigma->GetXaxis()->SetTitleSize(0.05);
	grRatioSigma->GetXaxis()->SetLabelSize(0.05);
	grRatioSigma->GetYaxis()->SetTitle("Sigma Ratio");
	grRatioSigma->GetYaxis()->SetTitleOffset(0.7);
	grRatioSigma->GetYaxis()->SetTitleSize(0.05);
	grRatioSigma->GetYaxis()->SetLabelSize(0.05);
	grRatioSigma->Draw("AP*");
	grRatioSigma->SetMarkerStyle(8);
	c4->Print("/home/federica/eclipse-workspace/test_beam/file_root/image_minvalue/Ratio_FltWavevsOrgWave_signal.png", "png");

	TCanvas *cerr=new TCanvas("cerr","cerr");
	cerr->Divide(1,2);
	cerr->cd(1);
	TGraph *g_errNoFlt=new TGraph(ch,chn,errNoFlt);
	g_errNoFlt->SetTitle("Min/rms for original wave");
	g_errNoFlt->GetXaxis()->SetTitle("Entries");
	g_errNoFlt->GetYaxis()->SetTitle("|mean/rms|");
	g_errNoFlt->GetYaxis()->SetRangeUser(2,6);
	g_errNoFlt->Draw("AP*");
	g_errNoFlt->SetMarkerStyle(8);
	cerr->cd(2);
	TGraph *g_errFlt=new TGraph(ch,chn,err);
	g_errFlt->SetTitle("Min/rms for filter wave");
	g_errFlt->GetXaxis()->SetTitle("Entries");
	g_errFlt->GetYaxis()->SetTitle("|mean/rms|");
	g_errFlt->Draw("AP*");
	g_errFlt->SetMarkerStyle(8);

	cerr->Print("/home/federica/eclipse-workspace/test_beam/file_root/image_minvalue/Err.png", "png");

}
