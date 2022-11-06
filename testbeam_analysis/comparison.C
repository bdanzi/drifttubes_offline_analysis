/*
 * comparison.C
 *
 *  Created on: 4 nov 2019
 *      Author: federica
 */
#include <TStyle.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TVector.h>
#include "TF1.h"
#include "TFitResult.h"
#include "TMinuit.h"
#include "TGraphErrors.h"
#include <fstream>
#include <sstream>
#include <string>
#include "TSpectrum.h"
#include "TAxis.h"
#include "TMath.h"
#include "TH1F.h"
#include "TFitResult.h"
#include <iostream>
#include <stdio.h>
#include <dirent.h>
#include <map>
#include <iterator>
#include <vector>

using namespace std;
void comparison(TString fldName="."){

    struct dirent *entry = nullptr;
    DIR *dp = nullptr;

    string line;

	std::vector<std::string> hmix;
	std::vector<float> minHV;
	std::vector<float> maxHV;
	std::vector<float> Chs[20];
	std::vector<float> Slope[20];
	std::vector<float> ErrSlope[20];
	std::vector<float> Const[20];
	std::vector<float> ErrConst[20];
	Float_t mean_slope[3]={0};
	Float_t meanE_slope[3]={0};
	Float_t mean_const[3]={0};
	Float_t meanE_const[3]={0};
	EColor colors[4]={kBlue,kRed,kGreen,kPink};


	int tmpCh;
	float tmpSlope, tmpConst,tmpErrSlope,tmpErrConst;
	int nGas=0;

    dp = opendir(fldName.Data());
    if (dp != nullptr) {
        while ((entry = readdir(dp)) && nGas<20) {
            TString tmpFN=entry->d_name;
			cout<<"----------------- reading file: "<<tmpFN<<endl;
			if (tmpFN.Contains("txt")) {
				TString tmpGain(tmpFN(0,5));
				hmix.push_back(tmpGain.Data());
				minHV.push_back(TString(tmpFN(6,4)).Atof());
				maxHV.push_back(TString(tmpFN(11,4)).Atof());
				cout<<tmpFN<<" -> "<<tmpGain<<" minHV "<<minHV.back()<<" maxHV "<<maxHV.back()<<endl;
				ifstream tmpf;//lettura file
				tmpFN=fldName+"/"+tmpFN;
				tmpf.open(tmpFN.Data());
				tmpCh=-1;
				tmpSlope=tmpConst=tmpErrSlope=tmpErrConst=0.0;
				while(getline(tmpf,line)/*tmpf.good()*/) {
					stringstream tmpL(line);
					tmpL>>tmpCh>>tmpSlope>>tmpErrSlope>>tmpConst>>tmpErrConst;
					cout<<tmpCh<<"  "<<tmpSlope<<endl;
					Chs[nGas].push_back(tmpCh);
					//cout<<" canali "<< Chs.back()<<endl;
					Slope[nGas].push_back(tmpSlope);
					//Slope[nGas]=tmpSlope;
					cout<<"Slope["<<nGas<<"]"<<Slope[nGas].back()<<endl;
					ErrSlope[nGas].push_back(tmpErrSlope);
					Const[nGas].push_back(tmpConst);
					ErrConst[nGas].push_back(tmpErrConst);
				}
				tmpf.close();
				++nGas;
				cout<<" ngas "<<nGas<<endl;
			}
		}
    }
    closedir(dp);


	gStyle->SetOptFit(1);
	gStyle->SetStatX(0.6);
	gStyle->SetStatY(0.9);
	gStyle->SetStatW(0.2);
	gStyle->SetStatH(0.2);

	TF1 *pol =new TF1("pol","pol0");
	pol->SetNpx(1000);


//slope
    TCanvas *cvSlops = new TCanvas("cvSlops","cvSlops",1500,1500);
    cvSlops->Divide(nGas/2+nGas%2,2);
    TGraphErrors **grSlope=new TGraphErrors* [20];
    for(int i=0;i<nGas;++i){
    	grSlope[i] = new TGraphErrors(Chs[i].size(),Chs[i].data(),Slope[i].data(),0,ErrSlope[i].data());
    	grSlope[i]->SetTitle(Form("Slope of gas %s",hmix[i].c_str()));
    	grSlope[i]->SetMarkerStyle(8);
    	cvSlops->cd(i+1);
    	grSlope[i]->Draw("AP");
    	grSlope[i]->Fit("pol");
    	grSlope[i]->GetXaxis()->SetRangeUser(-1,12);
    	grSlope[i]->GetXaxis()->SetTitle("Channel");
    	grSlope[i]->GetYaxis()->SetTitle("Slope");
    	grSlope[i]->GetYaxis()->SetTitleOffset(1.65);

		mean_slope[i]=pol->GetParameter(0);
		meanE_slope[i]=pol->GetParError(0);
		cout << "i"<<i<<" slope "<< mean_slope [i]<< endl;
    }
   cvSlops->Print("/lustrehome/bdanzi/offline_analysis/testbeam_analysis/Plots/histosTB_run_4.root/Slope.png", "png");
//   cvSlops->Print("/home/federica/eclipse-workspace/test_beam/Plot/Confronto_miscele/Slope.png", "png");
   //    cvSlops->Print("/home/federica/eclipse-workspace/test_beam/Plot/confronto_miscele_control/norm_loss/Slope.png", "png");


    //constant
    TCanvas *cvConst = new TCanvas("cvConst","cvConst",1500,1500);
    cvConst->Divide(nGas/2+nGas%2,2);
    TGraphErrors **grConst=new TGraphErrors* [20];
    for(int i=0;i<nGas;++i){
    	grConst[i] = new TGraphErrors(Chs[i].size(),Chs[i].data(),Const[i].data(),0,ErrConst[i].data());
    	grConst[i]->SetTitle(Form("Constant of gas %s",hmix[i].c_str()));
    	grConst[i]->SetMarkerStyle(8);
    	cvConst->cd(i+1);
    	grConst[i]->Draw("AP");
    	grConst[i]->Fit("pol");
    	grConst[i]->GetXaxis()->SetRangeUser(-1,12);
    	grConst[i]->GetXaxis()->SetTitle("Channel");
    	grConst[i]->GetYaxis()->SetTitle("Constant");

		mean_const[i]=pol->GetParameter(0);
		meanE_const[i]=pol->GetParError(0);

		cout<<" i "<< i <<" costante "<<mean_const[i] <<endl;
    }
       cvConst->Print("/lustrehome/bdanzi/offline_analysis/testbeam_analysis/Plots/histosTB_run_4.root/Constant.png", "png");
//	   cvConst->Print("/home/federica/eclipse-workspace/test_beam/Plot/Confronto_miscele/Constant.png", "png");
//       cvConst->Print("/home/federica/eclipse-workspace/test_beam/Plot/confronto_miscele_control/norm_loss/Constant.png", "png");

    float norm[3]={18.08,24.79,31.75};
    TCanvas *cExp=new TCanvas();
    cExp->Divide(nGas/2+nGas%2,2);
    TF1 **exp= new TF1* [20];
    for(int i=0;i<nGas;++i){
    exp[i]=new TF1(Form("exp%i",i),"expo",minHV[i],maxHV[i]);
    exp[i]->SetTitle(Form("Gain of gas %s",hmix[i].c_str()));
    exp[i]->SetParameter(0,mean_const[i]);
    exp[i]->SetParError(0, meanE_const[i]);
    exp[i]->SetParameter(1,mean_slope[i]);
    exp[i]->SetParError(1, meanE_slope[i]);
    exp[i]->SetLineColor(colors[i]);
    exp[i]->GetXaxis()->SetTitle("HV [V]");
    cExp->cd(i+1);
    exp[i]->Draw();
    }

const double ymax=6;//3.2;//13.0; //3.2; //3.2 e 1.8 sono i valori per gli histo normalizzati
const double ymin=0.0;//2.0;//1.8;


    TCanvas *cGain=new TCanvas("cGain","Canvas for Gains",1500,1500);
    TH1F *hGain=new TH1F("hGain","Gain vs HV for different gas mixtures",100,1490,1840);
    hGain->SetStats(0);
    hGain->SetMaximum(ymax);
    hGain->SetMinimum(ymin);
    cGain->cd();
    hGain->Draw();
    hGain->GetXaxis()->SetTitle("HV [V]");
    hGain->GetYaxis()->SetTitle("Gain (10^5)");
    exp[0]->Draw("SAME");
    exp[1]->Draw("SAME");
    exp[2]->Draw("SAME");

    //legend
    TLegend *legend = new TLegend(0.1,0.7,0.48,0.9);
    legend->SetHeader("Gas mixtures","C");
    for(int i=0;i<nGas;i++){
    legend->AddEntry(exp[i],hmix[i].c_str(),"l");
    legend->Draw();
    }
    cGain->Print("/lustrehome/bdanzi/offline_analysis/testbeam_analysis/Plots/histosTB_run_4.root/GainTot.png", "png");
  //  cGain->Print("/home/federica/eclipse-workspace/test_beam/Plot/Confronto_miscele/GainTot.png", "png");
 //   cGain->Print("/home/federica/eclipse-workspace/test_beam/Plot/confronto_miscele_control/norm_loss/GainTot.png", "png");




    /*for (Slope_it= Slope.begin(); Slope_it!=Slope.end(); ++Slope_it){
       cout << Slope_it->first << " => " << Slope_it->second[0] << '\n';
    	}*/
	}






