#include <stdio.h>
#include <dirent.h>
#include <Riostream.h>
#include <map>
#include <vector>
#include <iterator>
#include "TGraphErrors.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TMath.h"

void PlotChargeInteg(TString fldName=".",int cnvId=0, TString fOutName="") {///home/federica/eclipse-workspace/test_beam/Plot/curves.txt
    struct dirent *entry = nullptr;
    DIR *dp = nullptr;

    //scrittura di dati su file
	 ofstream myfile;
	  if (!fOutName.IsNull()) {
	    myfile.open (fOutName.Data());
	  }

	float e= 1.602e-19;
	float R=330;//impedenza camera
	float gEl=9.74; //guadagno elettronica
	float ergEl=0.08; //errore sul guadagno dell'elettronica
	float frequ_camp=0.8e-9;//larghezza temporale bin
	float loss=0.50;
	float convFact=(frequ_camp*1e-5/(e*R*gEl*loss))*2.0;
	float a=convFact/18.08;//90-10
	float b=convFact/24.79;//85-15
	float c=convFact/31.75;//80-20
	float d=convFact/16;//95-05?
	float convFactC[5]={1,a,b,c,d};
	cout<< " fattore di conversione "<< convFactC[1]<< endl;
	cout<< " fattore di conversione senza norm "<< convFact << endl;

	std::vector<float> hVs;
	std::vector<int> Chs;
	std::map<int, std::vector<float> > MPv;
	std::map<int, std::vector<float> > errMPv;
	std::map<int, std::vector<float> > Sgm;
	std::map<int, std::vector<float> > errSgm;
	
	int tmpCh;
	float tmpMPv, tmpErrMPv, tmpSgm, tmpErrSgm;
	
    dp = opendir(fldName.Data());
    if (dp != nullptr) {
        while ((entry = readdir(dp))) {
            TString tmpFN=entry->d_name;
			if (tmpFN.Contains("txt")) {
				TString tmpHVs(tmpFN(0,tmpFN.Sizeof()-5));
				cout<<tmpFN<<" -> "<<tmpHVs<<endl;
				hVs.push_back(tmpHVs.Atof());
				ifstream tmpf;
				tmpFN=fldName+"/"+tmpFN;
				tmpf.open(tmpFN.Data());
				tmpCh=-1;
				tmpMPv=tmpErrMPv=tmpSgm=tmpErrSgm=0.0;
				while(!tmpf.eof()) {
					tmpf>>tmpCh>>tmpMPv>>tmpErrMPv>>tmpSgm>>tmpErrSgm;
					//cout<<tmpCh<<" tmpMPv "<<tmpMPv<<endl;
					//Chs.push_back(tmpCh);
					//MPv[tmpCh].push_back(tmpMPv*convFact);//convfact per gli isto normalizzati
					MPv[tmpCh].push_back(tmpMPv*convFactC[cnvId]);
					//errMPv[tmpCh].push_back(convFact*TMath::Sqrt(tmpErrMPv*tmpErrMPv+((tmpMPv*tmpMPv*ergEl*ergEl)/(gEl*gEl))));
					errMPv[tmpCh].push_back((TMath::Sqrt(tmpErrMPv*tmpErrMPv+((tmpMPv*tmpMPv*ergEl*ergEl)/(gEl*gEl))))*convFactC[cnvId]);

					Sgm[tmpCh].push_back(tmpSgm);
					errSgm[tmpCh].push_back(tmpErrSgm);
					tmpCh++;
				}
				tmpf.close();
			}
		}
    }
    closedir(dp);
	
	const int n=9;
	Float_t ch[9]={0,1,2,3,4,5,9,10,11};
	Float_t slope[9]={0,0,0,0,0,0,0,0,0};
	Float_t constant[9]={0,0,0,0,0,0,0,0,0};
	Float_t constErr[9]={0};
	Float_t err[9]={0,0,0,0,0,0,0,0,0};
	gStyle->SetOptFit(1);
	gStyle->SetStatX(0.6);
	gStyle->SetStatY(0.9);
	gStyle->SetStatW(0.2);                
	gStyle->SetStatH(0.2);   

	TF1 *exp =new TF1("exp","expo");
	exp->SetNpx(1000);

	TCanvas *cv = new TCanvas("cv","cv",1500,1500);
	cv->Divide(3,3);
	TGraphErrors **grMPs = new TGraphErrors* [MPv.size()];
	std::map<int, std::vector<float> >::iterator MPv_it = MPv.begin();
	std::map<int, std::vector<float> >::iterator errMPv_it = errMPv.begin();
	for (int ich=0; ich<MPv.size()-1; ++ich) {
		grMPs[ich] = new TGraphErrors(hVs.size(),&hVs[0],&MPv_it->second[0],0,&errMPv_it->second[0]);
		grMPs[ich]->SetMarkerStyle(8);
		grMPs[ich]->SetTitle(Form("Gain ch%i",MPv_it->first));
		grMPs[ich]->GetXaxis()->SetTitle("HV [V]");
		grMPs[ich]->GetYaxis()->SetTitle("Gain (10^5)");
		grMPs[ich]->GetYaxis()->SetTitleSize(0.06);
		grMPs[ich]->GetYaxis()->SetTitleOffset(0.8);
		grMPs[ich]->GetXaxis()->SetTitleSize(0.05);		
		grMPs[ich]->GetXaxis()->SetLabelSize(0.05);		
		grMPs[ich]->GetYaxis()->SetLabelSize(0.05);		
		grMPs[ich]->Fit(exp);

		slope[ich]=exp->GetParameter(1);
		err[ich]=exp->GetParError(1);
		constant[ich]=exp->GetParameter(0);
		constErr[ich]=exp->GetParError(0);
		//cout<< " canale "<< ich<< " costante " << constant[ich] << endl;
		if(myfile.is_open()){
			myfile << ch[ich]<< "\t" << slope[ich] << "\t"<<err[ich]<<"\t" <<constant[ich] << "\t"<<constErr[ich]<<endl;
		}
		cv->cd(ich+1);
		grMPs[ich]->Draw("AP");
		MPv_it++;
		errMPv_it++;
	}
	cv->Print("/lustrehome/bdanzi/offline_analysis/testbeam_analysis/Plots/histosTB_run_4.root/gain_chargeinteg.png", "png");

	//grafico per la pendenza
	TCanvas *cSlope=new TCanvas();
	cSlope->cd();
	TGraphErrors *grSlope= new TGraphErrors(n,ch,slope,0,err);
	grSlope->SetTitle("Slope of exponential fit");
	grSlope->SetMarkerStyle(8);
	grSlope->GetXaxis()->SetTitle("Channel");
	grSlope->GetYaxis()->SetTitle("Slope of exponential fit");
	grSlope->GetYaxis()->SetRangeUser(0.14,0.008);
	grSlope->GetXaxis()->SetRangeUser(-1,12);
	grSlope->GetYaxis()->SetTitleOffset(1.25);
	grSlope->Draw("AP");
	cSlope->Print("/lustrehome/bdanzi/offline_analysis/testbeam_analysis/Plots/histosTB_run_4.root/slope_chargeinteg.png", "png");

	//grafico per la costante
	TCanvas *cConst=new TCanvas();
	cConst->cd();
	TGraphErrors *grConst=new TGraphErrors(n,ch,constant,0,constErr);
	grConst->SetTitle("Constant of exponential fit");
	grConst->SetMarkerStyle(8);
	grConst->GetXaxis()->SetTitle("Channel");
	grConst->GetYaxis()->SetTitle("Constant of exponential fit");
	grConst->GetYaxis()->SetTitleOffset(1.25);
	grConst->GetXaxis()->SetRangeUser(-1,12);
	grConst->Draw("AP");
	cConst->Print("/lustrehome/bdanzi/offline_analysis/testbeam_analysis/Plots/histosTB_run_4.root/Constant_chargeinteg.png", "png");


	if (myfile.is_open()) { myfile.close(); }
}
