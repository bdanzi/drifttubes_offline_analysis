/*
 * read_resultFT.C
 *
 *  Created on: 8 gen 2020
 *      Author: federica
 */
#include "stdio.h"
#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TText.h"
#include "TCanvas.h"
#include "TLegend.h"
#include <vector>
#include "TMath.h"
#include "TString.h"
#include <dirent.h>
#include <Riostream.h>
#include <map>
#include <vector>
#include "funcUtils.h"

using namespace std;

void read_resultFT(TString file="",TString outname="", TString  fOutName="",int m=0){
	//num bin su cui fai smooth=m*2+1
	ofstream myfile;
	  if (!fOutName.IsNull()) {
	    myfile.open (fOutName.Data());
	  }


	string line;

	float tmpGain,tmpEG,tmpfreq2,tmpEfreq,tmpGainTot,tmpEGtot,tmpPhase_diff,tmpEphase_diff,tmpPhase_diffTot,tmpEphase_diffTot;
	vector<float> gain;
	vector<float> gainTot;
	vector<float> freq;
	vector<float> Efreq;
	vector<float> Egain;
	vector<float> EgainTot;
	vector<float> phase_diff;
	vector<float> Ephase_diff;
	vector<float> phase_diffTot;
	vector<float> Ephase_diffTot;
	gROOT->SetStyle("Modern");
	gStyle->SetStatX(1);
	TCanvas *c_gain=new TCanvas("c0","c_gain",1500,1500);
	TCanvas *c_gainTot=new TCanvas("c1","c_gaintot",1500,1500);
	TCanvas *c_phase=new TCanvas("c2","c_phase",1500,1500);
	TCanvas *c_phaseTot=new TCanvas("c3","c_phaseTot",1500,1500);
	vector < pair <float,pair<float,float> > >vecF;
	pair<float,float> pairF;//seconda coppia
	pair <float,pair<float,float> > pairA;
	vector<float>gainSG;
	vector<float>phase_diffSG;
	ifstream tmpf;
	tmpf.open(file.Data());
	int size=0;
	while(getline(tmpf,line)) {
		stringstream tmpL(line);
		tmpL>>tmpfreq2>>tmpEfreq>>tmpGain>>tmpEG>>tmpGainTot>>tmpEGtot>>tmpPhase_diff>>tmpEphase_diff>>tmpPhase_diffTot>>tmpEphase_diffTot;
//		cout<<tmpfreq2<<"  "<<tmpEfreq<<" "<<tmpGain<<endl;
		freq.push_back(tmpfreq2*1e+6);
		size=freq.size()-2*m;
		Efreq.push_back(tmpEfreq);
		gain.push_back(tmpGain);
		gainSG=smooth(gain,m);
		Egain.push_back(tmpEG);
		gainTot.push_back(tmpGainTot);
		phase_diff.push_back(tmpPhase_diff);
		phase_diffSG=smooth(phase_diff,m);
		Ephase_diff.push_back(tmpEphase_diff);
		phase_diffTot.push_back(tmpPhase_diffTot);
		Ephase_diffTot.push_back(tmpEphase_diffTot);

		for (int i=0;i<gainSG.size();++i){
			pairF.first=gainSG[i];
//			cout<<"pairF.first "<<pairF.first<<" \t"<<"gain "<<gainSG[i]<<endl;
		}

		for (int i=0;i<phase_diffSG.size();++i){
			pairF.second=phase_diffSG[i]*TMath::RadToDeg();
		}

//		pairF.first=tmpGain;
//		pairF.second=tmpPhase_diff*TMath::RadToDeg();
		for(int i=0;i<size;++i){
			pairA.first=freq[i];
		}
		pairA.second=pairF;
		vecF.push_back(pairA);

		//pairA.first=tmpfreq2*1e+6;
	}
	tmpf.close();
//	for (int i=0;i<gainSG.size();++i){
//		cout<< gainSG[i]<< "\t";
//	}
//	cout<<size<<"\t"<<gainSG.size()<<"\t"<<gain.size()<<"\t"<<vecF.size()<<endl;

	sort(vecF.begin(), vecF.end());

//	for (vector<pair<float, pair<float, float > > >::iterator it = vecF.begin(); it < vecF.end();++it){
//		cout << it->first <<"\t"<< it->second.first <<"\t"<<  it->second.second<<endl;
//	}
	for (vector<pair<float, pair<float, float > > >::iterator it = vecF.begin(); it < vecF.end();++it){
		for (vector<pair<float, pair<float, float > > >::iterator it2 = it+1; it2 < vecF.end();){
			if(it->first==it2->first){
				vecF.erase(it2);
			}else{++it2;}
		}
	}



	for (vector<pair<float, pair<float, float > > >::iterator it = vecF.begin(); it < vecF.end();++it){
		if(myfile.is_open()){
//			if(it->first<2.e+8){
			myfile<<"\t"<<it->first<<"\t"<<it->second.first<<"\t"<<it->second.second<<endl;
//			}
		}
	}


	/***gain uscita catena-uscita buffer***/



	TGraphErrors *gr_gain= new TGraphErrors(freq.size(), &freq[0],&gain[0],&Efreq[0],&Egain[0]);
//	TGraph *gr_gain= new TGraph(freq.size(),&freq[0],&gain[0]);
	c_gain->cd();
	c_gain->SetLogx();
//  c_gain->SetLogx();
	gr_gain->SetMarkerStyle(7);
//	gr_gain->SetTitle(Form("Transfer function for frequencies %.1f MHz-%d MHz",frstRg[runId],lstRg[runId]));
	gr_gain->SetTitle("Gain for frequencies 0.1MHz-1500 MHz");
	gr_gain->GetYaxis()->SetRangeUser(-1,40);
//	gr_gain->GetXaxis()->SetRangeUser(frstRg[runId],lstRg[runId]);
	gr_gain->GetYaxis()->SetTitle("Gain [dB]");
	gr_gain->GetXaxis()->SetTitleOffset(1.2);
	gr_gain->GetXaxis()->SetTitle("Freq [MHz]");
	gr_gain->Draw("AP");

	/***gain uscita catena-uscita generatore***/

//	TGraphErrors *gr_gainTot= new TGraphErrors(freq.size(), &freq[0],&gainTot[0],&Efreq[0],&EgainTot[0]);
	TGraph *gr_gainTot= new TGraph(freq.size(),&freq[0],&gainTot[0]);
	c_gainTot->cd();
//  c_gainTot->SetLogx();
	gr_gainTot->SetMarkerStyle(7);
	gr_gainTot->SetTitle("Gain for frequencies 0.1MHz-1500 MHz of all chain");
	gr_gainTot->SetLineColor(kBlue);
	gr_gainTot->SetMarkerColor(kBlue);
	gr_gainTot->GetYaxis()->SetRangeUser(-1,100);
	gr_gainTot->GetYaxis()->SetTitle("Gain");
	gr_gainTot->GetXaxis()->SetTitle("Freq [MHz]");
	gr_gainTot->Draw("AP");


	/***fase uscita catena-uscita buffer***/
	TGraphErrors *gr_phase=new TGraphErrors(freq.size(), &freq[0],&phase_diff[0],&Efreq[0],&Ephase_diff[0]);
//	TGraph *gr_phase= new TGraph(freq.size(),&freq[0],&phase_diff[0]);
	c_phase->cd();
	gr_phase->SetMarkerStyle(7);
	gr_phase->SetTitle("Phase difference");
	//gr_phase->SetMarkerColor(kBlue);
	gr_phase->GetXaxis()->SetTitle("Freq [MHz]");
	gr_phase->GetYaxis()->SetTitle("Phase difference");
	gr_phase->Draw("AP");

	/***fase uscita catena-uscita generatore***/
	TGraphErrors *gr_phaseTot=new TGraphErrors(freq.size(), &freq[0],&phase_diffTot[0],&Efreq[0],&Ephase_diffTot[0]);
//	TGraph *gr_phase= new TGraph(freq.size(),&freq[0],&phase_diff[0]);
	c_phaseTot->cd();
	gr_phaseTot->SetMarkerStyle(7);
	gr_phaseTot->SetTitle("Phase difference for all chain");
	gr_phaseTot->GetXaxis()->SetTitle("Freq [MHz]");
	gr_phaseTot->GetYaxis()->SetTitle("Phase difference");
	gr_phaseTot->SetLineColor(kBlue);
	gr_phaseTot->Draw("AP");



	if (myfile.is_open()) { myfile.close(); }
    /**********************/
    //  SAVING
    /*********************/
    //rootpla
    TFile *saving_rootpla = new TFile (outname.Data(),"recreate");
    c_gain->Write("gr_gain");
    c_phase->Write("gr_phase");
    c_gainTot->Write("gr_gain_tot");
    c_phaseTot->Write("gr_phase_tot");
    saving_rootpla->Close();


}


