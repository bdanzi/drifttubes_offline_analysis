/*
 * read_FT.C
 *
 *  Created on: 5 dic 2019
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
using namespace std;

struct data0{
	float_t ampl0=0.0;
	float_t Eampl0=0.0;
	float_t phase0=0.0;
	float_t Ephase0=0.0;
	float_t freq0=0.0;
	float_t Efreq0=0.0;
};

struct data1{
	float_t ampl1=0.0;
	float_t Eampl1=0.0;
	float_t phase1=0.0;
	float_t Ephase1=0.0;
	float_t freq1=0.0;
	float_t Efreq1=0.0;
};

struct data2{
	float_t ampl2=0.0;
	float_t Eampl2=0.0;
	float_t phase2=0.0;
	float_t Ephase2=0.0;
	float_t freq2=0.0;
	float_t Efreq2=0.0;
};

void read_FT(TString fname="",TString outname=""/* int runId=0*/,TString fOutName=""){

	  ofstream myfile;
	  if (!fOutName.IsNull()) {
	    myfile.open (fOutName.Data());
	  }



	TObjArray *tx = fname.Tokenize(";");
	tx->Print();

	data0 spettro0;
	data1 spettro1;
	data2 spettro2;

	int ChEntries =0;
	float tmpGain=0.;
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

	float ampl_in,ampl_out,tmpEG,Eampl_out,Eampl_in,tmpEfreq,ampl_inTot,Eampl_inTot,tmpGainTot,tmpEGtot,tmpPhase_diff,tmpEphase_diff,tmpPhase_diffTot,tmpEphase_diffTot;
	TCanvas *c_gain=new TCanvas("c0","c_gain",1500,1500);
	TCanvas *c_gainTot=new TCanvas("c1","c_gaintot",1500,1500);
	TCanvas *c_phase=new TCanvas("c2","c_phase",1500,1500);
	TCanvas *c_phaseTot=new TCanvas("c3","c_phaseTot",1500,1500);

	float frstRg[8]={1,0.1,1,10,300,500,1,1000};
	int lstRg[8]={1,1,10,300,500,1000,1,1200};

	TFile *f_Spectrum = new TFile(((TObjString *)(tx->At(0)))->String().Data(),"read");
	if (!f_Spectrum->IsOpen()) {
		cout<<" Error opening Reference root file! "<<endl;
		return 0;
	}
	else {cout<< " ok, fatto "<< endl;}

/***tree***/
//	TTree *t = (TTree*)f_Spectrum->Get("tree");
//	TBranch *b_branch0=t->GetBranch("branch0");
//	b_branch0->SetAddress(&spettro0.ampl0);
//	TBranch *b_branch1=t->GetBranch("branch1");
//	b_branch1->SetAddress(&spettro1.ampl1);
//	TBranch *b_branch2=t->GetBranch("branch2");
//	b_branch2->SetAddress(&spettro2.ampl2);

//	NEntries=t->GetEntries();
// cout<<" entries "<< NEntries<<endl;
/***tree***/

/***chain***/
		TChain *t=new TChain("tree");

		for (Int_t i = 0; i < tx->GetEntries(); i++){
			cout << ((TObjString *)(tx->At(i)))->String() << endl;
			t->Add(((TObjString *)(tx->At(i)))->String());
		}

		TObjArray *fileElements=t->GetListOfFiles();
		ChEntries=t->GetEntries();
		t->SetBranchAddress("branch0",&spettro0);
		t->SetBranchAddress("branch1",&spettro1);
		t->SetBranchAddress("branch2",&spettro2);

	for (int i=5;i<ChEntries-1;++i){
		t->GetEntry(i);
		//Long64_t ientry = t->LoadTree(i);
		//cout<<"conteggio della entries "<< ientry<<endl;
		//if (ientry < 0) break;
		ampl_in = spettro1.ampl1;//t->GetLeaf("ampl1")->GetValue();
//		cout << "ampl_in "<< ampl_in<<endl;
		ampl_out = spettro0.ampl0;//t->GetLeaf("ampl0")->GetValue();
//		cout << "ampl_out "<< ampl_out<<endl;
		ampl_inTot=spettro2.ampl2;
		freq.push_back(spettro2.freq2/*t->GetLeaf("freq2")->GetValue()*/);
		tmpEfreq = spettro2.Efreq2;//t->GetLeaf("Efreq2")->GetValue();
//		cout << "tmpEfreq "<< tmpEfreq<<endl;
		Eampl_out=spettro0.Eampl0;//t->GetLeaf("Eampl0")->GetValue();
		Eampl_in=spettro1.Eampl1;//t->GetLeaf("Eampl1")->GetValue();
		Eampl_inTot=spettro2.Eampl2;
//		cout << "Eampl_out "<< Eampl_out<<endl;
//		cout << "Eampl_in "<< Eampl_in<<endl;
		tmpGain=20*TMath::Log10(ampl_out/ampl_in);
//		tmpEG=(1/ampl_in)*TMath::Sqrt((ampl_out*ampl_out/ampl_in*ampl_in)*(Eampl_in*Eampl_in)+(Eampl_out*Eampl_out));
		tmpEG=(20/TMath::Log(10))*TMath::Sqrt((Eampl_out*Eampl_out/ampl_out*ampl_out)+(Eampl_in*Eampl_in/ampl_in*ampl_in));
		gain.push_back(tmpGain);
		Egain.push_back(tmpEG);

		Efreq.push_back(tmpEfreq);
//		if(spettro2.Efreq2>0.5){cout<< spettro2.freq2<< " "<< tmpGain<<endl;}

//		tmpGainTot=20*TMath::Log10(ampl_out/ampl_inTot);
		tmpGainTot=(spettro1.ampl1/spettro2.ampl2);
		gainTot.push_back(20*TMath::Log10(tmpGainTot));
		tmpEGtot=(1/ampl_inTot)*TMath::Sqrt((ampl_out*ampl_out/ampl_inTot*ampl_inTot)*(Eampl_inTot*Eampl_inTot)+(Eampl_out*Eampl_out));
		EgainTot.push_back(tmpEGtot);
		tmpPhase_diff=spettro0.phase0-spettro1.phase1;
		tmpEphase_diff=spettro0.Ephase0-spettro1.Ephase1;
		phase_diff.push_back(spettro0.phase0-spettro1.phase1);//(t->GetLeaf("phase0")->GetValue()-t->GetLeaf("phase1")->GetValue());
		Ephase_diff.push_back(spettro0.Ephase0-spettro1.Ephase1);//(t->GetLeaf("Ephase0")->GetValue()-t->GetLeaf("Ephase1")->GetValue());

		tmpPhase_diffTot=spettro0.phase0-spettro2.phase2;
		tmpEphase_diffTot=spettro0.Ephase0-spettro2.Ephase2;
		phase_diffTot.push_back(spettro0.phase0-spettro2.phase2);//(t->GetLeaf("phase0")->GetValue()-t->GetLeaf("phase1")->GetValue());
		Ephase_diffTot.push_back(spettro0.Ephase0-spettro2.Ephase2);//(t->GetLeaf("Ephase0")->GetValue()-t->GetLeaf("Ephase1")->GetValue());
		//cout << spettro2.freq2<<"\t"<<spettro2.Efreq2<<"\t"<<tmpGain<<"\t"<<tmpEG<<"\t"<<tmpGainTot<<"\t"<<tmpEGtot<<endl;
		if (myfile.is_open()) {
//			myfile<<spettro2.freq2<<"\t"<<spettro2.Efreq2<<"\t"<<tmpGain<<"\t"<<tmpEG<<"\t"<<tmpGainTot<<"\t"<<tmpEGtot<<"\t"<<tmpPhase_diff<<"\t"<<tmpEphase_diff<<"\t"<<tmpPhase_diffTot<<"\t"<<tmpEphase_diffTot<<endl;
			myfile<<"\t"<<spettro2.freq2<<"\t"<<tmpGain<<"\t"<<tmpPhase_diff<<endl;

		}
	}

	/***gain uscita catena-uscita buffer***/
	TGraphErrors *gr_gain= new TGraphErrors(freq.size(), &freq[0],&gain[0],&Efreq[0],&Egain[0]);
//	TGraph *gr_gain= new TGraph(freq.size(),&freq[0],&gain[0]);
	c_gain->cd();
//  c_gain->SetLogx();
	gr_gain->SetMarkerStyle(7);
//	gr_gain->SetTitle(Form("Transfer function for frequencies %.1f MHz-%d MHz",frstRg[runId],lstRg[runId]));
	gr_gain->SetTitle("Gain for frequencies 0.1MHz-1500 MHz");
	gr_gain->GetYaxis()->SetRangeUser(-1,100);
//	gr_gain->GetXaxis()->SetRangeUser(frstRg[runId],lstRg[runId]);
	gr_gain->GetYaxis()->SetTitle("Gain");
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

    if (myfile.is_open()) { myfile.close(); }
}


