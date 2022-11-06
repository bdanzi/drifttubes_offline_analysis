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
#include "funcUtils.h"

using namespace std;

void Occupancy(TString file=""){

	gROOT->SetStyle("Modern");


	string line;
	ifstream tmpf;
	tmpf.open(file.Data());
	string s1,s2;
	int tmpEntry,tmpChannel;
	vector<int> channel;
	vector<int> entry;
	TH1F *h1=new TH1F("h1","Channel Occupancy",100,0,24);
	TH1F *h2=new TH1F("h2","Full events with no duplicates",1000,-1,2000);
	TH1F *h3=new TH1F("h3","Full events ",1000,-1,2000);
	TH1F *h4=new TH1F ("h4"," Multiplicity",20,0,20);
	int oldEntry;
	oldEntry=-1;
	int mult;
	mult=0;
	while(getline(tmpf,line)) {
		stringstream tmpL(line);
		tmpL>>s1>>tmpEntry>>s2>>tmpChannel;
		//		cout<<s1 <<" "<<tmpEntry<<" "<<s2<<" "<<tmpChannel<<endl;
		h1->Fill(tmpChannel);
		channel.push_back(tmpChannel);
		entry.push_back(tmpEntry);
		h3->Fill(tmpEntry);
		if(tmpEntry!=oldEntry){
			if(oldEntry!=-1)h4->Fill(mult);
			oldEntry=tmpEntry;
			mult=1;
		}
			else{mult++;}

	}
	h4->Fill(mult);

	vector<int> entry_nd=entry;
	int tmpentry_nd;
	sort(entry_nd.begin(),entry_nd.end());
	entry_nd.erase( unique( entry_nd.begin(), entry_nd.end() ), entry_nd.end() );
	cout<<" size of delete vector "<< entry_nd.size()<<" size of original vector "<<entry.size()<<endl;
	for(int i=0;i<entry_nd.size();++i){
		tmpentry_nd=entry_nd[i];
//		cout<<tmpentry_nd<<endl;
		h2->Fill(tmpentry_nd);
	}



	//______________________________________________Draw________________________________________________//
	TCanvas *c1=new TCanvas("c1","c1");
	c1->cd();
	h1->Draw();
	h1->GetXaxis()->SetTitle("Channel");
	h1->GetXaxis()->SetNdivisions(25);
	h1->GetYaxis()->SetTitle("Entries");

	TCanvas *c2=new TCanvas("c2","c2",1000,1000);
	TGraph *graph =new TGraph(entry.size(),entry.data(),channel.data());
	graph->SetMarkerStyle(8);
	graph->SetTitle("Event Distribution");
	graph->GetXaxis()->SetTitle("Event");
	graph->GetYaxis()->SetTitle("Channel");
	graph->Draw("AP");


	TCanvas *c3=new TCanvas("c3","c3",1000,1000);
	c3->cd();
	h2->Draw();
	h2->GetXaxis()->SetTitle("Full events");
	h2->GetYaxis()->SetTitle("Entries");


	TCanvas *c4=new TCanvas("c4","c4",1000,1000);
	c4->cd();
	h3->Draw();
	h3->GetXaxis()->SetTitle("Full events");
	h3->GetYaxis()->SetTitle("Entries");

	TCanvas *c5=new TCanvas();
	c5->cd();
	h4->Draw();
	h4->GetXaxis()->SetTitle("Multiplicity");
	h4->GetYaxis()->SetTitle("Entries");

	tmpf.close();
}
