#include "Riostream.h"
#include <TCanvas.h>
#include <TGraph.h>
//#include<TVector.h>
//#include "TF1.h"
//#include "TDirectory.h"

#include <map>
#include <vector>
#include <fstream>

void PlotNSelEv(int maxNSgn=9) {
	std::map<int, std::vector<float> > dataCont;
	for (int i=1; i<=maxNSgn; ++i) {
		ifstream iSelFl(Form("nSelEv-5000-sgm-%d.txt",i));
		int tmpCh, tmpNEv;
//		while (!iSelFl.eof()) {
		while (iSelFl>>tmpCh>>tmpNEv) {
//			cout<<tmpCh<<"\t"<<tmpNEv<<endl;
			dataCont[tmpCh].push_back(((float)tmpNEv)/5000.0);
		}
		iSelFl.close();
	}

	std::vector<float> nsgm;
	for (int isgm=1; isgm<=maxNSgn; ++isgm) {
		nsgm.push_back(isgm);
	}
	std::vector<TGraph *> grps;
	std::vector<TCanvas *> cvs;
	for (std::map<int, std::vector<float> >::iterator it = dataCont.begin(); it != dataCont.end(); ++it) {
		//cout<<"Ch "<<it->first;
		for (int ipt=0; ipt<it->second.size(); ++ipt) {
			//cout<<"\t"<<it->second[ipt];
		}
//		cout<<endl;
		cvs.push_back(new TCanvas(Form("CV-ch%d",it->first)));
		cvs.back()->cd();
		grps.push_back(new TGraph(maxNSgn,nsgm.data(),it->second.data()));
		grps.back()->SetMarkerStyle(8);
		grps.back()->SetTitle(Form("Ev selection for ch %d",it->first));
		grps.back()->Draw("AP");
	}

}
