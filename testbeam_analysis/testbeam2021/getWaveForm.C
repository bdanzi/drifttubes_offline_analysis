#include "/home/federica/eclipse-workspace/testBeam_CluCou/sw_test_beam/offline_analysis/testbeam_analysis/FindPeak-algo.C"
//TGraph* getWaveForm(int channel, int entry, const char* filename) {
	void getWaveForm(int channel, int entry, const char* filename) {

	auto gr = new TGraph();
	auto fl = TFile::Open(filename, "read");
	auto tree = static_cast<TTree*>(fl->Get("data"));
	int n = tree->GetEntries();
	WaveData* wd = new WaveData();
	tree->SetBranchAddress("x", &wd);
	tree->GetEntry(entry);
	Float_t *dataY=new Float_t [1024];
	std::vector<float> dataX;
	int NPeak;
	Int_t pkPos[100];
	Float_t pkHgt[100];
	int skipFstBin=1;
	for (auto point : wd->getX742Data()) {
		if (point.first == channel) {
			for (int i = 0; i < 1024; ++i) {
				gr->SetPoint(i, (double) i, -1*((double) point.second[i]-0.450));
				dataY[i]=-1*((double) point.second[i]-0.450);
				dataX.push_back((double) i);
			}
			break;
		}
	}

	NPeak = FindPeaks(1024,dataY,
			1.2e-3,pkPos,pkHgt);
	cout<<NPeak<<endl;
	gr->Draw("");
	for (int ipk=0; ipk<NPeak; ipk++){
		cout<<pkHgt[ipk]<<" "<<pkPos[ipk]<<endl;
		TMarker *tm = new  TMarker(dataX[pkPos[ipk]], pkHgt[ipk]/*+PeakBaseShift_1*/, 23);
		tm->SetMarkerSize(1.5);
		tm->SetMarkerColor(2);
		tm->Draw("same");
	}


	gr->SetMarkerStyle(7);
	fl->Close();
	delete wd;
	return gr;
}
