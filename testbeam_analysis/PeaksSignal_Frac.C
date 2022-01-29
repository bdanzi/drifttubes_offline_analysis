//#include "SignalCont.h"
//#include "ConvSignalCont.h"
//#include "PeakCont.h"
//#include "PeakData.h"
//#include "CellCont.h"
#include "dataUtils.h"

void PeaksSignalFrac(TString fName="", Int_t swId=11, Float_t timeRes=0.5e-9, Float_t noise=0.9e-3, Int_t cEv=-1)
{

	TH1F *hNCl = new TH1F ("hNCl","n Cluster",100,0,100);
	TH1F *hNEl = new TH1F ("hNEl","n Primary Electrons",100,0,100);
	TH1F *hTotChr = new TH1F ("hTotChr","Total Charge",100,0,5e+7);
	TH1F *hChrPerCl = new TH1F ("hChrPerCl","Total Charge over N. Cluster",100,0,5e+6);

	TH1F *hBsl_R = new TH1F ("hBsl_R","Base line - R",1000,-0.5,0.5);
	TH1F *hRms_R = new TH1F ("hRms_R","noise RMS - R",500,0,0.05);
	TH1F *hMaxV_R = new TH1F ("hMaxV_R","Max val - R",600,-0.02,0.1);
	TH1F *hMaxVN_R = new TH1F ("hMaxVN_R","Max val over base line - R",300,-0.02,0.1);
	TH1F *hInteg_R = new TH1F ("hInteg_R","Integral - R",1000,-10.,10.);//600,20.,80.
	TH1F *hIntegN_R = new TH1F ("hIntegN_R","Integral minius PDS - R",1000,-10.,10.);
	TH1F *hIntegPerCl_R = new TH1F ("hIntegPerCl_R","Integral per Cluster - R",1000,-10.,10.);//600,20.,80.
	TH1F *hIntegNPerCl_R = new TH1F ("hIntegNPerCl_R","Integral minius PDS per Cluster - R",1000,-10.,10.);
	TH1F *hIntegPerPk_R = new TH1F ("hIntegPerPk_R","Integral per Peak - R",1000,-10.,10.);//600,20.,80.
	TH1F *hIntegNPerPk_R = new TH1F ("hIntegNPerPk_R","Integral minius PDS per Peak - R",1000,-10.,10.);

	TH1F *hBsl_L = new TH1F ("hBsl_L","Base line - R",1000,-0.5,0.5);
	TH1F *hRms_L = new TH1F ("hRms_L","noise RMS - L",500,0,0.05);
	TH1F *hMaxV_L = new TH1F ("hMaxV_L","Max val - L",600,-0.02,0.1);
	TH1F *hMaxVN_L = new TH1F ("hMaxVN_L","Max val over base line - L",300,-0.02,0.1);
	TH1F *hInteg_L = new TH1F ("hInteg_L","Integral - L",1000,-10.,10.);//600,20.,80.
	TH1F *hIntegN_L = new TH1F ("hIntegN_L","Integral minius PDS - L",1000,-10.,10.);
	TH1F *hIntegPerCl_L = new TH1F ("hIntegPerCl_L","Integral per Cluster - L",1000,-10.,10.);//600,20.,80.
	TH1F *hIntegNPerCl_L = new TH1F ("hIntegNPerCl_L","Integral minius PDS per Cluster - L",1000,-10.,10.);
	TH1F *hIntegPerPk_L = new TH1F ("hIntegPerPk_L","Integral per Peak - L",1000,-10.,10.);//600,20.,80.
	TH1F *hIntegNPerPk_L = new TH1F ("hIntegNPerPk_L","Integral minius PDS per Peak - L",1000,-10.,10.);

    TFile simulfile(fName.Data());
    if (simulfile.IsZombie()) {
        std::cout<<"Error opening file"<<std::endl;
        exit(-1);
    }

    TString FlSgnName=fName;
    FlSgnName.ReplaceAll("SimOut","Signal");
    TFile sgnfile(FlSgnName.Data());
    if (sgnfile.IsZombie()) {
        std::cout<<"Error opening file"<<std::endl;
        exit(-1);
    }

    TString FlCnvSgnName=fName;
    FlCnvSgnName.ReplaceAll("SimOut","CnvSignal");
    TFile cvSgnfile(FlCnvSgnName.Data());
    if (cvSgnfile.IsZombie()) {
        std::cout<<"Error opening file"<<std::endl;
        exit(-1);
    }

    TString FlPkName=fName;
    FlPkName.ReplaceAll("SimOut","Peaks");
    TFile pkfile(FlPkName.Data());
    if (pkfile.IsZombie()) {
        std::cout<<"Error opening file"<<std::endl;
        exit(-1);
    }

    TTree * trdata = (TTree *) simulfile.Get("SimulData");
    CellCont *cellCont = new CellCont();
    trdata->SetBranchAddress("DataCell",&cellCont);
    Int_t nEvnt = trdata->GetEntries();

    TTree *trsignal = (TTree*) sgnfile.Get("SimulSignal");
    SignalCont *signalCont= new SignalCont();
    trsignal->SetBranchAddress("SignalDataBlock",&signalCont);

    TTree *trConvsignal = (TTree*) cvSgnfile.Get("ConvSignal");
    ConvSignalCont *convSignalCont= new ConvSignalCont();
    trConvsignal->SetBranchAddress("ConvolutedSignal",&convSignalCont);

    TTree *trPeakGen = (TTree *) pkfile.Get("GenSignalPeak");
    PeakCont *peakContSG = new PeakCont();
	trPeakGen->SetBranchAddress("PeaksData",&peakContSG);
	TTree *trPeakCvR = (TTree *) pkfile.Get("RgtSignalPeak");
	PeakCont *peakContCvR = new PeakCont();
	trPeakCvR->SetBranchAddress("PeaksData",&peakContCvR);
	TTree *trPeakCvL = (TTree *) pkfile.Get("LftSignalPeak");
	PeakCont *peakContCvL = new PeakCont();
	trPeakCvL->SetBranchAddress("PeaksData",&peakContCvL);

    trsignal->GetEntry(0);

    Int_t S_npoint = signalCont->getSignal(swId).GetNWpoints();
    cout<<"NPnt "<<S_npoint<<endl;
    Float_t *xArr=new Float_t [S_npoint];
    for (int i=0; i<S_npoint; ++i){
        xArr[i]=i*timeRes;
    }
    //    for (int is=0; is<100/*S_npoint*/; ++is) {
    //        cout<<xArr[is]<<endl;
    //    }

    cout<<"TotNev "<<nEvnt<<endl;

//	if (trConvsignal->GetEntries()<nEvnt) { nEvnt=trConvsignal->GetEntries(); }
	Int_t fstEv=0;
	if ( (cEv>=0) && (cEv<(nEvnt-1)) ) {
		fstEv=cEv;
		nEvnt=cEv+1;
	}


	for (int iEv=fstEv; iEv<nEvnt; ++iEv) {
		if (iEv%50==0) cout<<"IEv: " <<iEv<<endl;

		trdata->GetEntry(iEv);
		CellData &BlData = cellCont->getCell(swId);
		float nCl = BlData.getNCl();
		hNCl->Fill(nCl);
		hNEl->Fill(BlData.getNCharges());
//		cout<<"iev "<<iEv<<endl;
		float totChr=0.0;
		for (int iel=0; iel<BlData.getNCharges(); ++iel) {
			totChr+=BlData.getCharge(iel);
		}
		hTotChr->Fill(totChr);
		hChrPerCl->Fill(totChr/nCl);

	    trsignal->GetEntry(iEv);
	    trPeakCvR->GetEntry(iEv);
	    trPeakCvL->GetEntry(iEv);
	    SignalDataVl &SignalBlData = signalCont->getSignal(swId);
	    trConvsignal->GetEntry(iEv);
	    ConvSignalDataVl &CSignalBlData = convSignalCont->getSignal(swId);

	    trPeakGen->GetEntry(iEv);
		PeakData &pSG = peakContSG->getPeak(swId);
	    PeakData& pCvR = peakContCvR->getPeak(swId);
	    PeakData& pCvL = peakContCvL->getPeak(swId);

	    float nSG  = pSG.GetNpeaks();
	    float nCvR = pCvR.GetNpeaks();
	    float nCvL = pCvL.GetNpeaks();

		wave wsr, wsl;
	    wsr.fillWave(CSignalBlData.GetSignalR());
	    wsl.fillWave(CSignalBlData.GetSignalL());

	    hBsl_R->Fill(wsr.bsln);
		hRms_R->Fill(wsr.rms);
		hMaxV_R->Fill(wsr.max);
		hMaxVN_R->Fill(wsr.nMax());
		hInteg_R->Fill(wsr.integ);
		hIntegN_R->Fill(wsr.nInteg());
		hIntegPerCl_R->Fill(wsr.integ/nCl);
		hIntegNPerCl_R->Fill(wsr.nInteg()/nCl);
		hIntegPerPk_R->Fill(wsr.integ/nCvR);
		hIntegNPerPk_R->Fill(wsr.nInteg()/nCvR);

		hBsl_L->Fill(wsl.bsln);
		hRms_L->Fill(wsl.rms);
		hMaxV_L->Fill(wsl.max);
		hMaxVN_L->Fill(wsl.nMax());
		hInteg_L->Fill(wsl.integ);
		hIntegN_L->Fill(wsl.nInteg());
		hIntegPerCl_L->Fill(wsl.integ/nCl);
		hIntegNPerCl_L->Fill(wsl.nInteg()/nCl);
		hIntegPerPk_L->Fill(wsl.integ/nCvL);
		hIntegNPerPk_L->Fill(wsl.nInteg()/nCvL);
	}

	simulfile.Close();
	sgnfile.Close();
	cvSgnfile.Close();
	pkfile.Close();

	TCanvas *cv= new TCanvas();
	cv->Divide(2,2);
	cv->cd(1);
	hNCl->Draw();
	cv->cd(2);
	hNEl->Draw();
	cv->cd(3);
	hTotChr->Draw();
	cv->cd(4);
	hChrPerCl->Draw();

	TCanvas *cvR = new TCanvas();
	cvR->Divide(2,3);
	cvR->cd(1);
	hBsl_R->Draw();
	cvR->cd(2);
	hRms_R->Draw();
	cvR->cd(3);
	hMaxV_R->Draw();
	cvR->cd(4);
	hMaxVN_R->Draw();
	cvR->cd(5);
	hInteg_R->Draw();
	cvR->cd(6);
	hIntegN_R->Draw();

	TCanvas *cvL = new TCanvas();
	cvL->Divide(2,3);
	cvL->cd(1);
	hBsl_L->Draw();
	cvL->cd(2);
	hRms_L->Draw();
	cvL->cd(3);
	hMaxV_L->Draw();
	cvL->cd(4);
	hMaxVN_L->Draw();
	cvL->cd(5);
	hInteg_L->Draw();
	cvL->cd(6);
	hIntegN_L->Draw();

	TCanvas *cvRnrm = new TCanvas();
	cvRnrm->Divide(2,2);
	cvRnrm->cd(1);
	hIntegPerCl_R->Draw();
	cvRnrm->cd(2);
	hIntegNPerCl_R->Draw();
	cvRnrm->cd(3);
	hIntegPerPk_R->Draw();
	cvRnrm->cd(4);
	hIntegNPerPk_R->Draw();

	TCanvas *cvLnrm = new TCanvas();
	cvLnrm->Divide(2,2);
	cvLnrm->cd(1);
	hIntegPerCl_L->Draw();
	cvLnrm->cd(2);
	hIntegNPerCl_L->Draw();
	cvLnrm->cd(3);
	hIntegPerPk_L->Draw();
	cvLnrm->cd(4);
	hIntegNPerPk_L->Draw();

}
