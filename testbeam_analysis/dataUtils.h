//contenitore di forme d'onda
#ifndef DATAUTIL_H
#define DATAUTIL_H

#include <vector>
#include <map>

#include "TH1F.h"
#include "TH2F.h"
#include "TMath.h"

static int skipFstBin = 1; //525 for El Cal   //100 dati proto
static int skipLstBin = 1;  //475 for El Cal   //50	dati proto
static int fbin_rms=30;
static int fbin_bsl=30;
static float invfbin=1.0/((float)fbin_rms);

struct wave {
	std::vector<float> Y; //vettore dei dati.
	//variabili di appoggio
	float bsln; //baseline
	float rms;
	float integ;
	float max;
	float integInR;
	float maxInR;
	float min;
	float minInR;
	int   maxPos;
	int   maxInRPos;
	int   minPos;
	int   minInRPos;
	float sumX;
	float sumY;
	float sumXY;
	float sumX2;
	int   nPtReg;
	float regA;
	float regB;
	std::vector<float> deriv;
	std::vector<float> sderiv;

	wave() : bsln(0.0), rms(0.0), integ(0.0), max(-9999), integInR(0.0), maxInR(-9999), min(9999), minInR(9999)
			, maxPos(-1), maxInRPos(-1), minPos(-1), minInRPos(-1)
			, sumX(0.0), sumY(0.0), sumXY(0.0), sumX2(0.0), nPtReg(0), regA(0.0), regB(0.0) {
		Y.clear();
		deriv.clear();
		sderiv.clear();
	}/////inizializzazione

	void clear() {
		Y.clear();
		bsln=0.0;
		rms=0.0;
		integ=0.0;
		max=-9999;
		integInR=0.0;
		maxInR=-9999;
		min=9999;
		minInR=9999;
		maxPos=-1;
		maxInRPos=-1;
		minPos=-1;
		minInRPos=-1;
		sumX=0.0;
		sumY=0.0;
		sumXY=0.0;
		sumX2=0.0;
		regA=0.0;
		regB=0.0;
		deriv.clear();
		sderiv.clear();
	}

	int nPt() { return Y.size(); }
	int nPtInR() { return (Y.size()-(skipFstBin+skipLstBin)); }
	float nInteg() { return (integ-bsln*((float)nPt())); }//funzione float che deve ritornare l'integrale  a cui togliamo la baseline. ////
	float nIntegInR() { return (integInR-bsln*((float)nPtInR())); }//funzione float che deve ritornare l'integrale  a cui togliamo la baseline. ////
	float nMax() { return (max-bsln); }
	float nMaxInR() { return (maxInR-bsln); }
	float nMin() {return (min-bsln);}
	float nMinInR(){return (minInR-bsln);}

	float n1Integ() { return (integ-norm1(0.0,(float)(nPt()-1))); }//funzione float che deve ritornare l'integrale  a cui togliamo la baseline. ////
	float n1IntegInR() { return (integInR-norm1((float)skipFstBin,(float)(nPt()-1-skipLstBin))); }//funzione float che deve ritornare l'integrale  a cui togliamo la baseline. ////
	float n1Max() { return (max-(RegA()+RegB()*((float)maxPos))); }
	float n1MaxInR() { return (maxInR-(RegA()+RegB()*((float)maxInRPos))); }
	float n1Min() {return (min-(RegA()+RegB()*((float)minPos)));}
	float n1MinInR(){return (minInR-(RegA()+RegB()*((float)minInRPos)));}
	float norm1(float x1, float x2) { return (RegA()*(x2-x1)+0.5*RegB()*(x2*x2-x1*x1)); }

	float n1At(int i) { return (Y[i]-(RegA()+RegB()*((float)i))); }
	float nAt(int i) { return (Y[i]-bsln); }
	 
	float nnAt(int i, float cut=0.05) { //0.05 proto_cosmic
		if (nMaxInR()>cut) { return nAt(i); }
		else { return n1At(i); }
	}
	float nnInteg(float cut=0.05) {//0.05 proto_cosmic
		if (nMaxInR()>cut) { return nInteg(); }
		else { return n1Integ(); }
	}
	float nnIntegInR(float cut=0.05) {//0.05 proto_cosmic
		if (nMaxInR()>cut) { return nIntegInR(); }
		else { return n1IntegInR(); }
	}

	float RegA() {
		if (regA==0.0) {
		float N=nPtReg;//skipFstBin+skipLstBin;
		regA=(sumY*sumX2-sumX*sumXY)/(N*sumX2-sumX*sumX);
		}
		return regA;
	}
	float RegB() {
		if (regB==0.0) {
		float N=nPtReg;//skipFstBin+skipLstBin;
		regB=(N*sumXY-sumX*sumY)/(N*sumX2-sumX*sumX);
		}
		return regB;
	}
	//metodi per riempire la waveform
	//void fillWave(int nPt, __uint16_t *arr) {
	void fillWave(std::vector<__uint16_t> &arr, int maxDim=-1) {
		float tmpval;
		int nPt = arr.size();   //dimensione del vettore arr.
		if (maxDim!=-1 && maxDim<nPt) { nPt=maxDim; }
		clear();               //richiama la funzione sopra definita. 
		//riempimento della variabili di appoggio.
		for (int ipt=0; ipt<nPt; ++ipt) {
			tmpval=((1.0/65536.0)*arr[ipt]-0.5);
			tmpval*=-1;
			Y.push_back(tmpval);
			integ+=Y.back();            //l'integrale � in sostanza una somma, integ � inizializzata a zero e poi viene incrementata.
		
			//Rms somma Y nei primi 100 bin
			if (ipt==fbin_rms-1) {
				bsln=integ*invfbin;
				//rms=sqrt(rms*invfbin-bsln*bsln);
				
			}
			if (ipt==fbin_rms/*100*/) {
				rms+=(Y.back()-bsln)*(Y.back()-bsln);    //.back mi restituisce l'ultimo elemento del vettore.
			}

			if (Y.back()>max) { max=Y.back(); maxPos=ipt; }
			if (Y.back()<min) { min=Y.back(); minPos=ipt; }
			
			if (ipt>=skipFstBin && ipt<(nPt-skipLstBin)) {
				integInR+=Y.back();
				if (Y.back()>maxInR) { maxInR=Y.back(); maxInRPos=ipt; }
				if (Y.back()<minInR) { minInR=Y.back(); minInRPos=ipt; }
			} else {
				sumX+=ipt;
				sumX2+=ipt*ipt;
				sumY+=Y.back();
				sumXY+=((float)ipt)*Y.back();
				++nPtReg;
			}
			
		}//fine ciclo sui punti.
		for(int i=0;i<nPt;++i){
	        deriv.push_back((Y[(i+1)>(nPt-1)?(nPt-1):(i+1)]-Y[(i-1)<0?0:(i-1)])/2);
		}
		for(int i=0;i<nPt;++i){
	        sderiv.push_back((deriv[(i+1)>(nPt-1)?(nPt-1):(i+1)]-deriv[(i-1)<0?0:(i-1)])/2);
		}
		rms=sqrt(rms*invfbin);
	}//fine funzione riempimento.
	//void fillWave(int nPt, float *arr) {
	void fillWave(std::vector<float> &arr, int maxDim=-1, bool print=false) {
		int nPt=arr.size();
		if (maxDim!=-1 && maxDim<nPt) { nPt=maxDim; }
		clear();
		//definizione delle variabili d'appoggio
		for (int ipt=0; ipt<nPt; ++ipt) {
			Y.push_back(arr[ipt]);        
			integ+=Y.back();
			if (ipt<fbin_rms) {
				rms+=Y.back()*Y.back();
			}
			if (ipt==fbin_rms-1) {
				bsln=integ*invfbin;
				rms=sqrt(rms*invfbin-bsln*bsln);
			}
			if (Y.back()>max) { max=Y.back(); maxPos=ipt; }
			if (Y.back()<min) { min=Y.back(); minPos=ipt; }
			if (ipt>=skipFstBin && ipt<(nPt-skipLstBin)) {
				integInR+=Y.back();
				if (Y.back()>maxInR) { maxInR=Y.back(); maxInRPos=ipt; }
				if(Y.back()<minInR)  { minInR=Y.back(); minInRPos=ipt; }
			} else {
				sumX+=ipt;
				sumX2+=ipt*ipt;
				sumY+=Y.back();
				sumXY+=((float)ipt)*Y.back();
				++nPtReg;
			}
//			if (print) std::cout<<"ipt "<<ipt<<" pnt "<<Y.back()<<" integ "<<integ<<std::endl;
		}
		for(int i=0;i<nPt;++i){
	            deriv.push_back((Y[(i+1)>(nPt-1)?(nPt-1):(i+1)]-Y[(i-1)<0?0:(i-1)])/2);
		}
		for(int i=0;i<nPt;++i){
	        sderiv.push_back((deriv[(i+1)>(nPt-1)?(nPt-1):(i+1)]-deriv[(i-1)<0?0:(i-1)])/2);
		}
	}

	void addPnt(float &val, bool InRng=true) {
		Y.push_back(val);
		int ipt=nPt()-1;
		integ+=val;

		if (ipt<fbin_rms) {
			rms+=val*val;
		}
		if (ipt==fbin_rms-1) {
			bsln=integ*invfbin;
			rms=sqrt(rms*invfbin-bsln*bsln);
		}
		if (val>max) { max=val; maxPos=ipt; }
		if (val<min) { min=val; minPos=ipt; }
		if (ipt>=skipFstBin && InRng) {
			integInR+=Y.back();
			if (Y.back()>maxInR) { maxInR=Y.back(); maxInRPos=ipt; }
			if (Y.back()<minInR) { minInR=Y.back(); minInRPos=ipt; }
		} else {
			sumX+=ipt;
			sumX2+=ipt*ipt;
			sumY+=Y.back();
			sumXY+=((float)ipt)*Y.back();
			++nPtReg;
		}
	}
}; //fine della struct wave


typedef std::map<int,wave> WvCont;                    //typedef assegna un altro nome al tipo specificato.
typedef std::map<std::pair<int,int>,wave> diffWvCont; //serve per fare la differenza
	


//struttura per costruire istogrammi
struct hstPerCh { //istogrammi per tutti i canali dell'oscilloscopio.
//isto che verranno riempiti con wave filtrate
	TH1F *hSum;
	TH1F *hHPeaks;
	TH2F *hHNPeaks;
	TH1F *hTPeaks;
	TH1F *hTFstPeaks;
	TH1F *hNPeaks_1;
	TH1F *hBsl;
	TH1F *hInteg;
	TH1F *hIntegN;
	TH1F *hIntegInR;
	TH1F *hIntegNInR;
	TH1F *hIntegNInRC1;
	TH1F *hIntegNInRC2;
	TH1F *hIntegNInRC3;
	TH1F *hIntegNInRC4;

	TH1F *hNeventSignals;
	TH1F *hNPeaks;	


	TH1F *hRms;
	TH1F *hMaxV;
	TH1F *hMaxVN;
	TH1F *hMaxVInR;
	TH1F *hMaxVNInR;
 	TH1F *hMinV;
	TH1F *hMinVN;
	TH1F *hMinVInR;
	TH1F *hMinVNInR;
//histo che verranno riempiti con wave non filtrate 
	TH1F *hIntegInRoriginalW;
	TH1F *hIntegNInRoriginalW;
	TH1F *hMaxVoriginalW;
	TH1F *hMaxVNoriginalW;
	TH1F *hMaxVInRoriginalW;
	TH1F *hMaxVNInRoriginalW;
	TH1F *hMinVoriginalW;
	TH1F *hMinVNoriginalW;
	TH1F *hMinVInRoriginalW;
	TH1F *hMinVNInRoriginalW;
	TH1F *hder;
	TH1F *hchiFT;
	TH1F *hP0;
	TH1F *hRmsOriginalW;
	TH1F *hMaxVNSmooth;

	TH1F *hIntegNInRFullW;
	//Derivative study
	TH1F *hFirstDeriv;	
    TH1F *hSecDeriv;

	hstPerCh(int Ch=0, int SF=0) {
		//gRootDir->cd();
		//TDirectory fldch(Form("H-Ch%d",Ch),Form("folder for ch%d",Ch));
		//fldch.cd();
		if (SF==0) {//senza smooth
			hSum = new TH1F (Form("hSum_ch%d",Ch),Form("Sum amplitudes - Ch %d",Ch),1000,-1000,1000.); 
			hP0=new TH1F (Form("hP0_ch%d",Ch),Form("HP0 - Ch %d",Ch),1000,-1,1);
			hder = new TH1F (Form("hder_ch%d",Ch),Form("Hder - Ch %d",Ch),1000,-1e+7,1e+7);
			hchiFT = new TH1F (Form("hchiFT_ch%d",Ch),Form("HchiFT - Ch %d",Ch),1000,-10,10);
			hNPeaks = new TH1F (Form("hNPeaks_ch%d",Ch),Form("N Peaks found - Ch %d",Ch),100,-0.5,99.5);
			hNeventSignals = new TH1F (Form("hNeventSignals_ch%d",Ch),Form("N event Signals - Ch %d",Ch),2,-0.5,1.5); 

			hHPeaks = new TH1F (Form("hHPeaks_ch%d",Ch),Form("Height of Peaks found - Ch %d",Ch),500,0,0.3);
			hHNPeaks = new TH2F (Form("hHNPeaks_ch%d",Ch),Form("Height vs N of Peaks found - Ch %d",Ch),100,0,100,500,0,0.5);
			hTPeaks = new TH1F (Form("hTPeaks_ch%d",Ch),Form("Time of Peaks found - Ch %d",Ch),500,0,1000);
			hTFstPeaks = new TH1F (Form("hTFstPeaks_ch%d",Ch),Form("Time of First Peak found - Ch %d",Ch),500,0,1000);
			hNPeaks_1 = new TH1F (Form("hNPeaks_1_ch%d",Ch),Form("Time of Last Peak found - Ch %d",Ch),500,0,1000);


			//hFirstDeriv= new TH1F (Form("hFirstDeriv_ch%d",Ch),Form("First derivative- Ch %d",Ch),10000,-0.01,0.01);
			//hSecDeriv= new TH1F (Form("hSecDeriv_ch%d",Ch),Form("Second derivative- Ch %d",Ch),10000,-0.01,0.01);
	    


			hBsl = new TH1F (Form("hBsl_ch%d",Ch),Form("Base line - Ch %d",Ch),1000,-0.5,0.2); 
			hInteg = new TH1F (Form("hInteg_ch%d",Ch),Form("Integral - Ch %d",Ch),10000,-0.5,0.5);//600,20.,80.
			hIntegN = new TH1F (Form("hIntegN_ch%d",Ch),Form("Integral minius PDS - Ch %d",Ch),1000,-10.,10.);
			hIntegInR = new TH1F (Form("hIntegInR_ch%d",Ch),Form("Integral - Ch %d",Ch),10000,-0.01,1);
			hIntegNInR = new TH1F (Form("hIntegNInR_ch%d",Ch),Form("Integral minius PDS - Ch %d",Ch),1000,-10.,10.);
			hIntegNInRC1 = new TH1F (Form("hIntegNInRC1_ch%d",Ch),Form("Integral minius PDS Norm. on NPeak - Ch %d",Ch),1000,-10.,10.);
			hIntegNInRC2 = new TH1F (Form("hIntegNInRC2_ch%d",Ch),Form("Integral minius PDS Norm. on NPeak and loss - Ch %d",Ch),1000,-10.,10.);
			hIntegNInRC3 = new TH1F (Form("hIntegNInRC3_ch%d",Ch),Form("Integral minius PDS whitout norm - Ch %d",Ch),1000,-10.,20.);
			hIntegNInRC4 = new TH1F (Form("hIntegNInRC4_ch%d",Ch),Form("Integral minius PDS Norm. on loss - Ch %d",Ch),1000,-10.,20.);
			hIntegNInRFullW=new TH1F (Form("hIntegNInRFullW_ch%d",Ch),Form("Integral minius PDS for Full Waves - Ch %d",Ch),1000,-10.,10.);
			hRms = new TH1F (Form("hRms_ch%d",Ch),Form("noise RMS - Ch %d",Ch),500,0,0.005);
			//distribution of derivative for meg FE

			hMaxV = new TH1F (Form("hMaxV_ch%d",Ch),Form("Max val - Ch %d",Ch),600,-0.02,0.3);
			hMaxVN = new TH1F (Form("hMaxVN_ch%d",Ch),Form("Max val over base line - Ch %d",Ch),300,-0.02,0.3);
			hMaxVNSmooth= new TH1F (Form("hMaxVNSmooth_ch%d",Ch),Form("Max val Smooth over base line - Ch %d",Ch),300,-0.02,0.3);
			hMaxVInR = new TH1F (Form("hMaxVInR_ch%d",Ch),Form("Max val - Ch %d",Ch),600,-0.02,0.3);
			hMaxVNInR = new TH1F (Form("hMaxVNInR_ch%d",Ch),Form("Max val over base line - Ch %d",Ch),300,-0.02,0.3);

			hMinV = new TH1F (Form("hMinV_ch%d",Ch),Form("Min val - Ch %d",Ch),400,-0.02,0);
			hMinVN = new TH1F (Form("hMinVN_ch%d",Ch),Form("Min val over base line - Ch %d",Ch),400,-0.02,0);
			hMinVInR = new TH1F (Form("hMinVInR_ch%d",Ch),Form("Min val - Ch %d",Ch),400,-0.02,0);
			hMinVNInR = new TH1F (Form("hMinVNInR_ch%d",Ch),Form("Min val over base line - Ch %d",Ch),400,-0.02,0);

			hRmsOriginalW = new TH1F (Form("hRmsOriginalW_ch%d",Ch),Form("noise RMS OriginalW- Ch %d",Ch),500,0,0.05);
			hIntegInRoriginalW = new TH1F (Form("hIntegInRoriginalW_ch%d",Ch),Form("Integral - Ch %d",Ch),1000,-10.,10.);
			hIntegNInRoriginalW = new TH1F (Form("hIntegNInRoriginalW_ch%d",Ch),Form("Integral minius PDS - Ch %d",Ch),1000,-10.,10.);

			hMinVoriginalW = new TH1F (Form("hMinVoriginalW_ch%d",Ch),Form("Min val for Original Wave - Ch %d",Ch),400,-0.02,0);
			hMinVNoriginalW = new TH1F (Form("hMinVNoriginalW_ch%d",Ch),Form("Min val over base line for Original Wave - Ch %d",Ch),400,-0.02,0);
			hMinVInRoriginalW = new TH1F (Form("hMinVInRoriginalW_ch%d",Ch),Form("Min val for Original Wave - Ch %d",Ch),400,-0.02,0);
			hMinVNInRoriginalW = new TH1F (Form("hMinVNInRoriginalW_ch%d",Ch),Form("Min val over base line for Original Wave- Ch %d",Ch),400,-0.02,0);

			hMaxVoriginalW = new TH1F (Form("hMaxVoriginalW_ch%d",Ch),Form("Max val for Original Wave - Ch %d",Ch),600,-0.02,0.1);
			hMaxVNoriginalW = new TH1F (Form("hMaxVNoriginalW_ch%d",Ch),Form("Max val over base line for Original Wave - Ch %d",Ch),300,-0.02,0.1);
			hMaxVInRoriginalW = new TH1F (Form("hMaxVInRoriginalW_ch%d",Ch),Form("Max val for Original Wave - Ch %d",Ch),600,-0.02,0.1);
			hMaxVNInRoriginalW = new TH1F (Form("hMaxVNInRoriginalW_ch%d",Ch),Form("Max val over base line for Original Wave- Ch %d",Ch),300,-0.02,0.1);
			//////////////////////////////////////
			hNPeaks->GetYaxis()->SetTitle("Entries");
			hNeventSignals->GetYaxis()->SetTitle("Entries");
        	hHPeaks->GetYaxis()->SetTitle("Entries");
			hHPeaks->GetXaxis()->SetTitle("Height [V]");
        	hHNPeaks->GetYaxis()->SetTitle("Height of Peaks found");
			hHNPeaks->GetXaxis()->SetTitle("Number of Peaks found");
        	hTPeaks->GetYaxis()->SetTitle("Entries");
        	hTFstPeaks->GetYaxis()->SetTitle("Entries");
			hTFstPeaks->GetXaxis()->SetTitle("Time [ns]");
        	hNPeaks_1->GetYaxis()->SetTitle("Entries");
        	hBsl->GetYaxis()->SetTitle("Entries");
        	hInteg->GetYaxis()->SetTitle("Entries");
        	hIntegN->GetYaxis()->SetTitle("Entries");
        	hIntegInR->GetYaxis()->SetTitle("Entries");
        	hIntegNInR->GetYaxis()->SetTitle("Entries");
        	hIntegNInRC1->GetYaxis()->SetTitle("Entries");
			hIntegNInRC2->GetYaxis()->SetTitle("Entries");
			hIntegNInRC3->GetYaxis()->SetTitle("Entries");
			hIntegNInRC4->GetYaxis()->SetTitle("Entries");
			hRms->GetYaxis()->SetTitle("Entries");
			hMaxV->GetYaxis()->SetTitle("Entries");
			hMaxVN->GetYaxis()->SetTitle("Entries");
			hMaxVInR->GetYaxis()->SetTitle("Entries");
			hMaxVNInR->GetYaxis()->SetTitle("Entries");
			hMinV->GetYaxis()->SetTitle("Entries");
			hMinVN->GetYaxis()->SetTitle("Entries");
			hMinVInR->GetYaxis()->SetTitle("Entries");
			hMinVNInR->GetYaxis()->SetTitle("Entries");
	        hIntegInRoriginalW->GetYaxis()->SetTitle("Entries");
	        hIntegNInRoriginalW->GetYaxis()->SetTitle("Entries");
	        hMaxVoriginalW->GetYaxis()->SetTitle("Entries");
	        hMaxVNoriginalW->GetYaxis()->SetTitle("Entries");
	        hMaxVInRoriginalW->GetYaxis()->SetTitle("Entries");
	        hMaxVNInRoriginalW->GetYaxis()->SetTitle("Entries");
	        hMinVoriginalW->GetYaxis()->SetTitle("Entries");
	        hMinVNoriginalW->GetYaxis()->SetTitle("Entries");
	        hMinVInRoriginalW->GetYaxis()->SetTitle("Entries");
	        hMinVNInRoriginalW->GetYaxis()->SetTitle("Entries");
			hder->GetYaxis()->SetTitle("Entries");
			hchiFT->GetYaxis()->SetTitle("Entries");
			hP0->GetYaxis()->SetTitle("Entries");
			hRmsOriginalW->GetYaxis()->SetTitle("Entries");
			hMaxVNSmooth->GetYaxis()->SetTitle("Entries");


			/////////////////////////////////////////////////////////
			hNPeaks->GetXaxis()->SetTitle("N Peaks");
			hNPeaks->GetYaxis()->SetTitle("Entries");
			hBsl->GetXaxis()->SetTitle("Baseline [V]");
			hBsl->GetYaxis()->SetTitle("Entries");
			hIntegNInR-> GetXaxis()->SetTitle("Integral [V]");
			hIntegNInR-> GetYaxis()->SetTitle("Entries");
			hIntegNInRC1-> GetXaxis()->SetTitle("Integral [V]");
			hIntegNInRC1-> GetYaxis()->SetTitle("Entries");
			hIntegNInRC2-> GetXaxis()->SetTitle("Integral [V]");
			hIntegNInRC2-> GetYaxis()->SetTitle("Entries");
			hRms->GetXaxis()->SetTitle("Rms [V]");
			hRms->GetYaxis()->SetTitle("Entries");

			hMaxVNInR->GetXaxis()->SetTitle("Max_value [V]");
			hMaxVNInR->GetYaxis()->SetTitle("Entries");
			hMinVNInR->GetYaxis()->SetTitle("Entries");
			hMinVNInR->GetXaxis()->SetTitle("Min value [V]");

			hIntegNInRoriginalW-> GetXaxis()->SetTitle("Integral [mA]");
			hIntegNInRoriginalW-> GetYaxis()->SetTitle("Entries");

			hMinVNInRoriginalW->GetXaxis()->SetTitle("Min value [V]");
			hMinVNInRoriginalW->GetYaxis()->SetTitle("Entries");
			hMaxVNInRoriginalW->GetXaxis()->SetTitle("Min_value [V]");
			hMaxVNInRoriginalW->GetYaxis()->SetTitle("Entries");



		} else {//con smooth
			hBsl = new TH1F (Form("hBsl_SF%d_ch%d",SF,Ch),Form("Base line Smoothed Wave with SF %d - Ch %d",SF,Ch),1000,-0.5,0.5);
			hInteg = new TH1F (Form("hInteg_SF%d_ch%d",SF,Ch),Form("Integral Smoothed Wave with SF %d - Ch %d",SF,Ch),1000,-10.,10.);
			hIntegN = new TH1F (Form("hIntegN_SF%d_ch%d",SF,Ch),Form("Integral minius PDS Smoothed Wave with SF %d - Ch %d",SF,Ch),1000,-10.,10.);
			hIntegInR = new TH1F (Form("hIntegInR_SF%d_ch%d",SF,Ch),Form("Integral Smoothed Wave with SF %d - Ch %d",SF,Ch),1000,-10.,10.);
			hIntegNInR = new TH1F (Form("hIntegNInR_SF%d_ch%d",SF,Ch),Form("Integral minius PDS Smoothed Wave with SF %d - Ch %d",SF,Ch),1000,-10.,10.);
			hRms = new TH1F (Form("hRms_SF%d_ch%d",SF,Ch),Form("noise RMS Smoothed Wave with SF %d - Ch %d",SF,Ch),500,0,0.05);
			hMaxV = new TH1F (Form("hMaxV_SF%d_ch%d",SF,Ch),Form("max val Smoothed Wave with SF %d - Ch %d",SF,Ch),600,-0.02,0.1);
			hMaxVN = new TH1F (Form("hMaxVN_SF%d_ch%d",SF,Ch),Form("max val over base line Smoothed Wave with SF %d - Ch %d",SF,Ch),300,-0.02,0.1);//,1000,-0.02,0.02);
			hMaxVInR = new TH1F (Form("hMaxVInR_SF%d_ch%d",SF,Ch),Form("max val Smoothed Wave with SF %d - Ch %d",SF,Ch),600,-0.02,0.1);
			hMaxVNInR = new TH1F (Form("hMaxVNInR_SF%d_ch%d",SF,Ch),Form("max val over base line Smoothed Wave with SF %d - Ch %d",SF,Ch),300,-0.02,0.1);//,1000,-0.02,0.02);
		}
	}

	~hstPerCh() {
//		delete hBsl;
//		delete hInteg;
//		delete hIntegN;
//		delete hRms;
//		delete hMaxV;
//		delete hMaxVN;
	}
};

struct hstDiffCh {//non serve

	TH1F *hBsl;
	TH1F *hInteg;
	TH1F *hIntegN;
	TH1F *hRms;
	TH1F *hMaxV;
	TH1F *hMaxVN;

	hstDiffCh(int Ch1=0, int Ch2=1, bool ped=false) {
		//gRootDir->cd();
		//TDirectory fldch(Form("H-Ch%d",Ch),Form("folder for ch%d",Ch));
		//fldch.cd();
		std::string pfix="";
		if (ped) { pfix="_ped"; }
		hBsl = new TH1F (Form("hBsl_ch%d-ch%d%s",Ch1,Ch2,pfix.c_str()),Form("Base line - Ch %d -Ch %d",Ch1,Ch2),1000,-0.5,0.5);
		hInteg = new TH1F (Form("hInteg_ch%d-ch%d%s",Ch1,Ch2,pfix.c_str()),Form("Integral - Ch %d -Ch %d",Ch1,Ch2),600,20.,80.);
		hIntegN = new TH1F (Form("hIntegN_ch%d-ch%d%s",Ch1,Ch2,pfix.c_str()),Form("Integral minius PDS - Ch %d -Ch %d",Ch1,Ch2),1000,-10.,10.);
		hRms = new TH1F (Form("hRms_ch%d-ch%d%s",Ch1,Ch2,pfix.c_str()),Form("noise RMS - Ch %d -Ch %d",Ch1,Ch2),500,0,0.05);
		hMaxV = new TH1F (Form("hMaxV_ch%d-ch%d%s",Ch1,Ch2,pfix.c_str()),Form("max val - Ch %d -Ch %d",Ch1,Ch2),600,-0.02,0.1);
		hMaxVN = new TH1F (Form("hMaxVN_ch%d-ch%d%s",Ch1,Ch2,pfix.c_str()),Form("max val over base line - Ch %d -Ch %d",Ch1,Ch2),300,-0.02,0.1);
	}

	~hstDiffCh() {
//		delete hBsl;
//		delete hInteg;
//		delete hIntegN;
//		delete hRms;
//		delete hMaxV;
//		delete hMaxVN;
	}

};//fine della struct per gli istogrammi.

typedef std::map<std::pair<int,int>, hstDiffCh> hDiffCont;
//struttura per fare la trasformata di Fourier
struct wavefft {
	int nPt;
	double *real;
	double *img;
	double *ampl;
	double *phi;
	double *omega;

	wavefft(int npt=0) : nPt(npt) {   
		if (npt>0) {
			real = new double[npt];      //new alloca zone di memoria quando non si sa quanta memoria ci serve. Il suo complemento � delete
			img  = new double[npt];
			ampl = new double[npt];
			phi  = new double[npt];
			omega = new double[npt];
		}
		else { real = img = ampl = phi = omega =0x0; }//il puntatore punta a zero:0x0
	}

	~wavefft() {
		if (real != 0x0 ) { delete real; }
		if (img  != 0x0 ) { delete img; }
		if (ampl != 0x0 ) { delete ampl; }
		if (phi  != 0x0 ) { delete phi; }
		if (omega!= 0x0 ) { delete omega; }
	}

	void setSize(int npt=0) {
		if (npt>0 && npt!=nPt) {
			if (real != 0x0 ) { delete real; }
			if (img  != 0x0 ) { delete img; }
			if (ampl != 0x0 ) { delete ampl; }
			if (phi  != 0x0 ) { delete phi; }
			if (omega!= 0x0 ) { delete omega; }
			nPt=npt;
			real = new double[npt];
			img  = new double[npt];
			ampl = new double[npt];
			phi  = new double[npt];
			omega = new double[npt];
		}
	}

	void clear() {
		for (int iPt=0; iPt<nPt; ++iPt) {
			real[iPt]=0.0;
			img[iPt]=0.0;
			ampl[iPt]=0.0;
			phi[iPt]=0.0;
			omega[iPt]=0.0;
		}
	}
 
  void fillRM(int npt, double *Real, double *Img, double *Omega) {
    if (nPt!=npt) {setSize(npt);}
    for (int iPt=0; iPt<nPt; ++iPt) {
      real[iPt]=Real[iPt];
      img[iPt]=Img[iPt];
      ampl[iPt]=TMath::Sqrt(Real[iPt]*Real[iPt]+Img[iPt]*Img[iPt]);
      phi[iPt]=(Img[iPt]||Real[iPt])?TMath::ATan2(Img[iPt],Real[iPt]):0.0;
      omega[iPt]=Omega[iPt];
      //if (Omega!=0x0) { omega[iPt]=Omega[iPt];}
      //else { omega[iPt]=TMath::TwoPi()*((Double_t) (iPt+1))/tmax; }
    }
  }

  void fillAP(int npt, double *Ampl, double *Phi, double *Omega) {
    if (nPt!=npt) {setSize(npt);}
    for (int iPt=0; iPt<nPt; ++iPt) {
      ampl[iPt]=Ampl[iPt];
      phi[iPt]=Phi[iPt];
      omega[iPt]=Omega[iPt];
      real[iPt]=ampl[iPt]*TMath::Cos(phi[iPt]);
      img[iPt]=ampl[iPt]*TMath::Sin(phi[iPt]);
    }
  }

};//fine struttura wavefft.
typedef std::map<int,wavefft> fftCont; //mappa per le forme d'onda filtrate.


#endif
