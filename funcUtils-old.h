#ifndef FUNCUTIL_H
#define FUNCUTIL_H

#include "dataUtils.h"
#include "TVirtualFFT.h"
#include "TMath.h"

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

                    //////////////Define Functions//////////////

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

std::vector<float> smooth(std::vector<float> &arr, int sf=1) {
  int vsz=arr.size();
  float ssz=2*sf+1;
  std::vector<float> out(0);
  for (int i=sf; i<vsz-sf; ++i) {
    float tmpv=0;
    for (int j=-sf; j<=sf; ++j) {
      tmpv+=arr[i+j];
    }
    tmpv/=ssz;
    out.push_back(tmpv);
  }
  return out;
}

static double tmax=0.0;

void FFT(wave &tmpW, wavefft &tmpWfft) {
//double tmax = 1.0e-6;
   Int_t N = tmpW.nPt();

   TVirtualFFT *fft = TVirtualFFT::FFT(1,&N,"R2C");
   for (int ipt=0; ipt<N; ++ipt) { fft->SetPoint(ipt,tmpW.Y[ipt]); }
   fft->Transform();

   tmpWfft.setSize(N);
//   fft->GetPointsComplex(tmpWfft.real,tmpWfft.img);
   Int_t index;
   for (index=0; index<(N*0.5); ++index) {
     fft->GetPointComplex(index, tmpWfft.real[index], tmpWfft.img[index]);
 //     tmpWfft.real[index]*=1./TMath::Sqrt(2.*pi);
 //     tmpWfft.img[index]*=1./TMath::Sqrt(2.*pi);
     tmpWfft.ampl[index]=TMath::Sqrt(pow(tmpWfft.real[index],2)+pow(tmpWfft.img[index],2)); ///TMath::Sqrt(2.*pi);
//     tmpWfft.phi[index]=TMath::ATan(tmpWfft.img[index]/tmpWfft.real[index]);
     tmpWfft.phi[index]=(tmpWfft.img[index]||tmpWfft.real[index])?TMath::ATan2(tmpWfft.img[index],tmpWfft.real[index]):0.0;
 //    cout<<re[index]<<" "<<im[index]<<" "<<amp[index]<<endl;
     tmpWfft.omega[index] = TMath::TwoPi()*((Double_t) (index+1))/tmax;
   }
   for (index=N-1; index>=(N*0.5); --index) {
     fft->GetPointComplex(index, tmpWfft.real[index], tmpWfft.img[index]);
 //     tmpWfft.real[index]*=1./TMath::Sqrt(2.*pi);
 //     tmpWfft.img[index]*=1./TMath::Sqrt(2.*pi);
     tmpWfft.ampl[index]=TMath::Sqrt(pow(tmpWfft.real[index],2)+pow(tmpWfft.img[index],2)); ///TMath::Sqrt(2.*pi);
//     tmpWfft.phi[index]=TMath::ATan(tmpWfft.img[index]/tmpWfft.real[index]);
     tmpWfft.phi[index]=(tmpWfft.img[index]||tmpWfft.real[index])?TMath::ATan2(tmpWfft.img[index],tmpWfft.real[index]):0.0;
 //    cout<<re[index]<<" "<<im[index]<<" "<<amp[index]<<endl;
     tmpWfft.omega[index] = TMath::TwoPi()*((Double_t) (N-index+1))/tmax;
   }

}

void InverseFFT(Double_t *ReFFT, Double_t *ImFFT, Int_t N, wave &outW) {

   TVirtualFFT *ifft = TVirtualFFT::FFT(1,&N,"C2R");
   ifft->SetPointsComplex(ReFFT,ImFFT);
   ifft->Transform();

   for (Int_t i=0; i<N; i++) {
	   float tmpVal = ifft->GetPointReal(i)/((Double_t)N);
	   outW.addPnt(tmpVal);
   }

}

void filterWave(wavefft &inWfft, Double_t *ReFFT, Double_t *ImFFT, double minCut=50e+6, double maxCut=1500e+6){
	double aveMin=0.0;
	int nMin=0;
	double aveMax=0.0;
	int nMax=0;
	for (int ip=0; ip<inWfft.nPt/2; ++ip) {
		if(inWfft.omega[ip]>minCut*0.9 && inWfft.omega[ip]<minCut*1.1) {
			aveMin+=inWfft.ampl[ip];
			++nMin;
		}
		if(inWfft.omega[ip]>maxCut*0.9 && inWfft.omega[ip]<maxCut*1.1) {
			aveMax+=inWfft.ampl[ip];
			++nMax;
		}
	}
	aveMin/=((double)nMin);
	aveMax/=((double)nMax);
//	std::cout<<"AveMin "<<aveMin<<" aveMax "<<aveMax<<std::endl;
	for (int ip=0; ip<inWfft.nPt/2; ++ip) {
		if (inWfft.omega[ip]>minCut && inWfft.omega[ip]<maxCut) {
			double newAmp = (aveMax-aveMin)/(maxCut-minCut)*(inWfft.omega[ip]-minCut) + aveMin;
//			std::cout<<"NewAmp "<<newAmp<<std::endl;
			ReFFT[ip]=newAmp*TMath::Cos(inWfft.phi[ip]);
			ImFFT[ip]=newAmp*TMath::Sin(inWfft.phi[ip]);
//			ReFFT[ip]=inWfft.ampl[ip]*TMath::Cos(inWfft.phi[ip]);
//			ImFFT[ip]=inWfft.ampl[ip]*TMath::Sin(inWfft.phi[ip]);
		} else {
			ReFFT[ip]=inWfft.real[ip];
			ImFFT[ip]=inWfft.img[ip];
//			ReFFT[ip]=inWfft.ampl[ip]*TMath::Cos(inWfft.phi[ip]);
//			ImFFT[ip]=inWfft.ampl[ip]*TMath::Sin(inWfft.phi[ip]);
		}
	}
	aveMin=0.0;
	nMin=0;
	aveMax=0.0;
	nMax=0;
	for (int ip=inWfft.nPt-1; ip>=inWfft.nPt/2; --ip) {
		if(inWfft.omega[ip]>minCut*0.9 && inWfft.omega[ip]<minCut*1.1) {
			aveMin+=inWfft.ampl[ip];
			++nMin;
		}
		if(inWfft.omega[ip]>maxCut*0.9 && inWfft.omega[ip]<maxCut*1.1) {
			aveMax+=inWfft.ampl[ip];
			++nMax;
		}
	}
	aveMin/=((double)nMin);
	aveMax/=((double)nMax);
	for (int ip=inWfft.nPt-1; ip>=inWfft.nPt/2; --ip) {
		if (inWfft.omega[ip]>minCut && inWfft.omega[ip]<maxCut) {
			double newAmp = (aveMax-aveMin)/(maxCut-minCut)*(inWfft.omega[ip]-minCut) + aveMin;
			ReFFT[ip]=0.0;//newAmp*TMath::Cos(inWfft.phi[ip]);
			ImFFT[ip]=0.0;//newAmp*TMath::Sin(inWfft.phi[ip]);
//			ReFFT[ip]=inWfft.ampl[ip]*TMath::Cos(inWfft.phi[ip]);
//			ImFFT[ip]=inWfft.ampl[ip]*TMath::Sin(inWfft.phi[ip]);
		} else {
			ReFFT[ip]=0.0;//inWfft.real[ip];
			ImFFT[ip]=0.0;//inWfft.img[ip];
//			ReFFT[ip]=inWfft.ampl[ip]*TMath::Cos(inWfft.phi[ip]);
//			ImFFT[ip]=inWfft.ampl[ip]*TMath::Sin(inWfft.phi[ip]);
		}
	}

}

#endif
