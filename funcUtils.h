#ifndef FUNCUTIL_H
#define FUNCUTIL_H

#include "dataUtils.h"
#include "TVirtualFFT.h"
#include "TMath.h"
#include "iostream"
#include "TSpectrum.h"
#include "TComplex.h"


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

                    //////////////Define Functions//////////////

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///----------smooth------------
std::vector<float> smooth(float *arr, int vsz=0, int sf=1) {
  float ssz=2*sf+1;
  std::vector<float> out(0);
  for (int i=sf; i<vsz-sf; ++i) {
    float tmpv=0;//valore d'appoggio su cui fare media
    for (int j=-sf; j<=sf; ++j) {
      tmpv+=arr[i+j];
    }
    tmpv/=ssz;
    out.push_back(tmpv);
  }
  return out;
}

std::vector<float> smooth(std::vector<float> &arr, int sf=1) {
  return smooth(arr.data(),arr.size(),sf);
}
///////////////////////////////SG//////////////////////////
namespace smtSG {
  int norm_1=35;
  int norm_2=21;
  int norm_3=231;
  int norm_4=231;
  int norm_5=429;
  int norm_6=143;
  int norm_7=805;
  int C_1[5]={-3,12,17,12,-3};//k=2
  int C_2[7]={-2,3,6,7,6,3,-2};//k=2
  int C_3[9]={-21,14,39,54,59,54,39,14,-21};//k=2
  int C_4[7]={5,-30,75,131,75,-30,5}; //k=4
  int C_5[9]={15,-55,30,135,179,135,30,-55,15}; //k=4
  int C_6[13]={-11,0,9,16,21,24,25,24,21,16,9,0,-11};//k=2
  int C_7[23]={-42,-21,-2,15,30,43,54,63,70,75,78,79,78,75,70,63,54,43,30,15,-2,-21,-42};//k=2
}

void getSGpars(int **C, float &norm, int m=5, int k=2) {
    if (k==0) {
    	int *Cf = new int[m];
    	for (int i=0; i<m; ++i) {Cf[i]=1;}
    	*C=Cf;
    	norm=m; //media mobile.
    } else if(k==2||k==3){
    	if (m==5) {
    		*C=smtSG::C_1;
    		norm=smtSG::norm_1;
    	} else if (m==7) {
    		*C=smtSG::C_2;
    		norm=smtSG::norm_2;
    	} else if (m==9) {
    		*C=smtSG::C_3;
    		norm=smtSG::norm_3;
    	}else if (m==13) {
    		*C=smtSG::C_6;
    		norm=smtSG::norm_6;
    	}else if (m==23) {
    		*C=smtSG::C_7;
    		norm=smtSG::norm_7;
    	}
    } else if (k==4||k==5) {
    	if (m==9) {
    		*C=smtSG::C_5;
    		norm=smtSG::norm_5;
    	}else if (m==7) {
    		*C=smtSG::C_4;
    		norm=smtSG::norm_4;
    	}
    }
}
/*-----------------SMOOTHING SAVITZKY-GOLAY-------------------------------------------------------*/
std::vector<float> smoothSG(double *arr, int vd=0, int m=5, int k=2){ //k il grado della polinomiale.
	//    int vd=arr.size(); //dimensione dell'array (vector dimension)
	//    float n = (m-1)/2; //bin su cui vado a fare lo smoothing.
	int n = (m-1)/2;
	std::vector<float> outSG(0);

	int *C=0x0;
	float norm=1.0;
	getSGpars(&C,norm,m,k);

	for (int i=n;i<=vd-n;++i){//m 
		float tmpvSG=0.0;//variabile temporanea
		for(int j=-n;j<=n;++j){
			tmpvSG+= C[j+n]*arr[i+j];
		}
		tmpvSG=tmpvSG/norm;
		outSG.push_back(tmpvSG);
	}
	if (k==0) {delete [] C;}
	return outSG;
}
////////////////////////////////////////////////////////////-----------------
std::vector<float> smoothSG(float *arr, int vd=0, int m=5, int k=2){ //k il grado della polinomiale.
//    int vd=arr.size(); //dimensione dell'array (vector dimension)
//    float n = (m-1)/2; //bin su cui vado a fare lo smoothing.
    int n = (m-1)/2;
    std::vector<float> outSG(0);

    int *C=0x0;
    float norm=1.0;
    getSGpars(&C,norm,m,k);
//	std::cout<<"C "<<C<<" norm "<<norm<<std::endl;

	for (int i=n;i<=vd-n;++i){
		float tmpvSG=0.0;//variabile temporanea
		for(int j=-n;j<=n;++j){
			tmpvSG+= C[j+n]*arr[i+j];
		}
		tmpvSG=tmpvSG/norm;
		outSG.push_back(tmpvSG);
	}
	if (k==0) {delete [] C;}
	return outSG;
}

std::vector<float> smoothSG(std::vector<float> &arr, int m=5, int k=2){ //k il grado della polinomiale.
  return smoothSG(arr.data(),arr.size(),m,k);
}

//---------------------------

static double tmax=0.0;

void FFT(wave &tmpW, wavefft &tmpWfft) {
//double tmax = 1.0e-6;
   Int_t N = tmpW.nPt(); //Y.size

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
	   //outW.addPnt(tmpVal);
   }

}

void filterWave(wavefft &inWfft, Double_t *ReFFT, Double_t *ImFFT, double minCut=50e+6, double maxCut=1500e+6){//gli passi l'onda trasformata, perch� andremo ad agire sull'ampiezza, tagliando le frequenze.
	double aveMin=0.0;
	int nMin=0;
	double aveMax=0.0;
	int nMax=0;
	for (int ip=0; ip<inWfft.nPt/2; ++ip) {
		if(inWfft.omega[ip]>minCut*0.9 && inWfft.omega[ip]<minCut*1.1) {
			aveMin+=inWfft.ampl[ip];
			++nMin; //quante volte entri nell'if.
		}
		if(inWfft.omega[ip]>maxCut*0.9 && inWfft.omega[ip]<maxCut*1.1) {
			aveMax+=inWfft.ampl[ip];
			++nMax;
		}
	}
	aveMin/=((double)nMin);//media.
	aveMax/=((double)nMax);
	std::cout<<"AveMin "<<aveMin<<" aveMax "<<aveMax<<std::endl;
	for (int ip=0; ip<inWfft.nPt/2; ++ip) {
		if (inWfft.omega[ip]>minCut && inWfft.omega[ip]<maxCut) {
			double newAmp = (aveMax-aveMin)/(maxCut-minCut)*(inWfft.omega[ip]-minCut) + aveMin;//retta di taglio.
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

void filterWaveSM(wavefft &inWfft, Double_t *ReFFT, Double_t *ImFFT, int smtM=5, int smtK=2){ //gli passi sempre l'onda trasformata ma il "taglio" stavolta � lo smooth
	std::vector<float> fltAmpl = smoothSG(inWfft.ampl,inWfft.nPt/2,smtM,smtK); //la nuova ampiezza � quella delle forme d'onda con smooth
  //std::vector<float> fltPhi = smoothSG(inWfft.phi,inWfft.nPt/2,5,2);
	int n=(smtM-1)/2;
	int nPt=inWfft.nPt/2;
	for (int ip=0; ip<nPt; ++ip) {
	double newAmp=0.0;
//	double newPhi=0.0;
//			std::cout<<"NewAmp "<<newAmp<<std::endl;
		if (ip<n){ 
    newAmp=fltAmpl[0];  
    //newPhi=fltPhi[0];
    }
		else if (ip>=nPt-n){ 
    newAmp=fltAmpl.back(); 
   // newPhi=fltPhi.back();
    }
		else { newAmp=fltAmpl[ip-n]; 
    //newPhi=fltPhi[ip-n];
    }
		ReFFT[ip]=newAmp*TMath::Cos(/*newPhi*/inWfft.phi[ip]);
		ImFFT[ip]=newAmp*TMath::Sin(/*newPhi*/inWfft.phi[ip]); 
	}
	for (int ip=inWfft.nPt-1; ip>=inWfft.nPt/2; --ip) {
		ReFFT[ip]=0.0;//newAmp*TMath::Cos(inWfft.phi[ip]);
		ImFFT[ip]=0.0;//newAmp*TMath::Sin(inWfft.phi[ip]);
	}

}
////////////////// trasformata di fourier semplice, dove togliamo i primi bin. Serve a ripristinare la baseline.
void filterWaveBsl(wavefft &inWfft, Double_t *ReFFT, Double_t *ImFFT){
	for (int ip=0; ip<inWfft.nPt/2; ++ip) {
		if(ip>0){
			ReFFT[ip]=inWfft.ampl[ip]*TMath::Cos(inWfft.phi[ip]);
			ImFFT[ip]=inWfft.ampl[ip]*TMath::Sin(inWfft.phi[ip]);
		}
		else{
			ReFFT[ip]=0;
			ImFFT[ip]=0;     
		}
	}
	for (int ip=inWfft.nPt-1; ip>=inWfft.nPt/2; --ip) {
		ReFFT[ip]=0.0;//inWfft.real[ip];
		ImFFT[ip]=0.0;//inWfft.img[ip];
		//			ReFFT[ip]=inWfft.ampl[ip]*TMath::Cos(inWfft.phi[ip]);
		//			ImFFT[ip]=inWfft.ampl[ip]*TMath::Sin(inWfft.phi[ip]);
	}

}
//FILTRO CR filtro passa alto
 void filterWaveCR(wavefft &inWfft, Double_t *ReFFT, Double_t *ImFFT, double tau=1/50e+6){
//	wavefft FILTER;
  for (int ip=0; ip<inWfft.nPt/2; ++ip) {
      double filterAmp=inWfft.omega[ip]*tau/TMath::Sqrt(1.0+(inWfft.omega[ip]*tau)*(inWfft.omega[ip]*tau));
      double filterPhi= 0;//TMath::ATan2(1.0,(inWfft.omega[ip]*tau));	
      //cout<< " omega " << inWfft.omega[ip] << " ampiezza flt " << filterAmp << " filterPhi " << filterPhi << endl;
    	ReFFT[ip]=inWfft.ampl[ip]*filterAmp*TMath::Cos(inWfft.phi[ip]+filterPhi);
			ImFFT[ip]=inWfft.ampl[ip]*filterAmp*TMath::Sin(inWfft.phi[ip]+filterPhi);		 
      
      	}       
	for (int ip=inWfft.nPt-1; ip>=inWfft.nPt/2; --ip) {
			ReFFT[ip]=0.0;//newAmp*TMath::Cos(inWfft.phi[ip]);
			ImFFT[ip]=0.0;//newAmp*TMath::Sin(inWfft.phi[ip]);	
      //FILTER.real=ReFFT;
      //FILTER.img=ImFFT;
	}  
    
    //return FILTER;
}
//FILTRO RC filtro passa basso
void filterWaveRC(wavefft &inWfft, Double_t *ReFFT, Double_t *ImFFT, double tau=50e+6){
	for (int ip=0; ip<inWfft.nPt/2; ++ip) {
      double filterAmp=1.0/TMath::Sqrt(1.0+(inWfft.omega[ip]*tau)*(inWfft.omega[ip]*tau));
      double filterPhi= 0;//TMath::ATan2((inWfft.omega[ip]*tau),1.0);	
    	ReFFT[ip]=inWfft.ampl[ip]*filterAmp*TMath::Cos(inWfft.phi[ip]+filterPhi);
			ImFFT[ip]=inWfft.ampl[ip]*filterAmp*TMath::Sin(inWfft.phi[ip]+filterPhi);		
	}
	for (int ip=inWfft.nPt-1; ip>=inWfft.nPt/2; --ip) {
			ReFFT[ip]=0.0;//newAmp*TMath::Cos(inWfft.phi[ip]);
			ImFFT[ip]=0.0;//newAmp*TMath::Sin(inWfft.phi[ip]);	 
	}
}
//Notch filter
void filterNotch(wavefft &inWfft, Double_t *ReFFT, Double_t *ImFFT, double omega_0=50e+6, double width=1e+6, int type=0){
	TComplex omg0_2(omega_0,0);
  TComplex ratio(omega_0/width,0);
  ratio=omg0_2/ratio;
  omg0_2*=omg0_2; 
  TComplex s, den, num, outFlt;
  
  for (int ip=0; ip<inWfft.nPt/2; ++ip) {
      s(0,inWfft.omega[ip]);
      den = s*s+ratio*s+omg0_2;
      if(type==0){ num= omg0_2; }        //lowpass
      else if(type==1){ num= ratio*s;}   //bandpass
      else if(type==2){ num= s*s+omg0_2;}//Notch,bandstop
      else if(type==3){ num= s*s;}       //Highpass
      outFlt(inWfft.real[ip],inWfft.img[ip]);
      outFlt=(num/den)*outFlt;
    	ReFFT[ip]=outFlt.Re();
			ImFFT[ip]=outFlt.Im();		
	}
	for (int ip=inWfft.nPt-1; ip>=inWfft.nPt/2; --ip) {
			ReFFT[ip]=0.0;//newAmp*TMath::Cos(inWfft.phi[ip]);
			ImFFT[ip]=0.0;//newAmp*TMath::Sin(inWfft.phi[ip]);	
	}
}
//filtro cancellazione della coda.
void filterpole_zero(wavefft &inWfft, Double_t *ReFFT, Double_t *ImFFT, Double_t tau=50e+6, int A=1,int n=1){
	TComplex s,den,outFlt;
	TComplex t(tau,0);
	TComplex Gain(A,0);
	den = TComplex::Power(1.0+s*t,n+1);
	for (int ip=0; ip<inWfft.nPt/2; ++ip) {
		s(0,inWfft.omega[ip]);
		outFlt(inWfft.real[ip],inWfft.img[ip]);
		outFlt=Gain/den*outFlt;
		ReFFT[ip]=outFlt.Re();
		ImFFT[ip]=outFlt.Im();
	}
	for (int ip=inWfft.nPt-1; ip>=inWfft.nPt/2; --ip) {
		ReFFT[ip]=0.0;//newAmp*TMath::Cos(inWfft.phi[ip]);
		ImFFT[ip]=0.0;//newAmp*TMath::Sin(inWfft.phi[ip]);
	}
}

#endif
