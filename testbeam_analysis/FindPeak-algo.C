Int_t FindPeaks(Int_t jentry, Int_t skipFstBin, Int_t npt, Float_t *amplitude, Float_t rms, Float_t *fderiv, Float_t *sderiv, Int_t *pkPos, Float_t *pkHgt) {
  // In read_data.C you find this line:
  // FindPeaks(((wave)Waves_signal_1[channel]).nPtInR(),&((wave)Waves_signal_1[channel]).Y[skipFstBin],&((wave)Waves_signal_1[channel]).rms,&((wave)Waves_signal_1[channel]).deriv[skipFstBin],&((wave)Waves_signal_1[channel]).sderiv[skipFstBin],pkPos,pkHgt);
  
  Float_t sigd1 = (float)(rms/sqrt(2));
  Float_t sigd2 = (float)(rms/2);
  
  
  Int_t maxNPks = 120; //With maxNPks = 130 *** Break *** segmentation violation
  Int_t nPks=0;
  
  for (int ip=0; ip<(npt-64); ++ip) {
    //if (amplitude[ip]>(float)(rms) && fderiv[ip]< sigd1 && sderiv[ip]< sigd2  
    if (amplitude[ip]>(float)(2.5*rms) && (TMath::Abs(fderiv[ip-1]+fderiv[ip+1])< 0. || fderiv[ip]<0.) && sderiv[ip]< -sigd2/2  && (fderiv[ip-1]> 0. || fderiv[ip+1]< 0. )
	//if (amplitude[ip]>(float)(rms*0.66) && fderiv[ip]< sigd1 && sderiv[ip]< 0. && (amplitude[ip]-amplitude[ip-1])>(rms) && (amplitude[ip+1]-amplitude[ip])<(rms) && fderiv[ip-1]> sigd1 && fderiv[ip+1]< sigd1 && sderiv[ip-1]<sigd2 && sderiv[ip+1]<sigd2 
	) {
      
      pkHgt[nPks]= amplitude[ip];
      pkPos[nPks]= ip;
      cout << "i-th event: " << jentry << " Peak in time [ns]: " <<(ip+skipFstBin)*0.833 <<" Peaks: "<< nPks << " Ntot points "<< npt <<"\n";
      ++nPks;
      
    }
    if(nPks==maxNPks) {break;}
    
    
  }
  
  
  return nPks;
  
}
