Int_t FindPeaks(Int_t npt, Float_t *amplitude, Float_t rms, Float_t *fderiv, Float_t *sderiv, Int_t *pkPos, Float_t *pkHgt) {
//FindPeaks(((wave)Waves_signal_1[channel]).nPtInR(),&((wave)Waves_signal_1[channel]).Y[skipFstBin],&((wave)Waves_signal_1[channel]).rms,&((wave)Waves_signal_1[channel]).deriv[skipFstBin],&((wave)Waves_signal_1[channel]).sderiv[skipFstBin],pkPos,pkHgt);
  
  Float_t sigd1 = (float)(rms/sqrt(2));
  Float_t sigd2 = (float)(rms/2);


  Int_t  maxNPks = 100;
  Int_t nPks=0;

  for (int ip=0; ip<npt; ++ip) {
      
      if (amplitude[ip]>(float)(rms*0.66) && fderiv[ip]< sigd1 && sderiv[ip]< 0. && (amplitude[ip]-amplitude[ip-1])>rms && (amplitude[ip+1]-amplitude[ip])>rms && fderiv[ip-1]> sigd1 && fderiv[ip+1]< -sigd1 && sderiv[ip-1]<sigd2 && sderiv[ip+1]<sigd2   
			  ) {
          cout << "HELLO!"<<endl;
            pkHgt[nPks]= amplitude[ip];
            pkPos[nPks]= ip;
            cout << " i-th Bin: " <<ip <<" Peaks: "<< nPks << "\n";
            ++nPks;

      }
      if(nPks==maxNPks) {break;}
    
   
  }


  return nPks;

}
