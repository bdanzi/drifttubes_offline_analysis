Int_t FindPeaks(Int_t jentry, Int_t skipFstBin, Int_t channel, Int_t npt, Float_t *amplitude, Float_t rms, Float_t *fderiv, Float_t *sderiv, Int_t *pkPos, Float_t *pkHgt, Float_t timeRes, Float_t N_1, Float_t N_2, Float_t N_3, Float_t N_4, Float_t bslnTimeInterval, Int_t isChannel_1cm, Int_t isChannel_2cm) {
  // In read_data.C you find this line:
  // FindPeaks(((wave)Waves_signal_1[channel]).nPtInR(),&((wave)Waves_signal_1[channel]).Y[skipFstBin],&((wave)Waves_signal_1[channel]).rms,&((wave)Waves_signal_1[channel]).deriv[skipFstBin],&((wave)Waves_signal_1[channel]).sderiv[skipFstBin],pkPos,pkHgt);
  
  
  Float_t sigd1 = (float)(rms/sqrt(2));
  Float_t sigd2 = (float)(rms/2);
  
  
  Int_t maxNPks = 200; //With maxNPks = 130 *** Break *** segmentation violation
  Int_t nPks=0;
  Int_t difference;
  //cout <<"N1: " << N_1 << endl;
  //cout <<"N2: " << N_2 << endl;
  //cout <<"N3: " << N_3 << endl;
  //cout <<"N4: " << N_4 << endl;
  /* if (channel <= 10 && channel != 4 && channel != 0 && channel !=7){ //1cm tubes
     for (int ip=0; ip<(npt-424); ++ip) {
     // Condition with data presented in Workshop: if (amplitude[ip]>(float)(2.5*rms) && (TMath::Abs(fderiv[ip-1]+fderiv[ip+1])< 0. || fderiv[ip]<0.) && sderiv[ip]< -sigd2/2  && (fderiv[ip-1]> 0. || fderiv[ip+1]<0. )
     //if (amplitude[ip]>(float)(rms) && fderiv[ip]< sigd1 && sderiv[ip]< sigd2  
     
     if (amplitude[ip]>(float)(4*rms) && ( fderiv[ip]< sigd1/2) && sderiv[ip]< 0.  && ((amplitude[ip]-amplitude[ip-1])>(rms) || (amplitude[ip+1]-amplitude[ip])<(rms)) && (fderiv[ip-1]> sigd1 || fderiv[ip+1]< -sigd1 )
     //if (amplitude[ip]>(float)(rms*0.66) && fderiv[ip]< sigd1 && sderiv[ip]< 0. && (amplitude[ip]-amplitude[ip-1])>(rms) && (amplitude[ip+1]-amplitude[ip])<(rms) && fderiv[ip-1]> sigd1 && fderiv[ip+1]< sigd1 && sderiv[ip-1]<sigd2 && sderiv[ip+1]<sigd2 
     ) {
     pkHgt[nPks]= amplitude[ip];
     pkPos[nPks]= ip;
     //cout << "i-th event: " << jentry << " Peak in time [ns]: " <<(ip+skipFstBin)*0.833 <<" Peak in Bin: " <<(ip+skipFstBin) <<" Peaks: "<< nPks << " Amplitude of the Electron Peak "<< pkHgt[nPks] <<"\n";
     
     nPks=nPks+1;
     
     
     }
     if(nPks==maxNPks) {break;}
     
     
     }
     }*/
  
  // else if(channel == 0 || channel == 4 || channel == 7 || channel == 11 ){
  
  for (int ip=0; ip<(npt); ++ip) {
    // Condition with data presented in Workshop: if (amplitude[ip]>(float)(2.5*rms) && (TMath::Abs(fderiv[ip-1]+fderiv[ip+1])< 0. || fderiv[ip]<0.) && sderiv[ip]< -sigd2/2  && (fderiv[ip-1]> 0. || fderiv[ip+1]<0. )
    //if (amplitude[ip]>(float)(rms) && fderiv[ip]< sigd1 && sderiv[ip]< sigd2  
    
    //Original condition 
    //  if (amplitude[ip]>(float)(3*rms) && ( fderiv[ip]< sigd1/2) && sderiv[ip]< 0.  && ((amplitude[ip]-amplitude[ip-1])>(rms) || (amplitude[ip+1]-amplitude[ip])<(rms)) && (fderiv[ip-1]> sigd1 || fderiv[ip+1]< -sigd1 )
    
    // First attempt 28 October
    //if (amplitude[ip]>(float)(3*rms) && ( fderiv[ip]< sigd1/2) && sderiv[ip]< sigd2/2  && ((amplitude[ip]-amplitude[ip-1])>(rms) || (amplitude[ip+1]-amplitude[ip])<(rms)) && (fderiv[ip-1]> sigd1/2 || fderiv[ip+1]< -sigd1/2 )
    
    // Attempt 2 November 
    //if (amplitude[ip]>(float)(3*rms) && ( TMath::Abs(fderiv[ip])< sigd1/2) && sderiv[ip]< sigd2/2  && ((amplitude[ip]-amplitude[ip-1])>(rms) || (amplitude[ip+1]-amplitude[ip])<(rms)) && (fderiv[ip-1]> sigd1/2 || fderiv[ip+1]< sigd1/2 )
    
    // Attempt 2 November with variables 18.51 pm
    if (amplitude[ip]>(float)(N_1*rms) && (amplitude[ip] - (float) (amplitude[ip-1]+amplitude[ip+1])/2 > (float) N_2*rms) && ((fderiv[ip])< (float) N_3 * sigd1) && (fderiv[ip-1] > N_3 * sigd1 || fderiv[ip+1] < N_3 * sigd1) && sderiv[ip] < (float) N_4 *sigd2 
	
	) {
      
      if (isChannel_1cm && ((ip+skipFstBin)* timeRes <= (300+300)) && ((ip+skipFstBin)* timeRes >= (bslnTimeInterval+skipFstBin*timeRes))){
	//cout << "i-th event: " << jentry << " Peak in time [ns]"<< endl;
	if(nPks==0){
	  pkHgt[nPks]= amplitude[ip];
	  pkPos[nPks]= ip;
	  nPks=nPks+1;
	}
	else{
	  if(ip - pkPos[nPks-1] > 1){
	    pkHgt[nPks]= amplitude[ip];
	    pkPos[nPks]= ip;
	    nPks=nPks+1;
	  }
	  else{
	    pkHgt[nPks-1] = amplitude[ip];
	    pkPos[nPks-1] = ip;
	  }
	  
	}
	//cout << "i-th event: " << jentry << " Peak in time [ns]: " <<(ip+skipFstBin)*0.833 <<" Peak in Bin: " <<(ip+skipFstBin) <<" Peaks: "<< nPks << " Amplitude of the Electron Peak "<< pkHgt[nPks] <<"\n";
	
	
      }
      
      if (isChannel_2cm && ((ip+skipFstBin)* timeRes >= (bslnTimeInterval+skipFstBin*timeRes)) && ((ip+skipFstBin)*timeRes <= 800)){
	//cout << "i-th event: " << jentry << " Peak in time [ns]"<< endl;
	if(nPks==0){
	  pkHgt[nPks]= amplitude[ip];
	  pkPos[nPks]= ip;
	  nPks=nPks+1;
	}
	else{
	  if(ip - pkPos[nPks-1] > 1){
	    pkHgt[nPks]= amplitude[ip];
	    pkPos[nPks]= ip;
	    nPks=nPks+1;
	  }
	  else{
	    pkHgt[nPks-1] = amplitude[ip];
	    pkPos[nPks-1] = ip;
	  }
	  
	}
	//cout << "i-th event: " << jentry << " Peak in time [ns]: " <<(ip+skipFstBin)*0.833 <<" Peak in Bin: " <<(ip+skipFstBin) <<" Peaks: "<< nPks << " Amplitude of the Electron Peak "<< pkHgt[nPks] <<"\n";
	
	
      }
      
    }
    
    if(nPks==maxNPks) {
      // cout <<"\n";
      break;}
    
    
  }
  //}
  
  
  
  
  return nPks;
  
}
