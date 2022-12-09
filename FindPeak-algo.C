Int_t FindPeaks(Int_t jentry, Int_t skipFstBin, Int_t channel, Int_t npt, Float_t *amplitude, Float_t rms, Float_t *fderiv, Float_t *sderiv, Int_t *pkPos, Float_t *pkHgt, Float_t timeRes, Float_t N_1, Float_t N_2, Float_t N_3, Float_t N_4, Float_t bslnTimeInterval, Int_t isChannel_1cm, Int_t isChannel_2cm,Int_t isChannel_1p5cm) {
  // In read_data.C you find this line:
  // FindPeaks(((wave)Waves_signal_1[channel]).nPtInR(),&((wave)Waves_signal_1[channel]).Y[skipFstBin],&((wave)Waves_signal_1[channel]).rms,&((wave)Waves_signal_1[channel]).deriv[skipFstBin],&((wave)Waves_signal_1[channel]).sderiv[skipFstBin],pkPos,pkHgt);
  
  
  //Float_t sigd1 = (float)(rms/(sqrt(2)*timeRes));
  //Float_t sigd2 = (float)(rms/(2*timeRes*timeRes));
  Float_t sigd1 = (float)(rms/(sqrt(2)*timeRes));
  Float_t sigd2 = (float)(rms/(2*timeRes*timeRes));
  
  Int_t maxNPks = 200; //With maxNPks = 130 *** Break *** segmentation violation
  Int_t nPks=0;
  Float_t max_amplitude= -1.0;
  Int_t index_max = -1;
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
    
    //cout << " i-th event: " << jentry <<" channel "<<channel << " ip "<< ip <<  " Amplitude ip "<< amplitude[ip] << " First Derivative ip "<< fderiv[ip] << " Second Derivative ip "<< sderiv[ip]<<"\n";
	
    if (amplitude[ip]>(float)(N_1*rms) && (amplitude[ip] - (float) (amplitude[ip-1]+amplitude[ip+1])/2 > (float) N_2*rms) && ((abs(fderiv[ip])< (float) (2.0 * N_3 * sigd1)) || (fderiv[ip-1] >  (float) N_3 * sigd1 || fderiv[ip+1] < (float) (-1.0 *(float) N_3 * sigd1))) && sderiv[ip] < (float)( -1.0 * N_4 * sigd2 )
	
	) {
      
      //if (isChannel_1cm && ((ip+skipFstBin)* timeRes <= (300+300)) && ((ip+skipFstBin)* timeRes >= (bslnTimeInterval+skipFstBin*timeRes))){
        if (isChannel_1cm){
	if(nPks==0){
	  pkHgt[nPks]= amplitude[ip];
	  pkPos[nPks]= ip;
	  nPks=nPks+1;
	}
	else{
	  if(ip - pkPos[nPks-1] > 1){ // Separated Electrons
      if (max_amplitude > 0. && index_max > 0){
        pkHgt[nPks-1] = max_amplitude;
	      pkPos[nPks-1] = index_max;
      }
      max_amplitude = -1.0; 
      index_max = -1; 
	    pkHgt[nPks]= amplitude[ip];
	    pkPos[nPks]= ip;
	    nPks=nPks+1;
	  }
	  else{ // consecutive electrons
      if(amplitude[ip]> max_amplitude){ 
        if(amplitude[ip-1] > amplitude[ip]){
          max_amplitude = amplitude[ip-1]; 
          index_max = ip-1; 
        }
        else{
        max_amplitude = amplitude[ip]; 
        index_max = ip; 
        }
      }
      pkHgt[nPks-1] = amplitude[ip];
	    pkPos[nPks-1] = ip;
	  }
	  
	}
	
	//cout << "i-th event: " << jentry << " Peak in bin "<< pkPos[nPks-1] << " Peak in ns "<< pkPos[nPks-1]*timeRes<<endl;
      }
      
      //if (isChannel_2cm && ((ip+skipFstBin)* timeRes >= (bslnTimeInterval+skipFstBin*timeRes)) && ((ip+skipFstBin)*timeRes <= 800)){
       if (isChannel_2cm || isChannel_1p5cm){
	if(nPks==0){
	  pkHgt[nPks]= amplitude[ip];
	  pkPos[nPks]= ip;
	  nPks=nPks+1;
	}
	else{
	  if(ip - pkPos[nPks-1] > 1){ // if ip is No consecutive electrons
      if (max_amplitude > 0. && index_max > 0){ // If before there was a group of electrons
        pkHgt[nPks-1] = max_amplitude; // put in the previous Peak the maximum of the group of electrons
	      pkPos[nPks-1] = index_max;
      }
      max_amplitude = -1; // initialization to start an eventual new group of electrons
      index_max = -1; 
	    pkHgt[nPks]= amplitude[ip]; // put in the current Peak the new peak found
	    pkPos[nPks]= ip;
	    nPks=nPks+1;
      
	  }
	  else{ // If consecutive electrons
      if(amplitude[ip]> max_amplitude){  // If it is the maximum 
        if(amplitude[ip-1] > amplitude[ip]){ // Check if the first electron of the group had a higher amplitude
          max_amplitude = amplitude[ip-1]; 
          index_max = ip-1; // added
        }
        else{
        max_amplitude = amplitude[ip]; 
        index_max = ip; 
        }
      }
      pkHgt[nPks-1] = amplitude[ip]; // store in the same position the current consecutive electron
	    pkPos[nPks-1] = ip;
	  }
	  
	}
	//cout << "i-th event: " << jentry << " Peak in time [ns]: " <<(ip+skipFstBin)*0.833 <<" Peak in Bin: " <<(ip+skipFstBin) <<" Peaks: "<< nPks << " Amplitude of the Electron Peak "<< pkHgt[nPks] <<"\n";
	
	//cout << "i-th event: " << jentry << " Peak in bin "<< pkPos[nPks-1] << " Peak in ns "<< pkPos[nPks-1]*timeRes<<endl;
      }
      
    }
    
    if(nPks==maxNPks) {
      // cout <<"\n";
      break;}
    
    
  }
  //}
  
  
  
  
  return nPks;
  
}