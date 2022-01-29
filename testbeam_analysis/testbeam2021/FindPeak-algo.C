//______________________________________________________________________________
Int_t FindPeaks(Int_t npt, Float_t *amplitude, Float_t sig, Int_t *pkPos, Float_t *pkHgt) {

  Int_t nrise=12;
  Int_t checkUpTo=6;
  
  if(checkUpTo>nrise) { checkUpTo=nrise; }

//  Float_t sig=0.5e-3;
  Float_t *sigd = new Float_t[nrise];
  Float_t *sigd_2 = new Float_t[nrise];
  
/*  TH1F tmpWdist("tmpWdist","",1000,-0.5,0.5);
  for (int i=100; i<1000; ++i) { tmpWdist.Fill(amplitude[i]); }
  sig = tmpWdist.GetRMS();
  Float_t mean = tmpWdist.GetMean();
  
  cout<<"Signal Off "<<mean<<" noise "<<sig<<endl;
*/
  Float_t nSigCut=4.0;//5.0; //stava a 4
  for (int ir=checkUpTo; ir<=nrise; ++ir) {
    int irId = ir-1;
    sigd[irId]=2.4495*sig/((float)(2*ir+1));
    sigd_2[irId]=1.414*sigd[irId];//%0
    sigd_2[irId]*=nSigCut;
    sigd[irId]*=nSigCut;
  }
  sig*=1.414;
  float sig1=3.0*sig;
  sig*=nSigCut;

  Float_t *wave = new Float_t[npt];
  Float_t **riseMeas = new Float_t*[npt];// (r,nrise);
  Float_t **deri = new Float_t*[npt];// zeros(r,nrise);
  for (int ip=0; ip<npt; ++ip) {
    //wave[ip]=-1.0*(amplitude[ip]-mean);
	wave[ip]=amplitude[ip];
    riseMeas[ip]=new Float_t [nrise];
    deri[ip]=new Float_t [nrise];
    for (int ir=0; ir<nrise; ++ir) {
      riseMeas[ip][ir]=0.0;
      deri[ip][ir]=0.0;
    }
  }

  Int_t  maxNPks = 100;
//  Int_t pkPos[maxNPks];
//  Float_t pkHgt[maxNPks];
  Int_t nPks=0;
  Int_t last=-nrise;
  for (int i=nrise+1; i<npt-1; ++i) {
    for (int ir=nrise; ir>=checkUpTo; --ir) {
      int irId = ir-1;
      deri[i][irId] = (2.0*wave[i]-wave[i-ir]-wave[i-ir-1])/((float)(2*ir+1));
      riseMeas[i][irId] = (deri[i][irId]-deri[i-ir][irId]);
      
      if ( deri[i][irId]>=sigd[irId] && (deri[i][irId]-deri[i-ir][irId])>sigd_2[irId] && (wave[i]-wave[i-ir])>sig
    		  //&& (((wave[i]-wave[i-1])<sig && (wave[i]-wave[i+1])<sig) && ((wave[i]-wave[i-2])<sig && (wave[i]-wave[i+2])<sig))
			  ) {
//    	  bool skip=false;
//    	  for (int ick=1; ick<3; ++ick) {
//    		  if (fabs(wave[i]-wave[i-ick])<sig1) { skip=true; break; }
//    	  }
//        //%if globFlag
//    	  if (!skip) {
        if ((i-last)>ir /*%1*/) {
//            pkHgt[nPks]=-1.0*(wave[i]-mean);
            pkHgt[nPks]=/*-*/(wave[i]);
            pkPos[nPks]=i;
            ++nPks;
            last=i;
        } else {
//            pkHgt[nPks-1]=-1.0*(wave[i]-mean);
            pkHgt[nPks-1]=(wave[i]);
            pkPos[nPks-1]=i;
            last=i;
        }
//    	  }
      }
      if(nPks==maxNPks) {break;}
    }
    if(nPks==maxNPks) {break;}
  }

  delete [] sigd;
  delete [] sigd_2;
  delete [] wave;
  for (int ip=0; ip<npt; ++ip) {
    delete [] riseMeas[ip];
    delete [] deri[ip];
  }
  delete [] riseMeas;
  delete [] deri;
    
  return nPks;

}
