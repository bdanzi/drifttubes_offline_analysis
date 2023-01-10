Int_t ClusterizationFindPeaks(Float_t cut,Int_t *nElectrons_per_cluster,Int_t jentry, Int_t skipFstBin, Int_t channel, Int_t npt, Float_t rms, Int_t *pkPos_clust, Float_t *pkHgt_clust, Int_t *pkPos, Float_t *pkHgt, Int_t NPeak, Float_t timeRes,Int_t isChannel_1cm, Int_t isChannel_2cm, Int_t isChannel_1p5cm, Float_t scale_cut) {
  
  Int_t nPks_constant = NPeak;
  Int_t NPeak_clust = 0;
  Float_t difference_ns = 0.;
  Float_t max_amplitude= -1.0;
  Int_t index_max = -1;
  
  for(int i = 0; i < nPks_constant; ++i){
    pkPos_clust[i] = 0;
    pkHgt_clust[i] = 0.;
    nElectrons_per_cluster[i] = 1;
  }
  
  for (int ip = 0; ip < nPks_constant; ++ip) {
  //   if(isChannel_1cm){
  //   cut = 0.24*sqrt(pkPos[ip]*timeRes);
  // }
  
    cut = scale_cut*sqrt(pkPos[ip]*timeRes);
  
  
    if(NPeak_clust == 0){
      pkHgt_clust[NPeak_clust]= pkHgt[ip];
      pkPos_clust[NPeak_clust]= pkPos[ip];
      NPeak_clust = NPeak_clust + 1;
    }
    else{
      difference_ns = (pkPos[ip] - pkPos_clust[NPeak_clust-1])*timeRes;
      if(difference_ns > cut){ // Electrons in Different Clusters
	
  if (max_amplitude > 0. && index_max > 0){
	  pkHgt_clust[NPeak_clust-1] = max_amplitude;
	  pkPos_clust[NPeak_clust-1] = index_max;
	}
	max_amplitude = -1.0; 
	index_max = -1; 
	pkHgt_clust[NPeak_clust]= pkHgt[ip];
	pkPos_clust[NPeak_clust]= pkPos[ip];
	NPeak_clust = NPeak_clust + 1;
      }
      else{ // Electrons in the same Cluster
	if(pkHgt[ip]> pkHgt_clust[NPeak_clust-1]){
	    max_amplitude = pkHgt[ip]; 
	    index_max = pkPos[ip]; 
	  }
	  else{
	    max_amplitude = pkHgt_clust[NPeak_clust-1]; 
	    index_max = pkPos_clust[NPeak_clust-1]; 
	  }
	
	nElectrons_per_cluster[NPeak_clust-1] = nElectrons_per_cluster[NPeak_clust-1] + 1;
	pkHgt_clust[NPeak_clust-1] = max_amplitude;
	pkPos_clust[NPeak_clust-1] = index_max;
      }
      
    }
  }
  //for(int m = 0; m < NPeak_clust; m++){
  //	cout << "i-th event: " << jentry << " m Cluster index: "<< m << " Cluster Peak in time [ns]: "<<(pkPos_clust[m]+skipFstBin)*0.833 <<" channel: "<< channel <<" Cluster Peak in Bin: " <<(pkPos_clust[m]+skipFstBin) <<" Amplitude Cluster Peaks: "<< pkHgt_clust[m] <<" Peaks: "<< NPeak_clust << " Number of Electrons per Cluster "<< nElectrons_per_cluster[m] <<"\n";
  //}
  //cout << "i-th event: " << jentry << endl;
  return NPeak_clust;
  
  
}
