Int_t ClusterizationFindPeaks(Float_t cut,Int_t *nElectrons_per_cluster,Int_t jentry, Int_t skipFstBin, Int_t channel, Int_t npt, Float_t rms, Int_t *pkPos_clust, Float_t *pkHgt_clust, Int_t *pkPos, Float_t *pkHgt, Int_t NPeak_clust, Float_t timeRes) {
  
  
  Int_t nPks_constant = NPeak_clust;
  Int_t difference_ns = 0;
  
  bool CLUST = false;
  bool firstTime = true;
  Int_t j = 0;
  
  for(int ip=0;ip<nPks_constant;++ip){
    pkPos_clust[ip]=0;
    pkHgt_clust[ip]=0;
    nElectrons_per_cluster[ip] = 1;
  }
  
  for (int i=1; i < nPks_constant; ++i) {
    
    
    difference_ns = (pkPos[i] - pkPos[i-1])*timeRes;
    pkPos_clust[j] = pkPos[i-1];
    
    
    if(difference_ns > cut){	// separated electrons belonging to different cluster
      
      j = j + 1; // cluster index
      pkHgt_clust[j] = pkHgt[i];
      CLUST = false;
      firstTime = false;
    }
	else { // inside the cluster
	  
	  NPeak_clust = NPeak_clust - 1;
	  
	  if(CLUST){ //next association of electrons to the cluster
	    pkPos_clust[j] = pkPos[i];
	    pkHgt_clust[j] = pkHgt_clust[j]+ pkHgt[i];
	    if(difference_ns > timeRes){
	      nElectrons_per_cluster[j] = nElectrons_per_cluster[j] + 1;
		}
	  }
	  else { // first association of electrons to the cluster
	    
	    pkPos_clust[j] = pkPos[i]; // position of the cluster is associated to the last electron peak
	    pkHgt_clust[j] = pkHgt[i] + pkHgt[i-1];
	    if(difference_ns > timeRes){
	      nElectrons_per_cluster[j] = nElectrons_per_cluster[j] + 1;
	    }
	    CLUST = true;
	    
	  }
	}
    
    if((i == nPks_constant - 1) && (difference_ns > cut)){
      
      pkPos_clust[j]=pkPos[i];
      pkHgt_clust[j]=pkHgt[i];
    }
    
      
  }
  
  //for(int m = 0; m < NPeak_clust; m++){
  //	cout << "i-th event: " << jentry << " m Cluster index: "<< m << " Cluster Peak in time [ns]: "<<(pkPos_clust[m]+skipFstBin)*0.833 <<" channel: "<< channel <<" Cluster Peak in Bin: " <<(pkPos_clust[m]+skipFstBin) <<" Amplitude Cluster Peaks: "<< pkHgt_clust[m] <<" Peaks: "<< NPeak_clust << " Number of Electrons per Cluster "<< nElectrons_per_cluster[m] <<"\n";
  //}
  
  return NPeak_clust;
  
}
