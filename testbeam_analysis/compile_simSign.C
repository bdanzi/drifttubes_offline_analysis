#if !defined(__CINT__) || defined(__MAKECINT__)

#include <TTree.h>
#include <TFile.h>
#include <TString.h>
#include <TROOT.h>
//#include <string>
//#include <iostream>
#include <TSystem.h>
#include <TH2.h>
#include "TChain.h"
#include <stdlib.h>
#include <cstring>

#include "simulationSign.h"


#endif

using namespace std;


int main (int argc, char ** argv){

  Char_t *name;//[300];
  char fld[200];
  if (argc>=2) {
//	  if (argv[1][0]=='d') {
	  if (strcmp(argv[1],"d")==0) {
		  sprintf(fld,"/lustre/cms/store/user/taliercio/TestBeam/Drift");
	  } else {
		  sprintf(fld,"%s",argv[1]);
	  }
  }
  if (argc>=3) {
    name=argv[2];
  }

  int nEv=1;
  if (argc>=4) {
    nEv=atoi(argv[3]);
  }

  int nCh=1;
  if (argc>=5) {
    nCh=atoi(argv[4]);
  }
//    printf("%s\n",name);
//    return 0;
  DoSim(name,nEv,nCh);

}
