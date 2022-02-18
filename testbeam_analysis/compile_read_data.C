#if !defined(__CINT__) || defined(__MAKECINT__)

#include "read_data.h"
#include <TTree.h>
#include <TFile.h>
#include <TString.h>
#include <TROOT.h>
//#include <string>
#include <iostream>
#include <TSystem.h>
#include <TH2.h>
#include "TChain.h"
#include <stdlib.h>
#include <cstring>
#include <fstream>


#endif

using namespace std;


int main (int argc, char ** argv){

  Char_t name[300];
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
#ifndef _OSC
    sprintf(name,"%s/run_%d.root",fld,atoi(argv[2]));
#else
    sprintf(name,"%s/run-%05d.root",fld,atoi(argv[2]));
#endif
  } else {
#ifndef _OSC
    sprintf(name,"%s/run_126.root",fld);
#else
    sprintf(name,"%s/run-00000.root",fld);
#endif
  }

  int MidEv=0;
  if (argc>=4) {
	  MidEv=atoi(argv[3]);
  }

  int nEv=-1;
  if (argc>=5) {
    nEv=atoi(argv[4]);
  }


  bool evalCut=false;
  if (argc>=6) {
    evalCut=atoi(argv[5]);
  }

  TString fOutName="Full_Wave.txt";
  if (argc>=7) {
	  fOutName=argv[6];
  }



//    printf("%s\n",name);
//    return 0;
  TFile *file;
  TTree *tree;

  file = TFile::Open(name);

  cout << "Read file with name: " << name << endl;

  tree = (TTree*)file->Get("data");
  read_data test(tree);
  test.Loop(name,MidEv,nEv,evalCut,fOutName);

}
