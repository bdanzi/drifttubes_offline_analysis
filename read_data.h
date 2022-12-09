//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Mon Oct  8 15:07:32 2018 by ROOT version 6.10/09
// from TTree data/
// found on file: /lustre/cms/store/user/taliercio/TestBeam/Drift/run_126.root
//////////////////////////////////////////////////////////

#ifndef read_data_h
#define read_data_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

#include "WaveData/inc/WaveData.h"
#include "WaveData/inc/RunHeader.h"
#include "TObject.h"
//#include <vector>
#include <string>
#include <iomanip>
#include <stdio.h>
#include <fstream>
#include <iostream>
//#include <map>

#include "dataUtils.h"

using namespace std;

class read_data {
 public :
  TTree          *fChain;   //!pointer to the analyzed TTree or TChain
  Int_t           fCurrent; //!current Tree number in a TChain
  
  // Fixed size dimensions of array or collections stored in the TTree if any.
  static constexpr Int_t kMax_chdata_x742 = 18;
  
  WaveData *wd; 
  
  // Declaration of leaf types
  read_data(TTree *tree=0);
  virtual ~read_data();
  virtual Int_t    Cut(Long64_t entry);
  virtual Int_t    GetEntry(int entry);
  virtual Long64_t LoadTree(Long64_t entry);
  virtual void     Init(TTree *tree);
  virtual void     Loop(Char_t *name, Int_t MidEv=0,Int_t eventn=-1, float _gsample = 1.2, float N_1 = 3.0,float N_2= 1.0, float N_3 = 0.5, float N_4=0.5, float bslnTimeInterval = 25, int _dim = 1024, float _scale_cut = 0.20);
  virtual Bool_t   Notify();
  virtual void     Show(Long64_t entry = -1);
  
  double _tmax;
  bool isNov2021TestBeam;
  int dim;
  int nMaxCh;
  
 private:
  void doChDiff(int evN, std::pair<int,int> &chDiff, std::map<int, bool> &isFull, diffWvCont &diffWaves, WvCont &Waves, WvCont &Smt4Waves, int SF, hDiffCont &hDiff, hDiffCont &hDiffPed, TFile *theFile);
  
};

#endif

#ifdef read_data_cxx
read_data::read_data(TTree *tree) : fChain(0) 
{
  // if parameter tree is not specified (or zero), connect the file
  // used to generate this class and read the Tree.
  if (tree == 0) {
    TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("/lustre/cms/store/user/taliercio/TestBeam/Drift/run_126.root");
    if (!f || !f->IsOpen()) {
      f = new TFile("/lustre/cms/store/user/taliercio/TestBeam/Drift/run_trial.root");
      cout << "Reading run_trial.root that does not exist" << endl;
    }
    f->GetObject("data",tree);
  }
  
  
  // dim = 1024;
  // //_tmax = 853.3e-9;//1.0e-6;//853.3e-9;//1.0e-6;//853.3e-9;//1.0e-6;
  // _tmax = (float) 1/(_gsample*1.0e+9)*dim;
  // nMaxCh = 15;
  Init(tree);
}

read_data::~read_data()
{
  if (!fChain) return;
  delete fChain->GetCurrentFile();
}

Int_t read_data::GetEntry(int entry)
{
  // Read contents of entry.
  if (!fChain) return 0;
  return fChain->GetEntry(entry);
}
Long64_t read_data::LoadTree(Long64_t entry)
{
  // Set the environment to read one entry
  if (!fChain) return -5;
  Long64_t centry = fChain->LoadTree(entry);
  if (centry < 0) return centry;
  if (fChain->GetTreeNumber() != fCurrent) {
    fCurrent = fChain->GetTreeNumber();
    Notify();
  }
  return centry;
}

void read_data::Init(TTree *tree)
{
  // The Init() function is called when the selector needs to initialize
  // a new tree or chain. Typically here the branch addresses and branch
  // pointers of the tree will be set.
  // It is normally not necessary to make changes to the generated
  // code, but the routine can be extended by the user if needed.
  // Init() will be called many times when running on PROOF
  // (once per file to be processed).
  
  // Set branch addresses and branch pointers
  if (!tree) return;
  fChain = tree;
  fCurrent = -1;
  
  wd = new WaveData();
  fChain->SetBranchAddress("x", &wd);

  Notify();
}

Bool_t read_data::Notify()
{
  // The Notify() function is called when a new file is opened. This
  // can be either for a new TTree in a TChain or when when a new TTree
  // is started when using PROOF. It is normally not necessary to make changes
  // to the generated code, but the routine can be extended by the
  // user if needed. The return value is currently not used.
  
  return kTRUE;
}

void read_data::Show(Long64_t entry)
{
  // Print contents of entry.
  // If entry is not specified, print current entry
  if (!fChain) return;
  fChain->Show(entry);
}
Int_t read_data::Cut(Long64_t entry)
{
  // This function may be called from Loop.
  // returns  1 if entry is accepted.
  // returns -1 otherwise.
  return 1;
}
#endif // #ifdef read_data_cxx
