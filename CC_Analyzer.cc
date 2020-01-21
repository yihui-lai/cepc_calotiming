
#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TBrowser.h"
#include "TH2.h"
#include "TRandom.h"

void CC_Analyzer() {

  TFile *f = new TFile("testtest.root");
  TTree *t1 = (TTree*)f->Get("tree");
  float depositedEnergyTotal;
  t1->SetBranchAddress("depositedEnergyTotal",&depositedEnergyTotal);
  Int_t nentries = (Int_t)t1->GetEntries();
  for(Int_t i=0;i<nentries; i++) {
    t1->GetEntry(i);
    std::cout<<"total energy deposited is "<<depositedEnergyTotal<<std::endl;
  }

}
