
#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TBrowser.h"
#include "TH2.h"
#include "TRandom.h"

void CC_Analyzer() {

  TFile *f = new TFile("testtest.root");
  TTree *t1 = (TTree*)f->Get("tree");
  float depositedEnergyTotal,depositedEnergyWorld;
  float depositedEnergyTiming_f,depositedEnergyTiming_r;
  float depositedEnergyECAL_f,depositedEnergyECAL_r;

  t1->SetBranchAddress("depositedEnergyTotal",&depositedEnergyTotal);
  t1->SetBranchAddress("depositedEnergyWorld",&depositedEnergyWorld);
  t1->SetBranchAddress("depositedEnergyTiming_f",&depositedEnergyTiming_f);
  t1->SetBranchAddress("depositedEnergyTiming_r",&depositedEnergyTiming_r);
  t1->SetBranchAddress("depositedEnergyECAL_f",&depositedEnergyECAL_f);
  t1->SetBranchAddress("depositedEnergyECAL_r",&depositedEnergyECAL_r);

  Int_t nentries = (Int_t)t1->GetEntries();
  for(Int_t i=0;i<nentries; i++) {
    t1->GetEntry(i);
    std::cout<<endl<<"event number "<<i<<std::endl;
    std::cout<<"total energy deposited is "<<depositedEnergyTotal<<std::endl;
    std::cout<<"world energy deposited is "<<depositedEnergyWorld<<std::endl;
    std::cout<<"timing front energy deposited is "<<depositedEnergyTiming_f<<std::endl;
    std::cout<<"timing rear energy deposited is "<<depositedEnergyTiming_r<<std::endl;
    std::cout<<"ECAL front energy deposited is "<<depositedEnergyECAL_f<<std::endl;
    std::cout<<"ECAL read energy deposited is "<<depositedEnergyECAL_r<<std::endl;
  }

}
