
#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TBrowser.h"
#include "TH2.h"
#include "TRandom.h"


void CC_Ana(const char* inputfilename,const char* outputfilename) {

  TH1F *hTotalE = new TH1F("hTotalE","energy total world (GeV)",100,0.,100.);

  TFile *f = new TFile(inputfilename);
  TTree *t1 = (TTree*)f->Get("tree");
  float depositedEnergyTotal,depositedEnergyWorld;
  float depositedEnergyTiming_f,depositedEnergyTiming_r;
  float depositedEnergyECAL_f,depositedEnergyECAL_r;
  float depositedEnergyHCALAct,depositedEnergyHCALPas;
  float depositedEnergyEscapeWorld;

  t1->SetBranchAddress("depositedEnergyTotal",&depositedEnergyTotal);
  t1->SetBranchAddress("depositedEnergyWorld",&depositedEnergyWorld);
  t1->SetBranchAddress("depositedEnergyTiming_f",&depositedEnergyTiming_f);
  t1->SetBranchAddress("depositedEnergyTiming_r",&depositedEnergyTiming_r);
  t1->SetBranchAddress("depositedEnergyECAL_f",&depositedEnergyECAL_f);
  t1->SetBranchAddress("depositedEnergyECAL_r",&depositedEnergyECAL_r);
  t1->SetBranchAddress("depositedEnergyHCALAct",&depositedEnergyHCALAct);
  t1->SetBranchAddress("depositedEnergyHCALPas",&depositedEnergyHCALPas);
  t1->SetBranchAddress("depositedEnergyEscapeWorld",&depositedEnergyEscapeWorld);



  Int_t nentries = (Int_t)t1->GetEntries();
  for(Int_t i=0;i<nentries; i++) {
    t1->GetEntry(i);
    std::cout<<endl<<"event number "<<i<<std::endl;
    std::cout<<"total energy deposited is "<<depositedEnergyTotal<<std::endl;
    std::cout<<"world energy deposited is "<<depositedEnergyWorld<<std::endl;
    std::cout<<"timing front energy deposited is "<<depositedEnergyTiming_f<<std::endl;
    std::cout<<"timing rear energy deposited is "<<depositedEnergyTiming_r<<std::endl;
    std::cout<<"ECAL front energy deposited is "<<depositedEnergyECAL_f<<std::endl;
    std::cout<<"ECAL rear energy deposited is "<<depositedEnergyECAL_r<<std::endl;
    std::cout<<"HCAL act  energy deposited is "<<depositedEnergyHCALAct<<std::endl;
    std::cout<<"HCAL Pas  energy deposited is "<<depositedEnergyHCALPas<<std::endl;
    std::cout<<"escape energy deposited is "<<depositedEnergyEscapeWorld<<std::endl;
    float eee=depositedEnergyTiming_f+depositedEnergyTiming_r+depositedEnergyECAL_f+depositedEnergyECAL_r+depositedEnergyHCALAct+depositedEnergyHCALPas+depositedEnergyEscapeWorld;
    float fff=depositedEnergyTotal+depositedEnergyEscapeWorld;
    std::cout<<" sum in detectors is "<<eee<<std::endl;
    std::cout<<" deposited plus escaped is "<<fff<<std::endl;
    double ggg=20.-fff;
    std::cout<<" mystery is "<<ggg<<std::endl;
    hTotalE->Fill(depositedEnergyTotal);
  }

  f->Close();

  TFile * out = new TFile(outputfilename,"RECREATE");
  hTotalE->Write();
  out->Close();

}


void CC_Analyzer(bool debug) {
  if(debug) {  
    std::cout<<"running on temp.root"<<std::endl;
    CC_Ana("temp.root","temp_hists.root");
  }
  else {
    std::cout<<"running on normal files"<<std::endl;
  CC_Ana("pions_1GeV.root","pions_1GeV_hists.root");
  CC_Ana("pions_5GeV.root","pions_5GeV_hists.root");
  CC_Ana("pions_10GeV.root","pions_10GeV_hists.root");
  CC_Ana("pions_20GeV.root","pions_20GeV_hists.root");
  CC_Ana("pions_50GeV.root","pions_50GeV_hists.root");
  CC_Ana("pions_100GeV.root","pions_100GeV_hists.root");
  }
}
