#include <iostream>
#include <iomanip>
#include <locale>

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

#include "vector"
#include "vector"
using std::vector;

#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

TTree          *fChain;   //!pointer to the analyzed TTree or TChain               
Int_t           fCurrent; //!current Tree number in a TChain                       



void resolution(const char* inputfilename,double* aamean,double* aarms) {

  std::cout<<"file name is "<<inputfilename<<std::endl;
  TFile *f = new TFile(inputfilename);
  TH1F* nhist = static_cast<TH1F*>(f->Get("hTotalE")->Clone());
  Int_t imax = nhist->GetMaximumBin();
  Double_t amax = nhist->GetBinCenter(imax);
  std::cout<<"hist max at "<<amax<<std::endl;
  Double_t arms = nhist->GetRMS();
  std::cout<<"hist rms is "<<arms<<std::endl;
  TF1 *f1 = new TF1("f1","gaus",amax-2*arms,amax+2*arms);
  nhist->Fit("f1","R0");
  TF1 *fit=nhist->GetFunction("f1");
  Double_t p0= f1->GetParameter(0);
  Double_t p1= f1->GetParameter(1);
  Double_t p2= f1->GetParameter(2);
  std::cout<<"fit parameters are "<<p0<<" "<<p1<<" "<<p2<<std::endl;
  *aamean=p1;
  *aarms=p2;
}



void res() {

  //  gStyle->SetOptStat(111111);


  Int_t npoints=6;
  double aamean[npoints],aarms[npoints],rrres[npoints];
  double abc,dej;
  for(int j=0;j<npoints;j++){
    resolution("pions_50GeV_hists.root",&abc,&dej);
    aamean[j]=abc;aarms[j]=dej;
    rrres[j]=0;
    if(aamean[j]!=0) rrres[j]=aarms[j]/aamean[j];
    std::cout<<"jjjjjj "<<j<<" "<<aamean[j]<<" "<<rrres[j]<<std::endl;
  }
  

  
  auto Canvas= new TCanvas("Canvas","Canvas",200,10,700,500);
  Canvas->SetGrid();
  auto g = new TGraph(npoints,aamean,rrres);
  g->SetMarkerStyle(6);
  g->SetMarkerSize(6);
  g->Draw("ACP");
  
  


  return;


}
