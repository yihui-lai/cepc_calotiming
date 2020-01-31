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

  const int npoints=5;
  const char* filenames[npoints];
  filenames[0]="pions_5GeV_hists.root"; 
  filenames[1]="pions_10GeV_hists.root"; 
  filenames[2]="pions_20GeV_hists.root"; 
  filenames[3]="pions_50GeV_hists.root"; 
  filenames[4]="pions_100GeV_hists.root"; 
  


  double aamean[npoints],aarms[npoints],rrres[npoints];
  double abc,dej;
  for(int j=0;j<npoints;j++){
    resolution(filenames[j],&abc,&dej);
    aamean[j]=abc;aarms[j]=dej;
    rrres[j]=0;
    if(aamean[j]!=0) rrres[j]=aarms[j]/aamean[j];
    std::cout<<"jjjjjj "<<j<<" "<<aamean[j]<<" "<<rrres[j]<<std::endl;
  }
  

  
  auto Canvas= new TCanvas("Canvas","Canvas",200,10,700,500);

  auto g = new TGraph(npoints,aamean,rrres);
  TF1 *f2 = new TF1("f2","sqrt([0]*[0]/x+[1]*[1]/x/x+[2]*[2])");
  g->Fit("f2");
  gStyle->SetOptFit();
  g->SetMarkerColor(kBlue);
  g->SetMarkerStyle(21);
  g->SetMarkerSize(1.5);
  g->Draw("AP");
  
  


  return;


}
