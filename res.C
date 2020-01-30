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

void res() {

  gStyle->SetOptStat(111111);

  TString canvName = "Fig_";
  canvName += "haha";
  int W = 800;
  int H = 600;
  TCanvas* canv = new TCanvas(canvName,canvName,50,50,W,H);
  // references for T, B, L, R                                                  
  float T = 0.08*H;
  float B = 0.12*H;
  float L = 0.12*W;
  float R = 0.04*W;

  canv->SetFillColor(0);
  canv->SetBorderMode(0);
  canv->SetFrameFillStyle(0);
  canv->SetFrameBorderMode(0);
  canv->SetLeftMargin( L/W );
  canv->SetRightMargin( R/W );
  canv->SetTopMargin( T/H );
  canv->SetBottomMargin( B/H );
  canv->SetTickx(0);
  canv->SetTicky(0);

  float x1_l = 0.9;
  float y1_l = 0.80;
  float dx_l = 0.60;
  float dy_l = 0.1;
  float x0_l = x1_l-dx_l;
  float y0_l = y1_l-dy_l;

  TLegend *lgd = new TLegend(x0_l,y0_l,x1_l, y1_l);
  lgd->SetBorderSize(0); lgd->SetTextSize(0.04); lgd->SetTextFont(62); lgd->SetFillColor(\
											 0);





  TFile *f = new TFile("pions_50GeV_hists.root");
  // 
    TH1F* nhist = static_cast<TH1F*>(f->Get("hTotalE")->Clone());
    nhist->Fit("gaus");
    nhist->Draw("");



    lgd->AddEntry(nhist, "haha", "l");


    lgd->Draw();


  canv->Update();
  canv->RedrawAxis();
  canv->GetFrame()->Draw();


  return;


}
