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





void resolution(const char* inputfilename,const char* histname,double* aamean,double* aarms,const char* imagename) {
    
    auto Canvas= new TCanvas("Canvas","Canvas",200,10,700,500);
    std::cout<<"file name is "<<inputfilename<<std::endl;
    TFile *f = new TFile(inputfilename);
    TH1F* nhist = static_cast<TH1F*>(f->Get(histname)->Clone());
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

    nhist->GetXaxis()->SetTitle("Energy/Etrue");
    nhist->GetYaxis()->SetTitle("Events");
    
    Canvas->SaveAs(imagename,".png");
}



void res() {
    
    const int npoints=6;
    const char* elec_filenames[npoints];
    const char* filenames[npoints];
    const char* e_name[npoints];
    const char* p_name[npoints];
    
    elec_filenames[0]="electron_1GeV_hists.root";
    elec_filenames[1]="electron_2GeV_hists.root";
    elec_filenames[2]="electron_5GeV_hists.root";
    elec_filenames[3]="electron_10GeV_hists.root";
    elec_filenames[4]="electron_20GeV_hists.root";
    elec_filenames[5]="electron_50GeV_hists_good.root";
    
    filenames[0]="pion_1GeV_hists.root";
    filenames[1]="pion_2GeV_hists.root";
    filenames[2]="pion_5GeV_hists.root";
    filenames[3]="pion_10GeV_hists.root";
    filenames[4]="pion_20GeV_hists.root";
    filenames[5]="pion_50GeV_hists_good.root";
    
    e_name[0]="electron_1GeV_hists";
    e_name[1]="electron_2GeV_hists";
    e_name[2]="electron_5GeV_hists";
    e_name[3]="electron_10GeV_hists";
    e_name[4]="electron_20GeV_hists";
    e_name[5]="electron_50GeV_hists";
    
    p_name[0]="pion_1GeV_hists";
    p_name[1]="pion_2GeV_hists";
    p_name[2]="pion_5GeV_hists";
    p_name[3]="pion_10GeV_hists";
    p_name[4]="pion_20GeV_hists";
    p_name[5]="pion_50GeV_hists";
    
    double aatruemean[npoints];
    aatruemean[0]=1;
    aatruemean[1]=2;
    aatruemean[2]=5;
    aatruemean[3]=10;
    aatruemean[4]=20;
    aatruemean[5]=50;
    
    const int nhst=2;
    double aaamean[npoints][nhst],aarms[npoints][nhst],rrres[npoints][nhst];
    
    vector<string> hnam(nhst);
    hnam[0]="hestE"; //electron
    hnam[1]="hestE"; //pion
    
    
    
    double abc,dej;
    for(int k=0;k<nhst;k++) {
        for(int j=0;j<npoints;j++){
            std::cout<<"fitting "<<hnam[k]<<" at energy "<<aatruemean[j]<<std::endl;
            if (k==0) resolution(elec_filenames[j],hnam[k].c_str(),&abc,&dej,e_name[j]);
            if (k==1) resolution(filenames[j],hnam[k].c_str(),&abc,&dej,p_name[j]);
            aaamean[j][k]=abc;
            aarms[j][k]=dej;
            rrres[j][k]=0;
            if(abc!=0) rrres[j][k]=aarms[j][k]/abc;
        }
    }
    
    
    
    
    TF1 *f2 = new TF1("f2","sqrt([0]*[0]/x+[1]*[1]/x/x+[2]*[2])");
    TF1 *calice = new TF1("calice","sqrt(0.57*0.57/x+0.016*0.016)",5,100);  // arXiv:1507.05892 taken from pg 21 but needs to be checked
    
    double arrres[npoints];
    for (int k=0;k<npoints;k++) {arrres[k]=rrres[k][0];}
    auto g1 = new TGraph(npoints,aatruemean,arrres);
    g1->Fit("f2");
    
    
    for (int k=0;k<npoints;k++) {arrres[k]=rrres[k][1];}
    auto g2 = new TGraph(npoints,aatruemean,arrres);
    g2->Fit("f2");
    
    
    auto Canvas= new TCanvas("Canvas","Canvas",200,10,700,500);
    //  gStyle->SetOptStat(111111);
    //gStyle->SetOptFit();
    
    float x1_l = 0.9;
    float y1_l = 0.80;
    float dx_l = 0.60;
    float dy_l = 0.1;
    float x0_l = x1_l-dx_l;
    float y0_l = y1_l-dy_l;
    TLegend *lgd = new TLegend(x0_l-0.1,y0_l,x1_l, y1_l);
    lgd->SetBorderSize(0); lgd->SetTextSize(0.04); lgd->SetTextFont(62); lgd->SetFillColor(0);
    
    
    TH1 *frame = new TH1F("frame","",1000,0,120);
    frame->SetMinimum(0.);
    frame->SetMaximum(0.5);
    frame->SetStats(0);
    frame->GetXaxis()->SetTitle("true energy (GeV)");
    frame->GetXaxis()->SetTickLength(0.02);
    frame->GetXaxis()->SetLabelSize(0.03);
    frame->GetYaxis()->SetTitle("percent resolution");
    frame->GetYaxis()->SetLabelSize(0.03);
    frame->Draw("");
    
    calice->Draw("same");
    lgd->AddEntry(calice, "calice detector resolution", "l");
    
    g1->SetMarkerColor(kBlue);
    g1->SetMarkerStyle(21);
    g1->SetMarkerSize(1.5);
    g1->SetLineColor(kBlue);
    g1->Draw("P");
    lgd->AddEntry(g1, Form("electron scintillation energy: %g/sqrt(x) & %g/x & %g",0.001,1.089,0.06), "L");
    
    g2->SetMarkerColor(kGreen);
    g2->SetMarkerStyle(23);
    g2->SetMarkerSize(1.5);
    g2->SetLineColor(kGreen);
    g2->Draw("P");
    lgd->AddEntry(g2, Form("pion scintillation energy: %g/sqrt(x) & %g/x & %g",0.696,0.023,0.121), "L");
    
    
    lgd->Draw();
    Canvas->Print("Resolution.png",".png");
    
    
    return;
    
    
}
