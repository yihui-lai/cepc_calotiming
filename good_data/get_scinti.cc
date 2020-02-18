#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TBrowser.h"
#include "TH2.h"
#include "TRandom.h"

void CC_Ana_e(const char* inputfilename,const char* name,double energy, double* aamean, double* aarms) {
    
    TFile *f = new TFile(inputfilename);
    TTree *t1 = (TTree*)f->Get("tree");
    
    auto Canvas= new TCanvas("Canvas","Canvas",200,10,700,500);
    TH1F* nhist = new TH1F("nhist","nhist",300,0.95,1);
    
    t1->Draw(Form("(depositedEnergyTiming_f-depositedIonEnergyTiming_f + depositedEnergyTiming_r -depositedIonEnergyTiming_r + depositedEnergyECAL_f - depositedIonEnergyECAL_f + depositedEnergyECAL_r -depositedIonEnergyECAL_r) / %g>>nhist",energy));
    
    //t1->Draw(Form("( depositedEnergyECAL_f - depositedIonEnergyECAL_f + depositedEnergyECAL_r -depositedIonEnergyECAL_r) / %g>>nhist",energy));
    
    int imax = nhist->GetMaximumBin();
    double amax = nhist->GetBinCenter(imax);
    std::cout<<"hist max at "<<amax<<std::endl;
    double arms = nhist->GetRMS();
    std::cout<<"hist rms is "<<arms<<std::endl;
    TF1 *f1 = new TF1("f1","gaus",amax-2*arms,amax+2*arms);
    nhist->Fit("f1","R");
    
    TF1 *fit=nhist->GetFunction("f1");
    double p0= f1->GetParameter(0);
    double p1= f1->GetParameter(1);
    double p2= f1->GetParameter(2);
    std::cout<<"fit parameters are "<<p0<<" "<<p1<<" "<<p2<<std::endl;
    
    *aamean=p1;
    *aarms=p2;
    
    nhist->GetXaxis()->SetTitle("Energy/Etrue");
    nhist->GetYaxis()->SetTitle("Events");
    
    //cout<<name<<endl;
    
    Canvas->SaveAs( name );
}

void CC_Ana_p(const char* inputfilename,const char* name,double energy, double* aamean, double* aarms) {
    
    TFile *f = new TFile(inputfilename);
    TTree *t1 = (TTree*)f->Get("tree");
    
    auto Canvas= new TCanvas("Canvas","Canvas",200,10,700,500);
    TH1F* nhist = new TH1F("nhist","nhist",100,0.5,1.1);
    
    t1->Draw(Form("(depositedEnergyTiming_f-depositedIonEnergyTiming_f + depositedEnergyTiming_r -depositedIonEnergyTiming_r + depositedEnergyECAL_f - depositedIonEnergyECAL_f + depositedEnergyECAL_r -depositedIonEnergyECAL_r) / %g>>nhist",energy));
    
    //t1->Draw(Form("(depositedEnergyECAL_f - depositedIonEnergyECAL_f + depositedEnergyECAL_r -depositedIonEnergyECAL_r) / %g>>nhist",energy));
    
    int imax = nhist->GetMaximumBin();
    double amax = nhist->GetBinCenter(imax);
    std::cout<<"hist max at "<<amax<<std::endl;
    double arms = nhist->GetRMS();
    std::cout<<"hist rms is "<<arms<<std::endl;
    TF1 *f1 = new TF1("f1","gaus",amax-2*arms,amax+2*arms);
    nhist->Fit("f1","R");
    
    TF1 *fit=nhist->GetFunction("f1");
    double p0= f1->GetParameter(0);
    double p1= f1->GetParameter(1);
    double p2= f1->GetParameter(2);
    std::cout<<"fit parameters are "<<p0<<" "<<p1<<" "<<p2<<std::endl;
    
    *aamean=p1;
    *aarms=p2;
    
    nhist->GetXaxis()->SetTitle("Energy/Etrue");
    nhist->GetYaxis()->SetTitle("Events");
    
    //cout<<name<<endl;
    
    Canvas->SaveAs( name );
}

void get_scinti() {
    
    std::cout<<"running on normal files"<<std::endl;
    
    int npoints = 7;
    int nhst=2;
    const char* filenames_p[npoints];
    filenames_p[0]="pion_1GeV.root";
    filenames_p[1]="pion_2GeV.root";
    filenames_p[2]="pion_5GeV.root";
    filenames_p[3]="pion_10GeV.root";
    filenames_p[4]="pion_20GeV.root";
    filenames_p[5]="pion_50GeV_good.root";
    filenames_p[6]="pion_100GeV_good.root";
    const char* filenames_e[npoints];
    filenames_e[0]="electron_1GeV.root";
    filenames_e[1]="electron_2GeV.root";
    filenames_e[2]="electron_5GeV.root";
    filenames_e[3]="electron_10GeV.root";
    filenames_e[4]="electron_20GeV.root";
    filenames_e[5]="electron_50GeV_good.root";
    filenames_e[6]="electron_100GeV_good.root";
    
    const char* outname_p[npoints];
    outname_p[0]="s_pion1GeV";
    outname_p[1]="s_pion2GeV";
    outname_p[2]="s_pion5GeV";
    outname_p[3]="s_pion10GeV";
    outname_p[4]="s_pion20GeV";
    outname_p[5]="s_pion50GeV";
    outname_p[6]="s_pion100GeV";
    const char* outname_e[npoints];
    outname_e[0]="s_electron1GeV";
    outname_e[1]="s_electron2GeV";
    outname_e[2]="s_electron5GeV";
    outname_e[3]="s_electron10GeV";
    outname_e[4]="s_electron20GeV";
    outname_e[5]="s_electron50GeV";
    outname_e[6]="s_electron100GeV";
    
    
    double rrmean[npoints][nhst];
    double rrrms[npoints][nhst];
    double rrres[npoints][nhst];
    
    double arrres[npoints],aatruemean[npoints];
    aatruemean[0]=1;
    aatruemean[1]=2;
    aatruemean[2]=5;
    aatruemean[3]=10;
    aatruemean[4]=20;
    aatruemean[5]=50;
    aatruemean[6]=100;
    
    
    double mean,rms;
    double cali_r[npoints];
    
    for(int k=0;k<nhst;k++) {
        for(int j=0;j<npoints;j++){
            //std::cout<<"fitting "<<hnam[k]<<" at energy "<<aatruemean[j]<<std::endl;
            if (k==0) CC_Ana_p(filenames_p[j],outname_p[j], aatruemean[j] ,&mean,&rms);
            if (k==1){
                CC_Ana_e(filenames_e[j],outname_e[j], aatruemean[j] ,&mean,&rms);
                cali_r[j] = 1/mean;
            }
            rrmean[j][k]=mean;
            rrrms[j][k]=rms;
            rrres[j][k]=0;
            if(mean!=0) rrres[j][k]=rrrms[j][k]/mean;
        }
    }
    
    
    for(int k=0;k<nhst;k++) {
        for(int j=0;j<npoints;j++){
            rrmean[j][k]=rrmean[j][k]*cali_r[j];
            rrrms[j][k]=rrrms[j][k]*cali_r[j];
            cout<<rrres[j][k]<<" ";
            rrres[j][k]=0;
            if(rrmean[j][k]!=0) rrres[j][k]=rrrms[j][k]/rrmean[j][k];
            cout<<rrres[j][k]<<" "<<endl;
        }
    }
    
    
    /*
     CC_Ana_p("pion_1GeV.root","g_pion1GeV", 1 ,&mean,&rms);
     rrres[0][0] = rms;
     rrmean[0][0] = mean;
     CC_Ana_p("pion_2GeV.root","g_pion2GeV", 2 ,&mean,&rms);
     rrres[1][0] = rms;
     rrmean[1][0] = mean;
     CC_Ana_p("pion_5GeV.root","g_pion5GeV", 5 ,&mean,&rms);
     rrres[2][0] = rms;
     rrmean[2][0] = mean;
     CC_Ana_p("pion_10GeV.root","g_pion10GeV", 10 ,&mean,&rms);
     rrres[3][0] = rms;
     rrmean[3][0] = mean;
     CC_Ana_p("pion_20GeV.root","g_pion20GeV", 20 ,&mean,&rms);
     rrres[4][0] = rms;
     rrmean[4][0] = mean;
     CC_Ana_p("pion_50GeV_good.root","g_pion50GeV",50,&mean,&rms);
     rrres[5][0] = rms;
     rrmean[5][0] = mean;
     
     CC_Ana_e("electron_1GeV.root","g_electron1GeV",1,&mean,&rms);
     rrres[0][1] = rms;
     rrmean[0][1] = mean;
     CC_Ana_e("electron_2GeV.root","g_electron2GeV",2,&mean,&rms);
     rrres[1][1] = rms;
     rrmean[1][1] = mean;
     CC_Ana_e("electron_5GeV.root","g_electron5GeV",5,&mean,&rms);
     rrres[2][1] = rms;
     rrmean[2][1] = mean;
     CC_Ana_e("electron_10GeV.root","g_electron10GeV",10,&mean,&rms);
     rrres[3][1] = rms;
     rrmean[3][1] = mean;
     CC_Ana_e("electron_20GeV.root","g_electron20GeV",20,&mean,&rms);
     rrres[4][1] = rms;
     rrmean[4][1] = mean;
     CC_Ana_e("electron_50GeV_good.root","g_electron50GeV",50,&mean,&rms);
     rrres[5][1] = rms;
     rrmean[5][1] = mean;
     */
    
    
    TF1 *f_p = new TF1("f_p","sqrt([0]*[0]/x+[1]*[1]/x/x+[2]*[2])");
    TF1 *f_e = new TF1("f_e","sqrt([0]*[0]/x+[1]*[1]/x/x+[2]*[2])");
    TF1 *calice = new TF1("calice","sqrt(0.57*0.57/x+0.016*0.016)",5,100);  // arXiv:1507.05892 taken from pg 21 but needs to be checked
    
    
    //resolution
    for (int k=0;k<npoints;k++) {arrres[k]=rrres[k][0];}
    auto g1 = new TGraph(npoints,aatruemean,arrres);
    g1->Fit("f_p");
    double p_a= f_p->GetParameter(0);
    double p_b= f_p->GetParameter(1);
    double p_c= f_p->GetParameter(2);
    
    for (int k=0;k<npoints;k++) {arrres[k]=rrres[k][1];}
    auto g2 = new TGraph(npoints,aatruemean,arrres);
    g2->Fit("f_e");
    double e_a= f_e->GetParameter(0);
    double e_b= f_e->GetParameter(1);
    double e_c= f_e->GetParameter(2);
    
    
    auto Canvas= new TCanvas("Canvas","Canvas",200,10,700,500);
    //  gStyle->SetOptStat(111111);
    //gStyle->SetOptFit();
    
    float x1_l = 0.9;
    float y1_l = 0.80;
    float dx_l = 0.60;
    float dy_l = 0.1;
    float x0_l = x1_l-dx_l;
    float y0_l = y1_l-dy_l;
    TLegend *lgd = new TLegend(x0_l-0.1,y0_l-0.1,x1_l-0.2, y1_l+0.1);
    lgd->SetBorderSize(0); lgd->SetTextSize(0.04); lgd->SetTextFont(62); lgd->SetFillColor(0);
    
    
    TH1 *frame = new TH1F("frame","",1000,0,120);
    frame->SetMinimum(0.);
    frame->SetMaximum(0.2);
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
    lgd->AddEntry(g1, Form("#pi- scinti energy reso: #frac{%0.2g}{#sqrt{x}} & #frac{%0.2g}{x} & %0.2g",p_a,p_b,p_c), "L");
    
    g2->SetMarkerColor(kGreen);
    g2->SetMarkerStyle(23);
    g2->SetMarkerSize(1.5);
    g2->SetLineColor(kGreen);
    g2->Draw("P");
    lgd->AddEntry(g2, Form("e- scinti energy reso: #frac{%0.2g}{#sqrt{x}} & #frac{%0.2g}{x} & %0.2g",e_a,e_b,e_c), "L");
    
    lgd->Draw();
    Canvas->Print("Resolution.png",".png");
    
    
    
    //mean value
    double arrmean[npoints];
    for (int k=0;k<npoints;k++) {arrmean[k]=rrmean[k][0];}
    auto g_p_m = new TGraph(npoints,aatruemean,arrmean);
    for (int k=0;k<npoints;k++) {arrmean[k]=rrmean[k][1];}
    auto g_e_m = new TGraph(npoints,aatruemean,arrmean);
    
    
    auto Canvas_m= new TCanvas("Canvas_m","Canvas_m",200,10,700,500);
    TLegend *lgd_m = new TLegend(x0_l,y0_l,x1_l, y1_l);
    lgd_m->SetBorderSize(0); lgd_m->SetTextSize(0.04); lgd_m->SetTextFont(62); lgd_m->SetFillColor(0);
    TH1 *frame_m = new TH1F("frame_m","",1000,0,120);
    frame_m->SetMinimum(0.7);
    frame_m->SetMaximum(1.1);
    frame_m->SetStats(0);
    frame_m->GetXaxis()->SetTitle("true energy (GeV)");
    frame_m->GetXaxis()->SetTickLength(0.02);
    frame_m->GetXaxis()->SetLabelSize(0.03);
    frame_m->GetYaxis()->SetTitle("E_{scinti}/E_{true} ");
    frame_m->GetYaxis()->SetLabelSize(0.03);
    frame_m->Draw("");
    
    g_p_m->SetMarkerColor(kBlue);
    g_p_m->SetMarkerStyle(21);
    g_p_m->SetMarkerSize(1.5);
    g_p_m->SetLineColor(kBlue);
    g_p_m->Draw("CP");
    lgd_m->AddEntry(g_p_m, "#pi- scinti mean energy ");
    
    g_e_m->SetMarkerColor(kGreen);
    g_e_m->SetMarkerStyle(23);
    g_e_m->SetMarkerSize(1.5);
    g_e_m->SetLineColor(kGreen);
    g_e_m->Draw("CP");
    lgd_m->AddEntry(g_e_m, "e- scinti mean energy");
    
    lgd_m->Draw();
    Canvas_m->Print("mean_Energy.png",".png");
    
    
    return;
    
    
}
