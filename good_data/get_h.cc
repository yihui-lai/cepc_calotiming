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
    
    //Canvas->SaveAs( name );
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
    
    //Canvas->SaveAs( name );
}

void CC_ceren(const char* inputfilename,const char* name, double energy, double* aamean, double* aarms) {
    
    TFile *f = new TFile(inputfilename);
    TTree *t1 = (TTree*)f->Get("tree");
    
    auto Canvas= new TCanvas("Canvas","Canvas",200,10,700,500);
    TH1F* htest = new TH1F("htest","htest",100000,0,100000);
    t1->Draw("tot_phot_cer_ECAL_f + tot_phot_cer_ECAL_r >>htest ");
    int imax = htest->GetMaximumBin();
    double amax = htest->GetBinCenter(imax);
    int imin = htest->GetMinimumBin();
    double amin = htest->GetBinCenter(imin);
    std::cout<<"hist max at "<<amax<<std::endl;
    double arms = htest->GetRMS();
    std::cout<<"hist rms is "<<arms<<std::endl;
    
    TH1F* nhist = new TH1F("nhist","nhist",100,amax-3*arms,amax+3*arms);
    t1->Draw(" tot_phot_cer_ECAL_f + tot_phot_cer_ECAL_r >>nhist ");
    
    TF1 *f1 = new TF1("f1","gaus",amax-2*arms,amax+2*arms);
    nhist->Fit("f1","R");
    
    TF1 *fit=nhist->GetFunction("f1");
    double p0= f1->GetParameter(0);
    double p1= f1->GetParameter(1);
    double p2= f1->GetParameter(2);
    std::cout<<"fit parameters are "<<p0<<" "<<p1<<" "<<p2<<std::endl;
    
    *aamean=p1;
    *aarms=p2;
    
    nhist->GetXaxis()->SetTitle("Nceren");
    nhist->GetYaxis()->SetTitle("Events");
    
    //cout<<name<<endl;
    
    //Canvas->SaveAs( name );
}

void CC_Fem(const char* inputfilename,const char* name, double energy, double* aamean, double* aarms) {
    
    TFile *f = new TFile(inputfilename);
    TTree *t1 = (TTree*)f->Get("tree");
    
    auto Canvas= new TCanvas("Canvas","Canvas",200,10,700,500);
    TH1F* htest = new TH1F("htest","htest",300,0.1,1.1);
    
    t1->Draw(Form("(depositedElecEnergyECAL_f + depositedElecEnergyECAL_r + depositedElecEnergyTiming_f + depositedElecEnergyTiming_r) / %g>>htest",energy));
    
    int imax = htest->GetMaximumBin();
    double amax = htest->GetBinCenter(imax);
    int imin = htest->GetMinimumBin();
    double amin = htest->GetBinCenter(imin);
    std::cout<<"hist max at "<<amax<<std::endl;
    double arms = htest->GetRMS();
    std::cout<<"hist rms is "<<arms<<std::endl;
    
    TH1F* nhist = new TH1F("nhist","nhist",100,amax-3*arms,amax+3*arms);
    t1->Draw(Form("(depositedElecEnergyECAL_f + depositedElecEnergyECAL_r + depositedElecEnergyTiming_f + depositedElecEnergyTiming_r) / %g>>nhist",energy));
    
    TF1 *f1 = new TF1("f1","gaus",amax-2*arms,amax+2*arms);
    nhist->Fit("f1","R");
    
    TF1 *fit=nhist->GetFunction("f1");
    double p0= f1->GetParameter(0);
    double p1= f1->GetParameter(1);
    double p2= f1->GetParameter(2);
    std::cout<<"fit parameters are "<<p0<<" "<<p1<<" "<<p2<<std::endl;
    
    *aamean=p1;
    *aarms=p2;
    
    nhist->GetXaxis()->SetTitle("E_elec/Etrue");
    nhist->GetYaxis()->SetTitle("Events");
    
    //cout<<name<<endl;
    
    //Canvas->SaveAs( name );
}

void get_h() {
    
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
    outname_p[0]="f_pion1GeV";
    outname_p[1]="f_pion2GeV";
    outname_p[2]="f_pion5GeV";
    outname_p[3]="f_pion10GeV";
    outname_p[4]="f_pion20GeV";
    outname_p[5]="f_pion50GeV";
    outname_p[6]="f_pion100GeV";
    const char* outname_e[npoints];
    outname_e[0]="f_electron1GeV";
    outname_e[1]="f_electron2GeV";
    outname_e[2]="f_electron5GeV";
    outname_e[3]="f_electron10GeV";
    outname_e[4]="f_electron20GeV";
    outname_e[5]="f_electron50GeV";
    outname_e[6]="f_electron100GeV";
    
    
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
    double cali_ce[npoints],cali_sc[npoints];
    
    double SS[npoints],CC[npoints],Fem[npoints];
    double SS_reso[npoints],CC_reso[npoints],Fem_reso[npoints];
    
    //ceren
    for(int k=0;k<nhst;k++) {
        for(int j=0;j<npoints;j++){
            if (k==0) CC_ceren(filenames_p[j],outname_p[j], aatruemean[j] ,&mean,&rms);
            if (k==1){
                CC_ceren(filenames_e[j],outname_e[j], aatruemean[j] ,&mean,&rms);
                cali_ce[j] = mean;
            }
            rrmean[j][k]=mean;
            rrrms[j][k]=rms;
            rrres[j][k]=0;
            if(mean!=0) rrres[j][k]=rrrms[j][k]/mean;
        }
    }
    
    for(int k=0;k<nhst;k++) {
        for(int j=0;j<npoints;j++){
            rrmean[j][k]=rrmean[j][k]/cali_ce[j];
            rrrms[j][k]=rrrms[j][k]/cali_ce[j];
            cout<<rrmean[j][k]<<" ";
            rrres[j][k]=0;
            if(rrmean[j][k]!=0) rrres[j][k]=rrrms[j][k]/rrmean[j][k];
            cout<<rrres[j][k]<<" "<<endl;
        }
    }
    for(int l=0;l<npoints;l++) {CC[l]=rrmean[l][0];CC_reso[l]=rrres[l][0];
        cout<<"CC[l]: "<<CC[l]<<endl;}
    
    //sinti
    for(int k=0;k<nhst;k++) {
        for(int j=0;j<npoints;j++){
            if (k==0) CC_Ana_p(filenames_p[j],outname_p[j], aatruemean[j] ,&mean,&rms);
            if (k==1){
                CC_Ana_e(filenames_e[j],outname_e[j], aatruemean[j] ,&mean,&rms);
                cali_sc[j] = mean;
            }
            rrmean[j][k]=mean;
            rrrms[j][k]=rms;
            rrres[j][k]=0;
            if(mean!=0) rrres[j][k]=rrrms[j][k]/mean;
        }
    }
    
    
    for(int k=0;k<nhst;k++) {
        for(int j=0;j<npoints;j++){
            rrmean[j][k]=rrmean[j][k]/cali_sc[j];
            rrrms[j][k]=rrrms[j][k]/cali_sc[j];
            cout<<rrres[j][k]<<" ";
            rrres[j][k]=0;
            if(rrmean[j][k]!=0) rrres[j][k]=rrrms[j][k]/rrmean[j][k];
            cout<<rrres[j][k]<<" "<<endl;
        }
    }
    for(int l=0;l<npoints;l++) {SS[l]=rrmean[l][0];SS_reso[l]=rrres[l][0];
        cout<<"SS[l]: "<<SS[l]<<endl;}
    
    
    double cali_r[npoints];
    for(int k=0;k<nhst;k++) {
        for(int j=0;j<npoints;j++){
            if (k==0) CC_Fem(filenames_p[j],outname_p[j], aatruemean[j] ,&mean,&rms);
            if (k==1){
                CC_Fem(filenames_e[j],outname_e[j], aatruemean[j] ,&mean,&rms);
                cali_r[j] = mean;
            }
            rrmean[j][k]=mean;
            rrrms[j][k]=rms;
            rrres[j][k]=0;
            if(mean!=0) rrres[j][k]=rrrms[j][k]/mean;
        }
    }
    /*
     for(int k=0;k<nhst;k++) {
     for(int j=0;j<npoints;j++){
     rrmean[j][k]=rrmean[j][k];//*cali_r[j];
     rrrms[j][k]=rrrms[j][k];//*cali_r[j];
     rrres[j][k]=0;
     if(rrmean[j][k]!=0) rrres[j][k]=rrrms[j][k]/rrmean[j][k];
     
     }
     }
     */
    for(int l=0;l<npoints;l++) {
        Fem[l]=rrmean[l][0];Fem_reso[l]=rrres[l][0];
        cout<<"Fem[l]: "<<Fem[l]<<endl;
    }
    
    
    
    //really begin something
    
    double chi[npoints];//,d_ceren[npoints],d_scinti[npoints];
    double C_to_S[npoints];
    double h_S[npoints],h_C[npoints];
    double E_total[npoints],E_reso[npoints];
    for(int l=0;l<npoints;l++) {
        h_S[l] = (SS[l]-Fem[l])/(1-Fem[l]);
        h_C[l] = (CC[l]-Fem[l])/(1-Fem[l]);
        
        chi[l] = (1-h_S[l])/(1-h_C[l]);
        C_to_S[l] = CC[l]/SS[l];
        
        E_total[l] = (SS[l]-chi[l]*CC[l])/(1-chi[l]);
        E_reso[l] = (sqrt(SS_reso[l]*SS_reso[l]+chi[l]*CC_reso[l]*chi[l]*CC_reso[l]))/(1-chi[l]);

        //d_ceren[l] = cali_ce[l];
        //d_scinti[l] = cali_sc[l];
        cout<<"h_S[l]: "<<h_S[l]<<endl;
        cout<<"h_C[l]: "<<h_C[l]<<endl;
    }
    for (int k=0;k<npoints;k++) {cout<<C_to_S[k]<<" "<<SS[k]<<endl;}
    for (int k=0;k<npoints;k++) {cout<<" total E:"<<E_total[k]<<" "<<E_reso[k]<<endl;}
    
    
    //add functions here, change the way calculating E total and E reso
    
    for(int loop_f=0;loop_f<npoints;loop_f++){
        
        TFile *f = new TFile(filenames_p[loop_f]);
        TTree *t1 = (TTree*)f->Get("tree");
        
        int tot_phot_cer_ECAL_f,tot_phot_cer_ECAL_r;
        float depositedIonEnergyTiming_f,depositedIonEnergyTiming_r;
        float depositedIonEnergyECAL_f,depositedIonEnergyECAL_r;
        float depositedEnergyTiming_f,depositedEnergyTiming_r;
        float depositedEnergyECAL_f,depositedEnergyECAL_r;
        
        t1->SetBranchAddress("depositedEnergyTiming_f",&depositedEnergyTiming_f);
        t1->SetBranchAddress("depositedEnergyTiming_r",&depositedEnergyTiming_r);
        t1->SetBranchAddress("depositedEnergyECAL_f",&depositedEnergyECAL_f);
        t1->SetBranchAddress("depositedEnergyECAL_r",&depositedEnergyECAL_r);
        t1->SetBranchAddress("depositedIonEnergyTiming_f",&depositedIonEnergyTiming_f);
        t1->SetBranchAddress("depositedIonEnergyTiming_r",&depositedIonEnergyTiming_r);    t1->SetBranchAddress("depositedIonEnergyECAL_f",&depositedIonEnergyECAL_f);
        t1->SetBranchAddress("depositedIonEnergyECAL_r",&depositedIonEnergyECAL_r);
        
        t1->SetBranchAddress("tot_phot_cer_ECAL_f",&tot_phot_cer_ECAL_f);
        t1->SetBranchAddress("tot_phot_cer_ECAL_r",&tot_phot_cer_ECAL_r);
        
        const int nentries = t1->GetEntries();
        
        double C_E[nentries],S_E[nentries];
        for(Int_t ii=0;ii<nentries; ii++) {
            t1->GetEntry(ii);
            C_E[ii]=(tot_phot_cer_ECAL_f + tot_phot_cer_ECAL_r)/cali_ce[loop_f];
            S_E[ii]=(depositedEnergyTiming_f-depositedIonEnergyTiming_f + depositedEnergyTiming_r -depositedIonEnergyTiming_r + depositedEnergyECAL_f - depositedIonEnergyECAL_f + depositedEnergyECAL_r -depositedIonEnergyECAL_r) /aatruemean[loop_f]/cali_sc[loop_f];
            //cout<<C_E[ii]<<" "<<S_E[ii]<<endl;
        }
        
        auto g_fit = new TGraph(nentries,S_E,C_E);
        
        auto Canvas= new TCanvas("Canvas","Canvas",200,10,700,500);
        
        g_fit->GetXaxis()->SetRangeUser(0,1);
        g_fit->GetYaxis()->SetRangeUser(0,1);
        g_fit->GetXaxis()->SetTitle("S/E");
        g_fit->GetXaxis()->SetTickLength(0.02);
        g_fit->GetXaxis()->SetLabelSize(0.03);
        g_fit->GetYaxis()->SetTitle("C/E");
        g_fit->GetYaxis()->SetLabelSize(0.03);

        TF1 *f1 = new TF1("f1","[0]*x+[1]",0,1.2);
        g_fit->Fit("f1","R");
        
        TF1 *fit=g_fit->GetFunction("f1");
        double p0= f1->GetParameter(0);
        double p1= f1->GetParameter(1);
        
        g_fit->SetMarkerStyle(21);
        g_fit->SetMarkerSize(1.5);
        g_fit->Draw("AP");
        TH1F* nhist = new TH1F("nhist","nhist",400,0.2,1.8);
        double xx,yy;
        for(Int_t ii=0;ii<nentries; ii++) {
            t1->GetEntry(ii);
            yy=(tot_phot_cer_ECAL_f + tot_phot_cer_ECAL_r)/cali_ce[loop_f];
            xx=(depositedEnergyTiming_f-depositedIonEnergyTiming_f + depositedEnergyTiming_r -depositedIonEnergyTiming_r + depositedEnergyECAL_f - depositedIonEnergyECAL_f + depositedEnergyECAL_r -depositedIonEnergyECAL_r) /aatruemean[loop_f]/cali_sc[loop_f];
            nhist->Fill(1+(p0*xx-yy+p1)/sqrt(1+p0*p0));
        }
        //nhist->Draw();
        
        TF1 *ff1 = new TF1("ff1","gaus",0.2,1.8);
        nhist->Fit("ff1","R");
        double pp0= ff1->GetParameter(0);
        double pp1= ff1->GetParameter(1);
        double pp2= ff1->GetParameter(2);
        E_reso[loop_f] = pp2/pp1;
    }
    
    
    
    
    
    
    
    
    //plots not important
    
     double arrmean[npoints];
     float x1_l = 0.9;
     float y1_l = 0.80;
     float dx_l = 0.60;
     float dy_l = 0.1;
     float x0_l = x1_l-dx_l;
     float y0_l = y1_l-dy_l;
     
     for (int k=0;k<npoints;k++) {arrmean[k]=SS_reso[k];}
     auto g_s = new TGraph(npoints,aatruemean,arrmean);
     for (int k=0;k<npoints;k++) {arrmean[k]=CC_reso[k];}
     auto g_c = new TGraph(npoints,aatruemean,arrmean);
     
     for (int k=0;k<npoints;k++) {arrmean[k]=E_reso[k];}
     auto g_factor = new TGraph(npoints,aatruemean,arrmean);
    
    TF1 *f_s = new TF1("f_s","sqrt([0]*[0]/x+[1]*[1]/x/x+[2]*[2])");
    TF1 *f_c = new TF1("f_c","sqrt([0]*[0]/x+[1]*[1]/x/x+[2]*[2])");
    TF1 *f_factor = new TF1("f_factor","sqrt([0]*[0]/x+[1]*[1]/x/x+[2]*[2])");
    g_s->Fit("f_s");
    double g_s_a= f_s->GetParameter(0);
    double g_s_b= f_s->GetParameter(1);
    double g_s_c= f_s->GetParameter(2);
    g_c->Fit("f_c");
    double g_c_a= f_c->GetParameter(0);
    double g_c_b= f_c->GetParameter(1);
    double g_c_c= f_c->GetParameter(2);
    g_factor->Fit("f_factor");
    double g_factor_a= f_factor->GetParameter(0);
    double g_factor_b= f_factor->GetParameter(1);
    double g_factor_c= f_factor->GetParameter(2);
    
    
     //auto g_factor = new TGraph(npoints,arrmean,SS);
     //TF1 *ff = new TF1("ff","[0]+x*[1]+x*x*[2]+x*x*x*[3]",0.54,1);
     //g_factor->Fit("ff","R");
     
     auto Canvas= new TCanvas("Canvas","Canvas",200,10,700,500);
     TLegend *lgd = new TLegend(x0_l,y0_l,x1_l, y1_l);
     lgd->SetBorderSize(0); lgd->SetTextSize(0.04); lgd->SetTextFont(62); lgd->SetFillColor(0);
     TH1 *frame = new TH1F("frame","",1000,0,120);
     frame->SetMinimum(0);
     frame->SetMaximum(1);
     frame->SetStats(0);
     frame->GetXaxis()->SetTitle("true energy (GeV)");
     frame->GetXaxis()->SetTickLength(0.02);
     frame->GetXaxis()->SetLabelSize(0.03);
     frame->GetYaxis()->SetTitle("hadronic response ");
     frame->GetYaxis()->SetLabelSize(0.03);
     frame->Draw("");
     
     g_factor->SetMarkerColor(kRed);
     g_factor->SetMarkerStyle(21);
     g_factor->SetMarkerSize(1.5);
     g_factor->SetLineColor(kRed);
     g_factor->Draw("CP");
     frame->GetXaxis()->SetTitle("E_{in}");
     frame->GetYaxis()->SetTitle("resolution");
     lgd->AddEntry(g_factor, Form("S+C #frac{%0.2g}{#sqrt{x}} & #frac{%0.2g}{x} & %0.2g",g_factor_a,g_factor_b,g_factor_c), "L");
     
     
     g_s->SetMarkerColor(kBlue);
     g_s->SetMarkerStyle(21);
     g_s->SetMarkerSize(1.5);
     g_s->SetLineColor(kBlue);
     g_s->Draw("CP");
     lgd->AddEntry(g_s, Form("Scinti #frac{%0.2g}{#sqrt{x}} & #frac{%0.2g}{x} & %0.2g",g_s_a,g_s_b,g_s_c), "L");
     g_c->SetMarkerColor(kGreen);
     g_c->SetMarkerStyle(23);
     g_c->SetMarkerSize(1.5);
     g_c->SetLineColor(kGreen);
     g_c->Draw("CP");
     lgd->AddEntry(g_c, Form("Ceren #frac{%0.2g}{#sqrt{x}} & #frac{%0.2g}{x} & %0.2g",g_c_a,g_c_b,g_c_c), "L");
     
     lgd->Draw();
     //Canvas->Print("h.png",".png");
     
    
    
    return;
}
