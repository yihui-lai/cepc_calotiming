#include "CreateTree.hh"
#include <algorithm>

using namespace std ;

CreateTree* CreateTree::fInstance = NULL ;


// ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----


CreateTree::CreateTree (TString name)
{
  if( fInstance )
  {
    return ;
  }
  
  this -> fInstance = this ;
  this -> fname     = name ;
  this -> ftree     = new TTree (name,name) ;
  
  this -> GetTree ()->Branch ("Event", &this->Event, "Event/I") ;
  
  inputInitialPosition = new vector<float>(3,0.); 
  inputMomentum =        new vector<float>(4,0.); 
  this -> GetTree() -> Branch("inputInitialPosition", "vector<float>", &inputInitialPosition);
  this -> GetTree() -> Branch("inputMomentum",        "vector<float>", &inputMomentum);
  
  //integrated per longitudinal layer 
  this -> GetTree() -> Branch("depositedEnergyTotal",     &this->depositedEnergyTotal,         "depositedEnergyTotal/F");

  this -> GetTree() -> Branch("depositedEnergyAbsorber",      &this->depositedEnergyAbsorber,           "depositedEnergyAbsorber/F");

  this -> GetTree() -> Branch("depositedEnergyTiming_f",      &this->depositedEnergyTiming_f,           "depositedEnergyTiming_f/F");
  this -> GetTree() -> Branch("depositedEnergyTiming_r",      &this->depositedEnergyTiming_r,           "depositedEnergyTiming_r/F");

  this -> GetTree() -> Branch("depositedEnergyECAL_f",      &this->depositedEnergyECAL_f,           "depositedEnergyECAL_f/F");
  this -> GetTree() -> Branch("depositedEnergyECAL_r",      &this->depositedEnergyECAL_r,           "depositedEnergyECAL_r/F");

  this -> GetTree() -> Branch("depositedEnergyWorld",     &this->depositedEnergyWorld,         "depositedEnergyWorld/F");

  //single channels
  this -> GetTree() -> Branch("Edep_Timing_f_ch", &this->Edep_Timing_f_ch,    "Edep_Timing_f_ch[18]/F");
  this -> GetTree() -> Branch("Edep_Timing_r_ch", &this->Edep_Timing_r_ch,    "Edep_Timing_r_ch[18]/F");

  this -> GetTree() -> Branch("Edep_ECAL_f_ch",   &this->Edep_ECAL_f_ch,      "Edep_ECAL_f_ch[25]/F");
  this -> GetTree() -> Branch("Edep_ECAL_r_ch",   &this->Edep_ECAL_r_ch,      "Edep_ECAL_r_ch[25]/F");
  
/*
  this -> GetTree() -> Branch("E_dep_f",   &E_dep_f);
  this -> GetTree() -> Branch("E_dep_time_f",   &E_dep_time_f);
  this -> GetTree() -> Branch("E_dep_r",   &E_dep_r);
  this -> GetTree() -> Branch("E_dep_time_r",   &E_dep_time_r);

  this -> GetTree() -> Branch("time_prod_scint", &time_prod_scint);
  this -> GetTree() -> Branch("time_prod_cher",  &time_prod_cher);
  this -> GetTree() -> Branch("time_ext_scint",  &time_ext_scint);
  this -> GetTree() -> Branch("time_ext_cher", 	 &time_ext_cher);

  this -> GetTree() -> Branch("lambda_prod_scint", &lambda_prod_scint);
  this -> GetTree() -> Branch("lambda_prod_cher",  &lambda_prod_cher);
  this -> GetTree() -> Branch("lambda_ext_scint",  &lambda_ext_scint);
  this -> GetTree() -> Branch("lambda_ext_cher",   &lambda_ext_cher);
  
  this -> GetTree() -> Branch("time_prod_scint_ref", &time_prod_scint_ref);
  this -> GetTree() -> Branch("time_prod_cher_ref",  &time_prod_cher_ref);
  this -> GetTree() -> Branch("time_ext_scint_ref",  &time_ext_scint_ref);
  this -> GetTree() -> Branch("time_ext_cher_ref",   &time_ext_cher_ref);
    
  this -> GetTree() -> Branch("tot_phot_sci",        &this->tot_phot_sci,               "tot_phot_sci/I");
  this -> GetTree() -> Branch("tot_phot_cer",        &this->tot_phot_cer,               "tot_phot_cer/I");
  this -> GetTree() -> Branch("tot_gap_phot_sci",    &this->tot_gap_phot_sci,       "tot_gap_phot_sci/I");
  this -> GetTree() -> Branch("tot_gap_phot_cer",    &this->tot_gap_phot_cer,       "tot_gap_phot_cer/I");
  this -> GetTree() -> Branch("tot_det_phot_sci",    &this->tot_det_phot_sci,       "tot_det_phot_sci/I");
  this -> GetTree() -> Branch("tot_det_phot_cer",    &this->tot_det_phot_cer,       "tot_det_phot_cer/I");
  
  h_phot_sci_lambda = new TH1F("h_phot_sci_lambda","",1250,0.,1250.);
  h_phot_sci_time = new TH1F("h_phot_sci_time","",100000,0.,5000.);
  h_phot_sci_angleAtProduction = new TH1F("h_phot_sci_angleAtProduction","",2000,0.,360.);

  h_phot_cer_lambda = new TH1F("h_phot_cer_lambda","",1250, 0.,1250.);
  h_phot_cer_time = new TH1F("h_phot_cer_time","",100000,0.,5000.);
  h_phot_cer_angleAtProduction = new TH1F("h_phot_cer_angleAtProduction","",2000,0.,360.);
  
  h_phot_sci_gap_lambda = new TH1F("h_phot_sci_gap_lambda","",1250, 0.,1250.);
  h_phot_sci_gap_time = new TH1F("h_phot_sci_gap_time","",100000,0.,5000.);
  h_phot_sci_gap_angleAtProduction = new TH1F("h_phot_sci_gap_angleAtProduction","",2000,0.,360.);
  h_phot_sci_gap_angleWithSurfNormal = new TH1F("h_phot_sci_gap_angleWithSurfNormal","",2000,0.,360.);

  h_phot_cer_gap_lambda = new TH1F("h_phot_cer_gap_lambda","",1250,0.,1250.);
  h_phot_cer_gap_time = new TH1F("h_phot_cer_gap_time","",100000,0.,5000.);
  h_phot_cer_gap_angleAtProduction = new TH1F("h_phot_cer_gap_angleAtProduction","",2000,0.,360.);
  h_phot_cer_gap_angleWithSurfNormal = new TH1F("h_phot_cer_gap_angleWithSurfNormal","",2000,0.,360.);

  h_phot_sci_det_lambda = new TH1F("h_phot_sci_gap_lambda","",1250, 0.,1250.);
  h_phot_sci_det_time = new TH1F("h_phot_sci_det_time","",100000,0.,5000.);
  h_phot_sci_det_angleAtProduction = new TH1F("h_phot_sci_det_angleAtProduction","",2000,0.,360.);
  h_phot_sci_det_angleWithSurfNormal = new TH1F("h_phot_sci_det_angleWithSurfNormal","",2000,0.,360.);

  h_phot_cer_det_lambda = new TH1F("h_phot_cer_det_lambda","",1250,0.,1250.);
  h_phot_cer_det_time = new TH1F("h_phot_cer_det_time","",100000,0.,5000.);
  h_phot_cer_det_angleAtProduction = new TH1F("h_phot_cer_det_angleAtProduction","",2000,0.,360.);
  h_phot_cer_det_angleWithSurfNormal = new TH1F("h_phot_cer_det_angleWithSurfNormal","",2000,0.,360.);
  
*/
  
  this -> Clear() ;
}

// ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----



CreateTree::~CreateTree()
{}

// ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----



int CreateTree::Fill() 
{ 
  return this -> GetTree() -> Fill(); 
}

// ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----



bool CreateTree::Write(TFile * outfile)
{
  outfile -> cd();
  ftree -> Write();
  /*
  h_phot_sci_lambda->Write();
  h_phot_sci_time->Write();
  h_phot_sci_angleAtProduction -> Write();

  h_phot_cer_lambda->Write();
  h_phot_cer_time->Write();
  h_phot_cer_angleAtProduction -> Write();

  h_phot_sci_gap_lambda->Write();
  h_phot_sci_gap_time->Write();
  h_phot_sci_gap_angleAtProduction->Write();
  h_phot_sci_gap_angleWithSurfNormal->Write();

  h_phot_cer_gap_lambda->Write();
  h_phot_cer_gap_time->Write();
  h_phot_cer_gap_angleAtProduction->Write();
  h_phot_cer_gap_angleWithSurfNormal->Write();

  h_phot_sci_det_lambda->Write();
  h_phot_sci_det_time->Write();
  h_phot_sci_det_angleAtProduction->Write();
  h_phot_sci_det_angleWithSurfNormal->Write();

  h_phot_cer_det_lambda->Write();
  h_phot_cer_det_time->Write();
  h_phot_cer_det_angleAtProduction->Write();
  h_phot_cer_det_angleWithSurfNormal->Write();
  */
  return true ;
}

// ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----



void CreateTree::Clear()
{
  Event	= 0;
  
  depositedEnergyTotal = 0.;

  depositedEnergyAbsorber = 0.;
  depositedEnergyTiming_f = 0.;
  depositedEnergyTiming_r = 0.;
  depositedEnergyECAL_f = 0.;
  depositedEnergyECAL_r = 0.;

  for (int iBar = 0; iBar<18; iBar++)
  {
	Edep_Timing_f_ch[iBar] = 0.;
	Edep_Timing_r_ch[iBar] = 0.;
  }
  for (int iCh = 0; iCh<25; iCh++)
  {
	Edep_ECAL_f_ch[iCh] = 0.;
	Edep_ECAL_r_ch[iCh] = 0.;
  }

  depositedEnergyWorld = 0.;
  

  
  for (int i = 0 ; i < 3 ; ++i) 
  {
    inputInitialPosition -> at(i) = 0.;
  }
  for (int i = 0 ; i < 4 ; ++i) 
  {
    inputMomentum ->at(i) = 0.;
  }




  /*
  E_dep_f.clear();
  E_dep_time_f.clear();
  E_dep_r.clear();
  E_dep_time_r.clear();

  tot_phot_sci = 0;
  tot_phot_cer = 0;
  tot_gap_phot_sci = 0;
  tot_gap_phot_cer = 0;
  tot_det_phot_sci = 0;
  tot_det_phot_cer = 0;
  
  time_ext_cher.clear();
  time_ext_scint.clear();
  time_prod_cher.clear();
  time_prod_scint.clear();

  lambda_ext_cher.clear();
  lambda_ext_scint.clear();
  lambda_prod_cher.clear();
  lambda_prod_scint.clear();
  
  time_ext_cher_ref.clear();
  time_ext_scint_ref.clear();
  time_prod_cher_ref.clear();
  time_prod_scint_ref.clear();
*/
}
