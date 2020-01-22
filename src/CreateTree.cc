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


  this -> GetTree() -> Branch("inputTrackerX0",        		&this->inputTrackerX0,               	"inputTrackerX0/F");
  this -> GetTree() -> Branch("inputServiceAlmm",      	&this->inputServiceAlmm,            "inputServiceAlmm/F");
  this -> GetTree() -> Branch("inputTimingThick",      		&this->inputTimingThick,              "inputTimingThick/F");
  this -> GetTree() -> Branch("inputE1Thick",      		&this->inputE1Thick,               	"inputE1Thick/F");
  this -> GetTree() -> Branch("inputE2Thick",      		&this->inputE2Thick,               	"inputE2Thick/F");
  this -> GetTree() -> Branch("inputE1Width",      		&this->inputE1Width,              	"inputE1Width/F");
  this -> GetTree() -> Branch("inputTimingECAL_dist",      	&this->inputTimingECAL_dist,       "inputTimingECAL_dist/F");

  
  inputInitialPosition	= new vector<float>(3,0.); 
  inputMomentum		= new vector<float>(4,0.);
  primaryPosT1		= new vector<float>(3,0.); 
  primaryMomT1		= new vector<float>(4,0.); 
  primaryPosE1		= new vector<float>(3,0.); 
  primaryMomE1		= new vector<float>(4,0.);  

  this -> GetTree() -> Branch("inputInitialPosition",    "vector<float>", &inputInitialPosition);
  this -> GetTree() -> Branch("inputMomentum",        "vector<float>", &inputMomentum);
  this -> GetTree() -> Branch("primaryPosT1", 	"vector<float>", &primaryPosT1);
  this -> GetTree() -> Branch("primaryMomT1",	"vector<float>", &primaryMomT1);
  this -> GetTree() -> Branch("primaryPosE1", 	"vector<float>", &primaryPosE1);
  this -> GetTree() -> Branch("primaryMomE1",	"vector<float>", &primaryMomE1);  

  this -> GetTree() -> Branch("nTracksT1",        &this->nTracksT1,               "nTracksT1/I");
  this -> GetTree() -> Branch("nTracksT2",        &this->nTracksT2,               "nTracksT2/I");
  this -> GetTree() -> Branch("nTracksE1",        &this->nTracksE1,               "nTracksE1/I");
  this -> GetTree() -> Branch("nTracksE2",        &this->nTracksE2,               "nTracksE2/I");
  this -> GetTree() -> Branch("nTracksTRK",     &this->nTracksTRK,    	  "nTracksTRK[6]/F");


  //integrated per longitudinal layer 
  this -> GetTree() -> Branch("depositedEnergyTotal",     	&this->depositedEnergyTotal,        	"depositedEnergyTotal/F");
  this -> GetTree() -> Branch("depositedEnergyEscapeWorld",     	&this->depositedEnergyEscapeWorld,        	"depositedEnergyEscapeWorld/F");

  this -> GetTree() -> Branch("depositedEnergyTiming_f",  &this->depositedEnergyTiming_f,         "depositedEnergyTiming_f/F");
  this -> GetTree() -> Branch("depositedEnergyTiming_r",  &this->depositedEnergyTiming_r,         "depositedEnergyTiming_r/F");
  this -> GetTree() -> Branch("depositedEnergyECAL_f",     &this->depositedEnergyECAL_f,            "depositedEnergyECAL_f/F");
  this -> GetTree() -> Branch("depositedEnergyECAL_r",     &this->depositedEnergyECAL_r,            "depositedEnergyECAL_r/F");
  this -> GetTree() -> Branch("depositedEnergyHCAL",       &this->depositedEnergyHCAL,               "depositedEnergyHCAL/F");

  this -> GetTree() -> Branch("depositedEnergyWorld",     	&this->depositedEnergyWorld,         	"depositedEnergyWorld/F");

  //single channels
  this -> GetTree() -> Branch("Edep_Tracker_layer", 		&this->Edep_Tracker_layer,    		"Edep_Tracker_layer[6]/F");
  this -> GetTree() -> Branch("Edep_Timing_f_ch", 		&this->Edep_Timing_f_ch,    		"Edep_Timing_f_ch[18]/F");
  this -> GetTree() -> Branch("Edep_Timing_r_ch", 		&this->Edep_Timing_r_ch,    		"Edep_Timing_r_ch[18]/F");
  this -> GetTree() -> Branch("Edep_ECAL_f_ch",   		&this->Edep_ECAL_f_ch,      			"Edep_ECAL_f_ch[100]/F");
  this -> GetTree() -> Branch("Edep_ECAL_r_ch",   		&this->Edep_ECAL_r_ch,      			"Edep_ECAL_r_ch[100]/F");
  

  //Cerenkov photons
//  this -> GetTree() -> Branch("tot_phot_sci_Timing",        &this->tot_phot_sci_Timing,               "tot_phot_sci_Timing/I");
  this -> GetTree() -> Branch("tot_phot_cer_Timing_f",     &this->tot_phot_cer_Timing_f,            "tot_phot_cer_Timing_f/I");
  this -> GetTree() -> Branch("tot_phot_cer_Timing_r",     &this->tot_phot_cer_Timing_r,            "tot_phot_cer_Timing_r/I");
  this -> GetTree() -> Branch("tot_phot_cer_ECAL_f",        &this->tot_phot_cer_ECAL_f,               "tot_phot_cer_ECAL_f/I");
  this -> GetTree() -> Branch("tot_phot_cer_ECAL_r",        &this->tot_phot_cer_ECAL_r,               "tot_phot_cer_ECAL_r/I");
  this -> GetTree() -> Branch("tot_phot_cer_HCAL",        &this->tot_phot_cer_HCAL,               "tot_phot_cer_HCAL/I");


//  h_phot_sci_lambda = new TH1F("h_phot_sci_lambda","",1250,0.,1250.);
//  h_phot_sci_time = new TH1F("h_phot_sci_time","",100000,0.,5000.);
//  h_phot_sci_angleAtProduction = new TH1F("h_phot_sci_angleAtProduction","",2000,0.,360.);

  h_phot_cer_lambda_Timing_f = new TH1F("h_phot_cer_lambda_Timing_f","",1250, 0.,1250.);
  h_phot_cer_lambda_Timing_r = new TH1F("h_phot_cer_lambda_Timing_r","",1250, 0.,1250.);
  h_phot_cer_lambda_ECAL_f    = new TH1F("h_phot_cer_lambda_ECAL_f","",1250, 0.,1250.);
  h_phot_cer_lambda_ECAL_r    = new TH1F("h_phot_cer_lambda_ECAL_r","",1250, 0.,1250.);
  h_phot_cer_lambda_HCAL	   = new TH1F("h_phot_cer_lambda_HCAL","",1250, 0.,1250.);

//  h_phot_cer_time = new TH1F("h_phot_cer_time","",100000,0.,5000.);
//  h_phot_cer_angleAtProduction = new TH1F("h_phot_cer_angleAtProduction","",2000,0.,360.);

  /*
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
    
  this -> GetTree() -> Branch("tot_gap_phot_sci",    &this->tot_gap_phot_sci,       "tot_gap_phot_sci/I");
  this -> GetTree() -> Branch("tot_gap_phot_cer",    &this->tot_gap_phot_cer,       "tot_gap_phot_cer/I");
  this -> GetTree() -> Branch("tot_det_phot_sci",    &this->tot_det_phot_sci,       "tot_det_phot_sci/I");
  this -> GetTree() -> Branch("tot_det_phot_cer",    &this->tot_det_phot_cer,       "tot_det_phot_cer/I");
  

  
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


  h_phot_cer_lambda_Timing_f ->Write();
  h_phot_cer_lambda_Timing_r ->Write();
  h_phot_cer_lambda_ECAL_f    ->Write();
  h_phot_cer_lambda_ECAL_r    ->Write();
  h_phot_cer_lambda_HCAL	   ->Write();

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
  depositedEnergyEscapeWorld=0.;
  nTracksT1 = 0;
  nTracksT2 = 0;
  nTracksE1 = 0;
  nTracksE2 = 0;

  for (int iLayer = 0; iLayer<6; iLayer++)
  {
	nTracksTRK[iLayer] = 0;
  }

  depositedEnergyTiming_f = 0.;
  depositedEnergyTiming_r = 0.;
  depositedEnergyECAL_f = 0.;
  depositedEnergyECAL_r = 0.;
  depositedEnergyHCAL = 0.;

  tot_phot_cer_Timing_f = 0.;
  tot_phot_cer_Timing_r = 0.;
  tot_phot_cer_ECAL_f = 0.;
  tot_phot_cer_ECAL_r = 0.;
  tot_phot_cer_HCAL = 0.;

  for (int iLayer = 0; iLayer<6; iLayer++)
  {
	Edep_Tracker_layer[iLayer] = 0.;
  }

  for (int iBar = 0; iBar<18; iBar++)
  {
	Edep_Timing_f_ch[iBar] = 0.;
	Edep_Timing_r_ch[iBar] = 0.;
  }
  for (int iCh = 0; iCh<100; iCh++)
  {
	Edep_ECAL_f_ch[iCh] = 0.;
	Edep_ECAL_r_ch[iCh] = 0.;
  }

  depositedEnergyWorld = 0.;
  

  
  for (int i = 0 ; i < 3 ; ++i) 
  {
    inputInitialPosition -> at(i) = 0.;
    primaryPosT1 -> at(i) = 0.;
    primaryPosE1 -> at(i) = 0.;

  }
  for (int i = 0 ; i < 4 ; ++i) 
  {
    inputMomentum ->at(i) = 0.;
    primaryMomT1 -> at(i) = 0.;
    primaryMomE1 -> at(i) = 0.;
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
