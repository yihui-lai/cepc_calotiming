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
  this -> GetTree() -> Branch("depositedEnergyHCALAct",       &this->depositedEnergyHCALAct,               "depositedEnergyHCALAct/F");
  this -> GetTree() -> Branch("depositedEnergyHCALPas",       &this->depositedEnergyHCALPas,               "depositedEnergyHCALPas/F");

  this -> GetTree() -> Branch("depositedEnergyWorld",     	&this->depositedEnergyWorld,         	"depositedEnergyWorld/F");
  this -> GetTree() -> Branch("depositedEnergyTimingGap",     	&this->depositedEnergyTimingGap,         	"depositedEnergyTimingGap/F");
  this -> GetTree() -> Branch("depositedEnergyServices",     	&this->depositedEnergyServices,         	"depositedEnergyServices/F");
  this -> GetTree() -> Branch("depositedEnergyEcalGap",     	&this->depositedEnergyEcalGap,         	"depositedEnergyEcalGap/F");
  this -> GetTree() -> Branch("depositedEnergyEcalDet",     	&this->depositedEnergyEcalDet,         	"depositedEnergyEcalDet/F");
  this -> GetTree() -> Branch("depositedEnergySolenoid",     	&this->depositedEnergySolenoid,         	"depositedEnergySolenoid/F");

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




  h_phot_cer_lambda_Timing_f = new TH1F("h_phot_cer_lambda_Timing_f","",1250, 0.,1250.);
  h_phot_cer_lambda_Timing_r = new TH1F("h_phot_cer_lambda_Timing_r","",1250, 0.,1250.);
  h_phot_cer_lambda_ECAL_f    = new TH1F("h_phot_cer_lambda_ECAL_f","",1250, 0.,1250.);
  h_phot_cer_lambda_ECAL_r    = new TH1F("h_phot_cer_lambda_ECAL_r","",1250, 0.,1250.);
  h_phot_cer_lambda_HCAL	   = new TH1F("h_phot_cer_lambda_HCAL","",1250, 0.,1250.);
  h_phot_cer_parentID = new TH1F("h_phot_cer_parentID","",2000,-1000,1000);


  
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
  h_phot_cer_parentID->Write();

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
  depositedEnergyHCALAct = 0.;
  depositedEnergyHCALPas = 0.;

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
  depositedEnergyServices = 0.;
  depositedEnergyTimingGap = 0.;
  depositedEnergyEcalGap = 0.;
  depositedEnergyEcalDet = 0.;
  depositedEnergySolenoid = 0.;
  

  
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




}
