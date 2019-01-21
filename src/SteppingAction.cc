#include "SteppingAction.hh"
#include "TrackingAction.hh"
#include "DetectorConstruction.hh"
#include "TString.h"
#include "TRandom3.h"
#include "TCint.h"
#include "G4Step.hh"
#include "G4Track.hh"
#include "G4SteppingManager.hh"
#include <time.h>

#include <iostream>
#include <fstream>
#include <vector>
#include "TTree.h"

//long int CreateSeed();



using namespace std;
using namespace CLHEP;
/*
SteppingAction::SteppingAction (const string& configFileName)
{
  //---------------------------------------
  //------------- Parameters --------------
  //---------------------------------------
  
  ConfigFile config (configFileName) ;

  config.readInto(core_material, "core_material"); 

  if (core_material == 0)
  {
	  config.readInto(toy_ly,	"toy_ly");
	  config.readInto(toy_decay,	"toy_decay");
	  config.readInto(toy_rise,	"toy_rise");
  }
  

}*/



int to_int (string name)
{
  int Result ;             // int which will contain the result
  stringstream convert (name) ;
  string dummy ;           
  convert >> dummy ;       
  convert >> Result ;
  return Result ;
}


//---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- 


SteppingAction::SteppingAction (DetectorConstruction* detectorConstruction,
                                const G4int& scint, const G4int& cher):
  fDetectorConstruction(detectorConstruction),
  propagateScintillation(scint),
  propagateCerenkov(cher)
{}


// ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----


SteppingAction::~SteppingAction ()
{}


// ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----


void SteppingAction::UserSteppingAction (const G4Step * theStep)
{
 
  G4Track* theTrack = theStep->GetTrack () ;
  
  const G4ThreeVector& theTrackDirection = theTrack->GetMomentumDirection();
  const G4ThreeVector& theTrackVertexDirection = theTrack->GetVertexMomentumDirection();
  
  G4int trackID = theTrack->GetTrackID();
  TrackInformation* theTrackInfo = (TrackInformation*)(theTrack->GetUserInformation());
  G4ParticleDefinition* particleType = theTrack->GetDefinition () ;
  
  G4StepPoint * thePrePoint  = theStep->GetPreStepPoint () ;
  G4StepPoint * thePostPoint = theStep->GetPostStepPoint () ;
  const G4ThreeVector & thePrePosition  = thePrePoint->GetPosition () ;
  G4VPhysicalVolume * thePrePV  = thePrePoint->GetPhysicalVolume () ;
  G4VPhysicalVolume * thePostPV = thePostPoint->GetPhysicalVolume () ;
  G4String thePrePVName  = "" ; if ( thePrePV )  thePrePVName  = thePrePV  -> GetName () ;
  G4String thePostPVName = "" ; if ( thePostPV ) thePostPVName = thePostPV -> GetName () ;
  G4VSolid* thePreS = thePrePV->GetLogicalVolume()->GetSolid();
    
  G4int nStep = theTrack -> GetCurrentStepNumber();
  
//        cout << " step length = " << theStep->GetStepLength() << endl;
  //-------------
 
  // get position
  G4double global_x = thePrePosition.x()/mm;
  G4double global_y = thePrePosition.y()/mm;
  G4double global_z = thePrePosition.z()/mm;
  
  
  // optical photon

  if( particleType == G4OpticalPhoton::OpticalPhotonDefinition() )
  {
  /*
    G4String processName = theTrack->GetCreatorProcess()->GetProcessName();
        
    //----------------------------
    // count photons at production
    if( ( theTrack->GetLogicalVolumeAtVertex()->GetName().contains("core") ) &&
        (nStep == 1) && (processName == "Scintillation") )
    {

//       cout << " in log volume " << endl;
      //save only prompt photons
      if (thePrePoint->GetGlobalTime()/picosecond>2500) 	theTrack->SetTrackStatus(fKillTrackAndSecondaries);

      else if (thePrePVName == "corePV") 
      {
	CreateTree::Instance()->tot_phot_sci += 1;
//	CreateTree::Instance()->h_phot_sci_time   -> Fill( thePrePoint->GetGlobalTime()/ns );	//ATTENTION** error: in previous ntuple the time was in ns!
	CreateTree::Instance()->h_phot_sci_time   -> Fill( thePrePoint->GetGlobalTime()/picosecond );
	CreateTree::Instance()->h_phot_sci_angleAtProduction -> Fill( G4ThreeVector(0.,0.,1.).angle(theTrackVertexDirection) *57.2958);
        CreateTree::Instance()->h_phot_sci_lambda -> Fill( MyMaterials::fromEvToNm(theTrack->GetTotalEnergy()/eV) );
      }
           
      if( !propagateScintillation ) theTrack->SetTrackStatus(fKillTrackAndSecondaries);      

    }
    
        
    if( ( theTrack->GetLogicalVolumeAtVertex()->GetName().contains("core")) &&
        (nStep == 1) && (processName == "Cerenkov") )
    {
      //kill very long wavelengths
      if (MyMaterials::fromEvToNm(theTrack->GetTotalEnergy()/eV) > 1200)  theTrack->SetTrackStatus(fKillTrackAndSecondaries); 
      else if (thePrePVName == "corePV" && MyMaterials::fromEvToNm(theTrack->GetTotalEnergy()/eV) < 1200)
      {      
	CreateTree::Instance()->tot_phot_cer += 1;
        CreateTree::Instance()->h_phot_cer_time   -> Fill( thePrePoint->GetGlobalTime()/picosecond );
        CreateTree::Instance()->h_phot_cer_angleAtProduction -> Fill( G4ThreeVector(0.,0.,1.).angle(theTrackVertexDirection)*57.2958 );
        CreateTree::Instance()->h_phot_cer_lambda -> Fill( MyMaterials::fromEvToNm(theTrack->GetTotalEnergy()/eV) );

//	CreateTree::Instance()->time_prod_cher.push_back(thePrePoint->GetGlobalTime()/picosecond );
//	CreateTree::Instance()->lambda_prod_cher.push_back(MyMaterials::fromEvToNm(theTrack->GetTotalEnergy()/eV));
      }

      if( !propagateCerenkov ) theTrack->SetTrackStatus(fKillTrackAndSecondaries);      

    }
    
    
    //----------------------------
    // count photons at fiber exit
    
    if(  theTrack->GetLogicalVolumeAtVertex()->GetName().contains("core")  &&        processName == "Scintillation"  &&
        ((thePrePVName.contains("gapLayerPV") && thePostPVName.contains("gapPV")) || (thePrePVName.contains("gapLayerPV_ref") && thePostPVName.contains("gapPV_ref"))) )
    {
      CreateTree::Instance()->tot_gap_phot_sci += 1;
      if (thePrePVName.contains("gapLayerPV")) 
      {
	//uncomment line below to count photons at coupling gap
	CreateTree::Instance()->time_ext_scint.push_back(thePrePoint->GetGlobalTime()/picosecond );
//	CreateTree::Instance()->lambda_ext_scint.push_back(MyMaterials::fromEvToNm(theTrack->GetTotalEnergy()/eV));

       CreateTree::Instance()->h_phot_sci_gap_lambda -> Fill( MyMaterials::fromEvToNm(theTrack->GetTotalEnergy()/eV) );
       CreateTree::Instance()->h_phot_sci_gap_time   -> Fill( thePrePoint->GetGlobalTime()/picosecond );
       CreateTree::Instance()->h_phot_sci_gap_angleAtProduction -> Fill( (G4ThreeVector(0.,0.,1.).angle(theTrackVertexDirection))*57.2958 );
       CreateTree::Instance()->h_phot_sci_gap_angleWithSurfNormal -> Fill( (G4ThreeVector(0.,0.,1.).angle(theTrackDirection))*57.2958 );
      }
//      else if (thePrePVName == "gapLayerPV_ref") CreateTree::Instance()->time_ext_scint_ref.push_back(thePrePoint->GetGlobalTime()/picosecond );
      // if you do not want to kill a photon once it exits the fiber, comment here below
//      theTrack->SetTrackStatus(fKillTrackAndSecondaries);      
    }
    
    if( ( theTrack->GetLogicalVolumeAtVertex()->GetName().contains("core") ||
          theTrack->GetLogicalVolumeAtVertex()->GetName().contains("capillary") ||
          theTrack->GetLogicalVolumeAtVertex()->GetName().contains("cladding") ) &&
        (processName == "Cerenkov") &&
        ( (thePrePVName.contains("gapLayerPV") && thePostPVName.contains("gapPV")) || (thePrePVName.contains("gapLayerPV_ref") && thePostPVName.contains("gapPV_ref")) ) 
      )
    {
      CreateTree::Instance()->tot_gap_phot_cer += 1;
      // if you do not want to kill a photon once it exits the fiber, comment here below
//      theTrack->SetTrackStatus(fKillTrackAndSecondaries);
      
      if (thePrePVName.contains("gapLayerPV"))
      {
	//uncomment line below to count photons at coupling gap
	CreateTree::Instance()->time_ext_cher.push_back(thePrePoint->GetGlobalTime()/picosecond );
//	CreateTree::Instance()->lambda_ext_cher.push_back(MyMaterials::fromEvToNm(theTrack->GetTotalEnergy()/eV));

        CreateTree::Instance()->h_phot_cer_gap_lambda -> Fill( MyMaterials::fromEvToNm(theTrack->GetTotalEnergy()/eV) );
        CreateTree::Instance()->h_phot_cer_gap_time   -> Fill( thePrePoint->GetGlobalTime()/picosecond );
        CreateTree::Instance()->h_phot_cer_gap_angleAtProduction -> Fill( (G4ThreeVector(0.,0.,1.).angle(theTrackVertexDirection))*57.2958 );
	CreateTree::Instance()->h_phot_cer_gap_angleWithSurfNormal -> Fill( (G4ThreeVector(0.,0.,1.).angle(theTrackDirection))*57.2958 );
      } 
//      else if (thePrePVName == "gapLayerPV_ref") CreateTree::Instance()->time_ext_cher_ref.push_back(thePrePoint->GetGlobalTime()/picosecond );
      

    }
    
   
    
    
    //------------------------------
    // count photons at the detector
    
    if( ( theTrack->GetLogicalVolumeAtVertex()->GetName().contains("core") ) &&
        (processName == "Scintillation") &&
        (thePrePVName.contains("detLayerPV")) && (thePostPVName.contains("detPV")) )
    {
      CreateTree::Instance()->tot_det_phot_sci += 1;
      CreateTree::Instance()->h_phot_sci_det_lambda -> Fill( MyMaterials::fromEvToNm(theTrack->GetTotalEnergy()/eV) );
      CreateTree::Instance()->h_phot_sci_det_time   -> Fill( thePrePoint->GetGlobalTime()/picosecond );
      CreateTree::Instance()->h_phot_sci_det_angleAtProduction -> Fill( (G4ThreeVector(0.,0.,1.).angle(theTrackVertexDirection))*57.2958 );
      CreateTree::Instance()->h_phot_sci_det_angleWithSurfNormal -> Fill( (G4ThreeVector(0.,0.,1.).angle(theTrackDirection))*57.2958 );
      //comment line below if standard SiPM
//      CreateTree::Instance()->time_ext_scint.push_back(thePrePoint->GetGlobalTime()/picosecond );
      // if you do not want to kill a photon once it enters the detector, comment here below
      theTrack->SetTrackStatus(fKillTrackAndSecondaries);
    }
    
    if( ( theTrack->GetLogicalVolumeAtVertex()->GetName().contains("core") ||
          theTrack->GetLogicalVolumeAtVertex()->GetName().contains("capillary") ||
          theTrack->GetLogicalVolumeAtVertex()->GetName().contains("cladding") ) &&
        (processName == "Cerenkov") &&
        (thePrePVName.contains("detLayerPV")) && (thePostPVName.contains("detPV")) )
    {
      CreateTree::Instance()->tot_det_phot_cer += 1;

      CreateTree::Instance()->h_phot_cer_det_lambda -> Fill( MyMaterials::fromEvToNm(theTrack->GetTotalEnergy()/eV) );
      CreateTree::Instance()->h_phot_cer_det_time   -> Fill( thePrePoint->GetGlobalTime()/picosecond );
      CreateTree::Instance()->h_phot_cer_det_angleAtProduction -> Fill( (G4ThreeVector(0.,0.,1.).angle(theTrackVertexDirection))*57.2958 );
      CreateTree::Instance()->h_phot_cer_det_angleWithSurfNormal -> Fill( (G4ThreeVector(0.,0.,1.).angle(theTrackDirection))*57.2958 );
      //comment line below if standard SiPM
//      CreateTree::Instance()->time_ext_cher.push_back(thePrePoint->GetGlobalTime()/picosecond );
//    
      // if you do not want to kill a photon once it enters the detector, comment here below
      theTrack->SetTrackStatus(fKillTrackAndSecondaries);
    }
    
    
    
    
    
    if( (theTrack->GetLogicalVolumeAtVertex()->GetName().contains("core")) && (nStep == 1) )
    {    
      //----------------------------------------------------------
      // storing time, energy and position at gap with fast timing
      Photon ph;
      ph.position.SetX(global_x);
      ph.position.SetY(global_y);
      ph.position.SetZ(global_z);
      ph.direction.SetX(theTrack->GetVertexMomentumDirection().x());
      ph.direction.SetY(theTrack->GetVertexMomentumDirection().y());
      ph.direction.SetZ(theTrack->GetVertexMomentumDirection().z());
      ph.dist = (global_z/(0.5*fiber_length));
      ph.energy = theTrack->GetTotalEnergy()/eV;
      
      Fiber* fib = fDetectorConstruction -> GetFiber();
      std::map<int,Travel> trc = GetTimeAndProbability(ph,fib,theTrackInfo->GetParticleProdTime());
      
      for(unsigned int it = 0; it < CreateTree::Instance()->attLengths->size(); ++it)
      {
        int attLength = int( CreateTree::Instance()->attLengths->at(it) );
        
        if( trc[attLength].prob[0] < 1.E-09 ) theTrack->SetTrackStatus(fKillTrackAndSecondaries);      
        
        for(int it2 = 0; it2 < 3; ++it2)
        {
          CreateTree::Instance()->tot_gap_photFast_cer->at(it) += trc[attLength].prob[it2];
          
          //CreateTree::Instance()->h_photFast_cer_gap_lambda[attLength] -> Fill( MyMaterials::fromEvToNm(theTrack->GetTotalEnergy()/eV), trc[attLength].prob[it2] );
          //CreateTree::Instance()->h_photFast_cer_gap_E[attLength]      -> Fill( theTrack->GetTotalEnergy()/eV, trc[attLength].prob[it2] );
          //CreateTree::Instance()->h_photFast_cer_gap_time[attLength]   -> Fill( trc[attLength].time[it2], trc[attLength].prob[it2] );
        }
      }
    }
  */  
  } // optical photon
  
  
  // non optical photon
  else
  {
    //G4cout << ">>> begin non optical photon" << G4endl;
    
    G4double energy = theStep->GetTotalEnergyDeposit() - theStep->GetNonIonizingEnergyDeposit();
    if ( energy == 0. ) return;
    
    CreateTree::Instance() -> depositedEnergyTotal += energy/GeV;
    
    if( thePrePVName.contains("absorber") )
    {
      CreateTree::Instance()->depositedEnergyAbsorber += energy/GeV;
    }

    //timing
    if( thePrePVName.contains("corePV_front") )
    {
      CreateTree::Instance()->depositedEnergyTiming_f += energy/GeV;
      for (int iBar = 0; iBar<18; iBar++)
      {
	if (thePrePVName == Form("corePV_front_%d", iBar)) CreateTree::Instance()->Edep_Timing_f_ch[iBar] += energy/GeV;
      }
    }
    if( thePrePVName.contains("corePV_rear") )
    {
      CreateTree::Instance()->depositedEnergyTiming_r += energy/GeV;
      for (int iBar = 0; iBar<18; iBar++)
      {
	if (thePrePVName == Form("corePV_rear_%d", iBar)) CreateTree::Instance()->Edep_Timing_r_ch[iBar] += energy/GeV;
      }	
    }

    //ecal
    if( thePrePVName.contains("ecalCrystalP_f") )
    {
      CreateTree::Instance()->depositedEnergyECAL_f += energy/GeV;
      for (int iCh = 0; iCh<25; iCh++)
      {
	if (thePrePVName == Form("ecalCrystalP_f_%d", iCh)) CreateTree::Instance()->Edep_ECAL_f_ch[iCh] += energy/GeV;
      }
    }
    if( thePrePVName.contains("ecalCrystalP_r") )
    {
      CreateTree::Instance()->depositedEnergyECAL_r += energy/GeV;
      for (int iCh = 0; iCh<25; iCh++)
      {
	if (thePrePVName == Form("ecalCrystalP_r_%d", iCh)) CreateTree::Instance()->Edep_ECAL_r_ch[iCh] += energy/GeV;
      }
    }    


    if( thePrePVName.contains("world") )
    {
      CreateTree::Instance() -> depositedEnergyWorld += energy/GeV;
    }
    
    //G4cout << ">>> end non optical photon" << G4endl;
  } // non optical photon
  
  
  return ;
}

