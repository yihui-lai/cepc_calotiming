//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes, nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
// $Id: DetectorConstruction.cc, v 1.18 2010-10-23 19:27:38 gum Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- 
//---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- 

#include "DetectorConstruction.hh"
#include "G4MagneticField.hh"
#include "G4UniformMagField.hh"
#include "G4FieldManager.hh"
#include "G4TransportationManager.hh"
#include "G4MagIntegratorStepper.hh"
#include "G4Mag_UsualEqRhs.hh"
#include "G4ClassicalRK4.hh"
#include "G4ExplicitEuler.hh"
#include "G4ChordFinder.hh"
#include "G4EqMagElectricField.hh"
#include "G4PropagatorInField.hh"
#include "G4Material.hh"
#include "G4Box.hh"
#include "G4Trd.hh"
#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4PVReplica.hh"
#include "G4SubtractionSolid.hh"
#include "G4GeometryTolerance.hh"
#include "G4GeometryManager.hh"
#include "G4NistManager.hh"
#include "G4VisAttributes.hh"
#include "G4Colour.hh"
#include "G4SDManager.hh"

#include "G4UserLimits.hh"


#include "DetectorConstruction.hh"
#include <G4TransportationManager.hh>
#include <G4MagneticField.hh>
#include <G4UniformMagField.hh>
#include <G4FieldManager.hh>
#include "CreateTree.hh"
#include <algorithm>
#include <string>
#include <sstream>

using namespace CLHEP;



DetectorConstruction::DetectorConstruction (const string& configFileName)
{
  //---------------------------------------
  //------------- Parameters --------------
  //---------------------------------------
  
  ConfigFile config (configFileName) ;
  
  config.readInto(checkOverlaps,"checkOverlaps");
  
  config.readInto(world_material,	"world_material");
  config.readInto(fibre_length,  	"fibre_length");
  config.readInto(detector,      	"detector");
  
  config.readInto(core_radius_x,   "core_radius_x");
  config.readInto(core_radius_y,   "core_radius_y");
  config.readInto(core_material,   "core_material");
  config.readInto(core_rIndex,     "core_rIndex");
  config.readInto(core_absLength,  "core_absLength");
  
  config.readInto(gap_l,         "gap_l");
  config.readInto(gap_size_x,    "gap_size_x");
  config.readInto(gap_size_y,    "gap_size_y");
  config.readInto(gap_material,  "gap_material");
  
  config.readInto(det_l,         "det_l");
  config.readInto(det_size_x,    "det_size_x");
  config.readInto(det_size_y,    "det_size_y");
  config.readInto(det_material,  "det_material");
  
  config.readInto(depth,	"depth");
  config.readInto(cryst_dist,	"cryst_dist");
  config.readInto(abs_thick,	"abs_thick");


  config.readInto(ecal_material, 	"ecal_material");
  config.readInto(ecal_front_length,	"ecal_front_length");
  config.readInto(ecal_rear_length,	"ecal_rear_length");
  config.readInto(ecal_front_face,	"ecal_front_face");
  config.readInto(ecal_rear_face,	"ecal_rear_face");
  config.readInto(ecal_timing_distance,	"ecal_timing_distance");
  config.readInto(ecal_det_size,	"ecal_det_size");

  
  if (core_material == 0)
  {
	  config.readInto(toy_ly,	"toy_ly");
	  config.readInto(toy_decay,	"toy_decay");
	  config.readInto(toy_rise,	"toy_rise");
  }
  
  B_field_intensity = config.read<double>("B_field_intensity") * tesla ;
  
  expHall_x = 300.*cm;
  expHall_y = 300.*cm;
  expHall_z = 300.*cm;
  
  B_field_IsInitialized = false ;
  
  initializeMaterials();
}

//---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- 



DetectorConstruction::~DetectorConstruction ()
{
  delete stepLimit;
}

//---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- 



G4VPhysicalVolume* DetectorConstruction::Construct ()
{
  G4cout << ">>>>>> DetectorConstruction::Construct ()::begin <<<<<<" << G4endl ;
  
  
  //------------------------------------
  //------------- Geometry -------------
  //------------------------------------
  
  
  // The experimental Hall
  // ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----
  G4VSolid * worldS = new G4Box ("worldS", 0.5 * expHall_x, 0.5 * expHall_y, 0.5 * expHall_z) ;
  G4LogicalVolume * worldLV = new G4LogicalVolume (worldS, WoMaterial, "worldLV", 0, 0, 0) ;
  G4VPhysicalVolume * worldPV = new G4PVPlacement (0, G4ThreeVector (), worldLV, "worldPV", 0, false, 0, checkOverlaps) ;


  // the Tracker
  // ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----
  G4VSolid* trackerS;
  if (abs_thick > 0) 
  { 
	  trackerS = new G4Box  ("trackerS", core_radius_x*4, core_radius_y*4, 0.5*abs_thick);    
	  G4LogicalVolume* trackerLV = new G4LogicalVolume (trackerS, MyMaterials::Silicon(), "trackerLV") ;
	  new G4PVPlacement(0, G4ThreeVector(0.,0., -0.5 * (fibre_length + abs_thick) - 3*mm), trackerLV, "trackerPV", worldLV, false, 0, checkOverlaps) ;
  }

  
  // the pre-shower
  // ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----
  G4VSolid* absorberS;
  if (abs_thick > 0) 
  { 
	  absorberS = new G4Box  ("absorberS", core_radius_x*4, core_radius_y*4, 0.5*abs_thick);    
	  G4LogicalVolume* absorberLV = new G4LogicalVolume (absorberS, MyMaterials::Lead(), "absorberLV") ;
	  new G4PVPlacement(0, G4ThreeVector(0.,0., -0.5 * (fibre_length + abs_thick) - 3*mm), absorberLV, "absorberPV", worldLV, false, 0, checkOverlaps) ;
  }
  


  //*********************
  //  TIMING LAYER
  //*********************

  G4RotationMatrix* piRotTiming_front = new G4RotationMatrix;
  G4RotationMatrix* piRotTiming_rear  = new G4RotationMatrix;

//  piRotTiming->rotateY(-pointingAngle*rad*(iX + 0.5));
  piRotTiming_front->rotateX(90*deg);
  piRotTiming_rear->rotateY(90*deg);
//  piRotTiming->rotateY(180*rad);

  int NBARS = 18;
  float bar_spacing = core_radius_x*2;	
  // the crystals
  // ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----
  G4VSolid* coreS;
  coreS = new G4Box  ("coreS", core_radius_x, core_radius_y, 0.5*fibre_length) ;    
  G4LogicalVolume* coreLV = new G4LogicalVolume (coreS, CoMaterial, "coreLV");

  G4PVPlacement* corePV;

  for (int iBar = 0; iBar < NBARS; iBar++)
  { 
 	float x = iBar*bar_spacing - bar_spacing*NBARS/2 + bar_spacing/2;	

	//first layer
  	corePV = new G4PVPlacement(piRotTiming_front, G4ThreeVector(x, 0.,0.), coreLV, Form("corePV_front_%d", iBar), worldLV, false, 0, checkOverlaps) ;

	//second layer
  	corePV = new G4PVPlacement(piRotTiming_rear, G4ThreeVector(0, x, 6), coreLV, Form("corePV_rear_%d", iBar), worldLV, false, 0, checkOverlaps) ;
  }

  
  // fibre gap for photon counting
  // ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----  
  //end gaps
  G4VSolid* gapLayerS;
  G4VSolid* gapS;
  if (gap_size_x == -1)
  {
	gapLayerS = new G4Box  ("gapLayerS", core_radius_x, core_radius_y, 0.5*depth) ;
	gapS      = new G4Box  (     "gapS", core_radius_x, core_radius_y, 0.5*(gap_l-depth)) ;
  }
  else
  {
  	gapLayerS = new G4Box  ("gapLayerS", gap_size_x*0.5, gap_size_y*0.5, 0.5*depth) ;
  	gapS      = new G4Box  (     "gapS", gap_size_x*0.5, gap_size_y*0.5, 0.5*(gap_l-depth)) ;
  }
	
  
  G4RotationMatrix* piRot = new G4RotationMatrix;
  double theta = M_PI/4.*rad;
  double phi = 0;
  piRot->rotateZ(theta);
//  piRot->rotateX(phi);

  G4LogicalVolume* gapLayerLV = new G4LogicalVolume (gapLayerS, GaMaterial, "gapLayerLV") ;
  G4LogicalVolume* gapLV      = new G4LogicalVolume (gapS,      GaMaterial,      "gapLV") ;

  //new G4PVPlacement (0, G4ThreeVector (0., 0., 0.5*fibre_length+0.5*depth ),          gapLayerLV, "gapLayerPV", worldLV, false, 0, checkOverlaps) ;
 // new G4PVPlacement (0, G4ThreeVector (0., 0., 0.5*fibre_length+depth+0.5*(gap_l-depth)), gapLV,      "gapPV", worldLV, false, 0, checkOverlaps) ;
 


  
  // Si detector for photon counting
  // ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----
  G4VSolid * detLayerS = NULL;
  G4VSolid * detS      = NULL;
  G4LogicalVolume * detLayerLV = NULL;
  G4LogicalVolume * detLV      = NULL;
  if( detector )
  {
    if (det_size_x == -1)	//set the SiPM active area equal to the crystal surface
    {
      detLayerS = new G4Box ("detLayerS", core_radius_x*mm, core_radius_y*mm, 0.5*depth) ;
      detS      = new G4Box (     "detS", core_radius_x*mm, core_radius_y*mm, 0.5*(det_l-depth)) ;
    }
    else if (det_size_x > 0)
    {
      detLayerS = new G4Box ("detLayerS", det_size_x*0.5, det_size_y*0.5, 0.5*depth) ;
      detS      = new G4Box (     "detS", det_size_x*0.5, det_size_y*0.5, 0.5*(det_l-depth)) ;
    }

    detLayerLV = new G4LogicalVolume (detLayerS, DeMaterial, "detLayerLV") ;
    detLV      = new G4LogicalVolume (detS,      DeMaterial,      "detLV") ;
    
  //  new G4PVPlacement (0, G4ThreeVector (0., 0., 0.5*fibre_length+gap_l+0.5*depth),          detLayerLV, "detLayerPV", worldLV, false, 0, checkOverlaps) ;
  //  new G4PVPlacement (0, G4ThreeVector (0., 0., 0.5*fibre_length+gap_l+depth+0.5*(det_l-depth)), detLV,      "detPV", worldLV, false, 0, checkOverlaps) ;         
  }

   //SURFACE STATE IMPLEMENTATION

   //front face surface state (opposite to photodetector)
   G4Box* Front_skin 			= new G4Box("Front_skin", core_radius_x, core_radius_y, 0.5*depth);
   G4LogicalVolume* Front_skin_log 	= new G4LogicalVolume(Front_skin, MyMaterials::Air(), "Front_skin_log", 0, 0, 0);
//   G4VPhysicalVolume* Front_skin_phys 	= new G4PVPlacement(0, G4ThreeVector(0., 0., -(0.5*fibre_length+0.5*depth)) , Front_skin_log, "Front_skin_phys", worldLV, false, 0);

   //side surface state
   G4VSolid* dummySk = new G4Box ("dummySk", core_radius_x+depth, core_radius_y+depth, 0.5*fibre_length) ;
   G4VSolid* subSk   = new G4Box ("subSk", core_radius_x, core_radius_y, 0.51*fibre_length);
   
   G4SubtractionSolid* Side_skin 	= new G4SubtractionSolid("Side_skin", dummySk, subSk);				
   G4LogicalVolume* Side_skin_log 	= new G4LogicalVolume(Side_skin, MyMaterials::Air(), "Side_skin_log", 0, 0, 0);
  // G4VPhysicalVolume* Side_skin_phys 	= new G4PVPlacement(0, G4ThreeVector(0., 0., 0.), Side_skin_log, "Side_skin_phys", worldLV, false, 0);

   G4LogicalBorderSurface *CrystalFrontSkin   	= NULL;
   G4LogicalBorderSurface *CrystalSideSkin   	= NULL;
   G4OpticalSurface *OpCrystalSurface		= NULL;
	
   ///-------CRYSTAL SURFACE-------
   OpCrystalSurface = new G4OpticalSurface("crystal");
   OpCrystalSurface->SetType(dielectric_dielectric);
   OpCrystalSurface->SetModel(unified);




  //**********************************
  //  ELECTROMAGNETIC CALORIMETER
  //**********************************


  //double pointingAngle = 0.001976598*M_PI;	//~0.36° -- arctan (29.31/4720)
  double pointingAngle = 0.;	//~0.36° -- arctan (29.31/4720)
  double alveola_thickness = 0.2*mm;
  float total_ecal_length = ecal_front_length + ecal_rear_length;
  
  // ECAL solid
  G4Trd* ecalCrystalS_f = new G4Trd("ecalCrystalS_f",0.5*ecal_front_face, 0.5*ecal_rear_face, 0.5*ecal_front_face , 0.5*ecal_rear_face, 0.5*ecal_front_length);
  G4Trd* ecalCrystalS_r = new G4Trd("ecalCrystalS_r",0.5*ecal_front_face, 0.5*ecal_rear_face, 0.5*ecal_front_face , 0.5*ecal_rear_face, 0.5*ecal_rear_length);
//  G4Trd* ecalCrystalS_big = new G4Trd("ecalCrystalS_big",0.5*(ecal_front_face + alveola_thickness/2), 0.5*(ecal_rear_face + alveola_thickness/2), 0.5*(ecal_front_length + alveola_thickness/2) , 0.5*(ecal_rear_face+ alveola_thickness/2), 0.5*(total_ecal_length - 0.1*mm));
 // G4SubtractionSolid* ecalAlveolaS = new G4SubtractionSolid ("ecalAlveolaS", ecalCrystalS_big, ecalCrystalS);
  
  G4Box* ecalGapLayerS = new G4Box("ecalGapLayerS", ecal_det_size, ecal_det_size, 0.5*depth );
  G4Box* ecalGapS      = new G4Box("ecalGapS",      ecal_det_size, ecal_det_size, 0.5*(gap_l-depth) );
 
  G4Box* ecalDetLayerS = new G4Box("ecalDetLayerS", ecal_det_size, ecal_det_size, 0.5*depth);
  G4Box* ecalDetS      = new G4Box("ecalDetS",      ecal_det_size, ecal_det_size, 0.5*(det_l-depth));
  
  
  // ECAL logic
  G4LogicalVolume* ecalCrystalL_f = new G4LogicalVolume(ecalCrystalS_f, EcalMaterial, "ecalCrystalL_f", 0, 0, 0);
  G4LogicalVolume* ecalCrystalL_r = new G4LogicalVolume(ecalCrystalS_r, EcalMaterial, "ecalCrystalL_r", 0, 0, 0);
//  G4LogicalVolume* ecalAlveolaL = new G4LogicalVolume(ecalAlveolaS, MyMaterials::Air(), "ecalAlveolaL", 0, 0, 0);
  
  G4LogicalVolume* ecalGapLayerL = new G4LogicalVolume(ecalGapLayerS, GaMaterial, "ecalGapLayerL", 0, 0, 0);
  G4LogicalVolume* ecalGapL	 = new G4LogicalVolume(ecalGapS, GaMaterial, "ecalGapL", 0, 0, 0);
  
  G4LogicalVolume* ecalDetLayerL = new G4LogicalVolume(ecalDetLayerS, DeMaterial, "ecalDetLayerL", 0, 0, 0);
  G4LogicalVolume* ecalDetL      = new G4LogicalVolume(ecalDetS, DeMaterial, "ecalDetL", 0, 0, 0);

  
  // ECAL physical placement
  G4VPhysicalVolume* ecalCrystalP_f[25];
  G4VPhysicalVolume* ecalCrystalP_r[25];
  G4VPhysicalVolume* ecalAlveolaP[25];

  G4VPhysicalVolume* ecalGapLayerP[25];
  G4VPhysicalVolume* ecalGapP[25];
  
  G4VPhysicalVolume* ecalDetLayerP[25];
  G4VPhysicalVolume* ecalDetP[25];
  
  char name[60];

  G4double x_pos[25];
  G4double y_pos[25];

  int iCrystal;
  for (int iX = 0; iX < 5; iX ++)
  {
    for (int iY = 0; iY < 5; iY ++)
    {
      
      iCrystal = 5*iX + iY;
      x_pos[iCrystal] = iX*(ecal_rear_face+ecal_front_face + 2*alveola_thickness)/2. - ecal_front_face*2.5 + ecal_front_face/2;	// position the baricenter of crystals and then rotating them by
      y_pos[iCrystal] = iY*(ecal_rear_face+ecal_front_face + 2*alveola_thickness)/2. - ecal_front_face*2.5 + ecal_front_face/2;
      
      G4RotationMatrix* piRotEcal = new G4RotationMatrix;
      piRotEcal->rotateY(-pointingAngle*rad*(iX + 0.5));
      piRotEcal->rotateX(pointingAngle*rad*(iY + 0.5));

      cout << " x_pos [" <<iCrystal << "] = " << x_pos[iCrystal] << " :: y_pos[" << iCrystal << "] = " << y_pos[iCrystal] << " :: angle = [" << pointingAngle*iX << ", " << pointingAngle*iY << "] " << endl;
      
      sprintf(name, "ecalCrystalP_f_%d", iCrystal);
      ecalCrystalP_f[iCrystal] = new G4PVPlacement(piRotEcal, G4ThreeVector(x_pos[iCrystal], y_pos[iCrystal], ecal_timing_distance+ecal_front_length*0.5), ecalCrystalL_f, name, worldLV, false, 0);

      sprintf(name, "ecalCrystalP_r_%d", iCrystal);
      ecalCrystalP_r[iCrystal] = new G4PVPlacement(piRotEcal, G4ThreeVector(x_pos[iCrystal], y_pos[iCrystal], ecal_timing_distance + ecal_front_length + ecal_rear_length*0.5), ecalCrystalL_r, name, worldLV, false, 0);
  /*    
      sprintf(name, "ecalAlveolaP_%d", iCrystal);
      ecalAlveolaP[iCrystal] = new G4PVPlacement(piRotEcal, G4ThreeVector(x_pos[iCrystal], y_pos[iCrystal], ecal_timing_distance), ecalAlveolaL, name, worldLV, false, 0);

      sprintf(name,"ecalGapLayerP_%d",iCrystal);
      ecalGapLayerP[iCrystal] = new G4PVPlacement(piRotEcal, G4ThreeVector(x_pos_pmt[iCrystal], y_pos_pmt[iCrystal], (total_ecal_length + gap_l + ecal_timing_distance)/2), ecalGapLayerL, name, worldLV, false, 0);
      
      sprintf(name,"ecalGapP_%d",iCrystal);
      ecalGapP[iCrystal] = new G4PVPlacement(piRotEcal, G4ThreeVector(x_pos_pmt[iCrystal], y_pos_pmt[iCrystal], (total_ecal_length + win_l)*0.5 + gap_l + ecal_timing_distance), ecalGapL, name, worldLV, false, 0);

      sprintf(name,"ecalDetLayerP_%d", iCrystal);
      ecalDetLayerP[iCrystal] = new G4PVPlacement(piRotEcal, G4ThreeVector(x_pos_pmt[iCrystal], y_pos_pmt[iCrystal], (total_ecal_length + depth)*0.5 + gap_l + det_l + ecal_timing_distance), ecalDetLayerL, name, worldLV, false, 0);

      sprintf(name,"ecalDetP_%d", iCrystal);
      ecalDetP[iCrystal] = new G4PVPlacement(piRotEcal, G4ThreeVector(x_pos_pmt[iCrystal], y_pos_pmt[iCrystal], total_ecal_length*0.5 + gap_l + det_l + depth + ecal_timing_distance) , ecalDetL, name, worldLV, false, 0);
*/
    }
  }
  
  
  
  //-----------------------------------------------------
  //------------- Visualization attributes --------------
  //-----------------------------------------------------
  
  G4Colour white  (1.00, 1.00, 1.00);  // white
  G4Colour gray   (0.50, 0.50, 0.50);  // gray
  G4Colour black  (0.00, 0.00, 0.00);  // black
  G4Colour red    (1.00, 0.00, 0.00);  // red
  G4Colour green  (0.00, 1.00, 0.00);  // green
  G4Colour blue   (0.00, 0.00, 1.00);  // blue
  G4Colour cyan   (0.00, 1.00, 1.00);  // cyan
  G4Colour air    (0.90, 0.94, 1.00);  // cyan
  G4Colour magenta(1.00, 0.00, 1.00);  // magenta 
  G4Colour yellow (1.00, 1.00, 0.00);  // yellow
  G4Colour brass  (0.80, 0.60, 0.40);  // brass
  G4Colour brown  (0.70, 0.40, 0.10);  // brown
  
  G4VisAttributes* VisAttWorld = new G4VisAttributes(black);
  VisAttWorld -> SetVisibility(true) ;
  VisAttWorld -> SetForceWireframe(true) ;
  worldLV -> SetVisAttributes(VisAttWorld) ;
  
  G4VisAttributes* VisAttCore = new G4VisAttributes(green);
  VisAttCore -> SetVisibility(true);
  VisAttCore -> SetForceWireframe(true);
  coreLV -> SetVisAttributes(VisAttCore);

  G4VisAttributes* VisCrystalCore = new G4VisAttributes(cyan);
  VisCrystalCore -> SetVisibility(true);
  VisCrystalCore -> SetForceWireframe(true);
  ecalCrystalL_f -> SetVisAttributes(VisCrystalCore);

  G4VisAttributes* VisCrystalCore_r = new G4VisAttributes(blue);
  VisCrystalCore_r -> SetVisibility(true);
  VisCrystalCore_r -> SetForceWireframe(true);
  ecalCrystalL_r -> SetVisAttributes(VisCrystalCore_r);
  
  G4VisAttributes* VisAttGapLayer = new G4VisAttributes(gray);
  VisAttGapLayer -> SetVisibility(true);
  VisAttGapLayer -> SetForceWireframe(true);
  gapLayerLV -> SetVisAttributes(VisAttGapLayer);

  
  G4VisAttributes* VisAttGap = new G4VisAttributes(blue);
  VisAttGap -> SetVisibility(true);
  VisAttGap -> SetForceWireframe(true);
  gapLV -> SetVisAttributes(VisAttGap);

  
  
  if( detector )
  {
    G4VisAttributes* VisAttDetLayer = new G4VisAttributes(red);
    VisAttDetLayer -> SetVisibility(true);
    VisAttDetLayer -> SetForceWireframe(false);
    detLayerLV -> SetVisAttributes(VisAttDetLayer);
    
    G4VisAttributes* VisAttDet = new G4VisAttributes(gray);
    VisAttDet -> SetVisibility(true);
    VisAttDet -> SetForceWireframe(false);
    detLV -> SetVisAttributes(VisAttDet);
  }
  
  if (B_field_intensity > 0.1 * tesla) ConstructField () ; 
  
  
  
  // //-----------------------------------------------
  // //------------- Fast photon timing --------------
  // //-----------------------------------------------
  
  // std::vector<std::pair<double,double> > rIndVecCore;
  // std::vector<std::pair<double,double> > rIndVecClad;
  // std::vector<std::pair<double,double> > rIndVecAir;
  // std::vector<std::pair<double,double> > rIndVecGap;
  
  // G4MaterialPropertyVector* mpVec;
  
  // mpVec = ClMaterial->GetMaterialPropertiesTable()->GetProperty("RINDEX");
  // for(unsigned int it = 0; it < mpVec->GetVectorLength(); ++it)
  // {
  //   std::pair<double,double> dummy(mpVec->GetLowEdgeEnergy(it)/eV,(*mpVec)[it]);
  //   rIndVecCore.push_back(dummy);
  // }
  
  // mpVec = WoMaterial->GetMaterialPropertiesTable()->GetProperty("RINDEX");
  // for(unsigned int it = 0; it < mpVec->GetVectorLength(); ++it)
  // {
  //   std::pair<double,double> dummy(mpVec->GetLowEdgeEnergy(it)/eV,(*mpVec)[it]);
  //   std::pair<double,double> dummy2(mpVec->GetLowEdgeEnergy(it)/eV,fibre_cladRIndex);
  //   rIndVecAir.push_back(dummy);
  //   rIndVecClad.push_back(dummy2);
  // }
  
  // mpVec = GaMaterial->GetMaterialPropertiesTable()->GetProperty("RINDEX");
  // for(unsigned int it = 0; it < mpVec->GetVectorLength(); ++it)
  // {
  //   std::pair<double,double> dummy(mpVec->GetLowEdgeEnergy(it)/eV,(*mpVec)[it]);
  //   rIndVecGap.push_back(dummy);
  // }
  
  
  // fib = FiberInit(fibre_length,fibre_radius,CreateTree::Instance()->attLengths,rIndVecCore,rIndVecClad,rIndVecAir,rIndVecGap) ;
  
  
  
  G4cout << ">>>>>> DetectorConstruction::Construct ()::end <<< " << G4endl ;
  return worldPV ;
}

//---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- 



void DetectorConstruction::initializeMaterials ()
{
  //-----------------
  // define materials
  
  WoMaterial = NULL ;
  if      ( world_material == 1 ) WoMaterial = MyMaterials::Air () ;
  else
  {
    G4cerr << "<DetectorConstructioninitializeMaterials>: Invalid fibre world material specifier " << world_material << G4endl ;
    exit (-1) ;
  }
  G4cout << "Wo. material: "<< WoMaterial << G4endl ;
  
  
  CoMaterial = NULL ;
  if      ( core_material == 1 )  CoMaterial = MyMaterials::Quartz();
  else if ( core_material == 2 )  CoMaterial = MyMaterials::SiO2();
  else if ( core_material == 3 )  CoMaterial = MyMaterials::SiO2_Ce();
  else if ( core_material == 4 )  CoMaterial = MyMaterials::LuAG_Ce();
  else if ( core_material == 5 )  CoMaterial = MyMaterials::YAG_Ce();
  else if ( core_material == 6 )  CoMaterial = MyMaterials::LSO();
  else if ( core_material == 7 )  CoMaterial = MyMaterials::LYSO();
  else if ( core_material == 8 )  CoMaterial = MyMaterials::LuAG_undoped();
  else if ( core_material == 9 )  CoMaterial = MyMaterials::GAGG_Ce();
  else if ( core_material == 11 ) CoMaterial = MyMaterials::LuAG_Pr();
  else if ( core_material == 12 ) CoMaterial = MyMaterials::PbF2();
  else if ( core_material == 13 ) CoMaterial = MyMaterials::PlasticBC408();
  else if ( core_material == 14 ) CoMaterial = MyMaterials::PlasticBC418();
  else if ( core_material == 15)  CoMaterial = MyMaterials::ToyLSO(toy_ly, toy_decay, toy_rise);
  else
  {
    G4cerr << "<DetectorConstructioninitializeMaterials>: Invalid fibre clad material specifier " << core_material << G4endl ;
    exit (-1) ;
  }
  G4cout << "Co. material: "<< CoMaterial << G4endl ;

  EcalMaterial = NULL ;
  if      ( ecal_material == 1 )  EcalMaterial = MyMaterials::Quartz();
  else if ( ecal_material == 2 )  EcalMaterial = MyMaterials::SiO2();
  else if ( ecal_material == 3 )  EcalMaterial = MyMaterials::SiO2_Ce();
  else if ( ecal_material == 4 )  EcalMaterial = MyMaterials::LuAG_Ce();
  else if ( ecal_material == 5 )  EcalMaterial = MyMaterials::YAG_Ce();
  else if ( ecal_material == 6 )  EcalMaterial = MyMaterials::LSO();
  else if ( ecal_material == 7 )  EcalMaterial = MyMaterials::LYSO();
  else if ( ecal_material == 8 )  EcalMaterial = MyMaterials::LuAG_undoped();
  else if ( ecal_material == 9 )  EcalMaterial = MyMaterials::GAGG_Ce();
  else if ( ecal_material == 10 ) EcalMaterial = MyMaterials::LuAG_Pr();
  else if ( ecal_material == 11 ) EcalMaterial = MyMaterials::PbF2();
  else if ( ecal_material == 12 ) EcalMaterial = MyMaterials::PlasticBC408();
  else if ( ecal_material == 13 ) EcalMaterial = MyMaterials::PlasticBC418();
  else if ( ecal_material == 14 ) EcalMaterial = MyMaterials::PWO();
  else if ( ecal_material == 15)  EcalMaterial = MyMaterials::ToyLSO(toy_ly, toy_decay, toy_rise);
  else
  {
    G4cerr << "<DetectorConstructioninitializeMaterials>: Invalid fibre clad material specifier " << ecal_material << G4endl ;
    exit (-1) ;
  }
  G4cout << "ECAL material: "<< EcalMaterial << G4endl ;
  
  
  GaMaterial = NULL;
  if     ( gap_material == 1 ) GaMaterial = MyMaterials::Air();
  else if( gap_material == 2 ) GaMaterial = MyMaterials::OpticalGrease();
  else if( gap_material == 3 ) GaMaterial = MyMaterials::MeltMount168();
  else if( gap_material == 4 ) GaMaterial = MyMaterials::OpticalGrease155();
  else
  {
    G4cerr << "<DetectorConstructioninitializeMaterials>: Invalid gap material specifier " << gap_material << G4endl;
    exit(-1);
  }
  G4cout << "Gap material: " << gap_material << G4endl;
  
  
  DeMaterial = NULL;
  if	 ( det_material == 1 ) DeMaterial = MyMaterials::Silicon();
  else if( det_material == 2 ) DeMaterial = MyMaterials::Quartz();
  else if( det_material == 3 ) DeMaterial = MyMaterials::Air();
  else
  {
    G4cerr << "<DetectorConstructioninitializeMaterials>: Invalid detector material specifier " << det_material << G4endl;
    exit(-1);
  }
  G4cout << "Detector material: " << det_material << G4endl;
  
  
  
  //------------------
  // change properties
  
  if( core_absLength > 0 )
  {
    const G4int nEntries_ABS = 2;
    G4double PhotonEnergy_ABS[nEntries_ABS] = { 1.*eV, 10.*eV };
    G4double Absorption[nEntries_ABS] = { core_absLength*mm, core_absLength*mm };
    
    CoMaterial -> GetMaterialPropertiesTable() -> RemoveProperty("ABSLENGTH");
    CoMaterial -> GetMaterialPropertiesTable() -> AddProperty("ABSLENGTH",PhotonEnergy_ABS,Absorption,nEntries_ABS);
  }
  if( core_rIndex > 0 )
  {
    const G4int nEntries_RI = 2;
    G4double PhotonEnergy_RI[nEntries_RI] = { 1.*eV, 10.*eV };
    G4double RefractiveIndex[nEntries_RI] = { core_rIndex, core_rIndex };
    
    CoMaterial -> GetMaterialPropertiesTable() -> RemoveProperty("RINDEX");
    CoMaterial -> GetMaterialPropertiesTable() -> AddProperty("RINDEX",PhotonEnergy_RI,RefractiveIndex,nEntries_RI);
  }
  
}

//---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- 



void DetectorConstruction::ConstructField () 
{
  G4cout << ">>>>>> DetectorConstruction::ConstructField ()::begin <<<<<<" << G4endl ;
  
  static G4TransportationManager * trMgr = G4TransportationManager::GetTransportationManager () ; 
  
  // A field object is held by a field manager
  // Find the global Field Manager
  G4FieldManager * globalFieldMgr = trMgr->GetFieldManager () ;
  
  if( !B_field_IsInitialized )
  {
    // magnetic field parallel to the beam direction (w/ tilt)
    G4ThreeVector fieldVector(0.0522*B_field_intensity,0.0522*B_field_intensity,0.9973*B_field_intensity);   
    
    B_field = new G4UniformMagField (fieldVector) ; 
    globalFieldMgr->SetDetectorField (B_field) ;
    globalFieldMgr->CreateChordFinder (B_field) ;
    globalFieldMgr->GetChordFinder ()->SetDeltaChord (0.005 * mm) ;
    B_field_IsInitialized = true ;
  }
  
  G4cout << ">>>>>> DetectorConstruction::ConstructField ()::end <<< " << G4endl ;
  return ;
}

void DetectorConstruction::SetMaxStep(G4double maxStep)
{
  if ((stepLimit)&&(maxStep>0.)) stepLimit->SetMaxAllowedStep(maxStep);
}


// Initialization classes
//
/*
void DetectorConstruction::initializeSurface(G4OpticalSurface *mySurface, string surfaceType)
{
    if(surfaceType == "crystal") 
    {
//         cout << "CRISTALLO " << crystalSurfinish << endl;

        surfinish   	= crystalSurfinish;
        RefFile     	= cReffile;
        reflectivity    = cReflectivity;
        surrefind   	= cSurrefind;
        surtype     	= cSurtype;
        specularspike   = cSpecularspike;
        specularlobe    = cSpecularlobe;
        sigmaalpha  	= cSigmaalpha;
        backscatter     = cBackscatter;
        lambertian  	= cLambertian;
    }



    if(this->surfinish <= 5)
    {
        G4cout << "Using unified model." << G4endl;
        mySurface -> SetModel(unified);
        switch(this->surtype) {
        case 0:
            mySurface -> SetType(dielectric_metal);
            G4cout << "Surface type: dielectric_metal" << G4endl;
            break;
        case 1:
            mySurface -> SetType(dielectric_dielectric);
            G4cout << "Surface type: dielectric_dielectric" << G4endl;
            break;
        }
    }

    else if(this->surfinish > 5 && surfaceType == "wrapping") G4cout << "Value not allowed" << G4endl;

    else {
        G4cout << "Using LUT for surface treatment." << G4endl;
        mySurface -> SetModel(LUT);
        mySurface -> SetType(dielectric_LUT);
    }

    switch(this->surfinish) {
    case 0:
        mySurface -> SetFinish(polished);
        G4cout << "Surface finish: polished" << G4endl;
        break;
    case 1:
        mySurface -> SetFinish(polishedfrontpainted);
        G4cout << "Surface finish: polishedfrontpainted" << G4endl;
        break;
    case 2:
        mySurface -> SetFinish(polishedbackpainted);
        G4cout << "Surface finish: polishedbackpainted" << G4endl;
        break;
    case 3:
        mySurface -> SetFinish(ground);
        G4cout << "Surface finish: ground" << G4endl;
        break;
    case 4:
        mySurface -> SetFinish(groundfrontpainted);
        G4cout << "Surface finish: groundfrontpainted" << G4endl;
        break;
    case 5:
        mySurface -> SetFinish(groundbackpainted);
        G4cout << "Surface finish: groundbackpainted" << G4endl;
        break;
    case 17:
        mySurface -> SetFinish(polishedteflonair);
        G4cout << "Surface finish: polishedteflonair" << G4endl;
        break;
    case 18:
        mySurface -> SetFinish(polishedtioair);
        G4cout << "Surface finish: polishedtioair" << G4endl;
        break;
    case 26:
        mySurface -> SetFinish(etchedtioair);
        G4cout << "Surface finish: etchedtioair" << G4endl;
        break;
    case 34:
        mySurface -> SetFinish(groundtioair);
        G4cout << "Surface finish: groundtioair" << G4endl;
        break;
    case 36:
        mySurface -> SetFinish(polishedtyvekair);
        G4cout << "Surface finish: polishedtyvekair" << G4endl;
        break;
    default:
        G4cout << "Surface finish unkown!" << G4endl;
        exit(0);
    }
}*/


//
// reflectivity
//
/*
void DetectorConstruction::initializeReflectivitySurface(G4OpticalSurface *surface, string surfaceType)
{

    int NumRefl = 0;
    G4double EphotonRefl[1000];
    G4double Refl[1000];
    
    EphotonRefl[0] = 0.0001 * eV;
    EphotonRefl[1] = 1.0 * eV;
    EphotonRefl[2] = 4.08 * eV;
    Refl[0] = 0.0; // suppress photons with energy < 1eV (will not be detected)
    Refl[1] = this->crystal_reflectivity;
    Refl[2] = this->crystal_reflectivity;
    NumRefl = 100;

    G4cout << "Reflectivities as a function of the photon energy:" << G4endl;
    for(int i = 0; i < NumRefl; i++)         G4cout << i << "   " << EphotonRefl[i] << "   " << Refl[i] << G4endl;
    

    Ephoton[0] = 0.0001 * eV;
    Ephoton[1] = 1.0 * eV;
    Ephoton[2] = 4.08 * eV;
    G4double RefractiveIndex[3] = {this->surrefind, this->surrefind, this->surrefind};
    G4double SpecularLobe[3]    = {this->specularlobe, this->specularlobe, this->specularlobe};
    G4double SpecularSpike[3]   = {this->specularspike, this->specularspike, this->specularspike};
    G4double Backscatter[3]     = {this->backscatter, this->backscatter, this->backscatter};
    G4double Lambertian[3]      = {this->lambertian, this->lambertian, this->lambertian};
    G4MaterialPropertiesTable *myST = new G4MaterialPropertiesTable();
    G4cout << "Read from config-file: " << G4endl;
    G4cout << "Read SPECULARLOBECONSTANT  : " << SpecularLobe[0] << G4endl;
    G4cout << "Read SPECULARSPIKECONSTANT : " << SpecularSpike[0] << G4endl;
    G4cout << "Read BACKSCATTERCONSTANT   : " << Backscatter[0] << G4endl;
    G4cout << "Read LAMBERTIAN            : " << Lambertian[0] << G4endl;
    G4cout << "Read ref. index            : " << RefractiveIndex[0] << G4endl;

    myST->AddProperty("RINDEX", Ephoton, RefractiveIndex, 3);
    if(this->specularlobe >= 0) {
        G4cout << "Setting SPECULARLOBECONSTANT to : " << SpecularLobe[0] << G4endl;
        myST->AddProperty("SPECULARLOBECONSTANT",  Ephoton, SpecularLobe,    3);
    }
    if(this->specularspike >= 0) {
        G4cout << "Setting SPECULARSPIKECONSTANT to : " << SpecularSpike[0] << G4endl;
        myST->AddProperty("SPECULARSPIKECONSTANT", Ephoton, SpecularSpike,   3);
    }
    if(this->backscatter >= 0) {
        G4cout << "Setting BACKSCATTERCONSTANT to : " << Backscatter[0] << G4endl;
        myST->AddProperty("BACKSCATTERCONSTANT",   Ephoton, Backscatter,     3);
    }
    if(this->lambertian >= 0) {
        G4cout << "Setting LAMBERTIAN to : " << Lambertian[0] << G4endl;
        myST->AddProperty("LAMBERTIAN",            Ephoton, Lambertian,      3);
    }

//     setWrappingIndices(myST,1.82,"teflon");

    surface->SetMaterialPropertiesTable(myST);
    if(this->sigmaalpha >= 0) surface->SetSigmaAlpha(this->sigmaalpha);


}*/

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

//---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- 

