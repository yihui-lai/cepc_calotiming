#ifndef CreateTree_H
#define CreateTree_H 1

#include <iostream>
#include <vector>
#include "TString.h"

#include "TH2F.h"
#include "TH3F.h"
#include "TFile.h"
#include "TTree.h"
#include "TNtuple.h"

class CreateTree
{
private:
  
  TTree*  ftree ;
  TString fname ;
  
public:
  
  CreateTree (TString name);
  ~CreateTree () ;
  
  TTree*             GetTree() const { return ftree; };
  TString            GetName() const { return fname; };
  void               AddEnergyDeposit(int index, float deposit);
  void               AddScintillationPhoton(int index);
  void               AddCerenkovPhoton(int index);
  int                Fill();
  bool               Write(TFile *);
  void               Clear() ;
  
  static CreateTree* Instance() { return fInstance; } ;
  static CreateTree* fInstance;
  
  int Event;
  
  std::vector<float>* inputMomentum ; // Px Py Pz E
  std::vector<float>* inputInitialPosition ; // x, y, z

  //integrated energy in each longitudinal layer
  float depositedEnergyTotal;

  float depositedEnergyAbsorber;

  float depositedEnergyTiming_f;
  float depositedEnergyTiming_r;

  float depositedEnergyECAL_f;
  float depositedEnergyECAL_r;

  float depositedEnergyWorld;


  //energy deposit in each trasnversally segmented channel
  float Edep_Timing_f_ch[18];
  float Edep_Timing_r_ch[18];

  float Edep_ECAL_f_ch[25];
  float Edep_ECAL_r_ch[25];


  /*
  std::vector<float> E_dep_f;
  std::vector<float> E_dep_time_f;
  std::vector<float> E_dep_r;
  std::vector<float> E_dep_time_r;
  
  std::vector<float> time_ext_scint;
  std::vector<float> time_ext_cher;
  std::vector<float> time_prod_scint;
  std::vector<float> time_prod_cher;

  std::vector<float> lambda_ext_scint;
  std::vector<float> lambda_ext_cher;
  std::vector<float> lambda_prod_scint;
  std::vector<float> lambda_prod_cher;
  
  std::vector<float> time_ext_scint_ref;
  std::vector<float> time_ext_cher_ref;
  std::vector<float> time_prod_scint_ref;
  std::vector<float> time_prod_cher_ref;
  */


  /*


  int tot_phot_sci;
  int tot_phot_cer;
  int tot_latGap_phot_sci;
  int tot_latGap_phot_cer;
  int tot_gap_phot_sci;
  int tot_gap_phot_cer;
  int tot_det_phot_sci;
  int tot_det_phot_cer;
  
  TH1F* h_phot_sci_lambda;
  TH1F* h_phot_sci_E;
  TH1F* h_phot_sci_time;
  TH1F* h_phot_sci_angleAtProduction;
  TH1F* h_phot_cer_lambda;
  TH1F* h_phot_cer_E;
  TH1F* h_phot_cer_time;
  TH1F* h_phot_cer_angleAtProduction;
  
  TH1F* h_phot_sci_gap_lambda;
  TH1F* h_phot_sci_gap_E;
  TH1F* h_phot_sci_gap_time;
  TH1F* h_phot_sci_gap_angleAtProduction;
  TH1F* h_phot_sci_gap_angleWithSurfNormal;
  TH1F* h_phot_cer_gap_lambda;
  TH1F* h_phot_cer_gap_E;
  TH1F* h_phot_cer_gap_time;
  TH1F* h_phot_cer_gap_angleAtProduction;
  TH1F* h_phot_cer_gap_angleWithSurfNormal;

  TH1F* h_phot_sci_det_lambda;
  TH1F* h_phot_sci_det_E;
  TH1F* h_phot_sci_det_time;
  TH1F* h_phot_sci_det_angleAtProduction;
  TH1F* h_phot_sci_det_angleWithSurfNormal;
  TH1F* h_phot_cer_det_lambda;
  TH1F* h_phot_cer_det_E;
  TH1F* h_phot_cer_det_time;
  TH1F* h_phot_cer_det_angleAtProduction;
  TH1F* h_phot_cer_det_angleWithSurfNormal;
  */

};

#endif
