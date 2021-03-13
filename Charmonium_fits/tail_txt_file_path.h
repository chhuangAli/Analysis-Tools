#include "TString.h"

  TString path_tails_diretory = Form("/home/chunlu/cernbox/AnaSourceCode/anaJPsi_2015_2018PbPb/MC_analysis/determine_signal_tail_parameters/Tails_2015_2018");
//  TString path_tails_diretory = Form("/home/chunlu/cernbox/AnaSourceCode/anaJPsi_2015_2018PbPb/MC_analysis/determine_signal_tail_parameters/Tails_2015_2018/y_dep");

//  TString path_geant4_CB2_tail00 = Form("%s/Geant4_CB2_signal_tail_2.0to4.8.txt", path_tails_diretory.Data() );
//  TString path_geant4_CB2_tail01 = Form("%sGeant4_CB2_signal_tail_2.2to4.4.txt");
  TString path_geant3_CB2_tail00 = Form("%s/Geant3_CB2_signal_tail_2.2to4.5.txt", path_tails_diretory.Data() ); 
  TString path_geant3_CB2_tail01 = Form("%s/Geant3_CB2_signal_tail_2.4to4.7.txt", path_tails_diretory.Data() );
//  TString path_geant4_NA60_tail00 = Form("/home/chunlu/cernbox/AnaSourceCode/anaJPsi_2017pp/MC_analysis/signal_extract/Geant4_plots_and_values/Geant4_NA60_signal_tail_2.0to4.8.txt");
//  TString path_geant4_NA60_tail01 = Form("/home/chunlu/cernbox/AnaSourceCode/anaJPsi_2017pp/MC_analysis/signal_extract/Geant4_plots_and_values/Geant4_NA60_signal_tail_2.2to4.4.txt");
  TString path_geant3_NA60_tail00 = Form("%s/Geant3_NA60_signal_tail_2.2to4.5.txt", path_tails_diretory.Data());
  TString path_geant3_NA60_tail01 = Form("%s/Geant3_NA60_signal_tail_2.4to4.7.txt", path_tails_diretory.Data());

  int num_pT_bins = 17;

//  int num_pT_bins = 7; // attension 7 is the number of rapidity bins
//  int pt_bin =0;
//  double** values_geant4_CB2_tail_range00 = new double*[num_pT_bins];
//  double** values_geant4_CB2_tail_range01 = new double*[num_pT_bins];
  double** values_geant3_CB2_tail_range00 = new double*[num_pT_bins];
  double** values_geant3_CB2_tail_range01 = new double*[num_pT_bins];

//  double** values_geant4_NA60_tail_range00 = new double*[num_pT_bins];
//  double** values_geant4_NA60_tail_range01 = new double*[num_pT_bins];
  double** values_geant3_NA60_tail_range00 = new double*[num_pT_bins];
  double** values_geant3_NA60_tail_range01 = new double*[num_pT_bins];
