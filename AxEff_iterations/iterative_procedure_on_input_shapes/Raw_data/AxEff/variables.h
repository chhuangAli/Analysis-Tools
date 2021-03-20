#include "TString.h"
#include "TH1F.h"

  const int NUMOFFILES_LHC18q=130;
//  TString runList_location = Form("/home/chunlu/local/runList_LHC18q/muon_calo_pass2/LHC18q_runList.txt");
  const int NUMOFFILES_LHC18r=99;
//  const int NUMOFFILES = 228;
  const int NUMOFFILES = 229;
// TString LHC18r_runList_location = Form("/home/chunlu/local/runList_LHC18r/muon_calo_pass2/LHC18r_runList.txt");

  TString runList_location = Form("/home/chunlu/local/runList_LHC18q_r/muon_calo_pass2/LHC18q_r_runList.txt");
//  TString runList_location = Form("/home/chunlu/local/runList_LHC18q_r/muon_calo_pass2/LHC18q_r_runList_without_296752.txt");

  
  int run_number[NUMOFFILES_LHC18q] = {0};
  TString tstr_run_number[NUMOFFILES_LHC18q]={};

  int LHC18r_run_number[NUMOFFILES_LHC18r] = {0};
  TString tstr_LHC18r_run_number[NUMOFFILES_LHC18r]={};

  int LHC18q_r_run_number[NUMOFFILES] = {0};
  TString tstr_LHC18q_r_run_number[NUMOFFILES]={};  

  TString centralityBins_alicounter[] = {"m0","0_10","10_20","20_30","30_40","40_50","50_60","60_70","70_80","80_90","90_100"};
  const int MAXINDEX_CENTBINS = (int)(sizeof(centralityBins_alicounter)/sizeof(centralityBins_alicounter[0]));

  TString anaResult_file[NUMOFFILES];
  TH1F *histo_Ae_run; 

