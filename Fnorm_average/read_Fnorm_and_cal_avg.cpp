////// ********************************
// This macro is to read Fnorm obtained by other analyzers, and calculate the average Fnorm by taking weights of different CMUL events. 
////// ********************************


#include <iostream>
#include <fstream>
#include <cmath>

// ROOT includes
#include <TFile.h>
#include <TCanvas.h>
#include <TH1F.h>
#include <TH1D.h>
#include <TObjArray.h>
#include "TObjString.h"
#include "TString.h"
#include "TPaveText.h"
#include "TPaveStats.h"
#include "TLegend.h"

#include "TFitResult.h"
#include "TFitResultPtr.h"
#include "TF1.h"
#include "TPad.h"
#include "TMatrixDSym.h"
#include "TAttMarker.h"
#include "TLatex.h"
#include "TGaxis.h"
#include "TAxis.h"

//PWG includes
//#include "AliCounterCollection.h"

#include "MathTechnicalFunc.h"

double read_alicounter_obtain_events(TString, const int, TString, int *);
double read_Fnorm_cal_avg(TString, TString , TString , const int, TString, double*);
double read_Fnorm_cal_ave_err(TString *path_to_Fnorm, TString *name_Fnorm_histo, TString *analysis_run_files_location, const int *number_of_runs, TString *tstr_alicounter_path, double *Fnorm_fraction);
double mean_with_stat_weight(double *, double* , int);


int main()
{

//************************************
//
// The calculation of Fnorm for LHC15o
//
//************************************

  const int number_methods = 3;

  TString path_to_directory = "$HOME/local/anaedPbPbData2015/Fnorm/";
  TString path_to_file = "FNormRunByRun_15o.root";

  TString path_to_Fnorm[3];
  path_to_Fnorm[0] = Form("%s%s", path_to_directory.Data(), path_to_file.Data() );
  TString name_Fnorm_histo[3];
  name_Fnorm_histo[0] = "fFNormDirect_CINT7"; // This variable tells the Fnorm obtained from a given method. 

  TString analysis_run_files_location = Form("/home/chunlu/local/runList_LHC15o/runList.txt");
  const int number_of_runs = 137; 

  TString tstr_alicounter_path = "$HOME/local/anaedPbPbData2015/LHC15o/muon_calo_pass1/AOD229_CMULEvent_AnaResults/";


  double avg_Fnorm_LHC15o[number_methods] = {0}; 
  double avg_Fnorm_err_LHC15o[number_methods] = {0};
  double final_Fnorm_LHC15o = 0; 
  double final_Fnorm_err_LHC15o = 0;
  double Fnorm_fraction_LHC15o[3]={0}; 

//****************************
//
// Calculate the average of Fnorm by integrating over runs. This Fnorm is obtained from offline direct method.
//
//****************************

  avg_Fnorm_LHC15o[0] = read_Fnorm_cal_avg(path_to_Fnorm[0], name_Fnorm_histo[0], analysis_run_files_location, number_of_runs, tstr_alicounter_path, Fnorm_fraction_LHC15o );
  avg_Fnorm_err_LHC15o[0] = sqrt(Fnorm_fraction_LHC15o[2]) / Fnorm_fraction_LHC15o[1];

  std::cout << Fnorm_fraction_LHC15o[0] << " / " << Fnorm_fraction_LHC15o[1] << "+/-" << sqrt(Fnorm_fraction_LHC15o[2]) / Fnorm_fraction_LHC15o[1]  << std::endl;

//****************************
//
// Calculate the average of Fnorm by integrating over runs. This Fnorm is obtained from offline indirect method.
//
//****************************
  name_Fnorm_histo[0] = "fFNormIndirect_CINT7";
  avg_Fnorm_LHC15o[1] = read_Fnorm_cal_avg(path_to_Fnorm[0], name_Fnorm_histo[0], analysis_run_files_location, number_of_runs, tstr_alicounter_path, Fnorm_fraction_LHC15o );
  avg_Fnorm_err_LHC15o[1] = sqrt(Fnorm_fraction_LHC15o[2]) / Fnorm_fraction_LHC15o[1];
  std::cout << Fnorm_fraction_LHC15o[0] << " / " << Fnorm_fraction_LHC15o[1] << "+/-" << avg_Fnorm_err_LHC15o[1]   << std::endl;

//****************************
//
// Calculate an average of Fnorm by integrating over runs. This Fnorm is obtained from online method.
//
//****************************
  name_Fnorm_histo[0] = "fFNormOnline_C0V0M";
  avg_Fnorm_LHC15o[2] = read_Fnorm_cal_avg(path_to_Fnorm[0], name_Fnorm_histo[0], analysis_run_files_location, number_of_runs, tstr_alicounter_path, Fnorm_fraction_LHC15o );
  avg_Fnorm_err_LHC15o[2] = sqrt(Fnorm_fraction_LHC15o[2]) / Fnorm_fraction_LHC15o[1];
  std::cout << Fnorm_fraction_LHC15o[0] << " / " << Fnorm_fraction_LHC15o[1] << "+/-" << avg_Fnorm_err_LHC15o[2]  << std::endl;

//****************************
//
// Calculate an average of the three Fnorm obtained above. 
//
//****************************
  final_Fnorm_LHC15o = mean_with_stat_weight(avg_Fnorm_LHC15o, avg_Fnorm_err_LHC15o, number_methods);
  final_Fnorm_err_LHC15o = sqrt(1/(pow(1/avg_Fnorm_err_LHC15o[0],2) + pow(1/avg_Fnorm_err_LHC15o[1],2) + pow(1/avg_Fnorm_err_LHC15o[2],2)));

  std::cout << final_Fnorm_LHC15o << "+/-" << final_Fnorm_err_LHC15o << std::endl;

//******************************************
//
// The calculation of Fnorm for LHC18q. Using the same formulas which are shown above.
//
//******************************************
  double avg_Fnorm_LHC18q[number_methods] = {0};
  double avg_Fnorm_err_LHC18q[number_methods] = {0};
  double final_Fnorm_LHC18q = 0;
  double final_Fnorm_err_LHC18q = 0;
 
  path_to_directory = "$HOME/local/anaedPbPbData2018/Fnorm/";
  path_to_file = "FNormRunByRun_18q.root";
  path_to_Fnorm[1] = Form("%s%s", path_to_directory.Data(), path_to_file.Data() );

//****************************************************
//
// For Fnorm obtained from offline direct method
//
//****************************************************

  name_Fnorm_histo[1] = "fFNormDirect_CINT7";

//  analysis_run_files_location = Form("/home/chunlu/local/runList_LHC18q/muon_calo_pass2/LHC18q_runList.txt");
  analysis_run_files_location = Form("/home/chunlu/local/runList_LHC18q_r/muon_calo_pass3/LHC18q_runList.txt");
  const int number_of_runs_LHC18q = 130;

  tstr_alicounter_path = "$HOME/local/anaedPbPbData2018/LHC18q_r/muon_calo_pass3/AOD225_CMULEvent_AnaResults/";

  double Fnorm_fraction_LHC18q[3]={0};

  avg_Fnorm_LHC18q[0] = read_Fnorm_cal_avg(path_to_Fnorm[1], name_Fnorm_histo[1], analysis_run_files_location, number_of_runs_LHC18q, tstr_alicounter_path, Fnorm_fraction_LHC18q );
  avg_Fnorm_err_LHC18q[0] = sqrt(Fnorm_fraction_LHC18q[2]) / Fnorm_fraction_LHC18q[1];
  std::cout << Fnorm_fraction_LHC18q[0] << " / " << Fnorm_fraction_LHC18q[1] << "+/- " << sqrt(Fnorm_fraction_LHC18q[2]) / Fnorm_fraction_LHC18q[1] << std::endl;

//****************************************************
//
// For Fnorm obtained from offline indirect method
//
//****************************************************

  name_Fnorm_histo[1] = "fFNormIndirect_CINT7";
  avg_Fnorm_LHC18q[1] = read_Fnorm_cal_avg(path_to_Fnorm[1], name_Fnorm_histo[1], analysis_run_files_location, number_of_runs_LHC18q, tstr_alicounter_path, Fnorm_fraction_LHC18q );
  avg_Fnorm_err_LHC18q[1] = sqrt(Fnorm_fraction_LHC18q[2]) / Fnorm_fraction_LHC18q[1];
  std::cout << Fnorm_fraction_LHC18q[0] << " / " << Fnorm_fraction_LHC18q[1] << "+/- " << avg_Fnorm_err_LHC18q[1] << std::endl;

//****************************************************
//
// For Fnorm obtained from online method
//
//****************************************************

  name_Fnorm_histo[1] = "fFNormOnline_C0V0M";
  avg_Fnorm_LHC18q[2] = read_Fnorm_cal_avg(path_to_Fnorm[1], name_Fnorm_histo[1], analysis_run_files_location, number_of_runs_LHC18q, tstr_alicounter_path, Fnorm_fraction_LHC18q );
  avg_Fnorm_err_LHC18q[2] = sqrt(Fnorm_fraction_LHC18q[2]) / Fnorm_fraction_LHC18q[1];
  std::cout << Fnorm_fraction_LHC18q[0] << " / " << Fnorm_fraction_LHC18q[1] << "+/- " << avg_Fnorm_err_LHC18q[2] << std::endl;

  final_Fnorm_LHC18q = mean_with_stat_weight(avg_Fnorm_LHC18q, avg_Fnorm_err_LHC18q, number_methods);
  final_Fnorm_err_LHC18q = sqrt(1/(pow(1/avg_Fnorm_err_LHC18q[0],2) + pow(1/avg_Fnorm_err_LHC18q[1],2) + pow(1/avg_Fnorm_err_LHC18q[2],2)));
  std::cout << final_Fnorm_LHC18q << "+/-" << final_Fnorm_err_LHC18q << std::endl;

//******************************************
//
// The calculation of Fnorm for LHC18r. Using the same formulas which are shown above.
//
//******************************************

  path_to_file = "FNormRunByRun_18r.root";
  path_to_Fnorm[2] = Form("%s%s", path_to_directory.Data(), path_to_file.Data() );

  name_Fnorm_histo[2] = "fFNormDirect_CINT7";

//  analysis_run_files_location = Form("/home/chunlu/local/runList_LHC18r/muon_calo_pass2/LHC18r_runList.txt");
//  const int number_of_runs_LHC18r = 99;

  analysis_run_files_location = Form("/home/chunlu/local/runList_LHC18q_r/muon_calo_pass3/LHC18r_runList.txt");
  const int number_of_runs_LHC18r = 98;

  double avg_Fnorm_LHC18r[number_methods] = {0};
  double avg_Fnorm_err_LHC18r[number_methods] = {0};
  double Fnorm_fraction_LHC18r[3]={0}; // 
  double final_Fnorm_LHC18r = 0;
  double final_Fnorm_err_LHC18r=0;

//****************************************************
//
// For Fnorm obtained from offline direct method
//
//****************************************************

  avg_Fnorm_LHC18r[0] = read_Fnorm_cal_avg(path_to_Fnorm[2], name_Fnorm_histo[2], analysis_run_files_location, number_of_runs_LHC18r, tstr_alicounter_path, Fnorm_fraction_LHC18r );
  avg_Fnorm_err_LHC18r[0] = sqrt(Fnorm_fraction_LHC18r[2]) / Fnorm_fraction_LHC18r[1];
  std::cout << Fnorm_fraction_LHC18r[0] << " / " << Fnorm_fraction_LHC18r[1] << " +/- " << sqrt(Fnorm_fraction_LHC18r[2]) / Fnorm_fraction_LHC18r[1]  << std::endl;

//****************************************************
//
// For Fnorm obtained from offline indirect method
//
//****************************************************

  name_Fnorm_histo[2] = "fFNormIndirect_CINT7";
  avg_Fnorm_LHC18r[1] = read_Fnorm_cal_avg(path_to_Fnorm[2], name_Fnorm_histo[2], analysis_run_files_location, number_of_runs_LHC18r, tstr_alicounter_path, Fnorm_fraction_LHC18r );
  avg_Fnorm_err_LHC18r[1] = sqrt(Fnorm_fraction_LHC18r[2]) / Fnorm_fraction_LHC18r[1];
  std::cout << Fnorm_fraction_LHC18r[0] << " / " << Fnorm_fraction_LHC18r[1] << " +/- " << avg_Fnorm_err_LHC18r[1]  << std::endl;

//****************************************************
//
// For Fnorm obtained from online method
//
//****************************************************

  name_Fnorm_histo[2] = "fFNormOnline_C0V0M";
  avg_Fnorm_LHC18r[2] = read_Fnorm_cal_avg(path_to_Fnorm[2], name_Fnorm_histo[2], analysis_run_files_location, number_of_runs_LHC18r, tstr_alicounter_path, Fnorm_fraction_LHC18r );
  avg_Fnorm_err_LHC18r[2] = sqrt(Fnorm_fraction_LHC18r[2]) / Fnorm_fraction_LHC18r[1];
  std::cout << Fnorm_fraction_LHC18r[0] << " / " << Fnorm_fraction_LHC18r[1] << " +/- " << avg_Fnorm_err_LHC18r[2]  << std::endl;

  final_Fnorm_LHC18r = mean_with_stat_weight(avg_Fnorm_LHC18r, avg_Fnorm_err_LHC18r, number_methods);
  final_Fnorm_err_LHC18r = sqrt(1/(pow(1/avg_Fnorm_err_LHC18r[0],2) + pow(1/avg_Fnorm_err_LHC18r[1],2) + pow(1/avg_Fnorm_err_LHC18r[2],2)) );

  std::cout << final_Fnorm_LHC18r << "+/-" << final_Fnorm_err_LHC18r << std::endl;

//************************************
//
// calculate the average Fnorm over the three periods, LHC15o, LHC18q, and LHC18r
//
//************************************

  double Fnorm_total = 0;
  double Fnorm_total_error = 0;

  Fnorm_total = (final_Fnorm_LHC15o*Fnorm_fraction_LHC15o[1] + final_Fnorm_LHC18q*Fnorm_fraction_LHC18q[1] + final_Fnorm_LHC18r*Fnorm_fraction_LHC18r[1]);
  Fnorm_total = Fnorm_total / (Fnorm_fraction_LHC15o[1]+ Fnorm_fraction_LHC18q[1] + Fnorm_fraction_LHC18r[1]);

//  Fnorm_total_error = sqrt(pow(final_Fnorm_err_LHC15o*Fnorm_fraction_LHC15o[1],2) + pow(final_Fnorm_err_LHC18q*Fnorm_fraction_LHC18q[1],2) + pow(final_Fnorm_err_LHC18r*Fnorm_fraction_LHC18r[1],2));
  Fnorm_total_error = sqrt(pow(final_Fnorm_err_LHC15o, 2)* pow(Fnorm_fraction_LHC15o[1],2) + pow(final_Fnorm_err_LHC18q,2)*pow(Fnorm_fraction_LHC18q[1],2) + pow(final_Fnorm_err_LHC18r,2)*pow(Fnorm_fraction_LHC18r[1],2));
  Fnorm_total_error = Fnorm_total_error / (Fnorm_fraction_LHC15o[1]+ Fnorm_fraction_LHC18q[1] + Fnorm_fraction_LHC18r[1]);
  std::cout << "Fnorm total: " << Fnorm_total << " +/- " << Fnorm_total_error << std::endl;

  std::cout << "N_MB: " << Fnorm_total * (Fnorm_fraction_LHC15o[1]+ Fnorm_fraction_LHC18q[1] + Fnorm_fraction_LHC18r[1])  << std::endl;

  

//  Fnorm_total = (avg_Fnorm_LHC18q[0]*Fnorm_fraction_LHC18q[1] + avg_Fnorm_LHC18r[0]*Fnorm_fraction_LHC18r[1]);
//  Fnorm_total = Fnorm_total / (Fnorm_fraction_LHC18q[1] + Fnorm_fraction_LHC18r[1]);
//  std::cout << "Fnorm LHC18q+18r: " << Fnorm_total << std::endl;
 

// calculate the Fnorm directly from LHC15o to LHC18r
  Fnorm_total = 0;
  Fnorm_total = Fnorm_fraction_LHC15o[0] + Fnorm_fraction_LHC18q[0] + Fnorm_fraction_LHC18r[0];
  Fnorm_total = Fnorm_total / (Fnorm_fraction_LHC15o[1]+ Fnorm_fraction_LHC18q[1] + Fnorm_fraction_LHC18r[1]);

  Fnorm_total_error = sqrt(Fnorm_fraction_LHC15o[2] + Fnorm_fraction_LHC18q[2] + Fnorm_fraction_LHC18r[2]);
  Fnorm_total_error = Fnorm_total_error / (Fnorm_fraction_LHC15o[1]+ Fnorm_fraction_LHC18q[1] + Fnorm_fraction_LHC18r[1]);

  std::cout << "Fnorm total: " << Fnorm_total << " +/- " << Fnorm_total_error << std::endl;

//**************************************************
//
// calculate the Fnorm directly from LHC15o to LHC18r
//
//**************************************************
/*
  const int number_of_runs[3] = {137, 130, 99};
  TString analysis_run_files_path[3] = {
                                         Form("/home/chunlu/local/runList_LHC15o/runList.txt"),
                                         Form("/home/chunlu/local/runList_LHC18q/muon_calo_pass2/LHC18q_runList.txt"),
                                         Form("/home/chunlu/local/runList_LHC18r/muon_calo_pass2/LHC18r_runList.txt"),
                                        };
  TString tstr_alicounter_path_[3] = {"$HOME/local/anaedPbPbData2015/LHC15o/muon_calo_pass1/AOD197_CMULEvent_AnaResults_Pcorr/2015",
                                      "$HOME/local/anaedPbPbData2018/LHC18q_r/muon_calo_pass3/AOD225_CMULEvent_AnaResults/2018",
                                      "$HOME/local/anaedPbPbData2018/LHC18q_r/muon_calo_pass3/AOD225_CMULEvent_AnaResults/2018"
                                     };
  Fnorm_fraction_LHC18r

  read_Fnorm_cal_avg_error(path_to_Fnorm, name_Fnorm_histo, analysis_run_files_path, number_of_runs, tstr_alicounter_path, Fnorm_fraction_LHC18r)
*/

//***************************************************  
//
// obtain the CMUL events in a given centrality bin
//
//***************************************************
 
  int cent_bin[2] = {1, 9}; 
  double total_events_LHC15o = 0;
  double total_events_LHC18q = 0;
  double total_events_LHC18r = 0;
  double total_events=0;

  analysis_run_files_location = Form("/home/chunlu/local/runList_LHC15o/runList.txt");
  tstr_alicounter_path = "$HOME/local/anaedPbPbData2015/LHC15o/muon_calo_pass1/AOD229_CMULEvent_AnaResults/";
  total_events_LHC15o = read_alicounter_obtain_events(analysis_run_files_location, number_of_runs, tstr_alicounter_path, cent_bin);

  analysis_run_files_location = Form("/home/chunlu/local/runList_LHC18q/muon_calo_pass2/LHC18q_runList.txt");
  tstr_alicounter_path = "$HOME/local/anaedPbPbData2018/LHC18q_r/muon_calo_pass3/AOD225_CMULEvent_AnaResults/";
  total_events_LHC18q = read_alicounter_obtain_events(analysis_run_files_location, number_of_runs_LHC18q, tstr_alicounter_path, cent_bin);

  analysis_run_files_location = Form("/home/chunlu/local/runList_LHC18r/muon_calo_pass2/LHC18r_runList.txt");
  tstr_alicounter_path = "$HOME/local/anaedPbPbData2018/LHC18q_r/muon_calo_pass3/AOD225_CMULEvent_AnaResults/";
  total_events_LHC18r = read_alicounter_obtain_events(analysis_run_files_location, number_of_runs_LHC18r, tstr_alicounter_path, cent_bin);

  total_events = (total_events_LHC15o + total_events_LHC18q + total_events_LHC18r);

  std::cout << "centrality 0-90%" << std::endl;
  std::cout << "LHC15o: " << total_events_LHC15o << std::endl;
  std::cout << "LHC18q: " << total_events_LHC18q << std::endl;
  std::cout << "LHC18r: " << total_events_LHC18r << std::endl;
  std::cout << "LHC15o+18q+18r: " << (total_events_LHC15o+total_events_LHC18q+total_events_LHC18r) << std::endl;
 
  std::cout << "MB: " << Fnorm_total <<" x " << total_events << " = " << (Fnorm_total*total_events) << std::endl;

return 0;
}

double read_alicounter_obtain_events(TString analysis_run_files_location, const int number_of_runs, TString tstr_alicounter_path, int *cent_bin)
{
//********************************************
//
// This function returns the total number of the CMUL events, for a given range, for a given period. 
//
//********************************************

      int run_number[number_of_runs] = {0};
      TString tstr_run_number[number_of_runs]={};
      read_run_number_and_obtain_array(number_of_runs, analysis_run_files_location, run_number, tstr_run_number);
      //std::cout << LHC18q_run_number[0] << std::endl;

      int centrality_bins = 12;
      TString centralityBins_alicounter[] = {"m0","0_10","10_20","20_30","30_40","40_50","50_60","60_70","70_80","80_90","90_100","no_centrality"};
      TString triggerName = Form("CMUL7-B-NOPF-MUFAST");
      TString selection = Form("yes");
      TString selection_ = Form("no");
      TString anaResult_file[number_of_runs];
      TString counterName = Form("eventCounters");

      for(int iRun=0; iRun < number_of_runs; iRun++)  anaResult_file[iRun] = Form("%sPbPb_CMULResultsHistos_Run000%i.root", tstr_alicounter_path.Data(), run_number[iRun] );

      double **numOfevents_Cent_per_Run = new double*[number_of_runs];
      for(int iRun=0; iRun< number_of_runs; iRun++) numOfevents_Cent_per_Run[iRun] = new double[centrality_bins]; //  **numOfevents_Cent_per_Run   {NUMOFFILES  {..MAXINDEX_CENTBINS.. }  NUMOFFILES}

//      double **numOfrejevents_Cent_per_Run = new double*[number_of_runs];
//      for(int iRun=0; iRun< number_of_runs; iRun++) numOfrejevents_Cent_per_Run[iRun] = new double[centrality_bins]; //  **numOfevents_Cent_per_Run   {NUMOFFILES  {..MAXINDEX_CENTBINS.. }  NUMOFFILES}

      read_alicounter(number_of_runs, centrality_bins, run_number, centralityBins_alicounter, triggerName, selection, anaResult_file, counterName, numOfevents_Cent_per_Run);

      double numberOfevents_per_Run[number_of_runs] = {0}; // Must be initiallized, otherwise it obtaines random number in array. 
      double total_numberOfevents = 0;

      for(int iRun=0; iRun< number_of_runs; iRun++)
      {
         for(int iCent = cent_bin[0]; iCent < cent_bin[1]; iCent++)
         {
           numberOfevents_per_Run[iRun] += (numOfevents_Cent_per_Run[iRun][iCent]);
           total_numberOfevents += (numOfevents_Cent_per_Run[iRun][iCent]);
//           std::cout << numOfevents_Cent_per_Run[iRun][iCent] << " ";
         }
//         std::cout << ", " << numberOfevents_per_Run[iRun] << " ";
//         std::cout << std::endl;
      }

return total_numberOfevents;
}

double read_Fnorm_cal_avg(TString path_to_Fnorm, TString name_Fnorm_histo, TString analysis_run_files_location, const int number_of_runs, TString tstr_alicounter_path, double *Fnorm_fraction)
{
  // *************************************
  //
  // This function returns the average Fnorm integrating over the run numbers
  // 
  // *************************************

//      TString path_to_directory = "$HOME/local/anaedPbPbData2015/Fnorm/";
//      TString path_to_file = "Fnorm_offline_LHC15o.root";

//      TString path_to_Fnorm = Form("%s%s", path_to_directory.Data(), path_to_file.Data() );

      TFile *inputFile = new TFile(path_to_Fnorm,"READ");

      const int xi = 0;
      const int xf1 = 137;
      const int NUMOFBINS1 = (xf1-xi);
      const int xf2 = 130;
      const int NUMOFBINS2 = (xf2-xi);
      const int xf3 = 99;
      const int NUMOFBINS3 = (xf3-xi);
      TH1F *Fnorm;
      if( number_of_runs == 137)
      {
        Fnorm = new TH1F("Fnorm_15o","; run number;", NUMOFBINS1, xi, xf1); // Fnorm_offline_LHC15o.root has one more bin
      }
      else if( number_of_runs == 130)
      {
        Fnorm = new TH1F("Fnorm_18q","; run number;", NUMOFBINS2, xi, xf2);
      }
      else if( number_of_runs == 98 )
      {
        Fnorm = new TH1F("Fnorm_18r","; run number;", NUMOFBINS3, xi, xf3);
      }
      else 
      {
        std::cout << "Please check out the number of runs in LHC15o, LHC18q, or LHC18r" << std::endl;
        return 1;
      }

      Fnorm->Add( (TH1F*)inputFile->Get(name_Fnorm_histo.Data() ) );

      double avg_Fnorm = 0;
      double avg_Fnorm_err = 0;
  //
  //
  // read the CMUL evnets in alicounter 
  //
  //

//      const int number_of_runs = xf1;
//      TString analysis_run_files_location = Form("/home/chunlu/local/runList_LHC15o/runList.txt"); // environment variable is not working here. ie. $HOME
//      TString analysis_run_files_location = Form("/home/chunlu/local/runList_LHC15o/runList01.txt");
//      TString analysis_run_files_location = Form("/home/chunlu/cernbox/alice-muon-qa-reports/data/2018/LHC18q/muon_calo_pass1/runListGoodForQA.txt");

      int run_number[number_of_runs] = {0};
      TString tstr_run_number[number_of_runs]={};
      read_run_number_and_obtain_array(number_of_runs, analysis_run_files_location, run_number, tstr_run_number);
      //std::cout << LHC18q_run_number[0] << std::endl;

      int centrality_bins = 12;
      TString centralityBins_alicounter[] = {"m0","0_10","10_20","20_30","30_40","40_50","50_60","60_70","70_80","80_90","90_100","no_centrality"};
      TString triggerName = Form("CMUL7-B-NOPF-MUFAST");
      TString selection = Form("yes");
      TString selection_ = Form("no");
      TString anaResult_file[number_of_runs];
      TString counterName = Form("DimuonHistosUnlike/listOfObjects_CMUL7/eventCounters");

      for(int iRun=0; iRun < number_of_runs; iRun++)  anaResult_file[iRun] = Form("%sDimuonResultsHistos_Run000%i.root", tstr_alicounter_path.Data(), run_number[iRun] );

      double **numOfevents_Cent_per_Run = new double*[number_of_runs];
      for(int iRun=0; iRun< number_of_runs; iRun++) numOfevents_Cent_per_Run[iRun] = new double[centrality_bins]; //  **numOfevents_Cent_per_Run   {NUMOFFILES  {..MAXINDEX_CENTBINS.. }  NUMOFFILES}

      double **numOfEvents_Cent_per_Run = new double*[number_of_runs];
      for(int iRun=0; iRun< number_of_runs; iRun++) numOfEvents_Cent_per_Run[iRun] = new double[centrality_bins];
//      double **numOfrejevents_Cent_per_Run = new double*[number_of_runs];
//      for(int iRun=0; iRun< number_of_runs; iRun++) numOfrejevents_Cent_per_Run[iRun] = new double[centrality_bins]; //  **numOfevents_Cent_per_Run   {NUMOFFILES  {..MAXINDEX_CENTBINS.. }  NUMOFFILES}

      double numberOfevents_per_Run[number_of_runs] = {0}; // Must be initiallized, otherwise it obtaines random number in array. 
      double total_numberOfevents = 0;

      if( number_of_runs == 137)
      {
         for(int iTrig=0; iTrig<2; iTrig++)
         {
           if(iTrig == 0)
           {
             triggerName = Form("CMLL7-B-NOPF-MUFAST&CMUL7-B-NOPF-MUFAST");
             counterName = Form("DimuonHistosUnlike/listOfObjects_CMULandCMLL/eventCounters");
             read_alicounter(number_of_runs, centrality_bins, run_number, centralityBins_alicounter, triggerName, selection, anaResult_file, counterName, numOfevents_Cent_per_Run);
           }
           else if(iTrig == 1)
           {
             triggerName = Form("CMUL7-B-NOPF-MUFAST&!CMLL7-B-NOPF-MUFAST");
             counterName = Form("DimuonHistosUnlike/listOfObjects_CMULnoCMLL/eventCounters");
             read_alicounter(number_of_runs, centrality_bins, run_number, centralityBins_alicounter, triggerName, selection, anaResult_file, counterName, numOfEvents_Cent_per_Run);
           }
          
           for(int iRun=0; iRun< number_of_runs; iRun++)
           {
             for(int iCent = 1; iCent < (centrality_bins-2); iCent++)
             {
               if(iTrig == 0) 
               {
                 numberOfevents_per_Run[iRun] += (numOfevents_Cent_per_Run[iRun][iCent]);
//                 total_numberOfevents += (numOfevents_Cent_per_Run[iRun][iCent]);
    //           std::cout << numOfevents_Cent_per_Run[iRun][iCent] << " ";
               }
               else if(iTrig == 1)
               {
                 numberOfevents_per_Run[iRun] += (numOfEvents_Cent_per_Run[iRun][iCent]);
                 total_numberOfevents += (numOfevents_Cent_per_Run[iRun][iCent] + numOfEvents_Cent_per_Run[iRun][iCent]);
               }
             }
//         std::cout << ", " << numberOfevents_per_Run[iRun] << " ";
//         std::cout << std::endl;
            }
 
         }         

      }
      else
      {
        read_alicounter(number_of_runs, centrality_bins, run_number, centralityBins_alicounter, triggerName, selection, anaResult_file, counterName, numOfevents_Cent_per_Run); 
//      read_alicounter(number_of_runs, centrality_bins, run_number, centralityBins_alicounter, triggerName, selection_, anaResult_file, counterName, numOfrejevents_Cent_per_Run);
        for(int iRun=0; iRun< number_of_runs; iRun++)
        {
                      
           for(int iCent = 1; iCent < (centrality_bins-2); iCent++)
           {
             if(run_number[iRun] == 296977)
             {
                numOfevents_Cent_per_Run[iRun][iCent] = 0.;
             } 

             numberOfevents_per_Run[iRun] += (numOfevents_Cent_per_Run[iRun][iCent]);
             total_numberOfevents += (numOfevents_Cent_per_Run[iRun][iCent]);
//           std::cout << numOfevents_Cent_per_Run[iRun][iCent] << " ";
           }
//         std::cout << ", " << numberOfevents_per_Run[iRun] << " ";
//         std::cout << std::endl;
        }

      }


//      for(int iRun=0; iRun< number_of_runs; iRun++) std::cout << numberOfevents_per_Run[iRun] << std::endl;
        if( number_of_runs == 98 )
        {
          for(int iRun = 0; iRun< number_of_runs; iRun++)
          {
            if((iRun+1) < 43) 
            {
              avg_Fnorm += ( (Fnorm->GetBinContent(iRun+1)) * numberOfevents_per_Run[iRun] );
              avg_Fnorm_err += pow(Fnorm->GetBinError(iRun+1), 2) * pow(numberOfevents_per_Run[iRun], 2);
            }
            else 
            {
              avg_Fnorm += ( (Fnorm->GetBinContent(iRun+2)) * numberOfevents_per_Run[iRun] );
              avg_Fnorm_err += pow(Fnorm->GetBinError(iRun+2), 2) * pow(numberOfevents_per_Run[iRun], 2);
            }
//        std::cout << ( (Fnorm[0]->GetBinContent(iRun+1)) ) << std::endl;
//        std::cout << ( numberOfevents_per_Run[iRun] ) << std::endl;
          }
        }
        else 
        {
          for(int iRun = 0; iRun< number_of_runs; iRun++)
          {        
            avg_Fnorm += ( (Fnorm->GetBinContent(iRun+1)) * numberOfevents_per_Run[iRun] ); 
            avg_Fnorm_err += pow(Fnorm->GetBinError(iRun+1), 2) * pow(numberOfevents_per_Run[iRun], 2);
//        std::cout << ( (Fnorm[0]->GetBinContent(iRun+1)) ) << std::endl;
//        std::cout << ( numberOfevents_per_Run[iRun] ) << std::endl;
          }
        }
      
//      std::cout << avg_Fnorm << std::endl;
      Fnorm_fraction[0] = avg_Fnorm;
      Fnorm_fraction[1] = total_numberOfevents;
      Fnorm_fraction[2] = avg_Fnorm_err;

//      std::cout << avg_Fnorm << std::endl;

      avg_Fnorm = avg_Fnorm / total_numberOfevents; 

//      std::cout << total_numberOfevents << std::endl;
      std::cout << avg_Fnorm << std::endl;
      delete Fnorm;
      inputFile->Close();

      for(int iRun=0; iRun < number_of_runs; iRun++) delete []  numOfevents_Cent_per_Run[iRun];
      delete [] numOfevents_Cent_per_Run;

      for(int iRun=0; iRun < number_of_runs; iRun++) delete []  numOfEvents_Cent_per_Run[iRun];
      delete [] numOfEvents_Cent_per_Run;

return avg_Fnorm;
}

double read_Fnorm_cal_avg_error(TString *path_to_Fnorm, TString *name_Fnorm_histo, TString *analysis_run_files_location, const int *number_of_runs, TString *tstr_alicounter_path, double *Fnorm_fraction)
{
      TFile *inputFile[3];
      inputFile[0] = new TFile(path_to_Fnorm[0],"READ");
      inputFile[1] = new TFile(path_to_Fnorm[1],"READ");
      inputFile[2] = new TFile(path_to_Fnorm[2],"READ");  

      const int xi = 0;
      const int xf1 = 137;
      const int NUMOFBINS1 = (xf1-xi)+1;
      const int xf2 = 130;
      const int NUMOFBINS2 = (xf2-xi);
      const int xf3 = 99;
      const int NUMOFBINS3 = (xf3-xi);

      TH1F *Fnorm[3];
      Fnorm[0] = new TH1F("Fnorm_15o","; run number;", NUMOFBINS1, xi, xf1+1); // Fnorm_offline_LHC15o.root has one more bin
      Fnorm[1] = new TH1F("Fnorm_18q","; run number;", NUMOFBINS2, xi, xf2);
      Fnorm[2] = new TH1F("Fnorm_18r","; run number;", NUMOFBINS3, xi, xf3);

      Fnorm[0]->Add( (TH1F*)inputFile[0]->Get(name_Fnorm_histo[0].Data() ) );
      Fnorm[1]->Add( (TH1F*)inputFile[1]->Get(name_Fnorm_histo[1].Data() ) );
      Fnorm[2]->Add( (TH1F*)inputFile[2]->Get(name_Fnorm_histo[2].Data() ) );

      double avg_Fnorm = 0;
      double avg_Fnorm_err = 0;


      int run_number_LHC15o[number_of_runs[0]] = {0};
      TString tstr_run_number_LHC15o[number_of_runs[0]]={};
      read_run_number_and_obtain_array(number_of_runs[0], analysis_run_files_location[0], run_number_LHC15o, tstr_run_number_LHC15o);

      int run_number_LHC18q[number_of_runs[1]] = {0};
      TString tstr_run_number_LHC18q[number_of_runs[1]]={};
      read_run_number_and_obtain_array(number_of_runs[1], analysis_run_files_location[1], run_number_LHC18q, tstr_run_number_LHC18q);

      int run_number_LHC18r[number_of_runs[2]] = {0};
      TString tstr_run_number_LHC18r[number_of_runs[2]]={};
      read_run_number_and_obtain_array(number_of_runs[1], analysis_run_files_location[2], run_number_LHC18r, tstr_run_number_LHC18r);


      int centrality_bins = 12;
      TString centralityBins_alicounter[] = {"m0","0_10","10_20","20_30","30_40","40_50","50_60","60_70","70_80","80_90","90_100","no_centrality"};
      TString triggerName = Form("CMUL7-B-NOPF-MUFAST");
      TString selection = Form("yes");
      TString selection_ = Form("no");
      TString counterName = Form("eventCounters");

      TString anaResult_file_LHC15o[number_of_runs[0]];


      for(int iRun=0; iRun < number_of_runs[0]; iRun++)  anaResult_file_LHC15o[iRun] = Form("%sPbPb_CMULResultsHistos_Run000%i.root", tstr_alicounter_path[0].Data(), run_number_LHC15o[iRun] );

      double **numOfevents_Cent_per_Run = new double*[number_of_runs[0]];
      for(int iRun=0; iRun< number_of_runs[0]; iRun++) numOfevents_Cent_per_Run[iRun] = new double[centrality_bins];



      delete Fnorm[0];
      delete Fnorm[1];
      delete Fnorm[2];
      inputFile[0]->Close();
      inputFile[1]->Close();
      inputFile[2]->Close();

      for(int iRun=0; iRun < number_of_runs[0]; iRun++) delete []  numOfevents_Cent_per_Run[iRun];
      delete [] numOfevents_Cent_per_Run;

}

double mean_with_stat_weight(double *value, double* err_value, int number_methods)
{
// ****************************
//
// This function returns a mathmetical average by taking into account the statistical uncertainty as weights. 
//
// ****************************
	
  double mean_weighted=0;
  double numerator = 0;
  double denominator = 0;

  for(int i=0; i<number_methods; i++)
  {
   numerator += (value[i]/err_value[i]/err_value[i]);
   denominator += (1/err_value[i]/err_value[i]);
  }

  mean_weighted = numerator / denominator;

return mean_weighted;
}
