//--------------------------------------------------------------------------
// Base macro for analyzing the analysis task output file.
// Usage: `root-config --libs --cflags` -I${ALICE_ROOT}/include -L${ALICE_ROOT}/lib -lSTEERBase
// Input: 
// Output: 
//--------------------------------------------------------------------------

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
#include "Math/MinimizerOptions.h" //This contains the function:  ROOT::Math::MinimizerOptions::SetDefaultMaxFunctionCalls(50000);


//PWG includes
#include "AliCounterCollection.h"

#include "gsl/gsl_errno.h"
#include "TError.h"

//Self library
#include "PtCentRangeValues.h"
#include "FitFunction.h"
#include "FitFunctionTails.h"
#include "FitFunctionMaterials.h"
#include "MathTechnicalFunc.h"
#include "tail_txt_file_path.h"
#include "variables.h"
#include "read_access_histos.h"
#include "read_access_histos_y_dependence.h"
#include "read_access_histos_double_diff.h"

//#include "fit_params/fit_params_0to0p3.h"
//#include "fit_params/fit_params_0to0p3_cent60to90.h"
//#include "fit_params/fit_params_0p3to1.h"
//#include "fit_params/fit_params_1to2.h"
//#include "fit_params/fit_params_1to2_cent20to30.h"
//#include "fit_params/fit_params_2to3.h"
//#include "fit_params/fit_params_2to3_cent10to20.h"
//#include "fit_params/fit_params_2to3_cent40to50.h"
//#include "fit_params/fit_params_2to3_cent60to90.h"

//#include "fit_params/fit_params_3to4.h"
//#include "fit_params/fit_params_3to4_cent40to50.h"
//#include "fit_params/fit_params_4to5.h"
//#include "fit_params/fit_params_4to5_cent30to40.h"
//#include "fit_params/fit_params_4to5_cent60to90.h"
//#include "fit_params/fit_params_5to6.h"
//#include "fit_params/fit_params_5to6_cent20to30.h"
//#include "fit_params/fit_params_5to6_cent60to90.h"

//#include "fit_params/fit_params_6to7.h"
//#include "fit_params/fit_params_6to7_cent40to50.h"
//#include "fit_params/fit_params_6to7_cent50to60.h"
//#include "fit_params/fit_params_7to8.h"
//#include "fit_params/fit_params_7to8_cent20to30.h"
//#include "fit_params/fit_params_7to8_cent50to60.h"
//#include "fit_params/fit_params_7to8_cent60to90.h"

//#include "fit_params/fit_params_8to9.h"
//#include "fit_params/fit_params_8to9_cent20to30.h" // for cent0to10, 40to50ma also
//#include "fit_params/fit_params_8to9_cent30to40.h"
//#include "fit_params/fit_params_9to10.h" // use fit_params_10to11.h for cent60to90
//#include "fit_params/fit_params_9to10_cent10to20.h"
//#include "fit_params/fit_params_9to10_cent30to40.h" // for 40to50 also
//#include "fit_params/fit_params_10to11.h"
//#include "fit_params/fit_params_10to11_cent10to20.h"
//#include "fit_params/fit_params_10to11_cent20to30.h"
//#include "fit_params/fit_params_11to12.h"
//#include "fit_params/fit_params_11to12_cent10to20.h" // for cent20to30, 30to40 also
//#include "fit_params/fit_params_11to12_cent40to50.h"
//#include "fit_params/fit_params_12to15.h"
//#include "fit_params/fit_params_12to15_cent60to90.h"
//#include "fit_params/fit_params_15to20.h"
//#include "fit_params/fit_params_15to20_cent30to40.h"
//#include "fit_params/fit_params_15to20_cent60to90.h"

//#include "fit_params/fit_params_3p5to3p75_cent0to10.h"
//#include "fit_params/fit_params_2p5to2p75_cent10to20.h"
//#include "fit_params/fit_params_0p3to2_y3_3p5.h"
#include "fit_params/fit_params_2to4_y2p5_3.h"

struct fit_variable_results
{
//  const int x =5;
//  const int y =5;
//  const int z =5;
  double Njpsi_CB2VWG;

  int fit_valide;
  int fit_status;
  int fit_cov_matrix_status;

};

//struct fit_variable_results fitJPsi_Psi2S(int, double, double, TH1F*, TF1*, TF1*, TF1*, TF1*, TF1*,  struct paramOfFunc, double *, double *, double *, double *, double *, double *, double *,double *, TCanvas*, TCanvas*, struct fit_variable_results);
struct fit_variable_results fitJPsi_Psi2S(int, double, double, TH1F*, TF1*, TF1*, TF1*, TF1*, TF1*,  struct paramOfFunc, double *, double *, double *, double *, double *, double *, double *,double *, TCanvas*, struct fit_variable_results);
double calRMS(int ,double *, double );
double calWeightAverage(int , double *, double *);
double calWeightRMS(int , double *, double *);
//double calWeightRMS
//TString opt="IRS0"; //fitting parameters
TString opt="IRSL0";
//"R" Use the Range specified in the function range
//"S" The result of the fit is returned in the TFitResultPtr 
//"L" Likelihood
//"0" Do not plot the function after fit. Must be used because the global fit function is draw also so use this option to avoid drawing twice at the end

int main()
{
     ROOT::Math::MinimizerOptions::SetDefaultMaxFunctionCalls(60000);//This function allows the number of calling fitting up to 50000 
//     gErrorIgnoreLevel = kError+1;

     arrayAvailableCentBins[0] = CentMin;
     for(int i=1; i<nCentIndex; i++)  
     {
         arrayAvailableCentBins[i] = arrayAvailableCentBins[i-1] + CentInterval;
         std::cout << "arrayAvailableCentBins: " << arrayAvailableCentBins[i] << std::endl;
     }

//--------------------------------------------------------------------------------------
//
// Initialize the Histograms
//
//--------------------------------------------------------------------------------------  

    TH1F *histoRawInvM = new TH1F(Form("histoRaw%sInvM", dimuonChargeNames[1].Data() ), ";M_{#mu#mu} (GeV/c^{2});Nevents", NUMOFBINS,lowerLimit_histoInvMass,upperLimit_histoInvMass);
    histoRawInvM->Sumw2();

    TH1F *histoInvMass_Pt12to15 = new TH1F(Form("histoInvMass_Merged_Pt12to15"), Form("dimuon invariant mass distribution vs Pt12to15; M_{#mu#mu} (GeV/c^{2}); # of events per %g GeV", (upperLimit_histoInvMass-lowerLimit_histoInvMass)/NUMOFBINS ), NUMOFBINS,lowerLimit_histoInvMass,upperLimit_histoInvMass);
      histoInvMass_Pt12to15->Sumw2();
      histoInvMass_Pt12to15->SetAxisRange(2,5,"X");

    TH1F *histoInvMass_Pt15to20 = new TH1F(Form("histoInvMass_Merged_Pt15to20"), Form("dimuon invariant mass distribution vs Pt15to20; M_{#mu#mu} (GeV/c^{2}); Events per %g GeV",(upperLimit_histoInvMass-lowerLimit_histoInvMass)/NUMOFBINS ), NUMOFBINS,lowerLimit_histoInvMass,upperLimit_histoInvMass);
      histoInvMass_Pt15to20->Sumw2();

    
    TH1F *histoInvMass[numberOfPtBins_photon];
    for(int iPt =0; iPt < (numberOfPtBins_photon); iPt++) 
    {
         histoInvMass[iPt] = new TH1F(Form("histoInvMass_Pt%gto%g", arrayAvailablePtBins_photon[iPt], arrayAvailablePtBins_photon[iPt+1]), Form("; m_{#mu#mu} (GeV/c^{2}); Events per %g GeV", (upperLimit_histoInvMass-lowerLimit_histoInvMass)/NUMOFBINS ), NUMOFBINS, lowerLimit_histoInvMass, upperLimit_histoInvMass);
         histoInvMass[iPt]->Sumw2();
         histoInvMass[iPt]->SetAxisRange(2,5,"X");
    }

    double arrayAvailableYBins[] = {2.5, 2.75, 3.0, 3.25, 3.5, 3.75, 4.0};
    const int numberOfYBins = 6;

    TH1F *histoInvMass_Ydependence[numberOfYBins];
    for(int iY =0; iY < (numberOfYBins); iY++)
    {
         histoInvMass_Ydependence[iY] = new TH1F(Form("histoInvMass_Y%gto%g", arrayAvailableYBins[iY], arrayAvailableYBins[iY+1]), Form("; m_{#mu#mu} (GeV/c^{2}); Events per %g GeV", (upperLimit_histoInvMass-lowerLimit_histoInvMass)/NUMOFBINS ), NUMOFBINS, lowerLimit_histoInvMass, upperLimit_histoInvMass);
         histoInvMass_Ydependence[iY]->Sumw2();
         histoInvMass_Ydependence[iY]->SetAxisRange(2,5,"X");
    }

    double arrayPtBins_double_diff[] = {0.3,2,4,6,12,20};
    const int numberOfPtBins_double_diff = 5;
    double arrayYBins_double_diff[] = {2.5,3,3.5,4};
    const int numberOfYBins_double_diff = 3;
    TH1F ***histoInvMass_double_diff;    
    
    histoInvMass_double_diff = new TH1F**[numberOfPtBins_double_diff];
    for(int ipT =0; ipT < numberOfPtBins_double_diff; ipT++)
    {
      histoInvMass_double_diff[ipT] = new TH1F*[numberOfYBins_double_diff];
         std::cout << "ipT = 0" << std::endl;
    }

    for(int ipT =0; ipT < numberOfPtBins_double_diff; ipT++)
    {
      for(int iY =0; iY < numberOfYBins_double_diff; iY++)
      {
         histoInvMass_double_diff[ipT][iY] = new TH1F(Form("histoInvMass_pT%gto%g_Y%gto%g", arrayPtBins_double_diff[ipT], arrayPtBins_double_diff[ipT+1], arrayYBins_double_diff[iY], arrayYBins_double_diff[iY+1]), Form("; m_{#mu#mu} (GeV/c^{2}); Events per %g GeV", (upperLimit_histoInvMass-lowerLimit_histoInvMass)/NUMOFBINS ), NUMOFBINS, lowerLimit_histoInvMass, upperLimit_histoInvMass);
         histoInvMass_double_diff[ipT][iY]->Sumw2();
         histoInvMass_double_diff[ipT][iY]->SetAxisRange(2,5,"X");
      }

    }
 //--------------------------------------------------------------------------------------
 //
 // Read raw root file to be ready for the fitting
 //
 //--------------------------------------------------------------------------------------
    TFile *inputFile;    
    bool iMerge = true;
//    bool iY = true;
    bool iY = false;
    bool ipT_Y =true;
    int init_pt = 0; // indicate initial pT bin for mass spectra
    int init_y = 0;


    int init_cent = 1; // indicate initial cent bin for mass spectra
    int end_cent = 2; //indicate end cent bin for mass spectra

    if(iMerge)
    {
      for(int iFile=0; iFile<3; iFile++)
      {
        TString directory;
        TString name ;
        if(iFile ==0) 
        {
//          directory = "$HOME/local/anaedPbPbData2018/LHC18q_r/muon_calo_pass3/AOD225_CMULEvent_AnaResults_Pcorr/";
//          directory = "$HOME/local/anaedPbPbData2018/LHC18q_r/muon_calo_pass3/AOD225_CMULEvent_AnaResults/";
          directory = "$HOME/local/anaedPbPbData2018/LHC18q_r/muon_calo_pass3/AOD225_CMULEvent_AnaResults/double_differential_pT_y/";
          name = "DimuonResultsHistos_LHC18q.root";
        }
        else if(iFile ==1)
        {
//          directory = "$HOME/local/anaedPbPbData2018/LHC18q_r/muon_calo_pass3/AOD225_CMULEvent_AnaResults_Pcorr/";
//          directory = "$HOME/local/anaedPbPbData2018/LHC18q_r/muon_calo_pass3/AOD225_CMULEvent_AnaResults/";
          directory = "$HOME/local/anaedPbPbData2018/LHC18q_r/muon_calo_pass3/AOD225_CMULEvent_AnaResults/double_differential_pT_y/";
          name = "DimuonResultsHistos_LHC18r.root";
        }
        else if(iFile ==2)
        {
          directory = "$HOME/local/anaedPbPbData2015/LHC15o/muon_calo_pass1/AOD229_CMULEvent_AnaResults/double_differential_pT_y/";
//          directory = "$HOME/local/anaedPbPbData2015/LHC15o/muon_calo_pass1/AOD197_CMULEvent_AnaResults/";
          name = "DimuonResultsHistos_LHC15o.root";
        }
        else
        {
          break;
        }

        TString directory_name = Form("%s%s", directory.Data(), name.Data() );

        TFile *inputFile = new TFile(directory_name,"READ");

        if(iY)
        {
           std::cout << "iY test" << std::endl;
          read_access_histos_y_dependence(inputFile, init_y, init_cent, end_cent, histoInvMass, histoRawInvM);
          pt_bin = y_bin;
        }
        else if(ipT_Y)
        {
          read_access_histos_double_diff(inputFile, 0, init_cent, end_cent, histoInvMass_double_diff, histoRawInvM);
        }
        else
        {
          read_access_histos(inputFile, init_pt, init_cent, end_cent, histoInvMass, histoRawInvM);
        }
 
        inputFile->Close();
      }
    }
    else 
    {
        TFile *inputFile = new TFile(path_to_histos,"READ");

        read_access_histos(inputFile, init_pt, init_cent, end_cent, histoInvMass, histoRawInvM);

        inputFile->Close();
    }
//    inputFile_LHC18q_after_corr->Close();
//    inputFile_LHC18r->Close();  
  
    delete arrayAvailableCentBins; 
//-------------------------------------------------------------------------------------------------------------------------
//
//Fitting 
//
//-------------------------------------------------------------------------------------------------------------------------   

//  if(pt_bin == 12 || pt_bin ==13 || pt_bin ==14)
  if(pt_bin == 116)
//  if(pt_bin == 12)
//  if(pt_bin == 10)
  {
    histoInvMass[pt_bin]->Rebin();
    histoInvMass[pt_bin]->SetAxisRange(2,5,"X");
     nbin = histoInvMass[pt_bin]->GetNbinsX();
    _xmin = histoInvMass[pt_bin]->GetXaxis()->GetBinLowEdge(1);
    _xmax = histoInvMass[pt_bin]->GetXaxis()->GetBinUpEdge(nbin);
    binWidth = (_xmax-_xmin)/nbin;

    histoInvMass[pt_bin]->SetTitle(Form("; m_{#mu#mu} (GeV/c^{2}); Events per %g GeV",binWidth ));
  }
  else
  {
    nbin = histoInvMass[0]->GetNbinsX();
    _xmin = histoInvMass[0]->GetXaxis()->GetBinLowEdge(1);
    _xmax = histoInvMass[0]->GetXaxis()->GetBinUpEdge(nbin);
    binWidth = (_xmax-_xmin)/nbin;
  }


  struct paramOfFunc paramCB2VWG[NUMOFTAILPARAM_CB2];
  paramCB2VWG[0].strNameOfFitFunc = "fitFuncCB2VWG_pp13";
  paramCB2VWG[1].strNameOfFitFunc = "fitFuncCB2VWG_geant3";
//  paramCB2VWG[2].strNameOfFitFunc = "fitFuncCB2VWG_geant3";

  struct paramOfFunc paramCB2Pol2[NUMOFTAILPARAM_CB2];
  paramCB2Pol2[0].strNameOfFitFunc = "fitFuncCB2Pol_pp13";
  paramCB2Pol2[1].strNameOfFitFunc = "fitFuncCB2Pol_geant3";
//  paramCB2Pol2[2].strNameOfFitFunc = "fitFuncCB2Pol_geant3";  

  double Njpsi_CB2VWG[FITRANGE][NUMOFTAILPARAM_CB2][NUMOFCENT_BINS];
  double SigmaJpsi_CB2VWG[FITRANGE][NUMOFTAILPARAM_CB2][NUMOFCENT_BINS];
  double MassJpsi_CB2VWG[FITRANGE][NUMOFTAILPARAM_CB2][NUMOFCENT_BINS];
  double ChiJpsi_CB2VWG[FITRANGE][NUMOFTAILPARAM_CB2][NUMOFCENT_BINS];
  double ErrorOfNjpsi_CB2VWG[FITRANGE][NUMOFTAILPARAM_CB2][NUMOFCENT_BINS];
  double ErrorOfSigmaJpsi_CB2VWG[FITRANGE][NUMOFTAILPARAM_CB2][NUMOFCENT_BINS];
  double ErrorOfMassJpsi_CB2VWG[FITRANGE][NUMOFTAILPARAM_CB2][NUMOFCENT_BINS];
  double ErrorOfChiJpsi_CB2VWG[FITRANGE][NUMOFTAILPARAM_CB2][NUMOFCENT_BINS];

  double Njpsi_CB2Pol2[FITRANGE][NUMOFTAILPARAM_CB2][NUMOFCENT_BINS];
  double SigmaJpsi_CB2Pol2[FITRANGE][NUMOFTAILPARAM_CB2][NUMOFCENT_BINS];
  double MassJpsi_CB2Pol2[FITRANGE][NUMOFTAILPARAM_CB2][NUMOFCENT_BINS];
  double ChiJpsi_CB2Pol2[FITRANGE][NUMOFTAILPARAM_CB2][NUMOFCENT_BINS];
  double ErrorOfNjpsi_CB2Pol2[FITRANGE][NUMOFTAILPARAM_CB2][NUMOFCENT_BINS];
  double ErrorOfSigmaJpsi_CB2Pol2[FITRANGE][NUMOFTAILPARAM_CB2][NUMOFCENT_BINS];
  double ErrorOfMassJpsi_CB2Pol2[FITRANGE][NUMOFTAILPARAM_CB2][NUMOFCENT_BINS];
  double ErrorOfChiJpsi_CB2Pol2[FITRANGE][NUMOFTAILPARAM_CB2][NUMOFCENT_BINS];

  
  for (int iTail=0; iTail < NUMOFTAILPARAM_CB2; iTail++)
  {
    // initialize the struct for CB2 + VWG
    paramCB2VWG[iTail].binWidth = binWidth;
    paramCB2VWG[iTail].NumOfBackGroundFuncParams = numOfVWG_Parm;
    paramCB2VWG[iTail].NumOfJPsiSignalFuncParams = numOfCB2_Parm+1;

    paramCB2VWG[iTail].NameOfBackGroundFuncParams = new TString[numOfVWG_Parm];
    for (int iName=0; iName<numOfVWG_Parm; iName++) paramCB2VWG[iTail].NameOfBackGroundFuncParams[iName].Append(NameOfVWGFuncParams[iName].Data());

    paramCB2VWG[iTail].NameOfJPsiSignalFuncParams = new TString[numOfCB2_Parm+1];
    for (int iName=0; iName<numOfCB2_Parm+1; iName++) paramCB2VWG[iTail].NameOfJPsiSignalFuncParams[iName].Append(NameOfDoubleCB2FuncParams[iName].Data());
   
    paramCB2VWG[iTail].ParamJPsiFunc = new double[numOfCB2_Parm+1];

    // initialize the struct for CB2 + Pol2overPol3
    paramCB2Pol2[iTail].binWidth = binWidth;
    paramCB2Pol2[iTail].NumOfBackGroundFuncParams = numOfPol2overPol3_Parm;
    paramCB2Pol2[iTail].NumOfJPsiSignalFuncParams = numOfCB2_Parm+1;

    paramCB2Pol2[iTail].NameOfBackGroundFuncParams = new TString[numOfPol2overPol3_Parm];
    for (int iName=0; iName<numOfPol2overPol3_Parm; iName++) paramCB2Pol2[iTail].NameOfBackGroundFuncParams[iName].Append(NameOfPol2overPol3FuncParams[iName].Data());

    paramCB2Pol2[iTail].NameOfJPsiSignalFuncParams = new TString[numOfCB2_Parm+1];
    for (int iName=0; iName<numOfCB2_Parm+1; iName++) paramCB2Pol2[iTail].NameOfJPsiSignalFuncParams[iName].Append(NameOfDoubleCB2FuncParams[iName].Data());

    paramCB2Pol2[iTail].ParamJPsiFunc = new double[numOfCB2_Parm+1];

 }

  for(int iParam=0; iParam<numOfCB2_Parm; iParam++)
  {
    paramCB2VWG[0].ParamJPsiFunc[iParam] = Param_doubleCB2_ppdata[iParam];
    paramCB2VWG[1].ParamJPsiFunc[iParam] = Param_doubleCB2_geant3[iParam];
//    paramCB2VWG[2].ParamJPsiFunc[iParam] = Param_doubleCB2_geant3[iParam];
    paramCB2Pol2[0].ParamJPsiFunc[iParam] = Param_doubleCB2_ppdata[iParam];
    paramCB2Pol2[1].ParamJPsiFunc[iParam] = Param_doubleCB2_geant3[iParam];
//    paramCB2Pol2[2].ParamJPsiFunc[iParam] = Param_doubleCB2_geant3[iParam];
  }

   
  for(int iPt=0; iPt < num_pT_bins; iPt++) 
  {
//     values_geant4_CB2_tail_range00[iPt] = new double[numOfCB2_Parm];
//     values_geant4_CB2_tail_range01[iPt] = new double[numOfCB2_Parm];
     values_geant3_CB2_tail_range00[iPt] = new double[numOfCB2_Parm];
     values_geant3_CB2_tail_range01[iPt] = new double[numOfCB2_Parm];

//     values_geant4_NA60_tail_range00[iPt] = new double[numOfNA60_Parm];
//     values_geant4_NA60_tail_range01[iPt] = new double[numOfNA60_Parm];
     values_geant3_NA60_tail_range00[iPt] = new double[numOfNA60_Parm];
     values_geant3_NA60_tail_range01[iPt] = new double[numOfNA60_Parm];
  }


//  read_signal_tail_params(path_geant4_CB2_tail00, numOfCB2_Parm, num_pT_bins, values_geant4_CB2_tail_range00 );
//  read_signal_tail_params(path_geant4_CB2_tail01, numOfCB2_Parm, num_pT_bins, values_geant4_CB2_tail_range01 );
  read_signal_tail_params(path_geant3_CB2_tail00, numOfCB2_Parm, num_pT_bins, values_geant3_CB2_tail_range00 );
  read_signal_tail_params(path_geant3_CB2_tail01, numOfCB2_Parm, num_pT_bins, values_geant3_CB2_tail_range01 );

//  read_signal_tail_params(path_geant4_NA60_tail00, numOfNA60_Parm, num_pT_bins, values_geant4_NA60_tail_range00 );
//  read_signal_tail_params(path_geant4_NA60_tail01, numOfNA60_Parm, num_pT_bins, values_geant4_NA60_tail_range01 );
  read_signal_tail_params(path_geant3_NA60_tail00, numOfNA60_Parm, num_pT_bins, values_geant3_NA60_tail_range00 );
  read_signal_tail_params(path_geant3_NA60_tail01, numOfNA60_Parm, num_pT_bins, values_geant3_NA60_tail_range01 );

  double rangeDo=2.0;
  double rangeUp=5.0;
  const double FitRange[FITRANGE][2] = {{2.2,4.5},{2.4,4.7}}; // for Pb-Pb fitting
//  const double FitRange[FITRANGE][2] = {{2.0,4.8},{2.2,4.4}};  


  const int FIT_NUMBER_CB2 = FITRANGE * 2*  NUMOFTAILPARAM_CB2 * 1 ;
//  TCanvas *fitCanvas_background[FIT_NUMBER_CB2];
  TCanvas *fitCanvas_CB2background[FIT_NUMBER_CB2];

  for(int ifit=0; ifit<FIT_NUMBER_CB2; ifit++)
  {
   fitCanvas_CB2background[ifit] = new TCanvas(Form("c1_%i",ifit),"",200, 10, 800, 600);
//   fitCanvas_background[ifit] = new TCanvas(Form("c2_%i",ifit),"",200, 10, 800, 600);
  }
  int iCount1=0;


  struct fit_variable_results outcome_fit;

  for (int iRange=0; iRange < FITRANGE; iRange++ )
//  for (int iRange=1; iRange < 2; iRange++ )
  {

//    for (int iTail=2; iTail< 3; iTail++)
    for (int iTail=0; iTail< NUMOFTAILPARAM_CB2; iTail++)
    {     
//      paramCB2VWG[iTail].fitRange[iRange][0]=FitRange[iRange][0]; 
//      paramCB2VWG[iTail].fitRange[iRange][1]=FitRange[iRange][1];

      const int numOfVWGCb2_Param = numOfVWG_Parm + numOfCB2_Parm;      

    if(iTail==0 && iRange == 0) // Tail: pp13TeV. Fit Range: 2.0 4.8
    {
      //paramCB2VWG[iTail].strNameOfFitFunc = Form("fitFuncCB2VWG_pp13_%i_%i", iTail, iRange);
      for(int iParam=0; iParam<numOfCB2_Parm; iParam++) paramCB2VWG[iTail].ParamJPsiFunc[iParam] = Param_doubleCB2_ppdata[iParam];
      paramCB2VWG[iTail].paramsInitialVWG[0]= pp13_fit0_par0_CB2_VWG[0];
      paramCB2VWG[iTail].paramsLimitVWG[0][0] = pp13_fit0_par0_CB2_VWG[1];
      paramCB2VWG[iTail].paramsLimitVWG[0][1] = pp13_fit0_par0_CB2_VWG[2];
  
      //background function fit
      paramCB2VWG[iTail].paramsInitialVWG[1] = pp13_fit0_par1_CB2_VWG[0]; 
      paramCB2VWG[iTail].paramsLimitVWG[1][0] = pp13_fit0_par1_CB2_VWG[1];
      paramCB2VWG[iTail].paramsLimitVWG[1][1] = pp13_fit0_par1_CB2_VWG[2];

      paramCB2VWG[iTail].paramsInitialVWG[2] = pp13_fit0_par2_CB2_VWG[0];
      paramCB2VWG[iTail].paramsLimitVWG[2][0] = pp13_fit0_par2_CB2_VWG[1];
      paramCB2VWG[iTail].paramsLimitVWG[2][1] = pp13_fit0_par2_CB2_VWG[2];

      paramCB2VWG[iTail].paramsInitialVWG[3] = pp13_fit0_par3_CB2_VWG[0];
      paramCB2VWG[iTail].paramsLimitVWG[3][0] = pp13_fit0_par3_CB2_VWG[1];
      paramCB2VWG[iTail].paramsLimitVWG[3][1] = pp13_fit0_par3_CB2_VWG[2];
      
      paramCB2VWG[iTail].paramsInitialVWG[4] = pp13_fit0_par4_CB2_VWG[0];
      paramCB2VWG[iTail].paramsLimitVWG[4][0] = pp13_fit0_par4_CB2_VWG[1];
      paramCB2VWG[iTail].paramsLimitVWG[4][1] = pp13_fit0_par4_CB2_VWG[2];


      //background function + signal function fit. If the above range is not nice, adjust the range here.
      paramCB2VWG[iTail].paramsVarianceVWG[1] = pp13_fit0_par1_variance_CB2_VWG[0];
      paramCB2VWG[iTail].paramsLimitVarianceVWG[1][0] = pp13_fit0_par1_variance_CB2_VWG[1];
      paramCB2VWG[iTail].paramsLimitVarianceVWG[1][1] = pp13_fit0_par1_variance_CB2_VWG[2];

      paramCB2VWG[iTail].paramsVarianceVWG[2] = pp13_fit0_par2_variance_CB2_VWG[0];
      paramCB2VWG[iTail].paramsLimitVarianceVWG[2][0] = pp13_fit0_par2_variance_CB2_VWG[1];
      paramCB2VWG[iTail].paramsLimitVarianceVWG[2][1] = pp13_fit0_par2_variance_CB2_VWG[2];

      paramCB2VWG[iTail].paramsVarianceVWG[3] = pp13_fit0_par3_variance_CB2_VWG[0];
      paramCB2VWG[iTail].paramsLimitVarianceVWG[3][0] = pp13_fit0_par3_variance_CB2_VWG[1];
      paramCB2VWG[iTail].paramsLimitVarianceVWG[3][1] = pp13_fit0_par3_variance_CB2_VWG[2];

      paramCB2VWG[iTail].paramsVarianceVWG[4] = pp13_fit0_par4_variance_CB2_VWG[0];
      paramCB2VWG[iTail].paramsLimitVarianceVWG[4][0] = pp13_fit0_par4_variance_CB2_VWG[1];
      paramCB2VWG[iTail].paramsLimitVarianceVWG[4][1] = pp13_fit0_par4_variance_CB2_VWG[2];

      paramCB2VWG[iTail].mean_jpsi = pp13_fit0_mean_jpsi_CB2_VWG[0];
      paramCB2VWG[iTail].mean_LowerLimits_jpsi = pp13_fit0_mean_jpsi_CB2_VWG[1];
      paramCB2VWG[iTail].mean_UpperLimits_jpsi = pp13_fit0_mean_jpsi_CB2_VWG[2];

      paramCB2VWG[iTail].width_jpsi = pp13_fit0_width_jpsi_CB2_VWG[0];
      paramCB2VWG[iTail].width_LowerLimits_jpsi = pp13_fit0_width_jpsi_CB2_VWG[1];
      paramCB2VWG[iTail].width_UpperLimits_jpsi = pp13_fit0_width_jpsi_CB2_VWG[2];

      paramCB2VWG[iTail].Norm_psi2S = pp13_fit0_norm_psi2S_CB2_VWG[0];
      paramCB2VWG[iTail].LowerLimits_psi2S = pp13_fit0_norm_psi2S_CB2_VWG[1];
      paramCB2VWG[iTail].UpperLimits_psi2S = pp13_fit0_norm_psi2S_CB2_VWG[2];

      //paramCB2VWG[iTail].paramsLimitVWG[4][0] = 0;
      //paramCB2VWG[iTail].paramsLimitVWG[4][1] = 1000;
 
    }
    if(iTail==0 && iRange == 1) // Tail: pp13TeV. Fit Range: 2.4 4.4
    {
      //paramCB2VWG[iTail].strNameOfFitFunc = Form("fitFuncCB2VWG_pp13_%i_%i", iTail, iRange);
      for(int iParam=0; iParam<numOfCB2_Parm; iParam++) paramCB2VWG[iTail].ParamJPsiFunc[iParam] = Param_doubleCB2_ppdata[iParam];

      paramCB2VWG[iTail].paramsInitialVWG[0]= pp13_fit1_par0_CB2_VWG[0];
      paramCB2VWG[iTail].paramsLimitVWG[0][0] = pp13_fit1_par0_CB2_VWG[1];
      paramCB2VWG[iTail].paramsLimitVWG[0][1] = pp13_fit1_par0_CB2_VWG[2];

      paramCB2VWG[iTail].paramsInitialVWG[1] = pp13_fit1_par1_CB2_VWG[0];
      paramCB2VWG[iTail].paramsLimitVWG[1][0] = pp13_fit1_par1_CB2_VWG[1];
      paramCB2VWG[iTail].paramsLimitVWG[1][1] = pp13_fit1_par1_CB2_VWG[2];

      paramCB2VWG[iTail].paramsInitialVWG[2] = pp13_fit1_par2_CB2_VWG[0];
      paramCB2VWG[iTail].paramsLimitVWG[2][0] = pp13_fit1_par2_CB2_VWG[1];
      paramCB2VWG[iTail].paramsLimitVWG[2][1] = pp13_fit1_par2_CB2_VWG[2];

      paramCB2VWG[iTail].paramsInitialVWG[3] = pp13_fit1_par3_CB2_VWG[0];
      paramCB2VWG[iTail].paramsLimitVWG[3][0] = pp13_fit1_par3_CB2_VWG[1];
      paramCB2VWG[iTail].paramsLimitVWG[3][1] = pp13_fit1_par3_CB2_VWG[2];

      paramCB2VWG[iTail].paramsInitialVWG[4] = pp13_fit1_par4_CB2_VWG[0];
      paramCB2VWG[iTail].paramsLimitVWG[4][0] = pp13_fit1_par4_CB2_VWG[1];
      paramCB2VWG[iTail].paramsLimitVWG[4][1] = pp13_fit1_par4_CB2_VWG[2];

      paramCB2VWG[iTail].paramsVarianceVWG[1] = pp13_fit1_par1_variance_CB2_VWG[0];
      paramCB2VWG[iTail].paramsLimitVarianceVWG[1][0] = pp13_fit1_par1_variance_CB2_VWG[1];
      paramCB2VWG[iTail].paramsLimitVarianceVWG[1][1] = pp13_fit1_par1_variance_CB2_VWG[2];

      paramCB2VWG[iTail].paramsVarianceVWG[2] = pp13_fit1_par2_variance_CB2_VWG[0];
      paramCB2VWG[iTail].paramsLimitVarianceVWG[2][0] = pp13_fit1_par2_variance_CB2_VWG[1];
      paramCB2VWG[iTail].paramsLimitVarianceVWG[2][1] = pp13_fit1_par2_variance_CB2_VWG[2];

      paramCB2VWG[iTail].paramsVarianceVWG[3] = pp13_fit1_par3_variance_CB2_VWG[0];
      paramCB2VWG[iTail].paramsLimitVarianceVWG[3][0] = pp13_fit1_par3_variance_CB2_VWG[1];
      paramCB2VWG[iTail].paramsLimitVarianceVWG[3][1] = pp13_fit1_par3_variance_CB2_VWG[2];

      paramCB2VWG[iTail].paramsVarianceVWG[4] = pp13_fit1_par4_variance_CB2_VWG[0];
      paramCB2VWG[iTail].paramsLimitVarianceVWG[4][0] = pp13_fit1_par4_variance_CB2_VWG[1];
      paramCB2VWG[iTail].paramsLimitVarianceVWG[4][1] = pp13_fit1_par4_variance_CB2_VWG[2];

      paramCB2VWG[iTail].mean_jpsi = pp13_fit1_mean_jpsi_CB2_VWG[0];
      paramCB2VWG[iTail].mean_LowerLimits_jpsi = pp13_fit1_mean_jpsi_CB2_VWG[1];
      paramCB2VWG[iTail].mean_UpperLimits_jpsi = pp13_fit1_mean_jpsi_CB2_VWG[2];

      paramCB2VWG[iTail].width_jpsi = pp13_fit1_width_jpsi_CB2_VWG[0];
      paramCB2VWG[iTail].width_LowerLimits_jpsi = pp13_fit1_width_jpsi_CB2_VWG[1];
      paramCB2VWG[iTail].width_UpperLimits_jpsi = pp13_fit1_width_jpsi_CB2_VWG[2];

      paramCB2VWG[iTail].Norm_psi2S = pp13_fit1_norm_psi2S_CB2_VWG[0];
      paramCB2VWG[iTail].LowerLimits_psi2S = pp13_fit1_norm_psi2S_CB2_VWG[1];
      paramCB2VWG[iTail].UpperLimits_psi2S = pp13_fit1_norm_psi2S_CB2_VWG[2];

    }

     if(iTail==1 && iRange==0) 
    {
       //paramCB2VWG[iTail].strNameOfFitFunc = Form("fitFuncCB2VWG_geant4_%i_%i", iTail, iRange);
       for(int iParam=0; iParam<numOfCB2_Parm; iParam++) paramCB2VWG[iTail].ParamJPsiFunc[iParam] = values_geant3_CB2_tail_range00[pt_bin][iParam];
 
       paramCB2VWG[iTail].paramsInitialVWG[0]= geant3_fit0_par0_CB2_VWG[0];
       paramCB2VWG[iTail].paramsLimitVWG[0][0] = geant3_fit0_par0_CB2_VWG[1];
       paramCB2VWG[iTail].paramsLimitVWG[0][1] = geant3_fit0_par0_CB2_VWG[2];

       paramCB2VWG[iTail].paramsInitialVWG[1] = geant3_fit0_par1_CB2_VWG[0];
       paramCB2VWG[iTail].paramsLimitVWG[1][0] = geant3_fit0_par1_CB2_VWG[1];
       paramCB2VWG[iTail].paramsLimitVWG[1][1] = geant3_fit0_par1_CB2_VWG[2];

       paramCB2VWG[iTail].paramsInitialVWG[2] = geant3_fit0_par2_CB2_VWG[0];
       paramCB2VWG[iTail].paramsLimitVWG[2][0] = geant3_fit0_par2_CB2_VWG[1];
       paramCB2VWG[iTail].paramsLimitVWG[2][1] = geant3_fit0_par2_CB2_VWG[2];

       paramCB2VWG[iTail].paramsInitialVWG[3] = geant3_fit0_par3_CB2_VWG[0];
       paramCB2VWG[iTail].paramsLimitVWG[3][0] = geant3_fit0_par3_CB2_VWG[1];
       paramCB2VWG[iTail].paramsLimitVWG[3][1] = geant3_fit0_par3_CB2_VWG[2];

       paramCB2VWG[iTail].paramsInitialVWG[4] = geant3_fit0_par4_CB2_VWG[0];
       paramCB2VWG[iTail].paramsLimitVWG[4][0] = geant3_fit0_par4_CB2_VWG[1];
       paramCB2VWG[iTail].paramsLimitVWG[4][1] = geant3_fit0_par4_CB2_VWG[2];

       paramCB2VWG[iTail].mean_jpsi = geant3_fit0_mean_jpsi_CB2_VWG[0];
       paramCB2VWG[iTail].mean_LowerLimits_jpsi = geant3_fit0_mean_jpsi_CB2_VWG[1];
       paramCB2VWG[iTail].mean_UpperLimits_jpsi = geant3_fit0_mean_jpsi_CB2_VWG[2];

       paramCB2VWG[iTail].width_jpsi = geant3_fit0_width_jpsi_CB2_VWG[0];
       paramCB2VWG[iTail].width_LowerLimits_jpsi = geant3_fit0_width_jpsi_CB2_VWG[1];
       paramCB2VWG[iTail].width_UpperLimits_jpsi = geant3_fit0_width_jpsi_CB2_VWG[2];

       paramCB2VWG[iTail].Norm_psi2S = geant3_fit0_norm_psi2S_CB2_VWG[0];
       paramCB2VWG[iTail].LowerLimits_psi2S = geant3_fit0_norm_psi2S_CB2_VWG[1];
       paramCB2VWG[iTail].UpperLimits_psi2S = geant3_fit0_norm_psi2S_CB2_VWG[2];

     }
     if(iTail==1 && iRange==1) 
     {
       //paramCB2VWG[iTail].strNameOfFitFunc = Form("fitFuncCB2VWG_geant4_%i_%i", iTail, iRange);
       for(int iParam=0; iParam<numOfCB2_Parm; iParam++) paramCB2VWG[iTail].ParamJPsiFunc[iParam] = values_geant3_CB2_tail_range01[pt_bin][iParam];

       paramCB2VWG[iTail].paramsInitialVWG[0]= geant3_fit1_par0_CB2_VWG[0];
       paramCB2VWG[iTail].paramsLimitVWG[0][0] = geant3_fit1_par0_CB2_VWG[1];
       paramCB2VWG[iTail].paramsLimitVWG[0][1] = geant3_fit1_par0_CB2_VWG[2];

       paramCB2VWG[iTail].paramsInitialVWG[1] = geant3_fit1_par1_CB2_VWG[0];      
       paramCB2VWG[iTail].paramsLimitVWG[1][0] = geant3_fit1_par1_CB2_VWG[1];  
       paramCB2VWG[iTail].paramsLimitVWG[1][1] = geant3_fit1_par1_CB2_VWG[2];    

       paramCB2VWG[iTail].paramsInitialVWG[2] = geant3_fit1_par2_CB2_VWG[0];      
       paramCB2VWG[iTail].paramsLimitVWG[2][0] = geant3_fit1_par2_CB2_VWG[1];   
       paramCB2VWG[iTail].paramsLimitVWG[2][1] = geant3_fit1_par2_CB2_VWG[2];    

       paramCB2VWG[iTail].paramsInitialVWG[3] = geant3_fit1_par3_CB2_VWG[0];      
       paramCB2VWG[iTail].paramsLimitVWG[3][0] = geant3_fit1_par3_CB2_VWG[1];   
       paramCB2VWG[iTail].paramsLimitVWG[3][1] = geant3_fit1_par3_CB2_VWG[2];    
 
       paramCB2VWG[iTail].paramsInitialVWG[4] = geant3_fit1_par4_CB2_VWG[0];
       paramCB2VWG[iTail].paramsLimitVWG[4][0] = geant3_fit1_par4_CB2_VWG[1];
       paramCB2VWG[iTail].paramsLimitVWG[4][1] = geant3_fit1_par4_CB2_VWG[2];

       paramCB2VWG[iTail].paramsVarianceVWG[1] = geant3_fit1_par1_variance_CB2_VWG[0];// 5* (1.01616e+03);
       paramCB2VWG[iTail].paramsLimitVarianceVWG[1][0] = geant3_fit1_par1_variance_CB2_VWG[1];
       paramCB2VWG[iTail].paramsLimitVarianceVWG[1][1] = geant3_fit1_par1_variance_CB2_VWG[2];

       paramCB2VWG[iTail].paramsVarianceVWG[2] = geant3_fit1_par2_variance_CB2_VWG[0]; 
       paramCB2VWG[iTail].paramsLimitVarianceVWG[2][0] = geant3_fit1_par2_variance_CB2_VWG[1];
       paramCB2VWG[iTail].paramsLimitVarianceVWG[2][1] = geant3_fit1_par2_variance_CB2_VWG[2];

       paramCB2VWG[iTail].paramsVarianceVWG[3] = geant3_fit1_par3_variance_CB2_VWG[0];
       paramCB2VWG[iTail].paramsLimitVarianceVWG[3][0] = geant3_fit1_par3_variance_CB2_VWG[1];
       paramCB2VWG[iTail].paramsLimitVarianceVWG[3][1] = geant3_fit1_par3_variance_CB2_VWG[2];

       paramCB2VWG[iTail].mean_jpsi = geant3_fit1_mean_jpsi_CB2_VWG[0];
       paramCB2VWG[iTail].mean_LowerLimits_jpsi = geant3_fit1_mean_jpsi_CB2_VWG[1];
       paramCB2VWG[iTail].mean_UpperLimits_jpsi = geant3_fit1_mean_jpsi_CB2_VWG[2];

       paramCB2VWG[iTail].width_jpsi = geant3_fit1_width_jpsi_CB2_VWG[0];
       paramCB2VWG[iTail].width_LowerLimits_jpsi = geant3_fit1_width_jpsi_CB2_VWG[1];
       paramCB2VWG[iTail].width_UpperLimits_jpsi = geant3_fit1_width_jpsi_CB2_VWG[2];

       paramCB2VWG[iTail].Norm_psi2S = geant3_fit1_norm_psi2S_CB2_VWG[0];
       paramCB2VWG[iTail].LowerLimits_psi2S = geant3_fit1_norm_psi2S_CB2_VWG[1];
       paramCB2VWG[iTail].UpperLimits_psi2S = geant3_fit1_norm_psi2S_CB2_VWG[2];
 
     }
/*
     if(iTail==2 && iRange==0) 
     {
       //paramCB2VWG[iTail].strNameOfFitFunc = Form("fitFuncCB2VWG_geant3_%i_%i", iTail, iRange);
       for(int iParam=0; iParam<numOfCB2_Parm; iParam++) paramCB2VWG[iTail].ParamJPsiFunc[iParam] = values_geant3_CB2_tail_range00[pt_bin][iParam];
       
       paramCB2VWG[iTail].paramsInitialVWG1[0] = geant3_fit0_par0_CB2_VWG1[0];
       paramCB2VWG[iTail].paramsLimitVWG1[0][0] = geant3_fit0_par0_CB2_VWG1[1];
       paramCB2VWG[iTail].paramsLimitVWG1[0][1] = geant3_fit0_par0_CB2_VWG1[2];

       paramCB2VWG[iTail].paramsInitialVWG1[1] = geant3_fit0_par1_CB2_VWG1[0];         
       paramCB2VWG[iTail].paramsLimitVWG1[1][0] = geant3_fit0_par1_CB2_VWG1[1];
       paramCB2VWG[iTail].paramsLimitVWG1[1][1] = geant3_fit0_par1_CB2_VWG1[2];

       paramCB2VWG[iTail].paramsInitialVWG1[2] = geant3_fit0_par2_CB2_VWG1[0]; 
       paramCB2VWG[iTail].paramsLimitVWG1[2][0] = geant3_fit0_par2_CB2_VWG1[1]; 
       paramCB2VWG[iTail].paramsLimitVWG1[2][1] = geant3_fit0_par2_CB2_VWG1[2];

       paramCB2VWG[iTail].paramsInitialVWG1[3] = geant3_fit0_par3_CB2_VWG1[0];
       paramCB2VWG[iTail].paramsLimitVWG1[3][0] = geant3_fit0_par3_CB2_VWG1[1]; 
       paramCB2VWG[iTail].paramsLimitVWG1[3][1] =  geant3_fit0_par3_CB2_VWG1[2];
       //paramCB2VWG[iTail].paramsLimitVWG[4][0] = 0;
       //paramCB2VWG[iTail].paramsLimitVWG[4][1] = 100;

       paramCB2VWG[iTail].paramsVarianceVWG1[1] = geant3_fit0_par1_variance_CB2_VWG1[0];
       paramCB2VWG[iTail].paramsLimitVarianceVWG1[1][0] = geant3_fit0_par1_variance_CB2_VWG1[1];
       paramCB2VWG[iTail].paramsLimitVarianceVWG1[1][1] = geant3_fit0_par1_variance_CB2_VWG1[2];

       paramCB2VWG[iTail].paramsVarianceVWG1[2] = geant3_fit0_par2_variance_CB2_VWG1[0];
       paramCB2VWG[iTail].paramsLimitVarianceVWG1[2][0]=geant3_fit0_par2_variance_CB2_VWG1[1];
       paramCB2VWG[iTail].paramsLimitVarianceVWG1[2][1]=geant3_fit0_par2_variance_CB2_VWG1[2];

       paramCB2VWG[iTail].paramsVarianceVWG1[3] = geant3_fit0_par3_variance_CB2_VWG1[0];
       paramCB2VWG[iTail].paramsLimitVarianceVWG1[3][0]=geant3_fit0_par3_variance_CB2_VWG1[1];
       paramCB2VWG[iTail].paramsLimitVarianceVWG1[3][1]=geant3_fit0_par3_variance_CB2_VWG1[2];

       paramCB2VWG[iTail].mean_jpsi = geant3_fit0_mean_jpsi_CB2_VWG1[0];
       paramCB2VWG[iTail].mean_LowerLimits_jpsi = geant3_fit0_mean_jpsi_CB2_VWG1[1];
       paramCB2VWG[iTail].mean_UpperLimits_jpsi = geant3_fit0_mean_jpsi_CB2_VWG1[2];

       paramCB2VWG[iTail].width_jpsi = geant3_fit0_width_jpsi_CB2_VWG1[0];
       paramCB2VWG[iTail].width_LowerLimits_jpsi = geant3_fit0_width_jpsi_CB2_VWG1[1];
       paramCB2VWG[iTail].width_UpperLimits_jpsi = geant3_fit0_width_jpsi_CB2_VWG1[2];

       paramCB2VWG[iTail].Norm_psi2S = geant3_fit0_norm_psi2S_CB2_VWG1[0];
       paramCB2VWG[iTail].LowerLimits_psi2S = geant3_fit0_norm_psi2S_CB2_VWG1[1];
       paramCB2VWG[iTail].UpperLimits_psi2S = geant3_fit0_norm_psi2S_CB2_VWG1[2];

     }

     if(iTail==2 && iRange==1) // Tail: geant3. Fit Range: 2.2 4.4
     {
       //paramCB2VWG[iTail].strNameOfFitFunc = Form("fitFuncCB2VWG_geant3_%i_%i", iTail, iRange);
       for(int iParam=0; iParam<numOfCB2_Parm; iParam++) paramCB2VWG[iTail].ParamJPsiFunc[iParam] = values_geant3_CB2_tail_range01[pt_bin][iParam];

       paramCB2VWG[iTail].paramsInitialVWG1[0]= geant3_fit1_par0_CB2_VWG1[0];
       paramCB2VWG[iTail].paramsLimitVWG1[0][0] = geant3_fit1_par0_CB2_VWG1[1];
       paramCB2VWG[iTail].paramsLimitVWG1[0][1] = geant3_fit1_par0_CB2_VWG1[2];
       
       paramCB2VWG[iTail].paramsInitialVWG1[1] = geant3_fit1_par1_CB2_VWG1[0]; 
       paramCB2VWG[iTail].paramsLimitVWG1[1][0] = geant3_fit1_par1_CB2_VWG1[1];          
       paramCB2VWG[iTail].paramsLimitVWG1[1][1] = geant3_fit1_par1_CB2_VWG1[2];     

       paramCB2VWG[iTail].paramsInitialVWG1[2] = geant3_fit1_par2_CB2_VWG1[0]; 
       paramCB2VWG[iTail].paramsLimitVWG1[2][0] = geant3_fit1_par2_CB2_VWG1[1];         
       paramCB2VWG[iTail].paramsLimitVWG1[2][1] = geant3_fit1_par2_CB2_VWG1[2];      

       paramCB2VWG[iTail].paramsInitialVWG1[3] = geant3_fit1_par3_CB2_VWG1[0];         
       paramCB2VWG[iTail].paramsLimitVWG1[3][0] = geant3_fit1_par3_CB2_VWG1[1];        
       paramCB2VWG[iTail].paramsLimitVWG1[3][1] = geant3_fit1_par3_CB2_VWG1[2];       
       //paramCB2VWG[iTail].paramsLimitVWG[4][0] = 0;
       //paramCB2VWG[iTail].paramsLimitVWG[4][1] = 1000;

       paramCB2VWG[iTail].paramsVarianceVWG1[1] = geant3_fit1_par1_variance_CB2_VWG1[0];
       paramCB2VWG[iTail].paramsLimitVarianceVWG1[1][0]=geant3_fit1_par1_variance_CB2_VWG1[1];
       paramCB2VWG[iTail].paramsLimitVarianceVWG1[1][1]=geant3_fit1_par1_variance_CB2_VWG1[2];

       paramCB2VWG[iTail].paramsVarianceVWG1[2] = geant3_fit1_par2_variance_CB2_VWG1[0];
       paramCB2VWG[iTail].paramsLimitVarianceVWG1[2][0]=geant3_fit1_par2_variance_CB2_VWG1[1];
       paramCB2VWG[iTail].paramsLimitVarianceVWG1[2][1]=geant3_fit1_par2_variance_CB2_VWG1[2]; 

       paramCB2VWG[iTail].paramsVarianceVWG1[3] = geant3_fit1_par3_variance_CB2_VWG1[0];
       paramCB2VWG[iTail].paramsLimitVarianceVWG1[3][0]= geant3_fit1_par3_variance_CB2_VWG1[1];
       paramCB2VWG[iTail].paramsLimitVarianceVWG1[3][1]= geant3_fit1_par3_variance_CB2_VWG1[2];

       paramCB2VWG[iTail].mean_jpsi = geant3_fit1_mean_jpsi_CB2_VWG1[0];
       paramCB2VWG[iTail].mean_LowerLimits_jpsi = geant3_fit1_mean_jpsi_CB2_VWG1[1];
       paramCB2VWG[iTail].mean_UpperLimits_jpsi = geant3_fit1_mean_jpsi_CB2_VWG1[2];

       paramCB2VWG[iTail].width_jpsi = geant3_fit1_width_jpsi_CB2_VWG1[0];
       paramCB2VWG[iTail].width_LowerLimits_jpsi = geant3_fit1_width_jpsi_CB2_VWG1[1];
       paramCB2VWG[iTail].width_UpperLimits_jpsi = geant3_fit1_width_jpsi_CB2_VWG1[2];

       paramCB2VWG[iTail].Norm_psi2S = geant3_fit1_norm_psi2S_CB2_VWG1[0];
       paramCB2VWG[iTail].LowerLimits_psi2S = geant3_fit1_norm_psi2S_CB2_VWG1[1];
       paramCB2VWG[iTail].UpperLimits_psi2S = geant3_fit1_norm_psi2S_CB2_VWG1[2];

     }
*/
      for (int iCent=0; iCent < 1; iCent++)
      {

        TF1 *fitFuncCB2VWG = new TF1(Form("%s",paramCB2VWG[iTail].strNameOfFitFunc.Data()), fitFunctionDoubleCB2VWG, FitRange[iRange][0], FitRange[iRange][1], numOfVWGCb2_Param+1);

        TF1 *fitVWGreject = new TF1("fitVWGreject", BackgroundVWGreject_JPsi_Psi2S, FitRange[iRange][0], FitRange[iRange][1], numOfVWG_Parm); 
 
        TF1 *fitVWG = new TF1("fitVWG", BackgroundVWG, FitRange[iRange][0], FitRange[iRange][1], numOfVWG_Parm);

        TF1 *fitCB2_JPsi = new TF1("fitCB2_JPsi", CrystalBallExtended, FitRange[iRange][0], FitRange[iRange][1], numOfCB2_Parm);

        TF1 *fitCB2_Psi2S = new TF1("fitCB2_Psi2S", CrystalBallExtended, FitRange[iRange][0], FitRange[iRange][1], numOfCB2_Parm);


        outcome_fit = fitJPsi_Psi2S( iCent,
                      FitRange[iRange][0],
                      FitRange[iRange][1], 
//                      histoRawInvM,
                      histoInvMass_double_diff[pt_bin][y_bin],
//                      histoInvMass[pt_bin],
//                      histoInvMass_Pt12to15, 
                      fitVWGreject,
                      fitFuncCB2VWG,
                      fitVWG,
                      fitCB2_JPsi,
                      fitCB2_Psi2S,
                      paramCB2VWG[iTail],
                      Njpsi_CB2VWG[iRange][iTail],
                      SigmaJpsi_CB2VWG[iRange][iTail],
                      MassJpsi_CB2VWG[iRange][iTail],
                      ChiJpsi_CB2VWG[iRange][iTail],
                      ErrorOfNjpsi_CB2VWG[iRange][iTail],
                      ErrorOfSigmaJpsi_CB2VWG[iRange][iTail],
                      ErrorOfMassJpsi_CB2VWG[iRange][iTail],
                      ErrorOfChiJpsi_CB2VWG[iRange][iTail],
                      fitCanvas_CB2background[iCount1],
//                      fitCanvas_CB2,
//                      fitCanvas_background[iCount1],
                      outcome_fit);

                      std::cout << "num jpsi: " << outcome_fit.Njpsi_CB2VWG << std::endl;
                      std::cout << "fit validation: " << outcome_fit.fit_valide << std::endl;
                      std::cout << "fit status: " << outcome_fit.fit_status << std::endl;
                      std::cout << "fit covariant matrix: " << outcome_fit.fit_cov_matrix_status << std::endl;

          fit_valide_CB2VWG[iRange][iTail][0] = outcome_fit.fit_valide;
          fit_status_CB2VWG[iRange][iTail][0] = outcome_fit.fit_status;
          fit_cov_matrix_status_CB2VWG[iRange][iTail][0] = outcome_fit.fit_cov_matrix_status;

//        histoRawInvM_copy[iCount1]->Draw();
//        fitCB2_JPsi->Draw("SAME");
//        fitCB2_Psi2S->Draw("SAME");
//        fitVWG->Draw("SAME");
 //       fitFuncCB2VWG->Draw("SAME");

        iCount1 +=1;
//        std::cout << "iCount1: " << iCount1 << std::endl;

//        delete fitVWGreject;
//        delete fitFuncCB2VWG;
 //       delete fitVWG;
//        delete fitCB2;


      } // end of for (int iCent=0; iCent < 1; iCent++)  
    } // end of for (int iTail=0; iTail< NUMOFTAILPARAM_CB2; iTail++)
  }  // end of for (iRange=0; iRange < 2; iRange++ )

 
  for (int iRange=0; iRange < FITRANGE; iRange++ )
//    for (int iRange=0; iRange < 1; iRange++ )
   {
    for (int iTail=0; iTail< NUMOFTAILPARAM_CB2; iTail++)
//    for (int iTail=0; iTail< 1; iTail++)
    {     
//      paramCB2VWG[iTail].fitRange[iRange][0]=FitRange[iRange][0]; 
//      paramCB2VWG[iTail].fitRange[iRange][1]=FitRange[iRange][1];
     if(iTail==0 && iRange == 0)        //CB2 pol pp13 2.2 4.5
     {
       for(int iParam=0; iParam<numOfCB2_Parm; iParam++) paramCB2Pol2[iTail].ParamJPsiFunc[iParam] = Param_doubleCB2_ppdata[iParam];
       paramCB2Pol2[iTail].paramsInitialPol2overPol3[0] = pp13_fit0_par0_CB2_pol[0];
       paramCB2Pol2[iTail].paramsLimitPol2overPol3[0][0] = pp13_fit0_par0_CB2_pol[1];
       paramCB2Pol2[iTail].paramsLimitPol2overPol3[0][1] = pp13_fit0_par0_CB2_pol[2];

       paramCB2Pol2[iTail].paramsInitialPol2overPol3[1] = pp13_fit0_par1_CB2_pol[0];  
       paramCB2Pol2[iTail].paramsLimitPol2overPol3[1][0] = pp13_fit0_par1_CB2_pol[1];     
       paramCB2Pol2[iTail].paramsLimitPol2overPol3[1][1] =  pp13_fit0_par1_CB2_pol[2];     

       paramCB2Pol2[iTail].paramsInitialPol2overPol3[2] = pp13_fit0_par2_CB2_pol[0];         
       paramCB2Pol2[iTail].paramsLimitPol2overPol3[2][0] = pp13_fit0_par2_CB2_pol[1];   
       paramCB2Pol2[iTail].paramsLimitPol2overPol3[2][1] =  pp13_fit0_par2_CB2_pol[2];   

       paramCB2Pol2[iTail].paramsInitialPol2overPol3[3] = pp13_fit0_par3_CB2_pol[0];         
       paramCB2Pol2[iTail].paramsLimitPol2overPol3[3][0] = pp13_fit0_par3_CB2_pol[1];  
       paramCB2Pol2[iTail].paramsLimitPol2overPol3[3][1] =  pp13_fit0_par3_CB2_pol[2];  

       paramCB2Pol2[iTail].paramsInitialPol2overPol3[4] = pp13_fit0_par4_CB2_pol[0];         
       paramCB2Pol2[iTail].paramsLimitPol2overPol3[4][0] = pp13_fit0_par4_CB2_pol[1];   
       paramCB2Pol2[iTail].paramsLimitPol2overPol3[4][1] =  pp13_fit0_par4_CB2_pol[2];   

       paramCB2Pol2[iTail].paramsInitialPol2overPol3[5] = pp13_fit0_par5_CB2_pol[0];         
       paramCB2Pol2[iTail].paramsLimitPol2overPol3[5][0] = pp13_fit0_par5_CB2_pol[1];
       paramCB2Pol2[iTail].paramsLimitPol2overPol3[5][1] = pp13_fit0_par5_CB2_pol[2];

       paramCB2Pol2[iTail].paramsInitialPol2overPol3[6] = pp13_fit0_par6_CB2_pol[0];         
       paramCB2Pol2[iTail].paramsLimitPol2overPol3[6][0] = pp13_fit0_par6_CB2_pol[1];
       paramCB2Pol2[iTail].paramsLimitPol2overPol3[6][1] = pp13_fit0_par6_CB2_pol[2];

       paramCB2Pol2[iTail].paramsVariancePol2overPol3[1] = pp13_fit0_par1_variance_CB2_pol[0];
       paramCB2Pol2[iTail].paramsLimitVariancePol2overPol3[1][0] = pp13_fit0_par1_variance_CB2_pol[1];
       paramCB2Pol2[iTail].paramsLimitVariancePol2overPol3[1][1] = pp13_fit0_par1_variance_CB2_pol[2]; 

       paramCB2Pol2[iTail].paramsVariancePol2overPol3[2] = pp13_fit0_par2_variance_CB2_pol[0];
       paramCB2Pol2[iTail].paramsLimitVariancePol2overPol3[2][0] = pp13_fit0_par2_variance_CB2_pol[1];
       paramCB2Pol2[iTail].paramsLimitVariancePol2overPol3[2][1] = pp13_fit0_par2_variance_CB2_pol[2];

       paramCB2Pol2[iTail].paramsVariancePol2overPol3[3] = pp13_fit0_par3_variance_CB2_pol[0];
       paramCB2Pol2[iTail].paramsLimitVariancePol2overPol3[3][0] = pp13_fit0_par3_variance_CB2_pol[1];
       paramCB2Pol2[iTail].paramsLimitVariancePol2overPol3[3][1] = pp13_fit0_par3_variance_CB2_pol[2];

       paramCB2Pol2[iTail].paramsVariancePol2overPol3[4] = pp13_fit0_par4_variance_CB2_pol[0];
       paramCB2Pol2[iTail].paramsLimitVariancePol2overPol3[4][0] = pp13_fit0_par4_variance_CB2_pol[1];
       paramCB2Pol2[iTail].paramsLimitVariancePol2overPol3[4][1] = pp13_fit0_par4_variance_CB2_pol[2];

       paramCB2Pol2[iTail].paramsVariancePol2overPol3[5] = pp13_fit0_par5_variance_CB2_pol[0];
       paramCB2Pol2[iTail].paramsLimitVariancePol2overPol3[5][0] = pp13_fit0_par5_variance_CB2_pol[1];
       paramCB2Pol2[iTail].paramsLimitVariancePol2overPol3[5][1] = pp13_fit0_par5_variance_CB2_pol[2];

       paramCB2Pol2[iTail].paramsVariancePol2overPol3[6] = pp13_fit0_par6_variance_CB2_pol[0];
       paramCB2Pol2[iTail].paramsLimitVariancePol2overPol3[6][0] = pp13_fit0_par6_variance_CB2_pol[1];
       paramCB2Pol2[iTail].paramsLimitVariancePol2overPol3[6][1] = pp13_fit0_par6_variance_CB2_pol[2];

       paramCB2Pol2[iTail].mean_jpsi = pp13_fit0_mean_jpsi_CB2_pol[0];
       paramCB2Pol2[iTail].mean_LowerLimits_jpsi = pp13_fit0_mean_jpsi_CB2_pol[1];
       paramCB2Pol2[iTail].mean_UpperLimits_jpsi = pp13_fit0_mean_jpsi_CB2_pol[2];
 
       paramCB2Pol2[iTail].width_jpsi = pp13_fit0_width_jpsi_CB2_pol[0];
       paramCB2Pol2[iTail].width_LowerLimits_jpsi = pp13_fit0_width_jpsi_CB2_pol[1];
       paramCB2Pol2[iTail].width_UpperLimits_jpsi = pp13_fit0_width_jpsi_CB2_pol[2];

       paramCB2Pol2[iTail].Norm_psi2S = pp13_fit0_norm_psi2S_CB2_pol[0];
       paramCB2Pol2[iTail].LowerLimits_psi2S = pp13_fit0_norm_psi2S_CB2_pol[1];
       paramCB2Pol2[iTail].UpperLimits_psi2S = pp13_fit0_norm_psi2S_CB2_pol[2];       
     }
     if(iTail==0 && iRange == 1)        //CB2 pol pp13 2.4 4.7
     {
       for(int iParam=0; iParam<numOfCB2_Parm; iParam++) paramCB2Pol2[iTail].ParamJPsiFunc[iParam] = Param_doubleCB2_ppdata[iParam];
       paramCB2Pol2[iTail].paramsInitialPol2overPol3[0] = pp13_fit1_par0_CB2_pol[0];
       paramCB2Pol2[iTail].paramsLimitPol2overPol3[0][0] = pp13_fit1_par0_CB2_pol[1];
       paramCB2Pol2[iTail].paramsLimitPol2overPol3[0][1] = pp13_fit1_par0_CB2_pol[2];

       paramCB2Pol2[iTail].paramsInitialPol2overPol3[1] = pp13_fit1_par1_CB2_pol[0];
       paramCB2Pol2[iTail].paramsLimitPol2overPol3[1][0] = pp13_fit1_par1_CB2_pol[1];
       paramCB2Pol2[iTail].paramsLimitPol2overPol3[1][1] =  pp13_fit1_par1_CB2_pol[2];

       paramCB2Pol2[iTail].paramsInitialPol2overPol3[2] = pp13_fit1_par2_CB2_pol[0];
       paramCB2Pol2[iTail].paramsLimitPol2overPol3[2][0] = pp13_fit1_par2_CB2_pol[1];
       paramCB2Pol2[iTail].paramsLimitPol2overPol3[2][1] =  pp13_fit1_par2_CB2_pol[2];

       paramCB2Pol2[iTail].paramsInitialPol2overPol3[3] = pp13_fit1_par3_CB2_pol[0];
       paramCB2Pol2[iTail].paramsLimitPol2overPol3[3][0] = pp13_fit1_par3_CB2_pol[1];
       paramCB2Pol2[iTail].paramsLimitPol2overPol3[3][1] =  pp13_fit1_par3_CB2_pol[2];

       paramCB2Pol2[iTail].paramsInitialPol2overPol3[4] = pp13_fit1_par4_CB2_pol[0];
       paramCB2Pol2[iTail].paramsLimitPol2overPol3[4][0] = pp13_fit1_par4_CB2_pol[1];
       paramCB2Pol2[iTail].paramsLimitPol2overPol3[4][1] =  pp13_fit1_par4_CB2_pol[2];

       paramCB2Pol2[iTail].paramsInitialPol2overPol3[5] = pp13_fit1_par5_CB2_pol[0];
       paramCB2Pol2[iTail].paramsLimitPol2overPol3[5][0] = pp13_fit1_par5_CB2_pol[1];
       paramCB2Pol2[iTail].paramsLimitPol2overPol3[5][1] =  pp13_fit1_par5_CB2_pol[2];

       paramCB2Pol2[iTail].paramsInitialPol2overPol3[6] = pp13_fit1_par6_CB2_pol[0];
       paramCB2Pol2[iTail].paramsLimitPol2overPol3[6][0] = pp13_fit1_par6_CB2_pol[1];
       paramCB2Pol2[iTail].paramsLimitPol2overPol3[6][1] =  pp13_fit1_par6_CB2_pol[2];

       for(int iParam=1; iParam < numOfPol2overPol3_Parm; iParam++)
       {
           paramCB2Pol2[iTail].paramsVariancePol2overPol3[iParam] = 0;
           paramCB2Pol2[iTail].paramsLimitVariancePol2overPol3[iParam][0] = 0;
           paramCB2Pol2[iTail].paramsLimitVariancePol2overPol3[iParam][1] = 0;
       }
       //paramCB2Pol2[iTail].paramsLimitPol2overPol3[5][0] = -1000;
       //paramCB2Pol2[iTail].paramsLimitPol2overPol3[5][1] =  1000;
       //paramCB2Pol2[iTail].paramsLimitPol2overPol3[6][0] = -100;
       //paramCB2Pol2[iTail].paramsLimitPol2overPol3[6][1] =  100;

       paramCB2Pol2[iTail].mean_jpsi = pp13_fit1_mean_jpsi_CB2_pol[0];
       paramCB2Pol2[iTail].mean_LowerLimits_jpsi = pp13_fit1_mean_jpsi_CB2_pol[1];
       paramCB2Pol2[iTail].mean_UpperLimits_jpsi = pp13_fit1_mean_jpsi_CB2_pol[2];
 
       paramCB2Pol2[iTail].width_jpsi = pp13_fit1_width_jpsi_CB2_pol[0];
       paramCB2Pol2[iTail].width_LowerLimits_jpsi = pp13_fit1_width_jpsi_CB2_pol[1];
       paramCB2Pol2[iTail].width_UpperLimits_jpsi = pp13_fit1_width_jpsi_CB2_pol[2];

       paramCB2Pol2[iTail].Norm_psi2S = pp13_fit1_norm_psi2S_CB2_pol[0];
       paramCB2Pol2[iTail].LowerLimits_psi2S = pp13_fit1_norm_psi2S_CB2_pol[1];
       paramCB2Pol2[iTail].UpperLimits_psi2S = pp13_fit1_norm_psi2S_CB2_pol[2];  
     }

     if(iTail==1 && iRange==0)        
     {
       for(int iParam=0; iParam<numOfCB2_Parm; iParam++) paramCB2Pol2[iTail].ParamJPsiFunc[iParam] = values_geant3_CB2_tail_range00[pt_bin][iParam];

       paramCB2Pol2[iTail].paramsInitialPol2overPol3[0] = geant3_fit0_par0_CB2_pol[0];
       paramCB2Pol2[iTail].paramsLimitPol2overPol3[0][0] = geant3_fit0_par0_CB2_pol[1];
       paramCB2Pol2[iTail].paramsLimitPol2overPol3[0][1] = geant3_fit0_par0_CB2_pol[2];

       paramCB2Pol2[iTail].paramsInitialPol2overPol3[1] = geant3_fit0_par1_CB2_pol[0];
       paramCB2Pol2[iTail].paramsLimitPol2overPol3[1][0] = geant3_fit0_par1_CB2_pol[1];
       paramCB2Pol2[iTail].paramsLimitPol2overPol3[1][1] =  geant3_fit0_par1_CB2_pol[2];

       paramCB2Pol2[iTail].paramsInitialPol2overPol3[2] = geant3_fit0_par2_CB2_pol[0];
       paramCB2Pol2[iTail].paramsLimitPol2overPol3[2][0] = geant3_fit0_par2_CB2_pol[1];
       paramCB2Pol2[iTail].paramsLimitPol2overPol3[2][1] =  geant3_fit0_par2_CB2_pol[2];

       paramCB2Pol2[iTail].paramsInitialPol2overPol3[3] = geant3_fit0_par3_CB2_pol[0];
       paramCB2Pol2[iTail].paramsLimitPol2overPol3[3][0] = geant3_fit0_par3_CB2_pol[1];
       paramCB2Pol2[iTail].paramsLimitPol2overPol3[3][1] =  geant3_fit0_par3_CB2_pol[2];

       paramCB2Pol2[iTail].paramsInitialPol2overPol3[4] = geant3_fit0_par4_CB2_pol[0];
       paramCB2Pol2[iTail].paramsLimitPol2overPol3[4][0] = geant3_fit0_par4_CB2_pol[1];
       paramCB2Pol2[iTail].paramsLimitPol2overPol3[4][1] =  geant3_fit0_par4_CB2_pol[2];

       paramCB2Pol2[iTail].paramsInitialPol2overPol3[5] = geant3_fit0_par5_CB2_pol[0];
       paramCB2Pol2[iTail].paramsLimitPol2overPol3[5][0] = geant3_fit0_par5_CB2_pol[1];
       paramCB2Pol2[iTail].paramsLimitPol2overPol3[5][1] =  geant3_fit0_par5_CB2_pol[2];

       paramCB2Pol2[iTail].paramsInitialPol2overPol3[6] = geant3_fit0_par6_CB2_pol[0];
       paramCB2Pol2[iTail].paramsLimitPol2overPol3[6][0] = geant3_fit0_par6_CB2_pol[1];
       paramCB2Pol2[iTail].paramsLimitPol2overPol3[6][1] =  geant3_fit0_par6_CB2_pol[2];

       for(int iParam=1; iParam < numOfPol2overPol3_Parm; iParam++)
       {
           paramCB2Pol2[iTail].paramsVariancePol2overPol3[iParam] = 0;
           paramCB2Pol2[iTail].paramsLimitVariancePol2overPol3[iParam][0] = 0;
           paramCB2Pol2[iTail].paramsLimitVariancePol2overPol3[iParam][1] = 0;
       }
       //paramCB2Pol2[iTail].paramsLimitPol2overPol3[5][0] = -10000;
       //paramCB2Pol2[iTail].paramsLimitPol2overPol3[5][1] =  10000;
       //paramCB2Pol2[iTail].paramsLimitPol2overPol3[6][0] = -1000;
       //paramCB2Pol2[iTail].paramsLimitPol2overPol3[6][1] =  1000;

       paramCB2Pol2[iTail].mean_jpsi = geant3_fit0_mean_jpsi_CB2_pol[0];
       paramCB2Pol2[iTail].mean_LowerLimits_jpsi = geant3_fit0_mean_jpsi_CB2_pol[1];
       paramCB2Pol2[iTail].mean_UpperLimits_jpsi = geant3_fit0_mean_jpsi_CB2_pol[2];
 
       paramCB2Pol2[iTail].width_jpsi = geant3_fit0_width_jpsi_CB2_pol[0];
       paramCB2Pol2[iTail].width_LowerLimits_jpsi = geant3_fit0_width_jpsi_CB2_pol[1];
       paramCB2Pol2[iTail].width_UpperLimits_jpsi = geant3_fit0_width_jpsi_CB2_pol[2];

       paramCB2Pol2[iTail].Norm_psi2S = geant3_fit0_norm_psi2S_CB2_pol[0];
       paramCB2Pol2[iTail].LowerLimits_psi2S = geant3_fit0_norm_psi2S_CB2_pol[1];
       paramCB2Pol2[iTail].UpperLimits_psi2S = geant3_fit0_norm_psi2S_CB2_pol[2]; 
     }
     if(iTail==1 && iRange==1)        //CB2 Pol geant4 2.4 4.7
     {
       for(int iParam=0; iParam<numOfCB2_Parm; iParam++) paramCB2Pol2[iTail].ParamJPsiFunc[iParam] = values_geant3_CB2_tail_range01[pt_bin][iParam];

       paramCB2Pol2[iTail].paramsInitialPol2overPol3[0] = geant3_fit1_par0_CB2_pol[0];
       paramCB2Pol2[iTail].paramsLimitPol2overPol3[0][0] = geant3_fit1_par0_CB2_pol[1];
       paramCB2Pol2[iTail].paramsLimitPol2overPol3[0][1] = geant3_fit1_par0_CB2_pol[2];

       paramCB2Pol2[iTail].paramsInitialPol2overPol3[1] = geant3_fit1_par1_CB2_pol[0];
       paramCB2Pol2[iTail].paramsLimitPol2overPol3[1][0] = geant3_fit1_par1_CB2_pol[1];
       paramCB2Pol2[iTail].paramsLimitPol2overPol3[1][1] =  geant3_fit1_par1_CB2_pol[2];

       paramCB2Pol2[iTail].paramsInitialPol2overPol3[2] = geant3_fit1_par2_CB2_pol[0];
       paramCB2Pol2[iTail].paramsLimitPol2overPol3[2][0] = geant3_fit1_par2_CB2_pol[1];
       paramCB2Pol2[iTail].paramsLimitPol2overPol3[2][1] =  geant3_fit1_par2_CB2_pol[2];

       paramCB2Pol2[iTail].paramsInitialPol2overPol3[3] = geant3_fit1_par3_CB2_pol[0];
       paramCB2Pol2[iTail].paramsLimitPol2overPol3[3][0] = geant3_fit1_par3_CB2_pol[1];
       paramCB2Pol2[iTail].paramsLimitPol2overPol3[3][1] =  geant3_fit1_par3_CB2_pol[2];

       paramCB2Pol2[iTail].paramsInitialPol2overPol3[4] = geant3_fit1_par4_CB2_pol[0];
       paramCB2Pol2[iTail].paramsLimitPol2overPol3[4][0] = geant3_fit1_par4_CB2_pol[1];
       paramCB2Pol2[iTail].paramsLimitPol2overPol3[4][1] =  geant3_fit1_par4_CB2_pol[2];

       paramCB2Pol2[iTail].paramsInitialPol2overPol3[5] = geant3_fit1_par5_CB2_pol[0];
       paramCB2Pol2[iTail].paramsLimitPol2overPol3[5][0] = geant3_fit1_par5_CB2_pol[1];
       paramCB2Pol2[iTail].paramsLimitPol2overPol3[5][1] =  geant3_fit1_par5_CB2_pol[2];

       paramCB2Pol2[iTail].paramsInitialPol2overPol3[6] = geant3_fit1_par6_CB2_pol[0];
       paramCB2Pol2[iTail].paramsLimitPol2overPol3[6][0] = geant3_fit1_par6_CB2_pol[1];
       paramCB2Pol2[iTail].paramsLimitPol2overPol3[6][1] =  geant3_fit1_par6_CB2_pol[2];

       for(int iParam=1; iParam < numOfPol2overPol3_Parm; iParam++)
       {
           paramCB2Pol2[iTail].paramsVariancePol2overPol3[iParam] = 0;
           paramCB2Pol2[iTail].paramsLimitVariancePol2overPol3[iParam][0] = 0;
           paramCB2Pol2[iTail].paramsLimitVariancePol2overPol3[iParam][1] = 0;
       }
       //paramCB2Pol2[iTail].paramsLimitPol2overPol3[5][0] = -10;
       //paramCB2Pol2[iTail].paramsLimitPol2overPol3[5][1] =  10;
       //paramCB2Pol2[iTail].paramsLimitPol2overPol3[6][0] = -10;
       //paramCB2Pol2[iTail].paramsLimitPol2overPol3[6][1] =  10;

       paramCB2Pol2[iTail].mean_jpsi = geant3_fit1_mean_jpsi_CB2_pol[0];
       paramCB2Pol2[iTail].mean_LowerLimits_jpsi = geant3_fit1_mean_jpsi_CB2_pol[1];
       paramCB2Pol2[iTail].mean_UpperLimits_jpsi = geant3_fit1_mean_jpsi_CB2_pol[2];
 
       paramCB2Pol2[iTail].width_jpsi = geant3_fit1_width_jpsi_CB2_pol[0];
       paramCB2Pol2[iTail].width_LowerLimits_jpsi = geant3_fit1_width_jpsi_CB2_pol[1];
       paramCB2Pol2[iTail].width_UpperLimits_jpsi = geant3_fit1_width_jpsi_CB2_pol[2];

       paramCB2Pol2[iTail].Norm_psi2S = geant3_fit1_norm_psi2S_CB2_pol[0];
       paramCB2Pol2[iTail].LowerLimits_psi2S = geant3_fit1_norm_psi2S_CB2_pol[1];
       paramCB2Pol2[iTail].UpperLimits_psi2S = geant3_fit1_norm_psi2S_CB2_pol[2]; 
     }

/*
     if(iTail==2 && iRange==0)        //CB2 Pol geant3 2.0 4.8
     {
       for(int iParam=0; iParam<numOfCB2_Parm; iParam++) paramCB2Pol2[iTail].ParamJPsiFunc[iParam] = values_geant3_CB2_tail_range00[pt_bin][iParam];

       paramCB2Pol2[iTail].paramsLimitPol1overPol2[0][0] = geant3_fit0_par0_CB2_pol[1];
       paramCB2Pol2[iTail].paramsLimitPol1overPol2[0][1] = geant3_fit0_par0_CB2_pol[2];

       paramCB2Pol2[iTail].paramsInitialPol1overPol2[1] = geant3_fit0_par1_CB2_pol[0];
       paramCB2Pol2[iTail].paramsLimitPol1overPol2[1][0] = geant3_fit0_par1_CB2_pol[1];
       paramCB2Pol2[iTail].paramsLimitPol1overPol2[1][1] =  geant3_fit0_par1_CB2_pol[2];

       paramCB2Pol2[iTail].paramsInitialPol1overPol2[2] = geant3_fit0_par2_CB2_pol[0];
       paramCB2Pol2[iTail].paramsLimitPol1overPol2[2][0] = geant3_fit0_par2_CB2_pol[1];
       paramCB2Pol2[iTail].paramsLimitPol1overPol2[2][1] =  geant3_fit0_par2_CB2_pol[2];

       paramCB2Pol2[iTail].paramsInitialPol1overPol2[3] = geant3_fit0_par3_CB2_pol[0];
       paramCB2Pol2[iTail].paramsLimitPol1overPol2[3][0] = geant3_fit0_par3_CB2_pol[1];
       paramCB2Pol2[iTail].paramsLimitPol1overPol2[3][1] =  geant3_fit0_par3_CB2_pol[2];

       paramCB2Pol2[iTail].paramsInitialPol1overPol2[4] = geant3_fit0_par4_CB2_pol[0];
       paramCB2Pol2[iTail].paramsLimitPol1overPol2[4][0] = geant3_fit0_par4_CB2_pol[1];
       paramCB2Pol2[iTail].paramsLimitPol1overPol2[4][1] =  geant3_fit0_par4_CB2_pol[2];

       for(int iParam=1; iParam < numOfPol1overPol2_Parm; iParam++)
       {
           paramCB2Pol2[iTail].paramsVariancePol1overPol2[iParam] = 0;
           paramCB2Pol2[iTail].paramsLimitVariancePol1overPol2[iParam][0] = 0;
           paramCB2Pol2[iTail].paramsLimitVariancePol1overPol2[iParam][1] = 0;
       }
       //paramCB2Pol2[iTail].paramsLimitPol2overPol3[5][0] = -1000;
       //paramCB2Pol2[iTail].paramsLimitPol2overPol3[5][1] =  1000;
       //paramCB2Pol2[iTail].paramsLimitPol2overPol3[6][0] = -1000;
       //paramCB2Pol2[iTail].paramsLimitPol2overPol3[6][1] =  1000;
       paramCB2Pol2[iTail].mean_jpsi = geant3_fit0_mean_jpsi_CB2_pol[0];
       paramCB2Pol2[iTail].mean_LowerLimits_jpsi = geant3_fit0_mean_jpsi_CB2_pol[1];
       paramCB2Pol2[iTail].mean_UpperLimits_jpsi = geant3_fit0_mean_jpsi_CB2_pol[2];
 
       paramCB2Pol2[iTail].width_jpsi = geant3_fit0_width_jpsi_CB2_pol[0];
       paramCB2Pol2[iTail].width_LowerLimits_jpsi = geant3_fit0_width_jpsi_CB2_pol[1];
       paramCB2Pol2[iTail].width_UpperLimits_jpsi = geant3_fit0_width_jpsi_CB2_pol[2];

       paramCB2Pol2[iTail].Norm_psi2S = geant3_fit0_norm_psi2S_CB2_pol[0];
       paramCB2Pol2[iTail].LowerLimits_psi2S = geant3_fit0_norm_psi2S_CB2_pol[1];
       paramCB2Pol2[iTail].UpperLimits_psi2S = geant3_fit0_norm_psi2S_CB2_pol[2];
     }

     if(iTail==2 && iRange==1)
     {

       for(int iParam=0; iParam<numOfCB2_Parm; iParam++) paramCB2Pol2[iTail].ParamJPsiFunc[iParam] = values_geant3_CB2_tail_range01[pt_bin][iParam];

       paramCB2Pol2[iTail].paramsLimitPol1overPol2[0][0] = geant3_fit1_par0_CB2_pol[1];
       paramCB2Pol2[iTail].paramsLimitPol1overPol2[0][1] = geant3_fit1_par0_CB2_pol[2];

       paramCB2Pol2[iTail].paramsInitialPol1overPol2[1] = geant3_fit1_par1_CB2_pol[0];
       paramCB2Pol2[iTail].paramsLimitPol1overPol2[1][0] = geant3_fit1_par1_CB2_pol[1];
       paramCB2Pol2[iTail].paramsLimitPol1overPol2[1][1] =  geant3_fit1_par1_CB2_pol[2];

       paramCB2Pol2[iTail].paramsInitialPol1overPol2[2] = geant3_fit1_par2_CB2_pol[0];
       paramCB2Pol2[iTail].paramsLimitPol1overPol2[2][0] = geant3_fit1_par2_CB2_pol[1];
       paramCB2Pol2[iTail].paramsLimitPol1overPol2[2][1] =  geant3_fit1_par2_CB2_pol[2];

       paramCB2Pol2[iTail].paramsInitialPol1overPol2[3] = geant3_fit1_par3_CB2_pol[0];
       paramCB2Pol2[iTail].paramsLimitPol1overPol2[3][0] = geant3_fit1_par3_CB2_pol[1];
       paramCB2Pol2[iTail].paramsLimitPol1overPol2[3][1] =  geant3_fit1_par3_CB2_pol[2];

       paramCB2Pol2[iTail].paramsInitialPol1overPol2[4] = geant3_fit1_par4_CB2_pol[0];
       paramCB2Pol2[iTail].paramsLimitPol1overPol2[4][0] = geant3_fit1_par4_CB2_pol[1];
       paramCB2Pol2[iTail].paramsLimitPol1overPol2[4][1] =  geant3_fit1_par4_CB2_pol[2];

       for(int iParam=1; iParam < numOfPol1overPol2_Parm; iParam++)
       {
           paramCB2Pol2[iTail].paramsVariancePol1overPol2[iParam] = 0;
           paramCB2Pol2[iTail].paramsLimitVariancePol1overPol2[iParam][0] = 0;
           paramCB2Pol2[iTail].paramsLimitVariancePol1overPol2[iParam][1] = 0;
       }
       //paramCB2Pol2[iTail].paramsLimitPol2overPol3[5][0] = -100;
       //paramCB2Pol2[iTail].paramsLimitPol2overPol3[5][1] =  100;
       //paramCB2Pol2[iTail].paramsLimitPol2overPol3[6][0] = -100;
       //paramCB2Pol2[iTail].paramsLimitPol2overPol3[6][1] =  100;
       paramCB2Pol2[iTail].mean_jpsi = geant3_fit1_mean_jpsi_CB2_pol[0];
       paramCB2Pol2[iTail].mean_LowerLimits_jpsi = geant3_fit1_mean_jpsi_CB2_pol[1];
       paramCB2Pol2[iTail].mean_UpperLimits_jpsi = geant3_fit1_mean_jpsi_CB2_pol[2];
 
       paramCB2Pol2[iTail].width_jpsi = geant3_fit1_width_jpsi_CB2_pol[0];
       paramCB2Pol2[iTail].width_LowerLimits_jpsi = geant3_fit1_width_jpsi_CB2_pol[1];
       paramCB2Pol2[iTail].width_UpperLimits_jpsi = geant3_fit1_width_jpsi_CB2_pol[2];

       paramCB2Pol2[iTail].Norm_psi2S = geant3_fit1_norm_psi2S_CB2_pol[0];
       paramCB2Pol2[iTail].LowerLimits_psi2S = geant3_fit1_norm_psi2S_CB2_pol[1];
       paramCB2Pol2[iTail].UpperLimits_psi2S = geant3_fit1_norm_psi2S_CB2_pol[2];
     }
*/

      const int numOfPol2Cb2_Param = numOfPol2overPol3_Parm + numOfCB2_Parm;
    
//    for (int iCent=0; iCent < NUMOFCENT_BINS; iCent++)
      for (int iCent=0; iCent < 1; iCent++)
      {


        TF1 *fitFuncCB2Pol2 =  new TF1(paramCB2Pol2[iTail].strNameOfFitFunc.Data(),
                                       fitFunctionDoubleCB2POL2overPOL3,
                                       FitRange[iRange][0],
                                       FitRange[iRange][1],
                                       numOfPol2Cb2_Param+1);

        TF1 *fitPol2 = new TF1("fitPol2", 
                              BackgroundPOL2overPOL3,
                              FitRange[iRange][0],
                              FitRange[iRange][1],
                              numOfPol2overPol3_Parm);

        TF1 *fitPol2reject = new TF1("fitPol2reject", 
                              BackgroundPOL2overPOL3reject_JPsi_Psi2S,
                              FitRange[iRange][0],
                              FitRange[iRange][1],
                              numOfPol2overPol3_Parm);                          

        TF1 *fitCB2_JPsi = new TF1("fitCB2_JPsi", CrystalBallExtended, FitRange[iRange][0], FitRange[iRange][1], numOfCB2_Parm);

        TF1 *fitCB2_Psi2S = new TF1("fitCB2_Psi2S", CrystalBallExtended, FitRange[iRange][0], FitRange[iRange][1], numOfCB2_Parm);
     

/*
        outcome_fit = fitJPsi_Psi2S(iCent,
                FitRange[iRange][0],
                FitRange[iRange][1], 
//                histoRawInvM,
                histoInvMass[pt_bin],
//                histoInvMass_Pt12to15,
                fitPol2reject,
                fitFuncCB2Pol2,
                fitPol2,
                fitCB2_JPsi,
                fitCB2_Psi2S,
                paramCB2Pol2[iTail],
                Njpsi_CB2Pol2[iRange][iTail],
                SigmaJpsi_CB2Pol2[iRange][iTail],
                MassJpsi_CB2Pol2[iRange][iTail],
                ChiJpsi_CB2Pol2[iRange][iTail],
                ErrorOfNjpsi_CB2Pol2[iRange][iTail],
                ErrorOfSigmaJpsi_CB2Pol2[iRange][iTail],
                ErrorOfMassJpsi_CB2Pol2[iRange][iTail],
                ErrorOfChiJpsi_CB2Pol2[iRange][iTail],
                fitCanvas_CB2background[iCount1],
//                fitCanvas_background[iCount1],
                outcome_fit);

          fit_valide_CB2Pol[iRange][iTail][0] = outcome_fit.fit_valide;
          fit_status_CB2Pol[iRange][iTail][0] = outcome_fit.fit_status;
          fit_cov_matrix_status_CB2Pol[iRange][iTail][0] = outcome_fit.fit_cov_matrix_status;
*/
//        std::cout << "iCount1: " << iCount1 << std::endl;
        iCount1 +=1;

/*
        histoRawInvM->Draw();

        fitCB2_JPsi->Draw("SAME");
        fitCB2_Psi2S->Draw("SAME");
        fitPol2->Draw("SAME");
        fitFuncCB2Pol2->Draw("SAME");
*/



//        delete fitFuncCB2Pol2;
//        delete fitPol2;
//        delete fitCB2;
      } // end of for (int iCent=0; iCent < 1; iCent++)  
    } // end of for (int iTail=0; iTail< NUMOFTAILPARAM_CB2; iTail++)
  }  // end of for (iRange=0; iRange < 2; iRange++ )
 


  struct paramOfFunc paramNA60VWG[NUMOFTAILPARAM_NA60];  
  paramNA60VWG[0].strNameOfFitFunc = "fitFuncNA60VWG_geant3";
//  paramNA60VWG[1].strNameOfFitFunc = "fitFuncNA60VWG_geant3";

  struct paramOfFunc paramNA60Pol2[NUMOFTAILPARAM_NA60];
  paramNA60Pol2[0].strNameOfFitFunc = "fitFuncNA60Pol2_geant3";
//  paramNA60Pol2[1].strNameOfFitFunc = "fitFuncNA60Pol2_geant3";


  for (int iTail=0; iTail < NUMOFTAILPARAM_NA60; iTail++)
  {
    //NA60 + VWG
    paramNA60VWG[iTail].binWidth = binWidth;
    paramNA60VWG[iTail].NumOfBackGroundFuncParams = numOfVWG_Parm;
    paramNA60VWG[iTail].NumOfJPsiSignalFuncParams = numOfNA60_Parm+1;

    paramNA60VWG[iTail].NameOfBackGroundFuncParams = new TString[numOfVWG_Parm];
    for (int iName=0; iName<numOfVWG_Parm; iName++) paramNA60VWG[iTail].NameOfBackGroundFuncParams[iName].Append(NameOfVWGFuncParams[iName].Data());

    paramNA60VWG[iTail].NameOfJPsiSignalFuncParams = new TString[numOfNA60_Parm+1];
    for (int iName=0; iName<numOfNA60_Parm+1; iName++) paramNA60VWG[iTail].NameOfJPsiSignalFuncParams[iName].Append(NameOfDoubleNA60FuncParams[iName].Data());

    paramNA60VWG[iTail].ParamJPsiFunc = new double[numOfNA60_Parm+1];


    //NA60 + (Pol2 / Pol3)
    paramNA60Pol2[iTail].binWidth = binWidth;
    paramNA60Pol2[iTail].NumOfBackGroundFuncParams = numOfPol2overPol3_Parm;
    paramNA60Pol2[iTail].NumOfJPsiSignalFuncParams = numOfNA60_Parm+1;

    paramNA60Pol2[iTail].NameOfBackGroundFuncParams = new TString[numOfPol2overPol3_Parm];
    for (int iName=0; iName<numOfPol2overPol3_Parm; iName++) paramNA60Pol2[iTail].NameOfBackGroundFuncParams[iName].Append(NameOfPol2overPol3FuncParams[iName].Data());

    paramNA60Pol2[iTail].NameOfJPsiSignalFuncParams = new TString[numOfNA60_Parm+1];
    for (int iName=0; iName<numOfNA60_Parm+1; iName++) paramNA60Pol2[iTail].NameOfJPsiSignalFuncParams[iName].Append(NameOfDoubleNA60FuncParams[iName].Data());

    paramNA60Pol2[iTail].ParamJPsiFunc = new double[numOfNA60_Parm+1];
  }



  for(int iParam=0; iParam<numOfNA60_Parm; iParam++)
  {
    paramNA60VWG[0].ParamJPsiFunc[iParam] = Param_doubleNA60_geant3[iParam];
//    paramNA60VWG[1].ParamJPsiFunc[iParam] = Param_doubleNA60_geant3[iParam];
    paramNA60Pol2[0].ParamJPsiFunc[iParam] = Param_doubleNA60_geant3[iParam];
//    paramNA60Pol2[1].ParamJPsiFunc[iParam] = Param_doubleNA60_geant3[iParam];
  }


  double Njpsi_NA60VWG[FITRANGE][NUMOFTAILPARAM_NA60][NUMOFCENT_BINS];  
  double SigmaJpsi_NA60VWG[FITRANGE][NUMOFTAILPARAM_NA60][NUMOFCENT_BINS];
  double MassJpsi_NA60VWG[FITRANGE][NUMOFTAILPARAM_NA60][NUMOFCENT_BINS];
  double ChiJpsi_NA60VWG[FITRANGE][NUMOFTAILPARAM_NA60][NUMOFCENT_BINS];
  double ErrorOfNjpsi_NA60VWG[FITRANGE][NUMOFTAILPARAM_NA60][NUMOFCENT_BINS];
  double ErrorOfSigmaJpsi_NA60VWG[FITRANGE][NUMOFTAILPARAM_NA60][NUMOFCENT_BINS];
  double ErrorOfMassJpsi_NA60VWG[FITRANGE][NUMOFTAILPARAM_NA60][NUMOFCENT_BINS];
  double ErrorOfChiJpsi_NA60VWG[FITRANGE][NUMOFTAILPARAM_NA60][NUMOFCENT_BINS];

  double Njpsi_NA60Pol2[FITRANGE][NUMOFTAILPARAM_NA60][NUMOFCENT_BINS]; 
  double SigmaJpsi_NA60Pol2[FITRANGE][NUMOFTAILPARAM_NA60][NUMOFCENT_BINS];
  double MassJpsi_NA60Pol2[FITRANGE][NUMOFTAILPARAM_NA60][NUMOFCENT_BINS];
  double ChiJpsi_NA60Pol2[FITRANGE][NUMOFTAILPARAM_NA60][NUMOFCENT_BINS];
  double ErrorOfNjpsi_NA60Pol2[FITRANGE][NUMOFTAILPARAM_NA60][NUMOFCENT_BINS];
  double ErrorOfSigmaJpsi_NA60Pol2[FITRANGE][NUMOFTAILPARAM_NA60][NUMOFCENT_BINS];
  double ErrorOfMassJpsi_NA60Pol2[FITRANGE][NUMOFTAILPARAM_NA60][NUMOFCENT_BINS];
  double ErrorOfChiJpsi_NA60Pol2[FITRANGE][NUMOFTAILPARAM_NA60][NUMOFCENT_BINS];

  const int FIT_NUMBER_NA60 = FITRANGE * 2*  NUMOFTAILPARAM_NA60 * 1 ;

//  TCanvas *fitCanvas_background02[FIT_NUMBER_NA60];
  TCanvas *fitCanvas_NA60background[FIT_NUMBER_NA60];

  for(int ifit=0; ifit<FIT_NUMBER_NA60; ifit++)
  {
   fitCanvas_NA60background[ifit] = new TCanvas(Form("c3_%i",ifit),"",200, 10, 800, 600);
//   fitCanvas_background02[ifit] = new TCanvas(Form("c4_%i",ifit),"",200, 10, 800, 600);
  }  

  iCount1=0;
  for (int iRange=0; iRange < FITRANGE; iRange++ )
  {
    for (int iTail=0; iTail < NUMOFTAILPARAM_NA60; iTail++)
    {
      const int numOfVWGNA60_Param = numOfVWG_Parm + numOfNA60_Parm;
//      for (int iCent=0; iCent < NUMOFCENT_BINS; iCent++)
      
    if(iTail==0 && iRange == 0) // Geant4 2.0 - 4.8
    {
      for(int iParam=0; iParam<numOfNA60_Parm; iParam++) paramNA60VWG[iTail].ParamJPsiFunc[iParam] = values_geant3_NA60_tail_range00[pt_bin][iParam];

      paramNA60VWG[iTail].paramsInitialVWG[0]= geant3_fit0_par0_NA60_VWG[0];
      paramNA60VWG[iTail].paramsLimitVWG[0][0] = geant3_fit0_par0_NA60_VWG[1];
      paramNA60VWG[iTail].paramsLimitVWG[0][1] = geant3_fit0_par0_NA60_VWG[2];

      paramNA60VWG[iTail].paramsInitialVWG[1] = geant3_fit0_par1_NA60_VWG[0];
      paramNA60VWG[iTail].paramsLimitVWG[1][0] = geant3_fit0_par1_NA60_VWG[1];
      paramNA60VWG[iTail].paramsLimitVWG[1][1] = geant3_fit0_par1_NA60_VWG[2];

      paramNA60VWG[iTail].paramsInitialVWG[2] = geant3_fit0_par2_NA60_VWG[0];
      paramNA60VWG[iTail].paramsLimitVWG[2][0] = geant3_fit0_par2_NA60_VWG[1];
      paramNA60VWG[iTail].paramsLimitVWG[2][1] = geant3_fit0_par2_NA60_VWG[2];

      paramNA60VWG[iTail].paramsInitialVWG[3] = geant3_fit0_par3_NA60_VWG[0];
      paramNA60VWG[iTail].paramsLimitVWG[3][0] = geant3_fit0_par3_NA60_VWG[1];
      paramNA60VWG[iTail].paramsLimitVWG[3][1] = geant3_fit0_par3_NA60_VWG[2];

      paramNA60VWG[iTail].paramsInitialVWG[4] = geant3_fit0_par4_NA60_VWG[0];
      paramNA60VWG[iTail].paramsLimitVWG[4][0] = geant3_fit0_par4_NA60_VWG[1];
      paramNA60VWG[iTail].paramsLimitVWG[4][1] = geant3_fit0_par4_NA60_VWG[2];

      paramNA60VWG[iTail].paramsVarianceVWG[0] = geant3_fit0_par0_variance_NA60_VWG[0];
      paramNA60VWG[iTail].paramsLimitVarianceVWG[0][0] = geant3_fit0_par0_variance_NA60_VWG[1];
      paramNA60VWG[iTail].paramsLimitVarianceVWG[0][1] = geant3_fit0_par0_variance_NA60_VWG[2];

      paramNA60VWG[iTail].paramsVarianceVWG[1] = geant3_fit0_par1_variance_NA60_VWG[0];
      paramNA60VWG[iTail].paramsLimitVarianceVWG[1][0]=geant3_fit0_par1_variance_NA60_VWG[1];
      paramNA60VWG[iTail].paramsLimitVarianceVWG[1][1]=geant3_fit0_par1_variance_NA60_VWG[2];

      paramNA60VWG[iTail].paramsVarianceVWG[2] = geant3_fit0_par2_variance_NA60_VWG[0];
      paramNA60VWG[iTail].paramsLimitVarianceVWG[2][0]=geant3_fit0_par2_variance_NA60_VWG[1];
      paramNA60VWG[iTail].paramsLimitVarianceVWG[2][1]=geant3_fit0_par2_variance_NA60_VWG[2];

      paramNA60VWG[iTail].paramsVarianceVWG[3] = geant3_fit0_par3_variance_NA60_VWG[0];
      paramNA60VWG[iTail].paramsLimitVarianceVWG[3][0]=geant3_fit0_par3_variance_NA60_VWG[1];
      paramNA60VWG[iTail].paramsLimitVarianceVWG[3][1]=geant3_fit0_par3_variance_NA60_VWG[2];

      paramNA60VWG[iTail].paramsVarianceVWG[4] = geant3_fit0_par4_variance_NA60_VWG[0];
      paramNA60VWG[iTail].paramsLimitVarianceVWG[4][0]=geant3_fit0_par4_variance_NA60_VWG[1];
      paramNA60VWG[iTail].paramsLimitVarianceVWG[4][1]=geant3_fit0_par4_variance_NA60_VWG[2];


      paramNA60VWG[iTail].mean_jpsi = geant3_fit0_mean_jpsi_NA60_VWG[0];
      paramNA60VWG[iTail].mean_LowerLimits_jpsi = geant3_fit0_mean_jpsi_NA60_VWG[1];
      paramNA60VWG[iTail].mean_UpperLimits_jpsi = geant3_fit0_mean_jpsi_NA60_VWG[2];
 
      paramNA60VWG[iTail].width_jpsi = geant3_fit0_width_jpsi_NA60_VWG[0];
      paramNA60VWG[iTail].width_LowerLimits_jpsi = geant3_fit0_width_jpsi_NA60_VWG[1];
      paramNA60VWG[iTail].width_UpperLimits_jpsi = geant3_fit0_width_jpsi_NA60_VWG[2];

      paramNA60VWG[iTail].Norm_psi2S = geant3_fit0_norm_psi2S_NA60_VWG[0];
      paramNA60VWG[iTail].LowerLimits_psi2S = geant3_fit0_norm_psi2S_NA60_VWG[1];
      paramNA60VWG[iTail].UpperLimits_psi2S = geant3_fit0_norm_psi2S_NA60_VWG[2];
 
    }
    if(iTail==0 && iRange == 1)
    {

      for(int iParam=0; iParam<numOfNA60_Parm; iParam++) paramNA60VWG[iTail].ParamJPsiFunc[iParam] = values_geant3_NA60_tail_range01[pt_bin][iParam];

      paramNA60VWG[iTail].paramsInitialVWG[0]= geant3_fit1_par0_NA60_VWG[0];
      paramNA60VWG[iTail].paramsLimitVWG[0][0] = geant3_fit1_par0_NA60_VWG[1];
      paramNA60VWG[iTail].paramsLimitVWG[0][1] = geant3_fit1_par0_NA60_VWG[2];

      paramNA60VWG[iTail].paramsInitialVWG[1] = geant3_fit1_par1_NA60_VWG[0];
      paramNA60VWG[iTail].paramsLimitVWG[1][0] = geant3_fit1_par1_NA60_VWG[1];
      paramNA60VWG[iTail].paramsLimitVWG[1][1] = geant3_fit1_par1_NA60_VWG[2];

      paramNA60VWG[iTail].paramsInitialVWG[2] = geant3_fit1_par2_NA60_VWG[0];
      paramNA60VWG[iTail].paramsLimitVWG[2][0] = geant3_fit1_par2_NA60_VWG[1];
      paramNA60VWG[iTail].paramsLimitVWG[2][1] = geant3_fit1_par2_NA60_VWG[2];

      paramNA60VWG[iTail].paramsInitialVWG[3] = geant3_fit1_par3_NA60_VWG[0];
      paramNA60VWG[iTail].paramsLimitVWG[3][0] = geant3_fit1_par3_NA60_VWG[1];
      paramNA60VWG[iTail].paramsLimitVWG[3][1] = geant3_fit1_par3_NA60_VWG[2];

      paramNA60VWG[iTail].paramsInitialVWG[4] = geant3_fit1_par4_NA60_VWG[0];
      paramNA60VWG[iTail].paramsLimitVWG[4][0] = geant3_fit1_par4_NA60_VWG[1];
      paramNA60VWG[iTail].paramsLimitVWG[4][1] = geant3_fit1_par4_NA60_VWG[2];

      paramNA60VWG[iTail].paramsVarianceVWG[0] = geant3_fit1_par0_variance_NA60_VWG[0];
      paramNA60VWG[iTail].paramsLimitVarianceVWG[0][0] = geant3_fit1_par0_variance_NA60_VWG[1];
      paramNA60VWG[iTail].paramsLimitVarianceVWG[0][1] = geant3_fit1_par0_variance_NA60_VWG[2];

      paramNA60VWG[iTail].paramsVarianceVWG[1] = geant3_fit1_par1_variance_NA60_VWG[0];
      paramNA60VWG[iTail].paramsLimitVarianceVWG[1][0]=geant3_fit1_par1_variance_NA60_VWG[1];
      paramNA60VWG[iTail].paramsLimitVarianceVWG[1][1]=geant3_fit1_par1_variance_NA60_VWG[2];

      paramNA60VWG[iTail].paramsVarianceVWG[2] = geant3_fit1_par2_variance_NA60_VWG[0];
      paramNA60VWG[iTail].paramsLimitVarianceVWG[2][0]=geant3_fit1_par2_variance_NA60_VWG[1];
      paramNA60VWG[iTail].paramsLimitVarianceVWG[2][1]=geant3_fit1_par2_variance_NA60_VWG[2];

      paramNA60VWG[iTail].paramsVarianceVWG[3] = geant3_fit1_par3_variance_NA60_VWG[0];
      paramNA60VWG[iTail].paramsLimitVarianceVWG[3][0] = geant3_fit1_par3_variance_NA60_VWG[1];
      paramNA60VWG[iTail].paramsLimitVarianceVWG[3][1] = geant3_fit1_par3_variance_NA60_VWG[2];

      paramNA60VWG[iTail].paramsVarianceVWG[4] = geant3_fit1_par4_variance_NA60_VWG[0];
      paramNA60VWG[iTail].paramsLimitVarianceVWG[4][0] = geant3_fit1_par4_variance_NA60_VWG[1];
      paramNA60VWG[iTail].paramsLimitVarianceVWG[4][1] = geant3_fit1_par4_variance_NA60_VWG[2];


      paramNA60VWG[iTail].mean_jpsi = geant3_fit1_mean_jpsi_NA60_VWG[0];
      paramNA60VWG[iTail].mean_LowerLimits_jpsi = geant3_fit1_mean_jpsi_NA60_VWG[1];
      paramNA60VWG[iTail].mean_UpperLimits_jpsi = geant3_fit1_mean_jpsi_NA60_VWG[2];
 
      paramNA60VWG[iTail].width_jpsi = geant3_fit1_width_jpsi_NA60_VWG[0];
      paramNA60VWG[iTail].width_LowerLimits_jpsi = geant3_fit1_width_jpsi_NA60_VWG[1];
      paramNA60VWG[iTail].width_UpperLimits_jpsi = geant3_fit1_width_jpsi_NA60_VWG[2];

      paramNA60VWG[iTail].Norm_psi2S = geant3_fit1_norm_psi2S_NA60_VWG[0];
      paramNA60VWG[iTail].LowerLimits_psi2S = geant3_fit1_norm_psi2S_NA60_VWG[1];
      paramNA60VWG[iTail].UpperLimits_psi2S = geant3_fit1_norm_psi2S_NA60_VWG[2];
 
    }
/*
     if(iTail==1 && iRange==0)  //NA60 VWG geant3 2.0 4.8
    {
     
      for(int iParam=0; iParam<numOfNA60_Parm; iParam++) paramNA60VWG[iTail].ParamJPsiFunc[iParam] = values_geant3_NA60_tail_range00[pt_bin][iParam];
      paramNA60VWG[iTail].paramsInitialVWG1[0]= geant3_fit0_par0_NA60_VWG1[0];
      paramNA60VWG[iTail].paramsLimitVWG1[0][0] = geant3_fit0_par0_NA60_VWG1[1];
      paramNA60VWG[iTail].paramsLimitVWG1[0][1] = geant3_fit0_par0_NA60_VWG1[2];

      paramNA60VWG[iTail].paramsInitialVWG1[1] = geant3_fit0_par1_NA60_VWG1[0];
      paramNA60VWG[iTail].paramsLimitVWG1[1][0] = geant3_fit0_par1_NA60_VWG1[1];
      paramNA60VWG[iTail].paramsLimitVWG1[1][1] = geant3_fit0_par1_NA60_VWG1[2];

      paramNA60VWG[iTail].paramsInitialVWG1[2] = geant3_fit0_par2_NA60_VWG1[0];
      paramNA60VWG[iTail].paramsLimitVWG1[2][0] = geant3_fit0_par2_NA60_VWG1[1];
      paramNA60VWG[iTail].paramsLimitVWG1[2][1] = geant3_fit0_par2_NA60_VWG1[2];

      paramNA60VWG[iTail].paramsInitialVWG1[3] = geant3_fit0_par3_NA60_VWG1[0];
      paramNA60VWG[iTail].paramsLimitVWG1[3][0] = geant3_fit0_par3_NA60_VWG1[1];
      paramNA60VWG[iTail].paramsLimitVWG1[3][1] = geant3_fit0_par3_NA60_VWG1[2];


      paramNA60VWG[iTail].paramsVarianceVWG1[0] = geant3_fit0_par0_variance_NA60_VWG1[0];
      paramNA60VWG[iTail].paramsLimitVarianceVWG1[0][0] = geant3_fit0_par0_variance_NA60_VWG1[1];
      paramNA60VWG[iTail].paramsLimitVarianceVWG1[0][1] = geant3_fit0_par0_variance_NA60_VWG1[2];

      paramNA60VWG[iTail].paramsVarianceVWG1[1] = geant3_fit0_par1_variance_NA60_VWG1[0];
      paramNA60VWG[iTail].paramsLimitVarianceVWG1[1][0]=geant3_fit0_par1_variance_NA60_VWG1[1];
      paramNA60VWG[iTail].paramsLimitVarianceVWG1[1][1]=geant3_fit0_par1_variance_NA60_VWG1[2];

      paramNA60VWG[iTail].paramsVarianceVWG1[2] = geant3_fit0_par2_variance_NA60_VWG1[0];
      paramNA60VWG[iTail].paramsLimitVarianceVWG1[2][0]=geant3_fit0_par2_variance_NA60_VWG1[1];
      paramNA60VWG[iTail].paramsLimitVarianceVWG1[2][1]=geant3_fit0_par2_variance_NA60_VWG1[2];

      paramNA60VWG[iTail].paramsVarianceVWG1[3] = geant3_fit0_par3_variance_NA60_VWG1[0];
      paramNA60VWG[iTail].paramsLimitVarianceVWG1[3][0]=geant3_fit0_par3_variance_NA60_VWG1[1];
      paramNA60VWG[iTail].paramsLimitVarianceVWG1[3][1]=geant3_fit0_par3_variance_NA60_VWG1[2];
      //paramNA60VWG[iTail].paramsLimitVWG[4][0] = 0;
      //paramNA60VWG[iTail].paramsLimitVWG[4][1] = 100;

       paramNA60VWG[iTail].mean_jpsi = geant3_fit0_mean_jpsi_NA60_VWG1[0];
      paramNA60VWG[iTail].mean_LowerLimits_jpsi = geant3_fit0_mean_jpsi_NA60_VWG1[1];
      paramNA60VWG[iTail].mean_UpperLimits_jpsi = geant3_fit0_mean_jpsi_NA60_VWG1[2];
 
      paramNA60VWG[iTail].width_jpsi = geant3_fit0_width_jpsi_NA60_VWG1[0];
      paramNA60VWG[iTail].width_LowerLimits_jpsi = geant3_fit0_width_jpsi_NA60_VWG1[1];
      paramNA60VWG[iTail].width_UpperLimits_jpsi = geant3_fit0_width_jpsi_NA60_VWG1[2];

      paramNA60VWG[iTail].Norm_psi2S = geant3_fit0_norm_psi2S_NA60_VWG1[0];
      paramNA60VWG[iTail].LowerLimits_psi2S = geant3_fit0_norm_psi2S_NA60_VWG1[1];
      paramNA60VWG[iTail].UpperLimits_psi2S = geant3_fit0_norm_psi2S_NA60_VWG1[2];

     }

     if(iTail==1 && iRange==1)       //NA60 VWG geant3 2.2 4.4
     {
      for(int iParam=0; iParam<numOfNA60_Parm; iParam++) paramNA60VWG[iTail].ParamJPsiFunc[iParam] = values_geant3_NA60_tail_range01[pt_bin][iParam];

      paramNA60VWG[iTail].paramsInitialVWG1[0]= geant3_fit1_par0_NA60_VWG1[0];
      paramNA60VWG[iTail].paramsLimitVWG1[0][0] = geant3_fit1_par0_NA60_VWG1[1];
      paramNA60VWG[iTail].paramsLimitVWG1[0][1] = geant3_fit1_par0_NA60_VWG1[2];

      paramNA60VWG[iTail].paramsInitialVWG1[1] = geant3_fit1_par1_NA60_VWG1[0];
      paramNA60VWG[iTail].paramsLimitVWG1[1][0] = geant3_fit1_par1_NA60_VWG1[1];
      paramNA60VWG[iTail].paramsLimitVWG1[1][1] = geant3_fit1_par1_NA60_VWG1[2];

      paramNA60VWG[iTail].paramsInitialVWG1[2] = geant3_fit1_par2_NA60_VWG1[0];
      paramNA60VWG[iTail].paramsLimitVWG1[2][0] = geant3_fit1_par2_NA60_VWG1[1];
      paramNA60VWG[iTail].paramsLimitVWG1[2][1] = geant3_fit1_par2_NA60_VWG1[2];

      paramNA60VWG[iTail].paramsInitialVWG1[3] = geant3_fit1_par3_NA60_VWG1[0];
      paramNA60VWG[iTail].paramsLimitVWG1[3][0] = geant3_fit1_par3_NA60_VWG1[1];
      paramNA60VWG[iTail].paramsLimitVWG1[3][1] = geant3_fit1_par3_NA60_VWG1[2];

      paramNA60VWG[iTail].paramsVarianceVWG1[0] = geant3_fit1_par0_variance_NA60_VWG1[0];
      paramNA60VWG[iTail].paramsLimitVarianceVWG1[0][0] = geant3_fit1_par0_variance_NA60_VWG1[1];
      paramNA60VWG[iTail].paramsLimitVarianceVWG1[0][1] = geant3_fit1_par0_variance_NA60_VWG1[2];

      paramNA60VWG[iTail].paramsVarianceVWG1[1] = geant3_fit1_par1_variance_NA60_VWG1[0];
      paramNA60VWG[iTail].paramsLimitVarianceVWG1[1][0]=geant3_fit1_par1_variance_NA60_VWG1[1];
      paramNA60VWG[iTail].paramsLimitVarianceVWG1[1][1]=geant3_fit1_par1_variance_NA60_VWG1[2];

      paramNA60VWG[iTail].paramsVarianceVWG1[2] = geant3_fit1_par2_variance_NA60_VWG1[0];
      paramNA60VWG[iTail].paramsLimitVarianceVWG1[2][0]=geant3_fit1_par2_variance_NA60_VWG1[1];
      paramNA60VWG[iTail].paramsLimitVarianceVWG1[2][1]=geant3_fit1_par2_variance_NA60_VWG1[2];

      paramNA60VWG[iTail].paramsVarianceVWG1[3] = geant3_fit1_par3_variance_NA60_VWG1[0];
      paramNA60VWG[iTail].paramsLimitVarianceVWG1[3][0]=geant3_fit1_par3_variance_NA60_VWG1[1];
      paramNA60VWG[iTail].paramsLimitVarianceVWG1[3][1]=geant3_fit1_par3_variance_NA60_VWG1[2];
      //paramNA60VWG[iTail].paramsLimitVWG[4][0] = 0;
      //paramNA60VWG[iTail].paramsLimitVWG[4][1] = 100;

       paramNA60VWG[iTail].mean_jpsi = geant3_fit1_mean_jpsi_NA60_VWG1[0];
      paramNA60VWG[iTail].mean_LowerLimits_jpsi = geant3_fit1_mean_jpsi_NA60_VWG1[1];
      paramNA60VWG[iTail].mean_UpperLimits_jpsi = geant3_fit1_mean_jpsi_NA60_VWG1[2];
 
      paramNA60VWG[iTail].width_jpsi = geant3_fit1_width_jpsi_NA60_VWG1[0];
      paramNA60VWG[iTail].width_LowerLimits_jpsi = geant3_fit1_width_jpsi_NA60_VWG1[1];
      paramNA60VWG[iTail].width_UpperLimits_jpsi = geant3_fit1_width_jpsi_NA60_VWG1[2];

      paramNA60VWG[iTail].Norm_psi2S = geant3_fit1_norm_psi2S_NA60_VWG1[0];
      paramNA60VWG[iTail].LowerLimits_psi2S = geant3_fit1_norm_psi2S_NA60_VWG1[1];
      paramNA60VWG[iTail].UpperLimits_psi2S = geant3_fit1_norm_psi2S_NA60_VWG1[2];
     }
*/

      for (int iCent=0; iCent < 1; iCent++)
      {


        TF1 *fitFuncNA60VWG =  new TF1(paramNA60VWG[iTail].strNameOfFitFunc.Data(), 
                                       fitFunctionDoubleNA60VWG, 
                                       FitRange[iRange][0], 
                                       FitRange[iRange][1],
                                       numOfVWGNA60_Param+1); 

        TF1 *fitVWGreject = new TF1("fitVWGreject", 
                            BackgroundVWGreject_JPsi_Psi2S,
                            FitRange[iRange][0],
                            FitRange[iRange][1],
                            numOfVWG_Parm);

        TF1 *fitVWG = new TF1("fitVWG", 
                              BackgroundVWG, 
                              FitRange[iRange][0], 
                              FitRange[iRange][1], 
                              numOfVWG_Parm);

        TF1 *fitNA60_JPsi = new TF1("fitNA60_JPsi", NA60, FitRange[iRange][0], FitRange[iRange][1],  numOfNA60_Parm);

        TF1 *fitNA60_Psi2S = new TF1("fitNA60_Psi2S", NA60, FitRange[iRange][0], FitRange[iRange][1],  numOfNA60_Parm);

/*
        outcome_fit = fitJPsi_Psi2S(iCent,
                FitRange[iRange][0],
                FitRange[iRange][1], 
//                histoRawInvM,
                histoInvMass[pt_bin],
//                histoInvMass_Pt12to15,
                fitVWGreject,
                fitFuncNA60VWG,
                fitVWG,
                fitNA60_JPsi,
                fitNA60_Psi2S,
                paramNA60VWG[iTail],
                Njpsi_NA60VWG[iRange][iTail],
                SigmaJpsi_NA60VWG[iRange][iTail],
                MassJpsi_NA60VWG[iRange][iTail],
                ChiJpsi_NA60VWG[iRange][iTail],
                ErrorOfNjpsi_NA60VWG[iRange][iTail],
                ErrorOfSigmaJpsi_NA60VWG[iRange][iTail],
                ErrorOfMassJpsi_NA60VWG[iRange][iTail],
                ErrorOfChiJpsi_NA60VWG[iRange][iTail],
                fitCanvas_NA60background[iCount1],
//                fitCanvas_background02[iCount1],
                outcome_fit);

          fit_valide_NA60VWG[iRange][iTail][0] = outcome_fit.fit_valide;
          fit_status_NA60VWG[iRange][iTail][0] = outcome_fit.fit_status;
          fit_cov_matrix_status_NA60VWG[iRange][iTail][0] = outcome_fit.fit_cov_matrix_status;
*/
        iCount1+=1;

//        delete fitFuncNA60VWG;
//        delete fitVWG;
//        delete fitNA60;
     }
   }
  }


  for (int iRange=0; iRange < FITRANGE; iRange++ )
  {
    for (int iTail=0; iTail < NUMOFTAILPARAM_NA60; iTail++)
    {
     
     if(iTail==0 && iRange == 0)       
     {

      for(int iParam=0; iParam<numOfNA60_Parm; iParam++) paramNA60Pol2[iTail].ParamJPsiFunc[iParam] = values_geant3_NA60_tail_range00[pt_bin][iParam];

       paramNA60Pol2[iTail].paramsInitialPol2overPol3[0] = geant3_fit0_par0_NA60_pol[0];
       paramNA60Pol2[iTail].paramsLimitPol2overPol3[0][0] = geant3_fit0_par0_NA60_pol[1];
       paramNA60Pol2[iTail].paramsLimitPol2overPol3[0][1] = geant3_fit0_par0_NA60_pol[2];

       paramNA60Pol2[iTail].paramsInitialPol2overPol3[1] = geant3_fit0_par1_NA60_pol[0];
       paramNA60Pol2[iTail].paramsLimitPol2overPol3[1][0] = geant3_fit0_par1_NA60_pol[1];
       paramNA60Pol2[iTail].paramsLimitPol2overPol3[1][1] = geant3_fit0_par1_NA60_pol[2];

       paramNA60Pol2[iTail].paramsInitialPol2overPol3[2] = geant3_fit0_par2_NA60_pol[0];
       paramNA60Pol2[iTail].paramsLimitPol2overPol3[2][0] = geant3_fit0_par2_NA60_pol[1];
       paramNA60Pol2[iTail].paramsLimitPol2overPol3[2][1] = geant3_fit0_par2_NA60_pol[2];

       paramNA60Pol2[iTail].paramsInitialPol2overPol3[3] = geant3_fit0_par3_NA60_pol[0];
       paramNA60Pol2[iTail].paramsLimitPol2overPol3[3][0] = geant3_fit0_par3_NA60_pol[1];
       paramNA60Pol2[iTail].paramsLimitPol2overPol3[3][1] = geant3_fit0_par3_NA60_pol[2];

       paramNA60Pol2[iTail].paramsInitialPol2overPol3[4] = geant3_fit0_par4_NA60_pol[0];
       paramNA60Pol2[iTail].paramsLimitPol2overPol3[4][0] = geant3_fit0_par4_NA60_pol[1];
       paramNA60Pol2[iTail].paramsLimitPol2overPol3[4][1] = geant3_fit0_par4_NA60_pol[2];

       paramNA60Pol2[iTail].paramsInitialPol2overPol3[5] = geant3_fit0_par5_NA60_pol[0];
       paramNA60Pol2[iTail].paramsLimitPol2overPol3[5][0] = geant3_fit0_par5_NA60_pol[1];
       paramNA60Pol2[iTail].paramsLimitPol2overPol3[5][1] = geant3_fit0_par5_NA60_pol[2];

       paramNA60Pol2[iTail].paramsInitialPol2overPol3[6] = geant3_fit0_par6_NA60_pol[0];
       paramNA60Pol2[iTail].paramsLimitPol2overPol3[6][0] = geant3_fit0_par6_NA60_pol[1];
       paramNA60Pol2[iTail].paramsLimitPol2overPol3[6][1] =  geant3_fit0_par6_NA60_pol[2];

       paramNA60Pol2[iTail].mean_jpsi = geant3_fit0_mean_jpsi_NA60_pol[0];
       paramNA60Pol2[iTail].mean_LowerLimits_jpsi = geant3_fit0_mean_jpsi_NA60_pol[1];
       paramNA60Pol2[iTail].mean_UpperLimits_jpsi = geant3_fit0_mean_jpsi_NA60_pol[2];
  
       paramNA60Pol2[iTail].width_jpsi = geant3_fit0_width_jpsi_NA60_pol[0];
       paramNA60Pol2[iTail].width_LowerLimits_jpsi = geant3_fit0_width_jpsi_NA60_pol[1];
       paramNA60Pol2[iTail].width_UpperLimits_jpsi = geant3_fit0_width_jpsi_NA60_pol[2];

       paramNA60Pol2[iTail].Norm_psi2S = geant3_fit0_norm_psi2S_NA60_pol[0];
       paramNA60Pol2[iTail].LowerLimits_psi2S = geant3_fit0_norm_psi2S_NA60_pol[1];
       paramNA60Pol2[iTail].UpperLimits_psi2S = geant3_fit0_norm_psi2S_NA60_pol[2];
     }

     if(iTail==0 && iRange == 1)        //NA60Pol2_geant4_2.40_4.70
     {
      for(int iParam=0; iParam<numOfNA60_Parm; iParam++) paramNA60Pol2[iTail].ParamJPsiFunc[iParam] = values_geant3_NA60_tail_range01[pt_bin][iParam];

       paramNA60Pol2[iTail].paramsInitialPol2overPol3[0] = geant3_fit1_par0_NA60_pol[0];
       paramNA60Pol2[iTail].paramsLimitPol2overPol3[0][0] = geant3_fit1_par0_NA60_pol[1];
       paramNA60Pol2[iTail].paramsLimitPol2overPol3[0][1] = geant3_fit1_par0_NA60_pol[2];

       paramNA60Pol2[iTail].paramsInitialPol2overPol3[1] = geant3_fit1_par1_NA60_pol[0];
       paramNA60Pol2[iTail].paramsLimitPol2overPol3[1][0] = geant3_fit1_par1_NA60_pol[1];
       paramNA60Pol2[iTail].paramsLimitPol2overPol3[1][1] = geant3_fit1_par1_NA60_pol[2];

       paramNA60Pol2[iTail].paramsInitialPol2overPol3[2] = geant3_fit1_par2_NA60_pol[0];
       paramNA60Pol2[iTail].paramsLimitPol2overPol3[2][0] = geant3_fit1_par2_NA60_pol[1];
       paramNA60Pol2[iTail].paramsLimitPol2overPol3[2][1] = geant3_fit1_par2_NA60_pol[2];

       paramNA60Pol2[iTail].paramsInitialPol2overPol3[3] = geant3_fit1_par3_NA60_pol[0];
       paramNA60Pol2[iTail].paramsLimitPol2overPol3[3][0] = geant3_fit1_par3_NA60_pol[1];
       paramNA60Pol2[iTail].paramsLimitPol2overPol3[3][1] = geant3_fit1_par3_NA60_pol[2];

       paramNA60Pol2[iTail].paramsInitialPol2overPol3[4] = geant3_fit1_par4_NA60_pol[0];
       paramNA60Pol2[iTail].paramsLimitPol2overPol3[4][0] = geant3_fit1_par4_NA60_pol[1];
       paramNA60Pol2[iTail].paramsLimitPol2overPol3[4][1] = geant3_fit1_par4_NA60_pol[2];

       paramNA60Pol2[iTail].paramsInitialPol2overPol3[5] = geant3_fit1_par5_NA60_pol[0];
       paramNA60Pol2[iTail].paramsLimitPol2overPol3[5][0] = geant3_fit1_par5_NA60_pol[1];
       paramNA60Pol2[iTail].paramsLimitPol2overPol3[5][1] = geant3_fit1_par5_NA60_pol[2];

       paramNA60Pol2[iTail].paramsInitialPol2overPol3[6] = geant3_fit1_par6_NA60_pol[0];
       paramNA60Pol2[iTail].paramsLimitPol2overPol3[6][0] = geant3_fit1_par6_NA60_pol[1];
       paramNA60Pol2[iTail].paramsLimitPol2overPol3[6][1] = geant3_fit1_par6_NA60_pol[2];

       paramNA60Pol2[iTail].mean_jpsi = geant3_fit1_mean_jpsi_NA60_pol[0];
       paramNA60Pol2[iTail].mean_LowerLimits_jpsi = geant3_fit1_mean_jpsi_NA60_pol[1];
       paramNA60Pol2[iTail].mean_UpperLimits_jpsi = geant3_fit1_mean_jpsi_NA60_pol[2];
  
       paramNA60Pol2[iTail].width_jpsi = geant3_fit1_width_jpsi_NA60_pol[0];
       paramNA60Pol2[iTail].width_LowerLimits_jpsi = geant3_fit1_width_jpsi_NA60_pol[1];
       paramNA60Pol2[iTail].width_UpperLimits_jpsi = geant3_fit1_width_jpsi_NA60_pol[2];

       paramNA60Pol2[iTail].Norm_psi2S = geant3_fit1_norm_psi2S_NA60_pol[0];
       paramNA60Pol2[iTail].LowerLimits_psi2S = geant3_fit1_norm_psi2S_NA60_pol[1];
       paramNA60Pol2[iTail].UpperLimits_psi2S = geant3_fit1_norm_psi2S_NA60_pol[2];
     }
/*
     if(iTail==1 && iRange == 0)        //NA60Pol2_geant3_2.20_4.50
     {
       for(int iParam=0; iParam<numOfNA60_Parm; iParam++) paramNA60Pol2[iTail].ParamJPsiFunc[iParam] = values_geant3_NA60_tail_range00[pt_bin][iParam];
       paramNA60Pol2[iTail].paramsLimitPol1overPol2[0][0] = geant3_fit0_par0_NA60_pol[1];
       paramNA60Pol2[iTail].paramsLimitPol1overPol2[0][1] = geant3_fit0_par0_NA60_pol[2];

       paramNA60Pol2[iTail].paramsInitialPol1overPol2[1] = geant3_fit0_par1_NA60_pol[0];
       paramNA60Pol2[iTail].paramsLimitPol1overPol2[1][0] = geant3_fit0_par1_NA60_pol[1];
       paramNA60Pol2[iTail].paramsLimitPol1overPol2[1][1] = geant3_fit0_par1_NA60_pol[2];

       paramNA60Pol2[iTail].paramsInitialPol1overPol2[2] = geant3_fit0_par2_NA60_pol[0];
       paramNA60Pol2[iTail].paramsLimitPol1overPol2[2][0] = geant3_fit0_par2_NA60_pol[1];
       paramNA60Pol2[iTail].paramsLimitPol1overPol2[2][1] = geant3_fit0_par2_NA60_pol[2];

       paramNA60Pol2[iTail].paramsInitialPol1overPol2[3] = geant3_fit0_par3_NA60_pol[0];
       paramNA60Pol2[iTail].paramsLimitPol1overPol2[3][0] = geant3_fit0_par3_NA60_pol[1];
       paramNA60Pol2[iTail].paramsLimitPol1overPol2[3][1] = geant3_fit0_par3_NA60_pol[2];

       paramNA60Pol2[iTail].paramsInitialPol1overPol2[4] = geant3_fit0_par4_NA60_pol[0];
       paramNA60Pol2[iTail].paramsLimitPol1overPol2[4][0] = geant3_fit0_par4_NA60_pol[1];
       paramNA60Pol2[iTail].paramsLimitPol1overPol2[4][1] = geant3_fit0_par4_NA60_pol[2];
       //paramNA60Pol2[iTail].paramsLimitPol2overPol3[5][0] = -100;
       //paramNA60Pol2[iTail].paramsLimitPol2overPol3[5][1] =  100;
       //paramNA60Pol2[iTail].paramsLimitPol2overPol3[6][0] = -100;
       //paramNA60Pol2[iTail].paramsLimitPol2overPol3[6][1] =  100;
       paramNA60Pol2[iTail].mean_jpsi = geant3_fit0_mean_jpsi_NA60_pol[0];
       paramNA60Pol2[iTail].mean_LowerLimits_jpsi = geant3_fit0_mean_jpsi_NA60_pol[1];
       paramNA60Pol2[iTail].mean_UpperLimits_jpsi = geant3_fit0_mean_jpsi_NA60_pol[2];
  
       paramNA60Pol2[iTail].width_jpsi = geant3_fit0_width_jpsi_NA60_pol[0];
       paramNA60Pol2[iTail].width_LowerLimits_jpsi = geant3_fit0_width_jpsi_NA60_pol[1];
       paramNA60Pol2[iTail].width_UpperLimits_jpsi = geant3_fit0_width_jpsi_NA60_pol[2];

       paramNA60Pol2[iTail].Norm_psi2S = geant3_fit0_norm_psi2S_NA60_pol[0];
       paramNA60Pol2[iTail].LowerLimits_psi2S = geant3_fit0_norm_psi2S_NA60_pol[1];
       paramNA60Pol2[iTail].UpperLimits_psi2S = geant3_fit0_norm_psi2S_NA60_pol[2];
     }

     if(iTail==1 && iRange == 1)        //NA60Pol2_geant3_2.40_4.70
     {
       for(int iParam=0; iParam<numOfNA60_Parm; iParam++) paramNA60Pol2[iTail].ParamJPsiFunc[iParam] = values_geant3_NA60_tail_range01[pt_bin][iParam];

       paramNA60Pol2[iTail].paramsLimitPol1overPol2[0][0] = geant3_fit1_par0_NA60_pol[1];
       paramNA60Pol2[iTail].paramsLimitPol1overPol2[0][1] = geant3_fit1_par0_NA60_pol[2];

       paramNA60Pol2[iTail].paramsInitialPol1overPol2[1] = geant3_fit1_par1_NA60_pol[0];
       paramNA60Pol2[iTail].paramsLimitPol1overPol2[1][0] = geant3_fit1_par1_NA60_pol[1];
       paramNA60Pol2[iTail].paramsLimitPol1overPol2[1][1] = geant3_fit1_par1_NA60_pol[2];

       paramNA60Pol2[iTail].paramsInitialPol1overPol2[2] = geant3_fit1_par2_NA60_pol[0];
       paramNA60Pol2[iTail].paramsLimitPol1overPol2[2][0] = geant3_fit1_par2_NA60_pol[1];
       paramNA60Pol2[iTail].paramsLimitPol1overPol2[2][1] = geant3_fit1_par2_NA60_pol[2];

       paramNA60Pol2[iTail].paramsInitialPol1overPol2[3] = geant3_fit1_par3_NA60_pol[0];
       paramNA60Pol2[iTail].paramsLimitPol1overPol2[3][0] = geant3_fit1_par3_NA60_pol[1];
       paramNA60Pol2[iTail].paramsLimitPol1overPol2[3][1] = geant3_fit1_par3_NA60_pol[2];

       paramNA60Pol2[iTail].paramsInitialPol1overPol2[4] = geant3_fit1_par4_NA60_pol[0];
       paramNA60Pol2[iTail].paramsLimitPol1overPol2[4][0] = geant3_fit1_par4_NA60_pol[1];
       paramNA60Pol2[iTail].paramsLimitPol1overPol2[4][1] = geant3_fit1_par4_NA60_pol[2];
       //paramNA60Pol2[iTail].paramsLimitPol2overPol3[5][0] = -100;
       //paramNA60Pol2[iTail].paramsLimitPol2overPol3[5][1] =  100;
       //paramNA60Pol2[iTail].paramsLimitPol2overPol3[6][0] = -100;
       //paramNA60Pol2[iTail].paramsLimitPol2overPol3[6][1] =  100;
       paramNA60Pol2[iTail].mean_jpsi = geant3_fit1_mean_jpsi_NA60_pol[0];
       paramNA60Pol2[iTail].mean_LowerLimits_jpsi = geant3_fit1_mean_jpsi_NA60_pol[1];
       paramNA60Pol2[iTail].mean_UpperLimits_jpsi = geant3_fit1_mean_jpsi_NA60_pol[2];
  
       paramNA60Pol2[iTail].width_jpsi = geant3_fit1_width_jpsi_NA60_pol[0];
       paramNA60Pol2[iTail].width_LowerLimits_jpsi = geant3_fit1_width_jpsi_NA60_pol[1];
       paramNA60Pol2[iTail].width_UpperLimits_jpsi = geant3_fit1_width_jpsi_NA60_pol[2];

       paramNA60Pol2[iTail].Norm_psi2S = geant3_fit1_norm_psi2S_NA60_pol[0];
       paramNA60Pol2[iTail].LowerLimits_psi2S = geant3_fit1_norm_psi2S_NA60_pol[1];
       paramNA60Pol2[iTail].UpperLimits_psi2S = geant3_fit1_norm_psi2S_NA60_pol[2];
     }
*/
      const int numOfPol2NA60_Param = numOfPol2overPol3_Parm + numOfNA60_Parm;
//      for (int iCent=0; iCent < NUMOFCENT_BINS; iCent++)
      for (int iCent=0; iCent < 1; iCent++)
      {

        fitCanvas_NA60background[iCount1]->cd();
        TF1 *fitFuncNA60Pol2 =  new TF1(paramNA60Pol2[iTail].strNameOfFitFunc.Data(), 
                                       fitFunctionDoubleNA60POL2overPOL3, 
                                       FitRange[iRange][0], 
                                       FitRange[iRange][1],
                                       numOfPol2NA60_Param+1); 

        TF1 *fitPol2reject = new TF1("fitPol2reject", 
                              BackgroundPOL2overPOL3reject_JPsi_Psi2S,
                              FitRange[iRange][0],
                              FitRange[iRange][1],
                              numOfPol2overPol3_Parm); 

        TF1 *fitPol2 = new TF1("fitPol2", 
                              BackgroundPOL2overPOL3, 
                              FitRange[iRange][0], 
                              FitRange[iRange][1], 
                              numOfPol2overPol3_Parm);

        TF1 *fitNA60_JPsi = new TF1("fitNA60_JPsi", NA60, FitRange[iRange][0], FitRange[iRange][1], numOfNA60_Parm);

        TF1 *fitNA60_Psi2S = new TF1("fitNA60_Psi2S", NA60, FitRange[iRange][0], FitRange[iRange][1], numOfNA60_Parm);

/*
        outcome_fit = fitJPsi_Psi2S(iCent,
                FitRange[iRange][0],
                FitRange[iRange][1], 
//                histoRawInvM,
                histoInvMass[pt_bin],
//                histoInvMass_Pt12to15,
                fitPol2reject,
                fitFuncNA60Pol2,
                fitPol2,
                fitNA60_JPsi,
                fitNA60_Psi2S,
                paramNA60Pol2[iTail],
                Njpsi_NA60Pol2[iRange][iTail],
                SigmaJpsi_NA60Pol2[iRange][iTail],
                MassJpsi_NA60Pol2[iRange][iTail],
                ChiJpsi_NA60Pol2[iRange][iTail],
                ErrorOfNjpsi_NA60Pol2[iRange][iTail],
                ErrorOfSigmaJpsi_NA60Pol2[iRange][iTail],
                ErrorOfMassJpsi_NA60Pol2[iRange][iTail],
                ErrorOfChiJpsi_NA60Pol2[iRange][iTail],
                fitCanvas_NA60background[iCount1],
//                fitCanvas_background02[iCount1],
                outcome_fit);

          fit_valide_NA60Pol[iRange][iTail][0] = outcome_fit.fit_valide;
          fit_status_NA60Pol[iRange][iTail][0] = outcome_fit.fit_status;
          fit_cov_matrix_status_NA60Pol[iRange][iTail][0] = outcome_fit.fit_cov_matrix_status;
*/
                iCount1+=1;      

//        delete fitFuncNA60VWG;
//        delete fitPol2reject;
//        delete fitPol2;
//        delete fitNA60;
     }
   }
  }

  histoRawInvM->SetAxisRange(rangeDo, rangeUp,"X");

  TString saveName_2[3]={"FitMethodsHistogram.pdf(","FitMethodsHistogram.pdf","FitMethodsHistogram.pdf)"};
  for(int ifit=0; ifit<FIT_NUMBER_CB2; ifit++) 
  {
    if(ifit==0)   fitCanvas_CB2background[ifit]->Print(saveName_2[0].Data());
//    else if(ifit==(FIT_NUMBER_CB2-1))   fitCanvas[ifit]->Print(saveName_2[2].Data());
    else fitCanvas_CB2background[ifit]->Print(saveName_2[1].Data());
  }

  for(int ifit=0; ifit<FIT_NUMBER_NA60; ifit++) 
  {
//    if(ifit==0)       fitCanvas_NA60[ifit]  [ifit]->Print(saveName_2[0].Data());
    if(ifit==(FIT_NUMBER_NA60-1))   fitCanvas_NA60background[ifit]->Print(saveName_2[2].Data());
    else fitCanvas_NA60background[ifit]->Print(saveName_2[1].Data());
  }



  //For saving each fit background histogram.
  TString saveName_backgroundfit[3]={"FitBackGroundHistogram.pdf(","FitBackGroundHistogram.pdf","FitBackGroundHistogram.pdf)"};
/*  for(int ifit=0; ifit<FIT_NUMBER_CB2; ifit++) 
  {
    if(ifit==0)   fitCanvas_background[ifit]->Print(saveName_backgroundfit[0].Data());
//    else if( ifit==(FIT_NUMBER_CB2-1) ) fitCanvas_background[ifit]->Print(saveName_backgroundfit[2].Data());
    else fitCanvas_background[ifit]->Print(saveName_backgroundfit[1].Data());
  }
  for(int ifit=0; ifit<FIT_NUMBER_NA60; ifit++) 
  {
//    if(ifit==0)       fitCanvas_NA60[ifit]  [ifit]->Print(saveName_2[0].Data());
    if(ifit==(FIT_NUMBER_NA60-1))   fitCanvas_background02[ifit]->Print(saveName_backgroundfit[2].Data());
    else fitCanvas_background02[ifit]->Print(saveName_backgroundfit[1].Data());
  }
*/
  
//  histo_NJpsifitUncertainty->GetXaxis()->SetLabelOffset(999);
//  histo_NJpsifitUncertainty->GetXaxis()->SetLabelSize(0);

  TString fitFuncElement_CB2[NUMOFTAILPARAM_CB2] = {"CB2+pp13","CB2+Geant3"};
  TString fitFuncElement_NA60[NUMOFTAILPARAM_NA60] = {"NA60+Geant3"};
  TString BackGroundFunc[] = {"VWG","polRatio"};
  TString fitRange[FITRANGE] = {"2.2,4.5","2.4,4.7"};

/*
  double WVALUE[max] = {2,1,1,2,1,
                        1,2,1,1,2,
                        1,1,2,2,2,
                        2,2,2,2,2}; 
*/
/*
  double WVALUE[max] = {4,1,1,4,1,
                        1,4,1,1,4,
                        1,1,1,1,1,
                        1,1,1,1,1};
*/

  double WVALUE[] = {2,1,2,1,2,1,
                     2,1,1,1,1,1};

/*
  double WVALUE[max] = {2,1,1,2,1, 
                        1,2,1,1,2,
                        1,1,2,2,2,
                        2,2,2,2,2,
                        4,2,2,4,2,
                        2,4,4,4,4};
*/

  double numberOfJPsi[max] = {0};
  double error_numberOfJPsi[max]={0};
  double sigmaOfJPsi[max]={0};
  double error_sigmaOfJPsi[max]={0};
  double massOfJPsi[max]={0};


  for (int iRange=0; iRange < FITRANGE; iRange++ )
//  for (int iRange=0; iRange < 1; iRange++ )
  {
//    for (int iTail=0; iTail< 1; iTail++)
    for (int iTail=0; iTail< NUMOFTAILPARAM_CB2; iTail++)
    {       
//    for (int iCent=0; iCent < NUMOFCENT_BINS; iCent++)
      for (int iCent=0; iCent < 1; iCent++)
      {
       xaxis->SetBinLabel(iCount, Form("%s+%s_%s",BackGroundFunc[0].Data(),fitFuncElement_CB2[iTail].Data(), fitRange[iRange].Data()) );
       xaxis_fit_valide->SetBinLabel(iCount, Form("%s+%s_%s",BackGroundFunc[0].Data(),fitFuncElement_CB2[iTail].Data(), fitRange[iRange].Data()) );
       xaxis_fit_status->SetBinLabel(iCount, Form("%s+%s_%s",BackGroundFunc[0].Data(),fitFuncElement_CB2[iTail].Data(), fitRange[iRange].Data()) );
       xaxis_cov_matrix_status->SetBinLabel(iCount, Form("%s+%s_%s",BackGroundFunc[0].Data(),fitFuncElement_CB2[iTail].Data(), fitRange[iRange].Data()) );


       std::cout<<"count: " << iCount<<std::endl;
       histo_NJpsifitUncertainty->SetBinContent(iCount, Njpsi_CB2VWG[iRange][iTail][iCent]);
       histo_NJpsifitUncertainty->SetBinError(iCount, ErrorOfNjpsi_CB2VWG[iRange][iTail][iCent]);
       histo_SigmaJpsifitUncertainty->SetBinContent(iCount,SigmaJpsi_CB2VWG[iRange][iTail][iCent] );
       histo_SigmaJpsifitUncertainty->SetBinError(iCount,ErrorOfSigmaJpsi_CB2VWG[iRange][iTail][iCent]);
       histo_MassJpsifitUncertainty->SetBinContent(iCount, MassJpsi_CB2VWG[iRange][iTail][iCent]);
       histo_MassJpsifitUncertainty->SetBinError(iCount, ErrorOfMassJpsi_CB2VWG[iRange][iTail][iCent]);
       histo_ChiJpsifitUncertainty->SetBinContent(iCount, ChiJpsi_CB2VWG[iRange][iTail][iCent] );

//       if(fit_valide_CB2VWG[iRange][iTail][0]) histo_fit_valide->SetBinContent(iCount, 0);
//       else histo_fit_valide->SetBinContent(iCount, 1);
       histo_fit_valide->SetBinContent(iCount, fit_valide_CB2VWG[iRange][iTail][0]);
       histo_fit_status->SetBinContent(iCount, fit_status_CB2VWG[iRange][iTail][0]);
       histo_cov_matrix_status->SetBinContent(iCount, fit_cov_matrix_status_CB2VWG[iRange][iTail][0]);       

     std::cout<<"Num Jpsi: "<<Njpsi_CB2VWG[iRange][iTail][iCent];
     std::cout<<" +- "<<ErrorOfNjpsi_CB2VWG[iRange][iTail][iCent] << std::endl;

     std::cout<<"Mass Jpsi: "<<MassJpsi_CB2VWG[iRange][iTail][iCent];
     std::cout<<" +- "<<ErrorOfMassJpsi_CB2VWG[iRange][iTail][iCent] <<std::endl;

     std::cout<<"Sigma Jpsi: "<< SigmaJpsi_CB2VWG[iRange][iTail][iCent];
     std::cout<<" +- "<<ErrorOfSigmaJpsi_CB2VWG[iRange][iTail][iCent] <<std::endl;

     std::cout<<"Chi2: " << ChiJpsi_CB2VWG[iRange][iTail][iCent] <<std::endl;
     
     numberOfJPsi[iCount-1] = Njpsi_CB2VWG[iRange][iTail][iCent];
     error_numberOfJPsi[iCount-1] = ErrorOfNjpsi_CB2VWG[iRange][iTail][iCent];
     sigmaOfJPsi[iCount-1] = SigmaJpsi_CB2VWG[iRange][iTail][iCent];
     error_sigmaOfJPsi[iCount-1] = ErrorOfSigmaJpsi_CB2VWG[iRange][iTail][iCent];
     massOfJPsi[iCount-1]   = MassJpsi_CB2VWG[iRange][iTail][iCent];
     iCount++;

     total_numberOfJPsi += Njpsi_CB2VWG[iRange][iTail][iCent];
     total_error_numberOfJPsi += ErrorOfNjpsi_CB2VWG[iRange][iTail][iCent];
     total_sigmaOfJPsi  += SigmaJpsi_CB2VWG[iRange][iTail][iCent];
     total_massOfJPsi   += MassJpsi_CB2VWG[iRange][iTail][iCent];
      
//     std::cout<<"total Num Jpsi: " << total_numberOfJPsi << std::endl;
      } // end of for (int iCent=0; iCent < 1; iCent++)  
    } // end of for (int iTail=0; iTail< NUMOFTAILPARAM_CB2; iTail++)
  }  // end of for (iRange=0; iRange < 2; iRange++ )

for (int iRange=0; iRange < FITRANGE; iRange++ )
//  for (int iRange=0; iRange < 1; iRange++ )
  {
//    for (int iTail=0; iTail< 1; iTail++)
    for (int iTail=0; iTail< NUMOFTAILPARAM_CB2; iTail++)
    {       
//    for (int iCent=0; iCent < NUMOFCENT_BINS; iCent++)
      for (int iCent=0; iCent < 1; iCent++)
      {
            xaxis->SetBinLabel(iCount, Form("%s+%s_%s",BackGroundFunc[1].Data(),fitFuncElement_CB2[iTail].Data(), fitRange[iRange].Data()) );
            xaxis_fit_valide->SetBinLabel(iCount, Form("%s+%s_%s",BackGroundFunc[1].Data(),fitFuncElement_CB2[iTail].Data(), fitRange[iRange].Data()) );
            xaxis_fit_status->SetBinLabel(iCount, Form("%s+%s_%s",BackGroundFunc[1].Data(),fitFuncElement_CB2[iTail].Data(), fitRange[iRange].Data()) );
            xaxis_cov_matrix_status->SetBinLabel(iCount, Form("%s+%s_%s",BackGroundFunc[1].Data(),fitFuncElement_CB2[iTail].Data(), fitRange[iRange].Data()) );

       std::cout<<"count: " << iCount<<std::endl;
       histo_NJpsifitUncertainty->SetBinContent(iCount, Njpsi_CB2Pol2[iRange][iTail][iCent]);
       histo_NJpsifitUncertainty->SetBinError(iCount, ErrorOfNjpsi_CB2Pol2[iRange][iTail][iCent]);
       histo_SigmaJpsifitUncertainty->SetBinContent(iCount,SigmaJpsi_CB2Pol2[iRange][iTail][iCent] );
       histo_SigmaJpsifitUncertainty->SetBinError(iCount,ErrorOfSigmaJpsi_CB2Pol2[iRange][iTail][iCent]);
       histo_MassJpsifitUncertainty->SetBinContent(iCount, MassJpsi_CB2Pol2[iRange][iTail][iCent]);
       histo_MassJpsifitUncertainty->SetBinError(iCount, ErrorOfMassJpsi_CB2Pol2[iRange][iTail][iCent]);
       histo_ChiJpsifitUncertainty->SetBinContent(iCount, ChiJpsi_CB2Pol2[iRange][iTail][iCent] );

//       if(fit_valide_CB2Pol[iRange][iTail][0]) histo_fit_valide->SetBinContent(iCount, 0);
//       else histo_fit_valide->SetBinContent(iCount, 1);
       histo_fit_valide->SetBinContent(iCount, fit_valide_CB2Pol[iRange][iTail][0]);
       histo_fit_status->SetBinContent(iCount, fit_status_CB2Pol[iRange][iTail][0]);
       histo_cov_matrix_status->SetBinContent(iCount, fit_cov_matrix_status_CB2Pol[iRange][iTail][0]);       


     std::cout<<"Num Jpsi: "<<Njpsi_CB2Pol2[iRange][iTail][iCent];
     std::cout<<" +- "<<ErrorOfNjpsi_CB2Pol2[iRange][iTail][iCent] << std::endl;

     std::cout<<"Mass Jpsi: "<<MassJpsi_CB2Pol2[iRange][iTail][iCent];
     std::cout<<" +- "<<ErrorOfMassJpsi_CB2Pol2[iRange][iTail][iCent] <<std::endl;

     std::cout<<"Sigma Jpsi: "<< SigmaJpsi_CB2Pol2[iRange][iTail][iCent];
     std::cout<<" +- "<<ErrorOfSigmaJpsi_CB2Pol2[iRange][iTail][iCent] <<std::endl;

     std::cout<<"Chi2: " << ChiJpsi_CB2Pol2[iRange][iTail][iCent] <<std::endl;
 
     numberOfJPsi[iCount-1] = Njpsi_CB2Pol2[iRange][iTail][iCent];
     error_numberOfJPsi[iCount-1] = ErrorOfNjpsi_CB2Pol2[iRange][iTail][iCent];
     sigmaOfJPsi[iCount-1]=SigmaJpsi_CB2Pol2[iRange][iTail][iCent];
     error_sigmaOfJPsi[iCount-1] = ErrorOfSigmaJpsi_CB2Pol2[iRange][iTail][iCent];
     massOfJPsi[iCount-1]   = MassJpsi_CB2Pol2[iRange][iTail][iCent];

     iCount++;

     total_numberOfJPsi += Njpsi_CB2Pol2[iRange][iTail][iCent];
     total_error_numberOfJPsi += ErrorOfNjpsi_CB2Pol2[iRange][iTail][iCent];
     total_sigmaOfJPsi  += SigmaJpsi_CB2Pol2[iRange][iTail][iCent];
     total_massOfJPsi   += MassJpsi_CB2Pol2[iRange][iTail][iCent];      
     //std::cout<<"total Num Jpsi: " << total_numberOfJPsi << std::endl;
      } // end of for (int iCent=0; iCent < 1; iCent++)  
    } // end of for (int iTail=0; iTail< NUMOFTAILPARAM_CB2; iTail++)
  }  // end of for (iRange=0; iRange < 2; iRange++ )



 for (int iRange=0; iRange < FITRANGE; iRange++ )
//  for (int iRange=0; iRange < 1; iRange++ )
  {
//    for (int iTail=0; iTail< 1; iTail++)
    for (int iTail=0; iTail< NUMOFTAILPARAM_NA60; iTail++)
    {       
//    for (int iCent=0; iCent < NUMOFCENT_BINS; iCent++)
      for (int iCent=0; iCent < 1; iCent++)
      {
            xaxis->SetBinLabel(iCount, Form("%s+%s_%s",BackGroundFunc[0].Data(),fitFuncElement_NA60[iTail].Data(), fitRange[iRange].Data()) );
            xaxis_fit_valide->SetBinLabel(iCount, Form("%s+%s_%s",BackGroundFunc[0].Data(),fitFuncElement_NA60[iTail].Data(), fitRange[iRange].Data()) );
            xaxis_fit_status->SetBinLabel(iCount, Form("%s+%s_%s",BackGroundFunc[0].Data(),fitFuncElement_NA60[iTail].Data(), fitRange[iRange].Data()) );
            xaxis_cov_matrix_status->SetBinLabel(iCount, Form("%s+%s_%s",BackGroundFunc[0].Data(),fitFuncElement_NA60[iTail].Data(), fitRange[iRange].Data()) );

     std::cout<<"count: " << iCount<<std::endl;
     histo_NJpsifitUncertainty->SetBinContent(iCount, Njpsi_NA60VWG[iRange][iTail][iCent]);
     histo_NJpsifitUncertainty->SetBinError(iCount, ErrorOfNjpsi_NA60VWG[iRange][iTail][iCent]);
     histo_SigmaJpsifitUncertainty->SetBinContent(iCount,SigmaJpsi_NA60VWG[iRange][iTail][iCent] );
     histo_SigmaJpsifitUncertainty->SetBinError(iCount,ErrorOfSigmaJpsi_NA60VWG[iRange][iTail][iCent]);
     histo_MassJpsifitUncertainty->SetBinContent(iCount, MassJpsi_NA60VWG[iRange][iTail][iCent]);
     histo_MassJpsifitUncertainty->SetBinError(iCount, ErrorOfMassJpsi_NA60VWG[iRange][iTail][iCent]);
     histo_ChiJpsifitUncertainty->SetBinContent(iCount, ChiJpsi_NA60VWG[iRange][iTail][iCent] );

//       if(fit_valide_NA60VWG[iRange][iTail][0]) histo_fit_valide->SetBinContent(iCount, 0);
//       else histo_fit_valide->SetBinContent(iCount, 1);
       histo_fit_valide->SetBinContent(iCount, fit_valide_NA60VWG[iRange][iTail][0]);
       histo_fit_status->SetBinContent(iCount, fit_status_NA60VWG[iRange][iTail][0]);
       histo_cov_matrix_status->SetBinContent(iCount, fit_cov_matrix_status_NA60VWG[iRange][iTail][0]); 

     std::cout<<"Num Jpsi: "<<Njpsi_NA60VWG[iRange][iTail][iCent];
     std::cout<<" +- "<<ErrorOfNjpsi_NA60VWG[iRange][iTail][iCent] << std::endl;

     std::cout<<"Mass Jpsi: "<<MassJpsi_NA60VWG[iRange][iTail][iCent];
     std::cout<<" +- "<<ErrorOfMassJpsi_NA60VWG[iRange][iTail][iCent] <<std::endl;

     std::cout<<"Sigma Jpsi: "<< SigmaJpsi_NA60VWG[iRange][iTail][iCent];
     std::cout<<" +- "<<ErrorOfSigmaJpsi_NA60VWG[iRange][iTail][iCent] <<std::endl;

     std::cout<<"Chi2: " << ChiJpsi_NA60VWG[iRange][iTail][iCent] <<std::endl;

     numberOfJPsi[iCount-1] = Njpsi_NA60VWG[iRange][iTail][iCent];
     error_numberOfJPsi[iCount-1] = ErrorOfNjpsi_NA60VWG[iRange][iTail][iCent];
     sigmaOfJPsi[iCount-1]=SigmaJpsi_NA60VWG[iRange][iTail][iCent];
     error_sigmaOfJPsi[iCount-1] = ErrorOfSigmaJpsi_NA60VWG[iRange][iTail][iCent];
     massOfJPsi[iCount-1]   = MassJpsi_NA60VWG[iRange][iTail][iCent];

     iCount++;
     total_numberOfJPsi += Njpsi_NA60VWG[iRange][iTail][iCent];
     total_error_numberOfJPsi += ErrorOfNjpsi_NA60VWG[iRange][iTail][iCent];
     total_sigmaOfJPsi  += SigmaJpsi_NA60VWG[iRange][iTail][iCent];
     total_massOfJPsi   += MassJpsi_NA60VWG[iRange][iTail][iCent];

     std::cout<<"total Num Jpsi: " << total_numberOfJPsi << std::endl;
      } // end of for (int iCent=0; iCent < 1; iCent++)  
    } // end of for (int iTail=0; iTail< NUMOFTAILPARAM_CB2; iTail++)
  }  // end of for (iRange=0; iRange < 2; iRange++ )


for (int iRange=0; iRange < FITRANGE; iRange++ )
//  for (int iRange=0; iRange < 1; iRange++ )
  {
//    for (int iTail=0; iTail< 1; iTail++)
    for (int iTail=0; iTail< NUMOFTAILPARAM_NA60; iTail++)
    {       
//    for (int iCent=0; iCent < NUMOFCENT_BINS; iCent++)
      for (int iCent=0; iCent < 1; iCent++)
      {
       xaxis->SetBinLabel(iCount, Form("%s+%s_%s",BackGroundFunc[1].Data(),fitFuncElement_NA60[iTail].Data(), fitRange[iRange].Data()) );
       xaxis_fit_valide->SetBinLabel(iCount, Form("%s+%s_%s",BackGroundFunc[1].Data(),fitFuncElement_NA60[iTail].Data(), fitRange[iRange].Data()) );
            xaxis_fit_status->SetBinLabel(iCount, Form("%s+%s_%s",BackGroundFunc[1].Data(),fitFuncElement_NA60[iTail].Data(), fitRange[iRange].Data()) );
            xaxis_cov_matrix_status->SetBinLabel(iCount, Form("%s+%s_%s",BackGroundFunc[1].Data(),fitFuncElement_NA60[iTail].Data(), fitRange[iRange].Data()) );

     std::cout<<"count: " << iCount<<std::endl;
     histo_NJpsifitUncertainty->SetBinContent(iCount, Njpsi_NA60Pol2[iRange][iTail][iCent]);
     histo_NJpsifitUncertainty->SetBinError(iCount, ErrorOfNjpsi_NA60Pol2[iRange][iTail][iCent]);
     histo_SigmaJpsifitUncertainty->SetBinContent(iCount,SigmaJpsi_NA60Pol2[iRange][iTail][iCent] );
     histo_SigmaJpsifitUncertainty->SetBinError(iCount,ErrorOfSigmaJpsi_NA60Pol2[iRange][iTail][iCent]);
     histo_MassJpsifitUncertainty->SetBinContent(iCount, MassJpsi_NA60Pol2[iRange][iTail][iCent]);
     histo_MassJpsifitUncertainty->SetBinError(iCount, ErrorOfMassJpsi_NA60Pol2[iRange][iTail][iCent]);
     histo_ChiJpsifitUncertainty->SetBinContent(iCount, ChiJpsi_NA60Pol2[iRange][iTail][iCent] );

//       if(fit_valide_NA60Pol[iRange][iTail][0]) histo_fit_valide->SetBinContent(iCount, 0);
//       else histo_fit_valide->SetBinContent(iCount, 1);
       histo_fit_valide->SetBinContent(iCount, fit_valide_NA60Pol[iRange][iTail][0]);
       histo_fit_status->SetBinContent(iCount, fit_status_NA60Pol[iRange][iTail][0]);
       histo_cov_matrix_status->SetBinContent(iCount, fit_cov_matrix_status_NA60Pol[iRange][iTail][0]); 


     std::cout<<"Num Jpsi: "<<Njpsi_NA60Pol2[iRange][iTail][iCent];
     std::cout<<" +- "<<ErrorOfNjpsi_NA60Pol2[iRange][iTail][iCent] << std::endl;

     std::cout<<"Mass Jpsi: "<<MassJpsi_NA60Pol2[iRange][iTail][iCent];
     std::cout<<" +- "<<ErrorOfMassJpsi_NA60Pol2[iRange][iTail][iCent] <<std::endl;

     std::cout<<"Sigma Jpsi: "<< SigmaJpsi_NA60Pol2[iRange][iTail][iCent];
     std::cout<<" +- "<<ErrorOfSigmaJpsi_NA60Pol2[iRange][iTail][iCent] <<std::endl;

     std::cout<<"Chi2: " << ChiJpsi_NA60Pol2[iRange][iTail][iCent] <<std::endl;
     
     numberOfJPsi[iCount-1] = Njpsi_NA60Pol2[iRange][iTail][iCent];
     error_numberOfJPsi[iCount-1] = ErrorOfNjpsi_NA60Pol2[iRange][iTail][iCent];
     sigmaOfJPsi[iCount-1]=SigmaJpsi_NA60Pol2[iRange][iTail][iCent];
     error_sigmaOfJPsi[iCount-1] = ErrorOfSigmaJpsi_NA60Pol2[iRange][iTail][iCent];
     massOfJPsi[iCount-1]   = MassJpsi_NA60Pol2[iRange][iTail][iCent];
     iCount++;

     total_numberOfJPsi += Njpsi_NA60Pol2[iRange][iTail][iCent];
     total_error_numberOfJPsi += ErrorOfNjpsi_NA60Pol2[iRange][iTail][iCent];
     total_sigmaOfJPsi  += SigmaJpsi_NA60Pol2[iRange][iTail][iCent];
     total_massOfJPsi   += MassJpsi_NA60Pol2[iRange][iTail][iCent];
     std::cout<<"total Num Jpsi: " << total_numberOfJPsi << std::endl;
      } // end of for (int iCent=0; iCent < 1; iCent++)  
    } // end of for (int iTail=0; iTail< NUMOFTAILPARAM_CB2; iTail++)
  }  // end of for (iRange=0; iRange < 2; iRange++ )

  double averageNJPsi = total_numberOfJPsi/(iCount-1);
//  std::cout << "iCount: " << iCount << std::endl;
//  std::cout << "Average N of JPsi: "<< averageNJPsi << std::endl;
  double averageErrorNJPsi = total_error_numberOfJPsi/(iCount-1);

  double averageSigmaJPsi = total_sigmaOfJPsi/(iCount-1);
  double averageMassJPsi  = total_massOfJPsi/(iCount-1);


//  for(int i=0; i<(iCount-1); i++) std::cout << numberOfJPsi[i] << std::endl;

/*
 * 
 * RMS calculation
 *
 */

  double rms = calRMS((iCount-1), numberOfJPsi, averageNJPsi);
  double rms_sigma = calRMS((iCount-1), sigmaOfJPsi, averageSigmaJPsi);
  double rms_mass = calRMS((iCount-1), massOfJPsi, averageMassJPsi);

/*
 * 
 * weight NJpsi calculation
 *
 */
  double wAverageNJpsi = calWeightAverage((iCount-1), WVALUE, numberOfJPsi);
  std::cout << wAverageNJpsi << std::endl;     
  double wErrorNJpsi = calWeightAverage((iCount-1),  WVALUE,  error_numberOfJPsi);
   std::cout<< wErrorNJpsi << std::endl;

/*
 * 
 * weighted RMS calculation
 *
 */
  double wrms = calWeightRMS((iCount-1), WVALUE, numberOfJPsi );
  std::cout <<"weight rms:"<< wrms << std::endl;    

/*
 * 
 * weight NJpsi calculation
 *
 */

  double wAverageSigmaJpsi = calWeightAverage((iCount-1), WVALUE, sigmaOfJPsi);
  double wErrorSigmaJpsi = calWeightAverage((iCount-1),  WVALUE, error_sigmaOfJPsi);

  std::cout << "weighted mean sigma of jpsi: " << wAverageSigmaJpsi << "+/-" << wErrorSigmaJpsi << std::endl;


  output_text_file.open(Form("jpsi_number_CB2_VWG_2p2_4p5.txt"));
  output_text_file << numberOfJPsi[1] <<"\n";
  output_text_file << error_numberOfJPsi[1] << "\n";
//  output_text_file << wrms << "\n";
  output_text_file.close();

  histo_NJpsifitUncertainty->SetAxisRange(wAverageNJpsi-display_NJPsi_range,wAverageNJpsi+display_NJPsi_range,"Y");
  histo_SigmaJpsifitUncertainty->SetAxisRange(averageSigmaJPsi-display_Sigma_range, averageSigmaJPsi+display_Sigma_range,"Y");
  histo_MassJpsifitUncertainty->SetAxisRange(averageMassJPsi-display_Mass_range,averageMassJPsi+display_Mass_range,"Y");
  histo_ChiJpsifitUncertainty->SetAxisRange(0.0,6.0,"Y");

//  TCanvas *c_uncertainty=new TCanvas("c_uncertainty","",200, 10, 800, 600);
  TCanvas *c_uncertainty=new TCanvas("c_uncertainty","",0,0, 800, 1000);

//  xaxis->SetLabelSize(0.02);
  xaxis->SetLabelSize(0.04);

  c_uncertainty->SetWindowSize(1,1);
  c_uncertainty->Divide(1,4,0.00,0.00);
  c_uncertainty->cd(1);

  TVirtualPad *pad1 = c_uncertainty->GetPad(1);
  pad1->SetTopMargin(0.2);
  pad1->SetRightMargin(0.2);
  pad1->Draw();
  pad1->cd();

  histo_NJpsifitUncertainty->SetStats(false);
//  histo_NJpsifitUncertainty->SetTitleSize(0.03,"y");
  histo_NJpsifitUncertainty->SetLabelSize(0.09,"Y");
  histo_NJpsifitUncertainty->Draw();

  TLine *line_average = new TLine();
  line_average->DrawLine(0,wAverageNJpsi, (max+1), wAverageNJpsi);

  TLine *line_rms[2];
  line_rms[0] = new TLine();
  line_rms[1] = new TLine();
  line_rms[0]->SetLineStyle(7);
  line_rms[1]->SetLineStyle(7);
  line_rms[0]->DrawLine(0,wrms+wAverageNJpsi, (max+1), wrms+wAverageNJpsi);
  line_rms[1]->DrawLine(0,-1*wrms+wAverageNJpsi, (max+1), -1*wrms + wAverageNJpsi);

  TLine *line_2rms[2];
  line_2rms[0] = new TLine();
  line_2rms[1] = new TLine();
  line_2rms[0]->SetLineStyle(7);
  line_2rms[0]->SetLineColor(2);
  line_2rms[1]->SetLineStyle(7);
  line_2rms[1]->SetLineColor(2);
  line_2rms[0]->DrawLine(0,2*wrms+wAverageNJpsi, (max+1), 2*wrms+wAverageNJpsi);
  line_2rms[1]->DrawLine(0,-2*wrms+wAverageNJpsi, (max+1), -2*wrms+wAverageNJpsi);

//  TLegend *legend = new TLegend(0.60,0.80,0.80,0.90);
  TLegend *legend = new TLegend(0.70,0.10,0.90,0.2);
  legend->SetTextFont(72);
//  legend->SetTextSize(0.04);

  legend->AddEntry(line_average,"Average of Number of J/#psi","l");
  legend->AddEntry(line_rms[0],"RMS","l");
  legend->AddEntry(line_2rms[0],"2*RMS","l");  
  legend->Draw();

//"NDC"
  TPaveText* t2 = new TPaveText(0.35, 0.85, 0.85,1.0,"NDC");
  t2->AddText(0.,0.,Form("%.1f <p_{T}< %.1f GeV/c N_{J/#psi} = %.0f #pm %.0f (stat) #pm %.0f (sys)", arrayAvailablePtBins_photon[pt_bin], arrayAvailablePtBins_photon[pt_bin+1], wAverageNJpsi,wErrorNJpsi, wrms ));
//  t2->AddText(0.,0.,Form("12.0 <p_{T}< 15.0 GeV/c N_{J/#psi} = %.0f #pm %.0f (stat) #pm %.0f (sys)", wAverageNJpsi,wErrorNJpsi, wrms ));
  t2->AddText(0.,0.,Form("(stat)/N_{J/#psi} = %0.2f %; (sys)/N_{J/#psi} = %0.2f %", 100*wErrorNJpsi/wAverageNJpsi, 100*wrms/wAverageNJpsi) );
  t2->AddText(0.,0.,Form("average width %0.5f",averageSigmaJPsi) );
  t2->SetTextSize(0.05);
  t2->Draw("same");

  c_uncertainty->cd(2);

  TVirtualPad *pad2 = c_uncertainty->GetPad(2);
  pad2->SetRightMargin(0.2);
  pad2->Draw();
  pad2->cd();

  histo_SigmaJpsifitUncertainty->SetStats(false);
  histo_SigmaJpsifitUncertainty->Draw();

  TLine *line_rms_sigma[2];   
  line_rms_sigma[0] = new TLine();
  line_rms_sigma[1] = new TLine();
  line_rms_sigma[0]->SetLineStyle(7);
  line_rms_sigma[1]->SetLineStyle(7);
  line_rms_sigma[0]->DrawLine(0,rms_sigma+averageSigmaJPsi, (max+1), rms_sigma+averageSigmaJPsi);
  line_rms_sigma[1]->DrawLine(0,-1*rms_sigma+averageSigmaJPsi, (max+1), -1*rms_sigma+averageSigmaJPsi);

  TLegend *legend_sigma = new TLegend(0.70,0.80,0.80,0.90);
  legend_sigma->SetTextFont(72);
//  legend->SetTextSize(0.04);
  legend_sigma->AddEntry(line_rms_sigma[0],"#pm RMS","l");
  legend_sigma->Draw();

  c_uncertainty->cd(3);

  TVirtualPad *pad3 = c_uncertainty->GetPad(3);
  pad3->SetRightMargin(0.2);
  pad3->Draw();
  pad3->cd();

  histo_MassJpsifitUncertainty->SetStats(false);
  histo_MassJpsifitUncertainty->Draw();

  TLine *line_rms_mass[2];
  line_rms_mass[0] = new TLine();
  line_rms_mass[1] = new TLine();
  line_rms_mass[0]->SetLineStyle(7);
  line_rms_mass[1]->SetLineStyle(7);
  line_rms_mass[0]->DrawLine(0, rms_mass+averageMassJPsi, (max+1), rms_mass+averageMassJPsi );
  line_rms_mass[1]->DrawLine(0, -1*rms_mass+averageMassJPsi, (max+1), -1*rms_mass+averageMassJPsi );

  TLegend *legend_mass = new TLegend(0.70,0.10,0.80,0.20);
  legend_mass->SetTextFont(72);
//  legend->SetTextSize(0.04);
  legend_mass->AddEntry(line_rms_mass[1],"#pm RMS","l");
  legend_mass->Draw();

  c_uncertainty->cd(4);
  TVirtualPad *pad4 = c_uncertainty->GetPad(4);
  pad4->SetBottomMargin(0.4);
  pad4->SetRightMargin(0.2);
  pad4->Draw();
  pad4->cd();

  histo_ChiJpsifitUncertainty->SetStats(false);
  histo_ChiJpsifitUncertainty->SetMarkerStyle(20);
  histo_ChiJpsifitUncertainty->LabelsOption("v");
//  histo_ChiJpsifitUncertainty->SetMarkerSize(2);
  histo_ChiJpsifitUncertainty->Draw("P");

  TLine *line_chi2 = new TLine();
  line_chi2->DrawLine(0,1, (max+1), 1);

  //To do: add the boolean to reflect the fitting according to the fit status, co-matrix and valildity.
  c_uncertainty->Print("FitUncertainty.pdf");

  histo_fit_valide->SetStats(false);
  histo_fit_valide->SetMarkerStyle(20);
  histo_fit_valide->LabelsOption("v");

  histo_fit_status->SetStats(false);
  histo_fit_status->SetMarkerStyle(20);
  histo_fit_status->LabelsOption("v");

  histo_cov_matrix_status->SetStats(false);
  histo_cov_matrix_status->SetMarkerStyle(20);
  histo_cov_matrix_status->LabelsOption("v");


  TCanvas *c_fit_quality = new TCanvas("c_fit_quality","",0,0, 800, 1000);
  c_fit_quality->Divide(1,3,0.00,0.00);
  c_fit_quality->cd(1);
  histo_fit_valide->Draw("P");

  c_fit_quality->cd(2);
  histo_fit_status->Draw("P");

  c_fit_quality->cd(3);
  TVirtualPad *pad_fit_quality_3 = c_fit_quality->GetPad(3);
  pad_fit_quality_3->SetBottomMargin(0.4);
  pad_fit_quality_3->Draw();
  pad_fit_quality_3->cd();
  histo_cov_matrix_status->Draw("P");

  c_fit_quality->Print(save_fit_quality_canvas_name.Data());
//  TH1F *Ratio_InvMass = new TH1F("ratio_Mine_over_Mine_pT_SummedUp","Ratio of two invariant mass; M_{#mu#mu}; Ratio",400,0,10);
//  TH1F *Ratio_InvMass = new TH1F("ratio_Mine_over_hM_V0M_0_10_pT_INTEGRATED","Ratio of two invariant mass; M_{#mu#mu}; Ratio",400,0,10);
 
//  Laure_InvMass->SetAxisRange(1,6,"X");
//  Ratio_InvMass->Divide(InvariantMass_totalCMUL,Laure_InvMass);
//  Ratio_InvMass->SetAxisRange(0.9,1.05,"Y"); 

//  histoRawInvMCheck[0]->SetAxisRange(1,6,"X");
//  Ratio_InvMass->Divide( histoRawInvMCheck[0],histoRawInvM_Direct);
/*
  TString saveName[3] = {"CMUL_DimuonInvMass_Cent_Bins.pdf(","CMUL_DimuonInvMass_Cent_Bins.pdf","CMUL_DimuonInvMass_Cent_Bins.pdf)"};
  TCanvas *c1=new TCanvas("c1","",200, 10, 800, 600);
  c1->cd();
  TPaveText* t1[10];
  for(int iCent=0; iCent<(NUMOFCENT_BINS+1); iCent++)
  {
    t1[iCent] = new TPaveText(0.56, 0.59, 0.87, 0.87,"NDC");
//    if(iCent==9) t1[iCent]->AddText(0.,0.,Form("number of CMUL events: %i",totalEvents) );
//    else t1[iCent]->AddText(0.,0.,Form("number of CMUL events: %i",totalEvents_CentBins[iCent]) );
    t1[iCent]->AddText(0.,0.,Form("Pt range: %.1f _ %.1f",arrayPtBins[0],arrayPtBins[1]));
  }
   

  for(int iCent=0; iCent<(NUMOFCENT_BINS+1); iCent++)
  {
    if(iCent==NUMOFCENT_BINS) histoRawInvM->Draw();
    else  histoRawInvM->Draw();
    t1[iCent]->Draw();
    if(iCent==0)  c1->Print(saveName[0].Data()); 
    else if(iCent == (NUMOFCENT_BINS)) c1->Print(saveName[2].Data());
    else c1->Print(saveName[1].Data());
  }
*/
//  c1->Divide(1,2,0.001,0.001);
 // c1->cd(1);
//  histoRawInvMCheck[0]->Draw();
//  histoRawInvM_PtSummedUp->Draw();
//  histoRawInvM_Direct->Draw("same");
//  Laure_InvMass->Draw("same");
//  Ratio_InvMass->Draw();
//  histoRawInvMCheck->Draw();
//    InvariantMass_totalCMUL->Draw();
  //histoRawSinglePt->Draw();
  //histoRawSinglePt_direct->Draw("same");
//  Ratio_InvMass->Draw();
//  histoRawInvM_Direct->Draw();
//    Laure_InvMass->Draw("same");

//    c1->cd(2);
//    Ratio_InvMass->Draw();

//    c1->Print("CMUL_DimuonInvMass_Cent_0_10_01.pdf");

 
//  inputFile_CMUL->Close();
//  inputRootFile_2->Close();
//  inputFile_CMUL_direct->Close();
//  TH1D *thEvent[2];
    return 0;
}




struct fit_variable_results fitJPsi_Psi2S(int iNum, double RangeMin, double RangeMax, TH1F* histo, 
             TF1* backgroundFuncReject, TF1* fitFunc, TF1* bgFunc, TF1* JPsiFunc, TF1* Psi_2SFunc,
             struct paramOfFunc parameters,
             double *NumOfJPsi, double *SigmaOfJpsi, double *MassOfJpsi, double *ChiOfJpsi,
             double *ErrorOfNjpsi, double *ErrorOfSigmaJpsi, double *ErrorOfMassJpsi, double*ErrorOfChiJpsi, 
//             TCanvas *tca, TCanvas *tca_background, struct fit_variable_results output_fit_results)
             TCanvas *tca, struct fit_variable_results output_fit_results)
{
//parameters in this functions:
//INPUT
//

//RangeMin, RangeMax: The minimum and maximum of Fit range 
//histo: The dimuon invariant mass distribution
//backgroundFuncReject:
//fitFunc:
//bgFunc:
//JPsiFunc
//Psi_2SFunc
//struct paramOfFunc parameters
//
//OUTPUT
//double *NumOfJPsi
//double *SigmaOfJpsi
//double *MassOfJpsi
//double *ChiOfJpsi
//double *ErrorOfNjpsi
//double *ErrorOfSigmaJpsi
//double *ErrorOfMassJpsi
// double*ErrorOfChiJpsi

  int totalFuncParams = parameters.NumOfBackGroundFuncParams + parameters.NumOfJPsiSignalFuncParams;

  fitFunc->SetLineColor(4);

  for(int iParams = 0; iParams < parameters.NumOfBackGroundFuncParams; iParams++)
  { fitFunc->SetParName(iParams, (parameters.NameOfBackGroundFuncParams[iParams]).Data() );}
 
  for(int iParams = parameters.NumOfBackGroundFuncParams; iParams < totalFuncParams; iParams++)
  {fitFunc->SetParName(iParams, (parameters.NameOfJPsiSignalFuncParams[iParams-parameters.NumOfBackGroundFuncParams]).Data() );}

  double digitsCritirea = 1000.;

  std::cout << "maximum of histogram: " << histo->GetMaximum() << std::endl;

  if(parameters.NumOfBackGroundFuncParams==7)
  {
//for pol2/pol3 each parameter has the different ParLimits.
//at the begining do the large values initially.
    backgroundFuncReject->SetParameter(0, parameters.paramsInitialPol2overPol3[0]*histo->GetMaximum());
    backgroundFuncReject->SetParLimits(0, parameters.paramsLimitPol2overPol3[0][0], parameters.paramsLimitPol2overPol3[0][1]*histo->GetMaximum());

    for (int iParam = 1; iParam< ( parameters.NumOfBackGroundFuncParams); iParam++)
    {
      backgroundFuncReject->SetParameter(iParam,parameters.paramsInitialPol2overPol3[iParam]);
      backgroundFuncReject->SetParLimits(iParam,parameters.paramsLimitPol2overPol3[iParam][0],parameters.paramsLimitPol2overPol3[iParam][1]);
    }

//    backgroundFuncReject->SetParameter(2,1);
//    backgroundFuncReject->SetParLimits(2,parameters.paramsLimitPol1overPol2[2][0],parameters.paramsLimitPol1overPol2[2][1]); 
//    backgroundFuncReject->SetParameter(3,1);
//    backgroundFuncReject->SetParLimits(3,parameters.paramsLimitPol1overPol2[3][0],parameters.paramsLimitPol1overPol2[3][1]);  
//    backgroundFuncReject->SetParameter(4,1);
//   backgroundFuncReject->SetParLimits(4,parameters.paramsLimitPol1overPol2[4][0],parameters.paramsLimitPol1overPol2[4][1]); 
    //backgroundFuncReject->SetParameter(5,1);
    //backgroundFuncReject->SetParLimits(5,parameters.paramsLimitPol1overPol2[5][0],parameters.paramsLimitPol1overPol2[5][1]);
    //backgroundFuncReject->SetParameter(6,1);
    //backgroundFuncReject->SetParLimits(6,parameters.paramsLimitPol1overPol2[6][0],parameters.paramsLimitPol1overPol2[6][1]);  
    histo->Fit("fitPol2reject", Form("%s",opt.Data()) ); 

    double parPol2overPol3[7]={0.}; 
    backgroundFuncReject->GetParameters(&parPol2overPol3[0]);

//    tca_background->cd();
//    histo->Draw();
//    backgroundFuncReject->Draw("SAME");

    TPaveText* t1 = new TPaveText(0.56, 0.59, 0.87, 0.87,"NDC");
    t1->AddText(0.,0.,Form("%s_%.2f_%.2f",parameters.strNameOfFitFunc.Data(),RangeMin, RangeMax));
    if( parPol2overPol3[0] >digitsCritirea || parPol2overPol3[0]< (-1*digitsCritirea) ) t1->AddText(0.,0.,Form("Background Normalization factor: %.0f #pm %.0f ", parPol2overPol3[0], backgroundFuncReject->GetParError(0) ));
    else t1->AddText(0.,0.,Form("Background Normalization factor: %.3f #pm %.3f ", parPol2overPol3[0], backgroundFuncReject->GetParError(0) ));

    for(int i=1; i< parameters.NumOfBackGroundFuncParams;i++)
    {
      if( parPol2overPol3[i] >digitsCritirea || parPol2overPol3[i] < (-1*digitsCritirea)) t1->AddText(0.,0.,Form("a%i: %.0f #pm %.0f",i, parPol2overPol3[i], backgroundFuncReject->GetParError(i)));
      else t1->AddText(0.,0.,Form("a%i: %.3f #pm %.3f",i, parPol2overPol3[i], backgroundFuncReject->GetParError(i)));
    }
//    t1->Draw();

//    tca_background->Print("background_fit.pdf");

    fitFunc->SetParameter(0,(parPol2overPol3[0]+ parameters.paramsVariancePol2overPol3[0]  ));
    fitFunc->SetParLimits(0,parameters.paramsLimitPol2overPol3[0][0],parameters.paramsLimitPol2overPol3[0][1]*histo->GetMaximum());
    for (int iParam = 1; iParam< ( parameters.NumOfBackGroundFuncParams); iParam++)
    {
        fitFunc->SetParameter(iParam,parPol2overPol3[iParam] + parameters.paramsVariancePol2overPol3[iParam]);
        fitFunc->SetParLimits(iParam,
                              (parameters.paramsLimitPol2overPol3[iParam][0] + parameters.paramsLimitVariancePol2overPol3[iParam][0]),  
                              (parameters.paramsLimitPol2overPol3[iParam][1] + parameters.paramsLimitVariancePol2overPol3[iParam][1])  );
    }
    //fitFunc->SetParameter(2,parPol1overPol2[2]);
    //fitFunc->SetParLimits(2,parameters.paramsLimitPol1overPol2[2][0],parameters.paramsLimitPol1overPol2[2][1]); 
    //fitFunc->SetParameter(3,parPol1overPol2[3]);
    //fitFunc->SetParLimits(3,parameters.paramsLimitPol1overPol2[3][0],parameters.paramsLimitPol1overPol2[3][1]);  
    //fitFunc->SetParameter(4,parPol1overPol2[4]);
    //fitFunc->SetParLimits(4,parameters.paramsLimitPol1overPol2[4][0],parameters.paramsLimitPol1overPol2[4][1]); 
    //fitFunc->SetParameter(5,parPol2overPol3[5]);
    //fitFunc->SetParLimits(5,parameters.paramsLimitPol1overPol2[5][0],parameters.paramsLimitPol1overPol2[5][1]);
    //fitFunc->SetParameter(6,parPol2overPol3[6]);
    //fitFunc->SetParLimits(6,parameters.paramsLimitPol1overPol2[6][0],parameters.paramsLimitPol1overPol2[6][1]);

    fitFunc->SetParameter(parameters.NumOfBackGroundFuncParams, histo->GetMaximum());
    fitFunc->SetParLimits(parameters.NumOfBackGroundFuncParams, 0, parameters.paramsLimitPol2overPol3[0][1]*histo->GetMaximum()); 
 }

 
  if(parameters.NumOfBackGroundFuncParams==5 )
  {
//    std::cout<< "VWG fit: " <<std::endl;

//    backgroundFuncReject->SetParameter(0,(parameters.paramsLimitVWG1[0][1])*histo->GetMaximum());
    backgroundFuncReject->SetParameter(0,(parameters.paramsInitialVWG[0])*histo->GetMaximum());
    backgroundFuncReject->SetParLimits(0,parameters.paramsLimitVWG[0][0],parameters.paramsLimitVWG[0][1]*histo->GetMaximum());

    for (int iParam = 1; iParam< ( parameters.NumOfBackGroundFuncParams); iParam++)
    {
        backgroundFuncReject->SetParameter(iParam,parameters.paramsInitialVWG[iParam]);
        backgroundFuncReject->SetParLimits(iParam,parameters.paramsLimitVWG[iParam][0],parameters.paramsLimitVWG[iParam][1]); 
    }
/*
    backgroundFuncReject->SetParameter(1,parameters.paramsinitialVWG[1]);
    backgroundFuncReject->SetParLimits(1,parameters.paramsLimitVWG[1][0],parameters.paramsLimitVWG[1][1]);

    backgroundFuncReject->SetParameter(2,parameters.paramsinitialVWG[2]);
    backgroundFuncReject->SetParLimits(2,parameters.paramsLimitVWG[2][0],parameters.paramsLimitVWG[2][1]); 

    backgroundFuncReject->SetParameter(3,parameters.paramsinitialVWG[3]);
    backgroundFuncReject->SetParLimits(3,parameters.paramsLimitVWG[3][0],parameters.paramsLimitVWG[3][1]);  
*/
//    backgroundFuncReject->SetParameter(4,1);
//    backgroundFuncReject->SetParLimits(4,parameters.paramsLimitVWG[4][0],parameters.paramsLimitVWG[4][1]); 
//    backgroundFuncReject->SetParameter(5,1);
//    backgroundFuncReject->SetParLimits(5,parameters.paramsLimitVWG[5][0],parameters.paramsLimitVWG[5][1]);
    histo->Fit("fitVWGreject", Form("%s",opt.Data() )); 

    double parVWG[5]={0.}; 
    backgroundFuncReject->GetParameters(&parVWG[0]);

//    tca_background->cd();
//    histo->Draw();
//    backgroundFuncReject->Draw("SAME");

    TPaveText* t1 = new TPaveText(0.56, 0.59, 0.87, 0.87,"NDC");
    t1->AddText(0.,0.,Form("%s_%.2f_%.2f",parameters.strNameOfFitFunc.Data(),RangeMin, RangeMax));

    if( parVWG[0] >digitsCritirea ) t1->AddText(0.,0.,Form("Background Normalization factor: %.3f #pm %.3f ", parVWG[0], backgroundFuncReject->GetParError(0) ));
    else t1->AddText(0.,0.,Form("Background Normalization factor: %.3f #pm %.3f ", parVWG[0], backgroundFuncReject->GetParError(0) ));

    for(int i=1; i< parameters.NumOfBackGroundFuncParams;i++)
    {
      if( parVWG[i] >digitsCritirea || parVWG[i] < (-1*digitsCritirea) ) t1->AddText(0.,0.,Form("a%i: %.0f #pm %.0f",i, parVWG[i], backgroundFuncReject->GetParError(i)));
      else t1->AddText(0.,0.,Form("a%i: %.3f #pm %.3f",i, parVWG[i], backgroundFuncReject->GetParError(i)));
    }

//    t1->Draw();
//    tca_background->Print("background_fit.pdf");

    fitFunc->SetParameter(0,(parVWG[0] + (parameters.paramsVarianceVWG[0])*histo->GetMaximum() ));
    fitFunc->SetParLimits(0,
                          parameters.paramsLimitVWG[0][0]*histo->GetMaximum(),
                          parameters.paramsLimitVWG[0][1]*histo->GetMaximum() );

    for (int iParam = 1; iParam< ( parameters.NumOfBackGroundFuncParams); iParam++)
    {
        fitFunc->SetParameter(iParam, (parVWG[iParam] + parameters.paramsVarianceVWG[iParam]) );
        fitFunc->SetParLimits(iParam, 
                              (parameters.paramsLimitVWG[iParam][0] + parameters.paramsLimitVarianceVWG[iParam][0]),
                              (parameters.paramsLimitVWG[iParam][1] + parameters.paramsLimitVarianceVWG[iParam][1]) );
    }
//    fitFunc->SetParameter(2,parVWG[2]);
//    fitFunc->SetParLimits(2,parameters.paramsLimitVWG1[2][0],parameters.paramsLimitVWG1[2][1]); 
//    fitFunc->SetParameter(3,parVWG[3]);
//    fitFunc->SetParLimits(3,parameters.paramsLimitVWG1[3][0],parameters.paramsLimitVWG1[3][1]);  
//    fitFunc->SetParameter(4,parVWG[4]);
//    fitFunc->SetParLimits(4,parameters.paramsLimitVWG[4][0],parameters.paramsLimitVWG[4][1]); 

    fitFunc->SetParameter(parameters.NumOfBackGroundFuncParams, histo->GetMaximum());
    fitFunc->SetParLimits(parameters.NumOfBackGroundFuncParams, 0, (parameters.paramsLimitVWG[0][1])*histo->GetMaximum()); 
//    fitFunc->SetParameter(5,parVWG[5]);
//    fitFunc->SetParLimits(5,parameters.paramsLimitVWG[5][0],parameters.paramsLimitVWG[5][1]);
  //  histo[iNum]->Fit(parameters.strNameOfFitFunc.Data());
  }


  if(parameters.NumOfBackGroundFuncParams==2 )
  {
    backgroundFuncReject->SetParameter(0,(parameters.paramsLimitExp[0][1])*histo->GetMaximum());
    backgroundFuncReject->SetParLimits(0,parameters.paramsLimitExp[0][0],parameters.paramsLimitExp[0][1]*histo->GetMaximum());
    backgroundFuncReject->SetParameter(1,1);
    backgroundFuncReject->SetParLimits(1,parameters.paramsLimitExp[1][0],parameters.paramsLimitExp[1][1]);
    backgroundFuncReject->SetParameter(2,(parameters.paramsLimitExp[2][1])*histo->GetMaximum());
    backgroundFuncReject->SetParLimits(2,parameters.paramsLimitExp[2][0],parameters.paramsLimitExp[2][1]*histo->GetMaximum());
    backgroundFuncReject->SetParameter(3,1);
    backgroundFuncReject->SetParLimits(3,parameters.paramsLimitExp[3][0],parameters.paramsLimitExp[3][1]);
 
    histo->Fit("fitEXPreject", Form("%s",opt.Data() )); 
   
    double parExp[2]={0.}; 
    backgroundFuncReject->GetParameters(&parExp[0]);

//    tca_background->cd();
//    histo->Draw();
//    backgroundFuncReject->Draw("SAME");

    TPaveText* t1 = new TPaveText(0.56, 0.59, 0.87, 0.87,"NDC");
    t1->AddText(0.,0.,Form("%s_%.2f_%.2f",parameters.strNameOfFitFunc.Data(),RangeMin, RangeMax));

    if( parExp[0] >digitsCritirea || parExp[0] < (-1*digitsCritirea) ) t1->AddText(0.,0.,Form("Background Normalization factor: %.0f #pm %.0f ", parExp[0], backgroundFuncReject->GetParError(0) ));
    else t1->AddText(0.,0.,Form("Background Normalization factor: %.3f #pm %.3f ", parExp[0], backgroundFuncReject->GetParError(0) ));

    if( parExp[1] >digitsCritirea || parExp[1] < (-1*digitsCritirea) ) t1->AddText(0.,0.,Form("a1: %.0f #pm %.0f", parExp[1], backgroundFuncReject->GetParError(1)));
    else t1->AddText(0.,0.,Form("a1: %.3f #pm %.3f", parExp[1], backgroundFuncReject->GetParError(1)));
    t1->Draw();

    fitFunc->SetParameter(0,parExp[0]);
    fitFunc->SetParLimits(0,parameters.paramsLimitExp[0][0],parameters.paramsLimitExp[0][1]*histo->GetMaximum());

    fitFunc->SetParameter(1,parExp[1]);
    fitFunc->SetParLimits(1,parameters.paramsLimitExp[1][0],parameters.paramsLimitExp[1][1]);
/*
    fitFunc->SetParameter(2,parExp[2]);
    fitFunc->SetParLimits(2,parameters.paramsLimitExp[2][0],parameters.paramsLimitExp[2][1]*histo->GetMaximum());

    fitFunc->SetParameter(3,parExp[3]);
    fitFunc->SetParLimits(3,parameters.paramsLimitExp[3][0],parameters.paramsLimitExp[3][1]); 
*/
  }

/*  
  fitFunc->SetParameter(parameters.NumOfBackGroundFuncParams+1, 3.10);
  fitFunc->SetParLimits((parameters.NumOfBackGroundFuncParams+1),3.07,3.2);
  fitFunc->SetParameter(parameters.NumOfBackGroundFuncParams+2, 0.07);
  fitFunc->SetParLimits((parameters.NumOfBackGroundFuncParams+2),0.04,0.10);
*/

  fitFunc->SetParameter(parameters.NumOfBackGroundFuncParams+1, parameters.mean_jpsi);
  fitFunc->SetParLimits((parameters.NumOfBackGroundFuncParams+1), parameters.mean_LowerLimits_jpsi, parameters.mean_UpperLimits_jpsi);
  fitFunc->SetParameter(parameters.NumOfBackGroundFuncParams+2, parameters.width_jpsi);
  fitFunc->SetParLimits((parameters.NumOfBackGroundFuncParams+2), parameters.width_LowerLimits_jpsi, parameters.width_UpperLimits_jpsi);


  
  for(int iParams = (parameters.NumOfBackGroundFuncParams+3); iParams < (totalFuncParams-1); iParams++)
  {
    fitFunc->FixParameter(iParams, parameters.ParamJPsiFunc[iParams-parameters.NumOfBackGroundFuncParams]);
//    std::cout<<"param:" << iParams;
//    std::cout<<"\t"<<parameters.ParamJPsiFunc[iParams-parameters.NumOfBackGroundFuncParams] << std::endl;
  }

//Psi(2S) normalization parameter determination
//  fitFunc->SetParameter((totalFuncParams-1), (histo->GetMaximum()/10.0) );
//  fitFunc->SetParLimits((totalFuncParams-1), 0, histo->GetMaximum()/2.0 );

  fitFunc->SetParameter((totalFuncParams-1), (histo->GetMaximum()*parameters.Norm_psi2S) );
  fitFunc->SetParLimits((totalFuncParams-1), -1*(histo->GetMaximum()*parameters.LowerLimits_psi2S), histo->GetMaximum()*parameters.UpperLimits_psi2S );

  fitFunc->SetNpx(10000);
  
  double *parFct = new double[totalFuncParams];

  TFitResultPtr tfitResult;
  tfitResult = histo->Fit( parameters.strNameOfFitFunc.Data(),Form("%s",opt.Data() ));
  fitFunc->GetParameters(&parFct[0]);

  double *parFct_excludePsi2S = new double[totalFuncParams-1];
  for(int iParam=0; iParam<(totalFuncParams-1); iParam++) parFct_excludePsi2S[iParam] = parFct[iParam];

  double *parFct_Psi2S = new double[totalFuncParams-1];
  for(int iParam = 0; iParam<(totalFuncParams-1); iParam++) parFct_Psi2S[iParam] = parFct[iParam];

  parFct_Psi2S[parameters.NumOfBackGroundFuncParams] = parFct[totalFuncParams-1];
  parFct_Psi2S[parameters.NumOfBackGroundFuncParams+1] = parFct_Psi2S[parameters.NumOfBackGroundFuncParams+1]+PDGmassDiffPsi;
  parFct_Psi2S[parameters.NumOfBackGroundFuncParams+2] = parFct_Psi2S[parameters.NumOfBackGroundFuncParams+2]*PDGmassRatio;

//  for(int iParam=0; iParam<parameters.NumOfBackGroundFuncParams; iParam++)  cout <<"parameter values: " << parFct[iParam] << "\t";
//  cout<<endl;

//  for(int iParam=0; iParam < totalFuncParams; iParam++)  std::cout <<"parameter values: " << parFct[iParam] << "\t";
//  std::cout<<std::endl;
//  for(int iParam=0; iParam < (totalFuncParams-1); iParam++)  std::cout <<"parameter values: " << parFct_Psi2S[iParam] << "\t";
//  std::cout<<std::endl;

  int FitStatus;
  FitStatus = tfitResult;
  int statusFitCov= tfitResult->CovMatrixStatus();
  
  bgFunc->SetParameters(&parFct[0]);
  bgFunc->SetLineColor(2);
  bgFunc->SetLineStyle(2);

//  for(int iParam=parameters.NumOfBackGroundFuncParams; iParam<totalFuncParams; iParam++)  cout <<"parameter values: " << parFct[iParam] << "\t";
//  cout<<endl;

  JPsiFunc->SetParameters(&parFct[parameters.NumOfBackGroundFuncParams]);
  JPsiFunc->SetLineColor(2);

//  std::cout<< JPsiFunc->GetParameter(0) << std::endl;
//  std::cout<< JPsiFunc->GetParameter(1) << std::endl;
//  std::cout<< JPsiFunc->GetParameter(2) << std::endl;
//  std::cout<< JPsiFunc->GetParameter(3) << std::endl;
//  std::cout<< JPsiFunc->GetParameter(4) << std::endl;
//  std::cout<< JPsiFunc->GetParameter(5) << std::endl;
//  std::cout<< JPsiFunc->GetParameter(6) << std::endl;


  Psi_2SFunc->SetParameters(&parFct[parameters.NumOfBackGroundFuncParams]);
  Psi_2SFunc->SetParameter(0,parFct[totalFuncParams-1]);
  Psi_2SFunc->SetParameter(1,parFct[parameters.NumOfBackGroundFuncParams+1]+PDGmassDiffPsi);
  Psi_2SFunc->SetParameter(2,parFct[parameters.NumOfBackGroundFuncParams+2]*PDGmassRatio);

  Psi_2SFunc->SetLineColor(3);

//  std::cout<< Psi_2SFunc->GetParameter(0) << std::endl;
//  std::cout<< Psi_2SFunc->GetParameter(1) << std::endl;
//  std::cout<< Psi_2SFunc->GetParameter(2) << std::endl;
//  std::cout<< Psi_2SFunc->GetParameter(3) << std::endl;
//  std::cout<< Psi_2SFunc->GetParameter(4) << std::endl;
//  std::cout<< Psi_2SFunc->GetParameter(5) << std::endl;
//  std::cout<< Psi_2SFunc->GetParameter(6) << std::endl;

//  for(int iParam=1; iParam < (parameters.NumOfJPsiSignalFuncParams+1); iParam++) Psi_2SFunc->SetParameter(iParam+1,parFct[parameters.NumOfBackGroundFuncParams+iParam]);
  

  //Get covariant matrix of the Jpsi signal
//  TMatrixDSym covMat = tfitResult->GetCovarianceMatrix().GetSub(parameters.NumOfBackGroundFuncParams, (totalFuncParams-1), parameters.NumOfBackGroundFuncParams, (totalFuncParams-1));
  TMatrixDSym covMat_JPsi = tfitResult->GetCovarianceMatrix().GetSub(parameters.NumOfBackGroundFuncParams, (totalFuncParams-2), parameters.NumOfBackGroundFuncParams, (totalFuncParams-2));
  double* covMatArray_JPsi = covMat_JPsi.GetMatrixArray();
//  covMat_JPsi.Print("FULL");
//  std::cout << covMatArray_JPsi[0] << " " << covMatArray_JPsi[1] << " " << covMatArray_JPsi[2] << " " << covMatArray_JPsi[7] << std::endl;

  NumOfJPsi[iNum]= JPsiFunc->Integral(RangeMin, RangeMax) / parameters.binWidth;
  ErrorOfNjpsi[iNum] = JPsiFunc->IntegralError(RangeMin, RangeMax, &parFct_excludePsi2S[parameters.NumOfBackGroundFuncParams], &covMatArray_JPsi[0])/parameters.binWidth;  //giving signal parameters and signal covariance matrix
//  std::cout << "Err: " << parameters.Err_Njpsi[iNum] << "\n";
  //double Err_Njpsi = JPsiFunc->IntegralError(RangeMin, RangeMax, &parFct[parameters.NumOfBackGroundFuncParams], &covMatArray[0])/parameters.binWidth;  //giving signal parameters and signal covariance matrix
//  double chi2 = JPsiFunc->GetChisquare()/JPsiFunc->GetNDF();

  //Get covariant matrix of the Psi2S signal
  TMatrixDSym covMat_Psi2S = tfitResult->GetCovarianceMatrix().GetSub((parameters.NumOfBackGroundFuncParams+1), (totalFuncParams-1), (parameters.NumOfBackGroundFuncParams+1), (totalFuncParams-1));

  TMatrixDSym covMat_new_Psi2S(parameters.NumOfJPsiSignalFuncParams-1);
  covMat_new_Psi2S(0,0) = covMat_Psi2S((parameters.NumOfJPsiSignalFuncParams-2), (parameters.NumOfJPsiSignalFuncParams-2));

  for(int i=1; i<(parameters.NumOfJPsiSignalFuncParams-1); i++) covMat_new_Psi2S(0,i) = covMat_Psi2S((parameters.NumOfJPsiSignalFuncParams-2), (i-1));

  for(int i=1; i<(parameters.NumOfJPsiSignalFuncParams-1); i++) covMat_new_Psi2S(i,0) = covMat_Psi2S(i-1, 6);


  for(int i = 1; i<(parameters.NumOfJPsiSignalFuncParams-1); i++)
  {
    for(int j=1; j<(parameters.NumOfJPsiSignalFuncParams-1); j++) covMat_new_Psi2S(i,j) = covMat_Psi2S(i-1, j-1);
  }

  double* covMatArray_Psi2S = covMat_new_Psi2S.GetMatrixArray();

  double NumOfPsi_2S = 0;
  double ErrorOfNumOfPsi_2S = 0;

  NumOfPsi_2S = Psi_2SFunc->Integral(RangeMin, RangeMax) / parameters.binWidth;
  ErrorOfNumOfPsi_2S = Psi_2SFunc->IntegralError(RangeMin, RangeMax, &parFct_Psi2S[parameters.NumOfBackGroundFuncParams], &covMatArray_Psi2S[0])/parameters.binWidth;
 
//  covMat_new_Psi2S.Print("FULL");

  double chi2 = fitFunc->GetChisquare()/fitFunc->GetNDF();
  ChiOfJpsi[iNum]=chi2;

  double mean = tfitResult->Parameter(parameters.NumOfBackGroundFuncParams+1);
  MassOfJpsi[iNum] = mean;
  double err_mean = tfitResult->ParError(parameters.NumOfBackGroundFuncParams+1);
  ErrorOfMassJpsi[iNum] = err_mean;

  double sigma_JPsi = tfitResult->Parameter(parameters.NumOfBackGroundFuncParams+2);
  SigmaOfJpsi[iNum] = sigma_JPsi;
  double err_sigma_JPsi = tfitResult->ParError(parameters.NumOfBackGroundFuncParams+2);
  ErrorOfSigmaJpsi[iNum] = err_sigma_JPsi;

  bool Fitvalide = tfitResult->IsValid();

  double signal_over_background = JPsiFunc->Integral(mean-3.*sigma_JPsi, mean+3.*sigma_JPsi) / bgFunc->Integral(mean-3.*sigma_JPsi, mean+3.*sigma_JPsi);
  double significance = JPsiFunc->Integral(mean-3.*sigma_JPsi, mean+3.*sigma_JPsi) / sqrt(JPsiFunc->Integral(mean-3.*sigma_JPsi, mean+3.*sigma_JPsi)+ bgFunc->Integral(mean-3.*sigma_JPsi, mean+3.*sigma_JPsi));
  
  tca->cd();
//  tca->SetLogy();
  histo->SetStats(false);

  TLegend *legend = new TLegend(0.55,0.15,0.90,0.35);
  legend->SetTextFont(72);
//  legend->SetTextSize(0.04);  
  legend->AddEntry(histo,"Raw events","lpe");
  legend->AddEntry(JPsiFunc,"Fitted J/#psi events","lpe");
  legend->AddEntry(Psi_2SFunc,"Fitted #psi(2S) events","lpe");
  legend->AddEntry(bgFunc,"Fitted background events","lpe");
  legend->AddEntry(fitFunc,"Fitted J/#psi #psi (2S) and background events","lpe");

  TGaxis::SetMaxDigits(4); 

  histo->Draw();
  fitFunc->Draw("SAME");
  JPsiFunc->Draw("SAME");
  Psi_2SFunc->Draw("SAME");
  bgFunc->Draw("SAME");
  legend->Draw("SAME");

  TLatex *tlatex = new TLatex();
  tlatex->SetTextFont(42);
/*
  tlatex->DrawLatexNDC(0.15,0.34,Form("#scale[0.5]{%s }", tstr_system.Data() ));
  tlatex->DrawLatexNDC(0.15,0.27,Form("#scale[0.5]{%s GeV/c}  }", tstr_pt_bin.Data() ));
  tlatex->DrawLatexNDC(0.15,0.20,Form("#scale[0.5]{%s}", tstr_cent_bin.Data()));
  tlatex->DrawLatexNDC(0.15,0.14,"#scale[0.5]{2.5 < #it{y} < 4}");
*/

  tlatex->DrawLatexNDC(0.15,0.63,Form("#scale[0.5]{%s }", tstr_system.Data() ));
  tlatex->DrawLatexNDC(0.15,0.56,Form("#scale[0.5]{%s GeV/c}  }", tstr_pt_bin.Data() ));
  tlatex->DrawLatexNDC(0.15,0.49,Form("#scale[0.5]{%s}", tstr_cent_bin.Data() ));
  tlatex->DrawLatexNDC(0.15,0.42,Form("#scale[0.5]{%s}", tstr_y_bin.Data() ));


  TString nbevents = "N_{evts} = ";
  TString centrality = "centrality bin : ";
  TString nbJpsi = "N_{J/#psi} = ";
  TString nbprim = "N_{#psi'} = ";
  TString masse = "Mass = ";
  TString sigmaJPSI = "#sigma_{J/#psi} = ";
  TString sigmaPrim = "#sigma_{#psi'} = ";
  TString alphaCB = "#alpha_{CB} = ";
  TString nCB = "n_{CB} = ";
  TString alphaUpCB = "#alpha.up_{CB} = ";
  TString nUpCB = "n.up_{CB} = ";
  TString SoverB = "S/B_{#pm3#sigma} = ";
  TString signalPrim = "S/ #sqrt{S+B} #approx ";
  TString signal = "S/ #sqrt{S+B} = ";
  TString chi2FIT = "#chi^{2}/nDoF = ";
  TString ratio = "N_{#psi'}/N_{J/#psi} = ";


  TPaveText* t1 = new TPaveText(0.56, 0.40, 0.87, 0.85,"NDC");
  t1->AddText(0.,0.,Form("%s_%.2f_%.2f",parameters.strNameOfFitFunc.Data(),RangeMin, RangeMax));
//    t1->AddText(0.,0.,Form("%s %s",centrality.Data(),strCentBins.Data()) );
  t1->AddText(0.,0.,Form("%s%.0f #pm %.0f",nbJpsi.Data(),NumOfJPsi[iNum],ErrorOfNjpsi[iNum]));
  t1->AddText(0.,0.,Form("%s%.4f #pm %.4f",masse.Data(),mean,err_mean));
  t1->AddText(0.,0.,Form("%s%.4f #pm %.4f",sigmaJPSI.Data(), sigma_JPsi, err_sigma_JPsi));

  t1->AddText(0.,0.,Form("%s%.0f #pm %.0f",nbprim.Data(),NumOfPsi_2S, ErrorOfNumOfPsi_2S));
//  t1->AddText(0.,0.,Form("%s%.4f #pm %.4f",masse.Data(),mean,err_mean));
//  t1->AddText(0.,0.,Form("%s%.4f #pm %.4f",sigmaPrim.Data(), sigma_JPsi, err_sigma_));

  t1->AddText(0.,0.,Form("%s%.2f",chi2FIT.Data(), chi2));
  t1->AddText(0.,0.,Form("%s%.2f", SoverB.Data(), signal_over_background) );
  t1->AddText(0.,0.,Form("%s%.2f", signal.Data(), significance) );
  t1->AddText(0.,0.,Form("Fit Status = %i", FitStatus));
  if(Fitvalide){t1->AddText(0.,0.,"Fit Validity : good \n");}
  else{t1->AddText(0.,0.,"Fit Validity : bad");}
  t1->AddText(0.,0.,Form("Fit Cov Matrix Status = %i", statusFitCov));

  output_fit_results.Njpsi_CB2VWG = NumOfJPsi[iNum];
  if(Fitvalide) output_fit_results.fit_valide = 0;
  else output_fit_results.fit_valide = 1;
  output_fit_results.fit_cov_matrix_status = statusFitCov;
  output_fit_results.fit_status = FitStatus;

  if(parameters.NumOfJPsiSignalFuncParams == 8)
  { 
    double signal_parameters[5] ={0};
    for(int iParams = (parameters.NumOfBackGroundFuncParams+3); iParams < (totalFuncParams-1); iParams++)
    {
      signal_parameters[iParams-parameters.NumOfBackGroundFuncParams-3] = parameters.ParamJPsiFunc[iParams-parameters.NumOfBackGroundFuncParams];
    }
    t1->AddText(0.,0.,Form("t1: %.2f", signal_parameters[0]) );
    t1->AddText(0.,0.,Form("p1: %.2f", signal_parameters[1]) );
    t1->AddText(0.,0.,Form("t2: %.2f", signal_parameters[2]) );
    t1->AddText(0.,0.,Form("p1: %.2f", signal_parameters[3]) );
  }

  if(parameters.NumOfJPsiSignalFuncParams == 12)
  {
    double signal_parameters[9] ={0};
    for(int iParams = (parameters.NumOfBackGroundFuncParams+3); iParams < (totalFuncParams-1); iParams++)
    {
      signal_parameters[iParams-parameters.NumOfBackGroundFuncParams-3] = parameters.ParamJPsiFunc[iParams-parameters.NumOfBackGroundFuncParams];
    }
    t1->AddText(0.,0.,Form("p1: %.2f", signal_parameters[0]) );
    t1->AddText(0.,0.,Form("p2: %.2f", signal_parameters[1]) );
    t1->AddText(0.,0.,Form("p3: %.2f", signal_parameters[2]) );
    t1->AddText(0.,0.,Form("p4: %.2f", signal_parameters[3]) );
    t1->AddText(0.,0.,Form("p5: %.2f", signal_parameters[4]) );
    t1->AddText(0.,0.,Form("p6: %.2f", signal_parameters[5]) );
    t1->AddText(0.,0.,Form("t1: %.2f", signal_parameters[6]) );
    t1->AddText(0.,0.,Form("t2: %.2f", signal_parameters[7]) );
  }
  
   if(parameters.NumOfBackGroundFuncParams==7)
   {
     if( tfitResult->Parameter(0) >digitsCritirea ) t1->AddText(0.,0.,Form("Background Normalization factor: %.0f #pm %.0f ", tfitResult->Parameter(0),tfitResult->ParError(0) ));
     else t1->AddText(0.,0.,Form("Background Normalization factor: %.3f #pm %.3f ", tfitResult->Parameter(0),tfitResult->ParError(0) ));

     if(tfitResult->Parameter(1) >digitsCritirea )  t1->AddText(0.,0.,Form("a1: %.0f #pm %.0f; a2: %.3f #pm %.3f ", tfitResult->Parameter(1), tfitResult->ParError(1), tfitResult->Parameter(2),tfitResult->ParError(2) ));
     else if(tfitResult->Parameter(2) >digitsCritirea )  t1->AddText(0.,0.,Form("a1: %.3f #pm %.3f; a2: %.0f #pm %.0f ", tfitResult->Parameter(1), tfitResult->ParError(1), tfitResult->Parameter(2),tfitResult->ParError(2) ));
     else if(tfitResult->Parameter(1) >digitsCritirea && tfitResult->Parameter(2) >5000 )  t1->AddText(0.,0.,Form("a1: %.0f #pm %.0f; a2: %.0f #pm %.0f ", tfitResult->Parameter(1), tfitResult->ParError(1), tfitResult->Parameter(2),tfitResult->ParError(2) ));
     else t1->AddText(0.,0.,Form("a1: %.3f #pm %.3f; a2: %.5f #pm %.3f ", tfitResult->Parameter(1), tfitResult->ParError(1), tfitResult->Parameter(2),tfitResult->ParError(2) ));

     if(tfitResult->Parameter(3) >digitsCritirea )  t1->AddText(0.,0.,Form("b1: %.0f #pm %.0f; b2: %.3f #pm %.3f ", tfitResult->Parameter(3), tfitResult->ParError(3), tfitResult->Parameter(4),tfitResult->ParError(4) ));
     else if(tfitResult->Parameter(4) >digitsCritirea )  t1->AddText(0.,0.,Form("b1: %.3f #pm %.3f; b2: %.0f #pm %.0f ", tfitResult->Parameter(3), tfitResult->ParError(3), tfitResult->Parameter(4),tfitResult->ParError(4) ));
     else if(tfitResult->Parameter(3) >digitsCritirea && tfitResult->Parameter(4) >5000 )  t1->AddText(0.,0.,Form("b1: %.0f #pm %.0f; b2: %.0f #pm %.0f ", tfitResult->Parameter(3), tfitResult->ParError(3), tfitResult->Parameter(4),tfitResult->ParError(4) ));
     else t1->AddText(0.,0.,Form("b1: %.3f #pm %.3f; b2: %.3f #pm %.3f ", tfitResult->Parameter(3), tfitResult->ParError(3), tfitResult->Parameter(4),tfitResult->ParError(4) ));

     if( tfitResult->Parameter(5) >digitsCritirea ) t1->AddText(0.,0.,Form("b3: %.0f #pm %.0f; b4: %.3f #pm %.3f", tfitResult->Parameter(5), tfitResult->ParError(5), tfitResult->Parameter(6), tfitResult->ParError(6) ));
     else if( tfitResult->Parameter(6) >digitsCritirea ) t1->AddText(0.,0.,Form("b3: %.3f #pm %.3f; b4: %.0f #pm %.0f", tfitResult->Parameter(5), tfitResult->ParError(5), tfitResult->Parameter(6), tfitResult->ParError(6) ));
     else if(tfitResult->Parameter(5) >digitsCritirea && tfitResult->Parameter(6) >5000) t1->AddText(0.,0.,Form("b3: %.0f #pm %.0f; b4: %.0f #pm %.0f", tfitResult->Parameter(5), tfitResult->ParError(5), tfitResult->Parameter(6), tfitResult->ParError(6) ));
     else t1->AddText(0.,0.,Form("b3: %.3f #pm %.3f; b4: %.3f #pm %.3f", tfitResult->Parameter(5), tfitResult->ParError(5), tfitResult->Parameter(6), tfitResult->ParError(6) )); 
 


   }     
   else if(parameters.NumOfBackGroundFuncParams==5 )
   {
     if( tfitResult->Parameter(0) >digitsCritirea ) t1->AddText(0.,0.,Form("Background Normalization factor: %.0f #pm %.0f ", tfitResult->Parameter(0),tfitResult->ParError(0) ));
     else t1->AddText(0.,0.,Form("Background Normalization factor: %.3f #pm %.3f ", tfitResult->Parameter(0),tfitResult->ParError(0) ));
     
     if( tfitResult->Parameter(1)>digitsCritirea ) t1->AddText(0.,0.,Form("a1: %.0f #pm %.0f", tfitResult->Parameter(1), tfitResult->ParError(1) ));
     else t1->AddText(0.,0.,Form("a1: %.3f #pm %.3f", tfitResult->Parameter(1), tfitResult->ParError(1) ));

     if( tfitResult->Parameter(2)>digitsCritirea ) t1->AddText(0.,0.,Form("a2: %.0f #pm %.0f", tfitResult->Parameter(2), tfitResult->ParError(2) ));
     else t1->AddText(0.,0.,Form("a2: %.3f #pm %.3f", tfitResult->Parameter(2), tfitResult->ParError(2) ));

     if(tfitResult->Parameter(3) > digitsCritirea) t1->AddText(0.,0.,Form("a3: %.0f #pm %.0f", tfitResult->Parameter(3), tfitResult->ParError(3) ));
     else t1->AddText(0.,0.,Form("a3: %.3f #pm %.3f", tfitResult->Parameter(3), tfitResult->ParError(3) ));

     if(tfitResult->Parameter(4) >digitsCritirea ) t1->AddText(0.,0.,Form("a4: %.0f #pm %.0f", tfitResult->Parameter(4), tfitResult->ParError(4) ));
     else t1->AddText(0.,0.,Form("a4: %.3f #pm %.3f", tfitResult->Parameter(4), tfitResult->ParError(4) ));

   }
   else if(parameters.NumOfBackGroundFuncParams==2 )
   {
     if( tfitResult->Parameter(0) >digitsCritirea ) t1->AddText(0.,0.,Form("Background Normalization factor: %.0f #pm %.0f ", tfitResult->Parameter(0),tfitResult->ParError(0) ));
     else t1->AddText(0.,0.,Form("Background Normalization factor: %.3f #pm %.3f ", tfitResult->Parameter(0),tfitResult->ParError(0) ));

     if( tfitResult->Parameter(1)>digitsCritirea ) t1->AddText(0.,0.,Form("a1: %.0f #pm %.0f", tfitResult->Parameter(1), tfitResult->ParError(1) ));
     else t1->AddText(0.,0.,Form("a1: %.3f #pm %.3f", tfitResult->Parameter(1), tfitResult->ParError(1) ));
     
//    t1->AddText(0.,0.,Form("a2:%.3f ", tfitResult->Parameter(2) ));
   }    
  t1->Draw("same");
//  tca->Print(Form("DimuonJPsiSignal%s_%.2f_%.2f.pdf",parameters.strNameOfFitFunc.Data(),RangeMin, RangeMax));
std::cout<< Form("JPsiSignal%s_%.2f_%.2f",parameters.strNameOfFitFunc.Data(),RangeMin, RangeMax) << std::endl;

return output_fit_results;
}

double calRMS(int index, double *value, double average)
{
  double minus_square = 0;
  for(int i=0; i<index; i++) minus_square += (average - value[i])*(average - value[i]);
  double rms = sqrt( minus_square / index);
  std::cout << "RMS: " << rms << std::endl;
return rms;
}

double calWeightAverage(int index, double *weight, double *value)
{
  double term1 = 0.0;
  double term2 = 0.0;
  double average = 0.0;
  for(int i=0; i<index; i++) term1 += weight[i] * value[i];
  for(int i=0; i<index; i++) term2 += weight[i];
  
  average = term1/term2;
  return average;
}

double calWeightRMS(int index, double *weight, double *value)
{
  double term1 = 0.0;
  double term2 = 0.0;
  double term3 = 0.0;
  double wrms=0;

  for(int i=0; i<index; i++) term1 += ( weight[i]*value[i]*value[i] );
  for(int i=0; i<index; i++) term2 += weight[i];

  for(int i=0; i<index; i++) term3 += ( weight[i]*value[i] );
  
  wrms = sqrt( (term1/term2) - (term3/term2)*(term3/term2)  );

return wrms;
}

//  thEvent[0] = eventCounters->Get("run","trigger:CMUL7-B-NOPF-MUFAST");
//  Hlist->Add(thEvent[0]);
//  thEvent[0]->GetYaxis()->SetTitle("# of events");
//  thEvent[0]->SetTitle("number of CMUL7 events as a function of runs");
//  thEvent[0]->Draw("");

  
//  thEvent[1] = eventCounters->Get("run","selected:yes/trigger:CMUL7-B-NOPF-MUFAST");
//  Hlist->Add(thEvent[1]);
//  thEvent[1]->SetLineColor(kRed);
  //thEvent[1] = eventCounters->Get("run","selected:yes/trigger:CMUL7-B-NOPF-MUFAST");
//  thEvent[1]->GetYaxis()->SetTitle("# of events");
//  thEvent[1]->Draw("same");  

  
//  thEvent[2] = (TH1D*)thEvent[0]->Clone();
//  thEvent[2]->Divide(thEvent[1], thEvent[0],1,1);
//  Hlist->Add(thEvent[2]);
//  thEvent[2]->Draw();

//  TFile outputfile("anaFile.root","recreate");
//  Hlist->Write();
//  outputfile.Close();
//  c1->Print("Run_Events_CMUL.eps");
  

  


/*
  TH1D *thEventRun1_01 = eventCounters[0]->Get("run","selected:yes");
  thEventRun1_01->Draw("same");

  c1->cd(2);
  TH1D *thEventRun2_00 = eventCounters[1]->Get("run","selected:no");
  thEventRun2_00->GetYaxis()->SetTitle("# of events");
  thEventRun2_00->Draw();

  TH1D *thEventRun2_01 = eventCounters[1]->Get("run","selected:yes");
  thEventRun2_01->Draw("same");

  c1->cd(3);
  TH1D *thEventRun3_00 = eventCounters[2]->Get("run","selected:no");
  thEventRun3_00->GetYaxis()->SetTitle("# of events");
  thEventRun3_00->Draw();

  TH1D *thEventRun3_01 = eventCounters[2]->Get("run","selected:yes");
  thEventRun3_01->Draw("same");

  c1->cd(4);
  TH1D *thEventRun3_00 = eventCounters[3]->Get("run","selected:no");
  thEventRun3_00->GetYaxis()->SetTitle("# of events");
  thEventRun3_00->Draw();

  TH1D *thEventRun3_01 = eventCounters[3]->Get("run","selected:yes");
  thEventRun3_01->Draw("same");
  //eventCounters->Get("run","selected:no")->Draw();
  //eventCounters->Get("run","selected:yes")->Draw("same");
  //eventCounters->Get("trigger","selected:yes")->Draw();
  //std::cout << eventCounters->GetSum("selected:yes") << std::endl;
*/
  //get dimuon pt histogram
/*
  TObjArray *listOfHisto = dynamic_cast<TObjArray*> (file->FindObjectAny("listOfHisto"));
  if (!listOfHisto) {
    printf("error with listOfHisto, exit!\n");
    return;
  }


  new TCanvas();
  TH1F *hDimuonPt = dynamic_cast<TH1F *>( listOfHisto->FindObject("hDimuonPt") );
  if (hDimuonPt) {
    hDimuonPt->Draw();
  }
  */
  
