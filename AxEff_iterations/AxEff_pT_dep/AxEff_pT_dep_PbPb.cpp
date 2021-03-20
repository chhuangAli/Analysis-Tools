#include <iostream>
#include <fstream>
#include <cmath>

#include "TF1.h"
#include "TList.h"
#include "TH1.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TAttLine.h"
#include "TLine.h"
#include "TLegend.h"
#include "TObjString.h"
#include "TLatex.h"
#include "TPaveStats.h"
#include "TStyle.h"
#include "TString.h"

#include "AliCounterCollection.h"

//#include "PtCentRangeValues.h"
//#include "FitFunction.h"
#include "MathTechnicalFunc.h" 
//#include "variables.h" 

double cal_binomial_error(double k, double N);

int main()
{  
  const int NUMOFFILES_LHC18q=130;
  const int NUMOFFILES_LHC18r=98;
  const int LHC15o_NUMOFFILES = 137;

  const int NUMOFFILES = 228;
  const int NUMOFFILES_LHC15o_18q_18r = LHC15o_NUMOFFILES+NUMOFFILES_LHC18r+NUMOFFILES_LHC18q;
//  const int NUMOFFILES_LHC15o_18q_18r = LHC15o_NUMOFFILES;
//  const int NUMOFFILES_LHC15o_18q_18r = 1;

  TString runList_location = Form("/home/chunlu/local/runList_LHC18q_r/muon_calo_pass3/LHC18q_r_runList.txt");

  int LHC18q_r_run_number[NUMOFFILES] = {0};
  TString tstr_LHC18q_r_run_number[NUMOFFILES]={};

  TString centralityBins_alicounter[] = {"m0","0_10","10_20","20_30","30_40","40_50","50_60","60_70","70_80","80_90","90_100"};
  const int MAXINDEX_CENTBINS = (int)(sizeof(centralityBins_alicounter)/sizeof(centralityBins_alicounter[0]));

  TString anaResult_file[NUMOFFILES];


  //This function is in "MathTechnicalFunc.h" 
  read_run_number_and_obtain_array(NUMOFFILES, runList_location, LHC18q_r_run_number, tstr_LHC18q_r_run_number);

  double **numOfevents_Cent_per_Run = new double*[NUMOFFILES];
  for(int iRun=0; iRun< (NUMOFFILES); iRun++) numOfevents_Cent_per_Run[iRun] = new double[MAXINDEX_CENTBINS]; //  **numOfevents_Cent_per_Run   {NUMOFFILES  {..MAXINDEX_CENTBINS.. }  NUMOFFILES}

  double **numOfRejecEvents_Cent_per_Run = new double*[NUMOFFILES];
  for(int iRun=0; iRun< (NUMOFFILES); iRun++) numOfRejecEvents_Cent_per_Run[iRun] = new double[MAXINDEX_CENTBINS];

 
  for(int iRun=0; iRun < (NUMOFFILES); iRun++)
  {
      std::cout << LHC18q_r_run_number[iRun] << std::endl;
     anaResult_file[iRun] = Form("$HOME/local/anaedPbPbData2018/LHC18q_r/muon_calo_pass3/AOD225_CMULEvent_AnaResults/DimuonResultsHistos_Run000%i.root",LHC18q_r_run_number[iRun]);
      
  }

//*******************************************
//
// obtain the run number of LHC15o
//
//*******************************************

  TString LHC15o_runList_location = Form("/home/chunlu/local/anaedPbPbData2015/runList_LHC15o/runList.txt");

  int LHC15o_run_number[LHC15o_NUMOFFILES] = {0};
  TString tstr_LHC15o_run_number[LHC15o_NUMOFFILES]={};

  read_run_number_and_obtain_array(LHC15o_NUMOFFILES, LHC15o_runList_location, LHC15o_run_number, tstr_LHC15o_run_number);

  TString LHC15o_anaResult_file[LHC15o_NUMOFFILES];
//  for(int iRun=0; iRun < LHC15o_NUMOFFILES; iRun++)  LHC15o_anaResult_file[iRun] = Form("$HOME/local/anaedPbPbData2015/CMULEvent_Results/AnalysisResults_run000%i.root", LHC15o_run_number[iRun]);
  for(int iRun=0; iRun < LHC15o_NUMOFFILES; iRun++)  LHC15o_anaResult_file[iRun] = Form("$HOME/local/anaedPbPbData2015/LHC15o/muon_calo_pass1/AOD229_CMULEvent_AnaResults/DimuonResultsHistos_Run000%i.root", LHC15o_run_number[iRun]);

  double **LHC15o_numOfevents_Cent_per_Run = new double*[LHC15o_NUMOFFILES];
  for(int iRun=0; iRun< LHC15o_NUMOFFILES; iRun++) LHC15o_numOfevents_Cent_per_Run[iRun] = new double[MAXINDEX_CENTBINS];

  double **LHC15o_numOfEvents_Cent_per_Run = new double*[LHC15o_NUMOFFILES];
  for(int iRun=0; iRun< LHC15o_NUMOFFILES; iRun++) LHC15o_numOfEvents_Cent_per_Run[iRun] = new double[MAXINDEX_CENTBINS];

//*******************************************
//
// obtain the number of events in a given trigger in LHC18q and LHC18r
//
//*******************************************

  TString triggerName = Form("CMUL7-B-NOPF-MUFAST");
  TString selection = Form("yes");
  TString counter_name = Form("DimuonHistosUnlike/listOfObjects_CMUL7/eventCounters");

  read_alicounter(NUMOFFILES, MAXINDEX_CENTBINS, LHC18q_r_run_number, centralityBins_alicounter, triggerName, selection, anaResult_file, counter_name, numOfevents_Cent_per_Run);
  selection = Form("no");
  read_alicounter(NUMOFFILES, MAXINDEX_CENTBINS, LHC18q_r_run_number, centralityBins_alicounter, triggerName, selection, anaResult_file, counter_name, numOfRejecEvents_Cent_per_Run);

  selection = Form("yes");
  for (int iTrig=0; iTrig <2; iTrig++)
  {
    if(iTrig == 0)
    {
      triggerName = Form("CMLL7-B-NOPF-MUFAST&CMUL7-B-NOPF-MUFAST");
      counter_name = Form("DimuonHistosUnlike/listOfObjects_CMULandCMLL/eventCounters");
      read_alicounter(LHC15o_NUMOFFILES, MAXINDEX_CENTBINS, LHC15o_run_number, centralityBins_alicounter, triggerName, selection, LHC15o_anaResult_file, counter_name, LHC15o_numOfevents_Cent_per_Run);
    }
    else if(iTrig == 1)
    {
      triggerName = Form("CMUL7-B-NOPF-MUFAST&!CMLL7-B-NOPF-MUFAST");
      counter_name = Form("DimuonHistosUnlike/listOfObjects_CMULnoCMLL/eventCounters");
      read_alicounter(LHC15o_NUMOFFILES, MAXINDEX_CENTBINS, LHC15o_run_number, centralityBins_alicounter, triggerName, selection, LHC15o_anaResult_file, counter_name, LHC15o_numOfEvents_Cent_per_Run);
    }

  }


//**********************************************
//
// Calculate the total number of runs in LHC15o, LHC18q and LHC18r
//
//**********************************************

  double LHC15o_numOfevents_per_Run[LHC15o_NUMOFFILES] = {0.0};

  for(int iRun=0; iRun<LHC15o_NUMOFFILES; iRun++)
  {
    for(int iCent=1; iCent< (MAXINDEX_CENTBINS-1); iCent++)  
    {
        LHC15o_numOfevents_per_Run[iRun] += (LHC15o_numOfevents_Cent_per_Run[iRun][iCent] + LHC15o_numOfEvents_Cent_per_Run[iRun][iCent]); 
//        LHC15o_numOfevents_Cent_per_Run[iRun][iCent] = (LHC15o_numOfevents_Cent_per_Run[iRun][iCent] + LHC15o_numOfEvents_Cent_per_Run[iRun][iCent]);
//        numOfRejecEvents_per_Run[iRun] += numOfRejecEvents_Cent_per_Run[iRun][iCent];
    }
  }

  double total_events=0;
  for(int iRun=0; iRun<NUMOFFILES; iRun++) total_events += LHC15o_numOfevents_per_Run[iRun];
  std::cout << "total accept events in CMUL for LHC15o: " << total_events << std::endl; 

  double numOfevents_per_Run[NUMOFFILES] = {0.0};
//  double numOfRejecEvents_per_Run[NUMOFFILES]={0};

  for(int iRun=0; iRun<NUMOFFILES; iRun++)
  {
    for(int iCent=1; iCent < (MAXINDEX_CENTBINS-1); iCent++)  
    {
//        numOfevents_per_Run[iRun] += (numOfevents_Cent_per_Run[iRun][iCent]+numOfRejecEvents_Cent_per_Run[iRun][iCent]);
        numOfevents_per_Run[iRun] += (numOfevents_Cent_per_Run[iRun][iCent]);  
//        numOfRejecEvents_per_Run[iRun] += numOfRejecEvents_Cent_per_Run[iRun][iCent];
    }
  }
  
  total_events=0;
  for(int iRun=0; iRun<NUMOFFILES; iRun++) total_events += numOfevents_per_Run[iRun];
  std::cout << "total accept events in CMUL for LHC18q+r: " << total_events << std::endl;

  double numOfevents_CentBins[MAXINDEX_CENTBINS] = {0};
  for(int iCent=0; iCent<MAXINDEX_CENTBINS; iCent++) 
  {
       for (int iRun=0; iRun<NUMOFFILES; iRun++) numOfevents_CentBins[iCent] += numOfevents_Cent_per_Run[iRun][iCent]; //LHC18q+18r
       for (int iRun=0; iRun<LHC15o_NUMOFFILES; iRun++) numOfevents_CentBins[iCent] += LHC15o_numOfevents_Cent_per_Run[iRun][iCent]; // LHC15o
  } 
  
//******************************************
//
// Define physics variables 
//
//******************************************

  int pTNBins = 2000; float pTMinBin = 0; float pTMaxBin =20;
  double value_each_pTbin = (pTMaxBin - pTMinBin) / pTNBins;

  int invMassNBins = 400;
  double invMassMinBin = 0;
  double invMassMaxBin =10;

  const int arrayAvailableCentralityBins[] = {0,10,20,30,40,50,60,70,80,90,100};
  const int arrayTuningCentralityBins[] = {0,10,20,30,40,50,60,90};
//  const int CENTRALITY_CLASS_SELECTION = 2; // for 0-10% and 10-20%
   const int CENTRALITY_CLASS_SELECTION = 7; // for 0-10%, 10-20%, 20-30%, 30-40% and 40-50%, 50-60% and 60-90%
  const int MAX_CENTRALITY_CLASS_SELECTION = (int)(sizeof(arrayAvailableCentralityBins)/sizeof(arrayAvailableCentralityBins[0]));

  TString centrality_total_range = Form("_Cent%dto%d", arrayAvailableCentralityBins[0],arrayAvailableCentralityBins[CENTRALITY_CLASS_SELECTION] );

  const int periods = 2;

  TH1F *hPt_rec[MAX_CENTRALITY_CLASS_SELECTION][NUMOFFILES_LHC15o_18q_18r];
  TH1F *hPt_gen[MAX_CENTRALITY_CLASS_SELECTION][NUMOFFILES_LHC15o_18q_18r];
  

//  const int ELEMENT_FOR_AE_PT_BIN =13;
//  float Pt_bin_for_Ae[ELEMENT_FOR_AE_PT_BIN]= {0, 0.3, 1, 2, 3, 4, 5, 6, 8, 10, 12, 15, 20};
  const int ELEMENT_FOR_AE_PT_BIN =16;
  float Pt_bin_for_Ae[ELEMENT_FOR_AE_PT_BIN]= {0,0.3,1,2,3,4,5,6,7,8,9,10,11,12,15,20};
  float Pt_bin_for_Ae_2015[]= {0,0.3,1,2,3,4,5,6,7,8,9,10,12};
  TH1F *hM_rec[MAX_CENTRALITY_CLASS_SELECTION][ELEMENT_FOR_AE_PT_BIN][NUMOFFILES_LHC15o_18q_18r] = {{0},{0},{0}};
  TH1F *hM_gen[MAX_CENTRALITY_CLASS_SELECTION][ELEMENT_FOR_AE_PT_BIN][NUMOFFILES_LHC15o_18q_18r] = {{0},{0},{0}};

  const int CENTRALITY_TUNING = 7;


//*************************************************
//
// Create histograms for saving the computational AxEff results
//
//*************************************************

    TH1F* hAe_Pt_wo_weight = new TH1F(Form("hAe_Pt_wo_weight"),Form("Ae as function of Pt %s; Pt(GeV/c); A#epsilon", centrality_total_range.Data()), (ELEMENT_FOR_AE_PT_BIN-1), Pt_bin_for_Ae); // without run number weighting
        
    for(int iCent=0; iCent<(MAX_CENTRALITY_CLASS_SELECTION-1); iCent++ )
    {

          for(int iRun=0; iRun < NUMOFFILES_LHC15o_18q_18r; iRun++)
          {
                  hPt_rec[iCent][iRun] = new TH1F(Form("hPt_phy_%i_%i", iCent, iRun),"hPt_phy",pTNBins,pTMinBin,pTMaxBin);
                  hPt_rec[iCent][iRun]->SetStats(kTRUE);
                  hPt_rec[iCent][iRun]->GetXaxis()->SetTitle("p_{T} (GeV/c)");
                  hPt_rec[iCent][iRun]->Sumw2();

                  hPt_gen[iCent][iRun] = new TH1F(Form("hPt_Gen_%i_%i", iCent, iRun ),"hPt_Gen",pTNBins,pTMinBin,pTMaxBin);
                  hPt_gen[iCent][iRun]->SetStats(kTRUE);
                  hPt_gen[iCent][iRun]->GetXaxis()->SetTitle("p_{T} (GeV/c)");
                  hPt_gen[iCent][iRun]->Sumw2();
                 
          }
          for( int iPt=0; iPt< (ELEMENT_FOR_AE_PT_BIN-1); iPt++)
          {
                  for(int iRun=0; iRun < NUMOFFILES_LHC15o_18q_18r; iRun++)
                  {
                    hM_rec[iCent][iPt][iRun] = new TH1F(Form("hMrecUnlikeSign_Cent%dto%d_Pt%gto%g_%i", arrayAvailableCentralityBins[iCent], arrayAvailableCentralityBins[iCent+1], Pt_bin_for_Ae[iPt], Pt_bin_for_Ae[iPt+1], iRun),"", invMassNBins, invMassMinBin, invMassMaxBin);
                    hM_rec[iCent][iPt][iRun]->Sumw2();

                    hM_gen[iCent][iPt][iRun] = new TH1F(Form("hMgenUnlikeSign_Cent%dto%d_Pt%gto%g_%i", arrayAvailableCentralityBins[iCent], arrayAvailableCentralityBins[iCent+1], Pt_bin_for_Ae[iPt], Pt_bin_for_Ae[iPt+1], iRun),"", invMassNBins, invMassMinBin, invMassMaxBin);
                    hM_gen[iCent][iPt][iRun]->Sumw2();
                  }
          }
    }
  
//****************************************************
//
// Load the .root files where dimuon histograms are saved for the AxEff computation and create a output file for the computational results.
//
//****************************************************

  
    for(int iRun = 0; iRun < NUMOFFILES_LHC15o_18q_18r; iRun++)
    {
        TString inputSimFile;
        TString inputSimFile2;
        TFile *inputFile;
        TFile *inputFile2;
        if(iRun < LHC15o_NUMOFFILES)
        {
            inputSimFile = Form("$HOME/local/anaedPbPbSimData2016/LHC16e2_jpsi_ttree/CINTEvent_Results_cuts/CINTEvent_AnaResults/PbPbMC_ResultsHistos_Run%i.root", LHC15o_run_number[iRun]);
            inputSimFile2 = Form("$HOME/local/anaedPbPbSimData2016/LHC16e2_plus_jpsi_ttree/CINTEvent_Results_cuts/CINTEvent_AnaResults/PbPbMC_ResultsHistos_Run%i.root", LHC15o_run_number[iRun]);
         
            inputFile = TFile::Open( inputSimFile.Data() );
            inputFile2 = TFile::Open( inputSimFile2.Data() );

//          if(inputFile->IsOpen() & inputFile2->IsOpen())
//          {
//            std::cout << inputSimFile.Data() << std::endl; 
//            std::cout << inputSimFile2.Data() << std::endl; 
//          }
        }
        else if(iRun >= LHC15o_NUMOFFILES)
        {
//          std::cout << LHC18q_r_run_number[iRun-LHC15o_NUMOFFILES] << std::endl;
            inputSimFile = Form("$HOME/local/anaedPbPbSim2018/LHC19a2_jpsi_ttree/CINTEvent_Results_cuts/CINTEvent_AnaResults/PbPbMC_ResultsHistos_Run%i.root",LHC18q_r_run_number[iRun-LHC15o_NUMOFFILES]);

            inputFile = TFile::Open( inputSimFile.Data() );   
            if(inputFile->IsOpen())
            {
              std::cout << inputSimFile.Data() << std::endl; 
            }     
        }
        else 
        {
            std::cout << " Please check the iRun value you use;" << std::endl;
            break;
        }

//***************************************
//
// Load the dimuon quantitiy histograms
//
//***************************************
          
        for(int iCent=0; iCent<(MAX_CENTRALITY_CLASS_SELECTION-2); iCent++ )
        {

            TString centrality_range = Form("_Cent%dto%d", arrayAvailableCentralityBins[iCent],arrayAvailableCentralityBins[iCent+1] );
//          std::cout << centrality_range.Data() << std::endl;
            
            for(int ipT=0; ipT< (ELEMENT_FOR_AE_PT_BIN-1); ipT++)
            {
                TString pTrange = Form("_Pt%gto%g", Pt_bin_for_Ae[ipT], Pt_bin_for_Ae[ipT+1]);
//                    std::cout << pTrange.Data() << std::endl;
                hPt_rec[iCent][iRun]->Add( (TH1F*)inputFile->Get( Form("DimuonHistosRecJpsi/Pt/histoPt%s%s", centrality_range.Data(), pTrange.Data())) );
                if(iRun < LHC15o_NUMOFFILES) hPt_rec[iCent][iRun]->Add( (TH1F*)inputFile2->Get( Form( "DimuonHistosRecJpsi/Pt/histoPt%s%s", centrality_range.Data(), pTrange.Data())) );

                hPt_gen[iCent][iRun]->Add( (TH1F*)inputFile->Get( Form("DimuonHistosGenJpsi/Pt/histoPt%s%s", centrality_range.Data(), pTrange.Data() )) );
                if(iRun < LHC15o_NUMOFFILES) hPt_gen[iCent][iRun]->Add( (TH1F*)inputFile2->Get( Form("DimuonHistosGenJpsi/Pt/histoPt%s%s", centrality_range.Data(), pTrange.Data())) );


                hM_rec[iCent][ipT][iRun]->Add( (TH1F*)inputFile->Get( Form("DimuonHistosRecJpsi/InvMass/histoInvMass%s%s", centrality_range.Data(), pTrange.Data() )) );
                if(iRun < LHC15o_NUMOFFILES) hM_rec[iCent][ipT][iRun]->Add( (TH1F*)inputFile2->Get( Form("DimuonHistosRecJpsi/InvMass/histoInvMass%s%s", centrality_range.Data(), pTrange.Data() )) );
                hM_gen[iCent][ipT][iRun]->Add( (TH1F*)inputFile->Get( Form("DimuonHistosGenJpsi/InvMass/histoInvMass%s%s", centrality_range.Data(), pTrange.Data() )) );
                if(iRun < LHC15o_NUMOFFILES) hM_gen[iCent][ipT][iRun]->Add( (TH1F*)inputFile2->Get( Form("DimuonHistosGenJpsi/InvMass/histoInvMass%s%s", centrality_range.Data(), pTrange.Data() )) );
             }
        }
        inputFile->Close();
        if(iRun < LHC15o_NUMOFFILES) inputFile2->Close();
    } // end of for(int iRun = 0; iRun < NUMOFFILES_LHC15o_18q_18r; iRun++)
  

//***********************************
//
// Define variables for computational AxEff results
//
//***********************************

    double entries_rec_per_cent[MAX_CENTRALITY_CLASS_SELECTION] = {0};

    double Ae_per_run[NUMOFFILES_LHC15o_18q_18r]={0};
    double errorAe_per_run[NUMOFFILES_LHC15o_18q_18r]={0};
    double centrality_events[NUMOFFILES_LHC15o_18q_18r]={0};


    double meta_Ae[MAX_CENTRALITY_CLASS_SELECTION][NUMOFFILES_LHC15o_18q_18r][ELEMENT_FOR_AE_PT_BIN] = { {{0}} };
    double meta_errorAe[MAX_CENTRALITY_CLASS_SELECTION][NUMOFFILES_LHC15o_18q_18r][ELEMENT_FOR_AE_PT_BIN] =  { {{0}}  };
    double meta_numOfevents_per_Run[MAX_CENTRALITY_CLASS_SELECTION][NUMOFFILES_LHC15o_18q_18r][ELEMENT_FOR_AE_PT_BIN] =  { {{0}} };


    double numerator_total_Ae = 0;
    double denominator_total_events = 0;
    double numerator_total_errorAe =0;
    double avg_Ae[2] = {0};


    const int lower_pT_integral_limit = 0;
//  const int upper_pT_integral_limit = 13; // which is 12 GeV/c
    const int upper_pT_integral_limit = (ELEMENT_FOR_AE_PT_BIN-1);

    double entries_hPt_rec = 0;
    double entries_hPt_gen = 0;

//*********************************************
//
// Compute the AxEff = Nrec / Ngen as a funtion of runs, centrality and pT.
//
//*********************************************

    for(int iCent=0; iCent<(MAX_CENTRALITY_CLASS_SELECTION-2); iCent++)
    {

      for(int iRun=0; iRun<NUMOFFILES_LHC15o_18q_18r; iRun++)
      {
//         std::cout << "iRun " << iRun << std::endl;
         for(int ipT = lower_pT_integral_limit; ipT< upper_pT_integral_limit; ipT++) 
         {
             entries_rec_per_cent[iCent] += (hM_rec[iCent][ipT][iRun]->Integral(1,400));

             entries_hPt_rec = (hM_rec[iCent][ipT][iRun]->Integral(1,400));
             entries_hPt_gen = (hM_gen[iCent][ipT][iRun]->Integral(1,400));
//           std::cout << entries_hPt_rec << " / " << entries_hPt_gen << std::endl;

             if(entries_hPt_gen != 0)
             {                      
                 numerator_total_Ae = (entries_hPt_rec / entries_hPt_gen);
                 numerator_total_errorAe = cal_binomial_error( entries_hPt_rec, entries_hPt_gen );
//             if(ipT==0) std::cout << numerator_total_Ae << " +/- " << numerator_total_errorAe << std::endl;
                 if(iRun < LHC15o_NUMOFFILES)          
                 {
//               denominator_total_events = LHC15o_numOfevents_Cent_per_Run[iRun][iCent+1];
                     denominator_total_events = LHC15o_numOfevents_per_Run[iRun];
                 }
                 else
                 {
                     denominator_total_events = numOfevents_per_Run[iRun-LHC15o_NUMOFFILES];
                 }

              }   
              else if(entries_hPt_gen == 0 )
              {
                 numerator_total_Ae = 0;
                 numerator_total_errorAe = 0;
                 denominator_total_events = 0;
              }

              meta_Ae[iCent][iRun][ipT] = numerator_total_Ae;
              meta_errorAe[iCent][iRun][ipT] = numerator_total_errorAe;
              meta_numOfevents_per_Run[iCent][iRun][ipT] = denominator_total_events;

              numerator_total_Ae = 0;
              numerator_total_errorAe = 0;
              denominator_total_events = 0;
              entries_hPt_rec = 0;
              entries_hPt_gen = 0;
        }

      }
    }

    std::cout << "=================" << std::endl;

//  std::cout << meta_Ae[0][0][0] << " +/- " << meta_errorAe[0][0][0] << ", " << meta_numOfevents_per_Run[0][0][0]  << std::endl;

    double Ae_summation_numerator = 0;
    double errorAe_summation_numerator = 0;
    double total_effective_events = 0;

    double Ae_pT_dep_per_cent[MAX_CENTRALITY_CLASS_SELECTION][ELEMENT_FOR_AE_PT_BIN]={{0}};
    double errorAe_pT_dep_per_cent[MAX_CENTRALITY_CLASS_SELECTION][ELEMENT_FOR_AE_PT_BIN]={{0}};
    double events_pT_dep_per_cent[MAX_CENTRALITY_CLASS_SELECTION][ELEMENT_FOR_AE_PT_BIN]={{0}};

//*****************************************
//
// Compute AxEff as a function of pT and centrailty by integrating over runs
//
//*****************************************

    for(int iCent=0; iCent<(MAX_CENTRALITY_CLASS_SELECTION-2); iCent++)
    {
        std::cout << arrayAvailableCentralityBins[iCent] <<" - " <<  arrayAvailableCentralityBins[iCent+1] <<"---> " << std::endl;
        for(int ipT = lower_pT_integral_limit; ipT< upper_pT_integral_limit; ipT++) 
        {
            std::cout << Pt_bin_for_Ae[ipT] <<" - " <<  Pt_bin_for_Ae[ipT+1] <<"---> " ;
            double AxEff_per_run_cent_pT = 0;
            double AxEff_error_per_run_cent_pT = 0;
            double events_per_run_cent_pT = 0;
            for(int iRun=0; iRun<NUMOFFILES_LHC15o_18q_18r; iRun++)     
            {
                AxEff_per_run_cent_pT = meta_Ae[iCent][iRun][ipT];
                AxEff_error_per_run_cent_pT = meta_errorAe[iCent][iRun][ipT];
                events_per_run_cent_pT = meta_numOfevents_per_Run[iCent][iRun][ipT];

                Ae_summation_numerator += (AxEff_per_run_cent_pT* events_per_run_cent_pT );
                errorAe_summation_numerator += (AxEff_error_per_run_cent_pT * AxEff_error_per_run_cent_pT *events_per_run_cent_pT *events_per_run_cent_pT);
                total_effective_events += (events_per_run_cent_pT);
       
            }      
            Ae_pT_dep_per_cent[iCent][ipT] = Ae_summation_numerator / total_effective_events;
            errorAe_pT_dep_per_cent[iCent][ipT] = sqrt(errorAe_summation_numerator/(total_effective_events * total_effective_events)); 

            std::cout << Ae_pT_dep_per_cent[iCent][ipT] <<" +/- " << errorAe_pT_dep_per_cent[iCent][ipT] <<  std::endl;
            Ae_summation_numerator = 0;
            errorAe_summation_numerator = 0;
            total_effective_events = 0;
        }    
  }


  Ae_summation_numerator = 0;
  errorAe_summation_numerator = 0;
  total_effective_events = 0;

  TString outputFileName = Form("final_AxEff_tuned_one_2015_2018PbPb.root");
  TFile *outputFile = new TFile(outputFileName.Data(),"recreate");

//*******************************************
//
// Create histograms for the computational AxEff results
//
//*******************************************

  TH1F* hAe_Pt_cent_dep[9];
  for(int iCent=0; iCent < 9; iCent++)
  {
      hAe_Pt_cent_dep[iCent] = new TH1F(Form("hAe_Pt_cent%ito%i", arrayAvailableCentralityBins[iCent], arrayAvailableCentralityBins[iCent+1]),Form("Ae as function of Pt cent%ito%i; Pt(GeV/c); A#epsilon", arrayAvailableCentralityBins[iCent], arrayAvailableCentralityBins[iCent+1]), (ELEMENT_FOR_AE_PT_BIN-1), Pt_bin_for_Ae);     
      hAe_Pt_cent_dep[iCent]->Sumw2();
      hAe_Pt_cent_dep[iCent]->SetStats(0);
  }

  TH1F* hAe_Pt_merged_cent[3];
  hAe_Pt_merged_cent[0] = new TH1F(Form("hAe_Pt_cent%ito%i", arrayTuningCentralityBins[0], arrayTuningCentralityBins[2]),Form("Ae as function of p_{T} cent%ito%i; p_{T} (GeV/c); A#epsilon", arrayTuningCentralityBins[0], arrayTuningCentralityBins[2]), (ELEMENT_FOR_AE_PT_BIN-1), Pt_bin_for_Ae);
  hAe_Pt_merged_cent[0]->Sumw2();
  hAe_Pt_merged_cent[0]->SetStats(0);

  hAe_Pt_merged_cent[1] = new TH1F(Form("hAe_Pt_cent%ito%i", arrayTuningCentralityBins[2], arrayTuningCentralityBins[4]),Form("Ae as function of p_{T} cent%ito%i; p_{T} (GeV/c); A#epsilon", arrayTuningCentralityBins[2], arrayTuningCentralityBins[4]), (ELEMENT_FOR_AE_PT_BIN-1), Pt_bin_for_Ae);
  hAe_Pt_merged_cent[1]->Sumw2();
  hAe_Pt_merged_cent[1]->SetStats(0);

  hAe_Pt_merged_cent[2] = new TH1F(Form("hAe_Pt_cent%ito%i", arrayTuningCentralityBins[4], arrayTuningCentralityBins[7]),Form("Ae as function of p_{T} cent%ito%i; p_{T} (GeV/c); A#epsilon", arrayTuningCentralityBins[4], arrayTuningCentralityBins[7]), (ELEMENT_FOR_AE_PT_BIN-1), Pt_bin_for_Ae);
  hAe_Pt_merged_cent[2]->Sumw2();
  hAe_Pt_merged_cent[2]->SetStats(0);

  for(int ipT = lower_pT_integral_limit; ipT< upper_pT_integral_limit; ipT++) 
  {
    for(int iCent=0; iCent < 9; iCent++)
    {
      hAe_Pt_cent_dep[iCent]->SetBinContent(ipT+1, Ae_pT_dep_per_cent[iCent][ipT]);
      hAe_Pt_cent_dep[iCent]->SetBinError(ipT+1, errorAe_pT_dep_per_cent[iCent][ipT]);
    }

  }


  const double entries_raw_jpsi_per_cert[] = {373767,236991,153202, 86369.4, 50531.2, 27016.8, 14465.6, 6767.11, 3118.41}; // This is the weight for taking weight average on AxEff. 

//**************************************************
//
//  Compute AxEff as a function of pT for 0-20%
//
//**************************************************

  for(int ipT = lower_pT_integral_limit; ipT< upper_pT_integral_limit; ipT++) 
  {
    for(int iCent = 0; iCent < 2; iCent++)
    {      
      Ae_summation_numerator += (Ae_pT_dep_per_cent[iCent][ipT] * entries_raw_jpsi_per_cert[iCent]);
      errorAe_summation_numerator += (errorAe_pT_dep_per_cent[iCent][ipT]*errorAe_pT_dep_per_cent[iCent][ipT]*entries_raw_jpsi_per_cert[iCent]*entries_raw_jpsi_per_cert[iCent]);
      total_effective_events += entries_raw_jpsi_per_cert[iCent];
    }    
    hAe_Pt_merged_cent[0]->SetBinContent(ipT+1, Ae_summation_numerator/total_effective_events);
    hAe_Pt_merged_cent[0]->SetBinError(ipT+1, sqrt(errorAe_summation_numerator/(total_effective_events * total_effective_events)) );

    std::cout << hAe_Pt_merged_cent[0]->GetBinContent(ipT+1) << " +/- " << hAe_Pt_merged_cent[0]->GetBinError(ipT+1) << std::endl;    

    Ae_summation_numerator = 0;
    errorAe_summation_numerator = 0;
    total_effective_events = 0;
  }

//**************************************************
//
//  Compute AxEff as a function of pT for 20-40%
//
//**************************************************

  std::cout << "AxEff, 20-40%" << std::endl;
  for(int ipT = lower_pT_integral_limit; ipT< upper_pT_integral_limit; ipT++) 
  {
    for(int iCent = 2; iCent < 4; iCent++)
    {      
      Ae_summation_numerator += (Ae_pT_dep_per_cent[iCent][ipT] * entries_raw_jpsi_per_cert[iCent]);
      errorAe_summation_numerator += (errorAe_pT_dep_per_cent[iCent][ipT]*errorAe_pT_dep_per_cent[iCent][ipT]*entries_raw_jpsi_per_cert[iCent]*entries_raw_jpsi_per_cert[iCent]);
      total_effective_events += entries_raw_jpsi_per_cert[iCent];
    }    
    hAe_Pt_merged_cent[1]->SetBinContent(ipT+1, Ae_summation_numerator/total_effective_events);
    hAe_Pt_merged_cent[1]->SetBinError(ipT+1, sqrt(errorAe_summation_numerator/(total_effective_events * total_effective_events)) );

    std::cout << hAe_Pt_merged_cent[1]->GetBinContent(ipT+1) << " +/- " << hAe_Pt_merged_cent[1]->GetBinError(ipT+1) << std::endl;    

    Ae_summation_numerator = 0;
    errorAe_summation_numerator = 0;
    total_effective_events = 0;
  }

//**************************************************
//
//  Compute AxEff as a function of pT for 40-90%
//
//**************************************************

  Ae_summation_numerator = 0;
  errorAe_summation_numerator = 0;
  total_effective_events = 0;
  std::cout << "AxEff, 40-90%" << std::endl;
  for(int ipT = lower_pT_integral_limit; ipT< upper_pT_integral_limit; ipT++) 
  {
    for(int iCent = 4; iCent < (MAX_CENTRALITY_CLASS_SELECTION-2); iCent++)
    {      
      Ae_summation_numerator += (Ae_pT_dep_per_cent[iCent][ipT] * entries_raw_jpsi_per_cert[iCent]);
      errorAe_summation_numerator += (errorAe_pT_dep_per_cent[iCent][ipT]*errorAe_pT_dep_per_cent[iCent][ipT]*entries_raw_jpsi_per_cert[iCent]*entries_raw_jpsi_per_cert[iCent]);
      total_effective_events += entries_raw_jpsi_per_cert[iCent];
    }    
    hAe_Pt_merged_cent[2]->SetBinContent(ipT+1, Ae_summation_numerator/total_effective_events);
    hAe_Pt_merged_cent[2]->SetBinError(ipT+1, sqrt(errorAe_summation_numerator/(total_effective_events * total_effective_events)) );

    std::cout << hAe_Pt_merged_cent[2]->GetBinContent(ipT+1) << " +/- " << hAe_Pt_merged_cent[2]->GetBinError(ipT+1) << std::endl;
    
    Ae_summation_numerator = 0;
    errorAe_summation_numerator = 0;
    total_effective_events = 0;
  }

//**************************************************
//
//  Compute AxEff as a function of pT for 0-90%
//
//**************************************************  

  Ae_summation_numerator = 0;
  errorAe_summation_numerator = 0;
  total_effective_events = 0;

  std::cout << "AxEff, integrated over centrality" << std::endl;
  for(int ipT = lower_pT_integral_limit; ipT< upper_pT_integral_limit; ipT++) 
  {
    for(int iCent = 0; iCent < (MAX_CENTRALITY_CLASS_SELECTION-2); iCent++)
    {      
      Ae_summation_numerator += (Ae_pT_dep_per_cent[iCent][ipT] * entries_rec_per_cent[iCent]);
      errorAe_summation_numerator += (errorAe_pT_dep_per_cent[iCent][ipT]*errorAe_pT_dep_per_cent[iCent][ipT]*entries_rec_per_cent[iCent]*entries_rec_per_cent[iCent]);
      total_effective_events += entries_rec_per_cent[iCent];
    }    
    std::cout << Ae_summation_numerator/total_effective_events << " +/- " << sqrt(errorAe_summation_numerator/(total_effective_events * total_effective_events)) << std::endl;
    
    Ae_summation_numerator = 0;
    errorAe_summation_numerator = 0;
    total_effective_events = 0;
  }

//***********************************************
//
// Create canvas for histograms to draw/display.
//
//***********************************************

  TCanvas *c1 = new TCanvas("c1","",800,600);
  TLine *line_average = new TLine();

  TLatex *myt[3];
  myt[0] = new TLatex(5,0.123, Form("cent 20 to 40 %%, 2.5 < y < 4 ") );
  myt[0] ->SetTextFont(42);
  myt[0] ->SetTextSize(0.02);
  myt[0] ->SetTextColor(kRed);

  myt[1] = new TLatex(5,0.123, Form("cent 40 to 90 %%, 2.5 < y < 4 ") );
  myt[1] ->SetTextFont(42);
  myt[1] ->SetTextSize(0.02);
  myt[1] ->SetTextColor(kRed);


  hAe_Pt_merged_cent[0]->Draw();
  myt[0]->Draw("SAME");
  c1->Print("AxEff_pT_dep_cent20to40.pdf");

  hAe_Pt_merged_cent[1]->Draw();
  myt[1]->Draw("SAME");
  c1->Print("AxEff_pT_dep_cent40to90.pdf");

  outputFile->Write();
  delete outputFile;


return 0;
}


double cal_binomial_error(double Nrec, double Ngen)
{
   double error_efficiency = 0;

   if( Nrec > Ngen)
   {
     error_efficiency = 1 / Ngen;
   }
   else
   {
     error_efficiency = TMath::Max(1/Ngen, sqrt((Nrec/Ngen )*( (1-( Nrec / Ngen )) / Ngen )) );
   }

return error_efficiency;
}

