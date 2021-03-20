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

#include "PtCentRangeValues.h"
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
// TString LHC18r_runList_location = Form("/home/chunlu/local/runList_LHC18r/muon_calo_pass2/LHC18r_runList.txt");

  TString runList_location = Form("/home/chunlu/local/runList_LHC18q_r/muon_calo_pass3/LHC18q_r_runList.txt");
//  TString runList_location = Form("/home/chunlu/local/runList_LHC18q_r/muon_calo_pass2/LHC18q_r_runList_without_296752.txt");


  int LHC18q_r_run_number[NUMOFFILES] = {0};
  TString tstr_LHC18q_r_run_number[NUMOFFILES]={};

  TString centralityBins_alicounter[] = {"m0","0_10","10_20","20_30","30_40","40_50","50_60","60_70","70_80","80_90","90_100"};
  const int MAXINDEX_CENTBINS = (int)(sizeof(centralityBins_alicounter)/sizeof(centralityBins_alicounter[0]));

  TString anaResult_file[NUMOFFILES];
  TH1F *histo_Ae_run;


  //#include "MathTechnicalFunc.h" 
  read_run_number_and_obtain_array(NUMOFFILES, runList_location, LHC18q_r_run_number, tstr_LHC18q_r_run_number);

//  read_run_number_and_obtain_array(NUMOFFILES_LHC18r, LHC18r_runList_location, LHC18r_run_number, tstr_LHC18r_run_number);

  double **numOfevents_Cent_per_Run = new double*[NUMOFFILES];
  for(int iRun=0; iRun< (NUMOFFILES); iRun++) numOfevents_Cent_per_Run[iRun] = new double[MAXINDEX_CENTBINS]; //  **numOfevents_Cent_per_Run   {NUMOFFILES  {..MAXINDEX_CENTBINS.. }  NUMOFFILES}

  double **numOfRejecEvents_Cent_per_Run = new double*[NUMOFFILES];
  for(int iRun=0; iRun< (NUMOFFILES); iRun++) numOfRejecEvents_Cent_per_Run[iRun] = new double[MAXINDEX_CENTBINS];

 
  for(int iRun=0; iRun < (NUMOFFILES); iRun++)
  {
      std::cout << LHC18q_r_run_number[iRun] << std::endl;
//      anaResult_file[iRun] = Form("$HOME/local/anaedPbPbData2018/LHC18q_r/muon_calo_pass2/AOD211_CMULEvent_Results/AnalysisResults_run000%i.root",LHC18q_r_run_number[iRun]);
     anaResult_file[iRun] = Form("$HOME/local/anaedPbPbData2018/LHC18q_r/muon_calo_pass3/AOD225_CMULEvent_AnaResults/DimuonResultsHistos_Run000%i.root",LHC18q_r_run_number[iRun]);
      
  }

//
//
// obtain the run number of LHC15o
//
//
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

//
//
// obtain the number of events in a given trigger in LHC18q and LHC18r
//
//
  TString triggerName = Form("CMUL7-B-NOPF-MUFAST");
  TString selection = Form("yes");
  TString counter_name = Form("DimuonHistosUnlike/listOfObjects_CMUL7/eventCounters");

  read_alicounter(NUMOFFILES, MAXINDEX_CENTBINS, LHC18q_r_run_number, centralityBins_alicounter, triggerName, selection, anaResult_file, counter_name, numOfevents_Cent_per_Run);
  
  read_alicounter(NUMOFFILES, MAXINDEX_CENTBINS, LHC18q_r_run_number, centralityBins_alicounter, triggerName, selection, anaResult_file, counter_name, numOfRejecEvents_Cent_per_Run);


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

  double LHC15o_numOfevents_per_Run[LHC15o_NUMOFFILES] = {0.0};

//  for(int iRun=0; iRun<NUMOFFILES; iRun++)  for(int iCent=0; iCent<MAXINDEX_CENTBINS; iCent++) std::cout << numOfevents_Cent_per_Run[iRun][iCent] << std::endl;

  for(int iRun=0; iRun<LHC15o_NUMOFFILES; iRun++)
  {
    for(int iCent=0; iCent<MAXINDEX_CENTBINS; iCent++)  
    {
        std::cout << LHC15o_numOfevents_Cent_per_Run[iRun][iCent] << " " << LHC15o_numOfEvents_Cent_per_Run[iRun][iCent] << ", ";
        LHC15o_numOfevents_per_Run[iRun] += (LHC15o_numOfevents_Cent_per_Run[iRun][iCent] + LHC15o_numOfEvents_Cent_per_Run[iRun][iCent]); 
        LHC15o_numOfevents_Cent_per_Run[iRun][iCent] = (LHC15o_numOfevents_Cent_per_Run[iRun][iCent] + LHC15o_numOfEvents_Cent_per_Run[iRun][iCent]);
        std::cout << LHC15o_numOfevents_Cent_per_Run[iRun][iCent] << std::endl;
//        numOfRejecEvents_per_Run[iRun] += numOfRejecEvents_Cent_per_Run[iRun][iCent];
    }
  }

  double numOfevents_per_Run[NUMOFFILES] = {0.0};
  double numOfRejecEvents_per_Run[NUMOFFILES]={0};

  for(int iRun=0; iRun<NUMOFFILES; iRun++)
  {
    for(int iCent=0; iCent<MAXINDEX_CENTBINS; iCent++)  
    {
        numOfevents_per_Run[iRun] += numOfevents_Cent_per_Run[iRun][iCent]; 
        numOfRejecEvents_per_Run[iRun] += numOfRejecEvents_Cent_per_Run[iRun][iCent];
    }
  }
  
  double total_events=0;
  for(int iRun=0; iRun<NUMOFFILES; iRun++) total_events += numOfevents_per_Run[iRun];
  std::cout << "total accept events in CMUL: " << total_events << std::endl;

  double numOfevents_CentBins[MAXINDEX_CENTBINS] = {0};
  for(int iCent=0; iCent<MAXINDEX_CENTBINS; iCent++) 
  {
       for (int iRun=0; iRun<NUMOFFILES; iRun++) numOfevents_CentBins[iCent] += numOfevents_Cent_per_Run[iRun][iCent]; //LHC18q+18r
       for (int iRun=0; iRun<LHC15o_NUMOFFILES; iRun++) numOfevents_CentBins[iCent] += LHC15o_numOfevents_Cent_per_Run[iRun][iCent]; // LHC15o
  } 
  
  int NBins = 60; float MinBin = -5; float MaxBin = -2;
  double value_each_bin = (MaxBin - MinBin) / NBins;

  const int arrayAvailableCentralityBins[] = {0,10,20,30,40,50,60,70,80,90,100};
  const int arrayTuningCentralityBins[] = {0,10,20,30,40,50,60,90};
  const int CENTRALITY_CLASS_SELECTION = 2; // for 0-10% and 10-20%
  const int MAX_CENTRALITY_CLASS_SELECTION = (int)(sizeof(arrayAvailableCentralityBins)/sizeof(arrayAvailableCentralityBins[0]));

  TString centrality_total_range = Form("_Cent%dto%d", arrayAvailableCentralityBins[0],arrayAvailableCentralityBins[CENTRALITY_CLASS_SELECTION] );

  const int periods = 2;

  TH1F *hY_rec[MAX_CENTRALITY_CLASS_SELECTION][NUMOFFILES_LHC15o_18q_18r];
  TH1F *hY_gen[MAX_CENTRALITY_CLASS_SELECTION][NUMOFFILES_LHC15o_18q_18r];
  TH1F *hY_ratio[MAX_CENTRALITY_CLASS_SELECTION][NUMOFFILES_LHC15o_18q_18r];
  
  TH1F *hY_rec_sum = new TH1F("hY_rec_sum","; Y ; # of events", NBins, MinBin, MaxBin);
  TH1F *hY_gen_sum = new TH1F("hY_gen_sum",Form("%s; Y ; # of events", centrality_total_range.Data()), NBins, MinBin, MaxBin);

  TH1F *hY_gen_minus_rec[MAX_CENTRALITY_CLASS_SELECTION];
  TH1F *hY_gen_minus_rec_with_small_rec[MAX_CENTRALITY_CLASS_SELECTION];


//  const int ELEMENT_FOR_AE_PT_BIN =13;
//  float Pt_bin_for_Ae[ELEMENT_FOR_AE_PT_BIN]= {0, 0.3, 1, 2, 3, 4, 5, 6, 8, 10, 12, 15, 20};
  const int ELEMENT_FOR_AE_PT_BIN =16;
  float Pt_bin_for_Ae[ELEMENT_FOR_AE_PT_BIN]= {0,0.3,1,2,3,4,5,6,7,8,9,10,11,12,15,20};

  float y_bin[] = {-4., -3.75, -3.5, -3.25, -3., -2.75, -2.5};
  float y_bin_pos[] = {2.5, 2.75, 3, 3.25, 3.5, 3.75, 4};
  const int number_y_bin = 6;

  const int CENTRALITY_TUNING = 7;
  TH1F* hAe_cent_dep[CENTRALITY_TUNING];

  TH1F *hM_rec[MAX_CENTRALITY_CLASS_SELECTION][number_y_bin][NUMOFFILES_LHC15o_18q_18r];
  TH1F *hM_gen[MAX_CENTRALITY_CLASS_SELECTION][number_y_bin][NUMOFFILES_LHC15o_18q_18r];

  for(int iCent = 0; iCent < (CENTRALITY_TUNING); iCent++)
  {
      hAe_cent_dep[iCent] = new TH1F(Form("hAe_cent%ito%i", arrayTuningCentralityBins[iCent], arrayTuningCentralityBins[iCent+1]),Form("Ae as function of y cent%ito%i; y; A#epsilon", arrayTuningCentralityBins[iCent], arrayTuningCentralityBins[iCent+1]), (number_y_bin), y_bin);     
      hAe_cent_dep[iCent]->Sumw2();
  }
  
  

  TH1F* hAe_Pt_wo_weight = new TH1F(Form("hAe_Pt_wo_weight"),Form("Ae as function of Pt %s; Pt(GeV/c); A#epsilon", centrality_total_range.Data()), (ELEMENT_FOR_AE_PT_BIN-1), Pt_bin_for_Ae); // without run number weighting

  TH1I* hRun_rec_greater_gen[CENTRALITY_CLASS_SELECTION];
  hRun_rec_greater_gen[0] = new TH1I("hRun_Cent0to10","Runs with reconstructed jpsi greater than generated jpsi" ,20,0,20);
  hRun_rec_greater_gen[1] = new TH1I("hRun_Cent10to20","Runs with reconstructed jpsi greater than generated jpsi",20,0,20);
  
           
//  for(int iCent=0; iCent<CENTRALITY_CLASS_SELECTION; iCent++)
  
  
    for(int iCent=0; iCent<(MAX_CENTRALITY_CLASS_SELECTION-1); iCent++ )
    {
//          hPt_gen_minus_rec[iCent] =  new TH1F(Form("hPt_gen_minus_rec_%i", iCent), "#Delta Pt(Pt) (= P_{T}^{Gen} - P_{T}^{Rec}); Pt(GeV/c); #Delta Pt", pTNBins,pTMinBin,pTMaxBin);
//          hPt_gen_minus_rec[iCent]->Sumw2();
//          hPt_gen_minus_rec_with_small_rec[iCent] = new TH1F(Form("hPt_gen_minus_rec_with_small_rec_%i",iCent), "#Delta Pt(Pt) (= P_{T}^{Gen} - P_{T}^{Rec}); Pt(GeV/c); #Delta Pt", pTNBins,pTMinBin,pTMaxBin);
//          hPt_gen_minus_rec_with_small_rec[iCent]->Sumw2();

          for(int iRun=0; iRun < NUMOFFILES_LHC15o_18q_18r; iRun++)
          {
                  hY_rec[iCent][iRun] = new TH1F(Form("hY_phy_%i_%i", iCent, iRun),"hY_phy", NBins, MinBin, MaxBin);
                  hY_rec[iCent][iRun]->SetStats(kTRUE);
                  hY_rec[iCent][iRun]->GetXaxis()->SetTitle("y");
                  hY_rec[iCent][iRun]->Sumw2();

                  hY_gen[iCent][iRun] = new TH1F(Form("hY_Gen_%i_%i", iCent, iRun ),"hY_Gen", NBins, MinBin, MaxBin);
                  hY_gen[iCent][iRun]->SetStats(kTRUE);
                  hY_gen[iCent][iRun]->GetXaxis()->SetTitle("y") ;
                  hY_gen[iCent][iRun]->Sumw2();

                  hY_ratio[iCent][iRun] = new TH1F(Form("hY_ratio_%i_%i", iCent, iRun),"; y; ", NBins, MinBin, MaxBin);
                  hY_ratio[iCent][iRun]->Sumw2();
          }

          for( int iY=0; iY<number_y_bin; iY++)
          {
                for(int iRun=0; iRun < NUMOFFILES_LHC15o_18q_18r; iRun++)
                {
                  hM_rec[iCent][iY][iRun] = new TH1F(Form("hM_phy_%i_%i_%i", iCent, iY, iRun),"hM_phy", 400, 0, 10);
                  hM_rec[iCent][iY][iRun]->Sumw2();
                  hM_gen[iCent][iY][iRun] = new TH1F(Form("hM_gen_%i_%i_%i", iCent, iY, iRun),"hM_gen", 400, 0, 10);
                  hM_gen[iCent][iY][iRun]->Sumw2();
                }
          }
    }
  

//  for(int iCent=0; iCent<CENTRALITY_CLASS_SELECTION; iCent++)

  
    for(int iRun = 0; iRun < NUMOFFILES_LHC15o_18q_18r; iRun++)
    {
        TString inputSimFile;
        TString inputSimFile2;
        TFile *inputFile;
        TFile *inputFile2;
        if(iRun < LHC15o_NUMOFFILES)
        {
          inputSimFile = Form("$HOME/local/anaedPbPbSimData2016/LHC16e2_jpsi_ttree/CINTEvent_MCinput_weight_AnaResults/PbPbMC_ResultsHistos_Run%i.root", LHC15o_run_number[iRun]);
          inputSimFile2 = Form("$HOME/local/anaedPbPbSimData2016/LHC16e2_plus_jpsi_ttree/CINTEvent_MCinput_weight_AnaResults/PbPbMC_ResultsHistos_Run%i.root", LHC15o_run_number[iRun]);

          inputFile = TFile::Open( inputSimFile.Data() );
          inputFile2 = TFile::Open( inputSimFile2.Data() );

        }
        else if(iRun >= LHC15o_NUMOFFILES)
        {
//          std::cout << LHC18q_r_run_number[iRun-LHC15o_NUMOFFILES] << std::endl;
          inputSimFile = Form("$HOME/local/anaedPbPbSim2018/LHC19a2_jpsi_ttree/CINTEvent_MCinput_weight_AnaResults/PbPbMC_ResultsHistos_Run%i.root",LHC18q_r_run_number[iRun-LHC15o_NUMOFFILES]);
          inputFile = TFile::Open( inputSimFile.Data() );        
        }
        else 
        {
          std::cout << " Please check the iRun value you use;" << std::endl;
          break;
        }
          std::cout << inputSimFile.Data() << std::endl; 
          
//        TFile *inputFile = TFile::Open( inputSimFile.Data() );
//       for(int iRun=0; iRun<NUMOFFILES; iRun++)
        for(int iCent=0; iCent<(MAX_CENTRALITY_CLASS_SELECTION-2); iCent++ )
        {

                  TString centrality_range = Form("_Cent%dto%d", arrayAvailableCentralityBins[iCent],arrayAvailableCentralityBins[iCent+1] );
                  std::cout << centrality_range.Data() << std::endl;
            
//                  TString inputSimFile_2 = Form("$HOME/local/anaedPbPbSim2018/LHC19a2/AnalysisResults_run%i.root",run_number[iRun]);


//                  TFile *inputFile_LHC16e2_plus = TFile::Open( inputSimFile_2.Data() );
 
//                  TList *tlist_rec[2]; tlist_rec[0] = new TList(); tlist_rec[1] = new TList();
//                  tlist_rec[0] = (TList *)inputFile->Get("PhysicsINT7inMUON");
//                  if(iRun < LHC15o_NUMOFFILES) tlist_rec[1] = (TList *)inputFile2->Get("PhysicsINT7inMUON");
//                  tlist_rec[1] = (TList *)inputFile_LHC16e2_plus->Get("PhysicsINT7inMUON");

//                  TList *tlist_gen[2]; tlist_gen[0] = new TList(); tlist_gen[1] = new TList();
//                  tlist_gen[0] = (TList *)inputFile->Get("GeneratedINT7inMUON");
//                  if(iRun < LHC15o_NUMOFFILES) tlist_gen[1] = (TList *)inputFile2->Get("GeneratedINT7inMUON");
//                  tlist_gen[1] = (TList *)inputFile_LHC16e2_plus->Get("GeneratedINT7inMUON");



                  hY_rec[iCent][iRun]->Add( (TH1F*)inputFile->Get( Form("DimuonHistosRecJpsi/Rapidity/histoRapidity%s", centrality_range.Data())) );
                  if(iRun < LHC15o_NUMOFFILES) hY_rec[iCent][iRun]->Add( (TH1F*)inputFile2->Get( Form("DimuonHistosRecJpsi/Rapidity/histoRapidity%s", centrality_range.Data() )) );
                  hY_rec_sum->Add(hY_rec[iCent][iRun]);

                  hY_gen[iCent][iRun]->Add( (TH1F*)inputFile->Get( Form("DimuonHistosGenJpsi/Rapidity/histoRapidity%s", centrality_range.Data())) );
                  if(iRun < LHC15o_NUMOFFILES) hY_gen[iCent][iRun]->Add( (TH1F*)inputFile2->Get( Form("DimuonHistosGenJpsi/Rapidity/histoRapidity%s", centrality_range.Data())) );

                  hY_gen_sum->Add(hY_gen[iCent][iRun]);

                  hY_ratio[iCent][iRun]->Divide(hY_rec[iCent][iRun], hY_gen[iCent][iRun]); 

                  for( int iY=0; iY<number_y_bin; iY++)
                  {
                    TString y_range = Form("_Y%gto%g", y_bin_pos[iY], y_bin_pos[iY+1] );
//                    std::cout << y_range.Data() << std::endl;

                    hM_rec[iCent][iY][iRun]->Add( (TH1F*)inputFile->Get( Form("DimuonHistosRecJpsi/InvMass/histoInvMass%s%s", centrality_range.Data(), y_range.Data() )) );
                    if(iRun < LHC15o_NUMOFFILES) hM_rec[iCent][iY][iRun]->Add( (TH1F*)inputFile2->Get( Form("DimuonHistosRecJpsi/InvMass/histoInvMass%s%s", centrality_range.Data(), y_range.Data() )) );
                    hM_gen[iCent][iY][iRun]->Add( (TH1F*)inputFile->Get( Form("DimuonHistosGenJpsi/InvMass/histoInvMass%s%s", centrality_range.Data(), y_range.Data() )) );
                    if(iRun < LHC15o_NUMOFFILES) hM_gen[iCent][iY][iRun]->Add( (TH1F*)inputFile2->Get( Form("DimuonHistosGenJpsi/InvMass/histoInvMass%s%s", centrality_range.Data(), y_range.Data() )) );
                  }
                  //hPt_ratio[iRun]->SetAxisRange(0.0,display_pTMaxBin,"X"); hPt_ratio[iRun]->SetAxisRange(0.0, 1.62, "Y");
//                  std::cout << hPt_rec[iCent][iRun] << 
//                  delete tlist_rec[0];
//                  delete tlist_rec[1];
//                  delete tlist_gen[0];
//                  delete tlist_gen[1]; 
        }
     inputFile->Close();
     if(iRun < LHC15o_NUMOFFILES) inputFile2->Close();
    }
  

  double entries_rec_hM = 0;
  double entries_gen_hM = 0;

  for(int iRun=0; iRun < NUMOFFILES_LHC15o_18q_18r; iRun++)
  {
//     for(int iCent=0; iCent<(MAX_CENTRALITY_CLASS_SELECTION-2); iCent++ )
     for(int iCent=0; iCent< 1; iCent++ )
     {
//       for( int iY=0; iY<; iY++)
       for( int iY=0; iY<1; iY++)
       {
//           entries_rec_hM += hM_rec[iCent][iY][iRun]->Integral(1, 400);
//           entries_gen_hM += hM_gen[iCent][iY][iRun]->Integral(1, 400);
           entries_rec_hM += hY_rec[iCent][iRun]->Integral(1, NBins);
           entries_gen_hM += hY_gen[iCent][iRun]->Integral(1, NBins);
//           entries_rec_hM += hM_rec[iCent][iY][iRun]->GetEntries();
//           entries_gen_hM += hM_gen[iCent][iY][iRun]->GetEntries();
       }
     }
  }

  std::cout <<"total"<< entries_rec_hM << ", " << entries_gen_hM << std::endl;


  TCanvas *c1 = new TCanvas("c1","",200, 10, 800, 600);
  c1->SetLogy();
  hY_gen_sum->SetMarkerColor(2);
  hY_gen_sum->SetLineColor(2);
  hY_gen_sum->Draw();
  hY_rec_sum->Draw("SAME");
  c1->Update();
 
//  double error;
//  std::cout << "Gen Njpsi: " << hY_gen_sum->Integral(1, 200) << " +/- " << error << std::endl;
//  error=0;
//  std::cout << "Rec Njpsi: " << hY_rec_sum->Integral(1, 200) << " +/- " << error << std::endl;

  TPaveStats *ps2 = (TPaveStats*)c1->GetPrimitive("stats");
  ps2->SetOptStat(1);
  ps2->SetName("Nrec_gen_stats");
  TList *listOfLines_rec_gen = ps2->GetListOfLines();

  TLatex *myt_rec_gen[2];
  myt_rec_gen[0] = new TLatex(0,0, Form("red: generated number of J/psi") );
  myt_rec_gen[0] ->SetTextFont(42);
  myt_rec_gen[0] ->SetTextSize(0.04);
  myt_rec_gen[0] ->SetTextColor(kRed);
  listOfLines_rec_gen->Add(myt_rec_gen[0]);

  myt_rec_gen[1] = new TLatex(0,0, Form("blue: reconstructed number of J/psi") );
  myt_rec_gen[1] ->SetTextFont(42);
  myt_rec_gen[1] ->SetTextSize(0.04);
  myt_rec_gen[1] ->SetTextColor(4);
  listOfLines_rec_gen->Add(myt_rec_gen[1]);

  hY_gen_sum->SetStats(0);
  c1->Modified();
  c1->Print("NumberOfJpsi_Cent0to20.pdf");
  delete c1;


  histo_Ae_run = new TH1F("Ae_vs_run","; Run number; J/#psi A*#epsilon", NUMOFFILES_LHC15o_18q_18r+1, 0, NUMOFFILES_LHC15o_18q_18r+1);
  histo_Ae_run->SetStats(0);
  TAxis *xaxis = histo_Ae_run->GetXaxis();
  


  double entries_rec_per_run[MAX_CENTRALITY_CLASS_SELECTION][NUMOFFILES_LHC15o_18q_18r]={0};
  double entries_gen_per_run[MAX_CENTRALITY_CLASS_SELECTION][NUMOFFILES_LHC15o_18q_18r]={0};

  double Ae_per_run_cent[MAX_CENTRALITY_CLASS_SELECTION][NUMOFFILES_LHC15o_18q_18r]={0};
  double errorAe_per_run_cent[MAX_CENTRALITY_CLASS_SELECTION][NUMOFFILES_LHC15o_18q_18r]={0};

  double Ae_per_run[NUMOFFILES_LHC15o_18q_18r]={0};
  double errorAe_per_run[NUMOFFILES_LHC15o_18q_18r]={0};
  double centrality_events[NUMOFFILES_LHC15o_18q_18r]={0};


//
//
// Calculate the AxEff in each centrality bin per run
//
//

  entries_rec_hM = 0;
  entries_gen_hM = 0;

//  for(int iCent=0; iCent<(MAX_CENTRALITY_CLASS_SELECTION-2); iCent++)
  for(int iCent=0; iCent<1; iCent++)
  {
      for(int iRun=0; iRun<NUMOFFILES_LHC15o_18q_18r; iRun++)
      {
//        entries_rec_per_run[iCent][iRun] += hY_rec[iCent][iRun]->Integral(1, NBins); 
//        entries_gen_per_run[iCent][iRun] += hY_gen[iCent][iRun]->Integral(1, NBins); 
//          entries_rec_per_run[iCent][iRun] += hPt_rec[iCent][iRun]->Integral(1, 120); 
//          entries_gen_per_run[iCent][iRun] += hPt_gen[iCent][iRun]->Integral(1, 120); 

//        entries_rec_hM += hY_rec[iCent][iRun]->Integral(hY_rec[iCent][iRun]->FindBin(-2.75), hY_rec[iCent][iRun]->FindBin(-2.5)); 
//        entries_gen_hM += hY_gen[iCent][iRun]->Integral(hY_gen[iCent][iRun]->FindBin(-2.75), hY_gen[iCent][iRun]->FindBin(-2.5));
         entries_rec_hM += hY_rec[iCent][iRun]->Integral(1, NBins);
         entries_gen_hM += hY_gen[iCent][iRun]->Integral(1, NBins); 
      }
  }
  
  std::cout << "total j/psi: " << entries_rec_hM << " " << entries_gen_hM << std::endl;

  for(int iCent=0; iCent<(MAX_CENTRALITY_CLASS_SELECTION-2); iCent++)
  {
      for(int iRun=0; iRun<NUMOFFILES_LHC15o_18q_18r; iRun++)
      {
         Ae_per_run_cent[iCent][iRun] = (entries_rec_per_run[iCent][iRun]/entries_gen_per_run[iCent][iRun]);
         errorAe_per_run_cent[iCent][iRun] = cal_binomial_error(entries_rec_per_run[iCent][iRun],entries_gen_per_run[iCent][iRun]);
      }
   }

//
// AxEff integrated over pT and centrality run by run 
// Sum the AxEff over the centrality bin with CMUL event weights per run 
//
//

  for(int iCent=0; iCent<(MAX_CENTRALITY_CLASS_SELECTION-2); iCent++)  //0-90%
  {
      for(int iRun=0; iRun<NUMOFFILES_LHC15o_18q_18r; iRun++)  
      {
//        if(iRun < LHC15o_NUMOFFILES) 
//        {
           Ae_per_run[iRun] += Ae_per_run_cent[iCent][iRun] * entries_rec_per_run[iCent][iRun];
           errorAe_per_run[iRun] += errorAe_per_run_cent[iCent][iRun] * errorAe_per_run_cent[iCent][iRun]* entries_rec_per_run[iCent][iRun]*entries_rec_per_run[iCent][iRun];
           centrality_events[iRun] += entries_rec_per_run[iCent][iRun];
//        }
//        else 
//        {
//           Ae_per_run[iRun] += Ae_per_run_cent[iCent][iRun] * numOfevents_Cent_per_Run[iRun-LHC15o_NUMOFFILES][iCent];
//           errorAe_per_run[iRun] += errorAe_per_run_cent[iCent][iRun] * errorAe_per_run_cent[iCent][iRun]* numOfevents_Cent_per_Run[iRun-LHC15o_NUMOFFILES][iCent]*numOfevents_Cent_per_Run[iRun-LHC15o_NUMOFFILES][iCent];
//           centrality_events[iRun] += numOfevents_Cent_per_Run[iRun-LHC15o_NUMOFFILES][iCent];
//        }
      }
  }

  for(int iRun=0; iRun<NUMOFFILES_LHC15o_18q_18r; iRun++)
  {
     Ae_per_run[iRun] = Ae_per_run[iRun]/centrality_events[iRun];
     errorAe_per_run[iRun] = sqrt(errorAe_per_run[iRun]/(centrality_events[iRun]*centrality_events[iRun]));
  }

  for(int iRun=0; iRun<NUMOFFILES_LHC15o_18q_18r; iRun++)
  {
     if(iRun < LHC15o_NUMOFFILES) xaxis->SetBinLabel(iRun+1, Form("%s", tstr_LHC15o_run_number[iRun].Data()));
     else xaxis->SetBinLabel(iRun+1, Form("%s",tstr_LHC18q_r_run_number[iRun-LHC15o_NUMOFFILES].Data()));

     histo_Ae_run->SetBinContent(iRun+1, Ae_per_run[iRun]);
     histo_Ae_run->SetBinError(iRun+1, errorAe_per_run[iRun]);
  }

//  histo_Ae_run->SetAxisRange(0.06,0.165,"Y");

  double numerator_total_Ae = 0;
  double denominator_total_events = 0;
  double numerator_total_errorAe =0;
  double avg_Ae = 0;
  double avg_errorAe = 0;

  for(int iRun=0; iRun<NUMOFFILES_LHC15o_18q_18r; iRun++) 
  {
    if(iRun < LHC15o_NUMOFFILES) 
    {
        numerator_total_Ae += Ae_per_run[iRun]*LHC15o_numOfevents_per_Run[iRun];
        denominator_total_events += (LHC15o_numOfevents_per_Run[iRun]);
        numerator_total_errorAe += errorAe_per_run[iRun]*errorAe_per_run[iRun]*(LHC15o_numOfevents_per_Run[iRun])*(LHC15o_numOfevents_per_Run[iRun]);
    }
    else 
    {
        numerator_total_Ae += Ae_per_run[iRun] *(numOfevents_per_Run[iRun-LHC15o_NUMOFFILES]); 
        denominator_total_events += (numOfevents_per_Run[iRun-LHC15o_NUMOFFILES]);
        numerator_total_errorAe += errorAe_per_run[iRun]*errorAe_per_run[iRun]*(numOfevents_per_Run[iRun-LHC15o_NUMOFFILES])*(numOfevents_per_Run[iRun-LHC15o_NUMOFFILES]);
    }
      

  }
  avg_Ae = numerator_total_Ae / denominator_total_events;
  avg_errorAe = sqrt(numerator_total_errorAe/(denominator_total_events*denominator_total_events)); 
  
  std::cout << "average Ae: " << avg_Ae << "+/-" << avg_errorAe <<std::endl;  

//  std::cout << "average Ae: " << numerator_total_Ae / denominator_total_events << "+/-" << sqrt(numerator_total_errorAe/(denominator_total_events*denominator_total_events)) <<std::endl;  

//
//  
//Calculate the A*E as a function of pT without weight
//
//


  for(int iPt=0; iPt<(number_y_bin); iPt++)
  {
    TString strYBins;
//             strPtBins.Form("Pt%gto%g",arrayPtBinsPhotonProduction[iPt],arrayPtBinsPhotonProduction[iPt+1]);
    strYBins.Form("Y%gto%g",y_bin[iPt], y_bin[iPt+1]);
//             std::cout << strPtBins.Data() <<std::endl;             
    const double display_MinBin = y_bin[iPt];
    const double display_MaxBin = y_bin[iPt+1];

    const int MinBin_statBox = (int)display_MinBin;
    const int MaxBin_statBox = (int)display_MaxBin;

    int lower_limit_bin_to_be_integral = ((display_MinBin - MinBin) / value_each_bin )+ 1;
    int upper_limit_bin_to_be_integral = ((display_MaxBin - MinBin) / value_each_bin) + 1;

    double error_phy_for_integral = 0;  double error_MC_for_integral = 0; double error_for_integral=0;
    double Ae =0; double error_Ae=0;
    double entries_hY_gen = hY_gen_sum->IntegralAndError(lower_limit_bin_to_be_integral,upper_limit_bin_to_be_integral, error_phy_for_integral);
    double entries_hY_rec = hY_rec_sum->IntegralAndError(lower_limit_bin_to_be_integral,upper_limit_bin_to_be_integral, error_phy_for_integral);

/*
    Ae = entries_hPt_rec / entries_hPt_gen;
    error_Ae = cal_binomial_error( entries_hPt_rec, entries_hPt_gen );

    hAe_Pt_wo_weight->SetBinContent((iPt+1), Ae);
    hAe_Pt_wo_weight->SetBinError((iPt+1), error_Ae );
*/
  }


//
//  
// Calculate the A*E as a function of pT with weight integrated over centrality and run
// Following tripo loops deal with the run weights
//
//

  double meta_Ae[NUMOFFILES_LHC15o_18q_18r][CENTRALITY_CLASS_SELECTION][numberOfPtBinsPhotonProduction] = {0};
  double meta_errorAe[NUMOFFILES_LHC15o_18q_18r][CENTRALITY_CLASS_SELECTION][numberOfPtBinsPhotonProduction]={0};
  double meta_numOfevents_per_Run[NUMOFFILES_LHC15o_18q_18r][CENTRALITY_CLASS_SELECTION][numberOfPtBinsPhotonProduction] ={0};
  double meta_numOfruns[CENTRALITY_CLASS_SELECTION][ELEMENT_FOR_AE_PT_BIN] = {0};

  double entries_hPt_rec_cent_pt[CENTRALITY_CLASS_SELECTION][ELEMENT_FOR_AE_PT_BIN] = {0};
  double entries_hPt_gen_cent_pt[CENTRALITY_CLASS_SELECTION][ELEMENT_FOR_AE_PT_BIN] = {0};
  double Ae_cent_pt[CENTRALITY_CLASS_SELECTION][ELEMENT_FOR_AE_PT_BIN] = {0};
  double errorAe_cent_pt[CENTRALITY_CLASS_SELECTION][ELEMENT_FOR_AE_PT_BIN] = {0};
  
  for(int iCent = 0; iCent < CENTRALITY_CLASS_SELECTION; iCent++)
  { 
    int iCounter = 0;

    for(int iRun = 0; iRun < NUMOFFILES_LHC15o_18q_18r; iRun++)
    {
//      std::cout <<"Run number: " << run_number[iRun] << " ";
//      std::cout << std::endl;
//      for(int iPt=0; iPt<7; iPt++)
      for(int iPt=0; iPt<(number_y_bin); iPt++)
      {
             TString strYBins;
//             strPtBins.Form("Pt%gto%g",arrayPtBinsPhotonProduction[iPt],arrayPtBinsPhotonProduction[iPt+1]);
             strYBins.Form("Y%gto%g", y_bin[iPt], y_bin[iPt+1]);
//             std::cout << strPtBins.Data() <<std::endl;             
             const double display_MinBin = y_bin[iPt];
             const double display_MaxBin = y_bin[iPt+1];

             const int MinBin_statBox = (int)display_MinBin;
             const int MaxBin_statBox = (int)display_MaxBin;

//             int lower_limit_bin_to_be_integral = ((display_MinBin - MinBin) / value_each_bin )+ 1;
//             int upper_limit_bin_to_be_integral = ((display_MaxBin - MinBin) / value_each_bin) + 1;
             int lower_limit_bin_to_be_integral = hY_rec[iCent][iRun]->FindBin(y_bin[iPt]);
             int upper_limit_bin_to_be_integral = hY_rec[iCent][iRun]->FindBin(y_bin[iPt+1]);
//             std::cout <<"bin number: " << lower_limit_bin_to_be_integral << " " << upper_limit_bin_to_be_integral << std::endl;         

             //
             //Compute the binomial error in each pT bin
             //
             for(int iBin=1; iBin<= NBins; iBin++)
             {
                    if(hY_ratio[iCent][iRun]->GetBinContent(iBin) ==0 ) hY_ratio[iCent][iRun]->SetBinError(iBin,0 );
                    else hY_ratio[iCent][iRun]->SetBinError(iBin, cal_binomial_error(hY_rec[iCent][iRun]->GetBinContent(iBin), hY_gen[iCent][iRun]->GetBinContent(iBin)));
                    //std::cout << hPt_ratio->GetBinContent(iBin) <<"+/-"<< hPt_ratio->GetBinError(iBin) << std::endl;
             }

             double error_phy_for_integral = 0;  double error_MC_for_integral = 0; double error_for_integral=0;
             double entries_hY_rec = hY_rec[iCent][iRun]->IntegralAndError(lower_limit_bin_to_be_integral, upper_limit_bin_to_be_integral-1, error_phy_for_integral);
             double entries_hY_gen = hY_gen[iCent][iRun]->IntegralAndError(lower_limit_bin_to_be_integral, upper_limit_bin_to_be_integral-1, error_MC_for_integral); 

             
             if(entries_hY_gen < entries_hY_rec)
             {                
                iCounter += 1;
                if(iRun < LHC15o_NUMOFFILES)
                {
                  hRun_rec_greater_gen[iCent]->GetXaxis()->SetBinLabel(iCounter, Form("%i", LHC15o_run_number[iRun] ));
                  hRun_rec_greater_gen[iCent]->SetBinContent(iCounter,1);
                  std::cout << "run number with #recjpsi > #genjpsi: " << LHC15o_run_number[iRun] << std::endl;
                }
                else 
                {
                  hRun_rec_greater_gen[iCent]->GetXaxis()->SetBinLabel(iCounter, Form("%i", LHC18q_r_run_number[iRun-LHC15o_NUMOFFILES] ));
                  hRun_rec_greater_gen[iCent]->SetBinContent(iCounter,1);
                  std::cout << "run number with #recjpsi > #genjpsi: " << LHC18q_r_run_number[iRun] << std::endl;
                }
             }

             entries_hPt_rec_cent_pt[iCent][iPt] += entries_hY_rec; // sum over runs
             entries_hPt_gen_cent_pt[iCent][iPt] += entries_hY_gen;

             double Ae =0; double error_Ae=0;
//             if(entries_hPt_rec != 0 )
//             {
                if(entries_hY_gen != 0)
                {                      
                       Ae = entries_hY_rec / entries_hY_gen;
//                       Ae = hPt_rec[iCent][iRun]->IntegralAndError(lower_limit_bin_to_be_integral,upper_limit_bin_to_be_integral, error_phy_for_integral) / hPt_gen[iCent][iRun]->IntegralAndError(lower_limit_bin_to_be_integral,upper_limit_bin_to_be_integral, error_MC_for_integral);
//                       error_phy_for_integral=0; error_MC_for_integral=0;
                       error_Ae = cal_binomial_error( entries_hY_rec, entries_hY_gen );

                       if(iRun < LHC15o_NUMOFFILES)          
                       {
                           meta_numOfevents_per_Run[iRun][iCent][iPt] = LHC15o_numOfevents_Cent_per_Run[iRun][iCent+1];
//                             meta_numOfevents_per_Run[iRun][iCent][iPt] = entries_rec_per_run[iCent][iRun];
                       }
                       else  meta_numOfevents_per_Run[iRun][iCent][iPt] = numOfevents_Cent_per_Run[iRun-LHC15o_NUMOFFILES][iCent+1];

                       meta_numOfruns[iCent][iPt] += 1;
//                       meta_numOfevents_per_Run[iRun][iCent][iPt] = numOfevents_per_Run[iRun];
                }
                else if(entries_hY_gen == 0 | entries_hY_rec == 0)
                {
                  std::cout << LHC18q_r_run_number[iRun] << " " << y_bin[iPt] << " " << y_bin[iPt+1];
                  std::cout << " " << entries_hY_rec << " " << entries_hY_gen << std::endl;

                }
//             }
/*
             std::cout <<"Ae: "<< Ae;
             std::cout << "+/-"<< error_Ae;
             std::cout <<" number of events: " << meta_numOfevents_per_Run[iRun][iCent][iPt];
             std::cout << std::endl;
*/
//             meta_Ae[iRun][iCent][iPt] = Ae *numOfevents_per_Run[iRun];
//             meta_errorAe[iRun][iCent][iPt] = error_Ae * error_Ae *numOfevents_per_Run[iRun];

//             Ae = hM_rec[iCent][5-iPt][iRun]->GetEntries() / hM_gen[iCent][5-iPt][iRun]->GetEntries();
//             error_Ae = cal_binomial_error(hM_rec[iCent][5-iPt][iRun]->GetEntries(), hM_gen[iCent][5-iPt][iRun]->GetEntries());

             meta_Ae[iRun][iCent][iPt] = Ae;
             meta_errorAe[iRun][iCent][iPt] = error_Ae * error_Ae ;

             
//             std::cout << Ae * numOfevents_per_Run[iRun];
//             std::cout << " meta Ae: " << meta_Ae[iRun][iCent][iPt];
//             std::cout << " meta error Ae: " << meta_errorAe[iRun][iCent][iPt];
//             std::cout << std::endl;
      } //for(int iPt=0; iPt<7; iPt++)
//      std::cout << std::endl;
          
    } //for(int iRun=0; iRun<NUMOFFILES; iRun++)
       
  } //(int iCent=0; iCent<CENTRALITY_CLASS_SELECTION; iCent++)



  for(int iCent=0; iCent<CENTRALITY_CLASS_SELECTION; iCent++) 
  {
    for(int iPt=0; iPt<(ELEMENT_FOR_AE_PT_BIN-1); iPt++) 
    {
//       std::cout << meta_numOfruns[iCent][iPt] - NUMOFFILES << std::endl;
       Ae_cent_pt[iCent][iPt] = entries_hPt_rec_cent_pt[iCent][iPt] / entries_hPt_gen_cent_pt[iCent][iPt]; 
       errorAe_cent_pt[iCent][iPt] = cal_binomial_error(entries_hPt_rec_cent_pt[iCent][iPt], entries_hPt_gen_cent_pt[iCent][iPt]);
//       std::cout << entries_hPt_rec_cent_pt[iCent][iPt] << " " << entries_hPt_gen_cent_pt[iCent][iPt] << std::endl;
    }
    
    std::cout << std::endl;
  }

    TCanvas *c_run = new TCanvas("c_run","",200, 10, 800, 600);
    c_run->Update();
    TString c_run_save_name[3] = {"run_number_rec_greater_gen.pdf(","run_number_rec_greater_gen.pdf","run_number_rec_greater_gen.pdf)"};
    for(int iCent=0; iCent <CENTRALITY_CLASS_SELECTION; iCent++)
    {
          hRun_rec_greater_gen[iCent]->Draw();
          if(iCent == 0) c_run->Print(c_run_save_name[0].Data());
          else if(iCent == (CENTRALITY_CLASS_SELECTION-1) ) c_run->Print(c_run_save_name[2].Data());
          else c_run->Print(c_run_save_name[1].Data());
    }
    

    double Ae_summation_numerator=0;
    double errorAe_summation_numerator=0;
    double total_effective_events=0;
//    double mean_Ae_run_by_run[CENTRALITY_CLASS_SELECTION] = {0};
//    double mean_errorAe_run_by_run[CENTRALITY_CLASS_SELECTION] = {0};

    double mean_Ae_run_by_run[numberOfPtBinsPhotonProduction]= {0};
    double mean_errorAe_run_by_run[numberOfPtBinsPhotonProduction] = {0};

    double run_number_weight[NUMOFFILES_LHC15o_18q_18r][numberOfPtBinsPhotonProduction]={0};
    
    for(int iPt=0; iPt<(ELEMENT_FOR_AE_PT_BIN-1); iPt++)
    {
        for(int iRun=0; iRun<NUMOFFILES_LHC15o_18q_18r; iRun++) for(int iCent=0; iCent<CENTRALITY_CLASS_SELECTION; iCent++)  run_number_weight[iRun][iPt] += meta_numOfevents_per_Run[iRun][iCent][iPt];
    }
   

    double mean_Ae_add_cent_weight=0;
    double mean_errorAe_cent_weight=0;
    double total_effective_centrality_events = 0;

    double mean_Ae_add_cent_weight2=0;
    double mean_errorAe_cent_weight2=0;

  //weight for each run's efficiency with event numbers in various pT regions.
//  for(int iPt=0; iPt<7; iPt++)
//  for(int iCent=0; iCent < 1; iCent++)
//  for(int iCent = 0; iCent < (CENTRALITY_TUNING); iCent++)
  for(int iCent = 0; iCent < 2; iCent++)
  {

    std::cout << " Cent:" << iCent << std::endl;
    for(int iPt=0; iPt<(number_y_bin); iPt++)
//  for(int iPt=(ELEMENT_FOR_AE_PT_BIN-2); iPt<(ELEMENT_FOR_AE_PT_BIN-1); iPt++)
    {
          TString strYBins;
//        strPtBins.Form("Pt%gto%g",arrayPtBinsPhotonProduction[iPt],arrayPtBinsPhotonProduction[iPt+1]);
          strYBins.Form("Y%gto%g",y_bin[iPt], y_bin[iPt+1]);
//        for(int iCent=0; iCent<CENTRALITY_CLASS_SELECTION; iCent++)
        //  for(int iCent=0; iCent<1; iCent++)
//        {
          Ae_summation_numerator=0;
          errorAe_summation_numerator=0;
          total_effective_events=0;
          for(int iRun=0; iRun<NUMOFFILES_LHC15o_18q_18r; iRun++)
          {
              // Weights for run by run

//                   if(meta_Ae[iRun][iCent][iPt] == 1) continue;
                      if( meta_Ae[iRun][iCent][iPt] == 2 ) continue;
                      else
                      {
//                             Ae_summation_numerator[iCent] += meta_Ae[iRun][iCent][iPt]*meta_numOfevents_per_Run[iRun][iCent];
//                             total_effective_events[iCent] += meta_numOfevents_per_Run[iRun][iCent]; //for run number weighting
//                               Ae_summation_numerator += meta_Ae[iRun][iCent][iPt]*meta_numOfevents_per_Run[iRun][iCent];
//                               Ae_summation_numerator += meta_Ae[iRun][iCent][iPt]*numOfevents_per_Run[iRun];

//                               Ae_summation_numerator += meta_Ae[iRun][iCent][iPt] *numOfevents_per_Run[iRun];
//                               errorAe_summation_numerator += meta_errorAe[iRun][iCent][iPt] * numOfevents_per_Run[iRun];
//                               total_effective_events += meta_numOfevents_per_Run[iRun][iCent][iPt]; //for run number weighting

                               Ae_summation_numerator += meta_Ae[iRun][iCent][iPt] * meta_numOfevents_per_Run[iRun][iCent][iPt];
                               errorAe_summation_numerator += meta_errorAe[iRun][iCent][iPt] *meta_numOfevents_per_Run[iRun][iCent][iPt] *meta_numOfevents_per_Run[iRun][iCent][iPt];
                               total_effective_events += meta_numOfevents_per_Run[iRun][iCent][iPt]; //for run number weighting

//                               Ae_summation_numerator += meta_Ae[iRun][iCent][iPt] * run_number_weight[iRun][iPt];
//                               errorAe_summation_numerator += meta_errorAe[iRun][iCent][iPt] * run_number_weight[iRun][iPt];
//                               total_effective_events += run_number_weight[iRun][iPt];

                      }                   
       
             }
//             std::cout << Ae_summation_numerator << " " << errorAe_summation_numerator;
//             std::cout << std::endl;
             //           std::cout << Ae_summation_numerator[iCent] << " ";
             //           std::cout << total_effective_events[iCent] << std::endl;
             //           std::cout <<"Ae: " <<  Ae_summation_numerator[iCent] / total_effective_events[iCent] << std::endl;
//             Ae_summation_numerator[iCent] = Ae_summation_numerator[iCent] / total_effective_events[iCent];           
            mean_Ae_run_by_run[iPt] =  Ae_summation_numerator / total_effective_events; //mean_Ae is summed over runs in a give pT and centrality bin.
            mean_errorAe_run_by_run[iPt] = sqrt(errorAe_summation_numerator / (total_effective_events*total_effective_events) );
//             std::cout <<  Ae_summation_centrality << " ";
//             std::cout <<  total_effective_events[iCent] << std::endl;
//        }

        mean_Ae_add_cent_weight = 0;
        mean_errorAe_cent_weight =0;
        total_effective_centrality_events=0;

        mean_Ae_add_cent_weight2=0;
        mean_errorAe_cent_weight2 =0;

/*       
        for(int iCent=0; iCent<CENTRALITY_CLASS_SELECTION; iCent++) 
        {
              mean_Ae_add_cent_weight  += mean_Ae_run_by_run[iCent]*numOfevents_CentBins[iCent];
//              mean_errorAe_cent_weight += mean_errorAe_run_by_run[iCent]*numOfevents_CentBins[iCent];
              mean_errorAe_cent_weight += mean_errorAe_run_by_run[iCent] * mean_errorAe_run_by_run[iCent] * numOfevents_CentBins[iCent] * numOfevents_CentBins[iCent];
              total_effective_centrality_events += numOfevents_CentBins[iCent];

              mean_Ae_add_cent_weight2 += Ae_cent_pt[iCent][iPt]*numOfevents_CentBins[iCent];
              mean_errorAe_cent_weight2 += errorAe_cent_pt[iCent][iPt]*errorAe_cent_pt[iCent][iPt]*numOfevents_CentBins[iCent] * numOfevents_CentBins[iCent];
        }
*/

        std::cout <<" y region: " << strYBins.Data() << " ";
        std::cout <<"Ae: " << std::fixed << std::setprecision(6) <<  mean_Ae_run_by_run[iPt] <<" +/- "<< std::fixed << std::setprecision(6) << mean_errorAe_run_by_run[iPt] << std::endl;

//        hAe_Pt->SetBinContent((iPt+1), mean_Ae_run_by_run[iPt]);
//        hAe_Pt->SetBinError((iPt+1), mean_errorAe_run_by_run[iPt] ));

        hAe_cent_dep[iCent]->SetBinContent((iPt+1), mean_Ae_run_by_run[iPt]);
        hAe_cent_dep[iCent]->SetBinError((iPt+1), mean_errorAe_run_by_run[iPt] );

        hAe_Pt_wo_weight->SetBinContent((iPt+1), mean_Ae_add_cent_weight2 / total_effective_centrality_events);
        hAe_Pt_wo_weight->SetBinError((iPt+1), sqrt(mean_errorAe_cent_weight2 / (total_effective_centrality_events*total_effective_centrality_events)) );

//        std::cout << "Ae wo weight: " << hAe_Pt_wo_weight->GetBinContent(iPt+1);
        double relative_diff = 0;
//        relative_diff = (hAe_Pt->GetBinContent(iPt+1) - hAe_Pt_wo_weight->GetBinContent(iPt+1) ) / hAe_Pt_wo_weight->GetBinContent(iPt+1) ; 
//        std::cout << " " << relative_diff << std::endl;
    }
  }

//  std::cout << "where the bug is" << std::endl;
  TFile *outputFile = new TFile(Form("./AxEff_LHC15o_18q_18r.root"),"recreate");
  for(int iCent = 0; iCent < (CENTRALITY_TUNING); iCent++)
  {
    outputFile->Add(hAe_cent_dep[iCent]);
  }
  outputFile->Add(histo_Ae_run);

  outputFile->Write();

  TPaveText* t2 = new TPaveText(0.35, 0.15, 0.85,0.30,"NDC");
  t2->AddText(0.,0., "LHC15o, LHC18q, LHC18r" );
  t2->AddText(0.,0., "Centrality 0 - 90%" );
  t2->AddText(0.,0., "0 < p_{T} < 20 GeV/c" );
  t2->AddText(0.,0., "-4.0 < y < -2.5" );
  t2->AddText(0.,0., Form("#bar{A#epsilon}: %f #pm %f", avg_Ae, avg_errorAe ));
  t2->SetTextSize(0.02);

  TLine *line_average = new TLine();


  TCanvas *c2 = new TCanvas("c2","",200, 10, 800, 600);
//  hAe_Pt->SetMarkerColor(2);
//  hAe_Pt->SetLineColor(2);
  hAe_cent_dep[0]->Draw();
  c2->Update();

  TPaveStats *ps = (TPaveStats*)c2->GetPrimitive("stats");
  ps->SetY1NDC(0.7);
  ps->SetY2NDC(0.9);
  ps->SetOptStat(1);
  ps->SetName("mystats");
  TList *listOfLines = ps->GetListOfLines();

  TLatex *myt[3];
  myt[0] = new TLatex(0,0, Form("include run number weighting") );
  myt[0] ->SetTextFont(42);
  myt[0] ->SetTextSize(0.04);
  myt[0] ->SetTextColor(kRed);
//  listOfLines->Add(myt[0]);

  myt[1] = new TLatex(0,0, Form("include centrality weighting") );
  myt[1] ->SetTextFont(42);
  myt[1] ->SetTextSize(0.04);
  myt[1] ->SetTextColor(kRed);
//  listOfLines->Add(myt[1]);

  myt[2] = new TLatex(0,0, Form("include centrality weighting") );
  myt[2] ->SetTextFont(42);
  myt[2] ->SetTextSize(0.04);
  myt[2] ->SetTextColor(4);
//  listOfLines->Add(myt[2]);

  hAe_cent_dep[0]->SetStats(0);
  c2->Modified();

  hAe_Pt_wo_weight->Draw("same");
  hAe_Pt_wo_weight->SetStats(0);
  c2->Print("Ae_func_pT_or_run.pdf(");
  
  c2->cd();
  histo_Ae_run->Draw();
  line_average->DrawLine(0, avg_Ae, (NUMOFFILES_LHC15o_18q_18r), avg_Ae);
  t2->Draw("SAME");
  c2->Print("Ae_func_pT_or_run.pdf)");

  delete outputFile;
return 0;
}


//input: k and N. N is the sample size. k is the expectation value for the number of events.
//in this case: k is the number of reconstructed jpsi and N is the number of generated Jpsi
double cal_binomial_error(double k, double N)
{
 double error_efficiency = 0;

 if(  (1-(k/N) )<0 )
 {
    error_efficiency = 0;
 }
 else
 {
    error_efficiency = (1./N)*sqrt(k*(1-(k/N)) );
 }

return error_efficiency;
}
