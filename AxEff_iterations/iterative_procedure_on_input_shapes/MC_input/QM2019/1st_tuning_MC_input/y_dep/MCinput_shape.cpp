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


int main()
{

        const int arrayAvailableCentralityBins[] = {0,10,20,30,40,50,60,70,80,90,100};
        const int CENTRALITY_CLASS_SELECTION = 2; // for 0-10% and 10-20%
        const int MAX_CENTRALITY_CLASS_SELECTION = (int)(sizeof(arrayAvailableCentralityBins)/sizeof(arrayAvailableCentralityBins[0]));

//        const int ELEMENT_FOR_AE_PT_BIN =16;
//        float Pt_bin_for_Ae[ELEMENT_FOR_AE_PT_BIN]= {0,0.3,1,2,3,4,5,6,7,8,9,10,11,12,15,20};

        const int number_y_bin = 6;
        float y_bin[] = {-4., -3.75, -3.5, -3.25, -3., -2.75, -2.5};
        float y_bin_pos[] = {2.5, 2.75, 3, 3.25, 3.5, 3.75, 4};

        int invMassNBins = 400;
        double invMassMinBin = 0;
        double invMassMaxBin =10;

        int NBins = 60;
        double MinBin = -5;
        double MaxBin = -2;

        double value_each_bin = (MaxBin - MinBin) / NBins;

        TString inputSimFile;
        TString inputSimFile2;
        TString inputSimFile3;
        TString inputSimFile4;

        TFile *inputFile;
        TFile *inputFile2;
        TFile *inputFile3;
        TFile *inputFile4;

//        inputSimFile = Form("$HOME/local/anaedPbPbSimData2016/LHC16e2_with_centrality_class/AnalysisResults_MCembedding_LHC16e2.root");
//        inputSimFile2 = Form("$HOME/local/anaedPbPbSimData2016/LHC16e2_plus_with_centrality_class/AnalysisResults_MCembedding_LHC16e2_plus.root");

        inputSimFile = Form("$HOME/local/anaedPbPbSimData2016/LHC16e2_jpsi_ttree/CINTEvent_MCinput_weight_AnaResults/PbPbMC_ResultsHistos_LHC16e2.root");
        inputSimFile2 = Form("$HOME/local/anaedPbPbSimData2016/LHC16e2_plus_jpsi_ttree/CINTEvent_MCinput_weight_AnaResults/PbPbMC_ResultsHistos_LHC16e2_plus.root");
        inputFile = TFile::Open( inputSimFile.Data() );
        inputFile2 = TFile::Open( inputSimFile2.Data() );

//        inputSimFile3 = Form("$HOME/local/anaedPbPbSim2018/LHC19a2/AnalysisResults_MCembedding_LHC19a2.root");
        inputSimFile3 = ("$HOME/local/anaedPbPbSim2018/LHC19a2_jpsi_ttree/CINTEvent_MCinput_weight_AnaResults/PbPbMC_ResultsHistos_LHC18q.root");
        inputFile3 = TFile::Open( inputSimFile3.Data() );

        inputSimFile4 = ("$HOME/local/anaedPbPbSim2018/LHC19a2_jpsi_ttree/CINTEvent_MCinput_weight_AnaResults/PbPbMC_ResultsHistos_LHC18r.root");
        inputFile4 = TFile::Open( inputSimFile4.Data() );

        TH1F *h_gen[MAX_CENTRALITY_CLASS_SELECTION];
        TH1F *hM_gen[MAX_CENTRALITY_CLASS_SELECTION][number_y_bin];


        for(int iCent=0; iCent< MAX_CENTRALITY_CLASS_SELECTION-1; iCent++)
        {                                  
                  h_gen[iCent] = new TH1F(Form("hRapidityUnlikeSign_Cent%gto%g", arrayAvailableCentralityBins[iCent], arrayAvailableCentralityBins[iCent+1]),"", NBins, MinBin, MaxBin);
//                  h_gen[iCent] = new TH1F(Form("hRapidityUnlikeSign_Cent%gto%g", arrayAvailableCentralityBins[iCent], arrayAvailableCentralityBins[iCent+1]),"", number_y_bin, y_bin);

                  for( int iY=0; iY<number_y_bin; iY++)
                  {                                
                    hM_gen[iCent][iY] = new TH1F(Form("hM_gen_%i_%i", iCent, iY),"hM_gen", invMassNBins, invMassMinBin, invMassMaxBin);
                    hM_gen[iCent][iY]->Sumw2();                
                  }

                 
        }
        h_gen[MAX_CENTRALITY_CLASS_SELECTION-1] = new TH1F(Form("hRapidityUnlikeSign_Cent60to90"),"", NBins, MinBin, MaxBin);

        for(int iCent=0; iCent<(MAX_CENTRALITY_CLASS_SELECTION-2); iCent++ )
        {

                  TString centrality_range = Form("_Cent%dto%d", arrayAvailableCentralityBins[iCent],arrayAvailableCentralityBins[iCent+1] );
                  std::cout << centrality_range.Data() << std::endl;

//                  TString inputSimFile_2 = Form("$HOME/local/anaedPbPbSim2018/LHC19a2/AnalysisResults_run%i.root",run_number[iRun]);
//                  TFile *inputFile_LHC16e2_plus = TFile::Open( inputSimFile_2.Data() );

                  h_gen[iCent]->Add( (TH1F*)inputFile->Get( Form("DimuonHistosGenJpsi/Rapidity/histoRapidity%s", centrality_range.Data())) );
                  h_gen[iCent]->Add( (TH1F*)inputFile2->Get( Form("DimuonHistosGenJpsi/Rapidity/histoRapidity%s", centrality_range.Data())) );
                  h_gen[iCent]->Add( (TH1F*)inputFile3->Get( Form("DimuonHistosGenJpsi/Rapidity/histoRapidity%s", centrality_range.Data())) );
                  h_gen[iCent]->Add( (TH1F*)inputFile4->Get( Form("DimuonHistosGenJpsi/Rapidity/histoRapidity%s", centrality_range.Data())) );

//                  h_gen[iCent]->Add( (TH1F*)inputFile->Get( Form("DimuonHistosGenJpsi/GenJpsi_Y%s", centrality_range.Data())) );
//                  h_gen[iCent]->Add( (TH1F*)inputFile2->Get( Form("DimuonHistosGenJpsi/GenJpsi_Y%s", centrality_range.Data())) );
//                  h_gen[iCent]->Add( (TH1F*)inputFile3->Get( Form("DimuonHistosGenJpsi/GenJpsi_Y%s", centrality_range.Data())) );
//                  h_gen[iCent]->Add( (TH1F*)inputFile4->Get( Form("DimuonHistosGenJpsi/GenJpsi_Y%s", centrality_range.Data())) );

                  if(iCent >= 6)
                  {
                    h_gen[MAX_CENTRALITY_CLASS_SELECTION-1]->Add(h_gen[iCent]);
                  }

                 for( int iY=0; iY<number_y_bin; iY++)
                 {
                   TString y_range = Form("_Y%gto%g", y_bin_pos[iY], y_bin_pos[iY+1] );
                   hM_gen[iCent][iY]->Add( (TH1F*)inputFile->Get( Form("DimuonHistosGenJpsi/InvMass/histoInvMass%s%s", centrality_range.Data(), y_range.Data() )) );
                   hM_gen[iCent][iY]->Add( (TH1F*)inputFile2->Get( Form("DimuonHistosGenJpsi/InvMass/histoInvMass%s%s", centrality_range.Data(), y_range.Data() )) );
                   hM_gen[iCent][iY]->Add( (TH1F*)inputFile3->Get( Form("DimuonHistosGenJpsi/InvMass/histoInvMass%s%s", centrality_range.Data(), y_range.Data() )) );
                   hM_gen[iCent][iY]->Add( (TH1F*)inputFile4->Get( Form("DimuonHistosGenJpsi/InvMass/histoInvMass%s%s", centrality_range.Data(), y_range.Data() )) );                   
                 }
//                  hPt_ratio[iCent][iRun]->Divide(hPt_rec[iCent][iRun], hPt_gen[iCent][iRun]);
                  //hPt_ratio[iRun]->SetAxisRange(0.0,display_pTMaxBin,"X"); hPt_ratio[iRun]->SetAxisRange(0.0, 1.62, "Y");
//                  std::cout << hPt_rec[iCent][iRun] << 
        }


        TFile *outputFile = new TFile(Form("JpsiMCinputHistos.root"),"recreate");
        TDirectory *directoryJPsi;
        directoryJPsi = outputFile->mkdir( Form("JPsiMCinputHistos") );

        for(int iCent = 0; iCent < MAX_CENTRALITY_CLASS_SELECTION-2; iCent++)
        {
            directoryJPsi->cd();
            TH1F *jpsi_MCinput = new TH1F(Form("histoMCinput_cent%ito%i", arrayAvailableCentralityBins[iCent], arrayAvailableCentralityBins[iCent+1] ),"", number_y_bin, y_bin);
            jpsi_MCinput->Sumw2();
            
//            std::cout << arrayAvailableCentralityBins[iCent] << std::endl;
        }

        directoryJPsi->cd();
        TH1F *jpsi_MCinput = new TH1F(Form("histoMCinput_cent60to90"),"", number_y_bin, y_bin);
        jpsi_MCinput->Sumw2();
        
        for(int iCent=0; iCent< MAX_CENTRALITY_CLASS_SELECTION-2; iCent++)
        {
                   for(int iPt=0; iPt<(number_y_bin); iPt++)
                   {
                     TString strYBins;
//              strPtBins.Form("Pt%gto%g",arrayPtBinsPhotonProduction[iPt],arrayPtBinsPhotonProduction[iPt+1]);
                     strYBins.Form("Y%gto%g", y_bin[iPt], y_bin[iPt+1]);
                     std::cout << strYBins.Data() <<std::endl;             
                     const double display_MinBin = y_bin[iPt];
                     const double display_MaxBin = y_bin[iPt+1];

                     const int pTMinBin_statBox = (int)display_MinBin;
                     const int pTMaxBin_statBox = (int)display_MaxBin;

//                     int lower_limit_bin_to_be_integral = ((display_MinBin - MinBin) / value_each_bin )+ 0;
//                     int upper_limit_bin_to_be_integral = ((display_MaxBin - MinBin) / value_each_bin) + 0;

                     int lower_limit_bin_to_be_integral = h_gen[iCent]->FindBin( y_bin[iPt]);
                     int upper_limit_bin_to_be_integral = h_gen[iCent]->FindBin( y_bin[iPt+1]);

                     std::cout << lower_limit_bin_to_be_integral << " " << upper_limit_bin_to_be_integral << std::endl;
//                     std::cout << h_gen[iCent]->GetBinContent(h_gen[iCent]->FindBin( y_bin[iPt+1])) << std::endl;
                     double error_phy_for_integral = 0;  double error_MC_for_integral = 0; double error_for_integral=0;
                     
                     double entries_hPt_gen = h_gen[iCent]->IntegralAndError(lower_limit_bin_to_be_integral, upper_limit_bin_to_be_integral-1, error_MC_for_integral);
//                     double entries_hPt_gen = h_gen[iCent]->GetBinContent(iPt+1);
//                     error_MC_for_integral = h_gen[iCent]->GetBinError(iPt+1);

//                     double entries_hPt_gen = hM_gen[iCent][5-iPt]->GetEntries();
//                     error_MC_for_integral = sqrt(entries_hPt_gen);
//                     entries_hPt_gen = entries_hPt_gen / (y_bin[iPt] - y_bin[iPt+1]);
//                     std::cout << fabs(y_bin[iPt] - y_bin[iPt+1]) << std::endl;
                     std::cout << entries_hPt_gen << std::endl;
//                   double entries_hPt_rec = hPt_rec->IntegralAndError(lower_limit_bin_to_be_integral,upper_limit_bin_to_be_integral, error_phy_for_integral);
                     ((TH1F*) outputFile->Get(Form("JPsiMCinputHistos/histoMCinput_cent%ito%i", arrayAvailableCentralityBins[iCent], arrayAvailableCentralityBins[iCent+1]) ))->SetBinContent(iPt+1, entries_hPt_gen );
                     ((TH1F*) outputFile->Get(Form("JPsiMCinputHistos/histoMCinput_cent%ito%i", arrayAvailableCentralityBins[iCent], arrayAvailableCentralityBins[iCent+1]) ))->SetBinError(iPt+1, error_MC_for_integral );
                     
/*
    Ae = entries_hPt_rec / entries_hPt_gen;
    error_Ae = cal_binomial_error( entries_hPt_rec, entries_hPt_gen );

    hAe_Pt_wo_weight->SetBinContent((iPt+1), Ae);
    hAe_Pt_wo_weight->SetBinError((iPt+1), error_Ae );
*/
                      
                    }
         }
/*
         for(int iPt=0; iPt<( number_y_bin); iPt++)
         {                                         
                     const double display_MinBin = y_bin[iPt];
                     const double display_MaxBin = y_bin[iPt+1];

                     const int pTMinBin_statBox = (int)display_MinBin;
                     const int pTMaxBin_statBox = (int)display_MaxBin;

                     int lower_limit_bin_to_be_integral = h_gen[MAX_CENTRALITY_CLASS_SELECTION-1]->FindBin(y_bin[iPt]);
                     int upper_limit_bin_to_be_integral = h_gen[MAX_CENTRALITY_CLASS_SELECTION-1]->FindBin(y_bin[iPt+1]);

                     double error_phy_for_integral = 0;  double error_MC_for_integral = 0; double error_for_integral=0;
//                     double entries_gen = h_gen[MAX_CENTRALITY_CLASS_SELECTION-1]->IntegralAndError(lower_limit_bin_to_be_integral, upper_limit_bin_to_be_integral, error_MC_for_integral);
                     double entries_gen = hM_gen[MAX_CENTRALITY_CLASS_SELECTION-1][5-iPt]->GetEntries();
//                     error_MC_for_integral = sqrt(entries_gen);
//                     entries_gen = entries_gen / fabs(y_bin[iPt] - y_bin[iPt+1]);
                     ((TH1F*) outputFile->Get(Form("JPsiMCinputHistos/histoMCinput_cent60to90" )))->SetBinContent(iPt+1, entries_gen );
                     ((TH1F*) outputFile->Get(Form("JPsiMCinputHistos/histoMCinput_cent60to90" )))->SetBinError(iPt+1, error_MC_for_integral);
         }
*/
         outputFile->Write();

         inputFile->Close();
         inputFile2->Close();
         inputFile3->Close();
         inputFile4->Close();

         delete outputFile;

return 0;
}
