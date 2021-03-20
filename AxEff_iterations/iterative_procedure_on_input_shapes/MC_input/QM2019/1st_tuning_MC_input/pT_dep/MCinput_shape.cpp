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

        const int ELEMENT_FOR_AE_PT_BIN =16;
        float Pt_bin_for_Ae[ELEMENT_FOR_AE_PT_BIN]= {0,0.3,1,2,3,4,5,6,7,8,9,10,11,12,15,20};

        int invMassNBins = 400;
        double invMassMinBin = 0;
        double invMassMaxBin =10;

        int pTNBins = 2000;
        double pTMinBin = 0;
        double pTMaxBin =20;

        double value_each_pTbin = (pTMaxBin - pTMinBin) / pTNBins;

        TString inputSimFile;
        TString inputSimFile2;
        TString inputSimFile3;
        TString inputSimFile4;

        TFile *inputFile;
        TFile *inputFile2;
        TFile *inputFile3;
        TFile *inputFile4;

        inputSimFile = Form("$HOME/local/anaedPbPbSimData2016/LHC16e2_jpsi_ttree/CINTEvent_MCinput_weight_AnaResults/PbPbMC_ResultsHistos_LHC16e2.root");
        inputSimFile2 = Form("$HOME/local/anaedPbPbSimData2016/LHC16e2_plus_jpsi_ttree/CINTEvent_MCinput_weight_AnaResults/PbPbMC_ResultsHistos_LHC16e2_plus.root");
        inputFile = TFile::Open( inputSimFile.Data() );
        inputFile2 = TFile::Open( inputSimFile2.Data() );

        inputSimFile3 = Form("/home/chunlu/local/anaedPbPbSim2018/LHC19a2_jpsi_ttree/CINTEvent_MCinput_weight_AnaResults/PbPbMC_ResultsHistos_LHC18q.root");
        inputFile3 = TFile::Open( inputSimFile3.Data() );
 
        inputSimFile4 = Form("$HOME/local/anaedPbPbSim2018/LHC19a2_jpsi_ttree/CINTEvent_MCinput_weight_AnaResults/PbPbMC_ResultsHistos_LHC18r.root");
        inputFile4 = TFile::Open( inputSimFile4.Data() );

        TH1F *hpT_gen[MAX_CENTRALITY_CLASS_SELECTION];
        TH1F *hM_gen[MAX_CENTRALITY_CLASS_SELECTION][ELEMENT_FOR_AE_PT_BIN];

        for(int iCent=0; iCent< MAX_CENTRALITY_CLASS_SELECTION-2; iCent++)
        {                                  
               hpT_gen[iCent] = new TH1F(Form("hRapidityUnlikeSign_Cent%gto%g", arrayAvailableCentralityBins[iCent], arrayAvailableCentralityBins[iCent+1]),"", pTNBins, pTMinBin, pTMaxBin);
               hpT_gen[iCent]->Sumw2();
        }


        for(int iCent=0; iCent< MAX_CENTRALITY_CLASS_SELECTION-1; iCent++)
        {
                  for( int iPt=0; iPt< ELEMENT_FOR_AE_PT_BIN-1; iPt++)
                  {
                      hM_gen[iCent][iPt] = new TH1F(Form("hMUnlikeSign_Cent%dto%d_Pt%gto%g", arrayAvailableCentralityBins[iCent], arrayAvailableCentralityBins[iCent+1], Pt_bin_for_Ae[iPt], Pt_bin_for_Ae[iPt+1]),"", invMassNBins, invMassMinBin, invMassMaxBin);
                      hM_gen[iCent][iPt]->Sumw2();
                  }

        }



        for(int iCent=0; iCent<(MAX_CENTRALITY_CLASS_SELECTION-2); iCent++ )
        {

                  TString centrality_range = Form("_Cent%dto%d", arrayAvailableCentralityBins[iCent],arrayAvailableCentralityBins[iCent+1] );
                  std::cout << centrality_range.Data() << std::endl;


                  for(int ipT=0; ipT< (ELEMENT_FOR_AE_PT_BIN-1); ipT++)
                  {
                      TString pTrange = Form("_Pt%gto%g", Pt_bin_for_Ae[ipT], Pt_bin_for_Ae[ipT+1]);
                      std::cout <<  pTrange.Data() << std::endl;
                      hpT_gen[iCent]->Add( (TH1F*)inputFile->Get( Form("DimuonHistosGenJpsi/Pt/histoPt%s%s", centrality_range.Data(), pTrange.Data()  )) );
                      hpT_gen[iCent]->Add( (TH1F*)inputFile2->Get( Form("DimuonHistosGenJpsi/Pt/histoPt%s%s", centrality_range.Data(), pTrange.Data()  )) );
                      hpT_gen[iCent]->Add( (TH1F*)inputFile3->Get( Form("DimuonHistosGenJpsi/Pt/histoPt%s%s", centrality_range.Data(), pTrange.Data()  )) );
                      hpT_gen[iCent]->Add( (TH1F*)inputFile4->Get( Form("DimuonHistosGenJpsi/Pt/histoPt%s%s", centrality_range.Data(), pTrange.Data()  )) );

                      hM_gen[iCent][ipT]->Add( (TH1F*)inputFile->Get( Form("DimuonHistosGenJpsi/InvMass/histoInvMass%s%s", centrality_range.Data(), pTrange.Data()  )) );
                      hM_gen[iCent][ipT]->Add( (TH1F*)inputFile2->Get( Form("DimuonHistosGenJpsi/InvMass/histoInvMass%s%s", centrality_range.Data(), pTrange.Data()  )) );
                      hM_gen[iCent][ipT]->Add( (TH1F*)inputFile3->Get( Form("DimuonHistosGenJpsi/InvMass/histoInvMass%s%s", centrality_range.Data(), pTrange.Data()  )) );
                      hM_gen[iCent][ipT]->Add( (TH1F*)inputFile4->Get( Form("DimuonHistosGenJpsi/InvMass/histoInvMass%s%s", centrality_range.Data(), pTrange.Data()  )) );                    
                  }

        }

        std::cout << "pT integral " << hpT_gen[0]->Integral(1, (pTNBins-1)) << " " << std::endl;
        double var=0;
        double var2=0;
        for(int ipT=0; ipT< (ELEMENT_FOR_AE_PT_BIN-1); ipT++)
        {          
            var += hM_gen[0][ipT]->Integral(1,invMassNBins);
            var2 += hM_gen[0][ipT]->GetEntries();
        }
        std::cout << "invmass integral" << var  << ", " << var2 << std::endl;

        TFile *outputFile = new TFile(Form("JpsiMCinputHistos.root"),"recreate");
        TDirectory *directoryJPsi;
        directoryJPsi = outputFile->mkdir( Form("JPsiMCinputHistos") );

        for(int iCent = 0; iCent < MAX_CENTRALITY_CLASS_SELECTION-2; iCent++)
        {
            directoryJPsi->cd();
            TH1F *jpsi_MCinput = new TH1F(Form("histoMCinput_cent%ito%i", arrayAvailableCentralityBins[iCent], arrayAvailableCentralityBins[iCent+1] ),"", ELEMENT_FOR_AE_PT_BIN-1, Pt_bin_for_Ae);
            jpsi_MCinput->Sumw2();
            
//            std::cout << arrayAvailableCentralityBins[iCent] << std::endl;
        }

        directoryJPsi->cd();
        TH1F *jpsi_MCinput = new TH1F(Form("histoMCinput_cent60to90"),"",ELEMENT_FOR_AE_PT_BIN-1, Pt_bin_for_Ae);
        jpsi_MCinput->Sumw2();
        
        for(int iCent=0; iCent< MAX_CENTRALITY_CLASS_SELECTION-2; iCent++)
        {
                   for(int iPt=0; iPt<(ELEMENT_FOR_AE_PT_BIN-1); iPt++)
                   {
                     TString strPtBins;
//              strPtBins.Form("Pt%gto%g",arrayPtBinsPhotonProduction[iPt],arrayPtBinsPhotonProduction[iPt+1]);
                     strPtBins.Form("Pt%gto%g",Pt_bin_for_Ae[iPt],Pt_bin_for_Ae[iPt+1]);
                     std::cout << strPtBins.Data() <<std::endl;             
                     const double display_pTMinBin = Pt_bin_for_Ae[iPt];
                     const double display_pTMaxBin = Pt_bin_for_Ae[iPt+1];

                     const int pTMinBin_statBox = (int)display_pTMinBin;
                     const int pTMaxBin_statBox = (int)display_pTMaxBin;

//                     int lower_limit_bin_to_be_integral = ((display_pTMinBin - pTMinBin) / value_each_pTbin )+ 1;
//                     int upper_limit_bin_to_be_integral = ((display_pTMaxBin - pTMinBin) / value_each_pTbin) + 1;

                     int lower_limit_bin_to_be_integral = hpT_gen[iCent]->FindBin(Pt_bin_for_Ae[iPt]);
                     int upper_limit_bin_to_be_integral = hpT_gen[iCent]->FindBin(Pt_bin_for_Ae[iPt+1]);                   
                     std::cout << "lower bound pT: " << Pt_bin_for_Ae[iPt] << ", bin: " << hpT_gen[iCent]->FindBin(Pt_bin_for_Ae[iPt]) << std::endl;
                     std::cout << "upper bound pT: " << Pt_bin_for_Ae[iPt+1] << ", bin: " << hpT_gen[iCent]->FindBin(Pt_bin_for_Ae[iPt+1]) << std::endl;

                     double error_phy_for_integral = 0;  double error_MC_for_integral = 0; double error_for_integral=0;
                     
                     double entries_hPt_gen = hpT_gen[iCent]->IntegralAndError(lower_limit_bin_to_be_integral, (upper_limit_bin_to_be_integral-1), error_MC_for_integral);
//                     double entries_hPt_gen = hM_gen[iCent][iPt]->Integral(1,invMassNBins);
//                     error_MC_for_integral = sqrt(entries_hPt_gen);
                     std::cout << entries_hPt_gen << " ";
                     std::cout << hM_gen[iCent][iPt]->Integral(1,invMassNBins) << ", ";
                     std::cout << hM_gen[iCent][iPt]->GetEntries() <<std::endl;
//                     double entries_hPt_gen = hpT_gen[iCent][iPt]->IntegralAndError(pTMinBin, pTNBins, error_MC_for_integral);
//                     double entries_hPt_gen = hpT_gen[iCent]->GetBinContent(iPt+1);
                     entries_hPt_gen = entries_hPt_gen / (Pt_bin_for_Ae[iPt+1]-Pt_bin_for_Ae[iPt]);
//                     error_MC_for_integral = hpT_gen[iCent]->GetBinError(iPt+1);
//                   double entries_hPt_rec = hPt_rec->IntegralAndError(lower_limit_bin_to_be_integral,upper_limit_bin_to_be_integral, error_phy_for_integral);
                     ((TH1F*) outputFile->Get(Form("JPsiMCinputHistos/histoMCinput_cent%ito%i", arrayAvailableCentralityBins[iCent], arrayAvailableCentralityBins[iCent+1]) ))->SetBinContent(iPt+1, entries_hPt_gen );
                     ((TH1F*) outputFile->Get(Form("JPsiMCinputHistos/histoMCinput_cent%ito%i", arrayAvailableCentralityBins[iCent], arrayAvailableCentralityBins[iCent+1]) ))->SetBinError(iPt+1, error_MC_for_integral );
                     
                      
                    }
         }

         outputFile->Write();

         inputFile->Close();
         inputFile2->Close();
         inputFile3->Close();
         inputFile4->Close();

         delete outputFile;

return 0;
}
