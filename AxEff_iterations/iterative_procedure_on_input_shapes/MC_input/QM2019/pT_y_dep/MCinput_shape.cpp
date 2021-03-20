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


        int number_y_bin = 4;
//  float y_bin[] = {2.5, 3, 3.5, 4 };
        float y_bin[] = {-4, -3.5, -3, -2.5};

        const int number_pT_bin = 5;
        float pT_bin[] = {0.3,2,4,6,12,20};

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


        inputSimFile = Form("$HOME/local/anaedPbPbSimData2016/LHC16e2_jpsi_ttree/CINTEvent_MCinput_weight_AnaResults/PbPbMC_ResultsHistos_LHC16e2.root");
        inputSimFile2 = Form("$HOME/local/anaedPbPbSimData2016/LHC16e2_plus_jpsi_ttree/CINTEvent_MCinput_weight_AnaResults/PbPbMC_ResultsHistos_LHC16e2_plus.root");
        inputFile = TFile::Open( inputSimFile.Data() );
        inputFile2 = TFile::Open( inputSimFile2.Data() );


        inputSimFile3 = ("$HOME/local/anaedPbPbSim2018/LHC19a2_jpsi_ttree/CINTEvent_MCinput_weight_AnaResults/PbPbMC_ResultsHistos_LHC18q.root");
        inputFile3 = TFile::Open( inputSimFile3.Data() );

        inputSimFile4 = ("$HOME/local/anaedPbPbSim2018/LHC19a2_jpsi_ttree/CINTEvent_MCinput_weight_AnaResults/PbPbMC_ResultsHistos_LHC18r.root");
        inputFile4 = TFile::Open( inputSimFile4.Data() );

        TString outFile;
        outFile = Form("JpsiMCinputHistos.root");

        TH1F *h_gen[MAX_CENTRALITY_CLASS_SELECTION];
        TH1F *hY_gen[MAX_CENTRALITY_CLASS_SELECTION][number_y_bin];


        for(int iCent=0; iCent< 2; iCent++)
        {                                  
                  h_gen[iCent] = new TH1F(Form("hRapidityUnlikeSign_Cent%gto%g", arrayAvailableCentralityBins[iCent], arrayAvailableCentralityBins[iCent+1]),"", NBins, MinBin, MaxBin);

                  for( int i=0; i < (number_pT_bin-1); i++)
                  {                                
                    hY_gen[iCent][i] = new TH1F(Form("hY_gen_%i_%i", iCent, i),"hY_gen", (number_y_bin-1), y_bin);
                    hY_gen[iCent][i]->Sumw2();                
                  }

                 
        }
//        h_gen[MAX_CENTRALITY_CLASS_SELECTION-1] = new TH1F(Form("hRapidityUnlikeSign_Cent60to90"),"", NBins, MinBin, MaxBin);

        for(int iCent=0; iCent<2; iCent++ )
        {

                  TString centrality_range = Form("_Cent%dto%d", arrayAvailableCentralityBins[iCent],arrayAvailableCentralityBins[iCent+1] );
                  std::cout << centrality_range.Data() << std::endl;


                  if(iCent >= 6)
                  {
                    h_gen[MAX_CENTRALITY_CLASS_SELECTION-1]->Add(h_gen[iCent]);
                  }

                 for( int i=0; i < (number_pT_bin-1); i++)
                 {
                   TString range = Form("_Pt%gto%g", pT_bin[i], pT_bin[i+1] );
                   hY_gen[iCent][i]->Add( (TH1F*)inputFile->Get( Form("DimuonHistosGenJpsi/GenJpsi_Y%s%s", centrality_range.Data(), range.Data() )) );
                   hY_gen[iCent][i]->Add( (TH1F*)inputFile2->Get( Form("DimuonHistosGenJpsi/GenJpsi_Y%s%s", centrality_range.Data(), range.Data() )) );
                   hY_gen[iCent][i]->Add( (TH1F*)inputFile3->Get( Form("DimuonHistosGenJpsi/GenJpsi_Y%s%s", centrality_range.Data(), range.Data() )) );
                   hY_gen[iCent][i]->Add( (TH1F*)inputFile4->Get( Form("DimuonHistosGenJpsi/GenJpsi_Y%s%s", centrality_range.Data(), range.Data() )) );                   
                 }
        }


        TFile *outputFile = new TFile(outFile.Data(),"recreate");
        TDirectory *directoryJPsi;
        directoryJPsi = outputFile->mkdir( Form("JPsiMCinputHistos") );

        for(int iCent = 0; iCent < MAX_CENTRALITY_CLASS_SELECTION-2; iCent++)
        {
            for( int i=0; i < (number_pT_bin-1); i++)
            {
              TString range = Form("_Pt%gto%g", pT_bin[i], pT_bin[i+1] );
              directoryJPsi->cd();
              TH1F *jpsi_MCinput = new TH1F(Form("histoMCinput_cent%ito%i%s", arrayAvailableCentralityBins[iCent], arrayAvailableCentralityBins[iCent+1], range.Data() ),"", (number_y_bin-1), y_bin);
              jpsi_MCinput->Sumw2();
            }
            
//            std::cout << arrayAvailableCentralityBins[iCent] << std::endl;
        }

        directoryJPsi->cd();
        TH1F *jpsi_MCinput = new TH1F(Form("histoMCinput_cent60to90"),"", number_y_bin, y_bin);
        jpsi_MCinput->Sumw2();
        
        for(int iCent = 0; iCent < 2; iCent++)
        {
            for( int i=0; i < (number_pT_bin-1); i++)
            {
              TString range = Form("_Pt%gto%g", pT_bin[i], pT_bin[i+1] );
              ((TH1F*) outputFile->Get(Form("JPsiMCinputHistos/histoMCinput_cent%ito%i%s", arrayAvailableCentralityBins[iCent], arrayAvailableCentralityBins[iCent+1],range.Data()) ))->Add(hY_gen[iCent][i] );
            }
            
//            std::cout << arrayAvailableCentralityBins[iCent] << std::endl;
        }

        
         outputFile->Write();

         inputFile->Close();
         inputFile2->Close();
         inputFile3->Close();
         inputFile4->Close();

         delete outputFile;

return 0;
}
