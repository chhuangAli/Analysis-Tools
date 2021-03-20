#include <iostream>
#include <fstream>
#include "TString.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TLine.h"


void read_numbers_in_txt(TString, double, double *[3]);

int main()
{

  const int ELEMENT_FOR_AE_PT_BIN =16;
  float Pt_bin_for_Ae[ELEMENT_FOR_AE_PT_BIN]= {0, 0.3, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 15, 20};
  int arrayAvailableCentralityBins[] ={0,10,20};

//  TFile *inputFile = new TFile("$HOME/cernbox/AnaSourceCode/anaJPsi_2015_2018PbPb/AxEff/cent_dep_tuning_input_shapes/Raw_data/AxEff/AxEff_LHC15o_18q_18r.root","READ");
//  TFile *inputFile2 = new TFile("$HOME/cernbox/AnaSourceCode/anaJPsi_2015_2018PbPb/AxEff/cent_dep_tuning_input_shapes/Raw_data/AxEff/after_tuning_MCinput/AxEff_LHC15o_18q_18r.root","READ");

  TFile *inputFile = new TFile("$HOME/cernbox/AnaSourceCode/anaJPsi_2015_2018PbPb/AxEff/cent_dep_tuning_input_shapes/Raw_data/AxEff/after_tuning_MCinput/AxEff_LHC15o_18q_18r.root","READ");
  TFile *inputFile2 = new TFile("$HOME/cernbox/AnaSourceCode/anaJPsi_2015_2018PbPb/AxEff/cent_dep_tuning_input_shapes/Raw_data/AxEff/after_2nd_tuning_MCinput/pT_dep/AxEff_LHC15o_18q_18r.root","READ");

  const int max_cent_bin = 2;
  TH1F* hAe_Pt[max_cent_bin];
  for(int iCent=0; iCent< max_cent_bin; iCent++ )
  {
    TString centrality_range = Form("_cent%dto%d", arrayAvailableCentralityBins[iCent],arrayAvailableCentralityBins[iCent+1] );
    hAe_Pt[iCent] = new TH1F(Form("hAe_Pt%d", arrayAvailableCentralityBins[iCent]),Form("A#epsilon as function of p_{T}; p_{T}(GeV/c); A#epsilon" ), (ELEMENT_FOR_AE_PT_BIN-1), Pt_bin_for_Ae);
    hAe_Pt[iCent]->SetStats(kFALSE);
    hAe_Pt[iCent]->Add((TH1F *)inputFile->Get( Form("hAe_Pt%s", centrality_range.Data()) ));   
  }

  TH1F* hAe_Pt_tuned[max_cent_bin];
  for(int iCent=0; iCent< max_cent_bin; iCent++ )
  {
    TString centrality_range = Form("_cent%dto%d", arrayAvailableCentralityBins[iCent],arrayAvailableCentralityBins[iCent+1] );
    hAe_Pt_tuned[iCent] = new TH1F(Form("hAe_Pt_tuned_%d", arrayAvailableCentralityBins[iCent]),Form("A#epsilon as function of p_{T}; p_{T}(GeV/c); A#epsilon" ),(ELEMENT_FOR_AE_PT_BIN-1), Pt_bin_for_Ae);
    hAe_Pt_tuned[iCent]->SetStats(kFALSE);
    hAe_Pt_tuned[iCent]->Add((TH1F *)inputFile2->Get( Form("hAe_Pt%s", centrality_range.Data()) ));
    hAe_Pt_tuned[iCent]->SetMarkerColor(3);
    hAe_Pt_tuned[iCent]->SetLineColor(3);

  }

  TH1F* hAe_Pt_rela_diff[max_cent_bin];
  for(int iCent=0; iCent< max_cent_bin; iCent++ )
  {
    TString centrality_range = Form("_cent%dto%d", arrayAvailableCentralityBins[iCent],arrayAvailableCentralityBins[iCent+1] );
    hAe_Pt_rela_diff[iCent] = new TH1F(Form("hAe_Pt_rela_diff_%d", arrayAvailableCentralityBins[iCent]),Form("A#epsilon as function of p_{T}; p_{T}(GeV/c); \%" ),(ELEMENT_FOR_AE_PT_BIN-1), Pt_bin_for_Ae);
    hAe_Pt_rela_diff[iCent]->SetStats(kFALSE);
    hAe_Pt_rela_diff[iCent]->Sumw2();
    hAe_Pt_rela_diff[iCent]->Add(hAe_Pt_tuned[iCent]);
    hAe_Pt_rela_diff[iCent]->Add(hAe_Pt[iCent], -1.0);
    hAe_Pt_rela_diff[iCent]->Divide(hAe_Pt[iCent]);
    hAe_Pt_rela_diff[iCent]->Scale(100);
    hAe_Pt_rela_diff[iCent]->GetYaxis()->SetTitle("#frac{A#epsilon_{2}-A#epsilon_{1}}{A#epsilon_{1}} %");
  }

 

  TCanvas *c2 = new TCanvas("c2","",200, 10, 800, 600);
  hAe_Pt[0]->Draw();
  hAe_Pt_tuned[0]->Draw("SAME");
  c2->Print("Ae_pT_comparison_0to10.pdf");

  hAe_Pt_rela_diff[0]->Draw();
  TLine *line = new TLine(0,0,20,0);
  line->SetLineStyle(2);
  line->Draw("same");
  c2->Print("Ae_pT_rela_diff_0to10.pdf");

  hAe_Pt[1]->Draw();
  hAe_Pt_tuned[1]->Draw("SAME");
  c2->Print("Ae_pT_comparison_10to20.pdf");

  hAe_Pt_rela_diff[1]->Draw();
  line->Draw("same");
  c2->Print("Ae_pT_rela_diff_10to20.pdf");

  inputFile->Close();
  inputFile2->Close();

return 0;
}

void read_numbers_in_txt(TString path_to_file, double number_pT_bin, double *output_variable[3])
{
    std::ifstream input_file(path_to_file.Data() );
    if(input_file.is_open() )
    {
       std::cout << "Opening file: " << path_to_file.Data() << std::endl;
       for(int iPt=0; iPt < number_pT_bin; iPt++)
       {
            input_file >> output_variable[iPt][0] >> output_variable[iPt][1] >> output_variable[iPt][2];
            std::cout << "The values are read: " << output_variable[iPt][0] << " " << output_variable[iPt][1] << " " << output_variable[iPt][2] << std::endl;
       }
    }
    else std::cout << "Unable to open file: " << path_to_file.Data() << std::endl;
    input_file.close();
}
