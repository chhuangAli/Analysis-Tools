#include <iostream>
#include <fstream>

#include "TFile.h"
#include "TString.h"
#include "TCanvas.h"
#include "TH1.h"
#include "TGraphAsymmErrors.h"
#include "TH2F.h"
#include "TLegend.h"

void makeup_tgraph(TGraphAsymmErrors *, TGraphAsymmErrors *);
void drawing_histograms(TGraphAsymmErrors* , TGraphAsymmErrors* );

int main()
{
//  const int ELEMENT_FOR_AE_PT_BIN =16;
//  float pT_bin[ELEMENT_FOR_AE_PT_BIN]= {0,0.3,1,2,3,4,5,6,7,8,9,10,11,12,15,20};

//  int number_y_bin = 7;
//  float y_bin[] ={2.5,2.75,3,3.25,3.5,3.75,4};
//  float y_bin[] = {-4., -3.75, -3.5, -3.25, -3., -2.75, -2.5};
  int number_y_bin = 4;
//  float y_bin[] = {2.5, 3, 3.5, 4 };
  float y_bin[] = {-4, -3.5, -3, -2.5};

  const int number_pT_bin = 5;
  float pT_bin[] = {0.3,2,4,6,12,20};
  
  const int ELEMENT_FOR_CENT_BIN = 3;
  int cent_bin[ELEMENT_FOR_CENT_BIN] = {0,10,20};

  TString path_to_directory = Form("/home/chunlu/cernbox/AnaSourceCode/anaJPsi_2015_2018PbPb/AxEff/cent_dep_tuning_input_shapes/Raw_data/Fit_results/");
  TString Njpsi_file_name = "jpsi_number_ppTail_VWG_2p2_4p5.txt";
//  TString path_to_Njpsi;
   
//  TString ReadFileLocation = Form("$HOME/cernbox/AnaSourceCode/anaJPsi_2015_2018PbPb/AxEff/cent_dep_tuning_input_shapes/Raw_data/AxEff/AxEff_y/AxEff_LHC15o_18q_18r.root");
  TString ReadFileLocation = Form("$HOME/cernbox/AnaSourceCode/anaJPsi_2015_2018PbPb/AxEff/cent_dep_tuning_input_shapes/Raw_data/AxEff/AxEff_double_diff/AxEff_LHC15o_18q_18r.root");
  TFile *inputFile_AxEff = new TFile(ReadFileLocation,"READ");
  if(inputFile_AxEff->IsOpen()) std::cout << "read" << std::endl;

  ReadFileLocation = Form("/home/chunlu/cernbox/AnaSourceCode/anaJPsi_2015_2018PbPb/AxEff/cent_dep_tuning_input_shapes/Raw_data/AxEff/AxEff_double_diff/AxEff_double_diff.txt");
  std::ifstream inputfile(ReadFileLocation.Data());
  

  const int CENTRALITY_TUNING = 2;
  TH1F* hAe_cent_dep[CENTRALITY_TUNING][number_pT_bin];

  for(int iCent = 0; iCent < (CENTRALITY_TUNING); iCent++)
  {
      for(int iPt =0; iPt < (number_pT_bin); iPt++)
      {
        hAe_cent_dep[iCent][iPt] = new TH1F(Form("hAe_cent%ito%i_pT%gto%g", cent_bin[iCent], cent_bin[iCent+1], pT_bin[iPt], pT_bin[iPt+1]),Form("Ae as function of y cent%ito%i; y; A#epsilon", cent_bin[iCent], cent_bin[iCent+1]), (number_y_bin-1), y_bin);    
        hAe_cent_dep[iCent][iPt]->Sumw2();
      }
  }

  double AxEff[CENTRALITY_TUNING][number_pT_bin][number_y_bin] = {0};
  double AxEff_err[CENTRALITY_TUNING][number_pT_bin][number_y_bin] = {0};

  for(int iCent = 0; iCent < (CENTRALITY_TUNING); iCent++)
  {
    TString cent_range = Form("_cent%ito%i", cent_bin[iCent], cent_bin[iCent+1]);    
    for(int iPt =0; iPt < (number_pT_bin); iPt++)
    {
        TString pT_range = Form("_pT%gto%g", pT_bin[iPt], pT_bin[iPt+1]);
        std::cout << Form("hAe_cent%ito%i_pT%gto%g", cent_bin[iCent], cent_bin[iCent+1], pT_bin[iPt], pT_bin[iPt+1])  << std::endl;
        for(int iY =0; iY < (number_y_bin-1); iY++)
        {
//        hAe_cent_dep[iCent][iPt]->Add( ((TH1F*)inputFile_AxEff->Get(Form("hAe%s%s",cent_range.Data(), pT_range.Data() ) )) );

          inputfile >> AxEff[iCent][iPt][iY] >> AxEff_err[iCent][iPt][iY];
          hAe_cent_dep[iCent][iPt]->SetBinContent(iY+1, AxEff[iCent][iPt][iY]);
          hAe_cent_dep[iCent][iPt]->SetBinError(iY+1, AxEff_err[iCent][iPt][iY]);
          std::cout <<  AxEff[iCent][iPt][iY]  << " " << AxEff_err[iCent][iPt][iY] << std::endl;
          std::cout <<  hAe_cent_dep[iCent][iPt]->GetBinContent(iY+1)   << " " << hAe_cent_dep[iCent][iPt]->GetBinError(iY+1) << std::endl;
        }
    }
  }

  double Njpsi[2]={0};
  std::ofstream dNdpT_output_stream("d2NdpTdy_func_y_PbPb.txt");

  TGraphAsymmErrors *gr_dN_dy = new TGraphAsymmErrors( number_y_bin );
  TGraphAsymmErrors *gr_dN_dy_syst = new TGraphAsymmErrors( number_y_bin );

  TFile *outputFile = new TFile(Form("JpsiYieldsHistos.root"),"recreate");
  TDirectory *directoryJPsi;
  directoryJPsi = outputFile->mkdir( Form("JPsiYieldsHistos") );

  for(int iCent = 0; iCent < ELEMENT_FOR_CENT_BIN-1; iCent++)  
  {        
    for(int iPt =0; iPt < (number_pT_bin); iPt++)
    {
      directoryJPsi->cd();
      TH1F *jpsi_yields = new TH1F(Form("histoYield_cent%ito%i_Pt%gto%g", cent_bin[iCent], cent_bin[iCent+1], pT_bin[iPt], pT_bin[iPt+1] ),"", (number_y_bin-1), y_bin);
      jpsi_yields->Sumw2();   
      directoryJPsi->Add(hAe_cent_dep[iCent][iPt]);
    }
  }

  for(int iCent = 0; iCent < ELEMENT_FOR_CENT_BIN-1; iCent++)
//  for(int iCent = 2; iCent < 3; iCent++)
  {
    TString cent_range = Form("cent%ito%i/", cent_bin[iCent], cent_bin[iCent+1]);    
    for(int ipT=0; ipT < (number_pT_bin); ipT++)
    {
      for(int iY=0; iY < (number_y_bin-1); iY++)
      {
        TString pt_range = Form("pT%gto%g_", pT_bin[ipT], pT_bin[ipT+1]);
        TString y_range = Form("y%gto%g/", fabs(y_bin[iY+1]), fabs(y_bin[iY]) );    
        TString path_to_Njpsi = Form("%s%sdouble_diff/%s%s%s", path_to_directory.Data(), cent_range.Data(), pt_range.Data(), y_range.Data(), Njpsi_file_name.Data() );

        std::ifstream input_JPsi_number(path_to_Njpsi.Data() );
        if (input_JPsi_number.is_open())
        {
          input_JPsi_number >> Njpsi[0] >> Njpsi[1];
          input_JPsi_number.close();
//          std::cout << path_to_Njpsi.Data() << std::endl;
//      dNdpT_output_stream << Njpsi[0] << " " << Njpsi[1] << " " << Njpsi[2] << std::endl;

          double xVal =0;
          double xValerr = 0, xValSysterr=0;
          double yVal=0;
          double yValerr =0, yValSysterr =0;

          xVal = (y_bin[iY+1] + y_bin[iY]) / 2;
          xValerr = (y_bin[iY+1] - y_bin[iY]) / 2.;
          xValSysterr = 0.3*xValerr;
 
          yVal = Njpsi[0] / (y_bin[iY+1]-y_bin[iY]);
          yValerr = yVal * Njpsi[1] / Njpsi[0];

          yValSysterr =  yVal * Njpsi[2] / Njpsi[0];

          std::cout << y_bin[iY] << " " << y_bin[iY+1] << " " << (y_bin[iY+1]-y_bin[iY]) << std::endl;
          std::cout << yVal << " +/- " << yValerr << std::endl;

          gr_dN_dy->SetPoint(iY+1, xVal, yVal );
          gr_dN_dy->SetPointError(iY+1, xValerr, xValerr, yValerr, yValerr);

          ((TH1F*) outputFile->Get(Form("JPsiYieldsHistos/histoYield_cent%ito%i_Pt%gto%g", cent_bin[iCent], cent_bin[iCent+1], pT_bin[ipT], pT_bin[ipT+1] ) ))->SetBinContent(iY+1, yVal);
          ((TH1F*) outputFile->Get(Form("JPsiYieldsHistos/histoYield_cent%ito%i_Pt%gto%g", cent_bin[iCent], cent_bin[iCent+1], pT_bin[ipT], pT_bin[ipT+1] ) ))->SetBinError(iY+1, yValerr);

          gr_dN_dy_syst->SetPoint(iY+1, xVal, yVal );
          gr_dN_dy_syst->SetPointError(iY+1, xValSysterr, xValSysterr, yValSysterr, yValSysterr );

          std::cout << yVal << std::endl;
          dNdpT_output_stream << yVal << " " << gr_dN_dy->GetErrorY(iY+1) << " " << gr_dN_dy_syst->GetErrorY(iY+1) <<     std::endl;

        }
        else
        {
          std::cout << "Unable to open file: " << path_to_Njpsi.Data() << std::endl;
        }
    
      }
      ((TH1F*) outputFile->Get(Form("JPsiYieldsHistos/histoYield_cent%ito%i_Pt%gto%g", cent_bin[iCent], cent_bin[iCent+1], pT_bin[ipT], pT_bin[ipT+1] ) ))->Divide(hAe_cent_dep[iCent][ipT]);
    }
  } 
//  makeup_tgraph(gr_dN_dy, gr_dN_dy_syst);
//  drawing_histograms(gr_dN_dy, gr_dN_dy_syst);
 
  outputFile->Write();

  delete gr_dN_dy;
  delete gr_dN_dy_syst;
//  delete jpsi_yields;
  delete outputFile;
  inputFile_AxEff->Close();
  delete inputFile_AxEff;


return 0;
}

void makeup_tgraph(TGraphAsymmErrors *gr, TGraphAsymmErrors *gr_syst)
{
//    gr->SetLineColor(2);
//  gr->SetLineWidth(2);

//    gr_syst->SetFillStyle(1000);
//    gr_syst->SetFillColor(46);
//    gr_syst->SetLineWidth(2);
//    gr_syst->SetLineStyle(1);

    gr->SetLineWidth(1);//cynthia
    gr->SetMarkerStyle(20);//cynthia
    gr->SetMarkerSize(0.5);//cynthia
    gr->SetMarkerColor(kRed);//cynthia
    gr->SetLineColor(kRed);//cynthia

    gr_syst->SetFillStyle(1000);//cynthia
    gr_syst->SetFillColor(kRed-9);//cynthia
  //gr_dsigma_syst->SetLineStyle(1); //cynthia
  //    gr_dsigma_syst->SetLineWidth(2);    //cynthia
    gr_syst->SetLineColor(kRed);
}

void drawing_histograms(TGraphAsymmErrors *gr_dN_dpT, TGraphAsymmErrors *gr_dN_dpT_syst)
{
  TH2F *histo = new TH2F("histo","",100,-4,-2,100, 1.e4, 5.e5);
  histo->GetYaxis()->SetLabelSize(0.04);
  histo->GetYaxis()->SetTitleSize(0.04);
  histo->GetYaxis()->SetTitleOffset(1.15);
  histo->SetStats(0);
  histo->GetYaxis()->SetTitle("dN/d#it{y} ");
  histo->GetXaxis()->SetTitle("#it{y} ");

  TCanvas *tc_jpsi_pt = new TCanvas("c0","c0",600,500);
  tc_jpsi_pt->SetLogy();
  histo->Draw();
  gr_dN_dpT_syst->Draw("2SAME");
  gr_dN_dpT->Draw("PSAME");

  TLegend *leg2 = new TLegend(0.4,0.7,0.8,0.85);
  leg2->SetBorderSize(0);
  //leg->SetTextSize(gStyle->GetTextSize()*0.8);
  leg2->SetLineColor(0);
  leg2->SetLineStyle(1);
  leg2->SetLineWidth(1);
  leg2->SetFillColor(0);
  leg2->SetFillStyle(1001);
//  leg2->AddEntry(gr_dsigma, "ALICE Inclusive J/#psi -> #mu^{+}#mu^{-}, 2.5 < y < 4, Lint = 1.23 pb","LPE");
  leg2->AddEntry(gr_dN_dpT, "Statistic uncertainty","LPE");
//  leg2->AddEntry(gr_dN_dpT_syst,"Systematic uncertainty","F");
  leg2->Draw();

  TString output_plot_name = "dNdy_func_y_PbPb.pdf";
  tc_jpsi_pt->Print( output_plot_name.Data() );


  delete histo;
  delete tc_jpsi_pt;
}

