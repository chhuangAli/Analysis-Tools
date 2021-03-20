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
#include "TMath.h"


#include "TFitResult.h"
#include "Math/MinimizerOptions.h"


double pt_shape(double *x, double *par)
{
//  return par[0]*x[0] / TMath::Power(1. + TMath::Power(x[0]/par[1], par[2]), par[3] );
  return (par[0]*x[0] / TMath::Power(1. + TMath::Power(x[0]/par[1], par[2]), par[3] ) ) + (par[4]*TMath::Exp(-par[5]*x[0]));
}

double y_shape(double *x, double *par)
{
  return par[0]*TMath::Exp(-(1./2.)*TMath::Power( ( ( x[0]-par[1])/par[2]), 2 ) );
}

int main()
{

  const int ELEMENT_FOR_AE_PT_BIN =16;
  float pT_bin[ELEMENT_FOR_AE_PT_BIN]= {0,0.3,1,2,3,4,5,6,7,8,9,10,11,12,15,20};
  const int ELEMENT_FOR_CENT_BIN = 8;
  int cent_bin[ELEMENT_FOR_CENT_BIN] = {0,10,20,30,40,50,60,90};
  
//
// Access to raw yield of Jpsi
//

  TFile *outputFile = new TFile(Form("WeightRawToMCinputHistos.root"),"recreate");
  TDirectory *directoryJPsiYield;
  directoryJPsiYield = outputFile->mkdir( Form("RawYieldAndFit") );
  TH1F *jpsi_yield;
  TH1F *tf1_fill_random_jpsi_yield;
  TH1F *ratio_jpsi_yield_tf1;

  TDirectory *directoryJPsiMCinput;
  directoryJPsiMCinput = outputFile->mkdir( Form("MCinputAndFit") );

  TDirectory *directoryJPsiYieldRatio;
  directoryJPsiYieldRatio = outputFile->mkdir( Form("YieldRatio") );

  for(int iCent = 0; iCent < ELEMENT_FOR_CENT_BIN-1; iCent++)
  {
    directoryJPsiYield->cd();
    jpsi_yield = new TH1F(Form("histoYield_cent%ito%i", cent_bin[iCent], cent_bin[iCent+1] ),"", ELEMENT_FOR_AE_PT_BIN-1, pT_bin);
    jpsi_yield->Sumw2();

    tf1_fill_random_jpsi_yield = new TH1F(Form("histoYield_tf1_cent%ito%i", cent_bin[iCent], cent_bin[iCent+1] ),"", ELEMENT_FOR_AE_PT_BIN-1, pT_bin);
    tf1_fill_random_jpsi_yield->Sumw2();

    ratio_jpsi_yield_tf1 = new TH1F(Form("histo_ratio_Yield_tf1_cent%ito%i", cent_bin[iCent], cent_bin[iCent+1] ),"", ELEMENT_FOR_AE_PT_BIN-1, pT_bin);
    ratio_jpsi_yield_tf1->Sumw2();
  }

  TString ReadFileLocation = Form("$HOME/cernbox/AnaSourceCode/anaJPsi_2015_2018PbPb/AxEff/cent_dep_tuning_input_shapes/Raw_data/Jpsi_yield/JpsiYieldsHistos.root");  // before the tuning
//  TString ReadFileLocation = Form("$HOME/cernbox/AnaSourceCode/anaJPsi_2015_2018PbPb/AxEff/cent_dep_tuning_input_shapes/Raw_data/Jpsi_yield/1st_tuning/JpsiYieldsHistos.root");
  TFile *inputFile_RawYield = new TFile(ReadFileLocation,"READ");

  for(int iCent = 0; iCent < ELEMENT_FOR_CENT_BIN-1; iCent++)
  {
    ((TH1F*) outputFile->Get(Form("RawYieldAndFit/histoYield_cent%ito%i", cent_bin[iCent], cent_bin[iCent+1]) ))->Add( (TH1F*)inputFile_RawYield->Get(Form("JPsiYieldsHistos/histoYield_cent%ito%i", cent_bin[iCent], cent_bin[iCent+1])) );    
  }
 
  directoryJPsiYield->Add(jpsi_yield );

//
// Access to the MC input of Jpsi
//

  TH1F *jpsi_MCinput;
  TH1F *tf1_fill_random;
  TH1F *ratio_MCinput_tf1;
  for(int iCent = 0; iCent < ELEMENT_FOR_CENT_BIN-1; iCent++)
  {
    directoryJPsiMCinput->cd();
    jpsi_MCinput = new TH1F(Form("histoMCinput_cent%ito%i", cent_bin[iCent], cent_bin[iCent+1] ),"", ELEMENT_FOR_AE_PT_BIN-1, pT_bin);
    jpsi_MCinput->Sumw2();
  
    tf1_fill_random = new TH1F(Form("histoMCinput_tf1_cent%ito%i", cent_bin[iCent], cent_bin[iCent+1] ),"", ELEMENT_FOR_AE_PT_BIN-1, pT_bin);
    tf1_fill_random->Sumw2();

    ratio_MCinput_tf1 = new TH1F(Form("histo_ratio_MCinput_tf1_cent%ito%i", cent_bin[iCent], cent_bin[iCent+1] ),"", ELEMENT_FOR_AE_PT_BIN-1, pT_bin);
    ratio_MCinput_tf1->Sumw2();
  }

  ReadFileLocation = Form("$HOME/cernbox/AnaSourceCode/anaJPsi_2015_2018PbPb/AxEff/cent_dep_tuning_input_shapes/MC_input/JpsiMCinputHistos.root"); // before the tuning
//  ReadFileLocation = Form("$HOME/cernbox/AnaSourceCode/anaJPsi_2015_2018PbPb/AxEff/cent_dep_tuning_input_shapes/MC_input/1st_tuning_MC_input/pT_dep/JpsiMCinputHistos.root"); // first tuning
  TFile *inputFile_MCinput = new TFile(ReadFileLocation,"READ");


  for(int iCent = 0; iCent < ELEMENT_FOR_CENT_BIN-1; iCent++)
  {
    ((TH1F*) outputFile->Get(Form("MCinputAndFit/histoMCinput_cent%ito%i", cent_bin[iCent], cent_bin[iCent+1]) ))->Add( (TH1F*)inputFile_MCinput->Get(Form("JPsiMCinputHistos/histoMCinput_cent%ito%i", cent_bin[iCent], cent_bin[iCent+1])) );
  }

  directoryJPsiMCinput->Add( jpsi_MCinput );
//
// Set an arbitrary normalization factor
//

   double sum_points1=0;
   double sum_points2=0;

   for(int iCent = 0; iCent < ELEMENT_FOR_CENT_BIN-1; iCent++)
   {
       sum_points1 = 0;
       sum_points2 = 0;
       for (int ipT = 0; ipT < ELEMENT_FOR_AE_PT_BIN; ipT++)
       {
           sum_points1 += ((TH1F*) outputFile->Get(Form("MCinputAndFit/histoMCinput_cent%ito%i", cent_bin[iCent], cent_bin[iCent+1]) ))->GetBinContent(ipT+1);
           sum_points2 += ((TH1F*) outputFile->Get(Form("RawYieldAndFit/histoYield_cent%ito%i", cent_bin[iCent], cent_bin[iCent+1]) ))->GetBinContent(ipT+1);
       }
       ((TH1F*) outputFile->Get(Form("MCinputAndFit/histoMCinput_cent%ito%i", cent_bin[iCent], cent_bin[iCent+1]) ))->Scale( sum_points2/sum_points1 );
   }

//
// Make a ratio plot 
//
   TH1F *jpsi_RatioRawToMCinput;
   for(int iCent = 0; iCent < ELEMENT_FOR_CENT_BIN-1; iCent++)
   {
       directoryJPsiYieldRatio->cd();
       jpsi_RatioRawToMCinput = new TH1F(Form("histoRatioRawToMCinput_cent%ito%i", cent_bin[iCent], cent_bin[iCent+1] ),"", ELEMENT_FOR_AE_PT_BIN-1, pT_bin);
       jpsi_RatioRawToMCinput->Sumw2();
   }

  for(int iCent = 0; iCent < ELEMENT_FOR_CENT_BIN-1; iCent++)
  {
      ((TH1F*) outputFile->Get(Form("YieldRatio/histoRatioRawToMCinput_cent%ito%i", cent_bin[iCent], cent_bin[iCent+1]) ))->Divide(((TH1F*) outputFile->Get(Form("RawYieldAndFit/histoYield_cent%ito%i", cent_bin[iCent], cent_bin[iCent+1]) )), ((TH1F*) outputFile->Get(Form("MCinputAndFit/histoMCinput_cent%ito%i", cent_bin[iCent], cent_bin[iCent+1]) )));
  }

  directoryJPsiYieldRatio->Add( jpsi_RatioRawToMCinput);

//
// Prepare to fit the raw yield of Jpsi
//  

//  bool fit_raw = false;
  bool fit_raw = true;
 
  ROOT::Math::MinimizerOptions::SetDefaultMaxFunctionCalls(50000);
  
  const int iStep = 0;
  int iCent = 1;


  double par1 = 0;
  double par1_err = 0;
  double par2 = 0;
  double par2_err = 0;
  double par3 = 0;
  double par3_err = 0;
  double par4 = 0;
  double par4_err = 0;
  double par5 = 0;
  double par5_err = 0;
  double par6 = 0;
  double par6_err = 0;
  double chi2 = 0;

  if(fit_raw)
  {
    TString opt = "IRSL";
    TCanvas *c2 = new TCanvas("c2","",200, 10, 800, 600);
    c2->SetLogy();

    TF1 *tf_Nraw_pt = new TF1("tf_Nraw_pt", pt_shape, pT_bin[1], pT_bin[ELEMENT_FOR_AE_PT_BIN-1], 6);
    tf_Nraw_pt->SetParameter(0, 1000000);
    tf_Nraw_pt->SetParameter(1, 3.5);
    tf_Nraw_pt->SetParameter(2, 2.);
    tf_Nraw_pt->SetParameter(3, 3.9);

    tf_Nraw_pt->SetParameter(4, 100);
    tf_Nraw_pt->SetParameter(5, 1.0);

    tf_Nraw_pt->SetParLimits(0, 100000, 10000000);
    tf_Nraw_pt->SetParLimits(1, 0, 10);
    tf_Nraw_pt->SetParLimits(2, 0, 5);
    tf_Nraw_pt->SetParLimits(3, 0, 5);
 
    tf_Nraw_pt->SetParLimits(4, -100, 1e10);
    tf_Nraw_pt->SetParLimits(5, -5, 100);


    TFitResultPtr tfitResult;


    tfitResult = ((TH1F*) outputFile->Get(Form("RawYieldAndFit/histoYield_cent%ito%i", cent_bin[iCent], cent_bin[iCent+1]) ))->Fit(tf_Nraw_pt, Form("%s", opt.Data()));

    par1 = tfitResult->Parameter(0);
    par1_err = tfitResult->ParError(0);
    par2 = tfitResult->Parameter(1);
    par2_err = tfitResult->ParError(1);
    par3 = tfitResult->Parameter(2);
    par3_err = tfitResult->ParError(2);
    par4 = tfitResult->Parameter(3);
    par4_err = tfitResult->ParError(3);

    par5 = tfitResult->Parameter(4);
    par5_err = tfitResult->ParError(4);
    par6 = tfitResult->Parameter(5);
    par6_err = tfitResult->ParError(5);

    chi2 = tf_Nraw_pt->GetChisquare() / tf_Nraw_pt->GetNDF();

    TF1 *tf_Nraw_pt_fake = new TF1("tf_Nraw_pt", pt_shape, pT_bin[0], pT_bin[ELEMENT_FOR_AE_PT_BIN-1], 6);
    tf_Nraw_pt_fake->SetParameter(0, par1);
    tf_Nraw_pt_fake->SetParameter(1, par2);
    tf_Nraw_pt_fake->SetParameter(2, par3);
    tf_Nraw_pt_fake->SetParameter(3, par4);
    tf_Nraw_pt_fake->SetParameter(4, par5);
    tf_Nraw_pt_fake->SetParameter(5, par6);
    tf_Nraw_pt_fake->SetLineStyle(2);

    TPaveText* t1 = new TPaveText(0.2, 0.2, 0.4, 0.4, "NDC");
    t1->AddText(0.,0.,Form("A:%.2f #pm %.2f", par1, par1_err)  );
    t1->AddText(0.,0.,Form("B:%.4f #pm %.4f", par2, par2_err)  );
    t1->AddText(0.,0.,Form("n1:%.4f #pm %.4f", par3, par3_err)  );
    t1->AddText(0.,0.,Form("n2:%.4f #pm %.4f", par4, par4_err)  );
    t1->AddText(0.,0.,Form("C:%.4f #pm %.4f", par5, par5_err)  );
    t1->AddText(0.,0.,Form("n3:%.4f #pm %.4f", par6, par6_err)  );
    t1->AddText(0.,0.,Form("#chi^{2}/nDoF: %.2f",chi2));

    ((TH1F*) outputFile->Get(Form("RawYieldAndFit/histoYield_tf1_cent%ito%i", cent_bin[iCent], cent_bin[iCent+1]) ))->FillRandom("tf_Nraw_pt", 1e7);
    ((TH1F*) outputFile->Get(Form("RawYieldAndFit/histoYield_tf1_cent%ito%i", cent_bin[iCent], cent_bin[iCent+1]) ))->SetMarkerColor(2);
    ((TH1F*) outputFile->Get(Form("RawYieldAndFit/histoYield_tf1_cent%ito%i", cent_bin[iCent], cent_bin[iCent+1]) ))->SetLineColor(2);
    for(int iCent = 0; iCent < ELEMENT_FOR_CENT_BIN-1; iCent++)
    {
      for (int ipT = 0; ipT < (ELEMENT_FOR_AE_PT_BIN-1); ipT++)
      {
        double yval=0;
        double yval_err = 0;
        if(ipT == 0)
        {
          yval = ((TH1F*) outputFile->Get(Form("RawYieldAndFit/histoYield_tf1_cent%ito%i", cent_bin[iCent], cent_bin[iCent+1]) ))->GetBinContent(ipT+1);
          yval = yval/(pT_bin[ipT+1]-pT_bin[ipT]);
          yval_err = ((TH1F*) outputFile->Get(Form("RawYieldAndFit/histoYield_tf1_cent%ito%i", cent_bin[iCent], cent_bin[iCent+1]) ))->GetBinError(ipT+1);
          yval_err = yval*( yval_err/(((TH1F*) outputFile->Get(Form("RawYieldAndFit/histoYield_tf1_cent%ito%i", cent_bin[iCent], cent_bin[iCent+1]) ))->GetBinContent(ipT+1)));

          ((TH1F*) outputFile->Get(Form("RawYieldAndFit/histoYield_tf1_cent%ito%i", cent_bin[iCent], cent_bin[iCent+1]) ))->SetBinContent(ipT+1, yval);
          ((TH1F*) outputFile->Get(Form("RawYieldAndFit/histoYield_tf1_cent%ito%i", cent_bin[iCent], cent_bin[iCent+1]) ))->SetBinError(ipT+1, yval_err);
        }
        else 
        {
          yval = ((TH1F*) outputFile->Get(Form("RawYieldAndFit/histoYield_tf1_cent%ito%i", cent_bin[iCent], cent_bin[iCent+1]) ))->GetBinContent(ipT+1);
          yval = yval/(pT_bin[ipT+1]-pT_bin[ipT]);
          yval_err = ((TH1F*) outputFile->Get(Form("RawYieldAndFit/histoYield_tf1_cent%ito%i", cent_bin[iCent], cent_bin[iCent+1]) ))->GetBinError(ipT+1);
          yval_err = yval*( yval_err/(((TH1F*) outputFile->Get(Form("RawYieldAndFit/histoYield_tf1_cent%ito%i", cent_bin[iCent], cent_bin[iCent+1]) ))->GetBinContent(ipT+1)));

          ((TH1F*) outputFile->Get(Form("RawYieldAndFit/histoYield_tf1_cent%ito%i", cent_bin[iCent], cent_bin[iCent+1]) ))->SetBinContent(ipT+1, yval);
          ((TH1F*) outputFile->Get(Form("RawYieldAndFit/histoYield_tf1_cent%ito%i", cent_bin[iCent], cent_bin[iCent+1]) ))->SetBinError(ipT+1, yval_err);
        }

      }
    }

     double sum_points1=0;
     double sum_points2=0;

     for(int iCent = 0; iCent < ELEMENT_FOR_CENT_BIN-1; iCent++)
     {
       sum_points1 = 0;
       sum_points2 = 0;
       for (int ipT = 0; ipT < ELEMENT_FOR_AE_PT_BIN; ipT++)
       {
           sum_points1 += ((TH1F*) outputFile->Get(Form("RawYieldAndFit/histoYield_tf1_cent%ito%i", cent_bin[iCent], cent_bin[iCent+1]) ))->GetBinContent(ipT+1);
           sum_points2 += ((TH1F*) outputFile->Get(Form("RawYieldAndFit/histoYield_cent%ito%i", cent_bin[iCent], cent_bin[iCent+1]) ))->GetBinContent(ipT+1);
       }
       ((TH1F*) outputFile->Get(Form("RawYieldAndFit/histoYield_tf1_cent%ito%i", cent_bin[iCent], cent_bin[iCent+1]) ))->Scale( sum_points2/sum_points1 );

     }

    ((TH1F*) outputFile->Get(Form("RawYieldAndFit/histo_ratio_Yield_tf1_cent%ito%i", cent_bin[iCent], cent_bin[iCent+1]) ))->Add(((TH1F*) outputFile->Get(Form("RawYieldAndFit/histoYield_tf1_cent%ito%i", cent_bin[iCent], cent_bin[iCent+1]) )));
    ((TH1F*) outputFile->Get(Form("RawYieldAndFit/histo_ratio_Yield_tf1_cent%ito%i", cent_bin[iCent], cent_bin[iCent+1]) ))->Divide(((TH1F*) outputFile->Get(Form("RawYieldAndFit/histoYield_cent%ito%i", cent_bin[iCent], cent_bin[iCent+1]) )));

    ((TH1F*) outputFile->Get(Form("RawYieldAndFit/histoYield_tf1_cent%ito%i", cent_bin[iCent], cent_bin[iCent+1]) ))->Draw();
    ((TH1F*) outputFile->Get(Form("RawYieldAndFit/histoYield_cent%ito%i", cent_bin[iCent], cent_bin[iCent+1]) ))->Draw("SAME");
    tf_Nraw_pt->Draw("SAME");
    tf_Nraw_pt_fake->Draw("SAME");
    t1->Draw("SAME");  
    c2->Print(Form("FitAndWeight_RawJpsi_cent%ito%i.pdf(", cent_bin[iCent], cent_bin[iCent+1]));
    c2->SetLogy(0);
    ((TH1F*) outputFile->Get(Form("RawYieldAndFit/histo_ratio_Yield_tf1_cent%ito%i", cent_bin[iCent], cent_bin[iCent+1]) ))->Draw();
    TLine *line = new TLine(0,1,20,1);
    line->SetLineStyle(2);
    line->Draw("same");
    c2->Print(Form("FitAndWeight_RawJpsi_cent%ito%i.pdf)", cent_bin[iCent], cent_bin[iCent+1]));
    

    delete tf_Nraw_pt;
    delete t1;
    delete c2;
  }
  else
  {
//
// Prepare to fit the MC input of Jpsi
//
    TString opt = "IRSL";

    TCanvas *c3 = new TCanvas("c3","",200, 10, 800, 600);
    c3->SetLogy();

    double par0_rec_pt_0to10[] = {1030648.09, 2.8528, 2.8051, 2.4259, 215765.9365, 0.6850};
    double par0_rec_pt_10to20[] = {596170.22, 2.9441, 2.6666, 2.5048, 128583.4635, 0.6359};

    double weight0_pt=0;
    double weight0_y=0;
    double weight0=0;

    double weight1_pt=0;
    double weight1_y=0;
    double weight1=0;
/*
    if(iCent == 0)
    {
      weight0_pt = setWeight_Pt(, par0_rec_pt_0to10, par0_gen_pt_0to10);
      (TH1F*)inputFile_MCinput->Get(Form("JPsiMCinputHistos/histoMCinput_cent%ito%i", cent_bin[iCent], cent_bin[iCent+1]))-
    }
    else if (iCent == 1)
    {
      weight1_pt = setWeight_Pt(, par0_rec_pt_10to20, par0_gen_pt_10to20);
    }
*/

//    TF1 *tf_Ngen_pt = new TF1("tf_Ngen_pt", pt_shape, pT_bin[0], pT_bin[ELEMENT_FOR_AE_PT_BIN-1], 4);
     TF1 *tf_Ngen_pt = new TF1("tf_Ngen_pt", pt_shape, pT_bin[0], pT_bin[ELEMENT_FOR_AE_PT_BIN-1], 6);
    if(iStep == 1)
    {
      if(iCent == 0)
      {
        tf_Ngen_pt->SetParameter(0, par0_rec_pt_0to10[0]);
        tf_Ngen_pt->FixParameter(1, par0_rec_pt_0to10[1]);
        tf_Ngen_pt->FixParameter(2, par0_rec_pt_0to10[2]);
        tf_Ngen_pt->FixParameter(3, par0_rec_pt_0to10[3]);

        tf_Ngen_pt->SetParameter(4, par0_rec_pt_0to10[4]);
        tf_Ngen_pt->FixParameter(5, par0_rec_pt_0to10[5]);

        tf_Ngen_pt->SetParLimits(0, 100000, 5000000);
        tf_Ngen_pt->SetParLimits(4, -10000, 1000000);
      }
      else if (iCent == 1)
      {
        tf_Ngen_pt->SetParameter(0, par0_rec_pt_10to20[0]);
        tf_Ngen_pt->FixParameter(1, par0_rec_pt_10to20[1]);
        tf_Ngen_pt->FixParameter(2, par0_rec_pt_10to20[2]);
        tf_Ngen_pt->FixParameter(3, par0_rec_pt_10to20[3]);

        tf_Ngen_pt->SetParameter(4, par0_rec_pt_10to20[4]);
        tf_Ngen_pt->FixParameter(5, par0_rec_pt_10to20[5]);

        tf_Ngen_pt->SetParLimits(0, 100000, 5000000);
        tf_Ngen_pt->SetParLimits(4, -10000, 1000000);
      }

      
    }
    else 
    {
      tf_Ngen_pt->SetParameter(0, 1.00715e6);
      tf_Ngen_pt->SetParameter(1, 3.50274);
      tf_Ngen_pt->SetParameter(2, 1.93403);
      tf_Ngen_pt->SetParameter(3, 3.96363);

      tf_Ngen_pt->SetParameter(4, 1000);
      tf_Ngen_pt->SetParameter(5, 1.5);
//    tf_Ngen_pt->FixParameter(1, 3.50274);
//  tf_Ngen_pt->FixParameter(2, 1.93403);
//  tf_Ngen_pt->FixParameter(3, 3.96363);

      tf_Ngen_pt->SetParLimits(0, 100000, 5000000);
      tf_Ngen_pt->SetParLimits(1, 2, 200);  
      tf_Ngen_pt->SetParLimits(2, 0, 5);
      tf_Ngen_pt->SetParLimits(3, 0, 100);
 
      tf_Ngen_pt->SetParLimits(4, -10000, 100000);
      tf_Ngen_pt->SetParLimits(5, -5, 10);
    }


    TFitResultPtr tfitResult_MCinput;
    tfitResult_MCinput = ((TH1F*) outputFile->Get(Form("MCinputAndFit/histoMCinput_cent%ito%i", cent_bin[iCent], cent_bin[iCent+1]) ))->Fit(tf_Ngen_pt, Form("%s", opt.Data()));

    par1 = tfitResult_MCinput->Parameter(0);
    par1_err = tfitResult_MCinput->ParError(0);
    par2 = tfitResult_MCinput->Parameter(1);
    par2_err = tfitResult_MCinput->ParError(1);
    par3 = tfitResult_MCinput->Parameter(2);
    par3_err = tfitResult_MCinput->ParError(2);
    par4 = tfitResult_MCinput->Parameter(3);
    par4_err = tfitResult_MCinput->ParError(3);

    par5 = tfitResult_MCinput->Parameter(4);
    par5_err = tfitResult_MCinput->ParError(4);
    par6 = tfitResult_MCinput->Parameter(5);
    par6_err = tfitResult_MCinput->ParError(5);
    chi2 = tf_Ngen_pt->GetChisquare() / tf_Ngen_pt->GetNDF();

    TPaveText* t1_MCinput = new TPaveText(0.2, 0.2, 0.4, 0.4, "NDC");
    t1_MCinput->AddText(0.,0.,Form("A:%.2f #pm %.2f", par1, par1_err)  );
    t1_MCinput->AddText(0.,0.,Form("B:%.4f #pm %.4f", par2, par2_err)  );
    t1_MCinput->AddText(0.,0.,Form("n1:%.4f #pm %.4f", par3, par3_err)  );
    t1_MCinput->AddText(0.,0.,Form("n2:%.4f #pm %.4f", par4, par4_err)  );
    t1_MCinput->AddText(0.,0.,Form("C:%.4f #pm %.4f", par5, par5_err)  );
    t1_MCinput->AddText(0.,0.,Form("n3:%.4f #pm %.4f", par6, par6_err)  );
    t1_MCinput->AddText(0.,0.,Form("#chi^{2}/nDoF: %.2f",chi2));

    ((TH1F*) outputFile->Get(Form("MCinputAndFit/histoMCinput_tf1_cent%ito%i", cent_bin[iCent], cent_bin[iCent+1]) ))->FillRandom("tf_Ngen_pt",1e7);
    ((TH1F*) outputFile->Get(Form("MCinputAndFit/histoMCinput_tf1_cent%ito%i", cent_bin[iCent], cent_bin[iCent+1]) ))->SetMarkerColor(2);
    ((TH1F*) outputFile->Get(Form("MCinputAndFit/histoMCinput_tf1_cent%ito%i", cent_bin[iCent], cent_bin[iCent+1]) ))->SetLineColor(2);
    for(int iCent = 0; iCent < ELEMENT_FOR_CENT_BIN-1; iCent++)
    {
      for (int ipT = 0; ipT < (ELEMENT_FOR_AE_PT_BIN-1); ipT++)
      {
        double yval=0;
        double yval_err = 0;
        yval = ((TH1F*) outputFile->Get(Form("MCinputAndFit/histoMCinput_tf1_cent%ito%i", cent_bin[iCent], cent_bin[iCent+1]) ))->GetBinContent(ipT+1);
        yval = yval/(pT_bin[ipT+1]-pT_bin[ipT]);
        yval_err = ((TH1F*) outputFile->Get(Form("MCinputAndFit/histoMCinput_tf1_cent%ito%i", cent_bin[iCent], cent_bin[iCent+1]) ))->GetBinError(ipT+1);
        yval_err = yval*( yval_err/(((TH1F*) outputFile->Get(Form("MCinputAndFit/histoMCinput_tf1_cent%ito%i", cent_bin[iCent], cent_bin[iCent+1]) ))->GetBinContent(ipT+1)));

        ((TH1F*) outputFile->Get(Form("MCinputAndFit/histoMCinput_tf1_cent%ito%i", cent_bin[iCent], cent_bin[iCent+1]) ))->SetBinContent(ipT+1, yval);
        ((TH1F*) outputFile->Get(Form("MCinputAndFit/histoMCinput_tf1_cent%ito%i", cent_bin[iCent], cent_bin[iCent+1]) ))->SetBinError(ipT+1, yval_err);
        
      }
    }

     double sum_points1=0;
     double sum_points2=0;

     for(int iCent = 0; iCent < ELEMENT_FOR_CENT_BIN-1; iCent++)
     {
       sum_points1 = 0;
       sum_points2 = 0;
       for (int ipT = 0; ipT < ELEMENT_FOR_AE_PT_BIN; ipT++)
       {
           sum_points1 += ((TH1F*) outputFile->Get(Form("MCinputAndFit/histoMCinput_tf1_cent%ito%i", cent_bin[iCent], cent_bin[iCent+1]) ))->GetBinContent(ipT+1);
           sum_points2 += ((TH1F*) outputFile->Get(Form("MCinputAndFit/histoMCinput_cent%ito%i", cent_bin[iCent], cent_bin[iCent+1]) ))->GetBinContent(ipT+1);
       }
       ((TH1F*) outputFile->Get(Form("MCinputAndFit/histoMCinput_tf1_cent%ito%i", cent_bin[iCent], cent_bin[iCent+1]) ))->Scale( sum_points2/sum_points1 );
       
     }


    ((TH1F*) outputFile->Get(Form("MCinputAndFit/histo_ratio_MCinput_tf1_cent%ito%i", cent_bin[iCent], cent_bin[iCent+1]) ))->Add(((TH1F*) outputFile->Get(Form("MCinputAndFit/histoMCinput_tf1_cent%ito%i", cent_bin[iCent], cent_bin[iCent+1]) )));
    ((TH1F*) outputFile->Get(Form("MCinputAndFit/histo_ratio_MCinput_tf1_cent%ito%i", cent_bin[iCent], cent_bin[iCent+1]) ))->Divide(((TH1F*) outputFile->Get(Form("MCinputAndFit/histoMCinput_cent%ito%i", cent_bin[iCent], cent_bin[iCent+1]) )));
  
    ((TH1F*) outputFile->Get(Form("MCinputAndFit/histoMCinput_tf1_cent%ito%i", cent_bin[iCent], cent_bin[iCent+1]) ))->Draw();
    ((TH1F*) outputFile->Get(Form("MCinputAndFit/histoMCinput_tf1_cent%ito%i", cent_bin[iCent], cent_bin[iCent+1]) ))->SetAxisRange(1,1e7,"Y");
    ((TH1F*) outputFile->Get(Form("MCinputAndFit/histoMCinput_cent%ito%i", cent_bin[iCent], cent_bin[iCent+1]) ))->Draw("SAME");
    tf_Ngen_pt->Draw("SAME");
    t1_MCinput->Draw("SAME");
    c3->Print(Form("FitAndWeight_MCinputJpsi_cent%ito%i.pdf(", cent_bin[iCent], cent_bin[iCent+1]));
    c3->SetLogy(0);
    ((TH1F*) outputFile->Get(Form("MCinputAndFit/histo_ratio_MCinput_tf1_cent%ito%i", cent_bin[iCent], cent_bin[iCent+1]) ))->Draw();
     TLine *line = new TLine(0,1,20,1);
    line->SetLineStyle(2);
    line->Draw("same");
    c3->Print(Form("FitAndWeight_MCinputJpsi_cent%ito%i.pdf)", cent_bin[iCent], cent_bin[iCent+1]));
    delete tf_Ngen_pt;
    delete t1_MCinput;
  }

  outputFile->Write();


return 0;
}
