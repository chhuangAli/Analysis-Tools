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
  return par[0]*x[0] / TMath::Power(1. + TMath::Power(x[0]/par[1], par[2]), par[3] );
}

double y_shape(double *x, double *par)
{
  return par[0]*TMath::Exp(-(1./2.)*TMath::Power( ( ( x[0]-par[1])/par[2]), 2 ) );
}

double setWeight_Y(double Y, double *parRaw, double *parGen)
{

  TFormula *fRecEval = new TFormula("fRecEval","[0]*TMath::Exp(-0.5*TMath::Power( (x[0]-[1])/[2],2))");
//    TF1 *fRecEval = new TF1("fRecEval","[0]*x[0] / TMath::Power( (1 + [1]*x[0]*x[0]) , [2] )",0,15);
  fRecEval->SetParameter(0, parRaw[0] );
  fRecEval->SetParameter(1, parRaw[1] );
  fRecEval->SetParameter(2, parRaw[2] );


  TFormula *fGenEval = new TFormula("fGenEval","[0]*TMath::Exp(-0.5*TMath::Power( (x[0]-[1])/[2],2))");
  fGenEval->SetParameter(0, parGen[0] );
  fGenEval->SetParameter(1, parGen[1] );
  fGenEval->SetParameter(2, parGen[2] );

  double reco = fRecEval->Eval( Y );
  double simu = fGenEval->Eval( Y );

  if( (simu==0) || (reco/simu <= 0) )
  {
       return 1;
  }

  delete fRecEval;
  delete fGenEval;
  return reco/simu;
}


int main()
{

  const int ELEMENT_FOR_AE_PT_BIN =16;
  float pT_bin[ELEMENT_FOR_AE_PT_BIN]= {0,0.3,1,2,3,4,5,6,7,8,9,10,11,12,15,20};
  const int ELEMENT_FOR_CENT_BIN = 8;
  int cent_bin[ELEMENT_FOR_CENT_BIN] = {0,10,20,30,40,50,60,90};
  float y_bin[] = {-4., -3.75, -3.5, -3.25, -3., -2.75, -2.5};
  int number_y_bin = 7;
  
//
// Access to raw yield of Jpsi
//

  TFile *outputFile = new TFile(Form("WeightRawToMCinputHistos.root"),"recreate");
  TDirectory *directoryJPsiYield;
  directoryJPsiYield = outputFile->mkdir( Form("RawYieldAndFit") );
  TH1F *jpsi_yield;

  TDirectory *directoryJPsiMCinput;
  directoryJPsiMCinput = outputFile->mkdir( Form("MCinputAndFit") );

  TDirectory *directoryJPsiYieldRatio;
  directoryJPsiYieldRatio = outputFile->mkdir( Form("YieldRatio") );

  for(int iCent = 0; iCent < ELEMENT_FOR_CENT_BIN-1; iCent++)
  {
    directoryJPsiYield->cd();
    jpsi_yield = new TH1F(Form("histoYield_cent%ito%i", cent_bin[iCent], cent_bin[iCent+1] ),"", number_y_bin-1, y_bin);
    jpsi_yield->Sumw2();
  }

//  TString ReadFileLocation = Form("$HOME/cernbox/AnaSourceCode/anaJPsi_2015_2018PbPb/AxEff/cent_dep_tuning_input_shapes/Raw_data/Jpsi_yield/y_dep/JpsiYieldsHistos.root");
  TString ReadFileLocation = Form("$HOME/cernbox/AnaSourceCode/anaJPsi_2015_2018PbPb/AxEff/cent_dep_tuning_input_shapes/Raw_data/Jpsi_yield/y_dep/1st_tuning/JpsiYieldsHistos.root");
  TFile *inputFile_RawYield = new TFile(ReadFileLocation,"READ");

//  for(int iCent = 0; iCent < ELEMENT_FOR_CENT_BIN-1; iCent++)
  for(int iCent = 0; iCent < 2; iCent++) 
  {
    ((TH1F*) outputFile->Get(Form("RawYieldAndFit/histoYield_cent%ito%i", cent_bin[iCent], cent_bin[iCent+1]) ))->Add( (TH1F*)inputFile_RawYield->Get(Form("JPsiYieldsHistos/histoYield_cent%ito%i", cent_bin[iCent], cent_bin[iCent+1])) );    
  }
 
  directoryJPsiYield->Add(jpsi_yield );

//
// Access to the MC input of Jpsi
//

  TH1F *jpsi_MCinput;
  for(int iCent = 0; iCent < ELEMENT_FOR_CENT_BIN-1; iCent++)
  {
    directoryJPsiMCinput->cd();
    jpsi_MCinput = new TH1F(Form("histoMCinput_cent%ito%i", cent_bin[iCent], cent_bin[iCent+1] ),"", number_y_bin-1, y_bin);
    jpsi_MCinput->Sumw2();
  }

//  ReadFileLocation = Form("$HOME/cernbox/AnaSourceCode/anaJPsi_2015_2018PbPb/AxEff/cent_dep_tuning_input_shapes/MC_input/y_dep/JpsiMCinputHistos.root");
  ReadFileLocation = Form("$HOME/cernbox/AnaSourceCode/anaJPsi_2015_2018PbPb/AxEff/cent_dep_tuning_input_shapes/MC_input/1st_tuning_MC_input/y_dep/JpsiMCinputHistos.root");
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
       for (int iY = 0; iY < number_y_bin+1; iY++)
       {
           sum_points1 += ((TH1F*) outputFile->Get(Form("MCinputAndFit/histoMCinput_cent%ito%i", cent_bin[iCent], cent_bin[iCent+1]) ))->GetBinContent(iY+1);
           sum_points2 += ((TH1F*) outputFile->Get(Form("RawYieldAndFit/histoYield_cent%ito%i", cent_bin[iCent], cent_bin[iCent+1]) ))->GetBinContent(iY+1);
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
       jpsi_RatioRawToMCinput = new TH1F(Form("histoRatioRawToMCinput_cent%ito%i", cent_bin[iCent], cent_bin[iCent+1] ),"", number_y_bin-1, y_bin);
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

  bool fit_raw = false;
//  bool fit_raw = true;
 
  ROOT::Math::MinimizerOptions::SetDefaultMaxFunctionCalls(50000);

  
  int iCent = 1;

  double par1 = 0;
  double par1_err = 0;
  double par2 = 0;
  double par2_err = 0;
  double par3 = 0;
  double par3_err = 0;
  double par4 = 0;
  double par4_err = 0;
  double chi2 = 0;

  if(fit_raw)
  {
    TString opt = "IRSL";
    TCanvas *c2 = new TCanvas("c2","",200, 10, 800, 600);
//    c2->SetLogy();

    TF1 *tf_Nraw_y = new TF1("tf_Nraw_y", y_shape, y_bin[0], y_bin[number_y_bin], 3);
    tf_Nraw_y->SetParameter(0, 1000000);
    tf_Nraw_y->SetParameter(1, 0.1);
    tf_Nraw_y->SetParameter(2, 2.);
//    tf_Nraw_pt->SetParameter(3, 1.);

    tf_Nraw_y->SetParLimits(0, 0, 10000000);
    tf_Nraw_y->SetParLimits(1, -5, 5);
    tf_Nraw_y->SetParLimits(2, 0, 5);

    TFitResultPtr tfitResult;


    tfitResult = ((TH1F*) outputFile->Get(Form("RawYieldAndFit/histoYield_cent%ito%i", cent_bin[iCent], cent_bin[iCent+1]) ))->Fit(tf_Nraw_y, Form("%s", opt.Data()));


    par1 = tfitResult->Parameter(0);
    par1_err = tfitResult->ParError(0);
    par2 = tfitResult->Parameter(1);
    par2_err = tfitResult->ParError(1);
    par3 = tfitResult->Parameter(2);
    par3_err = tfitResult->ParError(2);
//    par4 = tfitResult->Parameter(3);
//    par4_err = tfitResult->ParError(3);
    chi2 = tf_Nraw_y->GetChisquare() / tf_Nraw_y->GetNDF();


    TPaveText* t1 = new TPaveText(0.6, 0.2, 0.8, 0.4, "NDC");
    t1->AddText(0.,0.,Form("A:%.2f #pm %.2f", par1, par1_err)  );
    t1->AddText(0.,0.,Form("n1:%.4f #pm %.4f", par2, par2_err)  );
    t1->AddText(0.,0.,Form("n2:%.4f #pm %.4f", par3, par3_err)  );
//    t1->AddText(0.,0.,Form("n2:%.4f #pm %.4f", par4, par4_err)  );
    t1->AddText(0.,0.,Form("#chi^{2}/nDoF: %.2f",chi2));

    ((TH1F*) outputFile->Get(Form("RawYieldAndFit/histoYield_cent%ito%i", cent_bin[iCent], cent_bin[iCent+1]) ))->Draw();
    tf_Nraw_y->Draw("SAME");
    t1->Draw("SAME");  
    c2->Print(Form("Fit_y_dep_RawJpsi_cent%ito%i.pdf", cent_bin[iCent], cent_bin[iCent+1]));

    delete tf_Nraw_y;
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
//    c3->SetLogy();

    TF1 *tf_Ngen_y = new TF1("tf_Ngen_y", y_shape, y_bin[0], y_bin[number_y_bin], 3);
    tf_Ngen_y->SetParameter(0, 1.009886e6);
//    tf_Ngen_y->SetParameter(1, 0);
//    tf_Ngen_y->SetParameter(2, 2.12568);
//    tf_Ngen_pt->SetParameter(3, 3.96363);
    if(iCent ==0)
    {
      tf_Ngen_y->FixParameter(1, -2.4306);
      tf_Ngen_y->FixParameter(2, 1.0798);
      
    }
    else if(iCent ==1)
    {
      tf_Ngen_y->FixParameter(1, -1.5252);
      tf_Ngen_y->FixParameter(2, 1.6015); 
    }

//    tf_Ngen_y->SetParLimits(0, 0, 10000000);
//    tf_Ngen_y->SetParLimits(1, -10, 10);  
//    tf_Ngen_y->SetParLimits(2, -5, 5);
//    tf_Ngen_pt->SetParLimits(3, 5, 300);

    TFitResultPtr tfitResult_MCinput;
    tfitResult_MCinput = ((TH1F*) outputFile->Get(Form("MCinputAndFit/histoMCinput_cent%ito%i", cent_bin[iCent], cent_bin[iCent+1]) ))->Fit(tf_Ngen_y, Form("%s", opt.Data()));

    par1 = tfitResult_MCinput->Parameter(0);
    par1_err = tfitResult_MCinput->ParError(0);
    par2 = tfitResult_MCinput->Parameter(1);
    par2_err = tfitResult_MCinput->ParError(1);
    par3 = tfitResult_MCinput->Parameter(2);
    par3_err = tfitResult_MCinput->ParError(2);
//    par4 = tfitResult_MCinput->Parameter(3);
//    par4_err = tfitResult_MCinput->ParError(3);
    chi2 = tf_Ngen_y->GetChisquare() / tf_Ngen_y->GetNDF();

    TPaveText* t1_MCinput = new TPaveText(0.6, 0.2, 0.8, 0.4, "NDC");
    t1_MCinput->AddText(0.,0.,Form("A:%.2f #pm %.2f", par1, par1_err)  );
    t1_MCinput->AddText(0.,0.,Form("n1:%.4f #pm %.4f", par2, par2_err)  );
    t1_MCinput->AddText(0.,0.,Form("n2:%.4f #pm %.4f", par3, par3_err)  );
//    t1_MCinput->AddText(0.,0.,Form("n2:%.4f #pm %.4f", par4, par4_err)  );
    t1_MCinput->AddText(0.,0.,Form("#chi^{2}/nDoF: %.2f",chi2));

  
    ((TH1F*) outputFile->Get(Form("MCinputAndFit/histoMCinput_cent%ito%i", cent_bin[iCent], cent_bin[iCent+1]) ))->Draw();
    tf_Ngen_y->Draw("SAME");
    t1_MCinput->Draw("SAME");
    c3->Print(Form("Fit_y_dep_MCinputJpsi_cent%ito%i.pdf", cent_bin[iCent], cent_bin[iCent+1]));
    delete tf_Ngen_y;
    delete t1_MCinput;
  }

  outputFile->Write();


return 0;
}
