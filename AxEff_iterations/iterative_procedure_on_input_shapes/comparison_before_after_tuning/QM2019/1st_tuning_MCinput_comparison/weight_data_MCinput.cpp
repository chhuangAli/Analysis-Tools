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
#include "TH2.h"
#include "TGraphAsymmErrors.h"

#include "TFitResult.h"
#include "Math/MinimizerOptions.h"

TString opt = "IRSL";

void myPadSetUp(TPad *currentPad, float currentLeft=0.11, float currentTop=0.04, float currentRight=0.04, float currentBottom=0.15);
void SetOptCanvas(TCanvas *canvas, Int_t divide=1);

double pt_shape(double *x, double *par)
{
  return par[0]*x[0] / pow(1.+pow(x[0]/par[1], par[2]), par[3] );
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

  TFile *outputFile = new TFile(Form("RawAndMCinputHistos.root"),"recreate");
  TDirectory *directoryJPsiYield;
  directoryJPsiYield = outputFile->mkdir( Form("RawYield") );
  TH1F *jpsi_yield;

  TDirectory *directoryJPsiMCinput;
  directoryJPsiMCinput = outputFile->mkdir( Form("MCinput") );

  TDirectory *directoryJPsiYieldRatio;
  directoryJPsiYieldRatio = outputFile->mkdir( Form("YieldRatio") );

  for(int iCent = 0; iCent < ELEMENT_FOR_CENT_BIN-1; iCent++)
  {
    directoryJPsiYield->cd();
    jpsi_yield = new TH1F(Form("histoYield_cent%ito%i", cent_bin[iCent], cent_bin[iCent+1] ),"", ELEMENT_FOR_AE_PT_BIN-1, pT_bin);
    jpsi_yield->Sumw2();
  }

//  TString ReadFileLocation = Form("$HOME/cernbox/AnaSourceCode/anaJPsi_2015_2018PbPb/AxEff/cent_dep_tuning_input_shapes/Raw_data/Jpsi_yield/JpsiYieldsHistos.root");
  TString ReadFileLocation = Form("$HOME/cernbox/AnaSourceCode/anaJPsi_2015_2018PbPb/AxEff/cent_dep_tuning_input_shapes/Raw_data/Jpsi_yield/1st_tuning/JpsiYieldsHistos.root");
  TFile *inputFile_RawYield = new TFile(ReadFileLocation,"READ");

  for(int iCent = 0; iCent < ELEMENT_FOR_CENT_BIN-1; iCent++)
  {
    ((TH1F*) outputFile->Get(Form("RawYield/histoYield_cent%ito%i", cent_bin[iCent], cent_bin[iCent+1]) ))->Add( (TH1F*)inputFile_RawYield->Get(Form("JPsiYieldsHistos/histoYield_cent%ito%i", cent_bin[iCent], cent_bin[iCent+1])) );    
  }
 
  directoryJPsiYield->Add(jpsi_yield );

//
// Access to the MC input of Jpsi
//

  TH1F *jpsi_MCinput;
  for(int iCent = 0; iCent < ELEMENT_FOR_CENT_BIN-1; iCent++)
  {
    directoryJPsiMCinput->cd();
    jpsi_MCinput = new TH1F(Form("histoMCinput_cent%ito%i", cent_bin[iCent], cent_bin[iCent+1] ),"", ELEMENT_FOR_AE_PT_BIN-1, pT_bin);
    jpsi_MCinput->Sumw2();
  }

  ReadFileLocation = Form("$HOME/cernbox/AnaSourceCode/anaJPsi_2015_2018PbPb/AxEff/cent_dep_tuning_input_shapes/MC_input/1st_tuning_MC_input/pT_dep/JpsiMCinputHistos.root");
  TFile *inputFile_MCinput = new TFile(ReadFileLocation,"READ");


  for(int iCent = 0; iCent < ELEMENT_FOR_CENT_BIN-1; iCent++)
  {
    ((TH1F*) outputFile->Get(Form("MCinput/histoMCinput_cent%ito%i", cent_bin[iCent], cent_bin[iCent+1]) ))->Add( (TH1F*)inputFile_MCinput->Get(Form("JPsiMCinputHistos/histoMCinput_cent%ito%i", cent_bin[iCent], cent_bin[iCent+1])) );
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
           sum_points1 += ((TH1F*) outputFile->Get(Form("MCinput/histoMCinput_cent%ito%i", cent_bin[iCent], cent_bin[iCent+1]) ))->GetBinContent(ipT+1);
           sum_points2 += ((TH1F*) outputFile->Get(Form("RawYield/histoYield_cent%ito%i", cent_bin[iCent], cent_bin[iCent+1]) ))->GetBinContent(ipT+1);
       }
       ((TH1F*) outputFile->Get(Form("MCinput/histoMCinput_cent%ito%i", cent_bin[iCent], cent_bin[iCent+1]) ))->Scale( sum_points2/sum_points1 );
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
      ((TH1F*) outputFile->Get(Form("YieldRatio/histoRatioRawToMCinput_cent%ito%i", cent_bin[iCent], cent_bin[iCent+1]) ))->Divide(((TH1F*) outputFile->Get(Form("RawYield/histoYield_cent%ito%i", cent_bin[iCent], cent_bin[iCent+1]) )), ((TH1F*) outputFile->Get(Form("MCinput/histoMCinput_cent%ito%i", cent_bin[iCent], cent_bin[iCent+1]) )));
  }

  directoryJPsiYieldRatio->Add( jpsi_RatioRawToMCinput);

  TGraphAsymmErrors *gr_dNdpT = new TGraphAsymmErrors( ELEMENT_FOR_AE_PT_BIN-1 );
  TGraphAsymmErrors *gr_N_jpsi_gen = new TGraphAsymmErrors( ELEMENT_FOR_AE_PT_BIN-1);


  double xVal=0, deltaxVal = 0, deltaxValSyst = 0;//cynthia
  int iCent = 1;

  for(int iPt = 0; iPt < ELEMENT_FOR_AE_PT_BIN-1; iPt++)
  {
    double dNdpT_over_Ae = 0;
    double dNdpT_over_Ae_err = 0;
    xVal = (pT_bin[iPt+1] + pT_bin[iPt])/2;
    deltaxVal = (pT_bin[iPt+1] - pT_bin[iPt])/2;

    dNdpT_over_Ae = ((TH1F*) outputFile->Get(Form("RawYield/histoYield_cent%ito%i", cent_bin[iCent], cent_bin[iCent+1]) ))->GetBinContent(iPt+1);
    dNdpT_over_Ae_err = ((TH1F*) outputFile->Get(Form("RawYield/histoYield_cent%ito%i", cent_bin[iCent], cent_bin[iCent+1]) ))->GetBinError(iPt+1);
    gr_dNdpT->SetPoint(iPt+1, xVal, dNdpT_over_Ae);
    gr_dNdpT->SetPointError(iPt+1, deltaxVal, deltaxVal, dNdpT_over_Ae_err, dNdpT_over_Ae_err);

    std::cout << dNdpT_over_Ae << " +/- " << dNdpT_over_Ae_err << std::endl;
  }
  
  for(int iPt = 0; iPt < ELEMENT_FOR_AE_PT_BIN-1; iPt++)
  {
//    int iCent = 1;

    double N_jpsi_gen = 0;
    double N_jpsi_gen_err = 0;

    xVal = (pT_bin[iPt+1] + pT_bin[iPt])/2;
    deltaxVal = (pT_bin[iPt+1] - pT_bin[iPt])/2;

    N_jpsi_gen = ((TH1F*) outputFile->Get(Form("MCinput/histoMCinput_cent%ito%i", cent_bin[iCent], cent_bin[iCent+1]) ))->GetBinContent(iPt+1); 
    N_jpsi_gen_err = ((TH1F*) outputFile->Get(Form("MCinput/histoMCinput_cent%ito%i", cent_bin[iCent], cent_bin[iCent+1]) ))->GetBinError(iPt+1);
    
    gr_N_jpsi_gen->SetPoint(iPt, xVal,  N_jpsi_gen);
    gr_N_jpsi_gen->SetPointError(iPt, deltaxVal, deltaxVal, N_jpsi_gen_err, N_jpsi_gen_err);
    std::cout << N_jpsi_gen << " +/- " << N_jpsi_gen_err << std::endl;
  }

  TGraphAsymmErrors *gr_ratio = new TGraphAsymmErrors( ELEMENT_FOR_AE_PT_BIN-1 );
  for(int iPt = 0; iPt < ELEMENT_FOR_AE_PT_BIN-1; iPt++)
  {
    double ratio = 0;
//    int iCent = 1;
    xVal = (pT_bin[iPt+1] + pT_bin[iPt])/2;

    ratio = ((TH1F*) outputFile->Get(Form("YieldRatio/histoRatioRawToMCinput_cent%ito%i", cent_bin[iCent], cent_bin[iCent+1]) ))->GetBinContent(iPt+1);
    gr_ratio->SetPoint(iPt, xVal,  ratio);
    std::cout << ratio << std::endl;
//    gr_ratio-
  }
  
  gr_ratio->SetName("gr_ratio");
  gr_ratio->SetMarkerStyle(20);

  gr_dNdpT->SetLineWidth(1);
  gr_dNdpT->SetMarkerStyle(20);
  gr_dNdpT->SetMarkerSize(1);
  gr_dNdpT->SetMarkerColor(kRed); 
  gr_dNdpT->SetLineColor(kRed); 

  //gr_dNdpT_syst->SetFillStyle(1000); 
  //gr_dNdpT_syst->SetFillColor(kRed-9); 
  //gr_dsigma_syst->SetLineStyle(1);  //cynthia
  //    gr_dsigma_syst->SetLineWidth(2);    //cynthia
  //gr_dNdpT_syst->SetLineColor(kRed);

  gr_N_jpsi_gen->SetLineWidth(1);//cynthia
  gr_N_jpsi_gen->SetMarkerStyle(20);//cynthia
  gr_N_jpsi_gen->SetMarkerSize(1);//cynthia
  gr_N_jpsi_gen->SetMarkerColor(kGray);//cynthia
  gr_N_jpsi_gen->SetLineColor(kGray);

  TCanvas *c2 = new TCanvas("c2","",200, 10, 800, 600);
  c2->Divide(1,2,0,0);
  SetOptCanvas(c2,2);


  TH2F *histo = new TH2F("histo","",100,0,20,100,1,2e6);
  histo->SetStats(0);
  //  histo->GetXaxis()->SetLabelSize(0.05);
  //histo->GetXaxis()->SetTitleSize(0.05);
  histo->GetXaxis()->SetTitleOffset(1.10);
  histo->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
  //histo->GetYaxis()->SetLabelSize(0.05);
  //histo->GetYaxis()->SetTitleSize(0.05);
  histo->GetYaxis()->SetTitleOffset(1.10);
  histo->GetYaxis()->SetTitle("dN/dp_{T}/A#epsilon");

  TH2F *hratio = new TH2F("hratio","",100,0,20,100, 0.7, 1.5);
  hratio->SetStats(0);
  hratio->GetXaxis()->SetTitleOffset(1.10);
  hratio->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
  hratio->GetYaxis()->SetTitleOffset(1.15);
  hratio->GetYaxis()->SetLabelSize(0.08);
  hratio->GetXaxis()->SetLabelSize(0.08);
  hratio->GetXaxis()->SetTitleSize(0.08);


  TPad *pad = (TPad*) c2->GetPad(1);
  pad->SetLogy();
  pad->cd();
  histo->Draw();
//  gr_dNdpT_syst->Draw("2SAME");
  gr_dNdpT->Draw("PSAME");
  gr_N_jpsi_gen->Draw("LPSAME");

  TLegend *leg2 = new TLegend(0.4,0.7,0.8,0.85);
  leg2->SetBorderSize(0);
  if(iCent == 0)
  {
    leg2->SetHeader("0 - 10\% centrality bin");
  }
  else if (iCent == 1)
  {
    leg2->SetHeader("10 - 20\% centrality bin");
  }
  //leg->SetTextSize(gStyle->GetTextSize()*0.8);
  leg2->SetLineColor(0);
  leg2->SetLineStyle(1);
  leg2->SetLineWidth(1);
  leg2->SetFillColor(0);
  leg2->SetFillStyle(1001);
  leg2->AddEntry(gr_dNdpT, "dN/dp_{T} * 1/A#epsilon","LPE");
//  leg2->AddEntry(gr_dNdpT_syst,"Systematic uncertainty","F");
  leg2->AddEntry(gr_N_jpsi_gen, " input MC (arbitrary norm.)","LPE");
  leg2->Draw();

  pad = (TPad*) c2->GetPad(2);
  pad->cd();
  hratio->Draw();
  gr_ratio->Draw("CPSAME");
  TLine *line = new TLine(0,1,20,1);
  line->SetLineStyle(2);
  line->Draw("same");

  TLegend *leg3 = new TLegend(0.6,0.85,0.8,0.95);
  leg3->SetBorderSize(0);
  //leg->SetTextSize(gStyle->GetTextSize()*0.8);
  leg3->SetLineColor(0);
  leg3->SetLineStyle(1);
  leg3->SetLineWidth(1);
  leg3->SetFillColor(0);
  leg3->SetFillStyle(1001);
  leg3->AddEntry(gr_ratio, "(dN/dp_{T} * 1/A#epsilon) /  input MC (arbitrary norm.)","LPE");
  leg3->Draw();

  c2->Print("dN_dpT_over_Ae.pdf");

/*
//
// Prepare to fit the raw yield of Jpsi
//  

  ROOT::Math::MinimizerOptions::SetDefaultMaxFunctionCalls(50000);

  TCanvas *c2 = new TCanvas("c2","",200, 10, 800, 600);
//  c2->SetLogy();

  double par1 = 0;
  double par1_err = 0;
  double par2 = 0;
  double par2_err = 0;
  double par3 = 0;
  double par3_err = 0;
  double par4 = 0;
  double par4_err = 0;
  double chi2 = 0;

  TF1 *tf_Nraw_pt = new TF1("tf_Nraw_pt", pt_shape, pT_bin[1], pT_bin[ELEMENT_FOR_AE_PT_BIN], 4);
  tf_Nraw_pt->SetParameter(0, 1000000);
  tf_Nraw_pt->SetParameter(1, 0.1);
  tf_Nraw_pt->SetParameter(2, 2.);
  tf_Nraw_pt->SetParameter(3, 1.);

  tf_Nraw_pt->SetParLimits(0, 100000, 10000000);
  tf_Nraw_pt->SetParLimits(1, 0, 10);
  tf_Nraw_pt->SetParLimits(2, 0, 5);


  TFitResultPtr tfitResult;

  int iCent = 0;

  tfitResult = ((TH1F*) outputFile->Get(Form("RawYieldAndFit/histoYield_cent%ito%i", cent_bin[iCent], cent_bin[iCent+1]) ))->Fit(tf_Nraw_pt, Form("%s", opt.Data()));

  par1 = tfitResult->Parameter(0);
  par1_err = tfitResult->ParError(0);
  par2 = tfitResult->Parameter(1);
  par2_err = tfitResult->ParError(1);
  par3 = tfitResult->Parameter(2);
  par3_err = tfitResult->ParError(2);
  par4 = tfitResult->Parameter(3);
  par4_err = tfitResult->ParError(3);
  chi2 = tf_Nraw_pt->GetChisquare() / tf_Nraw_pt->GetNDF();


  TPaveText* t1 = new TPaveText(0.2, 0.2, 0.4, 0.4, "NDC");
  t1->AddText(0.,0.,Form("A:%.2f #pm %.2f", par1, par1_err)  );
  t1->AddText(0.,0.,Form("B:%.2f #pm %.4f", par2, par2_err)  );
  t1->AddText(0.,0.,Form("n1:%.2f #pm %.4f", par3, par3_err)  );
  t1->AddText(0.,0.,Form("n2:%.2f #pm %.4f", par4, par4_err)  );
  t1->AddText(0.,0.,Form("#chi^{2}/nDoF: %.2f",chi2));

  ((TH1F*) outputFile->Get(Form("RawYieldAndFit/histoYield_cent%ito%i", cent_bin[iCent], cent_bin[iCent+1]) ))->Draw();
  tf_Nraw_pt->Draw("SAME");
  t1->Draw("SAME");  
  c2->Print(Form("FitAndWeight_RawJpsi_cent%ito%i.pdf", cent_bin[iCent], cent_bin[iCent+1]));

  delete tf_Nraw_pt;
  delete c2;
//
// Prepare to fit the MC input of Jpsi
//

  TCanvas *c3 = new TCanvas("c3","",200, 10, 800, 600);
  c3->SetLogy();

  TF1 *tf_Ngen_pt = new TF1("tf_Ngen_pt", pt_shape, pT_bin[0], pT_bin[ELEMENT_FOR_AE_PT_BIN], 4);
  tf_Ngen_pt->SetParameter(0, 2000000);
  tf_Ngen_pt->SetParameter(1, 1);
  tf_Ngen_pt->SetParameter(2, 1.);
  tf_Ngen_pt->SetParameter(3, 1);

  tf_Ngen_pt->SetParLimits(0, 100000, 10000000);
  tf_Ngen_pt->SetParLimits(1, 0, 10);  
  tf_Ngen_pt->SetParLimits(2, 0, 5);
//  tf_Ngen_pt->SetParLimits(3, 0, 10);


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
  chi2 = tf_Ngen_pt->GetChisquare() / tf_Ngen_pt->GetNDF();

  TPaveText* t1_MCinput = new TPaveText(0.2, 0.2, 0.4, 0.4, "NDC");
  t1_MCinput->AddText(0.,0.,Form("A:%.2f #pm %.2f", par1, par1_err)  );
  t1_MCinput->AddText(0.,0.,Form("B:%.2f #pm %.4f", par2, par2_err)  );
  t1_MCinput->AddText(0.,0.,Form("n1:%.2f #pm %.4f", par3, par3_err)  );
  t1_MCinput->AddText(0.,0.,Form("n2:%.2f #pm %.4f", par4, par4_err)  );
  t1_MCinput->AddText(0.,0.,Form("#chi^{2}/nDoF: %.2f",chi2));

  
  ((TH1F*) outputFile->Get(Form("MCinputAndFit/histoMCinput_cent%ito%i", cent_bin[iCent], cent_bin[iCent+1]) ))->Draw();
  tf_Ngen_pt->Draw("SAME");
  t1_MCinput->Draw("SAME");
  c3->Print("FitAndWeight_MCinputJpsi.pdf");
*/

  outputFile->Write();

//  delete tf_Ngen_pt;
//  delete t1_MCinput;
//  delete t1;

return 0;
}

void SetOptCanvas(TCanvas *canvas, Int_t divide) {

  if(!canvas) return;

  if (divide == 1){
    TPad *pad = (TPad*) canvas->GetPad(1);
    myPadSetUp(pad,0.11,0.04,0.04,0.15);
  }
  else   if(divide==2){
    TPad *pad = (TPad*) canvas->GetPad(1);
    myPadSetUp(pad,0.11,0.04,0.04,0.001);
    pad->SetPad(0.,0.3,1.0,1.0);
    pad = (TPad*) canvas->GetPad(2);
    myPadSetUp(pad,0.11,0.001,0.04,0.25);
    pad->SetPad(0.,0.,1.0,0.3);
  }
}

//____________________________________________________________
//cynthia some pad style
void myPadSetUp(TPad *currentPad, float currentLeft, float currentTop, float currentRight, float currentBottom){
  currentPad->SetLeftMargin(currentLeft);
  currentPad->SetTopMargin(currentTop);
  currentPad->SetRightMargin(currentRight);
  currentPad->SetBottomMargin(currentBottom);
  return;
}

