// - compute the dN/dpT * 1/Ae 
// - compute the N^gen / Ae

#include <iostream>
#include <fstream>

#include "TCanvas.h"
#include "TH1.h"
#include "TH2.h"
#include "TGraphAsymmErrors.h"
#include "TLegend.h"
#include "TString.h"
#include "TLine.h"
#include "TFile.h"

void read_values_stat_sysm(TString, int, double *[3]);
void read_values_stat(TString, int, double *[2]);
double error_propagation(int chose, double a, double delta_a, double b, double delta_b, double c);
void myPadSetUp(TPad *currentPad, float currentLeft=0.11, float currentTop=0.04, float currentRight=0.04, float currentBottom=0.15);
void SetOptCanvas(TCanvas *canvas, Int_t divide=1);

int main()
{
   const int number_pT_bin = 15;
   double dNdpT[number_pT_bin][3];
   double Ae[number_pT_bin][3];

   double N_jpsi_gen[] = {187513, 1.69512e+06, 3.74128e+06, 3.0311e+06, 1.74076e+06, 872872, 418124, 203066, 102026, 53106, 29023, 16240, 9719, 11549, 4081};

   TString path_dNoverdpT = Form("/home/chunlu/cernbox/AnaSourceCode/anaJPsi_2017pp/JPsi_dNdpT_computation/dNdpT_pT_in_pp_2017.txt");
   TString path_dNgendpT = Form("/home/chunlu/cernbox/AnaSourceCode/anaJPsi_2017pp/JPsi_dNdpT_computation/Njpsi_gen_pt.txt");
   TString path_dNrecdpT = Form("/home/chunlu/cernbox/AnaSourceCode/anaJPsi_2017pp/JPsi_dNdpT_computation/Njpsi_rec_pt.txt");
//   TString path_to_Ae = Form("/home/chunlu/cernbox/AnaSourceCode/anaJPsi_2017pp/JPsi_accep_effi_computation/Ae_pT_in_pp_2017.txt");
   TString path_to_Ae = Form("/home/chunlu/cernbox/AnaSourceCode/anaJPsi_2017pp/JPsi_accep_effi_computation/correction_MC_input_shape/2nd_accep_effi_computation/Ae2_pt_in_pp_2017.txt");

   TString input_file_location = "/home/chunlu/local/anaedppSim2017/LHC18c5a_Geant3_jpsi_from_ttree_for_iteration_v4/2017ppMC_MULUResultsHistos.root";
   TString object_name = "DimuonHistosGenJpsi/GenJpsi_Pt";
   
   TFile *inputFile = new TFile(input_file_location.Data(),"READ");
   if ( inputFile->IsOpen() ) std::cout << Form("read %s successfully", input_file_location.Data()); std::cout << std::endl;
 
   float pt_bin[] = {0, 0.3, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 15, 20};
   TH1F* his_dNgendpt = new TH1F("his_dNgendpt",";p_{T}^{gen} GeV/c ;", number_pT_bin, pt_bin);
   his_dNgendpt->Sumw2();
   his_dNgendpt->SetStats(false);
   his_dNgendpt->SetMarkerStyle(4);
   his_dNgendpt->Add( (TH1F*)inputFile->Get( object_name.Data() ) ) ;

   for(int iPt=0; iPt < number_pT_bin; iPt++)
   {
      std::cout << his_dNgendpt->GetBinContent(iPt+1) << std::endl;
      N_jpsi_gen[iPt] = his_dNgendpt->GetBinContent(iPt+1);
   }
  
   const int ELEMENT_FOR_AE_PT_BIN = 16;

   TGraphAsymmErrors *gr_dNdpT = new TGraphAsymmErrors(ELEMENT_FOR_AE_PT_BIN );
//   gr_dsigma->SetName("gr_dsigma");
   TGraphAsymmErrors *gr_dNdpT_syst = new TGraphAsymmErrors((ELEMENT_FOR_AE_PT_BIN) );
//   gr_dsigma_syst->SetName("gr_dsigma_syst");
   TGraphAsymmErrors *gr_N_jpsi_gen = new TGraphAsymmErrors(  number_pT_bin-1 );

  double **variable;
  variable = new double *[number_pT_bin];
  for (int i=0; i<number_pT_bin; i++) variable[i] = new double[3];

  read_values_stat_sysm(path_dNoverdpT, number_pT_bin, variable);

  float Pt_bin_for_Ae[ELEMENT_FOR_AE_PT_BIN]= {0, 0.3, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 15, 20};

  TH1F* his_dNdpT = new TH1F(Form("his_dndpT"), " ", (ELEMENT_FOR_AE_PT_BIN-1), Pt_bin_for_Ae);
  his_dNdpT->Sumw2();

  TGraphAsymmErrors *gr_ratio = new TGraphAsymmErrors( number_pT_bin-1 );
  gr_ratio->SetName("gr_ratio");
  gr_ratio->SetMarkerStyle(20);

  for (int iPt=0; iPt<number_pT_bin; iPt++)
  {
      dNdpT[iPt][0] = variable[iPt][0];
      dNdpT[iPt][1] = variable[iPt][1];
      dNdpT[iPt][2] = variable[iPt][2];
  }

  double **variable_;
  variable_ = new double *[number_pT_bin];
  for (int i=0; i<number_pT_bin; i++) variable_[i] = new double[2];

//  read_values_stat(path_to_Ae, number_pT_bin, variable_);
  read_values_stat_sysm(path_to_Ae, number_pT_bin, variable);

  for (int iPt=0; iPt<number_pT_bin; iPt++)
  {
      Ae[iPt][0] = variable[iPt][0];
      Ae[iPt][1] = variable[iPt][1];
      Ae[iPt][2] = variable[iPt][2];            
  }

  double dNdpT_over_Ae[number_pT_bin]={0};
  double dNdpT_over_Ae_err=0;
  double dNdpT_over_Ae_syst_err=0;
  double ratio[number_pT_bin]={0};
  double sum_data_points1=0;
  double sum_data_points2=0;

  for (int iPt=0; iPt<number_pT_bin; iPt++) 
  {
    dNdpT_over_Ae[iPt] = dNdpT[iPt][0] / Ae[iPt][1];
    his_dNdpT->SetBinContent(iPt+1, dNdpT_over_Ae[iPt]);
    sum_data_points1 += dNdpT_over_Ae[iPt];

    dNdpT_over_Ae_err = error_propagation(1,  dNdpT[iPt][0],  dNdpT[iPt][1], Ae[iPt][1], Ae[iPt][2], dNdpT_over_Ae[iPt]);
    his_dNdpT->SetBinError(iPt+1, dNdpT_over_Ae_err);
    dNdpT_over_Ae_syst_err = error_propagation(1,  dNdpT[iPt][0],  dNdpT[iPt][2], Ae[iPt][1], 0, dNdpT_over_Ae[iPt]);
    std::cout << "dN/dPt over Ae: " << dNdpT_over_Ae[iPt] << " +/- " << dNdpT_over_Ae_err << " +/- " << dNdpT_over_Ae_syst_err << std::endl;

    double xVal=0, deltaxVal = 0, deltaxValSyst = 0;
    xVal = (Pt_bin_for_Ae[iPt+1] + Pt_bin_for_Ae[iPt])/2;
    deltaxVal = (Pt_bin_for_Ae[iPt+1] - Pt_bin_for_Ae[iPt])/2;
    deltaxValSyst = 0.3*deltaxVal;
    gr_dNdpT->SetPoint(iPt+1, xVal, dNdpT_over_Ae[iPt]);
    gr_dNdpT_syst->SetPoint(iPt+1, xVal, dNdpT_over_Ae[iPt]);

    gr_dNdpT->SetPointError(iPt+1, deltaxVal, deltaxVal, dNdpT_over_Ae_err, dNdpT_over_Ae_err);
    gr_dNdpT_syst->SetPointError(iPt+1, deltaxValSyst, deltaxValSyst, dNdpT_over_Ae_syst_err, dNdpT_over_Ae_syst_err);
 
    N_jpsi_gen[iPt] = N_jpsi_gen[iPt] / (Pt_bin_for_Ae[iPt+1] - Pt_bin_for_Ae[iPt]);
    sum_data_points2 += N_jpsi_gen[iPt];

//    ratio[number_pT_bin] = dNdpT_over_Ae[iPt] / N_jpsi_gen[iPt];
//    gr_N_jpsi_gen->SetPoint(iPt+1, xVal,  N_jpsi_gen[iPt]);

//    gr_ratio->SetPoint(iPt+1, xVal, ratio);
 
  }

  double norm = sum_data_points1 / sum_data_points2;  
  for (int iPt=0; iPt<number_pT_bin; iPt++)
  {    
    double xVal=0, deltaxVal = 0;
    xVal = (Pt_bin_for_Ae[iPt+1] + Pt_bin_for_Ae[iPt])/2;
    N_jpsi_gen[iPt] = N_jpsi_gen[iPt] * norm;
    ratio[iPt] = dNdpT_over_Ae[iPt] / N_jpsi_gen[iPt];
    gr_N_jpsi_gen->SetPoint(iPt, xVal,  N_jpsi_gen[iPt]);
    gr_ratio->SetPoint(iPt, xVal, ratio[iPt]);

    std::cout << " yeild ratio: " << ratio[iPt] << std::endl;
  }
  

//******************************
//
// Histogram makeup
//
//*****************************

  gr_dNdpT->SetLineWidth(1);
  gr_dNdpT->SetMarkerStyle(20);
  gr_dNdpT->SetMarkerSize(1);
  gr_dNdpT->SetMarkerColor(kRed);
  gr_dNdpT->SetLineColor(kRed);
  
  gr_dNdpT_syst->SetFillStyle(1000);
  gr_dNdpT_syst->SetFillColor(kRed-9);
  gr_dNdpT_syst->SetLineColor(kRed);

  gr_N_jpsi_gen->SetLineWidth(1);
  gr_N_jpsi_gen->SetMarkerStyle(20);
  gr_N_jpsi_gen->SetMarkerSize(1);
  gr_N_jpsi_gen->SetMarkerColor(kGray);
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

  TH2F *hratio = new TH2F("hratio","",100,0,20,100, 0.9, 1.3);
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
  gr_dNdpT_syst->Draw("2SAME");
  gr_dNdpT->Draw("PSAME");
  gr_N_jpsi_gen->Draw("LPSAME");

  TLegend *leg2 = new TLegend(0.4,0.7,0.8,0.85);
  leg2->SetBorderSize(0);
  //leg->SetTextSize(gStyle->GetTextSize()*0.8);
  leg2->SetLineColor(0);
  leg2->SetLineStyle(1);
  leg2->SetLineWidth(1);
  leg2->SetFillColor(0);
  leg2->SetFillStyle(1001);
  leg2->AddEntry(gr_dNdpT, "dN/dp_{T} * 1/A#epsilon","LPE");
  leg2->AddEntry(gr_dNdpT_syst,"Systematic uncertainty","F");
  leg2->AddEntry(gr_N_jpsi_gen, " input MC (arbitrary norm.)","LPE");
  leg2->Draw();

  pad = (TPad*) c2->GetPad(2);
  pad->cd();
  hratio->Draw();
  gr_ratio->Draw("CPSAME");
  TLine *line = new TLine(0,1,20,1);
  line->SetLineStyle(2);
  line->Draw("same");

  TLegend *leg3 = new TLegend(0.2,0.7,0.4,0.85);
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

  for (int i=0; i<number_pT_bin; i++) delete [] variable[i];
  delete [] variable;

  for (int i=0; i<number_pT_bin; i++) delete [] variable_[i];
  delete [] variable_;
  

return 0;
}


//____________________________________________________________

void read_values_stat_sysm(TString path_to_file, int number_pT_bin, double *output_variable[3])
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

//____________________________________________________________

void read_values_stat(TString path_to_file, int number_pT_bin, double *output_variable[2])
{
    std::ifstream input_file(path_to_file.Data() );
    if(input_file.is_open() )
    {
       std::cout << "Opening file: " << path_to_file.Data() << std::endl; 
       for(int iPt=0; iPt < number_pT_bin; iPt++)
       {
            input_file >> output_variable[iPt][0] >> output_variable[iPt][1]; 
            std::cout << "The values are read: " << output_variable[iPt][0] << " " << output_variable[iPt][1] << std::endl;
       }
    }
    else std::cout << "Unable to open file: " << path_to_file.Data() << std::endl; 
    input_file.close();
}


//____________________________________________________________
double error_propagation(int chose, double a, double delta_a, double b, double delta_b, double c)
{
  double delta_c=0;
  if(chose == 1)
  {  //for c= a/b, delta_c is c's uncertainty
     delta_c = c* sqrt(pow(delta_a/a,2)+pow(delta_b/b,2));
  }
  else if(chose == 2)
  {//for c= a-b
    delta_c = sqrt(pow(delta_a,2)+pow(delta_b,2));
  }

return delta_c;
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

