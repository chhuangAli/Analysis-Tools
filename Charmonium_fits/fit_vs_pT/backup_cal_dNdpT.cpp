#include <iostream>
#include <fstream>

#include "TString.h"
#include "TCanvas.h"
#include "TH1.h"
#include "TGraphAsymmErrors.h"
#include "TH2F.h"
#include "TLegend.h"

void makeup_tgraph(TGraphAsymmErrors *, TGraphAsymmErrors *);
void drawing_histograms(TGraphAsymmErrors* , TGraphAsymmErrors*, TString );

int main()
{
  const int ELEMENT_FOR_AE_PT_BIN =16;
  float pT_bin[ELEMENT_FOR_AE_PT_BIN]= {0,0.3,1,2,3,4,5,6,7,8,9,10,11,12,15,20};
  TString cent_range = "cent40to90";
//  TString path_to_directory = Form("/home/chunlu/cernbox/AnaSourceCode/anaJPsi_2015_2018PbPb/fit_JPsi/fit_vs_pT/cent0to20/");
  TString path_to_directory = Form("/home/chunlu/cernbox/AnaSourceCode/anaJPsi_2015_2018PbPb/fit_JPsi/fit_vs_pT/%s/", cent_range.Data());
  TString Njpsi_file_name = "jpsi_number.txt";
//  TString path_to_Njpsi;

  double Njpsi[3]={0};
  TH1F *dNdpT = new TH1F("","", ELEMENT_FOR_AE_PT_BIN-1, pT_bin);
  std::ofstream dNdpT_output_stream(Form("dNdpT_func_pT_PbPb_%s.txt", cent_range.Data()));

  TGraphAsymmErrors *gr_dN_dpT = new TGraphAsymmErrors( ELEMENT_FOR_AE_PT_BIN );
  TGraphAsymmErrors *gr_dN_dpT_syst = new TGraphAsymmErrors( ELEMENT_FOR_AE_PT_BIN );

  for(int ipT=0; ipT < (ELEMENT_FOR_AE_PT_BIN-1); ipT++)
  {
    TString pT_range;
//    if(ipT == 2 || ipT == 3 || ipT==10 || ipT == 11 || ipT == 12 || ipT == 13 || ipT == 14)
//    {
//       pT_range = Form("pT%gto%g/reject_bad_cov_matrix/", pT_bin[ipT], pT_bin[ipT+1]);
//    }
//    else
//    {
       pT_range = Form("pT%gto%g/", pT_bin[ipT], pT_bin[ipT+1]);
//    }
    TString path_to_Njpsi = Form("%s%s%s", path_to_directory.Data(), pT_range.Data(), Njpsi_file_name.Data() );

    std::ifstream input_JPsi_number(path_to_Njpsi.Data() );
    if (input_JPsi_number.is_open())
    {
      input_JPsi_number >> Njpsi[0] >> Njpsi[1] >> Njpsi[2];
      input_JPsi_number.close();

//      dNdpT_output_stream << Njpsi[0] << " " << Njpsi[1] << " " << Njpsi[2] << std::endl;
      std::cout << Njpsi[0] << " " << Njpsi[1] << " " << Njpsi[2] << std::endl;

      double xVal =0;
      double xValerr = 0, xValSysterr=0;
      double yVal=0;
      double yValerr =0, yValSysterr =0;

      xVal = (pT_bin[ipT+1] + pT_bin[ipT]) / 2;
      xValerr = (pT_bin[ipT+1] - pT_bin[ipT]) / 2.;
      xValSysterr = 0.3*xValerr;

      yVal = Njpsi[0] / (pT_bin[ipT+1]-pT_bin[ipT]);
      yValerr = yVal * Njpsi[1] / Njpsi[0];

      yValSysterr =  yVal * Njpsi[2] / Njpsi[0];

      gr_dN_dpT->SetPoint(ipT+1, xVal, yVal );
      gr_dN_dpT->SetPointError(ipT+1, xValerr, xValerr, yValerr, yValerr);


      gr_dN_dpT_syst->SetPoint(ipT+1, xVal, yVal );
      gr_dN_dpT_syst->SetPointError(ipT+1, xValSysterr, xValSysterr, yValSysterr, yValSysterr );

      std::cout << yVal << std::endl;
      dNdpT_output_stream << yVal << " " << gr_dN_dpT->GetErrorY(ipT+1) << " " << gr_dN_dpT_syst->GetErrorY(ipT+1) << std::endl;

    }
    else
    {
      std::cout << "Unable to open file: " << path_to_Njpsi.Data() << std::endl;
    }
    
  }
 
  makeup_tgraph(gr_dN_dpT, gr_dN_dpT_syst);
  drawing_histograms(gr_dN_dpT, gr_dN_dpT_syst, cent_range);
 
  delete dNdpT;

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

void drawing_histograms(TGraphAsymmErrors *gr_dN_dpT, TGraphAsymmErrors *gr_dN_dpT_syst, TString cent_range)
{
  TH2F *histo = new TH2F("histo","",100,0,20,100,9.,5.e5);
  histo->GetYaxis()->SetLabelSize(0.04);
  histo->GetYaxis()->SetTitleSize(0.04);
  histo->GetYaxis()->SetTitleOffset(1.15);
  histo->SetStats(0);
  histo->SetTitle(Form("%s", cent_range.Data()));
  histo->GetYaxis()->SetTitle("dN/d#it{p}_{T} (/(GeV/#it{c}))");
  histo->GetXaxis()->SetTitle("#it{p}_{T} ((GeV/#it{c}))");

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
  leg2->AddEntry(gr_dN_dpT_syst,"Systematic uncertainty","F");
  leg2->Draw();

//  TString output_plot_name = "dNdpT_func_pT_PbPb_cent20to40.pdf";
  TString output_plot_name = Form("dNdpT_func_pT_PbPb_%s.pdf", cent_range.Data());
  tc_jpsi_pt->Print( output_plot_name.Data() );


  delete histo;
  delete tc_jpsi_pt;
}

