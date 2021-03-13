#include <iostream>
#include <fstream>

#include "TString.h"
#include "TCanvas.h"
#include "TH1.h"
#include "TGraphAsymmErrors.h"
#include "TH2F.h"
#include "TLegend.h"

void makeup_tgraph(TGraphAsymmErrors *, TGraphAsymmErrors *, int);
void drawing_histograms(TGraphAsymmErrors** , TGraphAsymmErrors**, int, TString*);

int main()
{
  const int ELEMENT_FOR_AE_PT_BIN =16;
  float pT_bin[ELEMENT_FOR_AE_PT_BIN]= {0,0.3,1,2,3,4,5,6,7,8,9,10,11,12,15,20};
  TString cent_range[3] = {"cent0to20", "cent20to40", "cent40to90"};
//  TString path_to_directory = Form("/home/chunlu/cernbox/AnaSourceCode/anaJPsi_2015_2018PbPb/fit_JPsi/fit_vs_pT/cent0to20/");
  TString path_to_directory;
// = Form("/home/chunlu/cernbox/AnaSourceCode/anaJPsi_2015_2018PbPb/fit_JPsi/fit_vs_pT/%s/", cent_range.Data());
  TString Njpsi_file_name = "jpsi_number.txt";
//  TString path_to_Njpsi;

  double Njpsi[3]={0};
  TH1F *dNdpT = new TH1F("","", ELEMENT_FOR_AE_PT_BIN-1, pT_bin);
  std::ofstream dNdpT_cent0to20_output_stream(Form("dNdpT_func_pT_PbPb_%s.txt", cent_range[0].Data()));
  std::ofstream dNdpT_cent20to40_output_stream(Form("dNdpT_func_pT_PbPb_%s.txt", cent_range[1].Data()));
  std::ofstream dNdpT_cent40to90_output_stream(Form("dNdpT_func_pT_PbPb_%s.txt", cent_range[2].Data()));

//  TGraphAsymmErrors *gr_dN_dpT = new TGraphAsymmErrors( ELEMENT_FOR_AE_PT_BIN );
//  TGraphAsymmErrors *gr_dN_dpT_syst = new TGraphAsymmErrors( ELEMENT_FOR_AE_PT_BIN );

  TGraphAsymmErrors *gr_dN_dpT[3];
  TGraphAsymmErrors *gr_dN_dpT_syst[3];

  int color_id = 2;   
  for(int iCent=0; iCent < 3; iCent++)
  {
    path_to_directory = Form("/home/chunlu/cernbox/AnaSourceCode/anaJPsi_2015_2018PbPb/fit_JPsi/fit_vs_pT/%s/", cent_range[iCent].Data());
    gr_dN_dpT[iCent] = new TGraphAsymmErrors( ELEMENT_FOR_AE_PT_BIN );
    gr_dN_dpT_syst[iCent]  = new TGraphAsymmErrors( ELEMENT_FOR_AE_PT_BIN );

    for(int ipT=0; ipT < (ELEMENT_FOR_AE_PT_BIN-1); ipT++)
    {
      TString pT_range;
      pT_range = Form("pT%gto%g/", pT_bin[ipT], pT_bin[ipT+1]);
      if(iCent == 0)
      {
        if(ipT == 2 || ipT == 3 || ipT==10 || ipT == 11 || ipT == 12 || ipT == 13 || ipT == 14)
        {
           pT_range = Form("pT%gto%g/reject_bad_cov_matrix/", pT_bin[ipT], pT_bin[ipT+1]);
        }
      }
//      else
//      {
//         pT_range = Form("pT%gto%g/", pT_bin[ipT], pT_bin[ipT+1]);
//      }

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

        gr_dN_dpT[iCent]->SetPoint(ipT+1, xVal, yVal );
        gr_dN_dpT[iCent]->SetPointError(ipT+1, xValerr, xValerr, yValerr, yValerr);

        gr_dN_dpT_syst[iCent]->SetPoint(ipT+1, xVal, yVal );
        gr_dN_dpT_syst[iCent]->SetPointError(ipT+1, xValSysterr, xValSysterr, yValSysterr, yValSysterr );

        std::cout << yVal << std::endl;
        if(iCent == 0 ) 
          dNdpT_cent0to20_output_stream << yVal << " " << gr_dN_dpT[iCent]->GetErrorY(ipT+1) << " " << gr_dN_dpT_syst[iCent]->GetErrorY(ipT+1) << std::endl; 
        else if(iCent == 1)
          dNdpT_cent20to40_output_stream << yVal << " " << gr_dN_dpT[iCent]->GetErrorY(ipT+1) << " " << gr_dN_dpT_syst[iCent]->GetErrorY(ipT+1) << std::endl; 
        else if(iCent == 2)
          dNdpT_cent40to90_output_stream << yVal << " " << gr_dN_dpT[iCent]->GetErrorY(ipT+1) << " " << gr_dN_dpT_syst[iCent]->GetErrorY(ipT+1) << std::endl;           
      }
      else
      {
        std::cout << "Unable to open file: " << path_to_Njpsi.Data() << std::endl;
      }
    
    }     
    makeup_tgraph(gr_dN_dpT[iCent], gr_dN_dpT_syst[iCent], color_id+iCent);
  }

  drawing_histograms(gr_dN_dpT, gr_dN_dpT_syst, 3, cent_range);
 
  delete dNdpT;

return 0;
}

void makeup_tgraph(TGraphAsymmErrors *gr, TGraphAsymmErrors *gr_syst, int color_id)
{
//    gr->SetLineColor(2);
//  gr->SetLineWidth(2);

//    gr_syst->SetFillStyle(1000);
//    gr_syst->SetFillColor(46);
//    gr_syst->SetLineWidth(2);
//    gr_syst->SetLineStyle(1);

    gr->SetLineWidth(1);
    gr->SetMarkerStyle(20);
    gr->SetMarkerSize(0.5);
//    gr->SetMarkerColor(kRed);
//    gr->SetLineColor(kRed);
    gr->SetMarkerColor(color_id);
    gr->SetLineColor(color_id);    

    gr_syst->SetFillStyle(1000);
    gr_syst->SetFillColorAlpha(color_id, 0.2);
    gr_syst->SetLineColor(color_id);

//    gr_syst->SetFillColor(kRed-9);
//    gr_syst->SetLineColor(kRed);

  //gr_dsigma_syst->SetLineStyle(1); 
  //    gr_dsigma_syst->SetLineWidth(2);
}

void drawing_histograms(TGraphAsymmErrors **gr_dN_dpT, TGraphAsymmErrors **gr_dN_dpT_syst, int num_cent, TString* cent_range)
{
  TH2F *histo = new TH2F("histo","",100,0,20,100,9.,5.e5);
  histo->GetYaxis()->SetLabelSize(0.04);
  histo->GetYaxis()->SetTitleSize(0.04);
  histo->GetYaxis()->SetTitleOffset(1.15);
  histo->SetStats(0);
//  histo->SetTitle(Form("%s", cent_range.Data()));
  histo->GetYaxis()->SetTitle("dN/d#it{p}_{T} (/(GeV/#it{c}))");
  histo->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");


  TCanvas *tc_jpsi_pt = new TCanvas("c0","c0",600,500);
  tc_jpsi_pt->SetLogy();

  TLegend *leg2 = new TLegend(0.4,0.7,0.8,0.85);
  leg2->SetBorderSize(0);
  //leg->SetTextSize(gStyle->GetTextSize()*0.8);
  leg2->SetLineColor(0);
  leg2->SetLineStyle(1);
  leg2->SetLineWidth(1);
  leg2->SetFillColor(0);
  leg2->SetFillStyle(1001);
//  leg2->AddEntry(gr_dsigma, "ALICE Inclusive J/#psi -> #mu^{+}#mu^{-}, 2.5 < y < 4, Lint = 1.23 pb","LPE");

  histo->Draw();
  for(int iCent = 0; iCent < num_cent; iCent++)
  {
    gr_dN_dpT_syst[iCent]->Draw("2SAME");
    gr_dN_dpT[iCent]->Draw("PSAME");
    leg2->AddEntry(gr_dN_dpT[iCent], Form("%s: stat. ", cent_range[iCent].Data() ),"LPE");
    leg2->AddEntry(gr_dN_dpT_syst[iCent], Form("syst. " ),"F");  
  }
//  gr_dN_dpT_syst->Draw("2SAME");
//  gr_dN_dpT->Draw("PSAME");
  leg2->Draw();


//  TString output_plot_name = "dNdpT_func_pT_PbPb_cent20to40.pdf";
  TString output_plot_name = Form("dNdpT_func_pT_all_cent_PbPb.pdf");
  tc_jpsi_pt->Print( output_plot_name.Data() );


  delete histo;
  delete tc_jpsi_pt;
}

