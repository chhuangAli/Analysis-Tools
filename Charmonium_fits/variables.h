#include <iostream>
#include <fstream>
#include "TH1.h"
#include "TAxis.h"
#include "TString.h"

//      TString path_to_histos_directory = "$HOME/local/anaedPbPbData2018/LHC18q_r/muon_calo_pass3/AOD225_CMULEvent_AnaResults_Pcorr/";
//      TString path_to_histos_directory = "$HOME/local/anaedPbPbData2018/LHC18q_r/muon_calo_pass3/AOD225_CMULEvent_AnaResults/";
//      TString histos_name = "2018PbPb_CMULResultsHistos_LHC18q.root";

//      TString path_to_histos_directory = "$HOME/local/anaedPbPbData2018/LHC18q_r/muon_calo_pass3/AOD216_CMULEvent_AnaResults_Pcorr/";
//      TString path_to_histos_directory = "$HOME/local/anaedPbPbData2018/LHC18q_r/muon_calo_pass3/AOD216_CMULEvent_AnaResults/";
      TString path_to_histos_directory = "$HOME/local/anaedPbPbData2018/LHC18q_r/muon_calo_pass3/AOD225_CMULEvent_AnaResults/";
      TString histos_name = "2018PbPb_CMULResultsHistos_LHC18r.root";

//     TString path_to_histos_directory = "$HOME/local/anaedPbPbData2015/LHC15o/muon_calo_pass1/AOD197_CMULEvent_AnaResults_Pcorr/";
//     TString path_to_histos_directory = "$HOME/local/anaedPbPbData2015/LHC15o/muon_calo_pass1/AOD197_CMULEvent_AnaResults/";
//     TString histos_name = "2015PbPb_CMULResultsHistos_LHC15o.root";

      TString path_to_histos = Form("%s%s",path_to_histos_directory.Data(), histos_name.Data() );
//      TString path_to_histos_corr="$HOME/local/anaedPbPbData2018/LHC18q_r/muon_calo_pass3/AOD225_CMULEvent_AnaResults_Pcorr/2018PbPb_CMULResultsHistos_LHC18q.root";
//      TString path_to_histos="$HOME/local/anaedPbPbData2018/LHC18q_r/muon_calo_pass3/AOD225_CMULEvent_AnaResults/2018PbPb_CMULResultsHistos_LHC18q.root";
  
      const int NUMOFBINS=400;
      const int numberOfDimuonCharge =3;
      TString dimuonChargeNames[numberOfDimuonCharge] = {"MMLike","Unlike","PPLike"};

      float lowerLimit_histoInvMass = 0;
      float upperLimit_histoInvMass = 10;
//  const int NUMOFPTBIN_UNDER12 = 13;

      int nbin;
      double _xmin;
      double _xmax;
      double binWidth;

      const int numOfCB2_Parm= 7;
//  const int numOfNA60_Parm= 11;
      const int numOfVWG_Parm= 5;
//  const int numOfPol1overPol2_Parm= 5;
//  const int numOfPol2overPol3_Parm= 7;
      const int numOfExp_Parm = 4;

//  const int FITRANGE = 2;

      const int NUMOFTAILPARAM_CB2=2;
      const int NUMOFTAILPARAM_NA60=1;

//      const int NUMOFBACKGROUNDFUNC=2;
      const int NUMOFBACKGROUNDFUNC=3; // including the event mixing fitting

      // 12 -> 11<pT<12
      // 13 -> 12<pT<15
      // 14 -> 15<pT<20
      // 15 -> 0<pT<20
      int pt_bin = 15; // indicate which pT bin to be used. It is only effect when fit the signal as a function of pT. 
      int y_bin = 0;
      
      TString tstr_pt_bin = Form("0 < #it{p_{T}} < 20"); //for tlatex on mass plot
      TString tstr_cent_bin = Form("centrality 0-10%%");
      TString tstr_system = Form("Pb-Pb, #sqrt{#it{s}_{NN}} = 5.02 TeV");
      TString tstr_y_bin = Form("2.5 < #it{y} < 4");
      

//      double display_NJPsi_range = 2000;

      std::ofstream output_text_file; // for output the number of jpsi with statistical and systematic uncertainty

      bool fit_valide_CB2VWG[FITRANGE][NUMOFTAILPARAM_CB2][1];
      int fit_status_CB2VWG[FITRANGE][NUMOFTAILPARAM_CB2][1];
      int fit_cov_matrix_status_CB2VWG[FITRANGE][NUMOFTAILPARAM_CB2][1];

      bool fit_valide_CB2Pol[FITRANGE][NUMOFTAILPARAM_CB2][1];
      int fit_status_CB2Pol[FITRANGE][NUMOFTAILPARAM_CB2][1];
      int fit_cov_matrix_status_CB2Pol[FITRANGE][NUMOFTAILPARAM_CB2][1];

      bool fit_valide_CB2Exp[FITRANGE][NUMOFTAILPARAM_CB2][1];
      int fit_status_CB2Exp[FITRANGE][NUMOFTAILPARAM_CB2][1];
      int fit_cov_matrix_status_CB2Exp[FITRANGE][NUMOFTAILPARAM_CB2][1];

      bool fit_valide_NA60VWG[FITRANGE][NUMOFTAILPARAM_NA60][1];
      int fit_status_NA60VWG[FITRANGE][NUMOFTAILPARAM_NA60][1];
      int fit_cov_matrix_status_NA60VWG[FITRANGE][NUMOFTAILPARAM_NA60][1];

      bool fit_valide_NA60Pol[FITRANGE][NUMOFTAILPARAM_NA60][1];
      int fit_status_NA60Pol[FITRANGE][NUMOFTAILPARAM_NA60][1];
      int fit_cov_matrix_status_NA60Pol[FITRANGE][NUMOFTAILPARAM_NA60][1];

      bool fit_valide_NA60Exp[FITRANGE][NUMOFTAILPARAM_NA60][1];
      int fit_status_NA60Exp[FITRANGE][NUMOFTAILPARAM_NA60][1];
      int fit_cov_matrix_status_NA60Exp[FITRANGE][NUMOFTAILPARAM_NA60][1];

      const int min=0;
  //max is equal to 18 (6:fit after mixed-event method)
      const int max = FITRANGE * NUMOFBACKGROUNDFUNC * (NUMOFTAILPARAM_CB2 + NUMOFTAILPARAM_NA60);


      TH1D *histo_NJpsifitUncertainty = new TH1D("NJpsifitUncerntainty",";;N_{J/#psi}",(max+1),min,(max+1));
      TH1D *histo_SigmaJpsifitUncertainty = new TH1D("SigmaJpsifitUncerntainty",";;#sigma_{J/#psi}",(max+1),min,(max+1));
      TH1D *histo_MassJpsifitUncertainty = new TH1D("MassJpsifitUncerntainty",";;M_{J/#psi}",(max+1),min,(max+1));
      TH1D *histo_ChiJpsifitUncertainty = new TH1D("ChiJpsifitUncerntainty",";;#chi^{2}/ndf",(max+1),min,(max+1));

      TAxis *xaxis = histo_ChiJpsifitUncertainty->GetXaxis();
      //  TAxis *xaxis = histo_NJpsifitUncertainty->GetXaxis();

      int iCount=1;
      double total_numberOfJPsi=0;
      double total_error_numberOfJPsi=0;
      double total_sigmaOfJPsi=0;
      double total_massOfJPsi=0;

      TH1I *histo_fit_valide = new TH1I("fit_valide",";;fit valid",(max+1),min,(max+1));
      TAxis *xaxis_fit_valide = histo_fit_valide->GetXaxis();

      TH1I *histo_fit_status = new TH1I("fit_status",";;fit status",(max+1),min,(max+1));
      TAxis *xaxis_fit_status = histo_fit_status->GetXaxis();

      TH1I *histo_cov_matrix_status = new TH1I("fit_cov_matrix_status",";;fit cov. matrix status",(max+1),min,(max+1));
      TAxis *xaxis_cov_matrix_status = histo_cov_matrix_status->GetXaxis();

      TString save_fit_quality_canvas_name = Form("fit_quality.pdf");


