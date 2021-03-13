#include <iostream>
//#include "PtCentRangeValues.h"
#include "TString.h"
#include "TH1F.h"
#include "TFile.h"
#include "read_access_histos_y_dependence.h"

void read_access_histos_y_dependence(TFile *inputFile, int init_y, int init_Cent, int end_Cent, TH1F **histo_y, TH1F *histo_integrated_y)
{

   /*
   * Pt bins initialize
   */
//  double arrayAvailablePtBins[]={0.0,0.3,1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0,10.0,11.0,12.0,13.0,14.0,15.0};
//  double arrayAvailablePtBins[]={0.0,1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0,10.0,11.0,12.0};
//  const int numberOfPtBins = 13;

//  double arrayAvailablePtBins_photon[]={0,0.3,1,2,3,4,5,6,7,8,9,10,11,12,15,20};
//  const int numberOfPtBins_photon = (int)(sizeof(arrayAvailablePtBins_photon)/sizeof(arrayAvailablePtBins_photon[0]));

//  double arrayPtBinsPhotonProduction[]={0,0.3,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25};
//  const int numberOfPtBinsPhotonProduction = (int)(sizeof(arrayPtBinsPhotonProduction)/sizeof(arrayPtBinsPhotonProduction[0]));

//  double arrayPtBinsPhotonProduction[]={0,0.3,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20};
//  const int numberOfPtBinsPhotonProduction = (int)(sizeof(arrayPtBinsPhotonProduction)/sizeof(arrayPtBinsPhotonProduction[0]));

    double arrayAvailableYBins[] = {2.5, 2.75, 3.0, 3.25, 3.5, 3.75, 4.0};
    const int numberOfYBins = 6;
  
//  const int numberOfPtBinsPhotonProduction=22;


  /*
   * Centrality bins initialize
   */
  int CentMin = 0;
  int CentMax = 90;
  int CentInterval = 10;
  int nCentIndex =1 + (CentMax - CentMin) / CentInterval;
  const int NUMOFCENT_BINS = (CentMax - CentMin) / CentInterval;
  int *arrayAvailableCentBins = new int[nCentIndex];

  arrayAvailableCentBins[0] = CentMin;
  for(int i=1; i<nCentIndex; i++)  arrayAvailableCentBins[i] = arrayAvailableCentBins[i-1] + CentInterval;

//    for(int iCent = 0; iCent < (nCentIndex-1); iCent++)
    for(int iCent = init_Cent; iCent < end_Cent; iCent++)
    {
      TString strCentBins;
      strCentBins.Form("Cent%ito%i", arrayAvailableCentBins[iCent], arrayAvailableCentBins[iCent+1]);
      std::cout << strCentBins.Data() << std::endl;

      for(int iY = init_y; iY < (numberOfYBins); iY++)
      { 
        TString strYBins;
        strYBins.Form("Y%gto%g", arrayAvailableYBins[iY], arrayAvailableYBins[iY+1]);
        std::cout << strYBins.Data() << std::endl;

        histo_integrated_y->Add((TH1F*) inputFile->Get( Form("DimuonHistosUnlike/InvMass/histoInvMass_%s_%s", strCentBins.Data(), strYBins.Data() ) ) );
//        histoRawInvM->Add((TH1F*) inputFile_MergedCMUL_->Get( Form("DimuonHistosUnlike/InvMass/histoInvMass_%s",strPtBins.Data() ) ) ); 
        histo_y[iY]->Add((TH1F*) inputFile->Get( Form("DimuonHistosUnlike/InvMass/histoInvMass_%s_%s", strCentBins.Data(), strYBins.Data() ) ) );
//        histoInvMass[iPt]->Add((TH1F*) inputFile_MergedCMUL->Get( Form("DimuonHistosUnlike/InvMass/histoInvMass_%s",strPtBins.Data() ) ) ); 

      }
/*
      histo_pt[numberOfPtBins_photon-3]->Add( (TH1F*) inputFile->Get( Form("DimuonHistosUnlike/InvMass/histoInvMass_%s_Pt12to13", strCentBins.Data() ) ) );
      histo_pt[numberOfPtBins_photon-3]->Add( (TH1F*) inputFile->Get( Form("DimuonHistosUnlike/InvMass/histoInvMass_%s_Pt13to14", strCentBins.Data() ) ) );
      histo_pt[numberOfPtBins_photon-3]->Add( (TH1F*) inputFile->Get( Form("DimuonHistosUnlike/InvMass/histoInvMass_%s_Pt14to15", strCentBins.Data() ) ) );

      histo_pt[numberOfPtBins_photon-2]->Add( (TH1F*) inputFile->Get( Form("DimuonHistosUnlike/InvMass/histoInvMass_%s_Pt15to16", strCentBins.Data() ) ) );
      histo_pt[numberOfPtBins_photon-2]->Add( (TH1F*) inputFile->Get( Form("DimuonHistosUnlike/InvMass/histoInvMass_%s_Pt16to17", strCentBins.Data() ) ) );
      histo_pt[numberOfPtBins_photon-2]->Add( (TH1F*) inputFile->Get( Form("DimuonHistosUnlike/InvMass/histoInvMass_%s_Pt17to18", strCentBins.Data() ) ) );
      histo_pt[numberOfPtBins_photon-2]->Add( (TH1F*) inputFile->Get( Form("DimuonHistosUnlike/InvMass/histoInvMass_%s_Pt18to19", strCentBins.Data() ) ) );
      histo_pt[numberOfPtBins_photon-2]->Add( (TH1F*) inputFile->Get( Form("DimuonHistosUnlike/InvMass/histoInvMass_%s_Pt19to20", strCentBins.Data() ) ) );

      histo_integrated_pt->Add( (TH1F*) inputFile->Get( Form("DimuonHistosUnlike/InvMass/histoInvMass_%s_Pt12to13", strCentBins.Data() ) ) );
      histo_integrated_pt->Add( (TH1F*) inputFile->Get( Form("DimuonHistosUnlike/InvMass/histoInvMass_%s_Pt13to14", strCentBins.Data() ) ) );
      histo_integrated_pt->Add( (TH1F*) inputFile->Get( Form("DimuonHistosUnlike/InvMass/histoInvMass_%s_Pt14to15", strCentBins.Data() ) ) );

      histo_integrated_pt->Add( (TH1F*) inputFile->Get( Form("DimuonHistosUnlike/InvMass/histoInvMass_%s_Pt15to16", strCentBins.Data() ) ) );
      histo_integrated_pt->Add( (TH1F*) inputFile->Get( Form("DimuonHistosUnlike/InvMass/histoInvMass_%s_Pt16to17", strCentBins.Data() ) ) );
      histo_integrated_pt->Add( (TH1F*) inputFile->Get( Form("DimuonHistosUnlike/InvMass/histoInvMass_%s_Pt17to18", strCentBins.Data() ) ) );
      histo_integrated_pt->Add( (TH1F*) inputFile->Get( Form("DimuonHistosUnlike/InvMass/histoInvMass_%s_Pt18to19", strCentBins.Data() ) ) );
      histo_integrated_pt->Add( (TH1F*) inputFile->Get( Form("DimuonHistosUnlike/InvMass/histoInvMass_%s_Pt19to20", strCentBins.Data() ) ) );
*/
    }


}
