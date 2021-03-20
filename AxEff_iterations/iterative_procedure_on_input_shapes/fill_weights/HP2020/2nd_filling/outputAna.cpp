#include <iostream>
#include <fstream>
#include <string>

#include "TCanvas.h"
#include "TH1F.h"
#include "TH2.h"
#include "TH1D.h"
#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"
#include "TObjArray.h"
#include "TLorentzVector.h"
#include "TFormula.h"

#include "AliCounterCollection.h"

#ifdef __MAKECINT__

#pragma link C++ class std::vector<double>

#endif

//#include "AliCounterCollection.h"
//The definition of aliCounterCollection is in -I${ALICE_INCLUDE}/include -L${ALICE_ROOT}/lib -lSTEERBase

using std::cout; using std::endl; using std::string;

//const int numberOfDimuonCharges = 3;
const int TypeOfJPsi = 2;
enum enumEventDimuonCharge {kMMLike,kUnlike,kPPLike};
//TString dimuonChargeNames[numberOfDimuonCharges] = {"MMLike","Unlike","PPLike"};
TString dimuonChargeNames[TypeOfJPsi] = {"GenJpsi","RecJpsi"};

TString strCentralityBins = "0,10,20,30,40,50,60,70,80,90,110";
//TString strCentralityBins = "0,20,40,90";
TString strPtBins = "0,0.3,1,2,3,4,5,6,7,8,9,10,11,12,15,20";
//TString strPtBins = "0,0.3,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20";
TString strRapidityBins ="2.5,2.75,3,3.25,3.5,3.75,4";
TString strPtBins_Double_Diff = "0.3,2,4,6,12,20";
TString strYBins_Double_Diff = "2.5,3,3.5,4";
//TString strPtBins = "0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25";

std::vector<double> StringToVector(TString);
int FindBinInVector( double , std::vector<double> vectorBins );
TString FindCentBin(double lPercentile);
void AnaAndMakeHistos(int );
double setWeight_Pt(double, double *, double *);
double setWeight_Y(double, double *, double *);

int main(int argc, char *argv[])
{
  string input;
  input = argv[1]; 
  int runNumber = std::stoi(input);
  AnaAndMakeHistos(runNumber);

return 0;
}


void AnaAndMakeHistos(int RunNumber)
{
  std::vector<double> vectorCentBins = StringToVector(strCentralityBins);
  int numberOfCentBins = vectorCentBins.size()-1;

  std::vector<double> vectorPtBins = StringToVector(strPtBins);
  int numberOfPtBins = vectorPtBins.size()-1;

  std::vector<double> vectorRapidityBins = StringToVector(strRapidityBins);
  int numberOfRapidityBins = vectorRapidityBins.size()-1;

  std::vector<double> vectorPtBins_Double_Diff = StringToVector(strPtBins_Double_Diff);
  int numberOfPtBins_Double_Diff = vectorPtBins_Double_Diff.size()-1;

  std::vector<double> vectorYBins_Double_Diff = StringToVector(strYBins_Double_Diff);
  int numberOfYBins_Double_Diff = vectorYBins_Double_Diff.size()-1;

  float y_bin_double_diff[] = {-4, -3.5, -3., -2.5};

  int number_pt_bin = 15;
  int number_y_bin = 6;
  float pt_bin[] = {0, 0.3, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 15, 20};
//  float y_bin[]= {2.5, 2.75, 3.0, 3.25, 3.5, 3.75, 4};
  float y_bin[] = {-4., -3.75, -3.5, -3.25, -3., -2.75, -2.5};
  //RunNumber = 295581;
  //const int RunNumber = 245145;

//  TString ReadFileLocation = Form("$HOME/local/anaedPbPbSim2018/LHC19a2_jpsi_ttree/CINTEvent_Results_cuts/AnalysisResults_run%d.root", RunNumber);
//  TString ReadFileLocation = Form("$HOME/local/anaedPbPbSimData2016/LHC16e2_plus_jpsi_ttree/CINTEvent_Results_cuts/AnalysisResults_run%d.root", RunNumber);
  TString ReadFileLocation = Form("$HOME/local/anaedPbPbSimData2016/LHC16e2_jpsi_ttree/CINTEvent_Results_cuts/AnalysisResults_run%d.root", RunNumber);

  //TString ReadFileLocation = Form("../anaPbPbTask.root");
//  TString ReadFileLocation = Form("$HOME/local/anaedppSim2017/LHC18c5a_Geant3_jpsi_ttree_v2/AnalysisResults_run%d.root", RunNumber);   
  //TString ReadFileLocation = Form("/afs/cern.ch/work/c/chhuang/private/2018_TaskPbPb_Results/LHC18r/muon_calo_pass2/AOD211_CMUL/AnalysisResults_run000%d.root",RunNumber);
  //TString ReadFileLocation = Form("$HOME/private/GRID_Task/analysis_PbPb_2018/GRID_analysis_task_CMUL/AnalysisResults.root");

//  TFile *outputFile = new TFile(Form("/home/chunlu/local/anaedPbPbSim2018/LHC19a2_jpsi_ttree/CINTEvent_Results_cuts/CINTEvent_MCinput_2weights_AnaResults/PbPbMC_ResultsHistos_Run%d.root",RunNumber),"recreate");
//  TFile *outputFile = new TFile(Form("/home/chunlu/local/anaedPbPbSimData2016/LHC16e2_plus_jpsi_ttree/CINTEvent_Results_cuts/CINTEvent_MCinput_2weights_AnaResults/PbPbMC_ResultsHistos_Run%d.root",RunNumber),"recreate");
  TFile *outputFile = new TFile(Form("/home/chunlu/local/anaedPbPbSimData2016/LHC16e2_jpsi_ttree/CINTEvent_Results_cuts/CINTEvent_MCinput_2weights_AnaResults/PbPbMC_ResultsHistos_Run%d.root",RunNumber),"recreate");

  AliCounterCollection *fEventCounters;
  fEventCounters = new AliCounterCollection("JPsiCounters");
  fEventCounters->AddRubric("JPsi","gen/rec");
  fEventCounters->AddRubric("run",1000000);
  fEventCounters->AddRubric("centrality", "m0/0_10/10_20/20_30/30_40/40_50/50_60/60_70/70_80/80_90/90_100");
//  fEventCounters->AddRubric("selected","yes/no");
  fEventCounters->Init();

  TDirectory *directoryDimuonCounter;
  TDirectory *directoryDimuon[TypeOfJPsi];
//  TDirectory *directoryMCinput[TypeOfJPsi];
  TDirectory *directoryInvMass[TypeOfJPsi];
  TDirectory *directoryPt[TypeOfJPsi];
  TDirectory *directoryRapidity[TypeOfJPsi];

  directoryDimuonCounter = outputFile -> mkdir( Form("DimuonCounter") );
  directoryDimuonCounter->Add(fEventCounters);
  for( int iType=0; iType < TypeOfJPsi; iType++ )
  {
    directoryDimuon[iType] = outputFile -> mkdir( Form("DimuonHistos%s", dimuonChargeNames[iType].Data()) );
    directoryInvMass[iType] = directoryDimuon[iType] -> mkdir( "InvMass" );
    directoryPt[iType] = directoryDimuon[iType] -> mkdir( "Pt" );
    directoryRapidity[iType] = directoryDimuon[iType] -> mkdir( "Rapidity" );
  }


  TH1F *histoWeightGen_20to30 = new TH1F("WeightGen_20to30","", 200,-5.,5. );
  TH1F *histoWeightGenPt_20to30 = new TH1F("WeightGenPt_20to30","", 200,-5.,5. );
  TH1F *histoWeightGenY_20to30 = new TH1F("WeightGenY_20to30","", 200,-5.,5. );

  TH1F *histoWeightGen_30to40 = new TH1F("WeightGen_30to40","", 200,-5.,5. );
  TH1F *histoWeightGenPt_30to40 = new TH1F("WeightGenPt_30to40","", 200,-5.,5. );
  TH1F *histoWeightGenY_30to40 = new TH1F("WeightGenY_30to40","", 200,-5.,5. );

  TH1F *histoWeightGen_40to90 = new TH1F("WeightGen_40to90","", 200,-5.,5. );
  TH1F *histoWeightGenPt_40to90 = new TH1F("WeightGenPt_40to90","", 200,-5.,5. );
  TH1F *histoWeightGenY_40to90 = new TH1F("WeightGenY_40to90","", 200,-5.,5. );


  directoryDimuon[0]->cd();  
//  TH1F *histoGenJpsi_Y = new TH1F("GenJpsi_Y","", number_y_bin, y_bin);
//  histoGenJpsi_Y->Sumw2();

  for(int iCent=0; iCent < numberOfCentBins; iCent++)
  {
    TH1F *histoGenJpsi_Pt = new TH1F(Form("GenJpsi_Pt_Cent%gto%g", vectorCentBins[iCent], vectorCentBins[iCent+1]),"", number_pt_bin, pt_bin);
    histoGenJpsi_Pt->Sumw2();
    TH1F *histoGenJpsi_Y = new TH1F(Form("GenJpsi_Y_Cent%gto%g", vectorCentBins[iCent], vectorCentBins[iCent+1]),"", number_y_bin, y_bin);
    histoGenJpsi_Y->Sumw2();
  }


//  directoryDimuon[0]->Add(fEventCounters);

  for(int iCent=0; iCent < numberOfCentBins; iCent++)
  {
    
    for(int iPt=0; iPt<numberOfPtBins_Double_Diff; iPt++)
    {
      TH1F *histoGenJpsi_YPt = new TH1F(Form("GenJpsi_Y_Cent%gto%g_Pt%gto%g", vectorCentBins[iCent], vectorCentBins[iCent+1], vectorPtBins_Double_Diff[iPt], vectorPtBins_Double_Diff[iPt+1]),"", 60,-5,-2);
      histoGenJpsi_YPt->Sumw2();
    }
  }


  directoryDimuon[1]->cd();
//  TH1F *histoRecJpsi_Pt = new TH1F("RecJpsi_Pt","", number_pt_bin, pt_bin);
//  histoRecJpsi_Pt->Sumw2();

  for(int iCent=0; iCent < numberOfCentBins; iCent++)
  {
    TH1F *histoRecJpsi_Pt = new TH1F(Form("RecJpsi_Pt_Cent%gto%g", vectorCentBins[iCent], vectorCentBins[iCent+1]),"", number_pt_bin, pt_bin);
    histoRecJpsi_Pt->Sumw2();
    TH1F *histoRecJpsi_Y = new TH1F(Form("RecJpsi_Y_Cent%gto%g", vectorCentBins[iCent], vectorCentBins[iCent+1]),"", number_y_bin, y_bin);
    histoRecJpsi_Y->Sumw2();
  }

  for(int iCent=0; iCent < numberOfCentBins; iCent++)  
  {
    for(int iPt=0; iPt< numberOfPtBins_Double_Diff; iPt++)
    {
      TH1F *histoRecJpsi_YPt = new TH1F(Form("RecJpsi_Y_Cent%gto%g_Pt%gto%g", vectorCentBins[iCent], vectorCentBins[iCent+1], vectorPtBins_Double_Diff[iPt], vectorPtBins_Double_Diff[iPt+1]),"", 60,-5,-2);
      histoRecJpsi_YPt->Sumw2();
    }
  }


  for(int iDimuonCharge =0; iDimuonCharge < TypeOfJPsi; iDimuonCharge++)
  {
//    for(int iCentrality =0; iCentrality < numberOfCentralityBins; iCentrality++)
//    {
        
        directoryInvMass[iDimuonCharge]->cd();
        TH1F *histoInvMass_PtSummedUp = new TH1F(Form("histoInvMass_PtSummedUp"),"",400,0,10);
        histoInvMass_PtSummedUp->Sumw2();

        TH1F *histoInvMass_Pt12toInf = new TH1F(Form("histoInvMass_Pt12toInf"),"",400,0,10);
        histoInvMass_Pt12toInf->Sumw2();

      for(int iCent=0; iCent < numberOfCentBins; iCent++)
      {

          for(int iPt =0; iPt < numberOfPtBins; iPt++)
          {
              directoryInvMass[iDimuonCharge]->cd();
              TH1F *histoInvMass = new TH1F(Form("histoInvMass_Cent%gto%g_Pt%gto%g",vectorCentBins[iCent], vectorCentBins[iCent+1], vectorPtBins[iPt],vectorPtBins[iPt+1]),"",400,0,10);
              histoInvMass->Sumw2();

              directoryPt[iDimuonCharge]->cd();
              TH1F *histoPt =new TH1F(Form("histoPt_Cent%gto%g_Pt%gto%g",vectorCentBins[iCent], vectorCentBins[iCent+1], vectorPtBins[iPt],vectorPtBins[iPt+1]),"",2000,0 ,vectorPtBins[numberOfPtBins]);
              histoPt->Sumw2();

//              directoryRapidity[iDimuonCharge]->cd();
//              TH1F *histoRapidity =new TH1F(Form("histoRapidity_Cent%gto%g",vectorCentBins[iCent], vectorCentBins[iCent+1]),"",60,-5,-2);
//              histoRapidity->Sumw2();
          }
          directoryRapidity[iDimuonCharge]->cd();
          TH1F *histoRapidity =new TH1F(Form("histoRapidity_Cent%gto%g",vectorCentBins[iCent], vectorCentBins[iCent+1]),"",60,-5,-2);
          histoRapidity->Sumw2();
          
          for(int iY=0; iY< numberOfRapidityBins; iY++)
          { 
              directoryInvMass[iDimuonCharge]->cd();
              TH1F *histoInvMass_Y = new TH1F(Form("histoInvMass_Cent%gto%g_Y%gto%g", vectorCentBins[iCent], vectorCentBins[iCent+1], vectorRapidityBins[iY],vectorRapidityBins[iY+1]),"",400,0,10);
              histoInvMass_Y->Sumw2();
           }

      }

//      directoryPt[iDimuonCharge]->cd();
//      TH1F *histoSinglePt = new TH1F(Form("histoSinglePt_Cent%gto%g",vectorCentralityBins[iCentrality],vectorCentralityBins[iCentrality+1]),"",2000,0,20);
//      histoSinglePt->Sumw2();
//    }
  }
 
   double par0_gen_pt_20to30[] = {573409.30, 3.4611, 1.9373, 3.9033, -3835.3278, 0.4426};
   double par0_rec_pt_20to30[] = {432475.11, 3.1684, 2.3165, 3.0399, 42655.9235, 0.5356};

   double par0_gen_pt_30to40[] = {314782.32, 3.4950, 1.9353, 3.9531, -1567.1548, 0.7354};
   double par0_rec_pt_30to40[] = {173027.21, 3.2066, 2.6152, 2.5931, 56462.1904, 0.5636};

   double par0_gen_pt_40to90[] = {381156.11, 3.5154, 1.9323, 3.9835, 377.7856, 0.3198};
   double par0_rec_pt_40to90[] = {264618.37, 4.3018, 1.8622, 4.4136, 12818.4470, 0.8794};

   double par0_gen_y_20to30[] = {2423655.77, -0.0006, 2.1270};
   double par0_rec_y_20to30[] = {2613800.09, 0.4396, 2.3426};

   double par0_gen_y_30to40[] = {1342126.63, 0.0241, 2.1336};
   double par0_rec_y_30to40[] = {632739.19, -1.8112, 1.6280};

   double par0_gen_y_40to90[] = {1579359.94, 0.0199, 2.1321};
   double par0_rec_y_40to90[] = {782849.79, -1.8012, 1.5408};


   double par1_gen_pt_20to30[] = {426479.52, 3.1707, 2.3160, 3.0437, 42198.9956, 0.5336};
   double par1_rec_pt_20to30[] = {434512.30, 3.1646, 2.3145, 3.0364, 39111.7936, 0.5325};

   double par1_gen_pt_30to40[] = {163783.05, 3.2142, 2.6121, 2.6063, 53408.0353, 0.5582};
   double par1_rec_pt_30to40[] = {171126.28, 3.2000, 2.6075, 2.5881, 52000.9069, 0.5631};

   double par1_gen_pt_40to90[] = {299855.45, 4.3207, 1.8491, 4.4664, 24477.4146, 2.6914};
   double par1_rec_pt_40to90[] = {262983.56, 4.2870, 1.8528, 4.4137, 10221.3701, 0.8872};

   double par1_gen_y_20to30[] = {2611994.68, 0.4369, 2.3408};
   double par1_rec_y_20to30[] = {5797987.55, 2.9636, 3.0737};

   double par1_gen_y_30to40[] = {613847.67, -1.8257, 1.6214};
   double par1_rec_y_30to40[] = {600670.07, -1.8167, 1.6815};

   double par1_gen_y_40to90[] = {748582.43, -1.8004, 1.5413};
   double par1_rec_y_40to90[] = {697106.11, -1.9658, -1.5040};



  TString TObjArrayName = Form("TTree_DimuonINT7inMUON;1");
  TString TTreeName = Form("eventsTree");
  
  TFile *anaFile = new TFile(ReadFileLocation,"READ");
  if ( anaFile->IsOpen() ) cout << "read anaFile successfully\n";

  TObjArray *foutput = dynamic_cast<TObjArray*> (anaFile->Get(TObjArrayName));
  TTree *tree = (TTree*)foutput->FindObject(TTreeName); 
   
  int nEntries = tree->GetEntries();
  cout << "num of entries: " << nEntries << endl;

    //TH1F *th1fCMSLEventCent = dynamic_cast<TH1F*> ( foutput->FindObject("th1fCMSLEventCent") );
    //TH1F *th1fCMSLEventCent = (TH1F*)foutput->FindObject("th1fCMSLEventCent");
  
//  std::vector<TLorentzVector> vectorRawMuon;
//  std::vector<TLorentzVector> vectorRawMuon2;
//  std::vector<TLorentzVector> *tempoVectorMuon =0;
//  std::vector<TLorentzVector> *tempoVectorMuon2 =0;    
  double mainEventCentrality;
  double eventCentrality;
  
//  TBranch *muon = tree->GetBranch("gen_dimuon");
//  muon->SetAddress(&vectorRawMuon);
//  muon->SetAddress(&tempoVectorMuon); 

  TBranch *centrality = tree->GetBranch("eventCentrality");
  centrality->SetAddress(&eventCentrality);

  std::cout << "set branch successfully" << std::endl;


  TLorentzVector lvDimuon, lvDimuon2, lvMuon, lvMuon2;
  int numJPsi=0;
  eventCentrality=0;
  int dimuonCharge = 0;

  for (int iRaw=0;iRaw<nEntries;iRaw++)
//for (int iRaw=1;iRaw<=nEntries;iRaw++)  
  {
    std::vector<TLorentzVector> vectorRawMuon;
    std::vector<TLorentzVector> vectorRawMuon2;
    std::vector<TLorentzVector> *tempoVectorMuon =0;
    std::vector<TLorentzVector> *tempoVectorMuon2 =0;

    TBranch *muon = tree->GetBranch("gen_dimuon");
    muon->SetAddress(&tempoVectorMuon);
  
    TBranch *muon2 = tree->GetBranch("rec_dimuon");
    muon2->SetAddress(&tempoVectorMuon2);

    tree->GetEntry(iRaw);
//    printf("processing event... %.0f%%%s", 100.*iRaw/nEntries, (iRaw < nEntries) ? "\r" : "\n");

//    if( tempoVectorMuon->empty() ) continue;
    TString sCentrality;
    sCentrality = FindCentBin(eventCentrality);
    int centralityBin = FindBinInVector(eventCentrality, vectorCentBins);
    //if( centralityBin==1 && eventCentrality <10.5 ) cout <<centralityBin <<" "<< eventCentrality << endl;
    if( centralityBin == -1 ) continue;
    if( centralityBin == 9 ) continue;
    //cout <<centralityBin <<" "<< eventCentrality << endl;
//    if( numberOfEventInPool[centralityBin] < maxNumberOfMixedEvent ) continue;

    vectorRawMuon.clear();
    dimuonCharge = 0;

// the weight for 0-10% and 10-20% centrality

    double weight0_pt=0;
    double weight0_y=0;
    double weight0=0;

    double weight1_pt=0;
    double weight1_y=0;
    double weight1=0;

    double weight2_pt=0;
    double weight2_y=0;
    double weight2=0;

    double weight_results = 0;

  
    for(int iVector=0;iVector<tempoVectorMuon->size();iVector++)
    {
      vectorRawMuon.push_back(tempoVectorMuon->at(iVector));
      lvDimuon = vectorRawMuon[iVector];
      Float_t rapidity = 0.5 * log( ( lvDimuon.E() + lvDimuon.Pz() ) / ( lvDimuon.E() - lvDimuon.Pz() ) );

      TString fillName;
      fillName = Form("JPsi:gen/run:%d/centrality:%s", RunNumber, sCentrality.Data());
      fEventCounters->Count(fillName.Data());
//      std::cout << fillName.Data() << std::endl;
//      if( ( rapidity > -2.5) || rapidity < -4) continue;
      ( (TH1F*)outputFile->Get(Form("DimuonHistos%s/InvMass/histoInvMass_PtSummedUp",dimuonChargeNames[dimuonCharge].Data())))->Fill(lvDimuon.M());
      if(lvDimuon.Pt() >= 12.) ( (TH1F*)outputFile->Get(Form("DimuonHistos%s/InvMass/histoInvMass_Pt12toInf",dimuonChargeNames[dimuonCharge].Data())))->Fill(lvDimuon.M());

      int ptBin = FindBinInVector(lvDimuon.Pt(),vectorPtBins);
      if((ptBin == -1) ) continue;

//      histoGenJpsi_Pt->Fill(lvDimuon.Pt());

      int rapidityBin = FindBinInVector(-1.*rapidity, vectorRapidityBins);
      if( (rapidityBin == -1)) continue;

      if(centralityBin == 2)
      {
        weight0_pt = setWeight_Pt(lvDimuon.Pt(), par0_rec_pt_20to30, par0_gen_pt_20to30)* setWeight_Pt(lvDimuon.Pt(), par1_rec_pt_20to30, par1_gen_pt_20to30);
        weight0_y = setWeight_Y(rapidity, par0_rec_y_20to30, par0_gen_y_20to30) * setWeight_Y(rapidity, par1_rec_y_20to30, par1_gen_y_20to30);
        weight0 = weight0_pt * weight0_y;
 
        histoWeightGenPt_20to30->Fill(weight0_pt);
        histoWeightGenY_20to30->Fill(weight0_y);
        histoWeightGen_20to30->Fill(weight0);

        weight_results = weight0;
//        std::cout << "weight value:  " << weight_results << "for " << vectorCentBins[centralityBin] <<"-- "<< vectorCentBins[centralityBin+1] << std::endl;

      }
      else if( centralityBin == 3)
      {
        weight1_pt = setWeight_Pt(lvDimuon.Pt(), par0_rec_pt_30to40, par0_gen_pt_30to40) * setWeight_Pt(lvDimuon.Pt(), par1_rec_pt_30to40, par1_gen_pt_30to40);
        weight1_y = setWeight_Y(rapidity, par0_rec_y_30to40, par0_gen_y_30to40) * setWeight_Y(rapidity, par1_rec_y_30to40, par1_gen_y_30to40);
        weight1 = weight1_pt * weight1_y;
 
        histoWeightGenPt_30to40->Fill(weight1_pt);
        histoWeightGenY_30to40->Fill(weight1_y);
        histoWeightGen_30to40->Fill(weight1);
 
        weight_results = weight1;
      }
      else if( centralityBin > 3)
      {
        weight2_pt = setWeight_Pt(lvDimuon.Pt(), par0_rec_pt_40to90, par0_gen_pt_40to90) * setWeight_Pt(lvDimuon.Pt(), par1_rec_pt_40to90, par1_gen_pt_40to90);
        weight2_y = setWeight_Y(rapidity, par0_rec_y_40to90, par0_gen_y_40to90) * setWeight_Y(rapidity, par1_rec_y_40to90, par1_gen_y_40to90);
        weight2 = weight2_pt * weight2_y;
 
        histoWeightGenPt_40to90->Fill(weight2_pt);
        histoWeightGenY_40to90->Fill(weight2_y);
        histoWeightGen_40to90->Fill(weight2);
 
        weight_results = weight2;
      }
      else 
      {
        weight0_pt = 1.0;
        weight0_y = 1.0;
        weight0 = 1.0;

        weight1_pt = 1.0;
        weight1_y = 1.0;
        weight1 = 1.0;

        weight2_pt = 1.0;
        weight2_y = 1.0;
        weight2 = 1.0;

        weight_results = 1;
      }


//      histoGenJpsi_Y->Fill(rapidity);

//      if(ptBin > 13) ( (TH1F*)outputFile->Get(Form("DimuonHistos%s/InvMass/histoInvMass_Pt12toInf",dimuonChargeNames[dimuonCharge].Data())))->Fill(lvDimuon.M());
      if( lvDimuon.Pt() < 20)
      {
        ( (TH1F*)outputFile->Get(Form("DimuonHistos%s/GenJpsi_Y_Cent%gto%g", dimuonChargeNames[dimuonCharge].Data(), vectorCentBins[centralityBin], vectorCentBins[centralityBin+1])))->Fill( rapidity, weight_results);
        ( (TH1F*)outputFile->Get(Form("DimuonHistos%s/InvMass/histoInvMass_Cent%gto%g_Y%gto%g",dimuonChargeNames[dimuonCharge].Data(), vectorCentBins[centralityBin], vectorCentBins[centralityBin+1], vectorRapidityBins[rapidityBin], vectorRapidityBins[rapidityBin+1])))->Fill(lvDimuon.M(), weight_results);
        ((TH1F*) outputFile->Get(Form("DimuonHistos%s/Rapidity/histoRapidity_Cent%gto%g",dimuonChargeNames[dimuonCharge].Data(),vectorCentBins[centralityBin], vectorCentBins[centralityBin+1])))->Fill( rapidity, weight_results);
      }

//     std::cout << "pt: " << lvDimuon.Pt() << std::endl;
            ( (TH1F*)outputFile->Get(Form("DimuonHistos%s/InvMass/histoInvMass_Cent%gto%g_Pt%gto%g",dimuonChargeNames[dimuonCharge].Data(),vectorCentBins[centralityBin], vectorCentBins[centralityBin+1], vectorPtBins[ptBin],vectorPtBins[ptBin+1])))->Fill(lvDimuon.M(), weight_results);

//          if( (lvDimuon.M() < 5) && (lvDimuon.M() > 2))
//          {
            ((TH1F*) outputFile->Get(Form("DimuonHistos%s/Pt/histoPt_Cent%gto%g_Pt%gto%g",dimuonChargeNames[dimuonCharge].Data(),vectorCentBins[centralityBin], vectorCentBins[centralityBin+1], vectorPtBins[ptBin],vectorPtBins[ptBin+1])))->Fill(lvDimuon.Pt(), weight_results);
            ((TH1F*) outputFile->Get(Form("DimuonHistos%s/GenJpsi_Pt_Cent%gto%g", dimuonChargeNames[dimuonCharge].Data(),vectorCentBins[centralityBin], vectorCentBins[centralityBin+1])))->Fill(lvDimuon.Pt(), weight_results);

//       std::cout << "filling the rapidity in double differential bins" << std::endl;
       int ptBin_Double_Diff = FindBinInVector(lvDimuon.Pt(), vectorPtBins_Double_Diff); 
       if((ptBin_Double_Diff == -1) ) continue;        
       ((TH1F*) outputFile->Get(Form("DimuonHistos%s/GenJpsi_Y_Cent%gto%g_Pt%gto%g", dimuonChargeNames[dimuonCharge].Data(),vectorCentBins[centralityBin], vectorCentBins[centralityBin+1], vectorPtBins_Double_Diff[ptBin_Double_Diff],vectorPtBins_Double_Diff[ptBin_Double_Diff+1] )))->Fill(rapidity, weight_results);

//      vectorRawMuon2.push_back(tempoVectorMuon->at(iVector));
//      cout << " tempo vector muon size: " << tempoVectorMuon->size() << endl;
    } // end of for(int iVector=0;iVector<tempoVectorMuon->size();iVector++)

    tempoVectorMuon->clear();
    vectorRawMuon.clear();

    vectorRawMuon2.clear();
    dimuonCharge = 1;
    for(int iVector=0;iVector<tempoVectorMuon2->size();iVector++)
    {
      vectorRawMuon2.push_back(tempoVectorMuon2->at(iVector));
      lvDimuon2 = vectorRawMuon2[iVector];
      Float_t rapidity = 0.5 * log( ( lvDimuon2.E() + lvDimuon2.Pz() ) / ( lvDimuon2.E() - lvDimuon2.Pz() ) );

//      std::cout << "pt: " << lvDimuon2.Pt() << std::endl;
//      std::cout << "rapidity: " << rapidity << std::endl;
//      if( ( rapidity > -2.5) || rapidity < -4) continue;
      ( (TH1F*)outputFile->Get(Form("DimuonHistos%s/InvMass/histoInvMass_PtSummedUp",dimuonChargeNames[dimuonCharge].Data())))->Fill(lvDimuon2.M());
      if(lvDimuon2.Pt() >= 12.) ( (TH1F*)outputFile->Get(Form("DimuonHistos%s/InvMass/histoInvMass_Pt12toInf",dimuonChargeNames[dimuonCharge].Data())))->Fill(lvDimuon2.M());
      int ptBin = FindBinInVector(lvDimuon2.Pt(),vectorPtBins);
      if((ptBin == -1) ) continue;

//      histoRecJpsi_Pt->Fill(lvDimuon2.Pt());
      int rapidityBin = FindBinInVector(-1.*rapidity, vectorRapidityBins);
      if( (rapidityBin == -1)) continue;

//      histoRecJpsi_Y->Fill(rapidity);
//      if(ptBin > 13) ( (TH1F*)outputFile->Get(Form("DimuonHistos%s/InvMass/histoInvMass_Pt12toInf",dimuonChargeNames[dimuonCharge].Data())))->Fill(lvDimuon2.M());
      if( lvDimuon2.Pt() < 20)
      {
         ( (TH1F*)outputFile->Get(Form("DimuonHistos%s/RecJpsi_Y_Cent%gto%g", dimuonChargeNames[dimuonCharge].Data(), vectorCentBins[centralityBin], vectorCentBins[centralityBin+1])))->Fill( rapidity, weight_results);
//         std::cout << "rec region, weight value: " <<  weight_results << std::endl;
         ( (TH1F*)outputFile->Get(Form("DimuonHistos%s/InvMass/histoInvMass_Cent%gto%g_Y%gto%g",dimuonChargeNames[dimuonCharge].Data(), vectorCentBins[centralityBin], vectorCentBins[centralityBin+1], vectorRapidityBins[rapidityBin], vectorRapidityBins[rapidityBin+1] )))->Fill(lvDimuon2.M(), weight_results);
         ((TH1F*) outputFile->Get(Form("DimuonHistos%s/Rapidity/histoRapidity_Cent%gto%g",dimuonChargeNames[dimuonCharge].Data(), vectorCentBins[centralityBin], vectorCentBins[centralityBin+1] )))->Fill( rapidity, weight_results);

      }

      ( (TH1F*)outputFile->Get(Form("DimuonHistos%s/InvMass/histoInvMass_Cent%gto%g_Pt%gto%g",dimuonChargeNames[dimuonCharge].Data(), vectorCentBins[centralityBin], vectorCentBins[centralityBin+1], vectorPtBins[ptBin],vectorPtBins[ptBin+1])))->Fill(lvDimuon2.M(), weight_results);

//          if( (lvDimuon.M() < 5) && (lvDimuon.M() > 2))
//          {
      ((TH1F*) outputFile->Get(Form("DimuonHistos%s/Pt/histoPt_Cent%gto%g_Pt%gto%g",dimuonChargeNames[dimuonCharge].Data(), vectorCentBins[centralityBin], vectorCentBins[centralityBin+1], vectorPtBins[ptBin],vectorPtBins[ptBin+1])))->Fill(lvDimuon2.Pt(), weight_results);     
      ((TH1F*) outputFile->Get(Form("DimuonHistos%s/RecJpsi_Pt_Cent%gto%g", dimuonChargeNames[dimuonCharge].Data(),vectorCentBins[centralityBin], vectorCentBins[centralityBin+1])))->Fill(lvDimuon2.Pt(), weight_results); 

       int ptBin_Double_Diff = FindBinInVector(lvDimuon2.Pt(), vectorPtBins_Double_Diff);
       if((ptBin_Double_Diff == -1) ) continue;        
       ((TH1F*) outputFile->Get(Form("DimuonHistos%s/RecJpsi_Y_Cent%gto%g_Pt%gto%g", dimuonChargeNames[dimuonCharge].Data(),vectorCentBins[centralityBin], vectorCentBins[centralityBin+1], vectorPtBins_Double_Diff[ptBin_Double_Diff],vectorPtBins_Double_Diff[ptBin_Double_Diff+1] )))->Fill(rapidity, weight_results);


    } // end of for(int iVector=0;iVector<tempoVectorMuon2->size();iVector++)
    tempoVectorMuon2->clear();
    vectorRawMuon2.clear();

  } //end for(int iRaw=0;iRaw<nEntries;iRaw++)
   

//   TObjArrayName = Form("TTree_rec_DimuonMULU;1");
//    TTreeName = Form("eventsTree2");
  
//    TObjArray *foutput2 = dynamic_cast<TObjArray*> (anaFile->Get(TObjArrayName));
//    TTree *tree2 = (TTree*)foutput2->FindObject(TTreeName); 
   
//    nEntries = tree->GetEntries();
//    cout << "num of entries of rec: " << nEntries << endl;
/*
   TBranch *muon2 = tree->GetBranch("rec_dimuon");
   muon2->SetAddress(&tempoVectorMuon2); 

   centrality = tree->GetBranch("eventCentrality");
   centrality->SetAddress(&eventCentrality);

  std::cout << "set branch successfully" << std::endl;

   dimuonCharge = 1;

  for (int iRaw=0;iRaw<nEntries;iRaw++)
//for (int iRaw=1;iRaw<=nEntries;iRaw++)  
  {
    tree->GetEntry(iRaw);
//    printf("processing event... %.0f%%%s", 100.*iRaw/nEntries, (iRaw < nEntries) ? "\r" : "\n");
//    if( tempoVectorMuon->empty() ) continue;

    int centralityBin = FindBinInVector(eventCentrality, vectorCentBins);
    //if( centralityBin==1 && eventCentrality <10.5 ) cout <<centralityBin <<" "<< eventCentrality << endl;
    if( centralityBin == -1 ) continue;
//  cout <<centralityBin <<" "<< eventCentrality << endl;
    
    vectorRawMuon2.clear();
  
    for(int iVector=0;iVector<tempoVectorMuon2->size();iVector++)
    {
      vectorRawMuon2.push_back(tempoVectorMuon2->at(iVector));
      lvDimuon2 = vectorRawMuon2[iVector];
      Float_t rapidity = 0.5 * log( ( lvDimuon2.E() + lvDimuon2.Pz() ) / ( lvDimuon2.E() - lvDimuon2.Pz() ) );

//      std::cout << "pt: " << lvDimuon2.Pt() << std::endl;
//      std::cout << "rapidity: " << rapidity << std::endl;
//      if( ( rapidity > -2.5) || rapidity < -4) continue;
      ( (TH1F*)outputFile->Get(Form("DimuonHistos%s/InvMass/histoInvMass_PtSummedUp",dimuonChargeNames[dimuonCharge].Data())))->Fill(lvDimuon2.M());
      if(lvDimuon2.Pt() >= 12.) ( (TH1F*)outputFile->Get(Form("DimuonHistos%s/InvMass/histoInvMass_Pt12toInf",dimuonChargeNames[dimuonCharge].Data())))->Fill(lvDimuon2.M());
      int ptBin = FindBinInVector(lvDimuon2.Pt(),vectorPtBins);
      if((ptBin == -1) ) continue;

      histoRecJpsi_Pt->Fill(lvDimuon2.Pt());
      int rapidityBin = FindBinInVector(-1.*rapidity, vectorRapidityBins);
      if( (rapidityBin == -1)) continue;
   
      histoRecJpsi_Y->Fill(rapidity);
//      if(ptBin > 13) ( (TH1F*)outputFile->Get(Form("DimuonHistos%s/InvMass/histoInvMass_Pt12toInf",dimuonChargeNames[dimuonCharge].Data())))->Fill(lvDimuon2.M());
      if( lvDimuon2.Pt() < 20)
      {
         ( (TH1F*)outputFile->Get(Form("DimuonHistos%s/InvMass/histoInvMass_Cent%gto%g_Y%gto%g",dimuonChargeNames[dimuonCharge].Data(), vectorCentBins[centralityBin], vectorCentBins[centralityBin+1], vectorRapidityBins[rapidityBin], vectorRapidityBins[rapidityBin+1] )))->Fill(lvDimuon2.M());
         ((TH1F*) outputFile->Get(Form("DimuonHistos%s/Rapidity/histoRapidity_Cent%gto%g",dimuonChargeNames[dimuonCharge].Data(), vectorCentBins[centralityBin], vectorCentBins[centralityBin+1] )))->Fill( rapidity);
         
      }

      ( (TH1F*)outputFile->Get(Form("DimuonHistos%s/InvMass/histoInvMass_Cent%gto%g_Pt%gto%g",dimuonChargeNames[dimuonCharge].Data(), vectorCentBins[centralityBin], vectorCentBins[centralityBin+1], vectorPtBins[ptBin],vectorPtBins[ptBin+1])))->Fill(lvDimuon2.M());

//          if( (lvDimuon.M() < 5) && (lvDimuon.M() > 2))
//          {
      ((TH1F*) outputFile->Get(Form("DimuonHistos%s/Pt/histoPt_Cent%gto%g_Pt%gto%g",dimuonChargeNames[dimuonCharge].Data(), vectorCentBins[centralityBin], vectorCentBins[centralityBin+1], vectorPtBins[ptBin],vectorPtBins[ptBin+1])))->Fill(lvDimuon2.Pt());


    }
  }
*/


  anaFile->Close();
  outputFile->Write();
  delete outputFile;
}

std::vector<double> StringToVector(TString strBins = "0,90")
{
  std::vector<double> vectorBins;
  TObjArray *objBins = strBins.Tokenize(",");
  TIter nextBin(objBins);
  TObjString *strBin;
  while ((strBin=(TObjString*)nextBin())) {
    TString currentBin = strBin->GetName();
    Double_t bin = currentBin.Atof();
    vectorBins.push_back(bin);
  }
  return vectorBins;
}

int FindBinInVector( double valueToBin =0, std::vector<double> vectorBins = std::vector<double>() )
{
  for(int iElement =0; iElement < vectorBins.size()-1; iElement++)
  {
    if( (valueToBin < vectorBins[iElement+1]) && (valueToBin >= vectorBins[iElement]) )
    return iElement;
  }
return -1;
}

TString FindCentBin(double lPercentile)
{
       TString sCentrality;
       if(lPercentile < 0.0 )
       {
           sCentrality = "m0"; 
       }else if(lPercentile <= 10.0 && lPercentile >=0.0)
       {
           sCentrality = "0_10";
           //cout << "sCentrality: " << sCentrality.Data() << endl;
       }else if(lPercentile <= 20.0 && lPercentile >10.0)
       {
           sCentrality = "10_20";
           //cout << "sCentrality: " << sCentrality.Data() << endl;
       }else if(lPercentile <= 30.0 && lPercentile >20.0)
       {
           sCentrality = "20_30";
           //cout << "sCentrality: " << sCentrality.Data() << endl;
       }else if(lPercentile <= 40.0 && lPercentile >30.0)
       {
           sCentrality = "30_40";
           //cout << "sCentrality: " << sCentrality.Data() << endl;
       }else if(lPercentile <= 50.0 && lPercentile >40.0)
        {
           sCentrality = "40_50";
           //cout << "sCentrality: " << sCentrality.Data() << endl;
       }else if(lPercentile <= 60.0 && lPercentile >50.0)
       {
           sCentrality = "50_60";
           //cout << "sCentrality: " << sCentrality.Data() << endl;
       }else if(lPercentile <= 70.0 && lPercentile >60.0)
       {
           sCentrality = "60_70";
           //cout << "sCentrality: " << sCentrality.Data() << endl;
       }else if(lPercentile <= 80.0 && lPercentile >70.0)
       {
            sCentrality = "70_80";
           //cout << "sCentrality: " << sCentrality.Data() << endl;
       }else if(lPercentile <= 90.0 && lPercentile >80.0)
       {
            sCentrality = "80_90";
            //cout << "sCentrality: " << sCentrality.Data() << endl;
       }else if(lPercentile <= 110.0 && lPercentile >90.0)
       {
             sCentrality = "90_100";
            //cout << "sCentrality: " << sCentrality.Data() << endl;
       }

       return sCentrality;
}

double setWeight_Pt(double Pt, double *parRaw, double *parGen)
{  

//  TFormula *fRecEval = new TFormula("fRecEval","[0]*x[0] / TMath::Power(1 + TMath::Power(x[0]/[1],[2]),[3] )");
  TFormula *fRecEval = new TFormula("fRecEval","[0]*x[0] / TMath::Power(1 + TMath::Power(x[0]/[1],[2]),[3] )+[4]*TMath::Exp(-1.0*[5]*x[0])");
//    TF1 *fRecEval = new TF1("fRecEval","[0]*x[0] / TMath::Power( (1 + [1]*x[0]*x[0]) , [2] )",0,15);
  fRecEval->SetParameter(0, parRaw[0] );
  fRecEval->SetParameter(1, parRaw[1] );
  fRecEval->SetParameter(2, parRaw[2] );
  fRecEval->SetParameter(3, parRaw[3] );
  
  fRecEval->SetParameter(4, parRaw[4] );
  fRecEval->SetParameter(5, parRaw[5] );
  
//  TFormula *fGenEval = new TFormula("fGenEval","[0]*x[0] / TMath::Power(1 + TMath::Power(x[0]/[1],[2]),[3] )");
  TFormula *fGenEval =new TFormula("fGenEval","[0]*x[0] / TMath::Power(1 + TMath::Power(x[0]/[1],[2]),[3] )+[4]*TMath::Exp(-1.0*[5]*x[0])");
  fGenEval->SetParameter(0, parGen[0] );
  fGenEval->SetParameter(1, parGen[1] );
  fGenEval->SetParameter(2, parGen[2] );
  fGenEval->SetParameter(3, parGen[3] );

  fGenEval->SetParameter(4, parGen[4] );
  fGenEval->SetParameter(5, parGen[5] );
  
  double reco = fRecEval->Eval( Pt );
  double simu = fGenEval->Eval( Pt );
  
  if( (simu==0) || (reco/simu <= 0) )
  {
       return 1;
  }
  
  delete fRecEval;
  delete fGenEval;
  return reco/simu;
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
