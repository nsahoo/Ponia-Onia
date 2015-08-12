//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//  macro to calculate HLT efficiency as a function of pT, eta
//  author: Niladribihari Sahoo, <nsahoo@cern.ch>
//  date: 8 august 2015, SAT 13:10
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include<iostream>
#include <iomanip> //setprecision                                                                                                                                                       
#include <string.h>
#include <TSystem.h>
#include <TString.h>
#include <TLatex.h>
#include <TChain.h>
#include <TApplication.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TGaxis.h>
#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include<cmath>
#include <stdlib.h>
#include "TEfficiency.h"


using namespace std;

void eff_9aug(TString dataset){

  gStyle->SetOptStat(10);
  TChain *ch = new TChain("rootuple/oniaTree");
  if (dataset=="JpsiToMuMu"){
    printf("accessing JpsiToMuMu files\n");
    ch->Add("Onia2MuMuRootuple-JpsiToMuMu.root");
  } else if(dataset=="Upsilon1sToMuMu"){
    printf("accessing Upsilon1sToMuMu files\n");
    ch->Add("Onia2MuMuRootuple-Upsilon1SToMuMu.root");
  }

  TTree *tr = ch;
  int nentries = tr->GetEntries();
  printf("total entries in oniaTree: %i \n", nentries);

  const unsigned int histoN = 19;

  TString hltBit[histoN] = {"passbit0",  "passbit1",  "passbit2",  "passbit3",  "passbit4",  "passbit5",  "passbit6",  "passbit7", "passbit8", "passbit9",
                            "passbit10", "passbit11", "passbit12", "passbit13", "passbit14", "passbit15", "passbit16", "passbit17", "passbit18"};

  char hltName[19][32] = {"HLT_Dimuon16_Jpsi_v*", "HLT_Dimuon20_Jpsi_v*", "HLT_Dimuon13_PsiPrime_v*", "HLT_Dimuon13_Upsilon_v*",
                          "HLT_Dimuon10_Jpsi_Barrel_v*", "HLT_Dimuon8_PsiPrime_Barrel_v*", "HLT_Dimuon8_Upsilon_Barrel_v*", "HLT_Dimuon0_Phi_Barrel_v*",
                          "HLT_Mu25_TkMu0_dEta18_Onia_v*", "HLT_Mu16_TkMu0_dEta18_Onia_v*", "HLT_Mu16_TkMu0_dEta18_Phi_v*",
                          "HLT_Mu7p5_L2Mu2_Jpsi_v*", "HLT_Mu7p5_L2Mu2_Upsilon_v*", "HLT_Mu7p5_Track2_Jpsi_v*", "HLT_Mu7p5_Track3p5_Jpsi_v*",
                          "HLT_Mu7p5_Track7_Jpsi_v*", "HLT_Mu7p5_Track2_Upsilon_v*", "HLT_Mu7p5_Track3p5_Upsilon_v*", "HLT_Mu7p5_Track7_Upsilon_v*"};



  TEfficiency *eff_pt[histoN] = {0};
  TEfficiency *eff_eta[histoN] = {0};

  TH1D *h_pT_total = new TH1D("h_pT_total","total vs P_{T}; P_{T}^{#mu#mu}[GeV]; Events/1",75,0,75);
  TH1D *h_eta_total = new TH1D("h_eta_total","total vs #eta; #eta^{#mu#mu}; Events/0.1",60,-3.,3.);
  TCanvas *c = new TCanvas("c","",800,600);
  tr->Draw("dimuon_p4.Pt() >> h_pT_total");
  h_pT_total->Draw();
  c->Print(Form("effplots/%s_h_total_pT.pdf",dataset.Data()));
  c->Update();
  tr->Draw("dimuon_p4.Rapidity() >> h_eta_total");
  h_eta_total->Draw();
  c->Print(Form("effplots/%s_h_total_eta.pdf",dataset.Data()));
  //---------------------------------------
  TH1D *h_pT_pass[histoN];
  TH1D *h_eta_pass[histoN];
  TCanvas *cv[histoN];

  for (unsigned int i=0; i<histoN; i++){
    printf("--------------------------------------\n");
    printf("HLT path: %s \n", hltName[i]);

    cv[i] = new TCanvas(Form("cv[%i]",i),Form("canvas_%i",i),800,600);
    h_pT_pass[i] = new TH1D(Form("h_pT_pass[%i]",i),"pass vs P_{T}; P_{T}^{#mu#mu}[GeV]; Events/1",75,0,75);
    h_eta_pass[i] = new TH1D(Form("h_eta_pass[%i]",i),"pass vs #eta; #eta^{#mu#mu}; Events/0.1",60,-3.,3.);
    tr->Draw(Form("dimuon_p4.Pt() >> h_pT_pass[%i]",i),hltBit[i]);

    TLatex *t1 = new TLatex();
    t1->SetNDC();
    //t1->SetTextFont(12);                                                                                                                                                                                  
    t1->SetTextColor(1);
    t1->SetTextSize(0.03);
    t1->SetTextAlign(12);
    double fixNDC = 0.;
    t1->DrawLatex(.35,.26+fixNDC,TString::Format("%s",hltName[i]));
    cv[i]->Print(Form("effplots/%s_h_pass_pT_%s.pdf",dataset.Data(),hltName[i]));
    cv[i]->Update();

    tr->Draw(Form("dimuon_p4.Rapidity() >> h_eta_pass[%i]",i),hltBit[i]);
    t1->DrawLatex(.35,.26+fixNDC,TString::Format("%s",hltName[i]));
    cv[i]->Print(Form("effplots/%s_h_pass_eta_%s.pdf",dataset.Data(),hltName[i]));
    cv[i]->Update();
    //----------------------------------------------------------------
    if(TEfficiency::CheckConsistency(*h_pT_pass[i],*h_pT_total)){                                                                                                                                             
      eff_pt[i] = new TEfficiency(*h_pT_pass[i],*h_pT_total);                                                                                                                                            
      eff_pt[i]->SetTitle("HLT efficiency vs P_{T}; P_{T}^{#mu#mu}[GeV] ; HLT efficiency");                                                                                                       
      eff_pt[i]->Draw("AP");                        
    }

    t1->DrawLatex(.35,.26+fixNDC,TString::Format("%s",hltName[i]));
    cv[i]->Print(Form("effplots/%s_eff_pT_%s.pdf",dataset.Data(),hltName[i]));
    cv[i]->Update();

    if(TEfficiency::CheckConsistency(*h_eta_pass[i],*h_eta_total)){
      eff_eta[i] = new TEfficiency(*h_eta_pass[i],*h_eta_total);
      eff_eta[i]->SetTitle("HLT efficiency vs #eta; #eta^{#mu#mu} ; HLT efficiency");
      eff_eta[i]->Draw("AP");
    }

    t1->DrawLatex(.35,.26+fixNDC,TString::Format("%s",hltName[i]));
    cv[i]->Print(Form("effplots/%s_eff_eta_%s.pdf",dataset.Data(),hltName[i]));

  }

}

int main(int argc, char** argv){

  if (argc == 2){
    TString dataset = argv[1];    
    printf("dataset used: %s \n", dataset.Data());
    eff_9aug(dataset);

  } else {

    printf("You haven't put the right number of arguments, please give the arguments like below..\n");
    printf("./eff9 [dataset]  \n");
    printf("dataset: JpsiToMuMu Upsilon1sToMuMu \n");
    
  }


}


