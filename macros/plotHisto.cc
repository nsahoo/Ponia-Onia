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


using namespace std;


void plotHisto(TString dataset, TString varPlot){


  gStyle->SetOptStat(10);

  TChain *ch = new TChain("rootuple/oniaTree");
  if (dataset=="MuOnia"){
    printf("accessing MuOnia files\n");
    ch->Add("MuOnia-2015B-data-v3-30jul-merged.root");
  } else if(dataset=="Charmonium"){
    printf("accessing Charmonium files\n");
    ch->Add("Charmonium-2015B-data-v3-30jul-merged.root");
  }


  TString type;
  TString pt_type;
  TString eta_type;
  if (varPlot=="1"){
    printf("plotting mu+ pt,eta\n");
    type="mup";
    //pt_type=="muonP_p4.Pt()";
    //eta_type=="muonP_p4.Eta()";
  } else if (varPlot=="2"){
    printf("plotting mu- pt,eta\n");
    type="mun";
    //pt_type=="muonN_p4.Pt()";
    //eta_type=="muonN_p4.Eta()";
  } else if (varPlot=="3"){
    printf("plotting dimuon pt, rapidity \n");
    type="dimu";
    //pt_type=="dimuon_p4.Pt()";
    //eta_type=="dimuon_p4.Rapidity()";
  } else if (varPlot=="4"){
    printf("plotting dimuon inv. mass \n");
    type="dimu";
  }

  TTree *tr = ch;
  int nentries = tr->GetEntries();
  printf("total entries in oniaTree: %i \n", nentries);


  const unsigned int histoN = 19;                                                                                                                                                              
  TH1D *h[histoN];                                                                                                                                                                             
  TH1D *g[histoN];
  TCanvas *cv[histoN];                               
  //TCanvas *cn[histoN];                                                                                                                                        
  TString hltBit[histoN] = {"passbit0",  "passbit1",  "passbit2",  "passbit3",  "passbit4",  "passbit5",  "passbit6",  "passbit7", "passbit8", "passbit9",                                             
                            "passbit10", "passbit11", "passbit12", "passbit13", "passbit14", "passbit15", "passbit16", "passbit17", "passbit18"};


  char hltName[19][32] = {"HLT_Dimuon16_Jpsi_v*", "HLT_Dimuon20_Jpsi_v*", "HLT_Dimuon13_PsiPrime_v*", "HLT_Dimuon13_Upsilon_v*", 
			  "HLT_Dimuon10_Jpsi_Barrel_v*", "HLT_Dimuon8_PsiPrime_Barrel_v*", "HLT_Dimuon8_Upsilon_Barrel_v*", "HLT_Dimuon0_Phi_Barrel_v*", 
			  "HLT_Mu25_TkMu0_dEta18_Onia_v*", "HLT_Mu16_TkMu0_dEta18_Onia_v*", "HLT_Mu16_TkMu0_dEta18_Phi_v*",
			  "HLT_Mu7p5_L2Mu2_Jpsi_v*", "HLT_Mu7p5_L2Mu2_Upsilon_v*", "HLT_Mu7p5_Track2_Jpsi_v*", "HLT_Mu7p5_Track3p5_Jpsi_v*",                                                       
			  "HLT_Mu7p5_Track7_Jpsi_v*", "HLT_Mu7p5_Track2_Upsilon_v*", "HLT_Mu7p5_Track3p5_Upsilon_v*", "HLT_Mu7p5_Track7_Upsilon_v*"};                                          

  double xmin[19] = {2., 2., 3., 8., 2., 3., 8., 0.5, 2., 2., 0.5, 2., 8., 2., 2., 2., 8., 8., 8.};
  double xmax[19] = {4., 4., 5., 15.,4., 5.,15., 2., 13., 13., 2., 4., 15.,4., 4., 4.,15., 15., 15.};
  double nbins[19] = {20.,20.,20.,70.,20.,20.,70.,15.,110.,110.,15.,20.,70.,20.,20.,20.,70.,70.,70.};   

 
  for(int i=0; i<histoN; i++){
    printf("--------------------------------------\n");
    printf("HLT path: %s \n", hltName[i]);
    //h[i] = new TH1D(Form("h[%i]",i),"; M(#mu^{+}#mu^{-}) [GeV/c^{2}]; Events per GeV",300,0,300);
    cv[i] = new TCanvas(Form("cv[%i]",i),Form("canvas_%i",i),800,600);
    //cn[i] = new TCanvas(Form("cn[%i]",i),Form("canvas_%i",i),800,600);
    if (varPlot=="1"){
      h[i] = new TH1D(Form("h[%i]",i),"; P_{T}^{#mu^{+}}[GeV]; Events/0.5",150,0,75);
      g[i] = new TH1D(Form("g[%i]",i),"; #eta^{#mu^{+}}; Events/0.1",60,-3.,3.);
      tr->Draw(Form("muonP_p4.Pt() >> h[%i]",i),hltBit[i]);
      tr->Draw(Form("muonP_p4.Eta() >> g[%i]",i),hltBit[i]);
    } else if (varPlot=="2"){
      h[i] = new TH1D(Form("h[%i]",i),"; P_{T}^{#mu^{-}}[GeV]; Events/0.5",150,0,75);
      g[i] = new TH1D(Form("g[%i]",i),"; #eta^{#mu^{-}}; Events/0.1",60,-3.,3.);
      tr->Draw(Form("muonN_p4.Pt() >> h[%i]",i),hltBit[i]);
      tr->Draw(Form("muonN_p4.Eta() >> g[%i]",i),hltBit[i]);
    } else if (varPlot=="3"){
      h[i] = new TH1D(Form("h[%i]",i),"; P_{T}^{#mu#mu}[GeV]; Events/0.5",150,0,75);
      g[i] = new TH1D(Form("g[%i]",i),"; #eta^{#mu#mu}; Events/0.1",60,-3.,3.);
      tr->Draw(Form("dimuon_p4.Pt() >> h[%i]",i),hltBit[i]);
      tr->Draw(Form("dimuon_p4.Rapidity() >> g[%i]",i),hltBit[i]);
    } else if (varPlot=="4"){
      h[i] = new TH1D(Form("h[%i]",i),"; M(#mu^{+}#mu^{-}) [GeV/c^{2}]; Events/0.1",nbins[i],xmin[i],xmax[i]);
      tr->Draw(Form("dimuon_p4.M() >> h[%i]",i),hltBit[i]);
    }


    if (varPlot=="1" || varPlot=="2" || varPlot=="3"){
    //tr->Draw(Form("dimuon_p4.Pt() >> h[%i]",i),hltBit[i]);
    //tr->Draw(Form("pt_type >> h[%i]",i),hltBit[i]);
    h[i]->SetLineColor(1);
    //h[i]->SetXMinimum(0.1);
    h[i]->Draw();
    h[i]->GetYaxis()->SetTitleOffset(1.4);
    //h[i]->SetAxisRange(0.1, 300.,"X");


    TLatex *t1 = new TLatex();
    t1->SetNDC();
    //t1->SetTextFont(12);
    t1->SetTextColor(1);
    t1->SetTextSize(0.03);
    t1->SetTextAlign(12);
    double fixNDC = 0.;
    t1->DrawLatex(.35,.86+fixNDC,TString::Format("%s",hltName[i]));
    t1->DrawLatex(.19,.92,TString::Format("CMS Preliminary"));
    t1->DrawLatex(.54,.92,TString::Format("13 TeV data"));
    ////cv[i]->SetLogx();
    ////cv[i]->SetLogy();
    //cv[i]->Print(Form("MuOnia/plots_%i.pdf",i));
    //cv[i]->Print(Form("Charmonium/plots_%i.pdf",i));
    cv[i]->Print(Form("plots-v3-30jul/%s/%s_pT_%s.pdf",dataset.Data(),type.Data(),hltName[i]));
    //cv[i]->Print(Form("plots-v3-30jul/%s_%s.png",dataset.Data(),hltName[i]));
    //----------------------------------------------------------------------------------------------------
    cv[i]->Update();
    //tr->Draw(Form("dimuon_p4.Rapidity() >> g[%i]",i),hltBit[i]);
    //tr->Draw(Form("pt_type >> h[%i]",i),hltBit[i]);                                                                                                                                                 
    g[i]->SetLineColor(1);
    //h[i]->SetXMinimum(0.1);                                                                                                                                                                       
    g[i]->Draw();
    g[i]->GetYaxis()->SetTitleOffset(1.4);
    //h[i]->SetAxisRange(0.1, 300.,"X");                                                                                                                                                             
    t1->DrawLatex(.35,.86+fixNDC,TString::Format("%s",hltName[i]));
    t1->DrawLatex(.19,.92,TString::Format("CMS Preliminary"));
    t1->DrawLatex(.54,.92,TString::Format("13 TeV data"));

    ////cv[i]->SetLogx();                                                                                                                                                                             
    ////cv[i]->SetLogy();                                                                                                                                                                           
    //cv[i]->Print(Form("MuOnia/plots_%i.pdf",i));                                                                                                                                                     
    //cv[i]->Print(Form("Charmonium/plots_%i.pdf",i));                                                                                                                                                
    //cv[i]->Update();
    cv[i]->Print(Form("plots-v3-30jul/%s/%s_eta_%s.pdf",dataset.Data(),type.Data(),hltName[i]));
    //cv[i]->Print(Form("plots-v3-30jul/%s_%s.png",dataset.Data(),hltName[i]));      

    delete t1;
    //delete t2;
    } else if (varPlot=="4"){
      h[i]->SetLineColor(1);
      //h[i]->SetXMinimum(0.1);                                                                 
      h[i]->Draw();
      h[i]->GetYaxis()->SetTitleOffset(1.4);
      //h[i]->SetAxisRange(0.1, 300.,"X");                                                                                                                                                  

      TLatex *t1 = new TLatex();
      t1->SetNDC();
      //t1->SetTextFont(12);                                                                                                                                                            
      t1->SetTextColor(1);
      t1->SetTextSize(0.03);
      t1->SetTextAlign(12);
      double fixNDC = 0.;
      t1->DrawLatex(.35,.86+fixNDC,TString::Format("%s",hltName[i]));
      t1->DrawLatex(.19,.92,TString::Format("CMS Preliminary"));
      t1->DrawLatex(.54,.92,TString::Format("13 TeV data"));
      cv[i]->Print(Form("plots-v3-30jul/%s/%s_mass_%s.pdf",dataset.Data(),type.Data(),hltName[i]));

      delete t1;
    }

  }
 
}

int main(int argc, char** argv){

  if (argc == 3){
    TString dataset = argv[1];
    TString varPlot = argv[2];

    printf("dataset used: %s \n", dataset.Data());
    printf("plotted var: %s \n", varPlot.Data());

    plotHisto(dataset,varPlot);

  } else {

    printf("You haven't put the right number of arguments, please give the arguments like below..\n");
    printf("./plotHisto [dataset] [varPlot] \n");
    printf("dataset: MuOnia Charmonium \n");
    printf("varPlot: 1 2 3 4\n");

  }
  

}
