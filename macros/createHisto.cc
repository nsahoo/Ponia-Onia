#include<TString.h>
#include <TLatex.h>


void createHisto(){

  gStyle->SetOptStat(10);
  TChain *ch = new TChain("rootuple/oniaTree");
  //ch->Add("small/Onia2MuMuRootuple-Run2015B-MuOnia.root");
  //ch->Add("MuOnia-2015B-data-*.root");
  ch->Add("Charmonium-2015B-data-*.root");
  TTree *tr = ch;

  int nentries = tr->GetEntries();
  printf("total entries in oniaTree: %i \n", nentries);


  const unsigned int histoN = 18;
  TH1D *h[histoN];
  TCanvas *cv[histoN];
  TString hltBit[histoN] = {"passbit0", "passbit1", "passbit2", "passbit3", "passbit4", "passbit5", "passbit6", "passbit7", "passbit8", "passbit9", 
			    "passbit10", "passbit11", "passbit12", "passbit13", "passbit14", "passbit15", "passbit16", "passbit17"};

  char hltName[18][32] = {"HLT_Dimuon16_Jpsi_v", "HLT_Dimuon20_Jpsi_v", "HLT_Dimuon13_PsiPrime_v", "HLT_Dimuon13_Upsilon_v", "HLT_Dimuon10_Jpsi_Barrel_v", 
                             "HLT_Dimuon8_PsiPrime_Barrel_v", "HLT_Dimuon8_Upsilon_Barrel_v", "HLT_Dimuon0_Phi_Barrel_v", "HLT_Mu25_TkMu0_dEta18_Onia_v", "HLT_Mu16_TkMu0_dEta18_Onia_v",
			     "HLT_Mu7p5_L2Mu2_Jpsi_v", "HLT_Mu7p5_L2Mu2_Upsilon_v", "HLT_Mu7p5_Track2_Jpsi_v", "HLT_Mu7p5_Track3p5_Jpsi_v", 
			     "HLT_Mu7p5_Track7_Jpsi_v", "HLT_Mu7p5_Track2_Upsilon_v", "HLT_Mu7p5_Track3p5_Upsilon_v", "HLT_Mu7p5_Track7_Upsilon_v"};

  for(int i=0; i<histoN; i++){
    printf("--------------------------------------\n");
    printf("HLT path: %s \n", hltName[i]);
    h[i] = new TH1D(Form("h[%i]",i),"; M(#mu^{+}#mu^{-}) [GeV/c^{2}]; Events per GeV",300,0,300);
    cv[i] = new TCanvas(Form("cv[%i]",i),Form("canvas_%i",i),800,600);

    tr->Draw(Form("dimuon_p4.M() >> h[%i]",i),hltBit[i]);
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
    cv[i]->SetLogx();
    cv[i]->SetLogy();
    //cv[i]->Print(Form("MuOnia/plots_%i.pdf",i));
    //cv[i]->Print(Form("Charmonium/plots_%i.pdf",i));

    delete t1;
  }


}
