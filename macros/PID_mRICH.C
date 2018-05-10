#include "TStyle.h"
#include "TString.h"
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TH1D.h"
#include "TH2D.h"
#include "string.h"

using namespace std;

void PID_mRICH(const string inputfile = "/work/eic/xusun/output/modular_rich/test/PID_Likelihood_test.root") 
{
  Int_t           pid;
  Float_t         px;
  Float_t         py;
  Float_t         pz;
  Float_t         vx;
  Float_t         vy;
  Float_t         vz;
  Float_t         theta;
  Float_t         phi;
  Int_t           nhit;
  Int_t           nhitAegl;
  Int_t           nhitPhoDet;
  Int_t           nelSbKCs;
  Int_t           nelGaAsP;
  Int_t           nelGaAs;
  Double_t        Lpion;
  Double_t        LKaon;
  Double_t        Lproton;
  
  TChain *fchain  = new TChain("LLTreeDst");
  fchain->AddFile(inputfile.c_str(),-1,"LLTreeDst");
  fchain->SetMakeClass(1);
  fchain->SetBranchAddress("pid", &pid);
  fchain->SetBranchAddress("px", &px);
  fchain->SetBranchAddress("py", &py);
  fchain->SetBranchAddress("pz", &pz);
  fchain->SetBranchAddress("vx", &vx);
  fchain->SetBranchAddress("vy", &vy);
  fchain->SetBranchAddress("vz", &vz);
  fchain->SetBranchAddress("theta", &theta);
  fchain->SetBranchAddress("phi", &phi);
  fchain->SetBranchAddress("nhit", &nhit);
  fchain->SetBranchAddress("nhitAegl", &nhitAegl);
  fchain->SetBranchAddress("nhitPhoDet", &nhitPhoDet);
  fchain->SetBranchAddress("nelSbKCs", &nelSbKCs);
  fchain->SetBranchAddress("nelGaAsP", &nelGaAsP);
  fchain->SetBranchAddress("nelGaAs", &nelGaAs);
  fchain->SetBranchAddress("Lpion", &Lpion);
  fchain->SetBranchAddress("LKaon", &LKaon);
  fchain->SetBranchAddress("Lproton", &Lproton);

  TH2D *h_likelihood_piKplus = new TH2D("h_likelihood_piKplus","h_likelihood_piKplus",65,2.5,15.5,450,-20.0,430.0); // pi to K
  TH2D *h_likelihood_pipplus = new TH2D("h_likelihood_pipplus","h_likelihood_pipplus",65,2.5,15.5,450,-20.0,430.0); // pi to p
  double counter_pipiplus = 0.0;
  double counter_piKplus = 0.0;
  double counter_pipplus = 0.0;
  TH1D *h_prob_pipiplus = new TH1D("h_prob_pipiplus","h_prob_pipiplus",65,2.5,15.5); 
  TH1D *h_prob_piKplus  = new TH1D("h_prob_piKplus","h_prob_piKplus",65,2.5,15.5); 
  TH1D *h_prob_pipplus  = new TH1D("h_prob_pipplus","h_prob_pipplus",65,2.5,15.5);

  TH2D *h_likelihood_piKminus = new TH2D("h_likelihood_piKminus","h_likelihood_piKminus",65,2.5,15.5,450,-20.0,430.0);
  TH2D *h_likelihood_pipminus = new TH2D("h_likelihood_pipminus","h_likelihood_pipminus",65,2.5,15.5,450,-20.0,430.0);
  double counter_pipiminus = 0.0;
  double counter_piKminus = 0.0;
  double counter_pipminus = 0.0;
  TH1D *h_prob_pipiminus = new TH1D("h_prob_pipiminus","h_prob_pipiminus",65,2.5,15.5);
  TH1D *h_prob_piKminus  = new TH1D("h_prob_piKminus","h_prob_piKminus",65,2.5,15.5);
  TH1D *h_prob_pipminus  = new TH1D("h_prob_pipminus","h_prob_pipminus",65,2.5,15.5);

  TH2D *h_likelihood_Kpiplus = new TH2D("h_likelihood_Kpiplus","h_likelihood_Kpiplus",65,2.5,15.5,450,-20.0,430.0); // K to pi
  TH2D *h_likelihood_Kpplus = new TH2D("h_likelihood_Kpplus","h_likelihood_Kpplus",65,2.5,15.5,450,-20.0,430.0); // K to p
  double counter_Kpiplus = 0.0;
  double counter_KKplus = 0.0;
  double counter_Kpplus = 0.0;
  TH1D *h_prob_Kpiplus = new TH1D("h_prob_Kpiplus","h_prob_Kpiplus",65,2.5,15.5);
  TH1D *h_prob_KKplus  = new TH1D("h_prob_KKplus","h_prob_KKplus",65,2.5,15.5);
  TH1D *h_prob_Kpplus  = new TH1D("h_prob_Kpplus","h_prob_Kpplus",65,2.5,15.5);

  TH2D *h_likelihood_Kpiminus = new TH2D("h_likelihood_Kpiminus","h_likelihood_Kpiminus",65,2.5,15.5,450,-20.0,430.0);
  TH2D *h_likelihood_Kpminus = new TH2D("h_likelihood_Kpminus","h_likelihood_Kpminus",65,2.5,15.5,450,-20.0,430.0);
  double counter_Kpiminus = 0.0;
  double counter_KKminus = 0.0;
  double counter_Kpminus = 0.0;
  TH1D *h_prob_Kpiminus = new TH1D("h_prob_Kpiminus","h_prob_Kpiminus",65,2.5,15.5);
  TH1D *h_prob_KKminus  = new TH1D("h_prob_KKminus","h_prob_KKminus",65,2.5,15.5);
  TH1D *h_prob_Kpminus  = new TH1D("h_prob_Kpminus","h_prob_Kpminus",65,2.5,15.5);
  
  TH2D *h_likelihood_ppiplus = new TH2D("h_likelihood_ppiplus","h_likelihood_ppiplus",65,2.5,15.5,450,-20.0,430.0); // p to pi
  TH2D *h_likelihood_pKplus = new TH2D("h_likelihood_pKplus","h_likelihood_pKplus",65,2.5,15.5,450,-20.0,430.0); // p to K
  double counter_ppplus = 0.0;
  double counter_ppiplus = 0.0;
  double counter_pKplus = 0.0;
  TH1D *h_prob_ppiplus = new TH1D("h_prob_ppiplus","h_prob_ppiplus",65,2.5,15.5);
  TH1D *h_prob_pKplus  = new TH1D("h_prob_pKplus","h_prob_pKplus",65,2.5,15.5);
  TH1D *h_prob_ppplus  = new TH1D("h_prob_ppplus","h_prob_ppplus",65,2.5,15.5);

  TH2D *h_likelihood_ppiminus = new TH2D("h_likelihood_ppiminus","h_likelihood_ppiminus",65,2.5,15.5,450,-20.0,430.0);
  TH2D *h_likelihood_pKminus = new TH2D("h_likelihood_pKminus","h_likelihood_pKminus",65,2.5,15.5,450,-20.0,430.0);
  double counter_ppminus = 0.0;
  double counter_ppiminus = 0.0;
  double counter_pKminus = 0.0;
  TH1D *h_prob_ppiminus = new TH1D("h_prob_ppiminus","h_prob_ppiminus",65,2.5,15.5);
  TH1D *h_prob_pKminus  = new TH1D("h_prob_pKminus","h_prob_pKminus",65,2.5,15.5);
  TH1D *h_prob_ppminus  = new TH1D("h_prob_ppminus","h_prob_ppminus",65,2.5,15.5);

  Int_t nEvents = fchain->GetEntries();
  for(int i_event = 0; i_event < nEvents; i_event++)
  {
    if(i_event%1000 == 0) cout<<"====== Processing "<< i_event <<" evnets; ====="<< endl;
    
    fchain->GetEntry(i_event);
    double mom = sqrt(px*px+py*py+pz*pz);

    if(pid == 211) 
    {
      h_likelihood_piKplus->Fill(mom,Lpion-LKaon);
      h_likelihood_pipplus->Fill(mom,Lpion-Lproton);
      if( (Lpion-LKaon) > 0.0 && (Lpion-Lproton) > 0.0) counter_pipiplus++;
      if( (LKaon-Lpion) > 0.0 && (LKaon-Lproton) > 0.0) counter_piKplus++;
      if( (Lproton-Lpion) > 0.0 && (Lproton-LKaon) > 0.0) counter_pipplus++;
    }
    else if(pid == -211)
    {
      h_likelihood_piKminus->Fill(mom,Lpion-LKaon);
      h_likelihood_pipminus->Fill(mom,Lpion-Lproton);
      if( (Lpion-LKaon) > 0.0 && (Lpion-Lproton) > 0.0) counter_pipiminus++;
      if( (LKaon-Lpion) > 0.0 && (LKaon-Lproton) > 0.0) counter_piKminus++;
      if( (Lproton-Lpion) > 0.0 && (Lproton-LKaon) > 0.0) counter_pipminus++;
    }
    else if(pid == 321) 
    {
      h_likelihood_Kpiplus->Fill(mom,LKaon-Lpion);
      h_likelihood_Kpplus->Fill(mom,LKaon-Lproton);
      if( (Lpion-LKaon) > 0.0 && (Lpion-Lproton) > 0.0) counter_Kpiplus++;
      if( (LKaon-Lpion) > 0.0 && (LKaon-Lproton) > 0.0) counter_KKplus++;
      if( (Lproton-Lpion) > 0.0 && (Lproton-LKaon) > 0.0) counter_Kpplus++;
    }
    else if(pid == -321)
    {
      h_likelihood_Kpiminus->Fill(mom,LKaon-Lpion);
      h_likelihood_Kpminus->Fill(mom,LKaon-Lproton);
      if( (Lpion-LKaon) > 0.0 && (Lpion-Lproton) > 0.0) counter_Kpiminus++;
      if( (LKaon-Lpion) > 0.0 && (LKaon-Lproton) > 0.0) counter_KKminus++;
      if( (Lproton-Lpion) > 0.0 && (Lproton-LKaon) > 0.0) counter_Kpminus++;
    }
    else if(pid == 2212) 
    {
      h_likelihood_ppiplus->Fill(mom,Lproton-Lpion);
      h_likelihood_pKplus->Fill(mom,Lproton-LKaon);
      if( (Lpion-LKaon) > 0.0 && (Lpion-Lproton) > 0.0) counter_ppiplus++;
      if( (LKaon-Lpion) > 0.0 && (LKaon-Lproton) > 0.0) counter_pKplus++;
      if( (Lproton-Lpion) > 0.0 && (Lproton-LKaon) > 0.0) counter_ppplus++;
    }
    else if(pid == -321)
    {
      h_likelihood_ppiminus->Fill(mom,Lproton-Lpion);
      h_likelihood_pKminus->Fill(mom,Lproton-LKaon);
      if( (Lpion-LKaon) > 0.0 && (Lpion-Lproton) > 0.0) counter_ppiminus++;
      if( (LKaon-Lpion) > 0.0 && (LKaon-Lproton) > 0.0) counter_pKminus++;
      if( (Lproton-Lpion) > 0.0 && (Lproton-LKaon) > 0.0) counter_ppminus++;
    }
  }

  double sum_piplus = counter_pipiplus + counter_piKplus + counter_pipplus;
  h_prob_pipiplus->SetBinContent(h_prob_pipiplus->FindBin(5.0),counter_pipiplus/sum_piplus);
  h_prob_piKplus->SetBinContent(h_prob_pipiplus->FindBin(5.0),counter_piKplus/sum_piplus);
  h_prob_pipplus->SetBinContent(h_prob_pipiplus->FindBin(5.0),counter_pipplus/sum_piplus);

  double sum_piminus = counter_pipiminus + counter_piKminus + counter_pipiminus;
  h_prob_pipiminus->SetBinContent(h_prob_pipiplus->FindBin(5.0),counter_pipiminus/sum_piminus);
  h_prob_piKminus->SetBinContent(h_prob_pipiplus->FindBin(5.0),counter_piKminus/sum_piminus);
  h_prob_pipminus->SetBinContent(h_prob_pipiplus->FindBin(5.0),counter_pipminus/sum_piminus);

  double sum_Kplus = counter_Kpiplus + counter_KKplus + counter_Kpplus;
  h_prob_Kpiplus->SetBinContent(h_prob_pipiplus->FindBin(5.0),counter_Kpiplus/sum_Kplus);
  h_prob_KKplus->SetBinContent(h_prob_pipiplus->FindBin(5.0),counter_KKplus/sum_Kplus);
  h_prob_Kpplus->SetBinContent(h_prob_pipiplus->FindBin(5.0),counter_Kpplus/sum_Kplus);

  double sum_Kminus = counter_Kpiminus + counter_KKminus + counter_Kpiminus;
  h_prob_Kpiminus->SetBinContent(h_prob_pipiplus->FindBin(5.0),counter_Kpiminus/sum_Kminus);
  h_prob_KKminus->SetBinContent(h_prob_pipiplus->FindBin(5.0),counter_KKminus/sum_Kminus);
  h_prob_Kpminus->SetBinContent(h_prob_pipiplus->FindBin(5.0),counter_Kpminus/sum_Kminus);

  double sum_pplus = counter_ppiplus + counter_pKplus + counter_ppplus;
  h_prob_ppiplus->SetBinContent(h_prob_pipiplus->FindBin(5.0),counter_ppiplus/sum_pplus);
  h_prob_pKplus->SetBinContent(h_prob_pipiplus->FindBin(5.0),counter_pKplus/sum_pplus);
  h_prob_ppplus->SetBinContent(h_prob_pipiplus->FindBin(5.0),counter_ppplus/sum_pplus);

  double sum_pminus = counter_ppiminus + counter_pKminus + counter_ppiminus;
  h_prob_ppiminus->SetBinContent(h_prob_pipiplus->FindBin(5.0),counter_ppiminus/sum_pminus);
  h_prob_pKminus->SetBinContent(h_prob_pipiplus->FindBin(5.0),counter_pKminus/sum_pminus);
  h_prob_ppminus->SetBinContent(h_prob_pipiplus->FindBin(5.0),counter_ppminus/sum_pminus);
  
  TFile *Fill_OutPut = new TFile("/work/eic/xusun/output/modular_rich/test/LogLikelihood.root","RECREATE");

  h_likelihood_piKplus->Write();
  h_likelihood_piKminus->Write();
  h_likelihood_pipplus->Write();
  h_likelihood_pipminus->Write();

  h_likelihood_Kpiplus->Write();
  h_likelihood_Kpiminus->Write();
  h_likelihood_Kpplus->Write();
  h_likelihood_Kpminus->Write();
  
  h_likelihood_ppiplus->Write();
  h_likelihood_ppiminus->Write();
  h_likelihood_pKplus->Write();
  h_likelihood_pKminus->Write();

  h_prob_pipiplus->Write();
  h_prob_piKplus->Write();
  h_prob_pipplus->Write();

  h_prob_pipiminus->Write();
  h_prob_piKminus->Write();
  h_prob_pipminus->Write();

  h_prob_Kpiplus->Write();
  h_prob_KKplus->Write();
  h_prob_Kpplus->Write();

  h_prob_Kpiminus->Write();
  h_prob_KKminus->Write();
  h_prob_Kpminus->Write();

  h_prob_ppiplus->Write();
  h_prob_pKplus->Write();
  h_prob_ppplus->Write();

  h_prob_ppiminus->Write();
  h_prob_pKminus->Write();
  h_prob_ppminus->Write();
  
  Fill_OutPut->Close();
}
