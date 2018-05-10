#include "TStyle.h"
#include "TString.h"
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TGraph.h"
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

  TH2D *h_likelihood_piplusKplus = new TH2D("h_likelihood_piplusKplus","h_likelihood_piplusKplus",65,2.5,15.5,450,-20.0,430.0); // pi to K
  TH2D *h_likelihood_piplusproton = new TH2D("h_likelihood_piplusproton","h_likelihood_piplusproton",65,2.5,15.5,450,-20.0,430.0); // pi to p
  double counter_pipluspiplus = 0.0;
  double counter_piplusKplus = 0.0;
  double counter_piplusproton = 0.0;
  TGraph *g_prob_pipluspiplus = new TGraph();
  TGraph *g_prob_piplusKplus  = new TGraph();
  TGraph *g_prob_piplusproton  = new TGraph();

  TH2D *h_likelihood_piminusKminus = new TH2D("h_likelihood_piminusKminus","h_likelihood_piminusKminus",65,2.5,15.5,450,-20.0,430.0);
  TH2D *h_likelihood_piminusantiproton = new TH2D("h_likelihood_piminusantiproton","h_likelihood_piminusantiproton",65,2.5,15.5,450,-20.0,430.0);
  double counter_piminuspiminus = 0.0;
  double counter_piminusKminus = 0.0;
  double counter_piminusantiproton = 0.0;
  TGraph *g_prob_piminuspiminus = new TGraph();
  TGraph *g_prob_piminusKminus  = new TGraph();
  TGraph *g_prob_piminusantiproton  = new TGraph();

  TH2D *h_likelihood_Kpluspiplus = new TH2D("h_likelihood_Kpluspiplus","h_likelihood_Kpluspiplus",65,2.5,15.5,450,-20.0,430.0); // K to pi
  TH2D *h_likelihood_Kplusproton = new TH2D("h_likelihood_Kplusproton","h_likelihood_Kplusproton",65,2.5,15.5,450,-20.0,430.0); // K to p
  double counter_Kpluspiplus = 0.0;
  double counter_KplusKplus = 0.0;
  double counter_Kplusproton = 0.0;
  TGraph *g_prob_Kpluspiplus = new TGraph();
  TGraph *g_prob_KplusKplus  = new TGraph();
  TGraph *g_prob_Kplusproton  = new TGraph();

  TH2D *h_likelihood_Kminuspiminus = new TH2D("h_likelihood_Kminuspiminus","h_likelihood_Kminuspiminus",65,2.5,15.5,450,-20.0,430.0);
  TH2D *h_likelihood_Kminusantiproton = new TH2D("h_likelihood_Kminusantiproton","h_likelihood_Kminusantiproton",65,2.5,15.5,450,-20.0,430.0);
  double counter_Kminuspiminus = 0.0;
  double counter_KminusKminus = 0.0;
  double counter_Kminusantiproton = 0.0;
  TGraph *g_prob_Kminuspiminus = new TGraph();
  TGraph *g_prob_KminusKminus  = new TGraph();
  TGraph *g_prob_Kminusantiproton  = new TGraph();
  
  TH2D *h_likelihood_protonpiplus = new TH2D("h_likelihood_protonpiplus","h_likelihood_protonpiplus",65,2.5,15.5,450,-20.0,430.0); // p to pi
  TH2D *h_likelihood_protonKplus = new TH2D("h_likelihood_protonKplus","h_likelihood_protonKplus",65,2.5,15.5,450,-20.0,430.0); // p to K
  double counter_protonpiplus = 0.0;
  double counter_protonKplus = 0.0;
  double counter_protonproton = 0.0;
  TGraph *g_prob_protonpiplus = new TGraph();
  TGraph *g_prob_protonKplus  = new TGraph();
  TGraph *g_prob_protonproton  = new TGraph();

  TH2D *h_likelihood_antiprotonpiminus = new TH2D("h_likelihood_antiprotonpiminus","h_likelihood_antiprotonpiminus",65,2.5,15.5,450,-20.0,430.0);
  TH2D *h_likelihood_antiprotonKminus = new TH2D("h_likelihood_antiprotonKminus","h_likelihood_antiprotonKminus",65,2.5,15.5,450,-20.0,430.0);
  double counter_antiprotonpiminus = 0.0;
  double counter_antiprotonKminus = 0.0;
  double counter_antiprotonantiproton = 0.0;
  TGraph *g_prob_antiprotonpiminus = new TGraph();
  TGraph *g_prob_antiprotonKminus  = new TGraph();
  TGraph *g_prob_antiprotonantiproton  = new TGraph();

  Int_t nEvents = fchain->GetEntries();
  for(int i_event = 0; i_event < nEvents; i_event++)
  {
    if(i_event%1000 == 0) cout<<"====== Processing "<< i_event <<" evnets; ====="<< endl;
    
    fchain->GetEntry(i_event);
    double mom = sqrt(px*px+py*py+pz*pz);

    if(pid == 211) 
    {
      h_likelihood_piplusKplus->Fill(mom,Lpion-LKaon);
      h_likelihood_piplusproton->Fill(mom,Lpion-Lproton);
      if( (Lpion-LKaon) > 0.0 && (Lpion-Lproton) > 0.0) counter_pipluspiplus++;
      if( (LKaon-Lpion) > 0.0 && (LKaon-Lproton) > 0.0) counter_piplusKplus++;
      if( (Lproton-Lpion) > 0.0 && (Lproton-LKaon) > 0.0) counter_piplusproton++;
    }
    else if(pid == -211)
    {
      h_likelihood_piminusKminus->Fill(mom,Lpion-LKaon);
      h_likelihood_piminusantiproton->Fill(mom,Lpion-Lproton);
      if( (Lpion-LKaon) > 0.0 && (Lpion-Lproton) > 0.0) counter_piminuspiminus++;
      if( (LKaon-Lpion) > 0.0 && (LKaon-Lproton) > 0.0) counter_piminusKminus++;
      if( (Lproton-Lpion) > 0.0 && (Lproton-LKaon) > 0.0) counter_piminusantiproton++;
    }
    else if(pid == 321) 
    {
      h_likelihood_Kpluspiplus->Fill(mom,LKaon-Lpion);
      h_likelihood_Kplusproton->Fill(mom,LKaon-Lproton);
      if( (Lpion-LKaon) > 0.0 && (Lpion-Lproton) > 0.0) counter_Kpluspiplus++;
      if( (LKaon-Lpion) > 0.0 && (LKaon-Lproton) > 0.0) counter_KplusKplus++;
      if( (Lproton-Lpion) > 0.0 && (Lproton-LKaon) > 0.0) counter_Kplusproton++;
    }
    else if(pid == -321)
    {
      h_likelihood_Kminuspiminus->Fill(mom,LKaon-Lpion);
      h_likelihood_Kminusantiproton->Fill(mom,LKaon-Lproton);
      if( (Lpion-LKaon) > 0.0 && (Lpion-Lproton) > 0.0) counter_Kminuspiminus++;
      if( (LKaon-Lpion) > 0.0 && (LKaon-Lproton) > 0.0) counter_KminusKminus++;
      if( (Lproton-Lpion) > 0.0 && (Lproton-LKaon) > 0.0) counter_Kminusantiproton++;
    }
    else if(pid == 2212) 
    {
      h_likelihood_protonpiplus->Fill(mom,Lproton-Lpion);
      h_likelihood_protonKplus->Fill(mom,Lproton-LKaon);
      if( (Lpion-LKaon) > 0.0 && (Lpion-Lproton) > 0.0) counter_protonpiplus++;
      if( (LKaon-Lpion) > 0.0 && (LKaon-Lproton) > 0.0) counter_protonKplus++;
      if( (Lproton-Lpion) > 0.0 && (Lproton-LKaon) > 0.0) counter_protonproton++;
    }
    else if(pid == -2212)
    {
      h_likelihood_antiprotonpiminus->Fill(mom,Lproton-Lpion);
      h_likelihood_antiprotonKminus->Fill(mom,Lproton-LKaon);
      if( (Lpion-LKaon) > 0.0 && (Lpion-Lproton) > 0.0) counter_antiprotonpiminus++;
      if( (LKaon-Lpion) > 0.0 && (LKaon-Lproton) > 0.0) counter_antiprotonKminus++;
      if( (Lproton-Lpion) > 0.0 && (Lproton-LKaon) > 0.0) counter_antiprotonantiproton++;
    }
  }

  double sum_piplus = counter_pipluspiplus + counter_piplusKplus + counter_piplusproton;
  g_prob_pipluspiplus->SetPoint(0,5.0,counter_pipluspiplus/sum_piplus);
  g_prob_pipluspiplus->SetName("g_prob_pipluspiplus");
  g_prob_piplusKplus->SetPoint(0,5.0,counter_piplusKplus/sum_piplus);
  g_prob_piplusKplus->SetName("g_prob_piplusKplus");
  g_prob_piplusproton->SetPoint(0,5.0,counter_piplusproton/sum_piplus);
  g_prob_piplusproton->SetName("g_prob_piplusproton");

  double sum_piminus = counter_piminuspiminus + counter_piminusKminus + counter_piminusantiproton;
  g_prob_piminuspiminus->SetPoint(0,5.0,counter_piminuspiminus/sum_piminus);
  g_prob_piminuspiminus->SetName("g_prob_piminuspiminus");
  g_prob_piminusKminus->SetPoint(0,5.0,counter_piminusKminus/sum_piminus);
  g_prob_piminusKminus->SetName("g_prob_piminusKminus");
  g_prob_piminusantiproton->SetPoint(0,5.0,counter_piminusantiproton/sum_piminus);
  g_prob_piminusantiproton->SetName("g_prob_piminusantiproton");

  double sum_Kplus = counter_Kpluspiplus + counter_KplusKplus + counter_Kplusproton;
  g_prob_Kpluspiplus->SetPoint(0,5.0,counter_Kpluspiplus/sum_Kplus);
  g_prob_Kpluspiplus->SetName("g_prob_Kpluspiplus");
  g_prob_KplusKplus->SetPoint(0,5.0,counter_KplusKplus/sum_Kplus);
  g_prob_KplusKplus->SetName("g_prob_KplusKplus");
  g_prob_Kplusproton->SetPoint(0,5.0,counter_Kplusproton/sum_Kplus);
  g_prob_Kplusproton->SetName("g_prob_Kplusproton");

  double sum_Kminus = counter_Kminuspiminus + counter_KminusKminus + counter_Kminuspiminus;
  g_prob_Kminuspiminus->SetPoint(0,5.0,counter_Kminuspiminus/sum_Kminus);
  g_prob_Kminuspiminus->SetName("g_prob_Kminuspiminus");
  g_prob_KminusKminus->SetPoint(0,5.0,counter_KminusKminus/sum_Kminus);
  g_prob_KminusKminus->SetName("g_prob_KminusKminus");
  g_prob_Kminusantiproton->SetPoint(0,5.0,counter_Kminusantiproton/sum_Kminus);
  g_prob_Kminusantiproton->SetName("g_prob_Kminusantiproton");

  double sum_pplus = counter_protonpiplus + counter_protonKplus + counter_protonproton;
  g_prob_protonpiplus->SetPoint(0,5.0,counter_protonpiplus/sum_pplus);
  g_prob_protonpiplus->SetName("g_prob_protonpiplus");
  g_prob_protonKplus->SetPoint(0,5.0,counter_protonKplus/sum_pplus);
  g_prob_protonKplus->SetName("g_prob_protonKplus");
  g_prob_protonproton->SetPoint(0,5.0,counter_protonproton/sum_pplus);
  g_prob_protonproton->SetName("g_prob_protonproton");

  double sum_pminus = counter_antiprotonpiminus + counter_antiprotonKminus + counter_antiprotonantiproton;
  g_prob_antiprotonpiminus->SetPoint(0,5.0,counter_antiprotonpiminus/sum_pminus);
  g_prob_antiprotonpiminus->SetName("g_prob_antiprotonpiminus");
  g_prob_antiprotonKminus->SetPoint(0,5.0,counter_antiprotonKminus/sum_pminus);
  g_prob_antiprotonKminus->SetName("g_prob_antiprotonKminus");
  g_prob_antiprotonantiproton->SetPoint(0,5.0,counter_antiprotonantiproton/sum_pminus);
  g_prob_antiprotonantiproton->SetName("g_prob_antiprotonantiproton");
  
  TFile *Fill_OutPut = new TFile("/work/eic/xusun/output/modular_rich/test/LogLikelihood.root","RECREATE");

  h_likelihood_piplusKplus->Write();
  h_likelihood_piminusKminus->Write();
  h_likelihood_piplusproton->Write();
  h_likelihood_piminusantiproton->Write();

  h_likelihood_Kpluspiplus->Write();
  h_likelihood_Kminuspiminus->Write();
  h_likelihood_Kplusproton->Write();
  h_likelihood_Kminusantiproton->Write();
  
  h_likelihood_protonpiplus->Write();
  h_likelihood_antiprotonpiminus->Write();
  h_likelihood_protonKplus->Write();
  h_likelihood_antiprotonKminus->Write();

  g_prob_pipluspiplus->Write();
  g_prob_piplusKplus->Write();
  g_prob_piplusproton->Write();

  g_prob_piminuspiminus->Write();
  g_prob_piminusKminus->Write();
  g_prob_piminusantiproton->Write();

  g_prob_Kpluspiplus->Write();
  g_prob_KplusKplus->Write();
  g_prob_Kplusproton->Write();

  g_prob_Kminuspiminus->Write();
  g_prob_KminusKminus->Write();
  g_prob_Kminusantiproton->Write();

  g_prob_protonpiplus->Write();
  g_prob_protonKplus->Write();
  g_prob_protonproton->Write();

  g_prob_antiprotonpiminus->Write();
  g_prob_antiprotonKminus->Write();
  g_prob_antiprotonantiproton->Write();
  
  Fill_OutPut->Close();
}
