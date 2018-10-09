#include "TStyle.h"
#include "TString.h"
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TH2D.h"
#include "TH1D.h"
#include "string.h"
#include "../include/mRICH.h"

using namespace std;

int get_indexMomentumP(double momentum)
{
  const double delta_p = (mRICH::mMomP_stop-mRICH::mMomP_start)/mRICH::mNumOfIndexMomentumP; 

  int index = -1;
  if(momentum > mRICH::mMomP_start && momentum < mRICH::mMomP_stop)
  {
    for(int i_index = 0; i_index < mRICH::mNumOfIndexMomentumP; ++i_index)
    {
      if(momentum > mRICH::mMomP_start+i_index*delta_p && momentum <= mRICH::mMomP_start+(i_index+1)*delta_p)
      {
	index = i_index;
      }
    }
  }

  return index;
}

void PID_mRICH(const string inputfile = "/work/eic/xusun/output/likelihood/PID_Likelihood_May10_2018.root") 
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
  double counter_pipluspiplus[mRICH::mNumOfIndexMomentumP];
  double counter_piplusKplus[mRICH::mNumOfIndexMomentumP];
  double counter_piplusproton[mRICH::mNumOfIndexMomentumP];
  double sum_piplus[mRICH::mNumOfIndexMomentumP];
  TH1D *h_prob_pipluspiplus = new TH1D("h_prob_pipluspiplus","h_prob_pipluspiplus",mRICH::mNumOfIndexMomentumP,mRICH::mMomP_start,mRICH::mMomP_stop);
  TH1D *h_prob_piplusKplus  = new TH1D("h_prob_piplusKplus","h_prob_piplusKplus",mRICH::mNumOfIndexMomentumP,mRICH::mMomP_start,mRICH::mMomP_stop);
  TH1D *h_prob_piplusproton  = new TH1D("h_prob_piplusproton","h_prob_piplusproton",mRICH::mNumOfIndexMomentumP,mRICH::mMomP_start,mRICH::mMomP_stop);

  TH2D *h_likelihood_piminusKminus = new TH2D("h_likelihood_piminusKminus","h_likelihood_piminusKminus",65,2.5,15.5,450,-20.0,430.0);
  TH2D *h_likelihood_piminusantiproton = new TH2D("h_likelihood_piminusantiproton","h_likelihood_piminusantiproton",65,2.5,15.5,450,-20.0,430.0);
  double counter_piminuspiminus[mRICH::mNumOfIndexMomentumP];
  double counter_piminusKminus[mRICH::mNumOfIndexMomentumP];
  double counter_piminusantiproton[mRICH::mNumOfIndexMomentumP];
  double sum_piminus[mRICH::mNumOfIndexMomentumP];
  TH1D *h_prob_piminuspiminus = new TH1D("h_prob_piminuspiminus","h_prob_piminuspiminus",mRICH::mNumOfIndexMomentumP,mRICH::mMomP_start,mRICH::mMomP_stop);
  TH1D *h_prob_piminusKminus  = new TH1D("h_prob_piminusKminus","h_prob_piminusKminus",mRICH::mNumOfIndexMomentumP,mRICH::mMomP_start,mRICH::mMomP_stop);
  TH1D *h_prob_piminusantiproton  = new TH1D("h_prob_piminusantiproton","h_prob_piminusantiproton",mRICH::mNumOfIndexMomentumP,mRICH::mMomP_start,mRICH::mMomP_stop);

  TH2D *h_likelihood_Kpluspiplus = new TH2D("h_likelihood_Kpluspiplus","h_likelihood_Kpluspiplus",65,2.5,15.5,450,-20.0,430.0); // K to pi
  TH2D *h_likelihood_Kplusproton = new TH2D("h_likelihood_Kplusproton","h_likelihood_Kplusproton",65,2.5,15.5,450,-20.0,430.0); // K to p
  double counter_Kpluspiplus[mRICH::mNumOfIndexMomentumP];
  double counter_KplusKplus[mRICH::mNumOfIndexMomentumP];
  double counter_Kplusproton[mRICH::mNumOfIndexMomentumP];
  double sum_Kplus[mRICH::mNumOfIndexMomentumP];
  TH1D *h_prob_Kpluspiplus = new TH1D("h_prob_Kpluspiplus","h_prob_Kpluspiplus",mRICH::mNumOfIndexMomentumP,mRICH::mMomP_start,mRICH::mMomP_stop);
  TH1D *h_prob_KplusKplus  = new TH1D("h_prob_KplusKplus","h_prob_KplusKplus",mRICH::mNumOfIndexMomentumP,mRICH::mMomP_start,mRICH::mMomP_stop);
  TH1D *h_prob_Kplusproton  = new TH1D("h_prob_Kplusproton","h_prob_Kplusproton",mRICH::mNumOfIndexMomentumP,mRICH::mMomP_start,mRICH::mMomP_stop);

  TH2D *h_likelihood_Kminuspiminus = new TH2D("h_likelihood_Kminuspiminus","h_likelihood_Kminuspiminus",65,2.5,15.5,450,-20.0,430.0);
  TH2D *h_likelihood_Kminusantiproton = new TH2D("h_likelihood_Kminusantiproton","h_likelihood_Kminusantiproton",65,2.5,15.5,450,-20.0,430.0);
  double counter_Kminuspiminus[mRICH::mNumOfIndexMomentumP];
  double counter_KminusKminus[mRICH::mNumOfIndexMomentumP];
  double counter_Kminusantiproton[mRICH::mNumOfIndexMomentumP];
  double sum_Kminus[mRICH::mNumOfIndexMomentumP];
  TH1D *h_prob_Kminuspiminus = new TH1D("h_prob_Kminuspiminus","h_prob_Kminuspiminus",mRICH::mNumOfIndexMomentumP,mRICH::mMomP_start,mRICH::mMomP_stop);
  TH1D *h_prob_KminusKminus  = new TH1D("h_prob_KminusKminus","h_prob_KminusKminus",mRICH::mNumOfIndexMomentumP,mRICH::mMomP_start,mRICH::mMomP_stop);
  TH1D *h_prob_Kminusantiproton  = new TH1D("h_prob_Kminusantiproton","h_prob_Kminusantiproton",mRICH::mNumOfIndexMomentumP,mRICH::mMomP_start,mRICH::mMomP_stop);
  
  TH2D *h_likelihood_protonpiplus = new TH2D("h_likelihood_protonpiplus","h_likelihood_protonpiplus",65,2.5,15.5,450,-20.0,430.0); // p to pi
  TH2D *h_likelihood_protonKplus = new TH2D("h_likelihood_protonKplus","h_likelihood_protonKplus",65,2.5,15.5,450,-20.0,430.0); // p to K
  double counter_protonpiplus[mRICH::mNumOfIndexMomentumP];
  double counter_protonKplus[mRICH::mNumOfIndexMomentumP];
  double counter_protonproton[mRICH::mNumOfIndexMomentumP];
  double sum_pplus[mRICH::mNumOfIndexMomentumP];
  TH1D *h_prob_protonpiplus = new TH1D("h_prob_protonpiplus","h_prob_protonpiplus",mRICH::mNumOfIndexMomentumP,mRICH::mMomP_start,mRICH::mMomP_stop);
  TH1D *h_prob_protonKplus  = new TH1D("h_prob_protonKplus","h_prob_protonKplus",mRICH::mNumOfIndexMomentumP,mRICH::mMomP_start,mRICH::mMomP_stop);
  TH1D *h_prob_protonproton  = new TH1D("h_prob_protonproton","h_prob_protonproton",mRICH::mNumOfIndexMomentumP,mRICH::mMomP_start,mRICH::mMomP_stop);

  TH2D *h_likelihood_antiprotonpiminus = new TH2D("h_likelihood_antiprotonpiminus","h_likelihood_antiprotonpiminus",65,2.5,15.5,450,-20.0,430.0);
  TH2D *h_likelihood_antiprotonKminus = new TH2D("h_likelihood_antiprotonKminus","h_likelihood_antiprotonKminus",65,2.5,15.5,450,-20.0,430.0);
  double counter_antiprotonpiminus[mRICH::mNumOfIndexMomentumP];
  double counter_antiprotonKminus[mRICH::mNumOfIndexMomentumP];
  double counter_antiprotonantiproton[mRICH::mNumOfIndexMomentumP];
  double sum_pminus[mRICH::mNumOfIndexMomentumP];
  TH1D *h_prob_antiprotonpiminus = new TH1D("h_prob_antiprotonpiminus","h_prob_antiprotonpiminus",mRICH::mNumOfIndexMomentumP,mRICH::mMomP_start,mRICH::mMomP_stop);
  TH1D *h_prob_antiprotonKminus  = new TH1D("h_prob_antiprotonKminus","h_prob_antiprotonKminus",mRICH::mNumOfIndexMomentumP,mRICH::mMomP_start,mRICH::mMomP_stop);
  TH1D *h_prob_antiprotonantiproton  = new TH1D("h_prob_antiprotonantiproton","h_prob_antiprotonantiproton",mRICH::mNumOfIndexMomentumP,mRICH::mMomP_start,mRICH::mMomP_stop);

  for(int i_mom = 0; i_mom < mRICH::mNumOfIndexMomentumP; ++i_mom)
  {
    counter_pipluspiplus[i_mom] = 0.0;
    counter_piplusKplus[i_mom] = 0.0;
    counter_piplusproton[i_mom] = 0.0;
    sum_piplus[i_mom] = 0.0;

    counter_piminuspiminus[i_mom] = 0.0;
    counter_piminusKminus[i_mom] = 0.0;
    counter_piminusantiproton[i_mom] = 0.0;
    sum_piminus[i_mom] = 0.0;

    counter_Kpluspiplus[i_mom] = 0.0;
    counter_KplusKplus[i_mom] = 0.0;
    counter_Kplusproton[i_mom] = 0.0;
    sum_Kplus[i_mom] = 0.0;

    counter_Kminuspiminus[i_mom] = 0.0;
    counter_KminusKminus[i_mom] = 0.0;
    counter_Kminusantiproton[i_mom] = 0.0;
    sum_Kminus[i_mom] = 0.0;

    counter_protonpiplus[i_mom] = 0.0;
    counter_protonKplus[i_mom] = 0.0;
    counter_protonproton[i_mom] = 0.0;
    sum_pplus[i_mom] = 0.0;

    counter_antiprotonpiminus[i_mom] = 0.0;
    counter_antiprotonKminus[i_mom] = 0.0;
    counter_antiprotonantiproton[i_mom] = 0.0;
    sum_pminus[i_mom] = 0.0;
  }

  Int_t nEvents = fchain->GetEntries();
  for(int i_event = 0; i_event < nEvents; i_event++)
  {
    if(i_event%1000 == 0) cout<<"====== Processing "<< i_event <<" evnets; ====="<< endl;
    
    fchain->GetEntry(i_event);
    double mom = sqrt(px*px+py*py+pz*pz);
    int index = get_indexMomentumP(mom);
    if(index < 0) continue;

    if(pid == 211) 
    {
      h_likelihood_piplusKplus->Fill(mom,Lpion-LKaon);
      h_likelihood_piplusproton->Fill(mom,Lpion-Lproton);
      sum_piplus[index]++;
      if( (Lpion-LKaon) > 0.0 && (Lpion-Lproton) > 0.0) counter_pipluspiplus[index]++;
      if( (LKaon-Lpion) > 0.0 && (LKaon-Lproton) > 0.0) counter_piplusKplus[index]++;
      if( (Lproton-Lpion) > 0.0 && (Lproton-LKaon) > 0.0) counter_piplusproton[index]++;
    }
    else if(pid == -211)
    {
      h_likelihood_piminusKminus->Fill(mom,Lpion-LKaon);
      h_likelihood_piminusantiproton->Fill(mom,Lpion-Lproton);
      sum_piminus[index]++;
      if( (Lpion-LKaon) > 0.0 && (Lpion-Lproton) > 0.0) counter_piminuspiminus[index]++;
      if( (LKaon-Lpion) > 0.0 && (LKaon-Lproton) > 0.0) counter_piminusKminus[index]++;
      if( (Lproton-Lpion) > 0.0 && (Lproton-LKaon) > 0.0) counter_piminusantiproton[index]++;
    }
    else if(pid == 321) 
    {
      h_likelihood_Kpluspiplus->Fill(mom,LKaon-Lpion);
      h_likelihood_Kplusproton->Fill(mom,LKaon-Lproton);
      sum_Kplus[index]++;
      if( (Lpion-LKaon) > 0.0 && (Lpion-Lproton) > 0.0) counter_Kpluspiplus[index]++;
      if( (LKaon-Lpion) > 0.0 && (LKaon-Lproton) > 0.0) counter_KplusKplus[index]++;
      if( (Lproton-Lpion) > 0.0 && (Lproton-LKaon) > 0.0) counter_Kplusproton[index]++;
    }
    else if(pid == -321)
    {
      h_likelihood_Kminuspiminus->Fill(mom,LKaon-Lpion);
      h_likelihood_Kminusantiproton->Fill(mom,LKaon-Lproton);
      sum_Kminus[index]++;
      if( (Lpion-LKaon) > 0.0 && (Lpion-Lproton) > 0.0) counter_Kminuspiminus[index]++;
      if( (LKaon-Lpion) > 0.0 && (LKaon-Lproton) > 0.0) counter_KminusKminus[index]++;
      if( (Lproton-Lpion) > 0.0 && (Lproton-LKaon) > 0.0) counter_Kminusantiproton[index]++;
    }
    else if(pid == 2212) 
    {
      h_likelihood_protonpiplus->Fill(mom,Lproton-Lpion);
      h_likelihood_protonKplus->Fill(mom,Lproton-LKaon);
      sum_pplus[index]++;
      if( (Lpion-LKaon) > 0.0 && (Lpion-Lproton) > 0.0) counter_protonpiplus[index]++;
      if( (LKaon-Lpion) > 0.0 && (LKaon-Lproton) > 0.0) counter_protonKplus[index]++;
      if( (Lproton-Lpion) > 0.0 && (Lproton-LKaon) > 0.0) counter_protonproton[index]++;
    }
    else if(pid == -2212)
    {
      h_likelihood_antiprotonpiminus->Fill(mom,Lproton-Lpion);
      h_likelihood_antiprotonKminus->Fill(mom,Lproton-LKaon);
      sum_pminus[index]++;
      if( (Lpion-LKaon) > 0.0 && (Lpion-Lproton) > 0.0) counter_antiprotonpiminus[index]++;
      if( (LKaon-Lpion) > 0.0 && (LKaon-Lproton) > 0.0) counter_antiprotonKminus[index]++;
      if( (Lproton-Lpion) > 0.0 && (Lproton-LKaon) > 0.0) counter_antiprotonantiproton[index]++;
    }
  }

  for(int i_mom = 0; i_mom < mRICH::mNumOfIndexMomentumP; ++i_mom)
  {
    if(sum_piplus[i_mom]> 0.0)
    {
      h_prob_pipluspiplus->SetBinContent(i_mom+1,counter_pipluspiplus[i_mom]/sum_piplus[i_mom]);
      h_prob_piplusKplus->SetBinContent(i_mom+1,counter_piplusKplus[i_mom]/sum_piplus[i_mom]);
      h_prob_piplusproton->SetBinContent(0,5.0,counter_piplusproton[i_mom]/sum_piplus[i_mom]);
    }

    if(sum_piminus[i_mom]> 0.0)
    {
      h_prob_piminuspiminus->SetBinContent(i_mom+1,counter_piminuspiminus[i_mom]/sum_piminus[i_mom]);
      h_prob_piminusKminus->SetBinContent(i_mom+1,counter_piminusKminus[i_mom]/sum_piminus[i_mom]);
      h_prob_piminusantiproton->SetBinContent(i_mom+1,counter_piminusantiproton[i_mom]/sum_piminus[i_mom]);
    }

    if(sum_Kplus[i_mom]> 0.0)
    {
      h_prob_Kpluspiplus->SetBinContent(i_mom+1,counter_Kpluspiplus[i_mom]/sum_Kplus[i_mom]);
      h_prob_KplusKplus->SetBinContent(i_mom+1,counter_KplusKplus[i_mom]/sum_Kplus[i_mom]);
      h_prob_Kplusproton->SetBinContent(i_mom+1,counter_Kplusproton[i_mom]/sum_Kplus[i_mom]);
    }

    if(sum_Kminus[i_mom]> 0.0)
    {
      h_prob_Kminuspiminus->SetBinContent(i_mom+1,counter_Kminuspiminus[i_mom]/sum_Kminus[i_mom]);
      h_prob_KminusKminus->SetBinContent(i_mom+1,counter_KminusKminus[i_mom]/sum_Kminus[i_mom]);
      h_prob_Kminusantiproton->SetBinContent(i_mom+1,counter_Kminusantiproton[i_mom]/sum_Kminus[i_mom]);
    }

    if(sum_pplus[i_mom]> 0.0)
    {
      h_prob_protonpiplus->SetBinContent(i_mom+1,counter_protonpiplus[i_mom]/sum_pplus[i_mom]);
      h_prob_protonKplus->SetBinContent(i_mom+1,counter_protonKplus[i_mom]/sum_pplus[i_mom]);
      h_prob_protonproton->SetBinContent(i_mom+1,counter_protonproton[i_mom]/sum_pplus[i_mom]);
    }

    if(sum_pminus[i_mom]> 0.0)
    {
      h_prob_antiprotonpiminus->SetBinContent(i_mom+1,counter_antiprotonpiminus[i_mom]/sum_pminus[i_mom]);
      h_prob_antiprotonKminus->SetBinContent(i_mom+1,counter_antiprotonKminus[i_mom]/sum_pminus[i_mom]);
      h_prob_antiprotonantiproton->SetBinContent(i_mom+1,counter_antiprotonantiproton[i_mom]/sum_pminus[i_mom]);
    }
  }
  
  TFile *Fill_OutPut = new TFile("/work/eic/xusun/output/likelihood/LogLikelihood_May10_2018.root","RECREATE");

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

  h_prob_pipluspiplus->Write();
  h_prob_piplusKplus->Write();
  h_prob_piplusproton->Write();

  h_prob_piminuspiminus->Write();
  h_prob_piminusKminus->Write();
  h_prob_piminusantiproton->Write();

  h_prob_Kpluspiplus->Write();
  h_prob_KplusKplus->Write();
  h_prob_Kplusproton->Write();

  h_prob_Kminuspiminus->Write();
  h_prob_KminusKminus->Write();
  h_prob_Kminusantiproton->Write();

  h_prob_protonpiplus->Write();
  h_prob_protonKplus->Write();
  h_prob_protonproton->Write();

  h_prob_antiprotonpiminus->Write();
  h_prob_antiprotonKminus->Write();
  h_prob_antiprotonantiproton->Write();
  
  Fill_OutPut->Close();
}
