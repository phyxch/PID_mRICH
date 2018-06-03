#include "TFile.h"
#include "TString.h"
#include "TTree.h"
#include "TH2D.h"
#include "TCanvas.h"
#include <cmath>
#include <cassert>
#include <string>
#include <iostream>
#include <utility>
#include <functional>
#include <stdexcept>
#include "draw.h"
#include "../include/mRICH.h"
using namespace std;

std::pair<int,std::string> get_particle(int pid)
{
  if(pid ==  211)   return std::pair<int,std::string>(pid,"piplus");
  if(pid ==  321)   return std::pair<int,std::string>(pid,"Kplus");
  if(pid ==  2212)  return std::pair<int,std::string>(pid,"proton");
  if(pid == -211)   return std::pair<int,std::string>(pid,"piminus");
  if(pid == -321)   return std::pair<int,std::string>(pid,"Kminus");
  if(pid == -2212)  return std::pair<int,std::string>(pid,"antiproton");
  else return std::pair<int,std::string>(-1,"undefined");
}

std::pair<std::string,std::string> get_misIdentifiedParticle(int pid)
{
  if(pid ==  211)   return std::pair<std::string,std::string>("Kplus","proton");
  if(pid ==  321)   return std::pair<std::string,std::string>("piplus","proton");
  if(pid ==  2212)  return std::pair<std::string,std::string>("piplus","Kplus");
  if(pid == -211)   return std::pair<std::string,std::string>("Kminus","antiproton");
  if(pid == -321)   return std::pair<std::string,std::string>("piminus","antiproton");
  if(pid == -2212)  return std::pair<std::string,std::string>("piminus","Kminus");
  else return std::pair<std::string,std::string>("undefined","undefined");
}

void plotQA_NSigma(int pid = 211, int rank = 0)
{
  string date = "May23_2018";
  // string input_likelihood = Form("/work/eic/xusun/output/probability/3mm/PID_prob_%s.root",date.c_str());
  string input_likelihood = Form("/work/eic/xusun/output/probability/PID_prob_%s.root",date.c_str());
  cout << "read in file: " << input_likelihood.c_str() << endl;
  TFile *File_InPut_Likelihood = TFile::Open(input_likelihood.c_str());
  assert(File_InPut_Likelihood);

  auto identifiedParticle = get_particle(pid);
  std::pair<std::string,std::string> misIdentifiedParticle = get_misIdentifiedParticle(pid);
  const int index_vx = 0;
  const int index_vy = 0;
  const int index_theta = 2;
  const int index_phi = 3;

  string key_likelihood;
  if(rank == 0) key_likelihood = Form("h_mLikelihoodDiff_%s_%s_vx_%d_vy_%d_theta_%d_phi_%d",identifiedParticle.second.c_str(),misIdentifiedParticle.first.c_str(),index_vx,index_vy,index_theta,index_phi);
  if(rank == 1) key_likelihood = Form("h_mLikelihoodDiff_%s_%s_vx_%d_vy_%d_theta_%d_phi_%d",identifiedParticle.second.c_str(),misIdentifiedParticle.second.c_str(),index_vx,index_vy,index_theta,index_phi);

  TH2D *h_mLikelihoodDiff = (TH2D*)File_InPut_Likelihood->Get(key_likelihood.c_str());
  assert(h_mLikelihoodDiff);

  TCanvas *c_likelihood = new TCanvas("c_likelihood","c_likelihood",10,10,800,800);
  c_likelihood->cd(1);
  c_likelihood->cd(1)->SetLeftMargin(0.15);
  c_likelihood->cd(1)->SetBottomMargin(0.15);
  c_likelihood->cd(1)->SetGrid(0,0);
  c_likelihood->cd(1)->SetTicks(1,1);

  string title = Form("log(%s)-log(%s)",identifiedParticle.second.c_str(),misIdentifiedParticle.first.c_str());
  h_mLikelihoodDiff->GetXaxis()->SetTitle("momentum (GeV/c)");
  h_mLikelihoodDiff->GetXaxis()->CenterTitle();
  h_mLikelihoodDiff->GetXaxis()->SetTitleOffset(1.25);
  h_mLikelihoodDiff->GetXaxis()->SetNdivisions(505);

  h_mLikelihoodDiff->GetYaxis()->SetTitle(title.c_str());
  h_mLikelihoodDiff->GetYaxis()->CenterTitle();
  h_mLikelihoodDiff->GetYaxis()->SetTitleOffset(1.25);
  h_mLikelihoodDiff->GetYaxis()->SetNdivisions(505);
  h_mLikelihoodDiff->Draw("colz");

  string output_likelihood;
  if(rank == 0) output_likelihood = Form("../figures/nSigma/c_mLikelihoodDiff_%s_%s_vx_%d_vy_%d_theta_%d_phi_%d.pdf",identifiedParticle.second.c_str(),misIdentifiedParticle.first.c_str(),index_vx,index_vy,index_theta,index_phi);
  if(rank == 1) output_likelihood = Form("../figures/nSigma/c_mLikelihoodDiff_%s_%s_vx_%d_vy_%d_theta_%d_phi_%d.pdf",identifiedParticle.second.c_str(),misIdentifiedParticle.second.c_str(),index_vx,index_vy,index_theta,index_phi);
  c_likelihood->SaveAs(output_likelihood.c_str());

  int numBin = h_mLikelihoodDiff->GetNbinsX();
  const int NumOfMomBin = numBin;
  TH1D *h_mLikelihoodDiff_Proj[NumOfMomBin];
  // cout << "NumOfMomBin = " << NumOfMomBin << endl;
  for(int i_mom = 0; i_mom < NumOfMomBin; ++i_mom)
  {
    string HistName = Form("h_mLikelihoodDiff_Proj_%d",i_mom);
    h_mLikelihoodDiff_Proj[i_mom] = (TH1D*)h_mLikelihoodDiff->ProjectionY(HistName.c_str(),i_mom+1,i_mom+1);
  }

  int nPad_x = 4;
  int nPad_y = 2;
  TCanvas *c_proj = new TCanvas("c_proj","c_proj",10,10,nPad_x*400,nPad_y*400);
  c_proj->Divide(nPad_x,nPad_y);
  for(int i_pad = 0; i_pad < NumOfMomBin; ++i_pad)
  {
    c_proj->cd(i_pad+1);
    c_proj->cd(i_pad+1)->SetLeftMargin(0.15);
    c_proj->cd(i_pad+1)->SetBottomMargin(0.15);
    c_proj->cd(i_pad+1)->SetGrid(0,0);
    c_proj->cd(i_pad+1)->SetTicks(1,1);
    h_mLikelihoodDiff_Proj[i_pad]->SetMarkerStyle(24);
    h_mLikelihoodDiff_Proj[i_pad]->SetMarkerSize(1.4);
    h_mLikelihoodDiff_Proj[i_pad]->Draw("pE same");

    double mean  = h_mLikelihoodDiff_Proj[i_pad]->GetMean();
    double width = h_mLikelihoodDiff_Proj[i_pad]->GetStdDev();

    TF1 *f_gaus = new TF1("f_gaus","gaus",-500,500);
    f_gaus->SetParameter(0,100);
    f_gaus->SetParameter(1,mean);
    f_gaus->SetParameter(2,width);
    f_gaus->SetRange(mean-5.0*width,mean+5.0*width);
    h_mLikelihoodDiff_Proj[i_pad]->Fit(f_gaus,"NR");
    f_gaus->SetLineColor(2);
    f_gaus->SetLineWidth(4);
    f_gaus->SetLineStyle(2);
    f_gaus->Draw("l same");

    string leg_mom = Form("[%1.1f,%1.1f GeV/c]",mRICH::mBin_MomP[i_pad]-mRICH::mDelta_MomP,mRICH::mBin_MomP[i_pad]+mRICH::mDelta_MomP);
    plotTopLegend(leg_mom.c_str(),0.45,0.6,0.05,1,0.0,42,1,1);

    double mean_diff = f_gaus->GetParameter(1);
    double width_diff = f_gaus->GetParameter(2);
    // string leg_likelihood = Form("#Deltaln(likelihood) = %3.1f",mean_diff);
    string leg_likelihood = Form("#Deltaln(likelihood) = %3.1f",mean);
    plotTopLegend(leg_likelihood.c_str(),0.45,0.5,0.05,1,0.0,42,1,1);

    // double sigma_diff = TMath::Sqrt(2.0*mean_diff);
    // double err_diff = width_diff/TMath::Sqrt(2.0*mean_diff);
    // string leg_nSigma = Form("N_{#sigma} = %3.1f",TMath::Sqrt(2.0*mean_diff));
    double sigma_diff = TMath::Sqrt(2.0*mean);
    double err_diff = width/TMath::Sqrt(2.0*mean);
    string leg_nSigma = Form("N_{#sigma} = %3.1f",TMath::Sqrt(2.0*mean));
    plotTopLegend(leg_nSigma.c_str(),0.45,0.4,0.05,1,0.0,42,1,1);
  }

  string output_proj;
  if(rank == 0) output_proj = Form("../figures/nSigma/c_mProj_%s_%s_vx_%d_vy_%d_theta_%d_phi_%d.pdf",identifiedParticle.second.c_str(),misIdentifiedParticle.first.c_str(),index_vx,index_vy,index_theta,index_phi);
  if(rank == 1) output_proj = Form("../figures/nSigma/c_mProj_%s_%s_vx_%d_vy_%d_theta_%d_phi_%d.pdf",identifiedParticle.second.c_str(),misIdentifiedParticle.second.c_str(),index_vx,index_vy,index_theta,index_phi);
  c_proj->SaveAs(output_proj.c_str());

  // string input_nsigma = Form("/work/eic/xusun/output/probability/3mm/PID_nSigma_%s.root",date.c_str());
  string input_nsigma = Form("/work/eic/xusun/output/probability/PID_nSigma_%s.root",date.c_str());
  cout << endl;
  cout << "read in file: " << input_nsigma.c_str() << endl;
  TFile *File_InPut_NSigma = TFile::Open(input_nsigma.c_str());
  assert(File_InPut_NSigma);

  string key_nsigma;
  if(rank == 0) key_nsigma = Form("h_mNSigma_%s_%s_vx_%d_vy_%d_theta_%d_phi_%d",identifiedParticle.second.c_str(),misIdentifiedParticle.first.c_str(),index_vx,index_vy,index_theta,index_phi);
  if(rank == 1) key_nsigma = Form("h_mNSigma_%s_%s_vx_%d_vy_%d_theta_%d_phi_%d",identifiedParticle.second.c_str(),misIdentifiedParticle.second.c_str(),index_vx,index_vy,index_theta,index_phi);

  TH1D *h_mNSigma = (TH1D*)File_InPut_NSigma->Get(key_nsigma.c_str());
  assert(h_mNSigma);
  for(int i_bin = 0; i_bin < h_mNSigma->GetNbinsX(); ++i_bin)
  {
    h_mNSigma->SetBinError(i_bin+1,0.0);
  }

  TCanvas *c_nsigma = new TCanvas("c_nsigma","c_nsigma",10,10,800,800);
  c_nsigma->cd(1);
  c_nsigma->cd(1)->SetLeftMargin(0.15);
  c_nsigma->cd(1)->SetBottomMargin(0.15);
  c_nsigma->cd(1)->SetGrid(0,0);
  c_nsigma->cd(1)->SetTicks(1,1);

  h_mNSigma->SetTitle("");
  h_mNSigma->SetStats(0);
  h_mNSigma->GetXaxis()->SetTitle("momentum (GeV/c)");
  h_mNSigma->GetXaxis()->CenterTitle();
  h_mNSigma->GetXaxis()->SetTitleOffset(1.25);
  h_mNSigma->GetXaxis()->SetNdivisions(505);

  h_mNSigma->GetYaxis()->SetTitle("N_{#sigma} #pi-K separation");
  h_mNSigma->GetYaxis()->CenterTitle();
  h_mNSigma->GetYaxis()->SetTitleOffset(1.25);
  h_mNSigma->GetYaxis()->SetNdivisions(505);
  h_mNSigma->GetYaxis()->SetRangeUser(-0.1,13.1);
  h_mNSigma->SetMarkerStyle(24);
  h_mNSigma->SetMarkerSize(1.2);
  h_mNSigma->SetMarkerColor(1);
  h_mNSigma->Draw("P");
  PlotLine(2.5,10.5,3.0,3.0,1,2,2);
  plotTopLegend("3#sigma",9.55,3.1,0.04,1,0.0,42,0,1);

  string output_nsigma;
  if(rank == 0) output_nsigma = Form("../figures/nSigma/c_nSigma_%s_%s_vx_%d_vy_%d_theta_%d_phi_%d.pdf",identifiedParticle.second.c_str(),misIdentifiedParticle.first.c_str(),index_vx,index_vy,index_theta,index_phi);
  if(rank == 1) output_nsigma = Form("../figures/nSigma/c_nSigma_%s_%s_vx_%d_vy_%d_theta_%d_phi_%d.pdf",identifiedParticle.second.c_str(),misIdentifiedParticle.second.c_str(),index_vx,index_vy,index_theta,index_phi);
  c_nsigma->SaveAs(output_nsigma.c_str());
}
