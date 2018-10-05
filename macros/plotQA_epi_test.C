#include "TFile.h"
#include "TString.h"
#include "TTree.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TLegend.h"
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

void plotQA_epi_test()
{
  string input_Oct01 = "/work/eic/xusun/output/probability/PID_nSigma_Oct01_2018.root";
  cout << "read in file: " << input_Oct01.c_str() << endl;
  TFile *File_InPut_Oct01 = TFile::Open(input_Oct01.c_str());
  assert(File_InPut_Oct01);

  string input_Oct03 = "/work/eic/xusun/output/probability/PID_nSigma_Oct03_2018.root";
  cout << "read in file: " << input_Oct03.c_str() << endl;
  TFile *File_InPut_Oct03 = TFile::Open(input_Oct03.c_str());
  assert(File_InPut_Oct03);

  string input_Oct04 = "/work/eic/xusun/output/probability/PID_nSigma_Oct04_2018.root";
  cout << "read in file: " << input_Oct04.c_str() << endl;
  TFile *File_InPut_Oct04 = TFile::Open(input_Oct04.c_str());
  assert(File_InPut_Oct04);

  const int index_vx = 0;
  const int index_vy = 0;
  const int index_theta = 0;
  const int index_phi = 0;

  string key_electron = "h_mNSigma_electron_piminus_vx_0_vy_0_theta_0_phi_0";
  string key_pion = "h_mNSigma_piminus_electron_vx_0_vy_0_theta_0_phi_0";

  TH1D *h_mNSigma_pion_Oct01 = (TH1D*)File_InPut_Oct01->Get(key_pion.c_str());
  assert(h_mNSigma_pion_Oct01);
  TH1D *h_mNSigma_electron_Oct01 = (TH1D*)File_InPut_Oct01->Get(key_electron.c_str());
  assert(h_mNSigma_electron_Oct01);

  TH1D *h_mNSigma_pion_Oct03 = (TH1D*)File_InPut_Oct03->Get(key_pion.c_str());
  assert(h_mNSigma_pion_Oct03);
  TH1D *h_mNSigma_electron_Oct03 = (TH1D*)File_InPut_Oct03->Get(key_electron.c_str());
  assert(h_mNSigma_electron_Oct03);

  TH1D *h_mNSigma_pion_Oct04 = (TH1D*)File_InPut_Oct04->Get(key_pion.c_str());
  assert(h_mNSigma_pion_Oct04);
  TH1D *h_mNSigma_electron_Oct04 = (TH1D*)File_InPut_Oct04->Get(key_electron.c_str());
  assert(h_mNSigma_electron_Oct04);

  TCanvas *c_nsigma = new TCanvas("c_nsigma","c_nsigma",10,10,800,800);
  c_nsigma->cd(1);
  c_nsigma->cd(1)->SetLeftMargin(0.15);
  c_nsigma->cd(1)->SetBottomMargin(0.15);
  c_nsigma->cd(1)->SetGrid(0,0);
  c_nsigma->cd(1)->SetTicks(1,1);

  TH1D *h_play = new TH1D("h_play","h_play",30,0,3.0);
  for(int i_bin = 0; i_bin < 30; ++i_bin)
  {
    h_play->SetBinContent(i_bin+1,-10.0);
    h_play->SetBinError(i_bin+1,1.0);
  }
  h_play->SetTitle("");
  h_play->SetStats(0);
  h_play->GetXaxis()->SetTitle("momentum (GeV/c)");
  h_play->GetXaxis()->CenterTitle();
  h_play->GetXaxis()->SetTitleOffset(1.25);
  h_play->GetXaxis()->SetNdivisions(505);
  h_play->GetXaxis()->SetRangeUser(1.0,3.0);

  h_play->GetYaxis()->SetTitle("N_{#sigma}");
  h_play->GetYaxis()->CenterTitle();
  h_play->GetYaxis()->SetTitleOffset(1.25);
  h_play->GetYaxis()->SetNdivisions(505);
  h_play->GetYaxis()->SetRangeUser(-0.1,9.1);
  h_play->Draw("pE");

  h_mNSigma_pion_Oct01->SetMarkerStyle(24);
  h_mNSigma_pion_Oct01->SetMarkerSize(1.5);
  h_mNSigma_pion_Oct01->SetMarkerColor(2);
  h_mNSigma_pion_Oct01->Draw("pE same");

  h_mNSigma_pion_Oct03->SetMarkerStyle(24);
  h_mNSigma_pion_Oct03->SetMarkerSize(1.5);
  h_mNSigma_pion_Oct03->SetMarkerColor(2);
  h_mNSigma_pion_Oct03->Draw("pE same");

  h_mNSigma_pion_Oct04->SetMarkerStyle(24);
  h_mNSigma_pion_Oct04->SetMarkerSize(1.5);
  h_mNSigma_pion_Oct04->SetMarkerColor(2);
  h_mNSigma_pion_Oct04->Draw("pE same");

  h_mNSigma_electron_Oct01->SetMarkerStyle(20);
  h_mNSigma_electron_Oct01->SetMarkerSize(1.5);
  h_mNSigma_electron_Oct01->SetMarkerColor(kGray+2);
  h_mNSigma_electron_Oct01->Draw("pE same");

  h_mNSigma_electron_Oct03->SetMarkerStyle(20);
  h_mNSigma_electron_Oct03->SetMarkerSize(1.5);
  h_mNSigma_electron_Oct03->SetMarkerColor(kGray+2);
  h_mNSigma_electron_Oct03->Draw("pE same");

  h_mNSigma_electron_Oct04->SetMarkerStyle(20);
  h_mNSigma_electron_Oct04->SetMarkerSize(1.5);
  h_mNSigma_electron_Oct04->SetMarkerColor(kGray+2);
  h_mNSigma_electron_Oct04->Draw("pE same");

  PlotLine(1.0,3.0,3.0,3.0,1,2,2);
  plotTopLegend("3#sigma",2.75,3.1,0.04,1,0.0,42,0,1);

  TLegend *leg = new TLegend(0.6,0.6,0.9,0.8);
  leg->SetFillColor(0);
  leg->SetBorderSize(0);
  leg->AddEntry(h_mNSigma_electron_Oct04,"e^{-}/#pi^{-}","P");
  leg->AddEntry(h_mNSigma_pion_Oct01,"#pi^{-}/e^{-}","P");
  leg->Draw("same");

  c_nsigma->SaveAs("../figures/nSigma/c_nSigma_epi_test.eps");
}
