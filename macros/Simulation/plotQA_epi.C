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
#include "../draw.h"
#include "../../include/mRICH.h"
#include <fstream>

using namespace std;

void plotQA_epi()
{
  string date = "Oct08_2018";
  string inputfile = Form("/work/eic/xusun/output/probability/PID_nSigma_%s.root",date.c_str());
  cout << "read in file: " << inputfile.c_str() << endl;
  TFile *File_InPut = TFile::Open(inputfile.c_str());
  assert(File_InPut);

  const int index_vx = 0;
  const int index_vy = 0;
  const int index_theta = 0;
  const int index_phi = 0;

  string key_electron = "h_mNSigma_electron_piminus_vx_0_vy_0_theta_0_phi_0";
  TH1D *h_mNSigma_electron = (TH1D*)File_InPut->Get(key_electron.c_str());
  assert(h_mNSigma_electron);

  string key_pion = "h_mNSigma_piminus_electron_vx_0_vy_0_theta_0_phi_0";
  TH1D *h_mNSigma_pion = (TH1D*)File_InPut->Get(key_pion.c_str());
  assert(h_mNSigma_pion);

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
  h_play->GetXaxis()->SetRangeUser(0.5,3.0);

  h_play->GetYaxis()->SetTitle("N_{#sigma}");
  h_play->GetYaxis()->CenterTitle();
  h_play->GetYaxis()->SetTitleOffset(1.25);
  h_play->GetYaxis()->SetNdivisions(505);
  h_play->GetYaxis()->SetRangeUser(-0.1,9.1);
  h_play->Draw("pE");

  h_mNSigma_electron->SetMarkerStyle(20);
  h_mNSigma_electron->SetMarkerSize(1.5);
  h_mNSigma_electron->SetMarkerColor(kGray+2);
  h_mNSigma_electron->SetLineColor(1);
  h_mNSigma_electron->Draw("pE same");

  h_mNSigma_pion->SetMarkerStyle(24);
  h_mNSigma_pion->SetMarkerSize(1.5);
  h_mNSigma_pion->SetMarkerColor(2);
  h_mNSigma_pion->SetLineColor(2);
  h_mNSigma_pion->Draw("pE same");

  PlotLine(0.5,3.0,3.0,3.0,1,2,2);
  plotTopLegend("3#sigma",2.75,3.1,0.04,1,0.0,42,0,1);

  TLegend *leg = new TLegend(0.6,0.6,0.7,0.8);
  leg->SetFillColor(0);
  leg->SetBorderSize(0);
  leg->AddEntry(h_mNSigma_electron,"e^{-}/#pi^{-}","P");
  leg->AddEntry(h_mNSigma_pion,"#pi^{-}/e^{-}","P");
  leg->Draw("same");

  c_nsigma->SaveAs("../../figures/nSigma/c_nSigma_epi.eps");

  TFile *File_OutPut = new TFile("./epi_mRICH.root","RECREATE");
  File_OutPut->cd();
  h_mNSigma_electron->Write();
  File_OutPut->Close();

  // save e/pi and pi/e separation into a txt file
  string output_epi = "./output_epi.txt";
  cout << "save e/pi separation power to " << output_epi.c_str() << endl;
  ofstream File_OutPut_epi(output_epi.c_str());
  int NumOfBin_epi = h_mNSigma_electron->GetNbinsX();
  File_OutPut_epi << "momentum    " << "N_sigma    " << "errors" << endl;
  for(int i_bin = 0; i_bin < NumOfBin_epi; ++i_bin)
  {
    float p = h_mNSigma_electron->GetBinCenter(i_bin+1);
    float nsigma = h_mNSigma_electron->GetBinContent(i_bin+1);
    float error = h_mNSigma_electron->GetBinError(i_bin+1);
    File_OutPut_epi << p << "    " << nsigma << "    " << error << endl;
  }
  File_OutPut_epi.close();

  string output_pie = "./output_pie.txt";
  cout << "save pi/e separation power to " << output_pie.c_str() << endl;
  ofstream File_OutPut_pie(output_pie.c_str());
  int NumOfBin_pie = h_mNSigma_pion->GetNbinsX();
  File_OutPut_pie << "momentum    " << "N_sigma    " << "errors" << endl;
  for(int i_bin = 0; i_bin < NumOfBin_pie; ++i_bin)
  {
    float p = h_mNSigma_pion->GetBinCenter(i_bin+1);
    float nsigma = h_mNSigma_pion->GetBinContent(i_bin+1);
    float error = h_mNSigma_pion->GetBinError(i_bin+1);
    File_OutPut_pie << p << "    " << nsigma << "    " << error << endl;
  }
  File_OutPut_pie.close();
}
