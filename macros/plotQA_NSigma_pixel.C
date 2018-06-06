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

void plotQA_NSigma_pixel(int pid = 211, int rank = 0)
{
  string date = "May27_2018";

  string input_1mm = Form("/work/eic/xusun/output/probability/1mm/PID_nSigma_%s.root",date.c_str());
  cout << "read in file: " << input_1mm.c_str() << endl;
  TFile *File_InPut_1mm = TFile::Open(input_1mm.c_str());
  assert(File_InPut_1mm);

  string input_2mm = Form("/work/eic/xusun/output/probability/2mm/PID_nSigma_%s.root",date.c_str());
  cout << "read in file: " << input_2mm.c_str() << endl;
  TFile *File_InPut_2mm = TFile::Open(input_2mm.c_str());
  assert(File_InPut_2mm);

  string input_3mm = Form("/work/eic/xusun/output/probability/3mm/PID_nSigma_%s.root",date.c_str());
  cout << "read in file: " << input_3mm.c_str() << endl;
  TFile *File_InPut_3mm = TFile::Open(input_3mm.c_str());
  assert(File_InPut_3mm);

  auto identifiedParticle = get_particle(pid);
  std::pair<std::string,std::string> misIdentifiedParticle = get_misIdentifiedParticle(pid);
  const int index_vx = 0;
  const int index_vy = 0;
  const int index_theta = 0;
  const int index_phi = 3;

  string key_nsigma;
  if(rank == 0) key_nsigma = Form("h_mNSigma_%s_%s_vx_%d_vy_%d_theta_%d_phi_%d",identifiedParticle.second.c_str(),misIdentifiedParticle.first.c_str(),index_vx,index_vy,index_theta,index_phi);
  if(rank == 1) key_nsigma = Form("h_mNSigma_%s_%s_vx_%d_vy_%d_theta_%d_phi_%d",identifiedParticle.second.c_str(),misIdentifiedParticle.second.c_str(),index_vx,index_vy,index_theta,index_phi);

  TH1D *h_mNSigma_1mm = (TH1D*)File_InPut_1mm->Get(key_nsigma.c_str());
  assert(h_mNSigma_1mm);

  TH1D *h_mNSigma_2mm = (TH1D*)File_InPut_2mm->Get(key_nsigma.c_str());
  assert(h_mNSigma_2mm);

  TH1D *h_mNSigma_3mm = (TH1D*)File_InPut_3mm->Get(key_nsigma.c_str());
  assert(h_mNSigma_3mm);
  for(int i_bin = 0; i_bin < h_mNSigma_1mm->GetNbinsX(); ++i_bin)
  {
    h_mNSigma_1mm->SetBinError(i_bin+1,0.0);
    h_mNSigma_2mm->SetBinError(i_bin+1,0.0);
    h_mNSigma_3mm->SetBinError(i_bin+1,0.0);
  }

  TCanvas *c_nsigma = new TCanvas("c_nsigma","c_nsigma",10,10,800,800);
  c_nsigma->cd(1);
  c_nsigma->cd(1)->SetLeftMargin(0.15);
  c_nsigma->cd(1)->SetBottomMargin(0.15);
  c_nsigma->cd(1)->SetGrid(0,0);
  c_nsigma->cd(1)->SetTicks(1,1);

  h_mNSigma_1mm->SetTitle("");
  h_mNSigma_1mm->SetStats(0);
  h_mNSigma_1mm->GetXaxis()->SetTitle("momentum (GeV/c)");
  h_mNSigma_1mm->GetXaxis()->CenterTitle();
  h_mNSigma_1mm->GetXaxis()->SetTitleOffset(1.25);
  h_mNSigma_1mm->GetXaxis()->SetNdivisions(505);

  h_mNSigma_1mm->GetYaxis()->SetTitle("N_{#sigma} #pi-K separation");
  h_mNSigma_1mm->GetYaxis()->CenterTitle();
  h_mNSigma_1mm->GetYaxis()->SetTitleOffset(1.25);
  h_mNSigma_1mm->GetYaxis()->SetNdivisions(505);
  h_mNSigma_1mm->GetYaxis()->SetRangeUser(-0.1,15.1);
  h_mNSigma_1mm->SetMarkerStyle(28);
  h_mNSigma_1mm->SetMarkerSize(1.4);
  h_mNSigma_1mm->SetMarkerColor(1);
  h_mNSigma_1mm->Draw("P");

  h_mNSigma_2mm->SetMarkerStyle(25);
  h_mNSigma_2mm->SetMarkerSize(1.4);
  h_mNSigma_2mm->SetMarkerColor(4);
  h_mNSigma_2mm->Draw("P same");

  h_mNSigma_3mm->SetMarkerStyle(20);
  h_mNSigma_3mm->SetMarkerSize(1.4);
  h_mNSigma_3mm->SetMarkerColor(2);
  h_mNSigma_3mm->Draw("P same");

  PlotLine(2.5,10.5,3.0,3.0,1,2,2);
  plotTopLegend("3#sigma",3.55,3.1,0.04,1,0.0,42,0,1);

  TLegend *leg = new TLegend(0.6,0.6,0.9,0.8);
  leg->SetFillColor(0);
  leg->SetBorderSize(0);
  leg->AddEntry(h_mNSigma_1mm,"1mm*1mm","P");
  leg->AddEntry(h_mNSigma_2mm,"2mm*2mm","P");
  leg->AddEntry(h_mNSigma_3mm,"3mm*3mm","P");
  leg->Draw("same");

  string output_nsigma;
  if(rank == 0) output_nsigma = Form("../figures/nSigma/c_nSigmaPixel_%s_%s_vx_%d_vy_%d_theta_%d_phi_%d.pdf",identifiedParticle.second.c_str(),misIdentifiedParticle.first.c_str(),index_vx,index_vy,index_theta,index_phi);
  if(rank == 1) output_nsigma = Form("../figures/nSigma/c_nSigmaPixel_%s_%s_vx_%d_vy_%d_theta_%d_phi_%d.pdf",identifiedParticle.second.c_str(),misIdentifiedParticle.second.c_str(),index_vx,index_vy,index_theta,index_phi);
  c_nsigma->SaveAs(output_nsigma.c_str());
}
