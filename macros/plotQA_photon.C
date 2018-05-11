#include "TFile.h"
#include "TString.h"
#include "TTree.h"
#include "TH3F.h"
#include "TH2F.h"
#include "TCanvas.h"
#include <cmath>
#include <cassert>
#include <string>
#include <iostream>
#include <utility>
#include <functional>
#include <stdexcept>
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

void plotQA_photon(int pid = 211)
{
  auto particle = get_particle(pid);
  cout << "particle.first = " << particle.first << ", particle.second = " << particle.second << endl;

  if(particle.first == -1) return;

  int runid = 4;

  string date = "May10_2018/";
  string input_dir = Form("/work/eic/xusun/output/modular_rich/%s",date.c_str());
  string inputfile = Form("%sout.%s.%d.root",input_dir.c_str(),particle.second.c_str(),runid);
  cout << "read in file: " << inputfile.c_str() << endl;
  TFile *File_InPut = TFile::Open(inputfile.c_str());

  assert(File_InPut);
  TTree* eic_rich = (TTree *) File_InPut->GetObjectChecked("eic_rich", "TTree");
  assert(eic_rich);

  eic_rich->SetAlias("PID",Form("1*%d",pid));

  TCanvas *c_mRICH = new TCanvas("c_mRICH","c_mRICH",10,10,800,800);
  c_mRICH->cd()->SetLeftMargin(0.15);
  c_mRICH->cd()->SetBottomMargin(0.15);
  c_mRICH->cd()->SetGrid(0,0);
  c_mRICH->cd()->SetTicks(1,1);

  TH3F *h_photon_mRICH = new TH3F("h_photon_mRICH","h_photon_mRICH",230,60.5,299.5,151,-75.5,75.5,151,-75.5,75.5); // 1mm*1mm pixel
  eic_rich->Draw("out_x:out_y:out_z>>h_photon_mRICH","mpid==PID && pid==0");

  h_photon_mRICH->SetTitle(particle.second.c_str());
  h_photon_mRICH->SetStats(0);
  h_photon_mRICH->GetXaxis()->SetTitle("out_z (mm)");
  h_photon_mRICH->GetXaxis()->CenterTitle();
  h_photon_mRICH->GetXaxis()->SetNdivisions(505);
  h_photon_mRICH->GetYaxis()->SetTitle("out_y (mm)");
  h_photon_mRICH->GetYaxis()->CenterTitle();
  h_photon_mRICH->GetYaxis()->SetNdivisions(505);
  h_photon_mRICH->GetZaxis()->SetTitle("out_x (mm)");
  h_photon_mRICH->GetZaxis()->CenterTitle();
  h_photon_mRICH->GetZaxis()->SetNdivisions(505);

  string outputfile = Form("../figures/photon_dist/mRICH_photon_%s.pdf",particle.second.c_str());
  cout << "save plot to: " << outputfile.c_str() << endl;
  c_mRICH->SaveAs(outputfile.c_str());

  TCanvas *c_sensor = new TCanvas("c_sensor","c_sensor",10,10,800,800);
  c_sensor->cd()->SetLeftMargin(0.15);
  c_sensor->cd()->SetBottomMargin(0.15);
  c_sensor->cd()->SetGrid(0,0);
  c_sensor->cd()->SetTicks(1,1);

  TH2F *h_photon_sensor = new TH2F("h_photon_sensor","h_photon_sensor",151,-75.5,75.5,151,-75.5,75.5); // 1mm*1mm pixel
  eic_rich->Draw("out_y:out_x>>h_photon_sensor","mpid==PID && pid==0 && out_z>200.0","colz");

  h_photon_sensor->SetTitle(particle.second.c_str());
  h_photon_sensor->SetStats(0);
  h_photon_sensor->GetXaxis()->SetTitle("out_x (mm)");
  h_photon_sensor->GetXaxis()->CenterTitle();
  h_photon_sensor->GetXaxis()->SetNdivisions(505);
  h_photon_sensor->GetYaxis()->SetTitle("out_y (mm)");
  h_photon_sensor->GetYaxis()->CenterTitle();
  h_photon_sensor->GetYaxis()->SetNdivisions(505);

  outputfile = Form("../figures/photon_dist/sensor_photon_%s.pdf",particle.second.c_str());
  cout << "save plot to: " << outputfile.c_str() << endl;
  c_sensor->SaveAs(outputfile.c_str());
}
