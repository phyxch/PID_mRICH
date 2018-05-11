#include "TFile.h"
#include "TString.h"
#include "TTree.h"
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

void plotQA_Vertex(int pid = 211)
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

  TCanvas *c_rich = new TCanvas("c_rich","c_rich",10,10,800,800);
  c_rich->cd()->SetLeftMargin(0.15);
  c_rich->cd()->SetBottomMargin(0.15);
  c_rich->cd()->SetGrid(0,0);
  c_rich->cd()->SetTicks(1,1);

  TH2F *h_vertex_rich = new TH2F("h_vertex_rich","h_vertex_rich",151,-75.5,75.5,151,-75.5,75.5);
  eic_rich->Draw("mvy:mvx>>h_vertex_rich","mpid==PID","colz");

  h_vertex_rich->SetTitle(particle.second.c_str());
  // h_vertex_rich->SetStats(0);
  h_vertex_rich->GetXaxis()->SetTitle("mvx (mm)");
  h_vertex_rich->GetXaxis()->CenterTitle();
  h_vertex_rich->GetYaxis()->SetTitle("mvy (mm)");
  h_vertex_rich->GetYaxis()->CenterTitle();
  h_vertex_rich->GetYaxis()->SetNdivisions(505);

  string outputfile = Form("../figures/vertex_dist/vertex_rich_%s.pdf",particle.second.c_str());
  cout << "save plot to: " << outputfile.c_str() << endl;
  c_rich->SaveAs(outputfile.c_str());

  assert(File_InPut);
  TTree* generated = (TTree *) File_InPut->GetObjectChecked("generated", "TTree");
  assert(generated);

  generated->SetAlias("PID",Form("1*%d",pid));

  TCanvas *c_generated = new TCanvas("c_generated","c_generated",10,10,800,800);
  c_generated->cd()->SetLeftMargin(0.15);
  c_generated->cd()->SetBottomMargin(0.15);
  c_generated->cd()->SetGrid(0,0);
  c_generated->cd()->SetTicks(1,1);

  TH2F *h_vertex_generated = new TH2F("h_vertex_generated","h_vertex_generated",151,-75.5,75.5,151,-75.5,75.5);
  generated->Draw("vy:vx>>h_vertex_generated","pid==PID","colz");

  h_vertex_generated->SetTitle(particle.second.c_str());
  // h_vertex_generated->SetStats(0);
  h_vertex_generated->GetXaxis()->SetTitle("vx (mm)");
  h_vertex_generated->GetXaxis()->CenterTitle();
  h_vertex_generated->GetYaxis()->SetTitle("vy (mm)");
  h_vertex_generated->GetYaxis()->CenterTitle();
  h_vertex_generated->GetYaxis()->SetNdivisions(505);

  outputfile = Form("../figures/vertex_dist/vertex_generated_%s.pdf",particle.second.c_str());
  cout << "save plot to: " << outputfile.c_str() << endl;
  c_generated->SaveAs(outputfile.c_str());
}
