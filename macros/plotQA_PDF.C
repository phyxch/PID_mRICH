#include "TFile.h"
#include "TString.h"
#include "TTree.h"
#include "TH3D.h"
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

void plotQA_PDF(int pid = 211)
{
  auto particle = get_particle(pid);
  cout << "particle.first = " << particle.first << ", particle.second = " << particle.second << endl;

  if(particle.first == -1) return;

  string date = "test/";
  string input_dir = Form("/work/eic/xusun/output/modular_rich/%s",date.c_str());
  string inputfile = Form("%sPDF_database_test.root",input_dir.c_str());
  cout << "read in file: " << inputfile.c_str() << endl;
  TFile *File_InPut = TFile::Open(inputfile.c_str());

  assert(File_InPut);
  string histQA = Form("h_photonDist_%s",particle.second.c_str());
  TH3D *h_PDF = (TH3D*)File_InPut->Get(histQA.c_str());
  assert(h_PDF);

  TCanvas *c_PDF = new TCanvas("c_PDF","c_PDF",10,10,800,800);
  c_PDF->cd()->SetLeftMargin(0.15);
  c_PDF->cd()->SetBottomMargin(0.15);
  c_PDF->cd()->SetGrid(0,0);
  c_PDF->cd()->SetTicks(1,1);

  h_PDF->SetTitle(particle.second.c_str());
  h_PDF->SetStats(0);

  h_PDF->GetXaxis()->SetTitle("out_x (mm)");
  h_PDF->GetXaxis()->CenterTitle();
  h_PDF->GetXaxis()->SetTitleOffset(1.25);
  h_PDF->GetXaxis()->SetNdivisions(505);

  h_PDF->GetYaxis()->SetTitle("out_y (mm)");
  h_PDF->GetYaxis()->CenterTitle();
  h_PDF->GetYaxis()->SetTitleOffset(1.25);
  h_PDF->GetYaxis()->SetNdivisions(505);

  h_PDF->GetZaxis()->SetTitle("generated momemtum (GeV/c)");
  h_PDF->GetZaxis()->CenterTitle();
  h_PDF->GetZaxis()->SetTitleOffset(1.25);
  h_PDF->GetZaxis()->SetNdivisions(505);

  h_PDF->Draw();

  string outputfile = Form("../figures/masshypo_dist/mass_hypotheses_%s.pdf",particle.second.c_str());
  cout << "save plot to: " << outputfile.c_str() << endl;
  c_PDF->SaveAs(outputfile.c_str());
}
